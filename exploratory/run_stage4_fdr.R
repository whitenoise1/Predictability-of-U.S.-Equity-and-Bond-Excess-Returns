# working_paper/run_stage4_fdr.R
# Stage 4 / S4.1 (P11 / P-D) — universe-wide scalar-CW BH/BY-FDR over the P10
# trial universe. The load-bearing multiple-testing test: does the bond 9.6%
# P10 acceptance survive, or is it data-snooping? RUN FROM REPO ROOT.
#
# Per evaluable cell the published OOS Clark-West point stat is referenced to a
# Student-t with df = n_eff - 1, n_eff = floor(n_oos / h) (PD4 "primary OOS"
# leg; identical to cw_overlap_robust()$p_t_eff but computed scalar-only from
# the stored cw / n_oos / horizon). One-sided (CW>0 => model beats benchmark).
#
# BH + BY (headline, dependence-robust) at q<0.05 over NESTED families:
#   self  (per-target) : equity ; bond_<maturity> (5 — each maturity a
#                        pre-registered target, [[feedback_per_target_factor_pools]])
#   pool  (pooled-tgt)  : equity ; pooled-bond (charges the 5-maturity snooping tax)
#   univ  (universe)    : all evaluable cells (PD4 literal "universe-wide", strictest)
#
# Bond has no retained forecasts => FDR only here. A raw-signal "looks alive"
# flag (CW p<0.01 AND r2_oos>0 cluster, per maturity) gates the S4.3 bond
# forecast-retention re-run + bond SPA (SPA and FDR are NOT nested — WS4 had 0
# FDR survivors yet a SPA-superior best cell, so we trigger on raw signal).
#
# Gates: G1 reproduce stored accept_p10 EXACTLY from r2/cw/sharpe/threshold;
#        no-NA-p stop guard (no silent denominator shrinkage); df>=1;
#        family-count asserts (3438 evaluable; bond 2874/276, equity 564/21).
suppressMessages({library(arrow); library(dplyr)})
source("R/overlap_robust_inference.R")   # bh_fdr, by_fdr, overlap_effective_n
source("R/target_relative_sharpe.R")     # accept_decision_rel

RUNDATE_IN <- "2026-06-12"
ALIVE_CLUSTER <- 3L                       # >= this many raw-alive cells => "cluster"

p10 <- read_parquet(sprintf("data/audit/p10_rescored_%s.parquet", RUNDATE_IN))
fam <- p10 |> filter(!is.na(accept_p10))  # evaluable = P10-acceptance basis
cat(sprintf("[load] %d cells; %d evaluable (accept_p10 non-NA)\n", nrow(p10), nrow(fam)))

# ---- G1: reproduce stored accept_p10 exactly -------------------------------
acc_rep <- as.logical(mapply(accept_decision_rel, fam$r2_oos, fam$cw, fam$sharpe, fam$threshold))
bad <- which(is.na(acc_rep) | acc_rep != fam$accept_p10)
if (length(bad)) stop(sprintf("[G1 FAIL] %d cells: recomputed accept_p10 != stored", length(bad)))
cat(sprintf("[G1 OK] accept_p10 reproduced on all %d evaluable cells\n", nrow(fam)))

# ---- family-count asserts (the headline denominator) -----------------------
ct <- fam |> group_by(target) |> summarise(n = n(), acc = sum(accept_p10), .groups = "drop")
cat("[counts] "); print(ct)
stopifnot(nrow(fam) == 3438L,
          ct$n[ct$target == "bond"]   == 2874L, ct$acc[ct$target == "bond"]   == 276L,
          ct$n[ct$target == "equity"] == 564L,  ct$acc[ct$target == "equity"] == 21L)

# ---- scalar p_t_eff + no-NA guard ------------------------------------------
if (any(is.na(fam$cw) | is.na(fam$n_oos) | is.na(fam$horizon)))
  stop("NA in cw / n_oos / horizon on an evaluable cell")
fam <- fam |> mutate(
  n_eff = vapply(seq_len(n()), function(i) overlap_effective_n(n_oos[i], horizon[i]), integer(1)),
  df    = pmax(1L, n_eff - 1L),
  p_cw  = pt(cw, df = df, lower.tail = FALSE))
if (any(is.na(fam$p_cw))) stop("NA p-value produced (would silently shrink the FDR denominator)")

# ---- nested FDR families ----------------------------------------------------
fam <- fam |>
  mutate(fam_self = if_else(target == "equity", "equity", maturity),
         fam_pool = if_else(target == "equity", "equity", "bond")) |>
  group_by(fam_self) |> mutate(q_bh_self = bh_fdr(p_cw), q_by_self = by_fdr(p_cw)) |> ungroup() |>
  group_by(fam_pool) |> mutate(q_bh_pool = bh_fdr(p_cw), q_by_pool = by_fdr(p_cw)) |> ungroup() |>
  mutate(q_bh_univ = bh_fdr(p_cw), q_by_univ = by_fdr(p_cw))

# ---- persist per-cell ------------------------------------------------------
out <- file.path("data/audit", sprintf("stage4_fdr_%s.parquet", Sys.Date()))
write_parquet(fam, out); cat(sprintf("[write] %s (%d rows)\n", out, nrow(fam)))

# ---- summaries: survivors per family unit ----------------------------------
summ <- function(df, qbh, qby, label) tibble(
  family = label, m = nrow(df), raw_cw_p05 = sum(df$p_cw < 0.05),
  accept_p10 = sum(df$accept_p10), bh_surv = sum(df[[qbh]] < 0.05),
  by_surv = sum(df[[qby]] < 0.05), min_p = signif(min(df$p_cw), 3),
  min_q_bh = signif(min(df[[qbh]]), 3), min_q_by = signif(min(df[[qby]]), 3))

self_tbl <- bind_rows(lapply(split(fam, fam$fam_self), summ,
                             qbh = "q_bh_self", qby = "q_by_self", label = NULL),
                      .id = "family")
pool_tbl <- bind_rows(lapply(split(fam, fam$fam_pool), summ,
                             qbh = "q_bh_pool", qby = "q_by_pool", label = NULL),
                      .id = "family")
univ_tbl <- summ(fam, "q_bh_univ", "q_by_univ", "universe")

cat("\n=== self / per-target families (equity ; bond_<maturity>) ===\n"); print(as.data.frame(self_tbl))
cat("\n=== pooled-target families (equity ; pooled-bond) ===\n");          print(as.data.frame(pool_tbl))
cat("\n=== pooled universe (strictest) ===\n");                            print(as.data.frame(univ_tbl))

# ---- raw-signal "looks alive" gate for S4.3 (bond) -------------------------
alive <- fam |> filter(target == "bond") |> group_by(maturity) |>
  summarise(n_raw_alive = sum(p_cw < 0.01 & r2_oos > 0),
            min_p = signif(min(p_cw), 3), .groups = "drop") |>
  mutate(cluster = n_raw_alive >= ALIVE_CLUSTER)
bond_looks_alive <- any(alive$cluster)
cat(sprintf("\n=== S4.3 raw-signal gate (CW p<0.01 & r2_oos>0; cluster = >=%d) ===\n", ALIVE_CLUSTER))
print(as.data.frame(alive))

# ---- headline ---------------------------------------------------------------
bond_self_surv <- sum(self_tbl$by_surv[self_tbl$family != "equity"])
cat("\n================= HEADLINE =================\n")
cat(sprintf("BOND  9.6%% (276 P10-accepts): BH survivors = self %d / pooled-bond %d / universe %d\n",
            sum(self_tbl$bh_surv[self_tbl$family != "equity"]),
            pool_tbl$bh_surv[pool_tbl$family == "bond"], univ_tbl$bh_surv))
cat(sprintf("BOND  BY survivors (headline): self(per-mat) %d / pooled-bond %d / universe %d\n",
            bond_self_surv, pool_tbl$by_surv[pool_tbl$family == "bond"],
            univ_tbl$by_surv))
cat(sprintf("EQUITY 3.7%% (21 P10-accepts): BY survivors self/pooled %d / universe %d\n",
            self_tbl$by_surv[self_tbl$family == "equity"], univ_tbl$by_surv))
cat(sprintf("\nS4.3 bond forecast-retention + bond SPA triggered? %s\n",
            if (bond_looks_alive) "YES — a bond maturity shows a raw-signal cluster" else "NO — no raw-signal cluster; bond closes FDR-only"))
cat("[run_stage4_fdr complete]\n")

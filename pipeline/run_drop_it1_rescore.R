# working_paper/run_drop_it1_rescore.R
# ITERATION-1 REMOVAL (GS+PCA cautionary) — re-score the manuscript universe with
# the OLS-class data-driven "iteration-1" trials DROPPED, from COMMITTED artifacts
# only. NO model is re-fit (C-PCR is dormant in production); every number here is a
# re-aggregation / re-correction of the existing p10 + Stage-4 inputs. RUN FROM REPO ROOT.
#
# it1 := scenario in {S2_bond, S2_equity} AND model_config in {C-OLS, C-COMB, C-CT}
#        (the raw-input OLS path the audit relegates to a cautionary specimen; see
#        METHODS_NOTE_PCR_CAUTIONARY.md and METHOD_AUDIT_GS_PCA.md §7 option B).
#
# Three legs, mirroring the committed drivers with the it1 filter inserted:
#   (1) P10 rates          — re-aggregate p10_rescored over the kept (non-it1) cells.
#   (2) BH/BY-FDR          — RECOMPUTE p.adjust within each (now smaller) family from
#                            the per-cell p_cw (stored q_* were computed on the full
#                            universe and MUST NOT be copied). Mirrors run_stage4_fdr.R.
#   (3) equity-S2 SPA      — drop the OLS-class (it1) equity cells from the retained
#                            forecasts, re-run Hansen SPA per (h x sample). Mirrors
#                            run_stage4_spa.R; G1 reproduces stored losses < 1e-8.
#
# Emits data/audit/{p10_rescored,stage4_fdr,stage4_equity_inference,stage4_spa}_drop-it1_<date>.parquet
# and prints a full BEFORE/AFTER scalar table for every figure/inline number the .Rmd cites.
suppressMessages({library(arrow); library(dplyr); library(tidyr); library(purrr)})
source("R/overlap_robust_inference.R")   # bh_fdr, by_fdr, overlap_effective_n, cw_overlap_robust, hodrick_1b_mz, nw_hac_mz
source("R/target_relative_sharpe.R")     # accept_decision_rel
source("R/predictive_metrics.R"); source("R/predictive_regression.R")
source("R/spa_test.R")
source("R/pit_query.R"); source("R/panel.R"); source("R/targets.R")

RUNDATE_IN <- "2026-06-12"
FDR_IN     <- "2026-06-13"                # committed stage4_fdr (for the BEFORE column)
OUT_TAG    <- sprintf("drop-it1_%s", Sys.Date())
AS_OF      <- as.Date("2026-04-30")
SPA_B      <- 2000L
KEY <- c("chapter","scenario","model_config","iteration","horizon",
         "window_kind","window_length","sample")

is_it1 <- function(scenario, model_config)
  scenario %in% c("S2_bond","S2_equity") & model_config %in% c("C-OLS","C-COMB","C-CT")

`%||%` <- function(a, b) if (length(a)) a else b
hr <- function() cat(strrep("-", 78), "\n")

# ============================================================================
# (1) P10 RATES — re-aggregate over kept (non-it1) cells
# ============================================================================
p10_all <- read_parquet(sprintf("data/audit/p10_rescored_%s.parquet", RUNDATE_IN))
ev <- p10_all |> filter(!is.na(accept_p10))                 # evaluable basis (same as the .Rmd)

# G1a: reproduce stored accept_p10 exactly (same gate as run_stage4_fdr.R)
acc_rep <- as.logical(mapply(accept_decision_rel, ev$r2_oos, ev$cw, ev$sharpe, ev$threshold))
if (any(is.na(acc_rep) | acc_rep != ev$accept_p10))
  stop("[G1a FAIL] recomputed accept_p10 != stored")
cat(sprintf("[G1a OK] accept_p10 reproduced on all %d evaluable cells\n", nrow(ev)))

ev <- ev |> mutate(it1 = is_it1(scenario, model_config))
kept <- ev |> filter(!it1)

# headline-count asserts (verified against artifacts before coding — no silent drift)
ct_all  <- ev   |> count(target, wt = NULL, name = "n")
acc_all <- ev   |> group_by(target) |> summarise(acc = sum(accept_p10), .groups = "drop")
ct_keep <- kept |> group_by(target) |> summarise(n = n(), acc = sum(accept_p10), .groups = "drop")
stopifnot(nrow(ev) == 3438L, sum(ev$it1) == 1632L, nrow(kept) == 1806L,
          ct_keep$n[ct_keep$target == "bond"]   == 1530L, ct_keep$acc[ct_keep$target == "bond"]   == 76L,
          ct_keep$n[ct_keep$target == "equity"] == 276L,  ct_keep$acc[ct_keep$target == "equity"] == 0L)
cat("[asserts OK] evaluable 3438 -> 1806 (it1=1632); bond 76/1530; equity 0/276\n")

# emit the filtered FULL p10 (incl. below-floor non-it1 rows) so the .Rmd setup
# can repoint p10_all and keep its own evaluable/ch2/ch3 filtering unchanged.
p10_drop <- p10_all |> filter(!is_it1(scenario, model_config))
write_parquet(p10_drop, sprintf("data/audit/p10_rescored_%s.parquet", OUT_TAG))

# rate helpers identical to the .Rmd (evaluable basis, accept := accept_p10)
rate <- function(d) if (nrow(d)) sprintf("%.1f", 100 * sum(d$accept_p10) / nrow(d)) else "NA"
rate_pair <- function(filt) {
  a <- ev   |> filter(filt(scenario, target, maturity, horizon, model_config))
  k <- kept |> filter(filt(scenario, target, maturity, horizon, model_config))
  c(before = rate(a), after = rate(k), n_before = nrow(a), n_after = nrow(k))
}

# best bond trial: DECISION = accepted + positive r^2, max CW (the headline replacement)
best_old <- ev   |> filter(target == "bond", accept_p10) |> arrange(desc(cw)) |> slice(1)
best_new <- kept |> filter(target == "bond", accept_p10, r2_oos > 0) |> arrange(desc(cw)) |> slice(1)

# ============================================================================
# (2) BH/BY-FDR — recompute p.adjust within each family on the kept set
# ============================================================================
mk_fdr <- function(df) {
  if (any(is.na(df$cw) | is.na(df$n_oos) | is.na(df$horizon)))
    stop("NA in cw / n_oos / horizon on an evaluable cell")
  df <- df |> mutate(
    n_eff = vapply(seq_len(n()), function(i) overlap_effective_n(n_oos[i], horizon[i]), integer(1)),
    df    = pmax(1L, n_eff - 1L),
    p_cw  = pt(cw, df = df, lower.tail = FALSE),
    fam_self = if_else(target == "equity", "equity", maturity),
    fam_pool = if_else(target == "equity", "equity", "bond"))
  if (any(is.na(df$p_cw))) stop("NA p-value produced (would silently shrink the FDR denominator)")
  df |>
    group_by(fam_self) |> mutate(q_bh_self = bh_fdr(p_cw), q_by_self = by_fdr(p_cw)) |> ungroup() |>
    group_by(fam_pool) |> mutate(q_bh_pool = bh_fdr(p_cw), q_by_pool = by_fdr(p_cw)) |> ungroup() |>
    mutate(q_bh_univ = bh_fdr(p_cw), q_by_univ = by_fdr(p_cw))
}
fdr_keep <- mk_fdr(kept)
write_parquet(fdr_keep, sprintf("data/audit/stage4_fdr_%s.parquet", OUT_TAG))

fdr_old <- read_parquet(sprintf("data/audit/stage4_fdr_%s.parquet", FDR_IN))  # committed (full universe)
minq <- function(df, col) signif(min(df[[col]]), 3)
nsurv <- function(df, col) sum(df[[col]] < 0.05)

fdr_row <- function(label, dq_old, dq_new, surv_old, surv_new)
  tibble(family = label, by_q_before = dq_old, by_q_after = dq_new,
         by_surv_before = surv_old, by_surv_after = surv_new)
fdr_tbl <- bind_rows(
  fdr_row("bond_1y (self)",
          minq(filter(fdr_old, fam_self == "bond_1y"), "q_by_self"),
          minq(filter(fdr_keep, fam_self == "bond_1y"), "q_by_self"),
          nsurv(filter(fdr_old, fam_self == "bond_1y"), "q_by_self"),
          nsurv(filter(fdr_keep, fam_self == "bond_1y"), "q_by_self")),
  fdr_row("pooled-bond",
          minq(filter(fdr_old, fam_pool == "bond"), "q_by_pool"),
          minq(filter(fdr_keep, fam_pool == "bond"), "q_by_pool"),
          nsurv(filter(fdr_old, fam_pool == "bond"), "q_by_pool"),
          nsurv(filter(fdr_keep, fam_pool == "bond"), "q_by_pool")),
  fdr_row("equity (self)",
          minq(filter(fdr_old, fam_self == "equity"), "q_by_self"),
          minq(filter(fdr_keep, fam_self == "equity"), "q_by_self"),
          nsurv(filter(fdr_old, fam_self == "equity"), "q_by_self"),
          nsurv(filter(fdr_keep, fam_self == "equity"), "q_by_self")),
  fdr_row("universe",
          minq(fdr_old, "q_by_univ"), minq(fdr_keep, "q_by_univ"),
          nsurv(fdr_old, "q_by_univ"), nsurv(fdr_keep, "q_by_univ")))

# BH (lenient) bond short-mat survivors — the .Rmd reports these as the contrast
bh_before <- c(b1 = nsurv(filter(fdr_old, fam_self == "bond_1y"), "q_bh_self"),
               pool = nsurv(filter(fdr_old, fam_pool == "bond"), "q_bh_pool"),
               univ = nsurv(fdr_old, "q_bh_univ"))
bh_after  <- c(b1 = nsurv(filter(fdr_keep, fam_self == "bond_1y"), "q_bh_self"),
               pool = nsurv(filter(fdr_keep, fam_pool == "bond"), "q_bh_pool"),
               univ = nsurv(fdr_keep, "q_bh_univ"))
# n_live_bond guard used in the .Rmd setup (line 111) — recompute for part C
nlive_before <- sum(fdr_old$p_cw < 0.01 & fdr_old$r2_oos > 0 & grepl("bond", fdr_old$fam_self))
nlive_after  <- sum(fdr_keep$p_cw < 0.01 & fdr_keep$r2_oos > 0 & grepl("bond", fdr_keep$fam_self))
by_any_after <- any(fdr_keep$q_by_univ < 0.05) | any(fdr_keep$q_by_self < 0.05) | any(fdr_keep$q_by_pool < 0.05)

# ============================================================================
# (3) EQUITY-S2 SPA — drop it1 from the retained forecasts, re-run Hansen SPA
# ============================================================================
fc  <- read_parquet(sprintf("data/audit/equity_s2_forecasts_%s.parquet", RUNDATE_IN)) |>
  filter(!is_it1(scenario, model_config))
p10e <- p10_all |> filter(target == "equity", grepl("S2", scenario), !is.na(accept_p10),
                          !is_it1(scenario, model_config))
n_cells <- nrow(distinct(fc, across(all_of(KEY))))
cat(sprintf("[SPA load] %d non-it1 forecast cells; %d non-it1 equity-S2 evaluable in p10\n",
            n_cells, nrow(p10e)))
stopifnot(n_cells == nrow(p10e), n_cells == 144L)           # 432 - 288 it1 = 144 shrinkage cells

panel  <- build_panel(AS_OF, c("SP500TR","DGS3MO"), "M", "forward")
r1_tbl <- tibble(as_of_date = panel$reference_date,
                 r1 = construct_equity_excess_return(panel, h = 1L, sp_col = "SP500TR"))

fc <- fc |> arrange(across(all_of(KEY)), as_of_date)
inf <- fc |> group_by(across(all_of(KEY))) |> group_split() |> map_dfr(function(d) {
  h  <- d$horizon[1L]
  r2 <- r2_oos(d$y_hat, d$y_realized, d$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(d$y_hat, d$y_realized, d$y_bench), L = floor(1.5 * h))$stat
  sh  <- signal_sharpe(d$y_hat, d$y_realized, h = h)
  cro <- cw_overlap_robust(d$y_hat, d$y_realized, d$y_bench, h)
  dr  <- left_join(d, r1_tbl, by = "as_of_date")
  if (anyNA(dr$r1)) stop("r1 join left NA for an as-of date (alignment bug)")
  hod <- hodrick_1b_mz(dr$y_realized, dr$y_hat, dr$r1, h)
  nw  <- nw_hac_mz(dr$y_realized, dr$y_hat, h)
  d[1L, KEY] |> bind_cols(tibble(
    n_oos = nrow(d), r2_oos = r2, cw = cw, sharpe = sh,
    n_eff = cro$n_eff, cw_overlap = cro$cw_overlap, p_normal = cro$p_normal,
    p_t_eff = cro$p_t_eff, beta_mz = hod$beta, p_hodrick = hod$p_hodrick, p_nwhac = nw$p_nwhac))
})
if (nrow(filter(inf, is.na(p_t_eff) | is.na(p_hodrick) | is.na(p_nwhac) | is.na(p_normal))))
  stop("NA inference on an evaluable equity cell")

# G1: exact reproduction of stored r2/cw/sharpe vs p10
chk <- inf |> select(all_of(KEY), r2_oos, cw, sharpe) |>
  inner_join(p10e |> select(all_of(KEY), r2_s = r2_oos, cw_s = cw, sh_s = sharpe), by = KEY)
if (nrow(chk) != nrow(inf)) stop(sprintf("G1: %d re-run cells failed to match a p10 row", nrow(inf) - nrow(chk)))
g1 <- chk |> summarise(dr2 = max(abs(r2_oos - r2_s)), dcw = max(abs(cw - cw_s)), dsh = max(abs(sharpe - sh_s)))
cat(sprintf("[G1 OK] %d cells; max|dr2|=%.2e max|dcw|=%.2e max|dsharpe|=%.2e\n", nrow(chk), g1$dr2, g1$dcw, g1$dsh))
stopifnot(g1$dr2 < 1e-8, g1$dcw < 1e-8, g1$dsh < 1e-8)

inf <- inf |> left_join(p10e |> select(all_of(KEY), threshold, accept_p10), by = KEY) |>
  mutate(q_bh_oos = bh_fdr(p_t_eff), q_by_oos = by_fdr(p_t_eff),
         q_bh_is  = bh_fdr(p_hodrick), q_by_is  = by_fdr(p_hodrick))
write_parquet(inf, sprintf("data/audit/stage4_equity_inference_%s.parquet", OUT_TAG))

spa_family <- function(h, smp, fc_sub, family) {
  sub <- fc_sub |> filter(horizon == h, sample == smp)
  per <- sub |> group_by(across(all_of(KEY))) |>
    group_map(~ tibble(as_of_date = .x$as_of_date,
                       d = (.x$y_realized - .x$y_bench)^2 - (.x$y_realized - .x$y_hat)^2))
  row <- function(m_orig, m_used, n_common, n_dropped, o, reason) tibble(
    family = family, horizon = h, sample = smp, m_orig = m_orig, m_used = m_used,
    n_common = n_common, n_dropped = n_dropped, q = 1 / h, B = SPA_B,
    T_spa = o$T_spa %||% NA_real_, best_z = o$best_z %||% NA_real_,
    p_lower = o$p_lower %||% NA_real_, p_consistent = o$p_consistent %||% NA_real_,
    p_upper = o$p_upper %||% NA_real_, reason = reason)
  if (length(per) < 1L) return(row(0L, 0L, 0L, 0L, list(), "no cells"))
  common <- sort(as.Date(reduce(map(per, ~ as.character(.x$as_of_date)), intersect)))
  if (length(common) < 10L)
    return(row(length(per), 0L, length(common), 0L, list(), sprintf("common axis < 10 (%d)", length(common))))
  mat   <- vapply(per, function(p) p$d[match(common, p$as_of_date)], numeric(length(common)))
  sds   <- apply(mat, 2L, stats::sd)
  ndrop <- sum(sds == 0)
  mat   <- mat[, sds > 0, drop = FALSE]
  if (ncol(mat) < 1L) return(row(length(per), 0L, length(common), ndrop, list(), "all columns benchmark-identical"))
  o <- spa_test(mat, q = 1 / h, B = SPA_B, seed = 11L)
  row(length(per), o$m, length(common), ndrop, o, "")
}
combos <- fc |> distinct(horizon, sample)
spa <- bind_rows(
  pmap_dfr(combos, function(horizon, sample) spa_family(horizon, sample, fc, "all_windows")),
  pmap_dfr(combos, function(horizon, sample)
    spa_family(horizon, sample, filter(fc, window_kind == "expanding"), "expanding_only")))
write_parquet(spa, sprintf("data/audit/stage4_spa_%s.parquet", OUT_TAG))
spa_hits <- spa |> filter(!is.na(p_consistent), p_consistent < 0.05)

# ============================================================================
# BEFORE / AFTER REPORT — every figure / inline scalar the .Rmd cites
# ============================================================================
cat("\n"); hr(); cat("DROP-IT1 RE-SCORE — BEFORE / AFTER\n"); hr()

cat("\n[1] Evaluable universe & headline acceptance (evaluable basis, P10 rule)\n")
cat(sprintf("    evaluable trials : 3438 -> 1806   (it1 dropped = 1632, 47.5%%)\n"))
cat(sprintf("    BOND   acceptance: %s%% (%d/%d) -> %s%% (%d/%d)\n",
            rate(filter(ev, target=="bond")),  sum(filter(ev,target=="bond")$accept_p10),  nrow(filter(ev,target=="bond")),
            rate(filter(kept,target=="bond")), sum(filter(kept,target=="bond")$accept_p10), nrow(filter(kept,target=="bond"))))
cat(sprintf("    EQUITY acceptance: %s%% (%d/%d) -> %s%% (%d/%d)\n",
            rate(filter(ev, target=="equity")),  sum(filter(ev,target=="equity")$accept_p10),  nrow(filter(ev,target=="equity")),
            rate(filter(kept,target=="equity")), sum(filter(kept,target=="equity")$accept_p10), nrow(filter(kept,target=="equity"))))

cat("\n[2] Bond acceptance by maturity (Fig 2; ch2 level pools)\n")
bymat <- function(d) d |> filter(target=="bond") |> group_by(maturity) |>
  summarise(n=n(), acc=sum(accept_p10), rate=sprintf("%.1f", 100*mean(accept_p10)), .groups="drop")
print(full_join(bymat(filter(ev, chapter=="ch2"))  |> select(maturity, before=rate, nb=n),
                bymat(filter(kept, chapter=="ch2")) |> select(maturity, after=rate, na=n), by="maturity") |> as.data.frame())

cat("\n[3] Scenario-2 by horizon (Part I.1 prose)\n")
sc_h <- function(d, sc, h) { x <- d |> filter(scenario==sc, horizon==h); sprintf("%s%% (%d/%d)", rate(x), sum(x$accept_p10), nrow(x)) }
for (sc in c("S2_bond","S2_equity")) for (h in c(1,3))
  cat(sprintf("    %-10s h=%d : %s -> %s\n", sc, h, sc_h(ev,sc,h), sc_h(kept,sc,h)))

cat("\n[4] Best bond trial (headline; decision = accepted + positive r2, max CW)\n")
cat(sprintf("    BEFORE (it1): %s %s h=%d %s-%s %s | cw=%.2f r2=%.3f\n",
            best_old$model_config, best_old$maturity, best_old$horizon, best_old$window_kind,
            best_old$window_length, best_old$sample, best_old$cw, best_old$r2_oos))
cat(sprintf("    AFTER (non-it1): %s %s h=%d %s-%s %s | cw=%.2f r2=%.3f\n",
            best_new$model_config, best_new$maturity, best_new$horizon, best_new$window_kind,
            best_new$window_length, best_new$sample, best_new$cw, best_new$r2_oos))

cat("\n[5] Stage-4 BH/BY-FDR (headline = BY, dependence-robust)\n")
print(as.data.frame(fdr_tbl))
cat(sprintf("    bond_1y family size: %d -> %d\n",
            nrow(filter(fdr_old, fam_self=="bond_1y")), nrow(filter(fdr_keep, fam_self=="bond_1y"))))
cat(sprintf("    BH (lenient) bond survivors  bond_1y/pooled/univ : %d/%d/%d -> %d/%d/%d\n",
            bh_before["b1"], bh_before["pool"], bh_before["univ"],
            bh_after["b1"],  bh_after["pool"],  bh_after["univ"]))
cat(sprintf("    ANY BY survivor anywhere after drop-it1? %s   (verdict %s)\n",
            by_any_after, if (by_any_after) "CHANGED — investigate" else "unchanged: 0 survivors"))
cat(sprintf("    n_live_bond guard (.Rmd line 111): %d -> %d\n", nlive_before, nlive_after))

cat("\n[6] Equity-S2 Hansen SPA per (h x sample); p_consistent = headline\n")
print(as.data.frame(spa |> arrange(family, horizon, sample) |>
  select(family, horizon, sample, m_orig, m_used, n_common, n_dropped, T_spa, p_consistent, p_upper, reason)))
cat(sprintf("    SPA families with p_consistent < 0.05: %d / %d evaluable\n",
            nrow(spa_hits), sum(spa$reason == "")))
cat(sprintf("    equity-S2 evaluable cells: 432 -> 144 (shrinkage only); P10 accepts: 21 -> 0\n"))

cat("\n[write] artifacts:\n")
for (f in sprintf("data/audit/%s_%s.parquet",
                  c("p10_rescored","stage4_fdr","stage4_equity_inference","stage4_spa"), OUT_TAG))
  cat("   ", f, "\n")
cat("[run_drop_it1_rescore complete]\n")

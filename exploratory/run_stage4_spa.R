# working_paper/run_stage4_spa.R
# Stage 4 / S4.2 (P11 / P-D) — overlap-robust CW + Hodrick-1B + Hansen (2005)
# SPA on the equity-S2 sub-family, using the RETAINED pointwise forecasts
# (data/audit/equity_s2_forecasts_<date>.parquet). RUN FROM REPO ROOT.
#
# Adapts working_paper/run_residual_inference.R (WS4). Differences:
#   * family is the full equity-S2 evaluable set (432 cells), not the post-2008
#     residual region; keys carry no target/maturity (equity-only).
#   * SPA is per (horizon x SAMPLE) family + an expanding-only robustness variant
#     per family (PD4 / WS4 decisions 5-6): SPA p-values are driven by common-
#     as-of-axis length far more than family membership, so each (h, sample) is
#     a clean homogeneous-window family; expanding-only gives the cleanest axis.
#   * heavy-shrinkage cells can collapse to the benchmark (y_hat==y_bench =>
#     d==0 constant column); such benchmark-identical cells cannot be the best
#     model, so they are DROPPED from the SPA matrix with an explicit n_dropped
#     count (tracked, never silent — [[feedback_no_silent_fails]]).
#
# h=12 equity is infeasible under U2 (expanding-only-below-floor) and was not
# retained -> only h in {1,3} appear. The lone WS4 SPA survivor lived at h=12,
# so the equity SPA is now structurally CONFIRMATORY of the near-null.
#
# Gates: G1 recompute (r2_oos, cw, sharpe) from forecasts, match the stored p10
#        equity-S2 cells to < 1e-8; no-NA-inference stop guard; common axis >= 10.
suppressMessages({library(arrow); library(dplyr); library(tidyr); library(purrr)})
source("R/predictive_metrics.R");      source("R/predictive_regression.R")
source("R/target_relative_sharpe.R");  source("R/overlap_robust_inference.R")
source("R/spa_test.R")
source("R/pit_query.R"); source("R/panel.R"); source("R/targets.R")

AS_OF      <- as.Date("2026-04-30")
RUNDATE_IN <- "2026-06-12"
KEY <- c("chapter","scenario","model_config","iteration","horizon",
         "window_kind","window_length","sample")
SPA_B <- 2000L

# ---- load ------------------------------------------------------------------
fc  <- read_parquet(sprintf("data/audit/equity_s2_forecasts_%s.parquet", RUNDATE_IN))
p10 <- read_parquet(sprintf("data/audit/p10_rescored_%s.parquet", RUNDATE_IN)) |>
  filter(target == "equity", grepl("S2", scenario), !is.na(accept_p10))
n_cells <- nrow(distinct(fc, across(all_of(KEY))))
cat(sprintf("[load] %d forecast rows / %d cells; %d equity-S2 evaluable in p10\n",
            nrow(fc), n_cells, nrow(p10)))
stopifnot(n_cells == nrow(p10))                          # 432 == 432

# ---- one-period equity excess return (Hodrick 1B innovation) ---------------
panel  <- build_panel(AS_OF, c("SP500TR","DGS3MO"), "M", "forward")
r1_tbl <- tibble(as_of_date = panel$reference_date,
                 r1 = construct_equity_excess_return(panel, h = 1L, sp_col = "SP500TR"))

# ---- per-cell recompute + overlap-robust inference -------------------------
fc <- fc |> arrange(across(all_of(KEY)), as_of_date)
inf <- fc |> group_by(across(all_of(KEY))) |> group_split() |> map_dfr(function(d) {
  h  <- d$horizon[1L]
  r2 <- r2_oos(d$y_hat, d$y_realized, d$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(d$y_hat, d$y_realized, d$y_bench),
                        L = floor(1.5 * h))$stat
  sh  <- signal_sharpe(d$y_hat, d$y_realized, h = h)
  cro <- cw_overlap_robust(d$y_hat, d$y_realized, d$y_bench, h)
  dr  <- left_join(d, r1_tbl, by = "as_of_date")
  if (anyNA(dr$r1)) stop("run_stage4_spa: r1 join left NA for an as-of date (alignment bug)")
  hod <- hodrick_1b_mz(dr$y_realized, dr$y_hat, dr$r1, h)
  nw  <- nw_hac_mz(dr$y_realized, dr$y_hat, h)
  d[1L, KEY] |> bind_cols(tibble(
    n_oos = nrow(d), r2_oos = r2, cw = cw, sharpe = sh,
    n_eff = cro$n_eff, cw_overlap = cro$cw_overlap, p_normal = cro$p_normal,
    p_t_eff = cro$p_t_eff, beta_mz = hod$beta, p_hodrick = hod$p_hodrick, p_nwhac = nw$p_nwhac))
})
bad <- inf |> filter(is.na(p_t_eff) | is.na(p_hodrick) | is.na(p_nwhac) | is.na(p_normal))
if (nrow(bad)) stop(sprintf("%d evaluable cell(s) produced NA inference; investigate before scoring", nrow(bad)))

# ---- G1: exact reproduction vs stored p10 cells ----------------------------
chk <- inf |> select(all_of(KEY), r2_oos, cw, sharpe) |>
  inner_join(p10 |> select(all_of(KEY), r2_s = r2_oos, cw_s = cw, sh_s = sharpe), by = KEY)
if (nrow(chk) != nrow(inf))
  stop(sprintf("G1: %d of %d re-run cells failed to match a p10 row", nrow(inf) - nrow(chk), nrow(inf)))
g1 <- chk |> summarise(dr2 = max(abs(r2_oos - r2_s)), dcw = max(abs(cw - cw_s)), dsh = max(abs(sharpe - sh_s)))
cat(sprintf("[G1] %d cells matched; max|Δr²|=%.2e max|Δcw|=%.2e max|Δsharpe|=%.2e\n",
            nrow(chk), g1$dr2, g1$dcw, g1$dsh))
stopifnot(g1$dr2 < 1e-8, g1$dcw < 1e-8, g1$dsh < 1e-8)

# ---- attach P10 accept + FDR across the equity-S2 family -------------------
# (overlap-robust L=h-1 p_t_eff — the PD4 variant available only on retained
#  forecasts — plus the in-sample Hodrick p; complements the S4.1 scalar FDR.)
inf <- inf |> left_join(p10 |> select(all_of(KEY), threshold, accept_p10), by = KEY) |>
  mutate(q_bh_oos = bh_fdr(p_t_eff),  q_by_oos = by_fdr(p_t_eff),
         q_bh_is  = bh_fdr(p_hodrick), q_by_is  = by_fdr(p_hodrick))

# ---- SPA per (horizon x sample); + expanding-only robustness variant -------
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
    return(row(length(per), 0L, length(common), 0L, list(),
               sprintf("common axis < 10 (%d)", length(common))))
  mat   <- vapply(per, function(p) p$d[match(common, p$as_of_date)], numeric(length(common)))
  sds   <- apply(mat, 2L, stats::sd)
  ndrop <- sum(sds == 0)                                  # benchmark-identical cells
  mat   <- mat[, sds > 0, drop = FALSE]
  if (ncol(mat) < 1L)
    return(row(length(per), 0L, length(common), ndrop, list(), "all columns benchmark-identical"))
  o <- spa_test(mat, q = 1 / h, B = SPA_B, seed = 11L)
  row(length(per), o$m, length(common), ndrop, o, "")
}
`%||%` <- function(a, b) if (length(a)) a else b
combos <- fc |> distinct(horizon, sample)
spa <- bind_rows(
  pmap_dfr(combos, function(horizon, sample) spa_family(horizon, sample, fc, "all_windows")),
  pmap_dfr(combos, function(horizon, sample)
    spa_family(horizon, sample, filter(fc, window_kind == "expanding"), "expanding_only")))

# ---- persist ---------------------------------------------------------------
write_parquet(inf, file.path("data/audit", sprintf("stage4_equity_inference_%s.parquet", Sys.Date())))
write_parquet(spa, file.path("data/audit", sprintf("stage4_spa_%s.parquet", Sys.Date())))

# ---- headline --------------------------------------------------------------
acc <- inf |> filter(accept_p10)
surv <- function(df) tibble(
  cells = nrow(df), cw_overlap_p_lt05 = sum(df$p_t_eff < 0.05, na.rm = TRUE),
  oos_q_bh_lt05 = sum(df$q_bh_oos < 0.05, na.rm = TRUE),
  oos_q_by_lt05 = sum(df$q_by_oos < 0.05, na.rm = TRUE),
  hodrick_p_lt05 = sum(df$p_hodrick < 0.05, na.rm = TRUE),
  is_q_bh_lt05 = sum(df$q_bh_is < 0.05, na.rm = TRUE),
  is_q_by_lt05 = sum(df$q_by_is < 0.05, na.rm = TRUE))
cat(sprintf("\n=== equity-S2 family: %d evaluable, %d P10-accepted ===\n", nrow(inf), nrow(acc)))
cat("\n-- survivors among P10-ACCEPTED --\n");  print(as.data.frame(surv(acc)))
cat("\n-- survivors among ALL evaluable --\n"); print(as.data.frame(surv(inf)))
cat("\n-- Hansen SPA per (horizon x sample); p_consistent = headline --\n")
print(as.data.frame(spa |> arrange(family, horizon, sample) |>
  select(family, horizon, sample, m_orig, m_used, n_common, n_dropped,
         T_spa, best_z, p_consistent, p_upper, reason)))
spa_hits <- spa |> filter(!is.na(p_consistent), p_consistent < 0.05)
cat(sprintf("\n================= HEADLINE =================\nequity-S2 SPA families with p_consistent < 0.05: %d / %d evaluable families\n",
            nrow(spa_hits), sum(spa$reason == "")))
if (nrow(spa_hits)) print(as.data.frame(spa_hits |> select(family, horizon, sample, T_spa, p_consistent)))
cat("[run_stage4_spa complete]\n")

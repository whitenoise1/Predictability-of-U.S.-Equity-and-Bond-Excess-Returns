# working_paper/run_residual_inference.R
# WS4 Stage B (P8) — overlap-robust inference + FDR + SPA on the residual cell
# family, using the RETAINED pointwise forecasts (Stage A:
# data/audit/residual_forecasts_<date>.parquet). RUN FROM REPO ROOT.
#
# Pipeline:
#   G1  exact reproduction: recompute (r2_oos, cw, sharpe) from the retained
#       forecasts; must equal the stored p10_rescored values to < 1e-8.
#   (A) OOS Clark-West, overlap-robust  -> p_t_eff (t(n_eff) reference).
#   (B) in-sample Mincer-Zarnowitz      -> Hodrick(1992) 1B p_hodrick (+ NW-HAC).
#   FDR Benjamini-Hochberg + Benjamini-Yekutieli across the evaluable family.
#   SPA Hansen (2005), per horizon, on the common as-of axis (q = 1/h).
# Writes data/audit/residual_inference_<date>.parquet (+ _spa) and prints the
# headline: of the P10-accepted residual cells, how many survive each filter.
suppressMessages({library(arrow); library(dplyr); library(tidyr); library(purrr)})
source("R/predictive_metrics.R");      source("R/predictive_regression.R")
source("R/target_relative_sharpe.R");  source("R/overlap_robust_inference.R")
source("R/spa_test.R")
source("R/pit_query.R"); source("R/panel.R"); source("R/targets.R")

AS_OF   <- as.Date("2026-04-30")
RUNDATE <- "2026-06-01"
KEY <- c("chapter","scenario","model_config","iteration","horizon",
         "window_kind","window_length","sample")
SPA_B <- 2000L

# ---- load ------------------------------------------------------------------
fc  <- read_parquet(sprintf("data/audit/residual_forecasts_%s.parquet", RUNDATE))
p10 <- read_parquet("data/audit/p10_rescored_2026-06-01.parquet") |>
  filter(target=="equity", grepl("S2", scenario), sample=="post_2008",
         horizon %in% c(3L,12L), model_config %in% c("C-OLS","C-COMB","C-CT"))
cat(sprintf("[load] %d forecast rows; %d residual-region cells in p10 (%d evaluable)\n",
            nrow(fc), nrow(p10), sum(!is.na(p10$accept_p10))))

# ---- one-period equity excess return (Hodrick 1B innovation) ---------------
panel  <- build_panel(AS_OF, c("SP500TR","DGS3MO"), "M", "forward")
r1_tbl <- tibble(as_of_date = panel$reference_date,
                 r1 = construct_equity_excess_return(panel, h = 1L, sp_col = "SP500TR"))

# ---- per-cell recompute + inference ----------------------------------------
fc <- fc |> arrange(across(all_of(KEY)), as_of_date)
cell_rows <- fc |> group_by(across(all_of(KEY))) |> group_split()
inf <- map_dfr(cell_rows, function(d) {
  h <- d$horizon[1L]
  r2 <- r2_oos(d$y_hat, d$y_realized, d$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(d$y_hat, d$y_realized, d$y_bench),
                        L = floor(1.5*h))$stat
  sh <- signal_sharpe(d$y_hat, d$y_realized, h = h)
  cro <- cw_overlap_robust(d$y_hat, d$y_realized, d$y_bench, h)
  dr  <- left_join(d, r1_tbl, by = "as_of_date")
  if (anyNA(dr$r1)) stop("run_residual_inference: r1 join left NA for an as-of date (alignment bug)")
  hod <- hodrick_1b_mz(dr$y_realized, dr$y_hat, dr$r1, h)   # fail loud on degeneracy
  nw  <- nw_hac_mz(dr$y_realized, dr$y_hat, h)
  d[1L, KEY] |>
    bind_cols(tibble(n_oos = nrow(d), r2_oos = r2, cw = cw, sharpe = sh,
                     n_eff = cro$n_eff, cw_overlap = cro$cw_overlap,
                     p_normal = cro$p_normal, p_t_eff = cro$p_t_eff,
                     beta_mz = hod$beta, p_hodrick = hod$p_hodrick, p_nwhac = nw$p_nwhac))
})
# no-silent-fail guard: every evaluable cell must yield finite inference — an NA
# p-value would otherwise shrink the FDR denominators / survivor tallies untracked.
bad <- inf |> filter(is.na(p_t_eff) | is.na(p_hodrick) | is.na(p_nwhac) | is.na(p_normal))
if (nrow(bad)) stop(sprintf("%d evaluable cell(s) produced NA inference (degenerate LRV/fit); investigate before scoring", nrow(bad)))

# ---- G1: exact reproduction vs stored p10 cells ----------------------------
chk <- inf |> select(all_of(KEY), r2_oos, cw, sharpe) |>
  inner_join(p10 |> select(all_of(KEY), r2_s = r2_oos, cw_s = cw, sh_s = sharpe), by = KEY)
if (nrow(chk) != nrow(inf))
  stop(sprintf("G1: %d of %d re-run cells failed to match a p10 row", nrow(inf)-nrow(chk), nrow(inf)))
g1 <- chk |> summarise(dr2 = max(abs(r2_oos-r2_s)), dcw = max(abs(cw-cw_s)), dsh = max(abs(sharpe-sh_s)))
cat(sprintf("[G1] %d cells matched; max|Δr²|=%.2e max|Δcw|=%.2e max|Δsharpe|=%.2e\n",
            nrow(chk), g1$dr2, g1$dcw, g1$dsh))
stopifnot(g1$dr2 < 1e-8, g1$dcw < 1e-8, g1$dsh < 1e-8)

# ---- attach P10 threshold + accept; FDR across the evaluable family --------
inf <- inf |> left_join(p10 |> select(all_of(KEY), threshold, accept_p10), by = KEY)
inf <- inf |> mutate(
  q_bh_oos = bh_fdr(p_t_eff),  q_by_oos = by_fdr(p_t_eff),
  q_bh_is  = bh_fdr(p_hodrick), q_by_is = by_fdr(p_hodrick))

# ---- SPA per horizon (common as-of axis; q = 1/h) --------------------------
# Two families per horizon: "all_windows" (locked design) and "expanding_only"
# — a robustness variant. Expanding cells share a homogeneous t_eval, so their
# common as-of axis is longer/cleaner than the all-windows intersection (which
# the late-starting rolling/iter_2 cells shrink), giving the h=12 SPA more room.
spa_one <- function(h, fc_sub, family) {
  per <- fc_sub |> filter(horizon == h) |> group_by(across(all_of(KEY))) |>
    group_map(~ tibble(as_of_date = .x$as_of_date,
                       d = (.x$y_realized-.x$y_bench)^2 - (.x$y_realized-.x$y_hat)^2))
  common <- sort(as.Date(reduce(map(per, ~ as.character(.x$as_of_date)), intersect)))
  if (length(per) < 1L || length(common) < 10L)  # record the skip explicitly
    return(tibble(family = family, horizon = h, m = length(per), n = length(common),
                  n_common = length(common), q = 1/h, B = SPA_B, T_spa = NA_real_,
                  best_z = NA_real_, p_lower = NA_real_, p_consistent = NA_real_,
                  p_upper = NA_real_, reason = sprintf("common axis < 10 (%d)", length(common))))
  mat <- vapply(per, function(p) p$d[match(common, p$as_of_date)], numeric(length(common)))
  out <- spa_test(mat, q = 1/h, B = SPA_B, seed = 11L)
  out$family <- family; out$horizon <- h; out$n_common <- length(common); out$reason <- ""; out
}
spa <- bind_rows(
  map_dfr(c(3L,12L), ~ spa_one(.x, fc,                                   "all_windows")),
  map_dfr(c(3L,12L), ~ spa_one(.x, filter(fc, window_kind=="expanding"), "expanding_only")))

# ---- persist ---------------------------------------------------------------
write_parquet(inf, sprintf("data/audit/residual_inference_%s.parquet", RUNDATE))
if (nrow(spa)) write_parquet(spa, sprintf("data/audit/residual_spa_%s.parquet", RUNDATE))

# ---- headline summary ------------------------------------------------------
acc <- inf |> filter(accept_p10)
cat(sprintf("\n=== residual family: %d evaluable cells, %d P10-accepted ===\n",
            nrow(inf), nrow(acc)))
surv <- function(df) tibble(
  cells              = nrow(df),
  cw_overlap_p_lt05  = sum(df$p_t_eff   < 0.05, na.rm=TRUE),
  oos_q_bh_lt05      = sum(df$q_bh_oos  < 0.05, na.rm=TRUE),
  oos_q_by_lt05      = sum(df$q_by_oos  < 0.05, na.rm=TRUE),
  hodrick_p_lt05     = sum(df$p_hodrick < 0.05, na.rm=TRUE),
  is_q_bh_lt05       = sum(df$q_bh_is   < 0.05, na.rm=TRUE),
  is_q_by_lt05       = sum(df$q_by_is   < 0.05, na.rm=TRUE))
cat("\n-- survivors among P10-ACCEPTED cells --\n");   print(as.data.frame(surv(acc)))
cat("\n-- survivors among ALL evaluable cells --\n");  print(as.data.frame(surv(inf)))
cat("\n-- by horizon (accepted) --\n")
print(acc |> group_by(horizon) |>
  summarise(n=n(), n_eff_med=median(n_eff), cw_med=round(median(cw),2),
            p_norm_med=round(median(p_normal),4), p_teff_med=round(median(p_t_eff),4),
            p_hod_med=round(median(p_hodrick),4), .groups="drop") |> as.data.frame())
cat("\n-- Hansen SPA per horizon (all_windows = locked; expanding_only = robustness) --\n")
if (nrow(spa)) print(as.data.frame(spa |> arrange(family, horizon) |>
  select(family, horizon, m, n_common, T_spa, best_z, p_lower, p_consistent, p_upper)))
cat(sprintf("\n[done] wrote residual_inference_%s.parquet\n", RUNDATE))

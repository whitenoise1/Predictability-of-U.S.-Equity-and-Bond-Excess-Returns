# working_paper/rerun_bond_forecasts.R
#
# Forecast-RETAINING re-run of the S2_bond SHRINKAGE class for ALL FIVE maturities
# (1/2/3/5/10Y). Levels pool (chapter ch2) only. Bonds never retained pointwise
# forecasts (only equity did); this closes that gap so the Hansen SPA test becomes
# possible per maturity (not just 1Y). The shrinkage fit is unchanged by the dormant
# C-PCR refactor, so this reproduces the committed drop-it1 bond cells (all 360;
# G1 enforced in the merge step).
#
# The harness recomputes inner PIT CV at every as-of date (~6 min for an h=1
# expanding full-sample cell), so the work is SLICED for parallelism: one process
# per (config x sample x maturity), each writing a slice file. Drive with env CFG,
# SMP and MAT.
#   CFG in {C-ENET,C-RIDGE,C-PLS}  SMP in {full,pre_2008,post_2008,post_2016}
#   MAT in {bond_1y,bond_2y,bond_3y,bond_5y,bond_10y}
# Then merge + G1 with merge_bond_forecasts.R. RUN FROM REPO ROOT.
suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R")
source("R/u6_pool_prep.R")

CFG <- Sys.getenv("CFG"); SMP <- Sys.getenv("SMP"); MAT <- Sys.getenv("MAT")
if (!nzchar(CFG) || !nzchar(SMP) || !nzchar(MAT)) stop("set env CFG (C-ENET/C-RIDGE/C-PLS), SMP (full/pre_2008/post_2008/post_2016) and MAT (bond_1y/2y/3y/5y/10y)")
stopifnot(CFG %in% c("C-ENET","C-RIDGE","C-PLS"),
          SMP %in% c("full","pre_2008","post_2008","post_2016"),
          MAT %in% c("bond_1y","bond_2y","bond_3y","bond_5y","bond_10y"))
AS_OF        <- as.Date("2026-04-30")
W_Z          <- 60L
ROLL_LENGTHS <- c(12L, 24L, 36L, 48L, 60L, 84L)
HORIZONS_RUN <- c(1L, 3L)                      # h=12 has no evaluable bond_1y cell
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
MIN_OOS_ROWS <- 24L
T_START_S2_CUTOFF <- as.Date("1990-01-31")
# Canonical 5-maturity grid (p1_chapter_2.R:47-51); slice runs one MAT at a time.
MATURITY_GRID <- list(
  bond_1y  = list(y_col = "DGS1",  duration = 0.97),
  bond_2y  = list(y_col = "DGS2",  duration = 1.93),
  bond_3y  = list(y_col = "DGS3",  duration = 2.85),
  bond_5y  = list(y_col = "DGS5",  duration = 4.50),
  bond_10y = list(y_col = "DGS10", duration = 8.00))[MAT]

# ---- setup (verbatim p1_chapter_2.R sections 1-6) --------------------------
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
targets_bd <- bond_target_grid(panel_levels, c(1L,3L,12L), MATURITY_GRID)
s1bf <- bond_theory_pool_ff1989(panel_levels); panel_levels <- s1bf$panel
s1eq <- equity_theory_pool_wg_crr(panel_levels)
dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
s2f <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
         forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
s2d <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
         forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
all_ratio <- unique(c(s1eq$ratio_features, s2f$ratio_features))
all_level <- unique(c(s1bf$level_features, s1eq$level_features, s2f$level_features))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
  apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)
S2_POOL <- c(s2f$ratio_features, s2f$level_features)
S2_DEEP <- live_cols(panel_levels, c(s2d$ratio_features, s2d$level_features))
fnn <- function(x, d) { i <- which(!is.na(x)); if (!length(i)) as.Date(NA) else d[i[1L]] }
t_start <- max(as.Date(sapply(S2_DEEP, function(c) as.character(fnn(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval0 <- seq(t_start, by = "120 months", length.out = 2L)[2L]
sb <- switch(SMP,
  full      = list(t_eval = t_eval0, t_end = AS_OF),
  pre_2008  = list(t_eval = t_eval0, t_end = SPLIT_2008 - 1L),
  post_2008 = list(t_eval = max(t_eval0, SPLIT_2008), t_end = AS_OF),
  post_2016 = list(t_eval = max(t_eval0, SPLIT_2016), t_end = AS_OF))
panel_cell <- panel_xfm |> select(reference_date, all_of(S2_POOL))

# ---- sweep this (config, sample) slice over h in {1,3} x 7 windows ----------
WINDOWS <- c(list(c("expanding", NA)), lapply(ROLL_LENGTHS, function(w) c("rolling", w)))
fc_store <- list(); t0 <- Sys.time()
if (sb$t_eval < sb$t_end) for (h in HORIZONS_RUN) {
  for (win in WINDOWS) {
    wk <- win[[1L]]; wl <- if (wk == "rolling") as.integer(win[[2L]]) else NA_integer_
    out <- tryCatch(walk_forward_fs(y = targets_bd[[MAT]][[as.character(h)]],
        panel = panel_cell, horizon = h, window_kind = wk,
        window_length = if (wk == "rolling") wl else NULL,
        t_start = t_start, t_eval_start = sb$t_eval, t_end = sb$t_end,
        model_config = CFG, iteration_kind = "iter_shrink_full",
        features = S2_POOL, target = "bond",
        enet_pool = if (CFG == "C-ENET") S2_POOL else character(0),
        shrink_pool = if (CFG %in% c("C-RIDGE","C-PLS")) S2_POOL else character(0)),
      error = function(e) list(forecasts = tibble()))
    fl <- apply_feasibility_floor(out$forecasts, h, CFG); wf <- fl$wf
    if (nrow(wf) < MIN_OOS_ROWS) next            # not evaluable -> not retained
    fc_store[[length(fc_store)+1L]] <- tibble(chapter = "ch2", scenario = "S2_bond",
      target = "bond", maturity = MAT, model_config = CFG, iteration = "iter_shrink_full",
      horizon = h, window_kind = wk, window_length = wl, sample = SMP,
      as_of_date = wf$as_of_date, y_hat = wf$y_hat, y_realized = wf$y_realized, y_bench = wf$y_bench)
  }
}
fc_all <- if (length(fc_store)) bind_rows(fc_store) else tibble()
out <- file.path("data/audit", sprintf("_bondfc_slice_%s_%s_%s.parquet", CFG, SMP, MAT))
write_parquet(fc_all, out)
cat(sprintf("[slice %s/%s/%s] %d cells, %d rows, %.0f s -> %s\n", CFG, SMP, MAT,
            if (nrow(fc_all)) nrow(distinct(fc_all, horizon, window_kind, window_length)) else 0L,
            nrow(fc_all), as.numeric(Sys.time()-t0, units="secs"), out))

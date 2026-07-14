# working_paper/rerun_s1_forecasts.R
#
# Forecast-RETAINING re-run of the SCENARIO-1 (classic-literature) pools, so the
# manuscript can plot predicted-vs-realized for the theory replications (Figure 5).
# S1 never retained pointwise forecasts (only S2 equity + bond_1y did).
#
# Scope: the four S1 scenarios, CLASSIC OLS-family configs only --
#   S1_bond_FF1989      (C-OLS, C-COMB)         x 5 maturities
#   S1_bond_CP          (C-CP)                  x 5 maturities
#   S1_equity_WGCRR     (C-OLS, C-COMB, C-CT)
#   S1_equity_WGCRR_aug (C-OLS, C-COMB, C-CT)   [WG+CRR + CPI/INDPRO surprises]
# C-ENET is excluded: it is not the "classic OLS" the figure is about, and in S1
# the original driver routes it onto the wide S2 pool (a quirk we do not surface).
#
# Faithfulness: the setup + per-cell call mirror scripts/p1_chapter_2.R (FF/CP/
# WGCRR) and scripts/p1_chapter_3.R (WGCRR_aug) verbatim, so the reproduced cells
# reproduce the committed S1 metrics. A G1 check asserts r2_oos and cw match the
# committed stage-4 FDR cells (<1e-6 / <1e-4) before anything is written.
# RUN FROM REPO ROOT.
suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R")
source("R/surprise.R"); source("R/u6_pool_prep.R")

AS_OF        <- as.Date("2026-04-30")
W_Z          <- 60L
ROLL_LENGTHS <- c(36L, 48L, 60L, 84L)          # rolling-12/24 excluded by the n_eff floor (not in committed grid)
HORIZONS     <- c(1L, 3L, 12L)                  # bond grid; equity caps at h<=3 (committed set has no equity h=12)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
MIN_OOS_ROWS <- 24L
MATURITY_GRID <- list(
  bond_1y  = list(y_col = "DGS1",  duration = 0.97),
  bond_2y  = list(y_col = "DGS2",  duration = 1.93),
  bond_3y  = list(y_col = "DGS3",  duration = 2.85),
  bond_5y  = list(y_col = "DGS5",  duration = 4.50),
  bond_10y = list(y_col = "DGS10", duration = 8.00))

# ---- panel + targets (verbatim p1_chapter_2.R sections 1-2) ------------------
cat("[1] Building panel ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR", u6_extra_factor_ids(meta)))
panel0 <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
# Equity target = SP500TR month-end TOTAL return (matches the committed corrected
# pipeline rerun_equity_corrected.R), NOT the default SHILLER_PRICE within-month
# average (the temporal-averaging artifact). S1 equity predictors exclude the raw
# price, so only the target source matters here; bonds are unaffected.
targets_eq <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel0, h = h, sp_col = "SP500TR")), as.character(HORIZONS))
targets_bd <- bond_target_grid(panel0, HORIZONS, MATURITY_GRID)

# ---- base pools (FF / CP / WGCRR) on the no-surprise panel -------------------
cat("[2] Building S1 pools ...\n")
s1ff   <- bond_theory_pool_ff1989(panel0); panel_b <- s1ff$panel
s1eq   <- equity_theory_pool_wg_crr(panel_b)
ff_lvl <- s1ff$level_features; cp_lvl <- CP_FACTOR_COLS
eq_rat <- s1eq$ratio_features;  eq_lvl <- s1eq$level_features
xfm_b <- panel_b |>
  apply_transform(yoy_log_diff,   factor_ids = eq_rat) |>
  apply_transform(rolling_zscore, factor_ids = eq_rat, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = unique(c(ff_lvl, cp_lvl, eq_lvl)),
                  window = W_Z, min_obs = W_Z)

# ---- augmented WG+CRR on the surprise panel (verbatim p1_chapter_3.R 1.5/3/4) -
surp <- build_surprise_panel(AS_OF, meta = meta)
panel_a <- panel0 |> left_join(surp, by = "reference_date") |>
  mutate(CPI_SURPRISE    = ifelse(!is.na(CPI_SURPRISE)    & abs(CPI_SURPRISE)    > 50, NA_real_, CPI_SURPRISE),
         INDPRO_SURPRISE = ifelse(!is.na(INDPRO_SURPRISE) & abs(INDPRO_SURPRISE) > 15, NA_real_, INDPRO_SURPRISE))
panel_a <- bond_theory_pool_ff1989(panel_a)$panel
s1eqA   <- equity_theory_pool_wg_crr_aug(panel_a)
eqA_rat <- s1eqA$ratio_features; eqA_lvl <- s1eqA$level_features
xfm_a <- panel_a |>
  apply_transform(yoy_log_diff,   factor_ids = eqA_rat) |>
  apply_transform(rolling_zscore, factor_ids = eqA_rat, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = eqA_lvl, window = W_Z, min_obs = W_Z)

# ---- helpers (verbatim p1_chapter_2.R sections 5/7/8) ------------------------
fnn <- function(x, d) { i <- which(!is.na(x)); if (!length(i)) as.Date(NA) else d[i[1L]] }
t_start_for <- function(xfm, cols)
  max(as.Date(sapply(cols, function(c) as.character(fnn(xfm[[c]], xfm$reference_date)))))
t_eval_from <- function(t_start) seq(t_start, by = "120 months", length.out = 2L)[2L]
sample_bounds <- function(lbl, t_start) { td <- t_eval_from(t_start); switch(lbl,
  full = list(t_eval = td, t_end = AS_OF),
  pre_2008  = list(t_eval = td, t_end = SPLIT_2008 - 1L),
  post_2008 = list(t_eval = max(td, SPLIT_2008), t_end = AS_OF),
  post_2016 = list(t_eval = max(td, SPLIT_2016), t_end = AS_OF)) }
SAMPLES <- c("full","pre_2008","post_2008","post_2016")
WINDOWS <- c(list(c("expanding", NA)), lapply(ROLL_LENGTHS, function(w) c("rolling", w)))

# forecast-retaining cell runner: same walk_forward_fs call as run_cell, returns wf + metrics
run_cell_fc <- function(y_h, panel_cell, h, wk, wl, t_start, t_eval, t_end, cfg, tg) {
  L_nw <- floor(1.5 * h)
  out <- tryCatch(walk_forward_fs(y = y_h, panel = panel_cell, horizon = h,
      window_kind = wk, window_length = if (wk == "rolling") wl else NULL,
      t_start = t_start, t_eval_start = t_eval, t_end = t_end,
      model_config = cfg, iteration_kind = NA_character_,
      features = setdiff(names(panel_cell), "reference_date"), target = tg,
      enet_pool = character(0), shrink_pool = character(0),
      pcr_pool = character(0), comb_pool = character(0)),
    error = function(e) list(forecasts = tibble()))
  fl <- apply_feasibility_floor(out$forecasts, h, cfg); wf <- fl$wf
  if (nrow(wf) < MIN_OOS_ROWS) return(NULL)
  list(wf = wf,
       r2 = r2_oos(wf$y_hat, wf$y_realized, wf$y_bench),
       cw = clark_west_stat(clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench), L = L_nw)$stat,
       sh = signal_sharpe(wf$y_hat, wf$y_realized, h = h))
}

# ---- spec list: scenario x configs x (maturity) ------------------------------
specs <- list()
for (m in names(MATURITY_GRID)) for (cfg in c("C-OLS","C-COMB"))
  specs[[length(specs)+1L]] <- list(sc="S1_bond_FF1989", tg="bond", mat=m, cfg=cfg,
    xfm=xfm_b, cols=ff_lvl, yL=targets_bd[[m]])
for (m in names(MATURITY_GRID))
  specs[[length(specs)+1L]] <- list(sc="S1_bond_CP", tg="bond", mat=m, cfg="C-CP",
    xfm=xfm_b, cols=cp_lvl, yL=targets_bd[[m]])
for (cfg in c("C-OLS","C-COMB","C-CT"))
  specs[[length(specs)+1L]] <- list(sc="S1_equity_WGCRR", tg="equity", mat=NA_character_, cfg=cfg,
    xfm=xfm_b, cols=c(eq_rat, eq_lvl), yL=targets_eq)
for (cfg in c("C-OLS","C-COMB","C-CT"))
  specs[[length(specs)+1L]] <- list(sc="S1_equity_WGCRR_aug", tg="equity", mat=NA_character_, cfg=cfg,
    xfm=xfm_a, cols=c(eqA_rat, eqA_lvl), yL=targets_eq)

# ---- sweep -------------------------------------------------------------------
cat(sprintf("[3] Sweeping %d specs x %d samples x %d windows x 3 horizons ...\n",
            length(specs), length(SAMPLES), length(WINDOWS)))
fc_store <- list(); mt_store <- list(); t0 <- Sys.time()
for (sp in specs) {
  panel_cell <- sp$xfm |> select(reference_date, all_of(sp$cols))
  t_start <- t_start_for(sp$xfm, sp$cols)
  for (h in (if (sp$tg == "bond") HORIZONS else setdiff(HORIZONS, 12L))) {
    yh <- sp$yL[[as.character(h)]]
    for (lbl in SAMPLES) {
      sb <- sample_bounds(lbl, t_start); if (sb$t_eval >= sb$t_end) next
      for (win in WINDOWS) {
        wk <- win[[1L]]; wl <- if (wk == "rolling") as.integer(win[[2L]]) else NA_integer_
        r <- run_cell_fc(yh, panel_cell, h, wk, wl, t_start, sb$t_eval, sb$t_end, sp$cfg, sp$tg)
        if (is.null(r)) next
        key <- tibble(chapter = if (sp$sc == "S1_equity_WGCRR_aug") "ch3" else "ch2",
          scenario = sp$sc, target = sp$tg, maturity = sp$mat, model_config = sp$cfg,
          iteration = NA_character_, horizon = h, window_kind = wk, window_length = wl, sample = lbl)
        fc_store[[length(fc_store)+1L]] <- bind_cols(key[rep(1, nrow(r$wf)), ], tibble(
          as_of_date = r$wf$as_of_date, y_hat = r$wf$y_hat,
          y_realized = r$wf$y_realized, y_bench = r$wf$y_bench))
        mt_store[[length(mt_store)+1L]] <- bind_cols(key, tibble(r2_oos = r$r2, cw = r$cw, sharpe = r$sh))
      }
    }
  }
}
fc <- bind_rows(fc_store); mt <- bind_rows(mt_store)
cat(sprintf("    %d evaluable cells, %d forecast rows, %.0f s\n",
            nrow(mt), nrow(fc), as.numeric(Sys.time()-t0, units="secs")))

# ---- G1: reproduced metrics must match committed stage-4 FDR cells -----------
cat("[4] G1 check against committed FDR cells ...\n")
fdr <- read_parquet("data/audit/stage4_fdr_drop-it1_2026-06-17.parquet")
keys <- c("scenario","target","maturity","model_config","horizon","window_kind","window_length","sample")
# The committed FDR table carries one row per (cell x q-family universe) AND spans
# scenarios/iterations beyond this S1 sweep; scope to the cells mt references, then
# collapse the identical-metric q-family duplicates to the cell grain.
fdr_cells <- fdr |> select(all_of(keys), r2_oos, cw) |> semi_join(mt, by = keys) |> distinct()
stopifnot(nrow(fdr_cells) == nrow(distinct(fdr_cells, across(all_of(keys)))))  # one (r2,cw) per cell
cmp <- mt |> inner_join(fdr_cells, by = keys, suffix = c("", "_ref"))
stopifnot(nrow(cmp) == nrow(mt))                       # every reproduced cell found in committed set
d_r2 <- max(abs(cmp$r2_oos - cmp$r2_oos_ref))
d_cw <- max(abs(cmp$cw - cmp$cw_ref))
cat(sprintf("    matched %d/%d cells | max|dr2|=%.2e  max|dcw|=%.2e\n", nrow(cmp), nrow(mt), d_r2, d_cw))
if (d_r2 >= 1e-6 || d_cw >= 1e-4) {                     # self-diagnose any residual gap
  write_parquet(mt, "data/audit/_s1_mt_debug.parquet")
  off <- cmp |> mutate(dr2 = abs(r2_oos - r2_oos_ref)) |> filter(dr2 >= 1e-6)
  cat(sprintf("    G1 GAP: %d/%d cells off; by scenario:\n", nrow(off), nrow(cmp)))
  print(as.data.frame(dplyr::count(off, scenario, model_config, target)))
}
stopifnot(d_r2 < 1e-6, d_cw < 1e-4)

# ---- persist -----------------------------------------------------------------
out <- file.path("data/audit", sprintf("s1_forecasts_%s.parquet", Sys.Date()))
write_parquet(fc, out)
cat(sprintf("[5] Wrote %s (%d rows)\n", out, nrow(fc)))
print(fc |> count(scenario, model_config))

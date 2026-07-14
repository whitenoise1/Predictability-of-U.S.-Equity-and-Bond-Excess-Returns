# Phase-3 gate: one-cell old-vs-new diff. Compares the dropped S2 OLS-class path
# (C-OLS on the GS+PCA-selected deep cohort) against the new C-PCR (PCR on the
# full pool) on representative cells, BEFORE rewriting the driver / full re-run.
# Mirrors p1_chapter_2.R sections 1-6 for the panel/pools/targets.
suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R"); source("R/targets.R")
source("R/predictive_regression.R"); source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/u6_pool_prep.R")

AS_OF <- as.Date("2026-04-30"); W_Z <- 60L
SPLIT_2008 <- as.Date("2008-01-01"); T_START_S2_CUTOFF <- as.Date("1990-01-31")
PCA_GS_DEFAULTS <- list(n_groups=10L, min_var_explained=0.85, min_distance=0.5,
                        min_t_value=1.645, hac_lag=0L, verbose=FALSE)
MATURITY_GRID <- list(bond_1y=list(y_col="DGS1", duration=0.97),
                      bond_10y=list(y_col="DGS10", duration=8.00))

# --- panel / targets / pools (faithful to the driver) ---
meta <- read_parquet("data/metadata/factor_metadata.parquet")
ids <- meta |> filter(public_private %in% c("PUBLIC","MIXED"), !family %in% "FORECASTS",
                      !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors <- unique(c(ids, "SHILLER_PRICE", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(as_of_date=AS_OF, factor_ids=all_factors, frequency="M", fill="forward") |>
  cp_factor_columns()
targets_eq <- setNames(lapply(c(1L,3L,12L), function(h) construct_equity_excess_return(panel_levels, h=h)),
                       c("1","3","12"))
targets_bd <- bond_target_grid(panel_levels, c(1L,3L,12L), MATURITY_GRID)
s1_bond_ff <- bond_theory_pool_ff1989(panel_levels); panel_levels <- s1_bond_ff$panel
s1_eq <- equity_theory_pool_wg_crr(panel_levels)
dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff=NULL,
                                 forecast_token=U6_FORECAST_TOKEN, dashboard_features=dash$dashboard_features)
s2_deep <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff=T_START_S2_CUTOFF,
                                 forecast_token=U6_FORECAST_TOKEN, dashboard_features=dash$dashboard_features)
all_ratio <- unique(c(s1_eq$ratio_features, s2_full$ratio_features))
all_level <- unique(c(s1_bond_ff$level_features, CP_FACTOR_COLS, s1_eq$level_features, s2_full$level_features))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff, factor_ids=all_ratio) |>
  apply_transform(rolling_zscore, factor_ids=all_ratio, window=W_Z, min_obs=W_Z) |>
  apply_transform(rolling_zscore, factor_ids=all_level, window=W_Z, min_obs=W_Z)
S2_POOL_COLS <- c(s2_full$ratio_features, s2_full$level_features)
S2_DEEP_COLS <- live_cols(panel_levels, c(s2_deep$ratio_features, s2_deep$level_features))
cat(sprintf("pools: full=%d deep=%d\n", length(S2_POOL_COLS), length(S2_DEEP_COLS)))

first_non_na <- function(x, d){ i<-which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for <- function(cols) max(as.Date(sapply(cols, function(c)
  as.character(first_non_na(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval_from <- function(ts) seq(ts, by="120 months", length.out=2L)[2L]
panel_with_ref <- function(cols) panel_xfm |> select(reference_date, all_of(cols))

# --- scorecard (replicates run_cell's metric block) ---
score <- function(wf_forecasts, h) {
  fl <- apply_feasibility_floor(wf_forecasts, h, "x"); wf <- fl$wf
  if (nrow(wf) < 24L) return(tibble(n_oos=nrow(wf), r2_oos=NA, cw=NA, sharpe=NA, accept=NA))
  r2 <- r2_oos(wf$y_hat, wf$y_realized, wf$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench), L=floor(1.5*h))$stat
  sh <- signal_sharpe(wf$y_hat, wf$y_realized, h=h)
  tibble(n_oos=nrow(wf), r2_oos=r2, cw=cw, sharpe=sh, accept=accept_decision(r2,cw,sh))
}

run_old <- function(yL, h, t_start, t_eval, t_end, target, panel_cell) {
  # iter_1_full_sample GS+PCA selection on the deep cohort, then C-OLS.
  slc <- panel_xfm |> filter(reference_date>=T_START_S2_CUTOFF, reference_date<=AS_OF)
  yfull <- yL[["1"]][panel_xfm$reference_date>=T_START_S2_CUTOFF & panel_xfm$reference_date<=AS_OF]
  sel <- do.call(pca_gs_select, c(list(panel=slc, target_y=yfull, predictor_cols=S2_DEEP_COLS), PCA_GS_DEFAULTS))
  wf <- walk_forward_fs(y=yL[[as.character(h)]], panel=panel_cell, horizon=h,
                        window_kind="expanding", t_start=t_start, t_eval_start=t_eval, t_end=t_end,
                        model_config="C-OLS", iteration_kind="iter_1_full_sample",
                        features=sel$selected_features, target=target)$forecasts
  list(score=score(wf, h), n_feat=length(sel$selected_features))
}
run_new <- function(yL, h, t_start, t_eval, t_end, target, panel_cell) {
  wf <- walk_forward_fs(y=yL[[as.character(h)]], panel=panel_cell, horizon=h,
                        window_kind="expanding", t_start=t_start, t_eval_start=t_eval, t_end=t_end,
                        model_config="C-PCR", target=target, pcr_pool=S2_POOL_COLS)$forecasts
  list(score=score(wf, h), n_obj=if(nrow(wf)) round(mean(wf$n_features)) else NA)
}

fmtn <- function(x) if (length(x)==0 || is.na(x)) "NA" else sprintf("%.4f", x)
cells <- list(
  list(tag="bond_1y h=1 full",      yL=targets_bd[["bond_1y"]],  h=1L, target="bond",   smp="full"),
  list(tag="bond_10y h=1 full",     yL=targets_bd[["bond_10y"]], h=1L, target="bond",   smp="full"),
  list(tag="equity h=3 post2008",   yL=targets_eq,               h=3L, target="equity", smp="post_2008"))

for (cl in cells) {
  ts_old <- t_start_for(S2_DEEP_COLS); ts_new <- t_start_for(S2_DEEP_COLS)  # both anchor ~1990 on deep cohort
  te <- t_eval_from(ts_old)
  bnd <- if (cl$smp=="post_2008") list(te=max(te, SPLIT_2008), end=AS_OF) else list(te=te, end=AS_OF)
  pc_old <- panel_with_ref(S2_DEEP_COLS); pc_new <- panel_with_ref(S2_POOL_COLS)
  o <- run_old(cl$yL, cl$h, ts_old, bnd$te, bnd$end, cl$target, pc_old)
  n <- run_new(cl$yL, cl$h, ts_new, bnd$te, bnd$end, cl$target, pc_new)
  cat(sprintf("\n=== %s ===\n", cl$tag))
  cat(sprintf("  OLD C-OLS (deep cohort, %d GS+PCA feats): n_oos=%s r2=%s cw=%s sharpe=%s accept=%s\n",
      o$n_feat, o$score$n_oos, fmtn(o$score$r2_oos), fmtn(o$score$cw), fmtn(o$score$sharpe), o$score$accept))
  cat(sprintf("  NEW C-PCR (full pool, ~%s avail cols):     n_oos=%s r2=%s cw=%s sharpe=%s accept=%s\n",
      n$n_obj, n$score$n_oos, fmtn(n$score$r2_oos), fmtn(n$score$cw), fmtn(n$score$sharpe), n$score$accept))
}

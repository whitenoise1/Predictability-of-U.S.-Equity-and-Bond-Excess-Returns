# working_paper/rerun_equity_corrected.R
#
# Corrected-target re-run of the EQUITY legs of Chapter 2 (no surprises) and
# Chapter 3 (SPF surprises), mirroring scripts/p1_chapter_2.R + p1_chapter_3.R
# exactly for equity, with TWO deliberate corrections:
#   (1) equity target = SP500TR (month-end TOTAL return) instead of SHILLER_PRICE
#       (within-month average, price-only) — fixes the temporal-averaging artifact
#       (working_paper/VERIFICATION_01_equity_headline.md).
#   (2) SP500TR and SHILLER_PRICE excluded from the Scenario-2 predictor pool
#       (the index's own level must not predict its own return). Original Ch2
#       silently included SHILLER_PRICE among the 49 candidates.
# Everything else (configs, iterations, FS, window/sample grid, acceptance rule)
# is identical to Ch2/Ch3. Canonical scripts + artifacts are left untouched.
#
# SMOKE=1 env -> expanding-window / full-sample only, drop iter_2, for a fast
# end-to-end wiring check.  RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))  # hybrid_feature_selection (pca_gs dep)
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/u6_pool_prep.R")                        # PD5/PD7 — canonical U6 pool prep

SMOKE      <- nzchar(Sys.getenv("SMOKE"))
ONLY_ITER2 <- nzchar(Sys.getenv("ONLY_ITER2"))           # run ONLY the slow iter_2 S2 cells
SKIP_ITER2 <- (nzchar(Sys.getenv("SKIP_ITER2")) || SMOKE) && !ONLY_ITER2  # iter_2 (per-T FS) is slow
# WS4 (P8) forecast-retaining mode: restrict the sweep to the residual region
# (equity S2, post_2008, h in {3,12}, configs {C-OLS,C-COMB,C-CT}; all iters) and
# RETAIN the pointwise walk-forward forecasts so overlap-robust SEs + FDR/SPA can
# be computed downstream. Identical cell construction to the full run -> the
# scalar re-run must reproduce the stored cells exactly (Stage-B G1 guard).
RESIDUAL_FC    <- nzchar(Sys.getenv("RESIDUAL_FC"))
RESID_HORIZONS <- c(3L, 12L); RESID_CFGS <- c("C-OLS", "C-COMB", "C-CT")
AS_OF        <- as.Date("2026-04-30")
W_Z          <- 60L
ROLL_LENGTHS <- c(12L, 24L, 36L, 48L, 60L, 84L)   # U1 (2026-06-06): front-dense, capped at 7yr
HORIZONS     <- c(1L, 3L, 12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
MIN_OOS_ROWS <- 24L
T_START_S2_CUTOFF <- as.Date("1990-01-31")
PCA_GS_DEFAULTS <- list(n_groups = 10L, min_var_explained = 0.85,
                        min_distance = 0.5, min_t_value = 1.645,
                        hac_lag = 0L, verbose = FALSE)
SURPRISE_COLS   <- c("RGDP_SURPRISE","CPI_SURPRISE","UNRATE_SURPRISE","INDPRO_SURPRISE")
EXCL_FROM_POOL  <- c("SP500TR","SHILLER_PRICE")
EQUITY_SP_COL   <- "SP500TR"

# ---- 1. panel (SP500TR is already a corpus series) -------------------------
cat("[1] panel ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |>
  pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR",
                        u6_extra_factor_ids(meta)))   # PD-D: + 18 SPF M0
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward")

# ---- 1.5 surprises (for Ch3) ----------------------------------------------
surp <- build_surprise_panel(AS_OF, meta = meta)
panel_levels <- panel_levels |> left_join(surp, by = "reference_date") |>
  mutate(CPI_SURPRISE    = ifelse(!is.na(CPI_SURPRISE) & abs(CPI_SURPRISE) > 50, NA_real_, CPI_SURPRISE),
         INDPRO_SURPRISE = ifelse(!is.na(INDPRO_SURPRISE) & abs(INDPRO_SURPRISE) > 15, NA_real_, INDPRO_SURPRISE))

# ---- 2. corrected equity target (month-end total return) -------------------
targets_eq <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h = h, sp_col = EQUITY_SP_COL)),
  as.character(HORIZONS))

# ---- 3. pools (TERM/DEFAULT needed by WG+CRR) ------------------------------
cat("[3] pools ...\n")
panel_levels <- bond_theory_pool_ff1989(panel_levels)$panel    # adds TERM_SPREAD/DEFAULT_SPREAD
s1_eq     <- equity_theory_pool_wg_crr(panel_levels)
s1_eq_aug <- equity_theory_pool_wg_crr_aug(panel_levels)
# U6 (PD-D/PD7): derive the 13 dashboard constructs, then build TWO purged pools
# per chapter — full 114-var (shrinkage) + deep cohort (OLS/PCA-GS). SP500TR and
# SHILLER_PRICE are purged from the predictor pool (index must not predict itself).
dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
purge <- function(s2) { s2$ratio_features <- setdiff(s2$ratio_features, EXCL_FROM_POOL)
                        s2$level_features <- setdiff(s2$level_features, EXCL_FROM_POOL); s2 }
mk_full <- function(extra) purge(build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features,
             additional_candidates = extra))
mk_deep <- function(extra) purge(build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features,
             additional_candidates = extra))
cols2 <- function(s2) c(s2$ratio_features, s2$level_features)
s2f_ch2 <- mk_full(character(0));  s2f_ch3 <- mk_full(SURPRISE_COLS)
s2d_ch2 <- mk_deep(character(0));  s2d_ch3 <- mk_deep(SURPRISE_COLS)
S2_POOL_FULL <- list(ch2 = cols2(s2f_ch2), ch3 = cols2(s2f_ch3))           # shrinkage
S2_POOL_DEEP <- list(ch2 = live_cols(panel_levels, cols2(s2d_ch2)),        # OLS/PCA-GS, live
                     ch3 = live_cols(panel_levels, cols2(s2d_ch3)))
cat(sprintf("    S2 full ch2=%d ch3=%d | deep(live) ch2=%d ch3=%d (SP500TR/SHILLER_PRICE purged; dropped %s)\n",
            length(S2_POOL_FULL$ch2), length(S2_POOL_FULL$ch3),
            length(S2_POOL_DEEP$ch2), length(S2_POOL_DEEP$ch3),
            paste(setdiff(cols2(s2d_ch3), S2_POOL_DEEP$ch3), collapse=", ")))

# ---- 4. transforms ---------------------------------------------------------
all_ratio <- unique(c(s1_eq$ratio_features, s1_eq_aug$ratio_features,
                      s2f_ch2$ratio_features, s2f_ch3$ratio_features))
all_level <- unique(c(s1_eq$level_features, s1_eq_aug$level_features,
                      s2f_ch2$level_features, s2f_ch3$level_features))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
  apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)

# ---- 5. helpers ------------------------------------------------------------
first_non_na <- function(x, d){ i <- which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for <- function(cols) max(as.Date(sapply(cols, function(c)
  as.character(first_non_na(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval_from <- function(ts) seq(ts, by = "120 months", length.out = 2L)[2L]
sample_bounds <- function(s, ts){ d <- t_eval_from(ts); switch(s,
  full=list(t_eval=d,t_end=AS_OF), pre_2008=list(t_eval=d,t_end=SPLIT_2008-1L),
  post_2008=list(t_eval=max(d,SPLIT_2008),t_end=AS_OF),
  post_2016=list(t_eval=max(d,SPLIT_2016),t_end=AS_OF)) }

# ---- 6. FS pre-cache per chapter -------------------------------------------
cat("[6] FS pre-cache ...\n")
# PCA-GS on the live deep cohort; no silent swallow (no-silent-fails).
pca_gs_call <- function(pool, slice, y) do.call(pca_gs_select,
  c(list(panel=slice, target_y=y, predictor_cols=pool), PCA_GS_DEFAULTS))
fs_cache <- list()
for (ch in c("ch2","ch3")) {
  pool <- S2_POOL_DEEP[[ch]]
  sl_full <- panel_xfm |> filter(reference_date>=as.Date("1990-01-31"), reference_date<=AS_OF)
  mfull <- panel_xfm$reference_date>=as.Date("1990-01-31") & panel_xfm$reference_date<=AS_OF
  fs_cache[[paste(ch,"iter_1_full_sample","0",sep="|")]] <-
    pca_gs_call(pool, sl_full, targets_eq[["1"]][mfull])$selected_features
  sl16 <- panel_xfm |> filter(reference_date>=SPLIT_2016, reference_date<=AS_OF)
  m16 <- panel_xfm$reference_date>=SPLIT_2016 & panel_xfm$reference_date<=AS_OF
  fs_cache[[paste(ch,"iter_3_equal_window_2016","0",sep="|")]] <-
    pca_gs_call(pool, sl16, targets_eq[["1"]][m16])$selected_features
  for (h in HORIZONS)
    fs_cache[[paste(ch,"iter_4_multi_horizon",as.character(h),sep="|")]] <-
      pca_gs_call(pool, sl_full, targets_eq[[as.character(h)]][mfull])$selected_features
  cat(sprintf("    %s iter_1=%d iter_3=%d iter_4[1,3,12]=[%d,%d,%d]\n", ch,
    length(fs_cache[[paste(ch,"iter_1_full_sample","0",sep="|")]]),
    length(fs_cache[[paste(ch,"iter_3_equal_window_2016","0",sep="|")]]),
    length(fs_cache[[paste(ch,"iter_4_multi_horizon","1",sep="|")]]),
    length(fs_cache[[paste(ch,"iter_4_multi_horizon","3",sep="|")]]),
    length(fs_cache[[paste(ch,"iter_4_multi_horizon","12",sep="|")]])))
}
features_for_s2 <- function(ch, iter, h){
  if (iter=="iter_shrink_full")    return(S2_POOL_FULL[[ch]])   # full pool (shrinkage)
  if (iter=="iter_2_walk_forward") return(S2_POOL_DEEP[[ch]])   # deep cohort (OLS walk-fwd FS)
  fs_cache[[paste(ch, iter, if(iter=="iter_4_multi_horizon") as.character(h) else "0", sep="|")]]
}

# ---- 7. cell runner --------------------------------------------------------
run_cell <- function(y_h, panel_cell, h, wk, wl, ts, te, tend, cfg, iter, feats,
                     scenario, chapter, predset_label, enet_pool=character(0),
                     shrink_pool=character(0)) {
  out <- tryCatch(walk_forward_fs(y=y_h, panel=panel_cell, horizon=h, window_kind=wk,
    window_length=if(wk=="rolling") wl else NULL, t_start=ts, t_eval_start=te, t_end=tend,
    model_config=cfg, iteration_kind=iter, features=feats, target="equity",
    pca_gs_args=PCA_GS_DEFAULTS, enet_pool=enet_pool, shrink_pool=shrink_pool),
    error=function(e) list(forecasts=tibble()))
  fl <- apply_feasibility_floor(out$forecasts, h, cfg)  # U2/U3 on realised n_train
  wf <- fl$wf
  base <- tibble(chapter=chapter, scenario=scenario, target="equity", maturity=NA_character_,
    predset_label=predset_label, model_config=cfg, iteration=iter, horizon=h,
    window_kind=wk, window_length=ifelse(wk=="rolling", as.integer(wl), NA_integer_), n_oos=nrow(wf),
    n_train_min=if (nrow(wf)) min(wf$n_train) else NA_integer_,
    n_train_med=if (nrow(wf)) as.integer(round(median(wf$n_train))) else NA_integer_)
  if (nrow(wf) < MIN_OOS_ROWS)
    return(bind_cols(base, tibble(sample=NA_character_, r2_oos=NA_real_, cw=NA_real_,
      sharpe=NA_real_, accept=NA,
      note=if (fl$n_pre >= MIN_OOS_ROWS) "below_neff_floor" else "too_few_oos_rows")))
  r2 <- r2_oos(wf$y_hat,wf$y_realized,wf$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(wf$y_hat,wf$y_realized,wf$y_bench), L=floor(1.5*h))$stat
  sh <- signal_sharpe(wf$y_hat,wf$y_realized,h=h)
  out_row <- bind_cols(base, tibble(sample=NA_character_, r2_oos=r2, cw=cw, sharpe=sh,
    accept=accept_decision(r2,cw,sh), note=""))
  if (RESIDUAL_FC || scenario == "S2_equity")       # retain pointwise forecasts (WS4 + PD4)
    attr(out_row, "fc") <- tibble(chapter=chapter, scenario=scenario, model_config=cfg,
      iteration=iter, horizon=h, window_kind=wk,
      window_length=ifelse(wk=="rolling", as.integer(wl), NA_integer_),
      as_of_date=wf$as_of_date, y_hat=wf$y_hat, y_realized=wf$y_realized, y_bench=wf$y_bench)
  out_row
}

# ---- 8. specs (equity only; ch2 + ch3) -------------------------------------
mk_specs <- function() {
  o <- list(); iters <- c("iter_1_full_sample","iter_2_walk_forward",
                          "iter_3_equal_window_2016","iter_4_multi_horizon")
  for (ch in c("ch2","ch3")) {
    s1f <- if (ch=="ch2") c(s1_eq$ratio_features, s1_eq$level_features)
           else           c(s1_eq_aug$ratio_features, s1_eq_aug$level_features)
    s1_scn <- if (ch=="ch2") "S1_equity_WGCRR" else "S1_equity_WGCRR_aug"
    if (!ONLY_ITER2) for (cfg in c("C-OLS","C-COMB","C-ENET","C-CT"))
      o[[length(o)+1L]] <- list(chapter=ch, scenario=s1_scn, cfg=cfg, iter=NA_character_,
        feats_by_h=setNames(rep(list(s1f),length(HORIZONS)),as.character(HORIZONS)),
        pool_cols=s1f)
    use_iters <- if (ONLY_ITER2) "iter_2_walk_forward"
                 else if (SKIP_ITER2) setdiff(iters, "iter_2_walk_forward") else iters
    # PD2/PD7: shrinkage (ENET/RIDGE/PLS) runs once on the full pool under
    # iter_shrink_full; OLS-class crosses the PCA-GS iters on the deep cohort.
    if (!ONLY_ITER2) for (cfg in c("C-ENET","C-RIDGE","C-PLS"))
      o[[length(o)+1L]] <- list(chapter=ch, scenario="S2_equity", cfg=cfg, iter="iter_shrink_full",
        feats_by_h=setNames(lapply(HORIZONS, function(h) features_for_s2(ch,"iter_shrink_full",h)),as.character(HORIZONS)),
        pool_cols=S2_POOL_FULL[[ch]], t_start_cols=S2_POOL_DEEP[[ch]])
    for (cfg in c("C-OLS","C-COMB","C-CT")) for (iter in use_iters)
      o[[length(o)+1L]] <- list(chapter=ch, scenario="S2_equity", cfg=cfg, iter=iter,
        feats_by_h=setNames(lapply(HORIZONS, function(h) features_for_s2(ch,iter,h)),as.character(HORIZONS)),
        pool_cols=S2_POOL_DEEP[[ch]], t_start_cols=S2_POOL_DEEP[[ch]])
  }
  o
}
specs <- mk_specs()
if (RESIDUAL_FC) specs <- Filter(function(s) s$scenario=="S2_equity" && s$cfg %in% RESID_CFGS, specs)

SAMPLES <- if (SMOKE) "full" else if (RESIDUAL_FC) "post_2008" else c("full","pre_2008","post_2008","post_2016")
ALL_WINDOWS <- if (SMOKE) list(c("expanding",NA)) else
  c(list(c("expanding",NA)), lapply(ROLL_LENGTHS, function(w) c("rolling",w)))
pwr <- function(cols) panel_xfm |> select(reference_date, all_of(cols))
cat(sprintf("[8] %d specs x %d windows x %d samples%s\n", length(specs),
            length(ALL_WINDOWS), length(SAMPLES), if(SMOKE) " [SMOKE]" else ""))

# ---- 9. sweep --------------------------------------------------------------
results <- list(); fc_store <- list(); t0 <- Sys.time()
for (spec in specs) {
  if (SMOKE && isTRUE(spec$iter=="iter_2_walk_forward")) next
  for (h in HORIZONS) {
    if (RESIDUAL_FC && !(h %in% RESID_HORIZONS)) next
    feats <- spec$feats_by_h[[as.character(h)]]
    if (length(feats)==0L) next
    panel_cell <- pwr(spec$pool_cols)
    ts <- t_start_for(if (!is.null(spec$t_start_cols)) spec$t_start_cols else spec$pool_cols)
    is_shrink <- spec$cfg %in% c("C-ENET","C-RIDGE","C-PLS")
    plabel <- if (is_shrink && spec$scenario=="S2_equity") "s2_full_pool"
              else paste(feats, collapse="+")
    for (sl in SAMPLES) {
      sb <- sample_bounds(sl, ts); if (sb$t_eval >= sb$t_end) next
      for (win in ALL_WINDOWS) {
        row <- run_cell(targets_eq[[as.character(h)]],
          panel_cell, h, win[[1]], if(win[[1]]=="rolling") as.integer(win[[2]]) else NA_integer_,
          ts, sb$t_eval, sb$t_end, spec$cfg, spec$iter, feats, spec$scenario, spec$chapter,
          plabel, enet_pool=if(spec$cfg=="C-ENET") spec$pool_cols else character(0),
          shrink_pool=if(spec$cfg %in% c("C-RIDGE","C-PLS")) spec$pool_cols else character(0))
        fc <- attr(row, "fc")                          # PD4: retain all S2_equity forecasts
        if (!is.null(fc)) { fc$sample <- sl; fc_store[[length(fc_store)+1L]] <- fc }
        row$sample <- sl; results[[length(results)+1L]] <- row
      }
    }
  }
  if (length(results) %% 400L == 0L)
    cat(sprintf("    ... %d cells, %.0f s\n", length(results), as.numeric(Sys.time()-t0,units="secs")))
}
res <- bind_rows(results)
cat(sprintf("[9] sweep done: %d cells, %.0f s\n", nrow(res), as.numeric(Sys.time()-t0,units="secs")))

# ---- 10. persist + summary -------------------------------------------------
suffix <- if (SMOKE) "_SMOKE" else if (RESIDUAL_FC) "_residual" else if (ONLY_ITER2) "_iter2only" else if (SKIP_ITER2) "_noiter2" else ""
if (SKIP_ITER2 && !SMOKE) cat("[note] iter_2 (per-T walk-forward FS) EXCLUDED this run — slow; run separately.\n")
out_path <- file.path("data/audit", sprintf("rerun_equity_corrected_%s%s.parquet", Sys.Date(), suffix))
write_parquet(res, out_path)
cat(sprintf("[10] wrote %s (%d rows)\n", out_path, nrow(res)))
if (length(fc_store)) {        # WS4 residual region (RESIDUAL_FC) OR PD4 equity-S2 retention
  fc_all  <- bind_rows(fc_store)
  fc_name <- if (RESIDUAL_FC) "residual_forecasts" else "equity_s2_forecasts"
  fc_path <- file.path("data/audit", sprintf("%s_%s.parquet", fc_name, Sys.Date()))
  write_parquet(fc_all, fc_path)
  cat(sprintf("[10b] wrote %s (%d forecast rows across %d cells)\n", fc_path, nrow(fc_all),
      nrow(distinct(fc_all, chapter, scenario, model_config, iteration, horizon,
                    window_kind, window_length, sample))))
}
cat("\n=== acceptance by chapter x scenario (corrected SP500TR target) ===\n")
print(res |> group_by(chapter, scenario) |>
  summarise(n_cells=n(), n_accept=sum(accept,na.rm=TRUE),
            pass_pct=round(100*sum(accept,na.rm=TRUE)/n(),1),
            mean_r2=round(mean(r2_oos,na.rm=TRUE),4),
            mean_cw=round(mean(cw,na.rm=TRUE),2), .groups="drop"), n=Inf)
cat(sprintf("\nTotal equity accept: %d / %d\n", sum(res$accept,na.rm=TRUE), nrow(res)))

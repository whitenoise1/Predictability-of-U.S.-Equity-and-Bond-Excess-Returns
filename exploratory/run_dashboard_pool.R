# working_paper/run_dashboard_pool.R  — M1 practitioner-dashboard sensitivity (PD3).
#
# REWRITTEN for P-D (2026-06-08). The standalone S3a/S3b/S3c dashboard scenarios
# are RETIRED: under U6 the 13 derived practitioner constructs fold into the S2
# pool (tagged dashboard_subset). This driver now produces the M1 with/without
# sensitivity COUNTERFACTUAL — the WITHOUT-dashboard arm — to be diffed against the
# WITH-dashboard cells already produced by the main runs (p1_chapter_2.R bond +
# rerun_equity_corrected.R equity).
#
# Key fact (PD3): dashboard constructs live ONLY in the full shrinkage pool
# (s2_full); the deep cohort excludes them (1990 cutoff vs ~2006 proxy start), so
# the OLS-class S2 cells are identical with/without dashboard. M1 is therefore a
# SHRINKAGE-ONLY question. This run = shrinkage class (C-ENET/C-RIDGE/C-PLS) on the
# 101-var NO-dashboard pool (= 114 - 13), ch2-style (no surprises), equity + 5 bond
# maturities, U1 grid + U2/U3 floors. The acceptance delta vs the WITH cells is the
# M1 result (computed downstream in the rescore/M1 report).
# RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/u6_pool_prep.R")

AS_OF        <- as.Date("2026-04-30"); W_Z <- 60L
ROLL_LENGTHS <- c(12L, 24L, 36L, 48L, 60L, 84L)   # U1
HORIZONS     <- c(1L, 3L, 12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
MIN_OOS_ROWS <- 24L
T_START_S2_CUTOFF <- as.Date("1990-01-31")
SHRINK_CONFIGS <- c("C-ENET", "C-RIDGE", "C-PLS")
EXCL_FROM_POOL <- c("SP500TR", "SHILLER_PRICE")
MATS <- list(bond_1y  = list(y_col="DGS1",  duration=0.97), bond_2y = list(y_col="DGS2", duration=1.93),
             bond_3y  = list(y_col="DGS3",  duration=2.85), bond_5y = list(y_col="DGS5", duration=4.50),
             bond_10y = list(y_col="DGS10", duration=8.00))

# ---- 1. panel (macro + SP500TR + 18 SPF M0; NO dashboard constructs) --------
cat("[1] panel ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
panel_levels <- bond_theory_pool_ff1989(panel_levels)$panel

# ---- 2. targets (equity = SP500TR month-end TR; 5 bond maturities) ----------
cat("[2] targets ...\n")
targets_eq <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h = h, sp_col = "SP500TR")), as.character(HORIZONS))
targets_bd <- bond_target_grid(panel_levels, HORIZONS, MATS)
TARGETS <- c(list(equity = targets_eq), targets_bd)

# ---- 3. NO-dashboard U6 pool (shrinkage) + deep-cohort t_start anchor --------
# Canonical U6 call WITHOUT dashboard_features: 101 = 114 - 13. ch2-style (no
# surprises). The deep cohort (for the t_start anchor) is dashboard-free already.
cat("[3] no-dashboard pool ...\n")
cols2 <- function(s2) c(s2$ratio_features, s2$level_features)
s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
                                 forecast_token = U6_FORECAST_TOKEN)       # unpurged
s2_deep <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
                                 forecast_token = U6_FORECAST_TOKEN)
# Per-target purge to MATCH each WITH-arm convention so the M1 delta isolates the
# 13 dashboard cols and nothing else: ch2 BOND pool is UNPURGED (114 → WITHOUT 101);
# rerun EQUITY pool purges SP500TR/SHILLER_PRICE (112 → WITHOUT 99).
SHRINK_POOL_BD <- cols2(s2_full)                             # bond: unpurged (101)
SHRINK_POOL_EQ <- setdiff(SHRINK_POOL_BD, EXCL_FROM_POOL)    # equity: purged (99)
DEEP_COLS      <- live_cols(panel_levels, cols2(s2_deep))    # t_start anchor (~1990)
cat(sprintf("    no-dashboard pool: bond=%d (unpurged) equity=%d (purged); deep anchor=%d live\n",
            length(SHRINK_POOL_BD), length(SHRINK_POOL_EQ), length(DEEP_COLS)))

# ---- 4. transforms ----------------------------------------------------------
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = s2_full$ratio_features) |>
  apply_transform(rolling_zscore, factor_ids = s2_full$ratio_features, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = s2_full$level_features, window = W_Z, min_obs = W_Z)

# ---- 5. helpers (identical conventions to ch2) ------------------------------
first_non_na <- function(x, d){ i <- which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for  <- function(cols) max(as.Date(sapply(cols, function(c)
  as.character(first_non_na(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval_from  <- function(ts) seq(ts, by = "120 months", length.out = 2L)[2L]
sample_bounds <- function(s, ts){ d <- t_eval_from(ts); switch(s,
  full=list(t_eval=d,t_end=AS_OF), pre_2008=list(t_eval=d,t_end=SPLIT_2008-1L),
  post_2008=list(t_eval=max(d,SPLIT_2008),t_end=AS_OF),
  post_2016=list(t_eval=max(d,SPLIT_2016),t_end=AS_OF)) }
TS <- t_start_for(DEEP_COLS)                          # shrinkage t_start anchor (~1990)
panel_cell <- panel_xfm |> select(reference_date, all_of(SHRINK_POOL_BD))   # superset (101)

# ---- 6. cell runner (shrinkage only) ----------------------------------------
run_cell <- function(y_h, h, wk, wl, te, tend, cfg, target, maturity) {
  pool <- if (target=="equity") SHRINK_POOL_EQ else SHRINK_POOL_BD          # per-target purge
  out <- tryCatch(walk_forward_fs(y = y_h, panel = panel_cell, horizon = h,
      window_kind = wk, window_length = if (wk=="rolling") wl else NULL,
      t_start = TS, t_eval_start = te, t_end = tend, model_config = cfg,
      iteration_kind = "iter_shrink_full",
      target = if (target=="equity") "equity" else "bond",
      enet_pool   = if (cfg=="C-ENET") pool else character(0),
      shrink_pool = if (cfg %in% c("C-RIDGE","C-PLS")) pool else character(0)),
    error = function(e) list(forecasts = tibble()))
  fl <- apply_feasibility_floor(out$forecasts, h, cfg); wf <- fl$wf
  base <- tibble(scenario = "M1_nodash", target = target, maturity = maturity,
    predset_kind = "s2_shrink_nodash", predset_label = "s2_full_pool_nodash",
    model_config = cfg, iteration = "iter_shrink_full", horizon = h, window_kind = wk,
    window_length = ifelse(wk=="rolling", as.integer(wl), NA_integer_), n_oos = nrow(wf),
    n_train_min = if (nrow(wf)) min(wf$n_train) else NA_integer_,
    n_train_med = if (nrow(wf)) as.integer(round(median(wf$n_train))) else NA_integer_)
  if (nrow(wf) < MIN_OOS_ROWS)
    return(bind_cols(base, tibble(sample=NA_character_, r2_oos=NA_real_, cw=NA_real_,
      sharpe=NA_real_, accept=NA,
      note=if (fl$n_pre >= MIN_OOS_ROWS) "below_neff_floor" else "too_few_oos_rows")))
  r2 <- r2_oos(wf$y_hat, wf$y_realized, wf$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench), L=floor(1.5*h))$stat
  sh <- signal_sharpe(wf$y_hat, wf$y_realized, h=h)
  bind_cols(base, tibble(sample=NA_character_, r2_oos=r2, cw=cw, sharpe=sh,
                         accept=accept_decision(r2,cw,sh), note=""))
}

# ---- 7. sweep ---------------------------------------------------------------
SAMPLES     <- c("full","pre_2008","post_2008","post_2016")
ALL_WINDOWS <- c(list(c("expanding",NA)), lapply(ROLL_LENGTHS, function(w) c("rolling",w)))
cat(sprintf("[7] sweep: %d targets x %d configs x %d h x %d win x %d samp (shrinkage only)\n",
            length(TARGETS), length(SHRINK_CONFIGS), length(HORIZONS), length(ALL_WINDOWS), length(SAMPLES)))
results <- list(); t0 <- Sys.time()
for (tname in names(TARGETS)) {
  maturity <- if (tname=="equity") NA_character_ else tname
  for (cfg in SHRINK_CONFIGS) for (h in HORIZONS) {
    y_h <- TARGETS[[tname]][[as.character(h)]]
    for (sl in SAMPLES) {
      sb <- sample_bounds(sl, TS); if (sb$t_eval >= sb$t_end) next
      for (win in ALL_WINDOWS) {
        row <- run_cell(y_h, h, win[[1]], if(win[[1]]=="rolling") as.integer(win[[2]]) else NA_integer_,
                        sb$t_eval, sb$t_end, cfg, tname, maturity)
        row$sample <- sl; results[[length(results)+1L]] <- row
      }
    }
  }
  cat(sprintf("    %-9s done (%d cells, %.0fs)\n", tname, length(results),
              as.numeric(Sys.time()-t0, units="secs")))
}
res <- bind_rows(results)

# ---- 8. persist + summary ---------------------------------------------------
out_path <- file.path("data/audit", sprintf("dashboard_m1_nodash_%s.parquet", Sys.Date()))
write_parquet(res, out_path)
cat(sprintf("[8] wrote %s (%d cells) — diff vs WITH-dashboard cells (ch2/rerun) for M1\n",
            out_path, nrow(res)))
evl <- res |> filter(note == "" | is.na(note))
cat(sprintf("\n=== M1 no-dashboard arm (absolute-0.5 accept; P10 rescore downstream) ===\n"))
print(res |> mutate(grp = ifelse(target=="equity","equity","bond")) |>
  group_by(grp, model_config) |>
  summarise(evaluable = sum(note=="" | is.na(note)),
            below_floor = sum(note=="below_neff_floor", na.rm=TRUE),
            too_few     = sum(note=="too_few_oos_rows", na.rm=TRUE),
            n_accept = sum(accept, na.rm=TRUE),
            mean_r2  = round(mean(r2_oos, na.rm=TRUE), 4), .groups="drop"), n=Inf)
cat(sprintf("\nTotal evaluable %d / %d; below_neff_floor %d; too_few %d; accept %d\n",
            nrow(evl), nrow(res), sum(res$note=="below_neff_floor", na.rm=TRUE),
            sum(res$note=="too_few_oos_rows", na.rm=TRUE), sum(res$accept, na.rm=TRUE)))

# run_attribution_all.R — uniform variable-importance attribution across ALL
# procedures (consideration #3). Generalizes run_m3_attribution.R (which decomposed
# only the bond shrinkage class) so every config is attributed the same comparable
# way — no procedure black-boxed. Snapshot |coef| at the final FEASIBLE as-of, per
# evaluable p10 cell, for the in-scope configs × pools × targets.
#
# Stage S-A: snapshot coefficients only (the raw |coef| layer). GS selection
# frequency = S-B; combine + figure = S-C. Decisions (2026-06-15): OLS-class S2 is
# attributed on its COMMITTED iters (iter_1/3/4 fixed PCA-GS selection, NOT a fresh
# iter_2 walk-forward — supersedes the earlier "full walk-forward frequency" lock);
# C-COMB importance = |univariate slope|; measure-first STOP before the full sweep.
#
# Scope = chapter "ch2" only (level pools; avoids the ch3 surprise-panel rebuild,
# per run_m3). Authorities reproduced verbatim: bond = scripts/p1_chapter_2.R;
# equity = working_paper/rerun_equity_corrected.R (purges SP500TR/SHILLER_PRICE,
# retains S2_equity forecasts). Fidelity gates (full pass): (a) equity C-RIDGE
# snapshot vs retained equity_s2_forecasts <1e-8; (b) one bond OLS-class cell
# snapshot == walk_forward_fs last forecast; (b') one equity OLS-class cell snapshot
# vs retained forecast <1e-8. No re-implementation: reuses R/feature_selection.R +
# R/u6_pool_prep.R + R/walk_forward.R (feasibility_threshold).
#
# Invocation: measure-first (default) prints the cell matrix + extrapolated cost and
# STOPS. RUN_FULL=1 runs the full snapshot sweep + gates + writes the artifact.

suppressMessages({library(arrow); library(dplyr); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/u6_pool_prep.R")

RUN_FULL <- nzchar(Sys.getenv("RUN_FULL"))
AS_OF        <- as.Date("2026-04-30")
W_Z          <- 60L
HORIZONS     <- c(1L, 3L, 12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
T_START_S2_CUTOFF <- as.Date("1990-01-31")
EXCL_FROM_POOL <- c("SP500TR", "SHILLER_PRICE")    # purged for EQUITY only
PCA_GS_DEFAULTS <- list(n_groups = 10L, min_var_explained = 0.85, min_distance = 0.5,
                        min_t_value = 1.645, hac_lag = 0L, verbose = FALSE)
MATURITY_GRID <- list(bond_1y=list(y_col="DGS1",duration=0.97),
                      bond_2y=list(y_col="DGS2",duration=1.93),
                      bond_3y=list(y_col="DGS3",duration=2.85),
                      bond_5y=list(y_col="DGS5",duration=4.50),
                      bond_10y=list(y_col="DGS10",duration=8.00))
SHRINK_CFGS <- c("C-ENET","C-RIDGE","C-PLS")
OUT <- sprintf("data/audit/attribution_coefs_%s.parquet", AS_OF)

# ---- 1. panel + targets + pools (faithful to p1_chapter_2.R / rerun_equity) -----
cat("[1] panel + pools ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |>
  pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()

targets_eq <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h = h, sp_col = "SP500TR")),
  as.character(HORIZONS))
targets_bd <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)

s1_bond_ff <- bond_theory_pool_ff1989(panel_levels); panel_levels <- s1_bond_ff$panel
s1_eq      <- equity_theory_pool_wg_crr(panel_levels)
dash       <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
s2_deep <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
cols2 <- function(s2) c(s2$ratio_features, s2$level_features)
purge <- function(v) setdiff(v, EXCL_FROM_POOL)

# Per-(scenario,target) feature universes. Bond unpurged; equity purged.
BOND_FULL <- cols2(s2_full)                                   # S2 shrinkage bond
EQ_FULL   <- purge(cols2(s2_full))                            # S2 shrinkage equity
BOND_DEEP <- live_cols(panel_levels, cols2(s2_deep))          # S2 OLS-class bond
EQ_DEEP   <- live_cols(panel_levels, purge(cols2(s2_deep)))   # S2 OLS-class equity
S1_FF     <- s1_bond_ff$level_features                        # S1 bond FF1989
S1_CP     <- CP_FACTOR_COLS                                   # S1 bond CP (5 fwd rates)
S1_EQ     <- c(s1_eq$ratio_features, s1_eq$level_features)    # S1 equity WG∪CRR (unpurged — matches rerun_equity; purge is S2-only there)

# ---- 2. transform (yoy-log-diff ratios + 60-mo z-score; union over all pools) ---
# z-score is per-column, so values match either authoring runner as long as each
# fitted column is z-scored with the same ratio/level classification it had there.
all_ratio <- unique(c(s1_eq$ratio_features, s2_full$ratio_features))
all_level <- unique(c(s1_bond_ff$level_features, CP_FACTOR_COLS, s1_eq$level_features,
                      s2_full$level_features, "SP500TR", "SHILLER_PRICE"))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
  apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)
ref <- panel_xfm$reference_date

first_non_na <- function(x,d){ i<-which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for  <- function(cols) max(as.Date(sapply(cols, function(c)
  as.character(first_non_na(panel_xfm[[c]], ref)))))
# t_start anchors per scenario (S2 anchors on the deep cohort, S1 on its own pool).
TS <- list(S2_bond = t_start_for(BOND_DEEP), S2_equity = t_start_for(EQ_DEEP),
           S1_bond_FF1989 = t_start_for(S1_FF), S1_bond_CP = t_start_for(S1_CP),
           S1_equity_WGCRR = t_start_for(S1_EQ))
t_end_for <- function(s) if (s == "pre_2008") SPLIT_2008 - 1L else AS_OF

# ---- 3. reconstruct the OLS-class iter_1/3/4 PCA-GS selections (fs_cache) --------
# Identical slices/targets to p1_chapter_2.R:179-209 + rerun_equity_corrected.R:125-144.
cat("[2] reconstructing fs_cache (OLS-class iter_1/3/4 PCA-GS) ...\n")
pca_gs_call <- function(pool, slice, y) do.call(pca_gs_select,
  c(list(panel = slice, target_y = y, predictor_cols = pool), PCA_GS_DEFAULTS))$selected_features
fs_cache <- list()
m_full <- ref >= T_START_S2_CUTOFF & ref <= AS_OF; sl_full <- panel_xfm[m_full, , drop=FALSE]
m_16   <- ref >= SPLIT_2016        & ref <= AS_OF; sl_16   <- panel_xfm[m_16,   , drop=FALSE]
for (tk in c("equity", names(MATURITY_GRID))) {
  pool <- if (tk == "equity") EQ_DEEP else BOND_DEEP
  yL   <- if (tk == "equity") targets_eq else targets_bd[[tk]]
  fs_cache[[paste(tk,"iter_1_full_sample","0",sep="|")]]      <- pca_gs_call(pool, sl_full, yL[["1"]][m_full])
  fs_cache[[paste(tk,"iter_3_equal_window_2016","0",sep="|")]] <- pca_gs_call(pool, sl_16, yL[["1"]][m_16])
  for (h in HORIZONS)
    fs_cache[[paste(tk,"iter_4_multi_horizon",as.character(h),sep="|")]] <- pca_gs_call(pool, sl_full, yL[[as.character(h)]][m_full])
}

# ---- 4. per-config coefficient extraction (the native per-variable attribution) -
coef_of <- function(cfg, res, feats) {
  if (cfg == "C-RIDGE")      { v <- abs(res$beta); names(v) <- feats; v }            # std coefs
  else if (cfg == "C-PLS")   { v <- sqrt(rowSums(res$diag$loadings^2)); names(v) <- feats; v }  # net loading
  else if (cfg == "C-ENET")  { b <- res$beta[-1L]                                    # strip intercept
                               stopifnot(length(b) == length(feats)); v <- abs(b); names(v) <- feats; v }
  else if (cfg == "C-COMB")  { pf <- res$diag$per_feature; setNames(abs(pf$slope), pf$feature) }  # |univariate slope|
  else { b <- res$beta; b <- b[setdiff(names(b), "(Intercept)")]; abs(b) }           # C-OLS/C-CT/C-CP: |γ|
}

# ---- 5. snapshot fit at the FINAL feasible as-of (exact sweep path) --------------
# mode: "shrink" (PIT per-as-of col availability), "ols_fixed" (S1 / C-CP / OLS-class
# iter_1/3/4 fixed feats), "ols_iter2" (per-as-of PCA-GS, dormant unless p10 carries
# iter_2 cells). Lands on the last as-of that clears feasibility_threshold (the cell's
# last surviving forecast) — the as-of whose coefs the metrics are built on.
snapshot <- function(mode, x, ts, target, y_vec, cfg, h, wk, wl, sample) {
  tend <- t_end_for(sample)
  for (T_idx in rev(which(ref <= tend))) {
    if (is.na(y_vec[T_idx])) next
    fit_hi <- T_idx - h; if (fit_hi < 1L) next
    fit_lo <- if (wk == "rolling") max(fit_hi - wl + 1L, which(ref >= ts)[1L]) else which(ref >= ts)[1L]
    if (is.na(fit_lo) || fit_lo > fit_hi) next
    rows <- fit_lo:fit_hi
    if (mode == "shrink") {
      xT <- panel_xfm[T_idx, x, drop = FALSE]; win <- panel_xfm[rows, x, drop = FALSE]
      ok <- vapply(x, function(cc) !is.na(xT[[cc]]) && !anyNA(win[[cc]]), logical(1))
      feats <- x[ok]; if (length(feats) < 2L) next
    } else if (mode == "ols_iter2") {
      feats <- tryCatch(pca_gs_call(x, panel_xfm[rows, , drop=FALSE], y_vec[rows]), error = function(e) character(0))
      if (length(feats) < 2L) next
    } else feats <- x                                    # ols_fixed
    X_all <- as.matrix(panel_xfm[rows, feats, drop = FALSE]); y_t <- y_vec[rows]
    keep  <- stats::complete.cases(X_all) & !is.na(y_t)
    if (sum(keep) < feasibility_threshold(h, cfg, length(feats))) next
    xT_pred <- as.matrix(panel_xfm[T_idx, feats, drop = FALSE])
    if (any(is.na(xT_pred))) next
    res <- fit_and_predict(cfg, X_all[keep, , drop = FALSE], y_t[keep], xT_pred,
                           target = target, horizon = h)
    if (is.null(res)) next
    coef <- coef_of(cfg, res, feats)
    return(list(coef = coef, as_of = ref[T_idx], y_hat = res$y_hat, n = sum(keep), p = length(feats)))
  }
  NULL
}

# Resolve a p10 cell row -> snapshot arguments (mode/feature-universe/t_start/target/y).
resolve <- function(cl) {
  tg <- cl$target; scn <- cl$scenario; cfg <- cl$model_config; iter <- cl$iteration; h <- cl$horizon
  y_vec <- if (tg == "bond") targets_bd[[cl$maturity]][[as.character(h)]] else targets_eq[[as.character(h)]]
  ts <- TS[[scn]]
  # Shrinkage class (incl. the S1_equity C-ENET leg) runs the is_shrink PIT per-as-of
  # availability path on its pool, regardless of scenario.
  if (cfg %in% SHRINK_CFGS) {
    pool <- switch(scn, S2_bond = BOND_FULL, S2_equity = EQ_FULL, S1_equity_WGCRR = S1_EQ)
    return(list(mode="shrink", x=pool, ts=ts, target=tg, y=y_vec))
  }
  if (scn %in% c("S2_bond","S2_equity")) {                 # OLS-class S2 (deep-cohort PCA-GS)
    if (iter == "iter_2_walk_forward")
      return(list(mode="ols_iter2", x=if (tg=="bond") BOND_DEEP else EQ_DEEP, ts=ts, target=tg, y=y_vec))
    tk  <- if (tg == "bond") cl$maturity else "equity"
    key <- paste(tk, iter, if (iter == "iter_4_multi_horizon") as.character(h) else "0", sep="|")
    return(list(mode="ols_fixed", x=fs_cache[[key]], ts=ts, target=tg, y=y_vec))
  }
  feats <- switch(scn, S1_bond_FF1989 = S1_FF, S1_bond_CP = S1_CP,  # S1 OLS-class fixed pools
                  S1_equity_WGCRR = S1_EQ)                          # (C-OLS/C-COMB/C-CT; C-ENET took the shrink branch above)
  stopifnot(length(feats) > 0L)                                    # fail loud on an unresolved scenario (no silent empty fit)
  list(mode="ols_fixed", x=feats, ts=ts, target=tg, y=y_vec)
}

snap_cell <- function(cl) { a <- resolve(cl)
  snapshot(a$mode, a$x, a$ts, a$target, a$y, cl$model_config, cl$horizon,
           cl$window_kind, cl$window_length, cl$sample) }

# ---- 6. enumerate evaluable in-scope p10 cells (ch2) -----------------------------
p10 <- read_parquet("data/audit/p10_rescored_2026-06-12.parquet") |>
  filter(!is.na(accept_p10), chapter == "ch2")
INSCOPE <- p10 |> filter(
  (scenario == "S1_bond_FF1989"  & model_config %in% c("C-OLS","C-COMB")) |
  (scenario == "S1_bond_CP"      & model_config == "C-CP") |
  (scenario == "S2_bond"         & model_config %in% c("C-OLS","C-COMB", SHRINK_CFGS)) |
  (scenario == "S1_equity_WGCRR" & model_config %in% c("C-OLS","C-COMB","C-ENET","C-CT")) |
  (scenario == "S2_equity"       & model_config %in% c("C-OLS","C-COMB","C-CT", SHRINK_CFGS)))
cat(sprintf("[3] in-scope evaluable ch2 cells: %d / %d evaluable\n", nrow(INSCOPE), nrow(p10)))

# ---- 7. measure-first STOP (default) ---------------------------------------------
if (!RUN_FULL) {
  cat("\n=== in-scope cell matrix (scenario × config × target) ===\n")
  INSCOPE |> count(scenario, model_config, target, name = "cells") |>
    arrange(scenario, target, model_config) |> as.data.frame() |> print(row.names = FALSE)
  cat("\n=== by iteration ===\n")
  INSCOPE |> count(scenario, model_config, iteration, name = "cells") |>
    arrange(scenario, model_config, iteration) |> as.data.frame() |> print(row.names = FALSE)

  cat("\n[timing] sampling up to 3 cells per (config-class) ...\n")
  classes <- list(
    shrink   = INSCOPE |> filter(model_config %in% SHRINK_CFGS),
    ols_s2   = INSCOPE |> filter(grepl("^S2", scenario), model_config %in% c("C-OLS","C-COMB","C-CT")),
    ols_s1   = INSCOPE |> filter(scenario %in% c("S1_bond_FF1989","S1_equity_WGCRR")),
    cp       = INSCOPE |> filter(scenario == "S1_bond_CP"))
  per_class <- sapply(names(classes), function(nm) {
    d <- classes[[nm]]; if (!nrow(d)) return(NA_real_)
    idx <- unique(round(seq(1, nrow(d), length.out = min(3L, nrow(d)))))
    t0 <- Sys.time(); for (i in idx) invisible(snap_cell(d[i, ]))
    as.numeric(Sys.time() - t0, units = "secs") / length(idx)
  })
  cat("\n=== mean wall per cell by class (s) ===\n")
  for (nm in names(per_class)) cat(sprintf("    %-8s %s  (n=%d)\n", nm,
    if (is.na(per_class[nm])) "  —  " else sprintf("%6.3f", per_class[nm]), nrow(classes[[nm]])))
  est <- sum(mapply(function(nm) { r <- per_class[nm]; if (is.na(r)) 0 else r * nrow(classes[[nm]]) },
                    names(classes)))
  cat(sprintf("\n>>> estimated full-sweep wall: ~%.0f s (~%.1f min) over %d cells.\n",
              est, est/60, nrow(INSCOPE)))
  cat(">>> STOP (measure-first). Re-run with RUN_FULL=1 to execute the full sweep + gates + write.\n")
  quit(save = "no", status = 0)
}

# ---- 8. full snapshot sweep ------------------------------------------------------
cat("[4] full snapshot sweep ...\n")
t0 <- Sys.time(); acc <- list(); n_null <- 0L
for (i in seq_len(nrow(INSCOPE))) {
  cl <- INSCOPE[i, ]; s <- snap_cell(cl)
  if (is.null(s)) { n_null <- n_null + 1L; next }
  acc[[length(acc)+1L]] <- tibble(
    target = cl$target, scenario = cl$scenario, model_config = cl$model_config,
    maturity = cl$maturity, iteration = cl$iteration, horizon = cl$horizon,
    window_kind = cl$window_kind, window_length = cl$window_length, sample = cl$sample,
    variable = names(s$coef), coef = as.numeric(s$coef), as_of = s$as_of, n = s$n, p = s$p)
  if (i %% 200L == 0L) cat(sprintf("    ... %d / %d cells, %.0f s\n", i, nrow(INSCOPE),
                                   as.numeric(Sys.time() - t0, units = "secs")))
}
attr_long <- bind_rows(acc)
cat(sprintf("[4] %d cells snapshotted, %d returned no feasible as-of (%.0f s)\n",
            nrow(distinct(attr_long, scenario, model_config, maturity, iteration, horizon,
                          window_kind, window_length, sample)), n_null,
            as.numeric(Sys.time() - t0, units = "secs")))
stopifnot(nrow(attr_long) > 0)
fam <- meta |> select(variable = factor_id, family) |> distinct()
attribution <- attr_long |> left_join(fam, by = "variable")

# ---- 9. fidelity gates -----------------------------------------------------------
fc <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
gate_against_fc <- function(cfg) {                       # (a) shrink, (b') equity OLS-class
  g <- fc |> filter(model_config == cfg, window_kind == "rolling") |>
    count(iteration, window_length, sample, horizon, sort = TRUE) |> slice(1)
  st <- fc |> filter(model_config == cfg, iteration == g$iteration[1], horizon == g$horizon[1],
                     window_kind == "rolling", window_length == g$window_length[1],
                     sample == g$sample[1]) |> filter(as_of_date == max(as_of_date))
  gs <- snap_cell(tibble(target="equity", scenario="S2_equity", model_config=cfg,
    maturity=NA_character_, iteration=g$iteration[1], horizon=g$horizon[1],
    window_kind="rolling", window_length=g$window_length[1], sample=g$sample[1]))
  d <- if (!is.null(gs) && nrow(st) == 1L) abs(gs$y_hat - st$y_hat) else NA_real_
  cat(sprintf("    GATE %-7s equity %s iter=%s rolling-%d %s @ %s : |Δŷ| = %.2e\n",
              cfg, cfg, g$iteration[1], g$window_length[1], g$sample[1],
              as.character(st$as_of_date), d))
  if (is.na(d) || d >= 1e-8) stop(sprintf("attribution fidelity gate FAILED for %s (<1e-8).", cfg))
}
cat("[5] fidelity gates ...\n")
gate_against_fc("C-RIDGE")                               # (a) proven shrink gate
gate_against_fc("C-OLS")                                 # (b') equity OLS-class vs retained

# (b) bond OLS-class self-consistency: snapshot == walk_forward_fs last forecast.
bcell <- INSCOPE |> filter(scenario == "S2_bond", model_config == "C-OLS",
                           iteration != "iter_2_walk_forward") |>
  arrange(window_kind, desc(window_length)) |> slice(1)
if (nrow(bcell) == 1L) {
  tk  <- bcell$maturity; key <- paste(tk, bcell$iteration,
    if (bcell$iteration == "iter_4_multi_horizon") as.character(bcell$horizon) else "0", sep="|")
  feats <- fs_cache[[key]]
  wf <- walk_forward_fs(y = targets_bd[[tk]][[as.character(bcell$horizon)]],
    panel = panel_xfm |> select(reference_date, all_of(feats)), horizon = bcell$horizon,
    window_kind = bcell$window_kind,
    window_length = if (bcell$window_kind == "rolling") bcell$window_length else NULL,
    t_start = TS$S2_bond, t_eval_start = TS$S2_bond, t_end = t_end_for(bcell$sample),
    model_config = "C-OLS", iteration_kind = bcell$iteration, features = feats, target = "bond")
  fl  <- apply_feasibility_floor(wf$forecasts, bcell$horizon, "C-OLS")$wf
  last_wf <- fl |> filter(as_of_date == max(as_of_date))
  gs <- snap_cell(bcell)
  d  <- if (!is.null(gs) && nrow(last_wf) == 1L) abs(gs$y_hat - last_wf$y_hat) else NA_real_
  cat(sprintf("    GATE C-OLS   bond %s %s iter=%s @ %s : |Δŷ| = %.2e\n", tk, bcell$window_kind,
              bcell$iteration, as.character(last_wf$as_of_date), d))
  if (is.na(d) || d >= 1e-8) stop("attribution fidelity gate FAILED for bond C-OLS self-consistency (<1e-8).")
} else cat("    GATE C-OLS   bond: no eligible cell found — SKIPPED (flagged).\n")

# ---- 10. persist + summary -------------------------------------------------------
write_parquet(attribution, OUT)
cat(sprintf("[6] wrote %s (%d rows; %d cells × variables)\n", OUT, nrow(attribution),
            n_distinct(paste(attribution$scenario, attribution$model_config, attribution$maturity,
                             attribution$iteration, attribution$horizon, attribution$window_kind,
                             attribution$window_length, attribution$sample))))
cat("\n=== coverage: cells attributed per (scenario, config, target) ===\n")
attribution |> distinct(scenario, model_config, target, maturity, iteration, horizon,
                        window_kind, window_length, sample) |>
  count(scenario, model_config, target, name = "cells") |>
  arrange(scenario, target, model_config) |> as.data.frame() |> print(row.names = FALSE)

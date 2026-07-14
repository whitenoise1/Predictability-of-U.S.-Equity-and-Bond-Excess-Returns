# run_m3_attribution.R — M3 under-the-hood attribution (manuscript Pass 2).
#
# The shrinkage `beta` field (C-RIDGE standardized coefficients, C-PLS loadings) was
# diagnostic-only during the P-D sweep and never persisted. This bounded SNAPSHOT
# re-run (user decision 2026-06-13) re-fits C-RIDGE and C-PLS at the FINAL feasible
# as-of date for the live/accepted shrinkage cells, using the *exact* sweep fit
# (`fit_and_predict`) and the verbatim per-as-of column filter from `walk_forward_fs`,
# and dumps native per-variable attribution to data/audit/m3_attribution_<date>.parquet.
#
# Fidelity gate: for an equity C-RIDGE cell (which retains pointwise forecasts in
# equity_s2_forecasts), the snapshot y_hat must reproduce the stored forecast to <1e-8,
# proving the snapshot path == the sweep path (bond uses identical code, no retained
# forecasts). No re-implementation: reuses R/feature_selection.R + R/u6_pool_prep.R.

suppressMessages({library(arrow); library(dplyr); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/u6_pool_prep.R")

AS_OF        <- as.Date("2026-04-30")
W_Z          <- 60L
HORIZONS     <- c(1L, 3L, 12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
T_START_S2_CUTOFF <- as.Date("1990-01-31")
EXCL_FROM_POOL <- c("SP500TR", "SHILLER_PRICE")   # purged for EQUITY only (index can't predict itself)
MATURITY_GRID <- list(bond_1y=list(y_col="DGS1",duration=0.97),
                      bond_2y=list(y_col="DGS2",duration=1.93),
                      bond_3y=list(y_col="DGS3",duration=2.85),
                      bond_5y=list(y_col="DGS5",duration=4.50),
                      bond_10y=list(y_col="DGS10",duration=8.00))
OUT <- sprintf("data/audit/m3_attribution_%s.parquet", AS_OF)

# ---- 1. panel + targets + pools (faithful to p1_chapter_2.R / rerun_equity) -----
cat("[1] panel ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |>
  pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()

targets_eq <- setNames(lapply(HORIZONS, function(h)        # sp_col="SP500TR": month-end TR
  construct_equity_excess_return(panel_levels, h = h, sp_col = "SP500TR")),
  as.character(HORIZONS))
targets_bd <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)

panel_levels <- bond_theory_pool_ff1989(panel_levels)$panel
invisible(equity_theory_pool_wg_crr(panel_levels))
dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
s2_deep <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
cols2  <- function(s2) c(s2$ratio_features, s2$level_features)
purge  <- function(v) setdiff(v, EXCL_FROM_POOL)
POOL <- list(bond = cols2(s2_full), equity = purge(cols2(s2_full)))   # bond unpurged; equity purged
cat(sprintf("    s2_full=%d  bond pool=%d  equity pool=%d\n",
            length(cols2(s2_full)), length(POOL$bond), length(POOL$equity)))

# ---- 2. transform (yoy-log-diff ratios + 60-mo z-score; same as runner) ---------
all_ratio <- s2_full$ratio_features
all_level <- unique(c(s2_full$level_features, "SP500TR", "SHILLER_PRICE"))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
  apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)
ref <- panel_xfm$reference_date

first_non_na <- function(x,d){ i<-which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for  <- function(cols) max(as.Date(sapply(cols, function(c)
  as.character(first_non_na(panel_xfm[[c]], ref)))))
TS <- list(bond = t_start_for(cols2(s2_deep)),
           equity = t_start_for(purge(cols2(s2_deep))))
t_end_for <- function(s) if (s == "pre_2008") SPLIT_2008 - 1L else AS_OF

# ---- 3. snapshot fit at the FINAL feasible as-of (exact sweep path) --------------
# Replicates walk_forward_fs's per-as-of column filter (lines 640-652) verbatim, then
# calls fit_and_predict (the identical sweep estimator). Returns the native per-variable
# attribution: |std coef| for ridge, sqrt(rowSums(loadings^2)) for PLS.
snapshot <- function(target, y_vec, cfg, h, wk, wl, sample) {
  pool <- POOL[[target]]; ts <- TS[[target]]; tend <- t_end_for(sample)
  for (T_idx in rev(which(ref <= tend))) {
    if (is.na(y_vec[T_idx])) next
    fit_hi <- T_idx - h; if (fit_hi < 1L) next
    fit_lo <- if (wk == "rolling") max(fit_hi - wl + 1L, which(ref >= ts)[1L]) else which(ref >= ts)[1L]
    if (is.na(fit_lo) || fit_lo > fit_hi) next
    rows <- fit_lo:fit_hi
    xT <- panel_xfm[T_idx, pool, drop = FALSE]; win <- panel_xfm[rows, pool, drop = FALSE]
    ok <- vapply(pool, function(cc) !is.na(xT[[cc]]) && !anyNA(win[[cc]]), logical(1))
    feats <- pool[ok]; if (length(feats) < 2L) next
    X_all <- as.matrix(panel_xfm[rows, feats, drop = FALSE]); y_t <- y_vec[rows]
    keep  <- stats::complete.cases(X_all) & !is.na(y_t)
    if (sum(keep) < 8L * h + 2L) next
    res <- fit_and_predict(cfg, X_all[keep, , drop = FALSE], y_t[keep],
                           as.matrix(panel_xfm[T_idx, feats, drop = FALSE]),
                           target = target, horizon = h)
    if (is.null(res)) next
    coef <- if (cfg == "C-RIDGE") abs(res$beta) else sqrt(rowSums(res$diag$loadings^2))
    names(coef) <- feats
    return(list(coef = coef, as_of = ref[T_idx], y_hat = res$y_hat,
                n = sum(keep), p = length(feats)))
  }
  NULL
}

# ---- 4. cells to attribute: accepted bond shrinkage + feasible equity shrinkage --
# Restrict to ch2 (level pools): the natural base for "what the shrinkage class leans
# on", and it avoids the ch3 surprise-augmented pool (which would need the surprise
# panel rebuilt). equity_s2_forecasts carries one row per (cell, chapter), so the
# fidelity gate also fixes chapter == "ch2" for a unique match.
p10 <- read_parquet("data/audit/p10_rescored_2026-06-12.parquet") |> filter(!is.na(accept_p10))
sh  <- p10 |> filter(model_config %in% c("C-RIDGE","C-PLS"), grepl("^S2", scenario),
                     chapter == "ch2")
cells <- bind_rows(
  sh |> filter(target == "bond", accept_p10 == 1),                 # bond: accepted (live)
  sh |> filter(target == "equity"))                                # equity: all feasible (accepts none)
cat(sprintf("[4] attributing %d ch2 cells (bond accepted=%d, equity feasible=%d)\n",
            nrow(cells), sum(cells$target=="bond"), sum(cells$target=="equity")))

acc <- list()
for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  yv <- if (cl$target == "bond") targets_bd[[cl$maturity]][[as.character(cl$horizon)]]
        else targets_eq[[as.character(cl$horizon)]]
  s <- snapshot(cl$target, yv, cl$model_config, cl$horizon, cl$window_kind,
                cl$window_length, cl$sample)
  if (is.null(s)) next
  acc[[length(acc)+1L]] <- tibble(target = cl$target, model_config = cl$model_config,
    variable = names(s$coef), coef = as.numeric(s$coef))
}
attr_long <- bind_rows(acc)
stopifnot(nrow(attr_long) > 0)

fam <- meta |> select(variable = factor_id, family) |> distinct()
attribution <- attr_long |>
  group_by(target, model_config, variable) |>
  summarise(mean_abs_coef = mean(coef), n_cells = n(), .groups = "drop") |>
  left_join(fam, by = "variable") |>
  group_by(target, model_config) |>
  arrange(desc(mean_abs_coef), .by_group = TRUE) |>
  mutate(rank = row_number()) |> ungroup()

# ---- 5. fidelity gate: equity C-RIDGE snapshot vs retained forecast (<1e-8) ------
fc <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
gate <- fc |> filter(model_config == "C-RIDGE", window_kind == "rolling") |>
  count(window_length, sample, sort = TRUE) |> slice(1)
yv <- targets_eq[["1"]]
gs <- snapshot("equity", yv, "C-RIDGE", 1L, "rolling", gate$window_length[1], gate$sample[1])
stored <- fc |> filter(model_config == "C-RIDGE", horizon == 1L, window_kind == "rolling",
                       window_length == gate$window_length[1], sample == gate$sample[1]) |>
  filter(as_of_date == max(as_of_date))
delta <- if (!is.null(gs) && nrow(stored) == 1L) abs(gs$y_hat - stored$y_hat) else NA_real_
cat(sprintf("[5] FIDELITY equity C-RIDGE rolling-%d %s @ %s : |Δy_hat| = %.2e\n",
            gate$window_length[1], gate$sample[1], as.character(stored$as_of_date), delta))
if (is.na(delta) || delta >= 1e-8)
  stop("M3 fidelity gate FAILED: snapshot does not reproduce the retained forecast (<1e-8).")

write_parquet(attribution, OUT)
cat(sprintf("[6] wrote %s (%d rows)\n", OUT, nrow(attribution)))
cat("\n=== top 8 per (target, config) ===\n")
attribution |> filter(rank <= 8) |>
  select(target, model_config, rank, variable, family, mean_abs_coef, n_cells) |>
  as.data.frame() |> print(row.names = FALSE)

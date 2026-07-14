# run_family_signal_its_oos.R — the a-priori(ITS) -> PIT-OOS family-signal exhibit.
# Procedure-FREE, univariate/marginal, aggregated to economic family. Replaces the
# removed (it1-contaminated) Figure 3 attribution heatmap with a coherent map of
# where the raw signal APPEARS to lie (in-sample) and what SURVIVES out of sample.
#
# Per raw series x in the S2 full pool, for one target/horizon:
#   ITS leg  = full-sample univariate in-sample R^2  of  y_{t+h} ~ x_t   (a-priori; LOOK-AHEAD — descriptive only)
#   OOS leg  = univariate one-regressor walk-forward (C-OLS) OOS R^2 vs the expanding historical mean
#              (the paper's own benchmark; computed by r2_oos over feasibility-floored forecasts)
# Aggregated to family by MEDIAN; per-series long table retained for IQR whiskers.
#
# Setup (panel / pools / transform / t_start) reproduced from run_attribution_all.R
# VERBATIM so the numbers sit on the same panel as the rest of the paper. No
# re-implementation: reuses walk_forward_fs / apply_feasibility_floor / r2_oos.
#
# Invocation: measure-first (default) builds the panel, times a few univariate
# walk-forwards, extrapolates the full cost, and STOPS. RUN_FULL=1 runs the full
# sweep + the G1 self-consistency gate + writes the artifact.

suppressMessages({library(arrow); library(dplyr); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/u6_pool_prep.R")

RUN_FULL          <- nzchar(Sys.getenv("RUN_FULL"))
AS_OF             <- as.Date("2026-04-30")
W_Z               <- 60L
HORIZONS          <- c(1L, 3L, 12L)
H_PRIMARY         <- 1L                       # decision: h=1 primary (where the apparent signal concentrates)
BOND_MATURITY     <- "bond_1y"                # decision: short-end maturity for the bond panel (the live region)
T_START_S2_CUTOFF <- as.Date("1990-01-31")
EVAL_WAIT_MONTHS  <- 120L                     # walk-forward waits 120 mo after t_start before first scored forecast
EXCL_FROM_POOL    <- c("SP500TR", "SHILLER_PRICE")    # purged for EQUITY only
MATURITY_GRID <- list(bond_1y=list(y_col="DGS1",duration=0.97),
                      bond_2y=list(y_col="DGS2",duration=1.93),
                      bond_3y=list(y_col="DGS3",duration=2.85),
                      bond_5y=list(y_col="DGS5",duration=4.50),
                      bond_10y=list(y_col="DGS10",duration=8.00))
OUT <- sprintf("data/audit/family_signal_its_oos_%s.parquet", AS_OF)
add_months <- function(d, n) seq(as.Date(d), by = paste(n, "months"), length.out = 2L)[2L]

# Family resolver + display labels — reused VERBATIM from scratch_fig3_s2_gs.R so the
# economic families reconcile exactly with the (old) Figure 3 grammar. Derived features
# (spreads, DASH_*, X_* cross-asset, forward rates) that have no raw meta family are
# routed here rather than dropped (no silent fail).
fam_of <- function(variable, family) dplyr::case_when(
  !is.na(family) ~ family,
  grepl("^FWD_", variable) | variable %in% c("TERM_SPREAD","DASH_CURVE") ~ "YIELD_CURVE",
  variable %in% c("DEFAULT_SPREAD","X_HYIG_MOM","DASH_DEFAULT") ~ "CREDIT_SPREADS",
  variable == "DASH_PMI_PROXY" ~ "GROWTH", variable == "DASH_REAL2Y" ~ "MONETARY_POLICY_RATES",
  grepl("^X_EURJPY", variable) ~ "FX_CURRENCY", grepl("^X_OILGOLD", variable) ~ "COMMODITIES",
  variable %in% c("X_STKBOND_VAL","X_STKBOND_MOM","X_EMDM_MOM","X_INTLUS_MOM") ~ "CROSS_ASSET",
  TRUE ~ "OTHER")
FAM_LABEL <- c(MONETARY_POLICY_RATES="Monetary policy / rates", YIELD_CURVE="Yield curve",
  CREDIT_SPREADS="Credit spreads", INFLATION="Inflation", GROWTH="Growth / activity",
  LABOR_MARKET="Labor market", FORECASTS="SPF forecasts", EQUITY_VOLATILITY="Equity volatility",
  EQUITY_SENTIMENT="Equity sentiment", EQUITY_VALUATION="Equity valuation",
  FINANCIAL_CONDITIONS_COMPOSITE="Financial conditions", POLICY_UNCERTAINTY_GEOPOLITICAL="Policy uncertainty",
  FX_CURRENCY="FX", COMMODITIES="Commodities", HOUSING="Housing",
  CONSUMPTION_WEALTH="Consumption / wealth", CROSS_ASSET="Cross-asset (pract.)")

# ---- 1. panel + targets + S2 full pool (faithful to run_attribution_all.R) -------
cat("[1] panel + pools ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |>
  pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", u6_extra_factor_ids(meta)))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()

targets_eq <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h = h, sp_col = "SP500TR")), as.character(HORIZONS))
targets_bd <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)

dash    <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
cols2  <- function(s2) c(s2$ratio_features, s2$level_features)
BOND_FULL <- cols2(s2_full)                       # S2 full pool (bond)
EQ_FULL   <- setdiff(cols2(s2_full), EXCL_FROM_POOL)   # S2 full pool (equity; purged)

# ---- 2. transform (yoy-log-diff ratios + 60-mo PIT z-score) ----------------------
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = s2_full$ratio_features) |>
  apply_transform(rolling_zscore, factor_ids = s2_full$ratio_features, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = s2_full$level_features, window = W_Z, min_obs = W_Z)
ref <- panel_xfm$reference_date

# Drop columns absent / all-NA in the transformed panel — REPORTED, never silently absorbed.
avail <- function(cols) {
  present <- intersect(cols, names(panel_xfm))
  nonempty <- present[vapply(present, function(c) any(!is.na(panel_xfm[[c]])), logical(1))]
  dropped <- setdiff(cols, nonempty)
  if (length(dropped)) cat(sprintf("    [avail] dropped %d/%d cols (absent or all-NA): %s\n",
                                    length(dropped), length(cols), paste(head(dropped, 8), collapse=", ")))
  nonempty
}
BOND_FULL <- avail(BOND_FULL); EQ_FULL <- avail(EQ_FULL)

T_START      <- T_START_S2_CUTOFF
T_EVAL_START <- add_months(T_START, EVAL_WAIT_MONTHS)
EVAL_MASK    <- ref >= T_EVAL_START & ref <= AS_OF   # (B) same-window ITS support
cat(sprintf("[1] pools: bond=%d series, equity=%d series | h=%d, bond=%s | eval %s..%s\n",
            length(BOND_FULL), length(EQ_FULL), H_PRIMARY, BOND_MATURITY,
            as.character(T_EVAL_START), as.character(AS_OF)))

# ---- 3. the two legs, per series -------------------------------------------------
its_r2 <- function(xcol, yvec, mask = NULL) {     # univariate in-sample R^2 (a-priori mirage)
  x <- panel_xfm[[xcol]]; ok <- is.finite(x) & is.finite(yvec)
  if (!is.null(mask)) ok <- ok & mask             # (B) restrict to the OOS eval window
  if (sum(ok) < 24L) return(c(r2 = NA_real_, n = sum(ok)))
  fit <- stats::lm(yvec[ok] ~ x[ok])
  c(r2 = summary(fit)$r.squared, n = sum(ok))
}
oos_r2_one <- function(xcol, yvec, target) {      # univariate C-OLS walk-forward OOS R^2 (honest)
  wf <- walk_forward_fs(
    y = yvec, panel = panel_xfm |> select(reference_date, all_of(xcol)), horizon = H_PRIMARY,
    window_kind = "expanding", t_start = T_START, t_eval_start = T_EVAL_START, t_end = AS_OF,
    model_config = "C-OLS", features = xcol, target = target)$forecasts
  if (!"y_hat" %in% names(wf) || nrow(wf) == 0L) return(c(r2 = NA_real_, n = 0))
  fl <- apply_feasibility_floor(wf, H_PRIMARY, "C-OLS")$wf
  if (!"y_hat" %in% names(fl)) return(c(r2 = NA_real_, n = 0))
  fl <- fl[is.finite(fl$y_hat), , drop = FALSE]
  if (nrow(fl) < 2L) return(c(r2 = NA_real_, n = nrow(fl)))
  c(r2 = r2_oos(fl$y_hat, fl$y_realized, fl$y_bench), n = nrow(fl))
}
y_for <- function(target) if (target == "bond") targets_bd[[BOND_MATURITY]][[as.character(H_PRIMARY)]] else targets_eq[[as.character(H_PRIMARY)]]

# ---- 4. measure-first STOP (default) ---------------------------------------------
if (!RUN_FULL) {
  cat("\n[timing] sampling 3 series per target for the OOS walk-forward ...\n")
  samp <- function(v) v[unique(round(seq(1, length(v), length.out = min(3L, length(v)))))]
  time_target <- function(target, pool) {
    yv <- y_for(target); s <- samp(pool)
    t0 <- Sys.time(); for (c in s) invisible(oos_r2_one(c, yv, target))
    per <- as.numeric(Sys.time() - t0, units = "secs") / length(s)
    t1 <- Sys.time(); for (c in s) invisible(its_r2(c, yv))
    cat(sprintf("    %-7s OOS ~%.3f s/series, ITS ~%.4f s/series (n=%d series)\n",
                target, per, as.numeric(Sys.time() - t1, units="secs")/length(s), length(pool)))
    per * length(pool)
  }
  est <- time_target("bond", BOND_FULL) + time_target("equity", EQ_FULL)
  cat(sprintf("\n>>> estimated full-sweep wall: ~%.0f s (~%.1f min) over %d series-targets.\n",
              est, est/60, length(BOND_FULL) + length(EQ_FULL)))
  cat(">>> STOP (measure-first). Re-run with RUN_FULL=1 to execute the full sweep + G1 gate + write.\n")
  quit(save = "no", status = 0)
}

# ---- 5. full sweep ---------------------------------------------------------------
cat("[2] full sweep (both legs, all series, both targets) ...\n")
fam_map <- meta |> select(variable = factor_id, family) |> distinct()
sweep_target <- function(target, pool) {
  yv <- y_for(target)
  rows <- lapply(pool, function(c) {
    it <- its_r2(c, yv); itw <- its_r2(c, yv, EVAL_MASK); oo <- oos_r2_one(c, yv, target)
    tibble(target = target, variable = c, its_r2 = it["r2"], its_n = it["n"],
           itswin_r2 = itw["r2"], itswin_n = itw["n"], oos_r2 = oo["r2"], oos_n = oo["n"])
  })
  bind_rows(rows)
}
t0 <- Sys.time()
per_series <- bind_rows(sweep_target("bond", BOND_FULL), sweep_target("equity", EQ_FULL)) |>
  left_join(fam_map, by = "variable") |>
  mutate(fam = fam_of(variable, family), family_label = recode(fam, !!!FAM_LABEL))
n_other <- sum(per_series$fam == "OTHER")
cat(sprintf("[2] %d series-rows in %.0f s; %d routed to OTHER\n", nrow(per_series),
            as.numeric(Sys.time() - t0, units = "secs"), n_other))
if (n_other > 0) cat(sprintf("    [OTHER] %s\n",
  paste(unique(per_series$variable[per_series$fam == "OTHER"]), collapse = ", ")))
stopifnot(nrow(per_series) > 0)

per_family <- per_series |> group_by(target, fam, family_label) |>
  summarise(n_series = sum(!is.na(its_r2) | !is.na(oos_r2)),
            its_med = median(its_r2, na.rm = TRUE), oos_med = median(oos_r2, na.rm = TRUE),
            itswin_med = median(itswin_r2, na.rm = TRUE),
            its_q25 = quantile(its_r2, .25, na.rm = TRUE), its_q75 = quantile(its_r2, .75, na.rm = TRUE),
            itswin_q25 = quantile(itswin_r2, .25, na.rm = TRUE), itswin_q75 = quantile(itswin_r2, .75, na.rm = TRUE),
            oos_q25 = quantile(oos_r2, .25, na.rm = TRUE), oos_q75 = quantile(oos_r2, .75, na.rm = TRUE),
            horizon = H_PRIMARY, bond_maturity = BOND_MATURITY, .groups = "drop")

# ---- 6. G1 gate — univariate C-OLS == single-feature C-COMB forecasts (<1e-8) -----
cat("[3] G1 gate (univariate C-OLS forecasts == single-feature C-COMB) ...\n")
g1col <- BOND_FULL[1]; yv <- y_for("bond")
wf_ols  <- walk_forward_fs(y = yv, panel = panel_xfm |> select(reference_date, all_of(g1col)),
  horizon = H_PRIMARY, window_kind = "expanding", t_start = T_START, t_eval_start = T_EVAL_START,
  t_end = AS_OF, model_config = "C-OLS", features = g1col, target = "bond")$forecasts
wf_comb <- walk_forward_fs(y = yv, panel = panel_xfm |> select(reference_date, all_of(g1col)),
  horizon = H_PRIMARY, window_kind = "expanding", t_start = T_START, t_eval_start = T_EVAL_START,
  t_end = AS_OF, model_config = "C-COMB", features = g1col, target = "bond")$forecasts
j <- inner_join(wf_ols |> select(as_of_date, yh_ols = y_hat),
                wf_comb |> select(as_of_date, yh_comb = y_hat), by = "as_of_date")
d <- max(abs(j$yh_ols - j$yh_comb), na.rm = TRUE)
cat(sprintf("    GATE univariate C-OLS vs C-COMB on %s: |Δŷ| = %.2e over %d as-ofs\n", g1col, d, nrow(j)))
if (!is.finite(d) || d >= 1e-8) stop("G1 univariate-equivalence gate FAILED (<1e-8).")

# ---- 7. persist + summary --------------------------------------------------------
write_parquet(per_series, sub("\\.parquet$", "_perseries.parquet", OUT))
write_parquet(per_family, OUT)
cat(sprintf("[4] wrote %s (%d family rows) + per-series table\n", OUT, nrow(per_family)))
cat("\n=== family medians (ITS -> OOS), by target ===\n")
per_family |> arrange(target, desc(its_med)) |>
  transmute(target, family_label, n_series, its_med = round(100*its_med,2), oos_med = round(100*oos_med,2)) |>
  as.data.frame() |> print(row.names = FALSE)

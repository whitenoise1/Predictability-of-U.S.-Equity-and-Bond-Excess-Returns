# working_paper/rescore_p10.R  — WS1.
#
# Re-score the existing Chapter-2 / Chapter-3 acceptance cells under the P10
# target-relative Sharpe rule (WORKING_PAPER_PLAN.md P10), with NO new
# forecasts: only the Sharpe THRESHOLD changes, so re-evaluating
#   accept = (R2>0.01) & (CW>1.645) & (sharpe > threshold_target)
# on the stored r2_oos / cw / sharpe suffices.
#
#   Equity cells -> corrected-target rerun (rerun_equity_corrected_..._full).
#   Bond  cells  -> off-by-h-corrected 2026-05-30 chapter artifacts (their
#                   stale averaged-target equity rows are discarded).
#
# The threshold per cell needs the cell's OOS window [t_eval, t_end]; t_end is
# fixed by `sample`, t_eval = t_start + 120m, and t_start depends only on the
# scenario's full pool (NOT the selected features) — so it is reconstructed
# exactly by rebuilding the pools the generating drivers used. Two guards:
#   (G1) with threshold = 0.5 the re-score must reproduce the stored accept
#        on every evaluable cell (validates the r2/cw/sharpe read + logic);
#   (G2) reconstructed t_eval must imply the stored n_oos for clean S1 cells
#        (validates the t_start reconstruction).
# RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/u6_pool_prep.R")                     # PD-D — U6 admissions + live_cols
source("R/target_relative_sharpe.R")          # WS0

AS_OF        <- as.Date("2026-04-30"); W_Z <- 60L
HORIZONS     <- c(1L, 3L, 12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
T_START_S2_CUTOFF <- as.Date("1990-01-31")
SURPRISE_COLS  <- c("RGDP_SURPRISE","CPI_SURPRISE","UNRATE_SURPRISE","INDPRO_SURPRISE")
EXCL_FROM_POOL <- c("SP500TR","SHILLER_PRICE")
MATURITY_GRID  <- list(
  bond_1y=list(y_col="DGS1",duration=0.97),  bond_2y=list(y_col="DGS2",duration=1.93),
  bond_3y=list(y_col="DGS3",duration=2.85),  bond_5y=list(y_col="DGS5",duration=4.50),
  bond_10y=list(y_col="DGS10",duration=8.00))

# ---- 1. panel (superset corpus + SP500TR + SHILLER_PRICE + surprises) ------
cat("[1] panel + surprises ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors  <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR", u6_extra_factor_ids(meta)))  # +SPF M0
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
surp <- build_surprise_panel(AS_OF, meta = meta)
panel_levels <- panel_levels |> left_join(surp, by = "reference_date") |>
  mutate(CPI_SURPRISE    = ifelse(!is.na(CPI_SURPRISE) & abs(CPI_SURPRISE) > 50, NA_real_, CPI_SURPRISE),
         INDPRO_SURPRISE = ifelse(!is.na(INDPRO_SURPRISE) & abs(INDPRO_SURPRISE) > 15, NA_real_, INDPRO_SURPRISE))
ff <- bond_theory_pool_ff1989(panel_levels); panel_levels <- ff$panel  # adds TERM/DEFAULT spreads
dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel  # PD-D dashboard constructs

# ---- 2. targets (equity = SP500TR month-end TR; bond per maturity) ---------
cat("[2] targets ...\n")
eq_target <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h = h, sp_col = "SP500TR")), as.character(HORIZONS))
bd_target <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)

# ---- 3. pools (equity: corrected/purged; bond: pre-SP500TR universe) -------
cat("[3] pools ...\n")
purge <- function(p){ p$ratio_features <- setdiff(p$ratio_features, EXCL_FROM_POOL)
                      p$level_features <- setdiff(p$level_features, EXCL_FROM_POOL); p }
cols <- function(p) c(p$ratio_features, p$level_features)
s1_eq     <- equity_theory_pool_wg_crr(panel_levels)
s1_eq_aug <- equity_theory_pool_wg_crr_aug(panel_levels)
# PD-D: under U6 the S2 cells anchor t_start on the canonical DEEP COHORT
# (first-non-NA <= 1990, built with the U6 call — SPF M0 admitted, dashboard
# excluded by the cutoff). Equity is purged (SP500TR/SHILLER_PRICE), bond is not —
# mirroring rerun_equity_corrected.R vs ch2/ch3. live_cols (WURGLER) is irrelevant
# to t_start (it drops an early-starting series, not the max-first-non-NA).
mk_deep <- function(extra) build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
             forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features,
             additional_candidates = extra)
s2_eq_ch2 <- purge(mk_deep(character(0)));   s2_eq_ch3 <- purge(mk_deep(SURPRISE_COLS))
s2_bd_ch2 <-       mk_deep(character(0));     s2_bd_ch3 <-       mk_deep(SURPRISE_COLS)
cat(sprintf("    deep cohorts (t_start anchor): s2_eq_ch2=%d s2_eq_ch3=%d s2_bd_ch2=%d s2_bd_ch3=%d\n",
            length(cols(s2_eq_ch2)), length(cols(s2_eq_ch3)),
            length(cols(s2_bd_ch2)), length(cols(s2_bd_ch3))))

# ---- 4. transforms (union of all pools) ------------------------------------
cat("[4] transforms ...\n")
all_ratio <- unique(c(s1_eq$ratio_features, s1_eq_aug$ratio_features,
                      s2_eq_ch2$ratio_features, s2_eq_ch3$ratio_features,
                      s2_bd_ch2$ratio_features, s2_bd_ch3$ratio_features))
all_level <- unique(c(s1_eq$level_features, s1_eq_aug$level_features,
                      ff$level_features, CP_FACTOR_COLS,
                      s2_eq_ch2$level_features, s2_eq_ch3$level_features,
                      s2_bd_ch2$level_features, s2_bd_ch3$level_features))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
  apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)

# ---- 5. window helpers (identical to the drivers) + per-scenario t_start ----
first_non_na <- function(x, d){ i <- which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for  <- function(cs) max(as.Date(sapply(cs, function(c)
  as.character(first_non_na(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval_from  <- function(ts) seq(ts, by = "120 months", length.out = 2L)[2L]
sample_bounds <- function(s, ts){ d <- t_eval_from(ts); switch(s,
  full=list(t_eval=d,t_end=AS_OF), pre_2008=list(t_eval=d,t_end=SPLIT_2008-1L),
  post_2008=list(t_eval=max(d,SPLIT_2008),t_end=AS_OF),
  post_2016=list(t_eval=max(d,SPLIT_2016),t_end=AS_OF)) }

TSTART <- list(
  S1_bond_FF1989      = t_start_for(ff$level_features),
  S1_bond_CP          = t_start_for(CP_FACTOR_COLS),
  S1_equity_WGCRR     = t_start_for(cols(s1_eq)),
  S1_equity_WGCRR_aug = t_start_for(cols(s1_eq_aug)),
  S2_bond_ch2         = t_start_for(cols(s2_bd_ch2)),
  S2_bond_ch3         = t_start_for(cols(s2_bd_ch3)),
  S2_equity_ch2       = t_start_for(cols(s2_eq_ch2)),
  S2_equity_ch3       = t_start_for(cols(s2_eq_ch3)))
cat("[5] t_start by scenario:\n"); for (k in names(TSTART)) cat(sprintf("    %-22s %s\n", k, TSTART[[k]]))

tstart_key <- function(scenario, chapter)
  if (scenario %in% c("S2_bond","S2_equity")) paste(scenario, chapter, sep = "_") else scenario

# ---- 6. PIT expanding median-Sharpe paths per (target/maturity, h) ----------
cat("[6] Sharpe paths ...\n")
ref <- panel_levels$reference_date
sr_paths <- list()
for (h in HORIZONS) sr_paths[[paste("equity", h, sep="|")]] <-
  expanding_sharpe_path(eq_target[[as.character(h)]], ref, h)
for (m in names(MATURITY_GRID)) for (h in HORIZONS)
  sr_paths[[paste(m, h, sep="|")]] <- expanding_sharpe_path(bd_target[[m]][[as.character(h)]], ref, h)

threshold_for <- function(scenario, chapter, target, maturity, h, sample) {
  ts <- TSTART[[ tstart_key(scenario, chapter) ]]
  if (is.null(ts)) return(NA_real_)
  sb <- sample_bounds(sample, ts)
  pk <- if (target == "equity") paste("equity", h, sep="|") else paste(maturity, h, sep="|")
  median_sr_threshold(sr_paths[[pk]], sb$t_eval, sb$t_end)
}

# ---- 7. load cells: equity (rerun, SP500TR) + bond (ch2/ch3) — P-D artifacts -
cat("[7] load cells ...\n")
keep <- c("chapter","scenario","target","maturity","predset_label","model_config",
          "iteration","horizon","window_kind","window_length","n_oos","sample",
          "r2_oos","cw","sharpe","accept","note")
eq <- read_parquet("data/audit/rerun_equity_corrected_2026-06-12_noiter2.parquet")
b2 <- read_parquet("data/audit/p1_chapter_2_evaluation_2026-06-11.parquet") |>
  filter(target == "bond") |> mutate(chapter = "ch2")
b3 <- read_parquet("data/audit/p1_chapter_3_evaluation_2026-06-11.parquet") |>
  filter(target == "bond") |> mutate(chapter = "ch3")
cells <- bind_rows(eq |> select(any_of(keep)), b2 |> select(any_of(keep)), b3 |> select(any_of(keep)))
cat(sprintf("    %d cells (equity %d + bond %d)\n", nrow(cells),
            sum(cells$target=="equity"), sum(cells$target=="bond")))

# ---- 8. compute thresholds per distinct window combo, then re-score ---------
cat("[8] re-score ...\n")
combos <- cells |> filter(!is.na(sample)) |>
  distinct(scenario, chapter, target, maturity, horizon, sample)
combos$threshold <- mapply(threshold_for, combos$scenario, combos$chapter, combos$target,
                           combos$maturity, combos$horizon, combos$sample)
res <- cells |> left_join(combos, by = c("scenario","chapter","target","maturity","horizon","sample")) |>
  mutate(accept_05  = mapply(accept_decision,     r2_oos, cw, sharpe),
         accept_p10 = mapply(accept_decision_rel, r2_oos, cw, sharpe, threshold))

# ---- 9. GUARD G1: threshold=0.5 reproduces the stored accept exactly --------
ev <- res |> filter(!is.na(accept))
mism <- sum(ev$accept_05 != ev$accept)
if (mism > 0L) stop(sprintf("[G1 FAIL] 0.5-rule re-score disagrees with stored accept on %d/%d evaluable cells",
                            mism, nrow(ev)))
cat(sprintf("[G1 OK] 0.5-rule reproduces stored accept on all %d evaluable cells\n", nrow(ev)))

# ---- 10. GUARD G2: reconstructed t_eval implies stored n_oos (clean S1) -----
cat("[10] G2 t_eval/n_oos sanity (C-OLS, expanding, full sample):\n")
g2 <- res |> filter(window_kind=="expanding", sample=="full", model_config=="C-OLS",
                    scenario %in% c("S1_bond_FF1989","S1_bond_CP","S1_equity_WGCRR","S1_equity_WGCRR_aug"),
                    !is.na(n_oos))
g2chk <- g2 |> distinct(scenario, chapter, target, maturity, horizon, n_oos) |>
  rowwise() |>
  mutate(t_eval = sample_bounds("full", TSTART[[tstart_key(scenario, chapter)]])$t_eval,
         pk     = if (target=="equity") paste("equity",horizon,sep="|") else paste(maturity,horizon,sep="|"),
         pred_n = { r <- if (target=="equity") eq_target[[as.character(horizon)]] else bd_target[[maturity]][[as.character(horizon)]]
                    sum(ref >= t_eval & ref <= AS_OF & !is.na(r)) },
         gap = pred_n - n_oos) |> ungroup()
print(g2chk |> select(scenario, chapter, target, maturity, horizon, t_eval, n_oos, pred_n, gap), n = Inf)
if (any(abs(g2chk$gap) > 12L))
  cat(sprintf("    [G2 WARN] %d S1 cells with |pred_n - n_oos| > 12 — inspect t_start reconstruction\n",
              sum(abs(g2chk$gap) > 12L)))

# ---- 11. persist + headline -------------------------------------------------
out <- file.path("data/audit", sprintf("p10_rescored_%s.parquet", Sys.Date()))
write_parquet(res, out); cat(sprintf("[11] wrote %s (%d rows)\n", out, nrow(res)))
write_parquet(combos, file.path("data/audit", sprintf("p10_thresholds_%s.parquet", Sys.Date())))

hl <- res |> filter(!is.na(accept_p10)) |>
  group_by(chapter, scenario, target) |>
  summarise(evaluable = n(), acc_05 = sum(accept_05), acc_p10 = sum(accept_p10),
            pct_05 = round(100*mean(accept_05),1), pct_p10 = round(100*mean(accept_p10),1),
            .groups = "drop") |> arrange(target, chapter, scenario)
cat("\n=== Acceptance: 0.5-floor vs P10 target-relative (evaluable basis) ===\n")
print(hl, n = Inf)
cat("\n=== Roll-up by target ===\n")
print(res |> filter(!is.na(accept_p10)) |> group_by(target) |>
      summarise(evaluable=n(), acc_05=sum(accept_05), acc_p10=sum(accept_p10),
                pct_05=round(100*mean(accept_05),1), pct_p10=round(100*mean(accept_p10),1), .groups="drop"))
cat("\n=== threshold ranges by target ===\n")
print(combos |> group_by(target) |>
      summarise(min=round(min(threshold,na.rm=TRUE),3), med=round(median(threshold,na.rm=TRUE),3),
                max=round(max(threshold,na.rm=TRUE),3), .groups="drop"))
cat("\n[rescore_p10 complete]\n")

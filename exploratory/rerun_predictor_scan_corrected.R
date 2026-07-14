# working_paper/rerun_predictor_scan_corrected.R  — WS3.
#
# Corrected-target rebuild of the cross-chapter predictor-value scan (synthesis
# inputs). Mirrors scripts/predictor_value_scan.R with three corrections:
#   (1) equity target = SP500TR month-end TR (not Shiller monthly-average);
#   (2) SP500TR + SHILLER_PRICE purged from the candidate universe (P3);
#   (3) accept scored under BOTH the 0.5 floor and the P10 target-relative rule.
# Part A (selection frequency) is read from the P10-re-scored cells
# (data/audit/p10_rescored_2026-06-01.parquet) — equity (corrected target) +
# bond in one source. RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/target_relative_sharpe.R")

AS_OF        <- as.Date("2026-04-30"); W_Z <- 60L
HORIZONS_UNI <- c(1L, 12L)
MIN_OOS_ROWS <- 24L
T_START_S2_CUTOFF <- as.Date("1990-01-31")
SURPRISE_COLS  <- c("RGDP_SURPRISE","CPI_SURPRISE","UNRATE_SURPRISE","INDPRO_SURPRISE")
EXCL_FROM_POOL <- c("SP500TR","SHILLER_PRICE")
MATURITY_GRID  <- list(
  bond_1y=list(y_col="DGS1",duration=0.97),  bond_2y=list(y_col="DGS2",duration=1.93),
  bond_3y=list(y_col="DGS3",duration=2.85),  bond_5y=list(y_col="DGS5",duration=4.50),
  bond_10y=list(y_col="DGS10",duration=8.00))

# ============================================================================
# PART A — selection frequency from the P10-re-scored cells.
# ============================================================================
cat("[A] selection frequency (from p10_rescored) ...\n")
rs <- read_parquet("data/audit/p10_rescored_2026-06-01.parquet")
freq_from <- function(d, accept_col) {
  gs <- d |> filter(.data[[accept_col]] %in% TRUE,
                    grepl("^S2", scenario),
                    iteration %in% c("iter_1_full_sample","iter_3_equal_window_2016","iter_4_multi_horizon"),
                    model_config != "C-ENET",
                    !is.na(predset_label), predset_label != "", predset_label != "s2_full_pool")
  if (nrow(gs) == 0L) return(tibble())
  gs |> mutate(feat = strsplit(predset_label, "[+]")) |> unnest(feat) |>
    mutate(target_grp = ifelse(target == "equity", "equity", "bond")) |>
    group_by(chapter, target_grp, feat) |>
    summarise(n_accept_cells = n(), mean_cw = mean(cw, na.rm=TRUE),
              mean_r2 = mean(r2_oos, na.rm=TRUE), mean_sh = mean(sharpe, na.rm=TRUE), .groups="drop")
}
freq_p10 <- freq_from(rs, "accept_p10") |> mutate(rule = "P10")
freq_05  <- freq_from(rs, "accept_05")  |> mutate(rule = "floor_0.5")
freq_all <- bind_rows(freq_p10, freq_05)
write_parquet(freq_all, file.path("data/audit", sprintf("predictor_selection_frequency_%s.parquet", Sys.Date())))
cat(sprintf("    wrote selection-frequency (%d rows; P10 equity feats=%d, bond feats=%d)\n",
            nrow(freq_all), sum(freq_p10$target_grp=="equity"), sum(freq_p10$target_grp=="bond")))

# ============================================================================
# PART B — univariate standalone scan on the CORRECTED equity target.
# ============================================================================
cat("[B] panel + transforms (mirror chapter 3, corrected target) ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors  <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR"))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
surp <- build_surprise_panel(AS_OF, meta = meta)
panel_levels <- panel_levels |> left_join(surp, by="reference_date") |>
  mutate(CPI_SURPRISE    = ifelse(!is.na(CPI_SURPRISE) & abs(CPI_SURPRISE) > 50, NA_real_, CPI_SURPRISE),
         INDPRO_SURPRISE = ifelse(!is.na(INDPRO_SURPRISE) & abs(INDPRO_SURPRISE) > 15, NA_real_, INDPRO_SURPRISE))
ff <- bond_theory_pool_ff1989(panel_levels); panel_levels <- ff$panel
s1_eq <- equity_theory_pool_wg_crr_aug(panel_levels)
s2    <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
                               additional_candidates = SURPRISE_COLS)
ratio_features <- unique(c(s1_eq$ratio_features, s2$ratio_features))
level_features <- unique(c(ff$level_features, CP_FACTOR_COLS, s1_eq$level_features, s2$level_features,
                           intersect(SURPRISE_COLS, names(panel_levels))))
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = ratio_features) |>
  apply_transform(rolling_zscore, factor_ids = ratio_features, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = level_features, window = W_Z, min_obs = W_Z)
candidates <- setdiff(unique(c(ratio_features, level_features)), EXCL_FROM_POOL)
candidates <- candidates[candidates %in% names(panel_xfm)]
cat(sprintf("    %d candidates (purged %s)\n", length(candidates), paste(EXCL_FROM_POOL, collapse="/")))

# targets (corrected equity) + per-target Sharpe paths
eq_target <- setNames(lapply(HORIZONS_UNI, function(h)
  construct_equity_excess_return(panel_levels, h=h, sp_col="SP500TR")), as.character(HORIZONS_UNI))
bd_target <- bond_target_grid(panel_levels, HORIZONS_UNI, MATURITY_GRID)
target_set <- c(list(equity=eq_target), bd_target)
ref <- panel_levels$reference_date
sr_paths <- list()
for (tk in names(target_set)) for (h in HORIZONS_UNI)
  sr_paths[[paste(tk,h,sep="|")]] <- expanding_sharpe_path(target_set[[tk]][[as.character(h)]], ref, h)

first_non_na <- function(x,d){ i<-which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for1 <- function(col) first_non_na(panel_xfm[[col]], panel_xfm$reference_date)
t_eval_from  <- function(ts) seq(ts, by="120 months", length.out=2L)[2L]

cat("[B] univariate walk-forward scan ...\n")
rows <- list(); t0 <- Sys.time()
for (fac in candidates) {
  ts <- t_start_for1(fac); if (is.na(ts)) next
  te <- t_eval_from(ts);   if (te >= AS_OF) next
  pc <- panel_xfm |> select(reference_date, all_of(fac))
  for (tk in names(target_set)) for (h in HORIZONS_UNI) {
    y_h <- target_set[[tk]][[as.character(h)]]
    out <- tryCatch(walk_forward_fs(y=y_h, panel=pc, horizon=h, window_kind="expanding",
      window_length=NULL, t_start=ts, t_eval_start=te, t_end=AS_OF, model_config="C-OLS",
      iteration_kind=NA_character_, features=fac, target=if(tk=="equity")"equity" else "bond",
      pca_gs_args=list(), enet_pool=character(0)), error=function(e) list(forecasts=tibble()))
    wf <- out$forecasts; if (nrow(wf) < MIN_OOS_ROWS) next
    r2 <- r2_oos(wf$y_hat,wf$y_realized,wf$y_bench)
    cw <- clark_west_stat(clark_west_pointwise(wf$y_hat,wf$y_realized,wf$y_bench), L=floor(1.5*h))$stat
    sh <- signal_sharpe(wf$y_hat,wf$y_realized,h=h)
    thr <- median_sr_threshold(sr_paths[[paste(tk,h,sep="|")]], te, AS_OF)
    rows[[length(rows)+1L]] <- tibble(factor=fac, target=tk, horizon=h, n_oos=nrow(wf),
      r2_oos=r2, cw=cw, sharpe=sh, threshold=thr,
      accept_05=accept_decision(r2,cw,sh), accept_p10=accept_decision_rel(r2,cw,sh,thr))
  }
}
uni <- bind_rows(rows)
write_parquet(uni, file.path("data/audit", sprintf("predictor_univariate_scan_%s.parquet", Sys.Date())))
cat(sprintf("    wrote univariate scan (%d rows, %.0f s)\n", nrow(uni), as.numeric(Sys.time()-t0,units="secs")))

# ---- GUARD: bond rows must reproduce the prior scan (target-invariant) ------
old <- tryCatch(read_parquet("data/audit/predictor_univariate_scan_2026-05-31.parquet"), error=function(e) NULL)
if (!is.null(old)) {
  cmp <- uni |> filter(target != "equity") |>
    inner_join(old |> select(factor,target,horizon, r2_old=r2_oos, cw_old=cw, sh_old=sharpe),
               by=c("factor","target","horizon")) |>
    mutate(dr2=abs(r2_oos-r2_old), dcw=abs(cw-cw_old), dsh=abs(sharpe-sh_old))
  cat(sprintf("[guard] bond reproduction: %d matched rows, max |Δr2|=%.4f |Δcw|=%.4f |Δsh|=%.4f\n",
              nrow(cmp), max(cmp$dr2,na.rm=TRUE), max(cmp$dcw,na.rm=TRUE), max(cmp$dsh,na.rm=TRUE)))
  eqcmp <- uni |> filter(target=="equity") |>
    inner_join(old |> filter(target=="equity") |> select(factor,horizon, r2_old=r2_oos),
               by=c("factor","horizon"))
  cat(sprintf("[guard] equity SHOULD differ: matched %d, mean |Δr2|=%.3f (expect large)\n",
              nrow(eqcmp), mean(abs(eqcmp$r2_oos-eqcmp$r2_old),na.rm=TRUE)))
}

# ---- report ----------------------------------------------------------------
cat("\n=== univariate standalone accepters, EQUITY h=1 (corrected target) ===\n")
print(uni |> filter(target=="equity", horizon==1) |> arrange(desc(cw)) |>
      mutate(across(c(r2_oos,cw,sharpe,threshold), \(x) round(x,3))) |> head(15), n=15)
cat("\n=== univariate acceptance counts per (target, horizon) ===\n")
print(uni |> group_by(target, horizon) |>
      summarise(n=n(), acc05=sum(accept_05,na.rm=TRUE), accP10=sum(accept_p10,na.rm=TRUE), .groups="drop") |>
      arrange(target, horizon), n=Inf)
cat("\n[corrected predictor scan complete]\n")

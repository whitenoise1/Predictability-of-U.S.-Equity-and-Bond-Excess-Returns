# working_paper/run_surprise_only.R  — WS2.
#
# Surprise-only scenario (P9): completes the level / surprise / level+surprise
# triad. Chapter 2 = LEVEL pool; Chapter 3 = LEVEL + SURPRISE (degraded — but
# confounded by the PCA+GS reshuffle when 3-4 surprises enter a ~48-candidate
# pool). This run isolates the surprises: a FIXED small pool of the SPF-based
# macro surprises, run S1-style (configs applied DIRECTLY, no PCA+GS feature
# selection, no FS iterations), per target. If surprises carry standalone
# multivariate value, it shows here without the reshuffle confound.
#
# Pool   : {CPI_SURPRISE, INDPRO_SURPRISE, UNRATE_SURPRISE} — the surprises that
#          clear the 1990 t_start filter. RGDP_SURPRISE excluded (ALFRED first
#          GDPC1 vintage 1991-12 → fails the cutoff; documented in p1_chapter_3.R).
# Target : equity (corrected SP500TR month-end TR) + 5 bond maturities.
# Configs: equity C-OLS/C-COMB/C-ENET/C-CT; bond C-OLS/C-COMB/C-ENET.
# Grid   : 3 horizons x 7 windows x 4 samples. Surprises get the level treatment
#          (60-month rolling z-score), as in Chapter 3.
# Scoring: stored accept = 0.5 floor (comparability) + accept_p10 = P10 rule.
# RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R"); source("R/surprise.R")
source("R/target_relative_sharpe.R")

AS_OF        <- as.Date("2026-04-30"); W_Z <- 60L
ROLL_LENGTHS <- c(12L,24L,36L,48L,60L,84L)   # U1 (2026-06-06): front-dense, capped at 7yr
HORIZONS     <- c(1L,3L,12L)
SPLIT_2008   <- as.Date("2008-01-01"); SPLIT_2016 <- as.Date("2016-01-01")
MIN_OOS_ROWS <- 24L
SURPRISE_POOL <- c("CPI_SURPRISE","INDPRO_SURPRISE","UNRATE_SURPRISE")
MATURITY_GRID <- list(
  bond_1y=list(y_col="DGS1",duration=0.97),  bond_2y=list(y_col="DGS2",duration=1.93),
  bond_3y=list(y_col="DGS3",duration=2.85),  bond_5y=list(y_col="DGS5",duration=4.50),
  bond_10y=list(y_col="DGS10",duration=8.00))

# ---- 1. panel + surprises (mirror Chapter 3) -------------------------------
cat("[1] panel + surprises ...\n")
meta <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC","MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors  <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR"))
panel_levels <- build_panel(AS_OF, all_factors, "M", "forward") |> cp_factor_columns()
surp <- build_surprise_panel(AS_OF, meta = meta)
panel_levels <- panel_levels |> left_join(surp, by = "reference_date") |>
  mutate(CPI_SURPRISE    = ifelse(!is.na(CPI_SURPRISE) & abs(CPI_SURPRISE) > 50, NA_real_, CPI_SURPRISE),
         INDPRO_SURPRISE = ifelse(!is.na(INDPRO_SURPRISE) & abs(INDPRO_SURPRISE) > 15, NA_real_, INDPRO_SURPRISE))
stopifnot(all(SURPRISE_POOL %in% names(panel_levels)))
fobs <- sapply(SURPRISE_POOL, function(c){ i<-which(!is.na(panel_levels[[c]])); as.character(panel_levels$reference_date[i[1]]) })
cat(sprintf("    surprise pool first-non-NA: %s\n", paste(SURPRISE_POOL, fobs, sep="=", collapse=" | ")))

# ---- 2. targets ------------------------------------------------------------
cat("[2] targets ...\n")
eq_target <- setNames(lapply(HORIZONS, function(h)
  construct_equity_excess_return(panel_levels, h=h, sp_col="SP500TR")), as.character(HORIZONS))
bd_target <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)

# ---- 3. transforms (surprises = level treatment: 60m z-score) --------------
cat("[3] transforms ...\n")
panel_xfm <- panel_levels |>
  apply_transform(rolling_zscore, factor_ids = SURPRISE_POOL, window = W_Z, min_obs = W_Z)

# ---- 4. window helpers -----------------------------------------------------
first_non_na <- function(x,d){ i<-which(!is.na(x)); if(!length(i)) as.Date(NA) else d[i[1]] }
t_start_for  <- function(cs) max(as.Date(sapply(cs, function(c)
  as.character(first_non_na(panel_xfm[[c]], panel_xfm$reference_date)))))
t_eval_from  <- function(ts) seq(ts, by="120 months", length.out=2L)[2L]
sample_bounds <- function(s, ts){ d<-t_eval_from(ts); switch(s,
  full=list(t_eval=d,t_end=AS_OF), pre_2008=list(t_eval=d,t_end=SPLIT_2008-1L),
  post_2008=list(t_eval=max(d,SPLIT_2008),t_end=AS_OF),
  post_2016=list(t_eval=max(d,SPLIT_2016),t_end=AS_OF)) }
T_START <- t_start_for(SURPRISE_POOL)
cat(sprintf("[4] surprise-pool t_start=%s  t_eval(full)=%s\n", T_START, t_eval_from(T_START)))

# ---- 5. per-cell runner (mirrors p1_chapter_2.R::run_cell) -----------------
run_cell <- function(y_h, panel_cell, h, wk, wl, ts, te, tend, cfg, target, maturity) {
  out <- tryCatch(walk_forward_fs(y=y_h, panel=panel_cell, horizon=h, window_kind=wk,
    window_length=if(wk=="rolling") wl else NULL, t_start=ts, t_eval_start=te, t_end=tend,
    model_config=cfg, iteration_kind=NA_character_, features=SURPRISE_POOL, target=target,
    pca_gs_args=list(), enet_pool=if(cfg=="C-ENET") SURPRISE_POOL else character(0)),
    error=function(e) list(forecasts=tibble()))
  fl <- apply_feasibility_floor(out$forecasts, h, cfg)  # U2/U3 on realised n_train
  wf <- fl$wf
  base <- tibble(scenario="S3_surprise", target=target, maturity=maturity,
    model_config=cfg, iteration=NA_character_, horizon=h, window_kind=wk,
    window_length=ifelse(wk=="rolling", as.integer(wl), NA_integer_), n_oos=nrow(wf),
    n_train_min=if (nrow(wf)) min(wf$n_train) else NA_integer_,
    n_train_med=if (nrow(wf)) as.integer(round(median(wf$n_train))) else NA_integer_)
  if (nrow(wf) < MIN_OOS_ROWS)
    return(bind_cols(base, tibble(sample=NA_character_, r2_oos=NA_real_, cw=NA_real_,
      sharpe=NA_real_, accept=NA,
      note=if (fl$n_pre >= MIN_OOS_ROWS) "below_neff_floor" else "too_few_oos_rows")))
  r2 <- r2_oos(wf$y_hat,wf$y_realized,wf$y_bench)
  cw <- clark_west_stat(clark_west_pointwise(wf$y_hat,wf$y_realized,wf$y_bench), L=floor(1.5*h))$stat
  sh <- signal_sharpe(wf$y_hat,wf$y_realized,h=h)
  bind_cols(base, tibble(sample=NA_character_, r2_oos=r2, cw=cw, sharpe=sh,
    accept=accept_decision(r2,cw,sh), note=""))
}

# ---- 6. spec sweep ---------------------------------------------------------
SAMPLES <- c("full","pre_2008","post_2008","post_2016")
ALL_WINDOWS <- c(list(c("expanding",NA)), lapply(ROLL_LENGTHS, function(w) c("rolling",w)))
panel_cell <- panel_xfm |> select(reference_date, all_of(SURPRISE_POOL))
targets_by <- c(list(equity=eq_target), bd_target)
cfgs_for  <- function(tg) if (tg=="equity") c("C-OLS","C-COMB","C-ENET","C-CT") else c("C-OLS","C-COMB","C-ENET")

cat("[5] sweep ...\n"); results <- list(); t0 <- Sys.time()
for (tk in names(targets_by)) {
  tg  <- if (tk=="equity") "equity" else "bond"
  mat <- if (tk=="equity") NA_character_ else tk
  yL  <- targets_by[[tk]]
  for (cfg in cfgs_for(tg)) for (h in HORIZONS) {
    for (sl in SAMPLES) {
      sb <- sample_bounds(sl, T_START); if (sb$t_eval >= sb$t_end) next
      for (win in ALL_WINDOWS) {
        row <- run_cell(yL[[as.character(h)]], panel_cell, h, win[[1]],
          if(win[[1]]=="rolling") as.integer(win[[2]]) else NA_integer_,
          T_START, sb$t_eval, sb$t_end, cfg, tg, mat)
        row$sample <- sl; results[[length(results)+1L]] <- row
      }
    }
  }
}
res <- bind_rows(results)
cat(sprintf("[6] sweep done: %d cells, %.0f s\n", nrow(res), as.numeric(Sys.time()-t0,units="secs")))

# ---- 7. P10 threshold + accept_p10 (PIT median Sharpe path) ----------------
cat("[7] P10 scoring ...\n")
ref <- panel_levels$reference_date
sr_paths <- list()
for (h in HORIZONS) sr_paths[[paste("equity",h,sep="|")]] <-
  expanding_sharpe_path(eq_target[[as.character(h)]], ref, h)
for (m in names(MATURITY_GRID)) for (h in HORIZONS)
  sr_paths[[paste(m,h,sep="|")]] <- expanding_sharpe_path(bd_target[[m]][[as.character(h)]], ref, h)
thr_for <- function(target, maturity, h, sample) {
  sb <- sample_bounds(sample, T_START)
  pk <- if (target=="equity") paste("equity",h,sep="|") else paste(maturity,h,sep="|")
  median_sr_threshold(sr_paths[[pk]], sb$t_eval, sb$t_end)
}
combos <- res |> filter(!is.na(sample)) |> distinct(target, maturity, horizon, sample)
combos$threshold <- mapply(thr_for, combos$target, combos$maturity, combos$horizon, combos$sample)
res <- res |> left_join(combos, by=c("target","maturity","horizon","sample")) |>
  mutate(accept_05 = accept, accept_p10 = mapply(accept_decision_rel, r2_oos, cw, sharpe, threshold))

# ---- 8. persist + report ---------------------------------------------------
out <- file.path("data/audit", sprintf("p1_surprise_only_%s.parquet", Sys.Date()))
write_parquet(res, out); cat(sprintf("[8] wrote %s (%d rows)\n", out, nrow(res)))

cat("\n=== S3 surprise-only acceptance by target x config (evaluable basis) ===\n")
print(res |> filter(!is.na(accept_05)) |> group_by(target, model_config) |>
  summarise(eval=n(), acc05=sum(accept_05), accP10=sum(accept_p10),
            pct05=round(100*mean(accept_05),1), pctP10=round(100*mean(accept_p10),1),
            mean_r2=round(mean(r2_oos),3), mean_cw=round(mean(cw),2), .groups="drop"), n=Inf)
cat("\n=== roll-up by target ===\n")
print(res |> filter(!is.na(accept_05)) |> group_by(target) |>
  summarise(eval=n(), acc05=sum(accept_05), accP10=sum(accept_p10),
            pct05=round(100*mean(accept_05),1), pctP10=round(100*mean(accept_p10),1),
            best_r2=round(max(r2_oos),3), best_cw=round(max(cw),2), best_sh=round(max(sharpe),2),
            .groups="drop"), n=Inf)
cat(sprintf("\n=== TRIAD (S2 data-driven baselines from p10_rescored, P10 rule) ===\n"))
rs <- tryCatch(read_parquet("data/audit/p10_rescored_2026-06-01.parquet"), error=function(e) NULL)
if (!is.null(rs)) {
  tri <- rs |> filter(scenario=="S2_equity"|scenario=="S2_bond", !is.na(accept_p10)) |>
    mutate(chap=ifelse(chapter=="ch2","LEVEL (ch2 S2)","LEVEL+SURPRISE (ch3 S2)")) |>
    group_by(target, chap) |> summarise(pctP10=round(100*mean(accept_p10),1), .groups="drop")
  print(tri, n=Inf)
}
cat("    SURPRISE-ONLY (this run, P10): see roll-up above.\n")
cat("\n[surprise-only complete]\n")

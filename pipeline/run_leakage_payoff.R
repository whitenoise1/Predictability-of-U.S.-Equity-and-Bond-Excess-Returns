# working_paper/pipeline/run_leakage_payoff.R
#
# P1 "leakage payoff" (vintage-only). See working_paper/notes/EXEC_GUIDE_leakage_payoff.md.
#
# Question: does substituting REVISED (final-vintage) macro data for point-in-time
# data — release timing held identical — inflate measured acceptance? (The vintage
# analog of the equity averaging artifact 52.6%->0.0%; the payoff of Annex A.2.)
#
# Minimal slice: shrinkage class {C-RIDGE,C-PLS,C-ENET} on the Scenario-2 full pool,
# equity (ch2, no surprises) + bond_1y (ch2 levels), h=1, EXPANDING, FULL sample.
# Two panels: PIT (build_panel) vs LEAKED (build_panel_leaked_vintage). The feature
# set, t_start, and window are held FIXED from the PIT panel; only the data vintage
# varies. Realized targets (SP500TR / DGS yields) are single-vintage market data —
# identical across panels — so the P10 threshold is identical and ONLY y_hat moves.
#
# Guards (no-silent-fails):
#   - leak-applied assertion: revised cols differ, non-revised cols byte-identical.
#   - G1: the PIT leg reproduces committed p10_rescored_drop-it1_2026-06-17.parquet
#     (ch2 shrinkage / h1 / expanding / full) to max|Δ|=0 -> proves only the panel
#     differs. Hard stop on mismatch.
# SMOKE=1 -> only C-RIDGE (fast wiring check). RUN FROM REPO ROOT.

suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})
source("R/pit_query.R"); source("R/panel.R"); source("R/panel_leaked.R")
source("R/transforms.R"); source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R"); source("R/feature_pools.R")
source("R/u6_pool_prep.R"); source("R/target_relative_sharpe.R")

SMOKE   <- nzchar(Sys.getenv("SMOKE"))
CFGS    <- if (SMOKE) "C-RIDGE" else c("C-RIDGE", "C-PLS", "C-ENET")
# Slice axes are env-overridable so the full/h1 run and the pre_2008/h{1,3}
# extension share one generator (each writes a distinct, self-tagged artifact).
SAMPLES  <- strsplit(Sys.getenv("SAMPLES",  "full,pre_2008"), ",")[[1]]
HORIZONS <- as.integer(strsplit(Sys.getenv("HORIZONS", "1,3"), ",")[[1]])
AS_OF      <- as.Date("2026-04-30")
SPLIT_2008 <- as.Date("2008-01-01")
W_Z     <- 60L
MIN_OOS_ROWS      <- 24L
T_START_S2_CUTOFF <- as.Date("1990-01-31")
EXCL_FROM_POOL    <- c("SP500TR", "SHILLER_PRICE")
MATURITY_GRID     <- list(bond_1y = list(y_col = "DGS1", duration = 0.97))
sample_bounds <- function(s, ts) {                 # mirrors the canonical drivers
  d <- seq(ts, by = "120 months", length.out = 2L)[2L]
  switch(s, full = list(t_eval = d, t_end = AS_OF),
            pre_2008 = list(t_eval = d, t_end = SPLIT_2008 - 1L),
            stop("unknown sample: ", s))
}

meta   <- read_parquet("data/metadata/factor_metadata.parquet")
non_fc <- meta |> filter(public_private %in% c("PUBLIC", "MIXED"),
                         !family %in% "FORECASTS", !grepl("^RTDSM_", factor_id)) |> pull(factor_id)
all_factors <- unique(c(non_fc, "SHILLER_PRICE", "SP500TR", u6_extra_factor_ids(meta)))

# ---- 1. panels: PIT control + vintage-only leaked treatment ----------------
cat("[1] building PIT + leaked panels ...\n")
pit_levels    <- build_panel(AS_OF, all_factors, "M", "forward")
leaked_levels <- build_panel_leaked_vintage(AS_OF, all_factors, meta = meta)
stopifnot(identical(dim(pit_levels), dim(leaked_levels)),
          identical(pit_levels$reference_date, leaked_levels$reference_date))

# leak-applied assertion (no-silent-fails)
shared      <- setdiff(intersect(names(pit_levels), names(leaked_levels)), "reference_date")
revised_ids <- meta$factor_id[meta$vintage_available & meta$family != "FORECASTS"]
col_differs <- function(c) !isTRUE(all.equal(pit_levels[[c]], leaked_levels[[c]]))
diff_cols   <- shared[vapply(shared, col_differs, logical(1))]
nonrev      <- setdiff(shared, revised_ids)
stopifnot(length(diff_cols) > 0L,                       # leak actually did something
          all(diff_cols %in% revised_ids),              # only revised cols moved
          !any(vapply(nonrev, col_differs, logical(1)))) # non-revised cols identical
cat(sprintf("    leak applied: %d/%d revised cols in panel differ; %d non-revised cols identical\n",
            length(diff_cols), length(intersect(revised_ids, shared)), length(nonrev)))

# ---- 2. pools + transforms (identical recipe per panel) --------------------
# Mirrors rerun_equity_corrected.R (equity, purge SP500TR+SHILLER_PRICE) and
# rerun_bond_forecasts.R (bond, no SP500TR present -> purge only SP500TR here,
# SHILLER_PRICE kept as in the canonical bond pool). t_start uses the DEEP cohort.
build_pools <- function(levels) {
  levels <- cp_factor_columns(levels)
  levels <- bond_theory_pool_ff1989(levels)$panel
  s1eq   <- equity_theory_pool_wg_crr(levels)
  dash   <- derive_u6_dashboard(levels); levels <- dash$panel
  mk <- function(cut) build_scenario_2_pool(levels, meta, t_start_cutoff = cut,
          forecast_token = U6_FORECAST_TOKEN, dashboard_features = dash$dashboard_features)
  s2f <- mk(NULL); s2d <- mk(T_START_S2_CUTOFF)
  s2f_cols <- c(s2f$ratio_features, s2f$level_features)
  s2d_cols <- c(s2d$ratio_features, s2d$level_features)
  all_ratio <- unique(c(s1eq$ratio_features, s2f$ratio_features))
  all_level <- unique(c(s1eq$level_features, s2f$level_features))
  xfm <- levels |>
    apply_transform(yoy_log_diff,   factor_ids = all_ratio) |>
    apply_transform(rolling_zscore, factor_ids = all_ratio, window = W_Z, min_obs = W_Z) |>
    apply_transform(rolling_zscore, factor_ids = all_level, window = W_Z, min_obs = W_Z)
  fnn    <- function(x, d) { i <- which(!is.na(x)); if (!length(i)) as.Date(NA) else d[i[1L]] }
  tstart <- function(cols) max(as.Date(sapply(cols,
              function(c) as.character(fnn(xfm[[c]], xfm$reference_date)))))
  # Equity purges the index's own level (can't predict itself); bond purges
  # nothing (SP500TR/SHILLER_PRICE are legitimate bond predictors) — exactly as
  # rerun_equity_corrected.R vs rerun_bond_forecasts.R.
  eq_full <- setdiff(s2f_cols, EXCL_FROM_POOL)
  bd_full <- s2f_cols
  eq_deep <- live_cols(levels, setdiff(s2d_cols, EXCL_FROM_POOL))
  bd_deep <- live_cols(levels, s2d_cols)
  list(levels = levels, xfm = xfm,
       eq_full = eq_full, bd_full = bd_full,
       t_start_eq = tstart(eq_deep), t_start_bd = tstart(bd_deep))
}
cat("[2] pools + transforms ...\n")
pit_p <- build_pools(pit_levels); lk_p <- build_pools(leaked_levels)
# Feature set / window must be IDENTICAL across panels (else the comparison is a
# confound, not a clean vintage swap). Hold the experiment fixed on the PIT defs.
stopifnot(identical(pit_p$eq_full, lk_p$eq_full), identical(pit_p$bd_full, lk_p$bd_full),
          identical(pit_p$t_start_eq, lk_p$t_start_eq),
          identical(pit_p$t_start_bd, lk_p$t_start_bd))
cat(sprintf("    pools fixed: eq_full=%d bd_full=%d | t_start eq=%s bd=%s\n",
            length(pit_p$eq_full), length(pit_p$bd_full),
            pit_p$t_start_eq, pit_p$t_start_bd))

# ---- 3. targets (per horizon; panel-invariant: single-vintage markets) -----
ref     <- pit_p$xfm$reference_date
bd_grid <- bond_target_grid(pit_p$levels, c(1L, 3L, 12L), MATURITY_GRID)[["bond_1y"]]
TARGETS <- list(
  equity  = list(kind = "equity", pool = pit_p$eq_full, t_start = pit_p$t_start_eq,
                 scenario = "S2_equity", maturity = NA_character_,
                 y = setNames(lapply(HORIZONS, function(h)
                       construct_equity_excess_return(pit_p$levels, h = h, sp_col = "SP500TR")),
                     as.character(HORIZONS))),
  bond_1y = list(kind = "bond", pool = pit_p$bd_full, t_start = pit_p$t_start_bd,
                 scenario = "S2_bond", maturity = "bond_1y",
                 y = setNames(lapply(HORIZONS, function(h) bd_grid[[as.character(h)]]),
                     as.character(HORIZONS))))

# ---- 4. cell runner (one shrinkage cell) -----------------------------------
run_one <- function(xfm, pool, y, t_start, t_eval, t_end, h, cfg, kind) {
  panel_cell <- xfm |> select(reference_date, all_of(pool))
  out <- tryCatch(walk_forward_fs(y = y, panel = panel_cell, horizon = h,
           window_kind = "expanding", window_length = NULL, t_start = t_start,
           t_eval_start = t_eval, t_end = t_end, model_config = cfg,
           iteration_kind = "iter_shrink_full", features = pool, target = kind,
           enet_pool   = if (cfg == "C-ENET") pool else character(0),
           shrink_pool = if (cfg %in% c("C-RIDGE", "C-PLS")) pool else character(0)),
         error = function(e) list(forecasts = tibble()))
  fl <- apply_feasibility_floor(out$forecasts, h, cfg); wf <- fl$wf
  if (nrow(wf) < MIN_OOS_ROWS)
    return(tibble(r2_oos = NA_real_, cw = NA_real_, sharpe = NA_real_,
                  n_oos = nrow(wf), n_train_med = NA_integer_))
  tibble(r2_oos = r2_oos(wf$y_hat, wf$y_realized, wf$y_bench),
         cw     = clark_west_stat(clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench),
                                  L = floor(1.5 * h))$stat,
         sharpe = signal_sharpe(wf$y_hat, wf$y_realized, h = h),
         n_oos  = nrow(wf), n_train_med = as.integer(round(median(wf$n_train))))
}

cat(sprintf("[4] sweep: %d cfg x %d targets x %d horizons x %d samples x 2 panels = %d cells ...\n",
            length(CFGS), length(TARGETS), length(HORIZONS), length(SAMPLES),
            length(CFGS) * length(TARGETS) * length(HORIZONS) * length(SAMPLES) * 2L))
rows <- list(); t0 <- Sys.time()
for (panel_lab in c("pit", "leaked")) {
  xfm <- if (panel_lab == "pit") pit_p$xfm else lk_p$xfm
  for (tnm in names(TARGETS)) { tg <- TARGETS[[tnm]]
    for (h in HORIZONS) { y_h <- tg$y[[as.character(h)]]
      for (smp in SAMPLES) {
        sb <- sample_bounds(smp, tg$t_start); if (sb$t_eval >= sb$t_end) next
        thr <- median_sr_threshold(expanding_sharpe_path(y_h, ref, h), sb$t_eval, sb$t_end)
        for (cfg in CFGS) {
          m <- run_one(xfm, tg$pool, y_h, tg$t_start, sb$t_eval, sb$t_end, h, cfg, tg$kind)
          rows[[length(rows) + 1L]] <- tibble(
            panel = panel_lab, target = tnm, scenario = tg$scenario, maturity = tg$maturity,
            model_config = cfg, horizon = h, window_kind = "expanding", sample = smp,
            r2_oos = m$r2_oos, cw = m$cw, sharpe = m$sharpe, threshold = thr,
            accept_p10 = accept_decision_rel(m$r2_oos, m$cw, m$sharpe, thr),
            n_oos = m$n_oos, n_train_med = m$n_train_med)
          cat(sprintf("    %-6s %-7s %-7s %-8s h%d r2=%+.4f cw=%+.2f sh=%+.3f acc=%s (%.0fs)\n",
                      panel_lab, tnm, cfg, smp, h, m$r2_oos, m$cw, m$sharpe,
                      accept_decision_rel(m$r2_oos, m$cw, m$sharpe, thr),
                      as.numeric(Sys.time() - t0, units = "secs")))
        }}}}
}
res <- bind_rows(rows)

# ---- 5. GUARD G1: PIT leg reproduces committed drop-it1 (ch2) --------------
cat("[5] G1: pin PIT leg to drop-it1 ...\n")
pit_res <- res |> filter(panel == "pit")
stopifnot(!any(is.na(pit_res$r2_oos)))     # every slice we ran must be evaluable
ref_cells <- read_parquet("data/audit/p10_rescored_drop-it1_2026-06-17.parquet") |>
  filter(chapter == "ch2", model_config %in% CFGS, window_kind == "expanding",
         sample %in% SAMPLES, horizon %in% HORIZONS,
         (scenario == "S2_equity" & is.na(maturity)) |
         (scenario == "S2_bond"   & maturity == "bond_1y")) |>
  transmute(target = if_else(scenario == "S2_equity", "equity", "bond_1y"),
            model_config, horizon, sample, ref_r2 = r2_oos, ref_cw = cw, ref_sh = sharpe,
            ref_acc = accept_p10)
g1 <- pit_res |>
  inner_join(ref_cells, by = c("target", "model_config", "horizon", "sample")) |>
  mutate(d_r2 = abs(r2_oos - ref_r2), d_cw = abs(cw - ref_cw), d_sh = abs(sharpe - ref_sh),
         acc_ok = accept_p10 == ref_acc)
stopifnot(nrow(g1) == nrow(pit_res),
          max(g1$d_r2, g1$d_cw, g1$d_sh) < 1e-6, all(g1$acc_ok))
cat(sprintf("    G1 PASS: %d PIT cells, max|Δ| r2/cw/sh = %.2e, accept matches\n",
            nrow(g1), max(g1$d_r2, g1$d_cw, g1$d_sh)))

# ---- 6. persist + headline -------------------------------------------------
tag      <- paste0("_", paste(SAMPLES, collapse = "-"), "_h", paste(HORIZONS, collapse = ""))
suffix   <- if (SMOKE) "_SMOKE" else ""
out_path <- file.path("data/audit", sprintf("leakage_payoff_%s%s%s.parquet", Sys.Date(), tag, suffix))
write_parquet(res, out_path)
cat(sprintf("[6] wrote %s (%d rows)\n", out_path, nrow(res)))

acc <- res |> group_by(panel, sample, horizon) |>
  summarise(n_cells = n(), n_accept = sum(accept_p10, na.rm = TRUE),
            mean_r2 = round(mean(r2_oos, na.rm = TRUE), 4),
            mean_sharpe = round(mean(sharpe, na.rm = TRUE), 3), .groups = "drop")
cat("\n=== acceptance: PIT (control) vs LEAKED (vintage-only) ===\n"); print(acc, n = Inf)
cat("\n=== per-cell PIT -> leaked shift ===\n")
print(res |> select(panel, target, sample, horizon, model_config, r2_oos, cw, sharpe, threshold, accept_p10) |>
        arrange(target, sample, horizon, model_config, panel), n = Inf)

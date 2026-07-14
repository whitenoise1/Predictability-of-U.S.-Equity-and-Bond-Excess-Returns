# scratch_pcr_props_setup.R
# -----------------------------------------------------------------------------
# Shared setup for the PCR (Research-method) shrinkage-property explorations.
# Mirrors scripts/p1_chapter_2.R steps 1-4 (panel -> bond targets -> S2 full
# pool -> transforms), then exposes build_cell(): the faithful per-as-of
# rectangular training block that the harness hands to C-PCR / the shrinkage
# class (walk_forward_fs lines 780-804). Caches the built panel/pool to RDS so
# the per-exploration scripts reuse it without re-running the heavy build.
#
# NOT production code — exploration scratch under working_paper/, per the
# scratch_*.R convention. Bond S2 only (the data-driven path the audit is about).
# -----------------------------------------------------------------------------
suppressMessages({library(arrow); library(dplyr); library(tibble); library(tidyr)})

.here <- function(p) p   # run from repo root (Rscript working_paper/scratch_*.R)
CACHE <- "working_paper/.scratch_pcr_setup_cache_v2.rds"   # v2: + corrected equity target
EQUITY_SP_COL <- "SP500TR"                                 # month-end TOTAL return (artifact-free)
EQUITY_EXCL   <- c("SP500TR", "SHILLER_PRICE")             # purge target price from equity pool

# ---- source the same chain chapter 2 uses --------------------------------
source("R/pit_query.R");  source("R/panel.R");   source("R/transforms.R")
source("R/targets.R");    source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")
invisible(capture.output(source("sources/factor_orthogonalization.R")))
source("R/feature_selection.R");  source("R/feature_pools.R")
source("R/u6_pool_prep.R")

# ---- locked constants (verbatim from p1_chapter_2.R) ---------------------
AS_OF             <- as.Date("2026-04-30")
W_Z               <- 60L
HORIZONS          <- c(1L, 3L, 12L)
T_START_S2_CUTOFF <- as.Date("1990-01-31")
MATURITY_GRID <- list(
  bond_1y  = list(y_col = "DGS1",  duration = 0.97),
  bond_2y  = list(y_col = "DGS2",  duration = 1.93),
  bond_3y  = list(y_col = "DGS3",  duration = 2.85),
  bond_5y  = list(y_col = "DGS5",  duration = 4.50),
  bond_10y = list(y_col = "DGS10", duration = 8.00))

build_setup <- function() {
  cat("[setup] building panel ...\n")
  meta <- read_parquet("data/metadata/factor_metadata.parquet")
  base_ids <- meta |>
    filter(public_private %in% c("PUBLIC", "MIXED"),
           !family %in% "FORECASTS",
           !grepl("^RTDSM_", factor_id)) |>
    pull(factor_id)
  all_factors <- unique(c(base_ids, "SHILLER_PRICE", "SP500TR", u6_extra_factor_ids(meta)))
  panel_levels <- build_panel(as_of_date = AS_OF, factor_ids = all_factors,
                              frequency = "M", fill = "forward") |>
    cp_factor_columns()

  cat("[setup] bond + corrected equity targets + S2 full pool ...\n")
  targets_bd <- bond_target_grid(panel_levels, HORIZONS, MATURITY_GRID)
  targets_eq <- setNames(lapply(HORIZONS, function(h)                       # SP500TR month-end
    construct_equity_excess_return(panel_levels, h = h, sp_col = EQUITY_SP_COL)),
    as.character(HORIZONS))
  dash <- derive_u6_dashboard(panel_levels); panel_levels <- dash$panel
  s2_full <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = NULL,
                                   forecast_token     = U6_FORECAST_TOKEN,
                                   dashboard_features = dash$dashboard_features)
  s2_deep <- build_scenario_2_pool(panel_levels, meta, t_start_cutoff = T_START_S2_CUTOFF,
                                   forecast_token     = U6_FORECAST_TOKEN,
                                   dashboard_features = dash$dashboard_features)

  cat("[setup] transforms on S2 pool features ...\n")
  ratio_feats <- s2_full$ratio_features
  level_feats <- s2_full$level_features
  panel_xfm <- panel_levels |>
    apply_transform(yoy_log_diff,   factor_ids = ratio_feats) |>
    apply_transform(rolling_zscore, factor_ids = ratio_feats, window = W_Z, min_obs = W_Z) |>
    apply_transform(rolling_zscore, factor_ids = level_feats, window = W_Z, min_obs = W_Z)

  S2_POOL_COLS <- c(s2_full$ratio_features, s2_full$level_features)        # full shrink pool
  deep_cohort  <- c(s2_deep$ratio_features, s2_deep$level_features)
  S2_DEEP_COLS <- live_cols(panel_levels, deep_cohort)                     # deep cohort (~1990)

  # Equity pool/anchor purge the target's own price (leakage); bonds keep full pool.
  S2_POOL_COLS_EQ <- setdiff(S2_POOL_COLS, EQUITY_EXCL)
  S2_DEEP_COLS_EQ <- setdiff(S2_DEEP_COLS, EQUITY_EXCL)
  cat(sprintf("[setup] panel %d x %d | S2 full = %d | deep cohort = %d live | equity pool = %d\n",
              nrow(panel_xfm), ncol(panel_xfm), length(S2_POOL_COLS), length(S2_DEEP_COLS),
              length(S2_POOL_COLS_EQ)))
  list(panel_xfm = panel_xfm, targets_bd = targets_bd, targets_eq = targets_eq,
       S2_POOL_COLS = S2_POOL_COLS, S2_DEEP_COLS = S2_DEEP_COLS,
       S2_POOL_COLS_EQ = S2_POOL_COLS_EQ, S2_DEEP_COLS_EQ = S2_DEEP_COLS_EQ)
}

if (file.exists(CACHE)) {
  cat("[setup] loading cache", CACHE, "\n"); .S <- readRDS(CACHE)
} else {
  .S <- build_setup(); saveRDS(.S, CACHE); cat("[setup] cached ->", CACHE, "\n")
}
panel_xfm       <- .S$panel_xfm
targets_bd      <- .S$targets_bd
targets_eq      <- .S$targets_eq
S2_POOL_COLS    <- .S$S2_POOL_COLS
S2_DEEP_COLS    <- .S$S2_DEEP_COLS
S2_POOL_COLS_EQ <- .S$S2_POOL_COLS_EQ
S2_DEEP_COLS_EQ <- .S$S2_DEEP_COLS_EQ

# ---- t_start helper (chapter 2: max first-non-NA over the anchor cols) ----
.first_non_na <- function(x, dates) { i <- which(!is.na(x)); if (!length(i)) as.Date(NA) else dates[i[1L]] }
t_start_for <- function(cols)
  max(as.Date(vapply(cols, function(c)
    as.character(.first_non_na(panel_xfm[[c]], panel_xfm$reference_date)), character(1))))

# ---- build_cell: faithful per-as-of training block for C-PCR / shrinkage ---
# Mirrors walk_forward_fs lines 745-807 for a SINGLE as-of T: per-as-of column
# availability (observed at x_T AND complete across the window) then complete
# cases -> the exact rectangular X_train pcr_core receives. T defaults to the
# last evaluable as-of (largest expanding block).
build_cell <- function(maturity = "bond_5y", h = 1L,
                       target = c("bond", "equity"),
                       pool_cols    = NULL,
                       t_start_cols = NULL,
                       window_kind  = c("expanding", "rolling"),
                       window_length = NULL, T_ref = NULL) {
  window_kind <- match.arg(window_kind); target <- match.arg(target)
  if (target == "equity") {
    y <- targets_eq[[as.character(h)]]
    if (is.null(pool_cols))    pool_cols    <- S2_POOL_COLS_EQ
    if (is.null(t_start_cols)) t_start_cols <- S2_DEEP_COLS_EQ
  } else {
    y <- targets_bd[[maturity]][[as.character(h)]]
    if (is.null(pool_cols))    pool_cols    <- S2_POOL_COLS
    if (is.null(t_start_cols)) t_start_cols <- S2_DEEP_COLS
  }
  ref <- panel_xfm$reference_date
  t_start <- t_start_for(t_start_cols)
  T_idx <- if (is.null(T_ref)) max(which(ref <= AS_OF & !is.na(y)))
           else which(ref == as.Date(T_ref))
  fit_hi <- T_idx - h
  fit_lo <- if (window_kind == "expanding") which(ref >= t_start)[1L]
            else max(fit_hi - window_length + 1L, which(ref >= t_start)[1L])
  stopifnot(!is.na(fit_lo), fit_lo <= fit_hi)
  rows <- fit_lo:fit_hi

  xT  <- panel_xfm[T_idx, pool_cols, drop = FALSE]
  win <- panel_xfm[rows,  pool_cols, drop = FALSE]
  col_ok <- vapply(pool_cols,
                   function(cc) !is.na(xT[[cc]]) && !anyNA(win[[cc]]), logical(1))
  cols <- pool_cols[col_ok]
  X_all <- as.matrix(panel_xfm[rows, cols, drop = FALSE]); y_t <- y[rows]
  keep  <- stats::complete.cases(X_all) & !is.na(y_t)

  list(X_train = X_all[keep, , drop = FALSE], y_train = y_t[keep],
       x_T = as.matrix(panel_xfm[T_idx, cols, drop = FALSE]),
       cols = cols, maturity = maturity, h = h,
       t_start = t_start, T_date = ref[T_idx],
       n_train = sum(keep), p = length(cols),
       n_pool = length(pool_cols), n_dropped_avail = length(pool_cols) - length(cols))
}

# ---- Supervised PCA (Bair et al. 2006) — shared by explore4 + battery --------
# screen raw features by univariate |t| vs y >= theta, PCA the screened block,
# regress on the LEADING m components (PIT train-only basis). NULL when <2 clear.
spca_core <- function(X_train, y_train, X_pred, theta, m) {
  n <- nrow(X_train)
  r  <- as.numeric(cor(X_train, y_train))
  tt <- abs(r) * sqrt((n - 2) / pmax(1 - r^2, 1e-12)); tt[!is.finite(tt)] <- 0
  keep <- which(tt >= theta); if (length(keep) < 2L) return(NULL)
  ctr <- colMeans(X_train)[keep]; scl <- apply(X_train[, keep, drop = FALSE], 2L, sd)
  scl[!is.finite(scl) | scl == 0] <- 1
  Zs <- sweep(sweep(X_train[, keep, drop = FALSE], 2L, ctr, "-"), 2L, scl, "/")
  sv <- svd(Zs); mm <- min(as.integer(m), ncol(Zs), n - 1L)
  V  <- sv$v[, seq_len(mm), drop = FALSE]; S <- Zs %*% V
  fit <- lm.fit(cbind(1, S), y_train); b <- fit$coefficients
  Zp  <- sweep(sweep(X_pred[, keep, drop = FALSE], 2L, ctr, "-"), 2L, scl, "/")
  beta_x <- numeric(ncol(X_train)); beta_x[keep] <- as.numeric(V %*% b[-1L])
  names(beta_x) <- colnames(X_train)
  list(y_hat = as.numeric(cbind(1, Zp %*% V) %*% b), n_screened = length(keep),
       m = mm, df = mm, keep = keep, V = V,
       var_share = sv$d[seq_len(mm)]^2 / sum(sv$d^2), beta_x = beta_x)
}
SPCA_THETA <- 1.645; SPCA_M <- 1:3
SPCA_CANDS <- lapply(SPCA_M, function(mm) list(theta = SPCA_THETA, m = as.integer(mm)))
spca_cv <- function(X, y, xT, h = 1L) {
  fp <- function(Xtr, ytr, xval) vapply(SPCA_CANDS, function(cd) {
    r <- spca_core(Xtr, ytr, xval, cd$theta, cd$m); if (is.null(r)) mean(ytr) else r$y_hat
  }, numeric(1))
  best <- pit_expanding_cv(X, y, h = h, seed = h + 2L, fold_predict = fp,
                           n_candidates = length(SPCA_CANDS))
  if (is.null(best)) return(NULL)
  cd <- SPCA_CANDS[[best]]; r <- spca_core(X, y, xT, cd$theta, cd$m)
  if (!is.null(r)) r$theta <- cd$theta
  r
}

cat("[setup] ready: build_cell(target=), spca_core/spca_cv, panel_xfm, targets_bd/eq, pools\n")

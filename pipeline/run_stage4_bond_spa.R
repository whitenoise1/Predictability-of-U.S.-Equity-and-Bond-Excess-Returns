# working_paper/run_stage4_bond_spa.R
# Hansen (2005) SPA on the S2_bond SHRINKAGE family, using the RETAINED bond
# forecasts produced by rerun_bond_forecasts.R. Mirrors run_stage4_spa.R but the
# SPA family is per (MATURITY x horizon x sample) — each bond tenor is its own
# series, exactly as the FDR per-maturity families are. Closes the Table-8 SPA
# "n/a" for bonds (they previously retained no pointwise forecasts).
#
# Gates: G1 — recompute (r2_oos, cw, sharpe) from forecasts and match the committed
# drop-it1 bond cells < 1e-8; no-NA stop guard; common axis >= 10 (else reason-tagged).
# RUN FROM REPO ROOT.
suppressMessages({library(arrow); library(dplyr); library(tidyr); library(purrr)})
source("R/predictive_metrics.R"); source("R/predictive_regression.R")
source("R/spa_test.R")

SPA_B  <- 2000L
G1_REF <- "data/audit/p10_rescored_drop-it1_2026-06-17.parquet"
S2_BOND_SHRINK <- c("C-ENET","C-RIDGE","C-PLS")
KEY <- c("chapter","scenario","maturity","model_config","iteration","horizon",
         "window_kind","window_length","sample")
`%||%` <- function(a, b) if (length(a)) a else b

# ---- locate the newest non-SMOKE bond forecast artifact --------------------
cands <- sort(list.files("data/audit", pattern = "^bond_s2_forecasts_[0-9-]+\\.parquet$",
                         full.names = TRUE), decreasing = TRUE)
if (!length(cands)) stop("no bond_s2_forecasts_<date>.parquet found — run rerun_bond_forecasts.R first")
fc <- read_parquet(cands[1]); cat(sprintf("[load] %s (%d rows)\n", cands[1], nrow(fc)))

# ---- per-cell recompute (G1 source) ----------------------------------------
fc <- fc |> arrange(across(all_of(KEY)), as_of_date)
inf <- fc |> group_by(across(all_of(KEY))) |> group_split() |> map_dfr(function(d) {
  h <- d$horizon[1L]
  d[1L, KEY] |> bind_cols(tibble(
    n_oos = nrow(d),
    r2_oos = r2_oos(d$y_hat, d$y_realized, d$y_bench),
    cw = clark_west_stat(clark_west_pointwise(d$y_hat, d$y_realized, d$y_bench), L = floor(1.5 * h))$stat,
    sharpe = signal_sharpe(d$y_hat, d$y_realized, h = h)))
})

# ---- G1: exact reproduction vs committed drop-it1 bond shrinkage cells ------
ref <- read_parquet(G1_REF) |>
  filter(chapter == "ch2", target == "bond", scenario == "S2_bond",
         model_config %in% S2_BOND_SHRINK, !is.na(accept_p10))
chk <- inf |> select(all_of(KEY), r2_oos, cw, sharpe) |>
  inner_join(ref |> select(all_of(KEY), r2_s = r2_oos, cw_s = cw, sh_s = sharpe), by = KEY)
cat(sprintf("[G1] inf cells=%d ; committed ref=%d ; matched=%d\n", nrow(inf), nrow(ref), nrow(chk)))
if (nrow(chk) != nrow(ref))
  stop(sprintf("[G1 FAIL] %d committed bond cells unmatched by retained forecasts", nrow(ref) - nrow(chk)))
g1 <- chk |> summarise(dr2 = max(abs(r2_oos - r2_s)), dcw = max(abs(cw - cw_s)), dsh = max(abs(sharpe - sh_s)))
cat(sprintf("[G1] max|dr2|=%.2e max|dcw|=%.2e max|dsharpe|=%.2e\n", g1$dr2, g1$dcw, g1$dsh))
if (!(g1$dr2 < 1e-8 && g1$dcw < 1e-8 && g1$dsh < 1e-8))
  stop("[G1 FAIL] bond forecasts do not reproduce committed scalars < 1e-8")
cat("[G1 OK] bond scalars reproduced < 1e-8\n")

# ---- SPA per (maturity x horizon x sample), pooling estimation windows ------
spa_family <- function(mat, h, smp, fc_sub) {
  sub <- fc_sub |> filter(maturity == mat, horizon == h, sample == smp)
  per <- sub |> group_by(across(all_of(KEY))) |>
    group_map(~ tibble(as_of_date = .x$as_of_date,
                       d = (.x$y_realized - .x$y_bench)^2 - (.x$y_realized - .x$y_hat)^2))
  row <- function(m_orig, m_used, n_common, n_dropped, o, reason) tibble(
    target = "bond", maturity = mat, horizon = h, sample = smp,
    m_orig = m_orig, m_used = m_used, n_common = n_common, n_dropped = n_dropped,
    q = 1 / h, B = SPA_B, T_spa = o$T_spa %||% NA_real_, best_z = o$best_z %||% NA_real_,
    p_lower = o$p_lower %||% NA_real_, p_consistent = o$p_consistent %||% NA_real_,
    p_upper = o$p_upper %||% NA_real_, reason = reason)
  if (length(per) < 1L) return(row(0L, 0L, 0L, 0L, list(), "no cells"))
  common <- sort(as.Date(reduce(map(per, ~ as.character(.x$as_of_date)), intersect)))
  if (length(common) < 10L)
    return(row(length(per), 0L, length(common), 0L, list(), sprintf("common axis < 10 (%d)", length(common))))
  mat_d <- vapply(per, function(p) p$d[match(common, p$as_of_date)], numeric(length(common)))
  sds   <- apply(mat_d, 2L, stats::sd)
  ndrop <- sum(sds == 0)
  mat_d <- mat_d[, sds > 0, drop = FALSE]
  if (ncol(mat_d) < 1L) return(row(length(per), 0L, length(common), ndrop, list(), "all columns benchmark-identical"))
  o <- spa_test(mat_d, q = 1 / h, B = SPA_B, seed = 11L)
  row(length(per), o$m, length(common), ndrop, o, "")
}
combos <- fc |> distinct(maturity, horizon, sample)
spa <- pmap_dfr(combos, function(maturity, horizon, sample) spa_family(maturity, horizon, sample, fc))

out <- file.path("data/audit", sprintf("stage4_spa_bond_%s.parquet", Sys.Date()))
write_parquet(spa, out)
cat(sprintf("[write] %s (%d families)\n", out, nrow(spa)))
cat("\n-- bond SPA (evaluable families, p_consistent = headline) --\n")
print(as.data.frame(spa |> filter(reason == "") |> arrange(p_consistent) |>
  select(maturity, horizon, sample, m_used, n_common, T_spa, p_consistent)), row.names = FALSE)
hits <- spa |> filter(reason == "", p_consistent < 0.05)
cat(sprintf("\nbond SPA families with p_consistent < 0.05: %d / %d evaluable\n",
            nrow(hits), sum(spa$reason == "")))
cat("[run_stage4_bond_spa complete]\n")

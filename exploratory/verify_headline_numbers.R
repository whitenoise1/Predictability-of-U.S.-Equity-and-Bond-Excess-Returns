# working_paper/verify_headline_numbers.R
#
# Independent verification of the equity headline numbers in
# chapter_2_3_research_note (AAII_BULLBEAR R2_OOS ~ 0.19, Sharpe ~ 1.43, and the
# "cluster of ten standalone equity predictors"). RUN FROM REPO ROOT.
#
# Design notes:
#  - Uses walk_forward_predict() DIRECTLY (not the walk_forward_fs wrapper that
#    predictor_value_scan.R uses) so this is an independent code path.
#  - Reproduces the per-factor transform by type (ratio = yoy_log_diff -> z;
#    level = trailing z only), matching the chapter-2/3 setup. The reconcile
#    step (Check 1) confirms the transforms match the stored scan to tolerance.
#  - Read-only: writes nothing to data/. Figures -> working_paper/figures/.

suppressMessages({
  library(arrow); library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})
source("R/pit_query.R"); source("R/panel.R"); source("R/transforms.R")
source("R/targets.R"); source("R/predictive_regression.R")
source("R/predictive_metrics.R"); source("R/walk_forward.R")

AS_OF  <- as.Date("2026-04-30")
W_Z    <- 60L
H      <- 1L
BURNIN_MONTHS <- 120L

# The 10 standalone accepters + 3 near-misses from PREDICTOR_VALUE_SUMMARY.md.
LEVEL_FACTORS <- c("AAII_BULLBEAR","SKEW","VIXCLS","ANFCI","NFCICREDIT","NFCI",
                   "NFCIRISK","NFCILEVERAGE","GACDFSA066MSFRBPHI",
                   "CORESTICKM159SFRBATL","AAA","ADS_INDEX")
RATIO_FACTORS <- c("CPIAUCSL")                       # only ratio among the 13
FACTORS <- c(LEVEL_FACTORS, RATIO_FACTORS)

# NBER recession months (post-1990) and the two acute-crisis windows.
nber_rec <- function(d) {
  (d >= "1990-07-01" & d <= "1991-03-31") |
  (d >= "2001-03-01" & d <= "2001-11-30") |
  (d >= "2007-12-01" & d <= "2009-06-30") |
  (d >= "2020-02-01" & d <= "2020-04-30")
}
is_gfc   <- function(d) d >= "2008-09-01" & d <= "2009-06-30"
is_covid <- function(d) d >= "2020-02-01" & d <= "2020-05-31"

# --- panel + transforms -----------------------------------------------------
cat("[setup] building panel ...\n")
need <- unique(c(FACTORS, "SHILLER_PRICE", "DGS3MO"))
panel_levels <- build_panel(as_of_date = AS_OF, factor_ids = need,
                            frequency = "M", fill = "forward")
panel_xfm <- panel_levels |>
  apply_transform(yoy_log_diff,   factor_ids = RATIO_FACTORS) |>
  apply_transform(rolling_zscore, factor_ids = RATIO_FACTORS, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = LEVEL_FACTORS, window = W_Z, min_obs = W_Z)
y_eq <- construct_equity_excess_return(panel_levels, h = H)   # PRICE-ONLY (Shiller price)

# --- helpers ----------------------------------------------------------------
run_wf <- function(fac) {
  ref <- panel_xfm$reference_date
  ts  <- ref[which(!is.na(panel_xfm[[fac]]))[1L]]
  te  <- seq(ts, by = "120 months", length.out = 2L)[2L]
  if (is.na(ts) || te >= AS_OF) return(NULL)
  pc <- panel_xfm |> select(reference_date, all_of(fac))
  wf <- walk_forward_predict(y = y_eq, panel = pc, horizon = H,
                             window_kind = "expanding",
                             t_start = ts, t_eval_start = te, t_end = AS_OF)
  wf
}

metrics <- function(wf, h = H) {
  if (is.null(wf) || nrow(wf) < 24L) return(NULL)
  f_T <- clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench)
  tibble(n = nrow(wf),
         r2 = r2_oos(wf$y_hat, wf$y_realized, wf$y_bench),
         cw = clark_west_stat(f_T, L = floor(1.5 * h))$stat,
         sharpe = signal_sharpe(wf$y_hat, wf$y_realized, h = h))
}

subsample_metrics <- function(wf, keep, h = H) {
  metrics(wf[keep, , drop = FALSE], h)
}

cssed <- function(wf) {
  e_b <- (wf$y_realized - wf$y_bench)^2
  e_m <- (wf$y_realized - wf$y_hat)^2
  tibble(as_of_date = wf$as_of_date, cssed = cumsum(e_b - e_m))
}

sharpe_decomp <- function(wf, h = H, cost = 0.001) {
  pos   <- sign(wf$y_hat)
  strat <- pos * wf$y_realized
  nz    <- pos != 0
  hit   <- mean(sign(wf$y_hat[nz]) == sign(wf$y_realized[nz]))
  ord   <- order(abs(strat), decreasing = TRUE)
  top5  <- sum(strat[ord][1:min(5, length(strat))]) / sum(strat)
  flips <- mean(abs(diff(pos)) > 0)
  dpos  <- c(abs(pos[1]), abs(diff(pos)))
  strat_net <- strat - cost * dpos
  net_sh <- mean(strat_net) / sd(strat_net) * sqrt(12 / h)
  tibble(hit_rate = hit, top5_share = top5, flip_rate = flips, sharpe_net10bp = net_sh)
}

is_r2 <- function(fac) {                                  # full-sample in-sample R^2
  d <- tibble(y = y_eq, x = panel_xfm[[fac]]) |> filter(!is.na(y), !is.na(x))
  if (nrow(d) < 24L) return(NA_real_)
  summary(lm(y ~ x, data = d))$r.squared
}

# --- run --------------------------------------------------------------------
cat("[run] walk-forward per factor ...\n")
wfs <- setNames(lapply(FACTORS, run_wf), FACTORS)
wfs <- wfs[!vapply(wfs, is.null, logical(1))]

full_tab <- bind_rows(lapply(names(wfs), function(f)
  metrics(wfs[[f]]) |> mutate(factor = f, .before = 1)))

# Check 1 — reconcile against the stored scan.
scan_path <- "data/audit/predictor_univariate_scan_2026-05-31.parquet"
recon <- NULL
if (file.exists(scan_path)) {
  scan <- read_parquet(scan_path) |> filter(target == "equity", horizon == 1L) |>
    select(factor, r2_scan = r2_oos, cw_scan = cw, sh_scan = sharpe, n_scan = n_oos)
  recon <- full_tab |> left_join(scan, by = "factor") |>
    mutate(dr2 = round(r2 - r2_scan, 4), dcw = round(cw - cw_scan, 3),
           dsh = round(sharpe - sh_scan, 3))
}

# Checks 2-4 — robustness decompositions.
robust <- bind_rows(lapply(names(wfs), function(f) {
  wf <- wfs[[f]]; d <- wf$as_of_date
  ex_gfc   <- subsample_metrics(wf, !is_gfc(d))
  ex_covid <- subsample_metrics(wf, !is_covid(d))
  ex_both  <- subsample_metrics(wf, !is_gfc(d) & !is_covid(d))
  exp_only <- subsample_metrics(wf, !nber_rec(d))
  sd <- sharpe_decomp(wf)
  tibble(factor = f,
         r2_full = full_tab$r2[full_tab$factor == f],
         r2_is   = is_r2(f),
         r2_exGFC = ex_gfc$r2, r2_exCOVID = ex_covid$r2, r2_exBoth = ex_both$r2,
         r2_expOnly = exp_only$r2,
         sh_full = full_tab$sharpe[full_tab$factor == f],
         sh_exBoth = ex_both$sharpe, sh_net10bp = sd$sharpe_net10bp,
         hit_rate = sd$hit_rate, top5_share = sd$top5_share, flip_rate = sd$flip_rate)
}))

# --- console report ---------------------------------------------------------
rnd <- function(df) df |> mutate(across(where(is.numeric), \(x) round(x, 3)))
cat("\n================ CHECK 1: reconcile vs stored scan ================\n")
if (!is.null(recon)) print(rnd(recon), n = Inf, width = Inf) else cat("scan parquet not found\n")
cat("\n================ CHECKS 2-4: robustness decomposition ================\n")
cat("r2_full = OOS R2 (full); r2_is = in-sample R2; ex* = OOS R2 dropping that window;\n")
cat("expOnly = NBER-expansions only; sh_* = signal Sharpe; flip_rate = sign-change freq.\n\n")
print(rnd(robust), n = Inf, width = Inf)

cat("\n================ effective OOS sample ================\n")
print(full_tab |> transmute(factor, n_oos = n, r2 = round(r2,3),
                            cw = round(cw,2), sharpe = round(sharpe,2)) |>
      arrange(desc(r2)), n = Inf)

# --- figures: CSSED staircase for the headline factors ----------------------
hl <- intersect(c("AAII_BULLBEAR","VIXCLS","NFCI","CPIAUCSL"), names(wfs))
cssed_df <- bind_rows(lapply(hl, function(f) cssed(wfs[[f]]) |> mutate(factor = f)))
p <- ggplot(cssed_df, aes(as_of_date, cssed)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  annotate("rect", xmin = as.Date("2008-09-01"), xmax = as.Date("2009-06-30"),
           ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "red") +
  annotate("rect", xmin = as.Date("2020-02-01"), xmax = as.Date("2020-05-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "red") +
  geom_line(linewidth = 0.5, colour = "#2c6fbb") +
  facet_wrap(~ factor, scales = "free_y") +
  labs(title = "Cumulative SSE difference (benchmark - model); rising = model beats hist. mean",
       subtitle = "Shaded: GFC (2008-09..2009-06) and COVID (2020-02..05). Staircase jumps at shading = crisis-concentrated.",
       x = NULL, y = "cumulative (e_bench^2 - e_model^2)") +
  theme_minimal(base_size = 10)
ggsave("working_paper/figures/cssed_headline_equity.pdf", p, width = 9, height = 6)
cat("\n[fig] working_paper/figures/cssed_headline_equity.pdf\n")
cat("[done]\n")

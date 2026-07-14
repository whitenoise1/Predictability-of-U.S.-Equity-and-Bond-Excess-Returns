#!/usr/bin/env Rscript
# scratch_prob_figures.R — PROTOTYPE figures for the three considerations.
# Read-only on artifacts; writes PNGs to working_paper/figures/scratch_*.png.
# Run from the PROJECT ROOT:  Rscript working_paper/scratch_prob_figures.R
#
# Produces:
#   scratch_fig1a_distributions.png  — #1: OOS R^2 and CW across all evaluable trials
#                                      (the near-null, shown; how hard forecasting is)
#   scratch_fig1b_forecast_bands.png — #1: forecast vs realized + error band for a
#                                      representative equity cell (thin wiggle in a
#                                      wide realized cloud)
#   scratch_fig3_attribution.png     — #3: native shrinkage weight, BOTH targets x
#                                      both estimators (revives old Tables 5/6 content)

suppressMessages({library(arrow); library(dplyr); library(tidyr); library(ggplot2)})
if (!exists("p10")) source("working_paper/explore_artifacts.R")

OUT  <- file.path(ROOT, "working_paper", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BLUE <- "#2c6fbb"; ORANGE <- "#e08214"; RED <- "#b2182b"

# tidytext-free "reorder within facet" helpers --------------------------------
reorder_within <- function(x, by, within) {
  key <- paste(x, within, sep = "\031")
  factor(key, levels = unique(key[order(within, by)]))
}
scale_x_reordered <- function(...) scale_x_discrete(labels = function(z) sub("\031.*$", "", z), ...)

# ============================================================================
# FIG 1a — distribution of the THREE acceptance criteria across evaluable trials
#   OOS R^2 (vs 0.01 floor) | Clark-West t (vs 1.645) | Sharpe vs target bar (vs 0)
# The Sharpe leg is target-relative: the bar is each target's own PIT-expanding
# median Sharpe (`threshold`), so we plot excess = sharpe - threshold (red = 0).
# ============================================================================
ev <- p10 %>% filter(!is.na(accept_p10)) %>%
  mutate(Target = recode(target, equity = "Equity", bond = "Bond"),
         exc_sr = sharpe - threshold)

L_R2 <- "OOS R-squared  (red = 0.01 floor)"
L_CW <- "Clark-West t-stat  (red = 1.645 bar)"
L_SR <- "Sharpe - target bar  (red = 0)"
# R^2 clipped at -0.5 (long left tail); CW and Sharpe full range
n_r2_below <- sum(ev$r2_oos < -0.5)
disp <- bind_rows(
  ev %>% filter(r2_oos >= -0.5) %>% transmute(Target, stat = L_R2, value = r2_oos),
  ev %>% transmute(Target, stat = L_CW, value = cw),
  ev %>% transmute(Target, stat = L_SR, value = exc_sr)) %>%
  mutate(stat = factor(stat, levels = c(L_R2, L_CW, L_SR)))
vl <- data.frame(stat = factor(c(L_R2, L_CW, L_SR), levels = c(L_R2, L_CW, L_SR)),
                 xint = c(0.01, 1.645, 0))

p1a <- ggplot(disp, aes(value, fill = Target)) +
  geom_histogram(bins = 70, alpha = 0.6, position = "identity", colour = NA) +
  geom_vline(data = vl, aes(xintercept = xint), colour = RED, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ stat, scales = "free", nrow = 1) +
  scale_fill_manual(values = c(Equity = BLUE, Bond = ORANGE)) +
  labs(x = NULL, y = "trials") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUT, "scratch_fig1a_distributions.png"), p1a, width = 11, height = 3.4, dpi = 150)
cat(sprintf("fig1a: %d trials | OOS R2 clipped at -0.5 (%d trials below, off-view) | CW range [%.2f, %.2f] | excess-Sharpe range [%.2f, %.2f]\n",
            nrow(ev), n_r2_below, min(ev$cw), max(ev$cw), min(ev$exc_sr,na.rm=TRUE), max(ev$exc_sr,na.rm=TRUE)))
cat(sprintf("       share R2>0.01=%.1f%% | share CW>1.645=%.1f%% | share Sharpe>bar=%.1f%%\n",
            100*mean(ev$r2_oos > 0.01), 100*mean(ev$cw > 1.645), 100*mean(ev$sharpe > ev$threshold, na.rm = TRUE)))

# ============================================================================
# FIG 1b — forecast vs realized with an error band, representative equity cell
# ============================================================================
cellkeys <- c("model_config","iteration","horizon","window_kind","window_length","sample")
pick <- eqfc %>% filter(chapter == "ch2") %>%
  count(across(all_of(cellkeys)), sort = TRUE) %>% slice(1)
cell <- eqfc %>% filter(chapter == "ch2") %>% semi_join(pick, by = cellkeys) %>% arrange(as_of_date)
rmse  <- sqrt(mean((cell$y_realized - cell$y_hat)^2))
r2oos <- 1 - sum((cell$y_realized - cell$y_hat)^2) / sum((cell$y_realized - cell$y_bench)^2)
lab <- sprintf("%s, h=%d, %s-%s, %s sample  |  OOS R2 = %.3f",
               pick$model_config, pick$horizon, pick$window_kind, pick$window_length,
               gsub("_", " ", pick$sample), r2oos)

p1b <- ggplot(cell, aes(as_of_date)) +
  geom_ribbon(aes(ymin = y_hat - rmse, ymax = y_hat + rmse), fill = BLUE, alpha = 0.18) +
  geom_point(aes(y = y_realized), colour = "grey35", size = 0.7, alpha = 0.8) +
  geom_line(aes(y = y_bench), colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  geom_line(aes(y = y_hat), colour = BLUE, linewidth = 0.6) +
  labs(x = NULL, y = "excess return (h-period)",
       subtitle = paste0("Forecast (blue) +/- RMSE band vs realized (grey); benchmark dotted.  ", lab)) +
  theme_minimal(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7), panel.grid.minor = element_blank())
ggsave(file.path(OUT, "scratch_fig1b_forecast_bands.png"), p1b, width = 9, height = 3.0, dpi = 150)
cat(sprintf("fig1b: %s | n=%d | RMSE=%.4f | realized sd=%.4f | band/scatter=%.2f\n",
            lab, nrow(cell), rmse, sd(cell$y_realized), rmse/sd(cell$y_realized)))

# ============================================================================
# FIG 3 — native shrinkage attribution, BOTH targets x both estimators
# ============================================================================
m3top <- m3 %>%
  group_by(target, model_config) %>%
  mutate(rel = 100 * mean_abs_coef / sum(mean_abs_coef)) %>%
  slice_max(rel, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(lab = paste(recode(target, equity = "Equity", bond = "Bond"), model_config, sep = " / "))

p3 <- ggplot(m3top, aes(x = reorder_within(variable, rel, lab), y = rel)) +
  geom_col(width = 0.7, fill = BLUE) +
  geom_text(aes(label = sprintf("%.1f", rel)), hjust = -0.2, size = 1.9) +
  coord_flip() +
  facet_wrap(~ lab, scales = "free_y", ncol = 2) +
  scale_x_reordered() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "Relative weight within fit (%)") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "scratch_fig3_attribution.png"), p3, width = 9, height = 5.4, dpi = 150)
cat("fig3: top-8 features per target x config (bond/equity x C-PLS/C-RIDGE)\n")

cat("\nWrote 3 PNGs to working_paper/figures/scratch_*.png\n")

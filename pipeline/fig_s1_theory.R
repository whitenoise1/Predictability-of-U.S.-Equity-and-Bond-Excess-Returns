#!/usr/bin/env Rscript
# SUPERSEDED 2026-07-13: Figure 5 is now drawn inline (vector) by the Rmd chunk
# `theory-bridge`, mirroring this script's panel grammar; this PNG generator is
# kept only for standalone use and is no longer called by build.sh.
# fig_s1_theory.R — Annex theory bridge (companion to fig:bridge).
# Predicted-vs-realized for the STRONGEST cell of each classic-literature theory
# pool, the Scenario-1 analogue of the data-driven bridge. Four panels (2x2):
#   equity  : Welch-Goyal u Chen-Roll-Ross   (+ surprise-augmented variant)
#   bond    : Fama-French 1989 term+default  /  Cochrane-Piazzesi forwards
# Each panel takes the single highest-pooled-OOS-R2 cell of that pool (>=24 OOS
# points) -- the best the canon manages anywhere in the sweep -- with the Mincer-
# Zarnowitz fit, the 45-degree perfect-forecast line, OOS R2 + n, and the bold P10
# acceptance verdict (every theory pool accepts ZERO under the rule). Mirrors the
# bridge .mkpanel grammar exactly. Reads the committed s1_forecasts artifact.
# RUN FROM REPO ROOT.
suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork)})

s1  <- read_parquet("data/audit/s1_forecasts_2026-06-20.parquet")
p10 <- read_parquet("data/audit/p10_rescored_drop-it1_2026-06-17.parquet")
keys <- c("scenario","target","maturity","model_config","horizon","window_kind","window_length","sample")
GRN <- "#1a9850"; RED <- "#b2182b"; BLUE <- "#2c6fbb"; OUT <- "working_paper/figures/fig_s1_theory.png"

# strongest cell (pooled OOS R2, n>=24) of one scenario; returns the forecast rows
best_cell <- function(scn) {
  sel <- s1 %>% filter(scenario == scn) %>% group_by(across(all_of(keys))) %>%
    summarise(n = n(), oos = 1 - sum((y_realized-y_hat)^2)/sum((y_realized-y_bench)^2),
              .groups = "drop") %>% filter(n >= 24) %>% arrange(desc(oos)) %>% slice(1)
  list(cell = s1 %>% semi_join(sel, by = keys), sel = sel)
}
# P10 accept verdict for the plotted cell (joins the re-scored acceptance artifact)
accept_of <- function(sel) {
  k <- intersect(keys, names(p10))
  v <- p10 %>% semi_join(sel, by = k) %>% pull(accept_p10)
  isTRUE(any(v))
}
pools <- list(
  list(scn = "S1_equity_WGCRR",     ttl = "Equity · Welch–Goyal ∪ Chen–Roll–Ross"),
  list(scn = "S1_equity_WGCRR_aug", ttl = "Equity · WG∪CRR + macro surprises"),
  list(scn = "S1_bond_FF1989",      ttl = "Bond · Fama–French (term + default)"),
  list(scn = "S1_bond_CP",          ttl = "Bond · Cochrane–Piazzesi (forwards)"))

mkpanel <- function(po) {
  bc <- best_cell(po$scn); cell <- bc$cell; sel <- bc$sel
  oos <- 1 - sum((cell$y_realized-cell$y_hat)^2)/sum((cell$y_realized-cell$y_bench)^2)
  co <- summary(lm(y_realized ~ y_hat, data = cell))$coefficients  # realized-on-forecast fit
  a <- co[1, 1]; b <- co[2, 1]; tb <- co[2, 3]                     # intercept, slope, slope t
  acc <- accept_of(sel)
  mat <- if (is.na(sel$maturity)) NULL else toupper(sub("bond_", "", sel$maturity))
  desc <- paste(c(sel$model_config, mat, paste0("h", sel$horizon),
                  if (sel$window_kind == "rolling") paste0("roll", sel$window_length) else "expand",
                  gsub("_", " ", sel$sample)), collapse = " · ")
  cx <- median(cell$y_hat);      rx <- quantile(abs(cell$y_hat - cx), 0.90)
  cy <- median(cell$y_realized); ry <- quantile(abs(cell$y_realized - cy), 0.90)
  qx <- cx + c(-rx, rx); qy <- cy + c(-ry, ry)
  ggplot(cell, aes(y_hat, y_realized)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, colour = RED, linetype = "dashed", linewidth = 0.5) +
    geom_point(colour = "grey35", size = 0.7, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, colour = BLUE, fill = BLUE, alpha = 0.15, linewidth = 0.6) +
    annotate("text", x = qx[1], y = qy[2], hjust = 0, vjust = 1, size = 2.6, family = "mono",
             colour = "grey15",
             label = sprintf("OOS R2 %.3f\nn  %d\na  %.4f\nb  %.2f\nt(b) %.2f", oos, nrow(cell), a, b, tb)) +
    annotate("text", x = qx[2], y = qy[2], hjust = 1, vjust = 1, size = 2.7, fontface = "bold",
             colour = ifelse(acc, GRN, RED), label = sprintf("%s", ifelse(acc, "ACCEPT", "REJECT"))) +
    coord_cartesian(xlim = qx, ylim = qy) +
    labs(x = "Predicted", y = "Realized", title = po$ttl, subtitle = desc) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), aspect.ratio = 1,
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, hjust = 0.5, colour = "grey15"))
}

ps <- lapply(pools, mkpanel)
# bonds on the top row, equities on the bottom
fig <- ((ps[[3]] | ps[[4]]) / (ps[[1]] | ps[[2]])) +
  plot_annotation(title = "Scenario 1 — academic theory",
                  theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))
ggsave(OUT, fig, width = 6.8, height = 6.7, dpi = 300)
cat("wrote", OUT, "\n")
# echo the plotted-cell summary for the caption / sign-off
for (po in pools) { bc <- best_cell(po$scn); s <- bc$sel
  cat(sprintf("  %-22s best: %s mat=%s h%d %s/%s %s  n=%d oos=%.3f accept=%s\n",
      po$scn, s$model_config, ifelse(is.na(s$maturity),"-",s$maturity), s$horizon,
      s$window_kind, ifelse(is.na(s$window_length),"-",s$window_length), s$sample,
      bc$cell %>% nrow(), 1 - sum((bc$cell$y_realized-bc$cell$y_hat)^2)/sum((bc$cell$y_realized-bc$cell$y_bench)^2),
      accept_of(bc$sel))) }

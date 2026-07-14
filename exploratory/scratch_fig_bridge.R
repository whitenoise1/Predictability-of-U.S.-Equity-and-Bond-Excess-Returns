#!/usr/bin/env Rscript
# scratch_fig_bridge.R — PROTOTYPE for the I->II forecasting-power bridge figure.
# Two Mincer-Zarnowitz panels on a common, equal-aspect scale:
#   left  = representative equity cell (the canonical near-null; forecast barely moves)
#   right = best equity cell by OOS R2 (a real-looking slope that dies under Stage-4 SPA)
# Also prints the pooled / population numbers the accompanying paragraph will cite.
# Reads committed equity_s2_forecasts + p10_rescored. SCRATCH; not wired into the .Rmd.

suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork)})
setwd(".")  # run from the parent-project root (the PIT database)
BLUE <- "#2c6fbb"; RED <- "#b2182b"

eqfc <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
p10  <- read_parquet("data/audit/p10_rescored_2026-06-12.parquet") |> filter(!is.na(accept_p10))
keys <- c("model_config","iteration","horizon","window_kind","window_length","sample")

mz_panel <- function(cell, ttl, sub) {
  fit <- lm(y_realized ~ y_hat, data = cell); co <- summary(fit)$coefficients
  oos <- 1 - sum((cell$y_realized-cell$y_hat)^2)/sum((cell$y_realized-cell$y_bench)^2)
  ann <- sprintf("beta = %.2f (t %.2f)\nR2 = %.3f\nOOS R2 = %.3f\nn = %d",
                 co[2,1], co[2,3], summary(fit)$r.squared, oos, nrow(cell))
  rng <- range(c(cell$y_hat, cell$y_realized))
  ggplot(cell, aes(y_hat, y_realized)) +
    geom_hline(yintercept=0, colour="grey85", linewidth=0.3) +
    geom_vline(xintercept=0, colour="grey85", linewidth=0.3) +
    geom_abline(slope=1, intercept=0, colour=RED, linetype="dashed", linewidth=0.5) +
    geom_point(colour="grey35", size=1.0, alpha=0.6) +
    geom_smooth(method="lm", formula=y~x, colour=BLUE, fill=BLUE, alpha=0.15, linewidth=0.7) +
    annotate("text", x=rng[1], y=rng[2], hjust=0, vjust=1, size=2.7, family="mono", label=ann) +
    coord_equal(xlim=rng, ylim=rng) +
    labs(x="Predicted excess return", y="Realized excess return", title=ttl, subtitle=sub) +
    theme_minimal(base_size=9) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=10),
          plot.subtitle=element_text(size=7, colour="grey35"))
}

# representative = most-populated equity-S2 cell
pick_r <- eqfc |> count(across(all_of(keys)), sort=TRUE) |> slice(1)
cell_r <- eqfc |> semi_join(pick_r, by=keys)
sub_r  <- sprintf("%s, h=%d, %s, %s", pick_r$model_config, pick_r$horizon, pick_r$window_kind, gsub("_"," ",pick_r$sample))

# best = highest OOS R2 among equity-S2 cells with n>=24
cs <- eqfc |> group_by(across(all_of(keys))) |>
  summarise(n=n(), oos=1-sum((y_realized-y_hat)^2)/sum((y_realized-y_bench)^2), .groups="drop") |>
  filter(n>=24) |> arrange(desc(oos))
pick_b <- cs |> slice(1)
cell_b <- eqfc |> semi_join(pick_b, by=keys)
sub_b  <- sprintf("%s, h=%d, %s-%s, %s", pick_b$model_config, pick_b$horizon, pick_b$window_kind, pick_b$window_length, gsub("_"," ",pick_b$sample))

pL <- mz_panel(cell_r, "Representative trial", sub_r)
pR <- mz_panel(cell_b, "Best trial (by OOS R-squared)", sub_b)
fig <- pL + pR +
  plot_annotation(
    title = "How hard is the forecast? Predicted vs realized excess returns",
    subtitle = "Blue = realized~forecast fit (95% CI); red dashed = perfect forecast (slope 1). Equal aspect: the forecast barely moves.",
    theme = theme(plot.title=element_text(face="bold", size=11), plot.subtitle=element_text(size=8, colour="grey35")))
ggsave("working_paper/figures/scratch_fig_bridge.png", fig, width=8.4, height=4.6, dpi=200)

# ---- numbers for the paragraph -------------------------------------------------
pooled <- eqfc  # all retained equity-S2 forecasts
oos_pooled <- 1 - sum((pooled$y_realized-pooled$y_hat)^2)/sum((pooled$y_realized-pooled$y_bench)^2)
ev <- p10 |> filter(target=="equity")
cat(sprintf("\nrepresentative: %s  | best: %s\n", sub_r, sub_b))
cat(sprintf("pooled equity OOS R2 (shape view, non-iid) = %.3f\n", oos_pooled))
cat(sprintf("equity trials: median OOS R2 = %.3f | %% beating 0.01 floor = %.1f%% | %% CW>1.645 = %.1f%%\n",
            median(ev$r2_oos,na.rm=TRUE), 100*mean(ev$r2_oos>0.01,na.rm=TRUE), 100*mean(ev$cw>1.645,na.rm=TRUE)))
allev <- p10
cat(sprintf("ALL trials (eq+bond): median OOS R2 = %.3f | %% beating 0.01 = %.1f%% | %% CW>1.645 = %.1f%% | n=%d\n",
            median(allev$r2_oos,na.rm=TRUE), 100*mean(allev$r2_oos>0.01,na.rm=TRUE),
            100*mean(allev$cw>1.645,na.rm=TRUE), nrow(allev)))
cat("wrote working_paper/figures/scratch_fig_bridge.png\n")

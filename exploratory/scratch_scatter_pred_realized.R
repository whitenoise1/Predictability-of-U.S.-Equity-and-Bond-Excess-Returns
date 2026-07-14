# scratch_scatter_pred_realized.R — EXPLORATORY (not for the paper yet).
# Realized vs predicted excess returns for the representative equity cell, with the
# Mincer-Zarnowitz regression realized = a + b*forecast (line + 95% CI), the 45-degree
# perfect-forecast reference, and the full regression summary (printed + annotated).
# Reads the committed equity_s2_forecasts artifact. Run from repo root.

suppressMessages({library(arrow); library(dplyr); library(ggplot2)})
BLUE <- "#2c6fbb"; RED <- "#b2182b"
OUT  <- "working_paper/figures/scratch_scatter_pred_realized.png"

eqfc <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
keys <- c("model_config","iteration","horizon","window_kind","window_length","sample")
pick <- eqfc |> count(across(all_of(keys)), sort = TRUE) |> slice(1)
cell <- eqfc |> semi_join(pick, by = keys)

fit <- lm(y_realized ~ y_hat, data = cell)
s   <- summary(fit)
co  <- s$coefficients
n   <- nrow(cell)
oos_r2 <- 1 - sum((cell$y_realized - cell$y_hat)^2) / sum((cell$y_realized - cell$y_bench)^2)
rmse   <- sqrt(mean((cell$y_realized - cell$y_hat)^2))

cat(sprintf("\nrepresentative cell: %s | %s | h=%d | %s-%s | %s sample | n=%d\n",
            pick$model_config, pick$iteration, pick$horizon, pick$window_kind,
            pick$window_length, pick$sample, n))
cat("\n=== Mincer-Zarnowitz: realized ~ forecast (OLS) ===\n"); print(co)
cat(sprintf("\nR^2 = %.4f | adj R^2 = %.4f | resid SE = %.4f on %d df | F = %.3f (p = %.3f)\n",
            s$r.squared, s$adj.r.squared, s$sigma, fit$df.residual,
            s$fstatistic[1], pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
cat(sprintf("OOS R^2 (vs benchmark mean) = %.4f | RMSE = %.4f | corr(forecast,realized) = %.3f\n",
            oos_r2, rmse, cor(cell$y_hat, cell$y_realized)))
cat(sprintf("perfect-forecast target: intercept=0, slope=1.  Joint? estimated (a,b)=(%.4f, %.4f)\n",
            co[1,1], co[2,1]))

ann <- sprintf(paste0(
  "alpha = %.4f  (SE %.4f, t %.2f, p %.2f)\n",
  "beta  = %.4f  (SE %.4f, t %.2f, p %.2f)\n",
  "R^2 = %.3f   adj R^2 = %.3f   n = %d\n",
  "resid SE = %.4f   corr = %.3f\n",
  "OOS R^2 = %.3f   RMSE = %.4f"),
  co[1,1], co[1,2], co[1,3], co[1,4],
  co[2,1], co[2,2], co[2,3], co[2,4],
  s$r.squared, s$adj.r.squared, n, s$sigma, cor(cell$y_hat, cell$y_realized), oos_r2, rmse)

rng <- range(c(cell$y_hat, cell$y_realized))
p <- ggplot(cell, aes(y_hat, y_realized)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = RED, linetype = "dashed", linewidth = 0.5) +
  geom_point(colour = "grey35", size = 1.1, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, colour = BLUE, fill = BLUE, alpha = 0.15, linewidth = 0.7) +
  annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1, size = 2.9, family = "mono", label = ann) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "Predicted h-period excess return", y = "Realized h-period excess return",
       title = "Realized vs predicted (representative equity cell)",
       subtitle = sprintf("%s, h=%d, %s, %s sample.  Blue = OLS realized~forecast (95%% CI); red dashed = perfect forecast (slope 1).",
                          pick$model_config, pick$horizon, pick$window_kind, gsub("_", " ", pick$sample))) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey35"))
ggsave(OUT, p, width = 6.2, height = 6.4, dpi = 160)
cat(sprintf("\nwrote %s\n", OUT))

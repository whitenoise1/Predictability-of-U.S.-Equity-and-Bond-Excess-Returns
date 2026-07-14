# scratch_error_dist.R — EXPLORATORY. Distribution of the forecasting error
# (realized - predicted) for the representative equity cell, overlaid against the
# benchmark error (realized - training mean) — does the forecast tighten the errors?
# Reads committed equity_s2_forecasts. Run from repo root.

suppressMessages({library(arrow); library(dplyr); library(tidyr); library(ggplot2)})
BLUE <- "#2c6fbb"; GREY <- "#8c8c8c"

e <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
k <- c("model_config","iteration","horizon","window_kind","window_length","sample")
pick <- e |> count(across(all_of(k)), sort = TRUE) |> slice(1)
cell <- e |> semi_join(pick, by = k)

cell <- cell |> mutate(err_model = y_realized - y_hat, err_bench = y_realized - y_bench)
rmse_m <- sqrt(mean(cell$err_model^2)); rmse_b <- sqrt(mean(cell$err_bench^2))
oos_r2 <- 1 - sum(cell$err_model^2) / sum(cell$err_bench^2)

cat(sprintf("\nrepresentative cell: %s | h=%d | %s | %s sample | n=%d\n",
            pick$model_config, pick$horizon, pick$window_kind, pick$sample, nrow(cell)))
cat(sprintf("forecast error : mean(bias)=%.4f  sd=%.4f  RMSE=%.4f\n",
            mean(cell$err_model), sd(cell$err_model), rmse_m))
cat(sprintf("benchmark error: mean(bias)=%.4f  sd=%.4f  RMSE=%.4f\n",
            mean(cell$err_bench), sd(cell$err_bench), rmse_b))
cat(sprintf("RMSE ratio (model/bench)=%.4f  =>  OOS R^2=%.4f  (error reduction = %.2f%%)\n",
            rmse_m/rmse_b, oos_r2, 100*(1 - rmse_m/rmse_b)))

long <- cell |>
  select(err_model, err_bench) |>
  pivot_longer(everything(), names_to = "kind", values_to = "err") |>
  mutate(kind = recode(kind, err_model = "Forecast error (realized - prediction)",
                             err_bench = "Benchmark error (realized - training mean)"))
ann <- sprintf(paste0("Forecast:  bias %.4f,  RMSE %.4f\n",
                      "Benchmark: bias %.4f,  RMSE %.4f\n",
                      "RMSE ratio %.3f  ->  OOS R^2 %.3f"),
               mean(cell$err_model), rmse_m, mean(cell$err_bench), rmse_b, rmse_m/rmse_b, oos_r2)

p <- ggplot(long, aes(err, fill = kind, colour = kind)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.2, size = 2.9, family = "mono", label = ann) +
  scale_fill_manual(values = c("Forecast error (realized - prediction)" = BLUE,
                               "Benchmark error (realized - training mean)" = GREY)) +
  scale_colour_manual(values = c("Forecast error (realized - prediction)" = BLUE,
                                 "Benchmark error (realized - training mean)" = GREY)) +
  labs(x = "Error (h-period excess return)", y = "density", fill = NULL, colour = NULL,
       title = "Distribution of the forecasting error (representative equity cell)",
       subtitle = sprintf("%s, h=%d, %s, %s sample, n=%d.  Forecast error vs benchmark error nearly coincide -> no precision gain.",
                          pick$model_config, pick$horizon, pick$window_kind, gsub("_"," ",pick$sample), nrow(cell))) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 7, colour = "grey35"))
ggsave("working_paper/figures/scratch_error_dist.png", p, width = 7.2, height = 4.0, dpi = 160)
cat("\nwrote working_paper/figures/scratch_error_dist.png\n")

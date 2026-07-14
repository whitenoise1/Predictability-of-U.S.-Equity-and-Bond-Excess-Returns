# scratch_error_dist_pooled.R — EXPLORATORY. Forecasting-error distribution POOLED
# across ALL equity ch2 retained forecasts (every cell x as-of), forecast error vs
# benchmark error. NB: pools across cells so the same month recurs under many
# configs/windows -> n is non-independent; this is a shape/aggregate view, not a test.
# Reads committed equity_s2_forecasts. Run from repo root.

suppressMessages({library(arrow); library(dplyr); library(tidyr); library(ggplot2)})
BLUE <- "#2c6fbb"; GREY <- "#8c8c8c"

e <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2") |>
  mutate(err_model = y_realized - y_hat, err_bench = y_realized - y_bench)
k <- c("model_config","iteration","horizon","window_kind","window_length","sample")
n_cells <- nrow(distinct(e, across(all_of(k))))

rmse_m <- sqrt(mean(e$err_model^2)); rmse_b <- sqrt(mean(e$err_bench^2))
oos_r2 <- 1 - sum(e$err_model^2) / sum(e$err_bench^2)
cat(sprintf("\nPOOLED equity ch2: %d forecast rows across %d cells\n", nrow(e), n_cells))
cat(sprintf("forecast error : bias=%.4f sd=%.4f RMSE=%.4f\n", mean(e$err_model), sd(e$err_model), rmse_m))
cat(sprintf("benchmark error: bias=%.4f sd=%.4f RMSE=%.4f\n", mean(e$err_bench), sd(e$err_bench), rmse_b))
cat(sprintf("RMSE ratio (model/bench)=%.4f => pooled OOS R^2=%.4f (error reduction %.2f%%)\n",
            rmse_m/rmse_b, oos_r2, 100*(1 - rmse_m/rmse_b)))
cat(sprintf("by horizon: %s\n", paste(sprintf("h=%d r=%.3f",
  sort(unique(e$horizon)),
  sapply(sort(unique(e$horizon)), function(h){ s<-e[e$horizon==h,]; 1 - sum(s$err_model^2)/sum(s$err_bench^2) })),
  collapse="  ")))

long <- e |> select(err_model, err_bench) |>
  pivot_longer(everything(), names_to = "kind", values_to = "err") |>
  mutate(kind = recode(kind, err_model = "Forecast error (realized - prediction)",
                             err_bench = "Benchmark error (realized - training mean)"))
ann <- sprintf(paste0("Forecast:  bias %.4f,  RMSE %.4f\n",
                      "Benchmark: bias %.4f,  RMSE %.4f\n",
                      "RMSE ratio %.4f  ->  pooled OOS R^2 %.4f"),
               mean(e$err_model), rmse_m, mean(e$err_bench), rmse_b, rmse_m/rmse_b, oos_r2)

p <- ggplot(long, aes(err, fill = kind, colour = kind)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.2, size = 2.9, family = "mono", label = ann) +
  scale_fill_manual(values = c("Forecast error (realized - prediction)" = BLUE,
                               "Benchmark error (realized - training mean)" = GREY)) +
  scale_colour_manual(values = c("Forecast error (realized - prediction)" = BLUE,
                                 "Benchmark error (realized - training mean)" = GREY)) +
  labs(x = "Error (h-period excess return)", y = "density", fill = NULL, colour = NULL,
       title = "Forecasting-error distribution — pooled across all equity trials",
       subtitle = sprintf("%s forecast rows across %d Scenario-2 equity cells (ch2). Forecast vs benchmark error coincide -> no aggregate precision gain.",
                          format(nrow(e), big.mark=","), n_cells)) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 7, colour = "grey35"))
ggsave("working_paper/figures/scratch_error_dist_pooled.png", p, width = 7.2, height = 4.0, dpi = 160)
cat("\nwrote working_paper/figures/scratch_error_dist_pooled.png\n")

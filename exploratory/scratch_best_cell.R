# scratch_best_cell.R — EXPLORATORY. Find the best (horizon, trial, estimator) by OOS
# R^2 and show the realized-vs-predicted Mincer-Zarnowitz scatter for the best EQUITY
# cell (bond has no retained pointwise forecasts -> scalars only, can't scatter).
# Reads committed p10_rescored + equity_s2_forecasts. Run from repo root.

suppressMessages({library(arrow); library(dplyr); library(ggplot2)})
BLUE <- "#2c6fbb"; RED <- "#b2182b"

# ---- 1. best combinations across the whole trial universe (p10) ------------------
p10 <- read_parquet("data/audit/p10_rescored_2026-06-12.parquet") |> filter(!is.na(accept_p10))
kp <- c("chapter","scenario","target","maturity","model_config","iteration",
        "horizon","window_kind","window_length","sample")
cat("=== top 12 cells by OOS R^2 (all evaluable trials) ===\n")
p10 |> arrange(desc(r2_oos)) |>
  transmute(target, maturity, model_config, iteration, h = horizon,
            win = ifelse(window_kind=="rolling", paste0("roll", window_length), "exp"),
            sample, r2_oos = round(r2_oos,3), cw = round(cw,2), sharpe = round(sharpe,2),
            accept = accept_p10) |>
  head(12) |> as.data.frame() |> print(row.names = FALSE)
best <- p10 |> arrange(desc(r2_oos)) |> slice(1)
cat(sprintf("\nBEST OVERALL: %s %s | %s | %s | h=%d | %s%s | %s | OOS R2=%.3f CW=%.2f Sharpe=%.2f (accept=%d)\n",
            best$target, ifelse(is.na(best$maturity),"",best$maturity), best$scenario,
            best$model_config, best$horizon, best$window_kind,
            ifelse(best$window_kind=="rolling", best$window_length, ""), best$sample,
            best$r2_oos, best$cw, best$sharpe, best$accept_p10))

# ---- 2. best EQUITY cell that has retained forecasts (scatter-able) ---------------
eqfc <- read_parquet("data/audit/equity_s2_forecasts_2026-06-12.parquet") |> filter(chapter == "ch2")
keys <- c("model_config","iteration","horizon","window_kind","window_length","sample")
cellstats <- eqfc |> group_by(across(all_of(keys))) |>
  summarise(n = n(),
            oos_r2 = 1 - sum((y_realized-y_hat)^2)/sum((y_realized-y_bench)^2),
            .groups = "drop") |>
  filter(n >= 24) |> arrange(desc(oos_r2))
cat("\n=== top 8 EQUITY S2 cells by OOS R^2 (retained forecasts) ===\n")
cellstats |> mutate(oos_r2 = round(oos_r2,3),
                    win = ifelse(window_kind=="rolling", paste0("roll",window_length),"exp")) |>
  transmute(model_config, iteration, h=horizon, win, sample, n, oos_r2) |>
  head(8) |> as.data.frame() |> print(row.names = FALSE)
pick <- cellstats |> slice(1)
cell <- eqfc |> semi_join(pick, by = keys)

fit <- lm(y_realized ~ y_hat, data = cell); s <- summary(fit); co <- s$coefficients; n <- nrow(cell)
oos_r2 <- 1 - sum((cell$y_realized-cell$y_hat)^2)/sum((cell$y_realized-cell$y_bench)^2)
rmse <- sqrt(mean((cell$y_realized-cell$y_hat)^2))
cat(sprintf("\nBEST EQUITY (scatter): %s | %s | h=%d | %s-%s | %s | n=%d\n",
            pick$model_config, pick$iteration, pick$horizon, pick$window_kind,
            pick$window_length, pick$sample, n))
cat("\n=== Mincer-Zarnowitz: realized ~ forecast (OLS) ===\n"); print(co)
cat(sprintf("R^2=%.4f adjR^2=%.4f residSE=%.4f F=%.3f (p=%.3f) | OOS R2=%.4f RMSE=%.4f corr=%.3f\n",
            s$r.squared, s$adj.r.squared, s$sigma, s$fstatistic[1],
            pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail=FALSE),
            oos_r2, rmse, cor(cell$y_hat, cell$y_realized)))

ann <- sprintf(paste0("alpha = %.4f  (t %.2f, p %.2f)\nbeta  = %.4f  (t %.2f, p %.2f)\n",
                      "R^2 = %.3f   n = %d\nOOS R^2 = %.3f   RMSE = %.4f   corr = %.3f"),
               co[1,1], co[1,3], co[1,4], co[2,1], co[2,3], co[2,4],
               s$r.squared, n, oos_r2, rmse, cor(cell$y_hat, cell$y_realized))
rng <- range(c(cell$y_hat, cell$y_realized))
p <- ggplot(cell, aes(y_hat, y_realized)) +
  geom_hline(yintercept = 0, colour="grey85", linewidth=0.3) +
  geom_vline(xintercept = 0, colour="grey85", linewidth=0.3) +
  geom_abline(slope=1, intercept=0, colour=RED, linetype="dashed", linewidth=0.5) +
  geom_point(colour="grey35", size=1.1, alpha=0.7) +
  geom_smooth(method="lm", formula=y~x, colour=BLUE, fill=BLUE, alpha=0.15, linewidth=0.7) +
  annotate("text", x=rng[1], y=rng[2], hjust=0, vjust=1, size=2.9, family="mono", label=ann) +
  coord_equal(xlim=rng, ylim=rng) +
  labs(x="Predicted h-period excess return", y="Realized h-period excess return",
       title="Realized vs predicted — BEST equity cell by OOS R-squared",
       subtitle=sprintf("%s, h=%d, %s-%s, %s sample.  Blue = OLS realized~forecast (95%% CI); red dashed = perfect forecast.",
                        pick$model_config, pick$horizon, pick$window_kind, pick$window_length, gsub("_"," ",pick$sample))) +
  theme_minimal(base_size=9) +
  theme(panel.grid.minor=element_blank(), plot.subtitle=element_text(size=7, colour="grey35"))
ggsave("working_paper/figures/scratch_scatter_best_equity.png", p, width=6.2, height=6.4, dpi=160)
cat("\nwrote working_paper/figures/scratch_scatter_best_equity.png\n")

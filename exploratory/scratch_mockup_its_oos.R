# MOCKUP ONLY — synthetic/illustrative data, NOT results.
# Two geometries for the combined a-priori(ITS) -> OOS family-signal exhibit.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

set.seed(1)
fams <- c("Yield curve","Policy uncertainty","SPF forecasts","Monetary policy / rates",
          "Inflation","Growth / activity","Equity volatility","Credit spreads",
          "Financial conditions","Cross-asset (pract.)","FX","Housing","Commodities",
          "Equity valuation","Labor market","Equity sentiment","Consumption / wealth")

# Synthetic in-sample marginal R% (concentrated) and PIT-OOS marginal R% (collapses, some <0)
mk <- function(its_hi, seed){
  set.seed(seed)
  its <- sort(runif(length(fams), 0.4, its_hi), decreasing = TRUE)
  oos <- pmax(-1.6, its*runif(length(fams), -0.15, 0.28) + rnorm(length(fams), 0, 0.25))
  tibble(family = fams, its = its, oos = oos)
}
dat <- bind_rows(bond = mk(7.5, 11), equity = mk(5.0, 22), .id = "target") |>
  mutate(target = factor(target, c("bond","equity"), c("Bond","Equity")))

ord <- dat |> group_by(family) |> summarise(m = mean(its)) |> arrange(m) |> pull(family)
dat$family <- factor(dat$family, levels = ord)

# ---- (1) DUMBBELL / ARROW ----
long <- dat |> pivot_longer(c(its, oos), names_to = "leg", values_to = "r2") |>
  mutate(leg = factor(leg, c("its","oos"), c("In-sample (a-priori)","PIT out-of-sample")))

g1 <- ggplot(dat) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(y = family, yend = family, x = its, xend = oos),
               colour = "grey60",
               arrow = arrow(length = unit(0.07, "inches"), type = "closed")) +
  geom_point(data = long, aes(x = r2, y = family, colour = leg), size = 2.4) +
  facet_wrap(~target) +
  scale_colour_manual(values = c("In-sample (a-priori)" = "#2c6fbb",
                                 "PIT out-of-sample" = "#c0392b"), name = NULL) +
  labs(x = "Marginal R² (%)", y = NULL,
       title = "Where the a-priori signal lives — and what survives out of sample",
       subtitle = "Arrow = in-sample → PIT-OOS decay. Dashed line = survival threshold (OOS R² = 0).",
       caption = "MOCKUP — synthetic illustrative data, not results.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.caption = element_text(colour = "grey50"))

ggsave("figures/scratch_mockup_dumbbell.png", g1, width = 10, height = 6.2, dpi = 150)

# ---- (2) ITS-vs-OOS SCATTER ----
g2 <- ggplot(dat, aes(its, oos)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = target), size = 2.6) +
  geom_text(aes(label = family), size = 2.5, hjust = 0, nudge_x = 0.07,
            check_overlap = TRUE, colour = "grey25") +
  facet_wrap(~target) +
  scale_colour_manual(values = c("Bond" = "#2c6fbb", "Equity" = "#c0392b"), guide = "none") +
  coord_cartesian(xlim = c(0, 8.2)) +
  labs(x = "In-sample marginal R² (%)  — a-priori",
       y = "PIT out-of-sample marginal R² (%)",
       title = "In-sample vs out-of-sample family signal",
       subtitle = "45° line = no decay; points hug y ≈ 0 (or below) = the a-priori signal evaporates.",
       caption = "MOCKUP — synthetic illustrative data, not results.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(colour = "grey50"))

ggsave("figures/scratch_mockup_scatter.png", g2, width = 10, height = 5.6, dpi = 150)
cat("done\n")

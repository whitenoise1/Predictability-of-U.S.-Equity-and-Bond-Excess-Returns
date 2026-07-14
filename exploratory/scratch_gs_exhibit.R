#!/usr/bin/env Rscript
# scratch_gs_exhibit.R — PROTOTYPE for the Part-I Gram-Schmidt admission exhibit (S2).
# SCRATCH / not wired into the .Rmd. Reads committed gs_selection_freq (admission
# frequency of each variable across the PCA-GS screen) and shows the S2 full-database
# screen behaviour: how often, run fresh each trial, the screen re-selects each variable.
# Admission freq is identical across C-OLS/C-COMB/C-CT (one shared screen) -> no proc axis.
# Two candidate looks for sign-off:
#   1  horizontal bar / lollipop, faceted bond | equity
#   2  compact colour heatmap-table companion to Figure 5 (same grammar)
# Labels reuse the Figure-5 annex scheme: <family>-<Annex A.1 row no.>.

suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork)})
setwd(".")  # run from the parent-project root (the PIT database)

g  <- read_parquet("data/audit/gs_selection_freq_2026-04-30.parquet")
uv <- read_parquet("data/audit/annex_used_variables_2026-06-14.parquet")
pal <- c("#eef3f7", "#cfe0ec", "#9dc0db", "#5b91c0", "#2c6fbb", "#1b4a86")

FAM_ABBR <- c(MONETARY_POLICY_RATES="RATE", YIELD_CURVE="CURV", CREDIT_SPREADS="CRED",
  INFLATION="INFL", GROWTH="GROW", LABOR_MARKET="LABR", FORECASTS="SPF",
  EQUITY_VOLATILITY="EVOL", EQUITY_SENTIMENT="ESNT", EQUITY_VALUATION="EVAL",
  FINANCIAL_CONDITIONS_COMPOSITE="FCI", POLICY_UNCERTAINTY_GEOPOLITICAL="POLU",
  FX_CURRENCY="FX", COMMODITIES="COMM", HOUSING="HOUS", CONSUMPTION_WEALTH="CONW",
  CROSS_ASSET="XAST", DASHBOARD_DERIVED="DASH")
lab_map <- uv %>% transmute(variable = code, ann_fam = family, ann_no = No)

# S2 (full-database) screen only; freq is config-invariant -> take C-OLS rows
d <- g %>% filter(grain == "by_scenario", scenario %in% c("S2_bond","S2_equity"),
                  model_config == "C-OLS", freq > 0) %>%
  left_join(lab_map, by = "variable") %>%
  mutate(ann_fam = coalesce(ann_fam, family),
         label = ifelse(is.na(ann_no),
                        ifelse(is.na(ann_fam), variable, paste0(FAM_ABBR[ann_fam], "-?")),
                        paste0(FAM_ABBR[ann_fam], "-", ann_no)),
         Target = recode(target, bond = "Bond", equity = "Equity"),
         pct = 100 * freq)

cat("S2 admission profile (label / variable / pct):\n")
print(d %>% arrange(Target, desc(pct)) %>% select(Target, label, variable, pct) %>% as.data.frame())

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 12))

# per-facet ordering: disambiguate shared labels across targets, strip suffix on axis
d <- d %>% mutate(ord_label = paste(Target, label, sep = "___"))
ord_levels <- d %>% arrange(Target, pct) %>% pull(ord_label)

# ---- Look 1: horizontal bars, faceted -----------------------------------------
p1 <- ggplot(d, aes(x = factor(ord_label, levels = ord_levels), y = pct)) +
  geom_col(width = 0.65, fill = "#2c6fbb") +
  geom_text(aes(label = sprintf("%.0f", pct)), hjust = -0.25, size = 2.8) +
  coord_flip() +
  facet_wrap(~ Target, scales = "free_y") +
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_y_continuous(limits = c(0, 112), expand = c(0,0)) +
  labs(x = NULL, y = "Admission frequency (% of S2 trials the GS screen selected it)",
       title = "Gram-Schmidt admission frequency (S2 full-database screen)") +
  base_theme + theme(panel.grid.major.y = element_blank())
ggsave("working_paper/figures/scratch_gs_bars.png", p1, width = 8.6, height = 4.2, dpi = 200)

# ---- Look 2: colour heatmap-table companion to Figure 5 ------------------------
tile_panel <- function(df, ttl) {
  df <- df %>% mutate(label = factor(label, levels = label[order(pct)]))
  ggplot(df, aes(x = "S2 screen", y = label, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = pct >= 55), size = 2.8) +
    scale_fill_stepsn(colours = pal, breaks = c(5,10,25,50,75), limits = c(0,100),
                      name = "admit %", guide = guide_coloursteps(barheight = 4, barwidth = 0.6)) +
    scale_colour_manual(values = c(`TRUE`="white", `FALSE`="grey25"), guide = "none") +
    labs(x = NULL, y = NULL, title = ttl) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8), plot.title = element_text(face="bold", size=11))
}
p2 <- tile_panel(filter(d, Target=="Equity"), "Equity") +
      tile_panel(filter(d, Target=="Bond"), "Bond") +
  plot_annotation(title = "Gram-Schmidt admission frequency (S2 screen) — heatmap-table look",
                  subtitle = "Share of S2 trials in which the screen, refit each trial, selected the variable. Labels: <family>-<Annex A.1 row no.>.")
ggsave("working_paper/figures/scratch_gs_heat.png", p2, width = 7.0, height = 4.6, dpi = 200)

cat("\nwrote scratch_gs_bars.png  +  scratch_gs_heat.png\n")

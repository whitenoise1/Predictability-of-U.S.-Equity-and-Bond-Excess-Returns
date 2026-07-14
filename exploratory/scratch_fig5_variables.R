#!/usr/bin/env Rscript
# scratch_fig5_variables.R — PROTOTYPE for Figure 5 (variable-level attribution).
# SCRATCH / not wired into the .Rmd. Reads committed attribution_all (variable-level
# effective-importance shares; the same source Figure 4 aggregates to families) and
# renders three density variants for sign-off:
#   A  true-share colour            (honest; shrinkage columns go pale = the blank-column risk)
#   B  per-column-normalised colour (each procedure's leader saturates; number = TRUE share)
#   C  concentrated-procedures-only fallback (drop the diffuse C-RIDGE/C-PLS columns)
# Grammar (palette/theme) matches the committed Figure 4 exactly.

suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork)})
setwd(".")  # run from the parent-project root (the PIT database)

al <- read_parquet("data/audit/attribution_all_2026-04-30.parquet")
uv <- read_parquet("data/audit/annex_used_variables_2026-06-14.parquet")
PROC_ORDER <- c("C-CP","C-OLS","C-COMB","C-CT","C-ENET","C-RIDGE","C-PLS")
pal <- c("#eef3f7", "#cfe0ec", "#9dc0db", "#5b91c0", "#2c6fbb", "#1b4a86")
K <- 6L

# --- annex-traceable label: <FAM-abbr>-<annex No> (forward rates: curve fallback) ---
FAM_ABBR <- c(MONETARY_POLICY_RATES="RATE", YIELD_CURVE="CURV", CREDIT_SPREADS="CRED",
  INFLATION="INFL", GROWTH="GROW", LABOR_MARKET="LABR", FORECASTS="SPF",
  EQUITY_VOLATILITY="EVOL", EQUITY_SENTIMENT="ESNT", EQUITY_VALUATION="EVAL",
  FINANCIAL_CONDITIONS_COMPOSITE="FCI", POLICY_UNCERTAINTY_GEOPOLITICAL="POLU",
  FX_CURRENCY="FX", COMMODITIES="COMM", HOUSING="HOUS", CONSUMPTION_WEALTH="CONW",
  CROSS_ASSET="XAST", DASHBOARD_DERIVED="DASH")
lab_map <- uv %>% transmute(variable = code, ann_fam = family, ann_no = No)
fwd_lab <- function(v) paste0("CURV-", sub("^FWD_", "f", v))
mk_label <- function(variable, ann_fam, ann_no) ifelse(
  is.na(ann_no), fwd_lab(variable), paste0(FAM_ABBR[ann_fam], "-", ann_no))
al <- al %>% left_join(lab_map, by = "variable") %>%
  mutate(ann_fam = ifelse(is.na(ann_fam), "YIELD_CURVE", ann_fam),
         label   = mk_label(variable, ann_fam, ann_no))

# union of per-column top-K variables, per target; rows grouped by family then share
union_rows <- function(tg, k = K) {
  al %>% filter(target == tg) %>% group_by(model_config) %>%
    arrange(desc(share)) %>% slice_head(n = k) %>% ungroup() %>%
    distinct(variable) %>% pull(variable)
}
panel_df <- function(tg, k = K) {
  rows <- union_rows(tg, k)
  d <- al %>% filter(target == tg, variable %in% rows) %>%
    group_by(model_config) %>% mutate(col_norm = 100 * share / max(share)) %>% ungroup() %>%
    mutate(proc = factor(model_config, levels = intersect(PROC_ORDER, unique(model_config))))
  fam_ord <- d %>% group_by(ann_fam) %>% summarise(ft = sum(share), .groups="drop") %>%
    arrange(ft) %>% pull(ann_fam)
  ord <- d %>% group_by(label, ann_fam) %>% summarise(t = sum(share), .groups = "drop") %>%
    mutate(ann_fam = factor(ann_fam, levels = fam_ord)) %>%
    arrange(ann_fam, t) %>% pull(label)
  d %>% mutate(label = factor(label, levels = ord))
}

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold", size = 8.5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7.5),
        plot.title  = element_text(face = "bold", size = 11),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.position = "right")

# ---- Variant A : true-share colour (stepped, identical to Figure 4) -------------
tile_true <- function(d, ttl) {
  ggplot(d, aes(proc, label, fill = share)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(share >= 0.5, sprintf("%.0f", share), ""),
                  colour = share >= 25), size = 2.4) +
    scale_fill_stepsn(colours = pal, breaks = c(2,5,10,20,40), limits = c(0,100),
                      name = "share %", guide = guide_coloursteps(barheight = 4, barwidth = 0.6)) +
    scale_colour_manual(values = c(`TRUE`="white", `FALSE`="grey25"), guide = "none") +
    labs(x = NULL, y = NULL, title = ttl) + base_theme
}
# ---- Variant B : per-column-normalised colour, TRUE share as label --------------
tile_norm <- function(d, ttl) {
  ggplot(d, aes(proc, label, fill = col_norm)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(share >= 0.5, sprintf("%.0f", share), ""),
                  colour = col_norm >= 55), size = 2.4) +
    scale_fill_gradientn(colours = pal, limits = c(0,100),
                         name = "rel. to\ncol. max", guide = guide_colourbar(barheight = 4, barwidth = 0.6)) +
    scale_colour_manual(values = c(`TRUE`="white", `FALSE`="grey25"), guide = "none") +
    labs(x = NULL, y = NULL, title = ttl) + base_theme
}

eq <- panel_df("equity"); bd <- panel_df("bond")

# Variant A
pA <- tile_true(eq, "Equity") / tile_true(bd, "Bond") +
  plot_annotation(title = "Figure 5 — variant A: true effective-importance share (%)",
                  subtitle = "Honest scale; C-RIDGE/C-PLS columns go pale because no variable exceeds ~3%.",
                  theme = theme(plot.title = element_text(face="bold", size=12)))
ggsave("working_paper/figures/scratch_fig5_A_trueshare.png", pA, width = 8.2, height = 11, dpi = 200)

# Variant B
pB <- tile_norm(eq, "Equity") / tile_norm(bd, "Bond") +
  plot_annotation(title = "Figure 5 — variant B: per-column normalised colour (number = true share %)",
                  subtitle = "Each procedure's leader saturates -> no blank columns; magnitude still readable in the numbers.",
                  theme = theme(plot.title = element_text(face="bold", size=12)))
ggsave("working_paper/figures/scratch_fig5_B_colnorm.png", pB, width = 8.2, height = 11, dpi = 200)

# ---- LOCKED candidate: polished standalone Variant B (the figure to integrate) --
tile_final <- function(d, ttl) {
  ggplot(d, aes(proc, label, fill = col_norm)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(share >= 0.5, sprintf("%.0f", share), ""),
                  colour = col_norm >= 55), size = 2.5) +
    scale_fill_gradientn(colours = pal, limits = c(0,100),
                         name = "colour:\nrank within\nprocedure",
                         guide = guide_colourbar(barheight = 4.2, barwidth = 0.6)) +
    scale_colour_manual(values = c(`TRUE`="white", `FALSE`="grey25"), guide = "none") +
    labs(x = NULL, y = NULL, title = ttl) + base_theme
}
pFinal <- tile_final(eq, "Equity") / tile_final(bd, "Bond") +
  plot_annotation(
    title = "Figure 5. Principal variables behind each procedure",
    subtitle = paste0("Cell number = true effective-importance share (%) of that variable in the procedure's fit; ",
                      "colour = its rank within the procedure (each column's leader saturates).\n",
                      "Rows = union of each procedure's top-6 variables, grouped by economic family. ",
                      "Labels are <family>-<Annex A.1 row no.> (e.g. POLU-112); the five C-CP forward rates are CURV-f*."),
    theme = theme(plot.title = element_text(face="bold", size=12),
                  plot.subtitle = element_text(size=8)))
ggsave("working_paper/figures/scratch_fig5_final.png", pFinal, width = 8.4, height = 11, dpi = 200)

# Variant C : concentrated procedures only (drop the diffuse shrinkage pair)
conc <- c("C-CP","C-OLS","C-COMB","C-CT","C-ENET")
ecc <- panel_df("equity") %>% filter(model_config %in% conc) %>% droplevels()
bcc <- panel_df("bond")   %>% filter(model_config %in% conc) %>% droplevels()
# recompute union restricted to concentrated procedures so rows are not shrinkage-driven
union_rows_c <- function(tg) al %>% filter(target==tg, model_config %in% conc) %>%
  group_by(model_config) %>% arrange(desc(share)) %>% slice_head(n=K) %>% ungroup() %>%
  distinct(variable) %>% pull(variable)
mk_c <- function(tg) {
  rows <- union_rows_c(tg)
  d <- al %>% filter(target==tg, model_config %in% conc, variable %in% rows) %>%
    mutate(proc = factor(model_config, levels = intersect(PROC_ORDER, conc)))
  fam_ord <- d %>% group_by(ann_fam) %>% summarise(ft=sum(share),.groups="drop") %>% arrange(ft) %>% pull(ann_fam)
  ord <- d %>% group_by(label, ann_fam) %>% summarise(t=sum(share),.groups="drop") %>%
    mutate(ann_fam=factor(ann_fam,levels=fam_ord)) %>% arrange(ann_fam,t) %>% pull(label)
  d %>% mutate(label = factor(label, levels = ord))
}
pC <- tile_true(mk_c("equity"), "Equity") / tile_true(mk_c("bond"), "Bond") +
  plot_annotation(title = "Figure 5 — variant C (fallback): concentrated procedures only",
                  subtitle = "C-RIDGE/C-PLS dropped (diffuse -> see Figure 4). Cleaner, but loses the shrinkage rows.",
                  theme = theme(plot.title = element_text(face="bold", size=12)))
ggsave("working_paper/figures/scratch_fig5_C_concentrated.png", pC, width = 6.6, height = 9.5, dpi = 200)

cat("wrote scratch_fig5_{A_trueshare,B_colnorm,C_concentrated}.png\n")
cat(sprintf("union rows: equity=%d bond=%d (K=%d)\n", nlevels(eq$label), nlevels(bd$label), K))
cat("\n--- label legend (family abbreviations) ---\n")
print(FAM_ABBR)
cat("\n--- sample label map (top union vars) ---\n")
print(al %>% filter(variable %in% union_rows("equity")) %>% distinct(variable, label, ann_no) %>% arrange(ann_no) %>% as.data.frame())

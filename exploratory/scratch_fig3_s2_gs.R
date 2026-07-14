#!/usr/bin/env Rscript
# scratch_fig3_s2_gs.R — PROTOTYPE: S2-only attribution with a GS selection column.
# SCRATCH. Rebuilds Figure 3 restricted to the S2 (theory-agnostic) pool, and adds
# a "GS" column = the Gram-Schmidt admission profile (selection stage, OLS-class only),
# sitting beside the OLS-class fits it feeds, separated from the shrinkage block that
# bypasses GS. C-CP (S1 forward rates) is excluded (theory pool, different tier).
# Sources: attribution_coefs (per-cell |std coef|, scenario-keyed) + gs_selection_freq.

suppressMessages({library(arrow); library(dplyr); library(ggplot2)})
setwd(".")  # run from the parent-project root (the PIT database)
co <- read_parquet("data/audit/attribution_coefs_2026-04-30.parquet")
g  <- read_parquet("data/audit/gs_selection_freq_2026-04-30.parquet")

fam_of <- function(variable, family) dplyr::case_when(
  !is.na(family) ~ family,
  grepl("^FWD_", variable) | variable %in% c("TERM_SPREAD","DASH_CURVE") ~ "YIELD_CURVE",
  variable %in% c("DEFAULT_SPREAD","X_HYIG_MOM","DASH_DEFAULT") ~ "CREDIT_SPREADS",
  variable == "DASH_PMI_PROXY" ~ "GROWTH", variable == "DASH_REAL2Y" ~ "MONETARY_POLICY_RATES",
  grepl("^X_EURJPY", variable) ~ "FX_CURRENCY", grepl("^X_OILGOLD", variable) ~ "COMMODITIES",
  variable %in% c("X_STKBOND_VAL","X_STKBOND_MOM","X_EMDM_MOM","X_INTLUS_MOM") ~ "CROSS_ASSET",
  TRUE ~ "OTHER")
FAM_LABEL <- c(MONETARY_POLICY_RATES="Monetary policy / rates", YIELD_CURVE="Yield curve",
  CREDIT_SPREADS="Credit spreads", INFLATION="Inflation", GROWTH="Growth / activity",
  LABOR_MARKET="Labor market", FORECASTS="SPF forecasts", EQUITY_VOLATILITY="Equity volatility",
  EQUITY_SENTIMENT="Equity sentiment", EQUITY_VALUATION="Equity valuation",
  FINANCIAL_CONDITIONS_COMPOSITE="Financial conditions", POLICY_UNCERTAINTY_GEOPOLITICAL="Policy uncertainty",
  FX_CURRENCY="FX", COMMODITIES="Commodities", HOUSING="Housing",
  CONSUMPTION_WEALTH="Consumption / wealth", CROSS_ASSET="Cross-asset (pract.)")
OLS3 <- c("C-OLS","C-COMB","C-CT")
ck <- c("target","scenario","model_config","maturity","iteration","horizon","window_kind","window_length","sample")

# ---- fit columns: S2 only, effective importance by family (combine logic, S2-filtered) ----
fitcol <- co %>% filter(grepl("^S2", scenario)) %>%
  mutate(cell = do.call(paste, c(across(all_of(ck)), sep="|")), fam = fam_of(variable, family))
ntot <- fitcol %>% distinct(target, model_config, cell) %>% count(target, model_config, name="n_total")
fitfam <- fitcol %>% group_by(target, model_config, variable, fam) %>%
  summarise(sum_coef = sum(abs(coef)), mean_present = mean(abs(coef)), .groups="drop") %>%
  left_join(ntot, by=c("target","model_config")) %>%
  mutate(eff = if_else(model_config %in% OLS3, sum_coef/n_total, mean_present)) %>%
  group_by(target, model_config, fam) %>% summarise(eff = sum(eff), .groups="drop") %>%
  group_by(target, model_config) %>% mutate(share = 100*eff/sum(eff)) %>% ungroup() %>%
  select(target, col = model_config, fam, share)

# ---- GS column: S2 admission frequency by family, normalized to a share -------------
gscol <- g %>% filter(grain=="by_scenario", scenario %in% c("S2_bond","S2_equity"),
                      model_config=="C-OLS", freq>0) %>%
  mutate(target = ifelse(scenario=="S2_bond","bond","equity"), fam = fam_of(variable, family)) %>%
  group_by(target, fam) %>% summarise(f = sum(freq), .groups="drop") %>%
  group_by(target) %>% mutate(share = 100*f/sum(f)) %>% ungroup() %>%
  mutate(col = "GS") %>% select(target, col, fam, share)

dat <- bind_rows(gscol, fitfam) %>%
  mutate(Target = recode(target, equity="Equity", bond="Bond"),
         fam_lab = recode(fam, !!!FAM_LABEL))
COL_ORDER <- c("GS","C-OLS","C-COMB","C-CT","C-ENET","C-RIDGE","C-PLS")
dat <- dat %>% mutate(col = factor(col, levels = COL_ORDER))
fam_ord <- dat %>% group_by(fam_lab) %>% summarise(t=sum(share),.groups="drop") %>% arrange(t) %>% pull(fam_lab)
dat <- dat %>% mutate(fam_lab = factor(fam_lab, levels = fam_ord))

pal <- c("#eef3f7","#cfe0ec","#9dc0db","#5b91c0","#2c6fbb","#1b4a86")
mk <- function(tg) {
  d <- dat %>% filter(Target==tg) %>% droplevels()
  ggplot(d, aes(col, fam_lab, fill = share)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label = ifelse(share>=0.5, sprintf("%.1f", share), ""), colour = share>=25), size=2.7) +
    scale_fill_stepsn(colours=pal, breaks=c(2,5,10,20,40), limits=c(0,100), name="share %",
                      guide=guide_coloursteps(barheight=4, barwidth=0.6)) +
    scale_colour_manual(values=c(`TRUE`="white",`FALSE`="grey25"), guide="none") +
    labs(x=NULL, y=NULL, title=tg) +
    theme_minimal(base_size=10) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(face="bold", size=8, angle=45, hjust=1),
          axis.text.y=element_text(size=8.5), plot.title=element_text(face="bold", size=11))
}
library(patchwork)
fig <- mk("Bond") / mk("Equity") +
  plot_annotation(title="Figure 3 (S2-only) with a GS selection column",
                  subtitle="GS = Gram-Schmidt admission share by family (selection stage, feeds the OLS-class only). C-OLS/C-COMB/C-CT = freq x |coef|; shrinkage = native weight. C-CP (S1) excluded.")
ggsave("working_paper/figures/scratch_fig3_s2_gs.png", fig, width=8.2, height=10, dpi=200)
cat("wrote scratch_fig3_s2_gs.png\n")
cat("\nGS vs C-OLS top families (the freq != weight check):\n")
print(dat %>% filter(col %in% c("GS","C-OLS")) %>% group_by(Target,col) %>% arrange(desc(share)) %>%
        slice_head(n=3) %>% ungroup() %>% mutate(share=round(share,1)) %>% as.data.frame())

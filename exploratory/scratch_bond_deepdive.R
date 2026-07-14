#!/usr/bin/env Rscript
# scratch_bond_deepdive.R — SCRATCH content for the Bond 1Y/2Y FDR-only deep-dive
# (Part II.4). FDR + committed artifacts only; NO bond SPA (RETAIN_BOND_REGION deferred).
# Sources: stage4_fdr_2026-06-13 (BY-q cascade), p10 is the same universe, and
# attribution_coefs_2026-04-30 (maturity-keyed variable importance).
# Produces four tables; prints them for sign-off before any .Rmd integration.

suppressMessages({library(arrow); library(dplyr)})
setwd(".")  # run from the parent-project root (the PIT database)
fdr <- read_parquet("data/audit/stage4_fdr_2026-06-13.parquet")
co  <- read_parquet("data/audit/attribution_coefs_2026-04-30.parquet")
MAT <- c("bond_1y","bond_2y","bond_3y","bond_5y","bond_10y")
matlab <- c(bond_1y="1Y", bond_2y="2Y", bond_3y="3Y", bond_5y="5Y", bond_10y="10Y")

cat("================ TABLE A — the multiplicity-tax cascade ================\n")
cat("Min dependence-robust BY-FDR q at each nesting level. None clears 0.05.\n\n")
tabA <- fdr %>% filter(fam_self %in% MAT) %>% group_by(fam_self) %>%
  summarise(`per-maturity`=round(min(q_by_self),3), `pooled-bond`=round(min(q_by_pool),3),
            `whole-universe`=round(min(q_by_univ),3), .groups="drop") %>%
  mutate(Maturity=matlab[fam_self]) %>% arrange(match(fam_self,MAT)) %>%
  select(Maturity, `per-maturity`, `pooled-bond`, `whole-universe`)
print(as.data.frame(tabA), row.names=FALSE)

cat("\n================ TABLE B — the borderline cells (top by raw CW p) ================\n")
cat("bond_1y / bond_2y cells with the smallest raw CW p; q at all three nesting levels.\n\n")
tabB <- fdr %>% filter(fam_self %in% c("bond_1y","bond_2y")) %>%
  arrange(p_cw) %>% group_by(fam_self) %>% slice_head(n=3) %>% ungroup() %>%
  transmute(Maturity=matlab[fam_self], Pool=sub("_"," ",scenario), Proc=model_config,
            h=horizon, Window=ifelse(window_kind=="rolling",paste0("roll-",window_length),"expand"),
            Sample=sub("_"," ",sample), n=n_oos, r2=round(r2_oos,3), CW=round(cw,2),
            Sharpe=round(sharpe,2), `raw p`=signif(p_cw,2),
            `q self`=round(q_by_self,3), `q pool`=round(q_by_pool,3), `q univ`=round(q_by_univ,3)) %>%
  distinct()
print(as.data.frame(tabB), row.names=FALSE)

cat("\n================ TABLE C — live-region map (raw signal, pre-multiplicity) ================\n")
cat("Cells with raw CW p<0.01 AND r2_oos>0 (where the apparent signal lives).\n\n")
live <- fdr %>% filter(p_cw<0.01, r2_oos>0)
cat("bond_1y/2y live cells:", sum(live$fam_self %in% c("bond_1y","bond_2y")),
    " | all-bond:", sum(grepl("bond",live$fam_self)), " | equity:", sum(live$fam_self=="equity"), "\n\n")
tabC <- live %>% filter(fam_self %in% c("bond_1y","bond_2y")) %>%
  mutate(Maturity=matlab[fam_self], Window=ifelse(window_kind=="rolling","rolling","expanding")) %>%
  count(Maturity, Window, Sample=sub("_"," ",sample)) %>%
  arrange(Maturity, desc(n))
print(as.data.frame(tabC), row.names=FALSE)
cat("\nlive cells by procedure (bond_1y/2y):\n")
print(live %>% filter(fam_self %in% c("bond_1y","bond_2y")) %>% count(model_config) %>% arrange(desc(n)) %>% as.data.frame(), row.names=FALSE)
cat("\nlive cells by horizon (bond_1y/2y):\n")
print(live %>% filter(fam_self %in% c("bond_1y","bond_2y")) %>% count(horizon) %>% as.data.frame(), row.names=FALSE)

cat("\n================ TABLE D — what carries the borderline (GS-selected procedures) ================\n")
cat("C-OLS/C-COMB effective importance (mean|coef| x admission freq) at bond_1y / bond_2y.\n\n")
for (m in c("bond_1y","bond_2y")) {
  sub <- co %>% filter(maturity==m, model_config %in% c("C-OLS","C-COMB"))
  ncell <- sub %>% distinct(scenario,model_config,iteration,horizon,window_kind,window_length,sample) %>% nrow()
  d <- sub %>% mutate(acoef=abs(coef)) %>% group_by(Variable=variable, Family=family) %>%
    summarise(`admit freq`=round(n()/ncell,2), `mean|coef|`=signif(mean(acoef),3),
              `eff. imp.`=signif(mean(acoef)*n()/ncell,3), .groups="drop") %>%
    arrange(desc(`eff. imp.`)) %>% head(6)
  cat("---", matlab[m], "---\n"); print(as.data.frame(d), row.names=FALSE); cat("\n")
}

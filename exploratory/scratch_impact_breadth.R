#!/usr/bin/env Rscript
# scratch_impact_breadth.R — PROTOTYPE (publication-grade) of the proposed Figure 3:
# two-axis attribution for the GS+OLS data-driven path (S2). Separates the selection
# axis from the coefficient axis (Research-repo / VIP method); does not multiply them.
#   x = Breadth = GS admission frequency (% of S2 trials the selector picks it)
#   y = Impact  = mean |standardized coef| when selected
# A genuine single predictor would sit top-right (often selected AND heavily weighted).

suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork); library(ggrepel)})
setwd(".")  # run from the parent-project root (the PIT database)
co <- read_parquet("data/audit/attribution_coefs_2026-04-30.parquet")
ck <- c("target","scenario","model_config","maturity","iteration","horizon","window_kind","window_length","sample")
BLUE <- "#2c6fbb"; RED <- "#b2182b"

ib <- function(tg) {
  s <- co %>% filter(scenario == paste0("S2_", tg), model_config == "C-OLS") %>%
    mutate(cell = do.call(paste, c(across(all_of(ck)), sep="|")))
  nt <- n_distinct(s$cell)
  s %>% group_by(variable) %>%
    summarise(breadth = 100*n_distinct(cell)/nt, impact = mean(abs(coef)), .groups="drop") %>%
    mutate(target = tg)
}

panel <- function(tg) {
  d <- ib(tg)
  xm <- max(d$breadth)*1.15; ym <- max(d$impact)*1.15
  ythr <- as.numeric(quantile(d$impact, 2/3))   # "heavily weighted" = top third by impact
  ggplot(d, aes(breadth, impact)) +
    # "genuine predictor" corner: often selected (>=50%) AND above-median impact
    annotate("rect", xmin=50, xmax=xm, ymin=ythr, ymax=ym, fill=RED, alpha=0.06) +
    annotate("text", x=(50+xm)/2, y=(ythr+ym)/2, label="genuine predictor\n(empty)",
             colour=RED, alpha=0.8, size=2.6, fontface="italic", lineheight=0.9) +
    geom_hline(yintercept=ythr, colour="grey88", linewidth=0.3) +
    geom_vline(xintercept=50, colour="grey88", linewidth=0.3) +
    geom_point(colour=BLUE, size=2, alpha=0.85) +
    ggrepel::geom_text_repel(aes(label=variable), size=2.4, colour="grey25",
                             min.segment.length=0, segment.size=0.2, max.overlaps=Inf,
                             box.padding=0.4, seed=1) +
    scale_x_continuous(limits=c(0, xm), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, ym), expand=c(0,0)) +
    labs(x="Breadth — GS admission frequency (% of S2 trials)",
         y="Impact — mean abs. std. coefficient when selected",
         title=recode(tg, bond="Bond", equity="Equity")) +
    theme_minimal(base_size=9) +
    theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=11),
          axis.title=element_text(size=8))
}

fig <- panel("bond") + panel("equity") +
  plot_annotation(
    title="Figure 3. No data-driven predictor is both reliably selected and heavily weighted",
    subtitle="GS+OLS path, Scenario 2. Breadth (x) = how often the Gram-Schmidt screen admits the variable; Impact (y) = its mean absolute\nstandardized coefficient when admitted. A genuine predictor would sit top-right (selected in a majority of trials AND in the top third by\nimpact); for both targets that corner is empty.",
    theme=theme(plot.title=element_text(face="bold", size=11), plot.subtitle=element_text(size=8, colour="grey35")))
ggsave("working_paper/figures/scratch_impact_breadth.png", fig, width=9, height=4.8, dpi=200)
cat("wrote scratch_impact_breadth.png\n")

# Family a-priori(ITS) -> PIT-OOS marginal-R2 dumbbell (real data).
# Annex exhibit. Design points:
#   (1) one colour for both measurements: ITS (in-sample, a-priori) = filled dot,
#       OOS (PIT out-of-sample) = same-colour border, white inside. The light-red
#       x<0 band + dashed zero line carry the "no-go / no skill" reading.
#   (2) the decay arrow almost touches the OOS dot, and is drawn only where the
#       ITS->OOS gap exceeds MIN_ARROW (short gaps get no arrow).
#   (3) overall title/subtitle live in the manuscript caption, not the graphic.
suppressMessages({library(arrow); library(dplyr); library(ggplot2); library(patchwork)})
AS_OF <- "2026-04-30"
GAP       <- 0.12  # marginal-R2 (%) the arrowhead leaves before the OOS dot (almost touching)
MIN_ARROW <- 0.3   # ITS->OOS gaps shorter than this get no arrow (below this the dots visually overlap; 0.6 arbitrarily split near-identical gaps, e.g. equity labour -0.603 vs yield curve -0.600)
DOT_COL   <- "#1b7837"  # one colour (dark green): ITS filled, OOS same-colour border + white inside
raw <- read_parquet(sprintf("data/audit/family_signal_its_oos_%s.parquet", AS_OF))

render <- function(its_col, its_q25, its_q75, out, sub) {
  pf <- raw |>
    mutate(its_use = .data[[its_col]] * 100, its_lo = .data[[its_q25]] * 100, its_hi = .data[[its_q75]] * 100,
           oos_med = oos_med * 100, oos_q25 = oos_q25 * 100, oos_q75 = oos_q75 * 100)
  panel <- function(tg, title, show_legend = TRUE) {
    d <- pf |> filter(target == tg) |> arrange(its_use) |>
      mutate(family_label = factor(family_label, levels = family_label),
             arrow_len = oos_med - its_use,
             xend_arrow = oos_med - sign(arrow_len) * GAP)
    pts <- bind_rows(
      transmute(d, family_label, x = its_use, leg = "In-sample R² (a-priori)"),
      transmute(d, family_label, x = oos_med, leg = "PIT out-of-sample R²")) |>
      mutate(leg = factor(leg, c("In-sample R² (a-priori)", "PIT out-of-sample R²")))
    ggplot(d) +
      annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#f7c9c9", alpha = 0.5) +
      annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#e3f3e3", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
      geom_segment(aes(y=family_label, yend=family_label, x=its_lo, xend=its_hi), colour=DOT_COL, alpha=.16, linewidth=1.4) +
      geom_segment(aes(y=family_label, yend=family_label, x=oos_q25, xend=oos_q75), colour=DOT_COL, alpha=.16, linewidth=1.4) +
      geom_segment(data = filter(d, abs(arrow_len) > MIN_ARROW),
                   aes(y=family_label, yend=family_label, x=its_use, xend=xend_arrow),
                   colour="grey45", arrow=arrow(length=unit(0.06,"inches"), type="closed"), linewidth=0.45) +
      geom_point(data = pts, aes(x, family_label, fill = leg), shape = 21, colour = DOT_COL, size = 2.6, stroke = 0.7) +
      scale_fill_manual(values = c("In-sample R² (a-priori)" = DOT_COL,
                                   "PIT out-of-sample R²" = "white"), name = NULL) +
      labs(x = "Marginal R² (%)", y = NULL, title = title) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = if (show_legend) "top" else "none", legend.justification = "center")
  }
  g <- panel("bond","Bond (1Y, h=1)", TRUE) / panel("equity","Equity (h=1)", FALSE)
  ggsave(out, g, width = 8.5, height = 9.5, dpi = 150)
}
render("its_med", "its_q25", "its_q75",
       "working_paper/figures/fig_family_dumbbell.png",
       "ITS = full-sample in-sample marginal R² (1990–2026); OOS = PIT walk-forward. The gap is overfit + regime decay.")
cat("done\n")

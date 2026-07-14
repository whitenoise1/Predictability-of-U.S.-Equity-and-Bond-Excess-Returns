# working_paper/merge_bond_forecasts.R
# Bind the per-(config,sample,maturity) bond forecast slices, enforce the G1 guard
# (recomputed r2/cw/sharpe must reproduce ALL committed drop-it1 ch2 bond
# shrinkage cells < 1e-8, every maturity), and write the consolidated artifact.
# STOPS on any mismatch (no silent inconsistency). RUN FROM REPO ROOT.
suppressMessages({library(arrow); library(dplyr); library(purrr)})
source("R/predictive_metrics.R"); source("R/predictive_regression.R")
G1_REF <- "data/audit/p10_rescored_drop-it1_2026-06-17.parquet"
SHRINK <- c("C-ENET","C-RIDGE","C-PLS")
KEY <- c("chapter","scenario","target","maturity","model_config","iteration",
         "horizon","window_kind","window_length","sample")

slices <- list.files("data/audit", pattern = "^_bondfc_slice_.*\\.parquet$", full.names = TRUE)
if (!length(slices)) stop("no slice files found — run rerun_bond_forecasts.R slices first")
fc <- bind_rows(lapply(slices, read_parquet))
cat(sprintf("[merge] %d slices -> %d forecast rows / %d cells\n", length(slices), nrow(fc),
            nrow(distinct(fc, across(all_of(KEY))))))

# recompute scalars per cell
inf <- fc |> arrange(across(all_of(KEY)), as_of_date) |>
  group_by(across(all_of(KEY))) |> group_split() |> map_dfr(function(d) {
    h <- d$horizon[1L]
    d[1L, KEY] |> bind_cols(tibble(
      r2_oos = r2_oos(d$y_hat, d$y_realized, d$y_bench),
      cw = clark_west_stat(clark_west_pointwise(d$y_hat, d$y_realized, d$y_bench), L = floor(1.5*h))$stat,
      sharpe = signal_sharpe(d$y_hat, d$y_realized, h = h)))
  })

# G1 vs committed drop-it1 ch2 bond_1y shrinkage cells
ref <- read_parquet(G1_REF) |>
  filter(chapter == "ch2", target == "bond", scenario == "S2_bond",
         model_config %in% SHRINK, !is.na(accept_p10))
chk <- inf |> select(all_of(KEY), r2_oos, cw, sharpe) |>
  inner_join(ref |> select(all_of(KEY), r2_s = r2_oos, cw_s = cw, sh_s = sharpe), by = KEY)
cat(sprintf("[G1] recomputed cells=%d ; committed ref=%d ; matched=%d\n", nrow(inf), nrow(ref), nrow(chk)))
if (nrow(chk) != nrow(ref)) {
  miss <- anti_join(ref |> select(all_of(KEY)), inf |> select(all_of(KEY)), by = KEY)
  cat("UNMATCHED committed cells:\n"); print(as.data.frame(miss))
  stop(sprintf("[G1 FAIL] %d committed bond_1y cells unmatched by the re-run", nrow(ref) - nrow(chk)))
}
g1 <- chk |> summarise(dr2 = max(abs(r2_oos - r2_s)), dcw = max(abs(cw - cw_s)), dsh = max(abs(sharpe - sh_s)))
cat(sprintf("[G1] max|dr2|=%.2e max|dcw|=%.2e max|dsharpe|=%.2e\n", g1$dr2, g1$dcw, g1$dsh))
if (!(g1$dr2 < 1e-8 && g1$dcw < 1e-8 && g1$dsh < 1e-8))
  stop("[G1 FAIL] re-run does not reproduce committed bond scalars < 1e-8 — STOP")
cat("[G1 OK] all bond shrinkage scalars reproduced < 1e-8\n")

out <- file.path("data/audit", sprintf("bond_s2_forecasts_%s.parquet", Sys.Date()))
write_parquet(fc, out)
cat(sprintf("[write] %s (%d rows)\n", out, nrow(fc)))
invisible(file.remove(slices)); cat(sprintf("[cleanup] removed %d slice files\n", length(slices)))
cat("[merge_bond_forecasts complete]\n")

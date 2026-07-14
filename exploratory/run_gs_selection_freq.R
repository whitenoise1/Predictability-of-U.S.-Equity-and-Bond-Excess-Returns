# run_gs_selection_freq.R — S-B: GS / admission frequency layer for the uniform
# attribution build. Under the committed-iters decision (2026-06-15, D1) there is
# no walk-forward GS to run: a variable's admission frequency is simply the
# fraction of an OLS-class group's cells whose final-as-of fit admitted it, read
# straight off the S-A snapshot artifact (attribution_coefs). This is the freq term
# in the OLS-class effective weight (freq × |std coef|) that S-C combines, and the
# "raw GS selection-frequency" panel the figure shows alongside.
#
# Scope = OLS-class only (C-OLS/C-COMB/C-CT/C-CP). Shrinkage (C-RIDGE/C-PLS/C-ENET)
# is weighted by |coef| directly, not by admission frequency, so it is excluded
# here. S1 fixed pools + C-CP admit every pool member (freq=1) — recorded
# explicitly so the panel is honest, not implied. Pure aggregation, no fitting.

suppressMessages({library(arrow); library(dplyr); library(tidyr)})
AS_OF     <- "2026-04-30"
OLS_CLASS <- c("C-OLS","C-COMB","C-CT","C-CP")
S1_FIXED  <- c("S1_bond_FF1989","S1_bond_CP","S1_equity_WGCRR")   # always-admit scenarios
IN  <- sprintf("data/audit/attribution_coefs_%s.parquet", AS_OF)
OUT <- sprintf("data/audit/gs_selection_freq_%s.parquet", AS_OF)

cell_keys <- c("target","scenario","model_config","maturity","iteration","horizon",
               "window_kind","window_length","sample")
a <- read_parquet(IN) |> filter(model_config %in% OLS_CLASS) |>
  mutate(cell_id = do.call(paste, c(across(all_of(cell_keys)), sep = "|")))

# Admission frequency at a given grouping grain: distinct cells where a variable was
# admitted / distinct cells in the group.
freq_at <- function(gk) {
  totals  <- a |> distinct(across(all_of(c(gk, "cell_id")))) |>
    count(across(all_of(gk)), name = "n_cells_total")
  present <- a |> distinct(across(all_of(c(gk, "variable", "family", "cell_id")))) |>
    count(across(all_of(c(gk, "variable", "family"))), name = "n_cells_present")
  present |> left_join(totals, by = gk) |>
    mutate(freq = n_cells_present / n_cells_total)
}
pooled <- freq_at(c("target","model_config")) |> mutate(scenario = NA_character_, grain = "pooled")
by_scn <- freq_at(c("target","model_config","scenario")) |> mutate(grain = "by_scenario")
gs_freq <- bind_rows(pooled, by_scn) |>
  select(grain, target, model_config, scenario, variable, family,
         n_cells_present, n_cells_total, freq) |>
  arrange(grain, target, model_config, scenario, desc(freq))

# ---- verification gates (no-silent-fails) ----------------------------------------
stopifnot(sum(is.na(gs_freq$freq)) == 0L,
          all(gs_freq$freq > 0 & gs_freq$freq <= 1 + 1e-9))
# Every variable in an S1 fixed pool / C-CP admits in 100% of its scenario's cells.
s1_bad <- by_scn |> filter(scenario %in% S1_FIXED, abs(freq - 1) > 1e-9)
if (nrow(s1_bad)) { print(as.data.frame(s1_bad)); stop("S-B: S1/CP variable with freq != 1 (fixed pool should always admit).") }
# n_cells_total per (target,config,scenario) must reconcile with the S-A cell counts.
sa_cells <- a |> distinct(across(all_of(c("target","model_config","scenario","cell_id")))) |>
  count(target, model_config, scenario, name = "sa_cells")
chk <- by_scn |> distinct(target, model_config, scenario, n_cells_total) |>
  left_join(sa_cells, by = c("target","model_config","scenario"))
stopifnot(all(chk$n_cells_total == chk$sa_cells))

write_parquet(gs_freq, OUT)
cat(sprintf("[S-B] wrote %s (%d rows: %d pooled + %d by-scenario)\n", OUT, nrow(gs_freq),
            nrow(pooled), nrow(by_scn)))
cat("\n=== reconciliation: cells per (target, config, scenario) ===\n")
chk |> arrange(target, model_config, scenario) |> as.data.frame() |> print(row.names = FALSE)
cat("\n=== top-5 admission frequency per pooled (target, config) ===\n")
pooled |> group_by(target, model_config) |> arrange(desc(freq), .by_group = TRUE) |>
  slice_head(n = 5) |> mutate(freq = round(freq, 3)) |>
  select(target, model_config, variable, family, n_cells_present, n_cells_total, freq) |>
  as.data.frame() |> print(row.names = FALSE)

# run_attribution_all_combine.R — S-C1: combine the S-A snapshot coefficients with
# the S-B admission frequencies into one uniform variable-importance metric across
# ALL procedures, normalize to a comparable share scale, write the combined artifact,
# and render prototype figures for review (NO .Rmd touched — that is S-C2).
#
# Effective-importance rule (locked 2026-06-15, decision 1):
#   OLS-class (C-OLS/C-COMB/C-CT):  freq × |std coef|  ==  mean |coef| over ALL the
#       group's cells (a non-selected cell contributes 0) — this folds GS selection in.
#   C-CP:                           |γ|  (the 5 forward rates always admit; freq = 1).
#   Shrinkage (C-RIDGE/C-PLS/C-ENET): native weight directly, NO frequency — ridge
#       |std coef|, PLS net loading sqrt(Σ loading²), ENET |coef| — mean over the cells
#       where the variable was available (availability is not a selection step).
# Each (procedure × target) is then normalized to shares (% within fit) so the four
# different native metrics sit on one comparable scale (decision 1: "one share scale").

suppressMessages({library(arrow); library(dplyr); library(tidyr)})
AS_OF <- "2026-04-30"
A_IN  <- sprintf("data/audit/attribution_coefs_%s.parquet", AS_OF)
F_IN  <- sprintf("data/audit/gs_selection_freq_%s.parquet", AS_OF)
OUT   <- sprintf("data/audit/attribution_all_%s.parquet", AS_OF)
OLS3  <- c("C-OLS","C-COMB","C-CT")
cell_keys <- c("target","scenario","model_config","maturity","iteration","horizon",
               "window_kind","window_length","sample")

# ---- 1. combine into the effective-importance metric (pooled procedure × target) -
a <- read_parquet(A_IN) |> mutate(cell_id = do.call(paste, c(across(all_of(cell_keys)), sep = "|")))
ntot <- a |> distinct(target, model_config, cell_id) |> count(target, model_config, name = "n_total")
gsf  <- read_parquet(F_IN) |> filter(grain == "pooled") |>
  select(target, model_config, variable, gs_freq = freq)

attribution_all <- a |>
  group_by(target, model_config, variable, family) |>
  summarise(sum_coef = sum(coef), mean_present = mean(coef), n_present = n_distinct(cell_id),
            .groups = "drop") |>
  left_join(ntot, by = c("target","model_config")) |>
  left_join(gsf,  by = c("target","model_config","variable")) |>
  mutate(metric_kind = case_when(model_config %in% OLS3 ~ "freq x |std coef|",
                                 model_config == "C-CP"   ~ "|gamma|",
                                 model_config == "C-RIDGE"~ "|std coef|",
                                 model_config == "C-PLS"  ~ "PLS net loading",
                                 model_config == "C-ENET" ~ "|coef|"),
         eff_weight  = if_else(model_config %in% OLS3, sum_coef / n_total, mean_present)) |>
  group_by(target, model_config) |>
  mutate(share = 100 * eff_weight / sum(eff_weight)) |>
  arrange(desc(share), .by_group = TRUE) |> mutate(rank = row_number()) |> ungroup() |>
  select(target, model_config, metric_kind, variable, family, rank, share, eff_weight,
         gs_freq, mean_present_coef = mean_present, n_present, n_total)

# ---- 2. coverage audit (no-silent-fails) -----------------------------------------
cov <- attribution_all |> count(target, model_config, name = "n_vars")
expect <- tibble(
  target = c(rep("bond",6), rep("equity",6)),
  model_config = c("C-OLS","C-COMB","C-CP","C-ENET","C-RIDGE","C-PLS",
                   "C-OLS","C-COMB","C-CT","C-ENET","C-RIDGE","C-PLS"))
miss <- anti_join(expect, cov, by = c("target","model_config"))
if (nrow(miss)) { print(as.data.frame(miss)); stop("S-C1 coverage: a (target,config) cell is unattributed.") }
stopifnot(nrow(cov) == 12L, all(cov$n_vars >= 1L),
          all(abs((attribution_all |> group_by(target,model_config) |>
                   summarise(s = sum(share), .groups="drop"))$s - 100) < 1e-6))
write_parquet(attribution_all, OUT)
cat(sprintf("[S-C1] wrote %s (%d rows; 12 procedure×target cells). Figure: working_paper/fig_attribution.R\n", OUT, nrow(attribution_all)))

cat("\n=== coverage (vars attributed per procedure × target) ===\n")
cov |> arrange(target, model_config) |> as.data.frame() |> print(row.names = FALSE)
cat("\n=== top-3 effective-importance share per procedure × target ===\n")
attribution_all |> filter(rank <= 3) |> mutate(share = round(share, 1)) |>
  select(target, model_config, variable, share, gs_freq) |> as.data.frame() |> print(row.names = FALSE)

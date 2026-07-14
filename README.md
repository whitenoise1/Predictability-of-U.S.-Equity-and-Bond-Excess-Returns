# Predictability of U.S. Equity and Bond Excess Returns

Research document, data and source code of the working paper: "Predictability of U.S. Equity and Bond Excess Returns: a Point-in-Time Story", Stefan Bolta (2026).

**Keywords:** return predictability; point-in-time data; macroeconomic factors; out-of-sample forecasting; multiple-testing / data snooping; equity and bond risk premia.
**JEL Classification:** G12, G14, G17, C53, C55, E44.

## Abstract

Can macroeconomic data forecast monthly U.S. stock and Treasury-bond excess returns? Our attempt at an old problem is a leakage-free, point-in-time (PIT) panel of 506 macroeconomic series, on which we run a wide, pre-specified search across three forecast horizons, five bond maturities (and the S&P500), seven estimation windows and four sub-samples. We score every forecast against the naive historical-average benchmark. The result is a near-null, near-symmetric across targets: the macro data clear our bar in only 0.0% of equity trials and 5.0% of bond trials. The classic finance theory pools accept exactly zero, and the modest data-driven pool that shows potential is dismantled by the dependence-robust false-discovery control that a search this large demands — save a single pre-2008 short-end-bond cluster that clears the snooping-robust SPA but still fails the dependence-robust false-discovery control and does not persist past 2008.

## Repository contents

| File | Description |
|---|---|
| `Predictability-of-U.S.-Equity-and-Bond-Excess-Returns.pdf` | The research document (working-paper build of 14 July 2026). |
| `data/` | Every data artifact the manuscript reads at render time — 17 Parquet files sorted by the paper's narrative: the trial acceptance universe (6,768 trials / 1,806 evaluable), Stage-4 FDR and SPA inference, the retained pointwise forecast paths (equity, all five bond maturities, and the theory pools), panel provenance and PIT audit, predictor-pool manifests, attribution, and the vintage-leakage counterfactual. See `data/README.md` for the full manifest (object names, dimensions, grain, and where each is used in the paper). |
| `pipeline/` | The manuscript-support pipeline: artifact generators (target-relative acceptance re-scoring, the drop-iteration-1 consolidation, Stage-4 multiple-testing inference, full-curve bond forecast retention, the corpus revision audit, the vintage-leakage payoff) and the figure scripts. |
| `exploratory/` | Scratch / one-off exploration scripts kept for provenance — not part of the manuscript build. |

## How to explore / reproduce

1. **Requirements:** R ≥ 4.0. The `data/` artifacts are self-contained Parquet files; no external download is needed to explore them.
2. **Install the required packages:**

```r
install.packages(c("arrow", "dplyr", "tidyr", "ggplot2", "knitr",
                   "kableExtra", "patchwork", "ggrepel"))
```

3. **Load any artifact directly:**

```r
library(arrow)
p10 <- read_parquet("data/03_trials_acceptance/p10_rescored_drop-it1_2026-06-17.parquet")
```

**Note on scope:** the manuscript renders exclusively from the committed artifacts in `data/` — no model is re-estimated at render time. The `pipeline/` scripts are the generators of those artifacts and are distributed for methodology transparency; they are written to run from the root of the parent project (a private point-in-time macro database holding the underlying 506-series raw vintage store), so they will not run standalone from this repository. Every number in the paper traces to a committed artifact via `stopifnot()` pins in the manuscript source, and the heavy computations (the ~14-hour full-curve bond forecast retention, the Stage-4 bootstrap) are distributed precomputed in `data/`.

## Main pipeline scripts

* **`rescore_p10.R`**: re-scores every walk-forward trial under the target-relative acceptance rule (the strategy Sharpe must beat the target's own PIT-expanding median Sharpe path).
* **`run_drop_it1_rescore.R`**: consolidates the data-driven path onto the shrinkage estimator class and produces the acceptance and Stage-4 inference artifacts the paper reports.
* **`rerun_bond_forecasts.R` / `rerun_equity_corrected.R` / `rerun_s1_forecasts.R`**: retain the pointwise forecast paths — all five bond maturities, the corrected month-end equity target, and the classic-literature theory pools.
* **`run_stage4_bond_spa.R`**: Hansen (2005) superior-predictive-ability test across the full bond curve.
* **`run_revision_audit.R`**: first-release versus latest-vintage revision audit across every genuinely-revised series (Annex A.2).
* **`run_leakage_payoff.R`**: the counterfactual quantifying what a revised-data (vintage-leaked) shortcut would buy at fixed release timing (Annex A.2).

## License

This project is distributed under the [GNU GPL v3](LICENSE) license.

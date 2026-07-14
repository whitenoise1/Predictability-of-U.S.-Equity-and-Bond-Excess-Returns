# Manuscript data inventory

Every data object read by `Predictability of U.S. Equity and Bond Excess Returns.Rmd` at render
time, copied here as an exploration/validation snapshot (copied **2026-07-13**) and sorted by the
paper's narrative. **Canonical copies live at the project root** (`data/audit/`,
`data/metadata/`) — the Rmd reads those, *not* these copies. If a root artifact is regenerated,
re-copy it here.

The manuscript loads 17 distinct files (the annex table re-reads
`annex_used_variables_2026-06-14.parquet` a second time as object `av`). All headline scalars in
the paper are pinned by `stopifnot()` guards in the Rmd setup chunk, so a silent drift in any of
these artifacts fails the render.

## 01_panel_provenance — Section 3 (Data), Annex A.2

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `factor_metadata.parquet` | `meta` | 506 × 14 | ingested series | corpus counts (506 = 308 SPF + 198), PIT-tier and frequency tables, series-to-feature ladder |
| `pit_quality_report_2026-06-09.parquet` | `aud` | 1518 × 6 | series × audit test | "the store is audited" scalars: 198/198 monotonicity, 440/440 vintage-present, 102 release-lag flags (93 RTDSM + 9 shutdown) |
| `revision_audit_2026-06-23.parquet` | `rev` | 130 × 10 | genuinely-revised series | revision noise-to-signal exhibit: median 0.29, 25% > 0.5, 53% upward (Annex A.2, §3 footnote) |
| `revision_gdp_example_2026-06-23.parquet` | `rgx` | 1 × 2 | — | 2008-Q4 real GDP first-print −3.8% vs latest −8.5% example |

## 02_pools_features — §2.5, §3 ladder, Annex A.3

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `s2_pool_manifest_2026-06-06.parquet` | `s2man` | 114 × 5 | Scenario-2 pool feature | S2 all-variables pool composition: 83 macro + 18 SPF M0 + 13 dashboard |
| `annex_used_variables_2026-06-14.parquet` | `usedv`, `av` | 119 × 8 | model feature | audited used-feature table (Annex A.3); feeds the 506 → 119 ladder (101 direct + 18 derived) |

## 03_trials_acceptance — Section 4 Part I (headline near-null)

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `p10_rescored_drop-it1_2026-06-17.parquet` | `p10_all` → `p10`, `ch2`, `ch3` | 6768 × 20 | trial (scenario × target × maturity × config × h × window × sample) | THE acceptance artifact: 1,806 evaluable trials (276 equity / 1,530 bond), P10 target-relative rule (`accept_p10`), all acceptance rates, per-maturity bond counts 47/14/11/4/0 |

## 04_forecast_paths — bridge figures, pooled OOS, SPA inputs

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `equity_s2_forecasts_2026-06-12.parquet` | `eqfc` | 68,472 × 12 | as-of date × trial (equity S2) | retained pointwise equity forecasts (it-1 configs filtered out at load); first scored as-of 2004-12-31 |
| `bond_s2_forecasts_2026-06-21.parquet` | `bdfc` | 57,040 × 14 | as-of date × trial (bond S2, all 5 maturities) | full-curve retained bond forecasts (G1-pinned re-run); bridge figure shows the 1Y |
| `s1_forecasts_2026-06-20.parquet` | `s1fc` | 94,946 × 14 | as-of date × trial (theory pools) | Annex theory bridge: WG∪CRR (±surprises), Fama–French 1989, Cochrane–Piazzesi; pool starts 1991-05 (CP), 1996-08 (FF), 2016-07 (WG∪CRR) |

## 05_inference_stage4 — Section 4 Part II (multiple testing)

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `stage4_fdr_drop-it1_2026-06-17.parquet` | `fdr` | 1806 × 32 | evaluable trial | BH/BY FDR q-values at three family levels (self/pool/universe); raw-alive bond breadth (56 trials) |
| `stage4_spa_drop-it1_2026-06-17.parquet` | `spa` | 16 × 15 | family × horizon × sample | Hansen SPA, equity + bond-1y families |
| `stage4_spa_bond_2026-06-21.parquet` | `spabd` | 40 × 16 | maturity × horizon × sample | SPA across the full bond curve; only pre-2008 bond-1y significant (h3 p=.038, h1 p=.042) |
| `stage4_equity_inference_drop-it1_2026-06-17.parquet` | `eqi` | 144 × 25 | equity-S2 trial | overlap-robust Clark–West (t on n/h df), Mincer–Zarnowitz w/ Hodrick SEs, equity FDR |

## 06_attribution_families — Part I.2 diagnostics

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `attribution_all_2026-04-30.parquet` | `attall` | 752 × 12 | config × variable | shrinkage diffuseness scalar (top single-variable share ≤ 3%), attribution figure |
| `family_signal_its_oos_2026-04-30.parquet` | `fam_sig` | 34 × 15 | target × indicator family | family dumbbell: equity 0/17 families positive OOS; bond 8/17 positive, 5 improve OOS |

## 07_leakage_counterfactual — §2.6 footnote, Annex A.2, Table 14

| File | Rmd object | Dims | One row per | Used for |
|---|---|---|---|---|
| `leakage_payoff_combined_2026-06-24.parquet` | `lkp` → `lkw`, `lkb` | 36 × 15 | trial × {pit, leaked} panel | vintage-leak payoff: 18 paired cells, 0 accept flips, 7/9 bond cells inflate (PLS Sharpe 0.82 → 1.19) |

## Static figures (derived, not loaded as data)

One pre-rendered PNG is included by the Rmd, built from an artifact above:

- `figures/fig_family_dumbbell.png` ← `family_signal_its_oos_2026-04-30.parquet` (`pipeline/fig_family_its_oos.R`)

(Figure 5, the Scenario-1 theory bridge, is drawn inline as vector by the Rmd chunk
`theory-bridge` from `s1_forecasts` + `p10_rescored_drop-it1`; the old PNG generator
`pipeline/fig_s1_theory.R` is superseded.)

## Validation status (2026-07-13)

All 17 files load cleanly with `arrow::read_parquet()`; dimensions above are as read. Setup-chunk
`stopifnot` pins reproduced outside the render: trial counts (6,768 total / 1,806 evaluable),
per-maturity bond acceptances, revision-audit scalars, leakage-payoff scalars, PIT-audit counts,
and the OOS-window start dates (S2 2004-12-31; CP 1991-05-31; FF 1996-08-31; WG∪CRR 2016-07-31;
last scored as-of 2026-03-31).

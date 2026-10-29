# Open aortic arch surgery in the Netherlands, 2013–2024: independent statistical analysis

This is a reproducible R pipeline for the nationwide NHR cohort of open aortic arch surgery in 15 Dutch centres. It covers national volume, case mix, early outcomes, failure to rescue, organ-protection practice, centre variation, long-term survival and survival relative to the general population.

The analysis was developed independently. It is based only on the extract's variable list, the NHR Handbook Cardiochirurgie (v22.2.1), the research questions and established methods. No earlier scripts or results were used.

| Document | Content |
|---|---|
| `statistical_analysis_plan.md` | Questions, endpoints, models, covariates, missing data, multiplicity, sensitivity analyses. Written before analysis. |
| `statistical_methods_report.md` | Why each method was chosen, the assumptions and checks, and how to interpret the results. |
| `results/analysis_report.md` | Generated on each run: every key number, read from the fitted models. |

> **Status.** The pipeline has not yet been run on the patient data. Nothing in this repository is a study result.

## Folder structure

```
aortic_arch_analysis/
├── R/00_config.R              paths, NHR codings, ranges, analysis settings (edit here)
├── R/functions.R              shared functions (model fitting, pooling, splines, survival tools)
├── 01_data_validation.R       structure, codings, ranges, duplicates, consistency, derived variables
├── 02_cohort_descriptives.R   flow diagram, Table 1, case-mix trends, missing data, multiple imputation
├── 03_volume_trends.R         annual counts and rates, Poisson/negative binomial, acute vs elective
├── 04_clinical_outcomes.R     risk-adjusted trends, standardised risks, centre variation, LOS, failure to rescue
├── 05_organ_protection.R      temperature, hypothermia, circulatory arrest; within/between-centre trends
├── 06_survival_analysis.R     Kaplan-Meier, Cox, PH checks, piecewise HRs, RMST, landmark, conditional survival
├── 07_population_survival.R   expected survival, Ederer II, Pohar-Perme, SMR, excess hazard
├── 08_sensitivity_analyses.R  prespecified sensitivity analyses
├── 09_tables_figures.R        Tables 1-3, figures, supplement workbook, consistency checks, report
├── run_all.R                  runs 01-09 in fresh R sessions
├── data/external/             CBS 37360ned (life tables) and 03759ned (population): public, CC BY 4.0
├── data/derived/              intermediate files (never committed)
└── results/                   all output (never committed)
```

## How to run

1. **R ≥ 4.3** with the packages `data.table`, `ggplot2`, `lme4`, `survival`, `mice`, `ordinal`, `clubSandwich`, `MASS`, `patchwork`, `scales`, `writexl`, `readxl` and `haven`. On Ubuntu they are available as `r-cran-*` packages; elsewhere, use `install.packages()`. Exact versions are written to `results/session_info.txt`.
2. The data file is read directly from your desktop: **`AORTABOOG_NL_schoon v3(1).xlsx`** (first sheet). The code looks in `Desktop` and `Bureaublad`, including a desktop redirected to OneDrive. The file is only read, never changed or copied. If the file name, sheet or location changes, edit `RAW_FILE_NAME` / `RAW_FILE_SHEET` at the top of `R/00_config.R`, or set the full path in the environment variable `ARCH_RAW_DATA`.
3. From this folder, run `Rscript run_all.R`.
4. **Review `results/01_validation/validation_report.md` first.** By default the pipeline stops if a coded variable contains values that are not in the NHR handbook (`strict_coding = TRUE`). This is deliberate: it is better to stop than to guess what an undocumented code means.

### Settings (environment variables or `R/00_config.R`)

| Variable | Default | Meaning |
|---|---|---|
| `ARCH_M` | 20 | number of imputations |
| `ARCH_NSIM` | 500 | simulation draws for standardised-risk intervals |
| `ARCH_CORES` | cores − 1 | parallel model fits across imputations (Linux/macOS) |
| `ARCH_GLMM_FALLBACK` | `bobyqa` | last-resort refit for GLMM convergence warnings (`none` = faster; warnings are still reported) |

**Runtime.** With m = 20 on a 4-core machine, expect several hours, mostly in `04_clinical_outcomes.R`.

## Main outputs (`results/`)

* `publication/`:
  * Table 1 (baseline characteristics by urgency), Table 2 (organ protection), Table 3 (early outcomes);
  * Figures 1–8: volume, case mix, adjusted trends, organ protection, centre variation, Kaplan–Meier curves, landmark curves, relative survival;
  * supplementary figures S1–S4;
  * `Tables_and_Supplement.xlsx`;
  * `consistency_checks_tables_figures_models.csv`.
* `01_validation/` … `08_sensitivity/`: every model, diagnostic and intermediate table as CSV.
* `analysis_report.md`, `session_info.txt`, `logs/`.

All tables and figures are produced from the fitted analyses. None contain values entered by hand.

## Data protection

* Patient-level data, derived data and results are excluded by `.gitignore`. Keep the extract on secured local storage.
* Centres are anonymised (A–O, ordered by volume) in every output. The key is stored only in `data/derived/centre_key.rds`.
* Outputs are aggregates. Before external sharing, check small cells against NHR publication rules.

## Points to confirm with the data provider

These are flagged in the analysis plan (§1.3) and are not assumed:

* that `bronrij` is a row number and there is no patient identifier;
* that `jaar` is the operation year, `leeftijd` is the age in completed years, `overltijd_dag` is the number of days from operation to death or the status check, and `opnameduur` is the number of days from operation to discharge;
* the rules behind the `*_clean`, `pathologie_groep`, `klasse`, `s1_via_bewijs`, `excl_reden`/`incl_cohort` and `overl120d` variables;
* the eGFR formula;
* whether `circ_arrest_duur` includes periods with antegrade cerebral perfusion;
* whether 2024 registrations are complete.

## Data sources

* NHR Handbook Cardiochirurgie v22.2.1 (codings).
* CBS StatLine 37360ned, *Levensverwachting; geslacht, leeftijd* (version 22 June 2026).
* CBS StatLine 03759ned, *Bevolking op 1 januari en gemiddeld* (exported 18 July 2026).

Both CBS tables are published under CC BY 4.0 (source: CBS).

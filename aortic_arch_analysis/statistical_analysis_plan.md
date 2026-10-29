# Statistical Analysis Plan

**Study:** Nationwide trends, early outcomes and long-term survival after open aortic arch surgery in the Netherlands, 2013–2024 (15 cardiac surgery centres, Netherlands Heart Registration [NHR]).

**Version:** 1.0, written before any analysis of the patient-level data.
**Basis:** the extract's variable list (60 columns), NHR Handbook Cardiochirurgie v22.2.1 (codings and definitions), the research questions, and established methods. No earlier analysis scripts, outputs or reported results were consulted.

**Reporting guidelines:** STROBE (observational cohorts), RECORD (routinely collected data), and the EJCTS/JTCVS statistical and data reporting guidelines (Hickey et al., 2015). Estimates are reported with 95% confidence intervals (CIs); p-values are secondary.

---

## 1. Data source and the variables available

### 1.1 Variables in the extract

| Domain | Variables |
|---|---|
| Identification | `bronrij`, `centrum`, `jaar` |
| Demographics | `leeftijd`, `geslacht`, `BMI` |
| Preoperative risk | `kreatinine`, `eGFR`, `LVEF`, `cardiochir_eerder`, `endocarditis`, `urgentie`, `EuroscoreII` |
| Operative | `ECC`, `circ_arrest`, `circ_arrest_clean`, `chir_benadering`, `CABG`, `klepchir`, `klepchir_clean`, `aklepchir`, `mklepchir`, `pklepchir`, `tklepchir`, `aortachir`, `aortachir_clean`, `aotachir_benadering`, `aorta_asc`, `aorta_asc_path`, `aortaboog_path`, `aorta_desc`, `aorta_desc_path`, `cardiochir_overig`, `pathologie_groep`, `klasse` |
| Perfusion | `min_temp`, `ECC_duur`, `AoX_duur`, `circ_arrest_duur` |
| In-hospital outcomes | `overleden_zh`, `opnameduur`, `ECMO_opname`, `resp_insuff`, `beademing`, `heropname_IC`, `CVA_restletsel`, `CVA_zonder_restletsel`, `nierfalen`, `gastroint_compl`, `vasc_compl_opname`, `delier`, `postop_MI` |
| Follow-up | `mort_status`, `overltijd_dag`, `overl30d`, `overl120d`, `overl1j` |
| Cohort selection | `s1_via_bewijs`, `excl_reden`, `incl_cohort` |

### 1.2 Codings taken from the NHR handbook (v22.2.1)

* Binary variables: 0 = no, 1 = yes, −1 = unknown (−1 is set to missing).
* `geslacht`: 1 = male, 2 = female, −1 = unknown.
* `urgentie`: 10 = elective, 20 = urgent, 30 = emergency, 40 = salvage, −1 = unknown.
* `ECC`: 0 = none, 10 = type unknown, 20 = conventional, 30 = mini-ECC.
* `chir_benadering`: 10 = sternotomy, 20 = other non-minimally invasive, 30 = mini-sternotomy, 40 = mini-thoracotomy, 50 = percutaneous, 60 = other minimally invasive.
* `aotachir_benadering` (named `aortachir_benadering` in the handbook): 0 = none, 10 = open, 20 = simple (T)EVAR, 21 = complex (T)EVAR, 30 = combined open + (T)EVAR. The handbook states that a Thoraflex/frozen elephant trunk arch is coded 30.
* `aorta_asc`: 0 = no, 10 = isolated root, 20 = ascending including root, 30 = ascending excluding root, 90 = location unknown.
* `aorta_desc`: 0 = no, 1 = location unknown, 10 = thoracic, 20 = abdominal, 30 = thoraco-abdominal.
* `*_path`: 0 = none, 10 = aneurysm, 20 = acute dissection (≤ 2 weeks), 30 = chronic dissection, 40 = infection, 91 = other (90 in handbook versions before 22.2), −1 = unknown.
* `gastroint_compl`: 0 = none, 1 = type unknown, 10 = bleeding, 90 = other.
* `vasc_compl_opname`: 0 = none, 10 = major, 20 = minor (VARC-2).
* `beademing`: invasive ventilation > 24 h directly after the intervention (binary).
* `resp_insuff`: respiratory insufficiency requiring reintubation.
* `nierfalen`: new renal replacement therapy, or peak creatinine > 177 µmol/L together with a doubling of the baseline value (STS).
* `CVA_restletsel`: stroke with a permanent deficit, diagnosed by a neurologist. `CVA_zonder_restletsel`: neurological event without a deficit at discharge, including TIA.
* `min_temp`: minimum core temperature (°C, range 10–40). `ECC_duur`, `AoX_duur`, `circ_arrest_duur`: minutes.
* `mort_status`: 0 = alive, 1 = dead, verified preferably against the Personal Records Database (BRP).
* `overleden_zh`: death during the index admission in the operating centre (transfers to another hospital are not counted).

### 1.3 Definitions that are ambiguous and must be confirmed (flagged, not assumed)

| # | Issue | How the pipeline handles it |
|---|---|---|
| A1 | **No patient identifier.** `bronrij` is a source-row number. Repeat operations on the same patient (for example the record created by a "new heart operation during admission") cannot be linked. | Exact and near-duplicate records are counted and reported; none are removed automatically. Independence of records is assumed, and this is stated as a limitation. |
| A2 | **No dates** (operation, birth, discharge, status date). `jaar` is assumed to be the year of operation, `leeftijd` the age in completed years at operation, `opnameduur` the number of days from operation to discharge, and `overltijd_dag` the number of days from operation to death or to the last mortality-status check. | Checked empirically (for example, `overl30d` = 1 should imply `overltijd_dag` ≤ 30). For population matching, operation is placed at mid-year and exact age at `leeftijd` + 0.5. |
| A3 | **Inclusion and exclusion criteria** were applied upstream (`incl_cohort`, `excl_reden`, `s1_via_bewijs`). The extract has no arch-procedure field (`aortaboog`), so arch involvement cannot be checked directly. | `incl_cohort` defines the cohort. All exclusion reasons are tabulated (flow diagram). The `aortaboog_path` distribution among included patients is reported as a plausibility check. |
| A4 | **`*_clean` versions** of `klepchir`, `aortachir` and `circ_arrest`; the cleaning rules are unknown. | The cleaned versions are used for the primary analysis, raw versus cleaned values are cross-tabulated, and the raw versions are used in a sensitivity analysis. |
| A5 | **`pathologie_groep` and `klasse`** are derived variables whose definitions are unknown. `klasse` may encode arch extent (hemiarch or total arch), but this is unverified. | A pathology group is derived independently from the segment pathology codes (§4.3) and cross-tabulated against `pathologie_groep`. `klasse` is described only and is not used in any model. |
| A6 | **The definition of "acute"** is not in the handbook. | Study definition provided by the investigators (§4.3). Alternatives in sensitivity analyses: emergency/salvage only (`urgentie` 30/40), and non-elective (`urgentie` 20/30/40). Urgency is also reported in four levels. |
| A7 | **Changes to the handbook over time.** Version 22.2.0 (operations from 1 July 2022) added `ECMO_opname`, `postop_MI` (replacing "perioperative MI", which had a different definition) and `delier`, and split the pathology codes (chronic dissection 30 and infection 40 were added; "other" moved from 90 to 91). | Each variable's missingness is computed by year. Variables that are structurally missing in earlier years are analysed only for the years in which they were recorded, and **are not tested for trends across the definitional break**. Aneurysm and chronic dissection are combined to keep pathology comparable over time. |
| A8 | **`circ_arrest_duur`** is defined as "minutes of DHCA", while `circ_arrest` explicitly includes arrest with antegrade cerebral perfusion. It is unclear whether the duration includes periods of cerebral perfusion. | Duration is reported only among patients with arrest, and the definitional uncertainty is noted. |
| A9 | **`eGFR`**: the estimating formula is unknown. `kreatinine` is the last preoperative creatinine in µmol/L (missing by design when the patient is on dialysis). | `eGFR` is used as provided. The extract has no dialysis indicator, so patients on dialysis cannot be identified. |
| A10 | **Censoring date.** Without the status date, the date of administrative censoring is unknown, so loss to follow-up cannot be distinguished from administrative censoring. | Median follow-up is estimated with the reverse Kaplan–Meier method. Follow-up times of patients alive at last contact are shown by operation year to reveal incomplete follow-up. |
| A11 | **`overl120d`** has no handbook definition. | Described only. |
| A12 | **Completeness of 2024 (and 2013)**, given registration lag. | Annual counts are inspected, and a sensitivity analysis excludes 2024. |

---

## 2. Research questions

**Primary question (Q1).** Has risk-adjusted 30-day mortality after open aortic arch surgery changed over 2013–2024?

**Secondary questions.**

* Q2. National volume: the trend in annual procedure counts and population-adjusted rates, overall and by urgency.
* Q3. Case mix and operative complexity: how patient characteristics and concomitant procedures changed over time.
* Q4. Temporal trends in in-hospital and 1-year mortality and in the major postoperative complications.
* Q5. Differences in early outcomes and survival between acute and non-acute surgery (study definition, §4.3).
* Q6. Organ protection: temporal trends in minimum temperature, depth of hypothermia, and the use and duration of circulatory arrest, together with between-centre and within-centre variation.
* Q7. Between-centre variation in risk-adjusted outcomes.
* Q8. Failure to rescue (FTR) and its trend over time.
* Q9. Short- and long-term survival, early versus late mortality, 30-day landmark and conditional survival.
* Q10. Observed survival compared with the expected survival of the Dutch general population (relative survival, SMR, excess mortality).

---

## 3. Endpoints

| Type | Endpoint | Variable / definition |
|---|---|---|
| **Primary** | 30-day all-cause mortality | `overl30d`. BRP-verified, it captures deaths after discharge and after transfer, and it does not depend on per-intervention attribution rules or discharge practice. |
| Secondary | In-hospital mortality | `overleden_zh` |
| Secondary | 1-year mortality | 1 − Kaplan–Meier survival at 365 days (handles censoring). Binary `overl1j` is used only among patients with ≥ 365 days of potential follow-up or who died. |
| Secondary | Stroke with a permanent deficit | `CVA_restletsel` |
| Secondary | Stroke/TIA without a deficit | `CVA_zonder_restletsel` |
| Secondary | Renal failure | `nierfalen` |
| Secondary | Prolonged ventilation (> 24 h) | `beademing` |
| Secondary | Reintubation for respiratory insufficiency | `resp_insuff` |
| Secondary | ICU readmission | `heropname_IC` |
| Secondary | Major complication (composite) | Any of: `CVA_restletsel`, `nierfalen`, `resp_insuff`, `gastroint_compl` ∈ {1, 10, 90}, `vasc_compl_opname` = 10. Every component was recorded throughout 2013–2024. Missing only when no component is positive and at least one is missing. |
| Secondary | Hospital length of stay | `opnameduur` (days) |
| Secondary | Failure to rescue | 30-day death among patients with ≥ 1 major complication |
| Secondary | All-cause mortality (long term) | `overltijd_dag`, `mort_status` |
| Exploratory | ECMO, postoperative MI, delirium | Recorded from mid-2022 only: descriptive, without trend tests |
| Practice | Minimum temperature, hypothermia category, circulatory arrest use and duration | `min_temp`, `circ_arrest_clean`, `circ_arrest_duur` |

**Hypothermia categories** follow the international consensus (Yan et al., *Ann Cardiothorac Surg* 2013): profound ≤ 14 °C, deep 14.1–20 °C, low-moderate 20.1–24 °C, high-moderate 24.1–28 °C, mild 28.1–34 °C, normothermia > 34 °C. As requested, these are collapsed into **deep (≤ 20 °C, including profound), moderate (20.1–28 °C), mild (28.1–34 °C) and normothermia (> 34 °C)**. The five-level consensus grouping is reported in a supplementary table.

---

## 4. Analysis population and derived variables

### 4.1 Cohort

All rows with `incl_cohort` indicating inclusion. Each excluded row is counted by `excl_reden`, and a flow diagram is produced as a table. No record is dropped silently: every exclusion, recode to missing and out-of-range value is logged in `results/01_validation/`.

### 4.2 Cleaning rules

1. The code −1 is set to missing wherever the handbook defines it as "unknown".
2. Values outside the handbook ranges are set to missing and logged, for example `min_temp` outside 10–40 °C, `LVEF` outside 1–99%, `ECC_duur` outside 0–1000 min, `AoX_duur` outside 0–600 min, `circ_arrest_duur` outside 0–150 min (values > 150 are flagged but **kept**, because long arrest times exist clinically; see the sensitivity analysis), `leeftijd` outside 18–100, and `BMI` outside 12–70.
3. Logical inconsistencies (for example arrest = 0 with duration > 0, `overl30d` = 1 with `overltijd_dag` > 30, or in-hospital death while `mort_status` = alive) are tabulated. The raw data are never altered, and analysis variables follow the rules documented in `01_data_validation.R`.

### 4.3 Derived variables

* **Pathology** (from `aorta_asc_path`, `aortaboog_path`, `aorta_desc_path`, by precedence): acute dissection (20 in any segment) > chronic dissection or aneurysm (30, 10) > infection (40) > other (90/91). Acute aortic dissection is a strong determinant of risk and practice.
* **Calendar year**: `jaar` centred at 2018. Era: 2013–2016, 2017–2020, 2021–2024 (for description only).
* **Urgency**: four levels (adjustment variable).
* **Acute versus non-acute (study definition)**: **Acute**: `urgentie` 30 (emergency) or 40 (salvage), *or* acute dissection (code 20) in `aorta_asc_path`, `aortaboog_path` or `aorta_desc_path`. **Non-acute**: `urgentie` 10 (elective) or 20 (urgent) without an acute classification, *or* aneurysm (code 10) in any segment without an acute classification, also when urgency is unknown. **Not classifiable**: no acute-dissection or aneurysm code and no usable urgency code. Conflicts: acute dissection is always acute (also with `urgentie` 10/20); `urgentie` 30/40 is always acute (also with aneurysm code 10). Note that `urgentie` 20 (urgent) belongs to the non-acute group. Patients who are not classifiable are counted and reported, and are excluded from acute-versus-non-acute comparisons. The classification is based on observed data and is not imputed.
* **Operative complexity**: concomitant CABG, any valve surgery, aortic root involvement (`aorta_asc` ∈ {10, 20}), descending aorta involvement (`aorta_desc` > 0), hybrid open + endovascular repair (`aotachir_benadering` = 30, frozen elephant trunk), and other cardiac surgery.
* **LVEF categories** (EuroSCORE II): > 50, 31–50, 21–30, ≤ 20%.

---

## 5. Statistical methods

### 5.1 Descriptives (Q3)

* Table 1 by acute status (study definition): n (%), mean (SD) or median [IQR], and the number missing. Groups are compared with **standardised mean differences** (no hypothesis tests in Table 1).
* Change in characteristics over time: the per-year change from simple regression on calendar year (logistic for binary variables, giving an OR per year; linear for continuous variables, giving the mean change per year), with cluster-robust CR2 standard errors (by centre) and Satterthwaite degrees of freedom. These are descriptive.
* Missing data: the percentage missing for each variable by year and by centre (heatmaps).

### 5.2 National volume (Q2)

* Annual counts overall, acute and non-acute. Rates per million inhabitants use CBS population figures supplied in `data/external/population_nl.csv` (CBS servers are not reachable from the analysis environment).
* Poisson GLM: count ~ year + offset(log population). Overdispersion is assessed with the Pearson χ²/df statistic and a likelihood ratio test of the negative binomial against the Poisson model (boundary-corrected p-value). The negative binomial model is used if overdispersion is present.
* Non-linearity: linear trend versus a natural spline with 2 df (with only 12 data points, more flexible models are not justified), compared by likelihood ratio test and AIC.
* The result is reported as the incidence rate ratio (IRR) per year, that is, the annual percentage change. Whether the trend differs by urgency is tested with a year × urgency interaction.
* Sensitivity: excluding 2024 (possible registration lag) and excluding 2020–2021 (COVID-19 pandemic).

### 5.3 Early outcomes (Q1, Q4, Q5, Q7)

**Model.** A mixed-effects logistic regression with a random intercept for centre (`lme4::glmer`). The Laplace approximation (`nAGQ = 1`) is used throughout, so the same estimation method applies to the random-slope sensitivity models, where adaptive quadrature is not available. Fitting starts from an `nAGQ = 0` solution with the `nloptwrap` optimiser, followed by restarts and a `bobyqa` refit if convergence checks fail (see the methods report, §3.5).

**Estimand.** The adjusted odds ratio per calendar year, interpreted as the average within-cohort temporal change after accounting for changes in case mix.

**Covariates (primary model M1).** These are preoperative patient characteristics that plausibly affect both the case mix referred over time and the outcomes (confounders of the year–outcome association): age (natural spline, 3 df), sex, eGFR (natural spline, 3 df), LVEF category, previous cardiac surgery, active endocarditis, urgency (elective / urgent / emergency + salvage) and pathology (acute dissection / aneurysm–chronic dissection / other).

**Not adjusted for in M1.** Operative variables (concomitant procedures, hybrid repair, perfusion strategy). These are part of how practice changed over time, that is, mediators of the year effect. Adjusting for them would estimate a "controlled direct effect", which is not the question. They are added in sensitivity model M2.

**EuroSCORE II** is not combined with its own components (double adjustment). It is used as the only covariate (logit-transformed) in sensitivity model M3.

**Model set:** M0 (year + centre), M1 (primary), M2 (M1 + operative complexity), M3 (year + logit EuroSCORE II + centre).

**Non-linearity in year.** The linear year effect is the prespecified primary estimand. Non-linearity is tested by a likelihood ratio test of a natural spline with 3 df against the linear term. Adjusted, standardised risks by year (g-computation over the observed cohort, using the spline model) are plotted, with 95% simulation intervals.

**Sparse events.** Events per model parameter (EPV) are computed. If EPV < 10, a reduced model (year, age, sex, urgency, pathology) is used and flagged.

**Acute versus non-acute (Q5).** The adjusted OR for acute versus non-acute (primary: adjusted for age, sex, eGFR, LVEF, previous cardiac surgery and endocarditis, without pathology because pathology is part of the definition; secondary: also adjusted for pathology). An acute × year interaction tests whether the trends differ, and stratified models are also fitted.

**Centre variation (Q7).** The random intercept SD with a profile-likelihood 95% CI, the latent-scale ICC (σ²/(σ² + π²/3)) and the median odds ratio. Caterpillar plots of the centre conditional modes, and funnel plots of the observed/expected ratio against the expected number of events, where the expected count comes from M1 fixed effects without the centre effect. Centres are **anonymised (A–O)** in every output.

**Length of stay.** A linear mixed model of log(LOS) among patients discharged alive, giving ratios of geometric means. Sensitivity: a cause-specific Cox model for discharge alive, with in-hospital death as the competing event, stratified by centre.

### 5.4 Failure to rescue (Q8)

* **Denominator**: patients with ≥ 1 major complication (composite, §3). **Numerator**: 30-day deaths among them. Alternative: in-hospital death.
* Before analysis, the code verifies that the denominator equals the number of patients with the composite = 1 and that the numerator is no larger than the 30-day deaths.
* Model: a mixed logistic model within the complicated patients, with year plus the M1 covariates plus the number of distinct complications, and a centre random intercept.
* **Conditioning problem**: complication status is affected by year and by risk factors, so FTR trends can be distorted by changes in who develops (or is recorded as having) complications. To address this, an annual decomposition is reported: overall mortality = P(complication) × FTR + P(no complication) × mortality without complication. Unconditional trends in complications and mortality are also reported. Patients who die intraoperatively or very early have little opportunity for a complication to be recorded, which is a limitation that is reported.

### 5.5 Organ protection (Q6)

* Population: included patients with ECC. Temperature models are also run within patients with circulatory arrest.
* `min_temp`: a linear mixed model, `min_temp ~ year + case mix + (1 + year | centre)`. A **random slope** captures **within-centre trends**. With 15 centres, variance components are imprecise: they are reported with profile CIs, and a likelihood ratio test of random slope against intercept only uses a 50:50 χ² mixture.
* Hypothermia category: a cumulative-logit (proportional odds) mixed model (`ordinal::clmm`) with a centre random intercept. The proportional odds assumption is checked with `nominal_test` on the corresponding fixed-effects `clm`.
* Arrest use: a mixed logistic model. Arrest duration among patients with arrest: a linear mixed model of log(duration).
* **Within- versus between-centre change**: a centre fixed-effects model (year + centre dummies) estimates the pure within-centre trend, with CR2 cluster-robust inference and Satterthwaite df, which is appropriate for 15 clusters. Centre-specific slopes (BLUPs) are plotted.
* Each practice model is reported unadjusted and adjusted for case mix (age, sex, urgency, pathology, descending aorta involvement, hybrid repair).

### 5.6 Long-term survival (Q9)

* Time origin: the operation. Time: `overltijd_dag`. Event: `mort_status` = 1.
* Kaplan–Meier estimates (overall, by urgency, by pathology, by era) at 30 days and at 1, 3, 5 and 10 years. Median follow-up by reverse Kaplan–Meier.
* Cox regression: urgency and year plus the M1 covariates, **stratified by centre** (centre-specific baseline hazards, avoiding a distributional assumption for 15 clusters). Sensitivities: a gamma frailty model and cluster-robust SEs.
* **Proportional hazards**: scaled Schoenfeld residuals (`cox.zph`), global and per-term tests, and plots. Non-proportional hazards are expected for urgency and pathology (high early hazard). They are handled with **piecewise HRs** (0–30 days, 30 days–1 year, > 1 year) via `survSplit`, and with **restricted mean survival time** (RMST) differences at 5 and 10 years.
* **Early versus late**: early death (≤ 30 days) and late death (> 30 days).
* **30-day landmark** (prespecified): patients alive at day 30, with time reset to day 30. Kaplan–Meier curves and Cox models are fitted as above.
* **Conditional survival**: Kaplan–Meier among patients alive at s = 30 days, 1, 2 and 5 years, estimating the probability of surviving a further 1 and 5 years.
* **Selection effects**: patients alive at a landmark are a selected subgroup. Depletion of high-risk patients (especially after acute surgery) means that landmark HRs and conditional survival describe survivors and do not estimate a causal effect of urgency or year. Survival by operation year is also affected by shorter follow-up for recent years, so era comparisons use fixed horizons (1-year survival) and landmark analyses.

### 5.7 Population-referenced survival (Q10)

* Requires a Dutch life table by sex, single year of age and calendar year (CBS StatLine 37360ned, death probability `qx`), supplied as `data/external/lifetable_nl.csv`. If the file is absent, this module documents the gap and stops without producing estimates.
* Expected hazard: h*(t) = −log(1 − q<sub>sex, age(t), year(t)</sub>), with attained age(t) = `leeftijd` + 0.5 + t/365.25 and calendar year(t) = `jaar` + 0.5 + t/365.25 (because there are no dates). Ages above the table maximum use the highest age group. Years beyond the table use the last available year (documented).
* **Ederer II** relative survival in intervals (monthly in year 1, yearly to year 10), with observed survival from the actuarial life-table method and expected survival averaged over the patients at risk at the start of each interval. **Pohar-Perme** net survival is a sensitivity analysis, because Ederer II does not estimate net survival when excess hazard varies with age.
* **SMR** = observed deaths / expected deaths over follow-up, with exact Poisson CIs, overall and by follow-up interval (0–30 days, 30 days–1 year, 1–5 years, > 5 years), urgency, sex, age group and era.
* **Excess mortality** per 1000 person-years by interval, and Poisson excess-hazard regression (Dickman et al., 2004) for urgency, age, sex, era and follow-up interval.
* Verification: expected survival for reference profiles is recomputed by hand, and there are checks that attained age and calendar year advance correctly with follow-up.

### 5.8 Missing data

* Missingness is described by variable, year and centre.
* **Multiple imputation by chained equations** (`mice`, m = 20). Predictive mean matching for continuous variables, logistic or polytomous regression for categorical variables. Predictors are all analysis covariates, year, centre (fixed effects, appropriate for 15 clusters), the outcomes (30-day death, in-hospital death, the major complication composite, and the event indicator with the Nelson–Aalen cumulative hazard for survival). **Outcomes are not imputed.** Patients with a missing outcome are excluded from that outcome's analysis, and their number is reported.
* Pooling uses Rubin's rules (implemented in code), and the pooled variance components come from the average across imputations.
* **All adjusted analyses are run twice**: with multiple imputation (prespecified primary) and as a complete-case analysis (per model, patients with every model variable observed). Both sets are reported in full, and the key estimates are placed side by side (script 10). Additional sensitivity analysis: δ-adjusted imputation (missing not at random) for the continuous covariate with the most missing values.

### 5.9 Multiplicity

* Primary endpoint: a two-sided α of 0.05, with no adjustment.
* The year trends of the nine secondary binary early outcomes (Q4) are adjusted with **Holm's method** (family-wise error), shown alongside unadjusted values.
* Exploratory analyses (subgroups, interactions, outcomes from mid-2022 onward) use **Benjamini–Hochberg** (false discovery rate) and are labelled hypothesis-generating.
* Estimates and CIs are interpreted ahead of significance.

### 5.10 Model checking

* GLMMs: convergence and singular-fit warnings are recorded for every model. Fits are refitted with alternative optimisers (`allFit`) when a warning occurs. Binned residual plots are produced for the primary model, and discrimination (C statistic) and calibration (slope) are reported for case-mix descriptions.
* Linear mixed models: residual and random-effect Q–Q plots.
* Cox: Schoenfeld residuals, plus a martingale-residual check of the functional form of age.
* Ordinal: the proportional odds test.

### 5.11 Sensitivity analyses (prespecified)

1. Complete-case versus multiple imputation (every adjusted analysis, see §5.8).
2. Adjustment sets M0, M1, M2 and M3.
3. Clustering: random intercept, random intercept + slope for year, and a centre fixed-effects GLM with CR2/Satterthwaite inference.
4. Excluding 2024, and excluding salvage operations.
5. Alternative definitions of acute: emergency/salvage only; non-elective (`urgentie` 20/30/40).
6. Raw instead of `*_clean` versions of `klepchir`, `aortachir` and `circ_arrest`.
7. The major complication composite extended to include prolonged ventilation and stroke without a deficit (FTR sensitivity).
8. In-hospital death as the FTR numerator.
9. Survival: frailty versus stratified versus robust-SE Cox, piecewise HRs, and RMST.
10. δ-adjusted multiple imputation (missing not at random).
11. Arrest duration with values > 150 min excluded.

### 5.12 What will not be done

* No data-driven selection of covariates or cut-points, and no testing of many specifications in search of significance.
* No causal language for the comparisons of year, urgency or centre. All estimates are adjusted associations in observational data. Residual confounding (for example, the extent of arch replacement is not captured) and changes in coding over time are important limitations.

---

## 6. Software and reproducibility

R (version recorded in `results/session_info.txt`); packages `lme4`, `survival`, `mice`, `ordinal`, `clubSandwich`, `sandwich`, `data.table`, `ggplot2`, `patchwork`, `writexl`. A fixed random seed. `run_all.R` runs scripts 01–09 in order from a clean session. Raw data are read-only, and derived data are written to `data/derived/` (excluded from version control).

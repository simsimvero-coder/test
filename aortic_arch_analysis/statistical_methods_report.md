# Statistical Methods Report

This report explains the analyses implemented in scripts `01`–`09`: why each method was chosen, which assumptions are checked and how, and how the results should be interpreted. It accompanies `statistical_analysis_plan.md` (what is analysed) and `results/analysis_report.md` (the numbers, generated automatically by `09_tables_figures.R`).

> **Status.** The pipeline has been executed end to end on a **synthetic test file** with the same structure and codings as the real NHR extract (see §12). It has **not** been run on the patient data, which were not available in the development environment. No clinical result is reported or implied here. The "robust versus uncertain" assessment in §11 is therefore a set of prespecified criteria to apply once the real results exist, not a verdict.

---

## 1. Design and data

* **Design.** A nationwide retrospective cohort drawn from a mandatory quality registry (NHR), with 15 centres and operations from 2013 to 2024. The unit of analysis is the registered intervention, and each row is assumed to be a separate patient (see limitation A1).
* **Codings** come from the NHR Handbook Cardiochirurgie v22.2.1. Every coded variable is checked against the handbook. With `strict_coding = TRUE` (the default), the pipeline stops if it finds undocumented values, rather than guessing what they mean.
* **Raw data are read-only.** All cleaning happens in memory, and derived data go to `data/derived/`, which is never committed. Every recode to missing, out-of-range value, duplicate and inconsistency is written to `results/01_validation/`.

### Ambiguities that the analysis does not resolve by assumption

| Issue | Consequence for interpretation |
|---|---|
| No patient identifier | Repeat operations on the same patient cannot be linked. Duplicates are reported but never removed silently. If repeat records exist, independence is violated and standard errors may be slightly too small. |
| No dates | The operation is placed at mid-year and exact age at `leeftijd` + 0.5 for population matching. Errors in expected survival are at most half a year of age and calendar time, which is negligible for 1–10-year expected survival. |
| Upstream inclusion (`incl_cohort`) | The cohort definition was made outside this pipeline and cannot be fully verified, because there is no arch-procedure field. `aortaboog_path` is used as a plausibility check. |
| `*_clean` variables, `pathologie_groep`, `klasse` | The cleaning and derivation rules are unknown. Pathology is re-derived from the segment codes, `klasse` is described only, and the raw versions are used in sensitivity analyses. |
| Handbook revisions (v22.2.0, July 2022) | `ECMO_opname`, `postop_MI` and `delier` exist only from mid-2022 and are not trend-tested. Pathology codes 30/40/91 were added, so aneurysm and chronic dissection are combined to keep the categories stable over time. |
| **Extent of arch replacement** (hemiarch versus total arch) is not an explicit variable | This is an important unmeasured determinant of risk and practice. Combined open plus endovascular repair (frozen elephant trunk) and descending-aorta involvement act as partial proxies. Residual confounding of time trends by changing arch extent is possible. |

---

## 2. Choice of the primary endpoint

**30-day all-cause mortality** (`overl30d`) is primary, for three reasons:

* The vital status is verified against the BRP, the national population register.
* It includes deaths after discharge or transfer. NHR in-hospital mortality excludes patients transferred to another hospital, and transfer practices differ by centre and over time.
* In-hospital death is attributed to a single intervention under specific rules.

In-hospital mortality and 1-year mortality are secondary.

---

## 3. Risk-adjusted outcome models (scripts 04 and 08)

### 3.1 Why mixed-effects logistic regression

Outcomes are clustered within centres, which differ in case mix, practice and outcomes. A **random centre intercept** accounts for within-centre correlation, so standard errors are correct for patient-level covariates and the year effect. It also estimates the between-centre variance needed for Question 7.

With 15 centres, the centre variance is estimated imprecisely. It is therefore reported with a **likelihood-ratio (profile) 95% CI**, computed as described in §3.5. The **year effect is cross-checked** with a centre fixed-effects model using **CR2 cluster-robust standard errors with Satterthwaite degrees of freedom** (Bell–McCaffrey / Pustejovsky–Tipton). This small-sample correction is appropriate for few clusters, whereas conventional sandwich estimators are anti-conservative with 15 clusters. A random-slope model, (1 + year | centre), checks whether centre-specific trends inflate the precision of the average year effect.

### 3.2 Estimand and covariates (causal reasoning)

**Estimand.** The adjusted odds ratio per calendar year. It describes the temporal change in outcome among patients with the same preoperative risk profile.

**Confounders adjusted for** (primary model M1). These are preoperative characteristics whose distribution can change over time (referral and selection) and which affect outcome: age, sex, renal function (eGFR), LVEF, previous cardiac surgery, active endocarditis, urgency and aortic pathology. Age and eGFR enter as natural cubic splines (3 df, knots at fixed quantiles), so no linearity is assumed.

**Mediators not adjusted for.** Operative strategy (concomitant procedures, frozen elephant trunk, cerebral protection and temperature) is a channel through which outcomes may change over time. Adjusting for these would remove part of the change of interest. They are added in sensitivity model M2, which estimates the change *not* explained by measured operative strategy. The difference between M1 and M2 is descriptive and is not a formal mediation analysis.

**EuroSCORE II** is a composite of variables already in M1. Including both would double-adjust, so it is used alone in model M3.

**The year effect is a description, not a causal effect.** Calendar time has no single intervention. The estimate summarises everything that changed over time: surgical technique, perioperative care, patient selection, coding and registration completeness.

### 3.3 Non-linearity, standardised risks and multiplicity

* A **restricted cubic spline** for year (4 knots fixed in advance at 2014, 2017, 2020 and 2023) is compared with the linear term using a **pooled multivariate Wald test** (the D1 statistic of Li, Raghunathan and Rubin, 1991) on the two non-linear coefficients. Fixing the knots in advance avoids searching for the best-fitting shape.
* **Standardised risks by year** use g-computation: for each year, every patient's predicted risk is computed with year set to that value (keeping their own covariates and centre), and the predictions are averaged. Intervals come from simulating the fixed effects from their sampling distribution in each imputed dataset and pooling the draws. Centre effects are held at their conditional modes, which slightly understates uncertainty.
* **Multiplicity.** The primary endpoint is not adjusted. The nine secondary year trends are Holm-adjusted (family-wise error control), and exploratory analyses (interactions, subgroup trends) are BH-adjusted (false discovery rate). Interpretation focuses on effect sizes and CIs.

### 3.4 Sparse events

**Events per parameter (EPV)** is computed for each outcome.

* EPV < 10: a prespecified reduced model is used (year, age, sex, urgency, pathology).
* EPV still < 5 after reduction: descriptive results only.

The number of centres with zero events is reported for each outcome. Outcomes recorded only from mid-2022 (ECMO, postoperative MI, delirium) are expected to fall into these categories.

### 3.5 Convergence diagnostics

**Fitting strategy for GLMMs.** A fast `nAGQ = 0` fit provides starting values for the Laplace fit (`nAGQ = 1`, optimiser `nloptwrap`).

**Judging convergence warnings.** lme4's absolute-gradient check is known to give false positives with large n. A warning is therefore judged by the **scaled gradient** `solve(Hessian, gradient)`, as recommended in `?lme4::convergence`.

* If the scaled gradient is ≥ 0.002, the optimiser is restarted from the current estimates with tight tolerances (up to 5 times).
* If the warning persists, the model is refitted with `bobyqa` (`CONFIG$glmm_fallback = "bobyqa"`, the default). This is slower but the most rigorous option.

On the test data, refits changed the log-odds coefficients by less than 0.005, with identical log-likelihoods.

**Reporting.** Every fit records its optimiser, warnings, singularity and scaled gradient. These are tabulated per outcome in `model_convergence_diagnostics.csv` and summarised in the report.

**Centre-SD intervals.** The profile CI is computed by root-finding on lme4's deviance function: at each value of the SD, the fixed effects are profiled out by penalised iteratively reweighted least squares (`nAGQ = 0`), and the 95% limits are where the deviance exceeds its minimum by 3.84. lme4's general `confint(method = "profile")` was too slow for routine use (more than 10 minutes per model).

### 3.6 Model checking

* Primary model: binned residuals against fitted probabilities (Figure S4), the C statistic (conditional and fixed-effects only) and the calibration slope. These describe how well the case-mix adjustment works. They are not a prediction model.
* Length of stay: residual and random-effect Q–Q plots.

### 3.7 Acute versus elective

Urgency is not a treatment, so the comparison is **prognostic**. ORs are shown with and without adjustment for pathology. Pathology lies upstream of urgency: acute dissection makes surgery non-elective. The non-elective OR adjusted for pathology therefore answers a different question (urgency within the same pathology) and is expected to be much smaller. Year × urgency interactions and stratified trends are exploratory (BH-adjusted).

### 3.8 Length of stay

A linear mixed model of log(LOS) among patients discharged alive gives ratios of geometric means. Restricting to survivors is a form of conditioning on a post-operative outcome. A **cause-specific Cox model for discharge alive**, with in-hospital death treated as the competing event (censored) and stratified by centre, is therefore reported alongside. A hazard ratio above 1 means faster discharge.

---

## 4. Failure to rescue (FTR)

* **Definition.** FTR = 30-day death among patients with ≥ 1 major complication. The complications are stroke with a deficit, renal failure, reintubation, a GI complication and a major vascular complication. All five were recorded throughout 2013–2024 under stable definitions.
* **Verified identities.** The denominator equals the number with the composite = 1. The numerator is no larger than all 30-day deaths. Deaths without a recorded complication are counted, because they are outside FTR by definition. The decomposition *mortality = P(C)·FTR + (1 − P(C))·mortality without C* holds exactly each year. The script stops if any of these identities fails.
* **The conditioning problem.** FTR conditions on a post-operative event that is itself influenced by year (better prevention or more complete recording) and by patient risk. Changes in FTR can therefore reflect changes in *who* has complications rather than better rescue. Three measures address this:
  1. The FTR model adjusts for preoperative confounders and the number of complications.
  2. The annual decomposition is reported, so changes in complication incidence and in rescue can be seen side by side.
  3. Unconditional trends in complications and mortality are interpreted together with FTR.

  Without time stamps, a complication that followed a fatal event cannot be distinguished from one that preceded it. Patients dying in the operating room have little opportunity to have a complication recorded.
* **Sensitivity analyses.** In-hospital death as the numerator; an extended complication set that adds ventilation > 24 h and stroke without a deficit.

---

## 5. National volume (script 03)

* Annual counts come from 12 data points. **Poisson regression** with an offset for the CBS average population gives incidence rate ratios (the annual percentage change).
* **Overdispersion** is checked with the Pearson χ²/df statistic and a boundary-corrected likelihood ratio test of the negative binomial against Poisson. The negative binomial model is used if that test is significant. Because the test has little power with 12 points, a **dispersion-scaled (quasi-Poisson-type) CI** is always reported as well.
* **Non-linearity**: a natural spline with 2 df compared with a linear term (likelihood ratio test and AIC). More flexible models are not identifiable with 12 points.
* Acute and elective volumes are modelled separately, and a year × urgency interaction tests whether their trends differ.
* Sensitivity: excluding 2024 (registration lag) and 2020–2021 (COVID-19).
* The population figures are from CBS StatLine 03759ned (average population).

---

## 6. Organ protection (script 05)

* **Population**: patients operated with cardiopulmonary bypass. Temperature models are also run within patients with circulatory arrest.
* **Minimum temperature**: a linear mixed model with a **random intercept and random slope for year by centre**. The fixed year coefficient is the average within-centre trend, and the slope variance shows whether centres changed practice at different speeds. The random-slope model is compared with the intercept-only model by a likelihood ratio test referred to a 50:50 χ²₁/χ²₂ mixture, which accounts for the variance being tested at its boundary of zero.
* **Within- versus between-centre change.** A model with centre fixed effects uses only within-centre information and is not biased if centres with different practice changed their share of the national volume over time. Its inference uses CR2 cluster-robust SEs with Satterthwaite df. Centre-specific slopes (BLUPs, which are shrunk towards the mean) are plotted.
* **Depth of hypothermia**: a cumulative-logit (proportional odds) mixed model with a centre random intercept (`ordinal::clmm`). The proportional odds assumption is tested with `nominal_test` on the fixed-effects model. If it fails, the binary model for deep hypothermia (≤ 20 °C) is the reference result.
* **Use of arrest**: a mixed logistic model (random intercept, and intercept + slope), with the median odds ratio and ICC.
* **Arrest duration**: a linear mixed model of log(minutes), giving a ratio per year. Durations over 150 min are kept, flagged, and excluded in a sensitivity analysis.
* Every practice model is reported **unadjusted** (the national trend in practice) and **adjusted for case mix** (age, sex, urgency, pathology, descending aorta involvement, frozen elephant trunk), that is, the trend in practice for comparable patients.

---

## 7. Survival (script 06)

* **Kaplan–Meier** with log–log CIs. **Median follow-up** by reverse Kaplan–Meier. Follow-up among patients alive at last contact is shown by operation year to reveal incomplete follow-up (the censoring date is not in the extract).
* **Cox regression stratified by centre.** Stratification allows each centre its own baseline hazard without a distributional assumption, which is robust with 15 centres. Sensitivity models: gamma frailty, and cluster-robust SEs. The robust SEs are expected to be too small with 15 clusters, and are reported only for comparison.
* **Proportional hazards**: scaled Schoenfeld residual tests (global and per term) and plots. The functional form of age is checked with martingale residuals.
* **Non-proportional hazards** are expected for urgency and acute dissection, because the excess hazard is concentrated early. They are handled with:
  1. **Piecewise Cox models** with separate HRs for 0–30 days, 30 days–1 year and > 1 year.
  2. **RMST differences** at 5 and 10 years. These are assumption-free summaries in years of life gained or lost within the horizon. They are calculated only if both groups are followed up to that horizon.
* **30-day landmark**: patients alive at day 30, with follow-up restarted at day 30. This removes the immortal-time and early-hazard components.
* **Conditional survival**: Kaplan–Meier among patients alive at s (0 days, 30 days, 1, 2 and 5 years), estimating survival over a further 1 and 5 years.
* **Selection effects (essential for interpretation).** Patients alive at a landmark are a selected subset. Non-elective patients who survive 30 days have passed a strong selection filter (depletion of susceptibles). A landmark HR near 1 does not mean urgency is harmless; it means that, among survivors, later risk is similar. Neither landmark HRs nor conditional survival are causal effects. Survival by era is limited by shorter follow-up for recent operations, so era comparisons use fixed horizons (1-year mortality) and landmarks rather than the full curves.

---

## 8. Population-referenced survival (script 07)

* **Life table**: CBS StatLine 37360ned, period life tables by sex, single year of age (0–98, 99+) and calendar year (1950–2025), using the annual death probability *q*. It is converted to a constant annual hazard *h = −log(1 − q)* within each age and year.
* **Time scales.** On a weekly grid during follow-up, attained age = age at operation + 0.5 + t and calendar year = operation year + 0.5 + t. Both are updated as follow-up proceeds. Ages ≥ 99 use the 99+ row, and years after 2025 use 2025 (the number of patients affected is logged).
  * **Verification**: for four reference profiles, 5-year expected survival is recomputed by hand from the life table, and the script stops if the difference is ≥ 0.002. On the test data, the difference was below 1e-15.
  * One caveat: the CBS table indexes age as "age on 31 December", whereas the cohort has age at operation. The resulting error of about half a year in age is ignored.
* **Ederer II relative survival** (as requested). Observed interval survival uses the actuarial method (censored patients count as at risk for half an interval). Expected interval survival is the mean individual expected survival among patients at risk at the start of each interval. Intervals are monthly in year 1, then yearly. Greenwood-type CIs use the log scale.
* **Pohar-Perme net survival** is the sensitivity analysis. It weights each patient by the inverse of their expected survival, so its estimate of net survival is unbiased even when excess hazard varies with age, whereas Ederer II is not. If the two diverge, Pohar-Perme is the better estimate of net survival.
* **SMR** = observed/expected deaths with exact Poisson CIs, overall and by group. By follow-up interval, the SMR is computed among patients entering each interval (a conditional SMR). **Excess mortality** is reported per 1000 person-years.
* **Excess-hazard regression**: a Poisson GLM with the link log(μ − d*) (Dickman et al., 2004), using aggregated person-time by follow-up interval, urgency, age group, sex and era. It gives excess hazard ratios.
* **Interpretation.** The general population includes people with cardiovascular disease, but it is not a matched control group. Patients selected for elective surgery may be healthier than average for their age (a healthy-selection effect, which biases the SMR downwards), while patients with aortic disease carry the burden of their disease (biasing it upwards).

---

## 9. Missing data

* **Description**: missingness by variable, year and centre (heatmaps), plus checks for variables that are structurally missing in some years.
* **Assumption**: missing at random (MAR) given the variables in the imputation model. This is plausible because missingness in registries depends mainly on centre and period, and both are included as predictors.
* **Method**: multiple imputation by chained equations (`mice`; m = 20 by default) with predictive mean matching for continuous and binary variables, and polytomous regression for urgency and pathology.
  * Predictors: all model covariates, operative variables, year, **centre as fixed effects** (feasible and appropriate with 15 centres), and the outcomes (30-day death, in-hospital death, the major complication composite, survival status and the Nelson–Aalen cumulative hazard).
  * **Outcome values are never analysed as imputed** ("impute, then delete"). Patients with a missing outcome are excluded from that outcome's model, and their number is reported.
* **Pooling**: Rubin's rules with Barnard–Rubin-type degrees of freedom and the fraction of missing information, all implemented in code. Multivariate tests use D1.
* **Diagnostics**: trace plots, observed versus imputed densities, and mice's logged events.
* **Sensitivity**:
  * Complete-case analysis. It is unbiased if missingness depends only on covariates and not on the outcome.
  * **δ-adjusted imputation (missing not at random).** The imputed values of the continuous covariate with the most missing values (eGFR or LVEF) are shifted by ±10 units, and the primary model is refitted. This shows how far the year effect depends on the MAR assumption.

---

## 10. Sensitivity analyses: design logic

Each sensitivity analysis targets a specific threat:

| Threat | Analysis |
|---|---|
| Missing-data mechanism | Complete-case analysis, δ-adjusted imputation |
| Confounding versus mediation | M0, M2, M3 |
| Misspecified centre structure, or few clusters | Random slope; fixed effects with CR2 |
| Incomplete registration in the latest year | Excluding 2024 |
| Extreme-risk outliers | Excluding salvage operations |
| Definition of "acute" | Emergency/salvage versus non-elective |
| Unknown cleaning rules | Raw versus `*_clean` variables |
| FTR definition | Extended complication set; in-hospital numerator |
| Non-proportional hazards | Piecewise Cox, RMST, landmark |
| Clustering in Cox models | Stratified, frailty and robust-SE models |
| Ederer II bias | Pohar-Perme |
| Implausibly long arrest times | Excluding durations > 150 min |

None of these are run to search for significance. All of them are reported (`08_sensitivity`, Figure S1).

---

## 11. How results should be interpreted

### Association versus causation

Every estimate in this study is an **adjusted association from observational data**. Calendar year, urgency and centre are not interventions.

* **Year effects** describe change over time after accounting for measured case mix. They cannot be attributed to any specific innovation.
* **Centre effects** describe variation remaining after case-mix adjustment. They may reflect unmeasured case mix (notably the extent of arch replacement), coding practice or quality of care. Conditional modes are shrunk towards the mean, and centres must not be ranked on them.
* **Practice trends** (temperature, arrest) describe what changed. Any link between practice and outcome would need a separate, causally designed analysis with an explicit treatment-assignment model. That analysis is not part of this pipeline.

### Prespecified criteria for calling a finding robust (to apply to the real results)

A finding is called **robust** if all of the following hold:

1. Its 95% CI excludes the null in the primary analysis (Holm-adjusted for secondary outcomes).
2. Its direction and approximate magnitude (an OR within about ±25% of the primary estimate on the log scale) persist across complete-case, M0/M2/M3, random-slope and fixed-effects/CR2 models, the exclusion of 2024, and the δ-adjusted imputation.
3. Its models converged without unresolved warnings and had adequate EPV.
4. It is not explained by a known coding change (§1).

A finding is **uncertain** if it fails any of these: in particular, if it depends on the adjustment set, on one period, on imputation assumptions, or on a reduced model with low EPV.

**Expected sources of uncertainty, whatever the data show:**

* centre-variance estimates (only 15 clusters);
* non-linearity tests (only 12 years);
* trends in the post-2022 outcomes (not tested);
* FTR trends (the conditioning problem);
* long-term survival for recent eras (immature follow-up);
* SMR in the first 30 days (very few expected deaths, giving very large and unstable ratios).

---

## 12. Execution and verification

* **Development environment**: R 4.3.3 on Linux with Ubuntu-packaged CRAN packages. Versions are written to `results/session_info.txt` on every run.
* **Test run**: `tests/run_synthetic_test.sh` generates a **synthetic** extract (5,303 rows: 5,003 included, 300 excluded, 15 centres, 2013–2024) with all 60 columns and NHR codings. It deliberately includes −1 codes, out-of-range values, duplicates, structurally missing post-2022 variables and a Dutch-format CSV (`;` separator, `,` decimal). All nine scripts were executed on it.
  * The checks built into the code passed: FTR identities, time-scale verification, and consistency between tables, figures and models.
  * The models recovered the effects programmed into the synthetic data, for example a temperature trend of about +0.45 °C per year and a centre SD for 30-day mortality of 0.30 that lies within the estimated profile CI.
  * These outputs exist only to verify the code and are **not** study results.
* **Runtime**: with m = 20 imputations, script 04 is the slowest step. Expect it to take a few hours on a 4-core machine with the default `bobyqa` fallback, and less with `ARCH_GLMM_FALLBACK=none` (at the cost of possibly retaining fits with residual gradient warnings, which are reported).
* **What remains to be done on the real data**:
  1. Run `01_data_validation.R` and review `validation_report.md`, coding violations, `raw_vs_clean`, pathology cross-tabulations and the inclusion table.
  2. Confirm the assumptions in §1 with the data provider.
  3. Run the full pipeline.
  4. Read `analysis_report.md` and the consistency checks.
  5. Apply the robustness criteria in §11.

---

## References

* Austin PC, Merlo J. Intermediate and advanced topics in multilevel logistic regression analysis. *Stat Med* 2017;36:3257–77.
* Bell RM, McCaffrey DF. Bias reduction in standard errors for linear regression with multi-stage samples. *Surv Methodol* 2002;28:169–81.
* Dickman PW, Sloggett A, Hills M, Hakulinen T. Regression models for relative survival. *Stat Med* 2004;23:51–64.
* Ederer F, Heise H. Instructions to IBM 650 programmers in processing survival computations. NCI, 1959.
* Hickey GL, et al. Statistical and data reporting guidelines for the EJCTS and ICVTS. *Eur J Cardiothorac Surg* 2015;48:180–93.
* Li KH, Raghunathan TE, Rubin DB. Large-sample significance levels from multiply imputed data. *JASA* 1991;86:1065–73.
* Perme MP, Stare J, Estève J. On estimation in relative survival. *Biometrics* 2012;68:113–20.
* Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust variance estimation. *J Bus Econ Stat* 2018;36:672–83.
* Silber JH, et al. Hospital and patient characteristics associated with death after surgery: a study of adverse occurrence and failure to rescue. *Med Care* 1992;30:615–29.
* van Buuren S. *Flexible Imputation of Missing Data*, 2nd ed. CRC Press, 2018.
* van Houwelingen HC. Dynamic prediction by landmarking in event history analysis. *Scand J Stat* 2007;34:450–68.
* Yan TD, et al. Consensus on hypothermia in aortic arch surgery. *Ann Cardiothorac Surg* 2013;2:163–8.

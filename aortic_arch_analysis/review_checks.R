# =============================================================================
# review_checks.R  (standalone; run after run_all.R)
# Checks and small extra analyses requested in the review of paper 1.
#   A. Acute vs non-acute odds ratios without EuroSCORE II (own sample and M1 sample)
#   B. Flow diagram counts
#   C. Missing arrest duration: all patients versus patients with circulatory arrest
#   D. Concomitant valve surgery versus aortic valve surgery (cross-tabulation)
#   E. Unrounded values behind the text and tables
#   F. Extended adjustment, one-year mortality: sample versus covariate effect
#   G. Survival, relative survival, SMR and RMST with 95% confidence intervals
#   H. Sensitivity analyses for in-hospital mortality and the primary composite
#   I. Practice per centre (ranges for the text of Figure 2D-F)
# Reads the pipeline output and the derived data; writes only aggregate numbers.
# Output: results/<method>/review_checks/*.csv and overzicht_review.txt
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival) })
SUB <- "review_checks"
options(width = 220)
log_msg("Start review_checks")

cohort <- load_derived("cohort.rds")
specs  <- load_derived("specs.rds")
cc     <- analysis_datasets()
D1_END <- CONFIG$d1y_last_year
OUTS <- c(comp_primary = "Primary composite outcome", dzh = "In-hospital death", d30 = "30-day death",
          d1y = "One-year death", stroke_def = "Stroke with permanent deficit",
          stroke_nodef = "Stroke or TIA without deficit", renal = "Renal failure", gi = "GI complication",
          vent24 = "Ventilation >24 h", reintub = "Reintubation", icu_readm = "ICU readmission",
          major_compl = "Major complication")
OUTS <- OUTS[names(OUTS) %in% names(cohort)]
safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })
rd <- function(sub, f) { p <- file.path(res_path(sub), f); if (file.exists(p)) fread(p) else { log_msg("  not found: ", p); NULL } }
out_txt <- file.path(out_dir(SUB), "overzicht_review.txt")
txt <- character(0)
say <- function(title, x = NULL) {
  txt <<- c(txt, "", paste0("==== ", title, " ===="))
  if (!is.null(x)) txt <<- c(txt, capture.output(print(x, row.names = FALSE, nrows = 500)))
}
fit_term <- function(f, dd, term) {
  r <- safe(fit_pooled(f, list(dd)), deparse(f)[1])
  if (is.null(r)) return(data.table(n = NA_integer_, events = NA_integer_, or = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  tt <- term; p <- r$pooled[r$pooled$term == tt]
  data.table(n = r$diag$nobs, events = r$diag$n_events, or = exp(p$estimate), lower = exp(p$lower),
             upper = exp(p$upper), p = p$p)
}
in_m1 <- function(d, v) complete.cases(d[, .SD, .SDcols = intersect(all.vars(build_formula(v, c(COVARS_M1, "acute"), specs)), names(d))])

# ---- A. Acute vs non-acute without EuroSCORE II -----------------------------------------------
log_msg("A. acute vs non-acute without EuroSCORE II")
NO_ES <- c("age", "female", "prior_cs")
A <- rbindlist(lapply(names(OUTS), function(v) {
  end <- if (v == "d1y") D1_END else max(CONFIG$study_years)
  d <- copy(cc[[1]])[!is.na(get(v)) & year <= end]
  f_m1 <- build_formula(v, c(COVARS_M1, "acute"), specs)
  f_no <- build_formula(v, c(NO_ES, "acute"), specs)
  dm <- d[in_m1(d, v)]
  rbind(cbind(outcome = v, label = OUTS[[v]], model = "M1 (with EuroSCORE II)", fit_term(f_m1, d, "acute")),
        cbind(outcome = v, label = OUTS[[v]], model = "Without EuroSCORE II, all with age/sex/redo", fit_term(f_no, d, "acute")),
        cbind(outcome = v, label = OUTS[[v]], model = "Without EuroSCORE II, M1 sample", fit_term(f_no, dm, "acute")))
}))
write_out(A, SUB, "A_acute_vs_nonacute_without_euroscore")
say("A. OR acuut vs niet-acuut: M1, zonder EuroSCORE II (eigen cohort), zonder EuroSCORE II (M1-cohort)",
    A[, .(outcome, model, n, events, OR = sprintf("%.2f (%.2f-%.2f)", or, lower, upper), p = signif(p, 3))])

# ---- B. Flow diagram -----------------------------------------------------------------------------
log_msg("B. flow diagram")
fl <- rd("02_descriptives", "flow_diagram.csv")
sp <- rd("06_survival", "survival_population.csv")
d1 <- cohort[!is.na(d1y) & year <= D1_END]
B2 <- data.table(step = c("Analysis cohort (acute or non-acute)", "  acute", "  non-acute",
                          "Complete covariates for M1 (age, sex, EuroSCORE II, previous cardiac surgery)",
                          "  missing EuroSCORE II", "  missing previous cardiac surgery",
                          "Operations up to one-year cut-off year (one-year mortality cohort)",
                          "  of which one-year status known"),
                 n = c(nrow(cohort), sum(cohort$acute == 1), sum(cohort$acute == 0),
                       sum(complete.cases(cohort[, .(age, female, euroscore2, prior_cs)])),
                       sum(is.na(cohort$euroscore2)), sum(is.na(cohort$prior_cs)),
                       sum(cohort$year <= D1_END), nrow(d1)))
write_out(B2, SUB, "B_flow_analysis_steps")
say("B1. Flow (selectie, uit 02_descriptives)", fl)
say("B2. Flow (analyse)", B2)
say("B3. Survival-populatie (uit 06_survival)", sp)

# ---- C. Missing arrest duration --------------------------------------------------------------------
log_msg("C. arrest duration")
C <- cohort[, .(n = .N,
  arrest_yes = sum(circ_arrest_b %in% 1), arrest_no = sum(circ_arrest_b %in% 0), arrest_unknown = sum(is.na(circ_arrest_b)),
  duration_missing_all_patients = sum(is.na(ca_min) | !(circ_arrest_b %in% 1)),
  duration_missing_among_arrest = sum(circ_arrest_b %in% 1 & is.na(ca_min)),
  pct_missing_among_arrest = round(100 * mean(is.na(ca_min[circ_arrest_b %in% 1])), 1),
  duration_recorded_without_arrest = sum(circ_arrest_b %in% 0 & !is.na(ca_min) & ca_min > 0)),
  by = .(group = fifelse(acute == 1, "Acute", "Non-acute"))]
C <- rbind(C, cohort[, .(group = "All", n = .N,
  arrest_yes = sum(circ_arrest_b %in% 1), arrest_no = sum(circ_arrest_b %in% 0), arrest_unknown = sum(is.na(circ_arrest_b)),
  duration_missing_all_patients = sum(is.na(ca_min) | !(circ_arrest_b %in% 1)),
  duration_missing_among_arrest = sum(circ_arrest_b %in% 1 & is.na(ca_min)),
  pct_missing_among_arrest = round(100 * mean(is.na(ca_min[circ_arrest_b %in% 1])), 1),
  duration_recorded_without_arrest = sum(circ_arrest_b %in% 0 & !is.na(ca_min) & ca_min > 0))])
miss_other <- cohort[, .(temp_min = sum(is.na(temp_min)), ecc_min = sum(is.na(ecc_min)), aox_min = sum(is.na(aox_min)))]
write_out(C, SUB, "C_arrest_duration_missing")
say("C1. Stilstandsduur: ontbrekend over alle patienten versus onder patienten met stilstand", C)
say("C2. Ontbrekend, overige perfusievariabelen (alle patienten)", miss_other)

# ---- D. Valve surgery -------------------------------------------------------------------------------
log_msg("D. valve")
vv <- intersect(c("conc_valve", "aortic_valve", "klepchir", "klepchir_clean", "aklepchir"), names(cohort))
D <- cohort[, .N, by = c("acute", vv)][order(-acute, -N)]
D_conf <- cohort[, .(aortic_valve_1_but_conc_valve_0 = sum(aortic_valve %in% 1 & conc_valve %in% 0),
                     aortic_valve_1_but_conc_valve_NA = sum(aortic_valve %in% 1 & is.na(conc_valve))),
                 by = .(group = fifelse(acute == 1, "Acute", "Non-acute"))]
write_out(D, SUB, "D_valve_crosstab")
say("D1. Kruistabel klepchirurgie (conc_valve = klepchir_clean; aortic_valve = aklepchir, ruw)", D)
say("D2. Tegenstrijdig: aortaklep = 1 maar klepchirurgie = 0 of onbekend", D_conf)

# ---- E. Unrounded values --------------------------------------------------------------------------
log_msg("E. unrounded values")
say("E1. Temperatuurverschil zonder vs met stilstand (onafgerond)", rd("12_group_comparisons", "temperature_difference_no_arrest_vs_arrest.csv"))
bt <- rd("16_paper1_extra", "B_year_trends_overall_and_by_group.csv")
if (!is.null(bt)) say("E2. Jaartrends per groep (onafgerond), samengesteld eindpunt en ziekenhuissterfte",
                      bt[outcome %in% c("comp_primary", "dzh", "gi", "stroke_def")])

# ---- F. Extended adjustment, one-year mortality ---------------------------------------------------------
log_msg("F. extended adjustment")
d1c <- copy(cc[[1]])[!is.na(d1y) & year <= D1_END]
f_m1 <- build_formula("d1y", COVARS_M1, specs)
f_ex <- build_formula("d1y", COVARS_EXTENDED, specs)
f_bo <- build_formula("d1y", union(COVARS_M1, COVARS_EXTENDED), specs)
ccf <- function(f) complete.cases(d1c[, .SD, .SDcols = intersect(all.vars(f), names(d1c))])
d1c[, `:=`(s_m1 = ccf(f_m1), s_ex = ccf(f_ex))]
F1 <- d1c[, .(n = .N, deaths = sum(d1y == 1), mortality_pct = round(100 * mean(d1y), 1)), by = .(in_M1 = s_m1, in_extended = s_ex)][order(-in_M1, -in_extended)]
F2 <- d1c[s_m1 == TRUE, lapply(.SD, function(x) sum(is.na(x))), .SDcols = intersect(COVARS_EXTENDED, names(d1c))]
both <- d1c[s_m1 & s_ex]
F3 <- rbind(cbind(model = "M1, own sample", fit_term(f_m1, d1c, "year_c")),
            cbind(model = "Extended, own sample", fit_term(f_ex, d1c, "year_c")),
            cbind(model = "M1, common sample", fit_term(f_m1, both, "year_c")),
            cbind(model = "Extended, common sample", fit_term(f_ex, both, "year_c")),
            cbind(model = "M1 + extended covariates", fit_term(f_bo, d1c, "year_c")))
write_out(F1, SUB, "F_extended_overlap"); write_out(F3, SUB, "F_extended_models")
say("F0. Covariaten extended model", data.table(covariates = paste(COVARS_EXTENDED, collapse = ", ")))
say("F1. Overlap M1-cohort en extended-cohort (eenjaarssterfte)", F1)
say("F2. Ontbrekende extended-covariaten binnen het M1-cohort", F2)
say("F3. OR per jaar: cohort- versus covariaateffect", F3)

# ---- G. Survival with confidence intervals -----------------------------------------------------------
log_msg("G. survival with CIs")
km <- rd("06_survival", "km_survival_estimates.csv")
if (!is.null(km)) say("G1. Kaplan-Meier met 95%-BI", km[analysis %in% c("overall", "urgency")])
rs <- rd("07_population", "relative_survival_ederer2.csv")
if (!is.null(rs)) say("G2. Relatieve overleving (Ederer II) met 95%-BI op 1, 3 en 5 jaar",
                      rs[abs(end - 1) < 1e-6 | end %in% c(3, 5), .(group, end, observed, expected, relative, lower, upper)])
sm <- rd("07_population", "smr_and_excess_by_followup_interval.csv")
if (!is.null(sm)) say("G3. SMR met 95%-BI", sm[!grepl("Age", group), .(group, interval, observed, expected, smr, lower, upper)])
say("G4. Conditionele overleving", rd("06_survival", "conditional_survival.csv"))
say("G5. RMST", rd("06_survival", "rmst_by_acute_status.csv"))

# ---- H. Sensitivity analyses, in-hospital mortality and composite --------------------------------------
log_msg("H. sensitivity")
sens <- rd("08_sensitivity", "sensitivity_year_effect.csv")
if (!is.null(sens)) say("H. Sensitiviteitsanalyses jaartrend (08_sensitivity)", sens)

# ---- I. Practice per centre -------------------------------------------------------------------------
log_msg("I. centre practice")
say("I. Praktijk per centrum (bereiken voor de tekst bij Figuur 2D-F)", rd("16_paper1_extra", "I_centre_practice_profile.csv"))

writeLines(txt, out_txt)
cat(txt, sep = "\n")
log_msg("Done. Overview in ", out_txt)

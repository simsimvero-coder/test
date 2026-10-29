# =============================================================================
# 08_sensitivity_analyses.R
# Prespecified sensitivity analyses (statistical_analysis_plan.md, 5.11) for the
# primary endpoint and key secondary endpoints. Survival sensitivity analyses
# (frailty, robust SE, piecewise hazards, RMST) are in 06_survival_analysis.R.
# Output: results/08_sensitivity/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("08_sensitivity_analyses")
log_msg("Missing-data method: ", method_label())
SUB <- "08_sensitivity"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()   # imputed data sets or observed cohort (complete case)
specs  <- load_derived("specs.rds")

rows <- list()
add <- function(analysis, outcome, res, term = "year_c", n = NA, events = NA, note = "") {
  if (!is.null(attr(res, "n"))) { n <- attr(res, "n"); events <- attr(res, "events") }
  sel <- term; r <- res[res$term == sel]
  rows[[length(rows) + 1]] <<- data.table(analysis = analysis, outcome = outcome, term = term,
                                          or = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper),
                                          p = r$p, n = n, events = events, note = note)
}
pooled_or <- function(formula, dl) {
  r <- fit_pooled(formula, dl)
  out <- r$pooled; attr(out, "n") <- r$diag$nobs; attr(out, "events") <- r$diag$n_events
  out
}
cc_or <- function(formula, d) {
  vars <- all.vars(lme4::nobars(formula))
  dd <- d[complete.cases(d[, .SD, .SDcols = intersect(vars, names(d))])]
  f <- fit_mixed(formula, dd, binomial())
  s <- summary(f$fit)$coefficients
  list(res = data.table(term = rownames(s), estimate = s[, 1], lower = s[, 1] - 1.96 * s[, 2],
                        upper = s[, 1] + 1.96 * s[, 2], p = s[, 4]), n = nrow(dd), events = sum(dd[[vars[1]]] == 1))
}
fe_cr2_or <- function(formula, d) {
  f <- update(lme4::nobars(formula), . ~ . + centrum)
  vars <- all.vars(f)
  dd <- d[complete.cases(d[, .SD, .SDcols = intersect(vars, names(d))])]
  fit <- glm(f, family = binomial, data = dd)
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = dd$centrum, test = "Satterthwaite")
  r <- ct[ct$Coef == "year_c", ]
  out <- data.table(term = "year_c", estimate = r$beta, lower = r$beta - qt(.975, r$df_Satt) * r$SE,
                    upper = r$beta + qt(.975, r$df_Satt) * r$SE, p = r$p_Satt)
  attr(out, "n") <- nrow(dd); attr(out, "events") <- sum(dd[[vars[1]]] == 1)
  out
}

for (v in c("comp_primary", "d30", "dzh", "major_compl")) {
  log_msg("Sensitivity analyses: ", v)
  dl <- lapply(imps, function(d) d[!is.na(get(v))])
  n <- nrow(dl[[1]]); ev <- sum(dl[[1]][[v]] == 1)
  f1 <- build_formula(v, COVARS_M1, specs)
  add(paste0("M1 primary (", method_label(), ", random intercept)"), v, pooled_or(f1, dl))
  # 1. complete-case (only meaningful as a contrast to the imputation analysis)
  if (CONFIG$missing_method == "imputation") {
    cc <- cc_or(f1, cohort[!is.na(get(v))])
    add("Complete-case analysis", v, cc$res, n = cc$n, events = cc$events)
  }
  # 2. adjustment sets
  add("M0: year + centre only", v, pooled_or(build_formula(v, character(0), specs), dl), n = n, events = ev)
  add("M2: M1 + operative complexity", v, pooled_or(build_formula(v, c(COVARS_M1, COVARS_OPERATIVE), specs), dl), n = n, events = ev)
  add("M3: year + logit EuroSCORE II only", v, pooled_or(build_formula(v, "logit_es2", specs), dl), n = n, events = ev)
  add("Extended adjustment (age, sex, eGFR, LVEF, redo, endocarditis, urgency, pathology)", v,
      pooled_or(build_formula(v, COVARS_EXTENDED, specs), dl), n = n, events = ev)
  add("Age, sex, redo without EuroSCORE II (avoids double adjustment)", v,
      pooled_or(build_formula(v, c("age", "female", "prior_cs"), specs), dl), n = n, events = ev)
  # 3. clustering
  add("Random intercept + random slope for year", v,
      pooled_or(build_formula(v, COVARS_M1, specs, random = "(1 + year_c | centrum)"), dl), n = n, events = ev)
  add(if (CONFIG$missing_method == "imputation") "Centre fixed effects, CR2 cluster-robust (imputation 1)" else "Centre fixed effects, CR2 cluster-robust", v, fe_cr2_or(f1, dl[[1]]), n = n, events = ev)
  # 4. restrictions
  add(paste("Operations from", CONFIG$trend_sens_start), v,
      pooled_or(f1, lapply(dl, function(d) d[year >= CONFIG$trend_sens_start])))
  add("Excluding 2024", v, pooled_or(f1, lapply(dl, function(d) d[year != 2024])),
      n = sum(dl[[1]]$year != 2024), events = sum(dl[[1]][year != 2024][[v]] == 1))
  add("Excluding salvage operations", v, pooled_or(f1, lapply(dl, function(d) d[urg4 != "Salvage"])),
      n = sum(dl[[1]]$urg4 != "Salvage"), events = sum(dl[[1]][urg4 != "Salvage"][[v]] == 1))
  # 5. acute = emergency/salvage (comparison term)
  f_ac <- build_formula(v, c(COVARS_M1, "acute"), specs)
  add("Acute (study definition): OR acute vs non-acute", v, pooled_or(f_ac, dl), term = "acute", n = n, events = ev)
  f_es <- build_formula(v, c(COVARS_M1, "emerg_salv"), specs)
  r_es <- pooled_or(f_es, dl)
  add("Alternative: emergency/salvage only: year effect", v, r_es, n = n, events = ev)
  add("Alternative: emergency/salvage only (ignoring pathology): OR vs rest", v, r_es, term = "emerg_salv", n = n, events = ev)
  f_ne <- build_formula(v, c(COVARS_M1, "nonelective"), specs)
  add("Alternative: non-elective (urgentie 20/30/40): OR vs elective", v, pooled_or(f_ne, dl), term = "nonelective", n = n, events = ev)
  # 6. raw instead of cleaned valve variable (operative model M2)
  if (all(c("klepchir") %in% names(cohort))) {
    dl_raw <- lapply(dl, function(d) copy(d)[, conc_valve := klepchir])
    add("M2 with raw klepchir instead of klepchir_clean", v,
        pooled_or(build_formula(v, c(COVARS_M1, COVARS_OPERATIVE), specs), dl_raw), n = n, events = ev)
  }
}

# 7. MNAR delta adjustment for the continuous covariate with most missing values
cand <- CONFIG$mnar_var
miss <- vapply(cand, function(x) sum(is.na(cohort[[x]])), integer(1))
dvar <- names(which.max(miss))
if (max(miss) > 0 && CONFIG$missing_method == "imputation") {   # MNAR delta applies to imputation only
  for (delta in c(-CONFIG$mnar_delta, CONFIG$mnar_delta)) {
    log_msg("MNAR: ", dvar, " delta ", delta)
    di <- make_imputations(cohort, m = CONFIG$n_imputations, delta = delta, delta_var = dvar)
    dl <- lapply(di, function(d) d[!is.na(comp_primary)])
    add(sprintf("MNAR: imputed %s shifted by %+g", dvar, delta), "comp_primary",
        pooled_or(build_formula("comp_primary", COVARS_M1, specs), dl),
        note = sprintf("%d missing values shifted", max(miss)))
  }
}

# 8. Primary composite with the other stroke definition
alt <- if (CONFIG$composite_stroke == "any") "stroke with deficit only" else "any stroke (incl. without deficit/TIA)"
dl <- lapply(imps, function(d) d[!is.na(comp_primary_altstroke)])
add(paste("Primary composite with", alt), "comp_primary",
    pooled_or(build_formula("comp_primary_altstroke", COVARS_M1, specs), dl))

# 9. Failure to rescue with the extended complication definition (in-hospital death)
dl <- lapply(imps, function(d) d[major_compl_ext == 1 & !is.na(dzh)])
f_ftr <- build_formula("dzh", c(COVARS_M1, "n_compl_new"), specs)
add("FTR: extended complication set (broader set + ventilation >24 h, stroke without deficit)", "ftr_dzh",
    pooled_or(f_ftr, dl))

sens <- rbindlist(rows, fill = TRUE)
write_out(sens, SUB, "sensitivity_year_effect")

# 10. Practice variables: raw circulatory arrest variable; arrest duration without >150 min
prac <- list()
base <- function(d) d[ecc_any == 1 | is.na(ecc_any)]
f_ca <- as.formula("arrest ~ year_c + (1 | centrum)")
for (src in c("circ_arrest_b", "circ_arrest")) {
  if (!src %in% names(cohort)) next
  dl <- lapply(imps, function(d) base(d)[!is.na(get(src))][, arrest := get(src)])
  fp <- fit_pooled(f_ca, dl); r <- fp$pooled[term == "year_c"]
  prac[[src]] <- data.table(analysis = paste("Circulatory arrest use, variable:", src), or_or_ratio = exp(r$estimate),
                            lower = exp(r$lower), upper = exp(r$upper), p = r$p, n = fp$diag$nobs)
}
for (excl in c(FALSE, TRUE)) {
  dl <- lapply(imps, function(d) base(d)[circ_arrest_b == 1 & ca_min > 0 & (!excl | ca_min <= 150)][, log_ca := log(ca_min)])
  fp <- fit_pooled(as.formula("log_ca ~ year_c + (1 + year_c | centrum)"), dl, family = NULL); r <- fp$pooled[term == "year_c"]
  prac[[paste0("ca", excl)]] <- data.table(analysis = paste("Arrest duration (ratio per year)", if (excl) "excluding >150 min" else "all"),
                                           or_or_ratio = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper),
                                           p = r$p, n = fp$diag$nobs)
}
write_out(rbindlist(prac), SUB, "sensitivity_practice")

# Summary figure data: forest of the year effect on the primary composite
saveRDS(sens, derived_method("sensitivity.rds"))
log_msg("Done.")

# =============================================================================
# 04_clinical_outcomes.R
# Early outcomes: risk-adjusted temporal trends (mixed-effects logistic models
# with a centre random intercept, pooled over imputations), non-linearity,
# standardised risks by year, acute vs elective, centre variation, length of
# stay, and failure to rescue.
# Output: results/04_outcomes/*, data/derived/outcomes.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("04_clinical_outcomes")
SUB <- "04_outcomes"

cohort <- load_derived("cohort.rds")
imps   <- load_derived("imputed.rds")
specs  <- load_derived("specs.rds")
dv     <- load_derived("validated.rds")
STRUCT <- attr(dv, "structural_missing")
rm(dv)
yrs <- CONFIG$study_years

OUT <- data.table(
  var   = c("d30", "dzh", "d1y", "stroke_def", "stroke_nodef", "renal", "vent24",
            "reintub", "icu_readm", "major_compl", "ecmo", "pmi", "delirium"),
  label = c("30-day mortality", "In-hospital mortality", "1-year mortality (complete follow-up)",
            "Stroke with permanent deficit", "Stroke/TIA without deficit", "Renal failure",
            "Ventilation > 24 h", "Reintubation (respiratory insufficiency)", "ICU readmission",
            "Major complication (composite)", "ECMO during admission",
            "Postoperative myocardial infarction", "Delirium"),
  role  = c("primary", rep("secondary", 9), rep("exploratory", 3)))
# Variables recorded only in later years: restrict and do not test trends
OUT[, first_year := vapply(var, function(v) {
  src <- c(ecmo = "ECMO_opname", pmi = "postop_MI", delirium = "delier")[v]
  if (!is.na(src) && src %in% STRUCT$variable) STRUCT[variable == src, first_year] else min(yrs)
}, numeric(1))]
OUT[var %in% c("ecmo", "pmi", "delirium"), role := "exploratory (recorded from first_year only)"]

# ---- 1. Crude rates ---------------------------------------------------------------------
crude <- rbindlist(lapply(seq_len(nrow(OUT)), function(i) {
  v <- OUT$var[i]; d <- cohort[year >= OUT$first_year[i] & !is.na(get(v))]
  grp <- function(by) d[, .(events = sum(get(v) == 1), n = .N), by = by]
  rbind(cbind(group_type = "overall", group = "All", d[, .(events = sum(get(v) == 1), n = .N)]),
        grp("urgency2")[!is.na(urgency2), .(group_type = "urgency", group = as.character(urgency2), events, n)],
        grp("era")[, .(group_type = "era", group = as.character(era), events, n)],
        grp("year")[, .(group_type = "year", group = as.character(year), events, n)],
        fill = TRUE)[, `:=`(outcome = v, label = OUT$label[i], pct = 100 * events / n,
                            n_missing = sum(is.na(cohort[year >= OUT$first_year[i]][[v]])))]
}))
crude[, c("lower", "upper") := as.list(binom.test(events, n)$conf.int * 100), by = seq_len(nrow(crude))]
write_out(crude, SUB, "crude_rates")

# Sparse events by centre
sparse <- rbindlist(lapply(OUT$var, function(v)
  cohort[!is.na(get(v)), .(events = sum(get(v) == 1), n = .N), by = centre_lab][, outcome := v]))
write_out(sparse[, .(n_centres = .N, centres_zero_events = sum(events == 0),
                     min_events = min(events), median_events = median(events)), by = outcome],
          SUB, "sparse_events_by_centre")

# ---- 2. Risk-adjusted models ---------------------------------------------------------------
model_rows <- list(); diag_rows <- list(); gcomp_rows <- list(); urg_rows <- list()
centre_rows <- list(); fits_keep <- list()

analyse_binary <- function(v, label, first_year, test_trend = TRUE, extended = FALSE) {
  dl <- lapply(imps, function(d) d[year >= first_year & !is.na(get(v))])
  y <- dl[[1]][[v]]
  covars <- COVARS_M1
  f1 <- build_formula(v, covars, specs)
  npar <- n_fixed_par(f1, dl[[1]])
  ev_pp <- epv(y, npar)
  reduced <- ev_pp < 10
  if (reduced) {
    covars <- c("age", "female", "urg3", "path3")
    f1 <- build_formula(v, covars, specs)
    ev_pp <- epv(y, n_fixed_par(f1, dl[[1]]))
    log_msg(v, ": EPV < 10 -> reduced model (EPV now ", round(ev_pp, 1), ")")
  }
  if (sum(y == 1) < 10 || ev_pp < 5) {
    log_msg(v, ": too few events (EPV < 5 in reduced model), models not fitted")
    return(list(model = data.table(outcome = v, label = label, n = length(y), events = sum(y == 1),
                                   epv = ev_pp, note = "too few events (EPV < 5): descriptive only")))
  }
  m1 <- fit_pooled(f1, dl)
  tab <- exp_table(m1$pooled)
  tab[, `:=`(outcome = v, label = label, model = if (reduced) "M1 reduced" else "M1",
             n = length(y), events = sum(y == 1), epv = ev_pp)]
  res <- list(model = tab,
              diag = cbind(outcome = v, model = "M1", m1$diag, epv = ev_pp, reduced = reduced),
              re_sd = m1$re_sd, fits = m1$fits, data = dl, formula = f1, covars = covars)
  if (!test_trend) return(res)

  # Non-linear calendar-year effect (restricted cubic spline, 4 fixed knots)
  fs <- build_formula(v, covars, specs, year_term = YEAR_SPLINE)
  ms <- fit_pooled(fs, dl)
  nl <- pool_wald_D1(lapply(ms$fits, function(f) fixef(f$fit)),
                     lapply(ms$fits, function(f) vcov(f$fit)), YEAR_NL_TERMS)
  res$nonlinear <- cbind(outcome = v, nl)
  res$diag <- rbind(res$diag, cbind(outcome = v, model = "M1 spline", ms$diag, epv = ev_pp, reduced = reduced))
  # Standardised (g-computation) risk by year from the spline model
  gy <- gcomp_mi(lapply(ms$fits, `[[`, "fit"), dl, yrs[yrs >= first_year])
  res$gcomp <- cbind(outcome = v, label = label, gy)
  res
}

for (i in seq_len(nrow(OUT))) {
  v <- OUT$var[i]
  log_msg("Outcome: ", v)
  r <- tryCatch(analyse_binary(v, OUT$label[i], OUT$first_year[i],
                               test_trend = !grepl("exploratory", OUT$role[i])),
                error = function(e) { log_msg("ERROR ", v, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  model_rows[[v]] <- r$model
  if (!is.null(r$diag)) diag_rows[[v]] <- r$diag
  if (!is.null(r$gcomp)) gcomp_rows[[v]] <- merge(r$gcomp, r$nonlinear[, .(outcome, p_nonlinear = p)], by = "outcome")
  if (!is.null(r$fits)) fits_keep[[v]] <- r[c("fits", "data", "formula", "covars", "re_sd")]
}
models <- rbindlist(model_rows, fill = TRUE)
write_out(models, SUB, "adjusted_models_M1_all_terms")
write_out(rbindlist(diag_rows, fill = TRUE), SUB, "model_convergence_diagnostics")
gcomp <- rbindlist(gcomp_rows, fill = TRUE)
write_out(gcomp, SUB, "standardised_risk_by_year")

# Year effect summary with multiplicity adjustment
yr <- models[term == "year_c", .(outcome, label, model, n, events, epv, or_per_year = ratio,
                                 lower = ratio_lower, upper = ratio_upper, p)]
yr <- merge(yr, OUT[, .(outcome = var, role)], by = "outcome")
yr[, p_holm := NA_real_]; yr[role == "secondary", p_holm := holm(p)]
yr[, p_bh := NA_real_]; yr[grepl("exploratory", role), p_bh := bh(p)]
nlin <- unique(gcomp[, .(outcome, p_nonlinear)])
yr <- merge(yr, nlin, by = "outcome", all.x = TRUE)
setorder(yr, role, outcome)
write_out(yr, SUB, "year_trend_summary")

# ---- 3. Acute vs elective --------------------------------------------------------------------
# OR non-elective vs elective, with and without pathology (pathology lies upstream
# of urgency); year x urgency interaction; stratified year trends.
KEY <- intersect(c("d30", "dzh", "major_compl", "stroke_def", "renal"), names(fits_keep))
covars_nourg <- setdiff(COVARS_M1, "urg3")
for (v in KEY) {
  dl <- fits_keep[[v]]$data
  f_a <- build_formula(v, c(setdiff(covars_nourg, "path3"), "nonelective"), specs)
  f_b <- build_formula(v, c(covars_nourg, "nonelective"), specs)
  f_i <- build_formula(v, c(covars_nourg, "nonelective", "year_c:nonelective"), specs)
  pa <- exp_table(fit_pooled(f_a, dl)$pooled, "nonelective")[, adjustment := "without pathology"]
  pb <- exp_table(fit_pooled(f_b, dl)$pooled, "nonelective")[, adjustment := "with pathology"]
  p_int <- fit_pooled(f_i, dl)$pooled
  strat <- rbindlist(lapply(c(0, 1), function(g) {
    dg <- lapply(dl, function(d) d[nonelective == g])
    ff <- build_formula(v, setdiff(covars_nourg, if (g == 0) "path3" else character(0)), specs)
    if (sum(dg[[1]][[v]] == 1) < 20) return(NULL)
    exp_table(fit_pooled(ff, dg)$pooled, "year_c")[, stratum := if (g == 1) "Non-elective" else "Elective"]
  }))
  urg_rows[[v]] <- list(
    or = rbind(pa, pb)[, outcome := v],
    interaction = p_int[term == "year_c:nonelective"][, outcome := v],
    strata = strat[, outcome := v])
}
write_out(rbindlist(lapply(urg_rows, `[[`, "or")), SUB, "acute_vs_elective_or")
write_out(rbindlist(lapply(urg_rows, `[[`, "interaction"))[, `:=`(ratio_of_or = exp(estimate), p_bh = bh(p))],
          SUB, "acute_vs_elective_year_interaction")
write_out(rbindlist(lapply(urg_rows, `[[`, "strata")), SUB, "year_trend_by_urgency_stratum")

# ---- 4. Centre variation ----------------------------------------------------------------------
for (v in intersect(c("d30", "dzh", "major_compl", "stroke_def", "renal", "vent24"), names(fits_keep))) {
  fk <- fits_keep[[v]]
  f1 <- fk$fits[[1]]$fit; d1 <- fk$data[[1]]
  # Profile-likelihood CI for the centre SD (complete-case fit; documented)
  cc <- cohort[year >= OUT[var == v, first_year] & !is.na(get(v))]
  cc <- cc[complete.cases(cc[, .SD, .SDcols = setdiff(fk$covars, character(0))])]
  prof <- tryCatch(profile_sd_ci(fk$formula, cc), error = function(e) c(sd = NA, lower = NA, upper = NA))
  s <- fk$re_sd
  centre_rows[[v]] <- data.table(outcome = v, n_complete_case = nrow(cc), centre_sd_pooled = s,
                                 centre_sd_cc_profile = prof[["sd"]],
                                 sd_profile_lower = prof[["lower"]], sd_profile_upper = prof[["upper"]],
                                 mor = mor_from_sd(s), mor_lower = mor_from_sd(prof[["lower"]]),
                                 mor_upper = mor_from_sd(prof[["upper"]]), icc_latent = icc_logit(s))
  # Caterpillar data (conditional modes, imputation 1) and funnel data (O/E)
  re <- as.data.frame(ranef(f1, condVar = TRUE))
  key <- unique(cohort[, .(centrum = as.character(centrum), centre_lab)])
  cat_d <- merge(data.table(centrum = as.character(re$grp), effect = re$condval, se = re$condsd), key, by = "centrum")
  cat_d[, `:=`(outcome = v, or = exp(effect), lower = exp(effect - 1.96 * se), upper = exp(effect + 1.96 * se))]
  d1[, expected := predict(f1, newdata = d1, re.form = NA, type = "response")]
  fun <- d1[, .(observed = sum(get(v) == 1), expected = sum(expected), n = .N), by = centre_lab]
  fun[, `:=`(outcome = v, oe = observed / expected)]
  centre_rows[[paste0(v, "_cat")]] <- cat_d[, centrum := NULL]
  centre_rows[[paste0(v, "_fun")]] <- fun
}
write_out(rbindlist(centre_rows[!grepl("_cat|_fun", names(centre_rows))]), SUB, "centre_variation")
write_out(rbindlist(centre_rows[grepl("_cat", names(centre_rows))]), SUB, "centre_effects_caterpillar")
write_out(rbindlist(centre_rows[grepl("_fun", names(centre_rows))]), SUB, "centre_observed_expected")

# ---- 5. Primary model diagnostics -----------------------------------------------------------------
if ("d30" %in% names(fits_keep)) {
  f1 <- fits_keep$d30$fits[[1]]$fit; d1 <- fits_keep$d30$data[[1]]
  p_hat <- fitted(f1); yobs <- getME(f1, "y")
  lp <- qlogis(predict(f1, type = "response", re.form = NA))
  cstat <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
  cal <- glm(yobs ~ lp, family = binomial)
  bins <- data.table(p = p_hat, r = yobs - p_hat)[, bin := cut(rank(p), 20, labels = FALSE)]
  binned <- bins[, .(mean_fitted = mean(p), mean_resid = mean(r), se = 2 * sqrt(mean(p * (1 - p)) / .N)), by = bin]
  write_out(binned, SUB, "primary_model_binned_residuals")
  write_out(data.table(c_statistic_conditional = cstat(p_hat, yobs),
                       c_statistic_fixed_only = cstat(lp, yobs),
                       calibration_slope_fixed_only = coef(cal)[2],
                       calibration_intercept = coef(cal)[1]), SUB, "primary_model_performance")
}

# ---- 6. Length of stay ------------------------------------------------------------------------
los_dl <- lapply(imps, function(d) d[dzh == 0 & !is.na(los)][, log_los := log(pmax(los, 0.5))])
f_los <- build_formula("log_los", COVARS_M1, specs)
los_fit <- fit_pooled(f_los, los_dl, family = NULL)
los_tab <- exp_table(los_fit$pooled)[, `:=`(outcome = "LOS (discharged alive)",
                                            interpretation = "ratio of geometric means")]
# Sensitivity: cause-specific Cox model for discharge alive (in-hospital death censored)
cox_est <- lapply(imps, function(d) {
  dd <- d[!is.na(los) & !is.na(dzh)][, disch := as.integer(dzh == 0)]
  f <- as.formula(paste("Surv(pmax(los, 0.5), disch) ~ year_c +",
                        paste(covar_terms(COVARS_M1, specs), collapse = " + "), "+ strata(centrum)"))
  fit <- coxph(f, data = dd)
  list(coef(fit), vcov(fit))
})
los_cox <- exp_table(rubin_pool(lapply(cox_est, `[[`, 1), lapply(cox_est, `[[`, 2)), "year_c")
los_cox[, `:=`(outcome = "Discharge alive (cause-specific HR; >1 = shorter stay)")]
write_out(rbind(los_tab, los_cox, fill = TRUE), SUB, "length_of_stay_models")
write_out(cohort[dzh == 0 & !is.na(los), .(n = .N, median = median(los), q1 = quantile(los, .25),
                                            q3 = quantile(los, .75)), by = year][order(year)],
          SUB, "length_of_stay_by_year")
write_out(rbind(los_fit$diag[, model := "LOS LMM"]), SUB, "length_of_stay_diagnostics")
r1 <- residuals(los_fit$fits[[1]]$fit)
grDevices::png(file.path(out_dir(SUB), "los_lmm_qq.png"), 1400, 700, res = 150)
par(mfrow = c(1, 2)); qqnorm(r1, main = "LOS model residuals"); qqline(r1)
qqnorm(ranef(los_fit$fits[[1]]$fit)$centrum[, 1], main = "Centre random effects"); qqline(ranef(los_fit$fits[[1]]$fit)$centrum[, 1])
grDevices::dev.off()

# ---- 7. Failure to rescue ------------------------------------------------------------------------
MAJOR <- c("stroke_def", "renal", "reintub", "gi", "vasc_major")
ftr_d <- cohort[major_compl == 1]
check <- data.table(
  item = c("patients with major complication (denominator)", "denominator with 30-day status known",
           "deaths within 30 days among denominator (numerator)",
           "all 30-day deaths in cohort", "30-day deaths without recorded major complication",
           "denominator equals sum(major_compl == 1)", "numerator <= all 30-day deaths"),
  value = c(nrow(ftr_d), sum(!is.na(ftr_d$d30)), sum(ftr_d$d30 == 1, na.rm = TRUE),
            sum(cohort$d30 == 1, na.rm = TRUE), sum(cohort$d30 == 1 & cohort$major_compl == 0, na.rm = TRUE),
            nrow(ftr_d) == sum(cohort$major_compl == 1, na.rm = TRUE),
            sum(ftr_d$d30 == 1, na.rm = TRUE) <= sum(cohort$d30 == 1, na.rm = TRUE)))
stopifnot(check$value[6] == 1, check$value[7] == 1)
write_out(check, SUB, "ftr_numerator_denominator_check")
# Decomposition by year: mortality = P(C) * FTR + P(no C) * mortality without C
decomp <- cohort[!is.na(major_compl) & !is.na(d30),
                 .(n = .N, p_complication = mean(major_compl),
                   ftr = mean(d30[major_compl == 1]), mort_no_compl = mean(d30[major_compl == 0]),
                   mortality = mean(d30)), by = year][order(year)]
decomp[, check_identity := p_complication * ftr + (1 - p_complication) * mort_no_compl - mortality]
stopifnot(all(abs(decomp$check_identity) < 1e-10))
write_out(decomp, SUB, "ftr_decomposition_by_year")
ftr_models <- list()
for (num in c("d30", "dzh")) {
  dl <- lapply(imps, function(d) d[major_compl == 1 & !is.na(get(num))])
  f <- build_formula(num, c(COVARS_M1, "n_compl"), specs)
  npar <- n_fixed_par(f, dl[[1]]); ev <- epv(dl[[1]][[num]], npar)
  if (ev < 10) f <- build_formula(num, c("age", "female", "urg3", "path3", "n_compl"), specs)
  fp <- fit_pooled(f, dl)
  ftr_models[[num]] <- exp_table(fp$pooled)[, `:=`(numerator = num, n = nrow(dl[[1]]),
                                                   events = sum(dl[[1]][[num]] == 1), epv = ev,
                                                   reduced = ev < 10, centre_sd = fp$re_sd)]
}
write_out(rbindlist(ftr_models), SUB, "ftr_models")

saveRDS(list(OUT = OUT, crude = crude, models = models, year = yr, gcomp = gcomp,
             centre = rbindlist(centre_rows[!grepl("_cat|_fun", names(centre_rows))]),
             caterpillar = rbindlist(centre_rows[grepl("_cat", names(centre_rows))]),
             funnel = rbindlist(centre_rows[grepl("_fun", names(centre_rows))]),
             urgency = rbindlist(lapply(urg_rows, `[[`, "or")),
             los = rbind(los_tab, los_cox, fill = TRUE), ftr = rbindlist(ftr_models), ftr_decomp = decomp),
        derived_path("outcomes.rds"))
log_msg("Done.")

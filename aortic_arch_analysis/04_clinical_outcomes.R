# =============================================================================
# 04_clinical_outcomes.R
# Early outcomes: risk-adjusted temporal trends (mixed-effects logistic models
# with a centre random intercept, pooled over imputations), non-linearity,
# standardised risks by year, acute vs non-acute, centre variation, length of
# stay, and failure to rescue.
# Output: results/04_outcomes/*, data/derived/outcomes.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("04_clinical_outcomes")
log_msg("Missing-data method: ", method_label())
SUB <- "04_outcomes"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()   # imputed data sets or observed cohort (complete case)
specs  <- load_derived("specs.rds")
dv     <- load_derived("validated.rds")
STRUCT <- attr(dv, "structural_missing")
rm(dv)
yrs <- CONFIG$study_years

# Primary endpoint (decided explicitly, continuity with paper 1): composite of in-hospital
# death, stroke, renal failure or GI complication. 30-day mortality is a key secondary endpoint.
stroke_lab <- if (CONFIG$composite_stroke == "any") "any stroke" else "stroke with deficit"
OUT <- data.table(
  var   = c("comp_primary", "d30", "dzh", "d1y", "stroke_def", "stroke_nodef", "renal", "gi", "vent24",
            "reintub", "icu_readm", "major_compl", "major_compl_new", "ecmo", "pmi", "delirium"),
  label = c(sprintf("Primary composite (in-hospital death, %s, renal failure, GI complication)", stroke_lab),
            "30-day mortality", "In-hospital mortality",
            sprintf("1-year mortality (operations %d-%d)", min(yrs), CONFIG$d1y_last_year),
            "Stroke with permanent deficit", "Stroke/TIA without deficit", "Renal failure",
            "GI complication", "Ventilation > 24 h", "Reintubation (respiratory insufficiency)", "ICU readmission",
            sprintf("Major complication (%s, renal failure, GI)",
                    if (CONFIG$major_stroke == "any") "any stroke" else "stroke with deficit"),
            "Major complication, broader definition (+ reintubation, major vascular)",
            "ECMO during admission", "Postoperative myocardial infarction", "Delirium"),
  role  = c("primary", rep("secondary", 11), "supplementary", rep("exploratory", 3)))
# Variables recorded only in later years: restrict and do not test trends
OUT[, first_year := vapply(var, function(v) {
  src <- c(ecmo = "ECMO_opname", pmi = "postop_MI", delirium = "delier")[v]
  if (!is.na(src) && src %in% STRUCT$variable) STRUCT[variable == src, first_year] else min(yrs)
}, numeric(1))]
OUT[var %in% c("ecmo", "pmi", "delirium"), role := "exploratory (recorded from first_year only)"]
# Last operation year analysed: 1-year mortality only for years with complete potential
# 1-year follow-up (avoids selective inclusion of early deaths from the last year)
OUT[, last_year := max(yrs)]
OUT[var == "d1y", last_year := CONFIG$d1y_last_year]

# ---- 1. Crude rates ---------------------------------------------------------------------
crude <- rbindlist(lapply(seq_len(nrow(OUT)), function(i) {
  v <- OUT$var[i]; d <- cohort[year >= OUT$first_year[i] & year <= OUT$last_year[i] & !is.na(get(v))]
  grp <- function(by) d[, .(events = sum(get(v) == 1), n = .N), by = by]
  rbind(cbind(group_type = "overall", group = "All", d[, .(events = sum(get(v) == 1), n = .N)]),
        grp("urgency2")[!is.na(urgency2), .(group_type = "acute status", group = as.character(urgency2), events, n)],
        grp("era")[, .(group_type = "era", group = as.character(era), events, n)],
        grp("year")[, .(group_type = "year", group = as.character(year), events, n)],
        fill = TRUE)[, `:=`(outcome = v, label = OUT$label[i], pct = 100 * events / n,
                            n_missing = sum(is.na(cohort[year >= OUT$first_year[i] & year <= OUT$last_year[i]][[v]])))]
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
model_rows <- list(); diag_rows <- list(); gcomp_rows <- list(); urg_rows <- list(); from_rows <- list()
centre_rows <- list(); fits_keep <- list()

analyse_binary <- function(v, label, first_year, test_trend = TRUE, extended = FALSE, last_year = max(yrs)) {
  dl_all <- lapply(imps, function(d) d[year >= first_year & year <= last_year & !is.na(get(v))])
  covars <- COVARS_M1
  f1 <- build_formula(v, covars, specs)
  dl <- complete_for(dl_all, f1)          # complete-case: rows with all model variables
  y <- dl[[1]][[v]]
  npar <- n_fixed_par(f1, dl[[1]])
  ev_pp <- epv(y, npar)
  reduced <- ev_pp < 10
  if (reduced) {
    covars <- COVARS_REDUCED
    f1 <- build_formula(v, covars, specs)
    dl <- complete_for(dl_all, f1); y <- dl[[1]][[v]]
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
             n = length(y), events = sum(y == 1), epv = ev_pp,
             n_with_outcome = nrow(dl_all[[1]]), missing_method = CONFIG$missing_method)]
  res <- list(model = tab,
              diag = cbind(outcome = v, model = "M1", m1$diag, epv = ev_pp, reduced = reduced),
              re_sd = m1$re_sd, fits = m1$fits, data = dl, formula = f1, covars = covars)
  if (!test_trend) return(res)

  # Sensitivity: year trend from CONFIG$trend_sens_start onwards (same model)
  if (first_year < CONFIG$trend_sens_start) {
    dls <- lapply(dl, function(d) d[year >= CONFIG$trend_sens_start])
    ys <- dls[[1]][[v]]
    if (sum(ys == 1) >= 10) {
      ms15 <- fit_pooled(f1, dls)
      res$from <- exp_table(ms15$pooled, "year_c")[, `:=`(outcome = v, label = label,
        from_year = CONFIG$trend_sens_start, n = length(ys), events = sum(ys == 1))]
    }
  }

  # Non-linear calendar-year effect (restricted cubic spline, 4 fixed knots)
  fs <- build_formula(v, covars, specs, year_term = YEAR_SPLINE)
  ms <- fit_pooled(fs, dl)
  nl <- pool_wald_D1(lapply(ms$fits, function(f) fixef(f$fit)),
                     lapply(ms$fits, function(f) vcov(f$fit)), YEAR_NL_TERMS)
  res$nonlinear <- cbind(outcome = v, nl)
  res$diag <- rbind(res$diag, cbind(outcome = v, model = "M1 spline", ms$diag, epv = ev_pp, reduced = reduced))
  # Standardised (g-computation) risk by year from the spline model
  gy <- gcomp_mi(lapply(ms$fits, `[[`, "fit"), dl, yrs[yrs >= first_year & yrs <= last_year])
  res$gcomp <- cbind(outcome = v, label = label, gy)
  res
}

for (i in seq_len(nrow(OUT))) {
  v <- OUT$var[i]
  log_msg("Outcome: ", v)
  r <- tryCatch(analyse_binary(v, OUT$label[i], OUT$first_year[i], last_year = OUT$last_year[i],
                               test_trend = !grepl("exploratory", OUT$role[i])),
                error = function(e) { log_msg("ERROR ", v, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  model_rows[[v]] <- r$model
  if (!is.null(r$diag)) diag_rows[[v]] <- r$diag
  if (!is.null(r$gcomp)) gcomp_rows[[v]] <- merge(r$gcomp, r$nonlinear[, .(outcome, p_nonlinear = p)], by = "outcome")
  if (!is.null(r$fits)) fits_keep[[v]] <- r[c("fits", "data", "formula", "covars", "re_sd")]
  if (!is.null(r$from)) from_rows[[v]] <- r$from
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
yr_from <- rbindlist(from_rows, fill = TRUE)
if (nrow(yr_from)) {
  yr_from <- merge(yr_from, OUT[, .(outcome = var, role)], by = "outcome")
  yr_from[, p_holm := NA_real_]; yr_from[role == "secondary", p_holm := holm(p)]
  write_out(yr_from[, .(outcome, label, role, from_year, n, events, or_per_year = ratio,
                        lower = ratio_lower, upper = ratio_upper, p, p_holm)],
            SUB, paste0("year_trend_summary_from", CONFIG$trend_sens_start))
}

# ---- 3. Acute vs non-acute (study definition) --------------------------------------------------------------------
# OR acute vs non-acute (adjusted for the M1 covariates; and additionally for pathology)
# and year x acute interaction, for the main endpoints.
KEY <- intersect(c("comp_primary", "d30", "dzh", "major_compl", "stroke_def", "renal"), names(fits_keep))
for (v in KEY) {
  dl <- fits_keep[[v]]$data
  f_a <- build_formula(v, c(COVARS_M1, "acute"), specs)
  f_b <- build_formula(v, c(COVARS_M1, "path3", "acute"), specs)
  f_i <- build_formula(v, c(COVARS_M1, "acute", "year_c:acute"), specs)
  pa <- exp_table(fit_pooled(f_a, dl)$pooled, "acute")[, adjustment := "without pathology"]
  pb <- exp_table(fit_pooled(f_b, dl)$pooled, "acute")[, adjustment := "with pathology"]
  p_int <- fit_pooled(f_i, dl)$pooled
  urg_rows[[v]] <- list(or = rbind(pa, pb)[, outcome := v],
                        interaction = p_int[term == "year_c:acute"][, outcome := v])
}
write_out(rbindlist(lapply(urg_rows, `[[`, "or")), SUB, "acute_vs_nonacute_or")
write_out(rbindlist(lapply(urg_rows, `[[`, "interaction"))[, `:=`(ratio_of_or = exp(estimate), p_bh = bh(p))],
          SUB, "acute_vs_nonacute_year_interaction")

# ---- 3b. Year trends per group (acute / non-acute) x start year (2013 / 2015) -----------------------------------
# For every outcome: adjusted OR per calendar year within the acute and the non-acute group,
# starting in the first study year and in CONFIG$trend_sens_start (1-year mortality ends in
# CONFIG$d1y_last_year). Same covariates (M1); reduced model if EPV < 10; not fitted if EPV < 5.
group_trend <- function(dl_base, v, label, group, start, end) {
  g <- if (group == "Acute") 1L else 0L
  dg0 <- lapply(dl_base, function(d) d[acute == g & year >= start & year <= end & !is.na(get(v))])
  f <- build_formula(v, COVARS_M1, specs); model <- "M1"
  dg <- complete_for(dg0, f); y <- dg[[1]][[v]]; ev <- sum(y == 1)
  ep <- if (ev > 0) epv(y, n_fixed_par(f, dg[[1]])) else 0
  if (ep < 10) { f <- build_formula(v, COVARS_REDUCED, specs); model <- "reduced"
                 dg <- complete_for(dg0, f); y <- dg[[1]][[v]]; ev <- sum(y == 1)
                 ep <- if (ev > 0) epv(y, n_fixed_par(f, dg[[1]])) else 0 }
  base <- data.table(outcome = v, label = label, group = group, start_year = start, end_year = end,
                     n = length(y), events = ev)
  if (ev < 10 || ep < 5)
    return(cbind(base, epv = ep, model = "not fitted (EPV < 5)", or_per_year = NA_real_,
                 lower = NA_real_, upper = NA_real_, p = NA_real_))
  fp <- fit_pooled(f, dg)
  e <- exp_table(fp$pooled, "year_c")
  cbind(base, epv = ep, model = model, or_per_year = e$ratio, lower = e$ratio_lower,
        upper = e$ratio_upper, p = e$p, n_convergence_warning = fp$diag$n_convergence_warning)
}
gt_rows <- list()
for (v in intersect(OUT[!grepl("exploratory", role), var], names(fits_keep))) {
  o <- OUT[var == v]
  starts <- unique(c(o$first_year, max(o$first_year, CONFIG$trend_sens_start)))
  for (grp in c("Acute", "Non-acute")) for (st in starts) {
    log_msg("Group trend: ", v, " / ", grp, " / from ", st)
    gt_rows[[length(gt_rows) + 1]] <- tryCatch(group_trend(fits_keep[[v]]$data, v, o$label, grp, st, o$last_year),
      error = function(e) data.table(outcome = v, label = o$label, group = grp, start_year = st,
                                     end_year = o$last_year, model = paste("error:", conditionMessage(e))))
  }
}
# Failure to rescue per group and start year (paper-1 complication definition; in-hospital and 30-day death)
for (num in c("dzh", "d30")) {
  lab <- sprintf("Failure to rescue (%s among patients with major complication)",
                 if (num == "dzh") "in-hospital death" else "30-day death")
  dl_f <- lapply(imps, function(d) d[major_compl == 1 & !is.na(get(num))][, ftr_y := get(num)])
  dl_f <- complete_for(dl_f, build_formula("ftr_y", COVARS_M1, specs))
  for (grp in c("Acute", "Non-acute")) for (st in unique(c(min(yrs), CONFIG$trend_sens_start))) {
    log_msg("Group trend: FTR ", num, " / ", grp, " / from ", st)
    gt_rows[[length(gt_rows) + 1]] <- tryCatch(
      group_trend(dl_f, "ftr_y", lab, grp, st, max(yrs))[, outcome := paste0("ftr_", num)],
      error = function(e) data.table(outcome = paste0("ftr_", num), label = lab, group = grp, start_year = st,
                                     end_year = max(yrs), model = paste("error:", conditionMessage(e))))
  }
}
group_trends <- rbindlist(gt_rows, fill = TRUE)
write_out(group_trends, SUB, "year_trend_by_group_and_start")
# Publication layout: one row per outcome, four columns (group x start year)
gt_fmt <- copy(group_trends)[, cell := ifelse(is.na(or_per_year), model,
                                              sprintf("%s (%s-%s); p=%s", formatC(or_per_year, format = "f", digits = 3),
                                                      formatC(lower, format = "f", digits = 3),
                                                      formatC(upper, format = "f", digits = 3), fmt_p(p)))]
gt_fmt[, column := paste0(group, ", ", start_year, "-", end_year)]
gt_wide <- dcast(gt_fmt, outcome + label ~ column, value.var = "cell")
gt_wide <- gt_wide[order(match(outcome, c(OUT$var, "ftr_dzh", "ftr_d30")))]
write_out(gt_wide, SUB, "year_trend_by_group_and_start_table")

# ---- 4. Centre variation ----------------------------------------------------------------------
for (v in intersect(c("comp_primary", "d30", "dzh", "major_compl", "stroke_def", "renal", "vent24"), names(fits_keep))) {
  fk <- fits_keep[[v]]
  f1 <- fk$fits[[1]]$fit; d1 <- fk$data[[1]]
  # Profile-likelihood CI for the centre SD (complete-case fit; documented)
  cc <- cohort[year >= OUT[var == v, first_year] & year <= OUT[var == v, last_year] & !is.na(get(v))]
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
if ("comp_primary" %in% names(fits_keep)) {
  f1 <- fits_keep$comp_primary$fits[[1]]$fit; d1 <- fits_keep$comp_primary$data[[1]]
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
los_from <- exp_table(fit_pooled(f_los, lapply(los_dl, function(d) d[year >= CONFIG$trend_sens_start]), family = NULL)$pooled, "year_c")[
  , `:=`(outcome = paste0("LOS (discharged alive), from ", CONFIG$trend_sens_start), interpretation = "ratio of geometric means")]
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
write_out(rbind(los_tab, los_from, los_cox, fill = TRUE), SUB, "length_of_stay_models")
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
# Main analysis (paper 1 definitions): denominator = patients with a major complication
# (stroke, renal failure, GI complication); numerator = in-hospital death.
# Additional analyses: 30-day death as numerator; broader complication definition.
ftr_defs <- data.table(
  denominator = c("major_compl", "major_compl", "major_compl_new", "major_compl_new"),
  n_var       = c("n_compl", "n_compl", "n_compl_new", "n_compl_new"),
  numerator   = c("dzh", "d30", "dzh", "d30"),
  role        = c("main (paper 1 definitions)", "additional: 30-day death",
                  "additional: broader complication definition", "additional: broader definition, 30-day death"))
checks <- list(); decomps <- list(); ftr_models <- list()
for (k in seq_len(nrow(ftr_defs))) {
  den <- ftr_defs$denominator[k]; num <- ftr_defs$numerator[k]; nv <- ftr_defs$n_var[k]
  ftr_d <- cohort[get(den) == 1]
  chk_k <- data.table(
    definition = ftr_defs$role[k], denominator = den, numerator = num,
    item = c("patients with complication (denominator)", "denominator with outcome known",
             "deaths among denominator (numerator)", "all deaths in cohort",
             "deaths without recorded complication",
             "denominator equals sum(complication == 1)", "numerator <= all deaths"),
    value = c(nrow(ftr_d), sum(!is.na(ftr_d[[num]])), sum(ftr_d[[num]] == 1, na.rm = TRUE),
              sum(cohort[[num]] == 1, na.rm = TRUE),
              sum(cohort[[num]] == 1 & cohort[[den]] == 0, na.rm = TRUE),
              nrow(ftr_d) == sum(cohort[[den]] == 1, na.rm = TRUE),
              sum(ftr_d[[num]] == 1, na.rm = TRUE) <= sum(cohort[[num]] == 1, na.rm = TRUE)))
  stopifnot(chk_k$value[6] == 1, chk_k$value[7] == 1)
  checks[[k]] <- chk_k
  # Decomposition by year: mortality = P(C) * FTR + P(no C) * mortality without C
  dc <- cohort[!is.na(get(den)) & !is.na(get(num)),
               .(n = .N, p_complication = mean(get(den)),
                 ftr = mean(get(num)[get(den) == 1]), mort_no_compl = mean(get(num)[get(den) == 0]),
                 mortality = mean(get(num))), by = year][order(year)]
  dc[, check_identity := p_complication * ftr + (1 - p_complication) * mort_no_compl - mortality]
  stopifnot(all(abs(dc$check_identity) < 1e-10))
  decomps[[k]] <- cbind(definition = ftr_defs$role[k], denominator = den, numerator = num, dc)
  dl0 <- lapply(imps, function(d) d[get(den) == 1 & !is.na(get(num))])
  f <- build_formula(num, c(COVARS_M1, nv), specs)
  dl <- complete_for(dl0, f)
  ev <- epv(dl[[1]][[num]], n_fixed_par(f, dl[[1]]))
  if (ev < 10) {
    f <- build_formula(num, c(COVARS_REDUCED, nv), specs)
    dl <- complete_for(dl0, f)
  }
  fp <- fit_pooled(f, dl)
  ftr_models[[k]] <- exp_table(fp$pooled)[, `:=`(definition = ftr_defs$role[k], denominator = den,
                                                 numerator = num, n = nrow(dl[[1]]),
                                                 events = sum(dl[[1]][[num]] == 1), epv = ev,
                                                 reduced = ev < 10, centre_sd = fp$re_sd,
                                                 period = "all years")]
  dls <- lapply(dl, function(d) d[year >= CONFIG$trend_sens_start])
  if (sum(dls[[1]][[num]] == 1) >= 10) {
    fps <- fit_pooled(f, dls)
    ftr_models[[paste0(k, "_from")]] <- exp_table(fps$pooled, "year_c")[, `:=`(
      definition = ftr_defs$role[k], denominator = den, numerator = num, n = nrow(dls[[1]]),
      events = sum(dls[[1]][[num]] == 1), epv = NA_real_, reduced = ev < 10, centre_sd = fps$re_sd,
      period = paste("from", CONFIG$trend_sens_start))]
  }
}
write_out(rbindlist(checks), SUB, "ftr_numerator_denominator_check")
decomp_all <- rbindlist(decomps)
write_out(decomp_all, SUB, "ftr_decomposition_by_year")
decomp <- decomp_all[definition == ftr_defs$role[1]]
write_out(rbindlist(ftr_models), SUB, "ftr_models")

saveRDS(list(OUT = OUT, crude = crude, models = models, year = yr, year_from = yr_from, gcomp = gcomp,
             group_trends = group_trends, group_trends_wide = gt_wide,
             centre = rbindlist(centre_rows[!grepl("_cat|_fun", names(centre_rows))]),
             caterpillar = rbindlist(centre_rows[grepl("_cat", names(centre_rows))]),
             funnel = rbindlist(centre_rows[grepl("_fun", names(centre_rows))]),
             urgency = rbindlist(lapply(urg_rows, `[[`, "or")),
             los = rbind(los_tab, los_cox, fill = TRUE), ftr = rbindlist(ftr_models), ftr_decomp = decomp,
             ftr_decomp_all = decomp_all),
        derived_method("outcomes.rds"))
log_msg("Done.")

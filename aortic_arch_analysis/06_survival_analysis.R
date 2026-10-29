# =============================================================================
# 06_survival_analysis.R
# Kaplan-Meier, Cox regression (stratified by centre; pooled over imputations),
# proportional-hazards diagnostics, piecewise hazards, RMST, early vs late
# mortality, 30-day landmark analysis and conditional survival.
# Output: results/06_survival/*, data/derived/survival.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("06_survival_analysis")
log_msg("Missing-data method: ", method_label())
SUB <- "06_survival"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()   # imputed data sets or observed cohort (complete case)
specs  <- load_derived("specs.rds")

prep_surv <- function(d) {
  d <- d[!is.na(time_d) & !is.na(status)]
  d[, t_days := pmax(time_d, 0.5)]           # day-0 deaths/censorings -> 0.5 day
  d[, t_years := t_days / 365.25]
  d[]
}
sv <- prep_surv(cohort)
write_out(data.table(item = c("cohort", "with follow-up time and status", "excluded (missing time or status)",
                              "follow-up time 0 days (set to 0.5 day)", "deaths"),
                     n = c(nrow(cohort), nrow(sv), nrow(cohort) - nrow(sv),
                           sum(cohort$time_d == 0, na.rm = TRUE), sum(sv$status))),
          SUB, "survival_population")
imps_s <- lapply(imps, prep_surv)
TIMES <- CONFIG$km_times_days

# ---- 1. Follow-up --------------------------------------------------------------------------
rk <- reverse_km_median(sv$t_years, sv$status)
fu <- data.table(measure = c("reverse KM follow-up Q1 (years)", "reverse KM median follow-up (years)", "reverse KM follow-up Q3 (years)",
                             "maximum follow-up (years)", "total person-years"),
                 value = c(rk, max(sv$t_years), sum(sv$t_years)))
write_out(fu, SUB, "follow_up_summary")
# Follow-up of patients alive at last contact, by operation year (detects incomplete follow-up)
write_out(sv[status == 0, .(n_alive = .N, median_fu_years = median(t_years), min_fu_years = min(t_years),
                            max_fu_years = max(t_years)), by = year][order(year)], SUB, "follow_up_alive_by_year")

# ---- 2. Kaplan-Meier --------------------------------------------------------------------------
km_fits <- list(
  overall = survfit(Surv(t_days, status) ~ 1, data = sv),
  urgency = survfit(Surv(t_days, status) ~ urgency2, data = sv),
  pathology = survfit(Surv(t_days, status) ~ path3, data = sv),
  era = survfit(Surv(t_days, status) ~ era, data = sv))
km_tab <- rbindlist(lapply(names(km_fits), function(n) cbind(analysis = n, km_at(km_fits[[n]], TIMES))))
km_tab[, time_label := c("30 days", "1 year", "3 years", "5 years", "10 years")[match(round(time, 2), round(TIMES, 2))]]
write_out(km_tab, SUB, "km_survival_estimates")
lr <- rbindlist(lapply(c("urgency2", "path3", "era"), function(v) {
  s <- survdiff(as.formula(paste("Surv(t_days, status) ~", v)), data = sv)
  data.table(variable = v, chisq = s$chisq, df = length(s$n) - 1, p = pchisq(s$chisq, length(s$n) - 1, lower.tail = FALSE))
}))
write_out(lr, SUB, "logrank_tests")
# 1-year mortality by operation year (KM; handles censoring)
km1 <- rbindlist(lapply(sort(unique(sv$year)), function(y) {
  s <- summary(survfit(Surv(t_days, status) ~ 1, data = sv[year == y]), times = 365.25, extend = TRUE)
  data.table(year = y, n = nrow(sv[year == y]), mortality_1y = 1 - s$surv, lower = 1 - s$upper, upper = 1 - s$lower,
             n_risk_1y = s$n.risk)
}))
write_out(km1, SUB, "one_year_mortality_by_year_km")

# Plot data: curves on a grid with risk tables
curve_data <- function(sf, max_years = 12) {
  s <- summary(sf, times = seq(0, max_years * 365.25, by = 7), extend = TRUE)
  data.table(strata = if (is.null(s$strata)) "Overall" else sub(".*=", "", as.character(s$strata)),
             years = s$time / 365.25, surv = s$surv, lower = s$lower, upper = s$upper, n_risk = s$n.risk)
}
risk_table <- function(sf, max_years = 10) {
  s <- summary(sf, times = (0:max_years) * 365.25, extend = TRUE)
  data.table(strata = if (is.null(s$strata)) "Overall" else sub(".*=", "", as.character(s$strata)),
             years = round(s$time / 365.25), n_risk = s$n.risk)
}
curves <- lapply(km_fits, curve_data); risks <- lapply(km_fits, risk_table)

# ---- 3. Early vs late mortality ------------------------------------------------------------------
el <- sv[, .(deaths_0_30d = sum(status == 1 & t_days <= 30),
             deaths_31_365d = sum(status == 1 & t_days > 30 & t_days <= 365.25),
             deaths_gt_1y = sum(status == 1 & t_days > 365.25), n = .N), by = urgency2]
el <- rbind(el, sv[, .(urgency2 = "All", deaths_0_30d = sum(status == 1 & t_days <= 30),
                       deaths_31_365d = sum(status == 1 & t_days > 30 & t_days <= 365.25),
                       deaths_gt_1y = sum(status == 1 & t_days > 365.25), n = .N)])
write_out(el, SUB, "early_late_deaths")

# ---- 4. Cox regression (pooled over imputations) ------------------------------------------------
cox_rhs <- function(covars, extra = character(0)) paste(c("year_c", covar_terms(covars, specs), extra), collapse = " + ")
pool_cox <- function(formula, dl, label) {
  fits <- par_lapply(dl, function(d) coxph(formula, data = d, x = FALSE, model = FALSE))
  est <- lapply(fits, coef); vc <- lapply(fits, vcov)
  exp_table(rubin_pool(est, vc))[, `:=`(model = label, n = fits[[1]]$n, events = fits[[1]]$nevent)]
}
f_main <- as.formula(paste("Surv(t_years, status) ~", cox_rhs(COVARS_M1), "+ strata(centrum)"))
cox_main <- pool_cox(f_main, imps_s, "Cox stratified by centre (primary)")
covars_ne <- c(setdiff(COVARS_M1, c("urg3", "path3")), "nonelective")
cox_ne_a <- pool_cox(as.formula(paste("Surv(t_years, status) ~", cox_rhs(covars_ne), "+ strata(centrum)")),
                     imps_s, "non-elective vs elective, without pathology")
cox_ne_b <- pool_cox(as.formula(paste("Surv(t_years, status) ~", cox_rhs(c(covars_ne, "path3")), "+ strata(centrum)")),
                     imps_s, "non-elective vs elective, with pathology")
# Sensitivity: gamma frailty; robust cluster SE
cox_frail <- pool_cox(as.formula(paste("Surv(t_years, status) ~", cox_rhs(COVARS_M1), "+ frailty(centrum)")),
                      imps_s, "Cox with gamma frailty for centre")
cox_rob <- pool_cox(as.formula(paste("Surv(t_years, status) ~", cox_rhs(COVARS_M1), "+ cluster(centrum)")),
                    imps_s, "Cox, robust SE clustered by centre")
cox_all <- rbind(cox_main, cox_ne_a, cox_ne_b, cox_frail, cox_rob, fill = TRUE)
write_out(cox_all, SUB, "cox_models")

# Proportional hazards diagnostics (imputation 1)
fit1 <- coxph(f_main, data = imps_s[[1]])
zph <- cox.zph(fit1, terms = TRUE)
write_out(data.table(term = rownames(zph$table), zph$table), SUB, "cox_ph_test_schoenfeld")
grDevices::png(file.path(out_dir(SUB), "cox_schoenfeld_residuals.png"), 2400, 1600, res = 200)
nt <- nrow(zph$table) - 1
par(mfrow = c(ceiling(nt / 3), 3), mar = c(4, 4, 2, 1))
for (i in seq_len(nt)) plot(zph[i], resid = FALSE, main = rownames(zph$table)[i])
grDevices::dev.off()
# Functional form of age: martingale residuals from a model without age
f_noage <- as.formula(paste("Surv(t_years, status) ~", cox_rhs(setdiff(COVARS_M1, "age")), "+ strata(centrum)"))
d_mart <- complete_for(imps_s[1], f_noage)[[1]][!is.na(age)]
mres <- residuals(coxph(f_noage, data = d_mart), type = "martingale")
p_mart <- ggplot(data.table(age = d_mart$age, r = mres), aes(age, r)) +
  geom_point(alpha = .1, size = .5) + geom_smooth(method = "loess", se = TRUE, colour = PAL[2]) +
  labs(x = "Age (years)", y = "Martingale residual", title = "Functional form of age (model without age)") + theme_pub()
save_fig(p_mart, SUB, "cox_martingale_age")

# ---- 5. Piecewise hazards (non-proportional hazards) ---------------------------------------------
cuts <- CONFIG$piecewise_cuts / 365.25
pw <- lapply(imps_s, function(d) {
  s <- survSplit(Surv(t_years, status) ~ ., data = as.data.frame(d), cut = cuts, episode = "interval")
  s$interval <- factor(s$interval, labels = c("0-30 days", "30 days-1 year", ">1 year"))
  as.data.table(s)
})
f_pw <- as.formula(paste("Surv(tstart, t_years, status) ~",
                         paste(c("nonelective:interval", "acute_dissection:interval", "year_c:interval",
                                 covar_terms(c("age", "female", "egfr", "lvef_cat", "prior_cs", "endocarditis_b"), specs)),
                               collapse = " + "), "+ strata(centrum)"))
cox_pw <- pool_cox(f_pw, pw, "piecewise Cox (0-30 d, 30 d-1 y, >1 y)")
write_out(cox_pw[grepl("interval", term)], SUB, "cox_piecewise_hazard_ratios")

# ---- 6. RMST -----------------------------------------------------------------------------------
rmst <- rbindlist(lapply(CONFIG$rmst_tau_years, function(tau) {
  sf <- survfit(Surv(t_years, status) ~ urgency2, data = sv)
  max_fu <- sv[, .(m = max(t_years)), by = urgency2]
  if (any(max_fu$m < tau)) return(data.table(tau_years = tau, note = "tau exceeds follow-up in a group"))
  r <- rmst_by(sf, tau)
  dif <- r$rmst[2] - r$rmst[1]; se <- sqrt(sum(r$se^2))
  data.table(tau_years = tau, rmst_elective = r$rmst[1], rmst_nonelective = r$rmst[2],
             difference = dif, lower = dif - 1.96 * se, upper = dif + 1.96 * se,
             p = 2 * pnorm(-abs(dif / se)))
}), fill = TRUE)
write_out(rmst, SUB, "rmst_by_urgency")

# ---- 7. 30-day landmark ----------------------------------------------------------------------
LM <- CONFIG$landmark_days
lm_prep <- function(d) d[t_days > LM][, `:=`(t_lm = (t_days - LM) / 365.25)]
lm_d <- lm_prep(sv)
write_out(data.table(item = c("alive and in follow-up beyond day 30", "deaths after day 30", "excluded (died or censored <= 30 d)"),
                     n = c(nrow(lm_d), sum(lm_d$status), nrow(sv) - nrow(lm_d))), SUB, "landmark_population")
lm_fits <- list(urgency = survfit(Surv(t_lm * 365.25, status) ~ urgency2, data = lm_d),
                pathology = survfit(Surv(t_lm * 365.25, status) ~ path3, data = lm_d),
                era = survfit(Surv(t_lm * 365.25, status) ~ era, data = lm_d))
lm_tab <- rbindlist(lapply(names(lm_fits), function(n) cbind(analysis = n, km_at(lm_fits[[n]], c(365.25, 5 * 365.25, 10 * 365.25)))))
write_out(lm_tab, SUB, "landmark30_km_estimates")
imps_lm <- lapply(imps_s, lm_prep)
cox_lm <- rbind(
  pool_cox(as.formula(paste("Surv(t_lm, status) ~", cox_rhs(COVARS_M1), "+ strata(centrum)")), imps_lm, "landmark 30 d: primary covariates"),
  pool_cox(as.formula(paste("Surv(t_lm, status) ~", cox_rhs(covars_ne), "+ strata(centrum)")), imps_lm, "landmark 30 d: non-elective vs elective"))
write_out(cox_lm, SUB, "landmark30_cox")
zph_lm <- cox.zph(coxph(as.formula(paste("Surv(t_lm, status) ~", cox_rhs(COVARS_M1), "+ strata(centrum)")), data = imps_lm[[1]]), terms = TRUE)
write_out(data.table(term = rownames(zph_lm$table), zph_lm$table), SUB, "landmark30_ph_test")
lm_curves <- lapply(lm_fits, curve_data)

# ---- 8. Conditional survival ---------------------------------------------------------------
s_vals <- c(0, 30, 365.25, 2 * 365.25, 5 * 365.25)
cond <- rbindlist(lapply(s_vals, function(s) {
  rbindlist(lapply(c("All", levels(sv$urgency2)), function(g) {
    d <- sv[t_days > s & (g == "All" | urgency2 == g)]
    if (nrow(d) < 20) return(NULL)
    sf <- survfit(Surv(t_days - s, status) ~ 1, data = d)
    r <- km_at(sf, c(365.25, 5 * 365.25))
    data.table(survived_days = s, group = g, n_at_s = nrow(d), additional_years = c(1, 5),
               cond_surv = r$surv, lower = r$lower, upper = r$upper, n_risk = r$n_risk)
  }))
}))
write_out(cond, SUB, "conditional_survival")

saveRDS(list(curves = curves, risks = risks, km = km_tab, km1 = km1, cox = cox_all, pw = cox_pw,
             lm_curves = lm_curves, lm_risks = lapply(lm_fits, risk_table), cox_lm = cox_lm,
             cond = cond, rmst = rmst, fu = fu), derived_method("survival.rds"))
log_msg("Done.")

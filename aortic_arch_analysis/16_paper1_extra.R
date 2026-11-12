# =============================================================================
# 16_paper1_extra.R
# Additional analyses for paper 1 that were first run interactively, now part of
# the pipeline, plus every number that 17_figures_paper1.R needs. Same models,
# covariates and helper functions as 04/06 (complete-case or imputation via
# analysis_datasets()).
#   A. Year x acute interaction, 12 outcomes (Benjamini-Hochberg)
#   B. Year trends overall and per group (from 2013 and from 2015), 12 outcomes + FTR
#   C. Acute vs non-acute adjusted odds ratios
#   D. Case-mix-standardised risk by year, overall and per group (g-computation)
#   E. Crude rates by outcome, group and period
#   F. Perfusion practice by year and group; hypothermia categories
#   G. Failure-to-rescue decomposition by year and group; deaths without major complication
#   H. Patients with complete vs incomplete covariates; % complete by year and group
#   I. Practice profile per centre (no outcomes)
#   J. Survival: piecewise and landmark Cox, log-rank tests, KM grid
# Output: results/<method>/16_paper1_extra/*.csv
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival) })
SUB <- "16_paper1_extra"
log_msg("Start 16_paper1_extra")

cohort <- load_derived("cohort.rds")
specs  <- load_derived("specs.rds")
cc     <- analysis_datasets()
D1 <- "d1y"; D1_END <- CONFIG$d1y_last_year
OUTS <- c(comp_primary = "Primary composite outcome", dzh = "In-hospital death", d30 = "30-day death",
          d1y = "One-year death", stroke_def = "Stroke with permanent deficit",
          stroke_nodef = "Stroke or TIA without deficit", renal = "Renal failure", gi = "GI complication",
          vent24 = "Ventilation >24 h", reintub = "Reintubation", icu_readm = "ICU readmission",
          major_compl = "Major complication")
OUTS <- OUTS[names(OUTS) %in% names(cohort)]
SECONDARY <- setdiff(names(OUTS), "comp_primary")
end_of <- function(v) if (v == D1) D1_END else max(CONFIG$study_years)
`%||%` <- function(a, b) if (is.null(a)) b else a
in_group <- function(d, g) if (is.null(g)) d else d[acute == g]
safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })
na_row <- data.table(ratio = NA_real_, ratio_lower = NA_real_, ratio_upper = NA_real_, p = NA_real_)

# Model data for outcome v, optional group (1 acute, 0 non-acute) and start year, with
# the reduced covariate set when events per variable < 10 (as in 04).
model_data <- function(v, g = NULL, start = min(CONFIG$study_years), extra = character(0)) {
  sub <- function(d) in_group(d[year >= start & year <= end_of(v) & !is.na(get(v))], g)
  f <- build_formula(v, c(COVARS_M1, extra), specs); dl <- complete_for(lapply(cc, sub), f); model <- "M1"
  y <- dl[[1]][[v]]
  if (sum(y == 1) > 0 && epv(y, n_fixed_par(f, dl[[1]])) < 10) {
    f <- build_formula(v, c(COVARS_REDUCED, extra), specs); dl <- complete_for(lapply(cc, sub), f); model <- "reduced"
  }
  list(f = f, dl = dl, model = model, n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1))
}
year_or <- function(v, g = NULL, start = min(CONFIG$study_years)) {
  m <- model_data(v, g, start)
  r <- if (m$events < 10) na_row else
    safe(exp_table(fit_pooled(m$f, m$dl)$pooled, "year_c"), paste(v, g, start)) %||% na_row
  data.table(outcome = v, label = OUTS[[v]], group = if (is.null(g)) "All" else if (g == 1) "Acute" else "Non-acute",
             from = start, to = end_of(v), model = m$model, n = m$n, events = m$events,
             or_per_year = r$ratio[1], lower = r$ratio_lower[1], upper = r$ratio_upper[1], p = r$p[1])
}
# ---- A. Year x acute interaction -------------------------------------------------------
log_msg("A. interaction")
inter <- rbindlist(lapply(names(OUTS), function(v) {
  m <- model_data(v, extra = c("acute", "year_c:acute"))
  p <- safe(fit_pooled(m$f, m$dl)$pooled, v)
  r <- if (is.null(p)) NULL else p[grepl(":", term) & grepl("acute", term) & grepl("year_c", term)]
  if (is.null(r) || !nrow(r)) return(data.table(outcome = v, label = OUTS[[v]], n = m$n, events = m$events))
  data.table(outcome = v, label = OUTS[[v]], n = m$n, events = m$events, ratio_of_or = exp(r$estimate[1]),
             lower = exp(r$lower[1]), upper = exp(r$upper[1]), p = r$p[1])
}), fill = TRUE)
inter[, p_bh := bh(p)]
write_out(inter, SUB, "A_year_by_acute_interaction")

# ---- B. Year trends overall and per group ----------------------------------------------
log_msg("B. year trends")
trends <- rbindlist(lapply(names(OUTS), function(v) rbindlist(c(
  list(year_or(v)),
  lapply(c(1L, 0L), function(g) rbindlist(lapply(unique(c(min(CONFIG$study_years), CONFIG$trend_sens_start)),
                                                 function(st) year_or(v, g, st))))))))
trends[group == "All" & outcome %in% SECONDARY, p_holm := holm(p)]
trends[group != "All", p_bh := bh(p), by = .(from)]
write_out(trends, SUB, "B_year_trends_overall_and_by_group")
ftr_dl <- lapply(cc, function(d) d[major_compl == 1 & !is.na(dzh)][, ftr := dzh])
ftr_tr <- rbindlist(lapply(list(NULL, 1L, 0L), function(g) {
  dd <- lapply(ftr_dl, function(d) if (is.null(g)) d else d[acute == g])
  f <- build_formula("ftr", c(COVARS_M1, if ("n_compl" %in% names(dd[[1]])) "n_compl"), specs)
  dl <- complete_for(dd, f)
  r <- safe(exp_table(fit_pooled(f, dl)$pooled, "year_c"), "FTR") %||% na_row
  data.table(group = if (is.null(g)) "All" else if (g == 1) "Acute" else "Non-acute", n = nrow(dl[[1]]),
             deaths = sum(dl[[1]]$ftr == 1), or_per_year = r$ratio[1], lower = r$ratio_lower[1], upper = r$ratio_upper[1], p = r$p[1])
}))
write_out(ftr_tr, SUB, "B_failure_to_rescue_trend")

# ---- C. Acute vs non-acute -------------------------------------------------------------
log_msg("C. acute vs non-acute")
avn <- rbindlist(lapply(names(OUTS), function(v) {
  m <- model_data(v, extra = "acute")
  r <- safe(exp_table(fit_pooled(m$f, m$dl)$pooled, "acute"), v) %||% na_row
  data.table(outcome = v, label = OUTS[[v]], n = m$n, events = m$events, or_acute = r$ratio[1],
             lower = r$ratio_lower[1], upper = r$ratio_upper[1], p = r$p[1])
}))
write_out(avn, SUB, "C_acute_vs_nonacute_or")

# ---- D. Standardised risk by year -------------------------------------------------------
log_msg("D. g-computation")
std <- rbindlist(lapply(intersect(c("dzh", "d30", "comp_primary"), names(OUTS)), function(v) rbindlist(lapply(list(NULL, 1L, 0L), function(g) {
  sub <- function(d) in_group(d[!is.na(get(v))], g)
  fs <- build_formula(v, COVARS_M1, specs, year_term = YEAR_SPLINE)
  dl <- complete_for(lapply(cc, sub), fs)
  ms <- safe(fit_pooled(fs, dl), paste("gcomp", v))
  if (is.null(ms)) return(NULL)
  obs <- dl[[1]][, .(observed = mean(get(v) == 1), n = .N), by = year]
  out <- gcomp_mi(lapply(ms$fits, `[[`, "fit"), dl, sort(unique(dl[[1]]$year)))
  merge(cbind(outcome = v, group = if (is.null(g)) "All" else if (g == 1) "Acute" else "Non-acute", out), obs, by = "year")
}))))
write_out(std, SUB, "D_standardised_risk_by_year")

# ---- E. Crude rates by group and period ------------------------------------------------
cohort[, period := cut(year, c(2012, 2016, 2020, 2024), c("2013-2016", "2017-2020", "2021-2024"))]
crude <- rbindlist(lapply(names(OUTS), function(v) {
  d <- cohort[!is.na(get(v)) & year <= end_of(v)]
  rbind(d[, .(group = "All", period = "All", n = .N, events = sum(get(v) == 1))],
        d[, .(period = "All", n = .N, events = sum(get(v) == 1)), by = .(group = as.character(urgency2))],
        d[, .(group = "All", n = .N, events = sum(get(v) == 1)), by = .(period = as.character(period))],
        d[, .(n = .N, events = sum(get(v) == 1)), by = .(group = as.character(urgency2), period = as.character(period))],
        use.names = TRUE)[, `:=`(outcome = v, label = OUTS[[v]], pct = 100 * events / n)]
}))
write_out(crude, SUB, "E_crude_rates")

# ---- F. Practice by year and group ------------------------------------------------------
pr_base <- cohort[ecc_any == 1 | is.na(ecc_any)]
practice <- function(d) d[, .(n = .N,
  deep = 100 * mean(hypo4 == "Deep (<=20)", na.rm = TRUE), moderate = 100 * mean(hypo4 == "Moderate (20.1-28)", na.rm = TRUE),
  mild = 100 * mean(hypo4 == "Mild (28.1-34)", na.rm = TRUE), normo = 100 * mean(hypo4 == "Normothermia (>34)", na.rm = TRUE),
  temp_missing_pct = 100 * mean(is.na(temp_min)),
  temp_median = as.numeric(median(temp_min, na.rm = TRUE)), temp_q1 = as.numeric(quantile(temp_min, .25, na.rm = TRUE)),
  temp_q3 = as.numeric(quantile(temp_min, .75, na.rm = TRUE)), arrest_pct = 100 * mean(circ_arrest_b == 1, na.rm = TRUE),
  ca_median = as.numeric(median(ca_min[circ_arrest_b == 1], na.rm = TRUE))), by = year]
prac <- rbind(cbind(group = "All", practice(pr_base)),
              cbind(group = "Acute", practice(pr_base[urgency2 == "Acute"])),
              cbind(group = "Non-acute", practice(pr_base[urgency2 == "Non-acute"])))[order(group, year)]
write_out(prac, SUB, "F_practice_by_year_and_group")

# ---- G. Failure-to-rescue decomposition --------------------------------------------------
fd <- cohort[!is.na(major_compl) & !is.na(dzh)]
decomp <- function(d) d[, .(n = .N, compl_pct = 100 * mean(major_compl == 1), ftr_pct = 100 * mean(dzh[major_compl == 1] == 1),
                            mort_no_compl_pct = 100 * mean(dzh[major_compl == 0] == 1), mortality_pct = 100 * mean(dzh == 1),
                            n_compl = sum(major_compl == 1), deaths_compl = sum(dzh == 1 & major_compl == 1),
                            deaths_no_compl = sum(dzh == 1 & major_compl == 0))]
ftr_year <- rbind(cbind(group = "All", fd[, decomp(.SD), by = year]),
                  cbind(group = "Acute", fd[urgency2 == "Acute", decomp(.SD), by = year]),
                  cbind(group = "Non-acute", fd[urgency2 == "Non-acute", decomp(.SD), by = year]))[order(group, year)]
write_out(ftr_year, SUB, "G_ftr_decomposition_by_year_and_group")
ftr_grp <- rbind(cbind(group = "All", decomp(fd)), fd[, decomp(.SD), by = .(group = as.character(urgency2))], fill = TRUE)
ftr_grp[, pct_deaths_without_major_compl := 100 * deaths_no_compl / (deaths_compl + deaths_no_compl)]
write_out(ftr_grp, SUB, "G_ftr_by_group")

# ---- H. Complete vs incomplete covariate data -----------------------------------------
cohort[, complete_covariates := complete.cases(cohort[, COVARS_M1, with = FALSE])]
cmp <- cohort[, .(n = .N, age_median = as.numeric(median(age, na.rm = TRUE)), female_pct = 100 * mean(female == 1, na.rm = TRUE),
                  acute_pct = 100 * mean(acute == 1), redo_pct = 100 * mean(prior_cs == 1, na.rm = TRUE),
                  dzh_pct = 100 * mean(dzh == 1, na.rm = TRUE), comp_pct = 100 * mean(comp_primary == 1, na.rm = TRUE)),
              by = complete_covariates]
cmp[, p_dzh := chisq.test(table(cohort$complete_covariates, cohort$dzh))$p.value]
write_out(cmp, SUB, "H_complete_vs_incomplete")
write_out(cohort[, .(n = .N, pct_complete = 100 * mean(complete_covariates)), by = .(year, group = as.character(urgency2))][order(group, year)],
          SUB, "H_pct_complete_by_year_and_group")
write_out(rbindlist(lapply(COVARS_M1, function(v) cohort[, .(variable = v, missing = sum(is.na(get(v))), pct = 100 * mean(is.na(get(v))))])),
          SUB, "H_missing_model_covariates")

# ---- I. Centre practice profile ------------------------------------------------------------
prof <- cohort[, .(n = .N, pct_acute = 100 * mean(urgency2 == "Acute", na.rm = TRUE),
                   temp_median = as.numeric(median(temp_min, na.rm = TRUE)), temp_missing_pct = 100 * mean(is.na(temp_min)),
                   pct_deep = 100 * mean(temp_min <= 20, na.rm = TRUE), pct_arrest = 100 * mean(circ_arrest_b == 1, na.rm = TRUE),
                   ca_median = as.numeric(median(ca_min[circ_arrest_b == 1 & ca_min > 0], na.rm = TRUE))),
               by = centre_lab][order(centre_lab)]
write_out(prof, SUB, "I_centre_practice_profile")

# ---- J. Survival -------------------------------------------------------------------------
log_msg("J. survival")
LM <- CONFIG$landmark_days
prep <- function(d) d[!is.na(time_d) & !is.na(status)][, `:=`(t_days = pmax(time_d, 0.5), t_years = pmax(time_d, 0.5) / 365.25)]
sv <- prep(copy(cohort)); svl <- lapply(cc, function(d) prep(copy(d)))
pool_cox <- function(formula, dl) {
  fits <- lapply(dl, function(d) coxph(formula, data = d))
  exp_table(rubin_pool(lapply(fits, coef), lapply(fits, vcov)))
}
rhs <- paste(covar_terms(COVARS_M1, specs), collapse = " + ")
pw <- lapply(svl, function(d) {
  s <- survSplit(Surv(t_years, status) ~ ., data = as.data.frame(d), cut = CONFIG$piecewise_cuts / 365.25, episode = "interval")
  s$interval <- factor(s$interval, labels = c("0-30 days", "30 days-1 year", ">1 year")); as.data.table(s) })
f_pw <- as.formula(paste("Surv(tstart, t_years, status) ~ acute:interval + year_c:interval +", rhs, "+ strata(centrum)"))
hr_pw <- pool_cox(f_pw, complete_for(pw, f_pw))[grepl("acute", term) & grepl("interval", term)]
hr_pw[, interval := sub(".*interval", "", term)]
lm_l <- lapply(svl, function(d) d[t_days > LM][, t_lm := (t_days - LM) / 365.25])
f_lm <- as.formula(paste("Surv(t_lm, status) ~ year_c + acute +", rhs, "+ strata(centrum)"))
hr_lm <- pool_cox(f_lm, complete_for(lm_l, f_lm))[term == "acute"][, interval := "landmark (30-day survivors)"]
lr_a <- survdiff(Surv(pmin(t_days, LM), as.integer(status == 1 & t_days <= LM)) ~ urgency2, data = sv)
lr_b <- survdiff(Surv((t_days - LM) / 365.25, status) ~ urgency2, data = sv[t_days > LM])
hr <- rbind(hr_pw, hr_lm, fill = TRUE)[, .(interval, hr = ratio, lower = ratio_lower, upper = ratio_upper, p)]
write_out(hr, SUB, "J_hazard_ratios_acute")
write_out(data.table(test = c("days 0-30", "after day 30 (landmark)"),
                     p = c(pchisq(lr_a$chisq, 1, lower.tail = FALSE), pchisq(lr_b$chisq, 1, lower.tail = FALSE))),
          SUB, "J_logrank")
km_grid <- function(fit, times, scale, part) {
  s <- summary(fit, times = times * scale, extend = TRUE)
  out <- data.table(part = part, time = s$time / scale, group = sub("^urgency2=", "", as.character(s$strata)),
                    n_risk = s$n.risk, surv = s$surv, lower = s$lower, upper = s$upper)
  out[n_risk < KM_MIN_RISK, c("surv", "lower", "upper") := NA_real_]
  out }
fa <- survfit(Surv(pmin(t_days, LM), as.integer(status == 1 & t_days <= LM)) ~ urgency2, data = sv, conf.type = "log-log")
fb <- survfit(Surv((t_days - LM) / 365.25, status) ~ urgency2, data = sv[t_days > LM], conf.type = "log-log")
write_out(rbind(km_grid(fa, seq(0, LM, 2), 1, "days 0-30"), km_grid(fb, seq(0, 10, .25), 1, "years after day 30")),
          SUB, "J_km_grid_landmark")
log_msg("Done.")

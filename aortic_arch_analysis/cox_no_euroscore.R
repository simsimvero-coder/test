# =============================================================================
# cox_no_euroscore.R  (standalone; run after run_all.R)
# Piecewise Cox (0-30 days, 30 days-1 year, >1 year) and 30-day landmark Cox for
# acute versus non-acute presentation, with and without EuroSCORE II:
#   1. with EuroSCORE II (primary model, as in Table 5)
#   2. without EuroSCORE II, in the same patients as model 1
#   3. without EuroSCORE II, in all patients with complete age, sex and previous
#      cardiac surgery
# All models: stratified by centre, adjusted for age, sex, previous cardiac surgery
# and calendar year (as in 16_paper1_extra.R, section J).
# Writes only aggregate numbers. Output: results/<method>/cox_no_euroscore/
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival) })
SUB <- "cox_no_euroscore"
options(width = 200)
log_msg("Start cox_no_euroscore")

specs <- load_derived("specs.rds")
cc    <- analysis_datasets()
LM    <- CONFIG$landmark_days
prep  <- function(d) d[!is.na(time_d) & !is.na(status)][, `:=`(t_days = pmax(time_d, 0.5), t_years = pmax(time_d, 0.5) / 365.25)]
svl   <- lapply(cc, function(d) prep(copy(d)))
pool_cox <- function(formula, dl) {
  fits <- lapply(dl, function(d) coxph(formula, data = d))
  r <- exp_table(rubin_pool(lapply(fits, coef), lapply(fits, vcov)))
  r[, `:=`(n = fits[[1]]$n, events = fits[[1]]$nevent)][]
}
rhs_with <- paste(covar_terms(COVARS_M1, specs), collapse = " + ")
rhs_no   <- paste(covar_terms(c("age", "female", "prior_cs"), specs), collapse = " + ")
# Patients with complete data for the primary model (incl. EuroSCORE II)
vars_m1 <- unique(c("year_c", "acute", "centrum", all.vars(as.formula(paste("~", rhs_with)))))
same <- function(d) d[complete.cases(d[, .SD, .SDcols = intersect(vars_m1, names(d))])]

split_pw <- function(d) {
  s <- survSplit(Surv(t_years, status) ~ ., data = as.data.frame(d), cut = CONFIG$piecewise_cuts / 365.25, episode = "interval")
  s$interval <- factor(s$interval, labels = c("0-30 days", "30 days-1 year", ">1 year")); as.data.table(s)
}
run <- function(rhs, dl, label) {
  f_pw <- as.formula(paste("Surv(tstart, t_years, status) ~ acute:interval + year_c:interval +", rhs, "+ strata(centrum)"))
  pw <- complete_for(lapply(dl, split_pw), f_pw)
  a <- pool_cox(f_pw, pw)[grepl("acute", term) & grepl("interval", term)]
  a[, interval := sub(".*interval", "", term)]
  f_lm <- as.formula(paste("Surv(t_lm, status) ~ year_c + acute +", rhs, "+ strata(centrum)"))
  lm_l <- complete_for(lapply(dl, function(d) d[t_days > LM][, t_lm := (t_days - LM) / 365.25]), f_lm)
  b <- pool_cox(f_lm, lm_l)[term == "acute"][, interval := "landmark (30-day survivors)"]
  n_pat <- nrow(complete_for(dl, as.formula(paste("~", rhs, "+ year_c + acute")))[[1]])
  cbind(model = label, patients = n_pat,
        rbind(a, b, fill = TRUE)[, .(interval, n_rows_or_patients = n, events, hr = ratio, lower = ratio_lower, upper = ratio_upper, p)])
}
res <- rbind(
  run(rhs_with, lapply(svl, same), "1. With EuroSCORE II (primary)"),
  run(rhs_no,   lapply(svl, same), "2. Without EuroSCORE II, same patients"),
  run(rhs_no,   svl,               "3. Without EuroSCORE II, all patients"))
write_out(res, SUB, "cox_acute_with_without_euroscore")
out <- c(sprintf("Cox acute vs non-acute, with and without EuroSCORE II - %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
         capture.output(print(res[, .(model, patients, interval, events, HR = sprintf("%.3f (%.3f-%.3f)", hr, lower, upper), p = signif(p, 3))], nrows = 50)))
f_out <- file.path(out_dir(SUB), "overzicht_cox.txt")
writeLines(out, f_out); cat(out, sep = "\n")
log_msg("Done. Overview in ", f_out)

# =============================================================================
# sens_strict_groups.R  (standalone; run after run_all.R)
# Paper 1, sensitivity analysis with a stricter group definition:
#   Acute    = urgentie 30 (emergency) or 40 (salvage), or acute dissection (code 20)
#              in any aortic segment (unchanged from the primary definition);
#   Elective = urgentie 10 (elective) without a dissection code (20) in any segment.
#              (Chronic dissection, code 30, and infection, code 40, are already
#              excluded from the cohort.)
#   Urgent operations (urgentie 20) and non-acute operations with unknown urgency
#   belong to neither group and are left out; their number is reported (A).
# Same models as the primary analysis:
#   A. group counts (primary versus strict definition), overall and per year;
#   B. crude outcomes per group;
#   C. acute versus elective: adjusted OR (M1, or the reduced set when events per
#      variable < 10), with and without EuroSCORE II in the same sample;
#   D. year trend per group, 2013-2024 and from 2018 (Holm over the 11 secondary
#      outcomes, per group and period);
#   E. piecewise Cox (0-30 days, 30 days-1 year, >1 year) and 30-day landmark Cox,
#      acute versus elective, with and without EuroSCORE II (same patients).
# Writes only aggregate numbers. Output: results/<method>/sens_strict_groups/
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival) })
SUB <- "sens_strict_groups"
options(width = 220)
log_msg("Start sens_strict_groups")

specs <- load_derived("specs.rds")
strict <- function(d) {
  d <- copy(d)
  d[, group_primary := as.character(urgency2)]
  d[, acute := fcase(acute %in% 1L, 1L, urgentie %in% 10 & !(seg_ad %in% TRUE), 0L, default = NA_integer_)]
  d[, urgency2 := factor(fcase(acute == 1L, "Acute", acute == 0L, "Elective"), c("Elective", "Acute"))]
  d
}
cohort_all <- strict(load_derived("cohort.rds"))
cc_all     <- lapply(analysis_datasets(), strict)
cohort     <- cohort_all[!is.na(acute)]
cc         <- lapply(cc_all, function(d) d[!is.na(acute)])
END    <- max(CONFIG$study_years); D1_END <- CONFIG$d1y_last_year
OUTS <- c(comp_primary = "Primary composite outcome", dzh = "In-hospital death", d30 = "30-day death",
          d1y = "One-year death", stroke_def = "Stroke with permanent deficit",
          stroke_nodef = "Stroke or TIA without deficit", renal = "Renal failure", gi = "GI complication",
          vent24 = "Ventilation >24 h", reintub = "Reintubation", icu_readm = "ICU readmission",
          major_compl = "Major complication")
OUTS <- OUTS[names(OUTS) %in% names(cohort)]
end_of <- function(v) if (v == "d1y") D1_END else END
safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })
row_of <- function(fit, term) {
  if (is.null(fit)) return(data.table(or = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  tt <- term; p <- fit$pooled[fit$pooled$term == tt]
  if (!nrow(p)) return(data.table(or = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  data.table(or = exp(p$estimate), lower = exp(p$lower), upper = exp(p$upper), p = p$p)
}
# M1, or the reduced set when events per variable < 10 (as in the primary analysis)
fit_rule <- function(v, covs, dsub, reduced_covs = COVARS_REDUCED) {
  f <- build_formula(v, covs, specs); model <- "M1"
  dl <- complete_for(dsub, f); y <- dl[[1]][[v]]
  if (sum(y == 1) > 0 && epv(y, n_fixed_par(f, dl[[1]])) < 10) {
    f <- build_formula(v, reduced_covs, specs); model <- "reduced"
    dl <- complete_for(dsub, f); y <- dl[[1]][[v]]
  }
  list(fit = if (sum(y == 1) >= 10) safe(fit_pooled(f, dl), v) else NULL, model = model, n = nrow(dl[[1]]), events = sum(y == 1))
}

# ---- A. Group counts ----------------------------------------------------------------------
A1 <- cohort_all[, .N, by = .(group_primary, group_strict = fifelse(is.na(acute), "Neither (urgent or unknown urgency)",
                                                                   as.character(urgency2)), urgentie)][order(group_primary, group_strict, urgentie)]
A2 <- dcast(cohort_all[, .N, by = .(year, group_strict = fifelse(is.na(acute), "Neither", as.character(urgency2)))],
            year ~ group_strict, value.var = "N", fill = 0)
write_out(A1, SUB, "A1_group_counts"); write_out(A2, SUB, "A2_group_counts_by_year")

# ---- B. Crude outcomes per group (incl. the excluded 'neither' group, for comparison) ----------
pct <- function(x) round(100 * mean(x == 1, na.rm = TRUE), 1)
B <- rbindlist(lapply(names(OUTS), function(v) cohort_all[!is.na(get(v)) & year <= end_of(v),
  .(outcome = v, n = .N, events = sum(get(v) == 1), pct = pct(get(v))),
  by = .(group = fifelse(is.na(acute), "Neither", as.character(urgency2)))]))
B0 <- cohort_all[, .(n = .N, age_median = median(age, na.rm = TRUE), female_pct = pct(female),
                     euroscore2_median = round(median(euroscore2, na.rm = TRUE), 1)),
                 by = .(group = fifelse(is.na(acute), "Neither", as.character(urgency2)))]
write_out(B, SUB, "B_crude_outcomes"); write_out(B0, SUB, "B0_case_mix")

# ---- C. Acute versus elective, with and without EuroSCORE II (same sample) ----------------------
C <- rbindlist(lapply(names(OUTS), function(v) {
  dsub <- lapply(cc, function(d) d[!is.na(get(v)) & year <= end_of(v)])
  f_m1 <- build_formula(v, c(COVARS_M1, "acute"), specs)
  f_no <- build_formula(v, c("age", "female", "prior_cs", "acute"), specs)
  dl <- complete_for(dsub, f_m1)
  rbind(cbind(outcome = v, label = OUTS[[v]], model = "With EuroSCORE II", n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1),
              row_of(safe(fit_pooled(f_m1, dl), v), "acute")),
        cbind(outcome = v, label = OUTS[[v]], model = "Without EuroSCORE II (same sample)", n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1),
              row_of(safe(fit_pooled(f_no, dl), v), "acute")))
}))
write_out(C, SUB, "C_acute_vs_elective")

# ---- D. Year trend per group --------------------------------------------------------------------
D <- rbindlist(lapply(c(min(CONFIG$study_years), 2018L), function(start) rbindlist(lapply(names(OUTS), function(v)
  rbindlist(lapply(c(1L, 0L), function(g) {
    dsub <- lapply(cc, function(d) d[acute == g & year >= start & year <= end_of(v) & !is.na(get(v))])
    r <- fit_rule(v, COVARS_M1, dsub)
    cbind(data.table(outcome = v, label = OUTS[[v]], group = if (g == 1) "Acute" else "Elective",
                     from = start, to = end_of(v), model = r$model, n = r$n, events = r$events), row_of(r$fit, "year_c"))
  }))))))
D[outcome != "comp_primary", p_holm := p.adjust(p, "holm"), by = .(group, from)]
write_out(D, SUB, "D_year_trend_per_group")

# ---- E. Piecewise and landmark Cox, acute versus elective -------------------------------------------
LM <- CONFIG$landmark_days
svl <- lapply(cc, function(d) d[!is.na(time_d) & !is.na(status)][, `:=`(t_days = pmax(time_d, 0.5), t_years = pmax(time_d, 0.5) / 365.25)])
pool_cox <- function(formula, dl) {
  fits <- lapply(dl, function(d) coxph(formula, data = d))
  r <- exp_table(rubin_pool(lapply(fits, coef), lapply(fits, vcov)))
  r[, `:=`(n = fits[[1]]$n, events = fits[[1]]$nevent)][]
}
rhs_with <- paste(covar_terms(COVARS_M1, specs), collapse = " + ")
rhs_no   <- paste(covar_terms(c("age", "female", "prior_cs"), specs), collapse = " + ")
vars_m1  <- unique(c("year_c", "acute", "centrum", all.vars(as.formula(paste("~", rhs_with)))))
same     <- function(d) d[complete.cases(d[, .SD, .SDcols = intersect(vars_m1, names(d))])]
split_pw <- function(d) {
  s <- survSplit(Surv(t_years, status) ~ ., data = as.data.frame(d), cut = CONFIG$piecewise_cuts / 365.25, episode = "interval")
  s$interval <- factor(s$interval, labels = c("0-30 days", "30 days-1 year", ">1 year")); as.data.table(s)
}
run_cox <- function(rhs, dl, label) {
  f_pw <- as.formula(paste("Surv(tstart, t_years, status) ~ acute:interval + year_c:interval +", rhs, "+ strata(centrum)"))
  pw <- complete_for(lapply(dl, split_pw), f_pw)
  a <- pool_cox(f_pw, pw)[grepl("acute", term) & grepl("interval", term)]
  a[, interval := sub(".*interval", "", term)]
  f_lm <- as.formula(paste("Surv(t_lm, status) ~ year_c + acute +", rhs, "+ strata(centrum)"))
  lm_l <- complete_for(lapply(dl, function(d) d[t_days > LM][, t_lm := (t_days - LM) / 365.25]), f_lm)
  b <- pool_cox(f_lm, lm_l)[term == "acute"][, interval := "landmark (30-day survivors)"]
  cbind(model = label, patients = nrow(dl[[1]]),
        rbind(a, b, fill = TRUE)[, .(interval, events, hr = ratio, lower = ratio_lower, upper = ratio_upper, p)])
}
E <- rbind(safe(run_cox(rhs_with, lapply(svl, same), "With EuroSCORE II"), "cox with"),
           safe(run_cox(rhs_no,   lapply(svl, same), "Without EuroSCORE II, same patients"), "cox without"))
if (!is.null(E)) write_out(E, SUB, "E_cox_acute_vs_elective")

# ---- Overview ------------------------------------------------------------------------------------
fmt <- function(d, ratio = "OR", col = "or") { d <- copy(d)
  d[, (ratio) := sprintf("%.3f (%.3f-%.3f)", get(col), lower, upper)][, p := signif(p, 3)][] }
out <- c(sprintf("Sensitivity analysis, strict groups (acute vs elective only) - %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
         "", "==== A1. Groepsindeling: primair versus strikt ====", capture.output(print(A1, nrows = 100)),
         "", "==== A2. Per jaar ====", capture.output(print(A2)),
         "", "==== B0. Case mix per groep ====", capture.output(print(B0)),
         "", "==== B. Ruwe uitkomsten per groep ====", capture.output(print(dcast(B, outcome ~ group, value.var = "pct"), nrows = 100)),
         "", "==== C. Acuut vs electief (OR), met en zonder EuroSCORE II (zelfde sample) ====",
         capture.output(print(fmt(C)[, .(outcome, model, n, events, OR, p)], nrows = 100)),
         "", "==== D. OR per jaar per groep ====",
         capture.output(print(fmt(D)[, .(from, group, outcome, to, model, n, events, OR, p, p_holm = signif(p_holm, 3))], nrows = 100)),
         "", "==== E. Cox acuut vs electief ====",
         if (is.null(E)) "Cox skipped" else capture.output(print(fmt(E, "HR", "hr")[, .(model, patients, interval, events, HR, p)], nrows = 50)))
f_out <- file.path(out_dir(SUB), "overzicht_strict_groups.txt")
writeLines(out, f_out); cat(out, sep = "\n")
log_msg("Done. Overview in ", f_out)

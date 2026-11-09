# =============================================================================
# 12_group_comparisons.R
# Additional tables for the manuscript, by group (acute / non-acute):
#  1. Table 1 with p-values (acute vs non-acute), incl. median age.
#  2. Table S1: 2013-2018 vs 2019-2024 within each group, with p-values.
#  3. Table 2 with p-values: temperature, hypothermia category, circulatory arrest.
#  4. Organ-protection trends per group: centre fixed effects with CR2
#     (temperature, arrest use), decomposition of the temperature trend into a
#     within-centre and a between-centre component, random-slope heterogeneity,
#     and the temperature difference without vs with circulatory arrest.
#  5. Centre ranges per group (median temperature, deep hypothermia, arrest use,
#     arrest duration).
#  6. Length of stay per group: medians and adjusted trends (all patients and
#     patients discharged alive).
#  7. Discharge alive without a major complication, per group.
# Complete-case analysis (observed data). Output: results/12_group_comparisons/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("12_group_comparisons")
SUB <- "12_group_comparisons"

cohort <- load_derived("cohort.rds")
specs  <- load_derived("specs.rds")
coh <- cohort[!is.na(urgency2)]
coh[, period := factor(fifelse(year <= 2018, "2013-2018", "2019-2024"), c("2013-2018", "2019-2024"))]
coh[, hypo3 := factor(fcase(temp_min <= 20, "Deep (<=20 C)", temp_min <= 28, "Moderate (20.1-28 C)",
                            temp_min > 28, "Above 28 C"), c("Deep (<=20 C)", "Moderate (20.1-28 C)", "Above 28 C"))]
coh[!(circ_arrest_b %in% 1), ca_min := NA]
coh[, discharge_free := fcase(dzh %in% 1 | major_compl %in% 1, 0L, dzh %in% 0 & major_compl %in% 0, 1L, default = NA_integer_)]
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ---- generic two-group comparison table ---------------------------------------------------
CONT <- c(age = "Age, years", bmi = "BMI, kg/m2", creat = "Creatinine, umol/L", egfr = "eGFR, mL/min/1.73m2",
          euroscore2 = "EuroSCORE II, %", temp_min = "Minimum temperature, C", ecc_min = "CPB time, min",
          aox_min = "Cross-clamp time, min", ca_min = "Circulatory arrest time, min (if arrest)")
BIN  <- c(female = "Female sex", prior_cs = "Previous cardiac surgery", endocarditis_b = "Active endocarditis",
          sternotomy = "Full sternotomy", asc_inv = "Ascending aorta replacement", root_inv = "Aortic root involved",
          desc_inv = "Descending aorta involved", hybrid_fet = "Combined open + endovascular (incl. FET)",
          conc_valve = "Concomitant valve surgery", aortic_valve = "Concomitant aortic valve surgery",
          conc_cabg = "Concomitant CABG", other_cardiac = "Other concomitant cardiac surgery",
          circ_arrest_b = "Circulatory arrest")
FAC  <- c(lvef_cat = "LVEF", urg4 = "Urgency", path3 = "Aortic pathology", hypo3 = "Hypothermia category")
compare_table <- function(d, g, levs) {
  a <- d[[g]] == levs[1]; b <- d[[g]] == levs[2]
  rows <- list()
  for (v in names(CONT)) {
    x <- d[[v]]
    p <- tryCatch(wilcox.test(x[a], x[b])$p.value, error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.table(characteristic = paste0(CONT[[v]], ", median [IQR]"),
                                           g1 = fmt_med(x[a]), g2 = fmt_med(x[b]), missing = sum(is.na(x)), p = fmt_p(p))
    if (v == "age") rows[[length(rows) + 1]] <- data.table(characteristic = "Age, years, mean (SD)",
                                                          g1 = fmt_mean(x[a]), g2 = fmt_mean(x[b]), missing = sum(is.na(x)),
                                                          p = fmt_p(tryCatch(t.test(x[a], x[b])$p.value, error = function(e) NA_real_)))
  }
  for (v in names(BIN)) {
    x <- d[[v]]
    tb <- table(factor(x[a | b], levels = 0:1), d[[g]][a | b])
    p <- tryCatch(suppressWarnings(chisq.test(tb)$p.value), error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.table(characteristic = BIN[[v]], g1 = fmt_npct(x[a]), g2 = fmt_npct(x[b]),
                                           missing = sum(is.na(x)), p = fmt_p(p))
  }
  for (v in names(FAC)) {
    f <- droplevels(as.factor(d[[v]]))
    tb <- table(f[a | b], d[[g]][a | b])
    p <- tryCatch(suppressWarnings(chisq.test(tb[rowSums(tb) > 0, , drop = FALSE])$p.value), error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.table(characteristic = FAC[[v]], g1 = "", g2 = "", missing = sum(is.na(f)), p = fmt_p(p))
    for (l in levels(f)) rows[[length(rows) + 1]] <- data.table(characteristic = paste0("    ", l),
                                                                g1 = fmt_npct(as.integer(f[a] == l)), g2 = fmt_npct(as.integer(f[b] == l)),
                                                                missing = NA_integer_, p = "")
  }
  out <- rbind(data.table(characteristic = "N", g1 = as.character(sum(a, na.rm = TRUE)), g2 = as.character(sum(b, na.rm = TRUE)),
                          missing = NA_integer_, p = ""), rbindlist(rows))
  setnames(out, c("g1", "g2"), levs); out
}

# 1. Table 1 with p-values
t1 <- compare_table(coh, "urgency2", c("Non-acute", "Acute"))
write_out(t1, SUB, "table1_with_p_acute_vs_nonacute")

# 2. Table S1: periods within each group
s1 <- rbindlist(lapply(c("Acute", "Non-acute"), function(gr)
  cbind(group = gr, compare_table(coh[urgency2 == gr], "period", c("2013-2018", "2019-2024")))), fill = TRUE)
write_out(s1, SUB, "tableS1_period_comparison_by_group")

# 3. Table 2 with p-values (organ protection, acute vs non-acute)
t2v <- c("temp_min", "hypo3", "circ_arrest_b", "ca_min", "ecc_min", "aox_min")
t2 <- t1[sapply(characteristic, function(x) any(startsWith(trimws(x), c(CONT[intersect(t2v, names(CONT))], BIN["circ_arrest_b"], FAC["hypo3"],
                                                                         levels(coh$hypo3)))))]
write_out(t2, SUB, "table2_with_p_acute_vs_nonacute")
# temperature by arrest use, per group
ta <- coh[!is.na(circ_arrest_b), .(n = .N, temp_median = fmt_med(temp_min), temp_n = sum(!is.na(temp_min))),
          by = .(urgency2, circulatory_arrest = fifelse(circ_arrest_b == 1, "yes", "no"))][order(urgency2, circulatory_arrest)]
write_out(ta, SUB, "temperature_by_arrest_use_and_group")

# ---- 4. Organ-protection trends per group --------------------------------------------------
rhs <- covar_terms(COVARS_M1, specs)
cr2 <- function(fit, d, term) {
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = d$centrum, test = "Satterthwaite")
  r <- ct[ct$Coef == term, ]
  data.table(estimate = r$beta, se = r$SE, df = r$df_Satt,
             lower = r$beta - qt(.975, r$df_Satt) * r$SE, upper = r$beta + qt(.975, r$df_Satt) * r$SE, p = r$p_Satt)
}
prac <- list(); dec <- list(); slopes <- list(); arrdiff <- list()
for (gr in c("Acute", "Non-acute")) for (st in unique(c(min(CONFIG$study_years), CONFIG$trend_sens_start))) {
  dg <- coh[urgency2 == gr & year >= st]
  lab <- paste0(gr, ", from ", st)
  # temperature: within-centre (FE) and total (no centre) trend; between = total - within
  dt <- dg[complete.cases(dg[, .(temp_min, year_c, age, female, logit_es2, prior_cs, centrum)])]
  f_w <- as.formula(paste("temp_min ~ year_c +", paste(c(rhs, "centrum"), collapse = " + ")))
  f_t <- as.formula(paste("temp_min ~ year_c +", paste(rhs, collapse = " + ")))
  w <- cr2(lm(f_w, dt), dt, "year_c"); t <- cr2(lm(f_t, dt), dt, "year_c")
  prac[[length(prac) + 1]] <- cbind(group = lab, outcome = "Minimum temperature (C/year)",
                                    model = "centre fixed effects, CR2", n = nrow(dt), w)
  dec[[length(dec) + 1]] <- data.table(group = lab, total_trend = t$estimate, total_lower = t$lower, total_upper = t$upper,
                                       total_p = t$p, within_centre = w$estimate, within_lower = w$lower, within_upper = w$upper,
                                       within_p = w$p, between_centre = t$estimate - w$estimate, n = nrow(dt))
  # arrest use
  da <- dg[complete.cases(dg[, .(circ_arrest_b, year_c, age, female, logit_es2, prior_cs, centrum)])]
  fa <- as.formula(paste("circ_arrest_b ~ year_c +", paste(c(rhs, "centrum"), collapse = " + ")))
  a <- tryCatch(cr2(glm(fa, binomial, da), da, "year_c"), error = function(e) NULL)
  if (!is.null(a)) prac[[length(prac) + 1]] <- cbind(group = lab, outcome = "Circulatory arrest use (OR/year)",
                                                     model = "centre fixed effects, CR2", n = nrow(da),
                                                     a[, `:=`(estimate = exp(estimate), lower = exp(lower), upper = exp(upper))])
  if (st == min(CONFIG$study_years)) {
    # random-slope heterogeneity of the temperature trend
    f0 <- as.formula(paste("temp_min ~ year_c +", paste(c(rhs, "(1 | centrum)"), collapse = " + ")))
    f1 <- as.formula(paste("temp_min ~ year_c +", paste(c(rhs, "(1 + year_c | centrum)"), collapse = " + ")))
    m0 <- lmer(f0, dt, REML = FALSE); m1 <- lmer(f1, dt, REML = FALSE)
    lr <- max(0, 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))))
    vc <- as.data.frame(VarCorr(m1))
    slopes[[length(slopes) + 1]] <- data.table(group = gr, sd_centre_slopes = vc$sdcor[vc$var1 == "year_c" & is.na(vc$var2)][1],
                                               lr = lr, p_mixture = 0.5 * pchisq(lr, 1, lower.tail = FALSE) + 0.5 * pchisq(lr, 2, lower.tail = FALSE))
    # temperature without vs with circulatory arrest (within-centre)
    dd <- dt[!is.na(circ_arrest_b)][, no_arrest := as.integer(circ_arrest_b == 0)]
    fd <- as.formula(paste("temp_min ~ no_arrest + year_c +", paste(c(rhs, "centrum"), collapse = " + ")))
    if (length(unique(dd$no_arrest)) == 2)
      arrdiff[[length(arrdiff) + 1]] <- cbind(group = gr, comparison = "no arrest vs arrest (C)", n = nrow(dd),
                                              n_no_arrest = sum(dd$no_arrest), cr2(lm(fd, dd), dd, "no_arrest"))
  }
}
write_out(rbindlist(prac, fill = TRUE), SUB, "organ_protection_trends_by_group_cr2")
write_out(rbindlist(dec), SUB, "temperature_trend_within_between_decomposition")
write_out(rbindlist(slopes), SUB, "temperature_random_slope_heterogeneity_by_group")
write_out(rbindlist(arrdiff), SUB, "temperature_difference_no_arrest_vs_arrest")

# ---- 5. Centre ranges per group ----------------------------------------------------------
cg <- coh[, .(n = .N, temp_median = as.numeric(median(temp_min, na.rm = TRUE)), temp_n = sum(!is.na(temp_min)),
              pct_deep = 100 * mean(temp_min <= 20, na.rm = TRUE), pct_arrest = 100 * mean(circ_arrest_b == 1, na.rm = TRUE),
              ca_median = as.numeric(median(ca_min, na.rm = TRUE))), by = .(urgency2, centre_lab)][order(urgency2, centre_lab)]
write_out(cg, SUB, "centre_practice_by_group")
rng <- cg[temp_n >= 10, .(centres = .N, temp_median_min = min(temp_median), temp_median_max = max(temp_median),
                          pct_deep_min = min(pct_deep), pct_deep_max = max(pct_deep))
          , by = urgency2]
rng2 <- cg[n >= 10, .(arrest_pct_min = min(pct_arrest, na.rm = TRUE), arrest_pct_max = max(pct_arrest, na.rm = TRUE),
                      ca_median_min = min(ca_median, na.rm = TRUE), ca_median_max = max(ca_median, na.rm = TRUE)), by = urgency2]
write_out(merge(rng, rng2, by = "urgency2"), SUB, "centre_practice_ranges_by_group")

# ---- 6. Length of stay per group -----------------------------------------------------------
los_desc <- coh[!is.na(los), .(n = .N, los_all = fmt_med(los), n_alive = sum(dzh == 0, na.rm = TRUE),
                               los_discharged_alive = fmt_med(los[dzh %in% 0])), by = urgency2]
write_out(los_desc, SUB, "length_of_stay_by_group")
los_tr <- list()
for (gr in c("Acute", "Non-acute")) for (pop in c("all patients", "discharged alive")) {
  d <- coh[urgency2 == gr & !is.na(los) & los > 0]
  if (pop == "discharged alive") d <- d[dzh %in% 0]
  d[, log_los := log(los)]
  f <- as.formula(paste("log_los ~ year_c +", paste(c(rhs, "(1 | centrum)"), collapse = " + ")), env = globalenv())
  r <- tryCatch(fit_pooled(f, list(d), family = NULL), error = function(e) NULL)
  if (is.null(r)) next
  e <- exp_table(r$pooled, "year_c")
  los_tr[[length(los_tr) + 1]] <- data.table(group = gr, population = pop, n = r$diag$nobs,
                                             ratio_per_year = e$ratio, lower = e$ratio_lower, upper = e$ratio_upper, p = e$p)
}
write_out(rbindlist(los_tr), SUB, "length_of_stay_trend_by_group")

# ---- 7. Discharge alive without a major complication ----------------------------------------
write_out(coh[, .(discharge_alive_without_major_complication = fmt_npct(discharge_free)), by = urgency2], SUB,
          "discharge_without_major_complication_by_group")

log_msg("Done.")

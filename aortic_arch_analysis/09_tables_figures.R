# =============================================================================
# 09_tables_figures.R
# Publication tables and figures generated directly from the fitted analyses,
# cross-checks of consistency between tables, figures and models, and an
# automatically generated results report.
# Output: results/publication/*, results/analysis_report.md
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
suppressPackageStartupMessages(library(patchwork))
start_log("09_tables_figures")
SUB <- "publication"
RES <- CONFIG$results_dir
rd <- function(sub, name) fread(file.path(res_path(sub), paste0(name, ".csv")))
has <- function(sub, name) file.exists(file.path(res_path(sub), paste0(name, ".csv")))
log_msg("Missing-data method: ", method_label())

cohort  <- load_derived("cohort.rds")
outc    <- readRDS(derived_method("outcomes.rds"))
organ   <- readRDS(derived_method("organ.rds"))
surv    <- readRDS(derived_method("survival.rds"))
vol     <- load_derived("volume.rds")
popf    <- derived_path("population.rds"); pop <- if (file.exists(popf)) readRDS(popf) else NULL
sensf   <- derived_method("sensitivity.rds"); sens <- if (file.exists(sensf)) readRDS(sensf) else NULL
checks  <- list()
chk <- function(what, ok, detail = "") checks[[length(checks) + 1]] <<- data.table(check = what, passed = isTRUE(ok), detail = detail)

# =============================== TABLES ==========================================================
# Table 1
t1 <- rd("02_descriptives", "table1_baseline_by_urgency")
setnames(t1, c("characteristic", "Overall", "Non-acute", "Acute", "missing", "smd"),
         c("Characteristic", "Overall", "Non-acute", "Acute", "Missing (n)", "SMD"), skip_absent = TRUE)
write_out(t1, SUB, "Table1_baseline_characteristics")
chk("Table 1 N (overall) + not classifiable = cohort size",
    as.integer(t1$Overall[1]) + sum(is.na(cohort$urgency2)) == nrow(cohort),
    sprintf("%s + %d vs %d", t1$Overall[1], sum(is.na(cohort$urgency2)), nrow(cohort)))

# Table 2
t2 <- copy(organ$t2)
write_out(t2, SUB, "Table2_organ_protection")
chk("Table 2 total N = patients with CPB (or CPB unknown)",
    t2[era == "Total", n] == nrow(cohort[ecc_any == 1 | is.na(ecc_any)]))

# Table 3: crude rates by urgency and era, adjusted year trend and urgency OR
cr <- outc$crude
fmt_rate <- function(e, n) sprintf("%d/%d (%.1f%%)", e, n, 100 * e / n)
t3 <- cr[group_type %in% c("overall", "acute status", "era"),
         .(cell = fmt_rate(events, n)), by = .(outcome, label, group)]
t3 <- dcast(t3, outcome + label ~ group, value.var = "cell")
yr <- outc$year[, .(outcome, adj_or_year = fmt_ci(or_per_year, lower, upper, 3), p_year = fmt_p(p),
                    p_holm = fmt_p(p_holm), p_nonlinear = fmt_p(p_nonlinear), model)]
ur <- outc$urgency[adjustment == "without pathology",
                   .(outcome, adj_or_acute = fmt_ci(ratio, ratio_lower, ratio_upper, 2))]
t3 <- merge(merge(t3, yr, by = "outcome", all.x = TRUE), ur, by = "outcome", all.x = TRUE)
t3 <- t3[match(outc$OUT$var, outcome)][!is.na(outcome)]
setcolorder(t3, intersect(c("label", "All", "Non-acute", "Acute", CONFIG$era_labels, "adj_or_year",
                            "p_year", "p_holm", "p_nonlinear", "adj_or_acute", "model", "outcome"), names(t3)))
write_out(t3, SUB, "Table3_early_outcomes")
# consistency: events in Table 3 equal events used in models (M1, imputation-independent outcome)
ev_tab <- cr[group_type == "overall", .(outcome, events)]
ev_mod <- unique(outc$models[term == "year_c", .(outcome, model_events = events)])
evc <- merge(ev_tab, ev_mod, by = "outcome")
if (CONFIG$missing_method == "imputation") {
  chk("Table 3 event counts equal events in fitted models", all(evc$events == evc$model_events),
      paste(evc[events != model_events, outcome], collapse = ","))
} else {
  chk("Complete-case models use no more events than Table 3", all(evc$model_events <= evc$events),
      paste(sprintf("%s: %d of %d events in model", evc$outcome, evc$model_events, evc$events), collapse = "; "))
}
chk("Table 3 30-day deaths equal cohort count",
    ev_tab[outcome == "d30", events] == sum(cohort$d30 == 1, na.rm = TRUE))

# Supplementary workbook
supp <- list()
add_supp <- function(name, sub, file) if (has(sub, file)) supp[[name]] <<- rd(sub, file)
add_supp("S1_flow", "02_descriptives", "flow_diagram")
add_supp("S2_missing_overall", "02_descriptives", "missing_overall")
add_supp("S3_casemix_trends", "02_descriptives", "characteristics_trend_per_year")
add_supp("S4_volume_models", "03_volume", "volume_trend_models")
add_supp("S5_volume_interaction", "03_volume", "volume_interaction_urgency")
add_supp("S6_outcome_models_all_terms", "04_outcomes", "adjusted_models_M1_all_terms")
add_supp("S7_year_trend_summary", "04_outcomes", "year_trend_summary")
add_supp("S8_acute_vs_nonacute", "04_outcomes", "acute_vs_nonacute_or")
add_supp("S9_acute_interaction", "04_outcomes", "acute_vs_nonacute_year_interaction")
add_supp("S10_centre_variation", "04_outcomes", "centre_variation")
add_supp("S11_los", "04_outcomes", "length_of_stay_models")
add_supp("S12_ftr_check", "04_outcomes", "ftr_numerator_denominator_check")
add_supp("S13_ftr_decomposition", "04_outcomes", "ftr_decomposition_by_year")
add_supp("S14_ftr_models", "04_outcomes", "ftr_models")
add_supp("S15_practice_models", "05_organ_protection", "practice_trend_models")
add_supp("S16_hypothermia_ordinal", "05_organ_protection", "hypothermia_ordinal_models")
add_supp("S17_random_slopes", "05_organ_protection", "random_slope_lrt_and_icc")
add_supp("S18_within_centre_FE", "05_organ_protection", "within_centre_fixed_effects_cr2")
add_supp("S19_km", "06_survival", "km_survival_estimates")
add_supp("S20_cox", "06_survival", "cox_models")
add_supp("S21_ph_test", "06_survival", "cox_ph_test_schoenfeld")
add_supp("S22_piecewise", "06_survival", "cox_piecewise_hazard_ratios")
add_supp("S23_rmst", "06_survival", "rmst_by_acute_status")
add_supp("S24_landmark_cox", "06_survival", "landmark30_cox")
add_supp("S25_conditional", "06_survival", "conditional_survival")
add_supp("S26_smr", "07_population", "smr_overall_and_by_group")
add_supp("S27_smr_interval", "07_population", "smr_and_excess_by_followup_interval")
add_supp("S28_excess_hazard", "07_population", "excess_hazard_regression")
add_supp("S29_timescale_check", "07_population", "time_scale_verification")
add_supp("S30_sensitivity", "08_sensitivity", "sensitivity_year_effect")
add_supp("S31_sens_practice", "08_sensitivity", "sensitivity_practice")
add_supp("S32_convergence", "04_outcomes", "model_convergence_diagnostics")
writexl::write_xlsx(c(list(Table1 = t1, Table2 = t2, Table3 = t3), supp),
                    file.path(out_dir(SUB), "Tables_and_Supplement.xlsx"))

# =============================== FIGURES =========================================================
# Figure 1: national volume
cnt <- vol$cnt
p1a <- ggplot(cnt[type != "All"], aes(factor(year), n, fill = type)) + geom_col(width = .75) +
  scale_fill_manual(values = PAL[c(1, 2)], name = NULL) +
  labs(x = NULL, y = "Procedures (n)", title = "A. Annual procedures") + theme_pub()
if (vol$has_pop) {
  fv <- vol$fitted
  p1b <- ggplot(cnt, aes(year, rate_per_million, colour = type)) + geom_point() +
    geom_line(data = fv, aes(year, fit, colour = series), inherit.aes = FALSE) +
    geom_ribbon(data = fv, aes(year, ymin = lower, ymax = upper, fill = series), inherit.aes = FALSE, alpha = .12) +
    scale_colour_manual(values = PAL[c(6, 2, 1)], breaks = c("All", "Acute", "Non-acute"), name = NULL) +
    scale_fill_manual(values = PAL[c(6, 2, 1)], breaks = c("All", "Acute", "Non-acute"), guide = "none") +
    scale_x_continuous(breaks = CONFIG$study_years) +
    labs(x = NULL, y = "Per million inhabitants", title = "B. Population-adjusted rate (fitted trend, 95% CI)") + theme_pub()
  save_fig(p1a / p1b, SUB, "Figure1_volume", 8, 8)
} else save_fig(p1a, SUB, "Figure1_volume", 8, 4.5)
chk("Figure 1 counts sum to cohort size (all years)", sum(cnt[type == "All", n]) == nrow(cohort))

# Figure 2: case-mix over time
an <- rd("02_descriptives", "characteristics_by_year")
anl <- melt(an, id.vars = "year", measure.vars = c("acute", "acute_dissection", "prior_cs", "hybrid_fet",
                                                   "desc_inv", "conc_valve", "conc_cabg", "female"))
lab2 <- c(acute = "Acute (study definition)", acute_dissection = "Acute dissection", prior_cs = "Previous cardiac surgery",
          hybrid_fet = "Open + endovascular (FET)", desc_inv = "Descending aorta involved",
          conc_valve = "Concomitant valve surgery", conc_cabg = "Concomitant CABG", female = "Female")
anl[, label := factor(lab2[as.character(variable)], lab2)]
p2 <- ggplot(anl, aes(year, 100 * value)) + geom_line(colour = PAL[1]) + geom_point(colour = PAL[1], size = 1.2) +
  facet_wrap(~label, ncol = 4) + scale_x_continuous(breaks = c(2013, 2018, 2024)) +
  labs(x = NULL, y = "%", title = "Patient characteristics and operative complexity by year") + theme_pub()
save_fig(p2, SUB, "Figure2_casemix_trends", 9, 5)

# Figure 3: adjusted (standardised) risk by year with crude points
gc <- outc$gcomp[outcome %in% c("d30", "dzh", "stroke_def", "renal", "major_compl", "vent24")]
crude_y <- outc$crude[group_type == "year" & outcome %in% unique(gc$outcome)][, year := as.integer(group)]
gc[, panel := factor(label, unique(label))]; crude_y[, panel := factor(label, levels(gc$panel))]
p3 <- ggplot(gc, aes(year, 100 * risk)) +
  geom_ribbon(aes(ymin = 100 * lower, ymax = 100 * upper), fill = PAL[1], alpha = .18) +
  geom_line(colour = PAL[1]) + geom_point(data = crude_y, aes(year, pct), colour = "grey35", size = 1.3) +
  facet_wrap(~panel, scales = "free_y", ncol = 3) + scale_x_continuous(breaks = c(2013, 2018, 2024)) +
  labs(x = NULL, y = "Risk (%)",
       title = "Risk-adjusted outcomes by year",
       subtitle = "Line: standardised risk (g-computation, spline in year, 95% interval); points: crude") + theme_pub()
save_fig(p3, SUB, "Figure3_adjusted_trends", 9, 6)

# Figure 4: organ protection
hy <- copy(organ$hypo_year)
p4a <- ggplot(hy, aes(factor(year), pct, fill = hypo4)) + geom_col(width = .8) +
  scale_fill_manual(values = c("#08306b", "#2a6ebb", "#9ecae1", "#edae49"), name = NULL) +
  labs(x = NULL, y = "%", title = "A. Lowest core temperature category") + theme_pub() +
  guides(fill = guide_legend(nrow = 2))
ann <- organ$annual
p4b <- ggplot(ann, aes(year, temp_median)) + geom_ribbon(aes(ymin = temp_q1, ymax = temp_q3), fill = PAL[1], alpha = .15) +
  geom_line(colour = PAL[1]) + geom_point(colour = PAL[1]) + scale_x_continuous(breaks = c(2013, 2018, 2024)) +
  labs(x = NULL, y = "Median (IQR), C", title = "B. Minimum core temperature") + theme_pub()
cy <- organ$centre_year
p4c <- ggplot(cy, aes(factor(year), centre_lab, fill = temp_median)) + geom_tile(colour = "white") +
  scale_fill_gradient(low = "#08306b", high = "#fdd49e", name = "Median C") +
  labs(x = NULL, y = "Centre", title = "C. Median minimum temperature by centre and year") + theme_pub() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = .5))
sl <- organ$slopes[grepl("^temp_min \\| CPB \\| adjusted", key)]
p4d <- ggplot(sl, aes(reorder(centre_lab, slope), slope)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), colour = PAL[1], size = .3) +
  coord_flip() + labs(x = "Centre", y = "C per year (BLUP)", title = "D. Within-centre temperature trend") + theme_pub()
save_fig((p4a | p4b) / (p4c | p4d), SUB, "Figure4_organ_protection", 11, 9)

# Figure 5: centre variation in 30-day mortality
ct <- outc$caterpillar[outcome == "d30"]; fu <- outc$funnel[outcome == "d30"]
p5a <- ggplot(ct, aes(reorder(centre_lab, effect), or)) + geom_hline(yintercept = 1, linetype = 2) +
  geom_pointrange(aes(ymin = lower, ymax = upper), colour = PAL[1], size = .3) + coord_flip() + scale_y_log10() +
  labs(x = "Centre", y = "Centre odds ratio (conditional mode, 95% interval)", title = "A. Centre effects, 30-day mortality") + theme_pub()
ex <- seq(max(1, min(fu$expected) * 0.8), max(fu$expected) * 1.2, length.out = 200)
lims <- data.table(expected = ex, lo95 = qpois(.025, ex) / ex, hi95 = qpois(.975, ex) / ex,
                   lo998 = qpois(.001, ex) / ex, hi998 = qpois(.999, ex) / ex)
p5b <- ggplot(fu, aes(expected, oe)) + geom_hline(yintercept = 1) +
  geom_line(data = lims, aes(expected, lo95), linetype = 2, colour = "grey40") +
  geom_line(data = lims, aes(expected, hi95), linetype = 2, colour = "grey40") +
  geom_line(data = lims, aes(expected, lo998), linetype = 3, colour = "grey40") +
  geom_line(data = lims, aes(expected, hi998), linetype = 3, colour = "grey40") +
  geom_point(colour = PAL[2]) + geom_text(aes(label = centre_lab), nudge_y = .06, size = 3) +
  labs(x = "Expected deaths (case-mix model)", y = "Observed / expected",
       title = "B. Funnel plot (95% and 99.8% Poisson limits)") + theme_pub()
save_fig(p5a | p5b, SUB, "Figure5_centre_variation", 11, 5)
chk("Funnel plot observed deaths sum to model events (30-day mortality)",
    sum(fu$observed) == outc$year[outcome == "d30", events])

# Figures 6-7: Kaplan-Meier and landmark curves with risk tables
km_plot <- function(cd, rk, title, cols = PAL) {
  p <- ggplot(cd, aes(years, surv, colour = strata, fill = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) + geom_step() +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + scale_x_continuous(breaks = 0:12) +
    scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, guide = "none") +
    labs(x = "Years after surgery", y = "Survival", title = title) + theme_pub()
  rt <- ggplot(rk, aes(years, strata, label = n_risk)) + geom_text(size = 2.8) +
    scale_x_continuous(breaks = 0:12, limits = range(cd$years)) + labs(x = NULL, y = NULL, title = "Number at risk") +
    theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), plot.title = element_text(size = 9))
  p / rt + plot_layout(heights = c(4, 1))
}
save_fig(km_plot(surv$curves$overall, surv$risks$overall, "Overall survival"), SUB, "Figure6a_KM_overall", 7, 6)
save_fig(km_plot(surv$curves$urgency, surv$risks$urgency, "Survival: acute vs non-acute"), SUB, "Figure6b_KM_acute_status", 7, 6)
save_fig(km_plot(surv$curves$pathology, surv$risks$pathology, "Survival by aortic pathology"), SUB, "Figure6c_KM_pathology", 7, 6)
save_fig(km_plot(surv$lm_curves$urgency, surv$lm_risks$urgency, "Survival among 30-day survivors (landmark at 30 days)"),
         SUB, "Figure7_landmark_acute_status", 7, 6)
chk("KM number at risk at time 0 equals survival population",
    surv$risks$overall[years == 0, n_risk] == nrow(cohort[!is.na(time_d) & !is.na(status)]))

# Figure 8: relative survival
if (!is.null(pop)) {
  rs <- pop$rs[group %in% c("All", "Non-acute", "Acute")]
  rsl <- rbind(rs[, .(group, years = end, value = observed, type = "Observed")],
               rs[, .(group, years = end, value = expected, type = "Expected (general population)")],
               rs[, .(group, years = end, value = relative, type = "Relative (Ederer II)")])
  p8a <- ggplot(rsl, aes(years, value, colour = type)) + geom_line() + facet_wrap(~group) +
    scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = PAL[c(1, 6, 2)], name = NULL) +
    labs(x = "Years after surgery", y = NULL, title = "A. Observed, expected and relative survival") + theme_pub()
  pp <- pop$pp[group %in% c("All", "Non-acute", "Acute")]
  p8b <- ggplot(pp, aes(years, net_survival, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) + geom_line() +
    scale_colour_manual(values = PAL[c(6, 1, 2)], name = NULL) + scale_fill_manual(values = PAL[c(6, 1, 2)], guide = "none") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Years after surgery", y = "Net survival", title = "B. Net survival (Pohar-Perme, sensitivity)") + theme_pub()
  si <- pop$smr_int[group %in% c("All", "Group: Non-acute", "Group: Acute")]
  si[, interval := factor(interval, c("0-30 days", "30 days-1 year", "1-5 years", ">5 years"))]
  p8c <- ggplot(si, aes(interval, smr, colour = group)) + geom_hline(yintercept = 1, linetype = 2) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(.5), size = .3) + scale_y_log10() +
    scale_colour_manual(values = PAL[c(6, 1, 2)], name = NULL) +
    labs(x = "Follow-up interval", y = "SMR (log scale)", title = "C. Standardised mortality ratio by interval") + theme_pub()
  save_fig(p8a / (p8b | p8c), SUB, "Figure8_relative_survival", 11, 9)
}

# Supplementary figures: sensitivity forest; FTR decomposition; conditional survival
if (!is.null(sens)) {
  s30 <- sens[outcome == "d30" & term == "year_c"]
  s30[, analysis := factor(analysis, rev(unique(analysis)))]
  ps <- ggplot(s30, aes(or, analysis)) + geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper), colour = PAL[1], size = .3) +
    labs(x = "Odds ratio per calendar year (95% CI)", y = NULL, title = "Sensitivity analyses: 30-day mortality trend") + theme_pub()
  save_fig(ps, SUB, "FigureS1_sensitivity_forest", 8, 5)
}
fd <- melt(outc$ftr_decomp, id.vars = "year", measure.vars = c("p_complication", "ftr", "mort_no_compl", "mortality"))
fd[, variable := factor(variable, c("p_complication", "ftr", "mort_no_compl", "mortality"),
                        c("Major complication", "Failure to rescue", "Mortality without complication", "30-day mortality"))]
pf <- ggplot(fd, aes(year, 100 * value)) + geom_line(colour = PAL[1]) + geom_point(colour = PAL[1], size = 1.2) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) + scale_x_continuous(breaks = c(2013, 2018, 2024)) +
  labs(x = NULL, y = "%", title = "Decomposition of 30-day mortality (crude)") + theme_pub()
save_fig(pf, SUB, "FigureS2_ftr_decomposition", 10, 3.5)
cs <- copy(surv$cond)[, s_lab := factor(round(survived_days), c(0, 30, 365, 730, 1826),
                                       c("0 d", "30 d", "1 y", "2 y", "5 y"))]
pc <- ggplot(cs, aes(s_lab, cond_surv, colour = group)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(.4), size = .3) +
  facet_wrap(~paste("Additional", additional_years, "year(s)")) + scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = PAL[c(6, 1, 2)], name = NULL) +
  labs(x = "Already survived", y = "Conditional survival", title = "Conditional survival") + theme_pub()
save_fig(pc, SUB, "FigureS3_conditional_survival", 8, 4.5)

# Model diagnostics figure (primary model binned residuals)
if (has("04_outcomes", "primary_model_binned_residuals")) {
  br <- rd("04_outcomes", "primary_model_binned_residuals")
  pb <- ggplot(br, aes(mean_fitted, mean_resid)) + geom_hline(yintercept = 0) +
    geom_line(aes(y = se), linetype = 2, colour = "grey50") + geom_line(aes(y = -se), linetype = 2, colour = "grey50") +
    geom_point(colour = PAL[1]) + labs(x = "Mean fitted probability (bin)", y = "Mean residual",
                                       title = "Binned residuals, primary model (30-day mortality)") + theme_pub()
  save_fig(pb, SUB, "FigureS4_binned_residuals", 6, 4)
}

# =============================== CONSISTENCY & REPORT ===========================================
cks <- rbindlist(checks)
write_out(cks, SUB, "consistency_checks_tables_figures_models")
if (!all(cks$passed)) log_msg("WARNING: consistency check(s) failed: ", paste(cks[passed == FALSE, check], collapse = "; "))

g <- function(x, d = 2) formatC(x, format = "f", digits = d)
yr <- outc$year; flow <- rd("02_descriptives", "flow_diagram")
d30 <- yr[outcome == "d30"]
vr <- vol$res
conv <- if (has("04_outcomes", "model_convergence_diagnostics")) rd("04_outcomes", "model_convergence_diagnostics") else NULL
pr <- organ$practice
tmp <- pr[outcome == "Minimum temperature (C)" & population == "CPB" & adjusted == TRUE]
report <- c(
  sprintf("# Analytical report (auto-generated): %s", method_label()), "",
  sprintf("Adjusted analyses use **%s**. Descriptive results (cohort, Table 1, volume) and population-referenced survival do not depend on this choice.", method_label()), "",
  sprintf("Generated %s. Input: `%s`. Every number below is read from the fitted analyses; no values were typed by hand.",
          format(Sys.time()), basename(CONFIG$raw_data_path)), "",
  "", "## Cohort", "",
  sprintf("* Records in extract: %d; analysis cohort: %d patients from %d centres (%d-%d).",
          flow[step == "Records in extract", n], nrow(cohort), uniqueN(cohort$centrum),
          min(cohort$year), max(cohort$year)),
  sprintf("* Acute (study definition): %d (%.1f%% of classifiable); not classifiable: %d; acute aortic dissection: %d (%.1f%%).",
          sum(cohort$acute == 1, na.rm = TRUE), 100 * mean(cohort$acute, na.rm = TRUE), sum(is.na(cohort$acute)),
          sum(cohort$acute_dissection == 1, na.rm = TRUE), 100 * mean(cohort$acute_dissection, na.rm = TRUE)),
  "", "## Volume", "",
  sprintf("* %s: %s model; IRR per year %s (%s), i.e. %+.1f%% per year; non-linearity p = %s.",
          vr$series, vr$model, g(vr$irr_per_year, 3), paste0(g(vr$irr_lower, 3), "-", g(vr$irr_upper, 3)),
          vr$annual_pct_change, fmt_p(vr$p_nonlinear)),
  "", "## Primary endpoint: 30-day mortality", "",
  with(outc$crude[group_type == "overall" & outcome == "d30"],
       sprintf("* Events (all patients with known 30-day status): %d of %d (%.1f%%).", events, n, pct)),
  sprintf("* Adjusted OR per calendar year (M1, %s, centre random intercept; %d patients, %d events in model): %s (95%% CI %s-%s), p = %s; test for non-linearity p = %s.",
          method_label(), d30$n, d30$events,
          g(d30$or_per_year, 3), g(d30$lower, 3), g(d30$upper, 3), fmt_p(d30$p), fmt_p(d30$p_nonlinear)),
  "", "## Secondary early outcomes (OR per year; Holm-adjusted p)", "",
  sprintf("* %s: %s (%s-%s); p = %s; Holm p = %s%s", yr[role == "secondary"]$label,
          g(yr[role == "secondary"]$or_per_year, 3), g(yr[role == "secondary"]$lower, 3), g(yr[role == "secondary"]$upper, 3),
          fmt_p(yr[role == "secondary"]$p), fmt_p(yr[role == "secondary"]$p_holm),
          ifelse(grepl("reduced", yr[role == "secondary"]$model), " (reduced model: EPV < 10)", "")),
  "", "## Organ protection", "",
  if (nrow(tmp)) sprintf("* Minimum temperature, adjusted change per year (random intercept + slope): %s C (%s to %s).",
                         g(tmp$estimate), g(tmp$lower), g(tmp$upper)) else "",
  sprintf("* %s, %s: OR per year %s (%s-%s).", organ$ordinal$outcome,
          ifelse(organ$ordinal$adjusted, "case-mix adjusted", "unadjusted"), g(organ$ordinal$estimate, 3),
          g(organ$ordinal$lower, 3), g(organ$ordinal$upper, 3)),
  "", "## Survival", "",
  sprintf("* %s: %s", surv$fu$measure, g(surv$fu$value)),
  sprintf("* KM survival (overall) at %s: %.1f%% (%.1f-%.1f).", surv$km[analysis == "overall"]$time_label,
          100 * surv$km[analysis == "overall"]$surv, 100 * surv$km[analysis == "overall"]$lower,
          100 * surv$km[analysis == "overall"]$upper),
  if (!is.null(pop)) c("", "## Population-referenced survival", "",
    sprintf("* SMR %s: %s (%s-%s).", pop$smr$group, g(pop$smr$smr), g(pop$smr$lower), g(pop$smr$upper)),
    sprintf("* Time-scale verification passed for %d reference profiles (max abs. difference %.1e).",
            nrow(pop$verify), max(pop$verify$abs_difference))) else "",
  "", "## Model diagnostics", "",
  if (!is.null(conv)) sprintf("* %s / %s: %d fits, %d failed, %d convergence warnings, %d singular.",
                              conv$outcome, conv$model, conv$n_fits, conv$n_failed, conv$n_convergence_warning, conv$n_singular) else "",
  "", "## Consistency checks", "",
  sprintf("* [%s] %s %s", ifelse(cks$passed, "x", " "), cks$check, ifelse(nzchar(cks$detail), paste0("(", cks$detail, ")"), "")))
writeLines(report, file.path(RES, CONFIG$missing_method, "analysis_report.md"))
write_session_info()
log_msg("Done.")

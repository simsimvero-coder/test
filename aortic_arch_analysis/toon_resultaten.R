# =============================================================================
# toon_resultaten.R — shows the main results in the R console, section by section,
# and writes the same overview to results/overzicht_resultaten.txt.
# Run after run_all.R, from the project folder:   source("toon_resultaten.R")
# =============================================================================

M <- Sys.getenv("ARCH_TOON_METHOD", "complete_case")
options(width = 250)
RES <- Sys.getenv("ARCH_RESULTS", "results")
out_file <- file.path(RES, "overzicht_resultaten.txt")
toon <- function(titel, f, kol = NULL, rij = NULL) {
  pad <- file.path(RES, f)
  cat("\n\n==== ", titel, " ====\n", sep = "")
  if (!file.exists(pad)) { cat("(niet beschikbaar: ", f, ")\n", sep = ""); return(invisible()) }
  x <- read.csv(pad, check.names = FALSE)
  if (!is.null(rij)) x <- x[eval(rij, x), , drop = FALSE]
  if (!is.null(kol)) x <- x[, intersect(kol, names(x)), drop = FALSE]
  num <- vapply(x, is.numeric, logical(1))
  rnd <- num & vapply(x, function(z) is.numeric(z) && any(abs(z - round(z)) > 1e-9, na.rm = TRUE), logical(1))
  x[rnd] <- lapply(x[rnd], function(z) ifelse(abs(z) >= 100, round(z, 1), signif(z, 3)))
  print(x, right = FALSE, row.names = FALSE)
}
mm <- function(sub, f) file.path(M, sub, f)

while (sink.number() > 0) sink()
sink(out_file, split = TRUE)
cat("Overzicht resultaten -", format(Sys.time()), "- methode:", M, "\n")

toon("1. Cohort (flowchart)", "02_descriptives/flow_diagram.csv")
toon("2a. Volume per jaar", "03_volume/annual_counts.csv")
toon("2b. Volumetrend per miljoen inwoners", "03_volume/volume_trend_models.csv",
     c("series", "analysis", "model", "irr_per_year", "irr_lower", "irr_upper", "annual_pct_change"))
toon("3. Tabel 1 (acuut vs niet-acuut, met p)", "12_group_comparisons/table1_with_p_acute_vs_nonacute.csv")
toon("4. Tabel S1 (2013-2018 vs 2019-2024 per groep)", "12_group_comparisons/tableS1_period_comparison_by_group.csv")
toon("5a. Tabel 2 (orgaanbescherming, met p)", "12_group_comparisons/table2_with_p_acute_vs_nonacute.csv")
toon("5b. Temperatuur met/zonder circulatiestilstand", "12_group_comparisons/temperature_by_arrest_use_and_group.csv")
toon("5c. Bereik per centrum", "12_group_comparisons/centre_practice_ranges_by_group.csv")
toon("6a. Trends werkwijze per groep (CR2)", "12_group_comparisons/organ_protection_trends_by_group_cr2.csv",
     c("group", "outcome", "n", "estimate", "lower", "upper", "p"))
toon("6b. Temperatuurtrend binnen/tussen centra", "12_group_comparisons/temperature_trend_within_between_decomposition.csv",
     c("group", "total_trend", "within_centre", "between_centre", "within_p"))
toon("6c. Spreiding temperatuurtrend tussen centra", "12_group_comparisons/temperature_random_slope_heterogeneity_by_group.csv")
toon("6d. Temperatuurverschil zonder vs met stilstand", "12_group_comparisons/temperature_difference_no_arrest_vs_arrest.csv",
     c("group", "n", "n_no_arrest", "estimate", "lower", "upper", "p"))
toon("7. Tabel 3 (uitkomsten)", mm("publication", "Table3_early_outcomes.csv"),
     c("outcome", "All", "Non-acute", "Acute", "2013-2016", "2021-2024", "adj_or_year", "p_year", "p_holm", "p_nonlinear"))
toon("8. Trends per groep, vanaf 2013 en 2015", mm("04_outcomes", "year_trend_by_group_and_start_table.csv"))
toon("9a. FTR", mm("04_outcomes", "ftr_models.csv"),
     c("definition", "period", "n", "events", "ratio", "ratio_lower", "ratio_upper", "p"), quote(term == "year_c"))
toon("9b. Ligduur", "12_group_comparisons/length_of_stay_by_group.csv")
toon("9c. Ligduurtrend", "12_group_comparisons/length_of_stay_trend_by_group.csv")
toon("9d. Levend ontslag zonder ernstige complicatie", "12_group_comparisons/discharge_without_major_complication_by_group.csv")
toon("10a. Overleving (Kaplan-Meier)", mm("06_survival", "km_survival_estimates.csv"),
     c("strata", "time_label", "surv", "lower", "upper", "n_risk"), quote(analysis == "urgency"))
toon("10b. Relatieve overleving (Ederer II)", "07_population/relative_survival_ederer2.csv",
     c("group", "end", "observed", "expected", "relative", "lower", "upper"),
     quote(end %in% c(1, 3, 5) & group %in% c("All", "Acute", "Non-acute")))
toon("10c. SMR", "07_population/smr_and_excess_by_followup_interval.csv",
     c("group", "interval", "observed", "expected", "smr", "lower", "upper"), quote(!grepl("Age", group)))
toon("11a. Werkwijze en uitkomst per centrum", mm("11_centre_practice", "centre_profile_practice_and_outcomes.csv"))
toon("11b. Werkwijze vs uitkomst op centrumniveau", mm("11_centre_practice", "centre_level_practice_vs_outcome.csv"),
     c("outcome", "practice", "spearman_rho", "p_spearman", "oe_ratio_per_sd", "lower", "upper", "p_weighted_regression"))
toon("12. Werkwijze binnen/tussen centra", mm("11_centre_practice", "within_between_centre_practice_models.csv"),
     c("outcome", "practice", "adjustment", "component", "or", "lower", "upper", "p"))
toon("13. Werkwijze vs uitkomst zonder centrum-effect", mm("11_centre_practice", "practice_outcome_without_centre_effect.csv"))
toon("14a. Cerebrale perfusie: beschrijving", mm("11_centre_practice", "sacp_descriptives_by_type.csv"))
toon("14b. Cerebrale perfusie: per jaar", mm("11_centre_practice", "sacp_type_by_year.csv"))
toon("14c. Cerebrale perfusie: type", mm("11_centre_practice", "sacp_models.csv"),
     c("analysis", "outcome", "term", "or", "lower", "upper", "p", "n", "events", "adjustment"))
toon("14d. SCP ja vs nee", mm("11_centre_practice", "sacp_any_vs_none.csv"))
toon("14e. SCP ja vs nee: overleving eerste jaar", mm("11_centre_practice", "sacp_any_vs_none_1year_survival.csv"))
toon("15a. Volume-uitkomst: sterfte per volumecategorie", mm("13_volume_outcome", "mortality_by_volume_band.csv"))
toon("15b. Volume-uitkomst: lineair (OR per 10 operaties/jaar)", mm("13_volume_outcome", "volume_linear_models.csv"),
     c("group", "period", "volume", "outcome", "n", "events", "vol_median", "or_per_10", "lower", "upper", "p",
       "or_per_10_random_intercept", "ri_p"))
toon("15c. Volume-uitkomst: spline", mm("13_volume_outcome", "volume_spline_tests.csv"))
toon("15d. Volume-uitkomst: drempel", mm("13_volume_outcome", "volume_threshold_models.csv"),
     c("group", "period", "volume", "outcome", "threshold", "threshold_ci_low", "threshold_ci_high",
       "or_per_10_below_threshold", "lower", "upper", "p_lr_vs_no_volume", "threshold_supported"))
toon("15e. Volume-uitkomst: drempel, cluster-bootstrap (centra)", mm("13_volume_outcome", "volume_threshold_bootstrap.csv"),
     c("group", "period", "outcome", "threshold", "threshold_boot_low", "threshold_boot_high",
       "or_per_10_below_threshold", "or_boot_low", "or_boot_high", "share_reps_hinge_better_than_linear"))
toon("15f. Volume-uitkomst: elbow, vroege uitkomsten", mm("13_volume_outcome", "elbow/elbow_early_outcomes.csv"),
     c("adjustment", "group", "period", "outcome", "centres", "n_units", "p_overall", "p_nonlinear",
       "elbow_last_significant_slope", "elbow_max_curvature"))
toon("15g. Volume-uitkomst: elbow, lange-termijnoverleving", mm("13_volume_outcome", "elbow/elbow_long_term_survival.csv"),
     c("adjustment", "group", "period", "centres", "n_units", "p_overall", "p_nonlinear",
       "elbow_last_significant_slope", "elbow_max_curvature"))

toon("16. Acuut vs niet-acuut (gecorrigeerde OR)", mm("04_outcomes", "acute_vs_nonacute_or.csv"))
toon("17. Acuut x jaar", mm("04_outcomes", "acute_vs_nonacute_year_interaction.csv"))
toon("18. Gestandaardiseerd risico per jaar", mm("04_outcomes", "standardised_risk_by_year.csv"),
     c("outcome", "year", "risk", "lower", "upper"), quote(outcome %in% c("comp_primary", "dzh", "d30")))
toon("19. Variatie tussen centra (MOR/ICC)", mm("04_outcomes", "centre_variation.csv"))
toon("20. FTR-ontleding per jaar", mm("04_outcomes", "ftr_decomposition_by_year.csv"),
     c("definition", "year", "n", "p_complication", "ftr", "mort_no_compl"))
toon("21. Sensitiviteit jaartrend", mm("08_sensitivity", "sensitivity_year_effect.csv"),
     c("analysis", "outcome", "term", "or", "lower", "upper", "p", "n"), quote(outcome %in% c("comp_primary", "dzh", "d30")))
toon("22. Cox-modellen", mm("06_survival", "cox_models.csv"), c("model", "term", "ratio", "ratio_lower", "ratio_upper", "p"),
     quote(term %in% c("year_c", "acute")))
toon("23. Cox stapsgewijs", mm("06_survival", "cox_piecewise_hazard_ratios.csv"), c("term", "ratio", "ratio_lower", "ratio_upper", "p"))
toon("24. Landmark 30 dagen", mm("06_survival", "landmark30_cox.csv"), c("model", "term", "ratio", "ratio_lower", "ratio_upper", "p"),
     quote(term %in% c("year_c", "acute")))
toon("25. Voorwaardelijke overleving", mm("06_survival", "conditional_survival.csv"))
toon("26. RMST", mm("06_survival", "rmst_by_acute_status.csv"))
toon("27. Excess hazard", "07_population/excess_hazard_regression.csv", c("term", "excess_hazard_ratio", "lower", "upper", "p"))
toon("28. Ontbrekende gegevens", "02_descriptives/missing_overall.csv")
toon("29. Code 30/40 per jaar", "02_descriptives/pathology_codes_30_40_by_year.csv")
toon("30. Case-mix trend per jaar", "02_descriptives/characteristics_trend_per_year.csv",
     c("label", "type", "n", "estimate", "lower", "upper", "p", "p_holm"))
toon("31. Sensitiviteit werkwijze", mm("08_sensitivity", "sensitivity_practice.csv"))
toon("32. Convergentie modellen", mm("04_outcomes", "model_convergence_diagnostics.csv"),
     c("outcome", "model", "n_fits", "n_failed", "n_convergence_warning", "n_singular", "n_events", "nobs"))
toon("33. 1-jaarssterfte ruw per groep", mm("15_additional", "one_year_mortality_crude.csv"))
toon("34. 1-jaarssterfte: sensitiviteit jaartrend", mm("15_additional", "sensitivity_one_year_mortality.csv"))
toon("35. Ontbrekend per jaar", mm("15_additional", "missing_by_year.csv"))
toon("36. Ontbrekend per centrum (anoniem)", mm("15_additional", "missing_by_centre_anonymised.csv"))
toon("37. Paper 2: volume aanvullend", mm("15_additional", "volume_additional_analyses.csv"))
toon("38. Paper 3: correctie operatieve complexiteit (OR, 95% BI)", mm("15_additional", "practice_outcome_operative_complexity_wide.csv"))
toon("39. Paper 3: zonder centra met >50% ontbrekende temperatuur", mm("15_additional", "practice_outcome_excluding_centres_missing_temperature.csv"),
     c("group", "practice", "outcome", "adjustment", "centres", "or", "lower", "upper", "p", "n"))
toon("40. Paper 3: stilstand in categorieen (vs <20 min)", mm("15_additional", "arrest_duration_categories.csv"))
toon("41. Paper 3: stilstand spline-toets", mm("15_additional", "arrest_duration_spline_tests.csv"))
toon("42. Paper 3: stilstand x temperatuur", mm("15_additional", "arrest_duration_by_temperature_interaction.csv"))
toon("43. Paper 3: stilstand -> nierfalen -> sterfte", mm("15_additional", "arrest_renal_death_attenuation.csv"))

toon("44. Paper 1: interactie jaar x acuut (BH)", mm("16_paper1_extra", "A_year_by_acute_interaction.csv"),
     c("outcome", "n", "events", "ratio_of_or", "lower", "upper", "p", "p_bh"))
toon("45. Paper 1: jaartrend totaal (Holm)", mm("16_paper1_extra", "B_year_trends_overall_and_by_group.csv"),
     c("outcome", "n", "events", "or_per_year", "lower", "upper", "p", "p_holm"), quote(group == "All"))
toon("46. Paper 1: jaartrend per groep (vanaf 2013 en 2015)", mm("16_paper1_extra", "B_year_trends_overall_and_by_group.csv"),
     c("outcome", "group", "from", "n", "events", "or_per_year", "lower", "upper", "p", "p_bh"), quote(group != "All"))
toon("47. Paper 1: FTR-trend", mm("16_paper1_extra", "B_failure_to_rescue_trend.csv"))
toon("48. Paper 1: acuut vs niet-acuut (OR)", mm("16_paper1_extra", "C_acute_vs_nonacute_or.csv"),
     c("outcome", "n", "events", "or_acute", "lower", "upper", "p"))
toon("49. Paper 1: FTR per groep (incl. % sterfte zonder ernstige complicatie)", mm("16_paper1_extra", "G_ftr_by_group.csv"))
toon("50. Paper 1: volledig vs onvolledig (covariaten)", mm("16_paper1_extra", "H_complete_vs_incomplete.csv"))
toon("51. Paper 1: % volledig per jaar en groep", mm("16_paper1_extra", "H_pct_complete_by_year_and_group.csv"))
toon("52. Paper 1: ontbrekende modelcovariaten", mm("16_paper1_extra", "H_missing_model_covariates.csv"))
toon("53. Paper 1: hazard ratios acuut (stapsgewijs, landmark)", mm("16_paper1_extra", "J_hazard_ratios_acute.csv"))
toon("54. Paper 1: log-rank (dag 0-30, na dag 30)", mm("16_paper1_extra", "J_logrank.csv"))

sink()
message("\nOverzicht ook opgeslagen in: ", out_file)

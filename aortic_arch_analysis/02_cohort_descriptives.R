# =============================================================================
# 02_cohort_descriptives.R
# Builds the analysis cohort (flow diagram), Table 1 by urgency, changes in
# patient characteristics over time, missing-data patterns, and the multiple
# imputations used by the adjusted analyses.
# Output: results/02_descriptives/*, data/derived/{cohort,imputed,specs}.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("02_cohort_descriptives")
SUB <- "02_descriptives"

dv <- load_derived("validated.rds")

# ---- 1. Cohort and flow --------------------------------------------------------------
flow <- data.table(step = "Records in extract", n = nrow(dv))
excl <- dv[included == FALSE, .N, by = .(reason = fifelse(is.na(excl_reden) | trimws(as.character(excl_reden)) == "",
                                                         "(no reason recorded)", as.character(excl_reden)))]
setorder(excl, -N)
flow <- rbind(flow, data.table(step = paste0("  excluded upstream: ", excl$reason), n = -excl$N))
step <- dv[included == TRUE]
flow <- rbind(flow, data.table(step = "Included according to incl_cohort", n = nrow(step)))
n_out_year <- step[is.na(jaar) | !jaar %in% CONFIG$study_years, .N]
step <- step[!is.na(jaar) & jaar %in% CONFIG$study_years]
flow <- rbind(flow, data.table(step = "  excluded: operation year missing or outside study period", n = -n_out_year))
# Explicit surgical selection (re-applied here; see CONFIG): open/hybrid arch surgery with ECC
drop_step <- function(step, cond, label) {
  cond[is.na(cond)] <- FALSE
  flow <<- rbind(flow, data.table(step = paste0("  excluded: ", label), n = -sum(cond)))
  step[!cond]
}
step <- drop_step(step, !(step$ECC %in% c(10, 20, 30)), "no cardiopulmonary bypass (ECC = 0) or ECC unknown")
step <- drop_step(step, step$aotachir_benadering %in% c(20, 21), "endovascular-only aortic repair ((T)EVAR, aotachir_benadering 20/21)")
step <- drop_step(step, step$chir_benadering %in% 50, "percutaneous approach (chir_benadering 50)")
if (!CONFIG$keep_unknown_approach)
  step <- drop_step(step, is.na(step$aotachir_benadering), "aortic surgical approach unknown")
step <- drop_step(step, step$p_chronic %in% TRUE, "chronic dissection (code 30) in any aortic segment")
step <- drop_step(step, step$p_infect %in% TRUE, "infection (code 40) in any aortic segment")
n_unknown_kept <- sum(is.na(step$aotachir_benadering))
n_dup <- step[drop_duplicate == TRUE, .N]
step <- step[drop_duplicate == FALSE]
flow <- rbind(flow, data.table(step = "  excluded: exact duplicate record", n = -n_dup))
flow <- rbind(flow, data.table(step = sprintf("Surgical cohort: open or hybrid arch surgery with ECC (incl. %d with unknown aortic approach)", n_unknown_kept),
                               n = nrow(step)))
write_out(step[, .N, by = .(aotachir_benadering, ECC, chir_benadering)][order(-N)], SUB, "selected_cohort_approach_ecc")

# Two data sets:
#  surgical_cohort - every surgically selected registration (incl. those not classifiable as
#                    acute/non-acute): national volume (total) and Figure 1.
#  cohort          - classified registrations (analysis cohort): all clinical analyses and the
#                    acute/non-acute volume series.
surgical_cohort <- copy(step)
cohort <- drop_step(step, is.na(step$acute), "not classifiable as acute or non-acute (clinical analyses only)")
flow <- rbind(flow, data.table(step = "Analysis cohort (classified as acute or non-acute)", n = nrow(cohort)))
write_out(flow, SUB, "flow_diagram")
log_msg("Surgical cohort: ", nrow(surgical_cohort), "; analysis cohort: ", nrow(cohort))

# Anonymised centre labels (A = largest, based on the surgical cohort). Key stays in data/derived.
cmap <- anonymise_centres(surgical_cohort$centrum)
surgical_cohort[, centre_lab := factor(cmap[centrum], levels = LETTERS[seq_along(cmap)])]
cohort[, centre_lab := factor(cmap[centrum], levels = LETTERS[seq_along(cmap)])]
surgical_cohort[, centrum := factor(centrum)]
cohort[, centrum := factor(centrum)]
saveRDS(data.table(centrum = names(cmap), label = cmap), derived_path("centre_key.rds"))
saveRDS(surgical_cohort, derived_path("surgical_cohort.rds"))

# Spline knots (fixed for all models)
specs <- list(age = spline_spec(cohort$age, 3), egfr = spline_spec(cohort$egfr, 3),
              year = spline_spec(cohort$year_c, 3))
saveRDS(specs, derived_path("specs.rds"))
saveRDS(cohort, derived_path("cohort.rds"))

# ---- 2. Table 1 by urgency --------------------------------------------------------------
t1_rows <- function(d, group) {
  g <- d[[group]] == levels(d[[group]])[2]
  grp_levels <- levels(d[[group]])
  sets <- c(list(Overall = rep(TRUE, nrow(d))),
            setNames(lapply(grp_levels, function(l) d[[group]] == l & !is.na(d[[group]])), grp_levels))
  row_bin <- function(v, label) {
    vals <- vapply(sets, function(s) fmt_npct(d[[v]][s]), character(1))
    data.table(characteristic = label, t(vals), missing = sum(is.na(d[[v]])), smd = smd(d[[v]], g))
  }
  row_cont <- function(v, label, f = fmt_med) {
    vals <- vapply(sets, function(s) f(d[[v]][s]), character(1))
    data.table(characteristic = label, t(vals), missing = sum(is.na(d[[v]])), smd = smd(d[[v]], g))
  }
  row_fac <- function(v, label) {
    f <- d[[v]]
    if (v == "lvef_cat") f <- factor(f, c("<=20%", "21-30%", "31-50%", ">50%"))
    head <- data.table(characteristic = label, t(setNames(rep("", length(sets)), names(sets))),
                       missing = sum(is.na(f)), smd = smd_factor(f, g))
    lv <- rbindlist(lapply(levels(f), function(l) {
      vals <- vapply(sets, function(s) fmt_npct(as.integer(f[s] == l)), character(1))
      data.table(characteristic = paste0("    ", l), t(vals), missing = NA_integer_, smd = NA_real_)
    }))
    rbind(head, lv)
  }
  out <- rbind(
    data.table(characteristic = "N", t(vapply(sets, function(s) as.character(sum(s)), character(1))),
               missing = NA_integer_, smd = NA_real_),
    data.table(characteristic = "Patient characteristics", t(setNames(rep("", length(sets)), names(sets))), missing = NA, smd = NA),
    row_cont("age", "Age, years, mean (SD)", fmt_mean),
    row_bin("female", "Female sex"),
    row_cont("bmi", "BMI, kg/m2, median [IQR]"),
    row_cont("creat", "Creatinine, umol/L, median [IQR]"),
    row_cont("egfr", "eGFR, mL/min/1.73m2, median [IQR]"),
    row_fac("lvef_cat", "LVEF"),
    row_bin("prior_cs", "Previous cardiac surgery"),
    row_bin("endocarditis_b", "Active endocarditis"),
    row_cont("euroscore2", "EuroSCORE II, %, median [IQR]"),
    row_fac("urg4", "Urgency"),
    row_fac("path3", "Aortic pathology"),
    data.table(characteristic = "Operative characteristics", t(setNames(rep("", length(sets)), names(sets))), missing = NA, smd = NA),
    row_bin("sternotomy", "Full sternotomy"),
    row_bin("asc_inv", "Ascending aorta replacement"),
    row_bin("root_inv", "Aortic root involved"),
    row_bin("desc_inv", "Descending aorta involved"),
    row_bin("hybrid_fet", "Combined open + endovascular (incl. frozen elephant trunk)"),
    row_bin("conc_valve", "Concomitant valve surgery"),
    row_bin("aortic_valve", "Concomitant aortic valve surgery"),
    row_bin("conc_cabg", "Concomitant CABG"),
    row_bin("other_cardiac", "Other concomitant cardiac surgery"),
    row_bin("ecc_any", "Cardiopulmonary bypass"),
    row_bin("circ_arrest_b", "Circulatory arrest"),
    row_cont("temp_min", "Minimum core temperature, C, median [IQR]"),
    row_cont("ecc_min", "CPB time, min, median [IQR]"),
    row_cont("aox_min", "Cross-clamp time, min, median [IQR]"),
    row_cont("ca_min", "Circulatory arrest time, min, median [IQR] (if arrest)")
  )
  out
}
cohort_ca <- copy(cohort); cohort_ca[circ_arrest_b != 1, ca_min := NA]
t1 <- t1_rows(cohort_ca[!is.na(urgency2)], "urgency2")
t1[, smd := round(smd, 3)]
write_out(t1, SUB, "table1_baseline_by_urgency")
log_msg("Table 1 written (", sum(is.na(cohort$urgency2)), " patients not classifiable as acute/non-acute: overall column only)")
write_out(cohort[, .N, by = acute_class][order(acute_class)], SUB, "acute_classification_counts")

# Table 1 by era (supplementary)
era_tab <- rbindlist(lapply(levels(cohort$era), function(e) {
  d <- cohort_ca[era == e]
  data.table(era = e, n = nrow(d), age = fmt_mean(d$age), female = fmt_npct(d$female),
             egfr = fmt_med(d$egfr), acute = fmt_npct(d$acute),
             acute_dissection = fmt_npct(d$acute_dissection), prior_cs = fmt_npct(d$prior_cs),
             euroscore2 = fmt_med(d$euroscore2), hybrid_fet = fmt_npct(d$hybrid_fet),
             desc_inv = fmt_npct(d$desc_inv), conc_cabg = fmt_npct(d$conc_cabg),
             conc_valve = fmt_npct(d$conc_valve), root_inv = fmt_npct(d$root_inv))
}))
write_out(era_tab, SUB, "characteristics_by_era")

# ---- 3. Change in characteristics over time ----------------------------------------------
trend_vars <- list(
  bin = c(female = "Female sex", prior_cs = "Previous cardiac surgery",
          endocarditis_b = "Active endocarditis", acute = "Acute (study definition)",
          emerg_salv = "Emergency/salvage", acute_dissection = "Acute aortic dissection",
          conc_cabg = "Concomitant CABG", conc_valve = "Concomitant valve surgery",
          root_inv = "Aortic root involved", desc_inv = "Descending aorta involved",
          hybrid_fet = "Combined open + endovascular", other_cardiac = "Other cardiac surgery",
          sternotomy = "Full sternotomy"),
  cont = c(age = "Age (years)", bmi = "BMI", egfr = "eGFR", lvef = "LVEF (%)",
           euroscore2 = "EuroSCORE II (%)"))
trend_one <- function(v, type, from = min(CONFIG$study_years)) {
  d <- cohort[!is.na(get(v)) & year >= from]
  if (nrow(d) < 30 || uniqueN(d[[v]]) < 2) return(NULL)
  f <- as.formula(paste(v, "~ year_c"))
  fit <- if (type == "bin") glm(f, binomial, data = d) else lm(f, data = d)
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = d$centrum, test = "Satterthwaite")
  r <- ct[2, ]
  est <- r$beta; se <- r$SE; dfs <- r$df_Satt
  lo <- est - qt(.975, dfs) * se; hi <- est + qt(.975, dfs) * se
  data.table(variable = v, label = if (type == "bin") trend_vars$bin[[v]] else trend_vars$cont[[v]],
             type = if (type == "bin") "OR per year" else "Mean change per year",
             n = nrow(d),
             estimate = if (type == "bin") exp(est) else est,
             lower = if (type == "bin") exp(lo) else lo,
             upper = if (type == "bin") exp(hi) else hi,
             df_satterthwaite = dfs, p = r$p_Satt)
}
trends <- rbind(rbindlist(lapply(names(trend_vars$bin), trend_one, type = "bin")),
                rbindlist(lapply(names(trend_vars$cont), trend_one, type = "cont")))
trends[, p_holm := holm(p)]
write_out(trends, SUB, "characteristics_trend_per_year")
# Sensitivity: trends from CONFIG$trend_sens_start onwards
trends_s <- rbind(rbindlist(lapply(names(trend_vars$bin), trend_one, type = "bin", from = CONFIG$trend_sens_start)),
                  rbindlist(lapply(names(trend_vars$cont), trend_one, type = "cont", from = CONFIG$trend_sens_start)))
trends_s[, p_holm := holm(p)]
write_out(trends_s, SUB, paste0("characteristics_trend_per_year_from", CONFIG$trend_sens_start))

# Annual proportions / medians for figures (used in 09)
annual <- cohort[, .(n = .N,
                     age_mean = mean(age, na.rm = TRUE), female = mean(female, na.rm = TRUE),
                     acute = mean(acute, na.rm = TRUE),
                     acute_dissection = mean(acute_dissection, na.rm = TRUE),
                     prior_cs = mean(prior_cs, na.rm = TRUE),
                     euroscore2_median = median(euroscore2, na.rm = TRUE),
                     hybrid_fet = mean(hybrid_fet, na.rm = TRUE), desc_inv = mean(desc_inv, na.rm = TRUE),
                     conc_cabg = mean(conc_cabg, na.rm = TRUE), conc_valve = mean(conc_valve, na.rm = TRUE),
                     root_inv = mean(root_inv, na.rm = TRUE)), by = year][order(year)]
write_out(annual, SUB, "characteristics_by_year")

# ---- 4. Missing data ---------------------------------------------------------------------
miss_vars <- c(age = "Age", female = "Sex", bmi = "BMI", egfr = "eGFR", lvef = "LVEF",
               prior_cs = "Previous cardiac surgery", endocarditis_b = "Endocarditis",
               urg4 = "Urgency", path3 = "Pathology", euroscore2 = "EuroSCORE II",
               conc_valve = "Valve surgery", hybrid_fet = "Open + endovascular",
               circ_arrest_b = "Circulatory arrest", temp_min = "Minimum temperature",
               ca_min = "Arrest duration", d30 = "30-day mortality", dzh = "In-hospital mortality",
               major_compl = "Major complication", vent24 = "Ventilation >24 h",
               icu_readm = "ICU readmission", los = "Length of stay", time_d = "Follow-up time",
               ecmo = "ECMO", pmi = "Postoperative MI", delirium = "Delirium")
mv <- names(miss_vars)
miss_overall <- data.table(variable = mv, label = unname(miss_vars),
                           n_missing = vapply(mv, function(v) sum(is.na(cohort[[v]])), integer(1)),
                           pct_missing = vapply(mv, function(v) round(100 * mean(is.na(cohort[[v]])), 1), numeric(1)))
write_out(miss_overall, SUB, "missing_overall")
miss_year <- melt(cohort[, lapply(.SD, function(x) 100 * mean(is.na(x))), by = year, .SDcols = mv],
                  id.vars = "year", variable.name = "variable", value.name = "pct")
miss_centre <- melt(cohort[, lapply(.SD, function(x) 100 * mean(is.na(x))), by = centre_lab, .SDcols = mv],
                    id.vars = "centre_lab", variable.name = "variable", value.name = "pct")
write_out(dcast(miss_year, variable ~ year, value.var = "pct"), SUB, "missing_by_year")
write_out(dcast(miss_centre, variable ~ centre_lab, value.var = "pct"), SUB, "missing_by_centre")
heat <- function(d, x, title) {
  d <- copy(d); d[, label := factor(miss_vars[as.character(variable)], rev(unname(miss_vars)))]
  ggplot(d, aes(x = factor(.data[[x]]), y = label, fill = pct)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(pct >= 1, sprintf("%.0f", pct), "")), size = 2.3) +
    scale_fill_gradient(low = "#f7fbff", high = "#d1495b", limits = c(0, 100), name = "% missing") +
    labs(x = NULL, y = NULL, title = title) + theme_pub() +
    theme(axis.text.x = element_text(angle = 0), panel.grid = element_blank())
}
save_fig(heat(miss_year, "year", "Missing data by operation year"), SUB, "missing_by_year", 8, 6.5)
save_fig(heat(miss_centre, "centre_lab", "Missing data by centre (anonymised)"), SUB, "missing_by_centre", 8, 6.5)

# ---- 5. Multiple imputation -------------------------------------------------------------
log_msg("Multiple imputation: m = ", CONFIG$n_imputations)
imps <- make_imputations(cohort, m = CONFIG$n_imputations)
mids <- attr(imps, "mids")
write_out(if (is.null(mids$loggedEvents)) data.table(note = "no logged events") else mids$loggedEvents,
          SUB, "imputation_logged_events")
write_out(data.table(variable = names(mids$method), method = mids$method,
                     n_missing = vapply(names(mids$method), function(v) sum(is.na(mids$data[[v]])), integer(1))),
          SUB, "imputation_methods")
# Diagnostics: observed vs imputed distributions; convergence traces
dens <- rbindlist(lapply(c("egfr", "lvef", "bmi", "euroscore2"), function(v) {
  if (sum(is.na(cohort[[v]])) == 0) return(NULL)
  rbind(data.table(variable = v, type = "Observed", value = cohort[[v]][!is.na(cohort[[v]])]),
        data.table(variable = v, type = "Imputed",
                   value = unlist(lapply(imps, function(d) d[[v]][is.na(cohort[[v]])]))))
}))
if (nrow(dens)) {
  p <- ggplot(dens, aes(value, colour = type)) + geom_density() +
    facet_wrap(~variable, scales = "free") + scale_colour_manual(values = PAL[1:2], name = NULL) +
    labs(title = "Observed vs imputed values", x = NULL, y = "Density") + theme_pub()
  save_fig(p, SUB, "imputation_density", 8, 5)
}
grDevices::png(file.path(out_dir(SUB), "imputation_convergence.png"), width = 1800, height = 1400, res = 200)
print(plot(mids))
grDevices::dev.off()
attr(imps, "mids") <- NULL
saveRDS(imps, derived_path("imputed.rds"))
saveRDS(mids, derived_path("mids.rds"))
log_msg("Done.")

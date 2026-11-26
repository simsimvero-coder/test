# =============================================================================
# sens_2018.R  (standalone; run after run_all.R)
# Sensitivity analysis restricted to operations from 2018 onwards, when
# registration was substantially more complete (EuroSCORE II missing < 9%,
# aortic pathology complete).
#   A. Adjusted odds ratio per calendar year, 12 outcomes, all patients and per
#      group (same model as the primary analysis: M1, or the reduced set when
#      events per variable < 10).
#   B. Acute versus non-acute, from 2018, with and without EuroSCORE II in the
#      same complete-case sample.
#   C. The procedures not classifiable as acute or non-acute (2013-2014): early
#      mortality compared with classified procedures, and the national year trend
#      with and without them (M0 year + centre, and M1).
# Writes only aggregate numbers. Output: results/<method>/sens_2018/
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages(library(data.table))
SUB <- "sens_2018"
options(width = 220)
log_msg("Start sens_2018")

cohort <- load_derived("cohort.rds")
specs  <- load_derived("specs.rds")
cc     <- analysis_datasets()
START  <- 2018L
D1_END <- CONFIG$d1y_last_year
OUTS <- c(comp_primary = "Primary composite outcome", dzh = "In-hospital death", d30 = "30-day death",
          d1y = "One-year death", stroke_def = "Stroke with permanent deficit",
          stroke_nodef = "Stroke or TIA without deficit", renal = "Renal failure", gi = "GI complication",
          vent24 = "Ventilation >24 h", reintub = "Reintubation", icu_readm = "ICU readmission",
          major_compl = "Major complication")
OUTS <- OUTS[names(OUTS) %in% names(cohort)]
end_of <- function(v) if (v == "d1y") D1_END else max(CONFIG$study_years)
safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })
sub_d <- function(d, v, g) { d <- d[year >= START & year <= end_of(v) & !is.na(get(v))]; if (is.null(g)) d else d[acute == g] }
row_of <- function(fit, term) {
  if (is.null(fit)) return(data.table(or = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  tt <- term; p <- fit$pooled[fit$pooled$term == tt]
  data.table(or = exp(p$estimate), lower = exp(p$lower), upper = exp(p$upper), p = p$p)
}

# ---- A. Year trend from 2018 ----------------------------------------------------------
A <- rbindlist(lapply(names(OUTS), function(v) rbindlist(lapply(list(NULL, 1L, 0L), function(g) {
  f <- build_formula(v, COVARS_M1, specs); model <- "M1"
  dl <- complete_for(lapply(cc, sub_d, v = v, g = g), f)
  y <- dl[[1]][[v]]
  if (sum(y == 1) > 0 && epv(y, n_fixed_par(f, dl[[1]])) < 10) {
    f <- build_formula(v, COVARS_REDUCED, specs); model <- "reduced"
    dl <- complete_for(lapply(cc, sub_d, v = v, g = g), f); y <- dl[[1]][[v]]
  }
  fit <- if (sum(y == 1) >= 10) safe(fit_pooled(f, dl), paste(v, g)) else NULL
  cbind(data.table(outcome = v, label = OUTS[[v]], group = if (is.null(g)) "All" else if (g == 1) "Acute" else "Non-acute",
                   from = START, to = end_of(v), model = model, n = nrow(dl[[1]]), events = sum(y == 1)),
        row_of(fit, "year_c"))
}))))
sec <- A$group == "All" & A$outcome != "comp_primary"
A[sec, p_holm := p.adjust(p, "holm")]
write_out(A, SUB, "A_year_trend_from_2018")

# ---- B. Acute vs non-acute from 2018, with and without EuroSCORE II -------------------------------
B <- rbindlist(lapply(names(OUTS), function(v) {
  f_m1 <- build_formula(v, c(COVARS_M1, "acute"), specs)
  f_no <- build_formula(v, c("age", "female", "prior_cs", "acute"), specs)
  dl <- complete_for(lapply(cc, sub_d, v = v, g = NULL), f_m1)          # same sample for both models
  rbind(cbind(outcome = v, label = OUTS[[v]], model = "With EuroSCORE II", n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1),
              row_of(safe(fit_pooled(f_m1, dl), v), "acute")),
        cbind(outcome = v, label = OUTS[[v]], model = "Without EuroSCORE II (same sample)", n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1),
              row_of(safe(fit_pooled(f_no, dl), v), "acute")))
}))
write_out(B, SUB, "B_acute_vs_nonacute_from_2018")

# ---- C. The 116 procedures not classifiable as acute or non-acute ------------------------------
# Surgical cohort = analysis cohort + unclassifiable procedures (2013-2014).
sc <- load_derived("surgical_cohort.rds")
sc[, centrum := as.character(centrum)]
sc[, klasse := fifelse(acute %in% 1, "Acute", fifelse(acute %in% 0, "Non-acute", "Not classifiable"))]
pct <- function(x) round(100 * mean(x == 1, na.rm = TRUE), 1)
C1 <- sc[year %in% 2013:2014, .(n = .N,
  in_hospital_pct = pct(dzh), d30_pct = pct(d30), composite_pct = pct(comp_primary), one_year_pct = pct(d1y),
  n_known_in_hospital = sum(!is.na(dzh)), euroscore2_median = round(median(euroscore2, na.rm = TRUE), 1),
  euroscore2_missing_pct = round(100 * mean(is.na(euroscore2)), 1), age_median = median(age, na.rm = TRUE)),
  by = .(year, klasse)][order(year, klasse)]
C2 <- rbind(sc[year %in% 2013:2014, .(set = "With unclassifiable", n = .N, in_hospital_pct = pct(dzh), d30_pct = pct(d30)), by = year],
            sc[year %in% 2013:2014 & klasse != "Not classifiable", .(set = "Analysis cohort", n = .N, in_hospital_pct = pct(dzh), d30_pct = pct(d30)), by = year])[order(year, set)]
sc_d <- rederive(sc)
trend_c <- function(v, covs, lab) {
  f <- build_formula(v, covs, specs)
  rbindlist(lapply(c("Analysis cohort", "With unclassifiable"), function(s) {
    d <- sc_d[!is.na(get(v)) & year <= end_of(v)]
    if (s == "Analysis cohort") d <- d[klasse != "Not classifiable"]
    dl <- complete_for(list(d), f)
    cbind(data.table(outcome = v, model = lab, set = s, n = nrow(dl[[1]]), events = sum(dl[[1]][[v]] == 1),
                     unclassifiable_in_model = sum(dl[[1]]$klasse == "Not classifiable")),
          row_of(safe(fit_pooled(f, dl), paste(v, s)), "year_c"))
  }))
}
C3 <- rbindlist(lapply(intersect(c("dzh", "d30", "comp_primary", "d1y"), names(OUTS)), function(v)
  rbind(trend_c(v, character(0), "M0: year + centre"), trend_c(v, COVARS_M1, "M1"))))
write_out(C1, SUB, "C1_unclassifiable_2013_2014"); write_out(C2, SUB, "C2_crude_mortality_2013_2014")
write_out(C3, SUB, "C3_year_trend_with_unclassifiable")

fmt <- function(d) d[, `:=`(OR = sprintf("%.3f (%.3f-%.3f)", or, lower, upper), p = signif(p, 3))][]
out <- c(sprintf("Sensitivity analysis from %d - %s", START, format(Sys.time(), "%Y-%m-%d %H:%M")),
         "", "==== A. OR per jaar vanaf 2018 ====",
         capture.output(print(fmt(copy(A))[, .(outcome, group, from, to, model, n, events, OR, p, p_holm = signif(p_holm, 3))], nrows = 100)),
         "", "==== B. Acuut vs niet-acuut vanaf 2018, met en zonder EuroSCORE II (zelfde cohort) ====",
         capture.output(print(fmt(copy(B))[, .(outcome, model, n, events, OR, p)], nrows = 100)),
         "", "==== C1. Niet-classificeerbaar vs acuut vs niet-acuut, 2013-2014 ====", capture.output(print(C1)),
         "", "==== C2. Ruwe sterfte 2013-2014 met en zonder niet-classificeerbaar ====", capture.output(print(C2)),
         "", "==== C3. OR per jaar met en zonder de niet-classificeerbare operaties ====",
         capture.output(print(fmt(copy(C3))[, .(outcome, model, set, n, events, unclassifiable_in_model, OR, p)], nrows = 100)))
f_out <- file.path(out_dir(SUB), "overzicht_2018.txt")
writeLines(out, f_out); cat(out, sep = "\n")
log_msg("Done. Overview in ", f_out)

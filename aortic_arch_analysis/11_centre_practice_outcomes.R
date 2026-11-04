# =============================================================================
# 11_centre_practice_outcomes.R
# Outcomes per centre related to organ-protection practice (minimum temperature,
# use of circulatory arrest, cross-clamp time, circulatory arrest time), and a
# subgroup analysis of selective cerebral perfusion (acute surgery, years with data).
#
# A. Centre profile: practice and crude / case-mix-adjusted outcomes per centre.
# B. Centre-level (ecological) association: Spearman correlation and weighted
#    regression of centre log(O/E) on centre practice (15 centres: exploratory).
# C. Patient-level within/between-centre ("hybrid", Mundlak) mixed models: each
#    practice variable split into the centre mean (between) and the patient's
#    deviation from it (within), adjusted for the primary covariates and year,
#    with a centre random intercept. Sensitivity: additionally adjusted for acute.
# D. Selective cerebral perfusion (cerebrale_perfusie_aanv / _duur), acute
#    patients with a recorded value: descriptives, adjusted comparison of
#    perfusion type, and perfusion duration.
# Cross-clamp and arrest durations partly reflect procedural complexity
# (confounding by indication); all results here are exploratory.
# Output: results/<method>/11_centre_practice/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("11_centre_practice_outcomes")
log_msg("Missing-data method: ", method_label())
SUB <- "11_centre_practice"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()
specs  <- load_derived("specs.rds")

OUTS <- c(comp_primary = "Primary composite", dzh = "In-hospital mortality", d30 = "30-day mortality",
          stroke_def = "Stroke with permanent deficit", renal = "Renal failure",
          major_compl = "Major complication")
OUTS <- OUTS[names(OUTS) %in% names(cohort)]

safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })

# ---- A. Centre profile ---------------------------------------------------------------------
prof <- cohort[, .(n = .N, pct_acute = 100 * mean(urgency2 == "Acute", na.rm = TRUE),
                   temp_median = as.numeric(median(temp_min, na.rm = TRUE)),
                   temp_missing_pct = 100 * mean(is.na(temp_min)),
                   pct_deep = 100 * mean(temp_min <= 20, na.rm = TRUE),
                   pct_arrest = 100 * mean(circ_arrest_b == 1, na.rm = TRUE),
                   aox_median = as.numeric(median(aox_min, na.rm = TRUE)),
                   ca_median = as.numeric(median(ca_min[circ_arrest_b == 1 & ca_min > 0], na.rm = TRUE))),
               by = centre_lab][order(centre_lab)]

oe_rows <- list()
for (v in names(OUTS)) {
  log_msg("Centre O/E: ", v)
  d <- complete_for(list(cohort[!is.na(get(v))]), build_formula(v, COVARS_M1, specs))[[1]]
  f <- safe(fit_mixed(build_formula(v, COVARS_M1, specs), d, binomial()), v)
  if (is.null(f) || is.null(f$fit)) next
  d[, expected := predict(f$fit, newdata = d, re.form = NA, type = "response")]
  re <- as.data.frame(ranef(f$fit, condVar = TRUE)); setDT(re)
  key <- unique(cohort[, .(centrum = as.character(centrum), centre_lab)])
  re <- merge(re[, .(centrum = as.character(grp), re = condval, re_se = condsd)], key, by = "centrum")
  o <- d[, .(n = .N, observed = sum(get(v) == 1), expected = sum(expected)), by = centre_lab]
  o <- merge(o, re[, .(centre_lab, re, re_se)], by = "centre_lab")
  o[, `:=`(outcome = v, crude_pct = 100 * observed / n, oe = observed / expected,
           centre_or = exp(re), centre_or_lower = exp(re - 1.96 * re_se), centre_or_upper = exp(re + 1.96 * re_se))]
  oe_rows[[v]] <- o
}
oe <- rbindlist(oe_rows)
write_out(prof, SUB, "centre_practice_profile")
write_out(oe[, .(outcome, centre_lab, n, observed, expected, crude_pct, oe, centre_or, centre_or_lower, centre_or_upper)],
          SUB, "centre_outcomes_observed_expected")
wide <- dcast(oe, centre_lab ~ outcome, value.var = "oe")
setnames(wide, setdiff(names(wide), "centre_lab"), paste0("oe_", setdiff(names(wide), "centre_lab")))
write_out(merge(prof, wide, by = "centre_lab"), SUB, "centre_profile_practice_and_outcomes")

# ---- B. Centre-level (ecological) association ------------------------------------------------
PRACT <- c(temp_median = "Median minimum temperature (C)", pct_arrest = "Circulatory arrest use (%)",
           aox_median = "Median cross-clamp time (min)", ca_median = "Median circulatory arrest time (min)",
           pct_deep = "Deep hypothermia <=20 C (%)")
eco <- rbindlist(lapply(names(OUTS), function(v) {
  o <- merge(oe[outcome == v], prof, by = "centre_lab", suffixes = c("", ".prof"))
  rbindlist(lapply(names(PRACT), function(x) {
    ok <- is.finite(o[[x]]) & o$observed > 0
    if (sum(ok) < 6) return(NULL)
    ct <- suppressWarnings(cor.test(o[[x]][ok], o$oe[ok], method = "spearman", exact = FALSE))
    z <- (o[[x]][ok] - mean(o[[x]][ok])) / sd(o[[x]][ok])
    w <- lm(log(o$oe[ok]) ~ z, weights = o$expected[ok])
    cf <- summary(w)$coefficients
    data.table(outcome = v, practice = PRACT[[x]], n_centres = sum(ok), spearman_rho = unname(ct$estimate),
               p_spearman = ct$p.value, oe_ratio_per_sd = exp(cf[2, 1]),
               lower = exp(cf[2, 1] - qt(.975, w$df.residual) * cf[2, 2]),
               upper = exp(cf[2, 1] + qt(.975, w$df.residual) * cf[2, 2]), p_weighted_regression = cf[2, 4])
  }))
}))
write_out(eco, SUB, "centre_level_practice_vs_outcome")

pd <- merge(oe[outcome %in% c("comp_primary", "dzh", "stroke_def")], prof, by = "centre_lab", suffixes = c("", ".prof"))
pd <- melt(pd, id.vars = c("centre_lab", "outcome", "oe", "expected"),
           measure.vars = c("temp_median", "pct_arrest", "aox_median", "ca_median"), variable.name = "practice", value.name = "value")
pd[, practice := factor(PRACT[as.character(practice)], PRACT[c("temp_median", "pct_arrest", "aox_median", "ca_median")])]
pd[, outcome := factor(OUTS[outcome], OUTS[c("comp_primary", "dzh", "stroke_def")])]
p <- ggplot(pd[is.finite(value)], aes(value, oe)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(aes(size = expected), colour = PAL[1], alpha = .7) + geom_text(aes(label = centre_lab), size = 2.4, vjust = -1) +
  facet_grid(outcome ~ practice, scales = "free") + scale_size_area(max_size = 4, guide = "none") +
  labs(x = NULL, y = "Observed / expected", title = "Centre practice and case-mix-adjusted outcome (each point = centre)") + theme_pub()
save_fig(p, SUB, "centre_practice_vs_oe", 11, 8)

# ---- C. Patient-level within/between-centre models --------------------------------------------
add_wb <- function(d, x, scale_w, scale_b) {
  d <- copy(d)
  d[, xb := mean(get(x), na.rm = TRUE), by = centrum]
  d[, `:=`(x_within = (get(x) - xb) / scale_w, x_between = xb / scale_b)]
  d
}
XS <- list(
  list(x = "temp_min", sw = 1, sb = 1, label = "Minimum temperature (per 1 C)", sub = function(d) d[!is.na(temp_min)]),
  list(x = "circ_arrest_b", sw = 1, sb = 0.1, label = "Circulatory arrest (within: yes vs no; between: per 10% points)",
       sub = function(d) d[!is.na(circ_arrest_b)]),
  list(x = "aox_min", sw = 30, sb = 30, label = "Cross-clamp time (per 30 min)", sub = function(d) d[!is.na(aox_min)]),
  list(x = "ca_min", sw = 10, sb = 10, label = "Circulatory arrest time (per 10 min; arrest patients)",
       sub = function(d) d[circ_arrest_b == 1 & ca_min > 0]))
hyb <- list()
for (xs in XS) for (v in names(OUTS)) for (adj in c("M1", "M1 + acute")) {
  log_msg("Within/between: ", xs$x, " -> ", v, " [", adj, "]")
  covs <- if (adj == "M1") COVARS_M1 else c(COVARS_M1, "acute")
  f <- as.formula(paste(v, "~ x_within + x_between +", paste(c("year_c", covar_terms(covs, specs), "(1 | centrum)"), collapse = " + ")),
                  env = globalenv())
  dl <- lapply(imps, function(d) add_wb(xs$sub(d[!is.na(get(v))]), xs$x, xs$sw, xs$sb))
  dl <- complete_for(dl, f)
  ev <- sum(dl[[1]][[v]] == 1)
  if (ev < 50) { log_msg("  too few events (", ev, ")"); next }
  r <- safe(fit_pooled(f, dl), paste(xs$x, v))
  if (is.null(r)) next
  e <- exp_table(r$pooled, c("x_within", "x_between"))
  hyb[[length(hyb) + 1]] <- e[, .(outcome = v, practice = xs$label, adjustment = adj,
                                  component = fifelse(term == "x_within", "within centre", "between centres"),
                                  or = ratio, lower = ratio_lower, upper = ratio_upper, p, n = nrow(dl[[1]]), events = ev)]
}
hyb <- rbindlist(hyb)
write_out(hyb, SUB, "within_between_centre_practice_models")

# ---- D. Selective cerebral perfusion -------------------------------------------------------
sc <- cohort[urgency2 == "Acute" & !is.na(sacp_cat)]
if (nrow(sc) < 30) {
  log_msg("Cerebral perfusion: ", nrow(sc), " acute patients with a recorded value - subgroup analysis skipped",
          if (all(is.na(cohort$cerebrale_perfusie_aanv))) " (columns cerebrale_perfusie_aanv/_duur absent or empty)" else "")
  write_out(data.table(note = "cerebral perfusion data absent or too few records", n = nrow(sc)), SUB, "sacp_not_analysed")
} else {
  log_msg("Cerebral perfusion: ", nrow(sc), " acute patients, years ", paste(range(sc$year), collapse = "-"))
  write_out(sc[, .N, by = .(year, sacp_cat)][order(year, sacp_cat)], SUB, "sacp_type_by_year")
  write_out(sc[, .N, by = .(centre_lab, sacp_group)][order(centre_lab, sacp_group)], SUB, "sacp_type_by_centre")
  desc <- sc[, c(list(n = .N, temp_median = fmt_med(temp_min), arrest = fmt_npct(circ_arrest_b),
                      ca_median = fmt_med(ca_min[circ_arrest_b == 1]), sacp_min_median = fmt_med(sacp_min)),
                 lapply(.SD, fmt_npct)), by = sacp_cat, .SDcols = names(OUTS)][order(sacp_cat)]
  write_out(desc, SUB, "sacp_descriptives_by_type")
  sacp_rows <- list()
  run_sacp <- function(label, data, term_rhs, keep) {
    for (v in names(OUTS)) {
      d0 <- data[!is.na(get(v))]
      for (covs in list(COVARS_M1, COVARS_REDUCED)) {
        f <- as.formula(paste(v, "~", term_rhs, "+", paste(c("year_c", covar_terms(covs, specs), "(1 | centrum)"), collapse = " + ")),
                        env = globalenv())
        dl <- complete_for(lapply(imps, function(d) d[urgency2 == "Acute" & bronrij %in% d0$bronrij]), f)
        ev <- sum(dl[[1]][[v]] == 1)
        if (epv(dl[[1]][[v]], n_fixed_par(f, dl[[1]])) >= 10 || identical(covs, COVARS_REDUCED)) break
      }
      if (ev < 20) { log_msg("  ", label, " / ", v, ": too few events (", ev, ")"); next }
      r <- safe(fit_pooled(f, dl), paste(label, v))
      if (is.null(r)) next
      e <- exp_table(r$pooled)[grepl(keep, term)]
      sacp_rows[[length(sacp_rows) + 1]] <<- e[, .(analysis = label, outcome = v, term, or = ratio, lower = ratio_lower,
                                                   upper = ratio_upper, p, n = nrow(dl[[1]]), events = ev,
                                                   adjustment = if (identical(covs, COVARS_M1)) "M1" else "reduced (EPV < 10)")]
    }
  }
  grp_n <- sc[, .N, by = sacp_group]
  ref <- as.character(grp_n[which.max(N), sacp_group])
  keep_groups <- as.character(grp_n[N >= 10, sacp_group])
  for (d in seq_along(imps)) {
    imps[[d]][, sacp_group_m := relevel(droplevels(factor(as.character(sacp_group), levels = keep_groups)), ref)]
    imps[[d]][, antegrade_bilateral := fifelse(sacp_cat == "Antegrade bilateral", 1L,
                                               fifelse(sacp_cat == "Antegrade unilateral", 0L, NA_integer_))]
    imps[[d]][, sacp_per10 := sacp_min / 10]
  }
  if (length(keep_groups) >= 2)
    run_sacp(paste0("Perfusion type (reference: ", ref, ")"), sc[sacp_group %in% keep_groups], "sacp_group_m", "^sacp_group_m")
  if (sc[sacp_cat %in% c("Antegrade unilateral", "Antegrade bilateral"), min(table(droplevels(sacp_cat)))] >= 10)
    run_sacp("Antegrade bilateral vs unilateral", sc[sacp_cat %in% c("Antegrade unilateral", "Antegrade bilateral")],
             "antegrade_bilateral", "^antegrade_bilateral$")
  dur <- sc[sacp_group == "Antegrade" & !is.na(sacp_min) & sacp_min > 0]
  if (nrow(dur) >= 50) {
    run_sacp("Antegrade perfusion duration (per 10 min)", dur, "sacp_per10", "^sacp_per10$")
    run_sacp("Antegrade perfusion duration (per 10 min), adjusted for minimum temperature", dur[!is.na(temp_min)],
             "sacp_per10 + temp_min", "^sacp_per10$|^temp_min$")
  }
  write_out(rbindlist(sacp_rows), SUB, "sacp_models")
  # Any selective cerebral perfusion vs none: crude (logistic) and adjusted (M1 + year + centre)
  if (sc[sacp_group == "None", .N] >= 10 && sc[sacp_group != "None", .N] >= 10) {
    for (d in seq_along(imps)) imps[[d]][, scp_any := fifelse(is.na(sacp_group), NA_integer_, as.integer(sacp_group != "None"))]
    any_rows <- rbindlist(lapply(c(names(OUTS), "stroke_any", "gi"), function(v) {
      if (!v %in% names(cohort)) return(NULL)
      x <- sc[!is.na(get(v))][, scp_any := as.integer(sacp_group != "None")]
      g <- summary(glm(as.formula(paste(v, "~ scp_any")), binomial, x))$coefficients["scp_any", ]
      f <- build_formula(v, c("scp_any", COVARS_M1), specs)
      dl <- complete_for(lapply(imps, function(d) d[urgency2 == "Acute" & bronrij %in% x$bronrij]), f)
      a <- safe(exp_table(fit_pooled(f, dl)$pooled, "scp_any"), paste("scp_any", v))
      data.table(outcome = v, scp_events = sum(x[scp_any == 1][[v]] == 1), scp_n = sum(x$scp_any == 1),
                 none_events = sum(x[scp_any == 0][[v]] == 1), none_n = sum(x$scp_any == 0),
                 or_crude = exp(g[1]), lower_crude = exp(g[1] - 1.96 * g[2]), upper_crude = exp(g[1] + 1.96 * g[2]), p_crude = g[4],
                 or_adjusted = if (is.null(a)) NA_real_ else a$ratio, lower_adjusted = if (is.null(a)) NA_real_ else a$ratio_lower,
                 upper_adjusted = if (is.null(a)) NA_real_ else a$ratio_upper, p_adjusted = if (is.null(a)) NA_real_ else a$p,
                 n_adjusted = nrow(dl[[1]]))
    }))
    write_out(any_rows, SUB, "sacp_any_vs_none")
  }
}

log_msg("Done.")

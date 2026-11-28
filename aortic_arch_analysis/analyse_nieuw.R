# =============================================================================
# analyse_nieuw.R - ONE file with all new analyses (no other new scripts needed)
#   Part 1  perfusion practice and volume per centre       -> results/centre_perfusion/
#   Part 2  paper 1, sensitivity: acute vs strictly elective -> results/sens_strict_groups/
#   Part 3  paper 2, volume-outcome, strict groups + Davies  -> results/volume/
# Needs the project folder (R/, data/derived from run_all.R; run_all.R need not be rerun).
# Use in RGui:  source(file.choose())   and choose this file.
# =============================================================================

# 1. Project folder: use the current folder if it is the project, otherwise ask for it
if (!file.exists(file.path("R", "00_config.R"))) {
  map <- if (.Platform$OS.type == "windows") utils::choose.dir(caption = "Kies de map aortic_arch_analysis") else
    readline("Pad naar de map aortic_arch_analysis: ")
  if (is.na(map) || !file.exists(file.path(map, "R", "00_config.R")))
    stop("Dit is niet de projectmap: R/00_config.R ontbreekt in ", map)
  setwd(map)
}
cat("Projectmap:", getwd(), "\n")

# 2. Settings (complete case is the primary analysis)
Sys.setenv(ARCH_METHODS = "complete_case", ARCH_MISSING = "complete_case", ARCH_VOLUME_BOOT = "500")

# 3. Data. The scripts do not read the Excel file themselves: they read the cleaned
#    cohort (data/derived/cohort.rds) that 01_data_validation.R and
#    02_cohort_descriptives.R make from the raw file named in R/00_config.R
#    ("AORTABOOG_NL_schoon v3(1).xlsx", found on the desktop, or the path in
#    Sys.setenv(ARCH_RAW_DATA = "...")). If that cohort is missing, it is made now.
source(file.path("R", "00_config.R"))
cat("Ruwe dataset volgens 00_config.R:", CONFIG$raw_data_path, "\n")
if (!file.exists(file.path(CONFIG$derived_dir, "cohort.rds"))) {
  cat("data/derived/cohort.rds ontbreekt: 01 en 02 worden nu gedraaid (leest de Excel in)\n")
  source("01_data_validation.R")
  source("02_cohort_descriptives.R")
}
cat("Gebruikt cohort: data/derived/cohort.rds, gemaakt op",
    format(file.mtime(file.path(CONFIG$derived_dir, "cohort.rds")), "%Y-%m-%d %H:%M"), "\n")

# 4. Package for the Davies test (volume.R skips only that test if it cannot be installed)
if (!requireNamespace("segmented", quietly = TRUE))
  try(install.packages("segmented"), silent = TRUE)

.status_nieuw <- c()
while (sink.number() > 0) sink()

# #############################################################################
# Part 1: perfusion practice per centre  (= centre_perfusion.R)
# #############################################################################
cat("\n\n#############  Part 1: perfusion practice per centre  #############\n"); .t0 <- Sys.time()
.status_nieuw["Part 1: perfusion practice per centre"] <- tryCatch({
# =============================================================================
# centre_perfusion.R  (standalone; run after run_all.R)
# Perfusion practice per centre, including centre volume, and comparison between
# centres. Centres are shown with their anonymised letter (centre_lab) and are
# ordered by volume only (never by practice or outcome); no ranking is made.
# Groups: All, Acute (urgentie 30/40 or acute dissection code) and Elective
# (urgentie 10 without dissection code), as in volume.R. Periods: 2013-2024 and
# 2018-2024 (complete registration).
#   A. Centre profile: number of operations, mean annual volume, % acute, and per
#      practice variable the centre median [IQR] or %, with the number recorded.
#   B. Between-centre comparison per practice variable:
#        - spread of centre values (minimum, quartiles, maximum; centres with >= 10
#          recorded values);
#        - heterogeneity test across centres (Kruskal-Wallis / chi-square);
#        - variation attributable to centre: ICC (and median odds ratio for binary
#          variables) from random-intercept models, crude and adjusted for case mix
#          (age, sex, logit EuroSCORE II, previous cardiac surgery, calendar year;
#          group All also for acute presentation);
#        - relation with centre volume: Spearman correlation across centres, and a
#          patient-level regression per 10 operations/year with CR2 standard errors
#          clustered by centre (Satterthwaite).
#   C. Funnel-plot data for binary practice variables (centre % against the
#      national %, with 95% and 99.8% control limits).
#   Figure: centre values per practice variable, centres ordered by volume.
# Writes only aggregate numbers. Output: results/<method>/centre_perfusion/
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(ggplot2); library(splines) })
SUB <- "centre_perfusion"
options(width = 250)
log_msg("Start centre_perfusion")
MIN_REC <- 10

specs <- load_derived("specs.rds")
co <- copy(load_derived("cohort.rds"))
co[, centrum := as.character(centrum)]
co[, grp := fcase(acute %in% 1L, "Acute", urgentie %in% 10 & !(seg_ad %in% TRUE), "Elective", default = NA_character_)]
co[, `:=`(
  deep        = fifelse(is.na(temp_min), NA_integer_, as.integer(temp_min <= 20)),
  moderate    = fifelse(is.na(temp_min), NA_integer_, as.integer(temp_min > 20 & temp_min <= 28)),
  mild_normo  = fifelse(is.na(temp_min), NA_integer_, as.integer(temp_min > 28)),
  ca_arrest   = fifelse(circ_arrest_b %in% 1 & ca_min > 0, ca_min, NA_real_),
  scp_ante    = fifelse(is.na(sacp_group), NA_integer_, as.integer(sacp_group == "Antegrade")),
  scp_bilat   = fifelse(sacp_cat %in% c("Antegrade unilateral", "Antegrade bilateral"),
                        as.integer(sacp_cat == "Antegrade bilateral"), NA_integer_),
  scp_min     = fifelse(sacp_group %in% "Antegrade", sacp_min, NA_real_))]
CONT <- c(temp_min = "Nadir temperature (C)", ca_arrest = "Circulatory arrest time (min, arrest patients)",
          scp_min = "Antegrade cerebral perfusion time (min)", ecc_min = "Cardiopulmonary bypass time (min)",
          aox_min = "Cross-clamp time (min)")
BIN  <- c(deep = "Deep hypothermia (<=20 C)", moderate = "Moderate hypothermia (20.1-28 C)",
          mild_normo = "Mild hypothermia or normothermia (>28 C)", circ_arrest_b = "Circulatory arrest",
          scp_ante = "Antegrade cerebral perfusion (of recorded)", scp_bilat = "Bilateral (of antegrade with side recorded)",
          hybrid_fet = "Combined open + endovascular (incl. FET)", root_inv = "Aortic root involved")
CONT <- CONT[names(CONT) %in% names(co)]; BIN <- BIN[names(BIN) %in% names(co)]
LAB  <- c(CONT, BIN)
safe <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })
PERIODS <- list("2013-2024" = range(CONFIG$study_years), "2018-2024" = c(2018L, max(CONFIG$study_years)))
GROUPS  <- c("All", "Acute", "Elective")
sel <- function(per, g) { r <- PERIODS[[per]]; d <- co[year >= r[1] & year <= r[2]]; if (g == "All") d else d[grp %in% g] }

prof_rows <- list(); cmp_rows <- list(); fun_rows <- list(); fig_rows <- list()
for (per in names(PERIODS)) for (g in GROUPS) {
  d <- sel(per, g); r <- PERIODS[[per]]; ny <- r[2] - r[1] + 1
  log_msg(per, " / ", g, ": ", nrow(d), " operations")
  cen <- d[, .(n = .N, vol_mean = round(.N / ny, 1), pct_acute = round(100 * mean(acute %in% 1), 1)), by = centre_lab]
  d <- merge(d, cen[, .(centre_lab, vol_mean)], by = "centre_lab")
  # ---- A. profile ----
  for (v in names(LAB)) {
    x <- d[!is.na(get(v)), .(rec = .N, val = if (v %in% names(CONT)) as.numeric(median(get(v))) else 100 * mean(get(v) == 1),
                             q1 = if (v %in% names(CONT)) as.numeric(quantile(get(v), .25)) else NA_real_,
                             q3 = if (v %in% names(CONT)) as.numeric(quantile(get(v), .75)) else NA_real_), by = centre_lab]
    x <- merge(cen, x, by = "centre_lab", all.x = TRUE)
    x[is.na(rec), rec := 0L]
    prof_rows[[length(prof_rows) + 1]] <- x[, .(period = per, group = g, centre_lab, n, vol_mean, pct_acute, variable = v,
      value = fifelse(rec == 0, "-", if (v %in% names(CONT)) sprintf("%.1f [%.1f-%.1f]", val, q1, q3) else sprintf("%.1f%%", val)),
      recorded = rec, missing_pct = round(100 * (1 - rec / n), 1), val)]
    fig_rows[[length(fig_rows) + 1]] <- x[rec >= MIN_REC, .(period = per, group = g, centre_lab, vol_mean, variable = v, val, rec)]
  }
  # ---- B. between-centre comparison ----
  rhs_adj <- paste(c("year_c", covar_terms(COVARS_M1, specs), if (g == "All") "acute"), collapse = " + ")
  for (v in names(LAB)) {
    dd <- d[!is.na(get(v))]
    keep <- dd[, .N, by = centre_lab][N >= MIN_REC, centre_lab]; dd <- dd[centre_lab %in% keep]
    if (length(keep) < 3) next
    isb <- v %in% names(BIN)
    cv <- dd[, .(val = if (isb) 100 * mean(get(v) == 1) else as.numeric(median(get(v))), vol_mean = vol_mean[1]), by = centre_lab]
    nat <- if (isb) 100 * mean(dd[[v]] == 1) else as.numeric(median(dd[[v]]))
    p_het <- if (isb) suppressWarnings(chisq.test(table(dd$centre_lab, dd[[v]]))$p.value) else
             kruskal.test(dd[[v]], factor(dd$centre_lab))$p.value
    sp <- suppressWarnings(cor.test(cv$val, cv$vol_mean, method = "spearman"))
    icc_of <- function(rhs, dat) safe({
      f <- as.formula(paste(v, "~", rhs, "+ (1 | centrum)"), env = globalenv())
      m <- if (isb) lme4::glmer(f, dat, binomial, control = lme4::glmerControl(optimizer = "bobyqa")) else lme4::lmer(f, dat)
      s2 <- as.numeric(lme4::VarCorr(m)$centrum)
      c(icc = if (isb) s2 / (s2 + pi^2 / 3) else s2 / (s2 + sigma(m)^2), mor = if (isb) exp(sqrt(2 * s2) * qnorm(.75)) else NA_real_)
    }, paste("ICC", per, g, v))
    da <- dd[complete.cases(dd[, intersect(all.vars(as.formula(paste("~", rhs_adj))), names(dd)), with = FALSE])]
    i0 <- icc_of("1", dd); i1 <- if (nrow(da) > 50) icc_of(rhs_adj, da) else NULL
    vr <- safe({
      dd[, vol10 := vol_mean / 10]
      m <- if (isb) glm(as.formula(paste(v, "~ vol10")), binomial, dd) else lm(as.formula(paste(v, "~ vol10")), dd)
      ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$centrum, test = "Satterthwaite"); ct <- ct[ct$Coef == "vol10", ]
      tq <- qt(.975, ct$df_Satt); e <- c(ct$beta, ct$beta - tq * ct$SE, ct$beta + tq * ct$SE)
      c(if (isb) exp(e) else e, ct$p_Satt)
    }, paste("volume", per, g, v))
    if (is.null(vr)) vr <- rep(NA_real_, 4)
    q <- quantile(cv$val, c(0, .25, .5, .75, 1))
    cmp_rows[[length(cmp_rows) + 1]] <- data.table(period = per, group = g, variable = v, label = LAB[[v]],
      type = if (isb) "%" else "median", patients = nrow(dd), centres = length(keep), national = nat,
      centre_min = q[1], centre_q1 = q[2], centre_median = q[3], centre_q3 = q[4], centre_max = q[5],
      p_heterogeneity = p_het, icc_crude = i0[1], icc_adjusted = i1[1], mor_crude = i0[2], mor_adjusted = i1[2],
      spearman_volume = unname(sp$estimate), p_spearman = sp$p.value,
      volume_effect = if (isb) "OR per 10 operations/year" else "difference per 10 operations/year",
      volume_est = vr[1], volume_lower = vr[2], volume_upper = vr[3], p_volume_cr2 = vr[4])
    # ---- C. funnel ----
    if (isb) {
      p0 <- mean(dd[[v]] == 1)
      fun_rows[[length(fun_rows) + 1]] <- dd[, .(period = per, group = g, variable = v, n = .N, pct = 100 * mean(get(v) == 1)), by = centre_lab][
        , `:=`(national = 100 * p0, lo95 = 100 * pmax(0, p0 - 1.96 * sqrt(p0 * (1 - p0) / n)), hi95 = 100 * pmin(1, p0 + 1.96 * sqrt(p0 * (1 - p0) / n)),
               lo998 = 100 * pmax(0, p0 - 3.09 * sqrt(p0 * (1 - p0) / n)), hi998 = 100 * pmin(1, p0 + 3.09 * sqrt(p0 * (1 - p0) / n)))][
        , position := fcase(pct > hi998, "above 99.8% limit", pct < lo998, "below 99.8% limit", pct > hi95, "above 95% limit",
                            pct < lo95, "below 95% limit", default = "within limits")][]
    }
  }
}
PR <- rbindlist(prof_rows); CMP <- rbindlist(cmp_rows); FUN <- rbindlist(fun_rows); FIG <- rbindlist(fig_rows)
write_out(PR[, !"val"], SUB, "A_centre_profile_long")
wide <- dcast(PR, period + group + centre_lab + n + vol_mean + pct_acute ~ variable, value.var = "value")
setcolorder(wide, c("period", "group", "centre_lab", "n", "vol_mean", "pct_acute", intersect(names(LAB), names(wide))))
wide <- wide[order(period, group, -vol_mean)]
write_out(wide, SUB, "A_centre_profile_wide")
write_out(dcast(PR, period + group + centre_lab ~ variable, value.var = "missing_pct"), SUB, "A_centre_missing_pct")
write_out(CMP, SUB, "B_between_centre_comparison")
write_out(FUN, SUB, "C_funnel_binary")

# ---- Figure: centre values per variable, centres ordered by volume ---------------------------------
for (per in names(PERIODS)) for (g in GROUPS) {
  x <- FIG[period == per & group == g]; if (!nrow(x)) next
  ord <- unique(x[, .(centre_lab, vol_mean)])[order(vol_mean)]
  x[, centre := factor(centre_lab, levels = ord$centre_lab)]
  nat <- CMP[period == per & group == g, .(variable, national)]
  x[, variable := factor(variable, levels = names(LAB), labels = LAB)]; nat[, variable := factor(variable, levels = names(LAB), labels = LAB)]
  p <- ggplot(x, aes(val, centre)) +
    geom_vline(data = nat, aes(xintercept = national), colour = "#C0392B", linetype = "22") +
    geom_point(aes(size = rec), colour = "#1B5E4B") + scale_size_area(max_size = 3.5, name = "Recorded") +
    facet_wrap(~ variable, scales = "free_x", ncol = 4) +
    labs(x = NULL, y = "Centre (ordered by mean annual volume, highest at top)",
         title = sprintf("%s, %s", g, per), caption = "Red dashed line: national value. Centres with fewer than 10 recorded values are not shown.") +
    theme_pub() + theme(plot.background = element_rect(fill = "white", colour = NA),
                        panel.grid.major.y = element_line(colour = "#E8F3EC"), axis.text.y = element_text(size = 8))
  save_fig(p, SUB, sprintf("centre_practice_%s_%s", tolower(g), sub("-", "_", per)), 14, 12)
}

# ---- Overview ------------------------------------------------------------------------------------------
f2 <- function(z) ifelse(is.na(z), NA, ifelse(abs(z) >= 100, round(z, 1), signif(z, 3)))
out <- c(sprintf("Perfusion practice per centre - %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
for (per in names(PERIODS)) for (g in GROUPS) {
  w <- wide[period == per & group == g]
  out <- c(out, sprintf("==== A. Centrumprofiel: %s, %s (centra op volume, hoog naar laag) ====", g, per),
           capture.output(print(w[, !c("period", "group")], row.names = FALSE)), "")
}
cm <- copy(CMP)[, (names(CMP)[vapply(CMP, is.numeric, logical(1))]) := lapply(.SD, f2), .SDcols = is.numeric]
out <- c(out, "==== B. Vergelijking tussen centra ====",
         capture.output(print(cm[, .(period, group, variable, type, patients, centres, national, centre_min, centre_q1, centre_median,
                                     centre_q3, centre_max, p_heterogeneity, icc_crude, icc_adjusted, mor_crude, mor_adjusted,
                                     spearman_volume, p_spearman, volume_est, volume_lower, volume_upper, p_volume_cr2)], nrows = 500, row.names = FALSE)),
         "", "==== C. Funnel: centra buiten de 99.8%-grenzen (binaire variabelen) ====",
         capture.output(print(FUN[grepl("99.8", position), .(period, group, variable, centre_lab, n, pct = round(pct, 1),
                                                           national = round(national, 1), position)], nrows = 500, row.names = FALSE)))
f_out <- file.path(out_dir(SUB), "overzicht_centre_perfusion.txt")
writeLines(out, f_out); cat(out, sep = "\n")
log_msg("Done. Overview in ", f_out)

"OK" }, error = function(e) paste("FOUT:", conditionMessage(e)))
while (sink.number() > 0) sink()
cat("\nPart 1: perfusion practice per centre: ", .status_nieuw["Part 1: perfusion practice per centre"], " (", round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1), " min)\n", sep = "")

# #############################################################################
# Part 2: paper 1 sensitivity, strict groups  (= sens_strict_groups.R)
# #############################################################################
cat("\n\n#############  Part 2: paper 1 sensitivity, strict groups  #############\n"); .t0 <- Sys.time()
.status_nieuw["Part 2: paper 1 sensitivity, strict groups"] <- tryCatch({
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

"OK" }, error = function(e) paste("FOUT:", conditionMessage(e)))
while (sink.number() > 0) sink()
cat("\nPart 2: paper 1 sensitivity, strict groups: ", .status_nieuw["Part 2: paper 1 sensitivity, strict groups"], " (", round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1), " min)\n", sep = "")

# #############################################################################
# Part 3: paper 2 volume-outcome  (= volume.R)
# #############################################################################
cat("\n\n#############  Part 3: paper 2 volume-outcome  #############\n"); .t0 <- Sys.time()
.status_nieuw["Part 3: paper 2 volume-outcome"] <- tryCatch({
# =============================================================================
# volume.R  (standalone; not part of run_all.R)
# Run from the project folder after the main pipeline has created data/derived:
#   source("volume.R")
# Paper 2: volume-outcome relation for open aortic arch surgery, 2018-2024
# (period with complete classification and registration). Same analysis cohort as
# paper 1. Groups for paper 2 (default, ARCH_VOLUME_GROUPS = "strict"):
#   Acute    = urgentie 30 (emergency) or 40 (salvage), or acute dissection (code 20)
#              in any aortic segment (identical to the acute group of paper 1);
#   Elective = urgentie 10 (elective) without a dissection code (20) in any segment.
#              Chronic dissection (30) and infection (40) are already excluded from the
#              cohort. Urgent operations (urgentie 20) and operations with unknown
#              urgency belong to neither group and are left out of the group analyses;
#              their number is reported (section 0).
#   ARCH_VOLUME_GROUPS = "original" restores acute / non-acute as in paper 1.
#
# Methods follow the volume-outcome literature (Kawczynski et al., Eur Heart J 2023;
# Heuts et al., Eur Respir J 2025; Govers et al., EJCTS 2026):
#   1. Volume per centre: (a) centre-year volume in the same group (primary,
#      patient-level exposure), (b) centre mean annual volume 2018-2024 (centre
#      characteristic, used for tertiles and as sensitivity), (c) total centre
#      volume (all arch operations), (d) previous-year volume.
#   2. Descriptive: outcomes by volume tertile of centres; adjusted OR/HR per tertile.
#   3. Linear: OR per 10 operations/year, CR2 cluster-robust (Satterthwaite) and
#      random-intercept sensitivity; within/between-centre decomposition (Mundlak).
#   4. Non-linearity: restricted cubic splines with 3 (primary) and 4 knots;
#      CR2 Wald tests of any volume effect and of non-linearity; AIC vs linear.
#      Davies test (segmented::davies.test) for a change in slope (breakpoint) of the
#      volume term in the adjusted logistic model; it treats patients as independent,
#      so it is a sensitivity check next to the cluster-robust spline tests.
#   5. Threshold: (a) elbow method on the case-mix-standardised risk curve (point of
#      maximum distance to the chord, as in Kawczynski/Heuts), (b) hinge model
#      (risk falls linearly up to k, flat above), both with 95% intervals from a
#      centre-level cluster bootstrap (whole centres resampled). A threshold is only
#      labelled "supported" if the overall volume effect is significant, the curve
#      decreases, and the bootstrap places it away from the edge of the data.
#   6. Number needed to treat: standardised risk at 10 operations/year versus at the
#      threshold (as in Kawczynski et al.).
#   Outcomes: in-hospital (primary), 30-day and one-year mortality, primary composite,
#   major complication, failure to rescue (adjusted for number of complications).
#   7. Long-term survival: Cox models (robust SE by centre; gamma frailty
#      sensitivity), linear, spline and tertiles.
#   8. Practice by threshold: centres below versus at/above the volume threshold X
#      (centre mean annual volume in the same group) compared on perfusion practice
#      (nadir temperature, deep hypothermia, circulatory arrest and its duration, CPB
#      and cross-clamp time), procedure, case mix and outcomes. X = the supported elbow
#      threshold for in-hospital mortality; if none is supported, the median centre
#      volume is used and flagged. Override with Sys.setenv(ARCH_VOLUME_THRESHOLD = "x")
#      (one value) or "acute,elective" (two values). Patient-level tests (Wilcoxon/chi-square)
#      and centre-adjusted regressions with CR2 standard errors (15 centres).
# Adjustment throughout: age (spline), sex, logit EuroSCORE II, previous cardiac
# surgery, calendar year. 15 centres: all centre-level inference is limited.
# Output: results/volume/*.csv, figures, and overzicht_volume.txt (summary also shown in the console)
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival); library(ggplot2) })
start_log("volume")
log_msg("Missing-data method: ", method_label())
SUB <- "volume"
P0 <- 2018; P1 <- max(CONFIG$study_years)
B_BOOT <- as.integer(Sys.getenv("ARCH_VOLUME_BOOT", "500"))
set.seed(CONFIG$seed)

cohort <- load_derived("cohort.rds")
surg   <- load_derived("surgical_cohort.rds")
specs  <- load_derived("specs.rds")
d0     <- copy(analysis_datasets()[[1]])
GROUPS <- Sys.getenv("ARCH_VOLUME_GROUPS", "strict")
GRPS   <- if (GROUPS == "original") c("Acute", "Non-acute") else c("Acute", "Elective")
grp_of <- function(d) if (GROUPS == "original") as.character(d$urgency2) else
  fcase(d$acute %in% 1L, "Acute", d$urgentie %in% 10 & !(d$seg_ad %in% TRUE), "Elective", default = NA_character_)
cohort <- copy(cohort)[, grp := grp_of(cohort)]
d0[, grp := grp_of(d0)]
log_msg("Groups (", GROUPS, "): ", paste(names(table(cohort$grp, useNA = "ifany")), table(cohort$grp, useNA = "ifany"), sep = " = ", collapse = "; "))
rhs    <- paste(c("year_c", covar_terms(COVARS_M1, specs)), collapse = " + ")
OUTS   <- c(dzh = "In-hospital mortality", d30 = "30-day mortality", comp_primary = "Primary composite",
            d1y = "One-year mortality", major_compl = "Major complication",
            ftr = "Failure to rescue (in-hospital death after major complication)")
d0[, ftr := fifelse(major_compl %in% 1, dzh, NA_integer_)]
OUTS   <- OUTS[names(OUTS) %in% names(d0)]
rhs_for <- function(v) if (v == "ftr") paste(rhs, "+ n_compl") else rhs
end_of <- function(v) if (v == "d1y") min(CONFIG$d1y_last_year, P1) else P1
safe   <- function(expr, what) tryCatch(expr, error = function(e) { log_msg("  skipped (", what, "): ", conditionMessage(e)); NULL })

# ---- 1. Volume measures --------------------------------------------------------------------
d0[, centrum_chr := as.character(centrum)]
vg <- cohort[!is.na(grp), .N, by = .(centrum_chr = as.character(centrum), year, grp)]
vt <- surg[, .(vol_total = .N), by = .(centrum_chr = as.character(centrum), year)]
vp <- copy(vg)[, year := year + 1L][, .(centrum_chr, year, grp, vol_prev = N)]
d0 <- merge(d0, vg[, .(centrum_chr, year, grp, vol_group = N)], by = c("centrum_chr", "year", "grp"), all.x = TRUE)
d0 <- merge(d0, vt, by = c("centrum_chr", "year"), all.x = TRUE)
d0 <- merge(d0, vp, by = c("centrum_chr", "year", "grp"), all.x = TRUE)
d0[is.na(vol_prev), vol_prev := 0]
cm <- vg[year >= P0 & year <= P1, .(vol_mean = sum(N) / (P1 - P0 + 1)), by = .(centrum_chr, grp)]
d0 <- merge(d0, cm, by = c("centrum_chr", "grp"), all.x = TRUE)
d0 <- d0[year >= P0 & year <= P1]
write_out(cm[order(grp, -vol_mean)], SUB, "centre_mean_annual_volume")
G0 <- cohort[year >= P0 & year <= P1, .(n = .N), by = .(grp = fifelse(is.na(grp), "Neither group", grp),
             paper1_group = as.character(urgency2), urgentie, dissection_code = seg_ad %in% TRUE)][order(grp, paper1_group, urgentie)]
write_out(G0, SUB, "group_definition_counts")

VOLS <- c(vol_group = "Centre-year volume, same group", vol_mean = "Centre mean annual volume 2018-2024",
          vol_total = "Centre-year volume, all arch surgery", vol_prev = "Previous-year volume, same group")

model_data <- function(gr, v) {
  d <- d0[grp == gr & year <= end_of(v) & !is.na(get(v))]
  d[complete.cases(d[, c(COVARS_M1, "year_c", "centrum", if (v == "ftr") "n_compl"), with = FALSE])]
}
cr2_term <- function(m, d, term) {
  r <- clubSandwich::coef_test(m, vcov = "CR2", cluster = d$centrum, test = "Satterthwaite")
  r <- r[r$Coef == term, ]
  data.table(or = exp(r$beta), lower = exp(r$beta - qt(.975, r$df_Satt) * r$SE),
             upper = exp(r$beta + qt(.975, r$df_Satt) * r$SE), p = r$p_Satt)
}
rcs_basis <- function(x, knots) {
  k <- length(knots); nrm <- (knots[k] - knots[1])^2
  b <- sapply(seq_len(k - 2), function(j) (pmax(x - knots[j], 0)^3 -
         pmax(x - knots[k - 1], 0)^3 * (knots[k] - knots[j]) / (knots[k] - knots[k - 1]) +
         pmax(x - knots[k], 0)^3 * (knots[k - 1] - knots[j]) / (knots[k] - knots[k - 1])) / nrm)
  m <- cbind(x, matrix(b, ncol = k - 2)); colnames(m) <- paste0("vs", seq_len(ncol(m))); m
}
knots_for <- function(x, nk) unique(as.numeric(quantile(x, if (nk == 3) c(.10, .50, .90) else c(.05, .35, .65, .95))))
add_basis <- function(d, x, kn) { B <- rcs_basis(x, kn); d[, (colnames(B)) := as.data.table(B)]; colnames(B) }
elbow <- function(x, y) {                       # maximum distance to the chord (Kneedle)
  xs <- (x - min(x)) / diff(range(x)); ys <- (y - min(y)) / max(diff(range(y)), 1e-12)
  chord <- ys[1] + (ys[length(ys)] - ys[1]) * xs
  x[which.max(abs(chord - ys))]
}

# ---- 2. Descriptive: tertiles of centres -----------------------------------------------------
tert <- cm[, .(centrum_chr, vol_mean, tertile = cut(rank(vol_mean, ties.method = "first"), 3, labels = c("T1 (low)", "T2", "T3 (high)"))),
           by = grp]
d0 <- merge(d0, tert[, .(centrum_chr, grp, tertile)], by = c("centrum_chr", "grp"), all.x = TRUE)
write_out(tert[, .(centres = .N, vol_min = min(vol_mean), vol_max = max(vol_mean)), by = .(grp, tertile)][order(grp, tertile)],
          SUB, "tertile_definitions")
tert_rows <- list()
for (gr in GRPS) for (v in names(OUTS)) {
  d <- model_data(gr, v); if (sum(d[[v]] == 1) < 20) next
  crude <- d[, .(n = .N, events = sum(get(v) == 1), pct = 100 * mean(get(v) == 1)), by = tertile][order(tertile)]
  d[, tert_num := as.integer(tertile) - 1L]
  m <- glm(as.formula(paste(v, "~ tertile +", rhs_for(v))), binomial, d)
  t2 <- safe(cr2_term(m, d, "tertileT2"), "T2"); t3 <- safe(cr2_term(m, d, "tertileT3 (high)"), "T3")
  mt <- glm(as.formula(paste(v, "~ tert_num +", rhs_for(v))), binomial, d); tr <- safe(cr2_term(mt, d, "tert_num"), "trend")
  tert_rows[[length(tert_rows) + 1]] <- cbind(group = gr, outcome = v, crude,
    or_vs_T1 = c(1, t2$or, t3$or), lower = c(NA, t2$lower, t3$lower), upper = c(NA, t2$upper, t3$upper), p = c(NA, t2$p, t3$p),
    or_per_tertile = tr$or, or_per_tertile_lower = tr$lower, or_per_tertile_upper = tr$upper, p_trend = tr$p)
}
write_out(rbindlist(tert_rows, fill = TRUE), SUB, "outcomes_by_volume_tertile")

# ---- 3-4. Linear, spline, within/between ----------------------------------------------------------
lin <- list(); spl <- list(); curves <- list(); wb <- list()
for (gr in GRPS) for (v in names(OUTS)) for (vv in names(VOLS)) {
  d <- model_data(gr, v); ev <- sum(d[[v]] == 1)
  if (ev < 20) { log_msg("  ", gr, "/", v, ": too few events (", ev, ")"); next }
  log_msg(gr, " / ", v, " / ", vv)
  d[, vol := get(vv)]; d[, vol10 := vol / 10]
  lab <- data.table(group = gr, outcome = v, volume = VOLS[[vv]], n = nrow(d), events = ev,
                    vol_p10 = quantile(d$vol, .1), vol_median = median(d$vol), vol_p90 = quantile(d$vol, .9))
  m_lin <- glm(as.formula(paste(v, "~ vol10 +", rhs_for(v))), binomial, d)
  a <- cr2_term(m_lin, d, "vol10")
  ri <- safe({ f <- fit_mixed(as.formula(paste(v, "~ vol10 +", rhs_for(v), "+ (1 | centrum)")), d, binomial())
               e <- summary(f$fit)$coefficients["vol10", ]; c(exp(e[1]), exp(e[1] - 1.96 * e[2]), exp(e[1] + 1.96 * e[2]), e[4]) }, "random intercept")
  if (is.null(ri)) ri <- rep(NA_real_, 4)
  lin[[length(lin) + 1]] <- cbind(lab, or_per_10 = a$or, lower = a$lower, upper = a$upper, p_cr2 = a$p,
                                  or_per_10_ri = ri[1], ri_lower = ri[2], ri_upper = ri[3], ri_p = ri[4], aic_linear = AIC(m_lin))
  for (nk in c(3, 4)) {
    kn <- knots_for(d$vol, nk); if (length(kn) < nk) next
    bn <- add_basis(d, d$vol, kn)
    m_s <- glm(as.formula(paste(v, "~", paste(bn, collapse = " + "), "+", rhs_for(v))), binomial, d)
    V <- clubSandwich::vcovCR(m_s, cluster = d$centrum, type = "CR2")
    w_all <- clubSandwich::Wald_test(m_s, constraints = clubSandwich::constrain_zero(bn), vcov = V, test = "HTZ")
    w_nl  <- clubSandwich::Wald_test(m_s, constraints = clubSandwich::constrain_zero(bn[-1]), vcov = V, test = "HTZ")
    lr_nl <- deviance(m_lin) - deviance(m_s)
    spl[[length(spl) + 1]] <- cbind(lab, knots = nk, knot_values = paste(round(kn, 1), collapse = "/"),
      p_overall_cr2 = w_all$p_val, p_nonlinear_cr2 = w_nl$p_val,
      p_nonlinear_lr = pchisq(lr_nl, nk - 2, lower.tail = FALSE), aic_spline = AIC(m_s), aic_linear = AIC(m_lin))
    if (nk == 3 && vv %in% c("vol_group", "vol_mean")) {
      grid <- seq(quantile(d$vol, .05), quantile(d$vol, .95), length.out = 40)
      std <- vapply(grid, function(g) { nd <- copy(d); add_basis(nd, rep(g, nrow(nd)), kn); mean(predict(m_s, nd, type = "response")) }, 0)
      curves[[length(curves) + 1]] <- data.table(group = gr, outcome = v, volume = VOLS[[vv]], vol = grid, std_risk = std)
    }
  }
  if (vv == "vol_group") {                                  # Mundlak: between (centre mean) vs within (yearly deviation)
    d[, `:=`(vb10 = vol_mean / 10, vw10 = (vol_group - vol_mean) / 10)]
    m_wb <- glm(as.formula(paste(v, "~ vb10 + vw10 +", rhs_for(v))), binomial, d)
    b <- cr2_term(m_wb, d, "vb10"); w <- cr2_term(m_wb, d, "vw10")
    wb[[length(wb) + 1]] <- cbind(lab[, .(group, outcome, n, events)],
      or_between = b$or, between_lower = b$lower, between_upper = b$upper, p_between = b$p,
      or_within = w$or, within_lower = w$lower, within_upper = w$upper, p_within = w$p)
  }
}
write_out(rbindlist(lin), SUB, "linear_models")
write_out(rbindlist(spl), SUB, "spline_tests")
write_out(rbindlist(wb), SUB, "within_between_centre")
CV <- rbindlist(curves); write_out(CV, SUB, "standardised_risk_curves")

# ---- 4b. Davies test for a breakpoint in the volume slope ---------------------------------------------
dav <- list()
if (!requireNamespace("segmented", quietly = TRUE)) {
  log_msg("  Davies test skipped: package 'segmented' not installed (install.packages(\"segmented\"))")
} else for (gr in GRPS) for (v in names(OUTS)) for (vv in c("vol_group", "vol_mean")) {
  d <- model_data(gr, v); ev <- sum(d[[v]] == 1); if (ev < 20) next
  d[, vol := get(vv)]
  nv <- uniqueN(d$vol); if (nv < 5) next
  m <- glm(as.formula(paste(v, "~ vol +", rhs_for(v))), binomial, d)
  r <- safe(segmented::davies.test(m, seg.Z = ~ vol, k = min(10L, nv - 2L)), paste("Davies", gr, v, vv))
  if (is.null(r)) next
  dav[[length(dav) + 1]] <- data.table(group = gr, outcome = v, volume = VOLS[[vv]], n = nrow(d), events = ev,
    distinct_volumes = nv, best_breakpoint = unname(r$statistic), p_davies = r$p.value)
}
DV <- rbindlist(dav); if (nrow(DV)) write_out(DV, SUB, "davies_test_breakpoint")

# ---- 5-6. Thresholds: elbow and hinge, with centre bootstrap; NNT ------------------------------------
hinge_fit <- function(y, X, vol, ks) {
  dev <- vapply(ks, function(k) suppressWarnings(glm.fit(cbind(X, pmin(vol, k)), y, family = binomial()))$deviance, 0)
  i <- which.min(dev)
  fh <- suppressWarnings(glm.fit(cbind(X, pmin(vol, ks[i])), y, family = binomial()))
  dl <- suppressWarnings(glm.fit(cbind(X, vol), y, family = binomial()))$deviance
  c(k = ks[i], or10 = exp(10 * tail(fh$coefficients, 1)), lr_vs_linear = dl - dev[i], edge = as.numeric(i %in% c(1, length(ks))))
}
elbow_fit <- function(y, X, vol, kn, grid) {
  Bm <- rcs_basis(vol, kn)
  f <- suppressWarnings(glm.fit(cbind(X, Bm), y, family = binomial()))
  cf <- f$coefficients; cf[is.na(cf)] <- 0
  std <- vapply(grid, function(g) mean(plogis(drop(cbind(X, rcs_basis(rep(g, length(vol)), kn)) %*% cf))), 0)
  c(elbow = elbow(grid, std), r10 = std[which.min(abs(grid - 10))], decreasing = as.numeric(std[length(std)] < std[1]), std = std)
}
thr <- list(); reps_all <- list()
for (gr in GRPS) for (v in intersect(c("dzh", "d30", "d1y", "comp_primary"), names(OUTS))) for (vv in c("vol_group", "vol_mean")) {
  d <- model_data(gr, v); if (sum(d[[v]] == 1) < 30) next
  log_msg("  thresholds ", gr, " / ", v, " / ", vv, " (", B_BOOT, " bootstrap replicates)")
  X <- model.matrix(as.formula(paste("~", rhs)), d); y <- d[[v]]; vol <- d[[vv]]
  ks <- unique(round(quantile(vol, seq(.10, .90, by = .02))))
  kn <- knots_for(vol, 3); grid <- seq(quantile(vol, .05), quantile(vol, .95), length.out = 40)
  h0 <- hinge_fit(y, X, vol, ks); e0 <- elbow_fit(y, X, vol, kn, grid)
  std0 <- e0[grep("^std", names(e0))]
  r_at <- function(x) std0[which.min(abs(grid - x))]
  nnt <- if (e0[["decreasing"]] == 1 && r_at(10) > r_at(e0[["elbow"]])) 1 / (r_at(10) - r_at(e0[["elbow"]])) else NA_real_
  idx <- split(seq_len(nrow(d)), as.character(d$centrum))
  reps <- t(vapply(seq_len(B_BOOT), function(b) {
    rows <- unlist(idx[sample(names(idx), length(idx), replace = TRUE)], use.names = FALSE)
    h <- tryCatch(hinge_fit(y[rows], X[rows, , drop = FALSE], vol[rows], ks), error = function(e) rep(NA_real_, 4))
    e <- tryCatch(elbow_fit(y[rows], X[rows, , drop = FALSE], vol[rows], kn, grid)[1:3], error = function(e) rep(NA_real_, 3))
    c(h, e)
  }, numeric(7)))
  colnames(reps) <- c("k", "or10", "lr_vs_linear", "edge", "elbow", "r10", "decreasing")
  r <- reps[stats::complete.cases(reps), , drop = FALSE]
  reps_all[[length(reps_all) + 1]] <- data.table(group = gr, outcome = v, volume = VOLS[[vv]], r)
  q <- function(x, p) if (length(x)) as.numeric(quantile(x, p)) else NA_real_
  thr[[length(thr) + 1]] <- data.table(group = gr, outcome = v, volume = VOLS[[vv]], n = nrow(d), events = sum(y == 1),
    centres = length(idx), replicates_ok = nrow(r),
    elbow = e0[["elbow"]], elbow_low = q(r[, "elbow"], .025), elbow_high = q(r[, "elbow"], .975),
    curve_decreasing = e0[["decreasing"]] == 1, share_reps_decreasing = mean(r[, "decreasing"] == 1),
    std_risk_at_10 = r_at(10), std_risk_at_elbow = r_at(e0[["elbow"]]), nnt_10_vs_elbow = nnt,
    hinge_threshold = h0[["k"]], hinge_low = q(r[r[, "edge"] == 0, "k"], .025), hinge_high = q(r[r[, "edge"] == 0, "k"], .975),
    hinge_or_per_10_below = h0[["or10"]], hinge_or_low = q(r[, "or10"], .025), hinge_or_high = q(r[, "or10"], .975),
    share_reps_hinge_at_edge = mean(r[, "edge"] == 1),
    share_reps_hinge_better_than_linear = mean(r[, "lr_vs_linear"] > 3.84 & r[, "edge"] == 0 & r[, "or10"] < 1))
}
TH <- rbindlist(thr)
sp <- rbindlist(spl)[knots == 3, .(group, outcome, volume, p_overall_cr2, p_nonlinear_cr2)]
TH <- merge(TH, sp, by = c("group", "outcome", "volume"), all.x = TRUE)
TH[, threshold_supported := p_overall_cr2 < .05 & curve_decreasing & share_reps_decreasing >= .95 &
                            share_reps_hinge_at_edge < .5]
write_out(TH, SUB, "thresholds_elbow_hinge_bootstrap")
write_out(rbindlist(reps_all), SUB, "thresholds_bootstrap_replicates")

# ---- 7. Long-term survival (Cox) ---------------------------------------------------------------------
sv <- d0[!is.na(time_d) & !is.na(status)][, t_years := pmax(time_d, .5) / 365.25]
cox_rows <- list(); cox_curves <- list()
for (gr in GRPS) for (vv in c("vol_group", "vol_mean")) {
  d <- sv[grp == gr]; d <- d[complete.cases(d[, c(COVARS_M1, "year_c", "centrum"), with = FALSE])]
  d[, vol := get(vv)]; d[, vol10 := vol / 10]
  fl <- coxph(as.formula(paste("Surv(t_years, status) ~ vol10 +", rhs, "+ cluster(centrum)")), d)
  s <- summary(fl)$coefficients["vol10", ]
  fr <- safe(coxph(as.formula(paste("Surv(t_years, status) ~ vol10 +", rhs, "+ frailty(centrum)")), d), "frailty")
  sf <- if (is.null(fr)) rep(NA_real_, 6) else summary(fr)$coefficients["vol10", ]
  kn <- knots_for(d$vol, 3); bn <- add_basis(d, d$vol, kn)
  fs <- coxph(as.formula(paste("Surv(t_years, status) ~", paste(bn, collapse = " + "), "+", rhs, "+ cluster(centrum)")), d)
  Vr <- fs$var; bi <- match(bn, names(coef(fs)))
  wald <- function(ii) { b <- coef(fs)[ii]; drop(t(b) %*% solve(Vr[ii, ii, drop = FALSE]) %*% b) }
  rows <- list(group = gr, volume = VOLS[[vv]], n = nrow(d), deaths = sum(d$status == 1),
               hr_per_10 = exp(s[["coef"]]), lower = exp(s[["coef"]] - 1.96 * s[["robust se"]]), upper = exp(s[["coef"]] + 1.96 * s[["robust se"]]),
               p_robust = s[["Pr(>|z|)"]],
               hr_per_10_frailty = exp(sf[1]), p_frailty = if (is.null(fr)) NA_real_ else sf[length(sf)],
               p_overall_spline_robust = pchisq(wald(bi), length(bi), lower.tail = FALSE),
               p_nonlinear_spline_robust = pchisq(wald(bi[-1]), length(bi) - 1, lower.tail = FALSE))
  if (vv == "vol_mean") {
    d[, tert_num := as.integer(tertile) - 1L]
    ft <- coxph(as.formula(paste("Surv(t_years, status) ~ tert_num +", rhs, "+ cluster(centrum)")), d)
    st <- summary(ft)$coefficients["tert_num", ]
    rows <- c(rows, hr_per_tertile = exp(st[["coef"]]), hr_tertile_lower = exp(st[["coef"]] - 1.96 * st[["robust se"]]),
              hr_tertile_upper = exp(st[["coef"]] + 1.96 * st[["robust se"]]), p_tertile_trend = st[["Pr(>|z|)"]])
    km <- survfit(Surv(t_years, status) ~ tertile, d)
    cox_curves[[length(cox_curves) + 1]] <- data.table(group = gr, tertile = sub("tertile=", "", rep(names(km$strata), km$strata)),
                                                       time = km$time, surv = km$surv, n_risk = km$n.risk)
  }
  cox_rows[[length(cox_rows) + 1]] <- as.data.table(rows)
}
write_out(rbindlist(cox_rows, fill = TRUE), SUB, "long_term_survival_cox")
KMT <- rbindlist(cox_curves); write_out(KMT[n_risk >= KM_MIN_RISK], SUB, "km_by_volume_tertile")


# ---- 8. Practice by threshold: centres < X versus >= X ------------------------------------------------
thr_env <- as.numeric(strsplit(Sys.getenv("ARCH_VOLUME_THRESHOLD", ""), ",")[[1]])
pick_x <- function(gr) {
  if (length(thr_env) == 2) return(list(x = thr_env[if (gr == "Acute") 1 else 2], source = "set by user"))
  if (length(thr_env) == 1 && !is.na(thr_env)) return(list(x = thr_env, source = "set by user"))
  t <- TH[group == gr & outcome == "dzh" & volume == VOLS[["vol_mean"]]]
  if (nrow(t) && isTRUE(t$threshold_supported)) return(list(x = t$elbow, source = "supported elbow, in-hospital mortality"))
  list(x = median(cm[grp == gr, vol_mean]), source = "median centre volume (no supported threshold)")
}
PRAC_CONT <- c(temp_min = "Nadir temperature (C)", ca_min = "Circulatory arrest duration (min, arrest patients)",
               ecc_min = "Cardiopulmonary bypass time (min)", aox_min = "Cross-clamp time (min)",
               age = "Age (years)", euroscore2 = "EuroSCORE II (%)")
PRAC_BIN  <- c(deep = "Deep hypothermia (<=20 C)", circ_arrest_b = "Circulatory arrest",
               root_inv = "Aortic root involved", desc_inv = "Descending aorta involved",
               hybrid_fet = "Combined open + endovascular (incl. FET)", conc_valve = "Concomitant valve surgery",
               acute_dissection = "Acute aortic dissection", female = "Female sex", prior_cs = "Previous cardiac surgery",
               dzh = "In-hospital mortality", d30 = "30-day mortality", d1y = "One-year mortality", comp_primary = "Primary composite",
               major_compl = "Major complication", ftr = "Failure to rescue")
d0[, deep := fifelse(is.na(temp_min), NA_integer_, as.integer(temp_min <= 20))]
if (!"acute_dissection" %in% names(d0) && "path3" %in% names(d0)) d0[, acute_dissection := as.integer(path3 == "Acute dissection")]
cmp_rows <- list(); xinfo <- list()
for (gr in GRPS) {
  px <- pick_x(gr); X <- px$x
  d <- d0[grp == gr]; d[, side := factor(fifelse(vol_mean >= X, "high", "low"), c("low", "high"))]
  cen <- unique(d[, .(centrum_chr, vol_mean, side)])
  xinfo[[gr]] <- data.table(group = gr, threshold = X, threshold_source = px$source,
                            centres_below = sum(cen$side == "low"), centres_above = sum(cen$side == "high"),
                            patients_below = sum(d$side == "low"), patients_above = sum(d$side == "high"))
  if (min(table(cen$side)) < 2) { log_msg("  ", gr, ": fewer than 2 centres on one side of X = ", round(X, 1), "; comparison skipped"); next }
  for (vname in intersect(names(PRAC_CONT), names(d))) {
    dd <- d[!is.na(get(vname))]; if (vname == "ca_min") dd <- dd[circ_arrest_b == 1 & ca_min > 0]
    if (uniqueN(dd$side) < 2) next
    s <- dd[, .(n = .N, median = as.numeric(median(get(vname))), q1 = as.numeric(quantile(get(vname), .25)), q3 = as.numeric(quantile(get(vname), .75))), by = side]
    pw <- suppressWarnings(wilcox.test(as.formula(paste(vname, "~ side")), dd)$p.value)
    m <- lm(as.formula(paste(vname, "~ side")), dd)
    r <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$centrum, test = "Satterthwaite"); r <- r[r$Coef == "sidehigh", ]
    cmp_rows[[length(cmp_rows) + 1]] <- data.table(group = gr, threshold = X, variable = PRAC_CONT[[vname]], type = "median [IQR]",
      below = sprintf("%.1f [%.1f-%.1f]", s[side == "low", median], s[side == "low", q1], s[side == "low", q3]),
      above = sprintf("%.1f [%.1f-%.1f]", s[side == "high", median], s[side == "high", q1], s[side == "high", q3]),
      n_below = s[side == "low", n], n_above = s[side == "high", n], p_patient = pw,
      difference_mean = r$beta, p_centre_adjusted_cr2 = r$p_Satt)
  }
  for (vname in intersect(names(PRAC_BIN), names(d))) {
    dd <- d[!is.na(get(vname))]; if (uniqueN(dd$side) < 2 || sum(dd[[vname]] == 1) < 5) next
    s <- dd[, .(n = .N, k = sum(get(vname) == 1)), by = side]
    pc <- suppressWarnings(chisq.test(table(dd$side, dd[[vname]]))$p.value)
    rhs_v <- if (vname %in% c("dzh", "d30", "d1y", "comp_primary", "major_compl", "ftr")) paste("+", rhs_for(vname)) else ""
    m <- glm(as.formula(paste(vname, "~ side", rhs_v)), binomial, dd)
    r <- tryCatch(cr2_term(m, dd, "sidehigh"), error = function(e) data.table(or = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
    cmp_rows[[length(cmp_rows) + 1]] <- data.table(group = gr, threshold = X, variable = PRAC_BIN[[vname]],
      type = if (rhs_v == "") "n (%); OR high vs low, unadjusted" else "n (%); OR high vs low, case-mix adjusted",
      below = sprintf("%d (%.1f)", s[side == "low", k], 100 * s[side == "low", k] / s[side == "low", n]),
      above = sprintf("%d (%.1f)", s[side == "high", k], 100 * s[side == "high", k] / s[side == "high", n]),
      n_below = s[side == "low", n], n_above = s[side == "high", n], p_patient = pc,
      or_high_vs_low = r$or, or_lower = r$lower, or_upper = r$upper, p_centre_adjusted_cr2 = r$p)
  }
}
write_out(rbindlist(xinfo), SUB, "threshold_used_for_comparison")
CMP <- rbindlist(cmp_rows, fill = TRUE); write_out(CMP, SUB, "practice_below_vs_above_threshold")

# ---- Figures --------------------------------------------------------------------------------------
GC <- setNames(c("#1B5E4B", "#5FAE88"), GRPS)
if (nrow(CV)) {
  cvp <- CV[volume == VOLS[["vol_group"]] & outcome %in% c("dzh", "d30", "d1y")]
  ths <- TH[volume == VOLS[["vol_group"]] & outcome %in% c("dzh", "d30", "d1y")]
  p <- ggplot(cvp, aes(vol, 100 * std_risk, colour = group)) + geom_line(linewidth = .9) +
    geom_vline(data = ths[threshold_supported == TRUE], aes(xintercept = elbow, colour = group), linetype = "22") +
    facet_wrap(~ outcome, scales = "free_y", labeller = as_labeller(OUTS)) + scale_colour_manual(values = GC, name = NULL) +
    labs(x = "Centre volume in the same group (operations per year)", y = "Case-mix-standardised risk (%)",
         caption = "Dashed line: elbow threshold, shown only where supported (see thresholds_elbow_hinge_bootstrap.csv)") + theme_pub()
  save_fig(p, SUB, "volume_outcome_spline_curves", 10, 4.5)
}
if (nrow(KMT)) {
  p <- ggplot(KMT[n_risk >= KM_MIN_RISK], aes(time, surv, colour = tertile)) + geom_step(linewidth = .7) + facet_wrap(~ group) +
    scale_colour_manual(values = c("#B5DAC6", "#5FAE88", "#1B5E4B"), name = "Centre volume") +
    scale_y_continuous(labels = scales::percent) + labs(x = "Years after surgery", y = "Survival") + theme_pub()
  save_fig(p, SUB, "survival_by_volume_tertile", 9, 4.5)
}
# ---- Summary in the console and in results/volume/overzicht_volume.txt --------------------------------
OUTD <- out_dir(SUB)
show <- function(title, f, cols = NULL) {
  x <- fread(file.path(OUTD, f)); if (!is.null(cols)) x <- x[, intersect(cols, names(x)), with = FALSE]
  num <- names(x)[vapply(x, is.numeric, logical(1))]
  x[, (num) := lapply(.SD, function(z) ifelse(abs(z) >= 100, round(z, 1), signif(z, 3))), .SDcols = num]
  cat("\n==== ", title, " ====\n", sep = ""); print(x, row.names = FALSE)
}
options(width = 250)
while (sink.number() > 0) sink()
sink(file.path(OUTD, "overzicht_volume.txt"), split = TRUE)
cat("Volume-outcome, ", P0, "-", P1, " - groups: ", paste(GRPS, collapse = " / "), " (", GROUPS, ") - ", format(Sys.time()), "\n", sep = "")
show("0. Groepsindeling (operaties 2018-2024)", "group_definition_counts.csv")
show("1. Volumetertielen (centra)", "tertile_definitions.csv")
show("2. Uitkomst per volumetertiel", "outcomes_by_volume_tertile.csv",
     c("group", "outcome", "tertile", "n", "events", "pct", "or_vs_T1", "lower", "upper", "p", "or_per_tertile", "p_trend"))
show("3. Lineair (OR per 10 operaties/jaar)", "linear_models.csv",
     c("group", "outcome", "volume", "n", "events", "vol_median", "or_per_10", "lower", "upper", "p_cr2", "or_per_10_ri", "ri_p"))
show("4. Spline / niet-lineariteit", "spline_tests.csv",
     c("group", "outcome", "volume", "knots", "knot_values", "p_overall_cr2", "p_nonlinear_cr2", "p_nonlinear_lr", "aic_spline", "aic_linear"))
if (file.exists(file.path(OUTD, "davies_test_breakpoint.csv")))
  show("4b. Davies-test (knikpunt in de volumehelling; patienten als onafhankelijk behandeld)", "davies_test_breakpoint.csv")
show("5. Drempel (elbow + hinge, cluster-bootstrap) en NNT", "thresholds_elbow_hinge_bootstrap.csv",
     c("group", "outcome", "volume", "elbow", "elbow_low", "elbow_high", "curve_decreasing", "share_reps_decreasing",
       "std_risk_at_10", "std_risk_at_elbow", "nnt_10_vs_elbow", "hinge_threshold", "hinge_low", "hinge_high",
       "hinge_or_per_10_below", "share_reps_hinge_better_than_linear", "p_overall_cr2", "p_nonlinear_cr2", "threshold_supported"))
show("6. Binnen- vs tussen-centrum", "within_between_centre.csv")
show("7. Overleving lange termijn (Cox)", "long_term_survival_cox.csv")
show("8a. Gebruikte drempel X per groep", "threshold_used_for_comparison.csv")
show("8b. Centra < X versus >= X: perfusie, operatie, case-mix, uitkomst", "practice_below_vs_above_threshold.csv",
     c("group", "variable", "below", "above", "p_patient", "difference_mean", "or_high_vs_low", "or_lower", "or_upper", "p_centre_adjusted_cr2"))
sink()
log_msg("Done. Results in ", OUTD)

"OK" }, error = function(e) paste("FOUT:", conditionMessage(e)))
while (sink.number() > 0) sink()
cat("\nPart 3: paper 2 volume-outcome: ", .status_nieuw["Part 3: paper 2 volume-outcome"], " (", round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1), " min)\n", sep = "")
# 6. Where the overviews are
cat("\n\n============ KLAAR ============\n")
print(.status_nieuw)
cat("\nOverzichten (stuur hiervan screenshots):\n")
for (f in list.files("results", pattern = "^overzicht_(centre_perfusion|strict_groups|volume)\\.txt$",
                     recursive = TRUE, full.names = TRUE))
  cat("  ", normalizePath(f), "\n")

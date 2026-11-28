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

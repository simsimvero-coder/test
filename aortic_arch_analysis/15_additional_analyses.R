# =============================================================================
# 15_additional_analyses.R
# Additional analyses requested for the three papers.
# Paper 1: A. one-year mortality trend, sensitivity analyses (operations up to
#             CONFIG$d1y_last_year); B. missing data by year and by centre (anonymised).
# Paper 2 (volume-outcome; same-group annual centre volume, as in 13):
#          C. prior-year volume; D. within-centre vs between-centre volume;
#          E. elective (non-acute) volume as predictor of acute outcomes;
#          F. complications and failure to rescue by volume;
#          G. one-year mortality and long-term survival by volume.
# Paper 3 (practice and outcome, patient level, without centre effect in the
#          estimate; cluster-robust CR1 SE by centre with t(G-1), as in 11 section E):
#          H. adjustment for operative complexity, each variable separately and all
#             together (the "All" group is also adjusted for acute vs non-acute);
#          I. sensitivity excluding centres with > 50% missing nadir temperature;
#          J. shape of the arrest-duration association (categories and spline);
#          K. arrest duration x nadir temperature interaction;
#          L. renal failure as an intermediate between arrest duration and death
#             (descriptive attenuation, not formal causal mediation).
# All analyses are exploratory. Output: results/<method>/15_additional/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("15_additional_analyses")
log_msg("Missing-data method: ", method_label())
SUB <- "15_additional"
suppressMessages(library(survival))

cohort <- load_derived("cohort.rds")
specs  <- load_derived("specs.rds")
d0     <- copy(analysis_datasets()[[1]])
d0[, centrum := as.character(centrum)]
G      <- uniqueN(d0$centrum)
m1     <- covar_terms(COVARS_M1, specs)
fmt    <- function(o, l, u) sprintf("%.2f (%.2f-%.2f)", o, l, u)

cl_term <- function(m, cluster, term, type = "CR1") {
  if (type == "CR2") {
    r <- clubSandwich::coef_test(m, vcov = "CR2", cluster = cluster, test = "Satterthwaite")
    r <- r[r$Coef == term, ]
    return(data.table(or = exp(r$beta), lower = exp(r$beta - qt(.975, r$df_Satt) * r$SE),
                      upper = exp(r$beta + qt(.975, r$df_Satt) * r$SE), p = r$p_Satt))
  }
  V <- sandwich::vcovCL(m, cluster = cluster, type = "HC1")
  b <- coef(m)[[term]]; se <- sqrt(V[term, term]); df <- uniqueN(cluster) - 1
  data.table(or = exp(b), lower = exp(b - qt(.975, df) * se), upper = exp(b + qt(.975, df) * se),
             p = 2 * pt(-abs(b / se), df))
}
glm_or <- function(f, x, term, type = "CR1") {
  m <- glm(f, binomial, x)
  if (!term %in% names(coef(m)) || is.na(coef(m)[[term]])) return(NULL)
  cbind(cl_term(m, x$centrum, term, type), n = nrow(x), events = sum(model.response(model.frame(m)) == 1))
}
cc_rows <- function(x, f) x[complete.cases(x[, .SD, .SDcols = intersect(all.vars(f), names(x))])]

# ---- A. One-year mortality: trend and sensitivity analyses --------------------------------------
log_msg("A. One-year mortality sensitivity analyses")
d1 <- d0[!is.na(d1y) & year <= CONFIG$d1y_last_year]
write_out(d1[, .(n = .N, deaths = sum(d1y == 1), mortality_pct = round(100 * mean(d1y), 1)), by = urgency2],
          SUB, "one_year_mortality_crude")
run1 <- function(lab, covs, dd = d1, rnd = "(1 | centrum)", t = "year_c") {
  r <- tryCatch(fit_pooled(build_formula("d1y", covs, specs, random = rnd), list(dd)),
                error = function(e) { log_msg("  ", lab, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(data.table(analysis = lab, term = t))
  p <- r$pooled[term == t]
  data.table(analysis = lab, term = t, or = exp(p$estimate), lower = exp(p$lower), upper = exp(p$upper),
             p = p$p, n = r$diag$nobs, events = r$diag$n_events)
}
fe1 <- {
  f <- update(lme4::nobars(build_formula("d1y", COVARS_M1, specs)), . ~ . + centrum)
  dd <- cc_rows(d1, f)
  cbind(data.table(analysis = "Centre fixed effects, CR2", term = "year_c"), glm_or(f, dd, "year_c", "CR2"))
}
sens1 <- rbindlist(list(
  run1("M1 primary", COVARS_M1),
  run1("M0: year + centre only", character(0)),
  run1("M2: M1 + operative complexity", c(COVARS_M1, COVARS_OPERATIVE)),
  run1("M3: year + logit EuroSCORE II only", "logit_es2"),
  run1("Extended adjustment", COVARS_EXTENDED),
  run1("Age, sex, redo without EuroSCORE II", c("age", "female", "prior_cs")),
  run1("Random intercept + random slope for year", COVARS_M1, rnd = "(1 + year_c | centrum)"),
  fe1,
  run1(paste("Operations from", CONFIG$trend_sens_start), COVARS_M1, dd = d1[year >= CONFIG$trend_sens_start]),
  run1("Excluding salvage operations", COVARS_M1, dd = d1[urg4 != "Salvage"]),
  run1("Acute vs non-acute (M1)", c(COVARS_M1, "acute"), t = "acute")), fill = TRUE)
write_out(sens1, SUB, "sensitivity_one_year_mortality")

# ---- B. Missing data by year and by centre --------------------------------------------------------
log_msg("B. Missing data by year and centre")
mb <- function(by) cohort[, .(n = .N,
  pct_missing_euroscore2 = round(100 * mean(is.na(euroscore2)), 1),
  pct_missing_temp_min = round(100 * mean(is.na(temp_min)), 1),
  pct_missing_pathology = round(100 * mean(is.na(path3)), 1),
  pct_missing_crossclamp = round(100 * mean(is.na(aox_min)), 1),
  pct_missing_arrest_duration_if_arrest = round(100 * mean(is.na(ca_min[circ_arrest_b %in% 1])), 1)), by = by]
write_out(mb("year")[order(year)], SUB, "missing_by_year")
mbc <- mb("centre_lab")[order(centre_lab)]
write_out(mbc, SUB, "missing_by_centre_anonymised")

# ---- Paper 2: volume -----------------------------------------------------------------------------------
log_msg("Paper 2: volume analyses")
vg <- cohort[!is.na(urgency2), .N, by = .(centrum = as.character(centrum), year, urgency2)]
vw <- dcast(vg, centrum + year ~ urgency2, value.var = "N", fill = 0)
setnames(vw, c("Acute", "Non-acute"), c("vol_acute", "vol_nonacute"), skip_absent = TRUE)
grid <- CJ(centrum = unique(vw$centrum), year = CONFIG$study_years)
vw <- merge(grid, vw, by = c("centrum", "year"), all.x = TRUE)
for (v in c("vol_acute", "vol_nonacute")) set(vw, which(is.na(vw[[v]])), v, 0L)
setorder(vw, centrum, year)
vw[, `:=`(prev_acute = shift(vol_acute), prev_nonacute = shift(vol_nonacute)), by = centrum]
dv <- merge(d0, vw, by = c("centrum", "year"), all.x = TRUE)
dv[, `:=`(vol = fifelse(urgency2 == "Acute", vol_acute, vol_nonacute),
          vol_prev = fifelse(urgency2 == "Acute", prev_acute, prev_nonacute))]
PERIODS <- list("2013-2024" = c(2013, 2024), "2018-2024" = c(2018, 2024))
rhs <- paste(c("year_c", m1), collapse = " + ")
vol_rows <- list()
addv <- function(lab, res, extra = list()) if (!is.null(res))
  vol_rows[[length(vol_rows) + 1]] <<- cbind(as.data.table(c(lab, extra)), res)
for (gr in c("Acute", "Non-acute")) for (pn in names(PERIODS)) {
  per <- PERIODS[[pn]]
  base <- dv[urgency2 == gr & year >= per[1] & year <= per[2]]
  for (v in c("dzh", "d30", "d1y")) {
    x <- base[!is.na(get(v))]
    if (v == "d1y") x <- x[year <= CONFIG$d1y_last_year]
    x <- cc_rows(x, as.formula(paste(v, "~", rhs)))
    if (sum(x[[v]] == 1) < 30) next
    lab <- list(group = gr, period = pn, outcome = v)
    x[, `:=`(vol10 = vol / 10, vprev10 = vol_prev / 10)]
    addv(c(lab, analysis = "Same-year volume (reference)"), glm_or(as.formula(paste(v, "~ vol10 +", rhs)), x, "vol10", "CR2"))
    xp <- x[!is.na(vprev10)]
    addv(c(lab, analysis = "Prior-year volume"), glm_or(as.formula(paste(v, "~ vprev10 +", rhs)), xp, "vprev10", "CR2"))
    x[, vmean10 := mean(vol10), by = centrum]; x[, vdev10 := vol10 - vmean10]
    fwb <- as.formula(paste(v, "~ vmean10 + vdev10 +", rhs))
    addv(c(lab, analysis = "Between centres (centre mean volume)"), glm_or(fwb, x, "vmean10", "CR2"))
    addv(c(lab, analysis = "Within centre (yearly deviation from centre mean)"), glm_or(fwb, x, "vdev10", "CR2"))
    if (gr == "Acute") {
      x[, `:=`(nonacute10 = vol_nonacute / 10, acute10 = vol_acute / 10)]
      addv(c(lab, analysis = "Elective (non-acute) volume"), glm_or(as.formula(paste(v, "~ nonacute10 +", rhs)), x, "nonacute10", "CR2"))
      addv(c(lab, analysis = "Elective volume, adjusted for acute volume"),
           glm_or(as.formula(paste(v, "~ nonacute10 + acute10 +", rhs)), x, "nonacute10", "CR2"))
    }
  }
  # F. complications and failure to rescue by volume
  x <- cc_rows(base[!is.na(major_compl)], as.formula(paste("major_compl ~", rhs)))
  x[, vol10 := vol / 10]
  addv(list(group = gr, period = pn, outcome = "major_compl", analysis = "Major complication by volume"),
       glm_or(as.formula(paste("major_compl ~ vol10 +", rhs)), x, "vol10", "CR2"))
  xf <- x[major_compl == 1 & !is.na(dzh)]
  if (sum(xf$dzh == 1) >= 20)
    addv(list(group = gr, period = pn, outcome = "ftr_dzh", analysis = "Failure to rescue by volume (adjusted for number of complications)"),
         glm_or(as.formula(paste("dzh ~ vol10 + n_compl +", rhs)), xf, "vol10", "CR2"))
  # G. long-term survival by volume (Cox, robust SE by centre)
  xs <- cc_rows(base[!is.na(time_d) & !is.na(status)], as.formula(paste("status ~", rhs)))
  xs[, vol10 := vol / 10]
  cs <- tryCatch(summary(coxph(as.formula(paste("Surv(time_d, status) ~ vol10 +", rhs, "+ cluster(centrum)")), xs)),
                 error = function(e) NULL)
  if (!is.null(cs))
    vol_rows[[length(vol_rows) + 1]] <- data.table(group = gr, period = pn, outcome = "overall survival",
      analysis = "Long-term survival by volume (Cox HR per 10, robust SE)", or = cs$conf.int["vol10", 1],
      lower = cs$conf.int["vol10", 3], upper = cs$conf.int["vol10", 4], p = cs$coefficients["vol10", "Pr(>|z|)"],
      n = cs$n, events = cs$nevent)
}
vol_extra <- rbindlist(vol_rows, fill = TRUE)
setnames(vol_extra, "or", "ratio_per_10")
write_out(vol_extra, SUB, "volume_additional_analyses")

# ---- Paper 3: practice and outcome ------------------------------------------------------------------
log_msg("Paper 3: practice-outcome analyses")
dp <- copy(d0)[, `:=`(temp_1C = temp_min, deep = as.integer(temp_min <= 20), arrest = circ_arrest_b,
                      aox_30 = aox_min / 30, ca_10 = fifelse(circ_arrest_b == 1 & ca_min > 0, ca_min / 10, NA_real_))]
PR <- c(temp_1C = "Minimum temperature (per 1 C warmer)", deep = "Deep hypothermia <=20 C (yes vs no)",
        arrest = "Circulatory arrest (yes vs no)", aox_30 = "Cross-clamp time (per 30 min)",
        ca_10 = "Circulatory arrest time (per 10 min; arrest patients)")
UIT <- c("comp_primary", "dzh", "d30", "d1y", "stroke_def", "renal", "major_compl")
OPS <- c(root = "root_inv", descending = "desc_inv", fet_hybrid = "hybrid_fet", valve = "conc_valve",
         cabg = "conc_cabg", other_cardiac = "other_cardiac")
ADJ <- c(list("M1" = character(0)), setNames(as.list(OPS), paste("M1 +", names(OPS))),
         list("M1 + all operative" = unname(OPS)))
rhs_p <- function(grp, extra) paste(c("year_c", m1, if (grp == "All") "acute", extra), collapse = " + ")
prep <- function(grp, p, v, dd = dp) {
  x <- dd[!is.na(get(v)) & !is.na(get(p))]
  if (v == "d1y") x <- x[year <= CONFIG$d1y_last_year]
  if (grp != "All") x <- x[urgency2 == grp]
  cc_rows(x, as.formula(paste(v, "~", rhs_p(grp, unname(OPS)))))
}
# H. operative complexity, separately and together (same patients in every column)
hrows <- list()
for (grp in c("All", "Acute", "Non-acute")) for (p in names(PR)) for (v in UIT) {
  x <- prep(grp, p, v)
  if (sum(x[[v]] == 1) < 20) next
  for (a in names(ADJ)) {
    r <- tryCatch(glm_or(as.formula(paste(v, "~", p, "+", rhs_p(grp, ADJ[[a]]))), x, p), error = function(e) NULL)
    if (!is.null(r)) hrows[[length(hrows) + 1]] <- cbind(data.table(group = grp, practice = PR[[p]], outcome = v, adjustment = a), r)
  }
}
hx <- rbindlist(hrows)
write_out(hx, SUB, "practice_outcome_operative_complexity")
hw <- dcast(hx[, .(group, practice, outcome, adjustment = factor(adjustment, names(ADJ)), est = fmt(or, lower, upper))],
            group + practice + outcome ~ adjustment, value.var = "est")
write_out(hw, SUB, "practice_outcome_operative_complexity_wide")

# I. excluding centres with > 50% missing nadir temperature
miss_c <- cohort[, .(miss = mean(is.na(temp_min))), by = .(centrum = as.character(centrum))]
keep <- miss_c[miss <= 0.5, centrum]
log_msg("  centres excluded for > 50% missing temperature: ", nrow(miss_c) - length(keep))
irows <- list()
for (grp in c("All", "Acute", "Non-acute")) for (p in names(PR)) for (v in UIT) {
  x <- prep(grp, p, v, dp[centrum %in% keep])
  if (sum(x[[v]] == 1) < 20) next
  for (a in c("M1", "M1 + all operative")) {
    r <- tryCatch(glm_or(as.formula(paste(v, "~", p, "+", rhs_p(grp, ADJ[[a]]))), x, p), error = function(e) NULL)
    if (!is.null(r)) irows[[length(irows) + 1]] <- cbind(data.table(group = grp, practice = PR[[p]], outcome = v, adjustment = a,
                                                                    centres = length(unique(x$centrum))), r)
  }
}
write_out(rbindlist(irows), SUB, "practice_outcome_excluding_centres_missing_temperature")

# J. shape of the arrest-duration association (arrest patients)
rcs4 <- function(x, kn) { k <- length(kn); nrm <- (kn[k] - kn[1])^2
  cbind(x, sapply(seq_len(k - 2), function(j) (pmax(x - kn[j], 0)^3 - pmax(x - kn[k - 1], 0)^3 * (kn[k] - kn[j]) / (kn[k] - kn[k - 1]) +
                                                pmax(x - kn[k], 0)^3 * (kn[k - 1] - kn[j]) / (kn[k] - kn[k - 1])) / nrm)) }
jcat <- list(); jspl <- list(); jcurve <- list()
CAT <- c(0, 20, 30, 40, 60, Inf); CATL <- c("<20", "20-29", "30-39", "40-59", ">=60")
for (grp in c("All", "Acute", "Non-acute")) for (v in c("comp_primary", "dzh", "d30", "d1y", "renal", "stroke_def")) {
  x <- prep(grp, "ca_10", v)
  if (sum(x[[v]] == 1) < 30) next
  x[, ca := ca_10 * 10]
  x[, ca_cat := relevel(factor(as.character(cut(ca, CAT, CATL, right = FALSE)), levels = CATL), "<20")]
  rr <- rhs_p(grp, unname(OPS))
  m <- glm(as.formula(paste(v, "~ ca_cat +", rr)), binomial, x)
  V <- sandwich::vcovCL(m, cluster = x$centrum, type = "HC1")
  for (lv in CATL[-1]) { tn <- paste0("ca_cat", lv); if (!tn %in% names(coef(m)) || is.na(coef(m)[tn])) next
    b <- coef(m)[[tn]]; se <- sqrt(V[tn, tn])
    jcat[[length(jcat) + 1]] <- data.table(group = grp, outcome = v, category = lv, n_in_category = sum(x$ca_cat == lv),
      events_in_category = sum(x[ca_cat == lv][[v]] == 1), or_vs_below20 = exp(b), lower = exp(b - qt(.975, G - 1) * se),
      upper = exp(b + qt(.975, G - 1) * se), p = 2 * pt(-abs(b / se), G - 1)) }
  kn <- as.numeric(quantile(x$ca, c(.05, .35, .65, .95)))
  B <- rcs4(x$ca, kn); colnames(B) <- paste0("cs", seq_len(ncol(B))); x[, (colnames(B)) := as.data.table(B)]
  ms <- glm(as.formula(paste(v, "~", paste(colnames(B), collapse = " + "), "+", rr)), binomial, x)
  Vs <- sandwich::vcovCL(ms, cluster = x$centrum, type = "HC1")
  nl <- colnames(B)[-1]; bn <- coef(ms)[nl]
  wald <- as.numeric(t(bn) %*% solve(Vs[nl, nl]) %*% bn) / length(nl)
  jspl[[length(jspl) + 1]] <- data.table(group = grp, outcome = v, n = nrow(x), events = sum(x[[v]] == 1),
    knots_min = paste(round(kn), collapse = "/"), p_nonlinear = pf(wald, length(nl), G - 1, lower.tail = FALSE))
  for (g in seq(10, 90, by = 5)) { nd <- copy(x); Bg <- rcs4(rep(g, nrow(nd)), kn); nd[, (colnames(B)) := as.data.table(Bg)]
    jcurve[[length(jcurve) + 1]] <- data.table(group = grp, outcome = v, arrest_min = g,
      standardised_risk = mean(predict(ms, nd, type = "response"))) }
}
write_out(rbindlist(jcat), SUB, "arrest_duration_categories")
write_out(rbindlist(jspl), SUB, "arrest_duration_spline_tests")
jc <- rbindlist(jcurve)
write_out(jc, SUB, "arrest_duration_standardised_risk")
if (nrow(jc)) {
  pj <- ggplot(jc[outcome %in% c("dzh", "renal", "comp_primary")], aes(arrest_min, 100 * standardised_risk, colour = group)) +
    geom_line(linewidth = .8) + facet_wrap(~ outcome, scales = "free_y") + scale_colour_manual(values = PAL[1:3], name = NULL) +
    labs(x = "Circulatory arrest time (min)", y = "Standardised risk (%)", title = "Circulatory arrest time and outcome") + theme_pub()
  save_fig(pj, SUB, "arrest_duration_standardised_risk", 9, 4)
}

# K. arrest duration x nadir temperature
krows <- list()
for (grp in c("All", "Acute", "Non-acute")) for (v in c("comp_primary", "dzh", "d30", "renal", "stroke_def", "major_compl")) {
  x <- prep(grp, "ca_10", v)
  x <- x[!is.na(temp_min)]
  if (sum(x[[v]] == 1) < 30) next
  x[, temp_c := temp_min - 24]
  m <- glm(as.formula(paste(v, "~ ca_10 * temp_c +", rhs_p(grp, unname(OPS)))), binomial, x)
  V <- sandwich::vcovCL(m, cluster = x$centrum, type = "HC1")
  bi <- coef(m)[["ca_10:temp_c"]]; sei <- sqrt(V["ca_10:temp_c", "ca_10:temp_c"])
  for (tt in c(20, 24, 28)) {
    L <- setNames(rep(0, length(coef(m))), names(coef(m))); L["ca_10"] <- 1; L["ca_10:temp_c"] <- tt - 24
    b <- sum(L * coef(m)); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    krows[[length(krows) + 1]] <- data.table(group = grp, outcome = v, n = nrow(x), events = sum(x[[v]] == 1),
      nadir_temperature = tt, or_per_10min_arrest = exp(b), lower = exp(b - qt(.975, G - 1) * se),
      upper = exp(b + qt(.975, G - 1) * se), interaction_ratio_per_1C = exp(bi),
      p_interaction = 2 * pt(-abs(bi / sei), G - 1))
  }
}
write_out(rbindlist(krows), SUB, "arrest_duration_by_temperature_interaction")

# L. renal failure as intermediate between arrest duration and death
lrows <- list()
for (grp in c("All", "Acute", "Non-acute")) for (v in c("dzh", "d30", "d1y")) {
  x <- prep(grp, "ca_10", v)
  x <- x[!is.na(renal)]
  if (sum(x[[v]] == 1) < 30) next
  rr <- rhs_p(grp, unname(OPS))
  a  <- glm_or(as.formula(paste(v, "~ ca_10 +", rr)), x, "ca_10")
  b  <- glm_or(as.formula(paste(v, "~ ca_10 + renal +", rr)), x, "ca_10")
  rn <- glm_or(as.formula(paste("renal ~ ca_10 +", rr)), x, "ca_10")
  rd <- glm_or(as.formula(paste(v, "~ renal + ca_10 +", rr)), x, "renal")
  lrows[[length(lrows) + 1]] <- data.table(group = grp, outcome = v, n = nrow(x),
    or_arrest_on_renal = rn$or, or_renal_on_death = rd$or,
    or_arrest_on_death = a$or, lower = a$lower, upper = a$upper,
    or_arrest_on_death_adj_renal = b$or, lower_adj = b$lower, upper_adj = b$upper,
    pct_attenuation_log_scale = 100 * (1 - log(b$or) / log(a$or)))
}
write_out(rbindlist(lrows), SUB, "arrest_renal_death_attenuation")

log_msg("Done.")

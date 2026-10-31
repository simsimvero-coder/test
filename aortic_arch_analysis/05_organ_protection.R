# =============================================================================
# 05_organ_protection.R
# Temperature management and circulatory arrest: temporal trends, between-centre
# variation and within-centre trends (random intercept + slope models; centre
# fixed-effects models with CR2 small-cluster inference).
# Output: results/05_organ_protection/*, data/derived/organ.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
suppressPackageStartupMessages(library(ordinal))
start_log("05_organ_protection")
log_msg("Missing-data method: ", method_label())
SUB <- "05_organ_protection"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()   # imputed data sets or observed cohort (complete case)
specs  <- load_derived("specs.rds")
CASEMIX <- COVARS_M1   # adjustment: age, sex, EuroSCORE II, previous cardiac surgery

# Population: patients operated with cardiopulmonary bypass
base <- function(d) d[ecc_any == 1 | is.na(ecc_any)]
write_out(data.table(group = c("cohort", "with CPB (or CPB unknown)", "CPB unknown", "no CPB"),
                     n = c(nrow(cohort), nrow(base(cohort)), sum(is.na(cohort$ecc_any)),
                           sum(cohort$ecc_any == 0, na.rm = TRUE))), SUB, "population")
coh <- base(cohort)

# ---- 1. Descriptives (Table 2 data) ---------------------------------------------------------
desc_by <- function(d, by) {
  d[, .(n = .N,
        arrest = fmt_npct(circ_arrest_b),
        temp_median = fmt_med(temp_min),
        deep = fmt_npct(as.integer(hypo4 == "Deep (<=20)")),
        moderate = fmt_npct(as.integer(hypo4 == "Moderate (20.1-28)")),
        mild = fmt_npct(as.integer(hypo4 == "Mild (28.1-34)")),
        normothermia = fmt_npct(as.integer(hypo4 == "Normothermia (>34)")),
        temp_missing = sum(is.na(temp_min)),
        arrest_min_median = fmt_med(ca_min[circ_arrest_b == 1]),
        arrest_gt150 = sum(ca_min > 150 & circ_arrest_b == 1, na.rm = TRUE),
        cpb_min_median = fmt_med(ecc_min), aox_min_median = fmt_med(aox_min)), by = by]
}
t2 <- rbind(cbind(stratum = "All", desc_by(coh, "era")),
            cbind(stratum = "Non-acute", desc_by(coh[urgency2 == "Non-acute"], "era")),
            cbind(stratum = "Acute", desc_by(coh[urgency2 == "Acute"], "era")),
            cbind(stratum = "All", era = "Total", desc_by(coh, NULL)), fill = TRUE)
write_out(t2, SUB, "table2_organ_protection_by_era")
write_out(coh[, .N, by = .(year, hypo5)][order(year, hypo5)], SUB, "hypothermia_consensus_5level_by_year")
annual <- coh[, .(n = .N, arrest = mean(circ_arrest_b, na.rm = TRUE),
                  temp_median = as.numeric(median(temp_min, na.rm = TRUE)),
                  temp_q1 = as.numeric(quantile(temp_min, .25, na.rm = TRUE)), temp_q3 = as.numeric(quantile(temp_min, .75, na.rm = TRUE)),
                  ca_median = as.numeric(median(ca_min[circ_arrest_b == 1], na.rm = TRUE))), by = year][order(year)]
write_out(annual, SUB, "annual_practice")
hypo_year <- coh[!is.na(hypo4), .N, by = .(year, hypo4)][, pct := 100 * N / sum(N), by = year]
write_out(hypo_year[order(year, hypo4)], SUB, "hypothermia_category_by_year")
centre_year <- coh[, .(n = .N, temp_median = as.numeric(median(temp_min, na.rm = TRUE)),
                       arrest = mean(circ_arrest_b, na.rm = TRUE)), by = .(centre_lab, year)]
write_out(centre_year[order(centre_lab, year)], SUB, "centre_year_practice")

# ---- helpers ------------------------------------------------------------------------------
rhs_terms <- function(adjust) if (adjust) covar_terms(CASEMIX, specs) else character(0)
mk_formula <- function(y, adjust, random) {
  as.formula(paste(y, "~", paste(c("year_c", rhs_terms(adjust), random), collapse = " + ")), env = globalenv())
}
lrt_mixture <- function(m0, m1) {
  lr <- max(0, 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))))
  data.table(lr = lr, p_mixture = 0.5 * pchisq(lr, 1, lower.tail = FALSE) + 0.5 * pchisq(lr, 2, lower.tail = FALSE))
}
vc_table <- function(fit) {
  v <- as.data.frame(VarCorr(fit)); setDT(v)
  v[, .(component = paste(grp, var1, ifelse(is.na(var2), "", var2)), sd_or_cor = sdcor)]
}
fe_cr2 <- function(y, d, adjust, family = NULL) {
  f <- as.formula(paste(y, "~", paste(c("year_c", rhs_terms(adjust), "centrum"), collapse = " + ")))
  d <- complete_for(list(d), f)[[1]]
  fit <- if (is.null(family)) lm(f, data = d) else glm(f, family = family, data = d)
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = d$centrum, test = "Satterthwaite")
  r <- ct[ct$Coef == "year_c", ]
  data.table(estimate = r$beta, se = r$SE, df = r$df_Satt,
             lower = r$beta - qt(.975, r$df_Satt) * r$SE, upper = r$beta + qt(.975, r$df_Satt) * r$SE, p = r$p_Satt)
}

practice <- list(); vcs <- list(); slopes <- list(); lrts <- list(); fe <- list(); cents <- list()

# ---- 2. Continuous outcomes: minimum temperature, arrest duration ----------------------------
cont_specs <- list(
  list(y = "temp_min", label = "Minimum temperature (C)", sub = function(d) d[!is.na(temp_min)], pop = "CPB", log = FALSE),
  list(y = "temp_min", label = "Minimum temperature (C)", sub = function(d) d[!is.na(temp_min) & circ_arrest_b == 1], pop = "circulatory arrest", log = FALSE),
  list(y = "log_ca", label = "Circulatory arrest duration (ratio per year)", sub = function(d) d[circ_arrest_b == 1 & ca_min > 0][, log_ca := log(ca_min)], pop = "circulatory arrest", log = TRUE))
for (cs in cont_specs) {
  for (adjust in c(FALSE, TRUE)) {
    key <- paste(cs$y, cs$pop, if (adjust) "adjusted" else "unadjusted", sep = " | ")
    log_msg(key)
    dl <- lapply(base_list <- lapply(imps, base), cs$sub)
    f_int <- mk_formula(cs$y, adjust, "(1 | centrum)")
    f_slp <- mk_formula(cs$y, adjust, "(1 + year_c | centrum)")
    dl <- complete_for(dl, f_slp)
    p_slp <- fit_pooled(f_slp, dl, family = NULL)
    est <- p_slp$pooled[term == "year_c"]
    if (cs$log) est[, `:=`(estimate = exp(estimate), lower = exp(lower), upper = exp(upper))]
    practice[[key]] <- cbind(outcome = cs$label, population = cs$pop, adjusted = adjust,
                             model = "random intercept + slope", period = "all years", n = nrow(dl[[1]]),
                             est[, .(estimate, lower, upper, p)], p_slp$diag[, .(n_convergence_warning, n_singular)])
    # Sensitivity: from CONFIG$trend_sens_start onwards
    dls <- lapply(dl, function(d) d[year >= CONFIG$trend_sens_start])
    p_s <- fit_pooled(f_slp, dls, family = NULL)
    es <- p_s$pooled[term == "year_c"]
    if (cs$log) es[, `:=`(estimate = exp(estimate), lower = exp(lower), upper = exp(upper))]
    practice[[paste(key, "from")]] <- cbind(outcome = cs$label, population = cs$pop, adjusted = adjust,
                             model = "random intercept + slope", period = paste("from", CONFIG$trend_sens_start),
                             n = nrow(dls[[1]]), es[, .(estimate, lower, upper, p)],
                             p_s$diag[, .(n_convergence_warning, n_singular)])
    # variance components, slope LRT, centre slopes and FE model (imputation 1)
    m_slp <- p_slp$fits[[1]]$fit
    m_int <- fit_mixed(f_int, dl[[1]])$fit
    vcs[[key]] <- cbind(analysis = key, vc_table(m_slp))
    icc <- { v <- as.data.frame(VarCorr(m_int)); v$vcov[1] / sum(v$vcov) }
    lrts[[key]] <- cbind(analysis = key, lrt_mixture(m_int, m_slp), icc_intercept_model = icc)
    re <- as.data.frame(ranef(m_slp, condVar = TRUE)); setDT(re)
    key_c <- unique(cohort[, .(grp = as.character(centrum), centre_lab)])
    sl <- merge(re[term == "year_c", .(grp = as.character(grp), dev = condval, se = condsd)], key_c, by = "grp")
    sl[, `:=`(key = key, slope = fixef(m_slp)["year_c"] + dev)]
    slopes[[key]] <- sl[, grp := NULL]
    fr <- fe_cr2(cs$y, dl[[1]], adjust)
    if (cs$log) fr[, `:=`(estimate = exp(estimate), lower = exp(lower), upper = exp(upper))]
    fe[[key]] <- cbind(analysis = key, model = "centre fixed effects, CR2", fr)
  }
}

# ---- 3. Binary: use of circulatory arrest -----------------------------------------------------
for (adjust in c(FALSE, TRUE)) {
  key <- paste("circ_arrest_b | CPB |", if (adjust) "adjusted" else "unadjusted")
  log_msg(key)
  dl <- lapply(imps, function(d) base(d)[!is.na(circ_arrest_b)])
  f_int <- mk_formula("circ_arrest_b", adjust, "(1 | centrum)")
  f_slp <- mk_formula("circ_arrest_b", adjust, "(1 + year_c | centrum)")
  dl <- complete_for(dl, f_slp)
  p_int <- fit_pooled(f_int, dl)
  p_slp <- fit_pooled(f_slp, dl)
  for (nm in c("random intercept", "random intercept + slope")) {
    pp <- if (nm == "random intercept") p_int else p_slp
    e <- exp_table(pp$pooled, "year_c")
    practice[[paste(key, nm)]] <- data.table(outcome = "Circulatory arrest use (OR per year)", population = "CPB",
                                             adjusted = adjust, model = nm, period = "all years", n = nrow(dl[[1]]),
                                             estimate = e$ratio, lower = e$ratio_lower, upper = e$ratio_upper,
                                             p = e$p, pp$diag[, .(n_convergence_warning, n_singular)])
  }
  dls <- lapply(dl, function(d) d[year >= CONFIG$trend_sens_start])
  ps <- fit_pooled(f_slp, dls); es <- exp_table(ps$pooled, "year_c")
  practice[[paste(key, "from")]] <- data.table(outcome = "Circulatory arrest use (OR per year)", population = "CPB",
                                               adjusted = adjust, model = "random intercept + slope",
                                               period = paste("from", CONFIG$trend_sens_start), n = nrow(dls[[1]]),
                                               estimate = es$ratio, lower = es$ratio_lower, upper = es$ratio_upper,
                                               p = es$p, ps$diag[, .(n_convergence_warning, n_singular)])
  m_int <- p_int$fits[[1]]$fit; m_slp <- p_slp$fits[[1]]$fit
  vcs[[key]] <- cbind(analysis = key, vc_table(m_slp))
  s <- attr(VarCorr(m_int)$centrum, "stddev")
  lrts[[key]] <- cbind(analysis = key, lrt_mixture(m_int, m_slp), icc_intercept_model = icc_logit(s), mor = mor_from_sd(s))
  fr <- fe_cr2("circ_arrest_b", dl[[1]], adjust, binomial())
  fe[[key]] <- cbind(analysis = key, model = "centre fixed effects, CR2 (OR)", fr[, `:=`(estimate = exp(estimate), lower = exp(lower), upper = exp(upper))])
  re <- as.data.frame(ranef(m_int, condVar = TRUE)); setDT(re)
  key_c <- unique(cohort[, .(grp = as.character(centrum), centre_lab)])
  cents[[key]] <- merge(re[, .(grp = as.character(grp), effect = condval, se = condsd)], key_c, by = "grp")[, key := key][, grp := NULL]
}

# ---- 4. Ordinal: hypothermia category -------------------------------------------------------------
# Ordinal mixed model with a usable covariance matrix. If the Hessian is not positive
# definite, refit with other optimisers; an imputation that still fails is left out
# of the Rubin pooling (counted and written to ordinal_fit_problems.csv).
ord_problems <- list()
fit_clmm_safe <- function(f, d) {
  fx <- function(m) names(coef(m))[!names(coef(m)) %in% names(m$alpha)]
  tries <- list(list(), list(control = clmm.control(method = "ucminf")),
                list(control = clmm.control(method = "nlminb", maxIter = 200, gradTol = 1e-6, innerCtrl = "noWarn")))
  for (k in seq_along(tries)) {
    m <- tryCatch(do.call(clmm, c(list(formula = f, data = d, link = "logit", nAGQ = 1), tries[[k]])),
                  error = function(e) NULL)
    if (is.null(m)) next
    v <- tryCatch(vcov(m), error = function(e) NULL)
    if (!is.null(v) && all(is.finite(diag(v)[fx(m)])) && all(diag(v)[fx(m)] > 0))
      return(list(fit = m, est = coef(m)[fx(m)], vc = v[fx(m), fx(m)], attempt = k))
  }
  NULL
}
pool_clmm <- function(fs, label) {
  ok <- !vapply(fs, is.null, logical(1))
  ord_problems[[label]] <<- data.table(model = label, imputations = length(fs), failed = sum(!ok),
                                       refitted = sum(vapply(fs[ok], function(x) x$attempt > 1, logical(1))))
  if (!any(ok)) return(NULL)
  if (any(!ok)) log_msg(label, ": ", sum(!ok), " of ", length(fs), " fits without usable covariance - left out of pooling")
  list(pooled = exp_table(rubin_pool(lapply(fs[ok], `[[`, "est"), lapply(fs[ok], `[[`, "vc")), "year_c"),
       fits = lapply(fs[ok], `[[`, "fit"))
}

ord_res <- list()
for (adjust in c(FALSE, TRUE)) {
  log_msg("hypothermia ordinal, adjusted = ", adjust)
  dl <- lapply(imps, function(d) base(d)[!is.na(hypo4)][, hypo4 := droplevels(hypo4)])
  f <- as.formula(paste("hypo4 ~", paste(c("year_c", rhs_terms(adjust), "(1 | centrum)"), collapse = " + ")))
  dl <- complete_for(dl, f)
  pr <- pool_clmm(par_lapply(dl, function(d) fit_clmm_safe(f, d)),
                  paste0("hypothermia ordinal, ", if (adjust) "adjusted" else "unadjusted", ", all years"))
  if (is.null(pr)) { log_msg("hypothermia ordinal: no usable fit, skipped"); next }
  pooled <- pr$pooled; fits <- pr$fits
  if (adjust) {
    dls <- lapply(dl, function(d) d[year >= CONFIG$trend_sens_start][, hypo4 := droplevels(hypo4)])
    ps2 <- pool_clmm(par_lapply(dls, function(d) fit_clmm_safe(f, d)),
                     paste0("hypothermia ordinal, adjusted, from ", CONFIG$trend_sens_start))
    if (!is.null(ps2)) { es <- ps2$pooled
    ord_res[[length(ord_res) + 1]] <- data.table(
      outcome = "Hypothermia category (cumulative OR per year; >1 = shift to warmer)", adjusted = TRUE,
      period = paste("from", CONFIG$trend_sens_start), n = nrow(dls[[1]]),
      estimate = es$ratio, lower = es$ratio_lower, upper = es$ratio_upper, p = es$p) }
  }
  ord_res[[length(ord_res) + 1]] <- data.table(
    outcome = "Hypothermia category (cumulative OR per year; >1 = shift to warmer)",
    adjusted = adjust, period = "all years", n = nrow(dl[[1]]), estimate = pooled$ratio, lower = pooled$ratio_lower,
    upper = pooled$ratio_upper, p = pooled$p,
    centre_sd = mean(vapply(fits, function(m) sqrt(as.numeric(VarCorr(m)$centrum)), numeric(1))))
  # proportional-odds check on the fixed-effects clm (imputation 1)
  f_clm <- as.formula(paste("hypo4 ~", paste(c("year_c", rhs_terms(adjust), "centrum"), collapse = " + ")))
  m_clm <- clm(f_clm, data = dl[[1]])
  nt <- tryCatch(nominal_test(m_clm, scope = ~ year_c), error = function(e) NULL)
  if (!is.null(nt)) write_out(cbind(term = rownames(nt), as.data.table(nt)), SUB,
                              paste0("proportional_odds_test_", if (adjust) "adjusted" else "unadjusted"))
}
# Sensitivity for PO: binary deep hypothermia (<=20 C) with mixed logistic model
dl <- lapply(imps, function(d) base(d)[!is.na(hypo4)][, deep := as.integer(hypo4 == "Deep (<=20)")])
dl <- complete_for(dl, mk_formula("deep", TRUE, "(1 | centrum)"))
pd <- exp_table(fit_pooled(mk_formula("deep", TRUE, "(1 | centrum)"), dl)$pooled, "year_c")
ord_res[[length(ord_res) + 1]] <- data.table(outcome = "Deep hypothermia <=20 C (OR per year, binary sensitivity)",
                                             adjusted = TRUE, period = "all years", n = nrow(dl[[1]]), estimate = pd$ratio,
                                             lower = pd$ratio_lower, upper = pd$ratio_upper, p = pd$p)

practice_tab <- rbindlist(practice, fill = TRUE)
write_out(practice_tab, SUB, "practice_trend_models")
write_out(rbindlist(ord_res, fill = TRUE), SUB, "hypothermia_ordinal_models")
if (length(ord_problems)) write_out(rbindlist(ord_problems), SUB, "ordinal_fit_problems")
write_out(rbindlist(vcs), SUB, "variance_components")
write_out(rbindlist(lrts, fill = TRUE), SUB, "random_slope_lrt_and_icc")
write_out(rbindlist(slopes), SUB, "centre_specific_slopes")
write_out(rbindlist(fe), SUB, "within_centre_fixed_effects_cr2")
write_out(rbindlist(cents), SUB, "centre_effects_arrest")

saveRDS(list(t2 = t2, annual = annual, hypo_year = hypo_year, centre_year = centre_year,
             practice = practice_tab, ordinal = rbindlist(ord_res, fill = TRUE),
             slopes = rbindlist(slopes), fe = rbindlist(fe), lrt = rbindlist(lrts, fill = TRUE)),
        derived_method("organ.rds"))
log_msg("Done.")

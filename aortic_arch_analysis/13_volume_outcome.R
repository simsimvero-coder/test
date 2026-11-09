# =============================================================================
# 13_volume_outcome.R
# Volume-outcome analysis: is a higher annual centre volume associated with lower
# mortality, and above which annual volume is the risk no longer lower?
# Separately for acute and non-acute patients (exposure: the centre's annual
# number of operations in that group), for 2013-2024 and 2018-2024.
#
# Exposure: centre-year volume (operations of the same group in that centre in
#   that calendar year). Sensitivity: total centre-year volume (all operations).
# Outcomes: in-hospital mortality (primary), 30-day mortality, primary composite.
# Adjustment: age (spline), sex, logit EuroSCORE II, previous cardiac surgery, year.
# Models (logistic, patients nested in centres):
#   a. linear in volume (OR per 10 operations/year), CR2 cluster-robust inference
#      by centre (volume effect includes between-centre differences);
#   b. restricted cubic spline in volume, with a Wald test of non-linearity and
#      case-mix-standardised risk across the observed volume range;
#   c. threshold ("hinge") model: risk changes linearly up to a volume k and is
#      flat above k. k is chosen by profile likelihood over the 10th-90th
#      percentile; its 95% interval contains every k within 3.84 deviance units.
#      A threshold is only reported as such if the hinge model fits better than a
#      model without volume (likelihood-ratio test).
#      The likelihood-ratio p-value and the profile interval treat patients as
#      independent and k as known, so they are anti-conservative (volume varies over
#      15 centres only). They are reported as descriptive; inference on the threshold
#      comes from (e).
#   d. sensitivity: random intercept for centre instead of CR2 (linear).
#   e. centre-level cluster bootstrap of the hinge model (whole centres resampled with
#      replacement, k re-estimated in every replicate) for in-hospital and 30-day
#      mortality with same-group volume: percentile intervals for k and for the OR per
#      10 operations below k, and how often the hinge beats a straight line.
# Exploratory; 15 centres, so centre-level volume varies over few clusters.
# Output: results/<method>/13_volume_outcome/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("13_volume_outcome")
log_msg("Missing-data method: ", method_label())
SUB <- "13_volume_outcome"

cohort  <- load_derived("cohort.rds")
surg    <- load_derived("surgical_cohort.rds")
imps    <- analysis_datasets()
specs   <- load_derived("specs.rds")
d0      <- imps[[1]]                       # complete case: observed cohort
OUTS    <- c(dzh = "In-hospital mortality", d30 = "30-day mortality", comp_primary = "Primary composite")
PERIODS <- list("2013-2024" = range(CONFIG$study_years), "2018-2024" = c(2018, max(CONFIG$study_years)))
rhs     <- paste(c("year_c", covar_terms(COVARS_M1, specs)), collapse = " + ")

# centre-year volumes
vg <- cohort[!is.na(urgency2), .N, by = .(centrum = as.character(centrum), year, urgency2)]
vt <- surg[, .(vol_total = .N), by = .(centrum = as.character(centrum), year)]
d0[, centrum_chr := as.character(centrum)]
d0 <- merge(d0, vg[, .(centrum_chr = centrum, year, urgency2, vol_group = N)], by = c("centrum_chr", "year", "urgency2"), all.x = TRUE)
d0 <- merge(d0, vt[, .(centrum_chr = centrum, year, vol_total)], by = c("centrum_chr", "year"), all.x = TRUE)
write_out(dcast(vg, centrum + year ~ urgency2, value.var = "N", fill = 0)[order(centrum, year)], SUB, "centre_year_volume")

cr2_term <- function(m, d, term) {
  r <- clubSandwich::coef_test(m, vcov = "CR2", cluster = d$centrum, test = "Satterthwaite")
  r <- r[r$Coef == term, ]
  data.table(or = exp(r$beta), lower = exp(r$beta - qt(.975, r$df_Satt) * r$SE),
             upper = exp(r$beta + qt(.975, r$df_Satt) * r$SE), p = r$p_Satt)
}
rcs_basis <- function(x, knots) {           # Harrell restricted cubic spline (k knots -> k-1 columns)
  k <- length(knots); nrm <- (knots[k] - knots[1])^2
  b <- sapply(seq_len(k - 2), function(j) (pmax(x - knots[j], 0)^3 -
         pmax(x - knots[k - 1], 0)^3 * (knots[k] - knots[j]) / (knots[k] - knots[k - 1]) +
         pmax(x - knots[k], 0)^3 * (knots[k - 1] - knots[j]) / (knots[k] - knots[k - 1])) / nrm)
  cbind(x, matrix(b, ncol = k - 2))
}

lin <- list(); spl <- list(); thr <- list(); curve <- list(); prof <- list(); desc <- list()
for (gr in c("Acute", "Non-acute")) for (pn in names(PERIODS)) for (vv in c("vol_group", "vol_total")) for (v in names(OUTS)) {
  per <- PERIODS[[pn]]
  d <- d0[urgency2 == gr & year >= per[1] & year <= per[2] & !is.na(get(v))]
  d <- d[complete.cases(d[, .(age, female, logit_es2, prior_cs, year_c, centrum)])]
  d[, vol := get(vv)]
  ev <- sum(d[[v]] == 1)
  lab <- list(group = gr, period = pn, volume = if (vv == "vol_group") "same-group volume" else "total centre volume",
              outcome = v, n = nrow(d), events = ev)
  if (ev < 30) { log_msg("  ", paste(unlist(lab[1:4]), collapse = " / "), ": too few events (", ev, ")"); next }
  log_msg(paste(unlist(lab[1:4]), collapse = " / "))
  if (v == "dzh" && vv == "vol_group")
    desc[[length(desc) + 1]] <- cbind(as.data.table(lab[1:2]),
      d[, .(n = .N, deaths = sum(dzh == 1), mortality_pct = round(100 * mean(dzh), 1)),
        by = .(volume_band = cut(vol, c(0, 10, 20, 30, 40, 60, Inf), c("1-10", "11-20", "21-30", "31-40", "41-60", ">60")))][order(volume_band)])
  # a. linear (per 10 operations/year), CR2; d. random-intercept sensitivity
  d[, vol10 := vol / 10]
  m_lin <- glm(as.formula(paste(v, "~ vol10 +", rhs)), binomial, d)
  a <- cr2_term(m_lin, d, "vol10")
  ri <- tryCatch({ f <- fit_mixed(as.formula(paste(v, "~ vol10 +", rhs, "+ (1 | centrum)")), d, binomial())
                   e <- summary(f$fit)$coefficients["vol10", ]; c(exp(e[1]), exp(e[1] - 1.96 * e[2]), exp(e[1] + 1.96 * e[2]), e[4]) },
                 error = function(e) rep(NA_real_, 4))
  lin[[length(lin) + 1]] <- cbind(as.data.table(lab), vol_min = min(d$vol), vol_median = as.numeric(median(d$vol)), vol_max = max(d$vol),
                                  or_per_10 = a$or, lower = a$lower, upper = a$upper, p = a$p,
                                  or_per_10_random_intercept = ri[1], ri_lower = ri[2], ri_upper = ri[3], ri_p = ri[4])
  # b. spline (4 knots at 5/35/65/95th percentiles of patient-level volume)
  kn <- unique(quantile(d$vol, c(.05, .35, .65, .95)))
  if (length(kn) >= 4) {
    B <- rcs_basis(d$vol, kn); colnames(B) <- paste0("vs", seq_len(ncol(B))); d[, (colnames(B)) := as.data.table(B)]
    m_spl <- glm(as.formula(paste(v, "~", paste(colnames(B), collapse = " + "), "+", rhs)), binomial, d)
    V <- clubSandwich::vcovCR(m_spl, cluster = d$centrum, type = "CR2")
    nl <- colnames(B)[-1]
    w <- clubSandwich::Wald_test(m_spl, constraints = clubSandwich::constrain_zero(nl), vcov = V, test = "HTZ")
    wa <- clubSandwich::Wald_test(m_spl, constraints = clubSandwich::constrain_zero(colnames(B)), vcov = V, test = "HTZ")
    spl[[length(spl) + 1]] <- cbind(as.data.table(lab), knots = paste(round(kn), collapse = "/"),
                                    p_overall_volume = wa$p_val, p_nonlinear = w$p_val)
    grid <- seq(ceiling(min(d$vol)), floor(max(d$vol)), length.out = 25)
    std <- sapply(grid, function(g) { nd <- copy(d); Bg <- rcs_basis(rep(g, nrow(nd)), kn); nd[, (colnames(B)) := as.data.table(Bg)]
                                      mean(predict(m_spl, nd, type = "response")) })
    curve[[length(curve) + 1]] <- data.table(group = gr, period = pn, volume_type = lab$volume, outcome = v, volume = grid,
                                             standardised_risk = std)
  }
  # c. hinge (threshold) model
  ks <- unique(round(quantile(d$vol, seq(.10, .90, by = .02))))
  m0 <- glm(as.formula(paste(v, "~", rhs)), binomial, d)
  dev <- sapply(ks, function(k) deviance(glm(as.formula(paste(v, "~ pmin(vol, ", k, ") +", rhs)), binomial, d)))
  kbest <- ks[which.min(dev)]
  mh <- glm(as.formula(paste(v, "~ pmin(vol, ", kbest, ") +", rhs)), binomial, d)
  slope <- cr2_term(mh, d, grep("pmin", names(coef(mh)), value = TRUE))
  lr <- deviance(m0) - min(dev)
  prof[[length(prof) + 1]] <- data.table(group = gr, period = pn, volume_type = lab$volume, outcome = v, k = ks, deviance = dev)
  thr[[length(thr) + 1]] <- cbind(as.data.table(lab), threshold = kbest,
                                  threshold_ci_low = min(ks[dev <= min(dev) + 3.84]), threshold_ci_high = max(ks[dev <= min(dev) + 3.84]),
                                  or_per_10_below_threshold = slope$or^10, lower = slope$lower^10, upper = slope$upper^10,
                                  p_slope_cr2 = slope$p, lr_vs_no_volume = lr,
                                  p_lr_vs_no_volume = pchisq(lr, 2, lower.tail = FALSE),   # 2 df: slope + estimated threshold
                                  lr_vs_linear = deviance(m_lin) - min(dev),
                                  threshold_at_grid_edge = kbest %in% range(ks),
                                  threshold_supported = pchisq(lr, 2, lower.tail = FALSE) < 0.05 & slope$or < 1 &
                                    !(kbest %in% range(ks)) & (deviance(m_lin) - min(dev)) > 3.84)
}
write_out(rbindlist(desc), SUB, "mortality_by_volume_band")
write_out(rbindlist(lin), SUB, "volume_linear_models")
write_out(rbindlist(spl), SUB, "volume_spline_tests")
write_out(rbindlist(thr), SUB, "volume_threshold_models")
write_out(rbindlist(prof), SUB, "volume_threshold_profile")
cv <- rbindlist(curve)
write_out(cv, SUB, "volume_standardised_risk_curves")
if (nrow(cv)) {
  p <- ggplot(cv[volume_type == "same-group volume"], aes(volume, 100 * standardised_risk, colour = period)) +
    geom_line(linewidth = .8) + facet_grid(outcome ~ group, scales = "free") +
    scale_colour_manual(values = PAL[1:2], name = NULL) +
    labs(x = "Centre volume in the same group (operations per year)", y = "Case-mix-standardised risk (%)",
         title = "Volume-outcome relation (restricted cubic spline)") + theme_pub()
  save_fig(p, SUB, "volume_outcome_curves", 9, 8)
}
# ---- e. cluster bootstrap of the threshold (centres resampled) ---------------------------------
B <- as.integer(Sys.getenv("ARCH_VOLUME_BOOT", "200"))
set.seed(CONFIG$seed)
hinge_fit <- function(y, X, vol, ks) {
  dev <- vapply(ks, function(k) suppressWarnings(glm.fit(cbind(X, pmin(vol, k)), y, family = binomial()))$deviance, numeric(1))
  i <- which.min(dev)
  fh <- suppressWarnings(glm.fit(cbind(X, pmin(vol, ks[i])), y, family = binomial()))
  dl <- suppressWarnings(glm.fit(cbind(X, vol), y, family = binomial()))$deviance
  c(k = ks[i], or10 = exp(10 * tail(fh$coefficients, 1)), lr_vs_linear = dl - dev[i], edge = as.numeric(i %in% c(1, length(ks))))
}
boot <- list(); boot_reps <- list()
for (gr in c("Acute", "Non-acute")) for (pn in names(PERIODS)) for (v in c("dzh", "d30")) {
  per <- PERIODS[[pn]]
  d <- d0[urgency2 == gr & year >= per[1] & year <= per[2] & !is.na(get(v))]
  d <- d[complete.cases(d[, .(age, female, logit_es2, prior_cs, year_c, centrum)])]
  if (sum(d[[v]] == 1) < 30) next
  log_msg("  bootstrap ", gr, " / ", pn, " / ", v, " (", B, " replicates)")
  X <- model.matrix(as.formula(paste("~", rhs)), d); y <- d[[v]]; vol <- d$vol_group
  ks <- unique(round(quantile(vol, seq(.10, .90, by = .02))))
  est <- hinge_fit(y, X, vol, ks)
  idx <- split(seq_len(nrow(d)), as.character(d$centrum))
  reps <- t(vapply(seq_len(B), function(b) {
    rows <- unlist(idx[sample(names(idx), length(idx), replace = TRUE)], use.names = FALSE)
    out <- tryCatch(hinge_fit(y[rows], X[rows, , drop = FALSE], vol[rows], ks), error = function(e) rep(NA_real_, 4))
    out
  }, numeric(4)))
  colnames(reps) <- c("k", "or10", "lr_vs_linear", "edge")
  ok <- stats::complete.cases(reps) & is.finite(reps[, "or10"])
  r <- reps[ok, , drop = FALSE]
  boot_reps[[length(boot_reps) + 1]] <- data.table(group = gr, period = pn, outcome = v, r)
  boot[[length(boot) + 1]] <- data.table(group = gr, period = pn, volume = "same-group volume", outcome = v,
    n = nrow(d), events = sum(y == 1), centres = length(idx), replicates_ok = nrow(r),
    threshold = est[["k"]], threshold_boot_low = quantile(r[, "k"], .025), threshold_boot_high = quantile(r[, "k"], .975),
    or_per_10_below_threshold = est[["or10"]],
    or_boot_low = quantile(r[, "or10"], .025), or_boot_high = quantile(r[, "or10"], .975),
    share_reps_or_below_1 = mean(r[, "or10"] < 1),
    share_reps_threshold_at_grid_edge = mean(r[, "edge"] == 1),
    share_reps_hinge_better_than_linear = mean(r[, "lr_vs_linear"] > 3.84 & r[, "edge"] == 0))
}
write_out(rbindlist(boot), SUB, "volume_threshold_bootstrap")
write_out(rbindlist(boot_reps), SUB, "volume_threshold_bootstrap_replicates")

log_msg("Done.")

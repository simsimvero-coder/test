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

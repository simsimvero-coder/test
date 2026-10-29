# =============================================================================
# 07_population_survival.R
# Observed vs expected survival using Dutch period life tables (CBS 37360ned):
# expected survival with attained age and calendar year updated during
# follow-up, Ederer II relative survival, Pohar-Perme net survival (sensitivity),
# standardised mortality ratios, excess mortality rates and excess-hazard
# regression (Dickman et al. 2004). Includes explicit time-scale verification.
# Output: results/07_population/*, data/derived/population.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("07_population_survival")
SUB <- "07_population"

# ---- 1. Life table ----------------------------------------------------------------------------
read_cbs_lifetable <- function(dir) {
  obs_f <- file.path(dir, "Observations.csv")
  if (!file.exists(obs_f)) return(NULL)
  obs <- fread(obs_f, sep = ";", dec = ",", colClasses = list(character = c("Measure", "Geslacht", "LeeftijdOp31December", "Perioden")))
  ages <- fread(file.path(dir, "LeeftijdOp31DecemberCodes.csv"), sep = ";", colClasses = list(character = "Identifier"))
  ages[, age := as.integer(sub("^([0-9]+).*", "\\1", Title))]
  lt <- obs[Measure == "M005112" & Geslacht %in% c("3000", "4000") & grepl("JJ00$", Perioden)]
  lt <- merge(lt, ages[, .(LeeftijdOp31December = Identifier, age)], by = "LeeftijdOp31December")
  lt[, `:=`(year = as.integer(substr(Perioden, 1, 4)), female = as.integer(Geslacht == "4000"),
            qx = as.numeric(Value))]
  lt[, .(year, female, age, qx)][order(year, female, age)]
}
lt <- read_cbs_lifetable(CONFIG$cbs_lifetable_dir)
if (is.null(lt) || !nrow(lt)) {
  msg <- paste("Life table (CBS 37360ned) not found in", CONFIG$cbs_lifetable_dir,
               "- population-referenced analyses were NOT performed.")
  writeLines(msg, file.path(out_dir(SUB), "NOT_PERFORMED.txt")); log_msg(msg)
  quit(save = "no", status = 0)
}
stopifnot(all(lt$qx >= 0 & lt$qx < 1), !anyNA(lt$qx))
lt_years <- range(lt$year); lt_maxage <- max(lt$age)
write_out(lt[year >= min(CONFIG$study_years) - 1], SUB, "lifetable_qx_used")
log_msg("Life table: years ", lt_years[1], "-", lt_years[2], ", ages 0-", lt_maxage, "+")
# Lookup array [sex, age, year] of annual hazards h = -log(1 - q)
yrs_lt <- lt_years[1]:lt_years[2]
H <- array(NA_real_, c(2, lt_maxage + 1, length(yrs_lt)))
H[cbind(lt$female + 1, lt$age + 1, lt$year - lt_years[1] + 1)] <- -log(1 - lt$qx)
stopifnot(!anyNA(H))

# ---- 2. Patients -------------------------------------------------------------------------------
cohort <- load_derived("cohort.rds")
d <- cohort[!is.na(time_d) & !is.na(status) & !is.na(age) & !is.na(female)]
write_out(data.table(item = c("cohort", "used (follow-up, status, age and sex known)"), n = c(nrow(cohort), nrow(d))),
          SUB, "population_used")
d[, t_years := pmax(time_d, 0.5) / 365.25]
d[, age_grp := cut(age, c(-Inf, 59, 69, Inf), c("<60", "60-69", ">=70"))]
# Assumptions (no dates in extract): operation at mid-year; exact age = completed age + 0.5
d[, `:=`(age0 = age + 0.5, cal0 = year + 0.5)]
n_beyond <- sum(d$cal0 + d$t_years > lt_years[2] + 1)
log_msg("Follow-up beyond last life-table year (", lt_years[2], "): ", n_beyond,
        " patients; last available year carried forward")

# Expected cumulative hazard on a weekly grid (age and calendar year updated each step)
STEP <- 1 / 52
tmax <- max(d$t_years)
grid <- seq(0, ceiling(tmax / STEP) * STEP, by = STEP)
exp_cumhaz <- function(age0, cal0, female, grid) {
  n <- length(age0); L <- matrix(0, length(grid), n)
  for (k in seq_len(length(grid) - 1)) {
    mid <- (grid[k] + grid[k + 1]) / 2
    a <- pmin(pmax(floor(age0 + mid), 0), lt_maxage)
    cy <- pmin(pmax(floor(cal0 + mid), lt_years[1]), lt_years[2])
    h <- H[cbind(female + 1, a + 1, cy - lt_years[1] + 1)]
    L[k + 1, ] <- L[k, ] + h * (grid[k + 1] - grid[k])
  }
  L
}
L <- exp_cumhaz(d$age0, d$cal0, d$female, grid)
# Lambda*_i(t) by linear interpolation on the grid
lam_at <- function(i, t) {
  t <- pmin(t, max(grid)); k <- pmin(floor(t / STEP) + 1, length(grid) - 1)
  frac <- (t - grid[k]) / STEP
  L[cbind(k, i)] * (1 - frac) + L[cbind(k + 1, i)] * frac
}

# ---- 3. Time-scale verification ------------------------------------------------------------------
verify <- rbindlist(lapply(list(c(65, 0, 2013), c(65, 1, 2013), c(75, 0, 2018), c(55, 1, 2020)), function(p) {
  a <- p[1]; f <- p[2]; y <- p[3]
  Li <- exp_cumhaz(a + 0.5, y + 0.5, f, grid[grid <= 5 + 1e-9])
  s_pipeline <- exp(-Li[nrow(Li), 1])
  # hand calculation: half-year pieces at (age, year) combinations over 5 years
  pieces <- seq(0.25, 4.75, by = 0.5)
  hh <- vapply(pieces, function(m) H[f + 1, min(floor(a + 0.5 + m), lt_maxage) + 1,
                                    min(floor(y + 0.5 + m), lt_years[2]) - lt_years[1] + 1], numeric(1))
  data.table(profile = sprintf("%s aged %d, operated %d", ifelse(f == 1, "woman", "man"), a, y),
             expected_5y_survival_pipeline = s_pipeline, expected_5y_survival_hand = exp(-sum(hh * 0.5)),
             attained_age_at_5y = a + 0.5 + 5, calendar_year_at_5y = y + 0.5 + 5)
}))
verify[, abs_difference := abs(expected_5y_survival_pipeline - expected_5y_survival_hand)]
write_out(verify, SUB, "time_scale_verification")
stopifnot(all(verify$abs_difference < 0.002))

# ---- 4. SMR ----------------------------------------------------------------------------------------
d[, E_total := lam_at(seq_len(.N), t_years)]
smr_group <- function(dd, label) cbind(group = label, smr_ci(sum(dd$status), sum(dd$E_total)), person_years = sum(dd$t_years))
intervals <- list(c(0, 30 / 365.25, "0-30 days"), c(30 / 365.25, 1, "30 days-1 year"),
                  c(1, 5, "1-5 years"), c(5, Inf, ">5 years"))
smr_int <- function(dd, label) rbindlist(lapply(intervals, function(iv) {
  a <- as.numeric(iv[1]); b <- as.numeric(iv[2])
  idx <- which(dd$t_years > a)
  if (!length(idx)) return(NULL)
  ii <- dd$row[idx]; tt <- dd$t_years[idx]
  e <- sum(lam_at(ii, pmin(tt, b)) - lam_at(ii, rep(a, length(ii))))
  o <- sum(dd$status[idx] == 1 & tt <= b)
  py <- sum(pmin(tt, b) - a)
  cbind(group = label, interval = iv[3], smr_ci(o, e), person_years = py,
        excess_per_1000py = 1000 * (o - e) / py,
        excess_lower = 1000 * ((if (o == 0) 0 else qchisq(.025, 2 * o) / 2) - e) / py,
        excess_upper = 1000 * (qchisq(.975, 2 * (o + 1)) / 2 - e) / py)
}))
d[, row := .I]
smr_overall <- rbind(smr_group(d, "All"),
                     rbindlist(lapply(levels(d$urgency2), function(g) smr_group(d[urgency2 == g], paste("Group:", g)))),
                     rbindlist(lapply(c(0, 1), function(f) smr_group(d[female == f], ifelse(f == 1, "Women", "Men")))),
                     rbindlist(lapply(levels(d$age_grp), function(g) smr_group(d[age_grp == g], paste("Age", g)))),
                     rbindlist(lapply(levels(d$era), function(g) smr_group(d[era == g], paste("Era", g)))),
                     rbindlist(lapply(levels(d$path3), function(g) smr_group(d[path3 == g], paste("Pathology:", g)))))
write_out(smr_overall, SUB, "smr_overall_and_by_group")
smr_by_int <- rbind(smr_int(d, "All"),
                    rbindlist(lapply(levels(d$urgency2), function(g) smr_int(d[urgency2 == g], paste("Group:", g)))),
                    rbindlist(lapply(levels(d$age_grp), function(g) smr_int(d[age_grp == g], paste("Age", g)))))
write_out(smr_by_int, SUB, "smr_and_excess_by_followup_interval")

# ---- 5. Ederer II relative survival -----------------------------------------------------------------
breaks <- c(seq(0, 1, by = 1 / 12), 2:10)
breaks <- breaks[breaks <= max(d$t_years) + 1e-9]
ederer2 <- function(dd) {
  out <- list(); cumP <- 1; cumPs <- 1; varsum <- 0
  for (j in seq_len(length(breaks) - 1)) {
    a <- breaks[j]; b <- breaks[j + 1]
    at <- dd[t_years > a]
    if (nrow(at) < 10) break
    dj <- sum(at$status == 1 & at$t_years <= b); wj <- sum(at$status == 0 & at$t_years <= b)
    neff <- nrow(at) - wj / 2
    p <- 1 - dj / neff
    ps <- mean(exp(-(lam_at(at$row, rep(b, nrow(at))) - lam_at(at$row, rep(a, nrow(at))))))
    cumP <- cumP * p; cumPs <- cumPs * ps
    varsum <- varsum + if (neff > dj) dj / (neff * (neff - dj)) else 0
    R <- cumP / cumPs; se_log <- sqrt(varsum)
    out[[j]] <- data.table(start = a, end = b, n_at_risk = nrow(at), deaths = dj, censored = wj,
                           observed = cumP, expected = cumPs, relative = R,
                           lower = R * exp(-1.96 * se_log), upper = R * exp(1.96 * se_log))
  }
  rbindlist(out)
}
rs <- rbind(cbind(group = "All", ederer2(d)),
            rbindlist(lapply(levels(d$urgency2), function(g) cbind(group = g, ederer2(d[urgency2 == g])))),
            rbindlist(lapply(levels(d$age_grp), function(g) cbind(group = paste("Age", g), ederer2(d[age_grp == g])))))
write_out(rs, SUB, "relative_survival_ederer2")

# ---- 6. Pohar-Perme net survival (sensitivity) ----------------------------------------------------------
pohar_perme <- function(dd) {
  ii <- dd$row; Tt <- dd$t_years; ev <- dd$status
  kmax <- min(length(grid) - 1, ceiling(min(10, max(Tt)) / STEP))
  cumh <- 0; cumv <- 0; out <- vector("list", kmax)
  for (k in seq_len(kmax)) {
    t0 <- grid[k]; t1 <- grid[k + 1]
    risk <- Tt > t0
    if (sum(risk) < 10) break
    w <- exp(L[k, ii[risk]])
    dN <- ev[risk] == 1 & Tt[risk] <= t1
    dLs <- L[k + 1, ii[risk]] - L[k, ii[risk]]
    sw <- sum(w)
    dh <- sum(w[dN]) / sw - sum(w * dLs) / sw
    cumh <- cumh + dh; cumv <- cumv + sum(w[dN]^2) / sw^2
    out[[k]] <- data.table(years = t1, net_survival = exp(-cumh),
                           lower = exp(-(cumh + 1.96 * sqrt(cumv))), upper = exp(-(cumh - 1.96 * sqrt(cumv))))
  }
  rbindlist(out)
}
pp <- rbind(cbind(group = "All", pohar_perme(d)),
            rbindlist(lapply(levels(d$urgency2), function(g) cbind(group = g, pohar_perme(d[urgency2 == g])))))
write_out(pp[abs(years * 4 - round(years * 4)) < 1e-6 | years %in% range(years)], SUB, "net_survival_pohar_perme_quarterly")

# ---- 7. Excess hazard regression (Dickman et al. 2004) -----------------------------------------------
fu_cuts <- c(0, 30 / 365.25, 1, 3, 5, 10)
fu_lab <- c("0-30 d", "30 d-1 y", "1-3 y", "3-5 y", "5-10 y")
split <- rbindlist(lapply(seq_len(length(fu_cuts) - 1), function(j) {
  a <- fu_cuts[j]; b <- fu_cuts[j + 1]
  dd <- d[t_years > a]
  if (!nrow(dd)) return(NULL)
  dd[, .(row, urgency2, age_grp, female, era, fu = fu_lab[j],
         y = pmin(t_years, b) - a, dth = as.integer(status == 1 & t_years <= b),
         dstar = lam_at(row, pmin(t_years, b)) - lam_at(row, rep(a, .N)))]
}))
agg <- split[!is.na(urgency2), .(y = sum(y), dth = sum(dth), dstar = sum(dstar)),
             by = .(fu, urgency2, age_grp, female, era)][y > 0]
agg[, fu := factor(fu, fu_lab)]
eh <- tryCatch({
  m <- glm(dth ~ fu + urgency2 + age_grp + female + era + offset(log(y)), data = agg,
           family = poisson(link = excess_link(agg$dstar)),
           start = c(log(sum(agg$dth) / sum(agg$y)), rep(0, 4 + 1 + 2 + 1 + 2)),
           control = glm.control(maxit = 100))
  ct <- summary(m)$coefficients
  data.table(term = rownames(ct), excess_hazard_ratio = exp(ct[, 1]),
             lower = exp(ct[, 1] - 1.96 * ct[, 2]), upper = exp(ct[, 1] + 1.96 * ct[, 2]),
             p = ct[, 4], converged = m$converged, deviance = m$deviance, df_resid = m$df.residual)
}, error = function(e) data.table(note = paste("excess hazard model failed:", conditionMessage(e))))
write_out(eh, SUB, "excess_hazard_regression")

saveRDS(list(rs = rs, pp = pp, smr = smr_overall, smr_int = smr_by_int, eh = eh, verify = verify,
             lt_years = lt_years), derived_path("population.rds"))
log_msg("Done.")

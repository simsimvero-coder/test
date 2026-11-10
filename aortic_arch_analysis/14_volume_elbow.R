# =============================================================================
# 14_volume_elbow.R
# Volume-outcome threshold by the "elbow of the curve" (method as in the HVS AV
# VSARR registry study, Eur J Cardiothorac Surg): centre-level analysis.
# Early outcomes: risk-adjusted in-hospital and 30-day mortality per centre against
# the centre's mean annual case volume (operations in the group per year of registry
# participation), restricted cubic spline with 3 knots, weights 1/variance.
# Long-term: adjusted HR per centre versus all other centres (Cox), variance = SE^2,
# same spline. Centres with >= 30 cases; sensitivity: small centres merged (1-10,
# 11-20, 21-29 cases). Adjustment: EuroSCORE II only, and age, sex, EuroSCORE II,
# redo and year. Elbow: volume beyond which the slope is no longer significantly
# negative, and the point of maximum curvature. Per-centre values are used only as
# data points; no centre is named or ranked.
# Output: results/<method>/13_volume_outcome/elbow/*
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
start_log("14_volume_elbow")
suppressMessages({ library(data.table); library(survival); library(ggplot2) })
OUTDIR <- file.path(res_path("13_volume_outcome"), "elbow")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
d0    <- copy(analysis_datasets()[[1]])
specs <- load_derived("specs.rds")
d0[, centrum := as.character(centrum)]
ADJ <- list("EuroSCORE II" = "logit_es2",
            "age, sex, EuroSCORE II, redo, year" = paste(c("year_c", covar_terms(COVARS_M1, specs)), collapse = " + "))
MIN_CASES <- 30
PERIODS <- list("2013-2024" = c(2013, 2024), "2018-2024" = c(2018, 2024))

rcs3 <- function(x, kn) {
  nl <- (pmax(x - kn[1], 0)^3 - pmax(x - kn[2], 0)^3 * (kn[3] - kn[1]) / (kn[3] - kn[2]) +
         pmax(x - kn[3], 0)^3 * (kn[2] - kn[1]) / (kn[3] - kn[2])) / (kn[3] - kn[1])^2
  cbind(v1 = x, v2 = nl)
}

fit_curve <- function(cd, yvar, wvar) {
  kn <- as.numeric(quantile(cd$vol, c(.10, .50, .90)))
  B  <- rcs3(cd$vol, kn)
  m  <- lm(cd[[yvar]] ~ B, weights = cd[[wvar]])
  b  <- coef(m); V <- vcov(m)
  grid <- seq(min(cd$vol), max(cd$vol), length.out = 200)
  Xg <- cbind(1, rcs3(grid, kn))
  fit <- drop(Xg %*% b); se <- sqrt(rowSums((Xg %*% V) * Xg))
  h  <- 1e-3 * diff(range(cd$vol))
  D  <- cbind(0, (rcs3(grid + h, kn) - rcs3(grid - h, kn)) / (2 * h))
  slope <- drop(D %*% b); slope_se <- sqrt(rowSums((D %*% V) * D))
  df <- m$df.residual
  sig <- slope + qt(.975, df) * slope_se < 0
  elbow_sig <- if (any(sig)) { i <- max(which(sig)); if (i < length(grid)) grid[i + 1] else NA_real_ } else NA_real_
  xs <- (grid - min(grid)) / diff(range(grid)); ys <- (fit - min(fit)) / diff(range(fit))
  chord <- ys[1] + (ys[length(ys)] - ys[1]) * xs
  elbow_knee <- grid[which.max(abs(chord - ys))]
  p_nonlin <- summary(m)$coefficients["Bv2", 4]
  p_overall <- anova(lm(cd[[yvar]] ~ 1, weights = cd[[wvar]]), m)[2, "Pr(>F)"]
  list(curve = data.table(vol = grid, fit = fit, lower = fit - qt(.975, df) * se, upper = fit + qt(.975, df) * se,
                          slope = slope, slope_lower = slope - qt(.975, df) * slope_se, slope_upper = slope + qt(.975, df) * slope_se),
       summary = data.table(n_units = nrow(cd), knots = paste(round(kn, 1), collapse = "/"),
                            p_overall = p_overall, p_nonlinear = p_nonlin,
                            elbow_last_significant_slope = elbow_sig, elbow_max_curvature = elbow_knee))
}

merge_small <- function(cd, outcome_cols) {
  big <- cd[n_total >= MIN_CASES]
  small <- cd[n_total < MIN_CASES]
  if (!nrow(small)) return(big)
  small[, band := cut(n_total, c(0, 10, 20, 30), c("merged 1-10", "merged 11-20", "merged 21-29"), include.lowest = TRUE)]
  rbind(big, small[, c(list(centrum = as.character(band[1]), vol = mean(vol), n_total = sum(n_total)),
                       lapply(.SD, sum)), by = band, .SDcols = outcome_cols][, band := NULL], fill = TRUE)
}

early <- list(); early_curves <- list(); late <- list(); late_curves <- list(); hr_rows <- list()
for (an in names(ADJ)) for (gr in c("Acute", "Non-acute")) for (pn in names(PERIODS)) {
  rhs <- ADJ[[an]]; tag <- if (an == "EuroSCORE II") "es2" else "full"
  per <- PERIODS[[pn]]
  dg  <- d0[urgency2 == gr & year >= per[1] & year <= per[2]]
  yrs <- merge(dg[, .(n = .N), by = centrum], d0[year >= per[1] & year <= per[2], .(years = uniqueN(year)), by = centrum], by = "centrum")
  yrs[, vol := n / years]
  for (v in c("dzh", "d30")) {
    d <- dg[!is.na(get(v))]
    d <- d[complete.cases(d[, .(age, female, logit_es2, prior_cs, year_c)])]
    pexp <- predict(glm(as.formula(paste(v, "~", rhs)), binomial, d), type = "response")
    d[, e := pexp]
    overall <- mean(d[[v]])
    cd <- d[, .(n_total = .N, obs = sum(get(v)), exp = sum(e), var_o = sum(e * (1 - e))), by = centrum]
    cd <- merge(cd, yrs[, .(centrum, vol)], by = "centrum")
    for (sens in c("centres >= 30 cases", "small centres merged")) {
      cc <- if (sens == "centres >= 30 cases") cd[n_total >= MIN_CASES] else merge_small(cd, c("obs", "exp", "var_o"))
      cc[, `:=`(crude_rate = obs / n_total, adj_rate = overall * obs / exp)]
      cc[, `:=`(w_adj = 1 / ((overall / exp)^2 * var_o),
                w_crude = n_total / pmax(crude_rate * (1 - crude_rate), 1 / n_total))]
      for (sc in if (an == names(ADJ)[1]) c("adjusted", "crude") else "adjusted") {
        f <- fit_curve(cc, if (sc == "adjusted") "adj_rate" else "crude_rate", if (sc == "adjusted") "w_adj" else "w_crude")
        lab <- data.table(adjustment = if (sc == "crude") "none" else an, group = gr, period = pn, outcome = v, scale = sc, centres = sens)
        early[[length(early) + 1]] <- cbind(lab, f$summary)
        early_curves[[length(early_curves) + 1]] <- cbind(lab, f$curve)
        if (sens == "centres >= 30 cases" && sc == "adjusted") {
          pts <- cbind(lab, cc[, .(centrum, vol, n_total, rate = adj_rate, w = w_adj)])
          p <- ggplot() + geom_ribbon(data = f$curve, aes(vol, ymin = 100 * lower, ymax = 100 * upper), alpha = .15) +
            geom_line(data = f$curve, aes(vol, 100 * fit)) +
            geom_point(data = pts, aes(vol, 100 * rate, size = w), alpha = .6) +
            geom_vline(xintercept = f$summary$elbow_last_significant_slope, linetype = 2) +
            scale_size_continuous(name = "Weight (1/variance)") +
            labs(x = "Annual case volume of the centre", y = "Risk-adjusted mortality (%)",
                 title = paste(gr, pn, if (v == "dzh") "in-hospital mortality" else "30-day mortality"), subtitle = paste("Adjusted for", an)) + theme_bw()
          ggsave(file.path(OUTDIR, paste0("early_", v, "_", tag, "_", gsub("-", "", gr), "_", gsub("-", "_", pn), ".png")), p, width = 7, height = 4.5, dpi = 300)
        }
      }
    }
  }
  ds <- dg[!is.na(time_d) & !is.na(status)]
  ds <- ds[complete.cases(ds[, .(age, female, logit_es2, prior_cs, year_c)])]
  hr <- rbindlist(lapply(unique(ds$centrum), function(cn) {
    ds[, this := as.integer(centrum == cn)]
    s <- summary(coxph(as.formula(paste("Surv(time_d, status) ~ this +", rhs)), ds))$coefficients["this", ]
    data.table(centrum = cn, n_total = sum(ds$this), deaths = sum(ds[this == 1]$status == 1),
               log_hr = s[["coef"]], se = s[["se(coef)"]])
  }))
  hr <- merge(hr, yrs[, .(centrum, vol)], by = "centrum")
  hr[, `:=`(hr = exp(log_hr), lower = exp(log_hr - 1.96 * se), upper = exp(log_hr + 1.96 * se), variance = se^2, w = 1 / se^2)]
  hr_rows[[length(hr_rows) + 1]] <- cbind(adjustment = an, group = gr, period = pn, hr)
  for (sens in c("centres >= 30 cases", "small centres merged")) {
    hh <- if (sens == "centres >= 30 cases") hr[n_total >= MIN_CASES] else {
      small <- hr[n_total < MIN_CASES]
      if (!nrow(small)) hr[n_total >= MIN_CASES] else {
        small[, band := cut(n_total, c(0, 10, 20, 30), c("merged 1-10", "merged 11-20", "merged 21-29"), include.lowest = TRUE)]
        mg <- rbindlist(lapply(split(small, small$band, drop = TRUE), function(s) {
          ids <- s$centrum
          ds[, this := as.integer(centrum %in% ids)]
          cf <- summary(coxph(as.formula(paste("Surv(time_d, status) ~ this +", rhs)), ds))$coefficients["this", ]
          data.table(centrum = as.character(s$band[1]), n_total = sum(s$n_total), vol = mean(s$vol),
                     log_hr = cf[["coef"]], se = cf[["se(coef)"]])
        }))
        mg[, w := 1 / se^2]
        rbind(hr[n_total >= MIN_CASES, .(centrum, n_total, vol, log_hr, se, w)], mg)
      }
    }
    f <- fit_curve(hh, "log_hr", "w")
    lab <- data.table(adjustment = an, group = gr, period = pn, outcome = "overall survival (HR vs all other centres)", centres = sens)
    s <- copy(f$summary)
    late[[length(late) + 1]] <- cbind(lab, s)
    late_curves[[length(late_curves) + 1]] <- cbind(lab, f$curve[, .(vol, hr = exp(fit), lower = exp(lower), upper = exp(upper), slope, slope_lower, slope_upper)])
    if (sens == "centres >= 30 cases") {
      p <- ggplot() + geom_ribbon(data = f$curve, aes(vol, ymin = exp(lower), ymax = exp(upper)), alpha = .15) +
        geom_line(data = f$curve, aes(vol, exp(fit))) +
        geom_point(data = hh, aes(vol, exp(log_hr), size = w), alpha = .6) +
        geom_hline(yintercept = 1, colour = "grey50") +
        geom_vline(xintercept = f$summary$elbow_last_significant_slope, linetype = 2) +
        scale_y_log10() + scale_size_continuous(name = "Weight (1/variance)") +
        labs(x = "Annual case volume of the centre", y = "Adjusted HR for death (vs all other centres)",
             title = paste(gr, pn, "long-term survival"), subtitle = paste("Adjusted for", an)) + theme_bw()
      ggsave(file.path(OUTDIR, paste0("late_survival_", tag, "_", gsub("-", "", gr), "_", gsub("-", "_", pn), ".png")), p, width = 7, height = 4.5, dpi = 300)
    }
  }
}
fwrite(rbindlist(early), file.path(OUTDIR, "elbow_early_outcomes.csv"))
fwrite(rbindlist(early_curves), file.path(OUTDIR, "elbow_early_curves.csv"))
fwrite(rbindlist(late), file.path(OUTDIR, "elbow_long_term_survival.csv"))
fwrite(rbindlist(late_curves), file.path(OUTDIR, "elbow_long_term_curves.csv"))
fwrite(rbindlist(hr_rows), file.path(OUTDIR, "hr_per_centre.csv"))
options(width = 250)
print(rbindlist(early)[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 3) else z)], row.names = FALSE)
print(rbindlist(late)[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 3) else z)], row.names = FALSE)
log_msg("Done. Output in ", OUTDIR)

# =============================================================================
# 03_volume_trends.R
# Annual national procedure counts and population-adjusted rates; Poisson vs
# negative binomial regression; non-linearity; acute vs non-acute; centre volume.
# Output: results/03_volume/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("03_volume_trends")
SUB <- "03_volume"

# National volume uses the SURGICAL cohort (all surgically selected registrations,
# including those not classifiable as acute/non-acute). The acute and non-acute series
# are its classified subsets (= the analysis cohort).
surg <- load_derived("surgical_cohort.rds")

# ---- population (CBS 03759ned, average population per year) --------------------------
read_cbs_population <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  rows <- grep('"Nederland";', lines, value = TRUE, fixed = TRUE, useBytes = TRUE)
  if (!length(rows)) return(NULL)
  parts <- lapply(strsplit(rows, ";", fixed = TRUE), function(x) gsub('"', "", x))
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", x)))
  dt <- data.table(year = as.integer(vapply(parts, `[`, "", 2)),
                   pop_jan1 = num(vapply(parts, `[`, "", 3)),
                   pop_avg  = num(vapply(parts, `[`, "", 8)))
  dt[]
}
pop <- read_cbs_population(CONFIG$population_path)
if (is.null(pop)) {
  log_msg("Population file not found: rates are not computed (counts only).")
} else {
  write_out(pop[year %in% CONFIG$study_years], SUB, "population_cbs_03759ned")
}

# ---- annual counts ---------------------------------------------------------------------
yrs <- CONFIG$study_years
cnt <- CJ(year = yrs, type = c("All", "Acute", "Non-acute"))
cnt[, n := mapply(function(y, t) {
  d <- surg[year == y]
  if (t == "All") nrow(d) else if (t == "Acute") sum(d$acute == 1, na.rm = TRUE)
  else sum(d$acute == 0, na.rm = TRUE)
}, year, type)]
n_unknown_urg <- surg[is.na(acute), .N, by = year][order(year)]   # not classifiable as acute/non-acute
write_out(n_unknown_urg, SUB, "counts_not_classifiable")
stopifnot(all(cnt[type == "All", n] >= cnt[type == "Acute", n] + cnt[type == "Non-acute", n]))
if (!is.null(pop)) {
  cnt <- merge(cnt, pop[, .(year, pop_avg)], by = "year", all.x = TRUE)
  if (anyNA(cnt$pop_avg)) log_msg("WARNING: population missing for some years: ",
                                  paste(unique(cnt[is.na(pop_avg), year]), collapse = ","))
  cnt[, rate_per_million := 1e6 * n / pop_avg]
}
cnt[, year_c := year - CONFIG$reference_year]
ac <- dcast(cnt, year ~ type, value.var = "n")
ac[, `Not classifiable` := All - Acute - `Non-acute`]
write_out(ac, SUB, "annual_counts")
if (!is.null(pop)) write_out(dcast(cnt, year ~ type, value.var = "rate_per_million"), SUB, "annual_rates_per_million")

# ---- count models ------------------------------------------------------------------------
fit_counts <- function(d, label) {
  use_off <- !is.null(pop) && all(!is.na(d$pop_avg))
  off <- if (use_off) " + offset(log(pop_avg / 1e6))" else ""
  f_lin <- as.formula(paste("n ~ year_c", off))
  f_ns  <- as.formula(paste("n ~ ns(year_c, df = 2)", off))
  pois <- glm(f_lin, family = poisson, data = d)
  disp <- sum(residuals(pois, "pearson")^2) / df.residual(pois)
  nb <- tryCatch(suppressWarnings(MASS::glm.nb(f_lin, data = d)), error = function(e) NULL)
  lrt_nb <- if (!is.null(nb)) max(0, 2 * (as.numeric(logLik(nb)) - as.numeric(logLik(pois)))) else NA
  p_nb <- if (!is.na(lrt_nb)) 0.5 * pchisq(lrt_nb, 1, lower.tail = FALSE) else NA
  use_nb <- !is.na(p_nb) && p_nb < 0.05
  fam_fit <- function(f) if (use_nb) suppressWarnings(MASS::glm.nb(f, data = d)) else glm(f, family = poisson, data = d)
  m_lin <- fam_fit(f_lin); m_ns <- fam_fit(f_ns)
  lr_nl <- 2 * (as.numeric(logLik(m_ns)) - as.numeric(logLik(m_lin)))
  p_nl <- pchisq(max(lr_nl, 0), 1, lower.tail = FALSE)
  b <- coef(m_lin)["year_c"]; se <- sqrt(vcov(m_lin)["year_c", "year_c"])
  res <- data.table(series = label, n_years = nrow(d), total = sum(d$n),
                    outcome_scale = if (use_off) "rate (offset: average population)" else "count",
                    pearson_dispersion = disp, lrt_nb_vs_poisson = lrt_nb, p_nb_vs_poisson = p_nb,
                    model = if (use_nb) "negative binomial" else "Poisson",
                    irr_per_year = exp(b), irr_lower = exp(b - 1.96 * se), irr_upper = exp(b + 1.96 * se),
                    annual_pct_change = 100 * (exp(b) - 1),
                    # quasi-Poisson-type interval: SE inflated by sqrt(max(1, dispersion)),
                    # because the NB test has little power with 12 annual points
                    irr_lower_dispersion_scaled = exp(b - 1.96 * se * sqrt(max(1, disp))),
                    irr_upper_dispersion_scaled = exp(b + 1.96 * se * sqrt(max(1, disp))),
                    p_trend = 2 * pnorm(-abs(b / se)),
                    aic_linear = AIC(m_lin), aic_spline = AIC(m_ns), p_nonlinear = p_nl)
  # fitted values (linear model) with CI on the response (rate) scale
  pr <- predict(m_lin, newdata = d, type = "link", se.fit = TRUE)
  scale <- if (use_off) 1e6 / d$pop_avg else 1
  fitted <- data.table(series = label, year = d$year,
                       fit = exp(pr$fit) * scale,
                       lower = exp(pr$fit - 1.96 * pr$se.fit) * scale,
                       upper = exp(pr$fit + 1.96 * pr$se.fit) * scale)
  list(res = res, fitted = fitted, model = m_lin)
}

series_fits <- lapply(c("All", "Acute", "Non-acute"), function(t) fit_counts(cnt[type == t], t))
vol_res <- rbindlist(lapply(series_fits, `[[`, "res"))
vol_fit <- rbindlist(lapply(series_fits, `[[`, "fitted"))

# Sensitivity: excluding 2024 (registration lag) and excluding 2020-2021 (COVID-19)
sens <- rbindlist(list(
  rbindlist(lapply(c("All", "Acute", "Non-acute"), function(t)
    fit_counts(cnt[type == t & year != 2024], t)$res))[, analysis := "excluding 2024"],
  rbindlist(lapply(c("All", "Acute", "Non-acute"), function(t)
    fit_counts(cnt[type == t & !year %in% 2020:2021], t)$res))[, analysis := "excluding 2020-2021"],
  rbindlist(lapply(c("All", "Acute", "Non-acute"), function(t)
    fit_counts(cnt[type == t & year >= CONFIG$trend_sens_start], t)$res))[, analysis := paste("from", CONFIG$trend_sens_start)]))
vol_res[, analysis := "primary"]
write_out(rbind(vol_res, sens), SUB, "volume_trend_models")
write_out(vol_fit, SUB, "volume_fitted_values")

# Does the trend differ between acute and non-acute? (year x group interaction)
d2 <- cnt[type != "All"]
off <- if (!is.null(pop)) " + offset(log(pop_avg / 1e6))" else ""
use_nb <- any(vol_res[series != "All", model == "negative binomial"])
mfun <- function(f) if (use_nb) suppressWarnings(MASS::glm.nb(f, data = d2)) else glm(f, family = poisson, data = d2)
m0 <- mfun(as.formula(paste("n ~ year_c + type", off)))
m1 <- mfun(as.formula(paste("n ~ year_c * type", off)))
lr <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
b_int <- coef(m1)["year_c:typeNon-acute"]; se_int <- sqrt(vcov(m1)["year_c:typeNon-acute", "year_c:typeNon-acute"])
write_out(data.table(test = "year x acute/non-acute interaction", lr_statistic = lr,
                     p = pchisq(lr, 1, lower.tail = FALSE),
                     ratio_of_irr_nonacute_vs_acute = exp(b_int),
                     lower = exp(b_int - 1.96 * se_int), upper = exp(b_int + 1.96 * se_int),
                     model = if (use_nb) "negative binomial" else "Poisson"),
          SUB, "volume_interaction_urgency")

# Diagnostics: Pearson residuals by year
diag <- rbindlist(lapply(series_fits, function(s)
  data.table(series = s$res$series, year = cnt[type == s$res$series, year],
             pearson_residual = residuals(s$model, "pearson"))))
write_out(diag, SUB, "volume_model_residuals")

# ---- centre volume -------------------------------------------------------------------------
cv <- surg[, .N, by = .(centre_lab, year)]
cv <- merge(CJ(centre_lab = levels(surg$centre_lab), year = yrs), cv, all.x = TRUE)
cv[is.na(N), N := 0L]
write_out(dcast(cv, centre_lab ~ year, value.var = "N"), SUB, "centre_annual_volume")
cv_sum <- cv[, .(total = sum(N), median_annual = median(N), min_annual = min(N), max_annual = max(N),
                 years_active = sum(N > 0)), by = centre_lab]
write_out(cv_sum, SUB, "centre_volume_summary")
saveRDS(list(cnt = cnt, fitted = vol_fit, res = vol_res, cv = cv, has_pop = !is.null(pop)),
        derived_path("volume.rds"))
log_msg("Done.")

# =============================================================================
# period_cr2.R  (standalone; run after run_all.R)
# Figure 2C: 2019-2024 versus 2013-2018, with cluster-robust inference.
# For each group (acute, non-acute) and each of four perfusion measures, one
# regression with period as the only explanatory variable and CR2 standard
# errors clustered by centre (Satterthwaite degrees of freedom):
#   circulatory arrest, deep hypothermia (<=20 C) -> logistic regression (OR)
#   nadir temperature, arrest duration            -> linear regression (difference)
# Deep hypothermia and nadir temperature: patients with recorded temperature.
# Arrest duration: patients with circulatory arrest.
# Writes only aggregate numbers. Output: results/period_cr2/
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages(library(data.table))
SUB <- "period_cr2"
options(width = 200)
d <- copy(load_derived("cohort.rds"))
d <- d[!is.na(urgency2)]
d[, `:=`(centrum = as.character(centrum), late = as.integer(year >= 2019),
         deep = fifelse(is.na(temp_min), NA_integer_, as.integer(temp_min <= 20)),
         ca_arrest = fifelse(circ_arrest_b %in% 1, ca_min, NA_real_))]
spec <- list(arrest = list(y = "circ_arrest_b", fam = "logistic"), deep = list(y = "deep", fam = "logistic"),
             temp = list(y = "temp_min", fam = "linear"), ca = list(y = "ca_arrest", fam = "linear"))
res <- rbindlist(lapply(c("Acute", "Non-acute"), function(g) rbindlist(lapply(names(spec), function(m) {
  s <- spec[[m]]; dd <- d[urgency2 == g & !is.na(get(s$y))]
  f <- as.formula(paste(s$y, "~ late"))
  fit <- if (s$fam == "logistic") glm(f, binomial, dd) else lm(f, dd)
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = dd$centrum, test = "Satterthwaite")
  r <- ct[ct$Coef == "late", ]; tq <- qt(.975, r$df_Satt)
  est <- r$beta; lo <- est - tq * r$SE; hi <- est + tq * r$SE
  if (s$fam == "logistic") { est <- exp(est); lo <- exp(lo); hi <- exp(hi) }
  summ <- if (s$fam == "logistic") dd[, .(v = round(100 * mean(get(s$y)), 1)), by = late] else dd[, .(v = round(median(get(s$y)), 1)), by = late]
  data.table(group = g, measure = m, model = s$fam, n = nrow(dd), centres = uniqueN(dd$centrum),
             value_2013_2018 = summ[late == 0, v], value_2019_2024 = summ[late == 1, v],
             effect = if (s$fam == "logistic") "OR" else "difference", estimate = est, lower = lo, upper = hi,
             df = r$df_Satt, p = r$p_Satt)
}))))
write_out(res, SUB, "figure2C_period_cr2")
out <- capture.output(print(res[, .(group, measure, model, n, value_2013_2018, value_2019_2024, effect,
                                     est = sprintf("%.3f (%.3f-%.3f)", estimate, lower, upper), df = round(df, 1), p = signif(p, 3))]))
writeLines(out, file.path(out_dir(SUB), "overzicht_period_cr2.txt")); cat(out, sep = "\n")

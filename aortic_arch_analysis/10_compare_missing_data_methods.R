# =============================================================================
# 10_compare_missing_data_methods.R
# Places the key adjusted estimates from the multiple-imputation analysis and the
# complete-case analysis side by side (outcome trends, acute vs non-acute,
# organ-protection trends, Cox models, sensitivity analyses).
# Requires scripts 04-06 (and optionally 08) to have run for both methods.
# Output: results/comparison_missing_data/*
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("10_compare_missing_data_methods")
SUB <- "comparison_missing_data"

need <- c("outcomes.rds", "organ.rds", "survival.rds")
files <- outer(need, MISSING_METHODS, function(n, m) mapply(derived_method, n, m))
if (!all(file.exists(files))) {
  stop("Results for both methods are needed. Missing: ", paste(basename(files[!file.exists(files)]), collapse = ", "),
       ". Run run_all.R (it runs scripts 04-09 for both methods).")
}
ld <- function(name, m) readRDS(derived_method(name, m))
LAB <- c(imputation = "MI", complete_case = "CC")

# side-by-side helper: x_mi and x_cc share key columns
side <- function(mi, cc, keys) {
  out <- merge(mi, cc, by = keys, all = TRUE, suffixes = c("_MI", "_CC"))
  if (all(c("estimate_MI", "estimate_CC") %in% names(out)))
    out[, ratio_CC_vs_MI := estimate_CC / estimate_MI]
  out[]
}

# ---- 1. Early outcomes: OR per year -------------------------------------------------------
yr <- lapply(MISSING_METHODS, function(m) ld("outcomes.rds", m)$year[
  , .(outcome, label, role, n, events, estimate = or_per_year, lower, upper, p, p_holm)])
t_year <- side(yr[[1]], yr[[2]][, !"label"][, !"role"], "outcome")
setcolorder(t_year, c("outcome", "label", "role"))
write_out(t_year, SUB, "outcome_year_trend_MI_vs_CC")

# ---- 1b. Outcome trends from CONFIG$trend_sens_start ------------------------------------------
yf <- lapply(MISSING_METHODS, function(m) {
  y <- ld("outcomes.rds", m)$year_from
  if (is.null(y) || !nrow(y)) return(NULL)
  y[term == "year_c", .(outcome, label, n, events, estimate = ratio, lower = ratio_lower, upper = ratio_upper, p)]
})
if (!any(vapply(yf, is.null, logical(1))))
  write_out(side(yf[[1]], yf[[2]][, !"label"], "outcome"), SUB, paste0("outcome_year_trend_from", CONFIG$trend_sens_start, "_MI_vs_CC"))

# ---- 1c. Trends by group and start year ---------------------------------------------------------
gt <- lapply(MISSING_METHODS, function(m) {
  y <- ld("outcomes.rds", m)$group_trends
  y[, .(outcome, label, group, start_year, end_year, n, events, estimate = or_per_year, lower, upper, p)]
})
write_out(side(gt[[1]], gt[[2]][, !"label"], c("outcome", "group", "start_year", "end_year")), SUB,
          "trends_by_group_and_start_MI_vs_CC")

# ---- 2. Acute vs non-acute ----------------------------------------------------------------------
ur <- lapply(MISSING_METHODS, function(m) ld("outcomes.rds", m)$urgency[
  , .(outcome, adjustment, estimate = ratio, lower = ratio_lower, upper = ratio_upper, p)])
write_out(side(ur[[1]], ur[[2]], c("outcome", "adjustment")), SUB, "acute_vs_nonacute_MI_vs_CC")

# ---- 3. Failure to rescue ------------------------------------------------------------------------
ft <- lapply(MISSING_METHODS, function(m) ld("outcomes.rds", m)$ftr[term == "year_c",
  .(definition, numerator, period, n, events, estimate = ratio, lower = ratio_lower, upper = ratio_upper, p)])
write_out(side(ft[[1]], ft[[2]], c("definition", "numerator", "period")), SUB, "ftr_year_trend_MI_vs_CC")

# ---- 4. Organ protection ---------------------------------------------------------------------------
pr <- lapply(MISSING_METHODS, function(m) {
  o <- ld("organ.rds", m)
  rbind(o$practice[, .(outcome, population, adjusted, model, period, n, estimate, lower, upper, p)],
        o$ordinal[, .(outcome, population = "CPB", adjusted, model = "ordinal / binary", period, n, estimate, lower, upper, p)])
})
write_out(side(pr[[1]], pr[[2]], c("outcome", "population", "adjusted", "model", "period")), SUB, "organ_protection_MI_vs_CC")

# ---- 5. Survival (Cox) -------------------------------------------------------------------------------
cx <- lapply(MISSING_METHODS, function(m) {
  s <- ld("survival.rds", m)
  rbind(s$cox[, .(model, term, n, events, estimate = ratio, lower = ratio_lower, upper = ratio_upper, p)],
        s$cox_lm[, .(model, term, n, events, estimate = ratio, lower = ratio_lower, upper = ratio_upper, p)])
})
cx <- lapply(cx, function(d) d[!grepl("^ns\\(", term)])
write_out(side(cx[[1]], cx[[2]], c("model", "term")), SUB, "cox_MI_vs_CC")

# ---- 6. Sensitivity analyses (if run) ---------------------------------------------------------------
sf <- vapply(MISSING_METHODS, function(m) file.exists(derived_method("sensitivity.rds", m)), logical(1))
if (all(sf)) {
  se <- lapply(MISSING_METHODS, function(m) {
    d <- ld("sensitivity.rds", m)[, .(analysis, outcome, term, n, events, estimate = or, lower, upper, p)]
    d[grepl("^M1 primary", analysis), analysis := "M1 primary"]
    d[grepl("^Centre fixed effects", analysis), analysis := "Centre fixed effects, CR2"]
    d
  })
  write_out(side(se[[1]], se[[2]], c("analysis", "outcome", "term")), SUB, "sensitivity_MI_vs_CC")
}

# ---- 7. Figure: OR per year by method ----------------------------------------------------------------
fd <- rbindlist(lapply(seq_along(MISSING_METHODS), function(i)
  yr[[i]][, method := c("Multiple imputation", "Complete case")[i]]))
fd[, label := factor(label, rev(unique(yr[[1]]$label)))]
p <- ggplot(fd, aes(estimate, label, colour = method)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = .5), size = .3) +
  scale_colour_manual(values = PAL[1:2], name = NULL) +
  labs(x = "Adjusted odds ratio per calendar year (95% CI)", y = NULL,
       title = "Temporal trends: multiple imputation vs complete-case analysis") + theme_pub()
save_fig(p, SUB, "year_trend_MI_vs_CC", 8, 5.5)

# ---- 8. Workbook and short summary ---------------------------------------------------------------------
sheets <- list(outcome_trends = t_year,
               trends_by_group = side(gt[[1]], gt[[2]][, !"label"], c("outcome", "group", "start_year", "end_year")),
               acute_vs_nonacute = side(ur[[1]], ur[[2]], c("outcome", "adjustment")),
               ftr = side(ft[[1]], ft[[2]], c("definition", "numerator", "period")),
               organ_protection = side(pr[[1]], pr[[2]], c("outcome", "population", "adjusted", "model", "period")),
               cox = side(cx[[1]], cx[[2]], c("model", "term")))
if (all(sf)) sheets$sensitivity <- side(se[[1]], se[[2]], c("analysis", "outcome", "term"))
writexl::write_xlsx(sheets, file.path(out_dir(SUB), "MI_vs_CC_comparison.xlsx"))

g <- function(x, d = 3) formatC(x, format = "f", digits = d)
tt <- t_year[!is.na(estimate_MI) & !is.na(estimate_CC)]
md <- c("# Multiple imputation vs complete-case analysis (auto-generated)", "",
        sprintf("Generated %s. Odds ratio per calendar year, primary model M1.", format(Sys.time())), "",
        "| Outcome | MI: n / events | MI: OR (95% CI) | CC: n / events | CC: OR (95% CI) |",
        "|---|---|---|---|---|",
        sprintf("| %s | %d / %d | %s (%s-%s) | %d / %d | %s (%s-%s) |", tt$label,
                tt$n_MI, tt$events_MI, g(tt$estimate_MI), g(tt$lower_MI), g(tt$upper_MI),
                tt$n_CC, tt$events_CC, g(tt$estimate_CC), g(tt$lower_CC), g(tt$upper_CC)),
        "", "Under complete-case analysis, patients with any missing model variable are excluded per model, so n and",
        "events are smaller. Material differences between the two columns suggest that missingness is related to",
        "the outcome or to calendar time, and should be discussed.")
writeLines(md, file.path(out_dir(SUB), "MI_vs_CC_summary.md"))
log_msg("Done.")

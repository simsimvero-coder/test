# =============================================================================
# run_all.R — runs the complete pipeline, each script in a fresh R process
# (clean session per step), stopping at the first failure.
#
# Usage (from this folder):   Rscript run_all.R
# Optional environment variables: ARCH_RAW_DATA, ARCH_RESULTS, ARCH_DERIVED,
# ARCH_M (imputations), ARCH_NSIM, ARCH_CORES, ARCH_GLMM_FALLBACK (see R/00_config.R).
# Scripts 04-09 run twice (multiple imputation and complete case); script 10 compares them.
# =============================================================================

if (!file.exists("R/00_config.R")) stop("Run this script from the project folder (aortic_arch_analysis/).")

required <- c("data.table", "ggplot2", "lme4", "survival", "splines", "mice", "ordinal",
              "clubSandwich", "MASS", "patchwork", "scales", "writexl", "readxl", "haven", "parallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the missing R packages first: ", paste(missing, collapse = ", "))

shared   <- c("01_data_validation.R", "02_cohort_descriptives.R", "03_volume_trends.R",
              "07_population_survival.R")
adjusted <- c("04_clinical_outcomes.R", "05_organ_protection.R", "06_survival_analysis.R",
              "08_sensitivity_analyses.R", "09_tables_figures.R")
# Adjusted analyses run twice: multiple imputation and complete-case analysis
steps <- c(lapply(shared, function(s) list(script = s, method = "imputation")),
           unlist(lapply(c("imputation", "complete_case"), function(m)
             lapply(adjusted, function(s) list(script = s, method = m))), recursive = FALSE),
           list(list(script = "10_compare_missing_data_methods.R", method = "imputation")))

rscript <- file.path(R.home("bin"), "Rscript")
timing <- data.frame(script = vapply(steps, `[[`, "", "script"),
                     method = vapply(steps, `[[`, "", "method"), seconds = NA_real_, status = NA_character_)
timing$method[timing$script %in% c(shared, "10_compare_missing_data_methods.R")] <- "-"
for (i in seq_along(steps)) {
  st <- steps[[i]]
  Sys.setenv(ARCH_MISSING = st$method)      # inherited by the child R process
  message("\n==== ", st$script, if (st$script %in% adjusted) paste0("  [", st$method, "]") else "", " ====")
  t0 <- Sys.time()
  rc <- system2(rscript, shQuote(st$script))   # system2 quotes the command itself
  timing$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  timing$status[i] <- if (rc == 0) "ok" else paste("failed (exit", rc, ")")
  if (rc != 0) { print(timing); stop("Pipeline stopped: ", st$script, " failed. See results/logs/.") }
}
print(timing)
res <- Sys.getenv("ARCH_RESULTS", "results")
write.csv(timing, file.path(res, "logs", "pipeline_timing.csv"), row.names = FALSE)
message("Pipeline complete. Reports: ", file.path(res, "imputation", "analysis_report.md"), ", ",
        file.path(res, "complete_case", "analysis_report.md"), " and ",
        file.path(res, "comparison_missing_data", "MI_vs_CC_summary.md"))

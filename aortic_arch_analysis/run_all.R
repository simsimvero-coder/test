# =============================================================================
# run_all.R — runs the complete pipeline, each script in a fresh R process
# (clean session per step), stopping at the first failure.
#
# Usage (from this folder):   Rscript run_all.R
# Optional environment variables: ARCH_RAW_DATA, ARCH_RESULTS, ARCH_DERIVED,
# ARCH_M (imputations), ARCH_NSIM, ARCH_CORES, ARCH_GLMM_FALLBACK (see R/00_config.R)
# =============================================================================

if (!file.exists("R/00_config.R")) stop("Run this script from the project folder (aortic_arch_analysis/).")

required <- c("data.table", "ggplot2", "lme4", "survival", "splines", "mice", "ordinal",
              "clubSandwich", "MASS", "patchwork", "scales", "writexl", "readxl", "haven", "parallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the missing R packages first: ", paste(missing, collapse = ", "))

scripts <- c("01_data_validation.R", "02_cohort_descriptives.R", "03_volume_trends.R",
             "04_clinical_outcomes.R", "05_organ_protection.R", "06_survival_analysis.R",
             "07_population_survival.R", "08_sensitivity_analyses.R", "09_tables_figures.R")
rscript <- file.path(R.home("bin"), "Rscript")
timing <- data.frame(script = scripts, seconds = NA_real_, status = NA_character_)
for (i in seq_along(scripts)) {
  message("\n==== ", scripts[i], " ====")
  t0 <- Sys.time()
  rc <- system2(rscript, scripts[i])
  timing$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  timing$status[i] <- if (rc == 0) "ok" else paste("failed (exit", rc, ")")
  if (rc != 0) { print(timing); stop("Pipeline stopped: ", scripts[i], " failed. See results/logs/.") }
}
print(timing)
res <- Sys.getenv("ARCH_RESULTS", "results")
write.csv(timing, file.path(res, "logs", "pipeline_timing.csv"), row.names = FALSE)
message("Pipeline complete. Report: ", file.path(res, "analysis_report.md"))

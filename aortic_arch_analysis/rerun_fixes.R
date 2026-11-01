# =============================================================================
# rerun_fixes.R — re-runs only the steps touched by the fixes after the full run:
#   - failure-to-rescue trends from 2015 (04b, both methods)
#   - survival tables: no estimates where < 10 patients remain at risk (06, both methods)
#   - sensitivity table: complete-case row (08, imputation)
# then the tables/figures (09) and the MI vs CC comparison (10).
# Run after the full pipeline has finished:   source("rerun_fixes.R")
# =============================================================================

if (!file.exists("R/00_config.R")) stop("Run this script from the project folder (aortic_arch_analysis/).")
rscript <- file.path(R.home("bin"), "Rscript")
steps <- list(c("04b_ftr_rerun.R", "imputation"), c("04b_ftr_rerun.R", "complete_case"),
              c("06_survival_analysis.R", "imputation"), c("06_survival_analysis.R", "complete_case"),
              c("08_sensitivity_analyses.R", "imputation"),
              c("09_tables_figures.R", "imputation"), c("09_tables_figures.R", "complete_case"),
              c("10_compare_missing_data_methods.R", "imputation"))
if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("caffeinate")))
  system2("caffeinate", c("-i", "-w", Sys.getpid()), wait = FALSE)
for (st in steps) {
  Sys.setenv(ARCH_MISSING = st[2])
  message("\n==== ", st[1], "  [", st[2], "] ====")
  rc <- system2(rscript, shQuote(st[1]))
  if (rc != 0) stop(st[1], " [", st[2], "] failed. See results/logs/.")
}
message("Re-run complete.")

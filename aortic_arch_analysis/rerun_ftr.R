# =============================================================================
# rerun_ftr.R — re-runs only the failure-to-rescue analysis for both missing-data
# methods, then the tables/figures (09) and the MI vs CC comparison (10).
# Run after the full pipeline has finished:   source("rerun_ftr.R")
# =============================================================================

if (!file.exists("R/00_config.R")) stop("Run this script from the project folder (aortic_arch_analysis/).")
rscript <- file.path(R.home("bin"), "Rscript")
steps <- list(c("04b_ftr_rerun.R", "imputation"), c("04b_ftr_rerun.R", "complete_case"),
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
message("FTR re-run complete.")

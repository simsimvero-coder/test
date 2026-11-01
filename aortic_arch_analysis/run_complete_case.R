# =============================================================================
# run_complete_case.R — runs the whole pipeline with the complete-case analysis
# only (no multiple imputation). Results: results/complete_case/.
# Usage, in R:   setwd(dirname(file.choose())); source("run_complete_case.R")
# Resume after a stop:   START <- "05"; source("run_complete_case.R")
# =============================================================================
Sys.setenv(ARCH_METHODS = "complete_case", ARCH_MISSING = "complete_case")
source("run_all.R")

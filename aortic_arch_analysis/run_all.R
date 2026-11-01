# =============================================================================
# run_all.R — runs the complete pipeline, each script in a fresh R process
# (clean session per step), stopping at the first failure.
#
# Usage (from this folder):   Rscript run_all.R
# Optional environment variables: ARCH_RAW_DATA, ARCH_RESULTS, ARCH_DERIVED,
# ARCH_M (imputations), ARCH_NSIM, ARCH_CORES, ARCH_GLMM_FALLBACK (see R/00_config.R).
# Scripts 04-09 run twice (multiple imputation and complete case); script 10 compares them.
#
# Resume after a stop (earlier results are kept and reused):
#   START <- "05"; source("run_all.R")              # from 05 [imputation] onwards
#   START <- "04 complete_case"; source("run_all.R") # from 04 [complete case] onwards
# =============================================================================

if (!file.exists("R/00_config.R")) stop("Run this script from the project folder (aortic_arch_analysis/).")

required <- c("data.table", "ggplot2", "lme4", "survival", "splines", "mice", "ordinal",
              "clubSandwich", "MASS", "patchwork", "scales", "writexl", "readxl", "haven", "parallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the missing R packages first: ", paste(missing, collapse = ", "))

shared   <- c("00_check_environment.R", "01_data_validation.R", "02_cohort_descriptives.R", "03_volume_trends.R",
              "07_population_survival.R")
adjusted <- c("04_clinical_outcomes.R", "05_organ_protection.R", "06_survival_analysis.R",
              "08_sensitivity_analyses.R", "09_tables_figures.R")
# Adjusted analyses run twice: multiple imputation and complete-case analysis
# ARCH_METHODS selects the methods; "complete_case" alone skips imputation (see run_complete_case.R)
methods <- trimws(strsplit(Sys.getenv("ARCH_METHODS", "imputation,complete_case"), ",")[[1]])
methods <- intersect(c("imputation", "complete_case"), methods)
if (!length(methods)) stop("ARCH_METHODS must contain imputation and/or complete_case")
steps <- c(lapply(shared, function(s) list(script = s, method = methods[1])),
           unlist(lapply(methods, function(m)
             lapply(adjusted, function(s) list(script = s, method = m))), recursive = FALSE),
           if (length(methods) == 2) list(list(script = "10_compare_missing_data_methods.R", method = "imputation")))

# Optional resume point (see top); removed afterwards so a later run starts from the beginning
start_at <- if (exists("START", envir = globalenv())) get("START", envir = globalenv()) else ""
if (exists("START", envir = globalenv())) rm("START", envir = globalenv())
first <- 1L
if (nzchar(start_at)) {
  parts <- strsplit(trimws(start_at), "[ ]+")[[1]]
  hit <- which(startsWith(vapply(steps, `[[`, "", "script"), parts[1]) &
               (length(parts) < 2 | vapply(steps, `[[`, "", "method") == parts[2]))
  if (!length(hit)) stop("START = '", start_at, "' does not match a script (e.g. \"05\" or \"04 complete_case\").")
  first <- hit[1]
  if (!file.exists(file.path(Sys.getenv("ARCH_DERIVED", "data/derived"), "cohort.rds")))
    stop("No earlier results found in this folder: run without START first.")
  message("Resuming at ", steps[[first]]$script, " [", steps[[first]]$method, "]; earlier steps are not re-run.")
}

# macOS: keep the computer awake while the pipeline runs (stops automatically afterwards)
if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("caffeinate")))
  system2("caffeinate", c("-i", "-w", Sys.getpid()), wait = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
timing <- data.frame(script = vapply(steps, `[[`, "", "script"),
                     method = vapply(steps, `[[`, "", "method"), seconds = NA_real_, status = NA_character_)
timing$method[timing$script %in% c(shared, "10_compare_missing_data_methods.R")] <- "-"
for (i in seq_along(steps)) {
  st <- steps[[i]]
  if (i < first) { timing$status[i] <- "skipped (START)"; next }
  Sys.setenv(ARCH_MISSING = st$method)      # inherited by the child R process
  message("\n==== ", st$script, if (st$script %in% adjusted) paste0("  [", st$method, "]") else "", " ====")
  t0 <- Sys.time()
  rc <- system2(rscript, shQuote(st$script))   # system2 quotes the command itself
  timing$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  timing$status[i] <- if (rc == 0) "ok" else paste("failed (exit", rc, ")")
  if (rc != 0) {
    print(timing)
    stop("Pipeline stopped: ", st$script, " failed. See results/logs/.\nAfter a fix, resume with:  START <- \"",
         substr(st$script, 1, 2), if (st$script %in% adjusted) paste0(" ", st$method) else "", "\"; source(\"run_all.R\")")
  }
}
print(timing)
res <- Sys.getenv("ARCH_RESULTS", "results")
write.csv(timing, file.path(res, "logs", "pipeline_timing.csv"), row.names = FALSE)
message("Pipeline complete. Reports: ", paste(file.path(res, methods, "analysis_report.md"), collapse = ", "),
        if (length(methods) == 2) paste(" and", file.path(res, "comparison_missing_data", "MI_vs_CC_summary.md")) else "")

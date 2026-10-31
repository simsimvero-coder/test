# =============================================================================
# 00_check_environment.R
# Run this first on the computer that will do the analysis (e.g. Intel Mac with
# R 4.2.3). It checks the R version, the required packages and their versions,
# whether the data file can be found and read, whether the expected columns are
# present, and whether model fitting and parallel processing work.
# It does not change anything. Output: console and results/logs/environment_check.txt
#
# Usage, from this folder:   Rscript 00_check_environment.R
# Install missing packages (CRAN binaries for your R version):
#   install.packages(c("data.table","ggplot2","lme4","survival","mice","ordinal",
#                      "clubSandwich","MASS","patchwork","scales","writexl","readxl","haven"))
# =============================================================================

source("R/00_config.R")
res <- character(0)
say <- function(status, what, detail = "") {
  line <- sprintf("[%-4s] %-45s %s", status, what, detail)
  message(line); res <<- c(res, line)
}

# ---- R and system ----------------------------------------------------------------------
rv <- getRversion()
say(if (rv >= "4.2.0") "OK" else "FAIL", "R version", paste(rv, "(minimum 4.2.0)"))
say("INFO", "Platform", paste(R.version$platform, "|", Sys.info()[["sysname"]], Sys.info()[["release"]]))
say("INFO", "CPU cores detected", parallel::detectCores())
say("INFO", "Cores used for model fits (CONFIG$cores)", CONFIG$cores)

# ---- packages ---------------------------------------------------------------------------
# Minimum versions: functions used by this pipeline exist from these versions on.
need <- c(data.table = "1.13.0", ggplot2 = "3.3.0", lme4 = "1.1-27", survival = "3.2-7",
          mice = "3.13.0", ordinal = "2019.12-10", clubSandwich = "0.5.0", MASS = "7.3-50",
          patchwork = "1.1.0", scales = "1.1.0", writexl = "1.3", readxl = "1.3.1", haven = "2.3.0",
          nloptr = "1.2.0", Matrix = "1.3-0")
ok_pkgs <- TRUE
for (p in names(need)) {
  if (!requireNamespace(p, quietly = TRUE)) { say("FAIL", paste("package", p), "not installed"); ok_pkgs <- FALSE; next }
  v <- utils::packageVersion(p)
  ok <- v >= need[[p]]
  say(if (ok) "OK" else "FAIL", paste("package", p), sprintf("%s (minimum %s)", v, need[[p]]))
  ok_pkgs <- ok_pkgs && ok
}

if (ok_pkgs) {
  source("R/functions.R")
  # ---- data file ---------------------------------------------------------------------------
  path <- CONFIG$raw_data_path
  if (file.exists(path)) {
    say("OK", "Data file found", path)
    d <- tryCatch(read_raw(path), error = function(e) { say("FAIL", "Data file readable", conditionMessage(e)); NULL })
    if (!is.null(d)) {
      say("OK", "Data file read", sprintf("%d rows, %d columns", nrow(d), ncol(d)))
      miss <- setdiff(EXPECTED_VARS, names(d))
      essential <- c("centrum", "jaar", "leeftijd", "geslacht", "urgentie", "incl_cohort",
                     "overl30d", "mort_status", "overltijd_dag")
      status <- if (any(essential %in% miss)) "FAIL" else if (length(miss)) "WARN" else "OK"
      say(status, "Expected columns present",
          if (length(miss)) paste("missing:", paste(miss, collapse = ", ")) else "all 60 present")
    }
  } else {
    say("FAIL", "Data file found", paste("not found; searched:", paste(attr(path, "searched"), collapse = " ; ")))
  }

  # ---- model fitting (lme4 example data, no patient data) ----------------------------------
  fit <- tryCatch(fit_mixed(cbind(incidence, size - incidence) ~ period + (1 | herd),
                            lme4::cbpp, binomial()), error = function(e) NULL)
  say(if (!is.null(fit) && !is.null(fit$fit)) "OK" else "FAIL", "Mixed-model fitting (glmer, nloptwrap)",
      if (!is.null(fit)) paste("optimizer:", fit$optimizer) else "")
  # ---- parallel ---------------------------------------------------------------------------
  pr <- tryCatch(par_lapply(as.list(1:4), function(i) i^2), error = function(e) NULL)
  say(if (identical(unlist(pr), c(1, 4, 9, 16))) "OK" else "WARN", "Parallel processing",
      if (.Platform$OS.type == "windows") "Windows: runs serially" else "forked (mclapply)")
  # ---- life table and population files -------------------------------------------------------
  say(if (file.exists(file.path(CONFIG$cbs_lifetable_dir, "Observations.csv"))) "OK" else "WARN",
      "CBS life table (37360ned)", CONFIG$cbs_lifetable_dir)
  say(if (file.exists(CONFIG$population_path)) "OK" else "WARN", "CBS population (03759ned)", CONFIG$population_path)
}

dir.create(file.path(CONFIG$results_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("Environment check", format(Sys.time())), res, "", capture.output(sessionInfo())),
           file.path(CONFIG$results_dir, "logs", "environment_check.txt"))
if (any(grepl("^\\[FAIL", res))) {
  message("\nOne or more checks FAILED - fix these before running run_all.R.")
  quit(save = "no", status = 1)
}
message("\nAll essential checks passed.")

# =============================================================================
# run_nieuw.R - runs the three new standalone scripts in one go
# (after run_all.R has created data/derived; run_all.R does not need to be rerun).
#   centre_perfusion.R    perfusion practice and volume per centre
#   sens_strict_groups.R  paper 1, sensitivity analysis acute vs strictly elective
#   volume.R              paper 2, volume-outcome with strict groups and Davies test
# Use in RGui: File > Source R code... > choose this file, or type
#   source(file.choose())
# =============================================================================

# 1. Project folder: use the current folder if it is the project, otherwise ask for it
if (!file.exists(file.path("R", "00_config.R"))) {
  map <- if (.Platform$OS.type == "windows") utils::choose.dir(caption = "Kies de map aortic_arch_analysis") else
    readline("Pad naar de map aortic_arch_analysis: ")
  if (is.na(map) || !file.exists(file.path(map, "R", "00_config.R")))
    stop("Dit is niet de projectmap: R/00_config.R ontbreekt in ", map)
  setwd(map)
}
cat("Projectmap:", getwd(), "\n")

# 2. Settings (complete case is the primary analysis)
Sys.setenv(ARCH_METHODS = "complete_case", ARCH_MISSING = "complete_case", ARCH_VOLUME_BOOT = "500")

# 3. Data. The scripts do not read the Excel file themselves: they read the cleaned
#    cohort (data/derived/cohort.rds) that 01_data_validation.R and
#    02_cohort_descriptives.R make from the raw file named in R/00_config.R
#    ("AORTABOOG_NL_schoon v3(1).xlsx", found on the desktop, or the path in
#    Sys.setenv(ARCH_RAW_DATA = "...")). If that cohort is missing, it is made now.
source(file.path("R", "00_config.R"))
cat("Ruwe dataset volgens 00_config.R:", CONFIG$raw_data_path, "\n")
if (!file.exists(file.path(CONFIG$derived_dir, "cohort.rds"))) {
  cat("data/derived/cohort.rds ontbreekt: 01 en 02 worden nu gedraaid (leest de Excel in)\n")
  source("01_data_validation.R")
  source("02_cohort_descriptives.R")
}
cat("Gebruikt cohort: data/derived/cohort.rds, gemaakt op",
    format(file.mtime(file.path(CONFIG$derived_dir, "cohort.rds")), "%Y-%m-%d %H:%M"), "\n")

# 4. Package for the Davies test (volume.R skips only that test if it cannot be installed)
if (!requireNamespace("segmented", quietly = TRUE))
  try(install.packages("segmented"), silent = TRUE)

# 5. Run the scripts; one failing script does not stop the others
scripts <- c("centre_perfusion.R", "sens_strict_groups.R", "volume.R")
status <- setNames(character(length(scripts)), scripts)
for (s in scripts) {
  cat("\n\n#############  ", s, "  #############\n")
  t0 <- Sys.time()
  status[s] <- tryCatch({ source(s); "OK" },
                        error = function(e) paste("FOUT:", conditionMessage(e)))
  while (sink.number() > 0) sink()
  cat("\n", s, ": ", status[s], " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)\n", sep = "")
}

# 6. Where the overviews are
cat("\n\n============ KLAAR ============\n")
print(status)
cat("\nOverzichten (stuur hiervan screenshots):\n")
for (f in list.files("results", pattern = "^overzicht_(centre_perfusion|strict_groups|volume)\\.txt$",
                     recursive = TRUE, full.names = TRUE))
  cat("  ", normalizePath(f), "\n")

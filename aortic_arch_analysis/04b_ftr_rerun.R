# =============================================================================
# 04b_ftr_rerun.R
# Re-runs only the failure-to-rescue part of 04_clinical_outcomes.R (section 7) and
# updates ftr_models.csv and the FTR entries in outcomes_<method>.rds.
# Use this after a fix to the FTR code, so the rest of 04 (hours with imputation)
# does not have to be re-run. The FTR code itself is read from 04_clinical_outcomes.R,
# so both scripts always run the same analysis.
# Usage: run rerun_ftr.R (both missing-data methods, then 09 and 10).
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("04b_ftr_rerun")
log_msg("Missing-data method: ", method_label())
SUB <- "04_outcomes"

cohort <- load_derived("cohort.rds")
imps   <- analysis_datasets()
specs  <- load_derived("specs.rds")

src <- readLines("04_clinical_outcomes.R")
i0 <- grep("^# ---- 7. Failure to rescue", src)
i1 <- grep('^write_out\\(rbindlist\\(ftr_models\\), SUB, "ftr_models"\\)', src)
stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0)
eval(parse(text = src[i0:i1]))

f <- derived_method("outcomes.rds")
if (file.exists(f)) {
  outc <- readRDS(f)
  outc$ftr <- rbindlist(ftr_models); outc$ftr_decomp <- decomp; outc$ftr_decomp_all <- decomp_all
  saveRDS(outc, f)
}
log_msg("FTR models: ", nrow(rbindlist(ftr_models)[term == "year_c"]), " year-trend rows written.")
log_msg("Done.")

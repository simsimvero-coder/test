# =============================================================================
# 00_config.R — project configuration
# All paths, codings (NHR Handboek Cardiochirurgie v22.2.1) and analysis
# settings live here. Edit this file, not the analysis scripts.
# Paths can be overridden with environment variables (used by the test run).
# =============================================================================

# ---- location of the patient-level file --------------------------------------
# The extract is read directly from the desktop (read-only; never copied or
# modified). Windows: R's "~" points to Documents, so the user profile is used,
# and desktops redirected to OneDrive ("OneDrive - <organisation>") are searched too.
RAW_FILE_NAME  <- "AORTABOOG_NL_schoon v3(1).xlsx"
# Excel sheet with the data. "auto" = the sheet containing the most expected
# variables (preferring "Alle_records", which also holds the excluded records and
# incl_cohort/excl_reden). A sheet name or number can be given instead.
RAW_FILE_SHEET <- "auto"

find_desktop_file <- function(file_name) {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  roots <- c(home, Sys.getenv("OneDrive"), Sys.getenv("OneDriveCommercial"),
             list.files(home, pattern = "^OneDrive", full.names = TRUE),
             # macOS: OneDrive under ~/Library/CloudStorage; iCloud-synced Desktop
             list.files(file.path(home, "Library", "CloudStorage"), pattern = "^OneDrive", full.names = TRUE),
             file.path(home, "Library", "Mobile Documents", "com~apple~CloudDocs"))
  roots <- unique(roots[nzchar(roots)])
  candidates <- unlist(lapply(roots, function(r) file.path(r, c("Desktop", "Bureaublad"), file_name)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(normalizePath(hit[1], winslash = "/"))
  # not found: return the first candidate so the error message shows where we looked
  out <- candidates[1]
  attr(out, "searched") <- candidates
  out
}

CONFIG <- list(
  # ---- input ----------------------------------------------------------------
  # Patient-level extract on the desktop (see above). Can be overridden with
  # the environment variable ARCH_RAW_DATA. Read-only.
  raw_data_path   = { p <- Sys.getenv("ARCH_RAW_DATA"); if (nzchar(p)) p else find_desktop_file(RAW_FILE_NAME) },
  raw_data_sheet  = RAW_FILE_SHEET,
  # CBS StatLine 37360ned (life tables) as downloaded (folder with Observations.csv)
  cbs_lifetable_dir = Sys.getenv("ARCH_CBS_DIR", "data/external/cbs_37360ned"),
  # CBS StatLine 03759ned (population on 1 January and average population), as
  # exported ("TabelIndeling" csv). Optional: volume rates are skipped if absent.
  population_path = Sys.getenv("ARCH_POPULATION",
                               "data/external/cbs_03759ned/03759ned_TabelIndeling.csv"),

  # ---- output ---------------------------------------------------------------
  results_dir = Sys.getenv("ARCH_RESULTS", "results"),
  derived_dir = Sys.getenv("ARCH_DERIVED", "data/derived"),

  # ---- general --------------------------------------------------------------
  seed            = 20261029L,
  study_years     = 2013:2024,
  reference_year  = 2018,          # calendar year is centred here
  era_breaks      = c(2012, 2016, 2020, 2024),
  era_labels      = c("2013-2016", "2017-2020", "2021-2024"),

  # values of incl_cohort that mean "included" (compared as character)
  inclusion_values = c("1", "TRUE", "True", "true", "ja", "Ja", "JA", "yes", "Yes"),
  use_clean_versions    = TRUE,    # klepchir_clean / aortachir_clean / circ_arrest_clean
  drop_exact_duplicates = FALSE,   # duplicates are always reported; never dropped silently
  strict_coding = TRUE,            # stop if a coded variable holds values outside the handbook

  # ---- surgical cohort selection (re-applied explicitly in 02) ----------------
  # Open or hybrid arch surgery with cardiopulmonary bypass (ECC 10/20/30);
  # excluded: endovascular-only aortic repair (aotachir_benadering 20/21),
  # percutaneous approach (chir_benadering 50), chronic dissection (code 30) or
  # infection (code 40) in any aortic segment. Records with unknown aortic approach
  # are KEPT (as in the original analysis) when they meet all other criteria;
  # set keep_unknown_approach = FALSE to exclude them. Counted in the flow diagram.
  keep_unknown_approach = TRUE,

  # ---- outcome definitions (paper 1 definitions are the main analysis) --------
  # Stroke counted in the primary composite (paper 1: with or without deficit):
  # "deficit" = CVA_restletsel only; "any" = CVA_restletsel or CVA_zonder_restletsel
  # (the latter includes TIA). The other choice is run as a sensitivity analysis.
  composite_stroke = "any",
  # Stroke counted in 'major complication' (and hence the FTR denominator).
  # Paper 1: stroke with residual deficit only.
  major_stroke = "deficit",
  # 1-year mortality (binary) is analysed for operation years up to this year only,
  # so that every included year has complete potential 1-year follow-up.
  d1y_last_year = 2023,
  # Sensitivity analysis: every temporal trend is re-estimated from this year on
  trend_sens_start = 2015,

  # ---- missing data ---------------------------------------------------------
  # Adjusted analyses (scripts 04, 05, 06, 08, 09) are run once per method;
  # run_all.R runs both. "imputation" = multiple imputation (mice);
  # "complete_case" = per model, patients with all model variables observed.
  missing_method = Sys.getenv("ARCH_MISSING", "imputation"),
  n_imputations = as.integer(Sys.getenv("ARCH_M", "20")),
  mnar_var      = "euroscore2",    # covariate shifted in the MNAR sensitivity analysis
  mnar_delta    = 2,               # shift of imputed values (EuroSCORE II, percentage points)

  # cores used to fit models on imputed data sets in parallel (forked; Linux/macOS)
  cores = as.integer(Sys.getenv("ARCH_CORES", max(1L, parallel::detectCores() - 1L))),

  # Last-resort refit for GLMMs that still show a convergence warning after
  # nloptwrap restarts: "bobyqa" (slow, most rigorous; default) or "none" (keep
  # the best restart, reported as a convergence warning in the diagnostics).
  glmm_fallback = Sys.getenv("ARCH_GLMM_FALLBACK", "bobyqa"),

  # ---- simulation-based intervals --------------------------------------------
  n_sim = as.integer(Sys.getenv("ARCH_NSIM", "500")),

  # ---- survival ---------------------------------------------------------------
  km_times_days  = c(30, 365.25, 3 * 365.25, 5 * 365.25, 10 * 365.25),
  landmark_days  = 30,
  piecewise_cuts = c(30, 365.25),
  rmst_tau_years = c(5, 10),

  # ---- reporting ----------------------------------------------------------------
  fig_width = 7, fig_height = 4.5, fig_dpi = 300
)

# Unknown code used throughout the NHR handbook
NHR_UNKNOWN <- -1

# -----------------------------------------------------------------------------
# Allowed codes per variable (NHR handbook v22.2.1, earlier versions where they
# differ: pathology 'overig' was 90 before v22.2.0 and 91 afterwards).
# -1 (unknown) is allowed everywhere it appears in the handbook and is set to NA.
# -----------------------------------------------------------------------------
BIN <- c(0, 1, -1)
CODES <- list(
  geslacht              = c(1, 2, -1),
  cardiochir_eerder     = BIN,
  endocarditis          = BIN,
  urgentie              = c(10, 20, 30, 40, -1),
  ECC                   = c(0, 10, 20, 30, -1),
  circ_arrest           = BIN,
  circ_arrest_clean     = BIN,
  chir_benadering       = c(10, 20, 30, 40, 50, 60, -1),
  CABG                  = BIN,
  klepchir              = BIN,
  klepchir_clean        = BIN,
  aklepchir             = BIN,
  mklepchir             = BIN,
  pklepchir             = BIN,
  tklepchir             = BIN,
  aortachir             = BIN,
  aortachir_clean       = BIN,
  aotachir_benadering   = c(0, 10, 20, 21, 30, -1),
  aorta_asc             = c(0, 10, 20, 30, 90, -1),
  aorta_asc_path        = c(0, 10, 20, 30, 40, 90, 91, -1),
  aortaboog_path        = c(0, 10, 20, 30, 40, 90, 91, -1),
  aorta_desc            = c(0, 1, 10, 20, 30, -1),
  aorta_desc_path       = c(0, 10, 20, 30, 40, 90, 91, -1),
  cardiochir_overig     = BIN,
  overleden_zh          = BIN,
  ECMO_opname           = BIN,
  resp_insuff           = BIN,
  beademing             = BIN,
  heropname_IC          = BIN,
  CVA_restletsel        = BIN,
  CVA_zonder_restletsel = BIN,
  nierfalen             = BIN,
  gastroint_compl       = c(0, 1, 10, 90, -1),
  vasc_compl_opname     = c(0, 10, 20, -1),
  mort_status           = BIN,
  overl30d              = BIN,
  overl120d             = BIN,
  overl1j               = BIN,
  delier                = BIN,
  postop_MI             = BIN
)

# Continuous variables: handbook ranges (outside -> NA, logged).
# 'hard' = physically impossible / outside handbook; flagged values are kept.
RANGES <- list(
  leeftijd         = c(18, 100),
  BMI              = c(12, 70),
  kreatinine       = c(1, 2000),
  eGFR             = c(1, 200),
  LVEF             = c(1, 99),
  EuroscoreII      = c(0.49, 100),
  min_temp         = c(10, 40),
  ECC_duur         = c(0, 1000),
  AoX_duur         = c(0, 600),
  circ_arrest_duur = c(0, 300),   # handbook max 150; 150-300 kept but flagged
  opnameduur       = c(0, 365),
  overltijd_dag    = c(0, 20 * 366)
)
FLAG_ABOVE <- list(circ_arrest_duur = 150)

# Variables expected in the extract (from the delivered variable list)
EXPECTED_VARS <- c(
  "bronrij", "centrum", "jaar", "leeftijd", "geslacht", "BMI", "kreatinine", "eGFR",
  "LVEF", "cardiochir_eerder", "endocarditis", "urgentie", "EuroscoreII", "ECC",
  "circ_arrest", "chir_benadering", "CABG", "klepchir", "aklepchir", "mklepchir",
  "pklepchir", "tklepchir", "aortachir", "aotachir_benadering", "aorta_asc",
  "aorta_asc_path", "aortaboog_path", "aorta_desc", "aorta_desc_path",
  "cardiochir_overig", "min_temp", "ECC_duur", "AoX_duur", "circ_arrest_duur",
  "overleden_zh", "ECMO_opname", "opnameduur", "resp_insuff", "beademing",
  "heropname_IC", "CVA_restletsel", "CVA_zonder_restletsel", "nierfalen",
  "gastroint_compl", "vasc_compl_opname", "mort_status", "overltijd_dag", "overl30d",
  "overl120d", "overl1j", "delier", "postop_MI", "klepchir_clean", "aortachir_clean",
  "circ_arrest_clean", "pathologie_groep", "s1_via_bewijs", "excl_reden",
  "incl_cohort", "klasse"
)

# Derived variables whose definitions are not in the handbook (described only)
UNDOCUMENTED_VARS <- c("bronrij", "pathologie_groep", "s1_via_bewijs", "excl_reden",
                       "incl_cohort", "klasse", "overl120d", "klepchir_clean",
                       "aortachir_clean", "circ_arrest_clean", "eGFR", "BMI",
                       "opnameduur", "overltijd_dag")

# Outcomes added in handbook v22.2.0 (operations from 1 July 2022)
LATE_ADDED_VARS <- c("ECMO_opname", "postop_MI", "delier")

# Covariates of the primary model, as specified by the investigators: age (spline),
# sex, EuroSCORE II (logit) and previous cardiac surgery (redo). Used in every
# adjusted analysis. NB: EuroSCORE II itself contains age, sex and previous
# cardiac surgery; these are therefore adjusted for partly twice (see methods report).
COVARS_M1 <- c("age", "female", "logit_es2", "prior_cs")
# Reduced model when events per parameter < 10
COVARS_REDUCED <- c("age", "female", "logit_es2")
# Broader set of the earlier analysis plan, kept as a sensitivity analysis
COVARS_EXTENDED <- c("age", "female", "egfr", "lvef_cat", "prior_cs", "endocarditis_b",
                     "urg3", "path3")
COVARS_OPERATIVE <- c("conc_cabg", "conc_valve", "root_inv", "desc_inv", "hybrid_fet",
                      "other_cardiac")

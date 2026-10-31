# =============================================================================
# 01_data_validation.R
# Reads the raw extract (read-only), checks structure, codings, ranges,
# duplicates and logical consistency, and creates analysis variables.
# Nothing is dropped here: every record is kept and every recode is logged.
# Output: results/01_validation/*, data/derived/validated.rds
# =============================================================================

source("R/00_config.R"); source("R/functions.R")
start_log("01_data_validation")
SUB <- "01_validation"

raw <- read_raw(CONFIG$raw_data_path)
log_msg("Read ", nrow(raw), " rows x ", ncol(raw), " columns from ", CONFIG$raw_data_path)
n_raw <- nrow(raw)

# ---- 1. Structure ----------------------------------------------------------------
missing_vars <- setdiff(EXPECTED_VARS, names(raw))
extra_vars   <- setdiff(names(raw), EXPECTED_VARS)
write_out(data.table(variable = c(missing_vars, extra_vars),
                     issue = c(rep("expected but absent", length(missing_vars)),
                               rep("present but not in variable list", length(extra_vars)))),
          SUB, "structure_check")
if (length(missing_vars)) log_msg("WARNING: absent variables: ", paste(missing_vars, collapse = ", "))
essential <- c("centrum", "jaar", "leeftijd", "geslacht", "urgentie", "incl_cohort",
               "overl30d", "mort_status", "overltijd_dag")
if (any(essential %in% missing_vars)) {
  stop("Essential variables absent: ", paste(intersect(essential, missing_vars), collapse = ", "))
}

# ---- 2. Variable inventory --------------------------------------------------------
inv <- rbindlist(lapply(names(raw), function(v) {
  x <- raw[[v]]
  xc <- suppressWarnings(to_code(x))
  data.table(variable = v, class = class(x)[1], n = length(x),
             n_na = sum(is.na(x) | (is.character(x) & trimws(x) == "")),
             n_minus1 = sum(xc == NHR_UNKNOWN, na.rm = TRUE),
             n_unique = uniqueN(x),
             min = if (is.numeric(xc) && any(!is.na(xc))) min(xc, na.rm = TRUE) else NA_real_,
             max = if (is.numeric(xc) && any(!is.na(xc))) max(xc, na.rm = TRUE) else NA_real_,
             examples = paste(head(unique(as.character(x)), 6), collapse = " | "),
             handbook_codes = if (v %in% names(CODES)) paste(CODES[[v]], collapse = ",") else
                              if (v %in% names(RANGES)) paste(RANGES[[v]], collapse = "-") else "",
             undocumented = v %in% UNDOCUMENTED_VARS)
}))
write_out(inv, SUB, "variable_inventory")

# Frequencies of every categorical / coded variable (including undocumented ones)
cat_vars <- intersect(c(names(CODES), "centrum", "jaar", "pathologie_groep", "klasse",
                        "s1_via_bewijs", "excl_reden", "incl_cohort"), names(raw))
freqs <- rbindlist(lapply(cat_vars, function(v) {
  t <- raw[, .N, by = .(value = as.character(get(v)))]
  t[, variable := v]; setorder(t, value); t
}))
write_out(freqs[, .(variable, value, N)], SUB, "value_frequencies")

# ---- 3. Codings -------------------------------------------------------------------
dv <- copy(raw)
dv[, row_id := .I]                       # internal row index (bronrij is kept as is)
code_issues <- list()
for (v in intersect(names(CODES), names(dv))) {
  x <- to_code(dv[[v]])
  nonnum <- !is.na(dv[[v]]) & is.na(x) & trimws(as.character(dv[[v]])) != ""
  bad <- !is.na(x) & !(x %in% CODES[[v]])
  if (any(nonnum | bad)) {
    code_issues[[v]] <- data.table(variable = v,
                                   value = as.character(dv[[v]][nonnum | bad]))[, .N, by = .(variable, value)]
  }
  x[bad] <- NA                           # only reached if strict_coding = FALSE
  x[x == NHR_UNKNOWN] <- NA
  set(dv, j = v, value = x)
}
code_issues <- rbindlist(code_issues)
write_out(if (nrow(code_issues)) code_issues else data.table(variable = character(), value = character(), N = integer()),
          SUB, "coding_violations")
if (nrow(code_issues) && CONFIG$strict_coding) {
  stop("Values outside the NHR handbook codes were found (see results/01_validation/",
       "coding_violations.csv). Check whether the extract was recoded, then update CODES ",
       "in R/00_config.R or set strict_coding = FALSE (values then become NA, logged).")
}

# Continuous variables: -1 -> NA, range checks
range_log <- list()
for (v in intersect(names(RANGES), names(dv))) {
  x <- to_code(dv[[v]])
  n_m1 <- sum(x == NHR_UNKNOWN, na.rm = TRUE)
  if (!v %in% c("overltijd_dag")) x[x == NHR_UNKNOWN] <- NA
  lo <- RANGES[[v]][1]; hi <- RANGES[[v]][2]
  out <- !is.na(x) & (x < lo | x > hi)
  flag <- if (!is.null(FLAG_ABOVE[[v]])) sum(!is.na(x) & x > FLAG_ABOVE[[v]] & !out) else 0L
  range_log[[v]] <- data.table(variable = v, range = paste(lo, hi, sep = "-"),
                               n_unknown_code = n_m1, n_out_of_range_set_NA = sum(out),
                               out_of_range_values = paste(head(sort(x[out]), 10), collapse = ","),
                               n_flagged_kept = flag)
  x[out] <- NA
  set(dv, j = v, value = x)
}
write_out(rbindlist(range_log), SUB, "range_checks")

# ---- 4. Identifiers, year, centre --------------------------------------------------
dv[, jaar := as.integer(to_code(jaar))]
dv[, centrum := as.character(centrum)]
id_check <- data.table(
  check = c("bronrij unique", "bronrij missing", "year outside study period",
            "year missing", "number of centres", "centre missing"),
  value = c(uniqueN(dv$bronrij) == nrow(dv), sum(is.na(dv$bronrij)),
            sum(!is.na(dv$jaar) & !dv$jaar %in% CONFIG$study_years), sum(is.na(dv$jaar)),
            uniqueN(dv$centrum[!is.na(dv$centrum)]), sum(is.na(dv$centrum))))
write_out(id_check, SUB, "identifier_checks")

# ---- 5. Duplicates ------------------------------------------------------------------
# No patient identifier exists; duplicates are judged on content.
content_vars <- setdiff(intersect(EXPECTED_VARS, names(dv)), c("bronrij"))
dv[, dup_exact := duplicated(dv, by = content_vars) | duplicated(dv, by = content_vars, fromLast = TRUE)]
near_vars <- intersect(c("centrum", "jaar", "leeftijd", "geslacht", "BMI", "kreatinine",
                         "LVEF", "urgentie", "EuroscoreII", "overltijd_dag", "mort_status"), names(dv))
dv[, dup_near := duplicated(dv, by = near_vars) | duplicated(dv, by = near_vars, fromLast = TRUE)]
dup_summary <- data.table(
  check = c("rows in exact-duplicate groups (all variables except bronrij)",
            "exact duplicates beyond first occurrence",
            "rows in near-duplicate groups (centre, year, age, sex, BMI, creatinine, LVEF, urgency, EuroSCORE II, follow-up, status)"),
  n = c(sum(dv$dup_exact), sum(duplicated(dv, by = content_vars)), sum(dv$dup_near)))
write_out(dup_summary, SUB, "duplicates")
dv[, drop_duplicate := CONFIG$drop_exact_duplicates & duplicated(dv, by = content_vars)]
log_msg("Exact duplicates beyond first: ", dup_summary$n[2],
        if (CONFIG$drop_exact_duplicates) " (will be excluded)" else " (kept; drop_exact_duplicates = FALSE)")

# ---- 6. Analysis variables -----------------------------------------------------------
pick <- function(clean, raw_v) {
  if (CONFIG$use_clean_versions && clean %in% names(dv)) dv[[clean]] else dv[[raw_v]]
}
seg_any <- function(code) {
  m <- cbind(dv$aorta_asc_path, dv$aortaboog_path, dv$aorta_desc_path)
  apply(m, 1, function(r) if (all(is.na(r))) NA else any(r %in% code))
}

dv[, `:=`(
  included   = as.character(incl_cohort) %in% CONFIG$inclusion_values,
  year       = jaar,
  year_c     = jaar - CONFIG$reference_year,
  era        = cut(jaar, CONFIG$era_breaks, CONFIG$era_labels),
  age        = leeftijd,
  female     = as.integer(geslacht == 2),
  bmi        = BMI,
  creat      = kreatinine,
  egfr       = eGFR,
  lvef       = LVEF,
  euroscore2 = EuroscoreII,
  prior_cs   = cardiochir_eerder,
  endocarditis_b = endocarditis,
  urg4 = factor(urgentie, c(10, 20, 30, 40), c("Elective", "Urgent", "Emergency", "Salvage")),
  urg3 = factor(fcase(urgentie == 10, "Elective", urgentie == 20, "Urgent",
                      urgentie %in% c(30, 40), "Emergency/salvage"),
                c("Elective", "Urgent", "Emergency/salvage")),
  nonelective = as.integer(urgentie %in% c(20, 30, 40)),
  emerg_salv  = as.integer(urgentie %in% c(30, 40))
)]
dv[is.na(urgentie), `:=`(nonelective = NA_integer_, emerg_salv = NA_integer_)]
# (urgency2 = acute vs non-acute group; defined below, after pathology)
# EuroSCORE II on the logit scale (adjustment variable; EuroSCORE II is a predicted probability in %)
dv[, logit_es2 := qlogis(pmin(pmax(euroscore2 / 100, 1e-4), 1 - 1e-4))]
dv[, lvef_cat := cut(lvef, c(0, 20, 30, 50, 100), c("<=20%", "21-30%", "31-50%", ">50%"))]
dv[, lvef_cat := relevel(lvef_cat, ">50%")]

# Pathology (precedence: acute dissection > chronic dissection/aneurysm > infection > other)
dv[, `:=`(p_acute = seg_any(20), p_aneur = seg_any(c(10, 30)),
          p_infect = seg_any(40), p_other = seg_any(c(90, 91)),
          p_chronic = seg_any(30))]
dv[, path5 := fcase(p_acute, "Acute dissection",
                    p_aneur & seg_any(30), "Chronic dissection",
                    p_aneur, "Aneurysm",
                    p_infect, "Infection",
                    p_other, "Other",
                    !is.na(p_acute), "None recorded",
                    default = NA_character_)]
# Chronic dissection and infection are excluded from the analysis cohort (02);
# within the cohort, path3 is therefore aneurysm / acute dissection / other.
dv[, path3 := factor(fcase(path5 == "Acute dissection", "Acute dissection",
                           path5 %in% c("Aneurysm", "Chronic dissection"), "Aneurysm",
                           path5 %in% c("Infection", "Other", "None recorded"), "Other"),
                     c("Aneurysm", "Acute dissection", "Other"))]
dv[, acute_dissection := as.integer(path3 == "Acute dissection")]

# Acute vs non-acute (study definition, specified by the investigators):
#   Acute:       urgentie 30 (spoed) or 40 (redding), OR acute dissection (code 20)
#                in aorta_asc_path, aortaboog_path or aorta_desc_path.
#   Non-acute:   urgentie 10 (electief) or 20 (urgent) without an acute
#                classification, OR aneurysm (code 10) in any segment without an
#                acute classification, also when urgency is unknown.
#   Not classifiable: no acute dissection or aneurysm code AND no usable urgency code.
# Conflicts: acute dissection -> always acute (also with urgentie 10/20);
#            urgentie 30/40 -> always acute (also with aneurysm code 10).
# NB: urgentie 20 (urgent) belongs to the NON-acute group.
dv[, `:=`(seg_ad = seg_any(20) %in% TRUE, seg_aneur10 = seg_any(10) %in% TRUE)]
dv[, acute := fcase(urgentie %in% c(30, 40) | seg_ad, 1L,
                    urgentie %in% c(10, 20) | seg_aneur10, 0L,
                    default = NA_integer_)]
dv[, acute_class := factor(fcase(acute == 1, "Acute", acute == 0, "Non-acute", default = "Not classifiable"),
                           c("Non-acute", "Acute", "Not classifiable"))]
dv[, urgency2 := factor(fcase(acute == 1, "Acute", acute == 0, "Non-acute"), c("Non-acute", "Acute"))]

# Operative characteristics
dv[, `:=`(
  conc_cabg   = CABG,
  conc_valve  = pick("klepchir_clean", "klepchir"),
  aortic_valve = aklepchir,
  aorta_any   = pick("aortachir_clean", "aortachir"),
  root_inv    = as.integer(aorta_asc %in% c(10, 20)),
  asc_inv     = as.integer(aorta_asc > 0),
  desc_inv    = as.integer(aorta_desc > 0),
  hybrid_fet  = as.integer(aotachir_benadering == 30),
  other_cardiac = cardiochir_overig,
  sternotomy  = as.integer(chir_benadering == 10),
  ecc_any     = as.integer(ECC > 0),
  circ_arrest_b = pick("circ_arrest_clean", "circ_arrest"),
  temp_min    = min_temp,
  ecc_min     = ECC_duur, aox_min = AoX_duur, ca_min = circ_arrest_duur
)]
for (v in c("root_inv", "asc_inv")) dv[is.na(aorta_asc), (v) := NA_integer_]
dv[is.na(aorta_desc), desc_inv := NA_integer_]
dv[is.na(aotachir_benadering), hybrid_fet := NA_integer_]
dv[is.na(chir_benadering), sternotomy := NA_integer_]
dv[is.na(ECC), ecc_any := NA_integer_]

# Hypothermia (Yan et al. 2013 consensus)
dv[, hypo5 := cut(temp_min, c(-Inf, 14, 20, 24, 28, 34, Inf),
                  c("Profound (<=14)", "Deep (14.1-20)", "Low-moderate (20.1-24)",
                    "High-moderate (24.1-28)", "Mild (28.1-34)", "Normothermia (>34)"))]
dv[, hypo4 := cut(temp_min, c(-Inf, 20, 28, 34, Inf),
                  c("Deep (<=20)", "Moderate (20.1-28)", "Mild (28.1-34)", "Normothermia (>34)"),
                  ordered_result = TRUE)]

# Outcomes
dv[, `:=`(
  d30 = overl30d, dzh = overleden_zh, d1y_flag = overl1j,
  stroke_def = CVA_restletsel, stroke_nodef = CVA_zonder_restletsel,
  renal = nierfalen, vent24 = beademing, reintub = resp_insuff, icu_readm = heropname_IC,
  ecmo = ECMO_opname, pmi = postop_MI, delirium = delier,
  gi = as.integer(gastroint_compl %in% c(1, 10, 90)),
  vasc_major = as.integer(vasc_compl_opname == 10),
  los = opnameduur,
  time_d = overltijd_dag, status = mort_status
)]
dv[is.na(gastroint_compl), gi := NA_integer_]
dv[is.na(vasc_compl_opname), vasc_major := NA_integer_]

composite <- function(cols) {
  m <- as.matrix(dv[, .SD, .SDcols = cols])
  apply(m, 1, function(r) if (any(r == 1, na.rm = TRUE)) 1L else if (any(is.na(r))) NA_integer_ else 0L)
}
# Any stroke (with or without deficit), with the same strict missing-data rule as the
# composites: 1 if either is positive; 0 only if both are known and negative; else NA
dv[, stroke_any := composite(c("stroke_def", "stroke_nodef"))]
# Stroke used in the composites (CONFIG$composite_stroke; the other choice is a sensitivity analysis)
dv[, stroke_comp := if (CONFIG$composite_stroke == "any") stroke_any else stroke_def]
dv[, stroke_alt  := if (CONFIG$composite_stroke == "any") stroke_def else stroke_any]

# PRIMARY ENDPOINT (paper 1): in-hospital death, stroke, renal failure or GI complication
dv[, comp_primary := composite(c("dzh", "stroke_comp", "renal", "gi"))]
dv[, comp_primary_altstroke := composite(c("dzh", "stroke_alt", "renal", "gi"))]

# MAJOR COMPLICATION (paper 1, main analysis and FTR denominator): stroke with residual
# deficit (CONFIG$major_stroke), renal failure, GI complication
dv[, stroke_major := if (CONFIG$major_stroke == "any") stroke_any else stroke_def]
MAJOR_COMPONENTS <- c("stroke_major", "renal", "gi")
dv[, major_compl := composite(MAJOR_COMPONENTS)]
dv[, n_compl := rowSums(as.matrix(.SD) == 1, na.rm = TRUE), .SDcols = MAJOR_COMPONENTS]
# Supplementary definitions: broader set (adds reintubation and major vascular complication),
# and extended set (adds ventilation > 24 h and stroke without deficit)
MAJOR_NEW <- c("stroke_def", "renal", "reintub", "gi", "vasc_major")
dv[, major_compl_new := composite(MAJOR_NEW)]
dv[, n_compl_new := rowSums(as.matrix(.SD) == 1, na.rm = TRUE), .SDcols = MAJOR_NEW]
dv[, major_compl_ext := composite(c(MAJOR_NEW, "vent24", "stroke_nodef"))]

# 1-year mortality: binary only when follow-up >= 365 days or death within 365 days
dv[, d1y := fifelse(!is.na(time_d) & !is.na(status),
                    fifelse(status == 1 & time_d <= 365, 1L,
                            fifelse(time_d >= 365, 0L, NA_integer_)), NA_integer_)]

# ---- 7. Logical consistency --------------------------------------------------------------
cc <- function(label, cond) data.table(check = label, n = sum(cond, na.rm = TRUE))
consistency <- rbindlist(list(
  cc("overl30d = 1 but overltijd_dag > 30", dv$d30 == 1 & dv$time_d > 30),
  cc("overl30d = 0 but died with overltijd_dag <= 30", dv$d30 == 0 & dv$status == 1 & dv$time_d <= 30),
  cc("overl30d = 1 but mort_status = 0", dv$d30 == 1 & dv$status == 0),
  cc("alive (mort_status = 0) with overltijd_dag < 30", dv$status == 0 & dv$time_d < 30),
  cc("overleden_zh = 1 but mort_status = 0", dv$dzh == 1 & dv$status == 0),
  cc("overleden_zh = 1 but overltijd_dag > opnameduur + 1", dv$dzh == 1 & dv$time_d > dv$los + 1),
  cc("overleden_zh = 1, opnameduur <= 30 but overl30d = 0", dv$dzh == 1 & dv$los <= 30 & dv$d30 == 0),
  cc("alive at discharge but overltijd_dag < opnameduur", dv$dzh == 0 & dv$time_d < dv$los),
  cc("overl1j = 1 but overltijd_dag > 366", dv$d1y_flag == 1 & dv$time_d > 366),
  cc("overl1j = 0 but died with overltijd_dag <= 365", dv$d1y_flag == 0 & dv$status == 1 & dv$time_d <= 365),
  cc("overl1j = 0 but follow-up < 365 days (alive)", dv$d1y_flag == 0 & dv$status == 0 & dv$time_d < 365),
  cc("overl120d = 1 but overltijd_dag > 120", dv$overl120d == 1 & dv$time_d > 120),
  cc("overltijd_dag negative", dv$time_d < 0),
  cc("overltijd_dag = 0", dv$time_d == 0),
  cc("circulatory arrest = 0 but arrest duration > 0", dv$circ_arrest_b == 0 & dv$ca_min > 0),
  cc("circulatory arrest = 1 but arrest duration = 0", dv$circ_arrest_b == 1 & dv$ca_min == 0),
  cc("no ECC but circulatory arrest = 1", dv$ecc_any == 0 & dv$circ_arrest_b == 1),
  cc("no ECC but ECC duration > 0", dv$ecc_any == 0 & dv$ecc_min > 0),
  cc("cross-clamp time > ECC time", dv$aox_min > dv$ecc_min),
  cc("both CVA with and without deficit = 1", dv$stroke_def == 1 & dv$stroke_nodef == 1),
  cc("aortaboog_path = 0 (no arch pathology) among included", dv$included & dv$aortaboog_path == 0),
  cc("aotachir_benadering = 0 (no aortic surgery coded) among included; not excluded", dv$included & dv$aotachir_benadering == 0),
  cc("aortaboog_path missing among included", dv$included & is.na(dv$aortaboog_path)),
  cc("included but excl_reden filled", dv$included & !is.na(dv$excl_reden) & trimws(as.character(dv$excl_reden)) != ""),
  cc("temperature recorded > 34 C with circulatory arrest", dv$circ_arrest_b == 1 & dv$temp_min > 34)
))
write_out(consistency, SUB, "consistency_checks")

# Acute classification: how each group arises (documents the conflict rules)
write_out(dv[, .N, by = .(acute_class, urgentie, acute_dissection_any_segment = seg_ad,
                          aneurysm10_any_segment = seg_aneur10)][order(acute_class, urgentie)],
          SUB, "acute_classification_detail")

# Raw vs clean versions
clean_pairs <- list(c("klepchir", "klepchir_clean"), c("aortachir", "aortachir_clean"),
                    c("circ_arrest", "circ_arrest_clean"))
clean_tab <- rbindlist(lapply(clean_pairs, function(p) {
  if (!all(p %in% names(dv))) return(NULL)
  t <- dv[, .N, by = .(raw = get(p[1]), clean = get(p[2]))]
  t[, pair := paste(p, collapse = " vs ")]; t
}))
write_out(clean_tab, SUB, "raw_vs_clean")

# Derived pathology vs supplied pathologie_groep; klasse description
if ("pathologie_groep" %in% names(dv))
  write_out(dv[, .N, by = .(path5_derived = path5, pathologie_groep)][order(path5_derived)], SUB, "pathology_derived_vs_supplied")
if ("klasse" %in% names(dv))
  write_out(dv[, .N, by = .(klasse, included)][order(klasse)], SUB, "klasse_description")
write_out(dv[, .N, by = .(incl_cohort, excl_reden, s1_via_bewijs)][order(-N)], SUB, "inclusion_exclusion")

# Structural missingness by year (included records): flags variables not recorded in some years
miss_vars <- intersect(EXPECTED_VARS, names(dv))
miss_year <- dv[included == TRUE, lapply(.SD, function(x) round(100 * mean(is.na(x)), 1)),
                by = jaar, .SDcols = miss_vars]
setorder(miss_year, jaar)
write_out(miss_year, SUB, "missing_pct_by_year_all_variables")
struct <- melt(miss_year, id.vars = "jaar", variable.name = "variable", value.name = "pct_missing")
struct_flag <- struct[, .(min_pct = min(pct_missing), max_pct = max(pct_missing),
                          years_ge90pct_missing = paste(jaar[pct_missing >= 90], collapse = ","),
                          first_year_lt90pct = suppressWarnings(min(jaar[pct_missing < 90]))),
                      by = variable][max_pct >= 90 & min_pct < 50]
write_out(struct_flag, SUB, "structurally_missing_variables")
STRUCTURAL <- struct_flag[, .(variable = as.character(variable), first_year = first_year_lt90pct)]

# ---- 8. Save -------------------------------------------------------------------------------
attr(dv, "structural_missing") <- STRUCTURAL
attr(dv, "major_components") <- MAJOR_COMPONENTS
saveRDS(dv, derived_path("validated.rds"))
stopifnot(nrow(dv) == n_raw)

# ---- 9. Report -------------------------------------------------------------------------------
md <- c("# Data validation report", "",
        sprintf("Generated %s from `%s`.", format(Sys.time()), basename(CONFIG$raw_data_path)), "",
        sprintf("* Rows read: **%d**; columns: %d.", n_raw, ncol(raw)),
        sprintf("* Expected variables absent: %s", if (length(missing_vars)) paste(missing_vars, collapse = ", ") else "none"),
        sprintf("* Extra variables: %s", if (length(extra_vars)) paste(extra_vars, collapse = ", ") else "none"),
        sprintf("* Included (`incl_cohort`): **%d**; not included: %d.", sum(dv$included), sum(!dv$included)),
        sprintf("* Centres: %d; years: %s.", uniqueN(dv$centrum), paste(range(dv$jaar, na.rm = TRUE), collapse = "-")),
        sprintf("* Coding violations: %d distinct values (see coding_violations.csv).", nrow(code_issues)),
        sprintf("* Exact duplicate rows beyond first occurrence: %d; rows in near-duplicate groups: %d.",
                dup_summary$n[2], dup_summary$n[3]),
        "", "## Logical consistency", "",
        paste0("* ", consistency$check, ": ", consistency$n),
        "", "## Variables structurally missing in some years", "",
        if (nrow(struct_flag)) paste0("* `", struct_flag$variable, "`: >=90% missing in ",
                                      struct_flag$years_ge90pct_missing,
                                      "; analysed (non-missing records) from ", struct_flag$first_year_lt90pct, " only.")
        else "* none",
        "", "## Out-of-range values set to missing", "",
        paste0("* `", rbindlist(range_log)$variable, "`: ", rbindlist(range_log)$n_out_of_range_set_NA,
               " (unknown code -1: ", rbindlist(range_log)$n_unknown_code, ")"),
        "", "Undocumented derived variables (described only): ",
        paste0("`", intersect(UNDOCUMENTED_VARS, names(dv)), "`", collapse = ", "))
writeLines(md, file.path(out_dir(SUB), "validation_report.md"))
log_msg("Validation complete.")

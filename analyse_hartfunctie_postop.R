###############################################################################
#  Ventrikelfunctie (LV / RV) vóór en na hartchirurgie
#  ---------------------------------------------------------------------------
#  1. Cohort geopereerd NA juni 2022: post-op TTE vs pre-op TTE (paired) en
#     absolute post-op waarden, gecorrigeerd voor tijd tot post-op TTE.
#  2. Idem gecorrigeerd voor leeftijd, geslacht, EuroSCORE II (+ tijd tot TTE).
#  3. Vergelijking cohort na juni 2022 vs vóór juni 2022.
#  4. Aortaklemtijd (AoX): effect op kamerfunctie + mediatie-analyses.
#  5. Aanvullende analyses (RV-dysfunctie, POAF, biomarkers, uitkomsten,
#     interrupted time series, multiple imputation).
#
#  Alle output komt in de map `output/` (CSV-tabellen + PDF-figuren).
#  Pas alleen het blok INSTELLINGEN aan.
###############################################################################

## ======================= INSTELLINGEN ======================================
data_bestand <- "Aortaboog 2014-2025 (1)_260412_AL (1).xlsx"   # .xlsx of .csv (bestaat het niet: kies-venster)
excel_sheet  <- 1
cutoff_datum <- as.Date("2022-07-01")   # "na juni 2022" = vanaf 1 juli 2022
ref_dag_tte  <- 5                       # rapporteer gecorrigeerde verandering bij TTE op dag X post-op
alleen_eerste_ok_per_patient <- TRUE    # bij meerdere OK's per pat_nr alleen de eerste
complicatie_NA_is_nee <- TRUE           # lege complicatie-cel = geen complicatie
n_boot_mediatie <- 2000                 # bootstrap-simulaties mediatie
run_mice <- TRUE                        # multiple imputation sensitiviteitsanalyse
out_dir  <- "output"
set.seed(2022)
## ===========================================================================

pkgs <- c("dplyr", "tidyr", "ggplot2", "readxl", "broom", "emmeans", "lme4",
          "mediation", "MASS", "tableone", "mice", "WeightIt", "cobalt", "sandwich")
nieuw <- setdiff(pkgs, rownames(installed.packages()))
if (length(nieuw)) try(install.packages(nieuw))
heeft <- function(p) requireNamespace(p, quietly = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(emmeans); library(lme4); library(splines)
})
select <- dplyr::select   # MASS maskeert anders dplyr::select
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pad <- normalizePath(out_dir)   # volledig pad, ook na source(chdir = TRUE)
schrijf <- function(x, naam) write.csv(x, file.path(out_dir, paste0(naam, ".csv")), row.names = FALSE)

###############################################################################
# 1. INLEZEN EN OPSCHONEN
###############################################################################
if (!file.exists(data_bestand)) data_bestand <- file.choose()
if (grepl("\\.csv$", data_bestand, ignore.case = TRUE)) {
  raw <- read.csv(data_bestand, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NULL", "-", "."))
} else {
  raw <- readxl::read_excel(data_bestand, sheet = excel_sheet, .name_repair = "minimal",
                            na = c("", "NA", "NULL", "-", "."))
}
raw <- as.data.frame(raw)

# 'AF' komt twee keer voor: eerste = ritme bij pre-op TTE, tweede = bij post-op TTE
af_idx <- which(trimws(names(raw)) == "AF")
if (length(af_idx) >= 1) names(raw)[af_idx[1]] <- "pre_af"
if (length(af_idx) >= 2) names(raw)[af_idx[2]] <- "post_af"

# Robuuste kolomnamen: kleine letters, alleen a-z0-9 en '_'
schoon_naam <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT", sub = "")
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(grepl("^[0-9]", x), paste0("x", x), x)
}
names(raw) <- make.unique(schoon_naam(names(raw)), sep = "_")

# AVG: direct identificerende gegevens verwijderen (niet nodig voor analyse)
raw <- raw %>% select(-any_of(c("geb_naam", "partnernaam", "postcode", "geb_datum")))

# Helpers ---------------------------------------------------------------------
parse_datum <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))   # Excel serienummer
  x <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(x)))
  for (f in c("%d-%m-%Y", "%Y-%m-%d", "%d/%m/%Y", "%d-%m-%y", "%d.%m.%Y")) {
    nieuw <- as.Date(x, format = f)
    out[is.na(out)] <- nieuw[is.na(out)]
  }
  out
}
parse_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", gsub("[<>\\s]", "", as.character(x), perl = TRUE))))
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x > 0)))
  s <- tolower(trimws(as.character(x)))
  ifelse(s %in% c("1", "ja", "j", "yes", "y", "true", "waar", "x"), 1L,
         ifelse(s %in% c("0", "nee", "n", "no", "false", "onwaar"), 0L, NA_integer_))
}

# Datums
datum_cols <- grep("datum", names(raw), value = TRUE)
raw[datum_cols] <- lapply(raw[datum_cols], parse_datum)

# Tekstkolommen die eigenlijk numeriek zijn (bv. decimale komma) omzetten
for (v in names(raw)) {
  if (is.character(raw[[v]])) {
    x <- raw[[v]]; nn <- !is.na(x) & x != ""
    if (sum(nn) > 0 && mean(!is.na(parse_num(x[nn]))) > 0.9) raw[[v]] <- parse_num(x)
  }
}

# Binaire variabelen ------------------------------------------------------------
comorb_vars <- c("diabetes", "chr_longziekte", "cardiochir_eerder", "endocarditis",
                 "krit_preop_toestand", "recent_mi", "thorac_aortachir", "dialyse",
                 "cva_eerder", "pci_eerder", "ecc", "circ_arrest", "cabg", "klepchir",
                 "aklepchir", "mklepchir", "pklepchir", "tklepchir", "aortachir",
                 "interv_afgebroken", "post_ic", "post_intubatie", "pre_af", "post_af",
                 "overleden_zkh", "x30d_mortaliteit", "nieuwe_interv", "ecmo_opname")
compl_vars <- c("postop_mi", "post_op_cardiale_dysfunctie", "arm_beenwond", "longinfectie",
                "urineweginfectie", "mediastinitis", "positieve_bloedkweek", "sepsis",
                "resp_insuff", "heropname_ic", "cva_restletsel", "cva_zonder_restletsel",
                "aki_zonder_dialyse", "nierfalen", "gastroint_compl", "ritmeprobleem",
                "delier", "longembolie", "dvt", "pneumothorax", "bloeding",
                "rethoracotomie", "tamponade")
comorb_vars <- intersect(comorb_vars, names(raw))
compl_vars  <- intersect(compl_vars, names(raw))
raw[comorb_vars] <- lapply(raw[comorb_vars], to01)
raw[compl_vars]  <- lapply(raw[compl_vars], function(x) {
  y <- to01(x); if (complicatie_NA_is_nee) y[is.na(y)] <- 0L; y })

# Echo-parameters: korte, consistente namen -------------------------------------
echo_map <- c(lvef_mod = "lvef_mod_biplane", lvef_teich = "lvef_teich",
              ivsd = "ivsd", lvesd = "lvesd", lvedd = "lvedd", lvpwd = "lvpwd",
              gls = "lv_gls_endo_peak_avg", rv_s = "s_snelheid_rv", tapse = "tapse",
              mv_ea = "mv_e_a", e_eprime = "e_e_gemiddelde", la_vol = "la_volume_biplane",
              la_strain = "la_ed_reservoir_ls_as", lvot_d = "lvot_diameter",
              lvot_vti = "lvot_v1_vti")
for (k in names(echo_map)) for (t in c("pre_", "post_")) {
  oud <- paste0(t, echo_map[[k]])
  if (oud %in% names(raw)) names(raw)[names(raw) == oud] <- paste0(t, k)
  else warning("Echo-kolom niet gevonden: ", oud)
}
# biomarkers: pre_hstnt, post_hstnt, pre_ck_mb, post_lactaat, ... (na schoon_naam)

###############################################################################
# 2. AFGELEIDE VARIABELEN
###############################################################################
med <- function(x) median(x, na.rm = TRUE)

d <- raw %>%
  filter(is.na(interv_afgebroken) | interv_afgebroken == 0) %>%
  arrange(pat_nr, interv_datum)
if (alleen_eerste_ok_per_patient) d <- d %>% group_by(pat_nr) %>% slice(1) %>% ungroup()

# Eenheden harmoniseren (echo: dimensies in cm, TAPSE in mm; lengte in cm)
for (t in c("pre_", "post_")) {
  for (k in c("ivsd", "lvesd", "lvedd", "lvpwd")) {
    v <- paste0(t, k); if (v %in% names(d) && isTRUE(med(d[[v]]) > 15)) d[[v]] <- d[[v]] / 10
  }
  v <- paste0(t, "lvot_d"); if (v %in% names(d) && isTRUE(med(d[[v]]) > 5)) d[[v]] <- d[[v]] / 10
  v <- paste0(t, "tapse");  if (v %in% names(d) && isTRUE(med(d[[v]]) < 5))  d[[v]] <- d[[v]] * 10
  v <- paste0(t, "gls");    if (v %in% names(d)) d[[v]] <- abs(d[[v]])   # |GLS|: hoger = beter
}
if (isTRUE(med(d$lengte) < 3)) d$lengte <- d$lengte * 100

d <- d %>%
  mutate(
    periode   = factor(ifelse(interv_datum >= cutoff_datum, "na_jun2022", "voor_jun2022"),
                       levels = c("voor_jun2022", "na_jun2022")),
    geslacht  = factor(geslacht),
    euroscore2 = euroscoreii,
    bsa       = sqrt(lengte * gewicht / 3600),                       # Mosteller, m2
    # tijd tot post-op TTE (dagen); herberekenen indien leeg
    post_dagen = coalesce(post_dagen, as.numeric(post_datum - interv_datum)),
    pre_dagen  = coalesce(pre_dagen,  as.numeric(interv_datum - pre_datum)),
    tdub       = log2(post_dagen + 1),        # coefficient = effect per verdubbeling van tijd
    aox30      = aox_duur / 30,               # per 30 min klemtijd
    ecc_zonder_klem = ecc_duur - aox_duur,    # reperfusie-/overige pomptijd
    log_post_hstnt  = log(post_hstnt),
    # LVEF: biplane (MOD) bij voorkeur, anders Teichholz
    pre_lvef  = coalesce(pre_lvef_mod,  pre_lvef_teich),
    post_lvef = coalesce(post_lvef_mod, post_lvef_teich),
    # Slagvolume (ml) = LVOT-oppervlak (cm2) x VTI (cm); geïndexeerd voor BSA
    pre_sv    = pi * (pre_lvot_d / 2)^2 * pre_lvot_vti,
    post_sv   = pi * (post_lvot_d / 2)^2 * post_lvot_vti,
    pre_svi   = pre_sv / bsa, post_svi = post_sv / bsa,
    # LV-massa (ASE, lineair, cm) en geïndexeerd
    pre_lvmi  = (0.8 * 1.04 * ((pre_ivsd + pre_lvedd + pre_lvpwd)^3 - pre_lvedd^3) + 0.6) / bsa,
    post_lvmi = (0.8 * 1.04 * ((post_ivsd + post_lvedd + post_lvpwd)^3 - post_lvedd^3) + 0.6) / bsa,
    # Afkapwaarden (ASE/EACVI 2015)
    pre_lv_dysf  = as.integer(pre_lvef < 50),  post_lv_dysf  = as.integer(post_lvef < 50),
    pre_gls_abn  = as.integer(pre_gls < 18),   post_gls_abn  = as.integer(post_gls < 18),
    pre_rv_dysf  = as.integer(pre_tapse < 17 | pre_rv_s < 9.5),
    post_rv_dysf = as.integer(post_tapse < 17 | post_rv_s < 9.5),
    # Type ingreep
    ingreep = case_when(
      aortachir == 1                                   ~ "Aorta(+/-)",
      cabg == 1 & klepchir == 1                        ~ "CABG+klep",
      cabg == 1 & (is.na(klepchir) | klepchir == 0)    ~ "Geisoleerde CABG",
      klepchir == 1                                    ~ "Klep",
      TRUE                                             ~ "Overig"),
    ingreep = relevel(factor(ingreep), ref = "Geisoleerde CABG"),
    off_pump = as.integer(ecc == 0)
  )

# Composiet complicaties
d$compl_any <- as.integer(rowSums(d[compl_vars], na.rm = TRUE) > 0)
major <- intersect(c("x30d_mortaliteit", "postop_mi", "cva_restletsel", "nierfalen",
                     "rethoracotomie", "post_op_cardiale_dysfunctie", "resp_insuff",
                     "sepsis", "mediastinitis"), names(d))
d$compl_major <- as.integer(rowSums(d[major], na.rm = TRUE) > 0)

# Parameters voor de analyses
lv_par   <- c("lvef", "lvef_mod", "lvef_teich", "gls", "lvedd", "lvesd", "svi", "lvot_vti")
rv_par   <- c("tapse", "rv_s")
dia_par  <- c("e_eprime", "mv_ea", "la_vol", "la_strain")
alle_par <- c(lv_par, rv_par, dia_par, "lvmi")
hoofd_par <- c("lvef", "gls", "tapse", "rv_s")   # primaire LV/RV functie-uitkomsten

for (p in alle_par) d[[paste0("delta_", p)]] <- d[[paste0("post_", p)]] - d[[paste0("pre_", p)]]

na  <- d %>% filter(periode == "na_jun2022")
voor <- d %>% filter(periode == "voor_jun2022")
cat("N totaal:", nrow(d), "| na juni 2022:", nrow(na), "| voor juni 2022:", nrow(voor), "\n")

###############################################################################
# 3. TABEL 1 (baseline per periode)
###############################################################################
t1_vars <- intersect(c("leeftijd_operatie", "geslacht", "bmi", "bsa", "kreatinine_gehalte",
                       "diabetes", "chr_longziekte", "dialyse", "cva_eerder", "pci_eerder",
                       "cardiochir_eerder", "recent_mi", "krit_preop_toestand", "endocarditis",
                       "euroscore2", "ingreep", "off_pump", "ecc_duur", "aox_duur", "min_temp",
                       "rbc", "pre_lvef", "pre_gls", "pre_tapse", "pre_rv_s", "pre_e_eprime",
                       "pre_la_strain", "pre_af", "pre_dagen", "post_dagen", "post_ic",
                       "post_intubatie"), names(d))
t1_cat <- intersect(c("geslacht", "diabetes", "chr_longziekte", "dialyse", "cva_eerder",
                      "pci_eerder", "cardiochir_eerder", "recent_mi", "krit_preop_toestand",
                      "endocarditis", "ingreep", "off_pump", "pre_af", "post_ic",
                      "post_intubatie"), t1_vars)
t1_nonnormal <- intersect(c("euroscore2", "kreatinine_gehalte", "ecc_duur", "aox_duur",
                            "rbc", "pre_dagen", "post_dagen"), t1_vars)
if (heeft("tableone")) {
tab1 <- tableone::CreateTableOne(vars = t1_vars, strata = "periode", data = d, factorVars = t1_cat,
                       addOverall = TRUE)
tab1_df <- print(tab1, nonnormal = t1_nonnormal, smd = TRUE, missing = TRUE,
                 showAllLevels = TRUE, printToggle = FALSE)
write.csv(tab1_df, file.path(out_dir, "tabel1_baseline_per_periode.csv"))
}

###############################################################################
# 4. COHORT NA JUNI 2022: PRE vs POST
###############################################################################
## 4a. Ongecorrigeerd, gepaard --------------------------------------------------
paired_tab <- function(dat, pars) {
  bind_rows(lapply(pars, function(p) {
    pre <- dat[[paste0("pre_", p)]]; post <- dat[[paste0("post_", p)]]
    ok <- complete.cases(pre, post)
    if (sum(ok) < 5) return(NULL)
    tt <- t.test(post[ok], pre[ok], paired = TRUE)
    wt <- suppressWarnings(wilcox.test(post[ok], pre[ok], paired = TRUE, exact = FALSE))
    dl <- post[ok] - pre[ok]
    tibble(parameter = p, n_paar = sum(ok),
           pre_mean = mean(pre[ok]), pre_sd = sd(pre[ok]),
           post_mean = mean(post[ok]), post_sd = sd(post[ok]),
           post_mean_alle = mean(post, na.rm = TRUE), n_post_alle = sum(!is.na(post)),
           delta_mean = mean(dl), ci_laag = tt$conf.int[1], ci_hoog = tt$conf.int[2],
           delta_pct = 100 * mean(dl) / mean(pre[ok]),
           delta_mediaan = median(dl), delta_q1 = quantile(dl, .25), delta_q3 = quantile(dl, .75),
           p_ttest = tt$p.value, p_wilcoxon = wt$p.value)
  }))
}
tab_paired <- paired_tab(na, alle_par)
tab_paired$p_holm <- p.adjust(tab_paired$p_ttest, "holm")
schrijf(tab_paired, "4a_paired_pre_post_na2022")

## 4b. Absoluut: % patiënten met afwijkende functie pre vs post (McNemar) --------
mcnemar_tab <- function(dat, vars) {
  bind_rows(lapply(vars, function(v) {
    pre <- dat[[paste0("pre_", v)]]; post <- dat[[paste0("post_", v)]]
    ok <- complete.cases(pre, post); if (sum(ok) < 5) return(NULL)
    tb <- table(factor(pre[ok], 0:1), factor(post[ok], 0:1))
    tibble(maat = v, n = sum(ok), pre_pct = 100 * mean(pre[ok]), post_pct = 100 * mean(post[ok]),
           nieuw_afwijkend = tb["0", "1"], hersteld = tb["1", "0"],
           p_mcnemar = tryCatch(mcnemar.test(tb)$p.value, error = function(e) NA))
  }))
}
schrijf(mcnemar_tab(na, c("lv_dysf", "gls_abn", "rv_dysf")), "4b_afwijkend_pre_post_na2022")

## 4c. Gecorrigeerd voor tijd tot post-op TTE (ANCOVA op verandering) ----------
#  delta = b0 + b1*(pre - gemiddelde pre) + b2*log2(dagen+1) [+ covariaten]
#  -> emmeans geeft gecorrigeerde gemiddelde verandering bij TTE op dag `ref_dag_tte`
#  -> b2 = verandering in herstel per verdubbeling van tijd tot TTE
covsets <- list(
  A_tijd             = character(0),
  B_tijd_lft_gesl_es2 = c("leeftijd_operatie", "geslacht", "euroscore2"),
  C_B_plus_echo_omst = c("leeftijd_operatie", "geslacht", "euroscore2", "post_ic", "post_intubatie"),
  D_lft_gesl_es2_zonder_tijd = c("leeftijd_operatie", "geslacht", "euroscore2")
)
covsets <- lapply(covsets, intersect, names(d))

fit_change <- function(dat, p, covs, met_tijd = TRUE, extra = NULL) {
  pre <- paste0("pre_", p); dl <- paste0("delta_", p)
  rhs <- c(pre, if (met_tijd) "tdub", covs, extra)
  dd <- dat[complete.cases(dat[, c(dl, rhs)]), ]
  if (nrow(dd) < length(rhs) + 10) return(NULL)
  f <- reformulate(rhs, response = dl)
  m <- lm(f, data = dd)
  at <- setNames(list(mean(dd[[pre]])), pre)
  if (met_tijd) at$tdub <- log2(ref_dag_tte + 1)
  em <- summary(emmeans(m, ~ 1, at = at))
  tijd <- if (met_tijd) tidy(m, conf.int = TRUE) %>% filter(term == "tdub") else NULL
  list(model = m, n = nrow(dd), em = em, tijd = tijd)
}
res_4c <- list()
for (p in c(lv_par, rv_par, dia_par)) for (cs in names(covsets)) {
  r <- fit_change(na, p, covsets[[cs]], met_tijd = !grepl("zonder_tijd", cs))
  if (is.null(r)) next
  res_4c[[length(res_4c) + 1]] <- tibble(
    parameter = p, model = cs, n = r$n,
    adj_delta = r$em$emmean, ci_laag = r$em$lower.CL, ci_hoog = r$em$upper.CL,
    p_delta = 2 * pt(-abs(r$em$emmean / r$em$SE), r$em$df),
    tijd_per_verdubbeling = if (!is.null(r$tijd)) r$tijd$estimate else NA,
    tijd_ci_laag = if (!is.null(r$tijd)) r$tijd$conf.low else NA,
    tijd_ci_hoog = if (!is.null(r$tijd)) r$tijd$conf.high else NA,
    tijd_p = if (!is.null(r$tijd)) r$tijd$p.value else NA)
}
tab_4c <- bind_rows(res_4c)
schrijf(tab_4c, "4c_gecorrigeerde_verandering_na2022")

## 4d. Absolute post-op waarden gecorrigeerd (zonder pre-waarde in model) -------
res_4d <- list()
for (p in hoofd_par) for (cs in c("A_tijd", "B_tijd_lft_gesl_es2")) {
  y <- paste0("post_", p); rhs <- c("tdub", covsets[[cs]])
  dd <- na[complete.cases(na[, c(y, rhs)]), ]; if (nrow(dd) < 15) next
  m <- lm(reformulate(rhs, y), data = dd)
  em <- summary(emmeans(m, ~ 1, at = list(tdub = log2(ref_dag_tte + 1))))
  res_4d[[length(res_4d) + 1]] <- tibble(parameter = p, model = cs, n = nrow(dd),
                                         post_adj = em$emmean, ci_laag = em$lower.CL, ci_hoog = em$upper.CL)
}
schrijf(bind_rows(res_4d), "4d_absolute_postop_waarden_gecorrigeerd")

## 4e. Niet-lineair herstel: natural spline voor tijd tot TTE (sensitiviteit) --
spline_tab <- bind_rows(lapply(hoofd_par, function(p) {
  dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  dd <- na[complete.cases(na[, c(dl, pre, "post_dagen")]), ]; if (nrow(dd) < 30) return(NULL)
  m_lin <- lm(reformulate(c(pre, "tdub"), dl), dd)
  m_spl <- lm(reformulate(c(pre, "ns(post_dagen, df = 3)"), dl), dd)
  tibble(parameter = p, n = nrow(dd), AIC_log2 = AIC(m_lin), AIC_spline = AIC(m_spl),
         p_nonlineair = anova(m_lin, m_spl)$`Pr(>F)`[2])
}))
schrijf(spline_tab, "4e_spline_vs_log_tijd")

## 4f. Mixed model (gebruikt ook patiënten met alleen pre- óf post-echo) --------
echo_kol <- intersect(paste0(rep(c("pre_", "post_"), each = length(hoofd_par)), hoofd_par), names(na))
lang <- na %>%
  select(interv_nr, leeftijd_operatie, geslacht, euroscore2, post_dagen, all_of(echo_kol)) %>%
  pivot_longer(all_of(echo_kol), names_to = c("moment", "par"),
               names_pattern = "^(pre|post)_(.*)$", values_to = "waarde") %>%
  mutate(moment = factor(moment, c("pre", "post")),
         # tijd-term alleen op post-meting (pre = baseline)
         tdub_post = ifelse(moment == "post", log2(post_dagen + 1) - log2(ref_dag_tte + 1), 0)) %>%
  filter(!is.na(waarde))
mm_tab <- bind_rows(lapply(hoofd_par, function(p) {
  dd <- lang %>% filter(par == p, !is.na(euroscore2), !is.na(geslacht), !is.na(leeftijd_operatie)) %>%
    filter(moment == "pre" | !is.na(post_dagen))
  if (n_distinct(dd$interv_nr) < 20) return(NULL)
  m <- lmer(waarde ~ moment + tdub_post + leeftijd_operatie + geslacht + euroscore2 + (1 | interv_nr),
            data = dd)
  ci <- confint(m, parm = "momentpost", method = "Wald")
  tibble(parameter = p, n_pat = n_distinct(dd$interv_nr), n_obs = nrow(dd),
         post_vs_pre_bij_refdag = fixef(m)["momentpost"], ci_laag = ci[1], ci_hoog = ci[2])
}))
schrijf(mm_tab, "4f_mixed_model_post_vs_pre")

## 4g. Figuren ------------------------------------------------------------------
pdf(file.path(out_dir, "fig_4_pre_post_na2022.pdf"), width = 8, height = 6)
for (p in hoofd_par) {
  pre <- paste0("pre_", p); post <- paste0("post_", p)
  dd <- na %>% filter(!is.na(.data[[pre]]), !is.na(.data[[post]])) %>%
    select(interv_nr, post_dagen, all_of(c(pre, post)))
  if (nrow(dd) < 5) next
  # spaghetti pre -> post
  pl <- dd %>% pivot_longer(c(all_of(pre), all_of(post)), names_to = "moment", values_to = "w") %>%
    mutate(moment = factor(sub("_.*", "", moment), c("pre", "post"))) %>%
    ggplot(aes(moment, w, group = interv_nr)) + geom_line(alpha = .25) +
    geom_boxplot(aes(group = moment), width = .3, outlier.shape = NA, fill = NA, colour = "steelblue") +
    labs(title = paste(p, ": pre-op vs post-op (na juni 2022)"), x = NULL, y = p) + theme_minimal()
  print(pl)
  # verandering vs tijd tot TTE (herstelcurve)
  pl2 <- ggplot(dd, aes(post_dagen, .data[[post]] - .data[[pre]])) +
    geom_hline(yintercept = 0, linetype = 2) + geom_point(alpha = .5) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 3, 5, 7, 14, 30, 60, 90)) +
    labs(title = paste("Verandering", p, "vs tijd tot post-op TTE"),
         x = "Dagen na OK tot TTE (log-schaal)", y = paste("Delta", p)) + theme_minimal()
  print(pl2)
}
dev.off()

###############################################################################
# 5. NA JUNI 2022 vs VOOR JUNI 2022
###############################################################################
adj_basis <- intersect(c("leeftijd_operatie", "geslacht", "euroscore2", "ingreep"), names(d))

## 5a. Klinische uitkomsten (beschikbaar in beide periodes) ---------------------
uitk_bin <- intersect(c("x30d_mortaliteit", "overleden_zkh", "compl_any", "compl_major",
                        "ritmeprobleem", "aki_zonder_dialyse", "nierfalen", "delier",
                        "rethoracotomie", "post_op_cardiale_dysfunctie", "heropname_ic",
                        "ecmo_opname"), names(d))
uitk_cont <- intersect(c("icu_los", "in_hospital_los"), names(d))

res_5a <- list()
for (y in uitk_bin) {
  if (sum(d[[y]], na.rm = TRUE) < 5) next
  for (adj in list(character(0), adj_basis)) {
    dd <- d[complete.cases(d[, c(y, "periode", adj)]), ]
    m <- glm(reformulate(c("periode", adj), y), family = binomial, data = dd)
    tt <- tidy(m, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald") %>%
      filter(term == "periodena_jun2022")
    res_5a[[length(res_5a) + 1]] <- tibble(uitkomst = y, type = "OR",
      model = ifelse(length(adj), "gecorrigeerd", "ruw"), n = nrow(dd),
      pct_voor = 100 * mean(dd[[y]][dd$periode == "voor_jun2022"]),
      pct_na = 100 * mean(dd[[y]][dd$periode == "na_jun2022"]),
      schatting = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high, p = tt$p.value)
  }
}
# LOS: scheef verdeeld -> log-lineair; exp(coef) = ratio van geometrisch gemiddelde
for (y in uitk_cont) for (adj in list(character(0), adj_basis)) {
  dd <- d[complete.cases(d[, c(y, "periode", adj)]) & d[[y]] >= 0, ]
  m <- lm(reformulate(c("periode", adj), paste0("log1p(", y, ")")), data = dd)
  tt <- tidy(m, conf.int = TRUE) %>% filter(term == "periodena_jun2022")
  res_5a[[length(res_5a) + 1]] <- tibble(uitkomst = y, type = "ratio (geom. gem.)",
    model = ifelse(length(adj), "gecorrigeerd", "ruw"), n = nrow(dd),
    pct_voor = median(dd[[y]][dd$periode == "voor_jun2022"]),
    pct_na = median(dd[[y]][dd$periode == "na_jun2022"]),
    schatting = exp(tt$estimate), ci_laag = exp(tt$conf.low), ci_hoog = exp(tt$conf.high), p = tt$p.value)
}
schrijf(bind_rows(res_5a), "5a_uitkomsten_na_vs_voor_jun2022")   # pct_voor/na = mediaan bij LOS

## 5b. Sensitiviteit: inverse probability of treatment weighting (IPTW) ----------
if (heeft("WeightIt") && heeft("sandwich")) {
  ps_vars <- intersect(c("leeftijd_operatie", "geslacht", "euroscore2", "ingreep", "bmi",
                         "diabetes", "chr_longziekte", "cardiochir_eerder"), names(d))
  dw <- d[complete.cases(d[, ps_vars]), ]
  w <- WeightIt::weightit(reformulate(ps_vars, "periode"), data = dw, method = "glm",
                          estimand = "ATE", stabilize = TRUE)
  dw$w <- w$weights
  if (requireNamespace("cobalt", quietly = TRUE)) {
    bal <- cobalt::bal.tab(w, un = TRUE, thresholds = c(m = .1))
    write.csv(as.data.frame(bal$Balance), file.path(out_dir, "5b_balans_IPTW.csv"))
  }
  res_5b <- bind_rows(lapply(uitk_bin, function(y) {
    dd <- dw[!is.na(dw[[y]]), ]; if (sum(dd[[y]]) < 5) return(NULL)
    m <- glm(reformulate("periode", y), family = quasibinomial, data = dd, weights = w)
    V <- sandwich::vcovHC(m, type = "HC0"); b <- coef(m)["periodena_jun2022"]
    se <- sqrt(V["periodena_jun2022", "periodena_jun2022"])
    tibble(uitkomst = y, OR_iptw = exp(b), ci_laag = exp(b - 1.96 * se), ci_hoog = exp(b + 1.96 * se),
           p = 2 * pnorm(-abs(b / se)))
  }))
  schrijf(res_5b, "5b_uitkomsten_IPTW")
}

## 5c. Echo-vergelijking tussen periodes (alleen als er echo's vóór juni 2022 zijn)
res_5c <- list()
for (p in hoofd_par) {
  dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  if (sum(!is.na(voor[[dl]])) < 10) { message("Te weinig echo-data voor juni 2022 voor: ", p); next }
  for (cs in c("A_tijd", "B_tijd_lft_gesl_es2")) {
    rhs <- c("periode", pre, "tdub", covsets[[cs]])
    dd <- d[complete.cases(d[, c(dl, rhs)]), ]
    m <- lm(reformulate(rhs, dl), dd)
    tt <- tidy(m, conf.int = TRUE) %>% filter(term == "periodena_jun2022")
    em <- summary(emmeans(m, ~ periode, at = list(tdub = log2(ref_dag_tte + 1))))
    res_5c[[length(res_5c) + 1]] <- tibble(parameter = p, model = cs, n = nrow(dd),
      adj_delta_voor = em$emmean[1], adj_delta_na = em$emmean[2],
      verschil_na_min_voor = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high, p = tt$p.value)
  }
}
if (length(res_5c)) schrijf(bind_rows(res_5c), "5c_echo_verandering_na_vs_voor")

## 5d. Interrupted time series (trend vs. sprong rond juni 2022) -----------------
#  Onderscheidt een echte verandering per juni 2022 van een al bestaande trend.
d <- d %>% mutate(maand = as.numeric(interv_datum - min(interv_datum, na.rm = TRUE)) / 30.44,
                  maand_na = pmax(0, as.numeric(interv_datum - cutoff_datum) / 30.44))
its_tab <- bind_rows(lapply(intersect(c("compl_major", "compl_any", "x30d_mortaliteit"), uitk_bin), function(y) {
  dd <- d[complete.cases(d[, c(y, "maand", adj_basis)]), ]
  if (sum(dd[[y]]) < 10) return(NULL)
  m <- glm(reformulate(c("maand", "periode", "maand_na", adj_basis), y), binomial, dd)
  tidy(m, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald") %>%
    filter(term %in% c("maand", "periodena_jun2022", "maand_na")) %>% mutate(uitkomst = y, .before = 1)
}))
schrijf(its_tab, "5d_interrupted_time_series")

###############################################################################
# 6. AORTAKLEMTIJD (AoX) EN KAMERFUNCTIE
###############################################################################
aox_adj <- c("leeftijd_operatie", "geslacht", "euroscore2", "tdub")
on_pump <- na %>% filter(!is.na(aox_duur), aox_duur > 0)

## 6a. Totaal effect van AoX (per 30 min) op verandering in functie ---------------
res_6a <- list()
for (p in hoofd_par) {
  pre <- paste0("pre_", p); dl <- paste0("delta_", p)
  for (spec in list(c("aox30"),                               # alleen AoX
                    c("aox30", "ingreep"),                    # + type ingreep
                    c("aox30", "ecc_zonder_klem"),            # klemtijd los van overige pomptijd
                    c("aox30", "min_temp"))) {
    spec <- intersect(spec, names(on_pump))
    rhs <- c(spec, pre, aox_adj)
    dd <- on_pump[complete.cases(on_pump[, c(dl, rhs)]), ]; if (nrow(dd) < 25) next
    m <- lm(reformulate(rhs, dl), dd)
    tt <- tidy(m, conf.int = TRUE) %>% filter(term == "aox30")
    res_6a[[length(res_6a) + 1]] <- tibble(parameter = p, extra_correctie = paste(setdiff(spec, "aox30"), collapse = "+"),
      n = nrow(dd), effect_per_30min_aox = tt$estimate, ci_laag = tt$conf.low,
      ci_hoog = tt$conf.high, p = tt$p.value)
  }
  # effectmodificatie: zijn slechte ventrikels kwetsbaarder voor lange klemtijd?
  rhs <- c("aox30", pre, aox_adj)
  dd <- on_pump[complete.cases(on_pump[, c(dl, rhs)]), ]
  if (nrow(dd) >= 25) {
    mi <- lm(as.formula(paste(dl, "~ aox30 *", pre, "+", paste(aox_adj, collapse = "+"))), dd)
    ti <- tidy(mi) %>% filter(grepl(":", term))
    res_6a[[length(res_6a) + 1]] <- tibble(parameter = p, extra_correctie = paste0("interactie aox30 x ", pre),
      n = nrow(dd), effect_per_30min_aox = ti$estimate, ci_laag = NA, ci_hoog = NA, p = ti$p.value)
  }
}
schrijf(bind_rows(res_6a), "6a_AoX_effect_op_kamerfunctie")

# Dosis-respons figuur
pdf(file.path(out_dir, "fig_6_AoX_dosis_respons.pdf"), width = 7, height = 5)
for (p in hoofd_par) {
  dl <- paste0("delta_", p)
  dd <- on_pump %>% filter(!is.na(.data[[dl]]))
  if (nrow(dd) < 10) next
  print(ggplot(dd, aes(aox_duur, .data[[dl]])) + geom_hline(yintercept = 0, linetype = 2) +
          geom_point(alpha = .5) + geom_smooth(method = "lm", formula = y ~ ns(x, 3)) +
          labs(x = "Aortaklemtijd (min)", y = paste("Delta", p), title = paste("AoX vs delta", p)) +
          theme_minimal())
}
dev.off()

## 6b. Mediatie-analyses ---------------------------------------------------------
#  M1: AoX -> myocardschade (log post-op hsTnT) -> delta kamerfunctie
#  M2: AoX -> post-op kamerfunctie -> klinische uitkomst (IC-duur, complicaties)
#  ACME = indirect (gemedieerd) effect, ADE = direct effect, per 30 min AoX.
#  Let op: causale interpretatie vereist o.a. geen ongemeten mediator-uitkomst confounding;
#  rapporteer ook de sensitiviteitsanalyse (medsens) voor continue modellen.
med_run <- function(dat, mediator, uitkomst, covs, label, y_binair = FALSE) {
  vars <- unique(c("aox30", mediator, uitkomst, covs))
  dd <- dat[complete.cases(dat[, vars]), ]
  if (nrow(dd) < 40) { message("Mediatie ", label, ": te weinig data (n=", nrow(dd), ")"); return(NULL) }
  m_m <- lm(reformulate(c("aox30", covs), mediator), dd)
  m_y <- if (y_binair) glm(reformulate(c("aox30", mediator, covs), uitkomst), binomial(link = "probit"), dd)
         else lm(reformulate(c("aox30", mediator, covs), uitkomst), dd)
  md <- mediation::mediate(m_m, m_y, treat = "aox30", mediator = mediator, boot = TRUE, sims = n_boot_mediatie)
  s <- summary(md)
  sens <- if (!y_binair) tryCatch(mediation::medsens(md, rho.by = 0.05, sims = 200)$err.cr.d[1], error = function(e) NA) else NA
  tibble(analyse = label, n = nrow(dd),
         ACME = s$d.avg, ACME_ci_l = s$d.avg.ci[1], ACME_ci_h = s$d.avg.ci[2], ACME_p = s$d.avg.p,
         ADE = s$z.avg, ADE_ci_l = s$z.avg.ci[1], ADE_ci_h = s$z.avg.ci[2], ADE_p = s$z.avg.p,
         totaal = s$tau.coef, totaal_p = s$tau.p,
         prop_gemedieerd = s$n.avg, prop_ci_l = s$n.avg.ci[1], prop_ci_h = s$n.avg.ci[2],
         rho_waarbij_ACME_0 = sens)
}
on_pump$log_icu_los <- log1p(on_pump$icu_los)
on_pump$log_hosp_los <- log1p(on_pump$in_hospital_los)
cov_basis <- c("leeftijd_operatie", "geslacht", "euroscore2")
res_6b <- list()
if (heeft("mediation")) for (p in hoofd_par) {
  pre <- paste0("pre_", p)
  res_6b[[length(res_6b) + 1]] <- med_run(on_pump, "log_post_hstnt", paste0("delta_", p),
                                          c(cov_basis, pre, "tdub"), paste0("AoX -> hsTnT -> delta_", p))
  for (y in c("log_icu_los", "log_hosp_los"))
    res_6b[[length(res_6b) + 1]] <- med_run(on_pump, paste0("post_", p), y,
                                            c(cov_basis, pre, "tdub"), paste0("AoX -> post_", p, " -> ", y))
  res_6b[[length(res_6b) + 1]] <- med_run(on_pump, paste0("post_", p), "compl_major",
                                          c(cov_basis, pre, "tdub"), paste0("AoX -> post_", p, " -> compl_major"),
                                          y_binair = TRUE)
}
schrijf(bind_rows(res_6b), "6b_mediatie_AoX")

###############################################################################
# 7. AANVULLENDE ANALYSES
###############################################################################
## 7a. Determinanten van nieuwe post-op RV-dysfunctie en LV-dysfunctie ----------
det_vars <- intersect(c("leeftijd_operatie", "geslacht", "euroscore2", "ingreep", "aox30",
                        "min_temp", "tdub", "post_intubatie", "tklepchir", "mklepchir"), names(na))
res_7a <- bind_rows(lapply(c("post_rv_dysf", "post_lv_dysf"), function(y) {
  pre <- sub("post_", "pre_", y)
  dd <- na %>% filter(.data[[pre]] == 0)               # alleen patiënten met normale pre-op functie
  dd <- dd[complete.cases(dd[, c(y, det_vars)]), ]
  if (sum(dd[[y]]) < 10) return(NULL)
  m <- glm(reformulate(det_vars, y), binomial, dd)
  tidy(m, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald") %>%
    filter(term != "(Intercept)") %>% mutate(uitkomst = paste0("nieuwe ", y), n = nrow(dd), events = sum(dd[[y]]), .before = 1)
}))
schrijf(res_7a, "7a_determinanten_nieuwe_dysfunctie")

## 7b. Post-op AF / ritmeprobleem: voorspellende waarde van LA-strain/volume -----
if ("ritmeprobleem" %in% names(na)) {
  res_7b <- bind_rows(lapply(c("pre_la_strain", "pre_la_vol", "pre_e_eprime"), function(x) {
    vars <- c("ritmeprobleem", x, "leeftijd_operatie", "geslacht", "euroscore2", "pre_af")
    dd <- na[complete.cases(na[, intersect(vars, names(na))]), ]
    if (nrow(dd) < 30 || sum(dd$ritmeprobleem) < 10) return(NULL)
    dd[[x]] <- as.numeric(scale(dd[[x]]))               # per SD
    m <- glm(reformulate(setdiff(intersect(vars, names(dd)), "ritmeprobleem"), "ritmeprobleem"), binomial, dd)
    tidy(m, exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald") %>% filter(term == x) %>%
      mutate(n = nrow(dd), events = sum(dd$ritmeprobleem), .before = 1)
  }))
  schrijf(res_7b, "7b_LA_parameters_en_POAF_OR_per_SD")
}

## 7c. Biomarkers (hsTnT, CK-MB, lactaat) vs post-op functie -------------------
bio <- intersect(c("post_hstnt", "post_ck_mb", "post_ck", "post_lactaat", "aox_duur", "ecc_duur",
                   "post_lvef", "post_gls", "post_tapse", "post_rv_s", "delta_lvef", "delta_gls",
                   "delta_tapse", "delta_rv_s"), names(na))
cor_mat <- cor(na[bio], use = "pairwise.complete.obs", method = "spearman")
write.csv(round(cor_mat, 2), file.path(out_dir, "7c_spearman_biomarkers_functie.csv"))
pdf(file.path(out_dir, "fig_7c_correlatie_heatmap.pdf"), width = 8, height = 7)
print(as.data.frame(as.table(cor_mat)) %>%
        ggplot(aes(Var1, Var2, fill = Freq)) + geom_tile() + geom_text(aes(label = round(Freq, 2)), size = 2.5) +
        scale_fill_gradient2(limits = c(-1, 1), low = "#2166ac", high = "#b2182b") +
        labs(x = NULL, y = NULL, fill = "rho") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

## 7d. Post-op functie als voorspeller van klinische uitkomst -------------------
res_7d <- list()
for (p in hoofd_par) {
  x <- paste0("post_", p)
  for (y in c("icu_los", "in_hospital_los")) {
    rhs <- c(x, "tdub", cov_basis); dd <- na[complete.cases(na[, c(y, rhs)]), ]; if (nrow(dd) < 30) next
    dd[[x]] <- as.numeric(scale(dd[[x]]))
    tt <- tidy(lm(reformulate(rhs, paste0("log1p(", y, ")")), dd), conf.int = TRUE) %>% filter(term == x)
    res_7d[[length(res_7d) + 1]] <- tibble(predictor = x, uitkomst = y, maat = "ratio per SD", n = nrow(dd),
      schatting = exp(tt$estimate), ci_laag = exp(tt$conf.low), ci_hoog = exp(tt$conf.high), p = tt$p.value)
  }
  for (y in intersect(c("compl_major", "aki_zonder_dialyse", "ritmeprobleem", "delier"), names(na))) {
    rhs <- c(x, "tdub", cov_basis); dd <- na[complete.cases(na[, c(y, rhs)]), ]
    if (nrow(dd) < 30 || sum(dd[[y]]) < 10) next
    dd[[x]] <- as.numeric(scale(dd[[x]]))
    tt <- tidy(glm(reformulate(rhs, y), binomial, dd), exponentiate = TRUE, conf.int = TRUE, conf.method = "Wald") %>% filter(term == x)
    res_7d[[length(res_7d) + 1]] <- tibble(predictor = x, uitkomst = y, maat = "OR per SD", n = nrow(dd),
      schatting = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high, p = tt$p.value)
  }
}
schrijf(bind_rows(res_7d), "7d_postop_functie_en_uitkomst")

## 7e. Subgroepen: verandering per type ingreep en on- vs off-pump -------------
res_7e <- list()
for (p in hoofd_par) for (g in c("ingreep", "off_pump")) {
  dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  rhs <- c(g, pre, "tdub", cov_basis); dd <- na[complete.cases(na[, c(dl, rhs)]), ]
  dd[[g]] <- factor(dd[[g]]); if (nrow(dd) < 30 || nlevels(dd[[g]]) < 2) next
  m <- lm(reformulate(rhs, dl), dd)
  em <- as.data.frame(summary(emmeans(m, as.formula(paste("~", g)), at = list(tdub = log2(ref_dag_tte + 1)))))
  names(em)[1] <- "groep"
  res_7e[[length(res_7e) + 1]] <- em %>% mutate(parameter = p, factor = g, groep = as.character(groep),
    p_globaal = anova(m)[g, "Pr(>F)"], .before = 1)
}
schrijf(bind_rows(res_7e), "7e_subgroepen_ingreep_offpump")

## 7f. Multiple imputation (sensitiviteit voor hoofdmodel B) ---------------------
if (run_mice) {
  imp_vars <- intersect(c("leeftijd_operatie", "geslacht", "euroscore2", "tdub", "aox_duur",
                          "ingreep", "post_ic", "post_intubatie", "post_hstnt",
                          paste0("pre_", hoofd_par), paste0("post_", hoofd_par)), names(na))
  imp <- mice::mice(na[imp_vars], m = 20, maxit = 10, printFlag = FALSE, seed = 2022)
  res_7f <- bind_rows(lapply(hoofd_par, function(p) {
    f <- as.formula(paste0("I(post_", p, " - pre_", p, ") ~ I(pre_", p, " - ", mean(na[[paste0("pre_", p)]], na.rm = TRUE),
                           ") + I(tdub - ", log2(ref_dag_tte + 1), ") + I(leeftijd_operatie - ",
                           mean(na$leeftijd_operatie, na.rm = TRUE), ") + geslacht + I(euroscore2 - ",
                           mean(na$euroscore2, na.rm = TRUE), ")"))
    fit <- mice::as.mira(lapply(seq_len(imp$m), function(i) lm(f, data = mice::complete(imp, i))))
    summary(mice::pool(fit), conf.int = TRUE) %>% filter(term == "(Intercept)") %>%
      transmute(parameter = p, adj_delta_MI = estimate, ci_laag = `2.5 %`, ci_hoog = `97.5 %`, p = p.value)
  }))
  schrijf(res_7f, "7f_multiple_imputation_model_B")
  # Intercept = gecorrigeerde verandering bij gemiddelde pre-waarde/lft/ES2, TTE op ref-dag, referentiegeslacht
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

###############################################################################
# 8. RESULTATEN DIRECT TONEN (console + RStudio-tabbladen + PDF-figuren)
###############################################################################
toon_resultaten <- function(open_tabs = TRUE, open_pdf = TRUE) {
  lees <- function(f) { f <- file.path(out_pad, paste0(f, ".csv"))
    if (file.exists(f)) read.csv(f, check.names = FALSE) else NULL }
  rond <- function(x) { num <- sapply(x, is.numeric)
    x[num] <- lapply(x[num], function(v) ifelse(abs(v) < 0.001 & v != 0, signif(v, 2), round(v, 3))); x }
  tabellen <- list(
    "4a Pre vs post (na jun 2022), ongecorrigeerd" =
      function(x) x[, c("parameter", "n_paar", "pre_mean", "post_mean", "delta_mean", "ci_laag", "ci_hoog", "p_ttest", "p_holm")],
    "4b Afwijkende functie pre vs post" = identity,
    "4c Verandering gecorrigeerd (dag ref_dag_tte)" =
      function(x) x[x$model %in% c("A_tijd", "B_tijd_lft_gesl_es2"),
                    c("parameter", "model", "n", "adj_delta", "ci_laag", "ci_hoog", "p_delta", "tijd_per_verdubbeling", "tijd_p")],
    "4d Absolute post-op waarden gecorrigeerd" = identity,
    "5a Uitkomsten na vs voor jun 2022" = identity,
    "5b Uitkomsten IPTW" = identity,
    "5c Echo-verandering na vs voor jun 2022" = identity,
    "6a Effect AoX (per 30 min) op kamerfunctie" = identity,
    "6b Mediatie AoX" =
      function(x) x[, intersect(c("analyse", "n", "ACME", "ACME_ci_l", "ACME_ci_h", "ACME_p", "ADE", "ADE_p", "prop_gemedieerd"), names(x))],
    "7f Multiple imputation (model B)" = identity)
  bestanden <- c("4a_paired_pre_post_na2022", "4b_afwijkend_pre_post_na2022",
                 "4c_gecorrigeerde_verandering_na2022", "4d_absolute_postop_waarden_gecorrigeerd",
                 "5a_uitkomsten_na_vs_voor_jun2022", "5b_uitkomsten_IPTW", "5c_echo_verandering_na_vs_voor",
                 "6a_AoX_effect_op_kamerfunctie", "6b_mediatie_AoX", "7f_multiple_imputation_model_B")
  old <- options(width = 200); on.exit(options(old))
  for (i in seq_along(bestanden)) {
    x <- lees(bestanden[i]); if (is.null(x) || !nrow(x)) next
    x <- rond(tabellen[[i]](x)); titel <- names(tabellen)[i]
    cat("\n==================", titel, "==================\n"); print(x, row.names = FALSE)
    if (open_tabs && interactive() && Sys.getenv("RSTUDIO") == "1") utils::View(x, title = substr(titel, 1, 25))
  }
  if (open_pdf && interactive()) for (f in list.files(out_pad, "\\.pdf$", full.names = TRUE)) browseURL(f)
  cat("\nAlle tabellen en figuren staan in:", out_pad, "\n")
}
toon_resultaten()
cat("Klaar. Later opnieuw bekijken zonder alles te herberekenen: toon_resultaten()\n")

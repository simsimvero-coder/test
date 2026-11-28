###############################################################################
#  VERSIE 3 - MYOCARDUITKOMSTEN: kloppend hart (warme perfusie) vs cardioplegie
#  ---------------------------------------------------------------------------
#  Drie definities van de blootstelling:
#    ITT = toegewezen techniek (perf_tech_toegewezen)       -> hoofdanalyse
#    AT  = ontvangen techniek (perf_tech_ontvangen)         -> as-treated
#    PP  = alleen patiënten die de toegewezen techniek kregen -> per-protocol
#  M0. Indeling, omgezette patiënten (crossovers)
#  M1. Echo: verandering LVEF, GLS, RV S', TAPSE, LA-strain, E/e' per techniek
#  M2. Biomarkers na OK: CK, CK-MB, hsTnT, lactaat (ratio geometrische gemiddelden)
#  M3. Klemtijd x techniek (dosis-respons) voor echo en biomarkers
#  M4. Klinische myocarduitkomsten (Firth + Bayes)
#  NB: NIET gecorrigeerd voor circulatiestilstand/klemtijd in M1/M2/M4: die zijn
#      onderdeel van de strategie (mediatoren) -> correctie zou het effect wegnemen.
#  Zelfstandig: alleen dit bestand + het Excelbestand in dezelfde map.
#  Gebruik: source(file.choose(), encoding = "UTF-8", chdir = TRUE)
###############################################################################

## ======================= INSTELLINGEN v3 ===================================
code_kloppend <- 1           # waarde in perf_tech_toegewezen/-ontvangen die "kloppend hart" betekent
n_boot      <- 1000          # bootstrap voor BI's van de verandering
n_boot_med  <- 0             # 0 = snelle Monte Carlo-BI (aanbevolen); >0 = bootstrap (zeer traag)
prior_sd_zwak       <- 1.5   # zwak-informatieve prior op log-OR (95%: OR 0,05-19)
prior_sd_sceptisch  <- 0.35  # sceptische prior op log-OR (95%: OR 0,5-2)
## ===========================================================================

###############################################################################
# 0. DATA INLEZEN EN OPSCHONEN (identiek aan analyse_hartfunctie_postop.R)
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
if (getOption("repos")[["CRAN"]] == "@CRAN@") options(repos = c(CRAN = "https://cloud.r-project.org"))
# Mac/Windows: altijd kant-en-klare binaries (compileren vanaf broncode faalt zonder gfortran)
if (length(nieuw)) try(install.packages(nieuw, type = if (.Platform$pkgType == "source") "source" else "binary"))
heeft <- function(p) requireNamespace(p, quietly = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(broom)
  library(emmeans); library(lme4); library(splines)
})
select <- dplyr::select   # MASS maskeert anders dplyr::select
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pad <- normalizePath(out_dir)   # volledig pad, ook na source(chdir = TRUE)
schrijf <- function(x, naam) if (!is.null(x) && NROW(x) > 0) write.csv(x, file.path(out_dir, paste0(naam, ".csv")), row.names = FALSE)

###############################################################################
# 1. INLEZEN EN OPSCHONEN
###############################################################################
if (!file.exists(data_bestand)) data_bestand <- file.choose()
if (grepl("\\.csv$", data_bestand, ignore.case = TRUE)) {
  raw <- read.csv(data_bestand, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NULL", "-", "."))
} else {
  # Alles als tekst inlezen: kolommen met gemengde datums/tekst gaan anders deels verloren.
  # Excel-datums komen dan binnen als serienummer ("44927"); die zet parse_datum() om.
  raw <- readxl::read_excel(data_bestand, sheet = excel_sheet, .name_repair = "minimal",
                            col_types = "text", na = c("", "NA", "NULL", "-", "."))
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
  x <- tolower(trimws(as.character(x)))
  out <- as.Date(rep(NA_character_, length(x)))
  # 1) Excel-serienummers als tekst, bv. "44927" of "44927.5"
  serie <- grepl("^[0-9]{5}(\\.[0-9]+)?$", x)
  out[serie] <- as.Date(floor(as.numeric(x[serie])), origin = "1899-12-30")
  # 2) Tekstdatums: tijd eraf, maandnamen -> nummer, scheidingstekens -> "-"
  y <- sub("[ t]+[0-9]{1,2}:[0-9]{2}.*$", "", x)
  maanden <- c("januari|jan" = "01", "februari|feb" = "02", "maart|mrt|maa|mar" = "03",
               "april|apr" = "04", "mei|may" = "05", "juni|jun" = "06", "juli|jul" = "07",
               "augustus|aug" = "08", "september|sept|sep" = "09", "oktober|okt|oct" = "10",
               "november|nov" = "11", "december|dec" = "12")
  for (m in names(maanden)) y <- gsub(paste0("\\b(", m, ")\\.?\\b"), maanden[[m]], y)
  y <- gsub("[ /.]+", "-", y)
  for (f in c("%d-%m-%Y", "%Y-%m-%d", "%d-%m-%y")) {
    nieuw <- suppressWarnings(as.Date(y, format = f))
    out[is.na(out)] <- nieuw[is.na(out)]
  }
  # onwaarschijnlijke datums (typefouten) -> NA
  out[!is.na(out) & (out < as.Date("1990-01-01") | out > as.Date("2035-12-31"))] <- NA
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

# ja(): 1 bij elke niet-lege waarde behalve nee/0 (bv. "hemiboog" telt als ja)
ja <- function(x) {
  if (is.null(x)) return(0L)
  y <- to01(x); s <- tolower(trimws(as.character(x)))
  y[is.na(y) & !is.na(x) & !(s %in% c("", "nee", "n", "no", "0", "geen"))] <- 1L
  y[is.na(y)] <- 0L; y
}

# Datums
datum_cols <- grep("datum", names(raw), value = TRUE)
raw_datum_tekst <- raw$interv_datum
raw[datum_cols] <- lapply(raw[datum_cols], parse_datum)
cat("\nOperatiedatum ingelezen voor", sum(!is.na(raw$interv_datum)), "van", nrow(raw), "operaties",
    "| bereik:", format(min(raw$interv_datum, na.rm = TRUE)), "t/m", format(max(raw$interv_datum, na.rm = TRUE)), "\n")
if (any(is.na(raw$interv_datum))) {
  cat("LET OP - operatiedatum niet herkend; voorbeelden van de waarden in Excel:\n")
  print(head(unique(raw_datum_tekst[is.na(raw$interv_datum)]), 15))
}

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
ontbreekt <- character(0)
for (k in names(echo_map)) for (t in c("pre_", "post_")) {
  oud <- paste0(t, echo_map[[k]])
  if (oud %in% names(raw)) names(raw)[names(raw) == oud] <- paste0(t, k)
  else ontbreekt <- c(ontbreekt, oud)
}
if (length(ontbreekt)) {
  cat("\nLET OP - deze echo-kolommen zijn niet gevonden (worden overgeslagen):\n ",
      paste(ontbreekt, collapse = ", "), "\nKolomnamen in je bestand na opschonen:\n")
  print(names(raw))
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

# Eenheden harmoniseren PER WAARDE (cm en mm kunnen door elkaar staan)
# echo-dimensies -> cm, TAPSE -> mm
omzet <- function(x, grens, factor, richting = ">", naam = "") {
  i <- !is.na(x) & (if (richting == ">") x > grens else x < grens)
  if (any(i)) cat("  eenheid omgerekend", naam, ":", sum(i), "waarde(n) x", factor, "\n")
  x[i] <- x[i] * factor; x
}
for (t in c("pre_", "post_")) {
  for (k in c("ivsd", "lvpwd"))  if (paste0(t, k) %in% names(d)) d[[paste0(t, k)]] <- omzet(d[[paste0(t, k)]], 3, 0.1, naam = paste0(t, k))
  for (k in c("lvesd", "lvedd")) if (paste0(t, k) %in% names(d)) d[[paste0(t, k)]] <- omzet(d[[paste0(t, k)]], 15, 0.1, naam = paste0(t, k))
  v <- paste0(t, "lvot_d"); if (v %in% names(d)) d[[v]] <- omzet(d[[v]], 5, 0.1, naam = v)
  v <- paste0(t, "tapse");  if (v %in% names(d)) d[[v]] <- omzet(d[[v]], 5, 10, "<", naam = v)
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
    # Type ingreep (aortaboog-cohort): boog / overige aorta, met of zonder CABG/klep
    bijkomend = ja(cabg) == 1 | ja(klepchir) == 1,
    ingreep = case_when(
      ja(aortaboog) == 1 & bijkomend   ~ "Boog+CABG/klep",
      ja(aortaboog) == 1               ~ "Boog",
      ja(aortachir) == 1 & bijkomend   ~ "Aorta+CABG/klep",
      ja(aortachir) == 1               ~ "Aorta overig",
      ja(cabg) == 1 & ja(klepchir) == 1 ~ "CABG+klep",
      ja(cabg) == 1                    ~ "CABG",
      ja(klepchir) == 1                ~ "Klep",
      TRUE                             ~ "Overig"),
    ingreep = factor(ingreep),
    ingreep = relevel(ingreep, ref = names(which.max(table(ingreep)))),   # grootste groep = referentie
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

# Variabelen zonder variatie (bv. iedereen dezelfde ingreep) kunnen niet in een model
varieert <- names(d)[sapply(d, function(x) length(unique(na.omit(x))) > 1)]
cat("Type ingreep:\n"); print(table(d$ingreep, d$periode, useNA = "ifany"))

na  <- d %>% filter(periode == "na_jun2022")
voor <- d %>% filter(periode == "voor_jun2022")
cat("N totaal:", nrow(d), "| na juni 2022:", nrow(na), "| voor juni 2022:", nrow(voor),
    "| datum onbekend:", sum(is.na(d$periode)), "\n")
if (nrow(na) < 5) stop("Er zijn (bijna) geen operaties na ", format(cutoff_datum),
                       " gevonden. Controleer de melding over de operatiedatum hierboven en stuur die door.")


out_dir <- "output_v3"; dir.create(out_dir, showWarnings = FALSE)
out_pad <- normalizePath(out_dir)
schrijf <- function(x, naam) if (!is.null(x) && NROW(x) > 0)
  write.csv(x, file.path(out_dir, paste0(naam, ".csv")), row.names = FALSE)
for (p in c("logistf", "lavaan")) if (!heeft(p))
  try(install.packages(p, type = if (.Platform$pkgType == "source") "source" else "binary"))
set.seed(2022)

prim_par <- c("gls", "rv_s", "la_strain")          # primair
sec_par  <- c("lvef", "tapse", "e_eprime")         # secundair
alle_v2  <- c(prim_par, sec_par)
gunstig_hoger <- c(gls = TRUE, rv_s = TRUE, la_strain = TRUE, lvef = TRUE, tapse = TRUE, e_eprime = FALSE)
ref_t <- log2(ref_dag_tte + 1)

dc <- d %>% filter(!is.na(periode)) %>%
  mutate(per = as.integer(periode == "na_jun2022"),
         es2 = euroscore2,
         log_icu = log1p(icu_los),
         n_compl = rowSums(across(all_of(compl_vars)), na.rm = TRUE),
         # DAOH-30: dagen levend en buiten het ziekenhuis in de eerste 30 dagen na OK
         # alleen de 30-dagenmortaliteit-kolom (overleden_dagen_postop is niet betrouwbaar
         # als levende patiënten daar een 0 of follow-up-duur hebben)
         dood30 = ja(x30d_mortaliteit) == 1,
         daoh30 = ifelse(dood30, 0, pmax(0, 30 - pmin(in_hospital_los, 30))))

rond <- function(x) { x <- as.data.frame(x); num <- sapply(x, is.numeric)
  x[num] <- lapply(x[num], function(v) ifelse(abs(v) < 0.001 & v != 0, signif(v, 2), round(v, 3))); x }
toon <- function(x, titel) { if (is.null(x) || !NROW(x)) return(invisible())
  cat("\n==================", titel, "==================\n"); print(rond(x), row.names = FALSE) }
bayes <- function(est, se, prior_sd) {          # normaal-normaal benadering (log-schaal)
  v <- 1 / (1 / prior_sd^2 + 1 / se^2); m <- v * est / se^2
  c(post = m, laag = m - 1.96 * sqrt(v), hoog = m + 1.96 * sqrt(v), p_kleiner_0 = pnorm(0, m, sqrt(v)))
}


firth_eigen <- function(f, data) {
  mf <- model.frame(f, data); y <- model.response(mf); X <- model.matrix(f, mf)
  b <- rep(0, ncol(X))
  for (it in 1:200) {
    pr <- plogis(drop(X %*% b)); w <- pr * (1 - pr); XW <- X * sqrt(w)
    Iinv <- solve(crossprod(XW)); h <- rowSums((XW %*% Iinv) * XW)
    stap <- drop(Iinv %*% crossprod(X, y - pr + h * (0.5 - pr)))
    if (max(abs(stap)) > 5) stap <- stap * 5 / max(abs(stap))
    b <- b + stap; if (max(abs(stap)) < 1e-8) break
  }
  pr <- plogis(drop(X %*% b)); Iinv <- solve(crossprod(X * sqrt(pr * (1 - pr))))
  se <- sqrt(diag(Iinv)); names(b) <- names(se) <- colnames(X)
  list(b = b, se = se, ci_l = b - 1.96 * se, ci_h = b + 1.96 * se, p = 2 * pnorm(-abs(b / se)), methode = "Firth (Wald-BI)")
}
firth <- function(f, data) {
  if (heeft("logistf")) {
    r <- tryCatch({ fit <- logistf::logistf(f, data = data)
      b <- coef(fit); se <- sqrt(diag(fit$var)); names(se) <- names(b)
      ci_l <- fit$ci.lower; ci_h <- fit$ci.upper; names(ci_l) <- names(ci_h) <- names(b)
      list(b = b, se = se, ci_l = ci_l, ci_h = ci_h, p = setNames(fit$prob, names(b)), methode = "Firth (logistf, profiel-BI)")
    }, error = function(e) NULL)
    if (!is.null(r)) return(r)
  }
  firth_eigen(f, data)
}

cat("\n>> Bezig met M0. indeling technieken ...\n"); flush.console()
###############################################################################
# M0. BLOOTSTELLING: ITT / AT / PP
###############################################################################
tech_uit <- function(v) {
  if (!v %in% names(dc)) return(rep(NA_character_, nrow(dc)))
  x <- suppressWarnings(as.numeric(as.character(dc[[v]])))
  ifelse(is.na(x), NA, ifelse(x == code_kloppend, "kloppend_hart", "cardioplegie"))
}
dc$tech_itt <- tech_uit("perf_tech_toegewezen")
dc$tech_at  <- tech_uit("perf_tech_ontvangen")
if (all(is.na(dc$tech_itt))) { dc$tech_itt <- ifelse(dc$periode == "na_jun2022", "kloppend_hart", "cardioplegie")
  cat("perf_tech_toegewezen ontbreekt: ITT op basis van datum.\n") }
if (all(is.na(dc$tech_at))) dc$tech_at <- dc$tech_itt
dc$tech_pp <- ifelse(dc$tech_itt == dc$tech_at, dc$tech_itt, NA)
for (v in c("tech_itt", "tech_at", "tech_pp")) dc[[v]] <- factor(dc[[v]], levels = c("cardioplegie", "kloppend_hart"))
dc$omgezet <- as.integer(dc$tech_itt == "kloppend_hart" & dc$tech_at == "cardioplegie")
cat("\nToegewezen (rijen) x ontvangen (kolommen):\n"); print(table(dc$tech_itt, dc$tech_at, useNA = "ifany"))
definities <- c(ITT = "tech_itt", AT = "tech_at", PP = "tech_pp")

## Omgezette patiënten vs kloppend hart zoals bedoeld
bes_vars <- intersect(c("leeftijd_operatie", "es2", "bmi", "aox_duur", "circ_arrest_duur", "ecc_duur",
                        "pre_lvef", "pre_gls", "pre_rv_s", "min_temp"), names(dc))
bes_bin  <- intersect(c("geslacht", "ingreep", "cardiochir_eerder", "krit_preop_toestand", "endocarditis"), names(dc))
kh <- dc %>% filter(tech_itt == "kloppend_hart")
tab_M0 <- bind_rows(
  lapply(bes_vars, function(v) tibble(variabele = v,
    kloppend_zoals_bedoeld = sprintf("%.1f [%.1f-%.1f]", median(kh[[v]][kh$omgezet == 0], na.rm = TRUE),
      quantile(kh[[v]][kh$omgezet == 0], .25, na.rm = TRUE), quantile(kh[[v]][kh$omgezet == 0], .75, na.rm = TRUE)),
    omgezet_naar_cp = sprintf("%.1f [%.1f-%.1f]", median(kh[[v]][kh$omgezet == 1], na.rm = TRUE),
      quantile(kh[[v]][kh$omgezet == 1], .25, na.rm = TRUE), quantile(kh[[v]][kh$omgezet == 1], .75, na.rm = TRUE)),
    p = tryCatch(suppressWarnings(wilcox.test(kh[[v]] ~ kh$omgezet)$p.value), error = function(e) NA))),
  lapply(bes_bin, function(v) { tb <- table(kh[[v]], kh$omgezet)
    if (nrow(tb) < 2 || ncol(tb) < 2) return(NULL)
    pc <- round(100 * prop.table(tb, 2), 0)
    tibble(variabele = paste0(v, ": ", rownames(tb)), kloppend_zoals_bedoeld = paste0(tb[, "0"], " (", pc[, "0"], "%)"),
           omgezet_naar_cp = paste0(tb[, "1"], " (", pc[, "1"], "%)"),
           p = tryCatch(fisher.test(tb)$p.value, error = function(e) NA)) }))
schrijf(tab_M0, "M0_omgezette_patienten")

###############################################################################
# M1. ECHO
###############################################################################
cat("\n>> Bezig met M1. echo ...\n"); flush.console()
echo_par <- c("lvef", "gls", "rv_s", "tapse", "la_strain", "e_eprime")
res_M1 <- list()
for (def in names(definities)) for (p in echo_par) {
  tv <- definities[[def]]; dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  dd <- dc[complete.cases(dc[, c(dl, pre, "tdub", "es2", tv)]), ]
  dd$tech <- dd[[tv]]
  if (nrow(dd) < 10 || min(table(dd$tech)) < 3) next
  m <- lm(as.formula(paste(dl, "~ tech +", pre, "+ tdub + es2")), dd)
  em <- as.data.frame(summary(emmeans(m, ~ tech, at = setNames(list(mean(dd[[pre]]), ref_t), c(pre, "tdub")))))
  tt <- tidy(m, conf.int = TRUE) %>% filter(term == "techkloppend_hart")
  pg <- pnorm(0, tt$estimate, tt$std.error); gh <- gunstig_hoger[[p]]
  res_M1[[length(res_M1) + 1]] <- tibble(definitie = def, parameter = p,
    n_cp = sum(dd$tech == "cardioplegie"), n_kh = sum(dd$tech == "kloppend_hart"),
    pre_gem_cp = mean(dd[[pre]][dd$tech == "cardioplegie"]), pre_gem_kh = mean(dd[[pre]][dd$tech == "kloppend_hart"]),
    delta_cp = em$emmean[1], delta_cp_ci_laag = em$lower.CL[1], delta_cp_ci_hoog = em$upper.CL[1],
    delta_kh = em$emmean[2], delta_kh_ci_laag = em$lower.CL[2], delta_kh_ci_hoog = em$upper.CL[2],
    relatief_cp_pct = 100 * em$emmean[1] / mean(dd[[pre]][dd$tech == "cardioplegie"]),
    relatief_kh_pct = 100 * em$emmean[2] / mean(dd[[pre]][dd$tech == "kloppend_hart"]),
    verschil_kh_min_cp = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high, p = tt$p.value,
    P_kh_gunstiger = if (gh) 1 - pg else pg)
}
tab_M1 <- bind_rows(res_M1); schrijf(tab_M1, "M1_echo_per_techniek")
# delta = gecorrigeerde verandering (procentpunten voor LVEF/GLS) bij TTE op ref-dag, bij gemiddelde pre-waarde
# relatief_pct = delta als % van de gemiddelde pre-op waarde in die groep

## M1b. NIEUWE dysfunctie: normaal vóór de OK, afwijkend erna
# LV: LVEF >=50 -> <50 ; GLS: |GLS| >=18 -> <18 ; RV: TAPSE>=17 en S'>=9,5 -> TAPSE<17 of S'<9,5
nieuw_def <- list(lv_dysf = "nieuwe LV-dysfunctie (LVEF >=50 -> <50%)",
                  gls_abn = "nieuw afwijkende GLS (>=18 -> <18%)",
                  rv_dysf = "nieuwe RV-dysfunctie (TAPSE<17 of S'<9,5)")
res_M1b <- list()
for (def in names(definities)) for (k in names(nieuw_def)) {
  tv <- definities[[def]]; pre <- paste0("pre_", k); post <- paste0("post_", k)
  dd <- dc[!is.na(dc[[pre]]) & dc[[pre]] == 0 & !is.na(dc[[post]]) & !is.na(dc[[tv]]), ]
  dd$tech <- dd[[tv]]; dd$y <- dd[[post]]
  if (nrow(dd) < 10 || min(table(dd$tech)) < 3) next
  ev <- tapply(dd$y, dd$tech, sum); nn <- table(dd$tech)
  ci_bin <- function(x, n) if (n > 0) binom.test(x, n)$conf.int * 100 else c(NA, NA)
  ccp <- ci_bin(ev[["cardioplegie"]], nn[["cardioplegie"]]); ckh <- ci_bin(ev[["kloppend_hart"]], nn[["kloppend_hart"]])
  r <- if (sum(dd$y) >= 2 && !anyNA(dd$es2)) tryCatch(firth(y ~ tech + es2, dd), error = function(e) NULL) else NULL
  kk <- "techkloppend_hart"
  res_M1b[[length(res_M1b) + 1]] <- tibble(definitie = def, uitkomst = nieuw_def[[k]],
    n_risico_cp = nn[["cardioplegie"]], events_cp = ev[["cardioplegie"]],
    pct_cp = 100 * ev[["cardioplegie"]] / nn[["cardioplegie"]], ci_cp = sprintf("%.0f-%.0f%%", ccp[1], ccp[2]),
    n_risico_kh = nn[["kloppend_hart"]], events_kh = ev[["kloppend_hart"]],
    pct_kh = 100 * ev[["kloppend_hart"]] / nn[["kloppend_hart"]], ci_kh = sprintf("%.0f-%.0f%%", ckh[1], ckh[2]),
    p_fisher = tryCatch(fisher.test(table(dd$tech, dd$y))$p.value, error = function(e) NA),
    OR_firth_es2 = if (!is.null(r)) exp(r$b[[kk]]) else NA,
    OR_ci_laag = if (!is.null(r)) exp(r$ci_l[[kk]]) else NA, OR_ci_hoog = if (!is.null(r)) exp(r$ci_h[[kk]]) else NA)
}
tab_M1b <- bind_rows(res_M1b); schrijf(tab_M1b, "M1b_nieuwe_dysfunctie_per_techniek")

###############################################################################
# M2. BIOMARKERS
###############################################################################
cat("\n>> Bezig met M2. biomarkers ...\n"); flush.console()
bio <- intersect(c("ck", "ck_mb", "hstnt", "lactaat"), sub("^post_", "", names(dc)[grepl("^post_", names(dc))]))
res_M2 <- list()
for (def in names(definities)) for (b in bio) for (met_pre in c(FALSE, TRUE)) {
  tv <- definities[[def]]; y <- paste0("post_", b); pre <- paste0("pre_", b)
  if (met_pre && !pre %in% names(dc)) next
  vars <- c(y, tv, "es2", if (met_pre) pre)
  dd <- dc[complete.cases(dc[, vars]) & dc[[y]] > 0, ]
  if (met_pre) dd <- dd[dd[[pre]] > 0, ]
  dd$tech <- dd[[tv]]
  if (nrow(dd) < 10 || min(table(dd$tech)) < 3 || length(unique(dd[[y]])) < 3) next
  f <- as.formula(paste0("log(", y, ") ~ tech + es2", if (met_pre) paste0(" + log(", pre, ")")))
  m <- lm(f, dd); tt <- tidy(m, conf.int = TRUE) %>% filter(term == "techkloppend_hart")
  res_M2[[length(res_M2) + 1]] <- tibble(definitie = def, marker = b,
    correctie = if (met_pre) "ES2 + pre-waarde" else "ES2",
    n_cp = sum(dd$tech == "cardioplegie"), n_kh = sum(dd$tech == "kloppend_hart"),
    mediaan_cp = median(dd[[y]][dd$tech == "cardioplegie"]), mediaan_kh = median(dd[[y]][dd$tech == "kloppend_hart"]),
    ratio_kh_vs_cp = exp(tt$estimate), ci_laag = exp(tt$conf.low), ci_hoog = exp(tt$conf.high), p = tt$p.value,
    p_wilcoxon = suppressWarnings(wilcox.test(dd[[y]] ~ dd$tech)$p.value),
    P_lager_bij_kh = pnorm(0, tt$estimate, tt$std.error))
}
tab_M2 <- bind_rows(res_M2); schrijf(tab_M2, "M2_biomarkers_per_techniek")

###############################################################################
# M3. KLEMTIJD x TECHNIEK (dosis-respons)
###############################################################################
cat("\n>> Bezig met M3. klemtijd x techniek ...\n"); flush.console()
res_M3 <- list()
uitk_M3 <- c(paste0("delta_", echo_par), paste0("log_post_", bio))
for (b in bio) dc[[paste0("log_post_", b)]] <- ifelse(dc[[paste0("post_", b)]] > 0, log(dc[[paste0("post_", b)]]), NA)
for (def in c("ITT", "AT")) for (y in uitk_M3) {
  tv <- definities[[def]]; p <- sub("^delta_", "", y)
  pre <- if (grepl("^delta_", y)) c(paste0("pre_", p), "tdub") else character(0)
  dd <- dc[complete.cases(dc[, c(y, "aox30", "es2", tv, pre)]) & !is.na(dc$aox_duur), ]
  dd$tech <- dd[[tv]]
  if (nrow(dd) < 15 || min(table(dd$tech)) < 5) next
  m <- lm(as.formula(paste(y, "~ aox30 * tech + es2", if (length(pre)) paste("+", paste(pre, collapse = " + ")) else "")), dd)
  tr <- as.data.frame(summary(emtrends(m, ~ tech, var = "aox30")))
  p_int <- tidy(m) %>% filter(term == "aox30:techkloppend_hart") %>% pull(p.value)
  res_M3[[length(res_M3) + 1]] <- tibble(definitie = def, uitkomst = y, techniek = as.character(tr$tech),
    n = as.vector(table(dd$tech)[as.character(tr$tech)]),
    effect_per_30min = if (grepl("^log_", y)) exp(tr$aox30.trend) else tr$aox30.trend,
    ci_laag = if (grepl("^log_", y)) exp(tr$lower.CL) else tr$lower.CL,
    ci_hoog = if (grepl("^log_", y)) exp(tr$upper.CL) else tr$upper.CL,
    schaal = if (grepl("^log_", y)) "ratio per 30 min" else "verschil per 30 min", p_interactie = p_int)
}
tab_M3 <- bind_rows(res_M3); schrijf(tab_M3, "M3_klemtijd_x_techniek")

###############################################################################
# M4. KLINISCHE MYOCARDUITKOMSTEN
###############################################################################
cat("\n>> Bezig met M4. klinische myocarduitkomsten ...\n"); flush.console()
uitk_M4 <- intersect(c("post_op_cardiale_dysfunctie", "postop_mi", "ecmo_opname", "ritmeprobleem"), names(dc))
res_M4 <- list()
for (def in names(definities)) for (y in uitk_M4) {
  tv <- definities[[def]]
  dd <- dc[complete.cases(dc[, c(y, tv, "es2")]), ]; dd$tech <- dd[[tv]]
  if (sum(dd[[y]]) < 3 || min(table(dd$tech)) < 3) next
  r <- firth(as.formula(paste(y, "~ tech + es2")), dd); k <- "techkloppend_hart"
  bz <- bayes(r$b[[k]], r$se[[k]], prior_sd_zwak); bs <- bayes(r$b[[k]], r$se[[k]], prior_sd_sceptisch)
  res_M4[[length(res_M4) + 1]] <- tibble(definitie = def, uitkomst = y, n = nrow(dd), events = sum(dd[[y]]),
    pct_cp = 100 * mean(dd[[y]][dd$tech == "cardioplegie"]), pct_kh = 100 * mean(dd[[y]][dd$tech == "kloppend_hart"]),
    OR = exp(r$b[[k]]), ci_laag = exp(r$ci_l[[k]]), ci_hoog = exp(r$ci_h[[k]]), p = r$p[[k]],
    P_OR_kl1_zwak = bz[["p_kleiner_0"]], P_OR_kl1_sceptisch = bs[["p_kleiner_0"]])
}
tab_M4 <- bind_rows(res_M4); schrijf(tab_M4, "M4_klinische_myocarduitkomsten")

###############################################################################
# M5. TIMING VAN DE ECHO'S (gevoeligheidsanalyses)
###############################################################################
cat("\n>> Bezig met M5. timing echo's ...\n"); flush.console()
dc$pre_afstand <- abs(dc$pre_dagen)                 # dagen tussen pre-op TTE en OK
dc$tpre <- log2(dc$pre_afstand + 1)

## M5a. Verdeling van de timing per techniek
mq <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
tab_M5a <- bind_rows(lapply(c("ITT", "AT"), function(def) {
  tv <- definities[[def]]
  bind_rows(lapply(c("pre_afstand", "post_dagen"), function(v) {
      dd <- dc[!is.na(dc[[v]]) & !is.na(dc[[tv]]), ]
      tibble(definitie = def, tijd = ifelse(v == "pre_afstand", "pre-op TTE -> OK (dagen)", "OK -> post-op TTE (dagen)"),
        n_cp = sum(dd[[tv]] == "cardioplegie"), cardioplegie = mq(dd[[v]][dd[[tv]] == "cardioplegie"]),
        n_kh = sum(dd[[tv]] == "kloppend_hart"), kloppend_hart = mq(dd[[v]][dd[[tv]] == "kloppend_hart"]),
        pct_binnen_venster_cp = 100 * mean(if (v == "pre_afstand") dd[[v]][dd[[tv]] == "cardioplegie"] <= 90
                                           else dd[[v]][dd[[tv]] == "cardioplegie"] <= 10),
        pct_binnen_venster_kh = 100 * mean(if (v == "pre_afstand") dd[[v]][dd[[tv]] == "kloppend_hart"] <= 90
                                           else dd[[v]][dd[[tv]] == "kloppend_hart"] <= 10),
        p_wilcoxon = tryCatch(suppressWarnings(wilcox.test(dd[[v]] ~ dd[[tv]])$p.value), error = function(e) NA))
    }))
}))
schrijf(tab_M5a, "M5a_timing_echo_per_techniek")
# venster: pre <= 90 dagen; post <= 10 dagen

## Varianten van het hoofdmodel
varianten <- list(
  hoofdmodel        = list(extra = character(0), filter = function(x) rep(TRUE, nrow(x))),
  plus_tijd_pre_TTE = list(extra = "tpre",       filter = function(x) rep(TRUE, nrow(x))),
  pre_TTE_le90d     = list(extra = character(0), filter = function(x) !is.na(x$pre_afstand) & x$pre_afstand <= 90),
  post_TTE_dag1_10  = list(extra = character(0), filter = function(x) !is.na(x$post_dagen) & x$post_dagen >= 1 & x$post_dagen <= 10),
  post_TTE_dag2_10  = list(extra = character(0), filter = function(x) !is.na(x$post_dagen) & x$post_dagen >= 2 & x$post_dagen <= 10),
  pre90_en_post1_10 = list(extra = character(0), filter = function(x) !is.na(x$pre_afstand) & x$pre_afstand <= 90 &
                                                                  !is.na(x$post_dagen) & x$post_dagen >= 1 & x$post_dagen <= 10))

## M5b. Verschil kloppend hart - cardioplegie (zoals M1) per variant
res_M5b <- list()
for (def in c("ITT", "AT")) for (p in echo_par) for (vn in names(varianten)) {
  tv <- definities[[def]]; dl <- paste0("delta_", p); pre <- paste0("pre_", p); va <- varianten[[vn]]
  dd <- dc[va$filter(dc), ]
  dd <- dd[complete.cases(dd[, c(dl, pre, "tdub", "es2", tv, va$extra)]), ]; dd$tech <- dd[[tv]]
  if (nrow(dd) < 10 || min(table(dd$tech)) < 3) next
  m <- lm(as.formula(paste(dl, "~ tech +", paste(c(pre, "tdub", "es2", va$extra), collapse = " + "))), dd)
  tt <- tidy(m, conf.int = TRUE) %>% filter(term == "techkloppend_hart"); gh <- gunstig_hoger[[p]]
  pg <- pnorm(0, tt$estimate, tt$std.error)
  res_M5b[[length(res_M5b) + 1]] <- tibble(definitie = def, parameter = p, variant = vn,
    n_cp = sum(dd$tech == "cardioplegie"), n_kh = sum(dd$tech == "kloppend_hart"),
    verschil_kh_min_cp = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high, p = tt$p.value,
    P_kh_gunstiger = if (gh) 1 - pg else pg)
}
tab_M5b <- bind_rows(res_M5b); schrijf(tab_M5b, "M5b_echo_verschil_timing_varianten")

## M5c. Herstelt het hart even snel bij beide technieken? (interactie techniek x tijd tot post-op TTE)
res_M5c <- list()
for (def in c("ITT", "AT")) for (p in echo_par) {
  tv <- definities[[def]]; dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  dd <- dc[complete.cases(dc[, c(dl, pre, "tdub", "es2", tv)]), ]; dd$tech <- dd[[tv]]
  if (nrow(dd) < 15 || min(table(dd$tech)) < 5) next
  m <- lm(as.formula(paste(dl, "~ tech * tdub +", pre, "+ es2")), dd)
  tr <- as.data.frame(summary(emtrends(m, ~ tech, var = "tdub")))
  res_M5c[[length(res_M5c) + 1]] <- tibble(definitie = def, parameter = p, techniek = as.character(tr$tech),
    herstel_per_verdubbeling_dagen = tr$tdub.trend, ci_laag = tr$lower.CL, ci_hoog = tr$upper.CL,
    p_interactie = tidy(m) %>% filter(term == "techkloppend_hart:tdub") %>% pull(p.value))
}
tab_M5c <- bind_rows(res_M5c); schrijf(tab_M5c, "M5c_herstel_over_tijd_per_techniek")

## M5d. Klemtijd x techniek (zoals M3) per variant, voor LVEF, GLS en RV S'
res_M5d <- list()
for (def in c("ITT", "AT")) for (p in c("lvef", "gls", "rv_s")) for (vn in names(varianten)) {
  tv <- definities[[def]]; dl <- paste0("delta_", p); pre <- paste0("pre_", p); va <- varianten[[vn]]
  dd <- dc[va$filter(dc), ]
  dd <- dd[complete.cases(dd[, c(dl, pre, "tdub", "es2", "aox30", tv, va$extra)]), ]; dd$tech <- dd[[tv]]
  if (nrow(dd) < 15 || min(table(dd$tech)) < 5) next
  m <- lm(as.formula(paste(dl, "~ aox30 * tech +", paste(c(pre, "tdub", "es2", va$extra), collapse = " + "))), dd)
  tr <- as.data.frame(summary(emtrends(m, ~ tech, var = "aox30")))
  res_M5d[[length(res_M5d) + 1]] <- tibble(definitie = def, parameter = p, variant = vn,
    n_cp = sum(dd$tech == "cardioplegie"), n_kh = sum(dd$tech == "kloppend_hart"),
    aox_effect_cp = tr$aox30.trend[1], cp_ci = sprintf("%.2f tot %.2f", tr$lower.CL[1], tr$upper.CL[1]),
    aox_effect_kh = tr$aox30.trend[2], kh_ci = sprintf("%.2f tot %.2f", tr$lower.CL[2], tr$upper.CL[2]),
    p_interactie = tidy(m) %>% filter(term == "aox30:techkloppend_hart") %>% pull(p.value))
}
tab_M5d <- bind_rows(res_M5d); schrijf(tab_M5d, "M5d_klemtijd_x_techniek_timing_varianten")

###############################################################################
# FIGUREN
###############################################################################
pdf(file.path(out_dir, "fig_v3_myocard.pdf"), width = 8, height = 5.5)
kl <- c(cardioplegie = "#B2182B", kloppend_hart = "#2166AC")
for (b in bio) {
  y <- paste0("post_", b); dd <- dc %>% filter(!is.na(.data[[y]]), .data[[y]] > 0, !is.na(tech_itt))
  if (nrow(dd) < 10) next
  print(ggplot(dd, aes(tech_itt, .data[[y]], colour = tech_itt)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .15, alpha = .6) + scale_y_log10() + scale_colour_manual(values = kl, guide = "none") +
    labs(x = NULL, y = paste(b, "na OK (log-schaal)"), title = paste(b, ": kloppend hart vs cardioplegie (ITT)")) +
    theme_minimal())
}
for (y in c("delta_lvef", "delta_gls", "log_post_ck_mb", "log_post_ck")) {
  if (!y %in% names(dc)) next
  dd <- dc %>% filter(!is.na(.data[[y]]), !is.na(aox_duur), !is.na(tech_at)); if (nrow(dd) < 10) next
  print(ggplot(dd, aes(aox_duur, .data[[y]], colour = tech_at)) + geom_point(alpha = .6) +
    geom_smooth(method = "lm", formula = y ~ x) + scale_colour_manual(values = kl) +
    labs(x = "Aortaklemtijd (min)", y = y, colour = "Ontvangen techniek",
         title = paste("Klemtijd vs", y, "per techniek (as-treated)")) + theme_minimal())
}
dev.off()

###############################################################################
# TONEN
###############################################################################
toon(tab_M0, "M0. Omgezette patiënten (toegewezen kloppend hart, kreeg cardioplegie) vs zoals bedoeld")
toon(tab_M1 %>% select(definitie, parameter, n_cp, pre_gem_cp, delta_cp, delta_cp_ci_laag, delta_cp_ci_hoog, relatief_cp_pct,
                        n_kh, pre_gem_kh, delta_kh, delta_kh_ci_laag, delta_kh_ci_hoog, relatief_kh_pct),
     "M1. Echo: gecorrigeerde verandering PER GROEP met 95%-BI (pre-waarde + tijd tot TTE + ES2; TTE op ref-dag)")
toon(tab_M1 %>% select(definitie, parameter, verschil_kh_min_cp, ci_laag, ci_hoog, p, P_kh_gunstiger),
     "M1. Verschil kloppend hart - cardioplegie")
toon(tab_M1b, "M1b. Nieuwe dysfunctie na de OK (alleen patiënten met normale waarde vóór de OK)")
toon(tab_M2, "M2. Biomarkers na OK: ratio geometrisch gemiddelde kloppend hart / cardioplegie")
toon(tab_M3, "M3. Effect klemtijd per 30 min per techniek")
toon(tab_M4, "M4. Klinische myocarduitkomsten (Firth + Bayes)")
toon(tab_M5a, "M5a. Timing van de echo's per techniek (mediaan [IQR]; venster: pre <= 90 d, post <= 10 d)")
toon(tab_M5b %>% filter(parameter %in% c("lvef", "gls", "rv_s", "la_strain")),
     "M5b. Verschil kloppend hart - cardioplegie per timing-variant (alle parameters in CSV)")
toon(tab_M5c, "M5c. Herstel per verdubbeling van tijd tot post-op TTE, per techniek")
toon(tab_M5d, "M5d. Klemtijd (per 30 min) x techniek per timing-variant")
if (nrow(tab_M4)) cat("Methode M4:", unique(sapply(list(firth(as.formula(paste(uitk_M4[1], "~ tech_itt + es2")), dc)), `[[`, "methode")), "\n")
cat("\nITT = toegewezen, AT = ontvangen, PP = alleen wie de toegewezen techniek kreeg.\n")
cat("Alle tabellen en figuren staan in:", out_pad, "\n")

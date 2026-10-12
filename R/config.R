# =============================================================================
# Instellingen voor de analyse van de cardiochirurgische dataset
#
# Pas dit bestand aan aan JOUW dataset. De overige scripts lezen alles hieruit.
# =============================================================================

# ── Bestand ─────────────────────────────────────────────────────────────────
# Naam van je Excelbestand. R zoekt het op het bureaublad (ook in OneDrive);
# je mag hier ook een volledig pad invullen, bv. "C:/Users/naam/Desktop/data.xlsx".
# Wordt het bestand niet gevonden, dan verschijnt een venster om het te kiezen.
DATA_BESTAND  <- "AORTABOOG_NL_schoon v3(1).xlsx"
EXCEL_BLAD    <- NULL       # NULL = automatisch het werkblad met de data kiezen; of naam/nummer
UITVOER_MAP   <- "output"

# ── Selectie & groepen ──────────────────────────────────────────────────────
# Alleen patiënten in het inclusiecohort analyseren? (incl_cohort = ja/1)
ALLEEN_INCLUSIECOHORT <- TRUE
INCLUSIE_KOLOM        <- "incl_cohort"

# Kolom waarop groepen vergeleken worden in de tabellen en overlevingscurves.
# NULL = alleen het totaal. "acuut_electief" is een afgeleide variabele (zie hieronder).
GROEP <- "acuut_electief"

# ── Indeling acuut vs. electief aneurysma ───────────────────────────────────
# Acuut:               pathologie (asc/boog/desc) = 20, of urgentie = 30/40
# Electief aneurysma:  pathologie (asc/boog/desc) = 10, of urgentie = 10/20
# Pathologie gaat altijd vóór urgentie. Heeft een patiënt in het ene segment
# code 20 en in het andere 10, dan telt hij als acuut.
ACUUT_ELECTIEF <- list(
  path_kolommen      = c("aorta_asc_path", "aortaboog_path", "aorta_desc_path"),
  path_acuut         = 20,
  path_electief      = 10,
  urgentie_kolom     = "urgentie",
  urgentie_acuut     = c(30, 40),
  urgentie_electief  = c(10, 20),
  # Urgentie gebruiken als de pathologie ontbreekt. FALSE = alleen op
  # pathologie indelen.
  urgentie_als_path_ontbreekt = TRUE,
  # De dataset bestaat uit aortaboogpatiënten, dus iedereen komt in aanmerking.
  # Bevat de dataset ook andere operaties, zet dit dan op "aortachir": dan
  # wordt urgentie alleen gebruikt bij patiënten met aortachirurgie.
  aortachir_kolom    = NULL,
  labels             = c(electief = "electief aneurysma", acuut = "acuut")
)
# Patiënten die niet als acuut of electief aneurysma ingedeeld kunnen worden
# (bv. geen aortachirurgie) uit de analyse laten? FALSE = wel in 'Totaal'.
ALLEEN_INGEDEELDE_PATIENTEN <- TRUE

# ── Overleving ──────────────────────────────────────────────────────────────
OVERLEVING_TIJD   <- "overltijd_dag"   # follow-up in dagen
OVERLEVING_STATUS <- "mort_status"     # overleden ja/nee (of 1/0)
OVERLEVING_MAX_DAGEN <- 5 * 365        # tot hoever de curve getoond wordt

# ── Variabelen ──────────────────────────────────────────────────────────────
# Numerieke variabelen met plausibel bereik c(min, max). Waarden daarbuiten
# worden ontbrekend (NA) gemaakt en vastgelegd in het opschoonlogboek.
# Gebruik NA voor 'geen grens'.
NUMERIEK <- list(
  leeftijd         = c(18, 100),
  BMI              = c(12, 70),
  kreatinine       = c(1, 2000),    # µmol/L
  eGFR             = c(0, 200),     # mL/min/1.73m²
  LVEF             = c(1, 99),      # %
  EuroscoreII      = c(0.49, 100),  # %
  min_temp         = c(10, 40),     # °C
  ECC_duur         = c(0, 1000),    # min
  AoX_duur         = c(0, 600),     # min
  circ_arrest_duur = c(0, 150),     # min
  opnameduur       = c(0, 365),     # dagen
  overltijd_dag    = c(0, 20000)
)

# Categorische variabelen (ja/nee-variabelen horen hier ook bij).
# Variabelen met alleen ja/nee, j/n, yes/no, 1/0, true/false worden
# automatisch gehercodeerd naar "nee"/"ja".
CATEGORISCH <- c(
  "centrum", "jaar", "geslacht", "cardiochir_eerder", "endocarditis", "urgentie",
  "ECC", "circ_arrest", "chir_benadering", "CABG", "klepchir", "aklepchir",
  "mklepchir", "pklepchir", "tklepchir", "aortachir", "aotachir_benadering",
  "aorta_asc", "aorta_asc_path", "aortaboog_path", "aorta_desc",
  "aorta_desc_path", "cardiochir_overig", "overleden_zh", "ECMO_opname",
  "resp_insuff", "beademing", "heropname_IC", "CVA_restletsel",
  "CVA_zonder_restletsel", "nierfalen", "gastroint_compl", "vasc_compl_opname",
  "mort_status", "overl30d", "overl120d", "overl1j", "delier", "postop_MI",
  "klepchir_clean", "aortachir_clean", "circ_arrest_clean", "pathologie_groep",
  "s1_via_bewijs", "excl_reden", "incl_cohort", "klasse"
)

# Extra hercodering per variabele: c("ruwe waarde" = "nieuwe waarde").
# Hoofdletters en spaties rondom worden bij het vergelijken genegeerd.
# Labels volgens het NHR-handboek cardiochirurgie (versie 22.2.1). De code blijft
# vooraan staan, zodat de afgeleide variabelen er nog mee kunnen rekenen.
# -1 (onbekend) wordt al bij het inlezen als ontbrekend behandeld.
nhr_path <- c("0" = "0 geen", "10" = "10 aneurysma", "20" = "20 acute dissectie",
              "30" = "30 chronische dissectie", "40" = "40 infectie", "91" = "91 overig")
HERCODERING <- list(
  geslacht = c("1" = "man", "2" = "vrouw", m = "man", man = "man", male = "man",
               v = "vrouw", vrouw = "vrouw", f = "vrouw", female = "vrouw"),
  urgentie = c("10" = "10 electief", "20" = "20 urgent", "30" = "30 spoed", "40" = "40 redding"),
  ECC = c("0" = "0 geen ECC", "10" = "10 ECC soort onbekend", "20" = "20 conventionele ECC",
          "30" = "30 mini-ECC"),
  chir_benadering = c("10" = "10 sternotomie", "20" = "20 overig niet minimaal-invasief",
                      "30" = "30 mini-sternotomie", "40" = "40 mini-thoracotomie",
                      "50" = "50 percutaan", "60" = "60 overig minimaal-invasief"),
  aotachir_benadering = c("0" = "0 geen", "10" = "10 open", "20" = "20 enkelvoudige (T)EVAR",
                          "21" = "21 complexe (T)EVAR", "30" = "30 combinatie open en (T)EVAR"),
  aorta_asc = c("0" = "0 nee", "10" = "10 aortawortel geisoleerd",
                "20" = "20 ascendens incl. wortel", "30" = "30 ascendens excl. wortel",
                "90" = "90 locatie onbekend"),
  aorta_desc = c("0" = "0 nee", "1" = "1 locatie onbekend", "10" = "10 alleen thoracaal",
                 "20" = "20 alleen abdominaal", "30" = "30 thoraco-abdominaal"),
  aorta_asc_path = nhr_path, aortaboog_path = nhr_path, aorta_desc_path = nhr_path,
  gastroint_compl = c("0" = "0 geen", "1" = "1 soort onbekend", "10" = "10 bloeding",
                      "90" = "90 overig"),
  vasc_compl_opname = c("0" = "0 geen", "10" = "10 major", "20" = "20 minor")
)

# ── Tabellen ────────────────────────────────────────────────────────────────
# Welke variabelen in welke tabel komen (volgorde = volgorde in de tabel).
TABELLEN <- list(
  "Tabel 1 - Preoperatieve kenmerken" = c(
    "leeftijd", "geslacht", "BMI", "kreatinine", "eGFR", "LVEF",
    "cardiochir_eerder", "endocarditis", "urgentie", "EuroscoreII", "centrum", "jaar"
  ),
  "Tabel 2 - Operatiekenmerken" = c(
    "chir_benadering", "CABG", "klepchir", "aklepchir", "mklepchir", "pklepchir",
    "tklepchir", "aortachir", "aotachir_benadering", "aorta_asc", "aorta_asc_path",
    "aortaboog_path", "aorta_desc", "aorta_desc_path", "cardiochir_overig",
    "ECC", "ECC_duur", "AoX_duur", "circ_arrest", "circ_arrest_duur", "min_temp"
  ),
  "Tabel 3 - Postoperatieve uitkomsten" = c(
    "overleden_zh", "overl30d", "overl120d", "overl1j", "opnameduur",
    "ECMO_opname", "resp_insuff", "beademing", "heropname_IC", "CVA_restletsel",
    "CVA_zonder_restletsel", "nierfalen", "gastroint_compl", "vasc_compl_opname",
    "delier", "postop_MI"
  )
)

# Bij ja/nee-variabelen alleen de rij "ja" tonen (compacter, zoals gebruikelijk).
TOON_ALLEEN_JA <- TRUE

# ── Statistiek ──────────────────────────────────────────────────────────────
ALPHA_NORMALITEIT <- 0.05   # Shapiro-Wilk: bepaalt gem. (SD) vs. mediaan [IQR]

# ── Trendanalyse over de jaren (R/trend.R) ──────────────────────────────────
TREND_JAAR <- "jaar"

# De dataset bestaat al uit aortaboogpatiënten, dus geen extra selectie.
# TRUE = alleen rijen waarin BOOG_KOLOM is ingevuld (en niet 0) meenemen.
# Uitsplitsing naar acuut / electief gaat via GROEP.
TREND_ALLEEN_BOOG_OK <- FALSE
BOOG_KOLOM <- "aortaboog_path"

# Mortaliteit: label = c(variabele, minimale follow-up in dagen). Patiënten die
# nog leven maar korter gevolgd zijn dan deze termijn tellen voor die maat niet
# mee (anders wordt de mortaliteit in recente jaren onderschat).
TREND_MORTALITEIT <- list(
  "Ziekenhuismortaliteit" = list(var = "overleden_zh", min_followup = 0),
  "30-dagenmortaliteit"   = list(var = "overl30d",     min_followup = 30),
  "1-jaarsmortaliteit"    = list(var = "overl1j",      min_followup = 365)
)
# Mortaliteitstrend ook gecorrigeerd voor deze variabele(n) (NULL = niet)
TREND_CORRECTIE <- c("EuroscoreII")

# Complicaties; meerdere variabelen samen = 'ja' als minstens één 'ja' is.
TREND_COMPLICATIES <- list(
  "Nierfalen"                     = "nierfalen",
  "Respiratoire insufficiëntie"   = "resp_insuff",
  "Gastro-intestinale complicatie" = "gastroint_compl",
  "CVA (met of zonder restletsel)" = c("CVA_restletsel", "CVA_zonder_restletsel")
)

# Failure to rescue: overlijden (FTR_UITKOMST) bij patiënten met minstens één
# van de complicaties uit TREND_COMPLICATIES.
FTR_UITKOMST <- "overleden_zh"

TREND_TEMPERATUUR  <- "min_temp"
TREND_CIRC_ARREST  <- "circ_arrest"

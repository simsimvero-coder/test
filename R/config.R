# =============================================================================
# Instellingen voor de analyse van de cardiochirurgische dataset
#
# Pas dit bestand aan aan JOUW dataset. De overige scripts lezen alles hieruit.
# =============================================================================

# ── Bestanden ───────────────────────────────────────────────────────────────
# Pad naar de dataset: .csv (scheidingsteken ; of ,), .xlsx, .xls of .sav (SPSS)
DATA_BESTAND  <- "data/voorbeeld_cardiochirurgie.csv"
EXCEL_BLAD    <- 1          # werkblad bij Excel
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
  kreatinine       = c(20, 2000),   # µmol/L
  eGFR             = c(0, 200),     # mL/min/1.73m²
  LVEF             = c(5, 90),      # %; is LVEF bij jou een categorie, verplaats hem dan
  EuroscoreII      = c(0, 100),     # %
  min_temp         = c(10, 40),     # °C
  ECC_duur         = c(0, 1200),    # min
  AoX_duur         = c(0, 600),     # min
  circ_arrest_duur = c(0, 240),     # min
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
HERCODERING <- list(
  geslacht = c(m = "man", man = "man", male = "man", "1" = "man",
               v = "vrouw", vrouw = "vrouw", f = "vrouw", female = "vrouw", "2" = "vrouw")
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

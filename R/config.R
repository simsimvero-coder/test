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
# NULL = alleen het totaal.
GROEP <- "pathologie_groep"

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

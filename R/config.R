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

# ── PRIMAIRE indeling: originele combinatie urgentie/pathologie ──────────────
# voorrang = "urgentie" (standaard):
#   acuut              = urgentie 30/40 (ALTIJD), of pathologie 20 (acute dissectie)
#                        in asc/boog/desc (ALTIJD, ook bij urgentie 10/20)
#   electief aneurysma = urgentie 10 (geen dissectie); urgentie 20 alleen bij
#                        pathologie 10 (aneurysma)
#   niet ingedeeld     = urgentie 20 zonder pathologie 10/20, of geen bekende
#                        urgentie (behalve acute dissectie)
# voorrang = "strikt" (gevoeligheidsanalyse): urgentie én pathologie moeten de
#   groep bevestigen (electief = urgentie 10/20 + pathologie 10 zonder 20;
#   acuut = urgentie 30/40 + pathologie 10 of 20); de rest wordt uitgesloten.
# voorrang = "pathologie": pathologie 20/10 bepaalt de groep; urgentie alleen
#   als de pathologie onbekend is.
ACUUT_ELECTIEF <- list(
  voorrang           = "urgentie",   # "urgentie", "strikt" of "pathologie"
  path_kolommen      = c("aorta_asc_path", "aortaboog_path", "aorta_desc_path"),
  path_acuut         = 20,
  path_electief      = 10,
  urgentie_kolom     = "urgentie",
  urgentie_acuut     = c(30, 40),
  urgentie_electief  = c(10, 20),
  # Alleen voor "urgentie"/"pathologie": urgentie gebruiken als de pathologie
  # ontbreekt. FALSE = alleen op pathologie indelen.
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

# Ook gepoolde analyses over acuut + electief samen ("Totaal") uitvoeren?
# FALSE = elke analyse alleen apart voor acuut en electief, zonder correctie
# voor urgentie. (De kolom "Totaal" in Tabel 1-3 en het totale volume per jaar
# blijven altijd staan, omdat die het hele cohort beschrijven.)
ANALYSE_TOTAAL <- FALSE

# ── Overleving ──────────────────────────────────────────────────────────────
OVERLEVING_TIJD   <- "overltijd_dag"   # follow-up in dagen
OVERLEVING_STATUS <- "mort_status"     # overleden ja/nee (of 1/0)
OVERLEVING_MAX_DAGEN <- 5 * 365        # KM-curve tot maximaal 5 jaar
COX_MAX_DAGEN <- 365                  # Cox-trend: primair overlijden binnen 1 jaar

# ── Variabelen ──────────────────────────────────────────────────────────────
# Numerieke variabelen; bereik is alleen documentatie, er wordt NIET op afgekapt.
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
              "30" = "30 chronische dissectie", "40" = "40 infectie", "91" = "91 overig",
              "90" = "90 overig (codering voor 2022)")
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
    "leeftijd", "geslacht", "BMI", "kreatinine", "eGFR", "LVEF", "LVEF_categorie",
    "cardiochir_eerder", "endocarditis", "urgentie", "EuroscoreII", "centrum", "jaar"
  ),
  "Tabel 2 - Operatiekenmerken" = c(
    "chir_benadering", "concomitante_ingreep", "concomitant_excl_ascendens", "aantal_concomitant",
    "ingreep_wortel", "ingreep_ascendens", "CABG", "klepchir_clean", "aklepchir", "mklepchir", "pklepchir",
    "tklepchir", "aortachir", "aotachir_benadering", "aorta_asc", "aorta_asc_path",
    "aortaboog_path", "aorta_desc", "aorta_desc_path", "descendens_ingreep",
    "cardiochir_overig",
    "ECC", "ECC_duur", "AoX_duur", "circ_arrest_clean", "circ_arrest_duur",
    "circ_arrest_duur_bij_arrest", "min_temp", "temp_categorie"
  ),
  "Tabel 3 - Postoperatieve uitkomsten" = c(
    "overleden_zh", "overl30d", "overl120d", "overl1j", "opnameduur",
    "ECMO_opname", "resp_insuff", "beademing", "heropname_IC", "CVA_restletsel",
    "CVA_zonder_restletsel", "nierfalen", "gastroint_compl", "vasc_compl_opname",
    "delier", "postop_MI", "event_free_ontslag"
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
# Operatiejaren waarin meer dan dit aandeel van de patiënten (niet overleden)
# nog geen volledige follow-up KON hebben (einde database) tellen voor die
# mortaliteitsmaat niet mee. In de NHR-data: 74% voor 1 jaar in 2024, tegen
# <= 15% verlies aan follow-up in de andere jaren.
FU_MAX_ONVOLLEDIG <- 0.50
TREND_MORTALITEIT <- list(
  "Ziekenhuismortaliteit" = list(var = "overleden_zh", min_followup = 0),
  "30-dagenmortaliteit"   = list(var = "overl30d",     min_followup = 30),
  "1-jaarsmortaliteit"    = list(var = "overl1j",      min_followup = 365)
)
# In elke klinische groep eerst het ruwe en dan het gecorrigeerde model
# (covariaten: UITKOMST_COVARIATEN hieronder). NULL = alleen ruw.
TREND_CORRECTIE <- c("leeftijd", "geslacht", "cardiochir_eerder")

# Complicaties; meerdere variabelen samen = 'ja' als minstens één 'ja' is.
TREND_COMPLICATIES <- list(
  "Nierfalen"                     = "nierfalen",
  "Respiratoire insufficiëntie"   = "resp_insuff",
  "Gastro-intestinale complicatie" = "gastroint_compl",
  "CVA (met of zonder restletsel)" = c("CVA_restletsel", "CVA_zonder_restletsel")
)

# Failure to rescue: overlijden (FTR_UITKOMST) bij patiënten met minstens één
# major complicatie (MAJOR_COMPLICATIES, zie hieronder).
FTR_UITKOMST <- "overleden_zh"

TREND_TEMPERATUUR  <- "min_temp"
TREND_CIRC_ARREST  <- "circ_arrest_clean"   # zelfde variabele als CIRC_ARREST_KOLOM

# ── Volume en incidentie per miljoen inwoners (R/incidentie.R) ──────────────
# Bevolking van Nederland op 1 januari (CBS StatLine, tabel 37296ned
# "Bevolking; kerncijfers"). CONTROLEER deze cijfers bij CBS voordat je ze
# publiceert en vul nieuwe jaren aan.
BEVOLKING_1JAN <- c(
  "2013" = 16779575, "2014" = 16829289, "2015" = 16900726, "2016" = 16979120,
  "2017" = 17081507, "2018" = 17181084, "2019" = 17282163, "2020" = 17407585,
  "2021" = 17475415, "2022" = 17590672, "2023" = 17811291, "2024" = 17942942,
  "2025" = 18044027
)
# TRUE = gemiddelde bevolking van het jaar (gemiddelde van 1 januari van dit
# en het volgende jaar); FALSE = bevolking op 1 januari.
BEVOLKING_GEMIDDELD <- TRUE

# Vanaf welk jaar de indeling acuut / electief volledig is (de groepstrends per
# miljoen inwoners worden vanaf dit jaar berekend; het totaal over alle jaren).
GROEPEN_VANAF_JAAR <- 2015

# ── Vergelijking tussen periodes (R/periode.R) ──────────────────────────────
# Tabel 1-3 per periode, apart voor elke groep en voor het totaal, met p-waarden.
PERIODES <- list("2013-2018" = 2013:2018, "2019-2024" = 2019:2024)

# ── Perfusie (R/perfusie.R) ─────────────────────────────────────────────────
# Grenzen (°C) voor de temperatuurcategorieën: <= 20 diep, 20-28 matig, > 28.
TEMP_GRENZEN <- c(20, 28)
# Variabele die bepaalt of er circulatiestilstand was (voor de duur bij
# patiënten met stilstand). circ_arrest_clean = "ja" bij een geregistreerde
# stilstand óf een ingevulde stilstandduur > 0 (die patiënten hebben ook de lage
# temperaturen van patiënten met stilstand).
CIRC_ARREST_KOLOM <- "circ_arrest_clean"
# Kolom met het centrum, voor de random-intercept- en fixed-effects-modellen.
CENTRUM_KOLOM <- "centrum"

# Delier wordt pas vanaf 2022 systematisch geregistreerd (daarvoor > 95% leeg).
DELIER_VANAF_JAAR <- 2022
# Aantal bootstraptrekkingen voor de mediatieanalyse (nierfalen). Meer = trager.
MEDIATIE_BOOTSTRAP <- 200
# Bootstraptrekkingen voor de gestandaardiseerde mortaliteit per jaar (R/trendverschil.R).
TRENDVERSCHIL_BOOTSTRAP <- 200

# ── Klinische uitkomsten en variatie tussen centra (R/uitkomsten.R) ─────────
# Samengestelde uitkomst "major complicaties": ja als minstens één van deze ja is
# (gastroint_compl: elke code > 0). Ook de noemer voor failure to rescue.
MAJOR_COMPLICATIES <- c("CVA_restletsel", "nierfalen", "gastroint_compl")

# Secundaire morbiditeitsuitkomsten: Benjamini-Hochberg-correctie over deze
# uitkomsten, apart per groep (acuut / electief).
SECUNDAIRE_UITKOMSTEN <- c(
  "CVA met restletsel"             = "CVA_restletsel",
  "Nierfalen"                      = "nierfalen",
  "Gastro-intestinale complicatie" = "gastroint_compl",
  "Respiratoire insufficiëntie"    = "resp_insuff",
  "Beademing > 24 uur"             = "beademing"
)

# Correctievariabelen in alle patientgebonden modellen; per groep 2 modellen
# (ongecorrigeerd + gecorrigeerd). Leeftijd wordt per 10 jaar ingevoerd.
# Primair zonder EuroSCORE II (geen correctie voor urgentie binnen de groep).
UITKOMST_COVARIATEN <- c("leeftijd", "geslacht", "cardiochir_eerder")

# Gevoeligheidsanalyse: dezelfde modellen, gecorrigeerd voor EuroSCORE II.
# Alleen EuroSCORE II, omdat leeftijd, geslacht en eerdere cardiochirurgie daar
# al in zitten (anders dubbele correctie). Uitvoer: gevoeligheid_euroscore.xlsx.
# NULL = overslaan.
GEVOELIGHEID_EUROSCORE_COVARIATEN <- c("EuroscoreII")
# Alleen jaren waarin EuroSCORE II vrijwel compleet is geregistreerd (vanaf 2018
# ontbreekt het bij <= 13% van de patiënten; in 2013-2017 bij 27-82%). In de
# vroege jaren zou een complete-case-analyse een selectie van patiënten zijn.
GEVOELIGHEID_EUROSCORE_JAREN <- 2018:2024
PERFUSIE_COVARIATEN <- UITKOMST_COVARIATEN
OVERLEVING_COVARIATEN <- UITKOMST_COVARIATEN

# Gevoeligheidsanalyse met centrum-fixed-effects en CR2 (traag; FALSE = overslaan)
UITKOMST_CR2 <- TRUE

# Relatieve overleving (Ederer II): CSV-bestand met de Nederlandse sterftekansen
# (CBS StatLine, tabel 37360ned "Levensverwachting; geslacht, leeftijd").
# Kolommen: jaar, geslacht (man/vrouw), leeftijd (0-99), sterftekans (kans om binnen
# een jaar te overlijden, 0-1). Zoek je het bestand op het bureaublad, geef dan
# alleen de naam. NULL = relatieve overleving overslaan.
# Gevalideerd tegen relsurv::rs.surv(method = "ederer2"): verschil <= 0,3 %-punt.
LEVENSTABEL_BESTAND <- "levenstabel_NL.csv"

# ── Benchmark per centrum: funnel plots (R/benchmark.R) ─────────────────────
# Verwachte sterfte (EuroSCORE II) herkalibreren naar het landelijke niveau,
# zodat O/E = 1 het landelijke niveau bij dezelfde case-mix betekent.
FUNNEL_HERKALIBREREN <- TRUE

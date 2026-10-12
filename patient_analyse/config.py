"""
Instellingen voor de analyse van de patiëntendataset.

Pas dit bestand aan zodat het past bij JOUW dataset. De rest van de code
leest alles hieruit; meestal hoef je verder niets te wijzigen.
"""

# ── Bestanden ────────────────────────────────────────────────────
# Pad naar de dataset (.csv, .xlsx of .xls). Kan ook via de opdrachtregel
# worden opgegeven: python -m patient_analyse.main --bestand pad/naar/data.csv
DATA_BESTAND = "data/voorbeeld_patienten.csv"

# Werkblad bij Excel-bestanden (naam of index; 0 = eerste blad)
EXCEL_WERKBLAD = 0

# Map waarin tabellen, figuren en het opschoonlogboek worden opgeslagen
UITVOER_MAP = "output"

# ── Kolommen ─────────────────────────────────────────────────────
# Kolom met het (gepseudonimiseerde!) patiëntnummer. None = geen ID-kolom.
ID_KOLOM = "patient_id"

# Kolom waarop groepen vergeleken worden (bv. behandeling, cohort, uitkomst).
# None = geen groepsvergelijking, alleen beschrijvende statistiek.
GROEP_KOLOM = "behandeling"

# Numerieke variabelen met hun plausibele bereik (min, max).
# Waarden buiten dit bereik worden als ontbrekend (NaN) gemarkeerd en
# vastgelegd in het opschoonlogboek. Gebruik None als er geen grens is.
NUMERIEK = {
    "leeftijd": (18, 110),
    "bmi": (12, 70),
    "systolische_bd": (60, 260),
    "diastolische_bd": (30, 160),
    "hba1c": (15, 200),          # mmol/mol
    "opnameduur_dagen": (0, 365),
}

# Categorische variabelen. Per variabele een optionele hercodering:
# links de ruwe waarde in kleine letters (hoofdletters en spaties rondom
# worden genegeerd), rechts de gewenste waarde. Waarden die niet in de
# hercodering staan blijven ongewijzigd. Lege cellen, 'onbekend' en '?'
# worden als ontbrekend behandeld.
CATEGORISCH = {
    "geslacht": {
        "m": "man", "man": "man", "male": "man",
        "v": "vrouw", "vrouw": "vrouw", "f": "vrouw", "female": "vrouw",
    },
    "behandeling": {},
    "roker": {"ja": "ja", "j": "ja", "yes": "ja", "1": "ja",
              "nee": "nee", "n": "nee", "no": "nee", "0": "nee"},
    "diagnose": {},
    "heropname_30d": {"ja": "ja", "1": "ja", "nee": "nee", "0": "nee"},
}

# Datumvariabelen (worden ingelezen als dag-maand-jaar).
DATUMS = ["opnamedatum"]

# ── Statistiek ───────────────────────────────────────────────────
# Significantieniveau voor de normaliteitstoets (Shapiro-Wilk) die bepaalt
# of een parametrische of non-parametrische toets wordt gebruikt.
ALPHA_NORMALITEIT = 0.05

# Minimum aantal waarnemingen per groep om een toets uit te voeren.
MIN_N_PER_GROEP = 3

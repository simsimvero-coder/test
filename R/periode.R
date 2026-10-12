# =============================================================================
# Vergelijking tussen periodes (bv. 2013-2018 vs. 2019-2024)
#
# Maakt Tabel 1-3 met de periode als groepsvariabele, apart voor elke groep
# (acuut / electief) en voor het totaal. Toetsen en weergave zijn dezelfde als
# in de hoofdtabellen (zie maak_tabel in R/analyse.R).
# =============================================================================

# Variabelen die in de console getoond worden (alles staat in resultaten.xlsx)
PERIODE_KERNVARIABELEN <- c("leeftijd", "geslacht", "EuroscoreII", "urgentie", "klepchir",
                            "klepchir_clean", "aorta_asc_path", "aortaboog_path",
                            "descendens_ingreep", "descendens_pathologie", "circ_arrest",
                            "min_temp", "overleden_zh", "overl30d")

periodevergelijking <- function(df, groep = GROEP, jaar_kolom = TREND_JAAR) {
  jaar <- suppressWarnings(as.numeric(as.character(df[[jaar_kolom]])))
  periode <- rep(NA_character_, nrow(df))
  for (naam in names(PERIODES)) periode[jaar %in% PERIODES[[naam]]] <- naam
  df$periode <- factor(periode, levels = names(PERIODES))
  if (any(is.na(df$periode)))
    message("LET OP - ", sum(is.na(df$periode)), " patiënten vallen buiten PERIODES en tellen niet mee")

  variabelen <- unique(unlist(TABELLEN))
  subsets <- if (!is.null(groep) && groep %in% names(df))
    split(df, droplevels(factor(df[[groep]]))) else list()
  subsets[["Totaal"]] <- df
  lapply(subsets, function(d) maak_tabel(d[!is.na(d$periode), ], variabelen, "periode"))
}

# Rijen van een tabel voor een selectie variabelen (incl. de categorierijen eronder)
selecteer_rijen <- function(tabel, variabelen) {
  kop <- sub(",.*$", "", tabel$Variabele)
  is_kop <- !grepl("^    ", tabel$Variabele)
  blok <- cumsum(is_kop)
  tabel[blok %in% blok[is_kop & kop %in% variabelen], ]
}

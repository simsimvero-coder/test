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
  if (ANALYSE_TOTAAL || !length(subsets)) subsets[["Totaal"]] <- df
  lapply(subsets, function(d) maak_tabel(d[!is.na(d$periode), ], variabelen, "periode"))
}

# Rijen van een tabel voor een selectie variabelen (incl. de categorierijen eronder)
selecteer_rijen <- function(tabel, variabelen) {
  kop <- sub(",.*$", "", tabel$Variabele)
  is_kop <- !grepl("^    ", tabel$Variabele)
  blok <- cumsum(is_kop)
  tabel[blok %in% blok[is_kop & kop %in% variabelen], ]
}

# Supplement bij de beschrijvende periode-tabellen:
# per klinische groep en eindpunt twee modellen (periode laat versus vroeg),
# zonder en met de 4 vooraf vastgelegde correctievariabelen.
periode_uitkomstmodellen <- function(df, groep = GROEP) {
  if (length(PERIODES) != 2 || !all(c(groep, TREND_JAAR, CENTRUM_KOLOM) %in% names(df))) return(NULL)
  jaar <- suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]])))
  per <- ifelse(jaar %in% PERIODES[[1]], names(PERIODES)[1],
                ifelse(jaar %in% PERIODES[[2]], names(PERIODES)[2], NA_character_))
  basis <- data.frame(groep = as.character(df[[groep]]),
                      periode = factor(per, levels = names(PERIODES)),
                      .centrum = factor(df[[CENTRUM_KOLOM]]))
  basis <- voeg_modelcovariaten_toe(basis, df)
  m <- c("Ziekenhuismortaliteit" = "overleden_zh", "30-dagenmortaliteit" = "overl30d",
         "1-jaarsmortaliteit" = "overl1j")
  m <- m[unname(m) %in% names(df)]
  if (!length(m)) return(NULL)
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  rijen <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (uitkomst in names(m)) {
    y <- code_positief(df[[m[[uitkomst]]]])
    termijn <- if (uitkomst == "30-dagenmortaliteit") 30 else
               if (uitkomst == "1-jaarsmortaliteit") 365 else 0
    if (termijn) y[!is.na(y) & !y & (is.na(tijd) | tijd < termijn)] <- NA
    for (gecorrigeerd in c(FALSE, TRUE)) {
      d <- basis; d$y <- y
      d <- droplevels(d[d$groep == g & !is.na(d$y) & !is.na(d$periode) &
                          !is.na(d$.centrum), , drop = FALSE])
      if (gecorrigeerd) d <- modelgevallen(d, "y")
      or <- lo <- hi <- pval <- NA_real_
      status <- "onvoldoende observaties, events of beide periodes"
      if (nrow(d) >= 30 && sum(d$y) >= 5 && sum(!d$y) >= 5 &&
          nlevels(d$periode) == 2) {
        term <- paste0("periode", levels(d$periode)[2])
        r <- tryCatch(logit_mixed(cov_formule("y", "periode", gecorrigeerd), d, term),
                      error = function(e) list(model = paste("mislukt:", conditionMessage(e))))
        status <- r$model
        if (!is.null(r$est) && is.finite(r$est) && is.finite(r$se)) {
          or <- exp(r$est); lo <- exp(r$est - 1.96 * r$se)
          hi <- exp(r$est + 1.96 * r$se); pval <- r$p
        }
      }
      rijen[[length(rijen) + 1]] <- data.frame(
        groep = g, uitkomst = uitkomst, vergelijking = paste(names(PERIODES)[2], "vs", names(PERIODES)[1]),
        correctie = if (gecorrigeerd) "gecorrigeerd" else "ongecorrigeerd",
        n = nrow(d), events = sum(d$y), OR = round(or, 3),
        BI95_onder = round(lo, 3), BI95_boven = round(hi, 3),
        p = fmt_p(pval), model_status = status)
    }
  }
  do.call(rbind, rijen)
}

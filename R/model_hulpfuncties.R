# Gedeelde covariaten voor TWEE modellen binnen iedere klinische groep.
# Model 1: ongecorrigeerd. Model 2: de covariaten uit UITKOMST_COVARIATEN in
# R/config.R (leeftijd wordt per 10 jaar ingevoerd).
MODEL_TERMS <- sub("^leeftijd$", "leeftijd_10", UITKOMST_COVARIATEN)

voeg_modelcovariaten_toe <- function(d, bron) {
  nodig <- UITKOMST_COVARIATEN
  ontbreekt <- setdiff(nodig, names(bron))
  if (length(ontbreekt)) {
    warning("Geen volledig gecorrigeerd model mogelijk; ontbrekende kolommen: ",
            paste(ontbreekt, collapse = ", "))
    for (k in ontbreekt) bron[[k]] <- rep(NA, nrow(bron))
  }
  for (k in nodig) {
    x <- bron[[k]]
    if (k == "leeftijd") d$leeftijd_10 <- suppressWarnings(as.numeric(as.character(x))) / 10
    else if (is.numeric(x)) d[[k]] <- x
    else d[[k]] <- factor(as.character(x))
  }
  d
}

modelgevallen <- function(d, y = "y", cov = MODEL_TERMS) {
  nodig <- intersect(c(y, cov), names(d))
  if (length(nodig) != length(unique(c(y, cov)))) return(d[FALSE, , drop = FALSE])
  droplevels(d[stats::complete.cases(d[, nodig, drop = FALSE]), , drop = FALSE])
}

# zonder: covariaten die in dit model niet thuishoren (bv. EuroSCORE II bij een
# vergelijking tussen urgentiegroepen, omdat urgentie in EuroSCORE II zit).
cov_formule <- function(respons, tijd, correctie = FALSE, extra = NULL, zonder = NULL) {
  stats::reformulate(c(tijd, extra, if (correctie) setdiff(MODEL_TERMS, zonder)), response = respons)
}

# Voert 'expr' uit met tijdelijk andere correctievariabelen (bv. alleen
# EuroSCORE II voor de gevoeligheidsanalyse) en zet daarna alles terug.
met_covariaten <- function(cov, expr) {
  oud <- list(UITKOMST_COVARIATEN = UITKOMST_COVARIATEN, MODEL_TERMS = MODEL_TERMS)
  on.exit(for (k in names(oud)) assign(k, oud[[k]], envir = globalenv()))
  assign("UITKOMST_COVARIATEN", cov, envir = globalenv())
  assign("MODEL_TERMS", sub("^leeftijd$", "leeftijd_10", cov), envir = globalenv())
  expr
}

# =============================================================================
# Afgeleide variabelen
# =============================================================================

# Deelt patiënten in als "acuut" of "electief aneurysma" (zie ACUUT_ELECTIEF in
# config.R). Voegt de kolommen acuut_electief en acuut_electief_basis toe.
maak_acuut_electief <- function(df, cfg = ACUUT_ELECTIEF) {
  # Leest de code vooraan, zodat zowel "30" als "30 spoed" werkt
  als_getal <- function(x) suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", as.character(x))))
  n <- nrow(df)

  path_k <- intersect(cfg$path_kolommen, names(df))
  if (length(path_k) < length(cfg$path_kolommen))
    message("LET OP - pathologiekolommen niet gevonden: ",
            paste(setdiff(cfg$path_kolommen, path_k), collapse = ", "))
  path <- sapply(path_k, function(k) als_getal(df[[k]]))
  path <- matrix(path, nrow = n)
  path_acuut    <- rowSums(path == cfg$path_acuut, na.rm = TRUE) > 0
  path_electief <- rowSums(path == cfg$path_electief, na.rm = TRUE) > 0

  urg <- if (cfg$urgentie_kolom %in% names(df)) als_getal(df[[cfg$urgentie_kolom]]) else rep(NA, n)
  urg_acuut    <- urg %in% cfg$urgentie_acuut
  urg_electief <- urg %in% cfg$urgentie_electief
  aorta <- if (length(cfg$aortachir_kolom) && cfg$aortachir_kolom %in% names(df))
    df[[cfg$aortachir_kolom]] %in% "ja" else rep(TRUE, n)
  mag_urgentie <- cfg$urgentie_als_path_ontbreekt & aorta & !path_acuut & !path_electief

  groep <- rep(NA_character_, n)
  basis <- rep(NA_character_, n)
  # Volgorde is van belang: pathologie overschrijft urgentie, acuut gaat vóór electief.
  sel <- mag_urgentie & urg_electief; groep[sel] <- "electief"; basis[sel] <- "urgentie (pathologie onbekend)"
  sel <- mag_urgentie & urg_acuut;    groep[sel] <- "acuut";    basis[sel] <- "urgentie (pathologie onbekend)"
  groep[path_electief] <- "electief"; basis[path_electief] <- "pathologie"
  groep[path_acuut]    <- "acuut";    basis[path_acuut]    <- "pathologie"

  # Pathologie en urgentie spreken elkaar tegen -> pathologie wint, wel melden
  conflict <- (path_acuut & urg_electief) | (!path_acuut & path_electief & urg_acuut)
  basis[conflict] <- ifelse(path_acuut[conflict],
                            "pathologie (overschrijft electieve urgentie)",
                            "pathologie (overschrijft acute urgentie)")

  df$acuut_electief <- factor(unname(cfg$labels[groep]), levels = unname(cfg$labels))
  df$acuut_electief_basis <- basis

  overzicht <- as.data.frame(table(groep = addNA(df$acuut_electief),
                                   basis = addNA(df$acuut_electief_basis)),
                             responseName = "n")
  overzicht <- overzicht[overzicht$n > 0, ]
  overzicht$groep <- as.character(overzicht$groep)
  overzicht$groep[is.na(overzicht$groep)] <- "niet ingedeeld"
  overzicht$basis <- as.character(overzicht$basis)
  overzicht$basis[is.na(overzicht$basis)] <- "-"

  # Kruistabel om de indeling te controleren: urgentie x pathologie x groep
  controle <- as.data.frame(table(
    urgentie = addNA(factor(urg)),
    pathologie = addNA(factor(ifelse(path_acuut & path_electief, "20 en 10",
                                ifelse(path_acuut, cfg$path_acuut,
                                ifelse(path_electief, cfg$path_electief, NA))))),
    groep = addNA(df$acuut_electief)), responseName = "n")
  controle <- controle[controle$n > 0, ]
  for (k in c("urgentie", "pathologie", "groep")) {
    controle[[k]] <- as.character(controle[[k]])
    controle[[k]][is.na(controle[[k]])] <- "ontbreekt / anders"
  }

  list(data = df, overzicht = overzicht, controle = controle)
}

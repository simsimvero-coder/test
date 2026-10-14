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
  geen_path <- !path_acuut & !path_electief
  alle_path_onbekend <- rowSums(!is.na(path)) == 0
  zet <- function(sel, g, tekst) { groep[sel] <<- g; basis[sel] <<- tekst }

  if (identical(cfg$voorrang, "strikt")) {
    # Strikt (onderzoeksprotocol):
    #   - acute dissectie (pathologie 20 in asc/boog/desc) = altijd acuut, ook bij
    #     urgentie 10/20 of onbekende urgentie
    #   - aneurysma (pathologie 10, geen 20): urgentie 10/20 = electief, 30/40 = acuut
    #   - pathologie niet ingevuld (alle 3 leeg): indeling op urgentie
    #     (10/20 = electief, 30/40 = acuut)
    #   - uitgesloten: alleen andere pathologie (0/30/40/90/91, geen 10/20) en
    #     patiënten zonder bruikbare urgentie (behalve acute dissectie)
    # Latere regels overschrijven eerdere.
    zet(!is.na(urg) & !urg_acuut & !urg_electief, NA, "uitgesloten: urgentie niet 10-40")
    zet(is.na(urg), NA, "uitgesloten: urgentie onbekend")
    zet(geen_path & !alle_path_onbekend, NA, "uitgesloten: andere pathologie (0/30/40/90/91)")
    zet(alle_path_onbekend & urg_electief, "electief", "pathologie niet ingevuld + urgentie 10/20")
    zet(alle_path_onbekend & urg_acuut, "acuut", "pathologie niet ingevuld + urgentie 30/40")
    zet(path_electief & !path_acuut & urg_electief, "electief", "pathologie 10 + urgentie 10/20")
    zet(path_electief & !path_acuut & urg_acuut, "acuut", "pathologie 10 + urgentie 30/40")
    zet(path_acuut & urg_acuut, "acuut", "pathologie 20 + urgentie 30/40")
    zet(path_acuut & !urg_acuut, "acuut", "pathologie 20 + urgentie 10/20 of onbekend")
  } else if (identical(cfg$voorrang, "urgentie")) {
    # Urgentie 30/40 is altijd acuut; een acute dissectie (pathologie 20) ook.
    # Volgorde is van belang: latere regels overschrijven eerdere.
    zet(mag_urgentie & urg_electief, "electief", "urgentie 10/20 (geen pathologie 10/20)")
    zet(path_electief & !urg_electief, "electief", "pathologie 10 (urgentie onbekend)")
    zet(path_electief & urg_electief, "electief", "pathologie 10 + urgentie 10/20")
    zet(path_acuut & !urg_acuut, "acuut", "pathologie 20 (urgentie 10/20 of onbekend)")
    zet(urg_acuut & geen_path, "acuut", "urgentie 30/40 (geen pathologie 10/20)")
    zet(urg_acuut & path_electief & !path_acuut, "acuut", "urgentie 30/40 (overschrijft pathologie 10)")
    zet(urg_acuut & path_acuut, "acuut", "urgentie 30/40 + pathologie 20")
  } else {
    # Pathologie gaat voor urgentie; urgentie alleen als de pathologie onbekend is.
    zet(mag_urgentie & urg_electief, "electief", "urgentie (pathologie onbekend)")
    zet(mag_urgentie & urg_acuut, "acuut", "urgentie (pathologie onbekend)")
    zet(path_electief, "electief", "pathologie")
    zet(path_acuut, "acuut", "pathologie")
    conflict <- (path_acuut & urg_electief) | (!path_acuut & path_electief & urg_acuut)
    basis[conflict] <- ifelse(path_acuut[conflict],
                              "pathologie (overschrijft electieve urgentie)",
                              "pathologie (overschrijft acute urgentie)")
  }

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
  overzicht <- overzicht[order(overzicht$groep == "niet ingedeeld", overzicht$groep, -overzicht$n), ]

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

# LVEF in de EuroSCORE II-categorieën. De NHR codeert beschrijvende waarden als
# vaste percentages (goed = 55, matig = 40, slecht = 25, zeer slecht = 20), dus
# een indeling in categorieën is informatiever dan mediaan [IQR].
maak_lvef_categorie <- function(df, kolom = "LVEF") {
  if (!kolom %in% names(df) || !is.numeric(df[[kolom]])) return(df)
  df$LVEF_categorie <- cut(df[[kolom]], breaks = c(-Inf, 20, 30, 50, Inf),
                           labels = c("zeer slecht (<=20%)", "slecht (21-30%)",
                                      "matig (31-50%)", "goed (>50%)"))
  df
}

# Ja/nee-variabelen voor de aorta descendens:
#  - descendens_ingreep:    aorta_desc > 0 (ingreep aan de descendens; bijna compleet)
#  - descendens_pathologie: aorta_desc_path > 0 (pathologie geregistreerd; vaker ontbrekend)
maak_descendens <- function(df) {
  code <- function(x) suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", as.character(x))))
  janee <- function(x) factor(ifelse(is.na(x), NA, ifelse(x > 0, "ja", "nee")), levels = c("nee", "ja"))
  if ("aorta_desc" %in% names(df)) df$descendens_ingreep <- janee(code(df$aorta_desc))
  if ("aorta_desc_path" %in% names(df)) df$descendens_pathologie <- janee(code(df$aorta_desc_path))
  df
}

# Temperatuurcategorie en duur van de circulatiestilstand bij patiënten met stilstand
maak_perfusie <- function(df) {
  if ("min_temp" %in% names(df) && is.numeric(df$min_temp)) {
    g <- TEMP_GRENZEN
    df$temp_categorie <- cut(df$min_temp, breaks = c(-Inf, g[1], g[2], Inf),
                             labels = c(sprintf("diep (<=%g C)", g[1]),
                                        sprintf("matig (%g-%g C)", g[1], g[2]),
                                        sprintf("> %g C", g[2])))
  }
  kolom <- if (CIRC_ARREST_KOLOM %in% names(df)) CIRC_ARREST_KOLOM else "circ_arrest"
  if (all(c(kolom, "circ_arrest_duur") %in% names(df))) {
    met_arrest <- df[[kolom]] %in% "ja" & !is.na(df$circ_arrest_duur) & df$circ_arrest_duur > 0
    df$circ_arrest_duur_bij_arrest <- ifelse(met_arrest, df$circ_arrest_duur, NA)
  }
  df
}

# Hoeveel geregistreerde patiënten per jaar en per centrum in elke groep vallen
# en hoeveel er worden uitgesloten. Uitsluiting die per jaar of centrum verschilt
# kan trends en de benchmark vertekenen.
selectie_overzicht <- function(df_alle, groep = GROEP) {
  g <- as.character(df_alle[[groep]]); g[is.na(g)] <- "uitgesloten"
  per <- function(kolom) {
    tab <- table(df_alle[[kolom]], factor(g, levels = c(unname(ACUUT_ELECTIEF$labels), "uitgesloten")))
    out <- data.frame(x = rownames(tab), geregistreerd = as.integer(rowSums(tab)), check.names = FALSE)
    for (k in colnames(tab)) out[[k]] <- as.integer(tab[, k])
    out$uitgesloten_pct <- round(100 * out$uitgesloten / out$geregistreerd, 1)
    names(out)[1] <- kolom
    if (kolom == CENTRUM_KOLOM) out <- out[order(suppressWarnings(as.numeric(out[[1]]))), ]
    out
  }
  list(jaar = per(TREND_JAAR), centrum = per(CENTRUM_KOLOM))
}

# Gevoeligheid voor de definitie: groepsgrootte en ruwe vroege mortaliteit per
# periode onder elke indeling (strikt / urgentie / pathologie).
vergelijk_definities <- function(df_alle, definities = c("strikt", "urgentie", "pathologie")) {
  jaar <- suppressWarnings(as.numeric(as.character(df_alle[[TREND_JAAR]])))
  periode <- rep(NA_character_, nrow(df_alle))
  for (p in names(PERIODES)) periode[jaar %in% PERIODES[[p]]] <- p
  pct <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA else round(100 * mean(x == "ja"), 1) }
  do.call(rbind, lapply(definities, function(def) {
    cfg <- ACUUT_ELECTIEF; cfg$voorrang <- def
    gr <- as.character(maak_acuut_electief(df_alle, cfg)$data$acuut_electief)
    sel <- expand.grid(periode = c(names(PERIODES), "alle jaren"), groep = unname(cfg$labels),
                       stringsAsFactors = FALSE)
    do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
      rij <- gr %in% sel$groep[i] & (sel$periode[i] == "alle jaren" | periode %in% sel$periode[i])
      data.frame(definitie = def, groep = sel$groep[i], periode = sel$periode[i], n = sum(rij),
                 ziekenhuismortaliteit_pct = pct(df_alle$overleden_zh[rij]),
                 mortaliteit_30d_pct = pct(df_alle$overl30d[rij]))
    }))
  }))
}

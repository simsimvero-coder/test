# =============================================================================
# Inlezen en opschonen
# =============================================================================

# NHR-handboek: -1 (of -1,0 / -1,00) = onbekend
ONTBREKEND_CODES <- c("", "na", "n/a", "nan", "onbekend", "unknown", "missing", "?", "-",
                      "-1", "-1.0", "-1,0", "-1.00", "-1,00")
JA_CODES  <- c("ja", "j", "yes", "y", "1", "true", "waar", "overleden", "dood")
NEE_CODES <- c("nee", "n", "no", "0", "false", "onwaar", "levend", "in leven", "alive")

# Zoekt het bestand: als volledig pad, in de werkmap, of op het bureaublad
# (Windows: Desktop/Bureaublad, ook binnen OneDrive; Mac/Linux: ~/Desktop).
zoek_bestand <- function(naam) {
  if (is.null(naam) || !nzchar(naam)) return(NULL)
  if (file.exists(naam)) return(normalizePath(naam))
  thuis <- c(Sys.getenv("USERPROFILE"), Sys.getenv("HOME"), path.expand("~"))
  thuis <- unique(thuis[nzchar(thuis)])
  onedrive <- unlist(lapply(thuis, function(h) list.files(h, pattern = "^OneDrive", full.names = TRUE)))
  mappen <- c(file.path(rep(c(thuis, onedrive), each = 2), c("Desktop", "Bureaublad")),
              file.path(dirname(path.expand("~")), "Desktop"))
  kandidaten <- file.path(unique(mappen), basename(naam))
  gevonden <- kandidaten[file.exists(kandidaten)]
  if (length(gevonden)) normalizePath(gevonden[1]) else NULL
}

lees_data <- function(pad, blad = 1) {
  if (!file.exists(pad)) stop("Bestand niet gevonden: ", pad)
  if (!tolower(tools::file_ext(pad)) %in% c("xlsx", "xls"))
    stop("Kies een Excelbestand (.xlsx of .xls): ", pad)
  # Alles als tekst inlezen; het opschonen zet de typen daarna goed
  df <- as.data.frame(readxl::read_excel(pad, sheet = blad, col_types = "text"))

  names(df) <- trimws(names(df))
  # Kolomnamen hoofdletter-ongevoelig koppelen aan de namen uit config.R
  bekend <- unique(c(names(NUMERIEK), CATEGORISCH))
  for (i in seq_along(names(df))) {
    match <- bekend[tolower(bekend) == tolower(names(df)[i])]
    if (length(match) == 1) names(df)[i] <- match
  }
  df
}

schoon_op <- function(df) {
  logboek <- list()
  noteer <- function(stap, kolom, rij, oud, nieuw) {
    if (length(rij) == 0) return(invisible())
    # rij + 1: regel 1 van het bestand is de kopregel
    logboek[[length(logboek) + 1]] <<- data.frame(
      stap = stap, kolom = kolom %||% NA_character_, rij_in_bestand = rij + 1,
      oude_waarde = as.character(oud), nieuwe_waarde = as.character(nieuw)
    )
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a

  df$.rij <- seq_len(nrow(df))
  for (k in setdiff(names(df), ".rij")) {
    x <- trimws(as.character(df[[k]]))
    x[tolower(x) %in% ONTBREKEND_CODES] <- NA
    df[[k]] <- x
  }

  # Lege en dubbele rijen
  leeg <- rowSums(!is.na(df[setdiff(names(df), ".rij")])) == 0
  noteer("lege rij verwijderd", NULL, df$.rij[leeg], NA, NA)
  df <- df[!leeg, ]
  dubbel <- duplicated(df[setdiff(names(df), ".rij")])
  noteer("dubbele rij verwijderd", NULL, df$.rij[dubbel], NA, NA)
  df <- df[!dubbel, ]

  ontbreekt <- setdiff(unique(c(names(NUMERIEK), CATEGORISCH)), names(df))
  if (length(ontbreekt)) message("LET OP - kolommen uit config.R niet in de data: ",
                                 paste(ontbreekt, collapse = ", "))
  onbekend <- setdiff(names(df), c(".rij", names(NUMERIEK), CATEGORISCH))
  if (length(onbekend)) message("LET OP - kolommen niet in config.R (blijven tekst): ",
                                paste(onbekend, collapse = ", "))

  # Numerieke variabelen
  for (k in intersect(names(NUMERIEK), names(df))) {
    ruw <- df[[k]]
    getal <- suppressWarnings(as.numeric(gsub(",", ".", ruw, fixed = TRUE)))
    fout <- is.na(getal) & !is.na(ruw)
    if (sum(!is.na(ruw)) > 0 && mean(fout[!is.na(ruw)]) > 0.5) {
      message("LET OP - '", k, "' is grotendeels geen getal; wordt als categorie behandeld. ",
              "Zet de variabele in config.R bij CATEGORISCH.")
      CATEGORISCH <<- union(CATEGORISCH, k)
      next
    }
    noteer("geen geldig getal", k, df$.rij[fout], ruw[fout], NA)

    grens <- NUMERIEK[[k]]
    buiten <- !is.na(getal) & ((!is.na(grens[1]) & getal < grens[1]) |
                               (!is.na(grens[2]) & getal > grens[2]))
    noteer(sprintf("buiten bereik [%s, %s]", grens[1], grens[2]), k,
           df$.rij[buiten], getal[buiten], NA)
    getal[buiten] <- NA
    df[[k]] <- getal
  }

  # Categorische variabelen
  for (k in intersect(CATEGORISCH, names(df))) {
    ruw <- df[[k]]
    sleutel <- tolower(ruw)
    nieuw <- ruw
    if (!is.null(HERCODERING[[k]])) {
      map <- HERCODERING[[k]]
      names(map) <- tolower(names(map))
      treffer <- sleutel %in% names(map)
      nieuw[treffer] <- map[sleutel[treffer]]
    }
    waarden <- unique(na.omit(tolower(nieuw)))
    is_janee <- length(waarden) > 0 && all(waarden %in% c(JA_CODES, NEE_CODES))
    if (is_janee) {
      nieuw <- ifelse(tolower(nieuw) %in% JA_CODES, "ja", "nee")
      nieuw[is.na(ruw)] <- NA
    }
    veranderd <- !is.na(ruw) & ruw != nieuw
    # ja/nee-hercodering (1 -> ja) is voorspelbaar; alleen echte wijzigingen loggen
    # Alleen tekstcorrecties loggen (bv. "M " -> man); NHR-codes -> labels zijn voorspelbaar
    veranderd <- veranderd & !grepl("^-?[0-9.]+$", ruw)
    if (!is_janee) noteer("hercodering", k, df$.rij[veranderd], ruw[veranderd], nieuw[veranderd])
    df[[k]] <- if (is_janee) factor(nieuw, levels = c("nee", "ja")) else factor(nieuw)
  }

  # Inclusie
  if (ALLEEN_INCLUSIECOHORT && INCLUSIE_KOLOM %in% names(df)) {
    excl <- is.na(df[[INCLUSIE_KOLOM]]) | df[[INCLUSIE_KOLOM]] != "ja"
    noteer("niet in inclusiecohort", INCLUSIE_KOLOM, df$.rij[excl], df[[INCLUSIE_KOLOM]][excl], NA)
    df <- df[!excl, ]
  }

  logboek <- if (length(logboek)) do.call(rbind, logboek) else
    data.frame(stap = character(), kolom = character(), rij_in_bestand = integer(),
               oude_waarde = character(), nieuwe_waarde = character())
  df$.rij <- NULL
  rownames(df) <- NULL
  list(data = df, logboek = logboek)
}

ontbrekende_waarden <- function(df) {
  n <- colSums(is.na(df))
  out <- data.frame(variabele = names(n), ontbrekend_n = n,
                    ontbrekend_pct = round(100 * n / nrow(df), 1))
  out[order(-out$ontbrekend_n), ]
}

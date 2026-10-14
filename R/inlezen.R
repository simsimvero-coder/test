# =============================================================================
# Inlezen en typeconversie (zonder verwijderen of inhoudelijk opschonen)
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

lees_data <- function(pad, blad = NULL) {
  if (!file.exists(pad)) stop("Bestand niet gevonden: ", pad)
  if (!tolower(tools::file_ext(pad)) %in% c("xlsx", "xls"))
    stop("Kies een Excelbestand (.xlsx of .xls): ", pad)

  # Zoek het werkblad en de kopregel met de meeste bekende kolomnamen
  # (een tabblad met toelichting of titelregels boven de tabel wordt zo overgeslagen)
  bekend <- tolower(unique(c(names(NUMERIEK), CATEGORISCH)))
  bladen <- readxl::excel_sheets(pad)
  kandidaten <- if (is.null(blad)) bladen else blad
  beste <- list(score = 0)
  for (b in kandidaten) {
    voorproef <- suppressMessages(readxl::read_excel(pad, sheet = b, col_names = FALSE,
                                                     n_max = 30, col_types = "text"))
    for (r in seq_len(nrow(voorproef))) {
      score <- sum(tolower(trimws(unlist(voorproef[r, ]))) %in% bekend)
      if (score > beste$score) beste <- list(score = score, blad = b, rij = r)
    }
  }
  if (beste$score == 0)
    stop("Geen bekende kolomnamen gevonden. Werkbladen in het bestand: ",
         paste(bladen, collapse = ", "), ". Controleer de kolomnamen in R/config.R.")
  cat(sprintf("   werkblad '%s', kopregel op rij %d (%d van %d bekende kolommen gevonden)\n",
              beste$blad, beste$rij, beste$score, length(bekend)))

  # Alles eerst als tekst lezen; alleen benodigde typen en code-labels worden omgezet
  df <- as.data.frame(readxl::read_excel(pad, sheet = beste$blad, skip = beste$rij - 1,
                                         col_types = "text"))
  attr(df, "kopregel") <- beste$rij

  names(df) <- trimws(names(df))
  # Kolomnamen hoofdletter-ongevoelig koppelen aan de namen uit config.R
  bekend <- unique(c(names(NUMERIEK), CATEGORISCH))
  for (i in seq_along(names(df))) {
    match <- bekend[tolower(bekend) == tolower(names(df)[i])]
    if (length(match) == 1) names(df)[i] <- match
  }
  # Alleen representatie omzetten voor statistische functies. GEEN rijen verwijderen,
  # duplicaten verwijderen, inclusie selecteren, of waarden op grens afkappen.
  # -1/onbekend -> NA conform geregistreerde missingscodering (wordt gerapporteerd).
  for (k in names(df)) {
    x <- trimws(as.character(df[[k]]))
    x[tolower(x) %in% ONTBREKEND_CODES] <- NA_character_
    if (k %in% names(NUMERIEK)) {
      getal <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
      onleesbaar <- !is.na(x) & is.na(getal)
      if (any(onleesbaar)) warning(sprintf("Kolom '%s': %d niet-numerieke waarden als NA geinterpreteerd (geen rijen verwijderd)",
                                            k, sum(onleesbaar)))
      df[[k]] <- getal
    } else if (k %in% CATEGORISCH) {
      if (!is.null(HERCODERING[[k]])) {
        codering <- HERCODERING[[k]]
        hit <- tolower(x) %in% tolower(names(codering))
        x[hit] <- unname(codering[match(tolower(x[hit]), tolower(names(codering)))])
      }
      geldig <- unique(tolower(x[!is.na(x)]))
      if (length(geldig) && all(geldig %in% c(JA_CODES, NEE_CODES))) {
        x <- ifelse(is.na(x), NA_character_, ifelse(tolower(x) %in% JA_CODES, "ja", "nee"))
        df[[k]] <- factor(x, levels = c("nee", "ja"))
      } else df[[k]] <- factor(x)
    } else df[[k]] <- x
  }
  df
}

ontbrekende_waarden <- function(df) {
  n <- colSums(is.na(df))
  out <- data.frame(variabele = names(n), ontbrekend_n = n,
                    ontbrekend_pct = if (nrow(df)) round(100 * n / nrow(df), 1) else NA_real_)
  out[order(-out$ontbrekend_n), ]
}

# Louter signaleren: verandert GEEN enkele rij of waarde.
controle_bron <- function(df) {
  leeg <- if (ncol(df)) rowSums(is.na(df)) == ncol(df) else rep(TRUE, nrow(df))
  sm <- data.frame(controle = c("ingelezen rijen", "exact dubbele rijen", "geheel lege rijen"),
                   aantal = c(nrow(df), sum(duplicated(df)), sum(leeg)))
  numeriek <- do.call(rbind, lapply(intersect(names(NUMERIEK), names(df)), function(k) {
    x <- df[[k]]; gr <- NUMERIEK[[k]]
    fout <- (!is.na(x)) & ((!is.na(gr[1]) & x < gr[1]) |
                             (!is.na(gr[2]) & x > gr[2]))
    data.frame(variabele = k, ontbrekend = sum(is.na(x)),
               buiten_documentatiebereik = sum(fout), rijen_niet_verwijderd = TRUE)
  }))
  list(samenvatting = sm, numeriek = numeriek)
}

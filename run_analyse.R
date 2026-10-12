# =============================================================================
# Analyse aortaboogchirurgie (NHR-dataset)
#
# Gebruik:
#   1. Open aortaboog.Rproj in RStudio (zet de werkmap goed).
#   2. Open dit bestand en klik op 'Source'.
#   3. Het Excelbestand (DATA_BESTAND in R/config.R) wordt op je bureaublad
#      gezocht; staat het ergens anders, dan kies je het in een venster.
# Vanaf de opdrachtregel:  Rscript run_analyse.R pad/naar/dataset.xlsx
# =============================================================================

if (getRversion() < "4.1.0")
  stop("Je R-versie (", getRversion(), ") is te oud. Installeer R 4.1 of nieuwer via https://cran.r-project.org")

# Werkmap = de map van dit script, zodat de bestanden in R/ gevonden worden
# (werkt met source(), de Source-knop in RStudio en Archief > Source in R.app)
script_map <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (!is.null(script_map) && file.exists(file.path(script_map, "R", "config.R"))) setwd(script_map)
if (!file.exists(file.path("R", "config.R")))
  stop("Map R/ niet gevonden. Start het script met source(\"pad/naar/run_analyse.R\") ",
       "of via de Source-knop, niet door de code in de console te plakken.")

benodigd <- c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl")
ontbreekt <- benodigd[!sapply(benodigd, requireNamespace, quietly = TRUE)]
if (length(ontbreekt)) stop("Installeer eerst: install.packages(c(\"",
                            paste(ontbreekt, collapse = "\", \""), "\"), repos = \"https://cloud.r-project.org\")")

source("R/config.R", encoding = "UTF-8")
source("R/inlezen.R", encoding = "UTF-8")
source("R/afgeleid.R", encoding = "UTF-8")
source("R/analyse.R", encoding = "UTF-8")
source("R/figuren.R", encoding = "UTF-8")
source("R/trend.R", encoding = "UTF-8")
source("R/incidentie.R", encoding = "UTF-8")
source("R/periode.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)
gezocht <- if (length(args)) args[1] else DATA_BESTAND
bestand <- zoek_bestand(gezocht)
if (is.null(bestand)) {
  if (!interactive()) stop("Excelbestand niet gevonden: ", gezocht,
                           "\nGebruik: Rscript run_analyse.R pad/naar/dataset.xlsx")
  cat("'", gezocht, "' niet gevonden op het bureaublad. Kies je Excelbestand...\n", sep = "")
  bestand <- file.choose()
}
fig_map <- file.path(UITVOER_MAP, "figuren")
dir.create(fig_map, recursive = TRUE, showWarnings = FALSE)

# 1. Inlezen ------------------------------------------------------------------
cat("1. Inlezen:", bestand, "\n")
ruw <- lees_data(bestand, EXCEL_BLAD)
cat(sprintf("   %d rijen, %d kolommen\n", nrow(ruw), ncol(ruw)))

# 2. Opschonen ----------------------------------------------------------------
cat("2. Opschonen\n")
res <- schoon_op(ruw)
df <- res$data
logboek <- res$logboek
cat(sprintf("   %d patiënten in de analyse, %d regels in het opschoonlogboek\n",
            nrow(df), nrow(logboek)))
indeling <- maak_acuut_electief(df)
df <- maak_descendens(maak_lvef_categorie(indeling$data))
if (nrow(logboek)) {
  telling <- sort(table(paste0(logboek$stap, ifelse(is.na(logboek$kolom), "", paste0(" (", logboek$kolom, ")")))),
                  decreasing = TRUE)
  for (n in names(telling)) cat(sprintf("     - %s: %d\n", n, telling[[n]]))
}
ontbr <- ontbrekende_waarden(df)

cat("   Indeling acuut / electief aneurysma:\n")
print(indeling$overzicht, row.names = FALSE)
df_alle <- df   # alle patiënten, ook niet ingedeeld (voor volume per jaar)
if (ALLEEN_INGEDEELDE_PATIENTEN) {
  niet <- is.na(df$acuut_electief)
  df <- df[!niet, ]
  cat(sprintf("   %d niet-ingedeelde patiënten uitgesloten; %d over voor de analyse\n",
              sum(niet), nrow(df)))
}

if (nrow(df) == 0) stop("Geen patiënten over voor de analyse. Controleer de meldingen hierboven.")

# 3. Analyse ------------------------------------------------------------------
cat("3. Analyse\n")
tabellen <- lapply(TABELLEN, maak_tabel, df = df, groep = GROEP)
oe_groep   <- if (!is.null(GROEP)) oe_ratio(df, GROEP)
oe_centrum <- oe_ratio(df, "centrum")
km <- overleving(df, GROEP)
tr <- trendanalyse(df, GROEP)
inc <- incidentie(df_alle, GROEP)
per <- periodevergelijking(df, GROEP)

excel <- c(
  setNames(tabellen, substr(sub(" - .*", "", names(tabellen)), 1, 31)),
  list("Indeling acuut-electief" = indeling$controle,
       "Volume per jaar" = inc$per_jaar, "Volume per miljoen - RR" = inc$rr,
       "Periode - electief" = per[["electief aneurysma"]], "Periode - acuut" = per[["acuut"]],
       "Periode - totaal" = per[["Totaal"]],
       "O-E per groep" = oe_groep, "O-E per centrum" = oe_centrum,
       "Overleving" = km$tabel,
       "Trend - toetsen" = tr$toetsen, "Trend - aantallen" = tr$aantal,
       "Trend - percentages" = tr$pct, "Trend - temperatuur" = tr$temp, "Ontbrekende waarden" = ontbr, "Opschoonlogboek" = logboek)
)
excel <- excel[!sapply(excel, is.null)]
writexl::write_xlsx(excel, file.path(UITVOER_MAP, "resultaten.xlsx"))
write.csv(df, file.path(UITVOER_MAP, "opgeschoonde_data.csv"), row.names = FALSE, na = "")

# 4. Figuren ------------------------------------------------------------------
cat("4. Figuren\n")
figuren <- c(
  fig_histogrammen(df, fig_map),
  fig_boxplots(df, GROEP, fig_map),
  fig_complicaties(df, GROEP, fig_map),
  fig_kaplan_meier(km, fig_map),
  fig_oe(oe_groep, fig_map),
  fig_oe(oe_centrum, fig_map),
  fig_ontbrekend(ontbr, fig_map),
  trend_figuren(tr, fig_map),
  fig_incidentie(inc, fig_map)
)

# Samenvatting ----------------------------------------------------------------
for (n in names(tabellen)) {
  cat("\n", n, "\n", sep = "")
  print(tabellen[[n]], row.names = FALSE, right = FALSE)
}
if (!is.null(km)) {
  cat("\nOverleving (Kaplan-Meier)", if (km$logrank_p != "") paste0("- log-rank p = ", km$logrank_p), "\n")
  print(km$tabel, row.names = FALSE)
}
if (!is.null(oe_groep)) { cat("\nGeobserveerde/verwachte ziekenhuismortaliteit\n"); print(oe_groep, row.names = FALSE) }

cat("\nVolume per jaar (incl. niet ingedeelde patiënten)\n")
print(inc$per_jaar, row.names = FALSE)
cat("\nVolume per miljoen inwoners - Poisson-regressie (rate ratio per jaar)\n")
print(inc$rr, row.names = FALSE)

for (n in names(per)) {
  cat("\nVergelijking periodes -", n, "(selectie; volledige tabel in resultaten.xlsx)\n")
  print(selecteer_rijen(per[[n]], PERIODE_KERNVARIABELEN), row.names = FALSE, right = FALSE)
}

cat(sprintf("\nTrendanalyse (%d boog-OK's, %s-%s)\n", nrow(tr$data),
            min(tr$data$.jaar), max(tr$data$.jaar)))
print(tr$toetsen[c("maat", "groep", "n", "effectmaat", "schatting", "BI_onder", "BI_boven", "p",
                   "schatting_gecorrigeerd", "p_gecorrigeerd", "p_trendverschil_tussen_groepen")], row.names = FALSE)

cat("\nKlaar. Resultaten in '", UITVOER_MAP, "/':\n", sep = "")
cat("  - resultaten.xlsx (tabellen, O/E, overleving, ontbrekende waarden, opschoonlogboek)\n")
cat("  - opgeschoonde_data.csv\n")
for (f in figuren) cat("  - figuren/", basename(f), "\n", sep = "")

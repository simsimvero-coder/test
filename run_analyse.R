# =============================================================================
# Analyse cardiochirurgische patiëntendataset
#
# Gebruik (vanuit de hoofdmap van het project):
#   Rscript run_analyse.R                         # bestand uit R/config.R
#   Rscript run_analyse.R data/mijn_dataset.xlsx  # ander bestand
# Of in RStudio: open dit bestand en klik op 'Source'.
# =============================================================================

benodigd <- c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl")
ontbreekt <- benodigd[!sapply(benodigd, requireNamespace, quietly = TRUE)]
if (length(ontbreekt)) stop("Installeer eerst: install.packages(c(\"",
                            paste(ontbreekt, collapse = "\", \""), "\"))")

source("R/config.R", encoding = "UTF-8")
source("R/inlezen.R", encoding = "UTF-8")
source("R/afgeleid.R", encoding = "UTF-8")
source("R/analyse.R", encoding = "UTF-8")
source("R/figuren.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)
bestand <- if (length(args)) args[1] else DATA_BESTAND
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
df <- indeling$data
if (nrow(logboek)) {
  telling <- sort(table(paste0(logboek$stap, ifelse(is.na(logboek$kolom), "", paste0(" (", logboek$kolom, ")")))),
                  decreasing = TRUE)
  for (n in names(telling)) cat(sprintf("     - %s: %d\n", n, telling[[n]]))
}
ontbr <- ontbrekende_waarden(df)

cat("   Indeling acuut / electief aneurysma:\n")
print(indeling$overzicht, row.names = FALSE)
if (ALLEEN_INGEDEELDE_PATIENTEN) {
  niet <- is.na(df$acuut_electief)
  df <- df[!niet, ]
  cat(sprintf("   %d niet-ingedeelde patiënten uitgesloten; %d over voor de analyse\n",
              sum(niet), nrow(df)))
}

# 3. Analyse ------------------------------------------------------------------
cat("3. Analyse\n")
tabellen <- lapply(TABELLEN, maak_tabel, df = df, groep = GROEP)
oe_groep   <- if (!is.null(GROEP)) oe_ratio(df, GROEP)
oe_centrum <- oe_ratio(df, "centrum")
km <- overleving(df, GROEP)

excel <- c(
  setNames(tabellen, substr(sub(" - .*", "", names(tabellen)), 1, 31)),
  list("Indeling acuut-electief" = indeling$controle,
       "O-E per groep" = oe_groep, "O-E per centrum" = oe_centrum,
       "Overleving" = km$tabel, "Ontbrekende waarden" = ontbr, "Opschoonlogboek" = logboek)
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
  fig_ontbrekend(ontbr, fig_map)
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

cat("\nKlaar. Resultaten in '", UITVOER_MAP, "/':\n", sep = "")
cat("  - resultaten.xlsx (tabellen, O/E, overleving, ontbrekende waarden, opschoonlogboek)\n")
cat("  - opgeschoonde_data.csv\n")
for (f in figuren) cat("  - figuren/", basename(f), "\n", sep = "")

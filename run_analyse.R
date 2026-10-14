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

# De JIT-compiler van R kan op oudere Macs crashen ("caught segfault"); uitzetten
# maakt het script iets trager maar voorkomt dat.
invisible(compiler::enableJIT(0))
options(width = 250)   # brede tabellen niet afbreken in de console

if (getRversion() < "4.1.0")
  stop("Je R-versie (", getRversion(), ") is te oud. Installeer R 4.1 of nieuwer via https://cran.r-project.org")

# Werkmap = de map van dit script, zodat de bestanden in R/ gevonden worden
# (werkt met source(), de Source-knop in RStudio en Archief > Source in R.app)
script_map <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (!is.null(script_map) && file.exists(file.path(script_map, "R", "config.R"))) setwd(script_map)
if (!file.exists(file.path("R", "config.R")))
  stop("Map R/ niet gevonden. Start het script met source(\"pad/naar/run_analyse.R\") ",
       "of via de Source-knop, niet door de code in de console te plakken.")

benodigd <- c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl", "clubSandwich")
ontbreekt <- benodigd[!sapply(benodigd, requireNamespace, quietly = TRUE)]
if (length(ontbreekt)) stop("Installeer eerst: install.packages(c(\"",
                            paste(ontbreekt, collapse = "\", \""), "\"), repos = \"https://cloud.r-project.org\", type = \"binary\")")

source("R/config.R", encoding = "UTF-8")
source("R/inlezen.R", encoding = "UTF-8")
source("R/afgeleid.R", encoding = "UTF-8")
source("R/analyse.R", encoding = "UTF-8")
source("R/model_hulpfuncties.R", encoding = "UTF-8")
source("R/figuren.R", encoding = "UTF-8")
source("R/trend.R", encoding = "UTF-8")
source("R/incidentie.R", encoding = "UTF-8")
source("R/periode.R", encoding = "UTF-8")
source("R/perfusie.R", encoding = "UTF-8")
source("R/uitkomsten.R", encoding = "UTF-8")
source("R/benchmark.R", encoding = "UTF-8")
source("R/arrest_uitkomst.R", encoding = "UTF-8")

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

# 2. Inhoudelijke selectie, zonder opschoonfunctie of automatische dataverwijdering
cat("2. Selectie en variabelen afleiden (geen opschoning)\n")
df <- ruw
broncontrole <- controle_bron(ruw)
print(broncontrole$samenvatting, row.names = FALSE)
cat(sprintf("   %d originele rijen; dubbele en lege rijen NIET verwijderd\n", nrow(df)))
if (ALLEEN_INCLUSIECOHORT) {
  if (!INCLUSIE_KOLOM %in% names(df))
    stop("De inclusiekolom ontbreekt: ", INCLUSIE_KOLOM,
         ". Zet ALLEEN_INCLUSIECOHORT op FALSE als de bron al uitsluitend het cohort bevat.")
  inbegrepen <- tolower(as.character(df[[INCLUSIE_KOLOM]])) %in% c("ja", "1", "true")
  cat(sprintf("   Expliciete cohortselectie op '%s': %d van %d rijen\n",
              INCLUSIE_KOLOM, sum(inbegrepen), nrow(df)))
  df <- df[inbegrepen, , drop = FALSE]
}
indeling <- maak_acuut_electief(df)
df <- maak_event_free(maak_perfusie(maak_descendens(maak_lvef_categorie(indeling$data))))
ontbr <- ontbrekende_waarden(df)

cat("   Indeling acuut / electief aneurysma:\n")
print(indeling$overzicht, row.names = FALSE)
df_alle <- df  # cohort inclusief niet-ingedeelde patiënten (volume en selectiecontrole)
selectie <- selectie_overzicht(df_alle, GROEP)
df <- df_alle
if (ALLEEN_INGEDEELDE_PATIENTEN) {
  niet <- is.na(df$acuut_electief)
  df <- df[!niet, ]
  cat(sprintf("   %d niet-ingedeelde patiënten uitgesloten; %d over voor de analyse\n",
              sum(niet), nrow(df)))
  if (length(unique(na.omit(as.character(df$acuut_electief)))) != 2)
    stop("Na de indeling zijn niet beide groepen (acuut en electief) aanwezig. Controleer ACUUT_ELECTIEF in R/config.R.")
}

if (nrow(df) == 0) stop("Geen patiënten over voor de analyse. Controleer de meldingen hierboven.")

# 3. Analyse ------------------------------------------------------------------
cat("3. Analyse\n"); flush.console()
tabellen <- lapply(TABELLEN, maak_tabel, df = df, groep = GROEP)
oe_groep   <- if (!is.null(GROEP)) oe_ratio(df, GROEP)
oe_centrum <- oe_ratio(df, "centrum")
km <- overleving(df, GROEP)
cox <- overleving_cox_trend(df, GROEP)
tr <- trendanalyse(df, GROEP)
inc <- incidentie(df_alle, GROEP)
per <- periodevergelijking(df, GROEP)
per_mod <- periode_uitkomstmodellen(df, GROEP)
cat("   Temperatuurmodellen (CR2) - dit kan op een oudere computer 5-15 minuten duren\n")
flush.console()
pa <- perfusie_analyse(df, GROEP)
pd <- perfusie_duur_trends(df, GROEP)
cat("   Uitkomstmodellen (mixed models en CR2) - ook dit kan enkele minuten duren\n"); flush.console()
if (requireNamespace("lme4", quietly = TRUE) &&
    inherits(try(lme4::glmer(cbind(1:3, 3:1) ~ 1 + (1 | g), family = binomial,
                             data = data.frame(g = factor(1:3))), silent = TRUE), "try-error"))
  message("LET OP - lme4 werkt niet goed samen met het geïnstalleerde pakket Matrix. Installeer een nieuwere Matrix:\n",
          "  install.packages(\"Matrix\", repos = \"https://cloud.r-project.org\", type = \"binary\")\n",
          "en start R opnieuw. Tot die tijd wordt MASS::glmmPQL als terugval gebruikt.")
ut <- uitkomst_trends(df, GROEP)
cv <- centrum_variatie(df, GROEP)
los <- opnameduur_analyse(df, GROEP)
rs <- relatieve_overleving(df, GROEP)
smr <- smr_analyse(df, GROEP)
bm <- benchmark_centra(df, GROEP)
au <- arrest_uitkomsten(df, GROEP)
cat("   Mediatieanalyse nierfalen (bootstrap) - kan enkele minuten duren\n"); flush.console()
med <- mediatie_nierfalen(df, GROEP)
cva_tab <- tabel_cva(df, GROEP)

# Gevoeligheidsanalyse: dezelfde patiëntgebonden modellen, gecorrigeerd voor EuroSCORE II
sens_es <- NULL
if (length(GEVOELIGHEID_EUROSCORE_COVARIATEN)) {
  es_jaren <- if (length(GEVOELIGHEID_EUROSCORE_JAREN)) GEVOELIGHEID_EUROSCORE_JAREN else
    sort(unique(suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]])))))
  df_es <- df[suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))) %in% es_jaren, ]
  es_periode <- sprintf("%d-%d", min(es_jaren), max(es_jaren))
  cat("   Gevoeligheidsanalyse: correctie voor", paste(GEVOELIGHEID_EUROSCORE_COVARIATEN, collapse = ", "),
      "- jaren", es_periode, "(n =", nrow(df_es), ")\n")
  flush.console()
  # Periodevergelijking 2013-2018 vs 2019-2024 kan niet binnen deze jaren; wordt overgeslagen.
  sens_es <- met_covariaten(GEVOELIGHEID_EUROSCORE_COVARIATEN, list(
    "Uitkomst - trends" = uitkomst_trends(df_es, GROEP),
    "Trend - toetsen" = trendanalyse(df_es, GROEP)$toetsen,
    "Overleving - Cox per groep" = overleving_cox_trend(df_es, GROEP)))
  sens_es <- c(list("Toelichting" = data.frame(toelichting = c(
    paste0("Gevoeligheidsanalyse: in de kolommen 'gecorrigeerd' is gecorrigeerd voor ",
           paste(GEVOELIGHEID_EUROSCORE_COVARIATEN, collapse = ", "), " in plaats van ",
           paste(UITKOMST_COVARIATEN, collapse = ", "), "."),
    paste0("Alleen operaties uit ", es_periode, " (EuroSCORE II vrijwel compleet); n = ", nrow(df_es), "."),
    "De ongecorrigeerde kolommen gaan over dezelfde jaren (vergelijk die met de gecorrigeerde kolommen)."))),
    Filter(Negate(is.null), sens_es))
  writexl::write_xlsx(sens_es, file.path(UITVOER_MAP, "gevoeligheid_euroscore.xlsx"))
}

excel <- c(
  setNames(tabellen, substr(sub(" - .*", "", names(tabellen)), 1, 31)),
  list("Indeling acuut-electief" = indeling$controle, "Indeling - redenen" = indeling$overzicht,
       "Selectie per jaar" = selectie$jaar, "Selectie per centrum" = selectie$centrum,
       "Datacontrole bron" = broncontrole$samenvatting,
       "Datacontrole numeriek" = broncontrole$numeriek,
       "Volume per jaar" = inc$per_jaar, "Volume per miljoen - RR" = inc$rr,
       "Periode - electief" = per[["electief aneurysma"]], "Periode - acuut" = per[["acuut"]],
       "Periode - totaal" = per[["Totaal"]],
       "Periode uitkomstmodellen" = per_mod, "Temp - modellen" = pa$modellen, "Perfusieduur - trends" = pd, "Temp - per centrum" = pa$centra,
       "Temp - mediaan centrum-jaar" = pa$per_centrum_jaar,
       "Uitkomst - trends" = ut, "Centra - orgaanbescherming" = cv$per_centrum,
       "Centra - bereik" = cv$bereik, "Centra - ICC" = cv$icc,
       "Temp zonder vs met arrest" = cv$zonder_met, "Opnameduur" = los,
       "Relatieve overleving" = rs, "SMR" = smr, "Benchmark - funnel" = bm, "Arrest - uitkomsten" = au, "Arrest - mediatie nierfalen" = med,
       "Temp - categorie alle" = pa$categorieen[["alle procedures"]],
       "Temp - categorie met arrest" = pa$categorieen[["met circulatiestilstand"]], "Temp - ontbrekend per jaar" = pa$ontbrekend_jaar,
       "Temp - ontbrekend per centrum" = pa$ontbrekend_centrum,
       "Temp - % compleet centrum-jaar" = pa$compleet_centrum_jaar,
       "O-E per groep" = oe_groep, "O-E per centrum" = oe_centrum,
       "Overleving" = km$tabel,
       "Trend - toetsen" = tr$toetsen, "Trend - aantallen" = tr$aantal,
       "Trend - percentages" = tr$pct, "Trend - temperatuur" = tr$temp,
       "Overleving - Cox per groep" = cox, "Ontbrekende waarden" = ontbr)
)
excel <- c(excel, setNames(cva_tab, substr(names(cva_tab), 1, 31)))
excel <- excel[!sapply(excel, is.null)]
writexl::write_xlsx(excel, file.path(UITVOER_MAP, "resultaten.xlsx"))
write.csv(df, file.path(UITVOER_MAP, "analyse_data.csv"), row.names = FALSE, na = "")

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
  fig_incidentie(inc, fig_map),
  fig_temperatuur_centra(pa, fig_map),
  fig_temperatuur_per_centrum_jaar(pa, fig_map),
  fig_funnels(bm, fig_map)
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
if (!is.null(cox)) { cat("\nCox-overlevingstrend (ongecorrigeerd en gecorrigeerd per groep)\n"); print(cox, row.names = FALSE) }
if (!is.null(oe_groep)) { cat("\nGeobserveerde/verwachte ziekenhuismortaliteit\n"); print(oe_groep, row.names = FALSE) }

cat("\nSelectie per jaar (indeling: ", ACUUT_ELECTIEF$voorrang, ")\n", sep = "")
print(selectie$jaar, row.names = FALSE)
cat("\nSelectie per centrum\n")
print(selectie$centrum, row.names = FALSE)

cat("\nVolume per jaar (incl. niet ingedeelde patiënten)\n")
print(inc$per_jaar, row.names = FALSE)
cat("\nVolume per miljoen inwoners - Poisson-regressie (rate ratio per jaar)\n")
print(inc$rr, row.names = FALSE)

for (n in names(per)) {
  cat("\nVergelijking periodes -", n, "(selectie; volledige tabel in resultaten.xlsx)\n")
  print(selecteer_rijen(per[[n]], PERIODE_KERNVARIABELEN), row.names = FALSE, right = FALSE)
}

if (!is.null(per_mod)) { cat("\nPeriode-uitkomstmodellen (per groep, ongecorrigeerd + gecorrigeerd)\n"); print(per_mod, row.names = FALSE) }
cat("\nMinimale temperatuur - ontbrekend per jaar\n")
print(pa$ontbrekend_jaar, row.names = FALSE)
cat("\nMinimale temperatuur - ontbrekend per centrum\n")
print(pa$ontbrekend_centrum, row.names = FALSE)
for (n in names(pa$categorieen)) {
  cat("\nMinimale temperatuur - categorieën per groep -", n, "\n")
  print(pa$categorieen[[n]], row.names = FALSE, right = FALSE)
}
cat("\nMinimale temperatuur - trendmodellen (°C per jaar; volledig in tabblad 'Temp - modellen')\n")
m <- pa$modellen
print(data.frame(populatie = m$populatie, groep = m$groep, correctie = m$correctie, model = m$model, n = m$n,
                 trend = sprintf("%.3f (%.3f tot %.3f)", m$verandering_per_jaar, m$BI95_onder, m$BI95_boven),
                 p = m$p, SD_centra = m$SD_slope_tussen_centra, PI95_centra = m$PI95_trend_centrum,
                 ICC = m$ICC_centrum,
                 p_lineariteit = ifelse(grepl("lineariteit", m$extra_toets), m$p_extra, ""),
                 p_verschil_centra = ifelse(grepl("centra", m$extra_toets), m$p_extra, "")),
      row.names = FALSE, right = FALSE)
cat("\nPerfusieduren - trend per jaar (minuten per jaar; random intercept per centrum en CR2)\n")
print(pd, row.names = FALSE)
cat("\nMinimale temperatuur - trend per centrum\n")
print(pa$centra, row.names = FALSE)

cat("\nVariatie tussen centra in orgaanbescherming - bereik over centra\n")
print(cv$bereik, row.names = FALSE)
cat("\nICC's (aandeel variatie tussen centra; ja/nee-uitkomsten op latente schaal)\n")
print(cv$icc, row.names = FALSE)
cat("\nMinimale temperatuur zonder vs. met circulatiestilstand (verschil binnen centra, CR2)\n")
print(cv$zonder_met, row.names = FALSE)
cat("\nUitkomsttrends - OR per kalenderjaar (", attr(ut, "noot"), ")\n", sep = "")
print(ut, row.names = FALSE)
if (!is.null(sens_es)) {
  cat("\nGevoeligheidsanalyse - uitkomsttrends ", es_periode, ", ongecorrigeerd en gecorrigeerd voor ",
      paste(GEVOELIGHEID_EUROSCORE_COVARIATEN, collapse = ", "), " (volledig in gevoeligheid_euroscore.xlsx)\n", sep = "")
  u <- sens_es[["Uitkomst - trends"]]
  print(u[, intersect(c("groep", "uitkomst", "n_ongecorrigeerd", "OR_ongecorrigeerd", "p_ongecorrigeerd", "n_gecorrigeerd", "events_gecorrigeerd", "OR_gecorrigeerd",
                        "BI95_gecorrigeerd", "p_gecorrigeerd", "q_BH_gecorrigeerd", "OR_CR2_gecorrigeerd",
                        "BI95_CR2_gecorrigeerd", "p_CR2_gecorrigeerd"), names(u))], row.names = FALSE)
}
cat("\nCirculatiestilstand en uitkomsten (exploratief; gecorrigeerd voor ", paste(UITKOMST_COVARIATEN, collapse = ", "), ")\n", sep = "")
print(au[au$correctie != "ongecorrigeerd", c("groep", "blootstelling", "uitkomst", "correctie", "n", "ruw_met_vs_zonder", "effectmaat", "schatting", "BI95", "p", "schatting_CR2", "p_CR2")], row.names = FALSE)
if (!is.null(med)) { cat("\nMediatie: duur circulatiestilstand (+10 min) -> nierfalen -> mortaliteit (risicoverschil in procentpunten)\n"); print(med, row.names = FALSE) }
cat("\nOpnameduur (ratio per jaar = geometrische verandering)\n")
print(los, row.names = FALSE)
if (!is.null(rs)) { cat("\nKaplan-Meier en relatieve overleving (Ederer II)\n"); print(rs, row.names = FALSE) }
if (!is.null(smr)) { cat("\nStandardised mortality ratio t.o.v. de Nederlandse bevolking (CBS-levenstabel)\n"); print(smr, row.names = FALSE) }

cat("\nBenchmark per centrum (funnel plots) - centra buiten de grenzen\n")
print(bm[bm$status != "binnen grenzen", c("groep", "uitkomst", "centrum", "n", "geobserveerd", "verwacht",
                                          "waarde", "eenheid", "grens95", "status")], row.names = FALSE)
cat("Landelijk niveau per uitkomst:\n")
print(unique(bm[, c("groep", "uitkomst", "landelijk")]), row.names = FALSE)

cat(sprintf("\nTrendanalyse (%d boog-OK's, %s-%s)\n", nrow(tr$data),
            min(tr$data$.jaar), max(tr$data$.jaar)))
print(tr$toetsen[c("maat", "groep", "n", "effectmaat", "schatting", "BI_onder", "BI_boven", "p",
                   "schatting_gecorrigeerd", "BI_onder_gecorrigeerd", "BI_boven_gecorrigeerd",
                   "n_gecorrigeerd", "p_gecorrigeerd", "p_trendverschil_tussen_groepen")], row.names = FALSE)

cat("\nKlaar. Resultaten in '", UITVOER_MAP, "/':\n", sep = "")
cat("  - resultaten.xlsx (tabellen, twee modellen per groep, overleving, ontbrekende waarden)\n")
cat("  - analyse_data.csv\n")
if (!is.null(sens_es)) cat("  - gevoeligheid_euroscore.xlsx (modellen gecorrigeerd voor EuroSCORE II)\n")
for (f in figuren) cat("  - figuren/", basename(f), "\n", sep = "")

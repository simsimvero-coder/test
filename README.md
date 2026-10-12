# Analyse aortaboogchirurgie (R, NHR-dataset)

R-code voor de analyse van aortaboogpatiënten uit de Nederlandse Hart Registratie
(NHR), met coderingen volgens het **NHR-handboek Cardiochirurgie, versie 22.2.1**.
De code werkt alleen met je eigen Excelbestand. Er zit geen test- of voorbeelddata in.

## Installatie (eenmalig)

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl"))
```

## Gebruik

1. Dubbelklik op `aortaboog.Rproj` (RStudio opent met de juiste werkmap).
2. Open `run_analyse.R` en klik op **Source**.
3. R zoekt **`AORTABOOG_NL_schoon v3(1).xlsx`** op je bureaublad (ook als dat in
   OneDrive staat). Wordt het niet gevonden, dan verschijnt een venster om het bestand te kiezen.
   Een andere bestandsnaam of een ander pad zet je in `R/config.R` bij `DATA_BESTAND`.

De resultaten komen in de map `output/` in de projectmap.

## Coderingen (NHR-handboek)

- **-1 = onbekend** wordt overal als ontbrekend behandeld.
- 0/1-variabelen worden `nee`/`ja`. Codes worden labels, bv. urgentie
  `10 electief`, `20 urgent`, `30 spoed`, `40 redding`, en pathologie
  `10 aneurysma`, `20 acute dissectie`, `30 chronische dissectie`, `40 infectie`, `91 overig`.
- Bij complicaties met subtypes (`gastroint_compl` 1/10/90, `vasc_compl_opname` 10/20) geldt in
  de trendanalyse en bij FTR **elke code > 0** als complicatie.

## Indeling acuut vs. electief aneurysma (`R/afgeleid.R`)

| Groep | Pathologie (`aorta_asc_path`, `aortaboog_path`, `aorta_desc_path`) | Urgentie (alleen als pathologie 10/20 ontbreekt) |
|---|---|---|
| acuut | 20 (acute dissectie) in minstens één segment | 30 spoed / 40 redding |
| electief aneurysma | 10 (aneurysma), nergens 20 | 10 electief / 20 urgent |

- Pathologie gaat altijd vóór urgentie. Code 20 en 10 in verschillende segmenten → acuut.
- Patiënten met alleen chronische dissectie (30), infectie (40) of overig (91) worden op urgentie ingedeeld.
- Het tabblad *Indeling acuut-electief* in `resultaten.xlsx` toont de kruistabel ter controle.

## Analyses en p-waarden

| Onderdeel | Wat | p-waarden |
|---|---|---|
| Tabel 1–3 | kenmerken vóór, tijdens en na de operatie; gem. (SD) of mediaan [IQR], n (%) | acuut vs. electief: Welch t / Mann-Whitney U / chi-kwadraat / Fisher exact |
| O/E-ratio | geobserveerde ziekenhuismortaliteit / verwachte (EuroSCORE II), per groep en centrum | O/E ≠ 1 (exacte Poisson); verschil tussen groepen gecorrigeerd voor EuroSCORE II |
| Overleving | Kaplan-Meier, 30 d / 120 d / 1 / 3 / 5 jaar | log-rank; verschil per tijdpunt |
| Trend 1 | aantal boog-OK's per jaar; aandeel acuut | Poisson (IRR/jaar); logistisch (OR/jaar) |
| Trend 2 | mortaliteit: ziekenhuis, 30 dagen, 1 jaar | OR/jaar, ook gecorrigeerd voor EuroSCORE II |
| Trend 3 | nierfalen, resp. insufficiëntie, gastro-intestinaal, CVA (met + zonder restletsel) | OR/jaar |
| Trend 4 | failure to rescue: ziekenhuismortaliteit bij ≥ 1 complicatie uit trend 3 | OR/jaar |
| Trend 5 | minimale temperatuur | lineaire regressie (°C/jaar) |
| Trend 6 | gebruik circulatiestilstand | OR/jaar |

Voor alle trends geldt daarnaast:
- een p-waarde voor acuut vs. electief **per jaar**;
- een p-waarde voor **verschil in trend tussen de groepen** (interactie jaar × groep).

Bij 30-dagen- en 1-jaarsmortaliteit tellen patiënten die nog leven maar korter gevolgd zijn
dan die termijn niet mee.

## Uitvoer (`output/`)

- `resultaten.xlsx`: Tabel 1–3, Indeling acuut-electief, O/E per groep en per centrum, Overleving,
  Trend (toetsen, aantallen, percentages, temperatuur), Ontbrekende waarden, Opschoonlogboek
- `opgeschoonde_data.csv`
- `figuren/*.png`: alle figuren tonen de bijbehorende p-waarden

## Bestanden

| Bestand | Inhoud |
|---|---|
| `run_analyse.R` | startscript |
| `R/config.R` | alle instellingen: bestand, variabelen, bereiken, NHR-labels, definities |
| `R/inlezen.R` | Excel inlezen en opschonen |
| `R/afgeleid.R` | indeling acuut / electief aneurysma |
| `R/analyse.R` | tabellen, O/E-ratio, overleving |
| `R/trend.R` | trendanalyse |
| `R/figuren.R` | figuren |

## Privacy

`.gitignore` voorkomt dat Excel-/CSV-bestanden en `output/` in git terechtkomen.

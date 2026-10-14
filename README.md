# Analyse aortaboogchirurgie (R, NHR-dataset)

R-code voor de analyse van aortaboogpatiënten uit de Nederlandse Hart Registratie
(NHR), met coderingen volgens het **NHR-handboek Cardiochirurgie, versie 22.2.1**.
De code werkt alleen met je eigen Excelbestand. Er zit geen test- of voorbeelddata in.

## Installatie (eenmalig)

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl", "clubSandwich", "lme4"),
                 repos = "https://cloud.r-project.org", type = "binary")
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

Standaard (`voorrang = "urgentie"` in `R/config.R`):

| Groep | Regel |
|---|---|
| acuut | pathologie 20 (acute dissectie) in asc/boog/desc, **ongeacht urgentie**; of urgentie 30/40 (ook bij pathologie 10) |
| electief aneurysma | pathologie 10 (zonder 20) met urgentie 10/20 of onbekende urgentie |
| geen pathologie 10/20 (niet ingevuld of andere code) | indeling op urgentie: 10/20 = electief, 30/40 = acuut |
| niet ingedeeld | geen pathologie 10/20 en geen bruikbare urgentie |

- Controleer de tabbladen *Selectie per jaar*, *Selectie per centrum* en *Selectie - definities*
  (groepsgrootte en ruwe mortaliteit onder elke definitie). Het landelijke volume per jaar telt altijd
  alle geregistreerde boogoperaties.
- Met `voorrang = "strikt"` (gevoeligheidsanalyse) moeten urgentie én pathologie de groep bevestigen:
  electief = urgentie 10/20 + pathologie 10 zonder 20, acuut = urgentie 30/40 + pathologie 10/20;
  alle andere patiënten worden uitgesloten.
- Met `voorrang = "pathologie"` bepaalt pathologie 20/10 de groep en telt urgentie alleen als de pathologie onbekend is.
- De basis van de indeling staat per patiënt in `acuut_electief_basis`. De reden van uitsluiting staat in *Indeling - redenen*; het tabblad *Indeling acuut-electief* geeft de kruistabel urgentie × pathologie × groep.

## Volume en incidentie (`R/incidentie.R`)

- **Volume per jaar**: totaal, acuut, electief en **niet ingedeeld** per jaar, inclusief de niet-ingedeelde patiënten.
- **Per miljoen inwoners**: met de CBS-bevolking (`BEVOLKING_1JAN` in `R/config.R`; standaard het gemiddelde van 1 januari van het jaar en het volgende jaar). **Controleer de bevolkingscijfers bij CBS StatLine (tabel 37296ned).**
- **Poisson-regressie** met log(bevolking) als offset: rate ratio per jaar (95%-BI, p) voor het totaal over alle jaren, en per groep vanaf `GROEPEN_VANAF_JAAR` (standaard 2015).
- Uitvoer: tabbladen *Volume per jaar* en *Volume per miljoen - RR*, plus figuur `volume_per_miljoen.png`.

## Vergelijking tussen periodes (`R/periode.R`)

Tabel 1–3 met de periode als groep (standaard 2013–2018 vs. 2019–2024, `PERIODES` in `R/config.R`),
apart voor electief, acuut en totaal, met p-waarden. Tabbladen *Periode - electief / acuut / totaal*.

Extra variabelen in Tabel 2: `klepchir_clean`, `descendens_ingreep` (aorta_desc > 0) en
`descendens_pathologie` (aorta_desc_path > 0).

## Perfusie (`R/perfusie.R`)

- `temp_categorie`: ≤20 °C (diep), 20–28 °C (matig), >28 °C (`TEMP_GRENZEN`).
- `circ_arrest_duur_bij_arrest`: duur bij patiënten met circulatiestilstand (`CIRC_ARREST_KOLOM`) en duur > 0.
- Alles voor **alle procedures** en **alleen met circulatiestilstand**, per groep en totaal:
  - ontbrekende temperatuur per jaar, per centrum en % compleet per centrum × jaar;
  - temperatuurcategorieën per groep met p-waarden;
  - **landelijke trend** (marginaal, CR2 geclusterd op centrum) met een exacte decompositie:
    landelijke trend = bijdrage **binnen** centra (praktijkverandering) + bijdrage **verschuiving tussen**
    centra (volume naar centra met een ander temperatuurniveau);
  - **tussen centra**: ICC (aandeel van de variatie tussen centra) en 95%-voorspellingsinterval voor de
    trend van een centrum;
  - trend (°C per jaar):
    1. random intercept per centrum, met toets op afwijking van lineariteit (natural spline, 3 df);
    2. **random slope** (trend per centrum): SD van de trends tussen centra en LRT of de trend verschilt tussen centra;
    3. centrum-fixed-effects met CR2 (Satterthwaite, `clubSandwich`);
    4. alles apart voor acuut en electief, zonder correctie voor urgentie;
    5. gevoeligheidsanalyse: alleen centrum-jaren met ≥80% complete temperatuur (`TEMP_MIN_COMPLEET`);
  - trend per centrum (eigen regressie met 95%-BI en schatting uit het random-slope-model) + figuur
    `temperatuur_trend_per_centrum.png`, en per centrum de jaarmedianen met eigen trend naast de
    landelijke trend (`temperatuur_per_centrum_per_jaar.png`).
- Tabbladen *Temp - …* in `resultaten.xlsx`.

**Gepoolde analyses**: standaard wordt elke analyse apart gedaan voor acuut en electief, zonder correctie
voor urgentie (`ANALYSE_TOTAAL <- FALSE`). Vergelijkingen tussen de groepen (p per jaar, verschil in trend)
staan bij beide groepen. De kolom "Totaal" in Tabel 1–3 en het totale volume per jaar blijven staan.

## Klinische uitkomsten en variatie tussen centra (`R/uitkomsten.R`)

Apart voor acuut en electief, **zonder correctie voor urgentie of EuroSCORE II**.

- **Uitkomsttrends** (OR per kalenderjaar): ziekenhuis-, 30-dagen- en 1-jaarsmortaliteit (alleen
  patiënten met volledige follow-up), major complicaties (`MAJOR_COMPLICATIES`), secundaire
  morbiditeit (`SECUNDAIRE_UITKOMSTEN`) en failure to rescue (ziekenhuismortaliteit bij ≥ 1 major
  complicatie).
  - primair: logistisch mixed model, random intercept per centrum, gecorrigeerd voor leeftijd (per 10
    jaar), geslacht en eerdere cardiochirurgie (`UITKOMST_COVARIATEN`; `lme4::glmer`, met
    `MASS::glmmPQL` als terugval; de kolom `model_primair` toont welk model gebruikt is);
  - gevoeligheid: centrum-fixed-effects met CR2 (Satterthwaite) (`UITKOMST_CR2`);
  - Benjamini–Hochberg (q-waarden) over de secundaire uitkomsten, per groep.
- **Variatie tussen centra**: per centrum mediane temperatuur, % diepe hypothermie, % circulatiestilstand
  en mediane duur; bereik over centra; ICC's (ja/nee-uitkomsten op de latente schaal, π²/3).
- **Temperatuur zonder vs. met circulatiestilstand**, verschil binnen centra (CR2).
- **Opnameduur**: alle patiënten en levend ontslagen, met trend (geometrische ratio per jaar).
- **Relatieve overleving (Ederer II)** op 1, 3 en 5 jaar, met KM-overleving en 95%-BI.
  Hiervoor is een levenstabel nodig (`LEVENSTABEL_BESTAND`, standaard `levenstabel_NL.csv` op het
  bureaublad): CBS StatLine tabel **37360ned** "Levensverwachting; geslacht, leeftijd", onderwerp
  *Sterftekans*, geslacht mannen en vrouwen, leeftijd 0–99 jaar, perioden 2013–2024 → downloaden als
  CSV. Kolommen voor jaar/perioden, geslacht, leeftijd en sterftekans worden automatisch herkend.

Tabbladen *Uitkomst - trends*, *Centra - …*, *Temp zonder vs met arrest*, *Opnameduur*, *Relatieve overleving*.

## Benchmark per centrum (`R/benchmark.R`)

Funnel plots apart voor acuut en electief:
- ziekenhuis- en 30-dagenmortaliteit als **O/E t.o.v. EuroSCORE II** (indirecte standaardisatie;
  standaard geherkalibreerd naar het landelijke niveau, `FUNNEL_HERKALIBREREN`), grenzen exact Poisson;
- major complicaties en failure to rescue als ruw percentage t.o.v. het landelijke percentage,
  grenzen exact binomiaal.
95%- en 99,8%-controlegrenzen; per centrum de status (binnen/boven/onder). Tabblad *Benchmark - funnel*,
figuren `funnel_acuut.png` en `funnel_electief_aneurysma.png`.

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

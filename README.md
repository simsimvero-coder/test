# Analyse cardiochirurgische patiëntendataset (R)

R-code om de dataset in te lezen, op te schonen, te beschrijven (Tabel 1–3 met
groepsvergelijking), de mortaliteit te vergelijken met EuroSCORE II (O/E-ratio)
en de overleving te analyseren (Kaplan-Meier).

## Installatie

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "survival", "writexl"))
# optioneel, voor SPSS-bestanden (.sav):
install.packages("haven")
```

## Gebruik

1. Zet je dataset in de map `data/` (`.csv`, `.xlsx`, `.xls` of `.sav`).
2. Pas `R/config.R` aan; minimaal `DATA_BESTAND` en `GROEP`.
3. Draai vanuit de projectmap:

```bash
Rscript run_analyse.R                          # of: Rscript run_analyse.R data/mijn_data.xlsx
```

Of open `run_analyse.R` in RStudio (met de projectmap als working directory) en klik op **Source**.

Om de code eerst uit te proberen met **fictieve** data met dezelfde kolommen:

```bash
Rscript R/voorbeelddata.R
Rscript run_analyse.R
```

## Instellingen (`R/config.R`)

| Instelling | Standaard | Betekenis |
|---|---|---|
| `GROEP` | `"acuut_electief"` | kolom voor groepsvergelijking (`NULL` = alleen totaal) |
| `ACUUT_ELECTIEF` | zie hieronder | codes voor de indeling acuut / electief aneurysma |
| `ALLEEN_INGEDEELDE_PATIENTEN` | `TRUE` | niet-ingedeelde patiënten uitsluiten |
| `ALLEEN_INCLUSIECOHORT` | `TRUE` | alleen rijen met `incl_cohort` = ja/1 |
| `NUMERIEK` | leeftijd, BMI, kreatinine, … | numerieke variabelen + plausibel bereik |
| `CATEGORISCH` | overige variabelen | categorisch; ja/nee/1/0 wordt automatisch `nee`/`ja` |
| `HERCODERING` | geslacht | extra hercoderingen per variabele |
| `TABELLEN` | Tabel 1–3 | welke variabelen in welke tabel |
| `OVERLEVING_TIJD/STATUS` | `overltijd_dag` / `mort_status` | voor Kaplan-Meier |

Kolomnamen worden hoofdletter-ongevoelig gekoppeld (`bmi` = `BMI`).
Staat `LVEF` bij jou als categorie (bv. goed/matig/slecht) in de data, dan
herkent de code dat en behandelt hem als categorisch (met een melding).

## Indeling acuut vs. electief aneurysma (`R/afgeleid.R`)

| Groep | Pathologie (`aorta_asc_path`, `aortaboog_path`, `aorta_desc_path`) | Urgentie |
|---|---|---|
| acuut | 20 in minstens één segment | 30 of 40 |
| electief aneurysma | 10 (en nergens 20) | 10 of 20 |

- **Pathologie gaat altijd vóór urgentie.** Een dissectie (20) met urgentie 20 is acuut;
  een aneurysma (10) met urgentie 30 is electief.
- Code 20 in het ene segment en 10 in het andere → acuut.
- Urgentie wordt alleen gebruikt als de pathologie ontbreekt én de patiënt aortachirurgie
  had (`urgentie_als_path_ontbreekt`).
- De basis van de indeling staat in de kolom `acuut_electief_basis`; het tabblad
  *Indeling acuut-electief* in `resultaten.xlsx` geeft de kruistabel urgentie × pathologie × groep
  om de indeling te controleren.

## Wat de code doet

| Stap | Bestand | Inhoud |
|---|---|---|
| Inlezen | `R/inlezen.R` | CSV (`;` of `,`), Excel, SPSS; komma als decimaalteken |
| Opschonen | `R/inlezen.R` | lege/dubbele rijen weg, tekst → getal, waarden buiten bereik → NA, ja/nee-codering gelijktrekken, selectie inclusiecohort. **Elke wijziging** komt in het opschoonlogboek. |
| Tabellen | `R/analyse.R` | gem. (SD) of mediaan [IQR] (Shapiro-Wilk), n (%); p-waarden met Welch t-toets / Mann-Whitney U / Welch ANOVA / Kruskal-Wallis / chi-kwadraat / Fisher exact |
| O/E-ratio | `R/analyse.R` | geobserveerde ziekenhuismortaliteit (`overleden_zh`) / verwachte (gem. EuroSCORE II), met exact 95%-BI, per groep en per centrum |
| Overleving | `R/analyse.R` | Kaplan-Meier, overleving na 30 d, 120 d, 1, 3 en 5 jaar, log-rank toets |
| Figuren | `R/figuren.R` | histogrammen, boxplots per groep, postoperatieve uitkomsten, Kaplan-Meier, O/E-forestplots, ontbrekende waarden |

## Uitvoer (`output/`)

- `resultaten.xlsx`: tabbladen Tabel 1–3, O/E per groep en per centrum, Overleving, Ontbrekende waarden, Opschoonlogboek
- `opgeschoonde_data.csv`
- `figuren/*.png`

## Privacy

Werk alleen met gepseudonimiseerde data en volg het datamanagementbeleid van je
instelling. `.gitignore` houdt `data/` (behalve de fictieve voorbeelddata) en
`output/` buiten git.

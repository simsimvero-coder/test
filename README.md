# Analyse patiëntendataset

Python-code om een patiëntendataset in te lezen, op te schonen, te analyseren
(Tabel 1 met groepsvergelijking) en te visualiseren.

## Installatie

```bash
pip install -r requirements.txt
```

## Snel proberen met fictieve data

```bash
python -m patient_analyse.voorbeelddata   # maakt data/voorbeeld_patienten.csv (fictief!)
python -m patient_analyse.main            # draait de hele analyse
```

## Eigen dataset gebruiken

1. Zet je bestand (`.csv`, `.xlsx` of `.xls`) in de map `data/`.
2. Pas `patient_analyse/config.py` aan:
   - `DATA_BESTAND` – pad naar je bestand
   - `ID_KOLOM` – kolom met (gepseudonimiseerd) patiëntnummer
   - `GROEP_KOLOM` – kolom om groepen te vergelijken, of `None`
   - `NUMERIEK` – numerieke variabelen met plausibel bereik `(min, max)`
   - `CATEGORISCH` – categorische variabelen met optionele hercodering
   - `DATUMS` – datumvariabelen (dag-maand-jaar)
3. Draai `python -m patient_analyse.main`
   (of `python -m patient_analyse.main --bestand data/mijn_data.xlsx --groep geen`).

Kolomnamen worden automatisch omgezet naar kleine letters met `_` in plaats van spaties.

## Wat de code doet

| Stap | Module | Inhoud |
|---|---|---|
| Inlezen | `inlezen.py` | CSV (`;` of `,`) en Excel; komma als decimaalteken wordt herkend |
| Opschonen | `inlezen.py` | lege en dubbele rijen weg, tekst → getal, waarden buiten bereik → ontbrekend, hercoderen van categorieën, datums controleren. **Iedere wijziging** komt in een opschoonlogboek. |
| Analyse | `analyse.py` | Tabel 1: gem. (SD) of mediaan [IQR] (Shapiro-Wilk bepaalt welke), n (%); p-waarden met Welch t-toets / Mann-Whitney U / ANOVA / Kruskal-Wallis / chi-kwadraat / Fisher exact; Spearman-correlaties |
| Figuren | `figuren.py` | histogrammen, boxplots per groep, staafdiagrammen, correlatie-heatmap, ontbrekende waarden |

## Uitvoer (map `output/`)

- `resultaten.xlsx` – tabbladen *Tabel 1*, *Ontbrekende waarden*, *Correlaties*, *Opschoonlogboek*
- `opgeschoonde_data.csv`
- `figuren/*.png`

## Privacy

Werk alleen met gepseudonimiseerde data en volg het datamanagementbeleid van je
instelling. `.gitignore` voorkomt dat bestanden in `data/` (behalve de fictieve
voorbeelddata) en `output/` in git terechtkomen.

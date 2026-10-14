# Wijziging: primaire definitie plus urgentiegebaseerde gevoeligheidsanalyse

## Primair: ongewijzigd

De originele functie `maak_acuut_electief()` en de configuratie van
`ACUUT_ELECTIEF` zijn ten opzichte van de vorige ZIP **exact ongewijzigd**.
De hoofdgroep is dezelfde combinatie van NHR-pathologie- en urgentiecodes
als in het originele project. Dat blijft de basis voor `resultaten.xlsx`.

## Gevoeligheidsanalyse: alleen urgentie

- Urgentie **10**: `electief (10)`.
- Urgentie **20**: `urgent (20)`, als eigen groep; deze wordt dus nooit bij
  electief gevoegd.
- Urgentie **30 en 40**: `spoed/redding (30/40)`.
- Andere/ontbrekende urgentie: geen gevoeligheidsgroep; wel apart geteld.

De gevoeligheidsanalyse gebruikt **niet** de pathologiecodes voor de
indeling. Er worden geen afzonderlijke groepen gemaakt voor chronische
dissectie of overige pathologie. Dat betekent niet dat de oorspronkelijke
patiënten met zulke codes uit het cohort worden verwijderd: dat zou de
primaire definitie veranderen. Zij kunnen daarom op basis van de urgentiecode
wel meetellen in de urgentie-only analyse.

De uitkomsten zijn gescheiden:

- `output/resultaten.xlsx`: primaire definitie plus indelingscontroles.
- `output/gevoeligheid_urgentie.xlsx`: aantallen en tabellen, E-versus-S-
  vergelijking, urgent apart, uitkomsttrends, periodes, Cox- en KM-overleving,
  temperatuurtrends, centrumvariatie en opnameduur.
- `output/figuren_urgentie/kaplan_meier.png`: overleving voor electief
  versus spoed/redding op basis van urgentie.

Bij patiëntgebonden regressies wordt eerst ongecorrigeerd en vervolgens
gecorrigeerd voor leeftijd, geslacht en eerdere cardiochirurgie (`UITKOMST_COVARIATEN`),
net als in de vorige ZIP. Ruwe aantallen en Kaplan-Meier-curves blijven
beschrijvend. Directe E-versus-S-effectmaten beschrijven associaties,
geen causale effecten.

De gevoeligheidsanalyse werkt vanaf het gehele **geïncludeerde cohort**
`df_alle`, vóór de selectie op de primaire indeling. Dit is belangrijk,
want anders zouden patiënten met een bruikbare urgentie maar zonder primaire
classificatie ten onrechte verdwijnen uit de gevoeligheidsgroep.

In de directe vergelijking spoed/redding versus electief wordt EuroSCORE II
niet als covariaat gebruikt: urgentie is een onderdeel van EuroSCORE II.

**Tests:** `tests/smoke_test.R` bevat extra tests met synthetische
urgentiewaarden. De volledige analyse en de smoke test zijn uitgevoerd op
synthetische testdata; validatie op de echte NHR-data volgt bij de eerste run.

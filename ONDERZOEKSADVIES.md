# Analysekeuzes — paper 1 en paper 2

## Paper 1: landelijke trends in Nederlandse aortaboogchirurgie

- Cohort: alle geregistreerde aortaboogoperaties die aan vooraf vastgelegde inclusiecriteria voldoen; rapporteer operatie- én patiëntaantallen en registratievolledigheid per jaar/centrum.
- Beschrijvend: jaarlijks volume, operatietypen, pathologie, urgentie, patiëntkenmerken, orgaanbescherming, temperaturen en klinische uitkomsten. Aantallen en percentages zijn ongewogen; geen patiëntgebonden risicocorrectie voor nationaal volume.
- Regressies: per acuut/electief, telkens (1) zonder covariaten en (2) met leeftijd per 10 jaar, geslacht en eerdere cardiochirurgie. Gevoeligheidsanalyse: (3) gecorrigeerd voor alleen EuroSCORE II (`gevoeligheid_euroscore.xlsx`). Rapporteer `n`, gebeurtenissen, 95%-BI, betrouwbaarheidsdiagnostiek en de gebruikte modelvariant voor beide.
- EuroSCORE II bevat reeds informatie over leeftijd, geslacht, eerdere cardiochirurgie en operatie-urgentie; daarom niet in het primaire model (geen correctie voor urgentie binnen de groep) en in de gevoeligheidsanalyse alleen EuroSCORE II (geen dubbele correctie). Beoordeel correlatie, non-lineariteit, calibration drift en kwetsbaarheid door beperkte aantallen events. Overweeg aanvullend een gevoeligheidsmodel met alleen EuroSCORE II en een met afzonderlijke patiëntfactoren.
- Per centrum: modellen met centrumcluster/centrumintercept, en landelijke versus binnen-centrum temperatuurtrends. Standaardisatie van jaartrends voor case-mix is beschrijvend: niet causaal.
- Meet registratieverschillen en ontbrekende waarden; `n` kan verschillen tussen ruwe en gecorrigeerde modellen (complete cases). Een methodologisch vervolg is multiple imputation indien passend.

## Groepsindeling en gevoeligheidsanalyse

**Primaire indeling** (`voorrang = "urgentie"`; urgentie gaat boven pathologie):

| Situatie | Groep |
|---|---|
| urgentie 30/40 | acuut (altijd) |
| pathologie 20 (acute dissectie) | acuut (altijd, ook bij urgentie 10/20 of onbekend) |
| urgentie 10, geen dissectie | electief aneurysma |
| urgentie 20 + pathologie 10 | electief aneurysma |
| urgentie 20 zonder pathologie 10/20; geen bekende urgentie (geen dissectie) | niet ingedeeld |

**Gevoeligheidsanalyse:** dezelfde modellen gecorrigeerd voor alleen EuroSCORE II,
beperkt tot 2018-2024 omdat EuroSCORE II in 2013-2017 bij 27-82% van de patiënten
ontbreekt (`output/gevoeligheid_euroscore.xlsx`).

Voor publicatie: rapporteer de indeling met een stroomdiagram (tabbladen
*Indeling - redenen*, *Selectie per jaar* en *Selectie per centrum*).

## Overlevingsanalyse

1. Maak operatie de tijd-nul; onderscheid overlijden van laatste bekende follow-updatum. Controleer de volledigheid en eventuele links-censoring/truncation.
2. Primaire paper 1-uitkomsten: ziekenhuismortaliteit, 30-dagen- en 1-jaarsmortaliteit. Houd recente cohorten met onvoldoende volledige observatieduur **niet als levend tot dag 365**; gebruik correcte censurering bij tijd-tot-event-analyses en adequate volledigheidschecks voor binaire uitkomsten.
3. Toon ongewogen Kaplan-Meier per groep met 95%-BI en `at risk`-tabel (bijvoorbeeld 30 dagen, 1 en 3 jaar); rapporteer schattingen niet langer dan de beschikbare follow-up. Evalueer betrouwbaarheid van de staart bij weinig patiënten at risk.
4. Cox-model met operatiejaar als predictor **binnen** iedere klinische groep, ongecorrigeerd en vervolgens gecorrigeerd voor de vier factoren. Voor een *vergelijkbare vroege trend* primair 1-jaars tijdshorizon (`COX_MAX_DAGEN <- 365`); vijfjaarscurves zijn beschrijvend en worden vooral door oudere operatiecohorten ondersteund.
5. Toets proportionele hazards met Schoenfeld-residuen. Bij overtreding: tijdsafhankelijke effecten of restricted mean survival time (RMST) over dezelfde tijdshorizon, in plaats van de HR als constant te interpreteren.
6. Relatieve overleving (Ederer II) als secundaire uitkomst met de CBS-levenstabel (`levenstabel_NL.csv`). De eigen implementatie is gevalideerd tegen `relsurv::rs.surv(method = "ederer2")`: verschil hooguit 0,3 procentpunt op 1, 3 en 5 jaar.

## Paper 2: perfusiestrategieën en uitkomsten

Formuleer dan afzonderlijk de causaliteitsvraag (effect van bijvoorbeeld minimale temperatuur, circulatiestilstand of cerebrale perfusie op 30-dagenmortaliteit/neurologische uitkomsten). Maak vooraf een causale DAG; gebruik alleen preoperatieve confounders in primaire adjustment, en modelleer operatietype/ziekte-ernst, centrum en jaartal op een klinisch verdedigbare manier. Vermijd correctie voor mediatoren zoals postoperatieve complicaties en behandel postoperatieve variabelen niet als baseline. Onderzoek non-lineariteit en overlap, eventueel propensity- of doubly robust methoden. Dat is **niet** hetzelfde als de beschrijvende paper-1-analyse.

## Validatie vóór publicatie

De ZIP bevat geen patiëntdataset. De code is uitgevoerd op synthetische testdata (volledige run en `tests/smoke_test.R`), niet op de echte NHR-data. Voer eerst `source("tests/smoke_test.R")` uit in RStudio, vervolgens het hoofdprogramma, en controleer alle modelconvergentie, missingness, effectmaten, tijdsnul en groepsdefinities. Verifieer NHR-labels en CBS-populaties voor publicatie. Gebruik STROBE/RECORD voor het onderzoeksverslag.

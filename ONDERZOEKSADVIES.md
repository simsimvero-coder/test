# Analysekeuzes — paper 1 en paper 2

## Paper 1: landelijke trends in Nederlandse aortaboogchirurgie

- Cohort: alle geregistreerde aortaboogoperaties die aan vooraf vastgelegde inclusiecriteria voldoen; rapporteer operatie- én patiëntaantallen en registratievolledigheid per jaar/centrum.
- Beschrijvend: jaarlijks volume, operatietypen, pathologie, urgentie, patiëntkenmerken, orgaanbescherming, temperaturen en klinische uitkomsten. Aantallen en percentages zijn ongewogen; geen patiëntgebonden risicocorrectie voor nationaal volume.
- Regressies: per acuut/electief, telkens (1) zonder covariaten en (2) met leeftijd per 10 jaar, geslacht, EuroSCORE II in procentpunten en eerdere cardiochirurgie. Rapporteer `n`, gebeurtenissen, 95%-BI, betrouwbaarheidsdiagnostiek en de gebruikte modelvariant voor beide.
- Deze vier covariaten zijn door de onderzoeker opgegeven; EuroSCORE II bevat reeds informatie over leeftijd, geslacht, eerdere cardiochirurgie en operatie-urgentie. Beoordeel correlatie, non-lineariteit, calibration drift en kwetsbaarheid door beperkte aantallen events. Overweeg aanvullend een gevoeligheidsmodel met alleen EuroSCORE II en een met afzonderlijke patiëntfactoren.
- Per centrum: modellen met centrumcluster/centrumintercept, en landelijke versus binnen-centrum temperatuurtrends. Standaardisatie van jaartrends voor case-mix is beschrijvend: niet causaal.
- Meet registratieverschillen en ontbrekende waarden; `n` kan verschillen tussen ruwe en gecorrigeerde modellen (complete cases). Een methodologisch vervolg is multiple imputation indien passend.

## Definitieve hoofdindeling én geplande gevoeligheidsanalyse

**Primaire analyse (ongewijzigd):** `ACUUT_ELECTIEF` met
`voorrang = "urgentie"` combineert de aortapathologiecodes 10 (aneurysma)
en 20 (acute dissectie) met urgentie 10/20/30/40. De oorspronkelijke
prioriteits- en fallbackregels blijven **letterlijk ongewijzigd**.
De hoofdresultaten blijven dus `acuut` versus `electief aneurysma`.

**Gevoeligheidsanalyse (urgentie uitsluitend):**

| Urgentie | Groep | Vergelijking |
|---|---|---|
| 10 | electief (10) | hoofdpaar E versus S |
| 20 | urgent (20) | **apart**; niet samengevoegd met electief |
| 30/40 | spoed/redding (30/40) | hoofdpaar E versus S |
| ontbrekend/anders | geen groep | wel opgenomen in het volledigheidsoverzicht |

Deze indeling gebruikt **geen pathologiecode** om een urgentiegroep toe te wijzen.
Er worden **geen afzonderlijke subgroepen** aangemaakt voor chronische dissectie
of overige pathologie. De oorspronkelijke pathologiegegevens worden niet
verwijderd: zij blijven onderdeel van de bron en van de ongewijzigde primaire
definitie. De gevoeligheidsanalyse is dus geen restrictie op pathologie.

De gevoeligheidsanalyse begint bij het cohort vóór uitsluiting van patiënten
zonder primaire groep; aantallen kunnen daarom verschillen. In
`output/gevoeligheid_urgentie.xlsx` staan de drie groepen, de directe E-versus-S-
uitkomstmodellen, afzonderlijke ongecorrigeerde/gecorrigeerde modellen per groep,
per-jaartrends, KM/Cox, temperatuur/perfusie en het urgente subcohort. In
`output/resultaten.xlsx` staan de primaire analyse en controles.

Voor publicatie: rapporteer beide definities met een stroomdiagram, en maak
helder dat de directe urgentievergelijking *associatief* is en geen causaal
effect van urgentie vaststelt.

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

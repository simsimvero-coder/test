# Synthetische mini-test (bevat GEEN patientdata). Start vanuit de projectmap:
# source("tests/smoke_test.R")
for (p in c("R/config.R", "R/inlezen.R", "R/afgeleid.R", "R/analyse.R",
            "R/model_hulpfuncties.R", "R/trend.R", "R/gevoeligheid_urgentie.R")) source(p, encoding = "UTF-8")
stopifnot(is.function(overleving_cox_trend), is.function(trendanalyse))
set.seed(123)
n <- 1200L
D <- data.frame(
  acuut_electief = factor(rep(c("acuut", "electief aneurysma"), each = n / 2)),
  jaar = factor(sample(2013:2024, n, replace = TRUE)),
  centrum = factor(sample(1:12, n, replace = TRUE)),
  leeftijd = rnorm(n, 66, 11), geslacht = factor(sample(c("man", "vrouw"), n, TRUE)),
  EuroscoreII = pmin(90, pmax(0.5, rgamma(n, 2, 0.2))),
  cardiochir_eerder = factor(sample(c("ja", "nee"), n, TRUE)),
  min_temp = rnorm(n, 23, 4), circ_arrest = factor(sample(c("ja", "nee"), n, TRUE)),
  nierfalen = factor(sample(c("ja", "nee"), n, TRUE, prob = c(.12, .88))),
  resp_insuff = factor(sample(c("ja", "nee"), n, TRUE, prob = c(.18, .82))),
  gastroint_compl = factor(sample(c("ja", "nee"), n, TRUE, prob = c(.10, .90))),
  CVA_restletsel = factor(sample(c("ja", "nee"), n, TRUE, prob = c(.06, .94))),
  CVA_zonder_restletsel = factor(sample(c("ja", "nee"), n, TRUE, prob = c(.07, .93)))
)
p <- plogis(-2.5 + .025 * (D$leeftijd - 65) +
              .12 * (as.integer(as.character(D$jaar)) - 2013))
D$mort_status <- factor(ifelse(runif(n) < p, "ja", "nee"), levels = c("nee", "ja"))
D$overltijd_dag <- ifelse(D$mort_status == "ja", runif(n, 2, 500), runif(n, 370, 2000))
D$overleden_zh <- factor(ifelse(D$mort_status == "ja" & D$overltijd_dag < 30, "ja", "nee"))
D$overl30d <- D$overleden_zh
D$overl1j <- factor(ifelse(D$mort_status == "ja" & D$overltijd_dag <= 365, "ja", "nee"))
D$overl120d <- factor(ifelse(D$mort_status == "ja" & D$overltijd_dag <= 120, "ja", "nee"))
tr <- trendanalyse(D, "acuut_electief")
stopifnot(all(c("schatting", "schatting_gecorrigeerd", "n_gecorrigeerd") %in% names(tr$toetsen)))
km <- overleving(D, "acuut_electief")
stopifnot(!is.null(km), max(km$tabel$dag) <= OVERLEVING_MAX_DAGEN)
cox <- overleving_cox_trend(D, "acuut_electief")
stopifnot(nrow(cox) == 4L, all(c("ongecorrigeerd", "gecorrigeerd") %in% cox$correctie))
cat("Synthetische mini-test afgerond. Geen echte data gebruikt.\n")

# Controleer exacte betekenis van primaire versus alternatieve definitie.
# Ongeacht pathologie wordt urgentie 20 een zelfstandige sensibiliteitsgroep;
# primaire combinatie urgentie/pathologie blijft ongewijzigd.
klein <- data.frame(
  urgentie = factor(c("10 electief", "20 urgent", "30 spoed", "40 redding",
                      "-1", "90")),
  aorta_asc_path = factor(c("10 aneurysma", "10 aneurysma", "10 aneurysma",
                            "20 acute dissectie", "20 acute dissectie", "30 chronische dissectie")),
  aortaboog_path = factor(rep("0 geen", 6)),
  aorta_desc_path = factor(rep("0 geen", 6)),
  jaar = factor(rep(2020, 6)),
  overleden_zh = factor(rep("nee", 6)),
  overl30d = factor(rep("nee", 6))
)
stopifnot(identical(ACUUT_ELECTIEF$voorrang, "urgentie"))
prim <- maak_acuut_electief(klein)$data
sens <- maak_urgentie_sens(prim)
stopifnot(identical(as.character(sens$acuut_electief[2]), "electief aneurysma"))
stopifnot(identical(as.character(sens$acuut_electief[3]), "acuut"))
stopifnot(identical(as.character(sens$acuut_electief[5]), "acuut"))
stopifnot(identical(as.character(sens$urgentie_sens[1:4]),
                    c("electief (10)", "urgent (20)", "acuut (30/40)",
                      "acuut (30/40)")))
stopifnot(all(is.na(sens$urgentie_sens[5:6])))
summ <- overzicht_urgentie_sens(sens)
stopifnot(sum(summ$overzicht$n[summ$overzicht$periode == "alle jaren"]) == 6L)
stopifnot(!any(grepl("chronisch|overig", levels(sens$urgentie_sens), ignore.case = TRUE)))

# Groepsdefinitie voor de echte regressiemodellen met synthetische gegevens.
D$urgentie <- factor(sample(c("10 electief", "20 urgent", "30 spoed", "40 redding"), n, TRUE))
D <- maak_urgentie_sens(D)
stopifnot(nlevels(D$urgentie_sens) == 3L)
stopifnot(sum(D$urgentie_sens == "urgent (20)", na.rm = TRUE) > 0)
stopifnot(sum(as.character(D$urgentie_sens) %in% c("electief (10)", "acuut (30/40)")) +
          sum(D$urgentie_sens == "urgent (20)", na.rm = TRUE) == n)
cat("Urgentie-gebaseerde gevoeligheidsindeling gecontroleerd (synthetische data).\n")

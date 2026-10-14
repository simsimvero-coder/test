# Synthetische mini-test (bevat GEEN patientdata). Start vanuit de projectmap:
# source("tests/smoke_test.R")
for (p in c("R/config.R", "R/inlezen.R", "R/afgeleid.R", "R/analyse.R",
            "R/model_hulpfuncties.R", "R/trend.R")) source(p, encoding = "UTF-8")
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

# Controleer de primaire indeling (urgentie gaat boven pathologie; acute
# dissectie is altijd acuut).
klein <- data.frame(
  urgentie = factor(c("10 electief", "20 urgent", "30 spoed", "40 redding",
                      "-1", "90", "20 urgent", "10 electief", "-1")),
  aorta_asc_path = factor(c("10 aneurysma", "10 aneurysma", "10 aneurysma",
                            "20 acute dissectie", "20 acute dissectie", "30 chronische dissectie",
                            "0 geen", "20 acute dissectie", "10 aneurysma")),
  aortaboog_path = factor(rep("0 geen", 9)),
  aorta_desc_path = factor(rep("0 geen", 9))
)
stopifnot(identical(ACUUT_ELECTIEF$voorrang, "urgentie"))
prim <- as.character(maak_acuut_electief(klein)$data$acuut_electief)
verwacht <- c("electief aneurysma", "electief aneurysma", "acuut", "acuut", "acuut",
              NA, NA, "acuut", NA)
stopifnot(identical(prim, verwacht))
cat("Primaire indeling acuut/electief gecontroleerd (synthetische data).\n")

# Mortaliteit met onvolledige follow-up: het laatste jaar (vrijwel geen
# 1-jaars-follow-up, uitkomst onbekend) mag niet meetellen.
y <- c(rep(TRUE, 2), rep(FALSE, 8), TRUE, rep(NA, 9))
tijd <- c(rep(10, 2), rep(400, 8), 20, rep(100, 9))
jaar <- rep(c(2023, 2024), each = 10)
z <- mortaliteit_fu(y, tijd, jaar, 365)
stopifnot(identical(attr(z, "jaren_uitgesloten"), 2024), all(is.na(z[jaar == 2024])),
          sum(z, na.rm = TRUE) == 2)
cat("Follow-upcorrectie mortaliteit gecontroleerd (synthetische data).\n")

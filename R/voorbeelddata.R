# =============================================================================
# Genereer een FICTIEVE voorbeelddataset met dezelfde kolommen als de echte
# dataset (alleen aortaboogpatiënten), om de code mee te testen. Alle waarden
# zijn willekeurig.
# Bevat bewust 'vuile' data (komma's, tikfouten, uitschieters, dubbele rijen).
#
# Gebruik:  Rscript R/voorbeelddata.R
# =============================================================================

maak_voorbeelddata <- function(n = 1500, seed = 2026) {
  set.seed(seed)
  janee <- function(p) sample(0:1, n, TRUE, c(1 - p, p))
  kies  <- function(x, p = NULL) sample(x, n, TRUE, p)

  jaar <- kies(2015:2023)
  pathologie <- kies(c("aneurysma", "dissectie type A"), c(.55, .45))
  dissectie <- pathologie == "dissectie type A"
  aortachir <- rep(1L, n)
  circ <- as.integer(runif(n) < ifelse(dissectie, .9, .6 + 0.03 * (jaar - 2015)))
  leeftijd <- round(rnorm(n, ifelse(dissectie, 61, 67), 10))
  euro <- round(pmin(80, rlnorm(n, log(ifelse(dissectie, 8, 2.5)), 0.8)), 2)
  overl_zh <- rbinom(n, 1, pmin(0.9, euro / 100 * 1.1))
  ecc_duur <- round(rnorm(n, 110 + 60 * aortachir, 35))
  tijd <- round(ifelse(overl_zh == 1, runif(n, 0, 30), rexp(n, 1 / 4000)))
  tijd <- pmin(tijd, round(runif(n, 300, 3650)))
  status <- as.integer(overl_zh == 1 | (tijd < 3650 & runif(n) < 0.15))

  df <- data.frame(
    centrum = kies(paste("Centrum", LETTERS[1:6])),
    jaar = jaar,
    leeftijd = leeftijd,
    geslacht = kies(c("M", "V"), c(.68, .32)),
    BMI = round(rnorm(n, 27, 4.2), 1),
    kreatinine = round(rlnorm(n, log(90), 0.3)),
    eGFR = round(pmin(120, rnorm(n, 72, 18))),
    LVEF = round(pmin(75, rnorm(n, 55, 10))),
    cardiochir_eerder = janee(.08),
    endocarditis = janee(.04),
    # 10 = electief, 20 = urgent, 30 = spoed, 40 = salvage
    urgentie = ifelse(dissectie, kies(c(20, 30, 40), c(.05, .8, .15)),
                      kies(c(10, 20, 30), c(.7, .25, .05))),
    EuroscoreII = euro,
    ECC = 1L,
    circ_arrest = circ,
    chir_benadering = kies(c("sternotomie", "ministernotomie", "thoracotomie"), c(.85, .1, .05)),
    CABG = janee(.12),
    klepchir = janee(.3),
    aklepchir = janee(.35), mklepchir = janee(.12), pklepchir = janee(.01), tklepchir = janee(.05),
    aortachir = aortachir,
    aotachir_benadering = ifelse(aortachir == 1, kies(c("open", "hybride"), c(.9, .1)), NA),
    aorta_asc = janee(.9),
    # 10 = aneurysma, 20 = dissectie/acuut, 30 = overig
    aorta_asc_path = ifelse(aortachir == 1, ifelse(dissectie, 20, kies(c(10, 30), c(.95, .05))), NA),
    aortaboog_path = ifelse(dissectie, 20, 10),
    aorta_desc = janee(.15),
    aorta_desc_path = NA,
    cardiochir_overig = janee(.1),
    min_temp = round(ifelse(circ == 1, rnorm(n, 22 + 0.5 * (jaar - 2015), 2.5), rnorm(n, 33, 1.5)), 1),
    ECC_duur = ecc_duur,
    AoX_duur = round(ecc_duur * runif(n, .55, .8)),
    circ_arrest_duur = ifelse(circ == 1, round(rnorm(n, 30, 10)), NA),
    overleden_zh = overl_zh,
    ECMO_opname = janee(.02),
    opnameduur = round(rgamma(n, 2.5, 0.3)),
    resp_insuff = janee(.08), beademing = janee(.1), heropname_IC = janee(.06),
    CVA_restletsel = janee(.02), CVA_zonder_restletsel = janee(.02), nierfalen = janee(.05),
    gastroint_compl = janee(.03), vasc_compl_opname = janee(.03),
    mort_status = status,
    overltijd_dag = tijd,
    overl30d = as.integer(status == 1 & tijd <= 30),
    overl120d = as.integer(status == 1 & tijd <= 120),
    overl1j = as.integer(status == 1 & tijd <= 365),
    delier = janee(.12), postop_MI = janee(.02),
    klepchir_clean = NA, aortachir_clean = aortachir, circ_arrest_clean = circ,
    pathologie_groep = pathologie,
    s1_via_bewijs = janee(.5),
    excl_reden = NA,
    incl_cohort = 1L,
    klasse = kies(c("I", "II", "III")),
    stringsAsFactors = FALSE
  )
  df$klepchir_clean <- df$klepchir
  df$aorta_desc_path[df$aorta_desc == 1] <- ifelse(dissectie[df$aorta_desc == 1], 20, 10)
  i <- sample(n, 15); df$aorta_asc_path[i] <- NA; df$aortaboog_path[i] <- NA   # pathologie onbekend

  # Exclusies
  excl <- sample(n, 60)
  df$incl_cohort[excl] <- 0L
  df$excl_reden[excl] <- sample(c("leeftijd < 18", "dubbele registratie", "geen toestemming"), 60, TRUE)

  # Bewust 'vuile' data
  i <- sample(n, 80); df$BMI[i] <- sub(".", ",", df$BMI[i], fixed = TRUE)   # komma als decimaal
  df$geslacht[sample(n, 10)] <- sample(c("man", "vrouw", "m ", "F", "onbekend"), 10, TRUE)
  df$leeftijd[sample(n, 3)] <- c(4, 132, 670)                               # onmogelijke leeftijden
  df$min_temp[sample(n, 2)] <- c(3.2, 330)                                  # invoerfouten
  df$kreatinine[sample(n, 2)] <- c("hemolytisch", "<20")                    # tekst in getal
  df$overleden_zh[sample(n, 5)] <- c("ja", "nee", "Ja", "NEE", "")           # gemengde codering
  for (k in c("BMI", "kreatinine", "eGFR", "LVEF", "EuroscoreII")) df[[k]][sample(n, 40)] <- NA
  rbind(df, df[sample(n, 5), ])                                             # dubbele rijen
}

dir.create("data", showWarnings = FALSE)
write.table(maak_voorbeelddata(), "data/voorbeeld_cardiochirurgie.csv", sep = ";",
            row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("Fictieve voorbeelddataset opgeslagen in data/voorbeeld_cardiochirurgie.csv\n")

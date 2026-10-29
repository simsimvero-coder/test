# =============================================================================
# tests/make_synthetic_data.R
#
# !!! SYNTHETIC DATA — FOR CODE TESTING ONLY — NOT REAL PATIENTS !!!
#
# Generates a fake extract with the same 60 column names and NHR codings as the
# real data, so the pipeline can be executed and debugged without patient data.
# Effect sizes are arbitrary. Nothing produced from this file is a study result.
# Includes deliberate anomalies (duplicates, -1 codes, out-of-range values,
# structurally missing post-2022 variables) to exercise the validation code.
# =============================================================================

suppressPackageStartupMessages(library(data.table))
set.seed(1)
out <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "tests/synthetic/synthetic_arch.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

n <- 5300
centres <- 101:115
cw <- prop.table(c(14, 12, 11, 10, 9, 8, 8, 7, 6, 5, 5, 4, 4, 3, 2))
years <- 2013:2024
yw <- prop.table(seq(300, 520, length.out = 12) * exp(rnorm(12, 0, 0.05)))

d <- data.table(bronrij = seq_len(n),
                centrum = sample(centres, n, TRUE, cw),
                jaar = sample(years, n, TRUE, yw))
ci <- match(d$centrum, centres)
ycen <- d$jaar - 2018
u_c <- rnorm(15, 0, 0.30)          # centre effect on mortality
t_c <- rnorm(15, 0, 2)             # centre effect on temperature
ts_c <- rnorm(15, 0, 0.25)         # centre-specific temperature slope

d[, leeftijd := pmin(pmax(round(rnorm(n, 63, 12)), 18), 92)]
d[, geslacht := sample(c(1, 2), n, TRUE, c(.66, .34))]
d[, BMI := round(rnorm(n, 26.5, 4.5), 1)]
d[, kreatinine := round(exp(rnorm(n, log(88), 0.3)))]
d[, eGFR := round(pmin(141 * pmin(kreatinine / 88.4 / 0.9, 1)^-0.411 * pmax(kreatinine / 88.4 / 0.9, 1)^-1.209 * 0.993^leeftijd, 150))]
d[, LVEF := sample(c(55, 60, 50, 45, 40, 35, 25, 20), n, TRUE, c(.45, .15, .12, .08, .08, .05, .05, .02))]
d[, cardiochir_eerder := rbinom(n, 1, .15)]
d[, endocarditis := rbinom(n, 1, .03)]

# pathology
ad <- rbinom(n, 1, plogis(-0.4 + 0.03 * ycen))
other <- rbinom(n, 1, .05)
arch_code <- ifelse(ad == 1, 20, ifelse(other == 1, ifelse(d$jaar >= 2022, 91, 90), 10))
chronic <- ad == 0 & other == 0 & d$jaar >= 2022 & runif(n) < .15
arch_code[chronic] <- 30
d[, aortaboog_path := arch_code]
d[, aorta_asc_path := ifelse(runif(n) < .9, arch_code, 0)]
d[, aorta_desc := ifelse(runif(n) < plogis(-1.8 + 0.12 * ycen), 10, 0)]
d[, aorta_desc_path := ifelse(aorta_desc > 0, arch_code, 0)]
d[, aorta_asc := sample(c(20, 30, 10, 0), n, TRUE, c(.35, .5, .05, .1))]
d[, aotachir_benadering := ifelse(aorta_desc > 0 & runif(n) < .8, 30, 10)]
d[, urgentie := ifelse(ad == 1, sample(c(30, 20, 40), n, TRUE, c(.8, .17, .03)),
                       sample(c(10, 20, 30), n, TRUE, c(.85, .1, .05)))]
d[, EuroscoreII := round(100 * plogis(-3.2 + 0.035 * (leeftijd - 63) + 0.9 * (urgentie >= 30) +
                                         0.5 * cardiochir_eerder + rnorm(n, 0, .4)), 2)]
d[, ECC := sample(c(20, 10, 30), n, TRUE, c(.85, .13, .02))]
d[, circ_arrest := rbinom(n, 1, .93)]
d[, circ_arrest_clean := circ_arrest]
d[sample(n, 40), circ_arrest_clean := 1 - circ_arrest_clean]
d[, chir_benadering := sample(c(10, 30, 20), n, TRUE, c(.95, .04, .01))]
d[, CABG := rbinom(n, 1, .1)]
d[, aklepchir := rbinom(n, 1, .25)]
d[, mklepchir := rbinom(n, 1, .03)]
d[, pklepchir := 0]; d[, tklepchir := rbinom(n, 1, .01)]
d[, klepchir := as.integer(aklepchir | mklepchir | tklepchir)]
d[, klepchir_clean := klepchir]
d[, aortachir := 1]; d[, aortachir_clean := 1]
d[, cardiochir_overig := rbinom(n, 1, .05)]

# perfusion: deep hypothermia gradually replaced by moderate
d[, min_temp := round(pmin(pmax(rnorm(n, 22 + t_c[ci] + (0.45 + ts_c[ci]) * ycen + 1.5 * (aorta_desc > 0), 2.5), 14), 36), 1)]
d[circ_arrest_clean == 0, min_temp := round(rnorm(.N, 32, 1.5), 1)]
d[, ECC_duur := round(rnorm(n, 210, 50))]
d[, AoX_duur := round(pmin(ECC_duur - 20, rnorm(n, 130, 35)))]
d[, circ_arrest_duur := ifelse(circ_arrest_clean == 1, round(pmax(rnorm(n, 40 + 1.2 * ycen, 15), 5)), 0)]

# outcomes
lp30 <- -3.0 + 0.04 * (d$leeftijd - 63) + 1.1 * ad + 0.6 * (d$urgentie == 40) +
  0.5 * d$cardiochir_eerder - 0.04 * ycen + u_c[ci] - 0.01 * (d$eGFR - 75)
d30 <- rbinom(n, 1, plogis(lp30))
d[, CVA_restletsel := rbinom(n, 1, plogis(-2.7 + 0.5 * ad + 0.6 * d30))]
d[, CVA_zonder_restletsel := rbinom(n, 1, .04)]
d[, nierfalen := rbinom(n, 1, plogis(-2.3 + 0.6 * ad + 1.2 * d30 - 0.015 * (eGFR - 75)))]
d[, resp_insuff := rbinom(n, 1, plogis(-2.4 + 0.8 * d30))]
d[, beademing := rbinom(n, 1, plogis(-0.9 + 0.9 * ad + 1 * d30 - 0.05 * ycen))]
d[, heropname_IC := rbinom(n, 1, .08)]
d[, gastroint_compl := sample(c(0, 10, 90, 1), n, TRUE, c(.95, .02, .025, .005))]
d[, vasc_compl_opname := sample(c(0, 10, 20), n, TRUE, c(.93, .04, .03))]
d[, ECMO_opname := rbinom(n, 1, .02 + .1 * d30)]
d[, postop_MI := rbinom(n, 1, .02)]
d[, delier := rbinom(n, 1, .15)]
late <- d$jaar < 2022 | (d$jaar == 2022 & runif(n) < .5)
d[late, `:=`(ECMO_opname = NA, postop_MI = NA, delier = NA)]

# survival: early hazard + Gompertz background; censoring at status check mid-2025
op_time <- d$jaar + runif(n)
pot_fu <- (2025.5 - op_time) * 365.25
t_early <- ifelse(d30 == 1, sample(0:30, n, TRUE), Inf)
bg <- rexp(n, 1) / (exp(-9.5 + 0.085 * (d$leeftijd + 5)) * (1 + 0.6 * ad)) * 365.25
t_death <- pmin(t_early, ifelse(d30 == 1, Inf, 31 + bg))
d[, mort_status := as.integer(t_death <= pot_fu)]
d[, overltijd_dag := round(pmin(t_death, pot_fu))]
d[, overl30d := as.integer(mort_status == 1 & overltijd_dag <= 30)]
d[, overl120d := as.integer(mort_status == 1 & overltijd_dag <= 120)]
d[, overl1j := as.integer(mort_status == 1 & overltijd_dag <= 365)]
d[, opnameduur := round(pmax(1, exp(rnorm(n, log(10), .5)) + 3 * ad))]
d[mort_status == 1 & overltijd_dag <= 30, opnameduur := pmin(opnameduur, overltijd_dag)]
d[, overleden_zh := as.integer(mort_status == 1 & overltijd_dag <= opnameduur)]

d[, pathologie_groep := fcase(aortaboog_path == 20, "Acute dissectie", aortaboog_path %in% c(10, 30), "Aneurysma", default = "Overig")]
d[, klasse := sample(c("Hemiarch", "Total arch", "Onbekend"), n, TRUE, c(.6, .35, .05))]
d[, s1_via_bewijs := sample(c("ja", "nee"), n, TRUE, c(.7, .3))]
d[, incl_cohort := 1L]; d[, excl_reden := ""]
ex <- sample(n, 300)
d[ex, incl_cohort := 0L]
d[ex, excl_reden := sample(c("Geen boogingreep", "Leeftijd < 18", "Geisoleerde TEVAR"), .N, TRUE)]

# missingness and unknown codes
mk <- function(v, p, code = -1) { i <- sample(n, round(p * n)); set(d, i, v, code) }
mk("BMI", .03, NA); mk("kreatinine", .03, -1); mk("eGFR", .04, NA); mk("LVEF", .06, -1)
mk("geslacht", .002, -1); mk("urgentie", .005, -1); mk("EuroscoreII", .05, -1)
mk("min_temp", .04, -1); mk("circ_arrest_duur", .03, -1); mk("cardiochir_eerder", .01, -1)
mk("beademing", .01, -1); mk("heropname_IC", .02, -1); mk("opnameduur", .01, NA)
d[sample(which(!late), 20), delier := -1]
d[sample(n, 2), min_temp := 5]            # out of range
d[sample(which(d$circ_arrest_clean == 1), 3), circ_arrest_duur := 170]   # long, flagged
dups <- d[sample(n, 3)]; dups[, bronrij := n + seq_len(.N)]
d <- rbind(d, dups)

setcolorder(d, c("bronrij", "centrum", "jaar", "leeftijd", "geslacht", "BMI", "kreatinine", "eGFR",
  "LVEF", "cardiochir_eerder", "endocarditis", "urgentie", "EuroscoreII", "ECC",
  "circ_arrest", "chir_benadering", "CABG", "klepchir", "aklepchir", "mklepchir",
  "pklepchir", "tklepchir", "aortachir", "aotachir_benadering", "aorta_asc",
  "aorta_asc_path", "aortaboog_path", "aorta_desc", "aorta_desc_path",
  "cardiochir_overig", "min_temp", "ECC_duur", "AoX_duur", "circ_arrest_duur",
  "overleden_zh", "ECMO_opname", "opnameduur", "resp_insuff", "beademing",
  "heropname_IC", "CVA_restletsel", "CVA_zonder_restletsel", "nierfalen",
  "gastroint_compl", "vasc_compl_opname", "mort_status", "overltijd_dag", "overl30d",
  "overl120d", "overl1j", "delier", "postop_MI", "klepchir_clean", "aortachir_clean",
  "circ_arrest_clean", "pathologie_groep", "s1_via_bewijs", "excl_reden",
  "incl_cohort", "klasse"))
fwrite(d, out, sep = ";", dec = ",")
cat("SYNTHETIC test data written to", out, "(", nrow(d), "rows )\n")

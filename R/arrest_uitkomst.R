# =============================================================================
# Circulatiestilstand en uitkomsten (exploratief, associatief)
#
# Apart voor acuut en electief. Twee blootstellingen:
#   1. circulatiestilstand ja vs nee (alle patiënten)
#   2. duur van de circulatiestilstand per 10 minuten (alleen bij stilstand)
# Uitkomsten: nierfalen, GI-complicatie, ziekenhuis-, 30-dagen- en
# 1-jaarsmortaliteit (logistisch, OR) en opnameduur (log-lineair, ratio).
# Modellen: random intercept per centrum, ongecorrigeerd en gecorrigeerd voor
# UITKOMST_COVARIATEN (leeftijd, geslacht, eerdere cardiochirurgie);
# gevoeligheid: centrum fixed effects met CR2.
# LET OP: confounding by indication (operatie-uitgebreidheid, ziekte-ernst) is
# niet uit te sluiten; geen causale interpretatie.
# =============================================================================

arrest_uitkomsten <- function(df, groep = GROEP) {
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE)
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, fu) mortaliteit_fu(code_positief(df[[var]]), tijd, df[[TREND_JAAR]], fu)
  uitkomsten <- list(
    "Nierfalen" = code_positief(df$nierfalen),
    "Gastro-intestinale complicatie" = code_positief(df$gastroint_compl),
    "Ziekenhuismortaliteit" = dood("overleden_zh", 0),
    "30-dagenmortaliteit" = dood("overl30d", 30),
    "1-jaarsmortaliteit" = dood("overl1j", 365),
    "Opnameduur (dagen)" = df$opnameduur)
  basis <- data.frame(groep = as.character(df[[groep]]), .centrum = factor(df[[CENTRUM_KOLOM]]),
                      arrest = as.integer(code_positief(df[[CIRC_ARREST_KOLOM]])),
                      duur10 = df$circ_arrest_duur_bij_arrest / 10)
  basis <- voeg_modelcovariaten_toe(basis, df)
  blootstellingen <- list("Circulatiestilstand ja vs nee" = "arrest",
                          "Duur circulatiestilstand (per 10 min)" = "duur10")
  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (bl in names(blootstellingen)) for (u in names(uitkomsten)) {
    term <- blootstellingen[[bl]]
    for (gecorrigeerd in c(FALSE, TRUE)) {
      d <- basis; d$y <- uitkomsten[[u]]
      d <- d[d$groep %in% g & !is.na(d$y) & !is.na(d[[term]]) & !is.na(d$.centrum), ]
      if (gecorrigeerd) d <- modelgevallen(d, "y")
      d <- droplevels(d)
      continu <- u == "Opnameduur (dagen)"
      if (continu) d$y <- log(d$y + 1)
      est <- c(NA, NA, NA, NA); est2 <- c(NA, NA, NA, NA); model <- "te weinig gegevens"
      genoeg <- nrow(d) >= 30 && nlevels(d$.centrum) >= 3 && length(unique(d[[term]])) > 1 &&
                (continu || (sum(d$y) >= 5 && sum(!d$y) >= 5))
      if (genoeg) {
        f <- cov_formule("y", term, gecorrigeerd)
        if (continu) {
          r <- tryCatch({
            m <- nlme::lme(f, random = ~ 1 | .centrum, data = d, method = "REML")
            t <- summary(m)$tTable[term, ]; model <- "lme"
            c(t["Value"], t["Std.Error"], t["p-value"])
          }, error = function(e) c(NA, NA, NA))
        } else {
          r <- logit_mixed(f, d, term); model <- r$model
          r <- if (!is.null(r$est)) c(r$est, r$se, r$p) else c(NA, NA, NA)
        }
        est <- c(exp(r[1]), exp(r[1] - 1.96 * r[2]), exp(r[1] + 1.96 * r[2]), r[3])
        if (heeft_cs) est2 <- tryCatch({
          f2 <- cov_formule("y", term, gecorrigeerd, ".centrum")
          m2 <- if (continu) lm(f2, data = d) else glm(f2, family = binomial, data = d)
          V <- clubSandwich::vcovCR(m2, cluster = d$.centrum, type = "CR2")
          ci <- suppressWarnings(clubSandwich::conf_int(m2, vcov = V, coefs = term))
          ct <- suppressWarnings(clubSandwich::coef_test(m2, vcov = V, test = "Satterthwaite", coefs = term))
          c(exp(ci$beta), exp(ci$CI_L), exp(ci$CI_U), ct[[grep("^p_", names(ct), value = TRUE)[1]]])
        }, error = function(e) c(NA, NA, NA, NA))
      }
      vergelijk <- if (term == "arrest" && nrow(d) && !continu)
        sprintf("%.1f%% vs %.1f%%", 100 * mean(d$y[d$arrest == 1]), 100 * mean(d$y[d$arrest == 0])) else
        if (term == "arrest" && nrow(d)) sprintf("%.0f vs %.0f dagen (mediaan)", median(exp(d$y[d$arrest == 1]) - 1),
                                                  median(exp(d$y[d$arrest == 0]) - 1)) else ""
      res[[length(res) + 1]] <- data.frame(
        groep = g, blootstelling = bl, uitkomst = u, correctie = if (gecorrigeerd) "gecorrigeerd" else "ongecorrigeerd",
        n = nrow(d), ruw_met_vs_zonder = vergelijk, effectmaat = if (continu) "ratio (geometrisch)" else "OR",
        schatting = round(est[1], 3), BI95 = sprintf("%.3f-%.3f", est[2], est[3]), p = fmt_p(est[4]),
        schatting_CR2 = round(est2[1], 3), BI95_CR2 = sprintf("%.3f-%.3f", est2[2], est2[3]), p_CR2 = fmt_p(est2[4]),
        model = model)
    }
  }
  out <- do.call(rbind, res); rownames(out) <- NULL
  out
}

# =============================================================================
# Circulatiestilstand en uitkomsten (exploratief, associatief)
#
# Apart voor acuut en electief. Twee blootstellingen:
#   1. circulatiestilstand ja vs nee (alle patiënten)
#   2. duur van de circulatiestilstand per 10 minuten (alleen bij stilstand)
# Uitkomsten: nierfalen, GI-complicatie, CVA met restletsel, CVA (met of zonder
# restletsel), delier (alleen vanaf DELIER_VANAF_JAAR, eerder niet geregistreerd),
# ziekenhuis-, 30-dagen- en 1-jaarsmortaliteit (logistisch, OR) en opnameduur
# (log-lineair, ratio).
# Modellen: random intercept per centrum; ongecorrigeerd, gecorrigeerd voor
# UITKOMST_COVARIATEN (leeftijd, geslacht, eerdere cardiochirurgie) en daarnaast
# ook voor minimale temperatuur; gevoeligheid: centrum fixed effects met CR2.
# Mediatie: in hoeverre loopt het verband tussen duur en mortaliteit via nierfalen?
# LET OP: confounding by indication (operatie-uitgebreidheid, ziekte-ernst) is
# niet uit te sluiten; geen causale interpretatie.
# =============================================================================

arrest_kern <- function(df, groep = GROEP) {
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE)
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, fu) mortaliteit_fu(code_positief(df[[var]]), tijd, df[[TREND_JAAR]], fu)
  uitkomsten <- list(
    "Nierfalen" = code_positief(df$nierfalen),
    "Gastro-intestinale complicatie" = code_positief(df$gastroint_compl),
    "CVA met restletsel" = code_positief(df$CVA_restletsel),
    "CVA (met of zonder restletsel)" = samengesteld(df, c("CVA_restletsel", "CVA_zonder_restletsel")),
    "Delier (vanaf registratie)" = {
      y <- if ("delier" %in% names(df)) code_positief(df$delier) else rep(NA, nrow(df))
      y[!(suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))) >= DELIER_VANAF_JAAR)] <- NA
      y },
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

arrest_uitkomsten <- function(df, groep = GROEP) {
  basis <- arrest_kern(df, groep)
  if (!"min_temp" %in% names(df)) return(basis)
  # Extra model: ook gecorrigeerd voor minimale temperatuur (patiënten met temperatuur)
  temp <- met_covariaten(c(UITKOMST_COVARIATEN, "min_temp"),
                         arrest_kern(df[!is.na(df$min_temp), ], groep))
  temp <- temp[temp$correctie == "gecorrigeerd", ]
  temp$correctie <- "gecorrigeerd + temperatuur"
  out <- rbind(basis, temp)
  out <- out[order(out$groep, match(out$blootstelling, unique(basis$blootstelling)),
                   match(out$uitkomst, unique(basis$uitkomst))), ]
  rownames(out) <- NULL
  out
}

# Mediatie door nierfalen (counterfactual, g-formule), alleen bij patiënten met
# circulatiestilstand. Blootstelling: 10 minuten langere stilstand dan werkelijk.
# Mediator: nierfalen. Uitkomst: mortaliteit. Logistische modellen met centrum
# fixed effects en UITKOMST_COVARIATEN; interactie duur x nierfalen toegestaan.
#   TE  = totaal effect (risicoverschil per 10 min)
#   NDE = direct effect (niet via nierfalen); NIE = indirect effect via nierfalen
#   proportie gemedieerd = NIE / TE. 95%-BI met bootstrap (MEDIATIE_BOOTSTRAP).
# Veronderstelt geen ongemeten confounding van duur-nierfalen-mortaliteit;
# exploratief.
mediatie_nierfalen <- function(df, groep = GROEP, B = MEDIATIE_BOOTSTRAP) {
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, fu) mortaliteit_fu(code_positief(df[[var]]), tijd, df[[TREND_JAAR]], fu)
  uitk <- list("Ziekenhuismortaliteit" = dood("overleden_zh", 0), "30-dagenmortaliteit" = dood("overl30d", 30),
               "1-jaarsmortaliteit" = dood("overl1j", 365))
  basis <- data.frame(groep = as.character(df[[groep]]), .centrum = factor(df[[CENTRUM_KOLOM]]),
                      a = df$circ_arrest_duur_bij_arrest / 10, m = as.numeric(code_positief(df$nierfalen)))
  basis <- voeg_modelcovariaten_toe(basis, df)
  schat <- function(d) {
    fm <- cov_formule("m", "a", TRUE, ".centrum")
    fy <- update(cov_formule("y", "a", TRUE, ".centrum"), . ~ . + m + a:m)
    mm <- suppressWarnings(glm(fm, binomial, d)); my <- suppressWarnings(glm(fy, binomial, d))
    pm <- function(aw) { nd <- d; nd$a <- aw; predict(mm, nd, type = "response") }
    py <- function(aw, mw) { nd <- d; nd$a <- aw; nd$m <- mw; predict(my, nd, type = "response") }
    EY <- function(a, astar) { p <- pm(astar); mean(p * py(a, 1) + (1 - p) * py(a, 0)) }
    a0 <- d$a; a1 <- d$a + 1
    te <- EY(a1, a1) - EY(a0, a0); nde <- EY(a1, a0) - EY(a0, a0); nie <- EY(a1, a1) - EY(a1, a0)
    c(te = te, nde = nde, nie = nie, prop = nie / te)
  }
  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (u in names(uitk)) {
    d <- basis; d$y <- as.numeric(uitk[[u]])
    d <- modelgevallen(d[d$groep %in% g & !is.na(d$a) & !is.na(d$m) & !is.na(d$y) & !is.na(d$.centrum), ], "y")
    d <- droplevels(d)
    if (nrow(d) < 100 || sum(d$y) < 10 || sum(d$m) < 10) next
    p <- tryCatch(schat(d), error = function(e) rep(NA, 4))
    set.seed(2024)
    bs <- t(replicate(B, tryCatch(schat(d[sample(nrow(d), replace = TRUE), ]), error = function(e) rep(NA, 4))))
    ci <- apply(bs, 2, quantile, c(.025, .975), na.rm = TRUE)
    f <- function(i, schaal = 100, dec = 2) sprintf(paste0("%.", dec, "f (%.", dec, "f tot %.", dec, "f)"),
                                                   schaal * p[i], schaal * ci[1, i], schaal * ci[2, i])
    res[[length(res) + 1]] <- data.frame(
      groep = g, uitkomst = u, n = nrow(d), nierfalen = sum(d$m), overleden = sum(d$y),
      totaal_effect_pp = f(1), direct_effect_pp = f(2), indirect_via_nierfalen_pp = f(3),
      proportie_gemedieerd_pct = f(4, 100, 0), bootstrap = B)
  }
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res); rownames(out) <- NULL
  out
}

# Kenmerken van patiënten met en zonder CVA met restletsel, apart per groep
# (Tabel 1- en 2-variabelen plus perfusie), met p-waarden.
tabel_cva <- function(df, groep = GROEP, cva = "CVA_restletsel") {
  if (!cva %in% names(df)) return(NULL)
  vars <- setdiff(unique(c(TABELLEN[[1]], TABELLEN[[2]], "overleden_zh")), c("centrum", "jaar"))
  vars <- c(vars, "jaar")
  out <- lapply(split(df, as.character(df[[groep]])), function(x) {
    x$CVA <- factor(ifelse(code_positief(x[[cva]]), "CVA met restletsel", "geen CVA"),
                    levels = c("geen CVA", "CVA met restletsel"))
    maak_tabel(x[!is.na(x$CVA), ], intersect(vars, names(x)), "CVA")
  })
  setNames(out, paste("CVA -", names(out)))
}

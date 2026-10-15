# =============================================================================
# Verschil in mortaliteitstrend tussen acuut en electief (revisie)
#
# 1. Formele interactietoets kalenderjaar x presentatie (gepoold model):
#      y ~ jaar_c * acuut + covariaten * acuut + random intercept per centrum x groep
#    De covariaten krijgen een eigen effect per groep en ieder centrum een eigen
#    intercept per groep, zodat het gepoolde model zo dicht mogelijk bij de
#    gestratificeerde primaire modellen blijft. Rapporteert de OR per jaar per
#    groep, de ratio van de OR's (acuut / electief) met 95%-BI en p-interactie.
#    Gevoeligheid: centrum x groep als fixed effects, CR2 (Satterthwaite) met
#    clustering op centrum.
# 2. Absolute, voor case mix gestandaardiseerde mortaliteit per kalenderjaar
#    (marginale standaardisatie / g-computatie): per groep een logistisch model
#    met jaar als factor + UITKOMST_COVARIATEN; voor ieder jaar wordt de kans
#    voorspeld voor ALLE patiënten van die groep alsof ze in dat jaar geopereerd
#    waren, en gemiddeld. Referentiepopulatie = de eigen groep (alle jaren).
# 3. Absolute trend: hetzelfde met een lineaire jaarterm; gestandaardiseerde
#    mortaliteit in het eerste en laatste jaar, de absolute verandering
#    (procentpunt) en het verschil in absolute verandering acuut - electief.
#    95%-BI voor 2 en 3 met bootstrap (patiënten, binnen groep; TRENDVERSCHIL_BOOTSTRAP).
# Mortaliteit met dezelfde follow-upregels als de primaire modellen (mortaliteit_fu).
# =============================================================================

mortaliteit_uitkomsten <- function(df) {
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, min_fu) {
    if (!var %in% names(df)) return(NULL)
    as.vector(mortaliteit_fu(code_positief(df[[var]]), tijd, df[[TREND_JAAR]], min_fu))
  }
  Filter(Negate(is.null), list("Ziekenhuismortaliteit" = dood("overleden_zh", 0),
                               "30-dagenmortaliteit" = dood("overl30d", 30),
                               "1-jaarsmortaliteit" = dood("overl1j", 365)))
}

trendverschil_basis <- function(df, groep) {
  basis <- data.frame(jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                      .centrum = factor(df[[CENTRUM_KOLOM]]),
                      groep = as.character(df[[groep]]))
  basis <- voeg_modelcovariaten_toe(basis, df)
  basis$jaar_c <- basis$jaar - round(mean(basis$jaar, na.rm = TRUE))
  basis$acuut <- as.numeric(basis$groep == "acuut")
  basis
}

fmt_bi <- function(est, lo, hi, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f tot %.", d, "f)"), est, lo, hi)

# ── 1. Interactietoets ──────────────────────────────────────────────────────
mortaliteit_interactie <- function(df, groep = GROEP) {
  if (is.null(groep) || !all(c("acuut", "electief aneurysma") %in% df[[groep]])) return(NULL)
  basis <- trendverschil_basis(df, groep)
  cov <- MODEL_TERMS
  rechts <- c("jaar_c", "acuut", "jaar_c:acuut", cov, if (length(cov)) paste0(cov, ":acuut"))
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE) && UITKOMST_CR2
  uit <- mortaliteit_uitkomsten(df)
  out <- do.call(rbind, lapply(names(uit), function(u) {
    d <- basis; d$y <- uit[[u]]
    d <- d[d$groep %in% c("acuut", "electief aneurysma") & !is.na(d$y) & !is.na(d$jaar) & !is.na(d$.centrum), ]
    d <- modelgevallen(d, "y")
    d$.centrum_g <- droplevels(interaction(d$.centrum, d$acuut))
    voortgang("interactie jaar x presentatie: %s (n = %d)", u, nrow(d))
    rij <- data.frame(uitkomst = u, n_acuut = sum(d$acuut == 1), events_acuut = sum(d$y[d$acuut == 1]),
                      n_electief = sum(d$acuut == 0), events_electief = sum(d$y[d$acuut == 0]),
                      model = "", OR_jaar_electief = "", OR_jaar_acuut = "",
                      ratio_OR_acuut_vs_electief = "", p_interactie = NA_real_,
                      CR2_ratio_OR = "", CR2_p_interactie = NA_real_, stringsAsFactors = FALSE)
    f <- stats::reformulate(c(rechts, "(1 | .centrum_g)"), response = "y")
    fit <- tryCatch({
      m <- suppressMessages(suppressWarnings(lme4::glmer(f, data = d, family = binomial,
                                                         control = lme4::glmerControl(optimizer = "bobyqa"))))
      list(b = lme4::fixef(m), V = as.matrix(vcov(m)), model = "glmer (Laplace)")
    }, error = function(e) tryCatch({
      m <- suppressMessages(MASS::glmmPQL(stats::reformulate(rechts, response = "y"), random = ~ 1 | .centrum_g,
                                          family = binomial, data = d, verbose = FALSE))
      list(b = nlme::fixef(m), V = as.matrix(vcov(m)), model = "glmmPQL (terugval)")
    }, error = function(e2) list(model = paste("mislukt:", conditionMessage(e2)))))
    rij$model <- fit$model
    if (!is.null(fit$b)) {
      b <- fit$b; V <- fit$V
      lin <- function(w) { est <- sum(w * b[names(w)]); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
      e <- lin(c(jaar_c = 1)); a <- lin(c(jaar_c = 1, "jaar_c:acuut" = 1)); i <- lin(c("jaar_c:acuut" = 1))
      ci <- function(x) fmt_bi(exp(x[1]), exp(x[1] - 1.96 * x[2]), exp(x[1] + 1.96 * x[2]))
      rij$OR_jaar_electief <- ci(e); rij$OR_jaar_acuut <- ci(a); rij$ratio_OR_acuut_vs_electief <- ci(i)
      rij$p_interactie <- 2 * pnorm(-abs(i[1] / i[2]))
    }
    if (heeft_cs && nlevels(d$.centrum) >= 3) tryCatch({
      f_fe <- stats::reformulate(c("0", ".centrum_g", setdiff(rechts, "acuut")), response = "y")
      m <- glm(f_fe, family = binomial, data = d)
      V <- clubSandwich::vcovCR(m, cluster = d$.centrum, type = "CR2")
      ct <- suppressWarnings(clubSandwich::conf_int(m, vcov = V, coefs = "jaar_c:acuut"))
      tt <- suppressWarnings(clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite", coefs = "jaar_c:acuut"))
      rij$CR2_ratio_OR <- fmt_bi(exp(ct$beta), exp(ct$CI_L), exp(ct$CI_U))
      rij$CR2_p_interactie <- tt[[grep("^p_", names(tt), value = TRUE)[1]]]
    }, error = function(e) message("CR2 interactie overgeslagen voor ", u, ": ", conditionMessage(e)))
    rij
  }))
  if (is.null(out)) return(NULL)
  for (k in c("p_interactie", "CR2_p_interactie")) out[[k]] <- vapply(out[[k]], fmt_p, character(1))
  attr(out, "noot") <- paste0("Gepoold model: jaar x presentatie, covariaten (", paste(UITKOMST_COVARIATEN, collapse = ", "),
                              ") per groep, random intercept per centrum x groep. Ratio < 1: sterkere daling bij acuut.")
  out
}

# ── 2 en 3. Gestandaardiseerde mortaliteit per jaar en absolute trend ───────
standaardiseer <- function(d, jaren) {
  cov <- MODEL_TERMS
  m_f <- suppressWarnings(glm(stats::reformulate(c("factor(jaar)", cov), response = "y"), family = binomial, data = d))
  per_jaar <- vapply(jaren, function(j) { nd <- d; nd$jaar <- j; mean(predict(m_f, nd, type = "response")) }, numeric(1))
  m_l <- suppressWarnings(glm(stats::reformulate(c("jaar", cov), response = "y"), family = binomial, data = d))
  rand <- vapply(range(jaren), function(j) { nd <- d; nd$jaar <- j; mean(predict(m_l, nd, type = "response")) }, numeric(1))
  c(per_jaar, rand)
}

mortaliteit_gestandaardiseerd <- function(df, groep = GROEP, B = TRENDVERSCHIL_BOOTSTRAP) {
  if (is.null(groep)) return(NULL)
  basis <- trendverschil_basis(df, groep)
  uit <- mortaliteit_uitkomsten(df)
  per_jaar <- list(); absoluut <- list()
  for (u in names(uit)) {
    d0 <- basis; d0$y <- as.numeric(uit[[u]])
    d0 <- modelgevallen(d0[!is.na(d0$y) & !is.na(d0$jaar) & !is.na(d0$groep), ], "y")
    res <- list()
    for (g in intersect(c("acuut", "electief aneurysma"), unique(d0$groep))) {
      d <- d0[d0$groep == g, ]
      jaren <- sort(unique(d$jaar))
      voortgang("gestandaardiseerde mortaliteit: %s - %s (n = %d, %d bootstraps)", g, u, nrow(d), B)
      est <- standaardiseer(d, jaren)
      set.seed(2024)
      boot <- t(vapply(seq_len(B), function(b) tryCatch(
        standaardiseer(d[sample.int(nrow(d), replace = TRUE), ], jaren),
        error = function(e) rep(NA_real_, length(est))), numeric(length(est))))
      res[[g]] <- list(est = est, boot = boot, jaren = jaren)
      k <- length(jaren)
      lo <- apply(boot[, 1:k, drop = FALSE], 2, quantile, 0.025, na.rm = TRUE)
      hi <- apply(boot[, 1:k, drop = FALSE], 2, quantile, 0.975, na.rm = TRUE)
      per_jaar[[length(per_jaar) + 1]] <- data.frame(
        uitkomst = u, groep = g, jaar = jaren,
        n = as.vector(table(factor(d$jaar, levels = jaren))),
        overleden = as.vector(tapply(d$y, factor(d$jaar, levels = jaren), sum)),
        ruw_pct = round(100 * as.vector(tapply(d$y, factor(d$jaar, levels = jaren), mean)), 1),
        gestandaardiseerd_pct = round(100 * est[1:k], 1),
        BI95_onder = round(100 * lo, 1), BI95_boven = round(100 * hi, 1))
      verschil <- 100 * (boot[, k + 2] - boot[, k + 1])
      absoluut[[length(absoluut) + 1]] <- data.frame(
        uitkomst = u, groep = g, periode = sprintf("%d-%d", min(jaren), max(jaren)), n = nrow(d),
        gestand_pct_eerste_jaar = fmt_bi(100 * est[k + 1], 100 * quantile(boot[, k + 1], .025, na.rm = TRUE),
                                         100 * quantile(boot[, k + 1], .975, na.rm = TRUE), 1),
        gestand_pct_laatste_jaar = fmt_bi(100 * est[k + 2], 100 * quantile(boot[, k + 2], .025, na.rm = TRUE),
                                          100 * quantile(boot[, k + 2], .975, na.rm = TRUE), 1),
        absolute_verandering_pp = fmt_bi(100 * (est[k + 2] - est[k + 1]), quantile(verschil, .025, na.rm = TRUE),
                                         quantile(verschil, .975, na.rm = TRUE), 1),
        verandering_pp_per_jaar = round(100 * (est[k + 2] - est[k + 1]) / (max(jaren) - min(jaren)), 2),
        stringsAsFactors = FALSE)
    }
    if (all(c("acuut", "electief aneurysma") %in% names(res))) {
      a <- res[["acuut"]]; e <- res[["electief aneurysma"]]
      ka <- length(a$jaren); ke <- length(e$jaren)
      dd <- 100 * ((a$est[ka + 2] - a$est[ka + 1]) - (e$est[ke + 2] - e$est[ke + 1]))
      db <- 100 * ((a$boot[, ka + 2] - a$boot[, ka + 1]) - (e$boot[, ke + 2] - e$boot[, ke + 1]))
      absoluut[[length(absoluut) + 1]] <- data.frame(
        uitkomst = u, groep = "verschil acuut - electief", periode = "", n = NA,
        gestand_pct_eerste_jaar = "", gestand_pct_laatste_jaar = "",
        absolute_verandering_pp = fmt_bi(dd, quantile(db, .025, na.rm = TRUE), quantile(db, .975, na.rm = TRUE), 1),
        verandering_pp_per_jaar = NA_real_, stringsAsFactors = FALSE)
    }
  }
  list(per_jaar = do.call(rbind, per_jaar), absoluut = do.call(rbind, absoluut))
}

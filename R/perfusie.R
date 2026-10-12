# =============================================================================
# Perfusie: minimale temperatuur over de jaren, rekening houdend met centrum
#
# Per populatie (alle procedures / alleen met circulatiestilstand) en per groep
# (electief, acuut, totaal):
#   1. Random intercept per centrum (nlme::lme)
#   2. Random slope: trend per centrum (1 + jaar | centrum); SD van de slopes en
#      toets of de trend tussen centra verschilt (LRT, 50:50-mengverdeling)
#   3. Centrum-fixed-effects met CR2 cluster-robuuste SE en Satterthwaite-df
#   4. Lineariteit: lineair vs. natural spline (3 df) voor jaar (LRT, ML)
#   5. Totaal: random slope en CR2 gecorrigeerd voor acuut/electief
#   6. Gevoeligheidsanalyse: alleen centrum-jaren met >= TEMP_MIN_COMPLEET compleet
# Daarnaast: ontbrekende temperatuur per jaar, per centrum en per centrum x jaar,
# temperatuurcategorieën en de trend per centrum (+ figuur).
# =============================================================================

TEMP_MIN_COMPLEET <- 0.80

mengverdeling_p <- function(L) 0.5 * pchisq(L, 1, lower.tail = FALSE) + 0.5 * pchisq(L, 2, lower.tail = FALSE)

perfusie_analyse <- function(df, groep = GROEP, uitkomst = TREND_TEMPERATUUR,
                             jaar_kolom = TREND_JAAR, centrum = CENTRUM_KOLOM) {
  if (!all(c(uitkomst, jaar_kolom, centrum) %in% names(df))) return(NULL)
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE)
  if (!heeft_cs)
    message("LET OP - pakket 'clubSandwich' niet geïnstalleerd; CR2-analyses overgeslagen. Installeer met: ",
            "install.packages(\"clubSandwich\", repos = \"https://cloud.r-project.org\", type = \"binary\")")

  basis <- data.frame(
    y = df[[uitkomst]],
    jaar = suppressWarnings(as.numeric(as.character(df[[jaar_kolom]]))),
    .centrum = factor(df[[centrum]]),
    groep = if (!is.null(groep) && groep %in% names(df)) as.character(df[[groep]]) else "Totaal",
    arrest = if ("circ_arrest_duur_bij_arrest" %in% names(df)) {
      kol <- if (CIRC_ARREST_KOLOM %in% names(df)) CIRC_ARREST_KOLOM else "circ_arrest"
      df[[kol]] %in% "ja"
    } else rep(NA, nrow(df))
  )
  basis$acuut <- as.integer(basis$groep == "acuut")
  basis <- basis[!is.na(basis$jaar) & !is.na(basis$.centrum), ]
  jaar_midden <- round(mean(basis$jaar))
  basis$jaar_c <- basis$jaar - jaar_midden

  populaties <- list("alle procedures" = basis, "met circulatiestilstand" = basis[basis$arrest %in% TRUE, ])

  # ── Ontbrekende temperatuur ────────────────────────────────────────────────
  ontbrekend <- function(d, per) {
    a <- aggregate(cbind(n = 1, ontbrekend = is.na(d$y)) ~ d[[per]], FUN = sum, data = d)
    names(a)[1] <- per
    a$ontbrekend_pct <- round(100 * a$ontbrekend / a$n, 1)
    a
  }
  ontbr_jaar <- do.call(rbind, lapply(names(populaties), function(p)
    cbind(populatie = p, ontbrekend(populaties[[p]], "jaar"))))
  ontbr_centrum <- do.call(rbind, lapply(names(populaties), function(p)
    cbind(populatie = p, ontbrekend(populaties[[p]], ".centrum"))))
  names(ontbr_centrum)[names(ontbr_centrum) == ".centrum"] <- "centrum"
  matrix_compleet <- function(d) {
    t <- tapply(!is.na(d$y), list(centrum = d$.centrum, jaar = d$jaar), mean)
    m <- as.data.frame.matrix(round(100 * t, 0))
    cbind(centrum = rownames(m), m)
  }
  ontbr_matrix <- matrix_compleet(populaties[["met circulatiestilstand"]])

  # ── Modellen ───────────────────────────────────────────────────────────────
  rij <- function(pop, g, model, d, est = NA, lo = NA, hi = NA, p = NA, sd_slope = NA,
                  p_extra = NA, extra_label = "")
    data.frame(populatie = pop, groep = g, model = model, n = nrow(d),
               n_centra = length(unique(d$.centrum)), verandering_per_jaar = round(est, 3),
               BI95_onder = round(lo, 3), BI95_boven = round(hi, 3), p = fmt_p(p),
               SD_slope_tussen_centra = round(sd_slope, 3), extra_toets = extra_label,
               p_extra = fmt_p(p_extra))
  lme_ctrl <- nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)
  t_ci <- function(t) { q <- qt(0.975, t["DF"]); c(t["Value"], t["Value"] - q * t["Std.Error"],
                                                   t["Value"] + q * t["Std.Error"], t["p-value"]) }

  model_ri <- function(d, pop, g) tryCatch({
    m <- nlme::lme(y ~ jaar_c, random = ~ 1 | .centrum, data = d, method = "REML")
    r <- t_ci(summary(m)$tTable["jaar_c", ])
    # Lineariteit: lineair vs. natural spline (ML)
    m_lin <- nlme::lme(y ~ jaar_c, random = ~ 1 | .centrum, data = d, method = "ML")
    m_spl <- nlme::lme(y ~ splines::ns(jaar_c, df = 3), random = ~ 1 | .centrum, data = d, method = "ML")
    p_lin <- anova(m_lin, m_spl)[["p-value"]][2]
    rij(pop, g, "random intercept (centrum)", d, r[1], r[2], r[3], r[4],
        p_extra = p_lin, extra_label = "afwijking van lineariteit (spline 3 df)")
  }, error = function(e) rij(pop, g, paste("random intercept - mislukt:", conditionMessage(e)), d))

  model_rs <- function(d, pop, g, label = "random slope (trend per centrum)", covariaat = FALSE) tryCatch({
    f <- if (covariaat) y ~ jaar_c + acuut else y ~ jaar_c
    m <- nlme::lme(f, random = ~ 1 + jaar_c | .centrum, data = d, method = "REML", control = lme_ctrl)
    m0 <- nlme::lme(f, random = ~ 1 | .centrum, data = d, method = "REML")
    r <- t_ci(summary(m)$tTable["jaar_c", ])
    sd_s <- as.numeric(nlme::VarCorr(m)["jaar_c", "StdDev"])
    L <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(m0)))
    rij(pop, g, label, d, r[1], r[2], r[3], r[4], sd_slope = sd_s,
        p_extra = mengverdeling_p(max(L, 0)), extra_label = "trend verschilt tussen centra (LRT)")
  }, error = function(e) rij(pop, g, paste(label, "- mislukt:", conditionMessage(e)), d))

  model_cr2 <- function(d, pop, g, label = "centrum fixed effects, CR2 (Satterthwaite)", covariaat = FALSE) {
    if (!heeft_cs) return(rij(pop, g, paste(label, "- clubSandwich ontbreekt"), d))
    tryCatch({
      m <- lm(if (covariaat) y ~ jaar + acuut + .centrum else y ~ jaar + .centrum, data = d)
      ci <- clubSandwich::conf_int(m, vcov = "CR2", cluster = d$.centrum, coefs = "jaar")
      ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = d$.centrum, test = "Satterthwaite", coefs = "jaar")
      rij(pop, g, label, d, ci$beta, ci$CI_L, ci$CI_U, ct[[grep("^p_", names(ct), value = TRUE)[1]]])
    }, error = function(e) rij(pop, g, paste(label, "- mislukt:", conditionMessage(e)), d))
  }

  modellen <- list(); centra <- list()
  for (pop in names(populaties)) {
    dp <- populaties[[pop]]
    # centrum-jaren met voldoende complete temperatuur (gevoeligheidsanalyse)
    compleet <- tapply(!is.na(dp$y), paste(dp$.centrum, dp$jaar), mean)
    goed <- names(compleet)[compleet >= TEMP_MIN_COMPLEET]
    subsets <- c(split(dp, dp$groep), list(Totaal = dp))
    for (g in names(subsets)) {
      d <- droplevels(subsets[[g]][!is.na(subsets[[g]]$y), ])
      if (nrow(d) < 20 || nlevels(d$.centrum) < 3) next
      d_s <- droplevels(d[paste(d$.centrum, d$jaar) %in% goed, ])
      modellen[[length(modellen) + 1]] <- rbind(
        model_ri(d, pop, g), model_rs(d, pop, g), model_cr2(d, pop, g),
        if (g == "Totaal" && length(unique(d$acuut)) > 1) rbind(
          model_rs(d, pop, g, "random slope, gecorrigeerd voor acuut/electief", covariaat = TRUE),
          model_cr2(d, pop, g, "CR2, gecorrigeerd voor acuut/electief", covariaat = TRUE)),
        model_rs(d_s, pop, g, sprintf("gevoeligheid: centrum-jaren >=%d%% compleet - random slope", 100 * TEMP_MIN_COMPLEET)),
        model_cr2(d_s, pop, g, sprintf("gevoeligheid: centrum-jaren >=%d%% compleet - CR2", 100 * TEMP_MIN_COMPLEET))
      )
      # Trend per centrum: eigen lineaire regressie (met 95%-BI) en BLUP uit random-slope-model
      if (g == "Totaal") {
        blup <- tryCatch({
          m <- nlme::lme(y ~ jaar_c, random = ~ 1 + jaar_c | .centrum, data = d, method = "REML", control = lme_ctrl)
          setNames(coef(m)[, "jaar_c"], rownames(coef(m)))
        }, error = function(e) NULL)
        centra[[pop]] <- do.call(rbind, lapply(levels(d$.centrum), function(c) {
          x <- d[d$.centrum == c, ]
          ols <- if (length(unique(x$jaar)) >= 3) {
            s <- coef(summary(lm(y ~ jaar, data = x)))["jaar", ]
            q <- qt(0.975, nrow(x) - 2); c(s[1], s[1] - q * s[2], s[1] + q * s[2], s[4])
          } else rep(NA, 4)
          data.frame(populatie = pop, centrum = c, n = nrow(x),
                     jaren = sprintf("%d-%d", min(x$jaar), max(x$jaar)),
                     mediaan_temp = round(median(x$y), 1),
                     trend_eigen = round(ols[1], 3), BI95_onder = round(ols[2], 3),
                     BI95_boven = round(ols[3], 3), p = fmt_p(ols[4]),
                     trend_random_slope_BLUP = if (is.null(blup)) NA else round(unname(blup[c]), 3))
        }))
      }
    }
  }

  # Temperatuurcategorieën per groep (met p-waarden), per populatie
  categorieen <- lapply(setNames(names(populaties), names(populaties)), function(pop) {
    idx <- if (pop == "alle procedures") rep(TRUE, nrow(df)) else {
      kol <- if (CIRC_ARREST_KOLOM %in% names(df)) CIRC_ARREST_KOLOM else "circ_arrest"
      df[[kol]] %in% "ja"
    }
    maak_tabel(df[idx, ], c("min_temp", "temp_categorie"), groep)
  })

  list(modellen = do.call(rbind, modellen), centra = do.call(rbind, centra),
       ontbrekend_jaar = ontbr_jaar, ontbrekend_centrum = ontbr_centrum,
       compleet_centrum_jaar = ontbr_matrix, categorieen = categorieen, jaar_midden = jaar_midden)
}

fig_temperatuur_centra <- function(pa, map) {
  if (is.null(pa$centra)) return(NULL)
  d <- pa$centra[pa$centra$populatie == "met circulatiestilstand" & !is.na(pa$centra$trend_eigen), ]
  if (!nrow(d)) d <- pa$centra[!is.na(pa$centra$trend_eigen), ]
  if (!nrow(d)) return(NULL)
  rs <- pa$modellen[pa$modellen$populatie == d$populatie[1] & pa$modellen$groep == "Totaal" &
                      pa$modellen$model == "random slope (trend per centrum)", ]
  d$centrum <- factor(d$centrum, levels = d$centrum[order(d$trend_eigen)])
  p <- ggplot(d, aes(trend_eigen, centrum)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    { if (nrow(rs)) geom_vline(xintercept = rs$verandering_per_jaar, colour = "#C0392B") } +
    geom_errorbarh(aes(xmin = BI95_onder, xmax = BI95_boven), height = 0.25, colour = KLEUR) +
    geom_point(aes(size = n), colour = KLEUR) +
    geom_point(aes(x = trend_random_slope_BLUP), shape = 4, size = 2.5, colour = "#C0392B") +
    labs(title = "Trend in minimale temperatuur per centrum",
         subtitle = paste0(d$populatie[1], "; punt + 95%-BI = eigen regressie per centrum, ",
                           "x = schatting random-slope-model\nrode lijn = landelijke trend (random slope)",
                           if (nrow(rs)) sprintf(": %.3f °C/jaar (%.3f tot %.3f); SD tussen centra %.3f; p (verschil centra) = %s",
                                                 rs$verandering_per_jaar, rs$BI95_onder, rs$BI95_boven,
                                                 rs$SD_slope_tussen_centra, rs$p_extra)),
         x = "Verandering in minimale temperatuur (°C per jaar)", y = "Centrum", size = "n") + thema
  bewaar(p, map, "temperatuur_trend_per_centrum", 10, 0.35 * nrow(d) + 2.5)
}

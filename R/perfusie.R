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
#   5. Alles apart voor acuut en electief, ONGEcorrigeerd + vier covariaten
#      (gepoold "Totaal" alleen als ANALYSE_TOTAAL = TRUE)
#   6. Gevoeligheidsanalyse: alleen centrum-jaren met >= TEMP_MIN_COMPLEET compleet
#   7. Landelijke trend (marginaal, zonder centrumcorrectie, CR2 geclusterd op
#      centrum) en een exacte decompositie daarvan:
#        landelijke trend = bijdrage binnen centra + bijdrage verschuiving tussen centra
#      (y = centrumgemiddelde + afwijking daarvan; beide delen apart op jaar geregresseerd)
#   8. Tussen centra: ICC (aandeel variatie tussen centra) en 95%-voorspellingsinterval
#      voor de trend van een centrum (random slope ± 1,96 x SD)
# Daarnaast: ontbrekende temperatuur per jaar, per centrum en per centrum x jaar,
# temperatuurcategorieën en de trend per centrum (+ figuur).
# =============================================================================

TEMP_MIN_COMPLEET <- 0.80

# CR2-variantie één keer berekenen en hergebruiken (dit is de traagste stap)
cr2_jaar <- function(m, cluster) {
  V <- clubSandwich::vcovCR(m, cluster = cluster, type = "CR2")
  ci <- suppressWarnings(clubSandwich::conf_int(m, vcov = V, coefs = "jaar"))
  ct <- suppressWarnings(clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite", coefs = "jaar"))
  c(ci$beta, ci$CI_L, ci$CI_U, ct[[grep("^p_", names(ct), value = TRUE)[1]]])
}

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
  basis <- voeg_modelcovariaten_toe(basis, df)
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

  # ── Modellen: steeds eerst ongecorrigeerd, daarna gecorrigeerd ─────────────
  rij <- function(pop, g, model, d, est = NA, lo = NA, hi = NA, p = NA, sd_slope = NA,
                  p_extra = NA, extra_label = "", icc = NA)
    data.frame(populatie = pop, groep = g,
               correctie = corr_status, model = model, n = nrow(d),
               n_centra = length(unique(d$.centrum)), verandering_per_jaar = round(est, 3),
               BI95_onder = round(lo, 3), BI95_boven = round(hi, 3), p = fmt_p(p),
               SD_slope_tussen_centra = round(sd_slope, 3),
               PI95_trend_centrum = if (is.na(sd_slope)) "" else
                 sprintf("%.3f tot %.3f", est - 1.96 * sd_slope, est + 1.96 * sd_slope),
               ICC_centrum = round(icc, 3), extra_toets = extra_label, p_extra = fmt_p(p_extra))
  lme_ctrl <- nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)
  t_ci <- function(t) { q <- qt(0.975, t["DF"]); c(t["Value"], t["Value"] - q * t["Std.Error"],
                                                   t["Value"] + q * t["Std.Error"], t["p-value"]) }

  model_ri <- function(d, pop, g) tryCatch({
    m <- nlme::lme(cov_formule("y", "jaar_c", gecorrigeerd), random = ~ 1 | .centrum, data = d, method = "REML")
    r <- t_ci(summary(m)$tTable["jaar_c", ])
    # Lineariteit: lineair vs. natural spline (ML)
    m_lin <- nlme::lme(cov_formule("y", "jaar_c", gecorrigeerd), random = ~ 1 | .centrum, data = d, method = "ML")
    m_spl <- nlme::lme(cov_formule("y", "splines::ns(jaar_c, df = 3)", gecorrigeerd), random = ~ 1 | .centrum, data = d, method = "ML")
    p_lin <- anova(m_lin, m_spl)[["p-value"]][2]
    vc <- as.numeric(nlme::VarCorr(m)[, "Variance"])
    rij(pop, g, "random intercept (centrum)", d, r[1], r[2], r[3], r[4],
        p_extra = p_lin, extra_label = "afwijking van lineariteit (spline 3 df)", icc = vc[1] / sum(vc))
  }, error = function(e) rij(pop, g, paste("random intercept - mislukt:", conditionMessage(e)), d))

  model_rs <- function(d, pop, g, label = "random slope (trend per centrum)") tryCatch({
    f <- cov_formule("y", "jaar_c", gecorrigeerd)
    m <- nlme::lme(f, random = ~ 1 + jaar_c | .centrum, data = d, method = "REML", control = lme_ctrl)
    m0 <- nlme::lme(f, random = ~ 1 | .centrum, data = d, method = "REML")
    r <- t_ci(summary(m)$tTable["jaar_c", ])
    sd_s <- as.numeric(nlme::VarCorr(m)["jaar_c", "StdDev"])
    L <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(m0)))
    rij(pop, g, label, d, r[1], r[2], r[3], r[4], sd_slope = sd_s,
        p_extra = mengverdeling_p(max(L, 0)), extra_label = "trend verschilt tussen centra (LRT)")
  }, error = function(e) rij(pop, g, paste(label, "- mislukt:", conditionMessage(e)), d))

  model_cr2 <- function(d, pop, g, label = "centrum fixed effects, CR2 (Satterthwaite)") {
    if (!heeft_cs) return(rij(pop, g, paste(label, "- clubSandwich ontbreekt"), d))
    tryCatch({
      m <- lm(cov_formule("y", "jaar", gecorrigeerd, ".centrum"), data = d)
      r <- cr2_jaar(m, d$.centrum)
      rij(pop, g, label, d, r[1], r[2], r[3], r[4])
    }, error = function(e) rij(pop, g, paste(label, "- mislukt:", conditionMessage(e)), d))
  }

  # Landelijke trend (marginaal) en exacte decompositie in binnen/tussen centra.
  # Alle drie met CR2 (Satterthwaite) geclusterd op centrum.
  cr2_slope <- function(f, d) cr2_jaar(lm(f, data = d), d$.centrum)
  model_landelijk <- function(d, pop, g) {
    labels <- c("LANDELIJKE trend (marginaal, alle patienten), CR2",
                "  = bijdrage BINNEN centra (praktijkverandering), CR2",
                "  + bijdrage verschuiving TUSSEN centra (volume naar centra met ander niveau), CR2")
    if (!heeft_cs) return(do.call(rbind, lapply(labels, function(l) rij(pop, g, paste(l, "- clubSandwich ontbreekt"), d))))
    tryCatch({
      d$.gem_centrum <- ave(d$y, d$.centrum)
      d$.binnen <- d$y - d$.gem_centrum
      r1 <- cr2_slope(cov_formule("y", "jaar", gecorrigeerd), d)
      r2 <- cr2_slope(cov_formule(".binnen", "jaar", gecorrigeerd), d)
      r3 <- cr2_slope(cov_formule(".gem_centrum", "jaar", gecorrigeerd), d)
      rbind(rij(pop, g, labels[1], d, r1[1], r1[2], r1[3], r1[4]),
            rij(pop, g, labels[2], d, r2[1], r2[2], r2[3], r2[4]),
            rij(pop, g, labels[3], d, r3[1], r3[2], r3[3], r3[4]))
    }, error = function(e) rij(pop, g, paste("landelijke trend - mislukt:", conditionMessage(e)), d))
  }

  modellen <- list(); centra <- list()
  for (pop in names(populaties)) {
    dp <- populaties[[pop]]
    # centrum-jaren met voldoende complete temperatuur (gevoeligheidsanalyse)
    compleet <- tapply(!is.na(dp$y), paste(dp$.centrum, dp$jaar), mean)
    goed <- names(compleet)[compleet >= TEMP_MIN_COMPLEET]
    subsets <- split(dp, dp$groep)
    if (ANALYSE_TOTAAL || length(subsets) < 2) subsets[["Totaal"]] <- dp
    for (g in names(subsets)) for (gecorrigeerd in c(FALSE, TRUE)) {
      corr_status <- if (gecorrigeerd) "gecorrigeerd" else "ongecorrigeerd"
      d <- droplevels(subsets[[g]][!is.na(subsets[[g]]$y), , drop = FALSE])
      if (gecorrigeerd) d <- modelgevallen(d, "y")
      if (nrow(d) < 20 || nlevels(d$.centrum) < 3) {
        modellen[[length(modellen) + 1]] <- rij(pop, g,
          "niet geschat (te weinig complete observaties / centra)", d)
        next
      }
      d_s <- droplevels(d[paste(d$.centrum, d$jaar) %in% goed, , drop = FALSE])
      cat(sprintf("   temperatuurmodellen: %s, %s, %s (n = %d) ...\n",
                  pop, g, corr_status, nrow(d)))
      flush.console()
      zelfde <- nrow(d_s) == nrow(d)
      gevoelig <- if (zelfde) {
        rij(pop, g, sprintf("gevoeligheid: centrum-jaren >=%d%% compleet - alle centrum-jaren voldoen, gelijk aan hoofdanalyse",
                            100 * TEMP_MIN_COMPLEET), d)
      } else rbind(
        model_rs(d_s, pop, g, sprintf("gevoeligheid: centrum-jaren >=%d%% compleet - random slope", 100 * TEMP_MIN_COMPLEET)),
        model_cr2(d_s, pop, g, sprintf("gevoeligheid: centrum-jaren >=%d%% compleet - CR2", 100 * TEMP_MIN_COMPLEET)))
      modellen[[length(modellen) + 1]] <- rbind(
        model_landelijk(d, pop, g),
        model_ri(d, pop, g),
        model_rs(d, pop, g, "GEMIDDELDE trend binnen centra: random slope (trend per centrum)"),
        model_cr2(d, pop, g, "GEMIDDELDE trend binnen centra: centrum fixed effects, CR2 (Satterthwaite)"),
        gevoelig
      )
      # Trend per centrum: eigen lineaire regressie (met 95%-BI) en BLUP uit random-slope-model
      if (!gecorrigeerd) {
        blup <- tryCatch({
          m <- nlme::lme(y ~ jaar_c, random = ~ 1 + jaar_c | .centrum, data = d, method = "REML", control = lme_ctrl)
          setNames(coef(m)[, "jaar_c"], rownames(coef(m)))
        }, error = function(e) NULL)
        centra[[paste(pop, g)]] <- do.call(rbind, lapply(levels(d$.centrum), function(c) {
          x <- d[d$.centrum == c, ]
          ols <- if (length(unique(x$jaar)) >= 3) {
            s <- coef(summary(lm(y ~ jaar, data = x)))["jaar", ]
            q <- qt(0.975, nrow(x) - 2); c(s[1], s[1] - q * s[2], s[1] + q * s[2], s[4])
          } else rep(NA, 4)
          data.frame(populatie = pop, groep = g, centrum = c, n = nrow(x),
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

  per_centrum_jaar <- do.call(rbind, lapply(names(populaties), function(pop) {
    d <- populaties[[pop]][!is.na(populaties[[pop]]$y), ]
    a <- aggregate(y ~ groep + .centrum + jaar, data = d, FUN = function(x) c(mediaan = median(x), n = length(x)))
    data.frame(populatie = pop, groep = a$groep, centrum = a$.centrum, jaar = a$jaar,
               mediaan = a$y[, "mediaan"], n = a$y[, "n"])
  }))

  list(modellen = do.call(rbind, modellen), centra = do.call(rbind, centra),
       per_centrum_jaar = per_centrum_jaar,
       ontbrekend_jaar = ontbr_jaar, ontbrekend_centrum = ontbr_centrum,
       compleet_centrum_jaar = ontbr_matrix, categorieen = categorieen, jaar_midden = jaar_midden)
}

fig_temperatuur_centra <- function(pa, map, populatie = "met circulatiestilstand") {
  if (is.null(pa$centra)) return(NULL)
  paden <- c()
  for (g in unique(pa$centra$groep)) {
    d <- pa$centra[pa$centra$populatie == populatie & pa$centra$groep == g & !is.na(pa$centra$trend_eigen), ]
    if (!nrow(d)) next
    rs <- pa$modellen[pa$modellen$populatie == populatie & pa$modellen$groep == g &
                        pa$modellen$model == "GEMIDDELDE trend binnen centra: random slope (trend per centrum)" &
                        pa$modellen$correctie == "ongecorrigeerd", ]
    land <- pa$modellen[pa$modellen$populatie == populatie & pa$modellen$groep == g &
                          grepl("^LANDELIJKE", pa$modellen$model) &
                          pa$modellen$correctie == "ongecorrigeerd", ]
    d$centrum <- factor(d$centrum, levels = d$centrum[order(d$trend_eigen)])
    p <- ggplot(d, aes(trend_eigen, centrum)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      { if (nrow(land)) geom_vline(xintercept = land$verandering_per_jaar, colour = "#C0392B") } +
      geom_errorbarh(aes(xmin = BI95_onder, xmax = BI95_boven), height = 0.25, colour = KLEUR) +
      geom_point(aes(size = n), colour = KLEUR) +
      geom_point(aes(x = trend_random_slope_BLUP), shape = 4, size = 2.5, colour = "#C0392B") +
      labs(title = paste("Trend in minimale temperatuur per centrum -", g),
           subtitle = paste0(populatie, "; landelijke en centrumlijnen ongecorrigeerd; punt + 95%-BI = eigen regressie per centrum, x = schatting random-slope-model",
                             if (nrow(land)) sprintf("\nrode lijn = landelijke trend %.3f °C/jaar (95%%-BI %.3f tot %.3f, p = %s)",
                                                     land$verandering_per_jaar, land$BI95_onder, land$BI95_boven, land$p),
                             if (nrow(rs)) sprintf("\nSD trend tussen centra %.3f; 95%%-voorspellingsinterval %s; p (verschil tussen centra) = %s",
                                                   rs$SD_slope_tussen_centra, rs$PI95_trend_centrum, rs$p_extra)),
           x = "Verandering in minimale temperatuur (°C per jaar)", y = "Centrum", size = "n") + thema
    paden <- c(paden, bewaar(p, map, paste0("temperatuur_trend_per_centrum_", gsub("[^a-z]+", "_", tolower(g))),
                             10, 0.35 * nrow(d) + 2.8))
  }
  paden
}

# Per centrum: jaarmediaan + eigen trendlijn, naast de landelijke trend (rood)
fig_temperatuur_per_centrum_jaar <- function(pa, map, populatie = "met circulatiestilstand") {
  paden <- c()
  for (g in unique(pa$centra$groep)) {
    d <- pa$per_centrum_jaar[pa$per_centrum_jaar$populatie == populatie &
                               (pa$per_centrum_jaar$groep == g | g == "Totaal"), ]
    if (g == "Totaal") {   # gepoold: jaarmediaan over beide groepen opnieuw bepalen is niet mogelijk; gewogen gemiddelde
      d <- aggregate(cbind(mediaan = mediaan * n, n) ~ centrum + jaar, data = d, FUN = sum)
      d$mediaan <- d$mediaan / d$n
    }
    if (!nrow(d)) next
    land <- pa$modellen[pa$modellen$populatie == populatie & pa$modellen$groep == g &
                          grepl("^LANDELIJKE", pa$modellen$model) &
                          pa$modellen$correctie == "ongecorrigeerd", ]
    c_tab <- pa$centra[pa$centra$populatie == populatie & pa$centra$groep == g, ]
    d <- d[d$centrum %in% c_tab$centrum, ]
    d$centrum <- factor(d$centrum, levels = c_tab$centrum[order(-c_tab$trend_eigen)])
    midden <- weighted.mean(d$jaar, d$n)
    niveau <- weighted.mean(d$mediaan, d$n)
    lab <- setNames(sprintf("%s: %+.2f °C/jaar", c_tab$centrum, c_tab$trend_eigen), c_tab$centrum)
    p <- ggplot(d, aes(jaar, mediaan)) +
      { if (nrow(land)) geom_abline(intercept = niveau - land$verandering_per_jaar * midden,
                                    slope = land$verandering_per_jaar, colour = "#C0392B", linewidth = 0.7) } +
      geom_point(aes(size = n), colour = KLEUR, alpha = 0.7) +
      geom_smooth(aes(weight = n), method = "lm", formula = y ~ x, se = FALSE, colour = KLEUR, linewidth = 0.8) +
      facet_wrap(~centrum, labeller = labeller(centrum = lab)) +
      scale_size_area(max_size = 3, guide = "none") +
      scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 2)) +
      labs(title = paste("Minimale temperatuur per centrum over de jaren -", g),
           subtitle = paste0(populatie, "; landelijke en centrumlijnen ongecorrigeerd; punten = jaarmediaan (grootte = n), blauwe lijn = eigen trend van het centrum",
                             if (nrow(land)) sprintf("\nrode lijn = landelijke trend %.3f °C/jaar (95%%-BI %.3f tot %.3f, p = %s)",
                                                     land$verandering_per_jaar, land$BI95_onder, land$BI95_boven, land$p)),
           x = NULL, y = "Mediane minimale temperatuur (°C)") + thema
    paden <- c(paden, bewaar(p, map, paste0("temperatuur_per_centrum_per_jaar_", gsub("[^a-z]+", "_", tolower(g))), 12, 9))
  }
  paden
}

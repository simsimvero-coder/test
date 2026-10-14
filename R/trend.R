# =============================================================================
# Trendanalyse over de jaren
#
# Binaire uitkomsten: % per jaar met exact 95%-BI; trend = logistische
#   regressie met jaar als continue variabele (OR per jaar).
# Aantallen:          Poisson-regressie (incidence rate ratio per jaar).
# Temperatuur:        mediaan [IQR] per jaar; trend = lineaire regressie
#   (verandering in °C per jaar).
# Modellen per groep: ongecorrigeerd + leeftijd/geslacht/EuroSCORE II/redo.
# Beschrijvende jaarlijkse aantallen blijven ongecorrigeerde frequenties.
# =============================================================================

# TRUE bij "ja" of bij een NHR-code > 0 (bv. gastroint_compl "10 bloeding")
is_ja <- function(x) {
  x <- as.character(x)
  code <- suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", x)))
  ifelse(is.na(x), NA, x == "ja" | (!is.na(code) & code > 0))
}

# Combineer ja/nee-variabelen: ja als minstens één ja, nee als alles nee.
combineer_ja <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(rep(NA, nrow(df)))
  m <- sapply(vars, function(v) is_ja(df[[v]]))
  m <- matrix(m, nrow = nrow(df))
  ifelse(rowSums(m, na.rm = TRUE) > 0, TRUE, ifelse(rowSums(!is.na(m)) == ncol(m), FALSE, NA))
}

groepsvector <- function(d, groep)
  if (!is.null(groep) && groep %in% names(d)) droplevels(factor(d[[groep]])) else NULL

# Likelihood-ratiotoets: verschilt de trend over de jaren tussen de groepen?
p_interactie <- function(y, jaar, g, family) {
  ok <- !is.na(y) & !is.na(jaar) & !is.na(g)
  if (is.null(g) || nlevels(droplevels(g[ok])) < 2) return("")
  x <- data.frame(y = y[ok], jaar = jaar[ok], g = droplevels(g[ok]))
  m0 <- suppressWarnings(glm(y ~ jaar + g, family = family, data = x))
  m1 <- suppressWarnings(glm(y ~ jaar * g, family = family, data = x))
  a <- anova(m0, m1, test = if (family$family == "gaussian") "F" else "LRT")
  fmt_p(a[2, ncol(a)])
}

# Zet groepen om naar een lijst van subsets inclusief "Totaal".
per_groep <- function(d, groep) {
  lijst <- if (!is.null(groep) && groep %in% names(d))
    split(d, droplevels(factor(d[[groep]]))) else list()
  c(lijst, list(Totaal = d))
}

pct_per_jaar <- function(d, y, maat, groep) {
  d$.y <- y
  out <- do.call(rbind, lapply(names(gs <- per_groep(d, groep)), function(g) {
    x <- gs[[g]]; x <- x[!is.na(x$.y) & !is.na(x$.jaar), ]
    if (!nrow(x)) return(NULL)
    do.call(rbind, lapply(sort(unique(x$.jaar)), function(j) {
      xj <- x$.y[x$.jaar == j]; k <- sum(xj); n <- length(xj)
      bi <- binom.test(k, n)$conf.int
      data.frame(maat = maat, groep = g, jaar = j, n = n, events = k,
                 pct = round(100 * k / n, 1), BI_onder = round(100 * bi[1], 1),
                 BI_boven = round(100 * bi[2], 1))
    }))
  }))
  # p-waarde verschil tussen de groepen binnen elk jaar (chi-kwadraat / Fisher)
  g <- groepsvector(d, groep)
  out$p_tussen_groepen <- ""
  if (!is.null(g)) for (j in unique(out$jaar)) {
    sel <- d$.jaar == j & !is.na(d$.y) & !is.na(g)
    if (!any(sel)) next
    out$p_tussen_groepen[out$jaar == j & out$groep == "Totaal"] <- fmt_p(toets_cat(factor(d$.y[sel]), g[sel])$p)
  }
  out
}

# rbind voor tabellen met (deels) verschillende kolommen; ontbrekende kolommen worden NA/""
rbind_vul <- function(...) {
  delen <- Filter(function(x) !is.null(x) && nrow(x), list(...))
  if (!length(delen)) return(NULL)
  kol <- unique(unlist(lapply(delen, names)))
  do.call(rbind, lapply(delen, function(x) {
    for (k in setdiff(kol, names(x))) x[[k]] <- if (any(sapply(delen, function(y) k %in% names(y) && is.numeric(y[[k]])))) NA_real_ else ""
    x[kol]
  }))
}

trend_binair <- function(d, y, maat, groep, correctie = NULL) {
  d$.y <- y
  do.call(rbind, lapply(names(gs <- per_groep(d, groep)), function(g) {
    x <- gs[[g]]; x <- x[!is.na(x$.y) & !is.na(x$.jaar), ]
    rij <- data.frame(maat = maat, groep = g, n = nrow(x), events = sum(x$.y),
                      effectmaat = "OR per jaar", schatting = NA, BI_onder = NA, BI_boven = NA,
                      p = "", schatting_gecorrigeerd = NA, p_gecorrigeerd = "",
                      n_gecorrigeerd = NA_integer_, BI_onder_gecorrigeerd = NA_real_,
                      BI_boven_gecorrigeerd = NA_real_, status_gecorrigeerd = "",
                      p_trendverschil_tussen_groepen = "")
    if (g == "Totaal") rij$p_trendverschil_tussen_groepen <-
      p_interactie(d$.y, d$.jaar, groepsvector(d, groep), binomial())
    if (sum(x$.y) < 2 || sum(!x$.y) < 2 || length(unique(x$.jaar)) < 2) {
      if (length(correctie)) {
        rij$n_gecorrigeerd <- nrow(modelgevallen(x, ".y"))
        rij$status_gecorrigeerd <- "onvoldoende events of kalenderjaren"
      }
      return(rij)
    }
    m <- suppressWarnings(glm(.y ~ .jaar, family = binomial, data = x))
    s <- coef(summary(m))[".jaar", ]
    ci <- exp(s[1] + c(-1.96, 1.96) * s[2])
    rij$schatting <- round(exp(s[1]), 3); rij$BI_onder <- round(ci[1], 3)
    rij$BI_boven <- round(ci[2], 3); rij$p <- fmt_p(s[4])
    if (length(correctie)) {
      xc <- modelgevallen(x, ".y")
      rij$n_gecorrigeerd <- nrow(xc)
      if (nrow(xc) <= length(MODEL_TERMS) + 5 || sum(xc$.y) < 5 || sum(!xc$.y) < 5) {
        rij$status_gecorrigeerd <- "te weinig complete observaties/events"
      } else tryCatch({
        mc <- suppressWarnings(glm(cov_formule(".y", ".jaar", TRUE),
                                    family = binomial, data = xc))
        sc <- coef(summary(mc))[".jaar", ]
        rij$schatting_gecorrigeerd <- round(exp(sc[1]), 3)
        rij$BI_onder_gecorrigeerd <- round(exp(sc[1] - 1.96 * sc[2]), 3)
        rij$BI_boven_gecorrigeerd <- round(exp(sc[1] + 1.96 * sc[2]), 3)
        rij$p_gecorrigeerd <- fmt_p(sc[4])
        rij$status_gecorrigeerd <- "OK"
      }, error = function(e) rij$status_gecorrigeerd <<- conditionMessage(e))
    }
    rij
  }))
}

trend_aantal <- function(d, groep) {
  jaren <- seq(min(d$.jaar, na.rm = TRUE), max(d$.jaar, na.rm = TRUE))
  gs <- per_groep(d, groep)
  per_jaar <- do.call(rbind, lapply(names(gs), function(g)
    data.frame(maat = "Aantal boog-OK's", groep = g, jaar = jaren,
               n = as.integer(table(factor(gs[[g]]$.jaar, levels = jaren))))))
  toets <- do.call(rbind, lapply(names(gs), function(g) {
    x <- per_jaar[per_jaar$groep == g, ]
    m <- glm(n ~ jaar, family = poisson, data = x)
    s <- coef(summary(m))["jaar", ]
    ci <- exp(s[1] + c(-1.96, 1.96) * s[2])
    data.frame(maat = "Aantal boog-OK's", groep = g, n = sum(x$n), events = NA,
               effectmaat = "IRR per jaar", schatting = round(exp(s[1]), 3),
               BI_onder = round(ci[1], 3), BI_boven = round(ci[2], 3), p = fmt_p(s[4]),
               schatting_gecorrigeerd = NA, p_gecorrigeerd = "", p_trendverschil_tussen_groepen = "")
  }))
  x <- per_jaar[per_jaar$groep != "Totaal", ]
  if (length(unique(x$groep)) > 1)
    toets$p_trendverschil_tussen_groepen[toets$groep == "Totaal"] <-
      p_interactie(x$n, x$jaar, factor(x$groep), poisson())
  list(per_jaar = per_jaar, toets = toets)
}

trend_numeriek <- function(d, y, maat, groep) {
  d$.y <- y
  gs <- per_groep(d, groep)
  per_jaar <- do.call(rbind, lapply(names(gs), function(g) {
    x <- gs[[g]]; x <- x[!is.na(x$.y) & !is.na(x$.jaar), ]
    if (!nrow(x)) return(NULL)
    do.call(rbind, lapply(sort(unique(x$.jaar)), function(j) {
      q <- quantile(x$.y[x$.jaar == j], c(.25, .5, .75))
      data.frame(maat = maat, groep = g, jaar = j, n = sum(x$.jaar == j),
                 mediaan = q[2], Q1 = q[1], Q3 = q[3], gemiddelde = round(mean(x$.y[x$.jaar == j]), 1))
    }))
  }))
  toets <- do.call(rbind, lapply(names(gs), function(g) {
    x <- gs[[g]]; x <- x[!is.na(x$.y) & !is.na(x$.jaar), ]
    rij <- data.frame(maat = maat, groep = g, n = nrow(x), events = NA,
                      effectmaat = "verandering per jaar", schatting = NA, BI_onder = NA,
                      BI_boven = NA, p = "", schatting_gecorrigeerd = NA, p_gecorrigeerd = "",
                      p_trendverschil_tussen_groepen = "",
                      n_gecorrigeerd = NA_integer_, BI_onder_gecorrigeerd = NA_real_,
                      BI_boven_gecorrigeerd = NA_real_, status_gecorrigeerd = "")
    if (g == "Totaal") rij$p_trendverschil_tussen_groepen <-
      p_interactie(d$.y, d$.jaar, groepsvector(d, groep), gaussian())
    if (nrow(x) < 3 || length(unique(x$.jaar)) < 2) return(rij)
    s <- coef(summary(lm(.y ~ .jaar, data = x)))[".jaar", ]
    rij$schatting <- round(s[1], 2)
    rij$BI_onder <- round(s[1] - 1.96 * s[2], 2); rij$BI_boven <- round(s[1] + 1.96 * s[2], 2)
    rij$p <- fmt_p(s[4])
    xc <- modelgevallen(x, ".y")
    rij$n_gecorrigeerd <- nrow(xc)
    if (nrow(xc) > length(MODEL_TERMS) + 5 && length(unique(xc$.jaar)) >= 2) {
      tryCatch({
        t <- coef(summary(lm(cov_formule(".y", ".jaar", TRUE), data = xc)))[".jaar", ]
        rij$schatting_gecorrigeerd <- round(t[1], 3)
        rij$BI_onder_gecorrigeerd <- round(t[1] - 1.96 * t[2], 3)
        rij$BI_boven_gecorrigeerd <- round(t[1] + 1.96 * t[2], 3)
        rij$p_gecorrigeerd <- fmt_p(t[4])
        rij$status_gecorrigeerd <- "OK"
      }, error = function(e) rij$status_gecorrigeerd <<- conditionMessage(e))
    } else rij$status_gecorrigeerd <- "te weinig complete observaties"
    rij
  }))
  # p-waarde verschil tussen de groepen binnen elk jaar (Mann-Whitney / Kruskal-Wallis)
  g <- groepsvector(d, groep)
  per_jaar$p_tussen_groepen <- ""
  if (!is.null(g)) for (j in unique(per_jaar$jaar)) {
    sel <- d$.jaar == j & !is.na(d$.y) & !is.na(g)
    if (!any(sel)) next
    per_jaar$p_tussen_groepen[per_jaar$jaar == j & per_jaar$groep == "Totaal"] <-
      fmt_p(toets_num(d$.y[sel], g[sel], normaal = FALSE)$p)
  }
  rownames(per_jaar) <- NULL
  list(per_jaar = per_jaar, toets = toets)
}

# Hoofdfunctie: voert de zes trendanalyses uit.
trendanalyse <- function(df, groep = GROEP) {
  d <- voeg_modelcovariaten_toe(df, df)
  if (TREND_ALLEEN_BOOG_OK) {
    if (!BOOG_KOLOM %in% names(d)) stop("Kolom '", BOOG_KOLOM, "' niet gevonden")
    boog <- !is.na(d[[BOOG_KOLOM]]) & !as.character(d[[BOOG_KOLOM]]) %in% c("0", "nee")
    d <- d[boog, ]
  }
  d$.jaar <- suppressWarnings(as.numeric(as.character(d[[TREND_JAAR]])))
  d <- d[!is.na(d$.jaar), ]
  tijd <- if (OVERLEVING_TIJD %in% names(d)) d[[OVERLEVING_TIJD]] else rep(Inf, nrow(d))

  # (1) aantallen
  aantal <- trend_aantal(d, groep)

  # (2) mortaliteit, (3) complicaties, (4) FTR, (6) circulatiestilstand
  binair <- list()
  if (!is.null(groep) && groep %in% names(d) && "acuut" %in% d[[groep]]) {
    # Verandert het aandeel acute operaties over de jaren? (alleen totaal)
    tb <- trend_binair(d, as.character(d[[groep]]) == "acuut", "Aandeel acuut", NULL)
    aantal$toets <- rbind_vul(aantal$toets, tb)
  }
  for (label in names(TREND_MORTALITEIT)) {
    m <- TREND_MORTALITEIT[[label]]
    if (!m$var %in% names(d)) next
    y <- is_ja(d[[m$var]])
    # levend en te kort gevolgd -> uitkomst onbekend
    if (m$min_followup > 0) y[!is.na(y) & !y & (is.na(tijd) | tijd < m$min_followup)] <- NA
    binair[[label]] <- list(y = y, onderdeel = "2. Mortaliteit", correctie = TREND_CORRECTIE)
  }
  complicatie_lijst <- lapply(TREND_COMPLICATIES, combineer_ja, df = d)
  for (label in names(complicatie_lijst))
    binair[[label]] <- list(y = complicatie_lijst[[label]], onderdeel = "3. Complicaties",
                            correctie = TREND_CORRECTIE)

  # Failure to rescue: noemer = patiënten met >= 1 major complicatie
  enige_compl <- if (exists("MAJOR_COMPLICATIES")) combineer_ja(d, MAJOR_COMPLICATIES) %in% TRUE else
    Reduce(`|`, lapply(complicatie_lijst, function(x) x %in% TRUE))
  if (FTR_UITKOMST %in% names(d)) {
    ftr <- is_ja(d[[FTR_UITKOMST]])
    ftr[!enige_compl] <- NA            # alleen patiënten met een complicatie
    binair[["Failure to rescue"]] <- list(y = ftr, onderdeel = "4. Failure to rescue", correctie = TREND_CORRECTIE)
  }
  if (TREND_CIRC_ARREST %in% names(d))
    binair[["Gebruik circulatiestilstand"]] <- list(y = is_ja(d[[TREND_CIRC_ARREST]]),
                                                    onderdeel = "6. Circulatiestilstand", correctie = TREND_CORRECTIE)

  pct <- do.call(rbind, lapply(names(binair), function(l) {
    r <- pct_per_jaar(d, binair[[l]]$y, l, groep); r$onderdeel <- binair[[l]]$onderdeel; r }))
  toetsen_bin <- do.call(rbind, lapply(names(binair), function(l)
    trend_binair(d, binair[[l]]$y, l, groep, binair[[l]]$correctie)))

  # (5) temperatuur
  temp <- if (TREND_TEMPERATUUR %in% names(d))
    trend_numeriek(d, d[[TREND_TEMPERATUUR]], "Minimale temperatuur (°C)", groep)

  # Kolommen van de gepaarde modellen ook bij niet-corrigeerbare volumetellingen tonen.
  # NA-velden bij count-modellen (geen zinvolle case-mixcorrectie van nationaal volume).
  sjabloon <- toetsen_bin
  for (k in setdiff(names(sjabloon), names(aantal$toets))) {
    if (is.numeric(sjabloon[[k]])) aantal$toets[[k]] <- rep(NA_real_, nrow(aantal$toets))
    else aantal$toets[[k]] <- rep("", nrow(aantal$toets))
  }
  for (k in setdiff(names(sjabloon), names(temp$toets))) {
    if (is.numeric(sjabloon[[k]])) temp$toets[[k]] <- rep(NA_real_, nrow(temp$toets))
    else temp$toets[[k]] <- rep("", nrow(temp$toets))
  }
  toetsen <- rbind(aantal$toets[names(sjabloon)], toetsen_bin[names(sjabloon)], temp$toets[names(sjabloon)])
  rownames(toetsen) <- NULL
  per_jaar_temp <- temp$per_jaar

  if (!ANALYSE_TOTAAL && !is.null(groep) && groep %in% names(d)) {
    # Vergelijkingen tussen de groepen staan op de rij "Totaal": overzetten naar
    # de groepsrijen en daarna de gepoolde (Totaal) analyses weglaten.
    verplaats <- function(tab, sleutel, kolom) {
      if (is.null(tab) || !kolom %in% names(tab)) return(tab)
      tot <- tab[tab$groep == "Totaal", ]
      waarde <- setNames(tot[[kolom]], do.call(paste, tot[sleutel]))
      k <- do.call(paste, tab[sleutel])
      tab[[kolom]] <- ifelse(k %in% names(waarde), waarde[k], tab[[kolom]])
      tab
    }
    pct <- verplaats(pct, c("maat", "jaar"), "p_tussen_groepen")
    pct <- pct[pct$groep != "Totaal", ]
    per_jaar_temp <- verplaats(per_jaar_temp, c("maat", "jaar"), "p_tussen_groepen")
    if (!is.null(per_jaar_temp)) per_jaar_temp <- per_jaar_temp[per_jaar_temp$groep != "Totaal", ]
    toetsen <- verplaats(toetsen, "maat", "p_trendverschil_tussen_groepen")
    toetsen$groep[toetsen$maat == "Aandeel acuut"] <- "acuut vs. electief"
    toetsen <- toetsen[toetsen$groep != "Totaal", ]
    rownames(pct) <- rownames(toetsen) <- NULL
  }
  list(data = d, aantal = aantal$per_jaar, pct = pct, temp = per_jaar_temp, toetsen = toetsen)
}

# ── Figuren ─────────────────────────────────────────────────────────────────
trend_figuren <- function(tr, map) {
  jaar_as <- scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 1))
  kleuren <- scale_colour_manual(values = c("acuut" = "#C0392B", "electief aneurysma" = "#2E4F82",
                                            "Totaal" = "grey35"), aesthetics = c("colour", "fill"))
  ondertitel <- function(maten) {
    t <- tr$toetsen[tr$toetsen$maat %in% maten & tr$toetsen$p != "", ]
    if (!nrow(t)) return(NULL)
    tekst <- sprintf("%s: %s %.2f, p = %s", t$groep, t$effectmaat, t$schatting, t$p)
    inter <- t$p_trendverschil_tussen_groepen[t$p_trendverschil_tussen_groepen != ""]
    if (length(inter)) tekst <- c(tekst, paste("verschil in trend tussen groepen: p =", inter[1]))
    paste(tekst, collapse = "\n")
  }
  # Facetlabel: naam + trend-p per groep
  facetlabel <- function(maten) {
    sapply(maten, function(m) {
      t <- tr$toetsen[tr$toetsen$maat == m, ]
      regel <- paste(sprintf("%s p=%s", sub("electief aneurysma", "electief", t$groep),
                             ifelse(t$p == "", "-", t$p)), collapse = " | ")
      inter <- t$p_trendverschil_tussen_groepen[t$p_trendverschil_tussen_groepen != ""]
      if (length(inter)) regel <- paste0(regel, "\nverschil in trend p=", inter[1])
      paste0(m, "\ntrend: ", regel)
    })
  }
  paden <- c()

  # 1. Aantallen
  a <- tr$aantal[tr$aantal$groep != "Totaal", ]
  p <- ggplot(a, aes(jaar, n, fill = groep)) + geom_col() + jaar_as + kleuren +
    labs(title = "1. Aantal boog-OK's per jaar",
         subtitle = paste(ondertitel("Aantal boog-OK's"), ondertitel("Aandeel acuut"), sep = "\n"),
         x = NULL, y = "Aantal operaties", fill = NULL) + thema
  paden <- c(paden, bewaar(p, map, "trend_1_aantallen", 9, 6.5))

  lijnplot <- function(onderdeel, titel, naam, hoog = 5.5, ncol = 3) {
    d <- tr$pct[tr$pct$onderdeel == onderdeel, ]
    if (!nrow(d)) return(NULL)
    labels <- facetlabel(unique(d$maat))
    d$maat <- factor(labels[d$maat], levels = labels)
    p <- ggplot(d, aes(jaar, pct, colour = groep, fill = groep)) +
      geom_ribbon(aes(ymin = BI_onder, ymax = BI_boven), alpha = 0.10, colour = NA) +
      geom_line(linewidth = 0.8) + geom_point(aes(size = n)) +
      facet_wrap(~maat, ncol = ncol, scales = "free_y") + jaar_as + kleuren +
      scale_size_area(max_size = 3, guide = "none") +
      labs(title = titel, subtitle = "lijn = % per jaar, band = exact 95%-BI, puntgrootte = aantal patiënten",
           x = NULL, y = "%", colour = NULL, fill = NULL) +
      thema + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    bewaar(p, map, naam, 11, hoog)
  }
  paden <- c(paden,
    lijnplot("2. Mortaliteit", "2. Mortaliteit per jaar", "trend_2_mortaliteit"),
    lijnplot("3. Complicaties", "3. Complicaties per jaar", "trend_3_complicaties", 7.5, 2),
    lijnplot("4. Failure to rescue", "4. Failure to rescue per jaar (mortaliteit bij complicatie)",
             "trend_4_failure_to_rescue", 5.5, 1),
    lijnplot("6. Circulatiestilstand", "6. Gebruik circulatiestilstand per jaar",
             "trend_6_circulatiestilstand", 5.5, 1))

  # 5. Temperatuur
  if (!is.null(tr$temp)) {
    d <- tr$data[!is.na(tr$data[[TREND_TEMPERATUUR]]), ]
    d$groep <- if (!is.null(GROEP) && GROEP %in% names(d)) as.character(d[[GROEP]]) else "Totaal"
    p <- ggplot(d, aes(factor(.jaar), .data[[TREND_TEMPERATUUR]], fill = groep)) +
      geom_boxplot(outlier.size = 0.6, alpha = 0.85) + kleuren +
      labs(title = "5. Minimale temperatuur per jaar", subtitle = ondertitel("Minimale temperatuur (°C)"),
           x = NULL, y = "Minimale temperatuur (°C)", fill = NULL) + thema
    paden <- c(paden, bewaar(p, map, "trend_5_temperatuur", 10, 6.5))
  }
  paden[!sapply(paden, is.null)]
}

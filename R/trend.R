# =============================================================================
# Trendanalyse over de jaren
#
# Binaire uitkomsten: % per jaar met exact 95%-BI; trend = logistische
#   regressie met jaar als continue variabele (OR per jaar).
# Aantallen:          Poisson-regressie (incidence rate ratio per jaar).
# Temperatuur:        mediaan [IQR] per jaar; trend = lineaire regressie
#   (verandering in °C per jaar).
# Alles per groep (acuut / electief) en voor het totaal.
# =============================================================================

is_ja <- function(x) ifelse(is.na(x), NA, as.character(x) == "ja")

# Combineer ja/nee-variabelen: ja als minstens één ja, nee als alles nee.
combineer_ja <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(rep(NA, nrow(df)))
  m <- sapply(vars, function(v) is_ja(df[[v]]))
  m <- matrix(m, nrow = nrow(df))
  ifelse(rowSums(m, na.rm = TRUE) > 0, TRUE, ifelse(rowSums(!is.na(m)) == ncol(m), FALSE, NA))
}

# Zet groepen om naar een lijst van subsets inclusief "Totaal".
per_groep <- function(d, groep) {
  lijst <- if (!is.null(groep) && groep %in% names(d))
    split(d, droplevels(factor(d[[groep]]))) else list()
  c(lijst, list(Totaal = d))
}

pct_per_jaar <- function(d, y, maat, groep) {
  d$.y <- y
  do.call(rbind, lapply(names(gs <- per_groep(d, groep)), function(g) {
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
}

trend_binair <- function(d, y, maat, groep, correctie = NULL) {
  d$.y <- y
  do.call(rbind, lapply(names(gs <- per_groep(d, groep)), function(g) {
    x <- gs[[g]]; x <- x[!is.na(x$.y) & !is.na(x$.jaar), ]
    rij <- data.frame(maat = maat, groep = g, n = nrow(x), events = sum(x$.y),
                      effectmaat = "OR per jaar", schatting = NA, BI_onder = NA, BI_boven = NA,
                      p = "", schatting_gecorrigeerd = NA, p_gecorrigeerd = "")
    if (sum(x$.y) < 2 || sum(!x$.y) < 2 || length(unique(x$.jaar)) < 2) return(rij)
    m <- suppressWarnings(glm(.y ~ .jaar, family = binomial, data = x))
    s <- coef(summary(m))[".jaar", ]
    ci <- exp(s[1] + c(-1.96, 1.96) * s[2])
    rij$schatting <- round(exp(s[1]), 3); rij$BI_onder <- round(ci[1], 3)
    rij$BI_boven <- round(ci[2], 3); rij$p <- fmt_p(s[4])
    correctie <- intersect(correctie, names(x))
    if (length(correctie)) {
      f <- reformulate(c(".jaar", paste0("`", correctie, "`")), ".y")
      mc <- suppressWarnings(glm(f, family = binomial, data = x))
      sc <- coef(summary(mc))[".jaar", ]
      rij$schatting_gecorrigeerd <- round(exp(sc[1]), 3)
      rij$p_gecorrigeerd <- fmt_p(sc[4])
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
               schatting_gecorrigeerd = NA, p_gecorrigeerd = "")
  }))
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
                      BI_boven = NA, p = "", schatting_gecorrigeerd = NA, p_gecorrigeerd = "")
    if (nrow(x) < 3 || length(unique(x$.jaar)) < 2) return(rij)
    s <- coef(summary(lm(.y ~ .jaar, data = x)))[".jaar", ]
    rij$schatting <- round(s[1], 2)
    rij$BI_onder <- round(s[1] - 1.96 * s[2], 2); rij$BI_boven <- round(s[1] + 1.96 * s[2], 2)
    rij$p <- fmt_p(s[4])
    rij
  }))
  rownames(per_jaar) <- NULL
  list(per_jaar = per_jaar, toets = toets)
}

# Hoofdfunctie: voert de zes trendanalyses uit.
trendanalyse <- function(df, groep = GROEP) {
  d <- df
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
  for (label in names(TREND_MORTALITEIT)) {
    m <- TREND_MORTALITEIT[[label]]
    if (!m$var %in% names(d)) next
    y <- is_ja(d[[m$var]])
    # levend en te kort gevolgd -> uitkomst onbekend
    y[!is.na(y) & !y & !is.na(tijd) & tijd < m$min_followup] <- NA
    binair[[label]] <- list(y = y, onderdeel = "2. Mortaliteit", correctie = TREND_CORRECTIE)
  }
  complicatie_lijst <- lapply(TREND_COMPLICATIES, combineer_ja, df = d)
  for (label in names(complicatie_lijst))
    binair[[label]] <- list(y = complicatie_lijst[[label]], onderdeel = "3. Complicaties")

  enige_compl <- Reduce(`|`, lapply(complicatie_lijst, function(x) x %in% TRUE))
  if (FTR_UITKOMST %in% names(d)) {
    ftr <- is_ja(d[[FTR_UITKOMST]])
    ftr[!enige_compl] <- NA            # alleen patiënten met een complicatie
    binair[["Failure to rescue"]] <- list(y = ftr, onderdeel = "4. Failure to rescue")
  }
  if (TREND_CIRC_ARREST %in% names(d))
    binair[["Gebruik circulatiestilstand"]] <- list(y = is_ja(d[[TREND_CIRC_ARREST]]),
                                                    onderdeel = "6. Circulatiestilstand")

  pct <- do.call(rbind, lapply(names(binair), function(l) {
    r <- pct_per_jaar(d, binair[[l]]$y, l, groep); r$onderdeel <- binair[[l]]$onderdeel; r }))
  toetsen_bin <- do.call(rbind, lapply(names(binair), function(l)
    trend_binair(d, binair[[l]]$y, l, groep, binair[[l]]$correctie)))

  # (5) temperatuur
  temp <- if (TREND_TEMPERATUUR %in% names(d))
    trend_numeriek(d, d[[TREND_TEMPERATUUR]], "Minimale temperatuur (°C)", groep)

  toetsen <- rbind(aantal$toets, toetsen_bin, temp$toets)
  rownames(toetsen) <- NULL
  list(data = d, aantal = aantal$per_jaar, pct = pct, temp = temp$per_jaar, toetsen = toetsen)
}

# ── Figuren ─────────────────────────────────────────────────────────────────
trend_figuren <- function(tr, map) {
  jaar_as <- scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 1))
  kleuren <- scale_colour_manual(values = c("acuut" = "#C0392B", "electief aneurysma" = "#2E4F82",
                                            "Totaal" = "grey35"), aesthetics = c("colour", "fill"))
  ondertitel <- function(maten) {
    t <- tr$toetsen[tr$toetsen$maat %in% maten & tr$toetsen$p != "", ]
    if (!nrow(t)) return(NULL)
    paste(sprintf("%s: %s %.2f, p = %s", t$groep, t$effectmaat, t$schatting, t$p), collapse = "\n")
  }
  paden <- c()

  # 1. Aantallen
  a <- tr$aantal[tr$aantal$groep != "Totaal", ]
  p <- ggplot(a, aes(jaar, n, fill = groep)) + geom_col() + jaar_as + kleuren +
    labs(title = "1. Aantal boog-OK's per jaar", subtitle = ondertitel("Aantal boog-OK's"),
         x = NULL, y = "Aantal operaties", fill = NULL) + thema
  paden <- c(paden, bewaar(p, map, "trend_1_aantallen", 9, 5.5))

  lijnplot <- function(onderdeel, titel, naam, hoog = 5.5, ncol = 3) {
    d <- tr$pct[tr$pct$onderdeel == onderdeel, ]
    if (!nrow(d)) return(NULL)
    d$maat <- factor(d$maat, levels = unique(d$maat))
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
    paden <- c(paden, bewaar(p, map, "trend_5_temperatuur", 10, 5.5))
  }
  paden[!sapply(paden, is.null)]
}

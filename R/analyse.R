# =============================================================================
# Beschrijvende tabellen, groepsvergelijkingen, O/E-ratio en overleving
# =============================================================================

fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

is_normaal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(unique(x)) < 3) return(FALSE)
  if (length(x) > 5000) { set.seed(1); x <- sample(x, 5000) }
  shapiro.test(x)$p.value >= ALPHA_NORMALITEIT
}

beschrijf_num <- function(x, normaal) {
  x <- x[!is.na(x)]
  if (!length(x)) return("-")
  if (normaal) return(sprintf("%.1f (%.1f)", mean(x), sd(x)))
  q <- quantile(x, c(.25, .5, .75))
  sprintf("%.1f [%.1f-%.1f]", q[2], q[1], q[3])
}

toets_num <- function(x, g, normaal) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- droplevels(factor(g[ok]))
  if (nlevels(g) < 2 || any(table(g) < 3)) return(list(p = NA, toets = ""))
  if (nlevels(g) == 2) {
    if (normaal) list(p = t.test(x ~ g)$p.value, toets = "Welch t-toets")
    else list(p = suppressWarnings(wilcox.test(x ~ g)$p.value), toets = "Mann-Whitney U")
  } else {
    if (normaal) list(p = oneway.test(x ~ g)$p.value, toets = "Welch ANOVA")
    else list(p = kruskal.test(x ~ g)$p.value, toets = "Kruskal-Wallis")
  }
}

toets_cat <- function(x, g) {
  tab <- table(x, g)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(list(p = NA, toets = ""))
  verwacht <- suppressWarnings(chisq.test(tab)$expected)
  if (any(verwacht < 5)) {
    p <- if (all(dim(tab) == 2)) fisher.test(tab)$p.value
         else fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
    return(list(p = p, toets = "Fisher exact"))
  }
  list(p = chisq.test(tab)$p.value, toets = "chi-kwadraat")
}

# Maakt een tabel in "Tabel 1"-stijl voor de opgegeven variabelen.
maak_tabel <- function(df, variabelen, groep = NULL) {
  variabelen <- setdiff(intersect(variabelen, names(df)), groep)
  heeft_groep <- !is.null(groep) && groep %in% names(df)
  if (heeft_groep) {
    g <- droplevels(factor(df[[groep]]))
    niveaus <- levels(g)
    kolomnamen <- sprintf("%s (n=%d)", niveaus, as.integer(table(g)))
  }

  rijen <- list()
  voeg_toe <- function(variabele, totaal, per_groep = NULL, p = NA, toets = "", ontbr = "") {
    r <- data.frame(Variabele = variabele, check.names = FALSE)
    r[[sprintf("Totaal (n=%d)", nrow(df))]] <- totaal
    if (heeft_groep) {
      for (i in seq_along(niveaus)) r[[kolomnamen[i]]] <- if (is.null(per_groep)) "" else per_groep[i]
      r[["p-waarde"]] <- fmt_p(p)
      r[["Toets"]] <- toets
    }
    r[["Ontbrekend"]] <- as.character(ontbr)
    rijen[[length(rijen) + 1]] <<- r
  }

  for (v in variabelen) {
    x <- df[[v]]
    ontbr <- sum(is.na(x))
    if (is.numeric(x)) {
      normaal <- is_normaal(x)
      label <- paste0(v, if (normaal) ", gem. (SD)" else ", mediaan [IQR]")
      per_groep <- if (heeft_groep) sapply(niveaus, function(n) beschrijf_num(x[g %in% n], normaal))
      t <- if (heeft_groep) toets_num(x, g, normaal) else list(p = NA, toets = "")
      voeg_toe(label, beschrijf_num(x, normaal), per_groep, t$p, t$toets, ontbr)
    } else {
      x <- droplevels(factor(x))
      t <- if (heeft_groep) toets_cat(x, g) else list(p = NA, toets = "")
      n_pct <- function(sel, niv) {
        geldig <- sum(!is.na(x[sel]))
        if (!geldig) return("-")
        n <- sum(x[sel] == niv, na.rm = TRUE)
        sprintf("%d (%.1f%%)", n, 100 * n / geldig)
      }
      janee <- identical(levels(x), c("nee", "ja")) || identical(levels(x), "ja") || identical(levels(x), "nee")
      if (TOON_ALLEEN_JA && janee) {
        per_groep <- if (heeft_groep) sapply(niveaus, function(n) n_pct(g %in% n, "ja"))
        voeg_toe(paste0(v, ", n (%)"), n_pct(TRUE, "ja"), per_groep, t$p, t$toets, ontbr)
      } else {
        voeg_toe(paste0(v, ", n (%)"), "", NULL, t$p, t$toets, ontbr)
        for (niv in levels(x)) {
          per_groep <- if (heeft_groep) sapply(niveaus, function(n) n_pct(g %in% n, niv))
          voeg_toe(paste0("    ", niv), n_pct(TRUE, niv), per_groep)
        }
      }
    }
  }
  do.call(rbind, rijen)
}

# Geobserveerde vs. verwachte (EuroSCORE II) mortaliteit, met exact 95%-BI.
oe_ratio <- function(df, per, uitkomst = "overleden_zh", verwacht = "EuroscoreII") {
  if (!all(c(per, uitkomst, verwacht) %in% names(df))) return(NULL)
  d <- df[!is.na(df[[uitkomst]]) & !is.na(df[[verwacht]]) & !is.na(df[[per]]), ]
  d$obs <- as.integer(d[[uitkomst]] == "ja")
  d$exp <- d[[verwacht]] / 100
  groepen <- split(d, droplevels(factor(d[[per]])))
  if (ANALYSE_TOTAAL) groepen[["Totaal"]] <- d
  out <- do.call(rbind, lapply(names(groepen), function(n) {
    x <- groepen[[n]]
    O <- sum(x$obs); E <- sum(x$exp)
    bi <- poisson.test(O)$conf.int / E
    data.frame(groep = n, n = nrow(x), overleden = O,
               geobserveerd_pct = round(100 * O / nrow(x), 1),
               verwacht_pct = round(100 * E / nrow(x), 1),
               OE_ratio = round(O / E, 2),
               BI95_onder = round(bi[1], 2), BI95_boven = round(bi[2], 2),
               p = fmt_p(poisson.test(O, E)$p.value))   # H0: O/E = 1
  }))
  names(out)[1] <- per
  # Verschil tussen groepen, gecorrigeerd voor het verwachte risico:
  # logistische regressie met logit(EuroSCORE II) als offset
  dd <- d[d$exp > 0 & d$exp < 1, ]
  dd$.g <- droplevels(factor(dd[[per]]))
  if (nlevels(dd$.g) > 1) {
    m0 <- glm(obs ~ offset(qlogis(exp)), family = binomial, data = dd)
    m1 <- glm(obs ~ .g + offset(qlogis(exp)), family = binomial, data = dd)
    attr(out, "p_tussen_groepen") <- anova(m0, m1, test = "LRT")[2, "Pr(>Chi)"]
    out$p_verschil_tussen_groepen <- fmt_p(attr(out, "p_tussen_groepen"))
  }
  out
}

# Kaplan-Meier blijft ONGECORRIGEERD, per klinische groep. Geen extrapolatie.
# Nog lopende follow-up wordt op administratieve analysehorizon gecensureerd.
overleving <- function(df, groep = NULL) {
  nodig <- c(OVERLEVING_TIJD, OVERLEVING_STATUS)
  if (!all(nodig %in% names(df))) return(NULL)
  d <- df[!is.na(df[[OVERLEVING_TIJD]]) & !is.na(df[[OVERLEVING_STATUS]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # Bronrijen selecteren met dezelfde logische index.
  sel <- !is.na(df[[OVERLEVING_TIJD]]) & !is.na(df[[OVERLEVING_STATUS]])
  d$.tijd <- pmin(df[[OVERLEVING_TIJD]][sel], OVERLEVING_MAX_DAGEN)
  d$.status <- as.integer(as.character(d[[OVERLEVING_STATUS]]) == "ja" &
                            d[[OVERLEVING_TIJD]] <= OVERLEVING_MAX_DAGEN)
  d <- d[is.finite(d$.tijd) & d$.tijd >= 0, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  heeft_groep <- !is.null(groep) && groep %in% names(d)
  if (heeft_groep) {
    d <- d[!is.na(d[[groep]]), , drop = FALSE]
    d$.groep <- droplevels(factor(d[[groep]]))
    if (!nrow(d)) return(NULL)
    fit <- survival::survfit(survival::Surv(.tijd, .status) ~ .groep, data = d)
    p <- if (nlevels(d$.groep) >= 2) fmt_p(stats::pchisq(
      survival::survdiff(survival::Surv(.tijd, .status) ~ .groep, data = d)$chisq,
      nlevels(d$.groep) - 1, lower.tail = FALSE)) else ""
  } else {
    d$.groep <- factor("Totaal")
    fit <- survival::survfit(survival::Surv(.tijd, .status) ~ 1, data = d)
    p <- ""
  }
  tijden <- c(30, 120, 365, 3 * 365, 5 * 365)
  tabellen <- lapply(split(d, d$.groep, drop = TRUE), function(x) {
    f <- survival::survfit(survival::Surv(.tijd, .status) ~ 1, data = x)
    # Alleen tijdstippen waarvoor daadwerkelijk follow-up aanwezig is.
    geldige <- tijden[tijden <= max(x$.tijd)]
    if (!length(geldige)) return(NULL)
    v <- summary(f, times = geldige, extend = FALSE)
    if (!length(v$time)) return(NULL)
    data.frame(groep = as.character(x$.groep[1]), dag = v$time,
               n_initieel = nrow(x), at_risk = v$n.risk,
               events_cumulatief = vapply(v$time, function(t) sum(x$.status[x$.tijd <= t]), integer(1)),
               overleving_pct = round(100 * v$surv, 1),
               BI95 = sprintf("%.1f-%.1f", 100 * v$lower, 100 * v$upper),
               laatste_observatie_dag = max(x$.tijd))
  })
  tabel <- do.call(rbind, Filter(Negate(is.null), tabellen))
  if (is.null(tabel)) return(NULL)
  tabel$p_logrank <- c(p, rep("", nrow(tabel) - 1))
  rownames(tabel) <- NULL
  list(fit = fit, data = d, tabel = tabel, logrank_p = p, heeft_groep = heeft_groep)
}

# Eerste paper: tijd-tot-overlijden met kalenderjaar van operatie als determinant.
# Twee Cox-modellen *binnen* acuut en electief (HR per operatief kalenderjaar).
# Geen correctie van Kaplan-Meier zelf. Clustervariantie per ziekenhuis.
overleving_cox_trend <- function(df, groep = GROEP) {
  nodig <- c(OVERLEVING_TIJD, OVERLEVING_STATUS, TREND_JAAR, CENTRUM_KOLOM, groep)
  if (!all(nodig %in% names(df))) return(NULL)
  d <- data.frame(groep = as.character(df[[groep]]),
                  jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                  centrum = factor(df[[CENTRUM_KOLOM]]),
                  tijd_bron = df[[OVERLEVING_TIJD]],
                  status_bron = as.character(df[[OVERLEVING_STATUS]]))
  d <- voeg_modelcovariaten_toe(d, df)
  d <- d[complete.cases(d[, c("groep", "jaar", "centrum", "tijd_bron", "status_bron")]) &
           d$tijd_bron >= 0, , drop = FALSE]
  d$tijd <- pmin(d$tijd_bron, COX_MAX_DAGEN)
  d$dood <- as.integer(d$status_bron == "ja" & d$tijd_bron <= COX_MAX_DAGEN)
  if (!nrow(d)) return(NULL)
  d$jaar_c <- d$jaar - round(mean(d$jaar))
  res <- list()
  for (g in sort(unique(d$groep))) for (adjusted in c(FALSE, TRUE)) {
    x <- droplevels(d[d$groep == g, , drop = FALSE])
    if (adjusted) x <- modelgevallen(x, "dood")
    resultaat <- list(hr = NA_real_, onder = NA_real_, boven = NA_real_,
                      p = NA_real_, ph_p = NA_real_, status = "")
    if (nrow(x) < 30 || sum(x$dood) < 10 || length(unique(x$jaar)) < 2) {
      resultaat$status <- "onvoldoende observaties/events/kalenderjaren"
    } else tryCatch({
      f <- cov_formule("survival::Surv(tijd, dood)", "jaar_c", adjusted, "cluster(centrum)")
      m <- survival::coxph(f, data = x, ties = "efron")
      b <- unname(coef(m)["jaar_c"])
      se <- sqrt(diag(vcov(m)))["jaar_c"]
      if (!is.finite(b) || !is.finite(se)) stop("niet-schatbare kalenderjaarcoefficient")
      resultaat$hr <- round(exp(b), 3)
      resultaat$onder <- round(exp(b - 1.96 * se), 3)
      resultaat$boven <- round(exp(b + 1.96 * se), 3)
      resultaat$p <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
      resultaat$ph_p <- tryCatch(as.numeric(survival::cox.zph(m)$table["jaar_c", "p"]),
                                  error = function(e) NA_real_)
      resultaat$status <- "OK"
    }, error = function(e) resultaat$status <<- paste("mislukt:", conditionMessage(e)))
    res[[length(res) + 1]] <- data.frame(
      groep = g, correctie = if (adjusted) "gecorrigeerd" else "ongecorrigeerd",
      n = nrow(x), overlijden = sum(x$dood), HR_per_jaar = resultaat$hr,
      BI95_onder = resultaat$onder, BI95_boven = resultaat$boven,
      p = fmt_p(resultaat$p), p_PH_kalenderjaar = fmt_p(resultaat$ph_p),
      status = resultaat$status)
  }
  do.call(rbind, res)
}

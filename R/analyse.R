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
  groepen[["Totaal"]] <- d
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
    out$p_verschil_tussen_groepen <- c(rep("", nrow(out) - 1), fmt_p(attr(out, "p_tussen_groepen")))
  }
  out
}

# Kaplan-Meier-overleving; geeft model, overlevingstabel en log-rank p terug.
overleving <- function(df, groep = NULL) {
  if (!all(c(OVERLEVING_TIJD, OVERLEVING_STATUS) %in% names(df))) return(NULL)
  d <- df[!is.na(df[[OVERLEVING_TIJD]]) & !is.na(df[[OVERLEVING_STATUS]]), ]
  d$.tijd <- d[[OVERLEVING_TIJD]]
  d$.status <- as.integer(d[[OVERLEVING_STATUS]] == "ja")
  heeft_groep <- !is.null(groep) && groep %in% names(d)
  if (heeft_groep) {
    d <- d[!is.na(d[[groep]]), ]
    d$.groep <- droplevels(factor(d[[groep]]))
    fit <- survival::survfit(survival::Surv(.tijd, .status) ~ .groep, data = d)
    p <- fmt_p(1 - pchisq(survival::survdiff(survival::Surv(.tijd, .status) ~ .groep, data = d)$chisq,
                          nlevels(d$.groep) - 1))
  } else {
    fit <- survival::survfit(survival::Surv(.tijd, .status) ~ 1, data = d)
    p <- ""
  }
  tijden <- c(30, 120, 365, 3 * 365, 5 * 365)
  s <- summary(fit, times = tijden, extend = TRUE)
  tabel <- data.frame(
    groep = if (heeft_groep) sub("^\\.groep=", "", as.character(s$strata)) else "Totaal",
    dag = s$time, at_risk = s$n.risk,
    overleving_pct = round(100 * s$surv, 1),
    BI95 = sprintf("%.1f-%.1f", 100 * s$lower, 100 * s$upper),
    se = s$std.err
  )
  if (heeft_groep) {
    # Verschil in overleving op elk tijdpunt (z-toets op Kaplan-Meier-schattingen;
    # bij > 2 groepen een Wald chi-kwadraattoets op gelijkheid van alle groepen)
    tabel$p_tijdpunt <- ""
    for (t in unique(tabel$dag)) {
      r <- tabel[tabel$dag == t & tabel$se > 0, ]
      if (nrow(r) < 2) next
      w <- 1 / r$se^2
      gem <- sum(w * r$overleving_pct / 100) / sum(w)
      chi <- sum(w * (r$overleving_pct / 100 - gem)^2)
      tabel$p_tijdpunt[tabel$dag == t][1] <- fmt_p(1 - pchisq(chi, nrow(r) - 1))
    }
    tabel$p_logrank <- c(p, rep("", nrow(tabel) - 1))
  }
  tabel$se <- NULL
  list(fit = fit, data = d, tabel = tabel, logrank_p = p, heeft_groep = heeft_groep)
}

# =============================================================================
# Klinische uitkomsten, variatie tussen centra, opnameduur en (relatieve) overleving
#
# Alles apart voor acuut en electief: ruwe en gecorrigeerde trends (4 covariaten).
#
# 1. Uitkomsttrends (OR per kalenderjaar) voor mortaliteit, major complicaties
#    (composiet), secundaire morbiditeit en failure to rescue:
#    - primair: logistisch mixed model met random intercept per centrum,
#      gecorrigeerd voor UITKOMST_COVARIATEN (lme4::glmer)
#    - gevoeligheid: centrum-fixed-effects + CR2 (Satterthwaite)
#    - Benjamini-Hochberg over SECUNDAIRE_UITKOMSTEN, per groep
#    30-dagen- en 1-jaarsmortaliteit alleen bij patiënten met volledige
#    follow-up (overleden, of in leven en >= 30 / 365 dagen gevolgd).
# 2. Variatie tussen centra in orgaanbescherming: per centrum mediane temperatuur,
#    % diepe hypothermie, % circulatiestilstand en duur; ICC's (lineair, of op de
#    latente schaal voor ja/nee-uitkomsten: var / (var + pi^2/3)).
# 3. Temperatuur zonder vs. met circulatiestilstand, binnen centra (CR2).
# 4. Opnameduur (alle patiënten en levend ontslagen) en trend (ratio per jaar).
# 5. Kaplan-Meier en relatieve overleving (Ederer II) op 1, 3 en 5 jaar.
# =============================================================================

code_positief <- function(x) {
  x <- as.character(x)
  code <- suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", x)))
  ifelse(is.na(x), NA, x == "ja" | (!is.na(code) & code > 0))
}

samengesteld <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(rep(NA, nrow(df)))
  m <- matrix(sapply(vars, function(v) code_positief(df[[v]])), nrow = nrow(df))
  ifelse(rowSums(m, na.rm = TRUE) > 0, TRUE, ifelse(rowSums(!is.na(m)) == ncol(m), FALSE, NA))
}

# Event-free discharge: levend ontslagen zonder major complicatie
# (MAJOR_COMPLICATIES). Onbekend als overlijden of een complicatie onbekend is
# en er geen bekende gebeurtenis is.
maak_event_free <- function(df) {
  if (!"overleden_zh" %in% names(df)) return(df)
  gebeurtenis <- samengesteld(df, c("overleden_zh", MAJOR_COMPLICATIES))
  df$event_free_ontslag <- factor(ifelse(is.na(gebeurtenis), NA, ifelse(gebeurtenis, "nee", "ja")),
                                  levels = c("nee", "ja"))
  df
}

# Logistisch model met random intercept per centrum. Primair lme4::glmer (Laplace);
# lukt dat niet (bv. een lme4/Matrix-versieconflict), dan MASS::glmmPQL als terugval.
logit_mixed <- function(f, d, term) {
  f_re <- update(f, . ~ . + (1 | .centrum))
  fout <- ""
  if (requireNamespace("lme4", quietly = TRUE)) {
    r <- tryCatch({
      m <- suppressMessages(suppressWarnings(
        lme4::glmer(f_re, data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa"))))
      out <- list(model = "glmer (Laplace)", var_centrum = as.numeric(lme4::VarCorr(m)$.centrum))
      if (!is.null(term)) { s <- coef(summary(m))[term, ]; out$est <- s[1]; out$se <- s[2]; out$p <- s[4] }
      out
    }, error = function(e) { fout <<- conditionMessage(e); NULL })
    if (!is.null(r)) return(r)
  }
  r <- tryCatch({
    m <- suppressMessages(MASS::glmmPQL(f, random = ~ 1 | .centrum, family = binomial, data = d, verbose = FALSE))
    out <- list(model = "glmmPQL (terugval)", var_centrum = as.numeric(nlme::VarCorr(m)[1, "Variance"]))
    if (!is.null(term)) { t <- summary(m)$tTable[term, ]; out$est <- t["Value"]; out$se <- t["Std.Error"]; out$p <- t["p-value"] }
    out
  }, error = function(e) list(model = "mislukt", fout = paste("mislukt:", fout, conditionMessage(e))))
  r
}

voortgang <- function(...) { cat("   ", sprintf(...), "\n", sep = ""); flush.console() }

# ── 1. Uitkomsttrends ───────────────────────────────────────────────────────
uitkomst_trends <- function(df, groep = GROEP) {
  if (!requireNamespace("lme4", quietly = TRUE))
    warning("lme4 ontbreekt: voor sommige modellen kan PQL worden gebruikt; interpreteer met voorzichtigheid")
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE) && UITKOMST_CR2
  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, min_fu) {
    y <- code_positief(df[[var]])
    y <- mortaliteit_fu(y, tijd, df[[TREND_JAAR]], min_fu)
    weg <- attr(y, "jaren_uitgesloten")
    if (length(weg)) message("   ", var, ": operatiejaar ", paste(weg, collapse = ", "),
                             " niet meegenomen (follow-up nog geen ", min_fu, " dagen)")
    y
  }
  major <- samengesteld(df, MAJOR_COMPLICATIES)
  ftr <- code_positief(df[[FTR_UITKOMST]]); ftr[!(major %in% TRUE)] <- NA
  uitkomsten <- list(
    list("Ziekenhuismortaliteit", "mortaliteit", dood("overleden_zh", 0)),
    list("30-dagenmortaliteit", "mortaliteit", dood("overl30d", 30)),
    list("1-jaarsmortaliteit", "mortaliteit", dood("overl1j", 365)),
    list("Major complicaties", "composiet", major),
    list("Event-free discharge", "composiet",
         if ("event_free_ontslag" %in% names(df)) code_positief(df$event_free_ontslag) else rep(NA, nrow(df)))
  )
  for (l in names(SECUNDAIRE_UITKOMSTEN)) if (SECUNDAIRE_UITKOMSTEN[[l]] %in% names(df))
    uitkomsten[[length(uitkomsten) + 1]] <- list(l, "secundair", code_positief(df[[SECUNDAIRE_UITKOMSTEN[[l]]]]))
  uitkomsten[[length(uitkomsten) + 1]] <- list("Failure to rescue", "FTR", ftr)
  basis <- data.frame(jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                      .centrum = factor(df[[CENTRUM_KOLOM]]),
                      groep = as.character(df[[groep]]))
  basis <- voeg_modelcovariaten_toe(basis, df)
  basis$jaar_c <- basis$jaar - round(mean(basis$jaar, na.rm = TRUE))
  fit_1 <- function(d, adjusted) {
    d <- if (adjusted) modelgevallen(d, "y") else d
    n <- nrow(d); events <- sum(d$y)
    a <- list(n = n, events = events, OR = NA_real_, BI95 = "", p = NA_real_, model = "")
    if (n < 30 || events < 5 || n - events < 5 || length(unique(d$jaar)) < 2) {
      a$model <- "onvoldoende observaties of events"; return(a)
    }
    f <- cov_formule("y", "jaar_c", adjusted)
    r <- tryCatch(logit_mixed(f, d, "jaar_c"), error = function(e)
      list(model = paste("mislukt:", conditionMessage(e))))
    a$model <- r$model
    if (!is.null(r$est) && is.finite(r$est) && is.finite(r$se)) {
      a$OR <- round(exp(r$est), 3)
      a$BI95 <- sprintf("%.3f-%.3f", exp(r$est - 1.96 * r$se), exp(r$est + 1.96 * r$se))
      a$p <- unname(r$p)
    } else if (!is.null(r$fout)) a$model <- paste(a$model, r$fout)
    a
  }
  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (u in uitkomsten) {
    d <- basis
    d$y <- u[[3]]
    d <- droplevels(d[d$groep == g & !is.na(d$y) & !is.na(d$jaar) & !is.na(d$.centrum), , drop = FALSE])
    voortgang("uitkomstmodel: %s - %s (n ruw = %d)", g, u[[1]], nrow(d))
    onge <- fit_1(d, FALSE)
    gec <- fit_1(d, TRUE)
    cr2 <- list(OR = NA_real_, BI95 = "", p = NA_real_)
    dc <- modelgevallen(d, "y")
    if (heeft_cs && nrow(dc) >= 30 && sum(dc$y) >= 5 && sum(!dc$y) >= 5 && nlevels(dc$.centrum) >= 3) {
      tryCatch({
        m <- glm(cov_formule("y", "jaar", TRUE, ".centrum"), family = binomial, data = dc)
        r <- cr2_jaar(m, dc$.centrum)
        cr2 <- list(OR = round(exp(r[1]), 3),
                    BI95 = sprintf("%.3f-%.3f", exp(r[2]), exp(r[3])), p = unname(r[4]))
      }, error = function(e) message("CR2 overgeslagen voor ", g, ": ", conditionMessage(e)))
    }
    res[[length(res) + 1]] <- data.frame(
      groep = g, uitkomst = u[[1]], type = u[[2]],
      n_ongecorrigeerd = onge$n, events_ongecorrigeerd = onge$events,
      pct_ongecorrigeerd = if (onge$n) round(100 * onge$events / onge$n, 1) else NA_real_,
      model_ongecorrigeerd = onge$model,
      OR_ongecorrigeerd = onge$OR, BI95_ongecorrigeerd = onge$BI95, p_ongecorrigeerd = onge$p,
      n_gecorrigeerd = gec$n, events_gecorrigeerd = gec$events,
      model_gecorrigeerd = gec$model,
      OR_gecorrigeerd = gec$OR, BI95_gecorrigeerd = gec$BI95, p_gecorrigeerd = gec$p,
      OR_CR2_gecorrigeerd = cr2$OR, BI95_CR2_gecorrigeerd = cr2$BI95,
      p_CR2_gecorrigeerd = cr2$p)
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  # BH apart over secundaire eindpunten per groep én correctiestatus.
  for (naam in c("p_ongecorrigeerd", "p_gecorrigeerd", "p_CR2_gecorrigeerd")) {
    qnaam <- sub("^p_", "q_BH_", naam)
    out[[qnaam]] <- NA_real_
    for (g in unique(out$groep)) {
      i <- which(out$groep == g & out$type == "secundair")
      out[[qnaam]][i] <- p.adjust(out[[naam]][i], method = "BH")
    }
  }
  for (naam in grep("^(p_|q_BH_)", names(out), value = TRUE))
    out[[naam]] <- vapply(out[[naam]], fmt_p, character(1))
  attr(out, "noot") <- paste0("Per groep eerst ongecorrigeerd, dan gecorrigeerd voor ",
                              paste(UITKOMST_COVARIATEN, collapse = ", "),
                              "; CR2 alleen gecorrigeerd; complete-case-aantallen apart")
  out
}

# ── 2. Variatie tussen centra in orgaanbescherming ──────────────────────────
centrum_variatie <- function(df, groep = GROEP) {
  heeft_lme4 <- requireNamespace("lme4", quietly = TRUE)
  arrest_kol <- if (CIRC_ARREST_KOLOM %in% names(df)) CIRC_ARREST_KOLOM else "circ_arrest"
  d0 <- data.frame(groep = as.character(df[[groep]]), .centrum = factor(df[[CENTRUM_KOLOM]]),
                   temp = df[[TREND_TEMPERATUUR]],
                   diep = ifelse(is.na(df[[TREND_TEMPERATUUR]]), NA, df[[TREND_TEMPERATUUR]] <= TEMP_GRENZEN[1]),
                   arrest = code_positief(df[[arrest_kol]]),
                   duur = df$circ_arrest_duur_bij_arrest)

  d0 <- voeg_modelcovariaten_toe(d0, df)
  per_centrum <- do.call(rbind, lapply(split(d0, list(d0$groep, d0$.centrum), drop = TRUE), function(x)
    data.frame(groep = x$groep[1], centrum = as.character(x$.centrum[1]), n = nrow(x),
               mediaan_temp = round(median(x$temp, na.rm = TRUE), 1),
               pct_diepe_hypothermie = round(100 * mean(x$diep, na.rm = TRUE), 1),
               pct_circ_arrest = round(100 * mean(x$arrest, na.rm = TRUE), 1),
               mediaan_duur_arrest = round(median(x$duur, na.rm = TRUE), 1))))
  rownames(per_centrum) <- NULL
  per_centrum <- per_centrum[order(per_centrum$groep, suppressWarnings(as.numeric(per_centrum$centrum))), ]

  bereik <- do.call(rbind, lapply(split(per_centrum, per_centrum$groep), function(x) {
    b <- function(v) sprintf("%.1f - %.1f", min(v, na.rm = TRUE), max(v, na.rm = TRUE))
    data.frame(groep = x$groep[1], mediaan_temp = b(x$mediaan_temp),
               pct_diepe_hypothermie = b(x$pct_diepe_hypothermie),
               pct_circ_arrest = b(x$pct_circ_arrest), mediaan_duur_arrest = b(x$mediaan_duur_arrest))
  }))

  icc_lin <- function(x, y) tryCatch({
    x <- x[!is.na(x[[y]]), ]
    m <- nlme::lme(as.formula(paste(y, "~ 1")), random = ~ 1 | .centrum, data = x, method = "REML")
    v <- as.numeric(nlme::VarCorr(m)[, "Variance"]); v[1] / sum(v)
  }, error = function(e) NA)
  icc_bin <- function(x, y) {
    x <- x[!is.na(x[[y]]), ]; x$y <- as.integer(x[[y]])
    r <- logit_mixed(y ~ 1, x, NULL)
    if (is.null(r$var_centrum)) NA else r$var_centrum / (r$var_centrum + pi^2 / 3)
  }
  icc <- do.call(rbind, lapply(split(d0, d0$groep), function(x) {
    on.exit(NULL)
    voortgang("ICC's: %s", x$groep[1])
    data.frame(groep = x$groep[1], row.names = NULL,
               ICC_min_temp = round(icc_lin(x, "temp"), 3),
               ICC_diepe_hypothermie_latent = round(icc_bin(x, "diep"), 3),
               ICC_circ_arrest_latent = round(icc_bin(x, "arrest"), 3),
               ICC_duur_arrest = round(icc_lin(x, "duur"), 3))
  }))

  # Temperatuur zonder vs. met circulatiestilstand (verschil binnen centra, CR2)
  zonder_met <- do.call(rbind, lapply(split(d0, d0$groep), function(x) {
    x <- droplevels(x[!is.na(x$temp) & !is.na(x$arrest), ])
    q <- function(v) sprintf("%.1f [%.1f-%.1f]", median(v), quantile(v, .25), quantile(v, .75))
    r <- if (requireNamespace("clubSandwich", quietly = TRUE) && sum(!x$arrest) >= 5) tryCatch({
      x$zonder <- as.integer(!x$arrest)
      m <- lm(temp ~ zonder + .centrum, data = x)
      V <- clubSandwich::vcovCR(m, cluster = x$.centrum, type = "CR2")
      ci <- clubSandwich::conf_int(m, vcov = V, coefs = "zonder")
      ct <- clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite", coefs = "zonder")
      c(ci$beta, ci$CI_L, ci$CI_U, ct[[grep("^p_", names(ct), value = TRUE)[1]]])
    }, error = function(e) rep(NA, 4)) else rep(NA, 4)
    xc <- modelgevallen(x, "temp")
    rc <- if (requireNamespace("clubSandwich", quietly = TRUE) && nrow(xc) > 20 &&
              sum(xc$arrest) >= 5 && sum(!xc$arrest) >= 5 && nlevels(xc$.centrum) >= 3) tryCatch({
      xc$zonder <- as.integer(!xc$arrest)
      m <- lm(cov_formule("temp", "zonder", TRUE, ".centrum"), data = xc)
      V <- clubSandwich::vcovCR(m, cluster = xc$.centrum, type = "CR2")
      ci <- clubSandwich::conf_int(m, vcov = V, coefs = "zonder")
      ct <- clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite", coefs = "zonder")
      c(ci$beta, ci$CI_L, ci$CI_U, ct[[grep("^p_", names(ct), value = TRUE)[1]]])
    }, error = function(e) rep(NA_real_, 4)) else rep(NA_real_, 4)
    data.frame(groep = x$groep[1], n_zonder = sum(!x$arrest), temp_zonder = q(x$temp[!x$arrest]),
               n_met = sum(x$arrest), temp_met = q(x$temp[x$arrest]),
               n_gecorrigeerd = nrow(xc),
               verschil_binnen_centra = round(r[1], 2),
               BI95 = sprintf("%.2f tot %.2f", r[2], r[3]), p = fmt_p(r[4]),
               verschil_gecorrigeerd = round(rc[1], 2),
               BI95_gecorrigeerd = sprintf("%.2f tot %.2f", rc[2], rc[3]),
               p_gecorrigeerd = fmt_p(rc[4]))
  }))

  list(per_centrum = per_centrum, bereik = bereik, icc = icc, zonder_met = zonder_met)
}

# ── 4. Opnameduur ───────────────────────────────────────────────────────────
opnameduur_analyse <- function(df, groep = GROEP) {
  if (!"opnameduur" %in% names(df)) return(NULL)
  d0 <- data.frame(groep = as.character(df[[groep]]), .centrum = factor(df[[CENTRUM_KOLOM]]),
                   jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                   los = df$opnameduur, levend = !code_positief(df$overleden_zh))
  d0 <- voeg_modelcovariaten_toe(d0, df)
  d0$jaar_c <- d0$jaar - round(mean(d0$jaar, na.rm = TRUE))
  do.call(rbind, lapply(split(d0, d0$groep), function(x) do.call(rbind,
    lapply(c("alle patiënten", "levend ontslagen"), function(pop) do.call(rbind,
      lapply(c(FALSE, TRUE), function(gecorrigeerd) {
        y <- x[!is.na(x$los) & !is.na(x$jaar) & !is.na(x$.centrum) &
                   (pop == "alle patiënten" | x$levend %in% TRUE), , drop = FALSE]
        if (gecorrigeerd) y <- modelgevallen(y, "los")
        q <- if (nrow(y)) sprintf("%.1f [%.1f-%.1f]", median(y$los), quantile(y$los, .25), quantile(y$los, .75)) else ""
        r <- tryCatch({
          if (nrow(y) < 20 || nlevels(droplevels(y$.centrum)) < 2) stop("te weinig observaties/centra")
          m <- nlme::lme(cov_formule("log(los + 1)", "jaar_c", gecorrigeerd),
                         random = ~ 1 | .centrum, data = droplevels(y), method = "REML")
          t <- summary(m)$tTable["jaar_c", ]
          c(exp(t["Value"]), exp(t["Value"] - 1.96 * t["Std.Error"]),
            exp(t["Value"] + 1.96 * t["Std.Error"]), t["p-value"])
        }, error = function(e) rep(NA_real_, 4))
        data.frame(groep = x$groep[1], populatie = pop,
                   correctie = if (gecorrigeerd) "gecorrigeerd" else "ongecorrigeerd",
                   n = nrow(y), mediaan_IQR = q, ratio_per_jaar = round(r[1], 3),
                   BI95 = sprintf("%.3f-%.3f", r[2], r[3]), p = fmt_p(r[4]))
      }))))))
}

# ── 5. Relatieve overleving (Ederer II) ─────────────────────────────────────
lees_levenstabel <- function() {
  if (is.null(LEVENSTABEL_BESTAND)) return(NULL)
  pad <- zoek_bestand(LEVENSTABEL_BESTAND)
  if (is.null(pad)) {
    message("LET OP - levenstabel '", LEVENSTABEL_BESTAND, "' niet gevonden; relatieve overleving overgeslagen. ",
            "Zie README voor het maken van dit bestand.")
    return(NULL)
  }
  kop <- readLines(pad, n = 1, warn = FALSE, encoding = "UTF-8")
  sep <- if (grepl(";", kop)) ";" else ","
  lt <- read.csv(pad, sep = sep, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8-BOM")
  # Kolommen herkennen, ook in een CBS StatLine-export ("Perioden", "Sterftekans (kans)", ...)
  kies <- function(patroon) { k <- grep(patroon, names(lt), ignore.case = TRUE)[1]
    if (is.na(k)) stop("Kolom '", patroon, "' niet gevonden in de levenstabel"); lt[[k]] }
  getal <- function(x) suppressWarnings(as.numeric(gsub(",", ".", sub("^\\D*(\\d+[.,]?\\d*).*$", "\\1", x))))
  lt <- data.frame(jaar = getal(kies("jaar|period")),
                   geslacht = ifelse(grepl("^(m|1)", tolower(trimws(kies("geslacht")))), "man", "vrouw"),
                   leeftijd = getal(kies("leeftijd")),
                   sterftekans = suppressWarnings(as.numeric(gsub(",", ".", kies("sterftekans")))))
  lt <- lt[complete.cases(lt), ]
  if (max(lt$sterftekans) > 1) lt$sterftekans <- lt$sterftekans / 1000   # per 1000 opgegeven
  lt$hazard <- -log(1 - pmin(lt$sterftekans, 0.999))
  lt
}

relatieve_overleving <- function(df, groep = GROEP, tijden = c(1, 3, 5)) {
  lt <- lees_levenstabel()
  if (is.null(lt) || !all(c(OVERLEVING_TIJD, OVERLEVING_STATUS, "leeftijd", "geslacht") %in% names(df))) return(NULL)
  d <- data.frame(groep = as.character(df[[groep]]), t = df[[OVERLEVING_TIJD]] / 365.25,
                  dood = as.integer(code_positief(df[[OVERLEVING_STATUS]])),
                  leeftijd = df$leeftijd, geslacht = as.character(df$geslacht),
                  jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))))
  d <- d[complete.cases(d) & d$geslacht %in% c("man", "vrouw"), ]
  sleutel <- function(j, s, a) paste(pmin(pmax(j, min(lt$jaar)), max(lt$jaar)), s,
                                     pmin(floor(a), max(lt$leeftijd)))
  h_tab <- setNames(lt$hazard, paste(lt$jaar, lt$geslacht, lt$leeftijd))
  stap <- 1 / 12
  grid <- seq(stap, max(tijden), by = stap)

  do.call(rbind, lapply(split(d, d$groep), function(x) {
    km <- survival::survfit(survival::Surv(t, dood) ~ 1, data = x)
    s_km <- summary(km, times = c(0, grid), extend = TRUE)
    s_obs <- s_km$surv[-1]
    # Ederer II: verwachte overleving per interval = gemiddelde onder wie nog at risk is
    s_exp <- numeric(length(grid)); cum <- 1
    for (k in seq_along(grid)) {
      start <- grid[k] - stap
      r <- x[x$t >= start, ]
      h <- h_tab[sleutel(r$jaar + floor(start), r$geslacht, r$leeftijd + start)]
      cum <- cum * mean(exp(-h * stap), na.rm = TRUE)
      s_exp[k] <- cum
    }
    idx <- match(round(tijden / stap), round(grid / stap))
    lo <- s_km$lower[-1][idx]; hi <- s_km$upper[-1][idx]
    data.frame(groep = x$groep[1], n = nrow(x), jaar = tijden,
               overleving_KM = round(100 * s_obs[idx], 1),
               KM_BI95 = sprintf("%.1f-%.1f", 100 * lo, 100 * hi),
               verwachte_overleving = round(100 * s_exp[idx], 1),
               relatieve_overleving = round(100 * s_obs[idx] / s_exp[idx], 1),
               RS_BI95 = sprintf("%.1f-%.1f", 100 * lo / s_exp[idx], 100 * hi / s_exp[idx]),
               ratio_tov_1_jaar_pct = round(100 * (s_obs[idx] / s_exp[idx]) / (s_obs[idx[1]] / s_exp[idx[1]]), 1))
  }))
}

# Standardised mortality ratio (SMR): waargenomen sterfte / verwachte sterfte in de
# algemene bevolking van dezelfde leeftijd, geslacht en kalenderjaar (CBS-levenstabel).
# Verwacht = som over patiënten van de cumulatieve populatiehazard tijdens hun eigen
# follow-up binnen het tijdvenster (maandstappen). 95%-BI exact (Poisson).
# Vensters: 0-1 jaar, 1-5 jaar (alleen patiënten die het eerste jaar overleefden)
# en 0-5 jaar.
smr_analyse <- function(df, groep = GROEP,
                        vensters = list("0-1 jaar" = c(0, 1), "1-5 jaar" = c(1, 5), "0-5 jaar" = c(0, 5))) {
  lt <- lees_levenstabel()
  if (is.null(lt) || !all(c(OVERLEVING_TIJD, OVERLEVING_STATUS, "leeftijd", "geslacht") %in% names(df))) return(NULL)
  d <- data.frame(groep = as.character(df[[groep]]), t = df[[OVERLEVING_TIJD]] / 365.25,
                  dood = as.integer(code_positief(df[[OVERLEVING_STATUS]])),
                  leeftijd = df$leeftijd, geslacht = as.character(df$geslacht),
                  jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))))
  d <- d[complete.cases(d) & d$geslacht %in% c("man", "vrouw"), ]
  h_tab <- setNames(lt$hazard, paste(lt$jaar, lt$geslacht, lt$leeftijd))
  hazard <- function(r, s) h_tab[paste(pmin(pmax(r$jaar + floor(s), min(lt$jaar)), max(lt$jaar)),
                                       r$geslacht, pmin(floor(r$leeftijd + s), max(lt$leeftijd)))]
  stap <- 1 / 12
  out <- do.call(rbind, lapply(split(d, d$groep), function(x) do.call(rbind, lapply(names(vensters), function(v) {
    a <- vensters[[v]][1]; b <- vensters[[v]][2]
    if (a > 0) x <- x[x$t > a, ]           # at risk bij begin van het venster
    O <- sum(x$dood == 1 & x$t <= b)
    E <- 0
    for (s in seq(a, b - stap, by = stap)) {
      r <- x[x$t > s, ]
      if (!nrow(r)) break
      E <- E + sum(hazard(r, s) * pmin(stap, r$t - s), na.rm = TRUE)
    }
    ci <- if (O > 0) c(qchisq(0.025, 2 * O) / 2, qchisq(0.975, 2 * (O + 1)) / 2) / E else c(0, qchisq(0.975, 2) / 2 / E)
    data.frame(groep = x$groep[1], venster = v, n = nrow(x), waargenomen = O, verwacht = round(E, 1),
               SMR = round(O / E, 2), BI95 = sprintf("%.2f-%.2f", ci[1], ci[2]),
               p = fmt_p(poisson.test(O, E)$p.value),
               extra_sterfte_per_1000_persoonsjaren = round(1000 * (O - E) / sum(pmin(x$t, b) - a), 1))
  }))))
  rownames(out) <- NULL
  out
}

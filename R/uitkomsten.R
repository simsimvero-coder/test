# =============================================================================
# Klinische uitkomsten, variatie tussen centra, opnameduur en (relatieve) overleving
#
# Alles apart voor acuut en electief, zonder correctie voor urgentie of EuroSCORE.
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

voortgang <- function(...) { cat("   ", sprintf(...), "\n", sep = ""); flush.console() }

# ── 1. Uitkomsttrends ───────────────────────────────────────────────────────
uitkomst_trends <- function(df, groep = GROEP) {
  heeft_lme4 <- requireNamespace("lme4", quietly = TRUE)
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE) && UITKOMST_CR2
  if (!heeft_lme4) message("LET OP - pakket 'lme4' ontbreekt; mixed models overgeslagen. Installeer met: ",
                           "install.packages(\"lme4\", repos = \"https://cloud.r-project.org\", type = \"binary\")")

  tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(Inf, nrow(df))
  dood <- function(var, min_fu) {
    y <- code_positief(df[[var]])
    y[!is.na(y) & !y & !is.na(tijd) & tijd < min_fu] <- NA   # levend maar te kort gevolgd
    y
  }
  major <- samengesteld(df, MAJOR_COMPLICATIES)
  ftr <- code_positief(df[[FTR_UITKOMST]]); ftr[!(major %in% TRUE)] <- NA

  uitkomsten <- list(
    list("Ziekenhuismortaliteit", "mortaliteit", dood("overleden_zh", 0)),
    list("30-dagenmortaliteit (volledige follow-up)", "mortaliteit", dood("overl30d", 30)),
    list("1-jaarsmortaliteit (volledige follow-up)", "mortaliteit", dood("overl1j", 365)),
    list(paste0("Major complicaties (", paste(MAJOR_COMPLICATIES, collapse = "/"), ")"), "composiet", major)
  )
  for (l in names(SECUNDAIRE_UITKOMSTEN)) if (SECUNDAIRE_UITKOMSTEN[[l]] %in% names(df))
    uitkomsten[[length(uitkomsten) + 1]] <- list(l, "secundair", code_positief(df[[SECUNDAIRE_UITKOMSTEN[[l]]]]))
  uitkomsten[[length(uitkomsten) + 1]] <- list("Failure to rescue (overleden in zh bij major complicatie)", "FTR", ftr)

  basis <- data.frame(jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                      .centrum = factor(df[[CENTRUM_KOLOM]]),
                      groep = as.character(df[[groep]]))
  cov <- intersect(UITKOMST_COVARIATEN, names(df))
  for (k in cov) basis[[k]] <- if (k == "leeftijd") df[[k]] / 10 else df[[k]]
  basis$jaar_c <- basis$jaar - round(mean(basis$jaar, na.rm = TRUE))

  rij <- function(g, u, d, ...) {
    a <- list(...)
    data.frame(groep = g, uitkomst = u[[1]], type = u[[2]], n = nrow(d), events = sum(d$y),
               pct = round(100 * mean(d$y), 1),
               OR_primair = a$or %||% NA, BI95_primair = a$bi %||% "", p_primair = a$p %||% NA,
               OR_CR2 = a$or2 %||% NA, BI95_CR2 = a$bi2 %||% "", p_CR2 = a$p2 %||% NA)
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  f_cov <- if (length(cov)) paste("+", paste(cov, collapse = " + ")) else ""

  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (u in uitkomsten) {
    d <- basis; d$y <- u[[3]]
    d <- d[d$groep %in% g & !is.na(d$y) & !is.na(d$jaar) & !is.na(d$.centrum), ]
    d <- droplevels(d[complete.cases(d[, c("y", cov), drop = FALSE]), ])
    if (nrow(d) < 30 || sum(d$y) < 5) { res[[length(res) + 1]] <- rij(g, u, d); next }
    voortgang("uitkomstmodel: %s - %s (n = %d)", g, u[[1]], nrow(d))
    a <- list()
    if (heeft_lme4) tryCatch({
      m <- lme4::glmer(as.formula(paste("y ~ jaar_c", f_cov, "+ (1 | .centrum)")), data = d, family = binomial,
                       control = lme4::glmerControl(optimizer = "bobyqa"))
      s <- coef(summary(m))["jaar_c", ]
      a$or <- round(exp(s[1]), 3); a$p <- s[4]
      a$bi <- sprintf("%.3f-%.3f", exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]))
    }, error = function(e) a$bi <<- paste("mislukt:", conditionMessage(e)))
    if (heeft_cs) tryCatch({
      m <- glm(as.formula(paste("y ~ jaar", f_cov, "+ .centrum")), data = d, family = binomial)
      r <- cr2_jaar(m, d$.centrum)
      a$or2 <- round(exp(r[1]), 3); a$bi2 <- sprintf("%.3f-%.3f", exp(r[2]), exp(r[3])); a$p2 <- r[4]
    }, error = function(e) a$bi2 <<- paste("mislukt:", conditionMessage(e)))
    res[[length(res) + 1]] <- do.call(rij, c(list(g, u, d), a))
  }
  out <- do.call(rbind, res)

  # Benjamini-Hochberg over de secundaire uitkomsten, per groep
  out$q_BH_primair <- NA; out$q_BH_CR2 <- NA
  for (g in unique(out$groep)) {
    sel <- out$groep == g & out$type == "secundair"
    out$q_BH_primair[sel] <- p.adjust(out$p_primair[sel], method = "BH")
    out$q_BH_CR2[sel] <- p.adjust(out$p_CR2[sel], method = "BH")
  }
  for (k in c("p_primair", "p_CR2", "q_BH_primair", "q_BH_CR2"))
    out[[k]] <- sapply(out[[k]], function(p) if (is.na(p)) "" else fmt_p(p))
  attr(out, "noot") <- paste0("Gecorrigeerd voor: ", paste(cov, collapse = ", "),
                              " (leeftijd per 10 jaar); OR per kalenderjaar; BH over: ",
                              paste(names(SECUNDAIRE_UITKOMSTEN), collapse = ", "))
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
  icc_bin <- function(x, y) if (!heeft_lme4) NA else tryCatch({
    x <- x[!is.na(x[[y]]), ]; x$.y <- as.integer(x[[y]])
    m <- lme4::glmer(.y ~ 1 + (1 | .centrum), data = x, family = binomial)
    v <- as.numeric(lme4::VarCorr(m)$.centrum); v / (v + pi^2 / 3)
  }, error = function(e) NA)
  icc <- do.call(rbind, lapply(split(d0, d0$groep), function(x) {
    voortgang("ICC's: %s", x$groep[1])
    data.frame(groep = x$groep[1],
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
    data.frame(groep = x$groep[1], n_zonder = sum(!x$arrest), temp_zonder = q(x$temp[!x$arrest]),
               n_met = sum(x$arrest), temp_met = q(x$temp[x$arrest]),
               verschil_binnen_centra = round(r[1], 2),
               BI95 = sprintf("%.2f tot %.2f", r[2], r[3]), p = fmt_p(r[4]))
  }))

  list(per_centrum = per_centrum, bereik = bereik, icc = icc, zonder_met = zonder_met)
}

# ── 4. Opnameduur ───────────────────────────────────────────────────────────
opnameduur_analyse <- function(df, groep = GROEP) {
  if (!"opnameduur" %in% names(df)) return(NULL)
  d0 <- data.frame(groep = as.character(df[[groep]]), .centrum = factor(df[[CENTRUM_KOLOM]]),
                   jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                   los = df$opnameduur, levend = !code_positief(df$overleden_zh))
  d0$jaar_c <- d0$jaar - round(mean(d0$jaar, na.rm = TRUE))
  do.call(rbind, lapply(split(d0, d0$groep), function(x) do.call(rbind, lapply(c("alle patiënten", "levend ontslagen"), function(pop) {
    y <- x[!is.na(x$los) & !is.na(x$jaar) & (pop == "alle patiënten" | x$levend %in% TRUE), ]
    q <- sprintf("%.1f [%.1f-%.1f]", median(y$los), quantile(y$los, .25), quantile(y$los, .75))
    tr <- tryCatch({   # log(opnameduur + 1): ratio per jaar (geometrisch)
      m <- nlme::lme(log(los + 1) ~ jaar_c, random = ~ 1 | .centrum, data = y, method = "REML")
      t <- summary(m)$tTable["jaar_c", ]
      c(exp(t["Value"]), exp(t["Value"] - 1.96 * t["Std.Error"]), exp(t["Value"] + 1.96 * t["Std.Error"]), t["p-value"])
    }, error = function(e) rep(NA, 4))
    data.frame(groep = x$groep[1], populatie = pop, n = nrow(y), mediaan_IQR = q,
               ratio_per_jaar = round(tr[1], 3), BI95 = sprintf("%.3f-%.3f", tr[2], tr[3]), p = fmt_p(tr[4]))
  }))))
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

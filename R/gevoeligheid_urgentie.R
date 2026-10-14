# =============================================================================
# Gevoeligheidsanalyse op basis van alléén urgentie (ONAFHANKELIJK van pathologie)
# 10 electief; 20 urgent (APART); 30/40 spoed/redding.
# 10 versus 30/40 wordt als afzonderlijke tweegroepsvergelijking onderzocht.
# Geen extra groepen voor chronische dissectie of overige pathologie.
# =============================================================================

# Directe vergelijking: spoed/redding versus electief, niet versus urgent.
# Een centrumintercept in beide logistische modellen; de patiëntcovariaten
# alléén in het gecorrigeerde model, ZONDER EuroSCORE II: urgentie is een
# onderdeel van EuroSCORE II, dus correctie daarvoor zou het te vergelijken
# verschil (deels) wegcorrigeren.
vergelijk_urgentie_es <- function(d) {
  if (!all(c("urgentie_sens", CENTRUM_KOLOM) %in% names(d))) return(NULL)
  d$spoed <- as.integer(as.character(d$urgentie_sens) == "spoed/redding (30/40)")
  d$.centrum <- factor(d[[CENTRUM_KOLOM]])
  d <- voeg_modelcovariaten_toe(d, d)
  tijd <- if (OVERLEVING_TIJD %in% names(d)) d[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(d))
  uitkomsten <- c("Ziekenhuismortaliteit" = "overleden_zh",
                 "30-dagenmortaliteit" = "overl30d",
                 "1-jaarsmortaliteit" = "overl1j")
  res <- list()
  for (label in names(uitkomsten)) {
    v <- uitkomsten[[label]]
    if (!v %in% names(d)) next
    y <- code_positief(d[[v]])
    termijn <- if (v == "overl30d") 30 else if (v == "overl1j") 365 else 0
    if (termijn > 0) y[!is.na(y) & !y & (is.na(tijd) | tijd < termijn)] <- NA
    for (adj in c(FALSE, TRUE)) {
      x <- d; x$y <- as.numeric(y)
      x <- x[!is.na(x$y) & !is.na(x$.centrum), , drop = FALSE]
      if (adj) x <- modelgevallen(x, "y", setdiff(MODEL_TERMS, "EuroscoreII"))
      est <- ci1 <- ci2 <- pv <- NA_real_; status <- "onvoldoende gegevens/events"
      if (nrow(x) >= 30 && sum(x$y) >= 5 && sum(1 - x$y) >= 5 &&
          length(unique(x$spoed)) == 2) {
        rr <- tryCatch(logit_mixed(cov_formule("y", "spoed", adj, zonder = "EuroscoreII"),
                                   droplevels(x), "spoed"),
                       error = function(e) list(model = conditionMessage(e)))
        status <- rr$model
        if (!is.null(rr$est) && is.finite(rr$est) && is.finite(rr$se)) {
          est <- exp(rr$est); ci1 <- exp(rr$est - 1.96 * rr$se)
          ci2 <- exp(rr$est + 1.96 * rr$se); pv <- rr$p
        }
      }
      res[[length(res) + 1]] <- data.frame(
        uitkomst = label, vergelijking = "spoed/redding (30/40) vs electief (10)",
        model = if (adj) "gecorrigeerd" else "ongecorrigeerd",
        n = nrow(x), events = sum(x$y), OR = round(est, 3),
        BI95_onder = round(ci1, 3), BI95_boven = round(ci2, 3),
        p = fmt_p(pv), model_status = status)
    }
  }
  if (!length(res)) return(NULL)
  do.call(rbind, res)
}

# Minimaal gemeten temperatuur: tweegroepsverschil in °C, met een random
# centrumintercept. Geen causale interpretatie van de vergelijking.
vergelijk_urgentie_temp <- function(d) {
  if (!all(c("min_temp", "urgentie_sens", CENTRUM_KOLOM) %in% names(d))) return(NULL)
  d$spoed <- as.integer(as.character(d$urgentie_sens) == "spoed/redding (30/40)")
  d$.centrum <- factor(d[[CENTRUM_KOLOM]])
  d$y <- d$min_temp
  d <- voeg_modelcovariaten_toe(d, d)
  met_arrest <- if (CIRC_ARREST_KOLOM %in% names(d))
    as.character(d[[CIRC_ARREST_KOLOM]]) == "ja" else rep(NA, nrow(d))
  res <- list()
  for (pop in c("alle procedures", "met circulatiestilstand"))
    for (adj in c(FALSE, TRUE)) {
      idx <- if (pop == "alle procedures") rep(TRUE, nrow(d)) else met_arrest %in% TRUE
      x <- droplevels(d[idx & !is.na(d$y) & !is.na(d$.centrum), , drop = FALSE])
      if (adj) x <- modelgevallen(x, "y", setdiff(MODEL_TERMS, "EuroscoreII"))
      eff <- onder <- boven <- pv <- NA_real_
      status <- "onvoldoende observaties/centra"
      if (nrow(x) >= 30 && length(unique(x$spoed)) == 2 &&
          nlevels(x$.centrum) >= 2) {
        r <- tryCatch({
          m <- nlme::lme(cov_formule("y", "spoed", adj, zonder = "EuroscoreII"), random = ~ 1 | .centrum,
                         data = x, method = "REML")
          t <- summary(m)$tTable["spoed", ]
          list(est = unname(t["Value"]), se = unname(t["Std.Error"]),
               p = unname(t["p-value"]), status = "lme (random centrumintercept)")
        }, error = function(e) list(status = conditionMessage(e)))
        status <- r$status
        if (!is.null(r$est) && is.finite(r$est) && is.finite(r$se)) {
          eff <- r$est; onder <- eff - 1.96 * r$se; boven <- eff + 1.96 * r$se
          pv <- r$p
        }
      }
      res[[length(res) + 1]] <- data.frame(
        populatie = pop,
        vergelijking = "spoed/redding (30/40) vs electief (10)",
        model = if (adj) "gecorrigeerd" else "ongecorrigeerd",
        n = nrow(x), verschil_C = round(eff, 3), BI95_onder = round(onder, 3),
        BI95_boven = round(boven, 3), p = fmt_p(pv), model_status = status)
    }
  do.call(rbind, res)
}

analyse_urgentie_gevoeligheid <- function(df_alle, map = UITVOER_MAP) {
  stopifnot(URGENTIE_SENS_KOLOM %in% names(df_alle))
  g <- URGENTIE_SENS_KOLOM
  label_e <- "electief (10)"
  label_u <- "urgent (20)"
  label_s <- "spoed/redding (30/40)"

  # Cruciaal: vertrek van df_alle vóór selectie op de *primaire* groepsdefinitie.
  # Anders raakt men ongezien geldige urgentiegevallen kwijt.
  d <- droplevels(df_alle[!is.na(df_alle[[g]]), , drop = FALSE])
  d_es <- droplevels(d[as.character(d[[g]]) %in% c(label_e, label_s), , drop = FALSE])
  d_u <- droplevels(d[as.character(d[[g]]) == label_u, , drop = FALSE])
  aantallen <- data.frame(
    categorie = c(label_e, label_u, label_s, "onbekende/andere urgentie", "totaal cohort"),
    n = c(sum(as.character(df_alle[[g]]) == label_e, na.rm = TRUE),
          sum(as.character(df_alle[[g]]) == label_u, na.rm = TRUE),
          sum(as.character(df_alle[[g]]) == label_s, na.rm = TRUE),
          sum(is.na(df_alle[[g]])), nrow(df_alle)))

  # Toets electief vs spoed/redding; urgent blijft afzonderlijk zichtbaar in
  # de driedelige analyses en wordt nooit samengevoegd met electief.
  # Voorkom onduidelijke foutmeldingen in zeer kleine datasets.
  if (!all(c(label_e, label_s) %in% as.character(d_es[[g]]))) {
    warning("Urgentiegevoeligheid: onvoldoende electieve of spoed/redding-patiënten. ",
            "Alleen de indelingstabellen zijn beschikbaar.")
    return(list(aantallen = aantallen, onderdelen = list()))
  }

  # Alleen relevante klinische tabellen; géén nieuwe pathologie-categorieën.
  klinische_vars <- lapply(TABELLEN, function(v) setdiff(
    v, c("aorta_asc_path", "aortaboog_path", "aorta_desc_path", "descendens_pathologie")))
  tab_es <- lapply(klinische_vars, function(v) maak_tabel(d_es, v, g))
  tab_urgent <- lapply(klinische_vars, function(v) {
    if (!nrow(d_u)) return(NULL)
    maak_tabel(d_u, v, groep = NULL)
  })

  # Modellen: voor ELKE groep eerst ongecorrigeerd, dan gecorrigeerd voor
  # leeftijd, geslacht, EuroSCORE II en eerdere cardiochirurgie.
  # De 2-groeps trends laten daarnaast p voor verschil in jaartend E vs S zien.
  # Rekenintensieve modules apart afvangen zodat overige uitkomsten blijven staan.
  veilig <- function(naam, fn) tryCatch(fn(), error = function(e) {
    warning("Urgentiegevoeligheid, onderdeel '", naam, "' mislukt: ", conditionMessage(e))
    NULL
  })
  tr3 <- veilig("trends drie urgentiegroepen", function() trendanalyse(d, g))
  tr_es <- veilig("trendverschil electief vs spoed/redding", function() trendanalyse(d_es, g))
  ut <- veilig("klinische uitkomsttrends", function() uitkomst_trends(d, g))
  per_mod <- veilig("periode-uitkomsten", function() periode_uitkomstmodellen(d, g))
  cox <- veilig("Cox overleving", function() overleving_cox_trend(d, g))
  km <- veilig("Kaplan-Meier drie urgentiegroepen", function() overleving(d, g))
  km_es <- veilig("Kaplan-Meier electief vs spoed/redding", function() overleving(d_es, g))
  perf <- veilig("perfusie en temperatuur", function() perfusie_analyse(d, g))
  los <- veilig("opnameduur", function() opnameduur_analyse(d, g))
  variatie <- veilig("centrumverschillen perfusie", function() centrum_variatie(d, g))
  oe <- veilig("O/E per urgentiegroep", function() oe_ratio(d, g))
  es_mod <- veilig("directe E vs S uitkomstmodellen", function() vergelijk_urgentie_es(d_es))
  es_temp <- veilig("directe E vs S temperatuurmodellen", function() vergelijk_urgentie_temp(d_es))
  periodetab <- veilig("periode per urgentiegroep", function() periodevergelijking(d, g))

  # Functioneel gescheiden uitvoer. Alle tabnamen zijn kort en uniek.
  onderdelen <- c(
    list("Aantallen" = aantallen,
         "Tabel 1 E vs S" = tab_es[[1]], "Tabel 2 E vs S" = tab_es[[2]],
         "Tabel 3 E vs S" = tab_es[[3]],
         "Tabel 1 urgent apart" = tab_urgent[[1]],
         "Tabel 2 urgent apart" = tab_urgent[[2]],
         "Tabel 3 urgent apart" = tab_urgent[[3]],
         "Outcome trends E U S" = ut,
         "Uitkomsten E vs S" = es_mod,
         "Temperatuur E vs S" = es_temp,
         "Periode modellen E U S" = per_mod,
         "Cox E U S" = cox,
         "Kaplan Meier E U S" = if (!is.null(km)) km$tabel else NULL,
         "Kaplan Meier E vs S" = if (!is.null(km_es)) km_es$tabel else NULL,
         "Temperatuur modellen" = if (!is.null(perf)) perf$modellen else NULL,
         "Temperatuur centrum" = if (!is.null(perf)) perf$centra else NULL,
         "Temperatuur categorieen" = if (!is.null(perf)) perf$categorieen[["alle procedures"]] else NULL,
         "Temp missing per jaar" = if (!is.null(perf)) perf$ontbrekend_jaar else NULL,
         "Temp missing per centrum" = if (!is.null(perf)) perf$ontbrekend_centrum else NULL,
         "Opnameduur E U S" = los,
         "Perfusie per centrum" = if (!is.null(variatie)) variatie$per_centrum else NULL,
         "Perfusie bereik" = if (!is.null(variatie)) variatie$bereik else NULL,
         "Perfusie ICC" = if (!is.null(variatie)) variatie$icc else NULL,
         "Temp arrest vergelijking" = if (!is.null(variatie)) variatie$zonder_met else NULL,
         "OE per urgentie" = oe,
         "Trend E U S" = if (!is.null(tr3)) tr3$toetsen else NULL,
         "Trend E vs S" = if (!is.null(tr_es)) tr_es$toetsen else NULL,
         "Aantallen per jaar" = if (!is.null(tr3)) tr3$aantal else NULL,
         "Percentages per jaar" = if (!is.null(tr3)) tr3$pct else NULL,
         "Temperatuur per jaar" = if (!is.null(tr3)) tr3$temp else NULL),
    setNames(periodetab, paste0("Periode ", seq_along(periodetab), " ",
                                if (!is.null(periodetab)) gsub("/", "-", names(periodetab), fixed = TRUE) else ""))
  )
  onderdelen <- Filter(function(x) !is.null(x) && is.data.frame(x), onderdelen)
  names(onderdelen) <- substr(names(onderdelen), 1, 31)
  if (anyDuplicated(names(onderdelen))) stop("Dubbele tabbladnamen in urgentiegevoeligheid")
  writexl::write_xlsx(onderdelen, file.path(map, "gevoeligheid_urgentie.xlsx"))

  # Figuur krijgt aparte map zodat primaire Kaplan-Meier behouden blijft.
  if (!is.null(km_es)) {
    sens_fig_map <- file.path(map, "figuren_urgentie")
    dir.create(sens_fig_map, showWarnings = FALSE, recursive = TRUE)
    veilig("KM-figuur E vs S", function() fig_kaplan_meier(km_es, sens_fig_map))
  }
  list(aantallen = aantallen, onderdelen = onderdelen)
}

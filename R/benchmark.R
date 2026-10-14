# =============================================================================
# Benchmark per centrum: funnel plots
#
# Apart voor acuut en electief. Per centrum wordt een uitkomst uitgezet tegen het
# volume (of het verwachte aantal events), met 95%- en 99,8%-controlegrenzen.
#
# - Mortaliteit (ziekenhuis, 30 dagen): risico-gecorrigeerd met EuroSCORE II
#   (indirecte standaardisatie, O/E). Omdat EuroSCORE II landelijk miscalibreert,
#   wordt de verwachte sterfte standaard geherkalibreerd naar het landelijke
#   niveau (FUNNEL_HERKALIBREREN): O/E = 1 betekent dan "gelijk aan het landelijke
#   niveau bij deze case-mix". Grenzen: exacte Poisson-kwantielen.
# - Major complicaties en failure to rescue: ruwe percentages tegen het landelijke
#   percentage van de groep. Grenzen: exacte binomiale kwantielen.
# =============================================================================

funnel_grenzen <- function(type, x, ref) {
  # x = verwacht aantal events (O/E) of aantal patiënten (proportie)
  p <- c(0.001, 0.025, 0.975, 0.999)
  if (type == "OE") sapply(p, function(q) qpois(q, x) / x) + 0
  else sapply(p, function(q) qbinom(q, x, ref) / x) + 0
}

benchmark_centra <- function(df, groep = GROEP) {
  arrest_tijd <- if (OVERLEVING_TIJD %in% names(df)) df[[OVERLEVING_TIJD]] else rep(NA_real_, nrow(df))
  dood <- function(var, min_fu) {
    y <- code_positief(df[[var]])
    mortaliteit_fu(y, arrest_tijd, df[[TREND_JAAR]], min_fu)
  }
  major <- samengesteld(df, MAJOR_COMPLICATIES)
  ftr <- code_positief(df[[FTR_UITKOMST]]); ftr[!(major %in% TRUE)] <- NA
  exp_risico <- if ("EuroscoreII" %in% names(df)) df$EuroscoreII / 100 else rep(NA, nrow(df))

  uitkomsten <- list(
    list(label = "Ziekenhuismortaliteit (O/E, EuroSCORE II)", type = "OE", y = dood("overleden_zh", 0)),
    list(label = "30-dagenmortaliteit (O/E, EuroSCORE II)", type = "OE", y = dood("overl30d", 30)),
    list(label = "Major complicaties (ruw)", type = "prop", y = major),
    list(label = "Failure to rescue (ruw)", type = "prop", y = ftr)
  )

  basis <- data.frame(groep = as.character(df[[groep]]), centrum = as.character(df[[CENTRUM_KOLOM]]),
                      e = exp_risico)
  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (u in uitkomsten) {
    d <- basis; d$y <- u$y
    d <- d[d$groep %in% g & !is.na(d$y) & !is.na(d$centrum), ]
    if (u$type == "OE") d <- d[!is.na(d$e), ]
    if (nrow(d) < 30) next
    O_tot <- sum(d$y); n_tot <- nrow(d)
    if (u$type == "OE") {
      kal <- if (FUNNEL_HERKALIBREREN) O_tot / sum(d$e) else 1
      d$e <- d$e * kal
    }
    ref <- O_tot / n_tot
    per <- do.call(rbind, lapply(split(d, d$centrum), function(x) {
      O <- sum(x$y); n <- nrow(x)
      if (u$type == "OE") { E <- sum(x$e); waarde <- O / E; gr <- funnel_grenzen("OE", E, ref) }
      else { E <- n * ref; waarde <- O / n; gr <- funnel_grenzen("prop", n, ref) }
      status <- if (waarde > gr[4]) "boven 99,8%-grens" else if (waarde > gr[3]) "boven 95%-grens" else
        if (waarde < gr[1]) "onder 99,8%-grens" else if (waarde < gr[2]) "onder 95%-grens" else "binnen grenzen"
      data.frame(groep = g, uitkomst = u$label, centrum = x$centrum[1], n = n, geobserveerd = O,
                 verwacht = round(E, 1),
                 waarde = round(if (u$type == "OE") waarde else 100 * waarde, 2),
                 eenheid = if (u$type == "OE") "O/E" else "%",
                 grens95 = sprintf("%.2f-%.2f", 0 + (if (u$type == "OE") 1 else 100) * gr[2], (if (u$type == "OE") 1 else 100) * gr[3]),
                 grens998 = sprintf("%.2f-%.2f", (if (u$type == "OE") 1 else 100) * gr[1], (if (u$type == "OE") 1 else 100) * gr[4]),
                 status = status)
    }))
    attr_ruw <- if (u$type == "OE") sprintf("landelijk O/E t.o.v. EuroSCORE II zonder herkalibratie: %.2f", O_tot / sum(d$e / kal)) else
      sprintf("landelijk: %.1f%%", 100 * ref)
    per$landelijk <- attr_ruw
    res[[length(res) + 1]] <- per
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out[order(out$groep, out$uitkomst, suppressWarnings(as.numeric(out$centrum))), ]
}

fig_funnels <- function(bm, map) {
  if (is.null(bm) || !nrow(bm)) return(NULL)
  paden <- c()
  for (g in unique(bm$groep)) {
    d <- bm[bm$groep == g, ]
    lijnen <- do.call(rbind, lapply(unique(d$uitkomst), function(u) {
      x <- d[d$uitkomst == u, ]
      oe <- x$eenheid[1] == "O/E"
      ref <- if (oe) 1 else sum(x$geobserveerd) / sum(x$n)
      as <- if (oe) x$verwacht else x$n
      grid <- seq(max(1, min(as) * 0.8), max(as) * 1.1, length.out = 200)
      gr <- funnel_grenzen(if (oe) "OE" else "prop", if (oe) grid else round(grid), ref)
      f <- if (oe) 1 else 100
      data.frame(uitkomst = u, x = grid, ref = f * ref, o998 = f * gr[, 1], o95 = f * gr[, 2],
                 b95 = f * gr[, 3], b998 = f * gr[, 4])
    }))
    d$x <- ifelse(d$eenheid == "O/E", d$verwacht, d$n)
    d$status_kort <- factor(ifelse(grepl("99,8", d$status), "buiten 99,8%",
                                   ifelse(grepl("95", d$status), "buiten 95%", "binnen grenzen")),
                            levels = c("binnen grenzen", "buiten 95%", "buiten 99,8%"))
    p <- ggplot() +
      geom_line(data = lijnen, aes(x, o998), linetype = "dotted", colour = "grey40") +
      geom_line(data = lijnen, aes(x, b998), linetype = "dotted", colour = "grey40") +
      geom_line(data = lijnen, aes(x, o95), linetype = "dashed", colour = "grey55") +
      geom_line(data = lijnen, aes(x, b95), linetype = "dashed", colour = "grey55") +
      geom_line(data = lijnen, aes(x, ref), colour = "grey20") +
      geom_point(data = d, aes(x, waarde, colour = status_kort), size = 2.5) +
      geom_text(data = d, aes(x, waarde, label = centrum), size = 2.6, vjust = -0.9) +
      scale_colour_manual(values = c("binnen grenzen" = KLEUR, "buiten 95%" = "#E67E22",
                                     "buiten 99,8%" = "#C0392B"), drop = FALSE) +
      facet_wrap(~uitkomst, scales = "free") +
      labs(title = paste("Benchmark per centrum -", g),
           subtitle = paste0("O/E: x = verwacht aantal (EuroSCORE II", if (FUNNEL_HERKALIBREREN) ", geherkalibreerd naar landelijk niveau",
                             "); %: x = aantal patiënten\nstreep = 95%-grens, stippel = 99,8%-grens; label = centrum"),
           x = "Verwacht aantal events / aantal patiënten", y = "O/E of %", colour = NULL) + thema
    paden <- c(paden, bewaar(p, map, paste0("funnel_", gsub("[^a-z]+", "_", tolower(g))), 12, 9))
  }
  paden
}

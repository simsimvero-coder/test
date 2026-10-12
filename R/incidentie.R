# =============================================================================
# Operatievolume per jaar en incidentie per miljoen inwoners
#
# - Aantallen per jaar: totaal, acuut, electief en niet ingedeeld
# - Operaties per miljoen inwoners (CBS-bevolking uit config.R)
# - Poisson-regressie met log(bevolking) als offset: rate ratio (RR) per jaar
#   voor het totaal (alle jaren) en per groep (vanaf GROEPEN_VANAF_JAAR)
# =============================================================================

bevolking_voor <- function(jaren) {
  b <- BEVOLKING_1JAN
  sapply(jaren, function(j) {
    nu <- b[as.character(j)]; volgend <- b[as.character(j + 1)]
    if (is.na(nu)) return(NA_real_)
    if (BEVOLKING_GEMIDDELD && !is.na(volgend)) (nu + volgend) / 2 else unname(nu)
  })
}

poisson_rr <- function(n, jaar, bevolking, label, periode) {
  ok <- !is.na(bevolking)
  d <- data.frame(n = n[ok], jaar = jaar[ok], bev = bevolking[ok])
  leeg <- data.frame(maat = label, periode = periode, n_totaal = sum(d$n),
                     n_eerste_jaar = NA, n_laatste_jaar = NA, verandering_pct = NA,
                     per_miljoen_eerste = NA, per_miljoen_laatste = NA,
                     RR_per_jaar = NA, BI95_onder = NA, BI95_boven = NA, p = "")
  if (nrow(d) < 2 || sum(d$n) == 0) return(leeg)
  m <- glm(n ~ jaar + offset(log(bev)), family = poisson, data = d)
  s <- coef(summary(m))["jaar", ]
  ci <- exp(s[1] + c(-1.96, 1.96) * s[2])
  eerste <- d[which.min(d$jaar), ]; laatste <- d[which.max(d$jaar), ]
  leeg$n_eerste_jaar <- eerste$n; leeg$n_laatste_jaar <- laatste$n
  leeg$verandering_pct <- round(100 * (laatste$n - eerste$n) / eerste$n, 1)
  leeg$per_miljoen_eerste <- round(1e6 * eerste$n / eerste$bev, 1)
  leeg$per_miljoen_laatste <- round(1e6 * laatste$n / laatste$bev, 1)
  leeg$RR_per_jaar <- round(exp(s[1]), 3)
  leeg$BI95_onder <- round(ci[1], 3); leeg$BI95_boven <- round(ci[2], 3)
  leeg$p <- fmt_p(s[4])
  leeg
}

# df_alle: alle patiënten, dus vóór het uitsluiten van niet-ingedeelde patiënten
incidentie <- function(df_alle, groep = GROEP, jaar_kolom = TREND_JAAR) {
  jaar <- suppressWarnings(as.numeric(as.character(df_alle[[jaar_kolom]])))
  g <- if (!is.null(groep) && groep %in% names(df_alle)) as.character(df_alle[[groep]]) else rep(NA, nrow(df_alle))
  g[is.na(g)] <- "niet ingedeeld"
  jaren <- seq(min(jaar, na.rm = TRUE), max(jaar, na.rm = TRUE))
  niveaus <- c(setdiff(unique(c(unname(ACUUT_ELECTIEF$labels), g)), "niet ingedeeld"), "niet ingedeeld")

  tel <- table(factor(jaar, levels = jaren), factor(g, levels = niveaus))
  per_jaar <- data.frame(jaar = jaren, totaal = as.integer(rowSums(tel)), check.names = FALSE)
  for (niv in niveaus) per_jaar[[niv]] <- as.integer(tel[, niv])
  per_jaar$bevolking <- round(bevolking_voor(jaren))
  per_jaar$totaal_per_miljoen <- round(1e6 * per_jaar$totaal / per_jaar$bevolking, 1)
  for (niv in setdiff(niveaus, "niet ingedeeld"))
    per_jaar[[paste0(niv, "_per_miljoen")]] <- round(1e6 * per_jaar[[niv]] / per_jaar$bevolking, 1)
  if (any(is.na(per_jaar$bevolking)))
    message("LET OP - geen bevolkingscijfer voor: ",
            paste(per_jaar$jaar[is.na(per_jaar$bevolking)], collapse = ", "),
            " (vul BEVOLKING_1JAN in config.R aan)")
  if (sum(is.na(jaar)))
    message("LET OP - ", sum(is.na(jaar)), " patiënten zonder jaar tellen niet mee in het volume per jaar")

  vanaf <- per_jaar$jaar >= GROEPEN_VANAF_JAAR
  periode_alle <- sprintf("%d-%d", min(jaren), max(jaren))
  periode_groep <- sprintf("%d-%d", max(min(jaren), GROEPEN_VANAF_JAAR), max(jaren))
  rr <- rbind(
    poisson_rr(per_jaar$totaal, per_jaar$jaar, per_jaar$bevolking, "Totaal (incl. niet ingedeeld)", periode_alle),
    if (periode_groep != periode_alle)
      poisson_rr(per_jaar$totaal[vanaf], per_jaar$jaar[vanaf], per_jaar$bevolking[vanaf],
                 "Totaal (incl. niet ingedeeld)", periode_groep),
    do.call(rbind, lapply(setdiff(niveaus, "niet ingedeeld"), function(niv)
      poisson_rr(per_jaar[[niv]][vanaf], per_jaar$jaar[vanaf], per_jaar$bevolking[vanaf], niv, periode_groep)))
  )
  list(per_jaar = per_jaar, rr = rr)
}

fig_incidentie <- function(inc, map) {
  pj <- inc$per_jaar
  kolommen <- grep("_per_miljoen$", names(pj), value = TRUE)
  d <- do.call(rbind, lapply(kolommen, function(k) {
    x <- data.frame(jaar = pj$jaar, per_miljoen = pj[[k]], reeks = sub("_per_miljoen$", "", k))
    if (k != "totaal_per_miljoen") x <- x[x$jaar >= GROEPEN_VANAF_JAAR, ]
    x
  }))
  d <- d[!is.na(d$per_miljoen), ]
  rr <- inc$rr[inc$rr$p != "", ]
  ondertitel <- paste(sprintf("%s (%s): RR per jaar %.3f (95%%-BI %.3f-%.3f), p = %s",
                              rr$maat, rr$periode, rr$RR_per_jaar, rr$BI95_onder, rr$BI95_boven, rr$p),
                      collapse = "\n")
  p <- ggplot(d, aes(jaar, per_miljoen, colour = reeks)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2) +
    scale_colour_manual(values = c("totaal" = "grey30", "acuut" = "#C0392B",
                                   "electief aneurysma" = "#2E4F82")) +
    scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 1)) +
    expand_limits(y = 0) +
    labs(title = "Aortaboogoperaties per miljoen inwoners", subtitle = ondertitel,
         x = NULL, y = "Operaties per miljoen inwoners", colour = NULL) + thema
  bewaar(p, map, "volume_per_miljoen", 10, 6.5)
}

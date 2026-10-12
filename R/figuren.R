# =============================================================================
# Figuren (ggplot2), opgeslagen als PNG
# =============================================================================

library(ggplot2)

thema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom")
KLEUR <- "#2E4F82"

bewaar <- function(plot, map, naam, breed = 10, hoog = 6) {
  pad <- file.path(map, paste0(naam, ".png"))
  ggsave(pad, plot, width = breed, height = hoog, dpi = 150, bg = "white")
  pad
}

lang_numeriek <- function(df, vars) {
  vars <- intersect(vars, names(df)[sapply(df, is.numeric)])
  if (!length(vars)) return(NULL)
  tidyr::pivot_longer(df[vars], dplyr::everything(), names_to = "variabele", values_to = "waarde") |>
    dplyr::filter(!is.na(waarde))
}

fig_histogrammen <- function(df, map) {
  d <- lang_numeriek(df, setdiff(names(NUMERIEK), OVERLEVING_TIJD))
  if (is.null(d)) return(NULL)
  p <- ggplot(d, aes(waarde)) +
    geom_histogram(bins = 30, fill = KLEUR, colour = "white") +
    facet_wrap(~variabele, scales = "free") +
    labs(title = "Verdeling numerieke variabelen", x = NULL, y = "Aantal patiënten") + thema
  bewaar(p, map, "histogrammen", 12, 8)
}

fig_boxplots <- function(df, groep, map) {
  if (is.null(groep) || !groep %in% names(df)) return(NULL)
  vars <- intersect(c("leeftijd", "BMI", "eGFR", "LVEF", "EuroscoreII", "ECC_duur",
                      "AoX_duur", "opnameduur"), names(df)[sapply(df, is.numeric)])
  d <- df[!is.na(df[[groep]]), c(groep, vars)]
  labels <- sapply(vars, function(v) {
    t <- toets_num(d[[v]], d[[groep]], is_normaal(d[[v]]))
    sprintf("%s\n%s p = %s", v, t$toets, fmt_p(t$p))
  })
  d <- tidyr::pivot_longer(d, -dplyr::all_of(groep), names_to = "variabele", values_to = "waarde")
  d$variabele <- factor(labels[d$variabele], levels = labels)
  p <- ggplot(d[!is.na(d$waarde), ], aes(.data[[groep]], waarde, fill = .data[[groep]])) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
    facet_wrap(~variabele, scales = "free_y") +
    labs(title = paste("Kenmerken per", groep), x = NULL, y = NULL, fill = NULL) +
    thema + theme(axis.text.x = element_blank())
  bewaar(p, map, paste0("boxplots_per_", groep), 12, 8)
}

fig_complicaties <- function(df, groep, map) {
  vars <- intersect(TABELLEN[[3]], names(df))
  vars <- vars[sapply(df[vars], function(x) is.factor(x) && "ja" %in% levels(x))]
  if (!length(vars)) return(NULL)
  g <- if (!is.null(groep) && groep %in% names(df)) as.character(df[[groep]]) else rep("Totaal", nrow(df))
  d <- do.call(rbind, lapply(vars, function(v) {
    sel <- !is.na(df[[v]]) & !is.na(g)
    agg <- tapply(df[[v]][sel] == "ja", g[sel], mean)
    data.frame(variabele = v, groep = names(agg), pct = 100 * as.numeric(agg))
  }))
  labels <- sapply(vars, function(v) {
    if (all(g == "Totaal")) return(v)
    sprintf("%s (p = %s)", v, fmt_p(toets_cat(df[[v]], g)$p))
  })
  d$variabele <- factor(labels[d$variabele], levels = rev(labels))
  p <- ggplot(d, aes(pct, variabele, fill = groep)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    labs(title = "Postoperatieve uitkomsten", x = "% van patiënten", y = NULL,
         fill = if (is.null(groep)) NULL else groep) +
    guides(fill = guide_legend(reverse = TRUE)) + thema
  bewaar(p, map, "postoperatieve_uitkomsten", 10, 0.45 * length(vars) * max(1, length(unique(g)) / 2) + 2)
}

fig_kaplan_meier <- function(km, map) {
  if (is.null(km)) return(NULL)
  s <- km$fit
  d <- data.frame(tijd = s$time, overleving = s$surv, onder = s$lower, boven = s$upper,
                  groep = if (km$heeft_groep) rep(sub("^\\.groep=", "", names(s$strata)), s$strata) else "Totaal")
  # startpunt (t = 0, overleving 100%) per groep toevoegen
  start <- data.frame(tijd = 0, overleving = 1, onder = 1, boven = 1, groep = unique(d$groep))
  d <- rbind(start, d)
  d <- d[d$tijd <= OVERLEVING_MAX_DAGEN, ]
  d <- d[order(d$groep, d$tijd), ]
  d$tijd_eind <- ave(d$tijd, d$groep, FUN = function(t) c(t[-1], OVERLEVING_MAX_DAGEN))
  p <- ggplot(d, aes(tijd / 365.25, 100 * overleving, colour = groep, fill = groep)) +
    geom_step(linewidth = 0.8) +
    geom_rect(aes(xmin = tijd / 365.25, xmax = tijd_eind / 365.25,
                  ymin = 100 * onder, ymax = 100 * boven), alpha = 0.12, colour = NA) +
    labs(title = "Kaplan-Meier overleving",
         subtitle = if (km$logrank_p != "") paste("log-rank p =", km$logrank_p),
         x = "Jaren na operatie", y = "Overleving (%)", colour = NULL, fill = NULL) +
    coord_cartesian(ylim = c(max(0, 100 * min(d$onder, na.rm = TRUE) - 5), 100)) + thema
  bewaar(p, map, "kaplan_meier", 9, 6)
}

fig_oe <- function(oe, map) {
  if (is.null(oe)) return(NULL)
  per <- names(oe)[1]
  oe[[per]] <- sprintf("%s (p = %s)", oe[[per]], oe$p)
  oe[[per]] <- factor(oe[[per]], levels = rev(oe[[per]]))
  p_groep <- attr(oe, "p_tussen_groepen")
  p <- ggplot(oe, aes(OE_ratio, .data[[per]])) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = BI95_onder, xmax = BI95_boven), height = 0.25, colour = KLEUR) +
    geom_point(aes(size = n), colour = KLEUR) +
    labs(title = "Geobserveerde / verwachte mortaliteit (EuroSCORE II)",
         subtitle = paste0("p per rij: O/E verschilt van 1 (exacte Poisson-toets)",
                           if (!is.null(p_groep)) paste0("\nverschil tussen ", per, " (gecorrigeerd voor EuroSCORE II): p = ",
                                                         fmt_p(p_groep))),
         x = "O/E-ratio (95%-BI)", y = NULL, size = "n") + thema
  bewaar(p, map, paste0("OE_ratio_per_", per), 8, 0.4 * nrow(oe) + 2)
}

fig_ontbrekend <- function(ontbr, map) {
  d <- ontbr[ontbr$ontbrekend_n > 0, ]
  if (!nrow(d)) return(NULL)
  d$variabele <- factor(d$variabele, levels = rev(d$variabele))
  p <- ggplot(d, aes(ontbrekend_pct, variabele)) + geom_col(fill = KLEUR) +
    labs(title = "Ontbrekende waarden", x = "% ontbrekend", y = NULL) + thema
  bewaar(p, map, "ontbrekende_waarden", 8, 0.25 * nrow(d) + 1.5)
}

# =============================================================================
# 17_figures_paper1.R
# All figures for paper 1, drawn directly from the pipeline output (no numbers
# typed in): NEJM-style, green accents for acute (dark) and non-acute (light).
# Needs 02, 03, 07 and 16 to have run. Each figure is written as PNG and TIFF
# (600 dpi), PDF and SVG to results/<method>/figures_paper1/.
#   Figure 1  Flow diagram                       Figure 6  Landmark survival at day 30
#   Figure 2  Annual volume and rate             Figure 7  Standardised risk by year and group
#   Figure 3  Acute vs non-acute outcomes         Figure 8  Perfusion practice by year and group
#   Figure 4  Year trends by group (+ S1 2015)    Figure 9  Practice per centre (no outcomes)
#   Figure 5  Survival vs general population      Figure 10 Complications and failure to rescue
# =============================================================================
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(survival) })
SUB <- "figures_paper1"
OUT <- out_dir(SUB)
rd <- function(sub, f) { p <- file.path(res_path(sub), f); if (!file.exists(p)) stop("missing input: ", p); fread(p, check.names = FALSE) }
FONT <- if (.Platform$OS.type == "windows") "sans" else "Liberation Sans"
ink <- "#1A1A1A"; grid_col <- "#E4E4E4"; shade <- "#F1F1EE"
GC <- c("Acute" = "#1B5E4B", "Non-acute" = "#5FAE88")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
fp <- function(p) vapply(p, function(v) if (is.na(v)) "" else if (v < .001) "<.001" else
  sub("^0", "", formatC(v, format = "f", digits = if (v < .01) 3 else 2)), "")
f2 <- function(x) formatC(x, format = "f", digits = 2)
th <- theme_classic(base_size = 8.5, base_family = FONT) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink, size = 7.5),
        axis.line = element_line(linewidth = .35, colour = ink), axis.ticks = element_line(linewidth = .35, colour = ink),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = .3), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 11), plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
        plot.margin = margin(6, 8, 4, 6))
save_all <- function(p, name, w, h) {
  f <- file.path(OUT, name)
  ggsave(paste0(f, ".png"), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(paste0(f, ".tiff"), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(paste0(f, ".svg"), p, width = w, height = h, units = "mm", device = grDevices::svg)
  log_msg("  saved ", name)
}
fig <- function(name, expr) tryCatch(expr, error = function(e) log_msg("FIGURE SKIPPED (", name, "): ", conditionMessage(e)))
YRS <- CONFIG$study_years
YR <- scale_x_continuous(breaks = YRS, labels = function(x) ifelse(x %in% YRS[seq(1, length(YRS), 3)], x, ""), expand = c(.03, .03))
early <- annotate("rect", xmin = min(YRS) - .4, xmax = CONFIG$trend_sens_start - .5, ymin = -Inf, ymax = Inf, fill = shade)
key_groups <- function(extra = NULL) {
  p <- ggplot() + xlim(0, 10) + ylim(0, 1) +
    annotate("segment", x = c(.2, 2.6), xend = c(.8, 3.2), y = .5, yend = .5, colour = GC[c("Acute", "Non-acute")], linewidth = .9) +
    annotate("point", x = c(.5, 2.9), y = .5, colour = GC[c("Acute", "Non-acute")], size = 1.6) +
    annotate("text", x = c(1, 3.4), y = .5, label = c("Acute", "Non-acute"), hjust = 0, size = 2.7, family = FONT) +
    annotate("rect", xmin = 5, xmax = 5.5, ymin = .3, ymax = .7, fill = shade) +
    annotate("text", x = 5.7, y = .5, label = sprintf("%d–%d: incomplete classification", min(YRS), CONFIG$trend_sens_start - 1),
             hjust = 0, size = 2.5, family = FONT, colour = "#666666")
  (if (is.null(extra)) p else p + extra) + theme_void() }
log_msg("Start 17_figures_paper1")
cohort <- load_derived("cohort.rds")

# ---- Figure 1: flow diagram ------------------------------------------------------------
fig("1", {
  fl <- rd("02_descriptives", "flow_diagram.csv")
  nice <- c("year missing" = "Year of surgery missing or outside study period", "no cardiopulmonary bypass" = "No extracorporeal circulation",
            "endovascular" = "Endovascular procedure", "percutaneous" = "Percutaneous approach",
            "no aortic surgery" = "No aortic surgery recorded", "approach unknown" = "Aortic approach unknown",
            "chronic dissection" = "Chronic dissection", "infection" = "Aortic infection", "duplicate" = "Duplicate record",
            "upstream" = "Excluded before selection")
  lab_of <- function(s) { k <- names(nice)[vapply(names(nice), function(z) grepl(z, s, fixed = TRUE), logical(1))]; if (length(k)) nice[[k[1]]] else trimws(sub("excluded:", "", s)) }
  n_rec <- fl$n[1]
  ex <- fl[grepl("excluded", step) & !grepl("not classifiable", step) & n != 0][, .(lab = vapply(step, lab_of, ""), n = -n)]
  n_surg <- fl[grepl("^Surgical cohort", step), n]
  n_ncl <- -fl[grepl("not classifiable", step), n]
  n_an <- fl[grepl("^Analysis cohort", step), n]
  sc <- load_derived("surgical_cohort.rds")
  ncl_y <- sc[is.na(acute), .N, by = year][order(year)]
  grp <- cohort[, .N, by = urgency2]; nsv <- cohort[!is.na(time_d) & !is.na(status), .N, by = urgency2]
  g_n <- function(t, g) t[urgency2 == g, N]
  fmt <- function(x) format(x, big.mark = ",")
  box <- function(x, y, w, h, fill = "white", col = ink, lw = .4) annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill, colour = col, linewidth = lw)
  tx <- function(x, y, lab, size = 2.9, face = "plain", hjust = .5, col = ink) annotate("text", x = x, y = y, label = lab, size = size, fontface = face, hjust = hjust, family = FONT, lineheight = 1.15, colour = col)
  arr <- function(x, y, xe, ye) annotate("segment", x = x, y = y, xend = xe, yend = ye, linewidth = .4, colour = ink, arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
  ln <- function(x, y, xe, ye) annotate("segment", x = x, y = y, xend = xe, yend = ye, linewidth = .4, colour = ink)
  hx <- max(11, 5 + 3.2 * nrow(ex))
  p <- ggplot() + xlim(0, 100) + ylim(0, 100) +
    box(32, 93, 46, 9) + tx(32, 94.6, sprintf("Open aortic arch procedures registered in the NHR,\n%d centres, %d–%d", uniqueN(cohort$centrum), min(YRS), max(YRS)), face = "bold") +
    tx(32, 89.8, paste("n =", fmt(n_rec), "records"), size = 2.8) +
    box(76, 73, 44, hx * .9, fill = shade) + tx(56.5, 73 + hx * .45 - 3, sprintf("Excluded (n = %s)", fmt(sum(ex$n))), face = "bold", hjust = 0) +
    tx(56.5, 73 - 2, paste(sprintf("%s (n = %s)", ex$lab, fmt(ex$n)), collapse = "\n"), size = 2.6, hjust = 0) +
    box(32, 54, 46, 8) + tx(32, 55.4, "Open aortic arch procedures", face = "bold") + tx(32, 51.6, paste("n =", fmt(n_surg)), size = 2.8) +
    box(76, 41, 44, 11, fill = shade) + tx(56.5, 43.2, sprintf("Not classifiable as acute or non-acute (n = %s)", fmt(n_ncl)), face = "bold", hjust = 0, size = 2.6) +
    tx(56.5, 38.6, if (nrow(ncl_y)) paste(sprintf("%d, n = %s", ncl_y$year, fmt(ncl_y$N)), collapse = "; ") else "", size = 2.6, hjust = 0) +
    box(32, 30, 46, 8, lw = .6) + tx(32, 31.4, "Analysis cohort", face = "bold") + tx(32, 27.6, paste("n =", fmt(n_an)), size = 2.8) +
    box(16, 12, 26, 11, fill = GC[["Non-acute"]], col = GC[["Non-acute"]]) + tx(16, 14, "Non-acute", face = "bold") + tx(16, 9.6, paste("n =", fmt(g_n(grp, "Non-acute"))), size = 2.8) +
    box(48, 12, 26, 11, fill = GC[["Acute"]], col = GC[["Acute"]]) + tx(48, 14, "Acute", face = "bold", col = "white") + tx(48, 9.6, paste("n =", fmt(g_n(grp, "Acute"))), size = 2.8, col = "white") +
    tx(16, 3, paste0("Survival follow-up available,\nn = ", fmt(g_n(nsv, "Non-acute"))), size = 2.4, col = "#555555") +
    tx(48, 3, paste0("Survival follow-up available,\nn = ", fmt(g_n(nsv, "Acute"))), size = 2.4, col = "#555555") +
    arr(32, 88.5, 32, 58) + ln(32, 73, 54, 73) + arr(32, 50, 32, 34) + ln(32, 41, 54, 41) +
    ln(32, 26, 32, 21.5) + ln(16, 21.5, 48, 21.5) + arr(16, 21.5, 16, 17.5) + arr(48, 21.5, 48, 17.5) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  save_all(p, "Figure01_flowchart", 180, 150)
})

# ---- Figure 2: annual volume ---------------------------------------------------------------
fig("2", {
  ac <- rd("03_volume", "annual_counts.csv"); rt <- rd("03_volume", "annual_rates_per_million.csv"); vm <- rd("03_volume", "volume_trend_models.csv")
  lab_rr <- function(s, a) { r <- vm[series == s & analysis == a][1]
    sprintf("RR %.3f (%.3f–%.3f)%s", r$irr_per_year, r$irr_lower_dispersion_scaled, r$irr_upper_dispersion_scaled, stars(r$p_trend_dispersion_scaled)) }
  c_na <- "#93C9AE"; c_nc <- "#DADAD5"
  XL <- c(min(ac$year) - .45, max(ac$year) + 3.3)
  ac[, `:=`(top_a = Acute, top_n = Acute + `Non-acute`, top_t = All)]
  g <- seq(min(ac$year), max(ac$year), length.out = 500)
  sm <- function(y) splinefun(ac$year, y, method = "monoH.FC")(g)
  bands <- data.table(x = g, a = sm(ac$top_a), n = sm(ac$top_n), t = sm(ac$top_t))
  ymax <- ceiling(max(ac$All) * 1.15 / 100) * 100
  lab <- rbind(ac[, .(year, y = Acute / 2, n = Acute, col = "white", face = "bold")],
               ac[, .(year, y = Acute + `Non-acute` / 2, n = `Non-acute`, col = "#174A3C", face = "plain")],
               ac[`Not classifiable` >= 15, .(year, y = Acute + `Non-acute` + `Not classifiable` / 2, n = `Not classifiable`, col = "#555555", face = "plain")])
  lab[, xl := fifelse(year == min(year), year + .3, fifelse(year == max(year), year - .3, as.numeric(year)))]
  leg <- data.table(y = ymax * c(.93, .86, .79), lab = c("Acute", "Non-acute", "Not classifiable"), fill = c(GC[["Acute"]], c_na, c_nc))
  pA <- ggplot() +
    geom_ribbon(data = bands, aes(x, ymin = n, ymax = t), fill = c_nc) + geom_ribbon(data = bands, aes(x, ymin = a, ymax = n), fill = c_na) +
    geom_ribbon(data = bands, aes(x, ymin = 0, ymax = a), fill = GC[["Acute"]]) +
    geom_line(data = bands, aes(x, a), colour = "white", linewidth = .5) + geom_line(data = bands, aes(x, n), colour = "white", linewidth = .5) +
    geom_text(data = lab, aes(xl, y, label = n), colour = lab$col, fontface = lab$face, size = 2.6, family = FONT) +
    geom_text(data = ac, aes(year, All, label = All), vjust = -.8, fontface = "bold", size = 2.8, family = FONT, colour = ink) +
    geom_tile(data = leg, aes(x = max(ac$year) + 1.15, y = y, fill = fill), width = .38, height = ymax * .045) + scale_fill_identity() +
    geom_text(data = leg, aes(x = max(ac$year) + 1.5, y = y, label = lab), hjust = 0, size = 3, family = FONT, colour = ink) +
    annotate("text", x = max(ac$year) + .96, y = ymax * .70, label = "Bold: annual total", hjust = 0, size = 2.6, family = FONT, colour = "#666666", fontface = "italic") +
    scale_x_continuous(breaks = ac$year, limits = XL, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 100), expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Procedures per year (n)", tag = "A") + th + theme(panel.grid.major.y = element_blank(), axis.line.x = element_blank())
  rl <- rbind(data.table(year = rt$year, series = "All procedures", rate = rt$All),
              data.table(year = rt$year, series = "Acute", rate = rt$Acute)[year >= CONFIG$trend_sens_start],
              data.table(year = rt$year, series = "Non-acute", rate = rt$`Non-acute`)[year >= CONFIG$trend_sens_start])
  endp <- rl[year == max(year)]
  sa <- paste("from", CONFIG$trend_sens_start)
  endp[, txt := c("All procedures" = lab_rr("All", "primary"), "Acute" = lab_rr("Acute", sa), "Non-acute" = lab_rr("Non-acute", sa))[as.character(series)]]
  yr_max <- ceiling(max(rl$rate) * 1.1 / 5) * 5
  setorder(endp, rate); endp[, ylab := rate]
  for (i in seq_len(nrow(endp))[-1]) endp$ylab[i] <- max(endp$ylab[i], endp$ylab[i - 1] + yr_max * .1)
  pB <- ggplot(rl, aes(year, rate, colour = series)) +
    annotate("rect", xmin = XL[1], xmax = CONFIG$trend_sens_start - .5, ymin = 0, ymax = yr_max, fill = shade) +
    annotate("text", x = (XL[1] + CONFIG$trend_sens_start - .5) / 2, y = yr_max * .97, label = "Incomplete\nclassification", vjust = 1, size = 2.5, family = FONT, colour = "#777777", lineheight = .95) +
    geom_line(aes(linetype = series), linewidth = .7) + geom_point(aes(shape = series), size = 1.9, fill = c_na, stroke = .7) +
    geom_text(data = endp, aes(x = max(year) + .55, y = ylab, label = series), hjust = 0, vjust = -.25, fontface = "bold", size = 3, family = FONT) +
    geom_text(data = endp, aes(x = max(year) + .55, y = ylab, label = txt), hjust = 0, vjust = 1.35, size = 2.6, family = FONT) +
    annotate("text", x = max(rl$year) + .55, y = yr_max * .07, label = "RR per year (95% CI)\n*P<.05, **P<.01, ***P<.001", hjust = 0, size = 2.4, family = FONT, colour = "#666666", lineheight = 1) +
    scale_colour_manual(values = c("All procedures" = "#3C3C3C", "Acute" = GC[["Acute"]], "Non-acute" = "#4E9C78")) +
    scale_linetype_manual(values = c("All procedures" = "solid", "Acute" = "solid", "Non-acute" = "22")) +
    scale_shape_manual(values = c("All procedures" = 16, "Acute" = 16, "Non-acute" = 22)) +
    scale_x_continuous(breaks = rt$year, limits = XL, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, yr_max), breaks = seq(0, yr_max, 5), expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = "Year of surgery", y = "Procedures per million inhabitants", tag = "B") + th + theme(panel.grid.major.y = element_blank())
  save_all(pA / pB + plot_layout(heights = c(1.05, 1)), "Figure02_volume", 180, 200)
})

# ---- Figure 3: acute vs non-acute -------------------------------------------------------
OUT_ORDER <- c("comp_primary", "dzh", "d30", "d1y", "stroke_def", "stroke_nodef", "renal", "gi", "vent24", "reintub", "icu_readm", "major_compl")
SHORT <- c(comp_primary = "Composite\noutcome", dzh = "In-hospital\ndeath", d30 = "30-day\ndeath", d1y = "One-year\ndeath",
           stroke_def = "Stroke with\ndeficit", stroke_nodef = "Stroke/TIA\nwithout\ndeficit", renal = "Renal\nfailure",
           gi = "GI\ncompli-\ncation", vent24 = "Ventilation\n>24 h", reintub = "Reintu-\nbation", icu_readm = "ICU\nreadmission",
           major_compl = "Major\ncompli-\ncation")
fig("3", {
  cr <- rd("16_paper1_extra", "E_crude_rates.csv")[period == "All" & group %in% names(GC) & outcome %in% OUT_ORDER]
  cr[, `:=`(outcome = factor(outcome, OUT_ORDER), group = factor(group, c("Non-acute", "Acute")))]
  ymax <- ceiling(max(cr$pct) * 1.12 / 5) * 5
  nlab <- cohort[, .N, by = urgency2]
  pA <- ggplot(cr, aes(outcome, pct, fill = group)) + geom_col(position = position_dodge(.78), width = .74) +
    geom_text(aes(label = formatC(pct, format = "f", digits = 1)), position = position_dodge(.78), vjust = -.45, size = 2.2, family = FONT) +
    scale_fill_manual(values = c("Non-acute" = "#93C9AE", "Acute" = GC[["Acute"]]),
                      labels = sprintf("%s (n = %s)", c("Non-acute", "Acute"), format(c(nlab[urgency2 == "Non-acute", N], nlab[urgency2 == "Acute", N]), big.mark = ","))) +
    scale_x_discrete(labels = SHORT) + scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 5), expand = c(0, 0)) +
    labs(x = NULL, y = "Patients (%)", tag = "A", fill = NULL) + th +
    theme(axis.text.x = element_text(size = 6.6, lineheight = .9), legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = 7.5), legend.key.size = unit(3.5, "mm"))
  orx <- rd("16_paper1_extra", "C_acute_vs_nonacute_or.csv")[outcome %in% c("comp_primary", "dzh", "d30", "d1y", "stroke_def", "renal", "major_compl")]
  orx[, outcome := factor(outcome, OUT_ORDER)]; setorder(orx, outcome); orx[, y := .N:1]
  orx[, lab := sprintf("%s (%s–%s)%s", f2(or_acute), f2(lower), f2(upper), stars(p))]
  xl <- c(min(.42, min(orx$lower) * .8), max(4.6, max(orx$upper) * 1.6))
  pB <- ggplot(orx) + geom_vline(xintercept = 1, linewidth = .35, colour = ink) +
    geom_segment(aes(x = lower, xend = upper, y = y, yend = y), linewidth = .5, colour = GC[["Acute"]]) +
    geom_point(aes(or_acute, y), shape = 22, size = 2.4, fill = GC[["Acute"]], colour = GC[["Acute"]]) +
    geom_text(aes(x = xl[1], y = y, label = label), hjust = 0, size = 2.6, family = FONT) +
    geom_text(aes(x = xl[2], y = y, label = lab), hjust = 1, size = 2.6, family = FONT) +
    annotate("text", x = xl, y = nrow(orx) + .9, label = c("Outcome", "Adjusted OR (95% CI)"), hjust = c(0, 1), size = 2.6, fontface = "bold", family = FONT) +
    annotate("text", x = c(.85, 1.18), y = -.15, label = c("← Lower in acute", "Higher in acute →"), hjust = c(1, 0), size = 2.3, family = FONT, colour = "#555555") +
    scale_x_log10(limits = xl, breaks = c(.5, .75, 1, 1.5, 2, 3), labels = c("0.5", "0.75", "1.0", "1.5", "2.0", "3.0"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-.5, nrow(orx) + 1.2), expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = "Odds ratio, acute versus non-acute (log scale)", y = NULL, tag = "B") + th +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  save_all(pA / pB + plot_layout(heights = c(1.15, 1)), "Figure03_acute_vs_nonacute", 180, 175)
})

# ---- Figure 4 (+ S1): year trends by group -------------------------------------------------
fig("4", {
  tr <- rd("16_paper1_extra", "B_year_trends_overall_and_by_group.csv")[group %in% names(GC)]
  LBL <- c(comp_primary = "Primary composite outcome", dzh = "In-hospital death", d30 = "30-day death", d1y = "One-year death*",
           stroke_def = "Stroke with permanent deficit", stroke_nodef = "Stroke or TIA without deficit", renal = "Renal failure",
           gi = "GI complication", vent24 = "Ventilation >24 h", reintub = "Reintubation", icu_readm = "ICU readmission", major_compl = "Major complication")
  sec <- list(Primary = "comp_primary", Mortality = c("dzh", "d30", "d1y"),
              Morbidity = c("stroke_def", "stroke_nodef", "renal", "gi", "vent24", "reintub", "icu_readm", "major_compl"))
  build <- function(start, name, title) {
    x <- tr[from == start]; rows <- list(); hdr <- list(); y <- 0
    for (s in names(sec)) {
      if (s != "Primary") { y <- y - 1.1; hdr[[length(hdr) + 1]] <- data.table(y = y, lab = s); y <- y - .3 }
      for (o in intersect(sec[[s]], x$outcome)) { y <- y - 1.15
        rows[[length(rows) + 1]] <- data.table(outcome = o, lab = LBL[[o]], ylab = y - .27, group = c("Acute", "Non-acute"), y = c(y, y - .55)); y <- y - .55 } }
    R <- merge(rbindlist(rows), x, by = c("outcome", "group")); H <- rbindlist(hdr)
    R[, `:=`(sig = !is.na(p) & p < .05,
             ci = ifelse(is.na(or_per_year), "too few events", sprintf("%s (%s–%s)", f2(or_per_year), f2(lower), f2(upper))),
             en = sprintf("%s/%s", format(events, big.mark = ","), format(n, big.mark = ",")))]
    X <- c(lab = .30, grp = .565, en = .785, ci = 1.31, p = 1.415); top <- max(R$y) + 1.5; bot <- min(R$y) - .8
    p <- ggplot(R) + annotate("rect", xmin = .80, xmax = 1.22, ymin = bot + .2, ymax = top - .9, fill = "#F4F7F5") +
      geom_vline(xintercept = 1, linewidth = .35, colour = ink) +
      geom_segment(data = R[!is.na(or_per_year)], aes(x = pmax(lower, .80), xend = pmin(upper, 1.22), y = y, yend = y, colour = group), linewidth = .5) +
      geom_point(data = R[!is.na(or_per_year)], aes(or_per_year, y, colour = group, fill = ifelse(sig, group, "none")), shape = 22, size = 2.1, stroke = .6) +
      scale_colour_manual(values = GC) + scale_fill_manual(values = c(GC, none = "white")) +
      geom_text(data = unique(R[, .(lab, ylab)]), aes(X["lab"], ylab, label = lab), hjust = 0, size = 2.45, family = FONT) +
      geom_text(data = H, aes(X["lab"], y, label = lab), hjust = 0, size = 2.55, fontface = "bold", family = FONT) +
      geom_text(aes(X["grp"], y, label = group, colour = group), hjust = 0, size = 2.25, family = FONT, fontface = "bold") +
      geom_text(aes(X["en"], y, label = en), hjust = 1, size = 2.25, family = FONT) +
      geom_text(aes(X["ci"], y, label = ci), size = 2.3, family = FONT) + geom_text(aes(X["p"], y, label = fp(p)), size = 2.3, family = FONT) +
      annotate("text", x = X, y = top - .2, label = c("Outcome", "Group", "Events/n", "OR per year (95% CI)", "P"), hjust = c(0, 0, 1, .5, .5), size = 2.45, fontface = "bold", family = FONT) +
      annotate("segment", x = .30, xend = 1.44, y = top - .75, yend = top - .75, linewidth = .3) +
      annotate("text", x = c(.985, 1.015), y = bot - .2, label = c("← Decreasing", "Increasing →"), hjust = c(1, 0), size = 2.2, colour = "#555555", family = FONT) +
      scale_x_continuous(limits = c(.30, 1.45), breaks = c(.8, .9, 1, 1.1, 1.2), labels = c("0.80", "0.90", "1.00", "1.10", "1.20"), expand = c(0, 0)) +
      scale_y_continuous(limits = c(bot - .6, top + .3), expand = c(0, 0)) + coord_cartesian(clip = "off") +
      labs(x = "Adjusted odds ratio per calendar year", y = NULL, title = title) + theme_classic(base_size = 8.5, base_family = FONT) +
      theme(text = element_text(colour = ink), legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.text.x = element_text(colour = ink, size = 7.5), axis.title.x = element_text(hjust = .63),
            plot.title = element_text(size = 8.5, face = "bold"), plot.margin = margin(4, 6, 4, 6))
    save_all(p, name, 180, 175) }
  build(min(YRS), "Figure04_trends_by_group", sprintf("Operations %d–%d", min(YRS), max(YRS)))
  build(CONFIG$trend_sens_start, "FigureS1_trends_by_group_from2015", sprintf("Operations %d–%d", CONFIG$trend_sens_start, max(YRS)))
})

# ---- Figure 5: survival relative to the general population --------------------------------
fig("5", {
  rs <- rd("07_population", "relative_survival_ederer2.csv")[group %in% names(GC)]
  rs <- rbind(rs[, .SD[1], by = group][, `:=`(end = 0, observed = 1, relative = 1, lower = 1, upper = 1)], rs, fill = TRUE)[end <= 10]
  rs[, group := factor(group, names(GC))]
  ylo <- floor(min(c(rs$observed, rs$lower), na.rm = TRUE) * 10) / 10
  endl <- rs[, .SD[.N], by = group]
  pA <- ggplot(rs, aes(end, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) +
    geom_line(aes(y = observed), linewidth = .75) + geom_line(aes(y = relative), linewidth = .6, linetype = "22") +
    geom_text(data = endl, aes(x = end + .15, y = relative, label = sprintf("%s  %.1f%%", group, 100 * relative)), hjust = 0, size = 2.5, fontface = "bold", family = FONT) +
    annotate("segment", x = c(.2, .2), xend = c(.8, .8), y = ylo + c(.08, .04), yend = ylo + c(.08, .04), linetype = c("solid", "22"), linewidth = c(.75, .6)) +
    annotate("text", x = .95, y = ylo + c(.08, .04), label = c("Observed (Kaplan–Meier)", "Relative (Ederer II), with 95% CI"), hjust = 0, size = 2.5, family = FONT) +
    scale_colour_manual(values = GC) + scale_fill_manual(values = GC) +
    scale_x_continuous(limits = c(0, max(rs$end) * 1.45), breaks = 0:10, expand = c(0, 0)) +
    scale_y_continuous(limits = c(ylo, 1.01), breaks = seq(ylo, 1, .1), labels = function(x) 100 * x, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = "Years since surgery", y = "Survival (%)", tag = "A") + th
  sm <- rd("07_population", "smr_and_excess_by_followup_interval.csv")[group %in% paste("Group:", names(GC)) & interval %in% c("0-30 days", "30 days-1 year", "1-5 years")]
  sm[, `:=`(group = factor(sub("Group: ", "", group), c("Non-acute", "Acute")), interval = factor(interval, c("0-30 days", "30 days-1 year", "1-5 years"), c("0–30 days", "30 days–\n1 year", "1–5 years")))]
  pB <- ggplot(sm, aes(interval, smr, fill = group)) + geom_hline(yintercept = 1, linewidth = .35, linetype = "22") +
    geom_col(position = position_dodge(.8), width = .75) +
    geom_errorbar(aes(ymin = pmax(lower, .5), ymax = upper), position = position_dodge(.8), width = .2, linewidth = .3, colour = ink) +
    geom_text(aes(y = upper, label = ifelse(smr >= 10, formatC(smr, format = "f", digits = 0), formatC(smr, format = "f", digits = 1))),
              position = position_dodge(.8), vjust = -.5, size = 2.4, family = FONT) +
    scale_fill_manual(values = c("Non-acute" = "#93C9AE", "Acute" = GC[["Acute"]])) +
    scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300), expand = expansion(mult = c(0, .12))) +
    labs(x = "Interval after surgery", y = "Standardised mortality ratio (log scale)", tag = "B") + th +
    theme(legend.position = c(.85, .85), legend.title = element_blank(), legend.key.size = unit(3, "mm"))
  save_all(pA | pB, "Figure05_population_survival", 180, 85)
})

# ---- Figure 6: landmark at day 30 ---------------------------------------------------------
fig("6", {
  LM <- CONFIG$landmark_days
  hr <- rd("16_paper1_extra", "J_hazard_ratios_acute.csv"); lr <- rd("16_paper1_extra", "J_logrank.csv")
  sv <- cohort[!is.na(time_d) & !is.na(status) & urgency2 %in% names(GC)]
  sv[, `:=`(t_days = pmax(time_d, .5), grp = factor(as.character(urgency2), names(GC)))]
  km_df <- function(fit) { s <- summary(fit, censored = TRUE)
    d <- data.table(t = s$time, surv = s$surv, lower = s$lower, upper = s$upper, n_risk = s$n.risk, grp = factor(sub("^grp=", "", as.character(s$strata)), names(GC)))
    rbind(data.table(t = 0, surv = 1, lower = 1, upper = 1, n_risk = NA_real_, grp = factor(names(GC), names(GC))), d)[is.na(n_risk) | n_risk >= KM_MIN_RISK] }
  risk_at <- function(fit, times) { s <- summary(fit, times = times, extend = TRUE)
    data.table(t = s$time, n = s$n.risk, grp = factor(sub("^grp=", "", as.character(s$strata)), names(GC))) }
  fitA <- survfit(Surv(pmin(t_days, LM), as.integer(status == 1 & t_days <= LM)) ~ grp, data = sv, conf.type = "log-log")
  fitB <- survfit(Surv((t_days - LM) / 365.25, status) ~ grp, data = sv[t_days > LM], conf.type = "log-log")
  kmA <- km_df(fitA); kmB <- km_df(fitB); XMAX <- min(10, floor(max(kmB$t))); kmB <- kmB[t <= XMAX]
  W1 <- 2.6; kmA[, x := t / LM * W1]; kmB[, x := W1 + t]
  yrs <- seq(0, XMAX, 2)
  rk <- rbind(risk_at(fitA, 0)[, x := 0], risk_at(fitB, yrs)[, x := W1 + t])
  ylo <- max(0, floor(10 * min(kmB$lower, kmA$lower, na.rm = TRUE) - .5) / 10); XR <- W1 + XMAX + 1.6
  step_rib <- function(km) km[order(grp, x), { n <- .N; if (n < 2) NULL else
    list(x = c(rbind(x[-n], x[-1])), lower = rep(lower[-n], each = 2), upper = rep(upper[-n], each = 2)) }, by = grp]
  eB <- kmB[, .SD[.N], by = grp][order(surv)]; if (diff(eB$surv) < .03) eB$surv <- mean(eB$surv) + c(-.015, .015)
  htxt <- function(r) sprintf("adjusted HR %s%s (95%% CI, %s–%s)", f2(r$hr), stars(r$p), f2(r$lower), f2(r$upper))
  txt <- paste0("Days 0–30: ", htxt(hr[interval == "0-30 days"]), "; log-rank P", fp(lr$p[1]),
                "\nAfter day 30: ", htxt(hr[grepl("landmark", interval)]), "; log-rank P", fp(lr$p[2]))
  txt <- gsub("P(\\.|[0-9])", "P=\\1", txt)
  p <- ggplot() + annotate("rect", xmin = 0, xmax = W1, ymin = ylo, ymax = 1, fill = "#F3F6F4") +
    geom_ribbon(data = step_rib(kmA), aes(x, ymin = lower, ymax = upper, fill = grp), alpha = .14) +
    geom_ribbon(data = step_rib(kmB), aes(x, ymin = lower, ymax = upper, fill = grp), alpha = .14) +
    geom_step(data = kmA, aes(x, surv, colour = grp), linewidth = .7) + geom_step(data = kmB, aes(x, surv, colour = grp), linewidth = .7) +
    geom_vline(xintercept = W1, linetype = "22", linewidth = .45, colour = ink) +
    annotate("text", x = W1 / 2, y = 1.012, label = sprintf("Days 0–%d", LM), vjust = 0, size = 2.7, fontface = "bold", family = FONT) +
    annotate("text", x = W1 + .15, y = 1.012, label = sprintf("Landmark at day %d: survival among %d-day survivors", LM, LM), hjust = 0, vjust = 0, size = 2.7, fontface = "bold", family = FONT) +
    geom_text(data = eB, aes(max(kmB$x) + .12, surv, label = grp, colour = grp), hjust = 0, size = 2.8, fontface = "bold", family = FONT) +
    annotate("text", x = W1 + .3, y = ylo + (1 - ylo) * .05, label = txt, hjust = 0, vjust = 0, size = 2.5, family = FONT, lineheight = 1.1) +
    scale_colour_manual(values = GC) + scale_fill_manual(values = GC) +
    scale_x_continuous(limits = c(0, XR), breaks = c(W1 * c(0, 10, 20, 30) / LM, W1 + yrs[-1]), labels = c("0", "10", "20", "30", yrs[-1]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ylo, 1.035), breaks = seq(ylo, 1, .1), labels = function(x) 100 * x, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = "Days after surgery (left) and years after day 30 (right)", y = "Survival (%)") + th
  rt <- ggplot(rk, aes(x, grp, label = formatC(n, big.mark = ",", format = "d"), colour = grp)) + geom_text(size = 2.5, family = FONT) +
    scale_colour_manual(values = GC) + scale_y_discrete(limits = rev(names(GC))) +
    scale_x_continuous(limits = c(0, XR), expand = c(0, 0)) + coord_cartesian(clip = "off") + labs(x = NULL, y = NULL, title = "No. at risk") +
    theme_void(base_size = 9, base_family = FONT) +
    theme(axis.text.y.left = element_text(colour = ink, size = 8, hjust = 1, margin = margin(r = 14)), plot.title = element_text(size = 8, face = "bold"),
          legend.position = "none", plot.margin = margin(0, 8, 4, 6))
  save_all(p / rt + plot_layout(heights = c(5, 1)), "Figure06_landmark30", 180, 115)
})

# ---- Figure 7: standardised risk by year and group ----------------------------------------
fig("7", {
  st <- rd("16_paper1_extra", "D_standardised_risk_by_year.csv")[group %in% names(GC)]
  st[, group := factor(group, names(GC))]
  pan <- function(o, ttl, tag) { d <- st[outcome == o]; top <- ceiling(max(100 * d$upper) / 5) * 5
    ggplot(d, aes(year, 100 * risk, colour = group, fill = group)) + early +
      geom_ribbon(aes(ymin = 100 * lower, ymax = 100 * upper), alpha = .18, colour = NA) + geom_line(linewidth = .8) + geom_point(size = 1.4) +
      scale_colour_manual(values = GC) + scale_fill_manual(values = GC) + YR +
      scale_y_continuous(limits = c(0, top), expand = c(0, 0)) + labs(x = NULL, y = "Standardised risk (%)", title = ttl, tag = tag) + th }
  save_all((pan("dzh", "In-hospital mortality", "A") | pan("d30", "30-day mortality", "B") | pan("comp_primary", "Primary composite outcome", "C")) /
             key_groups() + plot_layout(heights = c(1, .08)), "Figure07_standardised_by_group", 180, 80)
})

# ---- Figure 8: perfusion practice by year and group ----------------------------------------
fig("8", {
  pr <- rd("16_paper1_extra", "F_practice_by_year_and_group.csv")[group %in% names(GC)]
  pr[, group := factor(group, names(GC))]
  lev <- c("Normothermia (>34 °C)", "Mild (28.1–34 °C)", "Moderate (20.1–28 °C)", "Deep (≤20 °C)")
  hl <- melt(pr[, .(group, year, normo, mild, moderate, deep)], id.vars = c("group", "year"))
  hl[, variable := factor(variable, c("normo", "mild", "moderate", "deep"), lev)]
  qA <- ggplot(hl, aes(year, value, fill = variable)) + geom_col(width = .78) +
    geom_text(data = hl[variable == lev[4]], aes(y = value / 2, label = sprintf("%.0f", value)), size = 2, family = FONT, colour = "white") +
    facet_wrap(~ group, nrow = 1) + scale_fill_manual(values = setNames(c("#DCEDE3", "#A9D4BC", "#5FAE88", "#1B5E4B"), lev), breaks = rev(lev)) +
    scale_x_continuous(breaks = YRS, labels = function(x) ifelse(x %in% YRS[seq(1, length(YRS), 3)], x, ""), expand = c(.02, .02)) +
    scale_y_continuous(limits = c(0, 100.5), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = NULL, y = "Patients with recorded\nnadir temperature (%)", fill = NULL, tag = "A") + th +
    theme(legend.position = "right", legend.text = element_text(size = 7.5), legend.key.size = unit(3.5, "mm"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5, hjust = 0))
  ln <- function(y, ttl, lim, ylab, tag) ggplot(pr, aes(year, get(y), colour = group)) + early + geom_line(linewidth = .8) + geom_point(size = 1.4) +
    scale_colour_manual(values = GC) + YR + scale_y_continuous(limits = lim, expand = c(0, 0)) + labs(x = NULL, y = ylab, title = ttl, tag = tag) + th
  qB <- ln("temp_median", "Nadir temperature", c(floor(min(pr$temp_median, na.rm = TRUE)) - 2, ceiling(max(pr$temp_median, na.rm = TRUE)) + 2), "°C, median", "B")
  qC <- ln("arrest_pct", "Circulatory arrest", c(0, 100), "Patients (%)", "C")
  qD <- ln("ca_median", "Duration of circulatory arrest", c(0, ceiling(max(pr$ca_median, na.rm = TRUE) / 15 + 1) * 15), "Minutes, median", "D")
  save_all(qA / (qB | qC | qD) / key_groups() + plot_layout(heights = c(1.1, 1, .1)), "Figure08_practice_by_group", 180, 165)
})

# ---- Figure 9: practice per centre -----------------------------------------------------------
fig("9", {
  c9 <- rd("16_paper1_extra", "I_centre_practice_profile.csv")
  c9[, `:=`(centre = factor(centre_lab, rev(sort(centre_lab))), miss = temp_missing_pct > 50)]
  dotp <- function(v, ttl, lim, tag, ylab = TRUE, temp = FALSE) {
    p <- ggplot(c9, aes(get(v), centre)) + geom_segment(aes(x = lim[1], xend = get(v), yend = centre), colour = grid_col, linewidth = .4) +
      { if (temp) geom_point(aes(shape = miss), size = 2, colour = GC[["Acute"]], fill = "white", stroke = .6) else geom_point(size = 2, colour = GC[["Acute"]]) } +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) + scale_x_continuous(limits = lim, expand = c(0, 0)) +
      labs(x = NULL, y = if (ylab) "Centre" else NULL, title = ttl, tag = tag) + th +
      theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_col, linewidth = .3))
    if (!ylab) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p }
  nlab <- ggplot(c9, aes(1, centre, label = formatC(n, big.mark = ",", format = "d"))) + geom_text(size = 2.4, family = FONT, hjust = 1) +
    scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) + labs(title = "\nn") + theme_void(base_family = FONT) +
    theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 1), plot.margin = margin(6, 4, 4, 0))
  tl <- c(floor(min(c9$temp_median, na.rm = TRUE) / 5) * 5, ceiling(max(c9$temp_median, na.rm = TRUE) / 5) * 5)
  save_all((dotp("pct_acute", "Acute\npresentation, %", c(0, 100), "A") | nlab |
              dotp("temp_median", "Nadir temperature,\n°C (median)", tl, "B", FALSE, TRUE) |
              dotp("pct_deep", "Deep\nhypothermia, %", c(0, 100), "C", FALSE, TRUE) |
              dotp("pct_arrest", "Circulatory\narrest, %", c(0, 100), "D", FALSE) |
              dotp("ca_median", "Arrest duration,\nmin (median)", c(0, ceiling(max(c9$ca_median, na.rm = TRUE) / 20) * 20), "E", FALSE)) +
             plot_layout(widths = c(1, .25, 1, 1, 1, 1)), "Figure09_centre_practice", 180, 95)
})

# ---- Figure 10: complications and failure to rescue by year and group ------------------------
fig("10", {
  f <- rd("16_paper1_extra", "G_ftr_decomposition_by_year_and_group.csv")[group %in% names(GC)]
  f[, `:=`(group = factor(group, names(GC)), deaths = deaths_compl + deaths_no_compl)]
  spec <- list(compl_pct = c("n_compl", "n"), ftr_pct = c("deaths_compl", "n_compl"),
               mort_no_compl_pct = c("deaths_no_compl", "n_nc"), mortality_pct = c("deaths", "n"))
  f[, n_nc := n - n_compl]
  trend <- function(k) rbindlist(lapply(names(GC), function(gname) {
    d <- f[group == gname]; y <- d[[spec[[k]][1]]]; N <- d[[spec[[k]][2]]]
    m <- glm(cbind(y, N - y) ~ I(d$year - CONFIG$reference_year), family = binomial); g <- seq(min(YRS), max(YRS), length.out = 80)
    X <- cbind(1, g - CONFIG$reference_year); eta <- drop(X %*% coef(m)); se <- sqrt(rowSums((X %*% vcov(m)) * X)); cf <- summary(m)$coefficients[2, ]
    data.table(group = factor(gname, names(GC)), year = g, fit = 100 * plogis(eta), lo = 100 * plogis(eta - 1.96 * se), hi = 100 * plogis(eta + 1.96 * se),
               txt = sprintf("%s: OR %s (%s–%s); P%s", gname, f2(exp(cf[1])), f2(exp(cf[1] - 1.96 * cf[2])), f2(exp(cf[1] + 1.96 * cf[2])),
                             ifelse(cf[4] < .001, "<.001", paste0("=", fp(cf[4]))))) }))
  pan <- function(k, ttl, tag) { tr <- trend(k); lab <- unique(tr[, .(group, txt)])
    top <- ceiling(max(c(f[[k]], tr$hi), na.rm = TRUE) * 1.25 / 10) * 10
    ggplot() + early + geom_ribbon(data = tr, aes(year, ymin = lo, ymax = hi, fill = group), alpha = .15) +
      geom_line(data = f, aes(year, get(k), colour = group), linewidth = .3, alpha = .6) + geom_point(data = f, aes(year, get(k), colour = group), size = 1.3) +
      geom_line(data = tr, aes(year, fit, colour = group), linewidth = .9) +
      annotate("text", x = CONFIG$trend_sens_start - .3, y = top * c(.97, .89), label = lab$txt, hjust = 0, vjust = 1, size = 2.3, family = FONT, colour = GC[as.character(lab$group)]) +
      scale_colour_manual(values = GC) + scale_fill_manual(values = GC) + YR +
      scale_y_continuous(limits = c(0, top), expand = c(0, 0)) + labs(x = NULL, y = "Patients (%)", title = ttl, tag = tag) + th }
  k10 <- key_groups(list(annotate("segment", x = 8.6, xend = 9.1, y = .5, yend = .5, linewidth = .9, colour = ink),
                         annotate("text", x = 9.2, y = .5, label = "Crude trend", hjust = 0, size = 2.5, family = FONT)))
  save_all(((pan("compl_pct", "Major complication", "A") | pan("ftr_pct", "Failure to rescue", "B")) /
              (pan("mort_no_compl_pct", "Mortality without major complication", "C") | pan("mortality_pct", "In-hospital mortality", "D"))) / k10 +
             plot_layout(heights = c(1, 1, .09)), "Figure10_ftr_by_group", 180, 150)
})
log_msg("Done. Figures in ", OUT)

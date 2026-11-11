# Figure: 30-day landmark analysis, acute versus non-acute (NEJM style).
# Panel A: survival during the first 30 days. Panel B: survival among 30-day survivors,
# with time restarted at day 30. Adjusted hazard ratios are read from the output of
# 06_survival_analysis.R; curves and numbers at risk are computed from the cohort.
source("R/00_config.R"); source("R/functions.R")
suppressMessages({ library(data.table); library(survival); library(ggplot2); library(patchwork) })

SUB <- "06_survival"
OUT <- file.path(CONFIG$results_dir, "figures_paper1"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LM  <- CONFIG$landmark_days
FONT <- if (.Platform$OS.type == "windows") "sans" else "Liberation Sans"
c_ac <- "#1B5E4B"; c_na <- "#4E9C78"; ink <- "#1A1A1A"; grid_col <- "#E6E6E6"
COLS <- c("Non-acute" = c_na, "Acute" = c_ac)

cohort <- load_derived("cohort.rds")
sv <- cohort[!is.na(time_d) & !is.na(status) & urgency2 %in% names(COLS)]
sv[, `:=`(t_days = pmax(time_d, 0.5), grp = factor(as.character(urgency2), levels = names(COLS)))]

stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
fmt_p <- function(p) ifelse(p < .001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = 3)))
hr_txt <- function(r) sprintf("Adjusted HR %.2f%s\n(95%% CI, %.2f\u2013%.2f)", r$ratio, stars(r$p), r$ratio_lower, r$ratio_upper)

pw <- fread(file.path(res_path(SUB), "cox_piecewise_hazard_ratios.csv"))
lmc <- fread(file.path(res_path(SUB), "landmark30_cox.csv"))
hr_early <- pw[grepl("acute", term) & grepl("0-30", term)][1]
hr_lm <- lmc[grepl("acute vs non-acute", model) & grepl("^acute", term) & !grepl(":", term)][1]
if (!nrow(hr_early) || is.na(hr_early$ratio)) stop("Piecewise HR for 0-30 days not found; run 06_survival_analysis.R first")
if (!nrow(hr_lm) || is.na(hr_lm$ratio)) stop("Landmark HR not found; run 06_survival_analysis.R first")

km_df <- function(fit, scale) {
  s <- summary(fit, censored = TRUE)
  d <- data.table(t = s$time / scale, surv = s$surv, lower = s$lower, upper = s$upper, n_risk = s$n.risk,
                  grp = factor(sub("^grp=", "", as.character(s$strata)), levels = names(COLS)))
  d <- rbind(data.table(t = 0, surv = 1, lower = 1, upper = 1, n_risk = NA_real_, grp = factor(names(COLS), levels = names(COLS))), d)
  d[is.na(n_risk) | n_risk >= KM_MIN_RISK]
}
risk_at <- function(fit, times, scale) {
  s <- summary(fit, times = times * scale, extend = TRUE)
  data.table(t = s$time / scale, n = s$n.risk, grp = factor(sub("^grp=", "", as.character(s$strata)), levels = names(COLS)))
}

# ---- A: first 30 days (administrative censoring at day 30) ----
a <- copy(sv)[, `:=`(ta = pmin(t_days, LM), sa = as.integer(status == 1 & t_days <= LM))]
fitA <- survfit(Surv(ta, sa) ~ grp, data = a, conf.type = "log-log")
pA_lr <- 1 - pchisq(survdiff(Surv(ta, sa) ~ grp, data = a)$chisq, 1)
kmA <- km_df(fitA, 1); rkA <- risk_at(fitA, c(0, 10, 20, 30), 1)

# ---- B: landmark at day 30 ----
b <- sv[t_days > LM][, tb := (t_days - LM) / 365.25]
fitB <- survfit(Surv(tb, status) ~ grp, data = b, conf.type = "log-log")
pB_lr <- 1 - pchisq(survdiff(Surv(tb, status) ~ grp, data = b)$chisq, 1)
kmB <- km_df(fitB, 1); XMAX <- min(10, floor(max(kmB$t)))
kmB <- kmB[t <= XMAX]; rkB <- risk_at(fitB, seq(0, XMAX, 2), 1)
endB <- summary(fitB, times = 5, extend = TRUE)
est5 <- data.table(grp = sub("^grp=", "", as.character(endB$strata)), surv = endB$surv, lower = endB$lower,
                   upper = endB$upper, n_risk = endB$n.risk)

th <- theme_classic(base_size = 9, base_family = FONT) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink, size = 8),
        axis.line = element_line(linewidth = .4, colour = ink), axis.ticks = element_line(linewidth = .4, colour = ink),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = .3), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 12), plot.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(6, 10, 2, 6))
step_rib <- function(km) km[order(grp, t), {
  n <- .N
  if (n < 2) NULL else list(t = c(rbind(t[-n], t[-1])), lower = rep(lower[-n], each = 2), upper = rep(upper[-n], each = 2))
}, by = grp]
curve <- function(km, xlim, xbr, ylim, xlab, title, tag, txt) {
  ggplot(km, aes(t, surv, colour = grp, fill = grp)) +
    geom_ribbon(data = step_rib(km), aes(x = t, ymin = lower, ymax = upper, fill = grp), inherit.aes = FALSE, alpha = .12) +
    geom_step(linewidth = .7) +
    annotate("text", x = xlim[2] * .03, y = ylim[1] + diff(ylim) * .08, label = txt, hjust = 0, vjust = 0,
             size = 2.6, family = FONT, colour = ink, lineheight = 1.05) +
    scale_colour_manual(values = COLS) + scale_fill_manual(values = COLS) +
    scale_x_continuous(limits = xlim, breaks = xbr, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], .1), labels = function(x) 100 * x, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = xlab, y = "Survival (%)", title = title, tag = tag) + th
}
legend_lab <- function(p, km, x) {
  e <- km[, .SD[.N], by = grp][order(surv)]
  gap <- .028
  if (nrow(e) == 2 && diff(e$surv) < gap) { m <- mean(e$surv); e$surv <- c(m - gap / 2, m + gap / 2) }
  p + geom_text(data = e, aes(x = x, y = surv, label = grp), hjust = 0, size = 2.9, fontface = "bold", family = FONT)
}
risk_tab <- function(rk, xlim, xbr) {
  ggplot(rk, aes(t, grp, label = formatC(n, big.mark = ",", format = "d"), colour = grp)) +
    geom_text(size = 2.6, family = FONT) + scale_colour_manual(values = COLS) +
    scale_y_discrete(limits = rev(names(COLS))) +
    scale_x_continuous(limits = xlim, breaks = xbr, expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = "No. at risk") + theme_void(base_size = 9, base_family = FONT) +
    theme(axis.text.y = element_text(colour = ink, size = 8, hjust = 1),
          plot.title = element_text(size = 8, face = "bold", hjust = 0), legend.position = "none",
          plot.margin = margin(0, 10, 4, 6), axis.text.y.left = element_text(colour = ink, size = 8, hjust = 1, margin = margin(r = 16)))
}

xA <- c(0, LM); xbA <- seq(0, LM, 10)
pA <- curve(kmA, c(0, LM * 1.4), xbA, c(.7, 1), "Days since surgery", "Days 0–30", "A",
            paste0(hr_txt(hr_early), "\nLog-rank ", fmt_p(pA_lr)))
pA <- legend_lab(pA, kmA, LM * 1.02)
xB <- c(0, XMAX * 1.2); xbB <- seq(0, XMAX, 2)
pB <- curve(kmB, xB, xbB, c(.4, 1), "Years since day 30", "Among 30-day survivors (landmark at day 30)", "B",
            paste0(hr_txt(hr_lm), "\nLog-rank ", fmt_p(pB_lr)))
pB <- legend_lab(pB, kmB, max(kmB$t) + XMAX * .015)

fig <- (pA + pB + plot_layout(widths = c(1, 1.9))) /
       (risk_tab(rkA, c(0, LM * 1.4), xbA) + risk_tab(rkB, xB, xbB) + plot_layout(widths = c(1, 1.9))) +
       plot_layout(heights = c(4, 1))
W <- 180; H <- 105
ggsave(file.path(OUT, "Figure_landmark30.pdf"), fig, width = W, height = H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_landmark30.svg"), fig, width = W, height = H, units = "mm", device = grDevices::svg)
ggsave(file.path(OUT, "Figure_landmark30.png"), fig, width = W, height = H, units = "mm", dpi = 600)
ggsave(file.path(OUT, "Figure_landmark30.tiff"), fig, width = W, height = H, units = "mm", dpi = 600, compression = "lzw")

cat("\nLandmark population: ", nrow(b), " patients alive beyond day ", LM, ", ", sum(b$status), " deaths\n", sep = "")
cat("Excluded (died or censored on or before day ", LM, "): ", nrow(sv) - nrow(b), "\n", sep = "")
cat("Log-rank days 0-30: ", fmt_p(pA_lr), "; landmark: ", fmt_p(pB_lr), "\n", sep = "")
cat("Piecewise HR 0-30 d: ", gsub("\n", " ", hr_txt(hr_early)), "\nLandmark HR: ", gsub("\n", " ", hr_txt(hr_lm)), "\n", sep = "")
cat("Five-year survival after day 30:\n"); print(est5[, .(grp, surv = round(100 * surv, 1), lower = round(100 * lower, 1),
                                                           upper = round(100 * upper, 1), n_risk)])
message("Figuur opgeslagen in ", OUT)

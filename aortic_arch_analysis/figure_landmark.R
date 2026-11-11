# Figure: 30-day landmark analysis, acute versus non-acute (NEJM style), in one graph.
# Days 0-30 on the left; at day 30 (vertical line) the Kaplan-Meier curves restart at 100%
# among 30-day survivors. The time axis is stretched before day 30 so both parts are readable. Adjusted hazard ratios are read from the output of
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

# first 30 days (administrative censoring at day 30) ----
a <- copy(sv)[, `:=`(ta = pmin(t_days, LM), sa = as.integer(status == 1 & t_days <= LM))]
fitA <- survfit(Surv(ta, sa) ~ grp, data = a, conf.type = "log-log")
pA_lr <- 1 - pchisq(survdiff(Surv(ta, sa) ~ grp, data = a)$chisq, 1)
kmA <- km_df(fitA, 1); rkA <- risk_at(fitA, c(0, 10, 20, 30), 1)

# landmark at day 30
b <- sv[t_days > LM][, tb := (t_days - LM) / 365.25]
fitB <- survfit(Surv(tb, status) ~ grp, data = b, conf.type = "log-log")
pB_lr <- 1 - pchisq(survdiff(Surv(tb, status) ~ grp, data = b)$chisq, 1)
kmB <- km_df(fitB, 1); XMAX <- min(10, floor(max(kmB$t)))
kmB <- kmB[t <= XMAX]
endB <- summary(fitB, times = 5, extend = TRUE)
est5 <- data.table(grp = sub("^grp=", "", as.character(endB$strata)), surv = endB$surv, lower = endB$lower,
                   upper = endB$upper, n_risk = endB$n.risk)

W1 <- 2.6   # plotting width of days 0-30, in "years" units of the landmark part
xA <- function(d) d / LM * W1
xB <- function(y) W1 + y
yrs <- seq(2, 10, 2); yrs <- yrs[yrs - LM / 365.25 <= XMAX]
kmA[, x := xA(t)]; kmB[, x := xB(t)]
rkA <- rkA[t == 0][, x := xA(t)]
rkB0 <- risk_at(fitB, c(0, yrs - LM / 365.25), 1); rkB0[, x := xB(t)]
rk <- rbind(rkA[, .(x, n, grp)], rkB0[, .(x, n, grp)])
ylo <- max(0, floor(10 * min(kmB$lower, kmA$lower, na.rm = TRUE) - .5) / 10)
XR <- xB(XMAX) + 1.6

step_rib <- function(km) km[order(grp, x), {
  n <- .N
  if (n < 2) NULL else list(x = c(rbind(x[-n], x[-1])), lower = rep(lower[-n], each = 2), upper = rep(upper[-n], each = 2))
}, by = grp]
ends <- function(km) { e <- km[, .SD[.N], by = grp][order(surv)]
  if (diff(e$surv) < .03) { m <- mean(e$surv); e$surv <- c(m - .015, m + .015) }; e }
eB <- ends(kmB)
txt <- paste0("Days 0–30: ", gsub("\n", " ", hr_txt(hr_early)), "; log-rank ", fmt_p(pA_lr),
              "\nAfter day 30: ", gsub("\n", " ", hr_txt(hr_lm)), "; log-rank ", fmt_p(pB_lr))

th <- theme_classic(base_size = 9, base_family = FONT) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink, size = 8),
        axis.line = element_line(linewidth = .4, colour = ink), axis.ticks = element_line(linewidth = .4, colour = ink),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = .3), legend.position = "none",
        plot.margin = margin(14, 8, 2, 6))
brks <- c(xA(c(0, 10, 20, 30)), xB(yrs - LM / 365.25))
labs_x <- c("0", "10", "20", "30", as.character(yrs))
p <- ggplot() +
  annotate("rect", xmin = 0, xmax = W1, ymin = ylo, ymax = 1, fill = "#F3F6F4") +
  geom_ribbon(data = step_rib(kmA), aes(x, ymin = lower, ymax = upper, fill = grp), alpha = .14) +
  geom_ribbon(data = step_rib(kmB), aes(x, ymin = lower, ymax = upper, fill = grp), alpha = .14) +
  geom_step(data = kmA, aes(x, surv, colour = grp), linewidth = .7) +
  geom_step(data = kmB, aes(x, surv, colour = grp), linewidth = .7) +
  geom_vline(xintercept = W1, linetype = "22", linewidth = .45, colour = ink) +
  annotate("text", x = W1 / 2, y = 1.012, label = "Days 0–30", vjust = 0, size = 2.7, fontface = "bold", family = FONT) +
  annotate("text", x = W1 + .15, y = 1.012, label = "Landmark at day 30: survival among 30-day survivors",
           hjust = 0, vjust = 0, size = 2.7, fontface = "bold", family = FONT) +
  geom_text(data = eB, aes(max(kmB$x) + .12, surv, label = sprintf("%s", grp), colour = grp), hjust = 0, size = 2.8, fontface = "bold", family = FONT) +
  annotate("text", x = W1 + .3, y = ylo + (1 - ylo) * .05, label = txt, hjust = 0, vjust = 0, size = 2.5, family = FONT, lineheight = 1.1) +
  scale_colour_manual(values = COLS) + scale_fill_manual(values = COLS) +
  scale_x_continuous(limits = c(0, XR), breaks = brks, labels = labs_x, expand = c(0, 0)) +
  scale_y_continuous(limits = c(ylo, 1.035), breaks = seq(ylo, 1, .1), labels = function(x) 100 * x, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Time since surgery (days 0–30, then years)", y = "Survival (%)") + th

rt <- ggplot(rk, aes(x, grp, label = formatC(n, big.mark = ",", format = "d"), colour = grp)) +
  geom_text(size = 2.5, family = FONT) + scale_colour_manual(values = COLS) +
  scale_y_discrete(limits = rev(names(COLS))) +
  scale_x_continuous(limits = c(0, XR), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "No. at risk") + theme_void(base_size = 9, base_family = FONT) +
  theme(axis.text.y.left = element_text(colour = ink, size = 8, hjust = 1, margin = margin(r = 14)),
        plot.title = element_text(size = 8, face = "bold", hjust = 0), legend.position = "none",
        plot.margin = margin(0, 8, 4, 6))

fig <- p / rt + plot_layout(heights = c(5, 1))
W <- 180; H <- 115
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

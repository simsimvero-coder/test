# =============================================================================
# Perfusie: trend in minimale temperatuur, rekening houdend met centrum
#
# 1. Random-intercept-model (nlme::lme): min_temp ~ jaar + (1 | centrum)
# 2. Centrum-fixed-effects: lm(min_temp ~ jaar + centrum) met CR2 cluster-robuuste
#    standaardfouten (geclusterd op centrum) en Satterthwaite-vrijheidsgraden
#    (pakket clubSandwich). Bij weinig centra (hier 15) zijn deze
#    betrouwbaarheidsintervallen eerlijker dan die van het random-intercept-model.
# =============================================================================

temperatuur_modellen <- function(df, groep = GROEP, uitkomst = TREND_TEMPERATUUR,
                                 jaar_kolom = TREND_JAAR, centrum = CENTRUM_KOLOM) {
  if (!all(c(uitkomst, jaar_kolom, centrum) %in% names(df))) return(NULL)
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE)
  if (!heeft_cs)
    message("LET OP - pakket 'clubSandwich' niet geïnstalleerd; CR2-analyse overgeslagen. ",
            "Installeer met: install.packages(\"clubSandwich\", repos = \"https://cloud.r-project.org\", type = \"binary\")")

  subsets <- if (!is.null(groep) && groep %in% names(df))
    split(df, droplevels(factor(df[[groep]]))) else list()
  subsets[["Totaal"]] <- df

  rij <- function(g, model, d, est = NA, lo = NA, hi = NA, p = NA)
    data.frame(groep = g, model = model, n = nrow(d), n_centra = length(unique(d$.centrum)),
               verandering_per_jaar = round(est, 3), BI95_onder = round(lo, 3),
               BI95_boven = round(hi, 3), p = fmt_p(p))

  do.call(rbind, lapply(names(subsets), function(g) {
    d <- subsets[[g]]
    d <- data.frame(y = d[[uitkomst]], jaar = suppressWarnings(as.numeric(as.character(d[[jaar_kolom]]))),
                    .centrum = factor(d[[centrum]]))
    d <- droplevels(d[complete.cases(d), ])
    if (nrow(d) < 10 || nlevels(d$.centrum) < 2) return(NULL)

    # 1. Random intercept per centrum
    ri <- tryCatch({
      m <- nlme::lme(y ~ jaar, random = ~ 1 | .centrum, data = d, method = "REML")
      t <- summary(m)$tTable["jaar", ]
      q <- qt(0.975, t["DF"])
      rij(g, "random intercept (centrum)", d, t["Value"], t["Value"] - q * t["Std.Error"],
          t["Value"] + q * t["Std.Error"], t["p-value"])
    }, error = function(e) rij(g, "random intercept (centrum) - mislukt", d))

    # 2. Centrum-fixed-effects met CR2 (Satterthwaite)
    fe <- if (heeft_cs) tryCatch({
      m <- lm(y ~ jaar + .centrum, data = d)
      ci <- clubSandwich::conf_int(m, vcov = "CR2", cluster = d$.centrum, coefs = "jaar")
      ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = d$.centrum,
                                    test = "Satterthwaite", coefs = "jaar")
      pkol <- grep("^p_", names(ct), value = TRUE)[1]
      rij(g, "centrum fixed effects, CR2 (Satterthwaite)", d, ci$beta, ci$CI_L, ci$CI_U, ct[[pkol]])
    }, error = function(e) rij(g, "centrum fixed effects, CR2 - mislukt", d))
    else rij(g, "centrum fixed effects, CR2 - clubSandwich ontbreekt", d)

    rbind(ri, fe)
  }))
}

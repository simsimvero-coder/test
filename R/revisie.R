# =============================================================================
# Aanvullende analyses voor de revisie
#
# 1. cohort_definitie:      hoe het cohort in het bronbestand is afgebakend
#                           (incl_cohort, excl_reden, klasse, s1_via_bewijs, aortachir,
#                           aotachir_benadering = open / (T)EVAR / hybride, aorta_asc, aorta_desc).
# 2. pathologie_urgentie:   pathologie (hoogste prioriteit over asc/boog/desc) x urgentie
#                           x groep, en de verdeling per jaar.
# 3. niet_ingedeeld:        niet-ingedeelde patiënten per jaar en centrum, met reden.
# 4. model_vergelijking:    per uitkomst en groep alle modellen naast elkaar, met
#                           n en events: (a) random intercept (primair), (b) random
#                           intercept + random slope voor jaar, (c) centrum fixed effects
#                           + CR2 (Satterthwaite-vrijheidsgraden), plus toetsen op
#                           niet-lineariteit van jaar (natural spline 3 df en jaar als
#                           factor, likelihood-ratiotoets binnen model a).
# 5. ontbrekend_jaar_groep / ontbrekend_centrum_groep: % ontbrekend per variabele.
# 6. followup_overzicht:    wie zit in de overlevingsanalyse, mediane (potentiële)
#                           follow-up (reverse Kaplan-Meier), en per operatiejaar wie
#                           meetelt voor 30-dagen- en 1-jaarsmortaliteit.
# =============================================================================

code_getal <- function(x) suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", as.character(x))))

groep_tekst <- function(df, groep = GROEP) {
  g <- as.character(df[[groep]]); g[is.na(g)] <- "niet ingedeeld"; g
}

telling <- function(...) {
  t <- as.data.frame(table(..., useNA = "ifany"), responseName = "n")
  t <- t[t$n > 0, , drop = FALSE]
  for (k in setdiff(names(t), "n")) { t[[k]] <- as.character(t[[k]]); t[[k]][is.na(t[[k]])] <- "ontbreekt" }
  t
}

# ── 1. Cohortdefinitie ──────────────────────────────────────────────────────
cohort_definitie <- function(ruw, df_alle, groep = GROEP) {
  kol_ruw <- intersect(c("incl_cohort", "excl_reden", "klasse", "s1_via_bewijs"), names(ruw))
  bron <- do.call(rbind, lapply(kol_ruw, function(k) {
    t <- telling(waarde = ruw[[k]]); data.frame(bestand = "bronbestand (alle rijen)", variabele = k, t)
  }))
  kol_coh <- intersect(c("aortachir", "aortachir_clean", "aotachir_benadering", "chir_benadering", "aorta_asc",
                         "aortaboog_path", "aorta_desc", "ECC", "circ_arrest_clean", "klasse", "s1_via_bewijs"),
                       names(df_alle))
  coh <- do.call(rbind, lapply(kol_coh, function(k) {
    t <- telling(waarde = df_alle[[k]], groep = groep_tekst(df_alle, groep))
    data.frame(bestand = "inclusiecohort", variabele = k, t)
  }))
  rbind_vul(bron, coh)
}

# ── 2. Pathologie x urgentie ────────────────────────────────────────────────
pathologie_categorie <- function(df, cfg = ACUUT_ELECTIEF) {
  k <- intersect(cfg$path_kolommen, names(df))
  p <- matrix(sapply(k, function(v) code_getal(df[[v]])), nrow = nrow(df))
  heeft <- function(codes) rowSums(matrix(p %in% codes, nrow = nrow(df)), na.rm = TRUE) > 0
  out <- ifelse(rowSums(!is.na(p)) == 0, "niet ingevuld",
         ifelse(heeft(20), "acute dissectie",
         ifelse(heeft(30), "chronische dissectie",
         ifelse(heeft(10), "aneurysma",
         ifelse(heeft(40), "infectie",
         ifelse(heeft(c(90, 91)), "overig", "geen (alle 0)"))))))
  factor(out, levels = c("acute dissectie", "chronische dissectie", "aneurysma", "infectie",
                         "overig", "geen (alle 0)", "niet ingevuld"))
}

pathologie_urgentie <- function(df_alle, groep = GROEP) {
  pat <- pathologie_categorie(df_alle)
  urg <- as.character(df_alle[[ACUUT_ELECTIEF$urgentie_kolom]])
  g <- groep_tekst(df_alle, groep)
  kruis <- telling(pathologie = pat, urgentie = urg, groep = g)
  kruis <- kruis[order(kruis$groep, kruis$pathologie, kruis$urgentie), ]
  jaar <- telling(jaar = df_alle[[TREND_JAAR]], groep = g, pathologie = pat)
  tot <- ave(jaar$n, jaar$jaar, jaar$groep, FUN = sum)
  jaar$pct_binnen_groep_jaar <- round(100 * jaar$n / tot, 1)
  jaar <- jaar[order(jaar$groep, jaar$jaar, jaar$pathologie), ]
  list(kruistabel = kruis, per_jaar = jaar)
}

# ── 3. Niet-ingedeelde patiënten per jaar en centrum, met reden ─────────────
niet_ingedeeld <- function(df_alle, groep = GROEP) {
  niet <- is.na(df_alle[[groep]])
  reden <- df_alle$acuut_electief_basis
  reden[is.na(reden)] <- "niet ingedeeld: geen bruikbare urgentie of pathologie"
  per <- function(kolom) {
    x <- df_alle[[kolom]]
    tot <- table(x)
    t <- telling(x = x[niet], reden = reden[niet])
    names(t)[1] <- kolom
    t$geregistreerd <- as.integer(tot[t[[kolom]]])
    t$pct_van_geregistreerd <- round(100 * t$n / t$geregistreerd, 1)
    t[order(code_getal(t[[kolom]]), t$reden), ]
  }
  list(jaar = per(TREND_JAAR), centrum = per(CENTRUM_KOLOM))
}

# ── 4. Modelvergelijking en niet-lineariteit ────────────────────────────────
model_vergelijking <- function(df, groep = GROEP) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  heeft_cs <- requireNamespace("clubSandwich", quietly = TRUE)
  uitkomsten <- uitkomst_lijst(df, meldingen = FALSE)
  basis <- data.frame(jaar = suppressWarnings(as.numeric(as.character(df[[TREND_JAAR]]))),
                      .centrum = factor(df[[CENTRUM_KOLOM]]), groep = as.character(df[[groep]]))
  basis <- voeg_modelcovariaten_toe(basis, df)
  basis$jaar_c <- basis$jaar - round(mean(basis$jaar, na.rm = TRUE))
  ctrl <- lme4::glmerControl(optimizer = "bobyqa")
  glmer_ <- function(f, d) suppressMessages(suppressWarnings(lme4::glmer(f, data = d, family = binomial, control = ctrl)))
  ci <- function(b, se) sprintf("%.3f (%.3f-%.3f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
  lrt <- function(m0, m1) tryCatch(anova(m0, m1)[2, "Pr(>Chisq)"], error = function(e) NA_real_)
  res <- list()
  for (g in sort(unique(na.omit(basis$groep)))) for (u in uitkomsten) {
    d <- basis; d$y <- u[[3]]
    d <- d[d$groep == g & !is.na(d$y) & !is.na(d$jaar) & !is.na(d$.centrum), ]
    n_uitkomst <- nrow(d)
    d <- modelgevallen(d, "y")
    voortgang("modelvergelijking: %s - %s (n = %d)", g, u[[1]], nrow(d))
    rij <- data.frame(groep = g, uitkomst = u[[1]], type = u[[2]],
                      n_met_uitkomst = n_uitkomst, n_complete_case = nrow(d),
                      pct_complete_case = if (n_uitkomst) round(100 * nrow(d) / n_uitkomst, 1) else NA_real_,
                      events = sum(d$y), jaren = if (nrow(d)) sprintf("%d-%d", min(d$jaar), max(d$jaar)) else "",
                      A_random_intercept_OR = "", A_p = NA_real_, A_SD_centrum_intercept = NA_real_,
                      B_random_slope_OR = "", B_p = NA_real_, B_SD_centrum_slope = NA_real_, B_p_LRT_slope = NA_real_,
                      C_FE_CR2_OR = "", C_p = NA_real_, C_df_Satterthwaite = NA_real_,
                      p_nietlineair_spline3 = NA_real_, p_nietlineair_jaar_factor = NA_real_,
                      status = "", stringsAsFactors = FALSE)
    if (nrow(d) < 30 || sum(d$y) < 10 || sum(!d$y) < 10 || nlevels(droplevels(d$.centrum)) < 3) {
      rij$status <- "onvoldoende observaties of events"; res[[length(res) + 1]] <- rij; next
    }
    cov <- MODEL_TERMS
    st <- character()
    mA <- tryCatch(glmer_(stats::reformulate(c("jaar_c", cov, "(1 | .centrum)"), "y"), d), error = function(e) NULL)
    if (!is.null(mA)) {
      s <- coef(summary(mA))["jaar_c", ]
      rij$A_random_intercept_OR <- ci(s[1], s[2]); rij$A_p <- s[4]
      rij$A_SD_centrum_intercept <- round(attr(lme4::VarCorr(mA)$.centrum, "stddev")[1], 3)
      mS <- tryCatch(glmer_(stats::reformulate(c("splines::ns(jaar_c, 3)", cov, "(1 | .centrum)"), "y"), d), error = function(e) NULL)
      if (!is.null(mS)) rij$p_nietlineair_spline3 <- lrt(mA, mS)
      mF <- tryCatch(glmer_(stats::reformulate(c("factor(jaar)", cov, "(1 | .centrum)"), "y"), d), error = function(e) NULL)
      if (!is.null(mF)) rij$p_nietlineair_jaar_factor <- lrt(mA, mF)
    } else st <- c(st, "model A mislukt")
    mB <- tryCatch(glmer_(stats::reformulate(c("jaar_c", cov, "(1 + jaar_c | .centrum)"), "y"), d), error = function(e) NULL)
    if (!is.null(mB)) {
      s <- coef(summary(mB))["jaar_c", ]
      rij$B_random_slope_OR <- ci(s[1], s[2]); rij$B_p <- s[4]
      rij$B_SD_centrum_slope <- round(attr(lme4::VarCorr(mB)$.centrum, "stddev")[2], 3)
      if (!is.null(mA)) rij$B_p_LRT_slope <- mengverdeling_p(max(0, 2 * (as.numeric(logLik(mB)) - as.numeric(logLik(mA)))))
      if (lme4::isSingular(mB)) st <- c(st, "model B singulier (slopevariantie ~0)")
    } else st <- c(st, "model B mislukt")
    if (heeft_cs) tryCatch({
      m <- glm(stats::reformulate(c("jaar_c", cov, ".centrum"), "y"), family = binomial, data = d)
      V <- clubSandwich::vcovCR(m, cluster = d$.centrum, type = "CR2")
      ct <- suppressWarnings(clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite", coefs = "jaar_c"))
      cf <- suppressWarnings(clubSandwich::conf_int(m, vcov = V, coefs = "jaar_c"))
      rij$C_FE_CR2_OR <- sprintf("%.3f (%.3f-%.3f)", exp(cf$beta), exp(cf$CI_L), exp(cf$CI_U))
      rij$C_p <- ct[[grep("^p_", names(ct), value = TRUE)[1]]]
      rij$C_df_Satterthwaite <- round(ct[[grep("^df", names(ct), value = TRUE)[1]]], 1)
    }, error = function(e) st <<- c(st, paste("CR2 mislukt:", conditionMessage(e))))
    rij$status <- if (length(st)) paste(st, collapse = "; ") else "OK"
    res[[length(res) + 1]] <- rij
  }
  out <- do.call(rbind, res)
  for (k in grep("(^|_)p(_|$)", names(out), value = TRUE)) out[[k]] <- vapply(out[[k]], fmt_p, character(1))
  out
}

# ── 5. Ontbrekende waarden per jaar / centrum en groep ──────────────────────
ONTBREKEND_VARIABELEN <- c("leeftijd", "geslacht", "cardiochir_eerder", "urgentie", "aorta_asc_path",
                           "aortaboog_path", "aorta_desc_path", "EuroscoreII", "min_temp", "circ_arrest",
                           "circ_arrest_clean", "circ_arrest_duur", "ECC_duur", "AoX_duur",
                           "overleden_zh", "overl30d", "overl1j", "overltijd_dag", "mort_status",
                           "CVA_restletsel", "CVA_zonder_restletsel", "nierfalen", "gastroint_compl",
                           "resp_insuff", "beademing", "delier", "opnameduur")

ontbrekend_per <- function(df, kolom, groep = GROEP) {
  vars <- intersect(ONTBREKEND_VARIABELEN, names(df))
  g <- groep_tekst(df, groep)
  sleutel <- split(seq_len(nrow(df)), list(g, df[[kolom]]), drop = TRUE)
  out <- do.call(rbind, lapply(names(sleutel), function(s) {
    i <- sleutel[[s]]
    r <- data.frame(groep = g[i[1]], x = as.character(df[[kolom]][i[1]]), n = length(i))
    for (v in vars) r[[v]] <- round(100 * mean(is.na(df[[v]][i])), 1)
    r
  }))
  names(out)[2] <- kolom
  out[order(out$groep, code_getal(out[[kolom]])), ]
}

# ── 6. Follow-up ────────────────────────────────────────────────────────────
followup_overzicht <- function(df, groep = GROEP) {
  if (!all(c(OVERLEVING_TIJD, OVERLEVING_STATUS) %in% names(df))) return(NULL)
  g <- groep_tekst(df, groep)
  tijd <- df[[OVERLEVING_TIJD]]; status <- as.character(df[[OVERLEVING_STATUS]])
  jaar <- code_getal(df[[TREND_JAAR]])
  reden <- ifelse(is.na(tijd) & is.na(status), "follow-uptijd en vitale status ontbreken",
           ifelse(is.na(tijd), "follow-uptijd ontbreekt",
           ifelse(is.na(status), "vitale status ontbreekt",
           ifelse(tijd < 0, "negatieve follow-uptijd", "in overlevingsanalyse"))))
  wie <- telling(groep = g, reden = reden)
  # Mediane follow-up: reverse Kaplan-Meier (overlijden = gecensureerd)
  ok <- reden == "in overlevingsanalyse"
  dood <- status %in% "ja"
  mediaan <- do.call(rbind, lapply(split(which(ok), g[ok]), function(i) {
    rk <- survival::survfit(survival::Surv(tijd[i], !dood[i]) ~ 1)
    q <- quantile(rk, c(.25, .5, .75))$quantile
    data.frame(groep = g[i[1]], n = length(i), overleden = sum(dood[i]),
               mediane_follow_up_dagen_reverseKM = round(q[2]),
               IQR_dagen = sprintf("%.0f-%.0f", q[1], q[3]),
               mediane_geobserveerde_follow_up_levenden = round(median(tijd[i][!dood[i]])),
               max_follow_up_dagen = max(tijd[i]),
               n_levend_ge_1j = sum(!dood[i] & tijd[i] >= 365),
               n_levend_ge_5j = sum(!dood[i] & tijd[i] >= 1825))
  }))
  # Per operatiejaar: wie telt mee voor 30-dagen- en 1-jaarsmortaliteit
  regel <- function(var, min_fu) {
    y <- code_positief(df[[var]])
    yfu <- mortaliteit_fu(y, tijd, df[[TREND_JAAR]], min_fu)
    weg <- attr(yfu, "jaren_uitgesloten")
    kort <- !(y %in% TRUE) & (is.na(tijd) | tijd < min_fu)
    do.call(rbind, lapply(split(seq_len(nrow(df)), list(g, jaar), drop = TRUE), function(i) data.frame(
      uitkomst = var, groep = g[i[1]], jaar = jaar[i[1]], n = length(i),
      overleden = sum(y[i] %in% TRUE),
      levend_volledig_gevolgd = sum(!(y[i] %in% TRUE) & !kort[i]),
      onvolledige_follow_up = sum(kort[i]),
      jaar_uitgesloten = jaar[i[1]] %in% weg,
      in_analyse = sum(!is.na(yfu[i])))))
  }
  per_jaar <- rbind(regel("overl30d", 30), regel("overl1j", 365))
  per_jaar <- per_jaar[order(per_jaar$uitkomst, per_jaar$groep, per_jaar$jaar), ]
  list(wie = wie, mediaan = mediaan, per_jaar = per_jaar)
}

# ── 7. Dekking van de FTR-definitie ─────────────────────────────────────────
# Welk deel van de ziekenhuissterfte valt onder de FTR-definitie (overleden na
# een van de MAJOR_COMPLICATIES), en hoe verandert FTR bij een bredere definitie
# (alle geregistreerde postoperatieve complicaties in FTR_BREDE_COMPLICATIES)?
# De registratie bevat geen tijdstip van complicatie of overlijden: dat de
# complicatie aan het overlijden voorafging, kan niet worden vastgesteld.
FTR_BREDE_COMPLICATIES <- c(MAJOR_COMPLICATIES, "CVA_zonder_restletsel", "resp_insuff", "beademing",
                            "heropname_IC", "vasc_compl_opname", "ECMO_opname", "postop_MI", "delier")

ftr_dekking <- function(df, groep = GROEP) {
  dood <- code_positief(df[[FTR_UITKOMST]])
  major <- samengesteld(df, MAJOR_COMPLICATIES)
  breed_vars <- intersect(FTR_BREDE_COMPLICATIES, names(df))
  breed <- samengesteld(df, breed_vars)
  g <- groep_tekst(df, groep)
  do.call(rbind, lapply(split(seq_len(nrow(df)), g), function(i) {
    d <- dood[i]; m <- major[i]; b <- breed[i]
    nd <- sum(d %in% TRUE)
    data.frame(groep = g[i[1]], n = length(i), ziekenhuissterfte = nd,
               sterfte_met_major_complicatie = sum(d %in% TRUE & m %in% TRUE),
               pct_sterfte_gedekt_door_FTR = round(100 * sum(d %in% TRUE & m %in% TRUE) / nd, 1),
               sterfte_met_brede_complicatie = sum(d %in% TRUE & b %in% TRUE),
               pct_sterfte_gedekt_breed = round(100 * sum(d %in% TRUE & b %in% TRUE) / nd, 1),
               sterfte_zonder_geregistreerde_complicatie = sum(d %in% TRUE & b %in% FALSE),
               sterfte_complicatiestatus_onbekend = sum(d %in% TRUE & is.na(b)),
               FTR_major_pct = sprintf("%.1f (%d/%d)", 100 * mean(d[m %in% TRUE], na.rm = TRUE),
                                       sum(d[m %in% TRUE], na.rm = TRUE), sum(m %in% TRUE & !is.na(d))),
               FTR_breed_pct = sprintf("%.1f (%d/%d)", 100 * mean(d[b %in% TRUE], na.rm = TRUE),
                                       sum(d[b %in% TRUE], na.rm = TRUE), sum(b %in% TRUE & !is.na(d))),
               brede_definitie = paste(breed_vars, collapse = ", "))
  }))
}

# ── 8. Volume: overdispersie en standaardisatie voor leeftijd en geslacht ────
# Overdispersie: Pearson-dispersie van het Poissonmodel, en de RR per jaar uit een
# quasi-Poisson- en een negatief-binomiaal model.
volume_overdispersie <- function(per_jaar, vanaf = GROEPEN_VANAF_JAAR) {
  reeksen <- intersect(c("totaal", unname(ACUUT_ELECTIEF$labels)), names(per_jaar))
  do.call(rbind, lapply(reeksen, function(k) {
    d <- per_jaar[!is.na(per_jaar$bevolking), ]
    if (k != "totaal") d <- d[d$jaar >= vanaf, ]
    d$n <- d[[k]]
    mp <- glm(n ~ jaar + offset(log(bevolking)), family = poisson, data = d)
    mq <- glm(n ~ jaar + offset(log(bevolking)), family = quasipoisson, data = d)
    disp <- sum(residuals(mp, type = "pearson")^2) / mp$df.residual
    sq <- coef(summary(mq))["jaar", ]
    nb <- tryCatch({
      m <- suppressWarnings(MASS::glm.nb(n ~ jaar + offset(log(bevolking)), data = d))
      s <- coef(summary(m))["jaar", ]
      c(sprintf("%.3f (%.3f-%.3f)", exp(s[1]), exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2])), fmt_p(s[4]), round(m$theta, 1))
    }, error = function(e) c("", "", NA))
    data.frame(reeks = k, periode = sprintf("%d-%d", min(d$jaar), max(d$jaar)), jaren = nrow(d),
               dispersie_pearson = round(disp, 2),
               p_overdispersie = fmt_p(pchisq(disp * mp$df.residual, mp$df.residual, lower.tail = FALSE)),
               RR_quasipoisson = sprintf("%.3f (%.3f-%.3f)", exp(sq[1]), exp(sq[1] - qt(.975, mq$df.residual) * sq[2]),
                                         exp(sq[1] + qt(.975, mq$df.residual) * sq[2])),
               p_quasipoisson = fmt_p(sq[4]), RR_negbin = nb[1], p_negbin = nb[2], theta_negbin = nb[3])
  }))
}

# Bevolking per jaar, geslacht en leeftijd (CBS StatLine 7461BEV, "Bevolking op
# 1 januari"; kolommen geslacht, leeftijd, jaar/perioden en aantal, ';' of ',').
lees_bevolking_leeftijd <- function(bestand = BEVOLKING_LEEFTIJD_BESTAND) {
  pad <- zoek_bestand(bestand)
  if (is.null(pad)) return(NULL)
  eerste <- readLines(pad, n = 1, warn = FALSE)
  b <- utils::read.csv(pad, sep = if (grepl(";", eerste)) ";" else ",", check.names = FALSE,
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  kies <- function(p) names(b)[grep(p, names(b), ignore.case = TRUE)[1]]
  getal <- function(x) suppressWarnings(as.numeric(gsub(",", ".", gsub("[ .](?=\\d{3}\\b)", "", sub("^\\D*(\\d[0-9 .,]*).*$", "\\1", x), perl = TRUE))))
  out <- data.frame(
    geslacht = ifelse(grepl("^m", tolower(b[[kies("geslacht|sex")]])), "man",
               ifelse(grepl("^v|^f", tolower(b[[kies("geslacht|sex")]])), "vrouw", NA)),
    leeftijd = suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", b[[kies("leeftijd|age")]]))),
    jaar = suppressWarnings(as.numeric(substr(gsub("\\D", "", b[[kies("jaar|period|year")]]), 1, 4))),
    aantal = getal(b[[names(b)[ncol(b)]]]))
  out[!is.na(out$geslacht) & !is.na(out$leeftijd) & !is.na(out$jaar) & !is.na(out$aantal), ]
}

LEEFTIJDSKLASSEN <- c(-Inf, 49, 59, 69, 79, Inf)
LEEFTIJDSLABELS <- c("<50", "50-59", "60-69", "70-79", ">=80")

# Operaties per miljoen, direct gestandaardiseerd naar de bevolking van het
# laatste jaar (leeftijd x geslacht), en een Poissonmodel met leeftijd-geslacht-
# strata en log(bevolking per stratum) als offset: RR per jaar gecorrigeerd voor
# de vergrijzing.
volume_gestandaardiseerd <- function(df_alle, groep = GROEP, vanaf = GROEPEN_VANAF_JAAR) {
  bev <- lees_bevolking_leeftijd()
  if (is.null(bev) || !all(c("leeftijd", "geslacht") %in% names(df_alle))) return(NULL)
  bev$klasse <- cut(bev$leeftijd, LEEFTIJDSKLASSEN, LEEFTIJDSLABELS)
  bev <- aggregate(aantal ~ jaar + geslacht + klasse, bev, sum)
  pt <- data.frame(jaar = code_getal(df_alle[[TREND_JAAR]]),
                   geslacht = ifelse(grepl("^m", tolower(df_alle$geslacht)), "man",
                              ifelse(grepl("^v", tolower(df_alle$geslacht)), "vrouw", NA)),
                   klasse = cut(as.numeric(df_alle$leeftijd), LEEFTIJDSKLASSEN, LEEFTIJDSLABELS),
                   groep = groep_tekst(df_alle, groep))
  ontbreekt <- sum(!complete.cases(pt))
  pt <- pt[complete.cases(pt), ]
  standaard_jaar <- max(intersect(pt$jaar, bev$jaar))
  std <- bev[bev$jaar == standaard_jaar, c("geslacht", "klasse", "aantal")]
  names(std)[3] <- "standaard"
  reeksen <- c("totaal", intersect(unname(ACUUT_ELECTIEF$labels), pt$groep))
  per_jaar <- list(); modellen <- list()
  for (k in reeksen) {
    x <- if (k == "totaal") pt else pt[pt$groep == k, ]
    tel <- aggregate(n ~ jaar + geslacht + klasse, transform(x, n = 1), sum)
    d <- merge(bev[bev$jaar %in% unique(pt$jaar), ], tel, all.x = TRUE)
    d$n[is.na(d$n)] <- 0
    if (k != "totaal") d <- d[d$jaar >= vanaf, ]
    d <- merge(d, std)
    d$rate <- d$n / d$aantal
    jr <- do.call(rbind, lapply(split(d, d$jaar), function(z) data.frame(
      reeks = k, jaar = z$jaar[1], n = sum(z$n), bevolking = sum(z$aantal),
      ruw_per_miljoen = round(1e6 * sum(z$n) / sum(z$aantal), 1),
      gestandaardiseerd_per_miljoen = round(1e6 * sum(z$rate * z$standaard) / sum(z$standaard), 1))))
    per_jaar[[k]] <- jr
    m0 <- glm(n ~ jaar + offset(log(aantal)), family = quasipoisson, data = d)
    m1 <- glm(n ~ jaar + geslacht * klasse + offset(log(aantal)), family = quasipoisson, data = d)
    r <- function(m) { s <- coef(summary(m))["jaar", ]; q <- qt(.975, m$df.residual)
      c(sprintf("%.3f (%.3f-%.3f)", exp(s[1]), exp(s[1] - q * s[2]), exp(s[1] + q * s[2])), fmt_p(s[4])) }
    modellen[[k]] <- data.frame(reeks = k, periode = sprintf("%d-%d", min(d$jaar), max(d$jaar)),
                                RR_ruw_quasipoisson = r(m0)[1], p_ruw = r(m0)[2],
                                RR_leeftijd_geslacht = r(m1)[1], p_leeftijd_geslacht = r(m1)[2],
                                standaardbevolking = standaard_jaar,
                                patienten_zonder_leeftijd_of_geslacht = ontbreekt)
  }
  list(per_jaar = do.call(rbind, per_jaar), modellen = do.call(rbind, modellen))
}

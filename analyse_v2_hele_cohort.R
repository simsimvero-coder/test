###############################################################################
#  VERSIE 2 - analyses die beter werken met kleine aantallen
#  ---------------------------------------------------------------------------
#  A. Data-beschikbaarheid per parameter en periode
#  B. Pre->post verandering in het HELE cohort (periode als covariaat),
#     hoofdmodel: pre-waarde + tijd tot TTE + EuroSCORE II (+ periode)
#     -> primaire uitkomsten: GLS, RV S', LA-strain (Holm over deze 3)
#     -> bootstrap-BI's en mixed model als gevoeligheidsanalyse
#  C. Aortaklemtijd -> verandering kamerfunctie (hele cohort)
#  D. Mediatie met lavaan + FIML (gebruikt ook patiënten met ontbrekende data)
#  E. Klinische uitkomsten na vs voor juni 2022:
#     Firth-logistische regressie, DAOH-30, aantal complicaties (ordinaal),
#     plus (benaderde) Bayesiaanse kans op voordeel
#
#  Vereist: analyse_hartfunctie_postop.R in DEZELFDE map (voor de data-opschoning).
#  Gebruik:  source(file.choose(), encoding = "UTF-8", chdir = TRUE)
#  Output:   map output_v2/
###############################################################################

## ======================= INSTELLINGEN v2 ===================================
n_boot      <- 2000          # bootstrap voor BI's van de verandering
n_boot_med  <- 500           # bootstrap voor mediatie (lavaan); hoger = trager
prior_sd_zwak       <- 1.5   # zwak-informatieve prior op log-OR (95%: OR 0,05-19)
prior_sd_sceptisch  <- 0.35  # sceptische prior op log-OR (95%: OR 0,5-2)
## ===========================================================================

# ---- 0. Data-opschoning uit script 1 hergebruiken ----------------------------
v1 <- "analyse_hartfunctie_postop.R"
if (!file.exists(v1)) stop("Zet ", v1, " in dezelfde map als dit script.")
code <- readLines(v1, encoding = "UTF-8", warn = FALSE)
eind <- grep("^# 3\\. TABEL 1", code)[1] - 2
eval(parse(text = code[seq_len(eind)], encoding = "UTF-8"), envir = globalenv())

out_dir <- "output_v2"; dir.create(out_dir, showWarnings = FALSE)
out_pad <- normalizePath(out_dir)
schrijf <- function(x, naam) if (!is.null(x) && NROW(x) > 0)
  write.csv(x, file.path(out_dir, paste0(naam, ".csv")), row.names = FALSE)
for (p in c("logistf", "lavaan")) if (!heeft(p))
  try(install.packages(p, type = if (.Platform$pkgType == "source") "source" else "binary"))
set.seed(2022)

prim_par <- c("gls", "rv_s", "la_strain")          # primair
sec_par  <- c("lvef", "tapse", "e_eprime")         # secundair
alle_v2  <- c(prim_par, sec_par)
gunstig_hoger <- c(gls = TRUE, rv_s = TRUE, la_strain = TRUE, lvef = TRUE, tapse = TRUE, e_eprime = FALSE)
ref_t <- log2(ref_dag_tte + 1)

dc <- d %>% filter(!is.na(periode)) %>%
  mutate(per = as.integer(periode == "na_jun2022"),
         es2 = euroscore2,
         log_icu = log1p(icu_los),
         n_compl = rowSums(across(all_of(compl_vars)), na.rm = TRUE),
         # DAOH-30: dagen levend en buiten het ziekenhuis in de eerste 30 dagen na OK
         dood30 = ja(x30d_mortaliteit) == 1 |
                  (!is.na(overleden_dagen_postop) & overleden_dagen_postop <= 30),
         daoh30 = ifelse(dood30, 0, pmax(0, 30 - pmin(in_hospital_los, 30))))

rond <- function(x) { x <- as.data.frame(x); num <- sapply(x, is.numeric)
  x[num] <- lapply(x[num], function(v) ifelse(abs(v) < 0.001 & v != 0, signif(v, 2), round(v, 3))); x }
toon <- function(x, titel) { if (is.null(x) || !NROW(x)) return(invisible())
  cat("\n==================", titel, "==================\n"); print(rond(x), row.names = FALSE) }
bayes <- function(est, se, prior_sd) {          # normaal-normaal benadering (log-schaal)
  v <- 1 / (1 / prior_sd^2 + 1 / se^2); m <- v * est / se^2
  c(post = m, laag = m - 1.96 * sqrt(v), hoog = m + 1.96 * sqrt(v), p_kleiner_0 = pnorm(0, m, sqrt(v)))
}

###############################################################################
# A. DATA-BESCHIKBAARHEID
###############################################################################
tab_A <- bind_rows(lapply(alle_v2, function(p) dc %>% group_by(periode) %>%
  summarise(parameter = p, n_ok = n(), pre = sum(!is.na(.data[[paste0("pre_", p)]])),
            post = sum(!is.na(.data[[paste0("post_", p)]])),
            gepaard = sum(!is.na(.data[[paste0("delta_", p)]])), .groups = "drop"))) %>%
  select(parameter, periode, everything())
schrijf(tab_A, "A_beschikbaarheid")
toon(tab_A, "A. Beschikbare echo's per periode")
cat("\nTAPSE-waarden (mm) ter controle:\n"); print(summary(dc$pre_tapse)); print(summary(dc$post_tapse))

###############################################################################
# B. VERANDERING IN HET HELE COHORT
###############################################################################
cov_sets <- list(
  hoofd       = "es2",
  plus_lft_gs = c("es2", "leeftijd_operatie", "geslacht"),
  plus_echo   = c("es2", intersect(c("post_ic", "post_intubatie"), varieert)))

adj_change <- function(m, dd, pre, periode_waarde = NULL) {
  nd <- dd; nd[[pre]] <- mean(dd[[pre]]); nd$tdub <- ref_t
  if (!is.null(periode_waarde)) nd$periode <- factor(periode_waarde, levels(dd$periode))
  mean(predict(m, newdata = nd))
}
res_B <- list()
for (p in alle_v2) for (cs in names(cov_sets)) {
  dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  rhs <- c(pre, "tdub", cov_sets[[cs]], "periode")
  dd <- dc[complete.cases(dc[, c(dl, rhs)]), ]
  if (nrow(dd) < 15) next
  if (n_distinct(dd$periode) < 2) rhs <- setdiff(rhs, "periode")
  f <- reformulate(rhs, dl); m <- lm(f, dd)
  at <- setNames(list(mean(dd[[pre]]), ref_t), c(pre, "tdub"))
  em <- summary(emmeans(m, ~ 1, at = at, weights = "proportional"))
  boot_ci <- c(NA, NA)
  if (cs == "hoofd") {
    bs <- replicate(n_boot, { i <- sample(nrow(dd), replace = TRUE); di <- dd[i, ]
      tryCatch(adj_change(lm(f, di), di, pre), error = function(e) NA) })
    boot_ci <- quantile(bs, c(.025, .975), na.rm = TRUE)
  }
  pr <- if ("periode" %in% rhs) tidy(m, conf.int = TRUE) %>% filter(term == "periodena_jun2022") else NULL
  gh <- gunstig_hoger[[p]]
  res_B[[length(res_B) + 1]] <- tibble(
    parameter = p, primair = p %in% prim_par, model = cs, n = nrow(dd),
    n_na2022 = sum(dd$periode == "na_jun2022"),
    adj_delta = em$emmean, ci_laag = em$lower.CL, ci_hoog = em$upper.CL,
    boot_ci_laag = boot_ci[1], boot_ci_hoog = boot_ci[2],
    p = 2 * pt(-abs(em$emmean / em$SE), em$df),
    tijd_per_verdubbeling = coef(m)[["tdub"]],
    periode_verschil = if (!is.null(pr)) pr$estimate else NA,
    periode_ci_laag = if (!is.null(pr)) pr$conf.low else NA,
    periode_ci_hoog = if (!is.null(pr)) pr$conf.high else NA,
    periode_p = if (!is.null(pr)) pr$p.value else NA,
    # kans (vlakke prior) dat de verandering na jun 2022 gunstiger is dan ervoor
    P_na2022_gunstiger = if (!is.null(pr)) { pg <- pnorm(0, pr$estimate, pr$std.error)
      if (gh) 1 - pg else pg } else NA)
}
tab_B <- bind_rows(res_B)
i <- tab_B$model == "hoofd" & tab_B$primair
tab_B$p_holm_primair <- NA; tab_B$p_holm_primair[i] <- p.adjust(tab_B$p[i], "holm")
schrijf(tab_B, "B_verandering_hele_cohort")

## B2. Mixed model: ook patiënten met alleen een pre- of alleen een post-echo
mm_v2 <- bind_rows(lapply(alle_v2, function(p) {
  kol <- intersect(paste0(c("pre_", "post_"), p), names(dc))
  lg <- dc %>% select(interv_nr, es2, periode, post_dagen, all_of(kol)) %>%
    pivot_longer(all_of(kol), names_to = "moment", values_to = "w") %>%
    mutate(moment = factor(sub("_.*", "", moment), c("pre", "post")),
           tpost = ifelse(moment == "post", log2(post_dagen + 1) - ref_t, 0)) %>%
    filter(!is.na(w), !is.na(es2), moment == "pre" | !is.na(post_dagen))
  if (n_distinct(lg$interv_nr) < 20 || sum(lg$moment == "post") < 10) return(NULL)
  m <- suppressMessages(lmer(w ~ moment + tpost + es2 + periode + (1 | interv_nr), lg))
  ci <- confint(m, parm = "momentpost", method = "Wald")
  tibble(parameter = p, n_pat = n_distinct(lg$interv_nr), n_echo = nrow(lg),
         post_vs_pre = fixef(m)[["momentpost"]], ci_laag = ci[1], ci_hoog = ci[2])
}))
schrijf(mm_v2, "B2_mixed_model_hele_cohort")

###############################################################################
# C. AORTAKLEMTIJD (hele cohort, alleen on-pump)
###############################################################################
op <- dc %>% filter(!is.na(aox_duur), aox_duur > 0)
res_C <- list()
for (p in alle_v2) {
  dl <- paste0("delta_", p); pre <- paste0("pre_", p)
  for (extra in list(character(0), "circ_arrest_duur", "ecc_zonder_klem")) {
    extra <- intersect(extra, varieert)
    rhs <- c("aox30", extra, pre, "tdub", "es2", "periode")
    dd <- op[complete.cases(op[, c(dl, rhs)]), ]; if (nrow(dd) < 20) next
    if (n_distinct(dd$periode) < 2) rhs <- setdiff(rhs, "periode")
    m <- lm(reformulate(rhs, dl), dd)
    tt <- tidy(m, conf.int = TRUE) %>% filter(term == "aox30")
    p_nl <- NA
    if (!length(extra)) {
      m2 <- lm(reformulate(c("ns(aox_duur, 3)", setdiff(rhs, "aox30")), dl), dd)
      p_nl <- tryCatch(anova(m, m2)$`Pr(>F)`[2], error = function(e) NA)
    }
    res_C[[length(res_C) + 1]] <- tibble(parameter = p, extra_correctie = paste(extra, collapse = "+"),
      n = nrow(dd), effect_per_30min = tt$estimate, ci_laag = tt$conf.low, ci_hoog = tt$conf.high,
      p = tt$p.value, p_niet_lineair = p_nl)
  }
}
tab_C <- bind_rows(res_C); schrijf(tab_C, "C_klemtijd_hele_cohort")

###############################################################################
# D. MEDIATIE MET LAVAAN + FIML
###############################################################################
tab_D <- NULL
if (heeft("lavaan")) {
  dm <- op %>% transmute(aox = aox30, es2, per, tdub, log_icu, daoh30, tnt = log_post_hstnt,
                         across(any_of(c(paste0("pre_", prim_par), paste0("post_", prim_par),
                                         paste0("delta_", prim_par)))))
  run_med <- function(model, lbl) {
    fit <- tryCatch(suppressWarnings(lavaan::sem(model, data = dm, missing = "fiml", fixed.x = FALSE,
                    se = "bootstrap", bootstrap = n_boot_med)),
                    error = function(e) { message(lbl, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) return(NULL)
    pe <- lavaan::parameterEstimates(fit, boot.ci.type = "perc")
    pe <- pe[pe$label %in% c("a", "b", "c", "ind", "tot"), ]
    data.frame(analyse = lbl, pad = pe$label, schatting = pe$est, ci_laag = pe$ci.lower,
               ci_hoog = pe$ci.upper, p = pe$pvalue, n_gebruikt = lavaan::lavInspect(fit, "nobs"))
  }
  res_D <- list()
  for (p in prim_par[sapply(prim_par, function(p) all(paste0(c("pre_", "post_", "delta_"), p) %in% names(dm)))]) {
    for (y in c("log_icu", "daoh30")) {
      mod <- sprintf("post_%s ~ a*aox + pre_%s + tdub + es2 + per
                      %s ~ b*post_%s + c*aox + es2 + per
                      ind := a*b
                      tot := c + a*b", p, p, y, p)
      res_D[[length(res_D) + 1]] <- run_med(mod, sprintf("AoX -> post_%s -> %s", p, y))
    }
    mod <- sprintf("tnt ~ a*aox + es2 + per
                    delta_%s ~ b*tnt + c*aox + pre_%s + tdub + es2 + per
                    ind := a*b
                    tot := c + a*b", p, p)
    res_D[[length(res_D) + 1]] <- run_med(mod, sprintf("AoX -> hsTnT -> delta_%s", p))
  }
  tab_D <- bind_rows(res_D)
  schrijf(tab_D, "D_mediatie_lavaan_FIML")
}

###############################################################################
# E. KLINISCHE UITKOMSTEN NA vs VOOR JUNI 2022
###############################################################################
# Firth-logistische regressie (betrouwbaar bij weinig events). Gebruikt logistf als
# dat geïnstalleerd is (profiel-likelihood BI); anders een eigen implementatie (Wald-BI).
firth_eigen <- function(f, data) {
  mf <- model.frame(f, data); y <- model.response(mf); X <- model.matrix(f, mf)
  b <- rep(0, ncol(X))
  for (it in 1:200) {
    pr <- plogis(drop(X %*% b)); w <- pr * (1 - pr); XW <- X * sqrt(w)
    Iinv <- solve(crossprod(XW)); h <- rowSums((XW %*% Iinv) * XW)
    stap <- drop(Iinv %*% crossprod(X, y - pr + h * (0.5 - pr)))
    if (max(abs(stap)) > 5) stap <- stap * 5 / max(abs(stap))
    b <- b + stap; if (max(abs(stap)) < 1e-8) break
  }
  pr <- plogis(drop(X %*% b)); Iinv <- solve(crossprod(X * sqrt(pr * (1 - pr))))
  se <- sqrt(diag(Iinv)); names(b) <- names(se) <- colnames(X)
  list(b = b, se = se, ci_l = b - 1.96 * se, ci_h = b + 1.96 * se, p = 2 * pnorm(-abs(b / se)), methode = "Firth (Wald-BI)")
}
firth <- function(f, data) {
  if (heeft("logistf")) {
    r <- tryCatch({ fit <- logistf::logistf(f, data = data)
      b <- coef(fit); se <- sqrt(diag(fit$var)); names(se) <- names(b)
      ci_l <- fit$ci.lower; ci_h <- fit$ci.upper; names(ci_l) <- names(ci_h) <- names(b)
      list(b = b, se = se, ci_l = ci_l, ci_h = ci_h, p = setNames(fit$prob, names(b)), methode = "Firth (logistf, profiel-BI)")
    }, error = function(e) NULL)
    if (!is.null(r)) return(r)
  }
  firth_eigen(f, data)
}
term <- "periodena_jun2022"
uitk <- intersect(c("x30d_mortaliteit", "compl_major", "compl_any", "ritmeprobleem", "rethoracotomie",
                    "post_op_cardiale_dysfunctie", "nierfalen", "aki_zonder_dialyse", "delier",
                    "heropname_ic", "ecmo_opname"), names(dc))
res_E <- list()
for (y in uitk) for (adj in unique(list(c("es2"), c("es2", intersect("ingreep", varieert))))) {
  dd <- dc[complete.cases(dc[, c(y, "periode", adj)]), ]
  if (sum(dd[[y]]) < 3) next
  r <- firth(reformulate(c("periode", adj), y), dd)
  bz <- bayes(r$b[[term]], r$se[[term]], prior_sd_zwak); bs <- bayes(r$b[[term]], r$se[[term]], prior_sd_sceptisch)
  res_E[[length(res_E) + 1]] <- tibble(uitkomst = y, correctie = paste(adj, collapse = "+"), methode = r$methode,
    n = nrow(dd), events = sum(dd[[y]]),
    pct_voor = 100 * mean(dd[[y]][dd$periode == "voor_jun2022"]), pct_na = 100 * mean(dd[[y]][dd$periode == "na_jun2022"]),
    OR = exp(r$b[[term]]), ci_laag = exp(r$ci_l[[term]]), ci_hoog = exp(r$ci_h[[term]]), p = r$p[[term]],
    P_OR_kleiner1_zwakke_prior = bz[["p_kleiner_0"]], P_OR_kleiner1_sceptische_prior = bs[["p_kleiner_0"]],
    OR_post_sceptisch = exp(bs[["post"]]))
}
tab_E <- bind_rows(res_E); schrijf(tab_E, "E1_uitkomsten_Firth_Bayes")

## E2. DAOH-30 en aantal complicaties (meer power dan ja/nee-uitkomsten)
res_E2 <- list()
for (y in c("daoh30", "n_compl")) {
  dd <- dc[complete.cases(dc[, c(y, "periode", "es2")]), ]
  if (nrow(dd) < 20) next
  # gemiddeld verschil (gecorrigeerd voor ES2) met bootstrap-BI
  m <- lm(reformulate(c("periode", "es2"), y), dd)
  bs <- replicate(n_boot, { di <- dd[sample(nrow(dd), replace = TRUE), ]
    tryCatch(coef(lm(reformulate(c("periode", "es2"), y), di))[[term]], error = function(e) NA) })
  # ordinaal (proportional odds): OR > 1 = hogere waarde
  dd$y_ord <- if (y == "daoh30") cut(dd$daoh30, c(-Inf, 0, 15, 20, 25, 30), ordered_result = TRUE)
              else factor(pmin(dd$n_compl, 4), ordered = TRUE)
  po <- tryCatch(MASS::polr(y_ord ~ periode + es2, data = dd, Hess = TRUE), error = function(e) NULL)
  po_b <- po_se <- NA
  if (!is.null(po)) { po_b <- coef(po)[[term]]; po_se <- sqrt(diag(vcov(po)))[[term]] }
  res_E2[[length(res_E2) + 1]] <- tibble(uitkomst = y, n = nrow(dd),
    mediaan_voor = median(dd[[y]][dd$periode == "voor_jun2022"]), mediaan_na = median(dd[[y]][dd$periode == "na_jun2022"]),
    gem_verschil_adj = coef(m)[[term]], boot_ci_laag = quantile(bs, .025, na.rm = TRUE), boot_ci_hoog = quantile(bs, .975, na.rm = TRUE),
    p_wilcoxon = suppressWarnings(wilcox.test(dd[[y]] ~ dd$periode)$p.value),
    OR_ordinaal = exp(po_b), OR_ci_laag = exp(po_b - 1.96 * po_se), OR_ci_hoog = exp(po_b + 1.96 * po_se),
    P_na2022_gunstiger = if (is.na(po_b)) NA else { pg <- pnorm(0, po_b, po_se); if (y == "daoh30") 1 - pg else pg })
}
tab_E2 <- bind_rows(res_E2); schrijf(tab_E2, "E2_DAOH30_en_aantal_complicaties")

###############################################################################
# RESULTATEN TONEN
###############################################################################
toon(tab_B %>% filter(model == "hoofd") %>%
       select(parameter, primair, n, n_na2022, adj_delta, ci_laag, ci_hoog, boot_ci_laag, boot_ci_hoog,
              p, p_holm_primair, tijd_per_verdubbeling),
     "B. Verandering pre->post, HELE cohort (ES2 + tijd + periode), TTE op ref-dag")
toon(tab_B %>% filter(model == "hoofd") %>%
       select(parameter, periode_verschil, periode_ci_laag, periode_ci_hoog, periode_p, P_na2022_gunstiger),
     "B. Verschil in verandering na vs voor jun 2022")
toon(tab_B %>% filter(model != "hoofd") %>% select(parameter, model, n, adj_delta, ci_laag, ci_hoog, p),
     "B. Gevoeligheidsanalyses")
toon(mm_v2, "B2. Mixed model (alle beschikbare echo's)")
toon(tab_C, "C. Aortaklemtijd (per 30 min) -> verandering")
toon(tab_D, "D. Mediatie (lavaan, FIML, bootstrap)")
toon(tab_E %>% select(-methode), "E1. Uitkomsten na vs voor jun 2022 (Firth + Bayes)")
if (nrow(tab_E)) cat("Methode E1:", unique(tab_E$methode), "\n")
toon(tab_E2, "E2. DAOH-30 en aantal complicaties")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nAlle tabellen staan in:", out_pad, "\n")

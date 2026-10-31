# =============================================================================
# functions.R — shared helper functions for the arch surgery pipeline
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lme4)
  library(survival)
  library(splines)
})

# ---- logging & output --------------------------------------------------------

.LOG_FILE <- NULL

start_log <- function(script_name) {
  if (substr(script_name, 1, 2) %in% c("04", "05", "06", "08", "09"))
    script_name <- paste0(script_name, "_", CONFIG$missing_method)
  d <- out_dir("logs")
  .LOG_FILE <<- file.path(d, paste0(script_name, ".log"))
  cat(sprintf("[%s] start %s\n", format(Sys.time()), script_name), file = .LOG_FILE)
  set.seed(CONFIG$seed)
  invisible(.LOG_FILE)
}

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  message(msg)
  if (!is.null(.LOG_FILE)) cat(msg, "\n", file = .LOG_FILE, append = TRUE)
}

# Outputs of the adjusted analyses go to results/<missing_method>/<sub>;
# descriptive and population outputs (01, 02, 03, 07) are shared.
METHOD_SUBS <- c("04_outcomes", "05_organ_protection", "06_survival", "08_sensitivity", "publication")
MISSING_METHODS <- c("imputation", "complete_case")
res_path <- function(sub, method = CONFIG$missing_method) {
  if (sub %in% METHOD_SUBS) file.path(CONFIG$results_dir, method, sub) else file.path(CONFIG$results_dir, sub)
}
out_dir <- function(sub) {
  d <- res_path(sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

write_out <- function(x, sub, name) {
  f <- file.path(out_dir(sub), paste0(name, ".csv"))
  fwrite(as.data.table(x), f)
  invisible(f)
}

save_fig <- function(p, sub, name, width = CONFIG$fig_width, height = CONFIG$fig_height) {
  d <- out_dir(sub)
  ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height,
         dpi = CONFIG$fig_dpi)
  ggsave(file.path(d, paste0(name, ".pdf")), p, width = width, height = height)
  invisible(NULL)
}

derived_path <- function(name) {
  dir.create(CONFIG$derived_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(CONFIG$derived_dir, name)
}

# Derived objects of the adjusted analyses carry the method in their file name
derived_method <- function(name, method = CONFIG$missing_method) {
  derived_path(sub("\\.rds$", paste0("_", method, ".rds"), name))
}

# Data sets analysed: the m imputed data sets, or the observed cohort once
# (rows with missing model variables are then dropped model by model).
analysis_datasets <- function() {
  if (!CONFIG$missing_method %in% MISSING_METHODS) stop("Unknown missing_method: ", CONFIG$missing_method)
  if (CONFIG$missing_method == "imputation") return(load_derived("imputed.rds"))
  list(rederive(load_derived("cohort.rds")))
}

# Keep rows with all variables of a model formula observed (no-op after imputation,
# apart from missing outcomes)
complete_for <- function(datasets, formula) {
  vars <- setdiff(all.vars(formula), ".")
  lapply(datasets, function(d) {
    v <- intersect(vars, names(d))
    d[complete.cases(d[, .SD, .SDcols = v])]
  })
}

method_label <- function(method = CONFIG$missing_method) {
  c(imputation = sprintf("multiple imputation (m = %d)", CONFIG$n_imputations),
    complete_case = "complete-case analysis")[[method]]
}

load_derived <- function(name) {
  f <- derived_path(name)
  if (!file.exists(f)) stop("Derived file '", f, "' not found. Run the earlier scripts first.")
  readRDS(f)
}

theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11),
          strip.text = element_text(face = "bold"))
}

PAL <- c("#2a6ebb", "#d1495b", "#3d9970", "#edae49", "#6f4e9c", "#7f7f7f")

# ---- reading -----------------------------------------------------------------

read_raw <- function(path, sheet = CONFIG$raw_data_sheet) {
  if (!file.exists(path)) {
    searched <- attr(path, "searched")
    stop("Data file not found. Looked for:\n  ",
         paste(if (is.null(searched)) path else searched, collapse = "\n  "),
         "\nCheck the file name (RAW_FILE_NAME in R/00_config.R) or set the full path ",
         "with the environment variable ARCH_RAW_DATA.")
  }
  ext <- tolower(tools::file_ext(path))
  dt <- switch(ext,
    csv = , txt = , tsv = {
      first <- readLines(path, n = 1, warn = FALSE)
      semi <- grepl(";", first)
      fread(path, sep = if (semi) ";" else "auto", dec = if (semi) "," else ".",
            na.strings = c("", "NA", "NULL", "."), encoding = "UTF-8")
    },
    xlsx = , xls = {
      sh <- readxl::excel_sheets(path)
      if (is.null(sheet) || identical(sheet, "auto")) {
        # pick the sheet whose column names match the most expected variables
        ov <- vapply(sh, function(s) {
          nm <- tryCatch(names(readxl::read_excel(path, sheet = s, n_max = 0)), error = function(e) character(0))
          length(intersect(nm, EXPECTED_VARS))
        }, numeric(1))
        best <- sh[ov == max(ov)]
        sheet <- if ("Alle_records" %in% best) "Alle_records" else best[1]
        message("Excel sheets: ", paste(sprintf("%s (%d expected columns)", sh, ov), collapse = ", "),
                " -> reading '", sheet, "'")
      } else {
        message("Excel sheets: ", paste(sh, collapse = ", "), " -> reading '",
                if (is.numeric(sheet)) sh[sheet] else sheet, "'")
      }
      as.data.table(readxl::read_excel(path, sheet = sheet, guess_max = 1e6, na = c("", "NA")))
    },
    sav = as.data.table(haven::zap_labels(haven::read_sav(path))),
    dta = as.data.table(haven::zap_labels(haven::read_dta(path))),
    rds = as.data.table(readRDS(path)),
    stop("Unsupported file type: ", ext))
  dt
}

# Convert a column to numeric codes without losing information.
# Text yes/no/TRUE/FALSE are mapped to 1/0; anything else non-numeric -> NA and
# reported by the caller.
to_code <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  y <- trimws(as.character(x))
  y[tolower(y) %in% c("true", "ja", "yes", "y")] <- "1"
  y[tolower(y) %in% c("false", "nee", "no", "n")] <- "0"
  y <- gsub(",", ".", y, fixed = TRUE)
  suppressWarnings(as.numeric(y))
}

# ---- descriptive helpers ---------------------------------------------------------

fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
}
fmt_npct <- function(x) {
  n <- sum(x == 1, na.rm = TRUE); d <- sum(!is.na(x))
  sprintf("%d (%.1f%%)", n, 100 * n / max(d, 1))
}
fmt_med <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf("%.1f [%.1f-%.1f]", q[1], q[2], q[3])
}
fmt_mean <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

# Standardised mean difference (continuous or binary)
smd <- function(x, g) {
  g <- as.logical(g); ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NA_real_)
  m1 <- mean(x[g]); m0 <- mean(x[!g])
  if (all(x %in% c(0, 1))) {
    s <- sqrt((m1 * (1 - m1) + m0 * (1 - m0)) / 2)
  } else {
    s <- sqrt((var(x[g]) + var(x[!g])) / 2)
  }
  if (s == 0) return(0)
  (m1 - m0) / s
}

# SMD for a multi-category factor (Yang & Dalton 2012)
smd_factor <- function(f, g) {
  g <- as.logical(g); ok <- !is.na(f) & !is.na(g)
  f <- droplevels(as.factor(f[ok])); g <- g[ok]
  if (nlevels(f) < 2 || length(unique(g)) < 2) return(NA_real_)
  p1 <- prop.table(table(f[g]))[-1]; p0 <- prop.table(table(f[!g]))[-1]
  k <- length(p1); d <- p1 - p0
  S <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    S[i, j] <- if (i == j) (p1[i] * (1 - p1[i]) + p0[i] * (1 - p0[i])) / 2
               else -(p1[i] * p1[j] + p0[i] * p0[j]) / 2
  }
  out <- tryCatch(sqrt(as.numeric(t(d) %*% solve(S) %*% d)), error = function(e) NA_real_)
  out
}

# ---- splines with fixed knots ------------------------------------------------------
# Knots are computed once on the analysis cohort and reused in every model and
# every imputed data set, so spline bases are identical across fits.

spline_spec <- function(x, df) {
  x <- x[!is.na(x)]
  probs <- seq(0, 1, length.out = df + 1)[-c(1, df + 1)]
  list(knots = unname(quantile(x, probs)), boundary = range(x))
}

ns_term <- function(var, spec) {
  sprintf("ns(%s, knots = c(%s), Boundary.knots = c(%s))", var,
          paste(signif(spec$knots, 8), collapse = ", "),
          paste(signif(spec$boundary, 8), collapse = ", "))
}

# Map covariate names to model terms
covar_terms <- function(covars, specs) {
  vapply(covars, function(v) {
    if (v == "age")  return(ns_term("age",  specs$age))
    if (v == "egfr") return(ns_term("egfr", specs$egfr))
    v
  }, character(1), USE.NAMES = FALSE)
}

build_formula <- function(outcome, covars, specs, year_term = "year_c",
                          random = "(1 | centrum)") {
  rhs <- c(year_term, covar_terms(covars, specs))
  if (!is.null(random) && nzchar(random)) rhs <- c(rhs, random)
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")), env = globalenv())
}

year_spline_term <- function(specs) ns_term("year_c", specs$year)

# ---- model fitting ----------------------------------------------------------------

# Fit a (G)LMM, capturing warnings and convergence status.
# GLMM strategy (for speed without loss of accuracy): a fast nAGQ = 0 fit supplies
# starting values for the Laplace (nAGQ = 1) fit with the nloptwrap optimiser.
# lme4's gradient check is conservative at large n; following ?lme4::convergence,
# a warning is treated as a false positive when the scaled gradient
# (solve(Hessian, gradient)) is < 0.002. Otherwise the model is refitted with bobyqa.
scaled_gradient <- function(fit) {
  dd <- fit@optinfo$derivs
  if (is.null(dd) || is.null(dd$Hessian)) return(NA_real_)
  tryCatch(max(abs(solve(dd$Hessian, dd$gradient))), error = function(e) NA_real_)
}

fit_mixed <- function(formula, data, family = NULL, refit = TRUE) {
  warns <- character(0)
  catcher <- function(expr) withCallingHandlers(expr,
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) if (grepl("singular", conditionMessage(m))) invokeRestart("muffleMessage"))
  fitter <- function(opt, tight = FALSE) {
    oc <- if (opt == "nloptwrap") {
      if (tight) list(maxeval = 1e5, xtol_abs = 1e-10, ftol_abs = 1e-10) else list(maxeval = 1e5)
    } else list(maxfun = 2e5)
    catcher(
      if (is.null(family)) {
        lmer(formula, data = data, REML = FALSE,
             control = lmerControl(optimizer = opt, optCtrl = oc))
      } else {
        f0 <- glmer(formula, data = data, family = family, nAGQ = 0,
                    control = glmerControl(optimizer = "nloptwrap"))
        glmer(formula, data = data, family = family, nAGQ = 1,
              start = list(theta = getME(f0, "theta"), fixef = fixef(f0)),
              control = glmerControl(optimizer = opt, optCtrl = oc))
      })
  }
  first <- if (is.null(family)) "bobyqa" else "nloptwrap"
  fit <- tryCatch(fitter(first), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(fit = NULL, ok = FALSE, singular = NA, warnings = conditionMessage(fit),
                optimizer = first, scaled_gradient = NA_real_))
  }
  sg <- scaled_gradient(fit)
  conv_warn <- any(grepl("converge|gradient|Hessian", warns, ignore.case = TRUE)) &&
               !(isTRUE(sg < 0.002))
  opt_used <- first
  if (conv_warn && refit && !is.null(family)) {
    # restart nloptwrap (tight tolerances) from the current estimates, up to 5 times
    cur <- fit
    for (i in 1:5) {
      warns_i <- character(0)
      cur <- tryCatch(withCallingHandlers(
        update(cur, start = list(theta = getME(cur, "theta"), fixef = fixef(cur)),
               control = glmerControl(optimizer = "nloptwrap",
                                      optCtrl = list(maxeval = 1e5, xtol_abs = 1e-10, ftol_abs = 1e-10))),
        warning = function(w) { warns_i <<- c(warns_i, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) NULL)
      if (is.null(cur)) break
      sg_i <- scaled_gradient(cur)
      if (!any(grepl("converge|gradient|Hessian", warns_i, ignore.case = TRUE)) || isTRUE(sg_i < 0.002)) {
        fit <- cur; sg <- sg_i; conv_warn <- FALSE; opt_used <- sprintf("nloptwrap (%d restart%s)", i, if (i > 1) "s" else "")
        break
      }
    }
  }
  if (conv_warn && refit && !is.null(family) && CONFIG$glmm_fallback == "none") {
    # keep the best restart; it remains counted as a convergence warning
    if (!is.null(cur)) { fit <- cur; sg <- scaled_gradient(cur); opt_used <- "nloptwrap (restarts exhausted)" }
  } else if (conv_warn && refit) {
    warns2 <- warns; warns <- character(0)
    fit2 <- tryCatch(fitter("bobyqa"), error = function(e) NULL)
    sg2 <- if (!is.null(fit2)) scaled_gradient(fit2) else NA
    if (!is.null(fit2) && (!any(grepl("converge|gradient|Hessian", warns, ignore.case = TRUE)) || isTRUE(sg2 < 0.002))) {
      fit <- fit2; opt_used <- "bobyqa"; conv_warn <- FALSE; sg <- sg2
    } else {
      warns <- warns2
    }
  }
  list(fit = fit, ok = !conv_warn, singular = isSingular(fit), warnings = unique(warns),
       optimizer = opt_used, scaled_gradient = sg)
}

# Parallel lapply over imputed data sets (forked processes; serial on Windows)
par_lapply <- function(X, FUN) {
  cores <- CONFIG$cores
  if (.Platform$OS.type == "windows" || cores <= 1 || length(X) == 1) return(lapply(X, FUN))
  res <- parallel::mclapply(X, FUN, mc.cores = min(cores, length(X)), mc.preschedule = FALSE)
  err <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(err)) stop("Parallel fit failed: ", as.character(res[[which(err)[1]]]))
  res
}

# Rubin's rules for a list of coefficient vectors and covariance matrices
rubin_pool <- function(est_list, vcov_list) {
  m <- length(est_list)
  Q <- do.call(rbind, est_list)
  qbar <- colMeans(Q)
  U <- Reduce(`+`, lapply(vcov_list, function(v) diag(as.matrix(v)))) / m
  B <- if (m > 1) apply(Q, 2, var) else rep(0, ncol(Q))
  Tv <- U + (1 + 1 / m) * B
  r <- ifelse(U > 0, (1 + 1 / m) * B / U, 0)
  df <- ifelse(B > 0, (m - 1) * (1 + 1 / r)^2, Inf)
  se <- sqrt(Tv)
  tq <- qt(0.975, df)
  data.table(term = names(qbar), estimate = qbar, se = se, df = df,
             lower = qbar - tq * se, upper = qbar + tq * se,
             p = 2 * pt(abs(qbar / se), df, lower.tail = FALSE),
             fmi = ifelse(Tv > 0, (1 + 1 / m) * B / Tv, 0))
}

# Fit the same mixed model on a list of data sets and pool the fixed effects
fit_pooled <- function(formula, datasets, family = binomial()) {
  datasets <- complete_for(datasets, formula)
  fits <- par_lapply(datasets, function(d) fit_mixed(formula, d, family))
  ok <- vapply(fits, function(f) !is.null(f$fit), logical(1))
  if (!any(ok)) stop("All fits failed for: ", deparse(formula)[1])
  est <- lapply(fits[ok], function(f) fixef(f$fit))
  n_vcov_fallback <- 0L
  vc  <- lapply(fits[ok], function(f) withCallingHandlers(as.matrix(vcov(f$fit)), warning = function(w) {
    if (grepl("falling back", conditionMessage(w))) n_vcov_fallback <<- n_vcov_fallback + 1L
    invokeRestart("muffleWarning") }))
  pooled <- rubin_pool(est, vc)
  re_sd <- vapply(fits[ok], function(f) {
    s <- as.data.frame(VarCorr(f$fit)); s$sdcor[is.na(s$var2)][1]
  }, numeric(1))
  diag <- data.table(
    n_fits = length(fits), n_failed = sum(!ok),
    n_convergence_warning = sum(vapply(fits[ok], function(f) !f$ok, logical(1))),
    n_singular = sum(vapply(fits[ok], function(f) isTRUE(f$singular), logical(1))),
    n_vcov_from_RX = n_vcov_fallback,
    n_events = if (!is.null(family)) sum(getME(fits[ok][[1]]$fit, "y") == 1) else NA_integer_,
    nobs = nobs(fits[ok][[1]]$fit),
    warnings = paste(unique(unlist(lapply(fits, `[[`, "warnings"))), collapse = " | "))
  list(pooled = pooled, fits = fits[ok], re_sd = mean(re_sd), diag = diag)
}

# Odds-ratio style output from a pooled table
exp_table <- function(pooled, terms = NULL) {
  p <- copy(pooled)
  if (!is.null(terms)) p <- p[term %in% terms]
  p[, `:=`(ratio = exp(estimate), ratio_lower = exp(lower), ratio_upper = exp(upper))]
  p[]
}

# Median odds ratio and latent-scale ICC from a random-intercept SD
mor_from_sd <- function(s) exp(sqrt(2 * s^2) * qnorm(0.75))
icc_logit   <- function(s) s^2 / (s^2 + base::pi^2 / 3)

# Events per parameter for a binary outcome
epv <- function(y, n_par) {
  ev <- min(sum(y == 1, na.rm = TRUE), sum(y == 0, na.rm = TRUE))
  ev / n_par
}

# Number of fixed-effect parameters implied by a formula on a data set
n_fixed_par <- function(formula, data) {
  f <- lme4::nobars(formula)
  ncol(model.matrix(f, model.frame(f, data))) - 1
}

# ---- g-computation of standardised risks by year ----------------------------------
# Predicted risk for every patient with calendar year set to each value, averaged
# over the observed cohort (case mix and centre mix). Uncertainty: fixed effects
# drawn from their sampling distribution, centre effects held at their conditional
# modes (approximation; documented in the methods report).
gcomp_year <- function(fit, data, years, nsim = CONFIG$n_sim) {
  tt <- delete.response(terms(fit, fixed.only = TRUE))
  beta <- fixef(fit); V <- as.matrix(vcov(fit))
  draws <- MASS::mvrnorm(nsim, beta, V)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
  d0 <- model.frame(fit)
  out <- lapply(years, function(y) {
    nd <- copy(as.data.table(data)); nd[, year_c := y - CONFIG$reference_year]
    re_off <- predict(fit, newdata = nd, re.form = NULL, type = "link", allow.new.levels = TRUE) -
              predict(fit, newdata = nd, re.form = NA, type = "link")
    xl <- .xlevels_merMod(fit); xl <- xl[intersect(names(xl), all.vars(tt))]
    X <- model.matrix(tt, model.frame(tt, nd, na.action = na.pass, xlev = xl),
                      contrasts.arg = attr(getME(fit, "X"), "contrasts"))
    X <- X[, names(beta), drop = FALSE]
    pt <- mean(plogis(X %*% beta + re_off), na.rm = TRUE)
    sims <- apply(draws, 1, function(b) mean(plogis(X %*% b + re_off), na.rm = TRUE))
    data.table(year = y, risk = pt, sims = list(sims))
  })
  rbindlist(out)
}

.xlevels_merMod <- function(fit) {
  fr <- model.frame(fit)
  fac <- names(fr)[vapply(fr, is.factor, logical(1))]
  setNames(lapply(fac, function(v) levels(fr[[v]])), fac)
}

# ---- survival helpers ---------------------------------------------------------------

km_at <- function(sf, times) {
  s <- summary(sf, times = times, extend = TRUE)
  data.table(strata = if (is.null(s$strata)) "Overall" else as.character(s$strata),
             time = s$time, n_risk = s$n.risk, n_event_cum = NA_integer_,
             surv = s$surv, lower = s$lower, upper = s$upper)
}

# Reverse Kaplan-Meier median follow-up
reverse_km_median <- function(time, status) {
  sf <- survfit(Surv(time, 1 - status) ~ 1)
  q <- quantile(sf, probs = c(.25, .5, .75), conf.int = FALSE)
  as.numeric(q)
}

# RMST up to tau for each stratum (from survfit) with SE
rmst_by <- function(sf, tau) {
  tab <- summary(sf, rmean = tau)$table
  if (is.null(dim(tab))) tab <- t(as.matrix(tab))
  data.table(strata = if (nrow(tab) == 1 && is.null(rownames(tab))) "Overall" else rownames(tab),
             rmst = tab[, "rmean"], se = tab[, "se(rmean)"])
}

# Exact Poisson CI for an SMR
smr_ci <- function(obs, exp) {
  lo <- ifelse(obs == 0, 0, qchisq(0.025, 2 * obs) / 2) / exp
  hi <- qchisq(0.975, 2 * (obs + 1)) / 2 / exp
  data.table(observed = obs, expected = exp, smr = obs / exp, lower = lo, upper = hi)
}

# Custom link for Poisson excess-hazard regression (Dickman et al. 2004)
excess_link <- function(d_star) {
  linkfun  <- function(mu) log(pmax(mu - d_star, 1e-10))
  linkinv  <- function(eta) d_star + exp(eta)
  mu.eta   <- function(eta) exp(eta)
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                 valideta = valideta, name = "log(mu - d*)"), class = "link-glm")
}

# ---- misc ------------------------------------------------------------------------------

anonymise_centres <- function(centre, n_by_centre = NULL) {
  # Letters assigned in order of decreasing volume (ties broken by original id)
  tab <- as.data.table(table(centre = as.character(centre)))
  setorder(tab, -N, centre)
  map <- setNames(LETTERS[seq_len(nrow(tab))], tab$centre)
  map
}

holm <- function(p) p.adjust(p, "holm")
bh   <- function(p) p.adjust(p, "BH")

write_session_info <- function() {
  f <- file.path(out_dir("."), "session_info.txt")
  sink(f); print(sessionInfo()); cat("\nRun at:", format(Sys.time()), "\n"); sink()
  invisible(f)
}

# ---- multiple imputation ---------------------------------------------------------------

IMP_VARS_CONT <- c("age", "bmi", "egfr", "lvef", "euroscore2")
IMP_VARS_BIN  <- c("female", "prior_cs", "endocarditis_b", "conc_cabg", "conc_valve",
                   "root_inv", "desc_inv", "hybrid_fet", "other_cardiac")
IMP_VARS_FAC  <- c("urg4", "path3")
IMP_AUX       <- c("comp_primary", "d30", "dzh", "major_compl", "status", "H")   # outcomes: predictors only

# Recreate derived categories after imputation
rederive <- function(d) {
  d <- copy(d)
  d[, lvef_cat := relevel(cut(lvef, c(0, 20, 30, 50, 100), c("<=20%", "21-30%", "31-50%", ">50%")), ">50%")]
  d[, urg3 := factor(fcase(urg4 == "Elective", "Elective", urg4 == "Urgent", "Urgent",
                           urg4 %in% c("Emergency", "Salvage"), "Emergency/salvage"),
                     c("Elective", "Urgent", "Emergency/salvage"))]
  d[, nonelective := as.integer(urg4 != "Elective")]
  d[, emerg_salv := as.integer(urg4 %in% c("Emergency", "Salvage"))]
  # urgency2 / acute (acute vs non-acute) are kept as observed: not re-derived from imputed values
  d[, acute_dissection := as.integer(path3 == "Acute dissection")]
  d[, logit_es2 := qlogis(pmin(pmax(euroscore2 / 100, 1e-4), 1 - 1e-4))]
  d[]
}

# Multiple imputation by chained equations. Outcomes are used as predictors but
# their observed values are restored afterwards (never analysed as imputed).
# delta/delta_var: shift imputed values of one continuous variable (MNAR sensitivity).
make_imputations <- function(cohort, m = CONFIG$n_imputations, delta = 0, delta_var = NULL,
                             maxit = 10, seed = CONFIG$seed) {
  d <- copy(cohort)
  ok <- !is.na(d$time_d) & !is.na(d$status)
  d[, H := NA_real_]
  d[ok, H := mice::nelsonaalen(as.data.frame(d[ok]), time_d, status)]
  vars <- c(IMP_VARS_CONT, IMP_VARS_BIN, IMP_VARS_FAC, IMP_AUX, "year_c", "centrum")
  dd <- as.data.frame(d[, ..vars])
  dd$centrum <- factor(dd$centrum)
  meth <- mice::make.method(dd)
  meth[c(IMP_VARS_CONT, IMP_VARS_BIN, "H", "comp_primary", "d30", "dzh", "major_compl", "status")] <- "pmm"
  meth[IMP_VARS_FAC] <- "polyreg"
  meth[!vapply(dd, anyNA, logical(1))] <- ""
  post <- mice::make.post(dd)
  if (!is.null(delta_var) && delta != 0) {
    post[delta_var] <- sprintf("imp[[j]][, i] <- imp[[j]][, i] + %s", delta)
  }
  imp <- suppressWarnings(mice::mice(dd, m = m, method = meth, post = post, maxit = maxit,
                                     seed = seed, printFlag = FALSE))
  out <- lapply(seq_len(m), function(i) {
    ci <- as.data.table(mice::complete(imp, i))
    di <- copy(d)
    for (v in c(IMP_VARS_CONT, IMP_VARS_BIN, IMP_VARS_FAC)) set(di, j = v, value = ci[[v]])
    rederive(di)
  })
  attr(out, "mids") <- imp
  out
}

# ---- calendar-year restricted cubic spline (fixed knots) ---------------------------------
# Linear term year_c plus k-2 non-linear terms (Harrell's parameterisation), so
# non-linearity is tested by a Wald test of the non-linear terms only.
# Knots fixed a priori at 2014, 2017, 2020, 2023 (approx. 5th/35th/65th/95th
# percentiles of a 2013-2024 cohort with rising volume).
YEAR_KNOTS <- c(2014, 2017, 2020, 2023) - CONFIG$reference_year

rcs_nonlin <- function(x, knots = YEAR_KNOTS) {
  k <- length(knots); norm <- (knots[k] - knots[1])^2
  m <- vapply(seq_len(k - 2), function(j) {
    (pmax(x - knots[j], 0)^3 -
       pmax(x - knots[k - 1], 0)^3 * (knots[k] - knots[j]) / (knots[k] - knots[k - 1]) +
       pmax(x - knots[k], 0)^3 * (knots[k - 1] - knots[j]) / (knots[k] - knots[k - 1])) / norm
  }, numeric(length(x)))
  matrix(m, ncol = k - 2)
}
ynl1 <- function(x) rcs_nonlin(x)[, 1]
ynl2 <- function(x) rcs_nonlin(x)[, 2]
YEAR_SPLINE <- "year_c + ynl1(year_c) + ynl2(year_c)"
YEAR_NL_TERMS <- c("ynl1(year_c)", "ynl2(year_c)")

# Multivariate Wald test pooled over imputations (D1; Li, Raghunathan & Rubin 1991)
pool_wald_D1 <- function(est_list, vcov_list, terms) {
  m <- length(est_list); k <- length(terms)
  Q <- do.call(rbind, lapply(est_list, function(e) e[terms]))
  qbar <- colMeans(Q)
  Ubar <- Reduce(`+`, lapply(vcov_list, function(v) as.matrix(v)[terms, terms])) / m
  if (m == 1) {
    stat <- as.numeric(t(qbar) %*% solve(Ubar) %*% qbar) / k
    return(data.table(F = stat, df1 = k, df2 = Inf, p = pchisq(stat * k, k, lower.tail = FALSE)))
  }
  B <- cov(Q)
  r <- (1 + 1 / m) * sum(diag(B %*% solve(Ubar))) / k
  D1 <- as.numeric(t(qbar) %*% solve(Ubar) %*% qbar) / (k * (1 + r))
  v <- k * (m - 1)
  w <- if (v > 4) 4 + (v - 4) * (1 + (1 - 2 / v) / r)^2 else v * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  data.table(F = D1, df1 = k, df2 = w, p = pf(D1, k, w, lower.tail = FALSE))
}

# g-computation across imputations: pools point estimates and simulation draws
gcomp_mi <- function(fits, datasets, years, nsim_total = CONFIG$n_sim) {
  per <- max(20, ceiling(nsim_total / length(fits)))
  res <- lapply(seq_along(fits), function(i) gcomp_year(fits[[i]], datasets[[i]], years, per))
  rbindlist(lapply(seq_along(years), function(j) {
    pts <- vapply(res, function(r) r$risk[j], numeric(1))
    sims <- unlist(lapply(res, function(r) r$sims[[j]]))
    data.table(year = years[j], risk = mean(pts),
               lower = quantile(sims, .025), upper = quantile(sims, .975))
  }))
}

# Likelihood-ratio (profile) 95% CI for a single random-intercept SD in a GLMM.
# Uses lme4's deviance function (nAGQ = 0: fixed effects profiled out by PIRLS at
# each SD value) and root-finding for deviance = minimum + 3.84.
profile_sd_ci <- function(formula, data, family = binomial(), level = 0.95) {
  dfun <- glmer(formula, data = data, family = family, nAGQ = 0, devFunOnly = TRUE)
  crit <- qchisq(level, 1)
  opt <- optimize(function(t) dfun(t), c(0, 5))
  th <- opt$minimum; dmin <- opt$objective
  if (dfun(0) <= dmin + 1e-8) { th <- 0; dmin <- dfun(0) }
  g <- function(t) dfun(t) - dmin - crit
  lo <- if (th == 0 || g(0) < 0) 0 else uniroot(g, c(0, th))$root
  hi_end <- max(1, 3 * th); while (g(hi_end) < 0 && hi_end < 50) hi_end <- hi_end * 2
  hi <- uniroot(g, c(th, hi_end))$root
  c(sd = th, lower = lo, upper = hi)
}

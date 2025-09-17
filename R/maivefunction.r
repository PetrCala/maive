#' R code for MAIVE
#'
#' R package for MAIVE: "Spurious Precision in Meta-Analysis of Observational Research" by
#' Zuzana Irsova, Pedro Bom, Tomas Havranek, and Heiko Rachinger.
#'
#' @param dat Data frame with columns bs, sebs, Ns, study_id (optional).
#' @param method 1 FAT-PET, 2 PEESE, 3 PET-PEESE, 4 EK.
#' @param weight 0 no weights, 1 standard weights, 2 adjusted weights.
#' @param instrument 1 yes, 0 no.
#' @param studylevel Correlation at study level: 0 none, 1 fixed effects, 2 cluster.
#' @param SE SE estimator: 0 CR0 (Huber–White), 1 CR1 (Standard empirical correction),
#' 2 CR2 (Bias-reduced estimator), 3 wild bootstrap.
#' @param AR Anderson Rubin corrected CI for weak instruments (only for unweighted MAIVE versions
#' of PET, PEESE, PET-PEESE, not available for fixed effects): 0 no, 1 yes.
#'
#' @details Data \code{dat} can be imported from an Excel file via:
#' \code{dat <- read_excel("inputdata.xlsx")} or from a csv file via: \code{dat <- read.csv("inputdata.csv")}
#' It should contain:
#' \itemize{
#'   \item Estimates: bs
#'   \item Standard errors: sebs
#'   \item Number of observations: Ns
#'   \item Optional: study_id
#' }
#' Default option for MAIVE: MAIVE-PET-PEESE, unweighted, instrumented, cluster SE, wild bootstrap, AR.
#'
#' @return \itemize{
#'   \item beta: MAIVE meta-estimate
#'   \item SE: MAIVE standard error
#'   \item F-test: heteroskedastic robust F-test of the first step instrumented SEs
#'   \item beta_standard: point estimate from the method chosen
#'   \item SE_standard: standard error from the method chosen
#'   \item Hausman: Hausman type test: comparison between MAIVE and standard version
#'   \item Chi2: 5% critical value for Hausman test
#'   \item SE_instrumented: instrumented standard errors
#'   \item AR_CI: Anderson-Rubin confidence interval for weak instruments
#'   \item pub bias p-value: p-value of test for publication bias / p-hacking based on instrumented FAT
#'   \item egger_coef: Egger Coefficient (PET estimate)
#'   \item egger_se: Egger Standard Error (PET standard error)
#'   \item is_quadratic_fit: Details on quadratic selection and slope behaviour
#'   \item boot_result: Boot result
#'   \item slope_coef: Slope coefficient
#' }
#'
#' @export
maive <- function(dat, method, weight, instrument, studylevel, SE, AR) {
  # Manual wild cluster bootstrap function:
  # wild bootstrap
  # clustered (weights drawn per cluster, not per observation)
  # with Rademacher distribution (i.e., ±1 with 0.5 probability each)

  methods <- c("PET", "PEESE", "PET-PEESE", "EK")
  instrumented <- c("not instrumented", "instrumented")
  weighted <- c("no weights", "standardly weighted", "adjusted weights")
  studylevelcorrelation <- c("none", "study level dummies", "cluster")

  # cluster for studylevel >=2, dummy for studylevel 1 or 3
  cluster <- studylevel %/% 2L
  dummy <- studylevel %% 2L

  type_map <- c("CR0", "CR1", "CR2")
  type_choice <- if (SE == 3L) "CR0" else type_map[SE + 1L] # not needed if bootstrap

  # AR not available for EK, for fixed effects, and for weighted
  if (method == 4L || weight %in% c(1L, 2L) || instrument == 0L) AR <- 0L

  # extracting data from excel
  dat <- as.data.frame(dat)
  bs <- dat[, 1]
  M <- length(bs)
  sebs <- dat[, 2]
  Ns <- dat[, 3]

  if (ncol(dat) >= 4) {
    studyid <- dat[, 4]
  } else {
    studyid <- (1:M)
    dummy <- 0L
    cluster <- 0L
  }

  alpha_s <- 0.05

  # create Dummies from studyid
  build_dummy_matrix <- function(values) {
    if (exists("to.dummy", mode = "function")) {
      to.dummy(data.frame(studyid = values), "studyid")
    } else {
      f <- factor(values)
      mm <- stats::model.matrix(~ f - 1)
      colnames(mm) <- paste0("studyid_", levels(f))
      mm
    }
  }
  D <- build_dummy_matrix(studyid)
  D <- D - matrix(colMeans(D), nrow = M, ncol = ncol(D), byrow = TRUE)
  D <- D[, 1:(dim(D)[2] - 1)]

  # g=studyid if clustered and g=(1:M)' if not clustered (gives heteroskedastic robust SE)
  g <- if (cluster == 0L) seq_len(M) else studyid
  dat$g <- g

  # (1) Instrumenting the variances with sample size allowing for a constant and including dummies
  invNs <- 1 / Ns
  sebs2 <- sebs^2
  Xiv <- cbind(1, invNs)
  varreg1 <- lm(sebs2 ~ 0 + Xiv)
  dimiv <- 2
  if (varreg1$coefficients[1] < 0) {
    Xiv <- as.matrix(invNs)
    varreg1 <- lm(sebs2 ~ 0 + Xiv)
    dimiv <- 1
  }

  sebs2fit1 <- fitted(varreg1)

  # F-statistic of first step. heteroskedasticity and autocorrelation robust variance HAC
  F_hac <- (varreg1$coefficients[dimiv]^2 / vcovCR(varreg1, cluster = g, type = type_choice)[dimiv, dimiv])
  F_hac <- if (instrument == 0L) "NA" else round(F_hac, 3)

  # ---- weights (implemented via division so that WLS is OLS on transformed data)
  if (weight == 0) {
    w <- rep(1, M)
  } else if (weight == 1) {
    w <- sebs
  } else if (weight == 2) {
    w <- sqrt(sebs2fit1)
  } else {
    stop("Invalid weight")
  }

  # ---- PET regressors (x) and PEESE regressors (x2)
  if (instrument == 0L) {
    x <- sebs
    x2 <- sebs^2
  } else {
    x <- sqrt(sebs2fit1) # instrumented SE (MAIVE)
    x2 <- sebs2fit1 # instrumented variance
  }

  # ---- Build design matrices honoring weights and study FE choice
  y <- bs / w
  X <- cbind(1, x) / w
  X2 <- cbind(1, x2) / w

  # Baseline (for comparisons): inverse-variance weighted, not instrumented
  y0 <- bs / sebs
  X0 <- cbind(1, sebs) / sebs
  X20 <- cbind(1, sebs^2) / sebs

  # Add study FEs if requested
  if (dummy == 1L) {
    X <- cbind(X, D / w)
    X2 <- cbind(X2, D / w)
    cD <- cbind(1, D) / w # for EK intercept + FE
    X0 <- cbind(X0, D / sebs)
    X20 <- cbind(X20, D / sebs)
    cD0 <- cbind(1, D) / sebs
  } else {
    cD <- matrix(1 / w, ncol = 1)
    cD0 <- matrix(1 / sebs, ncol = 1)
  }

  # ---- FE objects
  wis0 <- 1 / (w^2)
  fe <- sum(bs * wis0) / sum(wis0)
  varfe <- 1 / sum(wis0)

  wis00 <- 1 / (sebs^2)
  fe0 <- sum(bs * wis00) / sum(wis00)
  varfe0 <- 1 / sum(wis00)

  # ---- WLS levels (for sigma_h^2 later)
  wlsreg <- lm(y ~ 0 + cD)
  wlsreg0 <- lm(y0 ~ 0 + cD0)

  # ---- PET/PEESE fits (MAIVE & baseline)
  fatpet <- lm(y ~ 0 + X) # PET
  peese <- lm(y ~ 0 + X2) # PEESE
  fatpet0 <- lm(y0 ~ 0 + X0) # baseline PET
  peese0 <- lm(y0 ~ 0 + X20) # baseline PEESE


  # ---- PET-PEESE selection (uses standard OLS SEs for the intercept test)
  quadratic_decision <- abs(coef(fatpet)[1] / sqrt(vcov(fatpet)[1, 1])) >
    qt(1 - alpha_s / 2, M - ncol(X))
  petpeese <- if (quadratic_decision) peese else fatpet

  is_quadratic_fit0 <- abs(coef(fatpet0)[1] / sqrt(vcov(fatpet0)[1, 1])) >
    qt(1 - alpha_s / 2, M - ncol(X0))
  petpeese0 <- if (is_quadratic_fit0) peese0 else fatpet0

  # ---- sigma_h^2
  Qfe0 <- sum(residuals(wlsreg)^2)
  sigh2hat0 <- max(0, M * ((Qfe0 / (M - ncol(model.matrix(wlsreg)) - 1)) - 1) / sum(wis0))
  sighhat0 <- sqrt(sigh2hat0)

  Qfe00 <- sum(residuals(wlsreg0)^2)
  sigh2hat00 <- max(0, M * ((Qfe00 / (M - ncol(model.matrix(wlsreg0)) - 1)) - 1) / sum(wis00))
  sighhat00 <- sqrt(sigh2hat00)

  # Endogenous Kink (EK) Threshold- MAIVE
  if (coef(petpeese)[1] > 1.96 * sighhat0) {
    a0 <- (coef(petpeese)[1] - 1.96 * sighhat0) * (coef(petpeese)[1] + 1.96 * sighhat0) /
      (2 * 1.96 * coef(petpeese)[1])
  } else {
    a0 <- 0
  }

  # Endogenous Kink (EK) Threshold - baseline
  if (coef(petpeese0)[1] > 1.96 * sighhat00) {
    a00 <- (coef(petpeese0)[1] - 1.96 * sighhat00) * (coef(petpeese0)[1] + 1.96 * sighhat00) /
      (2 * 1.96 * coef(petpeese0)[1])
  } else {
    a00 <- 0
  }

  # EK - MAIVE
  ek_structure <- "intercept"
  if (!is.na(a0) && a0 > min(x) && a0 < max(x)) {
    xx_w <- (x - a0) * (x > a0) / w
    ekreg <- lm(y ~ 0 + cD + xx_w)
    ek_structure <- "kink"
  } else if (!is.na(a0) && a0 < min(x)) {
    x_w <- x / w
    ekreg <- lm(y ~ 0 + cD + x_w)
    ek_structure <- "linear"
  } else {
    ekreg <- lm(y ~ 0 + cD)
  }
  ek <- coef(ekreg)[1]

  # EK - baseline
  if (a00 > min(sebs) && a00 < max(sebs)) {
    xx0_w <- (sebs - a00) * (sebs > a00) / sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + xx0_w)
  } else if (a00 < min(sebs)) {
    x0_w <- sebs / sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + x0_w)
  } else {
    ekreg0 <- lm(y0 ~ 0 + cD0)
  }
  ek0 <- coef(ekreg0)[1]

  # ---------- Helper: robust/clustered vs bootstrap inference for a given coefficient
  infer_coef <- function(model, coef_index, SE, data, cluster_var, type_choice) {
    if (SE == 3L) {
      boot <- manual_wild_cluster_boot_se(
        model = model,
        data = data,
        cluster_var = "g",
        B = 500
      )
      # If your bootstrap helper returns p-values per coefficient, prefer those:
      if (!is.null(boot$boot_se)) {
        se <- boot$boot_se[coef_index]
      } else {
        stop("Bootstrap helper must return boot_se.")
      }
      boot_result <- boot
    } else {
      V <- vcovCR(model, cluster = data[[cluster_var]], type = type_choice)
      se <- sqrt(V[coef_index, coef_index])
      boot_result <- NULL
    }
    b <- coef(model)[coef_index]
    t <- as.numeric(b / se)
    p <- 2 * pnorm(-abs(t))
    list(b = b, se = se, p = p, boot_result = boot_result)
  }

  # ---------- Map: which model is “selected” for reporting
  slope_info <- switch(as.character(method),
    "1" = list(
      type = "linear",
      coefficient = round(as.numeric(coef(fatpet)[2]), 3),
      detail = NULL
    ),
    "2" = list(
      type = "quadratic",
      coefficient = round(as.numeric(coef(peese)[2]), 3),
      detail = NULL
    ),
    "3" = if (identical(petpeese, peese)) {
      list(
        type = "quadratic",
        coefficient = round(as.numeric(coef(peese)[2]), 3),
        detail = NULL
      )
    } else {
      list(
        type = "linear",
        coefficient = round(as.numeric(coef(fatpet)[2]), 3),
        detail = NULL
      )
    },
    "4" = {
      if (ek_structure == "kink") {
        kink_effect <- round(as.numeric(tail(coef(ekreg), 1)), 3)
        kink_location <- as.numeric(round(a0, 3))
        list(
          type = "kinked",
          coefficient = list(
            kink_effect = kink_effect,
            kink_location = kink_location
          ),
          detail = list(
            kink_location = kink_location,
            kink_effect = kink_effect
          )
        )
      } else if (ek_structure == "linear") {
        list(
          type = "linear",
          coefficient = round(as.numeric(tail(coef(ekreg), 1)), 3),
          detail = NULL
        )
      } else {
        list(
          type = "linear",
          coefficient = 0,
          detail = NULL
        )
      }
    },
    stop("Invalid method")
  )

  slope_coef <- slope_info$coefficient

  is_quadratic_summary <- list(
    quadratic = switch(as.character(method),
      "1" = FALSE,
      "2" = TRUE,
      "3" = isTRUE(quadratic_decision),
      "4" = FALSE,
      stop("Invalid method")
    ),
    slope_type = slope_info$type,
    slope_detail = slope_info$detail
  )

  # ---------- PET (“Egger”) outputs that respect user settings
  egger_inf <- infer_coef(fatpet,
    coef_index = 2L, SE = SE, data = dat,
    cluster_var = "g", type_choice = type_choice
  )
  egger_coef <- round(as.numeric(egger_inf$b), 3)
  egger_se <- round(as.numeric(egger_inf$se), 3)
  pb_p <- round(egger_inf$p, 3)

  # ---------- SE for intercepts of chosen MAIVE vs baseline models
  get_se_int <- function(model) {
    inf <- infer_coef(model,
      coef_index = 1L, SE = SE, data = dat,
      cluster_var = "g", type_choice = type_choice
    )
    list(se = inf$se, boot_result = inf$boot_result)
  }

  # Map method -> (MAIVE model, Standard model)
  cfg_map <- list(
    "1" = list(maive = fatpet, std = fatpet0),
    "2" = list(maive = peese, std = peese0),
    "3" = list(maive = petpeese, std = petpeese0),
    "4" = list(maive = ekreg, std = ekreg0)
  )
  cfg <- cfg_map[[as.character(method)]]
  if (is.null(cfg)) stop("Invalid method")

  beta <- coef(cfg$maive)[1]
  se_ma <- get_se_int(cfg$maive)
  betase <- se_ma$se
  boot_result <- se_ma$boot_result

  beta0 <- coef(cfg$std)[1]
  se_st <- get_se_int(cfg$std)
  beta0se <- se_st$se

  # Hausman (conservative: denom = Var(beta_MAIVE))
  V_maive <- vcovCR(cfg$maive, cluster = g, type = type_choice)
  Hausman <- (beta - beta0)^2 / V_maive[1, 1]
  Chi2 <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)

  # Anderson and Rubin Confidence intervals
  get_ar_ci_res <- function() {
    # Either AR disabled or EK (no AR CI)
    if (AR != 1L || method == 4L || weight %in% c(1L, 2L)) {
      return(list(b0_CI = "NA", b1_CI = "NA"))
    }
    cfg_ar <- switch(as.character(method),
      "1" = list(model = fatpet, adjust_fun = PET_adjust),
      "2" = list(model = peese, adjust_fun = PEESE_adjust),
      "3" = if (identical(petpeese, peese)) {
        list(model = peese, adjust_fun = PEESE_adjust)
      } else {
        list(model = fatpet, adjust_fun = PET_adjust)
      }
    )
    do.call(
      compute_AR_CI_optimized,
      c(cfg_ar, list(bs = bs, sebs = sebs, invNs = invNs, g = g, type_choice = type_choice))
    )
  }
  ar_ci_res <- get_ar_ci_res()

  list(
    "beta"              = round(beta, 3),
    "SE"                = round(betase, 3),
    "F-test"            = F_hac,
    "beta_standard"     = round(beta0, 3),
    "SE_standard"       = round(beta0se, 3),
    "Hausman"           = round(Hausman, 3),
    "Chi2"              = round(Chi2, 3),
    "SE_instrumented"   = sqrt(sebs2fit1),
    "AR_CI"             = ar_ci_res$b0_CI,
    "pub bias p-value"  = pb_p,
    "egger_coef"        = egger_coef,
    "egger_se"          = egger_se,
    "is_quadratic_fit"  = is_quadratic_summary,
    "boot_result"       = boot_result, # NULL unless SE == 3 and requested for intercept
    "slope_coef"        = slope_coef
  )
}

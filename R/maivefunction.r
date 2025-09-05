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
#'   \item MAIVE meta-estimate and standard error
#'   \item Hausman type test: comparison between MAIVE and standard version
#'   \item When instrumenting: heteroskedastic robust F-test of the first step instrumented SEs
#'   \item p-value of test for publication bias / p-hacking based on instrumented FAT
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

  if (studylevel == 0) {
    cluster <- 0
    dummy <- 0
  } else if (studylevel == 1) {
    cluster <- 0
    dummy <- 1
  } else if (studylevel == 2) {
    cluster <- 1
    dummy <- 0
  } else if (studylevel == 3) {
    cluster <- 1
    dummy <- 1
  }
  type_map <- c("CR0", "CR1", "CR2")
  if (SE < 3) {
    type_choice <- type_map[SE + 1]
  } else if (SE == 3) { # not needed if bootstrap
    type_choice <- type_map[0 + 1]
  }

  # AR not available for EK, for fixed effects, and for weighted
  if (method == 4 | weight == 1 | weight == 2 | instrument == 0) {
    AR <- 0
  }

  # extracting data from excel
  dat <- as.data.frame(dat)
  bs <- dat[, 1]
  M <- length(bs)
  sebs <- dat[, 2]
  Ns <- dat[, 3]

  if (dim(dat)[2] == 4) {
    studyid <- dat[, 4]
  } else {
    studyid <- (1:M)
    dummy <- 0
    cluster <- 0
  }

  alpha_s <- 0.05

  # create Dummies from studyid
  df <- data.frame(studyid)
  D <- to.dummy(df, "studyid")
  D <- D - matrix(colMeans(D), nrow = M, ncol = size(D)[2], byrow = TRUE)
  D <- D[, 1:(dim(D)[2] - 1)]

  # g=studyid if clustered and g=(1:M)' if not clustered (gives heteroskedastic robust SE)
  if (cluster == 0) {
    g <- (1:M)
  } else if (cluster == 1) {
    g <- studyid
  }

  dat$g <- g

  # (1) Instrumenting the variances with sample size allowing for a constant and including dummies
  invNs <- 1 / Ns
  sebs2 <- sebs^2
  Xiv <- matrix(c(ones(M, 1)[, 1], invNs), nrow = M)
  varreg1 <- lm(sebs2 ~ 0 + Xiv)
  dimiv <- 2
  if (varreg1$coefficients[1] < 0) {
    Xiv <- invNs
    varreg1 <- lm(sebs2 ~ 0 + Xiv)
    dimiv <- 1
  }

  sebs2fit1 <- varreg1$fitted.values

  # F-statistic of first step. heteroskedasticity and autocorrelation robust variance HAC
  F_hac <- (varreg1$coefficients[dimiv]^2 / vcovCR(varreg1, cluster = g, type = type_choice)[dimiv, dimiv])

  # weight
  if (weight == 0) {
    w <- ones(M, 1)[, 1]
  } else if (weight == 1) {
    w <- sebs
  } else if (weight == 2) {
    w <- sebs2fit1^(1 / 2)
  }

  # instrument
  if (instrument == 0) {
    x <- sebs
    x2 <- sebs^2
    F_hac <- "NA"
  } else if (instrument == 1) {
    x <- sebs2fit1^(1 / 2)
    x2 <- sebs2fit1
    F_hac <- round(F_hac, 3)
  }

  # choose dependent variable and regressor
  y <- bs / w
  x <- x
  x2 <- x2
  X <- matrix(c(ones(M, 1)[, 1], x) / w, nrow = M)
  X_d <- matrix(c(X, D / w), nrow = M)
  X2 <- matrix(c(ones(M, 1)[, 1], x2) / w, nrow = M)
  X2_d <- matrix(c(X2, D / w), nrow = M)

  # baseline, i.e. chosen method, with chosen options of study-level correlation
  #  but with inverse-variance weighting and without instrumenting
  y0 <- bs / sebs
  x0 <- sebs
  x20 <- sebs^2
  X0 <- matrix(c(ones(M, 1)[, 1], x0) / sebs, nrow = M)
  X0_d <- matrix(c(X0, D / sebs), nrow = M)
  X20 <- matrix(c(ones(M, 1)[, 1], x20) / sebs, nrow = M)
  X20_d <- matrix(c(X20, D / sebs), nrow = M)

  if (dummy == 0) {
    X <- X
    X0 <- X0
    X2 <- X2
    X20 <- X20
    cD <- ones(M, 1)[, 1]
  } else if (dummy == 1) {
    X <- X_d
    X0 <- X0_d
    X2 <- X2_d
    X20 <- X20_d
    cD <- matrix(c(ones(M, 1)[, 1], D), nrow = M) # for EK stack constant and dummies
  }

  cD <- cD / w
  cD0 <- cD / sebs

  ones_w <- ones(M, 1)[, 1] / w
  ones_w0 <- ones(M, 1)[, 1] / sebs

  # Fixed effects (FE)
  wis0 <- 1 / (w^2)
  fe <- sum(bs * wis0) / sum(wis0)
  varfe <- 1 / sum(wis0)
  # baseline
  wis00 <- 1 / (sebs^2)
  fe0 <- sum(bs * wis00) / sum(wis00)
  varfe0 <- 1 / sum(wis00)

  # WLS
  wlsreg <- lm(y ~ 0 + cD)
  wls <- wlsreg$coefficients[1]
  wlsse <- sqrt(vcovCR(wlsreg, cluster = g, type = type_choice)[1, 1])
  # baseline
  wlsreg0 <- lm(y0 ~ 0 + cD0)
  wls0 <- wlsreg0$coefficients[1]
  wlsse0 <- sqrt(vcovCR(wlsreg0, cluster = g, type = type_choice)[1, 1])

  # FAT-PET - MAIVE
  fatpet <- lm(y ~ 0 + X)
  # FAT-PET - baseline case
  fatpet0 <- lm(y0 ~ 0 + X0)

  # PEESE - MAIVE
  peese <- lm(y ~ 0 + X2)
  # PEESE - baseline case
  peese0 <- lm(y0 ~ 0 + X20)


  # PET-PEESE - MAIVE
  is_quadratic_fit <- abs(fatpet$coefficients[1] / sqrt(vcovCR(fatpet, cluster = g, type = type_choice)[1, 1])) > qt(1 - alpha_s / 2, M - dim(X)[2] - 1)
  if (is_quadratic_fit) {
    petpeese <- peese
  } else {
    petpeese <- fatpet
  }
  # PET-PEESE - baseline case
  is_quadratic_fit0 <- abs(fatpet0$coefficients[1] / sqrt(vcovCR(fatpet0, cluster = g, type = type_choice)[1, 1])) > qt(1 - alpha_s / 2, M - dim(X0)[2] - 1)
  if (is_quadratic_fit0) {
    petpeese0 <- peese0
  } else {
    petpeese0 <- fatpet0
  }

  # True effect variance - MAIVE
  Qfe0 <- sum(wlsreg$residuals * wlsreg$residuals)
  sigh2hat0 <- max(0, M * ((Qfe0 / (M - dim(wlsreg$model)[2] - 1)) - 1) / sum(wis0))
  sighhat0 <- sqrt(sigh2hat0)
  # True effect variance - baseline
  Qfe00 <- sum(wlsreg0$residuals * wlsreg0$residuals)
  sigh2hat00 <- max(0, M * ((Qfe00 / (M - dim(wlsreg0$model)[2] - 1)) - 1) / sum(wis00))
  sighhat00 <- sqrt(sigh2hat00)

  # Endogenous Kink (EK) Threshold- MAIVE
  if (petpeese$coefficients[1] > 1.96 * sighhat0) {
    a0 <- (petpeese$coefficients[1] - 1.96 * sighhat0) * (petpeese$coefficients[1] + 1.96 * sighhat0) / (2 * 1.96 * petpeese$coefficients[1])
  } else {
    a0 <- 0
  }
  # Endogenous Kink (EK) Threshold - baseline
  if (petpeese0$coefficients[1] > 1.96 * sighhat00) {
    a00 <- (petpeese0$coefficients[1] - 1.96 * sighhat00) * (petpeese0$coefficients[1] + 1.96 * sighhat00) / (2 * 1.96 * petpeese0$coefficients[1])
  } else {
    a00 <- 0
  }

  # EK - MAIVE
  if (a0 > min(x) && a0 < max(x)) {
    xx_w <- (x - a0) * (x > a0) / w
    ekreg <- lm(y ~ 0 + cD + xx_w)
  } else if (a0 < min(x)) {
    x_w <- x / w
    ekreg <- lm(y ~ 0 + cD + x_w)
  } else if (a0 > max(x)) {
    ekreg <- lm(y ~ 0 + cD)
  }
  ek <- ekreg$coefficients[1]

  # EK - baseline
  if (a00 > min(x0) && a00 < max(x0)) {
    xx0_w <- (x0 - a00) * (x0 > a00) / sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + xx0_w)
  } else if (a00 < min(x0)) {
    x0_w <- x0 / sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + x0_w)
  } else if (a00 > max(x0)) {
    ekreg0 <- lm(y0 ~ 0 + cD0)
  }
  ek0 <- ekreg0$coefficients[1]

  # Anderson and Rubin Confidence intervals
  get_ar_ci_res <- function() {
    # Either AR disabled or EK (no AR CI)
    if (AR != 1 || method == 4 || weight == 1 || weight == 2) {
      return(list(b0_CI = "NA", b1_CI = "NA"))
    }

    cfg <- switch(as.character(method),
      "1" = list(model = fatpet, adjust_fun = PET_adjust),
      "2" = list(model = peese, adjust_fun = PEESE_adjust),
      "3" = if (identical(petpeese, peese)) {
        list(model = peese, adjust_fun = PEESE_adjust)
      } else {
        list(model = fatpet, adjust_fun = PET_adjust)
      },
      stop("Invalid method")
    )

    do.call(
      compute_AR_CI_optimized,
      c(cfg, list(
        bs = bs, sebs = sebs, invNs = invNs, g = g, type_choice = type_choice
      ))
    )
  }

  ar_ci_res <- get_ar_ci_res()

  "p-value of test for publication bias / p-hacking based on instrumented FAT"
  model_list <- list(fatpet, peese, petpeese, ekreg)
  selected_model <- model_list[[method]]
  pb_p <- summary(selected_model)$coefficients[2, 4]
  slope_coef <- round(as.numeric(summary(selected_model)$coefficients[2, 1]), 3)

  get_se <- function(model, SE, dat, g, type_choice) {
    if (isTRUE(SE == 3)) {
      boot <- manual_wild_cluster_boot_se(
        model = model,
        data = dat,
        cluster_var = "g",
        B = 500
      )
      list(se = boot$boot_se[1], boot_result = boot)
    } else {
      v11 <- vcovCR(model, cluster = g, type = type_choice)[1, 1]
      list(se = sqrt(v11), boot_result = NULL)
    }
  }

  v11 <- function(model, g, type_choice) vcovCR(model, cluster = g, type = type_choice)[1, 1]

  # Map method -> (MAIVE model, Standard model, labels)
  cfg_map <- list(
    "1" = list(maive = fatpet, std = fatpet0, maive_label = "MAIVE-FAT-PET", std_label = "Standard FAT-PET"),
    "2" = list(maive = peese, std = peese0, maive_label = "MAIVE-PEESE", std_label = "Standard PEESE"),
    "3" = list(maive = petpeese, std = petpeese0, maive_label = "MAIVE-PET-PEESE", std_label = "Standard PET-PEESE"),
    "4" = list(maive = ekreg, std = ekreg0, maive_label = "MAIVE-EK", std_label = "Standard EK")
  )

  cfg <- cfg_map[[as.character(method)]]
  if (is.null(cfg)) stop("Invalid method")

  # MAIVE beta & SE (robust or WCR bootstrap)
  beta <- cfg$maive$coefficients[1]
  se_ma <- get_se(cfg$maive, SE, dat, g, type_choice)
  betase <- se_ma$se
  boot_result <- se_ma$boot_result

  # Standard beta & SE
  beta0 <- cfg$std$coefficients[1]
  se_st <- get_se(cfg$std, SE, dat, g, type_choice)
  beta0se <- se_st$se

  # Hausman-type test (conservative: denom = Var(beta_MAIVE))
  Hausman <- (cfg$maive$coefficients[1] - cfg$std$coefficients[1])^2 /
    v11(cfg$maive, g, type_choice)
  Chi2 <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)

  list(
    "beta"              = round(beta, 3),
    "SE"                = round(betase, 3),
    "F-test"            = F_hac,
    "beta_standard"     = round(beta0, 3),
    "SE_standard"       = round(beta0se, 3),
    "Hausman"           = round(Hausman, 3),
    "Chi2"              = round(Chi2, 3),
    "SE_instrumented"   = sebs2fit1^(1 / 2),
    "AR_CI"             = ar_ci_res$b0_CI,
    "pub bias p-value"  = round(pb_p, 3),
    "is_quadratic_fit"  = is_quadratic_fit,
    "boot_result"       = boot_result, # NULL unless SE == 3
    "slope_coef"        = slope_coef
  )
}

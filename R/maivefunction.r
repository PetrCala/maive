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
#' @keywords internal
maive_validate_inputs <- function(dat, method, weight, instrument, studylevel, SE, AR) {
  dat <- as.data.frame(dat)
  if (ncol(dat) < 3) {
    stop("dat must contain at least three columns: bs, sebs, and Ns.")
  }

  scalar_int <- function(value, name) {
    if (length(value) != 1L || is.na(value)) {
      stop(sprintf("%s must be a single non-missing value.", name))
    }
    as.integer(value)
  }

  method <- scalar_int(method, "method")
  weight <- scalar_int(weight, "weight")
  instrument <- scalar_int(instrument, "instrument")
  studylevel <- scalar_int(studylevel, "studylevel")
  SE <- scalar_int(SE, "SE")
  AR <- scalar_int(AR, "AR")

  if (!method %in% 1:4) stop("method must be between 1 and 4.")
  if (!weight %in% 0:2) stop("weight must be 0, 1, or 2.")
  if (!instrument %in% 0:1) stop("instrument must be 0 or 1.")
  if (!studylevel %in% 0:3) stop("studylevel must be between 0 and 3.")
  if (!SE %in% 0:3) stop("SE must be between 0 and 3.")
  if (!AR %in% 0:1) stop("AR must be 0 or 1.")

  type_map <- c("CR0", "CR1", "CR2")
  type_choice <- if (SE == 3L) "CR0" else type_map[SE + 1L]

  if (method == 4L || weight %in% c(1L, 2L) || instrument == 0L) {
    AR <- 0L
  }

  list(
    dat = dat,
    method = method,
    weight = weight,
    instrument = instrument,
    studylevel = studylevel,
    SE = SE,
    AR = AR,
    type_choice = type_choice,
    alpha_s = 0.05
  )
}

#' @keywords internal
maive_build_dummy_matrix <- function(values) {
  if (exists("to.dummy", mode = "function")) {
    to.dummy(data.frame(studyid = values), "studyid")
  } else {
    f <- factor(values)
    mm <- stats::model.matrix(~ f - 1)
    colnames(mm) <- paste0("studyid_", levels(f))
    mm
  }
}

#' @keywords internal
maive_center_dummy_matrix <- function(values) {
  D <- maive_build_dummy_matrix(values)
  if (is.null(dim(D)) || ncol(D) == 0L) {
    return(matrix(0, nrow = length(values), ncol = 0))
  }
  centered <- sweep(D, 2, colMeans(D), "-")
  if (ncol(centered) <= 1L) {
    centered[, 0, drop = FALSE]
  } else {
    centered[, seq_len(ncol(centered) - 1L), drop = FALSE]
  }
}

#' @keywords internal
maive_prepare_data <- function(dat, studylevel) {
  bs <- dat[[1]]
  sebs <- dat[[2]]
  Ns <- dat[[3]]

  if (!is.numeric(bs) || !is.numeric(sebs) || !is.numeric(Ns)) {
    stop("bs, sebs, and Ns must be numeric.")
  }
  if (any(!is.finite(bs)) || any(!is.finite(sebs)) || any(!is.finite(Ns))) {
    stop("bs, sebs, and Ns must be finite.")
  }
  if (any(sebs <= 0)) {
    stop("sebs must be strictly positive.")
  }
  if (any(Ns <= 0)) {
    stop("Ns must be strictly positive.")
  }

  M <- length(bs)
  if (length(sebs) != M || length(Ns) != M) {
    stop("bs, sebs, and Ns must have the same length.")
  }

  cluster <- studylevel %/% 2L
  dummy <- studylevel %% 2L

  if (ncol(dat) >= 4) {
    studyid <- dat[[4]]
  } else {
    studyid <- seq_len(M)
    dummy <- 0L
    cluster <- 0L
  }

  D <- maive_center_dummy_matrix(studyid)
  g <- if (cluster == 0L) seq_len(M) else studyid
  dat$g <- g

  list(
    dat = dat,
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    M = M,
    studyid = studyid,
    dummy = dummy,
    cluster = cluster,
    g = g,
    D = D
  )
}

#' @keywords internal
maive_compute_variance_instrumentation <- function(sebs, Ns, g, type_choice, instrument) {
  invNs <- 1 / Ns
  sebs2 <- sebs^2
  Xiv <- cbind(1, invNs)
  varreg1 <- lm(sebs2 ~ 0 + Xiv)
  dimiv <- 2L
  if (varreg1$coefficients[1] < 0) {
    Xiv <- as.matrix(invNs)
    varreg1 <- lm(sebs2 ~ 0 + Xiv)
    dimiv <- 1L
  }

  sebs2fit1 <- fitted(varreg1)
  if (instrument == 0L) {
    F_hac <- "NA"
  } else {
    V <- vcovCR(varreg1, cluster = g, type = type_choice)
    F_hac <- unname(round(varreg1$coefficients[dimiv]^2 / V[dimiv, dimiv], 3))
  }

  list(invNs = invNs, sebs2fit1 = sebs2fit1, F_hac = F_hac)
}

#' @keywords internal
maive_compute_weights <- function(weight, sebs, sebs2fit1) {
  if (weight == 0L) {
    rep(1, length(sebs))
  } else if (weight == 1L) {
    sebs
  } else if (weight == 2L) {
    sqrt(sebs2fit1)
  } else {
    stop("Invalid weight option.")
  }
}

#' @keywords internal
maive_build_design_matrices <- function(bs, sebs, w, x, x2, D, dummy) {
  y <- bs / w
  X <- cbind(1, x) / w
  X2 <- cbind(1, x2) / w

  y0 <- bs / sebs
  X0 <- cbind(1, sebs) / sebs
  X20 <- cbind(1, sebs^2) / sebs

  if (dummy == 1L) {
    X <- cbind(X, D / w)
    X2 <- cbind(X2, D / w)
    cD <- cbind(1, D) / w
    X0 <- cbind(X0, D / sebs)
    X20 <- cbind(X20, D / sebs)
    cD0 <- cbind(1, D) / sebs
  } else {
    cD <- matrix(1 / w, ncol = 1)
    cD0 <- matrix(1 / sebs, ncol = 1)
  }

  list(
    y = y,
    X = X,
    X2 = X2,
    y0 = y0,
    X0 = X0,
    X20 = X20,
    cD = cD,
    cD0 = cD0,
    w = w,
    sebs = sebs,
    x = x
  )
}

#' @keywords internal
maive_fit_models <- function(design) {
  y <- design$y
  X <- design$X
  X2 <- design$X2
  y0 <- design$y0
  X0 <- design$X0
  X20 <- design$X20
  cD <- design$cD
  cD0 <- design$cD0

  list(
    fatpet = lm(y ~ 0 + X),
    peese = lm(y ~ 0 + X2),
    fatpet0 = lm(y0 ~ 0 + X0),
    peese0 = lm(y0 ~ 0 + X20),
    wlsreg = lm(y ~ 0 + cD),
    wlsreg0 = lm(y0 ~ 0 + cD0)
  )
}

#' @keywords internal
maive_select_petpeese <- function(fits, design, alpha_s) {
  M <- length(design$y)

  quad_stat <- abs(coef(fits$fatpet)[1] / sqrt(vcov(fits$fatpet)[1, 1]))
  quad_cutoff <- qt(1 - alpha_s / 2, M - ncol(design$X))
  quadratic_decision <- quad_stat > quad_cutoff
  petpeese <- if (quadratic_decision) fits$peese else fits$fatpet

  quad_stat0 <- abs(coef(fits$fatpet0)[1] / sqrt(vcov(fits$fatpet0)[1, 1]))
  quad_cutoff0 <- qt(1 - alpha_s / 2, M - ncol(design$X0))
  quadratic_decision0 <- quad_stat0 > quad_cutoff0
  petpeese0 <- if (quadratic_decision0) fits$peese0 else fits$fatpet0

  list(
    quadratic_decision = quadratic_decision,
    petpeese = petpeese,
    quadratic_decision0 = quadratic_decision0,
    petpeese0 = petpeese0
  )
}

#' @keywords internal
maive_compute_sigma_h <- function(fits, w, sebs) {
  M <- length(w)
  wis0 <- 1 / (w^2)
  Qfe0 <- sum(residuals(fits$wlsreg)^2)
  denom0 <- M - ncol(model.matrix(fits$wlsreg)) - 1
  sigh2hat0 <- max(0, M * ((Qfe0 / denom0) - 1) / sum(wis0))
  sighhat0 <- sqrt(sigh2hat0)

  wis00 <- 1 / (sebs^2)
  Qfe00 <- sum(residuals(fits$wlsreg0)^2)
  denom00 <- M - ncol(model.matrix(fits$wlsreg0)) - 1
  sigh2hat00 <- max(0, M * ((Qfe00 / denom00) - 1) / sum(wis00))
  sighhat00 <- sqrt(sigh2hat00)

  list(sighhat0 = sighhat0, sighhat00 = sighhat00)
}

#' @keywords internal
maive_fit_ek <- function(selection, design, sighats, method) {
  if (method != 4L) {
    return(list(ekreg = NULL, ekreg0 = NULL, structure = "linear", a0 = NA_real_, a00 = NA_real_))
  }

  intercept <- coef(selection$petpeese)[1]
  threshold <- 1.96 * sighats$sighhat0
  if (intercept > threshold) {
    a0 <- (intercept - threshold) * (intercept + threshold) / (2 * 1.96 * intercept)
  } else {
    a0 <- 0
  }

  intercept0 <- coef(selection$petpeese0)[1]
  threshold0 <- 1.96 * sighats$sighhat00
  if (intercept0 > threshold0) {
    a00 <- (intercept0 - threshold0) * (intercept0 + threshold0) / (2 * 1.96 * intercept0)
  } else {
    a00 <- 0
  }

  y <- design$y
  cD <- design$cD
  y0 <- design$y0
  cD0 <- design$cD0

  if (!is.na(a0) && a0 > min(design$x) && a0 < max(design$x)) {
    xx_w <- (design$x - a0) * (design$x > a0) / design$w
    ekreg <- lm(y ~ 0 + cD + xx_w)
    ek_structure <- "kink"
  } else if (!is.na(a0) && a0 < min(design$x)) {
    x_w <- design$x / design$w
    ekreg <- lm(y ~ 0 + cD + x_w)
    ek_structure <- "linear"
  } else {
    ekreg <- lm(y ~ 0 + cD)
    ek_structure <- "intercept"
  }

  if (a00 > min(design$sebs) && a00 < max(design$sebs)) {
    xx0_w <- (design$sebs - a00) * (design$sebs > a00) / design$sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + xx0_w)
  } else if (a00 < min(design$sebs)) {
    x0_w <- design$sebs / design$sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + x0_w)
  } else {
    ekreg0 <- lm(y0 ~ 0 + cD0)
  }

  list(ekreg = ekreg, ekreg0 = ekreg0, structure = ek_structure, a0 = a0, a00 = a00)
}

#' @keywords internal
maive_infer_coef <- function(model, coef_index, SE, data, cluster_var, type_choice) {
  if (SE == 3L) {
    boot <- manual_wild_cluster_boot_se(
      model = model,
      data = data,
      cluster_var = cluster_var,
      B = 500
    )
    if (is.null(boot$boot_se)) {
      stop("Bootstrap helper must return boot_se.")
    }
    se <- unname(boot$boot_se[coef_index])
    boot_result <- boot
  } else {
    V <- vcovCR(model, cluster = data[[cluster_var]], type = type_choice)
    se <- unname(sqrt(V[coef_index, coef_index]))
    boot_result <- NULL
  }
  b <- unname(coef(model)[coef_index])
  t <- as.numeric(b / se)
  p <- 2 * pnorm(-abs(t))
  list(b = b, se = se, p = p, boot_result = boot_result)
}

#' @keywords internal
maive_get_intercept_se <- function(model, SE, data, cluster_var, type_choice) {
  inf <- maive_infer_coef(model, 1L, SE, data, cluster_var, type_choice)
  list(se = inf$se, boot_result = inf$boot_result)
}

#' @keywords internal
maive_slope_information <- function(method, fits, selection, ek) {
  method_str <- as.character(method)
  if (method_str == "1") {
    return(list(type = "linear", coefficient = round(as.numeric(coef(fits$fatpet)[2]), 3), detail = NULL))
  }
  if (method_str == "2") {
    return(list(type = "quadratic", coefficient = round(as.numeric(coef(fits$peese)[2]), 3), detail = NULL))
  }
  if (method_str == "3") {
    if (identical(selection$petpeese, fits$peese)) {
      return(list(type = "quadratic", coefficient = round(as.numeric(coef(fits$peese)[2]), 3), detail = NULL))
    }
    return(list(type = "linear", coefficient = round(as.numeric(coef(fits$fatpet)[2]), 3), detail = NULL))
  }
  if (method_str == "4") {
    if (is.null(ek$ekreg)) {
      return(list(type = "linear", coefficient = 0, detail = NULL))
    }
    if (ek$structure == "kink") {
      kink_effect <- round(as.numeric(tail(coef(ek$ekreg), 1)), 3)
      kink_location <- as.numeric(round(ek$a0, 3))
      detail <- list(kink_location = kink_location, kink_effect = kink_effect)
      return(list(
        type = "kinked",
        coefficient = list(kink_effect = kink_effect, kink_location = kink_location),
        detail = detail
      ))
    }
    if (ek$structure == "linear") {
      slope <- round(as.numeric(tail(coef(ek$ekreg), 1)), 3)
      return(list(type = "linear", coefficient = slope, detail = NULL))
    }
    return(list(type = "linear", coefficient = 0, detail = NULL))
  }
  stop("Invalid method")
}

#' @keywords internal
maive_quadratic_summary <- function(method, selection, slope_info) {
  quadratic_flag <- switch(as.character(method),
    "1" = FALSE,
    "2" = TRUE,
    "3" = isTRUE(selection$quadratic_decision),
    "4" = FALSE,
    stop("Invalid method")
  )
  list(
    quadratic = quadratic_flag,
    slope_type = slope_info$type,
    slope_detail = slope_info$detail
  )
}

#' @keywords internal
maive_get_config <- function(method, fits, selection, ek) {
  switch(as.character(method),
    "1" = list(maive = fits$fatpet, std = fits$fatpet0),
    "2" = list(maive = fits$peese, std = fits$peese0),
    "3" = list(maive = selection$petpeese, std = selection$petpeese0),
    "4" = list(maive = ek$ekreg, std = ek$ekreg0),
    stop("Invalid method")
  )
}

#' @keywords internal
maive_compute_hausman <- function(beta, beta0, model, g, type_choice) {
  V_maive <- vcovCR(model, cluster = g, type = type_choice)
  unname((beta - beta0)^2 / V_maive[1, 1])
}

#' @keywords internal
maive_compute_ar_ci <- function(opts, fits, selection, prepared, invNs, type_choice) {
  if (opts$AR != 1L || opts$method == 4L || opts$weight %in% c(1L, 2L)) {
    return(list(b0_CI = "NA", b1_CI = "NA"))
  }

  cfg_ar <- switch(as.character(opts$method),
    "1" = list(model = fits$fatpet, adjust_fun = PET_adjust),
    "2" = list(model = fits$peese, adjust_fun = PEESE_adjust),
    "3" = if (identical(selection$petpeese, fits$peese)) {
      list(model = fits$peese, adjust_fun = PEESE_adjust)
    } else {
      list(model = fits$fatpet, adjust_fun = PET_adjust)
    },
    stop("Invalid method")
  )

  do.call(
    compute_AR_CI_optimized,
    c(cfg_ar, list(
      bs = prepared$bs,
      sebs = prepared$sebs,
      invNs = invNs,
      g = prepared$g,
      type_choice = type_choice
    ))
  )
}

#' @export
maive <- function(dat, method, weight, instrument, studylevel, SE, AR) {
  opts <- maive_validate_inputs(dat, method, weight, instrument, studylevel, SE, AR)
  prepared <- maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- maive_compute_variance_instrumentation(prepared$sebs, prepared$Ns, prepared$g, opts$type_choice, opts$instrument)

  w <- maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1)
  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1

  design <- maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- maive_fit_models(design)
  selection <- maive_select_petpeese(fits, design, opts$alpha_s)
  sighats <- maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- maive_fit_ek(selection, design, sighats, opts$method)

  slope_info <- maive_slope_information(opts$method, fits, selection, ek)
  slope_summary <- maive_quadratic_summary(opts$method, selection, slope_info)

  egger_inf <- maive_infer_coef(fits$fatpet, 2L, opts$SE, prepared$dat, "g", opts$type_choice)
  cfg <- maive_get_config(opts$method, fits, selection, ek)
  if (is.null(cfg$maive) || is.null(cfg$std)) {
    stop("Failed to identify models for the selected method.")
  }

  beta <- unname(coef(cfg$maive)[1])
  beta0 <- unname(coef(cfg$std)[1])

  se_ma <- maive_get_intercept_se(cfg$maive, opts$SE, prepared$dat, "g", opts$type_choice)
  se_std <- maive_get_intercept_se(cfg$std, opts$SE, prepared$dat, "g", opts$type_choice)

  hausman <- maive_compute_hausman(beta, beta0, cfg$maive, prepared$g, opts$type_choice)
  chi2 <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)

  ar_ci_res <- maive_compute_ar_ci(opts, fits, selection, prepared, instrumentation$invNs, opts$type_choice)

  list(
    "beta" = round(beta, 3),
    "SE" = round(as.numeric(se_ma$se), 3),
    "F-test" = instrumentation$F_hac,
    "beta_standard" = round(beta0, 3),
    "SE_standard" = round(as.numeric(se_std$se), 3),
    "Hausman" = round(hausman, 3),
    "Chi2" = round(chi2, 3),
    "SE_instrumented" = sqrt(instrumentation$sebs2fit1),
    "AR_CI" = ar_ci_res$b0_CI,
    "pub bias p-value" = round(egger_inf$p, 3),
    "egger_coef" = round(egger_inf$b, 3),
    "egger_se" = round(egger_inf$se, 3),
    "is_quadratic_fit" = slope_summary,
    "boot_result" = se_ma$boot_result,
    "slope_coef" = slope_info$coefficient
  )
}

# Issue #24: the levels first stage can fit a negative variance for an
# individual estimate. Those estimates must be excluded explicitly from every
# quantity built on the instrumented SE, with the count reported.

negative_variance_fixture <- function() {
  i <- 0:19
  data.frame(
    bs = 0.3 + 0.08 * sin(i) + 0.03 * (i %% 3),
    sebs = 0.05 + 0.005 * i,
    Ns = 100 + 37 * i
  )
}

levels_fitted_variance <- function(dat) {
  fit <- lm(I(sebs^2) ~ I(1 / Ns), data = dat)
  unname(fitted(fit))
}

run_quietly <- function(expr) {
  withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
}

collect_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

run_maive <- function(dat, ...) {
  run_quietly(maive(dat, ...))
}

default_args <- list(method = 3, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 1, first_stage = 0)

test_that("the fixture reproduces a negative fitted variance in the levels first stage", {
  dat <- negative_variance_fixture()
  fitted_variance <- levels_fitted_variance(dat)

  expect_lt(fitted_variance[1], 0)
  expect_equal(sum(fitted_variance <= 0), 1L)
  expect_equal(fitted_variance[1], -0.00335, tolerance = 1e-2)
})

test_that("levels first stage excludes the estimate explicitly and reports it", {
  dat <- negative_variance_fixture()
  captured <- collect_warnings(do.call(maive, c(list(dat), default_args)))
  result <- captured$value

  exclusion_warnings <- grep("non-positive variance", captured$warnings, value = TRUE)
  expect_length(exclusion_warnings, 1L)
  expect_match(exclusion_warnings, "1 estimate")
  expect_match(exclusion_warnings, "first_stage = 1")
  expect_false(any(grepl("NaNs produced", captured$warnings)))

  expect_identical(result$n_excluded, 1L)
  expect_identical(result$excluded_rows, 1L)

  expect_length(result$SE_instrumented, nrow(dat))
  expect_true(is.na(result$SE_instrumented[1]))
  expect_false(is.nan(result$SE_instrumented[1]))
  expect_true(all(is.finite(result$SE_instrumented[-1])))
  expect_equal(result$SE_instrumented[-1], sqrt(levels_fitted_variance(dat)[-1]), tolerance = 1e-10)

  expect_length(result$weights, nrow(dat))
  expect_true(is.na(result$weights[1]))
  expect_true(all(is.finite(result$weights[-1])))
})

test_that("the count matches the number of rows whose square root would be NaN", {
  dat <- negative_variance_fixture()
  fitted_variance <- levels_fitted_variance(dat)
  nan_rows <- sum(is.nan(suppressWarnings(sqrt(fitted_variance))))

  result <- do.call(run_maive, c(list(dat), default_args))

  expect_identical(result$n_excluded, nan_rows)
  expect_identical(result$excluded_rows, which(fitted_variance <= 0))
})

test_that("PET, PEESE, the F-test and the Hausman comparison use the same rows", {
  dat <- negative_variance_fixture()
  opts <- MAIVE:::normalize_maive_options(dat, 3, 0, 1, 0, 0, 1, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs, prepared$Ns, prepared$g, opts$type_choice, opts$instrument, opts$first_stage_type
  )
  exclusion <- run_quietly(MAIVE:::maive_apply_variance_exclusion(opts, prepared, instrumentation))
  prepared_kept <- exclusion$prepared
  instrumentation_kept <- exclusion$instrumentation

  expect_identical(exclusion$n_excluded, 1L)
  expect_identical(prepared_kept$M, nrow(dat) - 1L)
  expect_true(all(instrumentation_kept$sebs2fit1 > 0))

  w <- MAIVE:::maive_compute_weights(opts$weight, prepared_kept$sebs, instrumentation_kept$sebs2fit1, prepared_kept$studyid)
  x <- sqrt(instrumentation_kept$sebs2fit1)
  x2 <- instrumentation_kept$sebs2fit1
  design <- MAIVE:::maive_build_design_matrices(prepared_kept$bs, prepared_kept$sebs, w, x, x2, prepared_kept$D, prepared_kept$dummy)
  fits <- MAIVE:::maive_fit_models(design)
  selection <- MAIVE:::maive_select_petpeese(fits, design, opts$alpha_s, opts$SE, prepared_kept$dat, opts$type_choice)
  cfg <- MAIVE:::maive_get_config(opts$method, fits, selection, NULL)
  hausman_cfg <- MAIVE:::maive_get_hausman_models(opts$method, cfg, selection, design)

  n_kept <- nrow(dat) - 1L
  expect_identical(stats::nobs(fits$fatpet), n_kept)
  expect_identical(stats::nobs(fits$peese), n_kept)
  expect_identical(stats::nobs(fits$fatpet0), n_kept)
  expect_identical(stats::nobs(fits$peese0), n_kept)
  expect_identical(stats::nobs(hausman_cfg$maive), n_kept)
  expect_identical(stats::nobs(hausman_cfg$std), n_kept)
  expect_false(anyNA(design$x))

  # The reported F-test is the first stage refit on the same rows.
  first_stage_kept <- lm(I(sebs^2) ~ I(1 / Ns), data = dat[-1, ])
  V <- clubSandwich::vcovCR(first_stage_kept, cluster = seq_len(n_kept), type = "CR0")
  expected_f <- unname(coef(first_stage_kept)[2]^2 / V[2, 2])

  result <- do.call(run_maive, c(list(dat), default_args))
  expect_equal(result$`F-test`, expected_f, tolerance = 1e-10)

  # The Hausman statistic compares two fits on the same rows.
  V_iv <- clubSandwich::vcovCR(hausman_cfg$maive, cluster = prepared_kept$g, type = opts$type_choice)
  V_ols <- clubSandwich::vcovCR(hausman_cfg$std, cluster = prepared_kept$g, type = opts$type_choice)
  expected_hausman <- as.numeric((coef(hausman_cfg$maive)[1] - coef(hausman_cfg$std)[1])^2 / (V_iv[1, 1] - V_ols[1, 1]))
  expect_equal(as.numeric(result$Hausman), expected_hausman, tolerance = 1e-10)
  expect_true(is.finite(result$Hausman))
})

test_that("the estimate is the second stage on the kept rows with the full-sample fitted variances", {
  dat <- negative_variance_fixture()
  result <- do.call(run_maive, c(list(dat), default_args))

  kept <- dat[-1, ]
  fitted_kept <- levels_fitted_variance(dat)[-1]
  expected_beta <- if (identical(result$petpeese_selected, "PEESE")) {
    unname(coef(lm(kept$bs ~ fitted_kept))[1])
  } else {
    unname(coef(lm(kept$bs ~ sqrt(fitted_kept)))[1])
  }
  expect_equal(result$beta, expected_beta, tolerance = 1e-10)
})

test_that("an excluded estimate cannot move the estimate or silently disable the Hausman test", {
  dat <- negative_variance_fixture()
  result <- do.call(run_maive, c(list(dat), default_args))

  outlier <- dat
  outlier$bs[1] <- 100
  result_outlier <- do.call(run_maive, c(list(outlier), default_args))

  expect_identical(result_outlier$n_excluded, 1L)
  expect_identical(result_outlier$beta, result$beta)
  expect_identical(result_outlier$SE, result$SE)
  expect_identical(result_outlier$`F-test`, result$`F-test`)
  expect_true(is.finite(result_outlier$Hausman))
  expect_identical(result_outlier$Hausman, result$Hausman)
})

test_that("AR = 1 returns a finite interval after the exclusion", {
  dat <- negative_variance_fixture()
  result <- do.call(run_maive, c(list(dat), default_args))

  expect_true(is.numeric(result$AR_CI))
  expect_length(result$AR_CI, 2L)
  expect_true(all(is.finite(result$AR_CI)))
  expect_lt(result$AR_CI[1], result$beta)
  expect_gt(result$AR_CI[2], result$beta)
})

test_that("adjusted weights with AR do not abort", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$weight <- 2
  result <- expect_no_error(do.call(run_maive, c(list(dat), args)))

  expect_identical(result$n_excluded, 1L)
  expect_true(is.finite(result$beta))
  expect_true(is.finite(result$SE))
  expect_true(all(is.finite(result$AR_CI)))
  expect_true(is.na(result$weights[1]))
  expect_true(all(result$weights[-1] > 0))
})

test_that("EK does not abort", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$method <- 4
  result <- expect_no_error(do.call(run_maive, c(list(dat), args)))

  expect_identical(result$n_excluded, 1L)
  expect_true(is.finite(result$beta))
  expect_true(is.finite(result$SE))
  expect_true(result$ek_structure %in% c("kink", "linear", "intercept"))
})

test_that("the wild bootstrap runs on the kept rows", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$SE <- 3
  args$AR <- 0
  result <- expect_no_error(do.call(run_maive, c(list(dat), args)))

  expect_identical(result$n_excluded, 1L)
  expect_true(is.finite(result$beta))
  expect_true(is.finite(result$SE))
  expect_true(all(is.finite(result$egger_boot_ci)))
})

test_that("WAIVE excludes the same estimate", {
  dat <- negative_variance_fixture()
  result <- run_quietly(do.call(waive, c(list(dat), default_args)))

  expect_identical(result$n_excluded, 1L)
  expect_identical(result$excluded_rows, 1L)
  expect_true(is.na(result$weights[1]))
  expect_true(all(is.finite(result$weights[-1])))
  expect_true(is.finite(result$beta))
})

test_that("a study that loses its only estimate does not break fixed effects or clustering", {
  dat <- negative_variance_fixture()
  dat$study_id <- c(99, rep(1:4, length.out = nrow(dat) - 1L))
  args <- default_args
  args$studylevel <- 3
  args$SE <- 1
  result <- expect_no_error(do.call(run_maive, c(list(dat), args)))

  expect_identical(result$n_excluded, 1L)
  expect_true(is.finite(result$beta))
  expect_true(is.finite(result$SE))
  expect_true(is.finite(result$Hausman))
})

test_that("log first stage excludes nothing and does not warn", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$first_stage <- 1
  captured <- collect_warnings(do.call(maive, c(list(dat), args)))
  result <- captured$value

  expect_false(any(grepl("non-positive variance", captured$warnings)))
  expect_identical(result$n_excluded, 0L)
  expect_identical(result$excluded_rows, integer(0))
  expect_length(result$SE_instrumented, nrow(dat))
  expect_true(all(is.finite(result$SE_instrumented)))
  expect_true(all(is.finite(result$weights)))
})

test_that("a clean dataset excludes nothing", {
  dat <- read.csv(test_path("fixtures", "euro.csv"))
  captured <- collect_warnings(do.call(maive, c(list(dat), default_args)))
  result <- captured$value

  expect_false(any(grepl("non-positive variance", captured$warnings)))
  expect_identical(result$n_excluded, 0L)
  expect_identical(result$excluded_rows, integer(0))
  expect_true(all(is.finite(result$SE_instrumented)))
})

test_that("without instrumenting nothing is excluded and SE_instrumented carries NA, not NaN", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$instrument <- 0
  args$AR <- 0
  captured <- collect_warnings(do.call(maive, c(list(dat), args)))
  result <- captured$value

  expect_false(any(grepl("NaNs produced", captured$warnings)))
  expect_false(any(grepl("non-positive variance", captured$warnings)))
  expect_identical(result$n_excluded, 0L)
  expect_true(is.na(result$SE_instrumented[1]))
  expect_false(is.nan(result$SE_instrumented[1]))
  expect_true(all(is.finite(result$weights)))
})

test_that("adjusted weights without instrumenting still exclude the unusable estimate", {
  dat <- negative_variance_fixture()
  args <- default_args
  args$instrument <- 0
  args$weight <- 2
  args$AR <- 0
  result <- do.call(run_maive, c(list(dat), args))

  expect_identical(result$n_excluded, 1L)
  expect_true(is.na(result$weights[1]))
  expect_true(is.finite(result$beta))
})

test_that("too few usable estimates abort with a pointer to the log first stage", {
  dat <- data.frame(
    bs = seq(0.1, 1.1, length.out = 11),
    sebs = c(rep(0.01, 9), 0.5, 0.5),
    Ns = c(rep(10, 8), 20, 1000, 2000)
  )
  fitted_variance <- levels_fitted_variance(dat)
  expect_gt(sum(fitted_variance <= 0), nrow(dat) - 4L)

  expect_error(
    run_quietly(do.call(maive, c(list(dat), default_args))),
    "first_stage = 1"
  )
})

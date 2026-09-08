test_that("log first stage applies smearing retransformation", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 1
  )

  log_model <- lm(log(dat$sebs^2) ~ log(dat$Ns))
  smearing <- mean(exp(residuals(log_model)))
  sehat_manual <- exp(predict(log_model)) * smearing

  expect_equal(result$SE_instrumented, unname(sqrt(sehat_manual)), tolerance = 1e-10)

  manual_vcov <- clubSandwich::vcovCR(log_model, cluster = seq_len(nrow(dat)), type = "CR0")
  slope <- coef(log_model)[2]
  manual_F <- unname(slope^2 / manual_vcov[2, 2])
  expect_equal(result$`F-test`, manual_F)
})

test_that("first-stage F-test does not crash when Ns is constant (rank deficient)", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = rep(100, 6)
  )

  # With constant Ns, the instrument (1/Ns) has no variation so IV is not identified.
  # MAIVE should not error; it should fall back to instrument=0 and report F-test as "NA".
  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  expect_identical(result$`F-test`, "NA")
  expect_true(is.finite(as.numeric(result$beta)))
  expect_true(is.finite(as.numeric(result$SE)))
})

test_that("Hausman statistic uses difference-in-estimators variance", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 2,
    SE = 0,
    AR = 1,
    first_stage = 0
  )

  opts <- MAIVE:::normalize_maive_options(dat, 1, 0, 1, 2, 0, 1, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs,
    prepared$Ns,
    prepared$g,
    opts$type_choice,
    opts$instrument,
    opts$first_stage_type
  )
  w <- MAIVE:::maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)
  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1
  design <- MAIVE:::maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- MAIVE:::maive_fit_models(design)
  selection <- MAIVE:::maive_select_petpeese(fits, design, opts$alpha_s)
  sighats <- MAIVE:::maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- MAIVE:::maive_fit_ek(selection, design, sighats, opts$method)
  cfg <- MAIVE:::maive_get_config(opts$method, fits, selection, ek)

  V_iv <- clubSandwich::vcovCR(cfg$maive, cluster = prepared$g, type = opts$type_choice)
  V_ols <- clubSandwich::vcovCR(cfg$std, cluster = prepared$g, type = opts$type_choice)
  var_diff <- V_iv[1, 1] - V_ols[1, 1]
  expected <- (coef(cfg$maive)[1] - coef(cfg$std)[1])^2 / var_diff
  expected_value <- as.numeric(expected)

  expect_equal(as.numeric(result$Hausman), expected_value)
})

test_that("Hausman PET-PEESE uses MAIVE weights", {
  dat <- data.frame(
    bs = c(0.42, 0.5, 0.48, 0.55, 0.46, 0.53),
    sebs = c(0.18, 0.2, 0.19, 0.21, 0.22, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92)
  )

  result <- maive(
    dat = dat,
    method = 3,
    weight = 2,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  opts <- MAIVE:::normalize_maive_options(dat, 3, 2, 1, 0, 0, 0, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs,
    prepared$Ns,
    prepared$g,
    opts$type_choice,
    opts$instrument,
    opts$first_stage_type
  )
  w <- MAIVE:::maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)
  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1
  design <- MAIVE:::maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- MAIVE:::maive_fit_models(design)
  selection <- MAIVE:::maive_select_petpeese(fits, design, opts$alpha_s)
  sighats <- MAIVE:::maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- MAIVE:::maive_fit_ek(selection, design, sighats, opts$method)
  cfg <- MAIVE:::maive_get_config(opts$method, fits, selection, ek)
  hausman_cfg <- MAIVE:::maive_get_hausman_models(opts$method, cfg, selection, design)

  V_iv <- clubSandwich::vcovCR(hausman_cfg$maive, cluster = prepared$g, type = opts$type_choice)
  V_ols <- clubSandwich::vcovCR(hausman_cfg$std, cluster = prepared$g, type = opts$type_choice)
  var_diff <- V_iv[1, 1] - V_ols[1, 1]
  expected <- (coef(hausman_cfg$maive)[1] - coef(hausman_cfg$std)[1])^2 / var_diff
  expected_value <- as.numeric(expected)

  expect_equal(as.numeric(result$Hausman), expected_value)
})

first_stage_default_fixture <- function() {
  Ns <- c(50, 80, 120, 200, 300, 450, 600, 800, 1000, 1500, 2000, 3000)
  scale <- c(1.1, 0.9, 1.2, 0.8, 1.0, 1.3, 0.7, 1.05, 0.95, 1.15, 0.85, 1.0)
  sebs <- scale / sqrt(Ns)
  noise <- c(0.02, -0.01, 0.03, -0.02, 0.01, 0.00, -0.03, 0.02, -0.01, 0.01, 0.00, -0.02)
  data.frame(
    bs = 0.3 + 0.5 * sebs + noise,
    sebs = sebs,
    Ns = Ns,
    study_id = rep(1:6, each = 2)
  )
}

first_stage_default_fields <- c("beta", "SE", "F-test", "SE_instrumented", "Hausman", "egger_coef")

test_that("normalize_maive_options() defaults first_stage to log", {
  dat <- first_stage_default_fixture()
  base <- list(dat = dat, method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0)

  omitted <- do.call(MAIVE:::normalize_maive_options, base)
  expect_identical(omitted$first_stage, 1L)
  expect_identical(omitted$first_stage_type, "log")

  null_arg <- do.call(MAIVE:::normalize_maive_options, c(base, list(first_stage = NULL)))
  expect_identical(null_arg$first_stage, 1L)
  expect_identical(null_arg$first_stage_type, "log")

  levels_arg <- do.call(MAIVE:::normalize_maive_options, c(base, list(first_stage = 0)))
  expect_identical(levels_arg$first_stage, 0L)
  expect_identical(levels_arg$first_stage_type, "levels")

  levels_name <- do.call(MAIVE:::normalize_maive_options, c(base, list(first_stage = "levels")))
  expect_identical(levels_name$first_stage, 0L)
})

test_that("maive() uses the log first stage when first_stage is omitted", {
  dat <- first_stage_default_fixture()
  args <- list(dat = dat, method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0)

  default <- do.call(maive, args)
  log_stage <- do.call(maive, c(args, list(first_stage = 1)))
  levels_stage <- do.call(maive, c(args, list(first_stage = 0)))

  expect_equal(default[first_stage_default_fields], log_stage[first_stage_default_fields])
  expect_false(isTRUE(all.equal(default$SE_instrumented, levels_stage$SE_instrumented)))
  expect_false(isTRUE(all.equal(default$beta, levels_stage$beta)))
  expect_false(isTRUE(all.equal(default$`F-test`, levels_stage$`F-test`)))

  # Explicit NULL falls back to the same default
  null_stage <- do.call(maive, c(args, list(first_stage = NULL)))
  expect_equal(null_stage[first_stage_default_fields], log_stage[first_stage_default_fields])
})

test_that("waive() uses the log first stage when first_stage is omitted", {
  dat <- first_stage_default_fixture()
  args <- list(dat = dat, method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0)

  default <- do.call(waive, args)
  log_stage <- do.call(waive, c(args, list(first_stage = 1)))
  levels_stage <- do.call(waive, c(args, list(first_stage = 0)))

  expect_equal(default[first_stage_default_fields], log_stage[first_stage_default_fields])
  expect_equal(default$weights, log_stage$weights)
  expect_false(isTRUE(all.equal(default$SE_instrumented, levels_stage$SE_instrumented)))
  expect_false(isTRUE(all.equal(default$beta, levels_stage$beta)))
  expect_false(isTRUE(all.equal(default$weights, levels_stage$weights)))
})

test_that("first_stage = 0 reproduces the levels first stage from the paper", {
  dat <- first_stage_default_fixture()
  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  # The published specification regresses sebs^2 on a constant and 1/Ns and
  # refits through the origin when the intercept is negative.
  levels_model <- lm(I(dat$sebs^2) ~ I(1 / dat$Ns))
  if (coef(levels_model)[1] < 0) {
    levels_model <- lm(I(dat$sebs^2) ~ 0 + I(1 / dat$Ns))
  }
  expect_equal(result$SE_instrumented, unname(sqrt(fitted(levels_model))), tolerance = 1e-10)
})

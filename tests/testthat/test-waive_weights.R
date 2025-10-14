test_that("WAIVE weight construction is stable and normalized", {
  df <- data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1.1, 2.0, 2.9, 4.1, 4.8)
  )
  fit <- lm(y ~ x, data = df)
  weights <- MAIVE:::waive_compute_decay_weights(fit)

  expect_length(weights, nrow(df))
  expect_gt(min(weights), 0)
  expect_equal(mean(weights), 1, tolerance = 1e-12)
})

test_that("WAIVE matches MAIVE when first stage has zero residuals", {
  Ns <- c(50, 75, 100, 125, 150)
  se2 <- 0.1 + 2 / Ns
  dat <- data.frame(
    bs = c(0.2, 0.25, 0.22, 0.3, 0.28),
    sebs = sqrt(se2),
    Ns = Ns
  )

  maive_res <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )
  waive_res <- waive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  expect_equal(waive_res$beta, maive_res$beta)
  expect_equal(waive_res$SE, maive_res$SE)
})

test_that("WAIVE downweights suspiciously small residuals", {
  dat <- data.frame(
    bs = c(0.3, 0.45, 0.5, 0.6, 0.55),
    sebs = c(0.22, 0.19, 0.18, 0.24, 0.23),
    Ns = c(60, 70, 80, 90, 100)
  )

  opts <- MAIVE:::maive_validate_inputs(dat, 1, 0, 1, 0, 0, 0, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs,
    prepared$Ns,
    prepared$g,
    opts$type_choice,
    opts$instrument,
    opts$first_stage_type
  )
  decay <- MAIVE:::waive_compute_decay_weights(instrumentation$first_stage_model)
  base_w <- MAIVE:::maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1)
  waive_w <- base_w * sqrt(decay)

  down_idx <- which(decay < 1)
  expect_true(length(down_idx) > 0)
  expect_true(all(waive_w[down_idx] < base_w[down_idx]))
})

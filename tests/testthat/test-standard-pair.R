education_fixture <- function() {
  read.csv(test_path("fixtures", "education.csv"), stringsAsFactors = FALSE)
}

run_education <- function(method, weight = 0) {
  suppressWarnings(maive(
    education_fixture(),
    method = method, weight = weight, instrument = 1, studylevel = 2, SE = 2, AR = 0, first_stage = 0
  ))
}

conventional_pet <- function(dat) {
  y0 <- dat$bs / dat$sebs
  X0 <- cbind(1, dat$sebs) / dat$sebs
  model <- lm(y0 ~ 0 + X0)
  V <- clubSandwich::vcovCR(model, cluster = dat$study_id, type = "CR2")
  list(coef = unname(coef(model)[1]), se = unname(sqrt(V[1, 1])))
}

test_that("beta_standard and SE_standard at method 3 come from the same conventional fit", {
  skip_if_not_installed("clubSandwich")
  dat <- education_fixture()
  res <- run_education(method = 3)

  # The conventional PET-PEESE selection on this data is PET (t = 1.17 on the
  # intercept), so the pair must equal the conventional PET coefficient and SE.
  pet <- conventional_pet(dat)
  expect_equal(res$beta_standard, pet$coef)
  expect_equal(res$SE_standard, pet$se)

  # Same pair as methods 1 and 4, which already used the conventional fit
  res_pet <- run_education(method = 1)
  expect_equal(res$beta_standard, res_pet$beta_standard)
  expect_equal(res$SE_standard, res_pet$SE_standard)

  # Values on 0.2.5 for reference: beta_standard was -0.0594701613 at method 3
  # (from the auxiliary weighted model) while SE_standard was 0.0175153250.
  expect_equal(res$beta_standard, 0.020527711568525672, tolerance = 1e-10)
  expect_equal(res$SE_standard, 0.017515324984599116, tolerance = 1e-10)
})

test_that("the conventional pair at method 3 does not move with MAIVE weights", {
  skip_if_not_installed("clubSandwich")
  results <- lapply(0:3, function(weight) run_education(method = 3, weight = weight))

  betas <- vapply(results, function(r) r$beta_standard, numeric(1))
  ses <- vapply(results, function(r) r$SE_standard, numeric(1))

  expect_equal(betas, rep(betas[1], 4))
  expect_equal(ses, rep(ses[1], 4))

  # The MAIVE estimate itself still responds to the weighting scheme
  maive_betas <- vapply(results, function(r) r$beta, numeric(1))
  expect_true(length(unique(round(maive_betas, 6))) > 1)
})

test_that("the Hausman statistic at method 3 is unchanged by the beta_standard fix", {
  skip_if_not_installed("clubSandwich")
  # Reference values computed with MAIVE 0.2.5 on the same fixture and options
  reference <- c(NA_real_, NA_real_, 7.8151132630996862, NA_real_)

  for (weight in 0:3) {
    res <- run_education(method = 3, weight = weight)
    expected <- reference[weight + 1]
    if (is.na(expected)) {
      expect_true(is.na(res$Hausman), label = sprintf("weight %d Hausman is NA", weight))
    } else {
      expect_equal(res$Hausman, expected, tolerance = 1e-10)
    }
  }
})

test_that("methods 1, 2, and 4 keep their conventional pair", {
  skip_if_not_installed("clubSandwich")
  reference <- list(
    "1" = c(beta = 0.020527711568525672, se = 0.017515324984599116, hausman = 0.6145003251123321),
    "2" = c(beta = -0.019854274240421819, se = 0.026056698213352728, hausman = 0.4570732166050237),
    "4" = c(beta = 0.020527711568525672, se = 0.017515324984599116, hausman = 0.6145003251123321)
  )
  for (method in c(1, 2, 4)) {
    res <- run_education(method = method)
    ref <- reference[[as.character(method)]]
    expect_equal(res$beta_standard, unname(ref["beta"]), tolerance = 1e-10)
    expect_equal(res$SE_standard, unname(ref["se"]), tolerance = 1e-10)
    expect_equal(res$Hausman, unname(ref["hausman"]), tolerance = 1e-10)
  }
})

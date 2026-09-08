mapping_fixture <- function() {
  data.frame(
    my_est = c(0.50, 0.60, 0.40, 0.55, 0.45, 0.52, 0.48, 0.58, 0.42, 0.50),
    my_se = c(0.20, 0.18, 0.25, 0.22, 0.24, 0.19, 0.23, 0.17, 0.26, 0.21),
    my_n = c(80, 120, 95, 110, 90, 130, 100, 140, 85, 105),
    my_study = c("A", "A", "B", "B", "C", "C", "D", "D", "A", "B"),
    stringsAsFactors = FALSE
  )
}

test_that("maive() accepts custom column names via estimate/se/n/study_id", {
  custom <- mapping_fixture()

  res_custom <- suppressWarnings(maive(
    dat = custom,
    estimate = "my_est", se = "my_se", n = "my_n", study_id = "my_study",
    method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 2, AR = 0, first_stage = 0
  ))

  standard <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n, study_id = custom$my_study,
    stringsAsFactors = FALSE
  )
  res_standard <- suppressWarnings(maive(
    standard,
    method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 2, AR = 0, first_stage = 0
  ))

  expect_equal(res_custom$beta, res_standard$beta)
  expect_equal(res_custom$SE, res_standard$SE)
  expect_equal(res_custom$Hausman, res_standard$Hausman)
  expect_equal(res_custom$weights, res_standard$weights)
})

test_that("custom column names work without a study identifier", {
  custom <- mapping_fixture()[, c("my_est", "my_se", "my_n")]

  res <- suppressWarnings(maive(
    dat = custom,
    estimate = "my_est", se = "my_se", n = "my_n",
    method = 1, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 0, first_stage = 0
  ))

  expect_true(is.numeric(res$beta))
  expect_true(is.finite(res$beta))
})

test_that("completely empty rows are dropped before custom columns are resolved", {
  custom <- mapping_fixture()
  with_empty <- rbind(custom, data.frame(my_est = NA, my_se = NA, my_n = NA, my_study = NA))

  expect_message(
    res_empty <- suppressWarnings(maive(
      dat = with_empty,
      estimate = "my_est", se = "my_se", n = "my_n", study_id = "my_study",
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    )),
    "Removed 1 completely empty row"
  )
  res_clean <- suppressWarnings(maive(
    dat = custom,
    estimate = "my_est", se = "my_se", n = "my_n", study_id = "my_study",
    method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
  ))

  expect_equal(res_empty$beta, res_clean$beta)
})

test_that("a missing custom column is reported by its mapped name", {
  custom <- mapping_fixture()
  expect_error(
    maive(
      dat = custom,
      estimate = "not_there", se = "my_se", n = "my_n",
      method = 1, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 0, first_stage = 0
    ),
    "Missing required columns: not_there"
  )
})

test_that("positional study_id fallback warns and names the column", {
  custom <- mapping_fixture()
  four_col <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n,
    year = rep(c(2001, 2002), 5)
  )

  expect_warning(
    maive(four_col, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0),
    "using the fourth column \\('year'\\)"
  )
})

test_that("an explicit study_id suppresses the positional fallback warning", {
  custom <- mapping_fixture()
  four_col <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n,
    year = rep(c(2001, 2002), 5)
  )

  expect_no_warning(
    maive(
      four_col,
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0,
      study_id = "year"
    )
  )
})

test_that("a column named study_id is used silently wherever it sits", {
  custom <- mapping_fixture()
  five_col <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n,
    year = rep(c(2001, 2002), 5),
    study_id = custom$my_study,
    stringsAsFactors = FALSE
  )
  four_col <- five_col[, c("bs", "sebs", "Ns", "study_id")]

  expect_no_warning(
    res_five <- maive(five_col, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
  )
  res_four <- maive(four_col, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)

  expect_equal(res_five$beta, res_four$beta)
  expect_equal(res_five$SE, res_four$SE)
})

test_that("the positional fallback still drives clustering when accepted", {
  custom <- mapping_fixture()
  four_col <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n,
    grp = custom$my_study, stringsAsFactors = FALSE
  )
  named <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n,
    study_id = custom$my_study, stringsAsFactors = FALSE
  )

  res_fallback <- suppressWarnings(
    maive(four_col, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
  )
  res_named <- maive(named, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)

  expect_equal(res_fallback$SE, res_named$SE)
})

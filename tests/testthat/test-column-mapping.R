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

test_that("the positional fallback never picks a column mapped as estimate, se, or n", {
  custom <- mapping_fixture()
  # Identifier first, sample size fourth: the layout harmonised data usually has
  lit <- data.frame(
    study_label = custom$my_study, effect = custom$my_est, se = custom$my_se, n_obs = custom$my_n,
    stringsAsFactors = FALSE
  )
  named <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n, study_id = custom$my_study,
    stringsAsFactors = FALSE
  )

  expect_warning(
    res <- maive(
      lit,
      estimate = "effect", se = "se", n = "n_obs",
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    ),
    "using the only unmapped column \\('study_label'\\)"
  )
  res_named <- maive(named, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
  expect_equal(res$beta, res_named$beta)
  expect_equal(res$SE, res_named$SE)
})

test_that("the fourth column keeps priority when it is free", {
  custom <- mapping_fixture()
  wide <- data.frame(
    year = rep(c(2001, 2002), 5), effect = custom$my_est, se = custom$my_se,
    paper_id = custom$my_study, n_obs = custom$my_n,
    stringsAsFactors = FALSE
  )
  named <- data.frame(
    bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n, study_id = custom$my_study,
    stringsAsFactors = FALSE
  )

  expect_warning(
    res <- maive(
      wide,
      estimate = "effect", se = "se", n = "n_obs",
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    ),
    "using the fourth column \\('paper_id'\\)"
  )
  res_named <- maive(named, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
  expect_equal(res$SE, res_named$SE)
})

test_that("an ambiguous fallback asks for study_id instead of guessing", {
  custom <- mapping_fixture()
  wide <- data.frame(
    study_label = custom$my_study, effect = custom$my_est, se = custom$my_se, n_obs = custom$my_n,
    year = rep(c(2001, 2002), 5),
    stringsAsFactors = FALSE
  )

  expect_error(
    maive(
      wide,
      estimate = "effect", se = "se", n = "n_obs",
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    ),
    "already mapped.*study_label.*year|Several columns.*study_label"
  )
  expect_error(
    MAIVE:::resolve_maive_columns(wide, estimate = "effect", se = "se", n = "n_obs", studylevel = 2),
    "Pass study_id"
  )

  # Naming the identifier resolves it
  expect_no_warning(
    maive(
      wide,
      estimate = "effect", se = "se", n = "n_obs", study_id = "study_label",
      method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    )
  )

  # At studylevel = 0 no identifier is needed, so the ambiguity is skipped silently
  expect_no_warning(
    res <- maive(
      wide,
      estimate = "effect", se = "se", n = "n_obs",
      method = 1, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 0, first_stage = 0
    )
  )
  expect_true(is.finite(res$beta))
})

test_that("a positionally inferred id does not trigger the degrees-of-freedom rule at studylevel = 0", {
  # Six rows, one per study, with a year in column four: ran in 0.2.5, aborted since 0.2.6
  six <- data.frame(
    bs = c(0.50, 0.60, 0.40, 0.55, 0.45, 0.52),
    sebs = c(0.20, 0.18, 0.25, 0.22, 0.24, 0.19),
    Ns = c(80, 120, 95, 110, 90, 130),
    year = c(2001, 2002, 2003, 2004, 2005, 2006)
  )

  expect_warning(
    res <- maive(six, method = 1, weight = 0, instrument = 0, studylevel = 0, SE = 0, AR = 0, first_stage = 0),
    "using the fourth column \\('year'\\)"
  )
  expect_true(is.finite(res$beta))

  # Clustering alone does not fit one regressor per study either
  expect_warning(
    res_cluster <- maive(six, method = 1, weight = 0, instrument = 0, studylevel = 2, SE = 0, AR = 0, first_stage = 0),
    "using the fourth column \\('year'\\)"
  )
  expect_true(is.finite(res_cluster$beta))

  # Study dummies still need the rows
  expect_error(
    suppressWarnings(maive(six, method = 1, weight = 0, instrument = 0, studylevel = 1, SE = 0, AR = 0, first_stage = 0)),
    "Insufficient degrees of freedom: 6 observations with 6 unique studies[[:space:]]+requires at least 9 rows"
  )
})

test_that("studylevel > 0 without any study identifier warns that it has no effect", {
  custom <- mapping_fixture()
  three_col <- data.frame(bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n)

  expect_no_warning(
    res0 <- maive(three_col, method = 1, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 0, first_stage = 0)
  )
  for (level in 1:3) {
    expect_warning(
      res <- maive(three_col, method = 1, weight = 0, instrument = 1, studylevel = level, SE = 0, AR = 0, first_stage = 0),
      sprintf("studylevel = %d has no effect because the data has no study identifier", level)
    )
    expect_equal(res$beta, res0$beta)
  }

  expect_warning(
    waive(three_col, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0),
    "studylevel = 2 has no effect"
  )
})

test_that("a three-column frame with one column mapped twice reaches the fallback message", {
  custom <- mapping_fixture()
  three_col <- data.frame(bs = custom$my_est, sebs = custom$my_se, Ns = custom$my_n)

  # se and n both point at sebs, so Ns is the only unmapped column; the frame
  # has no fourth column to compare against
  expect_warning(
    maive(
      three_col,
      n = "sebs",
      method = 1, weight = 0, instrument = 0, studylevel = 2, SE = 0, AR = 0, first_stage = 0
    ),
    "using the only unmapped column \\('Ns'\\)"
  )
})

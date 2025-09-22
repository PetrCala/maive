strip_names <- function(x) {
  if (is.list(x)) {
    lapply(x, strip_names)
  } else {
    if (!is.null(names(x))) {
      names(x) <- NULL
    }
    x
  }
}

test_that("maive results match baseline fixtures", {
  baseline_path <- testthat::test_path("..", "fixtures", "baseline_results.rds")
  skip_if_not(file.exists(baseline_path), "Baseline results fixture not found")

  baseline <- readRDS(baseline_path)
  expect_true(is.list(baseline))
  expect_gt(length(baseline), 0L)

  dat1 <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  dat2 <- data.frame(
    bs = c(0.2, 0.25, 0.3, 0.28, 0.26, 0.27, 0.29),
    sebs = c(0.15, 0.17, 0.16, 0.18, 0.17, 0.16, 0.15),
    Ns = c(60, 70, 65, 75, 80, 78, 72)
  )

  scenarios <- list(
    s1 = list(data = dat1, args = list(method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 1)),
    s2 = list(data = dat1, args = list(method = 2, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 1)),
    s3 = list(data = dat1, args = list(method = 3, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 1)),
    s4 = list(data = dat1, args = list(method = 4, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 1)),
    s5 = list(data = dat2, args = list(method = 1, weight = 0, instrument = 0, studylevel = 0, SE = 0, AR = 0)),
    s6 = list(data = dat1, args = list(method = 1, weight = 1, instrument = 1, studylevel = 2, SE = 1, AR = 1)),
    s7 = list(data = dat1, args = list(method = 1, weight = 2, instrument = 1, studylevel = 1, SE = 2, AR = 1))
  )

  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    expect_true(scenario_name %in% names(baseline))

    actual <- do.call(maive, c(list(dat = scenario$data), scenario$args))
    expected <- baseline[[scenario_name]]

    expect_equal(
      strip_names(actual),
      strip_names(expected),
      tolerance = 1e-8
    )
  }
})

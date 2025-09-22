#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clubSandwich)
  library(sandwich)
})

source("R/boot.r")
source("R/ar.r")
source("R/maivefunction.r")

baseline_path <- "tests/fixtures/baseline_results.rds"
if (!file.exists(baseline_path)) {
  stop("Baseline results fixture not found. Run the fixture generation script before testing.")
}

baseline <- readRDS(baseline_path)

if (!is.list(baseline) || length(baseline) == 0L) {
  stop("Baseline fixture must be a non-empty list of results.")
}

strip_names <- function(x) {
  if (is.list(x)) {
    lapply(x, strip_names)
  } else {
    names(x) <- NULL
    x
  }
}

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

compare_results <- function(actual, expected, scenario) {
  actual_clean <- strip_names(actual)
  expected_clean <- strip_names(expected)
  are_equal <- identical(actual_clean, expected_clean)
  if (!are_equal) {
    stop(sprintf("Results changed for %s", scenario))
  }
}

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  res <- do.call(maive, c(list(dat = scenario$data), scenario$args))
  expected <- baseline[[scenario_name]]
  compare_results(res, expected, scenario_name)
}

cat("All MAIVE regression equivalence tests passed.\n")

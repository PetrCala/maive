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

summarise_differences <- function(actual, expected, tolerance = 1e-8, limit = 10L) {
  missing_marker <- new.env(parent = emptyenv())
  results <- character()
  hidden_count <- 0L

  format_path <- function(path) {
    if (length(path) == 0) {
      "<root>"
    } else {
      paste(path, collapse = " -> ")
    }
  }

  format_value <- function(x, max_items = 4L) {
    if (identical(x, missing_marker)) {
      "<missing>"
    } else if (is.null(x)) {
      "NULL"
    } else if (is.list(x)) {
      sprintf("list(%d)", length(x))
    } else if (length(x) == 0) {
      "[]"
    } else {
      vals <- x
      if (is.numeric(vals)) {
        vals <- format(signif(vals, 6), trim = TRUE, scientific = FALSE)
      } else if (is.logical(vals)) {
        vals <- ifelse(is.na(vals), "NA", ifelse(vals, "TRUE", "FALSE"))
      } else if (is.character(vals)) {
        vals <- ifelse(is.na(vals), "NA", sprintf("\"%s\"", vals))
      }
      if (length(vals) > max_items) {
        vals <- c(vals[seq_len(max_items)], sprintf("... (+%d)", length(vals) - max_items))
      }
      paste(vals, collapse = ", ")
    }
  }

  append_line <- function(line) {
    if (length(results) < limit) {
      results <<- c(results, line)
    } else {
      hidden_count <<- hidden_count + 1L
    }
  }

  format_diff_detail <- function(path, expected, actual) {
    base <- sprintf(
      "%s expected %s, actual %s",
      format_path(path),
      format_value(expected),
      format_value(actual)
    )

    if (is.numeric(expected) && is.numeric(actual) && length(expected) == 1 && length(actual) == 1) {
      if (!is.na(expected) && !is.na(actual) && !is.nan(expected) && !is.nan(actual)) {
        delta <- actual - expected
        if (!is.infinite(delta)) {
          delta_fmt <- format(signif(delta, 6), trim = TRUE, scientific = FALSE)
          base <- sprintf("%s (Δ %s)", base, delta_fmt)
        }
      }
    }

    base
  }

  values_equal <- function(a, b) {
    if (identical(a, b)) {
      TRUE
    } else if (is.numeric(a) && is.numeric(b)) {
      if (is.nan(a) && is.nan(b)) {
        TRUE
      } else if (is.infinite(a) && is.infinite(b) && identical(sign(a), sign(b))) {
        TRUE
      } else if (is.na(a) && is.na(b)) {
        TRUE
      } else if (!is.na(a) && !is.na(b)) {
        abs(a - b) <= tolerance
      } else {
        FALSE
      }
    } else {
      FALSE
    }
  }

  compare_values <- function(act, exp, path) {
    if (identical(act, exp)) {
      return()
    }

    if (is.null(act) || is.null(exp)) {
      append_line(sprintf(
        "%s expected %s, actual %s",
        format_path(path),
        format_value(exp),
        format_value(act)
      ))
      return()
    }

    if (is.list(act) || is.list(exp)) {
      if (!is.list(act) || !is.list(exp)) {
        append_line(sprintf(
          "%s type mismatch: expected %s, actual %s",
          format_path(path),
          paste(class(exp), collapse = "/"),
          paste(class(act), collapse = "/")
        ))
        return()
      }

      act_len <- length(act)
      exp_len <- length(exp)
      max_len <- max(act_len, exp_len)
      act_names <- names(act)
      exp_names <- names(exp)

      for (i in seq_len(max_len)) {
        key <- NULL
        if (!is.null(act_names) && i <= length(act_names) && nzchar(act_names[i])) {
          key <- act_names[i]
        } else if (!is.null(exp_names) && i <= length(exp_names) && nzchar(exp_names[i])) {
          key <- exp_names[i]
        } else {
          key <- paste0("[", i, "]")
        }

        act_item <- if (i <= act_len) act[[i]] else missing_marker
        exp_item <- if (i <= exp_len) exp[[i]] else missing_marker
        child_path <- c(path, key)

        if (identical(act_item, missing_marker)) {
          append_line(sprintf(
            "%s missing from actual; expected %s",
            format_path(child_path),
            format_value(exp_item)
          ))
        } else if (identical(exp_item, missing_marker)) {
          append_line(sprintf(
            "%s missing from expected; actual %s",
            format_path(child_path),
            format_value(act_item)
          ))
        } else {
          compare_values(act_item, exp_item, child_path)
        }
      }

      return()
    }

    if (!is.atomic(act) || !is.atomic(exp)) {
      append_line(sprintf(
        "%s type mismatch: expected %s, actual %s",
        format_path(path),
        paste(class(exp), collapse = "/"),
        paste(class(act), collapse = "/")
      ))
      return()
    }

    act_len <- length(act)
    exp_len <- length(exp)

    if (act_len != exp_len) {
      append_line(sprintf(
        "%s length mismatch: expected %d, actual %d",
        format_path(path),
        exp_len,
        act_len
      ))
    }

    len <- min(act_len, exp_len)
    if (len == 0) {
      return()
    }

    for (i in seq_len(len)) {
      child_path <- if (len == 1) path else c(path, paste0("[", i, "]"))
      a_val <- act[[i]]
      e_val <- exp[[i]]
      if (!values_equal(a_val, e_val)) {
        append_line(format_diff_detail(child_path, e_val, a_val))
      }
    }
  }

  compare_values(actual, expected, character())

  if (hidden_count > 0) {
    results <- c(results, sprintf("... %d additional differences not shown", hidden_count))
  }

  results
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

compare_results <- function(actual, expected, scenario, tolerance = 1e-8) {
  actual_clean <- strip_names(actual)
  expected_clean <- strip_names(expected)
  if (identical(actual_clean, expected_clean)) {
    return()
  }

  diff_lines <- summarise_differences(actual_clean, expected_clean, tolerance = tolerance)
  if (length(diff_lines) == 0) {
    return()
  }

  message_lines <- c(
    sprintf("Results changed for %s:", scenario),
    paste0("  - ", diff_lines)
  )
  stop(paste(message_lines, collapse = "\n"), call. = FALSE)
}

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  res <- do.call(maive, c(list(dat = scenario$data), scenario$args))
  expected <- baseline[[scenario_name]]
  compare_results(res, expected, scenario_name)
}

cat("All MAIVE regression equivalence tests passed.\n")

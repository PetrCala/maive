test_that("maive_build_dummy_matrix produces one-hot encoded matrix", {
  values <- c("study_a", "study_b", "study_a", "study_c")
  result <- MAIVE:::maive_build_dummy_matrix(values)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(values))
  expect_equal(colnames(result), c("studyid.study_a", "studyid.study_b", "studyid.study_c"))
  expect_equal(result[, "studyid.study_a"], c(1, 0, 1, 0))
  expect_equal(result[, "studyid.study_b"], c(0, 1, 0, 0))
  expect_equal(result[, "studyid.study_c"], c(0, 0, 0, 1))
  expect_null(attr(result, "assign"))
  expect_null(attr(result, "contrasts"))
})

test_that("maive_build_dummy_matrix respects factor levels", {
  factor_values <- factor(c("x", "y", "x"), levels = c("y", "x"))
  result <- MAIVE:::maive_build_dummy_matrix(factor_values)

  expect_equal(colnames(result), c("studyid.y", "studyid.x"))
  expect_equal(result[, "studyid.y"], c(0, 1, 0))
  expect_equal(result[, "studyid.x"], c(1, 0, 1))
})

test_that("maive_build_dummy_matrix matches varhandle::to.dummy when available", {
  skip_if_not_installed("varhandle")

  library(varhandle)

  values <- c("A", "B", "A", "C", "B")
  expected <- varhandle::to.dummy(data.frame(studyid = values), "studyid")
  actual <- MAIVE:::maive_build_dummy_matrix(values)

  expect_identical(actual, expected)
})

test_that("maive_build_dummy_matrix returns a zero-column matrix for a single level", {
  result <- MAIVE:::maive_build_dummy_matrix(c("only", "only", "only"))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, 0L))

  centered <- MAIVE:::maive_center_dummy_matrix(c("only", "only", "only"))
  expect_equal(dim(centered), c(3L, 0L))
})

test_that("a single-level study identifier no longer hits the raw contrasts error", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6, 0.52, 0.48),
    sebs = c(0.25, 0.2, 0.22, 0.27, 0.21, 0.24),
    Ns = c(50, 80, 65, 90, 70, 60),
    study_id = rep("only", 6),
    stringsAsFactors = FALSE
  )

  # Without clustering the fit proceeds with no study dummies
  for (studylevel in c(0, 1)) {
    res <- suppressWarnings(maive(dat, method = 3, weight = 0, instrument = 1, studylevel = studylevel, SE = 2, AR = 0, first_stage = 0))
    expect_true(is.finite(res$beta))
  }

  # With clustering, clubSandwich reports its own readable message
  for (studylevel in c(2, 3)) {
    msg <- tryCatch(
      {
        suppressWarnings(maive(dat, method = 3, weight = 0, instrument = 1, studylevel = studylevel, SE = 2, AR = 0, first_stage = 0))
        ""
      },
      error = function(e) conditionMessage(e)
    )
    expect_match(msg, "single cluster")
    expect_false(grepl("contrasts can be applied", msg))
  }
})

test_that("euro results match the 0.2.5 reference at all four studylevels", {
  euro <- read.csv(test_path("fixtures", "euro.csv"), stringsAsFactors = FALSE)
  ref <- read.csv(
    test_path("fixtures", "euro_v0.2.5_reference.csv"),
    stringsAsFactors = FALSE, colClasses = "character"
  )
  ref$studylevel <- as.integer(ref$studylevel)
  ref$method <- as.integer(ref$method)

  # The reference was written with 17 significant digits and is bit-identical
  # to this code on the machine that produced it. Across BLAS builds the last
  # two or three digits move, so compare numerically with a tight tolerance.
  ref_tolerance <- 1e-8

  for (studylevel in 0:3) {
    for (method in 1:4) {
      res <- suppressWarnings(maive(
        euro,
        method = method, weight = 0, instrument = 1, studylevel = studylevel, SE = 2, AR = 0, first_stage = 0
      ))
      sub <- ref[ref$studylevel == studylevel & ref$method == method, ]
      for (i in seq_len(nrow(sub))) {
        field <- sub$field[i]
        # beta_standard at method 3 changed deliberately in 0.2.6 (it now comes
        # from the conventional PET-PEESE fit); it is covered by its own test.
        if (method == 3L && field == "beta_standard") next
        got <- if (field == "SE_instrumented_sum") sum(res$SE_instrumented) else res[[field]]
        expected <- sub$value[i]
        label <- sprintf("studylevel %d, method %d, %s", studylevel, method, field)
        if (is.na(expected)) {
          expect_true(is.na(got), label = label)
        } else if (is.numeric(got)) {
          expect_equal(got, as.numeric(expected), tolerance = ref_tolerance, label = label)
        } else {
          expect_identical(as.character(got), expected, label = label)
        }
      }
    }
  }
})

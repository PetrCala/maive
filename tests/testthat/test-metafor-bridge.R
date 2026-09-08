skip_if_not_installed("metafor")

# Ten two-group comparisons from five studies. Passing data = keeps the raw
# columns (n1i, n2i, study) on the escalc frame alongside yi and vi.
smd_raw <- function() {
  data.frame(
    m1i = c(2.1, 2.5, 1.9, 2.4, 2.2, 2.6, 2.3, 2.0, 2.4, 2.1),
    sd1i = c(1.0, 1.2, 0.9, 1.1, 1.0, 1.3, 1.1, 1.0, 1.2, 0.9),
    n1i = c(40, 50, 30, 45, 60, 35, 55, 48, 70, 38),
    m2i = c(1.8, 2.0, 1.7, 2.1, 1.9, 2.2, 2.0, 1.8, 2.0, 1.9),
    sd2i = c(1.0, 1.1, 1.0, 1.0, 1.1, 1.2, 1.0, 1.1, 1.0, 1.0),
    n2i = c(42, 55, 33, 48, 58, 37, 52, 50, 68, 41),
    study = rep(c("a", "b", "c", "d", "e"), each = 2),
    stringsAsFactors = FALSE
  )
}

smd_escalc <- function() {
  metafor::escalc(
    measure = "SMD",
    m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i,
    data = smd_raw()
  )
}

hand_built <- function(dat) {
  data.frame(
    bs = as.numeric(dat$yi),
    sebs = sqrt(as.numeric(dat$vi)),
    Ns = dat$n1i + dat$n2i
  )
}

run_maive <- function(d) {
  suppressWarnings(maive(d, method = 3, weight = 0, instrument = 1, studylevel = 0, SE = 0, AR = 0, first_stage = 0))
}

test_that("a two-group SMD escalc object converts without any extra input", {
  dat <- smd_escalc()
  expect_null(dat$ni) # metafor stores the total only as an attribute here

  out <- maive_from_metafor(dat)

  expect_identical(names(out), c("bs", "sebs", "Ns"))
  expect_identical(out, hand_built(dat))
  expect_identical(out$Ns, as.numeric(attr(dat$yi, "ni")))
})

test_that("escalc and rma.uni entry points reproduce the hand-built frame and MAIVE results", {
  dat <- smd_escalc()
  fit <- metafor::rma(yi, vi, data = dat)

  from_escalc <- maive_from_metafor(dat)
  from_fit <- maive_from_metafor(fit)
  manual <- hand_built(dat)

  expect_identical(from_escalc, manual)
  expect_identical(from_fit, manual)

  res_manual <- run_maive(manual)
  res_fit <- run_maive(from_fit)
  expect_identical(res_fit$beta, res_manual$beta)
  expect_identical(res_fit$SE, res_manual$SE)
})

test_that("the variance cannot end up in the sebs slot", {
  dat <- smd_escalc()
  out <- maive_from_metafor(dat)

  expect_equal(out$sebs^2, as.numeric(dat$vi))
  expect_false(isTRUE(all.equal(out$sebs, as.numeric(dat$vi))))

  # The mistake the adapter prevents: vi pasted as sebs shifts the estimate
  wrong <- data.frame(bs = as.numeric(dat$yi), sebs = as.numeric(dat$vi), Ns = out$Ns)
  expect_false(isTRUE(all.equal(run_maive(wrong)$beta, run_maive(out)$beta)))
})

test_that("custom var.names on the escalc frame are honoured", {
  dat <- metafor::escalc(
    measure = "COR", ri = c(0.30, 0.40, 0.20, 0.35, 0.25, 0.45),
    ni = c(50, 60, 70, 80, 90, 100), var.names = c("effect", "variance")
  )
  out <- maive_from_metafor(dat)

  expect_identical(out$bs, as.numeric(dat$effect))
  expect_identical(out$sebs, sqrt(as.numeric(dat$variance)))
  expect_identical(out$Ns, c(50, 60, 70, 80, 90, 100))
})

test_that("an explicit ni argument overrides the object and accepts a column name", {
  dat <- smd_escalc()
  dat$total <- dat$n1i + dat$n2i + 10

  expect_identical(maive_from_metafor(dat, ni = "total")$Ns, as.numeric(dat$total))
  expect_identical(maive_from_metafor(dat, ni = rep(99, nrow(dat)))$Ns, rep(99, nrow(dat)))
  expect_error(maive_from_metafor(dat, ni = 1:3), "length 3")
})

test_that("study_id can be a vector or a column name", {
  dat <- smd_escalc()

  by_name <- maive_from_metafor(dat, study_id = "study")
  by_vector <- maive_from_metafor(dat, study_id = dat$study)

  expect_identical(by_name$study_id, dat$study)
  expect_identical(by_name, by_vector)

  # The identifier then drives clustering in maive() without a warning
  expect_no_warning(
    maive(by_name, method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
  )
})

test_that("rma.uni rows pass through the fit's subset and NA masks", {
  dat <- smd_escalc()
  dat$yi[3] <- NA
  fit <- suppressWarnings(metafor::rma(yi, vi, data = dat, subset = -10))

  expect_equal(fit$k, 8)
  out <- maive_from_metafor(fit, study_id = "study")
  out_vec <- maive_from_metafor(fit, study_id = dat$study)

  kept <- c(1, 2, 4, 5, 6, 7, 8, 9)
  expect_identical(out$bs, as.numeric(dat$yi[kept]))
  expect_identical(out$sebs, sqrt(as.numeric(dat$vi[kept])))
  expect_identical(out$Ns, as.numeric(dat$n1i + dat$n2i)[kept])
  expect_identical(out$study_id, dat$study[kept])
  expect_identical(out, out_vec)

  # A vector already aligned to the kept rows is accepted as is
  expect_identical(maive_from_metafor(fit, study_id = dat$study[kept])$study_id, dat$study[kept])
  expect_error(maive_from_metafor(fit, study_id = letters[1:4]), "expected 8 .* or 10")
})

test_that("missing effects or variances are dropped from escalc input with a message", {
  dat <- smd_escalc()
  dat$vi[2] <- NA

  expect_message(out <- maive_from_metafor(dat, study_id = "study"), "Dropped 1 row")
  expect_equal(nrow(out), 9)
  expect_identical(out$study_id, dat$study[-2])
  expect_identical(out$Ns, as.numeric(dat$n1i + dat$n2i)[-2])
})

test_that("objects without a recoverable sample size give a clear error", {
  dat <- metafor::escalc(measure = "GEN", yi = c(0.1, 0.2, 0.3, 0.4), vi = c(0.01, 0.02, 0.03, 0.04))
  expect_error(maive_from_metafor(dat), "No sample sizes found")
  expect_error(maive_from_metafor(dat), "never infers")

  fit <- metafor::rma(yi, vi, data = dat)
  expect_error(maive_from_metafor(fit), "No sample sizes found")

  # Supplying ni resolves it, with rows aligned to the fit
  expect_identical(maive_from_metafor(fit, ni = c(10, 20, 30, 40))$Ns, c(10, 20, 30, 40))
})

test_that("a plain data frame with yi and vi columns is accepted", {
  plain <- data.frame(yi = c(0.1, 0.2, 0.3), vi = c(0.01, 0.02, 0.03), ni = c(10, 20, 30))
  out <- maive_from_metafor(plain)
  expect_identical(out, data.frame(bs = c(0.1, 0.2, 0.3), sebs = sqrt(c(0.01, 0.02, 0.03)), Ns = c(10, 20, 30)))

  expect_error(maive_from_metafor(data.frame(a = 1:4)), "missing column")
})

test_that("rma.mv and other rma subclasses are refused rather than flattened", {
  dat <- smd_escalc()
  dat$id <- seq_len(nrow(dat))
  mv <- suppressWarnings(metafor::rma.mv(yi, vi, random = ~ 1 | study / id, data = dat))
  expect_s3_class(mv, "rma.mv")
  expect_error(maive_from_metafor(mv), "Only .*rma.uni.* fits are supported; got .*rma.mv")

  mh <- metafor::rma.mh(
    measure = "OR",
    ai = c(4, 6, 3, 5), bi = c(36, 44, 27, 40), ci = c(11, 29, 11, 19), di = c(29, 21, 19, 26)
  )
  expect_error(maive_from_metafor(mh), "got .*rma.mh")
})

test_that("unsupported inputs are rejected by class", {
  expect_error(maive_from_metafor(1:10), "must be an .*escalc.* data frame or an .*rma.uni")
  expect_error(maive_from_metafor(list(yi = 1, vi = 1)), "must be an")
})

test_that("trimfill fits are refused before the generic rma.uni gate", {
  dat <- smd_escalc()
  fit <- metafor::rma(yi, vi, data = dat)
  tf <- metafor::trimfill(fit)
  expect_s3_class(tf, "rma.uni.trimfill")
  expect_true(inherits(tf, "rma.uni")) # would pass the generic gate

  expect_error(maive_from_metafor(tf), "trimfill.* fits contain imputed studies")
  expect_error(maive_from_metafor(tf), "Pass the original[[:space:]]+.*rma.* fit")

  # The original fit still converts
  expect_identical(maive_from_metafor(fit), hand_built(dat))
})

test_that("a stale ni attribute after dplyr::arrange() falls through to n1i + n2i", {
  skip_if_not_installed("dplyr")
  dat <- smd_escalc()
  arranged <- dplyr::arrange(dat, yi)
  correct <- as.numeric(arranged$n1i + arranged$n2i)

  # dplyr reorders the rows but leaves the attribute in the original order
  stale <- as.numeric(attr(arranged$yi, "ni"))
  expect_length(stale, nrow(arranged))
  expect_false(identical(stale, correct))

  out <- maive_from_metafor(arranged)
  expect_identical(out$Ns, correct)
  expect_identical(out$bs, as.numeric(arranged$yi))

  # rma() copies the same stale vector into fit$ni, so the rma entry point
  # needs the same consistency check
  fit <- metafor::rma(yi, vi, data = arranged)
  expect_false(identical(as.numeric(fit$ni), correct))
  expect_identical(maive_from_metafor(fit)$Ns, correct)
})

test_that("a stale ni attribute after a base reorder is caught by the consistency check", {
  dat <- smd_escalc()
  ord <- order(dat$yi)
  reordered <- as.data.frame(dat)[ord, ]
  # Simulate an attribute that did not follow the rows
  attr(reordered$yi, "ni") <- as.numeric(dat$n1i + dat$n2i)
  correct <- as.numeric(reordered$n1i + reordered$n2i)
  expect_false(identical(as.numeric(attr(reordered$yi, "ni")), correct))

  expect_identical(maive_from_metafor(reordered)$Ns, correct)

  # A stale ni column is subject to the same rule
  with_col <- reordered
  with_col$ni <- as.numeric(dat$n1i + dat$n2i)
  expect_identical(maive_from_metafor(with_col)$Ns, correct)
})

test_that("an ni column or attribute that agrees with n1i + n2i is still used", {
  dat <- smd_escalc()
  totals <- as.numeric(dat$n1i + dat$n2i)
  expect_identical(maive_from_metafor(dat)$Ns, totals)

  # Agreement is only required where both are present
  partial <- as.data.frame(dat)
  partial$n2i[3] <- NA
  attr(partial$yi, "ni") <- totals
  expect_identical(maive_from_metafor(partial)$Ns, totals)

  # Without n1i/n2i there is nothing to check against, so the attribute wins
  alone <- as.data.frame(dat)[, c("yi", "vi")]
  attr(alone$yi, "ni") <- totals
  expect_identical(maive_from_metafor(alone)$Ns, totals)

  # The explicit ni argument is never second-guessed
  expect_identical(maive_from_metafor(dat, ni = rep(7, nrow(dat)))$Ns, rep(7, nrow(dat)))
})

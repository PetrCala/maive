# MAIVE 0.3.1

*Unreleased*

## Bug Fixes

* The levels first stage (`first_stage = 0`) can fit a negative variance for an individual estimate. Previously `sqrt()` turned it into `NaN`, `lm()` dropped the row from PET but not from PEESE, the F-test, or the Hausman comparison, so the reported statistics used different samples, the AR interval came back `NA`, adjusted-weight AR aborted, and EK aborted. Such estimates are now excluded explicitly from everything built on the instrumented standard error (weights, PET, PEESE, PET-PEESE, EK, the F-test, the Hausman comparison, the Anderson-Rubin intervals, the wild bootstrap), so every reported statistic uses the same rows. A warning names the count and points to `first_stage = 1` (the log first stage), which cannot fit a negative variance. The fitted variance is not floored: a floored value would enter PET as an almost perfectly precise observation and, under `weight = 2`, as a weight of order `1/eps` (#24).

## New Features

* `maive()` and `waive()` return `n_excluded` (the number of excluded estimates, 0 with the log first stage) and `excluded_rows` (their positions in the input after empty rows are dropped). `SE_instrumented` and `weights` keep the input length with `NA` (never `NaN`) in the excluded positions.

---


# MAIVE 0.3.0

*Released: 2026-09-03*

## New Features

* `maive_from_metafor()` converts a metafor `escalc` data frame or an `rma.uni` fit into the data frame `maive()` and `waive()` expect. It takes the standard error as `sqrt(vi)`, so metafor's sampling variance cannot land in the `sebs` slot by mistake (which silently shifts the estimate), and it resolves sample sizes from `ni`, the `ni` column, the `ni` attribute `escalc()` stamps on the effect sizes, or `n1i + n2i`, never from the variance. `rma.uni` rows are taken through the fit's `subset` and missing-value masks so effects, variances, sample sizes, and study identifiers stay aligned. `rma.mv`, `rma.glmm`, and other `rma` subclasses are refused with a message. metafor is a suggested dependency only.

---


# MAIVE 0.2.6

*Released: 2026-09-03*

## Bug Fixes

* Custom column arguments now work: `estimate`, `se`, `n`, and `study_id` are resolved before the data frame is validated, so the documented custom column example in the vignette runs (it previously failed with "Missing required columns: bs, sebs, Ns"). Completely empty rows are still dropped before resolution.
* The `seed` argument now reaches the wild bootstrap (`SE = 3`). Previously the helper always used seed 123, so different seeds gave identical bootstrap confidence intervals and `seed = NULL` did not use the current RNG state. The reported SE under `SE = 3` is the CR1 cluster-robust SE and is unchanged; only the bootstrap confidence intervals (for example `egger_boot_ci`) depend on the seed, and a non-default seed now changes them.
* A study identifier with a single level no longer fails with the raw "contrasts can be applied only to factors with 2 or more levels" error. The dummy matrix is empty in that case, so `studylevel = 0` and `1` run, and `studylevel = 2` and `3` surface clubSandwich's own message that clustering needs more than one cluster. Results on multi-study data are unchanged.
* `beta_standard` at `method = 3` (PET-PEESE) now comes from the same conventional fit as `SE_standard`. It was previously read from the auxiliary PET-PEESE model that uses MAIVE's own weights, so the returned pair mixed two regressions and `beta_standard` moved with `weight` while `SE_standard` did not. This is a deliberate change to the reported value at `method = 3`; methods 1, 2, and 4 are unaffected, and the Hausman statistic still uses the auxiliary pair and is unchanged.

## New Features

* `maive()` and `waive()` return `ek_structure` ("kink", "linear", or "intercept") for `method = 4`, so an intercept-only degenerate EK fit is identifiable directly rather than inferred from a zero slope coefficient. It is `NA` for other methods.
* When no `study_id` argument is given and no column is named `study_id`, using the fourth column as the study identifier now emits a warning naming the column. Pass `study_id = "<column>"` to confirm the mapping, or drop the column if it is not a study identifier. A column named `study_id` is used regardless of its position.

## Documentation

* The vignette's custom column example is evaluated at build time with a fixture that satisfies the degrees of freedom rule (at least the number of unique studies plus three rows).
* `maive()` documents the `study_id` fallback and the origin of `beta_standard` and `SE_standard`.

---


# MAIVE 0.2.5

*Released: 2026-08-05*

## Bug Fixes

* Return unrounded estimates from maive() and waive(); values were previously rounded to 3 decimals, which distorted downstream z-ratios, p-values, and confidence intervals and could zero out small standard errors

---


# MAIVE 0.2.4

*Released: 2026-02-04*

## Bug Fixes

* Easymeta.org links

---


# MAIVE 0.2.3

*Released: 2026-02-04*

## Internal

* Add warnings for weak instruments

---


# MAIVE 0.2.2

*Released: 2026-01-07*

## Internal

* Update the ar calculation to build the weighted residual correctly

---


# MAIVE 0.2.1

*Released: 2026-01-07*

## New Features

* Add the option to set an RNG seed at the highest function level

---


# MAIVE 0.2.0

*Released: 2026-01-07*

## New Features

* Add the ottawa conference slides link to strategic locations around the package
* Add an explicit citation file to the inst folder, update CLAUDE.md with instructions
* Add the funnel plot vignette; add vignette preview


## Bug Fixes

* Failing tests
* A couple of failing tests
* Makefile targets, update docs
* Add a missing funnel plot topic to pkgdown
* Funnel plot documentation


## Documentation

* Update the introduction vignette to include info on the available column mapping


## Internal

* Re-generate R docs
* Allow positional arguments in the main functions
* Add empty rows removal
* Add a validation module, update the main functions to utilize it
* Regenerate r docs
* Update references to WAIVE to highlight its more aggressive correction for phacking
* Update outdated paper references (2024 -> 2025, nature communications)
* Update the docs to feature links to easymeta.org
* Clean up more unused special characters
* Get rid of funnel_plot docstring special characters

---


# MAIVE 0.1.12

*Released: 2026-01-07*

## Bug Fixes

* Failing tests with a gt operator


## Internal

* Add an explicit function for generating a funnel plot

---


# MAIVE 0.1.11

*Released: 2025-12-18*

## Other Changes

* Disable instrumentation when Ns has no variation; avoid aliased-slope vcov indexing
* Guard first-stage F-test against rank-deficient vcovCR; add regression test

---


# MAIVE 0.1.10

*Released: 2025-12-02*

## Internal

* Update outdated cran files

---


# MAIVE 0.1.9

*Released: 2025-12-02*

## Bug Fixes

* Cran submission issues

---


# MAIVE 0.1.8

*Released: 2025-11-27*

## Bug Fixes

* Further issues in the ar calculation
* Ar SE usage


## Internal

* Fix outdated tests

---


# MAIVE 0.1.7

*Released: 2025-11-27*

## Internal

* Keep the ar tests more neutral
* Ar calculation - avoid banana projection

---


# MAIVE 0.1.6

*Released: 2025-11-26*

## New Features

* Add automatic news updates


## Bug Fixes

* Automatic news release


## Documentation

* Update release instructions docs

---


# MAIVE 0.0.4

## Initial CRAN submission

* Implemented MAIVE (Meta-Analysis Instrumental Variable Estimator) for addressing spurious precision in meta-analysis
* Core functions:
  * `maive()`: Main function implementing PET, PEESE, PET-PEESE, and Endogenous Kink (EK) methods
  * `waive()`: Robust extension with downweighting of spurious precision and outliers
* Features:
  * Instrumental variable approach using inverse sample sizes
  * Multiple weighting schemes (no weights, inverse-variance, MAIVE-adjusted, WAIVE)
  * Study-level correlation handling (fixed effects, clustering, or both)
  * Robust standard errors (CR0, CR1, CR2, wild bootstrap)
  * Anderson-Rubin confidence intervals for weak instruments
  * First-stage specification options (levels or log transformation)
  * Publication bias testing based on instrumented FAT
* Comprehensive test suite with 9 test files
* Documentation with examples and usage guidelines

# Convert metafor objects into MAIVE input data

Builds the data frame that
[`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md) and
[`waive()`](https://petrcala.github.io/MAIVE/reference/waive.md) expect
(columns `bs`, `sebs`, `Ns`, and optionally `study_id`) from a metafor
`escalc` data frame or an `rma.uni` model fit.

## Usage

``` r
maive_from_metafor(x, ni = NULL, study_id = NULL)
```

## Arguments

- x:

  An `escalc` data frame (or any data frame with `yi` and `vi` columns)
  or an `rma.uni` fit.

- ni:

  Optional sample sizes: a numeric vector, or the name of a column in
  `x` (or in the data stored with the fit). Overrides any sample sizes
  carried by `x`.

- study_id:

  Optional study identifiers: a vector, or the name of a column in `x`
  (or in the data stored with the fit). Vectors may have either the
  length of the original data or the number of rows the fit kept.

## Value

A data frame with columns `bs`, `sebs`, `Ns`, and, when identifiers were
supplied, `study_id`, ready for
[`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md) or
[`waive()`](https://petrcala.github.io/MAIVE/reference/waive.md).

## Details

The conversion takes the standard error as `sqrt(vi)`, so metafor's
sampling variance cannot end up in the `sebs` slot by mistake. Sample
sizes are never inferred from the variances; doing so would reintroduce
the spurious precision that MAIVE corrects for. They are resolved, in
order, from the `ni` argument, an `ni` column (or the fit's `ni`), the
`ni` attribute that `escalc()` stamps on the effect size column, and
finally the exact two-group total `n1i + n2i`. The column and the
attribute are only used when they agree with `n1i + n2i` wherever both
are present: the attribute does not follow a reorder such as
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)
and `rma()` copies it into the fit, so a stale vector of the right
length falls through to the two-group total instead of binding every
sample size to the wrong estimate. When none of these is available the
function stops and names what it needs.

For `rma.uni` fits the effect sizes, variances, sample sizes, and any
vector supplied through `ni` or `study_id` are taken through the same
`subset` and missing-value masks the fit applied, so the rows stay
aligned. Multivariate (`rma.mv`), GLMM (`rma.glmm`), and other `rma`
subclasses are refused rather than silently flattened, and so are
[`metafor::trimfill()`](https://wviechtb.github.io/metafor/reference/trimfill.html)
fits: their imputed studies would enter MAIVE as observed estimates, so
pass the original `rma()` fit.

## Examples

``` r
dat <- metafor::escalc(
  measure = "SMD",
  m1i = c(2.1, 2.5, 1.9, 2.4, 2.2, 2.6),
  sd1i = c(1.0, 1.2, 0.9, 1.1, 1.0, 1.3),
  n1i = c(40, 50, 30, 45, 60, 35),
  m2i = c(1.8, 2.0, 1.7, 2.1, 1.9, 2.2),
  sd2i = c(1.0, 1.1, 1.0, 1.0, 1.1, 1.2),
  n2i = c(42, 55, 33, 48, 58, 37)
)
maive_dat <- maive_from_metafor(dat)
head(maive_dat)
#>          bs      sebs  Ns
#> 1 0.2971772 0.2221442  82
#> 2 0.4321116 0.1976638 105
#> 3 0.2071041 0.2529363  63
#> 4 0.2834756 0.2085368  93
#> 5 0.2837736 0.1850656 118
#> 6 0.3166657 0.2372653  72

fit <- metafor::rma(yi, vi, data = dat)
identical(maive_from_metafor(fit), maive_dat)
#> [1] TRUE

result <- maive(maive_dat,
  method = 3, weight = 0, instrument = 1,
  studylevel = 0, SE = 0, AR = 0
)
#> Warning: Sample size (6) is small for IV estimation. Results may be unreliable. Consider
#> using instrument=0 for small samples.
```

# Getting Started with MAIVE

## Overview

MAIVE (Meta-Analysis Instrumental Variable Estimator) addresses a
fundamental problem in meta-analysis of observational research:
**spurious precision**. Traditional meta-analysis assigns more weight to
studies with lower standard errors, assuming higher precision. However,
in observational research, precision must be estimated and is vulnerable
to manipulation through practices like p-hacking to achieve statistical
significance.

For a guided, interactive workflow, visit <https://www.easymeta.org>.

This manipulation can invalidate:

- Inverse-variance weighting schemes
- Bias-correction methods like funnel plots
- Traditional publication bias corrections

MAIVE introduces an **instrumental variable approach** to limit bias
caused by spurious precision in meta-analysis.

## The Problem: Spurious Precision

In observational research, researchers can inadvertently or deliberately
manipulate their analyses to achieve statistically significant results.
This includes:

- Selective reporting of specifications
- Outcome switching
- Sample trimming
- Selective controls inclusion

These practices create **spuriously precise estimates** that appear more
reliable than they actually are. Traditional meta-analysis methods that
weight by inverse variance will overweight these manipulated studies,
leading to biased conclusions.

## The MAIVE Solution

MAIVE uses **instrumental variables** to correct for spurious precision:

1.  **First-stage regression**: Instruments the potentially manipulated
    standard errors using sample sizes (which researchers cannot easily
    manipulate)
2.  **Second-stage regression**: Uses the instrumented standard errors
    in meta-regression models

This approach provides:

- Robust meta-estimates that account for spurious precision
- Hausman-type tests comparing IV and OLS estimates
- Anderson-Rubin confidence intervals for weak instruments
- Publication bias tests based on instrumented standard errors

## Installation

``` r

# Install from CRAN (once published)
install.packages("MAIVE")

# Or install development version from GitHub
install.packages("devtools")
devtools::install_github("PetrCala/MAIVE")
```

``` r

library(MAIVE)
```

## Data Structure

The [`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md)
function accepts either the default column names or custom mappings:

| Column              | Default name | Description                         |
|---------------------|--------------|-------------------------------------|
| Estimate            | `bs`         | Primary estimates (effect sizes)    |
| Std. Error          | `sebs`       | Standard errors (must be \> 0)      |
| Sample size         | `Ns`         | Sample sizes (must be \> 0)         |
| Study ID (optional) | `study_id`   | Clustering/fixed effects identifier |

**Custom column names:** You can map your own column names using:

- `estimate` for the estimate column
- `se` for the standard error column
- `n` for the sample size column
- `study_id` for study identifiers (optional)

Example:

``` r

custom_dat <- data.frame(
  my_est = c(0.50, 0.60, 0.40, 0.55, 0.45, 0.52, 0.48, 0.58, 0.42, 0.50, 0.47, 0.56),
  my_se = c(0.20, 0.18, 0.25, 0.22, 0.24, 0.19, 0.23, 0.17, 0.26, 0.21, 0.24, 0.18),
  my_n = c(80, 120, 95, 110, 90, 130, 100, 140, 85, 105, 88, 125),
  my_study = c("A", "A", "B", "B", "C", "C", "D", "D", "A", "B", "C", "D")
)

result <- maive(
  dat = custom_dat,
  estimate = "my_est",
  se = "my_se",
  n = "my_n",
  study_id = "my_study",
  method = 3,
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)
#> Registered S3 method overwritten by 'clubSandwich':
#>   method    from    
#>   bread.mlm sandwich

result$beta
#> [1] 0.7136277
```

Validation rules still apply: required columns must be numeric,
non-missing, finite, and there must be at least 4 observations after
removing completely empty rows. With a study identifier, there must be
at least as many rows as unique studies plus three (the fixture above
has 12 rows for 4 studies).

**Study identifier fallback:** when `study_id` is not supplied, a column
named `study_id` is used wherever it sits. If there is no such column
and the data frame has four or more columns, the fourth column is used
as the study identifier and a warning names it. If column four is a
moderator or a year rather than a study identifier, drop it or name the
correct column explicitly, because the identifier drives the study
dummies and clustering at every `studylevel` other than 0.

**Coming from metafor:** if your effect sizes are already in a
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
data frame or an `rma()` fit,
[`maive_from_metafor()`](https://petrcala.github.io/MAIVE/reference/maive_from_metafor.md)
builds this data frame for you. It uses `sqrt(vi)` as the standard error
and takes the sample sizes from the object (`ni`, the `ni` attribute
that `escalc()` records, or `n1i + n2i`), never inferring them from the
variance. Rows of an `rma.uni` fit pass through the fit’s `subset` and
missing-value masks, so a `study_id` vector for the original data stays
aligned.

``` r

smd <- metafor::escalc(
  measure = "SMD",
  m1i = c(2.1, 2.5, 1.9, 2.4, 2.2, 2.6, 2.3, 2.0),
  sd1i = c(1.0, 1.2, 0.9, 1.1, 1.0, 1.3, 1.1, 1.0),
  n1i = c(40, 50, 30, 45, 60, 35, 55, 48),
  m2i = c(1.8, 2.0, 1.7, 2.1, 1.9, 2.2, 2.0, 1.8),
  sd2i = c(1.0, 1.1, 1.0, 1.0, 1.1, 1.2, 1.0, 1.1),
  n2i = c(42, 55, 33, 48, 58, 37, 52, 50)
)

maive_from_metafor(smd)
#>          bs      sebs  Ns
#> 1 0.2971772 0.2221442  82
#> 2 0.4321116 0.1976638 105
#> 3 0.2071041 0.2529363  63
#> 4 0.2834756 0.2085368  93
#> 5 0.2837736 0.1850656 118
#> 6 0.3166657 0.2372653  72
#> 7 0.2829629 0.1943881 107
#> 8 0.1885831 0.2025211  98
```

## Basic Usage

Let’s create a simple example dataset:

``` r

# Simulated meta-analysis data: 50 estimates from 10 studies, with reported
# precision driven by sample size (as the first stage assumes)
set.seed(123)
n_estimates <- 50
Ns <- sample(100:1000, n_estimates, replace = TRUE)

data <- data.frame(
  bs = rnorm(n_estimates, mean = 0.3, sd = 0.2),
  sebs = 2 / sqrt(Ns) * runif(n_estimates, min = 0.8, max = 1.2),
  Ns = Ns,
  study_id = rep(1:10, each = 5)
)

head(data)
#>          bs       sebs  Ns study_id
#> 1 0.2409857 0.09350867 514        1
#> 2 0.4790251 0.08013701 562        1
#> 3 0.4756267 0.11301618 278        1
#> 4 0.4643162 0.08107801 625        1
#> 5 0.4377281 0.12785564 294        1
#> 6 0.4107835 0.05867782 917        2
```

### Default MAIVE Estimation

The default MAIVE estimator uses PET-PEESE with instrumented standard
errors, no weights, cluster-robust standard errors, and wild bootstrap:

``` r

# Run MAIVE with defaults
result <- maive(
  dat = data,
  method = 3,      # PET-PEESE (default)
  weight = 0,      # No weights (default)
  instrument = 1,  # Instrument SEs (default)
  studylevel = 2,  # Cluster-robust (default)
  SE = 3,          # Wild bootstrap (default)
  AR = 1           # Anderson-Rubin CI (default)
)

# View key results
cat("MAIVE Estimate:", round(result$beta, 3), "\n")
#> MAIVE Estimate: 0.336
cat("MAIVE SE:", round(result$SE, 3), "\n")
#> MAIVE SE: 0.051
cat("Standard Estimate:", round(result$beta_standard, 3), "\n")
#> Standard Estimate: 0.348
cat("Hausman Test:", round(result$Hausman, 3), "\n")
#> Hausman Test: 0
cat("First-stage F-test:", round(result$`F-test`, 3), "\n")
#> First-stage F-test: 546.273
```

### Understanding the Output

The [`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md)
function returns a list with the following key elements (see
[`?maive`](https://petrcala.github.io/MAIVE/reference/maive.md) for the
full list):

- `beta`: MAIVE point estimate (corrected for spurious precision)
- `SE`: MAIVE standard error
- `beta_standard`, `SE_standard`: Conventional (non-IV) estimate and
  standard error for comparison
- `Hausman`: Hausman-type test comparing IV vs OLS estimates (high value
  suggests spurious precision)
- `F-test`: First-stage F-test of instrument strength
- `AR_CI`: Anderson-Rubin confidence interval (robust to weak
  instruments)
- `pub bias p-value`: p-value for the publication bias test based on the
  instrumented FAT
- `SE_instrumented`: Vector of instrumented standard errors, one per
  input row (`NA` for excluded estimates)
- `n_excluded`, `excluded_rows`: Number and positions of estimates
  excluded because the levels first stage fitted a non-positive variance
  for them (0 and `integer(0)` with the log first stage)
- `petpeese_selected`: Which of PET and PEESE was selected when
  `method = 3`
- `ek_structure`: Structure of the EK fit when `method = 4` (“kink”,
  “linear”, or “intercept”)

## Method Options

MAIVE supports multiple meta-regression methods:

### 1. FAT-PET (Precision-Effect Test)

``` r

result_pet <- maive(
  dat = data,
  method = 1,  # FAT-PET
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)

cat("PET Estimate:", round(result_pet$beta, 3), "\n")
```

### 2. PEESE (Precision-Effect Estimate with Standard Error)

``` r

result_peese <- maive(
  dat = data,
  method = 2,  # PEESE
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)

cat("PEESE Estimate:", round(result_peese$beta, 3), "\n")
```

### 3. PET-PEESE (Conditional Method)

PET-PEESE uses PET if the PET estimate is not significantly different
from zero, otherwise uses PEESE:

``` r

result_petpeese <- maive(
  dat = data,
  method = 3,  # PET-PEESE (default)
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)

cat("PET-PEESE Estimate:", round(result_petpeese$beta, 3), "\n")
```

### 4. Endogenous Kink (EK)

``` r

result_ek <- maive(
  dat = data,
  method = 4,  # EK
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 0  # AR not available for EK
)

cat("EK Estimate:", round(result_ek$beta, 3), "\n")
```

## Weighting Schemes

### No Weights (Default)

Unweighted regression, recommended when spurious precision is a concern:

``` r

result_noweight <- maive(
  dat = data,
  method = 3,
  weight = 0,  # No weights
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)
```

### Inverse-Variance Weights

Traditional meta-analysis weighting:

``` r

result_ivweight <- maive(
  dat = data,
  method = 3,
  weight = 1,  # Inverse-variance weights
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 0  # AR not available with weights
)
```

### MAIVE-Adjusted Weights

Uses instrumented standard errors for weighting:

``` r

result_maiveweight <- maive(
  dat = data,
  method = 3,
  weight = 2,  # MAIVE-adjusted weights
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)
```

## Study-Level Correlation

Control for study-level correlation when you have multiple estimates per
study:

``` r

# No study-level adjustment
result_none <- maive(data, method = 3, weight = 0, instrument = 1, 
                     studylevel = 0, SE = 0, AR = 1)

# Study fixed effects (demeaned)
result_fe <- maive(data, method = 3, weight = 0, instrument = 1,
                   studylevel = 1, SE = 1, AR = 0)  # AR not available with FE

# Cluster-robust standard errors
result_cluster <- maive(data, method = 3, weight = 0, instrument = 1,
                        studylevel = 2, SE = 3, AR = 1)

# Both fixed effects and clustering
result_both <- maive(data, method = 3, weight = 0, instrument = 1,
                     studylevel = 3, SE = 3, AR = 0)
```

## Standard Error Options

``` r

# CR0 (Huber-White)
result_cr0 <- maive(data, method = 3, weight = 0, instrument = 1,
                    studylevel = 2, SE = 0, AR = 1)

# CR1 (Standard empirical correction)
result_cr1 <- maive(data, method = 3, weight = 0, instrument = 1,
                    studylevel = 2, SE = 1, AR = 1)

# CR2 (Bias-reduced estimator)
result_cr2 <- maive(data, method = 3, weight = 0, instrument = 1,
                    studylevel = 2, SE = 2, AR = 1)

# Wild bootstrap (recommended, default)
result_boot <- maive(data, method = 3, weight = 0, instrument = 1,
                     studylevel = 2, SE = 3, AR = 1)
```

## First-Stage Specification

MAIVE allows two functional forms for the first-stage regression. The
log specification is the default; pass `first_stage = 0` to reproduce
the levels specification from the published paper.

### Log (Default)

Log-linear regression of log(sebs²) on log(Ns) with smearing
retransformation. The fitted variance is always positive, so no estimate
can drop out of the second stage:

``` r

result_log <- maive(data, method = 3, weight = 0, instrument = 1,
                    studylevel = 2, SE = 3, AR = 1, first_stage = 1)

cat("First-stage (log) F-test:", round(result_log$`F-test`, 3), "\n")
```

### Levels (Paper Specification)

Regresses variance (sebs²) on constant and 1/Ns. This is the
specification in Irsova et al. (2025):

``` r

result_levels <- maive(data, method = 3, weight = 0, instrument = 1,
                       studylevel = 2, SE = 3, AR = 1, first_stage = 0)

cat("First-stage (levels) F-test:", round(result_levels$`F-test`, 3), "\n")
```

The levels fit is an ordinary least squares regression, so an individual
fitted variance can be negative even when the intercept is not. Such an
estimate has no instrumented standard error and is excluded from
everything built on it (weights, PET, PEESE, PET-PEESE, EK, the F-test,
the Hausman comparison, the Anderson-Rubin intervals, the bootstrap), so
every reported statistic uses the same rows. A warning reports the
count, which is also returned as `n_excluded` with the positions in
`excluded_rows`. The log specification above, the default, cannot fit a
negative variance and never excludes an estimate.

## WAIVE: More Aggressive Correction

WAIVE (Weighted Adjusted Instrumental Variable Estimator) provides a
more aggressive correction for p-hacking and spurious precision by
combining variance instrumentation with exponential-decay downweighting
of studies with spurious precision or extreme outliers:

For technical details and methodology, see the [WAIVE
slides](https://meta-analysis.cz/waive_ottawa.pdf).

``` r

result_waive <- waive(
  dat = data,
  method = 3,
  weight = 0,
  instrument = 1,
  studylevel = 2,
  SE = 3,
  AR = 1
)

cat("WAIVE Estimate:", round(result_waive$beta, 3), "\n")
cat("WAIVE SE:", round(result_waive$SE, 3), "\n")
```

WAIVE is particularly useful when:

- You need a more aggressive correction for p-hacking beyond standard
  MAIVE
- You suspect extreme outliers in your data
- Standard errors may be severely manipulated
- You want automatic downweighting of both spuriously precise estimates
  and outliers

## Interpretation Guidelines

### Hausman Test

The Hausman test compares the MAIVE (IV) estimate with the standard
(OLS) estimate:

- **High value**: Suggests spurious precision is a problem; MAIVE
  estimate is preferred
- **Low value**: IV and OLS are similar; spurious precision may not be
  severe

### First-Stage F-test

Tests the strength of the instrument (inverse sample size):

- **F \> 10**: Strong instrument
- **F \< 10**: Weak instrument; use Anderson-Rubin CI

### Anderson-Rubin Confidence Interval

Provides inference robust to weak instruments. Always check this CI when
F-test is low.

### Publication Bias p-value

Tests for publication bias using instrumented FAT:

- **Low p-value**: Evidence of publication bias
- **High p-value**: Little evidence of publication bias

## References

Irsova, Z., Bom, P.R.D., Havranek, T., & Rachinger, H. (2025). Spurious
precision in meta-analysis of observational research. *Nature
Communications*, 16, 8454. <https://doi.org/10.1038/s41467-025-63261-0>

Keane, M., & Neal, T. (2023). Instrument strength in IV estimation and
inference: A guide to theory and practice. *Journal of Econometrics*,
235(2), 1625-1653.

## See Also

- [`?maive`](https://petrcala.github.io/MAIVE/reference/maive.md) for
  detailed parameter documentation
- [`?waive`](https://petrcala.github.io/MAIVE/reference/waive.md) for
  the robust WAIVE estimator
- Project website: <https://meta-analysis.cz/maive/>

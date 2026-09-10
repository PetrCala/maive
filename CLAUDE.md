# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

MAIVE (Meta-Analysis Instrumental Variable Estimator) is an R package
that implements instrumental variable approaches to limit bias caused by
spurious precision in meta-analysis. The package addresses the problem
where observational research precision can be manipulated through
p-hacking, invalidating traditional inverse-variance weighting and
bias-correction methods.

**Project Page**: <https://meta-analysis.cz/maive/> **GitHub**:
PetrCala/MAIVE

## Citation

**CRITICAL**: The Nature Communications (2025) paper is the canonical
reference for MAIVE and must be cited everywhere users see references to
the package or method.

**Official Citation**:

``` text
Irsova, Z., Bom, P.R.D., Havranek, T., & Rachinger, H. (2025).
Spurious precision in meta-analysis of observational research.
Nature Communications, 16, 8454. https://doi.org/10.1038/s41467-025-63261-0
```

**When updating references**:

- Always use the Nature Communications (2025) paper as the primary
  reference
- Update both `DESCRIPTION` file (Description field) and `inst/CITATION`
  file
- Ensure documentation, vignettes, and README files reference this paper
- Never reference working papers or preprints as the primary citation

## Development Commands

### Package Installation

``` r

# Install from GitHub
devtools::install_github("PetrCala/MAIVE")

# Install development dependencies
devtools::install_dev_deps()

# Load package in-place for development
pkgload::load_all()
```

### Running Tests

``` r
# Run all tests
pkgload::load_all()
testthat::test_dir("tests/testthat")

# Or from command line
Rscript -e 'pkgload::load_all(); testthat::test_dir("tests/testthat")'

# Run single test file
testthat::test_file("tests/testthat/test-maive_first_stage.R")
```

### Building Documentation

``` r

# Generate documentation from roxygen comments
devtools::document()

# View help for main function
?maive
```

### Releasing to CRAN

``` bash
# Ensure clean working directory, then bump version (auto-commits, tags, pushes)
make version-patch   # 0.1.4 -> 0.1.5 (bug fixes)
make version-minor   # 0.1.4 -> 0.2.0 (new features)
make version-major   # 0.1.4 -> 1.0.0 (breaking changes)
```

The tag push triggers the CRAN submission workflow automatically.

## Architecture

### Core Function: `maive()`

The main exported functions are
[`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md) and
[`waive()`](https://petrcala.github.io/MAIVE/reference/waive.md)
(defined in R/maivefunction.r), which perform meta-analysis using
instrumental variable estimation. Both functions share a common pipeline
via `maive_run_pipeline()`:

1.  **Input Validation** (`normalize_maive_options`): Validates and
    normalizes all input parameters and resolves column mappings
    (`estimate`, `se`, `n`, optional `study_id`)
2.  **Data Preparation** (`maive_prepare_data`): Prepares study-level
    data, handles clustering and fixed effects through dummy matrix
    centering
3.  **Variance Instrumentation**
    (`maive_compute_variance_instrumentation`): First-stage regression
    of variances on inverse sample sizes, supporting both levels and log
    transformations
4.  **Weighting** (`maive_compute_weights`): Computes weights based on
    weighting scheme (none, inverse-variance, MAIVE-adjusted, or WAIVE)
5.  **Shared Pipeline** (`maive_run_pipeline`): Common second-stage
    analysis used by both
    [`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md) and
    [`waive()`](https://petrcala.github.io/MAIVE/reference/waive.md)
6.  **Model Fitting** (`maive_fit_models`): Fits PET (Precision-Effect
    Test), PEESE (Precision-Effect Estimate with Standard Error), and
    combined models
7.  **Method Selection** (`maive_select_petpeese`, `maive_fit_ek`):
    Selects between PET/PEESE or fits Endogenous Kink (EK) models
8.  **Confidence Intervals** (`maive_compute_egger_ar_ci`): Computes
    Anderson-Rubin confidence intervals for weak instruments (handled in
    R/ar.r)
9.  **Bootstrap** (R/boot.r): Wild cluster bootstrap for standard errors
    when SE=3

### Data Structure

The [`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md)
function expects a data frame with:

- Column 1: `bs` - Primary estimates
- Column 2: `sebs` - Standard errors (must be \> 0)
- Column 3: `Ns` - Sample sizes (must be \> 0)
- Column 4 (optional): `study_id` - Study identification for
  clustering/fixed effects

### Key Parameters

- `method`: 1=PET, 2=PEESE, 3=PET-PEESE (default), 4=EK
- `weight`: 0=none (default), 1=inverse-variance, 2=MAIVE-adjusted
- `instrument`: 0=no, 1=yes (default) - whether to instrument standard
  errors
- `studylevel`: 0=none, 1=fixed effects, 2=cluster (default), 3=both
- `SE`: 0=CR0 (Huber-White), 1=CR1, 2=CR2, 3=wild bootstrap (default)
- `AR`: 0=no, 1=yes (default) - compute Anderson-Rubin CI (only for
  unweighted IV)
- `first_stage`: 1=log (default), 0=levels (published paper) -
  first-stage functional form

**Note**: WAIVE is available as a standalone
[`waive()`](https://petrcala.github.io/MAIVE/reference/waive.md)
function that provides a more aggressive correction for p-hacking and
spurious precision by downweighting both spuriously precise estimates
and outliers.

### Study-Level Correlation

The `studylevel` parameter encodes two dimensions via modular
arithmetic:

- `cluster = studylevel %/% 2` (integer division)
- `dummy = studylevel %% 2` (modulo)

Study fixed effects are **demeaned** so the intercept measures a grand
mean. This is implemented via `maive_center_dummy_matrix()` which
centers dummy variables and drops one column.

### Anderson-Rubin Confidence Intervals

The AR CI computation (R/ar.r) uses adaptive grid search with
performance optimizations:

- Grid resolution adapts to sample size and parameter uncertainty
- Pre-computes projection matrices for efficiency
- Disables AR for datasets \>5000 observations to avoid memory issues
- Uses tighter bounds (2.5 standard errors) for better performance
- Only available for unweighted IV estimators

### First-Stage Regression

Two functional forms for instrumenting variances (sebs²):

1.  **Log** (`first_stage=1`, default): Log-linear regression with
    smearing retransformation to handle heteroskedasticity; the fitted
    variance is always positive
2.  **Levels** (`first_stage=0`): Linear regression of sebs² on constant
    and 1/Ns; the specification in the published paper

The levels fit can give a negative fitted variance for an individual
row. `maive_apply_variance_exclusion()` drops such rows explicitly (no
flooring, see \#24) from the weights, the second stage, the F-test
(refit on the kept rows), Hausman, AR, and the bootstrap, warns with the
count, and returns it as `n_excluded`.

The log specification uses Duan’s smearing estimator to retransform
predictions back to levels.

### Output Structure

Returns a named list with:

- `beta`, `SE`: Main MAIVE point estimate and standard error
- `beta_standard`, `SE_standard`: Conventional (non-IV) estimate and
  standard error from the same fit
- `Hausman`: Hausman-type test statistic comparing IV vs OLS intercepts
- `F-test`: First-stage F-test (if instrumenting)
- `AR_CI`: Anderson-Rubin confidence interval (if AR=1 and conditions
  met)
- `SE_instrumented`: Instrumented standard errors vector (input length,
  `NA` for excluded estimates)
- `n_excluded`, `excluded_rows`: Estimates excluded because the levels
  first stage fitted a non-positive variance for them; the first stage
  is fitted on all rows, every downstream statistic uses the kept rows
  (#24)
- `pub bias p-value`: p-value for publication bias test based on
  instrumented FAT
- `petpeese_selected`, `ek_structure`: model selection metadata for
  method 3 and method 4

## Testing

Tests are organized in `tests/testthat/` with fixtures in
`tests/testthat/fixtures/`. Test files follow the naming convention
`test-{component}.R`.

Key test patterns:

- Tests access internal functions using `MAIVE:::function_name` syntax
- Tests manually reconstruct computation pipelines to verify individual
  components
- Fixtures contain datasets for regression testing (e.g.,
  `egger_ar_no_acceptance.csv`)

## Important Implementation Details

- Study fixed effects are always demeaned, so
  `maive_center_dummy_matrix()` drops one column after centering
- The Hausman test uses difference-in-estimators variance (conservative
  approach)
- AR CI is automatically disabled for weighted methods or non-IV
  estimation
- Wild bootstrap (SE=3) uses Rademacher multipliers at the cluster level
- All robust standard errors use the `clubSandwich` package with
  CR0/CR1/CR2 estimators
- Infinite bounds in AR computation are handled by the optimization in
  [`compute_AR_CI_optimized()`](https://petrcala.github.io/MAIVE/reference/compute_AR_CI_optimized.md)

## Dependencies

Core dependencies (see `Imports` in `DESCRIPTION`):

- `stats`: Base R functionality
- `cli`: User-facing messages, warnings, and errors
- `clubSandwich`: Cluster-robust variance estimation

Development dependencies:

- `testthat`: Testing framework
- `devtools`: Package development tools
- `pkgload`: Load package in development

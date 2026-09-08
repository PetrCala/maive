# MAIVE: Meta-Analysis Instrumental Variable Estimator

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/MAIVE)](https://CRAN.R-project.org/package=MAIVE)
[![R-CMD-check](https://github.com/PetrCala/MAIVE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PetrCala/MAIVE/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/PetrCala/MAIVE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/PetrCala/MAIVE?branch=main)
<!-- badges: end -->

**Spurious Precision in Meta-Analysis of Observational Research**  
by Zuzana Irsova, Pedro R. D. Bom, Tomas Havranek, and Heiko Rachinger

**Project Website**: <https://meta-analysis.cz/maive/>

---

## Overview

MAIVE addresses a fundamental problem in meta-analysis of observational research: **spurious precision**.

Traditional meta-analysis assigns more weight to studies with lower standard errors, assuming higher precision. However, in observational research, precision can be manipulated through p-hacking and other questionable research practices, invalidating:

- Inverse-variance weighting schemes
- Traditional bias-correction methods (funnel plots, trim-and-fill)
- Selection models for publication bias

MAIVE implements an **instrumental variable approach** to limit bias caused by spurious precision in meta-analysis.

## Installation

### From CRAN

```r
install.packages("MAIVE")
```

### Development version

```r
install.packages("devtools")
devtools::install_github("PetrCala/MAIVE")
```

### Load package

```r
library(MAIVE)
```

## Quick Start

```r
# Prepare your data
data <- data.frame(
  bs = c(...),        # Effect sizes
  sebs = c(...),      # Standard errors
  Ns = c(...),        # Sample sizes
  study_id = c(...)   # Study IDs (optional)
)

# Run MAIVE with defaults (PET-PEESE, instrumented SEs, no weights)
result <- maive(
  dat = data,
  method = 3,      # PET-PEESE (default)
  weight = 0,      # No weights (default)
  instrument = 1,  # Instrument SEs (default)
  studylevel = 2,  # Cluster-robust (default)
  SE = 3,          # Wild bootstrap (default)
  AR = 1           # Anderson-Rubin CI (default)
)

# View results
print(result$beta)        # MAIVE estimate
print(result$SE)          # Standard error
print(result$Hausman)     # Hausman test
print(result$`F-test`)    # First-stage F-test
```

## Data Structure

The `maive()` function expects a data frame with:

| Column | Label | Description |
|--------|-------|-------------|
| 1 | `bs` | Primary estimates (effect sizes) |
| 2 | `sebs` | Standard errors (must be > 0) |
| 3 | `Ns` | Sample sizes (must be > 0) |
| 4 | `study_id` | Study identification (optional, for clustering/fixed effects) |

Other column names can be mapped with the `estimate`, `se`, `n`, and `study_id`
arguments. A column named `study_id` is used as the study identifier wherever it
sits. If there is no such column and no `study_id` argument, the fourth column is
used with a warning that names it, so a moderator or year kept in column four
does not silently drive the study dummies and clustering. With a study identifier,
the data needs at least the number of unique studies plus three rows.

### Using metafor objects

If your effect sizes already live in a metafor `escalc()` data frame or an
`rma()` fit, convert them with `maive_from_metafor()` rather than rebuilding the
frame by hand. It takes `sqrt(vi)` as the standard error (pasting metafor's
variance into `sebs` raises no error but shifts the estimate) and reads sample
sizes from the object, never from the variance:

```r
dat <- metafor::escalc(measure = "SMD", m1i = ..., sd1i = ..., n1i = ...,
                       m2i = ..., sd2i = ..., n2i = ...)
result <- maive(maive_from_metafor(dat), method = 3, weight = 0, instrument = 1,
                studylevel = 0, SE = 3, AR = 1)
```

Rows of an `rma.uni` fit are taken through the fit's `subset` and missing-value
masks, so a `study_id` vector for the original data stays aligned. `rma.mv` and
`rma.glmm` fits are refused rather than flattened.

## Key Features

### Methods

- **PET** (Precision-Effect Test)
- **PEESE** (Precision-Effect Estimate with Standard Error)
- **PET-PEESE** (Conditional method, default)
- **EK** (Endogenous Kink)

### Weighting Schemes

- No weights (recommended when spurious precision is a concern)
- Inverse-variance weights
- MAIVE-adjusted weights (using instrumented SEs)
- **WAIVE** (more aggressive correction via `waive()` function - downweights spurious precision and outliers) - [Details](https://meta-analysis.cz/waive_ottawa.pdf)

### Robust Inference

- Study-level correlation (fixed effects, clustering, or both)
- Multiple SE estimators (CR0, CR1, CR2, wild bootstrap)
- Anderson-Rubin confidence intervals for weak instruments
- First-stage specification options (log transformation by default, or levels as in the published paper)

### Output

The function returns:

- MAIVE point estimate and standard error
- Standard (non-IV) estimate for comparison
- Hausman-type test statistic
- First-stage F-test of instrument strength
- Anderson-Rubin confidence interval
- Publication bias test p-value
- Instrumented standard errors

## Documentation

- **Getting Started Guide**: See `vignette("introduction")`
- **Function Reference**: `?maive` and `?waive`
- **Guided interactive workflow**: <https://www.easymeta.org>
- **Development Workflow**: `.github/DEVELOPMENT-WORKFLOW.md`
- **CRAN Submission**: `.github/CRAN-SUBMISSION.md`
- **Project Page**: <https://meta-analysis.cz/maive/>

## Example

```r
# Create example data
set.seed(123)
Ns <- sample(100:1000, 50, replace = TRUE)
data <- data.frame(
  bs = rnorm(50, mean = 0.3, sd = 0.2),
  sebs = 2 / sqrt(Ns) * runif(50, min = 0.8, max = 1.2),  # precision driven by sample size
  Ns = Ns,
  study_id = rep(1:10, each = 5)
)

# Run MAIVE
result <- maive(data, method = 3, weight = 0, instrument = 1, 
                studylevel = 2, SE = 3, AR = 1)

# Compare with standard estimate
cat("MAIVE Estimate:", result$beta, "\n")
cat("Standard Estimate:", result$beta_standard, "\n")
cat("Hausman Test:", result$Hausman, "\n")

# Use WAIVE for more aggressive correction (downweights spurious precision + outliers)
result_waive <- waive(data, method = 3, weight = 0, instrument = 1,
                      studylevel = 2, SE = 3, AR = 1)
cat("WAIVE Estimate:", result_waive$beta, "\n")
```

## Citation

If you use MAIVE in your research, please cite:

> Irsova, Z., Bom, P.R.D., Havranek, T., & Rachinger, H. (2025).
> Spurious precision in meta-analysis of observational research.
> Nature Communications, 16, 8454. <https://doi.org/10.1038/s41467-025-63261-0>

## References

Keane, M., & Neal, T. (2023). Instrument strength in IV estimation and inference: A guide to theory and practice. *Journal of Econometrics*, 235(2), 1625-1653. <https://doi.org/10.1016/j.jeconom.2022.12.009>

Tipton, E. (2015). Small sample adjustments for robust variance estimation with cluster-correlated data. *Psychological Methods*, 20(3), 375–389. <https://doi.org/10.1037/met0000019>

## Contributing

We welcome contributions! Please see our [GitHub repository](https://github.com/PetrCala/MAIVE) for:

- Bug reports and feature requests (use [Issues](https://github.com/PetrCala/MAIVE/issues))
- Code contributions (submit Pull Requests)
- Questions and discussions

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Authors

- **Zuzana Irsova** - Charles University, Prague
- **Pedro R. D. Bom** - University of Deusto, Bilbao  
- **Tomas Havranek** - Charles University, Prague; CEPR, London; METRICS, Stanford
- **Petr Cala** (Maintainer) - Charles University, Prague
- **Heiko Rachinger** - University of the Balearic Islands, Palma

---

**Questions?** Contact the maintainer or visit our [project website](https://meta-analysis.cz/maive/).

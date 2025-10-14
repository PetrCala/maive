# WAIVE: Weighted Adjusted Instrumental Variable Estimator

## Overview

WAIVE (Weighted Adjusted Instrumental Variable Estimator) is an extension of MAIVE that adds robust downweighting of potentially p-hacked or outlier studies. While MAIVE corrects for spurious precision by instrumenting standard errors, WAIVE goes further by applying exponential-decay weights in the second stage to reduce the influence of:

1. **Spuriously precise studies** - Studies with negative first-stage residuals (smaller variance than expected given sample size)
2. **Extreme outliers** - Studies with residuals far from typical values in either direction

## Motivation

In meta-analysis of observational research, some studies may report suspiciously precise estimates (p-hacking) or exhibit extreme behavior that doesn't follow the expected variance-sample size relationship. WAIVE identifies these studies through first-stage residuals and smoothly downweights their contribution to the meta-analytic estimate.

## Usage

WAIVE is available as a weighting option (`weight=3`) in the `maive()` function:

```r
library(MAIVE)

# Example data
dat <- data.frame(
  bs = c(0.5, 0.45, 0.55, 0.6),
  sebs = c(0.25, 0.2, 0.22, 0.27),
  Ns = c(50, 80, 65, 90)
)

# Run WAIVE-PET-PEESE
result <- maive(
  dat = dat,
  method = 3,        # PET-PEESE
  weight = 3,        # WAIVE weights
  instrument = 1,    # Use instrumented SEs
  studylevel = 0,    # No study-level clustering
  SE = 0,            # CR0 standard errors
  AR = 0,            # No Anderson-Rubin CI (disabled for weighted methods)
  first_stage = 0    # Levels specification
)
```

## Algorithm

### First Stage (Same as MAIVE)

WAIVE uses the same first-stage regression as MAIVE to instrument variances:

**Levels specification** (`first_stage=0`):
```
se² ~ β₀ + β₁(1/N)
```

**Log specification** (`first_stage=1`):
```
log(se²) ~ β₀ + β₁·log(N)
```

### Weight Computation (New in WAIVE)

After the first stage, WAIVE computes exponential-decay weights based on residuals:

```r
# 1. Extract residuals
nu <- residuals(first_stage_model)

# 2. Robust normalization using MAD
sigma <- 1.4826 * MAD(nu)
if (sigma == 0) sigma <- sd(nu) + 1e-12
z <- nu / sigma

# 3. Compute exponential-decay components
z_neg <- max(-z, 0)              # Penalize spurious precision
z_out <- max(|z| - 2, 0)         # Penalize extreme outliers (tau=2)

# 4. Apply exponential decay
w <- exp(-1.0 * z_neg - 0.25 * z_out²)

# 5. Floor to avoid zero leverage
w <- max(w, 0.05)

# 6. Normalize to mean 1
w <- w / mean(w)
```

### Interpretation

- **Negative residuals** (`z < 0`): Study is more precise than expected → downweighted linearly via `exp(-z_neg)`
- **Extreme residuals** (`|z| > 2`): Study is an outlier → downweighted quadratically via `exp(-0.25·z_out²)`
- **Weight floor** (0.05): Ensures no study is completely excluded

### Second Stage

The second-stage regression uses the instrumented standard errors (as in MAIVE) weighted by the WAIVE adjustment weights:

```r
# Combine instrumented variances with WAIVE weights
w_combined <- sqrt(sebs²_instrumented * w_waive)

# Second stage (weighted regression)
bs / w_combined ~ 1 + sebs_instrumented / w_combined
```

## Comparison with MAIVE

| Feature | MAIVE (weight=2) | WAIVE (weight=3) |
|---------|------------------|------------------|
| First-stage IV | ✓ | ✓ |
| Instrumented SEs | ✓ | ✓ |
| Robust downweighting | ✗ | ✓ |
| Outlier handling | ✗ | ✓ |
| Anderson-Rubin CI | ✗ (disabled) | ✗ (disabled) |

## Properties

1. **Consistency**: WAIVE reduces to MAIVE when all residuals are similar (no extreme values)
2. **Robustness**: Uses MAD for normalization, resistant to outliers
3. **Smoothness**: Exponential decay provides smooth downweighting rather than hard cutoffs
4. **Compatibility**: Works with all MAIVE features (methods, first-stage specs, clustering, bootstrap)

## When to Use WAIVE

Consider WAIVE when:
- You suspect p-hacking or selective reporting in your literature
- Your dataset contains extreme outliers or leverage points
- Standard MAIVE estimates are sensitive to influential studies
- You want a more robust meta-analytic estimate

WAIVE is particularly useful for meta-analyses of observational research where publication bias and questionable research practices are common.

## Technical Details

### First-Stage Specification

WAIVE supports both levels and log specifications:

- **Levels** (`first_stage=0`): Standard linear regression, negative intercept → intercept dropped
- **Log** (`first_stage=1`): Log-linear with Duan smearing retransformation for heteroskedasticity

The residuals used for weighting come from the chosen first-stage model.

### Anderson-Rubin Confidence Intervals

WAIVE automatically disables Anderson-Rubin CIs because:
1. AR CIs require unweighted estimation for theoretical validity
2. WAIVE weights break the AR projection geometry
3. This is consistent with MAIVE weighted methods (weight=1,2)

### Standard Errors

WAIVE supports all standard error options:
- `SE=0`: CR0 (Huber-White)
- `SE=1`: CR1 (standard correction)
- `SE=2`: CR2 (bias-reduced)
- `SE=3`: Wild cluster bootstrap

### Study-Level Effects

WAIVE works with all study-level correlation structures:
- `studylevel=0`: No clustering or fixed effects
- `studylevel=1`: Study fixed effects only
- `studylevel=2`: Study clustering only
- `studylevel=3`: Both fixed effects and clustering

## Examples

See `examples/waive_demo.R` for a complete demonstration including:
- Comparison of MAIVE vs WAIVE estimates
- Weight computation details
- Identification of downweighted studies

## Testing

Comprehensive tests are available in `tests/testthat/test-waive.R` covering:
- Weight computation correctness
- First-stage specification compatibility (levels/log)
- Downweighting of spurious precision
- Downweighting of extreme outliers
- Study clustering support
- AR CI disabling
- Weight floor enforcement
- Structural consistency with MAIVE

## References

WAIVE builds on the MAIVE methodology from:

> Irsova, Z., Bom, P., Havranek, T., & Rachinger, H. "Spurious Precision in Meta-Analysis of Observational Research"

For questions or issues, see: https://meta-analysis.cz/maive

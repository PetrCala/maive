# WAIVE Demo: Weighted Adjusted Instrumental Variable Estimator
# This example demonstrates how WAIVE extends MAIVE by downweighting
# studies with spurious precision or extreme outlier behavior

library(MAIVE)

# Create example data with some potentially p-hacked studies
dat <- data.frame(
  bs = c(0.5, 0.48, 0.52, 0.51, 0.49, 0.45, 0.55, 0.50),
  sebs = c(0.25, 0.22, 0.23, 0.24, 0.21, 0.05, 0.80, 0.26),
  # Study 6 is suspiciously precise (tiny SE with medium N)
  # Study 7 is an outlier (very imprecise)
  Ns = c(50, 60, 55, 58, 62, 45, 40, 52)
)

# Run standard MAIVE (weight=2)
result_maive <- maive(
  dat = dat,
  method = 3,        # PET-PEESE
  weight = 2,        # MAIVE adjusted weights
  instrument = 1,    # Use instrumented SEs
  studylevel = 0,    # No clustering
  SE = 0,            # CR0 standard errors
  AR = 0,            # No Anderson-Rubin CI
  first_stage = 0    # Levels specification
)

# Run WAIVE (weight=3)
result_waive <- maive(
  dat = dat,
  method = 3,        # PET-PEESE
  weight = 3,        # WAIVE weights
  instrument = 1,    # Use instrumented SEs
  studylevel = 0,    # No clustering
  SE = 0,            # CR0 standard errors
  AR = 0,            # No Anderson-Rubin CI (disabled for weighted methods)
  first_stage = 0    # Levels specification
)

# Compare results
cat("MAIVE estimate:", result_maive$beta, "\n")
cat("MAIVE SE:", result_maive$SE, "\n\n")

cat("WAIVE estimate:", result_waive$beta, "\n")
cat("WAIVE SE:", result_waive$SE, "\n\n")

cat("Difference:", abs(result_maive$beta - result_waive$beta), "\n")

# WAIVE with log first-stage
result_waive_log <- maive(
  dat = dat,
  method = 3,
  weight = 3,
  instrument = 1,
  studylevel = 0,
  SE = 0,
  AR = 0,
  first_stage = 1    # Log specification
)

cat("\nWAIVE with log first-stage:", result_waive_log$beta, "\n")

# Inspect how WAIVE weights are computed
cat("\n--- WAIVE Weighting Details ---\n")

# First stage regression
invNs <- 1 / dat$Ns
sebs2 <- dat$sebs^2
varreg <- lm(sebs2 ~ invNs)
nu <- residuals(varreg)

# Robust normalization
mad_value <- mad(nu, constant = 1.4826)
sigma <- if (mad_value == 0) sd(nu) + 1e-12 else mad_value
z <- nu / sigma

# Compute WAIVE weights
z_neg <- pmax(-z, 0)
z_out <- pmax(abs(z) - 2, 0)
w_waive <- exp(-1.0 * z_neg - 0.25 * z_out^2)
w_waive <- pmax(w_waive, 0.05)

cat("First-stage residuals (standardized):\n")
print(round(z, 2))
cat("\nWAIVE weights (before normalization):\n")
print(round(w_waive, 3))
cat("\nNotice: Study 6 (spuriously precise) gets z < 0 → downweighted\n")
cat("        Study 7 (outlier) gets |z| > 2 → downweighted\n")

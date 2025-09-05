#' @keywords internal
manual_wild_cluster_boot_se <- function(model, data, cluster_var, B = 500, seed = 123) {
  set.seed(seed)

  # Extract residuals and fitted values
  resids <- residuals(model)
  fitted_vals <- fitted(model)

  # Get cluster IDs
  clusters <- unique(as.character(data[[cluster_var]]))
  G <- length(clusters)

  # Coefficient names
  coef_names <- names(coef(model))
  k <- length(coef_names)

  # Matrix to store bootstrap coefficients
  boot_coefs <- matrix(NA, nrow = B, ncol = k)
  colnames(boot_coefs) <- coef_names

  # Loop over bootstrap replications
  for (b in 1:B) {
    # Draw Rademacher multipliers per cluster
    u_g <- sample(c(-1, 1), size = G, replace = TRUE)
    names(u_g) <- as.character(clusters)

    # Create bootstrap outcome
    data$y_boot <- fitted_vals + resids * u_g[as.character(data[[cluster_var]])]

    # Refit the same model formula on bootstrap sample
    form_boot <- update(formula(model), y_boot ~ .)
    fit_boot <- lm(form_boot, data = data)

    # Store bootstrap coefficients
    boot_coefs[b, ] <- coef(fit_boot)
  }

  # Compute bootstrap SEs
  boot_se <- apply(boot_coefs, 2, sd)

  # Compute bootstrap percentile CI
  alpha <- 0.05
  boot_ci <- t(apply(boot_coefs, 2, function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2))))
  colnames(boot_ci) <- c("lower", "upper")

  # Return list
  return(list(
    boot_se = boot_se,
    boot_ci = boot_ci,
    boot_coefs = boot_coefs
  ))
}

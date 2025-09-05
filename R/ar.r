#' @keywords internal
PET_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs

#' @keywords internal
PEESE_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs^2

#' @keywords internal
compute_AR_CI_optimized <- function(model, adjust_fun, bs, sebs, invNs, g, type_choice) {
  # Beta estimates and robust SEs
  beta0 <- model$coefficients[1]
  beta0se <- sqrt(vcovCR(model, cluster = g, type = type_choice)[1, 1])
  beta1 <- model$coefficients[2]
  beta1se <- sqrt(vcovCR(model, cluster = g, type = type_choice)[2, 2])

  # Adaptive grid resolution based on dataset size
  M <- length(bs)
  base_resolution <- min(50, max(20, round(sqrt(M))))

  # Tighter bounds for better performance
  l0 <- beta0 - 3 * beta0se
  u0 <- beta0 + 3 * beta0se
  l1 <- beta1 - 3 * beta1se
  u1 <- beta1 + 3 * beta1se

  # Adaptive grid resolution
  pr0 <- max(base_resolution, round(base_resolution * 10 / (u0 - l0)))
  pr1 <- max(base_resolution, round(base_resolution * 10 / (u1 - l1)))

  b0_grid <- seq(l0, u0, length.out = min(pr0, 200))
  b1_grid <- seq(l1, u1, length.out = min(pr1, 200))

  # Pre-compute matrices (reused across all grid points)
  ones_vec <- rep(1, M)
  Z <- cbind(ones_vec, invNs)
  ZtZ_inv <- solve(t(Z) %*% Z)
  PZ <- Z %*% ZtZ_inv %*% t(Z)
  MZ <- diag(M) - PZ

  # Pre-compute sebs powers for adjustment functions
  sebs_sq <- sebs^2

  # Vectorized computation using outer product
  compute_AR_stats_vectorized <- function(b0_vals, b1_vals) {
    n_b0 <- length(b0_vals)
    n_b1 <- length(b1_vals)
    K <- n_b0 * n_b1

    # Grid of (b0, b1) pairs — K columns
    grid <- expand.grid(b0 = b0_vals, b1 = b1_vals) # K rows

    # Replicate vectors into M × K matrices
    bs_mat <- matrix(rep(bs, times = K), nrow = M, ncol = K)
    se_mat <- if (identical(adjust_fun, PET_adjust)) {
      matrix(rep(sebs, times = K), nrow = M, ncol = K)
    } else {
      matrix(rep(sebs_sq, times = K), nrow = M, ncol = K)
    }
    b0_mat <- matrix(rep(grid$b0, each = M), nrow = M, ncol = K)
    b1_mat <- matrix(rep(grid$b1, each = M), nrow = M, ncol = K)

    # Vectorized adjustment (PET or PEESE)
    bs_star_mat <- bs_mat - b0_mat - b1_mat * se_mat # M × K

    # AR statistic pieces (M × M) %*% (M × K) -> (M × K)
    PZ_bs_star <- PZ %*% bs_star_mat
    MZ_bs_star <- MZ %*% bs_star_mat

    num <- colSums(bs_star_mat * PZ_bs_star) # length K
    denom <- colSums(bs_star_mat * MZ_bs_star) # length K
    denom[abs(denom) < 1e-10] <- 1e-10

    stats <- (M - 2) * num / denom # length K

    # Reshape back to n_b0 × n_b1 surface
    matrix(stats, nrow = n_b0, ncol = n_b1)
  }

  AR_stats <- compute_AR_stats_vectorized(b0_grid, b1_grid)

  # Check which are accepted
  AR_accept <- AR_stats < 5.99

  b0_CI_all <- b0_grid[rowSums(AR_accept) > 0]
  b1_CI_all <- b1_grid[colSums(AR_accept) > 0]

  b0_CI <- c(b0_CI_all[1], b0_CI_all[length(b0_CI_all)])
  b1_CI <- c(b1_CI_all[1], b1_CI_all[length(b1_CI_all)])

  list(
    b0_CI = round(b0_CI, 3),
    b1_CI = round(b1_CI, 3)
  )
}

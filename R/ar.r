#' @keywords internal
PET_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs

#' @keywords internal
PEESE_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs^2

#' @keywords internal
compute_AR_CI_optimized <- function(model, adjust_fun, bs, sebs, invNs, g, type_choice) {
  # Beta estimates and robust SEs
  beta0 <- model$coefficients[1]
  beta0se <- sqrt(clubSandwich::vcovCR(model, cluster = g, type = type_choice)[1, 1])
  beta1 <- model$coefficients[2]
  beta1se <- sqrt(clubSandwich::vcovCR(model, cluster = g, type = type_choice)[2, 2])

  M <- length(bs)

  # For extremely large datasets, disable AR computation to avoid memory issues
  if (M > 5000) {
    warning("Dataset too large for AR computation (", M, " observations). AR CI disabled.")
    return(list(b0_CI = "NA", b1_CI = "NA"))
  }

  # For very large datasets, use much smaller grids to avoid memory issues
  if (M > 1000) {
    # Very aggressive grid reduction for large datasets
    max_grid_size <- min(30, max(15, round(500 / sqrt(M))))
    base_resolution <- max_grid_size
  } else {
    base_resolution <- min(50, max(20, round(sqrt(M))))
  }

  # Tighter bounds for better performance
  l0 <- beta0 - 2.5 * beta0se # Even tighter bounds
  u0 <- beta0 + 2.5 * beta0se
  l1 <- beta1 - 2.5 * beta1se
  u1 <- beta1 + 2.5 * beta1se

  # Adaptive grid resolution with strict limits
  pr0 <- min(max(base_resolution, round(base_resolution * 5 / (u0 - l0))), 100)
  pr1 <- min(max(base_resolution, round(base_resolution * 5 / (u1 - l1))), 100)

  b0_grid <- seq(l0, u0, length.out = pr0)
  b1_grid <- seq(l1, u1, length.out = pr1)

  # Pre-compute matrices (reused across all grid points)
  ones_vec <- rep(1, M)
  Z <- cbind(ones_vec, invNs)
  ZtZ_inv <- solve(t(Z) %*% Z)
  PZ <- Z %*% ZtZ_inv %*% t(Z)
  MZ <- diag(M) - PZ

  # Pre-compute sebs powers for adjustment functions
  sebs_sq <- sebs^2

  # Memory-efficient computation: process in chunks to avoid large matrices
  compute_AR_stats_chunked <- function(b0_vals, b1_vals) {
    n_b0 <- length(b0_vals)
    n_b1 <- length(b1_vals)

    # Process in smaller chunks to avoid memory issues
    # For very large datasets, use even smaller chunks
    max_chunk_size <- if (M > 2000) 500 else 1000
    chunk_size <- min(max_chunk_size, n_b0 * n_b1)
    n_chunks <- ceiling((n_b0 * n_b1) / chunk_size)

    all_stats <- numeric(n_b0 * n_b1)

    for (chunk in 1:n_chunks) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_b0 * n_b1)
      chunk_size_actual <- end_idx - start_idx + 1

      # Get grid points for this chunk
      grid_indices <- start_idx:end_idx
      b0_chunk <- rep(b0_vals, each = n_b1)[grid_indices]
      b1_chunk <- rep(b1_vals, times = n_b0)[grid_indices]

      # Create matrices for this chunk only
      bs_mat <- matrix(rep(bs, times = chunk_size_actual), nrow = M, ncol = chunk_size_actual)
      se_mat <- if (identical(adjust_fun, PET_adjust)) {
        matrix(rep(sebs, times = chunk_size_actual), nrow = M, ncol = chunk_size_actual)
      } else {
        matrix(rep(sebs_sq, times = chunk_size_actual), nrow = M, ncol = chunk_size_actual)
      }
      b0_mat <- matrix(rep(b0_chunk, each = M), nrow = M, ncol = chunk_size_actual)
      b1_mat <- matrix(rep(b1_chunk, each = M), nrow = M, ncol = chunk_size_actual)

      # Vectorized adjustment (PET or PEESE)
      bs_star_mat <- bs_mat - b0_mat - b1_mat * se_mat

      # AR statistic pieces
      PZ_bs_star <- PZ %*% bs_star_mat
      MZ_bs_star <- MZ %*% bs_star_mat

      num <- colSums(bs_star_mat * PZ_bs_star)
      denom <- colSums(bs_star_mat * MZ_bs_star)
      denom[abs(denom) < 1e-10] <- 1e-10

      stats <- (M - 2) * num / denom
      all_stats[grid_indices] <- stats

      # Force garbage collection after each chunk to free memory
      if (chunk %% 5 == 0) {
        gc()
      }
    }

    # Reshape back to n_b0 × n_b1 surface
    matrix(all_stats, nrow = n_b0, ncol = n_b1)
  }

  AR_stats <- compute_AR_stats_chunked(b0_grid, b1_grid)

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

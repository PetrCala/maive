seed_fixture <- function() {
  data.frame(
    bs = c(
      0.3286, 0.1947, 0.3324, 0.4257, 0.3225, 0.4305, 0.3983, 0.3590, 0.5049, 0.4508,
      0.4818, 0.4675, 0.4254, 0.2748, 0.3756, 0.2952, 0.3108, 0.3137, 0.3220, 0.4265,
      0.2266, 0.2718, 0.3592, 0.3577, 0.4073, 0.3033, 0.1090, 0.2585, 0.4321, 0.2493
    ),
    sebs = c(
      0.0736, 0.1190, 0.0786, 0.1368, 0.1175, 0.0831, 0.1172, 0.0975, 0.1150, 0.0698,
      0.0799, 0.0669, 0.0641, 0.0875, 0.0694, 0.0711, 0.1186, 0.0894, 0.1253, 0.0868,
      0.0823, 0.0686, 0.0767, 0.0792, 0.0667, 0.0881, 0.1691, 0.0726, 0.1042, 0.0819
    ),
    Ns = c(
      270, 170, 343, 94, 160, 396, 157, 213, 225, 438,
      336, 349, 363, 403, 393, 452, 128, 398, 159, 418,
      443, 487, 451, 431, 339, 406, 61, 361, 278, 353
    ),
    study_id = rep(seq_len(10), each = 3)
  )
}

run_boot <- function(dat, seed) {
  suppressWarnings(maive(
    dat,
    method = 1, weight = 0, instrument = 1, studylevel = 2, SE = 3, AR = 0, first_stage = 0,
    seed = seed
  ))
}

test_that("seed reaches the wild bootstrap: different seeds change the bootstrap CI", {
  skip_if_not_installed("clubSandwich")
  dat <- seed_fixture()

  res_seed_1 <- run_boot(dat, seed = 1)
  res_seed_2 <- run_boot(dat, seed = 2)
  res_seed_1_again <- run_boot(dat, seed = 1)

  # The bootstrap CI is built from draw-dependent t quantiles
  expect_false(identical(res_seed_1$egger_boot_ci, res_seed_2$egger_boot_ci))
  expect_false(identical(res_seed_1$boot_result$boot_ci, res_seed_2$boot_result$boot_ci))

  # The same seed reproduces the CI exactly
  expect_identical(res_seed_1$egger_boot_ci, res_seed_1_again$egger_boot_ci)
  expect_identical(res_seed_1$boot_result$boot_ci, res_seed_1_again$boot_result$boot_ci)
  expect_identical(res_seed_1$boot_result$boot_t_stats, res_seed_1_again$boot_result$boot_t_stats)

  # The point estimate never depends on the seed
  expect_identical(res_seed_1$beta, res_seed_2$beta)
})

test_that("seed = NULL leaves the bootstrap on the current RNG state", {
  skip_if_not_installed("clubSandwich")
  set.seed(123)
  n_clusters <- 8
  cluster <- rep(seq_len(n_clusters), each = 5)
  x <- rnorm(length(cluster))
  y <- 1 + 0.5 * x + rnorm(n_clusters, sd = 0.3)[cluster] + rnorm(length(cluster), sd = 0.5)
  data <- data.frame(y = y, x = x, cluster = cluster)
  model <- lm(y ~ x, data = data)

  set.seed(99)
  boot_a <- manual_wild_cluster_boot_se(model, data, "cluster", B = 30, seed = NULL)
  set.seed(99)
  boot_b <- manual_wild_cluster_boot_se(model, data, "cluster", B = 30, seed = NULL)
  boot_c <- manual_wild_cluster_boot_se(model, data, "cluster", B = 30, seed = NULL)

  expect_identical(boot_a$boot_coefs, boot_b$boot_coefs)
  expect_false(identical(boot_b$boot_coefs, boot_c$boot_coefs))

  # NA is the normalized form of NULL used by maive_analyze()
  set.seed(99)
  boot_na <- manual_wild_cluster_boot_se(model, data, "cluster", B = 30, seed = NA_integer_)
  expect_identical(boot_na$boot_coefs, boot_a$boot_coefs)
})

test_that("maive_infer_coef reads the bootstrap seed option", {
  skip_if_not_installed("clubSandwich")
  set.seed(123)
  cluster <- rep(seq_len(8), each = 5)
  x <- rnorm(length(cluster))
  y <- 1 + 0.5 * x + rnorm(length(cluster), sd = 0.5)
  data <- data.frame(y = y, x = x, g = cluster)
  model <- lm(y ~ x, data = data)

  old <- getOption("MAIVE.bootstrap.seed")
  on.exit(options(MAIVE.bootstrap.seed = old), add = TRUE)

  options(MAIVE.bootstrap.seed = 7L)
  inf_7 <- MAIVE:::maive_infer_coef(model, 1L, 3L, data, "g", "CR0")
  options(MAIVE.bootstrap.seed = 8L)
  inf_8 <- MAIVE:::maive_infer_coef(model, 1L, 3L, data, "g", "CR0")
  options(MAIVE.bootstrap.seed = 7L)
  inf_7_again <- MAIVE:::maive_infer_coef(model, 1L, 3L, data, "g", "CR0")

  expect_identical(inf_7$boot_result$boot_coefs, inf_7_again$boot_result$boot_coefs)
  expect_false(identical(inf_7$boot_result$boot_coefs, inf_8$boot_result$boot_coefs))
})

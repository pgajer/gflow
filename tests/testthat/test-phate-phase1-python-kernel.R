test_that("PHATE Phase 1 python kernel mode matches reference K and P", {
  fixture.dir <- testthat::test_path("fixtures", "phate_phase1")
  X <- as.matrix(utils::read.csv(file.path(fixture.dir, "dense_36x5_X.csv"),
                                 header = FALSE))
  K.ref <- unname(as.matrix(utils::read.csv(file.path(fixture.dir, "dense_36x5_K_python.csv"),
                                            header = FALSE)))
  P.ref <- unname(as.matrix(utils::read.csv(file.path(fixture.dir, "dense_36x5_P_python.csv"),
                                            header = FALSE)))

  fit <- phate.core(
    X = X,
    k = 5,
    alpha.decay = 40,
    t = 3,
    kernel.mode = "python",
    compute.D.pot = FALSE,
    verbose = FALSE
  )

  expect_equal(fit$diagnostics$kernel_mode, "python")
  expect_equal(fit$diagnostics$kernel_threshold, 1e-4)
  expect_equal(diag(fit$K), rep(1, nrow(X)), tolerance = 1e-12)
  expect_equal(rowSums(fit$P), rep(1, nrow(X)), tolerance = 1e-12)
  expect_equal(mean(fit$K == 0), mean(K.ref == 0), tolerance = 1e-12)

  expect_equal(fit$K, K.ref, tolerance = 1e-6)
  expect_equal(fit$P, P.ref, tolerance = 1e-6)
})


test_that("PHATE default gflow kernel mode preserves legacy diagonal behavior", {
  X <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1,
    2, 2
  ), ncol = 2, byrow = TRUE)

  fit <- phate.core(
    X = X,
    k = 2,
    t = 2,
    kernel.mode = "gflow",
    compute.D.pot = FALSE
  )

  expect_equal(fit$diagnostics$kernel_mode, "gflow")
  expect_equal(diag(fit$K), rep(0, nrow(X)), tolerance = 1e-12)
  expect_equal(rowSums(fit$P), rep(1, nrow(X)), tolerance = 1e-12)
})

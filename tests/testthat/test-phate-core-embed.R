test_that("phate.core works across input modes and preserves key invariants", {
  set.seed(101)
  X <- matrix(rnorm(90), ncol = 3)
  n <- nrow(X)

  # X input with auto t
  fit.x <- phate.core(
    X = X,
    k = 8,
    t = "auto",
    t.max = 12,
    compute.D.pot = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit.x, "phate_core")
  expect_equal(fit.x$input_type, "X")
  expect_equal(dim(fit.x$P), c(n, n))
  expect_equal(dim(fit.x$Pt), c(n, n))
  expect_equal(dim(fit.x$U.pot), c(n, n))
  expect_equal(dim(fit.x$D.pot), c(n, n))
  expect_true(all(is.finite(fit.x$U.pot)))
  expect_true(all(fit.x$P >= 0))
  expect_equal(rowSums(fit.x$P), rep(1, n), tolerance = 1e-8)
  expect_true(fit.x$t >= 1 && fit.x$t <= 12)
  expect_true(fit.x$diagnostics$t_auto)
  expect_true(is.numeric(fit.x$diagnostics$vne))
  expect_equal(length(fit.x$diagnostics$vne), 12)

  # D input with fixed t
  D <- as.matrix(dist(X))
  fit.d <- phate.core(
    D = D,
    k = 7,
    t = 3,
    compute.D.pot = FALSE,
    verbose = FALSE
  )
  expect_equal(fit.d$input_type, "D")
  expect_equal(fit.d$t, 3)
  expect_null(fit.d$D.pot)
  expect_equal(dim(fit.d$P), c(n, n))
  expect_equal(rowSums(fit.d$P), rep(1, n), tolerance = 1e-8)

  # K input
  K <- fit.x$K
  fit.k <- phate.core(
    K = K,
    t = 2,
    compute.D.pot = FALSE,
    verbose = FALSE
  )
  expect_equal(fit.k$input_type, "K")
  expect_equal(fit.k$t, 2)
  expect_equal(dim(fit.k$P), c(n, n))
  expect_equal(rowSums(fit.k$P), rep(1, n), tolerance = 1e-8)

  # P input
  fit.p <- phate.core(
    P = fit.k$P,
    t = 2,
    compute.D.pot = FALSE,
    verbose = FALSE
  )
  expect_equal(fit.p$input_type, "P")
  expect_equal(fit.p$t, 2)
  expect_null(fit.p$K)
  expect_equal(dim(fit.p$P), c(n, n))
  expect_equal(dim(fit.p$U.pot), c(n, n))
  expect_equal(rowSums(fit.p$P), rep(1, n), tolerance = 1e-8)
})


test_that("phate.embed returns expected shape and diagnostics (metric MDS)", {
  set.seed(202)
  X <- matrix(rnorm(120), ncol = 4)
  core <- phate.core(
    X = X,
    k = 9,
    t = "auto",
    t.max = 10,
    compute.D.pot = FALSE,
    verbose = FALSE
  )

  emb <- phate.embed(
    core = core,
    ndim = 2,
    method = "metric_mds",
    verbose = FALSE
  )

  expect_s3_class(emb, "phate_embedding")
  expect_equal(emb$method, "metric_mds")
  expect_equal(nrow(emb$embedding), nrow(X))
  expect_equal(ncol(emb$embedding), 2)
  expect_true(all(is.finite(emb$embedding)))
  expect_true(is.finite(emb$stress))
  expect_true(emb$stress >= 0)
  expect_true(is.finite(emb$cor_spearman))
  expect_true(emb$cor_spearman >= -1 && emb$cor_spearman <= 1)
})


test_that("phate.embed supports non-metric MDS and direct D.pot input", {
  skip_if_not_installed("MASS")

  set.seed(303)
  X <- matrix(rnorm(105), ncol = 3)
  core <- phate.core(
    X = X,
    k = 8,
    t = 3,
    compute.D.pot = TRUE,
    verbose = FALSE
  )

  emb <- phate.embed(
    D.pot = core$D.pot,
    ndim = 2,
    method = "nonmetric_mds",
    maxit = 40,
    tol = 1e-3,
    verbose = FALSE
  )

  expect_s3_class(emb, "phate_embedding")
  expect_equal(emb$method, "nonmetric_mds")
  expect_equal(nrow(emb$embedding), nrow(X))
  expect_equal(ncol(emb$embedding), 2)
  expect_true(is.finite(emb$stress))
  expect_true(emb$stress >= 0)
  expect_true(is.finite(emb$stress_raw))
})


test_that("phate one-call wrapper returns coherent core/embed outputs", {
  set.seed(404)
  X <- matrix(rnorm(160), ncol = 4)

  fit <- phate(
    X = X,
    k = 10,
    t = "auto",
    t.max = 12,
    ndim = 2,
    embed.method = "metric_mds",
    compute.D.pot = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "phate")
  expect_s3_class(fit$core, "phate_core")
  expect_s3_class(fit$embed, "phate_embedding")
  expect_equal(fit$t, fit$core$t)
  expect_equal(nrow(fit$embedding), nrow(X))
  expect_equal(ncol(fit$embedding), 2)
  expect_equal(fit$embed$method, "metric_mds")
  expect_true(all(is.finite(fit$embedding)))
})

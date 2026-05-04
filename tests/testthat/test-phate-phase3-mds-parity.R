.phase3_fixture_D <- function() {
  phase2.dir <- testthat::test_path("fixtures", "phate_phase2")
  unname(as.matrix(utils::read.csv(
    file.path(phase2.dir, "dense_36x5_Dpot_t3_gamma1_python.csv"),
    header = FALSE
  )))
}


.phase3_procrustes_rmse <- function(X, Y) {
  X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  Y <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
  sv <- svd(t(Y) %*% X)
  Y.rot <- Y %*% sv$u %*% t(sv$v)
  scale.factor <- sum(sv$d) / sum(Y^2)
  sqrt(mean((X - scale.factor * Y.rot)^2))
}


test_that("PHATE Phase 3 classic MDS matches stats::cmdscale geometry", {
  D <- .phase3_fixture_D()

  emb <- phate.embed(
    D.pot = D,
    ndim = 2,
    method = "classic",
    verbose = FALSE
  )

  ref <- stats::cmdscale(stats::as.dist(D), k = 2, eig = TRUE, add = TRUE)
  Y.ref <- as.matrix(ref$points)

  expect_equal(emb$method, "classic")
  expect_equal(emb$diagnostics$backend, "stats::cmdscale")
  expect_lt(.phase3_procrustes_rmse(emb$embedding, Y.ref), 1e-10)
  expect_equal(as.matrix(stats::dist(emb$embedding)),
               as.matrix(stats::dist(Y.ref)),
               tolerance = 1e-10)
  expect_true(is.finite(emb$stress))
})


test_that("PHATE Phase 3 metric SMACOF matches smacof::smacofSym ratio backend", {
  skip_if_not_installed("smacof")

  D <- .phase3_fixture_D()
  classic <- phate.embed(D.pot = D, ndim = 2, method = "classic", verbose = FALSE)

  emb <- phate.embed(
    D.pot = D,
    ndim = 2,
    method = "metric",
    maxit = 200,
    tol = 1e-8,
    verbose = FALSE
  )

  ref <- smacof::smacofSym(
    delta = stats::as.dist(D),
    ndim = 2,
    type = "ratio",
    init = classic$embedding,
    itmax = 200,
    eps = 1e-8,
    verbose = FALSE
  )

  expect_equal(emb$method, "metric")
  expect_equal(emb$diagnostics$backend, "smacof::smacofSym")
  expect_equal(emb$diagnostics$smacof_type, "ratio")
  expect_equal(emb$stress_raw, ref$stress, tolerance = 1e-10)
  expect_lt(.phase3_procrustes_rmse(emb$embedding, ref$conf), 1e-8)
  expect_lte(emb$stress_raw, classic$stress + 1e-6)
})


test_that("PHATE Phase 3 nonmetric SMACOF matches smacof::smacofSym ordinal backend", {
  skip_if_not_installed("smacof")

  D <- .phase3_fixture_D()
  classic <- phate.embed(D.pot = D, ndim = 2, method = "classic", verbose = FALSE)

  emb <- phate.embed(
    D.pot = D,
    ndim = 2,
    method = "nonmetric",
    maxit = 200,
    tol = 1e-8,
    verbose = FALSE
  )

  ref <- smacof::smacofSym(
    delta = stats::as.dist(D),
    ndim = 2,
    type = "ordinal",
    init = classic$embedding,
    itmax = 200,
    eps = 1e-8,
    verbose = FALSE
  )

  expect_equal(emb$method, "nonmetric")
  expect_equal(emb$diagnostics$backend, "smacof::smacofSym")
  expect_equal(emb$diagnostics$smacof_type, "ordinal")
  expect_equal(emb$stress_raw, ref$stress, tolerance = 1e-10)
  expect_lt(.phase3_procrustes_rmse(emb$embedding, ref$conf), 1e-8)
  expect_gte(emb$cor_spearman, classic$cor_spearman - 1e-6)
})


test_that("PHATE Phase 3 deprecated MDS method aliases map to canonical names", {
  skip_if_not_installed("smacof")

  D <- .phase3_fixture_D()

  expect_warning(
    metric <- phate.embed(D.pot = D, ndim = 2, method = "metric_mds",
                          maxit = 20, verbose = FALSE),
    "deprecated"
  )
  expect_warning(
    nonmetric <- phate.embed(D.pot = D, ndim = 2, method = "nonmetric_mds",
                             maxit = 20, verbose = FALSE),
    "deprecated"
  )

  expect_equal(metric$method, "metric")
  expect_equal(metric$method_requested, "metric_mds")
  expect_equal(nonmetric$method, "nonmetric")
  expect_equal(nonmetric$method_requested, "nonmetric_mds")
})

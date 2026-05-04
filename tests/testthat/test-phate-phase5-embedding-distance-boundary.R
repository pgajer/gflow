.phase5_fixture <- function() {
  phase2.dir <- testthat::test_path("fixtures", "phate_phase2")
  list(
    U = unname(as.matrix(utils::read.csv(
      file.path(phase2.dir, "dense_36x5_U_t3_gamma1_python.csv"),
      header = FALSE
    ))),
    D = unname(as.matrix(utils::read.csv(
      file.path(phase2.dir, "dense_36x5_Dpot_t3_gamma1_python.csv"),
      header = FALSE
    )))
  )
}


test_that("PHATE Phase 5 labels potential embedding distance boundary", {
  fixture <- .phase5_fixture()

  emb <- phate.embed(
    D.pot = fixture$D,
    ndim = 2,
    method = "classic",
    verbose = FALSE
  )

  expect_equal(emb$diagnostics$distance_source, "D.pot")
  expect_equal(emb$diagnostics$embedding_distance_source, "potential")
  expect_equal(emb$diagnostics$embedding_distance_input, "D.pot")
  expect_equal(emb$diagnostics$embedding_distance_backend, "provided_distance")
  expect_equal(emb$diagnostics$mds_backend, "stats::cmdscale")
  expect_equal(emb$diagnostics$backend, "stats::cmdscale")
})


test_that("PHATE Phase 5 U.pot and D.pot routes preserve classic embedding geometry", {
  fixture <- .phase5_fixture()

  emb.d <- phate.embed(
    D.pot = fixture$D,
    ndim = 2,
    method = "classic",
    verbose = FALSE
  )
  emb.u <- phate.embed(
    U.pot = fixture$U,
    ndim = 2,
    method = "classic",
    verbose = FALSE
  )

  expect_equal(emb.u$diagnostics$embedding_distance_source, "potential")
  expect_equal(emb.u$diagnostics$embedding_distance_input, "U.pot")
  expect_equal(emb.u$diagnostics$embedding_distance_backend, "euclidean_rows")
  expect_equal(as.matrix(stats::dist(emb.u$embedding)),
               as.matrix(stats::dist(emb.d$embedding)),
               tolerance = 1e-8)
  expect_equal(emb.u$stress, emb.d$stress, tolerance = 1e-7)
})


test_that("PHATE Phase 5 MDS backend boundary preserves metric and nonmetric outputs", {
  skip_if_not_installed("smacof")

  fixture <- .phase5_fixture()

  metric <- phate.embed(
    D.pot = fixture$D,
    ndim = 2,
    method = "metric",
    maxit = 100,
    tol = 1e-7,
    verbose = FALSE
  )
  nonmetric <- phate.embed(
    D.pot = fixture$D,
    ndim = 2,
    method = "nonmetric",
    maxit = 100,
    tol = 1e-7,
    verbose = FALSE
  )

  expect_equal(metric$diagnostics$embedding_distance_source, "potential")
  expect_equal(metric$diagnostics$mds_backend, "smacof::smacofSym")
  expect_equal(metric$diagnostics$backend, "smacof::smacofSym")
  expect_equal(metric$diagnostics$smacof_type, "ratio")
  expect_true(all(is.finite(metric$embedding)))
  expect_true(is.finite(metric$stress_raw))

  expect_equal(nonmetric$diagnostics$embedding_distance_source, "potential")
  expect_equal(nonmetric$diagnostics$mds_backend, "smacof::smacofSym")
  expect_equal(nonmetric$diagnostics$backend, "smacof::smacofSym")
  expect_equal(nonmetric$diagnostics$smacof_type, "ordinal")
  expect_true(all(is.finite(nonmetric$embedding)))
  expect_true(is.finite(nonmetric$stress_raw))
})

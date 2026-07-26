test_that("archived fitted-model analysis APIs are not exported by gflow", {
  retired <- c(
    "lslope.test",
    "significant.vertices",
    "subject.neighborhood.stats",
    "gflow.smooth.spline"
  )
  expect_false(any(retired %in% getNamespaceExports("gflow")))
})

test_that("posterior local correlation consumes supplied field draws", {
  adj <- list(2L, c(1L, 3L), 2L)
  weights <- lapply(adj, function(x) rep(1, length(x)))
  draws <- cbind(
    c(0, 1, 2),
    c(0.1, 0.9, 2.1),
    c(-0.1, 1.1, 1.9)
  )

  result <- lcor.with.posterior(
    adj.list = adj,
    weight.list = weights,
    y.hat = c(0, 1, 2),
    Z.hat.samples = draws,
    verbose = FALSE
  )

  expect_s3_class(result, "lcor.posterior")
  expect_identical(result$mode, "R")
  expect_identical(result$n.samples, 3L)
  expect_identical(result$n.features, 1L)
})

test_that("protected basin and extrema hooks remain explicitly out of scope", {
  protected.files <- testthat::test_path(
    "..", "..", "R",
    c(
      "graph_endpoint_geometry.R",
      "extremality_utils.R",
      "compute_pextrema_nbhds.R",
      "compute_gfc_modulation.R",
      "vertex_density.R"
    )
  )
  expect_true(all(file.exists(protected.files)))
  suggests <- utils::packageDescription("gflow", fields = "Suggests")
  expect_match(suggests, "(^|,\\s*)gflowx(\\s|,|$)")
})

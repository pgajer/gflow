test_that("canonical local and flow-aware association APIs are exported", {
  exports <- getNamespaceExports("gflow")
  canonical <- c(
    "lcor",
    "lslope",
    "lslope.neighborhood",
    "lcor.with.posterior",
    "permutation.test.lcor",
    "gfcor",
    "gfassoc.membership",
    "gfassoc.polarity",
    "gfassoc.overlap",
    "gfassoc.deviation"
  )
  expect_true(all(canonical %in% exports))
})

test_that("shape-specific local correlation helpers are private", {
  exports <- getNamespaceExports("gflow")
  private <- c(
    "lcor.vector.matrix",
    "lcor.matrix.matrix"
  )
  expect_false(any(private %in% exports))
})

test_that("archived functional association implementations are absent", {
  archived <- c(
    "fassoc.test",
    "fassoc0.test",
    "fassoc1.test",
    "zofam.test",
    "fofam.test",
    "functional.association.test",
    "total.associations",
    "create.delta1.Delta1.df",
    "delta.indices",
    "compute.diff.profiles",
    ".coef.assoc0", ".coef.assoc1",
    ".print.assoc0", ".print.assoc1", ".print.fassoc",
    ".summary.assoc0", ".summary.assoc1", ".summary.fassoc",
    ".print.summary.assoc0", ".print.summary.assoc1",
    ".plot.assoc0", ".plot.assoc1",
    ".extract.derivatives", ".extract.assoc1.derivatives"
  )
  expect_false(any(vapply(archived, exists, logical(1),
                         envir = asNamespace("gflow"), inherits = FALSE)))
})

test_that("canonical local front doors preserve vertex grain", {
  adj <- list(2L, c(1L, 3L), 2L)
  weights <- lapply(adj, function(x) rep(1, length(x)))
  y <- c(1, 2, 4)
  z <- c(2, 3, 7)

  cor.result <- lcor(
    adj, weights, y, z,
    y.diff.type = "difference",
    z.diff.type = "difference"
  )
  slope.result <- lslope(
    adj, weights, y, z,
    y.diff.type = "difference",
    z.diff.type = "difference"
  )

  expect_length(cor.result, length(adj))
  expect_s3_class(cor.result, "lcor_vector_vector_result")
  expect_length(slope.result, length(adj))
  expect_type(slope.result, "double")
})

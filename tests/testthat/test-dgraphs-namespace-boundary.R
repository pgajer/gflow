test_that("gflow and dgraphs overlap only on protected APIs", {
  protected.overlap <- c(
    "detect.graph.endpoints",
    "detect.local.extrema",
    "vertices"
  )

  overlap <- intersect(
    getNamespaceExports("gflow"),
    getNamespaceExports("dgraphs")
  )

  expect_setequal(overlap, protected.overlap)
})

test_that("removed graph compatibility wrappers are absent", {
  removed <- c(
    "as_igraph",
    "graph.connected.components",
    "geodesic.core.endpoints",
    "graph.geodesic.distances",
    "shortest.path",
    "summarize.isometry.deviation",
    "build.iknn.graphs.and.selectk",
    "create.adaptive.radius.graph",
    "create.cknn.graph",
    "create.cmst.graph",
    "create.geodesic.iknn.graph",
    "create.iknn.graphs",
    "create.iterated.iknn.graphs",
    "create.mknn.graph",
    "create.mknn.graphs",
    "create.radius.graph",
    "create.single.iknn.graph",
    "create.sknn.graph",
    "compute.stability.metrics"
  )

  expect_false(any(removed %in% getNamespaceExports("gflow")))

  still.defined <- vapply(
    removed,
    exists,
    logical(1),
    envir = asNamespace("gflow"),
    inherits = FALSE
  )
  expect_false(any(still.defined))
})

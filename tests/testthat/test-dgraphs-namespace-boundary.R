test_that("gflow and dgraphs share only the vertices generic", {
  overlap <- intersect(
    getNamespaceExports("gflow"),
    getNamespaceExports("dgraphs")
  )

  expect_setequal(overlap, "vertices")
  expect_identical(gflow::vertices, dgraphs::vertices)
  expect_identical(
    getS3method("vertices", "feature.carriers", envir = asNamespace("dgraphs")),
    getFromNamespace("vertices.feature.carriers", "gflow")
  )
  expect_false(exists("detect.local.extrema", asNamespace("gflow"),
                      inherits = FALSE))
})

test_that("gflow does not register methods for dgraphs extrema classes", {
  registered <- getNamespaceInfo("gflow", "S3methods")
  expect_false(any(registered[, 2L] %in%
                     c("local_extrema", "summary.local_extrema")))
  expect_identical(
    getS3method("summary", "local_extrema"),
    getFromNamespace("summary.local_extrema", "dgraphs")
  )
  expect_identical(
    getS3method("print", "summary.local_extrema"),
    getFromNamespace("print.summary.local_extrema", "dgraphs")
  )
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

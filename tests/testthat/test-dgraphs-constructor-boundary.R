test_that("graph constructors are not exposed by gflow", {
  migrated <- c(
    "create.adaptive.radius.graph",
    "create.bi.kNN.chain.graph",
    "create.bipartite.graph",
    "create.chain.graph",
    "create.chain.graph.with.offset",
    "create.circular.graph",
    "create.cknn.graph",
    "create.cmst.graph",
    "create.complete.graph",
    "create.empty.graph",
    "create.geodesic.iknn.graph",
    "create.grid.graph",
    "create.iknn.graphs",
    "create.iterated.iknn.graphs",
    "create.maximal.packing",
    "create.mknn.graph",
    "create.mknn.graphs",
    "create.path.graph",
    "create.path.graph.series",
    "create.plm.graph",
    "create.radius.graph",
    "create.random.graph",
    "create.rknn.graph",
    "create.single.iknn.graph",
    "create.sknn.graph",
    "create.star.graph",
    "create.subgraph",
    "create.threshold.distance.graph",
    "generate.circle.graph",
    "join.graphs",
    "nerve.graph",
    "validate.maximal.packing",
    "verify.maximal.packing"
  )

  gflow.exports <- getNamespaceExports("gflow")
  expect_false(any(migrated %in% gflow.exports))

  still.defined <- vapply(
    migrated,
    exists,
    logical(1),
    envir = asNamespace("gflow"),
    inherits = FALSE
  )
  expect_false(any(still.defined))
})

test_that("migrated graph constructors are supplied by dgraphs", {
  expect_true(requireNamespace("dgraphs", quietly = TRUE))

  migrated <- c(
    "create.adaptive.radius.graph",
    "create.bi.kNN.chain.graph",
    "create.bipartite.graph",
    "create.chain.graph",
    "create.chain.graph.with.offset",
    "create.circular.graph",
    "create.cknn.graph",
    "create.cmst.graph",
    "create.complete.graph",
    "create.empty.graph",
    "create.geodesic.iknn.graph",
    "create.grid.graph",
    "create.iknn.graphs",
    "create.iterated.iknn.graphs",
    "create.maximal.packing",
    "create.mknn.graph",
    "create.mknn.graphs",
    "create.path.graph",
    "create.path.graph.series",
    "create.plm.graph",
    "create.radius.graph",
    "create.random.graph",
    "create.rknn.graph",
    "create.single.iknn.graph",
    "create.sknn.graph",
    "create.star.graph",
    "create.subgraph",
    "create.threshold.distance.graph",
    "generate.circle.graph",
    "join.graphs",
    "nerve.graph",
    "validate.maximal.packing",
    "verify.maximal.packing"
  )

  expect_true(all(migrated %in% getNamespaceExports("dgraphs")))
})

test_that("unmigrated public graph constructors do not remain in gflow", {
  retired <- c("create.hHN.graph", "create.intersection.graph")

  expect_false(any(retired %in% getNamespaceExports("gflow")))

  still.defined <- vapply(
    retired,
    exists,
    logical(1),
    envir = asNamespace("gflow"),
    inherits = FALSE
  )
  expect_false(any(still.defined))
})

test_that("retired graph constructor native routines are not registered", {
  dll <- getLoadedDLLs()[["gflow"]]
  expect_false(is.null(dll))

  call.names <- names(getDLLRegisteredRoutines(dll)[[".Call"]])
  retired <- c(
    "S_adaptive_radius_edges_ann",
    "S_compare_iknn_graph_builders",
    "S_create_geodesic_iknn_graph",
    "S_create_iknn_graphs",
    "S_create_maximal_packing",
    "S_create_mknn_graph",
    "S_create_mknn_graphs",
    "S_create_mst_completion_graph",
    "S_create_path_graph_plm",
    "S_create_path_graph_plus",
    "S_create_path_graph_series",
    "S_create_single_iknn_graph",
    "S_create_sknn_graph",
    "S_create_uniform_grid_graph",
    "S_join_graphs"
  )

  expect_false(any(retired %in% call.names))

  temporary.internal <- c("S_create_hHN_graph")
  expect_true(all(temporary.internal %in% call.names))
})

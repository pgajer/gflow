.sknn_edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}


test_that("create.sknn.graph constructs union kNN support", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(3, 0),
    c(7, 0)
  )

  g <- create.sknn.graph(X, k = 1)

  expect_s3_class(g, "sknn_graph")
  expect_equal(g$k, 1L)
  expect_equal(.sknn_edge_keys(g$edge_matrix), c("1-2", "2-3", "3-4"))
  expect_equal(g$edge_weight, c(1, 2, 4))
  expect_equal(g$n_components_before, 1L)
  expect_equal(g$n_components_after, 1L)
  expect_equal(g$n_mst_edges_added, 0L)
})


test_that("ANN neighbor backend matches exact sKNN on general-position data", {
  X <- rbind(
    c(0.0, 0.0),
    c(1.1, 0.2),
    c(2.4, 0.1),
    c(4.0, 1.5),
    c(6.2, 0.3),
    c(7.7, 1.9)
  )

  exact <- create.sknn.graph(X, k = 2, neighbor.method = "exact")
  ann <- create.sknn.graph(X, k = 2, neighbor.method = "ann", ann.eps = 0)

  expect_equal(ann$neighbor_method, "ann")
  expect_equal(ann$ann_eps, 0)
  expect_equal(.sknn_edge_keys(ann$edge_matrix), .sknn_edge_keys(exact$edge_matrix))
  expect_equal(ann$edge_weight, exact$edge_weight, tolerance = 1e-12)
  expect_equal(ann$knn_index, exact$knn_index)
})


test_that("ANN backend returns Euclidean edge lengths rather than squared ANN distances", {
  X <- rbind(
    c(0, 0),
    c(3, 4),
    c(10, 0)
  )

  g <- create.sknn.graph(X, k = 1, neighbor.method = "ann")

  expect_equal(.sknn_edge_keys(g$edge_matrix), c("1-2", "2-3"))
  expect_equal(g$edge_weight, c(5, sqrt(65)), tolerance = 1e-12)
})


test_that("component.mst adds the minimum number of inter-component bridges", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(30, 0),
    c(31, 0)
  )

  g <- create.sknn.graph(X, k = 1, connect.components = TRUE)

  expect_equal(g$n_components_before, 3L)
  expect_equal(g$n_components_after, 1L)
  expect_equal(g$n_mst_edges_added, 2L)
  expect_equal(.sknn_edge_keys(g$mst_edge_matrix), c("2-3", "4-5"))
  expect_equal(g$mst_edge_weight, c(9, 19))
  expect_equal(.sknn_edge_keys(g$edge_matrix), c("1-2", "2-3", "3-4", "4-5", "5-6"))
})


test_that("ANN backend keeps component.mst bridge search exact", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(30, 0),
    c(31, 0)
  )

  exact <- create.sknn.graph(X, k = 1, connect.components = TRUE)
  ann <- create.sknn.graph(X, k = 1, connect.components = TRUE,
                           neighbor.method = "ann")

  expect_equal(ann$n_components_before, exact$n_components_before)
  expect_equal(ann$n_components_after, 1L)
  expect_equal(.sknn_edge_keys(ann$mst_edge_matrix), c("2-3", "4-5"))
  expect_equal(ann$mst_edge_weight, c(9, 19), tolerance = 1e-12)
  expect_equal(.sknn_edge_keys(ann$edge_matrix), .sknn_edge_keys(exact$edge_matrix))
})


test_that("component.mst.ann uses sparse ANN bridges when they connect components", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(30, 0),
    c(31, 0)
  )

  exact <- create.sknn.graph(X, k = 1, connect.components = TRUE,
                             connect.method = "component.mst")
  ann.bridge <- create.sknn.graph(
    X, k = 1, connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 2,
    bridge.k.max = 2
  )

  expect_equal(ann.bridge$connect_method, "component.mst.ann")
  expect_equal(ann.bridge$bridge_method, "ann")
  expect_equal(ann.bridge$bridge_k_used, 2L)
  expect_false(ann.bridge$bridge_exact_fallback_used)
  expect_equal(ann.bridge$n_components_after, 1L)
  expect_equal(.sknn_edge_keys(ann.bridge$mst_edge_matrix), .sknn_edge_keys(exact$mst_edge_matrix))
  expect_equal(ann.bridge$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
})


test_that("component.mst.ann falls back to exact bridges when sparse candidates do not connect", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(30, 0),
    c(31, 0)
  )

  exact <- create.sknn.graph(X, k = 1, connect.components = TRUE,
                             connect.method = "component.mst")
  fallback <- create.sknn.graph(
    X, k = 1, connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 1,
    bridge.k.max = 1
  )

  expect_equal(fallback$bridge_method, "ann_then_exact")
  expect_true(fallback$bridge_exact_fallback_used)
  expect_equal(fallback$n_components_after, 1L)
  expect_equal(.sknn_edge_keys(fallback$mst_edge_matrix), .sknn_edge_keys(exact$mst_edge_matrix))
  expect_equal(fallback$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
})


test_that("global.mst unions sKNN with the full Euclidean MST", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(30, 0),
    c(31, 0)
  )

  component <- create.sknn.graph(X, k = 1, connect.components = TRUE,
                                 connect.method = "component.mst")
  global <- create.sknn.graph(X, k = 1, connect.components = TRUE,
                              connect.method = "global.mst")

  expect_equal(global$n_components_after, 1L)
  expect_equal(.sknn_edge_keys(global$edge_matrix), .sknn_edge_keys(component$edge_matrix))
  expect_equal(.sknn_edge_keys(global$mst_edge_matrix), c("2-3", "4-5"))
})


test_that("local geometric pruning removes redundant long sKNN edges", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  raw <- create.sknn.graph(X, k = 2)
  pruned <- create.sknn.graph(
    X,
    k = 2,
    prune.edges = TRUE,
    prune.tau = 1.01,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.sknn_edge_keys(raw$edge_matrix), c("1-2", "1-3", "2-3", "2-4", "3-4"))
  expect_equal(.sknn_edge_keys(pruned$edge_matrix), c("1-2", "2-3", "3-4"))
  expect_equal(pruned$n_edges_before_pruning, 5L)
  expect_equal(pruned$n_edges_after_pruning, 3L)
  expect_equal(pruned$n_pruned_edges, 2L)
  expect_equal(.sknn_edge_keys(as.matrix(pruned$pruned_edge_stats[, c("u", "v")])),
               c("1-3", "2-4"))
  expect_equal(pruned$pruned_edge_stats$edge_length, c(2, 2), tolerance = 1e-12)
  expect_equal(pruned$pruned_edge_stats$alt_path_length, c(2, 2), tolerance = 1e-12)
  expect_equal(pruned$pruned_edge_stats$path_edge_ratio, c(1, 1), tolerance = 1e-12)
  expect_equal(pruned$n_components_before, raw$n_components_before)
  expect_equal(pruned$n_components_after, raw$n_components_after)
})


test_that("local geometric pruning is opt-in and can omit edge stats", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0)
  )

  raw <- create.sknn.graph(X, k = 2)
  pruned <- create.sknn.graph(X, k = 2, prune.edges = TRUE, prune.tau = 1.01)

  expect_equal(raw$n_pruned_edges, 0L)
  expect_false(raw$prune_edges)
  expect_equal(raw$n_edges_before_pruning, raw$n_edges)
  expect_equal(raw$n_edges_after_pruning, raw$n_edges)
  expect_true(is.data.frame(raw$pruned_edge_stats))
  expect_equal(nrow(raw$pruned_edge_stats), 0L)

  expect_true(pruned$prune_edges)
  expect_equal(pruned$n_pruned_edges, 1L)
  expect_equal(nrow(pruned$pruned_edge_stats), 0L)
  expect_equal(.sknn_edge_keys(pruned$edge_matrix), c("1-2", "2-3"))
})


test_that("local geometric pruning preserves components before MST repair", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(10, 0),
    c(11, 0),
    c(12, 0)
  )

  pruned <- create.sknn.graph(X, k = 2, prune.edges = TRUE, prune.tau = 1.01)
  connected <- create.sknn.graph(
    X,
    k = 2,
    prune.edges = TRUE,
    prune.tau = 1.01,
    connect.components = TRUE
  )

  expect_equal(pruned$n_components_before, 2L)
  expect_equal(pruned$n_components_after, 2L)
  expect_equal(pruned$n_pruned_edges, 2L)
  expect_equal(.sknn_edge_keys(pruned$edge_matrix), c("1-2", "2-3", "4-5", "5-6"))

  expect_equal(connected$n_components_before, 2L)
  expect_equal(connected$n_components_after, 1L)
  expect_equal(connected$n_mst_edges_added, 1L)
  expect_equal(.sknn_edge_keys(connected$mst_edge_matrix), "3-4")
  expect_equal(.sknn_edge_keys(connected$edge_matrix),
               c("1-2", "2-3", "3-4", "4-5", "5-6"))
})


test_that("local geometric pruning works with ANN neighbor backend", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0),
    c(4, 0)
  )

  exact <- create.sknn.graph(X, k = 2, prune.edges = TRUE, prune.tau = 1.01)
  ann <- create.sknn.graph(
    X,
    k = 2,
    prune.edges = TRUE,
    prune.tau = 1.01,
    neighbor.method = "ann"
  )

  expect_equal(ann$neighbor_method, "ann")
  expect_equal(.sknn_edge_keys(ann$edge_matrix), .sknn_edge_keys(exact$edge_matrix))
  expect_equal(ann$n_pruned_edges, exact$n_pruned_edges)
})


test_that("create.sknn.graph validates inputs", {
  X <- matrix(1:6, ncol = 2)

  expect_error(create.sknn.graph(X, k = 0), "positive integer")
  expect_error(create.sknn.graph(X, k = 3), "smaller than nrow")
  expect_error(create.sknn.graph(X, k = 1, connect.components = NA),
               "connect.components")
  expect_error(create.sknn.graph(X, k = 1, neighbor.method = "bogus"),
               "'arg' should be one of")
  expect_error(create.sknn.graph(X, k = 1, neighbor.method = "ann", ann.eps = -1),
               "ann.eps")
  expect_error(create.sknn.graph(X, k = 1, neighbor.method = "ann", ann.eps = 0.1),
               "not yet supported")
  expect_error(create.sknn.graph(X, k = 1, bridge.k = 0),
               "bridge.k")
  expect_error(create.sknn.graph(X, k = 1, bridge.k = 1, bridge.k.max = 0),
               "bridge.k.max")
  expect_error(create.sknn.graph(X, k = 1, bridge.growth = 1),
               "bridge.growth")
  expect_error(create.sknn.graph(X, k = 1, prune.edges = NA),
               "prune.edges")
  expect_error(create.sknn.graph(X, k = 1, with.pruned.edge.stats = NA),
               "with.pruned.edge.stats")
  expect_error(create.sknn.graph(X, k = 1, prune.tau = 1),
               "prune.tau")
  expect_error(create.sknn.graph(X, k = 1, prune.local.k = 0),
               "prune.local.k")
  expect_error(create.sknn.graph(c("a", "b"), k = 1), "matrix or data frame")
})

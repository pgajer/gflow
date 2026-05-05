.edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}


test_that("create.mknn.graph supports ANN component MST repair with exact fallback", {
  X <- cbind(c(0, 1, 10, 11, 30, 31), 0)

  exact <- create.mknn.graph(
    X, k = 2,
    connect.components = TRUE,
    connect.method = "component.mst"
  )
  ann <- create.mknn.graph(
    X, k = 2,
    connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 2,
    bridge.k.max = 2
  )
  fallback <- create.mknn.graph(
    X, k = 2,
    connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 1,
    bridge.k.max = 1
  )

  expect_equal(ann$bridge_method, "ann")
  expect_false(ann$bridge_exact_fallback_used)
  expect_equal(fallback$bridge_method, "ann_then_exact")
  expect_true(fallback$bridge_exact_fallback_used)
  expect_equal(ann$n_components_after, 1L)
  expect_equal(fallback$n_components_after, 1L)
  expect_equal(.edge_keys(ann$mst_edge_matrix), .edge_keys(exact$mst_edge_matrix))
  expect_equal(.edge_keys(fallback$mst_edge_matrix), .edge_keys(exact$mst_edge_matrix))
  expect_equal(ann$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
  expect_equal(fallback$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
})


test_that("create.single.iknn.graph supports ANN component MST repair with exact fallback", {
  X <- cbind(c(0, 1, 10, 11, 30, 31), 0)

  exact <- create.single.iknn.graph(
    X, k = 1,
    connect.components = TRUE,
    connect.method = "component.mst",
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )
  ann <- create.single.iknn.graph(
    X, k = 1,
    connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 2,
    bridge.k.max = 2,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )
  fallback <- create.single.iknn.graph(
    X, k = 1,
    connect.components = TRUE,
    connect.method = "component.mst.ann",
    bridge.k = 1,
    bridge.k.max = 1,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )

  expect_equal(ann$bridge_method, "ann")
  expect_false(ann$bridge_exact_fallback_used)
  expect_equal(fallback$bridge_method, "ann_then_exact")
  expect_true(fallback$bridge_exact_fallback_used)
  expect_equal(ann$n_components_after_mst, 1L)
  expect_equal(fallback$n_components_after_mst, 1L)
  expect_equal(.edge_keys(ann$mst_edge_matrix), .edge_keys(exact$mst_edge_matrix))
  expect_equal(.edge_keys(fallback$mst_edge_matrix), .edge_keys(exact$mst_edge_matrix))
  expect_equal(ann$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
  expect_equal(fallback$mst_edge_weight, exact$mst_edge_weight, tolerance = 1e-12)
  expect_equal(ann$n_edges_in_pruned_graph, sum(lengths(ann$pruned_adj_list)) / 2)
})


test_that("create.single.iknn.graph supports local geometric pruning", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  raw <- create.single.iknn.graph(
    X,
    k = 3,
    prune.method = "none",
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )
  pruned <- create.single.iknn.graph(
    X,
    k = 3,
    prune.method = "local.geodesic",
    prune.tau = 1.01,
    threshold.percentile = 0,
    with.edge.pruning.stats = TRUE,
    pca.dim = NULL,
    verbose = FALSE
  )

  expect_equal(.edge_keys(.graph.edge.table(raw$pruned_adj_list, raw$pruned_weight_list)),
               c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4"))
  expect_equal(.edge_keys(.graph.edge.table(pruned$pruned_adj_list, pruned$pruned_weight_list)),
               c("1-2", "2-3", "3-4"))
  expect_equal(pruned$prune_method, "local.geodesic")
  expect_equal(pruned$n_edges_before_pruning, 6L)
  expect_equal(pruned$n_edges_after_pruning, 3L)
  expect_equal(pruned$n_pruned_edges, 3L)
  expect_equal(pruned$n_edges_before_mst, 3L)
  expect_equal(pruned$n_edges_after_mst, 3L)
  expect_equal(.edge_keys(as.matrix(pruned$pruned_edge_stats[, c("u", "v")])),
               c("1-3", "1-4", "2-4"))
})


test_that("component MST repair validates bridge controls", {
  X <- cbind(c(0, 1, 10, 11), 0)

  expect_error(
    create.mknn.graph(X, k = 2, bridge.k = 0),
    "bridge.k"
  )
  expect_error(
    create.single.iknn.graph(
      X, k = 1, connect.components = TRUE, bridge.growth = 1,
      max.path.edge.ratio.deviation.thld = 0,
      threshold.percentile = 0,
      pca.dim = NULL,
      verbose = FALSE
    ),
    "bridge.growth"
  )
  expect_error(
    create.single.iknn.graph(
      X, k = 1, connect.components = TRUE, knn.metric = "linf.simplex",
      max.path.edge.ratio.deviation.thld = 0,
      threshold.percentile = 0,
      pca.dim = NULL,
      verbose = FALSE
    ),
    "knn.metric = 'euclidean'"
  )
})

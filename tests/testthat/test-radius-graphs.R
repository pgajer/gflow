.radius_edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}


test_that("create.radius.graph matches the fixed-radius definition", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  g <- create.radius.graph(X, radius = 1.1)

  expect_s3_class(g, "radius_graph")
  expect_equal(g$graph_rule, "fixed.radius")
  expect_equal(g$radius, 1.1)
  expect_equal(.radius_edge_keys(g$edge_matrix), "1-2")
  expect_equal(g$edge_weight, 1, tolerance = 1e-12)
  expect_equal(g$n_components_before, 2L)
  expect_equal(g$n_components_after, 2L)
})


test_that("create.radius.graph can repair components by component MST", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  g <- create.radius.graph(X, radius = 1.1, connect.components = TRUE)

  expect_equal(g$n_components_before, 2L)
  expect_equal(g$n_components_after, 1L)
  expect_equal(g$n_mst_edges_added, 1L)
  expect_equal(.radius_edge_keys(g$mst_edge_matrix), "2-3")
  expect_equal(.radius_edge_keys(g$edge_matrix), c("1-2", "2-3"))
})


test_that("create.radius.graph supports local geometric pruning", {
  X <- matrix(0:3, ncol = 1)

  g <- create.radius.graph(
    X,
    radius = 2.1,
    prune.method = "local.geodesic",
    prune.tau = 1.01,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.radius_edge_keys(.graph.edge.table(g$raw_adj_list,
                                                   g$raw_weight_list)),
               c("1-2", "1-3", "2-3", "2-4", "3-4"))
  expect_equal(.radius_edge_keys(.graph.edge.table(g$pruned_adj_list,
                                                   g$pruned_weight_list)),
               c("1-2", "2-3", "3-4"))
  expect_equal(.radius_edge_keys(g$edge_matrix), c("1-2", "2-3", "3-4"))
  expect_equal(g$n_edges_before_pruning, 5L)
  expect_equal(g$n_edges_after_pruning, 3L)
  expect_equal(g$n_pruned_edges, 2L)
  expect_equal(g$prune_method, "local.geodesic")
  expect_equal(g$prune_tau, 1.01)
  expect_equal(g$prune_local_k, 2L)
  expect_equal(nrow(g$pruned_edge_stats), 2L)
  expect_equal(g$pruned_edge_stats$path_edge_ratio, c(1, 1), tolerance = 1e-12)
})


test_that("create.adaptive.radius.graph distinguishes max, min, and geomean rules", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  max.g <- create.adaptive.radius.graph(
    X, k.scale = 1, radius.rule = "max", radius.factor = 1.5
  )
  min.g <- create.adaptive.radius.graph(
    X, k.scale = 1, radius.rule = "min", radius.factor = 1.5
  )
  geomean.g <- create.adaptive.radius.graph(
    X, k.scale = 1, radius.rule = "geomean", radius.factor = 1.5
  )

  expect_s3_class(max.g, "adaptive_radius_graph")
  expect_equal(max.g$radius_search, "ann")
  expect_equal(max.g$sigma, c(1, 1, 2), tolerance = 1e-12)
  expect_equal(min.g$sigma, c(1, 1, 2), tolerance = 1e-12)
  expect_equal(geomean.g$sigma, c(1, 1, 2), tolerance = 1e-12)
  expect_equal(.radius_edge_keys(max.g$edge_matrix), c("1-2", "1-3", "2-3"))
  expect_equal(.radius_edge_keys(min.g$edge_matrix), "1-2")
  expect_equal(.radius_edge_keys(geomean.g$edge_matrix), c("1-2", "2-3"))
  expect_equal(max.g$radius_rule, "max")
  expect_equal(min.g$radius_rule, "min")
  expect_equal(geomean.g$radius_rule, "geomean")
})

test_that("ANN adaptive-radius search matches all-pairs reference", {
  set.seed(919)
  X <- cbind(runif(35), runif(35), 0.2 * runif(35)^2)

  for (rule in c("max", "min", "geomean")) {
    ann.g <- create.adaptive.radius.graph(
      X,
      k.scale = 4L,
      radius.factor = 1.35,
      radius.rule = rule,
      radius.search = "ann"
    )
    ref.g <- create.adaptive.radius.graph(
      X,
      k.scale = 4L,
      radius.factor = 1.35,
      radius.rule = rule,
      radius.search = "all.pairs"
    )

    expect_equal(ann.g$sigma, ref.g$sigma, tolerance = 1e-10)
    expect_equal(.radius_edge_keys(ann.g$edge_matrix),
                 .radius_edge_keys(ref.g$edge_matrix))
    expect_equal(ann.g$edge_weight, ref.g$edge_weight, tolerance = 1e-10)
    expect_equal(ann.g$raw_adj_list, ref.g$raw_adj_list)
    expect_equal(ann.g$radius_search, "ann")
    expect_equal(ref.g$radius_search, "all.pairs")
  }
})


test_that("create.cknn.graph is the continuous-kNN geomean adaptive-radius wrapper", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  cknn.g <- create.cknn.graph(X, k.scale = 1, delta = 1.5)
  adaptive.g <- create.adaptive.radius.graph(
    X, k.scale = 1, radius.rule = "geomean", radius.factor = 1.5
  )

  expect_s3_class(cknn.g, "cknn_graph")
  expect_s3_class(cknn.g, "adaptive_radius_graph")
  expect_equal(cknn.g$graph_rule, "continuous.knn")
  expect_equal(cknn.g$radius_rule, "geomean")
  expect_equal(cknn.g$radius_factor, 1.5)
  expect_equal(cknn.g$delta, 1.5)
  expect_equal(.radius_edge_keys(cknn.g$edge_matrix),
               .radius_edge_keys(adaptive.g$edge_matrix))
  expect_equal(cknn.g$edge_weight, adaptive.g$edge_weight, tolerance = 1e-12)
  expect_equal(graph.geodesic.distances(cknn.g),
               graph.geodesic.distances(adaptive.g),
               tolerance = 1e-12)
})


test_that("create.adaptive.radius.graph supports local geometric pruning", {
  X <- matrix(0:3, ncol = 1)

  g <- create.adaptive.radius.graph(
    X,
    k.scale = 2,
    radius.rule = "max",
    prune.method = "local.geodesic",
    prune.tau = 1.01,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.radius_edge_keys(.graph.edge.table(g$raw_adj_list,
                                                   g$raw_weight_list)),
               c("1-2", "1-3", "2-3", "2-4", "3-4"))
  expect_equal(.radius_edge_keys(.graph.edge.table(g$pruned_adj_list,
                                                   g$pruned_weight_list)),
               c("1-2", "2-3", "3-4"))
  expect_equal(.radius_edge_keys(g$edge_matrix), c("1-2", "2-3", "3-4"))
  expect_equal(g$n_edges_before_pruning, 5L)
  expect_equal(g$n_edges_after_pruning, 3L)
  expect_equal(g$n_pruned_edges, 2L)
  expect_equal(g$prune_method, "local.geodesic")
  expect_equal(g$prune_local_k, 2L)
  expect_equal(nrow(g$pruned_edge_stats), 2L)
})


test_that("radius graph constructors validate inputs", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  expect_error(create.radius.graph(X, radius = 0), "radius")
  expect_error(create.radius.graph(X, radius = Inf), "radius")
  expect_error(create.radius.graph(c("a", "b"), radius = 1), "matrix or data frame")
  expect_error(create.radius.graph(X, radius = 1, prune.method = "global"),
               "'arg' should be one of")
  expect_error(create.radius.graph(X, radius = 1, prune.tau = 1),
               "prune.tau")
  expect_error(create.radius.graph(X, radius = 1, prune.local.k = 3),
               "prune.local.k")
  expect_error(create.adaptive.radius.graph(X, k.scale = 0), "k.scale")
  expect_error(create.adaptive.radius.graph(X, k.scale = 3), "k.scale")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1, radius.factor = 0),
               "radius.factor")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1, radius.rule = "median"),
               "'arg' should be one of")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1, prune.tau = 1),
               "prune.tau")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1,
                                            with.pruned.edge.stats = NA),
               "with.pruned.edge.stats")
  expect_error(create.cknn.graph(X, k.scale = 1, delta = 0), "delta")
})

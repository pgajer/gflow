global_ratio_edges <- function(adj.list, weight.list) {
  .graph.edge.table(adj.list, weight.list)
}

global_ratio_graph <- function(n, from, to, weight) {
  .graph.from.edge.table(n, data.frame(
    from = as.integer(from),
    to = as.integer(to),
    weight = as.numeric(weight)
  ))
}

global_ratio_edge_keys <- function(edges) {
  if (!nrow(edges)) {
    return(character())
  }
  sprintf("%d-%d", edges$from, edges$to)
}

test_that("global.geodesic.ratio prunes a redundant whole-graph shortcut", {
  graph <- global_ratio_graph(
    4,
    from = c(1, 2, 3, 4, 1),
    to = c(2, 3, 4, 1, 3),
    weight = c(1, 1, 1, 1, sqrt(2))
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 1.5,
    path.edge.ratio.percentile = 0,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(
    global_ratio_edge_keys(global_ratio_edges(pruned$adj_list, pruned$weight_list)),
    c("1-2", "1-4", "2-3", "3-4")
  )
  expect_equal(pruned$n_edges_before_pruning, 5L)
  expect_equal(pruned$n_edges_after_pruning, 4L)
  expect_equal(pruned$n_pruned_edges, 1L)
  expect_equal(pruned$pruned_edge_stats$u, 1L)
  expect_equal(pruned$pruned_edge_stats$v, 3L)
  expect_equal(pruned$pruned_edge_stats$path_edge_ratio, sqrt(2), tolerance = 1e-12)
})

test_that("global.geodesic.ratio keeps edges without an alternative path", {
  graph <- global_ratio_graph(
    3,
    from = c(1, 2),
    to = c(2, 3),
    weight = c(1, 1)
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 2,
    path.edge.ratio.percentile = 0,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(
    global_ratio_edge_keys(global_ratio_edges(pruned$adj_list, pruned$weight_list)),
    c("1-2", "2-3")
  )
  expect_equal(pruned$n_pruned_edges, 0L)
  expect_s3_class(pruned$pruned_edge_stats, "data.frame")
  expect_equal(nrow(pruned$pruned_edge_stats), 0L)
})

test_that("global.geodesic.ratio respects candidate percentile filtering", {
  graph <- global_ratio_graph(
    4,
    from = c(1, 2, 3, 4, 1),
    to = c(2, 3, 4, 1, 3),
    weight = c(1, 1, 1, 1, sqrt(2))
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 1.5,
    path.edge.ratio.percentile = 1,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(
    global_ratio_edge_keys(global_ratio_edges(pruned$adj_list, pruned$weight_list)),
    c("1-2", "1-3", "1-4", "2-3", "3-4")
  )
  expect_equal(pruned$n_pruned_edges, 0L)
})

test_that("global.geodesic.ratio rechecks candidates sequentially", {
  graph <- global_ratio_graph(
    4,
    from = c(1, 2, 1, 2),
    to = c(2, 3, 4, 4),
    weight = c(1, 1, 2.4, 1.2)
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 3,
    path.edge.ratio.percentile = 0,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(
    global_ratio_edge_keys(global_ratio_edges(pruned$adj_list, pruned$weight_list)),
    c("1-2", "2-3", "2-4")
  )
  expect_equal(pruned$n_pruned_edges, 1L)
  expect_equal(pruned$pruned_edge_stats$u, 1L)
  expect_equal(pruned$pruned_edge_stats$v, 4L)
})

test_that("global.geodesic.ratio candidate ties are deterministic by from and to", {
  graph <- global_ratio_graph(
    4,
    from = c(1, 2, 3, 4, 1, 2),
    to = c(2, 3, 4, 1, 3, 4),
    weight = c(1, 1, 1, 1, sqrt(2), sqrt(2))
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 1.5,
    path.edge.ratio.percentile = 0,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(pruned$pruned_edge_stats$u, c(1L, 2L))
  expect_equal(pruned$pruned_edge_stats$v, c(3L, 4L))
})

test_that("global.geodesic.ratio stats-off keeps the existing result shape", {
  graph <- global_ratio_graph(
    4,
    from = c(1, 2, 3, 4, 1),
    to = c(2, 3, 4, 1, 3),
    weight = c(1, 1, 1, 1, sqrt(2))
  )

  pruned <- .prune.graph.global.geodesic.ratio(
    graph$adj_list, graph$weight_list,
    max.ratio.threshold = 1.5,
    path.edge.ratio.percentile = 0,
    with.pruned.edge.stats = FALSE
  )

  expect_equal(pruned$n_pruned_edges, 1L)
  expect_s3_class(pruned$pruned_edge_stats, "data.frame")
  expect_equal(nrow(pruned$pruned_edge_stats), 0L)
  expect_false(pruned$with_pruned_edge_stats)
})

test_that("global.geodesic.ratio validates malformed graph inputs", {
  graph <- global_ratio_graph(
    3,
    from = c(1, 2, 1),
    to = c(2, 3, 3),
    weight = c(1, 1, 2)
  )
  bad <- graph
  bad$adj_list[[2L]] <- bad$adj_list[[2L]][-1L]
  bad$weight_list[[2L]] <- bad$weight_list[[2L]][-1L]

  expect_error(
    .prune.graph.global.geodesic.ratio(
      bad$adj_list, bad$weight_list,
      max.ratio.threshold = 1.5,
      path.edge.ratio.percentile = 0
    ),
    "undirected graph"
  )
})

test_that("global.geodesic.ratio is routed through graph-family constructors", {
  X <- rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1),
    c(3, 0), c(4, 0), c(4, 1), c(3, 1)
  )

  graphs <- list(
    sknn = create.sknn.graph(
      X, k = 3, prune.method = "global.geodesic.ratio",
      path.edge.ratio.percentile = 0, with.pruned.edge.stats = TRUE,
      connect.components = TRUE
    ),
    mknn = create.mknn.graph(
      X, k = 3, prune.method = "global.geodesic.ratio",
      path.edge.ratio.percentile = 0, with.pruned.edge.stats = TRUE,
      connect.components = TRUE
    ),
    radius = create.radius.graph(
      X, radius = 1.5, prune.method = "global.geodesic.ratio",
      path.edge.ratio.percentile = 0, with.pruned.edge.stats = TRUE,
      connect.components = TRUE
    ),
    adaptive.radius = create.adaptive.radius.graph(
      X, k.scale = 3, prune.method = "global.geodesic.ratio",
      path.edge.ratio.percentile = 0, with.pruned.edge.stats = TRUE,
      connect.components = TRUE
    ),
    iknn = create.single.iknn.graph(
      X, k = 3, prune.method = "global.geodesic.ratio",
      path.edge.ratio.percentile = 0, with.edge.pruning.stats = TRUE,
      connect.components = TRUE, with.lifecycle.branches = TRUE,
      verbose = FALSE
    )
  )

  for (graph in graphs) {
    expect_equal(graph$prune_method, "global.geodesic.ratio")
    expect_equal(length(graph$pruned_adj_list), nrow(X))
    expect_equal(length(graph$repaired_pruned_adj_list), nrow(X))
    expect_equal(length(graph$repaired_pruned_weight_list), nrow(X))
    expect_s3_class(graph$pruned_edge_stats, "data.frame")
    expect_s3_class(graph$repaired_pruned_pruning$pruned_edge_stats, "data.frame")
  }
})

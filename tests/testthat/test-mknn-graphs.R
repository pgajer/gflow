test_that("create.mknn.graph uses neighbor-rank distances for higher-index vertices", {
  X <- matrix(c(seq_len(12), rep(0, 12)), ncol = 2)

  graph <- create.mknn.graph(X, k = 2)

  expect_s3_class(graph, "mknn_graph")
  expect_true(any(unlist(graph$adj_list, use.names = FALSE) > 4L))

  for (i in seq_along(graph$adj_list)) {
    neighbors <- graph$adj_list[[i]]
    weights <- graph$weight_list[[i]]

    expect_equal(length(weights), length(neighbors))
    if (length(neighbors) > 0) {
      expected <- abs(X[neighbors, 1] - X[i, 1])
      expect_equal(weights, expected)
    }
  }
})

.mknn_edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}

test_that("create.mknn.graph local pruning removes redundant line chords", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  raw <- create.mknn.graph(X, k = 3)
  pruned <- create.mknn.graph(
    X,
    k = 3,
    prune.method = "local.geodesic",
    prune.tau = 1.01,
    with.pruned.edge.stats = TRUE
  )
  pruned.edges <- .graph.edge.table(pruned$adj_list, pruned$weight_list)

  expect_equal(.mknn_edge_keys(.graph.edge.table(raw$adj_list, raw$weight_list)),
               c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4"))
  expect_equal(.mknn_edge_keys(pruned.edges), c("1-2", "2-3", "3-4"))
  expect_equal(pruned$n_edges_before_pruning, 6L)
  expect_equal(pruned$n_edges_after_pruning, 3L)
  expect_equal(pruned$n_pruned_edges, 3L)
  expect_equal(pruned$n_edges_before_mst, 3L)
  expect_equal(pruned$n_edges_after_mst, 3L)
  expect_equal(.mknn_edge_keys(as.matrix(pruned$pruned_edge_stats[, c("u", "v")])),
               c("1-3", "1-4", "2-4"))
  expect_equal(pruned$pruned_edge_stats$path_edge_ratio, c(1, 1, 1),
               tolerance = 1e-12)
})

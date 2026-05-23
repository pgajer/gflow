.local_cpp_edge_keys <- function(adj.list, weight.list = NULL) {
  edges <- if (is.null(weight.list)) {
    adj.list
  } else {
    .graph.edge.table(adj.list, weight.list)
  }
  if (is.null(edges) || nrow(edges) == 0L) {
    return(character(0))
  }
  sort(paste(edges[, 1L], edges[, 2L], sep = "-"))
}

.local_cpp_graph <- function(n, edges) {
  .graph.from.edge.table(n, data.frame(
    from = as.integer(edges[, 1L]),
    to = as.integer(edges[, 2L]),
    weight = as.numeric(edges[, 3L])
  ))
}

test_that("C++ local pruning removes line chords and reports diagnostics", {
  X <- matrix(0:3, ncol = 1)
  graph <- .local_cpp_graph(4L, rbind(
    c(1, 2, 1),
    c(2, 3, 1),
    c(3, 4, 1),
    c(1, 3, 2),
    c(2, 4, 2),
    c(1, 4, 3)
  ))

  pruned <- .prune.graph.local.geodesic(
    X,
    graph$adj_list,
    graph$weight_list,
    k = 3,
    prune.tau = 1.01,
    prune.local.k = 3,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.local_cpp_edge_keys(pruned$adj_list, pruned$weight_list),
               c("1-2", "2-3", "3-4"))
  expect_equal(pruned$n_edges_before_pruning, 6L)
  expect_equal(pruned$n_edges_after_pruning, 3L)
  expect_equal(pruned$n_pruned_edges, 3L)
  expect_equal(pruned$prune_tau, 1.01)
  expect_equal(pruned$prune_local_k, 3L)
  expect_true(pruned$with_pruned_edge_stats)
  expect_equal(pruned$pruned_edge_stats$u, c(1L, 1L, 2L))
  expect_equal(pruned$pruned_edge_stats$v, c(4L, 3L, 4L))
  expect_equal(pruned$pruned_edge_stats$edge_length, c(3, 2, 2),
               tolerance = 1e-12)
  expect_equal(pruned$pruned_edge_stats$alt_path_length, c(3, 2, 2),
               tolerance = 1e-12)
  expect_equal(pruned$pruned_edge_stats$path_edge_ratio, c(1, 1, 1),
               tolerance = 1e-12)
})

test_that("C++ local pruning leaves bridges and local-only alternatives alone", {
  X <- matrix(0:2, ncol = 1)
  bridge <- .local_cpp_graph(3L, rbind(
    c(1, 2, 1),
    c(2, 3, 1)
  ))

  bridge.pruned <- .prune.graph.local.geodesic(
    X, bridge$adj_list, bridge$weight_list,
    k = 1, prune.tau = 1.5, prune.local.k = 1
  )
  expect_equal(.local_cpp_edge_keys(bridge.pruned$adj_list,
                                    bridge.pruned$weight_list),
               c("1-2", "2-3"))
  expect_equal(bridge.pruned$n_pruned_edges, 0L)

  X.local <- rbind(c(0, 0), c(0.1, 0), c(10, 0), c(11, 0))
  locality <- .local_cpp_graph(4L, rbind(
    c(1, 2, 10),
    c(1, 3, 3),
    c(3, 4, 3),
    c(2, 4, 3)
  ))
  local.pruned <- .prune.graph.local.geodesic(
    X.local, locality$adj_list, locality$weight_list,
    k = 1, prune.tau = 1.05, prune.local.k = 1,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.local_cpp_edge_keys(local.pruned$adj_list,
                                    local.pruned$weight_list),
               c("1-2", "1-3", "2-4", "3-4"))
  expect_equal(nrow(local.pruned$pruned_edge_stats), 0L)
})

test_that("C++ local pruning is sequential", {
  X <- matrix(c(0, 1, 2), ncol = 1)
  graph <- .local_cpp_graph(3L, rbind(
    c(1, 2, 1),
    c(2, 3, 4),
    c(1, 3, 5)
  ))

  pruned <- .prune.graph.local.geodesic(
    X,
    graph$adj_list,
    graph$weight_list,
    k = 2,
    prune.tau = 1.5,
    prune.local.k = 2,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.local_cpp_edge_keys(pruned$adj_list, pruned$weight_list),
               c("1-2", "2-3"))
  expect_equal(pruned$pruned_edge_stats$u, 1L)
  expect_equal(pruned$pruned_edge_stats$v, 3L)
  expect_equal(pruned$pruned_edge_stats$alt_path_length, 5,
               tolerance = 1e-12)
})

test_that("C++ local pruning uses deterministic kNN tie handling", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(-1, 0),
    c(0, 0.5)
  )
  graph <- .local_cpp_graph(4L, rbind(
    c(1, 4, 10),
    c(1, 2, 5),
    c(2, 4, 5)
  ))

  pruned <- .prune.graph.local.geodesic(
    X,
    graph$adj_list,
    graph$weight_list,
    k = 2,
    prune.tau = 1.01,
    prune.local.k = 2,
    with.pruned.edge.stats = TRUE
  )

  expect_equal(.local_cpp_edge_keys(pruned$adj_list, pruned$weight_list),
               c("1-2", "2-4"))
  expect_equal(pruned$pruned_edge_stats[, c("u", "v")],
               data.frame(u = 1L, v = 4L))
})

test_that("C++ local pruning keeps stats optional without changing decisions", {
  X <- matrix(0:3, ncol = 1)
  graph <- .local_cpp_graph(4L, rbind(
    c(1, 2, 1),
    c(2, 3, 1),
    c(3, 4, 1),
    c(1, 3, 2),
    c(2, 4, 2),
    c(1, 4, 3)
  ))

  no.stats <- .prune.graph.local.geodesic(
    X, graph$adj_list, graph$weight_list,
    k = 3, prune.tau = 1.01, prune.local.k = 3,
    with.pruned.edge.stats = FALSE
  )
  with.stats <- .prune.graph.local.geodesic(
    X, graph$adj_list, graph$weight_list,
    k = 3, prune.tau = 1.01, prune.local.k = 3,
    with.pruned.edge.stats = TRUE
  )

  expect_false(no.stats$with_pruned_edge_stats)
  expect_s3_class(no.stats$pruned_edge_stats, "data.frame")
  expect_equal(nrow(no.stats$pruned_edge_stats), 0L)
  expect_equal(no.stats$n_pruned_edges, with.stats$n_pruned_edges)
  expect_equal(.local_cpp_edge_keys(no.stats$adj_list, no.stats$weight_list),
               .local_cpp_edge_keys(with.stats$adj_list, with.stats$weight_list))
})

test_that("C++ local pruning validates malformed graph payloads", {
  X <- matrix(0:2, ncol = 1)
  graph <- .local_cpp_graph(3L, rbind(
    c(1, 2, 1),
    c(2, 3, 1)
  ))

  expect_error(
    .prune.graph.local.geodesic(X, graph$adj_list[-1L], graph$weight_list,
                                k = 1, prune.local.k = 1),
    "length nrow"
  )

  bad <- graph
  bad$adj_list[[1L]] <- c(2L, 4L)
  bad$weight_list[[1L]] <- c(1, 1)
  expect_error(
    .prune.graph.local.geodesic(X, bad$adj_list, bad$weight_list,
                                k = 1, prune.local.k = 1),
    "outside 1..nrow"
  )

  bad <- graph
  bad$weight_list[[1L]][[1L]] <- 0
  expect_error(
    .prune.graph.local.geodesic(X, bad$adj_list, bad$weight_list,
                                k = 1, prune.local.k = 1),
    "positive"
  )

  bad <- graph
  bad$adj_list[[2L]] <- bad$adj_list[[2L]][bad$adj_list[[2L]] != 1L]
  bad$weight_list[[2L]] <- bad$weight_list[[2L]][seq_along(bad$adj_list[[2L]])]
  expect_error(
    .prune.graph.local.geodesic(X, bad$adj_list, bad$weight_list,
                                k = 1, prune.local.k = 1),
    "reciprocal"
  )

  bad <- graph
  bad$weight_list[[2L]][bad$adj_list[[2L]] == 1L] <- 2
  expect_error(
    .prune.graph.local.geodesic(X, bad$adj_list, bad$weight_list,
                                k = 1, prune.local.k = 1),
    "mismatched"
  )
})

test_that("graph-family local pruning paths use the shared C++ helper", {
  X <- matrix(0:3, ncol = 1)
  graphs <- list(
    sknn = create.sknn.graph(
      X, k = 3, prune.method = "local.geodesic",
      prune.tau = 1.01, with.pruned.edge.stats = TRUE
    ),
    mknn = create.mknn.graph(
      X, k = 3, prune.method = "local.geodesic",
      prune.tau = 1.01, with.pruned.edge.stats = TRUE
    ),
    iknn = create.single.iknn.graph(
      X, k = 3, prune.method = "local.geodesic",
      prune.tau = 1.01, with.lifecycle.branches = TRUE,
      verbose = FALSE
    ),
    radius = create.rknn.graph(
      X, type = "fixed", radius = 2.1, prune.method = "local.geodesic",
      prune.tau = 1.01, with.pruned.edge.stats = TRUE
    ),
    adaptive = create.rknn.graph(
      X, type = "adaptive.radius", k.scale = 2, radius.rule = "max",
      prune.method = "local.geodesic", prune.tau = 1.01,
      with.pruned.edge.stats = TRUE
    )
  )

  for (graph in graphs) {
    expect_equal(graph$prune_method, "local.geodesic")
    expect_true(graph$n_edges_before_pruning >= graph$n_edges_after_pruning)
    expect_true(graph$n_pruned_edges >= 0L)
    expect_equal(length(graph$repaired_pruned_adj_list), nrow(X))
    expect_equal(length(graph$repaired_pruned_weight_list), nrow(X))
    expect_s3_class(graph$repaired_pruned_pruning$pruned_edge_stats, "data.frame")
    expect_equal(graph$repaired_pruned_pruning$prune_tau, 1.01)
  }
})

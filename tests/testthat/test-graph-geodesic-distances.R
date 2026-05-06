graph_geodesic_edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}


test_that("graph.geodesic.distances dispatches on standard final graph payloads", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  graphs <- list(
    sknn = create.sknn.graph(X, k = 1),
    mknn = create.mknn.graph(X, k = 2, connect.components = TRUE),
    radius = create.radius.graph(X, radius = 1.1, connect.components = TRUE),
    adaptive = create.adaptive.radius.graph(X, k.scale = 1, connect.components = TRUE)
  )

  for (g in graphs) {
    expected <- shortest.path(g$adj_list, g$weight_list, seq_along(g$adj_list))
    observed <- graph.geodesic.distances(g)
    expect_equal(observed, expected, tolerance = 1e-12)
  }
})


test_that("graph.geodesic.distances uses adj_list final payload after local pruning", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  graphs <- list(
    sknn = create.sknn.graph(X, k = 2, prune.method = "local.geodesic",
                             prune.tau = 1.01),
    mknn = create.mknn.graph(X, k = 3, prune.method = "local.geodesic",
                             prune.tau = 1.01)
  )

  for (g in graphs) {
    expect_equal(graph_geodesic_edge_keys(.graph.edge.table(g$adj_list, g$weight_list)),
                 c("1-2", "2-3", "3-4"))
    expected <- shortest.path(g$adj_list, g$weight_list, seq_along(g$adj_list))
    observed <- graph.geodesic.distances(g)
    expect_equal(observed, expected, tolerance = 1e-12)
  }
})


test_that("graph constructors expose raw, pruned, and final lifecycle fields", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(10, 0)
  )

  graphs <- list(
    sknn = create.sknn.graph(X, k = 1, connect.components = TRUE),
    mknn = create.mknn.graph(X, k = 2, connect.components = TRUE),
    radius = create.radius.graph(X, radius = 1.1, connect.components = TRUE),
    adaptive = create.adaptive.radius.graph(X, k.scale = 1, connect.components = TRUE),
    geodesic_iknn = create.geodesic.iknn.graph(
      create.sknn.graph(X, k = 1, connect.components = TRUE),
      k = 1
    ),
    iknn = create.single.iknn.graph(
      X,
      k = 1,
      connect.components = TRUE,
      prune.method = "none",
      threshold.percentile = 0,
      pca.dim = NULL,
      verbose = FALSE
    )
  )

  for (g in graphs) {
    expect_true(all(c(
      "raw_adj_list", "raw_weight_list",
      "pruned_adj_list", "pruned_weight_list",
      "adj_list", "weight_list"
    ) %in% names(g)))
    expect_equal(length(g$raw_adj_list), length(g$raw_weight_list))
    expect_equal(length(g$pruned_adj_list), length(g$pruned_weight_list))
    expect_equal(length(g$adj_list), length(g$weight_list))
    expect_equal(length(g$adj_list), nrow(X))
  }
})


test_that("graph.geodesic.distances uses final adj_list payload for IkNN objects", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  g <- create.single.iknn.graph(
    X,
    k = 2,
    prune.method = "none",
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )

  expected <- shortest.path(g$adj_list, g$weight_list, seq_along(g$adj_list))
  observed <- graph.geodesic.distances(g)

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("graph.geodesic.distances ignores raw and pruned lifecycle fields", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )

  g <- create.single.iknn.graph(
    X,
    k = 3,
    prune.method = "local.geodesic",
    prune.tau = 1.01,
    threshold.percentile = 0,
    compute.full = TRUE,
    pca.dim = NULL,
    verbose = FALSE
  )
  expected <- shortest.path(g$adj_list, g$weight_list, seq_along(g$adj_list))

  # If the wrapper accidentally used raw_* or pruned_* fields, these mutations
  # would make the call fail. The final payload is adj_list/weight_list.
  g$raw_adj_list <- list("not", "a", "numeric", "payload")
  g$raw_weight_list <- list("not", "a", "numeric", "payload")
  g$pruned_adj_list <- list("not", "a", "numeric", "payload")
  g$pruned_weight_list <- list("not", "a", "numeric", "payload")

  expect_equal(graph.geodesic.distances(g), expected, tolerance = 1e-12)
})


test_that("graph.geodesic.distances supports selected vertices", {
  X <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 0),
    c(3, 0)
  )
  g <- create.sknn.graph(X, k = 1)

  observed <- graph.geodesic.distances(g, vertices = c(1, 4))

  expect_equal(dim(observed), c(2L, 2L))
  expect_equal(observed[1, 2], 3, tolerance = 1e-12)
})


test_that("graph.geodesic.distances validates class and payload", {
  X <- rbind(c(0, 0), c(1, 0), c(2, 0))
  g <- create.sknn.graph(X, k = 1)

  expect_error(graph.geodesic.distances(list(adj_list = g$adj_list,
                                             weight_list = g$weight_list)),
               "must inherit")

  bad <- g
  bad$adj_list <- NULL
  expect_error(graph.geodesic.distances(bad), "cannot be NULL")

  iknn <- create.single.iknn.graph(
    X,
    k = 1,
    prune.method = "none",
    threshold.percentile = 0,
    pca.dim = NULL,
    verbose = FALSE
  )
  iknn$adj_list <- NULL
  expect_error(graph.geodesic.distances(iknn), "adj_list")

  bad <- g
  bad$weight_list <- bad$weight_list[-1]
  expect_error(graph.geodesic.distances(bad), "same size")

  bad <- g
  bad$adj_list[[1]] <- c(2L, 99L)
  bad$weight_list[[1]] <- c(1, 1)
  expect_error(graph.geodesic.distances(bad), "valid 1-based")

  bad <- g
  bad$weight_list[[1]][[1]] <- -1
  expect_error(graph.geodesic.distances(bad), "non-negative")

  expect_error(graph.geodesic.distances(g, vertices = 0), "vertices")
})

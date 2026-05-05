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


test_that("create.sknn.graph validates inputs", {
  X <- matrix(1:6, ncol = 2)

  expect_error(create.sknn.graph(X, k = 0), "positive integer")
  expect_error(create.sknn.graph(X, k = 3), "smaller than nrow")
  expect_error(create.sknn.graph(X, k = 1, connect.components = NA),
               "connect.components")
  expect_error(create.sknn.graph(c("a", "b"), k = 1), "matrix or data frame")
})

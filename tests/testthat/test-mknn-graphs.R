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

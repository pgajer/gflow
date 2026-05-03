count.undirected.edges <- function(adj.list) {
    sum(lengths(adj.list)) / 2
}

test_that("MST completion includes boundary-distance ties", {
    X <- matrix(c(
        0, 0,
        1, 0,
        0, 1,
        1, 1
    ), ncol = 2, byrow = TRUE)

    graph <- create.cmst.graph(X, q.thld = 0.5, pca.dim = NULL, verbose = FALSE)

    expect_equal(count.undirected.edges(graph$mst_adj_list), 3)
    expect_equal(graph$cmst_distance_threshold, 1)
    expect_equal(attr(graph, "cmst_distance_threshold"), 1)
    expect_equal(count.undirected.edges(graph$cmst_adj_list), 4)
    expect_true(all(unlist(graph$cmst_weight_list, use.names = FALSE) <= 1 + 1e-12))
})

test_that("MST completion uses R type-7 quantile threshold", {
    X <- cbind(c(0, 1, 3, 6), 0)

    graph <- create.cmst.graph(X, q.thld = 0.5, pca.dim = NULL, verbose = FALSE)

    expect_equal(sort(graph$mst_edge_weights), c(1, 2, 3))
    expect_equal(graph$cmst_distance_threshold, 2)
    expect_equal(count.undirected.edges(graph$cmst_adj_list), 3)
})

test_that("MST completion handles duplicate points at zero threshold", {
    X <- matrix(c(
        0, 0,
        0, 0,
        0, 0,
        1, 0
    ), ncol = 2, byrow = TRUE)

    graph <- create.cmst.graph(X, q.thld = 0.5, pca.dim = NULL, verbose = FALSE)

    expect_equal(graph$cmst_distance_threshold, 0)
    expect_equal(count.undirected.edges(graph$mst_adj_list), 3)
    expect_equal(count.undirected.edges(graph$cmst_adj_list), 4)
})

test_that("MST completion default is a high-quantile local completion", {
    expect_equal(formals(create.cmst.graph)$q.thld, 0.9)
})

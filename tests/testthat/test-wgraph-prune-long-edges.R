test_that("wgraph.prune.long.edges returns valid 1-based adjacency indices", {
    graph <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    edge.lengths <- list(c(1, 2), c(1, 3), c(2, 3))

    result <- dgraphs::wgraph.prune.long.edges(
        graph,
        edge.lengths,
        alt.path.len.ratio.thld = 1.1,
        use.total.length.constraint = TRUE,
        verbose = FALSE
    )

    if (!inherits(result, "dgraph")) expect_named(result, c("adj_list", "edge_lengths_list", "path_lengths",
                           "edge_lengths"))
    adj <- if(inherits(result,"dgraph")) .test.graph.adj(result) else result$adj_list
    expect_equal(adj, graph)
    expect_true(all(vapply(adj, function(neighbors) {
        all(neighbors >= 1L & neighbors <= length(adj))
    }, logical(1))))
})

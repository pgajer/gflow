test_that("wgraph.prune.long.edges returns valid 1-based adjacency indices", {
    graph <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    edge.lengths <- list(c(1, 2), c(1, 3), c(2, 3))

    result <- wgraph.prune.long.edges(
        graph,
        edge.lengths,
        alt.path.len.ratio.thld = 1.1,
        use.total.length.constraint = TRUE,
        verbose = FALSE
    )

    expect_named(result, c("adj_list", "edge_lengths_list", "path_lengths",
                           "edge_lengths"))
    expect_equal(result$adj_list, graph)
    expect_true(all(vapply(result$adj_list, function(neighbors) {
        all(neighbors >= 1L & neighbors <= length(result$adj_list))
    }, logical(1))))
})

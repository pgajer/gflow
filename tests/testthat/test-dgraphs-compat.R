test_that("connected-components compatibility wrapper delegates to dgraphs", {
    graph <- list(c(2L), c(1L), integer())

    expect_warning(
        observed <- graph.connected.components(graph),
        "dgraphs::graph.connected.components"
    )
    expect_identical(
        observed,
        dgraphs::graph.connected.components(graph)
    )
})

test_that("as_igraph compatibility wrapper delegates to dgraphs", {
    basin.graph <- list(
        adjacency.list = list(2L, 1L),
        weight.list = list(1.5, 1.5),
        intersection.matrix = matrix(c(0, 3, 3, 0), 2L),
        basin.metadata = data.frame(
            label = c("a", "b"),
            type = c("minimum", "maximum"),
            size = c(2L, 3L),
            extremum.vertex = c(1L, 2L),
            extremum.value = c(-1, 1)
        )
    )

    expect_warning(
        observed <- as_igraph(basin.graph),
        "dgraphs::as_igraph"
    )
    expected <- dgraphs::as_igraph(basin.graph)

    expect_identical(igraph::as_edgelist(observed), igraph::as_edgelist(expected))
    expect_identical(igraph::vertex_attr(observed), igraph::vertex_attr(expected))
    expect_identical(igraph::edge_attr(observed), igraph::edge_attr(expected))
})

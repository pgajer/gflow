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

test_that("shortest-path compatibility wrappers delegate to dgraphs", {
    graph <- list(2L, c(1L, 3L), 2L)
    edge.lengths <- list(1, c(1, 2), 2)

    expect_warning(
        observed <- shortest.path(graph, edge.lengths, c(1L, 3L)),
        "dgraphs::shortest.path"
    )
    expect_identical(
        observed,
        dgraphs::shortest.path(graph, edge.lengths, c(1L, 3L))
    )

    graph.object <- structure(
        list(adj_list = graph, weight_list = edge.lengths),
        class = "sknn_graph"
    )
    expect_warning(
        observed <- graph.geodesic.distances(graph.object),
        "dgraphs::graph.geodesic.distances"
    )
    expect_identical(
        observed,
        dgraphs::graph.geodesic.distances(graph.object)
    )
})

test_that("endpoint compatibility wrapper delegates to dgraphs", {
    graph <- list(2L, c(1L, 3L), 2L)
    edge.lengths <- list(1, c(1, 1), 1)
    arguments <- list(
        adj.list = graph,
        weight.list = edge.lengths,
        use.approx.eccentricity = FALSE
    )

    expect_warning(
        observed <- do.call(geodesic.core.endpoints, arguments),
        "dgraphs::geodesic.core.endpoints"
    )
    expect_identical(
        observed,
        do.call(dgraphs::geodesic.core.endpoints, arguments)
    )
})

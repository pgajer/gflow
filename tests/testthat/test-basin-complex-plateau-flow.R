plateau.flow.graph <- function(n, edges) {
    adj.list <- vector("list", n)
    edge.length.list <- vector("list", n)
    for (row in seq_len(nrow(edges))) {
        from <- as.integer(edges[row, 1L])
        to <- as.integer(edges[row, 2L])
        edge.length <- edges[row, 3L]
        adj.list[[from]] <- c(adj.list[[from]], to)
        adj.list[[to]] <- c(adj.list[[to]], from)
        edge.length.list[[from]] <- c(
            edge.length.list[[from]],
            edge.length
        )
        edge.length.list[[to]] <- c(
            edge.length.list[[to]],
            edge.length
        )
    }
    for (vertex in seq_len(n)) {
        if (is.null(adj.list[[vertex]])) {
            adj.list[[vertex]] <- integer()
            edge.length.list[[vertex]] <- numeric()
        }
    }
    list(
        adj.list = adj.list,
        edge.length.list = edge.length.list
    )
}

plateau.flow.create <- function(graph,
                                field,
                                direction = "max",
                                edge.length.quantile.thld = 1,
                                long.edge.fallback = "allow_and_flag",
                                ...) {
    create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        field,
        method = "trajectory_flow",
        direction = direction,
        method.params = list(
            modulation = "CLOSEST",
            plateau.policy = "connected_exact",
            edge.length.quantile.thld = edge.length.quantile.thld,
            long.edge.fallback = long.edge.fallback,
            symmetric.seeding = FALSE,
            tie.breaking = FALSE,
            ...
        )
    )
}

test_that("a flow-terminal plateau is defined by adjacent field values", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 1,
                3, 4, 1,
                4, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- plateau.flow.create(
        graph,
        c(0, 2, 2, 2, 0),
        direction = "max"
    )

    expect_identical(object$status, "ok")
    expect_identical(nrow(object$basin.table), 1L)
    expect_identical(
        object$extrema$representative.vertex,
        2L
    )
    expect_identical(
        object$extrema$plateau.vertices[[1L]],
        c(2L, 3L, 4L)
    )
    expect_identical(
        object$extrema$n.plateau.vertices,
        3L
    )
    expect_true(all(object$assignment$root.vertex == 2L))
    expect_true(all(lengths(graph$adj.list)[c(2L, 3L, 4L)] > 1L))
})

test_that("a non-flow-terminal plateau selects the shortest improving edge", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 2,
                3, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- plateau.flow.create(
        graph,
        c(0, 1, 1, 3, 2.5),
        direction = "max"
    )
    forest <- get.basin.trajectory.forest(object, required = TRUE)

    expect_identical(forest$next.vertex$max[2:3], c(3L, 5L))
    expect_identical(forest$root.vertex$max[2:3], c(5L, 5L))
    expect_identical(
        forest$plateau$outgoing.edges$max$plateau_00000002$source.vertex,
        3L
    )
    expect_identical(
        forest$plateau$outgoing.edges$max$plateau_00000002$target.vertex,
        5L
    )
})

test_that("equal boundary lengths resolve by target then source vertex", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 5, 1,
                3, 4, 1,
                3, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- plateau.flow.create(
        graph,
        c(0, 1, 1, 2, 2),
        direction = "max"
    )
    forest <- get.basin.trajectory.forest(object, required = TRUE)
    selected <-
        forest$plateau$outgoing.edges$max$plateau_00000002

    expect_identical(selected$target.vertex, 4L)
    expect_identical(selected$source.vertex, 3L)

    source.tie.graph <- plateau.flow.graph(
        4L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 1,
                3, 4, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    source.tie <- plateau.flow.create(
        source.tie.graph,
        c(0, 1, 1, 2),
        direction = "max"
    )
    source.selected <- get.basin.trajectory.forest(
        source.tie,
        required = TRUE
    )$plateau$outgoing.edges$max$plateau_00000002

    expect_identical(source.selected$target.vertex, 4L)
    expect_identical(source.selected$source.vertex, 2L)
})

test_that("disconnected equal-valued vertices remain distinct extrema", {
    graph <- plateau.flow.graph(
        4L,
        matrix(
            c(
                1, 2, 1,
                1, 3, 1,
                2, 4, 1,
                3, 4, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- plateau.flow.create(
        graph,
        c(0, 1, 1, 0),
        direction = "max"
    )

    expect_setequal(
        object$extrema$representative.vertex,
        c(2L, 3L)
    )
    expect_true(all(object$extrema$n.plateau.vertices == 1L))
    expect_identical(
        object$provenance$construction$n.contracted.plateaus,
        NULL
    )
})

test_that("minimum flow contracts flow-terminal and other plateaus", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 2,
                3, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    object <- plateau.flow.create(
        graph,
        c(3, 2, 2, 0, 1),
        direction = "both"
    )
    forest <- get.basin.trajectory.forest(object, required = TRUE)

    expect_identical(forest$next.vertex$min[2:3], c(3L, 5L))
    expect_identical(forest$root.vertex$min[2:3], c(5L, 5L))
    expect_true(all(!is.na(
        object$assignment$next.vertex[
            object$assignment$direction == "min" &
                object$assignment$vertex %in% c(2L, 3L)
        ]
    )))
    expect_true(all(vapply(
        list(forest$next.vertex$max, forest$next.vertex$min),
        function(next.vertex) {
            roots <- gflow:::.trajectory.follow.roots(
                next.vertex,
                length(next.vertex) + 1L
            )
            length(roots) == length(next.vertex)
        },
        logical(1)
    )))
})

test_that("plateau long-edge fallback is expanded and auditable", {
    graph <- plateau.flow.graph(
        6L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                3, 4, 10,
                4, 5, 1,
                5, 6, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    field <- c(0, 1, 1, 2, 1, 0)
    allowed <- plateau.flow.create(
        graph,
        field,
        direction = "max",
        edge.length.quantile.thld = 0.5,
        long.edge.fallback = "allow_and_flag"
    )
    forbidden <- plateau.flow.create(
        graph,
        field,
        direction = "max",
        edge.length.quantile.thld = 0.5,
        long.edge.fallback = "forbid"
    )
    allowed.forest <- get.basin.trajectory.forest(
        allowed,
        required = TRUE
    )
    forbidden.forest <- get.basin.trajectory.forest(
        forbidden,
        required = TRUE
    )

    expect_identical(allowed.forest$next.vertex$max[2:3], c(3L, 4L))
    expect_true(allowed.forest$long.edge.fallback$next.up.used[[3L]])
    expect_identical(
        allowed.forest$long.edge.fallback$n.next.up,
        1L
    )
    expect_true(allowed.forest$long.edge.fallback$attention.required)

    expect_identical(forbidden.forest$root.vertex$max[2:3], c(2L, 2L))
    expect_true(
        forbidden.forest$long.edge.fallback$next.up.blocked[[2L]]
    )
    expect_identical(
        forbidden.forest$long.edge.fallback$n.next.up.blocked,
        1L
    )
})

test_that("connected-plateau CLOSEST replay is deterministic", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 1,
                3, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    first <- plateau.flow.create(
        graph,
        c(0, 1, 1, 2, 2),
        direction = "both"
    )
    second <- plateau.flow.create(
        graph,
        c(0, 1, 1, 2, 2),
        direction = "both"
    )

    expect_identical(first$membership, second$membership)
    expect_identical(first$assignment, second$assignment)
    expect_identical(
        first$trajectory.forest$next.vertex,
        second$trajectory.forest$next.vertex
    )
    expect_identical(
        first$trajectory.forest$root.vertex,
        second$trajectory.forest$root.vertex
    )
})

test_that("CLOSEST plateau policy is mandatory and seeded replay is stable", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                2, 4, 1,
                3, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    expect_error(
        create.basin.complex(
            graph$adj.list,
            graph$edge.length.list,
            c(0, 1, 1, 2, 2),
            method = "trajectory_flow",
            method.params = list(
                modulation = "CLOSEST",
                plateau.policy = "none"
            )
        ),
        "requires exact connected-plateau contraction",
        class = "gflow_basin_input_error"
    )

    arguments <- list(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        field = c(0, 1, 1, 2, 2),
        method = "trajectory_flow",
        direction = "both",
        method.params = list(
            modulation = "CLOSEST",
            plateau.policy = "connected_exact",
            edge.length.quantile.thld = 1,
            symmetric.seeding = FALSE,
            tie.breaking = TRUE,
            tie.seed = 17L
        )
    )
    first <- do.call(create.basin.complex, arguments)
    second <- do.call(create.basin.complex, arguments)

    expect_identical(first$membership, second$membership)
    expect_identical(first$assignment, second$assignment)
    expect_identical(
        first$extrema$plateau.vertices,
        second$extrema$plateau.vertices
    )
    expect_true(first$field$tie.policy$applied)
})

test_that("support refinement maps trajectories to contracted basin ids", {
    graph <- plateau.flow.graph(
        5L,
        matrix(
            c(
                1, 2, 1,
                2, 3, 1,
                3, 4, 1,
                4, 5, 1
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    raw <- plateau.flow.create(
        graph,
        c(0, 1, 1, 2, 3),
        direction = "both"
    )
    counts <- .basin.trajectory.counts(raw)

    expect_true(all(counts >= 1L))

    refined <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        c(0, 1, 1, 2, 3),
        method = "trajectory_flow",
        direction = "both",
        method.params = list(
            modulation = "CLOSEST",
            plateau.policy = "connected_exact",
            edge.length.quantile.thld = 1,
            symmetric.seeding = FALSE,
            tie.breaking = FALSE
        ),
        simplify.params = list(
            support.filter = list(
                enabled = TRUE,
                min.basin.size = 0L,
                min.trajectory.count = 1L,
                min.basin.mass = 0
            )
        )
    )

    expect_true(all(refined$basin.table$retained))
})

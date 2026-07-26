phase.c.graph <- function(n, edges) {
    adj.list <- vector("list", n)
    edge.length.list <- vector("list", n)
    for (i in seq_len(nrow(edges))) {
        from <- as.integer(edges[i, 1L])
        to <- as.integer(edges[i, 2L])
        edge.length <- edges[i, 3L]
        adj.list[[from]] <- c(adj.list[[from]], to)
        adj.list[[to]] <- c(adj.list[[to]], from)
        edge.length.list[[from]] <- c(
            edge.length.list[[from]],
            edge.length
        )
        edge.length.list[[to]] <- c(edge.length.list[[to]], edge.length)
    }
    for (vertex in seq_len(n)) {
        if (is.null(adj.list[[vertex]])) {
            adj.list[[vertex]] <- integer()
            edge.length.list[[vertex]] <- numeric()
        } else {
            order.index <- order(adj.list[[vertex]])
            adj.list[[vertex]] <- as.integer(
                adj.list[[vertex]][order.index]
            )
            edge.length.list[[vertex]] <-
                edge.length.list[[vertex]][order.index]
        }
    }
    list(
        adj.list = adj.list,
        edge.length.list = edge.length.list
    )
}

phase.c.path.fixture <- function() {
    graph <- phase.c.graph(
        7L,
        cbind(1:6, 2:7, rep(1, 6))
    )
    c(
        graph,
        list(
            field = c(0, 3, 1, 2.5, 1, 2, 0),
            vertex.mass = c(1, 2, 1, 3, 1, 2, 1)
        )
    )
}

phase.c.modulation.fixture <- function() {
    graph <- phase.c.graph(
        12L,
        matrix(
            c(
                1, 2, 0.34, 1, 3, 1.95, 2, 3, 2.64, 3, 4, 2.37,
                2, 5, 1.91, 4, 5, 1.75, 5, 6, 0.31, 6, 7, 0.99,
                2, 8, 0.31, 5, 8, 0.65, 7, 8, 0.91, 1, 9, 1.88,
                8, 9, 2.79, 1, 10, 0.41, 2, 10, 2.78, 3, 10, 1.19,
                4, 10, 1, 6, 10, 2.71, 9, 10, 1.32, 3, 11, 1.53,
                4, 11, 0.88, 5, 11, 1.23, 8, 11, 2.12,
                10, 11, 1.98, 3, 12, 2.08, 10, 12, 0.76,
                11, 12, 1.9
            ),
            ncol = 3L,
            byrow = TRUE
        )
    )
    c(
        graph,
        list(
            field = c(
                2.2, 3.58, 3.61, 0.83, 3.47, 1.56,
                2.6, 1.68, 3.52, 2.91, 3.81, 1.69
            ),
            density = c(
                0.25, 2.04, 2.31, 0.97, 0.28, 4.42,
                0.21, 0.2, 0.83, 0.44, 2.92, 0.17
            )
        )
    )
}

phase.c.membership.signature <- function(membership) {
    split(
        membership$basin.id,
        paste(membership$direction, membership$vertex, sep = ":")
    )
}

test_that("geodesic adapter preserves legacy overlapping support", {
    fixture <- phase.c.path.fixture()
    legacy <- compute.basins.of.attraction(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        edge.length.quantile.thld = 1,
        with.trajectories = TRUE
    )
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "geodesic_reachability",
        vertex.mass = fixture$vertex.mass,
        method.params = list(
            edge.length.quantile.thld = 1,
            store.trajectories = TRUE
        )
    )

    expect_identical(object$status, "ok")
    expect_identical(
        object$provenance$method.backend,
        "compute.basins.of.attraction"
    )
    expect_silent(gflow:::.validate.basin.complex(object))
    expect_equal(nrow(object$assignment), 2L * length(fixture$field))
    expect_true(all(
        object$assignment$assignment.status == "not_applicable"
    ))
    expect_true(all(is.na(object$assignment$assignment.weight)))

    canonical.max <- lapply(
        object$basin.table$raw.support.vertices[
            object$basin.table$type == "max"
        ],
        sort
    )
    canonical.min <- lapply(
        object$basin.table$raw.support.vertices[
            object$basin.table$type == "min"
        ],
        sort
    )
    legacy.max <- lapply(
        legacy$lmax_basins,
        function(basin) sort(as.integer(basin$basin_df[, 1L]))
    )
    legacy.min <- lapply(
        legacy$lmin_basins,
        function(basin) sort(as.integer(basin$basin_df[, 1L]))
    )
    expect_setequal(canonical.max, legacy.max)
    expect_setequal(canonical.min, legacy.min)

    weight.sum <- aggregate(
        membership.weight ~ vertex + direction,
        object$membership,
        sum
    )
    expect_equal(weight.sum$membership.weight, rep(1, nrow(weight.sum)))
    expect_equal(
        object$basin.table$raw.support.mass,
        unname(vapply(
            object$basin.table$raw.support.vertices,
            function(vertices) {
                sum(fixture$vertex.mass[vertices]) /
                    sum(fixture$vertex.mass)
            },
            numeric(1)
        ))
    )
})

test_that("geodesic deterministic assignment is explicit", {
    fixture <- phase.c.path.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "geodesic_reachability",
        method.params = list(
            edge.length.quantile.thld = 1,
            primary.assignment.policy =
                "largest_membership_then_extremum"
        )
    )

    expect_true(all(object$assignment$assignment.status == "assigned"))
    expect_true(all(object$assignment$assignment.weight == 1))
    expect_true(all(
        object$membership$is.primary %in% c(TRUE, FALSE)
    ))
    expect_equal(
        vapply(
            split(object$assignment, object$assignment$direction),
            nrow,
            integer(1)
        ),
        c(max = 7L, min = 7L)
    )
})

test_that("trajectory adapter preserves legacy membership and assignment", {
    fixture <- phase.c.path.fixture()
    legacy <- compute.gfc.trajectory(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 1,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 0L,
        min.n.trajectories = 0L,
        store.trajectories = TRUE,
        tie.breaking = FALSE
    )
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        method.params = list(
            modulation = "CLOSEST",
            edge.length.quantile.thld = 1
        )
    )

    expect_identical(object$status, "ok")
    expect_identical(
        object$provenance$method.backend,
        "compute.gfc.trajectory"
    )
    expect_silent(gflow:::.validate.basin.complex(object))
    max.ids <- object$basin.table$basin.id[
        object$basin.table$type == "max"
    ]
    min.ids <- object$basin.table$basin.id[
        object$basin.table$type == "min"
    ]
    expect_identical(
        match(
            object$assignment$basin.id[
                object$assignment$direction == "max"
            ],
            max.ids
        ),
        legacy$max.assignment
    )
    expect_identical(
        match(
            object$assignment$basin.id[
                object$assignment$direction == "min"
            ],
            min.ids
        ),
        legacy$min.assignment
    )
    expect_identical(
        lapply(
            seq_along(legacy$max.membership.all),
            function(vertex) {
                match(
                    object$membership$basin.id[
                        object$membership$direction == "max" &
                            object$membership$vertex == vertex
                    ],
                    max.ids
                )
            }
        ),
        lapply(legacy$max.membership.all, as.integer)
    )
    expect_identical(
        object$trajectory.forest$next.vertex$max,
        legacy$next.up
    )
})

test_that("trajectory direction subsetting returns complete requested rows", {
    fixture <- phase.c.path.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        direction = "max",
        method.params = list(edge.length.quantile.thld = 1)
    )

    expect_identical(unique(object$extrema$type), "max")
    expect_identical(unique(object$membership$direction), "max")
    expect_identical(unique(object$assignment$direction), "max")
    expect_equal(nrow(object$assignment), length(fixture$field))
})

test_that("all five trajectory modulations retain distinct signatures", {
    fixture <- phase.c.modulation.fixture()
    modulations <- c(
        "CLOSEST", "NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN"
    )
    signatures <- vapply(modulations, function(modulation) {
        object <- create.basin.complex(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field,
            method = "trajectory_flow",
            vertex.density = fixture$density,
            method.params = list(
                modulation = modulation,
                edge.length.quantile.thld = 1
            )
        )
        expect_identical(object$status, "ok")
        paste(
            object$assignment$basin.id,
            object$raw$legacy.object$next.up,
            collapse = "|"
        )
    }, character(1))

    expect_length(unique(signatures), length(modulations))
})

test_that("seeded trajectory tie breaking is local and reproducible", {
    fixture <- phase.c.path.fixture()
    tied.field <- c(0, 1, 1, 2, 1, 1, 0)
    set.seed(318L)
    seed.before <- .Random.seed
    first <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        tied.field,
        method = "trajectory_flow",
        method.params = list(
            edge.length.quantile.thld = 1,
            tie.breaking = TRUE,
            tie.seed = 44L
        )
    )
    expect_identical(.Random.seed, seed.before)
    second <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        tied.field,
        method = "trajectory_flow",
        method.params = list(
            edge.length.quantile.thld = 1,
            tie.breaking = TRUE,
            tie.seed = 44L
        )
    )
    expect_identical(.Random.seed, seed.before)
    expect_identical(first$membership, second$membership)
    expect_identical(first$assignment, second$assignment)
    expect_identical(
        first$field$construction.values,
        second$field$construction.values
    )
    expect_true(first$field$tie.policy$applied)
})

test_that("Phase C keeps absent mass distinct from zero-mass vertices", {
    fixture <- phase.c.path.fixture()
    no.mass <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        method.params = list(edge.length.quantile.thld = 1)
    )
    with.zero <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        vertex.mass = c(1, 0, 1, 0, 1, 0, 1),
        method.params = list(edge.length.quantile.thld = 1)
    )

    expect_true(all(is.na(no.mass$basin.table$raw.support.mass)))
    expect_true(all(is.finite(with.zero$basin.table$raw.support.mass)))
    expect_equal(
        sum(with.zero$basin.table$primary.support.mass[
            with.zero$basin.table$type == "max"
        ]),
        1
    )
})

test_that("successful-object validation enforces Phase C invariants", {
    fixture <- phase.c.path.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        method.params = list(edge.length.quantile.thld = 1)
    )

    duplicate.assignment <- object
    duplicate.assignment$assignment[2L, ] <-
        duplicate.assignment$assignment[1L, ]
    expect_error(
        gflow:::.validate.basin.complex(duplicate.assignment),
        class = "gflow_basin_schema_error"
    )

    bad.membership <- object
    bad.membership$membership$membership.weight[[1L]] <- 0.25
    expect_error(
        gflow:::.validate.basin.complex(bad.membership),
        class = "gflow_basin_schema_error"
    )

    bad.status <- object
    bad.status$assignment$assignment.status[[1L]] <- "assigned"
    bad.status$assignment$assignment.weight[[1L]] <- 0
    expect_error(
        gflow:::.validate.basin.complex(bad.status),
        class = "gflow_basin_schema_error"
    )
})

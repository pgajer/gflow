phase.e.graph <- function(n, edges) {
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
    list(
        adj.list = adj.list,
        edge.length.list = edge.length.list
    )
}

phase.e.ring <- function(n) {
    adj.list <- vector("list", n)
    edge.length.list <- vector("list", n)
    for (vertex in seq_len(n)) {
        neighbors <- unique(as.integer(
            ((vertex - 1L + c(-2L, -1L, 1L, 2L)) %% n) + 1L
        ))
        adj.list[[vertex]] <- neighbors
        edge.length.list[[vertex]] <- rep(1, length(neighbors))
    }
    list(
        adj.list = adj.list,
        edge.length.list = edge.length.list,
        field = sin(2 * pi * seq_len(n) / n) +
            0.25 * cos(6 * pi * seq_len(n) / n)
    )
}

phase.e.grid <- function() {
    edges <- matrix(
        c(
            1, 2, 1, 1, 4, 1, 2, 3, 1, 2, 5, 1,
            3, 6, 1, 4, 5, 1, 4, 7, 1, 5, 6, 1,
            5, 8, 1, 6, 9, 1, 7, 8, 1, 8, 9, 1
        ),
        ncol = 3L,
        byrow = TRUE
    )
    c(
        phase.e.graph(9L, edges),
        list(field = c(0, 1, 0, 1, 3, 1, 0, 1, 2))
    )
}

test_that("RTCB adapter preserves legacy supports and resolved parameters", {
    fixture <- phase.e.ring(40L)
    method.params <- list(
        edge.length.quantile.thld = 1,
        store.trajectories = TRUE,
        n.min = 12L,
        m.min = 0.25,
        q.min = 0.2,
        run.max = 1L,
        tau0 = 0.05,
        kappa = 1.2,
        k.max = 8L,
        h.max = 100L,
        d.path.max = 4,
        eta.step = 0.01,
        sink.prune = TRUE,
        max.paths.per.sink = 8L
    )
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "rtcb",
        method.params = method.params
    )
    legacy <- gflow:::compute.basins.of.attraction.rtcb(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        edge.length.quantile.thld = 1,
        with.trajectories = TRUE,
        n.min = 12L,
        m.min = 0.25,
        q.min = 0.2,
        run.max = 1L,
        tau0 = 0.05,
        kappa = 1.2,
        k.max = 8L,
        h.max = 100L,
        d.path.max = 4,
        eta.step = 0.01,
        sink.prune = TRUE,
        max.paths.per.sink = 8L
    )

    expect_identical(object$status, "ok")
    expect_silent(gflow:::.validate.basin.complex(object))
    expect_identical(
        object$provenance$method.backend,
        "compute.basins.of.attraction.rtcb"
    )
    expect_identical(
        object$provenance$construction$backend.parameters,
        legacy$rtcb.params
    )
    expect_true(all(
        object$assignment$assignment.status == "not_applicable"
    ))
    expect_true(all(is.na(object$assignment$assignment.weight)))
    expect_setequal(
        object$basin.table$raw.support.vertices[
            object$basin.table$type == "max"
        ],
        lapply(
            legacy$lmax_basins,
            function(basin) sort(as.integer(basin$basin_df[, 1L]))
        )
    )
    expect_setequal(
        object$basin.table$raw.support.vertices[
            object$basin.table$type == "min"
        ],
        lapply(
            legacy$lmin_basins,
            function(basin) sort(as.integer(basin$basin_df[, 1L]))
        )
    )
})

test_that("RTCB records defaults and nondefaults change support", {
    fixture <- phase.e.ring(40L)
    default <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "rtcb",
        method.params = list(edge.length.quantile.thld = 1)
    )
    nondefault <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "rtcb",
        method.params = list(
            edge.length.quantile.thld = 1,
            n.min = 8L,
            m.min = 0.2,
            run.max = 1L
        )
    )

    expect_identical(default$parameters$method.params$n.min, 20L)
    expect_identical(
        default$provenance$construction$assignment.policy,
        "none"
    )
    expect_false(identical(
        default$basin.table$raw.support.size,
        nondefault$basin.table$raw.support.size
    ))
})

test_that("cell-complex adapter preserves every required legacy structure", {
    fixture <- phase.e.grid()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "overlap_cell_complex"
    )
    legacy <- object$raw$legacy.object

    expect_identical(object$status, "ok")
    expect_silent(gflow:::.validate.basin.complex(object))
    expect_identical(
        class(object),
        c(
            "basin_complex_overlap_cell_complex",
            "basin_complex",
            "gflow_basin_complex",
            "list"
        )
    )
    expect_false(inherits(object, "basin_cx"))
    expect_identical(object$cell.complex$cells, legacy$cells)
    expect_identical(
        object$cell.complex$cluster.assignments,
        legacy$cluster_assignments
    )
    expect_identical(
        object$cell.complex$cluster.mappings,
        legacy$cluster_mappings
    )
    expect_identical(
        object$cell.complex$complex.graph,
        legacy$graph
    )
    expect_identical(
        object$cell.complex$simplified.field,
        as.numeric(legacy$harmonic_predictions)
    )
    expect_identical(
        object$raw$initial.basin.complex,
        legacy$initial_basin_cx
    )
    expect_identical(object$raw$merged.basins, legacy$basins)
    expect_true(all(
        object$assignment$assignment.status == "not_applicable"
    ))
})

test_that("cell-complex nondefault profile is scientifically distinct", {
    fixture <- phase.e.grid()
    default <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "overlap_cell_complex"
    )
    nondefault <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "overlap_cell_complex",
        method.params = list(
            basin.merge.overlap.thld = 0.7,
            min.asc.desc.cell.size.thld = 2L,
            min.asc.asc.cell.size.thld = 2L,
            min.desc.desc.cell.size.thld = 2L
        )
    )

    expect_false(identical(
        default$basin.table$raw.support.vertices,
        nondefault$basin.table$raw.support.vertices
    ))
    expect_false(identical(
        default$cell.complex$cell.summary,
        nondefault$cell.complex$cell.summary
    ))
})

test_that("Phase E validation rejects parameter and cell payload drift", {
    ring <- phase.e.ring(30L)
    rtcb <- create.basin.complex(
        ring$adj.list,
        ring$edge.length.list,
        ring$field,
        method = "rtcb",
        method.params = list(edge.length.quantile.thld = 1)
    )
    bad.rtcb <- rtcb
    bad.rtcb$provenance$construction$backend.parameters$n_min <- 1L
    expect_error(
        gflow:::.validate.basin.complex(bad.rtcb),
        class = "gflow_basin_schema_error"
    )

    grid <- phase.e.grid()
    cells <- create.basin.complex(
        grid$adj.list,
        grid$edge.length.list,
        grid$field,
        method = "overlap_cell_complex"
    )
    bad.cells <- cells
    bad.cells$cell.complex$simplified.field <- numeric()
    expect_error(
        gflow:::.validate.basin.complex(bad.cells),
        class = "gflow_basin_schema_error"
    )
})

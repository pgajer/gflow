phase.g.grid.fixture <- function() {
    adjacency <- list(
        c(2L, 4L),
        c(1L, 3L, 5L),
        c(2L, 6L),
        c(1L, 5L, 7L),
        c(2L, 4L, 6L, 8L),
        c(3L, 5L, 9L),
        c(4L, 8L),
        c(5L, 7L, 9L),
        c(6L, 8L)
    )
    list(
        adj.list = adjacency,
        edge.length.list = lapply(
            adjacency,
            function(vertices) rep(1, length(vertices))
        ),
        field = c(5.0, 3.0, 0.5, 2.5, 1.2, 0.7, 0.6, 2.2, 4.3)
    )
}

test_that("canonical adapters do not dispatch through legacy exports", {
    fixture <- phase.g.grid.fixture()
    testthat::local_mocked_bindings(
        compute.basins.of.attraction = function(...) {
            stop("legacy geodesic export called")
        },
        compute.gfc.trajectory = function(...) {
            stop("legacy trajectory export called")
        },
        create.basin.cx = function(...) {
            stop("legacy cell export called")
        },
        .package = "gflow"
    )

    geodesic <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "geodesic_reachability",
        direction = "max",
        method.params = list(edge.length.quantile.thld = 1)
    )
    trajectory <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        direction = "max",
        method.params = list(
            modulation = "CLOSEST",
            edge.length.quantile.thld = 1,
            tie.breaking = FALSE
        )
    )
    cells <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "overlap_cell_complex"
    )

    expect_identical(geodesic$status, "ok")
    expect_identical(trajectory$status, "ok")
    expect_identical(cells$status, "ok")
})

test_that("two-dimensional canonical migration preserves assignment contract", {
    fixture <- phase.g.grid.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        direction = "max",
        method.params = list(
            modulation = "CLOSEST",
            edge.length.quantile.thld = 1,
            tie.breaking = FALSE
        )
    )
    maxima <- object$basin.table$extremum.vertex[
        object$basin.table$type == "max"
    ]
    assignment <- object$assignment[
        object$assignment$direction == "max",
        ,
        drop = FALSE
    ]

    expect_setequal(maxima, c(1L, 9L))
    expect_equal(nrow(assignment), 9L)
    expect_true(all(assignment$assignment.status == "assigned"))
    expect_true(all(assignment$basin.id %in% object$basin.table$basin.id))
    expect_silent(gflow:::.validate.basin.complex(object))
})

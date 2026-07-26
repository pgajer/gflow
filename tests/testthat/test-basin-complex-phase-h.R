phase.h.fixture <- function() {
    list(
        adj.list = list(2L, c(1L, 3L), 2L),
        edge.length.list = list(1, c(1, 1), 1),
        field = c(3, 1, 2)
    )
}

phase.h.capture <- function(expr) {
    tryCatch(expr, error = identity)
}

test_that("retired constructors fail with classed canonical guidance", {
    fixture <- phase.h.fixture()
    cases <- list(
        geodesic = phase.h.capture(compute.basins.of.attraction(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field
        )),
        trajectory = phase.h.capture(compute.gfc.trajectory(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field
        )),
        flow.alias = phase.h.capture(compute.gfc.flow(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field
        )),
        cells = phase.h.capture(create.basin.cx(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field
        )),
        mixed = phase.h.capture(compute.gfc(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field
        ))
    )

    for (condition in cases) {
        expect_s3_class(condition, "gflow_basin_lifecycle_error")
        expect_match(conditionMessage(condition), "retired")
        expect_match(conditionMessage(condition), "create.basin.complex")
        expect_true(nzchar(condition$details$argument.translation))
    }
    expect_identical(
        cases$geodesic$details$canonical.method,
        "geodesic_reachability"
    )
    expect_identical(
        cases$trajectory$details$canonical.method,
        "trajectory_flow"
    )
    expect_identical(
        cases$cells$details$canonical.method,
        "overlap_cell_complex"
    )
})

test_that("compute.gfc lifecycle guidance follows legacy modulation", {
    fixture <- phase.h.fixture()
    geodesic <- phase.h.capture(compute.gfc(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        modulation = "GEODESIC"
    ))
    rtcb <- phase.h.capture(compute.gfc(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        modulation = "RTCB"
    ))

    expect_identical(
        geodesic$details$canonical.method,
        "geodesic_reachability"
    )
    expect_identical(rtcb$details$canonical.method, "rtcb")
})

test_that("unresolved compute.gfc.basins backend remains functional", {
    fixture <- phase.h.fixture()
    object <- compute.gfc.basins(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        modulation = "none",
        verbose = FALSE
    )

    expect_s3_class(object, "gfc_basins")
})

test_that("archived trajectory objects reconstruct canonically", {
    fixture <- phase.h.fixture()
    archived <- gflow:::.compute.gfc.trajectory.backend(
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
        tie.breaking = FALSE
    )
    converted <- as.basin.complex(
        archived,
        adj.list = fixture$adj.list,
        edge.length.list = fixture$edge.length.list,
        method.params = list(edge.length.quantile.thld = 1)
    )

    expect_s3_class(converted, "gflow_basin_complex")
    expect_identical(converted$method, "trajectory_flow")
    expect_identical(converted$raw$archived.object, archived)
    expect_identical(
        converted$provenance$conversion$source.class,
        "gfc.flow"
    )
})

test_that("archived geodesic and cell objects reconstruct canonically", {
    fixture <- phase.h.fixture()
    geodesic <- gflow:::.compute.basins.of.attraction.backend(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        edge.length.quantile.thld = 1
    )
    cells <- gflow:::.create.basin.cx.backend(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field
    )
    converted.geodesic <- as.basin.complex(
        geodesic,
        adj.list = fixture$adj.list,
        edge.length.list = fixture$edge.length.list,
        method.params = list(edge.length.quantile.thld = 1)
    )
    converted.cells <- as.basin.complex(
        cells,
        adj.list = fixture$adj.list,
        edge.length.list = fixture$edge.length.list
    )

    expect_identical(converted.geodesic$method, "geodesic_reachability")
    expect_identical(converted.cells$method, "overlap_cell_complex")
    expect_identical(converted.geodesic$raw$archived.object, geodesic)
    expect_identical(converted.cells$raw$archived.object, cells)
})

test_that("archived conversion reports missing graph or field", {
    fixture <- phase.h.fixture()
    geodesic <- gflow:::.compute.basins.of.attraction.backend(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field
    )
    legacy.mixed <- gflow:::.compute.gfc.legacy.backend(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        edge.length.quantile.thld = 1,
        min.basin.size = 0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        expand.basins = FALSE,
        modulation = "GEODESIC"
    )

    expect_error(
        as.basin.complex(geodesic),
        class = "gflow_basin_conversion_error"
    )
    expect_error(
        as.basin.complex(
            legacy.mixed,
            adj.list = fixture$adj.list,
            edge.length.list = fixture$edge.length.list
        ),
        class = "gflow_basin_conversion_error"
    )
    expect_s3_class(
        as.basin.complex(
            legacy.mixed,
            adj.list = fixture$adj.list,
            edge.length.list = fixture$edge.length.list,
            field = fixture$field,
            method.params = list(edge.length.quantile.thld = 1)
        ),
        "gflow_basin_complex"
    )
})

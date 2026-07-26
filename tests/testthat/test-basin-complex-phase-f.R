phase.f.fixture <- function() {
    list(
        adj.list = list(
            2L,
            c(1L, 3L),
            c(2L, 4L),
            c(3L, 5L),
            c(4L, 6L),
            c(5L, 7L),
            6L
        ),
        edge.length.list = list(
            1,
            c(1, 1.2),
            c(1.2, 0.8),
            c(0.8, 1.4),
            c(1.4, 0.7),
            c(0.7, 1.1),
            1.1
        ),
        field = c(0.2, 3.1, 1.2, 2.7, 1.1, 2.3, 0.1),
        mass = c(1, 2, 1, 3, 1, 2, 1)
    )
}

phase.f.legacy.geodesic <- function(fixture, ...) {
    args <- list(
        adj.list = fixture$adj.list,
        edge.length.list = fixture$edge.length.list,
        fitted.values = fixture$field,
        modulation = "GEODESIC",
        edge.length.quantile.thld = 1,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 0L,
        expand.basins = FALSE,
        with.trajectories = FALSE,
        verbose = FALSE
    )
    do.call(
        gflow:::.compute.gfc.legacy.backend,
        utils::modifyList(args, list(...))
    )
}

phase.f.canonical.geodesic <- function(fixture, simplify.params) {
    create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "geodesic_reachability",
        vertex.mass = fixture$mass,
        method.params = list(edge.length.quantile.thld = 1),
        simplify.params = simplify.params
    )
}

phase.f.retained.signature <- function(object) {
    table <- object$basin.table[object$basin.table$retained, ]
    result <- setNames(
        lapply(
            unclass(table$retained.support.vertices),
            function(vertices) sort(as.integer(vertices))
        ),
        paste(table$type, table$extremum.vertex, sep = ":")
    )
    result[order(names(result))]
}

phase.f.legacy.signature <- function(object, expanded = FALSE) {
    max.support <- if (expanded && !is.null(object$expanded.max.vertices.list)) {
        object$expanded.max.vertices.list
    } else {
        object$max.vertices.list
    }
    min.support <- if (expanded && !is.null(object$expanded.min.vertices.list)) {
        object$expanded.min.vertices.list
    } else {
        object$min.vertices.list
    }
    max.vertices <- object$summary$vertex[object$summary$type == "max"]
    min.vertices <- object$summary$vertex[object$summary$type == "min"]
    result <- c(
        setNames(
            unname(max.support),
            paste0("max:", max.vertices)
        ),
        setNames(
            unname(min.support),
            paste0("min:", min.vertices)
        )
    )
    result <- lapply(result, function(vertices) sort(as.integer(vertices)))
    result[order(names(result))]
}

test_that("refinement provenance records every stage explicitly", {
    fixture <- phase.f.fixture()
    object <- phase.f.canonical.geodesic(
        fixture,
        list(
            relative.value = list(
                enabled = TRUE,
                min.relative.value.max = 1.7,
                max.relative.value.min = 0.5
            ),
            support.filter = list(
                enabled = TRUE,
                min.basin.size = 2L
            )
        )
    )
    stages <- object$provenance$refinement.stages

    expect_identical(
        stages$stage,
        c(
            "relative.value",
            "maxima.clustering",
            "minima.clustering",
            "geometric.filter",
            "support.filter",
            "expansion"
        )
    )
    expect_identical(
        stages$status,
        c(
            "completed",
            "disabled",
            "disabled",
            "disabled",
            "completed",
            "disabled"
        )
    )
    expect_true(all(vapply(
        unclass(stages$summary.snapshot),
        is.data.frame,
        logical(1)
    )))
    expect_silent(gflow:::.validate.basin.complex(object))
})

test_that("relative-value stage matches legacy geodesic retention", {
    fixture <- phase.f.fixture()
    legacy <- phase.f.legacy.geodesic(
        fixture,
        apply.relvalue.filter = TRUE,
        min.rel.value.max = 1.7,
        max.rel.value.min = 0.5
    )
    canonical <- phase.f.canonical.geodesic(
        fixture,
        list(relative.value = list(
            enabled = TRUE,
            min.relative.value.max = 1.7,
            max.relative.value.min = 0.5
        ))
    )

    expect_identical(
        phase.f.retained.signature(canonical),
        phase.f.legacy.signature(legacy)
    )
})

test_that("maxima and minima clustering match legacy stage output", {
    fixture <- phase.f.fixture()
    profiles <- list(
        maxima = list(
            legacy = list(
                apply.maxima.clustering = TRUE,
                max.overlap.threshold = 0.8
            ),
            canonical = list(maxima.clustering = list(
                enabled = TRUE,
                overlap.threshold = 0.8
            ))
        ),
        minima = list(
            legacy = list(
                apply.minima.clustering = TRUE,
                min.overlap.threshold = 0.8
            ),
            canonical = list(minima.clustering = list(
                enabled = TRUE,
                overlap.threshold = 0.8
            ))
        )
    )
    for (profile in profiles) {
        legacy <- do.call(
            phase.f.legacy.geodesic,
            c(list(fixture = fixture), profile$legacy)
        )
        canonical <- phase.f.canonical.geodesic(
            fixture,
            profile$canonical
        )
        expect_identical(
            phase.f.retained.signature(canonical),
            phase.f.legacy.signature(legacy)
        )
    }
})

test_that("geometric filtering matches legacy on reviewed metrics", {
    fixture <- phase.f.fixture()
    legacy <- phase.f.legacy.geodesic(
        fixture,
        apply.geometric.filter = TRUE,
        p.mean.nbrs.dist.threshold = 0.9,
        p.mean.hopk.dist.threshold = 0.9,
        p.deg.threshold = 1,
        hop.k = 2L
    )
    canonical <- phase.f.canonical.geodesic(
        fixture,
        list(geometric.filter = list(
            enabled = TRUE,
            mean.neighbor.distance.percentile.max = 0.9,
            mean.hop.distance.percentile.max = 0.9,
            degree.percentile.max = 1,
            hop.k = 2L
        ))
    )

    expect_identical(
        phase.f.retained.signature(canonical),
        phase.f.legacy.signature(legacy)
    )
})

test_that("expansion matches the legacy nearest-basin stage", {
    fixture <- phase.f.fixture()
    legacy <- phase.f.legacy.geodesic(
        fixture,
        apply.relvalue.filter = TRUE,
        min.rel.value.max = 1.7,
        max.rel.value.min = 0.5,
        expand.basins = TRUE
    )
    canonical <- phase.f.canonical.geodesic(
        fixture,
        list(
            relative.value = list(
                enabled = TRUE,
                min.relative.value.max = 1.7,
                max.relative.value.min = 0.5
            ),
            expansion = list(enabled = TRUE)
        )
    )

    expect_identical(
        phase.f.retained.signature(canonical),
        phase.f.legacy.signature(legacy, expanded = TRUE)
    )
})

test_that("raw membership remains unchanged through refinement", {
    fixture <- phase.f.fixture()
    raw <- phase.f.canonical.geodesic(fixture, list())
    refined <- phase.f.canonical.geodesic(
        fixture,
        list(
            maxima.clustering = list(
                enabled = TRUE,
                overlap.threshold = 0.8
            ),
            minima.clustering = list(
                enabled = TRUE,
                overlap.threshold = 0.8
            )
        )
    )

    expect_identical(refined$membership, raw$membership)
    expect_true(any(!refined$basin.table$retained))
    expect_true(all(
        refined$basin.table$retained.support.size[
            !refined$basin.table$retained
        ] == 0L
    ))
})

test_that("support filtering updates assignment and conserves retained mass", {
    fixture <- phase.f.fixture()
    trajectory <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "trajectory_flow",
        vertex.mass = fixture$mass,
        method.params = list(edge.length.quantile.thld = 1),
        simplify.params = list(support.filter = list(
            enabled = TRUE,
            min.basin.size = 3L,
            min.trajectory.count = 0L,
            min.basin.mass = 0.2
        ))
    )

    expect_true(all(
        trajectory$basin.table$retained.support.size[
            trajectory$basin.table$retained
        ] >= 3L
    ))
    expect_true(all(
        trajectory$basin.table$retained.support.mass[
            trajectory$basin.table$retained
        ] >= 0.2
    ))
    expect_true(all(
        trajectory$assignment$assignment.weight %in% c(0, 1)
    ))
})

test_that("positive mass filtering requires explicit vertex mass", {
    fixture <- phase.f.fixture()
    expect_error(
        create.basin.complex(
            fixture$adj.list,
            fixture$edge.length.list,
            fixture$field,
            method = "geodesic_reachability",
            simplify.params = list(support.filter = list(
                enabled = TRUE,
                min.basin.mass = 0.1
            ))
        ),
        class = "gflow_basin_input_error"
    )
})

test_that("non-applicable refinement stages are explicit", {
    fixture <- phase.f.fixture()
    object <- create.basin.complex(
        fixture$adj.list,
        fixture$edge.length.list,
        fixture$field,
        method = "superlevel_merge_tree"
    )
    stages <- object$provenance$refinement.stages

    expect_identical(
        stages$status,
        c(
            "not_applicable",
            "not_applicable",
            "not_applicable",
            "not_applicable",
            "disabled",
            "not_applicable"
        )
    )
})

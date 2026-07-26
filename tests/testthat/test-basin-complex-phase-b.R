phase.b.graph <- function() {
    list(
        adj.list = list(c(2L, 3L), c(1L, 3L), c(1L, 2L), integer()),
        edge.length.list = list(c(1, 2), c(1, 1.5), c(2, 1.5), numeric()),
        field = c(0, 1, 0.5, 2)
    )
}

phase.b.create <- function(method = "trajectory_flow", ...) {
    graph <- phase.b.graph()
    create.basin.complex(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        field = graph$field,
        method = method,
        ...
    )
}

test_that("Phase B constructor returns canonical failed objects for every method", {
    methods <- c(
        "trajectory_flow",
        "superlevel_merge_tree",
        "geodesic_reachability",
        "rtcb",
        "overlap_cell_complex"
    )

    for (method in methods) {
        object <- phase.b.create(method)

        expect_identical(
            class(object),
            c(
                paste0("basin_complex_", method),
                "basin_complex",
                "gflow_basin_complex",
                "list"
            )
        )
        expect_identical(object$method, method)
        expect_identical(object$direction, "both")
        expect_identical(object$status, "failed")
        expect_identical(nrow(object$diagnostics), 1L)
        expect_identical(
            object$diagnostics$condition.class,
            "gflow_basin_backend_not_implemented"
        )
        expect_true(gflow:::.validate.basin.complex(object))
        expect_false(any(class(object) %in% c(
            "basin_cx",
            "gfc.flow",
            "gfc_basins",
            "basins_of_attraction"
        )))
    }
})

test_that("graph validation accepts disconnected and zero-length geometry", {
    graph <- phase.b.graph()
    graph$edge.length.list[[1L]][1L] <- 0
    graph$edge.length.list[[2L]][1L] <- 0

    object <- create.basin.complex(
        graph$adj.list,
        graph$edge.length.list,
        graph$field
    )

    expect_identical(object$n.vertices, 4L)
    expect_identical(object$graph.input$validation$n.edges, 3L)
    expect_identical(object$graph.input$validation$n.components, 2L)
    expect_identical(object$graph.input$validation$isolated.vertices, 4L)
    expect_true(object$graph.input$validation$has.zero.length.edges)
    expect_true(object$graph.input$validation$edge.length.symmetric)
})

test_that("graph validation rejects malformed undirected graphs", {
    graph <- phase.b.graph()

    expect_error(
        create.basin.complex(list(), list(), numeric()),
        class = "gflow_basin_input_error"
    )
    expect_error(
        create.basin.complex(graph$adj.list, graph$edge.length.list[-1L], graph$field),
        class = "gflow_basin_input_error"
    )

    bad <- graph
    bad$adj.list[[1L]] <- c(2L, 2L)
    bad$edge.length.list[[1L]] <- c(1, 1)
    expect_error(
        create.basin.complex(bad$adj.list, bad$edge.length.list, bad$field),
        "duplicate",
        class = "gflow_basin_input_error"
    )

    bad <- graph
    bad$adj.list[[1L]][1L] <- 1L
    expect_error(
        create.basin.complex(bad$adj.list, bad$edge.length.list, bad$field),
        "self-loop",
        class = "gflow_basin_input_error"
    )

    bad <- graph
    bad$adj.list[[2L]] <- 3L
    bad$edge.length.list[[2L]] <- 1.5
    expect_error(
        create.basin.complex(bad$adj.list, bad$edge.length.list, bad$field),
        "reverse edge",
        class = "gflow_basin_input_error"
    )

    bad <- graph
    bad$edge.length.list[[2L]][1L] <- 1.01
    expect_error(
        create.basin.complex(bad$adj.list, bad$edge.length.list, bad$field),
        "symmetry tolerance",
        class = "gflow_basin_input_error"
    )
    expect_s3_class(
        create.basin.complex(
            bad$adj.list,
            bad$edge.length.list,
            bad$field,
            graph.params = list(edge.length.symmetry.tolerance = 0.02)
        ),
        "basin_complex"
    )
})

test_that("field, mass, and density have separate strict contracts", {
    graph <- phase.b.graph()

    expect_error(
        create.basin.complex(
            graph$adj.list,
            graph$edge.length.list,
            c(graph$field[-1L], NA_real_)
        ),
        class = "gflow_basin_input_error"
    )
    expect_error(
        create.basin.complex(
            graph$adj.list,
            graph$edge.length.list,
            graph$field,
            vertex.mass = rep(0, 4L)
        ),
        "positive total",
        class = "gflow_basin_input_error"
    )
    expect_error(
        create.basin.complex(
            graph$adj.list,
            graph$edge.length.list,
            graph$field,
            vertex.density = c(1, 1, -1, 1)
        ),
        "nonnegative",
        class = "gflow_basin_input_error"
    )

    without.mass <- phase.b.create()
    expect_null(without.mass$field$vertex.mass.input)
    expect_null(without.mass$field$vertex.mass.normalized)
    expect_true(is.na(without.mass$field$vertex.mass.input.total))

    with.mass <- phase.b.create(
        vertex.mass = c(0, 2, 1, 1),
        vertex.density = c(10, 8, 4, 2)
    )
    expect_equal(with.mass$field$vertex.mass.normalized, c(0, 0.5, 0.25, 0.25))
    expect_identical(with.mass$field$vertex.mass.input.total, 4)
    expect_equal(with.mass$field$vertex.density, c(10, 8, 4, 2))
})

test_that("method parameters are resolved and validated by family", {
    expect_error(
        phase.b.create(method.params = list(unknown = 1)),
        "Unknown method.params",
        class = "gflow_basin_input_error"
    )
    expect_error(
        phase.b.create(
            method.params = list(modulation = "DENSITY")
        ),
        "vertex.density",
        class = "gflow_basin_input_error"
    )
    density.object <- phase.b.create(
        method.params = list(modulation = "DENSITY"),
        vertex.density = c(1, 2, 3, 4)
    )
    expect_identical(
        density.object$parameters$method.params$modulation,
        "DENSITY"
    )

    expect_error(
        phase.b.create(
            method.params = list(tie.breaking = TRUE)
        ),
        "tie.seed",
        class = "gflow_basin_input_error"
    )
    tied <- phase.b.create(
        method.params = list(tie.breaking = TRUE, tie.seed = 17L)
    )
    expect_identical(tied$field$tie.policy$seed, 17L)
    expect_false(tied$field$tie.policy$applied)

    rtcb <- phase.b.create(
        "rtcb",
        method.params = list(n.min = NULL)
    )
    expect_identical(rtcb$parameters$method.params$n.min, 20L)
    expect_error(
        phase.b.create(
            "superlevel_merge_tree",
            method.params = list(plateau.tolerance = 0.01)
        ),
        "exact plateau",
        class = "gflow_basin_input_error"
    )
    expect_error(
        phase.b.create("overlap_cell_complex", direction = "max"),
        "requires direction",
        class = "gflow_basin_input_error"
    )
})

test_that("refinement configuration is strict even while disabled", {
    expect_error(
        phase.b.create(
            simplify.params = list(
                relative.value = list(unknown = 1)
            )
        ),
        "Unknown",
        class = "gflow_basin_input_error"
    )
    expect_error(
        phase.b.create(
            "superlevel_merge_tree",
            simplify.params = list(
                relative.value = list(enabled = TRUE)
            )
        ),
        "not applicable",
        class = "gflow_basin_input_error"
    )
    expect_error(
        phase.b.create(
            "superlevel_merge_tree",
            simplify.params = list(
                relative.value = list(
                    enabled = FALSE,
                    min.relative.value.max = 2
                )
            )
        ),
        "disabled nondefault",
        class = "gflow_basin_input_error"
    )
    expect_error(
        phase.b.create(
            simplify.params = list(
                expansion = list(enabled = TRUE)
            )
        ),
        "not applicable",
        class = "gflow_basin_input_error"
    )
    expect_s3_class(
        phase.b.create(
            "geodesic_reachability",
            simplify.params = list(
                expansion = list(enabled = TRUE)
            )
        ),
        "basin_complex"
    )
})

test_that("canonical tables have stable columns and accessors", {
    object <- phase.b.create()

    expect_identical(
        names(get.basin.table(object)),
        names(gflow:::.empty.basin.table())
    )
    expect_identical(
        names(get.basin.membership(object)),
        names(gflow:::.empty.membership.table())
    )
    expect_identical(
        names(get.basin.assignment(object)),
        names(gflow:::.empty.assignment.table())
    )
    expect_identical(as.data.frame(object), get.basin.table(object))
    expect_identical(as.basin.complex(object), object)
    expect_identical(basins(object), list())

    expect_true(is.list(get.basin.trajectory.forest(object)))
    expect_null(get.basin.merge.tree(object))
    expect_null(get.basin.cells(object))
    expect_error(
        get.basin.merge.tree(object, required = TRUE),
        class = "gflow_basin_input_error"
    )

    tree <- phase.b.create("superlevel_merge_tree")
    expect_s3_class(get.basin.merge.tree(tree), "data.frame")
    expect_identical(nrow(get.basin.merge.tree(tree)), 0L)

    cells <- phase.b.create("overlap_cell_complex")
    expect_true(is.list(get.basin.cells(cells)))
    expect_named(
        get.basin.cells(cells),
        c(
            "cells",
            "basin.intersections",
            "complex.graph",
            "cell.summary",
            "cluster.assignments",
            "cluster.mappings",
            "simplified.field"
        )
    )
})

test_that("print, summary, plot, and conversion methods are stable", {
    object <- phase.b.create()

    expect_output(print(object), "<basin_complex>")
    summary.object <- summary(object)
    expect_s3_class(summary.object, "summary.basin_complex")
    expect_identical(summary.object$n.vertices, 4L)
    expect_identical(summary.object$n.components, 2L)
    expect_false(summary.object$has.vertex.mass)
    expect_output(print(summary.object), "Canonical Basin Complex Summary")

    plot.file <- tempfile(fileext = ".pdf")
    grDevices::pdf(plot.file)
    expect_invisible(plot(object))
    grDevices::dev.off()
    expect_gt(file.info(plot.file)$size, 0)
    unlink(plot.file)

    expect_error(
        as.basin.complex(list()),
        class = "gflow_basin_conversion_error"
    )
})

test_that("schema validation rejects corrupted canonical objects", {
    object <- phase.b.create()

    bad.class <- object
    class(bad.class) <- rev(class(bad.class))
    expect_error(
        gflow:::.validate.basin.complex(bad.class),
        "classes",
        class = "gflow_basin_schema_error"
    )

    bad.table <- object
    bad.table$basin.table$extra <- character()
    expect_error(
        gflow:::.validate.basin.complex(bad.table),
        "schema",
        class = "gflow_basin_schema_error"
    )

    bad.failure <- object
    bad.failure$diagnostics <- gflow:::.empty.diagnostics.table()
    expect_error(
        gflow:::.validate.basin.complex(bad.failure),
        "diagnostic",
        class = "gflow_basin_schema_error"
    )
})

test_that("Phase B does not redirect or dual-class legacy constructors", {
    graph <- phase.b.graph()
    legacy <- compute.basins.of.attraction(
        graph$adj.list,
        graph$edge.length.list,
        graph$field,
        edge.length.quantile.thld = 1
    )

    expect_s3_class(legacy, "basins_of_attraction")
    expect_false(inherits(legacy, "basin_complex"))
})

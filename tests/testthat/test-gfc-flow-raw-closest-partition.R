raw_closest_args <- function(adj.list, weight.list, y, ...) {
    overrides <- list(...)
    args <- list(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 1,
        long.edge.fallback = "flag",
        plateau.tolerance = 0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        symmetric.seeding = FALSE,
        tie.breaking = FALSE
    )
    args[names(overrides)] <- overrides
    do.call(compute.gfc.trajectory, args)
}

raw_path_graph <- function(n, weights = rep(1, n - 1L)) {
    adj <- lapply(seq_len(n), function(i) {
        as.integer(c(if (i > 1L) i - 1L, if (i < n) i + 1L))
    })
    weight <- lapply(seq_len(n), function(i) {
        as.numeric(c(if (i > 1L) weights[[i - 1L]], if (i < n) weights[[i]]))
    })
    list(adj = adj, weight = weight)
}

follow_test_roots <- function(next.edge) {
    vapply(seq_along(next.edge), function(start) {
        current <- start
        seen <- integer()
        while (!is.na(next.edge[[current]])) {
            if (current %in% seen) stop("cycle")
            seen <- c(seen, current)
            current <- next.edge[[current]]
        }
        current
    }, integer(1))
}

test_that("raw CLOSEST handles increasing and two-peak paths", {
    graph <- raw_path_graph(5)
    increasing <- raw_closest_args(graph$adj, graph$weight, 1:5)
    part <- closest.basin.partition(increasing, "max")
    expect_equal(part$next.edge, c(2L, 3L, 4L, 5L, NA_integer_))
    expect_equal(part$root.vertex, rep(5L, 5L))
    expect_equal(part$assignment, rep(1L, 5L))

    two.peak <- raw_closest_args(graph$adj, graph$weight, c(0, 2, 0, 3, 0))
    part <- closest.basin.partition(two.peak, "max")
    expect_equal(part$next.edge, c(2L, NA, 2L, NA, 4L))
    expect_equal(part$root.vertex, c(2L, 2L, 2L, 4L, 4L))
    expect_equal(part$assignment, c(1L, 1L, 1L, 2L, 2L))
})

test_that("connected maxima and interior shelves are contracted deterministically", {
    graph <- raw_path_graph(5)
    maximum <- raw_closest_args(graph$adj, graph$weight, c(0, 1, 3, 3, 1))
    expect_equal(maximum$plateau.component, c(1L, 2L, 3L, 3L, 4L))
    expect_equal(maximum$plateau.representative, c(1L, 2L, 3L, 5L))
    expect_equal(maximum$next.up, c(2L, 3L, NA, 3L, 4L))
    expect_equal(maximum$raw.max.root.vertex, rep(3L, 5L))

    shelf.graph <- raw_path_graph(4)
    shelf <- raw_closest_args(shelf.graph$adj, shelf.graph$weight, c(0, 2, 2, 3))
    expect_equal(shelf$plateau.component, c(1L, 2L, 2L, 3L))
    expect_equal(shelf$next.up, c(2L, 3L, 4L, NA))
    expect_equal(shelf$raw.max.root.vertex, rep(4L, 4L))
})

test_that("positive plateau tolerance contracts a valid near-equal component", {
    graph <- raw_path_graph(3)
    fit <- raw_closest_args(
        graph$adj, graph$weight, c(0, 0.04, 1),
        plateau.tolerance = 0.05
    )
    part <- closest.basin.partition(fit, "max")

    expect_equal(fit$plateau.component, c(1L, 1L, 2L))
    expect_equal(fit$plateau.representative, c(1L, 3L))
    expect_equal(fit$plateau.value, c(0.02, 1))
    expect_equal(part$next.edge, c(2L, 3L, NA_integer_))
    expect_equal(part$root.vertex, rep(3L, 3L))
    expect_equal(part$assignment, rep(1L, 3L))
})

test_that("CLOSEST uses length before steepness and stable IDs for ties", {
    adj <- list(c(3L, 2L), 1L, 1L)
    short <- raw_closest_args(adj, list(c(2, 1), 1, 2), c(0, 1, 10))
    expect_equal(short$next.up[[1]], 2L)

    tied <- raw_closest_args(adj, list(c(1, 1), 1, 1), c(0, 1, 10))
    expect_equal(tied$next.up[[1]], 2L)
})

test_that("raw CLOSEST partitions disconnected components completely", {
    adj <- list(2L, 1L, 4L, 3L)
    weight <- list(1, 1, 1, 1)
    fit <- raw_closest_args(adj, weight, c(0, 1, 2, 0))
    part <- closest.basin.partition(fit, "max")
    expect_equal(part$root.vertex, c(2L, 2L, 3L, 3L))
    expect_equal(part$assignment, c(1L, 1L, 2L, 2L))
    expect_false(anyNA(part$assignment))
})

test_that("raw CLOSEST is invariant to adjacency ordering", {
    adj.a <- list(c(2L, 3L), c(1L, 4L), c(1L, 4L), c(2L, 3L))
    wt.a <- list(c(1, 2), c(1, 1), c(2, 1), c(1, 1))
    adj.b <- lapply(adj.a, rev)
    wt.b <- lapply(wt.a, rev)
    y <- c(0, 1, 1, 2)
    a <- raw_closest_args(adj.a, wt.a, y)
    b <- raw_closest_args(adj.b, wt.b, y)
    expect_identical(a$plateau.component, b$plateau.component)
    expect_identical(a$next.up, b$next.up)
    expect_identical(a$raw.max.assignment, b$raw.max.assignment)
})

test_that("basin accessor agrees with direct root following", {
    graph <- raw_path_graph(7)
    fit <- raw_closest_args(graph$adj, graph$weight, c(0, 2, 2, 1, 3, 3, 0))
    max.part <- closest.basin.partition(fit, "max")
    min.part <- closest.basin.partition(fit, "min")
    expect_equal(max.part$root.vertex, follow_test_roots(fit$next.up))
    expect_equal(min.part$root.vertex, follow_test_roots(fit$next.down))
    expect_identical(max.part$assignment, fit$raw.max.assignment)
    expect_identical(min.part$assignment, fit$raw.min.assignment)
})

test_that("raw CLOSEST remains separate from postfiltered assignment", {
    graph <- raw_path_graph(5)
    fit <- compute.gfc.trajectory(
        graph$adj, graph$weight, c(0, 2, 0, 3, 0),
        modulation = "CLOSEST",
        edge.length.quantile.thld = 1,
        long.edge.fallback = "flag",
        plateau.tolerance = 0,
        min.basin.size = 10L,
        tie.breaking = FALSE
    )
    expect_false(anyNA(fit$raw.max.assignment))
    expect_true(anyNA(fit$max.assignment) ||
                !identical(fit$raw.max.assignment, fit$max.assignment))
    expect_identical(
        closest.basin.partition(fit, "max")$assignment,
        fit$raw.max.assignment
    )
})

test_that("CLOSEST rejects invalid graphs and scalar fields explicitly", {
    expect_error(
        raw_closest_args(list(2L, integer()), list(1, numeric()), c(0, 1)),
        "^invalid_graph:"
    )
    expect_error(
        raw_closest_args(list(1L), list(1), 0),
        "^invalid_graph:"
    )
    graph <- raw_path_graph(3)
    expect_error(
        raw_closest_args(graph$adj, graph$weight, c(0, NA, 1)),
        "^invalid_scalar_field:"
    )
    expect_error(
        raw_closest_args(
            graph$adj, graph$weight, c(0, 0.05, 0.1),
            plateau.tolerance = 0.06
        ),
        "^invalid_scalar_field: plateau tolerance chaining"
    )
    expect_error(
        raw_closest_args(graph$adj, graph$weight, c(0, 1, 2), tie.breaking = TRUE),
        "^unsupported_configuration:"
    )
    for (threshold in list(NA_real_, Inf, 0, 1.01, c(0.5, 0.9), "0.9")) {
        expect_error(
            raw_closest_args(
                graph$adj, graph$weight, c(0, 1, 2),
                edge.length.quantile.thld = threshold
            ),
            "^invalid_graph: edge.length.quantile.thld"
        )
    }
})

test_that("the direct gfc.flow API exposes both raw and legacy assignments", {
    graph <- raw_path_graph(4)
    fit <- raw_closest_args(graph$adj, graph$weight, c(0, 2, 2, 3))
    expect_s3_class(fit, "gfc.flow")
    expect_named(fit$raw.closest.forest, c(
        "method.id", "plateau.component", "plateau.representative",
        "plateau.value", "next.up", "next.down", "max.root.vertex",
        "min.root.vertex", "max.assignment", "min.assignment",
        "edge.length.quantile.thld", "edge.length.thld",
        "plateau.tolerance", "long.edge.fallback"
    ))
    expect_equal(fit$edge.length.thld, 1)
    expect_equal(fit$raw.closest.forest$method.id,
                 "closest_improving_neighbor_forest")
})

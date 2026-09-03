adaptive.extrema.chain <- function(...) {
    detect.adaptive.extrema(
        adj.list = list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L),
        weight.list = list(1, c(1, 1), c(1, 1), c(1, 1), 1),
        y = c(1, 3, 2, 5, 1),
        max.radius = 2,
        min.neighborhood.size = 2,
        ...
    )
}

test_that("adaptive detection preserves strict extrema and their neighborhoods", {
    maxima <- adaptive.extrema.chain()
    expect_identical(class(maxima), "gflow_local_extrema")
    expect_false(inherits(maxima, "local_extrema"))
    expect_identical(maxima$vertices, c(2L, 4L))
    expect_identical(maxima$values, c(3, 5))
    expect_identical(maxima$radii, c(1, 2))
    expect_identical(maxima$neighborhood_sizes, c(2L, 3L))
    expect_identical(maxima$labels, c("M2", "M1"))
    expect_identical(maxima$type, rep("Maximum", 2))
    expect_true(maxima$detect.maxima)
    expect_setequal(maxima$neighborhood_vertices[[1]], c(1L, 3L))
    expect_setequal(maxima$neighborhood_vertices[[2]], c(2L, 3L, 5L))
    expect_identical(maxima$graph_diameter, 4)

    minima <- adaptive.extrema.chain(detect.maxima = FALSE,
                                     custom.prefix = "low")
    expect_identical(minima$vertices, c(1L, 3L, 5L))
    expect_identical(minima$values, c(1, 2, 1))
    expect_identical(minima$radii, c(2, 1, 2))
    expect_identical(minima$neighborhood_sizes, rep(2L, 3))
    expect_identical(minima$labels, c("low1", "low3", "low2"))
    expect_identical(minima$type, rep("Minimum", 3))
    expect_false(minima$detect.maxima)

    for (x in list(maxima, minima)) {
        s <- summary(x)
        expect_identical(class(s), "summary.gflow_local_extrema")
        expect_identical(s$extrema_type, x$type[[1]])
        expect_identical(s$n_extrema, length(x$vertices))
        expect_identical(s$extrema_details$fn_value,
                         sort(x$values, decreasing = x$detect.maxima))
        expect_output(print(s), "Function values statistics:")
    }
})

test_that("empty results retain the requested extrema type", {
    for (detect.maxima in c(TRUE, FALSE)) {
        x <- detect.adaptive.extrema(
            list(2L, c(1L, 3L), 2L), list(1, c(1, 1), 1),
            rep(1, 3), max.radius = 2, min.neighborhood.size = 2,
            detect.maxima = detect.maxima
        )
        expected.type <- if (detect.maxima) "Maximum" else "Minimum"
        expect_identical(x$vertices, integer())
        expect_identical(x$is_maxima, logical())
        expect_identical(x$type, character())
        expect_identical(x$detect.maxima, detect.maxima)
        s <- summary(x)
        expect_identical(class(s), "summary.gflow_local_extrema")
        expect_identical(s$extrema_type, expected.type)
        expect_identical(s$n_extrema, 0L)
        expect_identical(s$extrema_details, data.frame())
        expect_output(print(s), paste("Extrema type:", expected.type))
        expect_message(plot(x), "No extrema to plot")
    }

    unknown <- structure(list(vertices = integer()),
                         class = "gflow_local_extrema")
    expect_identical(summary(unknown)$extrema_type, NA_character_)
    expect_error(adaptive.extrema.chain(detect.maxima = NA), "single logical")
})

test_that("vertex extraction includes the center exactly once when requested", {
    x <- adaptive.extrema.chain()
    original <- x
    for (i in seq_along(x$vertices)) {
        center <- x$vertices[[i]]
        neighbors <- x$neighborhood_vertices[[i]]
        label <- x$labels[[i]]
        expect_false(center %in% neighbors)
        expect_identical(gflow::vertices(x, label, include.center = FALSE),
                         neighbors)
        with.center <- gflow::vertices(x, label)
        expect_setequal(with.center, c(neighbors, center))
        expect_equal(sum(with.center == center), 1L)
        expect_length(with.center, x$neighborhood_sizes[[i]] + 1L)
        expect_identical(dgraphs::vertices(x, label), with.center)
    }
    expect_identical(x, original)

    x$neighborhood_vertices[[1]] <- c(1L, 2L, 3L, 2L)
    expect_identical(vertices(x, "M2"), c(1L, 2L, 3L))
    expect_identical(vertices(x, "M2", include.center = FALSE), c(1L, 3L))
    expect_error(vertices(x, "M2", include.center = NA), "single logical")
    expect_error(vertices(x, "missing"), "not found")
})

test_that("dgraphs objects retain their own detector and methods", {
    x <- dgraphs::detect.local.extrema(
        list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L),
        list(1, c(1, 1), c(1, 1), c(1, 1), 1),
        c(1, 3, 2, 5, 1), max.radius = 2, min.neighborhood.size = 2
    )
    expect_identical(class(x), "local_extrema")
    expect_identical(x$vertices, 4L)
    expect_identical(x$neighborhood_sizes, 4L)
    expect_identical(class(summary(x)), "summary.local_extrema")
    expect_setequal(gflow::vertices(x, "M1"), 2:5)
    expect_setequal(dgraphs::vertices(x, "M1", include.center = FALSE),
                    c(2L, 3L, 5L))
})

test_that("plot dispatch handles nonempty adaptive extrema", {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path)
    on.exit({
        grDevices::dev.off()
        unlink(path)
    }, add = TRUE)
    x <- adaptive.extrema.chain()
    expect_identical(plot(x), x)
})

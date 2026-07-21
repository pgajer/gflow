closest_test_args <- function(adj.list, weight.list, y, ...) {
    compute.gfc.trajectory(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 1,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        symmetric.seeding = FALSE,
        tie.breaking = FALSE,
        ...
    )
}

test_that("CLOSEST resolves equal edge lengths by target vertex ID", {
    adj.list <- list(
        c(3L, 2L),
        c(1L, 4L),
        c(1L, 5L),
        2L,
        3L
    )
    weight.list <- list(
        c(1, 1),
        c(1, 1),
        c(1, 1),
        1,
        1
    )
    y <- c(0, 1, 1.1, 3, 2.8)

    fit <- closest_test_args(adj.list, weight.list, y)

    expect_equal(fit$next.up, c(2L, 4L, 5L, NA, NA))
    expect_false(any(fit$long.edge.fallback$next.up.used))
    expect_identical(fit$long.edge.fallback$n.next.up, 0L)
})

test_that("CLOSEST long-edge fallback is selectable and auditable", {
    adj.list <- list(
        c(3L, 2L),
        1L,
        1L,
        5L,
        c(4L, 6L),
        5L
    )
    weight.list <- list(
        c(12, 10),
        10,
        12,
        1,
        c(1, 1),
        1
    )
    y <- c(0, 2, 3, 0, 1, 2)

    common <- list(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 0.25,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        symmetric.seeding = FALSE,
        tie.breaking = FALSE
    )

    flagged <- do.call(
        compute.gfc.trajectory,
        c(common, list(long.edge.fallback = "flag"))
    )
    allowed <- do.call(
        compute.gfc.trajectory,
        c(common, list(long.edge.fallback = "allow"))
    )
    forbidden <- do.call(
        compute.gfc.trajectory,
        c(common, list(long.edge.fallback = "forbid"))
    )

    expect_equal(flagged$next.up[1], 2L)
    expect_true(flagged$long.edge.fallback$next.up.used[1])
    expect_true(flagged$long.edge.fallback$attention.required)
    expect_gt(flagged$long.edge.fallback$n.next.up, 0L)
    expect_true(any(flagged$long.edge.fallback$trajectory.step.counts > 0L))
    expect_true(any(flagged$long.edge.fallback$max.basin.vertex.counts > 0L))

    expect_equal(allowed$next.up, flagged$next.up)
    expect_true(allowed$long.edge.fallback$next.up.used[1])
    expect_false(allowed$long.edge.fallback$attention.required)

    expect_true(is.na(forbidden$next.up[1]))
    expect_false(forbidden$long.edge.fallback$next.up.used[1])
    expect_true(forbidden$long.edge.fallback$next.up.blocked[1])
    expect_identical(forbidden$long.edge.fallback$n.next.up, 0L)
    expect_gt(forbidden$long.edge.fallback$n.next.up.blocked, 0L)
    expect_false(forbidden$long.edge.fallback$attention.required)
})

test_that("matrix CLOSEST flow forwards the long-edge fallback policy", {
    adj.list <- list(c(2L, 3L), 1L, 1L, 5L, 4L)
    weight.list <- list(c(10, 12), 10, 12, 1, 1)
    y <- c(0, 2, 3, 0, 1)
    Y <- cbind(y, y + 0.25)

    fits <- compute.gfc.flow.matrix(
        adj.list = adj.list,
        weight.list = weight.list,
        Y = Y,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 0.25,
        long.edge.fallback = "forbid",
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        symmetric.seeding = FALSE
    )

    expect_length(fits, 2L)
    expect_true(all(vapply(
        fits,
        function(fit) identical(fit$long.edge.fallback$policy, "forbid"),
        logical(1)
    )))
    expect_true(all(vapply(fits, function(fit) is.na(fit$next.up[1]), logical(1))))
})

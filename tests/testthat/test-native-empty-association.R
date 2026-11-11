test_that("empty retained basin supports preserve result shapes", {
    adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    lengths <- lapply(adjacency, function(x) rep(1, length(x)))
    y <- c(0, 3, 1, 2, 0)
    basins <- create.basin.complex(
        adjacency, lengths, y, method = "trajectory_flow", direction = "both",
        simplify.params = list(support.filter = list(
            enabled = TRUE, min.basin.size = 100
        ))
    )
    membership <- gfassoc.membership(basins, support.stage = "retained")
    expect_identical(membership$n_max_basins, 0L)
    expect_identical(membership$n_min_basins, 0L)
    expect_identical(membership$max_values, numeric())
    expect_identical(membership$min_values, numeric())
    for (name in c("max_membership", "min_membership", "cell_membership")) {
        expect_identical(membership[[name]], rep(list(numeric()), 5))
    }
    for (cells in membership$cell_indices) {
        expect_identical(dim(cells), c(0L, 2L))
    }
    polarity <- gfassoc.polarity(y, membership)
    expect_false(any(polarity$is_valid))
    overlap <- gfassoc.overlap(membership, membership)
    for (name in c("O_pp", "O_mm", "O_pm", "O_mp")) {
        expect_identical(dim(overlap[[name]]), c(0L, 0L))
        expect_length(overlap[[name]], 0L)
    }
    deviation <- gfassoc.deviation(overlap$O_pp)
    for (name in c("delta", "zeta", "expected")) {
        expect_identical(dim(deviation[[name]]), c(0L, 0L))
    }
})

test_that("partially uncovered supports retain nonempty membership weights", {
    adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    basins <- create.basin.complex(
        adjacency, lapply(adjacency, function(x) rep(1, length(x))),
        c(0, 3, 1, 2, 0), method = "trajectory_flow", direction = "both",
        simplify.params = list(support.filter = list(
            enabled = TRUE, min.basin.size = 3
        ))
    )
    membership <- gfassoc.membership(basins, support.stage = "retained")
    expect_identical(lengths(membership$min_membership), c(0L, 1L, 1L, 1L, 0L))
    expect_identical(membership$max_membership[[3L]], c(0.5, 0.5))
    expect_identical(membership$cell_membership[[1L]], numeric())
    expect_identical(membership$cell_membership[[5L]], numeric())
    expect_equal(sum(membership$cell_membership[[3L]]), 1)
})

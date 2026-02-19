make_ring_graph <- function(n, offsets = c(-2L, -1L, 1L, 2L), base_weight = 1) {
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    for (i in seq_len(n)) {
        nbrs <- unique(as.integer(((i - 1L + offsets) %% n) + 1L))
        adj.list[[i]] <- nbrs
        weight.list[[i]] <- rep(as.numeric(base_weight), length(nbrs))
    }
    list(adj = adj.list, w = weight.list)
}

test_that("RTCB basin constructor returns expected structure and invariants", {
    n <- 80L
    g <- make_ring_graph(n)
    x <- seq_len(n)
    y <- sin(2 * pi * x / n) + 0.2 * sin(8 * pi * x / n)

    basins <- gflow:::compute.basins.of.attraction.rtcb(
        adj.list = g$adj,
        weight.list = g$w,
        y = y,
        edge.length.quantile.thld = 1.0,
        with.trajectories = TRUE,
        n.min = 30L,
        q.min = 0.2,
        run.max = 2L,
        tau0 = 0.1,
        kappa = 1.5,
        k.max = 24L,
        h.max = 200L
    )

    expect_s3_class(basins, "basins_of_attraction")
    expect_equal(basins$n_vertices, n)
    expect_equal(length(basins$y), n)
    expect_true(length(basins$lmin_basins) >= 1L)
    expect_true(length(basins$lmax_basins) >= 1L)

    all_basins <- c(basins$lmin_basins, basins$lmax_basins)
    for (b in all_basins) {
        expect_true(is.numeric(b$vertex))
        expect_true(b$vertex >= 1 && b$vertex <= n)
        expect_true(is.matrix(b$basin_df))
        expect_equal(ncol(b$basin_df), 2L)
        expect_true(all(b$basin_df[, 1] >= 1 & b$basin_df[, 1] <= n))
        expect_true(all(b$basin_df[, 2] >= 0))
        expect_true(is.vector(b$terminal_extrema))
        expect_true(all(b$terminal_extrema >= 1 & b$terminal_extrema <= n))
        expect_equal(length(b$all_predecessors), n)
        expect_true(is.list(b$trajectory_sets))
    }
})

test_that("compute.gfc supports modulation = RTCB and records parameters", {
    n <- 90L
    g <- make_ring_graph(n)
    x <- seq_len(n)
    y <- sin(2 * pi * x / n) + 0.25 * cos(6 * pi * x / n)

    fit <- compute.gfc(
        adj.list = g$adj,
        edge.length.list = g$w,
        fitted.values = y,
        modulation = "RTCB",
        rtcb.params = list(
            n.min = 35L,
            q.min = 0.2,
            run.max = 2L,
            tau0 = 0.1
        ),
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        expand.basins = FALSE,
        with.trajectories = FALSE,
        verbose = FALSE
    )

    expect_s3_class(fit, "basins_of_attraction")
    expect_s3_class(fit$summary, "basin_summary")
    expect_true(nrow(fit$summary) >= 1L)
    expect_true(all(fit$summary$type %in% c("min", "max")))
    expect_identical(fit$parameters$modulation, "RTCB")
    expect_true(is.list(fit$parameters$rtcb.params))
    expect_true(nrow(fit$stage.history) >= 1L)
})

test_that("compute.gfc keeps trajectory-first behavior for non-RTCB modulations", {
    n <- 60L
    g <- make_ring_graph(n)
    y <- sin(seq_len(n) / 5)

    fit <- compute.gfc(
        adj.list = g$adj,
        edge.length.list = g$w,
        fitted.values = y,
        modulation = "CLOSEST",
        verbose = FALSE
    )

    expect_s3_class(fit, "gfc.flow")
})

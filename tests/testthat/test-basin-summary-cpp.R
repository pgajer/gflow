test_that("summary.basins_of_attraction C++ path matches R fallback", {
    set.seed(42)
    n <- 80L

    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    for (i in seq_len(n)) {
        nbrs <- c(
            ((i - 2L) %% n) + 1L,
            (i %% n) + 1L,
            ((i - 3L) %% n) + 1L,
            ((i + 1L) %% n) + 1L
        )
        nbrs <- unique(as.integer(nbrs))
        adj.list[[i]] <- nbrs

        # Positive, non-uniform edge lengths.
        w <- abs(rnorm(length(nbrs), mean = 1.0, sd = 0.25)) + 1e-3
        weight.list[[i]] <- as.numeric(w)
    }

    y <- rnorm(n)
    basins <- compute.basins.of.attraction(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        edge.length.quantile.thld = 1.0,
        with.trajectories = FALSE
    )

    s_r <- summary(
        basins,
        adj.list,
        weight.list,
        hop.k = 2L,
        use.cpp = FALSE
    )

    s_cpp <- summary(
        basins,
        adj.list,
        weight.list,
        hop.k = 2L,
        use.cpp = TRUE
    )

    expect_equal(names(s_cpp), names(s_r))
    expect_equal(nrow(s_cpp), nrow(s_r))
    expect_equal(s_cpp$label, s_r$label)
    expect_equal(s_cpp$type, s_r$type)
    expect_equal(s_cpp$vertex, s_r$vertex)
    expect_equal(s_cpp$value, s_r$value, tolerance = 1e-12)
    expect_equal(s_cpp$rel.value, s_r$rel.value, tolerance = 1e-12)
    expect_equal(s_cpp$hop.idx, s_r$hop.idx, tolerance = 1e-12)
    expect_equal(s_cpp$basin.size, s_r$basin.size)
    expect_equal(s_cpp$p.mean.nbrs.dist, s_r$p.mean.nbrs.dist, tolerance = 1e-12)
    expect_equal(s_cpp$p.mean.hopk.dist, s_r$p.mean.hopk.dist, tolerance = 1e-12)
    expect_equal(s_cpp$deg, s_r$deg, tolerance = 1e-12)
    expect_equal(s_cpp$p.deg, s_r$p.deg, tolerance = 1e-12)
})

test_that("summary.basins_of_attraction supports cached vertex metrics", {
    set.seed(7)
    n <- 60L

    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    for (i in seq_len(n)) {
        nbrs <- c(
            ((i - 2L) %% n) + 1L,
            (i %% n) + 1L,
            ((i + 1L) %% n) + 1L
        )
        nbrs <- unique(as.integer(nbrs))
        adj.list[[i]] <- nbrs
        weight.list[[i]] <- as.numeric(abs(rnorm(length(nbrs), 1, 0.2)) + 1e-3)
    }

    y <- rnorm(n)
    basins <- compute.basins.of.attraction(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        edge.length.quantile.thld = 1.0,
        with.trajectories = FALSE
    )

    adj.0 <- lapply(adj.list, function(x) as.integer(x - 1L))
    metrics.k2 <- .Call(
        "S_precompute_basin_vertex_metrics",
        adj.0,
        weight.list,
        as.integer(2L),
        PACKAGE = "gflow"
    )

    s_direct <- summary(
        basins,
        adj.list,
        weight.list,
        hop.k = 2L,
        use.cpp = TRUE
    )

    s_cached <- summary(
        basins,
        adj.list,
        weight.list,
        hop.k = 2L,
        vertex.metrics = metrics.k2,
        use.cpp = TRUE
    )

    expect_equal(s_cached, s_direct, tolerance = 1e-12)

    metrics.k3 <- .Call(
        "S_precompute_basin_vertex_metrics",
        adj.0,
        weight.list,
        as.integer(3L),
        PACKAGE = "gflow"
    )

    expect_error(
        summary(
            basins,
            adj.list,
            weight.list,
            hop.k = 2L,
            vertex.metrics = metrics.k3,
            use.cpp = TRUE
        ),
        "hop\\.k"
    )
})

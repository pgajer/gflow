test_that("compute.gfc handles empty geometric retention without error", {
    set.seed(11)
    n <- 70L

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
        weight.list[[i]] <- as.numeric(abs(rnorm(length(nbrs), 1, 0.2)) + 1e-3)
    }

    y <- rnorm(n)

    res <- compute.gfc(
        adj.list = adj.list,
        edge.length.list = weight.list,
        fitted.values = y,
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = TRUE,
        p.mean.nbrs.dist.threshold = 0.0,
        p.mean.hopk.dist.threshold = 0.0,
        min.basin.size = n + 1L,
        expand.basins = FALSE,
        verbose = FALSE
    )

    expect_s3_class(res$summary, "basin_summary")
    expect_equal(nrow(res$summary), 0L)
    expect_length(res$basins$lmin_basins, 0L)
    expect_length(res$basins$lmax_basins, 0L)
})

test_that("compute.gfc geometric filtering ignores p.deg.threshold", {
    set.seed(12)
    n <- 80L

    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    for (i in seq_len(n)) {
        nbrs <- c(
            ((i - 2L) %% n) + 1L,
            (i %% n) + 1L,
            ((i - 4L) %% n) + 1L,
            ((i + 2L) %% n) + 1L
        )
        nbrs <- unique(as.integer(nbrs))
        adj.list[[i]] <- nbrs
        weight.list[[i]] <- as.numeric(abs(rnorm(length(nbrs), 1, 0.25)) + 1e-3)
    }

    y <- rnorm(n)

    set.seed(101)
    res.low <- compute.gfc(
        adj.list = adj.list,
        edge.length.list = weight.list,
        fitted.values = y,
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = TRUE,
        p.mean.nbrs.dist.threshold = 0.95,
        p.mean.hopk.dist.threshold = 0.95,
        p.deg.threshold = 0.01,
        min.basin.size = 3L,
        expand.basins = FALSE,
        verbose = FALSE
    )

    set.seed(101)
    res.high <- compute.gfc(
        adj.list = adj.list,
        edge.length.list = weight.list,
        fitted.values = y,
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = TRUE,
        p.mean.nbrs.dist.threshold = 0.95,
        p.mean.hopk.dist.threshold = 0.95,
        p.deg.threshold = 0.99,
        min.basin.size = 3L,
        expand.basins = FALSE,
        verbose = FALSE
    )

    expect_equal(res.low$summary, res.high$summary, tolerance = 1e-12)
    expect_equal(res.low$stage.history, res.high$stage.history, tolerance = 1e-12)
})

test_that("compute.gfc GEODESIC summary uses harmonized schema", {
    n <- 64L
    graph <- generate.circle.graph(n, type = "uniform")
    y <- 2 + sin(seq(0, 2 * pi, length.out = n)) + 0.05 * cos(seq(0, 6 * pi, length.out = n))

    res <- compute.gfc(
        adj.list = graph$adj.list,
        edge.length.list = graph$weight.list,
        fitted.values = y,
        modulation = "GEODESIC",
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        expand.basins = TRUE,
        verbose = FALSE
    )

    required.cols <- c(
        "label", "vertex", "value", "rel.value", "type",
        "basin.size.raw", "basin.size.exp", "basin.size",
        "hop.index", "hop.idx",
        "p.mean.nbrs.dist", "p.mean.hopk.dist",
        "degree", "deg.percentile", "deg", "p.deg"
    )

    expect_true(all(required.cols %in% names(res$summary)))
    expect_equal(as.integer(res$summary$basin.size.raw), as.integer(res$summary$basin.size))
    expect_equal(res$summary$hop.index, res$summary$hop.idx, tolerance = 1e-12)
    expect_equal(res$summary$degree, res$summary$deg, tolerance = 1e-12)
    expect_equal(res$summary$deg.percentile, res$summary$p.deg, tolerance = 1e-12)
    expect_equal(
        res$summary$rel.value,
        res$summary$value / mean(res$basins$y),
        tolerance = 1e-12
    )

    exp.size.map <- numeric(0)
    if (!is.null(res$expanded.max.vertices.list) && length(res$expanded.max.vertices.list) > 0) {
        exp.size.map <- c(
            exp.size.map,
            setNames(
                vapply(res$expanded.max.vertices.list, length, integer(1)),
                names(res$expanded.max.vertices.list)
            )
        )
    }
    if (!is.null(res$expanded.min.vertices.list) && length(res$expanded.min.vertices.list) > 0) {
        exp.size.map <- c(
            exp.size.map,
            setNames(
                vapply(res$expanded.min.vertices.list, length, integer(1)),
                names(res$expanded.min.vertices.list)
            )
        )
    }

    mapped.exp <- exp.size.map[as.character(res$summary$label)]
    idx <- !is.na(mapped.exp)
    expect_true(any(idx))
    expect_equal(as.integer(res$summary$basin.size.exp[idx]), as.integer(mapped.exp[idx]))
})

test_that("compute.gfc trajectory modulation summary uses harmonized schema", {
    n <- 64L
    graph <- generate.circle.graph(n, type = "uniform")
    y <- 2 + sin(seq(0, 2 * pi, length.out = n)) + 0.07 * sin(seq(0, 8 * pi, length.out = n))

    res <- compute.gfc(
        adj.list = graph$adj.list,
        edge.length.list = graph$weight.list,
        fitted.values = y,
        modulation = "CLOSEST",
        edge.length.quantile.thld = 1.0,
        apply.relvalue.filter = FALSE,
        apply.maxima.clustering = FALSE,
        apply.minima.clustering = FALSE,
        apply.geometric.filter = FALSE,
        min.basin.size = 1L,
        verbose = FALSE
    )

    required.cols <- c(
        "label", "vertex", "value", "rel.value", "type",
        "basin.size.raw", "basin.size.exp", "basin.size",
        "hop.index", "hop.idx",
        "p.mean.nbrs.dist", "p.mean.hopk.dist",
        "degree", "deg.percentile", "deg", "p.deg"
    )

    expect_s3_class(res, "gfc.flow")
    expect_true(all(required.cols %in% names(res$summary)))
    expect_equal(as.integer(res$summary$basin.size.raw), as.integer(res$summary$basin.size))
    expect_equal(as.integer(res$summary$basin.size.exp), as.integer(res$summary$basin.size.raw))
    expect_equal(res$summary$hop.index, res$summary$hop.idx, tolerance = 1e-12)
    expect_equal(res$summary$degree, res$summary$deg, tolerance = 1e-12)
    expect_equal(res$summary$deg.percentile, res$summary$p.deg, tolerance = 1e-12)
    expect_equal(
        res$summary$rel.value,
        res$summary$value / mean(res$y),
        tolerance = 1e-12
    )
})

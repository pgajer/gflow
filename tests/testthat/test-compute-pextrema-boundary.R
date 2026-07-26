test_that("probabilistic extrema remain owned by gflow", {
    fit <- structure(
        list(
            optimal.fit = list(
                y = c(1, 2, 3),
                fitted.values = c(1, 2, 3),
                graph = list(
                    adj.list = list(2L, c(1L, 3L), 2L),
                    vertex.densities = c(1, 1, 1),
                    edge.densities = c(1, 1, 1, 1),
                    edge.length.list = list(1, c(1, 1), 1)
                )
            )
        ),
        class = "riem.dcx"
    )

    extrema <- compute.pextrema.nbhds(
        fit,
        extremp_quantile = 0,
        eff_degree_quantile = 0
    )

    expect_s3_class(extrema, "data.frame")
    expect_true(is.infinite(
        extrema$hop_extremp_radius[
            extrema$vertex == 3 & extrema$type == "max"
        ]
    ))
    expect_true(is.infinite(
        extrema$hop_extremp_radius[
            extrema$vertex == 1 & extrema$type == "min"
        ]
    ))
})

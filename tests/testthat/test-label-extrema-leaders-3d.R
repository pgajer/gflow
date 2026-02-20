test_that("label.extrema.leaders.3d validates required extrema columns", {
    S <- matrix(
        c(
            0, 0, 0,
            1, 0, 0,
            0, 1, 0
        ),
        ncol = 3,
        byrow = TRUE
    )

    bad.df <- data.frame(
        vertex = c(1L, 2L),
        label = c("M1", "m1"),
        stringsAsFactors = FALSE
    )

    expect_error(
        label.extrema.leaders.3d(S, bad.df),
        "extrema.df must have either 'is_max' column \\(0/1\\) or 'type' column \\('min'/'max'\\)"
    )
})

test_that("label.extrema.leaders.3d warns when no extrema remain after filtering", {
    skip_if_not_installed("rgl")

    old.opt <- options(rgl.useNULL = TRUE)
    on.exit(options(old.opt), add = TRUE)
    rgl::open3d(useNULL = TRUE)
    on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

    S <- matrix(
        c(
            0, 0, 0,
            1, 0, 0,
            0, 1, 0
        ),
        ncol = 3,
        byrow = TRUE
    )

    extrema.df <- data.frame(
        vertex = c(1L, 2L),
        label = c("m1", "m2"),
        type = c("min", "min"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        label.extrema.leaders.3d(S, extrema.df, extrema.type = "maxima"),
        "No extrema to plot"
    )
})

test_that("label.extrema.leaders.3d runs with rowname-based vertices and avoidance", {
    skip_if_not_installed("rgl")

    old.opt <- options(rgl.useNULL = TRUE)
    on.exit(options(old.opt), add = TRUE)
    rgl::open3d(useNULL = TRUE)
    on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

    set.seed(1)
    S <- matrix(rnorm(18), ncol = 3)
    rownames(S) <- paste0("v", seq_len(nrow(S)))

    extrema.df <- data.frame(
        vertex = c("v1", "v3", "v5"),
        label = c("M1", "m1", "M2"),
        type = c("max", "min", "max"),
        stringsAsFactors = FALSE
    )

    expect_invisible(
        label.extrema.leaders.3d(
            graph.3d = S,
            extrema.df = extrema.df,
            extrema.type = "both",
            offset = c(0, -0.2, 0),
            radial.frac = 0.12,
            z.frac = 0.015,
            avoid.vertices = c("v2", "v4"),
            avoid.weight = 1.2,
            lab.cex = 1.1
        )
    )
})

test_that("label.extrema.leaders.3d warns on inconsistent type/is_max and still runs", {
    skip_if_not_installed("rgl")

    old.opt <- options(rgl.useNULL = TRUE)
    on.exit(options(old.opt), add = TRUE)
    rgl::open3d(useNULL = TRUE)
    on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

    S <- matrix(
        c(
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        ),
        ncol = 3,
        byrow = TRUE
    )

    extrema.df <- data.frame(
        vertex = c(1L, 2L),
        label = c("M1", "m1"),
        type = c("max", "min"),
        is_max = c(0L, 0L),
        stringsAsFactors = FALSE
    )

    expect_warning(
        expect_invisible(label.extrema.leaders.3d(S, extrema.df)),
        "inconsistent values; using 'is_max'"
    )
})

test_that("shared local PCA chart matches anchor-centered R SVD contract", {
    set.seed(301)
    X <- matrix(stats::rnorm(11 * 4), nrow = 11)
    center <- X[3, ] + c(0.2, -0.1, 0.05, 0)
    chart.dim <- 3L

    centered <- sweep(X, 2L, center, "-")
    sv <- svd(centered, nu = 0L, nv = chart.dim)
    coords.r <- centered %*% sv$v[, seq_len(chart.dim), drop = FALSE]

    chart <- rcpp_local_pca_chart(
        X_support = X,
        center = center,
        chart_dim = chart.dim,
        center_mode = "anchor",
        dim_rule = "fixed",
        rebase_to_anchor = TRUE,
        orient_basis = FALSE
    )

    expect_equal(chart$singular.values, sv$d, tolerance = 1e-10)
    expect_equal(
        tcrossprod(chart$coordinates),
        tcrossprod(coords.r),
        tolerance = 1e-10
    )
})

test_that("shared local PCA chart matches SSRHE mean-centered rebasing", {
    set.seed(302)
    X <- cbind(
        stats::rnorm(12),
        stats::rnorm(12),
        0.25 * stats::rnorm(12)
    )
    anchor <- X[5, ]
    chart.dim <- 2L

    centered.mean <- sweep(X, 2L, colMeans(X), "-")
    sv <- svd(centered.mean, nu = 0L, nv = chart.dim)
    coords.r <- centered.mean %*% sv$v[, seq_len(chart.dim), drop = FALSE]
    coords.r <- sweep(coords.r, 2L, coords.r[5, ], "-")

    chart <- rcpp_local_pca_chart(
        X_support = X,
        center = anchor,
        chart_dim = chart.dim,
        center_mode = "mean",
        dim_rule = "fixed",
        rebase_to_anchor = TRUE,
        orient_basis = FALSE
    )

    expect_equal(chart$singular.values, sv$d, tolerance = 1e-10)
    expect_equal(
        tcrossprod(chart$coordinates),
        tcrossprod(coords.r),
        tolerance = 1e-10
    )
})

test_that("shared local PCA chart supports cumulative dimension and weights", {
    set.seed(303)
    t <- seq(-1, 1, length.out = 18)
    X <- cbind(t, 0.02 * t^2, 0.001 * stats::rnorm(length(t)))
    w <- seq(0.5, 2, length.out = length(t))

    chart <- rcpp_local_pca_chart(
        X_support = X,
        center = X[9, ],
        chart_dim = 0L,
        center_mode = "anchor",
        dim_rule = "eigen.cumulative",
        eigen_tolerance = 0.95,
        weights = w,
        rebase_to_anchor = TRUE,
        orient_basis = TRUE
    )

    expect_true(chart$chart.dim >= 1L)
    expect_true(chart$chart.dim <= ncol(X))
    expect_equal(ncol(chart$coordinates), chart$chart.dim)
    expect_equal(ncol(chart$basis), chart$chart.dim)
    expect_true(is.finite(chart$selected.variance.ratio))
})

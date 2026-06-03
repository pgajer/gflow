test_that("kernel.local.polynomial.cv fits ambient-coordinate local polynomials", {
    set.seed(1)
    x <- seq(0, 1, length.out = 24)
    X <- cbind(x, x^2)
    y <- 1 + 2 * x - 0.5 * x^2
    foldid <- rep(1:4, length.out = length(y))

    fit <- kernel.local.polynomial.cv(
        X, y,
        foldid = foldid,
        support.grid = c(8L, 12L),
        degree.grid = 1:2,
        kernel.grid = c("gaussian", "tricube"),
        coordinate.method = "coordinates"
    )

    expect_s3_class(fit, "kernel.local.polynomial.cv")
    expect_equal(fit$method.id, "kernel_local_polynomial_cv")
    expect_equal(fit$coordinate.method, "coordinates")
    expect_equal(nrow(fit$cv.table), 8L)
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(
        as.numeric(predict(fit)),
        as.numeric(fit$fitted.values),
        tolerance = 1e-12
    )
})

test_that("kernel.local.polynomial.cv supports local PCA auto chart dimension", {
    t <- seq(-1, 1, length.out = 22)
    X <- cbind(t, t^2, 0.001 * sin(seq_along(t)))
    y <- sin(t) + 0.25 * t
    foldid <- rep(1:5, length.out = length(y))

    fit <- kernel.local.polynomial.cv(
        X, y,
        foldid = foldid,
        support.grid = c(8L, 10L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = "auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator"
    )

    expect_s3_class(fit, "kernel.local.polynomial.cv")
    expect_equal(fit$coordinate.method, "local.pca")
    expect_true(isTRUE(fit$auto.chart.dim))
    expect_true(fit$chart.dim >= 1L)
    expect_true(fit$chart.dim < ncol(X))
    expect_equal(fit$auto.chart.support.metric, "both")
    expect_equal(fit$auto.chart.selection.metric, "operator")
    expect_equal(
        fit$auto.chart.dim.diagnostics$summary$support.metric,
        "coordinates"
    )
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("kernel.local.polynomial.cv handles underdetermined ambient designs", {
    set.seed(44)
    X <- matrix(stats::rnorm(36 * 40), nrow = 36)
    y <- sin(X[, 1]) + stats::rnorm(36, sd = 0.02)
    foldid <- rep(1:4, length.out = length(y))

    fit <- kernel.local.polynomial.cv(
        X, y,
        foldid = foldid,
        support.grid = 12L,
        degree.grid = 2L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates"
    )

    expect_s3_class(fit, "kernel.local.polynomial.cv")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(is.finite(fit$cv.table$cv.rmse.observed)))
})

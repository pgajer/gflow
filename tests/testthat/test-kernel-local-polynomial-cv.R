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

test_that("cached kernel neighborhoods match direct neighborhoods", {
    set.seed(87)
    X <- matrix(stats::rnorm(16 * 4), nrow = 16)
    y <- stats::rnorm(16)
    center <- X[3, ] + 0.1

    direct <- .klp.local.neighborhood(
        X.train = X,
        y.train = y,
        center = center,
        support.size = 9L,
        coordinate.method = "coordinates",
        chart.dim = ncol(X)
    )
    ordered <- .klp.local.order(
        X.train = X,
        center = center,
        support.size = 12L
    )
    cached <- .klp.local.neighborhood.from.order(
        X.train = X,
        y.train = y,
        center = center,
        ordered = ordered,
        support.size = 9L,
        coordinate.method = "coordinates",
        chart.dim = ncol(X)
    )

    expect_equal(cached$index, direct$index)
    expect_equal(cached$distances, direct$distances)
    expect_equal(cached$y, direct$y)
    expect_equal(cached$z, direct$z)
})

test_that("cached kernel designs match direct local polynomial designs", {
    set.seed(88)
    z <- matrix(stats::rnorm(10 * 4), nrow = 10)
    y <- stats::rnorm(10)
    weights <- stats::runif(10, min = 0.2, max = 1)
    cand <- expand.grid(
        support.size = 10L,
        degree = 1:2,
        kernel = c("gaussian", "tricube"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    cand$chart.dim <- c(3L, 3L, 3L, 3L)

    cache <- .klp.local.design.cache(z, cand, seq_len(nrow(cand)))
    key <- .klp.design.cache.key(2L, 3L)
    direct.design <- .malps.design.matrix(z[, seq_len(3L), drop = FALSE], 2L)

    expect_equal(cache[[key]], direct.design)
    expect_equal(
        .klp.fit.intercept.design(cache[[key]], y, weights),
        .klp.fit.intercept(z[, seq_len(3L), drop = FALSE], y, weights, 2L),
        tolerance = 1e-12
    )
})

test_that("kernel design feasibility precheck skips impossible designs", {
    set.seed(89)
    z <- matrix(stats::rnorm(8 * 40), nrow = 8)
    y <- stats::rnorm(8)
    weights <- stats::runif(8, min = 0.2, max = 1)
    cache <- new.env(parent = emptyenv())
    key <- .klp.design.cache.key(2L, 40L)

    fit <- .klp.fit.intercept.lazy(
        z = z,
        y = y,
        weights = weights,
        degree = 2L,
        chart.dim = 40L,
        design.cache = cache
    )

    expect_equal(.klp.design.ncol(2L, 40L), 861L)
    expect_false(exists(key, envir = cache, inherits = FALSE))
    expect_equal(fit, stats::weighted.mean(y, weights), tolerance = 1e-12)
})

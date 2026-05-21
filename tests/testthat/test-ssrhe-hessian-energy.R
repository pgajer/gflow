test_that("ssrhe.hessian.operator constructs an auditable sparse Hessian operator", {
    skip_if_not_installed("Matrix")

    grid <- expand.grid(x = seq(0, 1, length.out = 5),
                        y = seq(0, 1, length.out = 5))
    X <- as.matrix(grid)
    op <- ssrhe.hessian.operator(
        X = X,
        k = 12L,
        tangent.dim = 2L,
        stabilizer = TRUE
    )

    expect_s3_class(op, "ssrhe.hessian.operator")
    expect_true(inherits(op$A, "sparseMatrix"))
    expect_true(inherits(op$B, "sparseMatrix"))
    expect_true(inherits(op$BS, "sparseMatrix"))
    expect_equal(dim(op$A), c(nrow(X) * 3L, nrow(X)))
    expect_equal(dim(op$B), c(nrow(X), nrow(X)))
    expect_equal(op$nn.index[, 1], seq_len(nrow(X)))
    expect_equal(as.matrix(op$B), as.matrix(Matrix::crossprod(op$A)),
                 tolerance = 1e-10)
    expect_true(max(abs(as.matrix(op$B - Matrix::t(op$B)))) < 1e-10)

    expect_true(all(op$diagnostics$design.rank == op$diagnostics$design.ncol))
    expect_true(all(op$row.table$diagonal %in% c(0L, 1L)))
    expect_equal(unique(op$row.table$scale[op$row.table$diagonal == 1L]),
                 sqrt(2), tolerance = 1e-12)
})

test_that("ssrhe.hessian.operator annihilates affine functions but not quadratics", {
    skip_if_not_installed("Matrix")

    grid <- expand.grid(x = seq(0, 1, length.out = 6),
                        y = seq(0, 1, length.out = 6))
    X <- as.matrix(grid)
    op <- ssrhe.hessian.operator(
        X = X,
        k = 14L,
        tangent.dim = 2L,
        stabilizer = TRUE
    )

    one <- rep(1, nrow(X))
    linear <- 1 + 2 * X[, 1] - 0.5 * X[, 2]
    quadratic <- X[, 1]^2 + 0.25 * X[, 1] * X[, 2] - X[, 2]^2

    expect_lt(sqrt(sum((op$A %*% one)^2)), 1e-8)
    expect_lt(sqrt(sum((op$A %*% linear)^2)), 1e-8)
    expect_gt(sqrt(sum((op$A %*% quadratic)^2)), 1e-3)

    expect_lt(sqrt(sum((op$BS %*% one)^2)), 1e-8)
    expect_lt(sqrt(sum((op$BS %*% linear)^2)), 1e-7)
    expect_lt(sqrt(sum((op$BS %*% quadratic)^2)), 1e-6)
})

test_that("ssrhe.hessian.operator supports supplied neighbor indices", {
    skip_if_not_installed("Matrix")

    X <- cbind(seq(0, 1, length.out = 8), 0)
    nn <- matrix(NA_integer_, nrow = nrow(X), ncol = 4L)
    for (i in seq_len(nrow(X))) {
        ord <- order(abs(X[, 1] - X[i, 1]), seq_len(nrow(X)))
        nn[i, ] <- ord[seq_len(ncol(nn))]
    }

    op <- ssrhe.hessian.operator(
        X = X,
        k = 4L,
        tangent.dim = 1L,
        nn.index = nn,
        return.BS = FALSE
    )

    expect_equal(op$nn.index, nn)
    expect_equal(dim(op$A), c(nrow(X), nrow(X)))
    expect_equal(op$row.table$a, rep(1L, nrow(X)))
    expect_null(op$BS)
})

test_that("ssrhe.hessian.operator supports supplied variable-size supports", {
    skip_if_not_installed("Matrix")

    X <- cbind(seq(0, 1, length.out = 9), 0)
    support <- lapply(seq_len(nrow(X)), function(i) {
        ord <- order(abs(X[, 1] - X[i, 1]), seq_len(nrow(X)))
        ord[seq_len(min(nrow(X), 3L + (i %% 3L)))]
    })
    support <- Map(function(ids, i) unique(c(i, ids)), support, seq_along(support))

    op <- ssrhe.hessian.operator(
        X = X,
        tangent.dim = 1L,
        neighborhood.type = "supplied",
        support.index = support,
        return.BS = FALSE
    )

    expect_null(op$nn.index)
    expect_equal(op$neighborhoods$type, "supplied")
    expect_equal(op$diagnostics$k, lengths(support))
    expect_equal(length(op$support.index), nrow(X))
    expect_equal(dim(op$A), c(nrow(X), nrow(X)))
})

test_that("ssrhe.hessian.operator builds adaptive-radius local supports", {
    skip_if_not_installed("Matrix")

    set.seed(42)
    X <- cbind(seq(0, 1, length.out = 25),
               0.1 * sin(seq(0, 2 * pi, length.out = 25)))

    op <- ssrhe.hessian.operator(
        X = X,
        tangent.dim = 1L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 2L,
        radius.rule = "geomean",
        radius.factor = 1.05,
        min.support = 5L,
        max.support = 7L,
        support.topup = "nearest",
        return.BS = FALSE
    )

    expect_equal(op$parameters$neighborhood.type, "adaptive.radius")
    expect_equal(op$neighborhoods$adaptive.k.scale, 2L)
    expect_true(all(op$neighborhoods$support.size >= 5L))
    expect_true(all(op$neighborhoods$support.size <= 7L))
    expect_true(any(op$neighborhoods$n.topup > 0L))
    expect_equal(op$diagnostics$k, op$neighborhoods$support.size)
    expect_equal(dim(op$A), c(nrow(X), nrow(X)))
})

test_that("ssrhe.hessian.operator can choose local dimension by cumulative PCA variance", {
    skip_if_not_installed("Matrix")

    theta <- seq(0, 2 * pi, length.out = 20)
    X <- cbind(cos(theta), sin(theta))
    op <- ssrhe.hessian.operator(
        X = X,
        k = 8L,
        tangent.dim = NULL,
        tangent.dim.rule = "eigen.cumulative",
        eigen.tolerance = 0.8,
        return.BS = FALSE
    )

    expect_true(all(op$diagnostics$tangent.dim >= 1L))
    expect_true(all(op$diagnostics$tangent.dim <= 2L))
    expect_equal(ncol(op$B), nrow(X))
})

test_that("fit.ssrhe.hessian.regression matches dense fixed-lambda reference", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    lambda1 <- 0.3
    lambda2 <- 0.05

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda1 = lambda1,
        lambda2 = lambda2,
        stabilizer = TRUE
    )

    W <- diag(1, nrow(X))
    M <- W + lambda1 * as.matrix(fit$operator$B) +
        lambda2 * as.matrix(fit$operator$BS)
    ref <- as.vector(solve(M, y))

    expect_s3_class(fit, "ssrhe.hessian.fit")
    expect_equal(fit$fitted.values, ref, tolerance = 1e-9)
    expect_equal(fit$residuals, y - ref, tolerance = 1e-9)
    expect_equal(fit$energies$hessian,
                 as.numeric(crossprod(ref, as.matrix(fit$operator$B) %*% ref)),
                 tolerance = 1e-9)
})

test_that("fit.ssrhe.hessian.regression reproduces fully observed y when penalties are zero", {
    skip_if_not_installed("Matrix")

    X <- cbind(seq(0, 1, length.out = 8), 0)
    y <- cos(seq_along(X[, 1]))

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 4L,
        tangent.dim = 1L,
        lambda1 = 0,
        lambda2 = 0
    )

    expect_equal(fit$fitted.values, y, tolerance = 1e-12)
    expect_equal(fit$residuals, rep(0, length(y)), tolerance = 1e-12)
    expect_equal(fit$objective, 0, tolerance = 1e-12)
})

test_that("refit.ssrhe.hessian.regression reuses the operator for new responses and lambdas", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y1 <- X[, 1] + X[, 2]
    y2 <- cbind(a = sin(X[, 1]), b = cos(X[, 2]))

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y1,
        k = 10L,
        tangent.dim = 2L,
        lambda1 = 0.2,
        lambda2 = 0
    )
    refit <- refit.ssrhe.hessian.regression(
        fit,
        y.new = y2,
        lambda1 = 0.4
    )

    M <- diag(1, nrow(X)) + 0.4 * as.matrix(fit$operator$B)
    dense <- solve(M, y2)

    expect_s3_class(refit, "ssrhe.hessian.refit")
    expect_equal(refit$fitted.values, dense, tolerance = 1e-9)
    expect_equal(colnames(refit$fitted.values), colnames(y2))
    expect_identical(refit$operator, fit$operator)
})

test_that("fit.ssrhe.hessian.regression handles missing responses and observation weights", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- X[, 1]^2 - X[, 2]
    y[c(2, 7)] <- NA_real_
    weights <- rep(1, nrow(X))
    weights[c(5, 9)] <- 0.25

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda1 = 0.15,
        weights = weights,
        ridge = 1e-8
    )

    y.clean <- y
    y.clean[is.na(y.clean)] <- 0
    w.clean <- weights
    w.clean[is.na(y)] <- 0
    M <- diag(w.clean + 1e-8) + 0.15 * as.matrix(fit$operator$B)
    ref <- as.vector(solve(M, w.clean * y.clean))

    expect_equal(fit$weights, w.clean)
    expect_equal(fit$fitted.values, ref, tolerance = 1e-8)
    expect_true(all(is.na(fit$residuals[c(2, 7)])))
    expect_false(any(is.na(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.regression matches SSRHE semi-supervised label convention", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 4)))
    y.full <- sin(2 * pi * X[, 1]) + 0.4 * X[, 2]^2
    labeled <- rep(FALSE, nrow(X))
    labeled[c(1, 3, 6, 9, 12, 15, 18, 20)] <- TRUE
    y.semisupervised <- y.full
    y.semisupervised[!labeled] <- NA_real_
    lambda1 <- 0.2
    lambda2 <- 0.03

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y.semisupervised,
        k = 12L,
        tangent.dim = 2L,
        lambda1 = lambda1,
        lambda2 = lambda2,
        stabilizer = TRUE
    )

    lflag <- as.numeric(labeled)
    M <- diag(lflag) + lambda1 * as.matrix(fit$operator$B) +
        lambda2 * as.matrix(fit$operator$BS)
    ref <- as.vector(solve(M, lflag * y.full))

    expect_equal(fit$weights, lflag)
    expect_equal(fit$observed, labeled)
    expect_equal(fit$fitted.values, ref, tolerance = 1e-9)
    expect_equal(fit$residuals[labeled], y.full[labeled] - ref[labeled],
                 tolerance = 1e-9)
    expect_true(all(is.na(fit$residuals[!labeled])))
})

test_that("fit.ssrhe.hessian.regression.cv reuses the operator and selects a grid point", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    labeled <- rep(FALSE, nrow(X))
    labeled[seq(1, nrow(X), by = 2)] <- TRUE
    y[!labeled] <- NA_real_
    fold.id <- integer(nrow(X))
    fold.id[labeled] <- rep(1:3, length.out = sum(labeled))

    fit <- fit.ssrhe.hessian.regression.cv(
        X = X,
        y = y,
        k = 12L,
        tangent.dim = 2L,
        lambda1.grid = c(0.05, 0.2),
        lambda2.grid = c(0, 0.03),
        fold.id = fold.id,
        stabilizer = TRUE,
        loss = "mse",
        ridge = 1e-8
    )

    expect_s3_class(fit, "ssrhe.hessian.cv.fit")
    expect_s3_class(fit, "ssrhe.hessian.fit")
    expect_equal(nrow(fit$cv.table), 4L)
    expect_true(fit$selection$lambda1 %in% c(0.05, 0.2))
    expect_true(fit$selection$lambda2 %in% c(0, 0.03))
    expect_true(all(is.finite(fit$cv.table$cv.mean)))
    expect_equal(fit$lambda$lambda1, fit$selection$lambda1)
    expect_equal(fit$lambda$lambda2, fit$selection$lambda2)
    expect_equal(fit$fold.id, fold.id)
})

test_that("fit.ssrhe.hessian.regression.cv supports adaptive-radius neighborhoods", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2

    fit <- fit.ssrhe.hessian.regression.cv(
        X = X,
        y = y,
        tangent.dim = 2L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 3L,
        radius.rule = "geomean",
        radius.factor = 1.25,
        min.support = 8L,
        lambda1.grid = c(0.05, 0.2),
        lambda2.grid = c(0, 0.03),
        nfolds = 3L,
        stabilizer = TRUE,
        loss = "mse",
        ridge = 1e-8
    )

    expect_s3_class(fit, "ssrhe.hessian.cv.fit")
    expect_equal(fit$operator$parameters$neighborhood.type, "adaptive.radius")
    expect_true(all(fit$operator$neighborhoods$support.size >= 8L))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(fit$selection$lambda1 %in% c(0.05, 0.2))
    expect_true(fit$selection$lambda2 %in% c(0, 0.03))
})

test_that("fit.ssrhe.hessian.regression.cv rejects matrix responses for now", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- cbind(a = X[, 1], b = X[, 2])

    expect_error(
        fit.ssrhe.hessian.regression.cv(
            X = X,
            y = y,
            k = 10L,
            tangent.dim = 2L,
            lambda1.grid = c(0.1, 1)
        ),
        "one response vector"
    )
})

test_that("fit.ssrhe.hessian.regression requires BS when lambda2 is positive", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- X[, 1]

    expect_error(
        fit.ssrhe.hessian.regression(
            X = X,
            y = y,
            k = 10L,
            tangent.dim = 2L,
            lambda1 = 0.1,
            lambda2 = 0.1,
            stabilizer = FALSE
        ),
        "lambda2 > 0 requires stabilizer = TRUE"
    )
})

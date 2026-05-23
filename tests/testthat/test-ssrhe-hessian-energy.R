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
    expect_s3_class(op$local.diagnostics, "data.frame")
    expect_equal(nrow(op$local.diagnostics), nrow(X))
    expect_true(all(c("support.size", "design.rank.deficiency",
                      "pca.variance.explained", "chart.distortion",
                      "boundary.asymmetry", "curvature.bias.proxy") %in%
                        names(op$local.diagnostics)))
    expect_true(all(is.finite(op$local.diagnostics$chart.distortion)))
    expect_true(all(is.finite(op$local.diagnostics$boundary.asymmetry)))
    expect_true(all(op$local.diagnostics$design.rank.deficiency == 0L))
    expect_true(all(op$row.table$diagonal %in% c(0L, 1L)))
    expect_equal(unique(op$row.table$scale[op$row.table$diagonal == 1L]),
                 sqrt(2), tolerance = 1e-12)
})

test_that("SSRHE fit paths can skip local geometry diagnostics", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]

    op <- ssrhe.hessian.operator(
        X = X,
        k = 10L,
        tangent.dim = 2L,
        return.local.diagnostics = FALSE,
        return.timing = TRUE
    )
    expect_null(op$local.diagnostics)
    expect_s3_class(op$diagnostics, "data.frame")
    expect_true(inherits(op$B, "sparseMatrix"))
    expect_s3_class(op$timing, "data.frame")
    expect_true(all(c("phase", "elapsed.sec", "total.elapsed.sec",
                      "fraction.of.total") %in% names(op$timing)))
    expect_true(all(c("validation", "neighborhood", "native.operator",
                      "local.diagnostics", "sparse.assembly",
                      "output.finalization") %in% op$timing$phase))
    expect_true(all(op$timing$elapsed.sec >= 0))
    expect_true(all(op$timing$total.elapsed.sec > 0))
    expect_s3_class(op$native.timing, "data.frame")
    expect_true(all(c("phase", "elapsed.sec") %in% names(op$native.timing)))
    expect_true(all(c("native.input.conversion",
                      "native.support.preparation",
                      "native.local.matrix.setup",
                      "native.local.pca",
                      "native.design.construction",
                      "native.pseudoinverse",
                      "native.triplet.insertion",
                      "native.output.materialization") %in%
                        op$native.timing$phase))
    expect_true(all(is.finite(op$native.timing$elapsed.sec)))
    expect_true(all(op$native.timing$elapsed.sec >= 0))

    op.adaptive <- ssrhe.hessian.operator(
        X = X,
        tangent.dim = 2L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        min.support = 10L,
        derivative.order = 2L,
        return.local.diagnostics = FALSE,
        return.timing = TRUE
    )
    expect_s3_class(op.adaptive$native.timing, "data.frame")
    expect_s3_class(op.adaptive$neighborhoods$timing, "data.frame")
    expect_true(all(c("neighborhood.validation", "neighborhood.create.graph",
                      "neighborhood.initial.supports", "neighborhood.topup",
                      "neighborhood.truncate.reorder",
                      "neighborhood.final.validation") %in%
                        op.adaptive$neighborhoods$timing$phase))
    graph <- create.adaptive.radius.graph(
        X = X,
        k.scale = 4L,
        radius.factor = 1.25,
        radius.rule = "geomean",
        prune.method = "none",
        connect.components = FALSE
    )
    nearest <- .exact.knn.index(X, nrow(X) - 1L)
    legacy.support <- vector("list", nrow(X))
    for (i in seq_len(nrow(X))) {
        ids <- unique(c(i, graph$adj_list[[i]]))
        if (length(ids) < 10L) {
            add <- nearest[i, !nearest[i, ] %in% ids]
            need <- 10L - length(ids)
            ids <- unique(c(ids, add[seq_len(min(need, length(add)))]))
        }
        if (ids[[1L]] != i) {
            ids <- c(i, setdiff(ids, i))
        }
        legacy.support[[i]] <- as.integer(ids)
    }
    expect_equal(op.adaptive$support.index, legacy.support)

    fit.fast <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda1 = 0.2,
        lambda2 = 0,
        return.timing = TRUE
    )
    expect_null(fit.fast$operator$local.diagnostics)
    expect_s3_class(fit.fast$operator$timing, "data.frame")

    fit.audit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda1 = 0.2,
        lambda2 = 0,
        return.local.diagnostics = TRUE
    )
    expect_s3_class(fit.audit$operator$local.diagnostics, "data.frame")
    expect_equal(nrow(fit.audit$operator$local.diagnostics), nrow(X))
})

test_that("ssrhe.hessian.operator local solver backends match SVD on full-rank supports", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
                               y = seq(0, 1, length.out = 6)))
    op.svd <- ssrhe.hessian.operator(
        X = X,
        k = 16L,
        tangent.dim = 2L,
        local.solver = "svd",
        return.local.diagnostics = FALSE,
        return.timing = TRUE
    )
    op.qr <- ssrhe.hessian.operator(
        X = X,
        k = 16L,
        tangent.dim = 2L,
        local.solver = "qr",
        return.local.diagnostics = FALSE,
        return.timing = TRUE
    )
    op.ne <- ssrhe.hessian.operator(
        X = X,
        k = 16L,
        tangent.dim = 2L,
        local.solver = "normal.equations",
        return.local.diagnostics = FALSE,
        return.timing = TRUE
    )

    expect_equal(as.matrix(op.qr$A), as.matrix(op.svd$A), tolerance = 1e-9)
    expect_equal(as.matrix(op.ne$A), as.matrix(op.svd$A), tolerance = 1e-8)
    expect_true(all(op.qr$diagnostics$local.solver == "qr"))
    expect_true(all(op.ne$diagnostics$local.solver == "normal.equations"))
    expect_true(all(op.qr$diagnostics$solver.fallback == 0L))
    expect_true(all(op.ne$diagnostics$solver.fallback == 0L))
    expect_equal(op.qr$parameters$local.solver, "qr")
    expect_equal(op.ne$parameters$local.solver, "normal.equations")
})

test_that("ssrhe.hessian.operator local diagnostics flag chart geometry", {
    skip_if_not_installed("Matrix")

    set.seed(24)
    x <- seq(-1, 1, length.out = 28)
    X.line <- cbind(x, 0)
    op.line <- ssrhe.hessian.operator(
        X = X.line,
        k = 8L,
        tangent.dim = 1L,
        return.BS = FALSE
    )

    X.curve <- cbind(x, x^2)
    op.curve <- ssrhe.hessian.operator(
        X = X.curve,
        k = 8L,
        tangent.dim = 1L,
        return.BS = FALSE
    )

    expect_lt(max(op.line$local.diagnostics$chart.distortion), 1e-10)
    expect_lt(max(op.line$local.diagnostics$curvature.bias.proxy), 1e-10)
    expect_gt(stats::median(op.curve$local.diagnostics$curvature.bias.proxy), 0)
    expect_gt(stats::median(op.curve$local.diagnostics$chart.distortion), 0)
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

test_that("ssrhe.hessian.operator constructs order-3 symmetric tensor rows", {
    skip_if_not_installed("Matrix")

    grid <- expand.grid(x = seq(0, 1, length.out = 5),
                        y = seq(0, 1, length.out = 5))
    X <- as.matrix(grid)
    op <- ssrhe.hessian.operator(
        X = X,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        return.BS = FALSE
    )

    expect_s3_class(op, "ssrhe.hessian.operator")
    expect_equal(op$parameters$derivative.order, 3L)
    expect_equal(op$parity$derivative.order, 3L)
    expect_equal(dim(op$A), c(nrow(X) * 4L, nrow(X)))
    expect_equal(dim(op$B), c(nrow(X), nrow(X)))
    expect_null(op$BS)
    expect_true(all(op$row.table$derivative.order == 3L))
    expect_equal(sort(unique(op$row.table$component)), 1:4)
    expect_equal(sort(unique(op$row.table$c)), 1:2)

    one_center <- op$row.table[op$row.table$center == 1L, ]
    expect_equal(
        one_center[, c("a", "b", "c")],
        data.frame(a = c(1L, 1L, 1L, 2L),
                   b = c(1L, 1L, 2L, 2L),
                   c = c(1L, 2L, 2L, 2L))
    )
    expect_equal(one_center$derivative.scale, c(6, 2, 2, 6),
                 tolerance = 1e-12)
    expect_equal(one_center$tensor.multiplicity, c(1L, 3L, 3L, 1L))
    expect_equal(one_center$scale, c(6, 2 * sqrt(3), 2 * sqrt(3), 6),
                 tolerance = 1e-12)
})

test_that("order-3 SSRHE operator annihilates one-dimensional quadratics", {
    skip_if_not_installed("Matrix")

    x <- seq(-1, 1, length.out = 30)
    X <- matrix(x, ncol = 1)
    op <- ssrhe.hessian.operator(
        X = X,
        k = 10L,
        tangent.dim = 1L,
        derivative.order = 3L,
        return.BS = FALSE
    )

    one <- rep(1, length(x))
    linear <- 2 - 0.5 * x
    quadratic <- 0.25 + 2 * x - 3 * x^2
    cubic <- x^3

    expect_equal(dim(op$A), c(length(x), length(x)))
    expect_equal(unique(op$row.table$scale), 6, tolerance = 1e-12)
    expect_lt(sqrt(sum((op$A %*% one)^2)), 1e-8)
    expect_lt(sqrt(sum((op$A %*% linear)^2)), 1e-8)
    expect_lt(sqrt(sum((op$A %*% quadratic)^2)), 1e-7)
    expect_gt(sqrt(sum((op$A %*% cubic)^2)), 1e-3)
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

test_that("fit.ssrhe.hessian.regression matches order-3 dense fixed-lambda reference", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + X[, 1] * X[, 2] + 0.1 * X[, 2]^3
    lambda1 <- 0.2

    fit <- fit.ssrhe.hessian.regression(
        X = X,
        y = y,
        k = 12L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda1 = lambda1,
        lambda2 = 0
    )

    M <- diag(1, nrow(X)) + lambda1 * as.matrix(fit$operator$B)
    ref <- as.vector(solve(M, y))

    expect_s3_class(fit, "ssrhe.hessian.fit")
    expect_equal(fit$operator$parameters$derivative.order, 3L)
    expect_equal(fit$fitted.values, ref, tolerance = 1e-9)
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

test_that("ssrhe.support.grid builds compact adaptive-radius profiles", {
    grid <- ssrhe.support.grid(
        n = 30L,
        tangent.dim = 1L,
        derivative.order = 3L,
        max.candidates = 4L
    )

    expect_s3_class(grid, "data.frame")
    expect_named(grid, c("adaptive.k.scale", "min.support", "max.support"))
    expect_lte(nrow(grid), 4L)
    expect_true(all(grid$adaptive.k.scale >= 1L))
    expect_true(all(grid$adaptive.k.scale < 30L))
    expect_true(all(grid$min.support >= 2L))
    expect_true(all(grid$min.support <= 30L))
})

test_that("fit.ssrhe.hessian.regression.cv supports outer support CV", {
    skip_if_not_installed("Matrix")

    x <- seq(0, 1, length.out = 24)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)
    fold.id <- rep(1:3, length.out = length(y))
    support.grid <- data.frame(
        adaptive.k.scale = c(2L, 4L),
        min.support = c(5L, 8L)
    )

    fit <- fit.ssrhe.hessian.regression.cv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        support.selection = "cv",
        support.grid = support.grid,
        lambda1.grid = c(1e-4, 1e-2),
        lambda2.grid = 0,
        fold.id = fold.id,
        ridge = 1e-8
    )

    expect_s3_class(fit, "ssrhe.hessian.cv.fit")
    expect_equal(fit$support.selection, "cv")
    expect_equal(nrow(fit$support.cv.table), 2L)
    expect_true(all(fit$support.cv.table$status == "ok"))
    expect_true(fit$selected.support$min.support %in% support.grid$min.support)
    expect_equal(fit$operator$neighborhoods$min.support,
                 fit$selected.support$min.support)
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.regression.gcv selects a finite grid point", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2

    fit <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        k = 12L,
        tangent.dim = 2L,
        lambda1.grid = c(0.05, 0.2),
        lambda2.grid = c(0, 0.03),
        stabilizer = TRUE,
        ridge = 1e-8
    )

    expect_s3_class(fit, "ssrhe.hessian.gcv.fit")
    expect_s3_class(fit, "ssrhe.hessian.fit")
    expect_equal(nrow(fit$gcv.table), 4L)
    expect_true(fit$selection$lambda1 %in% c(0.05, 0.2))
    expect_true(fit$selection$lambda2 %in% c(0, 0.03))
    expect_true(all(is.finite(fit$gcv.table$gcv)))
    expect_true(all(is.finite(fit$gcv.table$trace.S)))
    expect_equal(fit$lambda$lambda1, fit$selection$lambda1)
    expect_equal(fit$lambda$lambda2, fit$selection$lambda2)
})

test_that("fit.ssrhe.hessian.regression.gcv traces shrink with stronger penalty", {
    skip_if_not_installed("Matrix")

    x <- seq(0, 1, length.out = 24)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)

    fit <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        min.support = 8L,
        lambda1.grid = c(1e-4, 1e-1),
        lambda2.grid = 0,
        ridge = 1e-8
    )

    ordered <- fit$gcv.table[order(fit$gcv.table$lambda1), ]
    expect_lte(ordered$trace.S[2L], ordered$trace.S[1L] + 1e-7)
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.regression.gcv supports Hutchinson trace estimates", {
    skip_if_not_installed("Matrix")

    x <- seq(0, 1, length.out = 18)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)

    exact <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        min.support = 8L,
        lambda1.grid = c(1e-4, 1e-2),
        lambda2.grid = 0,
        ridge = 1e-8,
        gcv.trace.method = "exact"
    )
    hutch <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        min.support = 8L,
        lambda1.grid = c(1e-4, 1e-2),
        lambda2.grid = 0,
        ridge = 1e-8,
        gcv.trace.method = "hutchinson",
        gcv.trace.n.probes = 600L,
        gcv.trace.seed = 19L
    )
    hutch.again <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        min.support = 8L,
        lambda1.grid = c(1e-4, 1e-2),
        lambda2.grid = 0,
        ridge = 1e-8,
        gcv.trace.method = "hutchinson",
        gcv.trace.n.probes = 600L,
        gcv.trace.seed = 19L
    )

    expect_equal(hutch$gcv.table$trace.S,
                 hutch.again$gcv.table$trace.S,
                 tolerance = 0)
    expect_equal(unique(hutch$gcv.table$trace.method), "hutchinson")
    expect_equal(unique(hutch$gcv.table$trace.n.probes), 600L)
    expect_true(all(is.finite(hutch$gcv.table$trace.se)))
    expect_equal(hutch$gcv.trace$seed, 19L)
    expect_equal(hutch$gcv.trace$method, "hutchinson")
    expect_lt(max(abs(hutch$gcv.table$trace.S - exact$gcv.table$trace.S)),
              1.25)
})

test_that("fit.ssrhe.hessian.regression.gcv supports outer support GCV", {
    skip_if_not_installed("Matrix")

    x <- seq(0, 1, length.out = 24)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)
    support.grid <- data.frame(
        adaptive.k.scale = c(2L, 4L),
        min.support = c(5L, 8L)
    )

    fit <- fit.ssrhe.hessian.regression.gcv(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        support.selection = "gcv",
        support.grid = support.grid,
        lambda1.grid = c(1e-4, 1e-2),
        lambda2.grid = 0,
        ridge = 1e-8
    )

    expect_s3_class(fit, "ssrhe.hessian.gcv.fit")
    expect_equal(fit$support.selection, "gcv")
    expect_equal(nrow(fit$support.gcv.table), 2L)
    expect_true(all(fit$support.gcv.table$status == "ok"))
    expect_true(fit$selected.support$min.support %in% support.grid$min.support)
    expect_equal(fit$operator$neighborhoods$min.support,
                 fit$selected.support$min.support)
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.regression.gcv rejects missing responses", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- X[, 1]
    y[2L] <- NA_real_

    expect_error(
        fit.ssrhe.hessian.regression.gcv(
            X = X,
            y = y,
            k = 10L,
            tangent.dim = 2L,
            lambda1.grid = c(0.1, 1)
        ),
        "fully observed"
    )
})

test_that("fit.ssrhe.hessian.l1.regression matches direct genlasso reference", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    lambda <- 0.08

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda.grid = lambda,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )

    ref.path <- genlasso::genlasso(
        y = y,
        D = as.matrix(fit$operator$A),
        svd = TRUE,
        maxsteps = 1000L
    )
    ref <- as.vector(stats::coef(ref.path, lambda = lambda)$beta)

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$fitted.values, ref, tolerance = 1e-8)
    expect_equal(fit$residuals, y - ref, tolerance = 1e-8)
    expect_equal(fit$energies$hessian.l1,
                 sum(abs(as.vector(fit$operator$A %*% ref))),
                 tolerance = 1e-8)
    expect_equal(fit$solver$backend, "genlasso")
})

test_that("fit.ssrhe.hessian.l1.regression matches direct genlasso reference for order 3", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2 +
        0.1 * X[, 1] * X[, 2]
    lambda <- 0.015

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda.grid = lambda,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )

    ref.path <- genlasso::genlasso(
        y = y,
        D = as.matrix(fit$operator$A),
        svd = TRUE,
        maxsteps = 1000L
    )
    ref <- as.vector(stats::coef(ref.path, lambda = lambda)$beta)

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$operator$parameters$derivative.order, 3L)
    expect_true(all(fit$operator$row.table$derivative.order == 3L))
    expect_equal(fit$fitted.values, ref, tolerance = 1e-8)
    expect_equal(fit$energies$hessian.l1,
                 sum(abs(as.vector(fit$operator$A %*% ref))),
                 tolerance = 1e-8)
    expect_equal(fit$solver$backend, "genlasso")
})

test_that("fit.ssrhe.hessian.l1.regression supports ADMM and row scaling diagnostics", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda.grid = 0.01,
        lambda.selection = "fixed",
        solver = "admm",
        row.scaling = "l2",
        admm.maxiter = 2000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$solver$backend, "admm")
    expect_equal(fit$solver$row.scaling, "l2")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(fit$diagnostics$scaled.row.norm$min, 1, tolerance = 1e-8)
    expect_equal(fit$diagnostics$scaled.row.norm$median, 1, tolerance = 1e-8)
    expect_equal(fit$diagnostics$scaled.row.norm$max, 1, tolerance = 1e-8)
    expect_true(is.list(fit$solver$admm[[1L]]))
})

test_that("ADMM fixed-lambda fit is close to genlasso on a small order-3 operator", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    lambda <- 0.01

    fit.gen <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda.grid = lambda,
        lambda.selection = "fixed",
        solver = "genlasso",
        maxsteps = 1000L
    )
    fit.admm <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda.grid = lambda,
        lambda.selection = "fixed",
        solver = "admm",
        row.scaling = "none",
        admm.maxiter = 5000L,
        admm.abstol = 1e-5,
        admm.reltol = 1e-4
    )

    expect_equal(fit.admm$fitted.values, fit.gen$fitted.values,
                 tolerance = 2e-3)
})

test_that("fit.ssrhe.hessian.l1.regression selects a CV lambda", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    fold.id <- rep(1:4, length.out = nrow(X))
    lambda.grid <- c(0.02, 0.08, 0.2)

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda.grid = lambda.grid,
        lambda.selection = "cv",
        fold.id = fold.id,
        loss = "mse",
        maxsteps = 1000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_true(fit$lambda %in% lambda.grid)
    expect_equal(fit$cv$lambda, fit$lambda)
    expect_equal(fit$fold.id, fold.id)
    expect_true(all(is.finite(fit$cv$mean.error)))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(dim(fit$beta.grid), c(nrow(X), length(lambda.grid)))
})

test_that("fit.ssrhe.hessian.l1.regression handles missing responses and weights", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- X[, 1]^2 - X[, 2]
    y[c(3, 9)] <- NA_real_
    weights <- rep(1, nrow(X))
    weights[c(5, 11)] <- 0.25

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 10L,
        tangent.dim = 2L,
        lambda.grid = 0.05,
        lambda.selection = "fixed",
        weights = weights,
        maxsteps = 1000L
    )

    expect_equal(fit$weights[c(3, 9)], c(0, 0))
    expect_equal(fit$weights[c(5, 11)], c(0.25, 0.25))
    expect_true(all(is.na(fit$residuals[c(3, 9)])))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.l1.regression selects a CV lambda for order 3", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2
    fold.id <- rep(1:5, length.out = nrow(X))
    lambda.grid <- c(0.001, 0.01, 0.05)

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 14L,
        tangent.dim = 2L,
        derivative.order = 3L,
        lambda.grid = lambda.grid,
        lambda.selection = "cv",
        fold.id = fold.id,
        loss = "mse",
        maxsteps = 1000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$operator$parameters$derivative.order, 3L)
    expect_true(fit$lambda %in% lambda.grid)
    expect_equal(fit$cv$lambda, fit$lambda)
    expect_true(all(is.finite(fit$cv$mean.error)))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("refit.ssrhe.hessian.l1.regression reuses the operator", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y1 <- X[, 1] + X[, 2]
    y2 <- sin(X[, 1]) - cos(X[, 2])

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y1,
        k = 10L,
        tangent.dim = 2L,
        lambda.grid = 0.04,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )
    refit <- refit.ssrhe.hessian.l1.regression(
        fit,
        y.new = y2,
        lambda.grid = 0.12,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )

    expect_s3_class(refit, "ssrhe.hessian.l1.refit")
    expect_equal(refit$lambda, 0.12)
    expect_identical(refit$operator, fit$operator)
    expect_true(all(is.finite(refit$fitted.values)))
})

test_that("fit.ssrhe.hessian.l1.regression supports adaptive-radius neighborhoods", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                               y = seq(0, 1, length.out = 4)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        tangent.dim = 2L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 3L,
        radius.rule = "geomean",
        radius.factor = 1.25,
        min.support = 8L,
        lambda.grid = 0.05,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$operator$parameters$neighborhood.type, "adaptive.radius")
    expect_true(all(fit$operator$neighborhoods$support.size >= 8L))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.l1.regression supports adaptive-radius order-3 neighborhoods", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
                               y = seq(0, 1, length.out = 5)))
    y <- sin(2 * pi * X[, 1]) + 0.25 * X[, 2]^2

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        tangent.dim = 2L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        adaptive.k.scale = 4L,
        radius.rule = "geomean",
        radius.factor = 1.25,
        min.support = 14L,
        lambda.grid = 0.01,
        lambda.selection = "fixed",
        maxsteps = 1000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$operator$parameters$neighborhood.type, "adaptive.radius")
    expect_equal(fit$operator$parameters$derivative.order, 3L)
    expect_true(all(fit$operator$neighborhoods$support.size >= 14L))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.l1.regression handles one-dimensional SSRHE operators", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("genlasso")

    x <- seq(0, 1, length.out = 24)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        k = 8L,
        tangent.dim = 1L,
        lambda.grid = c(0.001, 0.01),
        lambda.selection = "cv",
        nfolds = 3L,
        maxsteps = 1000L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_true(fit$solver$svd)
    expect_equal(fit$solver$representation, "dense")
    expect_true(fit$lambda %in% c(0.001, 0.01))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.ssrhe.hessian.l1.regression supports outer support CV", {
    skip_if_not_installed("Matrix")

    x <- seq(0, 1, length.out = 24)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)
    fold.id <- rep(1:3, length.out = length(y))
    support.grid <- data.frame(
        adaptive.k.scale = c(2L, 4L),
        min.support = c(5L, 8L)
    )

    fit <- fit.ssrhe.hessian.l1.regression(
        X = X,
        y = y,
        tangent.dim = 1L,
        derivative.order = 3L,
        neighborhood.type = "adaptive.radius",
        support.selection = "cv",
        support.grid = support.grid,
        lambda.grid = c(0.001, 0.01),
        lambda.selection = "cv",
        fold.id = fold.id,
        solver = "admm",
        row.scaling = "l2",
        admm.maxiter = 500L
    )

    expect_s3_class(fit, "ssrhe.hessian.l1.fit")
    expect_equal(fit$support.selection, "cv")
    expect_equal(nrow(fit$support.cv.table), 2L)
    expect_true(all(fit$support.cv.table$status == "ok"))
    expect_true(fit$selected.support$min.support %in% support.grid$min.support)
    expect_equal(fit$operator$neighborhoods$min.support,
                 fit$selected.support$min.support)
    expect_true(all(is.finite(fit$fitted.values)))
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

test_that("order-3 SSRHE L2 rejects supplemental stabilizer paths", {
    skip_if_not_installed("Matrix")

    X <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4),
                               y = seq(0, 1, length.out = 4)))
    y <- X[, 1]

    expect_error(
        ssrhe.hessian.operator(
            X = X,
            k = 12L,
            tangent.dim = 2L,
            derivative.order = 3L,
            stabilizer = TRUE
        ),
        "stabilizer"
    )
    expect_error(
        fit.ssrhe.hessian.regression(
            X = X,
            y = y,
            k = 12L,
            tangent.dim = 2L,
            derivative.order = 3L,
            lambda1 = 0.1,
            lambda2 = 0.1
        ),
        "lambda2/stabilizer"
    )
    expect_error(
        fit.ssrhe.hessian.regression.cv(
            X = X,
            y = y,
            k = 12L,
            tangent.dim = 2L,
            derivative.order = 3L,
            lambda1.grid = c(0.1, 0.2),
            lambda2.grid = c(0, 0.1)
        ),
        "lambda2.grid/stabilizer"
    )
})

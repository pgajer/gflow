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

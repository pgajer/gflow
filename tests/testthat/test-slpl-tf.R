test_that("slpl.tf.operator builds inclusive synchronization rows", {
    x <- seq(-1, 1, length.out = 7)
    X <- matrix(x, ncol = 1)
    op <- slpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    expect_s3_class(op, "slpl_tf_operator")
    expect_s3_class(op$lpl.operator, "lpl_tf_operator")
    expect_equal(dim(op$A_LPL), c(7L, 7L))
    expect_gt(nrow(op$C_sync), 0L)
    expect_true(all(op$sync.map.table$sync_map_self_inclusion == "inclusive"))
    expect_true(all(op$sync.map.table$lpl_residual_self_inclusion ==
                        "self_excluded"))
    expect_equal(op$sync.row.table$row.norm.final[
        op$sync.row.table$status == "ok"
    ], rep(1, nrow(op$C_sync)), tolerance = 1e-12)
    probes <- cbind(1, x, x^2)
    expect_lt(max(abs(as.matrix(op$C_sync %*% probes))), 1e-10)
    expect_lt(op$diagnostics$sync.polynomial.reproduction.error, 1e-10)
    expect_equal(op$diagnostics$sync.coverage.zero, 0L)
    expect_true(op$settings$support.metric %in% c("coordinates", "graph.geodesic"))
})

test_that("fit.slpl.tf with lambda2 zero nests fit.lpl.tf", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)
    op <- slpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    lambda1 <- 0.02
    lpl.fit <- fit.lpl.tf(
        y = y,
        operator = op$lpl.operator,
        lambda = lambda1,
        maxsteps = 500L
    )
    slpl.fit <- fit.slpl.tf(
        y = y,
        operator = op,
        lambda1 = lambda1,
        lambda2 = 0,
        maxsteps = 500L
    )
    expect_s3_class(slpl.fit, "slpl_tf")
    expect_equal(slpl.fit$fitted.values, lpl.fit$fitted.values,
                 tolerance = 1e-8)
    expect_equal(as.numeric(predict(slpl.fit)), slpl.fit$fitted.values)
    expect_equal(slpl.fit$lambda1, lambda1)
    expect_equal(slpl.fit$lambda2, 0)
})

test_that("fit.slpl.tf positive lambda2 reduces synchronization energy", {
    skip_if_not_installed("genlasso")
    x <- seq(-1, 1, length.out = 9)
    y <- abs(x) + 0.4 * exp(-80 * (x - 0.35)^2)
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    fit0 <- fit.slpl.tf(y = y, operator = op, lambda1 = 0.05,
                        lambda2 = 0, maxsteps = 500L)
    fit10 <- fit.slpl.tf(y = y, operator = op, lambda1 = 0.05,
                         lambda2 = 10, maxsteps = 500L)
    expect_lt(fit10$diagnostics$sync.energy.l2,
              fit0$diagnostics$sync.energy.l2)
    expect_true(all(is.finite(fit10$fitted.values)))
})

test_that("fit.slpl.tf positive lambda2 uses compressed quadratic design", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 10)
    y <- sin(2 * pi * x)
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 5L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    fit <- fit.slpl.tf(y = y, operator = op, lambda1 = 0.02,
                       lambda2 = 0.5, maxsteps = 500L)
    expect_equal(fit$solver$augmented.design.dim,
                 paste(length(y), length(y), sep = " x "))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("refit.slpl.tf reuses fixed operator and penalties", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 10)
    y <- x^2
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2"
    )
    fit <- fit.slpl.tf(y = y, operator = op, lambda1 = 0.02, lambda2 = 1,
                       maxsteps = 500L)
    refit <- refit.slpl.tf(fit, y = y + 0.1 * sin(2 * pi * x))
    expect_s3_class(refit, "slpl_tf")
    expect_identical(refit$operator, fit$operator)
    expect_equal(refit$lambda1, fit$lambda1)
    expect_equal(refit$lambda2, fit$lambda2)
    expect_error(predict(refit, newdata = matrix(0, ncol = 1)), "newdata")
})

test_that("fit.slpl.tf deterministic CV materializes fold and candidate tables", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    y <- sin(2 * pi * x)
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    fit <- fit.slpl.tf(
        y = y,
        operator = op,
        lambda.selection = "cv",
        lambda1.grid = c(0.01, 0.05),
        lambda2.grid = c(0, 1),
        cv.folds = 3L,
        maxsteps = 500L
    )
    expect_s3_class(fit, "slpl_tf")
    expect_equal(fit$foldid, rep(1:3, length.out = length(x)))
    expect_equal(nrow(fit$cv$candidate.table), 4L)
    expect_equal(sum(fit$cv$candidate.table$selected), 1L)
    expect_true(fit$lambda1 %in% c(0.01, 0.05))
    expect_true(fit$lambda2 %in% c(0, 1))
    expect_true(all(fit$cv$candidate.table$status == "ok"))
    expect_equal(dim(fit$cv$fold.errors), c(3L, 2L, 2L))
    expect_equal(fit$fold.source, "generated_deterministic")
})

test_that("fit.slpl.tf RMSE CV selects on the reported RMSE scale", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    y <- sin(2 * pi * x) + 0.25 * x
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    fit <- fit.slpl.tf(
        y = y,
        operator = op,
        lambda.selection = "cv",
        lambda1.grid = c(0.001, 0.01, 0.05),
        lambda2.grid = c(0, 0.5),
        cv.folds = 3L,
        cv.loss = "rmse",
        selection = "min",
        maxsteps = 500L
    )
    candidates <- fit$cv$candidate.table
    selected <- which(candidates$selected)
    expect_equal(fit$cv$loss, "rmse")
    expect_equal(selected, which.min(candidates$mean.error))
    expect_equal(fit$cv$selected.idx, selected)
    expect_equal(fit$cv$error.selected, candidates$mean.error[[selected]])
    expect_equal(fit$lambda1, candidates$lambda1[[selected]])
    expect_equal(fit$lambda2, candidates$lambda2[[selected]])
})

test_that("fit.slpl.tf one-SE RMSE selection and sync.lambda.grid alias are deterministic", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    y <- cos(2 * pi * x) + 0.2 * x^2
    op <- slpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    fit <- fit.slpl.tf(
        y = y,
        operator = op,
        lambda.selection = "cv",
        lambda1.grid = c(0.001, 0.01, 0.05),
        sync.lambda.grid = c(0, 0.5),
        cv.folds = 3L,
        cv.loss = "rmse",
        selection = "one.se",
        maxsteps = 500L
    )
    candidates <- fit$cv$candidate.table
    best <- which.min(candidates$mean.error)
    threshold <- candidates$mean.error[[best]] + candidates$se[[best]]
    eligible <- which(candidates$mean.error <= threshold)
    expected <- eligible[order(-candidates$lambda1[eligible],
                               -candidates$lambda2[eligible],
                               candidates$mean.error[eligible])[[1L]]]
    expect_equal(fit$lambda2.grid, c(0, 0.5))
    expect_equal(fit$cv$loss, "rmse")
    expect_equal(which(candidates$selected), expected)
    expect_equal(fit$cv$selected.idx, expected)
    expect_equal(fit$lambda1, candidates$lambda1[[expected]])
    expect_equal(fit$lambda2, candidates$lambda2[[expected]])
})

test_that("fit.slpl.tf operator grid records failed candidates", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    y <- x^2
    fit <- fit.slpl.tf(
        X = matrix(x, ncol = 1),
        y = y,
        lambda.selection = "cv",
        lambda1.grid = c(0.01, 0.05),
        lambda2.grid = c(0, 0.5),
        cv.folds = 3L,
        operator.grid = list(
            list(degree = 2L, support.type = "knn", support.size = 2L,
                 support.buffer = 0L, kernel = "gaussian"),
            list(degree = 2L, support.type = "knn", support.size = 4L,
                 support.buffer = 0L, kernel = "gaussian")
        ),
        row.normalize = "l2",
        local.solver = "svd",
        maxsteps = 500L
    )
    expect_s3_class(fit, "slpl_tf")
    expect_equal(nrow(fit$operator.selection$table), 2L)
    expect_true(any(fit$operator.selection$table$status == "failed"))
    expect_true(any(fit$operator.selection$table$status == "ok"))
    expect_equal(sum(fit$operator.selection$table$selected), 1L)
    expect_identical(fit$operator.selection$foldid,
                     rep(1:3, length.out = length(x)))
    expect_equal(fit$fold.source, "generated_deterministic")
    expect_equal(fit$operator.selection$fold.source, "generated_deterministic")
})

test_that("slpl.tf.operator supports graph-geodesic synchronization", {
    x <- seq(0, 1, length.out = 9)
    X <- matrix(x, ncol = 1)
    adj <- vector("list", length(x))
    wt <- vector("list", length(x))
    for (ii in seq_along(x)) {
        nb <- integer(0)
        ww <- numeric(0)
        if (ii > 1L) {
            nb <- c(nb, ii - 1L)
            ww <- c(ww, abs(x[ii] - x[ii - 1L]))
        }
        if (ii < length(x)) {
            nb <- c(nb, ii + 1L)
            ww <- c(ww, abs(x[ii + 1L] - x[ii]))
        }
        adj[[ii]] <- nb
        wt[[ii]] <- ww
    }
    op <- slpl.tf.operator(
        X,
        adj.list = adj,
        weight.list = wt,
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        support.metric = "graph.geodesic",
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "svd"
    )
    expect_s3_class(op, "slpl_tf_operator")
    expect_equal(op$settings$support.metric, "graph.geodesic")
    expect_gt(nrow(op$C_sync), 0L)
    expect_lt(op$diagnostics$sync.polynomial.reproduction.error, 1e-10)
})

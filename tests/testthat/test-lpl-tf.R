test_that("lpl.tf.operator reproduces the tiny 1D divided-difference row", {
    X <- matrix(0:3, ncol = 1)
    op <- lpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        anchor.index = 1L,
        kernel = "gaussian",
        row.normalize = "none"
    )
    expect_s3_class(op, "lpl_tf_operator")
    expect_equal(dim(op$A), c(1L, 4L))
    expect_equal(as.numeric(op$A[1, ]), c(1, -3, 3, -1), tolerance = 1e-10)
    probes <- cbind(1, X[, 1], X[, 1]^2)
    expect_lt(max(abs(as.matrix(op$A %*% probes))), 1e-10)
    expect_true(all(is.finite(op$row.table$rank.tolerance)))
    expect_match(op$row.table$rank.tolerance.rule,
                 "Machine\\$double\\.eps")
    expect_equal(op$diagnostics$operator.rank, 1L)
    expect_equal(op$diagnostics$operator.nullity, 3L)
})

test_that("lpl.tf.operator matches unequal-spacing 1D divided differences", {
    x <- c(0, 0.07, 0.19, 0.31, 0.55, 0.7, 0.83, 1)
    op <- lpl.tf.operator(
        matrix(x, ncol = 1),
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        anchor.index = 1L,
        kernel = "gaussian",
        row.normalize = "l2"
    )
    D <- matrix(0, nrow = 1L, ncol = length(x))
    for (i in seq_len(nrow(D))) {
        idx <- i:(i + 3L)
        co <- vapply(idx, function(k) {
            others <- setdiff(idx, k)
            1 / prod(x[k] - x[others])
        }, numeric(1))
        D[i, idx] <- co
    }
    D <- D / sqrt(rowSums(D^2))
    got <- as.matrix(op$A)
    row.err <- vapply(seq_len(nrow(D)), function(i) {
        1 - abs(sum(got[i, ] * D[i, ]))
    }, numeric(1))
    expect_lt(max(row.err), 1e-10)
    expect_equal(op$diagnostics$operator.nullity, length(x) - 1L)
})

test_that("lpl.tf.operator reproduces 2D quadratic probes", {
    X <- rbind(
        c(0.00, 0.00), c(0.18, 0.03), c(0.05, 0.25), c(0.31, 0.19),
        c(0.62, 0.12), c(0.18, 0.71), c(0.76, 0.64), c(0.94, 0.38)
    )
    op <- lpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 8L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2"
    )
    B <- cbind(
        1,
        X[, 1], X[, 2],
        X[, 1]^2, X[, 1] * X[, 2], X[, 2]^2
    )
    expect_true(all(op$row.table$status == "ok"))
    expect_lt(max(abs(as.matrix(op$A %*% B))), 1e-10)
    expect_lt(op$diagnostics$polynomial.residual, 1e-10)
    expect_equal(op$diagnostics$operator.nullity, 6L)
})

test_that("lpl.tf.operator qr solver reproduces 2D quadratic probes", {
    X <- rbind(
        c(0.00, 0.00), c(0.18, 0.03), c(0.05, 0.25), c(0.31, 0.19),
        c(0.62, 0.12), c(0.18, 0.71), c(0.76, 0.64), c(0.94, 0.38)
    )
    op <- lpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 8L,
        support.buffer = 0L,
        kernel = "gaussian",
        row.normalize = "l2",
        local.solver = "qr"
    )
    B <- cbind(
        1,
        X[, 1], X[, 2],
        X[, 1]^2, X[, 1] * X[, 2], X[, 2]^2
    )
    expect_true(all(op$row.table$status == "ok"))
    expect_true(all(op$row.table$solver.used == "qr"))
    expect_lt(max(abs(as.matrix(op$A %*% B))), 1e-10)
})

test_that("lpl.tf.operator drops rank-deficient requested-degree rows", {
    X <- cbind(seq(0, 1, length.out = 8), 0)
    op <- lpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 8L,
        support.buffer = 0L,
        kernel = "gaussian"
    )
    expect_equal(nrow(op$A), 0L)
    expect_true(all(op$row.table$status == "dropped"))
    expect_true(all(op$row.table$drop.reason == "rank_deficient_requested_degree"))
    expect_true(all(is.finite(op$row.table$rank.tolerance)))
    expect_match(op$row.table$rank.tolerance.rule[[1L]],
                 "Machine\\$double\\.eps")
})

test_that("lpl.tf.operator enforces phase-1 input restrictions", {
    X <- matrix(seq_len(10), ncol = 1)
    expect_error(
        lpl.tf.operator(X, degree = 2L, support.size = 4L,
                        exclude.self = FALSE),
        "exclude.self"
    )
    expect_error(
        lpl.tf.operator(X, degree = 2L, support.size = 4L,
                        anchor.coordinates = matrix(0, 1, 1)),
        "observed anchors"
    )
    expect_error(
        lpl.tf.operator(rbind(c(0, 0), c(0, 0)), duplicate.action = "error"),
        "Duplicate"
    )
})

test_that("lpl.tf.operator supports graph-geodesic support distances", {
    X <- rbind(c(0, 0), c(1, 0), c(2, 0), c(3, 0), c(0, 1), c(0, 2))
    adj <- list(
        c(2L, 5L), c(1L, 3L), c(2L, 4L),
        3L, c(1L, 6L), 5L
    )
    weights <- lapply(seq_along(adj), function(i) rep(1, length(adj[[i]])))

    op <- lpl.tf.operator(
        X,
        adj.list = adj,
        weight.list = weights,
        degree = 1L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        support.metric = "graph.geodesic",
        kernel = "gaussian"
    )

    expect_s3_class(op, "lpl_tf_operator")
    expect_equal(op$settings$support.metric, "graph.geodesic")
    expect_equal(op$settings$graph.summary$n.edges, 5L)
    expect_true(all(op$row.table$support.metric == "graph.geodesic"))
    expect_equal(op$supports[[4L]]$predictors, c(3L, 2L, 1L))
})

test_that("lpl.tf.operator accepts create.rknn.graph objects", {
    X <- cbind(seq(0, 1, length.out = 8), 0.2 * seq(0, 1, length.out = 8))
    graph <- create.rknn.graph(
        X,
        type = "adaptive.radius",
        k.scale = 2L,
        radius.factor = 1.25,
        radius.rule = "geomean",
        prune.method = "none",
        connect.components = TRUE,
        connect.method = "component.mst"
    )

    op <- lpl.tf.operator(
        X,
        graph = graph,
        degree = 1L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        support.metric = "auto"
    )

    expect_equal(op$settings$support.metric, "graph.geodesic")
    expect_equal(op$settings$graph.summary$constructor, "graph.final")
    expect_equal(op$settings$graph.summary$n.vertices, nrow(X))
})

test_that("lpl.tf.operator extracts final raw and pruned graph stages", {
    X <- cbind(seq(0, 1, length.out = 9), 0.1 * seq(0, 1, length.out = 9))
    graph <- create.rknn.graph(
        X,
        type = "adaptive.radius",
        k.scale = 2L,
        radius.factor = 1.25,
        radius.rule = "geomean",
        prune.method = "none",
        connect.components = TRUE,
        connect.method = "component.mst",
        graph.detail = "full"
    )

    for (stage in c("final", "raw", "pruned")) {
        op <- lpl.tf.operator(
            X,
            graph = graph,
            graph.stage = stage,
            degree = 1L,
            support.type = "knn",
            support.size = 4L,
            support.buffer = 0L,
            support.metric = "graph.geodesic"
        )
        expect_equal(op$settings$graph.stage, stage)
        expect_equal(op$settings$graph.summary$constructor,
                     paste0("graph.", stage))
        expect_equal(op$settings$support.metric, "graph.geodesic")
    }
})

test_that("lpl.tf.operator supports deterministic local PCA charts", {
    X <- rbind(
        c(0.00, 0.00), c(0.15, 0.01), c(0.31, 0.04), c(0.47, 0.10),
        c(0.63, 0.18), c(0.79, 0.27), c(0.92, 0.38), c(1.00, 0.50)
    )
    op1 <- lpl.tf.operator(
        X,
        degree = 1L,
        support.type = "knn",
        support.size = 5L,
        support.buffer = 0L,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        kernel = "gaussian"
    )
    op2 <- lpl.tf.operator(
        X,
        degree = 1L,
        support.type = "knn",
        support.size = 5L,
        support.buffer = 0L,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        kernel = "gaussian"
    )

    expect_equal(op1$settings$coordinate.method, "local.pca")
    expect_equal(op1$settings$chart.dim, 2L)
    expect_true(all(op1$row.table$coordinate.method == "local.pca"))
    expect_equal(as.matrix(op1$A), as.matrix(op2$A), tolerance = 1e-12)
    expect_true(all(op1$row.table$status == "ok"))
})

test_that("lpl.tf.operator knn keeps target in duplicate-coordinate support", {
    X <- matrix(c(0, 0, 0, 0, 0, 1, 2, 3), ncol = 1)
    op <- lpl.tf.operator(
        X,
        degree = 2L,
        support.type = "knn",
        support.size = 5L,
        support.buffer = 0L,
        anchor.index = 5L,
        kernel = "gaussian",
        duplicate.action = "keep"
    )
    support <- op$supports[[1L]]
    expect_true(5L %in% support$candidate)
    expect_false(5L %in% support$predictors)
    expect_equal(length(support$candidate), 5L)
    expect_equal(op$row.table$support.size, length(support$candidate))
    expect_equal(op$row.table$predictor.size, length(support$predictors))
})

test_that("lpl.tf.operator row normalization is complete-row normalization", {
    X <- matrix(0:3, ncol = 1)
    op.none <- lpl.tf.operator(
        X, degree = 2L, support.type = "knn", support.size = 4L,
        support.buffer = 0L, anchor.index = 1L, kernel = "gaussian",
        row.normalize = "none"
    )
    op.l2 <- lpl.tf.operator(
        X, degree = 2L, support.type = "knn", support.size = 4L,
        support.buffer = 0L, anchor.index = 1L, kernel = "gaussian",
        row.normalize = "l2"
    )
    raw <- as.numeric(op.none$A[1, ])
    normed <- raw / sqrt(sum(raw^2))
    expect_equal(as.numeric(op.l2$A[1, ]), normed, tolerance = 1e-12)
    expect_equal(op.l2$row.table$row.norm.final, 1, tolerance = 1e-12)
})

test_that("lpl.tf.operator reports l1 row-normalization metadata", {
    X <- matrix(0:3, ncol = 1)
    op.l1 <- lpl.tf.operator(
        X, degree = 2L, support.type = "knn", support.size = 4L,
        support.buffer = 0L, anchor.index = 1L, kernel = "gaussian",
        row.normalize = "l1"
    )
    raw <- c(1, -3, 3, -1)
    expect_equal(as.numeric(op.l1$A[1, ]), raw / sum(abs(raw)),
                 tolerance = 1e-12)
    expect_equal(op.l1$row.table$row.norm.raw, sum(abs(raw)),
                 tolerance = 1e-12)
    expect_equal(op.l1$row.table$row.norm.raw.l2, sqrt(sum(raw^2)),
                 tolerance = 1e-12)
    expect_equal(op.l1$row.table$row.norm.raw.l1, sum(abs(raw)),
                 tolerance = 1e-12)
    expect_equal(op.l1$row.table$row.norm.used, "l1")
    expect_equal(op.l1$row.table$row.norm.final, 1, tolerance = 1e-12)
})

test_that("fit.lpl.tf fixed lambda matches genlasso on fixed operator", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    X <- matrix(x, ncol = 1)
    y <- sin(2 * pi * x)
    op <- lpl.tf.operator(
        X, degree = 2L, support.type = "knn", support.size = 4L,
        support.buffer = 0L, kernel = "gaussian", row.normalize = "l2"
    )
    lambda <- 0.02
    fit <- fit.lpl.tf(
        y = y,
        operator = op,
        lambda = lambda,
        maxsteps = 500L
    )
    solver.op <- .lpl.tf.prepare.solver.operator(op)
    ref <- genlasso::genlasso(y = y, D = as.matrix(solver.op$operator$A),
                              maxsteps = 500L,
                              svd = TRUE)
    beta.ref <- as.vector(stats::coef(ref, lambda = lambda)$beta)
    expect_s3_class(fit, "lpl_tf")
    expect_equal(fit$fitted.values, beta.ref, tolerance = 1e-8)
    expect_equal(as.numeric(predict(fit)), fit$fitted.values)
    expect_equal(fit$lambda, lambda)
    expect_equal(fit$lambda.selection, "fixed")
    expect_equal(fit$diagnostics$lambda.boundary, "none")
    expect_true(is.na(fit$diagnostics$solver.converged))
    expect_equal(fit$diagnostics$solver.status,
                 "convergence_not_reported_by_backend")
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.lpl.tf CV uses supplied foldid deterministically", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 15)
    X <- matrix(x, ncol = 1)
    y <- cos(2 * pi * x)
    foldid <- rep(1:3, length.out = length(y))
    lambda.grid <- c(0.1, 0.01, 0.001)
    op <- lpl.tf.operator(
        X, degree = 2L, support.type = "knn", support.size = 4L,
        support.buffer = 0L, kernel = "gaussian", row.normalize = "l2"
    )
    fit1 <- fit.lpl.tf(
        y = y,
        operator = op,
        lambda.selection = "cv",
        lambda.grid = lambda.grid,
        foldid = foldid,
        cv.loss = "mse",
        maxsteps = 500L
    )
    fit2 <- fit.lpl.tf(
        y = y,
        operator = op,
        lambda.selection = "cv",
        lambda.grid = lambda.grid,
        foldid = foldid,
        cv.loss = "mse",
        maxsteps = 500L
    )
    expect_equal(fit1$fitted.values, fit2$fitted.values, tolerance = 1e-12)
    expect_equal(fit1$lambda, fit2$lambda)
    expect_equal(fit1$foldid, foldid)
    expect_equal(fit1$fold.source, "supplied")
    expect_true(fit1$lambda %in% lambda.grid)
    expect_true(fit1$diagnostics$lambda.boundary %in% c("lower", "upper", "none"))
    expect_true(all(c("mean.error", "se", "selected.idx") %in% names(fit1$cv)))
})

test_that("fit.lpl.tf requires explicit CV lambda grid in Phase 2", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    op <- lpl.tf.operator(
        matrix(x, ncol = 1), degree = 2L, support.type = "knn",
        support.size = 4L, support.buffer = 0L
    )
    expect_error(
        fit.lpl.tf(
            y = x,
            operator = op,
            lambda.selection = "cv",
            foldid = rep(1:3, length.out = length(x))
        ),
        "explicit 'lambda.grid'"
    )
})

test_that("fit.lpl.tf records operator-grid selection provenance", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 15)
    y <- cos(2 * pi * x)
    grid <- data.frame(
        support.type = "knn",
        support.size = c(3L, 4L, 5L),
        kernel = "gaussian",
        row.normalize = "l2",
        stringsAsFactors = FALSE
    )
    fit <- fit.lpl.tf(
        X = matrix(x, ncol = 1),
        y = y,
        degree = 2L,
        support.buffer = 0L,
        operator.grid = grid,
        lambda.selection = "cv",
        lambda.grid = c(0.1, 0.01, 0.001),
        foldid = rep(1:3, length.out = length(y)),
        maxsteps = 500L
    )
    sel <- fit$operator.selection
    expect_s3_class(fit, "lpl_tf")
    expect_true(is.list(sel))
    expect_equal(nrow(sel$table), 3L)
    expect_true(any(sel$table$status == "failed"))
    expect_true(any(sel$table$status == "ok"))
    expect_equal(sum(sel$table$selected), 1L)
    expect_true(sel$selected.candidate.id %in% sel$table$candidate.id)
    expect_true(is.finite(sel$score))
    expect_equal(sel$foldid, rep(1:3, length.out = length(y)))
    expect_equal(rownames(sel$table), as.character(seq_len(nrow(sel$table))))
})

test_that("fit.lpl.tf reports duplicate-row keying metadata", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 15)
    op <- lpl.tf.operator(
        matrix(x, ncol = 1), degree = 2L, support.type = "knn",
        support.size = 4L, support.buffer = 0L
    )
    fit <- fit.lpl.tf(y = sin(2 * pi * x), operator = op, lambda = 0.01,
                      maxsteps = 500L)
    expect_true(fit$diagnostics$solver.duplicate.rows.collapsed >= 0L)
    solver.op <- .lpl.tf.prepare.solver.operator(op)
    expect_equal(solver.op$diagnostics$duplicate.key.rule, "signif_format")
    expect_equal(solver.op$diagnostics$duplicate.key.digits, 14L)
})

test_that("refit.lpl.tf reuses operator and selected lambda", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 12)
    X <- matrix(x, ncol = 1)
    y <- x^2
    fit <- fit.lpl.tf(
        X = X,
        y = y,
        degree = 2L,
        support.type = "knn",
        support.size = 4L,
        support.buffer = 0L,
        lambda = 0.01,
        maxsteps = 500L
    )
    y2 <- y + 0.1 * sin(2 * pi * x)
    refit <- refit.lpl.tf(fit, y = y2)
    expect_s3_class(refit, "lpl_tf")
    expect_identical(refit$operator, fit$operator)
    expect_equal(refit$lambda, fit$lambda)
    expect_equal(refit$y, y2)
    expect_true(all(is.finite(refit$fitted.values)))
})

test_that("fit.lpl.tf validates Phase 2 prediction and response contracts", {
    skip_if_not_installed("genlasso")
    x <- seq(0, 1, length.out = 8)
    op <- lpl.tf.operator(
        matrix(x, ncol = 1), degree = 2L, support.type = "knn",
        support.size = 4L, support.buffer = 0L
    )
    fit <- fit.lpl.tf(y = x, operator = op, lambda = 0.01, maxsteps = 500L)
    expect_error(predict(fit, newdata = matrix(0, ncol = 1)), "newdata")
    expect_error(fit.lpl.tf(y = cbind(x, x), operator = op, lambda = 0.01),
                 "one response")
    expect_error(refit.lpl.tf(fit, y = x, reuse.lambda = FALSE),
                 "Provide 'lambda'")
})

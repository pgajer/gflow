#' Construct A Synchronized LPL-TF Fixed Operator
#'
#' Builds the fixed S-LPL-TF operator pair. The LPL component
#' \eqn{A_{\mathrm{LPL}}} is the accepted self-excluded LPL-TF residual
#' operator from \code{\link{lpl.tf.operator}}. The synchronization component
#' \eqn{C_{\mathrm{sync}}} compares inclusive local polynomial prediction maps
#' over overlapping supports.
#'
#' @inheritParams lpl.tf.operator
#' @param sync.row.normalize Synchronization-row normalization. Phase S2 uses
#'   \code{"l2"} by default.
#' @param sync.min.norm Minimum raw synchronization-row norm. Rows below this
#'   threshold are dropped with metadata.
#'
#' @return A list of class \code{"slpl_tf_operator"} containing
#'   \code{A_LPL}, \code{C_sync}, row metadata, diagnostics, settings, and the
#'   embedded \code{"lpl_tf_operator"}.
#' @export
slpl.tf.operator <- function(
    X,
    adj.list = NULL,
    weight.list = NULL,
    graph = NULL,
    graph.stage = "final",
    anchor.index = NULL,
    anchor.coordinates = NULL,
    degree = 2L,
    support.type = c("adaptive.radius", "knn", "fixed.radius"),
    support.size = NULL,
    radius = NULL,
    min.support = NULL,
    support.buffer = 3L,
    kernel = c("epanechnikov", "triangular", "gaussian", "tricube"),
    coordinate.method = c("coordinates", "local.pca"),
    chart.dim = NULL,
    support.metric = c("auto", "coordinates", "graph.geodesic"),
    row.normalize = c("l2", "none", "l1"),
    local.solver = c("auto", "normal.equations", "qr", "svd"),
    normal.equations.max.condition = 1e8,
    duplicate.action = c("keep", "error"),
    drop.rank.deficient = TRUE,
    sync.row.normalize = c("l2"),
    sync.min.norm = 1e-12,
    verbose = FALSE,
    ...) {

    support.type <- match.arg(support.type, c("adaptive.radius", "knn", "fixed.radius"))
    kernel <- match.arg(kernel, c("epanechnikov", "triangular", "gaussian", "tricube"))
    coordinate.method <- match.arg(coordinate.method, c("coordinates", "local.pca"))
    support.metric <- match.arg(support.metric, c("auto", "coordinates", "graph.geodesic"))
    row.normalize <- match.arg(row.normalize, c("l2", "none", "l1"))
    local.solver <- match.arg(local.solver, c("auto", "normal.equations", "qr", "svd"))
    duplicate.action <- match.arg(duplicate.action, c("keep", "error"))
    sync.row.normalize <- match.arg(sync.row.normalize, c("l2"))
    if (!is.numeric(sync.min.norm) || length(sync.min.norm) != 1L ||
        !is.finite(sync.min.norm) || sync.min.norm <= 0) {
        stop("'sync.min.norm' must be a positive finite scalar.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    args <- .lpl.tf.validate.operator.args(
        X = X,
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        graph.stage = graph.stage,
        anchor.index = anchor.index,
        anchor.coordinates = anchor.coordinates,
        degree = degree,
        support.type = support.type,
        support.size = support.size,
        radius = radius,
        min.support = min.support,
        support.buffer = support.buffer,
        kernel = kernel,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        support.metric = support.metric,
        exclude.self = TRUE,
        row.normalize = row.normalize,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        duplicate.action = duplicate.action,
        drop.rank.deficient = drop.rank.deficient,
        verbose = verbose
    )

    lpl.operator <- lpl.tf.operator(
        X = X,
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        graph.stage = graph.stage,
        anchor.index = anchor.index,
        anchor.coordinates = anchor.coordinates,
        degree = degree,
        support.type = support.type,
        support.size = support.size,
        radius = radius,
        min.support = min.support,
        support.buffer = support.buffer,
        kernel = kernel,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        support.metric = support.metric,
        exclude.self = TRUE,
        row.normalize = row.normalize,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        duplicate.action = duplicate.action,
        drop.rank.deficient = drop.rank.deficient,
        verbose = verbose
    )
    sync <- .slpl.tf.sync.operator.from.lpl(
        lpl.operator = lpl.operator,
        args = args,
        sync.row.normalize = sync.row.normalize,
        sync.min.norm = sync.min.norm
    )
    diagnostics <- .slpl.tf.operator.diagnostics(lpl.operator, sync, X = args$X)
    out <- list(
        A_LPL = lpl.operator$A,
        C_sync = sync$C,
        X = args$X,
        lpl.operator = lpl.operator,
        sync.map.table = sync$map.table,
        sync.row.table = sync$row.table,
        sync.coverage = sync$coverage,
        diagnostics = diagnostics,
        settings = modifyList(lpl.operator$settings, list(
            sync.row.normalize = sync.row.normalize,
            sync.min.norm = sync.min.norm,
            sync.map.self.inclusion = "inclusive",
            lpl.residual.self.inclusion = "self_excluded",
            prediction.coverage = "Q_u_equals_P_u"
        )),
        call = match.call()
    )
    class(out) <- c("slpl_tf_operator", "list")
    out
}

#' @export
print.slpl_tf_operator <- function(x, ...) {
    cat("Synchronized Local Polynomial Lifting Trend-Filtering operator\n")
    cat("  observations:", ncol(x$A_LPL), "\n")
    cat("  LPL rows:", nrow(x$A_LPL), "\n")
    cat("  synchronization rows:", nrow(x$C_sync), "\n")
    cat("  degree:", x$settings$degree, "\n")
    cat("  support type:", x$settings$support.type, "\n")
    cat("  support metric:", x$settings$support.metric, "\n")
    cat("  stacked operator rank:", x$diagnostics$stacked.operator.rank, "\n")
    cat("  stacked operator nullity:", x$diagnostics$stacked.operator.nullity, "\n")
    invisible(x)
}

#' Fit A Fixed-Operator Synchronized LPL-TF Model
#'
#' Fits the fixed-operator S-LPL-TF estimator
#' \deqn{\hat f = \arg\min_f \frac12\|y-f\|_2^2 +
#' \lambda_1\|A_{\mathrm{LPL}}f\|_1 +
#' \frac{\lambda_2}{2}\|C_{\mathrm{sync}}f\|_2^2.}
#'
#' @param X Optional coordinate matrix used to build
#'   \code{\link{slpl.tf.operator}} when \code{operator} is \code{NULL}.
#' @param y Numeric response vector.
#' @param operator Optional \code{"slpl_tf_operator"} object.
#' @param lambda1 Fixed LPL-TF \eqn{\ell_1} penalty, or a shortcut for
#'   \code{lambda1.grid} when \code{lambda.selection = "cv"}.
#' @param lambda2 Fixed quadratic synchronization penalty, or a shortcut for
#'   \code{lambda2.grid} when \code{lambda.selection = "cv"}.
#' @param lambda1.grid Optional nonnegative \eqn{\lambda_1} grid for CV.
#' @param lambda2.grid Optional nonnegative \eqn{\lambda_2} grid for CV.
#' @param sync.lambda.grid Deprecated alias for \code{lambda2.grid}.
#' @param lambda.selection \code{"fixed"} or \code{"cv"}. CV uses materialized
#'   folds and a deterministic Cartesian grid over \code{lambda1.grid} and
#'   \code{lambda2.grid}.
#' @param operator.grid Optional Phase-S3 operator candidate grid. Supply a data
#'   frame with one row per candidate or a list of named lists.
#' @param foldid Optional deterministic fold assignments for CV.
#' @param cv.folds Number of generated folds when \code{foldid = NULL}.
#' @param cv.loss Cross-validation loss, \code{"mse"}, \code{"rmse"}, or
#'   \code{"mae"}.
#' @param cv.seed Optional seed for reproducible generated folds.
#' @param cv.repeats Phase S3 supports only one CV repeat.
#' @param solver Current implementation supports only \code{"genlasso"}.
#' @param selection Selection rule, \code{"min"} or \code{"one.se"}. For the
#'   two-parameter one-standard-error rule, the most regularized eligible pair
#'   is chosen by decreasing \code{lambda1} and then decreasing \code{lambda2}.
#' @param maxsteps,minlam,approx,rtol,btol,eps Genlasso controls.
#' @param verbose Logical.
#' @param ... Additional arguments passed to \code{\link{slpl.tf.operator}} when
#'   the operator is built internally.
#'
#' @return A list of class \code{"slpl_tf"}.
#' @export
fit.slpl.tf <- function(
    X = NULL,
    y,
    operator = NULL,
    lambda1 = NULL,
    lambda2 = 0,
    lambda1.grid = NULL,
    lambda2.grid = NULL,
    sync.lambda.grid = NULL,
    lambda.selection = c("fixed", "cv"),
    operator.grid = NULL,
    foldid = NULL,
    cv.folds = 5L,
    cv.loss = c("mse", "rmse", "mae"),
    cv.seed = NULL,
    cv.repeats = 1L,
    solver = c("genlasso"),
    selection = c("min", "one.se"),
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    verbose = FALSE,
    ...) {

    if (!requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for fit.slpl.tf().", call. = FALSE)
    }
    solver <- match.arg(solver)
    lambda.selection <- match.arg(lambda.selection)
    cv.loss <- match.arg(cv.loss)
    selection <- match.arg(selection)
    cv.folds <- .validate.ssrhe.positive.integer(cv.folds, "cv.folds")
    cv.repeats <- .validate.ssrhe.positive.integer(cv.repeats, "cv.repeats")
    if (cv.repeats != 1L) {
        stop("Phase S3 fit.slpl.tf() supports cv.repeats = 1 only.",
             call. = FALSE)
    }
    if (!is.null(sync.lambda.grid)) {
        if (!is.null(lambda2.grid)) {
            stop("Specify only one of 'lambda2.grid' or 'sync.lambda.grid'.",
                 call. = FALSE)
        }
        lambda2.grid <- sync.lambda.grid
    }
    if (identical(lambda.selection, "fixed")) {
        if (!is.null(lambda1.grid)) {
            stop("'lambda1.grid' is only used when lambda.selection = 'cv'.",
                 call. = FALSE)
        }
        if (!is.null(lambda2.grid)) {
            stop("'lambda2.grid' is only used when lambda.selection = 'cv'.",
                 call. = FALSE)
        }
        lambda1.grid <- .validate.ssrhe.hessian.l1.lambda.grid(lambda1, "fixed")
        lambda2.grid <- .slpl.tf.validate.lambda.grid(lambda2, "lambda2")
    } else {
        if (!is.null(lambda1.grid) && !is.null(lambda1)) {
            stop("Specify only one of 'lambda1' or 'lambda1.grid' for CV.",
                 call. = FALSE)
        }
        if (!is.null(lambda2.grid) && !is.null(lambda2) && length(lambda2) > 1L) {
            stop("Specify only one of vector 'lambda2' or 'lambda2.grid' for CV.",
                 call. = FALSE)
        }
        if (is.null(lambda1.grid)) {
            lambda1.grid <- lambda1
        }
        if (is.null(lambda1.grid)) {
            stop("'lambda1.grid' is required when lambda.selection = 'cv'.",
                 call. = FALSE)
        }
        if (is.null(lambda2.grid)) {
            lambda2.grid <- if (length(lambda2) > 1L) {
                lambda2
            } else {
                c(0, 1e-3, 1e-2, 1e-1, 1, 10)
            }
        }
        lambda1.grid <- .validate.ssrhe.hessian.l1.lambda.grid(lambda1.grid,
                                                               "cv")
        lambda2.grid <- .slpl.tf.validate.lambda.grid(lambda2.grid,
                                                      "lambda2.grid")
    }
    if (identical(lambda.selection, "fixed") && length(lambda2.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda2 value.",
             call. = FALSE)
    }
    operator.args <- list(...)
    if (!is.null(operator.grid)) {
        if (!is.null(operator)) {
            stop("Specify 'operator.grid' only when 'operator' is NULL.",
                 call. = FALSE)
        }
        if (is.null(X)) {
            stop("'X' is required when 'operator.grid' is supplied.",
                 call. = FALSE)
        }
        return(.slpl.tf.fit.operator.grid(
            X = X,
            y = y,
            operator.grid = operator.grid,
            base.operator.args = operator.args,
            fit.args = list(
                lambda1.grid = lambda1.grid,
                lambda2.grid = lambda2.grid,
                lambda.selection = lambda.selection,
                foldid = foldid,
                cv.folds = cv.folds,
                cv.loss = cv.loss,
                cv.seed = cv.seed,
                cv.repeats = cv.repeats,
                solver = solver,
                selection = selection,
                maxsteps = maxsteps,
                minlam = minlam,
                approx = approx,
                rtol = rtol,
                btol = btol,
                eps = eps,
                verbose = verbose
            ),
            call = match.call()
        ))
    }
    operator <- .slpl.tf.prepare.fit.operator(
        X = X,
        operator = operator,
        ...
    )
    n <- ncol(operator$A_LPL)
    y <- .lpl.tf.validate.response(y, n)
    weights <- rep(1, n)
    fold.source <- if (identical(lambda.selection, "fixed")) {
        NA_character_
    } else if (is.null(foldid)) {
        if (is.null(cv.seed)) "generated_deterministic" else "generated_seeded"
    } else {
        "supplied"
    }
    foldid <- .lpl.tf.prepare.foldid(
        foldid = foldid,
        n = n,
        cv.folds = cv.folds,
        cv.seed = cv.seed,
        lambda.selection = lambda.selection
    )
    solver.args <- .validate.ssrhe.hessian.l1.solver.args(
        n.lambda = length(lambda1.grid),
        nfolds = cv.folds,
        fold.id = foldid,
        observed = rep(TRUE, n),
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = "none",
        verbose = verbose
    )
    loss.internal <- if (identical(cv.loss, "rmse")) "mse" else cv.loss
    cv <- NULL
    solver.operator <- .lpl.tf.prepare.solver.operator(operator$lpl.operator)
    if (identical(lambda.selection, "cv")) {
        cv <- .slpl.tf.cv(
            y = y,
            weights = weights,
            operator = operator,
            solver.lpl.operator = solver.operator$operator,
            lambda1.grid = lambda1.grid,
            lambda2.grid = lambda2.grid,
            fold.id = solver.args$fold.id,
            loss = loss.internal,
            selection = selection,
            solver.args = solver.args
        )
        cv <- .slpl.tf.decorate.cv(cv, requested.loss = cv.loss)
        lambda1 <- cv$lambda1
        lambda2 <- cv$lambda2
    } else {
        lambda1 <- lambda1.grid[[1L]]
        lambda2 <- lambda2.grid[[1L]]
    }
    start <- proc.time()[["elapsed"]]
    raw.fit <- .slpl.tf.fit.from.operator(
        operator = operator,
        solver.lpl.operator = solver.operator$operator,
        y = y,
        lambda1 = lambda1,
        lambda2 = lambda2,
        solver.args = solver.args
    )
    runtime <- proc.time()[["elapsed"]] - start
    fitted <- raw.fit$fitted.values
    residuals <- y - fitted
    lpl.values <- as.vector(operator$A_LPL %*% fitted)
    sync.values <- as.vector(operator$C_sync %*% fitted)
    diagnostics <- modifyList(raw.fit$diagnostics, list(
        fit.status = "ok",
        runtime.seconds = runtime,
        lambda1 = lambda1,
        lambda2 = lambda2,
        n.lpl.rows = nrow(operator$A_LPL),
        n.sync.rows = nrow(operator$C_sync),
        n.active.lpl.rows = sum(abs(lpl.values) > 1e-8),
        sync.energy.l2 = sum(sync.values^2),
        data.loss = 0.5 * sum(residuals^2),
        lpl.energy.l1 = sum(abs(lpl.values)),
        objective = 0.5 * sum(residuals^2) +
            lambda1 * sum(abs(lpl.values)) +
            0.5 * lambda2 * sum(sync.values^2),
        lambda1.boundary = .lpl.tf.lambda.boundary(lambda1.grid,
                                                   if (is.null(cv)) 1L else cv$lambda1.selected.idx),
        lambda2.boundary = .lpl.tf.lambda.boundary(lambda2.grid,
                                                   if (is.null(cv)) 1L else cv$lambda2.selected.idx),
        solver.duplicate.rows.collapsed =
            solver.operator$diagnostics$duplicate.rows.collapsed,
        solver.operator.rows = solver.operator$diagnostics$solver.rows,
        solver.operator.raw.rows = solver.operator$diagnostics$raw.rows
    ))
    out <- list(
        fitted.values = fitted,
        residuals = residuals,
        y = y,
        X = X,
        operator = operator,
        lambda1 = lambda1,
        lambda2 = lambda2,
        lambda1.grid = lambda1.grid,
        lambda2.grid = lambda2.grid,
        lambda.selection = lambda.selection,
        cv = cv,
        foldid = solver.args$fold.id,
        fold.source = fold.source,
        cv.loss = cv.loss,
        cv.loss.internal = loss.internal,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        selection = selection,
        beta.grid = raw.fit$beta.grid,
        path = raw.fit$path,
        solver = raw.fit$solver,
        diagnostics = diagnostics,
        call = match.call()
    )
    class(out) <- c("slpl_tf", "list")
    out
}

#' Refit A Fixed-Operator Synchronized LPL-TF Model
#'
#' @param object A \code{"slpl_tf"} fit.
#' @param y New numeric response vector.
#' @param lambda1,lambda2 Optional fixed penalties. Defaults reuse
#'   \code{object$lambda1} and \code{object$lambda2}.
#' @param reuse.lambda Logical. If \code{TRUE}, reuse omitted penalties.
#' @param verbose Logical.
#' @param ... Reserved.
#'
#' @return A refitted \code{"slpl_tf"} object.
#' @export
refit.slpl.tf <- function(object, y, lambda1 = NULL, lambda2 = NULL,
                          reuse.lambda = TRUE, verbose = FALSE, ...) {
    if (!inherits(object, "slpl_tf")) {
        stop("'object' must be an slpl_tf fit.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    if (is.null(lambda1)) {
        if (!isTRUE(reuse.lambda)) {
            stop("Provide 'lambda1' or set reuse.lambda = TRUE.", call. = FALSE)
        }
        lambda1 <- object$lambda1
    }
    if (is.null(lambda2)) {
        if (!isTRUE(reuse.lambda)) {
            stop("Provide 'lambda2' or set reuse.lambda = TRUE.", call. = FALSE)
        }
        lambda2 <- object$lambda2
    }
    fit.slpl.tf(
        X = object$X,
        y = y,
        operator = object$operator,
        lambda1 = lambda1,
        lambda2 = lambda2,
        solver = object$solver$requested %||% "genlasso",
        verbose = verbose
    )
}

#' Predict From A Fixed-Operator Synchronized LPL-TF Model
#'
#' @param object A \code{"slpl_tf"} fit.
#' @param newdata Must be \code{NULL} in Phase S2.
#' @param type Prediction type. Phase S2 supports only \code{"response"}.
#' @param ... Reserved.
#'
#' @return Fitted values at the training points.
#' @method predict slpl_tf
#' @export
predict.slpl_tf <- function(object, newdata = NULL, type = c("response"),
                            ...) {
    type <- match.arg(type)
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    if (!is.null(newdata)) {
        stop("Phase S2 predict() for slpl_tf supports training-point ",
             "prediction only; 'newdata' must be NULL.", call. = FALSE)
    }
    object$fitted.values
}

.slpl.tf.sync.operator.from.lpl <- function(lpl.operator, args,
                                            sync.row.normalize,
                                            sync.min.norm) {
    X <- args$X
    n <- nrow(X)
    powers <- .lpl.tf.monomial.powers(args$chart.dim, args$degree)
    design.columns <- nrow(powers)
    coordinate.dist.matrix <- as.matrix(stats::dist(X))
    support.dist.matrix <- if (identical(args$support.metric, "graph.geodesic")) {
        args$graph.distance.matrix
    } else {
        coordinate.dist.matrix
    }
    maps <- list()
    map.table <- list()
    map.id <- 0L
    for (aa in seq_along(args$anchor.index)) {
        anchor <- args$anchor.index[[aa]]
        support <- .slpl.tf.inclusive.support(
            target = anchor,
            distances = support.dist.matrix[anchor, ],
            support.type = args$support.type,
            support.size = args$support.size,
            radius = args$radius,
            min.support = args$min.support
        )
        chart <- .lpl.tf.local.chart(
            X = X,
            target = anchor,
            candidate = support$candidate,
            coordinate.method = args$coordinate.method,
            chart.dim = args$chart.dim
        )
        D <- .lpl.tf.design.matrix(chart$coordinates, powers)
        weights <- .slpl.tf.inclusive.kernel.weights(
            distances = support$distances,
            kernel = args$kernel
        )
        positive <- is.finite(weights) & weights > 0
        support.positive <- support$candidate[positive]
        weights.positive <- weights[positive]
        D.positive <- D[positive, , drop = FALSE]
        for (target in support$candidate) {
            map.id <- map.id + 1L
            local.idx <- match(target, support$candidate)
            target.basis <- .lpl.tf.design.matrix(
                matrix(chart$coordinates[local.idx, ], nrow = 1L),
                powers
            )
            fit <- .lpl.tf.prediction.weights(
                D = D.positive,
                target.basis = target.basis,
                weights = weights.positive,
                local.solver = args$local.solver,
                normal.equations.max.condition =
                    args$normal.equations.max.condition
            )
            h <- rep(NA_real_, n)
            status <- "dropped"
            if (fit$ok && fit$rank >= design.columns) {
                h <- numeric(n)
                h[support.positive] <- fit$h
                status <- "ok"
            }
            maps[[map.id]] <- list(
                map.id = map.id,
                anchor = anchor,
                target = target,
                support = support.positive,
                h = h,
                status = status
            )
            map.table[[map.id]] <- data.frame(
                map.id = map.id,
                anchor.index = anchor,
                target.index = target,
                support.indices = paste(support.positive, collapse = ","),
                support.size = length(support$candidate),
                positive.weight.support.size = length(support.positive),
                prediction.coverage.set = paste(support$candidate, collapse = ","),
                support.type = args$support.type,
                support.metric = args$support.metric,
                kernel = args$kernel,
                coordinate.method = args$coordinate.method,
                chart.dim = args$chart.dim,
                degree = args$degree,
                design.columns = design.columns,
                weighted.design.rank = fit$rank,
                condition.number = fit$condition.number,
                solver.used = fit$solver.used,
                fallback.used = fit$fallback.used,
                status = status,
                drop.reason = if (identical(status, "ok")) {
                    NA_character_
                } else {
                    fit$drop.reason %||% "rank_deficient_requested_degree"
                },
                sync_map_self_inclusion = "inclusive",
                lpl_residual_self_inclusion = "self_excluded",
                stringsAsFactors = FALSE
            )
        }
    }
    map.table <- do.call(rbind, map.table)
    row.out <- .slpl.tf.sync.rows.from.maps(
        maps = maps,
        n = n,
        sync.row.normalize = sync.row.normalize,
        sync.min.norm = sync.min.norm
    )
    feasible.maps <- Filter(function(m) identical(m$status, "ok"), maps)
    attempted.coverage <- tabulate(
        vapply(maps, `[[`, integer(1), "target"),
        nbins = n
    )
    feasible.coverage <- if (length(feasible.maps)) {
        tabulate(vapply(feasible.maps, `[[`, integer(1), "target"),
                 nbins = n)
    } else {
        integer(n)
    }
    list(
        C = row.out$C,
        maps = maps,
        map.table = map.table,
        row.table = row.out$row.table,
        coverage = list(
            attempted = attempted.coverage,
            feasible = feasible.coverage,
            zero = sum(feasible.coverage == 0L),
            one = sum(feasible.coverage == 1L),
            two.plus = sum(feasible.coverage >= 2L)
        )
    )
}

.slpl.tf.inclusive.support <- function(target, distances, support.type,
                                       support.size, radius, min.support) {
    n <- length(distances)
    ord <- order(distances, seq_len(n), na.last = NA)
    if (identical(support.type, "knn")) {
        candidate <- head(ord, support.size)
    } else if (identical(support.type, "fixed.radius")) {
        candidate <- which(distances <= radius)
        candidate <- candidate[order(distances[candidate], candidate)]
    } else {
        candidate <- head(ord, min.support)
    }
    list(
        target = target,
        candidate = as.integer(candidate),
        distances = as.double(distances[candidate]),
        radius = if (length(candidate)) max(distances[candidate]) else NA_real_
    )
}

.slpl.tf.inclusive.kernel.weights <- function(distances, kernel) {
    if (!length(distances)) return(numeric(0))
    positive.dist <- distances[is.finite(distances) & distances > 0]
    bandwidth <- if (length(positive.dist)) max(positive.dist) else 0
    if (!is.finite(bandwidth) || bandwidth <= 0) {
        u <- rep(0, length(distances))
    } else {
        u <- distances / (bandwidth + sqrt(.Machine$double.eps))
    }
    w <- switch(
        kernel,
        epanechnikov = pmax(0, 1 - u^2),
        triangular = pmax(0, 1 - u),
        tricube = ifelse(u < 1, (1 - u^3)^3, 0),
        gaussian = exp(-0.5 * u^2)
    )
    w[!is.finite(w)] <- 0
    as.double(w)
}

.slpl.tf.sync.rows.from.maps <- function(maps, n, sync.row.normalize,
                                         sync.min.norm) {
    rows <- list()
    row.table <- list()
    row.id <- 0L
    for (target in seq_len(n)) {
        candidates <- Filter(function(m) {
            identical(m$status, "ok") && identical(m$target, target)
        }, maps)
        if (length(candidates) < 2L) next
        for (ii in seq_len(length(candidates) - 1L)) {
            for (jj in (ii + 1L):length(candidates)) {
                raw <- candidates[[ii]]$h - candidates[[jj]]$h
                normed <- .lpl.tf.normalize.row(raw, sync.row.normalize,
                                                min.norm = sync.min.norm)
                row.id <- row.id + 1L
                status <- if (normed$ok) "ok" else "dropped"
                if (normed$ok) rows[[length(rows) + 1L]] <- normed$row
                row.table[[row.id]] <- data.frame(
                    sync.row.id = row.id,
                    target.index = target,
                    anchor.left = candidates[[ii]]$anchor,
                    anchor.right = candidates[[jj]]$anchor,
                    map.left = candidates[[ii]]$map.id,
                    map.right = candidates[[jj]]$map.id,
                    left.support = paste(candidates[[ii]]$support, collapse = ","),
                    right.support = paste(candidates[[jj]]$support, collapse = ","),
                    status = status,
                    drop.reason = if (normed$ok) NA_character_ else normed$drop.reason,
                    row.norm.raw = normed$row.norm.raw,
                    row.norm.raw.l2 = normed$row.norm.raw.l2,
                    row.norm.raw.l1 = normed$row.norm.raw.l1,
                    row.norm.used = normed$row.norm.used,
                    row.norm.final = normed$row.norm.final,
                    weight.before.normalization = 1,
                    post.normalization.weight = 1,
                    scalar.weight.effective = "none_equal_weight_unit_normalized",
                    sync_map_self_inclusion = "inclusive",
                    lpl_residual_self_inclusion = "self_excluded",
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    row.table <- if (length(row.table)) do.call(rbind, row.table) else {
        data.frame()
    }
    C <- .lpl.tf.assemble.sparse(rows, ncol = n)
    list(C = C, row.table = row.table)
}

.slpl.tf.operator.diagnostics <- function(lpl.operator, sync, X) {
    A <- lpl.operator$A
    C <- sync$C
    stacked <- rbind(A, C)
    stack.rank <- .lpl.tf.matrix.rank(stacked)
    list(
        n.observations = ncol(A),
        n.lpl.rows = nrow(A),
        n.sync.rows = nrow(C),
        attempted.sync.rows = nrow(sync$row.table),
        valid.sync.rows = nrow(C),
        dropped.sync.rows = if (nrow(sync$row.table)) {
            sum(sync$row.table$status != "ok")
        } else {
            0L
        },
        sync.coverage.zero = sync$coverage$zero,
        sync.coverage.one = sync$coverage$one,
        sync.coverage.two.plus = sync$coverage$two.plus,
        sync.rows.per.target.mean = if (nrow(sync$row.table)) {
            mean(tabulate(as.integer(sync$row.table$target.index),
                          nbins = ncol(A)))
        } else {
            0
        },
        sync.rows.per.n = nrow(C) / ncol(A),
        sync.row.norm.raw.summary = if (nrow(sync$row.table)) {
            summary(sync$row.table$row.norm.raw)
        } else {
            summary(numeric(0))
        },
        lpl.operator.rank = lpl.operator$diagnostics$operator.rank,
        lpl.operator.nullity = lpl.operator$diagnostics$operator.nullity,
        stacked.operator.rank = stack.rank,
        stacked.operator.nullity = if (is.na(stack.rank)) NA_integer_ else ncol(A) - stack.rank,
        lpl.polynomial.reproduction.error =
            lpl.operator$diagnostics$polynomial.residual,
        sync.polynomial.reproduction.error =
            .slpl.tf.polynomial.residual(C, lpl.operator, X)
    )
}

.slpl.tf.polynomial.residual <- function(C, lpl.operator, X) {
    settings <- lpl.operator$settings
    if (!identical(settings$coordinate.method, "coordinates")) {
        return(NA_real_)
    }
    probes <- .lpl.tf.design.matrix(
        X,
        .lpl.tf.monomial.powers(ncol(X), settings$degree)
    )
    if (nrow(C) == 0L || ncol(probes) == 0L) return(NA_real_)
    norm(as.matrix(C %*% probes), type = "F") /
        max(norm(probes, type = "F"), .Machine$double.eps)
}

.slpl.tf.prepare.fit.operator <- function(X, operator, ...) {
    if (!is.null(operator)) {
        if (!inherits(operator, "slpl_tf_operator")) {
            stop("'operator' must be an slpl_tf_operator object.", call. = FALSE)
        }
        if (is.null(operator$A_LPL) || is.null(operator$C_sync)) {
            stop("'operator' must contain A_LPL and C_sync.", call. = FALSE)
        }
        return(operator)
    }
    if (is.null(X)) {
        stop("Provide 'operator' or 'X' to fit.slpl.tf().", call. = FALSE)
    }
    slpl.tf.operator(X = X, ...)
}

.slpl.tf.fit.from.operator <- function(operator, solver.lpl.operator, y,
                                       lambda1, lambda2, solver.args) {
    n <- ncol(operator$A_LPL)
    path.info <- .slpl.tf.genlasso.path(
        operator = operator,
        solver.lpl.operator = solver.lpl.operator,
        y = y,
        weights = rep(1, n),
        lambda2 = lambda2,
        solver.args = solver.args
    )
    beta.grid <- .ssrhe.hessian.l1.coef.matrix(path.info$path, lambda1, n)
    if (any(!is.finite(beta.grid))) {
        stop("S-LPL-TF genlasso path produced non-finite fitted values.",
             call. = FALSE)
    }
    fitted <- as.vector(beta.grid[, 1L])
    list(
        fitted.values = fitted,
        beta.grid = beta.grid,
        path = path.info$path,
        diagnostics = path.info$diagnostics,
        solver = list(
            backend = "genlasso_augmented_design",
            requested = solver.args$solver,
            representation = path.info$D.info$representation,
            svd = path.info$D.info$svd,
            row.scaling = path.info$D.info$row.scaling,
            augmented.design.dim = paste(path.info$design.dim, collapse = " x "),
            augmented.response.length = path.info$response.length,
            admm = NULL
        )
    )
}

.slpl.tf.genlasso.path <- function(operator, solver.lpl.operator, y, weights,
                                   lambda2, solver.args) {
    n <- ncol(operator$A_LPL)
    D.info <- .ssrhe.hessian.l1.penalty.matrix(
        solver.lpl.operator$A,
        row.scaling = solver.args$row.scaling
    )
    diagnostics <- .ssrhe.hessian.l1.diagnostics(
        A = solver.lpl.operator$A,
        D.info = D.info
    )
    observed <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(observed)) {
        stop("At least one observed positive-weight response is required.",
             call. = FALSE)
    }
    obs <- which(observed)
    sw <- sqrt(weights[obs])
    if (lambda2 > 0 && nrow(operator$C_sync) > 0L) {
        w <- numeric(n)
        w[obs] <- weights[obs]
        y.clean <- numeric(n)
        y.clean[obs] <- y[obs]
        b <- w * y.clean
        Q <- Matrix::Diagonal(n, x = w) +
            lambda2 * Matrix::crossprod(operator$C_sync)
        Q <- as.matrix(Matrix::forceSymmetric(Q))
        chol.info <- tryCatch(
            list(R = chol(Q), status = "ok", message = NA_character_),
            error = function(e) {
                list(R = NULL, status = "failed",
                     message = conditionMessage(e))
            }
        )
        if (!identical(chol.info$status, "ok")) {
            stop("Positive-lambda2 compressed design Cholesky failed: ",
                 chol.info$message, call. = FALSE)
        }
        X.design <- chol.info$R
        y.design <- forwardsolve(t(chol.info$R), b)
        design.representation <- "compressed_quadratic_design"
    } else {
        X.obs <- as.matrix(Matrix::sparseMatrix(
            i = seq_along(obs),
            j = obs,
            x = sw,
            dims = c(length(obs), n),
            giveCsparse = TRUE
        ))
        X.design <- X.obs
        y.design <- sw * y[obs]
        design.representation <- "observed_identity_rows"
    }
    path <- suppressWarnings(genlasso::genlasso(
        y = y.design,
        X = X.design,
        D = D.info$D,
        approx = solver.args$approx,
        maxsteps = solver.args$maxsteps,
        minlam = solver.args$minlam,
        rtol = solver.args$rtol,
        btol = solver.args$btol,
        eps = solver.args$eps,
        verbose = solver.args$verbose,
        svd = D.info$svd
    ))
    list(
        path = path,
        D.info = D.info,
        diagnostics = diagnostics,
        design.dim = dim(X.design),
        response.length = length(y.design),
        n.observed = sum(observed),
        design.representation = design.representation
    )
}

.slpl.tf.cv <- function(y, weights, operator, solver.lpl.operator,
                        lambda1.grid, lambda2.grid, fold.id, loss,
                        selection, solver.args) {
    n <- length(y)
    folds <- sort(unique(fold.id[fold.id > 0L]))
    n1 <- length(lambda1.grid)
    n2 <- length(lambda2.grid)
    cv.errors <- array(NA_real_, dim = c(length(folds), n1, n2),
                       dimnames = list(
                           paste0("Fold", folds),
                           format(signif(lambda1.grid, 6), scientific = TRUE),
                           format(signif(lambda2.grid, 6), scientific = TRUE)
                       ))
    fold.table <- vector("list", length(folds) * n2)
    tt <- 0L
    for (kk in seq_along(lambda2.grid)) {
        lambda2 <- lambda2.grid[[kk]]
        for (ii in seq_along(folds)) {
            fold <- folds[[ii]]
            test <- which(fold.id == fold)
            train.weights <- weights
            train.weights[test] <- 0
            path.info <- tryCatch(
                .slpl.tf.genlasso.path(
                    operator = operator,
                    solver.lpl.operator = solver.lpl.operator,
                    y = y,
                    weights = train.weights,
                    lambda2 = lambda2,
                    solver.args = solver.args
                ),
                error = function(e) e
            )
            tt <- tt + 1L
            if (inherits(path.info, "error")) {
                fold.table[[tt]] <- data.frame(
                    lambda2.index = kk,
                    lambda2 = lambda2,
                    fold = fold,
                    status = "failed",
                    failure.message = conditionMessage(path.info),
                    stringsAsFactors = FALSE
                )
                next
            }
            beta <- .ssrhe.hessian.l1.coef.matrix(path.info$path,
                                                  lambda1.grid, n)
            beta.test <- beta[test, , drop = FALSE]
            finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
            target <- matrix(y[test], nrow = length(test), ncol = n1)
            w.test <- matrix(weights[test], nrow = length(test), ncol = n1)
            if (identical(loss, "mse")) {
                fold.error <- colSums(w.test * (target - beta.test)^2) /
                    colSums(w.test)
            } else {
                fold.error <- colSums(w.test * abs(target - beta.test)) /
                    colSums(w.test)
            }
            fold.error[!finite.pred] <- NA_real_
            cv.errors[ii, , kk] <- fold.error
            fold.table[[tt]] <- data.frame(
                lambda2.index = kk,
                lambda2 = lambda2,
                fold = fold,
                status = "ok",
                failure.message = NA_character_,
                stringsAsFactors = FALSE
            )
        }
    }
    fold.table <- do.call(rbind, fold.table)
    mean.error <- apply(cv.errors, c(2L, 3L), mean, na.rm = TRUE)
    mean.error[!is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, c(2L, 3L), function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    candidates <- expand.grid(
        lambda1.index = seq_along(lambda1.grid),
        lambda2.index = seq_along(lambda2.grid),
        KEEP.OUT.ATTRS = FALSE
    )
    candidates$lambda1 <- lambda1.grid[candidates$lambda1.index]
    candidates$lambda2 <- lambda2.grid[candidates$lambda2.index]
    candidates$mean.error <- mapply(function(i, j) mean.error[i, j],
                                    candidates$lambda1.index,
                                    candidates$lambda2.index)
    candidates$se <- mapply(function(i, j) se[i, j],
                            candidates$lambda1.index,
                            candidates$lambda2.index)
    candidates$failed.folds <- vapply(seq_len(nrow(candidates)), function(ii) {
        sum(!is.finite(cv.errors[, candidates$lambda1.index[[ii]],
                                candidates$lambda2.index[[ii]]]))
    }, integer(1))
    candidates$status <- ifelse(is.finite(candidates$mean.error), "ok", "failed")
    candidates$selected <- FALSE
    best.idx <- which.min(candidates$mean.error)
    if (!length(best.idx) || !is.finite(candidates$mean.error[[best.idx]])) {
        stop("All S-LPL-TF CV candidates failed for the supplied grids.",
             call. = FALSE)
    }
    if (identical(selection, "one.se")) {
        threshold <- candidates$mean.error[[best.idx]] + candidates$se[[best.idx]]
        eligible <- which(candidates$mean.error <= threshold)
        ord <- order(-candidates$lambda1[eligible],
                     -candidates$lambda2[eligible],
                     candidates$mean.error[eligible])
        selected.idx <- eligible[ord[[1L]]]
    } else {
        selected.idx <- best.idx
    }
    candidates$selected[[selected.idx]] <- TRUE
    list(
        fold.id = fold.id,
        lambda1.grid = lambda1.grid,
        lambda2.grid = lambda2.grid,
        fold.errors = cv.errors,
        mean.error = mean.error,
        se = se,
        candidate.table = candidates,
        fold.fit.table = fold.table,
        best.idx = best.idx,
        selected.idx = selected.idx,
        lambda1.selected.idx = candidates$lambda1.index[[selected.idx]],
        lambda2.selected.idx = candidates$lambda2.index[[selected.idx]],
        lambda1.min = candidates$lambda1[[best.idx]],
        lambda2.min = candidates$lambda2[[best.idx]],
        lambda1 = candidates$lambda1[[selected.idx]],
        lambda2 = candidates$lambda2[[selected.idx]],
        error.min = candidates$mean.error[[best.idx]],
        error.selected = candidates$mean.error[[selected.idx]],
        selection = selection,
        loss = loss
    )
}

.slpl.tf.decorate.cv <- function(cv, requested.loss) {
    if (is.null(cv) || !identical(requested.loss, "rmse")) {
        return(cv)
    }
    out <- cv
    out$fold.errors <- sqrt(out$fold.errors)
    out$mean.error <- apply(out$fold.errors, c(2L, 3L), mean, na.rm = TRUE)
    out$mean.error[!is.finite(out$mean.error)] <- Inf
    out$se <- apply(out$fold.errors, c(2L, 3L), function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    for (ii in seq_len(nrow(out$candidate.table))) {
        i <- out$candidate.table$lambda1.index[[ii]]
        j <- out$candidate.table$lambda2.index[[ii]]
        out$candidate.table$mean.error[[ii]] <- out$mean.error[i, j]
        out$candidate.table$se[[ii]] <- out$se[i, j]
    }
    out <- .slpl.tf.update.cv.selection(out)
    out$loss <- "rmse"
    out$candidate.table$loss <- "rmse"
    out$fold.fit.table$loss <- "rmse"
    out
}

.slpl.tf.update.cv.selection <- function(cv) {
    candidates <- cv$candidate.table
    candidates$selected <- FALSE
    best.idx <- which.min(candidates$mean.error)
    if (!length(best.idx) || !is.finite(candidates$mean.error[[best.idx]])) {
        stop("All S-LPL-TF CV candidates failed for the supplied grids.",
             call. = FALSE)
    }
    if (identical(cv$selection, "one.se")) {
        threshold <- candidates$mean.error[[best.idx]] + candidates$se[[best.idx]]
        eligible <- which(candidates$mean.error <= threshold)
        ord <- order(-candidates$lambda1[eligible],
                     -candidates$lambda2[eligible],
                     candidates$mean.error[eligible])
        selected.idx <- eligible[ord[[1L]]]
    } else {
        selected.idx <- best.idx
    }
    candidates$selected[[selected.idx]] <- TRUE
    cv$candidate.table <- candidates
    cv$best.idx <- best.idx
    cv$selected.idx <- selected.idx
    cv$lambda1.selected.idx <- candidates$lambda1.index[[selected.idx]]
    cv$lambda2.selected.idx <- candidates$lambda2.index[[selected.idx]]
    cv$lambda1.min <- candidates$lambda1[[best.idx]]
    cv$lambda2.min <- candidates$lambda2[[best.idx]]
    cv$lambda1 <- candidates$lambda1[[selected.idx]]
    cv$lambda2 <- candidates$lambda2[[selected.idx]]
    cv$error.min <- candidates$mean.error[[best.idx]]
    cv$error.selected <- candidates$mean.error[[selected.idx]]
    cv
}

.slpl.tf.validate.lambda.grid <- function(x, name) {
    .validate.ssrhe.lambda.grid(x, name)
}

.slpl.tf.fit.operator.grid <- function(X, y, operator.grid,
                                       base.operator.args, fit.args, call) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    y <- .lpl.tf.validate.response(y, nrow(X))
    candidates <- .lpl.tf.operator.grid.candidates(operator.grid)
    if (!length(candidates)) {
        stop("'operator.grid' must contain at least one candidate.",
             call. = FALSE)
    }
    fold.source <- if (identical(fit.args$lambda.selection, "fixed")) {
        NA_character_
    } else if (is.null(fit.args$foldid)) {
        if (is.null(fit.args$cv.seed)) {
            "generated_deterministic"
        } else {
            "generated_seeded"
        }
    } else {
        "supplied"
    }
    if (identical(fit.args$lambda.selection, "cv") &&
        is.null(fit.args$foldid)) {
        fit.args$foldid <- .lpl.tf.prepare.foldid(
            foldid = NULL,
            n = nrow(X),
            cv.folds = fit.args$cv.folds,
            cv.seed = fit.args$cv.seed,
            lambda.selection = fit.args$lambda.selection
        )
    }
    fit.list <- vector("list", length(candidates))
    table.list <- vector("list", length(candidates))
    for (ii in seq_along(candidates)) {
        candidate.args <- .lpl.tf.merge.operator.args(
            base.operator.args,
            candidates[[ii]]
        )
        op <- tryCatch(
            do.call(slpl.tf.operator, c(list(X = X), candidate.args)),
            error = function(e) e
        )
        if (inherits(op, "error")) {
            table.list[[ii]] <- .slpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "operator",
                message = conditionMessage(op)
            )
            next
        }
        if (nrow(op$A_LPL) < 1L) {
            table.list[[ii]] <- .slpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "operator",
                message = "operator_no_valid_lpl_rows",
                operator = op
            )
            next
        }
        candidate.fit.args <- fit.args
        if (identical(candidate.fit.args$lambda.selection, "fixed")) {
            candidate.fit.args$lambda1 <- candidate.fit.args$lambda1.grid[[1L]]
            candidate.fit.args$lambda2 <- candidate.fit.args$lambda2.grid[[1L]]
            candidate.fit.args$lambda1.grid <- NULL
            candidate.fit.args$lambda2.grid <- NULL
        }
        fit <- tryCatch(
            do.call(fit.slpl.tf, c(list(
                X = NULL,
                y = y,
                operator = op,
                operator.grid = NULL
            ), candidate.fit.args)),
            error = function(e) e
        )
        if (inherits(fit, "error")) {
            table.list[[ii]] <- .slpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "fit",
                message = conditionMessage(fit),
                operator = op
            )
            next
        }
        score <- .slpl.tf.operator.selection.score(fit)
        fit.list[[ii]] <- fit
        table.list[[ii]] <- .slpl.tf.operator.selection.row(
            candidate.id = ii,
            candidate.args = candidate.args,
            status = "ok",
            failure.stage = NA_character_,
            message = NA_character_,
            operator = op,
            fit = fit,
            score = score
        )
    }
    selection.table <- do.call(rbind, table.list)
    rownames(selection.table) <- NULL
    ok <- which(selection.table$status == "ok" &
                    is.finite(selection.table$selection.score))
    if (!length(ok)) {
        stop("All S-LPL-TF operator-grid candidates failed.", call. = FALSE)
    }
    selected.row <- ok[which.min(selection.table$selection.score[ok])]
    selection.table$selected <- seq_len(nrow(selection.table)) == selected.row
    best <- fit.list[[selected.row]]
    best$X <- X
    best$fold.source <- fold.source
    best$operator.selection <- list(
        selected.candidate.id = selected.row,
        score = selection.table$selection.score[selected.row],
        score.metric = if (identical(fit.args$lambda.selection, "cv")) {
            paste0("cv.", fit.args$cv.loss)
        } else {
            "objective"
        },
        table = selection.table,
        candidates = candidates,
        foldid = fit.args$foldid,
        fold.source = fold.source
    )
    best$call <- call
    best
}

.slpl.tf.operator.selection.score <- function(fit) {
    if (identical(fit$lambda.selection, "cv") && !is.null(fit$cv)) {
        return(fit$cv$error.selected %||%
                   fit$cv$candidate.table$mean.error[fit$cv$selected.idx])
    }
    fit$diagnostics$objective %||% Inf
}

.slpl.tf.operator.selection.row <- function(candidate.id, candidate.args,
                                            status, failure.stage, message,
                                            operator = NULL, fit = NULL,
                                            score = Inf) {
    get.arg <- function(name) {
        value <- candidate.args[[name]]
        if (is.null(value)) return(NA_character_)
        paste(value, collapse = ",")
    }
    get.setting <- function(name) {
        if (!is.null(operator) && !is.null(operator$settings[[name]])) {
            return(paste(operator$settings[[name]], collapse = ","))
        }
        get.arg(name)
    }
    data.frame(
        candidate.id = candidate.id,
        status = status,
        selected = FALSE,
        selection.score = score,
        lambda1 = if (is.null(fit)) NA_real_ else fit$lambda1,
        lambda2 = if (is.null(fit)) NA_real_ else fit$lambda2,
        lambda1.boundary = if (is.null(fit)) {
            NA_character_
        } else {
            fit$diagnostics$lambda1.boundary
        },
        lambda2.boundary = if (is.null(fit)) {
            NA_character_
        } else {
            fit$diagnostics$lambda2.boundary
        },
        sync.energy.l2 = if (is.null(fit)) {
            NA_real_
        } else {
            fit$diagnostics$sync.energy.l2
        },
        failure.stage = failure.stage,
        failure.message = message,
        degree = get.arg("degree"),
        support.type = get.arg("support.type"),
        support.size = get.arg("support.size"),
        min.support = get.arg("min.support"),
        support.buffer = get.arg("support.buffer"),
        kernel = get.arg("kernel"),
        coordinate.method = get.setting("coordinate.method"),
        chart.dim = get.setting("chart.dim"),
        support.metric = get.setting("support.metric"),
        n.lpl.rows = if (is.null(operator)) NA_integer_ else nrow(operator$A_LPL),
        n.sync.rows = if (is.null(operator)) NA_integer_ else nrow(operator$C_sync),
        stringsAsFactors = FALSE
    )
}

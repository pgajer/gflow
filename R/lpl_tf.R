#' Construct a Local Polynomial Lifting Trend-Filtering Operator
#'
#' Builds the Local Polynomial Lifting Trend Filtering (LPL-TF) analysis
#' operator from local polynomial prediction residuals. For each observed
#' target point \eqn{i}, the operator predicts
#' \eqn{f_i} from nearby values \eqn{f_j} using a local polynomial model and
#' stores the residual row
#' \deqn{r_i(f) = f_i - \sum_{j \in S_i} h_{ij} f_j.}
#'
#' @param X Numeric coordinate matrix with one row per observation.
#' @param adj.list,weight.list Optional supplied undirected graph adjacency and
#'   positive edge-length lists using 1-based vertex indices.
#' @param graph Optional graph object returned by \code{\link{create.rknn.graph}}
#'   or a list containing \code{adj.list}/\code{weight.list} or
#'   \code{adj_list}/\code{weight_list}.
#' @param graph.stage Graph stage to extract from \code{graph}: \code{"final"},
#'   \code{"raw"}, or \code{"pruned"}.
#' @param anchor.index Optional integer vector of observed anchor/target
#'   indices. The current implementation supports observed anchors only.
#' @param anchor.coordinates Reserved for later off-observation anchors; this
#'   implementation requires \code{NULL}.
#' @param degree Polynomial degree. The current implementation supports
#'   \code{0L}, \code{1L}, and \code{2L}.
#' @param support.type Coordinate support rule. \code{"knn"} uses the nearest
#'   \code{support.size} candidates including the target before self exclusion.
#'   \code{"adaptive.radius"} uses the nearest \code{min.support} positive
#'   predictors after self exclusion. \code{"fixed.radius"} uses all points
#'   within \code{radius}.
#' @param support.size Positive integer for \code{support.type = "knn"}.
#' @param radius Positive coordinate radius for
#'   \code{support.type = "fixed.radius"}.
#' @param min.support Minimum number of positive-weight predictors after self
#'   exclusion. If \code{NULL}, uses
#'   \code{choose(ncol(X) + degree, degree) + support.buffer}.
#' @param support.buffer Nonnegative integer added to the local polynomial
#'   design size when choosing default supports.
#' @param kernel Kernel used to weight local predictors by distance.
#' @param coordinate.method Coordinate chart used for the local polynomial
#'   design. \code{"coordinates"} uses centered ambient coordinates;
#'   \code{"local.pca"} uses a deterministic local PCA chart on the selected
#'   support.
#' @param chart.dim Chart dimension. For \code{"coordinates"}, this must be
#'   \code{NULL} or \code{ncol(X)}. For \code{"local.pca"}, \code{NULL}
#'   defaults to \code{ncol(X)}.
#' @param support.metric Support-distance rule. \code{"coordinates"} uses
#'   Euclidean coordinate distances; \code{"graph.geodesic"} uses shortest-path
#'   distances in the supplied graph. \code{"auto"} uses graph geodesics when a
#'   graph is supplied and coordinate distances otherwise.
#' @param exclude.self Logical. LPL-TF residual rows require \code{TRUE}.
#' @param row.normalize Complete-row normalization rule applied after setting
#'   target coefficient \code{+1} and predictor coefficients \code{-h_ij}.
#' @param local.solver Local linear algebra rule.
#' @param normal.equations.max.condition Maximum condition number for
#'   \code{local.solver = "auto"} to use normal equations.
#' @param duplicate.action Duplicate coordinate policy.
#' @param drop.rank.deficient Logical. The current implementation requires
#'   \code{TRUE}.
#' @param verbose Logical. Reserved for diagnostic messages.
#' @param ... Reserved.
#'
#' @details
#' The operator uses observed anchors, excludes the target point from every
#' predictor support, and drops requested-degree rank-deficient rows with
#' explicit metadata. There is no silent degree downgrade. Graph support
#' distances affect support selection and kernel bandwidths only; local
#' polynomial coordinates are controlled separately by \code{coordinate.method}.
#'
#' If \eqn{a_i} is the complete assembled residual row, then
#' \code{row.normalize = "l2"} uses \eqn{a_i / \|a_i\|_2}, and
#' \code{row.normalize = "l1"} uses \eqn{a_i / \|a_i\|_1}. Normalization is
#' applied to the whole row, not only to predictor coefficients.
#'
#' @return A list of class \code{"lpl_tf_operator"} containing the sparse
#'   operator matrix \code{A}, row metadata, supports, local design summaries,
#'   diagnostics, settings, and the matched call.
#' @export
lpl.tf.operator <- function(
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
    exclude.self = TRUE,
    row.normalize = c("l2", "none", "l1"),
    local.solver = c("auto", "normal.equations", "qr", "svd"),
    normal.equations.max.condition = 1e8,
    duplicate.action = c("keep", "error"),
    drop.rank.deficient = TRUE,
    verbose = FALSE,
    ...) {

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for lpl.tf.operator().",
             call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    support.type <- match.arg(support.type, c("adaptive.radius", "knn", "fixed.radius"))
    kernel <- match.arg(kernel, c("epanechnikov", "triangular", "gaussian", "tricube"))
    coordinate.method <- match.arg(coordinate.method, c("coordinates", "local.pca"))
    support.metric <- match.arg(support.metric, c("auto", "coordinates", "graph.geodesic"))
    row.normalize <- match.arg(row.normalize, c("l2", "none", "l1"))
    local.solver <- match.arg(local.solver, c("auto", "normal.equations", "qr", "svd"))
    duplicate.action <- match.arg(duplicate.action, c("keep", "error"))
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
        exclude.self = exclude.self,
        row.normalize = row.normalize,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        duplicate.action = duplicate.action,
        drop.rank.deficient = drop.rank.deficient,
        verbose = verbose
    )

    powers <- .lpl.tf.monomial.powers(args$chart.dim, args$degree)
    design.columns <- nrow(powers)
    min.predictors <- args$min.support
    coordinate.dist.matrix <- as.matrix(stats::dist(args$X))
    support.dist.matrix <- if (identical(args$support.metric, "graph.geodesic")) {
        args$graph.distance.matrix
    } else {
        coordinate.dist.matrix
    }

    supports <- vector("list", length(args$anchor.index))
    local.design <- vector("list", length(args$anchor.index))
    rows <- vector("list", length(args$anchor.index))
    row.table <- vector("list", length(args$anchor.index))
    prediction.weights <- vector("list", length(args$anchor.index))

    for (a in seq_along(args$anchor.index)) {
        target <- args$anchor.index[[a]]
        support <- .lpl.tf.coordinate.support(
            target = target,
            distances = support.dist.matrix[target, ],
            support.type = args$support.type,
            support.size = args$support.size,
            radius = args$radius,
            min.predictors = min.predictors,
            exclude.self = args$exclude.self
        )
        supports[[a]] <- support
        chart <- .lpl.tf.local.chart(
            X = args$X,
            target = target,
            candidate = support$candidate,
            coordinate.method = args$coordinate.method,
            chart.dim = args$chart.dim
        )
        X.chart <- matrix(0, nrow = nrow(args$X), ncol = args$chart.dim)
        X.chart[support$candidate, ] <- chart$coordinates
        row <- .lpl.tf.operator.row(
            X = X.chart,
            target = target,
            anchor.row = a,
            predictors = support$predictors,
            distances = support$distances,
            powers = powers,
            degree = args$degree,
            design.columns = design.columns,
            min.predictors = min.predictors,
            kernel = args$kernel,
            row.normalize = args$row.normalize,
            local.solver = args$local.solver,
            normal.equations.max.condition = args$normal.equations.max.condition,
            coordinate.method = args$coordinate.method,
            support.metric = args$support.metric
        )
        rows[[a]] <- row$row
        row.table[[a]] <- row$row.table
        prediction.weights[[a]] <- row$prediction.weights
        local.design[[a]] <- row$local.design
    }

    row.table <- do.call(rbind, row.table)
    ok <- row.table$status == "ok"
    A <- .lpl.tf.assemble.sparse(rows[ok], ncol = nrow(args$X))
    diagnostics <- .lpl.tf.operator.diagnostics(
        A = A,
        row.table = row.table,
        X = args$X,
        degree = args$degree,
        chart.dim = args$chart.dim
    )

    out <- list(
        A = A,
        row.table = row.table,
        supports = supports,
        prediction.weights = prediction.weights,
        local.design = local.design,
        diagnostics = diagnostics,
        settings = list(
            degree = args$degree,
            support.type = args$support.type,
            support.size = args$support.size,
            radius = args$radius,
            min.support = args$min.support,
            support.buffer = args$support.buffer,
            kernel = args$kernel,
            coordinate.method = args$coordinate.method,
            chart.dim = args$chart.dim,
            support.metric = args$support.metric,
            graph.stage = args$graph.stage,
            graph.summary = args$graph.summary,
            exclude.self = args$exclude.self,
            row.normalize = args$row.normalize,
            local.solver = args$local.solver,
            normal.equations.max.condition = args$normal.equations.max.condition,
            duplicate.action = args$duplicate.action,
            drop.rank.deficient = args$drop.rank.deficient
        ),
        call = match.call()
    )
    class(out) <- c("lpl_tf_operator", "list")
    out
}

#' @export
print.lpl_tf_operator <- function(x, ...) {
    cat("Local Polynomial Lifting Trend-Filtering operator\n")
    cat("  observations:", ncol(x$A), "\n")
    cat("  valid rows:", nrow(x$A), "\n")
    cat("  degree:", x$settings$degree, "\n")
    cat("  support type:", x$settings$support.type, "\n")
    cat("  row normalization:", x$settings$row.normalize, "\n")
    cat("  operator rank:", x$diagnostics$operator.rank, "\n")
    cat("  operator nullity:", x$diagnostics$operator.nullity, "\n")
    invisible(x)
}

#' Fit A Local Polynomial Lifting Trend Filter
#'
#' Fits the Phase-2 LPL-TF estimator with a fixed operator:
#' \deqn{\hat f = \arg\min_f \frac12\|y-f\|_2^2 +
#' \lambda\|A_{\mathrm{LPL}}f\|_1.}
#'
#' @param X Optional coordinate matrix used to build \code{\link{lpl.tf.operator}}
#'   when \code{operator} is \code{NULL}.
#' @param y Numeric response vector.
#' @param operator Optional \code{"lpl_tf_operator"} object.
#' @param adj.list,weight.list,graph Optional graph inputs passed to
#'   \code{\link{lpl.tf.operator}} when \code{operator} is \code{NULL}.
#' @param degree Polynomial degree passed to \code{\link{lpl.tf.operator}} when
#'   \code{operator} is \code{NULL}.
#' @param lambda.grid Optional nonnegative lambda grid.
#' @param lambda Optional fixed lambda shortcut. If supplied, it is used as the
#'   single fixed lambda and \code{lambda.selection} must be \code{"fixed"}.
#' @param lambda.selection \code{"cv"} or \code{"fixed"}. Phase 2 requires an
#'   explicit \code{lambda.grid} for CV so the candidate grid is not generated
#'   from the full response vector before fold scoring.
#' @param operator.grid Optional Phase-3 operator candidate grid. Supply a data
#'   frame with one row per candidate or a list of named lists. Candidate fields
#'   may include operator-construction arguments such as \code{degree},
#'   \code{support.type}, \code{support.size}, \code{min.support},
#'   \code{support.buffer}, \code{kernel}, \code{row.normalize}, and
#'   \code{local.solver}.
#' @param foldid Optional deterministic fold assignments for CV.
#' @param cv.folds Number of generated folds when \code{foldid = NULL}.
#' @param cv.loss Cross-validation loss, \code{"mse"}, \code{"rmse"}, or
#'   \code{"mae"}. \code{"rmse"} uses the same lambda ordering as MSE and reports
#'   square-rooted fold errors.
#' @param cv.seed Reserved for reproducible generated folds. Phase 2 generated
#'   folds are deterministic even when this is \code{NULL}.
#' @param cv.repeats Phase 2 supports only one CV repeat.
#' @param solver Phase 2 supports only \code{"genlasso"}.
#' @param selection Lambda selection rule, \code{"min"} or \code{"one.se"}.
#' @param n.lambda Number of generated lambdas when a genlasso path is used and
#'   \code{lambda.grid = NULL}.
#' @param maxsteps,minlam,approx,rtol,btol,eps Genlasso controls.
#' @param verbose Logical.
#' @param ... Additional arguments passed to \code{\link{lpl.tf.operator}} when
#'   the operator is built internally.
#'
#' @return A list of class \code{"lpl_tf"}.
#' @export
fit.lpl.tf <- function(
    X = NULL,
    y,
    operator = NULL,
    adj.list = NULL,
    weight.list = NULL,
    graph = NULL,
    degree = 2L,
    lambda.grid = NULL,
    lambda = NULL,
    lambda.selection = c("cv", "fixed"),
    operator.grid = NULL,
    foldid = NULL,
    cv.folds = 5L,
    cv.loss = c("mse", "rmse", "mae"),
    cv.seed = NULL,
    cv.repeats = 1L,
    solver = c("genlasso"),
    selection = c("min", "one.se"),
    n.lambda = 80L,
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    verbose = FALSE,
    ...) {

    if (!requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for fit.lpl.tf().", call. = FALSE)
    }
    lambda.selection <- match.arg(lambda.selection)
    cv.loss <- match.arg(cv.loss)
    solver <- match.arg(solver)
    selection <- match.arg(selection)
    n.lambda <- .validate.ssrhe.positive.integer(n.lambda, "n.lambda")
    cv.folds <- .validate.ssrhe.positive.integer(cv.folds, "cv.folds")
    cv.repeats <- .validate.ssrhe.positive.integer(cv.repeats, "cv.repeats")
    if (cv.repeats != 1L) {
        stop("Phase 2 fit.lpl.tf() supports cv.repeats = 1 only.",
             call. = FALSE)
    }
    if (!is.null(lambda)) {
        if (!is.null(lambda.grid)) {
            stop("Specify only one of 'lambda' or 'lambda.grid'.", call. = FALSE)
        }
        lambda.grid <- lambda
        lambda.selection <- "fixed"
    }
    if (identical(lambda.selection, "cv") && is.null(lambda.grid)) {
        stop("Phase 2 fit.lpl.tf() requires an explicit 'lambda.grid' ",
             "when lambda.selection = 'cv' to avoid response-dependent ",
             "full-data lambda-grid generation.", call. = FALSE)
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
        return(.lpl.tf.fit.operator.grid(
            X = X,
            y = y,
            operator.grid = operator.grid,
            base.operator.args = c(list(
                adj.list = adj.list,
                weight.list = weight.list,
                graph = graph,
                degree = degree
            ), operator.args),
            fit.args = list(
                lambda.grid = lambda.grid,
                lambda.selection = lambda.selection,
                foldid = foldid,
                cv.folds = cv.folds,
                cv.loss = cv.loss,
                cv.seed = cv.seed,
                cv.repeats = cv.repeats,
                solver = solver,
                selection = selection,
                n.lambda = n.lambda,
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
    operator <- do.call(.lpl.tf.prepare.fit.operator, c(list(
        X = X,
        operator = operator,
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        degree = degree
    ), operator.args))
    n <- ncol(operator$A)
    y <- .lpl.tf.validate.response(y, n)
    weights <- rep(1, n)
    lambda.grid <- .validate.ssrhe.hessian.l1.lambda.grid(lambda.grid,
                                                          lambda.selection)
    solver.operator <- .lpl.tf.prepare.solver.operator(operator)
    fold.source <- if (is.null(foldid)) {
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
        n.lambda = n.lambda,
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
    start <- proc.time()[["elapsed"]]
    raw.fit <- .lpl.tf.fit.from.operator(
        operator = solver.operator$operator,
        y = y,
        weights = weights,
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        n.lambda = n.lambda,
        fold.id = solver.args$fold.id,
        loss = loss.internal,
        selection = selection,
        solver.args = solver.args,
        verbose = verbose
    )
    runtime <- proc.time()[["elapsed"]] - start
    selected.idx <- raw.fit$selection$selected.index %||% 1L
    cv <- .lpl.tf.decorate.cv(raw.fit$cv, requested.loss = cv.loss)
    penalty.values <- as.vector(operator$A %*% raw.fit$fitted.values)
    active.tol <- 1e-8
    n.active <- sum(is.finite(penalty.values) & abs(penalty.values) > active.tol)
    diagnostics <- modifyList(raw.fit$diagnostics, list(
        fit.status = "ok",
        runtime.seconds = runtime,
        n.operator.rows = nrow(operator$A),
        n.active.rows = n.active,
        active.row.fraction = if (nrow(operator$A)) n.active / nrow(operator$A) else NA_real_,
        lambda.boundary = .lpl.tf.lambda.boundary(raw.fit$lambda.grid,
                                                  selected.idx),
        nonfinite.fitted.values = sum(!is.finite(raw.fit$fitted.values)),
        solver.converged = NA,
        solver.status = "convergence_not_reported_by_backend",
        fit.source = raw.fit$solver$backend,
        solver.duplicate.rows.collapsed =
            solver.operator$diagnostics$duplicate.rows.collapsed,
        solver.operator.rows = solver.operator$diagnostics$solver.rows,
        solver.operator.raw.rows = solver.operator$diagnostics$raw.rows,
        objective = raw.fit$objective,
        data.loss = raw.fit$energies$data.loss,
        l1.penalty = raw.fit$energies$hessian.l1
    ))
    out <- list(
        fitted.values = raw.fit$fitted.values,
        residuals = raw.fit$residuals,
        y = raw.fit$y,
        X = X,
        operator = operator,
        lambda = raw.fit$lambda,
        lambda.grid = raw.fit$lambda.grid,
        lambda.selection = lambda.selection,
        cv = cv,
        foldid = raw.fit$fold.id,
        fold.source = fold.source,
        cv.loss = cv.loss,
        cv.loss.internal = loss.internal,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        selection = selection,
        n.lambda = n.lambda,
        beta.grid = raw.fit$beta.grid,
        path = raw.fit$path,
        solver = raw.fit$solver,
        diagnostics = diagnostics,
        call = match.call()
    )
    class(out) <- c("lpl_tf", "list")
    out
}

#' Refit A Local Polynomial Lifting Trend Filter
#'
#' Reuses the fixed operator and selected lambda from an existing LPL-TF fit.
#'
#' @param object A \code{"lpl_tf"} object.
#' @param y New numeric response vector.
#' @param lambda Optional fixed lambda.
#' @param reuse.lambda Logical. If \code{TRUE} and \code{lambda = NULL}, reuse
#'   \code{object$lambda}.
#' @param verbose Logical.
#' @param ... Reserved.
#'
#' @return A refitted \code{"lpl_tf"} object.
#' @export
refit.lpl.tf <- function(object, y, lambda = NULL, reuse.lambda = TRUE,
                         verbose = FALSE, ...) {
    if (!inherits(object, "lpl_tf")) {
        stop("'object' must be an lpl_tf fit.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    if (is.null(lambda)) {
        if (!isTRUE(reuse.lambda)) {
            stop("Provide 'lambda' or set reuse.lambda = TRUE.", call. = FALSE)
        }
        lambda <- object$lambda
    }
    fit.lpl.tf(
        X = object$X,
        y = y,
        operator = object$operator,
        lambda = lambda,
        lambda.selection = "fixed",
        solver = object$solver$requested %||% "genlasso",
        verbose = verbose
    )
}

#' Predict From A Local Polynomial Lifting Trend Filter
#'
#' @param object A \code{"lpl_tf"} object.
#' @param newdata Must be \code{NULL} in Phase 2.
#' @param type Prediction type. Phase 2 supports only \code{"response"}.
#' @param allow.incomplete Reserved.
#' @param ... Reserved.
#'
#' @return Fitted values at the training points.
#' @method predict lpl_tf
#' @export
predict.lpl_tf <- function(object, newdata = NULL, type = c("response"),
                           allow.incomplete = FALSE, ...) {
    type <- match.arg(type)
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    if (!is.null(newdata)) {
        stop("Phase 2 predict() for lpl_tf supports training-point prediction only; ",
             "'newdata' must be NULL.", call. = FALSE)
    }
    if (!identical(type, "response")) {
        stop("Phase 2 predict() for lpl_tf supports type = 'response' only.",
             call. = FALSE)
    }
    object$fitted.values
}

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

.lpl.tf.prepare.solver.operator <- function(operator, digits = 14L) {
    A <- methods::as(operator$A, "dgCMatrix")
    raw.rows <- nrow(A)
    diagnostics <- list(
        raw.rows = raw.rows,
        solver.rows = raw.rows,
        duplicate.rows.collapsed = 0L,
        duplicate.groups = 0L,
        duplicate.key.rule = "signif_format",
        duplicate.key.digits = as.integer(digits)
    )
    if (raw.rows <= 1L) {
        operator$solver.row.map <- data.frame(
            solver.row.id = seq_len(raw.rows),
            source.rows = as.character(seq_len(raw.rows)),
            duplicate.count = rep(1L, raw.rows),
            stringsAsFactors = FALSE
        )
        return(list(operator = operator, diagnostics = diagnostics))
    }
    dense <- as.matrix(A)
    row.key <- apply(dense, 1L, function(z) {
        paste(format(signif(z, digits), digits = digits,
                     scientific = TRUE, trim = TRUE), collapse = "|")
    })
    unique.key <- unique(row.key)
    groups <- lapply(unique.key, function(key) which(row.key == key))
    duplicate.count <- vapply(groups, length, integer(1))
    if (all(duplicate.count == 1L)) {
        operator$solver.row.map <- data.frame(
            solver.row.id = seq_len(raw.rows),
            source.rows = as.character(seq_len(raw.rows)),
            duplicate.count = rep(1L, raw.rows),
            stringsAsFactors = FALSE
        )
        return(list(operator = operator, diagnostics = diagnostics))
    }
    collapsed <- do.call(rbind, lapply(groups, function(idx) {
        dense[idx[[1L]], ] * length(idx)
    }))
    solver.map <- data.frame(
        solver.row.id = seq_along(groups),
        source.rows = vapply(groups, paste, character(1), collapse = ","),
        duplicate.count = duplicate.count,
        stringsAsFactors = FALSE
    )
    operator$A <- methods::as(Matrix::Matrix(collapsed, sparse = TRUE), "dgCMatrix")
    operator$solver.row.map <- solver.map
    diagnostics$solver.rows <- nrow(operator$A)
    diagnostics$duplicate.rows.collapsed <- sum(duplicate.count - 1L)
    diagnostics$duplicate.groups <- sum(duplicate.count > 1L)
    list(operator = operator, diagnostics = diagnostics)
}

.lpl.tf.fit.from.operator <- function(operator,
                                      y,
                                      weights,
                                      lambda.grid,
                                      lambda.selection,
                                      n.lambda,
                                      fold.id,
                                      loss,
                                      selection,
                                      solver.args,
                                      verbose = FALSE) {
    if (is.null(operator$A)) {
        stop("operator must contain A.", call. = FALSE)
    }
    n <- ncol(operator$A)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("LPL-TF fitting currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    y.raw <- as.vector(y.info$Y)
    weights <- as.vector(W)
    y.clean <- y.raw
    y.clean[!is.finite(y.clean)] <- 0

    D.info <- .ssrhe.hessian.l1.penalty.matrix(
        operator$A,
        row.scaling = solver.args$row.scaling
    )
    diagnostics <- .ssrhe.hessian.l1.diagnostics(
        A = operator$A,
        D.info = D.info
    )
    if (is.null(lambda.grid)) {
        stop("lambda.grid is required for Phase 2 LPL-TF fitting.",
             call. = FALSE)
    }
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.",
             call. = FALSE)
    }

    cv <- NULL
    if (identical(lambda.selection, "cv")) {
        if (is.null(fold.id)) {
            fold.id <- .prepare.ssrhe.cv.folds(
                observed = as.vector(y.info$observed & weights > 0),
                nfolds = solver.args$nfolds,
                fold.id = NULL
            )
        }
        cv <- .lpl.tf.cv(
            y = y.clean,
            weights = weights,
            D = D.info$D,
            lambda.grid = lambda.grid,
            fold.id = fold.id,
            loss = loss,
            selection = selection,
            solver.args = solver.args,
            use.svd = D.info$svd
        )
        selected.idx <- cv$selected.idx
    } else {
        selected.idx <- 1L
    }

    path.info <- .lpl.tf.genlasso.path(
        y = y.clean,
        weights = weights,
        D = D.info$D,
        solver.args = solver.args,
        use.svd = D.info$svd
    )
    beta.grid <- .ssrhe.hessian.l1.coef.matrix(
        path.info$path,
        lambda.grid,
        n
    )
    if (any(!is.finite(beta.grid))) {
        stop("LPL-TF genlasso path produced non-finite fitted values.",
             call. = FALSE)
    }
    fitted <- beta.grid[, selected.idx]
    observed <- as.vector(y.info$observed & weights > 0)
    residuals <- rep(NA_real_, n)
    residuals[observed] <- y.raw[observed] - fitted[observed]
    data.loss <- 0.5 * sum(weights * (y.clean - fitted)^2)
    hessian.l1 <- sum(abs(as.vector(operator$A %*% fitted)))
    lambda <- lambda.grid[selected.idx]

    list(
        fitted.values = as.vector(fitted),
        residuals = residuals,
        y = y.raw,
        weights = weights,
        lambda = lambda,
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        objective = data.loss + lambda * hessian.l1,
        energies = list(data.loss = data.loss, hessian.l1 = hessian.l1),
        operator = operator,
        beta.grid = beta.grid,
        path = path.info$path,
        cv = cv,
        fold.id = fold.id,
        observed = observed,
        n.responses = 1L,
        diagnostics = diagnostics,
        solver = list(
            backend = "genlasso_design_matrix",
            requested = solver.args$solver,
            representation = D.info$representation,
            svd = D.info$svd,
            row.scaling = D.info$row.scaling,
            n.observed = path.info$n.observed,
            admm = NULL
        ),
        selection = list(
            rule = selection,
            loss = loss,
            selected.index = selected.idx,
            lambda = lambda
        )
    )
}

.lpl.tf.genlasso.path <- function(y, weights, D, solver.args,
                                  use.svd = FALSE) {
    observed <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(observed)) {
        stop("At least one observed positive-weight response is required.",
             call. = FALSE)
    }
    obs <- which(observed)
    sw <- sqrt(weights[obs])
    X.design <- as.matrix(Matrix::sparseMatrix(
        i = seq_along(obs),
        j = obs,
        x = sw,
        dims = c(length(obs), length(y)),
        giveCsparse = TRUE
    ))
    path <- suppressWarnings(genlasso::genlasso(
        y = sw * y[obs],
        X = X.design,
        D = D,
        approx = solver.args$approx,
        maxsteps = solver.args$maxsteps,
        minlam = solver.args$minlam,
        rtol = solver.args$rtol,
        btol = solver.args$btol,
        eps = solver.args$eps,
        verbose = solver.args$verbose,
        svd = use.svd
    ))
    list(path = path, n.observed = sum(observed))
}

.lpl.tf.cv <- function(y, weights, D, lambda.grid, fold.id, loss,
                       selection, solver.args, use.svd = FALSE) {
    n <- length(y)
    folds <- sort(unique(fold.id[fold.id > 0L]))
    cv.errors <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda.grid))
    rownames(cv.errors) <- paste0("Fold", folds)
    colnames(cv.errors) <- format(signif(lambda.grid, 6), scientific = TRUE)

    for (ii in seq_along(folds)) {
        fold <- folds[ii]
        test <- which(fold.id == fold)
        train.weights <- weights
        train.weights[test] <- 0
        path.info <- tryCatch(
            .lpl.tf.genlasso.path(
                y = y,
                weights = train.weights,
                D = D,
                solver.args = solver.args,
                use.svd = use.svd
            ),
            error = function(e) NULL
        )
        if (is.null(path.info)) {
            next
        }
        beta <- .ssrhe.hessian.l1.coef.matrix(path.info$path, lambda.grid, n)
        beta.test <- beta[test, , drop = FALSE]
        finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
        target <- matrix(y[test], nrow = length(test), ncol = length(lambda.grid))
        w.test <- matrix(weights[test], nrow = length(test), ncol = length(lambda.grid))
        if (identical(loss, "mse")) {
            fold.error <- colSums(w.test * (target - beta.test)^2) / colSums(w.test)
        } else {
            fold.error <- colSums(w.test * abs(target - beta.test)) / colSums(w.test)
        }
        fold.error[!finite.pred] <- NA_real_
        cv.errors[ii, ] <- fold.error
    }

    mean.error <- colMeans(cv.errors, na.rm = TRUE)
    mean.error[!is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    best.idx <- which.min(mean.error)
    if (!is.finite(mean.error[best.idx])) {
        stop("All LPL-TF CV fits failed for the supplied lambda grid.",
             call. = FALSE)
    }
    if (identical(selection, "one.se")) {
        threshold <- mean.error[best.idx] + se[best.idx]
        eligible <- which(mean.error <= threshold)
        selected.idx <- eligible[which.max(lambda.grid[eligible])]
    } else {
        selected.idx <- best.idx
    }
    list(
        fold.id = fold.id,
        lambda.grid = lambda.grid,
        fold.errors = cv.errors,
        mean.error = mean.error,
        se = se,
        best.idx = best.idx,
        selected.idx = selected.idx,
        lambda.min = lambda.grid[best.idx],
        lambda = lambda.grid[selected.idx],
        error.min = mean.error[best.idx],
        error.selected = mean.error[selected.idx],
        selection = selection,
        loss = loss
    )
}

.lpl.tf.fit.operator.grid <- function(X, y, operator.grid,
                                      base.operator.args, fit.args, call) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    y <- .lpl.tf.validate.response(y, nrow(X))
    candidates <- .lpl.tf.operator.grid.candidates(operator.grid)
    if (!length(candidates)) {
        stop("'operator.grid' must contain at least one candidate.",
             call. = FALSE)
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
            do.call(lpl.tf.operator, c(list(X = X), candidate.args)),
            error = function(e) e
        )
        if (inherits(op, "error")) {
            table.list[[ii]] <- .lpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "operator",
                message = conditionMessage(op)
            )
            next
        }
        if (nrow(op$A) < 1L) {
            table.list[[ii]] <- .lpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "operator",
                message = "operator_no_valid_rows",
                operator = op
            )
            next
        }
        fit <- tryCatch(
            do.call(fit.lpl.tf, c(list(
                X = NULL,
                y = y,
                operator = op,
                operator.grid = NULL
            ), fit.args)),
            error = function(e) e
        )
        if (inherits(fit, "error")) {
            table.list[[ii]] <- .lpl.tf.operator.selection.row(
                candidate.id = ii,
                candidate.args = candidate.args,
                status = "failed",
                failure.stage = "fit",
                message = conditionMessage(fit),
                operator = op
            )
            next
        }
        score <- .lpl.tf.operator.selection.score(fit)
        fit.list[[ii]] <- fit
        table.list[[ii]] <- .lpl.tf.operator.selection.row(
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
        stop("All LPL-TF operator-grid candidates failed.", call. = FALSE)
    }
    selected.row <- ok[which.min(selection.table$selection.score[ok])]
    selection.table$selected <- seq_len(nrow(selection.table)) == selected.row
    best <- fit.list[[selected.row]]
    best$X <- X
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
        foldid = fit.args$foldid
    )
    best$call <- call
    best
}

.lpl.tf.operator.grid.candidates <- function(operator.grid) {
    if (is.data.frame(operator.grid)) {
        out <- vector("list", nrow(operator.grid))
        for (ii in seq_len(nrow(operator.grid))) {
            row <- as.list(operator.grid[ii, , drop = FALSE])
            row <- row[!vapply(row, function(x) length(x) == 1L && is.na(x),
                               logical(1))]
            out[[ii]] <- row
        }
        return(out)
    }
    if (is.list(operator.grid)) {
        if (!length(operator.grid)) return(list())
        if (all(vapply(operator.grid, is.list, logical(1)))) {
            return(operator.grid)
        }
        return(list(operator.grid))
    }
    stop("'operator.grid' must be a data frame or a list.", call. = FALSE)
}

.lpl.tf.merge.operator.args <- function(base.args, candidate.args) {
    out <- base.args
    for (nm in names(candidate.args)) {
        out[[nm]] <- candidate.args[[nm]]
    }
    out
}

.lpl.tf.operator.selection.score <- function(fit) {
    if (identical(fit$lambda.selection, "cv") && !is.null(fit$cv)) {
        return(fit$cv$error.selected %||% fit$cv$mean.error[fit$cv$selected.idx])
    }
    fit$diagnostics$objective %||% Inf
}

.lpl.tf.operator.selection.row <- function(candidate.id, candidate.args,
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
        lambda = if (is.null(fit)) NA_real_ else fit$lambda,
        lambda.boundary = if (is.null(fit)) {
            NA_character_
        } else {
            fit$diagnostics$lambda.boundary
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
        graph.stage = get.setting("graph.stage"),
        row.normalize = get.arg("row.normalize"),
        local.solver = get.arg("local.solver"),
        n.operator.rows = if (is.null(operator)) NA_integer_ else nrow(operator$A),
        operator.rank = if (is.null(operator)) {
            NA_integer_
        } else {
            operator$diagnostics$operator.rank
        },
        operator.nullity = if (is.null(operator)) {
            NA_integer_
        } else {
            operator$diagnostics$operator.nullity
        },
        polynomial.residual = if (is.null(operator)) {
            NA_real_
        } else {
            operator$diagnostics$polynomial.residual
        },
        stringsAsFactors = FALSE
    )
}

.lpl.tf.prepare.fit.operator <- function(X, operator, adj.list, weight.list,
                                         graph, degree, ...) {
    if (!is.null(operator)) {
        if (!inherits(operator, "lpl_tf_operator")) {
            stop("'operator' must be an lpl_tf_operator object.", call. = FALSE)
        }
        if (is.null(operator$A)) {
            stop("'operator' must contain an operator matrix 'A'.", call. = FALSE)
        }
        return(operator)
    }
    if (is.null(X)) {
        stop("Provide 'operator' or 'X' to fit.lpl.tf().", call. = FALSE)
    }
    lpl.tf.operator(
        X = X,
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        degree = degree,
        ...
    )
}

.lpl.tf.validate.response <- function(y, n) {
    if (is.matrix(y) || is.data.frame(y)) {
        if (NCOL(y) != 1L) {
            stop("fit.lpl.tf() currently supports one response vector.",
                 call. = FALSE)
        }
        y <- as.vector(as.matrix(y))
    }
    if (!is.numeric(y) || length(y) != n || anyNA(y) ||
        any(!is.finite(y))) {
        stop("'y' must be a finite numeric vector of length ncol(operator$A).",
             call. = FALSE)
    }
    as.numeric(y)
}

.lpl.tf.prepare.foldid <- function(foldid, n, cv.folds, cv.seed,
                                   lambda.selection) {
    if (identical(lambda.selection, "fixed")) {
        return(NULL)
    }
    observed <- rep(TRUE, n)
    if (!is.null(foldid)) {
        return(.prepare.ssrhe.cv.folds(observed = observed,
                                       nfolds = cv.folds,
                                       fold.id = foldid))
    }
    cv.folds <- min(cv.folds, n)
    if (cv.folds < 2L) {
        stop("At least two observations are required for CV.", call. = FALSE)
    }
    if (!is.null(cv.seed)) {
        old.seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                               inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old.seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old.seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(cv.seed)
        foldid <- sample(rep(seq_len(cv.folds), length.out = n))
    } else {
        foldid <- rep(seq_len(cv.folds), length.out = n)
    }
    .prepare.ssrhe.cv.folds(observed = observed, nfolds = cv.folds,
                            fold.id = foldid)
}

.lpl.tf.decorate.cv <- function(cv, requested.loss) {
    if (is.null(cv) || !identical(requested.loss, "rmse")) {
        return(cv)
    }
    out <- cv
    out$fold.errors <- sqrt(out$fold.errors)
    out$mean.error <- colMeans(out$fold.errors, na.rm = TRUE)
    out$mean.error[!is.finite(out$mean.error)] <- Inf
    out$se <- apply(out$fold.errors, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    out$loss <- "rmse"
    out
}

.lpl.tf.lambda.boundary <- function(lambda.grid, selected.idx) {
    if (!length(lambda.grid) || is.na(selected.idx)) return(NA_character_)
    if (length(lambda.grid) == 1L) return("none")
    selected <- lambda.grid[selected.idx]
    if (!is.finite(selected)) return(NA_character_)
    if (isTRUE(all.equal(selected, min(lambda.grid, na.rm = TRUE)))) {
        return("lower")
    }
    if (isTRUE(all.equal(selected, max(lambda.grid, na.rm = TRUE)))) {
        return("upper")
    }
    "none"
}

.lpl.tf.validate.operator.args <- function(
    X, adj.list, weight.list, graph, graph.stage, anchor.index,
    anchor.coordinates, degree, support.type, support.size, radius,
    min.support, support.buffer, kernel, coordinate.method, chart.dim,
    support.metric, exclude.self, row.normalize, local.solver,
    normal.equations.max.condition, duplicate.action, drop.rank.deficient,
    verbose) {

    X <- as.matrix(X)
    storage.mode(X) <- "double"
    if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) < 1L ||
        ncol(X) < 1L || any(!is.finite(X))) {
        stop("'X' must be a finite numeric matrix with at least one row and column.",
             call. = FALSE)
    }
    duplicate.action <- match.arg(duplicate.action, c("keep", "error"))
    if (identical(duplicate.action, "error") && any(duplicated(as.data.frame(X)))) {
        stop("Duplicate coordinate rows were found in 'X'.", call. = FALSE)
    }
    if (!is.null(anchor.coordinates)) {
        stop("lpl.tf.operator() supports observed anchors only; ",
             "'anchor.coordinates' must be NULL.", call. = FALSE)
    }
    n <- nrow(X)
    graph.stage <- match.arg(graph.stage, c("final", "raw", "pruned"))
    graph.info <- .lpl.tf.prepare.graph(
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        graph.stage = graph.stage,
        n = n
    )
    if (is.null(anchor.index)) {
        anchor.index <- seq_len(n)
    } else {
        if (!is.numeric(anchor.index) || any(!is.finite(anchor.index)) ||
            any(anchor.index != as.integer(anchor.index))) {
            stop("'anchor.index' must be an integer vector of observed row indices.",
                 call. = FALSE)
        }
        anchor.index <- as.integer(anchor.index)
        if (!length(anchor.index) || any(anchor.index < 1L | anchor.index > n)) {
            stop("'anchor.index' must contain observed row indices in 'X'.",
                 call. = FALSE)
        }
    }
    degree <- .lpl.tf.scalar.integer(degree, "degree", min = 0L, max = 2L)
    support.type <- match.arg(support.type, c("adaptive.radius", "knn", "fixed.radius"))
    support.buffer <- .lpl.tf.scalar.integer(
        support.buffer, "support.buffer", min = 0L
    )
    coordinate.method <- match.arg(coordinate.method, c("coordinates", "local.pca"))
    if (is.null(chart.dim)) {
        chart.dim <- ncol(X)
    } else {
        chart.dim <- .lpl.tf.scalar.integer(chart.dim, "chart.dim", min = 1L,
                                            max = ncol(X))
    }
    if (identical(coordinate.method, "coordinates") &&
        chart.dim != ncol(X)) {
        stop("'chart.dim' must be NULL or ncol(X) when coordinate.method = 'coordinates'.",
             call. = FALSE)
    }
    design.columns <- choose(chart.dim + degree, degree)
    if (is.null(min.support)) {
        min.support <- as.integer(design.columns + support.buffer)
    } else {
        min.support <- .lpl.tf.scalar.integer(min.support, "min.support", min = 1L)
    }
    if (is.null(support.size)) {
        support.size <- if (identical(support.type, "knn")) {
            as.integer(min.support + 1L)
        } else {
            NA_integer_
        }
    } else {
        support.size <- .lpl.tf.scalar.integer(support.size, "support.size", min = 1L)
    }
    if (identical(support.type, "fixed.radius")) {
        if (is.null(radius) || !is.numeric(radius) || length(radius) != 1L ||
            !is.finite(radius) || radius <= 0) {
            stop("'radius' must be a positive finite scalar for fixed-radius supports.",
                 call. = FALSE)
        }
        radius <- as.double(radius)
    } else if (!is.null(radius)) {
        if (!is.numeric(radius) || length(radius) != 1L ||
            !is.finite(radius) || radius <= 0) {
            stop("'radius' must be NULL or a positive finite scalar.",
                 call. = FALSE)
        }
        radius <- as.double(radius)
    }
    kernel <- match.arg(kernel, c("epanechnikov", "triangular", "gaussian", "tricube"))
    support.metric <- match.arg(support.metric, c("auto", "coordinates", "graph.geodesic"))
    if (identical(support.metric, "auto")) {
        support.metric <- if (is.null(graph.info)) "coordinates" else "graph.geodesic"
    }
    graph.distance.matrix <- NULL
    graph.summary <- list(constructor = "none", n.vertices = n,
                          n.edges = NA_integer_)
    if (identical(support.metric, "graph.geodesic")) {
        if (is.null(graph.info)) {
            stop("support.metric = 'graph.geodesic' requires 'graph' or both ",
                 "'adj.list' and 'weight.list'.", call. = FALSE)
        }
        graph.distance.matrix <- .pttf.geometry.all.source.distances(
            graph.info$adj.list,
            graph.info$weight.list
        )
        if (any(!is.finite(graph.distance.matrix))) {
            stop("Graph geodesic support requires a connected graph.",
                 call. = FALSE)
        }
        graph.summary <- graph.info$summary
    }
    if (!isTRUE(exclude.self)) {
        stop("lpl.tf.operator() requires exclude.self = TRUE.",
             call. = FALSE)
    }
    row.normalize <- match.arg(row.normalize, c("l2", "none", "l1"))
    local.solver <- match.arg(local.solver, c("auto", "normal.equations", "qr", "svd"))
    if (!is.numeric(normal.equations.max.condition) ||
        length(normal.equations.max.condition) != 1L ||
        !is.finite(normal.equations.max.condition) ||
        normal.equations.max.condition <= 0) {
        stop("'normal.equations.max.condition' must be a positive finite scalar.",
             call. = FALSE)
    }
    if (!isTRUE(drop.rank.deficient)) {
        stop("lpl.tf.operator() requires drop.rank.deficient = TRUE.",
             call. = FALSE)
    }
    verbose <- isTRUE(verbose)
    list(
        X = X,
        anchor.index = anchor.index,
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
        graph.stage = graph.stage,
        graph.distance.matrix = graph.distance.matrix,
        graph.summary = graph.summary,
        exclude.self = TRUE,
        row.normalize = row.normalize,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        duplicate.action = duplicate.action,
        drop.rank.deficient = TRUE,
        verbose = verbose
    )
}

.lpl.tf.prepare.graph <- function(adj.list, weight.list, graph, graph.stage, n) {
    if (!is.null(graph) && (!is.null(adj.list) || !is.null(weight.list))) {
        stop("Specify either 'graph' or 'adj.list'/'weight.list', not both.",
             call. = FALSE)
    }
    if (is.null(graph) && is.null(adj.list) && is.null(weight.list)) {
        return(NULL)
    }
    if (is.null(graph)) {
        if (is.null(adj.list) || is.null(weight.list)) {
            stop("Graph-geodesic support requires both 'adj.list' and 'weight.list'.",
                 call. = FALSE)
        }
        source <- "supplied"
    } else {
        extracted <- .lpl.tf.extract.graph.stage(graph, graph.stage)
        adj.list <- extracted$adj.list
        weight.list <- extracted$weight.list
        source <- extracted$source
    }
    validated <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    if (length(validated$adj.list) != n) {
        stop("The graph vertex count must match nrow(X).", call. = FALSE)
    }
    edge.table <- .pttf.geometry.edge.table(validated$adj.list,
                                            validated$weight.list)
    list(
        adj.list = validated$adj.list,
        weight.list = validated$weight.list,
        summary = list(
            constructor = source,
            stage = graph.stage,
            n.vertices = length(validated$adj.list),
            n.edges = nrow(edge.table),
            edge.length.summary = summary(edge.table$length)
        )
    )
}

.lpl.tf.extract.graph.stage <- function(graph, graph.stage) {
    if (!is.list(graph)) {
        stop("'graph' must be a list-like graph object.", call. = FALSE)
    }
    keys <- switch(
        graph.stage,
        final = list(c("adj_list", "weight_list"),
                     c("adj.list", "weight.list")),
        raw = list(c("raw_adj_list", "raw_weight_list"),
                   c("raw.adj.list", "raw.weight.list")),
        pruned = list(c("pruned_adj_list", "pruned_weight_list"),
                      c("pruned.adj.list", "pruned.weight.list"))
    )
    for (pair in keys) {
        if (!is.null(graph[[pair[[1L]]]]) && !is.null(graph[[pair[[2L]]]])) {
            return(list(
                adj.list = graph[[pair[[1L]]]],
                weight.list = graph[[pair[[2L]]]],
                source = paste0("graph.", graph.stage)
            ))
        }
    }
    stop("Could not extract graph adjacency and weights for graph.stage = '",
         graph.stage, "'.", call. = FALSE)
}

.lpl.tf.scalar.integer <- function(x, name, min = -Inf, max = Inf) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != as.integer(x) || x < min || x > max) {
        stop("'", name, "' must be an integer scalar", call. = FALSE)
    }
    as.integer(x)
}

.lpl.tf.monomial.powers <- function(dim, degree) {
    out <- list()
    rec <- function(prefix, remaining.dim, remaining.total) {
        if (remaining.dim == 1L) {
            out[[length(out) + 1L]] <<- c(prefix, remaining.total)
            return(invisible(NULL))
        }
        for (e in 0:remaining.total) {
            rec(c(prefix, e), remaining.dim - 1L, remaining.total - e)
        }
    }
    for (total in 0:degree) {
        start <- length(out) + 1L
        rec(integer(0), dim, total)
        out[start:length(out)] <- rev(out[start:length(out)])
    }
    mat <- do.call(rbind, out)
    storage.mode(mat) <- "integer"
    mat
}

.lpl.tf.design.matrix <- function(Z, powers) {
    Z <- as.matrix(Z)
    D <- matrix(1, nrow(Z), nrow(powers))
    for (k in seq_len(nrow(powers))) {
        for (j in seq_len(ncol(powers))) {
            if (powers[k, j] != 0L) {
                D[, k] <- D[, k] * Z[, j]^powers[k, j]
            }
        }
    }
    D
}

.lpl.tf.coordinate.support <- function(target, distances, support.type,
                                       support.size, radius, min.predictors,
                                       exclude.self) {
    n <- length(distances)
    ord <- order(distances, seq_len(n), na.last = NA)
    if (isTRUE(exclude.self)) {
        ord.no.self <- ord[ord != target]
    } else {
        ord.no.self <- ord
    }
    if (identical(support.type, "knn")) {
        predictors <- head(ord.no.self, support.size - 1L)
        candidate <- c(target, predictors)
    } else if (identical(support.type, "fixed.radius")) {
        candidate <- which(distances <= radius)
        predictors <- setdiff(candidate, target)
        predictors <- predictors[order(distances[predictors], predictors)]
    } else {
        predictors <- head(ord.no.self, min.predictors)
        candidate <- c(target, predictors)
    }
    list(
        target = target,
        candidate = as.integer(candidate),
        predictors = as.integer(predictors),
        distances = as.double(distances[predictors]),
        radius = if (length(predictors)) max(distances[predictors]) else NA_real_
    )
}

.lpl.tf.local.chart <- function(X, target, candidate, coordinate.method,
                                chart.dim) {
    centered <- sweep(X[candidate, , drop = FALSE], 2L,
                      X[target, , drop = TRUE], "-")
    if (identical(coordinate.method, "coordinates")) {
        return(list(
            coordinates = centered,
            basis = diag(ncol(X)),
            singular.values = NA_real_
        ))
    }
    sv <- svd(centered, nu = 0L, nv = chart.dim)
    basis <- sv$v[, seq_len(chart.dim), drop = FALSE]
    list(
        coordinates = centered %*% basis,
        basis = basis,
        singular.values = sv$d
    )
}

.lpl.tf.kernel.weights <- function(distances, kernel) {
    if (!length(distances)) return(numeric(0))
    bandwidth <- max(distances[is.finite(distances)], 0)
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

.lpl.tf.operator.row <- function(X, target, anchor.row, predictors, distances,
                                 powers, degree, design.columns,
                                 min.predictors, kernel, row.normalize,
                                 local.solver,
                                 normal.equations.max.condition,
                                 coordinate.method, support.metric) {
    predictor.size <- length(predictors)
    weights <- .lpl.tf.kernel.weights(distances, kernel)
    positive <- is.finite(weights) & weights > 0
    predictors.positive <- predictors[positive]
    weights.positive <- weights[positive]
    row.base <- data.frame(
        row.id = anchor.row,
        target.index = target,
        anchor.index = target,
        status = "dropped",
        drop.reason = NA_character_,
        support.size = predictor.size + 1L,
        positive.weight.support.size = length(predictors.positive) + 1L,
        predictor.size = predictor.size,
        positive.weight.predictor.size = length(predictors.positive),
        degree = degree,
        design.columns = design.columns,
        rank = NA_integer_,
        condition.number = NA_real_,
        solver.used = NA_character_,
        fallback.used = NA,
        row.norm.raw = NA_real_,
        row.norm.raw.l2 = NA_real_,
        row.norm.raw.l1 = NA_real_,
        row.norm.used = NA_character_,
        row.norm.final = NA_real_,
        radius = if (length(distances)) max(distances) else NA_real_,
        kernel = kernel,
        coordinate.method = coordinate.method,
        support.metric = support.metric,
        stringsAsFactors = FALSE
    )
    if (length(predictors.positive) < min.predictors) {
        row.base$drop.reason <- "too_few_positive_predictors_after_self_exclusion"
        return(.lpl.tf.row.result(NULL, row.base, NULL, numeric(0), predictors.positive))
    }
    Z <- sweep(X[predictors.positive, , drop = FALSE], 2L,
               X[target, , drop = TRUE], "-")
    D <- .lpl.tf.design.matrix(Z, powers)
    target.basis <- .lpl.tf.design.matrix(
        matrix(0, nrow = 1L, ncol = ncol(X)), powers
    )
    fit <- .lpl.tf.prediction.weights(
        D = D,
        target.basis = target.basis,
        weights = weights.positive,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition
    )
    row.base$rank <- fit$rank
    row.base$condition.number <- fit$condition.number
    row.base$solver.used <- fit$solver.used
    row.base$fallback.used <- fit$fallback.used
    if (fit$rank < design.columns) {
        row.base$drop.reason <- "rank_deficient_requested_degree"
        return(.lpl.tf.row.result(NULL, row.base, D, numeric(0), predictors.positive))
    }
    if (!fit$ok) {
        row.base$drop.reason <- fit$drop.reason
        return(.lpl.tf.row.result(NULL, row.base, D, numeric(0), predictors.positive))
    }
    full.row <- numeric(nrow(X))
    full.row[target] <- 1
    full.row[predictors.positive] <- full.row[predictors.positive] - fit$h
    normed <- .lpl.tf.normalize.row(full.row, row.normalize)
    row.base$row.norm.raw <- normed$row.norm.raw
    row.base$row.norm.raw.l2 <- normed$row.norm.raw.l2
    row.base$row.norm.raw.l1 <- normed$row.norm.raw.l1
    row.base$row.norm.used <- normed$row.norm.used
    row.base$row.norm.final <- normed$row.norm.final
    if (!normed$ok) {
        row.base$drop.reason <- normed$drop.reason
        return(.lpl.tf.row.result(NULL, row.base, D, fit$h, predictors.positive))
    }
    row.base$status <- "ok"
    row.base$drop.reason <- NA_character_
    .lpl.tf.row.result(normed$row, row.base, D, fit$h, predictors.positive)
}

.lpl.tf.prediction.weights <- function(D, target.basis, weights, local.solver,
                                       normal.equations.max.condition) {
    sqrt.w <- sqrt(weights)
    Xw <- D * sqrt.w
    sv <- svd(Xw, nu = 0L, nv = 0L)$d
    tol <- max(dim(Xw)) * max(sv, 1) * .Machine$double.eps
    positive <- sv[sv > tol]
    rank <- length(positive)
    p <- ncol(D)
    condition.number <- if (rank < p || length(positive) < 1L) {
        Inf
    } else {
        max(positive) / min(positive)
    }
    solver.used <- local.solver
    fallback.used <- FALSE
    if (identical(local.solver, "auto")) {
        if (rank == p && is.finite(condition.number) &&
            condition.number <= normal.equations.max.condition) {
            solver.used <- "normal.equations"
        } else {
            solver.used <- "svd"
            fallback.used <- TRUE
        }
    }
    if (rank < p) {
        return(list(ok = FALSE, h = numeric(0), rank = rank,
                    condition.number = condition.number,
                    solver.used = solver.used,
                    fallback.used = fallback.used,
                    drop.reason = "rank_deficient_requested_degree"))
    }
    h <- tryCatch({
        if (identical(solver.used, "normal.equations")) {
            XtWX <- crossprod(D, weights * D)
            as.numeric(target.basis %*% solve(XtWX, t(D) * rep(weights, each = p)))
        } else if (identical(solver.used, "qr")) {
            coef.map <- qr.solve(Xw, diag(nrow(D)))
            as.numeric(target.basis %*% coef.map * sqrt.w)
        } else {
            sv.full <- svd(Xw)
            tol.full <- max(dim(Xw)) * max(sv.full$d, 1) * .Machine$double.eps
            inv.d <- ifelse(sv.full$d > tol.full, 1 / sv.full$d, 0)
            coef.map <- sv.full$v %*% (inv.d * crossprod(sv.full$u, diag(nrow(D))))
            as.numeric(target.basis %*% coef.map * sqrt.w)
        }
    }, error = function(e) NULL)
    if (is.null(h) || length(h) != nrow(D) || any(!is.finite(h))) {
        return(list(ok = FALSE, h = numeric(0), rank = rank,
                    condition.number = condition.number,
                    solver.used = solver.used,
                    fallback.used = fallback.used,
                    drop.reason = "linear_solve_failed"))
    }
    list(ok = TRUE, h = h, rank = rank,
         condition.number = condition.number,
         solver.used = solver.used, fallback.used = fallback.used,
         drop.reason = NA_character_)
}

.lpl.tf.normalize.row <- function(row, row.normalize, min.norm = 1e-12) {
    raw.l2 <- sqrt(sum(row^2))
    raw.l1 <- sum(abs(row))
    if (!is.finite(raw.l2) || raw.l2 < min.norm) {
        return(list(ok = FALSE, row = row, row.norm.raw = raw.l2,
                    row.norm.raw.l2 = raw.l2, row.norm.raw.l1 = raw.l1,
                    row.norm.used = row.normalize,
                    row.norm.final = NA_real_,
                    drop.reason = "near_zero_row_norm"))
    }
    final <- switch(
        row.normalize,
        none = row,
        l2 = row / raw.l2,
        l1 = {
            if (!is.finite(raw.l1) || raw.l1 < min.norm) {
                return(list(ok = FALSE, row = row, row.norm.raw = raw.l1,
                            row.norm.raw.l2 = raw.l2, row.norm.raw.l1 = raw.l1,
                            row.norm.used = row.normalize,
                            row.norm.final = NA_real_,
                            drop.reason = "near_zero_row_norm"))
            }
            row / raw.l1
        }
    )
    final.norm <- if (identical(row.normalize, "l1")) {
        sum(abs(final))
    } else {
        sqrt(sum(final^2))
    }
    raw.used <- switch(row.normalize, none = raw.l2, l2 = raw.l2, l1 = raw.l1)
    list(ok = TRUE, row = final, row.norm.raw = raw.used,
         row.norm.raw.l2 = raw.l2, row.norm.raw.l1 = raw.l1,
         row.norm.used = row.normalize,
         row.norm.final = final.norm, drop.reason = NA_character_)
}

.lpl.tf.row.result <- function(row, row.table, local.design, h, predictors) {
    if (!length(h)) h <- rep(NA_real_, length(predictors))
    list(
        row = row,
        row.table = row.table,
        local.design = local.design,
        prediction.weights = data.frame(
            predictor.index = as.integer(predictors),
            h = as.double(h),
            stringsAsFactors = FALSE
        )
    )
}

.lpl.tf.assemble.sparse <- function(rows, ncol) {
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) {
        return(Matrix::Matrix(0, nrow = 0L, ncol = ncol, sparse = TRUE))
    }
    mat <- do.call(rbind, rows)
    methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
}

.lpl.tf.operator.diagnostics <- function(A, row.table, X, degree, chart.dim) {
    rank <- .lpl.tf.matrix.rank(A)
    probes <- if (identical(chart.dim, ncol(X))) {
        .lpl.tf.design.matrix(X, .lpl.tf.monomial.powers(ncol(X), degree))
    } else {
        matrix(numeric(0), nrow = nrow(X), ncol = 0L)
    }
    poly.residual <- if (nrow(A) > 0L && ncol(probes) > 0L) {
        norm(as.matrix(A %*% probes), type = "F") /
            max(norm(probes, type = "F"), .Machine$double.eps)
    } else {
        NA_real_
    }
    list(
        n.observations = nrow(X),
        n.candidate.rows = nrow(row.table),
        n.valid.rows = nrow(A),
        n.dropped.rows = sum(row.table$status != "ok"),
        dropped.by.reason = table(row.table$drop.reason, useNA = "ifany"),
        operator.rank = rank,
        operator.nullity = ncol(A) - rank,
        row.norm.raw.summary = summary(row.table$row.norm.raw),
        row.norm.final.summary = summary(row.table$row.norm.final),
        condition.number.summary = summary(row.table$condition.number),
        predictor.size.summary = summary(row.table$predictor.size),
        positive.weight.predictor.size.summary =
            summary(row.table$positive.weight.predictor.size),
        polynomial.residual = poly.residual,
        polynomial.probe.degree = degree,
        polynomial.probe.chart.dim = chart.dim,
        polynomial.probe.columns = ncol(probes)
    )
}

.lpl.tf.matrix.rank <- function(A, max.dense = 500L) {
    if (nrow(A) == 0L) return(0L)
    if (ncol(A) > max.dense || nrow(A) > max.dense) return(NA_integer_)
    as.integer(Matrix::rankMatrix(as.matrix(A), tol = 1e-8)[1])
}

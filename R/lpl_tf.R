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
#'   defaults to \code{ncol(X)}.  The special value \code{"auto"} estimates a
#'   single chart dimension from observed-coordinate local PCA spectra only,
#'   without using responses, truth values, latent coordinates, or labels.
#' @param support.metric Support-distance rule. \code{"coordinates"} uses
#'   Euclidean coordinate distances; \code{"graph.geodesic"} uses shortest-path
#'   distances in the supplied graph. \code{"auto"} uses graph geodesics when a
#'   graph is supplied and coordinate distances otherwise.
#' @param auto.chart.support.metric Support system used when
#'   \code{chart.dim = "auto"}. \code{"coordinates"} uses Euclidean coordinate
#'   neighborhoods, \code{"operator"} uses the resolved operator support metric,
#'   and \code{"both"} computes both diagnostics side by side.
#' @param auto.chart.selection.metric Which auto chart-dimension diagnostic to
#'   use for the fitted operator when \code{auto.chart.support.metric = "both"}.
#'   The default \code{"coordinates"} preserves historical behavior.
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
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator"),
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
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
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
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
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
            requested.chart.dim = args$requested.chart.dim,
            auto.chart.dim = args$auto.chart.dim,
            auto.chart.dim.diagnostics = args$auto.chart.dim.diagnostics,
            auto.chart.support.metric = args$auto.chart.support.metric,
            auto.chart.selection.metric = args$auto.chart.selection.metric,
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
            D.fallback = D.info$D.sparse,
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
        D.fallback = D.info$D.sparse,
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
            svd = path.info$genlasso.svd.used,
            svd.requested = path.info$genlasso.svd.requested,
            svd.fallback.used = path.info$genlasso.svd.fallback.used,
            svd.fallback.message = path.info$genlasso.svd.fallback.message,
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
                                  use.svd = FALSE, D.fallback = NULL) {
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
    run.genlasso <- function(D, svd) {
        suppressWarnings(genlasso::genlasso(
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
            svd = svd
        ))
    }
    primary <- tryCatch(
        list(path = run.genlasso(D, use.svd), error = NULL),
        error = function(e) list(path = NULL, error = e)
    )
    fallback.used <- FALSE
    fallback.message <- NA_character_
    svd.used <- use.svd
    if (!is.null(primary$error) && isTRUE(use.svd) &&
        .lpl.tf.is.svd.backend.error(primary$error)) {
        fallback.message <- conditionMessage(primary$error)
        path <- run.genlasso(D.fallback %||% D, FALSE)
        fallback.used <- TRUE
        svd.used <- FALSE
    } else if (!is.null(primary$error)) {
        stop(conditionMessage(primary$error), call. = FALSE)
    } else {
        path <- primary$path
    }
    list(
        path = path,
        n.observed = sum(observed),
        genlasso.svd.requested = use.svd,
        genlasso.svd.used = svd.used,
        genlasso.svd.fallback.used = fallback.used,
        genlasso.svd.fallback.message = fallback.message
    )
}

.lpl.tf.cv <- function(y, weights, D, D.fallback = NULL, lambda.grid, fold.id,
                       loss, selection, solver.args, use.svd = FALSE) {
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
                use.svd = use.svd,
                D.fallback = D.fallback
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
        auto.chart.support.metric = get.setting("auto.chart.support.metric"),
        auto.chart.selection.metric = get.setting("auto.chart.selection.metric"),
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
    support.metric, auto.chart.support.metric, auto.chart.selection.metric,
    exclude.self, row.normalize, local.solver,
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
    support.metric <- match.arg(support.metric, c("auto", "coordinates", "graph.geodesic"))
    auto.chart.support.metric <- match.arg(
        auto.chart.support.metric,
        c("coordinates", "operator", "both")
    )
    auto.chart.selection.metric <- match.arg(
        auto.chart.selection.metric,
        c("coordinates", "operator")
    )
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
    requested.chart.dim <- chart.dim
    auto.chart.dim <- FALSE
    auto.chart.dim.diagnostics <- NULL
    if (identical(chart.dim, "auto")) {
        if (!identical(coordinate.method, "local.pca")) {
            stop("'chart.dim = \"auto\"' is only supported when ",
                 "coordinate.method = 'local.pca'.", call. = FALSE)
        }
        auto.chart.dim <- TRUE
        auto.chart <- .local.pca.auto.chart.dim.with.metric(
            X = X,
            support.size = support.size,
            min.support = min.support,
            degree = degree,
            operator.distance.matrix = if (identical(support.metric, "graph.geodesic")) {
                graph.distance.matrix
            } else {
                NULL
            },
            operator.support.metric = support.metric,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        chart.dim <- auto.chart$chart.dim
        auto.chart.dim.diagnostics <- auto.chart
    } else if (is.null(chart.dim)) {
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
        requested.chart.dim = requested.chart.dim,
        auto.chart.dim = auto.chart.dim,
        auto.chart.dim.diagnostics = auto.chart.dim.diagnostics,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
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
        rank.tolerance = NA_real_,
        rank.tolerance.rule = NA_character_,
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
    row.base$rank.tolerance <- fit$rank.tolerance
    row.base$rank.tolerance.rule <- fit$rank.tolerance.rule
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
    tol.rule <- "max(dim(W_half_D)) * max(singular_value_1, 1) * .Machine$double.eps"
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
                    rank.tolerance = tol,
                    rank.tolerance.rule = tol.rule,
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
                    rank.tolerance = tol,
                    rank.tolerance.rule = tol.rule,
                    condition.number = condition.number,
                    solver.used = solver.used,
                    fallback.used = fallback.used,
                    drop.reason = "linear_solve_failed"))
    }
    list(ok = TRUE, h = h, rank = rank,
         rank.tolerance = tol,
         rank.tolerance.rule = tol.rule,
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

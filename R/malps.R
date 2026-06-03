#' Fit Model-Averaged Local Polynomial Smoothing
#'
#' Fits the current model-averaged local polynomial smoother (MALPS) on a
#' supplied coordinate matrix.  MALPS fits local polynomial regressions around
#' observed anchor points and averages all positive-weight local predictions at
#' each training point.  The current implementation supports supplied
#' coordinates, observed or supplied coordinate anchors, fixed support
#' parameters, deterministic cross-validation and exact dense GCV over support
#' parameters, local PCA chart coordinates,
#' ordinary and weighted new-response refits, coordinate-space new-point
#' prediction, linear-smoother diagnostics, and fixed-selection bootstrap
#' uncertainty summaries.  Graph-geodesic supports are available for training
#' and refit; graph-geodesic new-point prediction is deferred.
#'
#' Current development benchmarks suggest using the plain MALPS defaults for
#' production-style fits unless a specific diagnostic motivates a non-default
#' option: \code{support.selection = "cv"}, \code{robust.iterations = 0L},
#' \code{model.weight.rule = "none"}, and observed anchors
#' (\code{anchor.coordinates = NULL}).  Exact GCV selection
#' (\code{support.selection = "gcv"}) is a useful faster option for
#' fixed-weight, non-robust MALPS fits, especially exploratory continuous
#' smoothing runs, but it is not yet the universal default.  Binary responses
#' are treated as numeric conditional-expectation targets; fitted values are not
#' clipped inside the smoother.  Grid anchors, empirical-Bayes shrinkage, robust
#' local residual reweighting, and conservative/smoothed GCV selection remain
#' experimental controls rather than default recommendations.
#'
#' @param X Numeric coordinate matrix with one row per observation.
#' @param y Numeric response vector with length \code{nrow(X)}.
#' @param graph Optional supported gflow graph object.  When
#'   \code{support.metric = "graph.geodesic"} or \code{"auto"} with a graph
#'   supplied, weighted shortest-path distances from this graph are used for
#'   support construction.
#' @param adj.list Optional 1-based undirected adjacency list.  Must be supplied
#'   together with \code{weight.list} when using a supplied graph payload.
#' @param weight.list Optional positive edge-length list parallel to
#'   \code{adj.list}.  Edge weights are interpreted as graph lengths.
#' @param graph.stage Graph lifecycle stage used when \code{graph} is a gflow
#'   graph object; see \code{\link{graph.geodesic.distances}}.
#' @param anchor.index Optional integer vector of observed rows used as local
#'   model anchors.  Defaults to all rows when \code{anchor.coordinates} is
#'   \code{NULL}.
#' @param anchor.coordinates Optional numeric coordinate matrix of off-sample
#'   local model anchors with the same number of columns as \code{X}.  This
#'   option is currently supported only with coordinate support distances, not
#'   with \code{support.metric = "graph.geodesic"}.
#' @param degree Local polynomial degree, one of \code{0L}, \code{1L}, or
#'   \code{2L}.
#' @param degree.grid Optional degree grid used when
#'   \code{support.selection} is \code{"cv"} or \code{"gcv"}.
#' @param support.type Support construction rule.  \code{"knn"} uses exactly
#'   \code{support.size} nearest observations including the anchor itself for
#'   full-data local fitting supports.  Prediction supports are then rebuilt
#'   from the induced local radius, so tied targets at the kNN boundary can be
#'   covered even when they were not part of the exact-size fitting support.
#'   \code{"adaptive.radius"} uses a radius large enough to include at least
#'   \code{min.support} observations.  \code{"fixed.radius"} uses the supplied
#'   \code{radius}.
#' @param support.size Number of nearest observations for
#'   \code{support.type = "knn"}.  If \code{NULL}, defaults to
#'   \code{min.support}.
#' @param support.grid Optional kNN support-size grid used when
#'   \code{support.selection} is \code{"cv"} or \code{"gcv"}.
#' @param radius Radius for \code{support.type = "fixed.radius"} and optional
#'   base radius for \code{support.type = "adaptive.radius"}.
#' @param radius.grid Optional radius grid used when
#'   \code{support.selection} is \code{"cv"} or \code{"gcv"}.
#' @param min.support Minimum positive-weight support size.  If \code{NULL},
#'   defaults to the local polynomial design size plus \code{support.buffer}.
#' @param min.support.grid Optional minimum-support grid used when
#'   \code{support.selection} is \code{"cv"} or \code{"gcv"}.
#' @param support.buffer Nonnegative integer added to the local design size when
#'   deriving the default \code{min.support}.
#' @param kernel Kernel used for local fitting and model averaging.
#' @param kernel.grid Optional kernel grid used when
#'   \code{support.selection} is \code{"cv"} or \code{"gcv"}.
#' @param model.weight.rule Optional per-anchor model-quality multiplier applied
#'   to the ordinary kernel averaging weights.  \code{"none"} preserves the
#'   kernel-only averaging from earlier phases.  \code{"condition"} downweights
#'   poorly conditioned local designs with
#'   \eqn{q_u^{\mathrm{cond}}=1/\max(\kappa_u,1)}.  \code{"support"} upweights
#'   anchors with larger positive fitting supports using support size divided by
#'   the maximum support size.  \code{"boundary"} downweights asymmetric
#'   boundary-like supports with \eqn{q_u^{\mathrm{bdry}}=1/(1+b_u)}.
#'   \code{"quality"} multiplies these three component weights.  For every
#'   non-\code{"none"} rule, \code{diagnostics$model.weight.raw} stores the
#'   selected unnormalized rule multiplier, while \code{model.weights} and
#'   \code{diagnostics$model.weight} store the final multipliers after
#'   median-one normalization over positive raw values.  Raw zero or non-finite
#'   multipliers are retained diagnostically but are replaced by a tiny positive
#'   final floor before prediction averaging, so model-quality rules reweight
#'   existing prediction support rather than removing anchors from support
#'   membership.  These weights affect prediction averaging only; they do not
#'   change local supports or local coefficient estimation.
#' @param duplicate.action How duplicate coordinate rows should be handled.
#'   \code{"keep"} allows duplicates, records duplicate diagnostics, and
#'   preserves observed-anchor self-inclusion.  \code{"error"} rejects duplicate
#'   rows.  Jittering duplicates is deliberately not implemented in Phase 1.
#' @param coordinate.method Coordinate system used for each local polynomial
#'   design.  \code{"coordinates"} uses supplied coordinates centered at the
#'   anchor.  \code{"local.pca"} uses a deterministic local PCA chart estimated
#'   from the anchor support.
#' @param chart.dim Local chart dimension for \code{coordinate.method =
#'   "local.pca"}.  If \code{NULL}, defaults to \code{ncol(X)}.  The special
#'   value \code{"auto"} estimates a single observed-data local PCA dimension
#'   without using responses, truth values, latent coordinates, or labels.  For
#'   \code{coordinate.method = "coordinates"}, \code{chart.dim} must be
#'   \code{NULL}.
#' @param support.metric Distance system used for support construction.
#'   \code{"coordinates"} uses Euclidean distances in \code{X}.
#'   \code{"graph.geodesic"} uses weighted shortest-path distances from
#'   \code{graph} or \code{adj.list}/\code{weight.list}.  \code{"auto"} uses
#'   graph geodesic distances when graph input is supplied and coordinate
#'   distances otherwise.
#' @param auto.chart.support.metric Support system used when
#'   \code{chart.dim = "auto"}. \code{"coordinates"} uses Euclidean coordinate
#'   neighborhoods, \code{"operator"} uses the resolved MALPS support metric,
#'   and \code{"both"} computes both diagnostics side by side.
#' @param auto.chart.selection.metric Which auto chart-dimension diagnostic to
#'   use for the fitted MALPS model when both diagnostics are available. The
#'   default \code{"coordinates"} preserves historical behavior.
#' @param support.selection \code{"fixed"} for a single support profile,
#'   \code{"cv"} to tune support parameters by cross-validation, or
#'   \code{"gcv"} to tune support parameters by exact dense generalized
#'   cross-validation.  The GCV path is cached across candidates and is exact
#'   for fixed-weight linear MALPS fits; it rejects robust residual
#'   reweighting.
#' @param foldid Optional integer fold assignments.  If supplied,
#'   \code{cv.repeats} is forced to one.
#' @param cv.folds Number of folds generated when \code{foldid = NULL}.
#' @param cv.loss Cross-validation loss.
#' @param cv.repeats Number of independent fold assignments generated when
#'   \code{foldid = NULL}.
#' @param cv.seed Optional seed for generated folds.
#' @param cv.one.se Logical; use a one-standard-error rule for support
#'   selection.
#' @param gcv.exact.max.n Maximum number of observations allowed for exact dense
#'   GCV support selection.  Larger datasets should use \code{"cv"} until a
#'   sparse or trace-estimated GCV path is implemented.
#' @param local.solver Local weighted least-squares solver.
#'   \code{"auto"} uses normal equations only for full-rank, well-conditioned
#'   local designs and otherwise falls back to SVD.
#' @param normal.equations.max.condition Maximum local design condition number
#'   for normal equations under \code{local.solver = "auto"}.
#' @param robust.iterations Number of Cleveland-style robust residual
#'   reweighting iterations applied inside each local polynomial fit.  The
#'   default \code{0L} preserves ordinary MALPS.  Positive values multiply the
#'   geometric fitting weights by Tukey bisquare residual weights recomputed
#'   from the local fit.
#' @param robust.tuning.constant Positive bisquare tuning constant.  The
#'   default \code{6} matches the usual LOWESS convention in which residuals are
#'   scaled by \code{6 * median(abs(residuals))}.
#' @param verbose Logical; reserved for future progress messages.
#' @param ... Reserved for future extensions.  Supplying unused arguments is an
#'   error.
#'
#' @return A list of class \code{"malps"} with fitted values, local model
#'   coefficients, supports, averaging weights, diagnostics, and selection
#'   metadata.
#' @export
fit.malps <- function(
    X,
    y,
    graph = NULL,
    adj.list = NULL,
    weight.list = NULL,
    graph.stage = "final",
    anchor.index = NULL,
    anchor.coordinates = NULL,
    degree = 2L,
    degree.grid = NULL,
    support.type = c("adaptive.radius", "knn", "fixed.radius"),
    support.size = NULL,
    support.grid = NULL,
    radius = NULL,
    radius.grid = NULL,
    min.support = NULL,
    min.support.grid = NULL,
    support.buffer = 3L,
    kernel = c("epanechnikov", "triangular", "gaussian", "tricube"),
    kernel.grid = NULL,
    model.weight.rule = c("none", "condition", "support", "boundary", "quality"),
    duplicate.action = c("keep", "error"),
    coordinate.method = c("coordinates", "local.pca"),
    chart.dim = NULL,
    support.metric = c("auto", "coordinates", "graph.geodesic"),
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator"),
    support.selection = c("fixed", "cv", "gcv"),
    foldid = NULL,
    cv.folds = 5L,
    cv.loss = c("rmse", "mae", "mse"),
    cv.repeats = 1L,
    cv.seed = NULL,
    cv.one.se = FALSE,
    gcv.exact.max.n = 1000L,
    local.solver = c("auto", "normal.equations", "qr", "svd"),
    normal.equations.max.condition = 1e8,
    robust.iterations = 0L,
    robust.tuning.constant = 6,
    verbose = FALSE,
    ...) {

    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in fit.malps(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    X <- .malps.validate.X(X)
    n <- nrow(X)
    m <- ncol(X)
    y <- .malps.validate.y(y, n, "y")
    degree <- .malps.validate.degree(degree)
    support.type <- match.arg(support.type)
    kernel <- match.arg(kernel)
    model.weight.rule <- match.arg(model.weight.rule)
    cv.loss <- match.arg(cv.loss)
    duplicate.action <- match.arg(duplicate.action)
    support.metric <- match.arg(support.metric)
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    graph.input.supplied <- !is.null(graph) || !is.null(adj.list) ||
        !is.null(weight.list)
    support.metric <- .malps.resolve.support.metric(
        support.metric = support.metric,
        graph.input.supplied = graph.input.supplied
    )
    graph.info <- if (identical(support.metric, "graph.geodesic")) {
        .malps.prepare.graph.info(
            n = n,
            graph = graph,
            adj.list = adj.list,
            weight.list = weight.list,
            graph.stage = graph.stage
        )
    } else {
        NULL
    }
    duplicate.info <- .malps.duplicate.row.info(X)
    if (duplicate.info$has.duplicate.rows &&
        identical(duplicate.action, "error")) {
        stop(sprintf(
            paste0("X contains %d duplicate coordinate row(s). ",
                   "Duplicate rows are allowed with duplicate.action = 'keep'; ",
                   "jittering duplicates is not implemented in Phase 1a."),
            duplicate.info$n.duplicate.rows
        ), call. = FALSE)
    }
    coordinate.method <- match.arg(coordinate.method)
    auto.chart.dim.diagnostics <- NULL
    if (identical(chart.dim, "auto")) {
        if (!identical(coordinate.method, "local.pca")) {
            stop("'chart.dim = \"auto\"' is only supported when ",
                 "coordinate.method = 'local.pca'.", call. = FALSE)
        }
        auto.support.size <- support.size
        if (is.null(auto.support.size) && !is.null(support.grid)) {
            auto.support.size <- min(as.integer(support.grid), na.rm = TRUE)
        }
        auto.min.support <- min.support
        if (is.null(auto.min.support) && !is.null(min.support.grid)) {
            auto.min.support <- min(as.integer(min.support.grid), na.rm = TRUE)
        }
        auto.degree <- max(c(degree, degree.grid %||% degree), na.rm = TRUE)
        auto.distance.matrix <- NULL
        if (identical(support.metric, "graph.geodesic")) {
            auto.distance.matrix <- shortest.path(
                graph.info$adj.list,
                graph.info$weight.list,
                seq_along(graph.info$adj.list)
            )
            if (!is.matrix(auto.distance.matrix) ||
                any(dim(auto.distance.matrix) != c(n, n)) ||
                any(!is.finite(auto.distance.matrix))) {
                stop("Graph-geodesic auto chart dimension requires a connected graph.",
                     call. = FALSE)
            }
        }
        auto.chart <- .local.pca.auto.chart.dim.with.metric(
            X = X,
            support.size = auto.support.size,
            min.support = auto.min.support,
            degree = auto.degree,
            operator.distance.matrix = auto.distance.matrix,
            operator.support.metric = support.metric,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        chart.dim <- auto.chart$chart.dim
        auto.chart.dim.diagnostics <- auto.chart
    }
    chart.dim <- .malps.validate.chart.dim(chart.dim, m, coordinate.method)
    support.selection <- match.arg(support.selection)
    local.solver <- match.arg(local.solver)
    support.buffer <- .malps.validate.nonnegative.integer(support.buffer,
                                                         "support.buffer")
    cv.folds <- .malps.validate.positive.integer(cv.folds, "cv.folds")
    cv.repeats <- .malps.validate.positive.integer(cv.repeats, "cv.repeats")
    if (!is.null(cv.seed)) {
        cv.seed <- .malps.validate.integer.scalar(cv.seed, "cv.seed")
    }
    if (!is.logical(cv.one.se) || length(cv.one.se) != 1L ||
        is.na(cv.one.se)) {
        stop("cv.one.se must be TRUE or FALSE.", call. = FALSE)
    }
    gcv.exact.max.n <- .malps.validate.positive.integer(
        gcv.exact.max.n, "gcv.exact.max.n"
    )
    normal.equations.max.condition <- .malps.validate.positive.scalar(
        normal.equations.max.condition, "normal.equations.max.condition"
    )
    robust.iterations <- .malps.validate.nonnegative.integer(
        robust.iterations, "robust.iterations"
    )
    robust.tuning.constant <- .malps.validate.positive.scalar(
        robust.tuning.constant, "robust.tuning.constant"
    )
    anchor.info <- .malps.prepare.anchors(
        X = X,
        anchor.index = anchor.index,
        anchor.coordinates = anchor.coordinates,
        support.metric = support.metric
    )
    anchor.index <- anchor.info$anchor.index
    anchors <- anchor.info$anchors
    anchor.labels <- anchor.info$anchor.labels
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.", call. = FALSE)
    }
    if (identical(support.selection, "gcv") && robust.iterations > 0L) {
        stop(
            paste(
                "support.selection = \"gcv\" requires robust.iterations = 0",
                "because exact MALPS GCV is defined for fixed-weight linear",
                "fits."
            ),
            call. = FALSE
        )
    }

    support.distances <- .malps.support.distance.matrix(
        X = X,
        anchors = anchors,
        anchor.index = anchor.index,
        support.metric = support.metric,
        graph.info = graph.info
    )
    distance.matrix <- support.distances$distance.matrix
    if (identical(support.selection, "cv")) {
        selection <- .malps.select.cv(
            X = X,
            y = y,
            anchor.index = anchor.index,
            anchors = anchors,
            anchor.labels = anchor.labels,
            distance.matrix = distance.matrix,
            support.type = support.type,
            degree = degree,
            degree.grid = degree.grid,
            support.size = support.size,
            support.grid = support.grid,
            radius = radius,
            radius.grid = radius.grid,
            min.support = min.support,
            min.support.grid = min.support.grid,
            support.buffer = support.buffer,
            kernel = kernel,
            kernel.grid = kernel.grid,
            model.weight.rule = model.weight.rule,
            coordinate.method = coordinate.method,
            chart.dim = chart.dim,
            support.metric = support.metric,
            foldid = foldid,
            cv.folds = cv.folds,
            cv.loss = cv.loss,
            cv.repeats = cv.repeats,
            cv.seed = cv.seed,
            cv.one.se = cv.one.se,
            local.solver = local.solver,
            normal.equations.max.condition = normal.equations.max.condition,
            robust.iterations = robust.iterations,
            robust.tuning.constant = robust.tuning.constant
        )
        params <- selection$selected.params
    } else if (identical(support.selection, "gcv")) {
        selection <- .malps.select.gcv(
            X = X,
            y = y,
            anchor.index = anchor.index,
            anchors = anchors,
            anchor.labels = anchor.labels,
            distance.matrix = distance.matrix,
            support.type = support.type,
            degree = degree,
            degree.grid = degree.grid,
            support.size = support.size,
            support.grid = support.grid,
            radius = radius,
            radius.grid = radius.grid,
            min.support = min.support,
            min.support.grid = min.support.grid,
            support.buffer = support.buffer,
            kernel = kernel,
            kernel.grid = kernel.grid,
            model.weight.rule = model.weight.rule,
            coordinate.method = coordinate.method,
            chart.dim = chart.dim,
            support.metric = support.metric,
            graph.info = graph.info,
            duplicate.action = duplicate.action,
            duplicate.info = duplicate.info,
            local.solver = local.solver,
            normal.equations.max.condition = normal.equations.max.condition,
            robust.tuning.constant = robust.tuning.constant,
            gcv.exact.max.n = gcv.exact.max.n
        )
        params <- selection$selected.params
    } else {
        params <- .malps.prepare.fixed.params(
            n = n,
            m = m,
            degree = degree,
            support.type = support.type,
            support.size = support.size,
            radius = radius,
            min.support = min.support,
            support.buffer = support.buffer,
            chart.dim = chart.dim,
            support.metric = support.metric,
            model.weight.rule = model.weight.rule,
            kernel = kernel
        )
        selection <- list(
            method = "fixed",
            selected.index = 1L,
            selected.params = params,
            loss.name = NULL,
            one.se.used = FALSE,
            fallback.method = NULL,
            n.valid = NA_integer_,
            cv.folds = NULL,
            cv.repeats.requested = NULL,
            cv.repeats.effective = NULL,
            foldid.source = NULL,
            seed = NULL,
            fold.message = NULL
        )
    }
    out <- .malps.fit.core(
        X = X,
        y = y,
        anchor.index = anchor.index,
        anchors = anchors,
        anchor.labels = anchor.labels,
        distance.matrix = distance.matrix,
        params = params,
        coordinate.method = coordinate.method,
        support.metric = support.metric,
        graph.info = graph.info,
        model.weight.rule = model.weight.rule,
        duplicate.action = duplicate.action,
        duplicate.info = duplicate.info,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant,
        selection = selection,
        call = match.call()
    )
    out$auto.chart.dim.diagnostics <- auto.chart.dim.diagnostics
    out$auto.chart.support.metric <- auto.chart.support.metric
    out$auto.chart.selection.metric <- auto.chart.selection.metric
    out
}

#' Predict From A MALPS Fit
#'
#' @param object A \code{"malps"} object from \code{\link{fit.malps}}.
#' @param newdata Optional numeric coordinate matrix with the same number of
#'   columns as \code{object$X}.  If \code{NULL}, returns stored training-point
#'   fitted values.  New-point prediction is currently implemented for
#'   coordinate-support fits only; it fails informatively for objects fitted
#'   with \code{support.metric = "graph.geodesic"} because new targets need a
#'   graph-attachment/geodesic-distance contract.
#' @param type Prediction type.  Currently only \code{"response"}.
#' @param allow.incomplete Logical; if \code{FALSE}, new-point prediction fails
#'   when any target has no positive averaging coverage.  If \code{TRUE}, those
#'   targets receive \code{NA_real_}.
#' @param ... Reserved for future extensions.
#'
#' @return Numeric vector of fitted values.
#' @export
predict.malps <- function(object, newdata = NULL, type = c("response"),
                          allow.incomplete = FALSE, ...) {
    if (!inherits(object, "malps")) {
        stop("object must be a 'malps' object.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in predict.malps(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    type <- match.arg(type)
    if (!identical(type, "response")) {
        stop("predict.malps() supports only type = 'response'.",
             call. = FALSE)
    }
    if (!is.logical(allow.incomplete) || length(allow.incomplete) != 1L ||
        is.na(allow.incomplete)) {
        stop("allow.incomplete must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(newdata)) {
        if (identical(object$support.metric, "graph.geodesic")) {
            stop(
                paste(
                    "predict.malps(newdata = ...) is not implemented for",
                    "fits using support.metric = 'graph.geodesic'; new targets",
                    "need a graph attachment/geodesic-distance contract."
                ),
                call. = FALSE
            )
        }
        newdata <- .malps.validate.newdata(newdata, ncol(object$X))
        return(.malps.predict.newdata(
            object = object,
            newdata = newdata,
            allow.incomplete = allow.incomplete
        ))
    }
    object$fitted.values
}

#' Refit A MALPS Smoother With A New Response
#'
#' Reuses the anchors, selected support profile, local supports, prediction
#' supports, coordinate method, kernel, degree, and local-solver controls from a
#' previous \code{\link{fit.malps}} result, then recomputes local polynomial
#' coefficients for a new response vector.  Observation weights, when supplied,
#' multiply the stored geometric fitting weights for local coefficient
#' estimation.  Supports and model-averaging weights remain fixed.  A weighted
#' refit requires every stored local support to retain at least one positive
#' effective fitting weight; sparse bootstrap-style weights can therefore fail
#' even when the global weight vector is not all zero.  If the original fit used
#' robust local residual reweighting, robust weights are recomputed for the new
#' response because they are response-dependent.
#'
#' @param object A \code{"malps"} object from \code{\link{fit.malps}}.
#' @param y Optional new numeric response vector with length
#'   \code{nrow(object$X)}.  If \code{NULL}, the original response is reused.
#' @param weights Optional nonnegative numeric observation weights with length
#'   \code{nrow(object$X)}.  Zero weights remove observations from local
#'   coefficient estimation but do not change the stored supports or averaging
#'   weights.  The refit fails if any local support has no positive effective
#'   fitting weight after multiplying by these observation weights.
#' @param reuse.selection Currently requires \code{TRUE}.
#' @param refit.local.coefficients Currently requires \code{TRUE}.
#' @param verbose Logical; reserved for future progress messages.
#' @param ... Reserved for future extensions.
#'
#' @return A \code{"malps"} object with the same support profile and new fitted
#'   values.
#' @export
refit.malps <- function(
    object,
    y = NULL,
    weights = NULL,
    reuse.selection = TRUE,
    refit.local.coefficients = TRUE,
    verbose = FALSE,
    ...) {

    if (!inherits(object, "malps")) {
        stop("object must be a 'malps' object.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in refit.malps(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    weights <- .malps.validate.refit.weights(weights, nrow(object$X))
    if (!isTRUE(reuse.selection)) {
        stop("refit.malps() currently requires reuse.selection = TRUE.",
             call. = FALSE)
    }
    if (!isTRUE(refit.local.coefficients)) {
        stop("refit.malps() currently requires refit.local.coefficients = TRUE.",
             call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.", call. = FALSE)
    }
    if (is.null(y)) {
        y <- object$y
    }
    y <- .malps.validate.y(y, nrow(object$X), "y")

    local.solver <- object$local.solver %||% "auto"
    normal.equations.max.condition <-
        object$normal.equations.max.condition %||% 1e8
    robust.iterations <- object$robust.iterations %||% 0L
    robust.tuning.constant <- object$robust.tuning.constant %||% 6
    refit <- .malps.refit.with.supports(
        object = object,
        y = y,
        weights = weights,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant,
        call = match.call()
    )
    refit$refit <- list(
        source.class = class(object),
        reuse.selection = TRUE,
        refit.local.coefficients = TRUE,
        weights = weights
    )
    refit
}

#' Construct The Conditional MALPS Smoother Matrix
#'
#' Constructs the exact training-point smoother matrix for a fitted
#' \code{\link{fit.malps}} object, conditional on the stored supports, local
#' charts, fitting weights, and averaging weights.  For a non-robust MALPS fit
#' or weighted refit with fixed supports, the returned matrix \eqn{S} satisfies
#' \eqn{\hat y = S y}.  If the original fit used cross-validation, this matrix
#' is conditional on the selected support profile; it does not account for the
#' response-dependence of the model-selection step.
#'
#' Robust residual reweighting makes the effective fitting weights depend on
#' the response.  By default this function therefore rejects robust MALPS fits.
#' Setting \code{allow.robust = TRUE} constructs the fixed-final-weight
#' linearization, which reproduces the stored fit but should not be interpreted
#' as the exact response-to-fit map.
#'
#' @param object A \code{"malps"} object from \code{\link{fit.malps}} or
#'   \code{\link{refit.malps}}.
#' @param max.n Maximum number of training observations for which a dense
#'   smoother matrix may be constructed.  Use \code{Inf} to disable this guard.
#' @param allow.robust Logical; allow fixed-final-weight linearization for fits
#'   with \code{robust.iterations > 0L}.
#' @param ... Reserved for future extensions.
#'
#' @return A dense numeric \eqn{n \times n} matrix with attributes describing
#'   whether the matrix is conditional on support selection and whether it used
#'   a robust fixed-weight linearization.
#' @export
malps.smoother.matrix <- function(object, max.n = 1000L,
                                  allow.robust = FALSE, ...) {
    if (!inherits(object, "malps")) {
        stop("object must be a 'malps' object.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in malps.smoother.matrix(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    if (!is.logical(allow.robust) || length(allow.robust) != 1L ||
        is.na(allow.robust)) {
        stop("allow.robust must be TRUE or FALSE.", call. = FALSE)
    }
    n <- nrow(object$X)
    if (!is.numeric(max.n) || length(max.n) != 1L || is.na(max.n) ||
        max.n < 1) {
        stop("max.n must be a positive numeric scalar or Inf.", call. = FALSE)
    }
    if (is.finite(max.n) && n > max.n) {
        stop(sprintf(
            paste0("Dense MALPS smoother matrix requested for n = %d, ",
                   "which exceeds max.n = %d."),
            n, as.integer(max.n)
        ), call. = FALSE)
    }
    robust.used <- (object$robust.iterations %||% 0L) > 0L
    if (robust.used && !allow.robust) {
        stop(
            paste(
                "malps.smoother.matrix() is exact only for fixed-weight",
                "linear MALPS fits. Robust residual weights are",
                "response-dependent; use allow.robust = TRUE only for the",
                "fixed-final-weight linearization."
            ),
            call. = FALSE
        )
    }

    n.anchors <- length(object$supports)
    S <- matrix(0, nrow = n, ncol = n)
    for (a in seq_len(n.anchors)) {
        fit.support <- object$supports[[a]]
        pred.support <- object$prediction.supports[[a]]
        if (!length(pred.support$index)) {
            next
        }
        fit.index <- fit.support$index
        fit.weights <- object$fit.weights[[a]]
        if (length(fit.weights) != length(fit.index)) {
            stop("Stored MALPS fit weights are inconsistent with supports.",
                 call. = FALSE)
        }
        Z.fit <- .malps.chart.coordinates(
            object$X[fit.index, , drop = FALSE],
            object$charts[[a]]
        )
        D.fit <- .malps.design.matrix(Z.fit, object$degree)
        coef.map <- .malps.wls.coefficient.map(D.fit, fit.weights)

        pred.index <- pred.support$index
        Z.pred <- .malps.chart.coordinates(
            object$X[pred.index, , drop = FALSE],
            object$charts[[a]]
        )
        D.pred <- .malps.design.matrix(Z.pred, object$degree)
        local.map <- D.pred %*% coef.map
        alpha <- object$averaging.weights[pred.index, a]
        for (r in seq_along(pred.index)) {
            S[pred.index[r], fit.index] <-
                S[pred.index[r], fit.index] + alpha[r] * local.map[r, ]
        }
    }
    attr(S, "conditional.on.selection") <-
        !is.null(object$selection) &&
        !identical(object$selection$method %||% "fixed", "fixed")
    attr(S, "robust.fixed.weight.linearization") <- robust.used
    attr(S, "degree") <- object$degree
    attr(S, "support.type") <- object$support.type
    attr(S, "support.metric") <- object$support.metric
    S
}

#' Compute MALPS Linear-Smoother Diagnostics
#'
#' Computes effective degrees of freedom, generalized cross-validation (GCV),
#' and analytic leave-one-out residual diagnostics from a MALPS smoother matrix.
#' These diagnostics are exact for fixed-support, fixed-weight linear MALPS
#' fits.  For cross-validated fits they are conditional on the selected support
#' profile.  Robust fits are rejected by default for the same reason described
#' in \code{\link{malps.smoother.matrix}}.
#'
#' The GCV score is computed as
#' \deqn{
#'   \mathrm{GCV} =
#'   \frac{n^{-1}\|y - S y\|_2^2}
#'        {(1 - \operatorname{tr}(S)/n)^2}.
#' }
#' The analytic leave-one-out residuals are computed as
#' \deqn{
#'   e_i^{\mathrm{loo}} =
#'   \frac{y_i - \hat y_i}{1 - S_{ii}}.
#' }
#'
#' @param object A \code{"malps"} object from \code{\link{fit.malps}} or
#'   \code{\link{refit.malps}}.
#' @param y Optional response vector.  Defaults to \code{object$y}.
#' @param smoother.matrix Optional precomputed matrix from
#'   \code{\link{malps.smoother.matrix}}.
#' @param include.loocv Logical; include analytic leave-one-out residuals and
#'   mean squared error.
#' @param max.n Maximum dense smoother size passed to
#'   \code{\link{malps.smoother.matrix}} when \code{smoother.matrix = NULL}.
#' @param allow.robust Logical; passed to \code{\link{malps.smoother.matrix}}.
#' @param ... Reserved for future extensions.
#'
#' @return A list with residual, fitted-value, EDF, GCV, and optional LOOCV
#'   diagnostics.
#' @export
malps.gcv <- function(object, y = NULL, smoother.matrix = NULL,
                      include.loocv = TRUE, max.n = 1000L,
                      allow.robust = FALSE, ...) {
    if (!inherits(object, "malps")) {
        stop("object must be a 'malps' object.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in malps.gcv(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    n <- nrow(object$X)
    if (is.null(y)) {
        y <- object$y
    }
    y <- .malps.validate.y(y, n, "y")
    if (is.null(smoother.matrix)) {
        smoother.matrix <- malps.smoother.matrix(
            object,
            max.n = max.n,
            allow.robust = allow.robust
        )
    } else if (!is.matrix(smoother.matrix) ||
               !identical(dim(smoother.matrix), c(n, n)) ||
               !is.numeric(smoother.matrix)) {
        stop("smoother.matrix must be a numeric nrow(X) by nrow(X) matrix.",
             call. = FALSE)
    }
    if (!is.logical(include.loocv) || length(include.loocv) != 1L ||
        is.na(include.loocv)) {
        stop("include.loocv must be TRUE or FALSE.", call. = FALSE)
    }

    fitted <- as.numeric(smoother.matrix %*% y)
    residuals <- y - fitted
    rss <- sum(residuals^2)
    mse <- mean(residuals^2)
    leverage <- diag(smoother.matrix)
    edf <- sum(leverage)
    denom <- 1 - edf / n
    gcv <- if (is.finite(denom) && abs(denom) > sqrt(.Machine$double.eps)) {
        mse / denom^2
    } else {
        NA_real_
    }

    out <- list(
        n = n,
        fitted.values = fitted,
        residuals = residuals,
        rss = rss,
        mse = mse,
        edf = edf,
        leverage = leverage,
        gcv = gcv,
        conditional.on.selection =
            isTRUE(attr(smoother.matrix, "conditional.on.selection")),
        robust.fixed.weight.linearization =
            isTRUE(attr(smoother.matrix,
                        "robust.fixed.weight.linearization"))
    )
    if (include.loocv) {
        loo.denom <- 1 - leverage
        invalid <- !is.finite(loo.denom) |
            abs(loo.denom) <= sqrt(.Machine$double.eps)
        loo.resid <- residuals / loo.denom
        loo.resid[invalid] <- NA_real_
        out$loocv.residuals <- loo.resid
        out$loocv.invalid <- which(invalid)
        out$loocv.mse <- if (any(invalid)) NA_real_ else mean(loo.resid^2)
    }
    class(out) <- c("malps_gcv", "list")
    out
}

#' Bootstrap Uncertainty For A Fixed-Selection MALPS Fit
#'
#' Reuses the supports, local charts, prediction supports, averaging weights,
#' selected degree, and local-solver controls from an existing
#' \code{\link{fit.malps}} object, then repeatedly calls
#' \code{\link{refit.malps}} with bootstrap-style case weights.  This is a
#' fixed-selection uncertainty diagnostic: it measures refit variability
#' conditional on the fitted MALPS support profile and does not rerun
#' cross-validation or rebuild supports in each replicate.
#'
#' Two weight generators are available.  \code{"bayesian"} draws positive
#' exponential weights and rescales them to have mean one.  This is the default
#' because it preserves positive local support weights and therefore avoids many
#' local-support degeneracies.  \code{"multinomial"} draws ordinary
#' nonparametric bootstrap counts with total count \code{nrow(object$X)}; this
#' can fail when a replicate assigns zero weight to all observations in a local
#' support, and such failures are recorded.
#'
#' @param object A \code{"malps"} object from \code{\link{fit.malps}} or
#'   \code{\link{refit.malps}}.
#' @param B Number of successful bootstrap replicates requested.
#' @param weight.type Bootstrap weight generator, one of \code{"bayesian"} or
#'   \code{"multinomial"}.
#' @param y Optional response vector.  Defaults to \code{object$y}.
#' @param probs Optional nonnegative sampling probabilities for
#'   \code{weight.type = "multinomial"}.
#' @param conf.level Pointwise interval level used in the returned summary.
#' @param seed Optional integer seed for reproducible bootstrap weights.
#' @param max.failures Maximum failed refit attempts allowed before stopping.
#' @param keep.weights Logical; store the generated bootstrap weights in the
#'   returned object.
#' @param verbose Logical; emit progress messages every ten successful
#'   replicates.
#' @param ... Reserved for future extensions.
#'
#' @return A list of class \code{"malps_bootstrap"} containing:
#' \itemize{
#'   \item \code{replicate.fitted.values}: an \eqn{n \times B} matrix of
#'     successful replicate fitted values;
#'   \item \code{replicate.weights}: an optional \eqn{n \times B} matrix of
#'     bootstrap weights when \code{keep.weights = TRUE};
#'   \item \code{weights}: a compatibility alias for \code{replicate.weights};
#'   \item \code{failures}: a data frame with columns \code{attempt},
#'     \code{failure}, and \code{message};
#'   \item \code{summary}: pointwise fitted values, bootstrap summaries, and
#'     empirical lower/upper intervals;
#'   \item bootstrap metadata such as \code{B.completed}, \code{n.failures},
#'     \code{conf.level}, and \code{weight.type}.
#' }
#' @export
bootstrap.malps <- function(object, B = 200L,
                            weight.type = c("bayesian", "multinomial"),
                            y = NULL, probs = NULL, conf.level = 0.95,
                            seed = NULL, max.failures = max(10L, B),
                            keep.weights = FALSE, verbose = FALSE, ...) {
    if (!inherits(object, "malps")) {
        stop("object must be a 'malps' object.", call. = FALSE)
    }
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments in bootstrap.malps(): ",
             paste(names(dots), collapse = ", "), call. = FALSE)
    }
    B <- .malps.validate.positive.integer(B, "B")
    weight.type <- match.arg(weight.type)
    n <- nrow(object$X)
    if (is.null(y)) {
        y <- object$y
    }
    y <- .malps.validate.y(y, n, "y")
    conf.level <- .malps.validate.proportion(conf.level, "conf.level",
                                             strict = TRUE)
    max.failures <- .malps.validate.nonnegative.integer(max.failures,
                                                        "max.failures")
    if (!is.logical(keep.weights) || length(keep.weights) != 1L ||
        is.na(keep.weights)) {
        stop("keep.weights must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(seed)) {
        seed <- .malps.validate.integer.scalar(seed, "seed")
        set.seed(seed)
    }
    probs <- .malps.validate.bootstrap.probs(probs, n, weight.type)

    fitted.reps <- matrix(NA_real_, nrow = n, ncol = B)
    weights.reps <- if (keep.weights) matrix(NA_real_, nrow = n, ncol = B) else NULL
    failures <- data.frame(
        attempt = integer(0L),
        failure = integer(0L),
        message = character(0L),
        stringsAsFactors = FALSE
    )
    completed <- 0L
    attempts <- 0L
    failed <- 0L
    while (completed < B) {
        attempts <- attempts + 1L
        weights <- .malps.bootstrap.weights(n, weight.type, probs)
        refit <- tryCatch(
            refit.malps(object, y = y, weights = weights),
            error = function(e) e
        )
        if (inherits(refit, "error")) {
            failed <- failed + 1L
            failures <- rbind(
                failures,
                data.frame(
                    attempt = attempts,
                    failure = failed,
                    message = conditionMessage(refit),
                    stringsAsFactors = FALSE
                )
            )
            if (failed > max.failures) {
                stop(sprintf(
                    paste0("bootstrap.malps() exceeded max.failures = %d ",
                           "before completing B = %d replicates. Last error: %s"),
                    max.failures, B, conditionMessage(refit)
                ), call. = FALSE)
            }
            next
        }
        completed <- completed + 1L
        fitted.reps[, completed] <- refit$fitted.values
        if (keep.weights) {
            weights.reps[, completed] <- weights
        }
        if (verbose && (completed %% 10L == 0L || completed == B)) {
            message(sprintf("bootstrap.malps(): completed %d/%d replicates",
                            completed, B))
        }
    }
    colnames(fitted.reps) <- paste0("replicate.", seq_len(B))
    if (keep.weights) {
        colnames(weights.reps) <- colnames(fitted.reps)
    }
    alpha <- (1 - conf.level) / 2
    lower.prob <- alpha
    upper.prob <- 1 - alpha
    point.mean <- rowMeans(fitted.reps)
    point.sd <- if (B > 1L) {
        apply(fitted.reps, 1L, stats::sd)
    } else {
        rep(0, n)
    }
    point.median <- apply(fitted.reps, 1L, stats::median)
    point.lower <- apply(fitted.reps, 1L, stats::quantile,
                         probs = lower.prob, names = FALSE,
                         type = 7L)
    point.upper <- apply(fitted.reps, 1L, stats::quantile,
                         probs = upper.prob, names = FALSE,
                         type = 7L)
    summary <- data.frame(
        index = seq_len(n),
        fitted = object$fitted.values,
        mean = point.mean,
        sd = point.sd,
        median = point.median,
        lower = point.lower,
        upper = point.upper,
        stringsAsFactors = FALSE
    )
    out <- list(
        fitted.values = object$fitted.values,
        replicate.fitted.values = fitted.reps,
        summary = summary,
        replicate.weights = weights.reps,
        weights = weights.reps,
        weight.type = weight.type,
        B.requested = B,
        B.completed = completed,
        attempts = attempts,
        failures = failures,
        n.failures = failed,
        conf.level = conf.level,
        seed = seed,
        y = y,
        call = match.call(),
        object.call = object$call
    )
    class(out) <- c("malps_bootstrap", "list")
    out
}

#' @export
print.malps <- function(x, ...) {
    cat("Model-averaged local polynomial smoothing fit\n")
    cat("  n:", length(x$y), "\n")
    cat("  anchors:", length(x$anchor.index), "\n")
    cat("  degree:", x$degree, "\n")
    cat("  support.metric:", x$support.metric %||% "coordinates", "\n")
    cat("  support.type:", x$support.type, "\n")
    cat("  kernel:", x$kernel, "\n")
    if (!is.null(x$diagnostics)) {
        cat("  min coverage:", x$diagnostics$min.coverage.count, "\n")
        cat("  finite predictions:", x$diagnostics$finite.predictions, "\n")
    }
    invisible(x)
}

#' @export
print.malps_bootstrap <- function(x, ...) {
    cat("MALPS fixed-selection bootstrap\n")
    cat("  B completed:", x$B.completed, "\n")
    cat("  weight.type:", x$weight.type, "\n")
    cat("  failures:", x$n.failures, "\n")
    cat("  conf.level:", x$conf.level, "\n")
    invisible(x)
}

.malps.fit.core <- function(X, y, anchor.index, anchors, anchor.labels,
                            distance.matrix, params,
                            coordinate.method, support.metric, graph.info,
                            model.weight.rule,
                            duplicate.action,
                            duplicate.info, local.solver,
                            normal.equations.max.condition,
                            robust.iterations, robust.tuning.constant,
                            selection, call,
                            fit.index = seq_len(nrow(X)),
                            target.index = seq_len(nrow(X)),
                            fail.on.no.coverage = TRUE) {
    model.weight.rule <- params$model.weight.rule %||%
        model.weight.rule %||% "none"
    supports <- .malps.build.supports(
        distance.matrix = distance.matrix,
        anchor.index = anchor.index,
        anchor.labels = anchor.labels,
        support.type = params$support.type,
        support.size = params$support.size,
        radius = params$radius,
        min.support = params$min.support,
        kernel = params$kernel,
        candidate.index = fit.index
    )
    charts <- .malps.build.charts(
        X = X,
        anchor.index = anchor.index,
        anchors = anchors,
        anchor.labels = anchor.labels,
        supports = supports,
        coordinate.method = coordinate.method,
        chart.dim = params$chart.dim
    )
    fits <- .malps.fit.local.models(
        X = X,
        y = y,
        anchor.index = anchor.index,
        supports = supports,
        charts = charts,
        degree = params$degree,
        observation.weights = NULL,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant
    )
    boundary.score <- .malps.boundary.score(
        X = X,
        anchor.index = anchor.index,
        anchors = anchors,
        supports = supports
    )
    model.weight.info <- .malps.model.quality.diagnostics(
        rule = model.weight.rule,
        condition.number = fits$model.condition.number,
        positive.support.size = vapply(
            supports,
            function(s) sum(s$weights > 0),
            integer(1L)
        ),
        boundary.score = boundary.score
    )
    model.weights <- model.weight.info$weight
    prediction.supports <- .malps.build.prediction.supports(
        distance.matrix = distance.matrix,
        anchor.index = anchor.index,
        anchor.labels = anchor.labels,
        supports = supports,
        kernel = params$kernel,
        model.weights = model.weights,
        target.index = target.index
    )
    averaged <- .malps.average.predictions(
        X = X,
        anchor.labels = anchor.labels,
        supports = prediction.supports,
        charts = charts,
        coefficients = fits$coefficients,
        degree = params$degree
    )
    target.uncovered <- intersect(averaged$uncovered.index, target.index)
    if (fail.on.no.coverage && length(target.uncovered) > 0L) {
        uncovered.label <- .malps.format.index.message(target.uncovered)
        stop(sprintf(
            paste0("MALPS fit has %d target(s) with no positive averaging ",
                   "coverage; min.coverage.count = %d; uncovered.index = %s."),
            length(target.uncovered),
            min(averaged$coverage.count[target.index]),
            uncovered.label
        ), call. = FALSE)
    }
    if (fail.on.no.coverage &&
        !all(is.finite(averaged$fitted.values[target.index]))) {
        stop("MALPS produced non-finite fitted values.", call. = FALSE)
    }

    diagnostics <- list(
        required.support = .malps.design.size(params$chart.dim, params$degree),
        chart.dim = params$chart.dim,
        support.metric = support.metric,
        model.weight.rule = model.weight.rule,
        model.weight = model.weights,
        model.weight.raw = model.weight.info$raw.weight,
        model.weight.floor = model.weight.info$floor,
        model.weight.condition = model.weight.info$condition.weight,
        model.weight.support = model.weight.info$support.weight,
        model.weight.boundary = model.weight.info$boundary.weight,
        model.weight.source.condition.number = fits$model.condition.number,
        model.weight.fixed.from.original = FALSE,
        graph.source = if (is.null(graph.info)) NA_character_ else graph.info$source,
        min.support = params$min.support,
        support.size = vapply(supports, function(s) length(s$index),
                              integer(1L)),
        positive.weight.support.size = vapply(
            supports,
            function(s) sum(s$weights > 0),
            integer(1L)
        ),
        radius = vapply(supports, function(s) s$radius, numeric(1L)),
        coverage.count = averaged$coverage.count,
        min.coverage.count = min(averaged$coverage.count[target.index]),
        n.no.coverage = length(target.uncovered),
        uncovered.index = target.uncovered,
        rank = fits$rank,
        condition.number = fits$condition.number,
        rank.deficient = fits$rank < .malps.design.size(params$chart.dim,
                                                        params$degree),
        solver.used = fits$solver.used,
        fallback.used = fits$fallback.used,
        model.condition.number = fits$model.condition.number,
        effective.fit.support.size = fits$effective.fit.support.size,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant,
        robust.used = robust.iterations > 0L,
        robust.downweighted.count = fits$robust.downweighted.count,
        robust.weight.min = fits$robust.weight.min,
        robust.weight.max = fits$robust.weight.max,
        boundary.score = boundary.score,
        finite.predictions = all(is.finite(averaged$fitted.values[target.index])),
        has.duplicate.rows = duplicate.info$has.duplicate.rows,
        n.duplicate.rows = duplicate.info$n.duplicate.rows,
        duplicate.groups = duplicate.info$duplicate.groups
    )
    out <- list(
        fitted.values = averaged$fitted.values,
        residuals = y - averaged$fitted.values,
        X = X,
        y = y,
        graph = graph.info,
        anchor.index = anchor.index,
        anchors = anchors,
        anchor.labels = anchor.labels,
        degree = params$degree,
        chart.dim = params$chart.dim,
        coordinate.method = coordinate.method,
        support.metric = support.metric,
        model.weight.rule = model.weight.rule,
        model.weights = model.weights,
        charts = charts,
        duplicate.action = duplicate.action,
        support.type = params$support.type,
        kernel = params$kernel,
        supports = supports,
        prediction.supports = prediction.supports,
        radii = diagnostics$radius,
        fit.weights = fits$fit.weights,
        robust.weights = fits$robust.weights,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant,
        case.weights = NULL,
        averaging.weights = averaged$averaging.weights,
        local.coefficients = fits$coefficients,
        design.columns = .malps.design.column.names(params$chart.dim,
                                                    params$degree),
        diagnostics = diagnostics,
        selection = selection,
        cv = if (identical(selection$method, "cv")) selection else NULL,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        call = call
    )
    class(out) <- c("malps", "list")
    out
}

.malps.prepare.fixed.params <- function(n, m, degree, support.type,
                                        support.size, radius, min.support,
                                        support.buffer, chart.dim,
                                        support.metric, model.weight.rule,
                                        kernel) {
    degree <- .malps.validate.degree(degree)
    required.support <- .malps.design.size(chart.dim, degree)
    if (is.null(min.support)) {
        min.support <- required.support + support.buffer
    }
    min.support <- .malps.validate.positive.integer(min.support, "min.support")
    if (min.support > n) {
        stop("min.support cannot exceed nrow(X).", call. = FALSE)
    }
    if (is.null(support.size)) {
        support.size <- min.support
    }
    if (identical(support.type, "knn")) {
        support.size <- .malps.validate.positive.integer(support.size,
                                                         "support.size")
        if (support.size > n) {
            stop("support.size cannot exceed nrow(X).", call. = FALSE)
        }
        if (support.size < min.support) {
            stop("support.size must be at least min.support.", call. = FALSE)
        }
    } else if (!is.null(support.size)) {
        support.size <- .malps.validate.positive.integer(support.size,
                                                         "support.size")
    } else {
        support.size <- NA_integer_
    }
    if (identical(support.type, "fixed.radius")) {
        radius <- .malps.validate.positive.scalar(radius, "radius")
    } else if (!is.null(radius)) {
        radius <- .malps.validate.positive.scalar(radius, "radius")
    } else {
        radius <- NA_real_
    }
    list(
        degree = degree,
        support.type = support.type,
        support.size = if (is.na(support.size)) NULL else support.size,
        min.support = min.support,
        radius = if (is.na(radius)) NULL else radius,
        chart.dim = chart.dim,
        support.metric = support.metric,
        model.weight.rule = model.weight.rule,
        kernel = kernel
    )
}

.malps.select.cv <- function(X, y, anchor.index, distance.matrix, support.type,
                             degree, degree.grid, support.size, support.grid,
                             radius, radius.grid, min.support,
                             min.support.grid, support.buffer, kernel,
                             kernel.grid, coordinate.method, chart.dim,
                             anchors, anchor.labels,
                             support.metric, model.weight.rule,
                             foldid, cv.folds, cv.loss,
                             cv.repeats, cv.seed, cv.one.se, local.solver,
                             normal.equations.max.condition,
                             robust.iterations, robust.tuning.constant) {
    n <- nrow(X)
    m <- ncol(X)
    folds <- .malps.prepare.folds(
        n = n,
        foldid = foldid,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        cv.seed = cv.seed
    )
    candidates <- .malps.make.candidate.table(
        n = n,
        m = chart.dim,
        degree = degree,
        degree.grid = degree.grid,
        support.type = support.type,
        support.size = support.size,
        support.grid = support.grid,
        radius = radius,
        radius.grid = radius.grid,
        min.support = min.support,
        min.support.grid = min.support.grid,
        support.buffer = support.buffer,
        kernel = kernel,
        kernel.grid = kernel.grid,
        support.metric = support.metric,
        model.weight.rule = model.weight.rule
    )
    total.folds <- length(folds$folds)
    candidate.rows <- vector("list", nrow(candidates))
    fold.rows <- list()
    for (i in seq_len(nrow(candidates))) {
        params <- .malps.candidate.params(candidates[i, , drop = FALSE])
        losses <- rep(NA_real_, total.folds)
        failures <- character(total.folds)
        fallback.rates <- rep(NA_real_, total.folds)
        median.condition <- rep(NA_real_, total.folds)
        median.positive.support <- rep(NA_real_, total.folds)
        median.coverage <- rep(NA_real_, total.folds)
        no.coverage <- rep(NA_integer_, total.folds)
        for (f in seq_along(folds$folds)) {
            test.index <- which(folds$folds[[f]]$foldid ==
                                  folds$folds[[f]]$fold)
            train.index <- setdiff(seq_len(n), test.index)
            fit <- tryCatch(
                .malps.fit.core(
                    X = X,
                    y = y,
                    anchor.index = anchor.index,
                    anchors = anchors,
                    anchor.labels = anchor.labels,
                    distance.matrix = distance.matrix,
                    params = params,
                    coordinate.method = coordinate.method,
                    support.metric = support.metric,
                    model.weight.rule = model.weight.rule,
                    graph.info = NULL,
                    duplicate.action = "keep",
                    duplicate.info = .malps.duplicate.row.info(X),
                    local.solver = local.solver,
                    normal.equations.max.condition =
                        normal.equations.max.condition,
                    robust.iterations = robust.iterations,
                    robust.tuning.constant = robust.tuning.constant,
                    selection = list(method = "cv-fold"),
                    call = NULL,
                    fit.index = train.index,
                    target.index = test.index,
                    fail.on.no.coverage = FALSE
                ),
                error = function(e) e
            )
            if (inherits(fit, "error")) {
                failures[f] <- conditionMessage(fit)
                next
            }
            no.coverage[f] <- fit$diagnostics$n.no.coverage
            if (fit$diagnostics$n.no.coverage > 0L) {
                failures[f] <- paste0(
                    "MALPS fold fit has ", fit$diagnostics$n.no.coverage,
                    " target(s) with no positive averaging coverage; ",
                    "uncovered.index = ",
                    .malps.format.index.message(
                        fit$diagnostics$uncovered.index
                    ), "."
                )
                next
            }
            pred <- fit$fitted.values[test.index]
            if (any(!is.finite(pred))) {
                failures[f] <- "MALPS fold fit produced non-finite predictions."
                next
            }
            losses[f] <- .malps.cv.loss(y[test.index], pred, cv.loss)
            failures[f] <- ""
            fallback.rates[f] <- mean(fit$diagnostics$fallback.used)
            median.condition[f] <- stats::median(
                fit$diagnostics$condition.number, na.rm = TRUE
            )
            median.positive.support[f] <- stats::median(
                fit$diagnostics$positive.weight.support.size, na.rm = TRUE
            )
            median.coverage[f] <- stats::median(
                fit$diagnostics$coverage.count[test.index], na.rm = TRUE
            )
        }
        valid <- is.finite(losses)
        loss.mean <- if (all(valid)) mean(losses) else Inf
        loss.se <- if (sum(valid) > 1L) {
            stats::sd(losses[valid]) / sqrt(sum(valid))
        } else {
            NA_real_
        }
        candidate.rows[[i]] <- data.frame(
            candidate.id = i,
            degree = params$degree,
            chart.dim = params$chart.dim,
            support.metric = params$support.metric,
            model.weight.rule = params$model.weight.rule,
            support.type = params$support.type,
            support.size = params$support.size %||% NA_integer_,
            min.support = params$min.support,
            radius = params$radius %||% NA_real_,
            kernel = params$kernel,
            robust.iterations = robust.iterations,
            robust.tuning.constant = robust.tuning.constant,
            loss.mean = loss.mean,
            loss.se = loss.se,
            n.valid.folds = sum(valid),
            n.valid.scores = sum(valid),
            n.fit.failures = sum(!valid),
            n.no.coverage = sum(no.coverage, na.rm = TRUE),
            fallback.rate = mean(fallback.rates, na.rm = TRUE),
            median.condition = stats::median(median.condition, na.rm = TRUE),
            median.positive.support = stats::median(median.positive.support,
                                                     na.rm = TRUE),
            median.coverage.count = stats::median(median.coverage,
                                                   na.rm = TRUE),
            stringsAsFactors = FALSE
        )
        fold.rows[[i]] <- data.frame(
            candidate.id = i,
            fold.instance = seq_along(losses),
            loss = losses,
            valid = valid,
            message = failures,
            stringsAsFactors = FALSE
        )
    }
    candidate.table <- do.call(rbind, candidate.rows)
    fold.table <- do.call(rbind, fold.rows)
    finite <- is.finite(candidate.table$loss.mean)
    if (!any(finite)) {
        stop("All MALPS CV candidates failed.", call. = FALSE)
    }
    selected.index <- .malps.select.cv.index(candidate.table, cv.one.se)
    selected.params <- .malps.candidate.params(candidates[selected.index, ,
                                                          drop = FALSE])
    list(
        method = "cv",
        candidate.table = candidate.table,
        fold.table = fold.table,
        selected.index = selected.index,
        selected.params = selected.params,
        loss.name = cv.loss,
        one.se.used = cv.one.se,
        fallback.method = NULL,
        n.valid = sum(finite),
        cv.folds = folds$cv.folds.effective,
        cv.folds.requested = folds$cv.folds.requested,
        cv.folds.effective = folds$cv.folds.effective,
            cv.repeats.requested = cv.repeats,
            cv.repeats.effective = folds$cv.repeats.effective,
        foldid.source = folds$foldid.source,
        seed = cv.seed,
        fold.message = folds$fold.message,
        foldid = folds$folds[[1L]]$foldid,
        foldid.list = lapply(folds$folds, `[[`, "foldid")
    )
}

.malps.select.gcv <- function(X, y, anchor.index, distance.matrix, support.type,
                              degree, degree.grid, support.size, support.grid,
                              radius, radius.grid, min.support,
                              min.support.grid, support.buffer, kernel,
                              kernel.grid, coordinate.method, chart.dim,
                              anchors, anchor.labels,
                              support.metric, model.weight.rule, graph.info,
                              duplicate.action, duplicate.info, local.solver,
                              normal.equations.max.condition,
                              robust.tuning.constant, gcv.exact.max.n) {
    n <- nrow(X)
    if (n > gcv.exact.max.n) {
        stop(sprintf(
            paste0("Exact dense MALPS GCV requested for n = %d, which exceeds ",
                   "gcv.exact.max.n = %d. Increase gcv.exact.max.n only for ",
                   "small reference runs, or use support.selection = \"cv\"."),
            n,
            gcv.exact.max.n
        ), call. = FALSE)
    }
    candidates <- .malps.make.candidate.table(
        n = n,
        m = chart.dim,
        degree = degree,
        degree.grid = degree.grid,
        support.type = support.type,
        support.size = support.size,
        support.grid = support.grid,
        radius = radius,
        radius.grid = radius.grid,
        min.support = min.support,
        min.support.grid = min.support.grid,
        support.buffer = support.buffer,
        kernel = kernel,
        kernel.grid = kernel.grid,
        support.metric = support.metric,
        model.weight.rule = model.weight.rule
    )
    candidate.rows <- vector("list", nrow(candidates))
    for (i in seq_len(nrow(candidates))) {
        params <- .malps.candidate.params(candidates[i, , drop = FALSE])
        fit <- tryCatch(
            .malps.fit.core(
                X = X,
                y = y,
                anchor.index = anchor.index,
                anchors = anchors,
                anchor.labels = anchor.labels,
                distance.matrix = distance.matrix,
                params = params,
                coordinate.method = coordinate.method,
                support.metric = support.metric,
                model.weight.rule = model.weight.rule,
                graph.info = graph.info,
                duplicate.action = duplicate.action,
                duplicate.info = duplicate.info,
                local.solver = local.solver,
                normal.equations.max.condition =
                    normal.equations.max.condition,
                robust.iterations = 0L,
                robust.tuning.constant = robust.tuning.constant,
                selection = list(method = "fixed"),
                call = NULL
            ),
            error = function(e) e
        )
        gcv.value <- edf <- mse <- loocv.mse <- NA_real_
        fallback.rate <- median.condition <- median.positive.support <-
            median.coverage <- NA_real_
        failure <- ""
        status <- "ok"
        if (inherits(fit, "error")) {
            status <- "fit_error"
            failure <- conditionMessage(fit)
        } else {
            gcv <- tryCatch(
                malps.gcv(
                    fit,
                    max.n = gcv.exact.max.n,
                    include.loocv = TRUE
                ),
                error = function(e) e
            )
            if (inherits(gcv, "error")) {
                status <- "gcv_error"
                failure <- conditionMessage(gcv)
            } else {
                gcv.value <- gcv$gcv
                edf <- gcv$edf
                mse <- gcv$mse
                loocv.mse <- gcv$loocv.mse %||% NA_real_
                if (!is.finite(gcv.value)) {
                    status <- "invalid_gcv"
                    failure <- "GCV was not finite."
                }
            }
            fallback.rate <- mean(fit$diagnostics$fallback.used)
            median.condition <- stats::median(
                fit$diagnostics$condition.number, na.rm = TRUE
            )
            median.positive.support <- stats::median(
                fit$diagnostics$positive.weight.support.size, na.rm = TRUE
            )
            median.coverage <- stats::median(
                fit$diagnostics$coverage.count, na.rm = TRUE
            )
        }
        candidate.rows[[i]] <- data.frame(
            candidate.id = i,
            degree = params$degree,
            chart.dim = params$chart.dim,
            support.metric = params$support.metric,
            model.weight.rule = params$model.weight.rule,
            support.type = params$support.type,
            support.size = params$support.size %||% NA_integer_,
            min.support = params$min.support,
            radius = params$radius %||% NA_real_,
            kernel = params$kernel,
            gcv = gcv.value,
            edf = edf,
            mse = mse,
            loocv.mse = loocv.mse,
            status = status,
            valid = identical(status, "ok") & is.finite(gcv.value),
            n.fit.failures = as.integer(!identical(status, "ok")),
            fallback.rate = fallback.rate,
            median.condition = median.condition,
            median.positive.support = median.positive.support,
            median.coverage.count = median.coverage,
            message = failure,
            stringsAsFactors = FALSE
        )
    }
    candidate.table <- do.call(rbind, candidate.rows)
    finite <- candidate.table$valid & is.finite(candidate.table$gcv)
    if (!any(finite)) {
        stop("All MALPS GCV candidates failed or produced non-finite GCV.",
             call. = FALSE)
    }
    selected.index <- candidate.table$candidate.id[finite][
        which.min(candidate.table$gcv[finite])
    ]
    selected.params <- .malps.candidate.params(candidates[selected.index, ,
                                                          drop = FALSE])
    list(
        method = "gcv",
        candidate.table = candidate.table,
        selected.index = selected.index,
        selected.params = selected.params,
        loss.name = "gcv",
        one.se.used = FALSE,
        fallback.method = NULL,
        n.valid = sum(finite),
        cv.folds = NULL,
        cv.repeats.requested = NULL,
        cv.repeats.effective = NULL,
        foldid.source = NULL,
        seed = NULL,
        fold.message = "Selected by exact dense MALPS GCV.",
        gcv.exact.max.n = gcv.exact.max.n
    )
}

.malps.validate.X <- function(X) {
    if (is.data.frame(X)) {
        X <- as.matrix(X)
    }
    if (!is.matrix(X) || !is.numeric(X)) {
        stop("X must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(X) < 1L || ncol(X) < 1L) {
        stop("X must have at least one row and one column.", call. = FALSE)
    }
    if (!all(is.finite(X))) {
        stop("X must contain only finite values.", call. = FALSE)
    }
    storage.mode(X) <- "double"
    X
}

.malps.validate.y <- function(y, n, name) {
    if (!is.numeric(y) || length(y) != n) {
        stop(name, " must be a numeric vector with length nrow(X).",
             call. = FALSE)
    }
    y <- as.numeric(y)
    if (!all(is.finite(y))) {
        stop(name, " must contain only finite values.", call. = FALSE)
    }
    y
}

.malps.validate.newdata <- function(newdata, m) {
    if (is.data.frame(newdata)) {
        newdata <- as.matrix(newdata)
    }
    if (is.null(dim(newdata)) && is.numeric(newdata) && m == 1L) {
        newdata <- matrix(newdata, ncol = 1L)
    }
    if (!is.matrix(newdata) || !is.numeric(newdata)) {
        stop("newdata must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(newdata) < 1L) {
        stop("newdata must have at least one row.", call. = FALSE)
    }
    if (ncol(newdata) != m) {
        stop("newdata must have the same number of columns as object$X.",
             call. = FALSE)
    }
    if (!all(is.finite(newdata))) {
        stop("newdata must contain only finite values.", call. = FALSE)
    }
    storage.mode(newdata) <- "double"
    newdata
}

.malps.prepare.graph.info <- function(n, graph, adj.list, weight.list,
                                      graph.stage) {
    has.graph.object <- !is.null(graph)
    has.payload <- !is.null(adj.list) || !is.null(weight.list)
    if (has.graph.object && has.payload) {
        stop("Supply either graph or adj.list/weight.list, not both.",
             call. = FALSE)
    }
    if (!has.graph.object && !has.payload) {
        return(NULL)
    }
    if (has.payload) {
        if (is.null(adj.list) || is.null(weight.list)) {
            stop("Both adj.list and weight.list are required for graph-geodesic MALPS supports.",
                 call. = FALSE)
        }
        validated <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
        if (length(validated$adj.list) != n) {
            stop("Graph vertex count must equal nrow(X).", call. = FALSE)
        }
        return(list(
            source = "supplied",
            stage = NULL,
            adj.list = validated$adj.list,
            weight.list = validated$weight.list
        ))
    }

    fields <- .graph.geodesic.fields(graph, stage = graph.stage)
    adj <- graph[[fields$adj]]
    weights <- graph[[fields$weight]]
    .validate.graph.geodesic.payload(adj, weights, fields)
    validated <- .validate.metric.graph.lowpass.graph(adj, weights)
    if (length(validated$adj.list) != n) {
        stop("Graph vertex count must equal nrow(X).", call. = FALSE)
    }
    list(
        source = paste0("graph.", fields$stage),
        stage = fields$stage,
        adj.list = validated$adj.list,
        weight.list = validated$weight.list
    )
}

.malps.resolve.support.metric <- function(support.metric, graph.input.supplied) {
    if (identical(support.metric, "auto")) {
        return(if (graph.input.supplied) "graph.geodesic" else "coordinates")
    }
    if (identical(support.metric, "graph.geodesic") && !graph.input.supplied) {
        stop(
            paste(
                "support.metric = 'graph.geodesic' requires graph input via",
                "graph or adj.list/weight.list."
            ),
            call. = FALSE
        )
    }
    support.metric
}

.malps.prepare.anchors <- function(X, anchor.index, anchor.coordinates,
                                   support.metric) {
    n <- nrow(X)
    m <- ncol(X)
    if (!is.null(anchor.index) && !is.null(anchor.coordinates)) {
        stop("Supply either anchor.index or anchor.coordinates, not both.",
             call. = FALSE)
    }
    if (!is.null(anchor.coordinates)) {
        anchors <- .malps.validate.anchor.coordinates(anchor.coordinates, m)
        if (identical(support.metric, "graph.geodesic")) {
            stop(
                paste(
                    "anchor.coordinates is not supported with",
                    "support.metric = 'graph.geodesic'; graph-geodesic",
                    "supports require observed graph vertices as anchors."
                ),
                call. = FALSE
            )
        }
        return(list(
            anchor.index = rep(NA_integer_, nrow(anchors)),
            anchors = anchors,
            anchor.labels = paste0("offsample.", seq_len(nrow(anchors))),
            anchor.type = "coordinates"
        ))
    }
    if (is.null(anchor.index)) {
        anchor.index <- seq_len(n)
    } else {
        anchor.index <- .malps.validate.anchor.index(anchor.index, n)
    }
    list(
        anchor.index = anchor.index,
        anchors = X[anchor.index, , drop = FALSE],
        anchor.labels = as.character(anchor.index),
        anchor.type = "observed"
    )
}

.malps.validate.anchor.coordinates <- function(anchor.coordinates, m) {
    if (is.data.frame(anchor.coordinates)) {
        anchor.coordinates <- as.matrix(anchor.coordinates)
    }
    if (is.null(dim(anchor.coordinates)) && is.numeric(anchor.coordinates) &&
        m == 1L) {
        anchor.coordinates <- matrix(anchor.coordinates, ncol = 1L)
    }
    if (!is.matrix(anchor.coordinates) || !is.numeric(anchor.coordinates)) {
        stop("anchor.coordinates must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(anchor.coordinates) < 1L) {
        stop("anchor.coordinates must have at least one row.", call. = FALSE)
    }
    if (ncol(anchor.coordinates) != m) {
        stop("anchor.coordinates must have the same number of columns as X.",
             call. = FALSE)
    }
    if (!all(is.finite(anchor.coordinates))) {
        stop("anchor.coordinates must contain only finite values.",
             call. = FALSE)
    }
    storage.mode(anchor.coordinates) <- "double"
    anchor.coordinates
}

.malps.support.distance.matrix <- function(X, anchors, anchor.index,
                                           support.metric,
                                           graph.info) {
    if (identical(support.metric, "coordinates")) {
        return(list(
            distance.matrix = .malps.distance.matrix(
                X, anchors
            ),
            metric = "coordinates",
            graph.finite = NA
        ))
    }
    if (is.null(graph.info)) {
        stop("graph.info is required for graph-geodesic support distances.",
             call. = FALSE)
    }
    if (any(is.na(anchor.index))) {
        stop("graph-geodesic support distances require observed anchor.index values.",
             call. = FALSE)
    }
    all.dist <- shortest.path(
        graph.info$adj.list,
        graph.info$weight.list,
        seq_along(graph.info$adj.list)
    )
    if (!is.matrix(all.dist) || any(dim(all.dist) != c(nrow(X), nrow(X)))) {
        stop("Graph geodesic distance computation returned an unexpected shape.",
             call. = FALSE)
    }
    if (any(!is.finite(all.dist))) {
        stop(
            paste(
                "MALPS graph-geodesic supports require finite shortest-path",
                "distances among all vertices; repair or connect the graph first."
            ),
            call. = FALSE
        )
    }
    list(
        distance.matrix = all.dist[, anchor.index, drop = FALSE],
        metric = "graph.geodesic",
        graph.finite = TRUE
    )
}

.malps.validate.chart.dim <- function(chart.dim, m, coordinate.method) {
    if (identical(coordinate.method, "coordinates")) {
        if (!is.null(chart.dim)) {
            stop("chart.dim is only used when coordinate.method = 'local.pca'.",
                 call. = FALSE)
        }
        return(m)
    }
    if (is.null(chart.dim)) {
        return(m)
    }
    chart.dim <- .malps.validate.positive.integer(chart.dim, "chart.dim")
    if (chart.dim > m) {
        stop("chart.dim cannot exceed ncol(X).", call. = FALSE)
    }
    chart.dim
}

.malps.validate.refit.weights <- function(weights, n) {
    if (is.null(weights)) {
        return(NULL)
    }
    if (!is.numeric(weights) || length(weights) != n) {
        stop("weights must be NULL or a numeric vector with length nrow(object$X).",
             call. = FALSE)
    }
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) {
        stop("weights must contain only finite values.", call. = FALSE)
    }
    if (any(weights < 0)) {
        stop("weights must be nonnegative.", call. = FALSE)
    }
    if (!any(weights > 0)) {
        stop("weights must contain at least one positive value.", call. = FALSE)
    }
    weights
}

.malps.validate.degree <- function(degree) {
    if (!is.numeric(degree) || length(degree) != 1L ||
        is.na(degree) || degree != as.integer(degree) ||
        !degree %in% 0:2) {
        stop("degree must be one of 0L, 1L, or 2L.", call. = FALSE)
    }
    as.integer(degree)
}

.malps.validate.anchor.index <- function(anchor.index, n) {
    if (!is.numeric(anchor.index) || any(is.na(anchor.index)) ||
        any(anchor.index != as.integer(anchor.index))) {
        stop("anchor.index must be an integer vector of row indices.",
             call. = FALSE)
    }
    anchor.index <- as.integer(anchor.index)
    if (!length(anchor.index)) {
        stop("anchor.index must contain at least one anchor.", call. = FALSE)
    }
    if (any(anchor.index < 1L | anchor.index > n)) {
        stop("anchor.index values must be between 1 and nrow(X).",
             call. = FALSE)
    }
    if (anyDuplicated(anchor.index)) {
        stop("anchor.index must not contain duplicates.", call. = FALSE)
    }
    anchor.index
}

.malps.validate.positive.integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        x != as.integer(x) || x < 1L) {
        stop(name, " must be a positive integer.", call. = FALSE)
    }
    as.integer(x)
}

.malps.validate.nonnegative.integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        x != as.integer(x) || x < 0L) {
        stop(name, " must be a nonnegative integer.", call. = FALSE)
    }
    as.integer(x)
}

.malps.validate.integer.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        x != as.integer(x)) {
        stop(name, " must be an integer scalar.", call. = FALSE)
    }
    as.integer(x)
}

.malps.validate.positive.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0) {
        stop(name, " must be a positive finite scalar.", call. = FALSE)
    }
    as.numeric(x)
}

.malps.validate.proportion <- function(x, name, strict = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
    if (ok) {
        ok <- if (strict) x > 0 && x < 1 else x >= 0 && x <= 1
    }
    if (!ok) {
        bound <- if (strict) "strictly between 0 and 1" else "between 0 and 1"
        stop(name, " must be ", bound, ".", call. = FALSE)
    }
    as.numeric(x)
}

.malps.validate.degree.grid <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.numeric(x) || !length(x) || any(is.na(x)) ||
        any(x != as.integer(x)) || any(!x %in% 0:2)) {
        stop(name, " must contain only 0L, 1L, and 2L.", call. = FALSE)
    }
    sort(unique(as.integer(x)))
}

.malps.validate.positive.integer.grid <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.numeric(x) || !length(x) || any(is.na(x)) ||
        any(x != as.integer(x)) || any(x < 1L)) {
        stop(name, " must contain positive integers.", call. = FALSE)
    }
    sort(unique(as.integer(x)))
}

.malps.validate.positive.numeric.grid <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.numeric(x) || !length(x) || any(is.na(x)) ||
        any(!is.finite(x)) || any(x <= 0)) {
        stop(name, " must contain positive finite values.", call. = FALSE)
    }
    sort(unique(as.numeric(x)))
}

.malps.validate.bootstrap.probs <- function(probs, n, weight.type) {
    if (is.null(probs)) {
        return(NULL)
    }
    if (!identical(weight.type, "multinomial")) {
        stop("probs is currently used only with weight.type = 'multinomial'.",
             call. = FALSE)
    }
    if (!is.numeric(probs) || length(probs) != n ||
        any(!is.finite(probs)) || any(probs < 0)) {
        stop("probs must be NULL or a nonnegative finite numeric vector with length nrow(object$X).",
             call. = FALSE)
    }
    total <- sum(probs)
    if (!is.finite(total) || total <= 0) {
        stop("probs must contain at least one positive value.", call. = FALSE)
    }
    as.numeric(probs / total)
}

.malps.bootstrap.weights <- function(n, weight.type, probs) {
    if (identical(weight.type, "bayesian")) {
        w <- stats::rexp(n)
        total <- sum(w)
        if (!is.finite(total) || total <= 0) {
            w[] <- 1
        } else {
            w <- w / mean(w)
        }
        return(as.numeric(w))
    }
    as.numeric(tabulate(
        sample.int(n, size = n, replace = TRUE, prob = probs),
        nbins = n
    ))
}

.malps.validate.kernel.grid <- function(x, default) {
    allowed <- c("epanechnikov", "triangular", "gaussian", "tricube")
    if (is.null(x)) {
        x <- default
    }
    x <- as.character(x)
    if (!length(x) || any(!x %in% allowed)) {
        stop("kernel.grid must contain supported kernel names.", call. = FALSE)
    }
    unique(x)
}

.malps.design.size <- function(m, degree) {
    choose(m + degree, degree)
}

.malps.prepare.folds <- function(n, foldid, cv.folds, cv.repeats, cv.seed) {
    if (!is.null(foldid)) {
        if (!is.numeric(foldid) || length(foldid) != n || any(is.na(foldid)) ||
            any(foldid != as.integer(foldid)) || any(foldid < 1L)) {
            stop("foldid must be a positive integer vector with length nrow(X).",
                 call. = FALSE)
        }
        foldid <- as.integer(foldid)
        folds.present <- sort(unique(foldid))
        fold.list <- lapply(folds.present, function(f) {
            list(repeat.index = 1L, fold = f, foldid = foldid)
        })
        msg <- if (cv.repeats != 1L) {
            "foldid supplied; cv.repeats forced to 1."
        } else {
            "foldid supplied."
        }
        return(list(
            folds = fold.list,
            foldid.source = "supplied",
            cv.folds.requested = cv.folds,
            cv.folds.effective = length(folds.present),
            cv.repeats.effective = 1L,
            fold.message = msg
        ))
    }
    if (cv.folds > n) {
        stop("cv.folds cannot exceed nrow(X) when foldid is not supplied.",
             call. = FALSE)
    }
    if (is.null(cv.seed)) {
        warning(
            "cv.seed is NULL; generated MALPS CV folds are not reproducible.",
            call. = FALSE
        )
    }
    old.seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv)
    } else {
        NULL
    }
    if (!is.null(cv.seed)) {
        set.seed(cv.seed)
    }
    on.exit({
        if (!is.null(cv.seed)) {
            if (is.null(old.seed)) {
                rm(".Random.seed", envir = .GlobalEnv)
            } else {
                assign(".Random.seed", old.seed, envir = .GlobalEnv)
            }
        }
    }, add = TRUE)
    fold.list <- list()
    k <- 0L
    for (r in seq_len(cv.repeats)) {
        ids <- sample(rep(seq_len(cv.folds), length.out = n))
        for (f in seq_len(cv.folds)) {
            k <- k + 1L
            fold.list[[k]] <- list(repeat.index = r, fold = f, foldid = ids)
        }
    }
    list(
        folds = fold.list,
        foldid.source = if (is.null(cv.seed)) "generated_unseeded" else "generated_seeded",
        cv.folds.requested = cv.folds,
        cv.folds.effective = cv.folds,
        cv.repeats.effective = cv.repeats,
        fold.message = if (is.null(cv.seed)) {
            "foldid generated without cv.seed; CV is not reproducible."
        } else {
            "foldid generated from cv.seed."
        }
    )
}

.malps.make.candidate.table <- function(n, m, degree, degree.grid,
                                        support.type, support.size,
                                        support.grid, radius, radius.grid,
                                        min.support, min.support.grid,
                                        support.buffer, kernel, kernel.grid,
                                        support.metric, model.weight.rule) {
    degree.grid <- .malps.validate.degree.grid(degree.grid, "degree.grid")
    if (is.null(degree.grid)) {
        degree.grid <- degree
    }
    support.grid <- .malps.validate.positive.integer.grid(support.grid,
                                                          "support.grid")
    min.support.grid <- .malps.validate.positive.integer.grid(
        min.support.grid, "min.support.grid"
    )
    radius.grid <- .malps.validate.positive.numeric.grid(radius.grid,
                                                         "radius.grid")
    kernel.grid <- .malps.validate.kernel.grid(kernel.grid, kernel)
    rows <- list()
    row.id <- 0L
    for (d in degree.grid) {
        required <- .malps.design.size(m, d)
        min.grid <- min.support.grid
        if (is.null(min.grid)) {
            min.grid <- if (is.null(min.support)) {
                .malps.default.support.grid(
                    required + .malps.support.buffer.safe(support.buffer),
                    n
                )
            } else {
                min.support
            }
        }
        min.grid <- sort(unique(as.integer(min.grid)))
        min.grid <- min.grid[min.grid >= required & min.grid <= n]
        if (!length(min.grid)) {
            next
        }
        for (ms in min.grid) {
            if (identical(support.type, "knn")) {
                sg <- support.grid
                if (is.null(sg)) {
                    sg <- if (is.null(support.size)) {
                        unique(pmin(n, ceiling(ms * c(1, 1.5, 2, 3))))
                    } else {
                        support.size
                    }
                }
                sg <- sg[sg >= ms & sg <= n]
                rg <- NA_real_
            } else if (identical(support.type, "adaptive.radius")) {
                sg <- NA_integer_
                rg <- if (is.null(radius.grid)) {
                    if (is.null(radius)) NA_real_ else radius
                } else {
                    radius.grid
                }
            } else {
                sg <- NA_integer_
                rg <- radius.grid
                if (is.null(rg)) {
                    rg <- radius
                }
                if (is.null(rg) || any(is.na(rg))) {
                    stop("radius or radius.grid is required for fixed.radius CV.",
                         call. = FALSE)
                }
            }
            for (ss in sg) {
                for (rr in rg) {
                    for (kk in kernel.grid) {
                        row.id <- row.id + 1L
                        rows[[row.id]] <- data.frame(
                            candidate.id = row.id,
                            degree = d,
                            chart.dim = m,
                            support.metric = support.metric,
                            model.weight.rule = model.weight.rule,
                            support.type = support.type,
                            support.size = ss,
                            min.support = ms,
                            radius = rr,
                            kernel = kk,
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    }
    if (!length(rows)) {
        stop("No valid MALPS CV candidates were generated.", call. = FALSE)
    }
    do.call(rbind, rows)
}

.malps.support.buffer.safe <- function(x) {
    .malps.validate.nonnegative.integer(x, "support.buffer")
}

.malps.default.support.grid <- function(base, n) {
    base <- max(1L, as.integer(ceiling(base)))
    sort(unique(as.integer(pmin(n, ceiling(base * c(1, 1.5, 2, 3, 4))))))
}

.malps.candidate.params <- function(row) {
    support.size <- row$support.size[[1L]]
    radius <- row$radius[[1L]]
    list(
        degree = as.integer(row$degree[[1L]]),
        chart.dim = as.integer(row$chart.dim[[1L]]),
        support.metric = as.character(row$support.metric[[1L]]),
        model.weight.rule = as.character(row$model.weight.rule[[1L]]),
        support.type = as.character(row$support.type[[1L]]),
        support.size = if (is.na(support.size)) NULL else as.integer(support.size),
        min.support = as.integer(row$min.support[[1L]]),
        radius = if (is.na(radius)) NULL else as.numeric(radius),
        kernel = as.character(row$kernel[[1L]])
    )
}

.malps.refit.with.supports <- function(object, y, weights, local.solver,
                                       normal.equations.max.condition,
                                       robust.iterations,
                                       robust.tuning.constant, call) {
    fits <- .malps.fit.local.models(
        X = object$X,
        y = y,
        anchor.index = object$anchor.index,
        supports = object$supports,
        charts = object$charts,
        degree = object$degree,
        observation.weights = weights,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        robust.iterations = robust.iterations,
        robust.tuning.constant = robust.tuning.constant
    )
    averaged <- .malps.average.predictions(
        X = object$X,
        anchor.labels = object$anchor.labels %||% object$anchor.index,
        supports = object$prediction.supports,
        charts = object$charts,
        coefficients = fits$coefficients,
        degree = object$degree
    )
    if (!all(is.finite(averaged$fitted.values))) {
        stop("MALPS refit produced non-finite fitted values.", call. = FALSE)
    }
    diagnostics <- object$diagnostics
    diagnostics$coverage.count <- averaged$coverage.count
    diagnostics$min.coverage.count <- min(averaged$coverage.count)
    diagnostics$n.no.coverage <- length(averaged$uncovered.index)
    diagnostics$uncovered.index <- averaged$uncovered.index
    diagnostics$rank <- fits$rank
    diagnostics$condition.number <- fits$condition.number
    diagnostics$rank.deficient <- fits$rank <
        .malps.design.size(object$chart.dim %||% ncol(object$X), object$degree)
    diagnostics$solver.used <- fits$solver.used
    diagnostics$fallback.used <- fits$fallback.used
    diagnostics$model.condition.number <- fits$model.condition.number
    diagnostics$effective.fit.support.size <- fits$effective.fit.support.size
    diagnostics$robust.iterations <- robust.iterations
    diagnostics$robust.tuning.constant <- robust.tuning.constant
    diagnostics$robust.used <- robust.iterations > 0L
    diagnostics$robust.downweighted.count <- fits$robust.downweighted.count
    diagnostics$robust.weight.min <- fits$robust.weight.min
    diagnostics$robust.weight.max <- fits$robust.weight.max
    diagnostics$finite.predictions <- all(is.finite(averaged$fitted.values))
    diagnostics$model.weight.fixed.from.original <- TRUE

    out <- object
    out$fitted.values <- averaged$fitted.values
    out$residuals <- y - averaged$fitted.values
    out$y <- y
    out$averaging.weights <- averaged$averaging.weights
    out$local.coefficients <- fits$coefficients
    out$fit.weights <- fits$fit.weights
    out$robust.weights <- fits$robust.weights
    out$robust.iterations <- robust.iterations
    out$robust.tuning.constant <- robust.tuning.constant
    out$case.weights <- weights
    out$diagnostics <- diagnostics
    out$local.solver <- local.solver
    out$normal.equations.max.condition <- normal.equations.max.condition
    out$call <- call
    class(out) <- c("malps", "list")
    out
}

.malps.predict.newdata <- function(object, newdata, allow.incomplete) {
    distance.matrix <- .malps.distance.matrix(newdata, object$anchors)
    prediction.supports <- .malps.build.prediction.supports(
        distance.matrix = distance.matrix,
        anchor.index = object$anchor.index,
        anchor.labels = object$anchor.labels %||% object$anchor.index,
        supports = object$supports,
        kernel = object$kernel,
        model.weights = object$model.weights %||% rep(1, length(object$supports)),
        target.index = seq_len(nrow(newdata))
    )
    averaged <- .malps.average.predictions(
        X = newdata,
        anchor.labels = object$anchor.labels %||% object$anchor.index,
        supports = prediction.supports,
        charts = object$charts,
        coefficients = object$local.coefficients,
        degree = object$degree
    )
    if (!allow.incomplete && length(averaged$uncovered.index)) {
        stop(sprintf(
            paste0("MALPS prediction has %d target(s) with no positive ",
                   "averaging coverage; uncovered.index = %s."),
            length(averaged$uncovered.index),
            .malps.format.index.message(averaged$uncovered.index)
        ), call. = FALSE)
    }
    averaged$fitted.values
}

.malps.cv.loss <- function(y, pred, loss) {
    resid <- y - pred
    if (identical(loss, "rmse")) {
        sqrt(mean(resid^2))
    } else if (identical(loss, "mae")) {
        mean(abs(resid))
    } else {
        mean(resid^2)
    }
}

.malps.select.cv.index <- function(candidate.table, cv.one.se) {
    finite <- is.finite(candidate.table$loss.mean)
    best <- which.min(candidate.table$loss.mean)
    if (!cv.one.se) {
        return(best)
    }
    threshold <- candidate.table$loss.mean[best] +
        ifelse(is.na(candidate.table$loss.se[best]), 0,
               candidate.table$loss.se[best])
    eligible <- candidate.table[finite & candidate.table$loss.mean <= threshold,
                                , drop = FALSE]
    ord <- order(
        eligible$degree,
        -eligible$min.support,
        -ifelse(is.na(eligible$support.size), eligible$min.support,
                eligible$support.size),
        -ifelse(is.na(eligible$radius), 0, eligible$radius),
        eligible$candidate.id
    )
    eligible$candidate.id[ord[1L]]
}

.malps.design.column.names <- function(m, degree) {
    names <- "1"
    if (degree >= 1L) {
        names <- c(names, paste0("z", seq_len(m)))
    }
    if (degree >= 2L) {
        for (i in seq_len(m)) {
            for (j in i:m) {
                names <- c(names, paste0("z", i, ".z", j))
            }
        }
    }
    names
}

.malps.design.matrix <- function(Z, degree) {
    Z <- as.matrix(Z)
    n <- nrow(Z)
    m <- ncol(Z)
    out <- matrix(1, nrow = n, ncol = 1L)
    if (degree >= 1L) {
        out <- cbind(out, Z)
    }
    if (degree >= 2L) {
        quad <- vector("list", m * (m + 1L) / 2L)
        k <- 0L
        for (i in seq_len(m)) {
            for (j in i:m) {
                k <- k + 1L
                quad[[k]] <- Z[, i] * Z[, j]
            }
        }
        out <- cbind(out, do.call(cbind, quad))
    }
    out
}

.malps.distance.matrix <- function(X, anchors) {
    n <- nrow(X)
    a <- nrow(anchors)
    out <- matrix(0, nrow = n, ncol = a)
    for (j in seq_len(a)) {
        diffs <- sweep(X, 2L, anchors[j, ], "-")
        out[, j] <- sqrt(rowSums(diffs^2))
    }
    out
}

.malps.duplicate.row.info <- function(X) {
    keys <- apply(X, 1L, function(row) {
        paste(format(row, digits = 17L), collapse = "\r")
    })
    groups <- split(seq_len(nrow(X)), keys)
    groups <- groups[vapply(groups, length, integer(1L)) > 1L]
    groups <- unname(lapply(groups, as.integer))
    n.duplicate.rows <- sum(vapply(groups, length, integer(1L)))
    list(
        has.duplicate.rows = length(groups) > 0L,
        n.duplicate.rows = as.integer(n.duplicate.rows),
        duplicate.groups = groups
    )
}

.malps.compact.kernel <- function(kernel) {
    kernel %in% c("epanechnikov", "triangular", "tricube")
}

.malps.kernel.weights <- function(distance, radius, kernel) {
    t <- distance / radius
    if (identical(kernel, "epanechnikov")) {
        w <- ifelse(t < 1, 1 - t^2, 0)
    } else if (identical(kernel, "triangular")) {
        w <- ifelse(t < 1, 1 - t, 0)
    } else if (identical(kernel, "tricube")) {
        w <- ifelse(t < 1, (1 - abs(t)^3)^3, 0)
    } else {
        w <- exp(-0.5 * t^2)
    }
    w[!is.finite(w)] <- 0
    as.numeric(w)
}

.malps.inflate.radius <- function(distance, kernel) {
    max.distance <- max(distance, 0)
    if (max.distance <= 0) {
        return(1)
    }
    if (.malps.compact.kernel(kernel)) {
        return(max.distance * (1 + sqrt(.Machine$double.eps)) +
                   .Machine$double.eps)
    }
    max.distance
}

.malps.build.supports <- function(distance.matrix, anchor.index, anchor.labels,
                                  support.type,
                                  support.size, radius, min.support, kernel,
                                  candidate.index = seq_len(nrow(distance.matrix))) {
    n.anchors <- length(anchor.labels)
    supports <- vector("list", n.anchors)
    n <- nrow(distance.matrix)
    candidate.index <- as.integer(candidate.index)
    if (length(candidate.index) < min.support) {
        stop("candidate.index has fewer rows than min.support.", call. = FALSE)
    }
    for (a in seq_len(n.anchors)) {
        d.full <- distance.matrix[, a]
        d <- d.full[candidate.index]
        ord.local <- order(d, candidate.index)
        rows.ordered <- candidate.index[ord.local]
        if (identical(support.type, "knn")) {
            if (support.size > length(candidate.index)) {
                stop("support.size cannot exceed the number of candidate rows.",
                     call. = FALSE)
            }
            anchor <- anchor.index[a]
            if (!is.na(anchor) && anchor %in% candidate.index) {
                idx <- c(anchor, rows.ordered[rows.ordered != anchor])
            } else {
                idx <- rows.ordered
            }
            idx <- idx[seq_len(support.size)]
            local.radius <- .malps.inflate.radius(d.full[idx], kernel)
        } else if (identical(support.type, "adaptive.radius")) {
            kth <- d[ord.local[min.support]]
            base.radius <- if (is.null(radius)) 0 else radius
            local.radius <- max(base.radius, kth)
            if (.malps.compact.kernel(kernel)) {
                local.radius <- local.radius * (1 + sqrt(.Machine$double.eps)) +
                    .Machine$double.eps
            }
            idx <- candidate.index[d <= local.radius]
        } else {
            local.radius <- radius
            idx <- candidate.index[d <= local.radius]
        }
        weights <- .malps.kernel.weights(d.full[idx], local.radius, kernel)
        keep <- weights > 0
        idx <- idx[keep]
        weights <- weights[keep]
        distances <- d.full[idx]
        if (length(idx) < min.support) {
            stop(sprintf(
                paste0("Anchor %s has only %d positive-weight support ",
                       "point(s), below min.support = %d."),
                anchor.labels[a], length(idx), min.support
            ), call. = FALSE)
        }
        supports[[a]] <- list(
            anchor = anchor.index[a],
            anchor.label = anchor.labels[a],
            index = as.integer(idx),
            distance = as.numeric(distances),
            radius = as.numeric(local.radius),
            weights = as.numeric(weights)
        )
    }
    supports
}

.malps.build.prediction.supports <- function(distance.matrix, anchor.index,
                                             anchor.labels,
                                             supports, kernel,
                                             model.weights = NULL,
                                             target.index = seq_len(nrow(distance.matrix))) {
    out <- vector("list", length(supports))
    target.index <- as.integer(target.index)
    if (is.null(model.weights)) {
        model.weights <- rep(1, length(supports))
    }
    if (!is.numeric(model.weights) || length(model.weights) != length(supports) ||
        any(!is.finite(model.weights)) || any(model.weights < 0)) {
        stop("model.weights must be a nonnegative finite numeric vector with one value per anchor.",
             call. = FALSE)
    }
    for (a in seq_along(supports)) {
        radius <- supports[[a]]$radius
        d <- distance.matrix[target.index, a]
        weights <- .malps.kernel.weights(d, radius, kernel) * model.weights[a]
        keep <- weights > 0
        idx <- target.index[keep]
        out[[a]] <- list(
            anchor = anchor.index[a],
            anchor.label = anchor.labels[a],
            index = as.integer(idx),
            distance = as.numeric(distance.matrix[idx, a]),
            radius = as.numeric(radius),
            weights = as.numeric(weights[keep])
        )
    }
    out
}

.malps.build.charts <- function(X, anchor.index, anchors, anchor.labels,
                                supports, coordinate.method,
                                chart.dim) {
    n.anchors <- length(supports)
    charts <- vector("list", n.anchors)
    m <- ncol(X)
    for (a in seq_len(n.anchors)) {
        center <- anchors[a, ]
        if (identical(coordinate.method, "coordinates")) {
            basis <- diag(m)
        } else {
            Z <- sweep(X[supports[[a]]$index, , drop = FALSE], 2L, center, "-")
            Xw <- Z * sqrt(supports[[a]]$weights)
            sv <- svd(Xw, nu = 0L, nv = chart.dim)
            basis <- sv$v[, seq_len(chart.dim), drop = FALSE]
            basis <- .malps.orient.pca.basis(basis)
        }
        charts[[a]] <- list(
            anchor = anchor.index[a],
            anchor.label = anchor.labels[a],
            center = as.numeric(center),
            basis = basis,
            coordinate.method = coordinate.method,
            chart.dim = ncol(basis)
        )
    }
    charts
}

.malps.orient.pca.basis <- function(basis) {
    for (j in seq_len(ncol(basis))) {
        i <- which.max(abs(basis[, j]))
        if (length(i) && basis[i, j] < 0) {
            basis[, j] <- -basis[, j]
        }
    }
    basis
}

.malps.chart.coordinates <- function(X, chart) {
    Z <- sweep(X, 2L, chart$center, "-")
    Z %*% chart$basis
}

.malps.fit.local.models <- function(X, y, anchor.index, supports, degree,
                                    charts, observation.weights = NULL, local.solver,
                                    normal.equations.max.condition,
                                    robust.iterations,
                                    robust.tuning.constant) {
    n.anchors <- length(supports)
    chart.dim <- charts[[1L]]$chart.dim
    p <- .malps.design.size(chart.dim, degree)
    if (is.null(observation.weights)) {
        observation.weights <- rep(1, nrow(X))
    }
    coefficients <- matrix(NA_real_, nrow = n.anchors, ncol = p)
    rank <- integer(n.anchors)
    condition.number <- rep(NA_real_, n.anchors)
    model.condition.number <- rep(NA_real_, n.anchors)
    solver.used <- character(n.anchors)
    fallback.used <- logical(n.anchors)
    effective.fit.support.size <- integer(n.anchors)
    fit.weights <- vector("list", n.anchors)
    robust.weights <- vector("list", n.anchors)
    robust.downweighted.count <- integer(n.anchors)
    robust.weight.min <- rep(NA_real_, n.anchors)
    robust.weight.max <- rep(NA_real_, n.anchors)
    for (a in seq_len(n.anchors)) {
        s <- supports[[a]]
        base.weights <- s$weights * observation.weights[s$index]
        positive <- base.weights > 0
        if (sum(positive) < 1L) {
            stop(sprintf(
                "MALPS local fit at anchor %s has no positive effective fit weights.",
                s$anchor.label %||% s$anchor
            ), call. = FALSE)
        }
        Z <- .malps.chart.coordinates(X[s$index, , drop = FALSE], charts[[a]])
        D <- .malps.design.matrix(Z, degree)
        y.local <- y[s$index]
        current.robust.weights <- rep(1, length(s$index))
        effective.weights <- base.weights
        solved <- .malps.solve.wls(
            D = D,
            y = y.local,
            weights = effective.weights,
            local.solver = local.solver,
            normal.equations.max.condition = normal.equations.max.condition
        )
        base.solved <- solved
        if (robust.iterations > 0L) {
            for (iter in seq_len(robust.iterations)) {
                residual <- y.local - as.numeric(D %*% solved$coefficients)
                next.robust.weights <- .malps.robust.bisquare.weights(
                    residual = residual,
                    tuning.constant = robust.tuning.constant,
                    eligible = base.weights > 0
                )
                if (is.null(next.robust.weights)) {
                    break
                }
                current.robust.weights <- next.robust.weights
                effective.weights <- base.weights * current.robust.weights
                if (!any(effective.weights > 0)) {
                    stop(sprintf(
                        paste0("MALPS robust local fit at anchor %s has no ",
                               "positive effective fit weights after residual ",
                               "reweighting."),
                        s$anchor.label %||% s$anchor
                    ), call. = FALSE)
                }
                solved <- .malps.solve.wls(
                    D = D,
                    y = y.local,
                    weights = effective.weights,
                    local.solver = local.solver,
                    normal.equations.max.condition =
                        normal.equations.max.condition
                )
            }
        }
        effective.fit.support.size[a] <- sum(effective.weights > 0)
        coefficients[a, ] <- solved$coefficients
        rank[a] <- solved$rank
        condition.number[a] <- solved$condition.number
        model.condition.number[a] <- base.solved$condition.number
        solver.used[a] <- solved$solver.used
        fallback.used[a] <- solved$fallback.used
        fit.weights[[a]] <- as.numeric(effective.weights)
        robust.weights[[a]] <- as.numeric(current.robust.weights)
        robust.downweighted.count[a] <- sum(
            current.robust.weights[positive] < 0.999
        )
        robust.weight.min[a] <- min(current.robust.weights[positive])
        robust.weight.max[a] <- max(current.robust.weights[positive])
    }
    colnames(coefficients) <- .malps.design.column.names(chart.dim, degree)
    list(
        coefficients = coefficients,
        rank = rank,
        condition.number = condition.number,
        model.condition.number = model.condition.number,
        solver.used = solver.used,
        fallback.used = fallback.used,
        effective.fit.support.size = effective.fit.support.size,
        fit.weights = fit.weights,
        robust.weights = robust.weights,
        robust.downweighted.count = robust.downweighted.count,
        robust.weight.min = robust.weight.min,
        robust.weight.max = robust.weight.max
    )
}

.malps.robust.bisquare.weights <- function(residual, tuning.constant,
                                           eligible = rep(TRUE, length(residual))) {
    eligible <- as.logical(eligible)
    if (!length(eligible) || length(eligible) != length(residual) ||
        !any(eligible)) {
        return(NULL)
    }
    scale <- stats::median(abs(residual[eligible]), na.rm = TRUE)
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
        return(NULL)
    }
    weights <- rep(1, length(residual))
    u <- residual[eligible] / (tuning.constant * scale)
    weights[eligible] <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
    weights[!is.finite(weights)] <- 0
    as.numeric(weights)
}

.malps.solve.wls <- function(D, y, weights, local.solver,
                             normal.equations.max.condition) {
    sqrt.w <- sqrt(weights)
    Xw <- D * sqrt.w
    yw <- y * sqrt.w
    p <- ncol(D)
    qr.obj <- qr(Xw)
    s <- svd(Xw, nu = 0L, nv = 0L)$d
    tol <- max(dim(Xw)) * max(s, 1) * .Machine$double.eps
    positive <- s[s > tol]
    rank <- length(positive)
    condition.number <- if (rank < p || length(positive) < 1L) {
        Inf
    } else {
        max(positive) / min(positive)
    }
    use.solver <- local.solver
    fallback <- FALSE
    if (identical(local.solver, "auto")) {
        if (rank == p && is.finite(condition.number) &&
            condition.number <= normal.equations.max.condition) {
            use.solver <- "normal.equations"
        } else {
            use.solver <- "svd"
            fallback <- TRUE
        }
    }
    if (identical(use.solver, "normal.equations")) {
        if (rank < p) {
            stop("normal.equations local solver encountered a rank-deficient design.",
                 call. = FALSE)
        }
        coef <- as.numeric(solve(crossprod(Xw), crossprod(Xw, yw)))
    } else if (identical(use.solver, "qr")) {
        coef <- as.numeric(qr.coef(qr.obj, yw))
        coef[is.na(coef)] <- 0
    } else {
        sv <- svd(Xw)
        tol <- max(dim(Xw)) * max(sv$d, 1) * .Machine$double.eps
        inv.d <- ifelse(sv$d > tol, 1 / sv$d, 0)
        coef <- as.numeric(sv$v %*% (inv.d * crossprod(sv$u, yw)))
    }
    list(
        coefficients = coef,
        rank = rank,
        condition.number = condition.number,
        solver.used = use.solver,
        fallback.used = fallback
    )
}

.malps.wls.coefficient.map <- function(D, weights) {
    if (!is.numeric(weights) || length(weights) != nrow(D) ||
        any(!is.finite(weights)) || any(weights < 0)) {
        stop("MALPS smoother matrix encountered invalid local fit weights.",
             call. = FALSE)
    }
    if (!any(weights > 0)) {
        stop("MALPS smoother matrix encountered a support with no positive fit weights.",
             call. = FALSE)
    }
    sqrt.w <- sqrt(weights)
    Xw <- D * sqrt.w
    sv <- svd(Xw)
    tol <- max(dim(Xw)) * max(sv$d, 1) * .Machine$double.eps
    inv.d <- ifelse(sv$d > tol, 1 / sv$d, 0)
    coef.map <- sv$v %*% (inv.d * t(sv$u))
    sweep(coef.map, 2L, sqrt.w, "*")
}

.malps.average.predictions <- function(X, supports, charts, coefficients, degree,
                                       anchor.labels = seq_len(length(supports))) {
    n <- nrow(X)
    n.anchors <- length(supports)
    numerator <- numeric(n)
    denominator <- numeric(n)
    coverage.count <- integer(n)
    averaging.weights <- matrix(0, nrow = n, ncol = n.anchors)
    for (a in seq_len(n.anchors)) {
        s <- supports[[a]]
        Z <- .malps.chart.coordinates(X[s$index, , drop = FALSE], charts[[a]])
        D <- .malps.design.matrix(Z, degree)
        pred <- as.numeric(D %*% coefficients[a, ])
        w <- s$weights
        numerator[s$index] <- numerator[s$index] + w * pred
        denominator[s$index] <- denominator[s$index] + w
        coverage.count[s$index] <- coverage.count[s$index] + as.integer(w > 0)
        averaging.weights[s$index, a] <- w
    }
    uncovered <- which(denominator <= 0 | !is.finite(denominator))
    fitted <- numerator / denominator
    if (length(uncovered)) {
        fitted[uncovered] <- NA_real_
    }
    positive <- denominator > 0 & is.finite(denominator)
    averaging.weights[positive, ] <- averaging.weights[positive, ,
                                                       drop = FALSE] /
        denominator[positive]
    colnames(averaging.weights) <- paste0("anchor.", anchor.labels)
    list(
        fitted.values = fitted,
        averaging.weights = averaging.weights,
        coverage.count = coverage.count,
        min.coverage.count = min(coverage.count),
        n.no.coverage = length(uncovered),
        uncovered.index = uncovered
    )
}

.malps.format.index.message <- function(index, max.show = 20L) {
    if (!length(index)) {
        return("integer(0)")
    }
    shown <- index[seq_len(min(length(index), max.show))]
    out <- paste(shown, collapse = ",")
    if (length(index) > max.show) {
        out <- paste0(out, ",...")
    }
    out
}

.malps.boundary.score <- function(X, anchor.index, anchors, supports) {
    score <- rep(NA_real_, length(supports))
    for (a in seq_along(supports)) {
        s <- supports[[a]]
        Z <- sweep(X[s$index, , drop = FALSE], 2L,
                   anchors[a, ], "-")
        wsum <- sum(s$weights)
        if (!is.finite(wsum) || wsum <= 0) {
            next
        }
        centroid <- colSums(Z * s$weights) / wsum
        numerator <- sqrt(sum(centroid^2))
        denom <- sqrt(sum(s$weights * rowSums(Z^2)) / wsum)
        if (is.finite(denom) && denom > 0) {
            score[a] <- numerator / denom
        }
    }
    score
}

.malps.model.quality.diagnostics <- function(rule, condition.number,
                                             positive.support.size,
                                             boundary.score) {
    n <- length(condition.number)
    if (!identical(length(positive.support.size), n) ||
        !identical(length(boundary.score), n)) {
        stop("MALPS model-quality diagnostics must have one value per anchor.",
             call. = FALSE)
    }

    condition.weight <- 1 / pmax(condition.number, 1)
    condition.weight[!is.finite(condition.weight) | condition.weight < 0] <- 0

    support.weight <- rep(1, n)
    max.support <- suppressWarnings(max(positive.support.size, na.rm = TRUE))
    if (is.finite(max.support) && max.support > 0) {
        support.weight <- positive.support.size / max.support
    }
    support.weight[!is.finite(support.weight) | support.weight < 0] <- 0

    boundary.weight <- 1 / (1 + pmax(boundary.score, 0))
    boundary.weight[!is.finite(boundary.weight) | boundary.weight < 0] <- 1

    weights <- switch(
        rule,
        none = rep(1, n),
        condition = condition.weight,
        support = support.weight,
        boundary = boundary.weight,
        quality = condition.weight * support.weight * boundary.weight,
        stop("Unknown model.weight.rule.", call. = FALSE)
    )
    weights[!is.finite(weights) | weights < 0] <- 0
    if (!any(weights > 0)) {
        weights[] <- 1
    }
    positive <- weights[weights > 0]
    scale <- stats::median(positive, na.rm = TRUE)
    if (is.finite(scale) && scale > 0) {
        weights <- weights / scale
    }
    floor <- if (identical(rule, "none")) {
        0
    } else {
        .Machine$double.eps
    }
    if (floor > 0) {
        weights[!is.finite(weights) | weights <= 0] <- floor
    }
    list(
        weight = as.numeric(weights),
        floor = floor,
        raw.weight = as.numeric(switch(
            rule,
            none = rep(1, n),
            condition = condition.weight,
            support = support.weight,
            boundary = boundary.weight,
            quality = condition.weight * support.weight * boundary.weight,
            stop("Unknown model.weight.rule.", call. = FALSE)
        )),
        condition.weight = as.numeric(condition.weight),
        support.weight = as.numeric(support.weight),
        boundary.weight = as.numeric(boundary.weight)
    )
}

.malps.model.quality.weights <- function(rule, condition.number,
                                         positive.support.size,
                                         boundary.score) {
    .malps.model.quality.diagnostics(
        rule = rule,
        condition.number = condition.number,
        positive.support.size = positive.support.size,
        boundary.score = boundary.score
    )$weight
}

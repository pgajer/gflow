#' Construct a Transported Graph Hessian Operator
#'
#' Builds an experimental transported-Hessian graph operator without fitting a
#' response. This is the phase-0 through phase-2 diagnostic layer for the
#' transported graph Hessian trend-filtering project: it constructs directed
#' edge differences, matches directions across a base dart, assembles the
#' second-difference operator, and returns diagnostics that make the transport
#' rule auditable.
#'
#' For a dart \eqn{u\to v}, the first difference is
#' \deqn{
#'   \delta_{u\to v}f = \frac{f(v)-f(u)}{\ell_{uv}}.
#' }
#' For a base dart \eqn{u\to u'} and a matched outgoing direction
#' \eqn{u\to v \leftrightarrow u'\to v'}, the transported Hessian row is
#' \deqn{
#'   \delta_{u'\to v'}f - \delta_{u\to v}.
#' }
#'
#' The reference rule, \code{transport.rule = "exact.coordinate"}, matches
#' directions with identical coordinate labels. Labels may be supplied directly
#' through \code{direction.labels}, or inferred from axis-aligned
#' \code{coordinates}. The first general metric rule,
#' \code{transport.rule = "local.embedding.soft"}, uses supplied coordinates to
#' softly match directions by angle and edge length. The edge-angle rules,
#' \code{"edge.angle.hard"} and \code{"edge.angle.soft"}, match directions by
#' comparing their angle relative to the transported base edge.
#' \code{transport.rule = "regression.gradient"} estimates local
#' least-squares gradients and compares gradient components across graph darts.
#' Those gradients can be computed either from supplied coordinates or from a
#' graph-derived local embedding built in a shared disk around each base dart.
#'
#' @param adj.list List of integer neighbor vectors using 1-based vertex
#'   indices. The graph must be undirected.
#' @param weight.list Optional list of positive edge lengths parallel to
#'   \code{adj.list}. If \code{NULL}, all edge lengths are one.
#' @param transport.rule Character scalar. \code{"exact.coordinate"} uses exact
#'   direction labels. \code{"local.embedding.soft"} uses coordinate direction
#'   vectors and softmax transport weights. \code{"edge.angle.hard"} and
#'   \code{"edge.angle.soft"} use coordinate-derived edge angles relative to the
#'   base dart. \code{"regression.gradient"} compares locally estimated gradient
#'   components.
#' @param coordinates Optional numeric matrix with one row per vertex. For
#'   \code{transport.rule = "exact.coordinate"}, coordinates must define
#'   axis-aligned graph edges so direction labels can be inferred. Coordinates
#'   are required by \code{"edge.angle.hard"} and \code{"edge.angle.soft"}.
#' @param direction.labels Optional list parallel to \code{adj.list}; each entry
#'   gives the exact direction label for the corresponding outgoing dart. When
#'   supplied, these labels override labels inferred from \code{coordinates}.
#' @param polynomial.probes Optional numeric vector or matrix with one row per
#'   vertex. If supplied, residual diagnostics
#'   \eqn{\|AP\|_F/\|P\|_F} are reported for the transported Hessian matrix
#'   \eqn{A}.
#' @param local.embedding.method Character scalar controlling the coordinates
#'   used by \code{transport.rule = "local.embedding.soft"}. \code{"auto"}
#'   uses supplied \code{coordinates} when available and otherwise tries
#'   \code{"grip.edge.kk"}.
#'   \code{"coordinates"} uses supplied coordinates directly.
#'   \code{"grip.edge.kk"} builds a local graph disk and prefers weighted GRIP
#'   followed by true edge-KK optimization when the installed \pkg{grip} exposes
#'   \code{grip.optimize.edge.kk.layout()} or its edge-isometric alias.
#'   \code{"mds.edge.kk"} uses classical MDS followed by the same edge-KK
#'   optimizer when available. \code{"cmdscale"} uses
#'   classical MDS only.
#' @param local.embedding.dim Positive integer embedding dimension for local
#'   graph-disk embeddings.
#' @param local.disk.hops Non-negative integer hop radius around each base dart
#'   used by graph-derived local embedding methods.
#' @param local.max.vertices Positive integer cap on local disk size. If a disk
#'   is larger, the nearest vertices by hop distance are retained.
#' @param return.sparse Logical. If \code{TRUE}, attach \pkg{Matrix} sparse
#'   matrices to the returned payloads.
#' @param tol Positive numeric tolerance for coordinate-axis direction
#'   inference.
#' @param soft.angle.scale,soft.length.scale Optional positive numeric scales
#'   for the soft transport angular and length score components. If \code{NULL},
#'   robust global defaults are estimated from candidate dart comparisons.
#' @param soft.bandwidth Positive numeric softmax bandwidth \eqn{\tau}. Smaller
#'   values make the soft rule closer to hard nearest-direction matching.
#' @param max.match.angle Optional non-negative numeric angle threshold in
#'   radians. For soft transport, candidate Hessian rows whose best match has
#'   angle larger than this threshold are dropped.
#' @param max.length.relative.error Optional non-negative numeric threshold for
#'   the best match's relative edge-length error,
#'   \eqn{|\ell_{\mathrm{matched}}-\ell_{\mathrm{direction}}|/
#'   \ell_{\mathrm{direction}}}. Candidate rows above the threshold are dropped.
#' @param min.match.margin Optional non-negative numeric threshold for the
#'   difference between the second-best and best soft-match scores. Candidate
#'   rows with smaller margins are dropped.
#' @param max.effective.matches Optional positive numeric threshold for
#'   \eqn{\exp(H)}, where \eqn{H} is the softmax transport entropy. Candidate
#'   rows with more diffuse matches are dropped.
#' @param match.threshold.rule Character scalar. \code{"none"} and
#'   \code{"fixed"} apply only the explicit fixed gates above.
#'   \code{"local.quantile"} also gates rows by local candidate-score and
#'   margin quantiles. \code{"local.robust.z"} also gates rows by robust
#'   local z-scores for best-match score and margin.
#' @param match.score.quantile,match.margin.quantile Quantile thresholds used
#'   by \code{match.threshold.rule = "local.quantile"}. The best score must
#'   lie at or below \code{match.score.quantile}; the best-vs-second margin
#'   must lie at or above \code{match.margin.quantile}.
#' @param min.best.score.z,min.margin.z Non-negative robust-z thresholds used
#'   by \code{match.threshold.rule = "local.robust.z"}. A larger
#'   best-score z-score means the best candidate is unusually low-scoring
#'   relative to local alternatives. A larger margin z-score means the
#'   best-vs-second separation is unusually large.
#' @param max.effective.match.fraction Optional positive numeric threshold for
#'   \eqn{\exp(H) / n_{\mathrm{candidate}}}. Candidate rows whose softmax
#'   mass is too diffuse across the local target directions are dropped.
#' @param edge.angle.scale,edge.length.scale Optional positive numeric scales
#'   for edge-angle matching. If \code{NULL}, robust global medians are
#'   estimated from candidate angle differences and length differences.
#' @param edge.angle.bandwidth Positive numeric softmax bandwidth used only by
#'   \code{transport.rule = "edge.angle.soft"}.
#' @param edge.angle.max.angle.difference Optional non-negative threshold, in
#'   radians, for accepting an edge-angle match.
#' @param edge.angle.max.length.relative.error Optional non-negative threshold
#'   for accepting an edge-angle match based on relative edge-length error.
#' @param gradient.coordinate.method Character scalar used only by
#'   \code{transport.rule = "regression.gradient"}. \code{"coordinates"} uses
#'   the supplied global \code{coordinates}. \code{"local.embedding"} estimates
#'   the two endpoint gradients of each base dart in one shared graph-derived
#'   local chart.
#' @param gradient.embedding.method Character scalar used by
#'   \code{gradient.coordinate.method = "local.embedding"}. The choices match
#'   the graph-derived local embedding backends used by
#'   \code{local.embedding.soft}.
#' @param gradient.embedding.dim Optional positive integer dimension for
#'   graph-derived regression-gradient local charts. If \code{NULL}, the
#'   supplied coordinate dimension is used when available, otherwise
#'   \code{local.embedding.dim}.
#' @param gradient.disk.hops Non-negative integer hop radius for graph-derived
#'   regression-gradient local charts.
#' @param gradient.max.vertices Positive integer cap on graph-derived
#'   regression-gradient local chart size.
#' @param gradient.disk.rule Character scalar controlling graph-derived
#'   regression-gradient disk construction. \code{"hops"} uses
#'   \code{gradient.disk.hops}. \code{"metric.diameter.fraction"} uses a
#'   two-center graph-geodesic disk whose metric radius is
#'   \code{gradient.disk.radius.fraction} times the graph diameter.
#'   \code{"metric.local.scale"} uses a two-center graph-geodesic disk whose
#'   metric radius is \code{gradient.disk.local.scale.multiplier} times the
#'   larger endpoint local scale.
#' @param gradient.disk.radius.fraction Positive numeric scalar used by
#'   \code{gradient.disk.rule = "metric.diameter.fraction"}.
#' @param gradient.disk.local.scale.method Character scalar controlling the
#'   per-vertex metric scale used by \code{gradient.disk.rule =
#'   "metric.local.scale"}. \code{"knn.distance"} uses the weighted graph
#'   distance to the \code{gradient.disk.local.scale.k}-th nearest vertex.
#'   \code{"median.incident.length"} and \code{"quantile.incident.length"} use
#'   incident edge-length summaries.
#' @param gradient.disk.local.scale.k Positive integer neighborhood index used
#'   by \code{gradient.disk.local.scale.method = "knn.distance"}.
#' @param gradient.disk.local.scale.quantile Numeric probability in
#'   \code{[0, 1]} used by \code{gradient.disk.local.scale.method =
#'   "quantile.incident.length"}.
#' @param gradient.disk.local.scale.multiplier Positive numeric scalar
#'   multiplier for \code{gradient.disk.rule = "metric.local.scale"}.
#' @param gradient.disk.min.vertices Non-negative integer lower target for
#'   metric disk coverage. Metric disks whose radius selects too few vertices
#'   are expanded to include at least this many nearest vertices, subject to
#'   \code{gradient.max.vertices}.
#' @param gradient.chart.selection Character scalar used only by graph-derived
#'   regression-gradient transport. \code{"fixed"} uses
#'   \code{gradient.embedding.method} and \code{gradient.disk.hops}.
#'   \code{"adaptive"} tries \code{gradient.embedding.candidates} crossed with
#'   \code{gradient.disk.hops.candidates} for each base dart and chooses the
#'   chart with the lowest graph-only diagnostic score.
#' @param gradient.embedding.candidates Character vector of graph-derived local
#'   embedding backends considered by \code{gradient.chart.selection =
#'   "adaptive"}.
#' @param gradient.disk.rule.candidates Character vector of disk construction
#'   rules considered by \code{gradient.chart.selection = "adaptive"}.
#' @param gradient.disk.hops.candidates Non-negative integer vector of hop
#'   radii considered by \code{gradient.chart.selection = "adaptive"}.
#' @param gradient.disk.radius.fraction.candidates Positive numeric vector of
#'   graph-diameter fractions considered by \code{gradient.chart.selection =
#'   "adaptive"}.
#' @param gradient.disk.local.scale.multiplier.candidates Positive numeric
#'   vector of local-scale multipliers considered by
#'   \code{gradient.chart.selection = "adaptive"}.
#' @param gradient.ridge Non-negative ridge parameter used by
#'   \code{transport.rule = "regression.gradient"} when inverting local
#'   weighted least-squares normal equations.
#'
#' @return A list of class \code{"transported.graph.hessian.operator"} with the
#'   canonical graph, directed-edge table, first-difference matrix, transported
#'   Hessian matrix, row metadata, dropped candidate rows, and diagnostics. For
#'   graph-derived regression-gradient transport with supplied
#'   \code{coordinates}, the embedding diagnostics also include synthetic
#'   chart-quality checks against those ambient coordinates.
#'
#' @export
transported.graph.hessian.operator <- function(adj.list,
                                               weight.list = NULL,
                                               transport.rule = c("exact.coordinate",
                                                                  "local.embedding.soft",
                                                                  "edge.angle.hard",
                                                                  "edge.angle.soft",
                                                                  "regression.gradient"),
                                               coordinates = NULL,
                                               direction.labels = NULL,
                                               polynomial.probes = NULL,
                                               local.embedding.method = c("auto",
                                                                          "coordinates",
                                                                          "grip.edge.kk",
                                                                          "mds.edge.kk",
                                                                          "cmdscale"),
                                               local.embedding.dim = 2L,
                                               local.disk.hops = 1L,
                                               local.max.vertices = 50L,
                                               return.sparse = TRUE,
                                               tol = 1e-8,
                                               soft.angle.scale = NULL,
                                               soft.length.scale = NULL,
                                               soft.bandwidth = 0.25,
                                               max.match.angle = NULL,
                                               max.length.relative.error = NULL,
                                               min.match.margin = NULL,
                                               max.effective.matches = NULL,
                                               match.threshold.rule = c("none",
                                                                        "fixed",
                                                                        "local.quantile",
                                                                        "local.robust.z"),
                                               match.score.quantile = 0.25,
                                               match.margin.quantile = 0.50,
                                               min.best.score.z = 1,
                                               min.margin.z = 0,
                                               max.effective.match.fraction = NULL,
                                               edge.angle.scale = NULL,
                                               edge.length.scale = NULL,
                                               edge.angle.bandwidth = 0.35,
                                               edge.angle.max.angle.difference = NULL,
                                               edge.angle.max.length.relative.error = NULL,
                                               gradient.coordinate.method = c("coordinates",
                                                                              "local.embedding"),
                                               gradient.embedding.method = c("grip.edge.kk",
                                                                             "mds.edge.kk",
                                                                             "cmdscale"),
                                               gradient.embedding.dim = NULL,
                                               gradient.disk.hops = local.disk.hops,
                                               gradient.max.vertices = local.max.vertices,
                                               gradient.disk.rule = c("hops",
                                                                      "metric.diameter.fraction",
                                                                      "metric.local.scale"),
                                               gradient.disk.radius.fraction = 0.10,
                                               gradient.disk.local.scale.method = c("knn.distance",
                                                                                    "median.incident.length",
                                                                                    "quantile.incident.length"),
                                               gradient.disk.local.scale.k = 8L,
                                               gradient.disk.local.scale.quantile = 0.75,
                                               gradient.disk.local.scale.multiplier = 2,
                                               gradient.disk.min.vertices = 0L,
                                               gradient.chart.selection = c("fixed",
                                                                            "adaptive"),
                                               gradient.embedding.candidates = c("cmdscale",
                                                                                 "mds.edge.kk"),
                                               gradient.disk.rule.candidates = c("hops",
                                                                                 "metric.diameter.fraction",
                                                                                 "metric.local.scale"),
                                               gradient.disk.hops.candidates = 1:5,
                                               gradient.disk.radius.fraction.candidates = c(0.05,
                                                                                            0.075,
                                                                                            0.10,
                                                                                            0.15,
                                                                                            0.20),
                                               gradient.disk.local.scale.multiplier.candidates = c(1,
                                                                                                   1.5,
                                                                                                   2,
                                                                                                   3,
                                                                                                   4),
                                               gradient.ridge = 1e-8) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for transported.graph.hessian.operator().",
             call. = FALSE)
    }
    transport.rule <- match.arg(transport.rule)
    local.embedding.method <- match.arg(local.embedding.method)
    gradient.coordinate.method <- match.arg(gradient.coordinate.method)
    gradient.embedding.method <- match.arg(gradient.embedding.method)
    gradient.disk.rule <- match.arg(gradient.disk.rule)
    gradient.disk.local.scale.method <- match.arg(gradient.disk.local.scale.method)
    gradient.chart.selection <- match.arg(gradient.chart.selection)
    gradient.embedding.candidates <- .transported.graph.hessian.validate.embedding.candidates(
        gradient.embedding.candidates,
        "gradient.embedding.candidates"
    )
    gradient.disk.rule.candidates <- .transported.graph.hessian.validate.disk.rule.candidates(
        gradient.disk.rule.candidates,
        "gradient.disk.rule.candidates"
    )
    gradient.disk.hops.candidates <- .transported.graph.hessian.validate.nonnegative.integer.vector(
        gradient.disk.hops.candidates,
        "gradient.disk.hops.candidates"
    )
    gradient.disk.radius.fraction <- .transported.graph.hessian.validate.positive.numeric.scalar(
        gradient.disk.radius.fraction,
        "gradient.disk.radius.fraction"
    )
    gradient.disk.radius.fraction.candidates <- .transported.graph.hessian.validate.positive.numeric.vector(
        gradient.disk.radius.fraction.candidates,
        "gradient.disk.radius.fraction.candidates"
    )
    gradient.disk.local.scale.k <- .validate.positive.integer.scalar(
        gradient.disk.local.scale.k,
        "gradient.disk.local.scale.k"
    )
    gradient.disk.local.scale.quantile <- .transported.graph.hessian.validate.probability.scalar(
        gradient.disk.local.scale.quantile,
        "gradient.disk.local.scale.quantile"
    )
    gradient.disk.local.scale.multiplier <- .transported.graph.hessian.validate.positive.numeric.scalar(
        gradient.disk.local.scale.multiplier,
        "gradient.disk.local.scale.multiplier"
    )
    gradient.disk.local.scale.multiplier.candidates <- .transported.graph.hessian.validate.positive.numeric.vector(
        gradient.disk.local.scale.multiplier.candidates,
        "gradient.disk.local.scale.multiplier.candidates"
    )
    if (!is.logical(return.sparse) || length(return.sparse) != 1L || is.na(return.sparse)) {
        stop("return.sparse must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("tol must be a finite positive numeric scalar.", call. = FALSE)
    }
    if (!is.numeric(gradient.ridge) || length(gradient.ridge) != 1L ||
        !is.finite(gradient.ridge) || gradient.ridge < 0) {
        stop("gradient.ridge must be a finite non-negative numeric scalar.",
             call. = FALSE)
    }
    gradient.disk.hops <- .transported.graph.hessian.validate.nonnegative.integer(
        gradient.disk.hops,
        "gradient.disk.hops"
    )
    gradient.max.vertices <- .validate.positive.integer.scalar(
        gradient.max.vertices,
        "gradient.max.vertices"
    )
    gradient.disk.min.vertices <- .transported.graph.hessian.validate.nonnegative.integer(
        gradient.disk.min.vertices,
        "gradient.disk.min.vertices"
    )
    match.threshold.rule <- match.arg(match.threshold.rule)
    soft.args <- .transported.graph.hessian.validate.soft.args(
        soft.angle.scale = soft.angle.scale,
        soft.length.scale = soft.length.scale,
        soft.bandwidth = soft.bandwidth,
        local.embedding.dim = local.embedding.dim,
        local.disk.hops = local.disk.hops,
        local.max.vertices = local.max.vertices,
        max.match.angle = max.match.angle,
        max.length.relative.error = max.length.relative.error,
        min.match.margin = min.match.margin,
        max.effective.matches = max.effective.matches,
        match.threshold.rule = match.threshold.rule,
        match.score.quantile = match.score.quantile,
        match.margin.quantile = match.margin.quantile,
        min.best.score.z = min.best.score.z,
        min.margin.z = min.margin.z,
        max.effective.match.fraction = max.effective.match.fraction
    )
    edge.angle.args <- .transported.graph.hessian.validate.edge.angle.args(
        edge.angle.scale = edge.angle.scale,
        edge.length.scale = edge.length.scale,
        edge.angle.bandwidth = edge.angle.bandwidth,
        edge.angle.max.angle.difference = edge.angle.max.angle.difference,
        edge.angle.max.length.relative.error = edge.angle.max.length.relative.error
    )

    graph <- .validate.graph.trend.filtering.graph(adj.list, weight.list)
    if (sum(lengths(graph$adj.list)) < 2L) {
        stop("The graph must contain at least one undirected edge.", call. = FALSE)
    }
    if (identical(local.embedding.method, "auto")) {
        local.embedding.method <- if (is.null(coordinates)) "grip.edge.kk" else "coordinates"
    }
    if (is.null(gradient.embedding.dim)) {
        gradient.embedding.dim <- if (!is.null(coordinates)) {
            ncol(as.matrix(coordinates))
        } else {
            soft.args$local.embedding.dim
        }
    }
    gradient.embedding.dim <- .validate.positive.integer.scalar(
        gradient.embedding.dim,
        "gradient.embedding.dim"
    )
    coordinates <- .transported.graph.hessian.coordinates(
        coordinates = coordinates,
        n.vertices = graph$n.vertices,
        required = (identical(transport.rule, "local.embedding.soft") &&
            identical(local.embedding.method, "coordinates")) ||
            transport.rule %in% c("edge.angle.hard", "edge.angle.soft") ||
            (identical(transport.rule, "regression.gradient") &&
             identical(gradient.coordinate.method, "coordinates"))
    )
    labels <- .transported.graph.hessian.direction.labels(
        graph = graph,
        coordinates = coordinates,
        direction.labels = direction.labels,
        require.axis.aligned = identical(transport.rule, "exact.coordinate"),
        tol = tol
    )
    darts <- .transported.graph.hessian.darts(graph, labels)
    first.difference <- .transported.graph.hessian.first.difference(darts, graph$n.vertices)
    transport <- switch(
        transport.rule,
        exact.coordinate = .transported.graph.hessian.exact.transport(darts),
        local.embedding.soft = .transported.graph.hessian.soft.transport(
            darts = darts,
            coordinates = coordinates,
            graph = graph,
            local.embedding.method = local.embedding.method,
            soft.args = soft.args
        ),
        edge.angle.hard = .transported.graph.hessian.edge.angle.transport(
            darts = darts,
            coordinates = coordinates,
            transport.rule = transport.rule,
            edge.angle.args = edge.angle.args
        ),
        edge.angle.soft = .transported.graph.hessian.edge.angle.transport(
            darts = darts,
            coordinates = coordinates,
            transport.rule = transport.rule,
            edge.angle.args = edge.angle.args
        ),
        regression.gradient = .transported.graph.hessian.regression.gradient.transport(
            darts = darts,
            graph = graph,
            coordinates = coordinates,
            ridge = gradient.ridge,
            gradient.args = list(
                coordinate.method = gradient.coordinate.method,
                embedding.method = gradient.embedding.method,
                embedding.dim = gradient.embedding.dim,
                disk.hops = gradient.disk.hops,
                max.vertices = gradient.max.vertices,
                disk.rule = gradient.disk.rule,
                disk.radius.fraction = gradient.disk.radius.fraction,
                disk.local.scale.method = gradient.disk.local.scale.method,
                disk.local.scale.k = gradient.disk.local.scale.k,
                disk.local.scale.quantile = gradient.disk.local.scale.quantile,
                disk.local.scale.multiplier = gradient.disk.local.scale.multiplier,
                disk.min.vertices = gradient.disk.min.vertices,
                chart.selection = gradient.chart.selection,
                embedding.candidates = gradient.embedding.candidates,
                disk.rule.candidates = gradient.disk.rule.candidates,
                disk.hops.candidates = gradient.disk.hops.candidates,
                disk.radius.fraction.candidates = gradient.disk.radius.fraction.candidates,
                disk.local.scale.multiplier.candidates =
                    gradient.disk.local.scale.multiplier.candidates,
                local.scales = .transported.graph.hessian.local.scales(
                    graph = graph,
                    method = gradient.disk.local.scale.method,
                    k = gradient.disk.local.scale.k,
                    prob = gradient.disk.local.scale.quantile
                ),
                graph.diameter = .transported.graph.hessian.graph.diameter(graph)
            )
        )
    )
    hessian <- .transported.graph.hessian.matrix(transport$row.table, graph$n.vertices)
    residuals <- .transported.graph.hessian.probe.residuals(
        hessian,
        polynomial.probes,
        graph$n.vertices
    )

    diagnostics <- list(
        summary = .transported.graph.hessian.summary(darts, transport),
        nullity.estimate = .estimate.graph.trend.filtering.nullity(hessian),
        polynomial.residuals = residuals
    )
    out <- list(
        graph = list(
            n.vertices = graph$n.vertices,
            n.darts = nrow(darts),
            adj.list = graph$adj.list,
            weight.list = graph$weight.list
        ),
        transport.rule = transport.rule,
        coordinates = coordinates,
        direction.labels = labels,
        darts = darts,
        first.difference = .graph.trend.filtering.matrix.payload(first.difference, return.sparse),
        hessian = .graph.trend.filtering.matrix.payload(hessian, return.sparse),
        penalty = .graph.trend.filtering.matrix.payload(hessian, return.sparse),
        transport = list(
            rule = transport.rule,
            row.table = transport$row.table,
            dropped.table = transport$dropped.table,
            summary = diagnostics$summary,
            scales = transport$scales,
            embedding.table = transport$embedding.table %||% data.frame(),
            gradient.diagnostics = transport$gradient.diagnostics %||% data.frame()
        ),
        row.table = transport$row.table,
        dropped.table = transport$dropped.table,
        diagnostics = diagnostics
    )
    out$penalty$kind <- "transported_hessian"
    out$penalty$order <- 1L
    class(out) <- c("transported.graph.hessian.operator", "list")
    out
}

.transported.graph.hessian.validate.soft.args <- function(soft.angle.scale,
                                                          soft.length.scale,
                                                          soft.bandwidth,
                                                          local.embedding.dim,
                                                          local.disk.hops,
                                                          local.max.vertices,
                                                          max.match.angle,
                                                          max.length.relative.error,
                                                          min.match.margin,
                                                          max.effective.matches,
                                                          match.threshold.rule,
                                                          match.score.quantile,
                                                          match.margin.quantile,
                                                          min.best.score.z,
                                                          min.margin.z,
                                                          max.effective.match.fraction) {
    for (nm in c("soft.angle.scale", "soft.length.scale")) {
        val <- get(nm)
        if (!is.null(val) &&
            (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0)) {
            stop(sprintf("%s must be NULL or a finite positive numeric scalar.", nm),
                 call. = FALSE)
        }
    }
    if (!is.numeric(soft.bandwidth) || length(soft.bandwidth) != 1L ||
        !is.finite(soft.bandwidth) || soft.bandwidth <= 0) {
        stop("soft.bandwidth must be a finite positive numeric scalar.", call. = FALSE)
    }
    .transported.graph.hessian.validate.optional.nonnegative(max.match.angle,
                                                             "max.match.angle")
    .transported.graph.hessian.validate.optional.nonnegative(max.length.relative.error,
                                                             "max.length.relative.error")
    .transported.graph.hessian.validate.optional.nonnegative(min.match.margin,
                                                             "min.match.margin")
    if (!is.null(max.effective.matches) &&
        (!is.numeric(max.effective.matches) || length(max.effective.matches) != 1L ||
         !is.finite(max.effective.matches) || max.effective.matches <= 0)) {
        stop("max.effective.matches must be NULL or a finite positive numeric scalar.",
             call. = FALSE)
    }
    for (nm in c("match.score.quantile", "match.margin.quantile")) {
        val <- get(nm)
        if (!is.numeric(val) || length(val) != 1L || !is.finite(val) ||
            val < 0 || val > 1) {
            stop(sprintf("%s must be a finite numeric scalar in [0, 1].", nm),
                 call. = FALSE)
        }
    }
    for (nm in c("min.best.score.z", "min.margin.z")) {
        val <- get(nm)
        if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val < 0) {
            stop(sprintf("%s must be a finite non-negative numeric scalar.", nm),
                 call. = FALSE)
        }
    }
    if (!is.null(max.effective.match.fraction) &&
        (!is.numeric(max.effective.match.fraction) ||
         length(max.effective.match.fraction) != 1L ||
         !is.finite(max.effective.match.fraction) ||
         max.effective.match.fraction <= 0)) {
        stop("max.effective.match.fraction must be NULL or a finite positive numeric scalar.",
             call. = FALSE)
    }
    local.embedding.dim <- .validate.positive.integer.scalar(local.embedding.dim,
                                                             "local.embedding.dim")
    local.disk.hops <- .transported.graph.hessian.validate.nonnegative.integer(
        local.disk.hops,
        "local.disk.hops"
    )
    local.max.vertices <- .validate.positive.integer.scalar(local.max.vertices,
                                                            "local.max.vertices")
    list(
        angle.scale = soft.angle.scale,
        length.scale = soft.length.scale,
        bandwidth = as.double(soft.bandwidth),
        local.embedding.dim = local.embedding.dim,
        local.disk.hops = local.disk.hops,
        local.max.vertices = local.max.vertices,
        max.match.angle = max.match.angle,
        max.length.relative.error = max.length.relative.error,
        min.match.margin = min.match.margin,
        max.effective.matches = max.effective.matches,
        match.threshold.rule = match.threshold.rule,
        match.score.quantile = as.double(match.score.quantile),
        match.margin.quantile = as.double(match.margin.quantile),
        min.best.score.z = as.double(min.best.score.z),
        min.margin.z = as.double(min.margin.z),
        max.effective.match.fraction = max.effective.match.fraction
    )
}

.transported.graph.hessian.validate.optional.nonnegative <- function(x, name) {
    if (!is.null(x) &&
        (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0)) {
        stop(sprintf("%s must be NULL or a finite non-negative numeric scalar.", name),
             call. = FALSE)
    }
    invisible(NULL)
}

.transported.graph.hessian.validate.edge.angle.args <- function(edge.angle.scale,
                                                                edge.length.scale,
                                                                edge.angle.bandwidth,
                                                                edge.angle.max.angle.difference,
                                                                edge.angle.max.length.relative.error) {
    for (nm in c("edge.angle.scale", "edge.length.scale")) {
        val <- get(nm)
        if (!is.null(val) &&
            (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0)) {
            stop(sprintf("%s must be NULL or a finite positive numeric scalar.", nm),
                 call. = FALSE)
        }
    }
    if (!is.numeric(edge.angle.bandwidth) || length(edge.angle.bandwidth) != 1L ||
        !is.finite(edge.angle.bandwidth) || edge.angle.bandwidth <= 0) {
        stop("edge.angle.bandwidth must be a finite positive numeric scalar.",
             call. = FALSE)
    }
    .transported.graph.hessian.validate.optional.nonnegative(
        edge.angle.max.angle.difference,
        "edge.angle.max.angle.difference"
    )
    .transported.graph.hessian.validate.optional.nonnegative(
        edge.angle.max.length.relative.error,
        "edge.angle.max.length.relative.error"
    )
    list(
        angle.scale = edge.angle.scale,
        length.scale = edge.length.scale,
        bandwidth = as.double(edge.angle.bandwidth),
        max.match.angle = edge.angle.max.angle.difference,
        max.length.relative.error = edge.angle.max.length.relative.error,
        min.match.margin = NULL,
        max.effective.matches = NULL,
        match.threshold.rule = "none",
        match.score.quantile = 0.25,
        match.margin.quantile = 0.50,
        min.best.score.z = 0,
        min.margin.z = 0,
        max.effective.match.fraction = NULL
    )
}

.transported.graph.hessian.validate.nonnegative.integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L ||
        is.na(x) || !is.finite(x) ||
        x != floor(x) || x < 0) {
        stop(sprintf("%s must be a non-negative integer scalar.", name),
             call. = FALSE)
    }
    as.integer(x)
}

.transported.graph.hessian.validate.nonnegative.integer.vector <- function(x, name) {
    if (!is.numeric(x) || !length(x) ||
        any(is.na(x)) || any(!is.finite(x)) ||
        any(x != floor(x)) || any(x < 0)) {
        stop(sprintf("%s must be a non-empty non-negative integer vector.", name),
             call. = FALSE)
    }
    unique(as.integer(x))
}

.transported.graph.hessian.validate.embedding.candidates <- function(x, name) {
    allowed <- c("grip.edge.kk", "mds.edge.kk", "cmdscale")
    if (!is.character(x) || !length(x) || any(is.na(x)) ||
        any(!nzchar(x)) || any(!x %in% allowed)) {
        stop(sprintf("%s must be a non-empty character vector with values from %s.",
                     name, paste(shQuote(allowed), collapse = ", ")),
             call. = FALSE)
    }
    unique(x)
}

.transported.graph.hessian.validate.disk.rule.candidates <- function(x, name) {
    allowed <- c("hops", "metric.diameter.fraction", "metric.local.scale")
    if (!is.character(x) || !length(x) || any(is.na(x)) ||
        any(!nzchar(x)) || any(!x %in% allowed)) {
        stop(sprintf("%s must be a non-empty character vector with values from %s.",
                     name, paste(shQuote(allowed), collapse = ", ")),
             call. = FALSE)
    }
    unique(x)
}

.transported.graph.hessian.validate.positive.numeric.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(sprintf("%s must be a finite positive numeric scalar.", name),
             call. = FALSE)
    }
    as.double(x)
}

.transported.graph.hessian.validate.probability.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || x > 1) {
        stop(sprintf("%s must be a finite numeric scalar in [0, 1].", name),
             call. = FALSE)
    }
    as.double(x)
}

.transported.graph.hessian.validate.positive.numeric.vector <- function(x, name) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x <= 0)) {
        stop(sprintf("%s must be a non-empty finite positive numeric vector.", name),
             call. = FALSE)
    }
    unique(as.double(x))
}

.transported.graph.hessian.coordinates <- function(coordinates, n.vertices,
                                                   required) {
    if (is.null(coordinates)) {
        if (isTRUE(required)) {
            stop("coordinates must be supplied for the selected transport rule.",
                 call. = FALSE)
        }
        return(NULL)
    }
    coordinates <- as.matrix(coordinates)
    storage.mode(coordinates) <- "double"
    if (nrow(coordinates) != n.vertices || ncol(coordinates) < 1L ||
        any(!is.finite(coordinates))) {
        stop("coordinates must be a finite numeric matrix with one row per vertex.",
             call. = FALSE)
    }
    coordinates
}

.transported.graph.hessian.direction.labels <- function(graph, coordinates,
                                                        direction.labels,
                                                        require.axis.aligned,
                                                        tol) {
    if (!is.null(direction.labels)) {
        if (!is.list(direction.labels) || length(direction.labels) != graph$n.vertices) {
            stop("direction.labels must be a list parallel to adj.list.", call. = FALSE)
        }
        out <- vector("list", graph$n.vertices)
        for (idx in seq_len(graph$n.vertices)) {
            labels <- direction.labels[[idx]]
            if (length(labels) != length(graph$adj.list[[idx]])) {
                stop("Each direction.labels entry must be parallel to the corresponding adjacency entry.",
                     call. = FALSE)
            }
            labels <- as.character(labels)
            if (anyNA(labels) || any(!nzchar(labels))) {
                stop("direction.labels entries must be non-empty character labels.",
                     call. = FALSE)
            }
            out[[idx]] <- labels
        }
        return(out)
    }

    if (is.null(coordinates)) {
        if (isTRUE(require.axis.aligned)) {
            stop("coordinates or direction.labels must be supplied for exact.coordinate transport.",
                 call. = FALSE)
        }
        return(lapply(graph$adj.list, function(nbrs) paste0("dir", seq_along(nbrs))))
    }

    out <- vector("list", graph$n.vertices)
    for (from in seq_len(graph$n.vertices)) {
        nbrs <- graph$adj.list[[from]]
        labels <- character(length(nbrs))
        for (kk in seq_along(nbrs)) {
            delta <- coordinates[nbrs[kk], , drop = TRUE] -
                coordinates[from, , drop = TRUE]
            active <- which(abs(delta) > tol)
            if (length(active) != 1L) {
                if (!isTRUE(require.axis.aligned)) {
                    labels[kk] <- paste0("dir", kk)
                    next
                }
                stop("exact.coordinate transport requires axis-aligned coordinates or explicit direction.labels.",
                     call. = FALSE)
            }
            labels[kk] <- sprintf("axis%d%s", active, if (delta[active] > 0) "+" else "-")
        }
        out[[from]] <- labels
    }
    out
}

.transported.graph.hessian.darts <- function(graph, labels) {
    dart <- integer()
    from <- integer()
    to <- integer()
    neighbor.index <- integer()
    length <- numeric()
    direction.label <- character()

    for (u in seq_len(graph$n.vertices)) {
        nbrs <- graph$adj.list[[u]]
        lens <- graph$weight.list[[u]]
        if (length(nbrs) < 1L) next
        ids <- seq_along(nbrs)
        dart <- c(dart, seq.int(length(dart) + 1L, length(dart) + length(nbrs)))
        from <- c(from, rep.int(u, length(nbrs)))
        to <- c(to, nbrs)
        neighbor.index <- c(neighbor.index, ids)
        length <- c(length, lens)
        direction.label <- c(direction.label, labels[[u]])
    }

    data.frame(
        dart = as.integer(dart),
        from = as.integer(from),
        to = as.integer(to),
        neighbor.index = as.integer(neighbor.index),
        length = as.double(length),
        direction.label = direction.label,
        stringsAsFactors = FALSE
    )
}

.transported.graph.hessian.first.difference <- function(darts, n.vertices) {
    ii <- rep(darts$dart, each = 2L)
    jj <- as.integer(rbind(darts$from, darts$to))
    xx <- as.double(rbind(-1 / darts$length, 1 / darts$length))
    Matrix::sparseMatrix(
        i = ii,
        j = jj,
        x = xx,
        dims = c(nrow(darts), n.vertices),
        giveCsparse = TRUE
    )
}

.transported.graph.hessian.exact.transport <- function(darts) {
    row.records <- list()
    dropped.records <- list()
    row <- 0L
    candidate <- 0L
    key <- paste(darts$from, darts$direction.label, sep = "::")
    by.key <- split(seq_len(nrow(darts)), key)

    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        local <- darts[darts$from == base$from, , drop = FALSE]
        for (local.idx in seq_len(nrow(local))) {
            candidate <- candidate + 1L
            direction <- local[local.idx, , drop = FALSE]
            match.key <- paste(base$to, direction$direction.label, sep = "::")
            matches <- by.key[[match.key]]
            status <- if (length(matches) < 1L) {
                "dropped_no_exact_match"
            } else if (length(matches) > 1L) {
                "matched_ambiguous"
            } else {
                "matched"
            }
            if (length(matches) < 1L) {
                dropped.records[[length(dropped.records) + 1L]] <- data.frame(
                    candidate = candidate,
                    base.dart = base$dart,
                    base.from = base$from,
                    base.to = base$to,
                    direction.dart = direction$dart,
                    direction.from = direction$from,
                    direction.to = direction$to,
                    direction.label = direction$direction.label,
                    match.status = status,
                    stringsAsFactors = FALSE
                )
                next
            }
            matched <- darts[matches[1L], , drop = FALSE]
            row <- row + 1L
            row.records[[length(row.records) + 1L]] <- data.frame(
                row = row,
                candidate = candidate,
                base.dart = base$dart,
                base.from = base$from,
                base.to = base$to,
                base.length = base$length,
                direction.dart = direction$dart,
                direction.from = direction$from,
                direction.to = direction$to,
                direction.length = direction$length,
                direction.label = direction$direction.label,
                matched.dart = matched$dart,
                matched.from = matched$from,
                matched.to = matched$to,
                matched.length = matched$length,
                matched.direction.label = matched$direction.label,
                transport.weight = 1,
                transport.entropy = 0,
                match.status = status,
                n.exact.matches = length(matches),
                stringsAsFactors = FALSE
            )
        }
    }

    row.table <- if (length(row.records)) {
        do.call(rbind, row.records)
    } else {
        data.frame()
    }
    dropped.table <- if (length(dropped.records)) {
        do.call(rbind, dropped.records)
    } else {
        data.frame()
    }
    list(
        row.table = row.table,
        dropped.table = dropped.table,
        n.candidates = candidate,
        scales = NULL
    )
}

.transported.graph.hessian.soft.transport <- function(darts, coordinates, graph,
                                                      local.embedding.method,
                                                      soft.args) {
    if (!identical(local.embedding.method, "coordinates")) {
        return(.transported.graph.hessian.soft.transport.local(
            darts = darts,
            graph = graph,
            local.embedding.method = local.embedding.method,
            soft.args = soft.args
        ))
    }
    dart.vectors <- coordinates[darts$to, , drop = FALSE] -
        coordinates[darts$from, , drop = FALSE]
    dart.norms <- sqrt(rowSums(dart.vectors^2))
    if (any(dart.norms <= 0)) {
        stop("coordinates must give nonzero direction vectors for every graph dart.",
             call. = FALSE)
    }
    scales <- .transported.graph.hessian.soft.scales(
        darts = darts,
        dart.vectors = dart.vectors,
        dart.norms = dart.norms,
        soft.args = soft.args
    )

    row.records <- list()
    dropped.records <- list()
    row <- 0L
    candidate <- 0L
    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        local.ids <- which(darts$from == base$from)
        target.ids <- which(darts$from == base$to)
        if (length(target.ids) < 1L) next
        for (direction.id in local.ids) {
            candidate <- candidate + 1L
            direction <- darts[direction.id, , drop = FALSE]
            scores <- .transported.graph.hessian.soft.scores(
                direction.id = direction.id,
                target.ids = target.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms,
                scales = scales
            )
            weights <- .transported.graph.hessian.softmax(scores$total, scales$bandwidth)
            reverse.scores <- .transported.graph.hessian.soft.scores(
                direction.id = target.ids[which.min(scores$total)],
                target.ids = local.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms,
                scales = scales
            )
            reciprocal.best <- local.ids[which.min(reverse.scores$total)] == direction.id
            row <- row + 1L
            entropy <- -sum(weights * log(pmax(weights, .Machine$double.xmin)))
            effective.matches <- exp(entropy)
            confidence <- .transported.graph.hessian.match.confidence(
                scores = scores,
                weights = weights,
                direction.length = direction$length,
                matched.lengths = darts$length[target.ids],
                entropy = entropy,
                effective.matches = effective.matches,
                reciprocal.best = reciprocal.best
            )
            gate <- .transported.graph.hessian.match.gate(confidence, soft.args)
            if (!gate$accepted) {
                row <- row - 1L
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = confidence,
                        drop.reason = gate$reason
                    )
                next
            }
            for (kk in seq_along(target.ids)) {
                matched <- darts[target.ids[kk], , drop = FALSE]
                row.records[[length(row.records) + 1L]] <- data.frame(
                    row = row,
                    candidate = candidate,
                    base.dart = base$dart,
                    base.from = base$from,
                    base.to = base$to,
                    base.length = base$length,
                    direction.dart = direction$dart,
                    direction.from = direction$from,
                    direction.to = direction$to,
                    direction.length = direction$length,
                    direction.label = direction$direction.label,
                    matched.dart = matched$dart,
                    matched.from = matched$from,
                    matched.to = matched$to,
                    matched.length = matched$length,
                    matched.direction.label = matched$direction.label,
                    transport.weight = weights[kk],
                    transport.entropy = entropy,
                    effective.matches = effective.matches,
                    effective.match.fraction = confidence$effective.match.fraction,
                    angle = scores$angle[kk],
                    length.difference = scores$length.difference[kk],
                    match.score = scores$total[kk],
                    is.best.match = kk == confidence$best.index,
                    n.match.candidates = confidence$n.match.candidates,
                    best.score.local.rank = confidence$best.score.local.rank,
                    best.score.local.quantile = confidence$best.score.local.quantile,
                    best.score.robust.z = confidence$best.score.robust.z,
                    margin.local.quantile = confidence$margin.local.quantile,
                    margin.robust.z = confidence$margin.robust.z,
                    best.match.score = confidence$best.match.score,
                    second.match.score = confidence$second.match.score,
                    best.match.margin = confidence$best.match.margin,
                    best.match.angle = confidence$best.match.angle,
                    best.length.relative.error = confidence$best.length.relative.error,
                    best.transport.weight = confidence$best.transport.weight,
                    reciprocal.best.match = confidence$reciprocal.best.match,
                    match.status = "soft_matched",
                    match.accepted = TRUE,
                    drop.reason = NA_character_,
                    n.exact.matches = NA_integer_,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    row.table <- if (length(row.records)) {
        do.call(rbind, row.records)
    } else {
        data.frame()
    }
    dropped.table <- if (length(dropped.records)) {
        do.call(rbind, dropped.records)
    } else {
        data.frame()
    }
    list(
        row.table = row.table,
        dropped.table = dropped.table,
        n.candidates = candidate,
        scales = scales,
        embedding.table = data.frame()
    )
}

.transported.graph.hessian.edge.angle.transport <- function(darts,
                                                            coordinates,
                                                            transport.rule,
                                                            edge.angle.args) {
    dart.vectors <- coordinates[darts$to, , drop = FALSE] -
        coordinates[darts$from, , drop = FALSE]
    dart.norms <- sqrt(rowSums(dart.vectors^2))
    if (any(!is.finite(dart.norms) | dart.norms <= 0)) {
        stop("coordinates must give nonzero direction vectors for every graph dart.",
             call. = FALSE)
    }
    scales <- .transported.graph.hessian.edge.angle.scales(
        darts = darts,
        dart.vectors = dart.vectors,
        dart.norms = dart.norms,
        edge.angle.args = edge.angle.args
    )
    hard <- identical(transport.rule, "edge.angle.hard")
    row.records <- list()
    dropped.records <- list()
    row <- 0L
    candidate <- 0L

    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        local.ids <- which(darts$from == base$from)
        target.ids <- which(darts$from == base$to)
        if (length(target.ids) < 1L) next
        for (direction.id in local.ids) {
            candidate <- candidate + 1L
            direction <- darts[direction.id, , drop = FALSE]
            scores <- .transported.graph.hessian.edge.angle.scores(
                base.idx = base.idx,
                direction.id = direction.id,
                target.ids = target.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms,
                scales = scales
            )
            if (!length(scores$total) || all(!is.finite(scores$total))) {
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = NULL,
                        drop.reason = "no_valid_edge_angle_match"
                    )
                next
            }
            if (hard) {
                weights <- rep(0, length(target.ids))
                weights[which.min(scores$total)] <- 1
            } else {
                weights <- .transported.graph.hessian.softmax(scores$total,
                                                              scales$bandwidth)
            }
            entropy <- -sum(weights * log(pmax(weights, .Machine$double.xmin)))
            effective.matches <- exp(entropy)
            best.target.id <- target.ids[which.min(scores$total)]
            reverse.scores <- .transported.graph.hessian.edge.angle.scores(
                base.idx = best.target.id,
                direction.id = best.target.id,
                target.ids = local.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms,
                scales = scales
            )
            reciprocal.best <- if (length(reverse.scores$total) &&
                                   any(is.finite(reverse.scores$total))) {
                local.ids[which.min(reverse.scores$total)] == direction.id
            } else {
                NA
            }
            confidence <- .transported.graph.hessian.match.confidence(
                scores = scores,
                weights = weights,
                direction.length = direction$length,
                matched.lengths = darts$length[target.ids],
                entropy = entropy,
                effective.matches = effective.matches,
                reciprocal.best = reciprocal.best
            )
            gate <- .transported.graph.hessian.match.gate(confidence, edge.angle.args)
            if (!gate$accepted) {
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = confidence,
                        drop.reason = gate$reason
                    )
                next
            }
            row <- row + 1L
            keep.idx <- if (hard) which(weights > 0) else seq_along(target.ids)
            for (kk in keep.idx) {
                matched <- darts[target.ids[kk], , drop = FALSE]
                row.records[[length(row.records) + 1L]] <- data.frame(
                    row = row,
                    candidate = candidate,
                    base.dart = base$dart,
                    base.from = base$from,
                    base.to = base$to,
                    base.length = base$length,
                    direction.dart = direction$dart,
                    direction.from = direction$from,
                    direction.to = direction$to,
                    direction.length = direction$length,
                    direction.label = direction$direction.label,
                    matched.dart = matched$dart,
                    matched.from = matched$from,
                    matched.to = matched$to,
                    matched.length = matched$length,
                    matched.direction.label = matched$direction.label,
                    transport.weight = weights[kk],
                    transport.entropy = entropy,
                    effective.matches = effective.matches,
                    effective.match.fraction = confidence$effective.match.fraction,
                    source.angle = scores$source.angle[kk],
                    target.angle = scores$target.angle[kk],
                    angle = scores$angle[kk],
                    length.difference = scores$length.difference[kk],
                    length.relative.error = scores$length.relative.error[kk],
                    match.score = scores$total[kk],
                    is.best.match = kk == confidence$best.index,
                    n.match.candidates = confidence$n.match.candidates,
                    best.score.local.rank = confidence$best.score.local.rank,
                    best.score.local.quantile = confidence$best.score.local.quantile,
                    best.score.robust.z = confidence$best.score.robust.z,
                    margin.local.quantile = confidence$margin.local.quantile,
                    margin.robust.z = confidence$margin.robust.z,
                    best.match.score = confidence$best.match.score,
                    second.match.score = confidence$second.match.score,
                    best.match.margin = confidence$best.match.margin,
                    best.match.angle = confidence$best.match.angle,
                    best.length.relative.error = confidence$best.length.relative.error,
                    best.transport.weight = confidence$best.transport.weight,
                    reciprocal.best.match = confidence$reciprocal.best.match,
                    match.status = if (hard) "edge_angle_hard_matched" else "edge_angle_soft_matched",
                    match.accepted = TRUE,
                    drop.reason = NA_character_,
                    n.exact.matches = NA_integer_,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    row.table <- if (length(row.records)) do.call(rbind, row.records) else data.frame()
    dropped.table <- if (length(dropped.records)) {
        do.call(rbind, dropped.records)
    } else {
        data.frame()
    }
    list(
        row.table = row.table,
        dropped.table = dropped.table,
        n.candidates = candidate,
        scales = c(scales, list(
            transport.rule = transport.rule,
            representation = "edge_angle_relative_to_base_dart"
        )),
        embedding.table = data.frame()
    )
}

.transported.graph.hessian.regression.gradient.transport <- function(darts,
                                                                     graph,
                                                                     coordinates,
                                                                     ridge,
                                                                     gradient.args) {
    if (identical(gradient.args$coordinate.method, "local.embedding")) {
        return(.transported.graph.hessian.regression.gradient.transport.local(
            darts = darts,
            graph = graph,
            coordinates = coordinates,
            ridge = ridge,
            gradient.args = gradient.args
        ))
    }

    gradients <- .transported.graph.hessian.gradient.coefficients(
        graph = graph,
        coordinates = coordinates,
        ridge = ridge
    )
    n.dim <- ncol(coordinates)
    row.records <- list()
    row <- 0L
    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        for (component in seq_len(n.dim)) {
            row <- row + 1L
            from.coef <- gradients$coefficients[[base$from]][[component]]
            to.coef <- gradients$coefficients[[base$to]][[component]]
            combined.vertices <- c(to.coef$vertex, from.coef$vertex)
            combined.values <- c(to.coef$value, -from.coef$value)
            collapsed <- rowsum(combined.values, combined.vertices, reorder = FALSE)
            vertices <- as.integer(rownames(collapsed))
            values <- as.double(collapsed[, 1L])
            keep <- is.finite(values) & abs(values) > sqrt(.Machine$double.eps)
            if (!any(keep)) {
                next
            }
            row.records[[length(row.records) + 1L]] <- data.frame(
                row = row,
                candidate = row,
                base.dart = base$dart,
                base.from = base$from,
                base.to = base$to,
                base.length = base$length,
                component = component,
                coefficient.vertex = vertices[keep],
                coefficient.value = values[keep],
                transport.weight = 1,
                transport.entropy = NA_real_,
                match.status = "regression_gradient",
                gradient.coordinate.method = "coordinates",
                gradient.embedding.method = NA_character_,
                embedding.backend = NA_character_,
                disk.size = NA_integer_,
                edge.stress = NA_real_,
                graph.distance.stress = NA_real_,
                local.design.rank.from = gradients$diagnostics$rank[base$from],
                local.design.rank.to = gradients$diagnostics$rank[base$to],
                local.condition.from = gradients$diagnostics$condition[base$from],
                local.condition.to = gradients$diagnostics$condition[base$to],
                n.local.from = gradients$diagnostics$n.local[base$from],
                n.local.to = gradients$diagnostics$n.local[base$to],
                gradient.ridge = ridge,
                stringsAsFactors = FALSE
            )
        }
    }
    row.table <- if (length(row.records)) do.call(rbind, row.records) else data.frame()
    list(
        row.table = row.table,
        dropped.table = data.frame(),
        n.candidates = nrow(darts) * n.dim,
        scales = list(
            gradient.ridge = ridge,
            gradient.coordinate.method = "coordinates",
            representation = "local_weighted_least_squares"
        ),
        embedding.table = data.frame(),
        gradient.diagnostics = gradients$diagnostics
    )
}

.transported.graph.hessian.regression.gradient.transport.local <- function(darts,
                                                                            graph,
                                                                            coordinates,
                                                                            ridge,
                                                                            gradient.args) {
    ambient.coordinates <- if (is.null(coordinates)) NULL else as.matrix(coordinates)
    n.dim <- gradient.args$embedding.dim
    row.records <- list()
    diagnostics.records <- list()
    embedding.records <- list()
    row <- 0L

    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        chart.candidates <- .transported.graph.hessian.regression.gradient.chart.candidates(
            base = base,
            darts = darts,
            graph = graph,
            ambient.coordinates = ambient.coordinates,
            ridge = ridge,
            gradient.args = gradient.args
        )
        selected.idx <- .transported.graph.hessian.select.gradient.chart(chart.candidates)
        oracle.idx <- .transported.graph.hessian.oracle.gradient.chart(chart.candidates)
        for (idx in seq_along(chart.candidates)) {
            rec <- chart.candidates[[idx]]$embedding.record
            rec$chart.selected <- idx == selected.idx
            rec$chart.oracle <- idx == oracle.idx
            rec$gradient.chart.selection <- gradient.args$chart.selection
            embedding.records[[length(embedding.records) + 1L]] <- rec
        }

        selected <- chart.candidates[[selected.idx]]
        from.diag <- selected$from.gradient$diagnostics
        to.diag <- selected$to.gradient$diagnostics
        from.diag$base.dart <- base$dart
        from.diag$base.from <- base$from
        from.diag$base.to <- base$to
        from.diag$endpoint.role <- "from"
        from.diag$gradient.coordinate.method <- "local.embedding"
        from.diag$gradient.embedding.method <- selected$embedding.method
        from.diag$gradient.disk.rule <- selected$disk.rule
        from.diag$gradient.disk.hops <- selected$disk.hops
        from.diag$gradient.disk.radius.fraction <- selected$disk.radius.fraction
        from.diag$gradient.disk.local.scale.multiplier <-
            selected$disk.local.scale.multiplier
        from.diag$gradient.disk.local.scale <- selected$disk$local.scale %||% NA_real_
        from.diag$embedding.backend <- selected$embed$backend.used
        from.diag$disk.size <- length(selected$disk$vertices)
        from.diag$chart.selection.score <- selected$selection.score
        to.diag$base.dart <- base$dart
        to.diag$base.from <- base$from
        to.diag$base.to <- base$to
        to.diag$endpoint.role <- "to"
        to.diag$gradient.coordinate.method <- "local.embedding"
        to.diag$gradient.embedding.method <- selected$embedding.method
        to.diag$gradient.disk.rule <- selected$disk.rule
        to.diag$gradient.disk.hops <- selected$disk.hops
        to.diag$gradient.disk.radius.fraction <- selected$disk.radius.fraction
        to.diag$gradient.disk.local.scale.multiplier <-
            selected$disk.local.scale.multiplier
        to.diag$gradient.disk.local.scale <- selected$disk$local.scale %||% NA_real_
        to.diag$embedding.backend <- selected$embed$backend.used
        to.diag$disk.size <- length(selected$disk$vertices)
        to.diag$chart.selection.score <- selected$selection.score
        diagnostics.records[[length(diagnostics.records) + 1L]] <- from.diag
        diagnostics.records[[length(diagnostics.records) + 1L]] <- to.diag

        for (component in seq_len(ncol(selected$embed$coordinates))) {
            row <- row + 1L
            from.coef <- selected$from.gradient$coefficients[[component]]
            to.coef <- selected$to.gradient$coefficients[[component]]
            combined.vertices <- c(to.coef$vertex, from.coef$vertex)
            combined.values <- c(to.coef$value, -from.coef$value)
            collapsed <- rowsum(combined.values, combined.vertices, reorder = FALSE)
            vertices <- as.integer(rownames(collapsed))
            values <- as.double(collapsed[, 1L])
            keep <- is.finite(values) & abs(values) > sqrt(.Machine$double.eps)
            if (!any(keep)) {
                next
            }
            row.records[[length(row.records) + 1L]] <- data.frame(
                row = row,
                candidate = row,
                base.dart = base$dart,
                base.from = base$from,
                base.to = base$to,
                base.length = base$length,
                component = component,
                coefficient.vertex = vertices[keep],
                coefficient.value = values[keep],
                transport.weight = 1,
                transport.entropy = NA_real_,
                match.status = "regression_gradient",
                gradient.coordinate.method = "local.embedding",
                gradient.embedding.method = selected$embedding.method,
                gradient.disk.rule = selected$disk.rule,
                gradient.disk.parameter = selected$disk.parameter,
                gradient.disk.hops = selected$disk.hops,
                gradient.disk.radius.fraction = selected$disk.radius.fraction,
                gradient.disk.local.scale.multiplier =
                    selected$disk.local.scale.multiplier,
                gradient.disk.local.scale = selected$disk$local.scale %||% NA_real_,
                gradient.disk.metric.radius = selected$disk$metric.radius %||% NA_real_,
                gradient.chart.selection = gradient.args$chart.selection,
                embedding.backend = selected$embed$backend.used,
                disk.size = length(selected$disk$vertices),
                disk.fraction = selected$locality$disk.fraction,
                disk.max.hop = selected$locality$disk.max.hop,
                disk.mean.hop = selected$locality$disk.mean.hop,
                disk.weighted.radius = selected$locality$disk.weighted.radius,
                disk.weighted.mean.radius = selected$locality$disk.weighted.mean.radius,
                edge.stress = selected$embed$edge.stress,
                graph.distance.stress = selected$embed$graph.distance.stress,
                chart.selection.score = selected$selection.score,
                chart.oracle = selected.idx == oracle.idx,
                ambient.affine.residual = selected$chart.diagnostics$ambient.affine.residual,
                ambient.linear.gradient.residual.mean =
                    selected$chart.diagnostics$ambient.linear.gradient.residual.mean,
                ambient.linear.gradient.residual.max =
                    selected$chart.diagnostics$ambient.linear.gradient.residual.max,
                local.design.rank.from = selected$from.gradient$diagnostics$rank,
                local.design.rank.to = selected$to.gradient$diagnostics$rank,
                local.condition.from = selected$from.gradient$diagnostics$condition,
                local.condition.to = selected$to.gradient$diagnostics$condition,
                n.local.from = selected$from.gradient$diagnostics$n.local,
                n.local.to = selected$to.gradient$diagnostics$n.local,
                gradient.ridge = ridge,
                stringsAsFactors = FALSE
            )
        }
    }

    row.table <- if (length(row.records)) do.call(rbind, row.records) else data.frame()
    gradient.diagnostics <- if (length(diagnostics.records)) {
        do.call(rbind, diagnostics.records)
    } else {
        data.frame()
    }
    embedding.table <- if (length(embedding.records)) {
        do.call(rbind, embedding.records)
    } else {
        data.frame()
    }
    list(
        row.table = row.table,
        dropped.table = data.frame(),
        n.candidates = nrow(darts) * n.dim,
        scales = list(
            gradient.ridge = ridge,
            gradient.coordinate.method = "local.embedding",
            gradient.embedding.method = gradient.args$embedding.method,
            gradient.chart.selection = gradient.args$chart.selection,
            gradient.embedding.candidates = paste(gradient.args$embedding.candidates,
                                                  collapse = ","),
            gradient.disk.rule = gradient.args$disk.rule,
            gradient.disk.radius.fraction = gradient.args$disk.radius.fraction,
            gradient.disk.local.scale.method = gradient.args$disk.local.scale.method,
            gradient.disk.local.scale.k = gradient.args$disk.local.scale.k,
            gradient.disk.local.scale.quantile = gradient.args$disk.local.scale.quantile,
            gradient.disk.local.scale.multiplier =
                gradient.args$disk.local.scale.multiplier,
            gradient.disk.min.vertices = gradient.args$disk.min.vertices,
            gradient.disk.rule.candidates = paste(gradient.args$disk.rule.candidates,
                                                  collapse = ","),
            gradient.disk.hops.candidates = paste(gradient.args$disk.hops.candidates,
                                                  collapse = ","),
            gradient.disk.radius.fraction.candidates =
                paste(gradient.args$disk.radius.fraction.candidates, collapse = ","),
            gradient.disk.local.scale.multiplier.candidates =
                paste(gradient.args$disk.local.scale.multiplier.candidates,
                      collapse = ","),
            representation = "shared_local_embedding_weighted_least_squares"
        ),
        embedding.table = embedding.table,
        gradient.diagnostics = gradient.diagnostics
    )
}

.transported.graph.hessian.regression.gradient.chart.candidates <- function(base,
                                                                            darts,
                                                                            graph,
                                                                            ambient.coordinates,
                                                                            ridge,
                                                                            gradient.args) {
    local.ids <- which(darts$from == base$from)
    target.ids <- which(darts$from == base$to)
    required <- unique(c(base$from, base$to, darts$to[local.ids], darts$to[target.ids]))
    if (identical(gradient.args$chart.selection, "adaptive")) {
        methods <- gradient.args$embedding.candidates
        disk.rules <- gradient.args$disk.rule.candidates
    } else {
        methods <- gradient.args$embedding.method
        disk.rules <- gradient.args$disk.rule
    }
    candidates <- list()
    for (method in methods) {
        for (disk.rule in disk.rules) {
            disk.parameters <- if (identical(disk.rule, "hops")) {
                if (identical(gradient.args$chart.selection, "adaptive")) {
                    gradient.args$disk.hops.candidates
                } else {
                    gradient.args$disk.hops
                }
            } else if (identical(disk.rule, "metric.diameter.fraction")) {
                if (identical(gradient.args$chart.selection, "adaptive")) {
                    gradient.args$disk.radius.fraction.candidates
                } else {
                    gradient.args$disk.radius.fraction
                }
            } else {
                if (identical(gradient.args$chart.selection, "adaptive")) {
                    gradient.args$disk.local.scale.multiplier.candidates
                } else {
                    gradient.args$disk.local.scale.multiplier
                }
            }
            for (disk.parameter in disk.parameters) {
                candidates[[length(candidates) + 1L]] <-
                    .transported.graph.hessian.regression.gradient.chart.candidate(
                        base = base,
                        graph = graph,
                        required = required,
                        ambient.coordinates = ambient.coordinates,
                        ridge = ridge,
                        embedding.method = method,
                        disk.rule = disk.rule,
                        disk.parameter = disk.parameter,
                        embedding.dim = gradient.args$embedding.dim,
                        max.vertices = gradient.args$max.vertices,
                        min.vertices = gradient.args$disk.min.vertices,
                        local.scales = gradient.args$local.scales,
                        graph.diameter = gradient.args$graph.diameter
                    )
            }
        }
    }
    candidates
}

.transported.graph.hessian.regression.gradient.chart.candidate <- function(base,
                                                                           graph,
                                                                           required,
                                                                           ambient.coordinates,
                                                                           ridge,
                                                                           embedding.method,
                                                                           disk.rule,
                                                                           disk.parameter,
                                                                           embedding.dim,
                                                                           max.vertices,
                                                                           min.vertices,
                                                                           local.scales,
                                                                           graph.diameter) {
    if (identical(disk.rule, "hops")) {
        disk <- .transported.graph.hessian.local.disk(
            graph = graph,
            seeds = c(base$from, base$to),
            required = required,
            hops = disk.parameter,
            max.vertices = max.vertices
        )
    } else if (identical(disk.rule, "metric.diameter.fraction")) {
        disk <- .transported.graph.hessian.metric.local.disk(
            graph = graph,
            seeds = c(base$from, base$to),
            required = required,
            radius.fraction = disk.parameter,
            graph.diameter = graph.diameter,
            min.vertices = min.vertices,
            max.vertices = max.vertices
        )
    } else {
        disk <- .transported.graph.hessian.metric.local.scale.disk(
            graph = graph,
            seeds = c(base$from, base$to),
            required = required,
            local.scales = local.scales,
            multiplier = disk.parameter,
            min.vertices = min.vertices,
            max.vertices = max.vertices
        )
    }
    embed <- .transported.graph.hessian.local.embedding(
        graph = graph,
        vertices = disk$vertices,
        method = embedding.method,
        dim = embedding.dim
    )
    local.coordinates <- matrix(NA_real_,
                                nrow = graph$n.vertices,
                                ncol = ncol(embed$coordinates))
    local.coordinates[disk$vertices, ] <- embed$coordinates
    from.gradient <- .transported.graph.hessian.gradient.coefficients.one(
        graph = graph,
        coordinates = local.coordinates,
        ridge = ridge,
        vertex = base$from,
        available.vertices = disk$vertices
    )
    to.gradient <- .transported.graph.hessian.gradient.coefficients.one(
        graph = graph,
        coordinates = local.coordinates,
        ridge = ridge,
        vertex = base$to,
        available.vertices = disk$vertices
    )
    chart.diagnostics <- .transported.graph.hessian.local.chart.diagnostics(
        embedded.coordinates = embed$coordinates,
        ambient.coordinates = ambient.coordinates,
        disk.vertices = disk$vertices,
        from.gradient = from.gradient,
        to.gradient = to.gradient
    )
    locality <- .transported.graph.hessian.chart.locality(
        graph = graph,
        disk = disk,
        embed = embed,
        seeds = c(base$from, base$to)
    )
    selection.score <- .transported.graph.hessian.gradient.chart.score(
        embed = embed,
        from.gradient = from.gradient,
        to.gradient = to.gradient,
        embedding.dim = embedding.dim,
        locality = locality
    )
    embedding.record <- data.frame(
        base.dart = base$dart,
        base.from = base$from,
        base.to = base$to,
        requested.method = embedding.method,
        gradient.disk.rule = disk.rule,
        gradient.disk.parameter = disk.parameter,
        gradient.disk.hops = if (identical(disk.rule, "hops")) as.integer(disk.parameter) else NA_integer_,
        gradient.disk.radius.fraction = if (identical(disk.rule, "metric.diameter.fraction")) as.double(disk.parameter) else NA_real_,
        gradient.disk.local.scale.multiplier = if (identical(disk.rule, "metric.local.scale")) as.double(disk.parameter) else NA_real_,
        gradient.disk.local.scale = disk$local.scale %||% NA_real_,
        gradient.disk.metric.radius = disk$metric.radius %||% NA_real_,
        backend.used = embed$backend.used,
        fallback.reason = embed$fallback.reason,
        disk.size = length(disk$vertices),
        disk.fraction = locality$disk.fraction,
        disk.max.hop = locality$disk.max.hop,
        disk.mean.hop = locality$disk.mean.hop,
        disk.weighted.radius = locality$disk.weighted.radius,
        disk.weighted.mean.radius = locality$disk.weighted.mean.radius,
        embedding.dim = ncol(embed$coordinates),
        edge.stress = embed$edge.stress,
        graph.distance.stress = embed$graph.distance.stress,
        local.design.rank.from = from.gradient$diagnostics$rank,
        local.design.rank.to = to.gradient$diagnostics$rank,
        local.condition.from = from.gradient$diagnostics$condition,
        local.condition.to = to.gradient$diagnostics$condition,
        n.local.from = from.gradient$diagnostics$n.local,
        n.local.to = to.gradient$diagnostics$n.local,
        chart.selection.score = selection.score,
        ambient.affine.residual = chart.diagnostics$ambient.affine.residual,
        ambient.linear.gradient.residual.mean =
            chart.diagnostics$ambient.linear.gradient.residual.mean,
        ambient.linear.gradient.residual.max =
            chart.diagnostics$ambient.linear.gradient.residual.max,
        stringsAsFactors = FALSE
    )
    list(
        embedding.method = embedding.method,
        disk.rule = disk.rule,
        disk.parameter = disk.parameter,
        disk.hops = if (identical(disk.rule, "hops")) as.integer(disk.parameter) else NA_integer_,
        disk.radius.fraction = if (identical(disk.rule, "metric.diameter.fraction")) as.double(disk.parameter) else NA_real_,
        disk.local.scale.multiplier = if (identical(disk.rule, "metric.local.scale")) as.double(disk.parameter) else NA_real_,
        disk = disk,
        embed = embed,
        from.gradient = from.gradient,
        to.gradient = to.gradient,
        chart.diagnostics = chart.diagnostics,
        locality = locality,
        selection.score = selection.score,
        embedding.record = embedding.record
    )
}

.transported.graph.hessian.gradient.chart.score <- function(embed,
                                                            from.gradient,
                                                            to.gradient,
                                                            embedding.dim,
                                                            locality) {
    edge <- embed$edge.stress
    dist <- embed$graph.distance.stress
    if (!is.finite(edge)) edge <- 1e6
    if (!is.finite(dist)) dist <- 1e6
    ranks <- c(from.gradient$diagnostics$rank, to.gradient$diagnostics$rank)
    rank.penalty <- mean(pmax(0, embedding.dim - ranks))
    conditions <- c(from.gradient$diagnostics$condition,
                    to.gradient$diagnostics$condition)
    finite.conditions <- conditions[is.finite(conditions)]
    condition.term <- if (length(finite.conditions)) {
        log1p(max(finite.conditions))
    } else {
        50
    }
    radius <- locality$disk.weighted.radius
    if (!is.finite(radius)) radius <- 0
    edge + dist + 10 * rank.penalty + 0.05 * condition.term +
        0.05 * locality$disk.fraction + 0.01 * radius
}

.transported.graph.hessian.select.gradient.chart <- function(candidates) {
    scores <- vapply(candidates, `[[`, numeric(1), "selection.score")
    scores[!is.finite(scores)] <- Inf
    which.min(scores)
}

.transported.graph.hessian.oracle.gradient.chart <- function(candidates) {
    oracle <- vapply(candidates, function(candidate) {
        candidate$chart.diagnostics$ambient.linear.gradient.residual.mean
    }, numeric(1))
    if (!any(is.finite(oracle))) {
        return(.transported.graph.hessian.select.gradient.chart(candidates))
    }
    oracle[!is.finite(oracle)] <- Inf
    which.min(oracle)
}

.transported.graph.hessian.chart.locality <- function(graph, disk, embed, seeds) {
    finite.hops <- disk$hop.distance[is.finite(disk$hop.distance)]
    seed.idx <- match(seeds, disk$vertices)
    weighted.to.seed <- NA_real_
    weighted.radius <- NA_real_
    if (!is.null(embed$distances) && all(!is.na(seed.idx))) {
        d <- embed$distances[, seed.idx, drop = FALSE]
        weighted.to.seed <- apply(d, 1L, min)
        finite.weighted <- weighted.to.seed[is.finite(weighted.to.seed)]
        if (length(finite.weighted)) {
            weighted.radius <- max(finite.weighted)
            weighted.to.seed <- mean(finite.weighted)
        } else {
            weighted.to.seed <- NA_real_
        }
    }
    list(
        disk.fraction = length(disk$vertices) / graph$n.vertices,
        disk.max.hop = if (length(finite.hops)) max(finite.hops) else NA_real_,
        disk.mean.hop = if (length(finite.hops)) mean(finite.hops) else NA_real_,
        disk.weighted.radius = weighted.radius,
        disk.weighted.mean.radius = weighted.to.seed
    )
}

.transported.graph.hessian.gradient.coefficients <- function(graph,
                                                             coordinates,
                                                             ridge) {
    coordinates <- as.matrix(coordinates)
    n.vertices <- graph$n.vertices
    n.dim <- ncol(coordinates)
    coefficients <- vector("list", n.vertices)
    diagnostics <- vector("list", n.vertices)
    for (vertex in seq_len(n.vertices)) {
        one <- .transported.graph.hessian.gradient.coefficients.one(
            graph = graph,
            coordinates = coordinates,
            ridge = ridge,
            vertex = vertex,
            available.vertices = seq_len(n.vertices)
        )
        coefficients[[vertex]] <- one$coefficients
        diagnostics[[vertex]] <- one$diagnostics
    }
    list(coefficients = coefficients, diagnostics = do.call(rbind, diagnostics))
}

.transported.graph.hessian.gradient.coefficients.one <- function(graph,
                                                                 coordinates,
                                                                 ridge,
                                                                 vertex,
                                                                 available.vertices) {
    coordinates <- as.matrix(coordinates)
    n.dim <- ncol(coordinates)
    available <- rep(FALSE, graph$n.vertices)
    available[available.vertices] <- TRUE
    vertex.coord <- coordinates[vertex, , drop = TRUE]
    nbrs <- graph$adj.list[[vertex]]
    keep <- nbrs >= 1L & nbrs <= graph$n.vertices & available[nbrs]
    if (length(nbrs)) {
        nbr.coords <- coordinates[nbrs, , drop = FALSE]
        keep <- keep & apply(nbr.coords, 1L, function(z) all(is.finite(z)))
    }
    if (!all(is.finite(vertex.coord))) {
        keep[] <- FALSE
    }
    nbrs <- nbrs[keep]
    edge.length <- graph$weight.list[[vertex]][keep]
    if (!length(nbrs)) {
        coefficients <- lapply(seq_len(n.dim), function(component) {
            data.frame(vertex = vertex, value = 0, stringsAsFactors = FALSE)
        })
        diagnostics <- data.frame(
            vertex = vertex,
            n.local = 0L,
            rank = 0L,
            condition = Inf,
            stringsAsFactors = FALSE
        )
        return(list(coefficients = coefficients, diagnostics = diagnostics))
    }
    design <- coordinates[nbrs, , drop = FALSE] -
        matrix(vertex.coord, nrow = length(nbrs),
               ncol = n.dim, byrow = TRUE)
    weights <- 1 / pmax(edge.length, sqrt(.Machine$double.eps))^2
    weighted.design <- design * sqrt(weights)
    normal <- crossprod(weighted.design) + diag(ridge, n.dim)
    rhs <- t(design) %*% diag(weights, nrow = length(weights))
    inv.normal <- try(solve(normal), silent = TRUE)
    if (inherits(inv.normal, "try-error")) {
        inv.normal <- MASS::ginv(normal)
    }
    mapping <- inv.normal %*% rhs
    sv <- svd(normal, nu = 0, nv = 0)$d
    coefficients <- lapply(seq_len(n.dim), function(component) {
        vals <- as.double(mapping[component, ])
        data.frame(
            vertex = c(vertex, nbrs),
            value = c(-sum(vals), vals),
            stringsAsFactors = FALSE
        )
    })
    diagnostics <- data.frame(
        vertex = vertex,
        n.local = length(nbrs),
        rank = qr(weighted.design)$rank,
        condition = if (length(sv) && min(sv) > 0) max(sv) / min(sv) else Inf,
        stringsAsFactors = FALSE
    )
    list(coefficients = coefficients, diagnostics = diagnostics)
}

.transported.graph.hessian.local.chart.diagnostics <- function(embedded.coordinates,
                                                               ambient.coordinates,
                                                               disk.vertices,
                                                               from.gradient,
                                                               to.gradient) {
    empty <- list(
        ambient.affine.residual = NA_real_,
        ambient.linear.gradient.residual.mean = NA_real_,
        ambient.linear.gradient.residual.max = NA_real_
    )
    if (is.null(ambient.coordinates)) return(empty)
    ambient.coordinates <- as.matrix(ambient.coordinates)
    if (!nrow(ambient.coordinates) || any(disk.vertices < 1L) ||
        any(disk.vertices > nrow(ambient.coordinates))) {
        return(empty)
    }
    ambient <- ambient.coordinates[disk.vertices, , drop = FALSE]
    embedded <- as.matrix(embedded.coordinates)
    if (nrow(ambient) != nrow(embedded) ||
        any(!is.finite(ambient)) || any(!is.finite(embedded))) {
        return(empty)
    }

    design <- cbind(intercept = 1, embedded)
    fit <- try(qr.coef(qr(design), ambient), silent = TRUE)
    ambient.affine.residual <- NA_real_
    if (!inherits(fit, "try-error") && all(is.finite(fit))) {
        fitted <- design %*% fit
        denom <- sqrt(sum((ambient - matrix(colMeans(ambient),
                                            nrow = nrow(ambient),
                                            ncol = ncol(ambient),
                                            byrow = TRUE))^2))
        if (!is.finite(denom) || denom <= 0) denom <- sqrt(sum(ambient^2))
        if (!is.finite(denom) || denom <= 0) denom <- 1
        ambient.affine.residual <- sqrt(sum((ambient - fitted)^2)) / denom
    }

    gradient.values <- numeric()
    for (component in seq_along(from.gradient$coefficients)) {
        from.coef <- from.gradient$coefficients[[component]]
        to.coef <- to.gradient$coefficients[[component]]
        for (ambient.dim in seq_len(ncol(ambient.coordinates))) {
            from.value <- sum(from.coef$value *
                                  ambient.coordinates[from.coef$vertex, ambient.dim])
            to.value <- sum(to.coef$value *
                                ambient.coordinates[to.coef$vertex, ambient.dim])
            gradient.values <- c(gradient.values, to.value - from.value)
        }
    }
    gradient.values <- gradient.values[is.finite(gradient.values)]
    list(
        ambient.affine.residual = ambient.affine.residual,
        ambient.linear.gradient.residual.mean = if (length(gradient.values)) {
            mean(abs(gradient.values))
        } else {
            NA_real_
        },
        ambient.linear.gradient.residual.max = if (length(gradient.values)) {
            max(abs(gradient.values))
        } else {
            NA_real_
        }
    )
}

.transported.graph.hessian.soft.transport.local <- function(darts, graph,
                                                            local.embedding.method,
                                                            soft.args) {
    row.records <- list()
    dropped.records <- list()
    embedding.records <- list()
    row <- 0L
    candidate <- 0L

    for (base.idx in seq_len(nrow(darts))) {
        base <- darts[base.idx, , drop = FALSE]
        local.ids <- which(darts$from == base$from)
        target.ids <- which(darts$from == base$to)
        if (length(target.ids) < 1L) next
        disk <- .transported.graph.hessian.local.disk(
            graph = graph,
            seeds = c(base$from, base$to),
            required = unique(c(base$from, base$to,
                                darts$to[local.ids], darts$to[target.ids])),
            hops = soft.args$local.disk.hops,
            max.vertices = soft.args$local.max.vertices
        )
        embed <- .transported.graph.hessian.local.embedding(
            graph = graph,
            vertices = disk$vertices,
            method = local.embedding.method,
            dim = soft.args$local.embedding.dim
        )
        embedding.records[[length(embedding.records) + 1L]] <- data.frame(
            base.dart = base$dart,
            base.from = base$from,
            base.to = base$to,
            requested.method = local.embedding.method,
            backend.used = embed$backend.used,
            fallback.reason = embed$fallback.reason,
            disk.size = length(disk$vertices),
            embedding.dim = ncol(embed$coordinates),
            edge.stress = embed$edge.stress,
            graph.distance.stress = embed$graph.distance.stress,
            stringsAsFactors = FALSE
        )
        local.map <- match(seq_len(graph$n.vertices), disk$vertices)
        local.coords <- embed$coordinates
        local.vec <- function(dart.id) {
            from <- local.map[darts$from[dart.id]]
            to <- local.map[darts$to[dart.id]]
            local.coords[to, , drop = TRUE] - local.coords[from, , drop = TRUE]
        }

        for (direction.id in local.ids) {
            candidate <- candidate + 1L
            direction <- darts[direction.id, , drop = FALSE]
            direction.vec <- local.vec(direction.id)
            direction.norm <- sqrt(sum(direction.vec^2))
            if (!is.finite(direction.norm) || direction.norm <= 0) {
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = NULL,
                        drop.reason = "invalid_direction_vector",
                        embedding.backend = embed$backend.used,
                        disk.size = length(disk$vertices)
                    )
                next
            }
            target.vecs <- t(vapply(target.ids, local.vec, numeric(ncol(local.coords))))
            target.norms <- sqrt(rowSums(target.vecs^2))
            valid <- is.finite(target.norms) & target.norms > 0
            if (!any(valid)) {
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = NULL,
                        drop.reason = "no_valid_target_vectors",
                        embedding.backend = embed$backend.used,
                        disk.size = length(disk$vertices)
                    )
                next
            }
            use.target.ids <- target.ids[valid]
            target.vecs <- target.vecs[valid, , drop = FALSE]
            target.norms <- target.norms[valid]
            scores <- .transported.graph.hessian.soft.scores.from.vectors(
                direction.vec = direction.vec,
                direction.norm = direction.norm,
                direction.length = direction$length,
                target.ids = use.target.ids,
                target.vecs = target.vecs,
                target.norms = target.norms,
                darts = darts,
                soft.args = soft.args
            )
            weights <- .transported.graph.hessian.softmax(scores$total, soft.args$bandwidth)
            entropy <- -sum(weights * log(pmax(weights, .Machine$double.xmin)))
            effective.matches <- exp(entropy)
            best.target.id <- use.target.ids[which.min(scores$total)]
            reverse.direction.vec <- local.vec(best.target.id)
            reverse.direction.norm <- sqrt(sum(reverse.direction.vec^2))
            reverse.vecs <- t(vapply(local.ids, local.vec, numeric(ncol(local.coords))))
            reverse.norms <- sqrt(rowSums(reverse.vecs^2))
            reverse.valid <- is.finite(reverse.norms) & reverse.norms > 0
            reciprocal.best <- NA
            if (is.finite(reverse.direction.norm) && reverse.direction.norm > 0 &&
                any(reverse.valid)) {
                reverse.scores <- .transported.graph.hessian.soft.scores.from.vectors(
                    direction.vec = reverse.direction.vec,
                    direction.norm = reverse.direction.norm,
                    direction.length = darts$length[best.target.id],
                    target.ids = local.ids[reverse.valid],
                    target.vecs = reverse.vecs[reverse.valid, , drop = FALSE],
                    target.norms = reverse.norms[reverse.valid],
                    darts = darts,
                    soft.args = soft.args
                )
                reciprocal.best <- local.ids[reverse.valid][which.min(reverse.scores$total)] ==
                    direction.id
            }
            confidence <- .transported.graph.hessian.match.confidence(
                scores = scores,
                weights = weights,
                direction.length = direction$length,
                matched.lengths = darts$length[use.target.ids],
                entropy = entropy,
                effective.matches = effective.matches,
                reciprocal.best = reciprocal.best
            )
            gate <- .transported.graph.hessian.match.gate(confidence, soft.args)
            if (!gate$accepted) {
                dropped.records[[length(dropped.records) + 1L]] <-
                    .transported.graph.hessian.dropped.soft.record(
                        candidate = candidate,
                        base = base,
                        direction = direction,
                        confidence = confidence,
                        drop.reason = gate$reason,
                        embedding.backend = embed$backend.used,
                        disk.size = length(disk$vertices)
                    )
                next
            }
            row <- row + 1L
            for (kk in seq_along(use.target.ids)) {
                matched <- darts[use.target.ids[kk], , drop = FALSE]
                row.records[[length(row.records) + 1L]] <- data.frame(
                    row = row,
                    candidate = candidate,
                    base.dart = base$dart,
                    base.from = base$from,
                    base.to = base$to,
                    base.length = base$length,
                    direction.dart = direction$dart,
                    direction.from = direction$from,
                    direction.to = direction$to,
                    direction.length = direction$length,
                    direction.label = direction$direction.label,
                    matched.dart = matched$dart,
                    matched.from = matched$from,
                    matched.to = matched$to,
                    matched.length = matched$length,
                    matched.direction.label = matched$direction.label,
                    transport.weight = weights[kk],
                    transport.entropy = entropy,
                    effective.matches = effective.matches,
                    effective.match.fraction = confidence$effective.match.fraction,
                    angle = scores$angle[kk],
                    length.difference = scores$length.difference[kk],
                    match.score = scores$total[kk],
                    is.best.match = kk == confidence$best.index,
                    n.match.candidates = confidence$n.match.candidates,
                    best.score.local.rank = confidence$best.score.local.rank,
                    best.score.local.quantile = confidence$best.score.local.quantile,
                    best.score.robust.z = confidence$best.score.robust.z,
                    margin.local.quantile = confidence$margin.local.quantile,
                    margin.robust.z = confidence$margin.robust.z,
                    best.match.score = confidence$best.match.score,
                    second.match.score = confidence$second.match.score,
                    best.match.margin = confidence$best.match.margin,
                    best.match.angle = confidence$best.match.angle,
                    best.length.relative.error = confidence$best.length.relative.error,
                    best.transport.weight = confidence$best.transport.weight,
                    reciprocal.best.match = confidence$reciprocal.best.match,
                    match.status = "soft_matched",
                    match.accepted = TRUE,
                    drop.reason = NA_character_,
                    n.exact.matches = NA_integer_,
                    embedding.backend = embed$backend.used,
                    disk.size = length(disk$vertices),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    row.table <- if (length(row.records)) do.call(rbind, row.records) else data.frame()
    embedding.table <- if (length(embedding.records)) {
        do.call(rbind, embedding.records)
    } else {
        data.frame()
    }
    dropped.table <- if (length(dropped.records)) {
        do.call(rbind, dropped.records)
    } else {
        data.frame()
    }
    list(
        row.table = row.table,
        dropped.table = dropped.table,
        n.candidates = candidate,
        scales = list(
            angle.scale = soft.args$angle.scale %||% NA_real_,
            length.scale = soft.args$length.scale %||% NA_real_,
            bandwidth = soft.args$bandwidth,
            angle.scale.source = if (is.null(soft.args$angle.scale)) "per.row.local" else "user",
            length.scale.source = if (is.null(soft.args$length.scale)) "per.row.local" else "user",
            local.embedding.method = local.embedding.method
        ),
        embedding.table = embedding.table
    )
}

.transported.graph.hessian.match.confidence <- function(scores, weights,
                                                        direction.length,
                                                        matched.lengths,
                                                        entropy,
                                                        effective.matches,
                                                        reciprocal.best) {
    best.index <- which.min(scores$total)
    ordered <- sort(scores$total, partial = seq_len(min(2L, length(scores$total))))
    second <- if (length(ordered) >= 2L) ordered[[2L]] else Inf
    best.length.difference <- abs(matched.lengths[best.index] - direction.length)
    best.length.relative.error <- best.length.difference /
        max(abs(direction.length), sqrt(.Machine$double.eps))
    finite.scores <- scores$total[is.finite(scores$total)]
    n.candidates <- length(scores$total)
    best.rank <- rank(scores$total, ties.method = "min")[best.index]
    best.score.local.quantile <- if (n.candidates > 0L) {
        mean(scores$total <= scores$total[best.index], na.rm = TRUE)
    } else {
        NA_real_
    }
    score.center <- if (length(finite.scores)) stats::median(finite.scores) else NA_real_
    score.scale <- .transported.graph.hessian.robust.scale(finite.scores)
    best.score.robust.z <- if (is.finite(score.center) && is.finite(score.scale) &&
                               score.scale > 0) {
        (score.center - scores$total[best.index]) / score.scale
    } else {
        NA_real_
    }
    gaps <- diff(sort(finite.scores))
    gaps <- gaps[is.finite(gaps)]
    margin <- second - scores$total[best.index]
    margin.local.quantile <- if (is.finite(margin) && length(gaps)) {
        mean(gaps <= margin)
    } else if (is.infinite(margin) && margin > 0) {
        1
    } else {
        NA_real_
    }
    gap.center <- if (length(gaps)) stats::median(gaps) else 0
    gap.scale <- .transported.graph.hessian.robust.scale(gaps)
    margin.robust.z <- if (is.infinite(margin) && margin > 0) {
        Inf
    } else if (is.finite(margin) && is.finite(gap.center) &&
               is.finite(gap.scale) && gap.scale > 0) {
        (margin - gap.center) / gap.scale
    } else {
        NA_real_
    }
    effective.match.fraction <- if (n.candidates > 0L) {
        effective.matches / n.candidates
    } else {
        NA_real_
    }
    list(
        best.index = best.index,
        n.match.candidates = n.candidates,
        best.score.local.rank = best.rank,
        best.score.local.quantile = best.score.local.quantile,
        best.score.robust.z = best.score.robust.z,
        margin.local.quantile = margin.local.quantile,
        margin.robust.z = margin.robust.z,
        best.match.score = scores$total[best.index],
        second.match.score = second,
        best.match.margin = margin,
        best.match.angle = scores$angle[best.index],
        best.length.relative.error = best.length.relative.error,
        best.transport.weight = weights[best.index],
        transport.entropy = entropy,
        effective.matches = effective.matches,
        effective.match.fraction = effective.match.fraction,
        reciprocal.best.match = if (length(reciprocal.best) == 1L) {
            as.logical(reciprocal.best)
        } else {
            NA
        }
    )
}

.transported.graph.hessian.robust.scale <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) {
        return(NA_real_)
    }
    scale <- stats::mad(x, constant = 1.4826)
    if (!is.finite(scale) || scale <= 0) {
        scale <- stats::IQR(x) / 1.349
    }
    if (!is.finite(scale) || scale <= 0) {
        scale <- stats::sd(x)
    }
    if (!is.finite(scale) || scale <= 0) {
        scale <- 1
    }
    scale
}

.transported.graph.hessian.match.gate <- function(confidence, soft.args) {
    if (!is.null(soft.args$max.match.angle) &&
        is.finite(confidence$best.match.angle) &&
        confidence$best.match.angle > soft.args$max.match.angle) {
        return(list(accepted = FALSE, reason = "best_angle_above_threshold"))
    }
    if (!is.null(soft.args$max.length.relative.error) &&
        is.finite(confidence$best.length.relative.error) &&
        confidence$best.length.relative.error > soft.args$max.length.relative.error) {
        return(list(accepted = FALSE, reason = "length_relative_error_above_threshold"))
    }
    if (!is.null(soft.args$min.match.margin) &&
        is.finite(confidence$best.match.margin) &&
        confidence$best.match.margin < soft.args$min.match.margin) {
        return(list(accepted = FALSE, reason = "match_margin_below_threshold"))
    }
    if (!is.null(soft.args$max.effective.matches) &&
        is.finite(confidence$effective.matches) &&
        confidence$effective.matches > soft.args$max.effective.matches) {
        return(list(accepted = FALSE, reason = "effective_matches_above_threshold"))
    }
    if (!is.null(soft.args$max.effective.match.fraction) &&
        is.finite(confidence$effective.match.fraction) &&
        confidence$effective.match.fraction > soft.args$max.effective.match.fraction) {
        return(list(accepted = FALSE,
                    reason = "effective_match_fraction_above_threshold"))
    }
    if (identical(soft.args$match.threshold.rule, "local.quantile")) {
        if (is.finite(confidence$best.score.local.quantile) &&
            confidence$best.score.local.quantile > soft.args$match.score.quantile) {
            return(list(accepted = FALSE,
                        reason = "best_score_quantile_above_threshold"))
        }
        if (is.finite(confidence$margin.local.quantile) &&
            confidence$margin.local.quantile < soft.args$match.margin.quantile) {
            return(list(accepted = FALSE,
                        reason = "margin_quantile_below_threshold"))
        }
    }
    if (identical(soft.args$match.threshold.rule, "local.robust.z")) {
        if (is.finite(confidence$best.score.robust.z) &&
            confidence$best.score.robust.z < soft.args$min.best.score.z) {
            return(list(accepted = FALSE,
                        reason = "best_score_robust_z_below_threshold"))
        }
        if (is.finite(confidence$margin.robust.z) &&
            confidence$margin.robust.z < soft.args$min.margin.z) {
            return(list(accepted = FALSE,
                        reason = "margin_robust_z_below_threshold"))
        }
    }
    list(accepted = TRUE, reason = NA_character_)
}

.transported.graph.hessian.dropped.soft.record <- function(candidate, base,
                                                           direction,
                                                           confidence,
                                                           drop.reason,
                                                           embedding.backend = NA_character_,
                                                           disk.size = NA_integer_) {
    if (is.null(confidence)) {
        confidence <- list(
            n.match.candidates = NA_integer_,
            best.score.local.rank = NA_real_,
            best.score.local.quantile = NA_real_,
            best.score.robust.z = NA_real_,
            margin.local.quantile = NA_real_,
            margin.robust.z = NA_real_,
            best.match.score = NA_real_,
            second.match.score = NA_real_,
            best.match.margin = NA_real_,
            best.match.angle = NA_real_,
            best.length.relative.error = NA_real_,
            best.transport.weight = NA_real_,
            transport.entropy = NA_real_,
            effective.matches = NA_real_,
            effective.match.fraction = NA_real_,
            reciprocal.best.match = NA
        )
    }
    data.frame(
        candidate = candidate,
        base.dart = base$dart,
        base.from = base$from,
        base.to = base$to,
        base.length = base$length,
        direction.dart = direction$dart,
        direction.from = direction$from,
        direction.to = direction$to,
        direction.length = direction$length,
        direction.label = direction$direction.label,
        match.status = "dropped_soft_confidence",
        drop.reason = drop.reason,
        n.match.candidates = confidence$n.match.candidates,
        best.score.local.rank = confidence$best.score.local.rank,
        best.score.local.quantile = confidence$best.score.local.quantile,
        best.score.robust.z = confidence$best.score.robust.z,
        margin.local.quantile = confidence$margin.local.quantile,
        margin.robust.z = confidence$margin.robust.z,
        best.match.score = confidence$best.match.score,
        second.match.score = confidence$second.match.score,
        best.match.margin = confidence$best.match.margin,
        best.match.angle = confidence$best.match.angle,
        best.length.relative.error = confidence$best.length.relative.error,
        best.transport.weight = confidence$best.transport.weight,
        transport.entropy = confidence$transport.entropy,
        effective.matches = confidence$effective.matches,
        effective.match.fraction = confidence$effective.match.fraction,
        reciprocal.best.match = confidence$reciprocal.best.match,
        embedding.backend = embedding.backend,
        disk.size = disk.size,
        stringsAsFactors = FALSE
    )
}

.transported.graph.hessian.soft.scores.from.vectors <- function(direction.vec,
                                                                direction.norm,
                                                                direction.length,
                                                                target.ids,
                                                                target.vecs,
                                                                target.norms,
                                                                darts,
                                                                soft.args) {
    cos.theta <- as.vector(target.vecs %*% direction.vec) /
        (target.norms * direction.norm)
    cos.theta <- pmin(1, pmax(-1, cos.theta))
    angle <- acos(cos.theta)
    length.difference <- abs(darts$length[target.ids] - direction.length)

    angle.scale <- soft.args$angle.scale
    if (is.null(angle.scale)) {
        positive.angle <- angle[is.finite(angle) & angle > sqrt(.Machine$double.eps)]
        angle.scale <- if (length(positive.angle)) stats::median(positive.angle) else 1
    }
    length.scale <- soft.args$length.scale
    if (is.null(length.scale)) {
        positive.length <- length.difference[
            is.finite(length.difference) &
                length.difference > sqrt(.Machine$double.eps)
        ]
        length.scale <- if (length(positive.length)) {
            stats::median(positive.length)
        } else {
            stats::median(darts$length[target.ids])
        }
    }
    if (!is.finite(angle.scale) || angle.scale <= 0) angle.scale <- 1
    if (!is.finite(length.scale) || length.scale <= 0) length.scale <- 1

    list(
        angle = angle,
        length.difference = length.difference,
        total = sqrt((angle / angle.scale)^2 +
                     (length.difference / length.scale)^2)
    )
}

.transported.graph.hessian.local.disk <- function(graph, seeds, required,
                                                  hops, max.vertices) {
    seeds <- unique(as.integer(seeds))
    required <- unique(as.integer(required))
    if (any(seeds < 1L | seeds > graph$n.vertices) ||
        any(required < 1L | required > graph$n.vertices)) {
        stop("Internal error: local disk vertices are out of graph range.",
             call. = FALSE)
    }

    hop.distance <- rep(Inf, graph$n.vertices)
    hop.distance[seeds] <- 0
    frontier <- seeds
    if (hops > 0L) {
        for (hh in seq_len(hops)) {
            nbrs <- unique(unlist(graph$adj.list[frontier], use.names = FALSE))
            nbrs <- nbrs[is.infinite(hop.distance[nbrs])]
            if (!length(nbrs)) break
            hop.distance[nbrs] <- hh
            frontier <- nbrs
        }
    }

    in.disk <- which(hop.distance <= hops)
    vertices <- unique(c(required, in.disk))
    if (length(vertices) > max.vertices && length(required) < length(vertices)) {
        optional <- setdiff(vertices, required)
        optional.order <- order(hop.distance[optional], optional)
        keep.optional <- head(optional[optional.order],
                              max(0L, max.vertices - length(required)))
        vertices <- unique(c(required, keep.optional))
    }
    vertices <- vertices[order(ifelse(is.finite(hop.distance[vertices]),
                                      hop.distance[vertices],
                                      hops + 1L),
                               vertices)]
    list(
        vertices = as.integer(vertices),
        hop.distance = hop.distance[vertices],
        metric.distance = rep(NA_real_, length(vertices)),
        metric.radius = NA_real_,
        disk.rule = "hops",
        disk.parameter = as.integer(hops)
    )
}

.transported.graph.hessian.metric.local.disk <- function(graph, seeds, required,
                                                         radius.fraction,
                                                         graph.diameter,
                                                         min.vertices,
                                                         max.vertices) {
    seeds <- unique(as.integer(seeds))
    required <- unique(as.integer(required))
    if (any(seeds < 1L | seeds > graph$n.vertices) ||
        any(required < 1L | required > graph$n.vertices)) {
        stop("Internal error: metric local disk vertices are out of graph range.",
             call. = FALSE)
    }
    if (!is.finite(graph.diameter) || graph.diameter <= 0) {
        graph.diameter <- max(unlist(graph$weight.list, use.names = FALSE))
    }
    radius <- as.double(radius.fraction) * graph.diameter
    metric.distance <- .transported.graph.hessian.multi.source.distances(graph, seeds)
    hop.distance <- .transported.graph.hessian.hop.distances(graph, seeds)
    in.disk <- which(metric.distance <= radius)
    vertices <- .transported.graph.hessian.metric.disk.vertices(
        required = required,
        in.disk = in.disk,
        metric.distance = metric.distance,
        hop.distance = hop.distance,
        min.vertices = min.vertices,
        max.vertices = max.vertices
    )
    list(
        vertices = as.integer(vertices),
        hop.distance = hop.distance[vertices],
        metric.distance = metric.distance[vertices],
        metric.radius = radius,
        local.scale = NA_real_,
        disk.rule = "metric.diameter.fraction",
        disk.parameter = as.double(radius.fraction)
    )
}

.transported.graph.hessian.metric.local.scale.disk <- function(graph, seeds, required,
                                                               local.scales,
                                                               multiplier,
                                                               min.vertices,
                                                               max.vertices) {
    seeds <- unique(as.integer(seeds))
    required <- unique(as.integer(required))
    if (any(seeds < 1L | seeds > graph$n.vertices) ||
        any(required < 1L | required > graph$n.vertices)) {
        stop("Internal error: local-scale metric disk vertices are out of graph range.",
             call. = FALSE)
    }
    seed.scales <- local.scales[seeds]
    seed.scales <- seed.scales[is.finite(seed.scales) & seed.scales > 0]
    local.scale <- if (length(seed.scales)) max(seed.scales) else {
        max(unlist(graph$weight.list, use.names = FALSE), na.rm = TRUE)
    }
    if (!is.finite(local.scale) || local.scale <= 0) {
        local.scale <- 1
    }
    radius <- as.double(multiplier) * local.scale
    metric.distance <- .transported.graph.hessian.multi.source.distances(graph, seeds)
    hop.distance <- .transported.graph.hessian.hop.distances(graph, seeds)
    in.disk <- which(metric.distance <= radius)
    vertices <- .transported.graph.hessian.metric.disk.vertices(
        required = required,
        in.disk = in.disk,
        metric.distance = metric.distance,
        hop.distance = hop.distance,
        min.vertices = min.vertices,
        max.vertices = max.vertices
    )
    list(
        vertices = as.integer(vertices),
        hop.distance = hop.distance[vertices],
        metric.distance = metric.distance[vertices],
        metric.radius = radius,
        local.scale = local.scale,
        disk.rule = "metric.local.scale",
        disk.parameter = as.double(multiplier)
    )
}

.transported.graph.hessian.metric.disk.vertices <- function(required,
                                                            in.disk,
                                                            metric.distance,
                                                            hop.distance,
                                                            min.vertices,
                                                            max.vertices) {
    vertices <- unique(c(required, in.disk))
    target <- min(max.vertices, max(length(required), as.integer(min.vertices)))
    if (length(vertices) < target) {
        nearest <- order(metric.distance, hop.distance, seq_along(metric.distance))
        nearest <- nearest[is.finite(metric.distance[nearest])]
        vertices <- unique(c(vertices, head(nearest, target)))
    }
    if (length(vertices) > max.vertices && length(required) < length(vertices)) {
        optional <- setdiff(vertices, required)
        optional.order <- order(metric.distance[optional],
                                hop.distance[optional],
                                optional)
        keep.optional <- head(optional[optional.order],
                              max(0L, max.vertices - length(required)))
        vertices <- unique(c(required, keep.optional))
    }
    vertices[order(metric.distance[vertices],
                   hop.distance[vertices],
                   vertices)]
}

.transported.graph.hessian.hop.distances <- function(graph, seeds) {
    hop.distance <- rep(Inf, graph$n.vertices)
    hop.distance[seeds] <- 0
    frontier <- unique(as.integer(seeds))
    step <- 0L
    while (length(frontier)) {
        step <- step + 1L
        nbrs <- unique(unlist(graph$adj.list[frontier], use.names = FALSE))
        nbrs <- nbrs[is.infinite(hop.distance[nbrs])]
        if (!length(nbrs)) break
        hop.distance[nbrs] <- step
        frontier <- nbrs
    }
    hop.distance
}

.transported.graph.hessian.multi.source.distances <- function(graph, seeds) {
    n <- graph$n.vertices
    distance <- rep(Inf, n)
    visited <- rep(FALSE, n)
    distance[seeds] <- 0
    repeat {
        available <- which(!visited)
        if (!length(available)) break
        idx <- available[which.min(distance[available])]
        if (!is.finite(distance[idx])) break
        visited[idx] <- TRUE
        nbrs <- graph$adj.list[[idx]]
        if (!length(nbrs)) next
        alt <- distance[idx] + graph$weight.list[[idx]]
        improve <- alt < distance[nbrs]
        if (any(improve)) {
            distance[nbrs[improve]] <- alt[improve]
        }
    }
    distance
}

.transported.graph.hessian.local.scales <- function(graph,
                                                    method,
                                                    k,
                                                    prob) {
    edge.lengths <- unlist(graph$weight.list, use.names = FALSE)
    fallback <- stats::median(edge.lengths[is.finite(edge.lengths) &
                                           edge.lengths > 0],
                              na.rm = TRUE)
    if (!is.finite(fallback) || fallback <= 0) fallback <- 1

    if (identical(method, "median.incident.length") ||
        identical(method, "quantile.incident.length")) {
        scales <- vapply(seq_len(graph$n.vertices), function(v) {
            w <- graph$weight.list[[v]]
            w <- w[is.finite(w) & w > 0]
            if (!length(w)) return(fallback)
            if (identical(method, "median.incident.length")) {
                stats::median(w)
            } else {
                as.double(stats::quantile(w, probs = prob, names = FALSE,
                                          type = 7))
            }
        }, numeric(1L))
        scales[!is.finite(scales) | scales <= 0] <- fallback
        return(scales)
    }

    if (!identical(method, "knn.distance")) {
        stop("Internal error: unknown local scale method.", call. = FALSE)
    }

    scales <- vapply(seq_len(graph$n.vertices), function(v) {
        distance <- .transported.graph.hessian.multi.source.distances(graph, v)
        distance <- sort(distance[is.finite(distance) & distance > 0])
        if (!length(distance)) return(fallback)
        distance[min(length(distance), as.integer(k))]
    }, numeric(1L))
    scales[!is.finite(scales) | scales <= 0] <- fallback
    scales
}

.transported.graph.hessian.graph.diameter <- function(graph) {
    max.distance <- 0
    for (seed in seq_len(graph$n.vertices)) {
        distance <- .transported.graph.hessian.multi.source.distances(graph, seed)
        finite <- distance[is.finite(distance)]
        if (length(finite)) {
            max.distance <- max(max.distance, finite)
        }
    }
    if (is.finite(max.distance) && max.distance > 0) {
        max.distance
    } else {
        1
    }
}

.transported.graph.hessian.local.embedding <- function(graph, vertices, method,
                                                       dim) {
    subgraph <- .transported.graph.hessian.local.subgraph(graph, vertices)
    distances <- .transported.graph.hessian.local.distance.matrix(
        subgraph$adj.list,
        subgraph$weight.list
    )

    if (identical(method, "cmdscale")) {
        coords <- .transported.graph.hessian.cmdscale.embedding(distances, dim)
        return(.transported.graph.hessian.embedding.result(
            coordinates = coords,
            backend.used = "cmdscale",
            fallback.reason = NA_character_,
            subgraph = subgraph,
            distances = distances
        ))
    }

    if (identical(method, "mds.edge.kk")) {
        return(.transported.graph.hessian.mds.edge.kk.embedding(
            subgraph = subgraph,
            distances = distances,
            dim = dim,
            fallback.reason = NA_character_
        ))
    }

    .transported.graph.hessian.grip.edge.kk.embedding(
        subgraph = subgraph,
        distances = distances,
        dim = dim
    )
}

.transported.graph.hessian.local.subgraph <- function(graph, vertices) {
    vertices <- as.integer(vertices)
    local.index <- match(seq_len(graph$n.vertices), vertices)
    adj.list <- vector("list", length(vertices))
    weight.list <- vector("list", length(vertices))
    for (idx in seq_along(vertices)) {
        v <- vertices[idx]
        keep <- !is.na(local.index[graph$adj.list[[v]]])
        adj.list[[idx]] <- as.integer(local.index[graph$adj.list[[v]][keep]])
        weight.list[[idx]] <- as.double(graph$weight.list[[v]][keep])
    }
    list(
        vertices = vertices,
        adj.list = adj.list,
        weight.list = weight.list,
        n.vertices = length(vertices)
    )
}

.transported.graph.hessian.local.distance.matrix <- function(adj.list,
                                                             weight.list) {
    n <- length(adj.list)
    distances <- matrix(Inf, n, n)
    diag(distances) <- 0
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        if (!length(nbrs)) next
        distances[i, nbrs] <- pmin(distances[i, nbrs], weight.list[[i]])
    }
    if (n >= 2L) {
        for (kk in seq_len(n)) {
            via.k <- outer(distances[, kk], distances[kk, ], "+")
            distances <- pmin(distances, via.k)
        }
    }
    distances
}

.transported.graph.hessian.cmdscale.embedding <- function(distances, dim) {
    n <- nrow(distances)
    dim <- min(as.integer(dim), max(1L, n - 1L))
    d <- distances
    if (any(!is.finite(d))) {
        finite <- d[is.finite(d) & d > 0]
        fill <- if (length(finite)) 2 * max(finite) else 1
        d[!is.finite(d)] <- fill
    }
    fit <- try(stats::cmdscale(stats::as.dist(d), k = dim, eig = FALSE, add = TRUE),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
        coords <- matrix(0, nrow = n, ncol = dim)
    } else if (is.list(fit)) {
        coords <- as.matrix(fit$points)
    } else {
        coords <- as.matrix(fit)
    }
    .transported.graph.hessian.pad.embedding(coords, dim)
}

.transported.graph.hessian.pad.embedding <- function(coords, dim) {
    coords <- as.matrix(coords)
    storage.mode(coords) <- "double"
    if (ncol(coords) < dim) {
        coords <- cbind(coords, matrix(0, nrow(coords), dim - ncol(coords)))
    }
    coords[, seq_len(dim), drop = FALSE]
}

.transported.graph.hessian.grip.edge.kk.embedding <- function(subgraph,
                                                              distances,
                                                              dim) {
    edge.kk <- .transported.graph.hessian.edge.kk.function()
    if (is.null(edge.kk)) {
        return(.transported.graph.hessian.mds.edge.kk.embedding(
            subgraph = subgraph,
            distances = distances,
            dim = dim,
            fallback.reason = "true grip edge-KK optimizer unavailable"
        ))
    }

    edge.dim <- max(2L, as.integer(dim))
    ns <- getNamespace("grip")
    layout.weighted <- get("grip.layout.weighted", envir = ns, inherits = FALSE)
    init <- try(layout.weighted(adj_list = subgraph$adj.list,
                                weight_list = subgraph$weight.list,
                                dim = edge.dim,
                                rounds = 20L,
                                final_rounds = 20L,
                                num_init = 4L,
                                seed = 1L,
                                disconnected = "error"),
                silent = TRUE)
    if (!inherits(init, "try-error")) {
        opt <- .transported.graph.hessian.run.edge.kk(
            edge.kk = edge.kk,
            coords = init,
            subgraph = subgraph,
            dim = edge.dim
        )
        if (!inherits(opt, "try-error") && is.list(opt) && !is.null(opt$coords)) {
            return(.transported.graph.hessian.embedding.result(
                coordinates = .transported.graph.hessian.pad.embedding(opt$coords, dim),
                backend.used = paste0("grip.layout.weighted+", opt$backend),
                fallback.reason = NA_character_,
                subgraph = subgraph,
                distances = distances
            ))
        }
    }

    .transported.graph.hessian.mds.edge.kk.embedding(
        subgraph = subgraph,
        distances = distances,
        dim = dim,
        fallback.reason = "grip.layout.weighted or edge-KK failed"
    )
}

.transported.graph.hessian.mds.edge.kk.embedding <- function(subgraph,
                                                             distances,
                                                             dim,
                                                             fallback.reason) {
    init <- .transported.graph.hessian.cmdscale.embedding(distances, dim)
    edge.kk <- .transported.graph.hessian.edge.kk.function()
    if (is.null(edge.kk)) {
        return(.transported.graph.hessian.embedding.result(
            coordinates = init,
            backend.used = "cmdscale",
            fallback.reason = fallback.reason %||% "true grip edge-KK optimizer unavailable",
            subgraph = subgraph,
            distances = distances
        ))
    }

    edge.dim <- max(2L, as.integer(dim))
    init.edge <- .transported.graph.hessian.pad.embedding(init, edge.dim)
    opt <- .transported.graph.hessian.run.edge.kk(
        edge.kk = edge.kk,
        coords = init.edge,
        subgraph = subgraph,
        dim = edge.dim
    )
    if (inherits(opt, "try-error") || !is.list(opt) || is.null(opt$coords)) {
        return(.transported.graph.hessian.embedding.result(
            coordinates = init,
            backend.used = "cmdscale",
            fallback.reason = fallback.reason %||% "grip edge-KK optimizer failed",
            subgraph = subgraph,
            distances = distances
        ))
    }

    .transported.graph.hessian.embedding.result(
        coordinates = .transported.graph.hessian.pad.embedding(opt$coords, dim),
        backend.used = paste0("cmdscale+", opt$backend),
        fallback.reason = fallback.reason,
        subgraph = subgraph,
        distances = distances
    )
}

.transported.graph.hessian.edge.kk.function <- function() {
    if (!requireNamespace("grip", quietly = TRUE)) return(NULL)
    ns <- getNamespace("grip")
    for (nm in c("grip.optimize.edge.kk.layout",
                 "grip.optimize.edge.isometric.layout")) {
        if (exists(nm, envir = ns, inherits = FALSE)) {
            return(list(
                name = nm,
                fun = get(nm, envir = ns, inherits = FALSE)
            ))
        }
    }
    NULL
}

.transported.graph.hessian.run.edge.kk <- function(edge.kk, coords, subgraph,
                                                   dim) {
    args <- list(
        coords = coords,
        adj_list = subgraph$adj.list,
        weight_list = subgraph$weight.list,
        dim = dim,
        stiffness_method = "uniform",
        density_mix_schedule = 1,
        scale_mode = "identity",
        max_iter = 16L,
        return_trace = FALSE,
        diagnostics = FALSE
    )
    formals.names <- names(formals(edge.kk$fun))
    if ("engine" %in% formals.names) {
        args$engine <- "cpp"
    }
    call.args <- if ("..." %in% formals.names) {
        args
    } else {
        args[names(args) %in% formals.names]
    }
    out <- try(do.call(edge.kk$fun, call.args), silent = TRUE)
    if (!inherits(out, "try-error") && is.list(out)) {
        out$backend <- edge.kk$name
    }
    out
}

.transported.graph.hessian.embedding.result <- function(coordinates,
                                                        backend.used,
                                                        fallback.reason,
                                                        subgraph,
                                                        distances) {
    coordinates <- as.matrix(coordinates)
    list(
        coordinates = coordinates,
        backend.used = backend.used,
        fallback.reason = fallback.reason %||% NA_character_,
        edge.stress = .transported.graph.hessian.edge.stress(
            coordinates,
            subgraph
        ),
        graph.distance.stress = .transported.graph.hessian.distance.stress(
            coordinates,
            distances
        ),
        distances = distances
    )
}

.transported.graph.hessian.edge.stress <- function(coordinates, subgraph) {
    residuals <- numeric()
    weights <- numeric()
    for (i in seq_along(subgraph$adj.list)) {
        nbrs <- subgraph$adj.list[[i]]
        keep <- nbrs > i
        if (!any(keep)) next
        nbrs <- nbrs[keep]
        wts <- subgraph$weight.list[[i]][keep]
        edist <- sqrt(rowSums((coordinates[nbrs, , drop = FALSE] -
                               matrix(coordinates[i, ],
                                      nrow = length(nbrs),
                                      ncol = ncol(coordinates),
                                      byrow = TRUE))^2))
        residuals <- c(residuals, edist - wts)
        weights <- c(weights, wts)
    }
    if (!length(residuals)) return(NA_real_)
    denom <- mean(weights)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    sqrt(mean(residuals^2)) / denom
}

.transported.graph.hessian.distance.stress <- function(coordinates, distances) {
    idx <- upper.tri(distances) & is.finite(distances) & distances > 0
    if (!any(idx)) return(NA_real_)
    edist <- as.matrix(stats::dist(coordinates))
    denom <- mean(distances[idx])
    if (!is.finite(denom) || denom <= 0) denom <- 1
    sqrt(mean((edist[idx] - distances[idx])^2)) / denom
}

.transported.graph.hessian.soft.scales <- function(darts, dart.vectors,
                                                   dart.norms, soft.args) {
    angles <- numeric()
    length.diffs <- numeric()
    for (base.idx in seq_len(nrow(darts))) {
        local.ids <- which(darts$from == darts$from[base.idx])
        target.ids <- which(darts$from == darts$to[base.idx])
        if (length(target.ids) < 1L) next
        for (direction.id in local.ids) {
            scores <- .transported.graph.hessian.raw.soft.scores(
                direction.id = direction.id,
                target.ids = target.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms
            )
            angles <- c(angles, scores$angle)
            length.diffs <- c(length.diffs, scores$length.difference)
        }
    }
    positive.angle <- angles[is.finite(angles) & angles > sqrt(.Machine$double.eps)]
    positive.length <- length.diffs[is.finite(length.diffs) &
                                        length.diffs > sqrt(.Machine$double.eps)]
    angle.scale <- soft.args$angle.scale
    if (is.null(angle.scale)) {
        angle.scale <- if (length(positive.angle)) stats::median(positive.angle) else 1
    }
    length.scale <- soft.args$length.scale
    if (is.null(length.scale)) {
        length.scale <- if (length(positive.length)) {
            stats::median(positive.length)
        } else {
            stats::median(darts$length)
        }
    }
    if (!is.finite(angle.scale) || angle.scale <= 0) angle.scale <- 1
    if (!is.finite(length.scale) || length.scale <= 0) length.scale <- 1
    list(
        angle.scale = as.double(angle.scale),
        length.scale = as.double(length.scale),
        bandwidth = soft.args$bandwidth,
        angle.scale.source = if (is.null(soft.args$angle.scale)) "estimated.global.median" else "user",
        length.scale.source = if (is.null(soft.args$length.scale)) "estimated.global.median" else "user"
    )
}

.transported.graph.hessian.raw.soft.scores <- function(direction.id, target.ids,
                                                       darts, dart.vectors,
                                                       dart.norms) {
    direction.vec <- dart.vectors[direction.id, , drop = TRUE]
    target.vecs <- dart.vectors[target.ids, , drop = FALSE]
    cos.theta <- as.vector(target.vecs %*% direction.vec) /
        (dart.norms[target.ids] * dart.norms[direction.id])
    cos.theta <- pmin(1, pmax(-1, cos.theta))
    angle <- acos(cos.theta)
    length.difference <- abs(darts$length[target.ids] - darts$length[direction.id])
    list(angle = angle, length.difference = length.difference)
}

.transported.graph.hessian.soft.scores <- function(direction.id, target.ids,
                                                   darts, dart.vectors,
                                                   dart.norms, scales) {
    raw <- .transported.graph.hessian.raw.soft.scores(
        direction.id = direction.id,
        target.ids = target.ids,
        darts = darts,
        dart.vectors = dart.vectors,
        dart.norms = dart.norms
    )
    total <- sqrt((raw$angle / scales$angle.scale)^2 +
                  (raw$length.difference / scales$length.scale)^2)
    list(
        angle = raw$angle,
        length.difference = raw$length.difference,
        total = total
    )
}

.transported.graph.hessian.edge.angle.scales <- function(darts,
                                                         dart.vectors,
                                                         dart.norms,
                                                         edge.angle.args) {
    angles <- numeric()
    length.diffs <- numeric()
    for (base.idx in seq_len(nrow(darts))) {
        local.ids <- which(darts$from == darts$from[base.idx])
        target.ids <- which(darts$from == darts$to[base.idx])
        if (length(target.ids) < 1L) next
        for (direction.id in local.ids) {
            raw <- .transported.graph.hessian.raw.edge.angle.scores(
                base.idx = base.idx,
                direction.id = direction.id,
                target.ids = target.ids,
                darts = darts,
                dart.vectors = dart.vectors,
                dart.norms = dart.norms
            )
            angles <- c(angles, raw$angle)
            length.diffs <- c(length.diffs, raw$length.difference)
        }
    }
    positive.angle <- angles[is.finite(angles) & angles > sqrt(.Machine$double.eps)]
    positive.length <- length.diffs[is.finite(length.diffs) &
                                        length.diffs > sqrt(.Machine$double.eps)]
    angle.scale <- edge.angle.args$angle.scale
    if (is.null(angle.scale)) {
        angle.scale <- if (length(positive.angle)) stats::median(positive.angle) else 1
    }
    length.scale <- edge.angle.args$length.scale
    if (is.null(length.scale)) {
        length.scale <- if (length(positive.length)) {
            stats::median(positive.length)
        } else {
            stats::median(darts$length)
        }
    }
    if (!is.finite(angle.scale) || angle.scale <= 0) angle.scale <- 1
    if (!is.finite(length.scale) || length.scale <= 0) length.scale <- 1
    list(
        angle.scale = as.double(angle.scale),
        length.scale = as.double(length.scale),
        bandwidth = edge.angle.args$bandwidth,
        angle.scale.source = if (is.null(edge.angle.args$angle.scale)) {
            "estimated.global.median"
        } else {
            "user"
        },
        length.scale.source = if (is.null(edge.angle.args$length.scale)) {
            "estimated.global.median"
        } else {
            "user"
        }
    )
}

.transported.graph.hessian.raw.edge.angle.scores <- function(base.idx,
                                                            direction.id,
                                                            target.ids,
                                                            darts,
                                                            dart.vectors,
                                                            dart.norms) {
    base.vec <- dart.vectors[base.idx, , drop = TRUE]
    base.norm <- dart.norms[base.idx]
    direction.vec <- dart.vectors[direction.id, , drop = TRUE]
    direction.norm <- dart.norms[direction.id]
    source.angle <- .transported.graph.hessian.vector.angle(
        base.vec,
        direction.vec,
        base.norm,
        direction.norm
    )
    target.vecs <- dart.vectors[target.ids, , drop = FALSE]
    target.norms <- dart.norms[target.ids]
    target.angle <- vapply(seq_along(target.ids), function(idx) {
        .transported.graph.hessian.vector.angle(
            base.vec,
            target.vecs[idx, , drop = TRUE],
            base.norm,
            target.norms[idx]
        )
    }, numeric(1L))
    angle <- abs(target.angle - source.angle)
    angle <- pmin(angle, 2 * pi - angle)
    length.difference <- abs(darts$length[target.ids] - darts$length[direction.id])
    length.relative.error <- length.difference /
        max(abs(darts$length[direction.id]), sqrt(.Machine$double.eps))
    list(
        source.angle = rep(source.angle, length(target.ids)),
        target.angle = target.angle,
        angle = angle,
        length.difference = length.difference,
        length.relative.error = length.relative.error
    )
}

.transported.graph.hessian.edge.angle.scores <- function(base.idx,
                                                         direction.id,
                                                         target.ids,
                                                         darts,
                                                         dart.vectors,
                                                         dart.norms,
                                                         scales) {
    raw <- .transported.graph.hessian.raw.edge.angle.scores(
        base.idx = base.idx,
        direction.id = direction.id,
        target.ids = target.ids,
        darts = darts,
        dart.vectors = dart.vectors,
        dart.norms = dart.norms
    )
    total <- sqrt((raw$angle / scales$angle.scale)^2 +
                  (raw$length.difference / scales$length.scale)^2)
    c(raw, list(total = total))
}

.transported.graph.hessian.vector.angle <- function(a, b, a.norm, b.norm) {
    if (!is.finite(a.norm) || !is.finite(b.norm) || a.norm <= 0 || b.norm <= 0) {
        return(NA_real_)
    }
    cos.theta <- sum(a * b) / (a.norm * b.norm)
    acos(pmin(1, pmax(-1, cos.theta)))
}

.transported.graph.hessian.softmax <- function(score, bandwidth) {
    z <- -(score^2) / (bandwidth^2)
    z <- z - max(z)
    w <- exp(z)
    w / sum(w)
}

.transported.graph.hessian.matrix <- function(row.table, n.vertices) {
    if (!nrow(row.table)) {
        return(Matrix::sparseMatrix(
            i = integer(),
            j = integer(),
            x = numeric(),
            dims = c(0L, n.vertices),
            giveCsparse = TRUE
        ))
    }
    if (all(c("coefficient.vertex", "coefficient.value") %in% names(row.table))) {
        return(Matrix::sparseMatrix(
            i = row.table$row,
            j = as.integer(row.table$coefficient.vertex),
            x = as.double(row.table$coefficient.value),
            dims = c(max(row.table$row), n.vertices),
            giveCsparse = TRUE
        ))
    }
    n.rows <- max(row.table$row)
    ii <- rep(row.table$row, each = 4L)
    jj <- as.integer(rbind(
        row.table$matched.to,
        row.table$matched.from,
        row.table$direction.to,
        row.table$direction.from
    ))
    weight <- row.table$transport.weight
    xx <- as.double(rbind(
        weight / row.table$matched.length,
        -weight / row.table$matched.length,
        -weight / row.table$direction.length,
        weight / row.table$direction.length
    ))
    Matrix::sparseMatrix(
        i = ii,
        j = jj,
        x = xx,
        dims = c(n.rows, n.vertices),
        giveCsparse = TRUE
    )
}

.transported.graph.hessian.probe.residuals <- function(hessian, polynomial.probes,
                                                       n.vertices) {
    if (is.null(polynomial.probes)) return(NULL)
    probes <- as.matrix(polynomial.probes)
    storage.mode(probes) <- "double"
    if (nrow(probes) != n.vertices || any(!is.finite(probes))) {
        stop("polynomial.probes must be a finite numeric vector or matrix with one row per vertex.",
             call. = FALSE)
    }
    if (is.null(colnames(probes))) {
        colnames(probes) <- paste0("probe", seq_len(ncol(probes)))
    }
    applied <- as.matrix(hessian %*% probes)
    numerator <- sqrt(sum(applied^2))
    denominator <- sqrt(sum(probes^2))
    per.column <- data.frame(
        probe = colnames(probes),
        residual = vapply(seq_len(ncol(probes)), function(idx) {
            den <- sqrt(sum(probes[, idx]^2))
            if (den == 0) return(NA_real_)
            sqrt(sum(applied[, idx]^2)) / den
        }, numeric(1)),
        stringsAsFactors = FALSE
    )
    list(
        overall = if (denominator == 0) NA_real_ else numerator / denominator,
        per.column = per.column
    )
}

.transported.graph.hessian.summary <- function(darts, transport) {
    row.table <- transport$row.table
    dropped.table <- transport$dropped.table
    n.rows <- if (nrow(row.table)) max(row.table$row) else 0L
    n.row.records <- nrow(row.table)
    n.dropped <- nrow(dropped.table)
    n.ambiguous <- if (n.row.records) {
        sum(row.table$match.status == "matched_ambiguous")
    } else {
        0L
    }
    mean.effective.matches <- if (n.row.records && "effective.matches" %in% names(row.table)) {
        mean(row.table$effective.matches[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.effective.match.fraction <- if (n.row.records &&
                                         "effective.match.fraction" %in% names(row.table)) {
        mean(row.table$effective.match.fraction[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.best.score.robust.z <- if (n.row.records &&
                                    "best.score.robust.z" %in% names(row.table)) {
        mean(row.table$best.score.robust.z[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.margin.robust.z <- if (n.row.records &&
                                "margin.robust.z" %in% names(row.table)) {
        mean(row.table$margin.robust.z[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.best.score.local.quantile <- if (n.row.records &&
                                          "best.score.local.quantile" %in% names(row.table)) {
        mean(row.table$best.score.local.quantile[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.margin.local.quantile <- if (n.row.records &&
                                      "margin.local.quantile" %in% names(row.table)) {
        mean(row.table$margin.local.quantile[!duplicated(row.table$row)])
    } else {
        NA_real_
    }
    mean.match.score <- if (n.row.records && "match.score" %in% names(row.table)) {
        mean(row.table$match.score)
    } else {
        NA_real_
    }
    data.frame(
        n.darts = nrow(darts),
        n.candidate.rows = transport$n.candidates,
        n.rows = n.rows,
        n.row.records = n.row.records,
        n.dropped = n.dropped,
        n.ambiguous = n.ambiguous,
        n.direction.labels = length(unique(darts$direction.label)),
        mean.transport.entropy = if (n.row.records) mean(row.table$transport.entropy) else NA_real_,
        mean.effective.matches = mean.effective.matches,
        mean.effective.match.fraction = mean.effective.match.fraction,
        mean.best.score.robust.z = mean.best.score.robust.z,
        mean.margin.robust.z = mean.margin.robust.z,
        mean.best.score.local.quantile = mean.best.score.local.quantile,
        mean.margin.local.quantile = mean.margin.local.quantile,
        mean.match.score = mean.match.score,
        stringsAsFactors = FALSE
    )
}

#' @export
print.transported.graph.hessian.operator <- function(x, ...) {
    cat("\nTransported Graph Hessian Operator\n")
    cat("==================================\n\n")
    cat(sprintf("Vertices: %d; darts: %d\n", x$graph$n.vertices, x$graph$n.darts))
    cat(sprintf("Transport rule: %s\n", x$transport.rule))
    summary <- x$diagnostics$summary
    cat(sprintf("Rows: %d; dropped candidates: %d\n",
                summary$n.rows, summary$n.dropped))
    cat(sprintf("Nullity estimate: %s\n",
                if (is.na(x$diagnostics$nullity.estimate)) {
                    "NA"
                } else {
                    as.character(x$diagnostics$nullity.estimate)
                }))
    invisible(x)
}

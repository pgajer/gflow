
#' Detect Graph Endpoints from a 3D Embedding
#'
#' @description
#' Detects terminal arm tips in a graph by:
#' \enumerate{
#'   \item computing graph endpoint scores from the 3D embedding,
#'   \item optionally smoothing those scores with the archived `gflowx`
#'     rdgraph estimator,
#'   \item finding graph-local maxima with [detect.local.extrema()],
#'   \item retaining only maxima that persist across neighborhood scales.
#' }
#'
#' The function works with terminal regions in an embedding; returned endpoints
#' need not be degree-1 vertices.
#'
#' @param adj.list Graph adjacency list (1-based vertex indices).
#' @param weight.list Edge-length list aligned with `adj.list`.
#' @param layout.3d Numeric matrix or data frame with 3 columns giving the 3D
#'   embedding coordinates for graph vertices.
#' @param neighborhood Character string. Either `"geodesic_k"` (default) or
#'   `"geodesic_radius"`.
#' @param k Integer vector of geodesic neighborhood sizes used when
#'   `neighborhood = "geodesic_k"`.
#' @param radius Positive numeric vector of geodesic radii used when
#'   `neighborhood = "geodesic_radius"`.
#' @param q Quantile used for `s.q`. Must lie in `[0, 1]`. Default is `0.1`.
#' @param neighbor.weighting Character string controlling neighbor weights for
#'   `m` and weighted pair quantiles. One of `"uniform"` (default),
#'   `"inverse_distance"`, or `"gaussian"`.
#' @param gaussian.sigma Optional positive numeric scale used when
#'   `neighbor.weighting = "gaussian"`.
#' @param min.neighborhood.size Minimum number of usable neighbors required after
#'   removing zero-length embedding vectors.
#' @param scale.aggregate How to aggregate scores across multiple neighborhood
#'   scales. One of `"mean"` (default), `"median"`, `"min"`, or `"max"`.
#' @param score.metric Which score to maximize. One of `"score"` (default),
#'   `"m"`, `"s.q"`, or `"s.min"`.
#' @param smooth Logical; if `TRUE`, smooths aggregated scores on the graph
#'   before detecting extrema.
#' @param fitted.model Optional fitted graph smoother from
#'   `gflowx::fit.rdgraph.regression()`. If supplied, smoothing reuses its
#'   spectral decomposition through `gflowx`.
#' @param smooth.fit.args Optional named list of arguments passed to
#'   `gflowx::fit.rdgraph.regression()` when `smooth = TRUE` and `fitted.model` is
#'   `NULL`. `X`, `y`, `adj.list`, `weight.list`, and `verbose.level` are filled
#'   automatically when not supplied.
#' @param smooth.refit.args Optional named list of arguments passed to
#'   `gflowx::refit.rdgraph.regression()` when smoothing aggregated score fields.
#' @param detect.max.radius Maximum radius passed to [detect.local.extrema()]
#'   when calling endpoint candidates from the final detection score.
#' @param detect.min.neighborhood.size Minimum neighborhood size passed to
#'   [detect.local.extrema()] when calling endpoint candidates.
#' @param scale.stability.radius Optional non-negative geodesic radius used when
#'   translating per-scale local maxima into multiscale support. If `NULL`
#'   (default), uses `detect.max.radius`. This allows smoothed maxima that shift
#'   slightly inward to remain supported by nearby raw multiscale maxima.
#' @param min.scale.stability Minimum fraction of scales for which a vertex must
#'   be a local maximum of the per-scale detection score to be retained as an
#'   endpoint. Default is `0.5`.
#' @param min.score.quantile Quantile filter applied to the final detection
#'   score before returning endpoints. Default is `0.8`, which suppresses weak
#'   local maxima while keeping high-scoring terminal tips.
#' @param max.endpoints Optional positive integer cap on returned endpoints.
#' @param verbose Logical; if `TRUE`, prints progress.
#'
#' @return A list with endpoint vertices, score diagnostics, optional smoothed
#'   scores, scale stability summaries, and the local-extrema fit used for
#'   endpoint calling.
#'
#' @export
detect.graph.endpoints <- function(adj.list,
                                   weight.list,
                                   layout.3d,
                                   neighborhood = c("geodesic_k", "geodesic_radius"),
                                   k = c(10L, 20L, 30L),
                                   radius = NULL,
                                   q = 0.10,
                                   neighbor.weighting = c("uniform", "inverse_distance", "gaussian"),
                                   gaussian.sigma = NULL,
                                   min.neighborhood.size = 3L,
                                   scale.aggregate = c("mean", "median", "min", "max"),
                                   score.metric = c("score", "m", "s.q", "s.min"),
                                   smooth = FALSE,
                                   fitted.model = NULL,
                                   smooth.fit.args = NULL,
                                   smooth.refit.args = NULL,
                                   detect.max.radius = 2,
                                   detect.min.neighborhood.size = 2L,
                                   scale.stability.radius = NULL,
                                   min.scale.stability = 0.5,
                                   min.score.quantile = 0.8,
                                   max.endpoints = NULL,
                                   verbose = FALSE) {
    neighborhood <- match.arg(neighborhood)
    neighbor.weighting <- match.arg(neighbor.weighting)
    scale.aggregate <- match.arg(scale.aggregate)
    score.metric <- match.arg(score.metric)

    if (!is.logical(smooth) || length(smooth) != 1L) {
        stop("'smooth' must be a scalar logical.")
    }
    if (!is.numeric(detect.max.radius) || length(detect.max.radius) != 1L ||
        !is.finite(detect.max.radius) || detect.max.radius <= 0) {
        stop("'detect.max.radius' must be a finite scalar > 0.")
    }
    if (!is.numeric(detect.min.neighborhood.size) ||
        length(detect.min.neighborhood.size) != 1L ||
        !is.finite(detect.min.neighborhood.size) ||
        detect.min.neighborhood.size < 1 ||
        detect.min.neighborhood.size != floor(detect.min.neighborhood.size)) {
        stop("'detect.min.neighborhood.size' must be a positive integer.")
    }
    if (!is.numeric(min.scale.stability) || length(min.scale.stability) != 1L ||
        !is.finite(min.scale.stability) ||
        min.scale.stability < 0 || min.scale.stability > 1) {
        stop("'min.scale.stability' must be a finite scalar in [0, 1].")
    }
    if (!is.null(scale.stability.radius)) {
        if (!is.numeric(scale.stability.radius) || length(scale.stability.radius) != 1L ||
            !is.finite(scale.stability.radius) || scale.stability.radius < 0) {
            stop("'scale.stability.radius' must be NULL or a finite scalar >= 0.")
        }
    }
    if (!is.null(min.score.quantile)) {
        if (!is.numeric(min.score.quantile) || length(min.score.quantile) != 1L ||
            !is.finite(min.score.quantile) ||
            min.score.quantile < 0 || min.score.quantile > 1) {
                stop("'min.score.quantile' must be NULL or a finite scalar in [0, 1].")
        }
    }
    if (!is.null(max.endpoints)) {
        if (!is.numeric(max.endpoints) || length(max.endpoints) != 1L ||
            !is.finite(max.endpoints) || max.endpoints < 1 ||
            max.endpoints != floor(max.endpoints)) {
            stop("'max.endpoints' must be NULL or a positive integer.")
        }
    }

    scores <- dgraphs::compute.graph.endpoint.scores(
        adj.list = adj.list,
        weight.list = weight.list,
        layout.3d = layout.3d,
        neighborhood = neighborhood,
        k = k,
        radius = radius,
        q = q,
        neighbor.weighting = neighbor.weighting,
        gaussian.sigma = gaussian.sigma,
        min.neighborhood.size = min.neighborhood.size,
        scale.aggregate = scale.aggregate,
        verbose = verbose
    )

    metric.map <- list(
        score = scores$score,
        m = scores$m,
        s.q = scores$s.q,
        s.min = scores$s.min
    )
    detection.score.raw <- as.numeric(metric.map[[score.metric]])

    smoothed.scores <- NULL
    smoothing.model <- fitted.model
    smoothing.refit <- NULL
    detection.score <- detection.score.raw
    stability.score.by.scale <- NULL
    stability.smoothing.refit <- NULL

    if (smooth) {
        if (!requireNamespace("gflowx", quietly = TRUE)) {
            stop(
                "smooth = TRUE requires the optional archive package 'gflowx'. ",
                "Install gflowx or provide externally smoothed endpoint scores.",
                call. = FALSE
            )
        }
        fit.rdgraph <- getExportedValue("gflowx", "fit.rdgraph.regression")
        refit.rdgraph <- getExportedValue("gflowx", "refit.rdgraph.regression")

        if (verbose) {
            message("Smoothing aggregated endpoint scores on the graph")
        }

        score.matrix <- cbind(
            s.min = scores$s.min,
            s.q = scores$s.q,
            m = scores$m,
            score = scores$score
        )

        if (is.null(smoothing.model)) {
            fit.args <- smooth.fit.args
            if (is.null(fit.args)) fit.args <- list()
            fit.args$X <- NULL
            fit.args$y <- NULL

            if (is.null(fit.args$k)) {
                fit.args$k <- max(
                    2L,
                    min(
                        30L,
                        as.integer(round(stats::median(pmax(lengths(adj.list), 1L))))
                    )
                )
            }
            if (is.null(fit.args$adj.list)) fit.args$adj.list <- adj.list
            if (is.null(fit.args$weight.list)) fit.args$weight.list <- weight.list
            if (is.null(fit.args$verbose.level)) fit.args$verbose.level <- 0L
            if (is.null(fit.args$compute.extremality)) fit.args$compute.extremality <- FALSE

            fit.call <- c(list(
                X = as.matrix(layout.3d),
                y = detection.score.raw
            ), fit.args)
            smoothing.model <- do.call(fit.rdgraph, fit.call)
        }

        refit.args <- smooth.refit.args
        if (is.null(refit.args)) refit.args <- list()
        refit.args$fitted.model <- NULL
        refit.args$y.new <- NULL
        if (is.null(refit.args$per.column.gcv)) refit.args$per.column.gcv <- TRUE
        if (is.null(refit.args$verbose)) refit.args$verbose <- FALSE

        refit.call <- c(list(
            fitted.model = smoothing.model,
            y.new = score.matrix
        ), refit.args)
        smoothing.refit <- do.call(refit.rdgraph, refit.call)

        smoothed.scores <- smoothing.refit$fitted.values
        smoothed.scores <- as.matrix(smoothed.scores)
        colnames(smoothed.scores) <- colnames(score.matrix)
        detection.score <- as.numeric(smoothed.scores[, score.metric])
    }

    by.scale.metric <- scores$by.scale[[score.metric]]
    stability.score.by.scale <- by.scale.metric
    if (smooth) {
        refit.args.scale <- refit.args
        refit.args.scale$fitted.model <- NULL
        refit.args.scale$y.new <- NULL
        refit.call.scale <- c(list(
            fitted.model = smoothing.model,
            y.new = by.scale.metric
        ), refit.args.scale)
        stability.smoothing.refit <- do.call(refit.rdgraph, refit.call.scale)
        stability.score.by.scale <- stability.smoothing.refit$fitted.values
        stability.score.by.scale <- as.matrix(stability.score.by.scale)
        colnames(stability.score.by.scale) <- colnames(by.scale.metric)
    }

    local.max.by.scale <- matrix(FALSE, nrow = nrow(by.scale.metric), ncol = ncol(by.scale.metric))
    local.max.strong.by.scale <- matrix(FALSE, nrow = nrow(by.scale.metric), ncol = ncol(by.scale.metric))
    colnames(local.max.by.scale) <- colnames(by.scale.metric)
    colnames(local.max.strong.by.scale) <- colnames(by.scale.metric)
    scale.max.threshold.prob <- if (is.null(min.score.quantile)) {
        0.90
    } else {
        max(0.90, as.numeric(min.score.quantile))
    }
    scale.max.thresholds <- rep(NA_real_, ncol(by.scale.metric))

    for (scale.idx in seq_len(ncol(by.scale.metric))) {
        y.scale <- stability.score.by.scale[, scale.idx]
        if (!any(is.finite(y.scale))) next

        y.detect <- .replace.nonfinite.endpoint.values(y.scale)
        ext.scale <- detect.local.extrema(
            adj.list = adj.list,
            weight.list = weight.list,
            y = y.detect,
            max.radius = detect.max.radius,
            min.neighborhood.size = detect.min.neighborhood.size,
            detect.maxima = TRUE,
            custom.prefix = "E"
        )
        if (length(ext.scale$vertices) > 0L) {
            local.max.by.scale[ext.scale$vertices, scale.idx] <- TRUE
            scale.max.thresholds[scale.idx] <- stats::quantile(
                y.scale[is.finite(y.scale)],
                probs = scale.max.threshold.prob,
                na.rm = TRUE
            )[[1]]
            strong.vertices <- ext.scale$vertices[
                y.scale[ext.scale$vertices] >= scale.max.thresholds[scale.idx]
            ]
            if (length(strong.vertices) > 0L) {
                local.max.strong.by.scale[strong.vertices, scale.idx] <- TRUE
            }
        }
    }

    stability.radius <- if (is.null(scale.stability.radius)) detect.max.radius else scale.stability.radius
    local.max.filtered.by.scale <- .suppress.graph.endpoint.maxima.by.scale(
        adj.list = adj.list,
        weight.list = weight.list,
        local.max.by.scale = local.max.strong.by.scale,
        score.by.scale = stability.score.by.scale,
        radius = stability.radius,
        prefer.cpp = TRUE
    )
    local.max.support.by.scale <- local.max.filtered.by.scale

    if (stability.radius > 0) {
        local.max.support.by.scale <- .compute.graph.endpoint.support.by.scale(
            adj.list = adj.list,
            weight.list = weight.list,
            local.max.by.scale = local.max.filtered.by.scale,
            radius = stability.radius,
            prefer.cpp = TRUE
        )
    }

    finite.scale.count <- rowSums(is.finite(by.scale.metric))
    scale.stability <- rowSums(local.max.support.by.scale) / pmax(finite.scale.count, 1L)
    scale.stability[finite.scale.count < 1L] <- NA_real_

    y.final <- .replace.nonfinite.endpoint.values(detection.score)
    local.maxima <- detect.local.extrema(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y.final,
        max.radius = detect.max.radius,
        min.neighborhood.size = detect.min.neighborhood.size,
        detect.maxima = TRUE,
        custom.prefix = "E"
    )

    local.max.vertices <- local.maxima$vertices
    final.threshold <- NA_real_
    if (!is.null(min.score.quantile)) {
        final.threshold <- stats::quantile(
            detection.score[is.finite(detection.score)],
            probs = min.score.quantile,
            na.rm = TRUE
        )[[1]]
    }

    endpoint.vertices <- local.max.vertices
    if (length(endpoint.vertices) > 0L) {
        keep <- rep(TRUE, length(endpoint.vertices))
        if (!is.null(min.score.quantile) && is.finite(final.threshold)) {
            keep <- keep & (detection.score[endpoint.vertices] >= final.threshold)
        }
        keep <- keep & !is.na(scale.stability[endpoint.vertices]) &
            (scale.stability[endpoint.vertices] >= min.scale.stability)
        endpoint.vertices <- endpoint.vertices[keep]
    }

    if (!is.null(max.endpoints) && length(endpoint.vertices) > max.endpoints) {
        ord <- order(detection.score[endpoint.vertices], decreasing = TRUE)
        endpoint.vertices <- endpoint.vertices[ord[seq_len(max.endpoints)]]
    }

    summary.df <- data.frame(
        vertex = seq_along(detection.score),
        s.min = scores$s.min,
        s.q = scores$s.q,
        m = scores$m,
        score = scores$score,
        detection.score.raw = detection.score.raw,
        detection.score = detection.score,
        scale.coverage = scores$scale.coverage,
        scale.stability = scale.stability,
        is.local.max = seq_along(detection.score) %in% local.max.vertices,
        is.endpoint = seq_along(detection.score) %in% endpoint.vertices,
        stringsAsFactors = FALSE
    )

    if (!is.null(smoothed.scores)) {
        summary.df$s.min.smooth <- smoothed.scores[, "s.min"]
        summary.df$s.q.smooth <- smoothed.scores[, "s.q"]
        summary.df$m.smooth <- smoothed.scores[, "m"]
        summary.df$score.smooth <- smoothed.scores[, "score"]
    }

    res <- list(
        endpoints = as.integer(endpoint.vertices),
        score.metric = score.metric,
        scores = scores,
        smooth = smooth,
        smoothing.model = smoothing.model,
        smoothing.refit = smoothing.refit,
        smoothed.scores = smoothed.scores,
        detection.score.raw = as.numeric(detection.score.raw),
        detection.score = as.numeric(detection.score),
        local.maxima = local.maxima,
        local.max.by.scale = local.max.by.scale,
        local.max.strong.by.scale = local.max.strong.by.scale,
        local.max.filtered.by.scale = local.max.filtered.by.scale,
        local.max.support.by.scale = local.max.support.by.scale,
        stability.score.by.scale = stability.score.by.scale,
        stability.smoothing.refit = stability.smoothing.refit,
        scale.stability = as.numeric(scale.stability),
        finite.scale.count = as.integer(finite.scale.count),
        scale.stability.radius = stability.radius,
        scale.max.threshold.prob = scale.max.threshold.prob,
        scale.max.thresholds = as.numeric(scale.max.thresholds),
        min.scale.stability = min.scale.stability,
        score.threshold = final.threshold,
        summary = summary.df
    )

    class(res) <- c("graph_endpoints", class(res))
    res
}

.compute.graph.endpoint.support.by.scale <- function(adj.list,
                                                     weight.list,
                                                     local.max.by.scale,
                                                     radius,
                                                     prefer.cpp = TRUE)
{
    local.max.by.scale <- as.matrix(local.max.by.scale)
    local.max.by.scale <- matrix(
        as.logical(local.max.by.scale),
        nrow = nrow(local.max.by.scale),
        ncol = ncol(local.max.by.scale),
        dimnames = dimnames(local.max.by.scale)
    )

    if (isTRUE(prefer.cpp) && exists("rcpp_graph_multi_source_support_by_scale", mode = "function")) {
        support <- rcpp_graph_multi_source_support_by_scale(
            adj_list = adj.list,
            weight_list = weight.list,
            local_max_by_scale = local.max.by.scale,
            radius = as.numeric(radius)
        )
        support <- as.matrix(support)
        dimnames(support) <- dimnames(local.max.by.scale)
        return(support)
    }

    graph.obj <- .build.graph.endpoint.igraph(adj.list, weight.list)
    support <- matrix(
        FALSE,
        nrow = nrow(local.max.by.scale),
        ncol = ncol(local.max.by.scale),
        dimnames = dimnames(local.max.by.scale)
    )

    for (scale.idx in seq_len(ncol(local.max.by.scale))) {
        maxima.idx <- which(local.max.by.scale[, scale.idx])
        if (length(maxima.idx) < 1L) next

        d.to.max <- igraph::distances(
            graph.obj,
            v = maxima.idx,
            to = seq_len(nrow(local.max.by.scale)),
            weights = igraph::E(graph.obj)$weight
        )
        if (is.null(dim(d.to.max))) {
            d.to.max <- matrix(d.to.max, nrow = 1L)
        }
        support[, scale.idx] <- apply(d.to.max, 2L, min) <= radius
    }

    support
}

.suppress.graph.endpoint.maxima.by.scale <- function(adj.list,
                                                     weight.list,
                                                     local.max.by.scale,
                                                     score.by.scale,
                                                     radius,
                                                     prefer.cpp = TRUE)
{
    local.max.by.scale <- as.matrix(local.max.by.scale)
    local.max.by.scale <- matrix(
        as.logical(local.max.by.scale),
        nrow = nrow(local.max.by.scale),
        ncol = ncol(local.max.by.scale),
        dimnames = dimnames(local.max.by.scale)
    )
    score.by.scale <- as.matrix(score.by.scale)

    if (!all(dim(local.max.by.scale) == dim(score.by.scale))) {
        stop("'local.max.by.scale' and 'score.by.scale' must have the same dimensions.")
    }

    if (isTRUE(prefer.cpp) && exists("rcpp_graph_greedy_maxima_suppression_by_scale", mode = "function")) {
        keep <- rcpp_graph_greedy_maxima_suppression_by_scale(
            adj_list = adj.list,
            weight_list = weight.list,
            local_max_by_scale = local.max.by.scale,
            score_by_scale = score.by.scale,
            radius = as.numeric(radius)
        )
        keep <- as.matrix(keep)
        dimnames(keep) <- dimnames(local.max.by.scale)
        return(keep)
    }

    graph.obj <- .build.graph.endpoint.igraph(adj.list, weight.list)
    keep <- matrix(
        FALSE,
        nrow = nrow(local.max.by.scale),
        ncol = ncol(local.max.by.scale),
        dimnames = dimnames(local.max.by.scale)
    )

    for (scale.idx in seq_len(ncol(local.max.by.scale))) {
        cand <- which(local.max.by.scale[, scale.idx])
        if (length(cand) < 1L) next
        ord <- order(score.by.scale[cand, scale.idx], cand, decreasing = TRUE, na.last = NA)
        cand <- cand[ord]
        suppressed <- rep(FALSE, nrow(local.max.by.scale))
        for (vertex in cand) {
            if (suppressed[vertex]) next
            keep[vertex, scale.idx] <- TRUE
            if (radius <= 0) next
            d <- as.numeric(
                igraph::distances(
                    graph.obj,
                    v = vertex,
                    to = cand,
                    weights = igraph::E(graph.obj)$weight
                )
            )
            suppressed[cand[is.finite(d) & d <= radius]] <- TRUE
        }
    }

    keep
}

.build.graph.endpoint.igraph <- function(adj.list, weight.list) {
    graph.obj <- dgraphs::convert.adjacency.to.edge.matrix(adj.list, weight.list)
    edge.matrix <- graph.obj$edge.matrix
    weights <- graph.obj$weights

    n.vertices <- length(adj.list)
    if (nrow(edge.matrix) < 1L) {
        g <- igraph::make_empty_graph(n = n.vertices, directed = FALSE)
        return(g)
    }

    g <- igraph::graph_from_edgelist(edge.matrix, directed = FALSE)
    if (igraph::vcount(g) < n.vertices) {
        g <- igraph::add_vertices(g, n = n.vertices - igraph::vcount(g))
    }
    g <- igraph::set_edge_attr(g, "weight", value = weights)
    g
}

.replace.nonfinite.endpoint.values <- function(y) {
    y <- as.numeric(y)
    if (all(is.finite(y))) return(y)

    finite.y <- y[is.finite(y)]
    if (length(finite.y) < 1L) {
        return(rep(-1, length(y)))
    }

    range.y <- range(finite.y)
    gap <- max(diff(range.y), 1)
    floor.value <- range.y[[1]] - gap - 1
    y[!is.finite(y)] <- floor.value
    y
}

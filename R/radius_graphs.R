.validate.numeric.data.matrix <- function(X) {
    if (!(is.matrix(X) || is.data.frame(X))) {
        stop("'X' must be a matrix or data frame.", call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X)) {
        stop("'X' must contain numeric data.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("'X' cannot contain NA, NaN, or Inf values.", call. = FALSE)
    }
    if (nrow(X) < 2L) {
        stop("'X' must contain at least two observations.", call. = FALSE)
    }
    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }
    X
}

.pairwise.radius.edges <- function(X, keep.edge) {
    n <- nrow(X)
    rows <- vector("list", n * (n - 1L) / 2L)
    cursor <- 0L
    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            d <- .euclidean.distance(X, i, j)
            if (keep.edge(i, j, d)) {
                cursor <- cursor + 1L
                rows[[cursor]] <- data.frame(from = i, to = j, weight = d)
            }
        }
    }
    if (!cursor) {
        return(data.frame(from = integer(), to = integer(), weight = numeric()))
    }
    do.call(rbind, rows[seq_len(cursor)])
}

.default.radius.prune.k <- function(adj.list) {
    n <- length(adj.list)
    degree <- lengths(adj.list)
    positive.degree <- degree[degree > 0L]
    if (!length(positive.degree)) {
        return(1L)
    }
    as.integer(min(n - 1L, max(1L, round(stats::median(positive.degree)))))
}

.finalize.radius.graph <- function(X,
                                   edges,
                                   connect.components,
                                   connect.method,
                                   bridge.k,
                                   bridge.k.max,
                                   bridge.growth,
                                   class,
                                   prune.method = "none",
                                   max.path.edge.ratio.deviation.thld = 0.1,
                                   path.edge.ratio.percentile = 0.5,
                                   prune.tau = 1.05,
                                   prune.local.k = NULL,
                                   prune.k = 1L,
                                   with.pruned.edge.stats = FALSE) {
    n <- nrow(X)
    graph <- .graph.from.edge.table(n, edges)
    raw.adj.list <- graph$adj_list
    raw.weight.list <- graph$weight_list
    prune.method <- .normalize.local.prune.method(prune.method)
    prune.controls <- .normalize.local.prune.controls(
        n, prune.k, prune.tau, prune.local.k, with.pruned.edge.stats
    )
    global.ratio.controls <- .normalize.global.ratio.prune.controls(
        max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile
    )
    pruning <- .prune.graph.by.method(
        X = X,
        adj.list = raw.adj.list,
        weight.list = raw.weight.list,
        k = prune.k,
        prune.method = prune.method,
        max.path.edge.ratio.deviation.thld =
            global.ratio.controls$max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = global.ratio.controls$path.edge.ratio.percentile,
        prune.tau = prune.controls$prune.tau,
        prune.local.k = prune.controls$prune.local.k,
        with.pruned.edge.stats = prune.controls$with.pruned.edge.stats
    )
    pruned.adj.list <- pruning$adj_list
    pruned.weight.list <- pruning$weight_list
    n.edges.before.mst <- pruning$n_edges_after_pruning
    bridge <- .augment.graph.with.component.mst(
        X = X,
        adj.list = pruned.adj.list,
        weight.list = pruned.weight.list,
        k = prune.k,
        connect.components = connect.components,
        connect.method = connect.method,
        bridge.k = bridge.k,
        bridge.k.max = bridge.k.max,
        bridge.growth = bridge.growth
    )
    edge.table <- .graph.edge.table(bridge$adj_list, bridge$weight_list)
    out <- list(
        adj_list = bridge$adj_list,
        weight_list = bridge$weight_list,
        edge_matrix = as.matrix(edge.table[, c("from", "to"), drop = FALSE]),
        edge_weight = as.numeric(edge.table$weight),
        n_vertices = n,
        n_edges = nrow(edge.table),
        raw_adj_list = raw.adj.list,
        raw_weight_list = raw.weight.list,
        pruned_adj_list = pruned.adj.list,
        pruned_weight_list = pruned.weight.list,
        n_edges_before_mst = n.edges.before.mst,
        n_edges_after_mst = nrow(edge.table),
        n_components_before = bridge$n_components_before,
        n_components_after = bridge$n_components_after,
        component_id_before = bridge$component_id_before,
        component_id_after = bridge$component_id_after,
        mst_edge_matrix = bridge$mst_edge_matrix,
        mst_edge_weight = bridge$mst_edge_weight,
        n_mst_edges_added = bridge$n_mst_edges_added,
        connect_components = bridge$connect_components,
        connect_method = bridge$connect_method,
        bridge_method = bridge$bridge_method,
        bridge_k = bridge$bridge_k,
        bridge_k_max = bridge$bridge_k_max,
        bridge_growth = bridge$bridge_growth,
        bridge_k_used = bridge$bridge_k_used,
        bridge_exact_fallback_used = bridge$bridge_exact_fallback_used,
        n_edges_before_pruning = pruning$n_edges_before_pruning,
        n_edges_after_pruning = pruning$n_edges_after_pruning,
        n_pruned_edges = pruning$n_pruned_edges,
        pruned_edge_stats = pruning$pruned_edge_stats,
        prune_method = prune.method,
        prune_tau = pruning$prune_tau,
        prune_local_k = pruning$prune_local_k,
        with_pruned_edge_stats = pruning$with_pruned_edge_stats
    )
    out <- .add.graph.lifecycle.branches(
        result = out,
        X = X,
        k = prune.k,
        raw.adj.list = out$raw_adj_list,
        raw.weight.list = out$raw_weight_list,
        pruned.adj.list = out$pruned_adj_list,
        pruned.weight.list = out$pruned_weight_list,
        connect.method = connect.method,
        bridge.k = bridge$bridge_k,
        bridge.k.max = bridge$bridge_k_max,
        bridge.growth = bridge$bridge_growth,
        prune.method = prune.method,
        max.path.edge.ratio.deviation.thld =
            global.ratio.controls$max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = global.ratio.controls$path.edge.ratio.percentile,
        prune.tau = prune.controls$prune.tau,
        prune.local.k = prune.controls$prune.local.k,
        with.pruned.edge.stats = prune.controls$with.pruned.edge.stats
    )
    class(out) <- c(class, "list")
    out
}

#' Compute a Fixed-Radius Graph
#'
#' @description
#' Creates an undirected Euclidean radius graph from a numeric data matrix.
#' Vertices `i` and `j` are adjacent exactly when their Euclidean distance is
#' at most `radius`.
#'
#' @param X Numeric matrix or data frame with observations in rows.
#' @param radius Positive numeric scalar. Edges are included when
#'   \eqn{\|x_i - x_j\|_2 \le radius}.
#' @param connect.components Logical scalar. If `TRUE`, add MST bridge edges so
#'   the final graph is connected whenever possible.
#' @param connect.method Character scalar. `"component.mst"` adds exact shortest
#'   inter-component bridges. `"component.mst.ann"` tries sparse ANN bridge
#'   candidates before automatic exact fallback. `"global.mst"` unions the graph
#'   with the full Euclidean MST.
#' @param bridge.k Integer scalar or `NULL`. Initial ANN bridge neighborhood
#'   size for `connect.method = "component.mst.ann"`.
#' @param bridge.k.max Integer scalar or `NULL`. Maximum ANN bridge neighborhood
#'   size before exact fallback.
#' @param bridge.growth Numeric scalar greater than 1. Multiplicative growth
#'   factor for ANN bridge neighborhoods.
#' @param prune.method Character scalar. `"none"` disables geometric pruning.
#'   `"local.geodesic"` applies the experimental local geometric pruning stage
#'   before optional MST connectivity repair. `"global.geodesic.ratio"` applies
#'   whole-graph geodesic-ratio pruning before optional MST connectivity repair.
#' @param max.path.edge.ratio.deviation.thld Numeric scalar in `[0, 0.2)`.
#'   For `prune.method = "global.geodesic.ratio"`, an edge may be removed when
#'   the shortest alternative path is at most
#'   `1 + max.path.edge.ratio.deviation.thld` times the direct edge length.
#' @param path.edge.ratio.percentile Numeric scalar in `[0, 1]`. For
#'   `prune.method = "global.geodesic.ratio"`, only edges at or above this edge
#'   length percentile are considered.
#' @param prune.tau Numeric scalar greater than 1. For local geometric pruning,
#'   an edge may be removed when a retained local alternative path is at most
#'   this multiplicative factor times the direct edge length.
#' @param prune.local.k Integer scalar or `NULL`. Number of nearest neighbors
#'   used to form local neighborhoods for pruning. For fixed-radius graphs,
#'   `NULL` uses the median positive raw graph degree, clipped to `[1, n - 1]`.
#' @param with.pruned.edge.stats Logical scalar. If `TRUE`, return a data frame
#'   with one row per locally pruned edge.
#'
#' @return A list of class `"radius_graph"` containing adjacency lists, edge
#'   weights, edge matrix, and component diagnostics. The final graph is stored
#'   in `adj_list`/`weight_list`; `raw_adj_list`/`raw_weight_list` store the
#'   fixed-radius graph before pruning and optional MST component repair;
#'   `pruned_adj_list`/`pruned_weight_list` store the graph after optional local
#'   geometric pruning and before MST component repair. The repaired lifecycle
#'   branches `raw_repaired_*`, `pruned_repaired_*`, and
#'   `repaired_pruned_*` allow direct comparison of native, prune-first, and
#'   repair-first graph geodesic geometries.
#'
#' @examples
#' X <- matrix(c(0, 1, 3), ncol = 1)
#' create.radius.graph(X, radius = 1.1)$edge_matrix
#'
#' @export
create.radius.graph <- function(X,
                                radius,
                                prune.method = c("none", "local.geodesic", "global.geodesic.ratio"),
                                max.path.edge.ratio.deviation.thld = 0.1,
                                path.edge.ratio.percentile = 0.5,
                                prune.tau = 1.05,
                                prune.local.k = NULL,
                                with.pruned.edge.stats = FALSE,
                                connect.components = FALSE,
                                connect.method = c("component.mst", "component.mst.ann", "global.mst"),
                                bridge.k = NULL,
                                bridge.k.max = NULL,
                                bridge.growth = 2) {
    X <- .validate.numeric.data.matrix(X)
    if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) ||
        radius <= 0) {
        stop("'radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (!is.logical(connect.components) || length(connect.components) != 1L ||
        is.na(connect.components)) {
        stop("'connect.components' must be TRUE or FALSE.", call. = FALSE)
    }
    connect.method <- match.arg(connect.method)
    edges <- .pairwise.radius.edges(
        X,
        keep.edge = function(i, j, d) d <= radius
    )
    graph <- .graph.from.edge.table(nrow(X), edges)
    prune.k <- .default.radius.prune.k(graph$adj_list)
    out <- .finalize.radius.graph(
        X, edges, connect.components, connect.method,
        bridge.k, bridge.k.max, bridge.growth, "radius_graph",
        prune.method = prune.method,
        max.path.edge.ratio.deviation.thld = max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = path.edge.ratio.percentile,
        prune.tau = prune.tau,
        prune.local.k = prune.local.k,
        prune.k = prune.k,
        with.pruned.edge.stats = with.pruned.edge.stats
    )
    out$radius <- as.numeric(radius)
    out$graph_rule <- "fixed.radius"
    out
}

#' Compute an Adaptive-Radius Graph
#'
#' @description
#' Creates an undirected adaptive Euclidean radius graph. Let \eqn{\sigma_i}
#' be the distance from observation \eqn{i} to its `k.scale`-th non-self nearest
#' neighbor. Vertices are adjacent when their Euclidean distance is at most
#' `radius.factor` times either `max(sigma_i, sigma_j)` or
#' `min(sigma_i, sigma_j)`.
#'
#' @param X Numeric matrix or data frame with observations in rows.
#' @param k.scale Positive integer smaller than `nrow(X)`. Defines local scale
#'   distances \eqn{\sigma_i}.
#' @param radius.factor Positive numeric scalar multiplying the adaptive radius.
#' @param radius.rule Character scalar. `"max"` uses
#'   \eqn{d_{ij} \le radius.factor \max(\sigma_i,\sigma_j)}; `"min"` uses
#'   \eqn{d_{ij} \le radius.factor \min(\sigma_i,\sigma_j)}.
#' @param prune.local.k Integer scalar or `NULL`. Number of nearest neighbors
#'   used to form local neighborhoods for pruning. For adaptive-radius graphs,
#'   `NULL` defaults to `k.scale`.
#' @inheritParams create.radius.graph
#'
#' @return A list of class `"adaptive_radius_graph"` containing adjacency
#'   lists, edge weights, edge matrix, local scale values, and component
#'   diagnostics. The final graph is stored in `adj_list`/`weight_list`;
#'   `raw_adj_list`/`raw_weight_list` store the adaptive-radius graph before
#'   pruning and optional MST component repair; and
#'   `pruned_adj_list`/`pruned_weight_list` store the graph after optional local
#'   geometric pruning and before MST component repair. The repaired lifecycle
#'   branches `raw_repaired_*`, `pruned_repaired_*`, and
#'   `repaired_pruned_*` allow direct comparison of native, prune-first, and
#'   repair-first graph geodesic geometries.
#'
#' @examples
#' X <- matrix(c(0, 1, 3), ncol = 1)
#' create.adaptive.radius.graph(X, k.scale = 1, radius.rule = "max")$edge_matrix
#'
#' @export
create.adaptive.radius.graph <- function(X,
                                         k.scale,
                                         radius.factor = 1,
                                         radius.rule = c("max", "min"),
                                         prune.method = c("none", "local.geodesic", "global.geodesic.ratio"),
                                         max.path.edge.ratio.deviation.thld = 0.1,
                                         path.edge.ratio.percentile = 0.5,
                                         prune.tau = 1.05,
                                         prune.local.k = NULL,
                                         with.pruned.edge.stats = FALSE,
                                         connect.components = FALSE,
                                         connect.method = c("component.mst", "component.mst.ann", "global.mst"),
                                         bridge.k = NULL,
                                         bridge.k.max = NULL,
                                         bridge.growth = 2) {
    X <- .validate.numeric.data.matrix(X)
    n <- nrow(X)
    if (!is.numeric(k.scale) || length(k.scale) != 1L || !is.finite(k.scale) ||
        k.scale != floor(k.scale) || k.scale < 1L || k.scale >= n) {
        stop("'k.scale' must be a positive integer smaller than nrow(X).",
             call. = FALSE)
    }
    if (!is.numeric(radius.factor) || length(radius.factor) != 1L ||
        !is.finite(radius.factor) || radius.factor <= 0) {
        stop("'radius.factor' must be a positive finite numeric scalar.",
             call. = FALSE)
    }
    if (!is.logical(connect.components) || length(connect.components) != 1L ||
        is.na(connect.components)) {
        stop("'connect.components' must be TRUE or FALSE.", call. = FALSE)
    }
    radius.rule <- match.arg(radius.rule)
    connect.method <- match.arg(connect.method)

    knn.index <- .exact.knn.index(X, as.integer(k.scale))
    sigma <- numeric(n)
    for (i in seq_len(n)) {
        j <- knn.index[i, as.integer(k.scale)]
        sigma[[i]] <- .euclidean.distance(X, i, j)
    }
    radius.fun <- switch(radius.rule,
                         max = function(a, b) max(a, b),
                         min = function(a, b) min(a, b))
    edges <- .pairwise.radius.edges(
        X,
        keep.edge = function(i, j, d) {
            d <= radius.factor * radius.fun(sigma[[i]], sigma[[j]])
        }
    )
    out <- .finalize.radius.graph(
        X, edges, connect.components, connect.method,
        bridge.k, bridge.k.max, bridge.growth, "adaptive_radius_graph",
        prune.method = prune.method,
        max.path.edge.ratio.deviation.thld = max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = path.edge.ratio.percentile,
        prune.tau = prune.tau,
        prune.local.k = prune.local.k,
        prune.k = as.integer(k.scale),
        with.pruned.edge.stats = with.pruned.edge.stats
    )
    out$k_scale <- as.integer(k.scale)
    out$radius_factor <- as.numeric(radius.factor)
    out$radius_rule <- radius.rule
    out$sigma <- sigma
    out$graph_rule <- "adaptive.radius"
    out
}

#' @export
print.radius_graph <- function(x, ...) {
    cat("Fixed-radius graph\n")
    cat("Number of vertices:", x$n_vertices, "\n")
    cat("Number of edges:", x$n_edges, "\n")
    cat("Radius:", x$radius, "\n")
    cat("Connected components before MST augmentation:", x$n_components_before, "\n")
    cat("Connected components after MST augmentation:", x$n_components_after, "\n")
    invisible(x)
}

#' @export
print.adaptive_radius_graph <- function(x, ...) {
    cat("Adaptive-radius graph\n")
    cat("Number of vertices:", x$n_vertices, "\n")
    cat("Number of edges:", x$n_edges, "\n")
    cat("k.scale:", x$k_scale, "\n")
    cat("Radius factor:", x$radius_factor, "\n")
    cat("Radius rule:", x$radius_rule, "\n")
    cat("Connected components before MST augmentation:", x$n_components_before, "\n")
    cat("Connected components after MST augmentation:", x$n_components_after, "\n")
    invisible(x)
}

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

.finalize.radius.graph <- function(X,
                                   edges,
                                   connect.components,
                                   connect.method,
                                   bridge.k,
                                   bridge.k.max,
                                   bridge.growth,
                                   class) {
    n <- nrow(X)
    graph <- .graph.from.edge.table(n, edges)
    n.edges.before.mst <- nrow(edges)
    bridge <- .augment.graph.with.component.mst(
        X = X,
        adj.list = graph$adj_list,
        weight.list = graph$weight_list,
        k = 1L,
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
        bridge_exact_fallback_used = bridge$bridge_exact_fallback_used
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
#'
#' @return A list of class `"radius_graph"` containing adjacency lists, edge
#'   weights, edge matrix, and component diagnostics.
#'
#' @examples
#' X <- matrix(c(0, 1, 3), ncol = 1)
#' create.radius.graph(X, radius = 1.1)$edge_matrix
#'
#' @export
create.radius.graph <- function(X,
                                radius,
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
    out <- .finalize.radius.graph(
        X, edges, connect.components, connect.method,
        bridge.k, bridge.k.max, bridge.growth, "radius_graph"
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
#' @inheritParams create.radius.graph
#'
#' @return A list of class `"adaptive_radius_graph"` containing adjacency
#'   lists, edge weights, edge matrix, local scale values, and component
#'   diagnostics.
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
        bridge.k, bridge.k.max, bridge.growth, "adaptive_radius_graph"
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

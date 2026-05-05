#' Compute a Symmetric k-Nearest Neighbor Graph
#'
#' @description
#' Creates an undirected symmetric k-nearest neighbor (sKNN) graph from a numeric
#' data matrix. An edge between vertices `i` and `j` is present when either
#' `j` is among the `k` nearest neighbors of `i`, or `i` is among the `k`
#' nearest neighbors of `j`.
#'
#' @details
#' With no nearest-neighbor ties, the sKNN edge rule is equivalent to
#' \deqn{d_{ij} \le \max(\sigma_i, \sigma_j),}
#' where \eqn{\sigma_i} is the distance from point \eqn{i} to its `k`-th
#' non-self nearest neighbor.
#'
#' If `connect.components = TRUE`, the graph is augmented with MST edges. The
#' default `connect.method = "component.mst"` first collapses each connected
#' component of the sKNN graph to one component vertex, weights each component
#' pair by the shortest data-space edge between them, and adds the minimum
#' spanning tree of this component graph. This adds exactly `m - 1` bridge
#' edges when the original sKNN graph has `m` connected components.
#' `connect.method = "component.mst.ann"` uses ANN nearest-neighbor bridge
#' candidates first and automatically falls back to exact component MST when
#' the sparse candidate graph does not connect all components. The alternative
#' `connect.method = "global.mst"` unions the sKNN graph with the full
#' Euclidean MST on the data points.
#'
#' `neighbor.method = "exact"` computes all pairwise distances and is the
#' deterministic reference implementation. `neighbor.method = "ann"` uses
#' the bundled ANN kd-tree for the Euclidean kNN search step; exact parity can
#' differ only in nearest-neighbor tie cases. Edge weights are always returned
#' as ordinary Euclidean distances. Exact `component.mst` repair uses exact
#' pairwise bridge distances; `component.mst.ann` uses ANN bridge candidates
#' before guaranteed exact fallback.
#'
#' @param X Numeric matrix or data frame with observations in rows.
#' @param k Integer scalar. Number of non-self nearest neighbors.
#' @param connect.components Logical scalar. If `TRUE`, add MST bridge edges
#'   so that the final graph is connected whenever possible.
#' @param connect.method Character scalar. `"component.mst"` adds the minimum
#'   number of shortest inter-component bridges. `"component.mst.ann"` tries
#'   sparse ANN bridge candidates before automatic exact fallback. `"global.mst"`
#'   unions the sKNN graph with the Euclidean MST on all vertices.
#' @param edge.weight Character scalar. Currently only `"distance"` is
#'   supported.
#' @param neighbor.method Character scalar. `"exact"` uses all pairwise
#'   distances for kNN construction. `"ann"` uses the bundled ANN kd-tree
#'   backend for the kNN construction itself.
#' @param ann.eps Non-negative numeric scalar reserved for approximate ANN
#'   search. Currently only `0` is supported.
#' @param bridge.k Integer scalar or `NULL`. Initial number of ANN neighbors
#'   to use for sparse inter-component bridge search when
#'   `connect.method = "component.mst.ann"`. Defaults to `max(k, 10)`, capped
#'   at `nrow(X) - 1`.
#' @param bridge.k.max Integer scalar or `NULL`. Maximum number of ANN
#'   neighbors to use before exact fallback. Defaults to
#'   `min(nrow(X) - 1, max(50, 8 * k, bridge.k))`.
#' @param bridge.growth Numeric scalar greater than 1. Multiplicative growth
#'   factor for bridge candidate neighborhoods.
#'
#' @return A list of class `"sknn_graph"` with adjacency lists, weights, edge
#'   matrix, kNN index matrix, component diagnostics, and any MST edges added.
#'
#' @examples
#' X <- rbind(c(0, 0), c(1, 0), c(10, 0), c(11, 0))
#' g <- create.sknn.graph(X, k = 1, connect.components = TRUE)
#' g$n_components_before
#' g$n_components_after
#' g$mst_edge_matrix
#'
#' @export
create.sknn.graph <- function(X,
                              k,
                              connect.components = FALSE,
                              connect.method = c("component.mst", "component.mst.ann", "global.mst"),
                              edge.weight = c("distance"),
                              neighbor.method = c("exact", "ann"),
                              ann.eps = 0,
                              bridge.k = NULL,
                              bridge.k.max = NULL,
                              bridge.growth = 2) {
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
    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    n <- nrow(X)
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
        k != floor(k) || k < 1L || k >= n) {
        stop("'k' must be a positive integer smaller than nrow(X).", call. = FALSE)
    }
    if (!is.logical(connect.components) || length(connect.components) != 1L ||
        is.na(connect.components)) {
        stop("'connect.components' must be TRUE or FALSE.", call. = FALSE)
    }
    connect.method <- match.arg(connect.method)
    edge.weight <- match.arg(edge.weight)
    neighbor.method <- match.arg(neighbor.method)
    if (!is.numeric(ann.eps) || length(ann.eps) != 1L || !is.finite(ann.eps) ||
        ann.eps < 0) {
        stop("'ann.eps' must be a finite non-negative numeric scalar.", call. = FALSE)
    }
    if (identical(neighbor.method, "ann") && !identical(as.numeric(ann.eps), 0)) {
        stop("'ann.eps' values greater than 0 are not yet supported.", call. = FALSE)
    }
    if (!is.numeric(bridge.growth) || length(bridge.growth) != 1L ||
        !is.finite(bridge.growth) || bridge.growth <= 1) {
        stop("'bridge.growth' must be a finite numeric scalar greater than 1.", call. = FALSE)
    }
    if (is.null(bridge.k)) {
        bridge.k <- min(n - 1L, max(as.integer(k), 10L))
    }
    if (!is.numeric(bridge.k) || length(bridge.k) != 1L || !is.finite(bridge.k) ||
        bridge.k != floor(bridge.k) || bridge.k < 1L || bridge.k >= n) {
        stop("'bridge.k' must be a positive integer smaller than nrow(X).", call. = FALSE)
    }
    if (is.null(bridge.k.max)) {
        bridge.k.max <- min(n - 1L, max(50L, 8L * as.integer(k), as.integer(bridge.k)))
    }
    if (!is.numeric(bridge.k.max) || length(bridge.k.max) != 1L ||
        !is.finite(bridge.k.max) || bridge.k.max != floor(bridge.k.max) ||
        bridge.k.max < bridge.k || bridge.k.max >= n) {
        stop("'bridge.k.max' must be an integer between bridge.k and nrow(X) - 1.",
             call. = FALSE)
    }
    bridge.k <- as.integer(bridge.k)
    bridge.k.max <- as.integer(bridge.k.max)

    connect.method.id <- switch(connect.method,
                                component.mst = 0L,
                                component.mst.ann = 2L,
                                global.mst = 1L)
    neighbor.method.id <- switch(neighbor.method,
                                 exact = 0L,
                                 ann = 1L)
    knn.index <- matrix(integer(0), nrow = 0L, ncol = 0L)
    if (identical(neighbor.method, "ann")) {
        raw.knn <- .Call(
            "S_kNN",
            X,
            as.integer(k + 1L),
            PACKAGE = "gflow"
        )$indices
        knn.index <- matrix(NA_integer_, nrow = n, ncol = as.integer(k))
        for (i in seq_len(n)) {
            row <- raw.knn[i, ]
            non.self <- row[row != (i - 1L)]
            if (length(non.self) < k) {
                stop("ANN kNN search returned too few non-self neighbors.", call. = FALSE)
            }
            knn.index[i, ] <- non.self[seq_len(k)] + 1L
        }
        storage.mode(knn.index) <- "integer"
    }
    bridge.knn.index <- matrix(integer(0), nrow = 0L, ncol = 0L)
    if (isTRUE(connect.components) && identical(connect.method, "component.mst.ann")) {
        raw.bridge <- .Call(
            "S_kNN",
            X,
            as.integer(bridge.k.max + 1L),
            PACKAGE = "gflow"
        )$indices
        bridge.knn.index <- matrix(NA_integer_, nrow = n, ncol = bridge.k.max)
        for (i in seq_len(n)) {
            row <- raw.bridge[i, ]
            non.self <- row[row != (i - 1L)]
            if (length(non.self) < bridge.k.max) {
                stop("ANN bridge search returned too few non-self neighbors.", call. = FALSE)
            }
            bridge.knn.index[i, ] <- non.self[seq_len(bridge.k.max)] + 1L
        }
        storage.mode(bridge.knn.index) <- "integer"
    }

    result <- .Call(
        "S_create_sknn_graph",
        X,
        as.integer(k),
        as.logical(connect.components),
        as.integer(connect.method.id),
        as.integer(neighbor.method.id),
        as.numeric(ann.eps),
        knn.index,
        bridge.knn.index,
        as.integer(bridge.k),
        as.integer(bridge.k.max),
        as.numeric(bridge.growth),
        PACKAGE = "gflow"
    )
    result$edge.weight <- edge.weight
    class(result) <- c("sknn_graph", "list")
    result
}

#' @export
print.sknn_graph <- function(x, ...) {
    cat("Symmetric kNN graph\n")
    cat("Number of vertices:", x$n_vertices, "\n")
    cat("Number of edges:", x$n_edges, "\n")
    cat("k:", x$k, "\n")
    cat("Connected components before MST augmentation:", x$n_components_before, "\n")
    cat("Connected components after MST augmentation:", x$n_components_after, "\n")
    if (isTRUE(x$connect_components)) {
        cat("MST bridge edges added:", x$n_mst_edges_added, "\n")
        cat("Connection method:", x$connect_method, "\n")
        cat("Bridge method:", x$bridge_method, "\n")
        if (identical(x$bridge_method, "ann")) {
            cat("Bridge k used:", x$bridge_k_used, "\n")
        }
        if (isTRUE(x$bridge_exact_fallback_used)) {
            cat("Exact bridge fallback used: TRUE\n")
        }
    }
    cat("Neighbor method:", x$neighbor_method, "\n")
    if (identical(x$neighbor_method, "ann")) {
        cat("ANN eps:", x$ann_eps, "\n")
    }
    invisible(x)
}

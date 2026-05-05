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
#' edges when the original sKNN graph has `m` connected components. The
#' alternative `connect.method = "global.mst"` unions the sKNN graph with the
#' full Euclidean MST on the data points.
#'
#' @param X Numeric matrix or data frame with observations in rows.
#' @param k Integer scalar. Number of non-self nearest neighbors.
#' @param connect.components Logical scalar. If `TRUE`, add MST bridge edges
#'   so that the final graph is connected whenever possible.
#' @param connect.method Character scalar. `"component.mst"` adds the minimum
#'   number of shortest inter-component bridges. `"global.mst"` unions the
#'   sKNN graph with the Euclidean MST on all vertices.
#' @param edge.weight Character scalar. Currently only `"distance"` is
#'   supported.
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
                              connect.method = c("component.mst", "global.mst"),
                              edge.weight = c("distance")) {
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

    connect.method.id <- switch(connect.method,
                                component.mst = 0L,
                                global.mst = 1L)

    result <- .Call(
        "S_create_sknn_graph",
        X,
        as.integer(k),
        as.logical(connect.components),
        as.integer(connect.method.id),
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
    }
    invisible(x)
}

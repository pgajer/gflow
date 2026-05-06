.graph.geodesic.fields <- function(graph) {
    if (inherits(graph, "IkNN") ||
        inherits(graph, "sknn_graph") ||
        inherits(graph, "mknn_graph") ||
        inherits(graph, "radius_graph") ||
        inherits(graph, "adaptive_radius_graph") ||
        inherits(graph, "geodesic_iknn_graph")) {
        return(list(
            adj = "adj_list",
            weight = "weight_list",
            convention = paste(
                "Supported gflow graph objects store the final graph in",
                "adj_list/weight_list. raw_* fields hold the native graph and",
                "pruned_* fields hold the post-pruning, pre-MST graph when applicable."
            )
        ))
    }
    stop(
        "'graph' must inherit from one of: IkNN, sknn_graph, mknn_graph, ",
        "radius_graph, adaptive_radius_graph, or geodesic_iknn_graph.",
        call. = FALSE
    )
}

.validate.graph.geodesic.payload <- function(adj.list, weight.list, fields) {
    if (is.null(adj.list) || is.null(weight.list)) {
        stop(
            "Final adjacency and weight lists cannot be NULL. Expected fields: ",
            fields$adj, " and ", fields$weight, ". ", fields$convention,
            call. = FALSE
        )
    }
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("Final adjacency and weight payloads must both be lists.", call. = FALSE)
    }
    if (length(adj.list) != length(weight.list)) {
        stop("Final adjacency and weight lists must have the same size.", call. = FALSE)
    }
    n <- length(adj.list)
    if (!n) {
        stop("Final graph payload must contain at least one vertex.", call. = FALSE)
    }
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        weights <- weight.list[[i]]
        if (!is.numeric(nbrs) || !is.numeric(weights)) {
            stop("Every final adjacency and weight-list entry must be numeric.",
                 call. = FALSE)
        }
        if (length(nbrs) != length(weights)) {
            stop("Final adjacency and weight entries must have matching lengths.",
                 call. = FALSE)
        }
        if (length(nbrs) && (!all(is.finite(nbrs)) ||
            any(nbrs != floor(nbrs)) || any(nbrs < 1L) || any(nbrs > n))) {
            stop("Final adjacency lists must contain valid 1-based vertex indices.",
                 call. = FALSE)
        }
        if (length(weights) && (any(!is.finite(weights)) || any(weights < 0))) {
            stop("Final edge weights must be finite non-negative lengths.",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}

#' Compute Graph Geodesic Distances from a gflow Graph Object
#'
#' @description
#' Extracts the final adjacency and edge-length lists from a supported gflow
#' graph object and computes weighted shortest-path distances using the package's
#' C++-backed [shortest.path()] implementation.
#'
#' @param graph A supported gflow graph object. Currently supported classes are
#'   `"IkNN"`, `"sknn_graph"`, `"mknn_graph"`, `"radius_graph"`,
#'   `"adaptive_radius_graph"`, and `"geodesic_iknn_graph"`.
#' @param vertices Integer vector of 1-based vertex indices. If `NULL`, compute
#'   distances among all vertices in the final graph.
#'
#' @return A numeric matrix of weighted shortest-path distances among
#'   `vertices`.
#'
#' @details
#' The wrapper intentionally derives the final graph payload from the object
#' class and the graph-constructor naming convention. Supported gflow graph
#' objects store the final graph in `adj_list` and `weight_list`. When available,
#' `raw_adj_list`/`raw_weight_list` contain the native graph before pruning, and
#' `pruned_adj_list`/`pruned_weight_list` contain the graph after pruning but
#' before optional MST component repair. This wrapper always uses the final
#' `adj_list`/`weight_list` payload. The extracted payload is checked for
#' non-NULL, equal-size adjacency and weight lists, valid 1-based neighbor
#' indices, and finite non-negative edge lengths before calling
#' [shortest.path()].
#'
#' @examples
#' X <- rbind(c(0, 0), c(1, 0), c(2, 0))
#' g <- create.sknn.graph(X, k = 1)
#' graph.geodesic.distances(g)
#'
#' @export
graph.geodesic.distances <- function(graph, vertices = NULL) {
    fields <- .graph.geodesic.fields(graph)
    adj.list <- graph[[fields$adj]]
    weight.list <- graph[[fields$weight]]
    .validate.graph.geodesic.payload(adj.list, weight.list, fields)

    if (is.null(vertices)) {
        vertices <- seq_along(adj.list)
    }
    if (!is.numeric(vertices) || any(!is.finite(vertices)) ||
        any(vertices != floor(vertices)) || any(vertices < 1L) ||
        any(vertices > length(adj.list))) {
        stop("'vertices' must be valid 1-based vertex indices.", call. = FALSE)
    }
    shortest.path(adj.list, weight.list, as.integer(vertices))
}

.normalize.graph.geodesic.stage <- function(stage) {
    match.arg(
        stage,
        c("final", "raw", "raw.repaired", "pruned", "pruned.repaired",
          "repaired.pruned")
    )
}

.graph.geodesic.fields <- function(graph, stage = "final") {
    if (inherits(graph, "IkNN") ||
        inherits(graph, "sknn_graph") ||
        inherits(graph, "mknn_graph") ||
        inherits(graph, "radius_graph") ||
        inherits(graph, "adaptive_radius_graph") ||
        inherits(graph, "cknn_graph") ||
        inherits(graph, "geodesic_iknn_graph")) {
        stage <- .normalize.graph.geodesic.stage(stage)
        fields <- switch(
            stage,
            final = c("adj_list", "weight_list"),
            raw = c("raw_adj_list", "raw_weight_list"),
            raw.repaired = c("raw_repaired_adj_list", "raw_repaired_weight_list"),
            pruned = c("pruned_adj_list", "pruned_weight_list"),
            pruned.repaired = c("pruned_repaired_adj_list", "pruned_repaired_weight_list"),
            repaired.pruned = c("repaired_pruned_adj_list", "repaired_pruned_weight_list")
        )
        return(list(
            adj = fields[[1L]],
            weight = fields[[2L]],
            stage = stage,
            convention = paste(
                "Supported gflow graph objects store lifecycle graph stages as",
                "raw, raw.repaired, pruned, pruned.repaired, repaired.pruned,",
                "and final."
            )
        ))
    }
    stop(
        "'graph' must inherit from one of: IkNN, sknn_graph, mknn_graph, ",
        "radius_graph, adaptive_radius_graph, cknn_graph, or ",
        "geodesic_iknn_graph.",
        call. = FALSE
    )
}

.validate.graph.geodesic.payload <- function(adj.list, weight.list, fields) {
    if (is.null(adj.list) || is.null(weight.list)) {
        stop(
            "Graph adjacency and weight lists cannot be NULL for stage '",
            fields$stage, "'. Expected fields: ",
            fields$adj, " and ", fields$weight, ". ", fields$convention,
            call. = FALSE
        )
    }
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("Graph adjacency and weight payloads must both be lists.", call. = FALSE)
    }
    if (length(adj.list) != length(weight.list)) {
        stop("Graph adjacency and weight lists must have the same size.", call. = FALSE)
    }
    n <- length(adj.list)
    if (!n) {
        stop("Graph payload must contain at least one vertex.", call. = FALSE)
    }
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        weights <- weight.list[[i]]
        if (!is.numeric(nbrs) || !is.numeric(weights)) {
            stop("Every graph adjacency and weight-list entry must be numeric.",
                 call. = FALSE)
        }
        if (length(nbrs) != length(weights)) {
            stop("Graph adjacency and weight entries must have matching lengths.",
                 call. = FALSE)
        }
        if (length(nbrs) && (!all(is.finite(nbrs)) ||
            any(nbrs != floor(nbrs)) || any(nbrs < 1L) || any(nbrs > n))) {
            stop("Graph adjacency lists must contain valid 1-based vertex indices.",
                 call. = FALSE)
        }
        if (length(weights) && (any(!is.finite(weights)) || any(weights < 0))) {
            stop("Graph edge weights must be finite non-negative lengths.",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}

#' Compute Graph Geodesic Distances from a gflow Graph Object
#'
#' @description
#' Extracts one lifecycle-stage adjacency and edge-length payload from a
#' supported gflow graph object and computes weighted shortest-path distances
#' using the package's C++-backed [shortest.path()] implementation.
#'
#' @param graph A supported gflow graph object. Currently supported classes are
#'   `"IkNN"`, `"sknn_graph"`, `"mknn_graph"`, `"radius_graph"`,
#'   `"adaptive_radius_graph"`, `"cknn_graph"`, and
#'   `"geodesic_iknn_graph"`.
#' @param vertices Integer vector of 1-based vertex indices. If `NULL`, compute
#'   distances among all vertices in the selected graph stage.
#' @param stage Character scalar selecting the graph lifecycle stage. Supported
#'   values are `"final"`, `"raw"`, `"raw.repaired"`, `"pruned"`,
#'   `"pruned.repaired"`, and `"repaired.pruned"`.
#'
#' @return A numeric matrix of weighted shortest-path distances among
#'   `vertices`.
#'
#' @details
#' The wrapper intentionally derives the requested graph payload from the object
#' class and the graph-constructor naming convention. Supported gflow graph
#' objects store `"final"` in `adj_list`/`weight_list`, `"raw"` in
#' `raw_adj_list`/`raw_weight_list`, `"raw.repaired"` in
#' `raw_repaired_adj_list`/`raw_repaired_weight_list`, `"pruned"` in
#' `pruned_adj_list`/`pruned_weight_list`, `"pruned.repaired"` in
#' `pruned_repaired_adj_list`/`pruned_repaired_weight_list`, and
#' `"repaired.pruned"` in
#' `repaired_pruned_adj_list`/`repaired_pruned_weight_list`. The extracted
#' payload is checked for non-NULL, equal-size adjacency and weight lists, valid
#' 1-based neighbor indices, and finite non-negative edge lengths before calling
#' [shortest.path()].
#'
#' @examples
#' X <- rbind(c(0, 0), c(1, 0), c(2, 0))
#' g <- dgraphs::create.sknn.graph(X, k = 1)
#' graph.geodesic.distances(g)
#' graph.geodesic.distances(g, stage = "raw")
#'
#' @export
graph.geodesic.distances <- function(graph, vertices = NULL, stage = "final") {
    fields <- .graph.geodesic.fields(graph, stage = stage)
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

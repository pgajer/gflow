# Internal dispatcher for temporary dgraphs compatibility wrappers.
.gflow.dgraphs.compat <- function(name, args) {
    .Deprecated(
        new = paste0("dgraphs::", name),
        package = "gflow"
    )
    do.call(
        getExportedValue("dgraphs", name),
        args,
        envir = parent.frame()
    )
}

#' Deprecated graph-infrastructure compatibility wrappers
#'
#' These functions contain no graph implementation. They temporarily preserve
#' selected historical `gflow::` entry points while delegating directly to
#' `dgraphs`.
#'
#' @section Lifecycle:
#' Superseded in `gflow 0.2.0` by the corresponding `dgraphs::` API. These
#' wrappers are scheduled for removal in `gflow 0.3.0`.
#'
#' @param gflow.graph A gflow-style basin graph.
#' @param include.vertex.attrs,include.edge.attrs Logical controls forwarded to
#'   [dgraphs::as_igraph()].
#' @return The value returned by the corresponding `dgraphs` function.
#' @name gflow-dgraphs-compat
NULL

#' @rdname gflow-dgraphs-compat
#' @export
as_igraph <- function(gflow.graph,
                      include.vertex.attrs = TRUE,
                      include.edge.attrs = TRUE) {
    .gflow.dgraphs.compat(
        "as_igraph",
        list(
            gflow.graph = gflow.graph,
            include.vertex.attrs = include.vertex.attrs,
            include.edge.attrs = include.edge.attrs
        )
    )
}

#' @rdname gflow-dgraphs-compat
#' @param adj.list An adjacency list.
#' @export
graph.connected.components <- function(adj.list) {
    .gflow.dgraphs.compat(
        "graph.connected.components",
        list(adj.list = adj.list)
    )
}

#' @rdname gflow-dgraphs-compat
#' @param weight.list A list of edge lengths aligned with `adj.list`.
#' @param core.quantile,endpoint.quantile Numeric endpoint-selection quantiles.
#' @param use.approx.eccentricity Logical; use landmark-based eccentricity.
#' @param n.landmarks Number of approximation landmarks.
#' @param max.endpoints Optional maximum number of endpoints.
#' @param seed Random seed used for landmark initialization.
#' @param verbose Logical; report backend progress.
#' @export
geodesic.core.endpoints <- function(adj.list,
                                    weight.list,
                                    core.quantile = 0.10,
                                    endpoint.quantile = 0.90,
                                    use.approx.eccentricity = TRUE,
                                    n.landmarks = 64L,
                                    max.endpoints = NULL,
                                    seed = 1L,
                                    verbose = FALSE) {
    .gflow.dgraphs.compat(
        "geodesic.core.endpoints",
        list(
            adj.list = adj.list,
            weight.list = weight.list,
            core.quantile = core.quantile,
            endpoint.quantile = endpoint.quantile,
            use.approx.eccentricity = use.approx.eccentricity,
            n.landmarks = n.landmarks,
            max.endpoints = max.endpoints,
            seed = seed,
            verbose = verbose
        )
    )
}

#' @rdname gflow-dgraphs-compat
#' @param graph A graph object or adjacency list, depending on the wrapper.
#' @param vertices Optional vertex indices.
#' @param stage Graph lifecycle stage.
#' @export
graph.geodesic.distances <- function(graph, vertices = NULL, stage = "final") {
    .gflow.dgraphs.compat(
        "graph.geodesic.distances",
        list(graph = graph, vertices = vertices, stage = stage)
    )
}

#' @rdname gflow-dgraphs-compat
#' @param edge.lengths Edge lengths aligned with `graph`.
#' @export
shortest.path <- function(graph, edge.lengths, vertices) {
    .gflow.dgraphs.compat(
        "shortest.path",
        list(graph = graph, edge.lengths = edge.lengths, vertices = vertices)
    )
}

#' @rdname gflow-dgraphs-compat
#' @param D.estimated,D.true Estimated and reference distance matrices.
#' @param scale Logical; calibrate estimated distances.
#' @param true.tol Tolerance defining nonzero reference distances.
#' @export
summarize.isometry.deviation <- function(D.estimated,
                                         D.true,
                                         scale = TRUE,
                                         true.tol = sqrt(.Machine$double.eps)) {
    .gflow.dgraphs.compat(
        "summarize.isometry.deviation",
        list(
            D.estimated = D.estimated,
            D.true = D.true,
            scale = scale,
            true.tol = true.tol
        )
    )
}

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

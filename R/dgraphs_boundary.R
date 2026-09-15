.dgraphs.public <- function(name) getExportedValue("dgraphs", name)

## Keep graph algorithms in dgraphs. Adapt only its public representation.
.dgraphs.object.api <- function() {
    "graph" %in% names(formals(dgraphs::graph.connected.components))
}

.dgraphs.components <- function(adj.list) {
    if (.dgraphs.object.api()) {
        dgraphs::graph.connected.components(.dgraphs.public("dgraph")(adj.list))
    } else dgraphs::graph.connected.components(adj.list)
}

.dgraphs.path <- function(adj.list, length.list, h) {
    if (.dgraphs.object.api()) {
        series <- dgraphs::create.path.graph(.dgraphs.public("dgraph")(adj.list, length.list), h.values = h)
        graph <- series[[1L]]$graph
        list(adj.list = .dgraphs.public("graph.adjacency")(graph),
             edge.length.list = .dgraphs.public("graph.lengths")(graph))
    } else dgraphs::create.path.graph(adj.list, length.list, h = h)
}

.dgraphs.embedding <- function(adj.list, weights.list, invert.weights, dim, method, verbose) {
    if (.dgraphs.object.api()) {
        graph <- .dgraphs.public("dgraph")(adj.list, edge.attributes =
            if (is.null(weights.list)) list() else list(association = weights.list))
        dgraphs::graph.embedding(graph,
            edge.attribute = if (is.null(weights.list)) NULL else "association",
            transform = if (invert.weights) "reciprocal" else "identity",
            dim = dim, method = method, verbose = verbose)
    } else dgraphs::graph.embedding(adj.list = adj.list, weights.list = weights.list,
        invert.weights = invert.weights, dim = dim, method = method, verbose = verbose)
}

.dgraphs.nerve <- function(covering) {
    graph <- dgraphs::nerve.graph(covering)
    if (inherits(graph, "dgraph")) {
        list(adjacency.list = .dgraphs.public("graph.adjacency")(graph),
             weights.list = .dgraphs.public("graph.edge.attribute")(graph, "overlap"))
    } else graph
}

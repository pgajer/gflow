
#' Computes Shortest Path Distances for a Selected Set of Vertices of a Graph
#'
#' This function computes the shortest paths between specified vertices in a weighted graph.
#' It serves as an R interface to a C++ implementation of the shortest path algorithm.
#'
#' @param graph A list of numeric vectors. Each vector represents a vertex and contains
#'              the indices of its neighboring vertices.
#' @param edge.lengths A list of numeric vectors. Each vector corresponds to a vertex in
#'                     `graph` and contains the lengths of the edges to its neighbors.
#' @param vertices A numeric vector of vertex indices for which to compute the shortest paths.
#'
#' @return A matrix where the element at position \code{(i,j)} represents the shortest path
#'         distance from \code{vertices[i]} to \code{vertices[j]}.
#'
#' @details The function first performs input validation, ensuring that `graph` and
#'          `edge.lengths` are properly formatted and consistent. It then converts the
#'          graph to 0-based indexing (as required by the C++ function) before calling
#'          the C++ implementation.
#'
#' @note The C++ function assumes 0-based indexing, but this R function accepts 1-based
#'       index inputs as is standard in R. The conversion is handled internally.
#'
#' @examples
#' graph <- list(c(2,3), c(1,3,4), c(1,2,4), c(2,3))
#' edge.lengths <- list(c(1,4), c(1,2,5), c(4,2,1), c(5,1))
#' vertices <- c(1,3,4)
#' result <- shortest.path(graph, edge.lengths, vertices)
#' print(result)
#'
#' @export
shortest.path <- function(graph, edge.lengths, vertices) {
    if (!is.list(graph) || !all(sapply(graph, is.numeric)))
        stop("graph must be a list of numeric vectors.")
    if (!is.list(edge.lengths) || !all(sapply(edge.lengths, is.numeric)))
        stop("edge.lengths must be a list of numeric vectors.")
    if (!all(sapply(seq_along(graph), function(i) length(graph[[i]]) == length(edge.lengths[[i]]))))
        stop("The structure of graph and edge.lengths do not match.")
    if (!is.numeric(vertices))
        stop("vertices must be a numeric vector.")

                                        # Check that vertices are integers within the range 1 to length(graph)
    if (!all(vertices == as.integer(vertices)) ||
        any(vertices < 1) ||
        any(vertices > length(graph)))
        stop("vertices must be integers within the range 1 to length(graph).")

    ## Converting graph to 0-based indexing
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))
    ## Converting vertices to 0-based indexing
    vertices.0based <- as.integer(vertices - 1)

    res <- .Call(S_shortest_path,
                 graph.0based,
                 edge.lengths,
                 vertices.0based)
    return(res)
}

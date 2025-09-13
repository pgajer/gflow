#' Find All Shortest Paths Within a Radius
#'
#' @description
#' Finds all shortest paths from a starting vertex to other vertices in a weighted graph,
#' limited by a maximum distance (radius). This function implements a modified version of
#' Dijkstra's algorithm optimized for local neighborhood exploration.
#'
#' @param adj.list A list of integer vectors, where each vector contains the indices
#'        of vertices adjacent to the vertex at that position in the list.
#'        For example, `adj.list[[1]]` contains the indices of vertices connected to vertex 1.
#'
#' @param weight.list A list of numeric vectors containing the weights of edges in the
#'        corresponding `adj.list`. `weight.list[[i]][j]` is the weight of the edge between
#'        vertex i and vertex `adj.list[[i]][j]`.
#'
#' @param start The index of the starting vertex (1-based indexing as is standard in R).
#'
#' @param radius The maximum allowed distance for paths from the start vertex.
#'        Paths with a total weight greater than this value will not be included in the results.
#'
#' @return A list containing:
#'   \item{paths}{A list of integer vectors, where each vector represents a path from
#'         the start vertex to another vertex in the graph. Each element of the vector
#'         is a vertex index.}
#'   \item{reachable.vertices}{An integer vector containing all vertices that are
#'         reachable from the start vertex within the given radius.}
#'   \item{vertex.to.path.map}{A matrix with three columns: "vertex", "path.index", and
#'         "total.weight". Each row maps a reachable vertex to its corresponding path in
#'         the paths list and includes the total weight (distance) to that vertex. This
#'         provides efficient O(1) lookup for finding the shortest path to any specific vertex.}
#'
#' @details
#' The function converts the R graph representation to a C++ graph, and then uses an efficient
#' implementation of a bounded Dijkstra's algorithm to find all shortest paths within the
#' specified radius. The algorithm works in two phases:
#'
#' 1. A bounded Dijkstra's algorithm to find distances and predecessors
#' 2. Path reconstruction to generate actual paths from the collected information
#'
#' The resulting paths are sorted in descending order of their total weights.
#'
#' @examples
#' # Create a simple graph with 5 vertices
#' adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 5), c(2, 5), c(3, 4))
#' weight.list <- list(c(1, 2), c(1, 1, 3), c(2, 1, 2), c(3, 1), c(2, 1))
#'
#' # Find all shortest paths within radius 3 from vertex 1
#' (res <- find.shortest.paths.within.radius(adj.list, weight.list, 1, 3))
#'
#' @seealso
#' For graph creation and manipulation, consider using the \code{igraph} package.
#'
#' @export
find.shortest.paths.within.radius <- function(adj.list, weight.list, start, radius) {

    if (!is.list(adj.list) || !all(sapply(adj.list, is.numeric))) {
        stop("adj.list must be a list of numeric vectors")
    }

    if (!is.list(weight.list) || !all(sapply(weight.list, is.numeric))) {
        stop("weight.list must be a list of numeric vectors")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(paste("Mismatch between adjacency list and weight list for vertex", i))
        }
    }

    if (!is.numeric(start) || length(start) != 1 || start < 1 || start > length(adj.list)) {
        stop(paste("start must be a single integer between 1 and", length(adj.list)))
    }

    if (!is.numeric(radius) || length(radius) != 1 || radius < 0) {
        stop("radius must be a single non-negative number")
    }

    ## Adjust start index to 0-based for C++
    start.0based <- as.integer(start - 1)

    ## Convert to 0-based indexing for all adjacency lists
    adj.list.0based <- lapply(adj.list, function(neighbors) as.integer(neighbors - 1))

    ## Call the C++ implementation
    result <- .Call(S_find_graph_paths_within_radius,
                    adj.list.0based,
                    weight.list,
                    start.0based,
                    radius)

    return(result)
}

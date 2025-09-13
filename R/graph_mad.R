
#' Calculate Median Absolute Deviation for Graph Vertices
#'
#' @description
#' Computes the Median Absolute Deviation (MAD) for each vertex in a graph based on its
#' neighborhood values. The MAD is calculated using the values of each vertex and its
#' immediate neighbors in the graph.
#'
#' @param graph A list where each element is a numeric vector containing the indices
#'   of neighboring vertices. The indices should be positive integers representing
#'   1-based vertex numbering (as is standard in R).
#' @param y A numeric vector containing the values associated with each vertex.
#'
#' @return A numeric vector of the same length as \code{y}, where each element is
#'   the MAD value for the corresponding vertex. For vertices with no neighbors,
#'   the MAD value will be 0.
#'
#' @details
#' For each vertex, the function:
#' \itemize{
#'   \item Considers the vertex's value and the values of its immediate neighbors
#'   \item Calculates the median of these values
#'   \item Computes the absolute deviations from this median
#'   \item Returns the median of these absolute deviations
#' }
#'
#' @note
#' \itemize{
#'   \item The graph structure should represent an undirected graph
#'   \item Vertex indices in the graph list should be valid (between 1 and length(y))
#'   \item Isolated vertices (those with no neighbors) will have a MAD value of 0
#' }
#'
#' @examples
#' # Create a simple chain graph with 3 vertices
#' g <- list(c(2), c(1,3), c(2))
#' values <- c(1, 10, 2)
#' graph.mad(g, values)
#'
#' @seealso
#' \code{\link{mad}} for the standard MAD calculation on vectors
#'
#' @export
graph.mad <- function(graph,y) {

    if (!is.list(graph)) {
        stop("graph must be a list.")
    }

    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }

    ## if (any(unlist(graph) > n.vertices) || any(unlist(graph) < 1)) {
    ##     stop("Graph contains invalid vertex indices")
    ## }

    n.vertices <- length(y)
    if (length(graph) != n.vertices) {
        stop("The sizes of graph and y are not consistent.")
    }

    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    .Call(S_graph_mad,
          graph.0based,
          as.numeric(y))
}

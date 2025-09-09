
#' Create an Intersection Graph
#'
#' Constructs an intersection graph from an input graph based on neighborhood intersections.
#' In the output graph, an edge exists between two nodes if their neighborhoods in the input
#' graph have a significant intersection.
#'
#' @param adj.list A graph adjacency list of integer vectors or a matrix of
#'     integers. Each vector/row represents the neighbors of a node in the input
#'     graph. If using vectors, they should be sorted in ascending order.
#' @param p.thld A numeric value between 0 and 1. Determines the threshold for
#'     considering an intersection significant. For nodes i and j, the threshold
#'     is calculated as: \code{thld = p.thld * min(length(adj.list[[i]]),
#'     length(adj.list[[j]]))}
#' @param n.itrs Integer. Number of iterations to perform (default is 1).
#'
#' @return A list of integer vectors representing the intersection graph. Each vector
#'         contains the indices of neighbors for the corresponding node in the new graph.
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input parameters.
#' 2. Converts input to 0-based indexing for C++ processing.
#' 3. Calls the C++ implementation to create the intersection graph.
#' 4. Returns the result as a list of integer vectors.
#'
#' For multiple iterations (n.itrs > 1), each iteration uses the result of the previous
#' iteration as its input graph.
#'
#' @examples
#' adj.list <- list(c(1,2), c(0,2,3), c(0,1), c(1))
#' result <- create.intersection.graph(adj.list, 0.5)
#' print(result)
#'
#' @export
create.intersection.graph <- function(adj.list, p.thld, n.itrs = 1) {
    # Convert matrix to list if needed
    if (is.matrix(adj.list)) {
        adj.list <- lapply(seq_len(nrow(adj.list)), function(i) as.integer(adj.list[i,]))
    } else if (!is.list(adj.list)) {
        stop("adj.list must be a list or matrix")
    } else if (!all(sapply(adj.list, is.numeric))) {
        stop("all elements of adj.list must be numeric")
    }

    ## Check for isolated vertices
    if (any(sapply(adj.list, function(x) length(x) == 0))) {
        warning("adj.list contains isolated vertices")
    }

    ## Check for NA values
    if (any(sapply(adj.list, function(x) any(is.na(x))))) {
        stop("adj.list contains NA values")
    }

    ## Check for integer values
    if (!all(sapply(adj.list, function(x) all(x == as.integer(x))))) {
        stop("all values in adj.list must be integers")
    }

    # Validate p.thld
    if (!is.numeric(p.thld) || length(p.thld) != 1 ||
        p.thld < 0 || p.thld > 1) {
        stop("p.thld must be a single numeric value between 0 and 1")
    }

    # Validate n.itrs
    if (!is.numeric(n.itrs) || length(n.itrs) != 1 ||
        n.itrs != as.integer(n.itrs) || n.itrs < 1) {
        stop("n.itrs must be a positive integer")
    }

    # Continue with processing...
}

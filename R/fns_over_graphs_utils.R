
#' Identifies vertices of a graph at which a numeric function is locally constant
#'
#' This function takes an adjacency list representing a graph, a numeric vector of function
#' values at each vertex, and a precision threshold. It identifies the vertices at which the
#' function is considered locally constant, meaning the absolute difference between the function
#' value at the vertex and its neighbors is within the specified precision threshold.
#'
#' @param adj.list A list representing the graph adjacency list. Each element of the list
#'                 should be an integer vector specifying the neighboring vertex indices
#'                 for a given vertex. The vertex indices should be 1-based.
#' @param y A numeric vector representing the function values at each vertex of the graph.
#'          The length of `y` should be equal to the number of vertices in the graph.
#' @param prec A numeric scalar specifying the precision threshold used to determine local
#'             constancy. If the absolute difference between the value at a vertex and any
#'             of its neighbors is greater than this threshold, the vertex is not considered
#'             locally constant.
#'
#' @return An integer vector containing the indices of the vertices at which the function
#'         is locally constant. The vertex indices are 1-based.
#'
#' @export
loc.const.vertices <- function(adj.list, y, prec = 1e-8) {
    if (!is.numeric(y)) {
        stop("y has to be a numeric vector")
    }

    if (length(adj.list) != length(y)) {
        cat("length(adj.list):", length(adj.list), "\n")
        cat("length(y):", length(y), "\n")
        stop("The length of y has to be the same as the number of vertices in the graph.")
    }

    if (!is.numeric(prec)) {
        stop("prec has to be a number")
    }

    if (prec <= 0) {
        stop("prec has to be a positive number")
    }

    ## Converting each component of adj.list to an integer vector
    adj.list <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call(S_loc_const_vertices, adj.list, y, prec)

    result + 1  # Converting back to 1-based indexing
}

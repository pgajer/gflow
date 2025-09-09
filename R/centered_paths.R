#' Compute Graph Path Data with Kernel Weights
#'
#' @description
#' Analyzes paths in a graph centered around a reference vertex, computing distances,
#' kernel weights, and associated values along these paths. The function identifies
#' both single paths and composite paths that meet the minimum size requirement,
#' with the reference vertex serving as a central point in the path structure.
#'
#' @param adj.list A list of integer vectors representing the graph's adjacency structure.
#'        Each element `i` contains the vertices adjacent to vertex `i` (1-based indexing).
#' @param weight.list A list of numeric vectors containing edge weights.
#'        Each element `i` contains weights corresponding to edges in `adj.list[[i]]`.
#' @param y A numeric vector of values associated with each vertex in the graph.
#' @param ref.vertex An integer specifying the reference vertex around which paths are constructed (1-based indexing).
#' @param bandwidth A positive numeric value specifying the maximum allowable path distance from the reference vertex.
#' @param dist.normalization.factor A numeric value between 0 and 1 for normalizing distances in kernel calculations (default: 1.01).
#' @param min.path.size An integer specifying the minimum number of vertices required in a valid path (default: 5).
#' @param diff.threshold An integer specifying the number of vertices after the ref vertex that two paths have to have different (set intersection is empty) to produce a composite path from these two paths. Default is 5. If out of valid range, will be auto-adjusted to the midpoint of the valid range.
#' @param kernel.type An integer specifying the kernel function type:
#'        - 1: Epanechnikov
#'        - 2: Triangular
#'        - 4: Laplace
#'        - 5: Normal
#'        - 6: Biweight
#'        - 7: Tricube (default)
#' @param verbose Logical indicating whether to print progress information. Default is FALSE.
#'
#' @return A list where each element represents a path and contains:
#'   \item{vertices}{Integer vector of path vertices (1-based indices)}
#'   \item{ref_vertex}{Integer indicating the reference vertex (1-based index)}
#'   \item{rel_center_offset}{Numeric value indicating relative position of reference vertex (0 = center, 0.5 = endpoint)}
#'   \item{total_weight}{Numeric value representing total path length}
#'   \item{x_path}{Numeric vector of cumulative distances along path from start}
#'   \item{w_path}{Numeric vector of kernel weights for each vertex}
#'   \item{y_path}{Numeric vector of y-values for path vertices}
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 5 vertices
#' adj.list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2), c(3))
#' weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1), c(1))
#' y <- c(1.5, 2.0, 0.5, 1.0, 1.5)
#'
#' # Find paths centered around vertex 2
#' paths <- get.path.data(adj.list, weight.list, y,
#'                       ref.vertex = 2, bandwidth = 2)
#' }
#' @note
#' - All vertex indices must be positive integers
#' - Edge weights must be non-negative
#' - The length of adj.list and weight.list must match
#' - The reference vertex must exist in the graph
#'
#' @seealso
#' The C++ implementation details can be found in the source code of S_get_path_data
#'
#' @export
get.path.data <- function(adj.list,
                         weight.list,
                         y,
                         ref.vertex,
                         bandwidth,
                         dist.normalization.factor = 1.01,
                         min.path.size = 5L,
                         diff.threshold = 5L,
                         kernel.type = 7L,
                         verbose = FALSE) {
    ## Input validation with specific error messages
    if (!is.list(adj.list)) {
        stop("adj.list must be a list of integer vectors")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list of numeric vectors")
    }
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    # Validate consistency of adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Mismatch in lengths for vertex %d: adj.list has %d entries but weight.list has %d",
                        i, length(adj.list[[i]]), length(weight.list[[i]])))
        }
    }

    # Validate y vector
    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop(sprintf("y must be a numeric vector of length %d (number of vertices)", length(adj.list)))
    }

    # Validate ref.vertex
    if (!is.numeric(ref.vertex) || length(ref.vertex) != 1) {
        stop("ref.vertex must be a single numeric value")
    }
    if (ref.vertex <= 0 || ref.vertex > length(adj.list)) {
        stop(sprintf("ref.vertex must be between 1 and %d", length(adj.list)))
    }

    # Validate bandwidth
    if (!is.numeric(bandwidth) || length(bandwidth) != 1) {
        stop("bandwidth must be a single numeric value")
    }
    if (bandwidth <= 0) {
        stop("bandwidth must be positive")
    }

    # Validate dist.normalization.factor
    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1) {
        stop("dist.normalization.factor must be a single numeric value")
    }
    if (dist.normalization.factor <= 1 || dist.normalization.factor >= 1.5) {
        stop("dist.normalization.factor must be between 1 and 1.5")
    }

    # min.path.size
    if (!is.numeric(min.path.size) || length(min.path.size) != 1) {
        stop("min.path.size must be a single numeric value")
    }
    if (min.path.size <= 2 || min.path.size > length(adj.list)) {
        stop(sprintf("min.path.size must be between 2 and %d", length(adj.list)))
    }

    ## diff.threshold - FIXED: Use floor() to ensure integer value for %d
    if (!is.numeric(diff.threshold) || length(diff.threshold) != 1) {
        stop("diff.threshold must be a single numeric value")
    }
    if (diff.threshold < 2 || diff.threshold > length(adj.list) / 2) {
        old_value <- diff.threshold
        diff.threshold <- as.integer(floor(mean(c(2, floor(length(adj.list) / 2)))))
        if (verbose) {
            message(sprintf("diff.threshold adjusted from %g to %d (valid range: 2 to %d)",
                          old_value, diff.threshold, floor(length(adj.list) / 2)))
        }
    }

    # Validate kernel.type
    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("kernel.type must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    ## Check verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ function
    result <- .Call("S_get_path_data",
                   adj.list.0based,
                   weight.list,
                   as.double(y),
                   as.integer(ref.vertex - 1),
                   as.double(bandwidth),
                   as.double(dist.normalization.factor),
                   as.integer(min.path.size),
                   as.integer(diff.threshold),
                   as.integer(kernel.type),
                   as.logical(verbose))
    return(result)
}

#' Create a Refined Graph with Uniformly Spaced Grid Vertices and creates local paths throught the specified ref vertex
#'
#' @description
#' Creates a refined version of an input graph by adding grid vertices (points) along
#' its edges. The grid vertices are placed to maintain approximately uniform spacing
#' throughout the graph structure. This function is particularly useful for tasks
#' that require a denser sampling of points along the graph edges, such as
#' graph-based interpolation or spatial analysis.
#'
#' Analyzes paths in a graph centered around a reference vertex, computing distances,
#' kernel weights, and associated values along these paths. The function identifies
#' both single paths and composite paths that meet the minimum size requirement,
#' with the reference vertex serving as a central point in the path structure.
#'
#' @param adj.list A list where each element \code{i} is an integer vector containing the
#'        indices of vertices adjacent to vertex \code{i}. Vertex indices must be 1-based.
#'        The graph structure must be undirected, meaning if vertex \code{j} appears in
#'        \code{adj.list[[i]]}, then vertex \code{i} must appear in \code{adj.list[[j]]}.
#'
#' @param weight.list A list matching the structure of adj.list, where each element
#'        contains the corresponding edge weights (typically distances or lengths).
#'        \code{weight.list[[i]][j]} should contain the weight of the edge between vertex \code{i}
#'        and vertex \code{adj.list[[i]][j]}. Weights must be positive numbers.
#'
#' @param grid.size A positive integer specifying the desired number of grid vertices
#'        to add. The actual number of added vertices may differ slightly from this
#'        target due to the distribution of edge lengths in the graph.
#'
#' @param y A numeric vector of values associated with each vertex in the graph.
#' @param ref.vertex An integer specifying the reference vertex around which paths are constructed (1-based indexing).
#' @param bandwidth A positive numeric value specifying the maximum allowable path distance from the reference vertex.
#' @param dist.normalization.factor A numeric value between 0 and 1 for normalizing distances in kernel calculations (default: 1.01).
#' @param min.path.size An integer specifying the minimum number of vertices required in a valid path (default: 5).
#' @param diff.threshold An integer specifying the number of vertices after the ref vertex that two paths have to have different (set intersection is empty) to produce a composite path from these two paths
#' @param kernel.type An integer specifying the kernel function type:
#'        - 1: Epanechnikov
#'        - 2: Triangular
#'        - 4: Laplace
#'        - 5: Normal
#'        - 6: Biweight
#'        - 7: Tricube (default)
#' @param verbose Logical indicating whether to print progress information. Default is FALSE.
#'
#' @return A list where each element represents a path and contains:
#'   \item{vertices}{Integer vector of path vertices (1-based indices)}
#'   \item{ref_vertex}{Integer indicating the reference vertex (1-based index)}
#'   \item{rel_center_offset}{Numeric value indicating relative position of reference vertex (0 = center, 0.5 = endpoint)}
#'   \item{total_weight}{Numeric value representing total path length}
#'   \item{x_path}{Numeric vector of cumulative distances along path from start}
#'   \item{w_path}{Numeric vector of kernel weights for each vertex}
#'   \item{y_path}{Numeric vector of y-values for path vertices}
#'
#' @export
ugg.get.path.data <- function(adj.list,
                              weight.list,
                              grid.size,
                              y,
                              ref.vertex,
                              bandwidth,
                              dist.normalization.factor = 1.01,
                              min.path.size = 5L,
                              diff.threshold = 5L,
                              kernel.type = 7L,
                              verbose = FALSE) {

    ## Input validation with specific error messages
    if (!is.list(adj.list)) {
        stop("adj.list must be a list of integer vectors")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list of numeric vectors")
    }
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    # Validate consistency of adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Mismatch in lengths for vertex %d: adj.list has %d entries but weight.list has %d",
                        i, length(adj.list[[i]]), length(weight.list[[i]])))
        }
    }

    ## Verify positive grid.size
    if (!is.numeric(grid.size) || length(grid.size) != 1 ||
        grid.size != floor(grid.size) || grid.size < 2) {
        stop("grid.size must be a positive integer not smaller than 2")
    }

    # Validate y vector
    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop(sprintf("y must be a numeric vector of length %d (number of vertices)", length(adj.list)))
    }

    # Validate ref.vertex
    if (!is.numeric(ref.vertex) || length(ref.vertex) != 1) {
        stop("ref.vertex must be a single numeric value")
    }
    if (ref.vertex <= 0 || ref.vertex > length(adj.list)) {
        stop(sprintf("ref.vertex must be between 1 and %d", length(adj.list)))
    }

    # Validate bandwidth
    if (!is.numeric(bandwidth) || length(bandwidth) != 1) {
        stop("bandwidth must be a single numeric value")
    }
    if (bandwidth <= 0) {
        stop("bandwidth must be positive")
    }

    # Validate dist.normalization.factor
    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1) {
        stop("dist.normalization.factor must be a single numeric value")
    }
    if (dist.normalization.factor <= 1 || dist.normalization.factor >= 1.5) {
        stop("dist.normalization.factor must be between 1 and 1.5")
    }

    # min.path.size
    if (!is.numeric(min.path.size) || length(min.path.size) != 1) {
        stop("min.path.size must be a single numeric value")
    }
    if (min.path.size <= 2 || min.path.size > length(adj.list)) {
        stop(sprintf("min.path.size must be between 2 and %d", length(adj.list)))
    }

    ## diff.threshold
    if (!is.numeric(diff.threshold) || length(diff.threshold) != 1) {
        stop("diff.threshold must be a single numeric value")
    }
    if (diff.threshold <= 2 || diff.threshold > length(adj.list) / 2) {
        stop(sprintf("diff.threshold must be between 2 and %d", length(adj.list) / 2))
    }

    # Validate kernel.type
    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("kernel.type must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    ## Check verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ function
    result <- .Call("S_ugg_get_path_data",
                   adj.list.0based,
                   weight.list,
                   as.integer(grid.size),
                   as.double(y),
                   as.integer(ref.vertex - 1),
                   as.double(bandwidth),
                   as.double(dist.normalization.factor),
                   as.integer(min.path.size),
                   as.integer(diff.threshold),
                   as.integer(kernel.type),
                   as.logical(verbose))
    return(result)
}

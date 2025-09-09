#' Compute Relative Deviations for Graph Edge Weights
#'
#' @description
#' This function calculates the relative deviation for each edge in a weighted graph
#' to identify potentially redundant edges. For each edge (i,k), it finds the best
#' alternative path (i,j,k) through a common neighbor j and calculates how well
#' the triangle inequality holds.
#'
#' @details
#' The relative deviation is defined as:
#' \deqn{\frac{|w(i,j) + w(j,k) - w(i,k)|}{w(i,j) + w(j,k)}}
#'
#' This measures the fraction of the indirect path length that differs from
#' the direct path. Values close to 0 indicate potential geometric redundancy,
#' meaning the edge (i,k) follows almost perfectly from the edges (i,j) and (j,k).
#'
#' The function is an R wrapper around the C++ implementation that performs
#' the actual calculations. It handles the conversion between R and C++ data structures
#' and returns the results in a format that's easy to work with in R.
#'
#' @param adj_list A list where \code{adj_list[[i]]} contains the indices of vertices
#'   adjacent to vertex i (using 1-based indexing).
#' @param weight_list A list where \code{weight_list[[i]][j]} contains the weight of
#'   the edge from vertex i to vertex \code{adj_list[[i]][j]}.
#'
#' @return A numeric matrix with four columns:
#' \describe{
#'   \item{\code{source}}{The source node ID of the edge.}
#'   \item{\code{target}}{The target node ID of the edge.}
#'   \item{\code{best_intermediate}}{The intermediate node that gives the minimum deviation.}
#'   \item{\code{rel_deviation}}{The minimum relative deviation found for the edge.}
#' }
#' Each row represents one edge in the graph.
#'
#' @examples
#' # Create a simple weighted graph
#' adj_list <- list(c(2, 3, 4), c(1, 3), c(1, 2, 4), c(1, 3))
#' weight_list <- list(c(1, 2, 3), c(1, 1), c(2, 1, 1), c(3, 1))
#'
#' # Compute relative deviations
#' rel_dev <- compute.edge.weight.rel.deviations(adj_list, weight_list)
#'
#' # Display edges with very small relative deviations
#' redundant <- rel_dev[rel_dev[, "rel_deviation"] < 1e-10, ]
#' if (nrow(redundant) > 0) {
#'   print("Potentially redundant edges:")
#'   print(redundant)
#' }
#'
#' @seealso
#' \code{\link{remove.redundant.edges}} for automatic redundant edge removal
#'
#' @export
compute.edge.weight.rel.deviations <- function(adj_list, weight_list) {
    # Input validation
    if (!is.list(adj_list)) {
        stop("'adj_list' must be a list where each element contains the neighbor indices for a vertex")
    }

    if (!is.list(weight_list)) {
        stop("'weight_list' must be a list where each element contains the edge weights corresponding to the adjacency list")
    }

    if (length(adj_list) != length(weight_list)) {
        stop("'adj_list' and 'weight_list' must have the same length (equal to the number of vertices)")
    }

    # Check for empty graph
    if (length(adj_list) == 0) {
        return(matrix(numeric(0), ncol = 4,
                      dimnames = list(NULL, c("source", "target", "best_intermediate", "rel_deviation"))))
    }

    # Check that each element in the adjacency list is a vector of integers
    adj_numeric <- all(sapply(adj_list, function(x) {
        length(x) == 0 || (is.numeric(x) || is.integer(x))
    }))
    if (!adj_numeric) {
        stop("All elements of 'adj_list' must be numeric vectors containing vertex indices")
    }

    # Check that each element in the weight list is a numeric vector
    weight_numeric <- all(sapply(weight_list, function(x) {
        length(x) == 0 || is.numeric(x)
    }))
    if (!weight_numeric) {
        stop("All elements of 'weight_list' must be numeric vectors containing edge weights")
    }

    # Check that the adjacency list and weight list have the same lengths for each vertex
    lengths_match <- all(mapply(function(a, w) length(a) == length(w),
                                adj_list, weight_list))
    if (!lengths_match) {
        stop("For each vertex, the number of neighbors in 'adj_list' must match the number of weights in 'weight_list'")
    }

    # Check for valid vertex indices
    n_vertices <- length(adj_list)
    valid_indices <- all(unlist(lapply(adj_list, function(neighbors) {
        if (length(neighbors) == 0) return(TRUE)
        all(neighbors >= 1 & neighbors <= n_vertices & neighbors == floor(neighbors))
    })))
    if (!valid_indices) {
        stop(paste("Vertex indices in 'adj_list' must be integers between 1 and", n_vertices))
    }

    # Convert adjacency list from 1-based to 0-based indexing for C++
    adj_list_0 <- lapply(adj_list, function(neighbors) {
        if (length(neighbors) == 0) return(integer(0))
        as.integer(neighbors - 1)
    })

    # Ensure weight list contains doubles for C++ compatibility
    weight_list_dbl <- lapply(weight_list, as.double)

    result <- .Call("C_compute_edge_weight_rel_deviations",
                    adj_list_0,
                    weight_list_dbl)

    # Ensure result has proper column names
    if (is.matrix(result) && ncol(result) == 4) {
        colnames(result) <- c("source", "target", "best_intermediate", "rel_deviation")
    }

    return(result)
}


#' Remove Redundant Edges from a Graph While Preserving Connectivity
#'
#' @description
#' This function removes edges that have been identified as redundant (relative deviation of 0)
#' from a graph represented by adjacency and weight lists. It prevents the removal of edges
#' that would cause any vertex to become isolated (have 0 neighbors).
#'
#' @details
#' Redundant edges are those where an alternative two-hop path through an intermediate vertex
#' provides exactly the same total weight as the direct edge. This function processes a matrix
#' of redundant edges (typically those with relative deviation = 0) and systematically removes
#' them from the graph, ensuring no vertex becomes completely disconnected in the process.
#'
#' The function treats the graph as undirected, removing each edge in both directions. If an
#' edge removal would cause any vertex to become isolated (having no neighbors), the function
#' stops with a detailed error message.
#'
#' @param rd0_mat A numeric matrix with at least 3 columns named "source", "target", and
#'   "best_intermediate", containing edges with relative deviation equal to 0.
#' @param adj_list A list where each element i contains the indices of vertices adjacent
#'   to vertex i.
#' @param weight_list A list where each element i contains the weights of edges from
#'   vertex i to its adjacent vertices, in the same order as the adjacency list.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{adj_list}}{The updated adjacency list after redundant edge removal.}
#'   \item{\code{weight_list}}{The updated weight list corresponding to the new adjacency structure.}
#' }
#'
#' @examples
#' # Create a sample graph
#' adj_list <- list(c(2, 3, 4), c(1, 3), c(1, 2, 4), c(1, 3))
#' weight_list <- list(c(1, 2, 3), c(1, 1), c(2, 1, 1), c(3, 1))
#'
#' # Compute relative deviations
#' rd_mat <- compute.edge.weight.rel.deviations(adj_list, weight_list)
#'
#' # Filter for edges with zero relative deviation (tolerance for numerical precision)
#' rd0_mat <- rd_mat[abs(rd_mat[, "rel_deviation"]) < 1e-10, , drop = FALSE]
#'
#' if (nrow(rd0_mat) > 0) {
#'   # Remove redundant edges while preserving connectivity
#'   new_graph <- remove.redundant.edges.with.rel.dev(rd0_mat, adj_list, weight_list)
#'
#'   # Display the updated graph structure
#'   print("Updated adjacency list:")
#'   print(new_graph$adj_list)
#' }
#'
#' @seealso
#' \code{\link{remove.redundant.edges}} for automatic redundant edge removal
#'
#' @export
remove.redundant.edges.with.rel.dev <- function(rd0_mat, adj_list, weight_list) {
    # Validate inputs
    if (!is.matrix(rd0_mat) && !is.data.frame(rd0_mat)) {
        stop("'rd0_mat' must be a matrix or data frame")
    }

    # Convert to matrix if data frame
    if (is.data.frame(rd0_mat)) {
        rd0_mat <- as.matrix(rd0_mat)
    }

    # Check for required columns
    required_cols <- c("source", "target", "best_intermediate")
    if (!all(required_cols %in% colnames(rd0_mat))) {
        stop("'rd0_mat' must have columns named: ", paste(required_cols, collapse = ", "))
    }

    if (!is.list(adj_list) || !is.list(weight_list)) {
        stop("'adj_list' and 'weight_list' must be lists")
    }

    if (length(adj_list) != length(weight_list)) {
        stop("'adj_list' and 'weight_list' must have the same length")
    }

    # Handle empty edge set
    if (nrow(rd0_mat) == 0) {
        message("No redundant edges to remove.")
        return(list(adj_list = adj_list, weight_list = weight_list))
    }

    # Work with copies of the lists
    adj_list_new <- adj_list
    weight_list_new <- weight_list

    # Track successfully removed edges
    removed_edges <- 0
    skipped_edges <- 0
    isolation_prevented <- 0

    # Process each redundant edge
    for (i in seq_len(nrow(rd0_mat))) {
        source_idx <- rd0_mat[i, "source"]
        target_idx <- rd0_mat[i, "target"]
        intermediate <- rd0_mat[i, "best_intermediate"]

        # Validate indices
        n_vertices <- length(adj_list_new)
        if (source_idx < 1 || source_idx > n_vertices ||
            target_idx < 1 || target_idx > n_vertices) {
            warning(sprintf("Invalid vertex indices in row %d: source=%d, target=%d (skipping)",
                           i, source_idx, target_idx))
            skipped_edges <- skipped_edges + 1
            next
        }

        # Check if removing this edge would isolate either vertex
        if (length(adj_list_new[[source_idx]]) <= 1) {
            isolation_prevented <- isolation_prevented + 1
            next  # Skip this edge to prevent isolation
        }

        if (length(adj_list_new[[target_idx]]) <= 1) {
            isolation_prevented <- isolation_prevented + 1
            next  # Skip this edge to prevent isolation
        }

        # Remove source->target edge
        target_pos <- which(adj_list_new[[source_idx]] == target_idx)
        if (length(target_pos) > 0) {
            adj_list_new[[source_idx]] <- adj_list_new[[source_idx]][-target_pos[1]]
            weight_list_new[[source_idx]] <- weight_list_new[[source_idx]][-target_pos[1]]

            # For undirected graph, also remove target->source edge
            source_pos <- which(adj_list_new[[target_idx]] == source_idx)
            if (length(source_pos) > 0) {
                adj_list_new[[target_idx]] <- adj_list_new[[target_idx]][-source_pos[1]]
                weight_list_new[[target_idx]] <- weight_list_new[[target_idx]][-source_pos[1]]
            }

            removed_edges <- removed_edges + 1
        } else {
            # Edge not found - might have been removed already
            skipped_edges <- skipped_edges + 1
        }
    }

    # Report summary
    message(sprintf("Processed %d redundant edges:", nrow(rd0_mat)))
    message(sprintf("  - Removed: %d edges", removed_edges))
    if (skipped_edges > 0) {
        message(sprintf("  - Skipped: %d edges (not found or already removed)", skipped_edges))
    }
    if (isolation_prevented > 0) {
        message(sprintf("  - Preserved: %d edges to prevent vertex isolation", isolation_prevented))
    }

    # Return the updated lists
    return(list(adj_list = adj_list_new, weight_list = weight_list_new))
}


#' Remove Redundant Edges from a Weighted Graph
#'
#' @description
#' Identifies and removes redundant edges from a weighted graph while ensuring no vertex
#' becomes isolated (has zero neighbors). An edge is considered redundant if its relative
#' deviation is effectively zero, meaning there exists an alternative path through another
#' vertex with exactly the same total weight.
#'
#' @details
#' This function serves as an R interface to a C++ implementation that:
#' \enumerate{
#'   \item Computes the relative deviation for each edge using
#'     \code{\link{compute.edge.weight.rel.deviations}}
#'   \item Identifies edges with relative deviation values extremely close to zero (< 1e-16)
#'   \item Removes these redundant edges while ensuring no vertex becomes isolated
#'   \item Returns the updated graph structure
#' }
#'
#' The function will stop with an error message if removing an edge would leave any vertex
#' with zero neighbors, providing details about which edge was being removed and which vertex
#' would become isolated.
#'
#' @param adj_list A list where each element \code{i} contains the indices of vertices adjacent
#'   to vertex \code{i} (using 1-based indexing).
#' @param weight_list A list where each element \code{i} contains the weights of edges from
#'   vertex \code{i} to its adjacent vertices, in the same order as \code{adj_list[[i]]}.
#' @param tolerance Numeric value specifying the threshold below which a relative deviation
#'   is considered zero. Default is 1e-16.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{adj_list}}{The updated adjacency list after redundant edge removal
#'     (using 1-based indexing).}
#'   \item{\code{weight_list}}{The updated weight list corresponding to the new
#'     adjacency structure.}
#' }
#'
#' @examples
#' # Create a sample graph with a redundant edge
#' # Edge 1-3 has the same total weight as the path 1-2-3
#' adj_list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight_list <- list(c(1, 3), c(1, 2), c(3, 2))
#'
#' # Remove redundant edges
#' simplified_graph <- remove.redundant.edges(adj_list, weight_list)
#'
#' # Display the updated graph structure
#' cat("Original edges:", sum(lengths(adj_list)) / 2, "\n")
#' cat("Remaining edges:", sum(lengths(simplified_graph$adj_list)) / 2, "\n")
#' print(simplified_graph$adj_list)
#'
#' @seealso
#' \code{\link{remove.redundant.edges.with.rel.dev}} for removing specific redundant edges
#'
#' @note
#' The implementation processes edges in the order they are computed by the relative
#' deviation function, not in any particular order of deviation magnitude. For more
#' control over which edges are removed, use \code{\link{remove.redundant.edges.with.rel.dev}}
#' with a custom filtered edge set.
#'
#' @export
remove.redundant.edges <- function(adj_list, weight_list, tolerance = 1e-16) {
    # Input validation
    if (!is.list(adj_list) || !is.list(weight_list)) {
        stop("'adj_list' and 'weight_list' must be lists")
    }

    if (length(adj_list) != length(weight_list)) {
        stop("'adj_list' and 'weight_list' must have the same length")
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance < 0) {
        stop("'tolerance' must be a single non-negative numeric value")
    }

    # Handle empty graph
    if (length(adj_list) == 0) {
        return(list(adj_list = adj_list, weight_list = weight_list))
    }

    # Convert adjacency list from 1-based to 0-based indexing for C++
    adj_list_0 <- lapply(adj_list, function(neighbors) {
        if (length(neighbors) == 0) return(integer(0))
        as.integer(neighbors - 1)
    })

    # Ensure weight list contains doubles for C++ compatibility
    weight_list_dbl <- lapply(weight_list, as.double)

    result <- .Call("C_remove_redundant_edges",
                    adj_list_0,
                    weight_list_dbl,
                    as.double(tolerance))

    return(result)
}

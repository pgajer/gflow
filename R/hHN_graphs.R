#' Create a k-hop Neighborhood (hHN) Graph
#'
#' Generates a k-hop neighborhood graph from an input graph represented by an
#' adjacency list and corresponding edge lengths. The k-hop neighborhood of a
#' vertex includes all vertices reachable within at most k hops, with edges
#' representing the shortest paths between vertices.
#'
#' @param graph A list of numeric vectors representing the adjacency list of
#'   the input graph. Each element \code{graph[[i]]} contains the indices of
#'   vertices adjacent to vertex i. Note: Uses 1-based indexing as standard in R.
#' @param edge.lengths A list of numeric vectors representing the edge weights
#'   of the input graph. \code{edge.lengths[[i]][j]} is the weight of the edge
#'   from vertex i to its j-th neighbor in \code{graph[[i]]}.
#' @param h An integer specifying the maximum number of hops to consider for
#'   neighborhood connectivity. Must be greater than or equal to 1.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{adj_list}{A list of integer vectors representing the adjacency list
#'     of the hHN graph. Indices use 1-based indexing.}
#'   \item{dist_list}{A list of numeric vectors representing the edge weights
#'     (shortest path distances) of the hHN graph.}
#' }
#'
#' The time complexity is O(n * (m + n log n)), where n is the number of
#' vertices and m is the number of edges in the original graph.
#'
#' @examples
#' # Create a simple 4-vertex graph
#' graph <- list(
#'   c(2, 3),      # Vertex 1 connects to vertices 2 and 3
#'   c(1, 3, 4),   # Vertex 2 connects to vertices 1, 3, and 4
#'   c(1, 2, 4),   # Vertex 3 connects to vertices 1, 2, and 4
#'   c(2, 3)       # Vertex 4 connects to vertices 2 and 3
#' )
#'
#' edge.lengths <- list(
#'   c(1, 2),      # Edge weights from vertex 1
#'   c(1, 1, 3),   # Edge weights from vertex 2
#'   c(2, 1, 1),   # Edge weights from vertex 3
#'   c(3, 1)       # Edge weights from vertex 4
#' )
#'
#' # Generate the 2-hop neighborhood graph
#' h <- 2
#' khn.graph <- create.hHN.graph(graph, edge.lengths, h)
#'
#' # Access the results
#' khn.graph$adj_list   # Adjacency list of the hHN graph
#' khn.graph$dist_list  # Shortest path distances in the hHN graph
#'
#' @seealso
#' \code{\link{bbmwd.over.hHN.graphs}} for computing BBMWD over hHN graphs
#'
#' @export
create.hHN.graph <- function(graph, edge.lengths, h) {
    # Input validation
    if (!is.list(graph)) {
        stop("'graph' must be a list", call. = FALSE)
    }

    if (!all(vapply(graph, is.numeric, logical(1)))) {
        stop("All elements of 'graph' must be numeric vectors", call. = FALSE)
    }

    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list", call. = FALSE)
    }

    if (!all(vapply(edge.lengths, is.numeric, logical(1)))) {
        stop("All elements of 'edge.lengths' must be numeric vectors", call. = FALSE)
    }

    if (length(graph) != length(edge.lengths)) {
        stop("'graph' and 'edge.lengths' must have the same length", call. = FALSE)
    }

    # Check that each vertex has matching adjacency and edge length vectors
    length_mismatch <- vapply(
        seq_along(graph),
        function(i) length(graph[[i]]) != length(edge.lengths[[i]]),
        logical(1)
    )

    if (any(length_mismatch)) {
        stop("The structure of 'graph' and 'edge.lengths' do not match. ",
             "Each vertex must have the same number of neighbors and edge weights.",
             call. = FALSE)
    }

    # Validate h parameter
    if (!is.numeric(h) || length(h) != 1) {
        stop("'h' must be a single numeric value", call. = FALSE)
    }

    h <- as.integer(h)
    if (is.na(h) || h < 1) {
        stop("'h' must be greater than or equal to 1", call. = FALSE)
    }

    # Convert graph to 0-based indexing for C++ compatibility
    graph_0based <- lapply(graph, function(x) as.integer(x - 1))

    # Call the C++ implementation
    result <- .Call("S_create_hHN_graph",
                    graph_0based,
                    edge.lengths,
                    h)

    return(result)
}


#' Calculate Bayesian Bootstrap Mean Wasserstein Distance over k-Hop Neighbor Graphs
#'
#' Computes the Bayesian Bootstrap Mean Wasserstein Distance (BBMWD) over
#' k-Hop Neighbor (hHN) graphs for various combinations of nearest neighbors (k)
#' and hop values (h).
#'
#' @param y A numeric vector of binary outcomes (0 or 1) associated with each
#'   data point. Length must match the number of vertices in the graphs.
#' @param IkNN.graphs A list of pre-computed intersection k-Nearest Neighbor
#'   graphs. Each element must be a list containing:
#'   \describe{
#'     \item{pruned_adj_list}{Adjacency list of the pruned graph}
#'     \item{pruned_dist_list}{Distance list of the pruned graph}
#'   }
#' @param k.values An integer vector specifying the k values to consider.
#'   Default is \code{10:30}.
#' @param hop.values An integer vector specifying the hop values to consider.
#'   Default is \code{2:15}.
#' @param n.BB Integer; the number of bootstrap iterations for Bayesian
#'   Bootstrap. Default is 100.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'   Default is \code{TRUE}.
#'
#' @return A nested list structure where:
#' \describe{
#'   \item{Outer list}{Indexed by k values from \code{k.values}}
#'   \item{Inner list}{Contains results for each hop value, with elements:
#'     \describe{
#'       \item{h}{The hop value}
#'       \item{khn.graph}{The k-Hop Neighbor graph (list with adj_list and dist_list)}
#'       \item{bbmwd}{The calculated Bayesian Bootstrap Mean Wasserstein Distance}
#'     }
#'   }
#' }
#'
#' @details
#' This function iterates over all combinations of k and hop values, creating
#' hHN graphs and computing BBMWD for each. Progress is displayed if
#' \code{verbose = TRUE}.
#'
#' @examples
#' \dontrun{
#' # Assuming IkNN.graphs is pre-computed
#' set.seed(123)
#' n <- 100
#' y <- sample(0:1, n, replace = TRUE)
#'
#' # Run BBMWD analysis over a subset of k and hop values
#' results <- bbmwd.over.hHN.graphs(
#'   y = y,
#'   IkNN.graphs = IkNN.graphs,
#'   k.values = 10:15,
#'   hop.values = 2:5,
#'   n.BB = 50
#' )
#'
#' # Access specific result
#' bbmwd_k10_h3 <- results[[10]][[3]]$bbmwd
#' }
#'
#' @seealso
#' \code{\link{create.hHN.graph}},
#'
#' @export
bbmwd.over.hHN.graphs <- function(y,
                                  IkNN.graphs,
                                  k.values = 10:30,
                                  hop.values = 2:15,
                                  n.BB = 100,
                                  verbose = TRUE) {

    # Input validation
    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    if (!all(y %in% c(0, 1))) {
        stop("'y' must contain only binary values (0 or 1)", call. = FALSE)
    }

    if (!is.list(IkNN.graphs)) {
        stop("'IkNN.graphs' must be a list", call. = FALSE)
    }

    if (!is.numeric(k.values) || !all(k.values == round(k.values))) {
        stop("'k.values' must be a vector of integers", call. = FALSE)
    }

    if (!is.numeric(hop.values) || !all(hop.values == round(hop.values))) {
        stop("'hop.values' must be a vector of integers", call. = FALSE)
    }

    if (!is.numeric(n.BB) || length(n.BB) != 1 || n.BB < 1) {
        stop("'n.BB' must be a positive integer", call. = FALSE)
    }

    # Initialize results list
    results <- list()

    # Progress tracking
    total_iterations <- length(k.values) * length(hop.values)
    current_iteration <- 0

    for (k in k.values) {
        # Check if k exists in IkNN.graphs
        if (is.null(IkNN.graphs[[k]])) {
            warning(sprintf("No graph found for k = %d, skipping", k),
                    call. = FALSE)
            next
        }

        IkNN.graph <- IkNN.graphs[[k]]

        # Validate IkNN.graph structure
        if (!all(c("pruned_adj_list", "pruned_dist_list") %in% names(IkNN.graph))) {
            stop(sprintf("IkNN.graphs[[%d]] must contain 'pruned_adj_list' and 'pruned_dist_list'", k),
                 call. = FALSE)
        }

        hop.results <- list()

        for (h in hop.values) {
            current_iteration <- current_iteration + 1

            if (verbose) {
                cat(sprintf("\rProcessing k=%d, h=%d (%d/%d)        ",
                           k, h, current_iteration, total_iterations))
            }

            # Create k-hop neighborhood graph
            khn.graph <- create.hHN.graph(
                graph = IkNN.graph$pruned_adj_list,
                edge.lengths = IkNN.graph$pruned_dist_list,
                h = h
            )

            # Calculate BBMWD
            ## kmean.cri <- graph.kernel.weighted.mean.bb.cri(
            ##     graph = khn.graph$adj_list,
            ##     edge.lengths = khn.graph$dist_list,
            ##     y = y,
            ##     n.BB = n.BB
            ## )

            # Store results
            res <- list(
                h = h,
                khn.graph = khn.graph
                ##bbmwd = kmean.cri$bbmwd
            )

            hop.results[[as.character(h)]] <- res
        }

        results[[as.character(k)]] <- hop.results
    }

    if (verbose) cat("\n")

    return(results)
}


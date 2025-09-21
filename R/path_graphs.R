#' Create a Path Graph with Limited Hop Distance
#'
#' @description
#' Constructs a path graph from an input graph by finding all reachable vertices
#' within a specified maximum number of hops. For each pair of vertices, stores
#' the shortest path information along with edge weights and hop counts.
#'
#' @usage
#' create.path.graph(graph, edge.lengths, h)
#'
#' @param graph A list where each element \code{i} is a numeric vector containing the
#'        indices of vertices adjacent to vertex \code{i}. Indices should be 1-based
#'        (following R convention). Each vertex \code{i} should have an entry in the
#'        list, even if it has no neighbors (empty vector).
#' @param edge.lengths A list with the same structure as \code{graph} where each
#'        element \code{i} contains the weights/lengths of the edges specified in
#'        \code{graph[[i]]}. Must contain numeric values > 0. The order of weights
#'        must correspond to the order of vertices in \code{graph[[i]]}.
#' @param h Integer >= 1 specifying the maximum number of hops allowed in the
#'        path graph. A hop is a single edge traversal.
#'
#' @return An S3 object of class \code{"path.graph"} containing:
#' \describe{
#'   \item{adj.list}{A list where \code{adj.list[[i]]} contains all vertices
#'         reachable from vertex \code{i} within \code{h} hops}
#'   \item{edge.length.list}{A list where \code{edge.length.list[[i]]} contains
#'         the shortest path lengths to each vertex in \code{adj.list[[i]]}}
#'   \item{hop.list}{A list where \code{hop.list[[i]]} contains the number of
#'         hops required to reach each vertex in \code{adj.list[[i]]}}
#'   \item{shortest.paths}{A list with components:
#'         \describe{
#'           \item{i}{Source vertex indices}
#'           \item{j}{Target vertex indices}
#'           \item{paths}{List of vertex sequences for each shortest path}
#'         }}
#' }
#'
#' @details
#' The function uses Dijkstra's algorithm to compute shortest paths from each
#' vertex to all reachable vertices within the specified hop limit. The algorithm
#' is implemented in C++ for performance. The resulting path graph can be used
#' for various network analyses where connectivity within a limited number of
#' steps is relevant.
#'
#' @note
#' \itemize{
#'   \item All vertex indices in input and output are 1-based (R convention)
#'   \item The function internally converts indices to 0-based for C++ processing
#'   \item Empty adjacency lists are allowed for isolated vertices
#'   \item The graph can be directed or undirected
#'   \item Self-loops are not allowed
#' }
#'
#' @seealso
#' \code{\link[igraph]{shortest_paths}} for alternative shortest path calculations,
#' \code{\link{create.path.graph.series}} for creating multiple path graphs efficiently
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 3 vertices
#' graph <- list(
#'   c(2),    # Vertex 1 connected to 2
#'   c(1, 3), # Vertex 2 connected to 1 and 3
#'   c(2)     # Vertex 3 connected to 2
#' )
#'
#' # Define edge lengths
#' edge.lengths <- list(
#'   c(1.0),      # Length of edge 1->2
#'   c(1.0, 2.0), # Lengths of edges 2->1 and 2->3
#'   c(2.0)       # Length of edge 3->2
#' )
#'
#' # Create path graph with maximum 2 hops
#' pg <- create.path.graph(graph, edge.lengths, h = 2)
#'
#' # Print the path graph
#' print(pg)
#'
#' # Get shortest path between vertices 1 and 3
#' path <- get.shortest.path(pg, 1, 3)
#' }
#'
#' @export
create.path.graph <- function(graph, edge.lengths, h) {
    ## Input validation
    if (!is.list(graph)) {
        stop("'graph' must be a list.", call. = FALSE)
    }

    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list.", call. = FALSE)
    }

    if (length(graph) == 0) {
        stop("'graph' must not be empty.", call. = FALSE)
    }

    if (length(graph) != length(edge.lengths)) {
        stop("'graph' and 'edge.lengths' must have the same length.", call. = FALSE)
    }

    # Validate each adjacency list and edge length list
    for (i in seq_along(graph)) {
        if (!is.numeric(graph[[i]]) && length(graph[[i]]) > 0) {
            stop(sprintf("graph[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (!is.numeric(edge.lengths[[i]]) && length(edge.lengths[[i]]) > 0) {
            stop(sprintf("edge.lengths[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop(sprintf("graph[[%d]] and edge.lengths[[%d]] must have the same length.", i, i),
                 call. = FALSE)
        }

        # Check for positive edge lengths
        if (length(edge.lengths[[i]]) > 0 && any(edge.lengths[[i]] <= 0)) {
            stop(sprintf("All edge lengths in edge.lengths[[%d]] must be positive.", i),
                 call. = FALSE)
        }

        # Check for valid vertex indices
        if (length(graph[[i]]) > 0) {
            invalid_idx <- graph[[i]] < 1 | graph[[i]] > length(graph)
            if (any(invalid_idx)) {
                stop(sprintf("Invalid vertex indices in graph[[%d]]: indices must be between 1 and %d.",
                     i, length(graph)), call. = FALSE)
            }

            # Check for self-loops
            if (i %in% graph[[i]]) {
                stop(sprintf("Self-loop detected at vertex %d.", i), call. = FALSE)
            }
        }
    }

    # Validate h
    if (!is.numeric(h) || length(h) != 1) {
        stop("'h' must be a single numeric value.", call. = FALSE)
    }

    h <- as.integer(h)
    if (h < 1) {
        stop("'h' must be at least 1.", call. = FALSE)
    }

    ## Convert graph to 0-based indexing for C++
    graph.0based <- lapply(graph, function(x) {
        if (length(x) == 0) integer(0) else as.integer(x - 1)
    })

    ## Call C++ implementation
    res <- .Call(S_create_path_graph_plus,
                 graph.0based,
                 edge.lengths,
                 h)

    ## Create and return path.graph object
    new.path.graph(
        adj.list = res$adj_list,
        edge.length.list = res$edge_length_list,
        hop.list = res$hop_list,
        shortest.paths = res$shortest_paths
    )
}


#' Constructor for path.graph S3 Class
#'
#' @description
#' Internal constructor function for creating path.graph objects.
#' This function should not be called directly by users.
#'
#' @param adj.list List of adjacency lists
#' @param edge.length.list List of edge lengths
#' @param hop.list List of hop counts
#' @param shortest.paths List containing shortest path information
#'
#' @return An object of class "path.graph"
#'
#' @keywords internal
new.path.graph <- function(adj.list, edge.length.list, hop.list, shortest.paths) {
    structure(
        list(
            adj.list = adj.list,
            edge.length.list = edge.length.list,
            hop.list = hop.list,
            shortest.paths = shortest.paths
        ),
        class = "path.graph"
    )
}


#' Get Shortest Path Between Two Vertices
#'
#' @description
#' Retrieves the shortest path between two vertices in a path graph object.
#'
#' @param pg A path.graph object created by \code{\link{create.path.graph}}
#' @param from Source vertex index (1-based)
#' @param to Target vertex index (1-based)
#'
#' @return A list containing:
#' \describe{
#'   \item{path}{Integer vector of vertex indices representing the path}
#'   \item{length}{Numeric value of the total path length}
#'   \item{hops}{Integer number of hops in the path}
#' }
#' Returns NULL if no path exists between the vertices.
#'
#' @examples
#' \dontrun{
#' # Create a path graph
#' pg <- create.path.graph(graph, edge.lengths, h = 3)
#'
#' # Get shortest path from vertex 1 to vertex 3
#' path_info <- get.shortest.path(pg, 1, 3)
#' if (!is.null(path_info)) {
#'   cat("Path:", path_info$path, "\n")
#'   cat("Length:", path_info$length, "\n")
#' }
#' }
#'
#' @export
get.shortest.path <- function(pg, from, to) {
    if (!inherits(pg, "path.graph")) {
        stop("'pg' must be a path.graph object.", call. = FALSE)
    }

    # Validate vertex indices
    n.vertices <- length(pg$adj.list)
    if (!is.numeric(from) || length(from) != 1 || from < 1 || from > n.vertices) {
        stop(sprintf("'from' must be a single integer between 1 and %d.", n.vertices),
             call. = FALSE)
    }

    if (!is.numeric(to) || length(to) != 1 || to < 1 || to > n.vertices) {
        stop(sprintf("'to' must be a single integer between 1 and %d.", n.vertices),
             call. = FALSE)
    }

    from <- as.integer(from)
    to <- as.integer(to)

    # Find matching path
    idx <- which(pg$shortest.paths$i == from & pg$shortest.paths$j == to)

    if (length(idx) == 0) {
        return(NULL)  # No path found
    }

    # Extract path information
    path <- pg$shortest.paths$paths[[idx]]

    # Calculate path length
    path.length <- 0
    if (length(path) > 1) {
        for (k in 1:(length(path) - 1)) {
            v1 <- path[k]
            v2 <- path[k + 1]
            # Find edge in adjacency list
            edge.idx <- which(pg$adj.list[[v1]] == v2)
            if (length(edge.idx) > 0) {
                path.length <- path.length + pg$edge.length.list[[v1]][edge.idx[1]]
            }
        }
    }

    list(
        path = path,
        length = path.length,
        hops = length(path) - 1
    )
}


#' Print Method for path.graph Objects
#'
#' @description
#' Print a concise summary of a path.graph object.
#'
#' @param x A path.graph object
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisible copy of x
#'
#' @examples
#' \dontrun{
#' pg <- create.path.graph(graph, edge.lengths, h = 2)
#' print(pg)
#' }
#'
#' @export
#' @method print path.graph
print.path.graph <- function(x, ...) {
    cat("Path graph object\n")
    cat("  Number of vertices:", length(x$adj.list), "\n")
    cat("  Number of stored paths:", length(x$shortest.paths$paths), "\n")

    # Count edges
    n.edges <- sum(sapply(x$adj.list, length))
    cat("  Number of edges in path graph:", n.edges, "\n")

    invisible(x)
}


#' Summary Method for path.graph Objects
#'
#' @description
#' Provides a detailed summary of a path.graph object including statistics
#' about paths and connectivity.
#'
#' @param object A path.graph object
#' @param ... Additional arguments (currently ignored)
#'
#' @return A list containing summary statistics (invisibly)
#'
#' @examples
#' \dontrun{
#' pg <- create.path.graph(graph, edge.lengths, h = 3)
#' summary(pg)
#' }
#'
#' @export
#' @method summary path.graph
summary.path.graph <- function(object, ...) {
    n.vertices <- length(object$adj.list)
    n.paths <- length(object$shortest.paths$paths)

    cat("Path Graph Summary\n")
    cat("==================\n")
    cat("Number of vertices:", n.vertices, "\n")
    cat("Number of stored paths:", n.paths, "\n")

    if (n.paths > 0) {
        path.lengths <- sapply(object$shortest.paths$paths, length)
        cat("\nPath statistics:\n")
        cat("  Average path length:", mean(path.lengths), "vertices\n")
        cat("  Min path length:", min(path.lengths), "vertices\n")
        cat("  Max path length:", max(path.lengths), "vertices\n")

        # Connectivity analysis
        degree <- sapply(object$adj.list, length)
        cat("\nConnectivity:\n")
        cat("  Average degree:", mean(degree), "\n")
        cat("  Isolated vertices:", sum(degree == 0), "\n")

        # Example paths
        cat("\nExample paths (first 3):\n")
        n.show <- min(3, n.paths)
        for (k in 1:n.show) {
            i <- object$shortest.paths$i[k]
            j <- object$shortest.paths$j[k]
            path <- object$shortest.paths$paths[[k]]
            cat(sprintf("  %d -> %d: %s\n", i, j, paste(path, collapse = " -> ")))
        }
    }

    # Return summary statistics invisibly
    stats <- list(
        n.vertices = n.vertices,
        n.paths = n.paths,
        avg.path.length = if (n.paths > 0) mean(sapply(object$shortest.paths$paths, length)) else NA,
        avg.degree = mean(sapply(object$adj.list, length))
    )

    invisible(stats)
}


#' Create a Series of Path Graphs with Different Hop Limits
#'
#' @description
#' Efficiently constructs multiple path graphs from an input graph, each with a
#' different maximum hop limit. This is more efficient than creating path graphs
#' separately as it reuses computations.
#'
#' @usage
#' create.path.graph.series(graph, edge.lengths, h.values)
#'
#' @param graph A list where each element \code{i} is a numeric vector containing the
#'        indices of vertices adjacent to vertex \code{i}. Indices should be 1-based.
#' @param edge.lengths A list with the same structure as \code{graph} where each
#'        element \code{i} contains the weights/lengths of the edges specified in
#'        \code{graph[[i]]}.
#' @param h.values Integer vector specifying the hop limits for which to compute
#'        path graphs. All values must be >= 1. Values will be sorted in ascending
#'        order internally.
#'
#' @return A list of path.graph objects, one for each value in \code{h.values},
#'         with class "path.graph.series". Each path graph has an additional
#'         attribute "h" storing its hop limit.
#'
#' @details
#' The function computes path graphs for multiple hop limits efficiently by
#' building upon results from smaller hop values. This is particularly useful
#' for analyzing how network connectivity changes with increasing hop limits.
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' graph <- list(
#'   c(2),    # Vertex 1 connected to 2
#'   c(1, 3), # Vertex 2 connected to 1 and 3
#'   c(2)     # Vertex 3 connected to 2
#' )
#' edge.lengths <- list(c(1.0), c(1.0, 2.0), c(2.0))
#'
#' # Create path graphs for hop limits 1, 2, and 3
#' pgs <- create.path.graph.series(graph, edge.lengths, h.values = c(1, 2, 3))
#'
#' # Access individual path graphs
#' pg.h2 <- pgs[[2]]  # Path graph with h=2
#' }
#'
#' @seealso \code{\link{create.path.graph}}, \code{\link{compare.paths}}
#'
#' @export
create.path.graph.series <- function(graph, edge.lengths, h.values) {
    ## Validate basic inputs using same checks as create.path.graph
    if (!is.list(graph)) {
        stop("'graph' must be a list.", call. = FALSE)
    }

    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list.", call. = FALSE)
    }

    if (length(graph) == 0) {
        stop("'graph' must not be empty.", call. = FALSE)
    }

    if (length(graph) != length(edge.lengths)) {
        stop("'graph' and 'edge.lengths' must have the same length.", call. = FALSE)
    }

    # Validate graph structure
    for (i in seq_along(graph)) {
        if (!is.numeric(graph[[i]]) && length(graph[[i]]) > 0) {
            stop(sprintf("graph[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (!is.numeric(edge.lengths[[i]]) && length(edge.lengths[[i]]) > 0) {
            stop(sprintf("edge.lengths[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop(sprintf("graph[[%d]] and edge.lengths[[%d]] must have the same length.", i, i),
                 call. = FALSE)
        }

        if (length(edge.lengths[[i]]) > 0 && any(edge.lengths[[i]] <= 0)) {
            stop(sprintf("All edge lengths in edge.lengths[[%d]] must be positive.", i),
                 call. = FALSE)
        }

        if (length(graph[[i]]) > 0) {
            invalid_idx <- graph[[i]] < 1 | graph[[i]] > length(graph)
            if (any(invalid_idx)) {
                stop(sprintf("Invalid vertex indices in graph[[%d]].", i), call. = FALSE)
            }
        }
    }

    ## Validate h.values
    if (!is.numeric(h.values) || length(h.values) == 0) {
        stop("'h.values' must be a non-empty numeric vector.", call. = FALSE)
    }

    if (any(h.values < 1)) {
        stop("All values in 'h.values' must be at least 1.", call. = FALSE)
    }

    if (anyDuplicated(h.values)) {
        warning("Duplicate values in 'h.values' will be removed.", call. = FALSE)
        h.values <- unique(h.values)
    }

    # Sort h.values for efficient computation
    h.values <- sort(h.values)

    ## Convert graph to 0-based indexing
    graph.0based <- lapply(graph, function(x) {
        if (length(x) == 0) integer(0) else as.integer(x - 1)
    })

    ## Call C++ function
    res <- .Call(S_create_path_graph_series,
                 graph.0based,
                 edge.lengths,
                 as.integer(h.values))

    ## Convert each element to a path.graph object with h attribute
    path.graphs <- mapply(function(pg, h) {
        pg.obj <- new.path.graph(
            adj.list = pg$adj_list,
            edge.length.list = pg$edge_length_list,
            hop.list = pg$hop_list,
            shortest.paths = pg$shortest_paths
        )
        attr(pg.obj, "h") <- h
        pg.obj
    }, res, h.values, SIMPLIFY = FALSE)

    class(path.graphs) <- c("path.graph.series", "list")
    path.graphs
}


#' Print Method for path.graph.series Objects
#'
#' @description
#' Print a summary of a series of path graphs.
#'
#' @param x A path.graph.series object
#' @param ... Additional arguments passed to print.path.graph
#'
#' @return Invisible copy of x
#'
#' @export
#' @method print path.graph.series
print.path.graph.series <- function(x, ...) {
    h.values <- sapply(x, attr, "h")
    cat("Path graph series\n")
    cat("  Number of graphs:", length(x), "\n")
    cat("  Hop limits:", paste(h.values, collapse = ", "), "\n")
    cat("  Number of vertices:", length(x[[1]]$adj.list), "\n\n")

    for (i in seq_along(x)) {
        cat(sprintf("Graph %d (h = %d):\n", i, h.values[i]))
        cat("  Stored paths:", length(x[[i]]$shortest.paths$paths), "\n")
    }

    invisible(x)
}


#' Compare Paths Across Different Hop Limits
#'
#' @description
#' Analyzes how the shortest path between two vertices changes as the hop limit varies.
#'
#' @param x A path.graph.series object
#' @param from Source vertex index (1-based)
#' @param to Target vertex index (1-based)
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{h}{Hop limit}
#'   \item{path_exists}{Logical indicating if a path exists}
#'   \item{path_length}{Total length of the shortest path}
#'   \item{n_hops}{Number of hops in the path}
#'   \item{path}{List column containing the path vertices}
#' }
#'
#' @examples
#' \dontrun{
#' pgs <- create.path.graph.series(graph, edge.lengths, h.values = c(1, 2, 3))
#' compare.paths(pgs, from = 1, to = 3)
#' }
#'
#' @export
compare.paths <- function(x, from, to) {
    if (!inherits(x, "path.graph.series")) {
        stop("'x' must be a path.graph.series object.", call. = FALSE)
    }

    h.values <- sapply(x, attr, "h")

    results <- data.frame(
        h = h.values,
        path_exists = logical(length(x)),
        path_length = numeric(length(x)),
        n_hops = integer(length(x)),
        stringsAsFactors = FALSE
    )

    # Add path as a list column
    results$path <- vector("list", length(x))

    for (i in seq_along(x)) {
        path.info <- get.shortest.path(x[[i]], from, to)
        results$path_exists[i] <- !is.null(path.info)

        if (!is.null(path.info)) {
            results$path_length[i] <- path.info$length
            results$n_hops[i] <- path.info$hops
            results$path[[i]] <- path.info$path
        } else {
            results$path_length[i] <- NA_real_
            results$n_hops[i] <- NA_integer_
        }
    }

    results
}


#' Find Minimum Hop Limit for Path Existence
#'
#' @description
#' Determines the minimum hop limit required for a path to exist between two vertices.
#'
#' @param x A path.graph.series object
#' @param from Source vertex index (1-based)
#' @param to Target vertex index (1-based)
#'
#' @return Integer minimum h value where path exists, or NULL if no path exists
#'         for any h value in the series.
#'
#' @examples
#' \dontrun{
#' pgs <- create.path.graph.series(graph, edge.lengths, h.values = 1:5)
#' min.h <- minh.limit(pgs, from = 1, to = 5)
#' if (!is.null(min.h)) {
#'   cat("Minimum hops needed:", min.h, "\n")
#' }
#' }
#'
#' @export
minh.limit <- function(x, from, to) {
    if (!inherits(x, "path.graph.series")) {
        stop("'x' must be a path.graph.series object.", call. = FALSE)
    }

    for (i in seq_along(x)) {
        if (!is.null(get.shortest.path(x[[i]], from, to))) {
            return(attr(x[[i]], "h"))
        }
    }

    NULL
}


#' Summary Method for path.graph.series Objects
#'
#' @description
#' Provides detailed statistics about a series of path graphs.
#'
#' @param object A path.graph.series object
#' @param ... Additional arguments (currently ignored)
#'
#' @return A data.frame with summary statistics (invisibly)
#'
#' @export
#' @method summary path.graph.series
summary.path.graph.series <- function(object, ...) {
    h.values <- sapply(object, attr, "h")
    n.vertices <- length(object[[1]]$adj.list)

    # Collect statistics
    stats <- data.frame(
        h = h.values,
        n_paths = sapply(object, function(pg) length(pg$shortest.paths$paths)),
        avg_path_length = sapply(object, function(pg) {
            if (length(pg$shortest.paths$paths) == 0) return(NA_real_)
            mean(sapply(pg$shortest.paths$paths, length))
        }),
        max_path_length = sapply(object, function(pg) {
            if (length(pg$shortest.paths$paths) == 0) return(NA_integer_)
            max(sapply(pg$shortest.paths$paths, length))
        }),
        stringsAsFactors = FALSE
    )

    cat("Path Graph Series Summary\n")
    cat("=========================\n")
    cat("Number of vertices:", n.vertices, "\n")
    cat("Hop limits:", paste(h.values, collapse = ", "), "\n\n")
    cat("Statistics by hop limit:\n")
    print(stats, row.names = FALSE)

    invisible(stats)
}


#' Create a Path Length Matrix Graph Structure
#'
#' @description
#' Creates a PLM (Path Length Matrix) graph structure from an adjacency list
#' representation of an undirected graph with edge weights. This structure is
#' optimized for computing path-based statistics and contains precomputed paths
#' of specified length.
#'
#' @usage
#' create.plm.graph(graph, edge.lengths, h)
#'
#' @param graph A list where each element \code{i} is a numeric vector containing the
#'        indices of vertices adjacent to vertex \code{i}. Vertex indices should be
#'        1-based. For undirected graphs, edges should be present in both directions.
#' @param edge.lengths A list of the same structure as \code{graph} where each element
#'        contains the weights/lengths of the corresponding edges in \code{graph}.
#' @param h An odd positive integer specifying the maximum path length (in hops) to be
#'        precomputed. Common values are 3, 5, or 7.
#'
#' @return An S3 object of class "path.graph.plm" containing:
#' \describe{
#'   \item{adj_list}{Adjacency list representation of the graph}
#'   \item{edge_length_list}{List of edge weights}
#'   \item{hop_list}{List of hop counts for paths}
#'   \item{shortest_paths}{Precomputed paths information}
#'   \item{vertex_paths}{Information about paths containing each vertex}
#'   \item{h}{The hop limit used}
#' }
#'
#' @details
#' The PLM (Path Length Matrix) structure is designed for efficient computation
#' of path-based network statistics. It precomputes all paths up to length \code{h}
#' and stores them in a format that allows rapid queries. The restriction to odd
#' values of \code{h} is often used in certain network analysis contexts where
#' paths of odd length have special significance.
#'
#' For undirected graphs, ensure that if edge (i,j) exists, then edge (j,i)
#' also exists in the input with the same weight.
#'
#' @note
#' This function is memory-intensive for large graphs or large values of \code{h},
#' as it stores all paths up to length \code{h}.
#'
#' @examples
#' \dontrun{
#' # Create a simple path graph with 3 vertices: 1 -- 2 -- 3
#' graph <- list(c(2), c(1, 3), c(2))
#' edge.lengths <- list(c(1), c(1, 1), c(1))
#'
#' # Create PLM graph with paths up to length 3
#' plm_graph <- create.plm.graph(graph, edge.lengths, h = 3)
#'
#' print(plm_graph)
#' }
#'
#' @seealso \code{\link{create.path.graph}} for standard path graph creation
#'
#' @export
create.plm.graph <- function(graph, edge.lengths, h) {
    ## Input validation
    if (!is.list(graph) || length(graph) == 0) {
        stop("'graph' must be a non-empty list.", call. = FALSE)
    }

    if (!is.list(edge.lengths) || length(edge.lengths) == 0) {
        stop("'edge.lengths' must be a non-empty list.", call. = FALSE)
    }

    if (length(graph) != length(edge.lengths)) {
        stop("'graph' and 'edge.lengths' must have the same length.", call. = FALSE)
    }

    # Validate each adjacency list
    for (i in seq_along(graph)) {
        if (!is.numeric(graph[[i]]) && length(graph[[i]]) > 0) {
            stop(sprintf("graph[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (!is.numeric(edge.lengths[[i]]) && length(edge.lengths[[i]]) > 0) {
            stop(sprintf("edge.lengths[[%d]] must be numeric or empty.", i), call. = FALSE)
        }

        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop(sprintf("graph[[%d]] and edge.lengths[[%d]] must have the same length.", i, i),
                 call. = FALSE)
        }

        if (length(edge.lengths[[i]]) > 0 && any(edge.lengths[[i]] <= 0)) {
            stop(sprintf("All edge lengths in edge.lengths[[%d]] must be positive.", i),
                 call. = FALSE)
        }

        if (length(graph[[i]]) > 0) {
            invalid_idx <- graph[[i]] < 1 | graph[[i]] > length(graph)
            if (any(invalid_idx)) {
                stop(sprintf("Invalid vertex indices in graph[[%d]]: indices must be between 1 and %d.",
                     i, length(graph)), call. = FALSE)
            }
        }
    }

    ## Validate h
    if (!is.numeric(h) || length(h) != 1 || !is.finite(h)) {
        stop("'h' must be a single finite numeric value.", call. = FALSE)
    }

    h <- as.integer(h)

    if (h < 1) {
        stop("'h' must be at least 1.", call. = FALSE)
    }

    if (h %% 2 == 0) {
        stop("'h' must be odd (1, 3, 5, ...).", call. = FALSE)
    }

    ## Check for undirected graph consistency
    for (i in seq_along(graph)) {
        for (j_idx in seq_along(graph[[i]])) {
            j <- graph[[i]][j_idx]
            # Check if reverse edge exists
            if (!(i %in% graph[[j]])) {
                warning(sprintf("Graph may not be undirected: edge %d->%d exists but %d->%d does not.",
                        i, j, j, i), call. = FALSE)
            }
        }
    }

    ## Convert graph to integer and 0-based indexing
    graph.0based <- lapply(graph, function(x) {
        if (length(x) == 0) integer(0) else as.integer(x - 1)
    })

    ## Call C++ implementation
    res <- .Call(S_create_path_graph_plm,
                 graph.0based,
                 edge.lengths,
                 h)

    res$h <- h
    class(res) <- "path.graph.plm"

    res
}


#' Print Method for path.graph.plm Objects
#'
#' @description
#' Print a summary of a PLM path graph object.
#'
#' @param x A path.graph.plm object
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisible copy of x
#'
#' @export
#' @method print path.graph.plm
print.path.graph.plm <- function(x, ...) {
    cat("PLM Path Graph\n")
    cat("  Number of vertices:", length(x$adj_list), "\n")
    cat("  Maximum path length (h):", x$h, "\n")
    cat("  Number of stored paths:", length(x$shortest_paths$paths), "\n")

    # Count edges
    n.edges <- sum(sapply(x$adj_list, length)) / 2  # Divide by 2 for undirected
    cat("  Number of edges:", n.edges, "\n")

    invisible(x)
}

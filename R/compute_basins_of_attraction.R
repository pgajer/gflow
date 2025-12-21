#' Compute Basins of Attraction for Local Extrema
#'
#' @description
#' This function identifies all local extrema in a weighted graph and computes
#' their basins of attraction. A basin of attraction for a local extremum
#' consists of all vertices reachable via monotone paths: ascending paths for
#' minima and descending paths for maxima.
#'
#' @details
#' The basin of attraction captures the region of influence for each local
#' extremum. For a local minimum, the basin includes all vertices from which
#' the function value strictly increases along every edge of the path leading
#' back to the minimum. Similarly, for a local maximum, the basin includes
#' vertices from which the function value strictly decreases along paths to
#' the maximum.
#'
#' The algorithm performs a breadth-first search from each local extremum,
#' exploring only those edges that maintain the monotonicity property. A vertex
#' belongs to the basin if there exists at least one monotone path connecting
#' it to the extremum. The search continues until no new vertices satisfying
#' the monotonicity constraint can be reached.
#'
#' Each basin structure returned by the function contains the extremum vertex
#' index, its function value, the maximum hop distance within the basin, a
#' complete mapping of all basin vertices to their hop distances from the
#' extremum, and a boundary map identifying vertices adjacent to but outside
#' the basin.
#'
#' The function uses 0-based indexing internally for C++ compatibility but
#' returns results with 1-based indexing following R conventions.
#'
#' @param adj.list A list of integer vectors representing the graph's adjacency
#'   structure. Element \code{i} contains the indices of vertices adjacent to
#'   vertex \code{i}. Indices should be 1-based.
#'
#' @param weight.list A list of numeric vectors containing edge weights.
#'   Element \code{i} contains weights corresponding to the edges in
#'   \code{adj.list[[i]]}. Must have the same structure as \code{adj.list}.
#'
#' @param y A numeric vector of function values at each vertex. The length
#'   must equal the number of vertices (i.e., \code{length(adj.list)}).
#'
#' @param edge.length.quantile.thld Numeric value in (0, 1] specifying the quantile
#'   of the graph's edge length distribution to use as a threshold for basin
#'   construction. This parameter prevents "basin jumping" pathology by excluding
#'   long edges that may cause trajectories to skip over intermediate critical points.
#'
#'   The parameter is always interpreted as a quantile. For example,
#'   \code{edge.length.quantile.thld = 0.75} means that only edges with length
#'   at or below the 75th percentile are used during basin construction; edges
#'   longer than this threshold are excluded.
#'
#'   Common values:
#'   \itemize{
#'     \item \code{0.50}: Use only edges at or below median length (conservative)
#'     \item \code{0.75}: Use edges at or below third quartile (moderate)
#'     \item \code{0.90}: Exclude only the longest 10\% of edges (liberal)
#'     \item \code{1.00}: Use all edges, no constraint (original behavior)
#'   }
#'
#'   During basin construction, an edge is traversed only if both the monotonicity
#'   condition is satisfied AND the edge length does not exceed the computed
#'   threshold. This geometric constraint is essential when the graph contains
#'   edges spanning large distances relative to local neighborhood scales.
#'
#' @note The parameter must satisfy \code{0 < edge.length.quantile.thld <= 1}.
#'   Values outside this range will trigger an error. To disable edge length
#'   filtering entirely, use \code{edge.length.quantile.thld = 1.0}.
#'
#' @param with.trajectories Set to TRUE for the function to return gradient trajectories.
#'
#' @return An object of class \code{"basins_of_attraction"} containing:
#'   \item{lmin_basins}{A list of basin structures for local minima. Each
#'     structure contains:
#'     \itemize{
#'       \item \code{vertex}: The vertex index (1-based)
#'       \item \code{value}: Function value at the vertex
#'       \item \code{hop_idx}: Maximum hop distance within the basin
#'       \item \code{basin_df}: Matrix with columns (vertex, hop_distance)
#'         for all basin members
#'       \item \code{basin_bd_df}: Matrix with columns (vertex, y_value)
#'         for boundary vertices
#'     }}
#'   \item{lmax_basins}{A list of basin structures for local maxima with
#'     the same structure as \code{lmin_basins}}
#'   \item{n_vertices}{Total number of vertices in the graph}
#'   \item{y}{Copy of the input function values}
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' adj.list <- list(
#'   c(2, 3),      # vertex 1 connects to 2, 3
#'   c(1, 3, 4),   # vertex 2 connects to 1, 3, 4
#'   c(1, 2, 4),   # vertex 3 connects to 1, 2, 4
#'   c(2, 3)       # vertex 4 connects to 2, 3
#' )
#'
#' weight.list <- list(
#'   c(1.0, 1.0),
#'   c(1.0, 1.0, 1.0),
#'   c(1.0, 1.0, 1.0),
#'   c(1.0, 1.0)
#' )
#'
#' # Function values with one minimum and one maximum
#' y <- c(5.0, 2.0, 3.0, 6.0)
#'
#' # Compute basins
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Examine structure
#' print(basins)
#'
#' # Access individual basins
#' if (length(basins$lmin_basins) > 0) {
#'   min_basin <- basins$lmin_basins[[1]]
#'   cat("Minimum at vertex:", min_basin$vertex, "\n")
#'   cat("Basin size:", nrow(min_basin$basin_df), "\n")
#' }
#' }
#'
#' @seealso \code{\link{summary.basins_of_attraction}} for generating
#'   summary statistics
#'
#' @export
compute.basins.of.attraction <- function(adj.list,
                                         weight.list,
                                         y,
                                         edge.length.quantile.thld = 0.9,
                                         with.trajectories = FALSE
                                         ) {
    # Input validation
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices")
    }

    ## Validate edge.length.quantile.thld
    if (!is.numeric(edge.length.quantile.thld) || length(edge.length.quantile.thld) != 1) {
        stop("edge.length.quantile.thld must be a single numeric value")
    }

    if (is.na(edge.length.quantile.thld)) {
        stop("edge.length.quantile.thld cannot be NA")
    }

    if (edge.length.quantile.thld <= 0 || edge.length.quantile.thld > 1) {
        stop("edge.length.quantile.thld must be in the range (0, 1]; got ",
             edge.length.quantile.thld)
    }

    if (edge.length.quantile.thld < 0.1) {
        warning("edge.length.quantile.thld is very small (", edge.length.quantile.thld,
                "), which may result in highly fragmented or disconnected basins")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ function
    result <- .Call(S_compute_basins_of_attraction,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.numeric(edge.length.quantile.thld),
                    as.logical(with.trajectories),
                    PACKAGE = "gflow")

    # Add metadata
    result$n_vertices <- length(y)
    result$y <- y

    # Set class for S3 method dispatch
    class(result) <- c("basins_of_attraction", "list")

    return(result)
}

#' Summarize Basins of Attraction
#'
#' @description
#' Generates a comprehensive summary data frame of basins of attraction,
#' including basin characteristics and local density metrics for each extremum.
#'
#' @details
#' This function processes the raw basin structures returned by
#' \code{compute.basins.of.attraction} and produces a summary table suitable
#' for analysis and reporting. The summary includes both basin-specific
#' properties and local geometric characteristics of the graph at each extremum.
#'
#' The density metrics provide information about the local sparsity of the
#' graph around each extremum. The \code{d1} measure represents the distance
#' to the nearest neighbor, while \code{mean.nbrs.dist} captures the average
#' distance to all neighbors. These are converted to percentile ranks
#' (\code{p.d1} and \code{p.mean.nbrs.dist}) to facilitate comparison across
#' vertices. Higher percentile values indicate that an extremum lies in a
#' sparser region of the graph.
#'
#' The hop-k distance measure provides complementary information about isolation
#' at a specified graph distance. By examining vertices at exactly \code{hop_k}
#' steps from each extremum, we assess whether the extremum lies in a locally
#' sparse or dense region at that specific scale. This measure is particularly
#' informative when the immediate neighborhood (captured by \code{mean.nbrs.dist})
#' does not fully characterize the local geometry. The percentile rank
#' \code{p.mean.hopk.dist} enables comparison of this extended neighborhood
#' isolation across all vertices.
#'
#' The basin size quantifies the extent of the extremum's region of influence.
#' Larger basins suggest more prominent features in the function landscape.
#' The hop index indicates the maximum graph distance from the extremum to
#' any vertex within its basin, providing a measure of basin elongation.
#'
#' Extrema are labeled systematically: minima receive lowercase labels
#' (\code{m1}, \code{m2}, ...) in order of increasing function value, while
#' maxima receive uppercase labels (\code{M1}, \code{M2}, ...) in order of
#' decreasing function value. This labeling convention facilitates identification
#' of the most significant extrema.
#'
#' @param object An object of class \code{"basins_of_attraction"} returned by
#'   \code{compute.basins.of.attraction}.
#' @param adj.list A list of integer vectors representing the graph adjacency structure.
#'   Element \code{i} contains the vertex indices of neighbors of vertex \code{i}.
#'   Must have length equal to the number of vertices in the graph.
#' @param edgelen.list A list of numeric vectors containing edge lengths.
#'   Element \code{i} contains the lengths of edges incident to vertex \code{i},
#'   parallel to \code{adj.list[[i]]}. Must have length equal to the number of
#'   vertices in the graph.
#' @param hop.k Integer specifying the hop distance for computing extended
#'   neighborhood isolation. Default is 2. Must be a positive integer.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"basin_summary"} with one row per
#'   local extremum and the following columns:
#'   \item{label}{Character label for the extremum (\code{m1}, \code{m2}, ...
#'     for minima; \code{M1}, \code{M2}, ... for maxima)}
#'   \item{vertex}{Vertex index (1-based)}
#'   \item{value}{Function value at the extremum}
#'   \item{rel.value}{Function value relative to the mean (\code{value / mean(y)})}
#'   \item{type}{Extremum type: \code{"min"} or \code{"max"}}
#'   \item{hop.idx}{Maximum hop distance within the basin}
#'   \item{basin.size}{Number of vertices in the basin}
#'   \item{p.mean.nbrs.dist}{Percentile rank of mean neighbor distance (density measure)}
#'   \item{p.mean.hopk.dist}{Percentile rank of mean hop-k distance (extended isolation measure)}
#'   \item{deg}{Degree of the vertex}
#  ' \item{p.deg}{Percentile rank of degree of the vertex}
#'
#' @examples
#' \dontrun{
#' # Compute basins (see compute.basins.of.attraction examples)
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Create edge length list (Euclidean distances for example)
#' edgelen.list <- lapply(adj.list, function(nbrs) {
#'   rep(1.0, length(nbrs))  # uniform edge lengths
#' })
#'
#' # Generate summary with default hop.k = 2
#' basin_summary <- summary(basins, adj.list, edgelen.list)
#' print(basin_summary)
#'
#' # Generate summary with custom hop distance
#' basin_summary_h3 <- summary(basins, adj.list, edgelen.list, hop.k = 3)
#'
#' # Identify most prominent maximum
#' max_basins <- subset(basin_summary, type == "max")
#' most_prominent <- max_basins[which.max(max_basins$basin.size), ]
#' cat("Most prominent maximum:", most_prominent$label, "\n")
#' cat("Basin size:", most_prominent$basin.size, "\n")
#'
#' # Find extrema in sparse regions (both immediate and extended neighborhoods)
#' sparse_extrema <- subset(basin_summary,
#'                          p.mean.nbrs.dist > 0.75 & p.mean.hopk.dist > 0.75)
#' print(sparse_extrema[, c("label", "type", "value",
#'                          "p.mean.nbrs.dist", "p.mean.hopk.dist")])
#' }
#'
#' @seealso \code{\link{compute.basins.of.attraction}} for computing the basins
#'
#' @export
summary.basins_of_attraction <- function(object, adj.list, edgelen.list, hop.k = 2, ...) {

    if (!inherits(object, "basins_of_attraction")) {
        stop("Input must be of class 'basins_of_attraction'")
    }

    if (length(adj.list) != object$n_vertices) {
        stop("Length of adj.list must equal number of vertices")
    }

    if (length(edgelen.list) != object$n_vertices) {
        stop("Length of edgelen.list must equal number of vertices")
    }

    # Validate that adj.list and edgelen.list are parallel structures
    for (i in seq_len(object$n_vertices)) {
        if (length(adj.list[[i]]) != length(edgelen.list[[i]])) {
            stop(sprintf("Mismatch at vertex %d: adj.list has %d neighbors but edgelen.list has %d lengths",
                        i, length(adj.list[[i]]), length(edgelen.list[[i]])))
        }
    }

    if (!is.numeric(hop.k) || length(hop.k) != 1 || hop.k < 1 || hop.k != as.integer(hop.k)) {
        stop("hop.k must be a positive integer")
    }
    hop.k <- as.integer(hop.k)

    y <- object$y
    n <- object$n_vertices

    # Compute distance metrics
    d1 <- numeric(n)
    mean.nbrs.dist <- numeric(n)
    deg <- numeric(n)

    for (i in seq_len(n)) {
        if (length(edgelen.list[[i]]) > 0) {
            d1[i] <- min(edgelen.list[[i]])
            mean.nbrs.dist[i] <- mean(edgelen.list[[i]])
            deg[i] <- length(edgelen.list[[i]])
        } else {
            d1[i] <- Inf
            mean.nbrs.dist[i] <- Inf
            deg[i] <- 0
        }
    }

    # Compute hop-k distance metric
    mean.hopk.dist <- numeric(n)

    # Compute hop-k neighborhoods for each vertex using BFS
    for (i in seq_len(n)) {
        mean.hopk.dist[i] <- compute_mean_hopk_distance(i, adj.list, edgelen.list, hop.k)
    }

    # Compute density percentiles
    ## p.d1 <- sapply(d1, function(x) sum(d1 <= x) / n)
    p.mean.nbrs.dist <- sapply(mean.nbrs.dist, function(x) sum(mean.nbrs.dist <= x) / n)
    p.mean.hopk.dist <- sapply(mean.hopk.dist, function(x) sum(mean.hopk.dist <= x) / n)
    p.deg <- sapply(deg, function(x) sum(deg >= x) / n)

    # Initialize results
    results <- data.frame(
        label = character(),
        vertex = integer(),
        value = numeric(),
        rel.value = numeric(),
        type = character(),
        hop.idx = numeric(),
        basin.size = integer(),
        ## p.d1 = numeric(),
        p.mean.nbrs.dist = numeric(),
        p.mean.hopk.dist = numeric(),
        deg = numeric(),
        p.deg = numeric(),
        stringsAsFactors = FALSE
    )

    # Process minima basins
    if (length(object$lmin_basins) > 0) {
        min.df <- data.frame(
            label = character(length(object$lmin_basins)),
            vertex = sapply(object$lmin_basins, function(b) b$vertex),
            value = sapply(object$lmin_basins, function(b) b$value),
            rel.value = sapply(object$lmin_basins, function(b) b$value / mean(y)),
            type = "min",
            hop.idx = sapply(object$lmin_basins, function(b) b$hop_idx),
            basin.size = sapply(object$lmin_basins, function(b) nrow(b$basin_df)),
            stringsAsFactors = FALSE
        )

        # Add density metrics
        ## min.df$p.d1 <- p.d1[min.df$vertex]
        min.df$p.mean.nbrs.dist <- p.mean.nbrs.dist[min.df$vertex]
        min.df$p.mean.hopk.dist <- p.mean.hopk.dist[min.df$vertex]
        min.df$deg = deg[min.df$vertex]
        min.df$p.deg = p.deg[min.df$vertex]

        # Sort by value and assign labels
        min.df <- min.df[order(min.df$value), ]
        min.df$label <- paste0("m", seq_len(nrow(min.df)))

        results <- rbind(results, min.df)
    }

    # Process maxima basins
    if (length(object$lmax_basins) > 0) {
        max.df <- data.frame(
            label = character(length(object$lmax_basins)),
            vertex = sapply(object$lmax_basins, function(b) b$vertex),
            value = sapply(object$lmax_basins, function(b) b$value),
            rel.value = sapply(object$lmax_basins, function(b) b$value / mean(y)),
            type = "max",
            hop.idx = sapply(object$lmax_basins, function(b) b$hop_idx),
            basin.size = sapply(object$lmax_basins, function(b) nrow(b$basin_df)),
            stringsAsFactors = FALSE
        )

        # Add density metrics
        ## max.df$p.d1 <- p.d1[max.df$vertex]
        max.df$p.mean.nbrs.dist <- p.mean.nbrs.dist[max.df$vertex]
        max.df$p.mean.hopk.dist <- p.mean.hopk.dist[max.df$vertex]
        max.df$deg = deg[max.df$vertex]
        max.df$p.deg = p.deg[max.df$vertex]

        # Sort by value and assign labels
        max.df <- max.df[order(-max.df$value), ]
        max.df$label <- paste0("M", seq_len(nrow(max.df)))

        results <- rbind(results, max.df)
    }

    rownames(results) <- NULL

    # Return with class for potential print/plot methods
    class(results) <- c("basin_summary", "data.frame")
    return(results)
}

## ============================================================================
## S3 Generics and Methods for Extracting Local Extrema Summaries
## ============================================================================

#' Extract Local Minima from Basin Summary
#'
#' @description
#' Returns the subset of a basin summary data frame containing only local minima.
#'
#' @param object An object of class \code{"basin_summary"} returned by
#'   \code{summary.basins_of_attraction}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"basin_summary"} containing only rows
#'   where \code{type == "min"}, preserving all columns from the original summary.
#'
#' @examples
#' \dontrun{
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' summ <- summary(basins, adj.list, edgelen.list)
#'
#' # Extract only local minima
#' min.summ <- lmin(summ)
#'
#' # Work with minima only
#' print(min.summ)
#' }
#'
#' @seealso \code{\link{lmax.basin_summary}}, \code{\link{summary.basins_of_attraction}}
#'
#' @export
lmin.basin_summary <- function(object, ...) {
    result <- object[object$type == "min", , drop = FALSE]
    rownames(result) <- NULL
    class(result) <- c("basin_summary", "data.frame")
    return(result)
}


#' Extract Local Maxima from Basin Summary
#'
#' @description
#' Returns the subset of a basin summary data frame containing only local maxima.
#'
#' @param object An object of class \code{"basin_summary"} returned by
#'   \code{summary.basins_of_attraction}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"basin_summary"} containing only rows
#'   where \code{type == "max"}, preserving all columns from the original summary.
#'
#' @examples
#' \dontrun{
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' summ <- summary(basins, adj.list, edgelen.list)
#'
#' # Extract only local maxima
#' max.summ <- lmax(summ)
#'
#' # Work with maxima only
#' print(max.summ)
#' }
#'
#' @seealso \code{\link{lmin.basin_summary}}, \code{\link{summary.basins_of_attraction}}
#'
#' @export
lmax.basin_summary <- function(object, ...) {
    result <- object[object$type == "max", , drop = FALSE]
    rownames(result) <- NULL
    class(result) <- c("basin_summary", "data.frame")
    return(result)
}


#' Extract Local Minima
#'
#' @description
#' Generic function for extracting local minima information from an object.
#'
#' @param object An object containing local extrema information.
#' @param ... Additional arguments passed to methods.
#'
#' @return Local minima information, format depends on the method.
#'
#' @seealso \code{\link{lmin.basin_summary}}
#'
#' @export
lmin <- function(object, ...) {
    UseMethod("lmin")
}


#' Extract Local Maxima
#'
#' @description
#' Generic function for extracting local maxima information from an object.
#'
#' @param object An object containing local extrema information.
#' @param ... Additional arguments passed to methods.
#'
#' @return Local maxima information, format depends on the method.
#'
#' @seealso \code{\link{lmax.basin_summary}}
#'
#' @export
lmax <- function(object, ...) {
    UseMethod("lmax")
}

#' Compute mean hop-k distance for a vertex
#'
#' @description
#' Helper function to compute the mean distance from a vertex to all vertices
#' at exactly hop-k distance in the graph.
#'
#' @details
#' We perform breadth-first search to identify all vertices at exactly \code{hop.k}
#' graph distance from the given vertex. For each such vertex found, we compute
#' the geodesic distance by summing edge lengths along the shortest path. The
#' function returns the mean of these distances. If no vertices exist at hop-k
#' distance, the function returns infinity to indicate isolation.
#'
#' @param vertex Integer vertex index (1-based)
#' @param adj.list List of integer vectors representing graph adjacency
#' @param edgelen.list List of edge length vectors parallel to adjacency list
#' @param hop.k Integer hop distance to query
#'
#' @return Numeric value representing mean distance to hop-k neighbors,
#'   or Inf if no vertices exist at that hop distance
#'
#' @keywords internal
compute_mean_hopk_distance <- function(vertex, adj.list, edgelen.list, hop.k) {
    n <- length(adj.list)

    # BFS initialization
    visited <- rep(FALSE, n)
    hop_distance <- rep(-1, n)
    predecessor <- rep(NA_integer_, n)

    # Queue for BFS: stores vertex indices
    queue <- integer(0)

    # Initialize with starting vertex
    visited[vertex] <- TRUE
    hop_distance[vertex] <- 0
    queue <- c(queue, vertex)

    # Perform BFS
    while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]

        # Stop if we've gone beyond hop.k
        if (hop_distance[current] >= hop.k) {
            next
        }

        # Explore neighbors
        neighbors <- adj.list[[current]]
        if (length(neighbors) > 0) {
            for (neighbor in neighbors) {
                if (!visited[neighbor]) {
                    visited[neighbor] <- TRUE
                    hop_distance[neighbor] <- hop_distance[current] + 1
                    predecessor[neighbor] <- current
                    queue <- c(queue, neighbor)
                }
            }
        }
    }

    # Find all vertices at exactly hop.k distance
    hopk_vertices <- which(hop_distance == hop.k)

    if (length(hopk_vertices) == 0) {
        return(Inf)
    }

    # Compute geodesic distances to hop-k vertices by tracing back paths
    distances <- numeric(length(hopk_vertices))

    for (i in seq_along(hopk_vertices)) {
        v <- hopk_vertices[i]
        total_dist <- 0
        current <- v

        # Trace back to source vertex
        while (!is.na(predecessor[current])) {
            prev <- predecessor[current]

            # Find edge length from prev to current
            # Need to find which neighbor index current is in prev's adjacency list
            neighbor_idx <- which(adj.list[[prev]] == current)
            if (length(neighbor_idx) > 0) {
                total_dist <- total_dist + edgelen.list[[prev]][neighbor_idx[1]]
            }

            current <- prev
        }

        distances[i] <- total_dist
    }

    return(mean(distances))
}

#' @export
print.basins_of_attraction <- function(x, ...) {
    cat("Basins of Attraction\n")
    cat("====================\n")
    cat(sprintf("Number of vertices: %d\n", x$n_vertices))
    cat(sprintf("Local minima: %d\n", length(x$lmin_basins)))
    cat(sprintf("Local maxima: %d\n", length(x$lmax_basins)))
    cat("\nUse summary() to generate detailed basin statistics\n")
    invisible(x)
}



#' Draw Basin Vertices in 3D Space
#'
#' Visualizes the vertices belonging to a specific basin from a basins object
#' in 3D using rgl graphics. Searches for the basin in either local minima or
#' local maxima basins and draws the constituent vertices as spheres.
#'
#' @param basins.obj List object containing basin analysis results with structure:
#'   \itemize{
#'     \item \code{basins$lmin_basins}: Named list of local minima basins
#'     \item \code{basins$lmax_basins}: Named list of local maxima basins
#'   }
#'   Each basin contains a \code{basin_df} data frame where the first column
#'   contains vertex indices.
#' @param basin.label Character string identifying the basin to draw. Must match
#'   a name in either \code{lmin_basins} or \code{lmax_basins}
#' @param radius Numeric radius for basin vertex spheres. Default is 0.15
#' @param col Character string specifying color for basin spheres. Default is "cyan"
#'
#' @details
#' The function expects a global \code{graph.3d} object containing 3D coordinates
#' for graph vertices. It searches for \code{basin.label} first in local minima
#' basins, then in local maxima basins. If the basin is not found in either,
#' an error is thrown.
#'
#' @return None. Function is called for its side effect of adding 3D graphics
#'   to the current rgl device.
#'
#' @examples
#' \dontrun{
#' # Assumes graph.3d and final.basins are defined
#' draw.basin(final.basins, "M1", radius = 0.2, col = "cyan")
#' }
#'
#' @export
draw.basin <- function(basins.obj, basin.label, radius = 0.15, col = "cyan") {

    basin <- NULL
    if (basin.label %in% names(basins.obj$basins$lmin_basins)) {
        basin <- basins.obj$basins$lmin_basins[[basin.label]]
    } else if (basin.label %in% names(basins.obj$basins$lmax_basins)) {
        basin <- basins.obj$basins$lmax_basins[[basin.label]]
    } else {
        stop(basin.label, " not found")
    }

    basin.vertices <- basin$basin_df[,1]
    spheres3d(graph.3d[basin.vertices,], radius = radius, col = col)
}


## ============================================================================
## S3 Generic and Method for Basin Extraction
## ============================================================================

#' Extract Basin of Attraction for a Specific Extremum
#'
#' @description
#' Retrieves the basin of attraction associated with a specified local extremum
#' from a basins_of_attraction object. The extremum can be identified either by
#' its vertex index or by its character label when a summary data frame is provided.
#'
#' @details
#' This method provides convenient access to individual basins within a
#' basins_of_attraction object. When working interactively with basin analysis
#' results, users often need to extract specific basins for further examination
#' or visualization. The function supports two identification methods to
#' accommodate different workflows.
#'
#' Direct vertex identification requires knowing the 1-based vertex index of the
#' local extremum. This approach is useful when the extremum of interest was
#' identified through other means, such as domain knowledge or visualization.
#'
#' Label-based identification uses the systematic labels assigned by
#' \code{summary.basins_of_attraction}: lowercase labels (m1, m2, ...) for minima
#' ordered by increasing function value, and uppercase labels (M1, M2, ...) for
#' maxima ordered by decreasing function value. This approach requires passing
#' the summary data frame and is particularly convenient when referencing extrema
#' by their relative prominence.
#'
#' The \code{type} argument controls what information is returned. The "full"
#' option returns the complete basin structure including the extremum vertex,
#' function value, hop index, basin membership with hop distances, and boundary
#' vertices. The "vertices" option returns only the vector of vertex indices
#' belonging to the basin, which is useful for subsetting other data structures
#' or for visualization purposes.
#'
#' @param object An object of class \code{"basins_of_attraction"} returned by
#'   \code{compute.basins.of.attraction}.
#' @param id Either an integer specifying the vertex index (1-based) of a local
#'   extremum, or a character string specifying the extremum label (e.g., "m1",
#'   "M2") when \code{summary.df} is provided.
#' @param summary.df Optional data frame of class \code{"basin_summary"} returned
#'   by \code{summary.basins_of_attraction}. Required when \code{id} is a
#'   character label. Default is \code{NULL}.
#' @param type Character string specifying the return type. Either \code{"full"}
#'   to return the complete basin structure, or \code{"vertices"} to return only
#'   the vector of basin vertex indices. Default is \code{"full"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Depending on \code{type}:
#'   \describe{
#'     \item{type = "full"}{A list containing the complete basin structure:
#'       \code{vertex} (extremum index), \code{value} (function value at extremum),
#'       \code{hop_idx} (maximum hop distance), \code{basin_df} (matrix of basin
#'       vertices and hop distances), and \code{basin_bd_df} (matrix of boundary
#'       vertices and function values).}
#'     \item{type = "vertices"}{An integer vector of 1-based vertex indices
#'       belonging to the basin.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins of attraction
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Extract basin by vertex index (vertex 42 is a local minimum)
#' basin.full <- basin(basins, id = 42)
#' basin.verts <- basin(basins, id = 42, type = "vertices")
#'
#' # Extract basin by label using summary
#' summ <- summary(basins, adj.list, edgelen.list)
#' basin.m1 <- basin(basins, id = "m1", summary.df = summ)
#' basin.M1.verts <- basin(basins, id = "M1", summary.df = summ, type = "vertices")
#'
#' # Use basin vertices for subsetting
#' verts <- basin(basins, id = "m1", summary.df = summ, type = "vertices")
#' y.in.basin <- y[verts]
#' }
#'
#' @seealso \code{\link{compute.basins.of.attraction}} for computing basins,
#'   \code{\link{summary.basins_of_attraction}} for generating basin summaries,
#'   \code{\link{find.basin.by.vertex}}, \code{\link{get.basin.vertex.from.label}}
#'
#' @export
basin.basins_of_attraction <- function(object, id,
                                       summary.df = NULL,
                                       type = c("full", "vertices"),
                                       ...) {
    type <- match.arg(type)

    ## Resolve id to vertex index and extremum type
    if (is.character(id)) {
        ## Label-based lookup: can use summary.df or infer type from label pattern
        extremum.type <- get.extrema.type.from.label(id)

        if (!is.null(summary.df)) {
            ## Use summary.df for lookup
            if (!inherits(summary.df, "basin_summary")) {
                stop("summary.df must be of class 'basin_summary' from summary.basins_of_attraction()")
            }

            label.idx <- which(summary.df$label == id)
            if (length(label.idx) == 0) {
                stop(sprintf("Label '%s' not found in summary.df. Available labels: %s",
                             id, paste(summary.df$label, collapse = ", ")))
            }

            vertex.id <- summary.df$vertex[label.idx]
            extremum.type <- summary.df$type[label.idx]
        } else {
            ## No summary.df: infer type from label and search
            if (is.null(extremum.type)) {
                stop(sprintf("Cannot determine extremum type from label '%s'. ", id),
                     "Provide summary.df or use a standard label format (m1, m2, ..., M1, M2, ...)")
            }
            vertex.id <- NULL  # Will search by matching label pattern
        }
    } else if (is.numeric(id)) {
        ## Direct vertex index lookup
        vertex.id <- as.integer(id)
        if (vertex.id < 1 || vertex.id > object$n_vertices) {
            stop(sprintf("Vertex index %d is out of range [1, %d]",
                         vertex.id, object$n_vertices))
        }
        extremum.type <- NULL  # Will be determined by searching both lists
    } else {
        stop("id must be either an integer vertex index or a character label")
    }

    ## Find the basin in the appropriate list
    basin.result <- NULL

    if (!is.null(vertex.id)) {
        ## Search by vertex id
        if (is.null(extremum.type) || extremum.type == "min") {
            for (i in seq_along(object$lmin_basins)) {
                if (object$lmin_basins[[i]]$vertex == vertex.id) {
                    basin.result <- object$lmin_basins[[i]]
                    break
                }
            }
        }

        if (is.null(basin.result) && (is.null(extremum.type) || extremum.type == "max")) {
            for (i in seq_along(object$lmax_basins)) {
                if (object$lmax_basins[[i]]$vertex == vertex.id) {
                    basin.result <- object$lmax_basins[[i]]
                    break
                }
            }
        }
    } else {
        ## Character label without summary.df: need summary to resolve
        stop("summary.df must be provided when id is a character label")
    }

    ## Check if basin was found
    if (is.null(basin.result)) {
        if (!is.null(vertex.id)) {
            stop(sprintf("Vertex %d is not a local extremum", vertex.id))
        } else {
            stop(sprintf("Basin with id '%s' not found", id))
        }
    }

    ## Return based on type
    if (type == "full") {
        return(basin.result)
    } else {
        ## type == "vertices"
        return(as.integer(basin.result$basin_df[, 1]))
    }
}


#' Extract Basin of Attraction
#'
#' @description
#' Generic function for extracting a basin of attraction from an object.
#'
#' @param object An object containing basin of attraction information.
#' @param ... Additional arguments passed to methods.
#'
#' @return Basin information, format depends on the method.
#'
#' @seealso \code{\link{basin.basins_of_attraction}}
#'
#' @export
basin <- function(object, ...) {
    UseMethod("basin")
}

#' Extract All Basins of Attraction
#'
#' @description
#' Retrieves vertex membership for all basins of attraction, returning a named
#' list where each element contains the vertex indices belonging to that basin.
#'
#' @details
#' This function provides a convenient way to extract all basin memberships
#' simultaneously. In the statistical Morse-Smale framework, basins of attraction
#' are not necessarily disjoint: a vertex may belong to multiple basins when
#' the function landscape contains saddle regions or when noise creates
#' ambiguous gradient directions. The list structure accommodates this
#' overlapping membership.
#'
#' The returned list is named using the systematic labels from the basin summary:
#' lowercase labels (m1, m2, ...) for minima ordered by increasing function value,
#' and uppercase labels (M1, M2, ...) for maxima ordered by decreasing function
#' value. The list ordering follows the summary data frame, with minima appearing
#' before maxima by default.
#'
#' This function is particularly useful for visualization tasks where basin
#' membership needs to be mapped to colors or symbols, for computing overlap
#' statistics between basins, or for subsetting data by basin membership.
#'
#' @param object An object of class \code{"basins_of_attraction"} returned by
#'   \code{compute.basins.of.attraction}.
#' @param summary.df A data frame of class \code{"basin_summary"} returned by
#'   \code{summary.basins_of_attraction}. Required for basin labeling.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list where each element is an integer vector of 1-based
#'   vertex indices belonging to that basin. Names correspond to basin labels
#'   from the summary (e.g., "m1", "m2", "M1", "M2").
#'
#' @examples
#' \dontrun{
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' summ <- summary(basins, adj.list, edgelen.list)
#'
#' # Get all basin memberships
#' all.basins <- basins(basins, summary.df = summ)
#'
#' # Access individual basins by name
#' m1.verts <- all.basins[["m1"]]
#' M1.verts <- all.basins[["M1"]]
#'
#' # Compute basin sizes
#' basin.sizes <- sapply(all.basins, length)
#'
#' # Find vertices in multiple basins (overlap)
#' all.verts <- unlist(all.basins)
#' overlap.verts <- unique(all.verts[duplicated(all.verts)])
#'
#' # Create basin membership indicator matrix
#' n <- length(y)
#' membership <- sapply(all.basins, function(b) seq_len(n) %in% b)
#' }
#'
#' @seealso \code{\link{basin.basins_of_attraction}} for extracting a single basin,
#'   \code{\link{compute.basins.of.attraction}} for computing basins,
#'   \code{\link{summary.basins_of_attraction}} for generating basin summaries
#'
#' @export
basins.basins_of_attraction <- function(object, summary.df, ...) {
    ## Validate inputs
    if (!inherits(object, "basins_of_attraction")) {
        stop("object must be of class 'basins_of_attraction'")
    }
    if (!inherits(summary.df, "basin_summary")) {
        stop("summary.df must be of class 'basin_summary' from summary.basins_of_attraction()")
    }

    ## Initialize result list
    n.basins <- nrow(summary.df)
    result <- vector("list", n.basins)
    names(result) <- summary.df$label

    ## Extract vertices for each basin
    for (i in seq_len(n.basins)) {
        result[[i]] <- basin.basins_of_attraction(
            object,
            id = summary.df$label[i],
            summary.df = summary.df,
            type = "vertices"
        )
    }

    return(result)
}


#' Extract All Basins of Attraction
#'
#' @description
#' Generic function for extracting all basins of attraction from an object.
#'
#' @param object An object containing basins of attraction information.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list of basin memberships, format depends on the method.
#'
#' @seealso \code{\link{basins.basins_of_attraction}}
#'
#' @export
basins <- function(object, ...) {
    UseMethod("basins")
}

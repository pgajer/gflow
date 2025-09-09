#' Perform Harmonic Smoothing on Graph Function Values
#'
#' @description
#' Applies harmonic smoothing to function values defined on vertices of a graph,
#' preserving values at the boundary of a specified region while smoothly
#' interpolating interior values. This function implements a discrete Laplace
#' equation solution using weighted averaging.
#'
#' @details
#' Harmonic smoothing preserves the overall shape of a function defined on a
#' graph while removing local fluctuations. It works by iteratively updating
#' interior vertex values as weighted averages of their neighbors until
#' convergence, while keeping boundary values fixed.
#'
#' The algorithm:
#' \enumerate{
#'   \item Identifies boundary vertices (vertices with neighbors outside the
#'         region or degree 1 vertices)
#'   \item Iteratively updates interior vertex values using edge-weighted averaging
#'   \item Continues until convergence or maximum iterations reached
#' }
#'
#' Edge weights are incorporated by using their inverse as weighting factors,
#' respecting the geometric structure of the graph.
#'
#' @param adj.list A list of integer vectors, where each vector contains indices
#'   of vertices adjacent to the corresponding vertex. Indices must be 1-based.
#' @param weight.list A list of numeric vectors containing weights of edges
#'   corresponding to adjacencies in \code{adj.list}.
#' @param values A numeric vector of function values defined at each vertex.
#' @param region.vertices An integer vector of vertex indices (1-based) defining
#'   the region to be smoothed. Boundary vertices will have fixed values.
#' @param max.iterations Integer scalar, the maximum number of relaxation
#'   iterations to perform. Default is 100.
#' @param tolerance Numeric scalar, the convergence threshold for value changes.
#'   Default is 1e-6.
#'
#' @return A numeric vector of the same length as \code{values}, with smoothed
#'   values within the specified region.
#'
#' @examples
#' \dontrun{
#' # Create a simple grid graph
#' grid.graph <- create.graph.from.grid(10, 10)
#'
#' # Create noisy function values
#' values <- sin(0.1 * seq_len(100)) + rnorm(100, 0, 0.1)
#'
#' # Define a region for smoothing (center of the grid)
#' region <- 35:65
#'
#' # Apply harmonic smoothing
#' smoothed.values <- perform.harmonic.smoothing(
#'   grid.graph$adj.list,
#'   grid.graph$weight.list,
#'   values,
#'   region,
#'   max.iterations = 200,
#'   tolerance = 1e-8
#' )
#'
#' # Plot original vs smoothed values
#' plot(values, type = "l", col = "gray")
#' lines(smoothed.values, col = "red")
#' }
#'
#' @seealso \code{\link{harmonic.smoother}} for smoothing with topology tracking,
#'   \code{\link{get.region.boundary}} for boundary vertex identification
#'
#' @export
perform.harmonic.smoothing <- function(adj.list,
                                       weight.list,
                                       values,
                                       region.vertices,
                                       max.iterations = 100,
                                       tolerance = 1e-6) {
    ## Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }

    if (!is.numeric(values) || !is.vector(values)) {
        stop("'values' must be a numeric vector")
    }

    if (length(values) != length(adj.list)) {
        stop("'values' must have the same length as 'adj.list'")
    }

    if (!is.numeric(region.vertices) || !is.vector(region.vertices)) {
        stop("'region.vertices' must be a numeric vector of vertex indices")
    }

    region.vertices <- as.integer(region.vertices)

    if (length(region.vertices) == 0) {
        stop("'region.vertices' must not be empty")
    }

    if (any(region.vertices < 1) || any(region.vertices > length(adj.list))) {
        stop("'region.vertices' must contain valid vertex indices (1 to ",
             length(adj.list), ")")
    }

    if (!is.numeric(max.iterations) || length(max.iterations) != 1) {
        stop("'max.iterations' must be a single numeric value")
    }

    max.iterations <- as.integer(max.iterations)

    if (max.iterations < 1) {
        stop("'max.iterations' must be a positive integer")
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0) {
        stop("'tolerance' must be a single positive numeric value")
    }

    ## Validate adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (!is.numeric(adj.list[[i]])) {
            stop("All elements of 'adj.list' must be numeric vectors")
        }

        if (!is.numeric(weight.list[[i]])) {
            stop("All elements of 'weight.list' must be numeric vectors")
        }

        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop("Length of adj.list[[", i, "]] and weight.list[[", i,
                 "]] must match")
        }

        if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > length(adj.list))) {
            stop("adj.list[[", i, "]] contains invalid vertex indices")
        }

        if (any(weight.list[[i]] <= 0)) {
            stop("All weights in weight.list[[", i, "]] must be positive")
        }
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call("C_perform_harmonic_smoothing",
                    adj.list.0based,
                    weight.list,
                    as.numeric(values),
                    as.integer(region.vertices),
                    as.integer(max.iterations),
                    as.numeric(tolerance))

    return(result)
}


#' Perform Harmonic Smoothing with Topology Tracking
#'
#' @description
#' Applies harmonic smoothing to function values defined on vertices of a graph
#' while tracking how the topological structure (local extrema and their basins)
#' evolves during the smoothing process. This helps identify the optimal level
#' of smoothing that reduces noise while preserving significant features.
#'
#' @details
#' This function extends standard harmonic smoothing by monitoring the evolution
#' of local extrema during the iterative process. It identifies a "sweet spot"
#' where the topological structure stabilizes, indicating that noise has been
#' removed without over-flattening important features.
#'
#' The algorithm:
#' \enumerate{
#'   \item Iteratively performs harmonic smoothing on interior vertices
#'   \item Periodically identifies local extrema and their basins
#'   \item Monitors the stability of the topological structure
#'   \item Identifies when the topological structure stabilizes
#' }
#'
#' The function returns comprehensive information about the smoothing process,
#' including all intermediate states and the detected stability point.
#'
#' @param adj.list A list of integer vectors, where each vector contains indices
#'     of vertices adjacent to the corresponding vertex. Indices must be
#'     1-based.
#' @param weight.list A list of numeric vectors containing weights of edges
#'     corresponding to adjacencies in \code{adj.list}.
#' @param values A numeric vector of function values defined at each vertex.
#' @param region.vertices An integer vector of vertex indices (1-based) defining
#'     the region to be smoothed. Boundary vertices will have fixed values.
#' @param max.iterations Integer scalar, the maximum number of relaxation
#'     iterations to perform. Default is 100.
#' @param tolerance Numeric scalar, the convergence threshold for value changes.
#'     Default is 1e-6.
#' @param record.frequency Integer scalar, how often to record states (every N
#'     iterations). Default is 1 (record every iteration).
#' @param stability.window Integer scalar, number of consecutive iterations to
#'     check for topological stability. Default is 3.
#' @param stability.threshold Numeric scalar in \eqn{[0,1]}, maximum allowed
#'     difference in topology to consider stable. Default is 0.05.
#'
#' @return A list of class \code{"harmonic_smoother"} containing:
#'   \item{harmonic_predictions}{Numeric vector of smoothed function values}
#'   \item{i_harmonic_predictions}{Matrix of function values at each recorded
#'     iteration (columns are iterations)}
#'   \item{i_basins}{List of matrices representing extrema at each iteration}
#'   \item{stable_iteration}{Integer indicating the iteration at which topology
#'     stabilized}
#'   \item{topology_differences}{Numeric vector of differences between consecutive
#'     recorded iterations}
#'
#' @examples
#' \dontrun{
#' # Create a simple grid graph
#' grid.graph <- create.graph.from.grid(10, 10)
#'
#' # Create noisy function values
#' values <- sin(0.1 * seq_len(100)) + rnorm(100, 0, 0.1)
#'
#' # Define a region for smoothing (center of the grid)
#' region <- 35:65
#'
#' # Apply harmonic smoothing with topology tracking
#' result <- harmonic.smoother(
#'   grid.graph$adj.list,
#'   grid.graph$weight.list,
#'   values,
#'   region,
#'   max.iterations = 200,
#'   tolerance = 1e-8,
#'   record.frequency = 5,  # Record every 5 iterations
#'   stability.window = 3,
#'   stability.threshold = 0.05
#' )
#'
#' # Get the smoothed values
#' smoothed.values <- result$harmonic_predictions
#'
#' # Plot original vs smoothed values
#' plot(values, type = "l", col = "gray")
#' lines(smoothed.values, col = "red")
#'
#' # Plot the evolution of topology differences
#' plot(result$topology_differences, type = "l",
#'      xlab = "Iteration", ylab = "Topology Difference")
#' abline(v = result$stable_iteration, col = "blue", lty = 2)
#' }
#'
#' @seealso \code{\link{perform.harmonic.smoothing}} for basic smoothing,
#'   \code{\link{plot.harmonic_smoother}} for visualization methods,
#'   \code{\link{summary.harmonic_smoother}} for summary statistics
#'
#' @export
harmonic.smoother <- function(adj.list,
                              weight.list,
                              values,
                              region.vertices,
                              max.iterations = 100,
                              tolerance = 1e-6,
                              record.frequency = 1,
                              stability.window = 3,
                              stability.threshold = 0.05) {
    ## Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }

    if (!is.numeric(values) || !is.vector(values)) {
        stop("'values' must be a numeric vector")
    }

    if (length(values) != length(adj.list)) {
        stop("'values' must have the same length as 'adj.list'")
    }

    if (!is.numeric(region.vertices) || !is.vector(region.vertices)) {
        stop("'region.vertices' must be a numeric vector of vertex indices")
    }

    region.vertices <- as.integer(region.vertices)

    if (length(region.vertices) == 0) {
        stop("'region.vertices' must not be empty")
    }

    if (any(region.vertices < 1) || any(region.vertices > length(adj.list))) {
        stop("'region.vertices' must contain valid vertex indices (1 to ",
             length(adj.list), ")")
    }

    if (!is.numeric(max.iterations) || length(max.iterations) != 1) {
        stop("'max.iterations' must be a single numeric value")
    }

    max.iterations <- as.integer(max.iterations)

    if (max.iterations < 1) {
        stop("'max.iterations' must be a positive integer")
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0) {
        stop("'tolerance' must be a single positive numeric value")
    }

    if (!is.numeric(record.frequency) || length(record.frequency) != 1) {
        stop("'record.frequency' must be a single numeric value")
    }

    record.frequency <- as.integer(record.frequency)

    if (record.frequency < 1) {
        stop("'record.frequency' must be a positive integer")
    }

    if (!is.numeric(stability.window) || length(stability.window) != 1) {
        stop("'stability.window' must be a single numeric value")
    }

    stability.window <- as.integer(stability.window)

    if (stability.window < 1) {
        stop("'stability.window' must be a positive integer")
    }

    if (!is.numeric(stability.threshold) || length(stability.threshold) != 1 ||
        stability.threshold < 0 || stability.threshold > 1) {
        stop("'stability.threshold' must be a single numeric value between 0 and 1")
    }

    ## Validate adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (!is.numeric(adj.list[[i]])) {
            stop("All elements of 'adj.list' must be numeric vectors")
        }

        if (!is.numeric(weight.list[[i]])) {
            stop("All elements of 'weight.list' must be numeric vectors")
        }

        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop("Length of adj.list[[", i, "]] and weight.list[[", i,
                 "]] must match")
        }

        if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > length(adj.list))) {
            stop("adj.list[[", i, "]] contains invalid vertex indices")
        }

        if (any(weight.list[[i]] <= 0)) {
            stop("All weights in weight.list[[", i, "]] must be positive")
        }
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call("C_harmonic_smoother",
                    adj.list.0based,
                    weight.list,
                    as.numeric(values),
                    as.integer(region.vertices),
                    as.integer(max.iterations),
                    as.numeric(tolerance),
                    as.integer(record.frequency),
                    as.integer(stability.window),
                    as.numeric(stability.threshold))

    class(result) <- "harmonic_smoother"
    return(result)
}


#' Print Method for Harmonic Smoother Results
#'
#' @description
#' Prints a concise summary of harmonic smoother results.
#'
#' @param x An object of class \code{"harmonic_smoother"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the input object.
#'
#' @seealso \code{\link{harmonic.smoother}}, \code{\link{summary.harmonic_smoother}}
#'
#' @export
print.harmonic_smoother <- function(x, ...) {
  cat("Harmonic Smoother Results\n")
  cat("------------------------\n")
  cat(sprintf("Stable iteration: %d\n", x$stable_iteration))
  cat(sprintf("Number of recorded iterations: %d\n", ncol(x$i_harmonic_predictions)))

  n_extrema <- if (length(x$i_basins) > 0) {
    nrow(x$i_basins[[length(x$i_basins)]])
  } else {
    0
  }

  cat(sprintf("Number of tracked extrema at final iteration: %d\n", n_extrema))
  cat("\nUse summary() for more detailed information.\n")
  invisible(x)
}


#' Summary Method for Harmonic Smoother Results
#'
#' @description
#' Provides detailed summary statistics for harmonic smoother results, including
#' the evolution of extrema counts and basin structure differences.
#'
#' @param object An object of class \code{"harmonic_smoother"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.harmonic_smoother"} containing:
#'   \item{stable_iteration}{The iteration at which topology stabilized}
#'   \item{iterations_recorded}{Total number of recorded iterations}
#'   \item{initial_extrema}{Number of extrema at the first iteration}
#'   \item{initial_maxima}{Number of maxima at the first iteration}
#'   \item{initial_minima}{Number of minima at the first iteration}
#'   \item{final_extrema}{Number of extrema at the final iteration}
#'   \item{final_maxima}{Number of maxima at the final iteration}
#'   \item{final_minima}{Number of minima at the final iteration}
#'   \item{extrema_reduction}{Reduction in number of extrema}
#'   \item{topology_diff_summary}{Summary statistics of topology differences}
#'
#' @seealso \code{\link{harmonic.smoother}}, \code{\link{print.summary.harmonic_smoother}}
#'
#' @export
summary.harmonic_smoother <- function(object, ...) {
  # Calculate evolution of extrema counts
  extrema_counts <- sapply(object$i_basins, nrow)
  maxima_counts <- sapply(object$i_basins, function(b) {
    if (nrow(b) > 0 && "is_max" %in% colnames(b)) {
      sum(b[, "is_max"] == 1)
    } else {
      0
    }
  })
  minima_counts <- extrema_counts - maxima_counts

  # Create result structure
  result <- list(
    stable_iteration = object$stable_iteration,
    iterations_recorded = length(object$i_basins),
    initial_extrema = if (length(extrema_counts) > 0) extrema_counts[1] else 0,
    initial_maxima = if (length(maxima_counts) > 0) maxima_counts[1] else 0,
    initial_minima = if (length(minima_counts) > 0) minima_counts[1] else 0,
    final_extrema = if (length(extrema_counts) > 0) extrema_counts[length(extrema_counts)] else 0,
    final_maxima = if (length(maxima_counts) > 0) maxima_counts[length(maxima_counts)] else 0,
    final_minima = if (length(minima_counts) > 0) minima_counts[length(minima_counts)] else 0,
    extrema_reduction = if (length(extrema_counts) > 0) {
      extrema_counts[1] - extrema_counts[length(extrema_counts)]
    } else {
      0
    },
    topology_diff_summary = if (length(object$topology_differences) > 0) {
      summary(object$topology_differences)
    } else {
      NULL
    }
  )

  class(result) <- "summary.harmonic_smoother"
  return(result)
}


#' Print Method for Summary of Harmonic Smoother Results
#'
#' @description
#' Prints the summary of harmonic smoother results in a formatted manner.
#'
#' @param x An object of class \code{"summary.harmonic_smoother"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the input object.
#'
#' @seealso \code{\link{summary.harmonic_smoother}}
#'
#' @export
print.summary.harmonic_smoother <- function(x, ...) {
  cat("Summary of Harmonic Smoother Results\n")
  cat("-----------------------------------\n")
  cat(sprintf("Stability detected at iteration: %d\n", x$stable_iteration))
  cat(sprintf("Total iterations recorded: %d\n\n", x$iterations_recorded))

  cat("Extrema Evolution:\n")
  cat(sprintf("  Initial: %d total (%d maxima, %d minima)\n",
              x$initial_extrema, x$initial_maxima, x$initial_minima))
  cat(sprintf("  Final:   %d total (%d maxima, %d minima)\n",
              x$final_extrema, x$final_maxima, x$final_minima))

  reduction_pct <- if (x$initial_extrema > 0) {
    100 * x$extrema_reduction / x$initial_extrema
  } else {
    0
  }

  cat(sprintf("  Reduction: %d extrema (%.1f%%)\n\n",
              x$extrema_reduction, reduction_pct))

  if (!is.null(x$topology_diff_summary)) {
    cat("Topology Difference Statistics:\n")
    print(x$topology_diff_summary)
  }

  invisible(x)
}


#' Plot Method for Harmonic Smoother Results
#'
#' @description
#' Creates various plots to visualize the results of harmonic smoothing with
#' topology tracking.
#'
#' @param x An object of class \code{"harmonic_smoother"}.
#' @param y Ignored (included for S3 method consistency).
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @param type Character string specifying the type of plot. Options are:
#'   \describe{
#'     \item{"topology"}{Evolution of topology differences (default)}
#'     \item{"extrema"}{Evolution of extrema counts}
#'     \item{"values"}{Original vs smoothed values}
#'   }
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun{
#' # After running harmonic.smoother()
#' result <- harmonic.smoother(adj.list, weight.list, values, region)
#'
#' # Plot topology evolution
#' plot(result, type = "topology")
#'
#' # Plot extrema counts
#' plot(result, type = "extrema")
#'
#' # Plot smoothed values
#' plot(result, type = "values")
#' }
#'
#' @seealso \code{\link{harmonic.smoother}}
#'
#' @export
plot.harmonic_smoother <- function(x, y = NULL, ..., type = c("topology", "extrema", "values")) {
  type <- match.arg(type)

  if (type == "topology") {
    # Plot the evolution of topology differences
    if (length(x$topology_differences) == 0) {
      warning("No topology differences to plot")
      return(invisible(x))
    }

    plot(x$topology_differences, type = "l",
         xlab = "Iteration", ylab = "Topology Difference",
         main = "Evolution of Topological Structure", ...)
    abline(v = x$stable_iteration, col = "blue", lty = 2)
    legend("topright", legend = "Stability point", col = "blue", lty = 2)

  } else if (type == "extrema") {
    # Plot the evolution of extrema counts
    if (length(x$i_basins) == 0) {
      warning("No basin information to plot")
      return(invisible(x))
    }

    extrema_counts <- sapply(x$i_basins, nrow)
    maxima_counts <- sapply(x$i_basins, function(b) {
      if (nrow(b) > 0 && "is_max" %in% colnames(b)) {
        sum(b[, "is_max"] == 1)
      } else {
        0
      }
    })
    minima_counts <- sapply(x$i_basins, function(b) {
      if (nrow(b) > 0 && "is_max" %in% colnames(b)) {
        sum(b[, "is_max"] == 0)
      } else {
        0
      }
    })

    iterations <- seq_along(extrema_counts)
    plot(iterations, extrema_counts, type = "l", col = "black",
         xlab = "Iteration", ylab = "Count",
         main = "Evolution of Extrema Counts", ...)
    lines(iterations, maxima_counts, col = "red")
    lines(iterations, minima_counts, col = "blue")
    abline(v = x$stable_iteration, col = "green", lty = 2)

    legend("topright",
           legend = c("Total", "Maxima", "Minima", "Stability point"),
           col = c("black", "red", "blue", "green"),
           lty = c(1, 1, 1, 2))

  } else if (type == "values") {
    # Plot original vs smoothed values
    n_values <- nrow(x$i_harmonic_predictions)

    if (is.null(y)) {
      y <- seq_len(n_values)
    }

    # Get original values (first iteration)
    original_values <- x$i_harmonic_predictions[, 1]

    plot(y, original_values, type = "l", col = "gray",
         xlab = "Index", ylab = "Value",
         main = "Original vs Smoothed Values", ...)
    lines(y, x$harmonic_predictions, col = "red", lwd = 2)

    legend("topright",
           legend = c("Original", "Smoothed"),
           col = c("gray", "red"),
           lty = 1,
           lwd = c(1, 2))
  }

  invisible(x)
}


#' Get Boundary Vertices of a Region in a Graph
#'
#' @description
#' Identifies the boundary vertices of a specified region in a graph, defined as
#' vertices in the region that have at least one neighbor outside the region.
#'
#' @details
#' This function determines which vertices in a specified region are positioned at
#' the boundary - meaning they have at least one neighbor that is not part of the
#' region. These boundary vertices are often treated differently in graph-based
#' algorithms, such as harmonic smoothing, where their values are typically fixed
#' as constraints.
#'
#' The boundary definition used matches the one implemented in the C++ harmonic
#' smoothing functions, focusing on vertices that have connections to the outside
#' of the region.
#'
#' @param adj.list A list of integer vectors, where each vector contains the indices
#'   of vertices adjacent to the corresponding vertex. Indices must be 1-based.
#' @param region An integer vector of vertex indices (1-based) defining the region
#'   for which to find boundary vertices.
#'
#' @return An integer vector containing the indices of the boundary vertices,
#'   which is a subset of the input \code{region} vector. Returns an empty
#'   integer vector if no boundary vertices exist.
#'
#' @examples
#' \dontrun{
#' # Create a simple grid graph adjacency list (5x5 grid)
#' create_grid_adj_list <- function(n_rows, n_cols) {
#'   n <- n_rows * n_cols
#'   adj_list <- vector("list", n)
#'   for (i in 1:n_rows) {
#'     for (j in 1:n_cols) {
#'       v <- (i-1) * n_cols + j
#'       neighbors <- numeric(0)
#'
#'       # Add neighbors (up, down, left, right)
#'       if (i > 1) neighbors <- c(neighbors, (i-2) * n_cols + j)  # up
#'       if (i < n_rows) neighbors <- c(neighbors, i * n_cols + j)  # down
#'       if (j > 1) neighbors <- c(neighbors, (i-1) * n_cols + (j-1))  # left
#'       if (j < n_cols) neighbors <- c(neighbors, (i-1) * n_cols + (j+1))  # right
#'
#'       adj_list[[v]] <- neighbors
#'     }
#'   }
#'   return(adj_list)
#' }
#'
#' # Create a 5x5 grid graph
#' grid_adj_list <- create_grid_adj_list(5, 5)
#'
#' # Define a region (central 3x3 subgrid)
#' central_region <- c(7:9, 12:14, 17:19)
#'
#' # Find boundary vertices of the central region
#' boundary <- get.region.boundary(grid_adj_list, central_region)
#' print(boundary)
#' # Expected output: c(7, 8, 9, 12, 14, 17, 18, 19)
#' # (all except the center vertex 13)
#' }
#'
#' @seealso \code{\link{perform.harmonic.smoothing}}, \code{\link{harmonic.smoother}}
#'
#' @export
get.region.boundary <- function(adj.list, region) {
    ## Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }

    if (!is.numeric(region) || !is.vector(region)) {
        stop("'region' must be a numeric vector of vertex indices")
    }

    region <- as.integer(region)

    if (length(region) == 0) {
        return(integer(0))
    }

    ## Check if region indices are valid
    if (any(region < 1) || any(region > length(adj.list))) {
        stop("'region' must contain valid vertex indices (1 to ",
             length(adj.list), ")")
    }

    ## Check for duplicates
    if (anyDuplicated(region)) {
        region <- unique(region)
        warning("Duplicate vertices in 'region' were removed")
    }

    ## Find boundary vertices
    boundary_vertices <- integer(0)

    for (v in region) {
        neighbors <- adj.list[[v]]

        if (!is.numeric(neighbors)) {
            stop("adj.list[[", v, "]] must be a numeric vector")
        }

        ## Check if vertex has any neighbor outside the region
        has_outside_neighbor <- FALSE
        for (neighbor in neighbors) {
            if (!(neighbor %in% region)) {
                has_outside_neighbor <- TRUE
                break
            }
        }

        ## If vertex has a neighbor outside the region, it's a boundary vertex
        if (has_outside_neighbor) {
            boundary_vertices <- c(boundary_vertices, v)
        }
    }

    return(boundary_vertices)
}

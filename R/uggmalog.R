#' Uniform Grid Graph Model-Averaged LOGistic regression (UGGMALOG)
#'
#' @description
#' Implements a sophisticated algorithm for analyzing weighted graphs using
#' local path logistic models with model averaging. The algorithm performs
#' the following main steps:
#' \enumerate{
#'   \item Computes graph diameter and determines bandwidth range
#'   \item Creates a uniform grid representation of the input graph
#'   \item For each candidate bandwidth:
#'     \itemize{
#'       \item Processes paths through grid vertices
#'       \item Fits local logistic models to path data
#'       \item Computes weighted predictions and errors
#'     }
#'   \item Determines optimal bandwidth based on cross-validation errors
#' }
#'
#' The algorithm uses weighted logistic regression on paths through the graph
#' to create local models, which are then combined using weighted averaging.
#' Model evaluation is performed using leave-one-out cross-validation with
#' Brier score errors.
#'
#' @param adj.list A list of integer vectors representing the adjacency list
#'   of the graph. Each element \code{i} contains the indices of vertices
#'   adjacent to vertex \code{i}. Uses 1-based indexing.
#' @param weight.list A list of numeric vectors representing the weights of
#'   the edges. Must have the same structure as \code{adj.list}.
#' @param y A numeric vector of observations at each vertex. Length must
#'   match the number of vertices in the graph.
#' @param best.models.coverage.factor Numeric scalar between 0.5 and 1.0.
#'   Controls the proportion of best models used in model averaging.
#'   Default: 0.9
#' @param min.bw.factor Numeric scalar. Minimum bandwidth factor relative
#'   to graph diameter. Must be positive. Default: 0.05
#' @param max.bw.factor Numeric scalar. Maximum bandwidth factor relative
#'   to graph diameter. Must be greater than \code{min.bw.factor}.
#'   Default: 0.5
#' @param n.bws Positive integer. Number of bandwidths to test between
#'   \code{min.bw.factor} and \code{max.bw.factor}. Default: 50
#' @param grid.size Positive integer. Size of the evaluation grid for
#'   predictions. Default: 100
#' @param start.vertex Positive integer. Index of the starting vertex
#'   (1-based). Must be between 1 and the number of vertices. Default: 1
#' @param snap.tolerance Positive numeric scalar. Tolerance for snapping
#'   distances to grid points. Default: 0.1
#' @param dist.normalization.factor Numeric scalar greater than 1. Factor
#'   for normalizing distances. Default: 1.01
#' @param min.path.size Positive integer. Minimum path size for distance
#'   calculations. Default: 5
#' @param diff.threshold Non-negative integer. Threshold for difference in
#'   path lengths. Default: 5
#' @param kernel.type Integer between 0 and 7. Type of kernel to use:
#'   \itemize{
#'     \item 0: Uniform kernel
#'     \item 1: Triangular kernel
#'     \item 2: Epanechnikov kernel
#'     \item 3: Quartic (biweight) kernel
#'     \item 4: Triweight kernel
#'     \item 5: Tricube kernel
#'     \item 6: Gaussian kernel
#'     \item 7: Cosine kernel (default)
#'   }
#' @param fit.quadratic Logical scalar. Whether to fit quadratic terms in
#'   the local models. Default: FALSE
#' @param max.iterations Positive integer. Maximum number of iterations for
#'   optimization. Default: 100
#' @param ridge.lambda Non-negative numeric scalar. Ridge regression
#'   parameter. Default: 0.0
#' @param tolerance Positive numeric scalar. Convergence tolerance for
#'   optimization. Default: 1e-8
#' @param verbose Logical scalar. Whether to print progress messages.
#'   Default: FALSE
#'
#' @return A list with class \code{"uggmalog"} containing:
#'   \item{candidate_bws}{Numeric vector of candidate bandwidths tested}
#'   \item{bw_predictions}{Numeric matrix of predictions for each bandwidth
#'     (rows: vertices, columns: bandwidths)}
#'   \item{mean_errors}{Numeric vector of mean cross-validation errors for
#'     each bandwidth}
#'   \item{opt_bw_idx}{Integer scalar. Index of the optimal bandwidth
#'     (1-based)}
#'   \item{predictions}{Numeric vector of predictions using the optimal
#'     bandwidth}
#'   \item{graph_diameter}{Numeric scalar. Computed diameter of the graph}
#'
#' @details
#' The algorithm is particularly useful for prediction on graphs where local
#' structure is important. The model averaging approach helps to reduce
#' overfitting and provides more stable predictions.
#'
#' The bandwidth selection is performed automatically using cross-validation,
#' choosing the bandwidth that minimizes the mean Brier score error.
#'
#' @note
#' This function requires compilation of C++ code. The adjacency list uses
#' R's standard 1-based indexing, which is internally converted to 0-based
#' indexing for the C++ implementation.
#'
#' @examples
#' \dontrun{
#' # Create a simple chain graph
#' n <- 10
#' adj <- vector("list", n)
#' weights <- vector("list", n)
#'
#' # Build chain: 1 -- 2 -- 3 -- ... -- n
#' for (i in 1:n) {
#'   adj[[i]] <- integer(0)
#'   weights[[i]] <- numeric(0)
#'
#'   if (i > 1) {
#'     adj[[i]] <- c(adj[[i]], i - 1)
#'     weights[[i]] <- c(weights[[i]], 1.0)
#'   }
#'
#'   if (i < n) {
#'     adj[[i]] <- c(adj[[i]], i + 1)
#'     weights[[i]] <- c(weights[[i]], 1.0)
#'   }
#' }
#'
#' # Generate some response data
#' set.seed(123)
#' y <- sin(seq(0, pi, length.out = n)) + rnorm(n, sd = 0.1)
#'
#' # Run UGGMALOG
#' result <- uggmalog(adj, weights, y, verbose = TRUE)
#'
#' # Print optimal bandwidth index
#' cat("Optimal bandwidth index:", result$opt_bw_idx, "\n")
#'
#' # More complex example: Grid graph
#' # Create a 5x5 grid graph
#' grid_size <- 5
#' n <- grid_size^2
#' adj <- vector("list", n)
#' weights <- vector("list", n)
#'
#' # Function to convert 2D coordinates to 1D index
#' coord_to_idx <- function(i, j) (i - 1) * grid_size + j
#'
#' # Build grid adjacency
#' for (i in 1:grid_size) {
#'   for (j in 1:grid_size) {
#'     idx <- coord_to_idx(i, j)
#'     adj[[idx]] <- integer(0)
#'     weights[[idx]] <- numeric(0)
#'
#'     # Add horizontal edges
#'     if (j > 1) {
#'       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i, j - 1))
#'       weights[[idx]] <- c(weights[[idx]], 1.0)
#'     }
#'     if (j < grid_size) {
#'       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i, j + 1))
#'       weights[[idx]] <- c(weights[[idx]], 1.0)
#'     }
#'
#'     # Add vertical edges
#'     if (i > 1) {
#'       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i - 1, j))
#'       weights[[idx]] <- c(weights[[idx]], 1.0)
#'     }
#'     if (i < grid_size) {
#'       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i + 1, j))
#'       weights[[idx]] <- c(weights[[idx]], 1.0)
#'     }
#'   }
#' }
#'
#' # Generate response based on distance from center
#' center <- (grid_size + 1) / 2
#' y <- numeric(n)
#' for (i in 1:grid_size) {
#'   for (j in 1:grid_size) {
#'     dist_from_center <- sqrt((i - center)^2 + (j - center)^2)
#'     y[coord_to_idx(i, j)] <- exp(-dist_from_center / 2) + rnorm(1, sd = 0.05)
#'   }
#' }
#'
#' # Run UGGMALOG with custom parameters
#' result <- uggmalog(
#'   adj.list = adj,
#'   weight.list = weights,
#'   y = y,
#'   n.bws = 30,
#'   kernel.type = 5,  # Tricube kernel
#'   fit.quadratic = TRUE,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
uggmalog <- function(adj.list,
                     weight.list,
                     y,
                     best.models.coverage.factor = 0.9,
                     min.bw.factor = 0.05,
                     max.bw.factor = 0.5,
                     n.bws = 50L,
                     grid.size = 100L,
                     start.vertex = 1L,
                     snap.tolerance = 0.1,
                     dist.normalization.factor = 1.01,
                     min.path.size = 5L,
                     diff.threshold = 5L,
                     kernel.type = 7L,
                     fit.quadratic = FALSE,
                     max.iterations = 100L,
                     ridge.lambda = 0.0,
                     tolerance = 1e-8,
                     verbose = FALSE) {

  # Input validation with informative error messages

  # Check list inputs
  if (!is.list(adj.list)) {
    stop("'adj.list' must be a list", call. = FALSE)
  }

  if (!is.list(weight.list)) {
    stop("'weight.list' must be a list", call. = FALSE)
  }

  n_vertices <- length(adj.list)

  if (n_vertices == 0) {
    stop("'adj.list' must not be empty", call. = FALSE)
  }

  if (length(weight.list) != n_vertices) {
    stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
  }

  # Check adjacency and weight list structure
  for (i in seq_along(adj.list)) {
    if (!is.numeric(adj.list[[i]]) && !is.integer(adj.list[[i]])) {
      stop(sprintf("Element %d of 'adj.list' must be numeric or integer", i),
           call. = FALSE)
    }

    if (!is.numeric(weight.list[[i]])) {
      stop(sprintf("Element %d of 'weight.list' must be numeric", i),
           call. = FALSE)
    }

    if (length(adj.list[[i]]) != length(weight.list[[i]])) {
      stop(sprintf("Elements %d of 'adj.list' and 'weight.list' must have the same length", i),
           call. = FALSE)
    }

    # Check for valid vertex indices
    if (length(adj.list[[i]]) > 0) {
      if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n_vertices)) {
        stop(sprintf("Invalid vertex indices in element %d of 'adj.list'. Indices must be between 1 and %d",
                     i, n_vertices), call. = FALSE)
      }

      if (any(weight.list[[i]] < 0)) {
        stop(sprintf("Negative weights found in element %d of 'weight.list'", i),
             call. = FALSE)
      }
    }
  }

  # Check response variable
  if (!is.numeric(y)) {
    stop("'y' must be a numeric vector", call. = FALSE)
  }

  if (length(y) != n_vertices) {
    stop(sprintf("Length of 'y' (%d) must match the number of vertices (%d)",
                 length(y), n_vertices), call. = FALSE)
  }

  if (any(is.na(y))) {
    stop("'y' contains NA values", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("'y' contains non-finite values", call. = FALSE)
  }

  # Validate numeric parameters
  if (!is.numeric(best.models.coverage.factor) || length(best.models.coverage.factor) != 1) {
    stop("'best.models.coverage.factor' must be a single numeric value", call. = FALSE)
  }

  if (best.models.coverage.factor <= 0.5 || best.models.coverage.factor > 1.0) {
    stop("'best.models.coverage.factor' must be greater than 0.5 and not greater than 1.0",
         call. = FALSE)
  }

  if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1 || min.bw.factor <= 0) {
    stop("'min.bw.factor' must be a single positive number", call. = FALSE)
  }

  if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1) {
    stop("'max.bw.factor' must be a single numeric value", call. = FALSE)
  }

  if (max.bw.factor <= min.bw.factor) {
    stop("'max.bw.factor' must be greater than 'min.bw.factor'", call. = FALSE)
  }

  # Validate integer parameters
  n.bws <- as.integer(n.bws)
  if (length(n.bws) != 1 || n.bws < 1) {
    stop("'n.bws' must be a positive integer", call. = FALSE)
  }

  grid.size <- as.integer(grid.size)
  if (length(grid.size) != 1 || grid.size < 1) {
    stop("'grid.size' must be a positive integer", call. = FALSE)
  }

  start.vertex <- as.integer(start.vertex)
  if (length(start.vertex) != 1 || start.vertex < 1 || start.vertex > n_vertices) {
    stop(sprintf("'start.vertex' must be an integer between 1 and %d", n_vertices),
         call. = FALSE)
  }

  if (!is.numeric(snap.tolerance) || length(snap.tolerance) != 1 || snap.tolerance <= 0) {
    stop("'snap.tolerance' must be a single positive number", call. = FALSE)
  }

  if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1) {
    stop("'dist.normalization.factor' must be a single numeric value", call. = FALSE)
  }

  if (dist.normalization.factor <= 1) {
    stop("'dist.normalization.factor' must be greater than 1", call. = FALSE)
  }

  min.path.size <- as.integer(min.path.size)
  if (length(min.path.size) != 1 || min.path.size < 1) {
    stop("'min.path.size' must be a positive integer", call. = FALSE)
  }

  diff.threshold <- as.integer(diff.threshold)
  if (length(diff.threshold) != 1 || diff.threshold < 0) {
    stop("'diff.threshold' must be a non-negative integer", call. = FALSE)
  }

  kernel.type <- as.integer(kernel.type)
  if (length(kernel.type) != 1 || kernel.type < 0 || kernel.type > 7) {
    stop("'kernel.type' must be an integer between 0 and 7", call. = FALSE)
  }

  if (!is.logical(fit.quadratic) || length(fit.quadratic) != 1) {
    stop("'fit.quadratic' must be a single logical value", call. = FALSE)
  }

  max.iterations <- as.integer(max.iterations)
  if (length(max.iterations) != 1 || max.iterations < 1) {
    stop("'max.iterations' must be a positive integer", call. = FALSE)
  }

  if (!is.numeric(ridge.lambda) || length(ridge.lambda) != 1 || ridge.lambda < 0) {
    stop("'ridge.lambda' must be a single non-negative number", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0) {
    stop("'tolerance' must be a single positive number", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value", call. = FALSE)
  }

  # Convert to 0-based indexing for C++
  adj.list.0based <- lapply(adj.list, function(x) {
    if (length(x) == 0) {
      integer(0)
    } else {
      as.integer(x - 1)
    }
  })

  # Ensure weight.list elements are numeric
  weight.list <- lapply(weight.list, as.numeric)

  # Call the C++ function
  if (verbose) {
    message("Starting UGGMALOG algorithm...")
    message(sprintf("Graph has %d vertices", n_vertices))
    message(sprintf("Testing %d bandwidths between %.3f and %.3f times graph diameter",
                    n.bws, min.bw.factor, max.bw.factor))
  }

    result <- .Call(S_uggmalog,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    best.models.coverage.factor,
                    min.bw.factor,
                    max.bw.factor,
                    n.bws,
                    grid.size,
                    start.vertex,
                    snap.tolerance,
                    dist.normalization.factor,
                    min.path.size,
                    diff.threshold,
                    kernel.type,
                    fit.quadratic,
                    max.iterations,
                    ridge.lambda,
                    tolerance,
                    verbose)

  # Add class to result
  class(result) <- c("uggmalog", "list")
  
  # Add some additional metadata
  attr(result, "n_vertices") <- n_vertices
  attr(result, "kernel_type") <- kernel.type
  attr(result, "fit_quadratic") <- fit.quadratic

  if (verbose) {
    message(sprintf("Optimal bandwidth: %.3f (index %d)",
                    result$candidate_bws[result$opt_bw_idx],
                    result$opt_bw_idx))
    message("UGGMALOG completed successfully")
  }

  return(result)
}

#' Print method for uggmalog objects
#'
#' @param x An object of class \code{"uggmalog"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#' @export
#' @method print uggmalog
print.uggmalog <- function(x, ...) {
  cat("UGGMALOG Results\n")
  cat("================\n\n")

  n_vertices <- attr(x, "n_vertices")
  kernel_names <- c("Uniform", "Triangular", "Epanechnikov", "Quartic",
                    "Triweight", "Tricube", "Gaussian", "Cosine")
  kernel_type <- attr(x, "kernel_type")

  cat(sprintf("Number of vertices: %d\n", n_vertices))
  cat(sprintf("Graph diameter: %.3f\n", x$graph_diameter))
  cat(sprintf("Kernel type: %s\n", kernel_names[kernel_type + 1]))
  cat(sprintf("Quadratic terms: %s\n", ifelse(attr(x, "fit_quadratic"), "Yes", "No")))
  cat(sprintf("Number of bandwidths tested: %d\n", length(x$candidate_bws)))
  cat(sprintf("Optimal bandwidth: %.3f (index %d)\n",
              x$candidate_bws[x$opt_bw_idx], x$opt_bw_idx))
  cat(sprintf("Mean CV error at optimal bandwidth: %.6f\n",
              x$mean_errors[x$opt_bw_idx]))

  invisible(x)
}

#' Summary method for uggmalog objects
#'
#' @param object An object of class \code{"uggmalog"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return A list with summary information, invisibly
#' @export
#' @method summary uggmalog
summary.uggmalog <- function(object, ...) {
  n_vertices <- attr(object, "n_vertices")

  # Find bandwidth with minimum error
  min_error_idx <- which.min(object$mean_errors)

  # Compute prediction statistics
  pred_stats <- summary(object$predictions)

  # Create summary list
  summ <- list(
    n_vertices = n_vertices,
    graph_diameter = object$graph_diameter,
    n_bandwidths = length(object$candidate_bws),
    bandwidth_range = range(object$candidate_bws),
    optimal_bandwidth = object$candidate_bws[object$opt_bw_idx],
    optimal_bandwidth_idx = object$opt_bw_idx,
    min_cv_error = object$mean_errors[min_error_idx],
    cv_error_range = range(object$mean_errors),
    prediction_summary = pred_stats
  )

  class(summ) <- "summary.uggmalog"
  summ
}

#' Print method for summary.uggmalog objects
#'
#' @param x An object of class \code{"summary.uggmalog"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#' @export
#' @method print summary.uggmalog
print.summary.uggmalog <- function(x, ...) {
  cat("Summary of UGGMALOG Results\n")
  cat("===========================\n\n")

  cat("Graph Information:\n")
  cat(sprintf("  Number of vertices: %d\n", x$n_vertices))
  cat(sprintf("  Graph diameter: %.3f\n", x$graph_diameter))
  cat("\n")

  cat("Bandwidth Selection:\n")
  cat(sprintf("  Number of bandwidths tested: %d\n", x$n_bandwidths))
  cat(sprintf("  Bandwidth range: [%.3f, %.3f]\n",
              x$bandwidth_range[1], x$bandwidth_range[2]))
  cat(sprintf("  Optimal bandwidth: %.3f (index %d)\n",
              x$optimal_bandwidth, x$optimal_bandwidth_idx))
  cat("\n")

  cat("Cross-validation Errors:\n")
  cat(sprintf("  Minimum CV error: %.6f\n", x$min_cv_error))
  cat(sprintf("  CV error range: [%.6f, %.6f]\n",
              x$cv_error_range[1], x$cv_error_range[2]))
  cat("\n")

  cat("Predictions:\n")
  print(x$prediction_summary)

  invisible(x)
}

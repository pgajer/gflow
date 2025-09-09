#' Graph Degree-0 LOWESS with Buffer Zone Cross-Validation
#'
#' @description
#' Performs graph-based degree-0 LOWESS with buffer zone cross-validation for bandwidth selection.
#' This method prevents spatial autocorrelation from biasing performance estimates by
#' creating buffer zones around test vertices during cross-validation.
#'
#' @details
#' This function implements a graph-based extension of LOWESS degree-0 (locally weighted average)
#' with spatially-stratified cross-validation using buffer zones. The algorithm:
#' \enumerate{
#'   \item Creates a maximal packing of vertices to serve as fold seed points
#'   \item Assigns all vertices to the nearest seed point to form spatially coherent folds
#'   \item For each fold, creates a buffer zone around test vertices to ensure spatial independence
#'   \item For each candidate bandwidth, performs cross-validation across the folds
#'   \item Selects the bandwidth with the lowest cross-validation error
#'   \item Fits the final model with the optimal bandwidth
#' }
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of numeric vectors with edge weights corresponding to adjacencies.
#'   Each element \code{weight.list[[i]][j]} is the weight of the edge from vertex i to
#'   \code{adj.list[[i]][j]}.
#' @param y A numeric vector of response values at each vertex in the graph.
#' @param min.bw.factor Minimum bandwidth as a factor of graph diameter. Default is 0.01.
#' @param max.bw.factor Maximum bandwidth as a factor of graph diameter. Default is 0.5.
#' @param n.bws Number of bandwidths to test. Default is 10.
#' @param log.grid Logical, whether to use logarithmic spacing for bandwidth grid. Default is TRUE.
#' @param kernel.type Type of kernel function for weighting. Default is 7 (tricube).
#' @param dist.normalization.factor Factor for normalizing distances in kernel weights. Default is 1.1.
#' @param use.uniform.weights Whether to use uniform weights instead of kernel weights
#' @param buffer.hops Number of hops for buffer zone around test vertices. Default is 2.
#' @param auto.buffer.hops Logical, whether to automatically determine optimal buffer size
#'   based on spatial autocorrelation analysis. Default is TRUE.
#' @param n.folds Number of cross-validation folds. Default is 5.
#' @param with.bw.predictions Logical, whether to compute and store predictions for all bandwidths.
#'   Default is FALSE.
#' @param precision Precision for bandwidth grid computation. Default is 1e-6.
#' @param verbose Logical, whether to print progress information. Default is FALSE.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{predictions}{Vector of optimal predictions for each vertex.}
#'   \item{bw_predictions}{List of prediction vectors for each bandwidth (if \code{with.bw.predictions = TRUE}).}
#'   \item{bw_errors}{Vector of cross-validation errors for each bandwidth.}
#'   \item{bws}{Vector of bandwidth values tested.}
#'   \item{opt_bw}{Optimal bandwidth value.}
#'   \item{opt_bw_idx}{Index of the optimal bandwidth in the bws vector.}
#'   \item{buffer_hops_used}{Number of hops used for the buffer zone (useful when auto.buffer.hops = TRUE).}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a graph and compute graph-based LOWESS with buffer zone cross-validation
#' adj.list <- list(c(2,3), c(1,3,4), c(1,2), c(2))
#' weight.list <- list(c(1,1), c(1,1,1), c(1,1), c(1))
#' y <- c(1.2, 2.3, 0.7, 1.5)
#'
#' # Run with automatic buffer hop determination
#' result <- graph.deg0.lowess.buffer.cv(
#'   adj.list, weight.list, y,
#'   buffer.hops = 1,
#'   auto.buffer.hops = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Check the optimal bandwidth and predictions
#' result$opt_bw
#' result$predictions
#'
#' # See what buffer hop distance was used
#' result$buffer_hops_used
#' }
#'
#' @export
graph.deg0.lowess.buffer.cv <- function(
                                        adj.list,
                                        weight.list,
                                        y,
                                        min.bw.factor = 0.01,
                                        max.bw.factor = 0.5,
                                        n.bws = 20,
                                        log.grid = TRUE,
                                        dist.normalization.factor = 1.1,
                                        use.uniform.weights = FALSE,
                                        buffer.hops = 2,
                                        auto.buffer.hops = TRUE,
                                        kernel.type = 7L,
                                        n.folds = 5,
                                        with.bw.predictions = FALSE,
                                        precision = 1e-6,
                                        verbose = FALSE
                                        ) {
    ## Input validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of 'y' must match the number of vertices in the graph")
    }

    if (!is.numeric(min.bw.factor) || min.bw.factor <= 0) {
        stop("'min.bw.factor' must be a positive number")
    }

    if (!is.numeric(max.bw.factor) || max.bw.factor <= min.bw.factor) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor'")
    }

    if (!is.numeric(n.bws) || n.bws < 1) {
        stop("'n.bws' must be a positive integer")
    }

    if (!is.numeric(buffer.hops) || buffer.hops < 0) {
        stop("'buffer.hops' must be a non-negative integer")
    }

    if (!is.logical(auto.buffer.hops)) {
        stop("'auto.buffer.hops' must be a logical value")
    }

    if (!is.numeric(n.folds) || n.folds < 2) {
        stop("'n.folds' must be an integer greater than 1")
    }

                                        # Convert to 0-based indexing for C++ code
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call(
        "S_graph_deg0_lowess_buffer_cv",
        adj.list.0based,
        weight.list,
        as.numeric(y),
        as.double(min.bw.factor),
        as.double(max.bw.factor),
        as.integer(n.bws),
        as.logical(log.grid),
        as.integer(kernel.type),
        as.double(dist.normalization.factor),
        as.logical(use.uniform.weights),
        as.integer(buffer.hops),
        as.logical(auto.buffer.hops),
        as.integer(n.folds),
        as.logical(with.bw.predictions),
        as.double(precision),
        as.logical(verbose)
    )

                                        # Add class attribute to result for potential S3 methods
    class(result) <- c("graph_deg0_lowess_buffer_cv", "list")

    return(result)
}

#' Print method for graph_deg0_lowess_buffer_cv objects
#'
#' @param x A graph_deg0_lowess_buffer_cv object
#' @param ... Additional arguments (not used)
#'
#' @export
print.graph_deg0_lowess_buffer_cv <- function(x, ...) {
    cat("Graph Degree-0 LOWESS with Buffer Zone Cross-Validation\n")
    cat("------------------------------------------------------\n")
    cat("Number of vertices:", length(x$predictions), "\n")
    cat("Optimal bandwidth:", format(x$opt_bw, digits = 4), "\n")
    cat("Buffer hops used:", x$buffer_hops_used, "\n")
    cat("Number of bandwidths tested:", length(x$bws), "\n")
    cat("Minimum cross-validation error:", format(min(x$bw_errors), digits = 4), "\n")
    cat("\n")
    cat("Use summary() for more detailed information.\n")
}

#' Summary method for graph_deg0_lowess_buffer_cv objects
#'
#' @param object A graph_deg0_lowess_buffer_cv object
#' @param ... Additional arguments (not used)
#'
#' @export
summary.graph_deg0_lowess_buffer_cv <- function(object, ...) {
    result <- list()

                                        # Basic information
    result$n_vertices <- length(object$predictions)
    result$n_bandwidths <- length(object$bws)
    result$opt_bw <- object$opt_bw
    result$opt_bw_idx <- object$opt_bw_idx
    result$buffer_hops_used <- object$buffer_hops_used

                                        # Error statistics
    result$min_error <- min(object$bw_errors)
    result$bw_errors <- object$bw_errors

                                        # Prediction statistics
    result$pred_mean <- mean(object$predictions)
    result$pred_sd <- sd(object$predictions)
    result$pred_range <- range(object$predictions)

                                        # Print summary
    cat("Summary of Graph Degree-0 LOWESS with Buffer Zone Cross-Validation\n")
    cat("--------------------------------------------------------------\n")
    cat("Number of vertices:", result$n_vertices, "\n")
    cat("Buffer hop distance used:", result$buffer_hops_used, "\n\n")

    cat("Bandwidth selection:\n")
    cat("  Optimal bandwidth:", format(result$opt_bw, digits = 4),
        "(index", result$opt_bw_idx, "of", result$n_bandwidths, "tested)\n")
    cat("  Cross-validation error:", format(result$min_error, digits = 6), "\n\n")

    cat("Prediction statistics:\n")
    cat("  Mean:", format(result$pred_mean, digits = 4), "\n")
    cat("  Standard deviation:", format(result$pred_sd, digits = 4), "\n")
    cat("  Range: [", format(result$pred_range[1], digits = 4), ", ",
        format(result$pred_range[2], digits = 4), "]\n", sep = "")

    invisible(result)
}

#' Plot method for graph_deg0_lowess_buffer_cv objects
#'
#' @param x A graph_deg0_lowess_buffer_cv object
#' @param type Plot type: "cv" for cross-validation errors, "predictions" for predictions
#' @param ... Additional arguments passed to plotting functions
#'
#' @export
plot.graph_deg0_lowess_buffer_cv <- function(x, type = "cv", ...) {
    if (type == "cv") {
                                        # Plot cross-validation errors vs bandwidths
        plot(x$bws, x$bw_errors, type = "b", log = "x",
             xlab = "Bandwidth", ylab = "Cross-validation error",
             main = "Bandwidth selection via Buffer Zone CV", ...)

                                        # Highlight optimal bandwidth
        points(x$opt_bw, x$bw_errors[x$opt_bw_idx], col = "red", pch = 19, cex = 1.5)
        abline(v = x$opt_bw, col = "red", lty = 2)

                                        # Add legend
        legend("topright", legend = c("CV errors", "Optimal bandwidth"),
               col = c("black", "red"), pch = c(1, 19), lty = c(1, 2))

    } else if (type == "predictions") {
                                        # Plot predictions vs vertex indices
        plot(seq_along(x$predictions), x$predictions, type = "p",
             xlab = "Vertex index", ylab = "Predicted value",
             main = "LOWESS Predictions with Buffer Zone CV", ...)
    } else {
        stop("Invalid plot type. Use 'cv' or 'predictions'.")
    }
}

#' Print method for summary.graph_deg0_lowess_buffer_cv objects
#'
#' @param x A summary.graph_deg0_lowess_buffer_cv object
#' @param ... Additional arguments (not used)
#'
#' @export
print.summary.graph_deg0_lowess_buffer_cv <- function(x, ...) {
  cat("Summary of Graph Degree-0 LOWESS with Buffer Zone Cross-Validation\n")
  cat("--------------------------------------------------------------\n")
  cat("Number of vertices:", x$n_vertices, "\n")
  cat("Buffer hop distance used:", x$buffer_hops_used, "\n\n")

  cat("Bandwidth selection:\n")
  cat("  Optimal bandwidth:", format(x$opt_bw, digits = 4),
      "(index", x$opt_bw_idx, "of", x$n_bandwidths, "tested)\n")
  cat("  Cross-validation error:", format(x$min_error, digits = 6), "\n\n")

  cat("Prediction statistics:\n")
  cat("  Mean:", format(x$pred_mean, digits = 4), "\n")
  cat("  Standard deviation:", format(x$pred_sd, digits = 4), "\n")
  cat("  Range: [", format(x$pred_range[1], digits = 4), ", ",
      format(x$pred_range[2], digits = 4), "]\n", sep = "")

  invisible(x)
}

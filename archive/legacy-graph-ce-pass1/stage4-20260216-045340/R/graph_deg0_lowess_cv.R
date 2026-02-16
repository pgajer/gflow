#' Graph Degree 0 LOWESS with Cross-Validation for Bandwidth Selection
#'
#' @description Performs local constant fitting (degree 0 LOWESS) on graph data
#' with spatially-stratified cross-validation for automatic bandwidth selection.
#'
#' @details This function implements a graph-based extension of LOWESS degree 0
#' (locally weighted average) with spatially-stratified cross-validation for
#' bandwidth selection. For each vertex, the algorithm:
#' \enumerate{
#'   \item Creates a maximal packing of vertices to serve as fold seed points
#'   \item Assigns all vertices to the nearest seed point to form spatially coherent folds
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
#' @param y A numeric vector of response values for each vertex in the graph.
#' @param min.bw.factor Numeric value specifying the minimum bandwidth as a fraction of
#'   graph diameter (default: 0.05).
#' @param max.bw.factor Numeric value specifying the maximum bandwidth as a fraction of
#'   graph diameter (default: 0.25).
#' @param n.bws Integer specifying the number of candidate bandwidths to evaluate (default: 10).
#' @param log.grid Logical indicating whether to use logarithmic spacing for bandwidth
#'   grid (default: TRUE).
#' @param kernel.type Integer specifying the kernel function for weighting vertices:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#' @param dist.normalization.factor Numeric factor for normalizing distances when calculating
#'   kernel weights (default: 1.1).
#' @param use.uniform.weights Whether to use uniform weights instead of kernel weights
#' @param n.folds Integer specifying the number of cross-validation folds (default: 5).
#' @param with.bw.predictions Logical indicating whether to return predictions for all
#'   bandwidths (default: FALSE).
#' @param precision Numeric value specifying the precision tolerance for bandwidth
#'   grid calculations (default: 0.001).
#' @param verbose Logical indicating whether to display progress information (default: FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{predictions}: Numeric vector of smoothed values using optimal bandwidth
#'   \item \code{bw_predictions}: List of numeric vectors with predictions for each bandwidth
#'     (only if with.bw.predictions=TRUE, otherwise empty)
#'   \item \code{bw_errors}: Numeric vector of cross-validation errors for each bandwidth
#'   \item \code{bws}: Numeric vector of bandwidths tested
#'   \item \code{opt_bw}: Numeric value of the optimal bandwidth
#'   \item \code{opt_bw_idx}: Integer index of the optimal bandwidth in the bws vector
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 100 vertices
#' n <- 100
#' set.seed(123)
#'
#' # Create a ring graph
#' adj.list <- vector("list", n)
#' weight.list <- vector("list", n)
#'
#' for (i in 1:n) {
#'   neighbors <- c(i-1, i+1)
#'   # Handle wrap-around for ring structure
#'   neighbors[neighbors == 0] <- n
#'   neighbors[neighbors == n+1] <- 1
#'
#'   adj.list[[i]] <- neighbors
#'   weight.list[[i]] <- rep(1, length(neighbors))
#' }
#'
#' # Generate response values with spatial pattern
#' y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
#'
#' # Apply graph_deg0_lowess with CV
#' result <- graph.deg0.lowess.cv(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   n.folds = 5,
#'   verbose = TRUE
#' )
#'
#' # Plot results
#' plot(y, type="p", col="gray", pch=16, main="Graph Degree 0 LOWESS with CV")
#' lines(result$predictions, col="red", lwd=2)
#' legend("topright", legend=c("Observed", "Fitted"),
#'        col=c("gray", "red"), pch=c(16, NA), lty=c(NA, 1), lwd=c(NA, 2))
#'
#' # Plot CV errors
#' plot(result$bws, result$bw_errors, type="b", log="x",
#'      xlab="Bandwidth", ylab="CV Error", main="Cross-Validation Error by Bandwidth")
#' abline(v=result$opt_bw, col="red", lty=2)
#' }
#'
#' @export
graph.deg0.lowess.cv <- function(adj.list,
                                 weight.list,
                                 y,
                                 min.bw.factor = 0.05,
                                 max.bw.factor = 0.5,
                                 n.bws = 10,
                                 log.grid = TRUE,
                                 kernel.type = 7L,
                                 dist.normalization.factor = 1.2,
                                 use.uniform.weights = FALSE,
                                 n.folds = 5,
                                 with.bw.predictions = FALSE,
                                 precision = 0.001,
                                 verbose = FALSE) {

    ## Basic parameter validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices in the graph")
    }

    if (!is.numeric(min.bw.factor) || min.bw.factor <= 0)
        stop("min.bw.factor must be positive")

    if (!is.numeric(max.bw.factor) || max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (n.bws < 1)
        stop("n.bws must be positive")

    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 1)
        stop("dist.normalization.factor must be greater than 1")

    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel.type' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    n.folds <- as.integer(n.folds)
    if (n.folds < 2) {
        stop("n.folds must be at least 2")
    }

    if (precision <= 0)
        stop("precision must be positive")

    # Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ implementation
    result <- .Call("S_graph_deg0_lowess_cv",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.numeric(min.bw.factor),
                    as.numeric(max.bw.factor),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.integer(kernel.type),
                    as.numeric(dist.normalization.factor),
                    as.logical(use.uniform.weights),
                    as.integer(n.folds),
                    as.logical(with.bw.predictions),
                    as.numeric(precision),
                    as.logical(verbose))

    # Add class attribute for potential method dispatch
    class(result) <- c("graph_deg0_lowess_cv", "list")
    
    return(result)
}

#' Plot method for graph_deg0_lowess_cv objects
#' 
#' @param x An object of class "graph_deg0_lowess_cv"
#' @param type Type of plot: "fit" for the fitted curve, "cv" for cross-validation errors,
#'   or "both" for a two-panel plot showing both
#' @param ... Additional arguments passed to plotting functions
#' 
#' @return Invisibly returns the input object
#' 
#' @export
plot.graph_deg0_lowess_cv <- function(x, type = "both", ...) {
    if (!inherits(x, "graph_deg0_lowess_cv")) {
        stop("Object must be of class 'graph_deg0_lowess_cv'")
    }
    
    if (type == "both") {
        par(mfrow = c(1, 2))
        plot.graph_deg0_lowess_cv(x, type = "fit", ...)
        plot.graph_deg0_lowess_cv(x, type = "cv", ...)
        par(mfrow = c(1, 1))
        return(invisible(x))
    }
    
    if (type == "fit") {
        y_observed <- seq_along(x$predictions)  # Using index as x-coordinate
        plot(y_observed, x$predictions, 
             type = "l", col = "red", lwd = 2,
             xlab = "Vertex Index", ylab = "Value",
             main = "Graph Degree 0 LOWESS Fit", ...)
        return(invisible(x))
    }
    
    if (type == "cv") {
        plot(x$bws, x$bw_errors, type = "b", log = "x",
             xlab = "Bandwidth", ylab = "CV Error",
             main = "Cross-Validation Error by Bandwidth", ...)
        abline(v = x$opt_bw, col = "red", lty = 2)
        return(invisible(x))
    }
    
    stop("Invalid plot type. Must be one of: 'fit', 'cv', or 'both'")
}

#' Print method for graph_deg0_lowess_cv objects
#' 
#' @param x An object of class "graph_deg0_lowess_cv"
#' @param ... Additional arguments (not used)
#' 
#' @return Invisibly returns the input object
#' 
#' @export
print.graph_deg0_lowess_cv <- function(x, ...) {
    cat("Graph Degree 0 LOWESS with Cross-Validation\n")
    cat("-------------------------------------------\n")
    cat("Optimal bandwidth:", x$opt_bw, "\n")
    cat("Minimum CV error:", min(x$bw_errors), "\n")
    cat("Number of bandwidths evaluated:", length(x$bws), "\n")
    cat("Bandwidth range:", min(x$bws), "to", max(x$bws), "\n")
    cat("\n")
    
    invisible(x)
}

#' Summary method for graph_deg0_lowess_cv objects
#' 
#' @param object An object of class "graph_deg0_lowess_cv"
#' @param ... Additional arguments (not used)
#' 
#' @return A list with summary statistics
#' 
#' @export
summary.graph_deg0_lowess_cv <- function(object, ...) {
    # Calculate summary statistics for predictions
    pred_stats <- list(
        min = min(object$predictions),
        max = max(object$predictions),
        mean = mean(object$predictions),
        median = median(object$predictions),
        sd = sd(object$predictions)
    )
    
    # Create data frame of bandwidths and errors
    bw_results <- data.frame(
        bandwidth = object$bws,
        cv_error = object$bw_errors
    )
    
    # Sort by CV error
    bw_results <- bw_results[order(bw_results$cv_error), ]
    
    result <- list(
        opt_bw = object$opt_bw,
        min_cv_error = min(object$bw_errors),
        prediction_stats = pred_stats,
        bandwidth_results = bw_results
    )
    
    class(result) <- "summary.graph_deg0_lowess_cv"
    return(result)
}

#' Print method for summary.graph_deg0_lowess_cv objects
#' 
#' @param x An object of class "summary.graph_deg0_lowess_cv"
#' @param ... Additional arguments (not used)
#' 
#' @return Invisibly returns the input object
#' 
#' @export
print.summary.graph_deg0_lowess_cv <- function(x, ...) {
    cat("Summary for Graph Degree 0 LOWESS with Cross-Validation\n")
    cat("------------------------------------------------------\n")
    cat("Optimal bandwidth:", x$opt_bw, "\n")
    cat("Minimum CV error:", x$min_cv_error, "\n\n")
    
    cat("Prediction statistics:\n")
    cat("  Min:", x$prediction_stats$min, "\n")
    cat("  Max:", x$prediction_stats$max, "\n")
    cat("  Mean:", x$prediction_stats$mean, "\n")
    cat("  Median:", x$prediction_stats$median, "\n")
    cat("  Standard deviation:", x$prediction_stats$sd, "\n\n")
    
    cat("Top 5 bandwidths by CV error:\n")
    print(head(x$bandwidth_results, 5))
    
    invisible(x)
}

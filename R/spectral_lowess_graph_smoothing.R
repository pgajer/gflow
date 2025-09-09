#' Iterative Spectral LOWESS Graph Smoothing
#'
#' Applies iterative spectral LOWESS (Locally Weighted Scatterplot Smoothing) to
#' graph-structured data. The algorithm iteratively smooths features by estimating
#' conditional expectations using spectral embeddings of local neighborhoods,
#' then reconstructs the graph from the smoothed data.
#'
#' @param adj.list A list of integer vectors representing the adjacency list of
#'   the graph. Each element \code{adj.list[[i]]} contains the indices of vertices
#'   adjacent to vertex \code{i}. Indices must be between 1 and \code{length(adj.list)}.
#'
#' @param weight.list A list of numeric vectors containing edge weights corresponding
#'   to the adjacencies. Each element \code{weight.list[[i]][j]} is the weight of
#'   the edge from vertex \code{i} to vertex \code{adj.list[[i]][j]}. Weights must
#'   be non-negative. Must have the same structure as \code{adj.list}.
#'
#' @param X A numeric matrix where rows represent samples/vertices and columns
#'   represent features. Must have \code{nrow(X) == length(adj.list)}. Missing
#'   values are not allowed.
#'
#' @param max.iterations Maximum number of iterations to perform. Must be a positive
#'   integer. Default: 10.
#'
#' @param convergence.threshold Threshold for convergence. The algorithm stops when
#'   the convergence metric falls below this value. Must be positive. Default: 1e-4.
#'
#' @param convergence.type Type of convergence criterion to use:
#'   \itemize{
#'     \item{1: Maximum absolute difference between successive iterations}
#'     \item{2: Mean absolute difference between successive iterations}
#'     \item{3: Maximum relative change between successive iterations}
#'   }
#'   Default: 1.
#'
#' @param k Number of nearest neighbors for k-NN graph construction. Must be a
#'   positive integer less than the number of vertices. Default: 10.
#'
#' @param pruning.thld Threshold for pruning edges in graph construction. Edges
#'   with weights below this threshold are removed. Must be non-negative. Default: 0.1.
#'
#' @param n.evectors Number of eigenvectors to use for spectral embedding. Must be
#'   a positive integer. Larger values capture more graph structure but increase
#'   computation time. Default: 8.
#'
#' @param n.bws Number of candidate bandwidths to consider for LOWESS smoothing.
#'   Must be a positive integer. Default: 10.
#'
#' @param log.grid Logical. If \code{TRUE}, use logarithmic spacing for the bandwidth
#'   grid; if \code{FALSE}, use linear spacing. Default: \code{TRUE}.
#'
#' @param min.bw.factor Factor for minimum bandwidth, multiplied by the graph diameter
#'   to determine the actual minimum bandwidth. Must be positive and less than
#'   \code{max.bw.factor}. Default: 0.05.
#'
#' @param max.bw.factor Factor for maximum bandwidth, multiplied by the graph diameter
#'   to determine the actual maximum bandwidth. Must be positive and greater than
#'   \code{min.bw.factor}. Default: 0.5.
#'
#' @param dist.normalization.factor Factor for normalizing distances in kernel weight
#'   calculation. Must be at least 1.1. Higher values result in more uniform weights.
#'   Default: 1.1.
#'
#' @param kernel.type Integer specifying the kernel function for weighting:
#'   \itemize{
#'     \item{7: Gaussian kernel (default)}
#'     \item{Other values may be supported depending on the C++ implementation}
#'   }
#'   Default: 7L.
#'
#' @param n.cleveland.iterations Number of robustness iterations for Cleveland's
#'   LOWESS algorithm. Must be a non-negative integer. Default: 1.
#'
#' @param compute.errors Logical. If \code{TRUE}, compute prediction errors for
#'   each iteration. May increase computation time. Default: \code{TRUE}.
#'
#' @param compute.scales Logical. If \code{TRUE}, compute bandwidth/scale information
#'   for each iteration. May increase computation time. Default: \code{TRUE}.
#'
#' @param switch.to.residuals.after Number of iterations to perform direct smoothing
#'   before switching to residual smoothing (boosting mode). Must be a non-negative
#'   integer. Set to 0 to use residual smoothing from the start. If \code{NULL},
#'   defaults to \code{max.iterations} (never switch). Default: \code{NULL}.
#'
#' @param verbose Logical. If \code{TRUE}, print progress information during iterations.
#'   Default: \code{FALSE}.
#'
#' @details
#' The algorithm performs the following steps iteratively:
#' \enumerate{
#'   \item Compute spectral embedding of the graph using eigenvectors of the
#'         normalized Laplacian
#'   \item For each vertex, estimate conditional expectations of features using
#'         LOWESS on the spectral embedding coordinates
#'   \item Reconstruct the graph from the smoothed features using k-NN with
#'         specified pruning
#'   \item Check convergence criterion and stop if threshold is met
#' }
#'
#' The boosting mode (residual smoothing) can be activated by setting
#' \code{switch.to.residuals.after} to a value less than \code{max.iterations}.
#' In this mode, the algorithm smooths residuals from previous iterations rather
#' than the data directly, which can improve convergence in some cases.
#'
#' @return A list of class \code{"spectral.lowess.result"} containing:
#'   \item{smoothed.graphs}{A list of length \code{iterations.performed}, where each
#'         element is a graph (represented as adjacency and weight lists) after
#'         the corresponding iteration}
#'   \item{smoothed.X}{A list of length \code{iterations.performed}, where each
#'         element is the smoothed data matrix after the corresponding iteration}
#'   \item{convergence.metrics}{A numeric vector of length \code{iterations.performed}
#'         containing the convergence metric at each iteration}
#'   \item{iterations.performed}{Integer indicating the number of iterations actually
#'         performed before convergence or reaching \code{max.iterations}}
#'   \item{used.boosting}{Logical indicating whether boosting (residual smoothing)
#'         was used in any iteration}
#'   \item{call}{The matched function call}
#'   \item{parameters}{A list containing all input parameters for reproducibility}
#'
#' @note
#' \itemize{
#'   \item The function requires a C++ implementation accessed via \code{.Call}
#'   \item Large graphs or high-dimensional data may require substantial memory
#'   \item The choice of convergence type affects both speed and quality of results
#'   \item Edge weights should typically be similarity measures (larger = more similar)
#' }
#'
#' @references
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. \emph{Journal of the American Statistical Association},
#' 74(368), 829-836.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#'
#' # Create a k-NN graph
#' library(FNN)
#' knn.result <- get.knn(X, k = 10)
#' adj.list <- lapply(1:n, function(i) knn.result$nn.index[i, ])
#' weight.list <- lapply(1:n, function(i) 1 / (1 + knn.result$nn.dist[i, ]))
#'
#' # Apply spectral LOWESS smoothing
#' result <- spectral.lowess.graph.smoothing(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   X = X,
#'   max.iterations = 5,
#'   convergence.threshold = 1e-3,
#'   verbose = TRUE
#' )
#'
#' # Examine results
#' cat("Iterations performed:", result$iterations.performed, "\n")
#' cat("Final convergence metric:",
#'     result$convergence.metrics[result$iterations.performed], "\n")
#'
#' # Compare original and smoothed data
#' X.smoothed <- result$smoothed.X[[result$iterations.performed]]
#' par(mfrow = c(1, 2))
#' image(X, main = "Original Data", xlab = "Sample", ylab = "Feature")
#' image(X.smoothed, main = "Smoothed Data", xlab = "Sample", ylab = "Feature")
#'
#' # Example with boosting
#' result.boost <- spectral.lowess.graph.smoothing(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   X = X,
#'   max.iterations = 10,
#'   switch.to.residuals.after = 3,
#'   verbose = TRUE
#' )
#'
#' # Compare convergence
#' plot(result$convergence.metrics, type = "b", col = "blue",
#'      xlab = "Iteration", ylab = "Convergence Metric",
#'      main = "Convergence Comparison")
#' lines(result.boost$convergence.metrics, type = "b", col = "red")
#' legend("topright", legend = c("Standard", "Boosting"),
#'        col = c("blue", "red"), lty = 1, pch = 1)
#' }
#'
#' @export
spectral.lowess.graph.smoothing <- function(adj.list,
                                            weight.list,
                                            X,
                                            max.iterations = 10,
                                            convergence.threshold = 1e-4,
                                            convergence.type = 1,
                                            k = 10,
                                            pruning.thld = 0.1,
                                            n.evectors = 8,
                                            n.bws = 10,
                                            log.grid = TRUE,
                                            min.bw.factor = 0.05,
                                            max.bw.factor = 0.5,
                                            dist.normalization.factor = 1.1,
                                            kernel.type = 7L,
                                            n.cleveland.iterations = 1,
                                            compute.errors = TRUE,
                                            compute.scales = TRUE,
                                            switch.to.residuals.after = NULL,
                                            verbose = FALSE) {

    # Store the call for reproducibility
    cl <- match.call()

    ## Input validation with informative error messages

    # Check graph structure
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("'adj.list' and 'weight.list' must be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }

    n.vertices <- length(adj.list)
    if (n.vertices < 2) {
        stop("Graph must have at least 2 vertices")
    }

    # Validate adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (!is.numeric(adj.list[[i]]) || !is.numeric(weight.list[[i]])) {
            stop("All elements of 'adj.list' and 'weight.list' must be numeric")
        }

        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Length mismatch at vertex %d: adj.list has %d neighbors, weight.list has %d weights",
                         i, length(adj.list[[i]]), length(weight.list[[i]])))
        }

        # Check for valid vertex indices
        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n.vertices)) {
                stop(sprintf("Invalid vertex indices in adj.list[[%d]]: indices must be between 1 and %d",
                             i, n.vertices))
            }

            # Check for non-negative weights
            if (any(weight.list[[i]] < 0)) {
                stop(sprintf("Negative weights found in weight.list[[%d]]: all weights must be non-negative",
                             i))
            }
        }
    }

    # Validate data matrix
    if (!is.matrix(X) || !is.numeric(X)) {
        stop("'X' must be a numeric matrix")
    }

    if (anyNA(X)) {
        stop("'X' contains missing values (NA/NaN), which are not allowed")
    }

    if (nrow(X) != n.vertices) {
        stop(sprintf("Number of vertices in graph (%d) does not match number of rows in X (%d)",
                     n.vertices, nrow(X)))
    }

    if (ncol(X) < 1) {
        stop("'X' must have at least one column (feature)")
    }

    # Validate numeric parameters
    if (!is.numeric(max.iterations) || length(max.iterations) != 1 ||
        max.iterations != round(max.iterations) || max.iterations < 1) {
        stop("'max.iterations' must be a single positive integer")
    }
    max.iterations <- as.integer(max.iterations)

    if (!is.numeric(convergence.threshold) || length(convergence.threshold) != 1 ||
        convergence.threshold <= 0) {
        stop("'convergence.threshold' must be a single positive number")
    }

    if (!is.numeric(convergence.type) || length(convergence.type) != 1 ||
        !convergence.type %in% c(1, 2, 3)) {
        stop("'convergence.type' must be 1 (max absolute diff), 2 (mean absolute diff), or 3 (relative change)")
    }
    convergence.type <- as.integer(convergence.type)

    if (!is.numeric(k) || length(k) != 1 || k != round(k) || k < 1 || k >= n.vertices) {
        stop(sprintf("'k' must be a single positive integer less than the number of vertices (%d)",
                     n.vertices))
    }
    k <- as.integer(k)

    if (!is.numeric(pruning.thld) || length(pruning.thld) != 1 || pruning.thld < 0) {
        stop("'pruning.thld' must be a single non-negative number")
    }

    if (!is.numeric(n.evectors) || length(n.evectors) != 1 ||
        n.evectors != round(n.evectors) || n.evectors < 1) {
        stop("'n.evectors' must be a single positive integer")
    }
    n.evectors <- as.integer(n.evectors)

    if (!is.numeric(n.bws) || length(n.bws) != 1 ||
        n.bws != round(n.bws) || n.bws < 1) {
        stop("'n.bws' must be a single positive integer")
    }
    n.bws <- as.integer(n.bws)

    if (!is.logical(log.grid) || length(log.grid) != 1 || is.na(log.grid)) {
        stop("'log.grid' must be a single logical value (TRUE or FALSE)")
    }

    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1 ||
        min.bw.factor <= 0) {
        stop("'min.bw.factor' must be a single positive number")
    }

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1 ||
        max.bw.factor <= 0) {
        stop("'max.bw.factor' must be a single positive number")
    }

    if (min.bw.factor >= max.bw.factor) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor'")
    }

    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1 ||
        dist.normalization.factor < 1.1) {
        stop("'dist.normalization.factor' must be a single number >= 1.1")
    }

    if (!is.numeric(kernel.type) || length(kernel.type) != 1 ||
        kernel.type != round(kernel.type)) {
        stop("'kernel.type' must be a single integer")
    }
    kernel.type <- as.integer(kernel.type)

    if (!is.numeric(n.cleveland.iterations) || length(n.cleveland.iterations) != 1 ||
        n.cleveland.iterations != round(n.cleveland.iterations) || n.cleveland.iterations < 0) {
        stop("'n.cleveland.iterations' must be a single non-negative integer")
    }
    n.cleveland.iterations <- as.integer(n.cleveland.iterations)

    if (!is.logical(compute.errors) || length(compute.errors) != 1 || is.na(compute.errors)) {
        stop("'compute.errors' must be a single logical value (TRUE or FALSE)")
    }

    if (!is.logical(compute.scales) || length(compute.scales) != 1 || is.na(compute.scales)) {
        stop("'compute.scales' must be a single logical value (TRUE or FALSE)")
    }

    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
        stop("'verbose' must be a single logical value (TRUE or FALSE)")
    }

    # Handle switch.to.residuals.after parameter
    if (is.null(switch.to.residuals.after)) {
        switch.to.residuals.after <- max.iterations
    } else {
        if (!is.numeric(switch.to.residuals.after) || length(switch.to.residuals.after) != 1 ||
            switch.to.residuals.after != round(switch.to.residuals.after) ||
            switch.to.residuals.after < 0) {
            stop("'switch.to.residuals.after' must be a single non-negative integer or NULL")
        }
        switch.to.residuals.after <- as.integer(switch.to.residuals.after)

        if (switch.to.residuals.after > max.iterations) {
            warning("'switch.to.residuals.after' is greater than 'max.iterations', ",
                    "boosting will not be used")
            switch.to.residuals.after <- max.iterations
        }
    }

    # Inform about boosting mode
    if (verbose && switch.to.residuals.after < max.iterations) {
        if (switch.to.residuals.after == 0) {
            message("Using boosting (residual smoothing) for all iterations")
        } else {
            message(sprintf("Using direct smoothing for %d iteration%s, then switching to boosting",
                            switch.to.residuals.after,
                            ifelse(switch.to.residuals.after == 1, "", "s")))
        }
    }

    # Convert to 0-based indexing for C++ code
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) > 0) {
            as.integer(x - 1)
        } else {
            integer(0)
        }
    })

    # Store parameters for output
    params <- list(
        max.iterations = max.iterations,
        convergence.threshold = convergence.threshold,
        convergence.type = convergence.type,
        k = k,
        pruning.thld = pruning.thld,
        n.evectors = n.evectors,
        n.bws = n.bws,
        log.grid = log.grid,
        min.bw.factor = min.bw.factor,
        max.bw.factor = max.bw.factor,
        dist.normalization.factor = dist.normalization.factor,
        kernel.type = kernel.type,
        n.cleveland.iterations = n.cleveland.iterations,
        compute.errors = compute.errors,
        compute.scales = compute.scales,
        switch.to.residuals.after = switch.to.residuals.after,
        verbose = verbose
    )

    # Call the C++ function
    if (verbose) {
        message("Starting spectral LOWESS graph smoothing...")
    }

    result <- tryCatch({
        .Call("S_spectral_lowess_graph_smoothing",
              adj.list.0based,
              weight.list,
              X,
              max.iterations,
              convergence.threshold,
              convergence.type,
              k,
              pruning.thld,
              n.evectors,
              n.bws,
              log.grid,
              min.bw.factor,
              max.bw.factor,
              dist.normalization.factor,
              kernel.type,
              n.cleveland.iterations,
              compute.errors,
              compute.scales,
              switch.to.residuals.after,
              verbose)
    }, error = function(e) {
        stop("Error in C++ code: ", conditionMessage(e))
    })
    
    # Add metadata to result
    result$call <- cl
    result$parameters <- params
    
    # Set class for S3 methods
    class(result) <- c("spectral.lowess.result", "list")
    
    if (verbose) {
        message(sprintf("Completed after %d iterations", result$iterations.performed))
    }

    return(result)
}

#' Print method for spectral.lowess.result objects
#'
#' @param x An object of class "spectral.lowess.result"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.spectral.lowess.result <- function(x, ...) {
    cat("Spectral LOWESS Graph Smoothing Result\n")
    cat("--------------------------------------\n")
    cat("Call:\n")
    print(x$call)
    cat("\nIterations performed:", x$iterations.performed, "\n")
    cat("Converged:", x$iterations.performed < x$parameters$max.iterations, "\n")
    if (x$iterations.performed > 0) {
        cat("Final convergence metric:",
            x$convergence.metrics[x$iterations.performed], "\n")
    }
    cat("Used boosting:", x$used.boosting, "\n")
    cat("\nGraph dimensions:\n")
    cat("  Vertices:", length(x$smoothed.graphs[[1]]$adj.list), "\n")
    cat("  Features:", ncol(x$smoothed.X[[1]]), "\n")
    invisible(x)
}

#' Summary method for spectral.lowess.result objects
#'
#' @param object An object of class "spectral.lowess.result"
#' @param ... Additional arguments (ignored)
#'
#' @return A summary object containing key information
#'
#' @export
summary.spectral.lowess.result <- function(object, ...) {
    structure(list(
        call = object$call,
        iterations.performed = object$iterations.performed,
        converged = object$iterations.performed < object$parameters$max.iterations,
        convergence.metrics = object$convergence.metrics,
        used.boosting = object$used.boosting,
        n.vertices = length(object$smoothed.graphs[[1]]$adj.list),
        n.features = ncol(object$smoothed.X[[1]]),
        parameters = object$parameters
    ), class = "summary.spectral.lowess.result")
}

#' Print method for summary.spectral.lowess.result objects
#'
#' @param x An object of class "summary.spectral.lowess.result"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.summary.spectral.lowess.result <- function(x, ...) {
    cat("Summary of Spectral LOWESS Graph Smoothing\n")
    cat("==========================================\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\nConvergence:\n")
    cat("  Iterations performed:", x$iterations.performed, "\n")
    cat("  Converged:", x$converged, "\n")
    cat("  Convergence metrics:",
        paste(round(x$convergence.metrics, 6), collapse = ", "), "\n")
    cat("\nMethod:\n")
    cat("  Boosting used:", x$used.boosting, "\n")
    cat("  Convergence type:",
        c("Max absolute difference", "Mean absolute difference",
          "Relative change")[x$parameters$convergence.type], "\n")
    cat("\nData dimensions:\n")
    cat("  Number of vertices:", x$n.vertices, "\n")
    cat("  Number of features:", x$n.features, "\n")
    cat("\nKey parameters:\n")
    cat("  k (nearest neighbors):", x$parameters$k, "\n")
    cat("  n.evectors:", x$parameters$n.evectors, "\n")
    cat("  pruning.thld:", x$parameters$pruning.thld, "\n")
    invisible(x)
}


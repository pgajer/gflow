#' Iterative Spectral LOWESS Graph Smoothing
#'
#' Apply iterative spectral LOWESS smoothing to a graph and its associated data matrix.
#' At each iteration, conditional expectations of features are estimated using 
#' spectral embedding of local neighborhoods, and a new graph is constructed 
#' from the smoothed data.
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of numeric vectors with edge weights corresponding to adjacencies.
#'   Each element \code{weight.list[[i]][j]} is the weight of the edge from vertex i to
#'   \code{adj.list[[i]][j]}.
#' @param X A numeric matrix where rows are samples and columns are features
#' @param max.iterations Maximum number of iterations to perform
#' @param convergence.threshold Threshold for convergence
#' @param convergence.type Type of convergence criteria:
#'        1 = maximum absolute difference
#'        2 = mean absolute difference
#'        3 = maximum relative change
#' @param k Number of nearest neighbors for kNN graph construction
#' @param pruning.thld Threshold for pruning edges in graph construction
#' @param n.evectors Number of eigenvectors for spectral embedding
#' @param n.bws Number of candidate bandwidths for LOWESS
#' @param log.grid Logical, whether to use logarithmic spacing for bandwidth grid
#' @param min.bw.factor Factor for minimum bandwidth (multiplied by graph diameter)
#' @param max.bw.factor Factor for maximum bandwidth (multiplied by graph diameter)
#' @param dist.normalization.factor Factor for normalizing distances in kernel weight calculation
#' @param kernel.type Type of kernel function (1 = Gaussian, 2 = Exponential, etc.)
#' @param n.cleveland.iterations Number of robustness iterations for Cleveland's algorithm
#' @param compute.errors Logical, whether to compute prediction errors
#' @param compute.scales Logical, whether to compute bandwidth/scale information
#' @param switch.to.residuals.after Number of iterations to perform direct smoothing before
#'        switching to residual smoothing (boosting mode). Default is max.iterations (never switch).
#'        Set to 0 to use residual smoothing from the start.
#' @param verbose Logical, whether to print progress information
#'
#' @return A list containing:
#'   \item{smoothed.graphs}{List of smoothed graphs at each iteration}
#'   \item{smoothed.X}{List of smoothed data matrices at each iteration}
#'   \item{convergence.metrics}{Numeric vector of convergence metrics at each iteration}
#'   \item{iterations.performed}{Number of iterations actually performed}
#'   \item{used.boosting}{Logical, whether boosting (residual smoothing) was used}
#'
#' @examples
#' \dontrun{
#' # Create a graph and data matrix
#' graph <- create.iknn.graph(X, k = 10, pruning.thld = 0.1)
#' 
#' # Apply spectral LOWESS graph smoothing - traditional approach
#' result1 <- spectral.lowess.graph.smoothing(
#'   adj.list = graph$adj_list,
#'   weight.list = graph$weight_list,
#'   X = X,
#'   max.iterations = 10,
#'   convergence.threshold = 1e-4,
#'   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
#'   k = 10,
#'   pruning.thld = 0.1,
#'   n.evectors = 8,
#'   verbose = TRUE
#' )
#' 
#' # Apply spectral LOWESS graph smoothing with boosting
#' result2 <- spectral.lowess.graph.smoothing(
#'   adj.list = graph$adj_list,
#'   weight.list = graph$weight_list,
#'   X = X,
#'   max.iterations = 10,
#'   convergence.threshold = 1e-4,
#'   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
#'   k = 10,
#'   pruning.thld = 0.1,
#'   n.evectors = 8,
#'   switch.to.residuals.after = 2,  # Switch to boosting after 2 iterations
#'   verbose = TRUE
#' )
#' 
#' # Access final smoothed data matrix
#' X.smoothed <- result2$smoothed.X[[length(result2$smoothed.X)]]
#' 
#' # Plot convergence metrics for both approaches
#' plot(result1$convergence.metrics, type = "b", col = "blue",
#'      xlab = "Iteration", ylab = "Convergence Metric")
#' lines(result2$convergence.metrics, type = "b", col = "red")
#' legend("topright", legend = c("Traditional", "Boosting"), 
#'        col = c("blue", "red"), lty = 1)
#' }
#'
#' @export
spectral.lowess.graph.smoothing <- function(
                                            adj.list,
                                            weight.list,
                                            X,
                                            max.iterations = 10,
                                            convergence.threshold = 1e-4,
                                            convergence.type = 1,  # MAX.ABSOLUTE.DIFF
                                            k = 10,
                                            pruning.thld = 0.1,
                                            n.evectors = 8,
                                            n.bws = 10,
                                            log.grid = TRUE,
                                            min.bw.factor = 0.05,
                                            max.bw.factor = 0.5,
                                            dist.normalization.factor = 1.1,
                                            kernel.type = 7L,  # Gaussian
                                            n.cleveland.iterations = 1,
                                            compute.errors = TRUE,
                                            compute.scales = TRUE,
                                            switch.to.residuals.after = NULL,
                                            verbose = FALSE
                                            ) {
    ## Input validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (!is.matrix(X) || !is.numeric(X)) {
        stop("'X' must be a numeric matrix")
    }

    if (length(adj.list) != nrow(X)) {
        stop("Number of vertices in graph (", length(adj.list),
             ") does not match number of samples in X (", nrow(X), ")")
    }

    if (!is.numeric(max.iterations) || max.iterations < 1) {
        stop("'max.iterations' must be a positive integer")
    }

    if (!is.numeric(convergence.threshold) || convergence.threshold <= 0) {
        stop("'convergence.threshold' must be a positive number")
    }

    if (!convergence.type %in% c(1, 2, 3)) {
        stop("'convergence.type' must be 1 (MAX.ABSOLUTE.DIFF), 2 (MEAN.ABSOLUTE.DIFF), or 3 (RELATIVE.CHANGE)")
    }
    
    # If switch.to.residuals.after is NULL, default to max.iterations (never switch)
    if (is.null(switch.to.residuals.after)) {
        switch.to.residuals.after <- max.iterations
    }
    
    if (!is.numeric(switch.to.residuals.after) || switch.to.residuals.after < 0) {
        stop("'switch.to.residuals.after' must be a non-negative integer")
    }
    
    # Warn if boosting will be used
    if (switch.to.residuals.after < max.iterations) {
        if (verbose) {
            if (switch.to.residuals.after == 0) {
                message("Using boosting (residual smoothing) for all iterations")
            } else {
                message("Using direct smoothing for ", switch.to.residuals.after, 
                        " iterations, then switching to boosting (residual smoothing)")
            }
        }
    }

    if (max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (dist.normalization.factor < 1.1)
        stop("dist.normalization.factor must be greater than or equal to 1.1")

    # Convert to 0-based indexing for C++ code
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call(
        "S_spectral_lowess_graph_smoothing",
        adj.list.0based,
        weight.list,
        X,
        as.integer(max.iterations),
        as.double(convergence.threshold),
        as.integer(convergence.type),
        as.integer(k),
        as.double(pruning.thld),
        as.integer(n.evectors),
        as.integer(n.bws),
        as.logical(log.grid),
        as.double(min.bw.factor),
        as.double(max.bw.factor),
        as.double(dist.normalization.factor),
        as.integer(kernel.type),
        as.integer(n.cleveland.iterations),
        as.logical(compute.errors),
        as.logical(compute.scales),
        as.integer(switch.to.residuals.after),
        as.logical(verbose)
    )

    return(result)
}

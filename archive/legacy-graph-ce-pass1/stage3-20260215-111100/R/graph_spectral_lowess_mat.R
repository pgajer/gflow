#' Matrix version of Spectral Graph Local Polynomial Regression
#'
#' @description
#' Performs locally weighted regression for multiple response variables simultaneously
#' on a graph structure using spectral embeddings. This function is more efficient than
#' applying the standard `graph_spectral_lowess` to each response variable separately, as
#' it computes the graph Laplacian and eigenvectors only once.
#'
#' @param adj.list List where each element contains indices of adjacent vertices
#' @param weight.list List where each element contains weights of adjacent edges
#' @param Y Matrix or list of numeric vectors where each column/element represents
#'          a response variable at each vertex
#' @param n.evectors Number of eigenvectors to use for spectral embedding
#' @param n.bws Number of candidate bandwidths to evaluate (default: 20)
#' @param log.grid Logical; whether to use logarithmic spacing for bandwidth grid (default: FALSE)
#' @param min.bw.factor Factor for minimum bandwidth as fraction of graph diameter (default: 0.01)
#' @param max.bw.factor Factor for maximum bandwidth as fraction of graph diameter (default: 0.5)
#' @param dist.normalization.factor Factor for normalizing distances (default: 1.1)
#' @param kernel.type Integer indicating kernel type: 0=gaussian, 1=epanechnikov, etc. (default: 7)
#' @param precision Precision threshold for numerical calculations (default: 1e-6)
#' @param n.cleveland.iterations Number of iterations for Cleveland's robust fitting (default: 3)
#' @param with.errors Logical; whether to compute prediction errors (default: FALSE)
#' @param with.scale Logical; whether to compute bandwidth scales (default: FALSE)
#' @param verbose Logical; whether to print progress information (default: TRUE)
#'
#' @return A list with components:
#'   \item{predictions}{Matrix where each column contains smoothed values for a response variable}
#'   \item{errors}{Matrix of prediction errors (if with.errors=TRUE)}
#'   \item{scale}{Matrix of local bandwidth scales (if with.scale=TRUE)}
#'
#' @details
#' This function efficiently handles multiple response variables by leveraging the fact that
#' the graph structure processing (Laplacian calculation, spectral decomposition) only needs
#' to be performed once. For each vertex, it identifies an appropriate neighborhood, creates
#' a spectral embedding, and fits robust linear models for each response variable within 
#' this embedding space.
#'
#' The optimal bandwidth is selected independently for each response variable at each vertex,
#' balancing between overfitting (small bandwidth) and underfitting (large bandwidth).
#'
#' @examples
#' \dontrun{
#' # Create sample graph
#' n <- 100
#' g <- make.example.graph(n)
#' 
#' # Create multiple response variables (e.g., 3 columns)
#' Y <- matrix(rnorm(n*3), ncol=3)
#' 
#' # Apply spectral lowess to all columns at once
#' result <- graph.spectral.lowess.mat(g$adj.list, g$weight.list, Y)
#' 
#' # Plot results for first response variable
#' plot(Y[,1], result$predictions[,1], main="Original vs. Smoothed (Var 1)")
#' abline(0, 1, col="red", lty=2)
#' }
#'
#' @export
graph.spectral.lowess.mat <- function(adj.list,
                                      weight.list,
                                      Y,
                                      n.evectors = 5,
                                      n.bws = 20,
                                      log.grid = FALSE,
                                      min.bw.factor = 0.01,
                                      max.bw.factor = 0.5,
                                      dist.normalization.factor = 1.1,
                                      kernel.type = 7L,
                                      precision = 1e-6,
                                      n.cleveland.iterations = 0L,
                                      with.errors = FALSE,
                                      with.scale = FALSE,
                                      verbose = TRUE) {
    
    # Input validation
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }
    
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }
    
    # Check Y and convert if necessary
    if (is.matrix(Y)) {
        if (nrow(Y) != length(adj.list)) {
            stop("Number of rows in Y must match number of vertices in the graph")
        }
    } else if (is.list(Y)) {
        if (length(Y) == 0) {
            stop("Y list cannot be empty")
        }
        
        if (length(Y[[1]]) != length(adj.list)) {
            stop("Length of each vector in Y must match number of vertices in the graph")
        }
        
        # Verify all elements have the same length
        lengths <- sapply(Y, length)
        if (any(lengths != lengths[1])) {
            stop("All vectors in Y list must have the same length")
        }
    } else {
        stop("Y must be a matrix or a list of numeric vectors")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call the C++ implementation
    result <- .Call("S_graph_spectral_lowess_mat",
                    adj.list.0based,
                    weight.list,
                    Y,
                    as.integer(n.evectors),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.numeric(min.bw.factor),
                    as.numeric(max.bw.factor),
                    as.numeric(dist.normalization.factor),
                    as.integer(kernel.type),
                    as.numeric(precision),
                    as.integer(n.cleveland.iterations),
                    as.logical(with.errors),
                    as.logical(with.scale),
                    as.logical(verbose))
    
    # Add column names if Y is a matrix with column names
    if (is.matrix(Y) && !is.null(colnames(Y))) {
        colnames(result$predictions) <- colnames(Y)
        if (!is.null(result$errors)) colnames(result$errors) <- colnames(Y)
        if (!is.null(result$scale)) colnames(result$scale) <- colnames(Y)
    }
    
    return(result)
}

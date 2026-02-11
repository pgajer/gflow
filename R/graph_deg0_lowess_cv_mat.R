#' Matrix Version of Graph-Based LOWESS with Cross-Validation
#'
#' @description
#' This function implements a matrix version of graph-based LOWESS (Locally Weighted 
#' Smoothing) of degree 0 with spatially-stratified cross-validation for bandwidth selection.
#' It processes multiple response variables simultaneously, finding the optimal bandwidth
#' for each response variable. This is more efficient than calling the vector version
#' for each response variable separately.
#'
#' @details
#' The algorithm performs the following steps:
#' \enumerate{
#'   \item Creates a maximal packing of vertices to serve as fold seed points
#'   \item Assigns all vertices to the nearest seed point to form spatially coherent folds
#'   \item For each candidate bandwidth, performs cross-validation across the folds for each response variable
#'   \item Selects the bandwidth with the lowest cross-validation error for each response variable
#'   \item Fits the final model with the optimal bandwidth for each response variable
#' }
#'
#' @param adj.list A list of integer vectors. Each vector contains the indices of vertices
#'        adjacent to the vertex at the corresponding list position.
#' @param weight.list A list of numeric vectors. Each vector contains the weights (distances)
#'        of edges to the adjacent vertices specified in adj.list.
#' @param Y A list of numeric vectors. Each vector represents a response variable measured
#'        at each vertex of the graph.
#' @param min.bw.factor Minimum bandwidth as a factor of graph diameter. The actual minimum
#'        bandwidth will be min.bw.factor * graph.diameter.
#' @param max.bw.factor Maximum bandwidth as a factor of graph diameter. The actual maximum
#'        bandwidth will be max.bw.factor * graph.diameter.
#' @param n.bws Number of bandwidths to test in the grid between min.bw and max.bw.
#' @param log.grid Logical, if TRUE, use logarithmic spacing for the bandwidth grid;
#'        if FALSE, use linear spacing.
#' @param kernel.type Integer specifying the kernel function to use.
#'        Possible values: 0=uniform, 1=triangular, 2=epanechnikov, 3=quartic, 
#'        4=triweight, 5=tricube, 6=gaussian, 7=cosine.
#' @param dist.normalization.factor Factor for normalizing distances in kernel weights.
#'        A typical value is 1.1, which ensures all normalized distances fall within
#'        the effective support of most kernel functions.
#' @param use.uniform.weights Whether to use uniform weights instead of kernel weights
#' @param n.folds Number of cross-validation folds.
#' @param with.bw.predictions Logical, if TRUE, compute and return predictions for all
#'        bandwidths; if FALSE, only return predictions for the optimal bandwidths.
#' @param precision Numeric precision parameter for binary search and comparisons.
#' @param verbose Logical, if TRUE, print progress information during computation.
#'
#' @return A list containing:
#' \describe{
#'   \item{predictions}{A list of numeric vectors, where \code{predictions[[j]][i]} is the smoothed
#'         value for the jth response variable at the ith vertex}
#'   \item{bw.predictions}{A nested list structure, where \code{bw.predictions[[j]][[bw.idx]]}
#'         contains predictions for the jth response variable using the bandwidth at index bw.idx.
#'         Only included if with.bw.predictions is TRUE.}
#'   \item{bw.errors}{A list of numeric vectors, where \code{bw.errors[[j]][bw.idx]} is the
#'         cross-validation error for the jth response variable using the bandwidth at index bw.idx}
#'   \item{bws}{Numeric vector of bandwidths used in cross-validation}
#'   \item{opt.bws}{Numeric vector of optimal bandwidths for each response variable}
#'   \item{opt.bw.idxs}{Integer vector of indices of the optimal bandwidths in the bws vector
#'         for each response variable (1-based indices)}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 3 response variables
#' adj.list <- list(c(2,3), c(1,3), c(1,2))
#' weight.list <- list(c(1,1), c(1,1), c(1,1))
#' Y <- list(c(1,2,3), c(4,5,6), c(7,8,9))
#' 
#' # Run the algorithm
#' result <- graph.deg0.lowess.cv.mat(
#'   adj.list, weight.list, Y,
#'   min.bw.factor = 0.1, max.bw.factor = 0.5,
#'   n.bws = 10, log.grid = TRUE, kernel.type = 2,
#'   dist.normalization.factor = 1.1, n.folds = 3,
#'   with.bw.predictions = FALSE, precision = 1e-6, verbose = TRUE
#' )
#' 
#' # Access results
#' result$predictions  # Smoothed values for each response variable
#' result$opt.bws      # Optimal bandwidths for each response variable
#' }
#'
#' @export
graph.deg0.lowess.cv.mat <- function(adj.list,
                                     weight.list,
                                     Y,
                                     min.bw.factor = 0.1,
                                     max.bw.factor = 0.5,
                                     n.bws = 10,
                                     log.grid = TRUE,
                                     kernel.type = 7L,
                                     dist.normalization.factor = 1.1,
                                     use.uniform.weights = FALSE,
                                     n.folds = 5,
                                     with.bw.predictions = FALSE,
                                     precision = 1e-6,
                                     verbose = FALSE) {
  
    ## Input validation
    if (!is.list(adj.list) || !all(sapply(adj.list, is.integer) | sapply(adj.list, is.numeric)))
        stop("adj.list must be a list of integer vectors")

    if (!is.list(weight.list) || !all(sapply(weight.list, is.numeric)))
        stop("weight.list must be a list of numeric vectors")

    if (length(adj.list) != length(weight.list))
        stop("adj.list and weight.list must have the same length")

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

        ## Verify all elements have the same length
        lengths <- sapply(Y, length)
        if (any(lengths != lengths[1])) {
            stop("All vectors in Y list must have the same length")
        }
    } else {
        stop("Y must be a matrix or a list of numeric vectors")
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call("S_graph_deg0_lowess_cv_mat",
                    adj.list.0based,
                    weight.list,
                    Y,
                    as.double(min.bw.factor),
                    as.double(max.bw.factor),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.integer(kernel.type),
                    as.double(dist.normalization.factor),
                    as.logical(use.uniform.weights),
                    as.integer(n.folds),
                    as.logical(with.bw.predictions),
                    as.double(precision),
                    as.logical(verbose))

    ## Properly name components for better usability
    names(result) <- c("predictions", "bw_predictions", "bw_errors", "bws", "opt_bw", "opt_bw_idx")

    ## Name response variables if they have names
    if (!is.null(names(Y))) {
        response.names <- names(Y)
        names(result$predictions) <- response.names
        names(result$bw_errors) <- response.names
        names(result$opt_bws) <- response.names
        names(result$opt_bw_idxs) <- response.names

        if (with.bw.predictions)
            names(result$bw_predictions) <- response.names
    }

    return(result)
}

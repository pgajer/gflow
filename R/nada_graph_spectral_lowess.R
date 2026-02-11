#' Non-adaptive Local Regression on Graphs Using Spectral Embedding
#'
#' @description Performs local regression on graph data using spectral embeddings
#' with adaptive bandwidth selection.
#'
#' @details This function implements a graph-based extension of LOWESS (Locally
#' Weighted Scatterplot Smoothing) that uses spectral embedding to transform
#' graph distances into a Euclidean space suitable for local linear regression.
#' For each vertex, the function:
#' \itemize{
#'   \item Finds all vertices within the maximum bandwidth radius
#'   \item Creates a local spectral embedding using graph Laplacian eigenvectors
#'   \item Fits weighted linear models at multiple candidate bandwidths
#'   \item Selects the optimal bandwidth based on leave-one-out cross-validation error
#'   \item Computes smoothed predictions using the optimal model
#' }
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list\[\[i\]\]} contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of numeric vectors with edge weights corresponding to adjacencies.
#'   Each element \code{weight.list\[\[i\]\]\[j\]} is the weight of the edge from vertex i to
#'   \code{adj.list\[\[i\]\]\[j\]}.
#' @param y A numeric vector of response values for each vertex in the graph.
#' @param n.evectors Integer specifying the number of eigenvectors to use in the spectral
#'   embedding (default: 5).
#' @param n.bws Integer specifying the number of candidate bandwidths to evaluate (default: 10).
#' @param log.grid Logical indicating whether to use logarithmic spacing for bandwidth
#'   grid (default: TRUE).
#' @param min.bw.factor Numeric value specifying the minimum bandwidth as a fraction of
#'   graph diameter (default: 0.05).
#' @param max.bw.factor Numeric value specifying the maximum bandwidth as a fraction of
#'   graph diameter (default: 0.25).
#' @param dist.normalization.factor Numeric factor for normalizing distances when calculating
#'   kernel weights (default: 1.0).
#' @param kernel.type Integer specifying the kernel function for weighting vertices:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#'        Default is 7.
#' @param precision Numeric value specifying the precision tolerance for binary search and
#'   optimization algorithms (default: 0.001).
#' @param n.cleveland.iterations Number of Cleveland's robustness iterations (default: 0)
#' @param verbose Logical indicating whether to display progress information (default: FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{predictions}: Numeric vector of smoothed values for each vertex
#'   \item \code{errors}: Numeric vector of leave-one-out cross-validation errors
#'   \item \code{scale}: Numeric vector of optimal bandwidths (local scales) for each vertex
#'   \item \code{graph.diameter}: Numeric scalar with computed graph diameter
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
#'   neighbors\[neighbors == 0\] <- n
#'   neighbors\[neighbors == n+1\] <- 1
#'
#'   adj.list\[\[i\]\] <- neighbors
#'   weight.list\[\[i\]\] <- rep(1, length(neighbors))
#' }
#'
#' # Generate response values with spatial pattern
#' y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
#'
#' # Apply spectral LOWESS
#' result <- nada.graph.spectral.lowess(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   n.evectors = 5,
#'   verbose = TRUE
#' )
#'
#' # Plot results
#' plot(y, type="l", col="gray", main="Graph Spectral LOWESS")
#' lines(result$predictions, col="red", lwd=2)
#' }
#'
#' @export
nada.graph.spectral.lowess <- function(adj.list,
                                       weight.list,
                                       y,
                                       n.evectors = 5,
                                       n.bws = 20,
                                       log.grid = TRUE,
                                       min.bw.factor = 0.05,
                                       max.bw.factor = 0.33,
                                       dist.normalization.factor = 1.1,
                                       kernel.type = 7L,
                                       precision = 0.001,
                                       n.cleveland.iterations = 0L,
                                       verbose = FALSE) {

    ## Basic parameter validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices in the graph")
    }

    ## for (i in seq_along(adj.list)) {
    ##   if (length(adj.list[[i]]) != length(weight.list[[i]])) {
    ##     stop(paste0("Adjacency and weight lists have different lengths at vertex ", i))
    ##   }
    ## }

    if (n.evectors < 2) {
        warning("n.evectors should be at least 2; setting to 2")
        n.evectors <- 2
    }

    if (n.bws < 1)
        stop("n.bws must be positive")

    if (!is.numeric(min.bw.factor))
        stop("min.bw.factor must be numeri")

    if (!is.numeric(max.bw.factor) || max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (dist.normalization.factor < 1.05)
        stop("dist.normalization.factor must be greater than or equal to 1.05")

    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel.type' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    if (precision <= 0)
        stop("precision must be positive")

    if(!n.cleveland.iterations %in% c(0L, 1L, 2L, 3L, 4L, 5L)) {
        stop("'n.cleveland.iterations' must be and integer between 0 and 5")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ implementation
    result <- .Call("S_nada_graph_spectral_lowess",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.integer(n.evectors),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.numeric(min.bw.factor),
                    as.numeric(max.bw.factor),
                    as.numeric(dist.normalization.factor),
                    as.integer(kernel.type),
                    as.numeric(precision),
                    as.integer(n.cleveland.iterations),
                    as.logical(verbose))

    return(result)
}

#' Degree 0 LOWESS (Locally Weighted Average) on Graphs
#'
#' @description Performs local constant fitting (degree 0 LOWESS) on graph data
#' using a fixed bandwidth.
#'
#' @details This function implements a simplified version of LOWESS for graph data
#' where only degree 0 local models (weighted averages) are fit. For each vertex,
#' the function:
#' \enumerate{
#'   \item Finds all vertices within the specified bandwidth radius
#'   \item Computes a weighted average of response values using kernel weights
#'   \item Returns the smoothed prediction
#' }
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of numeric vectors with edge weights corresponding to adjacencies.
#'   Each element \code{weight.list[[i]][j]} is the weight of the edge from vertex i to
#'   \code{adj.list[[i]][j]}.
#' @param y A numeric vector of response values for each vertex in the graph.
#' @param bandwidth A numeric value specifying the fixed bandwidth (radius) to use for
#'   all local neighborhoods.
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
#'   kernel weights (default: 1.0).
#' @param verbose Logical indicating whether to display progress information (default: FALSE).
#'
#' @return A numeric vector of smoothed predictions for each vertex.
#'
#' @examples
#' \dontrun{
#' # Create a simple ring graph
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
#' # Generate response values with spatial pattern plus noise
#' y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
#'
#' # Apply degree 0 LOWESS with different bandwidths
#' bw1 <- 5
#' bw2 <- 10
#' bw3 <- 20
#'
#' result1 <- graph.deg0.lowess(adj.list, weight.list, y, bw1)
#' result2 <- graph.deg0.lowess(adj.list, weight.list, y, bw2)
#' result3 <- graph.deg0.lowess(adj.list, weight.list, y, bw3)
#'
#' # Plot results
#' plot(y, type="p", col="gray", main="Graph Degree 0 LOWESS")
#' lines(result1, col="blue", lwd=2)
#' lines(result2, col="red", lwd=2)
#' lines(result3, col="green", lwd=2)
#' legend("topright", legend=c(paste0("BW=", bw1), 
#'                             paste0("BW=", bw2),
#'                             paste0("BW=", bw3)),
#'        col=c("blue", "red", "green"), lwd=2)
#' }
#'
#' @export
graph.deg0.lowess <- function(adj.list,
                              weight.list,
                              y,
                              bandwidth,
                              kernel.type = 7L,
                              dist.normalization.factor = 1.1,
                              verbose = FALSE) {

    ## Basic parameter validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices in the graph")
    }

    if (!is.numeric(bandwidth) || bandwidth <= 0) {
        stop("bandwidth must be a positive number")
    }

    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 0) {
        stop("dist.normalization.factor must be a positive number")
    }

    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel.type' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    # Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ implementation
    result <- .Call("S_graph_deg0_lowess",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.numeric(bandwidth),
                    as.integer(kernel.type),
                    as.numeric(dist.normalization.factor),
                    as.logical(verbose))

    return(result)
}

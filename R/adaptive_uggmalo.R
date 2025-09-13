#' Adaptive Uniform Grid Graph Model-Averaged Local Linear Regression
#'
#' @description
#' Performs model-averaged local linear regression on graph-structured data using an adaptive
#' uniform grid approach. This method combines graph structure preservation with local linear
#' modeling to capture both global and local patterns in the data.
#'
#' @details
#' The adaptive UGGMALO algorithm proceeds through several steps:
#' 1. Creates a uniform grid representation of the input graph
#' 2. Computes optimal bandwidths using cross-validation
#' 3. Fits local linear models along paths through the graph
#' 4. Combines predictions using model averaging
#'
#' The function supports optional bootstrap confidence intervals and permutation testing
#' for statistical inference.
#'
#' @param adj.list List of integer vectors. Each vector contains indices of vertices
#'        adjacent to the corresponding vertex. Indices should be 1-based.
#' @param weight.list List of numeric vectors. Each vector contains weights of edges
#'        corresponding to adjacencies in adj.list.
#' @param y Numeric vector of response values at each vertex.
#' @param min.path.size Integer. Minimum number of vertices required in valid paths.
#' @param n.grid.vertices Integer. Number of vertices in uniform grid representation.
#' @param n.bws Integer. Number of candidate bandwidths to evaluate.
#' @param min.bw.factor Numeric. Factor multiplied by graph diameter for minimum bandwidth.
#' @param max.bw.factor Numeric. Factor multiplied by graph diameter for maximum bandwidth.
#' @param max.iterations Integer. Maximum number of iterations for convergence.
#' @param precision Numeric. Precision threshold for numerical computations.
#' @param dist.normalization.factor Numeric. Factor for normalizing graph distances.
#' @param kernel.type Integer. Type of kernel function (1: Gaussian, 2: Triangular).
#' @param tolerance Numeric. Convergence tolerance for model fitting.
#' @param n.bb Integer. Number of bootstrap iterations (0 for no bootstrap).
#' @param cri.probability Numeric. Confidence level for bootstrap intervals (0-1).
#' @param n.perms Integer. Number of permutation test iterations (0 for no testing).
#' @param blending.coef Numeric. Blending coefficient for model averaging.
#' @param verbose Logical. Whether to print progress information.
#'
#' @return A list containing:
#' \describe{
#'   \item{graph.diameter}{Numeric. Computed diameter of input graph.}
#'   \item{grid.opt.bw}{Numeric vector. Optimal bandwidth for each grid vertex.}
#'   \item{predictions}{Numeric vector. Model-averaged predictions for original vertices.}
#'   \item{grid.predictions}{Numeric vector. Model-averaged predictions for grid vertices.}
#'   \item{bb.predictions}{Matrix. Bootstrap predictions (if n.bb > 0).}
#'   \item{cri.lower}{Numeric vector. Lower confidence bounds (if n.bb > 0).}
#'   \item{cri.upper}{Numeric vector. Upper confidence bounds (if n.bb > 0).}
#'   \item{null.predictions}{Matrix. Permutation test predictions (if n.perms > 0).}
#'   \item{p.values}{Numeric vector. Vertex-wise p-values (if n.perms > 0).}
#'   \item{effect.sizes}{Numeric vector. Effect sizes (if n.perms > 0).}
#'   \item{significant.vertices}{Logical vector. Significance indicators (if n.perms > 0).}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 3 vertices
#' adj.list <- list(c(2), c(1, 3), c(2))
#' weight.list <- list(c(1), c(1, 1), c(1))
#' y <- c(1, 2, 3)
#'
#' # Run basic analysis
#' result <- adaptive.uggmalo(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   min.path.size = 2,
#'   n.grid.vertices = 5,
#'   n.bws = 10
#' )
#'
#' # Run with bootstrap confidence intervals
#' result.boot <- adaptive.uggmalo(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   min.path.size = 2,
#'   n.grid.vertices = 5,
#'   n.bws = 10,
#'   n.bb = 100
#' )
#' }
#'
#' @export
adaptive.uggmalo <- function(adj.list,
                             weight.list,
                             y,
                             min.path.size,
                             n.grid.vertices,
                             n.bws,
                             min.bw.factor = 0.05,
                             max.bw.factor = 0.6,
                             max.iterations = 20,
                             precision = 0.0001,
                             dist.normalization.factor = 1.1,
                             kernel.type = 7L,
                             tolerance = 1e-6,
                             n.bb = 0L,
                             cri.probability = 0.95,
                             n.perms = 0L,
                             blending.coef = 0.1,
                             verbose = FALSE) {

    ## Input validation
    if (!is.list(adj.list) || !is.list(weight.list))
        stop("adj.list and weight.list must be lists")

    if (length(adj.list) != length(weight.list))
        stop("adj.list and weight.list must have the same length")

    if (length(y) != length(adj.list))
        stop("y must have the same length as adj.list")

    if (!is.numeric(y))
        stop("y must be numeric")

    if (min.path.size < 2)
        stop("min.path.size must be at least 2")

    if (n.grid.vertices < 1) ## <<--- this needs to be changed - the grid needs to cover the graph
        stop("n.grid.vertices must be positive")

    if (!is.numeric(max.iterations) || length(max.iterations) != 1 ||
        max.iterations != floor(max.iterations) ||
        max.iterations < 1) {
        stop("max.iterations must be an integer greater than 0")
    }

    if (n.bws < 1)
        stop("n.bws must be positive")

    if (!is.numeric(min.bw.factor))
        stop("min.bw.factor must be numeri")

    if (max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (precision <= 0 || precision > 0.1)
        stop("precision must be positive not greater than 0.1")

    if (dist.normalization.factor < 1.1)
        stop("dist.normalization.factor must be greater than or equal to 1.1")

    if (!is.numeric(kernel.type) || kernel.type < 0 || kernel.type > 7)
        stop("kernel.type must be an integer between 0 and 7")

    if (tolerance <= 0)
        stop("tolerance must be positive")

    if (n.bb < 0)
        stop("n.bb must be non-negative")

    if (cri.probability <= 0 || cri.probability >= 1)
        stop("cri.probability must be between 0 and 1")

    if (n.perms < 0)
        stop("n.perms must be non-negative")

    if (blending.coef < 0 || blending.coef > 1) {
        stop("blending.coef must be between 0 and 1")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call(S_adaptive_uggmalo,
                   adj.list.0based,
                   weight.list,
                   as.double(y),
                   as.integer(min.path.size),
                   as.integer(n.grid.vertices),
                   as.integer(n.bws),
                   as.double(min.bw.factor),
                   as.double(max.bw.factor),
                   as.integer(max.iterations),
                   as.double(precision),
                   as.double(dist.normalization.factor),
                   as.integer(kernel.type),
                   as.double(tolerance),
                   as.integer(n.bb),
                   as.double(cri.probability),
                   as.integer(n.perms),
                   as.double(blending.coef),
                   as.logical(verbose))

    result
}

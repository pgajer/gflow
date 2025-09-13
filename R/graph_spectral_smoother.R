#' Graph Spectral Smoother utilizing spectrum of graph Laplacian
#'
#' @description
#' Smooths a graph signal using a spectral smoothing technique based on the graph Laplacian.
#' It determines the optimal number of low-eigenvalue eigenvectors for representing
#' the graph function and returns the smoothed signal along with cross-validation errors.
#'
#' @param graph A list of integer vectors representing the graph's adjacency list.
#'              Each vector contains the indices of neighboring vertices (1-based).
#' @param edge.lengths A list of edge lengths. Structure should match that of graph.
#' @param y A numeric vector representing the graph signal to be smoothed.
#' @param weights A numeric vector of weights for each vertex. If NULL, uniform weights are used.
#' @param imputation.method Integer specifying the imputation method:
#'        0: LOCAL_MEAN_THRESHOLD (default)
#'        1: NEIGHBORHOOD_MATCHING
#'        2: ITERATIVE_NEIGHBORHOOD_MATCHING
#'        3: SUPPLIED_THRESHOLD
#'        4: GLOBAL_MEAN_THRESHOLD
#' @param max.iterations Integer. Maximum iterations for iterative matching method. Default is 10.
#' @param convergence.threshold Numeric. Convergence threshold for iterative matching. Default is 1e-6.
#' @param apply.binary.threshold Logical. Whether to apply binary thresholding. Default is TRUE.
#' @param binary.threshold Numeric. Threshold for binary classification (0-1). Default is 0.5.
#' @param kernel Integer specifying the kernel function:
#'        0: Constant, 1: Epanechnikov (default), 2: Triangular,
#'        3: Truncated Exponential, 4: Normal
#' @param dist.normalization.factor Numeric. Scaling factor for distance normalization. Default is 1.01.
#' @param n.CVs Integer. Number of cross-validation iterations. Default is 0.
#' @param n.CV.folds Integer. Number of cross-validation folds. Default is 10.
#' @param epsilon Numeric. Small positive constant for numerical stability. Default is 1e-10.
#' @param min.plambda Numeric. Lower bound on proportion of eigenvectors to use. Default is 0.01.
#' @param max.plambda Numeric. Upper bound on proportion of eigenvectors to use. Default is 0.20.
#' @param seed Integer. Seed for random number generation. Default is 0.
#' @param verbose Logical. Whether to print additional information. Default is FALSE.
#'
#' @return A list containing:
#'   \item{optimal.num.eigenvectors}{Optimal number of eigenvectors used for smoothing}
#'   \item{y.smoothed}{Smoothed signal}
#'   \item{cv.errors}{Cross-validation errors for each fold}
#'   \item{mean.cv.errors}{Mean cross-validation errors across all folds}
#'   \item{median.cv.errors}{Median cross-validation errors across all folds}
#'   \item{Cmean.cv.errors}{Mean cross-validation errors computed in C++}
#'   \item{min.plambda, max.plambda}{Input lambda values}
#'   \item{min.num.eigenvectors, max.num.eigenvectors}{Range of eigenvectors used}
#'   \item{evalues, evectors}{Eigenvalues and eigenvectors of the graph Laplacian}
#'   \item{low.pass.ys}{Low-pass filtered versions of y for different numbers of eigenvectors}
#'
#' @export
graph.spectral.smoother <- function(graph,
                                    edge.lengths,
                                    y,
                                    weights = NULL,
                                    imputation.method = 0,
                                    max.iterations = 10,
                                    convergence.threshold = 1e-6,
                                    apply.binary.threshold = TRUE,
                                    binary.threshold = 0.5,
                                    kernel = 1,
                                    dist.normalization.factor = 1.1,
                                    n.CVs = 0,
                                    n.CV.folds = 10,
                                    epsilon = 1e-10,
                                    min.plambda = 0.01,
                                    max.plambda = 0.20,
                                    seed = 0,
                                    verbose = FALSE) {

    # Input validation
    if (!is.list(graph) || !all(sapply(graph, is.numeric)))
        stop("graph must be a list of numeric vectors.")

    if (!is.list(edge.lengths) || !all(sapply(edge.lengths, is.numeric)))
        stop("edge.lengths must be NULL or a list of numeric vectors.")

    if (!is.numeric(y))
        stop("y must be a numeric vector.")

    if (is.null(weights)) {
        weights <- rep(1.0, length(y))
    } else if (!is.null(weights) && (!is.numeric(weights) || length(weights) != length(y)))
        stop("weights must be NULL or a numeric vector of the same length as y.")

    if (!is.numeric(imputation.method) || length(imputation.method) != 1 || !(imputation.method %in% 0:4))
        stop("imputation.method must be an integer between 0 and 4.")

    if (!is.numeric(max.iterations) || length(max.iterations) != 1 || max.iterations < 1)
        stop("max.iterations must be a positive integer.")

    if (!is.numeric(convergence.threshold) || length(convergence.threshold) != 1 || convergence.threshold <= 0)
        stop("convergence.threshold must be a positive numeric value.")

    if (!is.logical(apply.binary.threshold))
        stop("apply.binary.threshold must be a logical value.")

    if (!is.numeric(binary.threshold) || length(binary.threshold) != 1 || binary.threshold < 0 || binary.threshold > 1)
        stop("binary.threshold must be a numeric value between 0 and 1.")

    if (!is.numeric(kernel) || length(kernel) != 1 || !(kernel %in% 0:4))
        stop("kernel must be an integer between 0 and 4.")

    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1 || dist.normalization.factor < 1)
        stop("dist.normalization.factor must be a numeric value greater than or equal to 1.")

    if (!is.numeric(n.CVs) || length(n.CVs) != 1 || n.CVs < 0)
        stop("n.CVs must be a non-negative integer.")

    if (!is.numeric(n.CV.folds) || length(n.CV.folds) != 1 || n.CV.folds < 2)
        stop("n.CV.folds must be an integer greater than 1.")

    if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0)
        stop("epsilon must be a positive numeric value.")

    if (!is.numeric(min.plambda) || length(min.plambda) != 1 || min.plambda <= 0 || min.plambda >= 1)
        stop("min.plambda must be a numeric value between 0 and 1.")

    if (!is.numeric(max.plambda) || length(max.plambda) != 1 || max.plambda <= 0 || max.plambda > 1 || max.plambda <= min.plambda)
        stop("max.plambda must be a numeric value between min.plambda and 1.")

    if (!is.numeric(seed) || length(seed) != 1 || seed != round(seed) || seed < 0)
        stop("seed must be a non-negative integer.")

    if (!is.logical(verbose))
        stop("verbose must be a logical value.")

    n.vertices <- length(y)

    if (n.vertices != length(graph))
        stop("The lengths of graph and y must be the same.")
    if (n.vertices != length(edge.lengths))
        stop("The lengths of edge.lengths and y must be the same.")

    if (!all(sapply(seq_along(graph), function(i) length(graph[[i]]) == length(edge.lengths[[i]])))) {
        stop("The structure of graph and edge.lengths do not match.")
    }

    if (verbose) {
        min.num.eigenvectors <- as.integer(min.plambda * n.vertices)
        max.num.eigenvectors <- as.integer(max.plambda * n.vertices)
        cat("min.plambda:", min.plambda, "\tmin.num.eigenvectors:", min.num.eigenvectors, "\n")
        cat("max.plambda:", max.plambda, "\tmax.num.eigenvectors:", max.num.eigenvectors, "\n")
    }

    # Converting each component of adj.list to an integer vector (0-based indexing)
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    result <- .Call(S_graph_spectral_smoother,
                    graph.0based,
                    edge.lengths,
                    weights,
                    y,
                    as.integer(imputation.method),
                    as.integer(max.iterations),
                    as.numeric(convergence.threshold),
                    as.logical(apply.binary.threshold),
                    as.numeric(binary.threshold),
                    as.integer(kernel),
                    as.numeric(dist.normalization.factor),
                    as.integer(n.CVs),
                    as.integer(n.CV.folds),
                    as.numeric(epsilon),
                    as.numeric(min.plambda),
                    as.numeric(max.plambda),
                    as.integer(seed))

    cv.errors <- result$cv_errors
    mean.cv.errors <- apply(cv.errors, 1, mean, na.rm = TRUE)
    median.cv.errors <- apply(cv.errors, 1, median, na.rm = TRUE)

    list(
        optimal.num.eigenvectors = result$optimal_num_eigenvectors,
        y.smoothed = result$y_smoothed,
        cv.errors = cv.errors,
        mean.cv.errors = mean.cv.errors,
        median.cv.errors = median.cv.errors,
        Cmean.cv.errors = result$mean_cv_errors,
        min.plambda = min.plambda,
        max.plambda = max.plambda,
        min.num.eigenvectors = result$min_num_eigenvectors,
        max.num.eigenvectors = result$max_num_eigenvectors,
        evalues = result$evalues,
        evectors = result$evectors,
        low.pass.ys = result$low_pass_ys
    )
}

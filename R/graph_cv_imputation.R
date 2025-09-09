#' Cross-Validation Imputation on Graphs
#'
#' @description
#' Performs cross-validation imputation on graph-structured data using various methods.
#' This function supports both binary and continuous data, and offers different
#' imputation strategies including local mean thresholding, neighborhood matching,
#' and iterative approaches.
#'
#' @param test.set A numeric vector of indices for the test set vertices (1-based).
#' @param graph A list of numeric vectors representing the graph structure.
#'   Each element of the list corresponds to a vertex and contains the indices of its neighbors.
#' @param edge.lengths A list of edge lengths. The structure should match that of `graph`.
#' @param y A numeric vector of original vertex values.
#' @param y.binary A logical value indicating whether the data is binary (TRUE) or continuous (FALSE).
#' @param imputation.method A string or integer specifying the imputation method. Options are:
#'        - 0 or "local_mean_threshold": Uses the mean of y computed over the training vertices (default).
#'        - 1 or "neighborhood_matching": Uses a matching method based on local neighborhood statistics.
#'        - 2 or "iterative_neighborhood_matching": Uses an iterative version of the neighborhood matching method.
#'        - 3 or "supplied_threshold": Uses a user-supplied threshold value.
#'        - 4 or "global_mean_threshold": Uses the global mean of y across all vertices.
#' @param max.iterations An integer specifying the maximum number of iterations for iterative methods.
#' @param convergence.threshold A numeric value specifying the convergence threshold for iterative methods.
#' @param apply.binary.threshold A logical value indicating whether to apply binary thresholding to the results.
#' @param binary.threshold A numeric value between 0 and 1 specifying the threshold for binary classification.
#' @param kernel A character string specifying the kernel function for distance weighting.
#'   Options are "Box", "Triangular", "Epanechnikov", or "Gaussian".
#' @param dist.normalization.factor A numeric value greater than 1 used for normalizing distances.
#'
#' @return A numeric vector containing the imputed values for the test set vertices.
#'
#' @details
#' This function serves as an R interface to a C++ implementation of graph-based
#' imputation methods. It handles various input checks and data type conversions
#' before calling the underlying C++ function.
#'
#' The choice of imputation method, kernel function, and other parameters allows
#' for flexibility in addressing different types of graph-structured data and
#' imputation scenarios.
#'
#' @note
#' - The function assumes that the graph structure is consistent and that all
#'   vertex indices in test.set and training.set are valid.
#' - For binary data (y.binary = TRUE), the imputed values will be either 0 or 1.
#' - For continuous data, the imputed values may fall outside the range of the
#'   original y values, depending on the chosen method.
#'
#' @examples
#' \dontrun{
#' # Example with a small graph
#' graph <- list(c(2,3), c(1,3), c(1,2,4), c(3))
#' y <- c(1, 0, 1, 0)
#' test.set <- c(2, 4)
#' result <- cv.imputation(test.set, graph, y = y, y.binary = TRUE,
#'                         imputation.method = "LOCAL_MEAN_THRESHOLD")
#' print(result)
#' }
#'
#' @seealso
#' For more details on graph-based imputation methods, see the documentation
#' of the underlying C++ implementation.
#'
#' @export
cv.imputation <- function(test.set,
                          graph,
                          edge.lengths,
                          y,
                          y.binary,
                          imputation.method = 1,
                          max.iterations = 10,
                          convergence.threshold = 1e-6,
                          apply.binary.threshold = TRUE,
                          binary.threshold = 0.5,
                          kernel = "Epanechnikov",
                          dist.normalization.factor = 1.01) {

    ## Checking input types and values
    if (!is.numeric(test.set)) {
        stop("test.set must be a numeric vector.")
    }
    if (!all(test.set == floor(test.set))) {
        stop("test.set must contain integer values")
    }
    if (any(test.set < 1)) {
        stop("test.set must contain positive integers")
    }
    if (any(duplicated(test.set))) {
        stop("test.set must not have overlapping elements")
    }

    if (!is.list(graph)) stop("graph must be a list")
    if (!all(sapply(graph, is.numeric))) stop("all elements of graph must be numeric")

    if (!is.list(edge.lengths)) stop("edge.lengths must be a list")
    if (length(edge.lengths) > 0 && !all(sapply(edge.lengths, is.numeric))) stop("all elements of edge.lengths must be numeric")

    if (!is.numeric(y)) stop("y must be a numeric vector")
    if (length(y) != length(graph)) stop("length of y must match the number of vertices in graph")

    if (!is.logical(y.binary)) stop("y.binary must be a logical value")

    if (is.character(imputation.method)) {
        imputation.methods <- c("local_mean_threshold", "neighborhood_matching", "iterative_neighborhood_matching", "supplied_threshold", "global_mean_threshold")
        imputation.method <- match.arg(imputation.method, imputation.methods)
        imputation.method <- which(imputation.methods == imputation.method) - 1 ## C++ enum is 0-indexed
    } else if (!is.numeric(imputation.method) || length(imputation.method) != 1 || !(imputation.method %in% 0:4)) {
        stop("imputation.method must be an integer between 0 and 4 or one of the following strings: 'local_mean_threshold', 'neighborhood_matching', 'iterative_neighborhood_matching', 'supplied_threshold', 'global_mean_threshold'")
    }

    if (!is.numeric(max.iterations) || max.iterations < 1 || max.iterations != floor(max.iterations)) {
        stop("max.iterations must be a positive integer")
    }

    if (!is.numeric(convergence.threshold) || convergence.threshold <= 0) {
        stop("convergence.threshold must be a positive numeric value")
    }

    if (!is.logical(apply.binary.threshold)) stop("apply.binary.threshold must be a logical value")

    if (!is.numeric(binary.threshold) || binary.threshold < 0 || binary.threshold > 1) {
        stop("binary.threshold must be a numeric value between 0 and 1")
    }

    if (is.character(kernel)) {
        kernel <- match.arg(kernel, c("epanechnikov", "triangular", "truncated_exponential", "normal"))
        kernel <- switch(kernel,
                         epanechnikov = 1,
                         triangular = 2,
                         truncated_exponential = 3,
                         normal = 4)
    } else if (kernel != as.integer(kernel) || kernel < 1 || kernel > 4) {
        stop("kernel must be an integer between 1 and 4 or one of the following strings: 'epanechnikov', 'triangular', 'truncated_exponential', 'normal'")
    }

    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 1) {
        stop("dist.normalization.factor must be a numeric value greater than 1")
    }

    ## Convert graph indices to 0-based
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    result <- .Call("S_cv_imputation",
                    as.integer(test.set),
                    graph.0based,
                    edge.lengths,
                    as.double(y),
                    as.logical(y.binary),
                    as.integer(imputation.method),
                    as.integer(max.iterations),
                    as.double(convergence.threshold),
                    as.logical(apply.binary.threshold),
                    as.double(binary.threshold),
                    as.integer(kernel),
                    as.double(dist.normalization.factor))

    return(result)
}

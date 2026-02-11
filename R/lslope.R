## Local Slope (Asymmetric Association) Measures
##
## Functions for computing asymmetric local association measures between a
## directing function y and response function z defined on graph vertices.
## Unlike symmetric local correlation, these measures treat y as the directing
## variable and z as the response.
##
## This file provides a unified interface supporting vector-vector, vector-matrix,
## and matrix-matrix computations with automatic dispatch based on input types.

################################################################################
#
# Main Dispatcher: lslope.grad()
#
################################################################################

#' Gradient-Restricted Local Slope
#'
#' @description
#' Computes the local slope of z with respect to y along the gradient direction
#' of y at each vertex. This asymmetric measure captures "for each unit increase
#' in y along its steepest direction, how much does z change?"
#'
#' The function automatically dispatches to specialized implementations based on
#' input types: vector-vector, vector-matrix, or matrix-matrix.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights.
#'   Must have same structure as adj.list.
#' @param y Numeric vector or matrix of directing function values. For vectors,
#'   length must equal the number of vertices. For matrices, number of rows must
#'   equal the number of vertices.
#' @param z Numeric vector or matrix of response function values. Dimension
#'   requirements match those for y.
#' @param type Character scalar specifying the type of measure:
#'   \describe{
#'     \item{"slope"}{Raw ratio Delta z / Delta y (unbounded). Default.}
#'     \item{"normalized"}{Sigmoid-normalized ratio (bounded to \eqn{[-1,1]}).}
#'     \item{"sign"}{Sign of Delta z along gradient direction \eqn{{-1, 0, +1}}.}
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y:
#'   \describe{
#'     \item{"logratio"}{Log-ratios (compositional data). Default.}
#'     \item{"difference"}{Standard differences (continuous data).}
#'   }
#' @param z.diff.type Character scalar specifying edge difference type for z.
#' @param epsilon Numeric scalar for pseudocount in log-ratios.
#'   If 0 (default), computed adaptively.
#' @param sigmoid.alpha Numeric scalar for sigmoid normalization scale.
#'   If 0 (default), auto-calibrated from median absolute slope.
#' @param sigmoid.type Character scalar specifying sigmoid function type:
#'   \itemize{
#'     \item "tanh": Hyperbolic tangent (default)
#'     \item "arctan": Scaled arctangent (heavier tails)
#'     \item "algebraic": Algebraic sigmoid (simpler derivative)
#'   }
#' @param ascending Logical. If TRUE (default), use ascending gradient
#'   (gradient toward higher y values). If FALSE, use descending gradient.
#' @param instrumented Logical. If TRUE, return full diagnostic output.
#'   If FALSE (default), return only coefficients. Only applies to vector-vector case.
#' @param mc.cores Integer scalar specifying the number of cores for parallel
#'   computation in the matrix-matrix case. Default is 1 (sequential).
#'   Note that parallel::mclapply is not available on Windows.
#'
#' @return The return type depends on the input types:
#'
#'   \strong{Vector-vector (y and z both vectors):}
#'   If instrumented = FALSE, a numeric vector of local slope coefficients.
#'   If instrumented = TRUE, a list with class "lslope_gradient_result" containing
#'   diagnostic information (gradient neighbors, delta values, extremum flags).
#'
#'   \strong{Vector-matrix (y vector, z matrix):}
#'   An object of class "lslope_vector_matrix_result". A numeric matrix of
#'   dimension (n.vertices x n.columns) where column j contains local slope
#'   coefficients for \code{z[,j]} with respect to y.
#'
#'   \strong{Matrix-matrix (y and z both matrices):} An object of class
#'   "lslope_matrix_matrix_result". A 3D numeric array of dimension (n.vertices
#'   x ncol(y) x ncol(z)). Element \code{[v, i, j]} contains the local slope of
#'   \code{z[,j]} with respect to \code{y[,i]} at vertex v.
#'
#' @details
#' The gradient edge at vertex v is the edge from v to the neighbor u that
#' maximizes the difference Delta y = y(u) - y(v) (for ascending=TRUE) or
#' minimizes it (for ascending=FALSE). At local extrema of y, no gradient edge
#' exists and the coefficient is set to 0.
#'
#' @section Asymmetry:
#' Unlike symmetric local correlation lcor(y,z), the local slope is asymmetric:
#' lslope(z; y) is not equal to lslope(y; z) in general. The slope measures the
#' response of z to changes in y, treating y as the independent (directing)
#' variable.
#'
#' @section Dispatch Behavior:
#' \strong{Vector-Vector:} Computes local slope of a single z on a single y.
#'
#' \strong{Vector-Matrix:} Efficiently computes slopes of multiple z columns
#' with respect to a single directing function y. The gradient edges are
#' determined once from y and reused for all z columns.
#'
#' \strong{Matrix-Matrix:} Computes slopes for all (y_i, z_j) pairs. Since the
#' measure is asymmetric, all p x q pairs are computed (no symmetry optimization).
#'
#' @examples
#' \dontrun{
#' ## Simple example
#' adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
#' weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
#'
#' ## Case (a): Vector-vector
#' y <- c(0.1, 0.5, 0.9)  # sPTB prevalence
#' z <- c(0.8, 0.4, 0.1)  # Secondary phylotype abundance
#' result <- lslope(adj.list, weight.list, y, z, type = "normalized")
#'
#' ## Case (b): Vector-matrix (feature screening)
#' Z <- matrix(runif(300), nrow = 100, ncol = 50)
#' colnames(Z) <- paste0("ASV", 1:50)
#' result <- lslope(adj.list, weight.list, y, Z,
#'                  type = "normalized",
#'                  z.diff.type = "logratio")
#'
#' ## Case (c): Matrix-matrix (multiple directing functions)
#' Y <- matrix(runif(200), nrow = 100, ncol = 2)
#' colnames(Y) <- c("sPTB", "PTB")
#' result <- lslope(adj.list, weight.list, Y, Z,
#'                  type = "normalized",
#'                  mc.cores = 4)
#' }
#'
#' @seealso \code{\link{lslope.neighborhood}} for neighborhood-based regression,
#'   \code{\link{lcor}} for symmetric local correlation
#'
#' @export
lslope <- function(adj.list,
                   weight.list,
                   y,
                   z,
                   type = c("slope", "normalized", "sign"),
                   y.diff.type = c("logratio", "difference"),
                   z.diff.type = c("logratio", "difference"),
                   epsilon = 0,
                   sigmoid.alpha = 0,
                   sigmoid.type = c("tanh", "arctan", "algebraic"),
                   ascending = TRUE,
                   instrumented = FALSE,
                   mc.cores = 1L) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)
    sigmoid.type <- match.arg(sigmoid.type)

    ## Detect input types FIRST (before any conversion)
    y.is.matrix <- is.matrix(y) || is.data.frame(y)
    z.is.matrix <- is.matrix(z) || is.data.frame(z)

    ## Convert data frames to matrices
    if (is.data.frame(y)) y <- as.matrix(y)
    if (is.data.frame(z)) z <- as.matrix(z)

    ## Validate common inputs
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")
    if (!is.numeric(sigmoid.alpha) || length(sigmoid.alpha) != 1)
        stop("sigmoid.alpha must be a single numeric value")
    if (!is.logical(ascending) || length(ascending) != 1)
        stop("ascending must be a single logical value")

    n.vertices <- length(adj.list)
    if (length(weight.list) != n.vertices)
        stop("Length of weight.list must equal number of vertices")

    ## Validate y dimensions
    if (y.is.matrix) {
        if (!is.numeric(y)) stop("y must be numeric")
        if (nrow(y) != n.vertices)
            stop("nrow(y) must equal number of vertices")
    } else {
        if (!is.numeric(y)) stop("y must be a numeric vector")
        if (length(y) != n.vertices)
            stop("length(y) must equal number of vertices")
    }

    ## Validate z dimensions
    if (z.is.matrix) {
        if (!is.numeric(z)) stop("z must be numeric")
        if (nrow(z) != n.vertices)
            stop("nrow(z) must equal number of vertices")
    } else {
        if (!is.numeric(z)) stop("z must be a numeric vector")
        if (length(z) != n.vertices)
            stop("length(z) must equal number of vertices")
    }

    ## Dispatch based on input types
    if (!y.is.matrix && !z.is.matrix) {
        ## Case (a): vector-vector
        result <- lslope.vector.vector(
            adj.list,
            weight.list,
            y,
            z,
            type,
            y.diff.type,
            z.diff.type,
            epsilon,
            sigmoid.alpha,
            sigmoid.type,
            ascending,
            instrumented
        )
        return(result)

    } else if (!y.is.matrix && z.is.matrix) {
        ## Case (b): vector-matrix
        result <- lslope.vector.matrix(
            adj.list,
            weight.list,
            y,
            z,
            type,
            y.diff.type,
            z.diff.type,
            epsilon,
            sigmoid.alpha,
            ascending
        )
        return(result)

    } else if (y.is.matrix && z.is.matrix) {
        ## Case (c): matrix-matrix
        result <- lslope.matrix.matrix(
            adj.list,
            weight.list,
            y,
            z,
            type,
            y.diff.type,
            z.diff.type,
            epsilon,
            sigmoid.alpha,
            ascending,
            mc.cores
        )
        return(result)

    } else {
        ## y is matrix, z is vector: not supported
        stop("y is matrix but z is vector: the directing function y should be ",
             "a vector when z is a vector. If you need to compute slopes for ",
             "multiple directing functions against a single response, use ",
             "lslope(adj.list, weight.list, Y, cbind(z), ...)")
    }
}


################################################################################
#
# Case (a): Vector-Vector Implementation
#
################################################################################

#' @describeIn lslope Vector-vector implementation (internal)
#' @keywords internal
lslope.vector.vector <- function(adj.list,
                                      weight.list,
                                      y,
                                      z,
                                      type,
                                      y.diff.type,
                                      z.diff.type,
                                      epsilon,
                                      sigmoid.alpha,
                                      sigmoid.type,
                                      ascending,
                                      instrumented) {

    n.vertices <- length(adj.list)

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    if (instrumented) {
        result <- .Call(
            "S_lslope_gradient_instrumented",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(sigmoid.alpha),
            as.character(sigmoid.type),
            as.logical(ascending),
            PACKAGE = "gflow"
        )

        ## Convert gradient neighbors back to 1-based indexing
        ## C++ returns SIZE_MAX for invalid vertices (local extrema)
        ## We convert these to NA in R
        max.size.t <- .Machine$integer.max * 2 + 1
        result$gradient.neighbors <- ifelse(
            result$gradient.neighbors >= max.size.t - 1,
            NA_integer_,
            as.integer(result$gradient.neighbors + 1)
        )

        ## Add class for method dispatch
        class(result) <- c("lslope_gradient_result", "list")

        ## Add metadata
        attr(result, "type") <- type
        attr(result, "y.diff.type") <- y.diff.type
        attr(result, "z.diff.type") <- z.diff.type
        attr(result, "ascending") <- ascending
        attr(result, "n.vertices") <- n.vertices

        return(result)

    } else {
        result <- .Call(
            "S_lslope_gradient",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(z),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(sigmoid.alpha),
            as.character(sigmoid.type),
            as.logical(ascending),
            PACKAGE = "gflow"
        )

        return(result)
    }
}


################################################################################
#
# Case (b): Vector-Matrix Implementation
#
################################################################################

#' Local Slope Between a Vector and Matrix Columns (C++ Implementation)
#'
#' Compute vertex-level slope coefficients measuring the response of each column
#' of Z to a directing function y on a graph. This implementation uses optimized
#' C++ code with OpenMP parallelization.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#' @param weight.list List of numeric vectors containing edge weights.
#' @param y Numeric vector of directing function values (length = number of vertices).
#' @param Z Numeric matrix or data frame of response function values.
#'   Number of rows must equal the number of vertices.
#' @param type Character scalar specifying measure type: "normalized", "slope", or "sign".
#' @param y.diff.type Character scalar specifying edge difference type for y.
#' @param z.diff.type Character scalar specifying edge difference type for Z columns.
#' @param epsilon Numeric scalar for pseudocount in log-ratios (0 = adaptive).
#' @param sigmoid.alpha Numeric scalar for sigmoid normalization scale (0 = auto).
#' @param ascending Logical. If TRUE (default), use ascending gradient.
#' @param n.threads Integer specifying number of OpenMP threads (0 = default).
#'
#' @return An object of class "lslope_vector_matrix_result". A numeric matrix of
#'   dimension (n.vertices x n.columns) where column j contains local slope
#'   coefficients for \code{Z[,j]} with respect to y.
#'
#' @details
#' This function uses an optimized C++ implementation that:
#' \enumerate{
#'   \item Pre-computes gradient edges from y once
#'   \item Calibrates sigmoid alpha from the data (if type = "normalized")
#'   \item Processes all columns of Z in parallel using OpenMP
#' }
#'
#' The algorithm complexity is O(n*d + q*n) where n is the number of vertices,
#' d is average vertex degree, and q is the number of columns in Z. This is
#' more efficient than the R-level loop which has O(q*n*d) complexity.
#'
#' @seealso \code{\link{lslope}} for the unified interface
#'
#' @export
lslope.vector.matrix <- function(adj.list,
                                 weight.list,
                                 y,
                                 Z,
                                 type = c("normalized", "slope", "sign"),
                                 y.diff.type = c("difference", "logratio"),
                                 z.diff.type = c("difference", "logratio"),
                                 epsilon = 0,
                                 sigmoid.alpha = 0,
                                 ascending = TRUE,
                                 n.threads = 0L) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Input validation
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.numeric(y)) stop("y must be a numeric vector")

    if (is.data.frame(Z)) Z <- as.matrix(Z)
    if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix")

    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")
    if (!is.numeric(sigmoid.alpha) || length(sigmoid.alpha) != 1)
        stop("sigmoid.alpha must be a single numeric value")
    if (!is.logical(ascending) || length(ascending) != 1)
        stop("ascending must be a single logical value")
    if (!is.numeric(n.threads) || length(n.threads) != 1)
        stop("n.threads must be a single numeric value")

    n.vertices <- length(adj.list)
    n.columns <- ncol(Z)

    if (length(y) != n.vertices)
        stop(sprintf("Length of y (%d) must equal number of vertices (%d)",
                     length(y), n.vertices))
    if (nrow(Z) != n.vertices)
        stop(sprintf("nrow(Z) (%d) must equal number of vertices (%d)",
                     nrow(Z), n.vertices))
    if (length(weight.list) != n.vertices)
        stop(sprintf("Length of weight.list (%d) must equal number of vertices (%d)",
                     length(weight.list), n.vertices))

    ## Store column names for output
    z.col.names <- colnames(Z)
    if (is.null(z.col.names)) z.col.names <- paste0("Z", seq_len(n.columns))

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    result <- .Call(
        "S_lslope_vector_matrix",
        adj.list.0,
        weight.list,
        as.numeric(y),
        Z,  ## Matrix is passed directly (column-major like R expects)
        as.character(type),
        as.character(y.diff.type),
        as.character(z.diff.type),
        as.numeric(epsilon),
        as.numeric(sigmoid.alpha),
        as.logical(ascending),
        as.integer(n.threads),
        PACKAGE = "gflow"
    )

    ## Extract coefficient matrix and add column names
    coef.mat <- result$coefficients
    colnames(coef.mat) <- z.col.names

    ## Add class and attributes
    class(coef.mat) <- c("lslope_vector_matrix_result", "matrix", "array")
    attr(coef.mat, "type") <- type
    attr(coef.mat, "y.diff.type") <- y.diff.type
    attr(coef.mat, "z.diff.type") <- z.diff.type
    attr(coef.mat, "ascending") <- ascending
    attr(coef.mat, "n.vertices") <- n.vertices
    attr(coef.mat, "n.columns") <- n.columns
    attr(coef.mat, "sigmoid.alpha") <- result$sigmoid.alpha
    attr(coef.mat, "n.local.maxima") <- result$n.local.maxima
    attr(coef.mat, "n.local.minima") <- result$n.local.minima

    ## Store gradient structure as attributes for potential debugging
    attr(coef.mat, "gradient.neighbors") <- result$gradient.neighbors
    attr(coef.mat, "gradient.delta.y") <- result$gradient.delta.y
    attr(coef.mat, "is.local.extremum") <- result$is.local.extremum

    return(coef.mat)
}

#' Local Slope Between a Vector and Matrix Columns (serial S_lslope_gradient call implementation)
#'
#' Compute vertex-level slope coefficients measuring the response of each column
#' of Z to a directing function y on a graph.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#' @param weight.list List of numeric vectors containing edge weights.
#' @param y Numeric vector of directing function values (length = number of vertices).
#' @param Z Numeric matrix or data frame of response function values.
#'   Number of rows must equal the number of vertices.
#' @param type Character scalar specifying measure type: "normalized", "slope", or "sign".
#' @param y.diff.type Character scalar specifying edge difference type for y.
#' @param z.diff.type Character scalar specifying edge difference type for Z columns.
#' @param epsilon Numeric scalar for pseudocount in log-ratios (0 = adaptive).
#' @param sigmoid.alpha Numeric scalar for sigmoid normalization scale (0 = auto).
#' @param ascending Logical. If TRUE (default), use ascending gradient.
#'
#' @return An object of class "lslope_vector_matrix_result". A numeric matrix of
#'   dimension (n.vertices x n.columns) where column j contains local slope
#'   coefficients for \code{Z[,j]} with respect to y. The matrix has column names
#'   from Z (or "Z1", "Z2", ... if Z has no column names).
#'
#' @details
#' This function computes the local slope of each column of Z with respect to
#' a single directing function y. The gradient edges are determined once from y
#' and reused for all columns of Z, making this more efficient than repeated
#' calls to the vector-vector version.
#'
#' At vertex v, for each column j:
#' \deqn{lslope(Z_j; y)(v) = f(\Delta_{grad} Z_j / \Delta_{grad} y)}
#' where f is identity (slope), tanh (normalized), or sign, and
#' \eqn{\Delta_{grad}}
#' denotes the difference along the gradient edge of y.
#'
#' @seealso \code{\link{lslope}} for the unified interface
#'
#' @export
lslope.vector.matrix.R <- function(adj.list,
                                        weight.list,
                                        y,
                                        Z,
                                        type = c("normalized", "slope", "sign"),
                                        y.diff.type = c("difference", "logratio"),
                                        z.diff.type = c("difference", "logratio"),
                                        epsilon = 0,
                                        sigmoid.alpha = 0,
                                        ascending = TRUE) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Input validation
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.numeric(y)) stop("y must be a numeric vector")

    if (is.data.frame(Z)) Z <- as.matrix(Z)
    if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix")

    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")
    if (!is.numeric(sigmoid.alpha) || length(sigmoid.alpha) != 1)
        stop("sigmoid.alpha must be a single numeric value")
    if (!is.logical(ascending) || length(ascending) != 1)
        stop("ascending must be a single logical value")

    n.vertices <- length(adj.list)
    n.columns <- ncol(Z)

    if (length(y) != n.vertices)
        stop(sprintf("Length of y (%d) must equal number of vertices (%d)",
                     length(y), n.vertices))
    if (nrow(Z) != n.vertices)
        stop(sprintf("nrow(Z) (%d) must equal number of vertices (%d)",
                     nrow(Z), n.vertices))
    if (length(weight.list) != n.vertices)
        stop(sprintf("Length of weight.list (%d) must equal number of vertices (%d)",
                     length(weight.list), n.vertices))

    ## Store column names for output
    z.col.names <- colnames(Z)
    if (is.null(z.col.names)) z.col.names <- paste0("Z", seq_len(n.columns))

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Get gradient edges and extrema info from first call (instrumented)
    ## This establishes the gradient structure from y
    first.result <- .Call(
        "S_lslope_gradient_instrumented",
        adj.list.0,
        weight.list,
        as.numeric(y),
        as.numeric(Z[, 1]),
        as.character(type),
        as.character(y.diff.type),
        as.character(z.diff.type),
        as.numeric(epsilon),
        as.numeric(sigmoid.alpha),
        "tanh",
        as.logical(ascending),
        PACKAGE = "gflow"
    )

    ## Extract gradient structure (same for all columns since y is fixed)
    is.extremum <- first.result$is.local.extremum
    gradient.neighbors <- first.result$gradient.neighbors
    gradient.delta.y <- first.result$gradient.delta.y
    calibrated.sigmoid.alpha <- first.result$sigmoid.alpha

    ## Initialize result matrix
    coef.mat <- matrix(NA_real_, nrow = n.vertices, ncol = n.columns)
    coef.mat[, 1] <- first.result$vertex.coefficients

    ## Process remaining columns
    ## Since gradient edges depend only on y, we can use non-instrumented calls
    ## but we need to ensure consistent sigmoid.alpha across columns
    for (j in seq_len(n.columns)[-1]) {
        result <- .Call(
            "S_lslope_gradient",
            adj.list.0,
            weight.list,
            as.numeric(y),
            as.numeric(Z[, j]),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(calibrated.sigmoid.alpha),  # Use calibrated alpha for consistency
            "tanh",
            as.logical(ascending),
            PACKAGE = "gflow"
        )
        coef.mat[, j] <- result
    }

    ## Add column names
    colnames(coef.mat) <- z.col.names

    ## Add class and attributes
    class(coef.mat) <- c("lslope_vector_matrix_result", "matrix", "array")
    attr(coef.mat, "type") <- type
    attr(coef.mat, "y.diff.type") <- y.diff.type
    attr(coef.mat, "z.diff.type") <- z.diff.type
    attr(coef.mat, "ascending") <- ascending
    attr(coef.mat, "n.vertices") <- n.vertices
    attr(coef.mat, "n.columns") <- n.columns
    attr(coef.mat, "sigmoid.alpha") <- calibrated.sigmoid.alpha
    attr(coef.mat, "n.local.maxima") <- first.result$n.local.maxima
    attr(coef.mat, "n.local.minima") <- first.result$n.local.minima

    return(coef.mat)
}


################################################################################
#
# Case (c): Matrix-Matrix Implementation
#
################################################################################

#' Local Slope Between All Column Pairs of Two Matrices
#'
#' Compute vertex-level slope coefficients for all pairs of columns from
#' directing matrix Y and response matrix Z.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#' @param weight.list List of numeric vectors containing edge weights.
#' @param Y Numeric matrix or data frame of directing function values.
#'   Number of rows must equal the number of vertices.
#' @param Z Numeric matrix or data frame of response function values.
#'   Number of rows must equal the number of vertices.
#' @param type Character scalar specifying measure type: "normalized", "slope", or "sign".
#' @param y.diff.type Character scalar specifying edge difference type for Y columns.
#' @param z.diff.type Character scalar specifying edge difference type for Z columns.
#' @param epsilon Numeric scalar for pseudocount in log-ratios (0 = adaptive).
#' @param sigmoid.alpha Numeric scalar for sigmoid normalization scale (0 = auto).
#' @param ascending Logical. If TRUE (default), use ascending gradient.
#' @param mc.cores Integer specifying number of cores for parallel computation.
#'   Default is 1 (sequential). Note: parallelization via mclapply is not
#'   available on Windows.
#'
#' @return An object of class "lslope_matrix_matrix_result". A 3D numeric array
#'     of dimension (n.vertices x ncol(Y) x ncol(Z)). Element \code{[v, i, j]}
#'     contains the local slope of \code{Z[,j]} with respect to \code{Y[,i]} at
#'     vertex v.
#'
#' @details
#' This function computes local slopes for all \code{(Y_i, Z_j)} pairs. Since
#' the local slope is asymmetric (treating Y as the directing variable and Z as
#' the response), all ncol(Y) x ncol(Z) pairs must be computed; there is no
#' symmetry optimization.
#'
#' The computation can be parallelized via mc.cores when running on Unix-like
#' systems. On Windows, the function falls back to sequential execution.
#'
#' @section Memory Considerations:
#' For p columns in Y and q columns in Z, the result requires storage for
#' n.vertices x p x q coefficients. For large matrices, consider processing
#' in batches or using a subset of columns.
#'
#' @seealso \code{\link{lslope}} for the unified interface,
#'   \code{\link{lslope.vector.matrix}} for vector-matrix computation
#'
#' @export
lslope.matrix.matrix <- function(adj.list,
                                      weight.list,
                                      Y,
                                      Z,
                                      type = c("normalized", "slope", "sign"),
                                      y.diff.type = c("difference", "logratio"),
                                      z.diff.type = c("difference", "logratio"),
                                      epsilon = 0,
                                      sigmoid.alpha = 0,
                                      ascending = TRUE,
                                      mc.cores = 1L) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Convert data frames to matrices
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    if (is.data.frame(Z)) Z <- as.matrix(Z)

    ## Input validation
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix")
    if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix")
    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")
    if (!is.numeric(sigmoid.alpha) || length(sigmoid.alpha) != 1)
        stop("sigmoid.alpha must be a single numeric value")
    if (!is.logical(ascending) || length(ascending) != 1)
        stop("ascending must be a single logical value")
    if (!is.numeric(mc.cores) || length(mc.cores) != 1 || mc.cores < 1)
        stop("mc.cores must be a positive integer")

    mc.cores <- as.integer(mc.cores)
    n.vertices <- length(adj.list)
    n.cols.y <- ncol(Y)
    n.cols.z <- ncol(Z)

    if (nrow(Y) != n.vertices)
        stop(sprintf("nrow(Y) (%d) must equal number of vertices (%d)",
                     nrow(Y), n.vertices))
    if (nrow(Z) != n.vertices)
        stop(sprintf("nrow(Z) (%d) must equal number of vertices (%d)",
                     nrow(Z), n.vertices))
    if (length(weight.list) != n.vertices)
        stop(sprintf("Length of weight.list (%d) must equal number of vertices (%d)",
                     length(weight.list), n.vertices))

    ## Get column names
    y.col.names <- colnames(Y)
    if (is.null(y.col.names)) y.col.names <- paste0("Y", seq_len(n.cols.y))

    z.col.names <- colnames(Z)
    if (is.null(z.col.names)) z.col.names <- paste0("Z", seq_len(n.cols.z))

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Build list of all (i, j) pairs to compute
    ## Since lslope is asymmetric, we compute all p * q pairs
    pairs <- vector("list", n.cols.y * n.cols.z)
    idx <- 1
    for (i in seq_len(n.cols.y)) {
        for (j in seq_len(n.cols.z)) {
            pairs[[idx]] <- c(i, j)
            idx <- idx + 1
        }
    }
    n.pairs <- length(pairs)

    ## Define the worker function for a single pair
    compute.pair <- function(pair) {
        i <- pair[1]
        j <- pair[2]

        result <- .Call(
            "S_lslope_gradient",
            adj.list.0,
            weight.list,
            as.numeric(Y[, i]),
            as.numeric(Z[, j]),
            as.character(type),
            as.character(y.diff.type),
            as.character(z.diff.type),
            as.numeric(epsilon),
            as.numeric(sigmoid.alpha),
            "tanh",
            as.logical(ascending),
            PACKAGE = "gflow"
        )

        return(result)
    }

    ## Execute computation (parallel or sequential)
    if (mc.cores > 1 && .Platform$OS.type != "windows") {
        results <- parallel::mclapply(pairs, compute.pair, mc.cores = mc.cores)
    } else {
        results <- lapply(pairs, compute.pair)
    }

    ## Assemble results into 3D array
    coef.array <- array(NA_real_,
                        dim = c(n.vertices, n.cols.y, n.cols.z),
                        dimnames = list(NULL, y.col.names, z.col.names))

    for (k in seq_along(pairs)) {
        i <- pairs[[k]][1]
        j <- pairs[[k]][2]
        coef.array[, i, j] <- results[[k]]
    }

    ## Add class and attributes
    class(coef.array) <- c("lslope_matrix_matrix_result", "array")
    attr(coef.array, "type") <- type
    attr(coef.array, "y.diff.type") <- y.diff.type
    attr(coef.array, "z.diff.type") <- z.diff.type
    attr(coef.array, "ascending") <- ascending
    attr(coef.array, "n.vertices") <- n.vertices
    attr(coef.array, "n.cols.y") <- n.cols.y
    attr(coef.array, "n.cols.z") <- n.cols.z

    return(coef.array)
}


################################################################################
#
# Print and Summary Methods
#
################################################################################

#' Print Method for Gradient-Restricted Local Slope Results
#'
#' @param x An object of class "lslope_gradient_result"
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_gradient_result <- function(x, ...) {
    cat("Gradient-Restricted Local Slope Result\n")
    cat("=======================================\n\n")

    cat("Parameters:\n")
    cat("  Type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Ascending gradient:", attr(x, "ascending"), "\n\n")

    cat("Results:\n")
    cat("  Number of vertices:", x$n.vertices, "\n")
    cat("  Local maxima:", x$n.local.maxima, "\n")
    cat("  Local minima:", x$n.local.minima, "\n")
    cat("  Mean coefficient:", round(x$mean.coefficient, 4), "\n")
    cat("  Median coefficient:", round(x$median.coefficient, 4), "\n")

    if (attr(x, "type") == "normalized") {
        cat("  Sigmoid alpha:", round(x$sigmoid.alpha, 4), "\n")
    }

    cat("\nCoefficient distribution (non-extrema):\n")
    non.extrema <- x$vertex.coefficients[!x$is.local.extremum]
    if (length(non.extrema) > 0) {
        print(summary(non.extrema))
    }

    invisible(x)
}


#' Print Method for Vector-Matrix Local Slope Results
#'
#' @param x An object of class "lslope_vector_matrix_result"
#' @param digits Number of digits for printing
#' @param max.show Maximum number of columns to display
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_vector_matrix_result <- function(x, digits = 4, max.show = 10, ...) {
    cat("Local Slope Vector-Matrix Result\n")
    cat("=================================\n\n")

    cat("Parameters:\n")
    cat("  Type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Ascending gradient:", attr(x, "ascending"), "\n")
    cat("  Vertices:", attr(x, "n.vertices"),
        ", Columns:", attr(x, "n.columns"), "\n")
    cat("  Local extrema:", attr(x, "n.local.maxima") + attr(x, "n.local.minima"),
        "(", attr(x, "n.local.maxima"), "maxima,",
        attr(x, "n.local.minima"), "minima )\n")

    if (attr(x, "type") == "normalized") {
        cat("  Sigmoid alpha:", round(attr(x, "sigmoid.alpha"), 4), "\n")
    }
    cat("\n")

    ## Show summary statistics
    col.means <- colMeans(x)
    n.show <- min(length(col.means), max.show)

    cat("Column mean coefficients (first", n.show, "):\n")
    print(round(col.means[seq_len(n.show)], digits))

    if (length(col.means) > max.show) {
        cat("  ... and", length(col.means) - max.show, "more columns\n")
    }

    cat("\nOverall summary of column means:\n")
    print(summary(col.means))

    invisible(x)
}


#' Summary Method for Vector-Matrix Local Slope Results
#'
#' @param object An object of class "lslope_vector_matrix_result"
#' @param ... Additional arguments (ignored)
#' @export
summary.lslope_vector_matrix_result <- function(object, ...) {
    cat("Local Slope Vector-Matrix Summary\n")
    cat("==================================\n\n")

    col.means <- colMeans(object)
    col.medians <- apply(object, 2, median)

    cat("Column mean coefficients:\n")
    print(summary(col.means))

    cat("\nColumn median coefficients:\n")
    print(summary(col.medians))

    ## Top features by absolute mean
    cat("\nTop 5 columns by |mean coefficient|:\n")
    n.cols <- ncol(object)
    if (n.cols > 0) {
        top.idx <- order(abs(col.means), decreasing = TRUE)[seq_len(min(5, n.cols))]
        top.df <- data.frame(
            Column = colnames(object)[top.idx],
            Mean = round(col.means[top.idx], 4),
            Median = round(col.medians[top.idx], 4),
            PropPositive = round(colMeans(object[, top.idx, drop = FALSE] > 1e-10), 3)
        )
        print(top.df, row.names = FALSE)
    }

    invisible(object)
}


#' Print Method for Matrix-Matrix Local Slope Results
#'
#' @param x An object of class "lslope_matrix_matrix_result"
#' @param digits Number of digits for printing
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_matrix_matrix_result <- function(x, digits = 4, ...) {
    cat("Local Slope Matrix-Matrix Result\n")
    cat("=================================\n\n")

    cat("Parameters:\n")
    cat("  Type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Ascending gradient:", attr(x, "ascending"), "\n\n")

    dims <- dim(x)
    cat("Dimensions:", dims[1], "vertices x",
        dims[2], "Y-columns x", dims[3], "Z-columns\n")
    cat("Total pairs:", dims[2] * dims[3], "\n\n")

    ## Overall summary
    cat("Summary of all coefficients:\n")
    print(summary(as.vector(x)))

    ## Mean across vertices for each (Y, Z) pair
    pair.means <- apply(x, c(2, 3), mean)
    cat("\nMean coefficient matrix (Y-cols as rows, Z-cols as cols):\n")
    print(round(pair.means[seq_len(min(5, dims[2])),
                           seq_len(min(5, dims[3]))], digits))
    if (dims[2] > 5 || dims[3] > 5) {
        cat("  ... (showing first 5x5 block)\n")
    }

    invisible(x)
}


#' Summary Method for Matrix-Matrix Local Slope Results
#'
#' @param object An object of class "lslope_matrix_matrix_result"
#' @param ... Additional arguments (ignored)
#' @export
summary.lslope_matrix_matrix_result <- function(object, ...) {
    cat("Local Slope Matrix-Matrix Summary\n")
    cat("==================================\n\n")

    dims <- dim(object)

    cat("Array dimensions:\n")
    cat("  Vertices:", dims[1], "\n")
    cat("  Y columns:", dims[2], "\n")
    cat("  Z columns:", dims[3], "\n")
    cat("  Total pairs:", dims[2] * dims[3], "\n\n")

    ## Compute mean matrix
    pair.means <- apply(object, c(2, 3), mean)

    cat("Summary of vertex-averaged coefficients:\n")
    print(summary(as.vector(pair.means)))

    ## Strongest associations
    cat("\nTop 5 (Y, Z) pairs by |mean coefficient|:\n")
    abs.means <- abs(pair.means)
    top.idx <- order(abs.means, decreasing = TRUE)[1:min(5, length(abs.means))]
    top.coords <- arrayInd(top.idx, dim(pair.means))

    top.df <- data.frame(
        Y = dimnames(object)[[2]][top.coords[, 1]],
        Z = dimnames(object)[[3]][top.coords[, 2]],
        Mean = round(pair.means[top.idx], 4)
    )
    print(top.df, row.names = FALSE)

    invisible(object)
}


################################################################################
#
# Neighborhood Local Regression Coefficient (unchanged from original)
#
################################################################################

#' Neighborhood Local Regression Coefficient
#'
#' @description
#' Computes the local regression coefficient of z on y using all edges in each
#' vertex's neighborhood. This is the asymmetric analog of local correlation,
#' satisfying beta_loc(z; y) = lcor(y, z) * sd_loc(z) / sd_loc(y).
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#' @param weight.list List of numeric vectors containing edge weights.
#' @param y Numeric vector of directing function values.
#' @param z Numeric vector of response function values.
#' @param weight.type Character scalar specifying weighting scheme:
#'   \describe{
#'     \item{"derivative"}{Geometric weights (w_e = 1/length^2). Default.}
#'     \item{"unit"}{Equal weights (w_e = 1).}
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y.
#' @param z.diff.type Character scalar specifying edge difference type for z.
#' @param epsilon Numeric scalar for pseudocount in log-ratios.
#'
#' @return A list with class "lslope_neighborhood_result" containing:
#'   \describe{
#'     \item{vertex.coefficients}{Numeric vector of local regression coefficients.}
#'     \item{sd.y}{Numeric vector of local standard deviation of y.}
#'     \item{sd.z}{Numeric vector of local standard deviation of z.}
#'     \item{lcor}{Numeric vector of local correlation (for reference).}
#'     \item{mean.coefficient}{Mean coefficient across all vertices.}
#'     \item{median.coefficient}{Median coefficient across all vertices.}
#'   }
#'
#' @details
#' The neighborhood local regression coefficient at vertex v is:
#' \deqn{\beta_{loc}(z; y, w)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z}{\sum w_e (\Delta_e y)^2}}
#'
#' Unlike the gradient-restricted slope which uses only the single gradient edge,
#' this measure uses all edges in the neighborhood, providing a more robust but
#' less direction-specific measure of local association.
#'
#' The relationship beta = rho * sigma_z / sigma_y holds locally: the regression
#' coefficient equals the correlation times the ratio of local standard deviations.
#'
#' @examples
#' \dontrun{
#' result <- lslope.neighborhood(adj.list, weight.list, y, z,
#'                               weight.type = "derivative",
#'                               y.diff.type = "difference",
#'                               z.diff.type = "logratio")
#'
#' ## Verify relationship with local correlation
#' beta.from.lcor <- result$lcor * result$sd.z / result$sd.y
#' all.equal(result$vertex.coefficients, beta.from.lcor)
#' }
#'
#' @seealso \code{\link{lslope}} for gradient-restricted slope,
#'   \code{\link{lcor}} for symmetric local correlation
#'
#' @export
lslope.neighborhood <- function(adj.list,
                                weight.list,
                                y,
                                z,
                                weight.type = c("derivative", "unit"),
                                y.diff.type = c("difference", "logratio"),
                                z.diff.type = c("difference", "logratio"),
                                epsilon = 0) {

    ## Match arguments
    weight.type <- match.arg(weight.type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Validate inputs
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.numeric(y)) stop("y must be a numeric vector")
    if (!is.numeric(z)) stop("z must be a numeric vector")
    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")

    n.vertices <- length(adj.list)

    if (length(y) != n.vertices) stop("Length of y must equal number of vertices")
    if (length(z) != n.vertices) stop("Length of z must equal number of vertices")
    if (length(weight.list) != n.vertices)
        stop("Length of weight.list must equal number of vertices")

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    result <- .Call(
        "S_lslope_neighborhood",
        adj.list.0,
        weight.list,
        as.numeric(y),
        as.numeric(z),
        as.character(weight.type),
        as.character(y.diff.type),
        as.character(z.diff.type),
        as.numeric(epsilon),
        PACKAGE = "gflow"
    )

    ## Add class for method dispatch
    class(result) <- c("lslope_neighborhood_result", "list")

    ## Add metadata
    attr(result, "weight.type") <- weight.type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "n.vertices") <- n.vertices

    return(result)
}


#' Print Method for Neighborhood Local Slope Results
#'
#' @param x An object of class "lslope_neighborhood_result"
#' @param ... Additional arguments (ignored)
#' @export
print.lslope_neighborhood_result <- function(x, ...) {
    cat("Neighborhood Local Regression Coefficient Result\n")
    cat("=================================================\n\n")

    cat("Parameters:\n")
    cat("  Weight type:", attr(x, "weight.type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n\n")

    cat("Results:\n")
    cat("  Number of vertices:", attr(x, "n.vertices"), "\n")
    cat("  Mean coefficient:", round(x$mean.coefficient, 4), "\n")
    cat("  Median coefficient:", round(x$median.coefficient, 4), "\n")

    cat("\nCoefficient distribution:\n")
    print(summary(x$vertex.coefficients))

    cat("\nLocal correlation summary:\n")
    print(summary(x$lcor))

    invisible(x)
}

#' Local Correlation Coefficients on Graphs
#'
#' Compute vertex-level correlation coefficients measuring the alignment of
#' directional changes between functions defined on graph vertices. This function
#' provides a unified interface supporting vector-vector, vector-matrix, and
#' matrix-matrix computations with automatic dispatch based on input types.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights (lengths).
#'   Must have same structure as adj.list.
#' @param y Numeric vector, matrix, or data frame of function values.
#'   For vectors, length must equal the number of vertices.
#'   For matrices/data frames, number of rows must equal the number of vertices.
#' @param z Numeric vector, matrix, or data frame of function values.
#'   Dimension requirements match those for y.
#' @param type Character scalar specifying the weighting scheme:
#'   \describe{
#'     \item{"derivative"}{Geometric weights (w_e = 1/length^2). Normalizes by
#'       edge length squared, making the measure analogous to comparing derivatives.
#'       Appropriate when functions represent continuous quantities sampled at
#'       irregular spatial positions. This is the default.}
#'     \item{"unit"}{Equal weights (w_e = 1). Treats all edges equally regardless
#'       of length. Appropriate when edge lengths are comparable or when counting
#'       directional agreements without geometric normalization.}
#'     \item{"sign"}{Sign-based weighting. Uses only the sign of edge products,
#'       providing robustness to outliers at the cost of magnitude information.}
#'   }
#' @param y.diff.type Character scalar specifying edge difference type for y:
#'   \describe{
#'     \item{"difference"}{Standard differences: Delta_e f = f(u) - f(v).
#'       Appropriate for continuous data in Euclidean space.}
#'     \item{"logratio"}{Log-ratios: Delta_e f = log((f(u) + epsilon) / (f(v) + epsilon)).
#'       Appropriate for compositional data (relative abundances, proportions).
#'       The log transformation maps multiplicative changes to an additive scale,
#'       corresponding to the Aitchison distance on the simplex.}
#'   }
#' @param z.diff.type Character scalar specifying edge difference type for z.
#'   Same options as y.diff.type.
#' @param epsilon Numeric scalar for pseudocount in log-ratio transformations.
#'   If 0 (default), computed adaptively as 1e-6 times the minimum non-zero value.
#'   Only used when the corresponding diff.type is "logratio".
#' @param winsorize.quantile Numeric scalar for winsorization of edge differences.
#'   If 0 (default), no winsorization is applied.
#'   If positive (e.g., 0.025), clips edge differences to the \eqn{[q, 1-q]} percentile
#'   range for robustness against outliers.
#' @param instrumented Logical scalar. If TRUE and both y and z are vectors,
#'   returns additional diagnostic information including edge-level differences
#'   and weights. Ignored for matrix inputs. Default is FALSE.
#' @param mc.cores Integer scalar specifying the number of cores for parallel
#'   computation in the matrix-matrix case. Default is 1 (sequential).
#'   Note that parallel::mclapply is not available on Windows.
#'
#' @return The return type depends on the input types:
#'
#'   \strong{Vector-vector (y and z both vectors):}
#'   An object of class "lcor_vector_vector_result".
#'   If instrumented = FALSE, a numeric vector of correlation coefficients
#'   at each vertex, with values in \eqn{[-1, 1]}.
#'   If instrumented = TRUE, a list containing:
#'   \describe{
#'     \item{vertex.coefficients}{Numeric vector of local correlation coefficients.}
#'     \item{vertex.delta.y}{List of edge differences for y at each vertex.}
#'     \item{vertex.delta.z}{List of edge differences for z at each vertex.}
#'     \item{vertex.weights}{List of edge weights at each vertex.}
#'     \item{all.delta.y}{All y edge differences (pre-winsorization).}
#'     \item{all.delta.z}{All z edge differences (pre-winsorization).}
#'     \item{y.lower, y.upper}{Winsorization bounds for y.}
#'     \item{z.lower, z.upper}{Winsorization bounds for z.}
#'   }
#'
#'   \strong{Vector-matrix (y vector, z matrix) or matrix-vector (y matrix, z vector):}
#'   An object of class "lcor_vector_matrix_result".
#'   If instrumented = FALSE, a numeric matrix of dimension (n.vertices x n.columns)
#'   where column j contains local correlation coefficients between y and \code{z[,j]}.
#'   If instrumented = TRUE, a list containing:
#'   \describe{
#'     \item{column.coefficients}{Matrix (n.vertices x n.columns) of local
#'       correlations. Column j contains \code{lcor(y, z[,j])} at each vertex.}
#'     \item{y.lower, y.upper}{Scalar winsorization bounds for y edge differences.}
#'     \item{z.lower, z.upper}{Numeric vectors of length n.columns giving
#'       per-column winsorization bounds for z edge differences.}
#'   }
#'   For matrix-vector input (y matrix, z vector), the attribute "transposed"
#'   is set to TRUE, and column j corresponds to \code{lcor(y[,j], z)}.
#'
#'   \strong{Matrix-matrix (y and z both matrices):}
#'   A list with class "lcor_matrix_matrix_result" containing:
#'   \describe{
#'     \item{pair.coefficients}{List of numeric vectors, one per column pair.
#'       Each vector contains vertex-wise coefficients.}
#'     \item{mean.coefficients}{Numeric vector of mean coefficients per pair.}
#'     \item{pair.names}{Character vector of pair names in "col_y:col_z" format.}
#'     \item{pair.indices}{Two-column integer matrix of (i, j) pair indices.}
#'     \item{symmetric}{Logical indicating whether y and z are identical
#'       (only upper triangular pairs computed).}
#'   }
#'
#' @details
#' The local correlation coefficient at vertex v measures the alignment of
#' directional changes in y and z within v's neighborhood:
#'
#' \deqn{lcor(y,z)(v) = \frac{\sum_e w_e \Delta_e y \cdot \Delta_e z}
#'                           {\sqrt{\sum_e w_e (\Delta_e y)^2}
#'                            \sqrt{\sum_e w_e (\Delta_e z)^2}}}
#'
#' where the sum is over edges e incident to vertex v, and Delta_e f denotes the
#' edge difference (either standard or log-ratio). This correlation-style
#' normalization makes the coefficient scale-invariant and interpretable as the
#' cosine of the angle between gradient vectors in the appropriate geometry.
#'
#' @section Dispatch Behavior:
#'
#' The function automatically dispatches to specialized implementations based on
#' input types:
#'
#' \strong{Case (a) Vector-Vector:} Computes a single local correlation between
#' two functions. This is the fundamental operation underlying all other cases.
#'
#' \strong{Case (b) Vector-Matrix:} Efficiently computes correlations between one
#' vector y and each column of matrix Z. The implementation pre-computes y-dependent
#' quantities once and reuses them across columns, providing significant speedup
#' over repeated vector-vector calls.
#'
#' \strong{Case (c) Matrix-Vector:} Equivalent to case (b) with roles swapped.
#' The function internally calls the vector-matrix implementation with y and z
#' exchanged and diff.types appropriately swapped.
#'
#' \strong{Case (d) Matrix-Matrix:} Computes correlations between all pairs of
#' columns from y and z. When y and z are identical (checked via identical()),
#' only upper triangular pairs are computed to avoid redundancy. Parallel
#' computation via mc.cores can accelerate this case substantially.
#'
#' @section Edge Difference Types:
#'
#' The flexibility to specify different edge difference types enables appropriate
#' treatment of mixed data:
#'
#' \strong{Continuous vs. compositional:} Use "difference" for y (e.g., clinical
#' severity score) and "logratio" for z (e.g., bacterial relative abundances).
#'
#' \strong{Compositional vs. compositional:} Use "logratio" for both when
#' analyzing co-variation of relative abundances in microbiome data.
#'
#' \strong{Continuous vs. continuous:} Use "difference" for both when analyzing
#' standard numeric measurements.
#'
#' @section Geometric Interpretation:
#'
#' With derivative weighting, the local correlation coefficient converges to
#' cos(theta) where theta is the angle between gradient vectors (or their
#' log-ratio analogs for compositional data). Values near +1 indicate parallel
#' gradients (functions increase together), values near -1 indicate anti-parallel
#' gradients (one increases while the other decreases), and values near 0
#' indicate orthogonal gradients (no directional relationship).
#'
#' @examples
#' \dontrun{
#'
#' # Build a simple graph (triangle)
#' adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
#' weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
#'
#' # Case (a): Vector-vector
#' y <- c(1.0, 2.0, 3.0)
#' z <- c(2.0, 4.0, 5.0)
#' result <- lcor(adj.list, weight.list, y, z, type = "unit")
#' print(result$vertex.coefficients)
#'
#' # Case (b): Vector-matrix (feature screening)
#' # Compare response against multiple bacterial abundances
#' response <- rnorm(100)
#' abundances <- matrix(runif(100 * 50), nrow = 100, ncol = 50)
#' colnames(abundances) <- paste0("ASV", 1:50)
#'
#' # Assume graph is built from the data
#' result <- lcor(adj.list, weight.list,
#'                response, abundances,
#'                type = "derivative",
#'                y.diff.type = "difference",
#'                z.diff.type = "logratio")
#'
#' # Identify top features by absolute mean correlation
#' top.features <- order(abs(result$mean.coefficients), decreasing = TRUE)[1:10]
#' print(names(result$mean.coefficients)[top.features])
#'
#' # Case (c): Matrix-vector (same as b with roles swapped)
#' result <- lcor(adj.list, weight.list,
#'                abundances, response,
#'                type = "derivative",
#'                y.diff.type = "logratio",
#'                z.diff.type = "difference")
#'
#' # Case (d): Matrix-matrix (pairwise correlations)
#' # Compute local correlation tensor for compositional data
#' Z <- abundances[, 1:10]  # First 10 ASVs
#' result <- lcor(adj.list, weight.list,
#'                Z, Z,
#'                type = "unit",
#'                y.diff.type = "logratio",
#'                z.diff.type = "logratio",
#'                mc.cores = 4)
#'
#' # Result contains upper triangular pairs only since y == z
#' print(result$pair.names[1:5])
#'
#' # With winsorization for robustness
#' result <- lcor(adj.list, weight.list, y, z,
#'                type = "derivative",
#'                winsorize.quantile = 0.025)
#'
#' # With instrumented output for diagnostics
#' result <- lcor(adj.list, weight.list, y, z,
#'                type = "derivative",
#'                instrumented = TRUE)
#' hist(result$all.delta.y, main = "Distribution of y edge differences")
#' }
#'
#' @seealso
#' \code{\link{lcor.vector.matrix}} for direct vector-matrix computation,
#' \code{\link{lcor.matrix.matrix}} for direct matrix-matrix computation,
#' \code{\link{comono}} for co-monotonicity coefficients with ratio normalization
#'
#' @references
#' Pawel Gajer and Jacques Ravel (2025). From Global to Local Correlation: Geometric Decomposition of Statistical Inference. \emph{arXiv:2511.04599}.
#'
#' @export
lcor <- function(adj.list,
                 weight.list,
                 y,
                 z,
                 type = c("derivative", "unit", "sign"),
                 y.diff.type = c("difference", "logratio"),
                 z.diff.type = c("difference", "logratio"),
                 epsilon = 0,
                 winsorize.quantile = 0,
                 instrumented = FALSE,
                 mc.cores = 1L) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Detect input types FIRST
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
    if (!is.numeric(winsorize.quantile) || length(winsorize.quantile) != 1)
        stop("winsorize.quantile must be a single numeric value")

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
        ## Case a): vector-vector
        result <- lcor.vector.vector(adj.list, weight.list, y, z,
                                     type, y.diff.type, z.diff.type,
                                     epsilon, winsorize.quantile,
                                     instrumented)
        return(result)

    } else if (!y.is.matrix && z.is.matrix) {
        ## Case b): vector-matrix
        result <- lcor.vector.matrix(adj.list, weight.list, y, z,
                                     type, y.diff.type, z.diff.type,
                                     epsilon, winsorize.quantile)
        return(result)

    } else if (y.is.matrix && !z.is.matrix) {
        ## Case c): matrix-vector (transpose of case b)
        result <- lcor.vector.matrix(adj.list, weight.list,
                                     z, y,           ## swap y and z
                                     type,
                                     z.diff.type,    ## swap diff types
                                     y.diff.type,
                                     epsilon, winsorize.quantile)
        attr(result, "transposed") <- TRUE
        return(result)

    } else {
        ## Case d): matrix-matrix
        result <- lcor.matrix.matrix(adj.list, weight.list, y, z,
                                     type, y.diff.type, z.diff.type,
                                     epsilon, winsorize.quantile,
                                     mc.cores)
        return(result)
    }
}

## Internal helper for vector-vector case (extracts current logic)
lcor.vector.vector <- function(adj.list, weight.list, y, z,
                                type, y.diff.type, z.diff.type,
                                epsilon, winsorize.quantile,
                                instrumented) {
    n.vertices <- length(adj.list)
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    if (instrumented) {
        result <- .Call("S_lcor_instrumented",
                        adj.list.0,
                        weight.list,
                        as.numeric(y),
                        as.numeric(z),
                        type,
                        y.diff.type,
                        z.diff.type,
                        as.numeric(epsilon),
                        as.numeric(winsorize.quantile),
                        PACKAGE = "gflow")
    } else {
        result <- .Call("S_lcor",
                        adj.list.0,
                        weight.list,
                        as.numeric(y),
                        as.numeric(z),
                        type,
                        y.diff.type,
                        z.diff.type,
                        as.numeric(epsilon),
                        as.numeric(winsorize.quantile),
                        PACKAGE = "gflow")
    }

    class(result) <- c("lcor_vector_vector_result", "vector")
    attr(result, "type") <- type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "epsilon") <- epsilon
    attr(result, "winsorize.quantile") <- winsorize.quantile
    attr(result, "n.vertices") <- n.vertices

    return(result)
}

#' Local Correlation Between a Vector and Matrix Columns
#'
#' Compute vertex-level correlation coefficients measuring the alignment of
#' directional changes between a response vector y and each column of a
#' feature matrix Z on a graph.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights.
#'   Must have same structure as adj.list.
#' @param y Numeric vector of response function values (length = number of vertices).
#' @param Z Numeric matrix or data frame of feature function values.
#'   Number of rows must equal the number of vertices.
#' @param type Character scalar specifying weighting scheme:
#'   "derivative" (default), "unit", or "sign".
#' @param y.diff.type Character scalar specifying edge difference type for y:
#'   "difference" (default) or "logratio".
#' @param z.diff.type Character scalar specifying edge difference type for Z columns:
#'   "difference" (default) or "logratio".
#' @param epsilon Numeric scalar for pseudocount in log-ratios (0 = adaptive).
#' @param winsorize.quantile Numeric scalar for winsorization (0 = none).
#' @param instrumented Logical. If FALSE (default), returns coefficient matrix only.
#'   If TRUE, returns list with coefficients and winsorization bounds.
#'
#' @return Depends on the instrumented parameter:
#'
#'   If instrumented = FALSE: A numeric matrix of dimension (n.vertices x n.columns)
#'   where column j contains local correlation coefficients between y and \code{Z[,j]}.
#'   The matrix has class "lcor_vector_matrix_result" and column names from Z.
#'
#'   If instrumented = TRUE: A list with class "lcor_vector_matrix_result" containing:
#'   \describe{
#'     \item{column.coefficients}{Matrix (n.vertices x n.columns) of local
#'       correlations. Column j contains \code{lcor(y, Z[,j])} at each vertex.}
#'     \item{y.lower}{Scalar lower winsorization bound for y edge differences
#'       (-Inf if no winsorization).}
#'     \item{y.upper}{Scalar upper winsorization bound for y edge differences
#'       (+Inf if no winsorization).}
#'     \item{z.lower}{Numeric vector of length n.columns giving per-column
#'       lower winsorization bounds for z edge differences.}
#'     \item{z.upper}{Numeric vector of length n.columns giving per-column
#'       upper winsorization bounds for z edge differences.}
#'   }
#'
#' @details
#' The local correlation coefficient at vertex v for column j measures the
#' alignment of directional changes in y and \code{Z[,j]} within v's neighborhood:
#'
#' \deqn{lcor(y, z_j)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z_j}
#'                              {\sqrt{\sum w_e (\Delta_e y)^2}
#'                               \sqrt{\sum w_e (\Delta_e z_j)^2}}}
#'
#' This function is optimized for the case where one response vector y is
#' compared against many feature columns. The y-dependent quantities are
#' computed once and reused across all columns.
#'
#' @section Computational Efficiency:
#'
#' For q columns and a graph with n vertices and m edges:
#' \itemize{
#'   \item One-time setup: O(m) to pre-compute y-dependent quantities
#'   \item Per-column processing: O(m) to compute z edge differences
#'   \item Total: O(m + q*m)
#' }
#'
#' @examples
#' \dontrun{
#' library(gflow)
#'
#' # Feature screening: response vs. multiple abundances
#' result <- lcor.vector.matrix(
#'   adj.list, weight.list,
#'   response, abundances,
#'   type = "derivative",
#'   y.diff.type = "difference",
#'   z.diff.type = "logratio"
#' )
#'
#' # Result is a matrix; get summary statistics easily
#' col.means <- colMeans(result)
#' top.features <- order(abs(col.means), decreasing = TRUE)[1:10]
#'
#' # With instrumented output for diagnostics
#' result <- lcor.vector.matrix(
#'   adj.list, weight.list,
#'   response, abundances,
#'   winsorize.quantile = 0.025,
#'   instrumented = TRUE
#' )
#'
#' # Examine winsorization bounds
#' cat("Y bounds:", result$y.lower, "to", result$y.upper, "\n")
#' cat("Z bounds (first 5 cols):\n")
#' print(head(data.frame(lower = result$z.lower, upper = result$z.upper)))
#' }
#'
#' @seealso
#' \code{\link{lcor}} for the unified interface supporting all input combinations
#'
#' @export
lcor.vector.matrix <- function(adj.list,
                                weight.list,
                                y,
                                Z,
                                type = c("derivative", "unit", "sign"),
                                y.diff.type = c("difference", "logratio"),
                                z.diff.type = c("difference", "logratio"),
                                epsilon = 0,
                                winsorize.quantile = 0.025,
                                instrumented = FALSE) {

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
    if (!is.numeric(winsorize.quantile) || length(winsorize.quantile) != 1)
        stop("winsorize.quantile must be a single numeric value")
    if (!is.logical(instrumented) || length(instrumented) != 1)
        stop("instrumented must be a single logical value")

    n.vertices <- length(adj.list)

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
    if (is.null(z.col.names)) z.col.names <- paste0("Z", seq_len(ncol(Z)))

    ## Convert to 0-based indexing
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++ function
    result <- .Call(
        "S_lcor_vector_matrix",
        adj.list.0,
        weight.list,
        as.numeric(y),
        Z,
        as.character(type),
        as.character(y.diff.type),
        as.character(z.diff.type),
        as.numeric(epsilon),
        as.numeric(winsorize.quantile),
        as.logical(instrumented),
        PACKAGE = "gflow"
    )

    ## Process result based on instrumented flag
    if (!instrumented) {
        ## Result is a matrix - add column names and class
        colnames(result) <- z.col.names
        class(result) <- c("lcor_vector_matrix_result", "matrix", "array")
    } else {
        ## Result is a list - add column names to coefficient matrix and bounds
        colnames(result$column.coefficients) <- z.col.names
        names(result$z.lower) <- z.col.names
        names(result$z.upper) <- z.col.names
        class(result) <- c("lcor_vector_matrix_result", "list")
    }

    ## Add metadata as attributes
    attr(result, "type") <- type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "epsilon") <- epsilon
    attr(result, "winsorize.quantile") <- winsorize.quantile
    attr(result, "n.vertices") <- n.vertices
    attr(result, "n.columns") <- ncol(Z)
    attr(result, "instrumented") <- instrumented

    return(result)
}


#' Print Method for lcor_vector_matrix_result
#'
#' @param x An object of class "lcor_vector_matrix_result"
#' @param digits Number of digits for printing
#' @param max.show Maximum number of columns to display in summary
#' @param ... Additional arguments (ignored)
#' @export
print.lcor_vector_matrix_result <- function(x, digits = 4, max.show = 10, ...) {
    cat("Local Correlation Vector-Matrix Result\n")
    cat("=======================================\n\n")

    cat("Parameters:\n")
    cat("  Weighting type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Vertices:", attr(x, "n.vertices"),
        ", Columns:", attr(x, "n.columns"), "\n")
    cat("  Instrumented:", attr(x, "instrumented"), "\n\n")

    ## Get coefficient matrix
    if (attr(x, "instrumented")) {
        coef.mat <- x$column.coefficients
    } else {
        coef.mat <- x
    }

    ## Show summary statistics
    col.means <- colMeans(coef.mat)
    n.show <- min(length(col.means), max.show)

    cat("Column mean coefficients (first", n.show, "):\n")
    print(round(col.means[seq_len(n.show)], digits))

    if (length(col.means) > max.show) {
        cat("  ... and", length(col.means) - max.show, "more columns\n")
    }

    cat("\nOverall summary of column means:\n")
    print(summary(col.means))

    ## Show winsorization bounds if instrumented
    if (attr(x, "instrumented") && attr(x, "winsorize.quantile") > 0) {
        cat("\nWinsorization bounds:\n")
        cat("  Y: [", round(x$y.lower, digits), ", ",
            round(x$y.upper, digits), "]\n", sep = "")
        cat("  Z (range across columns): [",
            round(min(x$z.lower), digits), ", ",
            round(max(x$z.upper), digits), "]\n", sep = "")
    }

    invisible(x)
}


#' Summary Method for lcor_vector_matrix_result
#'
#' @param object An object of class "lcor_vector_matrix_result"
#' @param ... Additional arguments (ignored)
#' @export
summary.lcor_vector_matrix_result <- function(object, ...) {
    cat("Local Correlation Vector-Matrix Summary\n")
    cat("========================================\n\n")

    ## Get coefficient matrix
    if (attr(object, "instrumented")) {
        coef.mat <- object$column.coefficients
    } else {
        coef.mat <- object
    }

    col.means <- colMeans(coef.mat)
    col.medians <- apply(coef.mat, 2, median)

    cat("Column mean coefficients:\n")
    print(summary(col.means))

    cat("\nColumn median coefficients:\n")
    print(summary(col.medians))

    ## Top features by absolute mean
    cat("\nTop 5 columns by |mean coefficient|:\n")
    n.cols <- ncol(coef.mat)
    if (n.cols > 0) {
        top.idx <- order(abs(col.means), decreasing = TRUE)[seq_len(min(5, n.cols))]
        top.df <- data.frame(
            Column = colnames(coef.mat)[top.idx],
            Mean = round(col.means[top.idx], 4),
            Median = round(col.medians[top.idx], 4),
            PropPositive = round(colMeans(coef.mat[, top.idx, drop = FALSE] > 1e-10), 3)
        )
        print(top.df, row.names = FALSE)
    }

    invisible(object)
}

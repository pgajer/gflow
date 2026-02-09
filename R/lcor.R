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
#'       The log transformation maps multiplicative changes to an additive scale.}
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
#' @param hop.radius Integer scalar specifying the hop distance for local
#'   neighborhoods. Default is 1 (immediate neighbors). When hop.radius > 1,
#'   neighborhoods include all vertices reachable within hop.radius hops.
#'   For \code{type = "derivative"}, edge lengths are replaced by shortest path
#'   lengths (sum of edge lengths) among paths with at most hop.radius hops.
#'   For \code{type = "unit"} or \code{type = "sign"}, all k-hop edges are
#'   treated as having length 1.
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
#'   where column j contains local correlation coefficients between the vector
#'   and the j-th column of the matrix.
#'   If instrumented = TRUE, a list containing:
#'   \describe{
#'     \item{column.coefficients}{Matrix (n.vertices x n.columns) of local
#'       correlations.}
#'     \item{y.lower, y.upper}{Scalar winsorization bounds for the vector input.}
#'     \item{z.lower, z.upper}{Numeric vectors of length n.columns giving
#'       per-column winsorization bounds for the matrix input.}
#'   }
#'   For matrix-vector input, the attribute "transposed" is set to TRUE.
#'
#'   \strong{Matrix-matrix with identical(y, z) = TRUE (symmetric case):}
#'   An object of class "lcor_matrix_matrix_result" with attribute symmetric = TRUE.
#'   Only upper triangular pairs (i, j) with i < j are computed.
#'   If instrumented = FALSE, a numeric matrix of dimension (n.vertices x n.pairs)
#'   where n.pairs = ncol(y) * (ncol(y) - 1) / 2. Column names are in
#'   "col_i:col_j" format.
#'   If instrumented = TRUE, a list containing:
#'   \describe{
#'     \item{pair.coefficients}{Matrix (n.vertices x n.pairs) of local correlations.}
#'     \item{pair.names}{Character vector of pair names in "col_i:col_j" format.}
#'     \item{pair.indices}{Integer matrix (n.pairs x 2) of column index pairs.}
#'     \item{column.lower, column.upper}{Numeric vectors of winsorization bounds
#'       for each column (shared between y and z since they are identical).}
#'   }
#'
#'   \strong{Matrix-matrix with identical(y, z) = FALSE (asymmetric case):}
#'   An object of class "lcor_matrix_matrix_result" with attribute symmetric = FALSE.
#'   All pairs (i, j) for i in 1:ncol(y) and j in 1:ncol(z) are computed.
#'   If instrumented = FALSE, a 3D numeric array of dimension
#'   (n.vertices x ncol(y) x ncol(z)). Element \code{[v, i, j]} contains
#'   \code{lcor(y[,i], z[,j])} at vertex v.
#'   If instrumented = TRUE, a list containing:
#'   \describe{
#'     \item{pair.coefficients}{3D array (n.vertices x ncol(y) x ncol(z)) of
#'       local correlations.}
#'     \item{y.column.lower, y.column.upper}{Numeric vectors of winsorization
#'       bounds for each column of y.}
#'     \item{z.column.lower, z.column.upper}{Numeric vectors of winsorization
#'       bounds for each column of z.}
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
#' @section Hop Radius:
#'
#' By default, the local correlation at vertex v is computed using edges incident
#' to v (hop.radius = 1). For hop.radius > 1, the neighborhood expands to all
#' vertices within hop.radius hops of v. The computation is then performed on the
#' k-hop star graph centered at v. For \code{type = "derivative"}, k-hop edges are
#' weighted by the inverse squared shortest path length within the hop limit; for
#' \code{type = "unit"} and \code{type = "sign"}, all k-hop edges are treated as
#' having unit length.
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
                 mc.cores = 1L,
                 hop.radius = 1L) {

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
                                     instrumented, hop.radius)
        return(result)

    } else if (!y.is.matrix && z.is.matrix) {
        ## Case b): vector-matrix
        result <- lcor.vector.matrix(adj.list, weight.list, y, z,
                                     type, y.diff.type, z.diff.type,
                                     epsilon, winsorize.quantile,
                                     instrumented = FALSE,
                                     hop.radius = hop.radius)
        return(result)

    } else if (y.is.matrix && !z.is.matrix) {
        ## Case c): matrix-vector (transpose of case b)
        result <- lcor.vector.matrix(adj.list, weight.list,
                                     z, y,           ## swap y and z
                                     type,
                                     z.diff.type,    ## swap diff types
                                     y.diff.type,
                                     epsilon, winsorize.quantile,
                                     instrumented = FALSE,
                                     hop.radius = hop.radius)
        attr(result, "transposed") <- TRUE
        return(result)

    } else {
        ## Case d): matrix-matrix
        result <- lcor.matrix.matrix(adj.list, weight.list, y, z,
                                     type, y.diff.type, z.diff.type,
                                     epsilon, winsorize.quantile,
                                     instrumented = FALSE,
                                     mc.cores = mc.cores,
                                     hop.radius = hop.radius)
        return(result)
    }
}

## Internal helper to expand neighborhoods to hop radius
.prepare_lcor_hop_graph <- function(adj.list, weight.list, hop.radius, type) {
    if (!is.numeric(hop.radius) || length(hop.radius) != 1 ||
        hop.radius < 1 || hop.radius != as.integer(hop.radius)) {
        stop("hop.radius must be a positive integer")
    }

    hop.radius <- as.integer(hop.radius)
    if (hop.radius <= 1L) {
        return(list(adj.list = adj.list, weight.list = weight.list, hop.radius = hop.radius))
    }

    if (type == "derivative") {
        pg <- create.path.graph(adj.list, weight.list, h = hop.radius)
        adj.list <- pg$adj.list
        weight.list <- pg$edge.length.list
    } else {
        unit.weights <- lapply(adj.list, function(neighbors) {
            if (length(neighbors) > 0) rep(1.0, length(neighbors)) else numeric(0)
        })
        pg <- create.path.graph(adj.list, unit.weights, h = hop.radius)
        adj.list <- pg$adj.list
        weight.list <- lapply(adj.list, function(neighbors) {
            if (length(neighbors) > 0) rep(1.0, length(neighbors)) else numeric(0)
        })
    }

    return(list(adj.list = adj.list, weight.list = weight.list, hop.radius = hop.radius))
}

## Internal helper for vector-vector case (extracts current logic)
lcor.vector.vector <- function(adj.list, weight.list, y, z,
                                type, y.diff.type, z.diff.type,
                                epsilon, winsorize.quantile,
                                instrumented, hop.radius = 1L) {
    hop.graph <- .prepare_lcor_hop_graph(adj.list, weight.list, hop.radius, type)
    adj.list <- hop.graph$adj.list
    weight.list <- hop.graph$weight.list
    hop.radius <- hop.graph$hop.radius

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
    attr(result, "hop.radius") <- hop.radius

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
#' @param hop.radius Integer scalar specifying the hop distance for local
#'   neighborhoods. Default is 1 (immediate neighbors). When hop.radius > 1,
#'   neighborhoods include all vertices reachable within hop.radius hops.
#'   For \code{type = "derivative"}, edge lengths are replaced by shortest path
#'   lengths (sum of edge lengths) among paths with at most hop.radius hops.
#'   For \code{type = "unit"} or \code{type = "sign"}, all k-hop edges are
#'   treated as having length 1.
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
#' When hop.radius > 1, the neighborhood expands to all vertices within the
#' specified hop distance, and the computation is performed on the resulting
#' k-hop star graph.
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
                                instrumented = FALSE,
                                hop.radius = 1L) {

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

    hop.graph <- .prepare_lcor_hop_graph(adj.list, weight.list, hop.radius, type)
    adj.list <- hop.graph$adj.list
    weight.list <- hop.graph$weight.list
    hop.radius <- hop.graph$hop.radius

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
    attr(result, "hop.radius") <- hop.radius

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
    cat("  Hop radius:", attr(x, "hop.radius"), "\n")
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

################################################################################
#
# lcor.matrix.matrix() Implementation
#
# Computes local correlation between all pairs of columns from matrices y and z.
# - Symmetric case (y == z): computes upper triangular pairs only
# - Asymmetric case (y != z): computes all pairs
#
################################################################################


#' Local Correlation Between All Column Pairs of Two Matrices
#'
#' Compute vertex-level correlation coefficients between all pairs of columns
#' from matrices y and z. When y and z are identical, only upper triangular
#' pairs are computed to avoid redundancy.
#'
#' @param adj.list List of integer vectors containing 1-based vertex indices.
#'   Element i contains the neighbors of vertex i.
#' @param weight.list List of numeric vectors containing edge weights.
#'   Must have same structure as adj.list.
#' @param y Numeric matrix or data frame of function values.
#'   Number of rows must equal the number of vertices.
#' @param z Numeric matrix or data frame of function values.
#'   Number of rows must equal the number of vertices.
#' @param type Character scalar specifying weighting scheme:
#'   "derivative" (default), "unit", or "sign".
#' @param y.diff.type Character scalar specifying edge difference type for y columns:
#'   "difference" (default) or "logratio".
#' @param z.diff.type Character scalar specifying edge difference type for z columns:
#'   "difference" (default) or "logratio".
#' @param epsilon Numeric scalar for pseudocount in log-ratios (0 = adaptive).
#' @param winsorize.quantile Numeric scalar for winsorization (0 = none).
#' @param instrumented Logical. If FALSE (default), returns coefficient array only.
#'   If TRUE, returns list with coefficients and winsorization bounds.
#' @param mc.cores Integer specifying number of cores for parallel computation.
#'   Default is 1 (sequential). Note: parallelization via mclapply is not
#'   available on Windows.
#' @param hop.radius Integer scalar specifying the hop distance for local
#'   neighborhoods. Default is 1 (immediate neighbors). When hop.radius > 1,
#'   neighborhoods include all vertices reachable within hop.radius hops.
#'   For \code{type = "derivative"}, edge lengths are replaced by shortest path
#'   lengths (sum of edge lengths) among paths with at most hop.radius hops.
#'   For \code{type = "unit"} or \code{type = "sign"}, all k-hop edges are
#'   treated as having length 1.
#'
#' @return Depends on whether y and z are identical and the instrumented parameter:
#'
#'   \strong{Symmetric case (identical(y, z) is TRUE):}
#'
#'   If instrumented = FALSE: A numeric matrix of dimension (n.vertices x n.pairs)
#'   where n.pairs = ncol(y) * (ncol(y) - 1) / 2. Column k contains local
#'   correlation coefficients for the k-th upper triangular pair.
#'   Column names are in "col_i:col_j" format.
#'
#'   If instrumented = TRUE: A list containing:
#'   \describe{
#'     \item{pair.coefficients}{Matrix (n.vertices x n.pairs) of local correlations.}
#'     \item{pair.names}{Character vector of pair names in "col_i:col_j" format.}
#'     \item{pair.indices}{Integer matrix (n.pairs x 2) of column index pairs (i, j).}
#'     \item{column.lower}{Numeric vector of lower winsorization bounds for each
#'       column of y (= z).}
#'     \item{column.upper}{Numeric vector of upper winsorization bounds for each
#'       column of y (= z).}
#'   }
#'
#'   \strong{Asymmetric case (identical(y, z) is FALSE):}
#'
#'   If instrumented = FALSE: A 3D numeric array of dimension
#'   (n.vertices x ncol(y) x ncol(z)). Element \code{[v, i, j]} contains the local
#'   correlation coefficient \code{lcor(y[,i], z[,j])} at vertex v.
#'
#'   If instrumented = TRUE: A list containing:
#'   \describe{
#'     \item{pair.coefficients}{3D array (n.vertices x ncol(y) x ncol(z)) of
#'       local correlations.}
#'     \item{y.column.lower}{Numeric vector of lower winsorization bounds for
#'       each column of y.}
#'     \item{y.column.upper}{Numeric vector of upper winsorization bounds for
#'       each column of y.}
#'     \item{z.column.lower}{Numeric vector of lower winsorization bounds for
#'       each column of z.}
#'     \item{z.column.upper}{Numeric vector of upper winsorization bounds for
#'       each column of z.}
#'   }
#'
#'   In both cases, the result has class "lcor_matrix_matrix_result" and
#'   attribute "symmetric" indicating which case applies.
#'
#' @details
#' This function computes local correlations between all relevant pairs of
#' columns from y and z. For each pair (i, j), it computes:
#'
#' \deqn{lcor(y_i, z_j)(v) = \frac{\sum w_e \Delta_e y_i \cdot \Delta_e z_j}
#'                                {\sqrt{\sum w_e (\Delta_e y_i)^2}
#'                                 \sqrt{\sum w_e (\Delta_e z_j)^2}}}
#'
#' at each vertex v.
#'
#' When hop.radius > 1, the neighborhood expands to all vertices within the
#' specified hop distance, and the computation is performed on the resulting
#' k-hop star graph.
#'
#' @section Symmetric Detection:
#'
#' The function uses \code{identical(y, z)} to detect the symmetric case.
#' This checks for object identity, not just numerical equality. If you want
#' symmetric treatment for numerically equal but distinct objects, pass the
#' same object for both y and z.
#'
#' @section Parallelization:
#'
#' When mc.cores > 1, the function uses \code{parallel::mclapply} to distribute
#' pair computations across cores. This provides near-linear speedup for large
#' numbers of pairs. Note that mclapply uses forking, which is not available
#' on Windows; on Windows, the function falls back to sequential execution
#' regardless of mc.cores.
#'
#' @section Memory Considerations:
#'
#' For p columns in y and q columns in z, the result requires storage for
#' n.vertices * p * q coefficients in the asymmetric case, or
#' n.vertices * p * (p-1) / 2 in the symmetric case. For large matrices,
#' consider processing in batches.
#'
#' @examples
#' \dontrun{
#' library(gflow)
#'
#' # Symmetric case: local correlation tensor for compositional data
#' Z <- abundances[, 1:20]  # First 20 ASVs
#' result <- lcor.matrix.matrix(
#'   adj.list, weight.list,
#'   Z, Z,
#'   type = "unit",
#'   y.diff.type = "logratio",
#'   z.diff.type = "logratio",
#'   mc.cores = 4
#' )
#'
#' # Result is matrix with 20*19/2 = 190 pairs
#' dim(result)  # n.vertices x 190
#'
#' # Find pairs with strongest mean correlation
#' pair.means <- colMeans(result)
#' top.pairs <- order(abs(pair.means), decreasing = TRUE)[1:10]
#' print(colnames(result)[top.pairs])
#'
#' # Asymmetric case: correlations between two different feature sets
#' Y <- clinical.scores    # n x 5 matrix of clinical variables
#' Z <- abundances[, 1:50] # n x 50 matrix of abundances
#' result <- lcor.matrix.matrix(
#'   adj.list, weight.list,
#'   Y, Z,
#'   type = "derivative",
#'   y.diff.type = "difference",
#'   z.diff.type = "logratio"
#' )
#'
#' # Result is 3D array
#' dim(result)  # n.vertices x 5 x 50
#'
#' # Mean correlation between clinical variable 2 and all ASVs
#' mean.cors <- colMeans(result[, 2, ])
#' }
#'
#' @seealso
#' \code{\link{lcor}} for the unified interface,
#' \code{\link{lcor.vector.matrix}} for vector-matrix computation
#'
#' @export
lcor.matrix.matrix <- function(adj.list,
                                weight.list,
                                y,
                                z,
                                type = c("derivative", "unit", "sign"),
                                y.diff.type = c("difference", "logratio"),
                                z.diff.type = c("difference", "logratio"),
                                epsilon = 0,
                                winsorize.quantile = 0,
                                instrumented = FALSE,
                                mc.cores = 1L,
                                hop.radius = 1L) {

    ## Match arguments
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)

    ## Convert data frames to matrices
    if (is.data.frame(y)) y <- as.matrix(y)
    if (is.data.frame(z)) z <- as.matrix(z)

    ## Input validation
    if (!is.list(adj.list)) stop("adj.list must be a list")
    if (!is.list(weight.list)) stop("weight.list must be a list")
    if (!is.matrix(y) || !is.numeric(y)) stop("y must be a numeric matrix")
    if (!is.matrix(z) || !is.numeric(z)) stop("z must be a numeric matrix")
    if (!is.numeric(epsilon) || length(epsilon) != 1)
        stop("epsilon must be a single numeric value")
    if (!is.numeric(winsorize.quantile) || length(winsorize.quantile) != 1)
        stop("winsorize.quantile must be a single numeric value")
    if (!is.logical(instrumented) || length(instrumented) != 1)
        stop("instrumented must be a single logical value")
    if (!is.numeric(mc.cores) || length(mc.cores) != 1 || mc.cores < 1)
        stop("mc.cores must be a positive integer")

    mc.cores <- as.integer(mc.cores)
    n.vertices <- length(adj.list)
    n.cols.y <- ncol(y)
    n.cols.z <- ncol(z)

    if (nrow(y) != n.vertices)
        stop(sprintf("nrow(y) (%d) must equal number of vertices (%d)",
                     nrow(y), n.vertices))
    if (nrow(z) != n.vertices)
        stop(sprintf("nrow(z) (%d) must equal number of vertices (%d)",
                     nrow(z), n.vertices))
    if (length(weight.list) != n.vertices)
        stop(sprintf("Length of weight.list (%d) must equal number of vertices (%d)",
                     length(weight.list), n.vertices))

    hop.graph <- .prepare_lcor_hop_graph(adj.list, weight.list, hop.radius, type)
    adj.list <- hop.graph$adj.list
    weight.list <- hop.graph$weight.list
    hop.radius <- hop.graph$hop.radius

    ## Get column names
    y.col.names <- colnames(y)
    if (is.null(y.col.names)) y.col.names <- paste0("Y", seq_len(n.cols.y))

    z.col.names <- colnames(z)
    if (is.null(z.col.names)) z.col.names <- paste0("Z", seq_len(n.cols.z))

    ## Convert to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Detect symmetric case
    symmetric <- identical(y, z)

    ## Build list of pairs to compute
    if (symmetric) {
        ## Upper triangular pairs only: (1,2), (1,3), ..., (m-1, m)
        if (n.cols.y < 2) {
            stop("For symmetric case, y must have at least 2 columns")
        }
        pairs <- combn(seq_len(n.cols.y), 2, simplify = FALSE)
        n.pairs <- length(pairs)
        pair.names <- sapply(pairs, function(p) {
            paste0(y.col.names[p[1]], ":", y.col.names[p[2]])
        })
    } else {
        ## All pairs: (1,1), (1,2), ..., (p, q)
        pairs <- vector("list", n.cols.y * n.cols.z)
        pair.names <- character(n.cols.y * n.cols.z)
        idx <- 1
        for (i in seq_len(n.cols.y)) {
            for (j in seq_len(n.cols.z)) {
                pairs[[idx]] <- c(i, j)
                pair.names[idx] <- paste0(y.col.names[i], ":", z.col.names[j])
                idx <- idx + 1
            }
        }
        n.pairs <- length(pairs)
    }

    ## Define the worker function for a single pair
    compute.pair <- function(pair) {
        i <- pair[1]
        j <- pair[2]

        if (instrumented) {
            ## Call instrumented version to get bounds
            res <- .Call(
                "S_lcor_instrumented",
                adj.list.0,
                weight.list,
                as.numeric(y[, i]),
                as.numeric(z[, j]),
                as.character(type),
                as.character(y.diff.type),
                as.character(z.diff.type),
                as.numeric(epsilon),
                as.numeric(winsorize.quantile),
                PACKAGE = "gflow"
            )
            return(res)
        } else {
            ## Call non-instrumented version for coefficients only
            res <- .Call(
                "S_lcor",
                adj.list.0,
                weight.list,
                as.numeric(y[, i]),
                as.numeric(z[, j]),
                as.character(type),
                as.character(y.diff.type),
                as.character(z.diff.type),
                as.numeric(epsilon),
                as.numeric(winsorize.quantile),
                PACKAGE = "gflow"
            )
            return(res)
        }
    }

    ## Execute computation (parallel or sequential)
    if (mc.cores > 1 && .Platform$OS.type != "windows") {
        results <- parallel::mclapply(pairs, compute.pair, mc.cores = mc.cores)
    } else {
        results <- lapply(pairs, compute.pair)
    }

    ## Assemble results based on symmetric/asymmetric and instrumented flags
    if (symmetric) {
        result <- .assemble.symmetric.result(
            results, pairs, pair.names, n.vertices, n.cols.y,
            y.col.names, instrumented
        )
    } else {
        result <- .assemble.asymmetric.result(
            results, pairs, n.vertices, n.cols.y, n.cols.z,
            y.col.names, z.col.names, instrumented
        )
    }

    ## Add common attributes
    attr(result, "type") <- type
    attr(result, "y.diff.type") <- y.diff.type
    attr(result, "z.diff.type") <- z.diff.type
    attr(result, "epsilon") <- epsilon
    attr(result, "winsorize.quantile") <- winsorize.quantile
    attr(result, "n.vertices") <- n.vertices
    attr(result, "symmetric") <- symmetric
    attr(result, "instrumented") <- instrumented
    attr(result, "hop.radius") <- hop.radius

    class(result) <- c("lcor_matrix_matrix_result",
                       if (instrumented) "list" else class(result))

    return(result)
}


## Internal helper: assemble symmetric case results
.assemble.symmetric.result <- function(results, pairs, pair.names, n.vertices,
                                        n.cols, col.names, instrumented) {
    n.pairs <- length(pairs)

    if (!instrumented) {
        ## Simple case: just coefficient matrix
        coef.mat <- matrix(NA_real_, nrow = n.vertices, ncol = n.pairs)
        for (k in seq_len(n.pairs)) {
            coef.mat[, k] <- results[[k]]
        }
        colnames(coef.mat) <- pair.names
        return(coef.mat)

    } else {
        ## Instrumented: extract coefficients and bounds
        coef.mat <- matrix(NA_real_, nrow = n.vertices, ncol = n.pairs)
        pair.indices <- matrix(NA_integer_, nrow = n.pairs, ncol = 2)

        ## Track bounds for each column (may appear in multiple pairs)
        ## Use first occurrence for each column's bounds
        column.lower <- rep(NA_real_, n.cols)
        column.upper <- rep(NA_real_, n.cols)
        names(column.lower) <- col.names
        names(column.upper) <- col.names

        for (k in seq_len(n.pairs)) {
            res <- results[[k]]
            i <- pairs[[k]][1]
            j <- pairs[[k]][2]

            coef.mat[, k] <- res$vertex.coefficients
            pair.indices[k, ] <- c(i, j)

            ## Store bounds for columns i and j if not yet recorded
            ## y bounds correspond to column i, z bounds to column j
            if (is.na(column.lower[i])) {
                column.lower[i] <- res$y.lower
                column.upper[i] <- res$y.upper
            }
            if (is.na(column.lower[j])) {
                column.lower[j] <- res$z.lower
                column.upper[j] <- res$z.upper
            }
        }

        colnames(coef.mat) <- pair.names
        colnames(pair.indices) <- c("i", "j")

        return(list(
            pair.coefficients = coef.mat,
            pair.names = pair.names,
            pair.indices = pair.indices,
            column.lower = column.lower,
            column.upper = column.upper
        ))
    }
}


## Internal helper: assemble asymmetric case results
.assemble.asymmetric.result <- function(results, pairs, n.vertices,
                                         n.cols.y, n.cols.z,
                                         y.col.names, z.col.names,
                                         instrumented) {

    if (!instrumented) {
        ## Simple case: 3D array
        coef.array <- array(NA_real_,
                            dim = c(n.vertices, n.cols.y, n.cols.z),
                            dimnames = list(NULL, y.col.names, z.col.names))

        for (k in seq_along(pairs)) {
            i <- pairs[[k]][1]
            j <- pairs[[k]][2]
            coef.array[, i, j] <- results[[k]]
        }
        return(coef.array)

    } else {
        ## Instrumented: 3D array plus per-column bounds
        coef.array <- array(NA_real_,
                            dim = c(n.vertices, n.cols.y, n.cols.z),
                            dimnames = list(NULL, y.col.names, z.col.names))

        ## Track bounds for each column of y and z
        y.column.lower <- rep(NA_real_, n.cols.y)
        y.column.upper <- rep(NA_real_, n.cols.y)
        z.column.lower <- rep(NA_real_, n.cols.z)
        z.column.upper <- rep(NA_real_, n.cols.z)
        names(y.column.lower) <- y.col.names
        names(y.column.upper) <- y.col.names
        names(z.column.lower) <- z.col.names
        names(z.column.upper) <- z.col.names

        for (k in seq_along(pairs)) {
            res <- results[[k]]
            i <- pairs[[k]][1]
            j <- pairs[[k]][2]

            coef.array[, i, j] <- res$vertex.coefficients

            ## Store bounds if not yet recorded
            if (is.na(y.column.lower[i])) {
                y.column.lower[i] <- res$y.lower
                y.column.upper[i] <- res$y.upper
            }
            if (is.na(z.column.lower[j])) {
                z.column.lower[j] <- res$z.lower
                z.column.upper[j] <- res$z.upper
            }
        }

        return(list(
            pair.coefficients = coef.array,
            y.column.lower = y.column.lower,
            y.column.upper = y.column.upper,
            z.column.lower = z.column.lower,
            z.column.upper = z.column.upper
        ))
    }
}


#' Print Method for lcor_matrix_matrix_result
#'
#' @param x An object of class "lcor_matrix_matrix_result"
#' @param digits Number of digits for printing
#' @param ... Additional arguments (ignored)
#' @export
print.lcor_matrix_matrix_result <- function(x, digits = 4, ...) {
    cat("Local Correlation Matrix-Matrix Result\n")
    cat("=======================================\n\n")

    symmetric <- attr(x, "symmetric")
    instrumented <- attr(x, "instrumented")

    cat("Parameters:\n")
    cat("  Weighting type:", attr(x, "type"), "\n")
    cat("  Y difference type:", attr(x, "y.diff.type"), "\n")
    cat("  Z difference type:", attr(x, "z.diff.type"), "\n")
    cat("  Symmetric:", symmetric, "\n")
    cat("  Instrumented:", instrumented, "\n")
    cat("  Vertices:", attr(x, "n.vertices"), "\n")
    cat("  Hop radius:", attr(x, "hop.radius"), "\n\n")

    ## Get coefficient structure
    if (instrumented) {
        coeffs <- x$pair.coefficients
    } else {
        coeffs <- x
    }

    if (symmetric) {
        ## Matrix case
        n.pairs <- ncol(coeffs)
        cat("Pairs computed:", n.pairs, "(upper triangular)\n\n")

        pair.means <- colMeans(coeffs)
        cat("Summary of pair mean coefficients:\n")
        print(summary(pair.means))

        ## Top pairs
        cat("\nTop 5 pairs by |mean coefficient|:\n")
        top.idx <- order(abs(pair.means), decreasing = TRUE)[seq_len(min(5, n.pairs))]
        top.df <- data.frame(
            Pair = colnames(coeffs)[top.idx],
            Mean = round(pair.means[top.idx], digits)
        )
        print(top.df, row.names = FALSE)

    } else {
        ## 3D array case
        dims <- dim(coeffs)
        cat("Dimensions:", dims[1], "vertices x",
            dims[2], "y-columns x", dims[3], "z-columns\n")
        cat("Total pairs:", dims[2] * dims[3], "\n\n")

        ## Overall summary
        cat("Summary of all coefficients:\n")
        print(summary(as.vector(coeffs)))

        ## Mean across vertices for each (y, z) pair
        pair.means <- apply(coeffs, c(2, 3), mean)
        cat("\nMean coefficient matrix (y-cols as rows, z-cols as cols):\n")
        print(round(pair.means[seq_len(min(5, dims[2])),
                               seq_len(min(5, dims[3]))], digits))
        if (dims[2] > 5 || dims[3] > 5) {
            cat("  ... (showing first 5x5 block)\n")
        }
    }

    invisible(x)
}


#' Summary Method for lcor_matrix_matrix_result
#'
#' @param object An object of class "lcor_matrix_matrix_result"
#' @param ... Additional arguments (ignored)
#' @export
summary.lcor_matrix_matrix_result <- function(object, ...) {
    cat("Local Correlation Matrix-Matrix Summary\n")
    cat("========================================\n\n")

    symmetric <- attr(object, "symmetric")
    instrumented <- attr(object, "instrumented")

    if (instrumented) {
        coeffs <- object$pair.coefficients
    } else {
        coeffs <- object
    }

    if (symmetric) {
        pair.means <- colMeans(coeffs)
        pair.medians <- apply(coeffs, 2, median)

        cat("Pair statistics:\n")
        cat("  Number of pairs:", length(pair.means), "\n")
        cat("  Mean of pair means:", round(mean(pair.means), 4), "\n")
        cat("  SD of pair means:", round(sd(pair.means), 4), "\n\n")

        cat("Distribution of pair mean coefficients:\n")
        print(summary(pair.means))

        if (instrumented && attr(object, "winsorize.quantile") > 0) {
            cat("\nWinsorization bounds by column:\n")
            bounds.df <- data.frame(
                Column = names(object$column.lower),
                Lower = round(object$column.lower, 4),
                Upper = round(object$column.upper, 4)
            )
            print(head(bounds.df, 10), row.names = FALSE)
            if (length(object$column.lower) > 10) {
                cat("  ... and", length(object$column.lower) - 10, "more columns\n")
            }
        }

    } else {
        dims <- dim(coeffs)

        cat("Array dimensions:\n")
        cat("  Vertices:", dims[1], "\n")
        cat("  Y columns:", dims[2], "\n")
        cat("  Z columns:", dims[3], "\n")
        cat("  Total pairs:", dims[2] * dims[3], "\n\n")

        ## Compute mean matrix
        pair.means <- apply(coeffs, c(2, 3), mean)

        cat("Summary of vertex-averaged coefficients:\n")
        print(summary(as.vector(pair.means)))

        ## Strongest associations
        cat("\nTop 5 (y, z) pairs by |mean coefficient|:\n")
        abs.means <- abs(pair.means)
        top.idx <- order(abs.means, decreasing = TRUE)[1:min(5, length(abs.means))]
        top.coords <- arrayInd(top.idx, dim(pair.means))

        top.df <- data.frame(
            Y = dimnames(coeffs)[[2]][top.coords[, 1]],
            Z = dimnames(coeffs)[[3]][top.coords[, 2]],
            Mean = round(pair.means[top.idx], 4)
        )
        print(top.df, row.names = FALSE)

        if (instrumented && attr(object, "winsorize.quantile") > 0) {
            cat("\nY column winsorization bounds:\n")
            y.bounds <- data.frame(
                Column = names(object$y.column.lower),
                Lower = round(object$y.column.lower, 4),
                Upper = round(object$y.column.upper, 4)
            )
            print(head(y.bounds, 5), row.names = FALSE)

            cat("\nZ column winsorization bounds:\n")
            z.bounds <- data.frame(
                Column = names(object$z.column.lower),
                Lower = round(object$z.column.lower, 4),
                Upper = round(object$z.column.upper, 4)
            )
            print(head(z.bounds, 5), row.names = FALSE)
        }
    }

    invisible(object)
}

#' Vertex-local correlation heatmap for `lcor()` matrix outputs with pair-named columns
#'
#' Designed for `lcor()` results like your `lcor_matrix_matrix_result` where each row is
#' a vertex and each column is a feature-pair named like "A:B" (often storing only the
#' off-diagonal triangle).
#'
#' The function reconstructs the per-vertex p x q matrix (or p x p when symmetric),
#' returns a ComplexHeatmap object, and returns a `map` that links each entry of the
#' row vector back to (i, j) indices of the original X/Y column sets.
#'
#' @param lcor.out Output of `lcor()`: typically a numeric matrix (n.vertices x n.pairs).
#' @param vertex Vertex selector: integer index (1..n) or a character rowname.
#' @param x.names Optional colnames of the left matrix X (used for ordering + indices).
#' @param y.names Optional colnames of the right matrix Y (used for ordering + indices).
#'   If NULL and `symmetric=TRUE`, y.names is set to x.names.
#' @param pair.sep.candidates Character vector of candidate separators for pair names.
#'   For your object, ":" is the right one.
#' @param symmetric Logical or NULL. If NULL, uses attr(lcor.out, "symmetric") when present.
#' @param fill.diagonal Value to put on the diagonal when `symmetric=TRUE`.
#' @param keep.all.features If TRUE and x.names/y.names provided, keep those axes even if
#'   some features never appear in any pair (will yield NA rows/cols). If FALSE, drop
#'   absent features but warn.
#' @param heatmap.name Legend title.
#' @param cluster.rows,cluster.columns Passed to ComplexHeatmap::Heatmap().
#' @param draw If TRUE, draw immediately.
#' @param ... Passed to ComplexHeatmap::Heatmap().
#'
#' @return List with `mat`, `map`, and `heatmap`.
lcor.vertex.heatmap <- function(lcor.out,
                                vertex,
                                x.names = NULL,
                                y.names = NULL,
                                pair.sep.candidates = c(":", "__"),
                                symmetric = NULL,
                                fill.diagonal = 1,
                                keep.all.features = FALSE,
                                heatmap.name = "lcor",
                                cluster.rows = FALSE,
                                cluster.columns = FALSE,
                                draw = TRUE,
                                ...) {
  ## --- deps ---
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required (for colorRamp2).")
  }

  ## --- core checks ---
  if (!is.matrix(lcor.out) || !is.numeric(lcor.out)) {
    stop("`lcor.out` must be a numeric matrix (n.vertices x n.pairs).")
  }

  ## --- symmetric flag ---
  sym.attr <- isTRUE(attr(lcor.out, "symmetric"))
  if (is.null(symmetric)) symmetric <- sym.attr
  symmetric <- isTRUE(symmetric)

  ## --- vertex index ---
  rn <- rownames(lcor.out)
  if (is.character(vertex)) {
    if (is.null(rn)) stop("`vertex` is character but `lcor.out` has no rownames.")
    v.idx <- match(vertex, rn)
    if (is.na(v.idx)) stop("`vertex` not found in `rownames(lcor.out)`.")
  } else {
    v.idx <- as.integer(vertex)
  }
  if (length(v.idx) != 1L || is.na(v.idx) || v.idx < 1L || v.idx > nrow(lcor.out)) {
    stop("`vertex` must select exactly one row (1..n.vertices).")
  }

  ## --- pair colnames parsing ---
  cn <- colnames(lcor.out)
  if (is.null(cn)) stop("`lcor.out` must have colnames encoding feature pairs (e.g. 'A:B').")

  sep.use <- NULL
  parts.use <- NULL
  for (sep.try in pair.sep.candidates) {
    parts.try <- strsplit(cn, sep.try, fixed = TRUE)
    ok <- all(vapply(parts.try, length, integer(1)) == 2L)
    if (isTRUE(ok)) {
      sep.use <- sep.try
      parts.use <- parts.try
      break
    }
  }
  if (is.null(sep.use)) {
    stop("Could not parse pair colnames. Try setting `pair.sep.candidates` to the correct separator.")
  }

  left.name <- vapply(parts.use, function(z) z[[1]], character(1))
  right.name <- vapply(parts.use, function(z) z[[2]], character(1))
  vals <- as.double(lcor.out[v.idx, ])

  ## --- axis names + ordering ---
  if (symmetric && is.null(y.names)) y.names <- x.names

  if (symmetric) {
    present <- unique(c(left.name, right.name))

    ## choose axis order: prefer x.names if provided
    if (!is.null(x.names)) {
      missing <- setdiff(x.names, present)
      if (length(missing) > 0L && !keep.all.features) {
        warning("Dropping features absent from pair colnames: ", paste(missing, collapse = ", "))
      }
      if (keep.all.features) {
        feat <- unique(c(x.names, setdiff(present, x.names)))
      } else {
        feat <- intersect(x.names, present)
        feat <- unique(c(feat, setdiff(present, feat)))
      }
    } else {
      ## stable by appearance
      feat <- unique(c(left.name, right.name))
    }

    p <- length(feat)
    mat <- matrix(NA_real_, nrow = p, ncol = p, dimnames = list(feat, feat))

    i.idx <- match(left.name, feat)
    j.idx <- match(right.name, feat)

    ## fill both triangles
    mat[cbind(i.idx, j.idx)] <- vals
    mat[cbind(j.idx, i.idx)] <- vals
    diag(mat) <- fill.diagonal

    ## mapping back to original columns, if available
    x.col <- if (!is.null(x.names)) match(left.name, x.names) else NA_integer_
    y.col <- if (!is.null(x.names)) match(right.name, x.names) else NA_integer_

    map <- data.frame(
      k = seq_along(vals),
      x.name = left.name,
      y.name = right.name,
      i = i.idx,
      j = j.idx,
      x.col = x.col,
      y.col = y.col,
      value = vals,
      stringsAsFactors = FALSE
    )

  } else {
    ## cross-block case: left are X features, right are Y features
    x.feat.present <- unique(left.name)
    y.feat.present <- unique(right.name)

    if (!is.null(x.names)) {
      missing.x <- setdiff(x.names, x.feat.present)
      if (length(missing.x) > 0L && !keep.all.features) {
        warning("Dropping X features absent from pair colnames: ", paste(missing.x, collapse = ", "))
      }
      x.feat <- if (keep.all.features) unique(c(x.names, setdiff(x.feat.present, x.names))) else {
        tmp <- intersect(x.names, x.feat.present)
        unique(c(tmp, setdiff(x.feat.present, tmp)))
      }
    } else {
      x.feat <- unique(left.name)
    }

    if (!is.null(y.names)) {
      missing.y <- setdiff(y.names, y.feat.present)
      if (length(missing.y) > 0L && !keep.all.features) {
        warning("Dropping Y features absent from pair colnames: ", paste(missing.y, collapse = ", "))
      }
      y.feat <- if (keep.all.features) unique(c(y.names, setdiff(y.feat.present, y.names))) else {
        tmp <- intersect(y.names, y.feat.present)
        unique(c(tmp, setdiff(y.feat.present, tmp)))
      }
    } else {
      y.feat <- unique(right.name)
    }

    p <- length(x.feat)
    q <- length(y.feat)
    mat <- matrix(NA_real_, nrow = p, ncol = q, dimnames = list(x.feat, y.feat))

    i.idx <- match(left.name, x.feat)
    j.idx <- match(right.name, y.feat)
    mat[cbind(i.idx, j.idx)] <- vals

    x.col <- if (!is.null(x.names)) match(left.name, x.names) else NA_integer_
    y.col <- if (!is.null(y.names)) match(right.name, y.names) else NA_integer_

    map <- data.frame(
      k = seq_along(vals),
      x.name = left.name,
      y.name = right.name,
      i = i.idx,
      j = j.idx,
      x.col = x.col,
      y.col = y.col,
      value = vals,
      stringsAsFactors = FALSE
    )
  }

  ## --- heatmap ---
  v.label <- if (!is.null(rn)) rn[v.idx] else as.character(v.idx)
  col.fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

  hm <- ComplexHeatmap::Heatmap(
    mat,
    name = heatmap.name,
    col = col.fun,
    na_col = "grey90",
    cluster_rows = cluster.rows,
    cluster_columns = cluster.columns,
    column_title = paste0("local correlations @ vertex: ", v.label),
    ...
  )

  if (isTRUE(draw)) ComplexHeatmap::draw(hm)

  list(mat = mat, map = map, heatmap = hm)
}

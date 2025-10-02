#' Build Nerve Complex from k-Nearest Neighbors with Riemannian Structure
#'
#' @description
#' Constructs a simplicial complex from the nerve of a k-nearest neighbor
#' covering of a point cloud, equipped with a Riemannian metric derived from
#' neighborhood overlaps. The construction addresses a fundamental question in
#' topological data analysis: how can we extract both the shape (topology) and
#' the local density (geometry) of data through neighborhood intersection patterns?
#'
#' The answer emerges by treating neighborhoods as sets in a measure space. When
#' neighborhoods intersect, they form simplices; the measure of these intersections
#' determines geometric weights. This unified approach produces a structure that
#' supports both topological analysis (homology, persistence) and geometric
#' computation (diffusion, regression).
#'
#' @details
#' \subsection{The Construction Problem}{
#' Given points \eqn{X = \{x_1, \ldots, x_n\} \subset \mathbb{R}^d}, we seek a
#' discrete structure that captures:
#' \enumerate{
#'   \item \strong{Connectivity}: Which points are "close" in a neighborhood sense?
#'   \item \strong{Topology}: What is the shape of the data at different scales?
#'   \item \strong{Geometry}: How does local density vary across the space?
#' }
#'
#' The k-nearest neighbor nerve provides a unified answer. Each point \eqn{x_i}
#' defines a neighborhood \eqn{\hat{N}_k(x_i)} containing itself and its \eqn{k}
#' nearest neighbors. These neighborhoods cover \eqn{X}, and their intersection
#' pattern determines both structure (which simplices exist) and weights (how
#' important they are).
#' }
#'
#' \subsection{From Neighborhoods to Simplices}{
#' A simplex with vertices \eqn{\{i_0, \ldots, i_p\}} exists in the nerve if and
#' only if the corresponding neighborhoods have non-empty intersection:
#' \deqn{\hat{N}_k(x_{i_0}) \cap \cdots \cap \hat{N}_k(x_{i_p}) \neq \emptyset}
#'
#' This is the \emph{nerve theorem} from algebraic topology: the nerve captures
#' the intersection pattern of the covering. Remarkably, under mild conditions
#' (contractible intersections), the nerve is homotopy equivalent to the union
#' of the neighborhoods, meaning they share the same topological features.
#'
#' \strong{Example}: Consider three points with \eqn{k=2}. If all three
#' neighborhoods mutually intersect, we get a triangle (2-simplex). If only
#' pairwise intersections exist, we get three edges forming a triangle boundary.
#' If only one pair intersects, we get a single edge. The nerve records exactly
#' this combinatorial information.
#' }
#'
#' \subsection{From Intersections to Geometry}{
#' The combinatorial nerve lacks metric structure. To perform calculus (compute
#' gradients, diffusion, etc.), we need inner products on chains. These arise
#' naturally by introducing a measure \eqn{\mu} and computing intersection measures:
#'
#' \strong{Vertex mass}: \eqn{m_i = \mu(\hat{N}_k(x_i))} measures the "volume"
#' of the neighborhood
#'
#' \strong{Edge weight}: \eqn{m_{ij} = \mu(\hat{N}_k(x_i) \cap \hat{N}_k(x_j))}
#' measures the overlap
#'
#' \strong{Triangle weight}: \eqn{m_{ijk} = \mu(\hat{N}_k(x_i) \cap \hat{N}_k(x_j) \cap \hat{N}_k(x_k))}
#' measures the three-way overlap
#'
#' These weights form diagonal metric matrices that define inner products on
#' chains at each dimension. The resulting structure is a \emph{Riemannian
#' simplicial complex} supporting discrete differential geometry.
#' }
#'
#' \subsection{Measure Options}{
#' The choice of measure \eqn{\mu} determines how the geometry reflects the data:
#'
#' \strong{Counting Measure} (\code{use_counting_measure = TRUE}):
#' \itemize{
#'   \item Definition: \eqn{\mu(A) = |A|} (cardinality)
#'   \item Interpretation: All points are equally important
#'   \item Use when: Data are approximately uniformly sampled
#'   \item Effect: Weights depend only on overlap size, not density
#' }
#'
#' \strong{Density Measure} (\code{use_counting_measure = FALSE}):
#' \itemize{
#'   \item Definition: \eqn{\mu(A) = \sum_{x \in A} w(x)} with
#'     \eqn{w(x) = (\varepsilon + d_k(x))^{-\alpha}}
#'   \item Interpretation: Points in dense regions receive more weight
#'   \item Use when: Sampling is non-uniform or density information is important
#'   \item Effect: Weights reflect both overlap and local density
#' }
#'
#' The density measure uses k-nearest neighbor distances as robust density
#' surrogates, avoiding the curse of dimensionality that affects kernel methods.
#' }
#'
#' \subsection{Directed vs. Mutual Neighborhoods}{
#' The k-NN relation is asymmetric: \eqn{x_j} being among \eqn{x_i}'s neighbors
#' does not imply the reverse. This creates a choice:
#'
#' \strong{Directed k-NN} (\code{directed_knn = TRUE}):
#' \itemize{
#'   \item Use all directed neighbor relationships
#'   \item Results in denser complex with more simplices
#'   \item May include "one-way" connections at density boundaries
#'   \item Suitable when all local information should be preserved
#' }
#'
#' \strong{Mutual k-NN} (\code{directed_knn = FALSE}, \strong{recommended}):
#' \itemize{
#'   \item Require symmetric neighbor relationships
#'   \item Edge exists only if \eqn{i \in N_k(j)} AND \eqn{j \in N_k(i)}
#'   \item Results in sparser, more stable complex
#'   \item Filters spurious connections at density gradients
#'   \item Default choice for most applications
#' }
#' }
#'
#' \subsection{Choosing k}{
#' The parameter \eqn{k} controls the scale of the construction:
#'
#' \strong{Small k} (e.g., \eqn{k \in [3, 8]}):
#' \itemize{
#'   \item Captures fine-scale local structure
#'   \item May fragment into disconnected components
#'   \item Sensitive to noise and outliers
#'   \item Suitable for well-sampled, low-noise data
#' }
#'
#' \strong{Moderate k} (e.g., \eqn{k \in [10, 20]}):
#' \itemize{
#'   \item Balances local and global structure
#'   \item Usually produces connected complex
#'   \item Robust to moderate noise
#'   \item Recommended starting point
#' }
#'
#' \strong{Large k} (e.g., \eqn{k > 30}):
#' \itemize{
#'   \item Emphasizes global structure
#'   \item May smooth over important features
#'   \item Very dense complex with many simplices
#'   \item Computational cost increases
#' }
#'
#' \strong{Rule of thumb}: Start with \eqn{k \approx \log(n)} and adjust based
#' on the complexity of features you wish to capture. Examine the resulting
#' complex using \code{\link{riem.dcx.summary}} to check connectivity.
#' }
#'
#' \subsection{Implementation Notes}{
#' \itemize{
#'   \item k-NN computation uses the ANN library with kd-trees for efficiency
#'   \item Complexity: \eqn{O(n \log n)} for k-NN, \eqn{O(n k^2)} for edge construction
#'   \item Memory scales as \eqn{O(nk + m)} where \eqn{m} is the number of edges
#'   \item For \code{max_p > 2}, expect combinatorial growth in simplex count
#'   \item The function automatically reduces \code{max_p} if higher simplices don't exist
#'   \item Sparse input matrices (dgCMatrix) are supported and recommended for high dimensions
#' }
#' }
#'
#' @param X Numeric matrix of dimension \eqn{n \times d} where rows are points
#'   in \eqn{d}-dimensional space. Can be dense matrix or sparse matrix (dgCMatrix
#'   from Matrix package). Features should be appropriately scaled.
#'
#' @param y Numeric vector of length \eqn{n} giving response values, or \code{NULL}.
#'   If provided, enables subsequent signal processing operations (smoothing,
#'   regression). If \code{NULL}, constructs topology-only complex.
#'
#' @param k Integer number of nearest neighbors, satisfying \eqn{2 \le k < n}.
#'   Controls the scale of neighborhood coverage. Default is 10.
#'
#' @param max.p Integer maximum simplex dimension to construct, satisfying
#'   \eqn{1 \le \text{max.p} < n}. The complex will include simplices from
#'   dimension 0 (vertices) through \code{max.p}. Default is 2 (include triangles).
#'
#' @param use.counting.measure Logical. If \code{TRUE} (default), use counting
#'   measure giving equal weight to all points. If \code{FALSE}, use density-based
#'   weights derived from k-NN distances.
#'
#' @param directed.knn Logical. If \code{TRUE}, use directed k-NN (asymmetric).
#'   If \code{FALSE} (default, recommended), enforce mutual k-NN (symmetric).
#'
#' @param density.normalization Non-negative numeric. Target sum for normalized
#'   weights. If 0 (default), normalize to sum to \eqn{n}. Ignored when
#'   \code{use.counting.measure = TRUE}.
#'
#' @return An S4 object of class \code{"Rcpp_riem_dcx"} (external pointer) to
#'   the constructed complex. This object contains:
#'   \itemize{
#'     \item Simplex tables recording vertices of each simplex
#'     \item Metric matrices encoding geometry at each dimension
#'     \item Boundary operators for discrete differential calculus
#'     \item Hodge Laplacians for diffusion and spectral analysis
#'     \item Density distributions and response data (if provided)
#'   }
#'
#'   Use \code{\link{riem.dcx.summary}} to inspect the structure, and other
#'   package functions for analysis (heat diffusion, spectral decomposition, etc.).
#'
#' @examples
#' # Example 1: Two-dimensional noisy circle
#' # Demonstrates topology recovery from noisy samples
#' set.seed(42)
#' theta <- seq(0, 2*pi, length.out = 100)
#' X_circle <- cbind(cos(theta) + rnorm(100, sd=0.1),
#'                   sin(theta) + rnorm(100, sd=0.1))
#' y_circle <- sin(3*theta) + rnorm(100, sd=0.2)
#'
#' # Build nerve with moderate k
#' dcx_circle <- build.nerve.from.knn(
#'   X = X_circle,
#'   y = y_circle,
#'   k = 8,
#'   max.p = 2,
#'   use.counting.measure = TRUE
#' )
#'
#' # Inspect structure
#' summary(dcx_circle)
#' # Should see: ~100 vertices, ~800 edges (dense graph), many triangles
#'
#' \dontrun{
#' # Example 2: High-dimensional data with density weights
#' # Demonstrates density-aware construction
#' set.seed(123)
#' n <- 200
#' d <- 50
#'
#' # Generate data with varying density
#' X_hd <- matrix(rnorm(n * d), nrow=n, ncol=d)
#' # Add dense cluster
#' cluster_idx <- 1:50
#' X_hd[cluster_idx, 1:5] <- X_hd[cluster_idx, 1:5] * 0.3  # Shrink variance
#'
#' y_hd <- rowSums(X_hd[, 1:3]) + rnorm(n)
#'
#' dcx_hd <- build.nerve.from.knn(
#'   X = X_hd,
#'   y = y_hd,
#'   k = 15,
#'   max.p = 2,
#'   use.counting.measure = FALSE,  # Density-aware
#'   directed.knn = FALSE
#' )
#'
#' # Vertices in the dense cluster receive higher weights
#'
#' # Example 3: Sparse matrix for very high dimensions
#' # Demonstrates memory-efficient computation
#' library(Matrix)
#' n <- 500
#' d <- 1000
#'
#' # Sparse random projection
#' X_sparse <- rsparsematrix(n, d, density=0.1)
#' y_sparse <- rowSums(as.matrix(X_sparse[, 1:10])) + rnorm(n)
#'
#' dcx_sparse <- build.nerve.from.knn(
#'   X = X_sparse,  # Automatically handled
#'   y = y_sparse,
#'   k = 20,
#'   max.p = 1,     # Graph only (edges, no triangles)
#'   use.counting.measure = TRUE
#' )
#'
#' # Example 4: Topology-only (no response)
#' # Demonstrates pure topological analysis
#' # Create two noisy circles (annulus)
#' n_outer <- 100
#' n_inner <- 60
#' theta_outer <- seq(0, 2*pi, length.out=n_outer)
#' theta_inner <- seq(0, 2*pi, length.out=n_inner)
#'
#' X_annulus <- rbind(
#'   cbind(cos(theta_outer) + rnorm(n_outer, sd=0.05),
#'         sin(theta_outer) + rnorm(n_outer, sd=0.05)),
#'   cbind(0.5*cos(theta_inner) + rnorm(n_inner, sd=0.05),
#'         0.5*sin(theta_inner) + rnorm(n_inner, sd=0.05))
#' )
#'
#' dcx_annulus <- build.nerve.from.knn(
#'   X = X_annulus,
#'   y = NULL,      # No response
#'   k = 6,
#'   max.p = 2,
#'   use.counting.measure = TRUE
#' )
#'
#' # Can now compute homology to detect two circles
#' # (First Betti number should be 2)
#' }
#'
#' @references
#' Carlsson, G. (2009). Topology and data. \emph{Bulletin of the American
#' Mathematical Society}, 46(2), 255-308.
#'
#' Edelsbrunner, H., & Harer, J. (2010). \emph{Computational Topology:
#' An Introduction}. American Mathematical Society.
#'
#' Lim, L. H. (2020). Hodge Laplacians on graphs. \emph{SIAM Review},
#' 62(3), 685-715.
#'
#' @seealso
#' \code{\link{riem.dcx.summary}} for extracting structural information,
#' \code{\link{riem.dcx.iteration.step}} for density evolution and signal fitting,
#' \code{\link{riem.dcx.heat.diffusion}} for heat-based smoothing
#'
#' @export
build.nerve.from.knn <- function(
                                 X,
                                 y = NULL,
                                 k = 10,
                                 max.p = 2,
                                 use.counting.measure = TRUE,
                                 directed.knn = FALSE,
                                 density.normalization = 0.0
                                 ) {
    ## Validate inputs at R level
    if (!is.matrix(X) && !inherits(X, "Matrix")) {
        stop("X must be a numeric matrix or sparse Matrix (dgCMatrix)")
    }

    if (!is.null(y) && !is.numeric(y)) {
        stop("y must be a numeric vector or NULL")
    }

    if (!is.numeric(k) || length(k) != 1 || k != as.integer(k)) {
        stop("k must be a single integer")
    }

    if (!is.numeric(max.p) || length(max.p) != 1 || max.p != as.integer(max.p)) {
        stop("max.p must be a single integer")
    }

    if (!is.logical(use.counting.measure) || length(use.counting.measure) != 1) {
        stop("use.counting.measure must be a single logical value")
    }

    if (!is.logical(directed.knn) || length(directed.knn) != 1) {
        stop("directed.knn must be a single logical value")
    }

    if (!is.numeric(density.normalization) || length(density.normalization) != 1) {
        stop("density.normalization must be a single numeric value")
    }

    storage.mode(X) <- "double"

    ## Call C++ implementation via .Call()
    dcx <- .Call(
        "S_build_nerve_from_knn",
        X,
        as.double(y),
        as.integer(k + 1), # ANN returns self as a member of k-nearest neighbors, so we need to adjust for it adding 1
        as.integer(max.p),
        as.logical(use.counting.measure),
        as.logical(directed.knn),
        as.numeric(density.normalization)
    )

    ## Attach construction parameters as attributes
    attr(dcx, "construction_params") <- list(
        n_points = nrow(X),
        n_features = ncol(X),
        k = k,
        max_p = max.p,
        use_counting_measure = use.counting.measure,
        directed_knn = directed.knn,
        density_normalization = density.normalization,
        has_response = !is.null(y)
    )

    class(dcx) <- c("riem_dcx", "list")

    return(dcx)
}

#' Extract Summary Information from Riemannian Complex
#'
#' @description
#' Extracts structural and geometric information from a constructed nerve
#' complex for inspection, diagnostics, and validation.
#'
#' @param dcx A Riemannian discrete chain complex object created by
#'   \code{\link{build.nerve.from.knn}}
#'
#' @return A list with components:
#' \describe{
#'   \item{max_dimension}{Maximum dimension of simplices}
#'   \item{n_simplices}{Integer vector of simplex counts by dimension}
#'   \item{n_vertices}{Number of 0-simplices (vertices)}
#'   \item{n_edges}{Number of 1-simplices (edges)}
#'   \item{n_triangles}{Number of 2-simplices (triangles)}
#'   \item{has_response}{Whether a response vector is stored}
#'   \item{response_length}{Length of response (0 if none)}
#'   \item{construction_params}{Parameters used to construct the complex}
#' }
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 5), ncol=5)
#' y <- rowSums(X[,1:2]) + rnorm(100)
#' dcx <- build.nerve.from.knn(X, y, k=10, max.p=2)
#'
#' info <- riem.dcx.summary(dcx)
#' print(info)
#' }
#'
#' @export
riem.dcx.summary <- function(dcx) {
  if (!inherits(dcx, "riem_dcx")) {
    stop("dcx must be a riem_dcx object created by build.nerve.from.knn")
  }

  # Call C++ summary via .Call()
  cpp_info <- .Call(S_riem_dcx_summary, dcx, PACKAGE = "gflow")

  # Add construction parameters if available
  params <- attr(dcx, "construction_params")
  if (!is.null(params)) {
    cpp_info$construction_params <- params
  }

  # Set class for pretty printing
  class(cpp_info) <- c("riem_dcx_summary", "list")

  return(cpp_info)
}

#' @export
print.riem_dcx_summary <- function(x, ...) {
  cat("Riemannian Discrete Chain Complex Summary\n")
  cat("==========================================\n\n")

  cat("Structure:\n")
  cat(sprintf("  Maximum dimension: %d\n", x$max_dimension))
  cat(sprintf("  Vertices (0-simplices): %d\n", x$n_vertices))
  cat(sprintf("  Edges (1-simplices): %d\n", x$n_edges))
  if (x$max_dimension >= 2) {
    cat(sprintf("  Triangles (2-simplices): %d\n", x$n_triangles))
  }

  if (x$max_dimension >= 3) {
    cat(sprintf("  Higher simplices:\n"))
    for (p in 3:x$max_dimension) {
      cat(sprintf("    Dimension %d: %d\n", p, x$n_simplices[p + 1]))
    }
  }

  cat("\nResponse:\n")
  if (x$has_response) {
    cat(sprintf("  Response vector present (length %d)\n", x$response_length))
  } else {
    cat("  No response vector (topology-only complex)\n")
  }

  if (!is.null(x$construction_params)) {
    cat("\nConstruction Parameters:\n")
    params <- x$construction_params
    cat(sprintf("  k (nearest neighbors): %d\n", params$k))
    cat(sprintf("  Requested max dimension: %d\n", params$max_p))
    cat(sprintf("  Measure type: %s\n",
                if (params$use_counting_measure) "counting" else "density-based"))
    cat(sprintf("  Neighborhood type: %s k-NN\n",
                if (params$directed_knn) "directed" else "mutual"))
  }

  invisible(x)
}

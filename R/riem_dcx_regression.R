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
        as.numeric(density.normalization)
    )

    ## Attach construction parameters as attributes
    attr(dcx, "construction_params") <- list(
        n_points = nrow(X),
        n_features = ncol(X),
        k = k,
        max_p = max.p,
        use_counting_measure = use.counting.measure,
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
  }

  invisible(x)
}

#' Fit k-Nearest Neighbor Riemannian Graph Regression
#'
#' Constructs a 1-skeleton k-nearest neighbor complex and iteratively adapts
#' its Riemannian geometry to reflect response structure, enabling conditional
#' expectation estimation that respects both feature topology and response
#' boundaries.
#'
#' @param X Numeric matrix (n × d) of features, or sparse matrix (dgCMatrix).
#'   Each row represents an observation, each column a feature. Features
#'   should be appropriately scaled if they have different units.
#'
#' @param y Numeric vector of length n containing response values. Cannot
#'   contain NA or infinite values.
#'
#' @param k Integer scalar giving the number of nearest neighbors. Must satisfy
#'   2 ≤ k < n. Larger k produces smoother fits but may oversmooth fine-scale
#'   features. Typical values are in [5, 30]. If NULL, selects k via
#'   cross-validation (NOT YET IMPLEMENTED).
#'
#' @param use.counting.measure Logical scalar. If TRUE (default), uses uniform
#'   vertex weights (counting measure). If FALSE, uses distance-based weights
#'   inversely proportional to local k-NN density: w(x) = (ε + d_k(x))^(-α).
#'   Distance-based weights are useful when sampling density varies across
#'   the feature space.
#'
#' @param density.normalization Numeric scalar, non-negative. Specifies target
#'   sum for normalized vertex densities. If 0 (default), densities are
#'   normalized to sum to n (giving average vertex density = 1). If positive,
#'   densities sum to the specified value. Most users should use the default.
#'   Set to 1.0 for probability interpretation, or to a large value (e.g., 1e6)
#'   for very large graphs to avoid numerical underflow.
#'
#' @param t.diffusion Diffusion time parameter (non-negative). Controls how much
#'   vertex density evolves in each iteration. If 0 (default), automatically
#'   selected as 0.5/λ₂ where λ₂ is the spectral gap. Larger values produce
#'   more aggressive density updates; smaller values are more conservative.
#'   Typical range: [0.1/λ₂, 1.0/λ₂]. Set to 0 for automatic selection (recommended).
#'
#' @param beta.damping Damping parameter (non-negative). Controls the strength
#'   of the restoring force toward uniform distribution. If 0 (default),
#'   automatically selected as 0.1/t.diffusion to maintain a 10:1 ratio of
#'   diffusion to damping. Larger values prevent concentration but may over-smooth;
#'   smaller values allow more geometric structure but risk collapse. Set to 0
#'   for automatic selection (recommended).
#'
#' @param gamma.modulation Numeric scalar in [0.5, 2.0]. Exponent controlling
#'   sharpness of response boundaries in the geometry. Default 1.0 corresponds
#'   to Cauchy kernel: Γ(Δ) = (1 + Δ²/σ²)^(-γ). Larger γ creates sharper
#'   boundaries (more aggressive edge pruning at response discontinuities).
#'   Smaller γ creates gentler transitions. Values outside [0.5, 2] trigger
#'   a warning.
#'
#' @param n.eigenpairs Integer scalar in [10, n]. Number of Laplacian
#'   eigenpairs to compute for spectral filtering. Default 200. Larger values
#'   capture more frequency components but increase computation time. Must be
#'   less than the number of observations n.
#'
#' @param filter.type Character scalar. Spectral filter type for response
#'   smoothing. One of:
#'   \itemize{
#'     \item \code{"heat_kernel"} (default): f(λ) = exp(-ηλ), corresponding
#'           to solving the heat equation. Provides exponential decay of
#'           high frequencies.
#'     \item \code{"tikhonov"}: f(λ) = 1/(1 + ηλ), corresponding to Tikhonov
#'           regularization. Gentler high-frequency attenuation than heat kernel.
#'     \item \code{"cubic_spline"}: f(λ) = 1/(1 + ηλ²), minimizing second
#'           derivatives. Produces spline-like smoothness.
#'   }
#'
#' @param epsilon.y Numeric scalar, positive. Relative convergence threshold
#'   for response. Iteration stops when the relative change in fitted values
#'   ||ŷ^(ℓ) - ŷ^(ℓ-1)|| / ||ŷ^(ℓ-1)|| falls below this value. Default 1e-4.
#'   Smaller values require more iterations but give more precise convergence.
#'
#' @param epsilon.rho Numeric scalar, positive. Relative convergence threshold
#'   for densities. Iteration stops when max_p ||ρ_p^(ℓ) - ρ_p^(ℓ-1)|| / ||ρ_p^(ℓ-1)||
#'   falls below this value across all dimensions. Default 1e-4.
#'
#' @param max.iterations Integer scalar, positive. Maximum number of iterations.
#'   Default 50. Typical convergence occurs within 5-10 iterations for
#'   well-chosen parameters. If maximum iterations reached without convergence,
#'   a warning is issued.
#'
#' @param max.ratio.threshold Numeric in \eqn{[0, 0.2)}.
#'     Geometric pruning removes an edge \eqn{(i,j)} when there exists an
#'     alternative path between \eqn{i} and \eqn{j} whose path/edge length ratio
#'     minus 1.0 is \emph{less than} this threshold. This is a deviation
#'     threshold \eqn{\delta} in \eqn{[0, 0.2)}. Internally we compare the
#'     path-to-edge ratio R to \eqn{1 + \delta}. Default 0.1.
#'
#' @param threshold.percentile Numeric in \eqn{[0,1]}. Only edges with
#'     length above this percentile are considered for geometric pruning. Default 0.5.
#'
#' @return An object of class \code{c("knn.riem.fit", "riem.dcx")}, which is
#'   an external pointer to a fitted riem_dcx_t object. Use accessor functions
#'   to extract components:
#'   \itemize{
#'     \item \code{\link{fitted.knn.riem.fit}}: Extract fitted values
#'     \item \code{\link{residuals.knn.riem.fit}}: Extract residuals
#'     \item \code{\link{summary.knn.riem.fit}}: Print model summary
#'     \item \code{\link{predict.knn.riem.fit}}: Predict on new data (future)
#'   }
#'
#' @details
#' \strong{Algorithm Overview}
#'
#' The algorithm proceeds in two main phases:
#'
#' \emph{Phase I: Initialization}
#' \enumerate{
#'   \item Constructs 1-skeleton k-NN complex from features using ANN library
#'   \item Computes initial vertex and edge densities from neighborhood overlaps
#'   \item Builds mass matrices M₀ (diagonal) and M₁ (full, via Gram determinants)
#'   \item Assembles Hodge Laplacian L₀ = B₁ M₁⁻¹ B₁ᵀ M₀
#'   \item Performs initial spectral filtering: ŷ⁽⁰⁾ = V F_η(Λ) Vᵀ y
#' }
#'
#' \emph{Phase II: Iterative Refinement}
#'
#' Repeats until convergence or max.iterations reached:
#' \enumerate{
#'   \item \strong{Density diffusion}: Evolve vertex density via damped heat
#'         equation (I + t(L₀ + βI))ρ₀ = ρ₀_prev + tβu
#'   \item \strong{Edge density update}: Recompute edge densities from evolved
#'         vertex masses via neighborhood intersections
#'   \item \strong{Response-coherence modulation}: Reduce edge masses where
#'         response varies rapidly: ρ₁([i,j]) ← ρ₁([i,j]) · (1 + Δ²ᵢⱼ/σ²)⁻ᵞ
#'   \item \strong{Metric update}: Rebuild M₀, M₁ from modified densities
#'   \item \strong{Laplacian reassembly}: Recompute L₀ with updated metric
#'   \item \strong{Response smoothing}: Re-filter response using new geometry
#'   \item \strong{Convergence check}: Test if response and densities stabilized
#' }
#'
#' \strong{Geometric Interpretation}
#'
#' The method creates a Riemannian metric on the k-NN graph where:
#' \itemize{
#'   \item Edges within response-coherent regions have high mass (short distance)
#'   \item Edges crossing response boundaries have low mass (long distance)
#'   \item Diffusion is rapid within coherent regions, slow across boundaries
#'   \item Spectral filtering respects this geometry, producing smooth fits
#'         within regions while preserving discontinuities at boundaries
#' }
#'
#' \strong{Parameter Selection Guidelines}
#'
#' For most applications, default parameters work well. Consider tuning:
#' \itemize{
#'   \item \strong{k}: Primary tuning parameter. Start with k = ⌈log(n)⌉ and
#'         adjust based on cross-validation. Larger k smooths more.
#'   \item \strong{gamma.modulation}: Increase (up to 2.0) for sharper response
#'         boundaries; decrease (down to 0.5) for gentler transitions.
#'   \item \strong{filter.type}: "heat_kernel" (default) for general use;
#'         "tikhonov" for better linear trend preservation; "cubic_spline"
#'         for smoother results on chain-like graphs.
#' }
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item k-NN search: O(n log n) using kd-trees (ANN library)
#'   \item Edge construction: O(n k²) for intersection computations
#'   \item Per iteration: O(n k² + m n_eig) where m ≈ n k is edge count
#'   \item Total: O(n k² + iter × m n_eig), typically manageable for n ≤ 10,000
#' }
#'
#' @section Memory Usage:
#' \itemize{
#'   \item Feature matrix: O(n d)
#'   \item Complex storage: O(n k) for vertices and edges
#'   \item Mass matrices: O(n + m²) but M₁ is sparse with O(m k²) nonzeros
#'   \item Eigenpairs: O(n × n_eigenpairs)
#'   \item Peak usage: dominated by Laplacian operations, typically < 1GB for n = 5000
#' }
#'
#' @references
#' Van der Waerden, B. L. (1975). \emph{Algebra}. Springer-Verlag.
#' (Writing style inspiration)
#'
#' Lim, L. H. (2020). Hodge Laplacians on graphs. \emph{SIAM Review}, 62(3), 685-715.
#' (Hodge theory on simplicial complexes)
#'
#' Belkin, M., & Niyogi, P. (2003). Laplacian eigenmaps for dimensionality
#' reduction and data representation. \emph{Neural Computation}, 15(6), 1373-1396.
#' (Graph-based geometric methods)
#'
#' @examples
#' \dontrun{
#' # Example 1: Simple nonlinear regression
#' set.seed(123)
#' n <- 200
#' X <- matrix(rnorm(n * 3), n, 3)
#' y <- sin(2 * pi * X[,1]) + 0.5 * X[,2]^2 + rnorm(n, sd = 0.2)
#'
#' fit <- fit.knn.riem.graph.regression(X, y, k = 15)
#' y.hat <- fitted(fit)
#' plot(y, y.hat, asp = 1, main = "Fitted vs Actual")
#' abline(0, 1, col = "red")
#'
#' # Example 2: Regression with sharp boundaries
#' X <- matrix(runif(n * 2, -1, 1), n, 2)
#' y <- ifelse(X[,1]^2 + X[,2]^2 < 0.5, 1, 0) + rnorm(n, sd = 0.1)
#'
#' # Use larger gamma for sharper boundary detection
#' fit <- fit.knn.riem.graph.regression(X, y, k = 20, gamma.modulation = 1.5)
#' y.hat <- fitted(fit)
#'
#' # Example 3: Using sparse matrix input
#' library(Matrix)
#' X.sparse <- Matrix(X, sparse = TRUE)
#' fit <- fit.knn.riem.graph.regression(X.sparse, y, k = 15)
#' }
#'
#' @seealso
#' \code{\link{fitted.knn.riem.fit}} for extracting fitted values
#' \code{\link{residuals.knn.riem.fit}} for extracting residuals
#' \code{\link{summary.knn.riem.fit}} for model summaries
#'
#' @export
fit.knn.riem.graph.regression <- function(
    X,
    y,
    k = NULL,
    use.counting.measure = TRUE,
    density.normalization = 0,
    t.diffusion = 0,
    beta.damping = 0,
    gamma.modulation = 1.0,
    n.eigenpairs = 10,
    filter.type = c("heat_kernel", "tikhonov", "cubic_spline"),
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    max.iterations = 50,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.5,
    test.stage
) {
    # ==================== Feature Matrix Validation ====================

    # Check X type
    is.dense <- is.matrix(X)
    is.sparse <- inherits(X, "dgCMatrix")

    if (!is.dense && !is.sparse) {
        stop("X must be a numeric matrix or sparse matrix (dgCMatrix from Matrix package)")
    }

    # For dense matrix, check numeric
    if (is.dense && !is.numeric(X)) {
        stop("X must be numeric")
    }

    # Extract dimensions
    n <- nrow(X)
    d <- ncol(X)

    # Check dimensions
    if (is.null(n) || is.null(d) || n < 1 || d < 1) {
        stop("X must have valid dimensions (positive rows and columns)")
    }

    if (n < 10) {
        stop(sprintf("X must have at least 10 observations (got n=%d)", n))
    }

    # Check for NA/Inf in X
    if (is.dense) {
        if (anyNA(X)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X))) {
            stop("X cannot contain infinite values")
        }
        # Ensure storage mode is double
        storage.mode(X) <- "double"
    } else {
        # For sparse matrix, check x slot
        if (anyNA(X@x)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X@x))) {
            stop("X cannot contain infinite values")
        }
    }

    # ==================== Response Vector Validation ====================

    # Check y type and length
    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    if (length(y) != n) {
        stop(sprintf("Length of y (%d) must equal number of rows in X (%d)",
                     length(y), n))
    }

    # Check for NA/Inf in y
    if (anyNA(y)) {
        stop("y cannot contain NA values")
    }

    if (any(!is.finite(y))) {
        stop("y cannot contain infinite values")
    }

    # Ensure y is a vector (not matrix)
    y <- as.vector(y)

    # ==================== Parameter k Validation ====================

    if (is.null(k)) {
        stop("Automatic k selection not yet implemented.\n",
             "Please specify k manually. Recommended starting value: k = ",
             ceiling(log(n)), " (roughly log(n))")
    }

    if (!is.numeric(k) || length(k) != 1) {
        stop("k must be a single numeric value")
    }

    # Convert to integer
    k <- as.integer(k)

    # Check range
    if (k < 2) {
        stop(sprintf("k must be at least 2 (got k=%d)", k))
    }

    if (k >= n) {
        stop(sprintf("k must be less than n (got k=%d, n=%d)", k, n))
    }

    # ==================== Logical Parameter Validation ====================

    if (!is.logical(use.counting.measure) || length(use.counting.measure) != 1) {
        stop("use.counting.measure must be TRUE or FALSE")
    }

    if (is.na(use.counting.measure)) {
        stop("use.counting.measure cannot be NA")
    }

    # ==================== Numeric Parameter Validation ====================

    # density.normalization
    if (!is.numeric(density.normalization) || length(density.normalization) != 1) {
        stop("density.normalization must be a single numeric value")
    }

    if (is.na(density.normalization) || !is.finite(density.normalization)) {
        stop("density.normalization cannot be NA or infinite")
    }

    if (density.normalization < 0) {
        stop(sprintf("density.normalization must be non-negative (got %.3f)",
                     density.normalization))
    }

    # t.diffusion
    if (!is.numeric(t.diffusion) || length(t.diffusion) != 1) {
        stop("t.diffusion must be a single numeric value")
    }

    if (is.na(t.diffusion) || !is.finite(t.diffusion)) {
        stop("t.diffusion cannot be NA or infinite")
    }

    if (t.diffusion < 0) {
        stop(sprintf("t.diffusion must be non-negative (got %.3f; use 0 for auto)",
                     t.diffusion))
    }

    # beta.damping
    if (!is.numeric(beta.damping) || length(beta.damping) != 1) {
        stop("beta.damping must be a single numeric value")
    }

    if (is.na(beta.damping) || !is.finite(beta.damping)) {
        stop("beta.damping cannot be NA or infinite")
    }

    if (beta.damping < 0) {
        stop(sprintf("beta.damping must be non-negative (got %.3f; use 0 for auto)",
                     beta.damping))
    }

    # gamma.modulation
    if (!is.numeric(gamma.modulation) || length(gamma.modulation) != 1) {
        stop("gamma.modulation must be a single numeric value")
    }

    if (is.na(gamma.modulation) || !is.finite(gamma.modulation)) {
        stop("gamma.modulation cannot be NA or infinite")
    }

    if (gamma.modulation <= 0) {
        stop(sprintf("gamma.modulation must be positive (got %.3f)", gamma.modulation))
    }

    if (gamma.modulation < 0.5 || gamma.modulation > 2.0) {
        warning(sprintf(paste0("gamma.modulation=%.3f is outside recommended range [0.5, 2.0].\n",
                              "Values outside this range may produce unstable geometry."),
                       gamma.modulation))
    }

    # n.eigenpairs
    if (!is.numeric(n.eigenpairs) || length(n.eigenpairs) != 1) {
        stop("n.eigenpairs must be a single numeric value")
    }

    n.eigenpairs <- as.integer(n.eigenpairs)

    if (n.eigenpairs < 10) {
        stop(sprintf("n.eigenpairs must be at least 10 (got %d)", n.eigenpairs))
    }

    if (n.eigenpairs > n) {
        stop(sprintf("n.eigenpairs cannot exceed n (got %d, n=%d)", n.eigenpairs, n))
    }

    # Warn if n.eigenpairs seems too small
    if (n.eigenpairs < 50 && n > 100) {
        warning(sprintf(paste0("n.eigenpairs=%d may be too small for n=%d observations.\n",
                              "Consider n.eigenpairs >= 100 for better frequency resolution."),
                       n.eigenpairs, n))
    }

    # epsilon.y
    if (!is.numeric(epsilon.y) || length(epsilon.y) != 1) {
        stop("epsilon.y must be a single numeric value")
    }

    if (is.na(epsilon.y) || !is.finite(epsilon.y)) {
        stop("epsilon.y cannot be NA or infinite")
    }

    if (epsilon.y <= 0) {
        stop(sprintf("epsilon.y must be positive (got %.3e)", epsilon.y))
    }

    if (epsilon.y >= 0.1) {
        warning(sprintf("epsilon.y=%.3f is very large; convergence may be premature", epsilon.y))
    }

    # epsilon.rho
    if (!is.numeric(epsilon.rho) || length(epsilon.rho) != 1) {
        stop("epsilon.rho must be a single numeric value")
    }

    if (is.na(epsilon.rho) || !is.finite(epsilon.rho)) {
        stop("epsilon.rho cannot be NA or infinite")
    }

    if (epsilon.rho <= 0) {
        stop(sprintf("epsilon.rho must be positive (got %.3e)", epsilon.rho))
    }

    if (epsilon.rho >= 0.1) {
        warning(sprintf("epsilon.rho=%.3f is very large; convergence may be premature", epsilon.rho))
    }

    # max.iterations
    if (!is.numeric(max.iterations) || length(max.iterations) != 1) {
        stop("max.iterations must be a single numeric value")
    }

    max.iterations <- as.integer(max.iterations)

    if (max.iterations < 1) {
        stop(sprintf("max.iterations must be at least 1 (got %d)", max.iterations))
    }

    if (max.iterations > 1000) {
        warning(sprintf(paste0("max.iterations=%d is very large.\n",
                              "Typical convergence occurs within 5-10 iterations."),
                       max.iterations))
    }

    # max.ratio.threshold
    if (!is.numeric(max.ratio.threshold) || length(max.ratio.threshold) != 1) {
        stop("max.ratio.threshold must be numeric.")
    }
        
    if (max.ratio.threshold < 0 || max.ratio.threshold >= 0.2) {
        stop("max.ratio.threshold must be in [0, 0.2).")
    }
        
    # threshold.percentile
    if (!is.numeric(threshold.percentile) || length(threshold.percentile) != 1 ||
        threshold.percentile < 0 || threshold.percentile > 1) {
        stop("threshold.percentile must be in [0, 1].")
    }
        
    # ==================== String Parameter Validation ====================

    # Match filter.type
    filter.type <- match.arg(filter.type)

    # ==================== Call C++ Function ====================

    fit <- .Call(
        S_fit_knn_riem_graph_regression,
        X,
        as.double(y),
        as.integer(k),
        as.logical(use.counting.measure),
        as.double(density.normalization),
        as.double(t.diffusion),
        as.double(beta.damping),
        as.double(gamma.modulation),
        as.integer(n.eigenpairs),
        as.character(filter.type),
        as.double(epsilon.y),
        as.double(epsilon.rho),
        as.integer(max.iterations),
        as.double(max.ratio.threshold),
        as.double(threshold.percentile),
        as.integer(test.stage),
        PACKAGE = "gflow"
    )

    # ==================== Post-Processing ====================

    # Add class
    class(fit) <- c("knn.riem.fit", "riem.dcx")

    # Store call for reproducibility
    attr(fit, "call") <- match.call()

    # Store key parameters for summary()
    attr(fit, "n") <- n
    attr(fit, "d") <- d
    attr(fit, "k") <- k

    return(fit)
}

#' Extract Fitted Values from kNN Riemannian Regression Fit
#'
#' @param object Object of class \code{knn.riem.fit}
#' @param ... Additional arguments (currently unused)
#' @return Numeric vector of fitted values (length n)
#' @export
fitted.knn.riem.fit <- function(object, ...) {
    object$fitted.values
}

#' @export
residuals.knn.riem.fit <- function(object, ...) {
    object$residuals
}

#' @export
print.knn.riem.fit <- function(x, ...) {
    cat("k-NN Riemannian Graph Regression Fit\n\n")
    cat(sprintf("Data: n = %d observations\n", length(x$y)))
    cat(sprintf("k-NN parameter: k = %d\n\n", x$parameters$k))
    cat(sprintf("Convergence: %s\n", if(x$iteration$converged) "YES" else "NO"))
    cat(sprintf("Iterations: %d\n", x$iteration$n.iterations))
    invisible(x)
}

#' @export
summary.knn.riem.fit <- function(object, ...) {
    # Can still create rich summary without additional C++ calls
    # All data is in the list

    resid <- object$residuals

    structure(
        list(
            call = attr(object, "call"),
            n = length(object$y),
            k = object$parameters$k,
            converged = object$iteration$converged,
            n_iterations = object$iteration$n.iterations,
            residuals = resid,
            sigma = sd(resid),
            r.squared = 1 - var(resid) / var(object$y)
        ),
        class = "summary.knn.riem.fit"
    )
}

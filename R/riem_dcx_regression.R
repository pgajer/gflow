#' Fit k-Nearest Neighbor Riemannian Graph Regression
#'
#' Constructs a 1-skeleton k-nearest neighbor complex and iteratively adapts
#' its Riemannian geometry to reflect response structure, enabling conditional
#' expectation estimation that respects both feature topology and response
#' boundaries.
#'
#' @param X Numeric matrix (n x d) of features, or sparse matrix (dgCMatrix).
#'   Each row represents an observation, each column a feature. Features
#'   should be appropriately scaled if they have different units.
#'
#' @param y Numeric vector of length n containing response values. Cannot
#'   contain NA or infinite values.
#'
#' @param k Integer scalar giving the number of nearest neighbors. Must satisfy
#'   \eqn{2 \le k < n}. Larger k produces smoother fits but may oversmooth
#'   fine-scale features. Typical values are in the range 5 to 30. If NULL,
#'   selects k via cross-validation (NOT YET IMPLEMENTED).
#'
#' @param pca.dim Positive integer or `NULL`. If not `NULL` and `ncol(X) >
#'     pca.dim`, PCA is used to reduce to at most `pca.dim` components.
#'
#' @param variance.explained Numeric in \eqn{(0,1]} or `NULL`. If not `NULL`,
#'     choose the smallest number of PCs whose cumulative variance explained
#'     exceeds this threshold, capped by `pca.dim`.
#'
#' @param max.iterations Integer scalar, positive. Maximum number of iterations.
#'   Default 50. Typical convergence occurs within 5-10 iterations for
#'   well-chosen parameters. If maximum iterations reached without convergence,
#'   a warning is issued.
#'
#' @param n.eigenpairs Integer scalar in the range 10 to n. Number of Laplacian
#'   eigenpairs to compute for spectral filtering. Default 200. Larger values
#'   capture more frequency components but increase computation time. Must be
#'   less than the number of observations n.
#'
#' @param filter.type Character string specifying the spectral filter type for
#'   response smoothing. The filter determines how high-frequency (rapid variation)
#'   components are attenuated relative to low-frequency (smooth) components.
#'   The smoothing parameter is selected automatically via generalized cross-validation.
#'   Available options:
#'   \describe{
#'     \item{\code{"heat_kernel"}}{Exponential decay exp(-η·λ). General-purpose
#'       filter with smooth, progressive attenuation. Most commonly used and
#'       recommended as the default choice. Equivalent to heat diffusion for
#'       time η.}
#'     \item{\code{"tikhonov"}}{Rational decay 1/(1+η·λ). Corresponds to
#'       Tikhonov regularization. Gentler than heat kernel, preserves more
#'       mid-frequency detail and linear trends. Good when moderate smoothing
#'       is sufficient.}
#'     \item{\code{"cubic_spline"}}{Spline smoothness 1/(1+η·λ²). Minimizes
#'       second derivatives of the fitted function. Produces very smooth,
#'       spline-like results. Excellent when the response should vary gradually
#'       and continuously.}
#'     \item{\code{"gaussian"}}{Super-exponential decay exp(-η·λ²). The most
#'       aggressive smoothing among exponential filters. Produces extremely
#'       smooth results with minimal oscillations. Use when maximum smoothness
#'       is desired and fine-scale features are noise.}
#'     \item{\code{"exponential"}}{Intermediate decay exp(-η·√λ). Less aggressive
#'       than heat kernel. Maintains more detail in mid-frequency range while
#'       still providing meaningful smoothing. Useful for preserving moderate-scale
#'       features while removing high-frequency noise.}
#'     \item{\code{"butterworth"}}{Smooth cutoff 1/(1+(λ/η)⁴). Fourth-order
#'       rational filter providing clear frequency separation with reduced
#'       ringing artifacts compared to ideal low-pass filters. Good balance
#'       between sharp cutoff and smooth transition.}
#'   }
#'   Default: \code{"heat_kernel"}
#'
#' @param t.diffusion Diffusion time parameter (non-negative). Controls how much
#'   vertex density evolves in each iteration. If 0 (default), automatically
#'   selected as \eqn{0.5/\lambda_2} where \eqn{\lambda_2} is the spectral gap.
#'   Larger values produce more aggressive density updates; smaller values are
#'   more conservative. Typical range is \eqn{0.1/\lambda_2} to \eqn{1.0/\lambda_2}.
#'   Set to 0 for automatic selection (recommended).
#'
#' @param density.uniform.pull Restoring force strength (non-negative).
#'   Controls how strongly vertex densities are pulled back toward uniform
#'   distribution during diffusion, preventing excessive concentration. If 0
#'   (default), automatically selected as 0.1/t.diffusion to maintain a 10:1
#'   ratio of diffusion to restoration. Larger values prevent concentration but
#'   may over-smooth; smaller values allow more geometric structure but risk
#'   collapse. Set to 0 for automatic selection (recommended). Former beta.damping.
#'
#' @param response.penalty.exp Response penalty exponent that controls the rate
#'     of edge weight reduction as response variation increases, via the penalty
#'     function \eqn{\Gamma(\Delta) = (1 + \Delta^2/\sigma^2)^{-\gamma}}. Larger
#'     values create steeper penalties (sharper adaptation to response changes);
#'     smaller values create gentler penalties. For smooth responses use
#'     0.5-0.75; for discontinuous responses use 1.5-2.0. When set to 0,
#'     disables response-coherence modulation entirely (geometry remains fixed
#'     based on feature space only). When negative (e.g., -1), triggers
#'     automatic selection via first-iteration GCV evaluation over a default
#'     grid {0.05, 0.1, 0.2, ..., 0.9, 1.0, 1.2, ..., 2.0}. Automatic selection
#'     efficiently determines the optimal value by evaluating GCV after a single
#'     iteration for each candidate, then uses the selected value for the full
#'     iterative refinement. Recommended when the response smoothness
#'     characteristics are unknown or when exploring new datasets. When
#'     positive, must be less than 2.0 and Default 0.3.
#'
#' @param use.counting.measure Logical scalar. If TRUE (default), uses uniform
#'   vertex weights (counting measure). If FALSE, uses distance-based weights
#'   inversely proportional to local k-NN density: \eqn{w(x) = (\epsilon + d_k(x))^{-\alpha}}.
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
#' @param density.alpha Exponent alpha in \eqn{[1,2]} for distance-based
#'                     reference measure formula \eqn{\mu(x) = (\epsilon + d_k(x))^{-\alpha}}.
#'                     Larger values (near 2) create stronger density-adaptive
#'                     weighting, concentrating measure on densely sampled
#'                     regions. Smaller values (near 1) produce more uniform
#'                     weighting. Only used when use_counting_measure = FALSE.
#'                     Default: 1.5.
#'
#' @param density.epsilon Regularization parameter epsilon > 0 in reference
#'                       measure formula \eqn{\mu(x) = (\epsilon + d_k(x))^{-\alpha}}.
#'                       Prevents numerical issues when nearest neighbor
#'                       distances are very small. Should be small relative to
#'                       typical k-NN distances but large enough for stability.
#'                       Only used when use_counting_measure = FALSE.
#'                       Default: 1e-10.
#'
#' @param max.ratio.threshold Numeric in the range 0 to 0.2.
#'     Geometric pruning removes an edge (i,j) when there exists an
#'     alternative path between i and j whose path/edge length ratio
#'     minus 1.0 is less than this threshold. This is a deviation
#'     threshold in the range 0 to 0.2. Internally we compare the
#'     path-to-edge ratio R to \eqn{1 + \delta}. Default 0.1.
#'
#' @param threshold.percentile Numeric in the range 0 to 1. Only edges with
#'     length above this percentile are considered for geometric pruning.
#'     Default 0.5.
#'
#' @param epsilon.y Numeric scalar, positive. Relative convergence threshold
#'   for response. Iteration stops when the relative change in fitted values
#'   \eqn{||\hat{y}^{(\ell)} - \hat{y}^{(\ell-1)}|| / ||\hat{y}^{(\ell-1)}||}
#'   falls below this value. Default 1e-4. Smaller values require more iterations
#'   but give more precise convergence.
#'
#' @param epsilon.rho Numeric scalar, positive. Relative convergence threshold
#'   for densities. Iteration stops when
#'   \eqn{\max_p ||\rho_p^{(\ell)} - \rho_p^{(\ell-1)}|| / ||\rho_p^{(\ell-1)}||}
#'   falls below this value across all dimensions. Default 1e-4.
#'
#' @param test.stage Integer for internal testing. Set to -1 for normal operation
#'     (default). Values 0-10 stop execution at intermediate stages for debugging.
#'     Not intended for end users.
#'
#' @param verbose Logical; print progress and timing.
#'
#' @return An object of class \code{c("knn.riem.fit", "list")}, which is
#'   a list containing:
#'   \itemize{
#'     \item \code{fitted.values}: Numeric vector of fitted values (selected
#'           from iteration with minimum GCV)
#'     \item \code{residuals}: Numeric vector of residuals
#'     \item \code{optimal.iteration}: Integer indicating which iteration's
#'           fitted values were selected (1-indexed)
#'     \item \code{graph}: List with graph structure information
#'     \item \code{iteration}: List with convergence information
#'     \item \code{parameters}: List of input parameters
#'     \item \code{y}: Original response values
#'     \item \code{gcv}: List with GCV information for each iteration
#'     \item \code{density}: List with density history
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
#'   \item Builds mass matrices \eqn{M_0} (diagonal) and \eqn{M_1} (full, via
#'         Gram determinants)
#'   \item Assembles Hodge Laplacian \eqn{L_0 = B_1 M_1^{-1} B_1^T M_0}
#'   \item Performs initial spectral filtering: \eqn{\hat{y}^{(0)} = V F_\eta(\Lambda) V^T y}
#' }
#'
#' \emph{Phase II: Iterative Refinement}
#'
#' Repeats until convergence or max.iterations reached:
#' \enumerate{
#'   \item \strong{Density diffusion}: Evolve vertex density via damped heat
#'         equation \eqn{(I + t(L_0 + \beta I))\rho_0 = \rho_{0,prev} + t\beta u}
#'   \item \strong{Edge density update}: Recompute edge densities from evolved
#'         vertex masses via neighborhood intersections
#'   \item \strong{Response-coherence modulation}: Reduce edge masses where
#'         response varies rapidly: \eqn{\rho_1([i,j]) \leftarrow \rho_1([i,j]) \cdot (1 + \Delta_{ij}^2/\sigma^2)^{-\gamma}}
#'   \item \strong{Metric update}: Rebuild \eqn{M_0}, \eqn{M_1} from modified densities
#'   \item \strong{Laplacian reassembly}: Recompute \eqn{L_0} with updated metric
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
#'   \item \strong{k}: Primary tuning parameter. Start with \eqn{k = \lceil\log(n)\rceil}
#'         and adjust based on cross-validation. Larger k smooths more.
#'   \item \strong{response.penalty.exp}: Increase (up to 2.0) for sharper response
#'         boundaries; decrease (down to 0.5) for gentler transitions.
#'   \item \strong{filter.type}: "heat_kernel" (default) for general use;
#'         "tikhonov" for better linear trend preservation; "cubic_spline"
#'         for smoother results on chain-like graphs.
#' }
#'
#' \strong{Automatic Parameter Selection}
#'
#' Several key parameters are automatically selected if not specified:
#' \itemize{
#'   \item \strong{t.diffusion}: Set to \eqn{0.5/\lambda_2} where \eqn{\lambda_2}
#'         is the spectral gap of the initial Laplacian
#'   \item \strong{density.uniform.pull}: Set to \eqn{0.1/t} to maintain 10:1
#'         diffusion to damping ratio
#'   \item \strong{response.penalty.exp}: When set to negative value (e.g., -1),
#'         automatically selected via first-iteration GCV evaluation across a
#'         grid of candidate values. This efficiently determines the optimal
#'         response-coherence strength by evaluating GCV after one iteration
#'         for each candidate, avoiding the computational cost of full convergence
#'         for each value
#' }
#'
#' The first-iteration GCV approach for gamma selection provides 5-10x speedup
#' compared to full cross-validation while maintaining comparable selection quality.
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item k-NN search: \eqn{O(n \log n)} using kd-trees (ANN library)
#'   \item Edge construction: \eqn{O(n k^2)} for intersection computations
#'   \item Per iteration: \eqn{O(n k^2 + m n_{eig})} where \eqn{m \approx n k}
#'         is edge count
#'   \item Total: \eqn{O(n k^2 + iter \times m n_{eig})}, typically manageable
#'         for \eqn{n \le 10,000}
#' }
#'
#' @section Memory Usage:
#' \itemize{
#'   \item Feature matrix: \eqn{O(n d)}
#'   \item Complex storage: \eqn{O(n k)} for vertices and edges
#'   \item Mass matrices: \eqn{O(n + m^2)} but \eqn{M_1} is sparse with
#'         \eqn{O(m k^2)} nonzeros
#'   \item Eigenpairs: \eqn{O(n \times n_{eigenpairs})}
#'   \item Peak usage: dominated by Laplacian operations, typically < 1GB
#'         for n = 5000
#' }
#'
#' @note \strong{Graph Connectivity Requirement:}
#'   The k-nearest neighbor graph must be fully connected (form a single
#'   connected component). If the graph is disconnected, the function will
#'   terminate with an error. To explore connectivity structure before fitting:
#'   \itemize{
#'     \item Use \code{create.iknn.graphs()} to examine connectivity across
#'           multiple k values
#'     \item Use \code{create.single.iknn.graph()} to inspect the graph
#'           structure for a specific k
#'     \item Examine connected components in the returned graph objects
#'   }
#'
#'   Common solutions for disconnected graphs:
#'   \itemize{
#'     \item Increase k to add more edges between vertices
#'     \item Remove isolated outlier observations
#'     \item Verify that your distance metric is appropriate for the data
#'   }
#'
#' @references
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
#' fit <- fit.knn.riem.graph.regression(X, y, k = 20, response.penalty.exp = 1.5)
#' y.hat <- fitted(fit)
#'
#' # Example 3: Using sparse matrix input
#' library(Matrix)
#' X.sparse <- Matrix(X, sparse = TRUE)
#' fit <- fit.knn.riem.graph.regression(X.sparse, y, k = 15)
#'
#' # Example 4: Automatic response.penalty.exp selection
#' # When optimal gamma is unknown, use automatic selection
#' set.seed(456)
#' X <- matrix(rnorm(n * 2), n, 2)
#' # Mix of smooth and discontinuous regions
#' y <- ifelse(X[,1] > 0, sin(3*X[,2]), 2*X[,2]^2) + rnorm(n, sd = 0.1)
#'
#' # Let algorithm select optimal gamma
#' fit.auto <- fit.knn.riem.graph.regression(X, y, k = 20,
#'                                            response.penalty.exp = -1)
#' cat("Selected gamma:", fit.auto$parameters$response.penalty.exp, "\n")
#' cat("Optimal iteration:", fit.auto$optimal.iteration, "\n")
#' }
#'
#' @seealso
#' \code{\link{fitted}} for extracting fitted values,
#' \code{\link{residuals}} for extracting residuals,
#' \code{\link{summary}} for model summaries
#'
#' @export
fit.knn.riem.graph.regression <- function(
    X,
    y,
    k,
    pca.dim = 100,
    variance.explained = 0.99,
    max.iterations = 10,
    n.eigenpairs = 10,
    filter.type = "heat_kernel",
    t.diffusion = 0,
    density.uniform.pull = 0,
    response.penalty.exp = 0,
    use.counting.measure = TRUE,
    density.normalization = 0,
    density.alpha = 1.5,
    density.epsilon = 1e-10,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.5,
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    test.stage = -1,
    verbose = TRUE
) {
    ## ==================== Feature Matrix Validation ====================

    ## Check X type
    is.dense <- is.matrix(X)
    is.sparse <- inherits(X, "dgCMatrix")

    if (!is.dense && !is.sparse) {
        stop("X must be a numeric matrix or sparse matrix (dgCMatrix from Matrix package)")
    }

    ## For dense matrix, check numeric
    if (is.dense && !is.numeric(X)) {
        stop("X must be numeric")
    }

    ## Extract dimensions
    n <- nrow(X)
    d <- ncol(X)

    ## Check dimensions
    if (is.null(n) || is.null(d) || n < 1 || d < 1) {
        stop("X must have valid dimensions (positive rows and columns)")
    }

    if (n < 10) {
        stop(sprintf("X must have at least 10 observations (got n=%d)", n))
    }

    ## Check for NA/Inf in X
    if (is.dense) {
        if (anyNA(X)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X))) {
            stop("X cannot contain infinite values")
        }
        ## Ensure storage mode is double
        storage.mode(X) <- "double"
    } else {
        ## For sparse matrix, check x slot
        if (anyNA(X@x)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X@x))) {
            stop("X cannot contain infinite values")
        }
    }

    ## ==================== Response Vector Validation ====================

    ## Check y type and length
    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    if (length(y) != n) {
        stop(sprintf("Length of y (%d) must equal number of rows in X (%d)",
                     length(y), n))
    }

    ## Check for NA/Inf in y
    if (anyNA(y)) {
        stop("y cannot contain NA values")
    }

    if (any(!is.finite(y))) {
        stop("y cannot contain infinite values")
    }

    ## Ensure y is a vector (not matrix)
    y <- as.vector(y)

    ## ==================== Parameter k Validation ====================

    if (!is.numeric(k) || length(k) != 1) {
        stop("k must be a single numeric value")
    }

    ## Convert to integer
    k <- as.integer(k)

    ## Check range
    if (k < 2) {
        stop(sprintf("k must be at least 2 (got k=%d)", k))
    }

    if (k >= n) {
        stop(sprintf("k must be less than n (got k=%d, n=%d)", k, n))
    }

    ## Validate density.alpha
    if (!is.numeric(density.alpha) || length(density.alpha) != 1) {
        stop("density.alpha must be a single numeric value")
    }
    if (!is.finite(density.alpha) || density.alpha < 1.0 || density.alpha > 2.0) {
        stop("density.alpha must be finite and in [1, 2], got ", density.alpha)
    }

    ## Validate density.epsilon
    if (!is.numeric(density.epsilon) || length(density.epsilon) != 1) {
        stop("density.epsilon must be a single numeric value")
    }
    if (!is.finite(density.epsilon) || density.epsilon <= 0) {
        stop("density.epsilon must be finite and positive, got ", density.epsilon)
    }

    ## ==================== Logical Parameter Validation ====================

    if (!is.logical(use.counting.measure) || length(use.counting.measure) != 1) {
        stop("use.counting.measure must be TRUE or FALSE")
    }

    if (is.na(use.counting.measure)) {
        stop("use.counting.measure cannot be NA")
    }

    ## ==================== Numeric Parameter Validation ====================

    ## density.normalization
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

    ## t.diffusion
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

    ## density.uniform.pull
    if (!is.numeric(density.uniform.pull) || length(density.uniform.pull) != 1) {
        stop("density.uniform.pull must be a single numeric value")
    }

    if (is.na(density.uniform.pull) || !is.finite(density.uniform.pull)) {
        stop("density.uniform.pull cannot be NA or infinite")
    }

    if (density.uniform.pull < 0) {
        stop(sprintf("density.uniform.pull must be non-negative (got %.3f; use 0 for auto)",
                     density.uniform.pull))
    }

    ## response.penalty.exp
    if (!is.numeric(response.penalty.exp) || length(response.penalty.exp) != 1) {
        stop("response.penalty.exp must be a single numeric value")
    }

    if (is.na(response.penalty.exp) || !is.finite(response.penalty.exp)) {
        stop("response.penalty.exp cannot be NA or infinite")
    }

    if (response.penalty.exp > 2.0) {
        warning(sprintf("response.penalty.exp=%.3f must be less than 2.0.\n",
                        response.penalty.exp))
    }

    ## n.eigenpairs
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

    ## Warn if n.eigenpairs seems too small
    if (n.eigenpairs < 50 && n > 100) {
        warning(sprintf(paste0("n.eigenpairs=%d may be too small for n=%d observations.\n",
                              "Consider n.eigenpairs >= 100 for better frequency resolution."),
                       n.eigenpairs, n))
    }

    ## epsilon.y
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

    ## epsilon.rho
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

    ## max.iterations
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

    ## max.ratio.threshold
    if (!is.numeric(max.ratio.threshold) || length(max.ratio.threshold) != 1) {
        stop("max.ratio.threshold must be numeric.")
    }
        
    if (max.ratio.threshold < 0 || max.ratio.threshold >= 0.2) {
        stop("max.ratio.threshold must be in [0, 0.2).")
    }
        
    ## threshold.percentile
    if (!is.numeric(threshold.percentile) || length(threshold.percentile) != 1 ||
        threshold.percentile < 0 || threshold.percentile > 1) {
        stop("threshold.percentile must be in [0, 1].")
    }

    ## verbose
    if (!is.logical(verbose) || length(verbose) != 1)
        stop("verbose must be TRUE/FALSE.")

    ## pca.dim
    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim))
            stop("pca.dim must be a positive integer or NULL.")
    }

    ## variance.explained
    if (!is.null(variance.explained)) {
        if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
            variance.explained <= 0 || variance.explained > 1)
            stop("variance.explained must be in (0, 1], or NULL.")
    }

    ## PCA (optional)
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose) message("High-dimensional data detected. Performing PCA.")
        original_dim <- ncol(X)
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(
                X, variance.threshold = variance.explained, max.components = pca.dim
            )
            n_components <- pca_analysis$n.components
            if (verbose) {
                message(sprintf("Using %d PCs (explains %.2f%% variance)",
                                n_components, 100 * pca_analysis$variance.explained))
            }
            X <- pca.project(X, pca_analysis$pca.result, n_components)
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance
            )
        } else {
            if (verbose) message(sprintf("Projecting to first %d PCs", pca.dim))
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained
            )
        }
    }

        
    ## ==================== String Parameter Validation ====================

    ## Match filter.type
    filter.type <- match.arg(filter.type)

    ## ==================== Call C++ Function ====================

    fit <- .Call(
        S_fit_knn_riem_graph_regression,
        X,
        as.double(y),
        as.integer(k + 1L), # this is to account for the fact that ANN library is set up to return for query point as the first elements of the list of kNN's
        as.logical(use.counting.measure),
        as.double(density.normalization),
        as.double(t.diffusion),
        as.double(density.uniform.pull),
        as.double(response.penalty.exp),
        as.integer(n.eigenpairs),
        as.character(filter.type),
        as.double(epsilon.y),
        as.double(epsilon.rho),
        as.integer(max.iterations),
        as.double(max.ratio.threshold + 1.0),
        as.double(threshold.percentile),
        as.double(density.alpha),
        as.double(density.epsilon),
        as.integer(test.stage),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## ==================== Post-Processing ====================

    ## Add class
    class(fit) <- c("knn.riem.fit", "riem.dcx")

    ## Store call for reproducibility
    attr(fit, "call") <- match.call()

    ## Store key parameters for summary()
    attr(fit, "n") <- n
    attr(fit, "d") <- d
    attr(fit, "k") <- k
    if (!is.null(pca_info)) attr(fit, "pca") <- pca_info

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
summary.knn.riem.fit <- function(object, ...) {
    # Extract key information
    resid <- object$residuals
    n <- length(object$y)

    # Compute fit quality metrics
    rss <- sum(resid^2)
    tss <- sum((object$y - mean(object$y))^2)
    r_squared <- 1 - rss / tss
    adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - object$parameters$k - 1)
    rmse <- sqrt(mean(resid^2))
    mae <- mean(abs(resid))

    # Get optimal iteration info
    optimal_iter <- object$optimal.iteration
    total_iters <- object$iteration$n.iterations

    # Get GCV information
    gcv_optimal <- object$gcv$gcv.optimal[optimal_iter]
    gcv_initial <- object$gcv$gcv.optimal[1]
    gcv_improvement <- (gcv_initial - gcv_optimal) / gcv_initial * 100

    # Check if parameters were auto-selected
    params <- object$parameters
    auto_selected <- list(
        t.diffusion = (params$t.diffusion > 0),  # Was auto-selected if positive
        density.uniform.pull = (params$density.uniform.pull > 0),
        response.penalty.exp = attr(object, "gamma.auto.selected")  # Need to store this
    )

    # Graph structure info
    n_vertices <- object$graph$n.vertices
    n_edges <- object$graph$n.edges
    avg_degree <- 2 * n_edges / n_vertices

    # Convergence diagnostics
    final_response_change <- if (length(object$iteration$response.changes) > 0) {
        tail(object$iteration$response.changes, 1)
    } else {
        NA
    }

    structure(
        list(
            call = attr(object, "call"),

            # Data dimensions
            n = n,
            d = attr(object, "d"),

            # Model parameters
            k = params$k,
            response.penalty.exp = params$response.penalty.exp,
            filter.type = params$filter.type,
            n.eigenpairs = params$n.eigenpairs,

            # Auto-selection flags
            auto.selected = auto_selected,

            # Graph structure
            graph = list(
                n.vertices = n_vertices,
                n.edges = n_edges,
                avg.degree = avg_degree,
                density = n_edges / choose(n_vertices, 2)
            ),

            # Convergence
            converged = object$iteration$converged,
            n.iterations = total_iters,
            optimal.iteration = optimal_iter,
            final.response.change = final_response_change,

            # Fit quality
            residuals = resid,
            sigma = sd(resid),
            rmse = rmse,
            mae = mae,
            r.squared = r_squared,
            adj.r.squared = adj_r_squared,

            # GCV diagnostics
            gcv.optimal = gcv_optimal,
            gcv.initial = gcv_initial,
            gcv.improvement.pct = gcv_improvement,

            # Response characteristics
            y.range = range(object$y),
            y.mean = mean(object$y),
            y.sd = sd(object$y)
        ),
        class = "summary.knn.riem.fit"
    )
}

#' @export
print.summary.knn.riem.fit <- function(x, digits = 4, ...) {
    cat("\n")
    cat("========================================================\n")
    cat("k-NN Riemannian Graph Regression Summary\n")
    cat("========================================================\n\n")

    # Call
    cat("Call:\n")
    print(x$call)
    cat("\n")

    # Data and Model
    cat("Data and Model:\n")
    cat(sprintf("  Observations:        n = %d\n", x$n))
    cat(sprintf("  Features:            d = %d\n", x$d))
    cat(sprintf("  Neighbors:           k = %d\n", x$k))
    cat(sprintf("  Filter type:         %s\n", x$filter.type))
    cat(sprintf("  Eigenpairs:          %d\n", x$n.eigenpairs))
    cat("\n")

    # Parameter selection
    cat("Parameter Selection:\n")
    cat(sprintf("  Response penalty:    gamma = %.3f", x$response.penalty.exp))
    if (!is.null(x$auto.selected$response.penalty.exp) && x$auto.selected$response.penalty.exp) {
        cat(" [auto-selected]")
    }
    cat("\n")
    cat(sprintf("  Diffusion time:      t = %.4f", x$k))  # Need actual value
    if (!is.null(x$auto.selected$t.diffusion) && x$auto.selected$t.diffusion) {
        cat(" [auto-selected]")
    }
    cat("\n\n")

    # Graph structure
    cat("Graph Structure:\n")
    cat(sprintf("  Vertices:            %d\n", x$graph$n.vertices))
    cat(sprintf("  Edges:               %d\n", x$graph$n.edges))
    cat(sprintf("  Avg degree:          %.1f\n", x$graph$avg.degree))
    cat(sprintf("  Graph density:       %.4f\n", x$graph$density))
    cat("\n")

    # Convergence
    cat("Convergence:\n")
    cat(sprintf("  Status:              %s\n",
                if(x$converged) "CONVERGED" else "Max iterations reached"))
    cat(sprintf("  Total iterations:    %d\n", x$n.iterations))
    cat(sprintf("  Optimal iteration:   %d (min GCV)\n", x$optimal.iteration))
    if (!is.na(x$final.response.change)) {
        cat(sprintf("  Final change:        %.2e\n", x$final.response.change))
    }
    cat("\n")

    # Fit quality
    cat("Fit Quality:\n")
    cat(sprintf("  R-squared:           %.4f\n", x$r.squared))
    cat(sprintf("  Adj R-squared:       %.4f\n", x$adj.r.squared))
    cat(sprintf("  RMSE:                %.4f\n", x$rmse))
    cat(sprintf("  MAE:                 %.4f\n", x$mae))
    cat(sprintf("  Residual SD:         %.4f\n", x$sigma))
    cat("\n")

    # GCV diagnostics
    cat("GCV Diagnostics:\n")
    cat(sprintf("  Optimal GCV:         %.4e\n", x$gcv.optimal))
    cat(sprintf("  Initial GCV:         %.4e\n", x$gcv.initial))
    cat(sprintf("  Improvement:         %.1f%%\n", x$gcv.improvement.pct))
    cat("\n")

    # Residual summary
    cat("Residuals:\n")
    print(summary(x$residuals), digits = digits)
    cat("\n")

    # Response summary
    cat("Response:\n")
    cat(sprintf("  Range:               [%.3f, %.3f]\n",
                x$y.range[1], x$y.range[2]))
    cat(sprintf("  Mean:                %.3f\n", x$y.mean))
    cat(sprintf("  SD:                  %.3f\n", x$y.sd))

    cat("\n")
    cat("========================================================\n")

    invisible(x)
}

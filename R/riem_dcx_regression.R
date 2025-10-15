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
#'     vertex density evolves in each iteration. If 0 (default), automatically
#'     selected as \eqn{0.5/\lambda_2} where \eqn{\lambda_2} is the spectral
#'     gap. Larger values produce more aggressive density updates; smaller
#'     values are more conservative. Typical range is \eqn{0.1/\lambda_2} to
#'     \eqn{1.0/\lambda_2}. Set to 0 for automatic selection (recommended).
#'
#' @param density.uniform.pull Restoring force strength (non-negative). Controls
#'     how strongly vertex densities are pulled back toward uniform distribution
#'     during diffusion, preventing excessive concentration. If 0 (default),
#'     automatically selected as 0.1/t.diffusion to maintain a 10:1 ratio of
#'     diffusion to restoration. Larger values prevent concentration but may
#'     over-smooth; smaller values allow more geometric structure but risk
#'     collapse. Set to 0 for automatic selection (recommended). Former
#'     beta.damping.
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
#' @param t.scale.factor Diffusion scale factor (positive numeric). Controls the
#'     dimensionless product \eqn{t \cdot \lambda_2} when \code{t.diffusion = 0}
#'     (automatic selection), where \eqn{\lambda_2} is the spectral gap of the
#'     graph Laplacian. This parameter determines how aggressively density
#'     diffuses per iteration relative to the graph's intrinsic geometric scale.
#'     Conservative values \eqn{(0.1-0.2)} produce slow but stable density
#'     evolution, requiring more iterations but reducing risk of oscillation.
#'     Moderate values \eqn{(0.4-0.6)} balance speed and stability for typical
#'     applications. Aggressive values \eqn{(0.8-1.0)} accelerate convergence but may
#'     cause oscillatory behavior when geometric-response feedback is strong.
#'     The scale factor is ignored when \code{t.diffusion > 0} (manual
#'     specification). Typical range: \eqn{[0.1, 1.0]}. Default: 0.5
#'
#' @param beta.coef.factor Damping coefficient factor (positive numeric).
#'     Controls the product \eqn{\beta \cdot t} when \code{density.uniform.pull
#'     = 0} (automatic selection). This parameter determines the strength of the
#'     restoring force toward uniform distribution during each iteration step.
#'     The damping prevents excessive concentration of density in small regions
#'     while allowing genuine accumulation in well-connected areas. Light
#'     damping (0.05) permits strong concentration and is appropriate when the
#'     data distribution is expected to be highly non-uniform. Moderate damping
#'     (0.1) provides balanced geometric adaptation for most applications.
#'     Strong damping (0.2-0.3) prevents concentration aggressively and is
#'     appropriate when maintaining dispersed density is important or when the
#'     response varies smoothly without sharp boundaries. The coefficient factor
#'     is ignored when \code{density.uniform.pull > 0} (manual specification). A
#'     coefficient of 0.1 means the restoring force contributes approximately
#'     10\% as much as the identity term in the system matrix. Typical range:
#'     \eqn{[0.05, 0.3]}. Default: 0.1
#'
#' @param use.counting.measure Logical scalar. If TRUE, uses uniform vertex
#'     weights (counting measure). If FALSE, uses distance-based weights
#'     inversely proportional to local k-NN density: \eqn{w(x) = (\epsilon +
#'     d_k(x))^{-\alpha}}. Distance-based weights are useful when sampling
#'     density varies across the feature space. Default: FALSE.
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
    t.scale.factor = 0.5,
    beta.coef.factor = 0.1,
    use.counting.measure = FALSE,
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

    ## ==================== Parameter use.counting.measure Validation ====================

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

    ## t.scale.factor
    if (!is.numeric(t.scale.factor) || length(t.scale.factor) != 1) {
        stop("t.scale.factor must be a single numeric value")
    }

    if (t.scale.factor <= 0) {
        stop("t.scale.factor must be positive (got ", t.scale.factor,
             "). Typical range: [0.1, 1.0]")
    }

    if (t.scale.factor < 0.05) {
        warning("Very small t.scale.factor (", t.scale.factor,
                "): convergence will be extremely slow. Consider increasing to at least 0.1.")
    }

    if (t.scale.factor > 2.0) {
        warning("Very large t.scale.factor (", t.scale.factor,
                "): may cause oscillations. Consider reducing to below 1.0.")
    }

    ## beta.coef.factor
    if (!is.numeric(beta.coef.factor) || length(beta.coef.factor) != 1) {
        stop("beta.coef.factor must be a single numeric value")
    }

    if (beta.coef.factor <= 0) {
        stop("beta.coef.factor must be positive (got ", beta.coef.factor,
             "). Typical range: [0.05, 0.3]")
    }

    if (beta.coef.factor < 0.01) {
        warning("Very small beta.coef.factor (", beta.coef.factor,
                "): density may collapse. Consider increasing to at least 0.05.")
    }

    if (beta.coef.factor > 0.5) {
        warning("Very large beta.coef.factor (", beta.coef.factor,
                "): may over-suppress structure. Consider reducing to below 0.3.")
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
        as.double(t.scale.factor),
        as.double(beta.coef.factor),
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

#' Diagnostic Plots for k-NN Riemannian Graph Regression
#'
#' @param results.list List of fitted model objects from fit.knn.riem.graph.regression
#' @param k.values Vector of k values corresponding to results.list
#' @param response Vector of observed response values (for computing relative ranges)
#' @param main.title Optional main title for the entire plot grid
#' @param mark.optimal Logical; if TRUE, mark the k selected by GCV in red
#'
#' @return Invisibly returns a data frame with all diagnostic metrics
#' @export
k.diagnostics.plots <- function(results.list,
                                 k.values,
                                 response,
                                 main.title = "k-NN Bandwidth Diagnostics",
                                 mark.optimal = TRUE) {

  n.models <- length(results.list)
  if (length(k.values) != n.models) {
    stop("Length of k.values must match length of results.list")
  }

    ## Compute all diagnostics
    gcv.scores <- numeric(n.models)
    range.lower <- numeric(n.models)
    range.upper <- numeric(n.models)
    stability.corr <- rep(NA_real_, n.models)
    n.extrema <- numeric(n.models)

    y.mean <- mean(response)

    for (j in seq_len(n.models)) {
        fit <- results.list[[j]]

        ## GCV score
        gcv.scores[j] <- fit$gcv$gcv.optimal[fit$optimal.iteration]

        ## Fitted values range (winsorized)
        fitted.wins <- winsorize(fit$fitted.values, p = 0.02)
        range.lower[j] <- min(fitted.wins) / y.mean
        range.upper[j] <- max(fitted.wins) / y.mean

        ## Stability correlation
        if (j > 1) {
            fit.prev <- results.list[[j - 1]]
            stability.corr[j] <- cor(fit.prev$fitted.values, fit$fitted.values)
        }

        ## Number of local extrema
        extrema.res <- compute.extrema.hop.nbhds(
            fit$graph$adj.list,
            fit$graph$edge.length.list,
            fit$fitted.values
        )
        n.extrema[j] <- nrow(extrema.res$extrema_df)
    }

    ## Find optimal k by GCV
    optimal.idx <- which.min(gcv.scores)
    optimal.k <- k.values[optimal.idx]

    ## Set up 2x2 plotting grid
    old.par <- par(mfrow = c(2, 2),
                   mar = c(4, 4, 2.5, 1),
                   oma = c(0, 0, 2, 0))
    on.exit(par(old.par))

    ## Plot 1: GCV scores (expressed as deviation from minimum)
    gcv.min <- min(gcv.scores)
    gcv.relative <- (gcv.scores - gcv.min) / gcv.min * 100  # percent increase from minimum

    plot(k.values, gcv.relative,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "GCV increase from minimum (%)",
         main = "Cross-Validation Criterion",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
        abline(h = 0, col = "gray40", lty = 1)
        text(optimal.k, max(gcv.relative),
             labels = paste0("k = ", optimal.k),
             pos = 4, col = "firebrick", cex = 0.9)
    }

    ## Plot 2: Fitted values range
    y.range <- range(c(range.lower, range.upper))
    plot(k.values, range.lower,
         type = "n",
         ylim = y.range,
         xlab = "Neighborhood size k",
         ylab = "Relative range (normalized by mean response)",
         main = "Prediction Range Stability",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    ## Add shaded region
    polygon(c(k.values, rev(k.values)),
            c(range.lower, rev(range.upper)),
            col = adjustcolor("steelblue", alpha.f = 0.3),
            border = NA)

    lines(k.values, range.lower, type = "b", pch = 19, col = "steelblue")
    lines(k.values, range.upper, type = "b", pch = 19, col = "steelblue")
    abline(h = 1, col = "gray40", lwd = 1.5, lty = 1)

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Add legend
    legend("topright",
           legend = c("Lower bound", "Upper bound", "Mean response"),
           col = c("steelblue", "steelblue", "gray40"),
           lty = c(1, 1, 1),
           lwd = c(1, 1, 1.5),
           pch = c(19, 19, NA),
           bty = "n",
           cex = 0.8)

    ## Plot 3: Stability correlation
    plot(k.values, stability.corr,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "Correlation with k-1 fit",
         main = "Fit Stability Across k",
         las = 1,
         ylim = c(0, 1),
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Plot 4: Number of extrema
    plot(k.values, n.extrema,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "Count",
         main = "Geometric Complexity (Local Extrema)",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Add overall title
    mtext(main.title, outer = TRUE, cex = 1.2, font = 2)

    ## Return diagnostic data frame
    diagnostics.df <- data.frame(
        k = k.values,
        gcv = gcv.scores,
        range.lower = range.lower,
        range.upper = range.upper,
        stability.corr = stability.corr,
        n.extrema = n.extrema,
        is.optimal = k.values == optimal.k
    )

    invisible(diagnostics.df)
}

#' Extract Optimal k from Diagnostics Data Frame
#'
#' @param diag.df Data frame returned by k.diagnostics.plots
#'
#' @return A list with optimal k value and its index
#' @export
get.rcx.optimal.k <- function(diag.df) {
  if (!any(diag.df$is.optimal)) {
    stop("No optimal k found in diagnostics data frame")
  }

  optimal.idx <- which(diag.df$is.optimal)
  optimal.k <- diag.df$k[optimal.idx]

  list(
    k = optimal.k,
    idx = optimal.idx,
    gcv = diag.df$gcv[optimal.idx]
  )
}

#' Refit Riemannian Graph Regression Model with New Response(s)
#'
#' @description
#' Apply the learned spectral filter from a fitted model to new response
#' variable(s) using the same graph geometry. This is computationally
#' efficient as it reuses the cached eigendecomposition without rebuilding
#' the graph or recomputing the Laplacian.
#'
#' @param fitted.model A fitted model object from \code{fit.knn.riem.graph.regression}
#' @param y.new New response data. Can be:
#'   \itemize{
#'     \item A numeric vector of length n (single response)
#'     \item A numeric matrix of dimension n × p (p response variables)
#'   }
#'   where n must match the number of observations in the original fit.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{fitted.values}: Fitted values (vector if y.new is vector,
#'           matrix if y.new is matrix)
#'     \item \code{residuals}: Residuals (same shape as y.new)
#'     \item \code{n.responses}: Number of response variables (1 for vector,
#'           p for matrix)
#'   }
#'
#' @details
#' The function applies the learned spectral operator
#' \eqn{\hat{y} = V \text{diag}(F_\eta(\Lambda)) V^T y} where V and
#' \eqn{F_\eta(\Lambda)} are cached from the optimal iteration of the
#' original fit. For matrix input, the operation is applied column-wise,
#' fitting all response variables simultaneously using the same geometry.
#'
#' This is particularly useful for:
#' \itemize{
#'   \item Bootstrap resampling (fit multiple bootstrap responses)
#'   \item Permutation testing (fit permuted responses)
#'   \item Multivariate regression (multiple outcomes sharing feature space)
#'   \item Cross-validation (refit on different subsets)
#' }
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' fit <- fit.knn.riem.graph.regression(X, y, k = 15)
#'
#' # Single new response
#' y.new <- y + rnorm(length(y), sd = 0.1)
#' refit.single <- refit.knn.riem.graph.regression(fit, y.new)
#'
#' # Multiple responses (e.g., bootstrap)
#' Y.boot <- replicate(100, sample(y, replace = TRUE))
#' refit.multi <- refit.knn.riem.graph.regression(fit, Y.boot)
#' dim(refit.multi$fitted.values)  # n × 100
#' }
#'
#' @export
refit.knn.riem.graph.regression <- function(fitted.model, y.new) {

    ## Validate input
    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be a 'knn.riem.fit' object from fit.knn.riem.graph.regression()")
    }

    if (is.null(fitted.model$spectral)) {
        stop("Fitted model does not contain spectral component. ",
             "Cannot refit without cached eigendecomposition.")
    }

    ## Extract spectral components
    V <- fitted.model$spectral$eigenvectors
    f_lambda <- fitted.model$spectral$filtered.eigenvalues
    n <- nrow(V)

    ## Validate and prepare y.new
    is_matrix <- is.matrix(y.new)

    if (is_matrix) {
        if (nrow(y.new) != n) {
            stop(sprintf("Number of rows in y.new (%d) must match fitted model (%d)",
                         nrow(y.new), n))
        }
        n_responses <- ncol(y.new)
    } else {
        ## Convert vector to column matrix for uniform handling
        if (length(y.new) != n) {
            stop(sprintf("Length of y.new (%d) must match fitted model (%d)",
                         length(y.new), n))
        }
        y.new <- as.matrix(y.new, ncol = 1)
        n_responses <- 1
    }

    ## Apply spectral filter: Ŷ = V diag(F_η(λ)) V^T Y
    ## This formula works for both single and multiple responses
    Vt_Y <- crossprod(V, y.new)           # V^T Y: (m × n) × (n × p) = m × p
    filtered <- f_lambda * Vt_Y           # Element-wise: m × p (broadcasts f_lambda)
    Y.hat <- V %*% filtered               # V × filtered: (n × m) × (m × p) = n × p

    ## Compute residuals
    residuals <- y.new - Y.hat

    ## Return in appropriate format
    result <- list(
        fitted.values = if (n_responses == 1) as.vector(Y.hat) else Y.hat,
        residuals = if (n_responses == 1) as.vector(residuals) else residuals,
        n.responses = n_responses
    )

    class(result) <- c("knn.riem.refit", "list")
    return(result)
}

#' @export
print.knn.riem.refit <- function(x, ...) {
  cat("\nRefitted Riemannian Graph Regression\n")
  cat("====================================\n\n")

  if (x$n.responses == 1) {
    cat(sprintf("Single response variable (n = %d)\n", length(x$fitted.values)))
    cat(sprintf("RMSE: %.4f\n", sqrt(mean(x$residuals^2))))
    cat(sprintf("MAE:  %.4f\n", mean(abs(x$residuals))))
  } else {
    cat(sprintf("Multiple responses (n = %d, p = %d)\n",
                nrow(x$fitted.values), x$n.responses))
    rmse_by_col <- apply(x$residuals, 2, function(r) sqrt(mean(r^2)))
    cat(sprintf("Mean RMSE across responses: %.4f\n", mean(rmse_by_col)))
    cat(sprintf("RMSE range: [%.4f, %.4f]\n", min(rmse_by_col), max(rmse_by_col)))
  }

  invisible(x)
}

#' @export
summary.knn.riem.refit <- function(object, ...) {
  n_resp <- object$n.responses
  n_obs <- if (n_resp == 1) length(object$fitted.values) else nrow(object$fitted.values)

  if (n_resp == 1) {
    # Single response summary
    y_hat <- object$fitted.values
    resid <- object$residuals
    y_obs <- y_hat + resid  # Reconstruct observed values

    # Basic fit quality metrics
    rss <- sum(resid^2)
    tss <- sum((y_obs - mean(y_obs))^2)
    r_squared <- 1 - rss / tss
    adj_r_squared <- 1 - (1 - r_squared) * (n_obs - 1) / (n_obs - 2)
    rmse <- sqrt(mean(resid^2))
    mae <- mean(abs(resid))

    # Normalization factors
    y_range <- diff(range(y_obs))
    y_sd <- sd(y_obs)
    y_mean <- mean(y_obs)

    # Relative metrics
    nrmse_range <- if (y_range > 0) rmse / y_range else NA
    nrmse_sd <- if (y_sd > 0) rmse / y_sd else NA
    nmae_range <- if (y_range > 0) mae / y_range else NA

    # CV-RMSE (only if y is positive or mean is not near zero)
    cv_rmse <- if (abs(y_mean) > 1e-10) rmse / abs(y_mean) else NA

    result <- list(
      n = n_obs,
      n.responses = 1,

      # Absolute metrics
      rmse = rmse,
      mae = mae,
      r.squared = r_squared,
      adj.r.squared = adj_r_squared,

      # Relative metrics
      nrmse.range = nrmse_range,
      nrmse.sd = nrmse_sd,
      nmae.range = nmae_range,
      cv.rmse = cv_rmse,

      # Response characteristics
      residual.range = range(resid),
      fitted.range = range(y_hat),
      y.range = range(y_obs),
      y.mean = y_mean,
      y.sd = y_sd
    )

  } else {
    # Multiple responses summary
    Y_hat <- object$fitted.values
    Resid <- object$residuals
    Y_obs <- Y_hat + Resid

    # Compute metrics for each response
    rmse_vec <- apply(Resid, 2, function(r) sqrt(mean(r^2)))
    mae_vec <- apply(Resid, 2, function(r) mean(abs(r)))

    r_squared_vec <- sapply(1:n_resp, function(j) {
      rss <- sum(Resid[, j]^2)
      tss <- sum((Y_obs[, j] - mean(Y_obs[, j]))^2)
      1 - rss / tss
    })

    # Compute normalized metrics for each response
    nrmse_range_vec <- sapply(1:n_resp, function(j) {
      y_range <- diff(range(Y_obs[, j]))
      if (y_range > 0) rmse_vec[j] / y_range else NA
    })

    nrmse_sd_vec <- sapply(1:n_resp, function(j) {
      y_sd <- sd(Y_obs[, j])
      if (y_sd > 0) rmse_vec[j] / y_sd else NA
    })

    result <- list(
      n = n_obs,
      n.responses = n_resp,

      # Absolute metrics
      rmse.mean = mean(rmse_vec),
      rmse.sd = sd(rmse_vec),
      rmse.range = range(rmse_vec),

      mae.mean = mean(mae_vec),
      mae.sd = sd(mae_vec),
      mae.range = range(mae_vec),

      r.squared.mean = mean(r_squared_vec),
      r.squared.sd = sd(r_squared_vec),
      r.squared.range = range(r_squared_vec),

      # Normalized metrics (averaged across responses)
      nrmse.range.mean = mean(nrmse_range_vec, na.rm = TRUE),
      nrmse.sd.mean = mean(nrmse_sd_vec, na.rm = TRUE),

      # Per-response details
      per.response = data.frame(
        response = 1:n_resp,
        rmse = rmse_vec,
        mae = mae_vec,
        r.squared = r_squared_vec,
        nrmse.range = nrmse_range_vec,
        nrmse.sd = nrmse_sd_vec
      )
    )
  }

  class(result) <- "summary.knn.riem.refit"
  return(result)
}

#' @export
print.summary.knn.riem.refit <- function(x, digits = 4, ...) {
  cat("\n")
  cat("========================================================\n")
  cat("Refitted Riemannian Graph Regression Summary\n")
  cat("========================================================\n\n")

  if (x$n.responses == 1) {
    # Single response output
    cat("Fit Quality (Absolute):\n")
    cat(sprintf("  Observations:    n = %d\n", x$n))
    cat(sprintf("  RMSE:            %.4f\n", x$rmse))
    cat(sprintf("  MAE:             %.4f\n", x$mae))
    cat(sprintf("  R-squared:       %.4f\n", x$r.squared))
    cat(sprintf("  Adj. R-squared:  %.4f\n\n", x$adj.r.squared))

    cat("Fit Quality (Relative):\n")
    if (!is.na(x$nrmse.range)) {
      cat(sprintf("  NRMSE (range):   %.2f%%\n", x$nrmse.range * 100))
    }
    ## if (!is.na(x$nrmse.sd)) {
    ##   cat(sprintf("  NRMSE (std dev): %.2f%%\n", x$nrmse.sd * 100))
    ## }
    if (!is.na(x$nmae.range)) {
      cat(sprintf("  NMAE (range):    %.2f%%\n", x$nmae.range * 100))
    }
    ## if (!is.na(x$cv.rmse)) {
    ##   cat(sprintf("  CV-RMSE:         %.2f%%\n", x$cv.rmse * 100))
    ## }
    cat("\n")

    cat("Response Characteristics:\n")
    cat(sprintf("  Mean:            %.4f\n", x$y.mean))
    cat(sprintf("  Std Dev:         %.4f\n", x$y.sd))
    cat(sprintf("  Observed range:  [%.4f, %.4f]\n", x$y.range[1], x$y.range[2]))
    cat(sprintf("  Fitted range:    [%.4f, %.4f]\n", x$fitted.range[1], x$fitted.range[2]))
    cat(sprintf("  Residual range:  [%.4f, %.4f]\n", x$residual.range[1], x$residual.range[2]))

  } else {
    # Multiple responses output
    cat(sprintf("Number of responses:       p = %d\n", x$n.responses))
    cat(sprintf("Observations per response: n = %d\n\n", x$n))

    cat("Fit Quality - Absolute (across all responses):\n")
    cat(sprintf("  RMSE:       %.4f ± %.4f  [%.4f, %.4f]\n",
                x$rmse.mean, x$rmse.sd, x$rmse.range[1], x$rmse.range[2]))
    cat(sprintf("  MAE:        %.4f ± %.4f  [%.4f, %.4f]\n",
                x$mae.mean, x$mae.sd, x$mae.range[1], x$mae.range[2]))
    cat(sprintf("  R-squared:  %.4f ± %.4f  [%.4f, %.4f]\n\n",
                x$r.squared.mean, x$r.squared.sd,
                x$r.squared.range[1], x$r.squared.range[2]))

    cat("Fit Quality - Relative (averaged across responses):\n")
    if (!is.na(x$nrmse.range.mean)) {
      cat(sprintf("  NRMSE (range):   %.2f%%\n", x$nrmse.range.mean * 100))
    }
    if (!is.na(x$nrmse.sd.mean)) {
      cat(sprintf("  NRMSE (std dev): %.2f%%\n\n", x$nrmse.sd.mean * 100))
    }

    # Show per-response details if not too many
    if (x$n.responses <= 10) {
      cat("Per-Response Details:\n")
      print(x$per.response, digits = digits, row.names = FALSE)
    } else {
      cat(sprintf("Per-response details available in $per.response (use summary(fit)$per.response)\n"))
      cat("Showing first 5 and last 5 responses:\n\n")
      to_show <- rbind(
        head(x$per.response, 5),
        tail(x$per.response, 5)
      )
      print(to_show, digits = digits, row.names = FALSE)
    }
  }

  cat("\n")
  cat("========================================================\n")
  cat("Normalized Error Metrics:\n")
  cat("  NRMSE (range)   = RMSE / range(y)  [scale-free error]\n")
  ##cat("  NRMSE (std dev) = RMSE / sd(y)     [error relative to variability]\n")
  cat("  NMAE (range)    = MAE / range(y)   [scale-free absolute error]\n")
  ##cat("  CV-RMSE         = RMSE / |mean(y)| [coefficient of variation]\n")
  cat("\n")
  cat("Interpretation: Normalized metrics express error as a percentage\n")
  cat("of response scale. Values < 10% typically indicate excellent fit,\n")
  cat("10-20% good fit, > 30% may warrant investigation.\n")
  cat("\n")
  cat("Note: This summary reflects fit quality using a pre-learned\n")
  cat("      spectral operator from the original model. No iteration\n")
  cat("      or parameter optimization was performed.\n")

  invisible(x)
}

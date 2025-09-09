#' Spectral Filtering for Graph Signals
#'
#' @description
#' Implements a comprehensive framework for spectral filtering of signals defined on graph vertices.
#' Supports multiple types of graph Laplacians, spectral filters, and parameter configurations
#' to achieve different smoothing characteristics.
#'
#' @details
#' The general spectral filtering process follows these steps:
#' \enumerate{
#'   \item Constructs a specific graph Laplacian operator based on the selected Laplacian type
#'   \item Computes the eigendecomposition of this Laplacian (or its transformation)
#'   \item Projects the signal onto the eigenbasis (Graph Fourier Transform)
#'   \item Applies filter weights based on the eigenvalues and selected filter type
#'   \item Reconstructs smoothed signals for a grid of filter parameter values
#'   \item Selects the optimal parameter using Generalized Cross-Validation (GCV)
#' }
#'
#' The function provides extensive flexibility through different Laplacian constructions,
#' various spectral filter types, and parameter controls for localization and smoothness.
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph. Each
#'   element \code{adj.list[[i]]} contains the indices of vertices adjacent to vertex \code{i}.
#'   Uses 1-based indexing. Must be non-empty.
#' @param weight.list A list of numeric vectors where \code{weight.list[[i]]} contains the edge
#'   weights corresponding to the edges in \code{adj.list[[i]]}. Must have the same structure
#'   as \code{adj.list}.
#' @param y Numeric vector of length equal to the number of vertices containing the signal values
#'   to be filtered. NA values are not allowed.
#' @param laplacian.type Integer code specifying the type of graph Laplacian to use (default: 4).
#'   Valid values are 0-11. See Details section for descriptions of each type.
#' @param filter.type Integer code specifying the type of spectral filter to apply (default: 1).
#'   Valid values are 0-11. See Details section for descriptions of each type.
#' @param laplacian.power Positive integer specifying the power to which the Laplacian is raised
#'   (default: 3). Higher powers provide smoother filtering.
#' @param kernel.tau.factor Positive numeric value in (0,1] determining kernel bandwidth as a
#'   fraction of graph diameter (default: 0.05).
#' @param kernel.radius.factor Numeric value >= 1 multiplying the search radius when finding
#'   vertices within kernel range (default: 5.0).
#' @param kernel.type Integer code (0-8) specifying the kernel function type (default: 0 for
#'   Gaussian). See Details section.
#' @param kernel.adaptive Logical indicating whether to use locally adaptive kernel bandwidth
#'   (default: FALSE).
#' @param kernel.min.radius.factor Minimum radius factor in (0,1) as fraction of graph diameter
#'   (default: 0.05).
#' @param kernel.max.radius.factor Maximum radius factor in (0,1] as fraction of graph diameter
#'   (default: 0.99).
#' @param kernel.domain.min.size Positive integer specifying minimum number of vertices required
#'   within radius for kernel computations (default: 4).
#' @param kernel.precision Positive numeric value for radius determination precision (default: 1e-6).
#' @param n.evectors.to.compute Positive integer specifying number of Laplacian eigenpairs to
#'   compute. Default is \code{min(20, length(y)-2)}.
#' @param n.candidates Positive integer specifying number of filter parameter values to evaluate
#'   (default: 200).
#' @param log.grid Logical indicating whether to use logarithmically-spaced filter parameters
#'   (default: TRUE).
#' @param with.t.predictions Logical indicating whether to return smoothed signals for all
#'   parameter values (default: TRUE).
#' @param verbose Logical indicating whether to print progress information (default: FALSE).
#'
#' @section Laplacian Types:
#' \describe{
#'   \item{0 - STANDARD}{L = D - A, the combinatorial Laplacian. Provides basic smoothing that
#'     minimizes first differences across edges. Best for general-purpose smoothing on regularly
#'     structured graphs.}
#'   \item{1 - NORMALIZED}{L_norm = D^(-1/2) L D^(-1/2), the normalized Laplacian. Accounts for
#'     varying vertex degrees, giving more balanced smoothing on irregular graphs. Recommended
#'     for graphs with highly variable connectivity.}
#'   \item{2 - RANDOM_WALK}{L_rw = D^(-1) L, the random walk Laplacian. Similar to normalized
#'     Laplacian but with asymmetric normalization. Useful when modeling diffusion processes.}
#'   \item{3 - KERNEL}{L_kernel = D_kernel - W_kernel, a Laplacian constructed using distance-based
#'     kernel weights rather than adjacency weights. Provides smoother spectral response,
#'     especially on irregular or noisy graphs.}
#'   \item{4 - NORMALIZED_KERNEL}{Normalized version of the kernel Laplacian. Combines benefits
#'     of normalization and kernel-based edge weighting. Excellent for irregular graphs where
#'     geometric relationships matter.}
#'   \item{5 - ADAPTIVE_KERNEL}{Kernel Laplacian with locally adaptive bandwidth. Adjusts smoothing
#'     automatically based on local graph density. Best for graphs with highly variable density.}
#'   \item{6 - SHIFTED}{I - L, the shifted standard Laplacian. Inverts the spectrum to emphasize
#'     smooth components. Particularly effective for chain graphs and 1D signals.}
#'   \item{7 - SHIFTED_KERNEL}{I - L_kernel, the shifted kernel Laplacian. Combines kernel-based
#'     edge weighting with spectral inversion.}
#'   \item{8 - REGULARIZED}{L + \eqn{\epsilon}*I, adds small regularization to ensure positive definiteness.
#'     Helps with numerical stability in filtering operations.}
#'   \item{9 - REGULARIZED_KERNEL}{L_kernel + \eqn{\epsilon}*I, regularized version of kernel Laplacian.}
#'   \item{10 - MULTI_SCALE}{Weighted combination of kernel Laplacians at different scales.
#'     Captures both fine and coarse features simultaneously.}
#'   \item{11 - PATH}{Path Laplacians for specialized graph structures.}
#' }
#'
#' @section Filter Types:
#' \describe{
#'   \item{0 - HEAT}{\eqn{\exp(-t \lambda)}, classic heat kernel filter. Provides smooth decay across
#'     frequencies with more pronounced filtering at higher frequencies.}
#'   \item{1 - GAUSSIAN}{\eqn{\exp(-t \lambda^2)}), Gaussian spectral filter. More aggressive decay at higher
#'     frequencies than heat kernel. Produces very smooth results with minimal ringing.}
#'   \item{2 - NON_NEGATIVE}{\eqn{\exp(-t \max(\lambda,0))}, truncated heat kernel that only attenuates
#'     non-negative eigenvalues.}
#'   \item{3 - CUBIC_SPLINE}{\eqn{1/(1+t \lambda^2)}, filter that mimics cubic spline behavior. Minimizes
#'     second derivatives, producing the smoothest results while preserving linear trends.}
#'   \item{4 - EXPONENTIAL}{\eqn{\exp(-t \sqrt(\lambda))}, less aggressive decay than heat kernel.}
#'   \item{5 - MEXICAN_HAT}{\eqn{\lambda \exp(-t \lambda^2)}, band-pass filter that enhances mid-frequencies.}
#'   \item{6 - IDEAL_LOW_PASS}{1 for \eqn{\lambda < t}, 0 otherwise. Sharp cutoff filter.}
#'   \item{7 - BUTTERWORTH}{\eqn{1/(1+(\lambda/t)^(2n))}, smoother cutoff than ideal filter.}
#'   \item{8 - TIKHONOV}{\eqn{1/(1+t\lambda)}, first-order smoothing filter.}
#'   \item{9 - POLYNOMIAL}{\eqn{(1-\lambda/\lambda_{\max})^p} for \eqn{\lambda < \lambda_{\max}}, polynomial decay filter.}
#'   \item{10 - INVERSE_COSINE}{\eqn{\cos(\pi \lambda/(2\lambda_{\max}))}, smooth filter with cosine profile.}
#'   \item{11 - ADAPTIVE}{Data-driven filter that adapts to signal properties.}
#' }
#'
#' @section Kernel Types:
#' \describe{
#'   \item{0 - GAUSSIAN}{\eqn{\exp(-d^2/\tau^2)}, Classic bell curve, smooth decay from center}
#'   \item{1 - EXPONENTIAL}{\eqn{\exp(-d/\tau)}, Sharper peak, heavier tails than Gaussian}
#'   \item{2 - HEAT}{\eqn{\exp(\frac{-d^2}{4\tau})}, Similar to Gaussian but with different scaling}
#'   \item{3 - TRICUBE}{\eqn{(1-(d/\tau)^3)^3} for \eqn{d < \tau}, Compact support, smooth transition to zero}
#'   \item{4 - EPANECHNIKOV}{\eqn{1-(d/\tau)^2} for \eqn{d < \tau}, Parabolic shape, optimal in statistical sense}
#'   \item{5 - UNIFORM}{1 for \eqn{d < \tau}, 0 otherwise, Equal weighting within radius}
#'   \item{6 - TRIANGULAR}{\eqn{1-|d/\tau|} for \eqn{d < \tau}, Linear decay from center}
#'   \item{7 - QUARTIC}{\eqn{(1-(d/\tau)^2)^2} for \eqn{d < \tau}, Similar to Gaussian but with compact support}
#'   \item{8 - TRIWEIGHT}{\eqn{(1-(d/\tau)^2)^3} for \eqn{d < \tau}, Higher-order version of quartic}
#' }
#'
#' @return An object of class \code{"graph_spectral_filter"}, which is a list containing:
#' \describe{
#'   \item{evalues}{Numeric vector of eigenvalues of the Laplacian operator}
#'   \item{evectors}{Matrix of corresponding eigenvectors (columns)}
#'   \item{candidate_ts}{Numeric vector of filter parameter values tested}
#'   \item{gcv_scores}{Numeric vector of Generalized Cross-Validation scores}
#'   \item{opt_t_idx}{Integer index of the optimal parameter value}
#'   \item{predictions}{Numeric vector of smoothed signal at optimal parameter}
#'   \item{t_predictions}{Matrix of smoothed signals at each parameter (if requested)}
#'   \item{laplacian_type}{Factor indicating Laplacian type used}
#'   \item{filter_type}{Factor indicating filter type applied}
#'   \item{laplacian_power}{Integer power to which Laplacian was raised}
#'   \item{kernel_params}{List of kernel parameters used}
#'   \item{compute_time_ms}{Numeric computation time in milliseconds}
#'   \item{gcv_min_score}{Numeric minimum GCV score achieved}
#' }
#'
#' @references
#' Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., & Vandergheynst, P. (2013).
#' The emerging field of signal processing on graphs: Extending high-dimensional data
#' analysis to networks and other irregular domains. IEEE Signal Processing Magazine,
#' 30(3), 83-98.
#'
#' @seealso
#' \code{\link{plot.graph_spectral_filter}} for visualization,
#' \code{\link{predict.graph_spectral_filter}} for prediction
#'
#' @examples
#' \dontrun{
#' # Create a simple chain graph with 30 vertices (for CRAN check timing)
#' n <- 30
#' adj_list <- vector("list", n)
#' weight_list <- vector("list", n)
#'
#' # Build adjacency structure for chain
#' for (i in seq_len(n)) {
#'   adj_list[[i]] <- integer(0)
#'   weight_list[[i]] <- numeric(0)
#'
#'   if (i > 1) {
#'     adj_list[[i]] <- c(adj_list[[i]], i - 1L)
#'     weight_list[[i]] <- c(weight_list[[i]], 1.0)
#'   }
#'   if (i < n) {
#'     adj_list[[i]] <- c(adj_list[[i]], i + 1L)
#'     weight_list[[i]] <- c(weight_list[[i]], 1.0)
#'   }
#' }
#'
#' # Create a noisy signal
#' set.seed(123)
#' x <- seq(0, 1, length.out = n)
#' y <- sin(2 * pi * x) + rnorm(n, 0, 0.2)
#'
#' # Standard Laplacian with heat kernel filter
#' result1 <- graph.spectral.filter(
#'   adj.list = adj_list,
#'   weight.list = weight_list,
#'   y = y,
#'   laplacian.type = 0L,    # STANDARD
#'   filter.type = 0L,       # HEAT
#'   laplacian.power = 1L,
#'   n.candidates = 50       # Reduced for faster execution
#' )
#'
#' # Cubic spline-like smoothing
#' result2 <- graph.spectral.filter(
#'   adj.list = adj_list,
#'   weight.list = weight_list,
#'   y = y,
#'   laplacian.type = 0L,    # STANDARD
#'   filter.type = 3L,       # CUBIC_SPLINE
#'   laplacian.power = 2L,
#'   n.candidates = 50
#' )
#'
#' # Compare results
#' plot(x, y, pch = 16, col = "gray", main = "Graph Spectral Filtering",
#'      xlab = "Position", ylab = "Value")
#' lines(x, result1$predictions, col = "blue", lwd = 2)
#' lines(x, result2$predictions, col = "red", lwd = 2)
#' legend("topright",
#'        legend = c("Noisy data", "Heat kernel", "Cubic spline"),
#'        pch = c(16, NA, NA),
#'        lty = c(NA, 1, 1),
#'        col = c("gray", "blue", "red"),
#'        lwd = c(NA, 2, 2))
#'
#' # Print summary
#' summary(result1)
#' }
#'
#' @export
graph.spectral.filter <- function(adj.list,
                                  weight.list,
                                  y,
                                  laplacian.type = 4L,
                                  filter.type = 1L,
                                  laplacian.power = 3L,
                                  kernel.tau.factor = 0.05,
                                  kernel.radius.factor = 5.0,
                                  kernel.type = 0L,
                                  kernel.adaptive = FALSE,
                                  kernel.min.radius.factor = 0.05,
                                  kernel.max.radius.factor = 0.99,
                                  kernel.domain.min.size = 4L,
                                  kernel.precision = 1e-6,
                                  n.evectors.to.compute = min(c(20L, length(y) - 2L)),
                                  n.candidates = 200L,
                                  log.grid = TRUE,
                                  with.t.predictions = TRUE,
                                  verbose = FALSE) {

    # Input validation
    if (!is.list(adj.list) || length(adj.list) == 0L) {
        stop("'adj.list' must be a non-empty list of integer vectors")
    }

    if (!is.list(weight.list) || length(weight.list) != length(adj.list)) {
        stop("'weight.list' must be a list of same length as 'adj.list'")
    }

    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop("'y' must be numeric vector of length equal to number of vertices")
    }

    if (anyNA(y)) {
        stop("'y' contains NA values")
    }

    # Validate adjacency list structure
    for (i in seq_along(adj.list)) {
        if (!is.null(adj.list[[i]]) && length(adj.list[[i]]) > 0L) {
            if (!is.numeric(adj.list[[i]])) {
                stop(sprintf("'adj.list[[%d]]' must contain numeric values", i))
            }
            if (any(adj.list[[i]] < 1L) || any(adj.list[[i]] > length(adj.list))) {
                stop(sprintf("'adj.list[[%d]]' contains out-of-range vertex indices", i))
            }
            if (length(weight.list[[i]]) != length(adj.list[[i]])) {
                stop(sprintf("'weight.list[[%d]]' must have same length as 'adj.list[[%d]]'", i, i))
            }
            if (!is.numeric(weight.list[[i]])) {
                stop(sprintf("'weight.list[[%d]]' must contain numeric values", i))
            }
            if (any(weight.list[[i]] < 0)) {
                stop(sprintf("'weight.list[[%d]]' contains negative weights", i))
            }
        }
    }

    # Parameter validation with informative messages
    laplacian.type <- as.integer(laplacian.type)
    if (length(laplacian.type) != 1L || is.na(laplacian.type) ||
        laplacian.type < 0L || laplacian.type > 11L) {
        stop("'laplacian.type' must be an integer between 0 and 11")
    }

    filter.type <- as.integer(filter.type)
    if (length(filter.type) != 1L || is.na(filter.type) ||
        filter.type < 0L || filter.type > 11L) {
        stop("'filter.type' must be an integer between 0 and 11")
    }

    laplacian.power <- as.integer(laplacian.power)
    if (length(laplacian.power) != 1L || is.na(laplacian.power) ||
        laplacian.power < 1L) {
        stop("'laplacian.power' must be a positive integer")
    }

    kernel.tau.factor <- as.numeric(kernel.tau.factor)
    if (length(kernel.tau.factor) != 1L || is.na(kernel.tau.factor) ||
        kernel.tau.factor <= 0 || kernel.tau.factor > 1) {
        stop("'kernel.tau.factor' must be a number in (0, 1]")
    }

    kernel.radius.factor <- as.numeric(kernel.radius.factor)
    if (length(kernel.radius.factor) != 1L || is.na(kernel.radius.factor) ||
        kernel.radius.factor < 1) {
        stop("'kernel.radius.factor' must be >= 1")
    }

    kernel.type <- as.integer(kernel.type)
    if (length(kernel.type) != 1L || is.na(kernel.type) ||
        kernel.type < 0L || kernel.type > 8L) {
        stop("'kernel.type' must be an integer between 0 and 8")
    }

    if (!is.logical(kernel.adaptive) || length(kernel.adaptive) != 1L ||
        is.na(kernel.adaptive)) {
        stop("'kernel.adaptive' must be TRUE or FALSE")
    }

    kernel.min.radius.factor <- as.numeric(kernel.min.radius.factor)
    if (length(kernel.min.radius.factor) != 1L || is.na(kernel.min.radius.factor) ||
        kernel.min.radius.factor <= 0 || kernel.min.radius.factor >= 1) {
        stop("'kernel.min.radius.factor' must be a number in (0, 1)")
    }

    kernel.max.radius.factor <- as.numeric(kernel.max.radius.factor)
    if (length(kernel.max.radius.factor) != 1L || is.na(kernel.max.radius.factor) ||
        kernel.max.radius.factor <= 0 || kernel.max.radius.factor > 1) {
        stop("'kernel.max.radius.factor' must be a number in (0, 1]")
    }

    if (kernel.min.radius.factor >= kernel.max.radius.factor) {
        stop("'kernel.min.radius.factor' must be less than 'kernel.max.radius.factor'")
    }

    kernel.domain.min.size <- as.integer(kernel.domain.min.size)
    if (length(kernel.domain.min.size) != 1L || is.na(kernel.domain.min.size) ||
        kernel.domain.min.size < 1L) {
        stop("'kernel.domain.min.size' must be a positive integer")
    }

    kernel.precision <- as.numeric(kernel.precision)
    if (length(kernel.precision) != 1L || is.na(kernel.precision) ||
        kernel.precision <= 0) {
        stop("'kernel.precision' must be a positive number")
    }

    n.evectors.to.compute <- as.integer(n.evectors.to.compute)
    if (length(n.evectors.to.compute) != 1L || is.na(n.evectors.to.compute) ||
        n.evectors.to.compute < 1L || n.evectors.to.compute > length(y)) {
        stop("'n.evectors.to.compute' must be between 1 and number of vertices")
    }

    n.candidates <- as.integer(n.candidates)
    if (length(n.candidates) != 1L || is.na(n.candidates) || n.candidates < 1L) {
        stop("'n.candidates' must be a positive integer")
    }

    if (!is.logical(log.grid) || length(log.grid) != 1L || is.na(log.grid)) {
        stop("'log.grid' must be TRUE or FALSE")
    }

    if (!is.logical(with.t.predictions) || length(with.t.predictions) != 1L ||
        is.na(with.t.predictions)) {
        stop("'with.t.predictions' must be TRUE or FALSE")
    }

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("'verbose' must be TRUE or FALSE")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(v) {
        if (length(v) == 0L) {
            integer(0)
        } else {
            as.integer(v - 1L)
        }
    })

    # Call C++ implementation
    res <- .Call("S_graph_spectral_filter",
                 adj.list.0based,
                 weight.list,
                 as.numeric(y),
                 laplacian.type,
                 filter.type,
                 laplacian.power,
                 kernel.tau.factor,
                 kernel.radius.factor,
                 kernel.type,
                 kernel.adaptive,
                 kernel.min.radius.factor,
                 kernel.max.radius.factor,
                 kernel.domain.min.size,
                 kernel.precision,
                 n.evectors.to.compute,
                 n.candidates,
                 log.grid,
                 with.t.predictions,
                 verbose)
    
    # Define factor labels
    laplacian_names <- c("STANDARD", "NORMALIZED", "RANDOM_WALK", "KERNEL",
                         "NORMALIZED_KERNEL", "ADAPTIVE_KERNEL", "SHIFTED",
                         "SHIFTED_KERNEL", "REGULARIZED", "REGULARIZED_KERNEL",
                         "MULTI_SCALE", "PATH")

    filter_names <- c("HEAT", "GAUSSIAN", "NON_NEGATIVE", "CUBIC_SPLINE",
                      "EXPONENTIAL", "MEXICAN_HAT", "IDEAL_LOW_PASS",
                      "BUTTERWORTH", "TIKHONOV", "POLYNOMIAL",
                      "INVERSE_COSINE", "ADAPTIVE")

    kernel_names <- c("GAUSSIAN", "EXPONENTIAL", "HEAT", "TRICUBE",
                      "EPANECHNIKOV", "UNIFORM", "TRIANGULAR", "QUARTIC",
                      "TRIWEIGHT")

    # Convert to factors for better interpretability
    res$laplacian_type <- factor(res$laplacian_type,
                                  levels = 0:11,
                                  labels = laplacian_names)

    res$filter_type <- factor(res$filter_type,
                              levels = 0:11,
                              labels = filter_names)

    if (!is.null(res$kernel_params) && !is.null(res$kernel_params$kernel_type)) {
        res$kernel_params$kernel_type <- factor(res$kernel_params$kernel_type,
                                                 levels = 0:8,
                                                 labels = kernel_names)
    }

    # Add class
    class(res) <- c("graph_spectral_filter", "list")

    return(res)
}

#' Print Method for Graph Spectral Filter Results
#'
#' @param x An object of class \code{"graph_spectral_filter"}
#' @param digits Integer; number of significant digits to print (default: 4)
#' @param ... Further arguments passed to or from other methods
#'
#' @return Invisibly returns the input object
#'
#' @export
print.graph_spectral_filter <- function(x, digits = 4L, ...) {
    cat("Graph Spectral Filter Result\n")
    cat("============================\n")
    cat("  Vertices:                ", length(x$predictions), "\n")
    cat("  Eigenpairs computed:     ", length(x$evalues), "\n")
    cat("  Diffusion times tested:  ", length(x$candidate_ts), "\n")
    cat("  Laplacian type:          ", as.character(x$laplacian_type), "\n")
    cat("  Filter type:             ", as.character(x$filter_type), "\n")
    cat("  Laplacian power:         ", x$laplacian_power, "\n")

    if (!is.null(x$kernel_params) && !is.null(x$kernel_params$kernel_type)) {
        cat("  Kernel type:             ", as.character(x$kernel_params$kernel_type), "\n")
    }

    cat("  Optimal parameter index: ", x$opt_t_idx, "\n")
    cat("  Optimal parameter value: ",
        format(x$candidate_ts[x$opt_t_idx], digits = digits), "\n")
    cat("  Computation time:        ", format(x$compute_time_ms, digits = digits), " ms\n")
    cat("  GCV minimum score:       ", format(x$gcv_min_score, digits = digits), "\n")

    invisible(x)
}

#' Summary Method for Graph Spectral Filter Results
#'
#' @param object An object of class \code{"graph_spectral_filter"}
#' @param ... Further arguments passed to or from other methods
#'
#' @return An object of class \code{"summary.graph_spectral_filter"} containing
#'   key summary statistics
#'
#' @export
summary.graph_spectral_filter <- function(object, ...) {
    out <- list(
        n_vertices      = length(object$predictions),
        n_eigenpairs    = length(object$evalues),
        n_candidates    = length(object$candidate_ts),
        laplacian_type  = object$laplacian_type,
        filter_type     = object$filter_type,
        laplacian_power = object$laplacian_power,
        kernel_params   = object$kernel_params,
        opt_t_idx       = object$opt_t_idx,
        opt_t           = object$candidate_ts[object$opt_t_idx],
        compute_time_ms = object$compute_time_ms,
        gcv_min_score   = object$gcv_min_score,
        eigenvalue_range = range(object$evalues),
        gcv_score_range = range(object$gcv_scores)
    )

    class(out) <- "summary.graph_spectral_filter"
    out
}

#' Print Method for Summary of Graph Spectral Filter Results
#'
#' @param x An object of class \code{"summary.graph_spectral_filter"}
#' @param digits Integer; number of significant digits to print (default: 4)
#' @param ... Further arguments passed to or from other methods
#'
#' @return Invisibly returns the input object
#'
#' @export
print.summary.graph_spectral_filter <- function(x, digits = 4L, ...) {
    cat("Summary: Graph Spectral Filter Result\n")
    cat("=====================================\n")
    cat("Graph properties:\n")
    cat("  Number of vertices:      ", x$n_vertices, "\n")
    cat("  Eigenpairs computed:     ", x$n_eigenpairs, "\n")
    cat("\n")

    cat("Filter configuration:\n")
    cat("  Laplacian type:          ", as.character(x$laplacian_type), "\n")
    cat("  Filter type:             ", as.character(x$filter_type), "\n")
    cat("  Laplacian power:         ", x$laplacian_power, "\n")

    if (!is.null(x$kernel_params)) {
        cat("\nKernel parameters:\n")
        cat("  - tau_factor:            ",
            format(x$kernel_params$tau_factor, digits = digits), "\n")
        cat("  - radius_factor:         ",
            format(x$kernel_params$radius_factor, digits = digits), "\n")
        if (!is.null(x$kernel_params$kernel_type)) {
            cat("  - kernel_type:           ", as.character(x$kernel_params$kernel_type), "\n")
        }
        cat("  - adaptive:              ", x$kernel_params$adaptive, "\n")
    }

    cat("\nOptimization results:\n")
    cat("  Parameters tested:       ", x$n_candidates, "\n")
    cat("  Optimal parameter index: ", x$opt_t_idx, "\n")
    cat("  Optimal parameter value: ", format(x$opt_t, digits = digits), "\n")
    cat("  GCV minimum score:       ", format(x$gcv_min_score, digits = digits), "\n")
    cat("  GCV score range:         [",
        format(x$gcv_score_range[1], digits = digits), ", ",
        format(x$gcv_score_range[2], digits = digits), "]\n", sep = "")

    cat("\nSpectral information:\n")
    cat("  Eigenvalue range:        [",
        format(x$eigenvalue_range[1], digits = digits), ", ",
        format(x$eigenvalue_range[2], digits = digits), "]\n", sep = "")

    cat("\nComputation time:          ",
        format(x$compute_time_ms, digits = digits), " ms\n")

    invisible(x)
}

#' Plot Method for Graph Spectral Filter Results
#'
#' @param x An object of class \code{"graph_spectral_filter"}
#' @param which Integer or character vector specifying which plots to produce:
#'   1 = "gcv" (GCV scores), 2 = "spectrum" (eigenvalue spectrum),
#'   3 = "filter" (filter response), 4 = "signal" (original vs filtered signal)
#' @param main Main title for the plots
#' @param ... Further graphical parameters
#'
#' @return Invisibly returns the input object
#'
#' @export
plot.graph_spectral_filter <- function(x,
                                       which = c(1L, 4L),
                                       main = NULL,
                                       ...) {
    if (!inherits(x, "graph_spectral_filter")) {
        stop("'x' must be of class 'graph_spectral_filter'")
    }

    # Handle which parameter
    show <- rep(FALSE, 4L)
    if (is.numeric(which)) {
        show[which] <- TRUE
    } else if (is.character(which)) {
        show[match(which, c("gcv", "spectrum", "filter", "signal"), nomatch = 0L)] <- TRUE
    }

    # Plot 1: GCV scores
    if (show[1L] && length(x$gcv_scores) > 0L) {
        plot(x$candidate_ts, x$gcv_scores,
             type = "l",
             xlab = "Filter parameter (t)",
             ylab = "GCV score",
             main = if (is.null(main)) "GCV Scores vs Filter Parameter" else main[1L],
             ...)
        abline(v = x$candidate_ts[x$opt_t_idx], col = "red", lty = 2)
        points(x$candidate_ts[x$opt_t_idx], x$gcv_scores[x$opt_t_idx],
               col = "red", pch = 19)
    }

    # Plot 2: Eigenvalue spectrum
    if (show[2L] && length(x$evalues) > 0L) {
        plot(seq_along(x$evalues), x$evalues,
             type = "h",
             xlab = "Eigenvalue index",
             ylab = "Eigenvalue",
             main = if (is.null(main)) "Laplacian Eigenvalue Spectrum" else
                    if (length(main) >= 2) main[2L] else "Eigenvalue Spectrum",
             ...)
    }

    # Plot 3: Filter response
    if (show[3L] && length(x$evalues) > 0L) {
        t_opt <- x$candidate_ts[x$opt_t_idx]
        lambda_range <- range(x$evalues)
        lambda_seq <- seq(lambda_range[1], lambda_range[2], length.out = 100)

        # Compute filter response based on filter type
        if (x$filter_type == "HEAT") {
            response <- exp(-t_opt * lambda_seq)
        } else if (x$filter_type == "GAUSSIAN") {
            response <- exp(-t_opt * lambda_seq^2)
        } else if (x$filter_type == "CUBIC_SPLINE") {
            response <- 1 / (1 + t_opt * lambda_seq^2)
        } else if (x$filter_type == "TIKHONOV") {
            response <- 1 / (1 + t_opt * lambda_seq)
        } else {
            # Default to heat kernel
            response <- exp(-t_opt * lambda_seq)
        }

        plot(lambda_seq, response,
             type = "l",
             xlab = "Eigenvalue",
             ylab = "Filter response",
             main = if (is.null(main)) paste("Filter Response:", x$filter_type) else
                    if (length(main) >= 3) main[3L] else "Filter Response",
             ylim = c(0, 1),
             ...)
    }

    # Plot 4: Original vs filtered signal
    if (show[4L] && length(x$predictions) > 0L) {
        n <- length(x$predictions)
        plot(seq_len(n), x$predictions,
             type = "l",
             col = "blue",
             lwd = 2,
             xlab = "Vertex index",
             ylab = "Signal value",
             main = if (is.null(main)) "Filtered Signal" else
                    if (length(main) >= 4) main[4L] else "Filtered Signal",
             ...)
    }

    invisible(x)
}

#' Predict Method for Graph Spectral Filter
#'
#' @param object An object of class \code{"graph_spectral_filter"}
#' @param newdata Optional numeric vector of new signal values. If NULL, returns
#'   the filtered signal from the original fit
#' @param t Optional filter parameter value. If NULL, uses the optimal value
#'   from the original fit
#' @param ... Further arguments (currently ignored)
#'
#' @return Numeric vector of filtered signal values
#'
#' @export
predict.graph_spectral_filter <- function(object, newdata = NULL, t = NULL, ...) {
    if (!inherits(object, "graph_spectral_filter")) {
        stop("'object' must be of class 'graph_spectral_filter'")
    }

    # If no new data, return original predictions
    if (is.null(newdata)) {
        return(object$predictions)
    }

    # Validate new data
    if (!is.numeric(newdata) || length(newdata) != ncol(object$evectors)) {
        stop("'newdata' must be numeric vector of same length as original signal")
    }

    # Use optimal t if not specified
    if (is.null(t)) {
        t <- object$candidate_ts[object$opt_t_idx]
    } else {
        t <- as.numeric(t)
        if (length(t) != 1L || is.na(t) || t <= 0) {
            stop("'t' must be a positive number")
        }
    }

    # Apply spectral filtering
    # Project onto eigenbasis
    coeffs <- crossprod(object$evectors, newdata)

    # Apply filter based on type
    if (object$filter_type == "HEAT") {
        filter_weights <- exp(-t * object$evalues)
    } else if (object$filter_type == "GAUSSIAN") {
        filter_weights <- exp(-t * object$evalues^2)
    } else if (object$filter_type == "CUBIC_SPLINE") {
        filter_weights <- 1 / (1 + t * object$evalues^2)
    } else if (object$filter_type == "TIKHONOV") {
        filter_weights <- 1 / (1 + t * object$evalues)
    } else {
        # Default to heat kernel
        filter_weights <- exp(-t * object$evalues)
    }

    # Reconstruct
    filtered_coeffs <- coeffs * filter_weights
    predictions <- object$evectors %*% filtered_coeffs

    return(as.vector(predictions))
}

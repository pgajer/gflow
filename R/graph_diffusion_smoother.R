#' Extended Graph Diffusion Smoother
#'
#' @description
#' Performs graph diffusion smoothing on a graph signal using an iterative diffusion
#' process. This method is useful for denoising graph signals, interpolating missing
#' values, and enhancing patterns in network data while preserving the underlying
#' graph structure.
#'
#' @param graph A list of integer vectors where each element contains the indices
#'   of neighbors for the corresponding vertex (1-based indexing).
#' @param edge.lengths A list of numeric vectors where each element contains the
#'   edge lengths corresponding to the neighbors in \code{graph}.
#' @param y A numeric vector of signal values at each vertex.
#' @param weights An optional numeric vector of weights for each vertex. If
#'   \code{NULL}, all vertices are weighted equally. Default is \code{NULL}.
#' @param n.time.steps An integer specifying the number of diffusion steps.
#'   Default is 100.
#' @param step.factor A numeric value between 0 and 1 controlling the magnitude
#'   of each diffusion step. Smaller values lead to more gradual smoothing.
#'   Default is 0.1.
#' @param normalize An integer (0, 1, or 2) specifying the normalization method:
#'   \itemize{
#'     \item 0: No normalization
#'     \item 1: Range adjustment to \eqn{[0, 1]}
#'     \item 2: Mean adjustment to preserve the mean of the signal
#'   }
#'   Default is 0.
#' @param preserve.local.maxima A logical value. If \code{TRUE}, local maxima are
#'   preserved during diffusion by adjusting weights based on neighborhood values.
#'   Default is \code{FALSE}.
#' @param local.maximum.weight.factor A numeric value between 0 and 1 controlling
#'   the degree of preservation for local maxima. Default is 1.0.
#' @param preserve.local.extrema A logical value. If \code{TRUE}, both local
#'   maxima and minima are preserved. Cannot be \code{TRUE} simultaneously with
#'   \code{preserve.local.maxima}. Default is \code{FALSE}.
#' @param imputation.method An integer (0-4) specifying the imputation method:
#'   \itemize{
#'     \item 0: LOCAL_MEAN_THRESHOLD
#'     \item 1: NEIGHBORHOOD_MATCHING
#'     \item 2: ITERATIVE_NEIGHBORHOOD_MATCHING
#'     \item 3: SUPPLIED_THRESHOLD
#'     \item 4: GLOBAL_MEAN_THRESHOLD
#'   }
#'   Default is 0.
#' @param max.iterations A positive integer for the maximum number of iterations
#'   in iterative imputation methods. Default is 10.
#' @param convergence.threshold A positive numeric value for convergence criteria.
#'   Default is 1e-6.
#' @param apply.binary.threshold A logical value indicating whether to apply
#'   binary thresholding to the output. Default is \code{TRUE}.
#' @param binary.threshold A numeric value between 0 and 1 for binary thresholding.
#'   Default is 0.5.
#' @param kernel An integer (0-4) specifying the kernel type for distance weighting:
#'   \itemize{
#'     \item 0: Uniform kernel
#'     \item 1: Gaussian kernel
#'     \item 2: Laplacian kernel
#'     \item 3: Exponential kernel
#'     \item 4: Cauchy kernel
#'   }
#'   Default is 1.
#' @param dist.normalization.factor A numeric value >= 1 for distance normalization.
#'   Default is 1.01.
#' @param n.CVs A non-negative integer for the number of cross-validation runs.
#'   If 0, no cross-validation is performed. Default is 0.
#' @param n.CV.folds An integer > 1 for the number of folds in cross-validation.
#'   Default is 10.
#' @param epsilon A positive numeric value for numerical stability. Default is 1e-10.
#' @param seed An integer for the random number generator seed. Default is 0.
#' @param verbose A logical value. If \code{TRUE}, prints progress information.
#'   Default is \code{TRUE}.
#'
#' @return A list of class \code{"ext.gds"} containing:
#'   \item{y.traj}{A matrix where each column represents the signal at a time step.}
#'   \item{y.optimal}{The optimally smoothed signal (based on CV if performed).}
#'   \item{cv.errors}{A matrix of cross-validation errors (rows: time steps,
#'     columns: CV runs).}
#'   \item{mean.cv.errors}{Mean CV errors across runs for each time step.}
#'   \item{Rmean.cv.errors}{R-computed mean CV errors (for verification).}
#'   \item{Rmedian.cv.errors}{Median CV errors across runs.}
#'   \item{optimal.time.step}{The time step with minimum CV error.}
#'   \item{min.cv.error}{The minimum mean CV error.}
#'
#' @details
#' The graph diffusion process iteratively updates vertex values based on their
#' neighbors' values, weighted by edge lengths and kernel functions. The algorithm
#' can preserve local extrema to maintain important features while smoothing noise.
#'
#' Cross-validation uses a masking approach where test vertices don't influence
#' the diffusion but their values are predicted to evaluate performance.
#'
#' @examples
#' \dontrun{
#' # Create a simple chain graph
#' n <- 20
#' graph <- vector("list", n)
#' edge.lengths <- vector("list", n)
#' for(i in 1:(n-1)) {
#'   graph[[i]] <- c(graph[[i]], i+1)
#'   graph[[i+1]] <- c(graph[[i+1]], i)
#'   edge.lengths[[i]] <- c(edge.lengths[[i]], 1)
#'   edge.lengths[[i+1]] <- c(edge.lengths[[i+1]], 1)
#' }
#'
#' # Create noisy signal
#' true.signal <- sin(seq(0, 2*pi, length.out = n))
#' y <- true.signal + rnorm(n, 0, 0.2)
#'
#' # Apply diffusion smoothing
#' result <- ext.graph.diffusion.smoother(
#'   graph = graph,
#'   edge.lengths = edge.lengths,
#'   y = y,
#'   n.time.steps = 50,
#'   step.factor = 0.1,
#'   n.CVs = 5,
#'   verbose = TRUE
#' )
#'
#' # Plot results
#' plot(y, type = "b", col = "red", main = "Graph Diffusion Smoothing")
#' lines(result$y.optimal, type = "b", col = "blue")
#' legend("topright", c("Original", "Smoothed"), col = c("red", "blue"), lty = 1)
#' }
#'
#' @references
#' Coifman, R. R., & Lafon, S. (2006). Diffusion maps. Applied and computational
#' harmonic analysis, 21(1), 5-30.
#'
#' @seealso \code{\link{graph.diffusion.matrix.smoother}},
#'   \code{\link{instrumented.gds}}
#'
#' @export
ext.graph.diffusion.smoother <- function(graph,
                                        edge.lengths,
                                        y,
                                        weights = NULL,
                                        n.time.steps = 100,
                                        step.factor = 0.1,
                                        normalize = 0,
                                        preserve.local.maxima = FALSE,
                                        local.maximum.weight.factor = 1.0,
                                        preserve.local.extrema = FALSE,
                                        imputation.method = 0,
                                        max.iterations = 10,
                                        convergence.threshold = 1e-6,
                                        apply.binary.threshold = TRUE,
                                        binary.threshold = 0.5,
                                        kernel = 1,
                                        dist.normalization.factor = 1.01,
                                        n.CVs = 0,
                                        n.CV.folds = 10,
                                        epsilon = 1e-10,
                                        seed = 0,
                                        verbose = TRUE) {

    # Input validation with informative error messages
    if (!is.list(graph)) {
        stop("'graph' must be a list", call. = FALSE)
    }

    if (!all(vapply(graph, is.numeric, logical(1)))) {
        stop("All elements of 'graph' must be numeric vectors", call. = FALSE)
    }

    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list", call. = FALSE)
    }

    if (!all(vapply(edge.lengths, is.numeric, logical(1)))) {
        stop("All elements of 'edge.lengths' must be numeric vectors", call. = FALSE)
    }

    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    n.vertices <- length(y)

    if (n.vertices == 0) {
        stop("'y' must not be empty", call. = FALSE)
    }

    if (length(graph) != n.vertices) {
        stop("Length of 'graph' must equal length of 'y'", call. = FALSE)
    }

    if (length(edge.lengths) != n.vertices) {
        stop("Length of 'edge.lengths' must equal length of 'y'", call. = FALSE)
    }

    # Validate graph structure
    for (i in seq_along(graph)) {
        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop(sprintf("Mismatch between graph[[%d]] and edge.lengths[[%d]]", i, i),
                 call. = FALSE)
        }

        if (length(graph[[i]]) > 0) {
            if (any(graph[[i]] < 1) || any(graph[[i]] > n.vertices)) {
                stop(sprintf("Invalid vertex indices in graph[[%d]]", i), call. = FALSE)
            }

            if (any(!is.finite(graph[[i]]))) {
                stop(sprintf("Non-finite values in graph[[%d]]", i), call. = FALSE)
            }

            if (any(!is.finite(edge.lengths[[i]]))) {
                stop(sprintf("Non-finite values in edge.lengths[[%d]]", i), call. = FALSE)
            }

            if (any(edge.lengths[[i]] <= 0)) {
                stop(sprintf("Non-positive edge lengths in edge.lengths[[%d]]", i),
                     call. = FALSE)
            }
        }
    }

    # Validate weights
    if (is.null(weights)) {
        weights <- rep(1.0, n.vertices)
    } else {
        if (!is.numeric(weights)) {
            stop("'weights' must be numeric", call. = FALSE)
        }

        if (length(weights) != n.vertices) {
            stop("Length of 'weights' must equal length of 'y'", call. = FALSE)
        }

        if (any(!is.finite(weights))) {
            stop("'weights' contains non-finite values", call. = FALSE)
        }

        if (any(weights < 0)) {
            stop("'weights' must be non-negative", call. = FALSE)
        }
    }

    # Validate numeric parameters
    validate_positive_integer <- function(x, name, min_val = 1) {
        if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
            stop(sprintf("'%s' must be a single finite number", name), call. = FALSE)
        }
        if (x < min_val || round(x) != x) {
            stop(sprintf("'%s' must be an integer >= %d", name, min_val), call. = FALSE)
        }
    }

    validate_numeric_range <- function(x, name, min_val, max_val,
                                       inclusive_min = TRUE, inclusive_max = TRUE) {
        if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
            stop(sprintf("'%s' must be a single finite number", name), call. = FALSE)
        }

        min_ok <- if (inclusive_min) x >= min_val else x > min_val
        max_ok <- if (inclusive_max) x <= max_val else x < max_val

        if (!min_ok || !max_ok) {
            min_sym <- if (inclusive_min) "[" else "("
            max_sym <- if (inclusive_max) "]" else ")"
            stop(sprintf("'%s' must be in range %s%g, %g%s",
                         name, min_sym, min_val, max_val, max_sym), call. = FALSE)
        }
    }

    validate_positive_integer(n.time.steps, "n.time.steps")
    validate_numeric_range(step.factor, "step.factor", 0, 1,
                           inclusive_min = FALSE, inclusive_max = FALSE)

    if (!is.numeric(normalize) || length(normalize) != 1 ||
        !(normalize %in% 0:2)) {
        stop("'normalize' must be 0, 1, or 2", call. = FALSE)
    }

    if (!is.logical(preserve.local.maxima) || length(preserve.local.maxima) != 1) {
        stop("'preserve.local.maxima' must be TRUE or FALSE", call. = FALSE)
    }

    if (!is.logical(preserve.local.extrema) || length(preserve.local.extrema) != 1) {
        stop("'preserve.local.extrema' must be TRUE or FALSE", call. = FALSE)
    }

    if (preserve.local.maxima && preserve.local.extrema) {
        stop("'preserve.local.maxima' and 'preserve.local.extrema' cannot both be TRUE",
             call. = FALSE)
    }

    validate_numeric_range(local.maximum.weight.factor,
                           "local.maximum.weight.factor", 0, 1)

    if (!is.numeric(imputation.method) || length(imputation.method) != 1 ||
        !(imputation.method %in% 0:4)) {
        stop("'imputation.method' must be an integer between 0 and 4", call. = FALSE)
    }

    validate_positive_integer(max.iterations, "max.iterations")
    validate_numeric_range(convergence.threshold, "convergence.threshold",
                           0, Inf, inclusive_min = FALSE)

    if (!is.logical(apply.binary.threshold) || length(apply.binary.threshold) != 1) {
        stop("'apply.binary.threshold' must be TRUE or FALSE", call. = FALSE)
    }

    validate_numeric_range(binary.threshold, "binary.threshold", 0, 1)

    if (!is.numeric(kernel) || length(kernel) != 1 || !(kernel %in% 0:4)) {
        stop("'kernel' must be an integer between 0 and 4", call. = FALSE)
    }

    validate_numeric_range(dist.normalization.factor,
                           "dist.normalization.factor", 1, Inf)
    validate_positive_integer(n.CVs, "n.CVs", min_val = 0)
    validate_positive_integer(n.CV.folds, "n.CV.folds", min_val = 2)
    validate_numeric_range(epsilon, "epsilon", 0, Inf, inclusive_min = FALSE)

    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
        stop("'seed' must be a finite number", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be TRUE or FALSE", call. = FALSE)
    }

    # Check for valid y values
    if (any(!is.finite(y))) {
        stop("'y' contains non-finite values", call. = FALSE)
    }

    # Convert to appropriate types
    n.time.steps <- as.integer(n.time.steps)
    normalize <- as.integer(normalize)
    imputation.method <- as.integer(imputation.method)
    max.iterations <- as.integer(max.iterations)
    kernel <- as.integer(kernel)
    n.CVs <- as.integer(n.CVs)
    n.CV.folds <- as.integer(n.CV.folds)
    seed <- as.integer(seed)

    # Convert graph to 0-based indexing for C++
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    # Set n.cores to 1 (sequential) as specified in the original
    n.cores <- 1L

    res <- .Call("S_ext_graph_diffusion_smoother",
                 graph.0based,
                 edge.lengths,
                 weights,
                 y,
                 n.time.steps,
                 as.numeric(step.factor),
                 normalize,
                 preserve.local.maxima,
                 as.numeric(local.maximum.weight.factor),
                 preserve.local.extrema,
                 imputation.method,
                 max.iterations,
                 as.numeric(convergence.threshold),
                 apply.binary.threshold,
                 as.numeric(binary.threshold),
                 kernel,
                 as.numeric(dist.normalization.factor),
                 n.CVs,
                 n.CV.folds,
                 as.numeric(epsilon),
                 seed,
                 n.cores,
                 verbose)

    # Post-process results
    result <- list(
        y.traj = res$y_traj,
        y.optimal = res$y_optimal,
        cv.errors = res$cv_errors,
        mean.cv.errors = res$mean_cv_errors,
        Rmean.cv.errors = if (n.CVs > 0) {
            apply(res$cv_errors, 1, mean, na.rm = TRUE)
        } else {
            numeric(0)
        },
        Rmedian.cv.errors = if (n.CVs > 0) {
            apply(res$cv_errors, 1, median, na.rm = TRUE)
        } else {
            numeric(0)
        },
        optimal.time.step = res$optimal_time_step,
        min.cv.error = res$min_cv_error
    )

    class(result) <- c("ext.gds", "list")
    result
}

#' Instrumented Graph Diffusion Smoother
#'
#' @description
#' Performs graph diffusion smoothing with comprehensive instrumentation and 
#' adaptive step sizes. This function collects detailed performance metrics
#' and allows fine-tuned control over the diffusion process.
#'
#' @param graph List of integer vectors representing adjacency structure.
#' @param edge.lengths List of numeric vectors containing edge lengths.
#' @param y Numeric vector of initial vertex values.
#' @param y.true Numeric vector of ground truth values for evaluation.
#' @param n.time.steps Integer number of diffusion steps.
#' @param base.step.factor Initial step size for all vertices.
#' @param use.pure.laplacian Logical; if \code{TRUE}, uses uniform weights.
#'   Default is \code{FALSE}.
#' @param ikernel Integer kernel type (1-3). Default is 1.
#' @param kernel.scale Scale parameter for kernels. Default is 1.0.
#' @param dist.normalization.factor Not used in current implementation.
#' @param increase.factor Multiplier for increasing step size. Default is 1.1.
#' @param decrease.factor Multiplier for decreasing step size. Default is 0.8.
#' @param oscillation.factor Multiplier when oscillating. Default is 0.5.
#' @param min.step Minimum allowed step size. Default is 0.01.
#' @param max.step Maximum allowed step size. Default is 2.0.
#'
#' @return An object of class \code{"instrumented.gds"} containing detailed
#'   metrics and trajectories.
#'
#' @examples
#' \dontrun{
#' # Example with simple graph
#' graph <- list(c(2), c(1,3), c(2))
#' edge.lengths <- list(c(1), c(1,1), c(1))
#' y.true <- c(0, 1, 0)
#' y.noisy <- y.true + rnorm(3, 0, 0.1)
#'
#' result <- instrumented.gds(
#'   graph = graph,
#'   edge.lengths = edge.lengths,
#'   y = y.noisy,
#'   y.true = y.true,
#'   n.time.steps = 50,
#'   base.step.factor = 0.5
#' )
#' }
#'
#' @export
instrumented.gds <- function(graph,
                            edge.lengths,
                            y,
                            y.true,
                            n.time.steps,
                            base.step.factor,
                            use.pure.laplacian = FALSE,
                            ikernel = 1,
                            kernel.scale = 1.0,
                            dist.normalization.factor = 1.01,
                            increase.factor = 1.1,
                            decrease.factor = 0.8,
                            oscillation.factor = 0.5,
                            min.step = 0.01,
                            max.step = 2.0) {

    # Type validation
    if (!is.list(graph)) {
        stop("'graph' must be a list", call. = FALSE)
    }

    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list", call. = FALSE)
    }

    if (!is.numeric(y)) {
        stop("'y' must be numeric", call. = FALSE)
    }

    if (!is.numeric(y.true)) {
        stop("'y.true' must be numeric", call. = FALSE)
    }

    # Length checks
    n.vertices <- length(y)

    if (n.vertices == 0) {
        stop("'y' must not be empty", call. = FALSE)
    }

    if (length(y.true) != n.vertices) {
        stop("'y' and 'y.true' must have the same length", call. = FALSE)
    }

    if (length(graph) != n.vertices) {
        stop("'graph' length must match 'y' length", call. = FALSE)
    }

    if (length(edge.lengths) != n.vertices) {
        stop("'edge.lengths' length must match 'y' length", call. = FALSE)
    }

    # Graph structure validation
    for (i in seq_along(graph)) {
        if (!is.numeric(graph[[i]])) {
            stop(sprintf("graph[[%d]] must be numeric", i), call. = FALSE)
        }

        if (!is.numeric(edge.lengths[[i]])) {
            stop(sprintf("edge.lengths[[%d]] must be numeric", i), call. = FALSE)
        }

        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop(sprintf("Mismatch between graph[[%d]] and edge.lengths[[%d]]", i, i),
                 call. = FALSE)
        }

        if (length(graph[[i]]) > 0) {
            if (any(graph[[i]] < 1 | graph[[i]] > n.vertices)) {
                stop(sprintf("Invalid vertex indices in graph[[%d]]", i), call. = FALSE)
            }

            if (any(edge.lengths[[i]] <= 0)) {
                stop(sprintf("Edge lengths must be positive at vertex %d", i),
                     call. = FALSE)
            }

            if (any(duplicated(graph[[i]]))) {
                stop(sprintf("Duplicate neighbors found in graph[[%d]]", i),
                     call. = FALSE)
            }
        }
    }

    # Parameter validation
    if (!is.numeric(n.time.steps) || length(n.time.steps) != 1 ||
        n.time.steps <= 0 || round(n.time.steps) != n.time.steps) {
        stop("'n.time.steps' must be a positive integer", call. = FALSE)
    }

    if (!is.numeric(base.step.factor) || length(base.step.factor) != 1 ||
        base.step.factor <= 0) {
        stop("'base.step.factor' must be a positive number", call. = FALSE)
    }

    if (!is.logical(use.pure.laplacian) || length(use.pure.laplacian) != 1) {
        stop("'use.pure.laplacian' must be TRUE or FALSE", call. = FALSE)
    }

    if (!is.numeric(ikernel) || length(ikernel) != 1 || !(ikernel %in% 1:3)) {
        stop("'ikernel' must be 1, 2, or 3", call. = FALSE)
    }

    if (!is.numeric(kernel.scale) || length(kernel.scale) != 1 ||
        kernel.scale <= 0) {
        stop("'kernel.scale' must be a positive number", call. = FALSE)
    }

    if (!is.numeric(dist.normalization.factor) ||
        length(dist.normalization.factor) != 1 ||
        dist.normalization.factor <= 1) {
        stop("'dist.normalization.factor' must be greater than 1", call. = FALSE)
    }

    # Step size factor validation
    if (!is.numeric(increase.factor) || length(increase.factor) != 1 ||
        increase.factor <= 1) {
        stop("'increase.factor' must be greater than 1", call. = FALSE)
    }

    if (!is.numeric(decrease.factor) || length(decrease.factor) != 1 ||
        decrease.factor <= 0 || decrease.factor >= 1) {
        stop("'decrease.factor' must be between 0 and 1", call. = FALSE)
    }

    if (!is.numeric(oscillation.factor) || length(oscillation.factor) != 1 ||
        oscillation.factor <= 0 || oscillation.factor >= 1) {
        stop("'oscillation.factor' must be between 0 and 1", call. = FALSE)
    }

    if (!is.numeric(min.step) || length(min.step) != 1 || min.step <= 0) {
        stop("'min.step' must be a positive number", call. = FALSE)
    }

    if (!is.numeric(max.step) || length(max.step) != 1 || max.step <= min.step) {
        stop("'max.step' must be greater than 'min.step'", call. = FALSE)
    }

    # Check for non-finite values
    if (any(!is.finite(y))) {
        stop("'y' contains non-finite values", call. = FALSE)
    }

    if (any(!is.finite(y.true))) {
        stop("'y.true' contains non-finite values", call. = FALSE)
    }

    # Convert graph to 0-based indexing
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    # Call C++ implementation
    result <- .Call("S_instrumented_gds",
                    graph.0based,
                    edge.lengths,
                    y,
                    y.true,
                    as.integer(n.time.steps),
                    as.double(base.step.factor),
                    use.pure.laplacian,
                    as.integer(ikernel),
                    as.double(kernel.scale),
                    as.double(increase.factor),
                    as.double(decrease.factor),
                    as.double(oscillation.factor),
                    as.double(min.step),
                    as.double(max.step))
    
    # Add class
    class(result) <- c("instrumented.gds", "list")
    result
}

#' Plot Instrumented GDS Results
#'
#' @description
#' Creates diagnostic plots for instrumented graph diffusion smoother results
#' using base R graphics. Supports both multi-panel and single-panel layouts.
#'
#' @param x An object of class \code{"instrumented.gds"}.
#' @param metrics Character vector of metrics to plot. Options include "snr",
#'   "mad", "energy.ratio". Default is c("snr", "mad", "energy.ratio").
#' @param layout Character string specifying the layout type:
#'   \itemize{
#'     \item "multi": Creates separate panels for each metric (default)
#'     \item "single": Plots all metrics in a single panel
#'   }
#' @param normalize Logical; if TRUE and layout="single", normalizes all metrics
#'   to \eqn{[0,1]} range for comparison. Default is FALSE.
#' @param main Main title for the plot. Default depends on layout.
#' @param xlab Label for x-axis. Default is "Iteration" for single layout.
#' @param ylab Label for y-axis. Default depends on layout and normalization.
#' @param col Colors for each metric. Can be a single color or a vector of colors.
#'   Default uses a color palette based on the number of metrics.
#' @param lty Line types for each metric. Default is 1 (solid) for multi-panel,
#'   or different line types for single-panel.
#' @param lwd Line width. Default is 2.
#' @param legend.pos Position of legend for single layout. Options include
#'   "topright", "topleft", "bottomright", "bottomleft", "top", "bottom",
#'   "left", "right", "center", or "none". Default is "topright".
#' @param cex.main Character expansion for main title. Default is 1.2.
#' @param cex.lab Character expansion for axis labels. Default is 1.
#' @param cex.axis Character expansion for axis text. Default is 0.9.
#' @param mar.panel Margins for each panel (multi layout only).
#'   Default is c(4, 4, 2, 1).
#' @param oma Outer margins (multi layout only). Default is c(0, 0, 3, 0).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns NULL. Called for side effect of creating plots.
#'
#' @details
#' This function supports two layout modes:
#' \itemize{
#'   \item Multi-panel layout: Creates vertically stacked panels with one panel
#'     per metric, each with its own y-axis scale.
#'   \item Single-panel layout: Plots all metrics on the same axes, optionally
#'     normalized to facilitate comparison.
#' }
#'
#' The choice of layout depends on your analysis needs:
#' \itemize{
#'   \item Use "multi" when metrics have very different scales or when you want
#'     to see detailed patterns in each metric.
#'   \item Use "single" when you want to compare the timing and relative changes
#'     across metrics.
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' graph <- list(c(2), c(1,3), c(2))
#' edge.lengths <- list(c(1), c(1,1), c(1))
#' y.true <- c(0, 1, 0)
#' y.noisy <- y.true + rnorm(3, 0, 0.1)
#'
#' # Run instrumented GDS
#' result <- instrumented.gds(
#'   graph = graph,
#'   edge.lengths = edge.lengths,
#'   y = y.noisy,
#'   y.true = y.true,
#'   n.time.steps = 50,
#'   base.step.factor = 0.5
#' )
#'
#' # Default multi-panel plot
#' plot(result)
#'
#' # Single panel with all metrics
#' plot(result, layout = "single")
#'
#' # Single panel with normalization
#' plot(result, layout = "single", normalize = TRUE)
#'
#' # Custom colors for multi-panel
#' plot(result, col = c("blue", "red", "darkgreen"), lwd = 3)
#'
#' # Plot only specific metrics in single panel
#' plot(result, metrics = c("snr", "mad"), layout = "single",
#'      col = c("navy", "darkred"), lty = c(1, 2))
#' }
#'
#' @seealso \code{\link{instrumented.gds}}
#'
#' @export
plot.instrumented.gds <- function(x,
                                  metrics = c("snr", "mad", "energy.ratio"),
                                  layout = c("multi", "single"),
                                  normalize = FALSE,
                                  main = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  col = NULL,
                                  lty = NULL,
                                  lwd = 2,
                                  legend.pos = "topright",
                                  cex.main = 1.2,
                                  cex.lab = 1,
                                  cex.axis = 0.9,
                                  mar.panel = c(4, 4, 2, 1),
                                  oma = c(0, 0, 3, 0),
                                  ...) {

    # Check that x is of correct class
    if (!inherits(x, "instrumented.gds")) {
        stop("'x' must be an object of class 'instrumented.gds'", call. = FALSE)
    }

    # Validate layout
    layout <- match.arg(layout)

    # Validate metrics
    valid_metrics <- c("snr", "mad", "energy.ratio")
    metrics <- match.arg(metrics, valid_metrics, several.ok = TRUE)
    n.metrics <- length(metrics)

    if (n.metrics == 0) {
        stop("At least one metric must be selected", call. = FALSE)
    }

    # Extract data
    iterations <- seq_along(x$snr_trajectory) - 1

    # Prepare data list
    data_list <- list()
    metric_labels <- list()
    metric_names <- character(0)

    if ("snr" %in% metrics) {
        if (is.null(x$snr_trajectory)) {
            warning("'snr_trajectory' not found in x", call. = FALSE)
        } else {
            data_list[["snr"]] <- x$snr_trajectory
            metric_labels[["snr"]] <- "SNR"
            metric_names <- c(metric_names, "snr")
        }
    }

    if ("mad" %in% metrics) {
        if (is.null(x$mean_absolute_deviation)) {
            warning("'mean_absolute_deviation' not found in x", call. = FALSE)
        } else {
            data_list[["mad"]] <- x$mean_absolute_deviation
            metric_labels[["mad"]] <- "Mean Absolute Deviation"
            metric_names <- c(metric_names, "mad")
        }
    }

    if ("energy.ratio" %in% metrics) {
        if (is.null(x$energy_ratio)) {
            warning("'energy_ratio' not found in x", call. = FALSE)
        } else {
            data_list[["energy.ratio"]] <- x$energy_ratio
            metric_labels[["energy.ratio"]] <- "Energy Ratio"
            metric_names <- c(metric_names, "energy.ratio")
        }
    }

    # Check if we have any data to plot
    if (length(data_list) == 0) {
        stop("No valid metric data found in x", call. = FALSE)
    }

    # Update number of metrics based on available data
    n.metrics <- length(data_list)

    # Set default titles based on layout
    if (is.null(main)) {
        main <- "Graph Diffusion Smoother Performance Metrics"
    }

    # Branch based on layout
    if (layout == "multi") {
        .plot_multi_panel(iterations, data_list, metric_labels, metric_names,
                         main, col, lty, lwd, cex.main, cex.lab, cex.axis,
                         mar.panel, oma, ...)
    } else {
        .plot_single_panel(iterations, data_list, metric_labels, metric_names,
                          normalize, main, xlab, ylab, col, lty, lwd,
                          legend.pos, cex.main, cex.lab, cex.axis, ...)
    }

    invisible(NULL)
}


# Internal function for multi-panel plot
.plot_multi_panel <- function(iterations, data_list, metric_labels, metric_names,
                             main, col, lty, lwd, cex.main, cex.lab, cex.axis,
                             mar.panel, oma, ...) {

    n.metrics <- length(data_list)

    # Set up colors and line types if not provided
    if (is.null(col)) {
        col <- c("blue", "red", "darkgreen", "purple", "orange", "brown")[seq_len(n.metrics)]
    } else {
        col <- rep_len(col, n.metrics)
    }

    if (is.null(lty)) {
        lty <- rep(1, n.metrics)
    } else {
        lty <- rep_len(lty, n.metrics)
    }

    # Save current par settings
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    # Set up multi-panel layout
    par(mfrow = c(n.metrics, 1), mar = mar.panel, oma = oma)

    # Plot each metric
    for (i in seq_along(metric_names)) {
        metric <- metric_names[i]
        y_data <- data_list[[metric]]
        y_label <- metric_labels[[metric]]

        # Create plot
        plot(iterations, y_data,
             type = "l",
             col = col[i],
             lty = lty[i],
             lwd = lwd,
             xlab = if (i == n.metrics) "Iteration" else "",
             ylab = y_label,
             cex.lab = cex.lab,
             cex.axis = cex.axis,
             las = 1,
             ...)

        # Add grid
        grid(col = "gray90", lty = "dotted")

        # Redraw line on top of grid
        lines(iterations, y_data,
              col = col[i],
              lty = lty[i],
              lwd = lwd)

        # Add metric name in top-left corner
        mtext(y_label, side = 3, line = 0.5, adj = 0, cex = cex.lab * 0.9, font = 2)
    }

    # Add main title
    mtext(main, outer = TRUE, cex = cex.main, line = 1, font = 2)
}


# Internal function for single-panel plot
.plot_single_panel <- function(iterations, data_list, metric_labels, metric_names,
                              normalize, main, xlab, ylab, col, lty, lwd,
                              legend.pos, cex.main, cex.lab, cex.axis, ...) {

    n.metrics <- length(data_list)

    # Set defaults for single panel
    if (is.null(xlab)) {
        xlab <- "Iteration"
    }

    if (is.null(ylab)) {
        ylab <- if (normalize) "Value (normalized)" else "Value"
    }

    # Set up colors and line types
    if (is.null(col)) {
        col <- rainbow(n.metrics)
    } else {
        col <- rep_len(col, n.metrics)
    }

    if (is.null(lty)) {
        lty <- seq_len(n.metrics)  # Different line types for distinction
    } else {
        lty <- rep_len(lty, n.metrics)
    }

    # Get the label names for the legend
    legend_labels <- unlist(metric_labels[metric_names])

    # Normalize data if requested
    plot_data <- data_list
    if (normalize) {
        plot_data <- lapply(data_list, function(y) {
            rng <- range(y, na.rm = TRUE)
            if (rng[2] - rng[1] > .Machine$double.eps) {
                (y - rng[1]) / (rng[2] - rng[1])
            } else {
                rep(0.5, length(y))
            }
        })
    }

    # Find y-axis range
    y_range <- range(unlist(plot_data), na.rm = TRUE)

    # Create base plot
    plot(iterations, plot_data[[1]],
         type = "n",
         xlim = range(iterations),
         ylim = y_range,
         xlab = xlab,
         ylab = ylab,
         main = main,
         cex.main = cex.main,
         cex.lab = cex.lab,
         cex.axis = cex.axis,
         las = 1,
         ...)

    # Add grid
    grid(col = "gray90", lty = "dotted")

    # Plot each metric
    for (i in seq_along(plot_data)) {
        lines(iterations, plot_data[[i]],
              col = col[i],
              lty = lty[i],
              lwd = lwd)
    }

    # Add legend
    if (!is.null(legend.pos) && legend.pos != "none") {
        legend(legend.pos,
               legend = legend_labels,
               col = col,
               lty = lty,
               lwd = lwd,
               bg = "white",
               box.lty = 1,
               box.col = "gray70",
               cex = cex.lab * 0.9)
    }
}

#' Graph Diffusion Smoother
#'
#' @description
#' A simplified interface for graph diffusion smoothing with cross-validation.
#' This function provides a more accessible API while maintaining full functionality.
#'
#' @param adj.list A list where each element contains neighbor indices.
#' @param weight.list A list where each element contains edge weights.
#' @param y A numeric vector of values to smooth.
#' @param n.time.steps Number of diffusion steps. Default is 100.
#' @param step.factor Step size factor (0, 1). Default is 0.5.
#' @param binary.threshold Threshold for binary classification. Default is 0.5.
#' @param ikernel Kernel type (0-4). Default is 1.
#' @param dist.normalization.factor Distance normalization > 1. Default is 1.1.
#' @param n.CVs Number of CV runs. Default is 0.
#' @param n.CV.folds Number of CV folds. Default is 10.
#' @param epsilon Numerical tolerance. Default is 1e-10.
#' @param verbose Print progress information. Default is FALSE.
#' @param seed Random seed. Default is 0.
#'
#' @return A list containing smoothing results and CV information.
#'
#' @examples
#' \dontrun{
#' # Create simple chain graph
#' adj.list <- list(2, c(1,3), c(2,4), c(3,5), 4)
#' weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
#' y <- c(1, 0, 1, 0, 1) + rnorm(5, 0, 0.1)
#'
#' result <- graph.diffusion.smoother(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   n.time.steps = 20,
#'   n.CVs = 3,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
graph.diffusion.smoother <- function(adj.list,
                                    weight.list,
                                    y,
                                    n.time.steps,
                                    step.factor,
                                    binary.threshold = 0.5,
                                    ikernel = 1,
                                    dist.normalization.factor = 1.1,
                                    n.CVs = 0,
                                    n.CV.folds = 10,
                                    epsilon = 1e-10,
                                    verbose = FALSE,
                                    seed = 0) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
    }

    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    n.vertices <- length(y)

    if (length(adj.list) != n.vertices) {
        stop("'adj.list' must have the same length as 'y'", call. = FALSE)
    }

    # Validate graph structure
    for (i in seq_along(adj.list)) {
        if (!is.numeric(adj.list[[i]])) {
            stop(sprintf("adj.list[[%d]] must be numeric", i), call. = FALSE)
        }

        if (!is.numeric(weight.list[[i]])) {
            stop(sprintf("weight.list[[%d]] must be numeric", i), call. = FALSE)
        }

        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("adj.list[[%d]] and weight.list[[%d]] must have same length",
                         i, i), call. = FALSE)
        }

        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1 | adj.list[[i]] > n.vertices)) {
                stop(sprintf("Invalid indices in adj.list[[%d]]", i), call. = FALSE)
            }

            if (any(weight.list[[i]] <= 0)) {
                stop(sprintf("Weights must be positive in weight.list[[%d]]", i),
                     call. = FALSE)
            }
        }
    }

    # Validate numeric parameters
    if (!is.numeric(n.time.steps) || n.time.steps < 1) {
        stop("'n.time.steps' must be a positive integer", call. = FALSE)
    }

    if (!is.numeric(step.factor) || step.factor <= 0 || step.factor >= 1) {
        stop("'step.factor' must be in (0, 1)", call. = FALSE)
    }

    if (!is.numeric(binary.threshold)) {
        stop("'binary.threshold' must be numeric", call. = FALSE)
    }

    if (!is.numeric(ikernel) || !(ikernel %in% 0:4)) {
        stop("'ikernel' must be an integer between 0 and 4", call. = FALSE)
    }

    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 1.0) {
        stop("'dist.normalization.factor' must be greater than 1.0", call. = FALSE)
    }

    if (!is.numeric(n.CVs) || n.CVs < 0) {
        stop("'n.CVs' must be a non-negative integer", call. = FALSE)
    }

    if (!is.numeric(n.CV.folds) || n.CV.folds < 2) {
        stop("'n.CV.folds' must be an integer >= 2", call. = FALSE)
    }

    if (!is.numeric(epsilon) || epsilon <= 0) {
        stop("'epsilon' must be positive", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be TRUE or FALSE", call. = FALSE)
    }

    if (!is.numeric(seed)) {
        stop("'seed' must be numeric", call. = FALSE)
    }

    # Convert to integers
    n.time.steps <- as.integer(n.time.steps)
    ikernel <- as.integer(ikernel)
    n.CVs <- as.integer(n.CVs)
    n.CV.folds <- as.integer(n.CV.folds)
    seed <- as.integer(seed)

    # Convert to 0-based indexing
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ implementation
    result <- .Call(
        "S_graph_diffusion_smoother",
        adj.list.0based,
        weight.list,
        y,
        n.time.steps,
        step.factor,
        binary.threshold,
        ikernel,
        dist.normalization.factor,
        n.CVs,
        n.CV.folds,
        epsilon,
        verbose,
        seed)

    # Convert trajectory to matrix
    if (!is.null(result$y_traj)) {
        result$y_traj <- as.matrix(as.data.frame(result$y_traj))
        colnames(result$y_traj) <- NULL
    }

    result
}

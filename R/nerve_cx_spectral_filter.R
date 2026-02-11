#' Apply Spectral Filtering to a Nerve Complex
#'
#' @description
#' Performs spectral filtering on a signal defined on the vertices of a nerve complex.
#' This method extends graph spectral filtering to incorporate higher-order relationships
#' captured by the simplicial structure of the nerve complex. The filtering uses a full
#' Laplacian that combines information from all dimensions of simplices, weighted according
#' to the specified dimension weights.
#'
#' @param complex A nerve complex object created by \code{\link{create.nerve.complex}}.
#'   Must contain a valid complex pointer and dimension information.
#' @param y Numeric vector of function values at vertices of the complex. Length must
#'   equal the number of vertices in the complex.
#' @param laplacian.type Character string specifying the type of Laplacian. One of:
#'   \describe{
#'     \item{\code{"STANDARD"}}{Standard combinatorial Laplacian (default)}
#'     \item{\code{"NORMALIZED"}}{Normalized Laplacian}
#'     \item{\code{"RANDOM_WALK"}}{Random walk Laplacian}
#'     \item{\code{"SHIFTED"}}{Shifted Laplacian (I - L)}
#'     \item{\code{"REGULARIZED"}}{Regularized Laplacian (L + epsilon*I)}
#'   }
#' @param filter.type Character string specifying the type of spectral filter. One of:
#'   \describe{
#'     \item{\code{"HEAT"}}{Heat kernel filter \eqn{\exp(-t\lambda)} (default)}
#'     \item{\code{"GAUSSIAN"}}{Gaussian filter \eqn{\exp(-t\lambda^2)}}
#'     \item{\code{"NON_NEGATIVE"}}{Non-negative heat kernel filter}
#'     \item{\code{"CUBIC_SPLINE"}}{Cubic spline filter \eqn{1/(1+t\lambda^2)}}
#'     \item{\code{"EXPONENTIAL"}}{Exponential filter}
#'     \item{\code{"MEXICAN_HAT"}}{Mexican hat wavelet filter}
#'     \item{\code{"IDEAL_LOW_PASS"}}{Ideal low-pass filter}
#'     \item{\code{"BUTTERWORTH"}}{Butterworth filter}
#'     \item{\code{"TIKHONOV"}}{Tikhonov filter}
#'     \item{\code{"POLYNOMIAL"}}{Polynomial filter}
#'     \item{\code{"INVERSE_COSINE"}}{Inverse cosine filter}
#'     \item{\code{"ADAPTIVE"}}{Adaptive filter}
#'   }
#' @param laplacian.power Positive integer specifying the power to which the Laplacian
#'   is raised. Default is 1. Higher powers produce smoother results.
#' @param dim.weights Numeric vector of non-negative weights for each dimension's
#'   contribution to the full Laplacian. Length should be at least max.dimension + 1.
#'   Default is a vector of 1s. The i-th element weights the i-1 dimensional simplices.
#' @param kernel.params Named list of kernel parameters used by certain Laplacian types:
#'   \describe{
#'     \item{\code{tau.factor}}{Factor for kernel bandwidth (0 < tau.factor <= 1)}
#'     \item{\code{radius.factor}}{Factor for neighborhood radius (>= 1)}
#'     \item{\code{kernel.type}}{Integer code for kernel type (0-8)}
#'   }
#' @param n.evectors Positive integer specifying the number of eigenvectors to compute.
#'   Default is 100. Set to 0 to compute all eigenvectors (may be slow for large complexes).
#' @param n.candidates Positive integer specifying the number of filter parameter values
#'   to evaluate for automatic selection via GCV. Default is 100.
#' @param log.grid Logical indicating whether to use logarithmic spacing for the parameter
#'   grid (TRUE, default) or linear spacing (FALSE).
#' @param with.t.predictions Logical indicating whether to return predictions for all
#'   tested parameter values (TRUE) or only the optimal value (FALSE, default).
#' @param verbose Logical indicating whether to print progress information during
#'   computation. Default is FALSE.
#'
#' @return An object of class \code{nerve_cx_spectral_filter} containing:
#'   \item{predictions}{Numeric vector of smoothed function values at the optimal parameter}
#'   \item{optimal_parameter}{The selected optimal filter parameter value}
#'   \item{gcv_score}{The GCV score at the optimal parameter}
#'   \item{all_parameters}{Vector of all tested parameter values}
#'   \item{all_gcv_scores}{Vector of GCV scores for each parameter}
#'   \item{compute_time_ms}{Computation time in milliseconds}
#'   \item{method}{List containing the input method parameters}
#'   \item{t_predictions}{(Optional) Matrix of predictions for all parameter values if
#'     \code{with.t.predictions = TRUE}}
#'
#' @details
#' The function implements a comprehensive framework for spectral filtering on nerve complexes.
#' The general process involves:
#' \enumerate{
#'   \item Construction of a weighted Laplacian operator incorporating all simplex dimensions
#'   \item Eigendecomposition of the Laplacian
#'   \item Application of the specified spectral filter in the frequency domain
#'   \item Automatic parameter selection using Generalized Cross-Validation (GCV)
#' }
#'
#' The dimension weights allow flexible control over how much each simplex dimension
#' contributes to the overall smoothing. Setting higher weights for higher dimensions
#' incorporates more global structure into the filtering.
#'
#' @section Parameter Selection:
#' The optimal filter parameter is selected automatically using leave-one-out
#' Generalized Cross-Validation (GCV). This avoids overfitting and typically produces
#' good results without manual tuning.
#'
#' @seealso
#' \code{\link{create.nerve.complex}} for creating nerve complexes,
#' \code{\link{set.complex.function.values}} for setting function values,
#' \code{\link{plot.nerve_cx_spectral_filter}} for visualization
#'
#' @examples
#' \dontrun{
#' # Generate 2D points
#' set.seed(123)
#' coords <- matrix(runif(200), ncol = 2)
#'
#' # Create a smooth function with noise
#' f_true <- function(x) sin(2*pi*x[1]) * cos(2*pi*x[2])
#' y_true <- apply(coords, 1, f_true)
#' y_noisy <- y_true + rnorm(length(y_true), 0, 0.2)
#'
#' # Create nerve complex
#' complex <- create.nerve.complex(coords, k = 8, max.dim = 2)
#' complex <- set.complex.function.values(complex, y_noisy)
#'
#' # Apply spectral filtering with default parameters
#' result <- nerve.cx.spectral.filter(complex, y_noisy)
#'
#' # Apply with custom parameters emphasizing higher dimensions
#' result2 <- nerve.cx.spectral.filter(
#'   complex, y_noisy,
#'   laplacian.type = "NORMALIZED",
#'   filter.type = "GAUSSIAN",
#'   dim.weights = c(1.0, 0.5, 0.25),
#'   n.candidates = 150,
#'   verbose = TRUE
#' )
#'
#' # Compare results
#' mse1 <- mean((result$predictions - y_true)^2)
#' mse2 <- mean((result2$predictions - y_true)^2)
#' cat("MSE (default):", mse1, "\n")
#' cat("MSE (custom):", mse2, "\n")
#' }
#'
#' @export
nerve.cx.spectral.filter <- function(complex,
                                   y,
                                   laplacian.type = "STANDARD",
                                   filter.type = "HEAT",
                                   laplacian.power = 1,
                                   dim.weights = NULL,
                                   kernel.params = list(tau.factor = 0.01,
                                                      radius.factor = 3.0,
                                                      kernel.type = 0),
                                   n.evectors = 100,
                                   n.candidates = 100,
                                   log.grid = TRUE,
                                   with.t.predictions = FALSE,
                                   verbose = FALSE) {

    # Input validation
    if (!inherits(complex, "nerve_complex")) {
        stop("'complex' must be a nerve_complex object", call. = FALSE)
    }

    if (!is.numeric(y) || anyNA(y)) {
        stop("'y' must be a numeric vector without NA values", call. = FALSE)
    }

    if (length(y) != complex$n_vertices) {
        stop(sprintf("Length of 'y' (%d) must match number of vertices in complex (%d)",
                     length(y), complex$n_vertices), call. = FALSE)
    }

    # Validate laplacian.type
    laplacian.type <- match.arg(laplacian.type,
                               c("STANDARD", "NORMALIZED", "RANDOM_WALK",
                                 "SHIFTED", "REGULARIZED"))

    # Validate filter.type
    filter.type <- match.arg(filter.type,
                           c("HEAT", "GAUSSIAN", "NON_NEGATIVE", "CUBIC_SPLINE",
                             "EXPONENTIAL", "MEXICAN_HAT", "IDEAL_LOW_PASS",
                             "BUTTERWORTH", "TIKHONOV", "POLYNOMIAL",
                             "INVERSE_COSINE", "ADAPTIVE"))

    # Validate numeric parameters
    if (!is.numeric(laplacian.power) || length(laplacian.power) != 1 ||
        laplacian.power < 1 || laplacian.power != as.integer(laplacian.power)) {
        stop("'laplacian.power' must be a positive integer", call. = FALSE)
    }
    laplacian.power <- as.integer(laplacian.power)

    # Set default dimension weights if not provided
    if (is.null(dim.weights)) {
        dim.weights <- rep(1.0, complex$max_dimension + 1)
    } else {
        if (!is.numeric(dim.weights) || any(dim.weights < 0) || anyNA(dim.weights)) {
            stop("'dim.weights' must be a numeric vector of non-negative values", call. = FALSE)
        }
        if (length(dim.weights) < complex$max_dimension + 1) {
            stop(sprintf("'dim.weights' must have at least %d elements",
                        complex$max_dimension + 1), call. = FALSE)
        }
    }

    # Validate kernel parameters
    if (!is.list(kernel.params)) {
        stop("'kernel.params' must be a list", call. = FALSE)
    }

    # Set defaults for kernel parameters
    default_kernel_params <- list(tau.factor = 0.01, radius.factor = 3.0, kernel.type = 0)
    kernel.params <- utils::modifyList(default_kernel_params, kernel.params)

    # Validate individual kernel parameters
    if (!is.numeric(kernel.params$tau.factor) || length(kernel.params$tau.factor) != 1 ||
        kernel.params$tau.factor <= 0 || kernel.params$tau.factor > 1) {
        stop("'kernel.params$tau.factor' must be a single value in (0, 1]", call. = FALSE)
    }

    if (!is.numeric(kernel.params$radius.factor) || length(kernel.params$radius.factor) != 1 ||
        kernel.params$radius.factor < 1) {
        stop("'kernel.params$radius.factor' must be a single value >= 1", call. = FALSE)
    }

    if (!is.numeric(kernel.params$kernel.type) || length(kernel.params$kernel.type) != 1 ||
        kernel.params$kernel.type < 0 || kernel.params$kernel.type > 8 ||
        kernel.params$kernel.type != as.integer(kernel.params$kernel.type)) {
        stop("'kernel.params$kernel.type' must be an integer between 0 and 8", call. = FALSE)
    }

    # Validate other parameters
    if (!is.numeric(n.evectors) || length(n.evectors) != 1 || n.evectors < 0 ||
        n.evectors != as.integer(n.evectors)) {
        stop("'n.evectors' must be a non-negative integer", call. = FALSE)
    }
    n.evectors <- as.integer(n.evectors)

    if (!is.numeric(n.candidates) || length(n.candidates) != 1 || n.candidates < 1 ||
        n.candidates != as.integer(n.candidates)) {
        stop("'n.candidates' must be a positive integer", call. = FALSE)
    }
    n.candidates <- as.integer(n.candidates)

    if (!is.logical(log.grid) || length(log.grid) != 1) {
        stop("'log.grid' must be TRUE or FALSE", call. = FALSE)
    }

    if (!is.logical(with.t.predictions) || length(with.t.predictions) != 1) {
        stop("'with.t.predictions' must be TRUE or FALSE", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be TRUE or FALSE", call. = FALSE)
    }

    # Map string types to integer codes for C++ function
    laplacian.type.map <- c(
        "STANDARD" = 0,
        "NORMALIZED" = 1,
        "RANDOM_WALK" = 2,
        "SHIFTED" = 3,
        "REGULARIZED" = 4
    )

    filter.type.map <- c(
        "HEAT" = 0,
        "GAUSSIAN" = 1,
        "NON_NEGATIVE" = 2,
        "CUBIC_SPLINE" = 3,
        "EXPONENTIAL" = 4,
        "MEXICAN_HAT" = 5,
        "IDEAL_LOW_PASS" = 6,
        "BUTTERWORTH" = 7,
        "TIKHONOV" = 8,
        "POLYNOMIAL" = 9,
        "INVERSE_COSINE" = 10,
        "ADAPTIVE" = 11
    )

    # Transform kernel.params names to match C++ function
    kernel.params.c <- list(
        tau_factor = kernel.params$tau.factor,
        radius_factor = kernel.params$radius.factor,
        kernel_type = as.integer(kernel.params$kernel.type)
    )

    # Call the C++ function
    if (verbose) {
        message("Starting nerve complex spectral filtering...")
        message(sprintf("  Laplacian type: %s", laplacian.type))
        message(sprintf("  Filter type: %s", filter.type))
        message(sprintf("  Number of candidates: %d", n.candidates))
    }

    result <- .Call("S_nerve_cx_spectral_filter",
                    complex$complex_ptr,
                    as.numeric(y),
                    as.integer(laplacian.type.map[laplacian.type]),
                    as.integer(filter.type.map[filter.type]),
                    laplacian.power,
                    as.numeric(dim.weights),
                    kernel.params.c,
                    n.evectors,
                    n.candidates,
                    log.grid,
                    with.t.predictions,
                    verbose)

    # Add method information to result
    result$method <- list(
        laplacian.type = laplacian.type,
        filter.type = filter.type,
        laplacian.power = laplacian.power,
        dim.weights = dim.weights,
        kernel.params = kernel.params,
        n.evectors = n.evectors,
        n.candidates = n.candidates,
        log.grid = log.grid
    )

    # Set class
    class(result) <- "nerve_cx_spectral_filter"

    if (verbose) {
        message(sprintf("Filtering complete. Optimal parameter: %.4e", result$optimal_parameter))
    }

    return(result)
}


#' Plot Nerve Complex Spectral Filter Results
#'
#' @description
#' Creates visualizations for the results of \code{\link{nerve.cx.spectral.filter}}.
#' Supports multiple plot types including parameter selection curves, prediction
#' comparisons, and residual analysis.
#'
#' @param x A \code{nerve_cx_spectral_filter} object returned by
#'   \code{\link{nerve.cx.spectral.filter}}.
#' @param type Character string specifying the plot type. One of:
#'   \describe{
#'     \item{\code{"parameters"}}{GCV scores vs. filter parameter values (default)}
#'     \item{\code{"predictions"}}{Original vs. smoothed values scatter plot}
#'     \item{\code{"residuals"}}{Residuals vs. original values}
#'   }
#' @param y Original function values. Required for \code{"predictions"} and
#'   \code{"residuals"} plot types. Must have the same length as the predictions.
#' @param main Character string giving the main title for the plot. If \code{NULL}
#'   (default), an appropriate title is chosen based on the plot type.
#' @param xlab Character string giving the x-axis label. If \code{NULL} (default),
#'   an appropriate label is chosen based on the plot type.
#' @param ylab Character string giving the y-axis label. If \code{NULL} (default),
#'   an appropriate label is chosen based on the plot type.
#' @param col Color specification for the plot. Default is "black" for most plots,
#'   with the optimal parameter highlighted in red for the parameters plot.
#' @param pch Plotting character. Default is 1 (open circles).
#' @param ... Additional graphical parameters passed to the underlying plot functions.
#'
#' @return No return value. Function is called for its side effect of creating a plot.
#'
#' @details
#' The \code{"parameters"} plot shows the GCV score as a function of the filter
#' parameter on a log-log scale, with the optimal parameter marked in red. This
#' helps visualize the parameter selection process and assess whether the minimum
#' is well-defined.
#'
#' The \code{"predictions"} plot shows a scatter plot of original vs. smoothed values
#' with a diagonal reference line. Points close to the diagonal indicate good agreement.
#'
#' The \code{"residuals"} plot shows residuals (original - smoothed) vs. original values
#' with a horizontal reference line at zero. This helps identify systematic patterns
#' in the filtering errors.
#'
#' @seealso
#' \code{\link{nerve.cx.spectral.filter}} for the main filtering function,
#' \code{\link{summary.nerve_cx_spectral_filter}} for numerical summaries
#'
#' @examples
#' \dontrun{
#' # Create example data and apply filtering
#' coords <- matrix(runif(200), ncol = 2)
#' complex <- create.nerve.complex(coords, k = 5, max.dim = 2)
#' y <- sin(coords[,1] * 2*pi) * cos(coords[,2] * 2*pi) + rnorm(100, 0, 0.1)
#' complex <- set.complex.function.values(complex, y)
#' result <- nerve.cx.spectral.filter(complex, y)
#'
#' # Create different plot types
#' par(mfrow = c(1, 3))
#' plot(result, type = "parameters")
#' plot(result, type = "predictions", y = y)
#' plot(result, type = "residuals", y = y)
#' }
#'
#' @method plot nerve_cx_spectral_filter
#' @export
plot.nerve_cx_spectral_filter <- function(x, type = c("parameters", "predictions", "residuals"),
                                        y = NULL, main = NULL, xlab = NULL, ylab = NULL,
                                        col = "black", pch = 1, ...) {
    # Match plot type
    type <- match.arg(type)

    # Check if y is required and provided
    if (type %in% c("predictions", "residuals")) {
        if (is.null(y)) {
            stop(sprintf("Original function values 'y' must be provided for plot type '%s'",
                        type), call. = FALSE)
        }
        if (!is.numeric(y) || length(y) != length(x$predictions)) {
            stop("'y' must be a numeric vector with the same length as predictions",
                 call. = FALSE)
        }
    }

    # Create the requested plot
    switch(type,
        parameters = {
            # Set default labels if not provided
            if (is.null(main)) main <- "GCV Score vs. Filter Parameter"
            if (is.null(xlab)) xlab <- "Filter Parameter (log scale)"
            if (is.null(ylab)) ylab <- "GCV Score (log scale)"

            # Create the plot
            plot(x$all_parameters, x$all_gcv_scores,
                 log = "xy", type = "b", col = col, pch = pch,
                 main = main, xlab = xlab, ylab = ylab, ...)

            # Mark optimal parameter
            opt_idx <- which.min(x$all_gcv_scores)
            points(x$all_parameters[opt_idx], x$all_gcv_scores[opt_idx],
                   col = "red", pch = 19, cex = 1.5)
            abline(v = x$optimal_parameter, col = "red", lty = 2)

            # Add legend
            legend("topright",
                   legend = c("GCV scores", "Optimal"),
                   col = c(col, "red"),
                   pch = c(pch, 19),
                   lty = c(1, NA),
                   bty = "n")
        },

        predictions = {
            # Set default labels if not provided
            if (is.null(main)) main <- "Original vs. Smoothed Values"
            if (is.null(xlab)) xlab <- "Original Values"
            if (is.null(ylab)) ylab <- "Smoothed Values"

            # Create scatter plot
            plot(y, x$predictions, col = col, pch = pch,
                 main = main, xlab = xlab, ylab = ylab, ...)

            # Add diagonal reference line
            abline(0, 1, col = "red", lty = 2)

            # Add correlation info
            cor_val <- cor(y, x$predictions)
            legend("bottomright",
                   legend = sprintf("Correlation: %.3f", cor_val),
                   bty = "n")
        },

        residuals = {
            # Set default labels if not provided
            if (is.null(main)) main <- "Residual Plot"
            if (is.null(xlab)) xlab <- "Original Values"
            if (is.null(ylab)) ylab <- "Residuals (Original - Smoothed)"

            # Calculate residuals
            residuals <- y - x$predictions

            # Create residual plot
            plot(y, residuals, col = col, pch = pch,
                 main = main, xlab = xlab, ylab = ylab, ...)

            # Add horizontal line at zero
            abline(h = 0, col = "red", lty = 2)

            # Add summary statistics
            rmse <- sqrt(mean(residuals^2))
            legend("topright",
                   legend = c(sprintf("RMSE: %.3f", rmse),
                            sprintf("Mean: %.3e", mean(residuals))),
                   bty = "n")
        }
    )

    invisible(NULL)
}


#' Print Method for nerve_cx_spectral_filter Objects
#'
#' @description
#' Prints a concise summary of nerve complex spectral filtering results.
#'
#' @param x A \code{nerve_cx_spectral_filter} object returned by
#'   \code{\link{nerve.cx.spectral.filter}}.
#' @param digits Integer indicating the number of significant digits to print
#'   for numeric values. Default is 4.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} is returned invisibly.
#'
#' @seealso
#' \code{\link{summary.nerve_cx_spectral_filter}} for more detailed summary,
#' \code{\link{nerve.cx.spectral.filter}} for the main function
#'
#' @method print nerve_cx_spectral_filter
#' @export
print.nerve_cx_spectral_filter <- function(x, digits = 4, ...) {
    cat("Nerve Complex Spectral Filter Results\n")
    cat("-------------------------------------\n")
    cat("Laplacian type:      ", x$method$laplacian.type, "\n")
    cat("Filter type:         ", x$method$filter.type, "\n")
    cat("Laplacian power:     ", x$method$laplacian.power, "\n")
    cat("Optimal parameter:   ", format(x$optimal_parameter, digits = digits), "\n")
    cat("GCV score:           ", format(x$gcv_score, digits = digits + 2), "\n")
    cat("Computation time:    ", format(x$compute_time_ms, digits = digits), "ms\n")
    cat("Parameters tested:   ", length(x$all_parameters), "\n")

    # Show dimension weights if not all equal to 1
    if (!all(x$method$dim.weights == 1)) {
        cat("Dimension weights:   ",
            paste(format(x$method$dim.weights, digits = 2), collapse = ", "), "\n")
    }

    cat("\nUse summary() for detailed statistics\n")
    cat("Use plot() to visualize results\n")

    invisible(x)
}


#' Summary Method for nerve_cx_spectral_filter Objects
#'
#' @description
#' Provides a comprehensive summary of nerve complex spectral filtering results,
#' including statistics about the predictions, parameter optimization, and
#' computational details.
#'
#' @param object A \code{nerve_cx_spectral_filter} object returned by
#'   \code{\link{nerve.cx.spectral.filter}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{summary.nerve_cx_spectral_filter} containing:
#'   \item{method}{List of method parameters used}
#'   \item{optimal_parameter}{The selected optimal filter parameter}
#'   \item{gcv_score}{The GCV score at the optimal parameter}
#'   \item{compute_time_ms}{Computation time in milliseconds}
#'   \item{predictions}{Summary statistics of the smoothed values}
#'   \item{parameters}{Summary of the parameter search}
#'   \item{improvement}{GCV improvement from worst to best parameter}
#'
#' @seealso
#' \code{\link{print.summary.nerve_cx_spectral_filter}} for printing,
#' \code{\link{nerve.cx.spectral.filter}} for the main function
#'
#' @method summary nerve_cx_spectral_filter
#' @export
summary.nerve_cx_spectral_filter <- function(object, ...) {
    # Calculate summary statistics for predictions
    pred_stats <- list(
        min = min(object$predictions),
        q1 = quantile(object$predictions, 0.25),
        median = median(object$predictions),
        mean = mean(object$predictions),
        q3 = quantile(object$predictions, 0.75),
        max = max(object$predictions),
        sd = sd(object$predictions)
    )

    # Summary of parameters and GCV scores
    param_stats <- list(
        min_param = min(object$all_parameters),
        max_param = max(object$all_parameters),
        n_params = length(object$all_parameters),
        min_gcv = min(object$all_gcv_scores),
        max_gcv = max(object$all_gcv_scores),
        gcv_range = diff(range(object$all_gcv_scores))
    )

    # Calculate improvement percentage
    improvement_pct <- 100 * (param_stats$max_gcv - param_stats$min_gcv) / param_stats$max_gcv

    # Create summary object
    result <- list(
        method = object$method,
        optimal_parameter = object$optimal_parameter,
        gcv_score = object$gcv_score,
        compute_time_ms = object$compute_time_ms,
        predictions = pred_stats,
        parameters = param_stats,
        improvement = improvement_pct
    )

    class(result) <- "summary.nerve_cx_spectral_filter"
    return(result)
}


#' Print Method for summary.nerve_cx_spectral_filter Objects
#'
#' @description
#' Prints a formatted summary of nerve complex spectral filtering results.
#'
#' @param x A \code{summary.nerve_cx_spectral_filter} object.
#' @param digits Integer indicating the number of significant digits to print.
#'   Default is 3.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} is returned invisibly.
#'
#' @method print summary.nerve_cx_spectral_filter
#' @export
print.summary.nerve_cx_spectral_filter <- function(x, digits = 3, ...) {
    cat("Summary of Nerve Complex Spectral Filter Results\n")
    cat("================================================\n\n")

    # Method information
    cat("Method Configuration:\n")
    cat("  Laplacian type:     ", x$method$laplacian.type, "\n")
    cat("  Filter type:        ", x$method$filter.type, "\n")
    cat("  Laplacian power:    ", x$method$laplacian.power, "\n")
    cat("  Dimension weights:  ",
        paste(format(x$method$dim.weights, digits = 2), collapse = ", "), "\n")
    cat("  Eigenvectors used:  ", x$method$n.evectors, "\n\n")

    # Optimization results
    cat("Parameter Optimization:\n")
    cat("  Optimal parameter:  ", format(x$optimal_parameter, digits = digits + 1), "\n")
    cat("  GCV score:          ", format(x$gcv_score, digits = digits + 3), "\n")
    cat("  Parameter range:    [", format(x$parameters$min_param, digits = digits + 1),
        ", ", format(x$parameters$max_param, digits = digits + 1), "]\n")
    cat("  Parameters tested:  ", x$parameters$n_params, "\n")
    cat("  GCV improvement:    ", format(x$improvement, digits = digits), "%\n\n")

    # Prediction statistics
    cat("Smoothed Values Summary:\n")
    cat("   Min.   1st Qu.    Median      Mean   3rd Qu.      Max.       SD\n")
    cat(sprintf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                x$predictions$min, x$predictions$q1, x$predictions$median,
                x$predictions$mean, x$predictions$q3, x$predictions$max,
                x$predictions$sd))

    cat("\nComputation time: ", format(x$compute_time_ms, digits = digits + 1), " ms\n")

    invisible(x)
}


#' Compare Graph vs. Nerve Complex Spectral Filtering
#'
#' @description
#' Compares the performance of graph-based spectral filtering (using only the
#' 1-skeleton) with nerve complex spectral filtering (incorporating higher-order
#' simplicial information). This function quantifies the improvement gained by
#' utilizing the full simplicial structure.
#'
#' @param complex A nerve complex object created by \code{\link{create.nerve.complex}}.
#' @param y Numeric vector of observed function values at vertices.
#' @param y.true Optional numeric vector of true function values for MSE calculation.
#'   If provided, must have the same length as \code{y}.
#' @param laplacian.type Character string specifying the Laplacian type.
#'   See \code{\link{nerve.cx.spectral.filter}} for options.
#' @param filter.type Character string specifying the filter type.
#'   See \code{\link{nerve.cx.spectral.filter}} for options.
#' @param laplacian.power Positive integer for the Laplacian power.
#' @param dim.weights.complex Numeric vector of dimension weights for complex filtering.
#'   If \code{NULL}, uses exponentially decreasing weights.
#' @param graph.only.weights Numeric vector of dimension weights for graph filtering.
#'   If \code{NULL}, uses weight 1 for dimension 0 and 0 for others.
#' @param kernel.params List of kernel parameters. See \code{\link{nerve.cx.spectral.filter}}.
#' @param n.evectors Number of eigenvectors to compute.
#' @param n.candidates Number of candidate parameter values to test.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list of class \code{nerve_cx_comparison} containing:
#'   \item{graph_predictions}{Smoothed values using graph-only filtering}
#'   \item{complex_predictions}{Smoothed values using full complex filtering}
#'   \item{graph_gcv}{GCV score for graph-only filtering}
#'   \item{complex_gcv}{GCV score for complex filtering}
#'   \item{gcv_improvement_pct}{Percentage improvement in GCV score}
#'   \item{mse_graph}{MSE for graph filtering (if \code{y.true} provided)}
#'   \item{mse_complex}{MSE for complex filtering (if \code{y.true} provided)}
#'   \item{mse_improvement_pct}{Percentage improvement in MSE (if \code{y.true} provided)}
#'   \item{graph_result}{Full result object from graph filtering}
#'   \item{complex_result}{Full result object from complex filtering}
#'
#' @details
#' This comparison function runs two separate filtering operations:
#' \enumerate{
#'   \item Graph-only filtering: Uses only the 1-skeleton (edges) of the complex
#'   \item Full complex filtering: Incorporates all simplex dimensions
#' }
#'
#' The improvement metrics help quantify the benefit of using higher-order
#' simplicial information. Positive improvement percentages indicate that
#' the complex filtering performs better than graph-only filtering.
#'
#' @examples
#' \dontrun{
#' # Generate test data
#' set.seed(123)
#' coords <- matrix(runif(200), ncol = 2)
#' f_true <- function(x) sin(2*pi*x[1]) * cos(2*pi*x[2])
#' y_true <- apply(coords, 1, f_true)
#' y_noisy <- y_true + rnorm(length(y_true), 0, 0.2)
#'
#' # Create nerve complex
#' complex <- create.nerve.complex(coords, k = 8, max.dim = 2)
#' complex <- set.complex.function.values(complex, y_noisy)
#'
#' # Compare filtering approaches
#' comparison <- compare.graph.vs.nerve.cx.filtering(
#'   complex, y_noisy, y_true,
#'   dim.weights.complex = c(1.0, 0.5, 0.25),
#'   verbose = TRUE
#' )
#'
#' # Print improvement
#' cat("GCV improvement:", comparison$gcv_improvement_pct, "%\n")
#' cat("MSE improvement:", comparison$mse_improvement_pct, "%\n")
#' }
#'
#' @export compare.graph.vs.nerve.cx.filtering
compare.graph.vs.nerve.cx.filtering <- function(
    complex,
    y,
    y.true = NULL,
    laplacian.type = "STANDARD",
    filter.type = "HEAT",
    laplacian.power = 1,
    dim.weights.complex = NULL,
    graph.only.weights = NULL,
    kernel.params = list(tau.factor = 0.01, radius.factor = 3.0),
    n.evectors = 100,
    n.candidates = 100,
    verbose = FALSE) {

    # Input validation
    if (!inherits(complex, "nerve_complex")) {
        stop("'complex' must be a nerve_complex object", call. = FALSE)
    }

    if (!is.null(y.true)) {
        if (!is.numeric(y.true) || length(y.true) != length(y)) {
            stop("'y.true' must be a numeric vector with the same length as 'y'", call. = FALSE)
        }
    }

    # Set default dimension weights if not provided
    if (is.null(dim.weights.complex)) {
        # Exponentially decreasing weights for higher dimensions
        dim.weights.complex <- 2^(-(0:complex$max_dimension))
    }

    # For graph-only, set weights for dimensions > 0 to zero
    if (is.null(graph.only.weights)) {
        graph.only.weights <- numeric(complex$max_dimension + 1)
        graph.only.weights[1] <- 1  # Only use 0-dimensional (vertex) information
    }

    # Ensure weight vectors have correct length
    min_length <- complex$max_dimension + 1
    if (length(dim.weights.complex) < min_length) {
        dim.weights.complex <- c(dim.weights.complex,
                               rep(0, min_length - length(dim.weights.complex)))
    }
    if (length(graph.only.weights) < min_length) {
        graph.only.weights <- c(graph.only.weights,
                              rep(0, min_length - length(graph.only.weights)))
    }

    # Run graph-only filtering
    if (verbose) {
        message("\nRunning graph-only spectral filtering...")
    }

    graph_result <- nerve.cx.spectral.filter(
        complex = complex,
        y = y,
        laplacian.type = laplacian.type,
        filter.type = filter.type,
        laplacian.power = laplacian.power,
        dim.weights = graph.only.weights,
        kernel.params = kernel.params,
        n.evectors = n.evectors,
        n.candidates = n.candidates,
        verbose = verbose
    )

    # Run full complex filtering
    if (verbose) {
        message("\nRunning full nerve complex spectral filtering...")
    }

    complex_result <- nerve.cx.spectral.filter(
        complex = complex,
        y = y,
        laplacian.type = laplacian.type,
        filter.type = filter.type,
        laplacian.power = laplacian.power,
        dim.weights = dim.weights.complex,
        kernel.params = kernel.params,
        n.evectors = n.evectors,
        n.candidates = n.candidates,
        verbose = verbose
    )

    # Calculate improvement in GCV score
    gcv_improvement_pct <- 100 * (graph_result$gcv_score - complex_result$gcv_score) /
                          graph_result$gcv_score

    # Initialize MSE-related metrics
    mse_graph <- NULL
    mse_complex <- NULL
    mse_improvement_pct <- NULL

    # Calculate MSE if true values are provided
    if (!is.null(y.true)) {
        mse_graph <- mean((graph_result$predictions - y.true)^2)
        mse_complex <- mean((complex_result$predictions - y.true)^2)

        if (mse_graph > 0) {
            mse_improvement_pct <- 100 * (mse_graph - mse_complex) / mse_graph
        }
    }

    # Prepare results
    result <- list(
        graph_predictions = graph_result$predictions,
        complex_predictions = complex_result$predictions,
        graph_gcv = graph_result$gcv_score,
        complex_gcv = complex_result$gcv_score,
        gcv_improvement_pct = gcv_improvement_pct,
        mse_graph = mse_graph,
        mse_complex = mse_complex,
        mse_improvement_pct = mse_improvement_pct,
        graph_result = graph_result,
        complex_result = complex_result
    )

    class(result) <- "nerve_cx_comparison"

    if (verbose) {
        message("\nComparison complete:")
        message(sprintf("  GCV improvement: %.1f%%", gcv_improvement_pct))
        if (!is.null(mse_improvement_pct)) {
            message(sprintf("  MSE improvement: %.1f%%", mse_improvement_pct))
        }
    }

    return(result)
}

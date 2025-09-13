#' Model-Averaged Locally Weighted Scatterplot Smoothing (MABILO)
#'
#' @description
#' Implements MABILO algorithm for robust local regression, extending LOWESS by incorporating
#' model averaging and Bayesian bootstrap for uncertainty quantification. The algorithm uses
#' symmetric k-hop neighborhoods and kernel-weighted averaging for predictions.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates (response values).
#' @param y.true Optional numeric vector of true y values for error calculation.
#' @param k.min Minimum number of neighbors on each side (positive integer).
#' @param k.max Maximum number of neighbors on each side. If NULL, defaults to min((n-2)/4, max(3*k.min, 10)).
#' @param n.bb Number of Bayesian bootstrap iterations (non-negative integer). Set to 0 to skip bootstrap.
#' @param p Probability level for credible intervals (between 0 and 1).
#' @param kernel.type Integer; kernel type for weight calculation:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 3: Truncated exponential
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'          \item 8: cosine
#'          \item 9: hyperbolic
#'          \item 10: constant
#'        }
#'        Default is 7.
#' @param dist.normalization.factor Positive number for distance normalization (default: 1.1).
#' @param epsilon Small positive number for numerical stability (default: 1e-10).
#' @param verbose Logical; if TRUE, prints progress information.
#'
#' @return A list of class "mabilo" containing:
#' \itemize{
#'   \item k_values - Vector of tested k values
#'   \item opt_k - Optimal k value for model averaging
#'   \item opt_k_idx - Index of optimal k value
#'   \item k_mean_errors - Mean LOOCV errors for each k
#'   \item k_mean_true_errors - Mean true errors if y.true provided
#'   \item ma_predictions - Model-averaged predictions using optimal k
#'   \item k_predictions - Model-averaged predictions for all k values
#'   \item bb_predictions - Central location of bootstrap estimates (if n.bb > 0)
#'   \item cri_L - Lower bounds of credible intervals (if n.bb > 0)
#'   \item cri_U - Upper bounds of credible intervals (if n.bb > 0)
#'   \item x_sorted - Input x values sorted in ascending order
#'   \item y_sorted - y values sorted corresponding to x_sorted
#'   \item y_true_sorted - y.true values sorted corresponding to x_sorted
#'   \item k_min - Minimum neighborhood size used
#'   \item k_max - Maximum neighborhood size used
#' }
#'
#' @details
#' The function automatically sorts input data by x values. For each point, it uses
#' k-hop neighborhoods (k points on each side when available) rather than k-nearest
#' neighbors, providing more symmetric neighborhoods. The optimal k is selected by
#' minimizing mean LOOCV errors.
#'
#' The Bayesian bootstrap analysis (when n.bb > 0) provides uncertainty
#' quantification through credible intervals computed at the specified
#' probability level p. Note that setting n.bb > 0 increases computation time
#' proportionally, as the algorithm must be run n.bb times with different weight
#' configurations.
#'
#' ## Bayesian Bootstrap
#'
#' The Bayesian bootstrap, introduced by Rubin (1981), is a variant of the classical
#' bootstrap that generates smooth posterior distributions. Instead of resampling with
#' replacement (which gives discrete weights n_i/n where n_i is the number of times
#' observation i is selected), the Bayesian bootstrap assigns continuous weights to
#' each observation.
#'
#' Specifically, for each bootstrap iteration:
#'
#' 1. Generate weights (w_1, ..., w_n) from a Dirichlet(1, 1, ..., 1) distribution
#' 2. These weights sum to 1 and provide a smooth reweighting of the data
#' 3. Compute the MABILO estimate using these weights
#'
#' This approach has several advantages:
#' - Provides smooth posterior distributions rather than discrete ones
#' - Every observation contributes to each bootstrap sample (with varying weights)
#' - Naturally incorporates uncertainty in a Bayesian framework
#' - Often produces less variable estimates than classical bootstrap
#'
#' The credible intervals computed from Bayesian bootstrap samples can be interpreted
#' as Bayesian posterior intervals under a noninformative prior. See Rubin (1981)
#' "The Bayesian Bootstrap" and Lo (1987) "A Large Sample Study of the Bayesian
#' Bootstrap" for theoretical details.
#'
#' @references
#' Rubin, D.B. (1981). The Bayesian Bootstrap. The Annals of Statistics, 9(1), 130-134.
#'
#' Lo, A.Y. (1987). A Large Sample Study of the Bayesian Bootstrap. The Annals of
#' Statistics, 15(1), 360-375.
#'
#' Newton, M.A. & Raftery, A.E. (1994). Approximate Bayesian Inference with the
#' Weighted Likelihood Bootstrap. Journal of the Royal Statistical Society, Series B,
#' 56(1), 3-48.
#'
#'
#' @examples
#' # Basic usage
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#' fit <- mabilo(x, y, k.min = 3, k.max = 10)
#' plot(x, y)
#' lines(x, fit$predictions, col = "red")
#'
#' # With Bayesian bootstrap
#' fit_bb <- mabilo(x, y, k.min = 3, k.max = 10, n.bb = 100, p = 0.95)
#' lines(x, fit_bb$cri_L, col = "blue", lty = 2)
#' lines(x, fit_bb$cri_U, col = "blue", lty = 2)
#'
#' @export
mabilo <- function(x,
                   y,
                   y.true = NULL,
                   k.min = max(3, as.integer(0.05 * length(x))),
                   k.max = NULL,
                   n.bb = 100,
                   p = 0.95,
                   kernel.type = 7L,
                   dist.normalization.factor = 1.1,
                   epsilon = 1e-10,
                   verbose = FALSE) {

    # Input validation for x and y
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("x and y must be numeric vectors")
    }

    n <- length(x)

    if (length(y) != n) {
        stop("x and y must have the same length")
    }

    if (n < 10) {
        stop("Need at least 10 observations")
    }

    # Check y.true if provided
    if (!is.null(y.true)) {
        if (!is.numeric(y.true)) {
            stop("y.true must be numeric")
        }
        if (length(y.true) != n) {
            stop("y.true must have the same length as x and y")
        }
    } else {
        y.true <- numeric(0)  # Empty vector for C call
    }

    # Validate k.min
    if (!is.numeric(k.min) || length(k.min) != 1) {
        stop("k.min must be a single numeric value")
    }
    k.min <- as.integer(round(k.min))
    if (k.min < 1) {
        stop("k.min must be a positive integer")
    }

    # Set and validate k.max
    if (is.null(k.max)) {
        k.max <- min(as.integer((n-2)/4), max(3*k.min, 10))
    } else {
        if (!is.numeric(k.max) || length(k.max) != 1) {
            stop("k.max must be a single numeric value")
        }
        k.max <- as.integer(round(k.max))
    }

    if (k.max <= k.min) {
        stop("k.max must be greater than k.min")
    }

    if (k.max >= (n-1)/2) {
        stop(sprintf("k.max (%d) is too large for data size (%d). k.max must be less than (n-1)/2 to ensure multiple windows.",
                     k.max, n))
    }

    # Validate n.bb
    if (!is.numeric(n.bb) || length(n.bb) != 1) {
        stop("n.bb must be a single numeric value")
    }
    n.bb <- as.integer(round(n.bb))
    if (n.bb < 0) {
        stop("n.bb must be a non-negative integer")
    }

    # Validate p
    if (!is.numeric(p) || length(p) != 1) {
        stop("p must be a single numeric value")
    }
    if (p <= 0 || p >= 1) {
        stop("p must be between 0 and 1 (exclusive)")
    }

    # Validate kernel.type
    if (!is.numeric(kernel.type) || length(kernel.type) != 1) {
        stop("kernel.type must be a single numeric value")
    }
    kernel.type <- as.integer(round(kernel.type))
    if (kernel.type < 1 || kernel.type > 10) {
        stop("kernel.type must be an integer between 1 and 10")
    }

    # Validate dist.normalization.factor
    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1) {
        stop("dist.normalization.factor must be a single numeric value")
    }
    if (dist.normalization.factor < 1) {
        stop("dist.normalization.factor must greater than or equal to 1")
    }

    # Validate epsilon
    if (!is.numeric(epsilon) || length(epsilon) != 1) {
        stop("epsilon must be a single numeric value")
    }
    if (epsilon <= 0) {
        stop("epsilon must be positive")
    }

    # Validate verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a single logical value (TRUE/FALSE)")
    }

    # Sort data by x if necessary
    if (!identical(x, sort(x))) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        if (length(y.true) > 0) {
            y.true <- y.true[ord]
        }
    }

    # Call the C function
    result <- .Call(S_mabilo,
                    as.double(x),
                    as.double(y),
                    as.double(y.true),
                    as.integer(k.min),
                    as.integer(k.max),
                    as.integer(n.bb),
                    as.double(p),
                    as.integer(kernel.type),
                    as.double(dist.normalization.factor),
                    as.double(epsilon),
                    as.logical(verbose))

    # Add sorted data and parameters to result
    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- if (length(y.true) > 0) y.true else NULL
    result$k_min <- k.min
    result$k_max <- k.max

    class(result) <- "mabilo"

    return(result)
}


#' Plot Method for Mabilo Objects
#'
#' @description
#' Generates diagnostic and visualization plots for mabilo objects. The function
#' supports different types of plots including fitted values with optional credible
#' intervals, error diagnostics, and residual analyses.
#'
#' @param x An object of class "mabilo"
#' @param type Character string specifying the type of plot to create (partial matching supported):
#'   * "fit": Plot fitted values against x values (e.g., "f" or "fi")
#'   * "diagnostic": Plot error diagnostics for different k values (e.g., "d" or "dia")
#'   * "residuals": Plot residuals against x values (e.g., "res")
#'   * "residuals.hist": Plot histogram of residuals (e.g., "rh" or "residuals.h")
#' @param title Plot title (default: "")
#' @param xlab Label for x-axis (default: "")
#' @param ylab Label for y-axis (default: "")
#' @param with.y.true Logical; if TRUE and true values are available, adds them to the plot
#'   (default: TRUE)
#' @param with.pts Logical; if TRUE, adds original data points to the plot
#'   (default: FALSE)
#' @param with.CrI Logical; if TRUE and bootstrap results available, adds credible intervals
#'   (default: TRUE)
#' @param true.col Color for true values line (default: "red")
#' @param ma.col Color for predictions (default: "blue")
#' @param pts.col Color for data points (default: "gray60")
#' @param with.predictions.pts Logical; if TRUE, adds points at prediction values
#'   (default: FALSE)
#' @param predictions.pts.col Color for prediction points (default: "blue")
#' @param predictions.pts.pch Point character for prediction points (default: 20)
#' @param CrI.as.polygon Logical; if TRUE, draws credible intervals as a polygon
#'   (default: TRUE)
#' @param CrI.polygon.col Color for credible interval polygon (default: "gray85")
#' @param CrI.line.col Color for credible interval lines (default: "gray10")
#' @param CrI.line.lty Line type for credible interval lines (default: 2)
#' @param ylim Numeric vector of length 2 giving y-axis limits (default: NULL)
#' @param legend.cex Character expansion factor for legend (default: 0.8)
#' @param k Specific k value to plot predictions for (default: NULL, uses optimal k)
#' @param with.legend Logical; if TRUE, displays legend on the plot (default: TRUE)
#' @param legend.pos Character string specifying legend position (default: "topright")
#' @param ... Additional arguments passed to plotting functions
#'
#' @details
#' The function produces different types of plots:
#'
#' For type = "fit":
#' * Plots the fitted values against x values
#' * Optionally includes true values and data points
#' * Shows credible intervals if bootstrap results available
#' * Supports customization of colors and styles
#'
#' For type = "diagnostic":
#' * Shows LOOCV error diagnostics for different k values
#' * Indicates optimal k value with vertical line
#' * Shows true errors if available
#'
#' For type = "residuals" and "residuals.hist":
#' * Visualizes model residuals in various ways
#' * Includes reference lines and density curves
#'
#' @return
#' Invisibly returns NULL. The function is called for its side effect of producing plots.
#'
#' @examples
#' # Generate example data
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#'
#' # Fit mabilo model
#' fit <- mabilo(x, y, n.bb = 100)
#'
#' # Basic fit plot with credible intervals
#' plot(fit)
#'
#' # Diagnostic plot showing errors
#' plot(fit, type = "diagnostic")
#'
#' # Residual plot
#' plot(fit, type = "residuals")
#'
#' @seealso
#' \code{\link{mabilo}} for fitting the mabilo model
#'
#' @importFrom graphics plot points lines polygon matplot matlines hist par mtext legend abline
#' @importFrom stats density shapiro.test qnorm quantile sd
#'
#' @export
plot.mabilo <- function(x, type = "fit",
                        title = "", xlab = "", ylab = "",
                        with.y.true = TRUE,
                        with.pts = FALSE,
                        with.CrI = TRUE,
                        true.col = "red",
                        ma.col = "blue",
                        pts.col = "gray60",
                        with.predictions.pts = FALSE,
                        predictions.pts.col = "blue",
                        predictions.pts.pch = 20,
                        CrI.as.polygon = TRUE,
                        CrI.polygon.col = "gray85",
                        CrI.line.col = "gray10",
                        CrI.line.lty = 2,
                        ylim = NULL,
                        legend.cex = 0.8,
                        k = NULL,
                        with.legend = TRUE,
                        legend.pos = "topright",
                        ...) {
    if (!inherits(x, "mabilo")) {
        stop("Input must be a 'mabilo' object")
    }

    # Define valid plot types
    valid_types <- c("fit", "diagnostic", "residuals", "residuals.hist")

    # First check for exact match
    if (type %in% valid_types) {
        # Use exact match
        type_match <- type
    } else {
        # If no exact match, try partial matching
        type_match <- grep(paste0("^", type), valid_types, value = TRUE)

        # Handle matching results
        if (length(type_match) == 0) {
            stop("Invalid plot type. Use 'fit', 'diagnostic', 'residuals' or 'residuals.hist'")
        } else if (length(type_match) > 1) {
            stop(sprintf("Ambiguous type. '%s' matches multiple types: %s",
                        type, paste(type_match, collapse = ", ")))
        }

        # Use the matched type
        type <- type_match
    }

    # Validate k if provided
    if (!is.null(k)) {
        if (!is.numeric(k) || length(k) != 1) {
            stop("k must be a single numeric value")
        }
        if (k < x$k_min || k > x$k_max) {
            stop(sprintf("k must be between %d and %d", x$k_min, x$k_max))
        }
        # Convert k to index
        k_idx <- which(x$k_values == k)
        if (length(k_idx) == 0) {
            stop(sprintf("k=%d not found in k_values", k))
        }
    }

    res <- x
    switch(type,
           "fit" = {
               mabilo.plot.fit(res, title, xlab, ylab, with.y.true,
                               with.pts, with.CrI, true.col, ma.col, pts.col,
                               with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                               CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                               CrI.line.lty, ylim, legend.cex, k = k,
                               with.legend = with.legend,
                               legend.pos = legend.pos,
                               ...)
           },
           "diagnostic" = {
               mabilo.plot.diagnostic(res, title, xlab, ylab,
                                      ma.col, true.col, legend.cex,
                                      legend.pos, ...)
           },
           "residuals" = {
               mabilo.plot.residuals(res, title, xlab, ylab, ...)
           },
           "residuals.hist" = {
               mabilo.plot.residuals.hist(res, title, xlab, ylab, ...)
           }
    )
    invisible(NULL)
}

# Helper function for fit plots
mabilo.plot.fit <- function(res, title, xlab, ylab, with.y.true,
                            with.pts, with.CrI, true.col, ma.col, pts.col,
                            with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                            CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                            CrI.line.lty, ylim, legend.cex, k = NULL,
                            with.legend = TRUE,
                            legend.pos = "topright", ...) {

    # Select appropriate predictions based on k
    if (!is.null(k)) {
        k_idx <- which(res$k_values == k)
        predictions <- res$k_predictions[[k_idx]]
        pred_label <- sprintf("Predictions (k=%d)", k)
    } else {
        predictions <- res$predictions
        pred_label <- "Predictions"
    }

    # Calculate y-limits if not provided
    if (is.null(ylim)) {
        ylim.data <- c(res$y_sorted, predictions)
        if (with.CrI && !is.null(res$cri_L)) {
            ylim.data <- c(ylim.data, res$cri_L, res$cri_U)
        }
        ylim <- range(ylim.data, na.rm = TRUE)
    }

    # Initialize plot
    plot(res$x_sorted, res$y_sorted, type = "n",
         las = 1, ylim = ylim, xlab = xlab, ylab = ylab,
         main = title, ...)

    # Add credible intervals if available and requested
    has_cri <- with.CrI && !is.null(res$cri_L) && !is.null(res$cri_U)
    if (has_cri) {
        if (CrI.as.polygon) {
            polygon(c(res$x_sorted, rev(res$x_sorted)),
                   c(res$cri_L, rev(res$cri_U)),
                   col = CrI.polygon.col, border = NA)
        } else {
            lines(res$x_sorted, res$cri_L, col = CrI.line.col, lty = CrI.line.lty)
            lines(res$x_sorted, res$cri_U, col = CrI.line.col, lty = CrI.line.lty)
        }
    }

    # Add predictions
    lines(res$x_sorted, predictions, col = ma.col, ...)

    # Add original data points if requested
    if (with.pts) {
        points(res$x_sorted, res$y_sorted, col = pts.col, ...)
    }

    # Add prediction points if requested
    if (with.predictions.pts) {
        points(res$x_sorted, predictions,
              col = predictions.pts.col,
              pch = predictions.pts.pch)
    }

    # Add true values if available and requested
    if (with.y.true && !is.null(res$y_true_sorted) && length(res$y_true_sorted) == length(res$x_sorted)) {
        lines(res$x_sorted, res$y_true_sorted, col = true.col)
    }

    # Add legend
    legend_items <- c()
    legend_cols <- c()
    legend_ltys <- c()
    legend_pchs <- c()

    if (with.pts) {
        legend_items <- c(legend_items, "Data Points")
        legend_cols <- c(legend_cols, pts.col)
        legend_ltys <- c(legend_ltys, NA)
        legend_pchs <- c(legend_pchs, 1)
    }

    legend_items <- c(legend_items, pred_label)
    legend_cols <- c(legend_cols, ma.col)
    legend_ltys <- c(legend_ltys, 1)
    legend_pchs <- c(legend_pchs, NA)

    if (has_cri) {
        legend_items <- c(legend_items, "95% Credible Interval")
        legend_cols <- c(legend_cols, if(CrI.as.polygon) CrI.polygon.col else CrI.line.col)
        legend_ltys <- c(legend_ltys, if(CrI.as.polygon) 1 else CrI.line.lty)
        legend_pchs <- c(legend_pchs, NA)
    }

    if (with.y.true && !is.null(res$y_true_sorted)) {
        legend_items <- c(legend_items, "True Values")
        legend_cols <- c(legend_cols, true.col)
        legend_ltys <- c(legend_ltys, 1)
        legend_pchs <- c(legend_pchs, NA)
    }

    if (with.legend) {
        legend(legend.pos, legend = legend_items,
               col = legend_cols, lty = legend_ltys, pch = legend_pchs,
               bg = "white", inset = 0.05, cex = legend.cex)
    }
}

# Helper function for diagnostic plots
mabilo.plot.diagnostic <- function(res, title, xlab, ylab, ma.col, true.col, legend.cex, legend.pos = "topright", ...) {
    # Prepare error data
    err_data <- matrix(res$k_mean_errors, ncol = 1)
    labels <- "LOOCV Errors"
    cols <- ma.col
    pchs <- 1

    # Add true errors if available
    if (!is.null(res$k_mean_true_errors)) {
        err_data <- cbind(err_data, res$k_mean_true_errors)
        labels <- c(labels, "True Errors")
        cols <- c(cols, true.col)
        pchs <- c(pchs, 2)
    }

    # Plot errors
    matplot(res$k_values, err_data, type = 'b',
            xlab = if(xlab == "") "k value" else xlab,
            ylab = if(ylab == "") "Error" else ylab,
            main = if(title == "") "Error Diagnostic Plot" else title,
            col = cols,
            pch = pchs, las = 1,
            ...)

    # Add optimal k value
    abline(v = res$opt_k, col = ma.col, lty = 2)
    mtext(sprintf("Optimal k = %d", res$opt_k),
          side = 3, line = 0.25, at = res$opt_k, col = ma.col)

    # Add legend
    legend(legend.pos, legend = labels,
           col = cols,
           pch = pchs,
           lty = 1, bg = "white", inset = 0.05, cex = legend.cex)
}

# Helper function for residual plots
mabilo.plot.residuals <- function(res, title, xlab, ylab, ...) {
    residuals <- res$y_sorted - res$predictions

    plot(res$x_sorted, residuals,
         xlab = if(xlab == "") "x" else xlab,
         ylab = if(ylab == "") "Residuals" else ylab,
         main = if(title == "") "Residuals Plot" else title,
         las = 1, ...)

    abline(h = 0, lty = 2)
}

# Helper function for residual histograms
mabilo.plot.residuals.hist <- function(res, title, xlab, ylab, ...) {
    residuals <- res$y_sorted - res$predictions

    hist(residuals,
         xlab = if(xlab == "") "Residuals" else xlab,
         ylab = if(ylab == "") "Frequency" else ylab,
         main = if(title == "") "Histogram of Residuals" else title,
         probability = TRUE,
         ...)

    lines(density(residuals), col = "blue")
    abline(v = 0, lty = 2)
}

#' Compute Summary Statistics for Mabilo Objects
#'
#' @description
#' Computes comprehensive summary statistics for mabilo fits including model parameters,
#' fit statistics, error analysis, and diagnostic information for model-averaged predictions.
#' Also includes Bayesian bootstrap statistics when available.
#'
#' @param object A 'mabilo' object
#' @param quantiles Numeric vector of probabilities for quantile computations
#' @param ... Additional arguments (currently unused)
#'
#' @return A 'summary.mabilo' object containing:
#' \itemize{
#'   \item model_info: Basic information about the model fit
#'   \item fit_stats: Prediction accuracy metrics
#'   \item k_error_stats: Error statistics for different k values
#'   \item residual_stats: Residual analysis results
#'   \item true_error_stats: Statistics comparing to true values (if available)
#'   \item bootstrap_stats: Bayesian bootstrap statistics (if available)
#' }
#'
#' @importFrom stats shapiro.test qnorm quantile sd
#'
#' @export
summary.mabilo <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    if (!inherits(object, "mabilo")) {
        stop("Input must be a 'mabilo' object")
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("quantiles must be numeric values between 0 and 1")
    }

    # Check for missing values
    if (any(is.na(object$y_sorted)) || any(is.na(object$predictions))) {
        warning("Missing values detected in fit results")
    }

    # Calculate residuals
    residuals <- object$y_sorted - object$predictions

    # Basic model information
    model_info <- list(
        n_observations = length(object$x_sorted),
        optimal_k = object$opt_k,
        k_range = range(object$k_values),
        min_x = min(object$x_sorted),
        max_x = max(object$x_sorted),
        range_x = diff(range(object$x_sorted))
    )

    # Fit statistics
    fit_stats <- list(
        mse = mean(residuals^2),
        rmse = sqrt(mean(residuals^2)),
        mae = mean(abs(residuals)),
        median_ae = median(abs(residuals))
    )

    # True error statistics if true values are available
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted)) {
        true_residuals <- object$y_true_sorted - object$predictions
        true_mse <- mean(true_residuals^2)
        equiv_normal_sd <- sqrt(true_mse)
        normalized_errors <- true_residuals / equiv_normal_sd
        normality_test <- shapiro.test(normalized_errors)
        error_quantiles <- quantile(normalized_errors,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
        theoretical_quantiles <- qnorm(c(0.025, 0.25, 0.5, 0.75, 0.975))

        true_error_stats <- list(
            true_mse = true_mse,
            true_rmse = sqrt(true_mse),
            true_mae = mean(abs(true_residuals)),
            true_median_ae = median(abs(true_residuals)),
            equiv_normal_sd = equiv_normal_sd,
            normalized_error_stats = list(
                mean = mean(normalized_errors),
                sd = sd(normalized_errors),
                quantiles = error_quantiles
            ),
            theoretical_normal_quantiles = theoretical_quantiles,
            shapiro_test = list(
                statistic = normality_test$statistic,
                p_value = normality_test$p_value
            ),
            relative_efficiency = mean(abs(true_residuals)) /
                                (sqrt(2/pi) * equiv_normal_sd)
        )
    }

    # Error statistics for different k values
    k_error_stats <- list(
        mean_errors = object$k_mean_errors,
        min_error = min(object$k_mean_errors),
        optimal_k_error = object$k_mean_errors[object$opt_k_idx]
    )

    # Residual statistics
    residual_stats <- list(
        mean = mean(residuals),
        sd = sd(residuals),
        quantiles = quantile(residuals, probs = quantiles)
    )

    # Bootstrap statistics if available
    bootstrap_stats <- NULL
    if (!is.null(object$bb_predictions)) {
        bootstrap_stats <- list(
            mean_prediction = object$bb_predictions,
            lower_ci = object$cri_L,
            upper_ci = object$cri_U,
            ci_width = object$cri_U - object$cri_L,
            mean_ci_width = mean(object$cri_U - object$cri_L)
        )
    }

    # Create summary object
    result <- list(
        call = object$call,
        model_info = model_info,
        fit_stats = fit_stats,
        k_error_stats = k_error_stats,
        residual_stats = residual_stats,
        true_error_stats = true_error_stats,
        bootstrap_stats = bootstrap_stats
    )

    class(result) <- "summary.mabilo"
    return(result)
}

#' Print Summary Statistics for Mabilo Fits
#'
#' @description
#' Formats and displays comprehensive summary statistics for mabilo fits.
#' Includes model parameters, fit statistics, error analysis, diagnostic information,
#' and Bayesian bootstrap results when available.
#'
#' @param x A 'summary.mabilo' object from summary.mabilo
#' @param digits Number of significant digits for numerical output (default: 4)
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.summary.mabilo <- function(x, digits = 4, ...) {
    ## Helper function for separator lines
    hr <- function() cat(paste0(rep("\u2500", 80), collapse = ""), "\n")

    # Helper function for section headers
    section_header <- function(text) {
        cat("\n")
        hr()
        cat(text, "\n")
        hr()
    }

    # Main title
    cat("\n")
    section_header("MABILO (Model-Averaged Local Weighted Smoothing) SUMMARY")

    # Model Information
    cat("\nModel Information:\n")
    cat(sprintf("Number of observations:         %d\n", x$model_info$n_observations))
    cat(sprintf("Optimal k:                      %d\n", x$model_info$optimal_k))
    cat(sprintf("k range:                        [%d, %d]\n",
                x$model_info$k_range[1], x$model_info$k_range[2]))
    cat(sprintf("X range:                        [%.3f, %.3f]\n",
                x$model_info$min_x, x$model_info$max_x))
    cat(sprintf("X span:                         %.3f\n", x$model_info$range_x))

    # Fit Statistics
    section_header("FIT STATISTICS")
    cat(sprintf("MSE:                            %.4f\n", x$fit_stats$mse))
    cat(sprintf("RMSE:                           %.4f\n", x$fit_stats$rmse))
    cat(sprintf("MAE:                            %.4f\n", x$fit_stats$mae))
    cat(sprintf("Median AE:                      %.4f\n", x$fit_stats$median_ae))

    # True Error Statistics (if available)
    if (!is.null(x$true_error_stats)) {
        section_header("TRUE ERROR STATISTICS")
        cat(sprintf("True MSE:                       %.4f\n", x$true_error_stats$true_mse))
        cat(sprintf("True RMSE:                      %.4f\n", x$true_error_stats$true_rmse))
        cat(sprintf("True MAE:                       %.4f\n", x$true_error_stats$true_mae))
        cat(sprintf("Equivalent Normal SD:           %.4f\n", x$true_error_stats$equiv_normal_sd))
        cat(sprintf("Relative efficiency:            %.4f\n", x$true_error_stats$relative_efficiency))
        cat(sprintf("Shapiro-Wilk p-value:          %.4f\n",
                   x$true_error_stats$shapiro_test$p_value))
    }

    # k Error Statistics
    section_header("K ERROR STATISTICS")
    cat(sprintf("Minimum Error:                   %.4f\n", x$k_error_stats$min_error))
    cat(sprintf("Optimal k Error:                 %.4f\n", x$k_error_stats$optimal_k_error))

    # Residual Statistics
    section_header("RESIDUAL STATISTICS")
    cat(sprintf("Mean:                            %.4f\n", x$residual_stats$mean))
    cat(sprintf("Standard Deviation:              %.4f\n", x$residual_stats$sd))
    cat("\nQuantiles:\n")
    print(format(data.frame(
        Quantile = names(x$residual_stats$quantiles),
        Value = x$residual_stats$quantiles
    ), digits = 4), row_names = FALSE)

    # Bootstrap Statistics (if available)
    if (!is.null(x$bootstrap_stats)) {
        section_header("BAYESIAN BOOTSTRAP STATISTICS")
        cat(sprintf("Mean CI Width:                  %.4f\n", x$bootstrap_stats$mean_ci_width))
        cat("\nCredible Interval Summary:\n")
        ci_summary <- summary(x$bootstrap_stats$ci_width)
        print(ci_summary, digits = 4)
    }

    cat("\n")  # Final newline for spacing
    invisible(x)
}

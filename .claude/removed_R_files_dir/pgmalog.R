

#' Adaptive Neighborhood Size Graph Path Linear Model for Univariate Data
#'
#' @description
#' Implements adaptive neighborhood selection for univariate graph path linear models (PGMALOG)
#' using path doubly weighted linear models (PWWLM). The function automatically determines
#' optimal neighborhood sizes both globally and locally for each point.
#'
#' @param x Numeric vector of predictor values, will be automatically sorted if not in order
#' @param y Numeric vector of response values corresponding to x
#' @param y.true Optional numeric vector of true response values for error calculation
#' @param use.median Logical; if TRUE uses median instead of mean for central tendency
#' @param h.min Minimum neighborhood size (must be even and >= 4)
#' @param h.max Maximum neighborhood size (must be even and <= n - 2, where n is length of x)
#' @param p Probability level for credible intervals (between 0 and 1)
#' @param n.bb Number of bootstrap iterations (>= 0, 0 skips bootstrap)'
#' @param kernel.type Integer; kernel type for weight calculation:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#'        Default is 7.
#' @param n.cores Number of cores for parallel processing
#' @param dist.normalization.factor Distance normalization factor (between 1 and (h-1)/2)
#' @param epsilon Numerical stability parameter (> 0)
#' @param seed Random seed for reproducibility
#' @param verbose Logical; if TRUE prints progress messages
#'
#' @details
#' The function constructs a chain graph from sorted x values and applies path linear
#' model estimation with adaptive neighborhood selection. It performs:
#' 1. Chain graph construction from ordered x values
#' 2. Path linear model estimation for different neighborhood sizes
#' 3. Cross-validation for optimal neighborhood selection
#' 4. Optional bootstrap for uncertainty quantification
#'
#' The function uses path distances in the chain graph and applies double weighting:
#' first by path distance, then by position within paths.
#'
#' @return An object of class "pgmalog" containing:
#' \itemize{
#'   \item h_values - Vector of tested neighborhood sizes
#'   \item cv_errors - Cross-validation errors for each h
#'   \item opt_h - Optimal neighborhood size
#'   \item opt_graph_adj_list - Adjacency list for optimal graph
#'   \item opt_graph_edge_lengths - Edge lengths for optimal graph
#'   \item opt_predictions - Conditional expectations using optimal h
#'   \item opt_bb_predictions - Bootstrap estimates (if n.bb > 0)
#'   \item opt_local_predictions - locally optimized estimates
#'   \item opt_ci_lower - Lower credible bounds (if n.bb > 0)
#'   \item opt_ci_upper - Upper credible bounds (if n.bb > 0)
#'   \item true_error - Mean absolute error (if y.true provided)
#'   \item x_sorted - Sorted x values
#'   \item y_sorted - Y values corresponding to sorted x
#'   \item y_true_sorted - True y values corresponding to sorted x (if provided)
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#'
#' # Fit model with default parameters
#' fit <- univariate.pgmalog(x, y)
#'
#' # Fit model with custom parameters and bootstrapping
#' fit2 <- univariate.pgmalog(
#'   x, y,
#'   h.min = 5,
#'   h.max = 21,
#'   n.bb = 1000,
#'   use.median = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{plot.pgmalog}} for plotting methods
#'
#' @references
#' Add relevant references here
#'
#' @export
upgmalog <- function(x, y, y.true = NULL,
                    use.median = TRUE,
                    h.min = 4,
                    h.max = min(30, length(x) - 2),
                    p = 0.95,
                    n.bb = 50,
                    bb.max.distance.deviation = 1,
                    n.CVs = 100,
                    n.CV.folds = 10,
                    seed = 0,
                    kernel.type = 7L,
                    n.cores = 1,
                    dist.normalization.factor = 1.01,
                    epsilon = 1e-15,
                    verbose = TRUE) {
    ## Check x and y inputs
    if (!is.numeric(x)) stop("x must be numeric")
    if (!is.numeric(y)) stop("y must be numeric")
    if (length(x) != length(y)) stop("x and y must have the same length")
    if (length(x) < 2) stop("x must have at least 2 points")
    if (any(is.na(x)) || any(is.na(y))) stop("x and y cannot contain NA values")
    if (any(is.infinite(x)) || any(is.infinite(y))) stop("x and y cannot contain infinite values")

    # Check if x is sorted and sort if necessary
    if (!identical(x, sort(x))) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        if (!is.null(y.true)) y.true <- y.true[ord]
    }

    # Check y.true if provided
    if (!is.null(y.true)) {
        if (!is.numeric(y.true)) stop("y.true must be numeric")
        if (length(y.true) != length(x)) stop("y.true must have the same length as x")
        if (any(is.na(y.true))) stop("y.true cannot contain NA values")
        if (any(is.infinite(y.true))) stop("y.true cannot contain infinite values")
    }

    # Check use.median
    if (!is.logical(use.median)) stop("use.median must be logical")
    if (length(use.median) != 1) stop("use.median must be a single value")

    # Check h.min and h.max
    if (!is.numeric(h.min) || h.min != as.integer(h.min)) stop("h.min must be an integer")
    if (!is.numeric(h.max) || h.max != as.integer(h.max)) stop("h.max must be an integer")
    if (h.min < 2) stop("h.min must be at least 2")
    if (h.max <= h.min) stop("h.max must be greater than h.min")
    if (h.max > length(x)) stop("h.max cannot exceed the number of points")
    if (h.min %% 2 == 1) stop("h.min must be even")
    if (h.max %% 2 == 1) stop("h.max must be even")
    if (h.min < 3) stop("h.min must be at least 3")
    if (h.max > length(x) - 2) stop("h.max cannot exceed (n_vertices - 2) points")

    # Check p
    if (!is.numeric(p)) stop("p must be numeric")
    if (length(p) != 1) stop("p must be a single value")
    if (p <= 0 || p >= 1) stop("p must be between 0 and 1")

    # Check n.bb
    if (!is.numeric(n.bb) || n.bb != as.integer(n.bb)) stop("n.bb must be an integer")
    if (n.bb < 0) stop("n.bb must be non-negative")

    # Check bb.max.distance.deviation
    if (!is.numeric(bb.max.distance.deviation) || bb.max.distance.deviation != as.integer(bb.max.distance.deviation)) stop("bb.max.distance.deviation must be an integer")
    if (bb.max.distance.deviation < -1) stop("bb.max.distance.deviation must be greater or equal to -1.")

    ## Validate n.CVs
    if (!is.numeric(n.CVs) || length(n.CVs) != 1 || n.CVs < 0 || n.CVs %% 1 != 0) {
        stop("n.CVs must be a non-negative integer.")
    }

    ## Validate n.CV.folds
    if (!is.numeric(n.CV.folds) || length(n.CV.folds) != 1 ||
        n.CV.folds <= 0 || n.CV.folds %% 1 != 0) {
        stop("n.CV.folds must be a positive integer.")
    }
    if (n.CV.folds > length(y)) {
        stop("n.CV.folds cannot be larger than the number of vertices.")
    }

    ## Validate seed
    if (!is.numeric(seed) || length(seed) != 1 || seed %% 1 != 0) {
        stop("seed must be an integer.")
    }

    # Check kernel.type
    if (!is.numeric(kernel.type) || kernel.type != as.integer(kernel.type)) stop("kernel.type must be an integer")
    if (!(kernel.type %in% 1:10)) stop("kernel.type must be between 1 and 10")

    # Check n.cores
    if (!is.numeric(n.cores) || n.cores != as.integer(n.cores)) stop("n.cores must be an integer")
    if (n.cores <= 0) stop("n.cores must be positive")

    # Check dist.normalization.factor
    if (!is.numeric(dist.normalization.factor)) stop("dist.normalization.factor must be numeric")
    if (length(dist.normalization.factor) != 1) stop("dist.normalization.factor must be a single value")
    if (dist.normalization.factor <= 1) stop("dist.normalization.factor must be greater than 1")

    # Check epsilon
    if (!is.numeric(epsilon)) stop("epsilon must be numeric")
    if (length(epsilon) != 1) stop("epsilon must be a single value")
    if (epsilon <= 0) stop("epsilon must be positive")

    ## Call the C++ function
    result <- .Call("S_upgmalog",
                    as.double(x),
                    as.double(y),
                    if (is.null(y.true)) double() else as.double(y.true),
                    as.logical(use.median),
                    as.integer(h.min),
                    as.integer(h.max),
                    as.double(p),
                    as.integer(n.bb),
                    as.integer(bb.max.distance.deviation),
                    as.integer(n.CVs),
                    as.integer(n.CV.folds),
                    as.integer(seed),
                    as.integer(kernel.type),
                    as.integer(n.cores),
                    as.double(dist.normalization.factor),
                    as.double(epsilon),
                    as.logical(verbose))

    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- y.true
    result$h_min <- h.min
    result$h_max <- h.max
    result$max_h_index <- (h.max - h.min) / 2

    class(result) <- "upgmalog"

    return(result)
}


#' Plot Method for upgmalog Objects
#'
#' @description
#' Creates visualizations for graph path linear model (UPGMALOG) results, including fitted values,
#' diagnostics, and residual analysis. Supports multiple types of conditional expectation
#' estimates.
#'
#' @param x An object of class "upgmalog"
#' @param type Character string specifying plot type:
#'        - "fit": fitted values and intervals
#'        - "diagnostic": cross-validation error analysis
#'        - "residuals": residual scatter plot
#'        - "residuals_hist": residual distribution
#' @param title Plot title (default: "")
#' @param xlab X-axis label (default: "")
#' @param ylab Y-axis label (default: "")
#' @param predictions.type Type of conditional expectation to plot:
#'        - "predictions": global optimal h estimate (default)
#'        - "bb.predictions": bootstrap-based estimate
#'        - "local.predictions": locally adaptive estimate
#' @param with.y.true Include true values if available (default: TRUE)
#' @param with.pts Show original data points (default: FALSE)
#' @param with.CrI Show credible intervals (default: TRUE)
#' @param y.true.col Color for true values (default: "red")
#' @param predictions.col Color for conditional expectation (default: "blue")
#' @param with.predictions.pts Show conditional expectation points (default: FALSE)
#' @param predictions.pts.col Point color for conditional expectation (default: "blue")
#' @param predictions.pts.pch Point character for conditional expectation (default: 20)
#' @param CrI.as.polygon Show credible intervals as polygon (default: TRUE)
#' @param CrI.polygon.col Credible interval polygon color (default: "gray95")
#' @param CrI.line.col Credible interval line color (default: "gray")
#' @param CrI.line.lty Credible interval line type (default: 2)
#' @param ylim Y-axis limits (default: NULL for automatic)
#' @param ... Additional graphical parameters
#'
#' @return NULL (invisibly)
#' @export
plot.upgmalog <- function(x, type = "fit",
                     title = "", xlab = "", ylab = "",
                     predictions.type = "predictions",
                     with.y.true = TRUE,
                     with.pts = FALSE,
                     with.CrI = TRUE,
                     y.true.col = "red",
                     predictions.col = "blue",
                     with.predictions.pts = FALSE,
                     predictions.pts.col = "blue",
                     predictions.pts.pch = 20,
                     CrI.as.polygon = TRUE,
                     CrI.polygon.col = "gray85",
                     CrI.line.col = "gray10",
                     CrI.line.lty = 2,
                     ylim = NULL,
                     ...) {

    if (!inherits(x, "upgmalog")) {
        stop("Input must be a 'upgmalog' object")
    }

    if (!predictions.type %in% c("predictions", "bb.predictions", "local.predictions")) {
        stop("predictions.type must be one of: 'predictions', 'bb.predictions', 'local.predictions'")
    }

    res <- x

    switch(type,
           "fit" = {
               plot.fit(res, title, xlab, ylab, predictions.type, with.y.true,
                       with.pts, with.CrI, y.true.col, predictions.col,
                       with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                       CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                       CrI.line.lty, ylim, ...)
           },
           "diagnostic" = {
               plot.diagnostic(res)
           },
           "residuals" = {
               plot.residuals(res, title, xlab, ylab, ...)
           },
           "residuals_hist" = {
               plot.residuals.hist(res, title, xlab, ylab, ...)
           },
           stop("Invalid plot type. Use 'fit', 'diagnostic', 'residuals' or 'residuals.hist'")
    )

    invisible(NULL)
}

# Internal plotting functions
plot.fit <- function(res, title, xlab, ylab, predictions.type, with.y.true,
                    with.pts, with.CrI, y.true.col, predictions.col,
                    with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                    CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                    CrI.line.lty, ylim, ...) {

    # Select appropriate conditional expectation based on type
    predictions <- switch(predictions.type,
                    "predictions" = res$opt_predictions,
                    "bb.predictions" = {
                        if (is.null(res$opt_bb_predictions)) {
                            stop("Bootstrap estimates not available")
                        }
                        res$opt_bb_predictions
                    },
                    "local.predictions" = {
                        if (is.null(res$opt_local_predictions)) {
                            stop("Local estimates not available")
                        }
                        res$opt_local_predictions
                    })

    # Calculate y-limits if not provided
    if (is.null(ylim)) {
        ylim_data <- c(res$y_sorted, predictions)
        if (with.CrI && !is.null(res$opt_ci_lower)) {
            ylim_data <- c(ylim_data, res$opt_ci_lower, res$opt_ci_upper)
        }
        ylim <- range(ylim_data, na.rm = TRUE)
    }

    # Initialize plot
    plot(res$x_sorted, res$y_sorted, type = "n",
         las = 1, ylim = ylim, xlab = xlab, ylab = ylab,
         main = title, ...)

    # Add credible intervals if requested and available
    # Only show for global or bootstrap estimates
    if (with.CrI && !is.null(res$opt_ci_lower) &&
        predictions.type %in% c("predictions", "bb.predictions")) {
        if (CrI.as.polygon) {
            polygon(c(res$x_sorted, rev(res$x_sorted)),
                   c(res$opt_ci_lower, rev(res$opt_ci_upper)),
                   col = CrI.polygon.col, border = NA)
        } else {
            matlines(res$x_sorted, cbind(res$opt_ci_lower, res$opt_ci_upper),
                    col = CrI.line.col, lty = CrI.line.lty)
        }
    }

    # Add conditional expectation line
    lines(res$x_sorted, predictions, col = predictions.col, ...)

    # Add original data points if requested
    if (with.pts) {
        points(res$x_sorted, res$y_sorted, ...)
    }

    # Add conditional expectation points if requested
    if (with.predictions.pts) {
        points(res$x_sorted, predictions,
               col = predictions.pts.col, pch = predictions.pts.pch)
    }

    # Add true values if available and requested
    if (with.y.true && !is.null(res$y_true_sorted)) {
        lines(res$x_sorted, res$y_true_sorted, col = y.true.col)
    }

    # Add legend
    legend_items <- character(0)
    legend_cols <- character(0)
    legend_ltys <- numeric(0)
    legend_pchs <- numeric(0)

    if (with.pts) {
        legend_items <- c(legend_items, "Data Points")
        legend_cols <- c(legend_cols, "black")
        legend_ltys <- c(legend_ltys, NA)
        legend_pchs <- c(legend_pchs, 1)
    }

    legend_items <- c(legend_items,
                     switch(predictions.type,
                            "predictions" = "Global Estimate",
                            "bb.predictions" = "Bootstrap Estimate",
                            "local.predictions" = "Local Estimate"))
    legend_cols <- c(legend_cols, predictions.col)
    legend_ltys <- c(legend_ltys, 1)
    legend_pchs <- c(legend_pchs, NA)

    if (with.y.true && !is.null(res$y_true_sorted)) {
        legend_items <- c(legend_items, "True Values")
        legend_cols <- c(legend_cols, y.true.col)
        legend_ltys <- c(legend_ltys, 1)
        legend_pchs <- c(legend_pchs, NA)
    }

    if (length(legend_items) > 0) {
        legend("topright", legend = legend_items,
               col = legend_cols, lty = legend_ltys, pch = legend_pchs,
               bg = "white")
    }
}

# Internal function for diagnostic plots
plot.diagnostic <- function(res) {
    if (!is.null(res$h_values) && !is.null(res$h_cv_errors)) {
        plot(res$h_values, res$h_cv_errors, type = 'b',
             las = 1, xlab = "Neighborhood Size (h)",
             ylab = "Cross-validation Error",
             main = "Cross-validation Diagnostic Plot")
        abline(v = res$opt_h, col = "gray", lty = 2)
        mtext(sprintf("Optimal h = %d", res$opt_h),
              side = 3, line = 0.25, at = res$opt_h)
    } else {
        stop("Diagnostic plot requires h_values and cv_errors in the result object")
    }
}

# Internal function for residual plots
plot.residuals <- function(res, title, xlab, ylab, predictions.type = "predictions", ...) {
    predictions <- switch(predictions.type,
                    "predictions" = res$opt_predictions,
                    "bb.predictions" = res$opt_bb_predictions,
                    "local.predictions" = res$opt_local_predictions)

    if (is.null(predictions)) {
        stop(sprintf("%s estimates not available", predictions.type))
    }

    residuals <- res$y_sorted - predictions

    plot(res$x_sorted, residuals, type = "p",
         xlab = if (xlab == "") "x" else xlab,
         ylab = if (ylab == "") "Residuals" else ylab,
         main = if (title == "") "Residual Plot" else title,
         las = 1,
         ...)
    abline(h = 0, col = "gray", lty = 2)
}

plot.residuals.hist <- function(res, title, xlab, ylab, predictions.type = "predictions", ...) {
    predictions <- switch(predictions.type,
                    "predictions" = res$opt_predictions,
                    "bb.predictions" = res$opt_bb_predictions,
                    "local.predictions" = res$opt_local_predictions)

    if (is.null(predictions)) {
        stop(sprintf("%s estimates not available", predictions.type))
    }

    residuals <- res$y_sorted - predictions
    hist(residuals, br = 100, col = "red", las = 1,
         xlab = if (ylab == "") "Residuals" else xlab,
         main = if (title == "") "Residual Plot" else title,
         ...)
}

#' Generate Summary Statistics for Graph Path Linear Model Fits
#'
#' @description
#' Produces comprehensive summary statistics for graph path linear model (UPGMALOG) fits,
#' including model parameters, fit quality measures, cross-validation metrics,
#' residual analysis, and optional bootstrap statistics. When true values are available,
#' computes additional error metrics and normal distribution equivalence statistics.
#'
#' @param object A 'upgmalog' object from adaptive.nbhd.size.univariate.upgmalog
#' @param quantiles Numeric vector of probabilities for residual quantiles,
#'        default: c(0, 0.25, 0.5, 0.75, 1)
#' @param ... Additional arguments (currently unused)
#'
#' @details
#' Calculates multiple sets of statistics:
#' - Basic model information (sample size, optimal h, data range)
#' - Fit quality measures (MSE, RMSE, MAE)
#' - Cross-validation metrics
#' - Residual analysis
#' - True error analysis (when y.true is available)
#' - Bootstrap statistics (when bootstrap was performed)
#'
#' For true error analysis, computes equivalent normal standard deviation
#' and tests for normality using Shapiro-Wilk test.
#'
#' @return A 'summary.upgmalog' object containing:
#' \describe{
#'   \item{model_info}{Basic model parameters and dimensions}
#'   \item{fit_stats}{Fit quality measures}
#'   \item{cv_stats}{Cross-validation statistics}
#'   \item{residual_stats}{Residual analysis}
#'   \item{true_error_stats}{True error analysis (if available)}
#'   \item{bootstrap_info}{Bootstrap statistics (if available)}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#'
#' # Fit model and generate summary
#' fit <- adaptive.nbhd.size.univariate.upgmalog(x, y)
#' summary(fit)
#'
#' # With custom quantiles
#' summary(fit, quantiles = c(0.1, 0.5, 0.9))
#' }
#'
#' @seealso
#' \code{\link{adaptive.nbhd.size.univariate.upgmalog}},
#' \code{\link{plot.upgmalog}}
#'
#' @export
summary.upgmalog <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    if (!inherits(object, "upgmalog")) {
        stop("Input must be a 'upgmalog' object")
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("quantiles must be numeric values between 0 and 1")
    }
    if (any(is.na(object$y_sorted)) || any(is.na(object$opt_predictions))) {
        warning("Missing values detected in fit results")
    }
    if (length(object$h_cv_errors) != length(object$h_values)) {
        warning("Inconsistent cross-validation results detected")
    }

    # Calculate residuals
    residuals <- object$y_sorted - object$opt_predictions

    # Basic model information
    model_info <- list(
        n_observations = length(object$x_sorted),
        optimal_h = object$opt_h,
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

    # True error statistics and normal error equivalence
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted)) {
        true_residuals <- object$y_true_sorted - object$opt_predictions

        # Calculate various error metrics
        true_mse <- mean(true_residuals^2)
        true_rmse <- sqrt(true_mse)
        true_mae <- mean(abs(true_residuals))
        true_median_ae <- median(abs(true_residuals))

        # Calculate equivalent normal standard deviation
        # This is the SD that would give the same MSE if the errors were normally distributed
        equiv_normal_sd <- sqrt(true_mse)

        # Calculate normalized errors (should be approximately N(0,1) if normal)
        normalized_errors <- true_residuals / equiv_normal_sd

        # Test for normality using Shapiro-Wilk test
        normality_test <- shapiro.test(normalized_errors)

        # Calculate error quantiles for QQ plot comparison
        error_quantiles <- quantile(normalized_errors, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
        theoretical_quantiles <- qnorm(c(0.025, 0.25, 0.5, 0.75, 0.975))

        true_error_stats <- list(
            true_mse = true_mse,
            true_rmse = true_rmse,
            true_mae = true_mae,
            true_median_ae = true_median_ae,
            equiv_normal_sd = equiv_normal_sd,
            normalized_error_stats = list(
                mean = mean(normalized_errors),
                sd = sd(normalized_errors),
                quantiles = error_quantiles
            ),
            theoretical_normal_quantiles = theoretical_quantiles,
            shapiro_test = list(
                statistic = normality_test$statistic,
                p_value = normality_test$p.value
            ),
            relative_efficiency = true_mae / (sqrt(2/pi) * equiv_normal_sd)  # Relative to normal efficiency
        )
    }

    # Cross-validation statistics
    cv_stats <- list(
        ##mean_cv_error = object$mean_cv_error,
        min_cv_error = min(object$h_cv_errors),
        optimal_h_cv_error = object$h_cv_errors[which(object$h_values == object$opt_h)],
        h_range = range(object$h_values)
    )

    # Residual statistics
    residual_stats <- list(
        mean = mean(residuals),
        sd = sd(residuals),
        quantiles = quantile(residuals, probs = quantiles)
    )

    # Bootstrap information (if available)
    bootstrap_info <- NULL
    if (!is.null(object$opt_bb_predictions)) {
        bootstrap_info <- list(
            n_bootstrap = length(object$opt_bb_predictions),
            has_cri = !is.null(object$opt_ci_lower) && !is.null(object$opt_ci_upper),
            mean_cri_width = if (!is.null(object$opt_ci_upper) && !is.null(object$opt_ci_lower))
                mean(object$opt_ci_upper - object$opt_ci_lower) else NA
        )
    }

    # Create summary object
    result <- list(
        call = object$call,
        model_info = model_info,
        fit_stats = fit_stats,
        cv_stats = cv_stats,
        residual_stats = residual_stats,
        true_error_stats = true_error_stats,
        bootstrap_info = bootstrap_info
    )

    class(result) <- "summary.upgmalog"
    return(result)
}


#' Print Summary Statistics for Graph Path Linear Model Fits
#'
#' @description
#' Formats and displays summary statistics for UPGMALOG fits in a structured,
#' easy-to-read format. Output includes model parameters, fit statistics,
#' error analysis, and diagnostic information.
#'
#' @param x A 'summary.upgmalog' object from summary.upgmalog
#' @param digits Number of significant digits for numerical output (default: 4)
#' @param ... Additional arguments (currently unused)
#'
#' @details
#' Displays results in organized sections:
#' - Model Information
#' - Fit Statistics
#' - True Error Statistics (if available)
#' - Cross-validation Statistics
#' - Residual Analysis
#' - Bootstrap Information (if available)
#'
#' Each section uses clear formatting with headers and appropriate spacing
#' for readability.
#'
#' @return Returns x invisibly, following R print method convention
#'
#' @examples
#' \dontrun{
#' fit <- adaptive.nbhd.size.univariate.upgmalog(x, y)
#' fit_summary <- summary(fit)
#' print(fit_summary)
#'
#' # With more decimal places
#' print(fit_summary, digits = 6)
#' }
#'
#' @seealso
#' \code{\link{summary.upgmalog}}
#'
#' @export
print.summary.upgmalog <- function(x, digits = 4, ...) {
    ## Helper function to create separator lines
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
    section_header("ADAPTIVE NEIGHBORHOOD SIZE K-MEANS REGRESSION SUMMARY")

    # Model Information
    cat("\nModel Information:\n")
    cat(sprintf("- Number of observations:      %d\n", x$model_info$n_observations))
    cat(sprintf("- Optimal neighborhood size:   %d\n", x$model_info$optimal_h))
    cat(sprintf("- X range:                     [%.3f, %.3f]\n",
                x$model_info$min_x, x$model_info$max_x))
    cat(sprintf("- X span:                      %.3f\n", x$model_info$range_x))

    # Fit Statistics
    section_header("FIT STATISTICS")
    cat(sprintf("MSE:                           %.4f\n", x$fit_stats$mse))
    cat(sprintf("RMSE:                          %.4f\n", x$fit_stats$rmse))
    cat(sprintf("MAE:                           %.4f\n", x$fit_stats$mae))
    cat(sprintf("Median AE:                     %.4f\n", x$fit_stats$median_ae))

    # True Error Statistics (if available)
    if (!is.null(x$true_error_stats)) {
        section_header("TRUE ERROR STATISTICS")
        cat(sprintf("True MSE:                      %.4f\n", x$true_error_stats$true_mse))
        cat(sprintf("True RMSE:                     %.4f\n", x$true_error_stats$true_rmse))
        cat(sprintf("True MAE:                      %.4f\n", x$true_error_stats$true_mae))
        cat(sprintf("True Median AE:                %.4f\n", x$true_error_stats$true_median_ae))
        cat(sprintf("Equivalent Normal SD:          %.4f\n", x$true_error_stats$equiv_normal_sd))
        cat("\nNormalized Error Analysis:\n")
        cat(sprintf("- Mean:                        %.4f\n",
                   x$true_error_stats$normalized_error_stats$mean))
        cat(sprintf("- SD:                          %.4f\n",
                   x$true_error_stats$normalized_error_stats$sd))
        cat(sprintf("- Shapiro-Wilk p-value:       %.4f\n",
                   x$true_error_stats$shapiro_test$p_value))
        cat(sprintf("- Relative efficiency:         %.4f\n",
                   x$true_error_stats$relative_efficiency))

        cat("\nQuantile Comparison (Normalized Errors vs Normal):\n")
        quant_table <- data.frame(
            Probability = c("2.5%", "25%", "50%", "75%", "97.5%"),
            Observed = x$true_error_stats$normalized_error_stats$quantiles,
            Theoretical = x$true_error_stats$theoretical_normal_quantiles
        )
        print(format(quant_table, digits = 4), row.names = FALSE)
    }

    # Cross-validation Statistics
    section_header("CROSS-VALIDATION STATISTICS")
    if (!is.null(x$cv_stats$mean_cv_error)) {
        cat(sprintf("Mean CV Error:                 %.4f\n", x$cv_stats$mean_cv_error))
    }
    cat(sprintf("Minimum CV Error:              %.4f\n", x$cv_stats$min_cv_error))
    cat(sprintf("Optimal h CV Error:            %.4f\n", x$cv_stats$optimal_h_cv_error))
    cat(sprintf("h range:                       [%d, %d]\n",
                x$cv_stats$h_range[1], x$cv_stats$h_range[2]))

    # Residual Statistics
    section_header("RESIDUAL STATISTICS")
    cat(sprintf("Mean:                          %.4f\n", x$residual_stats$mean))
    cat(sprintf("Standard Deviation:            %.4f\n", x$residual_stats$sd))
    cat("\nQuantiles:\n")
    quant_df <- data.frame(
        Quantile = names(x$residual_stats$quantiles),
        Value = x$residual_stats$quantiles
    )
    print(format(quant_df, digits = 4), row.names = FALSE)

    # Bootstrap Information
    if (!is.null(x$bootstrap_info)) {
        section_header("BOOTSTRAP INFORMATION")
        cat(sprintf("Number of bootstrap samples:    %d\n", x$bootstrap_info$n_bootstrap))
        cat(sprintf("Credible intervals:            %s\n",
                   if(x$bootstrap_info$has_cri) "Available" else "Not available"))
        if (x$bootstrap_info$has_cri) {
            cat(sprintf("Mean CI width:                %.4f\n",
                       x$bootstrap_info$mean_cri_width))
        }
    }

    cat("\n")  # Final newline for spacing
    invisible(x)
}



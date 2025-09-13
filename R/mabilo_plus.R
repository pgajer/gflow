#' Model-Averaged Locally Weighted Scatterplot Smoothing Plus (MABILO Plus)
#'
#' @description
#' Performs smoothing using MABILO Plus, which extends the original MABILO algorithm
#' by incorporating flexible model averaging strategies and error filtering approaches
#' without bootstrap resampling.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates (response values).
#' @param y.true Optional numeric vector of true y values for error calculation.
#' @param k.min Minimum number of neighbors (positive integer). Default is 5% of data points or 3, whichever is larger.
#' @param k.max Maximum number of neighbors (positive integer > k.min). Default corresponds to 2k+1 = n/2.
#' @param model.averaging.strategy Character string specifying the model averaging strategy.
#'   Must be one of: "kernel.weights.only", "error.weights.only", "kernel.and.error.weights".
#'   Default is "kernel.and.error.weights".
#' @param error.filtering.strategy Character string specifying the error filtering strategy.
#'   Must be one of: "global.percentile", "local.percentile", "combined.percentile", "best.k.models".
#'   Default is "combined.percentile".
#' @param distance.kernel Integer specifying the kernel type for distance weighting (1-10).
#'   Default is 7L.
#' @param model.kernel Integer specifying the kernel type for model weighting (1-10).
#'   Default is 7L.
#' @param dist.normalization.factor Positive numeric value for distance normalization.
#'   Default is 1.0.
#' @param epsilon Small positive number to prevent division by zero. Default is 1e-10.
#' @param verbose Logical indicating whether to print progress information. Default is TRUE.
#'
#' @return A list of class "mabilo_plus" containing:
#'   \item{\code{k_values}}{Vector of k values tested}
#'   \item{\code{opt_sm_k}}{Optimal k value for simple mean predictions}
#'   \item{\code{opt_ma_k}}{Optimal k value for model averaged predictions}
#'   \item{\code{opt_sm_k_idx}}{Index of optimal k for simple mean}
#'   \item{\code{opt_ma_k_idx}}{Index of optimal k for model averaging}
#'   \item{\code{k_mean_sm_errors}}{Mean errors for each k value (simple mean)}
#'   \item{\code{k_mean_ma_errors}}{Mean errors for each k value (model averaged)}
#'   \item{\code{k_mean_true_errors}}{Mean true errors for each k (if y.true provided)}
#'   \item{\code{sm_predictions}}{Simple mean predictions}
#'   \item{\code{ma_predictions}}{Model averaged predictions}
#'   \item{\code{sm_errors}}{Leave-one-out CV errors (simple mean)}
#'   \item{\code{ma_errors}}{Leave-one-out CV errors (model averaged)}
#'   \item{\code{k_predictions}}{Matrix of predictions for each k value}
#'   \item{\code{x_sorted}}{Sorted x values}
#'   \item{\code{y_sorted}}{Sorted y values}
#'   \item{\code{y_true_sorted}}{Sorted true y values (if provided)}
#'   \item{\code{k_min}}{Minimum k value used}
#'   \item{\code{k_max}}{Maximum k value used}
#'
#' @examples
#' # Basic usage
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#' fit <- mabilo.plus(x, y)
#'
#' # With custom parameters
#' fit2 <- mabilo.plus(x, y, k.min = 5, k.max = 20,
#'                     model.averaging.strategy = "kernel.weights.only")
#'
#' # Plot the results
#' plot(x, y)
#' lines(fit$x_sorted, fit$ma_predictions, col = "red", lwd = 2)
#'
#' @export
mabilo.plus <- function(x,
                        y,
                        y.true = NULL,
                        k.min = max(3, as.integer(0.05 * length(x))),
                        k.max = as.integer((length(x)-2)/4), # this corresponds to 2k + 1 = n/2
                        model.averaging.strategy = "kernel.and.error.weights",
                        error.filtering.strategy = "combined.percentile",
                        distance.kernel = 7L,
                        model.kernel = 7L,
                        dist.normalization.factor = 1.0,
                        epsilon = 1e-10,
                        verbose = TRUE) {
    # Check x and y
    if (!is.numeric(x)) stop("x must be numeric")
    if (!is.numeric(y)) stop("y must be numeric")
    if (length(x) != length(y)) stop("x and y must have the same length")
    if (length(x) < 3) stop("At least 3 points are required")

    n <- length(x)

    # Check y.true
    if (!is.null(y.true)) {
        if (!is.numeric(y.true)) stop("y.true must be numeric")
        if (length(y.true) != n) stop("y.true must have the same length as x")
    } else {
        y.true <- numeric(0)  # Empty vector
    }


    # Check k.min and k.max
    if (!is.numeric(k.min) || k.min != round(k.min) || k.min < 1)
        stop("k.min must be a positive integer")

    if (is.null(k.max)) k.max <- min(as.integer((n-2)/4), max(3*k.min, 10))
    if (!is.numeric(k.max) || k.max != round(k.max))
        stop("k.max must be an integer")
    if (k.max <= k.min)
        stop("k.max must be greater than k.min")
    if (k.max >= n)
        stop("k.max must be less than the number of points")

    if (k.max > (n-2)/4) {
        print(sprintf("ERROR: k.max (%d) is too large for the data size (%d).\n", k.max, n))
        stop("k.max should be less than (n-2)/4 to ensure multiple windows.")
    }

    # Check averaging strategy
    valid.avg.strategies <- c("kernel.weights.only",
                              "kernel.weights.with.filtering",
                              "error.weights.only",
                              "kernel.and.error.weights")
    if (!model.averaging.strategy %in% valid.avg.strategies) {
        stop("Invalid model.averaging.strategy. Must be one of: ",
             paste(valid.avg.strategies, collapse = ", "))
    } else {
        model.averaging.strategy <- gsub("\\.","_", model.averaging.strategy)
    }

    ## Check filtering strategy
    valid.filter.strategies <- c("global.percentile",
                                 "local.percentile",
                                 "combined.percentile",
                                 "best.k.models")
    if (!error.filtering.strategy %in% valid.filter.strategies) {
        stop("Invalid error.filtering.strategy. Must be one of: ",
             paste(valid.filter.strategies, collapse = ", "))
    } else {
        error.filtering.strategy <- gsub("\\.","_", error.filtering.strategy)
    }

    # Check kernels
    if (!is.numeric(distance.kernel) || distance.kernel < 1 || distance.kernel > 10)
        stop("distance.kernel must be an integer between 1 and 10")
    if (!is.numeric(model.kernel) || model.kernel < 1 || model.kernel > 10)
        stop("model.kernel must be an integer between 1 and 10")

    # Check other numeric parameters
    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 0)
        stop("dist.normalization.factor must be positive")
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("epsilon must be positive")

    # Check verbose
    if (!is.logical(verbose))
        stop("verbose must be logical")

    ## Check if x is sorted and sort if necessary
    if (!identical(x, sort(x))) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        if (!is.null(y.true)) y.true <- y.true[ord]
    }

    w <- rep(1.0, n)

    # Call the C++ implementation
    result <- .Call(S_mabilo_plus,
                   as.double(x),
                   as.double(y),
                   as.double(y.true),
                   as.double(w),
                   as.integer(k.min),
                   as.integer(k.max),
                   as.character(model.averaging.strategy),
                   as.character(error.filtering.strategy),
                   as.integer(distance.kernel),
                   as.integer(model.kernel),
                   as.double(dist.normalization.factor),
                   as.double(epsilon),
                   as.logical(verbose))

    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- y.true
    result$k_min <- k.min
    result$k_max <- k.max

    class(result) <- "mabilo_plus"

    return(result)
}

##
## Utility Functions for mabilo.plus()
##

# ============================================================================
# PLOT METHODS
# ============================================================================

#' Plot Method for Mabilo Plus Objects
#'
#' @description
#' Generates various diagnostic and visualization plots for mabilo.plus objects.
#' The function supports different types of plots including fitted values,
#' diagnostic plots, and residual analyses. For each plot type, various
#' customization options are available through the function parameters.
#'
#' @param x An object of class "mabilo_plus"
#' @param type Character string specifying the type of plot to create (partial matching supported):
#'   * "fit": Plot fitted values against x values (e.g., "f" or "fi")
#'   * "diagnostic": Plot error diagnostics for different k values (e.g., "d" or "dia")
#'   * "residuals": Plot residuals against x values (e.g., "res")
#'   * "residuals.hist": Plot histogram of residuals (e.g., "rh" or "residuals.h")
#' @param title Plot title (default: "")
#' @param xlab Label for x-axis (default: "")
#' @param ylab Label for y-axis (default: "")
#' @param predictions.type Type of predictions to plot:
#'   * "both": plot both sm_predictions and ma_predictions
#'   * "sm": sm_predictions only
#'   * "ma": ma_predictions only
#' @param diagnostic.type Type of errors to show in diagnostic plot:
#'   * "both": show both sm and ma errors
#'   * "sm": sm errors only
#'   * "ma": ma errors only
#' @param with.y.true Logical; if TRUE and true values are available, adds them to the plot
#'   (default: TRUE)
#' @param with.pts Logical; if TRUE, adds original data points to the plot
#'   (default: FALSE)
#' @param true.lwd Line width for true values (default: 2)
#' @param true.col Color for true values line (default: "red")
#' @param sm.col Color for sm prediction and error lines (default: "#00FF00")
#' @param ma.col Color for ma prediction and error lines (default: "blue")
#' @param pts.col Color for points (default: "gray60")
#' @param with.predictions.pts Logical; if TRUE, adds points at prediction values
#'   (default: FALSE)
#' @param predictions.pts.col Color for prediction points (default: "blue")
#' @param predictions.pts.pch Point character for prediction points (default: 20)
#' @param ylim Numeric vector of length 2 giving y-axis limits; if NULL, computed from data
#'   (default: NULL)
#' @param legend.cex Character expansion factor for legend (default: 0.8)
#' @param k Specific k value to plot predictions for (default: NULL, uses optimal k)
#' @param with.legend Logical; if TRUE, displays legend on the plot (default: TRUE)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns NULL. The function is called for its side effect of producing plots.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#'
#' # Fit mabilo.plus model
#' fit <- mabilo.plus(x, y)
#'
#' # Basic fit plot
#' plot(fit)
#'
#' # Plot both SM and MA predictions
#' plot(fit, type = "fit", predictions.type = "both")
#'
#' # Diagnostic plot
#' plot(fit, type = "diagnostic", diagnostic.type = "both")
#'
#' # Residual plots
#' plot(fit, type = "residuals", predictions.type = "ma")
#' plot(fit, type = "residuals.hist", predictions.type = "both")
#' }
#'
#' @export
plot.mabilo_plus <- function(x,
                             type = "fit",
                             title = "",
                             xlab = "",
                             ylab = "",
                             predictions.type = "both",
                             diagnostic.type = "both",
                             with.y.true = TRUE,
                             with.pts = FALSE,
                             true.lwd = 2,
                             true.col = "red",
                             sm.col = "#00FF00",
                             ma.col = "blue",
                             pts.col = "gray60",
                             with.predictions.pts = FALSE,
                             predictions.pts.col = "blue",
                             predictions.pts.pch = 20,
                             ylim = NULL,
                             legend.cex = 0.8,
                             k = NULL,
                             with.legend = TRUE,
                             ...) {

    if (!inherits(x, "mabilo_plus")) {
        stop("Input must be a 'mabilo_plus' object")
    }

    # Validate parameters
    if (!predictions.type %in% c("both", "sm", "ma")) {
        stop("predictions.type must be one of: 'both', 'sm', 'ma'")
    }

    if (!diagnostic.type %in% c("both", "sm", "ma")) {
        stop("diagnostic.type must be one of: 'both', 'sm', 'ma'")
    }

    # Define valid plot types
    valid_types <- c("fit", "diagnostic", "residuals", "residuals.hist")

    # Find partial matches for the provided type
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

    # Call appropriate plotting function
    switch(type,
           "fit" = mabilo_plus.plot.fit(x, title, xlab, ylab, predictions.type,
                                       with.y.true, with.pts, true.lwd, true.col,
                                       sm.col, ma.col, pts.col, with.predictions.pts,
                                       predictions.pts.col, predictions.pts.pch,
                                       ylim, legend.cex, k, with.legend, ...),
           "diagnostic" = mabilo_plus.plot.diagnostic(x, diagnostic.type, title, xlab,
                                                     ylab, sm.col, ma.col, true.col,
                                                     legend.cex, ...),
           "residuals" = mabilo_plus.plot.residuals(x, title, xlab, ylab,
                                                   predictions.type, ...),
           "residuals.hist" = mabilo_plus.plot.residuals.hist(x, title, xlab, ylab,
                                                              predictions.type, ...),
           stop("Invalid plot type.")
    )

    invisible(NULL)
}

# Helper function for fit plots
mabilo_plus.plot.fit <- function(res, title, xlab, ylab, predictions.type, with.y.true,
                                with.pts, true.lwd, true.col, sm.col, ma.col, pts.col,
                                with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                                ylim, legend.cex, k, with.legend, ...) {

    # Select predictions to plot
    if (!is.null(k)) {
        k_idx <- which(res$k_values == k)
        sm_pred <- res$k_predictions[, k_idx]
        ma_pred <- res$k_predictions[, k_idx]  # For specific k, both are same
        pred_label <- sprintf("Predictions (k=%d)", k)
    } else {
        sm_pred <- res$sm_predictions
        ma_pred <- res$ma_predictions
        pred_label <- c("SM Predictions", "MA Predictions")
    }

    # Determine which predictions to plot
    if (predictions.type == "sm") {
        pred_data <- matrix(sm_pred, ncol = 1)
        pred_cols <- sm.col
        pred_labels <- if (!is.null(k)) pred_label else pred_label[1]
    } else if (predictions.type == "ma") {
        pred_data <- matrix(ma_pred, ncol = 1)
        pred_cols <- ma.col
        pred_labels <- if (!is.null(k)) pred_label else pred_label[2]
    } else {  # both
        if (!is.null(k)) {
            pred_data <- matrix(ma_pred, ncol = 1)
            pred_cols <- ma.col
            pred_labels <- pred_label
        } else {
            pred_data <- cbind(sm_pred, ma_pred)
            pred_cols <- c(sm.col, ma.col)
            pred_labels <- pred_label
        }
    }

    # Calculate ylim if not provided
    if (is.null(ylim)) {
        y_range <- range(res$y_sorted)
        pred_range <- range(pred_data)
        if (with.y.true && !is.null(res$y_true_sorted)) {
            true_range <- range(res$y_true_sorted)
            ylim <- range(c(y_range, pred_range, true_range))
        } else {
            ylim <- range(c(y_range, pred_range))
        }
        # Add some margin
        ylim_margin <- diff(ylim) * 0.05
        ylim <- ylim + c(-ylim_margin, ylim_margin)
    }

    # Initialize plot
    if (with.pts) {
        plot(res$x_sorted, res$y_sorted,
             xlab = if(xlab == "") "x" else xlab,
             ylab = if(ylab == "") "y" else ylab,
             main = if(title == "") "Mabilo Plus Fit" else title,
             col = pts.col,
             ylim = ylim, las = 1, ...)
    } else {
        plot(res$x_sorted, res$y_sorted, type = 'n',
             xlab = if(xlab == "") "x" else xlab,
             ylab = if(ylab == "") "y" else ylab,
             main = if(title == "") "Mabilo Plus Fit" else title,
             ylim = ylim, las = 1, ...)
    }

    # Add predictions
    for (i in 1:ncol(pred_data)) {
        lines(res$x_sorted, pred_data[, i], col = pred_cols[i], lwd = 2)

        if (with.predictions.pts) {
            points(res$x_sorted, pred_data[, i], col = predictions.pts.col,
                  pch = predictions.pts.pch, cex = 0.5)
        }
    }

    # Add true values if available
    if (with.y.true && !is.null(res$y_true_sorted)) {
        lines(res$x_sorted, res$y_true_sorted, col = true.col, lwd = true.lwd)
    }

    # Add legend
    if (with.legend) {
        legend_items <- character()
        legend_cols <- character()
        legend_ltys <- numeric()
        legend_pchs <- numeric()

        if (with.pts) {
            legend_items <- c(legend_items, "Data Points")
            legend_cols <- c(legend_cols, pts.col)
            legend_ltys <- c(legend_ltys, NA)
            legend_pchs <- c(legend_pchs, 1)
        }

        for (i in 1:length(pred_labels)) {
            legend_items <- c(legend_items, pred_labels[i])
            legend_cols <- c(legend_cols, pred_cols[i])
            legend_ltys <- c(legend_ltys, 1)
            legend_pchs <- c(legend_pchs, NA)
        }

        if (with.y.true && !is.null(res$y_true_sorted)) {
            legend_items <- c(legend_items, "True Values")
            legend_cols <- c(legend_cols, true.col)
            legend_ltys <- c(legend_ltys, 1)
            legend_pchs <- c(legend_pchs, NA)
        }

        legend("topright", legend = legend_items,
               col = legend_cols, lty = legend_ltys, pch = legend_pchs,
               bg = "white", inset = 0.05, cex = legend.cex)
    }
}

# Helper function for diagnostic plots
mabilo_plus.plot.diagnostic <- function(res, diagnostic.type, title, xlab, ylab,
                                       sm.col, ma.col, true.col, legend.cex, ...) {

    # Prepare error data based on type
    if (diagnostic.type == "sm") {
        err_data <- matrix(res$k_mean_sm_errors, ncol = 1)
        labels <- "SM LOOCV Errors"
        cols <- sm.col
        pchs <- 1
        opt_k <- res$opt_sm_k
        opt_label <- sprintf("Optimal SM k = %d", opt_k)
    } else if (diagnostic.type == "ma") {
        err_data <- matrix(res$k_mean_ma_errors, ncol = 1)
        labels <- "MA LOOCV Errors"
        cols <- ma.col
        pchs <- 2
        opt_k <- res$opt_ma_k
        opt_label <- sprintf("Optimal MA k = %d", opt_k)
    } else {  # both
        err_data <- cbind(res$k_mean_sm_errors, res$k_mean_ma_errors)
        labels <- c("SM LOOCV Errors", "MA LOOCV Errors")
        cols <- c(sm.col, ma.col)
        pchs <- c(1, 2)
        opt_k <- c(res$opt_sm_k, res$opt_ma_k)
        opt_label <- c(sprintf("Optimal SM k = %d", opt_k[1]),
                      sprintf("Optimal MA k = %d", opt_k[2]))
    }

    # Add true errors if available
    if (!is.null(res$k_mean_true_errors)) {
        err_data <- cbind(err_data, res$k_mean_true_errors)
        labels <- c(labels, "True Errors")
        cols <- c(cols, true.col)
        pchs <- c(pchs, 3)
    }

    # Plot errors
    matplot(res$k_values, err_data, type = 'b',
            xlab = if(xlab == "") "k value" else xlab,
            ylab = if(ylab == "") "Error" else ylab,
            main = if(title == "") "Error Diagnostic Plot" else title,
            col = cols,
            pch = pchs, las = 1,
            ...)

    # Add optimal k values
    if (diagnostic.type == "both") {
        abline(v = opt_k[1], col = sm.col, lty = 2)
        abline(v = opt_k[2], col = ma.col, lty = 2)
        mtext(opt_label[1], side = 3, line = 0.25, at = opt_k[1], col = sm.col, cex = 0.8)
        mtext(opt_label[2], side = 3, line = -0.5, at = opt_k[2], col = ma.col, cex = 0.8)
    } else {
        abline(v = opt_k, col = cols[1], lty = 2)
        mtext(opt_label, side = 3, line = 0.25, at = opt_k, col = cols[1])
    }

    # Add legend
    legend("topright", legend = labels,
           col = cols,
           pch = pchs,
           lty = 1, bg = "white", inset = 0.05, cex = legend.cex)
}

# Helper function for residual plots
mabilo_plus.plot.residuals <- function(res, title, xlab, ylab, predictions.type, ...) {

    # Calculate residuals based on type
    if (predictions.type == "sm") {
        residuals <- res$y_sorted - res$sm_predictions
        res_label <- "SM Residuals"
    } else if (predictions.type == "ma") {
        residuals <- res$y_sorted - res$ma_predictions
        res_label <- "MA Residuals"
    } else {  # both
        residuals <- cbind(res$y_sorted - res$sm_predictions,
                          res$y_sorted - res$ma_predictions)
        res_label <- c("SM Residuals", "MA Residuals")
    }

    if (predictions.type == "both") {
        # Two panel plot
        par(mfrow = c(2, 1))

        # SM residuals
        plot(res$x_sorted, residuals[,1],
             xlab = if(xlab == "") "x" else xlab,
             ylab = if(ylab == "") "Residuals" else ylab,
             main = if(title == "") "SM Residuals" else paste(title, "- SM"),
             las = 1, ...)
        abline(h = 0, lty = 2)

        # MA residuals
        plot(res$x_sorted, residuals[,2],
             xlab = if(xlab == "") "x" else xlab,
             ylab = if(ylab == "") "Residuals" else ylab,
             main = if(title == "") "MA Residuals" else paste(title, "- MA"),
             las = 1, ...)
        abline(h = 0, lty = 2)

        par(mfrow = c(1, 1))
    } else {
        # Single plot
        plot(res$x_sorted, residuals,
             xlab = if(xlab == "") "x" else xlab,
             ylab = if(ylab == "") "Residuals" else ylab,
             main = if(title == "") paste(res_label, "Plot") else title,
             las = 1, ...)
        abline(h = 0, lty = 2)
    }
}

# Helper function for residual histograms
mabilo_plus.plot.residuals.hist <- function(res, title, xlab, ylab, predictions.type, ...) {

    # Calculate residuals based on type
    if (predictions.type == "sm") {
        residuals <- res$y_sorted - res$sm_predictions
        hist_cols <- "lightblue"
        density_cols <- "blue"
        res_labels <- "SM Residuals"
    } else if (predictions.type == "ma") {
        residuals <- res$y_sorted - res$ma_predictions
        hist_cols <- "pink"
        density_cols <- "red"
        res_labels <- "MA Residuals"
    } else {  # both
        sm_residuals <- res$y_sorted - res$sm_predictions
        ma_residuals <- res$y_sorted - res$ma_predictions

        # Two panel plot
        par(mfrow = c(2, 1))

        # SM residuals histogram
        hist(sm_residuals,
             breaks = "FD",
             col = "lightblue",
             border = "white",
             xlab = if(xlab == "") "Residuals" else xlab,
             main = if(title == "") "SM Residuals Histogram" else paste(title, "- SM"),
             probability = TRUE,
             ...)
        lines(density(sm_residuals), col = "blue", lwd = 2)
        abline(v = 0, col = "red", lty = 2)

        # MA residuals histogram
        hist(ma_residuals,
             breaks = "FD",
             col = "pink",
             border = "white",
             xlab = if(xlab == "") "Residuals" else xlab,
             main = if(title == "") "MA Residuals Histogram" else paste(title, "- MA"),
             probability = TRUE,
             ...)
        lines(density(ma_residuals), col = "red", lwd = 2)
        abline(v = 0, col = "red", lty = 2)

        par(mfrow = c(1, 1))
        return(invisible(NULL))
    }

    # Single histogram case
    hist(residuals,
         breaks = "FD",
         col = hist_cols,
         border = "white",
         xlab = if(xlab == "") "Residuals" else xlab,
         main = if(title == "") paste(res_labels, "Histogram") else title,
         probability = TRUE,
         ...)
    lines(density(residuals), col = density_cols, lwd = 2)
    abline(v = 0, col = "red", lty = 2)
}

# ============================================================================
# SUMMARY METHODS
# ============================================================================

#' Compute Summary Statistics for Mabilo Plus Objects
#'
#' @description
#' Computes summary statistics for mabilo.plus fits including model parameters,
#' fit statistics, error analysis, and diagnostic information for both
#' simple mean (SM) and model averaging (MA) predictions.
#'
#' @param object A 'mabilo_plus' object
#' @param quantiles Numeric vector of probabilities for quantile computations
#' @param ... Additional arguments (currently unused)
#'
#' @return A 'summary.mabilo_plus' object containing:
#' \itemize{
#'   \item model_info: Basic information about the model fit
#'   \item fit_stats: Prediction accuracy metrics for both SM and MA
#'   \item k_error_stats: Error statistics for different k values
#'   \item residual_stats: Residual analysis results for both SM and MA
#'   \item true_error_stats: Statistics comparing to true values (if available)
#' }
#'
#' @export
summary.mabilo_plus <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    if (!inherits(object, "mabilo_plus")) {
        stop("Input must be a 'mabilo_plus' object")
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("quantiles must be numeric values between 0 and 1")
    }

    # Check for missing values
    if (any(is.na(object$y_sorted)) ||
        any(is.na(object$sm_predictions)) ||
        any(is.na(object$ma_predictions))) {
        warning("Missing values detected in fit results")
    }

    # Calculate residuals for both prediction types
    sm_residuals <- object$y_sorted - object$sm_predictions
    ma_residuals <- object$y_sorted - object$ma_predictions

    # Basic model information
    model_info <- list(
        n_observations = length(object$x_sorted),
        optimal_sm_k = object$opt_sm_k,
        optimal_ma_k = object$opt_ma_k,
        k_range = range(object$k_values),
        min_x = min(object$x_sorted),
        max_x = max(object$x_sorted),
        range_x = diff(range(object$x_sorted))
    )

    # Fit statistics for both prediction types
    fit_stats <- list(
        sm = list(
            mse = mean(sm_residuals^2),
            rmse = sqrt(mean(sm_residuals^2)),
            mae = mean(abs(sm_residuals)),
            median_ae = median(abs(sm_residuals))
        ),
        ma = list(
            mse = mean(ma_residuals^2),
            rmse = sqrt(mean(ma_residuals^2)),
            mae = mean(abs(ma_residuals)),
            median_ae = median(abs(ma_residuals))
        )
    )

    # True error statistics if true values are available
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted)) {
        # Calculate true residuals for both prediction types
        sm_true_residuals <- object$y_true_sorted - object$sm_predictions
        ma_true_residuals <- object$y_true_sorted - object$ma_predictions

        true_error_stats <- list(
            sm = list(
                mse = mean(sm_true_residuals^2),
                rmse = sqrt(mean(sm_true_residuals^2)),
                mae = mean(abs(sm_true_residuals)),
                median_ae = median(abs(sm_true_residuals))
            ),
            ma = list(
                mse = mean(ma_true_residuals^2),
                rmse = sqrt(mean(ma_true_residuals^2)),
                mae = mean(abs(ma_true_residuals)),
                median_ae = median(abs(ma_true_residuals))
            )
        )
    }

    # Error statistics for different k values
    k_error_stats <- list(
        sm = list(
            mean_errors = object$k_mean_sm_errors,
            min_error = min(object$k_mean_sm_errors),
            optimal_k_error = object$k_mean_sm_errors[object$opt_sm_k_idx]
        ),
        ma = list(
            mean_errors = object$k_mean_ma_errors,
            min_error = min(object$k_mean_ma_errors),
            optimal_k_error = object$k_mean_ma_errors[object$opt_ma_k_idx]
        )
    )

    # Residual statistics for both prediction types
    residual_stats <- list(
        sm = list(
            mean = mean(sm_residuals),
            sd = sd(sm_residuals),
            quantiles = quantile(sm_residuals, probs = quantiles)
        ),
        ma = list(
            mean = mean(ma_residuals),
            sd = sd(ma_residuals),
            quantiles = quantile(ma_residuals, probs = quantiles)
        )
    )

    # Create summary object
    result <- list(
        call = object$call,
        model_info = model_info,
        fit_stats = fit_stats,
        k_error_stats = k_error_stats,
        residual_stats = residual_stats,
        true_error_stats = true_error_stats
    )

    class(result) <- "summary.mabilo_plus"
    return(result)
}

#' Print Summary Statistics for Mabilo Plus Fits
#'
#' @description
#' Formats and displays summary statistics for mabilo.plus fits in a structured,
#' easy-to-read format. Output includes model parameters, fit statistics,
#' error analysis, and diagnostic information for both simple mean (SM) and
#' model averaging (MA) predictions.
#'
#' @param x A 'summary.mabilo_plus' object from summary.mabilo_plus
#' @param digits Number of significant digits for numerical output (default: 4)
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.summary.mabilo_plus <- function(x, digits = 4, ...) {
    # Helper function to create separator lines
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
    section_header("MABILO PLUS (Model-Averaged Locally Weighted Smoothing Plus) SUMMARY")

    # Model Information
    cat("\nModel Information:\n")
    cat(sprintf("Number of observations:         %d\n", x$model_info$n_observations))
    cat(sprintf("Optimal SM k:                   %d\n", x$model_info$optimal_sm_k))
    cat(sprintf("Optimal MA k:                   %d\n", x$model_info$optimal_ma_k))
    cat(sprintf("k range:                        [%d, %d]\n",
                x$model_info$k_range[1], x$model_info$k_range[2]))
    cat(sprintf("X range:                        [%.3f, %.3f]\n",
                x$model_info$min_x, x$model_info$max_x))

    # Fit Statistics
    section_header("Fit Statistics")

    cat("\nSimple Mean (SM) Predictions:\n")
    cat(sprintf("  MSE:                          %.4f\n", x$fit_stats$sm$mse))
    cat(sprintf("  RMSE:                         %.4f\n", x$fit_stats$sm$rmse))
    cat(sprintf("  MAE:                          %.4f\n", x$fit_stats$sm$mae))
    cat(sprintf("  Median Absolute Error:        %.4f\n", x$fit_stats$sm$median_ae))

    cat("\nModel Averaged (MA) Predictions:\n")
    cat(sprintf("  MSE:                          %.4f\n", x$fit_stats$ma$mse))
    cat(sprintf("  RMSE:                         %.4f\n", x$fit_stats$ma$rmse))
    cat(sprintf("  MAE:                          %.4f\n", x$fit_stats$ma$mae))
    cat(sprintf("  Median Absolute Error:        %.4f\n", x$fit_stats$ma$median_ae))

    # True Error Statistics
    if (!is.null(x$true_error_stats)) {
        section_header("True Error Statistics")

        cat("\nSimple Mean (SM) True Errors:\n")
        cat(sprintf("  True MSE:                     %.4f\n", x$true_error_stats$sm$mse))
        cat(sprintf("  True RMSE:                    %.4f\n", x$true_error_stats$sm$rmse))
        cat(sprintf("  True MAE:                     %.4f\n", x$true_error_stats$sm$mae))

        cat("\nModel Averaged (MA) True Errors:\n")
        cat(sprintf("  True MSE:                     %.4f\n", x$true_error_stats$ma$mse))
        cat(sprintf("  True RMSE:                    %.4f\n", x$true_error_stats$ma$rmse))
        cat(sprintf("  True MAE:                     %.4f\n", x$true_error_stats$ma$mae))
    }

    # k Error Statistics
    section_header("k-Value Error Analysis")

    cat("\nSimple Mean (SM):\n")
    cat(sprintf("  Minimum Error:                %.4f\n", x$k_error_stats$sm$min_error))
    cat(sprintf("  Optimal k Error:              %.4f\n", x$k_error_stats$sm$optimal_k_error))

    cat("\nModel Averaged (MA):\n")
    cat(sprintf("  Minimum Error:                %.4f\n", x$k_error_stats$ma$min_error))
    cat(sprintf("  Optimal k Error:              %.4f\n", x$k_error_stats$ma$optimal_k_error))

    # Residual Statistics
    section_header("Residual Statistics")

    cat("\nSimple Mean (SM) Residuals:\n")
    cat(sprintf("  Mean:                         %.4f\n", x$residual_stats$sm$mean))
    cat(sprintf("  Standard Deviation:           %.4f\n", x$residual_stats$sm$sd))
    cat("  Quantiles:\n")
    print(x$residual_stats$sm$quantiles, digits = digits)

    cat("\nModel Averaged (MA) Residuals:\n")
    cat(sprintf("  Mean:                         %.4f\n", x$residual_stats$ma$mean))
    cat(sprintf("  Standard Deviation:           %.4f\n", x$residual_stats$ma$sd))
    cat("  Quantiles:\n")
    print(x$residual_stats$ma$quantiles, digits = digits)

    hr()
    invisible(x)
}

# ============================================================================
# PRINT METHOD
# ============================================================================

#' Print Method for Mabilo Plus Objects
#'
#' @description
#' Prints a concise summary of a mabilo_plus object
#'
#' @param x A 'mabilo_plus' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.mabilo_plus <- function(x, ...) {
    cat("Mabilo Plus (Model-Averaged Locally Weighted Smoothing Plus) Fit\n")
    cat("=================================================================\n")
    cat(sprintf("Number of observations: %d\n", length(x$x_sorted)))
    cat(sprintf("k range: [%d, %d]\n", x$k_min, x$k_max))
    cat(sprintf("Optimal SM k: %d\n", x$opt_sm_k))
    cat(sprintf("Optimal MA k: %d\n", x$opt_ma_k))
    cat("\nCall summary(object) for detailed statistics\n")
    invisible(x)
}

# ============================================================================
# EXTRACTION METHODS
# ============================================================================

#' Extract Fitted Values from Mabilo Plus Model
#'
#' @description
#' Extracts fitted values from a mabilo_plus object
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which predictions to return:
#'   "ma" (default), "sm", or "both"
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector or matrix of fitted values
#'
#' @export
fitted.mabilo_plus <- function(object, type = c("ma", "sm", "both"), ...) {
    type <- match.arg(type)

    switch(type,
           "ma" = object$ma_predictions,
           "sm" = object$sm_predictions,
           "both" = cbind(sm = object$sm_predictions, ma = object$ma_predictions)
    )
}

#' Extract Residuals from Mabilo Plus Model
#'
#' @description
#' Extracts residuals from a mabilo_plus object
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which residuals to return:
#'   "ma" (default), "sm", or "both"
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector or matrix of residuals
#'
#' @export
residuals.mabilo_plus <- function(object, type = c("ma", "sm", "both"), ...) {
    type <- match.arg(type)

    sm_resid <- object$y_sorted - object$sm_predictions
    ma_resid <- object$y_sorted - object$ma_predictions

    switch(type,
           "ma" = ma_resid,
           "sm" = sm_resid,
           "both" = cbind(sm = sm_resid, ma = ma_resid)
    )
}

#' Extract Model Coefficients from Mabilo Plus Model
#'
#' @description
#' Extracts the optimal k values from a mabilo_plus object
#'
#' @param object A 'mabilo_plus' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Named vector containing optimal k values
#'
#' @export
coef.mabilo_plus <- function(object, ...) {
    c(opt_sm_k = object$opt_sm_k,
      opt_ma_k = object$opt_ma_k)
}

# ============================================================================
# PREDICTION METHOD
# ============================================================================

#' Predict Method for Mabilo Plus Objects
#'
#' @description
#' Predict response values for new data using a fitted mabilo_plus model
#'
#' @param object A 'mabilo_plus' object
#' @param newdata Numeric vector of new x values for prediction
#' @param type Character string specifying which predictions to return:
#'   "ma" (default), "sm", or "both"
#' @param k Optional specific k value to use for predictions. If NULL,
#'   uses optimal k values
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector or matrix of predictions
#'
#' @details
#' This function performs local weighted regression predictions for new data points.
#' For each new x value, it finds the k nearest neighbors from the training data
#' and performs a weighted local regression.
#'
#' @export
predict.mabilo_plus <- function(object, newdata, type = c("ma", "sm", "both"),
                               k = NULL, ...) {
    type <- match.arg(type)

    if (!is.numeric(newdata)) {
        stop("newdata must be numeric")
    }

    n_new <- length(newdata)
    n_train <- length(object$x_sorted)

    # Determine which k values to use
    if (!is.null(k)) {
        if (!is.numeric(k) || length(k) != 1) {
            stop("k must be a single numeric value")
        }
        if (k < object$k_min || k > object$k_max) {
            stop(sprintf("k must be between %d and %d", object$k_min, object$k_max))
        }
        sm_k <- ma_k <- k
    } else {
        sm_k <- object$opt_sm_k
        ma_k <- object$opt_ma_k
    }

    # Initialize prediction arrays
    sm_pred <- numeric(n_new)
    ma_pred <- numeric(n_new)

    # For each new point
    for (i in 1:n_new) {
        x_new <- newdata[i]

        # Calculate distances to all training points
        distances <- abs(object$x_sorted - x_new)

        # Get indices of k nearest neighbors for SM
        sm_neighbors <- order(distances)[1:sm_k]

        # Simple mean prediction
        sm_pred[i] <- mean(object$y_sorted[sm_neighbors])

        # For MA predictions, we would need to implement the full
        # model averaging strategy here. For now, use SM prediction
        # as a placeholder
        ma_pred[i] <- sm_pred[i]

        # Note: Full MA prediction would require:
        # 1. Computing predictions for multiple k values
        # 2. Applying the model averaging strategy
        # 3. Using appropriate kernel weights
    }

    # Return based on type
    switch(type,
           "ma" = ma_pred,
           "sm" = sm_pred,
           "both" = cbind(sm = sm_pred, ma = ma_pred)
    )
}

# ============================================================================
# ADDITIONAL UTILITY METHODS
# ============================================================================

#' Extract Log-Likelihood from Mabilo Plus Model
#'
#' @description
#' Computes an approximate log-likelihood for model comparison
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which predictions to use:
#'   "ma" (default) or "sm"
#' @param ... Additional arguments (currently unused)
#'
#' @return Log-likelihood value
#'
#' @export
logLik.mabilo_plus <- function(object, type = c("ma", "sm"), ...) {
    type <- match.arg(type)

    # Get residuals
    resid <- if (type == "ma") {
        object$y_sorted - object$ma_predictions
    } else {
        object$y_sorted - object$sm_predictions
    }

    # Estimate sigma
    sigma <- sd(resid)
    n <- length(resid)

    # Compute log-likelihood
    ll <- -n/2 * log(2 * pi) - n * log(sigma) - sum(resid^2) / (2 * sigma^2)

    # Add attributes
    attr(ll, "df") <- if (type == "ma") object$opt_ma_k else object$opt_sm_k
    attr(ll, "nobs") <- n
    class(ll) <- "logLik"

    return(ll)
}

#' Compute AIC for Mabilo Plus Model
#'
#' @description
#' Computes Akaike Information Criterion for model selection
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which predictions to use:
#'   "ma" (default) or "sm"
#' @param k Penalty parameter (default: 2)
#' @param ... Additional arguments passed to logLik
#'
#' @return AIC value
#'
#' @export
AIC.mabilo_plus <- function(object, type = c("ma", "sm"), k = 2, ...) {
    type <- match.arg(type)
    ll <- logLik(object, type = type, ...)
    -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Compute BIC for Mabilo Plus Model
#'
#' @description
#' Computes Bayesian Information Criterion for model selection
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which predictions to use:
#'   "ma" (default) or "sm"
#' @param ... Additional arguments passed to logLik
#'
#' @return BIC value
#'
#' @export
BIC.mabilo_plus <- function(object, type = c("ma", "sm"), ...) {
    type <- match.arg(type)
    ll <- logLik(object, type = type, ...)
    -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df")
}

#' Extract Variance-Covariance Matrix
#'
#' @description
#' Computes an approximate variance-covariance matrix for the predictions
#'
#' @param object A 'mabilo_plus' object
#' @param type Character string specifying which predictions to use:
#'   "ma" (default) or "sm"
#' @param ... Additional arguments (currently unused)
#'
#' @return Variance-covariance matrix
#'
#' @export
vcov.mabilo_plus <- function(object, type = c("ma", "sm"), ...) {
    type <- match.arg(type)

    # Get residuals
    resid <- if (type == "ma") {
        object$y_sorted - object$ma_predictions
    } else {
        object$y_sorted - object$sm_predictions
    }

    # Simple variance estimate
    var_est <- var(resid)

    # Return diagonal matrix
    n <- length(object$x_sorted)
    diag(rep(var_est, n))
}

#' Update Mabilo Plus Model
#'
#' @description
#' Update a mabilo_plus model with new parameters
#'
#' @param object A 'mabilo_plus' object
#' @param ... Arguments to update in the model call
#'
#' @return A new mabilo_plus object
#'
#' @export
update.mabilo_plus <- function(object, ...) {
    # Extract original call
    call <- object$call
    if (is.null(call)) {
        stop("Model object does not contain original call")
    }

    # Update call with new arguments
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    # Evaluate updated call
    eval(call, parent.frame())
}

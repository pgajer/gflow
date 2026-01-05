#' Model-Averaged Binary Locally-Weighted Logistic Smoothing (MABILOG)
#'
#' @description
#' Implements MABILOG algorithm for robust local logistic regression on binary data,
#' extending the MABILO framework by incorporating model averaging with local
#' logistic regression and Bayesian bootstrap for uncertainty quantification.
#' The algorithm uses symmetric k-hop neighborhoods and kernel-weighted averaging
#' for predictions.
#'
#' @param x Numeric vector of x coordinates (predictors).
#' @param y Numeric vector of binary response values (0 or 1).
#' @param y.true Optional numeric vector of true probabilities for error calculation.
#' @param max.iterations Maximum number of iterations for logistic regression convergence
#'   (default: 100).
#' @param ridge.lambda Ridge regression penalty parameter for stabilization
#'   (default: 0.002).
#' @param max.beta Maximum allowed absolute value for regression coefficients
#'   (default: 100.0).
#' @param tolerance Convergence tolerance for logistic regression
#'   (default: 1e-8).
#' @param k.min Minimum number of neighbors on each side (positive integer).
#'   Default is 5 percent of data points or 3, whichever is larger.
#' @param k.max Maximum number of neighbors on each side.
#'   Default is \code{(n-1)/2 - 1} where n is the number of data points.
#' @param n.bb Number of Bayesian bootstrap iterations (non-negative integer).
#'   Set to 0 to skip bootstrap (default: 100).
#' @param p Probability level for credible intervals (between 0 and 1, default: 0.95).
#' @param distance.kernel Integer specifying kernel type for distance weighting (1-10).
#'   Common choices: 1 = Uniform, 2 = Triangular, 3 = Epanechnikov,
#'   4 = Quartic, 5 = Tricube, 6 = Gaussian, 7 = Cosine (default: 1).
#' @param dist.normalization.factor Positive factor for distance normalization - must be greater than 1
#'   (default: 1.1).
#' @param verbose Logical; if TRUE, prints progress information (default: FALSE).
#'
#' @details
#' The MABILOG algorithm extends MABILO to binary classification problems by:
#' \itemize{
#'   \item Using local weighted logistic regression instead of linear regression
#'   \item Incorporating ridge regularization for numerical stability
#'   \item Implementing coefficient constraints to prevent extreme predictions
#'   \item Providing probability estimates rather than continuous predictions
#' }
#'
#' The algorithm automatically sorts input data by x values. For each point,
#' it uses k-hop neighborhoods (k points on each side when available) rather
#' than k-nearest neighbors, providing more symmetric neighborhoods. The optimal
#' k is selected by minimizing mean leave-one-out cross-validation error.
#'
#' The Bayesian bootstrap analysis (when n.bb > 0) provides uncertainty
#' quantification through credible intervals computed at the specified
#' probability level p.
#'
#' @return A list of class "mabilog" containing:
#'   \item{\code{k_values}}{Vector of tested k values}
#'   \item{\code{opt_k}}{Optimal k value selected by minimizing LOOCV error}
#'   \item{\code{opt_k_idx}}{Index of optimal k in k_values vector}
#'   \item{\code{k_mean_errors}}{Mean LOOCV errors for each k value}
#'   \item{\code{k_mean_true_errors}}{Mean true errors for each k (if y.true provided)}
#'   \item{\code{predictions}}{Model-averaged probability predictions using optimal k}
#'   \item{\code{errors}}{Leave-one-out cross-validation errors}
#'   \item{\code{k_predictions}}{List of predictions for each k value}
#'   \item{\code{bb_predictions}}{Bootstrap mean predictions (if n.bb > 0)}
#'   \item{\code{cri_L}}{Lower credible interval bounds (if n.bb > 0)}
#'   \item{\code{cri_U}}{Upper credible interval bounds (if n.bb > 0)}
#'   \item{\code{x_sorted}}{Input x values (sorted)}
#'   \item{\code{y_sorted}}{Input y values (sorted by x)}
#'   \item{\code{y_true_sorted}}{True probabilities (sorted by x, if provided)}
#'   \item{\code{k_min}}{Minimum k value used}
#'   \item{\code{k_max}}{Maximum k value used}
#'
#' @examples
#' # Generate binary data with smooth probability structure
#' set.seed(42)
#' x <- seq(0, 10, length.out = 200)
#' true_prob <- plogis(sin(x) - 0.5)
#' y <- rbinom(length(x), 1, true_prob)
#'
#' # Fit MABILOG model
#' fit <- mabilog(x, y, y.true = true_prob, k.min = 5, k.max = 20)
#'
#' # Plot results
#' plot(x, y, col = c("red", "blue")[y + 1], pch = 19, cex = 0.5)
#' lines(fit$x_sorted, fit$predictions, lwd = 2)
#' lines(x, true_prob, col = "green", lty = 2)
#' legend("topright", c("Data (y=0)", "Data (y=1)", "Fitted", "True"),
#'        col = c("red", "blue", "black", "green"),
#'        pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 2))
#'
#' # With bootstrap confidence intervals
#' \donttest{
#' fit_bb <- mabilog(x, y, k.min = 5, k.max = 20, n.bb = 100)
#' plot(x, y, col = c("red", "blue")[y + 1], pch = 19, cex = 0.5)
#' polygon(c(fit_bb$x_sorted, rev(fit_bb$x_sorted)),
#'         c(fit_bb$cri_L, rev(fit_bb$cri_U)),
#'         col = "gray80", border = NA)
#' lines(fit_bb$x_sorted, fit_bb$predictions, lwd = 2)
#' }
#'
#' @seealso
#' \code{\link{mabilo}} for continuous response regression,
#' \code{\link{predict.mabilog}} for predictions on new data,
#' \code{\link{plot.mabilog}} for diagnostic plots
#'
#' @export
mabilog <- function(x,
                    y,
                    y.true = NULL,
                    max.iterations = 100,
                    ridge.lambda = 0.002,
                    max.beta = 100.0,
                    tolerance = 1e-8,
                    k.min = max(3, as.integer(0.05 * length(x))),
                    k.max = as.integer((length(x)-1)/2) - 1, # this corresponds to 2k + 1 =  n - 1
                    n.bb = 100,
                    p = 0.95,
                    distance.kernel = 1,
                    dist.normalization.factor = 1.1,
                    verbose = FALSE) {

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

    if (is.null(k.max)) k.max <- as.integer((length(x)-1)/2) - 1
    if (!is.numeric(k.max) || k.max != round(k.max))
        stop("k.max must be an integer")
    if (k.max <= k.min)
        stop("k.max must be greater than k.min")
    if (k.max >= n)
        stop("k.max must be less than the number of points")

    if (k.max >= (n - 1) / 2) {
        print(sprintf("ERROR: k.max (%d) is too large for the data size (%d).\n", k.max, n))
        stop("k.max should be less than (n -1) / 2 to ensure multiple windows. This is equivalent to 2k + 1 < n.")
    }

    ## Check n.bb
    if (!is.numeric(n.bb) || n.bb != as.integer(n.bb)) stop("n.bb must be an integer")
    if (n.bb < 0) stop("n.bb must be non-negative")

    ## Check p
    if (!is.numeric(p)) stop("p must be numeric")
    if (length(p) != 1) stop("p must be a single value")
    if (p <= 0 || p >= 1) stop("p must be between 0 and 1")

    ## Check kernels
    if (!is.numeric(distance.kernel) || distance.kernel < 1 || distance.kernel > 10)
        stop("distance.kernel must be an integer between 1 and 10")

    ## Check other numeric parameters
    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 1)
        stop("dist.normalization.factor must be greater than 1")

    if (!is.numeric(max.iterations) || max.iterations <= 0 || max.iterations != round(max.iterations))
        stop("max.iterations must be a positive integer")

    if (!is.numeric(ridge.lambda) || ridge.lambda < 0)
        stop("ridge.lambda must be non-negative")

    if (!is.numeric(max.beta) || max.beta <= 0)
        stop("max.beta must be positive")

    if (!is.numeric(tolerance) || tolerance <= 0)
        stop("tolerance must be positive")

    ## Check verbose
    if (!is.logical(verbose))
        stop("verbose must be logical")

    ## Check if x is sorted and sort if necessary
    if (!identical(x, sort(x))) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        if (!is.null(y.true)) y.true <- y.true[ord]
    }

    ## Call the C++ implementation
    result <- .Call(S_mabilog,
                    as.double(x),
                    as.double(y),
                    as.double(y.true),
                    as.integer(max.iterations),
                    as.double(ridge.lambda),
                    as.double(max.beta),
                    as.double(tolerance),
                    as.integer(k.min),
                    as.integer(k.max),
                    as.integer(n.bb),
                    as.double(p),
                    as.integer(distance.kernel),
                    as.double(dist.normalization.factor),
                    as.logical(verbose))

    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- y.true
    result$k_min <- k.min
    result$k_max <- k.max

    class(result) <- "mabilog"

    return(result)
}

#' Complete Utility Functions for mabilog()
#'
#' This file contains all utility functions adapted for the mabilog() function,
#' including plot, print, summary, predict, fitted, residuals, coef, and other S3 methods.

# ============================================================================
# PLOT METHODS
# ============================================================================

#' Plot Method for Mabilog Objects
#'
#' @description
#' Generates diagnostic and visualization plots for mabilog objects. The function
#' supports different types of plots including fitted values with optional credible
#' intervals, error diagnostics, and residual analyses.
#'
#' @param x An object of class "mabilog"
#' @param type Character string specifying the type of plot to create (partial matching supported):
#'   * "fit": Plot fitted values against x values (e.g., "f" or "fi")
#'   * "diagnostic": Plot error diagnostics for different k values (e.g., "d" or "dia")
#'   * "residuals": Plot residuals against x values (e.g., "res")
#'   * "residuals.hist": Plot histogram of residuals (e.g., "rh" or "residuals.h")
#'   * "roc": ROC curve for binary classification (e.g., "ro")
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
#' @param pts.col Color for data points when y=0 and y=1 (default: c("red", "darkgreen"))
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
#' @return Invisibly returns NULL. The function is called for its side effect of producing plots.
#'
#' @export
plot.mabilog <- function(x, type = "fit",
                         title = "", xlab = "", ylab = "",
                         with.y.true = TRUE,
                         with.pts = FALSE,
                         with.CrI = TRUE,
                         true.col = "red",
                         ma.col = "blue",
                         pts.col = c("red", "darkgreen"),
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

    if (!inherits(x, "mabilog")) {
        stop("Input must be a 'mabilog' object")
    }

    # Define valid plot types
    valid_types <- c("fit", "diagnostic", "residuals", "residuals.hist", "roc")

    # Find partial matches for the provided type
    type_match <- grep(paste0("^", type), valid_types, value = TRUE)

    # Handle matching results
    if (length(type_match) == 0) {
        stop("Invalid plot type. Use 'fit', 'diagnostic', 'residuals', 'residuals.hist', or 'roc'")
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

    switch(type,
           "fit" = mabilog.plot.fit(x, title, xlab, ylab, with.y.true,
                                   with.pts, with.CrI, true.col, ma.col, pts.col,
                                   with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                                   CrI.as.polygon, CrI.polygon.col, CrI.line.col, CrI.line.lty,
                                   ylim, legend.cex, k, with.legend, legend.pos, ...),
           "diagnostic" = mabilog.plot.diagnostic(x, title, xlab, ylab, ma.col, true.col,
                                                 legend.cex, legend.pos, ...),
           "residuals" = mabilog.plot.residuals(x, title, xlab, ylab, ...),
           "residuals.hist" = mabilog.plot.residuals.hist(x, title, xlab, ylab, ...),
           "roc" = mabilog.plot.roc(x, title, xlab, ylab, ma.col, legend.cex, ...),
           stop("Invalid plot type.")
    )

    invisible(NULL)
}

# Helper function for fit plots
mabilog.plot.fit <- function(res, title, xlab, ylab, with.y.true,
                            with.pts, with.CrI, true.col, ma.col, pts.col,
                            with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                            CrI.as.polygon, CrI.polygon.col, CrI.line.col, CrI.line.lty,
                            ylim, legend.cex, k, with.legend, legend.pos, ...) {

    # Select appropriate predictions based on k
    if (!is.null(k)) {
        k_idx <- which(res$k_values == k)
        predictions <- res$k_predictions[[k_idx]]
        pred_label <- sprintf("Predictions (k=%d)", k)
    } else {
        predictions <- res$predictions
        pred_label <- "Predictions"
    }

    # Check for bootstrap results
    has_cri <- !is.null(res$bb_predictions) && !is.null(res$cri_L) && !is.null(res$cri_U)

    # Calculate y-limits if not provided
    if (is.null(ylim)) {
        ylim <- c(-0.1, 1.1)  # Standard range for probabilities
    }

    # Initialize plot
    plot(res$x_sorted, res$y_sorted, type = "n",
         las = 1, ylim = ylim,
         xlab = if(xlab == "") "x" else xlab,
         ylab = if(ylab == "") "Probability" else ylab,
         main = if(title == "") "Mabilog Fit" else title, ...)

    # Add credible interval first if available
    if (with.CrI && has_cri) {
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
    lines(res$x_sorted, predictions, col = ma.col, lwd = 2, ...)

    # Add original data points if requested
    if (with.pts) {
        # Color points by y value
        point_colors <- ifelse(res$y_sorted == 0, pts.col[1], pts.col[2])
        points(res$x_sorted, res$y_sorted, col = point_colors, pch = 19, cex = 0.5, ...)
    }

    # Add prediction points if requested
    if (with.predictions.pts) {
        points(res$x_sorted, predictions,
              col = predictions.pts.col,
              pch = predictions.pts.pch, cex = 0.5)
    }

    # Add true values if available and requested
    if (with.y.true && !is.null(res$y_true_sorted) && length(res$y_true_sorted) == length(res$x_sorted)) {
        lines(res$x_sorted, res$y_true_sorted, col = true.col, lty = 2, lwd = 2)
    }

    # Add legend
    if (with.legend) {
        legend_items <- c()
        legend_cols <- c()
        legend_ltys <- c()
        legend_pchs <- c()

        if (with.pts) {
            legend_items <- c(legend_items, "y = 0", "y = 1")
            legend_cols <- c(legend_cols, pts.col[1], pts.col[2])
            legend_ltys <- c(legend_ltys, NA, NA)
            legend_pchs <- c(legend_pchs, 19, 19)
        }

        legend_items <- c(legend_items, pred_label)
        legend_cols <- c(legend_cols, ma.col)
        legend_ltys <- c(legend_ltys, 1)
        legend_pchs <- c(legend_pchs, NA)

        if (has_cri && with.CrI) {
            ci_label <- sprintf("%.0f%% Credible Interval",
                               100 * (as.numeric(res$call$p) %||% 0.95))
            legend_items <- c(legend_items, ci_label)
            legend_cols <- c(legend_cols, if(CrI.as.polygon) CrI.polygon.col else CrI.line.col)
            legend_ltys <- c(legend_ltys, if(CrI.as.polygon) 1 else CrI.line.lty)
            legend_pchs <- c(legend_pchs, NA)
        }

        if (with.y.true && !is.null(res$y_true_sorted)) {
            legend_items <- c(legend_items, "True Probabilities")
            legend_cols <- c(legend_cols, true.col)
            legend_ltys <- c(legend_ltys, 2)
            legend_pchs <- c(legend_pchs, NA)
        }

        legend(legend.pos, legend = legend_items,
               col = legend_cols, lty = legend_ltys, pch = legend_pchs,
               bg = "white", inset = 0.05, cex = legend.cex)
    }
}

# Helper function for diagnostic plots
mabilog.plot.diagnostic <- function(res, title, xlab, ylab, ma.col, true.col,
                                   legend.cex, legend.pos = "topright", ...) {
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
mabilog.plot.residuals <- function(res, title, xlab, ylab, ...) {
    # For binary data, use deviance residuals
    residuals <- sign(res$y_sorted - res$predictions) *
                 sqrt(-2 * (res$y_sorted * log(pmax(res$predictions, 1e-10)) +
                           (1 - res$y_sorted) * log(pmax(1 - res$predictions, 1e-10))))

    plot(res$x_sorted, residuals,
         xlab = if(xlab == "") "x" else xlab,
         ylab = if(ylab == "") "Deviance Residuals" else ylab,
         main = if(title == "") "Deviance Residuals Plot" else title,
         las = 1, ...)

    abline(h = 0, lty = 2)

    # Add loess smooth
    lo <- loess(residuals ~ res$x_sorted, span = 0.75)
    lines(res$x_sorted, predict(lo), col = "red", lwd = 2)
}

# Helper function for residual histograms
mabilog.plot.residuals.hist <- function(res, title, xlab, ylab, ...) {
    # Deviance residuals
    residuals <- sign(res$y_sorted - res$predictions) *
                 sqrt(-2 * (res$y_sorted * log(pmax(res$predictions, 1e-10)) +
                           (1 - res$y_sorted) * log(pmax(1 - res$predictions, 1e-10))))

    hist(residuals,
         xlab = if(xlab == "") "Deviance Residuals" else xlab,
         ylab = if(ylab == "") "Frequency" else ylab,
         main = if(title == "") "Histogram of Deviance Residuals" else title,
         probability = TRUE,
         breaks = "FD",
         ...)

    # Add normal density for comparison
    x_seq <- seq(min(residuals), max(residuals), length.out = 100)
    lines(x_seq, dnorm(x_seq, mean = mean(residuals), sd = sd(residuals)),
          col = "blue", lwd = 2)

    abline(v = 0, lty = 2)
}

# Helper function for ROC plots
mabilog.plot.roc <- function(res, title, xlab, ylab, ma.col, legend.cex, ...) {
    if (is.null(res$y_sorted) || is.null(res$predictions)) {
        stop("Cannot plot ROC curve: missing y values or predictions")
    }

    # Calculate ROC curve
    thresholds <- seq(0, 1, by = 0.01)
    tpr <- numeric(length(thresholds))
    fpr <- numeric(length(thresholds))

    for (i in seq_along(thresholds)) {
        pred_class <- as.numeric(res$predictions >= thresholds[i])
        tp <- sum(pred_class == 1 & res$y_sorted == 1)
        fp <- sum(pred_class == 1 & res$y_sorted == 0)
        tn <- sum(pred_class == 0 & res$y_sorted == 0)
        fn <- sum(pred_class == 0 & res$y_sorted == 1)

        tpr[i] <- tp / (tp + fn)
        fpr[i] <- fp / (fp + tn)
    }

    # Calculate AUC
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)

    # Plot
    plot(fpr, tpr, type = "l",
         xlab = if(xlab == "") "False Positive Rate" else xlab,
         ylab = if(ylab == "") "True Positive Rate" else ylab,
         main = if(title == "") sprintf("ROC Curve (AUC = %.3f)", auc) else title,
         col = ma.col, lwd = 2, las = 1, ...)

    # Add diagonal reference line
    abline(a = 0, b = 1, lty = 2, col = "gray")

    # Add AUC to plot
    text(0.7, 0.3, sprintf("AUC = %.3f", auc), cex = legend.cex)
}

# ============================================================================
# SUMMARY METHODS
# ============================================================================

#' Compute Summary Statistics for Mabilog Objects
#'
#' @description
#' Computes comprehensive summary statistics for mabilog fits including model parameters,
#' fit statistics, error analysis, and diagnostic information for logistic regression.
#' Also includes Bayesian bootstrap statistics when available.
#'
#' @param object A 'mabilog' object
#' @param quantiles Numeric vector of probabilities for quantile computations
#' @param ... Additional arguments (currently unused)
#'
#' @return A 'summary.mabilog' object containing:
#' \itemize{
#'   \item \code{model_info}: Basic information about the model fit
#'   \item \code{fit_stats}: Prediction accuracy metrics
#'   \item \code{k_error_stats}: Error statistics for different k values
#'   \item \code{residual_stats}: Residual analysis results
#'   \item \code{classification_stats}: Classification performance metrics
#'   \item \code{true_error_stats}: Statistics comparing to true values (if available)
#'   \item \code{bootstrap_stats}: Bayesian bootstrap statistics (if available)
#' }
#'
#' @export
summary.mabilog <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    if (!inherits(object, "mabilog")) {
        stop("Input must be a 'mabilog' object")
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("quantiles must be numeric values between 0 and 1")
    }

    # Check for missing values
    if (any(is.na(object$y_sorted)) || any(is.na(object$predictions))) {
        warning("Missing values detected in fit results")
    }

    # Calculate deviance residuals
    residuals <- sign(object$y_sorted - object$predictions) *
                 sqrt(-2 * (object$y_sorted * log(pmax(object$predictions, 1e-10)) +
                           (1 - object$y_sorted) * log(pmax(1 - object$predictions, 1e-10))))

    # Basic model information
    model_info <- list(
        n_observations = length(object$x_sorted),
        n_positive = sum(object$y_sorted == 1),
        n_negative = sum(object$y_sorted == 0),
        optimal_k = object$opt_k,
        k_range = range(object$k_values),
        min_x = min(object$x_sorted),
        max_x = max(object$x_sorted),
        range_x = diff(range(object$x_sorted))
    )

    # Fit statistics
    # Log-likelihood
    ll <- sum(object$y_sorted * log(pmax(object$predictions, 1e-10)) +
              (1 - object$y_sorted) * log(pmax(1 - object$predictions, 1e-10)))

    # Deviance
    deviance <- -2 * ll

    # Brier score
    brier_score <- mean((object$predictions - object$y_sorted)^2)

    fit_stats <- list(
        log_likelihood = ll,
        deviance = deviance,
        brier_score = brier_score,
        mean_prediction = mean(object$predictions),
        sd_prediction = sd(object$predictions)
    )

    # Classification statistics (using 0.5 threshold)
    pred_class <- as.numeric(object$predictions >= 0.5)
    confusion <- table(Predicted = pred_class, Actual = object$y_sorted)

    tp <- sum(pred_class == 1 & object$y_sorted == 1)
    tn <- sum(pred_class == 0 & object$y_sorted == 0)
    fp <- sum(pred_class == 1 & object$y_sorted == 0)
    fn <- sum(pred_class == 0 & object$y_sorted == 1)

    classification_stats <- list(
        confusion_matrix = confusion,
        accuracy = (tp + tn) / length(object$y_sorted),
        sensitivity = tp / (tp + fn),
        specificity = tn / (tn + fp),
        precision = tp / (tp + fp),
        f1_score = 2 * tp / (2 * tp + fp + fn)
    )

    # True error statistics if available
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted)) {
        # Mean squared error between predicted and true probabilities
        true_mse <- mean((object$predictions - object$y_true_sorted)^2)
        true_rmse <- sqrt(true_mse)
        true_mae <- mean(abs(object$predictions - object$y_true_sorted))

        # KL divergence
        kl_div <- mean(object$y_true_sorted * log(pmax(object$y_true_sorted /
                                                       pmax(object$predictions, 1e-10), 1e-10)) +
                      (1 - object$y_true_sorted) * log(pmax((1 - object$y_true_sorted) /
                                                            pmax(1 - object$predictions, 1e-10), 1e-10)))

        true_error_stats <- list(
            true_mse = true_mse,
            true_rmse = true_rmse,
            true_mae = true_mae,
            kl_divergence = kl_div
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
        classification_stats = classification_stats,
        true_error_stats = true_error_stats,
        bootstrap_stats = bootstrap_stats
    )

    class(result) <- "summary.mabilog"
    return(result)
}

#' Print Summary Statistics for Mabilog Fits
#'
#' @description
#' Formats and displays comprehensive summary statistics for mabilog fits.
#' Includes model parameters, fit statistics, error analysis, diagnostic information,
#' and Bayesian bootstrap results when available.
#'
#' @param x A 'summary.mabilog' object from summary.mabilog
#' @param digits Number of significant digits for numerical output (default: 4)
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.summary.mabilog <- function(x, digits = 4, ...) {
    # Helper function for separator lines
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
    section_header("MABILOG (Model-Averaged Binary Locally-Weighted Logistic Smoothing) SUMMARY")

    # Model Information
    cat("\nModel Information:\n")
    cat(sprintf("Number of observations:         %d (Positive: %d, Negative: %d)\n",
                x$model_info$n_observations,
                x$model_info$n_positive,
                x$model_info$n_negative))
    cat(sprintf("Optimal k:                      %d\n", x$model_info$optimal_k))
    cat(sprintf("k range:                        [%d, %d]\n",
                x$model_info$k_range[1], x$model_info$k_range[2]))
    cat(sprintf("X range:                        [%.3f, %.3f]\n",
                x$model_info$min_x, x$model_info$max_x))
    cat(sprintf("X span:                         %.3f\n", x$model_info$range_x))

    # Fit Statistics
    section_header("FIT STATISTICS")
    cat(sprintf("Log-likelihood:                 %.4f\n", x$fit_stats$log_likelihood))
    cat(sprintf("Deviance:                       %.4f\n", x$fit_stats$deviance))
    cat(sprintf("Brier Score:                    %.4f\n", x$fit_stats$brier_score))
    cat(sprintf("Mean Prediction:                %.4f\n", x$fit_stats$mean_prediction))
    cat(sprintf("SD of Predictions:              %.4f\n", x$fit_stats$sd_prediction))

    # Classification Statistics
    section_header("CLASSIFICATION PERFORMANCE (Threshold = 0.5)")
    cat("\nConfusion Matrix:\n")
    print(x$classification_stats$confusion_matrix)
    cat(sprintf("\nAccuracy:                       %.4f\n", x$classification_stats$accuracy))
    cat(sprintf("Sensitivity (Recall):           %.4f\n", x$classification_stats$sensitivity))
    cat(sprintf("Specificity:                    %.4f\n", x$classification_stats$specificity))
    cat(sprintf("Precision:                      %.4f\n", x$classification_stats$precision))
    cat(sprintf("F1 Score:                       %.4f\n", x$classification_stats$f1_score))

    # True Error Statistics (if available)
    if (!is.null(x$true_error_stats)) {
        section_header("TRUE PROBABILITY ERROR STATISTICS")
        cat(sprintf("True MSE:                       %.4f\n", x$true_error_stats$true_mse))
        cat(sprintf("True RMSE:                      %.4f\n", x$true_error_stats$true_rmse))
        cat(sprintf("True MAE:                       %.4f\n", x$true_error_stats$true_mae))
        cat(sprintf("KL Divergence:                  %.4f\n", x$true_error_stats$kl_divergence))
    }

    # k Error Statistics
    section_header("k-VALUE ERROR ANALYSIS")
    cat(sprintf("Minimum Error:                  %.4f\n", x$k_error_stats$min_error))
    cat(sprintf("Optimal k Error:                %.4f\n", x$k_error_stats$optimal_k_error))

    # Residual Statistics
    section_header("DEVIANCE RESIDUAL STATISTICS")
    cat(sprintf("Mean:                           %.4f\n", x$residual_stats$mean))
    cat(sprintf("Standard Deviation:             %.4f\n", x$residual_stats$sd))
    cat("Quantiles:\n")
    print(x$residual_stats$quantiles, digits = digits)

    # Bootstrap Statistics (if available)
    if (!is.null(x$bootstrap_stats)) {
        section_header("BAYESIAN BOOTSTRAP STATISTICS")
        cat(sprintf("Mean CI Width:                  %.4f\n", x$bootstrap_stats$mean_ci_width))
        cat("Note: Full bootstrap predictions and intervals available in object\n")
    }

    hr()
    invisible(x)
}

# ============================================================================
# PRINT METHOD
# ============================================================================

#' Print Method for Mabilog Objects
#'
#' @description
#' Prints a concise summary of a mabilog object
#'
#' @param x A 'mabilog' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.mabilog <- function(x, ...) {
    cat("Mabilog (Model-Averaged Binary Locally-Weighted Logistic Smoothing) Fit\n")
    cat("======================================================================\n")
    cat(sprintf("Number of observations: %d\n", length(x$x_sorted)))
    cat(sprintf("Response: %d zeros, %d ones\n",
                sum(x$y_sorted == 0), sum(x$y_sorted == 1)))
    cat(sprintf("k range: [%d, %d]\n", x$k_min, x$k_max))
    cat(sprintf("Optimal k: %d\n", x$opt_k))
    if (!is.null(x$bb_predictions)) {
        cat(sprintf("Bootstrap iterations: %d\n",
                   as.numeric(x$call$n.bb) %||% 100))
    }
    cat("\nCall summary(object) for detailed statistics\n")
    invisible(x)
}

# ============================================================================
# EXTRACTION METHODS
# ============================================================================

#' Extract Fitted Values from Mabilog Model
#'
#' @description
#' Extracts fitted probability values from a mabilog object
#'
#' @param object A 'mabilog' object
#' @param type Character string specifying which values to return:
#'   "response" for probabilities (default), "link" for logit scale
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector of fitted values
#'
#' @export
fitted.mabilog <- function(object, type = c("response", "link"), ...) {
    type <- match.arg(type)

    if (type == "response") {
        return(object$predictions)
    } else {
        # Convert to logit scale
        p <- pmax(pmin(object$predictions, 1 - 1e-10), 1e-10)
        return(log(p / (1 - p)))
    }
}

#' Extract Residuals from Mabilog Model
#'
#' @description
#' Extracts residuals from a mabilog object
#'
#' @param object A 'mabilog' object
#' @param type Character string specifying residual type:
#'   "deviance" (default), "pearson", "response", or "working"
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector of residuals
#'
#' @export
residuals.mabilog <- function(object, type = c("deviance", "pearson", "response", "working"), ...) {
    type <- match.arg(type)

    y <- object$y_sorted
    mu <- object$predictions

    # Protect against extreme values
    mu <- pmax(pmin(mu, 1 - 1e-10), 1e-10)

    switch(type,
        "deviance" = {
            sign(y - mu) * sqrt(-2 * (y * log(mu) + (1 - y) * log(1 - mu)))
        },
        "pearson" = {
            (y - mu) / sqrt(mu * (1 - mu))
        },
        "response" = {
            y - mu
        },
        "working" = {
            # Working residuals for logistic regression
            (y - mu) / (mu * (1 - mu))
        }
    )
}

#' Extract Model Coefficients from Mabilog Model
#'
#' @description
#' Extracts the optimal k value and model information from a mabilog object
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Named vector containing optimal k value and other parameters
#'
#' @export
coef.mabilog <- function(object, ...) {
    c(opt_k = object$opt_k,
      k_min = object$k_min,
      k_max = object$k_max,
      n_obs = length(object$x_sorted))
}

# ============================================================================
# PREDICTION METHOD
# ============================================================================

#' Predict Method for Mabilog Objects
#'
#' @description
#' Predict probability values for new data using a fitted mabilog model
#'
#' @param object A 'mabilog' object
#' @param newdata Numeric vector of new x values for prediction
#' @param type Character string specifying the type of prediction:
#'   "response" for probabilities (default), "link" for logit scale
#' @param k Optional specific k value to use for predictions. If NULL,
#'   uses optimal k value
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector of predictions
#'
#' @details
#' This function performs local weighted logistic regression predictions for new data points.
#' For each new x value, it finds the k nearest neighbors from the training data
#' and performs a weighted logistic regression.
#'
#' Note: This is a simplified implementation. Full local logistic regression
#' would require iterative fitting at each prediction point.
#'
#' @export
predict.mabilog <- function(object, newdata, type = c("response", "link"),
                           k = NULL, ...) {
    type <- match.arg(type)

    if (!is.numeric(newdata)) {
        stop("newdata must be numeric")
    }

    n_new <- length(newdata)
    n_train <- length(object$x_sorted)

    # Determine which k value to use
    if (!is.null(k)) {
        if (!is.numeric(k) || length(k) != 1) {
            stop("k must be a single numeric value")
        }
        if (k < object$k_min || k > object$k_max) {
            stop(sprintf("k must be between %d and %d", object$k_min, object$k_max))
        }
    } else {
        k <- object$opt_k
    }

    # Initialize prediction array
    predictions <- numeric(n_new)

    # For each new point
    for (i in 1:n_new) {
        x_new <- newdata[i]

        # Calculate distances to all training points
        distances <- abs(object$x_sorted - x_new)

        # Get indices of k nearest neighbors
        neighbors <- order(distances)[1:k]

        # Get neighbor data
        x_neighbors <- object$x_sorted[neighbors]
        y_neighbors <- object$y_sorted[neighbors]

        # Simple approach: weighted average of neighbor y values
        # (A full implementation would fit local logistic regression)
        weights <- exp(-distances[neighbors]^2 / (2 * median(distances[neighbors])^2))
        weights <- weights / sum(weights)

        # Weighted average (simple approximation)
        predictions[i] <- sum(weights * y_neighbors)

        # Ensure valid probability
        predictions[i] <- pmax(pmin(predictions[i], 1 - 1e-10), 1e-10)
    }

    # Return based on type
    if (type == "response") {
        return(predictions)
    } else {
        # Convert to logit scale
        return(log(predictions / (1 - predictions)))
    }
}

# ============================================================================
# ADDITIONAL UTILITY METHODS
# ============================================================================

#' Extract Log-Likelihood from Mabilog Model
#'
#' @description
#' Computes the log-likelihood for the fitted model
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Log-likelihood value with attributes
#'
#' @export
logLik.mabilog <- function(object, ...) {
    # Ensure valid probabilities
    p <- pmax(pmin(object$predictions, 1 - 1e-10), 1e-10)

    # Compute log-likelihood
    ll <- sum(object$y_sorted * log(p) + (1 - object$y_sorted) * log(1 - p))

    # Add attributes
    attr(ll, "df") <- object$opt_k  # Approximate degrees of freedom
    attr(ll, "nobs") <- length(object$y_sorted)
    class(ll) <- "logLik"

    return(ll)
}

#' Compute AIC for Mabilog Model
#'
#' @description
#' Computes Akaike Information Criterion for model selection
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments passed to logLik
#' @param k Penalty parameter (default: 2)
#'
#' @return AIC value
#'
#' @export
AIC.mabilog <- function(object, ..., k = 2) {
    ll <- logLik(object, ...)
    -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Compute BIC for Mabilog Model
#'
#' @description
#' Computes Bayesian Information Criterion for model selection
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments passed to logLik
#'
#' @return BIC value
#'
#' @export
BIC.mabilog <- function(object, ...) {
    ll <- logLik(object, ...)
    -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df")
}

#' Compute Deviance for Mabilog Model
#'
#' @description
#' Computes the deviance of the fitted model
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Deviance value
#'
#' @export
deviance.mabilog <- function(object, ...) {
    -2 * as.numeric(logLik(object, ...))
}

#' Extract Variance-Covariance Matrix
#'
#' @description
#' Computes an approximate variance-covariance matrix for the predictions
#' Note: This is a simplified approximation for binomial variance
#'
#' @param object A 'mabilog' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Diagonal variance-covariance matrix
#'
#' @export
vcov.mabilog <- function(object, ...) {
    # Binomial variance approximation
    p <- object$predictions
    var_est <- p * (1 - p)

    # Return diagonal matrix
    n <- length(object$x_sorted)
    diag(var_est)
}

#' Update Mabilog Model
#'
#' @description
#' Update a mabilog model with new parameters
#'
#' @param object A 'mabilog' object
#' @param ... Arguments to update in the model call
#'
#' @return A new mabilog object
#'
#' @export
update.mabilog <- function(object, ...) {
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


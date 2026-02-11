#' Grid-based Model Averaged Bandwidth Logistic Regression
#'
#' @description
#' Performs model-averaged bandwidth logistic regression using local polynomial fitting.
#' The function implements a flexible approach to binary regression by fitting local
#' models at points of a uniform grid spanning the range of predictor values. Multiple
#' local models are combined to produce robust predictions. It supports both linear and
#' quadratic local models, automatic bandwidth selection via cross-validation, and
#' various kernel types for weight calculation.
#'
#' @param x Numeric vector of predictor variables. Must not contain missing, infinite,
#'   or non-numeric values.
#' @param y Binary numeric vector (0 or 1) of response variables. Must be the same
#'   length as \code{x} and contain only 0s and 1s.
#' @param grid.size Integer; number of points in the uniform grid where local models
#'   are centered. Must be at least 2. Larger values provide smoother predictions but
#'   increase computation time. Default is 200.
#' @param fit.quadratic Logical; whether to include quadratic terms in the local models.
#'   If \code{TRUE}, local quadratic logistic regression is used; if \code{FALSE},
#'   local linear logistic regression is used. Default is \code{FALSE}.
#' @param pilot.bandwidth Numeric; bandwidth for local fitting. If less than or equal to 0,
#'   bandwidth is automatically selected using cross-validation. The bandwidth controls
#'   the size of the local neighborhood used for fitting. Default is -1 (automatic selection).
#' @param kernel Integer; kernel type for weight calculation:
#'   \describe{
#'     \item{1}{Epanechnikov kernel}
#'     \item{2}{Triangular kernel}
#'     \item{4}{Laplace kernel}
#'     \item{5}{Normal (Gaussian) kernel}
#'     \item{6}{Biweight kernel}
#'     \item{7}{Tricube kernel (default)}
#'   }
#' @param min.points Integer or \code{NULL}; minimum number of points required for local fitting.
#'   If \code{NULL}, automatically set to 3 for linear models or 4 for quadratic models.
#'   Must be at least 3 for linear models and 4 for quadratic models to ensure identifiability.
#'   Default is \code{NULL}.
#' @param cv.folds Integer; number of cross-validation folds for bandwidth selection.
#'   Must be at least 3 and cannot exceed the number of observations. Higher values
#'   provide more stable bandwidth selection but increase computation time. Default is 5.
#' @param n.bws Integer; number of bandwidths to try in automatic selection. Must be
#'   at least 2. More bandwidths provide finer resolution but increase computation time.
#'   Default is 50.
#' @param min.bw.factor Numeric; minimum bandwidth factor relative to data range. Must be
#'   positive and less than 1. Controls the smallest bandwidth to consider as a fraction
#'   of the x-range. Default is 0.05.
#' @param max.bw.factor Numeric; maximum bandwidth factor relative to data range. Must be
#'   greater than \code{min.bw.factor}. Controls the largest bandwidth to consider as a
#'   fraction of the x-range. Default is 0.9.
#' @param max.iterations Integer; maximum number of iterations for local fitting algorithm.
#'   Must be positive. Default is 100.
#' @param ridge.lambda Numeric; ridge parameter for numerical stability in local fitting.
#'   Must be positive. Larger values provide more stability but may introduce bias.
#'   Default is 1e-6.
#' @param tolerance Numeric; convergence tolerance for local fitting algorithm. Must be
#'   positive. Smaller values provide more precise convergence but may increase iterations.
#'   Default is 1e-8.
#' @param with.bw.predictions Logical; whether to return predictions for all bandwidth
#'   values tried during cross-validation. If \code{TRUE}, allows examination of how
#'   predictions vary with bandwidth. Default is \code{TRUE}.
#'
#' @return A list of class \code{"magelog"} containing:
#'   \describe{
#'     \item{\code{x.grid}}{Numeric vector of uniform grid points where local models are centered}
#'     \item{\code{predictions}}{Numeric vector of predicted probabilities at original x points}
#'     \item{\code{bw.grid.predictions}}{Matrix of predictions at grid points for each bandwidth
#'       (only if \code{with.bw.predictions = TRUE}). Rows correspond to grid points, columns
#'       to different bandwidths}
#'     \item{\code{mean.brier.errors}}{Numeric vector of cross-validation Brier scores for each bandwidth}
#'     \item{\code{opt.brier.bw.idx}}{Integer index of optimal bandwidth minimizing Brier score}
#'     \item{\code{bws}}{Numeric vector of bandwidth values tried}
#'     \item{\code{fit.info}}{List containing fitting parameters used}
#'     \item{\code{x}}{Original predictor values (sorted if necessary)}
#'     \item{\code{y}}{Original response values (sorted to match x)}
#'   }
#'
#' @details
#' The function fits local logistic regression models centered at points of a uniform
#' grid spanning the range of x values. Local models are fit using kernel-weighted
#' maximum likelihood estimation with optional ridge regularization for numerical stability.
#'
#' The bandwidth parameter controls the size of the local neighborhood used for each
#' local fit. When \code{pilot.bandwidth <= 0}, the function automatically selects an
#' optimal bandwidth using k-fold cross-validation to minimize the Brier score.
#'
#' Local models can be either linear or quadratic (controlled by \code{fit.quadratic}).
#' For numerical stability, each local fit requires a minimum number of points: at least
#' 3 for linear models and 4 for quadratic models.
#'
#' Predictions at the original x points are obtained by linear interpolation from the
#' grid-based predictions. This grid-based approach provides computational efficiency
#' and naturally smooth prediction curves.
#'
#' The function uses compiled C++ code for computational efficiency, particularly
#' important when performing cross-validation over multiple bandwidth values.
#'
#' @examples
#' \dontrun{
#' # Generate example data with logistic relationship
#' set.seed(123)
#' n <- 200
#' x <- seq(0, 1, length.out = n)
#' # True probability function
#' p_true <- 1/(1 + exp(-(x - 0.5)*10))
#' y <- rbinom(n, 1, p_true)
#'
#' # Fit model with automatic bandwidth selection
#' fit <- magelog(x, y, grid.size = 100, fit.quadratic = FALSE, cv.folds = 5)
#'
#' # Plot results
#' plot(x, y, pch = 19, col = adjustcolor("black", 0.5),
#'      xlab = "x", ylab = "Probability",
#'      main = "Local Logistic Regression")
#' # Add true probability curve
#' lines(x, p_true, col = "gray", lwd = 2, lty = 2)
#' # Add fitted curve
#' lines(fit$x.grid, fit$bw.grid.predictions[, fit$opt.brier.bw.idx],
#'       col = "red", lwd = 2)
#' legend("topleft", c("True probability", "Fitted curve", "Data"),
#'        col = c("gray", "red", "black"), lty = c(2, 1, NA),
#'        pch = c(NA, NA, 19), lwd = c(2, 2, NA))
#'
#' # Examine bandwidth selection
#' plot(fit$bws, fit$mean.brier.errors, type = "b",
#'      xlab = "Bandwidth", ylab = "Cross-validation Brier Score",
#'      main = "Bandwidth Selection")
#' abline(v = fit$bws[fit$opt.brier.bw.idx], col = "red", lty = 2)
#' }
#'
#' @seealso
#' \code{\link{predict.magelog}} for making predictions on new data,
#'
#' @references
#' Loader, C. (1999). Local Regression and Likelihood. Springer-Verlag.
#'
#' Fan, J. and Gijbels, I. (1996). Local Polynomial Modelling and Its Applications.
#' Chapman & Hall.
#'
#' @importFrom stats approx
#' @export
magelog <- function(x,
                    y,
                    grid.size = 200,
                    fit.quadratic = FALSE,
                    pilot.bandwidth = -1,
                    kernel = 7L,
                    min.points = NULL,
                    cv.folds = 5L,
                    n.bws = 50L,
                    min.bw.factor = 0.05,
                    max.bw.factor = 0.9,
                    max.iterations = 100L,
                    ridge.lambda = 1e-6,
                    tolerance = 1e-8,
                    with.bw.predictions = TRUE) {

    ## ====================================================================
    ## Input validation
    ## ====================================================================

    ## Basic type checks
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    if (!is.numeric(y)) stop("'y' must be a numeric vector")

    ## Convert to numeric to handle integer input
    x <- as.numeric(x)
    y <- as.numeric(y)

    ## Check for binary values
    if (!all(y %in% c(0, 1))) {
        stop("'y' must contain only binary values (0 and 1)")
    }

    ## Length and missing value checks
    n <- length(x)
    if (n != length(y)) stop("'x' and 'y' must have the same length")
    if (n == 0) stop("Input vectors cannot be empty")
    if (anyNA(x) || any(!is.finite(x))) {
        stop("'x' contains NA, NaN, or infinite values")
    }
    if (anyNA(y) || any(!is.finite(y))) {
        stop("'y' contains NA, NaN, or infinite values")
    }

    ## Parameter validation
    grid.size <- as.integer(grid.size)
    if (grid.size < 2) stop("'grid.size' must be at least 2")

    if (!is.logical(fit.quadratic) || length(fit.quadratic) != 1) {
        stop("'fit.quadratic' must be a single logical value")
    }

    if (!is.logical(with.bw.predictions) || length(with.bw.predictions) != 1) {
        stop("'with.bw.predictions' must be a single logical value")
    }

    ## Set minimum points based on model type
    required.min.points <- if (fit.quadratic) 4L else 3L
    if (is.null(min.points)) {
        min.points <- required.min.points
    } else {
        min.points <- as.integer(min.points)
        if (min.points < required.min.points) {
            stop(sprintf("'min.points' must be at least %d for %s model",
                        required.min.points,
                        if(fit.quadratic) "quadratic" else "linear"))
        }
    }

    ## Check dataset size
    if (n < min.points) {
        stop(sprintf("Dataset must contain at least %d observations", min.points))
    }

    ## Numeric parameter checks
    if (!is.numeric(pilot.bandwidth) || length(pilot.bandwidth) != 1) {
        stop("'pilot.bandwidth' must be a single numeric value")
    }

    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1) {
        stop("'min.bw.factor' must be a single numeric value")
    }
    if (min.bw.factor <= 0) stop("'min.bw.factor' must be positive")
    if (min.bw.factor >= 1) stop("'min.bw.factor' must be less than 1")

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1) {
        stop("'max.bw.factor' must be a single numeric value")
    }
    if (max.bw.factor <= min.bw.factor) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor'")
    }

    if (!is.numeric(ridge.lambda) || length(ridge.lambda) != 1) {
        stop("'ridge.lambda' must be a single numeric value")
    }
    if (ridge.lambda <= 0) stop("'ridge.lambda' must be positive")

    if (!is.numeric(tolerance) || length(tolerance) != 1) {
        stop("'tolerance' must be a single numeric value")
    }
    if (tolerance <= 0) stop("'tolerance' must be positive")

    ## Integer parameter checks
    kernel <- as.integer(kernel)
    if (!kernel %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel' must be one of: 1 (Epanechnikov), 2 (Triangular), ",
             "4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    cv.folds <- as.integer(cv.folds)
    if (cv.folds < 3) stop("'cv.folds' must be at least 3")
    if (cv.folds > n) stop("'cv.folds' cannot exceed the number of observations")

    n.bws <- as.integer(n.bws)
    if (n.bws < 2) stop("'n.bws' must be at least 2")

    max.iterations <- as.integer(max.iterations)
    if (max.iterations <= 0) stop("'max.iterations' must be positive")

    ## ====================================================================
    ## Call C++ implementation
    ## ====================================================================

    # Call the compiled C++ function
    res <- .Call("S_magelog",
                 x,
                 y,
                 grid.size,
                 fit.quadratic,
                 pilot.bandwidth,
                 kernel,
                 min.points,
                 cv.folds,
                 n.bws,
                 min.bw.factor,
                 max.bw.factor,
                 max.iterations,
                 ridge.lambda,
                 tolerance,
                 with.bw.predictions)

    ## ====================================================================
    ## Prepare return object
    ## ====================================================================

    # Add original data to results
    res$x <- x
    res$y <- y

    # Add fitting information
    res$fit.info <- list(
        grid.size = grid.size,
        fit.quadratic = fit.quadratic,
        kernel = kernel,
        min.points = min.points,
        cv.folds = cv.folds,
        ridge.lambda = ridge.lambda,
        tolerance = tolerance
    )

    # Assign class
    class(res) <- "magelog"

    return(res)
}

#' Predict Method for magelog Objects
#'
#' @description
#' Obtains predictions from a fitted magelog model. Can predict at the original
#' data points or at new x values.
#'
#' @param object A fitted model object of class \code{"magelog"}
#' @param newdata Numeric vector of new x values for prediction. If \code{NULL}
#'   (default), predictions are returned for the original x values.
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"response"}}{Predicted probabilities (default)}
#'     \item{\code{"logit"}}{Predicted values on the logit scale}
#'   }
#' @param ... Additional arguments (currently ignored)
#'
#' @return A numeric vector of predictions
#'
#' @details
#' Predictions are obtained by linear interpolation from the grid-based predictions
#' computed during model fitting. For x values outside the range of the original data,
#' predictions are extrapolated using the nearest grid point values (i.e., constant
#' extrapolation).
#'
#' @examples
#' \dontrun{
#' # Fit model
#' fit <- magelog(x, y)
#'
#' # Predictions at original points
#' pred <- predict(fit)
#'
#' # Predictions at new points
#' x_new <- seq(min(x), max(x), length.out = 50)
#' pred_new <- predict(fit, newdata = x_new)
#' }
#'
#' @export
predict.magelog <- function(object, newdata = NULL, type = c("response", "logit"), ...) {
    type <- match.arg(type)

    # Use original x if newdata not provided
    if (is.null(newdata)) {
        if (type == "response") {
            return(object$predictions)
        } else {
            # Convert to logit scale
            p <- object$predictions
            # Handle boundary cases
            p[p == 0] <- .Machine$double.eps
            p[p == 1] <- 1 - .Machine$double.eps
            return(log(p / (1 - p)))
        }
    }

    # Validate newdata
    if (!is.numeric(newdata)) stop("'newdata' must be numeric")
    newdata <- as.numeric(newdata)
    if (anyNA(newdata) || any(!is.finite(newdata))) {
        stop("'newdata' contains NA, NaN, or infinite values")
    }

    # Get optimal predictions on grid
    opt_idx <- object$opt.brier.bw.idx
    grid_pred <- object$bw.grid.predictions[, opt_idx]

    # Interpolate to new x values
    pred <- approx(x = object$x.grid, y = grid_pred, xout = newdata,
                   method = "linear", rule = 2)$y

    if (type == "logit") {
        # Convert to logit scale
        pred[pred == 0] <- .Machine$double.eps
        pred[pred == 1] <- 1 - .Machine$double.eps
        pred <- log(pred / (1 - pred))
    }

    return(pred)
}

#' Print Method for magelog Objects
#'
#' @description
#' Prints a summary of a fitted magelog model.
#'
#' @param x A fitted model object of class \code{"magelog"}
#' @param digits Integer; number of digits to display for numeric values
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.magelog <- function(x, digits = 4, ...) {
    cat("\nModel-Averaged Bandwidth Logistic Regression (magelog)\n")
    cat(rep("-", 50), "\n", sep = "")

    cat("Number of observations:", length(x$x), "\n")
    cat("Grid size:", x$fit.info$grid.size, "\n")
    cat("Model type:", if(x$fit.info$fit.quadratic) "Quadratic" else "Linear", "\n")
    cat("Kernel:", c("Epanechnikov", "Triangular", "", "Laplace",
                     "Normal", "Biweight", "Tricube")[x$fit.info$kernel], "\n")

    cat("\nBandwidth selection:\n")
    cat("  Method: ", x$fit.info$cv.folds, "-fold cross-validation\n", sep = "")
    cat("  Bandwidths tried:", length(x$bws), "\n")
    cat("  Optimal bandwidth:", round(x$bws[x$opt.brier.bw.idx], digits), "\n")
    cat("  Optimal CV Brier score:", round(x$mean.brier.errors[x$opt.brier.bw.idx], digits), "\n")

    invisible(x)
}

#' Summary Method for magelog Objects
#'
#' @description
#' Produces a summary of a fitted magelog model including goodness-of-fit statistics.
#'
#' @param object A fitted model object of class \code{"magelog"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return An object of class \code{"summary.magelog"} containing summary statistics
#'
#' @export
summary.magelog <- function(object, ...) {
    # Calculate Brier score on training data
    brier_score <- mean((object$predictions - object$y)^2)

    # Calculate log-likelihood
    p <- object$predictions
    p[p == 0] <- .Machine$double.eps
    p[p == 1] <- 1 - .Machine$double.eps
    log_lik <- sum(object$y * log(p) + (1 - object$y) * log(1 - p))

    # Calculate AIC (2 parameters per grid point as rough approximation)
    # This is a heuristic since local regression doesn't have fixed df
    approx_df <- 2 * object$fit.info$grid.size / object$bws[object$opt.brier.bw.idx]
    aic <- -2 * log_lik + 2 * approx_df

    res <- list(
        n = length(object$x),
        grid.size = object$fit.info$grid.size,
        model.type = if(object$fit.info$fit.quadratic) "Quadratic" else "Linear",
        kernel = c("Epanechnikov", "Triangular", "", "Laplace",
                   "Normal", "Biweight", "Tricube")[object$fit.info$kernel],
        optimal.bandwidth = object$bws[object$opt.brier.bw.idx],
        cv.brier.score = object$mean.brier.errors[object$opt.brier.bw.idx],
        training.brier.score = brier_score,
        log.likelihood = log_lik,
        aic = aic,
        approx.df = approx_df
    )

    class(res) <- "summary.magelog"
    return(res)
}

#' Print Method for summary.magelog Objects
#'
#' @description
#' Prints the summary of a magelog model.
#'
#' @param x A summary object of class \code{"summary.magelog"}
#' @param digits Integer; number of digits to display for numeric values
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.summary.magelog <- function(x, digits = 4, ...) {
    cat("\nSummary of Model-Averaged Bandwidth Logistic Regression\n")
    cat(rep("=", 55), "\n", sep = "")

    cat("\nModel Information:\n")
    cat("  Observations:", x$n, "\n")
    cat("  Grid size:", x$grid.size, "\n")
    cat("  Model type:", x$model.type, "\n")
    cat("  Kernel:", x$kernel, "\n")

    cat("\nBandwidth Selection:\n")
    cat("  Optimal bandwidth:", round(x$optimal.bandwidth, digits), "\n")
    cat("  CV Brier score:", round(x$cv.brier.score, digits), "\n")

    cat("\nGoodness of Fit:\n")
    cat("  Training Brier score:", round(x$training.brier.score, digits), "\n")
    cat("  Log-likelihood:", round(x$log.likelihood, digits), "\n")
    cat("  AIC:", round(x$aic, digits), "\n")
    cat("  Approx. degrees of freedom:", round(x$approx.df, 1), "\n")

    invisible(x)
}

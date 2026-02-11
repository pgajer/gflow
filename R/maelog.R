#' Model-Averaged Local Logistic Regression
#'
#' @description
#' Performs model-averaged local logistic regression using kernel-weighted maximum
#' likelihood estimation. The function implements bandwidth selection via cross-validation
#' and supports both linear and quadratic local polynomial models.
#'
#' @param x Numeric vector of predictor variables (explanatory variable).
#' @param y Binary response vector containing only 0 or 1 values.
#' @param fit.quadratic Logical; if \code{TRUE}, fits local quadratic models;
#'   if \code{FALSE} (default), fits local linear models.
#' @param pilot.bandwidth Numeric; the bandwidth parameter for local fitting.
#'   If less than or equal to 0 (default = -1), bandwidth is selected automatically
#'   using cross-validation or LOOCV approximation.
#' @param kernel Integer specifying the kernel function for weight calculation:
#'   \describe{
#'     \item{1}{Epanechnikov kernel}
#'     \item{2}{Triangular kernel}
#'     \item{4}{Laplace (double exponential) kernel}
#'     \item{5}{Gaussian (normal) kernel}
#'     \item{6}{Biweight (quartic) kernel}
#'     \item{7}{Tricube kernel (default)}
#'   }
#' @param min.points Integer; minimum number of points required in each local
#'   neighborhood for fitting. If \code{NULL} (default), automatically set to 3 for
#'   linear models or 4 for quadratic models. Must be at least 3 for linear and
#'   4 for quadratic models.
#' @param cv.folds Integer; number of folds for cross-validation when selecting
#'   bandwidth. If 0 (default), uses leave-one-out cross-validation (LOOCV)
#'   approximation. Must be between 0 and \code{length(x)}.
#' @param n.bws Integer; number of bandwidth candidates to evaluate during
#'   automatic selection. Default is 50. Must be at least 2.
#' @param min.bw.factor Numeric between 0 and 1; minimum bandwidth as a fraction
#'   of the data range. Default is 0.05. The minimum bandwidth tested is
#'   \code{min.bw.factor * diff(range(x))}.
#' @param max.bw.factor Numeric greater than \code{min.bw.factor}; maximum
#'   bandwidth as a fraction of the data range. Default is 0.9. The maximum
#'   bandwidth tested is \code{max.bw.factor * diff(range(x))}.
#' @param max.iterations Integer; maximum number of iterations for the local
#'   likelihood maximization algorithm. Default is 100. Must be positive.
#' @param ridge.lambda Numeric; ridge penalty parameter for numerical stability
#'   in the local regression. Default is 1e-6. Must be positive.
#' @param tolerance Numeric; convergence tolerance for the iterative fitting
#'   algorithm. Default is 1e-8. Must be positive.
#' @param with.errors Logical; if \code{TRUE}, computes and returns standard
#'   errors for the predictions. Default is \code{FALSE}.
#' @param with.bw.predictions Logical; if \code{TRUE} (default), returns the
#'   matrix of predictions for all bandwidth values tested during selection.
#'   Set to \code{FALSE} to save memory when only the optimal predictions are needed.
#'
#' @details
#' The function implements a local likelihood approach to logistic regression,
#' where a weighted logistic model is fit in a neighborhood around each point.
#' The weights are determined by a kernel function and bandwidth parameter.
#'
#' \subsection{Model Specification}{
#' At each point \eqn{x_0}, the local log-odds are modeled as:
#' \itemize{
#'   \item Linear: \eqn{\log\left(\frac{p(x)}{1-p(x)}\right) = \beta_0 + \beta_1(x - x_0)}
#'   \item Quadratic: \eqn{\log\left(\frac{p(x)}{1-p(x)}\right) = \beta_0 + \beta_1(x - x_0) + \beta_2(x - x_0)^2}
#' }
#' where \eqn{p(x) = P(Y = 1|X = x)}.
#' }
#'
#' \subsection{Bandwidth Selection}{
#' When \code{pilot.bandwidth <= 0}, the function automatically selects the
#' bandwidth by minimizing the cross-validation error. The candidate bandwidths
#' are logarithmically spaced between \code{min.bw.factor} and \code{max.bw.factor}
#' times the data range.
#' }
#'
#' \subsection{Kernel Functions}{
#' All kernel functions \eqn{K(u)} are defined on the interval \eqn{[-1, 1]} and
#' are zero outside this interval:
#' \itemize{
#'   \item Epanechnikov: \eqn{K(u) = \frac{3}{4}(1 - u^2)}
#'   \item Triangular: \eqn{K(u) = 1 - |u|}
#'   \item Laplace: \eqn{K(u) = \frac{1}{2}\exp(-|u|)}
#'   \item Gaussian: \eqn{K(u) = \frac{1}{\sqrt{2\pi}}\exp(-\frac{u^2}{2})}
#'   \item Biweight: \eqn{K(u) = \frac{15}{16}(1 - u^2)^2}
#'   \item Tricube: \eqn{K(u) = \frac{70}{81}(1 - |u|^3)^3}
#' }
#' }
#'
#' @return A list of class \code{"maelog"} containing:
#' \item{predictions}{Numeric vector of predicted probabilities at the input points.}
#' \item{errors}{Numeric vector of standard errors if \code{with.errors = TRUE};
#'   \code{NULL} otherwise.}
#' \item{opt.bw}{The selected optimal bandwidth if \code{pilot.bandwidth <= 0};
#'   the input bandwidth otherwise.}
#' \item{candidate.bandwidths}{Numeric vector of bandwidth values tested during
#'   selection if \code{pilot.bandwidth <= 0}; \code{NULL} otherwise.}
#' \item{mean.errors}{Numeric vector of cross-validation errors corresponding to
#'   \code{candidate.bandwidths} if \code{pilot.bandwidth <= 0}; \code{NULL} otherwise.}
#' \item{bw.predictions}{Matrix where column \code{i} contains predictions using
#'   \code{candidate.bandwidths[i]} if \code{with.bw.predictions = TRUE} and
#'   \code{pilot.bandwidth <= 0}; \code{NULL} otherwise.}
#' \item{x}{The input predictor values.}
#' \item{y}{The input response values.}
#' \item{call}{The matched function call.}
#' \item{kernel}{The kernel type used.}
#' \item{fit.quadratic}{Whether quadratic terms were included.}
#'
#' @seealso
#' \code{\link{predict.maelog}} for predictions at new points,
#' \code{\link{plot.maelog}} for diagnostic plots
#'
#' @references
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. \emph{Journal of the American Statistical Association}, 74(368), 829-836.
#'
#' Fan, J., & Gijbels, I. (1996). \emph{Local Polynomial Modelling and Its
#' Applications}. Chapman & Hall/CRC.
#'
#' Loader, C. (1999). \emph{Local Regression and Likelihood}. Springer.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' x <- seq(0, 1, length.out = n)
#' true.prob <- plogis(10 * (x - 0.5))  # Logistic function
#' y <- rbinom(n, 1, true.prob)
#'
#' # Fit model with automatic bandwidth selection
#' fit1 <- maelog(x, y)
#'
#' # Plot results
#' plot(x, y, col = rgb(0, 0, 0, 0.5), pch = 16,
#'      main = "Local Logistic Regression",
#'      xlab = "x", ylab = "Probability")
#' lines(x, fit1$predictions, col = "blue", lwd = 2)
#' lines(x, true.prob, col = "red", lwd = 2, lty = 2)
#' legend("topleft", c("Fitted", "True"),
#'        col = c("blue", "red"), lwd = 2, lty = c(1, 2))
#'
#' # Fit with quadratic local models
#' fit2 <- maelog(x, y, fit.quadratic = TRUE, cv.folds = 5)
#'
#' # Fit with fixed bandwidth
#' fit3 <- maelog(x, y, pilot.bandwidth = 0.1)
#'
#' # Compare different kernels
#' fit.tricube <- maelog(x, y, kernel = 7)
#' fit.gauss <- maelog(x, y, kernel = 5)
#'
#' # Larger example with cross-validation
#' n <- 1000
#' x <- runif(n, -2, 2)
#' prob <- plogis(2 * sin(pi * x))
#' y <- rbinom(n, 1, prob)
#'
#' # Compare linear vs quadratic with 10-fold CV
#' fit.linear <- maelog(x, y, fit.quadratic = FALSE, cv.folds = 10)
#' fit.quad <- maelog(x, y, fit.quadratic = TRUE, cv.folds = 10)
#'
#' # Plot bandwidth selection results
#' par(mfrow = c(1, 2))
#' plot(fit.linear$candidate.bandwidths, fit.linear$mean.errors,
#'      type = "l", xlab = "Bandwidth", ylab = "CV Error",
#'      main = "Linear Model")
#' abline(v = fit.linear$opt.bw, col = "red", lty = 2)
#'
#' plot(fit.quad$candidate.bandwidths, fit.quad$mean.errors,
#'      type = "l", xlab = "Bandwidth", ylab = "CV Error",
#'      main = "Quadratic Model")
#' abline(v = fit.quad$opt.bw, col = "red", lty = 2)
#' }
#'
#' @export
maelog <- function(x,
                   y,
                   fit.quadratic = FALSE,
                   pilot.bandwidth = -1,
                   kernel = 7L,
                   min.points = NULL,
                   cv.folds = 0L,
                   n.bws = 50L,
                   min.bw.factor = 0.05,
                   max.bw.factor = 0.9,
                   max.iterations = 100L,
                   ridge.lambda = 1e-6,
                   tolerance = 1e-8,
                   with.errors = FALSE,
                   with.bw.predictions = TRUE) {

    # Store the call
    cl <- match.call()

    # Basic input validation
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector")
    }
    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector")
    }

    # Convert to numeric and check binary constraint
    x <- as.double(x)
    y <- as.double(y)

    if (!all(y %in% c(0, 1))) {
        stop("'y' must contain only binary values (0 and 1)")
    }

    # Check lengths and missing values
    n <- length(x)
    if (n != length(y)) {
        stop("'x' and 'y' must have the same length")
    }
    if (n == 0) {
        stop("Input vectors cannot be empty")
    }
    if (anyNA(x) || any(!is.finite(x))) {
        stop("'x' contains NA, NaN, or infinite values")
    }
    if (anyNA(y) || any(!is.finite(y))) {
        stop("'y' contains NA, NaN, or infinite values")
    }

    # Validate logical parameters
    if (!is.logical(fit.quadratic) || length(fit.quadratic) != 1) {
        stop("'fit.quadratic' must be a single logical value")
    }
    if (!is.logical(with.errors) || length(with.errors) != 1) {
        stop("'with.errors' must be a single logical value")
    }
    if (!is.logical(with.bw.predictions) || length(with.bw.predictions) != 1) {
        stop("'with.bw.predictions' must be a single logical value")
    }

    # Set and validate min.points
    required.min.points <- if (fit.quadratic) 4L else 3L
    if (is.null(min.points)) {
        min.points <- required.min.points
    } else {
        if (!is.numeric(min.points) || length(min.points) != 1) {
            stop("'min.points' must be a single numeric value")
        }
        min.points <- as.integer(min.points)
        if (min.points < required.min.points) {
            stop(sprintf("'min.points' must be at least %d for %s model",
                        required.min.points,
                        if (fit.quadratic) "quadratic" else "linear"))
        }
    }

    # Check dataset size
    if (n < min.points) {
        stop(sprintf("Dataset must contain at least %d points", min.points))
    }

    # Validate numeric parameters
    if (!is.numeric(pilot.bandwidth) || length(pilot.bandwidth) != 1) {
        stop("'pilot.bandwidth' must be a single numeric value")
    }
    pilot.bandwidth <- as.double(pilot.bandwidth)

    # Validate kernel
    if (!is.numeric(kernel) || length(kernel) != 1) {
        stop("'kernel' must be a single numeric value")
    }
    kernel <- as.integer(kernel)
    if (!kernel %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel' must be one of: 1 (Epanechnikov), 2 (Triangular), ",
             "4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    # Validate cross-validation parameters
    if (!is.numeric(cv.folds) || length(cv.folds) != 1) {
        stop("'cv.folds' must be a single numeric value")
    }
    cv.folds <- as.integer(cv.folds)
    if (cv.folds < 0) {
        stop("'cv.folds' must be non-negative")
    }
    if (cv.folds > n) {
        stop("'cv.folds' cannot exceed the number of observations")
    }

    if (!is.numeric(n.bws) || length(n.bws) != 1) {
        stop("'n.bws' must be a single numeric value")
    }
    n.bws <- as.integer(n.bws)
    if (n.bws < 2) {
        stop("'n.bws' must be at least 2")
    }

    # Validate bandwidth factors
    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1) {
        stop("'min.bw.factor' must be a single numeric value")
    }
    if (min.bw.factor <= 0 || min.bw.factor >= 1) {
        stop("'min.bw.factor' must be between 0 and 1 (exclusive)")
    }

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1) {
        stop("'max.bw.factor' must be a single numeric value")
    }
    if (max.bw.factor <= min.bw.factor) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor'")
    }
    if (max.bw.factor > 1) {
        warning("'max.bw.factor' > 1 may lead to oversmoothing")
    }

    # Validate algorithm parameters
    if (!is.numeric(max.iterations) || length(max.iterations) != 1) {
        stop("'max.iterations' must be a single numeric value")
    }
    max.iterations <- as.integer(max.iterations)
    if (max.iterations <= 0) {
        stop("'max.iterations' must be positive")
    }

    if (!is.numeric(ridge.lambda) || length(ridge.lambda) != 1) {
        stop("'ridge.lambda' must be a single numeric value")
    }
    ridge.lambda <- as.double(ridge.lambda)
    if (ridge.lambda <= 0) {
        stop("'ridge.lambda' must be positive")
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1) {
        stop("'tolerance' must be a single numeric value")
    }
    tolerance <- as.double(tolerance)
    if (tolerance <= 0) {
        stop("'tolerance' must be positive")
    }

    # Call the C++ implementation
    res <- .Call("S_maelog",
                 x,
                 y,
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
                 with.errors,
                 with.bw.predictions)

    # Add input data and metadata to results
    res$x <- x
    res$y <- y
    res$call <- cl
    res$kernel <- kernel
    res$fit.quadratic <- fit.quadratic

    # Set class for S3 methods
    class(res) <- "maelog"

    return(res)
}

#' Predict Method for maelog Objects
#'
#' @description
#' Obtains predictions from a fitted \code{maelog} object at new predictor values.
#'
#' @param object A fitted object of class \code{"maelog"}.
#' @param newdata Numeric vector of new predictor values at which to make predictions.
#'   If missing, predictions at the original data points are returned.
#' @param type Character string specifying the type of prediction:
#'   \code{"response"} for predicted probabilities (default) or
#'   \code{"link"} for linear predictors on the logit scale.
#' @param se.fit Logical; if \code{TRUE}, returns standard errors of predictions
#'   (only available if the model was fit with \code{with.errors = TRUE}).
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' If \code{se.fit = FALSE}, a numeric vector of predictions.
#' If \code{se.fit = TRUE}, a list with components:
#' \item{fit}{Vector of predictions.}
#' \item{se.fit}{Vector of standard errors (if available).}
#'
#' @seealso \code{\link{maelog}}
#'
#' @examples
#' \dontrun{
#' # Fit model
#' set.seed(123)
#' x <- seq(0, 1, length.out = 100)
#' y <- rbinom(100, 1, plogis(10 * (x - 0.5)))
#' fit <- maelog(x, y, with.errors = TRUE)
#'
#' # Predictions at new points
#' x.new <- seq(0.2, 0.8, by = 0.1)
#' pred <- predict(fit, x.new)
#'
#' # Predictions with standard errors
#' pred.se <- predict(fit, x.new, se.fit = TRUE)
#' }
#' @export
predict.maelog <- function(object, newdata, type = c("response", "link"),
                          se.fit = FALSE, ...) {
    type <- match.arg(type)

    if (missing(newdata)) {
        # Return fitted values
        fit <- if (type == "response") {
            object$predictions
        } else {
            qlogis(object$predictions)
        }

        if (se.fit) {
            if (is.null(object$errors)) {
                warning("Standard errors not available; model was fit with with.errors = FALSE")
                return(list(fit = fit, se.fit = NULL))
            }
            return(list(fit = fit, se.fit = object$errors))
        }
        return(fit)
    }

    # Validate newdata
    if (!is.numeric(newdata)) {
        stop("'newdata' must be numeric")
    }
    newdata <- as.double(newdata)
    if (any(!is.finite(newdata))) {
        stop("'newdata' contains non-finite values")
    }

    # For new data, would need to call C++ prediction function
    # This is a placeholder - actual implementation would call C++
    stop("Prediction at new points not yet implemented")
}

#' Plot Method for maelog Objects
#'
#' @description
#' Produces diagnostic plots for fitted \code{maelog} objects.
#'
#' @param x A fitted object of class \code{"maelog"}.
#' @param which Integer specifying which plot to produce:
#'   \describe{
#'     \item{1}{Fitted values vs predictor with data points (default)}
#'     \item{2}{Cross-validation error vs bandwidth (if bandwidth was selected)}
#'     \item{3}{Residuals vs fitted values}
#'     \item{4}{QQ-plot of Pearson residuals}
#'   }
#' @param main Character string for plot title. Default title is used if \code{""}.
#' @param ... Additional graphical parameters passed to the plotting functions.
#'
#' @return
#' Invisibly returns the input object.
#'
#' @note
#' Plot 2 (bandwidth selection) is only available when the model was fitted
#' with automatic bandwidth selection (i.e., \code{pilot.bandwidth <= 0}).
#'
#' To display multiple plots, use \code{par(mfrow)} or call \code{plot} multiple times:
#' \preformatted{
#' par(mfrow = c(2, 2))
#' for (i in 1:4) plot(fit, which = i)
#' }
#'
#' @seealso \code{\link{maelog}}, \code{\link{summary.maelog}}
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- rbinom(n, 1, plogis(10 * (x - 0.5)))
#'
#' # Fit with fixed bandwidth
#' fit1 <- maelog(x, y, pilot.bandwidth = 0.1)
#'
#' # Plot fitted values
#' plot(fit1)
#'
#' \dontrun{
#' # Plot residuals
#' plot(fit1, which = 3)
#'
#' # Plot QQ plot
#' plot(fit1, which = 4)
#'
#' # Fit with automatic bandwidth selection
#' fit2 <- maelog(x, y, pilot.bandwidth = -1, n.bws = 30)
#'
#' # Show all plots
#' par(mfrow = c(2, 2))
#' plot(fit2, which = 1)  # Fitted values
#' plot(fit2, which = 2)  # Bandwidth selection
#' plot(fit2, which = 3)  # Residuals
#' plot(fit2, which = 4)  # QQ plot
#' par(mfrow = c(1, 1))
#' }
#'
#' @importFrom grDevices rgb
#' @importFrom graphics abline lines mtext par plot
#' @importFrom stats lowess plogis qlogis qqline qqnorm
#' @export
plot.maelog <- function(x, which = 1L, main = "", ...) {
    if (!inherits(x, "maelog")) {
        stop("'x' must be a maelog object")
    }

    # Validate which
    if (!is.numeric(which) || length(which) != 1) {
        stop("'which' must be a single integer")
    }
    which <- as.integer(which)
    if (which < 1L || which > 4L) {
        stop("'which' must be an integer between 1 and 4")
    }

    # Set default titles
    default.mains <- c("Fitted Values", "Bandwidth Selection",
                      "Residuals vs Fitted", "Normal Q-Q Plot")
    if (main == "") {
        main <- default.mains[which]
    }

    # Plot 1: Fitted values
    if (which == 1L) {
        plot(x$x, x$y, col = rgb(0, 0, 0, 0.5), pch = 16,
             xlab = "x", ylab = "Probability",
             main = main, ylim = c(0, 1), ...)
        o <- order(x$x)
        lines(x$x[o], x$bw_predictions[o], col = "blue", lwd = 2)
    }

    # Plot 2: CV error vs bandwidth
    else if (which == 2L) {
        if (is.null(x$candidate.bandwidths)) {
            stop("Plot 2 requires automatic bandwidth selection (pilot.bandwidth <= 0)")
        }
        plot(x$candidate.bandwidths, x$mean.errors,
             type = "b", pch = 19,
             xlab = "Bandwidth", ylab = "Cross-validation Error",
             main = main, ...)
        abline(v = x$opt.bw, col = "red", lty = 2)
        mtext(sprintf("Optimal: %.3f", x$opt.bw),
              side = 3, line = 0.5, col = "red", cex = 0.8)
    }

    # Plot 3: Residuals vs fitted
    else if (which == 3L) {
        residuals <- x$y - x$predictions
        plot(x$predictions, residuals,
             xlab = "Fitted Values", ylab = "Residuals",
             main = main, xlim = c(0, 1),
             pch = 16, col = rgb(0, 0, 0, 0.5), ...)
        abline(h = 0, col = "red", lty = 2)

        # Only add lowess line if we have enough unique x values and no issues
        if (length(unique(x$predictions)) > 10) {
            # Check for finite values before calling lowess
            finite.idx <- is.finite(x$predictions) & is.finite(residuals)
            if (sum(finite.idx) > 10) {
                tryCatch({
                    lo <- lowess(x$predictions[finite.idx], residuals[finite.idx])
                    lines(lo, col = "blue", lwd = 2)
                }, error = function(e) {
                    # Silently skip the lowess line if it fails
                })
            }
        }
    }

    # Plot 4: QQ plot of Pearson residuals
    else if (which == 4L) {
        pred.bounded <- pmax(pmin(x$predictions, 1 - 1e-10), 1e-10)
        pearson.resid <- (x$y - x$predictions) /
                        sqrt(pred.bounded * (1 - pred.bounded))
        finite.resid <- pearson.resid[is.finite(pearson.resid)]

        if (length(finite.resid) == 0) {
            stop("No finite Pearson residuals to plot")
        }

        if (length(finite.resid) < length(pearson.resid)) {
            warning(sprintf("%d non-finite residuals omitted from QQ plot",
                           length(pearson.resid) - length(finite.resid)))
        }
        qqnorm(finite.resid, main = main,
               ylab = "Pearson Residuals", ...)
        qqline(finite.resid, col = "red")
    }

    invisible(x)
}

#' Summary Method for maelog Objects
#'
#' @description
#' Produces a summary of a fitted \code{maelog} object.
#'
#' @param object A fitted object of class \code{"maelog"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns a list containing summary statistics. The function is called
#' primarily for its side effect of printing the summary.
#'
#' @seealso \code{\link{maelog}}
#'
#' @export
summary.maelog <- function(object, ...) {
    if (!inherits(object, "maelog")) {
        stop("'object' must be a maelog object")
    }

    cat("\nModel-Averaged Local Logistic Regression\n")
    cat("\nCall:\n")
    print(object$call)

    cat("\nNumber of observations:", length(object$x))
    cat("\nResponse: ", sum(object$y), " successes, ",
        sum(1 - object$y), " failures\n", sep = "")

    cat("\nModel type:", if (object$fit.quadratic) "Quadratic" else "Linear")

    kernel.names <- c("Epanechnikov", "Triangular", "", "Laplace",
                     "Gaussian", "Biweight", "Tricube")
    cat("\nKernel:", kernel.names[object$kernel])

    cat("\nBandwidth:", formatC(object$opt.bw, digits = 4))
    if (!is.null(object$candidate.bandwidths)) {
        cat(" (selected from", length(object$candidate.bandwidths), "candidates)")
    }

    # Compute deviance
    pred <- pmax(pmin(object$predictions, 1 - 1e-15), 1e-15)
    deviance <- -2 * sum(object$y * log(pred) + (1 - object$y) * log(1 - pred))
    null.deviance <- -2 * sum(object$y * log(mean(object$y)) +
                             (1 - object$y) * log(1 - mean(object$y)))

    cat("\n\nNull deviance:    ", formatC(null.deviance, digits = 4),
        " on ", length(object$y) - 1, " degrees of freedom", sep = "")
    cat("\nResidual deviance:", formatC(deviance, digits = 4))

    # McFadden's pseudo R-squared
    pseudo.r2 <- 1 - deviance / null.deviance
    cat("\nPseudo R-squared: ", formatC(pseudo.r2, digits = 4), "\n", sep = "")

    invisible(list(call = object$call,
                  n = length(object$x),
                  kernel = kernel.names[object$kernel],
                  bandwidth = object$opt.bw,
                  deviance = deviance,
                  null.deviance = null.deviance,
                  pseudo.r2 = pseudo.r2))
}

#' Print Method for maelog Objects
#'
#' @description
#' Prints a brief summary of a fitted \code{maelog} object.
#'
#' @param x A fitted object of class \code{"maelog"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input object.
#'
#' @seealso \code{\link{maelog}}, \code{\link{summary.maelog}}
#'
#' @export
print.maelog <- function(x, ...) {
    cat("\nModel-Averaged Local Logistic Regression\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nBandwidth:", formatC(x$opt.bw, digits = 4), "\n")
    invisible(x)
}

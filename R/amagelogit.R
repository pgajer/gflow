#' Grid-based Model Averaged Bandwidth Logistic Regression
#'
#' @description
#' Performs model-averaged bandwidth logistic regression using local polynomial fitting.
#' The function implements a flexible approach to binary regression by fitting local
#' models at points of a uniform grid spanning the range of predictor values. Multiple
#' local models are combined to produce robust predictions. It supports both linear and
#' quadratic local models, automatic bandwidth selection, and various kernel types
#' for weight calculation.
#'
#' @param x Numeric vector of predictor variables
#' @param y Binary vector (0 or 1) of response variables
#' @param grid.size Integer; number of points in the uniform grid where local models
#'        are centered. Default is 200
#' @param fit.quadratic Logical; whether to include quadratic terms in the local models.
#'        Default is FALSE
#' @param pilot.bandwidth Numeric; bandwidth for local fitting. If <= 0, bandwidth
#'        is automatically selected. Default is -1
#' @param kernel Integer; kernel type for weight calculation:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#' @param min.points Integer; minimum number of points required for local fitting.
#'        Default is automatically set based on fit.quadratic (4 for quadratic, 3 for linear)
#' @param cv.folds Integer; number of cross-validation folds for bandwidth selection.
#'        If 0, LOOCV approximation is used. Default is 0
#' @param n.bws Integer; number of bandwidths to try in automatic selection. Default is 50
#' @param min.bw.factor Numeric; minimum bandwidth factor relative to data range. Default is 0.05
#' @param max.bw.factor Numeric; maximum bandwidth factor relative to data range. Default is 0.9
#' @param max.iterations Integer; maximum number of iterations for local fitting. Default is 100
#' @param ridge.lambda Numeric; ridge parameter for local fitting. Default is 1e-6
#' @param tolerance Numeric; convergence tolerance for local fitting. Default is 1e-8
#' @param with.bw.predictions Logical; whether to return predictions for all bandwidths. Default is TRUE
#'
#' @return A list containing:
#' \itemize{
#'   \item x.grid: Uniform grid points where local models are centered
#'   \item predictions: Predicted probabilities at original x points
#'   \item bw.grid.predictions: Matrix of predictions at grid points for each bandwidth
#'   \item mean.brier.errors: Cross-validation errors for each bandwidth
#'   \item opt.brier.bw.idx: Index of optimal bandwidth
#'   \item bws: Vector of tried bandwidths
#'   \item fit.info: List of fitting parameters used
#'   \item x: Original predictor values
#'   \item y: Original response values
#' }
#'
#' @details
#' The function fits local logistic regression models centered at points of a uniform
#' grid spanning the range of x values. Local models are fit using kernel-weighted
#' maximum likelihood. The bandwidth determines the size of the local neighborhood.
#' When pilot.bandwidth <= 0, the function automatically selects a bandwidth using
#' cross-validation.
#'
#' The local models can be either linear or quadratic (controlled by fit.quadratic).
#' For numerical stability, the minimum number of points in each local fit is automatically
#' set to 3 for linear and 4 for quadratic models, but can be overridden with min.points.
#'
#' Predictions at the original x points are obtained by linear interpolation from the
#' grid-based predictions. The grid approach provides computational efficiency and
#' smooth prediction curves.
#'
#' @examples
#' \dontrun{
#' x <- seq(0, 1, length.out = 100)
#' p <- 1/(1 + exp(-(x - 0.5)*10))
#' y <- rbinom(100, 1, p)
#' fit <- amagelogit(x, y, grid.size = 200, fit.quadratic = TRUE, cv.folds = 5)
#' plot(x, y)
#' lines(fit$x.grid, fit$bw.grid.predictions[,fit$opt.brier.bw.idx], col = "red")
#' }
#'
#' @export
amagelogit <- function(x,
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
                      with.bw.predictions = TRUE ##  parallel = FALSE,  verbose = FALSE
                      ) {

    ## Basic input checks
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    if (!is.numeric(y)) stop("'y' must be a numeric vector")

    ## Convert and check binary values
    y <- as.numeric(y)
    if (!all(y %in% c(0, 1))) stop("'y' must contain only binary values (0 and 1)")

    ## Check vector lengths and missing values
    n <- length(x)
    if (n != length(y)) stop("'x' and 'y' must have the same length")
    if (n == 0) stop("Input vectors cannot be empty")
    if (anyNA(x) || any(is.infinite(x))) stop("'x' contains NA or infinite values")
    if (anyNA(y) || any(is.infinite(y))) stop("'y' contains NA or infinite values")

    ## Check grid.size
    grid.size <- as.integer(grid.size)
    if (grid.size < 2) stop("'grid.size' must be at least 2")

    ## Check logical parameters
    if (!is.logical(fit.quadratic)) stop("'fit.quadratic' must be logical")
    if (!is.logical(with.bw.predictions)) stop("'with.bw.predictions' must be logical")

    ## Set min.points based on model type if not provided
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

    ## Check if dataset has enough points
    if (n < min.points) {
        stop(sprintf("Input vectors must contain at least %d points", min.points))
    }

    ## Check numeric parameters
    if (!is.numeric(pilot.bandwidth)) stop("'pilot.bandwidth' must be numeric")
    if (!is.numeric(min.bw.factor)) stop("'min.bw.factor' must be numeric")
    if (!is.numeric(max.bw.factor)) stop("'max.bw.factor' must be numeric")
    if (!is.numeric(ridge.lambda)) stop("'ridge.lambda' must be numeric")
    if (!is.numeric(tolerance)) stop("'tolerance' must be numeric")

    ## Check and convert integer parameters
    kernel <- as.integer(kernel)
    if (!kernel %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    cv.folds <- as.integer(cv.folds)
    if (cv.folds < 3) stop("'cv.folds' must be at least 3")
    if (cv.folds > n) stop("'cv.folds' cannot exceed the dataset size")

    n.bws <- as.integer(n.bws)
    if (n.bws < 2) stop("'n.bws' must be at least 2")

    max.iterations <- as.integer(max.iterations)
    if (max.iterations <= 0) stop("'max.iterations' must be positive")

    ## Check bandwidth factors
    if (min.bw.factor <= 0) stop("'min.bw.factor' must be positive")
    if (min.bw.factor >= 1) stop("'min.bw.factor' must less than 1")
    if (max.bw.factor <= min.bw.factor) stop("'max.bw.factor' must be greater than min.bw.factor")

    ## Additional numeric parameter checks
    if (ridge.lambda <= 0) stop("'ridge.lambda' must be positive")
    if (tolerance <= 0) stop("'tolerance' must be positive")


    ## Check if x is sorted and sort if necessary
    ## if (!identical(x, sort(x))) {
    ##     ord <- order(x)
    ##     x <- x[ord]
    ##     y <- y[ord]
    ## }

    # Make the call to C++ function
    res <- .Call("S_amagelogit",
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

    res$x <- x
    res$y <- y

    ## result$x_sorted <- x
    ## result$y_sorted <- y

    return(res)
}

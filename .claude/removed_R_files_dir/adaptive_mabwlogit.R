#' Adaptive Bandwidth Local Logistic Regression
#'
#' @description
#' Performs local logistic regression with adaptive bandwidth selection based on local complexity.
#' The function first estimates a global pilot bandwidth using cross-validation, then adapts
#' the bandwidth locally based on the complexity of the fitted surface.
#'
#' @param x Numeric vector of predictor variables. Must be sorted in ascending order.
#' @param y Numeric vector of binary response variables (0 or 1).
#' @param c.min Numeric. Minimum allowed adaptation factor for bandwidth adjustment. Default is 0.2.
#' @param c.max Numeric. Maximum allowed adaptation factor for bandwidth adjustment. Default is 5.0.
#' @param power Numeric. Power for the adaptation formula. Default is -0.2.
#' @param kernel Integer; kernel type for weight calculation:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#'        Default is 7.
#' @param min.points Integer. Minimum number of points required in local neighborhood.
#'        Default is max(6, as.integer(0.05 * length(x))).
#' @param max.iterations Integer. Maximum number of iterations for logistic regression fitting.
#'        Default is 100.
#' @param ridge.lambda Numeric. Ridge regularization parameter. Default is 1e-6.
#' @param tolerance Numeric. Convergence tolerance for model fitting. Default is 1e-8.
#'
#' @return A list of class "mabilog" containing:
#' \itemize{
#'   \item predictions: Initial predictions using pilot bandwidth
#'   \item apredictions: Matrix of predictions using adaptive bandwidths
#'   \item baseline.bw: Initial pilot bandwidth
#'   \item var.bw: Vector of variable bandwidths
#'   \item beta1s: Linear coefficients
#'   \item beta2s: Quadratic coefficients
#'   \item n.points.used: Number of points used in each local fit
#'   \item used.knn: Logical vector indicating where k-NN fallback was used
#'   \item fit.info: List of fitting parameters and information
#' }
#'
#' @details
#' The adaptive bandwidth formula is:
#' \eqn{h(x) = h_0 \times \min\{\max\{(|\beta_2(x)|/\tilde{\beta}_2)^{power}, c_{min}\}, c_{max}\}}
#' where \eqn{h_0} is the pilot bandwidth selected by cross-validation,
#' \eqn{\beta_2(x)} is the local quadratic coefficient, and \eqn{\tilde{\beta}_2} is
#' the median of \eqn{|\beta_2(x)|}.
#'
#' @examples
#' x <- sort(runif(100))
#' y <- rbinom(100, 1, 0.5)
#' fit <- adaptive.mabwlogit(x, y)
#'
#' @export
adaptive.mabwlogit <- function(x,
                               y,
                               c.min = 0.2,
                               c.max = 5.0,
                               power = -0.2,
                               kernel = 7L,
                               min.points = max(6, as.integer(0.05 * length(x))),
                               max.iterations = 100,
                               ridge.lambda = 1e-6,
                               tolerance = 1e-8) {

    # Input validation
    if (!is.numeric(x)) stop("x must be numeric")
    if (!is.numeric(y)) stop("y must be numeric")
    if (length(x) != length(y)) stop("x and y must have the same length")
    if (length(x) < 3) stop("At least 3 points are required")

    # Validate adaptation parameters
    if (!is.numeric(c.min) || c.min <= 0) stop("c.min must be positive")
    if (!is.numeric(c.max) || c.max <= c.min) stop("c.max must be greater than c.min")
    if (!is.numeric(power)) stop("power must be numeric")

    ## Validate other parameters
    kernel <- as.integer(kernel)
    if (!kernel %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }
    if (!is.numeric(min.points) || min.points < 3)
        stop("min.points must be at least 3")
    if (!is.numeric(max.iterations) || max.iterations <= 0)
        stop("max.iterations must be positive")
    if (!is.numeric(ridge.lambda) || ridge.lambda < 0)
        stop("ridge.lambda must be non-negative")
    if (!is.numeric(tolerance) || tolerance <= 0)
        stop("tolerance must be positive")

    # Call the C++ implementation
    result <- .Call("S_adaptive_mabwlogit",
                   as.double(x),
                   as.double(y),
                   as.double(c.min),
                   as.double(c.max),
                   as.double(power),
                   as.integer(kernel),
                   as.integer(min.points),
                   as.integer(max.iterations),
                   as.double(ridge.lambda),
                   as.double(tolerance))

    return(result)
}

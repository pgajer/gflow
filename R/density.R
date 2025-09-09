#' Kernel Density Estimation on a Uniform Grid
#'
#' @description
#' Estimates probability density using kernel density estimation on a uniform grid.
#' Similar to base R's density() but uses a fixed uniform grid for estimation.
#'
#' @param x Numeric vector of observations
#' @param grid.size Integer specifying the number of points in the estimation grid (default: 512)
#' @param poffset Numeric value between 0 and 1 specifying the proportion of data range
#'        to add as padding on each end (default: 0.1)
#' @param bw Bandwidth for kernel estimation. If <= 0, bandwidth is automatically selected
#'        (default: 0)
#' @param kernel Integer specifying kernel type:
#'        1 = Epanechnikov,
#'        2 = Triangular,
#'        4 = Laplace,
#'        5 = Normal,
#'        6 = Biweight,
#'        7 = Tricube
#' @param verbose Logical indicating whether to print progress information (default: FALSE)
#'
#' @return A list containing:
#' \itemize{
#'   \item density - Vector of density estimates at grid points
#'   \item bw - Bandwidth used for estimation
#'   \item bw_auto_selected - Logical indicating if bandwidth was automatically selected
#'   \item offset - Actual offset added to data range
#'   \item start - Starting point of the grid
#'   \item end - End point of the grid
#' }
#'
#' @examples
#' data <- rnorm(100)
#' result <- gdensity(data, grid.size = 200, poffset = 0.1, bw = 0, kernel = 5)
#' plot(result$x, result$density, type = "l")
#'
#' @export
gdensity <- function(x,
                     grid.size = 512L,
                     poffset = 0.1,
                     bw = 0,
                     kernel = 5L,
                     verbose = FALSE) {
    # Check if x is numeric and has data
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector")
    }
    if (length(x) < 1) {
        stop("'x' must not be empty")
    }
    if (any(is.na(x))) {
        stop("'x' contains missing values")
    }
    if (any(is.infinite(x))) {
        stop("'x' contains infinite values")
    }

    # Check grid.size
    if (!is.numeric(grid.size) || length(grid.size) != 1) {
        stop("'grid.size' must be a single numeric value")
    }
    if (grid.size != round(grid.size)) {
        stop("'grid.size' must be an integer")
    }
    if (grid.size < 10) {
        stop("'grid.size' must be at least 10")
    }

    # Check poffset
    if (!is.numeric(poffset) || length(poffset) != 1) {
        stop("'poffset' must be a single numeric value")
    }
    if (poffset < 0 || poffset > 0.5) {
        stop("'poffset' must be between 0 and 0.5")
    }

    # Check bandwidth
    if (!is.numeric(bw) || length(bw) != 1) {
        stop("'bw' must be a single numeric value")
    }

    # Check kernel
    if (!is.numeric(kernel) || length(kernel) != 1) {
        stop("'kernel' must be a single numeric value")
    }
    if (!kernel %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    # Check verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value")
    }

    result <- .Call("S_estimate_local_density_over_grid",
                   as.double(x),
                   as.integer(grid.size),
                   as.double(poffset),
                   as.double(bw),
                   as.integer(kernel),
                   as.logical(verbose))

    result$x <- seq(result$start, result$end, length.out = grid.size)

    return(result)
}

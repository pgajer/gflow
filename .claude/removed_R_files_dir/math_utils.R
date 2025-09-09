#' Calculate Integrals Using Different Methods
#'
#' @description
#' Computes numerical integrals for multiple y-series against a common x-series
#' using either the trapezoidal rule or R's built-in integration method.
#'
#' @param x.series Numeric vector of x coordinates. Must contain at least 2 elements.
#' @param y.list Named list of numeric vectors, each containing y coordinates
#'              corresponding to x.series. Each vector must have the same length
#'              as x.series.
#' @param method Character string specifying the integration method.
#'              Options are "trapezoid" or "integrate" (default: "trapezoid")
#'
#' @return Named numeric vector containing the computed integrals.
#'         Names correspond to the names in y.list.
#'
#' @details
#' The function supports two integration methods:
#' \itemize{
#'   \item trapezoid: Uses trapezoidal rule for numerical integration
#'   \item integrate: Uses R's built-in integrate function with spline interpolation
#' }
#'
#' @examples
#' x <- seq(0, 10, 0.1)
#' y.list <- list(
#'   series1 = sin(x),
#'   series2 = cos(x)
#' )
#' # Using trapezoid method
#' calculate.integrals(x, y.list, method = "trapezoid")
#' # Using R's integrate function
#' calculate.integrals(x, y.list, method = "integrate")
#'
#' @importFrom stats integrate splinefun
#' @export
calculate.integrals <- function(x.series, y.list, method = c("trapezoid", "integrate")) {
    # Parameter validation
    method <- match.arg(method)

    # Input validation and conversion
    if (!is.list(y.list) && !is.matrix(y.list) && !is.data.frame(y.list)) {
        stop("y must be a named list, data frame, or matrix")
    }

    # Convert matrix to data frame
    if (is.matrix(y.list)) {
        ## if (is.null(colnames(y.list))) {
        ##     stop("Matrix input must have column names")
        ## }
        y.list <- as.data.frame(y.list)
    }

    # Convert data frame to list if necessary
    if (is.data.frame(y.list)) {
        y.list <- as.list(y.list)
    }

    if (!is.numeric(x.series) || length(x.series) < 2) {
        stop("x.series must be a numeric vector with at least 2 elements")
    }

    # Validate each y series
    ## lapply(names(y.list), function(name) {
    ##     y <- y.list[[name]]
    ##     if (!is.numeric(y) || length(y) != length(x.series)) {
    ##         stop(sprintf("y series '%s' must be numeric and same length as x.series", name))
    ##     }
    ## })

    # Initialize vector to store integrals
    integrals <- numeric(length(y.list))
    names(integrals) <- names(y.list)

    if (method == "trapezoid") {
        # Compute integrals using trapezoidal rule
        for(i in seq_along(y.list)) {
            y <- y.list[[i]]
            x <- x.series[seq(length(y))]
            integrals[i] <- -sum(diff(x) * (y[-1] + y[-length(y)])/2)
        }
    } else if (method == "integrate") {
        # Compute integrals using R's integrate function with spline interpolation
        for(i in seq_along(y.list)) {
            y <- y.list[[i]]
            x <- x.series[seq(length(y))]
            # Create spline function for interpolation
            spline.fn <- splinefun(x, y, method = "natural")
            # Integrate using R's built-in function
            integral.result <- try(
                integrate(spline.fn, min(x), max(x)),
                silent = TRUE
            )
            if (inherits(integral.result, "try-error")) {
                warning(sprintf("Integration failed for series '%s'. Returning NA.",
                              names(y.list)[i]))
                integrals[i] <- NA
            } else {
                integrals[i] <- integral.result$value
            }
        }
    }

    return(integrals)
}


#' Compute Bayesian Bootstrap Integrals for GAM Fits
#'
#' @description
#' Performs Bayesian bootstrap integration by fitting GAMs with random weights
#' and computing the integral of the resulting fits. Uses parallel processing
#' for improved performance.
#'
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates (response values)
#' @param n.bb Integer. Number of Bayesian bootstrap iterations (default: 1000)
#'
#' @return Numeric vector of length n.bb containing computed integral values
#'
#' @details
#' The function performs the following steps:
#' 1. Generates Dirichlet weights for each bootstrap iteration
#' 2. Fits a GAM using cubic regression splines
#' 3. Computes the integral of the fitted spline
#'
#' @note
#' Parallel backend should be registered before calling this function.
#' Requires the packages mgcv, foreach, and doParallel.
#'
#' @examples
#' \dontrun{
#' library(doParallel)
#' registerDoParallel(cores = 4)
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x)
#' integrals <- bb.integrals(x, y, n.bb = 100)
#' }
#'
#' @importFrom mgcv gam predict.gam
#' @importFrom foreach foreach %dopar%
#' @importFrom stats runif integrate splinefun
#'
bb.integrals <- function(x, y, n.bb = 1000) {
    # Input validation
    if (length(x) != length(y)) {
        stop("x and y must have the same length")
    }
    if (any(is.na(c(x, y)))) {
        stop("x and y cannot contain NA values")
    }
    if (any(is.infinite(c(x, y)))) {
        stop("x and y cannot contain infinite values")
    }

    x.min <- min(x)
    x.max <- max(x)

    # Check if parallel backend is registered
    if (!foreach::getDoParRegistered()) {
        stop("Parallel backend not registered. Please call registerDoParallel() first.")
    }

    bb.ints <- foreach(i = seq(n.bb),
                      .combine = c,
                      .packages = c("mgcv", "MCMCpack")) %dopar% {
        # Generate Dirichlet weights
        lambda <- MCMCpack::rdirichlet(1, rep(1, length(x)))

        # Fit GAM with try-catch
        bb.fit <- try(mgcv::gam(y ~ s(x, bs = "cr"),
                               weights = lambda,
                               method = "REML",  # more stable than GCV
                               select = TRUE),   # helps with numerical stability
                     silent = TRUE)

        if (inherits(bb.fit, "try-error")) {
            return(NA)
        }

        # Predict and create spline with error checking
        bb.Ey <- try(predict(bb.fit, type = "response"), silent = TRUE)

        if (inherits(bb.Ey, "try-error")) {
            return(NA)
        }

        # Compute integral
        spline.fn <- try(splinefun(x, bb.Ey, method = "natural"), silent = TRUE)

        if (inherits(spline.fn, "try-error")) {
            return(NA)
        }

        integral.result <- try(
            integrate(spline.fn, x.min, x.max),
            silent = TRUE
        )

        if (inherits(integral.result, "try-error")) {
            NA
        } else {
            integral.result$value
        }
    }

    return(bb.ints)
}

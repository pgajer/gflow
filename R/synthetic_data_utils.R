##
## Synthetic data generating functions
##

## --------------------------------------------------------------------------------------------------------------
##
## functions generating synthetic data in dimension 1
##
## --------------------------------------------------------------------------------------------------------------

#' Generate X-knot Positions for Gaussian Mixture Components
#'
#' @param strategy Character string specifying the positioning strategy. Must be one of:
#'   * "uniform": equally spaced components
#'   * "random": uniformly random positions
#'   * "jittered": perturbed uniform spacing
#'   * "offset": sequential spacing with sigma-relative gaps
#' @param n.components Integer specifying the number of Gaussian components
#' @param x.range Numeric vector of length 2 specifying the range \eqn{[\text{min}, \text{max}]} for x values
#' @param edge.buffer Numeric between 0 and 0.5 specifying the buffer from edges as fraction of range
#' @param jitter.fraction Numeric between 0 and 1 specifying maximum jitter as fraction of local spacing
#'
#' @return Numeric vector of length n.components containing x-coordinates for Gaussian components
#' @keywords internal
generate_x_knots <- function(strategy, n.components, x.range, edge.buffer, jitter.fraction) {
    x.width <- diff(x.range)
    edge.width <- edge.buffer * x.width

    # Special case for single component
    if (n.components == 1) {
        switch(strategy,
               "uniform" = {
                   return(mean(x.range))  # Center for single component
               },
               "random" = {
                   return(runif(1,
                                min = x.range[1] + edge.width,
                                max = x.range[2] - edge.width))
               },
               "jittered" = {
                   return(mean(x.range))  # Center for single component
               },
               "offset" = {
                   return(0)  # Start at origin for offset strategy
               },
               stop("Invalid x.knots.strategy")
        )
    }

    switch(strategy,
        "uniform" = {
            x.knots <- numeric(n.components)
            x.knots[1] <- x.range[1] + edge.width
            for (j in 2:n.components) {
                x.knots[j] <- x.knots[j-1] + (x.width - 2 * edge.width)/(n.components - 1)
            }
            x.knots
        },
        "random" = {
            sort(runif(n.components,
                      min = x.range[1] + edge.width,
                      max = x.range[2] - edge.width))
        },
        "jittered" = {
            # Start with uniform spacing
            x.knots <- numeric(n.components)
            x.knots[1] <- x.range[1] + edge.width
            for (j in 2:n.components) {
                x.knots[j] <- x.knots[j-1] + (x.width - 2 * edge.width)/(n.components - 1)
            }

            # Apply jitter
            for (j in 1:n.components) {
                if (j == 1) {
                    x.left <- x.knots[1] - x.range[1]
                    x.right <- x.knots[2] - x.knots[1]
                    jitter.range <- c(-x.left, jitter.fraction * x.right)
                } else if (j == n.components) {
                    x.left <- x.knots[j] - x.knots[j-1]
                    x.right <- x.range[2] - x.knots[j]
                    jitter.range <- c(-jitter.fraction * x.left, x.right)
                } else {
                    x.left <- x.knots[j] - x.knots[j-1]
                    x.right <- x.knots[j+1] - x.knots[j]
                    jitter.range <- c(-jitter.fraction * x.left,
                                    jitter.fraction * x.right)
                }
                x.knots[j] <- x.knots[j] + runif(1, min = jitter.range[1],
                                                max = jitter.range[2])
            }
            x.knots
        },
        "offset" = {
            # Note: This strategy ignores x.range and edge.buffer
            # The calling function should handle x.range adjustment
            return(NULL)  # Signal to handle in main function
        },
        stop("Invalid x.knots.strategy")
    )
}

#' Generate Standard Deviations for Gaussian Mixture Components
#'
#' @param strategy Character string specifying the SD generation strategy. Must be one of:
#'   "fixed" (constant SD), "adaptive" (based on neighbor distances),
#'   "stoch.adaptive" (randomized adaptive), or "custom" (user-specified)
#' @param x.knots Numeric vector of component x-coordinates
#' @param sd.knot Numeric specifying the fixed SD value when strategy = "fixed"
#' @param sd.factor Numeric between 0 and 1 specifying the fraction of neighbor distance to use for SD
#' @param sd.random.range Numeric vector of length 2 specifying \eqn{[\text{min}, \text{max}]} multipliers for random SD
#' @param custom.sd.knots Optional numeric vector of custom SD values
#'
#' @return Numeric vector of length length(x.knots) containing SD values for each component
#' @keywords internal
generate_sd_knots <- function(strategy, x.knots, sd.knot, sd.factor,
                            sd.random.range = c(0.5, 2), custom.sd.knots = NULL) {
    n.knots <- length(x.knots)

    switch(strategy,
        "fixed" = {
            rep(sd.knot, n.knots)
        },
        "adaptive" = {
            sd.vals <- numeric(n.knots)
            for (i in 1:n.knots) {
                left.dist <- if (i > 1) x.knots[i] - x.knots[i-1] else Inf
                right.dist <- if (i < n.knots) x.knots[i+1] - x.knots[i] else Inf
                sd.vals[i] <- sd.factor * min(left.dist, right.dist)
            }
            sd.vals
        },
        "stoch.adaptive" = {
            sd.vals <- numeric(n.knots)
            for (i in 1:n.knots) {
                left.dist <- if (i > 1) x.knots[i] - x.knots[i-1] else Inf
                right.dist <- if (i < n.knots) x.knots[i+1] - x.knots[i] else Inf
                base_sd <- sd.factor * min(left.dist, right.dist)
                sd.vals[i] <- runif(1,
                                  min = base_sd * sd.random.range[1],
                                  max = base_sd * sd.random.range[2])
            }
            sd.vals
        },
        "custom" = {
            if (length(custom.sd.knots) != n.knots) {
                stop("custom.sd.knots must have same length as number of components")
            }
            custom.sd.knots
        },
        stop("Invalid sd.strategy")
    )
}

#' Generate Amplitudes for Gaussian Mixture Components
#'
#' @param strategy Character string specifying the amplitude generation strategy. Must be one of:
#'   "mixed" (positive and negative), "positive" (only positive), or "custom" (user-specified)
#' @param n.components Integer specifying the number of Gaussian components
#' @param y.range Numeric vector of length 2 specifying the range \eqn{[\text{min}, \text{max}]} for y values
#' @param y.min.fraction Numeric between 0 and 1 specifying minimum absolute amplitude as fraction of range
#' @param custom.y.knots Optional numeric vector of custom amplitude values
#'
#' @return Numeric vector of length n.components containing amplitude values for each component
#' @keywords internal
generate_y_knots <- function(strategy, n.components, y.range, y.min.fraction,
                           custom.y.knots = NULL) {
    y.width <- diff(y.range)
    abs.y.min <- y.min.fraction * y.width
    abs.y.max <- abs(y.range[2])

    switch(strategy,
        "mixed" = {
            y.knots <- numeric(n.components)
            for (j in seq(n.components)) {
                if (runif(1) > 0.5) {
                    y.knots[j] <- runif(1, min = abs.y.min, max = abs.y.max)
                } else {
                    y.knots[j] <- runif(1, min = -abs.y.max, max = -abs.y.min)
                }
            }
            y.knots
        },
        "positive" = {
            runif(n.components, min = abs.y.min, max = abs.y.max)
        },
        "custom" = {
            if (length(custom.y.knots) != n.components) {
                stop("custom.y.knots must have same length as number of components")
            }
            custom.y.knots
        },
        stop("Invalid y.strategy")
    )
}

#' Generate Synthetic Data from a Mixture of Gaussians
#'
#' Creates synthetic data by generating a mixture of Gaussian components with
#' flexible positioning, scaling, and noise characteristics. This function provides
#' extensive control over the component placement, standard deviations, amplitudes,
#' and noise properties, making it suitable for testing and benchmarking smoothing
#' methods. Supports both absolute and sigma-relative positioning strategies.
#'
#' @param n.points Integer specifying the number of data points to generate
#' @param n.components Integer specifying the number of Gaussian components
#' @param x.range Numeric vector c(min, max) specifying the range for x values.
#'   Ignored when x.knots.strategy = "offset"
#' @param y.range Numeric vector c(min, max) specifying the range for y values
#'
#' @param x.knots.strategy Character string specifying how to position components.
#'   Must be one of:
#'   * "uniform": equally spaced components
#'   * "random": uniformly random positions
#'   * "jittered": perturbed uniform spacing
#'   * "offset": sequential spacing with sigma-relative gaps
#' @param jitter.fraction Numeric between 0 and 0.5 specifying maximum jitter as
#'   fraction of local spacing when x.knots.strategy = "jittered"
#' @param edge.buffer Numeric between 0 and 0.5 specifying minimum distance from
#'   range edges as fraction of total range (not used for "offset" strategy)
#'
#' @param sigma Numeric specifying the scale parameter when using offset strategy.
#'   Also informs the error distribution scale. Default is 1.5
#' @param xmin.factor Numeric specifying minimum spacing between components in
#'   units of sigma when x.knots.strategy = "offset". Default is 3
#' @param xmax.factor Numeric specifying maximum spacing between components in
#'   units of sigma when x.knots.strategy = "offset". Default is 7
#' @param x.offset.factor Numeric factor to extend x-range beyond components
#'   when x.knots.strategy = "offset". Default is 2
#'
#' @param sd.strategy Character string specifying how to generate component SDs.
#'   Must be one of:
#'   * "fixed": constant SD for all components
#'   * "adaptive": SD based on distance to neighbors
#'   * "stoch.adaptive": randomized adaptive SDs
#'   * "custom": user-specified SD values
#' @param sd.knot Numeric specifying the SD value when sd.strategy = "fixed"
#' @param sd.factor Numeric between 0 and 1 specifying the fraction of neighbor
#'   distance to use for SD in adaptive strategies
#' @param sd.random.range Numeric vector c(min, max) specifying multipliers for
#'   random SD generation in "stoch.adaptive" strategy
#' @param custom.sd.knots Optional numeric vector of length n.components
#'   specifying custom SD values
#'
#' @param y.strategy Character string specifying how to generate component amplitudes.
#'   Must be one of:
#'   * "mixed": mix of positive and negative amplitudes
#'   * "positive": only positive amplitudes
#'   * "custom": user-specified amplitudes
#' @param y.min.fraction Numeric between 0 and 1 specifying minimum absolute
#'   amplitude as fraction of y range
#' @param custom.y.knots Optional numeric vector of length n.components
#'   specifying custom amplitude values
#'
#' @param noise.fraction Numeric specifying noise standard deviation as fraction
#'   of y range
#' @param error.distribution Character string specifying noise distribution.
#'   Must be one of:
#'   * "norm": Gaussian noise
#'   * "laplace": Laplace noise
#'
#' @param seed Optional integer for random number generation
#' @param verbose Logical indicating whether to print generation details
#'
#' @section Positioning Strategies:
#' The function supports four strategies for positioning components along the x-axis:
#'
#' \describe{
#'   \item{\strong{uniform}}{Components are equally spaced within x.range with edge buffers}
#'   \item{\strong{random}}{Components are randomly positioned within x.range with edge buffers}
#'   \item{\strong{jittered}}{Components start on a uniform grid then are randomly perturbed.
#'     The jitter.fraction parameter controls the maximum displacement while maintaining order}
#'   \item{\strong{offset}}{Components are sequentially positioned with random gaps between
#'     xmin.factor*sigma and xmax.factor*sigma. The first component starts at 0, and the
#'     x-range is automatically determined}
#' }
#'
#' @section Jittering Details:
#' When using jittered positioning, components maintain their relative order:
#' \itemize{
#'   \item Interior components can move up to jitter.fraction times the distance
#'     to their neighbors
#'   \item Edge components have full freedom toward the range boundaries
#'   \item jitter.fraction must be \eqn{\le 0.5} to prevent crossing
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{x}{Numeric vector of x coordinates}
#'   \item{y}{Numeric vector of noisy y values}
#'   \item{y.true}{Numeric vector of true y values (without noise)}
#'   \item{x.knots}{Numeric vector of component x positions}
#'   \item{y.knots}{Numeric vector of component amplitudes}
#'   \item{sd.knots}{Numeric vector of component standard deviations}
#'   \item{parameters}{List of all input parameters used}
#' }
#'
#' @examples
#' # Basic usage with default parameters
#' data <- get.gaussian.mixture()
#'
#' # Custom mixture with specific parameters
#' data <- get.gaussian.mixture(
#'   n.points = 200,
#'   n.components = 3,
#'   x.knots.strategy = "jittered",
#'   sd.strategy = "stoch.adaptive",
#'   y.strategy = "mixed",
#'   noise.fraction = 0.1
#' )
#'
#' # Using sigma-scaled offset positioning
#' data <- get.gaussian.mixture(
#'   n.components = 4,
#'   x.knots.strategy = "offset",
#'   sigma = 2,
#'   xmin.factor = 2,
#'   xmax.factor = 5,
#'   sd.strategy = "adaptive"
#' )
#'
#' # Example with fixed standard deviations and custom amplitudes
#' data <- get.gaussian.mixture(
#'   n.components = 2,
#'   x.knots.strategy = "uniform",
#'   sd.strategy = "custom",
#'   custom.sd.knots = c(0.05, 0.15),
#'   y.strategy = "custom",
#'   custom.y.knots = c(1, -0.5)
#' )
#'
#' @export
get.gaussian.mixture <- function(
    # Core parameters
    n.points = 100,
    n.components = 5,
    x.range = c(0, 1),
    y.range = c(-1, 1),

    # Component positioning parameters
    x.knots.strategy = "uniform",
    jitter.fraction = 0.25,
    edge.buffer = 0.1,

    # Sigma-scaled positioning parameters (for offset strategy)
    sigma = 1.5,
    xmin.factor = 3,
    xmax.factor = 7,
    x.offset.factor = 2,

    # Component characteristics
    sd.strategy = "fixed",
    sd.knot = 0.1,
    sd.factor = 0.5,
    sd.random.range = c(0.5, 2),
    custom.sd.knots = NULL,

    y.strategy = "mixed",
    y.min.fraction = 0.1,
    custom.y.knots = NULL,

    # Noise parameters
    noise.fraction = 0.15,
    error.distribution = "norm",

    # Control parameters
    seed = NULL,
    verbose = FALSE) {

    # Parameter validation
    if (!is.null(seed)) set.seed(seed)

    # Validate numeric parameters
    if (!is.numeric(n.points) || n.points <= 0)
        stop("n.points must be positive")
    if (!is.numeric(n.components) || n.components <= 0)
        stop("n.components must be positive")
    if (!is.numeric(jitter.fraction) || jitter.fraction < 0 || jitter.fraction > 0.5)
        stop("jitter.fraction must be between 0 and 0.5")
    if (!is.numeric(edge.buffer) || edge.buffer < 0 || edge.buffer >= 0.5)
        stop("edge.buffer must be between 0 and 0.5")
    if (!is.numeric(y.min.fraction) || y.min.fraction <= 0 || y.min.fraction >= 1)
        stop("y.min.fraction must be between 0 and 1")
    if (!is.numeric(noise.fraction) || noise.fraction < 0)
        stop("noise.fraction must be non-negative")

    # Validate sigma-related parameters
    if (!is.numeric(sigma) || sigma <= 0)
        stop("sigma must be positive")
    if (!is.numeric(xmin.factor) || xmin.factor <= 0)
        stop("xmin.factor must be positive")
    if (!is.numeric(xmax.factor) || xmax.factor <= xmin.factor)
        stop("xmax.factor must be greater than xmin.factor")
    if (!is.numeric(x.offset.factor) || x.offset.factor <= 0)
        stop("x.offset.factor must be positive")

    # Validate ranges
    if (!is.numeric(x.range) || length(x.range) != 2)
        stop("x.range must be a numeric vector of length 2")
    if (x.range[2] <= x.range[1])
        stop("x.range[2] must be greater than x.range[1]")

    if (!is.numeric(y.range) || length(y.range) != 2)
        stop("y.range must be a numeric vector of length 2")
    if (y.range[2] <= y.range[1])
        stop("y.range[2] must be greater than y.range[1]")

    # Validate sd parameters
    if (!is.numeric(sd.knot) || sd.knot <= 0)
        stop("sd.knot must be positive")
    if (!is.numeric(sd.factor) || sd.factor <= 0 || sd.factor > 1)
        stop("sd.factor must be between 0 and 1")
    if (!is.numeric(sd.random.range) || length(sd.random.range) != 2 ||
        sd.random.range[1] <= 0 || sd.random.range[2] <= sd.random.range[1])
        stop("sd.random.range must be c(min, max) with 0 < min < max")

    # Validate strategies
    if (!x.knots.strategy %in% c("uniform", "random", "jittered", "offset"))
        stop("x.knots.strategy must be one of: 'uniform', 'random', 'jittered', 'offset'")
    if (!sd.strategy %in% c("fixed", "adaptive", "stoch.adaptive", "custom"))
        stop("sd.strategy must be one of: 'fixed', 'adaptive', 'stoch.adaptive', 'custom'")
    if (!y.strategy %in% c("mixed", "positive", "custom"))
        stop("y.strategy must be one of: 'mixed', 'positive', 'custom'")
    if (!error.distribution %in% c("norm", "laplace"))
        stop("error.distribution must be one of: 'norm', 'laplace'")

    # Validate custom inputs if specified
    if (sd.strategy == "custom") {
        if (is.null(custom.sd.knots))
            stop("custom.sd.knots must be provided when sd.strategy is 'custom'")
        if (length(custom.sd.knots) != n.components)
            stop("custom.sd.knots must have length equal to n.components")
        if (any(custom.sd.knots <= 0))
            stop("all custom.sd.knots must be positive")
    }

    if (y.strategy == "custom") {
        if (is.null(custom.y.knots))
            stop("custom.y.knots must be provided when y.strategy is 'custom'")
        if (length(custom.y.knots) != n.components)
            stop("custom.y.knots must have length equal to n.components")
    }

    # Generate component x positions
    if (x.knots.strategy == "offset") {
        # Sigma-scaled offset positioning
        x.knots <- numeric(n.components)
        x.knots[1] <- 0  # First knot at origin
        for (j in 2:n.components) {
            x.knots[j] <- runif(1,
                               min = x.knots[j - 1] + xmin.factor * sigma,
                               max = x.knots[j - 1] + xmax.factor * sigma)
        }
        # Update x.range based on generated knots
        x.range <- c(min(x.knots) - x.offset.factor * sigma,
                     max(x.knots) + x.offset.factor * sigma)
    } else {
        # Use existing strategies
        x.knots <- generate_x_knots(x.knots.strategy, n.components, x.range,
                                   edge.buffer, jitter.fraction)
    }

    # Generate x grid
    x <- seq(x.range[1], x.range[2], length.out = n.points)

    # Generate component amplitudes
    y.knots <- generate_y_knots(y.strategy, n.components, y.range,
                               y.min.fraction, custom.y.knots)

    # Generate standard deviations
    sd.vals <- generate_sd_knots(sd.strategy, x.knots, sd.knot, sd.factor,
                                sd.random.range, custom.sd.knots)

    # Generate Gaussian mixture
    y.true <- numeric(n.points)
    for (i in 1:n.points) {
        for (j in 1:n.components) {
            y.true[i] <- y.true[i] + y.knots[j] *
                exp(-((x[i] - x.knots[j])^2) / (2 * sd.vals[j]^2))
        }
    }

    # Add noise
    y.width <- diff(y.range)
    noise.sd <- noise.fraction * y.width

    y <- if (error.distribution == "norm") {
        y.true + rnorm(n.points, 0, noise.sd)
    } else {
        y.true + rlaplace(n.points, location = 0, scale = noise.sd)
    }

    # Prepare return structure
    result <- list(
        x = x,
        y = y,
        y.true = y.true,
        x.knots = x.knots,
        y.knots = y.knots,
        sd.knots = sd.vals,
        parameters = list(
            n.points = n.points,
            n.components = n.components,
            x.range = x.range,
            y.range = y.range,
            x.knots.strategy = x.knots.strategy,
            jitter.fraction = jitter.fraction,
            edge.buffer = edge.buffer,
            sigma = sigma,
            xmin.factor = xmin.factor,
            xmax.factor = xmax.factor,
            x.offset.factor = x.offset.factor,
            sd.strategy = sd.strategy,
            sd.knot = sd.knot,
            sd.factor = sd.factor,
            sd.random.range = sd.random.range,
            y.strategy = y.strategy,
            y.min.fraction = y.min.fraction,
            noise.fraction = noise.fraction,
            error.distribution = error.distribution,
            seed = seed
        )
    )

    if (verbose) {
        cat("Generated Gaussian mixture with:\n")
        cat("  Components:", n.components, "\n")
        cat("  X strategy:", x.knots.strategy, "\n")
        if (x.knots.strategy == "offset") {
            cat("  Sigma:", sigma, "\n")
            cat("  X range:", round(x.range[1], 3), "to", round(x.range[2], 3), "\n")
        }
        cat("  SD strategy:", sd.strategy, "\n")
        cat("  Y strategy:", y.strategy, "\n")
        cat("  Noise distribution:", error.distribution, "\n")
    }

    return(result)
}

#' Generate a 1D Gaussian Mixture with Direct Knot Specification
#'
#' A convenience function that generates a 1D Gaussian mixture by directly specifying
#' component positions and amplitudes. This function provides a simpler interface to
#' create Gaussian mixtures when exact component locations are known.
#'
#' @param n.points Integer. The number of points in the output grid. Default is 100.
#' @param x.knots Numeric vector. The x-coordinates of the Gaussian components' centers.
#'   Default is c(-5, 0, 10).
#' @param y.knots Numeric vector. The amplitudes of the Gaussian components.
#'   Must have the same length as x.knots. Default is c(5, 8, 2.5).
#' @param sd.knot Numeric. The standard deviation for all Gaussian components.
#'   Default is 1.0.
#' @param x.offset Numeric. The range offset to extend beyond the minimum and
#'   maximum x.knots. Default is 3.
#' @param add.noise Logical. Whether to add noise to the generated mixture.
#'   Default is FALSE.
#' @param noise.fraction Numeric. If add.noise is TRUE, the noise standard deviation
#'   as a fraction of the y range. Default is 0.1.
#' @param out.dir Character string or NULL. If provided, specifies the directory
#'   to save the results. Default is NULL (results not saved).
#' @param verbose Logical. Whether to print the output file path if saving.
#'   Default is FALSE.
#'
#' @return A list containing:
#'   \item{x}{Numeric vector. The x-coordinates of the generated points.}
#'   \item{y}{Numeric vector. The y-coordinates of the generated Gaussian mixture.}
#'   \item{y.true}{Numeric vector. The true y-coordinates without noise (if noise added).}
#'   \item{x.knots}{Numeric vector. The x-coordinates of the components.}
#'   \item{y.knots}{Numeric vector. The amplitudes of the components.}
#'   \item{sd.knots}{Numeric vector. The standard deviations of the components.}
#'   \item{file.name}{Character string. The file name if results were saved, NULL otherwise.}
#'   \item{params}{List. All input parameters used to generate the mixture.}
#'
#' @details
#' This function provides a simplified interface for generating Gaussian mixtures
#' when the user wants to specify exact component locations and amplitudes. Since
#' \code{get.gaussian.mixture} doesn't support custom x positions directly, this
#' function generates the mixture by evaluating the Gaussian components at the
#' specified positions. It's particularly useful for creating test cases or
#' reproducing specific mixture configurations.
#'
#' Note: Unlike \code{get.gaussian.mixture}, this function directly computes the
#' mixture values rather than using positioning strategies.
#'
#' @seealso \code{\link{get.gaussian.mixture}} for more flexible mixture generation options
#'
#' @examples
#' # Generate mixture with default parameters
#' result1 <- generate.1d.gaussian.mixture()
#' plot(result1$x, result1$y, type = 'l', main = "Default 1D Gaussian Mixture")
#'
#' # Generate mixture with custom parameters
#' result2 <- generate.1d.gaussian.mixture(
#'   n.points = 200,
#'   x.knots = c(-10, 0, 5, 15),
#'   y.knots = c(2, 5, 8, 3),
#'   sd.knot = 0.8,
#'   x.offset = 5
#' )
#'
#' # Generate mixture with noise
#' result3 <- generate.1d.gaussian.mixture(
#'   x.knots = c(-2, 2),
#'   y.knots = c(1, -1),
#'   add.noise = TRUE,
#'   noise.fraction = 0.15
#' )
#'
#' @export
generate.1d.gaussian.mixture <- function(n.points = 100,
                                         x.knots = c(-5, 0, 10),
                                         y.knots = c(5, 8, 2.5),
                                         sd.knot = 1.0,
                                         x.offset = 3,
                                         add.noise = FALSE,
                                         noise.fraction = 0.1,
                                         out.dir = NULL,
                                         verbose = FALSE) {

    ## Input validation
    if (length(x.knots) != length(y.knots)) {
        stop("x.knots and y.knots must have the same length")
    }
    if (n.points <= 0 || sd.knot <= 0 || x.offset <= 0) {
        stop("n.points, sd.knot, and x.offset must be positive")
    }
    if (!is.null(out.dir) && !dir.exists(out.dir)) {
        stop("Specified output directory does not exist")
    }

    ## Calculate x.range from knots and offset
    x.range <- c(min(x.knots) - x.offset, max(x.knots) + x.offset)

    ## Generate x grid
    x <- seq(x.range[1], x.range[2], length.out = n.points)

    ## Calculate y.range to accommodate the specified y.knots
    y.margin <- 0.1 * max(abs(y.knots))  # 10% margin
    y.range <- c(min(c(y.knots, 0)) - y.margin, max(c(y.knots, 0)) + y.margin)

    ## Since get.gaussian.mixture doesn't support custom x.knots,
    ## we'll compute the mixture directly
    n.components <- length(x.knots)
    y.true <- numeric(n.points)

    ## Generate Gaussian mixture by direct computation
    for (i in 1:n.points) {
        for (j in 1:n.components) {
            y.true[i] <- y.true[i] + y.knots[j] *
                exp(-((x[i] - x.knots[j])^2) / (2 * sd.knot^2))
        }
    }

    ## Add noise if requested
    if (add.noise) {
        y.width <- diff(y.range)
        noise.sd <- noise.fraction * y.width
        y <- y.true + rnorm(n.points, 0, noise.sd)
    } else {
        y <- y.true
    }

    ## Prepare results in the expected format
    results <- list(
        x = x,
        y = y,
        y.true = y.true,
        x.knots = x.knots,
        y.knots = y.knots,
        sd.knots = rep(sd.knot, n.components),
        file.name = NULL,
        params = list(
            n.points = n.points,
            x.knots = x.knots,
            y.knots = y.knots,
            sd.knot = sd.knot,
            x.offset = x.offset,
            add.noise = add.noise,
            noise.fraction = noise.fraction,
            out.dir = out.dir
        )
    )

    ## Save results if out.dir is specified
    if (!is.null(out.dir)) {
        ## Generate a unique identifier
        unique.id <- paste(
            paste(round(x.knots, 2), collapse = "_"),
            paste(round(y.knots, 2), collapse = "_"),
            sd.knot,
            sep = "_"
        )

        file.name <- sprintf("gaussian_mixture_1d_n%d_k%d_sd%.1f_%s.rda",
                             n.points,
                             length(x.knots),
                             sd.knot,
                             unique.id)

        results$file.name <- file.name

        file.path <- file.path(out.dir, file.name)
        save(results, file = file.path)
        if (verbose) {
            cat("Results saved to:", file.path, "\n")
        }
    }

    return(results)
}

#' Generate Binary Sample from Gaussian Mixture
#'
#' Transforms a Gaussian mixture into a probability function using logistic transformation
#' and generates binary samples based on these probabilities. This is useful for creating
#' synthetic binary response data with complex non-linear relationships.
#'
#' @param gaussian.mixture.result List or character. Either the result from
#'   \code{generate.1d.gaussian.mixture} or \code{get.gaussian.mixture}, or a path
#'   to a file containing such results.
#' @param n.sample.points Integer. Number of points to sample. Default is 50.
#' @param sample.method Character. Method for selecting sample points:
#'   \itemize{
#'     \item "uniform": Random uniform sampling within x range (default)
#'     \item "regular": Equally spaced points
#'     \item "custom": User-specified x locations (requires x.sample parameter)
#'   }
#' @param x.sample Numeric vector. Custom x locations when sample.method = "custom".
#' @param use.true.y Logical. If TRUE, use y.true (without noise) if available.
#'   Default is TRUE.
#' @param transform.method Character. Method to transform y values to probabilities:
#'   \itemize{
#'     \item "logit": Normalize to logit.range then apply inverse logit (default)
#'     \item "direct": Normalize directly to \eqn{[0,1]}
#'   }
#' @param logit.range Numeric vector of length 2. When transform.method = "logit",
#'   the range to normalize y values before applying inverse logit. Default is c(-3, 3).
#' @param out.dir Character or NULL. Directory to save results. If NULL, results
#'   are not saved. Default is NULL.
#' @param verbose Logical. Whether to print file path when saving. Default is FALSE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing:
#'   \item{x}{Numeric vector. The x-coordinates of the sampled points.}
#'   \item{y.binary}{Integer vector. The binary sample (0 or 1).}
#'   \item{y.prob}{Numeric vector. The probability of Y=1 for each point.}
#'   \item{y.smooth}{Numeric vector. The interpolated y values before transformation.}
#'   \item{gaussian.mixture}{List. The original Gaussian mixture result.}
#'   \item{params}{List. Parameters used for generation.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the smooth function from the Gaussian mixture
#'   \item Selects sample points according to sample.method
#'   \item Interpolates y values at sample points
#'   \item Transforms y values to probabilities using specified method
#'   \item Generates binary outcomes using binomial sampling
#' }
#'
#' For the logit transformation method, y values are first normalized to logit.range,
#' then the inverse logit function is applied: p = 1/(1 + exp(-y_normalized)).
#' This creates a smooth probability function that approaches 0 and 1 asymptotically.
#'
#' @examples
#' # First generate a Gaussian mixture
#' gm <- generate.1d.gaussian.mixture(
#'   x.knots = c(-2, 0, 2),
#'   y.knots = c(-1, 2, -1.5)
#' )
#'
#' # Generate binary sample with default settings
#' binary.data <- generate.binary.sample.from.gaussian.mixture(gm)
#'
#' # Plot the results
#' plot(gm$x, gm$y, type = 'l', col = 'blue',
#'      main = "Gaussian Mixture and Binary Sample")
#' points(binary.data$x, binary.data$y.binary, pch = 19,
#'        col = ifelse(binary.data$y.binary == 1, "green", "red"))
#' lines(binary.data$x, binary.data$y.prob, col = 'orange', lwd = 2)
#'
#' # Generate regular grid sample with wider logit range
#' binary.regular <- generate.binary.sample.from.gaussian.mixture(
#'   gm,
#'   n.sample.points = 100,
#'   sample.method = "regular",
#'   logit.range = c(-5, 5)  # Steeper probability transitions
#' )
#'
#' @importFrom stats approx rbinom runif
#' @export
generate.binary.sample.from.gaussian.mixture <- function(gaussian.mixture.result,
                                                        n.sample.points = 50,
                                                        sample.method = c("uniform", "regular", "custom"),
                                                        x.sample = NULL,
                                                        use.true.y = TRUE,
                                                        transform.method = c("logit", "direct"),
                                                        logit.range = c(-3, 3),
                                                        out.dir = NULL,
                                                        verbose = FALSE,
                                                        seed = NULL) {

    # Set seed if provided
    if (!is.null(seed)) set.seed(seed)

    # Match arguments
    sample.method <- match.arg(sample.method)
    transform.method <- match.arg(transform.method)

    # Validate inputs
    if (length(logit.range) != 2 || logit.range[1] >= logit.range[2]) {
        stop("logit.range must be a numeric vector of length 2 with min < max")
    }

    if (sample.method == "custom" && is.null(x.sample)) {
        stop("x.sample must be provided when sample.method = 'custom'")
    }

    if (!is.null(out.dir) && !dir.exists(out.dir)) {
        stop("Specified output directory does not exist")
    }

    # Load results if a file path is provided
    if (is.character(gaussian.mixture.result)) {
        if (!file.exists(gaussian.mixture.result)) {
            stop("File not found: ", gaussian.mixture.result)
        }
        # Load into a new environment to avoid conflicts
        env <- new.env()
        load(gaussian.mixture.result, envir = env)
        # Try to find the results object
        obj.names <- ls(env)
        if ("results" %in% obj.names) {
            gaussian.mixture.result <- env$results
        } else if (length(obj.names) == 1) {
            gaussian.mixture.result <- get(obj.names[1], envir = env)
        } else {
            stop("Cannot determine which object to use from loaded file")
        }
    }

    # Validate gaussian.mixture.result structure
    required.fields <- c("x", if(use.true.y && "y.true" %in% names(gaussian.mixture.result)) "y.true" else "y")
    if (!is.list(gaussian.mixture.result) || !all(required.fields %in% names(gaussian.mixture.result))) {
        stop("gaussian.mixture.result must be a list containing 'x' and 'y' (or 'y.true') components")
    }

    # Extract x and y from the Gaussian mixture result
    x.gaussian <- gaussian.mixture.result$x
    y.gaussian <- if (use.true.y && "y.true" %in% names(gaussian.mixture.result)) {
        gaussian.mixture.result$y.true
    } else {
        gaussian.mixture.result$y
    }

    # Generate sample points based on method
    if (sample.method == "uniform") {
        x.sample <- sort(runif(n.sample.points, min = min(x.gaussian), max = max(x.gaussian)))
    } else if (sample.method == "regular") {
        x.sample <- seq(min(x.gaussian), max(x.gaussian), length.out = n.sample.points)
    } else if (sample.method == "custom") {
        if (any(x.sample < min(x.gaussian) | x.sample > max(x.gaussian))) {
            warning("Some x.sample values are outside the range of x.gaussian")
        }
        x.sample <- sort(x.sample)
        n.sample.points <- length(x.sample)
    }

    # Interpolate y values for the sample points
    y.smooth <- approx(x.gaussian, y.gaussian, xout = x.sample, rule = 2)$y

    # Transform to probabilities
    if (transform.method == "logit") {
        # Normalize to logit.range then apply inverse logit
        y.normalized <- (y.smooth - min(y.smooth)) / (max(y.smooth) - min(y.smooth))
        y.normalized <- y.normalized * (logit.range[2] - logit.range[1]) + logit.range[1]
        y.prob <- 1 / (1 + exp(-y.normalized))
    } else if (transform.method == "direct") {
        # Normalize directly to [0, 1]
        y.prob <- (y.smooth - min(y.smooth)) / (max(y.smooth) - min(y.smooth))
    }

    # Generate binary sample
    y.binary <- rbinom(n.sample.points, size = 1, prob = y.prob)

    # Prepare results
    results <- list(
        x = x.sample,
        y.binary = y.binary,
        y.prob = y.prob,
        y.smooth = y.smooth,
        gaussian.mixture = gaussian.mixture.result,
        params = list(
            n.sample.points = n.sample.points,
            sample.method = sample.method,
            use.true.y = use.true.y,
            transform.method = transform.method,
            logit.range = logit.range,
            seed = seed
        )
    )

    # Save results if out.dir is specified
    if (!is.null(out.dir)) {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        file.name <- sprintf("binary_sample_n%d_%s_%s.rda",
                            n.sample.points,
                            sample.method,
                            timestamp)
        file.path <- file.path(out.dir, file.name)
        save(results, file = file.path)
        results[["file.name"]] <- file.name
        results[["file.path"]] <- file.path
        if (verbose) {
            cat("Results saved to:", file.path, "\n")
        }
    }

    return(results)
}

#' Compute Value of a Single Skewed Gaussian
#'
#' This function calculates the value of a skewed Gaussian distribution at a given point.
#'
#' @description
#' The skewed Gaussian distribution is defined by the following formula:
#'
#' f(x | xi, omega, alpha) = (2 / omega) * phi((x - xi) / omega) * Phi(alpha * ((x - xi) / omega))
#'
#' where phi is the standard normal density and Phi is the standard normal
#' cumulative distribution function.
#'
#' @details
#' The parameters of the distribution are:
#' * xi (equivalent to mu): location parameter
#' * omega (equivalent to sigma): scale parameter
#' * alpha: shape parameter controlling the skewness
#'
#' @param x Numeric. The point at which to evaluate the skewed Gaussian.
#' @param mu Numeric. The location parameter (equivalent to xi in the formula).
#' @param sigma Numeric. The scale parameter (equivalent to omega in the formula).
#' @param alpha Numeric. The shape parameter controlling the skewness.
#'
#' @return Numeric. The value of the skewed Gaussian distribution at point x.
#'
#' @examples
#' skewed.gaussian(0, mu = 0, sigma = 1, alpha = 2)
#'
#' @export
skewed.gaussian <- function(x, mu, sigma, alpha) {
    ## Compute the standard normal density (phi)
    phi_val <- dnorm(x, mean = mu, sd = sigma)
    ## Compute the standard normal cumulative distribution (Phi)
    Phi_val <- pnorm(alpha * (x - mu) / sigma)
    ## Compute the skewed Gaussian density
    skewed_val <- (2 / sigma) * phi_val * Phi_val

    return(skewed_val)
}

#' Generate a Spline Function
#'
#' This function creates a spline function based on a given set of points. It interpolates these points
#' using cubic splines, providing a smooth curve through the data.
#'
#' @param X A numeric matrix containing the coordinates of the points.
#' @param y A numeric vector containing the y-coordinates of the points. If not provided,
#'          random values will be generated.
#'
#' @return A spline function that interpolates the given points. The function takes a numeric
#'         value (x-coordinate) as input and returns the interpolated y-coordinate.
#'
#' @details The function first checks that the lengths of `x` and `y` are equal and then sorts
#'          these vectors based on the `x` values. It then creates a cubic spline function
#'          using the `splinefun` function in R, with the method set to "fmm" for flexible
#'          monotone splines. If `y` is not provided, random values are generated to create
#'          the spline.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' x <- seq(0, 10, length.out = 10)
#' y <- sin(x)  # Example y values
#' spline.fn <- synthetic.spline(x, y)
#' point <- 5  # An x-coordinate
#' value <- spline.fn(point)  # Evaluate the spline function at this point
#' }
#' @export
synthetic.xD.spline <- function(X, y = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  ## Check if y is provided, if not generate random y
  if (is.null(y)) {
    y <- runif(n)
  }

  ## Check if length of y matches the number of rows in X
  if (length(y) != n) {
    stop("Length of y must match the number of rows in X.")
  }

  ## List to store spline functions for each dimension
  spline.fns <- vector("list", p)

  ## Generate spline function for each dimension
  for (i in 1:p) {
    ## Extract the ith column of X
    x.i <- X[, i]

    ## Sort x.i and y based on x.i values
    order.indices <- order(x.i)
    x.i.sorted <- x.i[order.indices]
    y.sorted <- y[order.indices]

    ## Create spline function for this dimension
    spline.fns[[i]] <- splinefun(x.i.sorted, y.sorted, method = "fmm")
  }

  ## Returning a list of spline functions
  return(spline.fns)
}



#' Generate a Synthetic Smooth Function in One Dimension
#'
#' This function constructs a synthetic smooth function using a series of Gaussian distributions
#' in one-dimensional space. It is centered at given points with standard deviations determined
#' by the distances to neighboring points. The function can handle both positive and negative
#' Gaussian distributions.
#'
#' @param x       Values at which the function is to be evaluated.
#' @param x.knot A numeric vector of the centers of the Gaussian functions.
#' @param y.knot An optional numeric vector of response values at the locations specified in `x.knot`.
#'          If `y.knot` is not provided, random values will be generated.
#' @param sd.knot Standard deviations of the Gaussians.
#'
#' @return A function representing the synthetic smooth function. This function takes a numeric
#'         value as input and returns a numeric value, representing the function's output at that point.
#'
#' @details For each point in `x`, the function calculates the standard
#'     deviation of the Gaussian distribution as half of the minimum distance to
#'     its neighbors. It then forms a synthetic function by combining these
#'     Gaussian distributions. The function can generate smooth functions that
#'     take also negative values by multiplying Gaussian components by a vector
#'     of ±1 values.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' x <- seq(0, 1, length.out = 10)
#' synthetic.fn <- generate.synthetic.function(x)
#' point <- 0.5 # A point in 1D space
#' value <- synthetic.fn(point) # Evaluate the function at the given point
#' }
#' @export
synthetic.mixture.of.gaussians <- function(x, x.knot, y.knot = NULL, sd.knot = NULL) {

    n.knots <- length(x.knot)

    ## Check if y is provided, if not generate random y
    if (is.null(y.knot)) {
        y.knot <- runif(n)
    }

    ## Sort x.knot and y.knot based on x.knot values
    order.indices <- order(x.knot)
    x.knot <- x.knot[order.indices]
    y.knot <- y.knot[order.indices]

    if (!is.null(sd.knot)) {
        if (length(sd.knot) == length(y.knot)) {
            sd.vals <- sd.knot
        } else {
            sd.vals <- rep(sd.knot, length(y.knot))
        }
    } else {
        ## Calculate standard deviations for Gaussian functions
        sd.vals <- numeric(n)
        for (i in 1:n.knots) {
            left.dist <- if (i > 1) x.knot[i] - x.knot[i - 1] else Inf
            right.dist <- if (i < n.knots) x.knot[i + 1] - x.knot[i] else Inf
            sd.vals[i] <- min(left.dist, right.dist) / 2
        }
    }

    ## Define the synthetic function
    synthetic.function <- function(x.val) {
        result <- 0
        for (i in 1:n.knots) {
            result <- result + y.knot[i] * exp(-((x.val - x.knot[i])^2) / (2 * sd.vals[i]^2))
        }
        return(result)
    }

    return(synthetic.function(x))
}

#' Generate a Synthetic Smooth Function in Higher Dimensions
#'
#' This function creates a synthetic smooth function based on a set of points in a multi-dimensional space.
#' It utilizes Gaussian functions centered at each point with standard deviations determined by the minimum
#' distance to other points. The function allows for both positive and negative Gaussians.
#'
#' @param X A numeric matrix where each row represents a point in a multi-dimensional space (R^d).
#' @param y An optional numeric vector of response values corresponding to each row in `X`.
#'          If not provided, random values are generated.
#'
#' @return A function that represents the synthetic smooth function.
#'         This function takes a numeric vector (a point in R^d) as input and returns a numeric value.
#'
#' @details The function calculates the standard deviation for each Gaussian centered at a point in `X`
#'          as half the minimum distance to any other point in `X`. It generates a random vector of +1 or -1
#'          to multiply with the Gaussians, allowing the inclusion of negative values. The output is a function
#'          that, when evaluated at any given point in R^d, provides the corresponding value of
#'          the synthetic smooth function.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' X <- matrix(runif(20), ncol = 2) # Random 2D points
#' synthetic_fn <- generate_synthetic_function_higher_dim(X)
#' point <- c(0.5, 0.5) # A point in 2D space
#' value <- synthetic_fn(point) # Evaluate the function at the given point
#' }
#' @export
generate.synthetic.function.higher.dim <- function(X, y = NULL) {
  n <- nrow(X)

  ## Check if y is provided, if not generate random y
  if (is.null(y)) {
    y <- runif(n)
  }

  ## Calculate standard deviations for Gaussian functions
  sd.vals <- numeric(n)
  for (i in 1:n) {
    distances <- apply(X, 1, function(point) sqrt(sum((point - X[i,])^2)))
    distances <- distances[distances > 0]  ## exclude distance from point to itself
    sd.vals[i] <- if (length(distances) > 0) min(distances) / 2 else 0
  }

  ## Generate random signs
  signs <- sample(c(-1, 1), n, replace = TRUE)

  ## Define the synthetic function
  synthetic.function <- function(x.val) {
    result <- 0
    for (i in 1:n) {
      gauss.val <- exp(-sum((x.val - X[i,])^2) / (2 * sd.vals[i]^2))
      result <- result + y[i] * signs[i] * gauss.val
    }
    return(result)
  }

  return(synthetic.function)
}



#' Removes Close Neighbors
#'
#' This function takes a vector of positions (`x`) and a minimum distance
#' (`min.dist`), and returns a modified vector of knots where no two knots are
#' closer than the specified minimum distance.
#'
#' The function first sorts the input vector of knots. It then iterates through
#' the sorted knots, keeping only those that are at least `min.dist` away from
#' the previously accepted knot. This ensures that the output vector of knots
#' satisfies the distance requirement.
#'
#'
#' @param x         A numeric vector.
#' @param min.dist A positive numeric value specifying the minimum allowable
#'     distance between any two consecutive elements of x.
#'
#' @return A modified version of x such that no two consecutive elements are closer than `min.dist`.
#'
#' @examples
#' \dontrun{
#' x <- c(0.1, 0.2, 0.4, 0.7, 1.0)
#' min.dist <- 0.25
#' remove.close.neighbors(x, min.dist)
#' }
#' @export
remove.close.neighbors <- function(x, min.dist) {
  if (min.dist <= 0) {
    stop("min.dist must be positive")
  }

  # Sort the knots
  x.sorted <- sort(x)

  # Initialize the vector for the result
  x.filtered <- x.sorted[1]

  for (knot in x.sorted[-1]) {
    if (min(abs(x.filtered - knot)) >= min.dist) {
      x.filtered <- c(x.filtered, knot)
    }
  }

  return(x.filtered)
}

#' Bump Function
#'
#' This function creates a bump function, which is a smooth function with compact support.
#' The bump function is commonly used in mathematical analysis, particularly in the construction
#' of partitions of unity. The function reaches its peak at the `offset` and smoothly decreases
#' to zero at the boundary of its support.
#'
#' @param x A numeric vector or a single numeric value where the bump function is evaluated.
#' @param offset The center (offset) of the bump function (defaults to 0).
#'               The function reaches its maximum value at this point.
#' @param h The radius of the support of the bump function (defaults to 1).
#'          The function is zero outside the interval \eqn{[offset - h, offset + h]}.
#' @param q The power to which the Gaussian function is raised (defaults to 4).
#'
#' @return A numeric vector or a single numeric value representing the value(s) of the bump function
#'         at the input `x`. The values are in the range \eqn{[0, 1]}, with the function smoothly approaching
#'         zero at the boundaries of its support.
#'
#' @details The bump function is defined as exp(-1 / (h - (x - offset)^2)) for |x - offset| < h, and 0 otherwise.
#'          It is infinitely differentiable and is used in situations where a smooth function with compact
#'          support is needed.
#'
#' @examples
#' \dontrun{
#' x <- seq(-2, 2, length.out = 100)
#' plot(x, bump.fn(x), type = "l")
#' }
#' @export
bump.fn <- function(x, offset = 0, h = 1, q = 4) {
    n <- length(x)
    h2 <- h^2

    if ( n == 1 ) {
        ret <- 0
        if ( abs(x - offset) < h ) {
            ret <- exp( -1 / (h2 - (x - offset)^2) )
            ret <- ret^q
        }
    } else {
        ret <- numeric(n)
        for ( i in seq(n) ) {
            if ( abs(x[i] - offset) < h ) {
                ret[i] <- exp( -1 / (h2 - (x[i] - offset)^2) )
                ret[i] <- ret[i]^q
            }
        }
    }

    return(ret)
}

#' Generalized Gaussian Function
#'
#' Computes the values of a generalized Gaussian function, which allows for adjusting the shape of the curve
#' by altering the power of the absolute difference term.
#'
#' @param x A numeric vector or a single numeric value where the generalized Gaussian function is evaluated.
#' @param offset The center of the Gaussian function (defaults to 0).
#' @param h The scale parameter of the Gaussian function (defaults to 1).
#' @param p The power to which the absolute difference is raised (defaults to 2).
#'
#' @return A numeric vector or a single numeric value representing the value(s) of the generalized Gaussian
#'         function at the input `x`. The function's shape changes depending on the value of `p`.
#'
#' @details The generalized Gaussian function is defined as exp(-((abs(x - offset))^p) / h^p).
#'          The parameter `p` allows for adjusting the shape of the Gaussian curve:
#'          a value of `p` = 2 gives the standard Gaussian, `p` < 2 gives a curve with heavier tails,
#'          and `p` > 2 gives a curve with lighter tails. This function can be useful in various
#'          statistical applications where different shapes of Gaussian-like curves are required.
#'
#' @examples
#' \dontrun{
#' x <- seq(-5, 5, length.out = 100)
#' plot(x, ggaussian(x, p = 2), type = "l")  # Standard Gaussian
#' plot(x, ggaussian(x, p = 1), type = "l")  # Laplace distribution (heavier tails)
#' plot(x, ggaussian(x, p = 3), type = "l")  # Lighter tails than Gaussian
#' }
#' @export
ggaussian <- function(x, offset = 0, h = 1, p = 2) {
    n <- length(x)
    h2 <- h^p

    if ( n == 1 ) {
        ret <- 0
        ret <- exp( -(abs(x - offset))^p / h2 )
    } else {
        ret <- numeric(n)
        for ( i in seq(n) ) {
            ret[i] <- exp( -(abs(x[i] - offset))^p / h2 )
        }
    }

    return(ret)
}


#' q-Exponential Gaussian Function
#'
#' This function computes the values of a modified Gaussian function, known as a q-Gaussian,
#' which raises the standard Gaussian to the power of `q`. This modification can be used to
#' adjust the "flatness" or "peakiness" of the Gaussian curve.
#'
#' @param x A numeric vector or a single numeric value where the q-Gaussian function is evaluated.
#' @param offset The center of the q-Gaussian function (defaults to 0).
#' @param h The standard deviation of the Gaussian part of the function (defaults to 1).
#' @param q The power to which the Gaussian function is raised (defaults to 4).
#'
#' @return A numeric vector or a single numeric value representing the value(s) of the q-Gaussian function
#'         at the input `x`. The function alters the shape of a standard Gaussian curve based on the value of `q`.
#'
#' @details The q-Gaussian function is defined as (exp(-((x - offset)^2) / h2))^q.
#'          The parameter `q` allows for adjusting the shape of the Gaussian curve:
#'          higher values of `q` make the curve peakier, while lower values make it flatter.
#'          This function can be useful in various statistical and signal processing applications
#'          where a modified Gaussian shape is required.
#'
#' @examples
#' \dontrun{
#' x <- seq(-5, 5, length.out = 100)
#' plot(x, qgaussian(x, q = 4), type = "l")
#' }
#' @export
qexp.gaussian <- function(x, offset = 0, h = 1, q = 4) {
    n <- length(x)
    h2 <- h^2

    if ( n == 1 ) {
        ret <- 0
        ret <- exp( -q * (x - offset)^2 / h2 )
    } else {
        ret <- numeric(n)
        for ( i in seq(n) ) {
            ret[i] <- exp( -q * (x[i] - offset)^2 / h2 )
        }
    }

    return(ret)
}

#' Creates Left Asymmetric Bump Function
#'
#' This function creates an asymmetric bump function, which is a modification of the standard bump function
#' with asymmetric behavior to the left of the offset. It is smooth with compact support, primarily used
#' in mathematical analysis and applications requiring non-symmetric smooth functions.
#'
#' @param x A numeric vector or a single numeric value where the bump function is evaluated.
#' @param offset The center (offset) of the bump function (defaults to 0).
#'               The function behaves differently to the left and right of this point.
#' @param h The radius of the support of the bump function (defaults to 1).
#'          The function is zero outside the interval \eqn{[offset - h, offset + h]}.
#' @param q The power to which the Gaussian function to the right of the offset is raised (defaults to 1).
#'
#' @return A numeric vector or a single numeric value representing the value(s) of the asymmetric bump function
#'         at the input `x`. The function smoothly approaches zero at the boundaries of its support and has
#'         different behavior to the left of the offset.
#'
#' @details The asymmetric bump function is defined differently to the left and right of the offset.
#'          For `x < offset` and `offset - x < h`, it is defined as exp(-1 / h - (x - offset)^2).
#'          For `|x - offset| < h`, it follows the standard bump function exp(-1 / (h - (x - offset)^2)).
#'          It is zero outside the specified interval, maintaining smooth transitions and compact support.
#'
#' @examples
#' \dontrun{
#' x <- seq(-2, 2, length.out = 100)
#' plot(x, AsymmetricBumpFunction(x), type = "l")
#' }
#' @export
left.asymmetric.bump.fn <- function(x, offset = 0, h = 1, q = 1) {

    n <- length(x)
    h2 <- h^2

    if ( n == 1 ) {
        ret <- 0
        if ( x < offset ) {
            ret <- exp( -1 / h2 - (x - offset)^2 )
        } else if ( abs(x - offset) < h ) {
            ret <- exp( -1 / (h2 - (x - offset)^2) )
            ret <- ret^q
        }
    } else {
        ret <- numeric(n)
        for ( i in seq(n) ) {
            if ( x[i] < offset ) {
                ret[i] <- exp( -1 / h2 - (x[i] - offset)^2 )
            } else if ( x[i] - offset < h ) {
                ret[i] <- exp( -1 / (h2 - (x[i] - offset)^2) )
                ret[i] <- ret[i]^q
            }
        }
    }

    return(ret)
}

#' Creates Right Asymmetric Bump Function
#'
#' This function creates an asymmetric bump function, which is a modification of the standard bump function
#' with asymmetric behavior to the right of the offset. It is smooth with compact support, primarily used
#' in mathematical analysis and applications requiring non-symmetric smooth functions.
#'
#' @param x      A numeric vector or a single numeric value at which the bump function is evaluated.
#' @param offset The center (offset) of the bump function (defaults to 0).
#'               The function behaves differently to the right and right of this point.
#' @param h The radius of the support of the bump function (defaults to 1).
#'          The function is zero outside the interval \eqn{[offset - h, offset + h]}.
#' @param q The power to which the Gaussian function to the left of the offset is raised (defaults to 1).
#'
#' @return A numeric vector or a single numeric value representing the value(s) of the asymmetric bump function
#'         at the input `x`. The function smoothly approaches zero at the boundaries of its support and has
#'         different behavior to the right of the offset.
#'
#' @details The asymmetric bump function is defined differently to the right and right of the offset.
#'          For `x < offset` and `offset - x < h`, it is defined as exp(-1 / h - (x - offset)^2).
#'          For `|x - offset| < h`, it follows the standard bump function exp(-1 / (h - (x - offset)^2)).
#'          It is zero outside the specified interval, maintaining smooth transitions and compact support.
#'
#' @examples
#' \dontrun{
#' x <- seq(-2, 2, length.out = 100)
#' plot(x, AsymmetricBumpFunction(x), type = "l")
#' }
#' @export
right.asymmetric.bump.fn <- function(x, offset = 0, h = 1, q = 1) {

    n <- length(x)
    h2 <- h^2

    if ( n == 1 ) {
        ret <- 0
        if ( x >= offset ) {
            ret <- exp( -1 / h2 - (x - offset)^2 )
        } else if ( x < offset && offset - x < h ) {
            ret <- exp( -1 / (h2 - (x - offset)^2) )
            ret <- ret^q
        }
    } else {
        ret <- numeric(n)
        for ( i in seq(n) ) {
            if ( x[i] >= offset ) {
                ret[i] <- exp( -1 / h2 - (x[i] - offset)^2 )
            } else if ( x[i] < offset && offset - x[i] < h ) {
                ret[i] <- exp( -1 / (h2 - (x[i] - offset)^2) )
                ret[i] <- ret[i]^q
            }
        }
    }

    return(ret)
}

#' Generates the values of a synthetic 1d spline function over a uniform grid in the pre-specified range.
#'
#' This function generates a synthetic 1D spline with a specified number of local maxima.
#' The synthetic function is evaluated on a uniform grid over a specified x range.
#' The local minima are strategically placed to control the shape of the function.
#'
#' @param n.lmax Number of local maxima
#' @param x.min Minimum x value. Default 0.
#' @param x.max Maximum x value. Default 10.
#' @param y.min Minimum y value. Default 1.
#' @param y.max Maximum y value. Default 5.
#' @param n.grid Number of grid points. Default 400.
#' @param p.offset Fraction of interval between local maxima for local minima.
#'   Default 0.1.
#' @param alpha Shape parameter for beta distribution. Default 0.5.
#' @param method Spline method. Default 'natural'.
#'
#' @return List with components:
#' \itemize{
#'   \item \code{x}: Locations of grid points
#'   \item \code{y}: Function values at grid points
#'   \item \code{x.lmax}: Locations of local maxima
#'   \item \code{y.lmax}: Values of local maxima
#' }
#'
#' @examples
#' \dontrun{
#' synth <- synthetic.1D.spline(n.lmax = 5)
#' plot(synth$x, synth$y, type = "l")
#' points(synth$x.lmax, synth$y.lmax, col = "red")
#' }
#' @export
synthetic.1D.spline <- function(n.lmax,
                                x.min = 0,
                                x.max = 10,
                                y.min = 1,
                                y.max = 5,
                                n.grid = 400,
                                p.offset = 0.1,
                                alpha = 0.5,
                                method = "natural") {

    ## Generate local maxima locations
    x.lmax <- runif(n.lmax, x.min, x.max)
    x.lmax <- sort(x.lmax)

    ## Generate local maxima values
    y.lmax <- runif(n.lmax, y.min, y.max)

    ## Calculate local minima
    x.lmin <- numeric(n.lmax - 1)
    for(i in 1:(n.lmax-1)) {
        dx <- p.offset * (x.lmax[i+1] - x.lmax[i])
        int <- x.lmax[i] + dx
        int.end <- x.lmax[i+1] - dx
        l <- int.end - int
        x.lmin[i] <- int + rbeta(1, alpha, alpha) * l
    }

    ## Generate local minima values
    y.lmin <- runif(n.lmax-1, -y.max, -y.min)

    ## Create spline
    x.grid <- seq(x.min, x.max, length.out = n.grid)

    ## Create spline function
    x <- c(x.lmax, x.lmin)
    y <- c(y.lmax, y.lmin)
    o <- order(x)
    x <- x[o]
    y <- y[o]

    spline.fn <- splinefun(x, y, method = "natural")

    ## Evaluate spline at x grid
    y.grid <- spline.fn(x.grid)

    ## Return
    list(x = x.grid,
         y = y.grid,
         x.lmax = x.lmax,
         y.lmax = y.lmax)
}

#' Creates a Partition of Unity in 1D
#'
#' This function creates a partition of unity in one dimension using bump functions.
#' It is designed to cover a specified interval with a set of overlapping bump functions,
#' ensuring that the sum of these functions at any point in the interval equals one.
#'
#' @param x.center A numeric vector specifying the center points of the bump functions.
#' @param x.min The minimum value of the interval over which the partition is created (defaults to 0).
#' @param x.max The maximum value of the interval over which the partition is created (defaults to 10).
#' @param n.grid The number of points in the grid over the interval (defaults to 400).
#' @param C A scaling factor to control the spread of the bump functions. Defaults to 2.
#' @param q A q parameter in q-Gaussian
#' @param compact.supp Set to TRUE to use bump functions with compact support.
#'
#' @return A list with the following components
#' \itemize{
#'   \item \code{U}: A matrix with `n.grid` rows and `length(x.center)` columns, representing
#'         the values of the bump functions at each grid point. The sum of the column values at
#'         each row should be normalized to 1 to ensure a true partition of unity.
#'   \item \code{U.before}: U before normalization.
#'   \item \code{x.grid}: A uniform grid over the inverval \eqn{[x.min, x.max]}.
#'   \item \code{x.centers}: Sorted in the ascending order x.centers.
#'   \item \code{breaks}: Mid points between x.centers.
#' }
#'
#' @details The function uses asymmetric bump functions at the ends of the interval and
#'          symmetric bump functions for internal points. The spread of each bump function
#'          is determined by the distance to its neighboring centers, scaled by `C`.
#'          The function checks for valid numeric input and requires at least two center points.
#'
#' @examples
#' \dontrun{
#' x.center <- seq(1, 9, by = 2)
#' partition <- partition.of.unity.1D(x.center)
#' plot(seq(0, 10, length.out = 400), partition[,1], type = "l")
#' }
#' @export
partition.of.unity.1D <- function(x.center, x.min = 0, x.max = 10, n.grid = 400, C = 2, q = 3, compact.supp = FALSE) {

    n <- length(x.center)

    if ( n <= 1 ) {
        stop("x.center has to have at least 2 elements")
    }

    if ( !all(is.numeric(x.center)) ) {
        stop("x.center has to have a numeric value")
    }

    if ( !is.numeric(x.min) ) {
        stop("x.min has to have a numeric value")
    }

    if ( !is.numeric(x.max) ) {
        stop("x.max has to have a numeric value")
    }

    if ( !is.numeric(n.grid) ) {
        stop("n.grid has to have a numeric value")
    }

    ## Making sure the elements of x.center are sorted
    x.center <- sort(x.center)

    ## Identifying break points - points between x centers
    breaks <- numeric(n - 1)
    for ( i in 1:(n - 1) ) {
        breaks[i] <- x.center[i] + (x.center[i + 1] - x.center[i]) / 2
    }

    x.grid <- seq(x.min, x.max, length = n.grid)

    ## Creating the initial memebers of the partition of unity
    U <- matrix(nrow = n.grid, ncol = n)

    if ( compact.supp ) {
        U[,1] <- left.asymmetric.bump.fn(x.grid, x.center[1], C * (breaks[1] - x.center[1]), q = q)
        U[,n] <- right.asymmetric.bump.fn(x.grid, x.center[n], C * (x.center[n] - breaks[n - 1]), q = q)
        for ( i in 2:(n - 1) ) {
            h <- C * max(c(breaks[i] - x.center[i], x.center[i] - breaks[i - 1])) # choosing the larger of the distances to the mid point
            U[,i] <- bump.fn(x.grid, x.center[i], h, q = q)
        }
    } else {
        U[,1] <- qexp.gaussian(x.grid, x.center[1], C * (breaks[1] - x.center[1]), q)
        U[,n] <- qexp.gaussian(x.grid, x.center[n], C * (x.center[n] - breaks[n - 1]), q)
        for ( i in 2:(n - 1) ) {
            h <- C * max(c(breaks[i] - x.center[i], x.center[i] - breaks[i - 1])) # choosing the larger of the distances to the mid point
            U[,i] <- qexp.gaussian(x.grid, x.center[i], h, q)
        }
    }

    U.before <- U

    ## Normalizing U so that the sum of all functions sums to 1
    U <- U / rowSums(U)

    list(U = U,
         x.grid = x.grid,
         U.before = U.before,
         x.centers = x.center,
         breaks = breaks)
}

#' Bi-Gaussian defined as Gaussian with sigma.left for x <= mu and a Gaussian
#' with sigma.right for x > mu
#'
#' @param x           A vector of x-values.
#' @param mu          The location of the global maximum of the bi-gaussian function.
#' @param sigma.left  The standard deviation of the left part of the bi-gaussian function.
#' @param sigma.right The standard deviation of the right part of the bi-gaussian function.
bi.gaussian <- function(x, mu, sigma.left, sigma.right) {

    names(x) <- seq(length(x))

    x.right <- x[x >= mu]
    y.right <- exp(-(x.right - mu)^2/sigma.right^2)
    names(y.right) <- names(x.right)

    x.left <- x[x < mu]
    y.left <- exp(-(x.left - mu)^2/sigma.left^2)
    names(y.left) <- names(x.left)

    y = c(y.left, y.right)

    y[names(x)]
}

#' Creates a Synthetic Function with Specified Number of Local Maxima using bi-Gaussians.
#'
#' Generates a synthetic function with specified number of local maxima using bi-Gaussian functions
#' and a partition of unity. This function is useful for creating complex, non-linear synthetic data
#' for analysis and testing.
#'
#' @param n.lmax Number of local maxima.
#' @param x.lmax Locations of local maxima.
#' @param y.lmax Values of local maxima.
#' @param x.min Minimum x value for the grid (defaults to 0).
#' @param x.max Maximum x value for the grid (defaults to 10).
#' @param y.min Minimum y value for local maxima (defaults to 1).
#' @param y.max Maximum y value for local maxima (defaults to 5).
#' @param min.dist A positive numeric value specifying the minimum allowable
#'     distance between any two consecutive elements of x.
#' @param C.min.dist A real number between 0 and 1 used to specify min.dist if
#'     it is NULL. We use the formula min.dist = C.min.dist * (x.max - x.min) /
#'     n.lmax.
#' @param max.itr The maximal number of iterations for finding x.lmin so that
#'     the distance between consecutive elements is no less than min.dist.
#'
#' @param n.grid Number of grid points (defaults to 400).
#' @param C Scaling factor for the partition of unity (defaults to 2).
#' @param q A power parameter of q-Gaussian.
#' @param compact.support Set to TRUE to use partition of unity with the components
#'     with compact support and unbound support otherwise.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{x}: Vector of x locations for the grid points.
#'     \item \code{y}: Vector of function values at each grid point.
#'     \item \code{x.lmax}: Vector of x locations of the local maxima.
#'     \item \code{y.lmax}: Vector of y values of the local maxima.
#'   }
#'
#' @details The function generates a set of local maxima points and their corresponding y-values.
#'          It then constructs a partition of unity and a matrix of bi-Gaussian functions centered
#'          at these maxima. The bi-Gaussian function is defined as a Gaussian with sigma.left for
#'          x <= mu and a Gaussian with sigma.right for x > mu. The final synthetic function is
#'          the sum of these bi-Gaussians weighted by the partition of unity.
#'
#' @export
bi.gaussian.mixture <- function(n.lmax,
                                x.lmax = NULL,
                                y.lmax = NULL,
                                x.min = 0,
                                x.max = 10,
                                y.min = 1,
                                y.max = 5,
                                min.dist = NULL,
                                C.min.dist = 0.5,
                                max.itr = 1000,
                                n.grid = 400,
                                C = 1.5,
                                q = 3,
                                compact.support = FALSE) {
    if ( !is.numeric(x.min) ) {
        stop("x.min has to have a numeric value")
    }

    if ( !is.numeric(x.max) ) {
        stop("x.max has to have a numeric value")
    }

    if ( !is.numeric(n.grid) ) {
        stop("n.grid has to have a numeric value")
    }

    if ( x.max <= x.min ) {
        stop("x.max has to be greater than x.min")
    }

    if ( as.integer(n.lmax) != n.lmax ) {
        stop("n.lmax has to be an integer")
    }

    if ( n.lmax < 1 ) {
        stop("n.lmax has to be greater than 0")
    }

    if ( n.grid < 10 ) {
        stop("n.grid has to be greater than 9")
    }

    if ( is.null(min.dist) ) {
        min.dist <- C.min.dist * (x.max - x.min) / n.lmax
    }

    ## Generate local maxima locations
    if ( is.null(x.lmax) ) {

        x.lmax <- runif(n.lmax, x.min, x.max)
        x.lmax <- sort(x.lmax)

        itr <- 0
        while( itr < max.itr && min(diff(x.lmax)) < min.dist ) {
            x.lmax <- runif(n.lmax, x.min, x.max)
            x.lmax <- sort(x.lmax)
            itr <- itr + 1
        }

        if ( itr == max.itr ) {
            stop("Could not find x.lmax satisfying the min.dist condition")
        }
    } else {
        n.lmax <- length(x.lmax)
    }

    ## Generate local maxima values
    if ( is.null(y.lmax) ) {
        y.lmax <- runif(n.lmax, y.min, y.max)
    }

    ## Creating a partition of unity
    pU.res <- partition.of.unity.1D(x.lmax, x.min, x.max, n.grid, C, q, compact.support)
    x.grid <- pU.res$x.grid
    breaks <- pU.res$breaks
    U <- pU.res$U

    ## Creating a matrix of bi-Gaussians centered at the locations of the local
    ## maxima.
    biGaussians <- matrix(nrow = n.grid, ncol = n.lmax)

    ## The first and the last bi-Gaussian is set a bit differently as there are
    ## not break points to the left of the first local maximum and to the right
    ## of the last local maximum, so I am setting them up below.

    ## First bi-Gaussian
    sigma.left <- x.lmax[1] - x.min
    sigma.right <- breaks[1] - x.lmax[1]
    biGaussians[,1] <- y.lmax[1] * ggaussian(x.grid, offset = x.lmax[1], h = sigma.right)

    ## Last bi-Gaussian
    sigma.left <- x.lmax[n.lmax] - breaks[n.lmax - 1]
    sigma.right <- x.max - x.lmax[n.lmax]
    biGaussians[,n.lmax] <- y.lmax[n.lmax] * ggaussian(x.grid, offset = x.lmax[n.lmax], h = sigma.left)

    if ( n.lmax > 2 ) {
        for ( i in 2:(n.lmax - 1) ) {
            sigma.left <- x.lmax[i] - breaks[i - 1]
            sigma.right <- breaks[i] - x.lmax[i]
            biGaussians[,i] <- y.lmax[i] * bi.gaussian(x.grid, x.lmax[i], sigma.left, sigma.right)
        }
    }

    ## Creating the values of the synthetic function
    y.grid <- rowSums(biGaussians * U)

    list(x = x.grid,
         y = y.grid,
         x.lmax = x.lmax,
         y.lmax = y.lmax,
         U = U,
         min.dist = min.dist)
}

## --------------------------------------------------------------------------------------------------------------
##
## functions generating synthetic data in dimensions higher than 1
##
## --------------------------------------------------------------------------------------------------------------

#' Evaluates radial q-exponential Gaussian on set of points
#'
#' Evaluates radial q-exponential Gaussian cetnered at 'center' with standard deviation sigma
#' and power q on a set of points X.
#'
#' @param X A matrix of data frame of points at which the radia q-exponential Gaussian is to be evaluated.
#' @param center The center of the radial q-exponential Gaussian.
#' @param sigma  The standard deviation of the radial q-exponential Gaussian.
#' @param q      The power to which the radial q-exponential Gaussian is raised.
#'
#' @return Returns a vector of values of the radia q-exponential Gaussian at the points of X.
radial.qexp.gaussian <- function(X, center, sigma = 1, q = 1) {

    # Validate sigma
    if (!is.numeric(sigma) || sigma <= 0) {
        stop("sigma must be a positive numeric value.")
    }

    # Validate q
    if (!is.numeric(q) || q <= 0) {
        stop("q must be a positive numeric value.")
    }

    sigma2 <- sigma^2
    apply(X, 1, function(x) exp(-q * sum((x - center)^2) / sigma2))
}

#' Creates a Partition of Unity in dimensions greater than 1
#'
#' This function creates a partition of unity (PoU) over a subset S of \(R^d)\),
#' where d is the number of columns of S. The components of the partition of
#' unity are constructed from q-exponential Gaussian radial functions. If a matrix or data
#' frame, X.centers, of the centers of the components of the PoU, the centers of
#' the compoenents will be places at these points. If X.centers is NULL, then
#' the n.comp centers will be selected from the points of S.
#'
#' @param S A matrix or data frame of points over which the partition of unity will be evaluated.
#' @param X.centers A matrix or data frame of the centers of the components of the PoU.
#' @param n.comp The number of components of the PoU. Used when X.centers is NULL. In this case a random set of n.comp rows of S is taken to be X.centers.
#' @param C A scaling factor to control the spread of the bump functions. Defaults to 2.
#' @param q A q parameter in radial q-exponential Gaussian functions. Controls the steepness of the radial Gaussian function.
#'
#' @return A list with three components:
#' \itemize{
#'   \item \code{U}: A matrix of the values of each component of the PoU over the points of S. The return matrix has nrow(S) rows and n.comp columns.
#'   \item \code{nn.d}: A vector of distances to the nearest neighbor center.
#'   \item \code{nn.i}: A vector of indices of the nearest neighbor.
#' }
partition.of.unity.xD <- function(S,
                                  X.centers = NULL,
                                  n.comp = 3,
                                  C = 2,
                                  q = 3) {

    if ( !is.matrix(S) && !is.data.frame(S) ) {
        stop("S has to be  matrix of a data frame")
    }

    ## Check that S is numeric and that all entries are numeric (no NAs)
    if ( any(!is.numeric(S)) ) {
        stop("All elements of S have to be numeric. Detected non-numeric values.")
    }

    dim <- ncol(S)
    n.samples <- nrow(S)

    ## If X.centers is not NULL, check that its number of columns is the same as S
    ## If it is, set n.comp to the number of rows of X.centers
    if ( !is.null(X.centers) ) {

        if ( any(!is.numeric(X.centers)) ) {
            stop("All elements of X.centers have to be numeric. Detected non-numeric values.")
        }

        if ( !is.matrix(X.centers) && !is.data.frame(X.centers) ) {
            stop("X.centers has to be  matrix of a data frame")
        }

        if ( ncol(X.centers) != ncol(S) ) {
            stop("X.centers and S have to have the same number of columns. ncol(S):", ncol(S), " ncol(X.centers):",ncol(X.centers))
        }

        n.comp <- nrow(X.centers)

    } else {
        ## Selecting a set of n.comp centers in S
        ii <- sample(n.samples, size = n.comp)
        X.centers <- S[ii,]
    }

    ## Validate n.comp
    if (!is.numeric(n.comp) || n.comp <= 0 || n.comp != round(n.comp)) {
        stop("n.comp must be a positive integer.")
    }

    ## Validate C
    if (!is.numeric(C) || C <= 0) {
        stop("C must be a positive numeric value.")
    }

    ## Validate q
    if (!is.numeric(q) || q <= 0) {
        stop("q must be a positive numeric value.")
    }

    ## Computing the distances between the centers
    nn <- get.knn(X.centers, k = 1)
    nn.i <- nn$nn.index[,1]
    nn.d <- nn$nn.dist[,1]

    ## Creating the initial memebers of the partition of unity
    U <- matrix(nrow = n.samples, ncol = n.comp)
    for ( i in seq(n.comp) ) {
        sigma <- C * 0.5 * nn.d[i]
        u <- radial.qexp.gaussian(S, as.numeric(X.centers[i,]), sigma, q)
        u[!is.finite(u)] <- 0
        U[,i] <- u
    }

    U.before <- U

    ## Normalizing U so that the sum of all functions sums to 1
    U <- U / rowSums(U)

    list(U = U,
         U.before = U.before,
         X.centers = X.centers,
         nn.i = nn.i,
         nn.d = nn.d)
}


#' Creates a Synthetic Function with Specified Number of Local Maxima over a subset of R^d, d > 1.
#'
#' Generates a synthetic function with specified number of local maxima using a
#' mixture of Gaussians and a partition of unity. This function is useful for
#' creating complex, non-linear synthetic data for analysis and testing.
#'
#' @param S A matrix or data frame of points over which a the values of the constructed synthetic function will be evaluated.
#' @param X.lmax A matrix or data frame of the localtion of local maxima.
#' @param n.lmax The number of local maxima. Used when X.lmax is NULL. If this is the case, a random set of n.lmax rows of S is taken to be X.lmax.
#' @param y.lmax A vector of values the function being created supposet to have at the local maxima. A function that will be constructed is designed to have the values at the local maxima specified by y.lmax, but there are no warranty that it is going to be the case.
#' @param y.min Minimum y value for local maxima (defaults to 1).
#' @param y.max Maximum y value for local maxima (defaults to 5).
#' @param C A scaling factor to control the spread of the bump functions. Defaults to 2.
#' @param q A q parameter in radial q-exponential Gaussian functions. Controls the steepness of the radial Gaussian function.
#' @param with.partition.of.unity Set to TRUE to use a partition of unity.
#'
#' @details The function generates a set of local maxima and their corresponding
#'     y-values. It then constructs a partition of unity and a matrix of radial
#'     Gaussian functions centered at these local maxima. The final synthetic
#'     function is the sum of these radial Gaussians weighted by the partition
#'     of unity.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{X.lmax}: A matrix of the locations of the local maxima.
#'     \item \code{y.lmax}: A vector of function values at the local maxima.
#'     \item \code{y}: A vector of values of the created synthetic function.
#'   }
#'
#' @export
gaussian.mixture.xD <- function(S,
                                X.lmax = NULL,
                                n.lmax = 3,
                                y.lmax = NULL,
                                y.min = 1,
                                y.max = 5,
                                C = 2,
                                q = 1,
                                with.partition.of.unity = FALSE) {
                                ##min.dist = NULL,
                                ##C.min.dist = 0.5,

    ## Checking if S is a matrix of a data frame
    if ( !is.matrix(S) && !is.data.frame(S) ) {
        stop("S has to be  matrix of a data frame")
    }

    ## Checking that S is numeric and that all entries are numeric (no NAs)
    if ( any(!is.numeric(S)) ) {
        stop("All elements of S have to be numeric. Detected non-numeric values.")
    }

    dim <- ncol(S)
    n.samples <- nrow(S)

    ## If X.lmax is not NULL, check that its number of columns is the same as S
    ## If it is, set n.lmax to the number of rows of X.lmax
    if ( !is.null(X.lmax) ) {

        if ( any(!is.numeric(X.lmax)) ) {
            stop("All elements of X.lmax have to be numeric. Detected non-numeric values.")
        }

        if ( !is.matrix(X.lmax) && !is.data.frame(X.lmax) ) {
            stop("X.lmax has to be  matrix of a data frame")
        }

        if ( ncol(X.lmax) != ncol(S) ) {
            stop("X.lmax and S have to have the same number of columns. ncol(S):", ncol(S), " ncol(X.lmax):",ncol(X.lmax))
        }

        n.lmax <- nrow(X.lmax)

    } else {
        ## Selecting a set of n.lmax centers in S
        ii <- sample(n.samples, size = n.lmax)
        X.lmax <- S[ii,]
    }

    ## Validating n.lmax
    if (!is.numeric(n.lmax) || n.lmax <= 0 || n.lmax != round(n.lmax)) {
        stop("n.lmax must be a positive integer.")
    }

    ## Validating C
    if (!is.numeric(C) || C <= 0) {
        stop("C must be a positive numeric value.")
    }

    ## Validating q
    if (!is.numeric(q) || q <= 0) {
        stop("q must be a positive numeric value.")
    }

    ## Setting min.dist
    ## if ( is.null(min.dist) ) {
    ##     min.dist <- C.min.dist * (x.max - x.min) / n.lmax
    ## }

    ## Generating local maxima values
    if ( is.null(y.lmax) ) {
        y.lmax <- runif(n.lmax, y.min, y.max)
    }

    ## Generating local maxima values
    if ( is.null(y.lmax) ) {
        y.lmax <- runif(n.lmax, y.min, y.max)
    }

    ## Creating a partition of unity
    UoP.res <- partition.of.unity.xD(S, X.lmax, n.lmax, C, q)
    U <- UoP.res$U
    nn.d <- UoP.res$nn.d
    nn.i <- UoP.res$nn.i

    ## Creating a matrix of the values of the Gaussians centered at the locations of the local maxima.
    gaussians <- matrix(nrow = nrow(S), ncol = n.lmax)
    for ( i in seq(n.lmax) ) {
        gaussians[,i] <- y.lmax[i] * radial.qexp.gaussian(S, X.lmax[i,], sigma = nn.d[i], q = q)
    }

    ## Creating the values of the synthetic function
    if (with.partition.of.unity) {
        y <- rowSums(gaussians * U)
    } else {
        y <- rowSums(gaussians)
    }

    list(y = y,
         X.lmax = X.lmax,
         y.lmax = y.lmax,
         U = U)
         ##min.dist = min.dist)
}


#' Gaussian function in arbitrary dimension
#'
#' @param x Numeric vector of input values
#' @param mu Numeric vector of means (same length as x)
#' @param sigma Numeric vector of standard deviations (same length as x)
#' @param C Amplitude/scaling factor
#' @return Numeric value of the Gaussian function
#' @export
gaussian.xD <- function(x, mu, sigma, C) {
    C * exp(-(sum((x - mu)^2 / sigma^2)))
}

#' Gradient of the Gaussian function in arbitrary dimension
#'
#' @param x Numeric vector of input values
#' @param mu Numeric vector of means (same length as x)
#' @param sigma Numeric vector of standard deviations (same length as x)
#' @param C Amplitude/scaling factor
#' @return Numeric vector of gradient values
#' @export
gaussian.xD.grad <- function(x, mu, sigma, C) {
    ## Create symbolic variables for x
    x_sym <- paste0("x", seq(x))
    x_sym <- sapply(x_sym, as.symbol)

    ## Create the symbolic expression for the Gaussian function
    gaussian_expr <- substitute(gaussian.xD(x_sym, mu, sigma, C),
                                list(x_sym = x_sym, mu = mu, sigma = sigma, C = C))

    ## Compute the symbolic gradient
    grad_expr <- deriv(gaussian_expr, x_sym, function.arg = TRUE)

    ## Create a function from the symbolic gradient expression
    grad_func <- eval(grad_expr)

    ## Evaluate the gradient function at the given point
    grad_vals <- grad_func(x)

    return(grad_vals)
}

#' Generate Points on a Circle
#'
#' This function creates a set of points along a circle. The points can be either
#' equally spaced or randomly distributed along the circumference. Optional Gaussian
#' or Laplace noise can be added to the radial component.
#'
#' @param n An integer specifying the desired number of points to generate.
#' @param radius A positive numeric value specifying the radius of the circle (default is 1).
#' @param noise A non-negative numeric value specifying the level of noise to add to the points (default is 0.1).
#' @param type A character string, either "uniform" for equally spaced points or "random" for randomly distributed points along the circle.
#' @param noise.type A character string, either "normal" for Gaussian noise or "laplace" for Laplace (double exponential) noise.
#' @param seed An integer for the random seed. Default is NULL.
#'
#' @return A data frame with two columns, `x` and `y`, containing the coordinates of the generated points.
#'
#' @examples
#' # Generate 100 equally spaced points on a circle with radius 2
#' df <- generate.circle.data(100, radius = 2)
#'
#' # Generate 50 random points with Laplace noise
#' df_noisy <- generate.circle.data(50, noise = 0.1, type = "random", noise.type = "laplace")
#'
#' @export
generate.circle.data <- function(n,
                                 radius = 1,
                                 noise = 0.1,
                                 type = "random",
                                 noise.type = "laplace",
                                 seed = NULL) {
    ## Parameter checks
    if (!is.numeric(n) || n <= 0 || n != round(n))
        stop("n must be a positive integer.")
    if (!is.numeric(radius) || radius <= 0)
        stop("radius must be a positive number.")
    if (!is.numeric(noise) || noise < 0)
        stop("noise must be a non-negative number.")
    noise.type <- match.arg(noise.type, c("normal", "laplace"))
    type <- match.arg(type, c("uniform", "random"))

    ## Set seed if provided
    if (!is.null(seed)) set.seed(seed)

    # Generate angles
    if (type == "uniform") {
        angles <- seq(0, 2 * pi, length.out = n + 1)[-1]
    } else {
        angles <- sort(stats::runif(n, min = 0, max = 2 * pi))
    }

    ## Generating noise
    eps <- numeric(n)
    if (noise > 0) {
        if (noise.type == "laplace") {
            eps <- rlaplace(n, location = 0, scale = noise / sqrt(2))
        } else {
            eps <- stats::rnorm(n, mean = 0, sd = noise)
        }
    }

    ## Calculating coordinates
    x.noise.free <- radius * cos(angles)
    y.noise.free <- radius * sin(angles)
    x <- (radius + eps) * cos(angles)
    y <- (radius + eps) * sin(angles)

    data.frame(x = x,
               y = y,
               x.noise.free = x.noise.free,
               y.noise.free = y.noise.free,
               angles)
}

#' Generate Noisy Circle Data in Higher Dimensions
#'
#' This function generates data points from a noisy circle embedded in a higher-dimensional space.
#'
#' @param n An integer. The number of data points to generate.
#' @param dim An integer greater than 2. The dimension of the space to embed the circle in.
#' @param radius A positive number. The radius of the circle.
#' @param noise A non-negative number. The scale of the noise to add to the data.
#' @param type A string, either "uniform" or "random". Determines how angles are generated.
#' @param noise.type A string, either "normal" or "laplace". Determines the distribution of the noise.
#'
#' @return A matrix with n rows and dim columns. Each row represents a point in R^dim.
#'
#' @examples
#' # Generate 100 points on a noisy circle in 3D space
#' data <- generate.noisy.circle.embedding(100, dim = 3, radius = 1, noise = 0.1)
#'
#' @export
generate.noisy.circle.embedding <- function(n,
                                            dim = 3,
                                            radius = 1,
                                            noise = 0,
                                            type = "random",
                                            noise.type = "laplace") {
    if (!is.numeric(n) || n <= 0 || n != round(n))
        stop("n must be a positive integer.")

    if (!is.numeric(dim) || dim < 2 || dim != round(dim))
        stop("dim must be an integer greater than 2.")

    if (!is.numeric(radius) || radius <= 0)
        stop("radius must be a positive number.")

    if (!is.numeric(noise) || noise < 0)
        stop("noise must be a non-negative number.")

    noise.type <- match.arg(noise.type, c("normal", "laplace"))
    type <- match.arg(type, c("uniform", "random"))

    ## Generating angles
    if (type == "uniform") {
        angles <- seq(0, 2 * pi, length.out = n + 1)[-1]  # Exclude the last point to avoid duplication
    } else {
        angles <- sort(runif(n, min = 0, max = 2 * pi))
    }

    ## Generating noise
    m <- dim - 2
    eps <- NA
    eps2 <- NA
    if (noise > 0) {
        if (noise.type == "laplace") {
            eps <- rlaplace(n, location = 0, scale = noise / sqrt(2))
            ## Noise in other dimensions
            if (m > 0) {
                eps2 <- matrix(rlaplace(n * m, location = 0, scale = noise / sqrt(2)), nrow = n, ncol = m)
            }
        } else {
            eps <- rnorm(n, 0, noise)
            ## Noise in other dimensions
            if (m > 0) {
                eps2 <- matrix(rnorm(n * m, 0, noise), nrow = n, ncol = m)
            }
        }
    } else {
        eps2 <- matrix(0, nrow = n, ncol = m)
    }

    ## Calculating coordinates
    x <- (radius + eps) * cos(angles)
    y <- (radius + eps) * sin(angles)

    cbind(x, y, eps2)
}




#' Generate Equally Spaced Points on a Circle without any noise.
#'
#' This function creates a set of points that are equally spaced along a circle. The first and last point are the same.
#'
#' @param n The desired number of points to generate.
#' @param radius The radius of the circle (default is 1).
#'
#' @return A data frame with two columns, `x` and `y`, containing the coordinates of the generated points.
generate.circle <- function(n, radius = 1) {

    angles <- seq(0, 2 * pi, length.out = n + 1)
    angles <- angles[-1]

    x <- radius * cos(angles)
    y <- radius * sin(angles)

    data.frame(x = x, y = y)
}


#' Torus knot
#'
#' @param n The desired number of points to generate.
#' @param k Torus knot first parameter.
#' @param m Torus knot second parameter.
#'
v1.torus.knot <- function(n, k, m) {

    angles <- seq(0, 2 * pi, length.out = n)

    x <- sin(angles) + k * cos(k * angles)
    y <- cos(angles) - k * sin(k * angles)
    z <- -sin(m * angles)

    cbind(x, y, z)
}

torus.knot <- function(n, k, m) {
  u <- seq(0, 2 * pi, length.out = n + 1)[-1]  # Angles along the main circle
  v <- seq(0, 2 * pi, length.out = n)          # Angles along the tube

  x <- (cos(u) + k * cos(k * u)) * cos(m * v)
  y <- (cos(u) + k * cos(k * u)) * sin(m * v)
  z <- -sin(k * u)

  cbind(x, y, z)
}

#' Generate a Trefoil Knot in 3D
#'
#' @param n The desired number of points to generate.
#' @param scale The scale of the knot.
#' @param type  A character vector: "uniform" or "random".
#'
generate.trefoil.knot <- function(n = 100, scale = 1, type = "uniform") {

    if ( type == "uniform" ) {
        t <- seq(0, 2 * pi, length.out = n)
    } else {
        t <- runif(n, min = 0, max = 2 * pi)
    }

    x <- scale * sin(t) + 2 * sin(2*t)
    y <- scale * cos(t) - 2 * cos(2*t)
    z <- scale * -sin(3*t)

    cbind(x, y, z)
}

#' Generates clustered data (diagonal covariance)
#'
#' @param n.clusters    Integer, number of clusters.
#' @param cluster.sizes Integer vector of length n.clusters with sizes per cluster.
#' @param dimensions    Integer, number of dimensions (columns).
#' @param means         Optional numeric matrix \code{[n.clusters x dimensions]} of cluster means.
#' @param covariances   Optional list of length n.clusters; each element either
#'                      (i) a diagonal covariance matrix \code{(dimensions x dimensions)},
#'                      or (ii) a numeric vector of variances of length `dimensions`.
#'                      If NULL, uses 0.5 * I.
#' @param separation    Numeric, separation scale used when `means` is NULL.
#' @return Numeric matrix with `sum(cluster.sizes)` rows and `dimensions` columns.
generate.clustered.data <- function(n.clusters,
                                    cluster.sizes,
                                    dimensions = 2,
                                    means = NULL,
                                    covariances = NULL,
                                    separation = 2) {
  ## --- Error checking ---
  if (length(cluster.sizes) != n.clusters) {
    stop("Number of cluster sizes must match n.clusters.")
  }
  if (!is.null(means)) {
    if (!is.matrix(means) || nrow(means) != n.clusters || ncol(means) != dimensions) {
      stop("`means` must be a numeric matrix with nrow = n.clusters and ncol = dimensions.")
    }
  }
  if (!is.null(covariances)) {
    if (!is.list(covariances) || length(covariances) != n.clusters) {
      stop("`covariances` must be a list of length n.clusters (or NULL).")
    }
  }

  ## --- Preallocate result ---
  N <- sum(cluster.sizes)
  data <- matrix(0, nrow = N, ncol = dimensions)
  start.index <- 1L

  ## --- Helper to normalize covariance to sd vector (sqrt of diagonal variances) ---
  cov_to_sd <- function(cov_i) {
    if (is.null(cov_i)) {
      rep(sqrt(0.5), dimensions)
    } else if (is.vector(cov_i) && length(cov_i) == dimensions) {
      if (any(cov_i < 0)) stop("Covariance (variance vector) must be non-negative.")
      sqrt(cov_i)
    } else if (is.matrix(cov_i) && all(dim(cov_i) == c(dimensions, dimensions))) {
      if (!isTRUE(all(cov_i == diag(diag(cov_i))))) {
        stop("Only diagonal covariance supported: off-diagonal entries must be zero.")
      }
      v <- diag(cov_i)
      if (any(v < 0)) stop("Diagonal covariance must be positive semidefinite.")
      sqrt(v)
    } else {
      stop("Each covariance must be either a length-`dimensions` variance vector or a diagonal matrix.")
    }
  }

  ## --- Generate each cluster ---
  for (i in seq_len(n.clusters)) {
    n_i <- as.integer(cluster.sizes[i])

    # Means
    center <- if (is.null(means)) {
      runif(dimensions, -separation, separation) * i
    } else {
      means[i, ]
    }

    # Standard deviations per dimension
    sds <- if (is.null(covariances)) {
      rep(sqrt(0.5), dimensions)  # default: 0.5 * I
    } else {
      cov_to_sd(covariances[[i]])
    }

    # Sample: Z ~ N(0, I), then scale and shift: X = center + Z * sds
    Z <- matrix(rnorm(n_i * dimensions), nrow = n_i, ncol = dimensions)
    cluster.data <- sweep(Z, 2, sds, `*`)
    cluster.data <- sweep(cluster.data, 2, center, `+`)

    idx <- start.index:(start.index + n_i - 1L)
    data[idx, ] <- cluster.data
    start.index <- start.index + n_i
  }

  data
}

#' Tests graph Morse-Smale complex construct (its persistence and reconstruction) on synthetic datasets
#'
#' @param n.datasets  The number of datasets to generate.
#' @param n.pts       The number of points the synthetic cloud dataset suppose to have.
#' @param n.side      If the data is to be a sample from some grid, then n.side is the number of points on the axis of each dimension, thus the total number of grid points is n.size^dim.
#' @param dim         The dimension of R^dim in which data is constructed.
#' @param type        The method used for the construction.
#' @param Ks          A vector of positive integer values specifying the range of k values.
#' @param n.lmax      The number of gaussians to use.
#' @param y.min       The mimimum value of the range of y values.
#' @param y.max       The maximum value of the range of y values.
#' @param out.dir     The output directory where the generated data is stored.
#' @param verbose     Set it to TRUE, to print progress messages.
#'
#' @return A list containing the length of the longest stretch of isomorphic MS graphs for each generated dataset
#' @export
test.graph.MS.cx.on.synth.cloud.data <- function(n.datasets,
                                                 n.pts,
                                                 n.side = NULL,
                                                 dim = 2,
                                                 type = "grid",
                                                 Ks,
                                                 n.lmax = 10,
                                                 y.min = 1,
                                                 y.max = 15,
                                                 out.dir,
                                                 verbose = TRUE
                                                 ) {
    if (verbose) {
        routine.ptm <- proc.time()
    }

    if (!file.exists(out.dir))
        dir.create(out.dir)

    length.of.longest.stretch.of.isomorphic.MS.graphs.list <- list()
    for (i in seq(n.datasets)) {
        if (verbose) {
            cat("\ri:",i)
        }

        ##
        ## Creating a dataset
        ##
        X <- create.ENPs.grid.2D(n = n.side, x1.range = c(0,1), x2.range = c(0,1), f=0)
        ## plot(X, pch=".", las = 1)
        gm.res <- gaussian.mixture.xD(X,
                                      X.lmax = NULL,
                                      n.lmax = n.lmax,
                                      y.lmax = NULL,
                                      y.min = y.min,
                                      y.max = y.max,
                                      C = 2,
                                      q = 1,
                                      with.partition.of.unity = FALSE)
        ## plot2d.cont(X, gm.res$y, use.layout = TRUE)
        y <- gm.res$y

        ## Creating the reference graph and the associated MS graph
        ref.k <- 4
        ref.graph <- IW.kNN.graph(X, ref.k)$adj_list
        ref.MS.graph <- graph.MS.cx(ref.graph, y)$MS_graph

        ##
        ## Sampling from the domain of the function and restricting the function to that sample
        ##
        ii <- sample(nrow(X), size = n.pts)
        init.X <- X
        init.y <- y
        X <- X[ii,]
        y <- y[ii]

        pers.res <- morse.smale.graph.persistence(X,
                                                  y,
                                                  Ks = Ks,
                                                  ref.MS.graph = ref.MS.graph,
                                                  allow.disconnected.graph = TRUE,
                                                  verbose = FALSE)

        ## names(pers.res)
        ## [1] "nerve.graph"
        ## [2] "MS.graph"
        ## [3] "are.consecutive.graphs.isomorphic"
        ## [4] "first.index.of.longest.stretch.of.isomorphic.MS.graphs"
        ## [5] "length.of.longest.stretch.of.isomorphic.MS.graphs"
        ## [6] "is.isomorphic.with.ref.MS.graph"

        ## pers.res$first.index.of.longest.stretch.of.isomorphic.MS.graphs
        length.of.longest.stretch.of.isomorphic.MS.graphs.list[[i]] <- pers.res$length.of.longest.stretch.of.isomorphic.MS.graphs

        file <- paste0(out.dir,"/rgrid_",i,".rda")
        save(init.X, init.y, X, y, ii, pers.res, file = file)
    }

    if (verbose) {
        txt <- sprintf("\rTotal elapsed time")
         elapsed.time(routine.ptm, txt, with.brackets = FALSE)
    }

    length.of.longest.stretch.of.isomorphic.MS.graphs.list
}

#' Generate a 1D Circular Gaussian Mixture
#'
#' This function generates a smooth, periodic function over a circle by creating
#' a mixture of Gaussian distributions. The function ensures continuity and
#' smoothness at all points, including the 0/2pi connection point.
#'
#' @param n.points An integer specifying the number of points to generate along
#'   the circle. Must be positive. Default is 100.
#' @param x.knots A numeric vector specifying the positions of the Gaussian
#'   centers on the circle, in radians. Values are automatically wrapped to
#'   \code{[0, 2pi]}. Must have the same length as y.knots.
#' @param y.knots A numeric vector specifying the heights of the Gaussian
#'   centers. Must have the same length as x.knots.
#' @param sd.knot A positive numeric value specifying the standard deviation
#'   for all Gaussian distributions. Default is 0.5.
#' @param out.dir An optional character string specifying the directory to save
#'   the results. If NULL (default), results are not saved to a file.
#'
#' @return A list containing the following elements:
#'   \item{x}{A numeric vector of x values (angles in radians) from 0 to 2*pi.}
#'   \item{y}{A numeric vector of corresponding y values of the mixture function.}
#'   \item{file.name}{A character string with the name of the file where results
#'     are saved (if out.dir is specified).}
#'   \item{params}{A list of all input parameters used to generate the mixture.}
#'
#' @examples
#' # Generate a simple circular Gaussian mixture
#' result <- generate.1d.circular.gaussian.mixture(
#'   n.points = 200,
#'   x.knots = c(0, pi, 1.5*pi),
#'   y.knots = c(5, 8, 2.5),
#'   sd.knot = 0.3
#' )
#'
#' # Plot the result
#' plot(result$x, result$y, type = "l", xlab = "Angle (radians)", ylab = "Value")
#'
#' @export
generate.1d.circular.gaussian.mixture <- function(n.points = 100,
                                                  x.knots = c(0, pi, 1.5*pi),
                                                  y.knots = c(5, 8, 2.5),
                                                  sd.knot = 0.5,
                                                  out.dir = NULL) {
    ## Input validation
    if (length(x.knots) != length(y.knots)) {
        stop("x.knots and y.knots must have the same length")
    }
    if (n.points <= 0 || sd.knot <= 0) {
        stop("n.points and sd.knot must be positive")
    }
    if (!is.null(out.dir) && !dir.exists(out.dir)) {
        stop("Specified output directory does not exist")
    }

    ## Normalize x.knots to [0, 2pi] interval
    x.knots <- x.knots %% (2*pi)

    ## Set up grid
    x.grid <- seq(0, 2*pi, length.out = n.points)

    ## Generate synthetic data
    y.smooth <- circular.synthetic.mixture.of.gaussians(x.grid, x.knots, y.knots, sd.knot = sd.knot)

    ## Generate a unique identifier
    unique.id <- paste(
        paste(round(x.knots, 2), collapse = "_"),
        paste(y.knots, collapse = "_"),
        sd.knot,
        sep = "_"
    )

    ## Create the file name
    file.name <- sprintf("circular_gaussian_mixture_1d_n%d_k%d_sd%.2f_%s.rda",
                         n.points,
                         length(x.knots),
                         sd.knot,
                         unique.id)

    ## Prepare results
    results <- list(
        x = x.grid,
        y = y.smooth,
        file.name = file.name,
        params = list(
            n.points = n.points,
            x.knots = x.knots,
            y.knots = y.knots,
            sd.knot = sd.knot,
            out.dir = out.dir
        )
    )

    ## Save results if out.dir is specified
    if (!is.null(out.dir)) {
        file.path <- file.path(out.dir, file.name)
        save(results, file = file.path)
        cat("Results saved to:", file.path, "\n")
    }

    return(results)
}

#' Create a Synthetic Mixture of Circular Gaussians
#'
#' This function generates a mixture of periodic Gaussian distributions on a
#' circle. It ensures that the resulting function is smooth and continuous at
#' all points, including the 0/2pi connection point.
#'
#' @param x A numeric vector of angles (in radians) at which to evaluate the
#'   mixture function.
#' @param x.knot A numeric vector specifying the centers of the Gaussian
#'   distributions on the circle, in radians.
#' @param y.knot A numeric vector specifying the heights of the Gaussian
#'   distributions. Must have the same length as x.knot.
#' @param sd.knot A positive numeric value specifying the standard deviation
#'   for all Gaussian distributions in the mixture.
#'
#' @return A numeric vector of the same length as x, containing the values of
#'   the Gaussian mixture function evaluated at the input angles.
#'
#' @details
#' The function uses a periodic Gaussian distribution to ensure smoothness and
#' continuity around the entire circle. It calculates the shortest distance
#' between points on the circle in both clockwise and counterclockwise
#' directions, ensuring proper wrapping at the 0/2pi point.
#'
#' @examples
#' # Generate a simple circular Gaussian mixture
#' x <- seq(0, 2*pi, length.out = 100)
#' y <- circular.synthetic.mixture.of.gaussians(
#'   x,
#'   x.knot = c(0, pi, 1.5*pi),
#'   y.knot = c(5, 8, 2.5),
#'   sd.knot = 0.3
#' )
#'
#' # Plot the result
#' plot(x, y, type = "l", xlab = "Angle (radians)", ylab = "Value")
#'
#' @export
circular.synthetic.mixture.of.gaussians <- function(x, x.knot, y.knot, sd.knot) {
    n.knots <- length(x.knot)

    ## Define the periodic Gaussian function
    periodic.gaussian <- function(x, mu, sigma) {
        d1 <- (x - mu) %% (2*pi)
        d2 <- (mu - x) %% (2*pi)
        d <- pmin(d1, d2)
        return(exp(-(d^2) / (2 * sigma^2)))
    }

    ## Define the synthetic function
    synthetic.function <- function(x.val) {
        result <- 0
        for (i in 1:n.knots) {
            result <- result + y.knot[i] * periodic.gaussian(x.val, x.knot[i], sd.knot)
        }
        return(result)
    }

    return(sapply(x, synthetic.function))
}


#' Generate Random Partition of an Integer
#'
#' This function generates a random partition of an integer n into a fixed number
#' of parts (a), where each part is at least m. It uses a sampling-based approach
#' to generate the partition.
#'
#' @param n A positive integer to be partitioned
#' @param a The number of parts in the partition
#' @param m The minimum value for each part
#'
#' @return A numeric vector of length a containing the partition, where each element
#'         is greater than or equal to m and the sum equals n. Returns NULL if the
#'         partition is impossible (when n < a * m).
#'
#' @details
#' The algorithm works by first calculating r = n - (a * m), then randomly sampling
#' (a-1) indices from 1 to r to create partition points. The differences between
#' consecutive partition points (plus m) become the values in the final partition.
#'
#' @examples
#' generate.partition(20, 4, 3)  # Generates a partition of 20 into 4 parts, each >= 3
#' generate.partition(10, 3, 2)  # Generates a partition of 10 into 3 parts, each >= 2
#'
#' @export
generate.partition <- function(n, a, m) {
  # Input validation
    if (!is.numeric(n) || !is.numeric(a) || !is.numeric(m)) {
        stop("All arguments must be numeric")
    }
    if (n <= 0 || a <= 0 || m <= 0) {
        stop("All arguments must be positive")
    }
    if (floor(n) != n || floor(a) != a || floor(m) != m) {
        stop("All arguments must be integers")
    }

    ## Check if partition is possible
    if (n < a * m) {
        return(NULL)
    }

    ## Calculate remainder to be partitioned
    r <- n - (a * m)

    ## If r is 0, return vector of all m
    if (r == 0) {
        return(rep(m, a))
    }

    ## Sample a-1 indices from 1 to r
    if (a - 1 > r) {
        # If we need more samples than available, use all available and pad with repeats
        indices <- sort(c(1:r, sample(1:r, a - 1 - r, replace = TRUE)))
    } else {
        indices <- sort(sample(1:r, a-1))
    }

    ## Initialize result vector
    result <- numeric(a)

    ## Calculate first part
    result[1] <- indices[1]

    ## Calculate middle parts
    if (a > 2) {
        for (i in 2:(a-1)) {
            result[i] <- indices[i] - indices[i-1]
        }
    }

    ## Calculate last part
    result[a] <- r - indices[a-1]

    ## Add minimum value m to each part
    result <- result + m

    return(result)
}




#' Transform Continuous Values to Probabilities via Normalization and Logistic Function
#'
#' This function converts continuous-valued input into probabilities through a two-step process:
#' 1. Min-max normalization to a specified range
#' 2. Logistic transformation to convert normalized values to probabilities
#'
#' The transformation is useful for converting arbitrary continuous values (like smoothed
#' data or similarity scores) into valid probabilities that can be used for binary sampling
#' or probabilistic modeling.
#'
#' @param y Numeric vector of continuous values to be transformed. Must contain more than
#'   one unique finite value.
#' @param y.min Numeric scalar indicating the minimum value for the normalization range.
#'   Must be less than y.max. Default is -3.
#' @param y.max Numeric scalar indicating the maximum value for the normalization range.
#'   Must be greater than y.min. Default is 3.
#'
#' @return A numeric vector of the same length as y, containing probabilities in the
#'   range (0,1).
#'
#' @details
#' The function first applies min-max normalization to scale the input values to the
#' range \code{[y.min, y.max]}. It then applies the inverse logit (logistic) function
#' \code{1/(1 + exp(-x))} to convert these normalized values to probabilities.
#'
#' The default range \code{[-3, 3]} is chosen because these values, when transformed by the
#' logistic function, result in probabilities of approximately 0.05 and 0.95
#' respectively, providing a good spread of probability values while avoiding
#' extreme probabilities too close to 0 or 1.
#'
#' @examples
#' # Transform smoothed values into probabilities
#' y.smooth <- c(-1, 0, 0.5, 1, 2)
#' probs <- normalize.and.inv.logit.fn(y.smooth)
#'
#' # Use custom range for normalization
#' probs.wide <- normalize.and.inv.logit.fn(y.smooth, y.min = -5, y.max = 5)
#'
#' @section
#' Throws an error if:
#' * Input contains non-finite values
#' * Input contains fewer than two unique values
#' * y.min is greater than or equal to y.max
#'
#' @seealso
#' \code{\link{rbinom}} for generating binary samples from the resulting probabilities
#'
#' @export
normalize.and.inv.logit.fn <- function(y, y.min = -3, y.max = 3) {
    # Input validation
    if (!all(is.finite(y))) {
        stop("Input contains non-finite values")
    }
    if (length(unique(y)) <= 1) {
        stop("Input must contain more than one unique value")
    }
    if (y.min >= y.max) {
        stop("y.min must be less than y.max")
    }

    # Proceed with transformation
    y <- (y - min(y)) / (max(y) - min(y)) * (y.max - y.min) + y.min
    y <- 1 / (1 + exp(-y))
    return(y)
}




## create.gaussian.mixture <- function(x1 = 0.25, y1 = 0.25,  ## location of first Gaussian
##                                     x2 = 0.75, y2 = 0.75,    ## location of second Gaussian
##                                     A1 = 1.0,                ## amplitude of first Gaussian
##                                     A2 = 1.0,                ## amplitude of second Gaussian
##                                     sigma = 0.2) {           ## shared standard deviation
##     ## Input validation
##     if (sigma <= 0) stop("sigma must be positive")
##     if (A1 <= 0 || A2 <= 0) stop("amplitudes must be positive")

##     ## Create list of functions to return
##     result <- list()

##     ## Function for the Gaussian mixture
##     result$f <- function(x, y) {
##         ## First Gaussian
##         g1 <- A1 * exp(-((x - x1)^2 + (y - y1)^2)/(2 * sigma^2))
##         ## Second Gaussian
##         g2 <- A2 * exp(-((x - x2)^2 + (y - y2)^2)/(2 * sigma^2))
##         ## Sum of Gaussians
##         return(g1 + g2)
##     }

##     ## Function for the gradient
##     result$gradient <- function(x, y) {
##         ## First Gaussian and its derivatives
##         g1 <- A1 * exp(-((x - x1)^2 + (y - y1)^2)/(2 * sigma^2))
##         dg1_dx <- g1 * (-1/sigma^2) * (x - x1)
##         dg1_dy <- g1 * (-1/sigma^2) * (y - y1)

##         ## Second Gaussian and its derivatives
##         g2 <- A2 * exp(-((x - x2)^2 + (y - y2)^2)/(2 * sigma^2))
##         dg2_dx <- g2 * (-1/sigma^2) * (x - x2)
##         dg2_dy <- g2 * (-1/sigma^2) * (y - y2)

##         ## Return gradient vector (partial derivatives with respect to x and y)
##         return(c(dx = dg1_dx + dg2_dx,
##                 dy = dg1_dy + dg2_dy))
##     }

##     return(result)
## }


#' Generate 2D Gaussian Mixture
#'
#' Creates a 2D Gaussian mixture function with flexible component placement and characteristics.
#' Returns an object of class "gaussian_mixture" with associated plotting methods.
#'
#' @param n.components Integer specifying the number of Gaussian components
#' @param x.range Numeric vector c(min, max) specifying the range for x values
#' @param y.range Numeric vector c(min, max) specifying the range for y values
#' @param centers.strategy Character string specifying how to position components:
#'   * "random": uniformly random positions (default)
#'   * "grid": regular grid arrangement
#'   * "custom": user-specified centers
#' @param custom.centers Matrix with n.components rows and 2 columns (x, y coordinates)
#'   when centers.strategy = "custom"
#' @param amplitudes.strategy Character string specifying amplitude generation:
#'   * "mixed": positive and negative amplitudes (default)
#'   * "positive": only positive amplitudes
#'   * "custom": user-specified amplitudes
#' @param custom.amplitudes Numeric vector of length n.components when
#'   amplitudes.strategy = "custom"
#' @param amplitude.range Numeric vector c(min, max) for amplitude magnitudes
#' @param sd.strategy Character string specifying how to generate component SDs:
#'   * "fixed": constant SD for all components (default)
#'   * "random": random SDs within range
#'   * "custom": user-specified SDs
#' @param sd.value Numeric specifying SD when sd.strategy = "fixed"
#' @param sd.range Numeric vector c(min, max) when sd.strategy = "random"
#' @param custom.sds Numeric vector or 2x2 matrix. If vector, creates isotropic
#'   Gaussians; if matrix, specifies covariance for each component
#' @param correlation Numeric between -1 and 1 specifying correlation between
#'   x and y for all components (ignored if custom.sds is a matrix)
#' @param edge.buffer Numeric between 0 and 0.5 specifying buffer from edges
#' @param seed Optional integer for random number generation
#' @param verbose Logical indicating whether to print generation details
#'
#' @return An object of class "gaussian_mixture" containing:
#' \itemize{
#'   \item \code{f}: Function f(x, y) that evaluates the mixture at given coordinates
#'   \item \code{n.components}: Number of components
#'   \item \code{centers}: Matrix of component centers (n.components x 2)
#'   \item \code{amplitudes}: Vector of component amplitudes
#'   \item \code{covariances}: List of 2x2 covariance matrices for each component
#'   \item \code{x.range}: X range used
#'   \item \code{y.range}: Y range used
#'   \item \code{parameters}: List of all parameters used
#' }
#'
#' @examples
#' \dontrun{
#' # Generate random 2D mixture
#' gm1 <- get.gaussian.mixture.2d(n.components = 3)
#' plot(gm1)
#'
#' # Generate grid arrangement with custom amplitudes
#' gm2 <- get.gaussian.mixture.2d(
#'   n.components = 4,
#'   centers.strategy = "grid",
#'   amplitudes.strategy = "custom",
#'   custom.amplitudes = c(1, -0.5, 0.8, -0.3)
#' )
#'
#' # Custom centers with anisotropic Gaussians
#' gm3 <- get.gaussian.mixture.2d(
#'   n.components = 2,
#'   centers.strategy = "custom",
#'   custom.centers = matrix(c(0.3, 0.3, 0.7, 0.7), ncol = 2, byrow = TRUE),
#'   correlation = 0.5
#' )
#' }
#' @export
get.gaussian.mixture.2d <- function(
    n.components = 3,
    x.range = c(0, 1),
    y.range = c(0, 1),
    centers.strategy = c("random", "grid", "custom"),
    custom.centers = NULL,
    amplitudes.strategy = c("mixed", "positive", "custom"),
    custom.amplitudes = NULL,
    amplitude.range = c(0.5, 2),
    sd.strategy = c("fixed", "random", "custom"),
    sd.value = 0.1,
    sd.range = c(0.05, 0.2),
    custom.sds = NULL,
    correlation = 0,
    edge.buffer = 0.1,
    seed = NULL,
    verbose = FALSE) {

    # Set seed if provided
    if (!is.null(seed)) set.seed(seed)

    # Match arguments
    centers.strategy <- match.arg(centers.strategy)
    amplitudes.strategy <- match.arg(amplitudes.strategy)
    sd.strategy <- match.arg(sd.strategy)

    # Validate inputs
    if (n.components <= 0) stop("n.components must be positive")
    if (length(x.range) != 2 || x.range[2] <= x.range[1])
        stop("x.range must be c(min, max) with min < max")
    if (length(y.range) != 2 || y.range[2] <= y.range[1])
        stop("y.range must be c(min, max) with min < max")
    if (abs(correlation) > 1) stop("correlation must be between -1 and 1")
    if (edge.buffer < 0 || edge.buffer >= 0.5)
        stop("edge.buffer must be between 0 and 0.5")

    # Generate centers
    centers <- generate_2d_centers(centers.strategy, n.components, x.range, y.range,
                                   edge.buffer, custom.centers)

    # Generate amplitudes
    amplitudes <- generate_amplitudes_2d(amplitudes.strategy, n.components,
                                         amplitude.range, custom.amplitudes)

    # Generate covariance matrices
    covariances <- generate_covariances_2d(sd.strategy, n.components, sd.value,
                                           sd.range, custom.sds, correlation)

    # Create the mixture function
    f <- function(x, y) {
        result <- 0
        for (i in 1:n.components) {
            center <- centers[i, ]
            cov.inv <- solve(covariances[[i]])
            diff <- c(x - center[1], y - center[2])
            exponent <- -0.5 * sum(diff * (cov.inv %*% diff))
            result <- result + amplitudes[i] * exp(exponent)
        }
        result
    }

    # Create the return object
    gm <- list(
        f = f,
        n.components = n.components,
        centers = centers,
        amplitudes = amplitudes,
        covariances = covariances,
        x.range = x.range,
        y.range = y.range,
        parameters = list(
            centers.strategy = centers.strategy,
            amplitudes.strategy = amplitudes.strategy,
            amplitude.range = amplitude.range,
            sd.strategy = sd.strategy,
            sd.value = sd.value,
            sd.range = sd.range,
            correlation = correlation,
            edge.buffer = edge.buffer,
            seed = seed
        )
    )

    # Set class
    class(gm) <- c("gaussian_mixture", "list")

    if (verbose) {
        cat("Generated 2D Gaussian mixture with:\n")
        cat("  Components:", n.components, "\n")
        cat("  Centers strategy:", centers.strategy, "\n")
        cat("  Amplitudes strategy:", amplitudes.strategy, "\n")
        cat("  SD strategy:", sd.strategy, "\n")
    }

    return(gm)
}

#' Generate 2D component centers
#' @keywords internal
generate_2d_centers <- function(strategy, n.components, x.range, y.range,
                                edge.buffer, custom.centers) {
    x.width <- diff(x.range)
    y.width <- diff(y.range)
    x.buffer <- edge.buffer * x.width
    y.buffer <- edge.buffer * y.width

    if (strategy == "custom") {
        if (is.null(custom.centers))
            stop("custom.centers must be provided when strategy is 'custom'")
        if (!is.matrix(custom.centers) || nrow(custom.centers) != n.components ||
            ncol(custom.centers) != 2)
            stop("custom.centers must be a matrix with n.components rows and 2 columns")
        return(custom.centers)
    } else if (strategy == "random") {
        centers <- cbind(
            runif(n.components, x.range[1] + x.buffer, x.range[2] - x.buffer),
            runif(n.components, y.range[1] + y.buffer, y.range[2] - y.buffer)
        )
    } else if (strategy == "grid") {
        # Determine grid dimensions
        n.x <- ceiling(sqrt(n.components))
        n.y <- ceiling(n.components / n.x)

        # Create grid
        x.seq <- seq(x.range[1] + x.buffer, x.range[2] - x.buffer, length.out = n.x)
        y.seq <- seq(y.range[1] + y.buffer, y.range[2] - y.buffer, length.out = n.y)

        grid <- expand.grid(x = x.seq, y = y.seq)
        centers <- as.matrix(grid[1:n.components, ])
    }

    colnames(centers) <- c("x", "y")
    return(centers)
}

#' Generate amplitudes for 2D mixture
#' @keywords internal
generate_amplitudes_2d <- function(strategy, n.components, amplitude.range,
                                   custom.amplitudes) {
    if (strategy == "custom") {
        if (is.null(custom.amplitudes))
            stop("custom.amplitudes must be provided when strategy is 'custom'")
        if (length(custom.amplitudes) != n.components)
            stop("custom.amplitudes must have length n.components")
        return(custom.amplitudes)
    } else if (strategy == "positive") {
        return(runif(n.components, amplitude.range[1], amplitude.range[2]))
    } else if (strategy == "mixed") {
        amplitudes <- runif(n.components, amplitude.range[1], amplitude.range[2])
        signs <- sample(c(-1, 1), n.components, replace = TRUE)
        return(amplitudes * signs)
    }
}

#' Generate covariance matrices for 2D mixture
#' @keywords internal
generate_covariances_2d <- function(strategy, n.components, sd.value, sd.range,
                                    custom.sds, correlation) {
    covariances <- vector("list", n.components)

    if (strategy == "custom" && !is.null(custom.sds)) {
        if (is.matrix(custom.sds)) {
            # Single covariance matrix for all components
            if (nrow(custom.sds) != 2 || ncol(custom.sds) != 2)
                stop("custom.sds as matrix must be 2x2")
            for (i in 1:n.components) {
                covariances[[i]] <- custom.sds
            }
        } else if (is.vector(custom.sds)) {
            # Vector of standard deviations
            if (length(custom.sds) != n.components)
                stop("custom.sds as vector must have length n.components")
            for (i in 1:n.components) {
                sd <- custom.sds[i]
                cov.xy <- correlation * sd * sd
                covariances[[i]] <- matrix(c(sd^2, cov.xy, cov.xy, sd^2), 2, 2)
            }
        }
    } else if (strategy == "fixed") {
        cov.xy <- correlation * sd.value * sd.value
        cov.matrix <- matrix(c(sd.value^2, cov.xy, cov.xy, sd.value^2), 2, 2)
        for (i in 1:n.components) {
            covariances[[i]] <- cov.matrix
        }
    } else if (strategy == "random") {
        for (i in 1:n.components) {
            sd <- runif(1, sd.range[1], sd.range[2])
            cov.xy <- correlation * sd * sd
            covariances[[i]] <- matrix(c(sd^2, cov.xy, cov.xy, sd^2), 2, 2)
        }
    }

    return(covariances)
}

#' Plot method for gaussian_mixture objects
#'
#' Creates a 2x2 layout with contour plot, 3D surface, heat map, and filled contours.
#'
#' @param x A gaussian_mixture object
#' @param n.grid Number of grid points in each dimension. Default is 50.
#' @param main Main title for the plots. Default is "2D Gaussian Mixture"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisible NULL
#' @importFrom graphics axis contour image par persp points filled.contour
#' @importFrom grDevices heat.colors terrain.colors
#' @exportS3Method
plot.gaussian_mixture <- function(x, n.grid = 50, main = "2D Gaussian Mixture", ...) {
    # Generate evaluation grid
    x.seq <- seq(x$x.range[1], x$x.range[2], length.out = n.grid)
    y.seq <- seq(x$y.range[1], x$y.range[2], length.out = n.grid)

    # Evaluate function on grid
    z <- outer(x.seq, y.seq, Vectorize(function(xi, yi) x$f(xi, yi)))

    # Set up 2x2 layout
    old.par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old.par))

    # 1. Contour plot
    contour(x.seq, y.seq, z, main = paste(main, "- Contours"),
            xlab = "x", ylab = "y", nlevels = 15)
    points(x$centers[,1], x$centers[,2], pch = 19, col = "red", cex = 1.5)

    # 2. 3D perspective plot
    persp(x.seq, y.seq, z, theta = 45, phi = 30,
          main = paste(main, "- 3D Surface"),
          xlab = "x", ylab = "y", zlab = "f(x,y)",
          col = "lightblue", shade = 0.5)

    # 3. Heat map
    image(x.seq, y.seq, z, main = paste(main, "- Heat Map"),
          xlab = "x", ylab = "y", col = heat.colors(50))
    contour(x.seq, y.seq, z, add = TRUE, nlevels = 10)
    points(x$centers[,1], x$centers[,2], pch = 19, col = "black", cex = 1.5)

    # 4. Filled contour plot
    filled.contour(x.seq, y.seq, z,
                   main = paste(main, "- Filled Contours"),
                   xlab = "x", ylab = "y",
                   color.palette = terrain.colors,
                   plot.axes = {
                       graphics::axis(1)
                       graphics::axis(2)
                       graphics::points(x$centers[,1], x$centers[,2],
                                       pch = 19, col = "black", cex = 1.5)
                   })

    invisible(NULL)
}

#' Create 3D interactive plot of gaussian_mixture
#'
#' @param x A gaussian_mixture object
#' @param n.grid Number of grid points in each dimension
#' @param col Surface color
#' @param alpha Surface transparency
#' @param add.contours Whether to add contour lines at the base
#' @param ... Additional arguments passed to rgl functions
#'
#' @return Invisible surface object
#' @export
plot3d.gaussian_mixture <- function(x,
                                   n.grid = 50,
                                   col = "skyblue",
                                   alpha = 0.7,
                                   add.contours = TRUE,
                                   ...) {
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package 'rgl' is required for 3D plotting. Please install it.")
    }

    # Generate evaluation grid
    x.seq <- seq(x$x.range[1], x$x.range[2], length.out = n.grid)
    y.seq <- seq(x$y.range[1], x$y.range[2], length.out = n.grid)

    # Evaluate function on grid
    z <- matrix(0, n.grid, n.grid)
    for (i in 1:n.grid) {
        for (j in 1:n.grid) {
            z[i, j] <- x$f(x.seq[i], y.seq[j])
        }
    }

    # Clear previous plot
    rgl::open3d()

    # Add contours at base if requested
    if (add.contours) {
        base.z <- min(z) - 0.1 * diff(range(z))
        rgl::surface3d(x.seq, y.seq, matrix(base.z, n.grid, n.grid),
                      color = "white", alpha = 0.3)

        contour.lines <- grDevices::contourLines(x.seq, y.seq, z, nlevels = 10)
        for (line in contour.lines) {
            rgl::lines3d(line$x, line$y, rep(base.z, length(line$x)),
                        col = "darkgray", lwd = 2)
        }
    }

    # Add main surface
    surf <- rgl::surface3d(x.seq, y.seq, z, color = col, alpha = alpha, ...)

    # Add axes
    rgl::axes3d()
    rgl::title3d(xlab = "x", ylab = "y", zlab = "f(x,y)")

    # Add center points as spheres
    for (i in 1:x$n.components) {
        rgl::spheres3d(x$centers[i, 1], x$centers[i, 2],
                      x$f(x$centers[i, 1], x$centers[i, 2]),
                      radius = 0.02, col = "red")
    }

    # Set viewing angle
    rgl::view3d(theta = 45, phi = 30)

    invisible(surf)
}

#' Print method for gaussian_mixture objects
#'
#' @param x A gaussian_mixture object
#' @param ... Additional arguments (unused)
#'
#' @return Invisible x
#' @exportS3Method
print.gaussian_mixture <- function(x, ...) {
    cat("2D Gaussian Mixture\n")
    cat("===================\n")
    cat("Number of components:", x$n.components, "\n")
    cat("X range: [", x$x.range[1], ",", x$x.range[2], "]\n")
    cat("Y range: [", x$y.range[1], ",", x$y.range[2], "]\n")
    cat("\nComponent details:\n")
    for (i in 1:x$n.components) {
        cat(sprintf("  Component %d:\n", i))
        cat(sprintf("    Center: (%.3f, %.3f)\n", x$centers[i,1], x$centers[i,2]))
        cat(sprintf("    Amplitude: %.3f\n", x$amplitudes[i]))
        eig <- eigen(x$covariances[[i]])$values
        cat(sprintf("    Std devs: %.3f, %.3f\n", sqrt(eig[1]), sqrt(eig[2])))
    }
    cat("\nUse plot() for visualization\n")
    invisible(x)
}

#' Summary method for gaussian_mixture objects
#'
#' @param object A gaussian_mixture object
#' @param ... Additional arguments (unused)
#'
#' @return A summary list
#' @exportS3Method
summary.gaussian_mixture <- function(object, ...) {
    structure(
        list(
            n.components = object$n.components,
            x.range = object$x.range,
            y.range = object$y.range,
            centers = object$centers,
            amplitudes = object$amplitudes,
            total.amplitude = sum(abs(object$amplitudes)),
            parameters = object$parameters
        ),
        class = "summary.gaussian_mixture"
    )
}

#' Print summary of gaussian_mixture
#'
#' @param x A summary.gaussian_mixture object
#' @param ... Additional arguments (unused)
#'
#' @return Invisible x
#' @exportS3Method
print.summary.gaussian_mixture <- function(x, ...) {
    cat("Summary of 2D Gaussian Mixture\n")
    cat("==============================\n")
    cat("Components:", x$n.components, "\n")
    cat("Total absolute amplitude:", round(x$total.amplitude, 3), "\n")
    cat("X range: [", x$x.range[1], ",", x$x.range[2], "]\n")
    cat("Y range: [", x$y.range[1], ",", x$y.range[2], "]\n")
    cat("\nGeneration parameters:\n")
    cat("  Centers strategy:", x$parameters$centers.strategy, "\n")
    cat("  Amplitudes strategy:", x$parameters$amplitudes.strategy, "\n")
    cat("  SD strategy:", x$parameters$sd.strategy, "\n")
    invisible(x)
}

#' Sample and evaluate points from a Gaussian mixture
#'
#' Generates a dataset by sampling points from a domain and evaluating
#' a Gaussian mixture function at those points, with optional noise.
#'
#' @param mixture A gaussian_mixture object from get.gaussian.mixture.2d()
#' @param n Number of points to sample
#' @param sampling.method Method for sampling: "random", "grid", or "stratified"
#' @param noise.sd Standard deviation of additive Gaussian noise
#' @param seed Random seed for reproducibility (default: NULL)
#' @param ... Additional parameters
#'
#' @return A list of class "gaussian_mixture_data" containing:
#'   \item{X}{Matrix of sampled coordinates (n x 2)}
#'   \item{y}{Vector of noisy function values}
#'   \item{y.true}{Vector of true function values}
#'   \item{mixture}{The original gaussian_mixture object}
#'
#' @export
sample.gaussian_mixture <- function(mixture,
                                   n = 500,
                                   sampling.method = c("random", "grid", "stratified"),
                                   noise.sd = 0.1,
                                   seed = NULL,
                                   ...) {

    if (!inherits(mixture, "gaussian_mixture")) {
        stop("Input must be a gaussian_mixture object")
    }

    # Set seed if provided
    if (!is.null(seed)) set.seed(seed)

    sampling.method <- match.arg(sampling.method)

    # Generate sampling points based on method
    X <- generate_sampling_points(n, mixture$x.range, mixture$y.range, sampling.method)

    # Evaluate the mixture at these points
    y.true <- apply(X, 1, function(pt) mixture$f(pt[1], pt[2]))

    # Add noise
    y <- y.true + rnorm(length(y.true), 0, noise.sd)

    # Return structured result
    result <- list(
        X = X,
        y = y,
        y.true = y.true,
        mixture = mixture,
        noise.sd = noise.sd,
        sampling.method = sampling.method
    )

    class(result) <- c("gaussian_mixture_data", "list")
    return(result)
}

#' Generate sampling points using different methods
#'
#' @param n Number of points to generate
#' @param x.range Numeric vector of length 2 giving x-axis bounds
#' @param y.range Numeric vector of length 2 giving y-axis bounds
#' @param method Sampling method: "random", "grid", or "stratified"
#' @return Matrix of coordinates (n x 2)
#' @keywords internal
generate_sampling_points <- function(n, x.range, y.range, method = c("random", "grid", "stratified")) {
    method <- match.arg(method)

    if (method == "random") {
        # Random uniform sampling
        X <- matrix(nrow = n, ncol = 2)
        X[, 1] <- runif(n, x.range[1], x.range[2])
        X[, 2] <- runif(n, y.range[1], y.range[2])

    } else if (method == "grid") {
        # Regular grid sampling
        side_length <- ceiling(sqrt(n))
        x_seq <- seq(x.range[1], x.range[2], length.out = side_length)
        y_seq <- seq(y.range[1], y.range[2], length.out = side_length)
        full_grid <- expand.grid(x = x_seq, y = y_seq)

        # If n is less than full grid size, sample from grid
        if (n < nrow(full_grid)) {
            idx <- sample(nrow(full_grid), n)
            X <- as.matrix(full_grid[idx, ])
        } else {
            X <- as.matrix(full_grid[1:n, ])
        }

    } else if (method == "stratified") {
        # Stratified random sampling
        cells_per_side <- ceiling(sqrt(n))
        X <- matrix(nrow = n, ncol = 2)

        cell_width <- (x.range[2] - x.range[1]) / cells_per_side
        cell_height <- (y.range[2] - y.range[1]) / cells_per_side

        count <- 1
        for (i in 1:cells_per_side) {
            for (j in 1:cells_per_side) {
                if (count <= n) {
                    X[count, 1] <- x.range[1] + (i - 1 + runif(1)) * cell_width
                    X[count, 2] <- y.range[1] + (j - 1 + runif(1)) * cell_height
                    count <- count + 1
                }
            }
        }

        # Fill any remaining points randomly
        if (count <= n) {
            remaining <- n - count + 1
            X[count:n, 1] <- runif(remaining, x.range[1], x.range[2])
            X[count:n, 2] <- runif(remaining, y.range[1], y.range[2])
        }
    }

    return(X)
}

#' Plot Method for Gaussian Mixture Data Objects
#'
#' Provides various visualization options for data sampled from Gaussian mixture models.
#' Designed to work with objects created by \code{\link{sample.gaussian_mixture}}.
#'
#' @param x A gaussian_mixture_data object created by \code{\link{sample.gaussian_mixture}}
#' @param type Character string specifying the type of plot:
#'   \describe{
#'     \item{\code{"scatter2d"}}{2D scatter plot colored by function values (default)}
#'     \item{\code{"scatter3d"}}{3D scatter plot showing (x, y, z) points}
#'     \item{\code{"surface"}}{Interactive 3D surface with data points (requires \code{rgl})}
#'     \item{\code{"surface3d"}}{Static 3D surface using base graphics}
#'     \item{\code{"contour"}}{Contour plot of the underlying mixture function}
#'     \item{\code{"heatmap"}}{Heatmap visualization of the mixture function}
#'     \item{\code{"components"}}{Individual Gaussian components and full mixture}
#'     \item{\code{"nerve"}}{2D scatter with nerve complex overlay (requires \code{nerve.complex})}
#'     \item{\code{"compare"}}{Compare true values, noisy data, and estimates}
#'   }
#' @param nerve.complex Optional nerve complex object. Required only for \code{type = "nerve"}.
#' @param y.estimate Optional vector of estimated y values for comparison plots.
#'   Must have same length as data points. Used with \code{type = "compare"}.
#' @param grid_size Integer specifying grid resolution for function plots
#'   (contour, heatmap, surface). Default is 50.
#' @param knn_edges Logical. For \code{type = "scatter3d"}, whether to create
#'   k-nearest neighbor edges between points. If \code{TRUE}, connects each point
#'   to its k nearest neighbors, creating a graph structure that can help visualize
#'   local data relationships. If \code{FALSE} (default), points are displayed
#'   without any connecting edges. Requires the \code{FNN} package when \code{TRUE}.
#' @param k_neighbors Integer. For \code{type = "scatter3d"} with \code{knn_edges = TRUE},
#'   specifies the number of nearest neighbors to connect for each point. Default is 5.
#'   The actual number used will be the minimum of this value and n-1, where n is
#'   the number of data points. Larger values create denser graphs that may obscure
#'   the 3D structure, while smaller values show only the most local connections.
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the gaussian_mixture_data object
#'
#' @details
#' This method visualizes data sampled from Gaussian mixture models. The input object
#' must contain:
#' \itemize{
#'   \item \code{X}: Matrix of sampled points (n x 2)
#'   \item \code{y}: Vector of function values (possibly with noise)
#'   \item \code{y.true}: Vector of true function values (optional)
#'   \item \code{mixture}: The underlying gaussian_mixture object
#' }
#'
#' The plotting domain is automatically determined from the mixture object's
#' \code{x.range} and \code{y.range}.
#'
#' @note
#' The \code{"components"} plot type requires that the mixture object contains
#' component information (centers, amplitudes, covariances). The \code{"compare"}
#' plot type works best when \code{y.true} is available in the data object.
#'
#' @examples
#' \dontrun{
#' # Create a 2D Gaussian mixture
#' gm <- get.gaussian.mixture.2d(n.components = 3, sd.value = 0.1)
#'
#' # Sample points from the mixture
#' data <- sample(gm, n = 500, sampling.method = "random", noise.sd = 0.05)
#'
#' # Basic 2D scatter plot
#' plot(data)
#'
#' # Contour plot of underlying function
#' plot(data, type = "contour", nlevels = 20)
#'
#' # View individual components
#' plot(data, type = "components")
#'
#' # Compare true vs noisy values
#' plot(data, type = "compare")
#'
#' # 3D visualization (requires rgl)
#' plot(data, type = "surface", grid_size = 75)
#' }
#'
#' @seealso
#' \code{\link{sample.gaussian_mixture}} for creating data objects,
#' \code{\link{get.gaussian.mixture.2d}} for creating 2D mixtures
#'
#' @exportS3Method
plot.gaussian_mixture_data <- function(x,
                                      type = c("scatter2d", "scatter3d", "surface",
                                               "surface3d", "contour", "heatmap",
                                               "components", "nerve", "compare"),
                                      nerve.complex = NULL,
                                      y.estimate = NULL,
                                      grid_size = 50,
                                      knn_edges = FALSE,
                                      k_neighbors = 5,
                                      ...) {
    type <- match.arg(type)

    # Extract domain from mixture object
    x_range <- x$mixture$x.range
    y_range <- x$mixture$y.range

    switch(type,
           "scatter2d" = {
               # Simple 2D scatter plot colored by function values
               col_palette <- grDevices::heat.colors(50)
               color_breaks <- cut(x$y, 50, labels = FALSE)
               plot(x$X, col = col_palette[color_breaks], pch = 16,
                    xlim = x_range, ylim = y_range, ...)
           },

           "scatter3d" = {
               n_points <- nrow(x$X)

               if (knn_edges && requireNamespace("FNN", quietly = TRUE)) {
                   ## Create k-nearest neighbor graph
                   k <- min(k_neighbors, n_points - 1)
                   knn_result <- FNN::get.knn(x$X, k = k)

                   ## Build adjacency list
                   adj_list <- lapply(1:n_points, function(i) {
                       as.integer(knn_result$nn.index[i,])
                   })

                   ## Build weight list (using distances)
                   weight_list <- lapply(1:n_points, function(i) {
                       knn_result$nn.dist[i,]
                   })
               } else {
                   ## No edges
                   if (knn_edges && !requireNamespace("FNN", quietly = TRUE)) {
                       warning("FNN package required for knn_edges. Plotting without edges.")
                   }
                   adj_list <- lapply(1:n_points, function(i) integer(0))
                   weight_list <- lapply(1:n_points, function(i) numeric(0))
               }

               ## Create ggraph object
               g <- ggraph(adj_list, weight_list)

               ## Use plot.graph.3d
               plot.graph.3d(
                   x = g,
                   z = x$y,
                   layout = x$X,
                   z.point.size = 0.5,
                   base.plane = knn_edges,  # Show base only if edges exist
                   conn.points = knn_edges,  # Connect only if edges exist
                   vertical.lines = knn_edges,
                   ...
               )
           },

           "surface" = {
               # Interactive 3D surface with data points
               if (requireNamespace("rgl", quietly = TRUE)) {
                   # Create grid for the surface
                   x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
                   y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)

                   # Evaluate mixture on grid
                   z <- matrix(0, grid_size, grid_size)
                   for (i in 1:grid_size) {
                       for (j in 1:grid_size) {
                           z[i, j] <- x$mixture$f(x_seq[i], y_seq[j])
                       }
                   }

                   # Clear any existing rgl window
                   rgl::open3d()

                   # Plot surface
                   rgl::surface3d(x_seq, y_seq, z, col = "lightblue", alpha = 0.7)

                   # Add data points
                   rgl::points3d(x$X[,1], x$X[,2], x$y, col = "red", size = 3)

                   # Add axes
                   rgl::axes3d()
                   rgl::title3d(xlab = "X1", ylab = "X2", zlab = "y")
               } else {
                   warning("The rgl package is required for interactive 3D surface plots")
               }
           },

           "surface3d" = {
               # Static 3D surface using base graphics
               x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
               y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)

               # Evaluate mixture on grid
               z <- outer(x_seq, y_seq, Vectorize(function(xi, yi) x$mixture$f(xi, yi)))

               # Create perspective plot
               persp(x_seq, y_seq, z,
                     theta = 45, phi = 30,
                     xlab = "X1", ylab = "X2", zlab = "f(x,y)",
                     col = "lightblue", shade = 0.5,
                     main = "Gaussian Mixture Surface", ...)
           },

           "contour" = {
               # Contour plot of the mixture function
               x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
               y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)

               # Evaluate mixture on grid
               z <- outer(x_seq, y_seq, Vectorize(function(xi, yi) x$mixture$f(xi, yi)))

               # Create contour plot
               contour(x_seq, y_seq, z,
                       xlab = "X1", ylab = "X2",
                       main = "Mixture Function Contours", ...)

               # Add data points
               points(x$X, pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.3))
           },

           "heatmap" = {
               # Heatmap of the mixture function
               x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
               y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)

               # Evaluate mixture on grid
               z <- outer(x_seq, y_seq, Vectorize(function(xi, yi) x$mixture$f(xi, yi)))

               # Create heatmap
               image(x_seq, y_seq, z,
                     xlab = "X1", ylab = "X2",
                     main = "Mixture Function Heatmap",
                     col = heat.colors(50), ...)

               # Add contour lines
               contour(x_seq, y_seq, z, add = TRUE, col = "gray30")

               # Add data points
               points(x$X, pch = 19, cex = 0.3, col = "black")
           },

           "components" = {
               # Plot individual Gaussian components
               if (!is.null(x$mixture$centers) && !is.null(x$mixture$amplitudes) &&
                   !is.null(x$mixture$covariances)) {

                   x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
                   y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)

                   n_components <- x$mixture$n.components
                   n_rows <- ceiling(sqrt(n_components + 1))

                   # Set up plot layout
                   old_par <- par(no.readonly = TRUE)
                   on.exit(par(old_par))
                   par(mfrow = c(n_rows, n_rows))

                   # Plot each component
                   for (i in 1:n_components) {
                       # Create single-component function
                       component_f <- function(x_pt, y_pt) {
                           center <- x$mixture$centers[i, ]
                           amp <- x$mixture$amplitudes[i]
                           cov.inv <- solve(x$mixture$covariances[[i]])
                           diff <- c(x_pt - center[1], y_pt - center[2])
                           exponent <- -0.5 * sum(diff * (cov.inv %*% diff))
                           amp * exp(exponent)
                       }

                       z_comp <- outer(x_seq, y_seq, Vectorize(component_f))
                       image(x_seq, y_seq, z_comp,
                             main = paste("Component", i),
                             xlab = "X1", ylab = "X2",
                             col = heat.colors(50), ...)
                   }

                   # Plot full mixture
                   z_full <- outer(x_seq, y_seq,
                                   Vectorize(function(xi, yi) x$mixture$f(xi, yi)))
                   image(x_seq, y_seq, z_full,
                         main = "Full Mixture",
                         xlab = "X1", ylab = "X2",
                         col = heat.colors(50), ...)
               } else {
                   warning("Component information not available in mixture object")
               }
           },

           "nerve" = {
               # Plot with nerve complex edges
               if (is.null(nerve.complex)) {
                   stop("nerve.complex parameter is required for type='nerve'")
               }

               # Plot the data points
               col_palette <- grDevices::heat.colors(50)
               color_breaks <- cut(x$y, 50, labels = FALSE)
               plot(x$X, col = col_palette[color_breaks], pch = 16,
                    xlim = x_range, ylim = y_range,
                    xlab = "X1", ylab = "X2",
                    main = "Data with Nerve Complex", ...)

               # Add the nerve complex edges if igraph is available
               if (requireNamespace("igraph", quietly = TRUE)) {
                   # Extract edges from the nerve complex
                   edges <- igraph::as_edgelist(nerve.complex$skeleton)
                   for (i in 1:nrow(edges)) {
                       v1 <- edges[i, 1]
                       v2 <- edges[i, 2]
                       lines(rbind(x$X[v1,], x$X[v2,]), col = "black", lwd = 0.5)
                   }
               } else {
                   warning("igraph package required to plot nerve complex edges")
               }
           },

           "compare" = {
               # Compare true values, noisy data, and estimates
               old_par <- par(no.readonly = TRUE)
               on.exit(par(old_par))

               # Determine layout based on whether estimates are provided
               if (!is.null(y.estimate)) {
                   par(mfrow = c(2, 2))
               } else {
                   par(mfrow = c(1, 2))
               }

               # Use y.true if available, otherwise use y
               y_true <- if (!is.null(x$y.true)) x$y.true else x$y

               # Common color palette
               col_palette <- grDevices::heat.colors(50)

               # Plot 1: True smooth function
               if (!is.null(x$y.true)) {
                   color_breaks <- cut(y_true, 50, labels = FALSE)
                   plot(x$X, col = col_palette[color_breaks], pch = 16,
                        xlim = x_range, ylim = y_range,
                        main = "True Function",
                        xlab = "X1", ylab = "X2", ...)
               }

               # Plot 2: Noisy data
               color_breaks <- cut(x$y, 50, labels = FALSE)
               plot(x$X, col = col_palette[color_breaks], pch = 16,
                    xlim = x_range, ylim = y_range,
                    main = "Noisy Data",
                    xlab = "X1", ylab = "X2", ...)

               # Plot 3 & 4: Estimates and error (if provided)
               if (!is.null(y.estimate)) {
                   # Check length
                   if (length(y.estimate) != nrow(x$X)) {
                       warning("y.estimate length doesn't match data points")
                   } else {
                       # Plot estimates
                       color_breaks <- cut(y.estimate, 50, labels = FALSE)
                       plot(x$X, col = col_palette[color_breaks], pch = 16,
                            xlim = x_range, ylim = y_range,
                            main = "Estimates",
                            xlab = "X1", ylab = "X2", ...)

                       # Plot error (true - estimate)
                       error <- y_true - y.estimate
                       error_palette <- grDevices::colorRampPalette(
                           c("blue", "white", "red"))(50)
                       error_breaks <- cut(error, 50, labels = FALSE)
                       plot(x$X, col = error_palette[error_breaks], pch = 16,
                            xlim = x_range, ylim = y_range,
                            main = "Error (True - Estimate)",
                            xlab = "X1", ylab = "X2", ...)

                       # Add RMSE to title
                       rmse <- sqrt(mean(error^2))
                       graphics::title(sub = paste("RMSE:", round(rmse, 4)), line = 0.5)
                   }
               }
           }
    )

    invisible(x)
}

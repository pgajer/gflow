
#' Generate Random Variates from the Laplace Distribution
#'
#' This function generates random variates from the Laplace distribution
#' (also known as the double exponential distribution) with specified
#' location and scale parameters.
#'
#' @param n An integer specifying the number of observations to generate.
#'   Must be a positive integer.
#' @param location A numeric value specifying the location parameter \eqn{\mu}
#'   of the distribution. Default is 0.
#' @param scale A numeric value specifying the scale parameter b of the
#'   distribution. Must be positive. Default is 1.
#' @param seed An optional integer specifying the random seed. If NULL
#'   (the default), a random seed will be used.
#'
#' @return A numeric vector of length \code{n} containing the generated
#'   random variates.
#'
#' @details The probability density function of the Laplace distribution is:
#'   \deqn{f(x | \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right)}
#'   where \eqn{\mu} is the location parameter and b > 0 is the scale parameter.
#'
#' @note This function uses a C++ implementation for efficient random
#'   number generation. The underlying algorithm uses the inverse transform
#'   sampling method.
#'
#' @examples
#' # Generate 1000 random variates from the standard Laplace distribution
#' x <- rlaplace(1000)
#'
#' # Generate 500 random variates with location 2 and scale 3, using a specific seed
#' y <- rlaplace(500, location = 2, scale = 3, seed = 12345)
#'
#' @export
rlaplace <- function(n,
                     location = 0,
                     scale = 1,
                     seed = NULL) {

    if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n != round(n) || n <= 0) {
        stop("'n' must be a positive integer")
    }

    if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
        stop("'location' must be a finite numeric value")
    }

    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
        stop("'scale' must be a positive finite numeric value")
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed != round(seed)) {
            stop("'seed' must be a single integer or NULL")
        }
        # Convert to integer only if not NULL
        seed <- as.integer(seed)
    }

    .Call("S_rlaplace",
          as.integer(n),
          as.double(location),
          as.double(scale),
          seed)  # Pass seed as-is (NULL or integer vector)
}

#' Laplace Distribution Density Function
#'
#' Calculates the density function for the Laplace distribution.
#'
#' @param x Vector of quantiles.
#' @param location The location parameter \eqn{\mu}. Default is 0.
#' @param scale The scale parameter b. Must be positive. Default is 1.
#' @param log Logical; if TRUE, the log density is returned. Default is FALSE.
#'
#' @return A vector of density values (or log-density if log = TRUE).
#'
#' @examples
#' x <- seq(-5, 5, by = 0.1)
#' y <- dlaplace(x, location = 0, scale = 1)
#' plot(x, y, type = "l", main = "Laplace Density")
#'
#' @export
dlaplace <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("scale must be positive")

  log_dens <- -log(2 * scale) - abs(x - location) / scale

  if (log) return(log_dens)
  else return(exp(log_dens))
}

#' Laplace Distribution Function
#'
#' Calculates the cumulative distribution function for the Laplace distribution.
#'
#' @param q Vector of quantiles.
#' @param location The location parameter \eqn{\mu}. Default is 0.
#' @param scale The scale parameter b. Must be positive. Default is 1.
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]},
#'        otherwise, \eqn{P[X > x]}.
#' @param log.p Logical; if TRUE, probabilities p are given as \eqn{\log(p)}. Default is FALSE.
#'
#' @return A vector of cumulative probabilities.
#'
#' @examples
#' x <- seq(-5, 5, by = 0.1)
#' y <- plaplace(x, location = 0, scale = 1)
#' plot(x, y, type = "l", main = "Laplace CDF")
#'
#' @export
plaplace <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("scale must be positive")

  p <- 0.5 + 0.5 * sign(q - location) * (1 - exp(-abs(q - location) / scale))

  if (!lower.tail) p <- 1 - p

  if (log.p) return(log(p))
  else return(p)
}

#' Laplace Distribution Quantile Function
#'
#' Calculates the quantile function for the Laplace distribution.
#'
#' @param p Vector of probabilities.
#' @param location The location parameter \eqn{\mu}. Default is 0.
#' @param scale The scale parameter b. Must be positive. Default is 1.
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]},
#'        otherwise, \eqn{P[X > x]}.
#' @param log.p Logical; if TRUE, probabilities p are given as \eqn{\log(p)}. Default is FALSE.
#'
#' @return A vector of quantiles.
#'
#' @examples
#' p <- seq(0, 1, by = 0.1)
#' q <- qlaplace(p, location = 0, scale = 1)
#' plot(p, q, type = "l", main = "Laplace Quantile Function")
#'
#' @export
qlaplace <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("scale must be positive")

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  q <- location - scale * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))

  return(q)
}

#' Produces a random sample from the uniform distribution over the (K-1)-dimensional simplex
#'
#' @param K The size of vector of non-negative numbers that sum up to 1.
#'
#' @return A numeric vector of non-negative numbers that sum to 1.
runif.simplex <- function(K)
{
    lambda <- numeric(K)

    out <- .C(C_runif_simplex,
             as.integer(K),
             lambda=as.double(lambda))

    out$lambda
}

#' Uniformly Sample Points from a Unit Sphere Surface
#'
#' Generates uniform random samples from the surface of a unit sphere in (dim+1)-dimensional
#' space using the Gaussian normalization method. For example, dim=2 generates points on
#' the surface of a 3D unit ball (2-sphere).
#'
#' @param n.samples Integer > 0 indicating the number of samples to generate.
#' @param dim Integer >= 0 indicating the intrinsic dimension of the sphere.
#'            The sphere will be embedded in (dim+1)-dimensional space.
#'            E.g., dim=2 for points on a regular 3D sphere surface.
#'
#' @return A matrix of dimension \code{n.samples x (dim+1)} where each row represents
#'         a point uniformly sampled from the surface of the unit sphere. All points
#'         have unit Euclidean norm.
#'
#' @details
#' The function uses the fact that normalizing a vector of independent standard normal
#' variables results in a uniform distribution on the sphere surface.
#'
#' @references
#' Muller, M. E. (1959). A note on a method for generating points uniformly on
#' n-dimensional spheres. Communications of the ACM, 2(4), 19-20.
#'
#' @examples
#' # Generate 10 points on a 2D circle (1-sphere)
#' points_circle <- runif.sphere(10, 1)
#'
#' # Generate 100 points on a 3D sphere surface (2-sphere)
#' points_sphere <- runif.sphere(100, 2)
#'
#' @export
runif.sphere <- function(n.samples, dim) {
  ## Initialize the samples matrix
  samples <- matrix(0, nrow = n.samples, ncol = dim + 1)

  for(i in 1:n.samples) {
    ## Generate n independent standard Gaussian random variables
    point <- rnorm(dim + 1)
    ## Normalize the point to make it lie on the surface of the unit ball
    point <- point / sqrt(sum(point^2))
    samples[i,] <- point
  }

  return(samples)
}

#' Generate Uniform Random Sample from n-dimensional Torus
#'
#' @description
#' Generates points uniformly distributed on an n-dimensional torus. Each point on the
#' torus is represented by its embedding in 2n-dimensional Euclidean space, where each
#' pair of coordinates represents a circle in a different dimension.
#'
#' @param n.pts numeric; Number of points to generate
#' @param dim numeric; Dimension of the torus (default = 1, which gives points on a circle)
#'
#' @return A matrix with n.pts rows and 2*dim columns. Each row represents a point
#' on the torus, with pairs of columns representing the (x,y) coordinates of each
#' circular component.
#'
#' @details
#' The n-dimensional torus is the product of n circles. For each dimension, a random
#' angle is generated uniformly in \code{[0, 2pi]}, and its sine and cosine give the
#' coordinates of the point's projection onto that circular component. The resulting
#' points are uniformly distributed with respect to the natural (product) measure
#' on the torus.
#'
#' @examples
#' # Generate 100 points on a circle (1-dimensional torus)
#' circle_points <- runif.torus(100, dim = 1)
#'
#' # Generate 100 points on a 2-dimensional torus
#' torus_points <- runif.torus(100, dim = 2)
#'
#' @export
runif.torus <- function(n.pts, dim = 1) {
    # Input validation
    if(!is.numeric(n.pts) || !is.numeric(dim)) {
        stop("Both n.pts and dim must be numeric")
    }
    if(n.pts <= 0 || dim <= 0) {
        stop("Both n.pts and dim must be positive")
    }
    if(round(dim) != dim) {
        stop("dim must be a positive integer")
    }

    # Generate random angles for each dimension
    angles <- matrix(2 * pi * runif(n.pts * dim), nrow = n.pts, ncol = dim)

    # Initialize result matrix
    result <- matrix(0, nrow = n.pts, ncol = 2 * dim)

    # Fill in coordinates
    for(i in 1:dim) {
        result[, 2*i - 1] <- cos(angles[, i])  # x-coordinate for dimension i
        result[, 2*i] <- sin(angles[, i])      # y-coordinate for dimension i
    }

    # Column names for clarity
    colnames(result) <- paste0(
        rep(c("x", "y"), dim),
        rep(1:dim, each = 2)
    )

    return(result)
}

#' Sample points from empirical distribution with optional linear approximation
#'
#' This function implements the inverse transform sampling method to generate random
#' samples from an empirical distribution. The algorithm consists of the following steps:
#'
#' 1. Density Estimation:
#'    - Constructs a histogram of the input data using nbins equally spaced bins
#'    - The histogram counts represent the empirical density
#'
#' 2. PDF Normalization:
#'    - Normalizes the density values to create a proper PDF that integrates to 1
#'    - For histogram with equal bin widths dx: pdf = counts / (sum(counts) * dx)
#'
#' 3. CDF Computation:
#'    - Computes the cumulative distribution function (CDF) from the normalized PDF
#'    - For histogram with equal bin widths: \code{cdf[i] = sum(pdf[1:i]) * dx}
#'    - The CDF is a step function by nature, jumping at bin midpoints
#'
#' 4. Inverse Transform Sampling:
#'    - Generates n uniform random numbers u ~ U(0,1)
#'    - For each u, finds x such that F(x) = u, where F is the CDF
#'    - Two methods available for finding x:
#'      a) Step Function (use.linear.approximation=FALSE):
#'         - Returns the first x value where CDF \eqn{\ge} u
#'         - Preserves the discrete nature of the histogram
#'      b) Linear Interpolation (use.linear.approximation=TRUE):
#'         - Linearly interpolates between surrounding x values
#'         - Smooths the distribution but may not perfectly represent histogram
#'
#' @param n     Integer. Number of sample points to generate.
#' @param y     Numeric vector. Raw data points from which to estimate distribution.
#' @param nbins Integer. Number of bins to use when estimating density.
#'              More bins provide finer resolution but may introduce noise.
#' @param use.linear.approximation Logical. If TRUE, uses linear interpolation
#'              between CDF points. If FALSE, uses step function approach.
#'
#' @return      Numeric vector of length n containing samples from the distribution
#'              estimated from y.
#'
#' @examples
#' # Generate bimodal test data
#' y <- c(rnorm(1000, 0.3, 0.1), rnorm(1000, 0.7, 0.1))
#'
#' # Sample using step function (default)
#' samples1 <- sample.from.empirical.distribution(1000, y)
#'
#' # Sample using linear approximation
#' samples2 <- sample.from.empirical.distribution(1000, y, use.linear.approximation=TRUE)
#'
#' # Compare distributions
#' par(mfrow=c(1,3))
#' hist(y, main="Original Data")
#' hist(samples1, main="Step Function Sampling")
#' hist(samples2, main="Linear Interpolation Sampling")
#'
#' @importFrom graphics hist
#' @export
sample.from.empirical.distribution <- function(n,
                                               y,
                                               nbins = 100,
                                               use.linear.approximation = FALSE) {
    if (!is.numeric(n) || n <= 0 || n != round(n)) {
        stop("n must be a positive integer")
    }
    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }

    # Estimate density from raw data
    h <- hist(y, breaks=nbins, plot=FALSE)
    density.values <- h$counts
    x <- h$mids

    # Normalize density to create proper PDF
    dx <- diff(x[1:2])
    density.normalized <- density.values / (sum(density.values) * dx)

    # Compute CDF
    cdf <- c(0, cumsum(density.normalized) * dx)

    # Generate uniform random numbers
    u <- runif(n)

    if (use.linear.approximation) {
        # Linear interpolation between CDF points
        samples <- sapply(u, function(u.val) {
            idx <- max(which(cdf <= u.val))
            if (idx == length(x)) return(x[idx])

            alpha <- (u.val - cdf[idx]) / (cdf[idx + 1] - cdf[idx])
            x[idx] + alpha * (x[idx + 1] - x[idx])
        })
    } else {
        # Step function (no interpolation)
        samples <- x[sapply(u, function(u.val) {
            which.max(cdf >= u.val)
        })]
    }

    return(samples)
}

#' Sample points from an empirical probability distribution
#'
#' This function generates random samples from a probability distribution defined by
#' discrete density values over a grid. It uses inverse transform sampling:
#' 1. Normalizes the input density values to create a proper PDF
#' 2. Computes the cumulative distribution function (CDF)
#' 3. Generates uniform random numbers
#' 4. Uses inverse CDF method to transform uniform samples to the target distribution
#'
#' @param n     Integer. Number of sample points to generate.
#' @param y     Numeric vector. Non-negative values representing the density/histogram
#'              heights at each grid point. Does not need to be normalized.
#' @param x     Optional numeric vector. Grid points corresponding to y values.
#'              Must be strictly increasing. If not provided, a uniform grid on \code{[0,1]}
#'              is created with length(y) points.
#'
#' @return      Numeric vector of length n containing samples from the distribution
#'              defined by the empirical density y over grid x.
#'
#' @details
#' The function performs the following steps:
#' 1. If x is not provided, creates a uniform grid on \code{[0,1]}
#' 2. Normalizes y to create proper PDF by ensuring y*dx integrates to 1
#' 3. Computes CDF through cumulative trapezoidal integration
#' 4. Generates n uniform random numbers
#' 5. Uses binary search to find inverse CDF values
sample.from.1d.density <- function(n, y, x = NULL) {
    if (!is.numeric(n) || n <= 0 || n != round(n)) {
        stop("n must be a positive integer")
    }
    if (!is.numeric(y) || any(y < 0)) {
        stop("y must be a non-negative numeric vector")
    }

    # Create or validate grid
    if (is.null(x)) {
        x <- seq(0, 1, length.out = length(y))
    } else {
        if (length(x) != length(y)) {
            stop("x and y must have the same length")
        }
        if (!all(diff(x) > 0)) {
            stop("x must be strictly increasing")
        }
    }

    # Compute grid spacing (can be non-uniform)
    dx <- diff(x)

    # Normalize y to create proper PDF
    # For non-uniform grid, need to account for dx in normalization
    y.normalized <- y / sum(y[-length(y)] * dx)

    # Compute CDF through cumulative trapezoidal integration
    cdf <- c(0, cumsum(y.normalized[-length(y)] * dx))

    # Generate uniform random numbers
    u <- runif(n)

    # Find inverse CDF values through binary search
    samples <- sapply(u, function(u.val) {
        idx <- max(which(cdf <= u.val))
        if (idx == length(x)) return(x[idx])

        # Linear interpolation between grid points
        alpha <- (u.val - cdf[idx]) / (cdf[idx + 1] - cdf[idx])
        x[idx] + alpha * (x[idx + 1] - x[idx])
    })

    return(samples)
}




#' Estimates the mode of a numeric vector
#'
#' @param x A numeric vector.
#' @param min.size The minimal number of finite values within x. If length(x) < min.size, then the result is NA.
#' @param ... Parameters to be passed to density()
#' @export
mode.1D <- function(x, min.size=10, ...) {
    res <- NA
    x <- x[is.finite(x)]
    if ( length(x) > min.size )
    {
        d <- density(x, ...)
        res <- d$x[which.max(d$y)]
    }

    res
}

#' Entropy of a vector of non-negative values
#'
#' @param x    A vector of non-negative values
#' @param base The base of the logarithm in the entropy formula.
#'
#' @return   entropy with log 'base'
#' @export
entropy <- function(x, base = 10) {
    x <- x[x>0]
    x <- x / sum(x)
    -sum( x * log2(x) ) / log2(base)
}

#' Computes Shannon Evenness of a discrete probability distribution
#'
#' Computes normalized Shannon entropy (evenness) of a numeric vector defined as entropy(x) / logB(length(x)), where B is the base of the logarithm.
#'
#' @param x Numeric vector of counts or abundances
#' @param base Logarithm base (default: 2)
#' @return A number in \eqn{[0,1]} quantifying evenness (1 = perfectly even)
#' @export
evenness <- function(x, base = 2) {
  x <- x[x > 0]
  if (length(x) <= 1) return(0)

  p <- x / sum(x)
  H <- -sum(p * log(p, base = base))
  C <- log(length(p), base = base)

  E <- H / C
  return(E)
}


#' Min-Max Normalization
#'
#' Scales a numeric vector to a specified range \eqn{[\text{y.min}, \text{y.max}]}.
#' This transformation shifts and rescales the data, so the minimum value
#' corresponds to y.min and the maximum value to y.max.
#'
#' @param x A numeric vector to be normalized.
#' @param y.min The minimum value in the normalized range (default is 0).
#' @param y.max The maximum value in the normalized range (default is 1).
#'
#' @return A numeric vector where the original values of \code{x} have been
#' scaled to fall within the range \code{[y.min, y.max]}.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' minmax.normalize(x, 0, 1)  # Normalizes x to the range \eqn{[0, 1]}
#'
#' @export
minmax.normalize <- function(x, y.min = 0, y.max = 1) {
    if (!is.numeric(x)) {
        stop("Input must be a numeric vector.")
    }

    x.min <- min(x)
    x.max <- max(x)

    if (x.min == x.max) {
        return(rep(y.min, length(x)))
    }

    s <- (y.max - y.min) / (x.max - x.min)
    y <- s * (x - x.min) + y.min

    y
}

#' Inverse logit transformation x -> 1 / (1 + exp(-x)) from real numbers into the unit interval \((0,1)\)
#'
#' @param x  A numerical vector.
#'
#' @return Inverse logit transformed data.
#' @export
inv.logit <- function(x) {
    1 / (1 + exp(-x))
}

#' Normalize and Apply Inverse Logit Transformation
#'
#' First performs Min-Max normalization on a numeric vector to scale it within a specified range,
#' and then applies the inverse logit function to each element.
#'
#' The Min-Max normalization scales the data linearly between `y.min` and `y.max`. The inverse
#' logit function (logistic function) then transforms these normalized values into probabilities,
#' mapping any real-valued number into the range (0, 1).
#'
#' @param x A numeric vector to be transformed.
#' @param y.min The minimum value of the range for Min-Max normalization (default is -3).
#' @param y.max The maximum value of the range for Min-Max normalization (default is 3).
#'
#' @return A numeric vector where each original value in \code{x} has been first normalized
#' and then transformed via the inverse logit function.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(10)
#' normalize.and.inv.logit(x)  # Normalizes x and applies inverse logit
#' }
#' @export
normalize.and.inv.logit <- function(x, y.min = -3, y.max = 3) {
    y <- minmax.normalize(x, y.min, y.max)
    y <- inv.logit(y)
    y
}



#' Empirical Cumulative Distribution Function (ECDF)
#'
#' Calculates the empirical cumulative distribution function (ECDF) for a given vector of data.
#'
#' The ECDF is a step function that represents the cumulative distribution of a sample of data.
#' It is defined as the proportion of observations in the sample that are less than or equal to
#' a given value. This function calculates the ECDF for each unique value in the input vector.
#'
#' @param x A numeric vector of data.
#'
#' @return A numeric vector representing the ECDF values for each unique value in the input vector.
#'         The length of the return value is the same as the number of unique values in the input vector.
#'         The names attribute of the return value is set to "ecdf".
#'
#' @note The input vector is sorted internally to compute the ECDF. If the input vector contains
#'       missing values (NA), they will be removed before computing the ECDF.
#'
#' @examples
#' x <- c(1.0, 2.0, 3.0, 4.0, 5.0)
#' result <- ecdf.cpp(x)
#' print(result)
#'
#' @export
ecdf.cpp <- function(x) {
  if (!is.numeric(x)) {
    stop("Argument 'x' must be a numeric vector.")
  }

  if ( length(x) < 2 ) {
      stop("x has to have at least two elements.")
  }

  result <- .Call(S_ecdf, x)

  return(result)
}

#' Calculates K-Nearest Neighbor Weighted Mean
#'
#' This function calculates the weighted mean of the response variable `y` for each observation in `X`
#' based on its k-nearest neighbors. The weights are determined by applying a kernel function to the
#' distances between each observation and its neighbors.
#'
#' @param X A matrix or data frame containing the predictor variables.
#' @param y A vector containing the response variable.
#' @param k An integer specifying the number of nearest neighbors to consider.
#' @param kernel A character string specifying the kernel function to use for weighting the neighbors.
#'               Default is "epanechnikov". Other options include "uniform", "triangular", "biweight",
#'               "triweight", "cosine", "gaussian", and "rank".
#' @param nn.factor A numeric value greater than 1 used to adjust the kernel bandwidth. Default is 1.01.
#'
#' @return A vector of weighted means of the response variable for each observation in `X`.
#'
#' @details The function first finds the k-nearest neighbors for each observation in `X` using the
#'          `get.knn()` function. It then calculates the distances between each observation and its
#'          neighbors and applies the specified kernel function to the distances to obtain weights.
#'          The weights are normalized to sum to 1 for each observation. Finally, the weighted mean
#'          of the response variable is calculated for each observation using its neighbors' weights.
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 2)
#' y <- rnorm(50)
#' weighted_means <- knn.weighted.mean(X, y, k = 5)
#'
#' @importFrom FNN get.knn
#'
#' @export
knn.weighted.mean <- function(X, y, k, kernel = "epanechnikov", nn.factor = 1.01) {

    nn <- get.knn(X, k = k)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist
    nn.i <- cbind(seq(nrow(X)), nn.i)
    nn.d <- cbind(rep(0,nrow(X)), nn.d)

    nn.r <- nn.factor * nn.d[,k+1]
    rw <- row.weighting(nn.d, nn.r, kernel) # applying kernel to NN distances with kernel bandwidth nn.r
    nn.w <- rw$nn.w
    nn.w <- row.TS.norm(nn.w) # Total Sum row normalization

    Ey <- sapply(seq(nrow(nn.i)), function(i) weighted.mean(y[nn.i[i,]], w = nn.w[i,]) )

    Ey
}

#' Returns the number of times a numeric vector changes the sign between
#' consecutive points
#'
#' @param x A numeric vector
n.zeros <- function(x)
{
    n <- length(x)

    n.zeros <- 0
    zeros.ii <- c()
    for ( i in 1:(n-1) )
    {
        if ( x[i]*x[i+1] < 0 )
        {
            n.zeros <- n.zeros + 1
            zeros.ii <- c(zeros.ii, i)
        }
    }

    list(n.zeros=n.zeros,
         zeros.ii=zeros.ii)
}

#' Returns a matrix with two columns. The first column is the index of x where
#' \code{x[i] = x[i+1]} and the second column is the largest k such that \code{x[i] = x[i + k]}
#'
#' @param x A numeric vector.
#' @export
persistent.values <- function(x)
{
    stopifnot(is.numeric(x))

    n <- length(x)

    pers.mat <- c()
    i <- 1
    while ( i < n )
    {
        if ( i+1 <= n && x[i] == x[i+1] )
        {
            k <- 1
            while ( i+k <= n && x[i] == x[i+k] )
            {
                k <- k + 1
            }
            pers.mat <- rbind(pers.mat, c(i, k-1))
            i <- i + k
        } else {
            i <- i + 1
        }
    }

    pers.mat
}

#' Returns a matrix with two columns. The first column is the index of x where
#' x\[i\] is not NA and the second column is the largest k such that x\[i\] == x\[i + k\]
#'
#' @param x A numeric vector.
#'
non.NA.values <- function(x)
{
    n <- length(x)

    start.len.mat <- c()
    i <- 1
    while ( i <= n )
    {
        if ( !is.na(x[i]) )
        {
            k <- 1
            while ( i+k <= n && !is.na(x[i+k]) && x[i] == x[i+k] )
            {
                k <- k + 1
            }
            start.len.mat <- rbind(start.len.mat, c(i, k-1))
            i <- i + k
        } else {
            i <- i + 1
        }
    }

    start.len.mat
}

#' Returns a matrix with two columns. The first column is the index of x where
#' x\[i\] is not 0 and the second column is the largest k such that x\[i\] == x\[i + k\]
#'
#' @param x A numeric vector.
#'
non.zero.values <- function(x)
{
    stopifnot(is.finite(x))

    n <- length(x)

    start.len.mat <- c()
    i <- 1
    while ( i <= n )
    {
        if ( x[i] > 0 )
        {
            k <- 1
            while ( i+k <= n && x[i+k] > 0 && x[i] == x[i+k] )
            {
                k <- k + 1
            }
            start.len.mat <- rbind(start.len.mat, c(i, k-1))
            i <- i + k
        } else {
            i <- i + 1
        }
    }

    start.len.mat
}

#' Proportion of cases in a binary variable y
#'
#' @param y  A binary variable with 0/1 values.
#' @export
p.cases <- function(y)
{
    f <- table(y)

    ret <- f["1"]/sum(f)

    ret[[1]]
}


#' Computes a cross product between two 3D vectors
#'
#' This routine calculates the cross product between two vectors in \code{R^3}.
#'
#' @param x A 3D vector.
#' @param y A 3D vector.
#'
#' @return The cross product of x and y.
#'
#' @examples
#' \dontrun{
#' x <- c(0,1,3)
#' y <- c(2,3,4)
#' z <- cross.prod(x, y)
#' }
#' @export
cross.prod <- function(x, y)
{
    create3D <- function(x) utils::head(c(x, rep(0, 3)), 3)

    x <- create3D(x)
    y <- create3D(y)

    j <- function(i) (i-1) %% 3+1

    i=1:3

    return (x[j(i+1)]*y[j(i+2)] - x[j(i+2)]*y[j(i+1)])
}

#' The norm of a vector
#'
#' @param x A numeric vector.
#' @export
vector.norm <- function(x)
{
    sqrt(sum(x^2))
}


#' Computes the angle between two 3D vectors in radians
#'
#' @param v A 3D vector.
#' @param w A 3D vector.
#' @param method The method to be used for the angle estimate. Possible choices are: "arcsine" and "arctan".
#'
#' @return The angle between two 3D vectors in radians.
#'
#' @examples
#' \dontrun{
#' v <- c(0,1,3)
#' w <- c(2,3,4)
#' angle <- angle.3D(v, w)
#' }
#' @export
## Based on the Wikipedia article:  https://en.wikipedia.org/wiki/Great-circle_distance
angle.3D <- function(v, w, method="arcsine") {
    stopifnot(length(v)==length(w))
    stopifnot(is.numeric(v), is.numeric(w))
    stopifnot(all(is.finite(v), is.finite(w)))
    stopifnot(method %in% c("arcsine", "arctan"))

    d <- NA
    if (method == "arcsine") {
        Ed <- sqrt(sum((v-w)^2)) # the Euclidean distance between v and w
        d <- 2*asin(Ed/2)
    } else {
        cp <- cross.prod(v, w)
        d <- sqrt(sum(cp^2)) / sum(v * w)
        d <- atan(d)
    }

    d
}

## ----------------------------------------------------------------------------------------------------
##
## Liear Projection
##
## ----------------------------------------------------------------------------------------------------

#' Projects a Vector onto a Subspace in R^n
#'
#' This function projects a given vector in R^n onto a subspace spanned by
#' unit vectors. The subspace is defined by a matrix where each column is a
#' unit vector spanning the subspace.
#'
#' @param x A numeric vector of length n representing the vector to be projected.
#' @param U A matrix of dimensions n x m, where each column is a unit vector
#'   spanning the subspace. The matrix should have full column rank.
#'
#' @return A numeric vector representing the projection of x onto the subspace
#'   spanned by the columns of U.
#'
#' @examples
#' u1 <- c(1, 0, 0)
#' u2 <- c(0, 1, 0)
#' u3 <- c(0, 0, 1)
#' U <- matrix(c(u1, u2, u3), nrow=3)
#' x <- c(2, 3, 4)
#' project.onto.subspace(x, U)
#'
#' @export
project.onto.subspace <- function(x, U) {

    ## Compute the projection matrix P
    P <- U %*% solve(t(U) %*% U) %*% t(U)

    ## Compute the projection of x onto the subspace
    p <- P %*% x

    return(p)
}

## ----------------------------------------------------------------------------------------------------
##
## Distance Measures
##
## ----------------------------------------------------------------------------------------------------

#' Calculates a Symmetric Density-Associated Distance
#'
#' This function computes a custom, symmetric distance metric between two
#' points `p` and `q` based on the provided density function `density.func`.
#' The Euclidean distance between `p` and `q` is normalized by the average
#' of the densities at these points, making the distance symmetric.
#'
#' Any function that is non-negative over all points of the given set can be
#' used as a density function in this function. This allows for scaling of the
#' "density function" so that it affects the Euclidean distance more and more.
#' The function transform.ED() can be used to scale any "density function".
#'
#' @param p A numeric vector representing coordinates of the first point.
#' @param q A numeric vector representing coordinates of the second point.
#' @param p.density The density value at point p.
#' @param q.density The density value at point q.
#' @param density.func A function representing the density at a given point in the space (optional, default = NULL).
#'
#' @return The calculated density-associated distance between `p` and `q`.
#' @examples
#' density.func <- function(x) { return(dnorm(x, mean=0, sd=1)) }
#' p <- c(0, 0)
#' q <- c(1, 1)
#' p.density <- density.func(p)
#' q.density <- density.func(q)
#' Rdensity.distance(p, q, p.density, q.density, density.func)
#'
#' @export
Rdensity.distance <- function(p, q, p.density, q.density, density.func = NULL) {

    ## Calculate Euclidean distance between p and q
    euclidean.distance <- sqrt(sum((p - q)^2))

    ## Calculate densities at p and q using the given density function
    ## p.density <- density.func(p)
    ## q.density <- density.func(q)

    ## Calculate the average of the densities at p and q
    mean.density <- (p.density + q.density) / 2

    ## Return the density-normalized distance
    return(euclidean.distance / mean.density)
}


## ----------------------------------------------------------------------------------------------------
##
## Numerical differentiation
##
## ----------------------------------------------------------------------------------------------------

#' Second-order accurate method of derivative of a function estimate
#'
#' @description
#' Computes the derivative using the second-order accurate finite difference method:
#' \deqn{f'(t_i) \approx \frac{-f(t_{i+2}) + 8f(t_{i+1}) - 8f(t_{i-1}) + f(t_{i-2})}{12\Delta t}}
#'
#' @details
#' It is assumed that y is defined over a uniform grid. The method uses:
#' - Forward difference for the first point
#' - Central difference for the second and second-to-last points
#' - Backward difference for the last point
#' - Five-point stencil for interior points
#'
#' @param y  A numeric vector of response values defined over a uniform grid.
#' @param dx The distance between consecutive points of the grid over which y is defined.
#'
#' @return A numeric vector of the same length as y containing the derivative estimates.
#'
#' @examples
#' \dontrun{
#' # Example with a quadratic function
#' x <- seq(0, 2*pi, length.out = 100)
#' y <- sin(x)
#' dx <- x[2] - x[1]
#' dy <- derivative.second.order.method(y, dx)
#'
#' # Compare with analytical derivative
#' dy_true <- cos(x)
#' plot(x, dy, type = "l", col = "blue", main = "Numerical vs Analytical Derivative")
#' lines(x, dy_true, col = "red", lty = 2)
#' legend("topright", c("Numerical", "Analytical"), col = c("blue", "red"), lty = c(1, 2))
#' }
#' @export
derivative.second.order.method <- function(y, dx)
{
    n <- length(y)
    dy <- numeric(n)
    dy[1] <- (y[2] - y[1]) / dx
    dy[2] <- (y[3] - y[1]) / (2 * dx)
    dy[n-1] <- (y[n] - y[n-2]) / (2 * dx)
    dy[n] <- (y[n] - y[n-1]) / dx
    ddx <- 12 * dx
    for ( i in 3:(n-2) ) {
        ## f'(t_i) approx (-f(t_{i+2}) + 8*f(t_{i+1}) - 8*f(t_{i-1}) + f(t_{i-2})) / (12*dt)
        dy[i] <- (-y[i+2] + 8*y[i+1] - 8*y[i-1] + y[i-2]) / ddx
    }
    dy
}

## ----------------------------------------------------------------------------------------------------
##
## Transformation functions
##
## ----------------------------------------------------------------------------------------------------

#' Linear transformation of 'x' values so that the min(x) is sent to ymin and max(x) to ymax
#'
#' @param x  Input values
#' @param ymin The value to which min(x) is sent.
#' @param ymax The value to which max(x) is sent.
#' @param xrange The minimum of x and the maximum of x that is to be mapped to ymin and ymax, respectively.
#' @return A numeric vector of the same length as 'x' with values linearly transformed to the range \code{[ymin, ymax]}.
#' @export
scale_to_range <- function(x, ymin, ymax, xrange = NULL) {

    if( !is.numeric(x) ) {
        stop("x has to be numeric.")
    }

    if( !is.numeric(ymin) ) {
        stop("ymin has to be numeric.")
    }

    if( !is.numeric(ymax) ) {
        stop("ymax has to be numeric.")
    }

    if (!is.null(xrange)) {
        if ( !is.numeric(xrange) || length(xrange) != 2 ) {
            stop("xrange, if provided, must be numeric and contain exactly two elements.")
        }
        xmin <- xrange[1]
        xmax <- xrange[2]
    } else {
        xmin <- min(x)
        xmax <- max(x)
    }

    if( xmin == xmax ) { # Avoid division by zero if all x values are identical
        return(rep((ymin + ymax) / 2, length(x)))
    }

    slope <- (ymax - ymin) / (xmax - xmin)

    y <- slope * (x - xmin) + ymin # Correctly applying element-wise multiplication

    return(y)
}

#' Subtracts the centroid from each point of the given dataset
#'
#' This function centers a dataset by subtracting a centroid from each point (row).
#' If no centroid is provided, it computes the centroid as the column means of X.
#'
#' @param X A matrix or data.frame specifying a set of points, where each row represents
#'   a point and each column represents a dimension.
#' @param centroid A numeric vector specifying the centroid to subtract from each point.
#'   If NULL (default), the centroid is computed as the column means of X.
#'   Must have length equal to ncol(X).
#' @return A numeric matrix of the same dimensions as X, with the centroid subtracted
#'   from each row. The returned matrix represents the centered dataset.
#' @examples
#' # Example 1: Auto-compute centroid
#' X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#' centroid.shift(X)
#'
#' # Example 2: Provide custom centroid
#' X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#' custom_centroid <- c(2, 5)
#' centroid.shift(X, centroid = custom_centroid)
#' @export
centroid.shift <- function(X, centroid = NULL) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    if ( is.null(centroid) ) {
        ## Calculate the centroid
        centroid <- colMeans(X)
    } else if ( length(centroid) != ncol(X) ) {
        stop("centroid length has to be the same as the number of columns of X")
    }

    ## Center the points by subtracting the centroid
    cX <- sweep(X, 2, centroid, FUN = "-")

    return(cX)
}

#' Maps points from the interior of a unit disk onto the unit sphere of the same dimension as the disk
#'
#' This function transforms points from the interior of an n-dimensional unit disk (points with
#' L2 norm < 1) onto the surface of an (n+1)-dimensional unit sphere using a stereographic-like
#' projection. Points at the center of the disk map to the "north pole" of the sphere.
#'
#' @param X A matrix or data frame specifying a set of points within the unit disk, where each
#'   row represents a point and each column represents a dimension. All points must have
#'   L2 norm strictly less than 1.
#' @return A numeric matrix with nrow(X) rows and ncol(X)+1 columns, where each row represents
#'   a point on the unit sphere in (n+1)-dimensional space. The returned points all have
#'   L2 norm equal to 1.
#' @details The mapping uses the following transformation:
#'   \itemize{
#'     \item For a point p with radius \eqn{r = ||p||}, the angle \eqn{\phi = \pi r}
#'     \item The sphere coordinates are: \eqn{(\sin(\phi) p/r, \cos(\phi))}
#'     \item Points at the origin (r = 0) map to (0, 0, ..., 0, 1)
#'   }
#' @examples
#' # Example 1: Map 2D disk points to 3D sphere
#' X <- matrix(c(0.5, 0.3, -0.2, 0.4, 0, 0), nrow = 3, ncol = 2, byrow = TRUE)
#' sphere_points <- disk.to.sphere(X)
#' # Verify all points are on unit sphere
#' apply(sphere_points, 1, function(x) sqrt(sum(x^2)))
#'
#' # Example 2: Center point maps to north pole
#' origin <- matrix(c(0, 0), nrow = 1)
#' disk.to.sphere(origin)  # Returns (0, 0, 1)
#' @export
disk.to.sphere <- function(X) {
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }
    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }
    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }
    if ( ncol(X) < 2 ) {
        stop("X has to have at least two columns.")
    }
    dist.to.0 <- apply(X, 1, function(x) sqrt(sum(x^2)))
    if ( !all(dist.to.0 < 1) ) {
        stop("All points of X must have L2 norm less than 1.")
    }
    sX <- matrix(nrow = nrow(X), ncol = ncol(X) + 1)
    for ( i in seq(nrow(X)) ) {
        p <- as.numeric(X[i,])
        r <- sqrt(sum(p^2))
        if ( r > 0 ) {
            x <- p / r
            phi <- pi * r
            sX[i,] <- c(sin(phi) * x, cos(phi))
        } else {
            # At origin, map to north pole: (0, 0, ..., 0, 1)
            sX[i,] <- c(rep(0, ncol(X)), 1)
        }
    }
    return(sX)
}

#' Perform Lp Radial Transformation
#'
#' This function applies the Lp radial transformation to each row of the input matrix or data frame.
#' The transformation is defined as x -> r^p * x / r, where r is the L2 norm of x.
#'
#' @param X A matrix or data frame containing numeric values to be transformed.
#' @param p The Lp exponent. Must be greater than 0. Default is 0.5.
#'
#' @return A matrix or data frame of the same dimensions as the input, with the Lp radial transformation applied to each row.
#'
#' @examples
#' X <- matrix(runif(20), ncol = 2)
#' transformed_X <- radial.Lp.transform(X, p = 0.5)
#'
#' @export
radial.Lp.transform <- function(X, p = 0.5) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    if ( p <= 0 ) {
        stop("p has to be greater than 1")
    }

    for ( i in seq(nrow(X)) ) {
        q <- as.numeric(X[i,])
        r <- sqrt(sum(q^2))
        X[i,] <- r^p * q / r
    }

    return(X)
}

#' Restrict Points to L-infinity Unit Sphere
#'
#' This function restricts the points represented by rows in the input matrix or data frame
#' to the L-infinity unit sphere. If the maximum absolute value of a row exceeds 1,
#' the row is scaled down to have a maximum absolute value of 1.
#'
#' @param X A matrix or data frame containing numeric values to be restricted.
#'
#' @return A matrix or data frame of the same dimensions as the input, with each row
#'         restricted to the L-infinity unit sphere.
#'
#' @examples
#' X <- matrix(runif(20, -2, 2), ncol = 2)
#' restricted_X <- restric.to.Linf.unit.sphere(X)
#'
#' @export
restric.to.Linf.unit.sphere <- function(X) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    for ( i in seq(nrow(X)) ) {
        q <- as.numeric(X[i,])
        m <- max(q)
        if ( m > 1 ) {
            X[i,] <- q / m
        }
    }

    return(X)
}

#' Applies a Linear Transformation and Tangent Function
#'
#' This function applies a linear transformation to the input variable 'x',
#' scaling it so that the range of the transformed 'x' falls within the interval
#' \code{[-pi/2 + offset, pi/2 - offset]}. The function then applies a tangent
#' function to the result of the transformation.
#'
#' @param x The input variable to undergo transformation.
#' @param offset A value representing the desired offset of the transformed variable range from -pi/2 and pi/2.
#'
#' @return Returns the tangent of the linearly transformed 'x'.
#'
#' @examples
#' \dontrun{
#' x <- 0.001
#' shifted.tan(x)
#' }
#' @export
shifted.tan <- function(x, offset=0.1)
{
    x0 <- min(x)
    x1 <- max(x)
    y0 <- -pi/2 + offset
    C <- (pi - 2*offset) / (x1 - x0)
    y <- C * (x - x0) + y0
    tan(y)
}


#' Applies a Scaled Log Transformation
#'
#' This function first applies a scaled log transformation, C*log(x), to the
#' input variable 'x'. The scaling is performed such that the range of the
#' transformed variable 'x' lies within the interval \code{[-pi/2, pi/2]}.
#' Subsequently, a tangent function is applied to the result of the
#' transformation.
#'
#' @param x     The input variable to undergo transformation.
#' @param scale A value representing the desired offset of the transformed variable range from -pi/2 and pi/2.
#'
#' @return Returns the tangent of the scaled log of x.
#'
#' @examples
#' \dontrun{
#' x <- 0.1
#' scaled.log.tan(x)
#' }
#' @export
scaled.log.tan <- function(x, scale=0.1)
{
    x <- scale * log(x)

    stopifnot(max(x) < pi/2)
    stopifnot(min(x) > -pi/2)

    tan(x)
}

#' Restricts the range of value of a numeric vector to \code{[xmin, xmax]}
#'
#' @param x     A numeric vector.
#' @param xmin  The minimum value of the clamped x.
#' @param xmax  The maximum value of the clamped x.
#'
#' @return A vector with all values below xmin set to xmin and all value above xmax set to xmax.
#' @export
clamp <- function(x, xmin, xmax) {

    clamped.x <- pmax(pmin(x, xmax), xmin)

  clamped.x
}

#' Returns robust log transform of a non-negative vector
#'
#' This function applies a log transformation to positive values and normalizes them by their mean,
#' while keeping zero values as zero. This creates a scale-invariant transformation that preserves
#' the relative differences between positive values.
#'
#' @param x A numeric vector of non-negative numbers.
#' @return A numeric vector of the same length as x, where positive values are log-transformed
#'   and normalized by the mean of log-transformed values, and zero values remain zero.
#' @details The transformation process:
#'   \enumerate{
#'     \item Positive values are log-transformed: \code{log(x[i])} for \code{x[i] > 0}
#'     \item The mean of log-transformed values is calculated: \code{M = mean(log(x[x > 0]))}
#'     \item Log-transformed values are divided by M for normalization
#'     \item Zero values remain unchanged as 0
#'   }
#' @note This function assumes all input values are non-negative. Negative values will cause
#'   unexpected behavior as they are not checked.
#' @examples
#' # Example 1: Basic usage
#' x <- c(0, 1, 10, 100, 1000)
#' robust.log.transform(x)
#'
#' # Example 2: Scale invariance property
#' x1 <- c(0, 1, 2, 5, 10)
#' x2 <- c(0, 10, 20, 50, 100)  # x1 * 10
#' # The relative differences are preserved
#' robust.log.transform(x1)
#' robust.log.transform(x2)
#' @export
robust.log.transform <- function(x) {
    idx <- x > 0
    x[idx] <- log(x[idx])
    M <- mean(x[idx])
    x[idx] <- x[idx] / M

    return(x)
}

#' Returns robust transform of a non-negative vector
#'
#' This function normalizes positive values by their geometric mean while keeping zero values
#' as zero. The geometric mean is more robust to outliers than the arithmetic mean, making
#' this transformation useful for data with skewed distributions.
#'
#' @param x A numeric vector of non-negative numbers.
#' @return A numeric vector of the same length as x, where positive values are divided by
#'   the geometric mean of all positive values, and zero values remain zero.
#' @details The transformation process:
#'   \enumerate{
#'     \item The geometric mean of positive values is calculated: \code{D = exp(mean(log(x[x > 0])))}
#'     \item Each positive value is divided by \code{D: x[i] / D for x[i] > 0}
#'     \item Zero values remain unchanged as 0
#'   }
#'   The geometric mean of the transformed positive values will be 1.
#' @note This function assumes all input values are non-negative. Negative values will cause
#'   unexpected behavior as they are not checked.
#' @examples
#' # Example 1: Basic usage
#' x <- c(0, 1, 4, 16, 64)
#' robust.transform(x)
#' # Geometric mean of positive transformed values equals 1
#'
#' # Example 2: Comparison with arithmetic mean normalization
#' x <- c(0, 1, 2, 100)  # Outlier present
#' # Robust (geometric mean) normalization
#' robust.transform(x)
#' # Compare to arithmetic mean normalization
#' x_arith <- x
#' x_arith[x > 0] <- x[x > 0] / mean(x[x > 0])
#' x_arith  # More affected by the outlier
#' @export
robust.transform <- function(x) {
    idx <- x > 0
    D <- exp(mean(log(x[idx])))
    x[idx] <- x[idx] / D

    return(x)
}

#' Determine Optimal Number of Principal Components Based on Variance Explained
#'
#' This helper function performs PCA and determines the optimal number of principal
#' components needed to explain a specified percentage of variance in the data.
#'
#' @param X A numeric matrix to analyze.
#' @param variance.threshold The percentage of variance to be explained (between 0 and 1).
#' @param max.components Maximum number of components to consider (default is NULL, which means all columns).
#' @param center Logical indicating whether the variables should be centered (default is TRUE).
#' @param scale Logical indicating whether the variables should be scaled (default is FALSE).
#'
#' @return A list containing:
#'   \item{n.components}{Number of components needed to meet the variance threshold}
#'   \item{variance.explained}{Percentage of variance explained by the selected components}
#'   \item{cumulative.variance}{Vector of cumulative variance explained for all components}
#'   \item{pca.result}{The full PCA result object}
#'
#' @importFrom stats prcomp
#' @keywords internal
pca.optimal.components <- function(X, variance.threshold = 0.99, max.components = NULL,
                                  center = TRUE, scale = FALSE) {
    # Set maximum components to consider
    if (is.null(max.components) || max.components > ncol(X)) {
        max.components <- ncol(X)
    }

    # Perform PCA
    pca_result <- prcomp(X, center = center, scale. = scale)

    # Calculate total variance and cumulative variance explained
    total.variance <- sum(pca_result$sdev^2)
    cumulative.variance <- cumsum(pca_result$sdev^2) / total.variance

    # Find the smallest number of components exceeding the variance threshold
    n.components <- which(cumulative.variance >= variance.threshold)[1]
    if (is.na(n.components)) {
        n.components <- length(cumulative.variance)
    }

    # Limit to max.components
    n.components <- min(n.components, max.components)

    # Return results
    list(
        n.components = n.components,
        variance.explained = cumulative.variance[n.components],
        cumulative.variance = cumulative.variance,
        pca.result = pca_result
    )
}

#' Project Data onto Principal Components
#'
#' This helper function projects a data matrix onto its first n principal components.
#'
#' @param X A numeric matrix to be projected.
#' @param pca.result A PCA result object from prcomp().
#' @param n.components The number of principal components to retain.
#'
#' @return A matrix with the same number of rows as X but with n.components columns,
#'         representing the projection of X onto its first n.components principal components.
#'         Row names from the original matrix are preserved.
#'
#' @keywords internal
pca.project <- function(X, pca.result, n.components) {
    # Store row names before projection
    row_names <- rownames(X)

    # Extract the projection onto the first n.components principal components
    projection <- pca.result$x[, 1:n.components, drop = FALSE]

    # Assign original row names to the projection if they exist
    if (!is.null(row_names)) {
        rownames(projection) <- row_names
    }

    # Return the projection
    return(projection)
}

#' Perform PCA Variance Analysis
#'
#' This function performs Principal Component Analysis (PCA) on the input data
#' and analyzes the variance explained by different numbers of principal components.
#'
#' @param X A numeric matrix or data frame with at least two columns.
#' @param pc.90 Threshold for explaining 90% of variance (default: 0.90).
#' @param pc.95 Threshold for explaining 95% of variance (default: 0.95).
#' @param pc.97.5 Threshold for explaining 97.5% of variance (default: 0.975).
#' @param pc.99 Threshold for explaining 99% of variance (default: 0.99).
#' @param max.components Maximum number of components to analyze (default: 300).
#'
#' @return A list containing:
#'   \item{variance_explained}{Vector of cumulative variance explained}
#'   \item{pc_90}{Number of PCs explaining >90% variance}
#'   \item{pc_95}{Number of PCs explaining >95% variance}
#'   \item{pc_97_5}{Number of PCs explaining >97.5% variance}
#'   \item{pc_99}{Number of PCs explaining >99% variance}
#'   \item{X_pca90}{PCA projection for 90% variance threshold}
#'   \item{X_pca95}{PCA projection for 95% variance threshold}
#'   \item{X_pca97_5}{PCA projection for 97.5% variance threshold}
#'   \item{X_pca99}{PCA projection for 99% variance threshold}
#'   \item{auc}{Area under the curve of variance explained}
#'   \item{pca_result}{Full PCA result object}
#'
#' @examples
#' data <- matrix(rnorm(1000*50), ncol=50)
#' result <- pca.variance.analysis(data)
#' print(result$pc.95)
#'
#' @export
pca.variance.analysis <- function(X, pc.90 = 0.90, pc.95 = 0.95, pc.97.5 = 0.975, pc.99 = 0.99, max.components = 300) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    if (ncol(X) < 2) {
        stop("X must have at least two columns")
    }
    if (!all(sapply(list(pc.90, pc.95, pc.97.5, pc.99), is.numeric)) ||
        !all(sapply(list(pc.90, pc.95, pc.97.5, pc.99), function(x) x > 0 && x < 1))) {
        stop("All threshold parameters must be numeric values between 0 and 1")
    }
    if (!is.numeric(max.components) || max.components <= 0 || max.components != as.integer(max.components) ) {
        stop("max.components must be a positive integer not exceeding the number of columns in X")
    }

    if (max.components > ncol(X)) {
        max.components <- ncol(X)
    }

    ## PCA
    X.pca.res <- prcomp(X)

    ## Calculate total variance
    total.variance <- sum(X.pca.res$sdev^2)

    ## Calculate cumulative variance explained
    variance.explained <- cumsum(X.pca.res$sdev^2) / total.variance

    ## Function to find the smallest number of PCs exceeding a threshold
    find.smallest.pc <- function(variance.explained, threshold) {
        which(variance.explained > threshold)[1]
    }

    ## Find the number of PCs for each threshold
    pc.90.value <- find.smallest.pc(variance.explained, pc.90)
    pc.95.value <- find.smallest.pc(variance.explained, pc.95)
    pc.97.5.value <- find.smallest.pc(variance.explained, pc.97.5)
    pc.99.value <- find.smallest.pc(variance.explained, pc.99)

    ## Calculate area under the curve
    auc <- sum(variance.explained[1:max.components])

    ## Create PCA projections
    create.pca.projection <- function(n.components) {
        pca.components <- X.pca.res$rotation[, 1:n.components]
        as.matrix(X) %*% pca.components
    }

    list(
        variance.explained = variance.explained[1:max.components],
        pc.90 = pc.90.value,
        pc.95 = pc.95.value,
        pc.97.5 = pc.97.5.value,
        pc.99 = pc.99.value,
        X.pca90 = create.pca.projection(pc.90.value),
        X.pca95 = create.pca.projection(pc.95.value),
        X.pca97.5 = create.pca.projection(pc.97.5.value),
        X.pca99 = create.pca.projection(pc.99.value),
        auc = auc,
        pca.result = X.pca.res
    )
}


#' Calculate the Overlap Coefficient Between Two Numeric Vectors
#'
#' @description
#' Computes the Overlap Coefficient (also known as Szymkiewicz-Simpson coefficient)
#' between two numeric vectors treated as sets. The coefficient is defined as the size
#' of the intersection divided by the size of the smaller set. Duplicates are removed
#' before calculation.
#'
#' @param x A numeric vector representing the first set
#' @param y A numeric vector representing the second set
#'
#' @return A numeric value between 0 and 1, where:
#'   \itemize{
#'     \item 1 indicates that one set is a subset of the other
#'     \item 0 indicates the sets are disjoint
#'     \item Values between 0 and 1 indicate partial overlap
#'   }
#'
#' @examples
#' overlap.coefficient(c(1, 2, 3), c(1, 2, 3))     # Returns 1
#' overlap.coefficient(c(1, 2), c(1, 2, 3, 4))     # Returns 1
#' overlap.coefficient(c(1, 2, 2, 3), c(1, 2, 4))  # Returns 0.67
#' overlap.coefficient(c(1, 2), c(3, 4))           # Returns 0
#' overlap.coefficient(numeric(0), c(1, 2))        # Returns 0
#'
#' @export
#'
#' @seealso
#' Other set similarity metrics:
#' \code{\link[stats]{dist}} for general distance measures
#'
#' @references
#' Simpson, G.G. (1960) Notes on the measurement of faunal resemblance.
#' American Journal of Science, 258-A: 300-311.
overlap.coefficient <- function(x, y) {
    # Input validation
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Both inputs must be numeric vectors")
    }
    if (length(x) == 0 || length(y) == 0) {
        return(0)  # Empty sets have no overlap
    }

    # Find intersection size using unique elements
    intersection.size <- length(intersect(unique(x), unique(y)))

    # Find minimum size of the two sets (after removing duplicates)
    min.size <- min(length(unique(x)), length(unique(y)))

    # Calculate overlap coefficient
    return(intersection.size / min.size)
}

#' Calculate Overlap Coefficients Between Basins of Attraction
#'
#' @description
#' Computes a matrix of overlap coefficients between basins of attraction (BoA)
#' identified by local maxima indices in two variables. Only includes BoAs that meet
#' the minimum size requirement.
#'
#' @param x Named integer vector where values are indices of local maxima and names
#'        are IDs of points in the basin for variable X
#' @param y Named integer vector where values are indices of local maxima and names
#'        are IDs of points in the basin for variable Y
#' @param min.BoA.size Minimum number of elements required for a basin to be included
#' @param x.labels Named vector of labels for basins in x. Names must match unique values in x.
#' @param y.labels Named vector of labels for basins in y. Names must match unique values in y.
#'
#' @return A matrix where rows correspond to BoAs from x and columns to BoAs from y,
#'         containing their pairwise overlap coefficients. Only includes BoAs meeting
#'         the size threshold. Row and column names are taken from x.labels and y.labels
#'         if provided, otherwise from the local maxima indices.
#'
#' @examples
#' ids <- paste0("id", 1:500)
#' x <- c(rep(1, 200), rep(2, 300))
#' names(x) <- ids
#' y <- c(rep(1, 150), rep(2, 350))
#' names(y) <- ids
#' x.labels <- c("1" = "Peak A", "2" = "Peak B")
#' y.labels <- c("1" = "Max 1", "2" = "Max 2")
#' boa.overlap(x, y, min.BoA.size = 100, x.labels = x.labels, y.labels = y.labels)
#'
#' @export
boa.overlap <- function(x, y, min.BoA.size, x.labels = NULL, y.labels = NULL) {
   # Input validation
   if (!identical(names(x), names(y))) {
       stop("x and y must have identical names")
   }
   if (!is.integer(x) && !is.numeric(x) || !is.integer(y) && !is.numeric(y)) {
       stop("x and y must be integer or numeric vectors")
   }
   if (!is.numeric(min.BoA.size) || min.BoA.size < 1) {
       stop("min.BoA.size must be a positive number")
   }

   # Validate labels if provided
   x.lmax <- as.character(unique(x))
   y.lmax <- as.character(unique(y))

   ##  if (!is.null(x.labels)) {
   ##     if (!all(x.lmax %in% names(x.labels))) {
   ##         stop("x.labels must contain names matching all unique values in x")
   ##     }
   ## }

   ## if (!is.null(y.labels)) {
   ##     if (!all(y.lmax %in% names(y.labels))) {
   ##         stop("y.labels must contain names matching all unique values in y")
   ##     }
   ## }

   # Find basins and their sizes for x
   x.basins <- lapply(x.lmax, function(lmax) names(x)[x == lmax])
   x.sizes <- sapply(x.basins, length)

   # Find basins and their sizes for y
   y.basins <- lapply(y.lmax, function(lmax) names(y)[y == lmax])
   y.sizes <- sapply(y.basins, length)

   # Filter basins by minimum size
   x.valid <- x.sizes >= min.BoA.size
   y.valid <- y.sizes >= min.BoA.size

   if (!any(x.valid) || !any(y.valid)) {
       stop("No basins meet the minimum size requirement")
   }

   x.basins <- x.basins[x.valid]
   y.basins <- y.basins[y.valid]

   # Create matrix of overlap coefficients
   n.x <- sum(x.valid)
   n.y <- sum(y.valid)

   result <- matrix(0, nrow = n.x, ncol = n.y)

   # Calculate overlap coefficients
   for (i in 1:n.x) {
       for (j in 1:n.y) {
           intersection.size <- length(intersect(x.basins[[i]], y.basins[[j]]))
           min.size <- min(length(x.basins[[i]]), length(y.basins[[j]]))
           result[i, j] <- intersection.size / min.size
       }
   }

   # Set row and column names using labels if provided, otherwise use lmax indices
   valid.x.lmax <- x.lmax[x.valid]
   valid.y.lmax <- y.lmax[y.valid]

   if (!is.null(x.labels)) {
       rownames(result) <- x.labels[valid.x.lmax]
   } else {
       rownames(result) <- valid.x.lmax
   }

   if (!is.null(y.labels)) {
       colnames(result) <- y.labels[valid.y.lmax]
   } else {
       colnames(result) <- valid.y.lmax
   }

   return(result)
}


#' Format Numbers in a Matrix with Specified Decimal Places
#'
#' @description
#' Formats all numeric elements in a matrix to have a specified number of decimal places,
#' ensuring a consistent display with leading and trailing zeros as needed.
#'
#' @param x Input matrix (numeric)
#' @param digits Number of digits to round to (default = 2)
#' @param nsmall Minimum number of digits to display after decimal point (default = 2)
#' @param ... Additional arguments (currently unused)
#'
#' @return A character matrix with all numbers formatted according to specifications
#'
#' @examples
#' m <- matrix(c(0.7, 1.235, 0.1, 0.8876), nrow = 2)
#' matrix.format(m)               # default: 2 digits
#' matrix.format(m, digits = 3)   # 3 digits
#' matrix.format(m, nsmall = 3)   # at least 3 decimal places
#'
#' @export
matrix.format <- function(x, digits = 2, nsmall = 2, ...) {
   # Input validation
   if (!is.matrix(x)) {
       stop("Input must be a matrix")
   }
   if (!is.numeric(digits) || digits < 0) {
       stop("digits must be a non-negative number")
   }
   if (!is.numeric(nsmall) || nsmall < 0) {
       stop("nsmall must be a non-negative number")
   }

   # Preserve row and column names
   rnames <- rownames(x)
   cnames <- colnames(x)

   # Format the matrix
   result <- apply(x, 1:2, function(val) format(round(val, digits), nsmall = nsmall))

   # Restore row and column names
   rownames(result) <- rnames
   colnames(result) <- cnames

   return(result)
}

#' Perform Robust Z-Score Normalization
#'
#' @description
#' Normalizes data using robust statistics by subtracting the median and dividing by the Median
#' Absolute Deviation (MAD). This method is more resistant to outliers compared to traditional
#' z-score normalization that uses mean and standard deviation.
#'
#' @param x A numeric vector to be normalized
#'
#' @return A numeric vector of the same length as the input containing robust z-scores.
#'         Missing values (NA) in the input will result in NA in the output at the same positions.
#'
#' @details
#' The function implements robust z-score normalization using the formula:
#' z = (x - median(x)) / MAD(x)
#'
#' The MAD is scaled by a factor of 1.4826 to make it consistent with the standard deviation
#' when the data is normally distributed.
#'
#' @note
#' If the MAD is zero, the function will return a vector of zeros and issue a warning.
#'
#' @examples
#' # Basic usage
#' data <- c(1, 2, 3, 100, 4, 5, 6)  # Note the outlier
#' robust.z.normalize(data)
#'
#' # Handling missing values
#' data.with.na <- c(1, NA, 3, 100, 4, NA, 6)
#' robust.z.normalize(data.with.na)
#'
#' @export
robust.z.normalize <- function(x) {
  # Check if input is numeric
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }

  # Remove NA values for calculations
  x.clean <- na.omit(x)

  # Calculate median
  x.median <- median(x.clean)

  # Calculate MAD
  # Using constant c=1.4826 for consistency with normal distribution
  mad.value <- mad(x.clean, constant = 1.4826)

  # Check if MAD is zero to avoid division by zero
  if (mad.value == 0) {
    warning("MAD is zero, returning zeros to avoid division by zero")
    return(rep(0, length(x)))
  }

  # Calculate robust z-scores
  z.scores <- (x - x.median) / mad.value

  return(z.scores)
}

#' Compute Standardized Effect Sizes with Bayesian Inference
#'
#' @description
#' Computes standardized effect sizes and associated Bayesian statistics for comparing
#' multiple methods against a reference method. The standardization is performed using
#' the standard deviation of the reference method. The function provides point estimates
#' (mean and median effects), interval estimates (credible intervals and highest density
#' intervals), and posterior probabilities for different hypotheses.
#'
#' @param bb.integrals List of Bayesian bootstrap integral values, where each element
#'        contains bootstrap samples for a method. The list must be named with method
#'        identifiers.
#' @param reference_idx Index of the reference method in the list (default: 1)
#'
#' @details
#' The function computes several Bayesian statistics:
#' - Effect sizes are standardized by dividing differences by the reference method's SD
#' - 95% credible intervals show the central 95% of the effect size distribution
#' - 95% HDI (Highest Density Interval) shows the most probable 95% of effect sizes
#' - Posterior probabilities are computed for:
#'   * Method having smaller effect than reference (prob_smaller)
#'   * Method having larger effect than reference (prob_larger)
#'   * Effect being practically equivalent to reference (within +/-0.1 SD)
#'
#' @return A data frame where each row represents a method (named by method identifiers)
#'         and columns contain:
#' \itemize{
#'   \item mean_effect: Mean standardized effect size
#'   \item median_effect: Median standardized effect size
#'   \item ci_lower: Lower bound of 95% credible interval
#'   \item ci_upper: Upper bound of 95% credible interval
#'   \item hdi_lower: Lower bound of 95% highest density interval
#'   \item hdi_upper: Upper bound of 95% highest density interval
#'   \item prob_smaller: Probability of smaller effect than reference
#'   \item prob_larger: Probability of larger effect than reference
#'   \item prob_practical_equiv: Probability of practical equivalence
#' }
#'
#' @examples
#' # Create example data
#' bb.integrals <- list(
#'   method1 = rnorm(100, mean = 0, sd = 1),
#'   method2 = rnorm(100, mean = 0.5, sd = 1),
#'   method3 = rnorm(100, mean = -0.3, sd = 1.2)
#' )
#' 
#' results <- compute.bayesian.effects(bb.integrals)
#'
#' # View results for specific method
#' results["method2", ]
#'
#' # Compare probabilities
#' results[, c("prob_smaller", "prob_larger", "prob_practical_equiv")]
#'
#' @export
compute.bayesian.effects <- function(bb.integrals, reference_idx = 1) {
   ref_samples <- unlist(bb.integrals[[reference_idx]])
   ref_sd <- sd(ref_samples)
   n_methods <- length(bb.integrals)

   # Initialize results matrix
   result_df <- data.frame(
       mean_effect = numeric(n_methods),
       median_effect = numeric(n_methods),
       ci_lower = numeric(n_methods),
       ci_upper = numeric(n_methods),
       hdi_lower = numeric(n_methods),
       hdi_upper = numeric(n_methods),
       prob_smaller = numeric(n_methods),
       prob_larger = numeric(n_methods),
       prob_practical_equiv = numeric(n_methods)
   )

   # Set row names to method names
   rownames(result_df) <- names(bb.integrals)

   # Compute statistics for each method
   for(i in seq_len(n_methods)) {
       current_samples <- unlist(bb.integrals[[i]])
       effect_samples <- (current_samples - ref_samples)/ref_sd

       ci <- quantile(effect_samples, c(0.025, 0.975))
       # Use HDInterval if available, otherwise use quantile
       if (requireNamespace("HDInterval", quietly = TRUE)) {
         hdi <- HDInterval::hdi(effect_samples, 0.95)
       } else {
         hdi <- quantile(effect_samples, c(0.025, 0.975))  # fallback to CI
       }

       result_df[i, ] <- c(
           mean(effect_samples),
           median(effect_samples),
           ci[1],
           ci[2],
           hdi[1],
           hdi[2],
           mean(effect_samples < 0),
           mean(effect_samples > 0),
           mean(abs(effect_samples) < 0.1)
       )
   }

   return(result_df)
}

#' Compute Pairwise Bayes Factors for Method Comparisons
#'
#' @description
#' Computes pairwise Bayes factors comparing all methods using Savage-Dickey density ratios.
#' Each Bayes factor represents evidence for one method versus another.
#'
#' @param bb.integrals A named list of Bayesian bootstrap integral values, where each element
#'   contains numeric samples representing the performance metric for a method. Names of the
#'   list elements will be used as row and column names in the output matrix.
#' @return A square matrix of Bayes factors where \code{BF[i,j]} represents the evidence for
#'   method i being better than method j. Values > 1 indicate evidence for method i, values < 1
#'   indicate evidence for method j. The diagonal contains 1s (method compared to itself).
#' @details The function uses the Savage-Dickey density ratio to compute Bayes factors:
#'   \itemize{
#'     \item For each pair of methods, it computes the difference in their performance samples
#'     \item Estimates the posterior density of these differences using kernel density estimation
#'     \item Computes the prior density at 0 assuming a normal distribution
#'     \item The Bayes factor is the ratio: posterior density at 0 / prior density at 0
#'     \item When density estimation fails, uses extreme values (1000 or 0.001) based on mean difference
#'   }
#' @note
#'   \itemize{
#'     \item Uses Sheather-Jones bandwidth selection for robust density estimation
#'     \item Interpretation: BF > 10 (strong evidence), BF > 3 (moderate evidence), BF > 1 (weak evidence)
#'     \item The function assumes that higher values in bb.integrals indicate better performance
#'   }
#' @examples
#' # Example: Compare three methods
#' set.seed(123)
#' bb.integrals <- list(
#'   method_A = rnorm(1000, mean = 0.5, sd = 0.1),
#'   method_B = rnorm(1000, mean = 0.52, sd = 0.1),
#'   method_C = rnorm(1000, mean = 0.48, sd = 0.1)
#' )
#'
#' bf_matrix <- compute.pairwise.bayes.factors(bb.integrals)
#' print(bf_matrix)
#'
#' # Interpret results
#' # bf_matrix[1,2] > 1 suggests evidence for method_A over method_B
#' # bf_matrix[2,3] > 1 suggests evidence for method_B over method_C
#' @references
#' Wagenmakers, E. J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010).
#' Bayesian hypothesis testing for psychologists: A tutorial on the Savage-Dickey method.
#' Cognitive Psychology, 60(3), 158-189.
#' @export
compute.pairwise.bayes.factors <- function(bb.integrals) {
    n_methods <- length(bb.integrals)
    bf_matrix <- matrix(NA, n_methods, n_methods)
    rownames(bf_matrix) <- colnames(bf_matrix) <- names(bb.integrals)

    for(i in 1:n_methods) {
        for(j in 1:n_methods) {
            if(i != j) {
                samples_i <- unlist(bb.integrals[[i]])
                samples_j <- unlist(bb.integrals[[j]])
                diff_samples <- samples_i - samples_j

                # Compute density with larger bandwidth for more stable estimation
                posterior_density <- try(density(diff_samples,
                                              bw = "SJ",    # Sheather-Jones bandwidth
                                              n = 1024,     # More points for better resolution
                                              from = min(diff_samples) - sd(diff_samples),
                                              to = max(diff_samples) + sd(diff_samples)),
                                      silent = TRUE)

                if(!inherits(posterior_density, "try-error")) {
                    # Compute prior density at 0
                    prior_density_0 <- dnorm(0, mean(diff_samples), sd(diff_samples))

                    # Get posterior density at 0
                    posterior_at_0 <- approx(posterior_density$x,
                                           posterior_density$y,
                                           xout = 0)$y

                    # If density estimation succeeded
                    if(!is.na(posterior_at_0)) {
                        bf_matrix[i,j] <- posterior_at_0/prior_density_0
                    } else {
                        # If methods are very different, use extreme value
                        mean_diff <- mean(diff_samples)
                        if(mean_diff > 0) {
                            bf_matrix[i,j] <- 1000  # Strong evidence for i > j
                        } else {
                            bf_matrix[i,j] <- 0.001 # Strong evidence for j > i
                        }
                    }
                }
            }
        }
    }

    # Set diagonal to 1 (comparing method with itself)
    diag(bf_matrix) <- 1

    return(bf_matrix)
}

#' Perform Hierarchical Bayesian Analysis of Method Differences
#'
#' @description
#' Fits a hierarchical Bayesian model to compare multiple methods, estimating method-specific
#' means and variances while accounting for between-method variability. The model provides
#' posterior distributions for all pairwise method differences and probabilities of superiority.
#'
#' @param bb.integrals A named list of Bayesian bootstrap integral values, where each element
#'   contains numeric samples representing the performance metric for a method. All elements
#'   must have the same length.
#' @return A stanfit object containing:
#'   \itemize{
#'     \item \code{mu}: Posterior samples of method-specific means
#'     \item \code{sigma}: Posterior samples of method-specific standard deviations
#'     \item \code{tau}: Posterior samples of between-method variability
#'     \item \code{global_mu}: Posterior samples of the overall mean across methods
#'     \item \code{diff}: Matrix of pairwise differences (mu\[i\] - mu\[j\])
#'     \item \code{prob_diff}: Matrix of probabilities that method i is worse than method j
#'   }
#' @details The hierarchical model structure:
#'   \itemize{
#'     \item Method means: \code{mu[k] ~ Normal(global_mu, tau)}
#'     \item Observations: \code{y[n,k] ~ Normal(mu[k], sigma[k])}
#'     \item Priors: \code{global_mu ~ Normal(0, 10)}, \code{tau ~ Cauchy(0, 2.5)},
#'           \code{sigma ~ Cauchy(0, 2.5)}
#'     \item The model pools information across methods through the hierarchical structure
#'   }
#' @note
#'   \itemize{
#'     \item Requires the 'rstan' package to be installed
#'     \item Uses 4 chains with 2000 iterations (1000 warmup, 1000 sampling)
#'     \item \code{prob_diff\[i,j\]} represents P(mu\[i\] < mu\[j\]), so values < 0.5 indicate method i is better
#'     \item Check convergence using \code{rstan::check_hmc_diagnostics()} and \code{summary()}
#'   }
#' @examples
#' \dontrun{
#' # Example: Compare three methods with hierarchical model
#' set.seed(123)
#' bb.integrals <- list(
#'   method_A = rnorm(100, mean = 0.5, sd = 0.1),
#'   method_B = rnorm(100, mean = 0.52, sd = 0.12),
#'   method_C = rnorm(100, mean = 0.48, sd = 0.08)
#' )
#'
#' # Fit hierarchical model
#' fit <- analyze.hierarchical.differences(bb.integrals)
#'
#' # Check convergence
#' rstan::check_hmc_diagnostics(fit)
#'
#' # Extract posterior summaries
#' print(fit, pars = c("mu", "sigma", "tau"))
#'
#' # Get pairwise probabilities
#' prob_matrix <- rstan::extract(fit, "prob_diff")$prob_diff
#' apply(prob_matrix, c(2,3), mean)  # Mean probabilities
#' }
#' @seealso
#'   \code{\link[rstan]{stan}} for model fitting,
#'   \code{\link[rstan]{extract}} for extracting posterior samples,
#'   \code{\link[rstan]{check_hmc_diagnostics}} for convergence diagnostics
#' @references
#' Gelman, A., & Hill, J. (2007). Data analysis using regression and
#' multilevel/hierarchical models. Cambridge University Press.
#' @export
analyze.hierarchical.differences <- function(bb.integrals) {
    # Prepare data for Stan
    n_methods <- length(bb.integrals)
    n_samples <- length(bb.integrals[[1]])

    # Convert to matrix
    y <- sapply(bb.integrals, unlist)

    # Stan model code for hierarchical analysis
    stan_code <- "
    data {
        int<lower=0> N;
        int<lower=0> K;
        matrix[N, K] y;
    }
    parameters {
        vector[K] mu;
        vector<lower=0>[K] sigma;
        real<lower=0> tau;
        real global_mu;
    }
    model {
        // Priors
        global_mu ~ normal(0, 10);
        tau ~ cauchy(0, 2.5);

        mu ~ normal(global_mu, tau);
        sigma ~ cauchy(0, 2.5);

        // Likelihood
        for(k in 1:K)
            y[,k] ~ normal(mu[k], sigma[k]);
    }
    generated quantities {
        matrix[K,K] diff;
        matrix[K,K] prob_diff;
        for(i in 1:K) {
            for(j in 1:K) {
                diff[i,j] = mu[i] - mu[j];
                prob_diff[i,j] = normal_cdf(0, mu[i] - mu[j],
                    sqrt(sigma[i]^2 + sigma[j]^2));
            }
        }
    }
    "

    # Fit model using rstan
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("Package 'rstan' is required for this function. Please install it with install.packages('rstan')")
    }
    fit <- rstan::stan(model_code = stan_code,
                data = list(N = n_samples,
                           K = n_methods,
                           y = y),
                chains = 4,
                iter = 2000)

    return(fit)
}

#' Format Matrix of Bayes Factors into Simplified Notation
#'
#' @description
#' Converts a matrix of Bayes factors into a more readable format using simplified
#' notation for extreme values and comparison symbols. This is particularly useful
#' for creating publication-ready tables with LaTeX math symbols.
#'
#' @param x A numeric matrix or data frame of Bayes factors, typically from pairwise
#'   method comparisons where \code{x[i,j]} represents evidence for method i vs method j
#' @param thresholds Named numeric vector of thresholds for different symbols
#'        (default: \code{c(strong = 1000, moderate = 10, weak = 3)}). Values are used
#'        for both directions (e.g., BF > 10 and BF < 1/10)
#' @param digits Number of decimal places for non-extreme values (default: 2)
#' @param ... Additional arguments (currently unused)
#' @return A character matrix with the same dimensions and names as the input, where numeric
#'         values are replaced with symbolic notation:
#'         \itemize{
#'           \item \code{"$\\gg$"}: BF \eqn{\ge} strong threshold (very strong evidence)
#'           \item \code{"$>$"}: weak \eqn{\le} BF < moderate (weak to moderate evidence)
#'           \item \code{"$<$"}: 1/moderate < BF \eqn{\le} 1/weak (weak to moderate evidence against)
#'           \item \code{"$\\ll$"}: BF \eqn{\le} 1/strong (very strong evidence against)
#'           \item Numeric value: 1/weak < BF < weak (inconclusive)
#'           \item \code{"1"}: Diagonal elements (method compared to itself)
#'         }
#' @details The function maps Bayes factors to LaTeX symbols based on evidence strength:
#'   \itemize{
#'     \item Strong evidence (BF \eqn{\ge} 1000 or BF \eqn{\le} 0.001) uses double symbols
#'     \item Moderate evidence (BF \eqn{\ge} 10 or BF \eqn{\le} 0.1) also uses double symbols
#'     \item Weak evidence (BF \eqn{\ge} 3 or BF \eqn{\le} 0.333) uses single symbols
#'     \item Values between thresholds are displayed numerically
#'   }
#' @note
#'   \itemize{
#'     \item The output contains LaTeX math symbols suitable for knitr/RMarkdown documents
#'     \item The current implementation uses the same symbol ($\\gg$) for both strong
#'           and moderate evidence, which may be confusing
#'     \item For interpretation: symbols point toward the better method
#'   }
#' @examples
#' # Create example Bayes factor matrix
#' bf.matrix <- matrix(c(1, 15.2, 0.05, 1200,
#'                       0.066, 1, 0.002, 8.5,
#'                       20.1, 500, 1, 2.1,
#'                       0.0008, 0.118, 0.476, 1),
#'                     nrow = 4, byrow = TRUE)
#' rownames(bf.matrix) <- colnames(bf.matrix) <- c("A", "B", "C", "D")
#'
#' # Format with default thresholds
#' bayes.factors.symbols(bf.matrix)
#'
#' # Format with custom thresholds
#' bayes.factors.symbols(bf.matrix,
#'                      thresholds = c(strong = 100, moderate = 10, weak = 3),
#'                      digits = 1)
#'
#' # Use in knitr/RMarkdown for LaTeX output
#' # knitr::kable(bayes.factors.symbols(bf.matrix), escape = FALSE)
#' @seealso
#'   \code{\link{compute.pairwise.bayes.factors}} for generating Bayes factor matrices
#' @export
bayes.factors.symbols <- function(x,
                                  thresholds = c(strong = 1000, moderate = 10, weak = 3),
                                  digits = 2, ...) {
   bf.matrix <- x  # Use standard S3 first parameter name

   # Convert to matrix if data frame
   if(is.data.frame(bf.matrix)) {
       bf.matrix <- as.matrix(bf.matrix)
   }

   # Validate input
   if(!is.numeric(bf.matrix)) {
       stop("Input must be numeric")
   }

   # Create output matrix
   result <- matrix("", nrow = nrow(bf.matrix), ncol = ncol(bf.matrix))
   dimnames(result) <- dimnames(bf.matrix)

   # Ensure thresholds have names and sort in descending order
   if(is.null(names(thresholds))) {
       names(thresholds) <- c("strong", "moderate", "weak")[1:length(thresholds)]
   }
   thresholds <- sort(thresholds, decreasing = TRUE)

   # Extract threshold values with defaults
   strong_thresh <- if("strong" %in% names(thresholds)) thresholds["strong"] else 1000
   moderate_thresh <- if("moderate" %in% names(thresholds)) thresholds["moderate"] else 10
   weak_thresh <- if("weak" %in% names(thresholds)) thresholds["weak"] else 3

   # Function to format a single value
   format.value <- function(x) {
       if(is.na(x)) return("NA")

       # Handle diagonal
       if(abs(x - 1) < .Machine$double.eps) return("1")

       if(x >= strong_thresh) {
           return("$\\gg$")
       } else if(x >= moderate_thresh) {
           return("$\\gg$")
       } else if(x >= weak_thresh) {
           return("$>$")
       } else if(x <= 1/strong_thresh) {
           return("$\\ll$")
       } else if(x <= 1/moderate_thresh) {
           return("$\\ll$")
       } else if(x <= 1/weak_thresh) {
           return("$<$")
       } else {
           return(sprintf(paste0("%.", digits, "f"), x))
       }
   }

   # Apply formatting to each element
   for(i in 1:nrow(bf.matrix)) {
       for(j in 1:ncol(bf.matrix)) {
           result[i,j] <- format.value(bf.matrix[i,j])
       }
   }

   return(result)
}

#' Finds inflection points of a vector
#'
#' @description
#' Identifies inflection points in a numeric vector by computing the second derivative
#' using either the MAGELO (Monotonic Averaging with Gaussian Error and Linear Order) or
#' MABILO (Monotonic Averaging with Binomial Error and Linear Order) smoothing methods.
#' Inflection points are where the second derivative equals zero.
#'
#' @param y A numeric vector of values for which to find inflection points.
#' @param x Optional numeric vector of x-coordinates corresponding to y values.
#'   If NULL (default), uses indices 1:length(y). Must be the same length as y.
#' @param method Character string specifying the smoothing method to use. Either "magelo"
#'   (Gaussian error model) or \code{\link{mabilo}} (Binomial error model). Default is "magelo".
#'   Only the first element is used if a vector is provided.
#' @return A list containing:
#'   \itemize{
#'     \item \code{y}: The sorted y values (reordered if x was provided)
#'     \item \code{r0}: Result from the initial smoothing of y
#'     \item \code{r1}: Result from smoothing the first derivative
#'     \item \code{r2}: Result from smoothing the second derivative
#'     \item \code{infl.pts}: Integer indices of inflection points in the sorted data
#'     \item \code{infl.pts.vals}: The y values at the inflection points
#'   }
#' @details The function uses a three-step process:
#'   \enumerate{
#'     \item Smooths the original data using the specified method
#'     \item Computes and smooths the first derivative (differences of smoothed values)
#'     \item Computes and smooths the second derivative
#'     \item Finds zeros of the second derivative using root-finding
#'   }
#'   If x values are provided, the data is first sorted by x before analysis.
#' @note
#'   \itemize{
#'     \item Requires the 'rootSolve' package for finding zeros of the second derivative
#'     \item Also requires either the 'magelo' or \code{\link{mabilo}} function to be available
#'     \item The inflection point indices are floored to integers, which may introduce slight imprecision
#'     \item When x is NULL, the function sorts y values, which changes the interpretation of results
#'   }
#' @examples
#' \dontrun{
#' # Example 1: Find inflection points in a sigmoid-like curve
#' x <- seq(-5, 5, length.out = 100)
#' y <- 1 / (1 + exp(-x)) + 0.1 * sin(2*x)  # Sigmoid with oscillation
#'
#' result <- find.inflection.pts(y, x, method = "magelo")
#'
#' # Plot the results
#' plot(x, y, type = "l", main = "Inflection Points")
#' points(x[result$infl.pts], result$infl.pts.vals, col = "red", pch = 19)
#'
#' # Example 2: Using indices as x-coordinates
#' y <- c(1, 2, 5, 10, 15, 18, 19, 19.5, 19.8, 20)  # Growth curve
#' result <- find.inflection.pts(y, method = "mabilo")
#' print(result$infl.pts)
#' }
#' @seealso
#'   \code{\link[rootSolve]{uniroot.all}} for root finding,
#'   \code{magelo} and \code{mabilo} for the smoothing methods
#' @export
find.inflection.pts <- function(y, x = NULL, method = c("magelo","mabilo")) {

    if (!is.null(x)) {
        if (length(x) != length(y)) stop("x and y have to be of the same length")
        o <- order(x)
        x <- x[o]
        y <- y[o]
    } else {
        x <- seq(y)
        y <- sort(y)
    }

    r0 <- r1 <- r2 <- NULL
    if (method == "magelo") {
        r0 <- magelo(x, y)
        r1 <- magelo(r0$xgrid[-length(r0$xgrid)], diff(r0$gpredictions))
        r2 <- magelo(r1$xgrid[-length(r1$xgrid)], diff(r1$gpredictions))
        r2fn <- approxfun(r2$xg, r2$Eyg)
    } else if (method == "mabilo") {
        n <- length(x)
        r0 <- mabilo(x, y)

        x1 <- x[-n]
        r1 <- mabilo(x1, diff(r0$predictions))

        x2 <- x1[-length(x1)]
        r2 <- mabilo(x2, diff(r1$predictions))
        r2fn <- approxfun(x2, r2$predictions)
    }

    if (!requireNamespace("rootSolve", quietly = TRUE)) {
      stop("Package 'rootSolve' is required for this function. Please install it with install.packages('rootSolve')")
    }
    infl.pts <- rootSolve::uniroot.all(r2fn, interval=range(r2$xg))
    infl.pts <- floor(infl.pts)
    infl.pts.vals <- y[infl.pts]

    list(y=y,
         r0=r0,
         r1=r1,
         r2=r2,
         infl.pts=infl.pts,
         infl.pts.vals=infl.pts.vals)
}

#' Calculate Second Derivative of Gaussian Function
#'
#' @description
#' This function computes the second derivative of the unnormalized Gaussian function:
#' \deqn{f(x) = \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)}{f(x) = exp(-((x - mu)^2)/(2 * sd^2))}
#' with respect to x. The second derivative identifies inflection points of the Gaussian
#' curve, which occur at \eqn{\mu \pm \sigma}{mu +/- sd}.
#'
#' @param x Numeric vector of x values at which to evaluate the second derivative
#' @param mu Mean parameter (location) of the Gaussian function
#' @param sd Standard deviation parameter (scale) of the Gaussian function. Must be positive.
#' @return Numeric vector of second derivative values with the same length as x. The second
#'   derivative equals zero at the inflection points (x = mu +/- sd), is negative between
#'   them (concave down), and positive outside them (concave up).
#' @details The second derivative is computed using the formula:
#'   \deqn{f''(x) = \frac{(z^2 - 1) \cdot \exp(-z^2/2)}{\sigma^2}}{f''(x) = (z^2 - 1) * exp(-z^2/2) / sd^2}
#'   where \eqn{z = (x - \mu)/\sigma}{z = (x - mu)/sd} is the standardized value.
#'
#'   Key properties:
#'   \itemize{
#'     \item Zero at x = mu +/- sd (inflection points)
#'     \item Maximum at x = mu (where f''(mu) = -1/sd^2)
#'     \item Approaches 0 as \eqn{x \to \pm \infty}
#'   }
#' @note This computes the second derivative of the unnormalized Gaussian (without the
#'   \eqn{1/(\sigma\sqrt{2\pi})}{1/(sd*sqrt(2*pi))} normalization factor). For the
#'   normalized version, multiply the result by \eqn{1/(\sigma\sqrt{2\pi})}{1/(sd*sqrt(2*pi))}.
#' @examples
#' # Example 1: Basic usage
#' x <- seq(-5, 5, length.out = 100)
#' result <- gaussian.second.derivative(x, mu = 0, sd = 1)
#'
#' # Example 2: Verify inflection points at mu +/- sd
#' mu <- 2
#' sd <- 1.5
#' x_inflection <- c(mu - sd, mu + sd)
#' gaussian.second.derivative(x_inflection, mu, sd)  # Should be ~0
#'
#' # Example 3: Visualize the Gaussian and its second derivative
#' x <- seq(-4, 4, length.out = 200)
#' f <- exp(-(x^2)/2)  # Original function (mu=0, sd=1)
#' f_double_prime <- gaussian.second.derivative(x, mu = 0, sd = 1)
#'
#' par(mfrow = c(2, 1))
#' plot(x, f, type = "l", main = "Gaussian Function", ylab = "f(x)")
#' abline(v = c(-1, 1), col = "red", lty = 2)  # Inflection points
#' plot(x, f_double_prime, type = "l", main = "Second Derivative", ylab = "f''(x)")
#' abline(h = 0, col = "gray")
#' abline(v = c(-1, 1), col = "red", lty = 2)  # Zero crossings
#' @seealso
#'   \code{\link{dnorm}} for the normalized Gaussian density function,
#'   \code{\link{find.inflection.pts}} for finding inflection points in general curves
#' @export
gaussian.second.derivative <- function(x, mu, sd) {
    # Standardized difference
    z <- (x - mu)/sd

    # Base exponential term
    exp.term <- exp(-z^2/2)

    # Calculate second derivative
    # Using the chain rule and product rule twice
    # The result is: (z^2 - 1) * exp(-z^2/2) / sd^2
    second.deriv <- (z^2 - 1) * exp.term / sd^2

    return(second.deriv)
}

#' Calculate Jaccard index between two sets represented as numeric vectors
#'
#' @description
#' The Jaccard index (also known as Jaccard similarity coefficient) is defined as the ratio
#' of the size of the intersection of the two sets to the size of their union. Formally:
#' \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|}}{J(A, B) = |A  intersect  B| / |A  union  B|}
#' The index ranges from 0 (no overlap) to 1 (identical sets).
#'
#' @param vec.1 A numeric vector representing the first set. Duplicate values are
#'   automatically removed to form a proper set.
#' @param vec.2 A numeric vector representing the second set. Duplicate values are
#'   automatically removed to form a proper set.
#' @return A numeric value between 0 and 1 representing the Jaccard index, where:
#'   \itemize{
#'     \item 0 indicates no overlap between the sets
#'     \item 1 indicates the sets are identical
#'     \item Values between 0 and 1 indicate partial overlap
#'   }
#'   Returns NaN if both sets are empty (0/0 case).
#' @details The function first converts the input vectors to sets by removing duplicates
#'   using \code{unique()}, then computes the intersection and union using R's built-in
#'   set operations. This ensures that the calculation follows proper set theory principles
#'   where each element is counted only once.
#' @note
#'   \itemize{
#'     \item The function treats the inputs as sets, so duplicate values are ignored
#'     \item For comparing multisets (where duplicates matter), consider other similarity measures
#'     \item The Jaccard distance can be computed as 1 - Jaccard index
#'     \item NA values in the input vectors are preserved and treated as distinct elements
#'   }
#' @examples
#' # Example 1: Basic usage with duplicates
#' samples.1 <- c(1, 2, 5, 5)
#' samples.2 <- c(2, 2, 5, 7)
#' jaccard.index(samples.1, samples.2)  # Returns 0.5 (intersection: {2,5}, union: {1,2,5,7})
#'
#' # Example 2: Identical sets
#' jaccard.index(c(1, 2, 3), c(3, 2, 1))  # Returns 1
#'
#' # Example 3: No overlap
#' jaccard.index(c(1, 2, 3), c(4, 5, 6))  # Returns 0
#'
#' # Example 4: One empty set
#' jaccard.index(c(1, 2, 3), numeric(0))  # Returns 0
#'
#' # Example 5: Comparing sample indices or cluster memberships
#' cluster1_samples <- c(1, 3, 5, 7, 9)
#' cluster2_samples <- c(2, 4, 5, 7, 8)
#' similarity <- jaccard.index(cluster1_samples, cluster2_samples)
#' print(paste("Cluster overlap:", round(similarity * 100, 1), "%"))
#' @references
#' Jaccard, P. (1912). The distribution of the flora in the alpine zone.
#' New Phytologist, 11(2), 37-50.
#' @seealso
#'   \code{\link{intersect}} for set intersection,
#'   \code{\link{union}} for set union,
#'   \code{\link{setdiff}} for set difference
#' @export
jaccard.index <- function(vec.1, vec.2) {
  ## Convert input vectors to sets by removing duplicates
  set.1 <- unique(vec.1)
  set.2 <- unique(vec.2)

  ## Compute the intersection and union
  intersection.size <- length(intersect(set.1, set.2))
  union.size <- length(union(set.1, set.2))

  ## Return the Jaccard index
  intersection.size / union.size
}


#' Sigmoidal Function with Lower and Upper Thresholds
#'
#' Creates a sigmoidal function that smoothly transitions from approximately 0 to
#' approximately 1 between two specified threshold values. This is useful for
#' creating soft thresholds in data processing or for probability-based filtering.
#'
#' @param x Numeric vector or scalar. The input values to transform.
#' @param lower_threshold Numeric scalar. The lower threshold value where the
#'        function begins its transition from 0 to 1. For values well below this
#'        threshold, the function returns values close to 0.
#' @param upper_threshold Numeric scalar. The upper threshold value where the
#'        function completes its transition from 0 to 1. For values well above this
#'        threshold, the function returns values close to 1.
#' @param steepness Numeric scalar. Controls the steepness of the transition between
#'        thresholds. Higher values create a sharper transition. Default is 10.
#'
#' @return A numeric vector of the same length as `x` with values ranging from
#'         approximately 0 to approximately 1.
#'
#' @examples
#' # Define thresholds
#' q90 <- 0.1  # Lower threshold (e.g., 90th percentile)
#' q95 <- 0.15 # Upper threshold (e.g., 95th percentile)
#'
#' # Apply to a sequence of values
#' x_seq <- seq(0, 0.3, length.out = 100)
#' y_values <- thresholded.sigmoid(x_seq, q90, q95)
#'
#' # Plot the result
#' plot(x_seq, y_values, type = "l",
#'      xlab = "Input Value", ylab = "Transformed Value",
#'      main = "Thresholded Sigmoid Function")
#' abline(v = c(q90, q95), lty = 2, col = c("blue", "red"))
#'
#' @export
thresholded.sigmoid <- function(x, lower_threshold, upper_threshold, steepness = 10) {
  # Calculate midpoint between thresholds
  midpoint <- (lower_threshold + upper_threshold) / 2

  # Calculate scale factor to control steepness of transition
  # Higher steepness values create sharper transitions
  scale_factor <- steepness / (upper_threshold - lower_threshold)

  # Apply logistic function
  1 / (1 + exp(-scale_factor * (x - midpoint)))
}


#' Winsorized Z-score normalization
#'
#' @param data A numeric matrix or data frame where rows are samples and columns are features
#' @param limits A numeric vector of length 2 specifying the lower and upper proportion
#'        of values to be winsorized (default: c(0.01, 0.01) for 1% at each tail)
#' @return A matrix of winsorized and Z-scored values
#' @examples
#' # Example with random data
#' set.seed(123)
#' example.data <- matrix(rnorm(100, 5, 2), ncol=5)
#' # Add some outliers
#' example.data[1,1] <- 25
#' example.data[2,3] <- -15
#' normalized.data <- winsorize.zscore(example.data)
#' @export
winsorize.zscore <- function(data, limits = c(0.01, 0.01)) {
    ## Convert data to matrix if it's a data frame
    if(is.data.frame(data)) {
        data <- as.matrix(data)
    }

    ## Initialize the result matrix
    winsorized.data <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    colnames(winsorized.data) <- colnames(data)

    ## Process each column separately
    for(col in 1:ncol(data)) {
        ## Extract column
        x <- data[, col]

        ## Calculate quantiles for winsorization
        lower.quantile <- quantile(x, limits[1], na.rm = TRUE)
        upper.quantile <- quantile(x, 1 - limits[2], na.rm = TRUE)

        ## Winsorize the data
        winsorized.x <- pmin(pmax(x, lower.quantile), upper.quantile)

        ## Z-score the winsorized data
        mean.x <- mean(winsorized.x, na.rm = TRUE)
        sd.x <- sd(winsorized.x, na.rm = TRUE)

        ## Handle the case where sd is 0 (constant feature)
        if(sd.x == 0) sd.x <- 1e-10

        winsorized.data[, col] <- (winsorized.x - mean.x) / sd.x
    }

    return(winsorized.data)
}

#' Robust Z-score normalization using median and MAD
#'
#' @param data A numeric matrix or data frame where rows are samples and columns are features
#' @param scale.factor Scaling factor for MAD (default: 1.4826 to make MAD consistent with
#'        standard deviation for normal distributions)
#' @return A matrix of robust Z-scored values
#' @examples
#' # Example with random data
#' set.seed(123)
#' example.data <- matrix(rnorm(100, 5, 2), ncol=5)
#' # Add some outliers
#' example.data[1,1] <- 25
#' example.data[2,3] <- -15
#' normalized.data <- robust.zscore(example.data)
#' @export
robust.zscore <- function(data, scale.factor = 1.4826) {
    ## Convert data to matrix if it's a data frame
    if(is.data.frame(data)) {
        data <- as.matrix(data)
    }

    ## Initialize the result matrix
    robust.data <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    colnames(robust.data) <- colnames(data)

    ## Process each column separately
    for(col in 1:ncol(data)) {
        ## Extract column
        x <- data[, col]

        ## Calculate median
        med.x <- median(x, na.rm = TRUE)

        ## Calculate MAD
        mad.x <- median(abs(x - med.x), na.rm = TRUE) * scale.factor

        ## Handle the case where MAD is 0 (constant feature)
        if(mad.x == 0) mad.x <- 1e-10

        ## Calculate robust Z-scores
        robust.data[, col] <- (x - med.x) / mad.x
    }

    return(robust.data)
}

#' Quantize Continuous Variable (Internal Function)
#'
#' This is a placeholder for the quantize.cont.var function that is referenced
#' but not defined in the original code. This function should be implemented
#' based on your specific requirements.
#'
#' @param x Numeric vector to quantize.
#' @param method Quantization method ("uniform" or "quantile").
#' @param wins.p Winsorization parameter.
#' @param round Whether to round endpoints.
#' @param dig.lab Number of digits for labels.
#' @param start Start value for color range.
#' @param end End value for color range.
#' @param n.levels Number of quantization levels.
#'
#' @return List containing x.cat and x.col.tbl
#'
#' @keywords internal
#' @noRd
quantize.cont.var <- function(x, method = "uniform", wins.p = 0.02,
                             round = FALSE, dig.lab = 2, start = 1/6,
                             end = 0, n.levels = 10) {

    x_original <- x  # Keep original for later

    if (method == "uniform") {
        # Winsorize if requested
        if (wins.p > 0 && wins.p < 0.5) {
            q_low <- quantile(x, wins.p, na.rm = TRUE)
            q_high <- quantile(x, 1 - wins.p, na.rm = TRUE)

            n_low <- sum(x < q_low, na.rm = TRUE)
            n_high <- sum(x > q_high, na.rm = TRUE)

            x[x < q_low] <- q_low
            x[x > q_high] <- q_high
        }
        # Create uniform breaks
        breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                     length.out = n.levels + 1)
    } else if (method == "quantile") {
        # Create quantile breaks
        breaks <- quantile(x, probs = seq(0, 1, length.out = n.levels + 1),
                          na.rm = TRUE)
        breaks <- unique(breaks)
    }

    if (round) {
        breaks <- round(breaks, dig.lab)
    }

    # CRITICAL FIX: Extend breaks slightly to ensure all values are included
    breaks[1] <- breaks[1] - .Machine$double.eps * abs(breaks[1]) - 1e-10
    breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps * abs(breaks[length(breaks)]) + 1e-10

    # Cut the variable - use the winsorized values
    x.cat <- cut(x, breaks = breaks, include.lowest = TRUE, dig.lab = dig.lab)

    # Create color table
    cols <- grDevices::rainbow(length(levels(x.cat)), start = start, end = end)
    x.col.tbl <- cols
    names(x.col.tbl) <- levels(x.cat)

    return(list(x.cat = x.cat, x.col.tbl = x.col.tbl))
}

#' Find Local Minima in Distance Sequence
#'
#' Identifies local minima in a sequence of edit distances, which can indicate
#' stable graph configurations.
#'
#' @param distances Numeric vector of distances.
#' @param k.values Optional numeric vector of corresponding k values. If NULL,
#'   indices are used.
#' @param threshold Numeric value for the minimum prominence of a local minimum
#'   to be considered significant (default: 0, all local minima).
#' @param na.rm Logical. Should NA values be handled? If TRUE, NA values are
#'   interpolated; if FALSE (default), function fails if NA values are present.
#' @param sort.by.prominence Logical. Should results be sorted by prominence?
#'   If TRUE (default), results are sorted with most prominent first.
#'   If FALSE, results are returned in order of occurrence.
#'
#' @return A data frame with columns:
#'   \item{index}{Integer indices of local minima}
#'   \item{k}{Corresponding k values (or indices if k.values is NULL)}
#'   \item{distance}{Distance values at local minima}
#'   \item{prominence}{Prominence of each local minimum}
#'
#' @details
#' A local minimum is defined as a point that is smaller than both its
#' neighbors. Prominence is calculated as the minimum height the function
#' must climb to reach a higher value.
#'
#' If \code{na.rm = TRUE}, NA values are linearly interpolated before
#' finding minima. Minima at originally NA positions are excluded from results.
#'
#' @examples
#' # Create sample data with clear minima
#' k <- 1:20
#' dist <- sin(k / 3) + 0.1 * rnorm(20)
#'
#' # Find local minima
#' minima <- find.local.minima(dist, k)
#' print(minima)
#'
#' # Plot with minima marked
#' plot(k, dist, type = "b")
#' points(minima$k, minima$distance, col = "red", pch = 19, cex = 1.5)
#'
#' # Example with NA handling
#' dist_na <- dist
#' dist_na[c(5, 15)] <- NA
#' minima_na <- find.local.minima(dist_na, k, na.rm = TRUE)
#'
#' @export
find.local.minima <- function(distances,
                              k.values = NULL,
                              threshold = 0,
                              na.rm = FALSE,
                              sort.by.prominence = TRUE) {

  # Input validation
  if (!is.numeric(distances)) {
    stop("distances must be a numeric vector")
  }

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a non-negative scalar")
  }

  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("na.rm must be a single logical value")
  }

  if (!is.logical(sort.by.prominence) || length(sort.by.prominence) != 1) {
    stop("sort.by.prominence must be a single logical value")
  }

  n <- length(distances)
  if (n < 3) {
    return(data.frame(
      index = integer(0),
      k = numeric(0),
      distance = numeric(0),
      prominence = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(k.values)) {
    k.values <- seq_along(distances)
  } else {
    if (!is.numeric(k.values)) {
      stop("k.values must be numeric")
    }
    if (length(k.values) != n) {
      stop("k.values must have the same length as distances")
    }
  }

  # Handle NA values
  original.na <- is.na(distances)
  if (any(original.na)) {
    if (!na.rm) {
      stop("distances contains NA values. Set na.rm = TRUE to handle them.")
    }

    # Linear interpolation for NA values
    distances <- approx(x = which(!original.na),
                       y = distances[!original.na],
                       xout = 1:n,
                       method = "linear",
                       rule = 2)$y  # rule = 2 uses nearest value for extrapolation
  }

  # Vectorized approach to find local minima
  # Create shifted versions for comparison
  left <- c(Inf, distances[-n])
  right <- c(distances[-1], Inf)

  # Find where current value is less than both neighbors
  is.local.min <- distances < left & distances < right

  # Don't consider originally NA positions as minima
  if (na.rm && any(original.na)) {
    is.local.min[original.na] <- FALSE
  }

  # Calculate prominence for each local minimum
  prominence <- rep(NA_real_, n)
  min.indices <- which(is.local.min)

  if (length(min.indices) > 0) {
    # Vectorized prominence calculation
    for (idx in min.indices) {
      # Find barriers on both sides
      left.barrier <- find.barrier(distances, idx, direction = "left")
      right.barrier <- find.barrier(distances, idx, direction = "right")

      prominence[idx] <- min(left.barrier - distances[idx],
                            right.barrier - distances[idx])
    }
  }

  # Filter by threshold
  valid <- is.local.min & !is.na(prominence) & prominence >= threshold

  result <- data.frame(
    index = which(valid),
    k = k.values[valid],
    distance = distances[valid],
    prominence = prominence[valid],
    stringsAsFactors = FALSE
  )

  # Sort by prominence if requested
  if (sort.by.prominence && nrow(result) > 0) {
    result <- result[order(result$prominence, decreasing = TRUE), ]
  }

  rownames(result) <- NULL
  return(result)
}

#' Find Barrier Height for Prominence Calculation
#'
#' Internal function to find the barrier height in a given direction
#' from a local minimum.
#'
#' @param distances Numeric vector of distances
#' @param idx Index of the local minimum
#' @param direction Either "left" or "right"
#'
#' @return The barrier height (maximum value before reaching a higher value
#'   than the minimum, or the maximum in that direction if no higher value exists)
#'
#' @keywords internal
#' @noRd
find.barrier <- function(distances, idx, direction = c("left", "right")) {
  direction <- match.arg(direction)
  n <- length(distances)
  min.value <- distances[idx]

  if (direction == "left") {
    if (idx == 1) return(min.value)

    # Search leftward
    search.range <- (idx - 1):1
    higher.idx <- which(distances[search.range] > min.value)

    if (length(higher.idx) == 0) {
      # No higher value found, return max of entire left side
      return(max(distances[1:(idx - 1)]))
    } else {
      # Found higher value, return max between minimum and that point
      boundary <- search.range[min(higher.idx)]
      return(max(distances[boundary:idx]))
    }
  } else {  # direction == "right"
    if (idx == n) return(min.value)

    # Search rightward
    search.range <- (idx + 1):n
    higher.idx <- which(distances[search.range] > min.value)

    if (length(higher.idx) == 0) {
      # No higher value found, return max of entire right side
      return(max(distances[(idx + 1):n]))
    } else {
      # Found higher value, return max between minimum and that point
      boundary <- search.range[min(higher.idx)]
      return(max(distances[idx:boundary]))
    }
  }
}

#' Finds local minima in a vector (internal/simplified version)
#'
#' Internal helper function that identifies local minima in a numeric vector
#' and returns the corresponding k values at those positions.
#'
#' @param x A numeric vector in which to find local minima
#' @param k.values A numeric vector of k values corresponding to positions in x
#'
#' @return A numeric vector containing the k values at positions where local minima occur.
#'   Returns an empty numeric vector if x has fewer than 3 elements or no local minima exist.
#'
#' @details
#' A local minimum is defined as a point where the value is less than both its
#' immediate neighbors. The function uses second differences to identify these points.
#'
#' @keywords internal
#' @noRd
internal.find.local.minima <- function(x, k.values) {
  # Input validation
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }
  if (!is.numeric(k.values)) {
    stop("k.values must be a numeric vector")
  }
  if (length(x) != length(k.values)) {
    stop("x and k.values must have the same length")
  }

  # Need at least 3 points to have a local minimum
  if (length(x) < 3) {
    return(numeric(0))
  }

  # Find local minima
  # A point is a local minimum if x[i-1] > x[i] < x[i+1]
  n <- length(x)
  is.min <- logical(n)

  for (i in 2:(n-1)) {
    is.min[i] <- x[i] < x[i-1] && x[i] < x[i+1]
  }

  k.values[is.min]
}

#' Pearson Weighted Correlation Coefficient
#'
#' Computes the Pearson correlation coefficient between two numeric vectors with weights.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as \code{x}.
#' @param w A numeric vector of non-negative weights of the same length as \code{x}.
#'
#' @return A single numeric value representing the weighted Pearson correlation coefficient.
#'
#' @details
#' This function computes the weighted Pearson correlation coefficient using C code
#' for efficiency. The weights must be non-negative and the three vectors must have
#' the same length.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' w <- runif(100)
#' pearson.wcor(x, y, w)
#'
#' @export
pearson.wcor <- function(x, y, w) {
    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(w)) {
        stop("x, y, and w must be numeric vectors")
    }

    n <- length(x)

    if (length(y) != n || length(w) != n) {
        stop("x, y, and w must have the same length")
    }

    if (any(w < 0)) {
        stop("weights must be non-negative")
    }

    cc <- 0
    out <- .C(C_pearson_wcor,
             as.double(x),
             as.double(y),
             as.double(w),
             as.integer(n),
             cc = as.double(cc))
    out$cc
}

#' Bayesian Bootstrap Credible Intervals for Pearson Weighted Correlation Coefficients
#'
#' Computes Bayesian bootstrap credible interval estimates of Pearson weighted
#' correlation coefficients for nearest neighbor data.
#'
#' @param nn.y1 A matrix of y1 values corresponding to nearest neighbors.
#'        Dimensions: \code{ng} x \code{K}, where \code{ng} is the number of groups
#'        and \code{K} is the number of nearest neighbors.
#' @param nn.y2 A matrix of y2 values corresponding to nearest neighbors.
#'        Same dimensions as \code{nn.y1}.
#' @param nn.i A matrix of nearest neighbor indices (1-based indexing).
#'        Same dimensions as \code{nn.y1}.
#' @param nn.w A matrix of nearest neighbor weights. Same dimensions as \code{nn.y1}.
#' @param nx The number of elements of x over which y1 and y2 are defined.
#' @param n.BB The number of Bayesian bootstrap iterations. Default is 1000.
#' @param alpha The significance level for credible intervals. Default is 0.05.
#'
#' @return A 2 x \code{ng} matrix where the first row contains the lower bounds
#'         and the second row contains the upper bounds of the
#'         \code{(1-alpha)*100\%} credible intervals for each group.
#'
#' @details
#' This internal function uses Bayesian bootstrap to compute credible intervals
#' for weighted correlation coefficients in a nearest neighbor context. The function
#' calls C code for computational efficiency. Note that the C function expects
#' 0-based indexing, so the R indices are adjusted by subtracting 1.
#'
#' @keywords internal
pearson.wcor.BB.qCrI <- function(nn.y1, nn.y2, nn.i, nn.w, nx, n.BB = 1000, alpha = 0.05) {
    # Input validation
    if (!is.matrix(nn.y1) || !is.matrix(nn.y2) || !is.matrix(nn.i) || !is.matrix(nn.w)) {
        stop("nn.y1, nn.y2, nn.i, and nn.w must be matrices")
    }

    K <- ncol(nn.w)
    ng <- nrow(nn.w)

    # Check dimensions consistency
    if (ncol(nn.y1) != K || nrow(nn.y1) != ng ||
        ncol(nn.y2) != K || nrow(nn.y2) != ng ||
        ncol(nn.i) != K || nrow(nn.i) != ng) {
        stop("All input matrices must have the same dimensions")
    }

    if (n.BB < 1) {
        stop("n.BB must be at least 1")
    }

    if (alpha <= 0 || alpha >= 1) {
        stop("alpha must be between 0 and 1")
    }

    lwcor.CI <- numeric(2 * ng)
    out <- .C(C_pearson_wcor_BB_qCrI,
             as.double(t(nn.y1)),
             as.double(t(nn.y2)),
             as.integer(t(nn.i - 1)),  # Convert to 0-based indexing for C
             as.double(t(nn.w)),
             as.integer(K),
             as.integer(ng),
             as.integer(nx),
             as.integer(n.BB),
             as.double(alpha),
             lwcor.CI = as.double(lwcor.CI))

    lwcor.CI <- matrix(out$lwcor.CI, nrow = 2, ncol = ng, byrow = FALSE)
    lwcor.CI
}

#' Remove outliers from a state space
#'
#' Removes outliers from a state space (SS) and optionally from a variable
#' defined over it using k-nearest neighbor distances.
#'
#' @param S A numeric matrix or data frame representing the state space,
#'   where rows are observations and columns are dimensions.
#' @param y An optional numeric vector of length \code{nrow(S)} representing
#'   a variable defined over the state space. Default is NULL.
#' @param p Numeric between 0 and 1. The percentile threshold for outlier removal.
#'   Default value of 0.98 removes samples whose distance to the nearest neighbor
#'   is greater than the 98th percentile of all nearest neighbor distances.
#' @param dist.factor A positive numeric constant used for outlier selection.
#'   For method "dist.factor", points with ratio of furthest to nearest neighbor
#'   distance exceeding this value are considered outliers. For method
#'   "diff.dist.factor", points with relative jump in neighbor distances
#'   exceeding this value are considered outliers. Default is 100.
#' @param K Integer specifying the number of nearest neighbors to use for
#'   outlier identification. Must be positive. Default is 30.
#' @param method Character string specifying the outlier detection method when K > 1.
#'   Options are:
#'   \itemize{
#'     \item "diff.dist.factor" (default): Uses the maximum jump in consecutive
#'       neighbor distances relative to the median jump.
#'     \item "dist.factor": Uses the ratio of the K-th to 1st neighbor distance.
#'     \item Any other value: Uses a threshold-based approach on all neighbor distances.
#'   }
#'
#' @return A list containing:
#'   \item{S.q}{The filtered state space with outliers removed}
#'   \item{y.q}{The filtered y variable (NULL if y was NULL)}
#'   \item{nn.d}{The K-nearest neighbor distance matrix from the original data}
#'   \item{d.thld}{The distance threshold used (p-th percentile)}
#'   \item{idx}{Logical vector indicating which observations were kept (TRUE) or
#'     removed as outliers (FALSE)}
#'
#' @importFrom FNN get.knn
#' @importFrom stats quantile median
#'
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(123)
#' S <- matrix(rnorm(200), ncol = 2)
#' y <- rnorm(100)
#'
#' # Remove outliers
#' result <- rm.SS.outliers(S, y, p = 0.95)
#'
#' # Check how many points were removed
#' sum(!result$idx)
#' }
#'
#' @seealso \code{\link[FNN]{get.knn}} for K-nearest neighbor computation
#' @export
rm.SS.outliers <- function(S, y = NULL, p = 0.98, dist.factor = 100, K = 30,
                          method = "diff.dist.factor") {
    # Validate inputs
    if (!is.matrix(S) && !is.data.frame(S)) {
        stop("S must be a matrix or data frame")
    }
    if (!is.null(y)) {
        if (!is.numeric(y)) {
            stop("y must be numeric or NULL")
        }
        if (length(y) != nrow(S)) {
            stop("Length of y must equal number of rows in S")
        }
    }
    stopifnot(is.numeric(K))
    stopifnot(as.integer(K) == K)
    stopifnot(K > 0)
    stopifnot(is.numeric(p))
    stopifnot(p > 0 & p <= 1)
    stopifnot(is.numeric(dist.factor))
    stopifnot(dist.factor > 0)

    # Get K nearest neighbors
    nn <- FNN::get.knn(S, k = K)
    nn.d <- nn$nn.dist
    d.nn <- nn.d[, 1]
    d.thld <- quantile(d.nn, probs = c(p))[[1]]

    # Determine outliers based on method
    if (K == 1) {
        idx <- d.nn < d.thld
    } else if (method == "dist.factor") {
        d.rat <- nn.d[, ncol(nn.d)] / d.nn
        idx <- d.rat < dist.factor
    } else if (method == "diff.dist.factor") {
        jump <- apply(nn.d, 1, function(x) max(diff(x)))
        jump.median <- median(jump)
        rel.jump <- jump / jump.median
        idx <- rel.jump < dist.factor
    } else {
        idx <- rep(TRUE, nrow(nn.d))
        for (i in seq_len(nrow(nn.d))) {
            d <- nn.d[i, ]
            d1 <- d[1]
            dd <- diff(d)
            idx1 <- d1 >= d.thld || sum(dd >= d.thld) > 0
            if (idx1) {
                idx[i] <- FALSE
            }
        }
    }

    # Filter data
    S.q <- S[idx, , drop = FALSE]
    y.q <- NULL
    if (!is.null(y)) {
        y.q <- y[idx]
    }

    list(S.q = S.q,
         y.q = y.q,
         nn.d = nn.d,
         d.thld = d.thld,
         idx = idx)
}

#' Apply Floor to Function Values
#'
#' @description
#' Sets all values below threshold to the threshold value.
#' Simple, explicit operation for enforcing lower bounds.
#'
#' @param y Numeric vector
#' @param floor.at Lower bound threshold
#' @param verbose Logical; print diagnostic info
#'
#' @return Numeric vector with floor applied
#'
#' @examples
#' y <- c(-0.1, 0.2, 0.5)
#' apply.floor(y, floor.at = 0)
#'
#' @export
apply.floor <- function(y, floor.at, verbose = FALSE) {
  n.floored <- sum(y < floor.at, na.rm = TRUE)

  if (verbose && n.floored > 0) {
    cat("Flooring", n.floored, "values to", floor.at, "\n")
    cat("  Original range: [", min(y, na.rm = TRUE), ",", max(y, na.rm = TRUE), "]\n")
  }

  y[!is.na(y) & y < floor.at] <- floor.at

  if (verbose && n.floored > 0) {
    cat("  New range: [", min(y, na.rm = TRUE), ",", max(y, na.rm = TRUE), "]\n")
  }

  return(y)
}


#' Apply Ceiling to Function Values
#'
#' @description
#' Sets all values above threshold to the threshold value.
#' Simple, explicit operation for enforcing upper bounds.
#'
#' @param y Numeric vector
#' @param ceiling.at Upper bound threshold
#' @param verbose Logical; print diagnostic info
#'
#' @return Numeric vector with ceiling applied
#'
#' @examples
#' y <- c(0.2, 0.5, 1.2)
#' apply.ceiling(y, ceiling.at = 1.0)
#'
#' @export
apply.ceiling <- function(y, ceiling.at, verbose = FALSE) {
  n.ceilinged <- sum(y > ceiling.at, na.rm = TRUE)

  if (verbose && n.ceilinged > 0) {
    cat("Ceiling", n.ceilinged, "values to", ceiling.at, "\n")
    cat("  Original range: [", min(y, na.rm = TRUE), ",", max(y, na.rm = TRUE), "]\n")
  }

  y[!is.na(y) & y > ceiling.at] <- ceiling.at

  if (verbose && n.ceilinged > 0) {
    cat("  New range: [", min(y, na.rm = TRUE), ",", max(y, na.rm = TRUE), "]\n")
  }

  return(y)
}


#' Left Winsorize (for completeness, though you may not use it)
#'
#' @description
#' Replaces extreme values on the left tail with p-th quantile.
#'
#' @param y Numeric vector
#' @param p Proportion for left tail (default 0.01)
#' @param verbose Logical; print diagnostic info
#'
#' @return Numeric vector with left tail winsorized
#'
#' @export
left.winsorize <- function(y, p = 0.01, verbose = FALSE) {
  if (p < 0 || p >= 0.25) {
    stop("p must be between 0 and 0.25")
  }
  if (p == 0) return(y)

  q.left <- quantile(y, probs = p, na.rm = TRUE)
  n.winsorized <- sum(y < q.left, na.rm = TRUE)

  if (verbose) {
    cat("Left winsorization (p =", p, "):\n")
    cat("  Threshold:", q.left, "\n")
    cat("  Values winsorized:", n.winsorized, "\n")
  }

  y[!is.na(y) & y < q.left] <- q.left

  return(y)
}


#' Right Winsorize
#'
#' @description
#' Replaces extreme values on the right tail with (1-p)-th quantile.
#'
#' @param y Numeric vector
#' @param p Proportion for right tail (default 0.01)
#' @param verbose Logical; print diagnostic info
#'
#' @return Numeric vector with right tail winsorized
#'
#' @examples
#' y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' right.winsorize(y, p = 0.1)
#'
#' @export
right.winsorize <- function(y, p = 0.01, verbose = FALSE) {
  if (p < 0 || p >= 0.25) {
    stop("p must be between 0 and 0.25")
  }
  if (p == 0) return(y)
  if (all(is.na(y))) return(y)

  q.right <- quantile(y, probs = 1 - p, na.rm = TRUE)
  n.winsorized <- sum(y > q.right, na.rm = TRUE)

  if (verbose) {
    cat("Right winsorization (p =", p, "):\n")
    cat("  Threshold:", q.right, "\n")
    cat("  Values winsorized:", n.winsorized, "\n")
  }

  y[!is.na(y) & y > q.right] <- q.right

  return(y)
}

#' Winsorize a numeric vector
#'
#' Replaces extreme values in a numeric vector with less extreme values.
#' Values below the p-th percentile are set to the p-th percentile,
#' and values above the (1-p)-th percentile are set to the (1-p)-th percentile.
#'
#' @param x A numeric vector to be winsorized
#' @param p The proportion of data to be winsorized on each tail (default is 0.01).
#'   Must be between 0 and 0.25.
#' @param verbose Set to TRUE to report the new min/max values.
#'
#' @return A numeric vector of the same length as \code{x} with extreme values replaced.
#'   NA values in the input are preserved in the output.
#'
#' @seealso \code{winsorize.zscore} for robust z-score winsorization
#'
#' @examples
#' # Simple example
#' x <- 1:10
#' winsorize(x, p = 0.2)
#'
#' # With normally distributed data
#' set.seed(123)
#' y <- rnorm(100)
#' y.wins <- winsorize(y, p = 0.05)
#'
#' # Compare ranges
#' range(y)
#' range(y.wins)
#'
#' @importFrom stats quantile
#' @export
winsorize <- function(x, p = 0.01, verbose = FALSE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }
  if (!is.numeric(p) || length(p) != 1) {
    stop("p must be a single numeric value")
  }
  if (p <= 0 || p >= 0.25) {
    stop("p must be greater than 0 and less than 0.25")
  }

  # Handle edge cases
  if (length(x) == 0) return(x)

  if (all(is.na(x))) {
    warning("All values in x are NA - no winsorization performed")
    return(x)
  }

  q <- quantile(x, probs = c(p, 1 - p), na.rm = TRUE)

  n.left <- sum(x < q[1], na.rm = TRUE)
  n.right <- sum(x > q[2], na.rm = TRUE)

  if (verbose) {
    cat("Two-sided winsorization (p =", p, "):\n")
    cat("  Left threshold:", q[1], "(", n.left, "values)\n")
    cat("  Right threshold:", q[2], "(", n.right, "values)\n")
  }

  x[!is.na(x) & x < q[1]] <- q[1]
  x[!is.na(x) & x > q[2]] <- q[2]

  return(x)
}

#' Break Ties with Adaptive Noise
#'
#' @description
#' Adds minimal noise to break ties in function values.
#' Automatically adapts noise scale to the data range.
#'
#' @param y Numeric vector
#' @param noise.scale Relative noise as fraction of range (default 1e-10)
#' @param min.abs.noise Minimum absolute noise magnitude (default 1e-12)
#' @param preserve.bounds Keep min/max exactly (default TRUE)
#' @param seed Random seed for reproducibility (default NULL)
#' @param verbose Logical; print diagnostic info
#'
#' @return Numeric vector with ties broken
#'
#' @details
#' For values at global min: perturbs upward only
#' For values at global max: perturbs downward only
#' For interior values: perturbs symmetrically
#'
#' @export
break.ties <- function(y,
                      noise.scale = 1e-10,
                      min.abs.noise = 1e-12,
                      preserve.bounds = TRUE,
                      seed = NULL,
                      verbose = FALSE) {

  if (!is.null(seed)) set.seed(seed)

  if (any(is.na(y))) {
    stop("y contains NA values")
  }

  # Get range
  y.min <- min(y)
  y.max <- max(y)
  y.range <- y.max - y.min

  if (y.range == 0) {
    stop("y is constant - cannot break ties")
  }

  # Check for ties
  tied.values <- unique(y[duplicated(y)])
  n.tied.values <- length(tied.values)

  if (n.tied.values == 0) {
    if (verbose) cat("No ties to break\n")
    return(y)
  }

  # Determine noise magnitude
  noise.magnitude <- max(y.range * noise.scale, min.abs.noise)

  if (verbose) {
    cat("Breaking ties:\n")
    cat("  Tied values:", n.tied.values, "\n")
    cat("  Noise magnitude:", format(noise.magnitude, scientific = TRUE), "\n")
    cat("  As % of range:", format(100 * noise.magnitude / y.range, digits = 4), "%\n")
  }

  # Process each tied value
  for (val in tied.values) {
    tied.indices <- which(y == val)
    n.tied <- length(tied.indices)

    is.global.min <- abs(val - y.min) < 1e-14
    is.global.max <- abs(val - y.max) < 1e-14

    # Generate perturbations
    if (preserve.bounds && is.global.min) {
      # Keep first at exact min, perturb others up
      perturbations <- c(0, sort(runif(n.tied - 1, 0, noise.magnitude)))

    } else if (preserve.bounds && is.global.max) {
      # Keep first at exact max, perturb others down
      perturbations <- c(0, -sort(runif(n.tied - 1, 0, noise.magnitude), decreasing = TRUE))

    } else {
      # Symmetric perturbation
      perturbations <- runif(n.tied, -noise.magnitude/2, noise.magnitude/2)
      perturbations <- perturbations - mean(perturbations)
    }

    y[tied.indices] <- val + perturbations
  }

  if (verbose) {
    cat("  Final unique values:", length(unique(y)), "/", length(y), "\n")
  }

  return(y)
}

#' Prepare Binary Outcome Conditional Expectation for Extrema Detection
#'
#' @description
#' Standard pipeline for `E[Y|X]` where `Y` is binary with values in `{0, 1}`.
#' Enforces `[0, 1]` bounds, applies right winsorization, and breaks ties.
#'
#' @param y.hat Estimated conditional expectation vector
#' @param p.right Right tail proportion for winsorization (default 0.01)
#' @param apply.right.winsorization Logical; whether to winsorize right tail
#' @param noise.scale Noise scale for tie breaking (default 1e-10)
#' @param seed Random seed (default 123)
#' @param verbose Logical; print diagnostics
#'
#' @return Numeric vector ready for extrema detection
#'
#' @examples
#' \dontrun{
#' # After fitting
#' sptb.cond.exp <- sptb.cond.exp.fit$predictions
#' sptb.cond.exp.clean <- prepare.binary.cond.exp(
#'   sptb.cond.exp,
#'   p.right = 0.01,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
prepare.binary.cond.exp <- function(y.hat,
                                   p.right = 0.01,
                                   apply.right.winsorization = TRUE,
                                   noise.scale = 1e-10,
                                   seed = 123,
                                   verbose = FALSE) {

  if (verbose) {
    cat("\n=== PREPARING BINARY CONDITIONAL EXPECTATION ===\n")
    cat("Original range: [", min(y.hat), ",", max(y.hat), "]\n")
  }

  # Step 1: Floor at 0 (conditional expectation cannot be negative)
  y.clean <- apply.floor(y.hat, floor.at = 0, verbose = verbose)

  # Step 2: Right winsorization (if requested)
  if (apply.right.winsorization) {
    y.clean <- right.winsorize(y.clean, p = p.right, verbose = verbose)
  }

  # Step 3: Ceiling at 1 (conditional expectation cannot exceed 1)
  # Note: This might do nothing if data is left-skewed
  y.clean <- apply.ceiling(y.clean, ceiling.at = 1, verbose = verbose)

  # Step 4: Break ties
  y.clean <- break.ties(y.clean,
                       noise.scale = noise.scale,
                       preserve.bounds = TRUE,
                       seed = seed,
                       verbose = verbose)

  if (verbose) {
    cat("Final range: [", min(y.clean), ",", max(y.clean), "]\n")
    cat("=== PREPARATION COMPLETE ===\n\n")
  }

  return(y.clean)
}


#' Prepare Continuous Outcome for Extrema Detection
#'
#' @description
#' Standard pipeline for continuous outcomes.
#' Applies two-sided winsorization and breaks ties.
#'
#' @param y Continuous outcome vector
#' @param p.winsorize Proportion for two-sided winsorization (default 0.01)
#' @param apply.winsorization Logical; whether to winsorize (default TRUE)
#' @param noise.scale Noise scale for tie breaking (default 1e-10)
#' @param seed Random seed (default 123)
#' @param verbose Logical; print diagnostics
#'
#' @return Numeric vector ready for extrema detection
#'
#' @examples
#' \dontrun{
#' y.clean <- prepare.continuous.outcome(
#'   y.raw,
#'   p.winsorize = 0.025,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
prepare.continuous.outcome <- function(y,
                                      p.winsorize = 0.01,
                                      apply.winsorization = TRUE,
                                      noise.scale = 1e-10,
                                      seed = 123,
                                      verbose = FALSE) {

  if (verbose) {
    cat("\n=== PREPARING CONTINUOUS OUTCOME ===\n")
    cat("Original range: [", min(y), ",", max(y), "]\n")
  }

  y.clean <- y

  # Step 1: Two-sided winsorization
  if (apply.winsorization) {
    y.clean <- winsorize(y.clean, p = p.winsorize, verbose = verbose)
  }

  # Step 2: Break ties
  y.clean <- break.ties(y.clean,
                       noise.scale = noise.scale,
                       preserve.bounds = TRUE,
                       seed = seed,
                       verbose = verbose)

  if (verbose) {
    cat("Final range: [", min(y.clean), ",", max(y.clean), "]\n")
    cat("=== PREPARATION COMPLETE ===\n\n")
  }

  return(y.clean)
}

#' Prepare for Extrema Detection (Fully Customizable)
#'
#' @description
#' Flexible wrapper that allows manual specification of all preprocessing steps.
#' Use this when you need full control or have already applied some steps.
#'
#' @param y Function values (possibly already preprocessed)
#' @param floor.at Optional floor value (NULL = no flooring)
#' @param ceiling.at Optional ceiling value (NULL = no ceiling)
#' @param left.winsorize.p Optional left winsorization proportion (NULL = skip)
#' @param right.winsorize.p Optional right winsorization proportion (NULL = skip)
#' @param break.ties.after Logical; break ties at end (default TRUE)
#' @param noise.scale Noise scale for tie breaking (default 1e-10)
#' @param seed Random seed (default 123)
#' @param verbose Logical; print diagnostics
#'
#' @return Numeric vector ready for extrema detection
#'
#' @examples
#' \dontrun{
#' # Your current workflow:
#' sptb.cond.exp <- sptb.cond.exp.fit$predictions
#' sptb.cond.exp <- ifelse(sptb.cond.exp < 0, 0, sptb.cond.exp)  # manual floor
#' sptb.cond.exp <- right.winsorize(sptb.cond.exp)
#'
#' # Then just break ties:
#' sptb.cond.exp.clean <- prepare.for.extrema.detection(
#'   sptb.cond.exp,
#'   floor.at = NULL,  # Already done manually
#'   ceiling.at = NULL,  # Not needed
#'   break.ties.after = TRUE,
#'   verbose = FALSE
#' )
#'
#' # Or let the function do everything:
#' sptb.cond.exp.clean <- prepare.for.extrema.detection(
#'   sptb.cond.exp.fit$predictions,
#'   floor.at = 0,
#'   right.winsorize.p = 0.01,
#'   ceiling.at = 1,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
prepare.for.extrema.detection <- function(y,
                                         floor.at = NULL,
                                         ceiling.at = NULL,
                                         left.winsorize.p = NULL,
                                         right.winsorize.p = NULL,
                                         break.ties.after = TRUE,
                                         noise.scale = 1e-10,
                                         seed = 123,
                                         verbose = FALSE) {

  if (verbose) {
    cat("\n=== CUSTOM PREPARATION PIPELINE ===\n")
    cat("Original range: [", min(y), ",", max(y), "]\n")
  }

  y.clean <- y

  # Apply operations in logical order

  # 1. Floor
  if (!is.null(floor.at)) {
    y.clean <- apply.floor(y.clean, floor.at = floor.at, verbose = verbose)
  }

  # 2. Left winsorization
  if (!is.null(left.winsorize.p)) {
    y.clean <- left.winsorize(y.clean, p = left.winsorize.p, verbose = verbose)
  }

  # 3. Right winsorization
  if (!is.null(right.winsorize.p)) {
    y.clean <- right.winsorize(y.clean, p = right.winsorize.p, verbose = verbose)
  }

  # 4. Ceiling
  if (!is.null(ceiling.at)) {
    y.clean <- apply.ceiling(y.clean, ceiling.at = ceiling.at, verbose = verbose)
  }

  # 5. Break ties
  if (break.ties.after) {
    y.clean <- break.ties(y.clean,
                         noise.scale = noise.scale,
                         preserve.bounds = TRUE,
                         seed = seed,
                         verbose = verbose)
  }

  if (verbose) {
    cat("Final range: [", min(y.clean), ",", max(y.clean), "]\n")
    cat("=== PREPARATION COMPLETE ===\n\n")
  }

  return(y.clean)
}

#' Break Ties in Compositional Data Matrix
#'
#' Resolves duplicate rows in a relative abundance or compositional data matrix
#' by adding small, locally-informed perturbations that respect the structure of
#' nearby communities. The perturbation is proportional to local variance estimated
#' from a neighborhood defined either by k-nearest neighbors or by distance radius.
#'
#' @param rel.abund.mat Numeric matrix (n x p) of relative abundances where rows
#'   represent samples and columns represent taxa/features. Rows should sum to 1
#'   (or close to 1 within numerical precision). Can contain zeros.
#' @param neighborhood.method Character string specifying how to define the
#'   neighborhood for variance estimation. Either \code{"knn"} for k-nearest
#'   neighbors or \code{"radius"} for distance-based threshold. Default is
#'   \code{"knn"}.
#' @param neighborhood.size Integer specifying the number of nearest neighbors
#'   to use when \code{neighborhood.method = "knn"}. Default is 20.
#' @param neighborhood.radius Numeric value specifying the distance threshold
#'   when \code{neighborhood.method = "radius"}. Only samples within this
#'   distance are included in the neighborhood. Default is 0.01.
#' @param distance.metric Character string specifying the distance metric to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"euclidean"}: L² norm, \eqn{d(x,y) = \sqrt{\sum(x_i - y_i)^2}}
#'     \item \code{"manhattan"}: L¹ norm, \eqn{d(x,y) = \sum|x_i - y_i|}
#'     \item \code{"chebyshev"}: L∞ norm, \eqn{d(x,y) = \max_i |x_i - y_i|}
#'     \item \code{"bray.curtis"}: \eqn{d(x,y) = \sum|x_i - y_i| / \sum(x_i + y_i)}
#'   }
#'   Default is \code{"euclidean"}.
#' @param noise.scale Numeric value controlling the magnitude of perturbation.
#'   The noise added to each taxon is drawn from \eqn{N(0, \sigma_i^2 \cdot
#'   \text{noise.scale}^2)} where \eqn{\sigma_i} is the local standard deviation
#'   for taxon i. Use very small values (e.g., 1e-10 to 1e-8) to minimally
#'   affect downstream analyses. Default is 1e-10.
#' @param min.neighborhood.size Integer specifying the minimum number of neighbors
#'   required for valid variance estimation when using \code{neighborhood.method =
#'   "radius"}. If fewer neighbors are found within the radius, the function falls
#'   back to k-NN with this value. Default is 5.
#' @param seed Integer seed for random number generation to ensure reproducibility.
#'   If \code{NULL}, no seed is set. Default is \code{NULL}.
#' @param verbose Logical indicating whether to print diagnostic information
#'   including number of duplicates, neighborhood statistics, and warnings.
#'   Default is \code{FALSE}.
#'
#' @return A numeric matrix with the same dimensions as \code{rel.abund.mat}
#'   where duplicate rows have been perturbed to be unique. Rows are guaranteed
#'   to sum to 1 (within numerical precision) and all values are non-negative.
#'
#' @details
#' The function identifies duplicate rows and perturbs them by adding noise
#' proportional to the local compositional variance. For each duplicate:
#' \enumerate{
#'   \item A neighborhood is identified using either k-NN or distance radius
#'   \item Local variance is estimated for each taxon from the neighborhood
#'   \item Gaussian noise scaled by local variance is added to each taxon
#'   \item Negative values are set to zero
#'   \item The composition is re-normalized to sum to 1
#' }
#'
#' For taxa with zero local variance, the function uses a fallback strategy:
#' first trying the minimum positive local variance scaled by 0.1, then falling
#' back to global variance if necessary.
#'
#' When \code{neighborhood.method = "radius"} and too few neighbors are found
#' (< \code{min.neighborhood.size}), the function automatically falls back to
#' k-NN with \code{min.neighborhood.size} neighbors and reports this in verbose
#' mode.
#'
#' @section Distance Metric Guidance:
#' \itemize{
#'   \item \strong{Manhattan (L¹)}: Sensitive to cumulative differences across
#'     all taxa. Range: 0 to 2 for compositions.
#'   \item \strong{Euclidean (L²)}: Standard geometric distance, penalizes large
#'     single-taxon differences. Range: 0 to sqrt(2) for compositions.
#'   \item \strong{Chebyshev (L∞)}: Focuses on maximum single-taxon difference.
#'     Useful when you want uniform bounds on all taxa. Range: 0 to 1 for
#'     compositions.
#'   \item \strong{Bray-Curtis}: Normalized dissimilarity common in ecology,
#'     scale-invariant. Range: 0 to 1.
#' }
#'
#' @section Choosing Parameters:
#' \itemize{
#'   \item For high-dimensional sparse data (e.g., microbiome): Use larger
#'     \code{neighborhood.size} (30-50) or larger \code{neighborhood.radius}
#'     (0.05-0.1) to ensure sufficient neighbors.
#'   \item For dense regions (e.g., samples dominated by single taxon): Use
#'     smaller \code{neighborhood.radius} (0.001-0.01) to capture local structure.
#'   \item Start with \code{noise.scale = 1e-10} and increase only if duplicates
#'     remain after perturbation.
#' }
#'
#' @examples
#' # Simulate compositional data with duplicates
#' set.seed(42)
#' n <- 100
#' p <- 50
#' X <- matrix(rexp(n * p), nrow = n, ncol = p)
#' X <- X / rowSums(X)  # normalize to compositions
#'
#' # Introduce some duplicates
#' X[5, ] <- X[1, ]
#' X[10, ] <- X[1, ]
#' X[20, ] <- X[15, ]
#'
#' sum(duplicated(X))  # Should be 3
#'
#' # Break ties using k-NN with Euclidean distance
#' X.dedup.knn <- break.composition.ties(
#'   rel.abund.mat = X,
#'   neighborhood.method = "knn",
#'   neighborhood.size = 20,
#'   distance.metric = "euclidean",
#'   noise.scale = 1e-10,
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' sum(duplicated(X.dedup.knn))  # Should be 0
#'
#' # Break ties using radius method with Manhattan distance
#' X.dedup.radius <- break.composition.ties(
#'   rel.abund.mat = X,
#'   neighborhood.method = "radius",
#'   neighborhood.radius = 0.01,
#'   distance.metric = "manhattan",
#'   noise.scale = 1e-10,
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' # For L. iners dominated samples (restrictive neighborhood)
#' X.dedup.linf <- break.composition.ties(
#'   rel.abund.mat = X,
#'   neighborhood.method = "radius",
#'   neighborhood.radius = 0.01,
#'   distance.metric = "chebyshev",  # L∞ ensures all taxa within bounds
#'   noise.scale = 1e-10,
#'   min.neighborhood.size = 10,
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' @seealso
#' \code{\link{duplicated}} for identifying duplicate rows,
#' \code{\link{prcomp}} for PCA-based alternative approaches
#'
#' @export
break.composition.ties <- function(
  rel.abund.mat,
  neighborhood.method = c("knn", "radius"),
  neighborhood.size = 20,
  neighborhood.radius = 0.01,
  distance.metric = c("euclidean", "manhattan", "chebyshev", "bray.curtis"),
  noise.scale = 1e-10,
  min.neighborhood.size = 5,
  seed = NULL,
  verbose = FALSE
) {
  # rel.abund.mat: n x p matrix of relative abundances (rows sum to 1)
  # neighborhood.method: "knn" for k-nearest neighbors, "radius" for distance threshold
  # neighborhood.size: k for knn method
  # neighborhood.radius: distance threshold for radius method
  # distance.metric: "euclidean" (L2), "manhattan" (L1), "chebyshev" (L∞), or "bray.curtis"
  # noise.scale: multiplier for perturbation magnitude
  # min.neighborhood.size: minimum neighbors required for valid variance estimate

  if (!is.null(seed)) set.seed(seed)

  neighborhood.method <- match.arg(neighborhood.method)
  distance.metric <- match.arg(distance.metric)

  n <- nrow(rel.abund.mat)
  p <- ncol(rel.abund.mat)

  # Identify duplicate rows
  dup.rows <- duplicated(rel.abund.mat) | duplicated(rel.abund.mat, fromLast = TRUE)
  n.dup <- sum(dup.rows)

  if (verbose) {
    cat("Total samples:", n, "\n")
    cat("Duplicate samples:", n.dup, "\n")
    cat("Neighborhood method:", neighborhood.method, "\n")
    cat("Distance metric:", distance.metric, "\n")
  }

  if (n.dup == 0) {
    if (verbose) cat("No duplicates found\n")
    return(rel.abund.mat)
  }

  # Distance function
  compute.distance <- function(x, y, metric) {
    switch(metric,
      euclidean = sqrt(sum((x - y)^2)),
      manhattan = sum(abs(x - y)),
      chebyshev = max(abs(x - y)),
      bray.curtis = {
        # Bray-Curtis dissimilarity: sum(|x_i - y_i|) / sum(x_i + y_i)
        sum(abs(x - y)) / sum(x + y)
      }
    )
  }

  result <- rel.abund.mat
  neighborhood.stats <- list()

  for (i in which(dup.rows)) {
    # Compute distances to all other samples
    distances <- apply(rel.abund.mat, 1, function(x) {
      compute.distance(x, rel.abund.mat[i, ], distance.metric)
    })

    # Select neighborhood based on method
    if (neighborhood.method == "knn") {
      # k nearest neighbors (excluding self at distance 0)
      neighbor.idx <- order(distances)[2:(neighborhood.size + 1)]
      max.dist <- max(distances[neighbor.idx])

    } else {  # radius method
      # All samples within radius (excluding self)
      neighbor.idx <- which(distances > 0 & distances <= neighborhood.radius)
      max.dist <- neighborhood.radius

      # Check if we have enough neighbors
      if (length(neighbor.idx) < min.neighborhood.size) {
        if (verbose) {
          cat(sprintf("Sample %d: Only %d neighbors within radius %.4f, ",
                     i, length(neighbor.idx), neighborhood.radius))
        }

        # Fallback to knn with min.neighborhood.size
        neighbor.idx <- order(distances)[2:(min.neighborhood.size + 1)]
        max.dist <- max(distances[neighbor.idx])

        if (verbose) {
          cat(sprintf("using %d-NN fallback (max dist: %.4f)\n",
                     min.neighborhood.size, max.dist))
        }
      }
    }

    neighborhood <- rel.abund.mat[neighbor.idx, , drop = FALSE]

    if (verbose) {
      neighborhood.stats[[as.character(i)]] <- list(
        n.neighbors = length(neighbor.idx),
        max.distance = max.dist,
        mean.distance = mean(distances[neighbor.idx])
      )
    }

    # Estimate local variance for each taxon
    local.sd <- apply(neighborhood, 2, sd)

    # For taxa with zero local variance, use fallback
    zero.var.taxa <- local.sd == 0
    if (any(zero.var.taxa)) {
      min.positive.sd <- min(local.sd[!zero.var.taxa], na.rm = TRUE)
      if (is.finite(min.positive.sd)) {
        local.sd[zero.var.taxa] <- min.positive.sd * 0.1
      } else {
        # Ultimate fallback: use global variance
        global.sd <- apply(rel.abund.mat, 2, sd)
        local.sd[zero.var.taxa] <- global.sd[zero.var.taxa] * 0.1
      }
    }

    # Add scaled noise proportional to local variance
    noise <- rnorm(p, mean = 0, sd = local.sd * noise.scale)

    # Apply perturbation
    perturbed <- rel.abund.mat[i, ] + noise

    # Ensure non-negativity
    perturbed[perturbed < 0] <- 0

    # Re-normalize to sum to 1 (preserve compositional constraint)
    if (sum(perturbed) > 0) {
      result[i, ] <- perturbed / sum(perturbed)
    } else {
      # Extremely rare case: add uniform tiny noise
      result[i, ] <- (rel.abund.mat[i, ] + 1e-15) / (1 + p * 1e-15)
    }
  }

  # Verify all duplicates are resolved
  remaining.dup <- sum(duplicated(result))

  if (verbose) {
    cat("\nNeighborhood statistics:\n")
    for (sample.id in names(neighborhood.stats)) {
      stats <- neighborhood.stats[[sample.id]]
      cat(sprintf("  Sample %s: %d neighbors, max dist: %.4f, mean dist: %.4f\n",
                 sample.id, stats$n.neighbors, stats$max.distance, stats$mean.distance))
    }
    cat("\nRemaining duplicates after perturbation:", remaining.dup, "\n")
    if (remaining.dup > 0) {
      cat("WARNING: Some duplicates remain. Consider increasing noise.scale\n")
    }
  }

  return(result)
}

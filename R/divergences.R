
#' Estimate Relative Entropy Using Nearest Neighbor Distance Ratio
#'
#' @description
#' This function estimates the relative entropy (Kullback-Leibler divergence) between two datasets
#' using the nearest neighbor distance ratio method. The intuition behind this approach is that
#' if two distributions are similar, the ratio of distances to nearest neighbors within a set
#' and to the other set should be close to 1. If the distributions differ, this ratio will deviate from 1.
#'
#' The algorithm works by:
#' 1. Finding the k-nearest neighbors for each point within its own set and in the other set.
#' 2. Computing the ratios of these distances.
#' 3. Using the log of these ratios to estimate the relative entropy.
#'
#' This method is particularly useful in high-dimensional spaces where traditional density
#' estimation techniques may fail due to the curse of dimensionality.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#' @param k The number of nearest neighbors to consider (default is 1).
#'
#' @return A numeric value representing the estimated relative entropy.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- nn.distance.ratio.estimator(X, Y)
#' print(result)
#'
#' @references
#' Wang, Q., Kulkarni, S. R., & Verdu, S. (2009). Divergence estimation for multidimensional densities
#' via k-nearest-neighbor distances. IEEE Transactions on Information Theory, 55(5), 2392-2405.
init.version.nn.distance.ratio.estimator <- function(X, Y, k = 1) {
  n <- nrow(X)

  # Find k+1 nearest neighbors for each point in X and Y
  nn.X <- get.knn(X, k = k + 1)
  nn.Y <- get.knn(Y, k = k + 1)

  # Compute distances to k-th nearest neighbor within each set
  d.X <- nn.X$nn.dist[, k]
  d.Y <- nn.Y$nn.dist[, k]

  # Compute distances to k-th nearest neighbor in the other set
  d.XY <- get.knn(rbind(X, Y), k = k)$nn.dist[1:n, k]
  d.YX <- get.knn(rbind(Y, X), k = k)$nn.dist[1:n, k]

  # Estimate relative entropy
  relative.entropy <- mean(log(d.XY / d.X)) + mean(log(d.YX / d.Y))

  return(relative.entropy)
}

#' Estimate Relative Entropy Using Nearest Neighbor Distance Ratio
#'
#' @description
#' This function estimates the relative entropy (Kullback-Leibler divergence) between two datasets
#' using the nearest neighbor distance ratio method. It includes safeguards against zero distances
#' and potential numerical instabilities.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#' @param k The number of nearest neighbors to consider (default is 1).
#' @param eps A small constant to add to distances to avoid numerical issues (default is NULL, which means it will be automatically determined).
#' @param eps.factor A small factor used to compute eps when eps is NULL. eps is set to the smallest non-zero distance multiplied by this factor (default is 1e-8).
#'
#' @return A numeric value representing the estimated relative entropy.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- nn.distance.ratio.estimator(X, Y)
#' print(result)
#'
#' @references
#' Wang, Q., Kulkarni, S. R., & Verdu, S. (2009). Divergence estimation for multidimensional densities
#' via k-nearest-neighbor distances. IEEE Transactions on Information Theory, 55(5), 2392-2405.
#' @export
nn.distance.ratio.estimator <- function(X, Y, k = 1, eps = NULL, eps.factor = 1e-8) {

    n <- nrow(X)

    ## Finding k+1 nearest neighbors for each point in X and Y
    nn.X <- FNN::get.knn(X, k = k + 1)
    nn.Y <- FNN::get.knn(Y, k = k + 1)

    ## Computing distances to k-th nearest neighbor within each set
    d.X <- nn.X$nn.dist[, k]
    d.Y <- nn.Y$nn.dist[, k]

    ## Computing distances to k-th nearest neighbor in the other set
    d.XY <- FNN::get.knn(rbind(X, Y), k = k)$nn.dist[1:n, k]
    d.YX <- FNN::get.knn(rbind(Y, X), k = k)$nn.dist[1:n, k]

    ## Determining eps if not provided
    if (is.null(eps)) {
        min.non.zero <- min(c(d.X[d.X > 0], d.Y[d.Y > 0], d.XY[d.XY > 0], d.YX[d.YX > 0]))
        eps <- min.non.zero * eps.factor  # A small fraction of the smallest non-zero distance
    }

    ## Adding eps to all distances
    d.X <- d.X + eps
    d.Y <- d.Y + eps
    d.XY <- d.XY + eps
    d.YX <- d.YX + eps

    ## Estimating relative entropy
    relative.entropy <- mean(log(d.XY / d.X)) + mean(log(d.YX / d.Y))

    ## return(max(0, relative.entropy))  # Ensure non-negativity
    return(relative.entropy)
}


nn.divergence <- function(X, Y, k = 1, eps = NULL, eps.factor = 1e-8) {

    n <- nrow(X)

    ## Computing distances to k-th nearest neighbor within each set
    d.X <- FNN::get.knn(X, k)$nn.dist[, k]

    ## Computing distances to k-th nearest neighbor in the other set
    d.YX <- FNN::get.knnx(Y, X, k)$nn.dist[, k]

    ## Determining eps if not provided
    if (is.null(eps)) {
        min.non.zero <- min(c(d.X[d.X > 0], d.YX[d.YX > 0]))
        eps <- min.non.zero * eps.factor  # A small fraction of the smallest non-zero distance
    }

    ## Adding eps to all distances
    d.X <- d.X + eps
    d.YX <- d.YX + eps

    ## Estimating relative entropy
    relative.entropy <- mean(log(d.YX / d.X)) + log((n - 1) / n)

    return(relative.entropy)
}


#' Compute Wasserstein Distance Between Two Datasets
#'
#' @description
#' This function calculates the Wasserstein distance (also known as Earth Mover's Distance)
#' between two datasets. The Wasserstein distance can be intuitively understood as the minimum
#' "cost" of transforming one distribution into another, where cost is measured as the amount
#' of probability mass that needs to be moved, multiplied by the distance it needs to be moved.
#'
#' Imagine each distribution as a pile of earth, and the Wasserstein distance as the minimum
#' amount of work required to transform one pile into the other. This makes it particularly
#' useful for comparing distributions with different supports or when KL-divergence might be
#' undefined or infinite.
#'
#' The algorithm uses the transport package to compute the distance efficiently. It's worth
#' noting that while this isn't a direct proxy for relative entropy, it provides a robust
#' measure of dissimilarity between distributions.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#'
#' @return A numeric value representing the Wasserstein distance between X and Y.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- wasserstein.distance(X, Y)
#' print(result)
#'
#' @references
#' Villani, C. (2008). Optimal transport: old and new (Vol. 338). Springer Science & Business Media.
#' @export
wasserstein.distance <- function(X, Y) {
  # Ensure X and Y are matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # Create point patterns from the data
  # Assuming equal weights for all points
  X_pp <- transport::pp(X)
  Y_pp <- transport::pp(Y)

  # Compute Wasserstein distance
  w <- transport::wasserstein(X_pp, Y_pp, p = 2)

  return(w)
}

#' Compute Wasserstein Divergence Between Two Point Sets
#'
#' This function calculates a modified Wasserstein divergence between two point sets X and Y.
#' It uses k-nearest neighbors in X to define local neighborhoods, and computes the
#' Wasserstein distance between corresponding neighborhoods in X and Y.
#'
#' @param X A numeric matrix where each row represents a point in n-dimensional space.
#' @param Y A numeric matrix with the same dimensions as X, representing the second point set.
#' @param k An integer specifying the number of nearest neighbors to consider.
#'
#' @return A numeric value representing the sum of Wasserstein distances between local neighborhoods.
#'
#' @details
#' The function first checks the validity of inputs, then uses k-nearest neighbors to define
#' local neighborhoods in X. For each neighborhood, it computes the 2-Wasserstein distance
#' between the corresponding points in X and Y, using Euclidean ground distance.
#' The final divergence is the sum of these local Wasserstein distances.
#'
#' @note
#' This function requires the 'FNN' package for k-nearest neighbor calculations and
#' the 'transport' package for Wasserstein distance computations.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(1000 * 3), nrow = 1000, ncol = 3)
#' Y <- matrix(rnorm(1000 * 3), nrow = 1000, ncol = 3)
#' div <- wasserstein.divergence(X, Y, k = 10)
#' }
#'
#' @importFrom FNN get.knn
#' @importFrom transport wasserstein
#'
#' @export
wasserstein.divergence <- function(X, Y, k) {

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

    if (!is.matrix(Y)) {
        Y <- try(as.matrix(Y), silent = TRUE)
        if (inherits(Y, "try-error")) {
            stop("Y must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(Y)) {
        stop("Y must contain numeric values")
    }

    if (any(is.na(Y)) || any(is.infinite(Y))) {
        stop("Y cannot contain NA, NaN, or Inf values")
    }

    if (nrow(X) != nrow(Y)) {
        stop("X and Y have to have the same number of rows.")
    }

    if (!is.numeric(k) || k != round(k) || k < 1) {
        stop("k has to be a positive integer.")
    }

    nn <- FNN::get.knn(X, k)
    nn.i <- nn$nn.index
    dW <- 0
    for (i in seq(nrow(X))) {
        X_subset <- X[nn.i[i,], , drop = FALSE]
        Y_subset <- Y[nn.i[i,], , drop = FALSE]

        ## Creating cost matrix
        cost_matrix <- as.matrix(dist(rbind(X_subset, Y_subset), method = "euclidean"))

        ## Creating uniform weights
        wx <- rep(1/nrow(X_subset), nrow(X_subset))
        wy <- rep(1/nrow(Y_subset), nrow(Y_subset))

        dW <- dW + transport::wasserstein(wx, wy, costm = cost_matrix, p = 2)
    }

    return(dW / nrow(X))
}

#' Compute Angular Wasserstein Index between two point sets
#'
#' This function calculates the Angular Wasserstein Index I_W(X, Y) between two point sets X and Y.
#' It computes the angular distribution of k-nearest neighbors for each point and then
#' calculates the Wasserstein distance between these angular distributions.
#'
#' @param X A numeric matrix where each row represents a point in n-dimensional space.
#' @param Y A numeric matrix with the same dimensions as X, representing the second point set.
#' @param k An integer specifying the number of nearest neighbors to consider.
#'
#' @return A numeric value representing the Angular Wasserstein Index.
#'
#' @importFrom FNN get.knn
#' @importFrom transport wasserstein
#' @export
angular.wasserstein.index <- function(X, Y, k) {

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

    if (!is.matrix(Y)) {
        Y <- try(as.matrix(Y), silent = TRUE)
        if (inherits(Y, "try-error")) {
            stop("Y must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(Y)) {
        stop("Y must contain numeric values")
    }

    if (any(is.na(Y)) || any(is.infinite(Y))) {
        stop("Y cannot contain NA, NaN, or Inf values")
    }

    if (nrow(X) != nrow(Y)) {
        stop("X and Y have to have the same number of rows.")
    }

    if (!is.numeric(k) || k != round(k) || k < 2) {
        stop("k has to be an integer greater than 1.")
    }

    ## Compute Angles Between Vectors and a Reference Vector
    ##
    ## This function calculates the angles between a set of vectors and a reference vector.
    ## It handles the edge case of zero vectors and a zero reference vector.
    ##
    ## @param vectors A matrix where each row represents a vector.
    ## @param reference A vector to which the angles are computed.
    ##
    ## @return A numeric vector of angles (in radians) between each row of 'vectors'
    ##         and the 'reference' vector.
    ##
    ## @details
    ## The function computes the angle between each row of 'vectors' and the 'reference'
    ## vector using the dot product method. It handles several edge cases:
    ##
    ## 1. Zero vectors in 'vectors':
    ##    - If a row in 'vectors' is a zero vector (all elements are zero), the angle
    ##      is set to pi/2 (90 degrees). This assumes that a zero vector is
    ##      perpendicular to any other vector.
    ##
    ## 2. Zero reference vector:
    ##    - If 'reference' is a zero vector, the function returns NA for all angles,
    ##      as the angle is undefined in this case.
    ##
    ## 3. Numerical precision:
    ##    - To account for potential floating-point precision issues, the cosine
    ##      values are clamped to the range [-1, 1] before computing the arccosine.
    ##
    ## @note
    ## - The choice of pi/2 for zero vectors is arbitrary and can be adjusted based
    ##   on the specific needs of the application.
    ## - Consider logging or warning when zero vectors are encountered, as they may
    ##   indicate unexpected data or preprocessing issues.
    ##
    ## @examples
    ## vectors <- matrix(c(1, 0, 0,
    ##                     0, 1, 0,
    ##                     0, 0, 0), nrow = 3, byrow = TRUE)
    ## reference <- c(1, 1, 0)
    ## angles <- compute.angles(vectors, reference)
    ## print(angles)
    ## # Expected output: c(pi/4, pi/4, pi/2)
    ##
    ## @seealso
    ## \code{\link{acos}}, \code{\link{sqrt}}
    ##
    ## vector is a (k - 1)-by-ncol(X) matrix, whereas 'reference' is a vector of length ncol(X)
    compute.angles <- function(vectors, reference) {
        dot.products <- vectors %*% reference
        vector_norms <- sqrt(rowSums(vectors^2))
        ref_norm <- sqrt(sum(reference^2))

        ## Handle zero vectors
        zero_vectors <- vector_norms == 0
        zero_ref <- ref_norm == 0

        if (zero_ref) {
            ## If reference is zero, return NA or a default angle
            return(rep(NA_real_, nrow(vectors)))
        }

        ## Compute cosine of angles
        cos_angles <- rep(0, nrow(vectors))

        ## For non-zero vectors, compute as usual
        non_zero <- !zero_vectors
        cos_angles[non_zero] <- dot.products[non_zero] / (vector_norms[non_zero] * ref_norm)

        ## Clamp values to [-1, 1] to account for potential numerical issues
        cos_angles <- pmin(pmax(cos_angles, -1), 1)

        ## Compute angles
        angles <- acos(cos_angles)

        ## Optionally, assign a specific value for zero vectors (e.g., 0 or pi/2)
        angles[zero_vectors] <- pi/2  ## or another chosen value

        return(angles)
    }

    ## Rest of the function remains the same
    nn <- FNN::get.knn(X, k)
    nn.i <- nn$nn.index

    total.dist <- 0

    for (i in 1:nrow(X)) {
        neighbors.X <- X[nn.i[i, ], ] - matrix(X[i, ], nrow = k, ncol = ncol(X), byrow = TRUE)
        neighbors.Y <- Y[nn.i[i, ], ] - matrix(Y[i, ], nrow = k, ncol = ncol(Y), byrow = TRUE)

        reference.X <- neighbors.X[k, ]
        reference.Y <- neighbors.Y[k, ]

        angles.X <- compute.angles(neighbors.X[-k, ], reference.X)
        angles.Y <- compute.angles(neighbors.Y[-k, ], reference.Y)

        if (any(is.na(angles.X))) {
            cat("angles.X\n")
            print(angles.X)
            stop("i:",i," angles.X has NAs")
        }

        W <- 0
        out <- .C(C_wasserstein_distance_1D,
                  as.double(angles.X),
                  as.double(angles.Y),
                  as.integer(k-1),
                  W=as.double(W))

        total.dist <- total.dist + out$W
    }

    return(total.dist / nrow(X))
}


cpp.angular.wasserstein.index <- function(X, Y, k) {

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

    if (!is.matrix(Y)) {
        Y <- try(as.matrix(Y), silent = TRUE)
        if (inherits(Y, "try-error")) {
            stop("Y must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(Y)) {
        stop("Y must contain numeric values")
    }

    if (any(is.na(Y)) || any(is.infinite(Y))) {
        stop("Y cannot contain NA, NaN, or Inf values")
    }

    if (nrow(X) != nrow(Y)) {
        stop("X and Y have to have the same number of rows.")
    }

    if (!is.numeric(k) || k != round(k) || k < 2) {
        stop("k has to be an integer greater than 1.")
    }

    ## Performance warning for large matrices
    if (nrow(X) * ncol(X) > 1e6) {
        warning("Large matrices detected. Computation may take a while.")
    }

    result <- .Call(S_angular_wasserstein_index, X, Y, as.integer(k + 1))

    return(result)
}

# Energy Distance
#' Compute Energy Distance Between Two Datasets
#'
#' @description
#' This function calculates the energy distance between two datasets. The energy distance
#' is a statistical distance between the distributions of random vectors, which is zero
#' if and only if the distributions are identical.
#'
#' Intuitively, the energy distance can be thought of as a measure of the "tension" between
#' two distributions. If we imagine the points in each distribution as particles with the
#' same charge, the energy distance is analogous to the potential energy of the system.
#' When the distributions are identical, the "tension" or "energy" is minimized.
#'
#' The energy distance has several desirable properties:
#' 1. It's simple to compute, requiring only pairwise distances between points.
#' 2. It's applicable to data of any dimension.
#' 3. It doesn't require density estimation, making it suitable for high-dimensional data.
#'
#' While not a direct proxy for relative entropy, the energy distance provides a robust
#' measure of distributional differences and can be used in similar contexts.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#'
#' @return A numeric value representing the energy distance between X and Y.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- energy.distance(X, Y)
#' print(result)
#'
#' @references
#' Szekely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of statistics based on distances.
#' Journal of statistical planning and inference, 143(8), 1249-1272.
#' @export
energy.distance <- function(X, Y) {
  n <- nrow(X)
  m <- ncol(X)

  # Ensure X and Y are matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # Compute pairwise distances
  XX <- as.matrix(dist(X))
  YY <- as.matrix(dist(Y))
  XY <- as.matrix(dist(rbind(X, Y)))[(n+1):(2*n), 1:n]

  # Compute energy distance
  ed <- 2 * mean(XY) - mean(XX) - mean(YY)

  return(ed)
}

#' Estimate Relative Entropy Using Entropy Difference
#'
#' @description
#' This function estimates the relative entropy between two datasets by computing
#' the difference between the joint entropy and the average of individual entropies.
#' The intuition behind this method comes from the relationship between mutual
#' information and entropy:
#'
#' I(X;Y) = H(X) + H(Y) - H(X,Y) = KL(P(X,Y) || P(X)P(Y))
#'
#' Where I(X;Y) is the mutual information, H() is entropy, and KL() is the Kullback-Leibler divergence.
#'
#' The algorithm works by:
#' 1. Discretizing the continuous data into bins.
#' 2. Estimating the entropy of X, Y, and their union separately.
#' 3. Computing the difference to approximate the relative entropy.
#'
#' This method provides a rough approximation of the relative entropy and can be
#' useful when dealing with high-dimensional data where direct density estimation is challenging.
#' However, it's sensitive to the choice of binning and may not be as accurate as some other methods.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#' @param num.bins The number of bins to use for discretization (default is 10).
#'
#' @return A numeric value representing the estimated relative entropy.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- entropy.difference(X, Y)
#' print(result)
#'
#' @seealso
#'   \code{\link[infotheo]{discretize}} for the discretization method,
#' @references
#' Cover, T. M., & Thomas, J. A. (2006). Elements of information theory. John Wiley & Sons.
#' @export
entropy.difference <- function(X, Y, num.bins = 10) {
    ## Check if package is available
    if (!requireNamespace("infotheo", quietly = TRUE)) {
        stop("Package 'infotheo' is required. Please install it with install.packages('infotheo')")
    }
    n <- ncol(X)
    ## Estimate entropy for X, Y, and their union
    entropy.X <- sum(sapply(1:n, function(i) {
        entropy(infotheo::discretize(X[,i], disc = "equalwidth", nbins = num.bins))
    }))
    entropy.Y <- sum(sapply(1:n, function(i) {
        entropy(infotheo::discretize(Y[,i], disc = "equalwidth", nbins = num.bins))
    }))
    entropy.union <- sum(sapply(1:n, function(i) {
        entropy(infotheo::discretize(c(X[,i], Y[,i]), disc = "equalwidth", nbins = num.bins))
    }))
    ## Estimate relative entropy
    relative.entropy <- entropy.union - 0.5 * (entropy.X + entropy.Y)
    return(relative.entropy)
}

#' Estimate Mutual Information Between Dataset Membership and Features
#'
#' @description
#' This function estimates the mutual information between the dataset membership (X or Y)
#' and the feature values. In the context of comparing two datasets, this mutual information
#' is equivalent to the Jensen-Shannon divergence, which is closely related to relative entropy.
#'
#' Intuitively, mutual information measures how much knowing the features tells us about
#' which dataset a point came from, and vice versa. If the datasets are very different,
#' knowing the features will give us a lot of information about which dataset a point
#' belongs to, resulting in high mutual information.
#'
#' The algorithm works by:
#' 1. Combining X and Y into a single dataset with labels.
#' 2. Discretizing the continuous data.
#' 3. Computing the mutual information between each feature and the dataset labels.
#' 4. Summing these mutual information values.
#'
#' This method is particularly useful when you want to understand which features contribute
#' most to the difference between the datasets, as you can look at the mutual information
#' for each feature separately.
#'
#' @param X A matrix or data frame representing the first dataset, where rows are observations
#'   and columns are features. Must have the same number of columns as Y.
#' @param Y A matrix or data frame representing the second dataset, where rows are observations
#'   and columns are features. Must have the same number of columns as X.
#' @param num.bins The number of bins to use for discretization (default is 10).
#' @return A numeric value representing the estimated mutual information in nats (natural units).
#'   Higher values indicate greater differences between the datasets. Returns 0 when the
#'   datasets have identical distributions.
#' @details
#'   The mutual information I(F;D) between features F and dataset label D is computed as:
#'   \deqn{I(F;D) = \sum_{f,d} p(f,d) \log\frac{p(f,d)}{p(f)p(d)}}{I(F;D) = sum p(f,d) * log(p(f,d)/(p(f)*p(d)))}
#'
#'   The function uses the infotheo package's discretization method (equal frequency binning
#'   by default) to handle continuous features. The total mutual information is the sum
#'   across all features, which assumes feature independence given the dataset label.
#' @note
#'   \itemize{
#'     \item Requires the 'infotheo' package to be installed
#'     \item The discretization step can lose information, especially for highly continuous data
#'     \item Results depend on the discretization method and number of bins used
#'     \item For high-dimensional data, consider using only the most informative features
#'     \item The assumption of summing MI across features may overestimate total information
#'       if features are correlated
#'   }
#' @examples
#' # Example 1: Datasets with different means
#' set.seed(123)
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- mutual.information(X, Y)
#' print(paste("MI between datasets:", round(result, 4), "nats"))
#'
#' # Example 2: Identical datasets (MI should be ~0)
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000), ncol = 2)
#' result <- mutual.information(X, Y)
#' print(paste("MI for identical distributions:", round(result, 4)))
#'
#' # Example 3: Feature-wise contribution
#' X <- matrix(rnorm(500), ncol = 5)
#' Y <- X
#' Y[,3] <- Y[,3] + 2  # Only change feature 3
#' # Compute MI for each feature separately
#' \dontrun{
#' feature_mi <- sapply(1:5, function(i) {
#'   mutual.information(X[,i,drop=FALSE], Y[,i,drop=FALSE])
#' })
#' barplot(feature_mi, names.arg = paste("Feature", 1:5),
#'         main = "Feature-wise MI Contribution")
#' }
#' @references
#' Cover, T. M., & Thomas, J. A. (2006). Elements of information theory. John Wiley & Sons.
#'
#' Lin, J. (1991). Divergence measures based on the Shannon entropy.
#' IEEE Transactions on Information Theory, 37(1), 145-151.
#' @seealso
#'   \code{\link[infotheo]{mutinformation}} for the underlying MI computation,
#'   \code{\link[infotheo]{discretize}} for the discretization method,
#'   \code{total.variation.distance} for an alternative dataset comparison metric
#' @export
mutual.information <- function(X, Y, num.bins = 10) {
  # Check if package is available
  if (!requireNamespace("infotheo", quietly = TRUE)) {
    stop("Package 'infotheo' is required. Please install it with install.packages('infotheo')")
  }
  n <- nrow(X)
  m <- ncol(X)
  # Create a dataset with labels
  data <- rbind(cbind(X, rep(0, n)), cbind(Y, rep(1, n)))
  # Discretize only the feature columns (not the label column)
  disc.features <- infotheo::discretize(data[, 1:m], disc = "equalwidth", nbins = num.bins)
  # Keep the label column as is (it's already discrete: 0 or 1)
  disc.data <- cbind(disc.features, data[, m+1])
  # Compute mutual information
  mi <- sapply(1:m, function(i) infotheo::mutinformation(disc.data[, i], disc.data[, m+1]))
  return(sum(mi))
}

# Classifier-based Approach
#' Estimate Divergence Using a Classifier-based Approach
#'
#' @description
#' This function estimates the divergence between two datasets using a classifier-based approach.
#' The intuition behind this method is that if two datasets are very different, it should be
#' easy to train a classifier to distinguish between them, resulting in a low classification error.
#'
#' The algorithm works by:
#' 1. Combining X and Y into a single dataset with labels.
#' 2. Training a logistic regression model with elastic net regularization to classify the points.
#' 3. Using the negative log-likelihood of the model as a proxy for divergence.
#'
#' This approach has several advantages:
#' 1. It can handle high-dimensional data well, especially with the elastic net regularization.
#' 2. It provides a flexible framework that can be adapted to different types of data by
#'    choosing appropriate classification algorithms.
#' 3. The regularization parameter alpha allows for tuning between L1 and L2 regularization,
#'    providing control over feature selection and multicollinearity handling.
#'
#' While not a direct measure of relative entropy, this method provides a robust way to
#' quantify the dissimilarity between datasets, especially in high-dimensional spaces.
#'
#' @param X A matrix or data frame representing the first dataset.
#' @param Y A matrix or data frame representing the second dataset.
#' @param alpha The elastic net mixing parameter, with 0 <= alpha <= 1. alpha=1 is the lasso penalty,
#'        and alpha=0 is the ridge penalty.
#'
#' @return A numeric value representing the estimated divergence.
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol = 2)
#' Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
#' result <- classifier.based.divergence(X, Y)
#' print(result)
#'
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear
#' models via coordinate descent. Journal of statistical software, 33(1), 1.
#'
#' @export
classifier.based.divergence <- function(X, Y, alpha = 0.5) {
  n <- nrow(X)
  m <- ncol(X)

  # Create a dataset with labels
  data <- rbind(X, Y)
  labels <- c(rep(0, n), rep(1, n))

  # Check if glmnet is available
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for this function. Please install it with install.packages('glmnet')")
  }
  
  # Fit a logistic regression model with elastic net regularization
  cv.fit <- glmnet::cv.glmnet(data, labels, family = "binomial", alpha = alpha)

  # Use the negative log-likelihood as a proxy for divergence
  divergence <- cv.fit$cvm[which.min(cv.fit$cvm)]

  return(divergence)
}

#' Calculate Jensen-Shannon Divergence Between Two Probability Mass Functions
#'
#' @description
#' Computes the Jensen-Shannon divergence (JSD) between two discrete probability
#' mass functions. The JSD is a symmetrized and smoothed version of the
#' Kullback-Leibler divergence. It is defined as:
#' JSD(P||Q) = 1/2 * D(P||M) + 1/2 * D(Q||M), where M = 1/2 * (P + Q)
#' and D denotes the Kullback-Leibler divergence.
#'
#' @param p Numeric vector representing the first probability mass function.
#'          Will be normalized to sum to 1. Must not contain negative values.
#' @param q Numeric vector representing the second probability mass function.
#'          Will be normalized to sum to 1. Must not contain negative values.
#'
#' @return A numeric value representing the Jensen-Shannon divergence in bits
#'         (using log base 2). The value is always non-negative and bounded
#'         above by 1.
#'
#' @examples
#' # Calculate JSD between two simple distributions
#' p <- c(0.4, 0.6)
#' q <- c(0.5, 0.5)
#' jensen.shannon.divergence(p, q)
#'
#' # Calculate JSD between distributions of different lengths
#' p <- c(0.3, 0.7)
#' q <- c(0.2, 0.3, 0.5)
#' jensen.shannon.divergence(p, q)  # p will be padded with 0 to match q's length
#'
#' # Calculate JSD between two categorical distributions
#' p <- c(0.2, 0.3, 0.5)
#' q <- c(0.1, 0.4, 0.5)
#' jensen.shannon.divergence(p, q)
#'
#' @details
#' The function implements the following behavior:
#' \itemize{
#'   \item If input vectors have different lengths, the shorter vector is padded
#'         with zeros to match the length of the longer vector
#'   \item Input vectors are automatically normalized to sum to 1
#'   \item Negative probabilities are not allowed
#'   \item Zero probabilities are handled appropriately by only computing
#'         the divergence over positive probability values
#' }
#'
#' When vectors of different lengths are provided, the function effectively treats
#' the shorter distribution as having zero probability for the additional categories
#' present in the longer distribution. This allows comparison of distributions
#' with different support sizes while preserving the mathematical properties of
#' the Jensen-Shannon divergence.
#'
#' @references
#' Lin, J. (1991). Divergence measures based on the Shannon entropy.
#' IEEE Transactions on Information Theory, 37(1), 145-151.
#'
#' @export
jensen.shannon.divergence <- function(p, q) {
  ## Pad shorter vector with zeros
  if (length(p) != length(q)) {
    max_length <- max(length(p), length(q))
    if (length(p) < max_length) {
      p <- c(p, rep(0, max_length - length(p)))
    } else {
      q <- c(q, rep(0, max_length - length(q)))
    }
  }

  ## Normalize vectors to sum to 1 (in case they weren't already)
  p <- p / sum(p)
  q <- q / sum(q)

  ## Input validation for negative values
  if (any(p < 0) || any(q < 0)) {
    stop("Probability vectors cannot contain negative values")
  }

  ## Calculate the midpoint distribution
  m <- (p + q) / 2

  ## Helper function for Kullback-Leibler divergence
  kl.divergence <- function(x, y) {
    ## Handle zero probabilities to avoid log(0)
    ## Only consider indices where x is positive
    indices <- x > 0
    sum(x[indices] * log2(x[indices] / y[indices]))
  }

  ## Calculate Jensen-Shannon divergence
  jsd <- (kl.divergence(p, m) + kl.divergence(q, m)) / 2
  return(jsd)
}

#' Calculate Kullback-Leibler Divergence Between Two Probability Mass Functions
#'
#' @description
#' Computes the Kullback-Leibler divergence (KL divergence) between two discrete
#' probability mass functions. The KL divergence is a measure of the information
#' lost when \code{Q} is used to approximate \code{P}. Note that KL divergence is not symmetric:
#' \eqn{\text{KL}(P||Q) \ne \text{KL}(Q||P)} in general.
#'
#' @param p Numeric vector representing the first probability mass function \code{P}.
#'          Will be normalized to sum to 1. Must not contain negative values.
#' @param q Numeric vector representing the second probability mass function \code{Q}.
#'          Will be normalized to sum to 1. Must not contain negative values.
#' @param zero.handling Character string specifying how to handle zeros in \code{Q}.
#'          Must be one of \code{"strict"} (error if \code{Q[i] = 0} where \code{P[i] > 0}),
#'          \code{"exclude"} (compute only where \code{Q > 0}), or
#'          \code{"pseudo"} (add small constant to \code{Q}). Defaults to \code{"pseudo"}.
#' @param epsilon Numeric value specifying the small constant to add when
#'          \code{zero.handling = "pseudo"}. Defaults to \code{1e-10}.
#'
#' @return A numeric value representing the Kullback-Leibler divergence in bits
#'         (using log base 2). The value is always non-negative, and equals
#'         infinity if there exists an \code{i} where \code{P[i] > 0} and \code{Q[i] = 0}
#'         (when \code{zero.handling = "strict"}).
#'
#' @examples
#' # Calculate KL divergence between two simple distributions
#' p <- c(0.4, 0.6)
#' q <- c(0.5, 0.5)
#' kullback.leibler.divergence(p, q)
#'
#' # Using different zero handling methods
#' p <- c(0.2, 0.8, 0)
#' q <- c(0.3, 0.6, 0.1)
#' kullback.leibler.divergence(p, q, zero.handling = "exclude")
#' kullback.leibler.divergence(p, q, zero.handling = "pseudo")
#'
#' # Calculate KL divergence between distributions of different lengths
#' p <- c(0.3, 0.7)
#' q <- c(0.2, 0.3, 0.5)
#' kullback.leibler.divergence(p, q)  # p will be padded with 0 to match q's length
#'
#' @details
#' The function implements the following behavior:
#' \itemize{
#'   \item If input vectors have different lengths, the shorter vector is padded
#'         with zeros to match the length of the longer vector
#'   \item Input vectors are automatically normalized to sum to 1
#'   \item Negative probabilities are not allowed
#' }
#' The implementation provides three ways to handle zeros in the \code{Q} distribution:
#' \itemize{
#'   \item \code{"strict"}: Throws an error if \code{Q[i] = 0} where \code{P[i] > 0}
#'   \item \code{"exclude"}: Computes divergence only over indices where \code{Q > 0}
#'   \item \code{"pseudo"}: Adds small constant \code{epsilon} to \code{Q} (and renormalizes)
#' }
#'
#' When vectors of different lengths are provided, the function treats the shorter
#' distribution as having zero probability for the additional categories present
#' in the longer distribution. This is particularly important when using
#' "strict" zero handling, as padding with zeros could trigger errors if P
#' is the shorter vector.
#'
#' @references
#' Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
#' The Annals of Mathematical Statistics, 22(1), 79-86.
#'
#' @export
kullback.leibler.divergence <- function(p, q, zero.handling = "pseudo", epsilon = 1e-10) {
    # Validate inputs
    if (!is.numeric(p) || !is.numeric(q)) {
        stop("Input vectors must be numeric")
    }
    if (any(is.infinite(p)) || any(is.infinite(q)) || any(is.na(p)) || any(is.na(q))) {
        stop("Input vectors contain infinite or NA values")
    }
    if (!is.numeric(epsilon) || epsilon <= 0) {
        stop("epsilon must be a positive number")
    }
    if (!zero.handling %in% c("strict", "exclude", "pseudo")) {
        stop('zero.handling must be one of "strict", "exclude", or "pseudo"')
    }

    # Handle different lengths
    if (length(p) != length(q)) {
        max_length <- max(length(p), length(q))
        if (length(p) < max_length) {
            warning("p is shorter than q; padding p with zeros")
            p <- c(p, rep(0, max_length - length(p)))
        } else {
            warning("q is shorter than p; padding q with zeros")
            q <- c(q, rep(0, max_length - length(q)))
        }

        if (zero.handling == "strict") {
            warning("Using 'strict' zero handling with padded zeros may raise errors")
        }
    }

    # Normalize vectors to sum to 1
    p <- p / sum(p)
    q <- q / sum(q)

    # Input validation for negative values
    if (any(p < 0) || any(q < 0)) {
        stop("Probability vectors cannot contain negative values")
    }

    # Handle zeros according to specified method
    if (zero.handling == "strict") {
        if (any(q == 0 & p > 0)) {
            stop("Q has zero probability where P is non-zero")
        }
        indices <- p > 0
    } else if (zero.handling == "exclude") {
        indices <- p > 0 & q > 0
    } else { # zero.handling == "pseudo"
        # Add epsilon to q and renormalize
        q.mod <- q + epsilon
        q.mod <- q.mod / sum(q.mod)
        q <- q.mod
        indices <- p > 0
    }

    # Calculate KL divergence
    sum(p[indices] * log2(p[indices] / q[indices]))
}

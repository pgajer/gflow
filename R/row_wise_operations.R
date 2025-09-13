#' Row-wise weighted mean (internal function)
#'
#' @description
#' Internal helper function for MAGELO that computes weighted means row-wise across
#' a matrix of y values using corresponding weight matrices. Used primarily for
#' computing predictions at grid points.
#'
#' @param nn.y A matrix of y values where rows correspond to grid points and columns
#'   to nearest neighbors
#' @param nn.w A matrix of weights with same dimensions as nn.y. Weights should sum
#'   to 1 for each row
#' @param max.K Integer vector of length nrow(nn.y) indicating the last column index
#'   with non-zero weight for each row (1-based indexing)
#'
#' @return Numeric vector of length nrow(nn.y) containing weighted means
#'
#' @keywords internal
#' @noRd
row.weighted.mean <- function(nn.y, nn.w, max.K) {
    ng <- nrow(nn.y)
    K <- ncol(nn.y)

    stopifnot(nrow(nn.w)==ng)
    stopifnot(ncol(nn.w)==K)

    Eyg <- numeric(ng)

    out <- .C(C_columnwise_wmean,
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             Eyg=as.double(Eyg))

    return(out$Eyg)
}

#' Bayesian bootstrap row-wise weighted means (internal function)
#'
#' @description
#' Internal helper function that generates Bayesian bootstrap samples of weighted means
#' by resampling weights according to a Dirichlet distribution. Used for uncertainty
#' quantification in MAGELO.
#'
#' @param nn.y Matrix of y values (rows: grid points, columns: nearest neighbors)
#' @param nn.w Matrix of weights with same dimensions as nn.y
#' @param max.K Integer vector indicating last non-zero weight column for each row
#' @param n.BB Number of Bayesian bootstrap iterations (default: 100)
#'
#' @return Matrix of dimensions nrow(nn.y) x n.BB containing bootstrap samples
#'   of weighted means
#'
#' @keywords internal
#' @noRd
row.weighted.mean.BB <- function(nn.y, nn.w, max.K, n.BB=100) {
    ng <- nrow(nn.y)
    K <- ncol(nn.y)

    stopifnot(nrow(nn.w)==ng)
    stopifnot(ncol(nn.w)==K)

    stopifnot(is.numeric(n.BB))
    stopifnot(n.BB>=1)

    BB.Eyg <- numeric(ng * n.BB)

    out <- .C(C_columnwise_wmean_BB,
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(n.BB),
             BB.Eyg=as.double(BB.Eyg))

    BB.Eyg <- matrix(out$BB.Eyg, nrow=ng, ncol=n.BB, byrow=FALSE)

    return(BB.Eyg)
}

#' Bayesian bootstrap credible intervals using quantiles (internal function)
#'
#' @description
#' Internal helper function that computes quantile-based credible intervals for
#' weighted means using Bayesian bootstrap. This is the primary method used by
#' MAGELO for constructing uncertainty bands.
#'
#' @param y.binary Logical; if TRUE, clips predictions to \code{[0,1]} interval
#' @param nn.y Matrix of y values (rows: grid points, columns: nearest neighbors)
#' @param nn.w Matrix of weights with same dimensions as nn.y
#' @param max.K Integer vector indicating last non-zero weight column for each row
#' @param n.BB Number of Bayesian bootstrap iterations (default: 100)
#' @param alpha Significance level for credible intervals (default: 0.05 for 95% CI)
#'
#' @return Matrix with 2 rows (lower and upper bounds) and ncol = nrow(nn.y),
#'   containing the credible interval bounds at each grid point
#'
#' @keywords internal
#' @noRd
row.weighted.mean.BB.qCrI <- function(y.binary, nn.y, nn.w, max.K, n.BB=100, alpha=0.05) {
    ng <- nrow(nn.y)
    K <- ncol(nn.y)

    stopifnot(nrow(nn.w)==ng)
    stopifnot(ncol(nn.w)==K)

    stopifnot(is.numeric(n.BB))
    stopifnot(n.BB>=1)

    Eyg.qCI <- numeric(ng * 2)

    out <- .C(C_columnwise_wmean_BB_qCrI,
             as.integer(y.binary),
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(n.BB),
             as.double(alpha),
             Eyg.qCI=as.double(Eyg.qCI))

    Eyg.qCI <- matrix(out$Eyg.qCI, nrow=2, ncol=ng, byrow=FALSE)

    return(Eyg.qCI)
}

#' Bayesian bootstrap credible intervals version 1 (internal function)
#'
#' @description
#' Internal helper function implementing an alternative method for computing
#' credible intervals based on bootstrap variance estimation. This version
#' uses a different approach than the quantile-based method.
#'
#' @param Eyg Numeric vector of pre-computed weighted means at grid points
#' @param nn.y Matrix of y values (rows: grid points, columns: nearest neighbors)
#' @param nn.w Matrix of weights with same dimensions as nn.y
#' @param max.K Integer vector indicating last non-zero weight column for each row
#' @param n.BB Number of Bayesian bootstrap iterations (default: 100)
#'
#' @return Numeric vector of length nrow(nn.y) containing half-widths of
#'   credible intervals (to be added/subtracted from Eyg)
#'
#' @keywords internal
#' @noRd
row.weighted.mean.BB.CI.v1 <- function(Eyg, nn.y, nn.w, max.K, n.BB=100) {
    ng <- nrow(nn.y)
    K <- ncol(nn.y)

    stopifnot(length(Eyg)==ng)
    stopifnot(nrow(nn.w)==ng)
    stopifnot(ncol(nn.w)==K)

    stopifnot(is.numeric(n.BB))
    stopifnot(n.BB>=1)

    Ey.CI <- numeric(ng)

    out <- .C(C_columnwise_wmean_BB_CrI_1,
             as.double(Eyg),
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(n.BB),
             Ey.CI=as.double(Ey.CI))

    return( out$Ey.CI )
}

#' Bayesian bootstrap credible intervals version 2 (internal function)
#'
#' @description
#' Internal helper function implementing the variance-based method for computing
#' credible intervals. This is the version actually used by MAGELO when
#' get.predictions.CrI = FALSE in degree 0 regression.
#'
#' @param Eyg Numeric vector of pre-computed weighted means at grid points
#' @param nn.y Matrix of y values (rows: grid points, columns: nearest neighbors)
#' @param nn.w Matrix of weights with same dimensions as nn.y
#' @param max.K Integer vector indicating last non-zero weight column for each row
#' @param n.BB Number of Bayesian bootstrap iterations (default: 100)
#'
#' @return Numeric vector of length nrow(nn.y) containing half-widths of
#'   credible intervals (to be added/subtracted from Eyg)
#'
#' @keywords internal
#' @noRd
row.weighted.mean.BB.CI.v2 <- function(Eyg, nn.y, nn.w, max.K, n.BB=100) {
    ng <- nrow(nn.y)
    K <- ncol(nn.y)

    stopifnot(length(Eyg)==ng)
    stopifnot(nrow(nn.w)==ng)
    stopifnot(ncol(nn.w)==K)

    stopifnot(is.numeric(n.BB))
    stopifnot(n.BB>=1)

    Ey.CI <- numeric(ng)

    out <- .C(C_columnwise_wmean_BB_CrI_2,
             as.double(Eyg),
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(n.BB),
             Ey.CI=as.double(Ey.CI))

    return( out$Ey.CI )
}


#' Compute weighted means for multiple response variables (internal function)
#'
#' @description
#' Internal helper function for MAGELO that computes row-wise weighted means for
#' multiple response variables simultaneously. This is used when applying MAGELO
#' to matrix-valued outcomes, where each column of Y represents a different
#' response variable (e.g., multiple permutations or bootstrap samples).
#'
#' @param Y Numeric matrix where rows correspond to observations (same length as
#'   the original x used to create nn matrices) and columns represent different
#'   response variables to be smoothed
#' @param nn.i Integer matrix of nearest neighbor indices (rows: grid points,
#'   columns: k nearest neighbors). Uses 1-based indexing
#' @param nn.w Numeric matrix of weights with same dimensions as nn.i. Weights
#'   should sum to 1 for each row
#' @param max.K Integer vector of length nrow(nn.i) indicating the last column
#'   index with non-zero weight for each row (1-based indexing)
#'
#' @return Matrix of dimensions nrow(nn.i) x ncol(Y) containing weighted mean
#'   predictions at each grid point for each response variable
#'
#' @details
#' This function efficiently computes weighted means for multiple response variables
#' in a single C call, which is more efficient than repeatedly calling
#' row.weighted.mean() for each column of Y. The function handles the conversion
#' to 0-based indexing for the C implementation and transposes matrices as needed.
#'
#' @keywords internal
#' @noRd
matrix.weighted.means <- function(Y, nn.i, nn.w, max.K) {
    nrY <- nrow(Y)
    ncY <- ncol(Y)

    Tnn.i <- t(nn.i-1)
    Tnn.w <- t(nn.w)
    nrTnn <- nrow(Tnn.i)
    ncTnn <- ncol(Tnn.i) # number of grid points

    Eyg <- numeric(ncTnn * ncY)

    out <- .C(C_matrix_wmeans,
             as.double(Y),
             as.integer(nrY),
             as.integer(ncY),
             as.integer(Tnn.i),
             as.double(Tnn.w),
             as.integer(nrTnn),
             as.integer(ncTnn),
             as.integer(max.K-1),
             Eyg=as.double(Eyg))

    matrix(out$Eyg, nrow=ncTnn, ncol=ncY, byrow=FALSE)
}

#' A generic interface for evaluating a kernel over a numeric vector
#'
#' @param x           A numeric vector.
#' @param kernel.str  The name (character string) of a kernel that will be evaluated over x.
#' @export
kernel.fn <- function(x, kernel.str="normal")
{
    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    ikernel <- pmatch(kernel.str, kernels)

    nx <- length(x)
    y <- numeric(nx)

    x[is.na(x)] <- 1

    out <- .C("C_kernel",
             as.integer(ikernel),
             as.double(x),
             as.integer(nx),
             y=as.double(y))

    return(out$y)
}

#' Triangular kernel
#'
#' The function equals to 1 - abs(x/bw) within (-bw, bw) and 0 otherwise
#'
#' @param x  A numeric vector.
#' @param bw A bandwidth numeric parameter.
#' @export
triangular.kernel <- function(x, bw=1)
{
    stopifnot(is.numeric(x))

    nx <- length(x)
    y <- numeric(nx)
    max.K <- -1 # index of the last element with value < bw. If no element has
               # value < bw, max.K in C will retain -1 value. In R it will be
               # changed to 0.

    out <- .C("C_triangular_kernel_with_stop",
             as.double(x),
             as.integer(nx),
             as.double(bw),
             max.K=as.integer(max.K),
             y=as.double(y))

    list(y=out$y,
         max.K=out$max.K + 1)
}


#' Epanechnikov kernel function
#'
#' The function equals to 1-t^2 for t in (-1,1) and 0 otherwise; t = abs(x/bw)
#'
#' @param x  A numeric vector.
#' @param bw A bandwidth numeric parameter.
#' @export
epanechnikov.kernel <- function(x, bw=1)
{
    stopifnot(is.numeric(x))

    nx <- length(x)
    y <- numeric(nx)
    max.K <- -1 # index of the last element with value < bw. If no element has
               # value < bw, max.K in C will retain -1 value. In R it will be
               # changed to 0.

    out <- .C("C_epanechnikov_kernel_with_stop",
             as.double(x),
             as.integer(nx),
             as.double(bw),
             max.K=as.integer(max.K),
             y=as.double(y))

    list(y=out$y,
         max.K=out$max.K + 1)
}


#'
#' Truncated exponential kernel
#'
#' The function equals to (1-t)*exp(-t) for t within (-1, 1) and 0 otherwise; t = abs(x/bw)
#'
#' @param x  A numeric vector.
#' @param bw A bandwidth numeric parameter.
#' @export
tr.exponential.kernel <- function(x, bw=1)
{
    stopifnot(is.numeric(x))

    nx <- length(x)
    y <- numeric(nx)
    max.K <- -1 # index of the last element with value < bw. If no element has
               # value < bw, max.K in C will retain -1 value. In R it will be
               # changed to 0.

    out <- .C("C_tr_exponential_kernel_with_stop",
             as.double(x),
             as.integer(nx),
             as.double(bw),
             max.K=as.integer(max.K),
             y=as.double(y))

    list(y=out$y,
         max.K=out$max.K + 1)
}


#' Generates kernel defined weights of the rows of an input numeric matrix
#'
#' @param x          A numeric matrix.
#' @param bws        A numeric vector of bandwidths.
#' @param kernel.str A character string indicating what kernel to apply to rows of x.
#'
#' ALERT: Values of max.K can be 0 if there is a row of x, where the weight of
#' the first elements is 0. Thus \code{max.K[i] = 0} indicates that there are no
#' non-zero elements in that row.
#'
#' @export
row.weighting <- function(x, bws, kernel.str="epanechnikov")
{
    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    ikernel <- pmatch(kernel.str, kernels)

    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))

    nr <- nrow(x)
    nc <- ncol(x)

    stopifnot(is.finite(bws))
    stopifnot(length(bws)==nr)

    w <- numeric(nr * nc)
    max.K <- numeric(nr)

    out <- .C("C_columnwise_weighting",
             as.double(t(x)),
             as.integer(nc),
             as.integer(nr),
             as.integer(ikernel),
             as.double(bws),
             max.K=as.integer(max.K),
             w=as.double(w))

    w <- t(matrix(out$w, nrow=nc, ncol=nr, byrow=FALSE))

    max.K <- out$max.K + 1 # max.K is the __index__ of the last non-zero weight in
                          # each row of x; since it is calculated in C, it is
                          # 0-based, so we need to add 1 for 1-based indexing of R
    list(nn.w=w,
         max.K=max.K)
}


#' Normalizes a distance matrix
#'
#' normalize.dist(d, min.K, bw) divides the rows of the input distance matrix,
#' d, by 'bw' if the distance of the 'min.K'-th element of the given row, i, is
#' <bw, if it is not, then all elements are divided by the distance of the
#' min.K-st element. The radius, r, is bw in the first case and \code{d[i,minK]} in the
#' second.
#'
#' @param d      A numeric matrix.
#' @param min.K  The minimum __number_ of elements of each row of d after normalization with normalized distance < 1.
#' @param bw     A normalization non-negative constant.
#'
#' @return A list with two components: nd (normalized distances), r (radii).
#'
#' ALERT: This routine makes sense only in the situation when ncol(d) > minK. Thus,
#' the user has to make sure this condition is satisfied before calling
#' normalize_dist().
#'
#' @export
normalize.dist <- function(d, min.K, bw)
{
    nr <- nrow(d)
    nc <- ncol(d)

    ##stopifnot( min.K+1 < nc )
    if( nc - 1 < min.K  )
    {
        stop(paste0("min.K+1 < nc:  min.K=",min.K," nc=",nc))
    }

    nd <- matrix(0, nrow=nc, ncol=nr)
    r <- numeric(nr)

    out <- .C("C_normalize_dist", ## this C function normalizes the columns of the input matrix
             as.double(t(d)),
             as.integer(nc),   ## number of rows of t(d)
             as.integer(nr),   ## number of columns of t(d)
             as.integer(min.K),
             as.double(bw),
             nd=as.double(nd),
             r=as.double(r))

    nd <- t(matrix(out$nd, nrow=nc, ncol=nr, byrow=FALSE))

    list(nd=nd, r=out$r)
}

#' Gets bandwidths defined as radii of linear model's disks of support
#'
#' Finds a bandwidth for each column of d, such that \code{bws[i] = bw}
#' if \code{d[iminK + ir] < bw} and \code{bws[i] = d[iminK + 1 + ir] = d[minK + ir]}, otherwise.
#'
#' @param d      A numeric matrix nn.d of distances to NNs.
#' @param min.K  The minimum number of elements of each row with weights > 0.
#' @param bw     A bandwidth.
#'
#' @return A vector of bandwidths.
get.bws <- function(d, min.K, bw)
{
    nr <- nrow(d)
    nc <- ncol(d)

    stopifnot( min.K <= nc )

    bws <- numeric(nr)

    out <- .C("C_get_bws",
             as.double(t(d)),
             as.integer(nc),    # number of rows of t(d)
             as.integer(nr),    # number of columns of t(d)
             as.integer(min.K), # minimal number of elements with weights > 0
             as.double(bw),
             bws=as.double(bws))

    out$bws
}


#' Gets bandwidths defined as radia of linear model's disks of support
#'
#' The difference between this function and get.bws() is that in the former,
#' min.K is a constant, which is replaced by a vector minK in here. This
#' function is not meant to be accessible by user and serves as an R interface
#' to get_bws_with_minK_a() for unit testing this C function.
#'
#' @param d      A numeric matrix.
#'
#' @param minK   A vector such that minK\[i\] is the required minimum number of
#'               elements in the row of the i-th column of d with weights > 0.
#'
#' @param bw     A bandwidth.
#'
#' @return A vector of bandwidths.
#' @export
get.bws.with.minK <- function(d, minK, bw)
{
    nr <- nrow(d)
    nc <- ncol(d)

    stopifnot(length(minK)==nr)

    bws <- numeric(nr)

    out <- .C("C_get_bws_with_minK_a",
             as.double(t(d)),
             as.integer(nc),    # number of rows of t(d)
             as.integer(nr),    # number of columns of t(d)
             as.integer(minK),  # minimal number of elements with weights > 0
             as.double(bw),
             bws=as.double(bws))

    out$bws
}

#' Row-wise total sum normalization
#'
#' @param x  A matrix.
#'
#' @return A matrix obtained from x by dividing the rows of x by the sum of the
#'     given row's elements, if the sum is not 0.
#' @export
row.TS.norm <- function(x)
{
    nr <- nrow(x)
    nc <- ncol(x)

    nx <- numeric(nr * nc)

    out <- .C("C_columnwise_TS_norm",
             as.double(t(x)),
             as.integer(nc),
             as.integer(nr),
             nx=as.double(nx))

    nx <- matrix(out$nx, nrow=nr, ncol=nc, byrow=TRUE)

    return(nx)
}

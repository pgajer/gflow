#' Arcsine square-root (angular) transform for proportions
#'
#' @description
#' Applies the arcsine square-root (angular) transformation
#' \deqn{y = \arcsin(\sqrt{x})}
#' to proportion data \eqn{x \in [0,1]}. This transformation is often used as a
#' variance-stabilizing map for proportion-like data and expands the scale near
#' 0 and 1.
#'
#' @param x Numeric vector, matrix, or array of proportions in \eqn{[0,1]}.
#' @param clip Logical; if TRUE (default), values are clipped to \eqn{[0,1]}
#'   before transforming.
#' @param eps Nonnegative numeric; if \code{clip = TRUE}, values in
#'   \eqn{[-eps, 1+eps]} are clipped to \eqn{[0,1]}, otherwise an error is thrown.
#' @return An object of the same shape as \code{x} containing transformed values
#'   in radians.
#'
#' @examples
#' x <- c(0, 0.1, 0.5, 0.9, 1)
#' arcsin.sqrt.transform(x)
#'
#' X <- matrix(runif(12), nrow = 3)
#' arcsin.sqrt.transform(X)
#'
#' @export
arcsin.sqrt.transform <- function(x, clip = TRUE, eps = 0) {
    ## Input validation
    if (!is.numeric(x)) {
        stop("x must be numeric.")
    }
    if (!is.logical(clip) || length(clip) != 1L) {
        stop("clip must be a logical scalar.")
    }
    eps <- as.numeric(eps)
    if (!is.finite(eps) || eps < 0) {
        stop("eps must be a nonnegative finite number.")
    }

    ## Handle missing values
    y <- x

    if (clip) {
        ## Allow small numerical spillover controlled by eps
        bad.low <- which(is.finite(y) & (y < -eps))
        bad.high <- which(is.finite(y) & (y > 1 + eps))
        if (length(bad.low) > 0L || length(bad.high) > 0L) {
            stop("x contains values outside [0,1] beyond the allowed eps tolerance.")
        }
        y[is.finite(y) & y < 0] <- 0
        y[is.finite(y) & y > 1] <- 1
    } else {
        if (any(is.finite(y) & (y < 0 | y > 1))) {
            stop("x must be in [0,1] when clip = FALSE.")
        }
    }

    asin(sqrt(y))
}

#' Centered Log-Ratio (CLR) Transform for Compositional Data
#'
#' @description
#' Applies the centered log-ratio (CLR) transform to a numeric matrix of compositional
#' measurements (e.g., relative abundances). For each row, the CLR transform is:
#' \deqn{\mathrm{clr}(x_j) = \log(x_j + \epsilon) - \frac{1}{p}\sum_{k=1}^p \log(x_k + \epsilon),}
#' where \eqn{p} is the number of columns (features) and \code{epsilon} is a pseudocount.
#'
#' This is a convenience wrapper around \code{\link{clr.transform.matrix}} that supports
#' optional subsetting of rows and columns.
#'
#' @param X A numeric matrix or a data.frame coercible to a numeric matrix.
#'   Rows correspond to samples/vertices; columns correspond to features.
#' @param rows Optional integer vector of row indices (1-based) selecting which rows to transform.
#'   If NULL (default), all rows are used.
#' @param cols Optional integer vector of column indices (1-based) selecting which columns to transform.
#'   If NULL (default), all columns are used.
#' @param pseudo.count Positive numeric scalar. Pseudocount added to \code{X} prior to log transform.
#'   Default is \code{1e-6}.
#' @param na.rm Logical. If TRUE (default), rows containing any non-finite values after subsetting
#'   are removed prior to computing CLR. If FALSE, an error is thrown when non-finite values are present.
#'
#' @return A numeric matrix of CLR values with \code{length(rows)} rows (after any \code{na.rm}
#'   filtering) and \code{length(cols)} columns. Row and column names are preserved when present.
#'
#' @examples
#' \dontrun{
#' ## CLR over all rows/columns
#' clr.X <- clr.transform(phi.zmb, pseudo.count = 1e-6)
#'
#' ## CLR over a subset of vertices (rows)
#' clr.disk <- clr.transform(phi.zmb, rows = M1.disk.res$vertices)
#' }
#'
#' @seealso \code{\link{clr.transform.matrix}}
#'
#' @export
clr.transform <- function(X,
                          rows = NULL,
                          cols = NULL,
                          pseudo.count = 1e-6,
                          na.rm = TRUE) {

    if (is.data.frame(X)) {
        X <- as.matrix(X)
    }
    if (!is.matrix(X) || !is.numeric(X)) {
        stop("X must be a numeric matrix or a data.frame coercible to a numeric matrix.")
    }

    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric value.")
    }

    nr <- nrow(X)
    nc <- ncol(X)

    if (is.null(rows)) {
        rows <- seq_len(nr)
    } else {
        rows <- as.integer(rows)
        rows <- rows[!is.na(rows)]
        rows <- rows[rows >= 1L & rows <= nr]
        if (length(rows) == 0L) stop("rows is empty after filtering to valid indices.")
    }

    if (is.null(cols)) {
        cols <- seq_len(nc)
    } else {
        cols <- as.integer(cols)
        cols <- cols[!is.na(cols)]
        cols <- cols[cols >= 1L & cols <= nc]
        if (length(cols) == 0L) stop("cols is empty after filtering to valid indices.")
    }

    X.sub <- X[rows, cols, drop = FALSE]

    if (isTRUE(na.rm)) {
        ok <- apply(X.sub, 1, function(z) all(is.finite(z)))
        X.sub <- X.sub[ok, , drop = FALSE]
    } else {
        if (any(!is.finite(X.sub))) stop("Non-finite values detected in X after subsetting; set na.rm=TRUE to drop rows.")
    }

    clr.transform.matrix(X.sub, pseudo.count = pseudo.count)
}


#' Centered Log-Ratio (CLR) Transform for a Numeric Matrix
#'
#' @description
#' Computes the centered log-ratio (CLR) transform of a numeric matrix \code{X}.
#' CLR is computed row-wise:
#' \deqn{\mathrm{clr}(x_j) = \log(x_j + \epsilon) - \frac{1}{p}\sum_{k=1}^p \log(x_k + \epsilon),}
#' where \eqn{p} is the number of columns and \code{epsilon} is a pseudocount.
#'
#' This function performs the core CLR computation and expects a numeric matrix input.
#'
#' @param X A numeric matrix with rows representing samples/vertices and columns representing features.
#' @param pseudo.count Positive numeric scalar. Pseudocount added to \code{X} prior to log transform.
#'   Default is \code{1e-6}.
#' @param check Logical. If TRUE (default), validates that \code{X} is numeric and finite. If FALSE,
#'   assumes inputs are valid (useful for internal performance).
#'
#' @return A numeric matrix of the same dimensions as \code{X}, containing CLR values.
#' Row and column names are preserved when present.
#'
#' @examples
#' \dontrun{
#' X.disk <- phi.zmb[M1.disk.res$vertices, , drop = FALSE]
#' clr.disk <- clr.transform.matrix(X.disk, pseudo.count = 1e-6)
#' }
#'
#' @export
clr.transform.matrix <- function(X,
                                 pseudo.count = 1e-6,
                                 check = TRUE) {

    if (isTRUE(check)) {
        if (!is.matrix(X) || !is.numeric(X)) {
            stop("X must be a numeric matrix.")
        }
        if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
            stop("pseudo.count must be a single positive numeric value.")
        }
        if (any(!is.finite(X))) {
            stop("X contains non-finite values; filter or impute before calling clr.transform.matrix().")
        }
    }

    log.X <- log(X + pseudo.count)
    row.center <- rowMeans(log.X)
    out <- sweep(log.X, 1, row.center, FUN = "-")

    ## Preserve dimnames
    dimnames(out) <- dimnames(X)

    out
}

#' Construct an SSRHE-Style Local Hessian Energy Operator
#'
#' Constructs the discrete local Hessian operator used by the semi-supervised
#' regularization framework of Kim, Steinke, and Hein (2009). The constructor
#' is intended as a package-facing, auditable R/C++ port of the reference
#' behavior: it returns the sparse row operator \eqn{A}, the quadratic energy
#' matrix \eqn{B=A^\top A}, and optionally the supplemental stabilizer
#' \eqn{B_S}.
#'
#' @param X Numeric matrix with one observation per row. Distances for the
#'   default neighborhood search are Euclidean distances in these coordinates.
#' @param k Integer number of nearest neighbors per point, including the point
#'   itself. This follows the SSRHE Matlab convention.
#' @param tangent.dim Integer tangent dimension used for the local PCA chart.
#'   Required when \code{tangent.dim.rule = "fixed"}. If \code{NULL} and
#'   \code{tangent.dim.rule = "eigen.cumulative"}, the dimension is selected
#'   locally from the PCA variance ratio.
#' @param nn.index Optional integer matrix of neighbor indices, with
#'   \code{nrow(X)} rows and \code{k} columns. Each row must contain its center
#'   vertex. If \code{NULL}, Euclidean k-NN including self is computed in C++.
#' @param tangent.dim.rule Either \code{"fixed"} or \code{"eigen.cumulative"}.
#' @param eigen.tolerance Cumulative local PCA variance threshold used when
#'   \code{tangent.dim.rule = "eigen.cumulative"} and \code{tangent.dim = NULL}.
#' @param stabilizer Logical. If \code{TRUE}, also construct the supplemental
#'   stabilizer matrix described in the SSRHE supplement.
#' @param pinv.tol Nonnegative tolerance multiplier for local pseudoinverses.
#' @param return.A Logical. If \code{TRUE}, return \eqn{A} as a sparse matrix.
#' @param return.B Logical. If \code{TRUE}, return \eqn{B=A^\top A}.
#' @param return.BS Logical. If \code{TRUE} and \code{stabilizer = TRUE}, return
#'   the supplemental stabilizer \eqn{B_S}.
#' @param return.sparse Logical. If \code{TRUE}, attach \pkg{Matrix} sparse
#'   matrices in addition to raw triplets.
#' @param verbose Logical. If \code{TRUE}, print a short native construction
#'   message.
#'
#' @details
#' For each point \eqn{x_i}, the operator constructs a local PCA chart from the
#' \code{k} nearest neighbors, projects the neighborhood into
#' \eqn{\mathbb R^m}, and centers the local coordinates at \eqn{x_i}. A local
#' quadratic model is then fitted with the intercept fixed at \eqn{f_i}:
#' \deqn{
#'   f(x_j) - f(x_i)
#'   \approx
#'   \sum_{a \le b} h_{ab} z_{ja} z_{jb}
#'   +
#'   \sum_a g_a z_{ja}.
#' }
#' The local linear map from function values to the quadratic coefficients is
#' the Matlab-compatible
#' \deqn{
#'   \mathrm{RegMat} = X^+ - X^+ \mathrm{IndMat},
#' }
#' where \eqn{X} is the reduced quadratic-plus-linear design and
#' \eqn{\mathrm{IndMat}} copies the center value into all local rows. Diagonal
#' Hessian components are scaled by \eqn{\sqrt{2}} and off-diagonal components
#' by \eqn{1}, so that the row energy matches the reference Hessian energy.
#'
#' The returned \eqn{A} is the stacked sparse matrix of these local Hessian
#' coefficient rows. The \eqn{\ell_2} SSRHE penalty is
#' \deqn{
#'   f^\top B f = \|Af\|_2^2,\quad B=A^\top A.
#' }
#' Exposing \eqn{A} is useful for auditability and for future
#' \eqn{\ell_1}-style variants based on \eqn{\|Af\|_1}.
#'
#' @return A list of class \code{"ssrhe.hessian.operator"} containing:
#'   \itemize{
#'     \item \code{A}: sparse local Hessian row operator, if requested;
#'     \item \code{B}: sparse quadratic energy matrix \eqn{A^\top A}, if requested;
#'     \item \code{BS}: optional supplemental stabilizer matrix;
#'     \item \code{A.triplet}, \code{BS.triplet}: raw sparse triplets;
#'     \item \code{nn.index}: neighbor index matrix;
#'     \item \code{row.table}: one row per Hessian component row of \eqn{A};
#'     \item \code{diagnostics}: local design rank, condition, and PCA summaries;
#'     \item \code{parity}: notes documenting reference-behavior conventions.
#'   }
#'
#' @references
#' Kim, K. I., Steinke, F., and Hein, M. (2009). Semi-supervised regression
#' using Hessian energy with an application to semi-supervised dimensionality
#' reduction. \emph{Neural Computation}.
#'
#' @export
ssrhe.hessian.operator <- function(
    X,
    k,
    tangent.dim,
    nn.index = NULL,
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    stabilizer = FALSE,
    pinv.tol = sqrt(.Machine$double.eps),
    return.A = TRUE,
    return.B = TRUE,
    return.BS = stabilizer,
    return.sparse = TRUE,
    verbose = FALSE) {

    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    k <- .validate.ssrhe.positive.integer(k, "k")
    if (k < 2L || k > n) {
        stop("k must be between 2 and nrow(X).", call. = FALSE)
    }

    tangent.dim.rule <- match.arg(tangent.dim.rule)
    if (missing(tangent.dim) || is.null(tangent.dim)) {
        if (identical(tangent.dim.rule, "fixed")) {
            stop("tangent.dim is required when tangent.dim.rule = 'fixed'.",
                 call. = FALSE)
        }
        tangent.dim.cpp <- 0L
    } else {
        tangent.dim.cpp <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        if (tangent.dim.cpp > min(k, ncol(X))) {
            stop("tangent.dim must be no larger than min(k, ncol(X)).",
                 call. = FALSE)
        }
    }

    eigen.tolerance <- .validate.ssrhe.numeric.scalar(eigen.tolerance, "eigen.tolerance")
    if (eigen.tolerance <= 0 || eigen.tolerance > 1) {
        stop("eigen.tolerance must be in (0, 1].", call. = FALSE)
    }
    pinv.tol <- .validate.ssrhe.numeric.scalar(pinv.tol, "pinv.tol")
    if (pinv.tol < 0) {
        stop("pinv.tol must be nonnegative.", call. = FALSE)
    }
    stabilizer <- isTRUE(stabilizer)
    return.A <- isTRUE(return.A)
    return.B <- isTRUE(return.B)
    return.BS <- isTRUE(return.BS)
    return.sparse <- isTRUE(return.sparse)
    verbose <- isTRUE(verbose)

    if ((return.A || return.B || return.BS || return.sparse) &&
        !requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required to attach sparse SSRHE operators.",
             call. = FALSE)
    }
    if (!return.sparse && (return.A || return.B || return.BS)) {
        stop("return.sparse must be TRUE when returning A, B, or BS matrices.",
             call. = FALSE)
    }

    nn.index.cpp <- .validate.ssrhe.nn.index(nn.index, n, k)

    raw <- .Call(
        "S_ssrhe_hessian_operator",
        X,
        as.integer(k),
        as.integer(tangent.dim.cpp),
        nn.index.cpp,
        tangent.dim.rule,
        as.double(eigen.tolerance),
        as.logical(stabilizer),
        as.double(pinv.tol),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    raw$row.table <- as.data.frame(raw$row.table, stringsAsFactors = FALSE)
    raw$diagnostics <- as.data.frame(raw$diagnostics, stringsAsFactors = FALSE)

    A <- .ssrhe.triplet.to.sparse(raw$A.triplet)
    B <- Matrix::crossprod(A)
    BS <- NULL
    if (stabilizer && !is.null(raw$BS.triplet)) {
        BS <- .ssrhe.triplet.to.sparse(raw$BS.triplet)
        BS <- (BS + Matrix::t(BS)) / 2
    }

    out <- list(
        A = if (return.A) A else NULL,
        B = if (return.B) B else NULL,
        BS = if (return.BS && stabilizer) BS else NULL,
        A.triplet = raw$A.triplet,
        BS.triplet = raw$BS.triplet,
        nn.index = raw$nn.index,
        row.table = raw$row.table,
        diagnostics = raw$diagnostics,
        parity = raw$parity,
        parameters = raw$parameters
    )

    class(out) <- c("ssrhe.hessian.operator", "list")
    out
}

.validate.ssrhe.X <- function(X) {
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    storage.mode(X) <- "double"
    if (!is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L) {
        stop("X must be a numeric matrix with at least two rows and one column.",
             call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("X must contain only finite values.", call. = FALSE)
    }
    X
}

.validate.ssrhe.positive.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 ||
        abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a positive integer scalar.", name), call. = FALSE)
    }
    as.integer(round(x))
}

.validate.ssrhe.numeric.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) {
        stop(sprintf("%s must be a finite numeric scalar.", name), call. = FALSE)
    }
    as.double(x)
}

.validate.ssrhe.nn.index <- function(nn.index, n, k) {
    if (is.null(nn.index)) {
        return(NULL)
    }
    if (!is.matrix(nn.index)) {
        nn.index <- as.matrix(nn.index)
    }
    storage.mode(nn.index) <- "integer"
    if (!identical(dim(nn.index), c(n, k))) {
        stop("nn.index must be an nrow(X) by k integer matrix.", call. = FALSE)
    }
    if (any(is.na(nn.index)) || any(nn.index < 1L) || any(nn.index > n)) {
        stop("nn.index must contain integer indices in 1:nrow(X).", call. = FALSE)
    }
    has.center <- vapply(seq_len(n), function(i) any(nn.index[i, ] == i), logical(1))
    if (!all(has.center)) {
        stop("Each row of nn.index must contain its center vertex.", call. = FALSE)
    }
    nn.index
}

.ssrhe.triplet.to.sparse <- function(triplet) {
    Matrix::sparseMatrix(
        i = triplet$i,
        j = triplet$j,
        x = triplet$x,
        dims = triplet$dim
    )
}

#' @export
print.ssrhe.hessian.operator <- function(x, ...) {
    n <- if (!is.null(x$B)) ncol(x$B) else x$A.triplet$dim[[2]]
    nrow.A <- x$A.triplet$dim[[1]]
    cat("SSRHE Hessian operator\n")
    cat("  vertices:", n, "\n")
    cat("  operator rows:", nrow.A, "\n")
    cat("  k:", x$parameters$k, "\n")
    cat("  tangent.dim.rule:", x$parameters$tangent.dim.rule, "\n")
    cat("  stabilizer:", isTRUE(x$parameters$stabilizer), "\n")
    invisible(x)
}

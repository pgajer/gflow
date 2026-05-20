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

#' Fit SSRHE-Style Hessian-Energy Regression
#'
#' Fits the \eqn{\ell_2} Hessian-energy regularized estimator associated with
#' \code{\link{ssrhe.hessian.operator}}. This is a direct SSRHE-style
#' comparator for Hessian smoothing on point clouds; it is not a replacement
#' for \code{\link{fit.rdgraph.regression}} and is distinct from
#' \eqn{\ell_1}-adaptive graph trend filtering.
#'
#' @inheritParams ssrhe.hessian.operator
#' @param y Numeric response vector of length \code{nrow(X)} or numeric matrix
#'   with \code{nrow(X)} rows. Matrix columns are fit as separate responses.
#'   \code{NA} values are allowed and are treated as unobserved by setting the
#'   corresponding observation weight to zero.
#' @param lambda1 Nonnegative Hessian-energy penalty multiplier for
#'   \eqn{f^\top B f = \|Af\|_2^2}.
#' @param lambda2 Nonnegative supplemental-stabilizer multiplier for
#'   \eqn{f^\top B_S f}. If positive, \code{stabilizer} must be \code{TRUE}.
#' @param weights Optional nonnegative observation weights. May be \code{NULL},
#'   a vector of length \code{nrow(X)}, or a matrix with the same dimensions as
#'   \code{y}.
#' @param ridge Nonnegative diagonal ridge added to the linear system for
#'   numerical stabilization.
#'
#' @details
#' For each response column, this function solves
#' \deqn{
#'   (W + \lambda_1 B + \lambda_2 B_S + \epsilon I)\hat f = Wy,
#' }
#' where \eqn{W} is the diagonal matrix of observation weights and
#' \eqn{\epsilon} is \code{ridge}. Fully observed data with
#' \code{lambda1 = lambda2 = ridge = 0} therefore reproduce the observed
#' response exactly. Missing responses are excluded from the data-fit term by
#' setting their weights to zero. The Matlab SSRHE semi-supervised convention is
#' represented by a \code{0}/\code{1} labeled-indicator \code{weights} vector, or
#' equivalently by setting unlabeled responses to \code{NA}: labeled vertices
#' contribute unit diagonal data-fit terms and unlabeled vertices contribute no
#' data-fit term.
#'
#' The lower-level operator is matched to the Kim--Steinke--Hein SSRHE Matlab
#' construction: self-including kNN neighborhoods, local PCA charts, a
#' fixed-intercept local quadratic fit, doubled diagonal Hessian components,
#' and optional supplemental stabilizer. The package exposes both \eqn{A} and
#' \eqn{B=A^\top A}; the fitted \eqn{\ell_2} estimator uses \eqn{B}.
#'
#' @return A list of class \code{"ssrhe.hessian.fit"} containing fitted values,
#'   residuals, input response, weights, lambda parameters, objective/energy
#'   diagnostics, the reused \code{operator}, solver metadata, and the call.
#'
#' @references
#' Kim, K. I., Steinke, F., and Hein, M. (2009). Semi-supervised regression
#' using Hessian energy with an application to semi-supervised dimensionality
#' reduction. \emph{Neural Computation}.
#'
#' @export
fit.ssrhe.hessian.regression <- function(
    X,
    y,
    k,
    tangent.dim,
    lambda1,
    lambda2 = 0,
    weights = NULL,
    nn.index = NULL,
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    stabilizer = lambda2 > 0,
    pinv.tol = sqrt(.Machine$double.eps),
    ridge = 0,
    return.A = TRUE,
    verbose = FALSE) {

    X <- .validate.ssrhe.X(X)
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    stabilizer <- isTRUE(stabilizer)
    if (lambda2 > 0 && !stabilizer) {
        stop("lambda2 > 0 requires stabilizer = TRUE.", call. = FALSE)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        stabilizer = stabilizer,
        pinv.tol = pinv.tol,
        return.A = return.A,
        return.B = TRUE,
        return.BS = stabilizer,
        return.sparse = TRUE,
        verbose = verbose
    )

    out <- .fit.ssrhe.hessian.from.operator(
        operator = operator,
        y = y,
        lambda1 = lambda1,
        lambda2 = lambda2,
        weights = weights,
        ridge = ridge,
        verbose = verbose
    )
    out$X <- X
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.fit", "list")
    out
}

#' Refit SSRHE-Style Hessian-Energy Regression
#'
#' Reuses the operator from \code{\link{fit.ssrhe.hessian.regression}} to fit
#' new responses or new fixed penalty weights without rebuilding local PCA
#' neighborhoods or Hessian-energy matrices.
#'
#' @param fitted.model A \code{"ssrhe.hessian.fit"} object.
#' @param y.new New numeric response vector or matrix with one row per vertex.
#'   If \code{NULL}, the original response is reused.
#' @inheritParams fit.ssrhe.hessian.regression
#'
#' @return A list of class \code{"ssrhe.hessian.refit"}.
#' @export
refit.ssrhe.hessian.regression <- function(fitted.model,
                                           y.new = NULL,
                                           lambda1 = fitted.model$lambda$lambda1,
                                           lambda2 = fitted.model$lambda$lambda2,
                                           weights = NULL,
                                           ridge = fitted.model$lambda$ridge,
                                           verbose = FALSE) {
    if (!inherits(fitted.model, "ssrhe.hessian.fit")) {
        stop("fitted.model must be a 'ssrhe.hessian.fit' object.", call. = FALSE)
    }
    if (is.null(fitted.model$operator) ||
        !inherits(fitted.model$operator, "ssrhe.hessian.operator")) {
        stop("fitted.model does not contain a reusable SSRHE operator.",
             call. = FALSE)
    }
    if (is.null(y.new)) {
        y.new <- fitted.model$y
    }
    if (is.null(weights)) {
        weights <- fitted.model$weights
    }
    out <- .fit.ssrhe.hessian.from.operator(
        operator = fitted.model$operator,
        y = y.new,
        lambda1 = lambda1,
        lambda2 = lambda2,
        weights = weights,
        ridge = ridge,
        verbose = verbose
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.refit", "list")
    out
}

.fit.ssrhe.hessian.from.operator <- function(operator,
                                             y,
                                             lambda1,
                                             lambda2,
                                             weights,
                                             ridge,
                                             verbose = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE regression solves.",
             call. = FALSE)
    }
    if (is.null(operator$B)) {
        stop("operator must contain B. Rebuild with return.B = TRUE.",
             call. = FALSE)
    }
    n <- ncol(operator$B)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    Y.raw <- y.info$Y
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    Y.clean <- Y.raw
    Y.clean[!is.finite(Y.clean)] <- 0
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    if (lambda2 > 0 && is.null(operator$BS)) {
        stop("lambda2 > 0 requires an operator with BS. Rebuild with stabilizer = TRUE.",
             call. = FALSE)
    }

    n.responses <- ncol(Y.clean)
    Y.hat <- matrix(NA_real_, nrow = n, ncol = n.responses)
    residuals <- matrix(NA_real_, nrow = n, ncol = n.responses)
    data.loss <- hessian.energy <- stabilizer.energy <- objective <- numeric(n.responses)
    solve.method <- character(n.responses)

    common.weights <- n.responses == 1L || all(vapply(seq_len(n.responses), function(j) {
        identical(W[, j], W[, 1L])
    }, logical(1)))

    if (common.weights) {
        solved <- .solve.ssrhe.hessian.system(
            B = operator$B,
            BS = operator$BS,
            weights = W[, 1L],
            Y = Y.clean,
            lambda1 = lambda1,
            lambda2 = lambda2,
            ridge = ridge,
            verbose = verbose
        )
        Y.hat[,] <- solved$fitted
        solve.method[] <- solved$method
    } else {
        for (j in seq_len(n.responses)) {
            solved <- .solve.ssrhe.hessian.system(
                B = operator$B,
                BS = operator$BS,
                weights = W[, j],
                Y = Y.clean[, j, drop = FALSE],
                lambda1 = lambda1,
                lambda2 = lambda2,
                ridge = ridge,
                verbose = verbose
            )
            Y.hat[, j] <- solved$fitted[, 1L]
            solve.method[j] <- solved$method
        }
    }

    for (j in seq_len(n.responses)) {
        observed <- y.info$observed[, j] & W[, j] > 0
        residuals[observed, j] <- Y.raw[observed, j] - Y.hat[observed, j]
        residuals[!observed, j] <- NA_real_
        data.loss[j] <- 0.5 * sum(W[, j] * (Y.clean[, j] - Y.hat[, j])^2)
        hessian.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$B %*% Y.hat[, j]))
        if (!is.null(operator$BS)) {
            stabilizer.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$BS %*% Y.hat[, j]))
        } else {
            stabilizer.energy[j] <- 0
        }
        objective[j] <- data.loss[j] +
            0.5 * lambda1 * hessian.energy[j] +
            0.5 * lambda2 * stabilizer.energy[j]
    }

    if (!is.null(y.info$col.names)) {
        colnames(Y.hat) <- y.info$col.names
        colnames(residuals) <- y.info$col.names
        names(data.loss) <- y.info$col.names
        names(hessian.energy) <- y.info$col.names
        names(stabilizer.energy) <- y.info$col.names
        names(objective) <- y.info$col.names
    }

    single <- n.responses == 1L
    list(
        fitted.values = if (single) as.vector(Y.hat) else Y.hat,
        residuals = if (single) as.vector(residuals) else residuals,
        y = if (single) as.vector(Y.raw) else Y.raw,
        weights = if (single) as.vector(W) else W,
        lambda = list(lambda1 = lambda1, lambda2 = lambda2, ridge = ridge),
        objective = if (single) objective[1L] else objective,
        energies = list(
            data.loss = if (single) data.loss[1L] else data.loss,
            hessian = if (single) hessian.energy[1L] else hessian.energy,
            stabilizer = if (single) stabilizer.energy[1L] else stabilizer.energy
        ),
        operator = operator,
        n.responses = n.responses,
        solver = list(method = if (length(unique(solve.method)) == 1L) solve.method[1L] else solve.method),
        observed = if (single) as.vector(y.info$observed & W > 0) else y.info$observed & W > 0
    )
}

.solve.ssrhe.hessian.system <- function(B, BS, weights, Y,
                                        lambda1, lambda2, ridge,
                                        verbose = FALSE) {
    n <- nrow(B)
    if (length(weights) != n) stop("Internal error: weights length mismatch.", call. = FALSE)
    system <- lambda1 * B
    if (lambda2 > 0) {
        system <- system + lambda2 * BS
    }
    if (ridge > 0) {
        system <- system + Matrix::Diagonal(n, ridge)
    }
    if (any(weights > 0)) {
        system <- system + Matrix::Diagonal(n, weights)
    }
    rhs <- weights * Y

    fit <- tryCatch({
        as.matrix(Matrix::solve(system, rhs))
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse solve failed; retrying with dense solve().")
        }
        tryCatch(
            solve(as.matrix(system), as.matrix(rhs)),
            error = function(e2) {
                stop("SSRHE linear solve failed: ", conditionMessage(e2),
                     call. = FALSE)
            }
        )
    })
    list(fitted = fit, method = "linear.solve")
}

.prepare.ssrhe.response.matrix <- function(y, n, name) {
    is.matrix.input <- is.matrix(y) || inherits(y, "Matrix") ||
        (is.data.frame(y) && ncol(y) > 1L)
    if (is.matrix.input) {
        Y <- if (inherits(y, "Matrix")) as.matrix(y) else as.matrix(y)
        if (nrow(Y) != n) stop(sprintf("nrow(%s) must be %d.", name, n), call. = FALSE)
        col.names <- colnames(Y)
    } else {
        if (length(y) != n) stop(sprintf("%s must have length %d.", name, n), call. = FALSE)
        Y <- matrix(y, ncol = 1L)
        col.names <- NULL
    }
    storage.mode(Y) <- "double"
    if (any(is.infinite(Y))) {
        stop(sprintf("%s cannot contain infinite values.", name), call. = FALSE)
    }
    observed <- is.finite(Y)
    list(Y = Y, observed = observed, col.names = col.names)
}

.prepare.ssrhe.weight.matrix <- function(weights, y.info, n) {
    p <- ncol(y.info$Y)
    if (is.null(weights)) {
        W <- matrix(1, nrow = n, ncol = p)
    } else if (is.matrix(weights) || inherits(weights, "Matrix") || is.data.frame(weights)) {
        W <- if (inherits(weights, "Matrix")) as.matrix(weights) else as.matrix(weights)
        if (!identical(dim(W), dim(y.info$Y))) {
            stop("weights matrix must have the same dimensions as y.", call. = FALSE)
        }
    } else {
        if (length(weights) != n) stop("weights vector must have length nrow(X).", call. = FALSE)
        W <- matrix(as.double(weights), nrow = n, ncol = p)
    }
    storage.mode(W) <- "double"
    if (any(!is.finite(W)) || any(W < 0)) {
        stop("weights must be finite and nonnegative.", call. = FALSE)
    }
    W[!y.info$observed] <- 0
    W
}

.validate.ssrhe.nonnegative.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
        stop(sprintf("%s must be a finite nonnegative numeric scalar.", name),
             call. = FALSE)
    }
    as.double(x)
}

#' @export
print.ssrhe.hessian.fit <- function(x, ...) {
    cat("SSRHE Hessian regression fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda1:", format(x$lambda$lambda1, digits = 4), "\n")
    cat("  lambda2:", format(x$lambda$lambda2, digits = 4), "\n")
    cat("  objective:", paste(format(x$objective, digits = 4), collapse = ", "), "\n")
    invisible(x)
}

#' @export
print.ssrhe.hessian.refit <- function(x, ...) {
    cat("SSRHE Hessian regression refit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda1:", format(x$lambda$lambda1, digits = 4), "\n")
    cat("  lambda2:", format(x$lambda$lambda2, digits = 4), "\n")
    invisible(x)
}

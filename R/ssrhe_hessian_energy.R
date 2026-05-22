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
#'   itself, used when \code{neighborhood.type = "knn"}. This follows the
#'   SSRHE Matlab convention. For adaptive-radius neighborhoods, \code{k} may be
#'   omitted when \code{adaptive.k.scale} is supplied.
#' @param tangent.dim Integer tangent dimension used for the local PCA chart.
#'   Required when \code{tangent.dim.rule = "fixed"}. If \code{NULL} and
#'   \code{tangent.dim.rule = "eigen.cumulative"}, the dimension is selected
#'   locally from the PCA variance ratio.
#' @param nn.index Optional integer matrix of neighbor indices, with
#'   \code{nrow(X)} rows and \code{k} columns. Each row must contain its center
#'   vertex. If \code{NULL}, Euclidean k-NN including self is computed in C++.
#' @param neighborhood.type Local support rule. \code{"knn"} uses the original
#'   rectangular self-including kNN neighborhoods. \code{"adaptive.radius"}
#'   builds variable-size supports from \code{\link{create.adaptive.radius.graph}}.
#'   \code{"supplied"} uses \code{support.index} directly.
#' @param support.index Optional list of integer vectors, one per row of
#'   \code{X}. Each element gives a variable-size local support and must contain
#'   its center vertex.
#' @param adaptive.k.scale Integer local-scale k used by
#'   \code{\link{create.adaptive.radius.graph}} when
#'   \code{neighborhood.type = "adaptive.radius"}.
#' @param radius.rule,radius.factor Adaptive-radius graph parameters passed to
#'   \code{\link{create.adaptive.radius.graph}}.
#' @param min.support Optional minimum local support size for adaptive-radius
#'   supports. Defaults to the local quadratic design size plus
#'   \code{support.buffer}, clamped to \code{nrow(X)}.
#' @param max.support Optional maximum local support size. Oversized adaptive
#'   supports are trimmed to the closest \code{max.support} vertices while
#'   keeping the center vertex.
#' @param support.buffer Nonnegative integer added to the local quadratic design
#'   size when choosing the default \code{min.support}.
#' @param support.topup How undersized adaptive-radius supports are enlarged.
#'   \code{"nearest"} appends ambient nearest neighbors. \code{"none"} leaves
#'   supports unchanged and lets the C++ operator reject undersized supports.
#' @param tangent.dim.rule Either \code{"fixed"} or \code{"eigen.cumulative"}.
#' @param eigen.tolerance Cumulative local PCA variance threshold used when
#'   \code{tangent.dim.rule = "eigen.cumulative"} and \code{tangent.dim = NULL}.
#' @param derivative.order Integer derivative order of the local SSRHE-style
#'   operator. \code{2L} is the original local Hessian energy. \code{3L} is an
#'   experimental third-derivative energy whose rows estimate unique symmetric
#'   third-derivative tensor components with factorial and tensor-multiplicity
#'   scaling.
#' @param stabilizer Logical. If \code{TRUE}, also construct the supplemental
#'   stabilizer matrix described in the SSRHE supplement. Currently supported
#'   only for \code{derivative.order = 2L}.
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
#' For each point \eqn{x_i}, the operator constructs a local PCA chart from a
#' self-including local support, projects the support into
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
#' With \code{derivative.order = 3L}, the reduced local design uses cubic
#' monomials first, followed by quadratic and linear monomials. The returned
#' rows estimate unique symmetric third-derivative components
#' \eqn{\partial_{abc} f} with scale
#' \deqn{
#'   s_{abc}=m_{abc}\sqrt{\mu_{abc}},
#' }
#' where \eqn{m_{abc}} is the factorial monomial-to-derivative multiplier
#' (\eqn{6} for \eqn{aaa}, \eqn{2} for \eqn{aab}, and \eqn{1} for three
#' distinct indices) and \eqn{\mu_{abc}} is the ordered-tensor multiplicity
#' (\eqn{1}, \eqn{3}, or \eqn{6}). Thus \eqn{\|Af\|_2^2} approximates the
#' squared Frobenius norm of the full symmetric third-derivative tensor.
#'
#' The returned \eqn{A} is the stacked sparse matrix of these local derivative
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
#'     \item \code{nn.index}: neighbor index matrix for rectangular kNN
#'       neighborhoods, if used;
#'     \item \code{support.index}: list of actual local supports used by the
#'       backend;
#'     \item \code{neighborhoods}: support-construction diagnostics;
#'     \item \code{row.table}: one row per Hessian component row of \eqn{A};
#'     \item \code{diagnostics}: local design rank, condition, and PCA summaries;
#'     \item \code{local.diagnostics}: one row per center vertex with support
#'       size, design rank/condition, PCA variance, chart-distortion,
#'       boundary-asymmetry, and curvature-bias proxy diagnostics;
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
    k = NULL,
    tangent.dim,
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = FALSE,
    pinv.tol = sqrt(.Machine$double.eps),
    return.A = TRUE,
    return.B = TRUE,
    return.BS = stabilizer,
    return.sparse = TRUE,
    verbose = FALSE) {

    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    neighborhood.type <- match.arg(neighborhood.type)
    radius.rule <- match.arg(radius.rule)
    support.topup <- match.arg(support.topup)

    tangent.dim.rule <- match.arg(tangent.dim.rule)
    if (missing(tangent.dim) || is.null(tangent.dim)) {
        if (identical(tangent.dim.rule, "fixed")) {
            stop("tangent.dim is required when tangent.dim.rule = 'fixed'.",
                 call. = FALSE)
        }
        tangent.dim.cpp <- 0L
    } else {
        tangent.dim.cpp <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        if (tangent.dim.cpp > ncol(X)) {
            stop("tangent.dim must be no larger than ncol(X).",
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
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    stabilizer <- isTRUE(stabilizer)
    if (derivative.order == 3L && stabilizer) {
        stop("stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }
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

    neighborhood <- .prepare.ssrhe.neighborhood(
        X = X,
        k = k,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
        tangent.dim.rule = tangent.dim.rule,
        derivative.order = derivative.order
    )

    if (!missing(tangent.dim) && !is.null(tangent.dim) &&
        tangent.dim.cpp > min(min(neighborhood$support.size), ncol(X))) {
        stop("tangent.dim must be no larger than the smallest local support size and ncol(X).",
             call. = FALSE)
    }

    raw <- .Call(
        "S_ssrhe_hessian_operator",
        X,
        as.integer(neighborhood$k),
        as.integer(tangent.dim.cpp),
        neighborhood$nn.index,
        neighborhood$support.index,
        tangent.dim.rule,
        as.double(eigen.tolerance),
        as.integer(derivative.order),
        as.logical(stabilizer),
        as.double(pinv.tol),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    raw$row.table <- as.data.frame(raw$row.table, stringsAsFactors = FALSE)
    raw$diagnostics <- as.data.frame(raw$diagnostics, stringsAsFactors = FALSE)
    raw$local.diagnostics <- .ssrhe.local.geometry.diagnostics(
        X = X,
        support.index = raw$support.index,
        diagnostics = raw$diagnostics
    )

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
        support.index = raw$support.index,
        neighborhoods = neighborhood$metadata,
        row.table = raw$row.table,
        diagnostics = raw$diagnostics,
        local.diagnostics = raw$local.diagnostics,
        parity = raw$parity,
        parameters = raw$parameters
    )
    out$parameters$neighborhood.type <- neighborhood$type
    out$parameters$k <- neighborhood$k
    out$parameters$adaptive.k.scale <- neighborhood$adaptive.k.scale
    out$parameters$radius.rule <- neighborhood$radius.rule
    out$parameters$radius.factor <- neighborhood$radius.factor
    out$parameters$derivative.order <- derivative.order

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

.validate.ssrhe.derivative.order <- function(x) {
    x <- .validate.ssrhe.positive.integer(x, "derivative.order")
    if (!x %in% c(2L, 3L)) {
        stop("derivative.order must be either 2L or 3L.", call. = FALSE)
    }
    x
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

.validate.ssrhe.support.index <- function(support.index, n) {
    if (is.null(support.index)) {
        return(NULL)
    }
    if (is.matrix(support.index)) {
        support.index <- lapply(seq_len(nrow(support.index)), function(i) {
            support.index[i, !is.na(support.index[i, ])]
        })
    }
    if (!is.list(support.index) || length(support.index) != n) {
        stop("support.index must be a list with length nrow(X).", call. = FALSE)
    }
    out <- vector("list", n)
    for (i in seq_len(n)) {
        ids <- as.integer(support.index[[i]])
        if (length(ids) < 2L || any(is.na(ids)) || any(ids < 1L) || any(ids > n)) {
            stop("Each support.index element must contain at least two valid vertex indices.",
                 call. = FALSE)
        }
        ids <- unique(ids)
        if (!i %in% ids) {
            stop("Each support.index element must contain its center vertex.",
                 call. = FALSE)
        }
        out[[i]] <- ids
    }
    out
}

.ssrhe.design.ncol <- function(m, derivative.order = 2L) {
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    q2 <- m * (m + 1L) / 2L
    if (derivative.order == 2L) {
        q2 + m
    } else {
        q3 <- m * (m + 1L) * (m + 2L) / 6L
        q3 + q2 + m
    }
}

.ssrhe.default.min.support <- function(tangent.dim,
                                       tangent.dim.rule,
                                       ambient.dim,
                                       support.buffer,
                                       n,
                                       derivative.order = 2L) {
    support.buffer <- .validate.ssrhe.nonnegative.integer(support.buffer,
                                                         "support.buffer")
    if (!is.null(tangent.dim)) {
        m <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
    } else if (identical(tangent.dim.rule, "eigen.cumulative")) {
        m <- ambient.dim
    } else {
        stop("tangent.dim is required to choose the default min.support.",
             call. = FALSE)
    }
    min(n, .ssrhe.design.ncol(m, derivative.order = derivative.order) + support.buffer)
}

#' Build Candidate Support Profiles for SSRHE Adaptive-Radius Tuning
#'
#' Constructs a compact grid of adaptive-radius support profiles for use with
#' \code{support.selection = "cv"} in SSRHE fitting functions. Each row gives an
#' \code{adaptive.k.scale} value and a requested \code{min.support}. The grid is
#' intentionally small by default because support selection nests operator
#' construction inside response cross-validation.
#'
#' @param n Number of observations.
#' @param tangent.dim Tangent dimension used by the SSRHE local polynomial
#'   design.
#' @param derivative.order SSRHE derivative order, currently \code{2L} or
#'   \code{3L}.
#' @param support.buffer Nonnegative integer added to the local design size when
#'   constructing the smallest candidate support.
#' @param max.candidates Maximum number of candidate support profiles to return.
#'
#' @return A data frame with columns \code{adaptive.k.scale},
#'   \code{min.support}, and \code{max.support}. \code{max.support} is
#'   \code{NA} by default, meaning no truncation.
#' @export
ssrhe.support.grid <- function(n,
                               tangent.dim,
                               derivative.order = 2L,
                               support.buffer = 2L,
                               max.candidates = 8L) {
    n <- .validate.ssrhe.positive.integer(n, "n")
    tangent.dim <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    support.buffer <- .validate.ssrhe.nonnegative.integer(support.buffer,
                                                         "support.buffer")
    max.candidates <- .validate.ssrhe.positive.integer(max.candidates,
                                                      "max.candidates")
    if (n < 3L) {
        stop("n must be at least 3.", call. = FALSE)
    }
    base <- min(n - 1L,
                .ssrhe.design.ncol(tangent.dim, derivative.order) +
                    support.buffer)
    candidates <- unique(as.integer(round(c(
        base,
        base + c(3L, 5L, 10L),
        1.5 * base,
        2 * base,
        3 * base,
        4 * base,
        floor(0.15 * n),
        floor(0.25 * n)
    ))))
    candidates <- sort(unique(pmax(2L, pmin(candidates, n - 1L))))
    if (length(candidates) > max.candidates) {
        keep <- unique(round(seq(1, length(candidates),
                                 length.out = max.candidates)))
        candidates <- candidates[keep]
    }
    adaptive.k.scale <- pmax(2L, floor(candidates / 3L))
    adaptive.k.scale <- pmin(adaptive.k.scale, n - 1L)
    unique(data.frame(
        adaptive.k.scale = as.integer(adaptive.k.scale),
        min.support = as.integer(candidates),
        max.support = NA_integer_
    ))
}

.validate.ssrhe.support.grid <- function(support.grid,
                                         n,
                                         tangent.dim,
                                         derivative.order,
                                         support.buffer,
                                         max.candidates) {
    if (is.null(support.grid)) {
        return(ssrhe.support.grid(
            n = n,
            tangent.dim = tangent.dim,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = max.candidates
        ))
    }
    if (!is.data.frame(support.grid)) {
        stop("support.grid must be a data frame.", call. = FALSE)
    }
    required <- c("adaptive.k.scale", "min.support")
    if (!all(required %in% names(support.grid))) {
        stop("support.grid must contain adaptive.k.scale and min.support columns.",
             call. = FALSE)
    }
    grid <- support.grid[, unique(c(required, intersect("max.support", names(support.grid)))),
                         drop = FALSE]
    if (!"max.support" %in% names(grid)) {
        grid$max.support <- NA_integer_
    }
    grid$adaptive.k.scale <- as.integer(round(grid$adaptive.k.scale))
    grid$min.support <- as.integer(round(grid$min.support))
    grid$max.support <- ifelse(is.na(grid$max.support), NA_integer_,
                               as.integer(round(grid$max.support)))
    bad <- is.na(grid$adaptive.k.scale) | grid$adaptive.k.scale < 1L |
        grid$adaptive.k.scale >= n |
        is.na(grid$min.support) | grid$min.support < 2L |
        grid$min.support > n |
        (!is.na(grid$max.support) &
             (grid$max.support < grid$min.support | grid$max.support > n))
    if (any(bad)) {
        stop("support.grid contains invalid adaptive.k.scale/min.support/max.support values.",
             call. = FALSE)
    }
    grid <- unique(grid)
    if (nrow(grid) > max.candidates) {
        grid <- grid[seq_len(max.candidates), , drop = FALSE]
    }
    rownames(grid) <- NULL
    grid
}

.ssrhe.support.diagnostics <- function(operator) {
    support.size <- operator$neighborhoods$support.size
    topup <- operator$neighborhoods$n.topup
    truncated <- operator$neighborhoods$n.truncated
    data.frame(
        support.size.min = min(support.size),
        support.size.median = stats::median(support.size),
        support.size.max = max(support.size),
        support.topup.median = stats::median(topup, na.rm = TRUE),
        support.topup.max = max(topup, na.rm = TRUE),
        support.truncated.max = max(truncated, na.rm = TRUE),
        operator.rows = operator$A.triplet$dim[[1L]],
        operator.cols = operator$A.triplet$dim[[2L]],
        stringsAsFactors = FALSE
    )
}

.ssrhe.local.geometry.diagnostics <- function(X, support.index, diagnostics) {
    n <- nrow(X)
    if (!is.list(support.index) || length(support.index) != n) {
        stop("Internal error: invalid SSRHE support.index.", call. = FALSE)
    }
    out <- diagnostics
    support.size <- out$k
    if (is.null(support.size)) {
        support.size <- lengths(support.index)
    }
    out$support.size <- support.size
    out$design.rank.deficiency <- pmax(0L, out$design.ncol - out$design.rank)
    out$pca.variance.explained <- out$local.variance.ratio
    out$pca.discarded.variance.ratio <- pmax(0, 1 - out$local.variance.ratio)

    chart.distortion <- rep(NA_real_, n)
    boundary.asymmetry <- rep(NA_real_, n)
    curvature.bias.proxy <- out$pca.discarded.variance.ratio

    for (center in seq_len(n)) {
        ids <- as.integer(support.index[[center]])
        if (!length(ids) || !center %in% ids) next
        local <- X[ids, , drop = FALSE]
        centered <- sweep(local, 2L, colMeans(local), "-")
        svd.local <- tryCatch(svd(centered), error = function(e) NULL)
        if (is.null(svd.local) || !length(svd.local$d)) next
        m <- as.integer(out$tangent.dim[center])
        m <- max(1L, min(m, ncol(svd.local$v), length(svd.local$d)))
        coords <- centered %*% svd.local$v[, seq_len(m), drop = FALSE]
        base.local <- match(center, ids)
        coords <- sweep(coords, 2L, coords[base.local, ], "-")

        radius <- sqrt(max(rowSums(coords^2), na.rm = TRUE))
        centroid.norm <- sqrt(sum(colMeans(coords)^2))
        boundary.asymmetry[center] <- if (is.finite(radius) && radius > 0) {
            centroid.norm / radius
        } else {
            0
        }

        if (length(ids) >= 3L) {
            d.ambient <- stats::dist(local)
            d.chart <- stats::dist(coords)
            d.ambient <- as.numeric(d.ambient)
            d.chart <- as.numeric(d.chart)
            ok <- is.finite(d.ambient) & is.finite(d.chart) & d.ambient > 0
            denom <- mean(d.ambient[ok])
            chart.distortion[center] <- if (any(ok) && is.finite(denom) && denom > 0) {
                sqrt(mean((d.chart[ok] - d.ambient[ok])^2)) / denom
            } else {
                0
            }
        } else {
            chart.distortion[center] <- 0
        }
    }

    out$chart.distortion <- chart.distortion
    out$boundary.asymmetry <- boundary.asymmetry
    out$curvature.bias.proxy <- curvature.bias.proxy
    out
}

.validate.ssrhe.nonnegative.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 ||
        abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a nonnegative integer scalar.", name),
             call. = FALSE)
    }
    as.integer(round(x))
}

.prepare.ssrhe.neighborhood <- function(X,
                                        k,
                                        nn.index,
                                        neighborhood.type,
                                        support.index,
                                        adaptive.k.scale,
                                        radius.rule,
                                        radius.factor,
                                        min.support,
                                        max.support,
                                        support.buffer,
                                        support.topup,
                                        tangent.dim,
                                        tangent.dim.rule,
                                        derivative.order = 2L) {
    n <- nrow(X)
    if (!is.null(support.index) &&
        identical(neighborhood.type, "knn") &&
        is.null(nn.index) && is.null(k)) {
        neighborhood.type <- "supplied"
    }

    if (identical(neighborhood.type, "knn")) {
        if (is.null(k)) {
            if (!is.null(nn.index)) {
                k <- ncol(as.matrix(nn.index))
            } else {
                stop("k is required when neighborhood.type = 'knn'.",
                     call. = FALSE)
            }
        }
        k <- .validate.ssrhe.positive.integer(k, "k")
        if (k < 2L || k > n) {
            stop("k must be between 2 and nrow(X).", call. = FALSE)
        }
        nn.index <- .validate.ssrhe.nn.index(nn.index, n, k)
        support.size <- rep(k, n)
        return(list(
            type = "knn",
            k = k,
            nn.index = nn.index,
            support.index = NULL,
            support.size = support.size,
            adaptive.k.scale = NA_integer_,
            radius.rule = NA_character_,
            radius.factor = NA_real_,
            metadata = list(
                type = "knn",
                support.size = support.size,
                n.topup = rep(0L, n),
                n.truncated = rep(0L, n)
            )
        ))
    }

    if (identical(neighborhood.type, "supplied")) {
        support.index <- .validate.ssrhe.support.index(support.index, n)
        support.size <- lengths(support.index)
        return(list(
            type = "supplied",
            k = 0L,
            nn.index = NULL,
            support.index = support.index,
            support.size = support.size,
            adaptive.k.scale = NA_integer_,
            radius.rule = NA_character_,
            radius.factor = NA_real_,
            metadata = list(
                type = "supplied",
                support.size = support.size,
                n.topup = rep(0L, n),
                n.truncated = rep(0L, n)
            )
        ))
    }

    if (!identical(neighborhood.type, "adaptive.radius")) {
        stop("Unknown SSRHE neighborhood type.", call. = FALSE)
    }
    if (is.null(adaptive.k.scale)) {
        if (is.null(k)) {
            stop("adaptive.k.scale is required when neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        adaptive.k.scale <- k
    }
    adaptive.k.scale <- .validate.ssrhe.positive.integer(adaptive.k.scale,
                                                        "adaptive.k.scale")
    if (adaptive.k.scale >= n) {
        stop("adaptive.k.scale must be smaller than nrow(X).", call. = FALSE)
    }
    radius.factor <- .validate.ssrhe.numeric.scalar(radius.factor, "radius.factor")
    if (radius.factor <= 0) {
        stop("radius.factor must be positive.", call. = FALSE)
    }
    if (is.null(min.support)) {
        min.support <- .ssrhe.default.min.support(
            tangent.dim = tangent.dim,
            tangent.dim.rule = tangent.dim.rule,
            ambient.dim = ncol(X),
            support.buffer = support.buffer,
            n = n,
            derivative.order = derivative.order
        )
    } else {
        min.support <- .validate.ssrhe.positive.integer(min.support,
                                                       "min.support")
    }
    if (min.support < 2L || min.support > n) {
        stop("min.support must be between 2 and nrow(X).", call. = FALSE)
    }
    if (!is.null(max.support)) {
        max.support <- .validate.ssrhe.positive.integer(max.support,
                                                       "max.support")
        if (max.support < min.support || max.support > n) {
            stop("max.support must be between min.support and nrow(X).",
                 call. = FALSE)
        }
    }

    graph <- create.adaptive.radius.graph(
        X = X,
        k.scale = adaptive.k.scale,
        radius.factor = radius.factor,
        radius.rule = radius.rule,
        prune.method = "none",
        connect.components = FALSE
    )
    nearest <- NULL
    if (identical(support.topup, "nearest") || !is.null(max.support)) {
        nearest <- .exact.knn.index(X, n - 1L)
    }
    support.index <- vector("list", n)
    n.topup <- integer(n)
    n.truncated <- integer(n)
    for (i in seq_len(n)) {
        ids <- unique(c(i, graph$adj_list[[i]]))
        if (length(ids) < min.support && identical(support.topup, "nearest")) {
            add <- nearest[i, !nearest[i, ] %in% ids]
            need <- min.support - length(ids)
            ids <- unique(c(ids, add[seq_len(min(need, length(add)))]))
            n.topup[[i]] <- max(0L, length(ids) - length(unique(c(i, graph$adj_list[[i]]))))
        }
        if (!is.null(max.support) && length(ids) > max.support) {
            non.center <- setdiff(ids, i)
            d <- rowSums((t(t(X[non.center, , drop = FALSE]) - X[i, ]))^2)
            keep <- non.center[order(d, non.center)][seq_len(max.support - 1L)]
            n.truncated[[i]] <- length(ids) - max.support
            ids <- c(i, keep)
        } else if (ids[[1L]] != i) {
            ids <- c(i, setdiff(ids, i))
        }
        support.index[[i]] <- as.integer(ids)
    }
    support.index <- .validate.ssrhe.support.index(support.index, n)
    support.size <- lengths(support.index)
    list(
        type = "adaptive.radius",
        k = 0L,
        nn.index = NULL,
        support.index = support.index,
        support.size = support.size,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        metadata = list(
            type = "adaptive.radius",
            graph = graph,
            adaptive.k.scale = adaptive.k.scale,
            radius.rule = radius.rule,
            radius.factor = radius.factor,
            min.support = min.support,
            max.support = max.support,
            support.topup = support.topup,
            support.size = support.size,
            n.topup = n.topup,
            n.truncated = n.truncated,
            sigma = graph$sigma
        )
    )
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
    cat("  neighborhood.type:", x$parameters$neighborhood.type %||% "knn", "\n")
    if (identical(x$parameters$neighborhood.type, "knn")) {
        cat("  k:", x$parameters$k, "\n")
    } else if (!is.null(x$neighborhoods$support.size)) {
        cat("  support size:",
            paste(range(x$neighborhoods$support.size), collapse = "-"),
            "\n")
    }
    cat("  tangent.dim.rule:", x$parameters$tangent.dim.rule, "\n")
    cat("  derivative.order:", x$parameters$derivative.order %||% 2L, "\n")
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
#'   Currently supported only for \code{derivative.order = 2L}.
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
    k = NULL,
    tangent.dim,
    lambda1,
    lambda2 = 0,
    weights = NULL,
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = lambda2 > 0,
    pinv.tol = sqrt(.Machine$double.eps),
    ridge = 0,
    return.A = TRUE,
    verbose = FALSE) {

    X <- .validate.ssrhe.X(X)
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    stabilizer <- isTRUE(stabilizer)
    if (lambda2 > 0 && !stabilizer) {
        stop("lambda2 > 0 requires stabilizer = TRUE.", call. = FALSE)
    }
    if (derivative.order == 3L && (lambda2 > 0 || stabilizer)) {
        stop("lambda2/stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
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

#' Select SSRHE Hessian Regression Penalties by Label Cross-Validation
#'
#' Fits \code{\link{fit.ssrhe.hessian.regression}} over a grid of fixed
#' \code{lambda1}/\code{lambda2} values using cross-validation on observed
#' labels. This is intended for semi-supervised SSRHE use: validation folds are
#' formed only from entries that are observed and have positive data-fit weight.
#'
#' @inheritParams fit.ssrhe.hessian.regression
#' @param lambda1.grid Nonnegative numeric vector of Hessian-energy penalty
#'   candidates.
#' @param lambda2.grid Nonnegative numeric vector of supplemental-stabilizer
#'   penalty candidates. Use \code{0} to omit the supplemental stabilizer from
#'   selection.
#' @param nfolds Number of validation folds over observed positive-weight
#'   labels. Ignored when \code{fold.id} is supplied.
#' @param fold.id Optional integer vector of length \code{nrow(X)} assigning
#'   observed positive-weight labels to validation folds. Nonpositive or
#'   \code{NA} entries are ignored.
#' @param loss Validation loss, currently \code{"mse"} or \code{"mae"}.
#' @param selection Selection rule. \code{"min"} chooses the smallest mean
#'   validation loss. \code{"one.se"} chooses the largest total penalty among
#'   candidates within one standard error of the minimum.
#' @param support.selection Support-profile selection rule.
#'   \code{"rule"} uses the supplied \code{adaptive.k.scale},
#'   \code{min.support}, and \code{max.support}. \code{"cv"} is currently
#'   supported for \code{neighborhood.type = "adaptive.radius"} and chooses
#'   among rows of \code{support.grid} by outer response cross-validation.
#' @param support.grid Optional data frame of support profiles with columns
#'   \code{adaptive.k.scale}, \code{min.support}, and optional
#'   \code{max.support}. If \code{NULL}, \code{\link{ssrhe.support.grid}} builds
#'   a compact default grid.
#' @param support.cv.max.candidates Maximum number of support profiles to try
#'   when \code{support.selection = "cv"}.
#'
#' @details
#' For each validation fold, the held-out labels are removed from the data-fit
#' term by setting their weights to zero. The fitted values are then scored only
#' on those held-out labels. The final returned fit is refit with the selected
#' penalties using all observed positive-weight labels.
#'
#' With \code{support.selection = "cv"}, this function performs an outer
#' support-profile selection loop. For each candidate adaptive-radius support
#' profile, it constructs a fresh SSRHE operator, runs the usual
#' \code{lambda1.grid}/\code{lambda2.grid} cross-validation with the same fold
#' assignments, and selects the support profile with the smallest selected CV
#' error. This can substantially increase runtime because local operator
#' construction and lambda CV are nested.
#'
#' The current implementation supports a single response vector. Matrix-response
#' penalty selection should be performed column-by-column.
#'
#' @return A list of class \code{"ssrhe.hessian.cv.fit"} and
#'   \code{"ssrhe.hessian.fit"} containing the final fit plus
#'   \code{cv.table}, \code{fold.id}, and \code{selection} diagnostics.
#' @export
fit.ssrhe.hessian.regression.cv <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda1.grid,
    lambda2.grid = 0,
    weights = NULL,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = any(lambda2.grid > 0),
    pinv.tol = sqrt(.Machine$double.eps),
    ridge = 0,
    return.A = TRUE,
    support.selection = c("rule", "cv"),
    support.grid = NULL,
    support.cv.max.candidates = 8L,
    verbose = FALSE) {

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE CV regression.",
             call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("fit.ssrhe.hessian.regression.cv currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda1.grid <- .validate.ssrhe.lambda.grid(lambda1.grid, "lambda1.grid")
    lambda2.grid <- .validate.ssrhe.lambda.grid(lambda2.grid, "lambda2.grid")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    support.selection <- match.arg(support.selection)
    neighborhood.type <- match.arg(neighborhood.type)
    support.cv.max.candidates <- .validate.ssrhe.positive.integer(
        support.cv.max.candidates, "support.cv.max.candidates")
    stabilizer <- isTRUE(stabilizer)
    if (any(lambda2.grid > 0) && !stabilizer) {
        stop("positive lambda2.grid values require stabilizer = TRUE.",
             call. = FALSE)
    }
    if (derivative.order == 3L && (any(lambda2.grid > 0) || stabilizer)) {
        stop("lambda2.grid/stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }
    if (identical(support.selection, "cv")) {
        if (!identical(neighborhood.type, "adaptive.radius")) {
            stop("support.selection = 'cv' is currently supported only for neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        tangent.dim.grid <- if (missing(tangent.dim) || is.null(tangent.dim)) {
            ncol(X)
        } else {
            .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        }
        support.grid <- .validate.ssrhe.support.grid(
            support.grid = support.grid,
            n = n,
            tangent.dim = tangent.dim.grid,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = support.cv.max.candidates
        )
        fold.id.cv <- .prepare.ssrhe.cv.folds(
            observed = as.vector(y.info$observed & W > 0),
            nfolds = nfolds,
            fold.id = fold.id
        )
        candidate.fits <- vector("list", nrow(support.grid))
        candidate.rows <- vector("list", nrow(support.grid))
        for (ii in seq_len(nrow(support.grid))) {
            max.support.ii <- support.grid$max.support[ii]
            if (is.na(max.support.ii)) max.support.ii <- NULL
            elapsed <- system.time({
                cand <- tryCatch(
                    fit.ssrhe.hessian.regression.cv(
                        X = X,
                        y = y,
                        k = k,
                        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
                        lambda1.grid = lambda1.grid,
                        lambda2.grid = lambda2.grid,
                        weights = weights,
                        nfolds = nfolds,
                        fold.id = fold.id.cv,
                        loss = loss,
                        selection = selection,
                        support.selection = "rule",
                        nn.index = nn.index,
                        neighborhood.type = "adaptive.radius",
                        support.index = support.index,
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        radius.rule = radius.rule,
                        radius.factor = radius.factor,
                        min.support = support.grid$min.support[ii],
                        max.support = max.support.ii,
                        support.buffer = support.buffer,
                        support.topup = support.topup,
                        tangent.dim.rule = tangent.dim.rule,
                        eigen.tolerance = eigen.tolerance,
                        derivative.order = derivative.order,
                        stabilizer = stabilizer,
                        pinv.tol = pinv.tol,
                        ridge = ridge,
                        return.A = return.A,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            })[["elapsed"]]
            if (inherits(cand, "error")) {
                candidate.rows[[ii]] <- data.frame(
                    candidate = ii,
                    status = "error",
                    adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                    min.support = support.grid$min.support[ii],
                    max.support = support.grid$max.support[ii],
                    cv.mean = Inf,
                    lambda1 = NA_real_,
                    lambda2 = NA_real_,
                    runtime.sec = unname(elapsed),
                    message = conditionMessage(cand),
                    stringsAsFactors = FALSE
                )
            } else {
                candidate.fits[[ii]] <- cand
                diag <- .ssrhe.support.diagnostics(cand$operator)
                candidate.rows[[ii]] <- cbind(
                    data.frame(
                        candidate = ii,
                        status = "ok",
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        min.support = support.grid$min.support[ii],
                        max.support = support.grid$max.support[ii],
                        cv.mean = cand$selection$cv.mean,
                        cv.se = cand$selection$cv.se,
                        lambda1 = cand$selection$lambda1,
                        lambda2 = cand$selection$lambda2,
                        runtime.sec = unname(elapsed),
                        message = "",
                        stringsAsFactors = FALSE
                    ),
                    diag
                )
            }
        }
        support.cv.table <- do.call(rbind, candidate.rows)
        if (!any(is.finite(support.cv.table$cv.mean))) {
            stop("All SSRHE support-CV candidates failed.", call. = FALSE)
        }
        selected.support.idx <- which.min(support.cv.table$cv.mean)
        final <- candidate.fits[[selected.support.idx]]
        final$support.selection <- "cv"
        final$support.cv.table <- support.cv.table
        final$selected.support <- support.grid[selected.support.idx, , drop = FALSE]
        final$selection$support.index <- selected.support.idx
        final$selection$support.cv.mean <- support.cv.table$cv.mean[selected.support.idx]
        attr(final, "call") <- match.call()
        class(final) <- c("ssrhe.hessian.cv.fit", "ssrhe.hessian.fit", "list")
        return(final)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = stabilizer,
        pinv.tol = pinv.tol,
        return.A = return.A,
        return.B = TRUE,
        return.BS = stabilizer,
        return.sparse = TRUE,
        verbose = verbose
    )

    fold.id <- .prepare.ssrhe.cv.folds(
        observed = as.vector(y.info$observed & W > 0),
        nfolds = nfolds,
        fold.id = fold.id
    )
    folds <- sort(unique(fold.id[fold.id > 0]))
    grid <- expand.grid(lambda1 = lambda1.grid,
                        lambda2 = lambda2.grid,
                        KEEP.OUT.ATTRS = FALSE)
    fold.loss <- matrix(NA_real_, nrow = nrow(grid), ncol = length(folds))
    colnames(fold.loss) <- paste0("fold", folds)

    Y <- as.vector(y.info$Y)
    for (g in seq_len(nrow(grid))) {
        for (ff in seq_along(folds)) {
            validation <- fold.id == folds[ff]
            train.weights <- as.vector(W)
            train.weights[validation] <- 0
            fold.fit <- tryCatch(
                .fit.ssrhe.hessian.from.operator(
                    operator = operator,
                    y = Y,
                    lambda1 = grid$lambda1[g],
                    lambda2 = grid$lambda2[g],
                    weights = train.weights,
                    ridge = ridge,
                    verbose = verbose
                ),
                error = function(e) NULL
            )
            if (!is.null(fold.fit)) {
                err <- fold.fit$fitted.values[validation] - Y[validation]
                wv <- as.vector(W)[validation]
                if (identical(loss, "mse")) {
                    fold.loss[g, ff] <- sum(wv * err^2) / sum(wv)
                } else {
                    fold.loss[g, ff] <- sum(wv * abs(err)) / sum(wv)
                }
            }
        }
    }

    cv.mean <- rowMeans(fold.loss, na.rm = TRUE)
    cv.se <- apply(fold.loss, 1L, stats::sd, na.rm = TRUE) /
        sqrt(rowSums(is.finite(fold.loss)))
    cv.mean[!is.finite(cv.mean)] <- Inf
    cv.se[!is.finite(cv.se)] <- Inf
    cv.table <- data.frame(
        lambda1 = grid$lambda1,
        lambda2 = grid$lambda2,
        cv.mean = cv.mean,
        cv.se = cv.se,
        n.successful.folds = rowSums(is.finite(fold.loss)),
        total.penalty = grid$lambda1 + grid$lambda2,
        stringsAsFactors = FALSE
    )
    cv.table <- cbind(cv.table, as.data.frame(fold.loss, check.names = FALSE))
    if (!any(is.finite(cv.table$cv.mean))) {
        stop("All SSRHE cross-validation fits failed.", call. = FALSE)
    }
    min.idx <- which.min(cv.table$cv.mean)
    if (identical(selection, "one.se")) {
        threshold <- cv.table$cv.mean[min.idx] + cv.table$cv.se[min.idx]
        eligible <- which(cv.table$cv.mean <= threshold)
        selected.idx <- eligible[which.max(cv.table$total.penalty[eligible])]
    } else {
        selected.idx <- min.idx
    }

    final <- .fit.ssrhe.hessian.from.operator(
        operator = operator,
        y = Y,
        lambda1 = cv.table$lambda1[selected.idx],
        lambda2 = cv.table$lambda2[selected.idx],
        weights = as.vector(W),
        ridge = ridge,
        verbose = verbose
    )
    final$X <- X
    final$support.selection <- "rule"
    final$selected.support <- data.frame(
        adaptive.k.scale = operator$parameters$adaptive.k.scale,
        min.support = operator$neighborhoods$min.support %||% NA_integer_,
        max.support = operator$neighborhoods$max.support %||% NA_integer_
    )
    final$cv.table <- cv.table
    final$fold.id <- fold.id
    final$selection <- list(
        rule = selection,
        loss = loss,
        selected.index = selected.idx,
        min.index = min.idx,
        lambda1 = cv.table$lambda1[selected.idx],
        lambda2 = cv.table$lambda2[selected.idx],
        cv.mean = cv.table$cv.mean[selected.idx],
        cv.se = cv.table$cv.se[selected.idx],
        nfolds = length(folds)
    )
    attr(final, "call") <- match.call()
    class(final) <- c("ssrhe.hessian.cv.fit", "ssrhe.hessian.fit", "list")
    final
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

.validate.ssrhe.lambda.grid <- function(x, name) {
    if (!is.numeric(x) || !length(x) ||
        any(is.na(x)) || any(!is.finite(x)) || any(x < 0)) {
        stop(sprintf("%s must be a nonempty finite nonnegative numeric vector.",
                     name), call. = FALSE)
    }
    sort(unique(as.double(x)))
}

.prepare.ssrhe.cv.folds <- function(observed, nfolds, fold.id) {
    n <- length(observed)
    if (is.null(fold.id)) {
        nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
        observed.idx <- which(observed)
        if (length(observed.idx) < 2L) {
            stop("At least two observed positive-weight labels are required for CV.",
                 call. = FALSE)
        }
        nfolds <- min(nfolds, length(observed.idx))
        out <- integer(n)
        out[observed.idx] <- rep(seq_len(nfolds), length.out = length(observed.idx))
    } else {
        if (length(fold.id) != n) {
            stop("fold.id must have length nrow(X).", call. = FALSE)
        }
        out <- as.integer(fold.id)
        out[is.na(out) | out < 1L | !observed] <- 0L
    }
    folds <- sort(unique(out[out > 0L]))
    if (length(folds) < 2L) {
        stop("At least two nonempty validation folds are required.", call. = FALSE)
    }
    for (ff in folds) {
        n.validation <- sum(out == ff)
        n.training <- sum(observed & out != ff)
        if (n.validation < 1L || n.training < 1L) {
            stop("Each validation fold must leave at least one training label.",
                 call. = FALSE)
        }
    }
    out
}

#' Fit SSRHE-Style Hessian L1 Regression
#'
#' Fits an \eqn{\ell_1}-adaptive SSRHE Hessian comparator using the row
#' operator \eqn{A} from \code{\link{ssrhe.hessian.operator}}:
#' \deqn{
#'   \widehat f_\lambda =
#'   \arg\min_f
#'   \left\{
#'     \frac{1}{2}\sum_i w_i(y_i-f_i)^2
#'     + \lambda \|Af\|_1
#'   \right\}.
#' }
#'
#' @inheritParams ssrhe.hessian.operator
#' @param y Numeric response vector of length \code{nrow(X)}. Matrix responses
#'   are not yet supported for the \eqn{\ell_1} path.
#' @param lambda.grid Optional nonnegative lambda grid. For
#'   \code{lambda.selection = "fixed"}, this must contain exactly one value. For
#'   \code{lambda.selection = "cv"}, \code{NULL} builds a default grid from the
#'   fitted full-data generalized-lasso path.
#' @param lambda.selection \code{"cv"} for observed-label K-fold
#'   cross-validation or \code{"fixed"} for a supplied fixed lambda.
#' @param weights Optional nonnegative observation weights. Missing responses
#'   automatically receive zero weight.
#' @param n.lambda Number of default lambda candidates when
#'   \code{lambda.grid = NULL} and \code{lambda.selection = "cv"}.
#' @param nfolds Number of validation folds over observed positive-weight
#'   labels. Ignored when \code{fold.id} is supplied.
#' @param fold.id Optional integer fold assignments of length \code{nrow(X)}.
#' @param loss Validation loss, currently \code{"mse"} or \code{"mae"}.
#' @param selection Selection rule. \code{"min"} chooses the smallest mean
#'   validation loss. \code{"one.se"} chooses the largest lambda within one
#'   standard error of the minimum.
#' @param support.selection Support-profile selection rule.
#'   \code{"rule"} uses the supplied \code{adaptive.k.scale},
#'   \code{min.support}, and \code{max.support}. \code{"cv"} is currently
#'   supported for \code{neighborhood.type = "adaptive.radius"} with
#'   \code{lambda.selection = "cv"} and chooses among rows of
#'   \code{support.grid} by an outer response cross-validation loop.
#' @param support.grid Optional data frame of support profiles with columns
#'   \code{adaptive.k.scale}, \code{min.support}, and optional
#'   \code{max.support}. If \code{NULL}, \code{\link{ssrhe.support.grid}} builds
#'   a compact default grid.
#' @param support.cv.max.candidates Maximum number of support profiles to try
#'   when \code{support.selection = "cv"}.
#' @param solver Solver backend. \code{"genlasso"} uses the generalized-lasso
#'   path backend. \code{"admm"} uses a fixed-lambda ADMM solver for each lambda
#'   value. \code{"auto"} tries \code{"genlasso"} first and falls back to ADMM
#'   when path extraction produces an error or non-finite fitted values.
#' @param row.scaling Optional row scaling for the penalty matrix before
#'   solving. \code{"l2"} scales nonzero rows of \eqn{A} to unit Euclidean norm.
#' @param admm.rho,admm.maxiter,admm.abstol,admm.reltol ADMM controls used when
#'   \code{solver = "admm"} or when \code{solver = "auto"} falls back to ADMM.
#' @param maxsteps,minlam,approx,rtol,btol,eps,verbose Controls passed to
#'   \pkg{genlasso}.
#'
#' @details
#' This function exposes the SSRHE local polynomial derivative estimator as an
#' \eqn{\ell_1} penalty, rather than as the original SSRHE \eqn{\ell_2}
#' Hessian-energy penalty used by
#' \code{\link{fit.ssrhe.hessian.regression}}. With
#' \code{derivative.order = 2L}, \eqn{A} contains local quadratic/Hessian rows.
#' With \code{derivative.order = 3L}, \eqn{A} contains local cubic
#' third-derivative tensor rows with the same symmetric component scaling used
#' by \code{\link{ssrhe.hessian.operator}}. The existing \eqn{\ell_2} fit
#' solves a linear system involving \eqn{B=A^\top A}. This function instead
#' keeps the rows of \eqn{A} and penalizes their absolute values. It is
#' therefore closer in spirit to graph trend filtering, while retaining SSRHE's
#' local-PCA derivative construction.
#'
#' The default solver backend is \pkg{genlasso}. For larger or numerically
#' fragile third-order operators, \code{solver = "admm"} fits the requested
#' lambda values directly without computing a generalized-lasso path, and
#' \code{solver = "auto"} falls back to ADMM when the path backend fails or
#' returns non-finite fitted values. Cross-validation removes held-out labels
#' from the data-fit term by setting their weights to zero, then scores
#' predictions on those held-out labels.
#'
#' With \code{support.selection = "cv"}, each adaptive-radius support candidate
#' builds a fresh SSRHE operator and runs the usual lambda CV with shared fold
#' assignments. The selected fit is the candidate with the smallest selected CV
#' error. This is useful when the default support-size rule is too rigid, but it
#' can be much slower than tuning lambda for one fixed operator.
#'
#' @return A list of class \code{"ssrhe.hessian.l1.fit"} containing fitted
#'   values, residuals, selected lambda, the SSRHE operator, generalized-lasso
#'   path metadata, lambda-grid fitted values, and CV diagnostics when
#'   requested.
#' @export
fit.ssrhe.hessian.l1.regression <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda.grid = NULL,
    lambda.selection = c("cv", "fixed"),
    weights = NULL,
    n.lambda = 40L,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    pinv.tol = sqrt(.Machine$double.eps),
    solver = c("genlasso", "admm", "auto"),
    row.scaling = c("none", "l2"),
    admm.rho = 1,
    admm.maxiter = 2000L,
    admm.abstol = 1e-4,
    admm.reltol = 1e-3,
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    support.selection = c("rule", "cv"),
    support.grid = NULL,
    support.cv.max.candidates = 8L,
    verbose = FALSE) {

    solver <- match.arg(solver)
    if (!identical(solver, "admm") &&
        !requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for fit.ssrhe.hessian.l1.regression().",
             call. = FALSE)
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for fit.ssrhe.hessian.l1.regression().",
             call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("fit.ssrhe.hessian.l1.regression currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda.selection <- match.arg(lambda.selection)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    support.selection <- match.arg(support.selection)
    neighborhood.type <- match.arg(neighborhood.type)
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    support.cv.max.candidates <- .validate.ssrhe.positive.integer(
        support.cv.max.candidates, "support.cv.max.candidates")
    solver.args <- .validate.ssrhe.hessian.l1.solver.args(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        observed = as.vector(y.info$observed & W > 0),
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        verbose = verbose
    )
    lambda.grid <- .validate.ssrhe.hessian.l1.lambda.grid(
        lambda.grid,
        lambda.selection
    )
    if (identical(support.selection, "cv")) {
        if (!identical(neighborhood.type, "adaptive.radius")) {
            stop("support.selection = 'cv' is currently supported only for neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        if (!identical(lambda.selection, "cv")) {
            stop("support.selection = 'cv' requires lambda.selection = 'cv'.",
                 call. = FALSE)
        }
        tangent.dim.grid <- if (missing(tangent.dim) || is.null(tangent.dim)) {
            ncol(X)
        } else {
            .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        }
        support.grid <- .validate.ssrhe.support.grid(
            support.grid = support.grid,
            n = n,
            tangent.dim = tangent.dim.grid,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = support.cv.max.candidates
        )
        fold.id.cv <- solver.args$fold.id %||% .prepare.ssrhe.cv.folds(
            observed = as.vector(y.info$observed & W > 0),
            nfolds = nfolds,
            fold.id = NULL
        )
        candidate.fits <- vector("list", nrow(support.grid))
        candidate.rows <- vector("list", nrow(support.grid))
        for (ii in seq_len(nrow(support.grid))) {
            max.support.ii <- support.grid$max.support[ii]
            if (is.na(max.support.ii)) max.support.ii <- NULL
            elapsed <- system.time({
                cand <- tryCatch(
                    fit.ssrhe.hessian.l1.regression(
                        X = X,
                        y = y,
                        k = k,
                        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
                        lambda.grid = lambda.grid,
                        lambda.selection = "cv",
                        weights = weights,
                        n.lambda = n.lambda,
                        nfolds = nfolds,
                        fold.id = fold.id.cv,
                        loss = loss,
                        selection = selection,
                        support.selection = "rule",
                        nn.index = nn.index,
                        neighborhood.type = "adaptive.radius",
                        support.index = support.index,
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        radius.rule = radius.rule,
                        radius.factor = radius.factor,
                        min.support = support.grid$min.support[ii],
                        max.support = max.support.ii,
                        support.buffer = support.buffer,
                        support.topup = support.topup,
                        tangent.dim.rule = tangent.dim.rule,
                        eigen.tolerance = eigen.tolerance,
                        derivative.order = derivative.order,
                        pinv.tol = pinv.tol,
                        solver = solver,
                        row.scaling = row.scaling,
                        admm.rho = admm.rho,
                        admm.maxiter = admm.maxiter,
                        admm.abstol = admm.abstol,
                        admm.reltol = admm.reltol,
                        maxsteps = maxsteps,
                        minlam = minlam,
                        approx = approx,
                        rtol = rtol,
                        btol = btol,
                        eps = eps,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            })[["elapsed"]]
            if (inherits(cand, "error")) {
                candidate.rows[[ii]] <- data.frame(
                    candidate = ii,
                    status = "error",
                    adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                    min.support = support.grid$min.support[ii],
                    max.support = support.grid$max.support[ii],
                    cv.error = Inf,
                    lambda = NA_real_,
                    runtime.sec = unname(elapsed),
                    message = conditionMessage(cand),
                    stringsAsFactors = FALSE
                )
            } else {
                candidate.fits[[ii]] <- cand
                diag <- .ssrhe.support.diagnostics(cand$operator)
                candidate.rows[[ii]] <- cbind(
                    data.frame(
                        candidate = ii,
                        status = "ok",
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        min.support = support.grid$min.support[ii],
                        max.support = support.grid$max.support[ii],
                        cv.error = cand$cv$error.selected,
                        cv.se = cand$cv$se[cand$cv$selected.idx],
                        lambda = cand$lambda,
                        runtime.sec = unname(elapsed),
                        message = "",
                        stringsAsFactors = FALSE
                    ),
                    diag
                )
            }
        }
        support.cv.table <- do.call(rbind, candidate.rows)
        if (!any(is.finite(support.cv.table$cv.error))) {
            stop("All SSRHE Hessian L1 support-CV candidates failed.",
                 call. = FALSE)
        }
        selected.support.idx <- which.min(support.cv.table$cv.error)
        out <- candidate.fits[[selected.support.idx]]
        out$support.selection <- "cv"
        out$support.cv.table <- support.cv.table
        out$selected.support <- support.grid[selected.support.idx, , drop = FALSE]
        out$selection$support.index <- selected.support.idx
        out$selection$support.cv.error <- support.cv.table$cv.error[selected.support.idx]
        attr(out, "call") <- match.call()
        class(out) <- c("ssrhe.hessian.l1.fit", "list")
        return(out)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = FALSE,
        pinv.tol = pinv.tol,
        return.A = TRUE,
        return.B = FALSE,
        return.BS = FALSE,
        return.sparse = TRUE,
        verbose = verbose
    )

    out <- .fit.ssrhe.hessian.l1.from.operator(
        operator = operator,
        y = as.vector(y.info$Y),
        weights = as.vector(W),
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        n.lambda = n.lambda,
        fold.id = solver.args$fold.id,
        loss = loss,
        selection = selection,
        solver.args = solver.args,
        verbose = verbose
    )
    out$X <- X
    out$support.selection <- "rule"
    out$selected.support <- data.frame(
        adaptive.k.scale = operator$parameters$adaptive.k.scale,
        min.support = operator$neighborhoods$min.support %||% NA_integer_,
        max.support = operator$neighborhoods$max.support %||% NA_integer_
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.l1.fit", "list")
    out
}

#' Refit SSRHE-Style Hessian L1 Regression
#'
#' Reuses the \eqn{A} operator from
#' \code{\link{fit.ssrhe.hessian.l1.regression}} to fit a new response or a new
#' lambda grid without rebuilding local neighborhoods or local PCA Hessian
#' rows.
#'
#' @param fitted.model A \code{"ssrhe.hessian.l1.fit"} object.
#' @param y.new Optional new numeric response vector. If \code{NULL}, the
#'   original response is reused.
#' @inheritParams fit.ssrhe.hessian.l1.regression
#'
#' @return A list of class \code{"ssrhe.hessian.l1.refit"}.
#' @export
refit.ssrhe.hessian.l1.regression <- function(
    fitted.model,
    y.new = NULL,
    lambda.grid = fitted.model$lambda,
    lambda.selection = c("fixed", "cv"),
    weights = NULL,
    n.lambda = 40L,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    solver = c("genlasso", "admm", "auto"),
    row.scaling = c("none", "l2"),
    admm.rho = 1,
    admm.maxiter = 2000L,
    admm.abstol = 1e-4,
    admm.reltol = 1e-3,
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    verbose = FALSE) {

    if (!inherits(fitted.model, "ssrhe.hessian.l1.fit")) {
        stop("fitted.model must be a 'ssrhe.hessian.l1.fit' object.",
             call. = FALSE)
    }
    if (is.null(fitted.model$operator) ||
        !inherits(fitted.model$operator, "ssrhe.hessian.operator")) {
        stop("fitted.model does not contain a reusable SSRHE operator.",
             call. = FALSE)
    }
    if (is.null(fitted.model$operator$A)) {
        stop("fitted.model operator does not contain A.", call. = FALSE)
    }
    if (is.null(y.new)) {
        y.new <- fitted.model$y
    }
    n <- length(fitted.model$fitted.values)
    y.info <- .prepare.ssrhe.response.matrix(y.new, n, "y.new")
    if (ncol(y.info$Y) != 1L) {
        stop("refit.ssrhe.hessian.l1.regression currently supports one response vector.",
             call. = FALSE)
    }
    if (is.null(weights)) {
        weights <- fitted.model$weights
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda.selection <- match.arg(lambda.selection)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    solver <- match.arg(solver)
    if (!identical(solver, "admm") &&
        !requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for the selected SSRHE Hessian L1 solver.",
             call. = FALSE)
    }
    solver.args <- .validate.ssrhe.hessian.l1.solver.args(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        observed = as.vector(y.info$observed & W > 0),
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        verbose = verbose
    )
    lambda.grid <- .validate.ssrhe.hessian.l1.lambda.grid(
        lambda.grid,
        lambda.selection
    )

    out <- .fit.ssrhe.hessian.l1.from.operator(
        operator = fitted.model$operator,
        y = as.vector(y.info$Y),
        weights = as.vector(W),
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        n.lambda = n.lambda,
        fold.id = solver.args$fold.id,
        loss = loss,
        selection = selection,
        solver.args = solver.args,
        verbose = verbose
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.l1.refit", "ssrhe.hessian.l1.fit", "list")
    out
}

.fit.ssrhe.hessian.l1.from.operator <- function(operator,
                                                y,
                                                weights,
                                                lambda.grid,
                                                lambda.selection,
                                                n.lambda,
                                                fold.id,
                                                loss,
                                                selection,
                                                solver.args,
                                                verbose = FALSE) {
    if (is.null(operator$A)) {
        stop("operator must contain A. Rebuild with return.A = TRUE.",
             call. = FALSE)
    }
    n <- ncol(operator$A)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("SSRHE Hessian L1 fitting currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    y.raw <- as.vector(y.info$Y)
    weights <- as.vector(W)
    y.clean <- y.raw
    y.clean[!is.finite(y.clean)] <- 0

    D.info <- .ssrhe.hessian.l1.penalty.matrix(
        operator$A,
        row.scaling = solver.args$row.scaling
    )
    diagnostics <- .ssrhe.hessian.l1.diagnostics(
        A = operator$A,
        D.info = D.info
    )
    backend <- solver.args$solver
    path.info <- NULL
    path <- NULL
    if (!identical(backend, "admm")) {
        path.info <- tryCatch(
            .fit.ssrhe.hessian.l1.path(
                y = y.clean,
                weights = weights,
                D = D.info$D,
                solver.args = solver.args,
                use.svd = D.info$svd
            ),
            error = function(e) {
                if (identical(backend, "auto")) {
                    structure(list(error = e), class = "ssrhe.l1.path.error")
                } else {
                    stop(e)
                }
            }
        )
        if (!inherits(path.info, "ssrhe.l1.path.error")) {
            path <- path.info$path
        }
        if (is.null(lambda.grid)) {
            if (inherits(path.info, "ssrhe.l1.path.error")) {
                stop("lambda.grid is required when solver = 'auto' falls back before a genlasso path is available.",
                     call. = FALSE)
            }
            lambda.grid <- .ssrhe.hessian.l1.default.lambda.grid(path, n.lambda)
        }
    } else if (is.null(lambda.grid)) {
        stop("lambda.grid is required when solver = 'admm'.", call. = FALSE)
    }
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.",
             call. = FALSE)
    }

    cv <- NULL
    if (identical(lambda.selection, "cv")) {
        if (is.null(fold.id)) {
            fold.id <- .prepare.ssrhe.cv.folds(
                observed = as.vector(y.info$observed & weights > 0),
                nfolds = solver.args$nfolds,
                fold.id = NULL
            )
        }
        cv <- .ssrhe.hessian.l1.cv(
            y = y.clean,
            weights = weights,
            D = D.info$D,
            lambda.grid = lambda.grid,
            fold.id = fold.id,
            loss = loss,
            selection = selection,
            solver.args = solver.args,
            use.svd = D.info$svd
        )
        selected.idx <- cv$selected.idx
    } else {
        selected.idx <- 1L
    }

    fit.source <- backend
    beta.grid <- NULL
    if (!identical(backend, "admm") &&
        !inherits(path.info, "ssrhe.l1.path.error")) {
        beta.grid <- .ssrhe.hessian.l1.coef.matrix(path, lambda.grid, n)
    }
    if (identical(backend, "admm") ||
        inherits(path.info, "ssrhe.l1.path.error") ||
        any(!is.finite(beta.grid))) {
        if (identical(backend, "genlasso")) {
            stop("SSRHE Hessian L1 path produced non-finite fitted values.",
                 call. = FALSE)
        }
        beta.grid <- .ssrhe.hessian.l1.admm.grid(
            y = y.clean,
            weights = weights,
            D = D.info$D.sparse,
            lambda.grid = lambda.grid,
            solver.args = solver.args
        )
        fit.source <- if (identical(backend, "auto")) "admm_fallback" else "admm"
    }
    fitted <- beta.grid[, selected.idx]
    observed <- as.vector(y.info$observed & weights > 0)
    residuals <- rep(NA_real_, n)
    residuals[observed] <- y.raw[observed] - fitted[observed]
    data.loss <- 0.5 * sum(weights * (y.clean - fitted)^2)
    hessian.l1 <- sum(abs(as.vector(operator$A %*% fitted)))
    lambda <- lambda.grid[selected.idx]

    list(
        fitted.values = as.vector(fitted),
        residuals = residuals,
        y = y.raw,
        weights = weights,
        lambda = lambda,
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        objective = data.loss + lambda * hessian.l1,
        energies = list(data.loss = data.loss, hessian.l1 = hessian.l1),
        operator = operator,
        beta.grid = beta.grid,
        path = path,
        cv = cv,
        fold.id = fold.id,
        observed = observed,
        n.responses = 1L,
        diagnostics = diagnostics,
        solver = list(
            backend = fit.source,
            requested = backend,
            representation = D.info$representation,
            svd = D.info$svd,
            row.scaling = D.info$row.scaling,
            n.observed = if (is.null(path.info) ||
                inherits(path.info, "ssrhe.l1.path.error")) {
                sum(observed)
            } else {
                path.info$n.observed
            },
            admm = if (grepl("admm", fit.source, fixed = TRUE)) {
                attr(beta.grid, "admm")
            } else {
                NULL
            }
        ),
        selection = list(
            rule = selection,
            loss = loss,
            selected.index = selected.idx,
            lambda = lambda
        )
    )
}

.ssrhe.hessian.l1.penalty.matrix <- function(A,
                                             row.scaling = c("none", "l2")) {
    row.scaling <- match.arg(row.scaling)
    A.sparse <- methods::as(A, "dgCMatrix")
    row.scale <- rep(1, nrow(A.sparse))
    if (identical(row.scaling, "l2")) {
        row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
        nz <- is.finite(row.norm) & row.norm > 0
        row.scale[nz] <- 1 / row.norm[nz]
        A.sparse <- Matrix::Diagonal(x = row.scale) %*% A.sparse
        A.sparse <- methods::as(A.sparse, "dgCMatrix")
    }
    use.svd <- nrow(A) >= ncol(A)
    list(
        D = if (use.svd) as.matrix(A.sparse) else A.sparse,
        D.sparse = A.sparse,
        svd = use.svd,
        representation = if (use.svd) "dense" else "sparse",
        row.scaling = row.scaling,
        row.scale = row.scale
    )
}

.ssrhe.hessian.l1.diagnostics <- function(A, D.info) {
    A.sparse <- methods::as(A, "dgCMatrix")
    scaled <- D.info$D.sparse
    row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
    scaled.row.norm <- sqrt(Matrix::rowSums(scaled^2))
    list(
        nrow = nrow(A.sparse),
        ncol = ncol(A.sparse),
        nnzero = Matrix::nnzero(A.sparse),
        density = Matrix::nnzero(A.sparse) / (nrow(A.sparse) * ncol(A.sparse)),
        row.norm = list(
            min = suppressWarnings(min(row.norm, na.rm = TRUE)),
            median = stats::median(row.norm),
            max = suppressWarnings(max(row.norm, na.rm = TRUE)),
            zero = sum(!is.finite(row.norm) | row.norm == 0)
        ),
        scaled.row.norm = list(
            min = suppressWarnings(min(scaled.row.norm, na.rm = TRUE)),
            median = stats::median(scaled.row.norm),
            max = suppressWarnings(max(scaled.row.norm, na.rm = TRUE)),
            zero = sum(!is.finite(scaled.row.norm) | scaled.row.norm == 0)
        ),
        row.scaling = D.info$row.scaling,
        representation = D.info$representation,
        svd = D.info$svd
    )
}

.fit.ssrhe.hessian.l1.path <- function(y, weights, D, solver.args,
                                       use.svd = FALSE) {
    observed <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(observed)) {
        stop("At least one observed positive-weight response is required.",
             call. = FALSE)
    }
    if (all(observed) && all(abs(weights - 1) < sqrt(.Machine$double.eps))) {
        path <- suppressWarnings(genlasso::genlasso(
            y = y,
            D = D,
            approx = solver.args$approx,
            maxsteps = solver.args$maxsteps,
            minlam = solver.args$minlam,
            rtol = solver.args$rtol,
            btol = solver.args$btol,
            eps = solver.args$eps,
            verbose = solver.args$verbose,
            svd = use.svd
        ))
    } else {
        obs <- which(observed)
        sw <- sqrt(weights[obs])
        X.design <- as.matrix(Matrix::sparseMatrix(
            i = seq_along(obs),
            j = obs,
            x = sw,
            dims = c(length(obs), length(y)),
            giveCsparse = TRUE
        ))
        path <- suppressWarnings(genlasso::genlasso(
            y = sw * y[obs],
            X = X.design,
            D = D,
            approx = solver.args$approx,
            maxsteps = solver.args$maxsteps,
            minlam = solver.args$minlam,
            rtol = solver.args$rtol,
            btol = solver.args$btol,
            eps = solver.args$eps,
            verbose = solver.args$verbose,
            svd = use.svd
        ))
    }
    list(path = path, n.observed = sum(observed))
}

.ssrhe.hessian.l1.default.lambda.grid <- function(path, n.lambda) {
    lambda.path <- path$lambda
    lambda.path <- lambda.path[is.finite(lambda.path) & lambda.path >= 0]
    lambda.max <- suppressWarnings(max(lambda.path, na.rm = TRUE))
    if (!is.finite(lambda.max) || lambda.max <= 0) {
        lambda.max <- 1
    }
    n.lambda <- max(2L, as.integer(n.lambda))
    lambda.positive <- lambda.path[lambda.path > 0]
    complete.path <- isTRUE(path$completepath)
    lambda.min <- if (complete.path || !length(lambda.positive)) {
        lambda.max * 1e-4
    } else {
        max(min(lambda.positive) * (1 + 1e-8), lambda.max * 1e-4)
    }
    grid <- exp(seq(log(lambda.max), log(lambda.min), length.out = n.lambda - 1L))
    if (complete.path) grid <- c(grid, 0)
    unique(grid)
}

.ssrhe.hessian.l1.coef.matrix <- function(path, lambda.grid, n) {
    beta <- tryCatch(
        as.matrix(stats::coef(path, lambda = lambda.grid)$beta),
        error = function(e) NULL
    )
    if (is.null(beta) ||
        !identical(dim(beta), c(as.integer(n), length(lambda.grid)))) {
        beta <- vapply(lambda.grid, function(lambda) {
            tryCatch(
                as.vector(stats::coef(path, lambda = lambda)$beta),
                error = function(e) rep(NA_real_, n)
            )
        }, numeric(n))
    }
    beta
}

.ssrhe.soft.threshold <- function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
}

.ssrhe.hessian.l1.admm <- function(y, weights, D, lambda, solver.args) {
    D <- methods::as(D, "dgCMatrix")
    n <- length(y)
    m <- nrow(D)
    weights <- as.numeric(weights)
    y <- as.numeric(y)
    rho <- solver.args$admm.rho
    AtA <- Matrix::crossprod(D)
    system <- AtA * rho + Matrix::Diagonal(n, x = weights)
    ridge <- max(1e-10, sqrt(.Machine$double.eps) * max(1, mean(Matrix::diag(system))))
    system <- system + Matrix::Diagonal(n, x = rep(ridge, n))
    factor <- Matrix::Cholesky(system, LDL = FALSE, Imult = 0)
    wy <- weights * y
    beta <- rep(0, n)
    z <- rep(0, m)
    u <- rep(0, m)
    converged <- FALSE
    iter <- 0L
    primal <- dual <- Inf
    for (iter in seq_len(solver.args$admm.maxiter)) {
        z.old <- z
        rhs <- wy + rho * as.vector(Matrix::crossprod(D, z - u))
        beta <- as.vector(Matrix::solve(factor, rhs))
        Dbeta <- as.vector(D %*% beta)
        z <- .ssrhe.soft.threshold(Dbeta + u, lambda / rho)
        u <- u + Dbeta - z
        primal <- sqrt(sum((Dbeta - z)^2))
        dual <- rho * sqrt(sum(as.vector(Matrix::crossprod(D, z - z.old))^2))
        eps.pri <- sqrt(m) * solver.args$admm.abstol +
            solver.args$admm.reltol * max(sqrt(sum(Dbeta^2)), sqrt(sum(z^2)))
        eps.dual <- sqrt(n) * solver.args$admm.abstol +
            solver.args$admm.reltol *
            sqrt(sum(as.vector(rho * Matrix::crossprod(D, u))^2))
        if (primal <= eps.pri && dual <= eps.dual) {
            converged <- TRUE
            break
        }
    }
    attr(beta, "admm") <- list(
        iterations = iter,
        converged = converged,
        primal.residual = primal,
        dual.residual = dual,
        rho = rho
    )
    beta
}

.ssrhe.hessian.l1.admm.grid <- function(y, weights, D, lambda.grid,
                                        solver.args) {
    admm.info <- vector("list", length(lambda.grid))
    beta <- vapply(seq_along(lambda.grid), function(ii) {
        fit <- .ssrhe.hessian.l1.admm(
            y = y,
            weights = weights,
            D = D,
            lambda = lambda.grid[ii],
            solver.args = solver.args
        )
        admm.info[[ii]] <<- attr(fit, "admm")
        as.vector(fit)
    }, numeric(length(y)))
    attr(beta, "admm") <- admm.info
    beta
}

.ssrhe.hessian.l1.cv <- function(y, weights, D, lambda.grid, fold.id, loss,
                                 selection, solver.args, use.svd = FALSE) {
    n <- length(y)
    folds <- sort(unique(fold.id[fold.id > 0L]))
    cv.errors <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda.grid))
    rownames(cv.errors) <- paste0("Fold", folds)
    colnames(cv.errors) <- format(signif(lambda.grid, 6), scientific = TRUE)

    for (ii in seq_along(folds)) {
        fold <- folds[ii]
        test <- which(fold.id == fold)
        train.weights <- weights
        train.weights[test] <- 0
        if (identical(solver.args$solver, "admm")) {
            beta <- tryCatch(
                .ssrhe.hessian.l1.admm.grid(
                    y = y,
                    weights = train.weights,
                    D = D,
                    lambda.grid = lambda.grid,
                    solver.args = solver.args
                ),
                error = function(e) NULL
            )
        } else {
            path.info <- tryCatch(
                .fit.ssrhe.hessian.l1.path(
                    y = y,
                    weights = train.weights,
                    D = D,
                    solver.args = solver.args,
                    use.svd = use.svd
                ),
                error = function(e) NULL
            )
            if (is.null(path.info)) {
                beta <- NULL
            } else {
                beta <- .ssrhe.hessian.l1.coef.matrix(path.info$path, lambda.grid, n)
            }
            if (identical(solver.args$solver, "auto") &&
                (is.null(beta) || any(!is.finite(beta)))) {
                beta <- tryCatch(
                    .ssrhe.hessian.l1.admm.grid(
                        y = y,
                        weights = train.weights,
                        D = methods::as(D, "dgCMatrix"),
                        lambda.grid = lambda.grid,
                        solver.args = solver.args
                    ),
                    error = function(e) NULL
                )
            }
        }
        if (is.null(beta)) {
            next
        }
        beta.test <- beta[test, , drop = FALSE]
        finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
        target <- matrix(y[test], nrow = length(test), ncol = length(lambda.grid))
        w.test <- matrix(weights[test], nrow = length(test), ncol = length(lambda.grid))
        if (identical(loss, "mse")) {
            fold.error <- colSums(w.test * (target - beta.test)^2) / colSums(w.test)
        } else {
            fold.error <- colSums(w.test * abs(target - beta.test)) / colSums(w.test)
        }
        fold.error[!finite.pred] <- NA_real_
        cv.errors[ii, ] <- fold.error
    }

    mean.error <- colMeans(cv.errors, na.rm = TRUE)
    mean.error[!is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    best.idx <- which.min(mean.error)
    if (!is.finite(mean.error[best.idx])) {
        stop("All SSRHE Hessian L1 CV fits failed for the supplied lambda grid.",
             call. = FALSE)
    }
    if (identical(selection, "one.se")) {
        threshold <- mean.error[best.idx] + se[best.idx]
        eligible <- which(mean.error <= threshold)
        selected.idx <- eligible[which.max(lambda.grid[eligible])]
    } else {
        selected.idx <- best.idx
    }
    list(
        fold.id = fold.id,
        lambda.grid = lambda.grid,
        fold.errors = cv.errors,
        mean.error = mean.error,
        se = se,
        best.idx = best.idx,
        selected.idx = selected.idx,
        lambda.min = lambda.grid[best.idx],
        lambda = lambda.grid[selected.idx],
        error.min = mean.error[best.idx],
        error.selected = mean.error[selected.idx],
        selection = selection,
        loss = loss
    )
}

.validate.ssrhe.hessian.l1.lambda.grid <- function(lambda.grid,
                                                   lambda.selection) {
    if (is.null(lambda.grid)) {
        if (identical(lambda.selection, "fixed")) {
            stop("lambda.grid is required when lambda.selection = 'fixed'.",
                 call. = FALSE)
        }
        return(NULL)
    }
    lambda.grid <- .validate.ssrhe.lambda.grid(lambda.grid, "lambda.grid")
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.",
             call. = FALSE)
    }
    lambda.grid
}

.validate.ssrhe.hessian.l1.solver.args <- function(n.lambda,
                                                   nfolds,
                                                   fold.id,
                                                   observed,
                                                   maxsteps,
                                                   minlam,
                                                   approx,
                                                   rtol,
                                                   btol,
                                                   eps,
                                                   solver = c("genlasso", "admm", "auto"),
                                                   row.scaling = c("none", "l2"),
                                                   admm.rho = 1,
                                                   admm.maxiter = 2000L,
                                                   admm.abstol = 1e-4,
                                                   admm.reltol = 1e-3,
                                                   verbose) {
    n.lambda <- .validate.ssrhe.positive.integer(n.lambda, "n.lambda")
    nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
    maxsteps <- .validate.ssrhe.positive.integer(maxsteps, "maxsteps")
    solver <- match.arg(solver)
    row.scaling <- match.arg(row.scaling)
    minlam <- .validate.ssrhe.nonnegative.scalar(minlam, "minlam")
    rtol <- .validate.ssrhe.nonnegative.scalar(rtol, "rtol")
    btol <- .validate.ssrhe.nonnegative.scalar(btol, "btol")
    eps <- .validate.ssrhe.nonnegative.scalar(eps, "eps")
    admm.rho <- .validate.ssrhe.numeric.scalar(admm.rho, "admm.rho")
    if (!is.finite(admm.rho) || admm.rho <= 0) {
        stop("admm.rho must be a finite positive numeric scalar.",
             call. = FALSE)
    }
    admm.maxiter <- .validate.ssrhe.positive.integer(admm.maxiter, "admm.maxiter")
    admm.abstol <- .validate.ssrhe.nonnegative.scalar(admm.abstol, "admm.abstol")
    admm.reltol <- .validate.ssrhe.nonnegative.scalar(admm.reltol, "admm.reltol")
    if (!is.null(fold.id)) {
        fold.id <- .prepare.ssrhe.cv.folds(observed, nfolds, fold.id)
    }
    list(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        maxsteps = maxsteps,
        minlam = minlam,
        approx = isTRUE(approx),
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        verbose = isTRUE(verbose)
    )
}

#' @export
print.ssrhe.hessian.l1.fit <- function(x, ...) {
    cat("SSRHE Hessian L1 regression fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda.selection:", x$lambda.selection, "\n")
    cat("  lambda:", format(x$lambda, digits = 4), "\n")
    cat("  objective:", format(x$objective, digits = 4), "\n")
    if (!is.null(x$cv)) {
        cat("  CV", x$cv$loss, ":",
            format(x$cv$error.selected, digits = 4), "\n")
    }
    invisible(x)
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

#' @export
print.ssrhe.hessian.cv.fit <- function(x, ...) {
    cat("SSRHE Hessian regression CV fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  selected lambda1:", format(x$selection$lambda1, digits = 4), "\n")
    cat("  selected lambda2:", format(x$selection$lambda2, digits = 4), "\n")
    cat("  CV", x$selection$loss, ":",
        format(x$selection$cv.mean, digits = 4), "\n")
    invisible(x)
}

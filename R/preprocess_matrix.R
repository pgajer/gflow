#' Preprocess a Feature Matrix by Log-Shift, Winsorization, and Robust Z-Scoring
#'
#' @description
#' Applies a per-feature log transform with an additive pseudocount (log-shift),
#' followed by winsorization and robust z-scoring, to produce a matrix suitable
#' for geometry-based graph construction (e.g., kNN or ikNN graphs).
#'
#' @details
#' The preprocessing steps are applied column-wise:
#' \enumerate{
#'   \item \code{S} is coerced to a numeric matrix; non-finite values are replaced
#'         with \code{0}.
#'   \item Negative values are clamped to \code{0} (useful for assays reporting
#'         background-subtracted signals such as net MFI).
#'   \item A pseudocount \code{pc[j]} is added per feature and a natural log is taken:
#'         \code{log(x + pc[j])}.
#'   \item Values are winsorized using \code{winsorize.vec()}.
#'   \item Values are robustly standardized using \code{robust.zscore.vec()}.
#' }
#'
#' The default pseudocount rule \code{"half.min.pos"} sets \code{pc[j]} to one half
#' of the minimum strictly positive finite value in column \code{j}; if a column
#' has no positive finite values, \code{pc[j]} defaults to \code{1}.
#'
#' @param S Numeric matrix (or coercible to a numeric matrix) with samples in rows
#'   and features in columns.
#' @param winsor.probs Numeric vector of length 2 giving lower and upper
#'   winsorization probabilities passed to \code{winsorize.vec()}.
#' @param pseudocount Pseudocount specification. Either \code{"half.min.pos"} for
#'   a per-feature pseudocount as described, or a single positive numeric value
#'   used for all features.
#'
#' @return A numeric matrix of the same dimensions as \code{S} containing the
#'   transformed values.
#'
#' @seealso \code{\link{preprocess.matrix.asinh}},
#'   \code{\link{preprocess.matrix.normscores}},
#'   \code{\link{select.features.by.floor.mass}}
#'
#' @examples
#' ## Toy example with zeros and negatives
#' set.seed(1)
#' S <- matrix(rnorm(20), nrow = 5)
#' S[1, 1] <- -3
#' S[2, 2] <- 0
#' X <- preprocess.matrix.logshift(S)
#' dim(X)
#'
#' @export
preprocess.matrix.logshift <- function(S,
                                       winsor.probs = c(0.01, 0.99),
                                       pseudocount = "half.min.pos") {
    ## ---- validate inputs (CRAN-safe defensive checks) ----
    if (missing(S) || is.null(S)) stop("`S` must be a non-null matrix-like object.")
    X <- tryCatch(as.matrix(S), error = function(e) NULL)
    if (is.null(X)) stop("`S` must be coercible to a matrix via `as.matrix()`.")
    if (!is.numeric(X)) {
        suppressWarnings(storage.mode(X) <- "double")
        if (!is.numeric(X)) stop("`S` must be numeric (or coercible to numeric).")
    }
    if (length(dim(X)) != 2L) stop("`S` must be a 2D matrix.")
    if (ncol(X) < 1L) stop("`S` must have at least one column.")
    if (nrow(X) < 1L) stop("`S` must have at least one row.")

    if (!is.numeric(winsor.probs) || length(winsor.probs) != 2L ||
        any(!is.finite(winsor.probs))) {
        stop("`winsor.probs` must be a finite numeric vector of length 2.")
    }
    if (winsor.probs[1] < 0 || winsor.probs[2] > 1 || winsor.probs[1] >= winsor.probs[2]) {
        stop("`winsor.probs` must satisfy 0 <= p1 < p2 <= 1.")
    }

    if (is.character(pseudocount)) {
        if (length(pseudocount) != 1L || is.na(pseudocount) ||
            pseudocount != "half.min.pos") {
            stop("`pseudocount` character value must be \"half.min.pos\".")
        }
    } else if (is.numeric(pseudocount)) {
        if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0) {
            stop("`pseudocount` numeric value must be a single finite number > 0.")
        }
    } else {
        stop("`pseudocount` must be either \"half.min.pos\" or a single numeric > 0.")
    }

    ## ---- core transform ----
    X[!is.finite(X)] <- 0
    X <- pmax(X, 0)  ## clamp negative net MFI to 0

    ## per-feature pseudocounts
    pc <- rep(1.0, ncol(X))
    if (!is.null(colnames(X))) names(pc) <- colnames(X)

    if (is.character(pseudocount) && pseudocount == "half.min.pos") {
        for (j in seq_len(ncol(X))) {
            v <- X[, j]
            vpos <- v[v > 0 & is.finite(v)]
            pc[j] <- if (length(vpos) > 0L) 0.5 * min(vpos) else 1.0
            ## safety: enforce strictly positive pseudocount
            if (!is.finite(pc[j]) || pc[j] <= 0) pc[j] <- 1.0
        }
    } else {
        pc[] <- pseudocount
    }

    for (j in seq_len(ncol(X))) X[, j] <- log(X[, j] + pc[j])
    for (j in seq_len(ncol(X))) X[, j] <- winsorize.vec(X[, j], probs = winsor.probs)
    for (j in seq_len(ncol(X))) X[, j] <- robust.zscore.vec(X[, j])

    X
}

#' Preprocess a Feature Matrix by Scaled asinh, Winsorization, and Robust Z-Scoring
#'
#' @description
#' Applies a per-feature scaled inverse hyperbolic sine (\code{asinh}) transform,
#' followed by winsorization and robust z-scoring, to produce a matrix suitable
#' for geometry-based graph construction (e.g., kNN or ikNN graphs).
#'
#' @details
#' The preprocessing steps are applied column-wise:
#' \enumerate{
#'   \item \code{S} is coerced to a numeric matrix; non-finite values are replaced
#'         with \code{0}.
#'   \item A per-feature scale \code{aa[j]} is selected:
#'         \itemize{
#'           \item \code{a = "median.pos"}: \code{aa[j]} is the median of strictly
#'                 positive finite values in column \code{j}; if none exist, \code{aa[j] = 1}.
#'           \item \code{a} numeric scalar: the same positive scale is used for all features.
#'         }
#'   \item The scaled transform is applied: \code{asinh(x / aa[j])}.
#'         (Unlike log transforms, \code{asinh} can handle negative values.)
#'   \item Values are winsorized using \code{winsorize.vec()}.
#'   \item Values are robustly standardized using \code{robust.zscore.vec()}.
#' }
#'
#' @param S Numeric matrix (or coercible to a numeric matrix) with samples in rows
#'   and features in columns.
#' @param a Scale specification for the \code{asinh} transform. Either \code{"median.pos"}
#'   for per-feature median-positive scaling, or a single positive numeric value used
#'   for all features.
#' @param winsor.probs Numeric vector of length 2 giving lower and upper
#'   winsorization probabilities passed to \code{winsorize.vec()}.
#'
#' @return A numeric matrix of the same dimensions as \code{S} containing the
#'   transformed values.
#'
#' @seealso \code{\link{preprocess.matrix.logshift}},
#'   \code{\link{preprocess.matrix.normscores}},
#'   \code{\link{select.features.by.floor.mass}}
#'
#' @examples
#' ## Toy example with negatives, zeros, and non-finite values
#' set.seed(1)
#' S <- matrix(rnorm(20), nrow = 5)
#' S[1, 1] <- -10
#' S[2, 2] <- 0
#' S[3, 3] <- NA_real_
#' X <- preprocess.matrix.asinh(S)
#' dim(X)
#'
#' @export
preprocess.matrix.asinh <- function(S,
                                    a = "median.pos",
                                    winsor.probs = c(0.01, 0.99)) {
    ## ---- validate inputs (CRAN-safe defensive checks) ----
    if (missing(S) || is.null(S)) stop("`S` must be a non-null matrix-like object.")
    X <- tryCatch(as.matrix(S), error = function(e) NULL)
    if (is.null(X)) stop("`S` must be coercible to a matrix via `as.matrix()`.")
    if (!is.numeric(X)) {
        suppressWarnings(storage.mode(X) <- "double")
        if (!is.numeric(X)) stop("`S` must be numeric (or coercible to numeric).")
    }
    if (length(dim(X)) != 2L) stop("`S` must be a 2D matrix.")
    if (ncol(X) < 1L) stop("`S` must have at least one column.")
    if (nrow(X) < 1L) stop("`S` must have at least one row.")

    if (!is.numeric(winsor.probs) || length(winsor.probs) != 2L ||
        any(!is.finite(winsor.probs))) {
        stop("`winsor.probs` must be a finite numeric vector of length 2.")
    }
    if (winsor.probs[1] < 0 || winsor.probs[2] > 1 || winsor.probs[1] >= winsor.probs[2]) {
        stop("`winsor.probs` must satisfy 0 <= p1 < p2 <= 1.")
    }

    if (is.character(a)) {
        if (length(a) != 1L || is.na(a) || a != "median.pos") {
            stop("`a` character value must be \"median.pos\".")
        }
    } else if (is.numeric(a)) {
        if (length(a) != 1L || !is.finite(a) || a <= 0) {
            stop("`a` numeric value must be a single finite number > 0.")
        }
    } else {
        stop("`a` must be either \"median.pos\" or a single numeric > 0.")
    }

    ## ---- core transform ----
    X[!is.finite(X)] <- 0

    ## choose scale per feature
    aa <- rep(1.0, ncol(X))
    if (!is.null(colnames(X))) names(aa) <- colnames(X)

    if (is.character(a) && a == "median.pos") {
        for (j in seq_len(ncol(X))) {
            v <- X[, j]
            vpos <- v[v > 0 & is.finite(v)]
            aa[j] <- if (length(vpos) > 0L) stats::median(vpos) else 1.0
            if (!is.finite(aa[j]) || aa[j] <= 0) aa[j] <- 1.0
        }
    } else {
        aa[] <- a
    }

    ## safety: avoid division by zero or non-finite scale
    aa[!is.finite(aa) | aa <= 0] <- 1.0

    for (j in seq_len(ncol(X))) X[, j] <- asinh(X[, j] / aa[j])
    for (j in seq_len(ncol(X))) X[, j] <- winsorize.vec(X[, j], probs = winsor.probs)
    for (j in seq_len(ncol(X))) X[, j] <- robust.zscore.vec(X[, j])

    X
}

#' Preprocess a Feature Matrix by Normal Scores, Winsorization, and Robust Z-Scoring
#'
#' @description
#' Produces a geometry-friendly feature matrix by applying, column-wise:
#' winsorization in the original scale, rank-based normal scores
#' (approximately standard normal marginal distributions), and robust z-scoring.
#'
#' @details
#' The preprocessing steps are applied column-wise:
#' \enumerate{
#'   \item \code{S} is coerced to a numeric matrix; non-finite values are replaced
#'         with \code{0}.
#'   \item Optionally, negative values are clamped to \code{0} (useful for assays
#'         with background-subtracted signals).
#'   \item Values are winsorized via \code{winsorize.vec()} to limit extreme tails
#'         prior to ranking.
#'   \item Rank-based normal scores are computed via \code{normal.scores.vec()}.
#'   \item Values are robustly standardized using \code{robust.zscore.vec()}.
#' }
#'
#' This transform is often effective when features have heterogeneous scales,
#' heavy tails, or outliers, and when the downstream geometry should be less
#' sensitive to marginal distributions.
#'
#' @param S Numeric matrix (or coercible to a numeric matrix) with samples in rows
#'   and features in columns.
#' @param winsor.probs Numeric vector of length 2 giving lower and upper
#'   winsorization probabilities passed to \code{winsorize.vec()}.
#' @param clamp.neg.to.zero Logical; if \code{TRUE}, negative entries are set to \code{0}
#'   before winsorization and ranking.
#'
#' @return A numeric matrix of the same dimensions as \code{S} containing the
#'   transformed values.
#'
#' @seealso \code{\link{normal.scores.vec}}, \code{\link{winsorize.vec}},
#'   \code{\link{robust.zscore.vec}}, \code{\link{preprocess.matrix.logshift}},
#'   \code{\link{preprocess.matrix.asinh}}
#'
#' @examples
#' set.seed(1)
#' S <- matrix(rnorm(20), nrow = 5)
#' S[1, 1] <- -5
#' S[2, 2] <- 100
#' X <- preprocess.matrix.normscores(S)
#' dim(X)
#'
#' @export
preprocess.matrix.normscores <- function(S,
                                        winsor.probs = c(0.01, 0.99),
                                        clamp.neg.to.zero = TRUE) {
    ## ---- validate inputs (CRAN-safe defensive checks) ----
    if (missing(S) || is.null(S)) stop("`S` must be a non-null matrix-like object.")
    X <- tryCatch(as.matrix(S), error = function(e) NULL)
    if (is.null(X)) stop("`S` must be coercible to a matrix via `as.matrix()`.")
    if (!is.numeric(X)) {
        suppressWarnings(storage.mode(X) <- "double")
        if (!is.numeric(X)) stop("`S` must be numeric (or coercible to numeric).")
    }
    if (length(dim(X)) != 2L) stop("`S` must be a 2D matrix.")
    if (ncol(X) < 1L) stop("`S` must have at least one column.")
    if (nrow(X) < 1L) stop("`S` must have at least one row.")

    if (!is.numeric(winsor.probs) || length(winsor.probs) != 2L ||
        any(!is.finite(winsor.probs))) {
        stop("`winsor.probs` must be a finite numeric vector of length 2.")
    }
    if (winsor.probs[1] < 0 || winsor.probs[2] > 1 || winsor.probs[1] >= winsor.probs[2]) {
        stop("`winsor.probs` must satisfy 0 <= p1 < p2 <= 1.")
    }

    if (!is.logical(clamp.neg.to.zero) || length(clamp.neg.to.zero) != 1L ||
        is.na(clamp.neg.to.zero)) {
        stop("`clamp.neg.to.zero` must be a single non-missing logical value.")
    }

    ## ---- core transform ----
    X[!is.finite(X)] <- 0
    if (clamp.neg.to.zero) X <- pmax(X, 0)

    ## winsorize in raw space before ranking
    for (j in seq_len(ncol(X))) X[, j] <- winsorize.vec(X[, j], probs = winsor.probs)

    for (j in seq_len(ncol(X))) X[, j] <- normal.scores.vec(X[, j])

    ## after normal scores, robust scale is usually unnecessary, but safe
    for (j in seq_len(ncol(X))) X[, j] <- robust.zscore.vec(X[, j])

    X
}


#' Winsorize a Numeric Vector
#'
#' @description
#' Winsorizes a numeric vector by clamping values below/above specified
#' quantiles to the corresponding quantile values.
#'
#' @param x Numeric vector.
#' @param probs Numeric vector of length 2 giving lower and upper probabilities.
#'   Must satisfy \code{0 <= p1 < p2 <= 1}.
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @examples
#' x <- c(-100, rnorm(100), 100)
#' y <- winsorize.vec(x, probs = c(0.01, 0.99))
#' range(x); range(y)
#'
#' @export
winsorize.vec <- function(x, probs = c(0.01, 0.99)) {
    ## ---- validate inputs ----
    if (!is.numeric(x)) {
        suppressWarnings(x <- as.numeric(x))
        if (!is.numeric(x)) stop("`x` must be numeric (or coercible to numeric).")
    }
    if (!is.numeric(probs) || length(probs) != 2L || any(!is.finite(probs))) {
        stop("`probs` must be a finite numeric vector of length 2.")
    }
    if (probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
        stop("`probs` must satisfy 0 <= p1 < p2 <= 1.")
    }

    if (length(x) == 0L) return(x)
    if (all(!is.finite(x))) return(x)

    xf <- x[is.finite(x)]
    if (length(xf) == 0L) return(x)

    qs <- stats::quantile(xf, probs = probs, na.rm = TRUE, type = 7)
    x[x < qs[1]] <- qs[1]
    x[x > qs[2]] <- qs[2]
    x
}


#' Robust Z-Score Standardization for a Numeric Vector
#'
#' @description
#' Standardizes a numeric vector using a robust location and scale estimator:
#' the median for location and MAD (median absolute deviation) for scale.
#' If MAD is too small or non-finite, falls back to the standard deviation.
#' If both scale estimates are degenerate, returns a vector of zeros.
#'
#' @param x Numeric vector.
#' @param eps Small positive threshold used to detect degenerate scale estimates.
#'
#' @return A numeric vector of the same length as \code{x}, standardized to
#'   approximately zero median and unit robust scale.
#'
#' @examples
#' x <- c(rep(1, 10), 100)
#' z <- robust.zscore.vec(x)
#' median(z); stats::mad(z)
#'
#' @export
robust.zscore.vec <- function(x, eps = 1e-12) {
    ## ---- validate inputs ----
    if (!is.numeric(x)) {
        suppressWarnings(x <- as.numeric(x))
        if (!is.numeric(x)) stop("`x` must be numeric (or coercible to numeric).")
    }
    if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
        stop("`eps` must be a single finite numeric value > 0.")
    }

    if (length(x) == 0L) return(x)

    m <- stats::median(x, na.rm = TRUE)
    s <- stats::mad(x, center = m, constant = 1.4826, na.rm = TRUE)

    if (!is.finite(s) || s < eps) s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s < eps) return(x * 0)

    (x - m) / s
}

#' Select Features by Floor-Mass Criterion
#'
#' @description
#' Computes, for each feature (column), the fraction of finite values that are
#' equal (within tolerance) to that feature's minimum (the "floor mass").
#' Returns a logical vector indicating which features pass a maximum floor-mass
#' threshold. This is useful for identifying features dominated by a detection
#' floor / censoring (many values piled at the minimum), which can distort
#' distance-based geometry and kNN/ikNN graph construction.
#'
#' @param X Numeric matrix (or coercible to a numeric matrix) with samples in rows
#'   and features in columns.
#' @param floor.mass.max Numeric scalar in \eqn{[0,1]}; features with floor mass
#'   greater than this value are flagged for removal.
#' @param tol Numeric scalar > 0. Tolerance for deciding whether a value is equal
#'   to the column minimum.
#'
#' @return A list with components:
#' \describe{
#'   \item{keep}{Logical vector of length \code{ncol(X)} indicating which features
#'     satisfy \code{floor.mass <= floor.mass.max}. Features with undefined floor
#'     mass (no finite values) are set to \code{FALSE}.}
#'   \item{floor.mass}{Named numeric vector of floor-mass values for each feature.}
#' }
#'
#' @seealso \code{\link{floor_mass}}, \code{\link{preprocess.matrix.logshift}},
#'   \code{\link{preprocess.matrix.asinh}}, \code{\link{preprocess.matrix.normscores}}
#'
#' @examples
#' X <- matrix(c(0,0,0,1,  2,2,2,2), nrow = 4)
#' colnames(X) <- c("floor.heavy", "constant")
#' select.features.by.floor.mass(X, floor.mass.max = 0.6)
#'
#' @rawNamespace export(select.features.by.floor.mass)
select.features.by.floor.mass <- function(X,
                                         floor.mass.max = 0.60,
                                         tol = 1e-12) {
    ## ---- validate inputs (CRAN-safe defensive checks) ----
    if (missing(X) || is.null(X)) stop("`X` must be a non-null matrix-like object.")
    M <- tryCatch(as.matrix(X), error = function(e) NULL)
    if (is.null(M)) stop("`X` must be coercible to a matrix via `as.matrix()`.")
    if (!is.numeric(M)) {
        suppressWarnings(storage.mode(M) <- "double")
        if (!is.numeric(M)) stop("`X` must be numeric (or coercible to numeric).")
    }
    if (length(dim(M)) != 2L) stop("`X` must be a 2D matrix.")
    if (nrow(M) < 1L) stop("`X` must have at least one row.")
    if (ncol(M) < 1L) stop("`X` must have at least one column.")

    if (!is.numeric(floor.mass.max) || length(floor.mass.max) != 1L ||
        !is.finite(floor.mass.max)) {
        stop("`floor.mass.max` must be a single finite numeric value in [0, 1].")
    }
    if (floor.mass.max < 0 || floor.mass.max > 1) {
        stop("`floor.mass.max` must lie in [0, 1].")
    }

    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("`tol` must be a single finite numeric value > 0.")
    }

    ## ---- compute floor mass per feature ----
    fm <- vapply(seq_len(ncol(M)),
                 function(j) floor_mass(M[, j], tol = tol),
                 numeric(1))

    if (!is.null(colnames(M))) names(fm) <- colnames(M)

    keep <- fm <= floor.mass.max
    keep[is.na(keep)] <- FALSE

    list(keep = keep, floor.mass = fm)
}


#' Floor Mass of a Numeric Vector
#'
#' @description
#' Computes the fraction of finite values that are equal (within tolerance) to
#' the minimum finite value of the vector. This quantifies the extent to which
#' observations pile up at a lower detection bound or floor.
#'
#' @param x Numeric vector.
#' @param tol Numeric scalar > 0. Tolerance for deciding whether a value is equal
#'   to the minimum.
#'
#' @return A numeric scalar in \eqn{[0,1]} giving the fraction of finite values within
#'   \code{tol} of the minimum finite value. Returns \code{NA_real_} if \code{x}
#'   contains no finite values.
#'
#' @examples
#' floor_mass(c(0, 0, 0, 1, 2))
#' floor_mass(c(NA, NaN, Inf, -Inf))
#'
#' @export
floor_mass <- function(x, tol = 1e-12) {
    ## ---- validate inputs ----
    if (!is.numeric(x)) {
        suppressWarnings(x <- as.numeric(x))
        if (!is.numeric(x)) stop("`x` must be numeric (or coercible to numeric).")
    }
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("`tol` must be a single finite numeric value > 0.")
    }

    xf <- x[is.finite(x)]
    if (length(xf) == 0L) return(NA_real_)

    mn <- min(xf)
    mean(abs(xf - mn) <= tol)
}


#' Rank-Based Normal Scores for a Numeric Vector
#'
#' @description
#' Converts a numeric vector to rank-based normal scores (Blom-style),
#' using average ranks for ties and mapping ranks to quantiles of the standard
#' normal distribution. Non-finite values are returned as \code{NA_real_}.
#'
#' @details
#' For finite entries, ranks \code{r} are mapped to \code{u = (r - 0.5) / n} and
#' then transformed via \code{stats::qnorm(u)}. If there are fewer than two finite
#' values, the function returns a vector of zeros (matching the input length).
#'
#' @param x Numeric vector (or coercible to numeric).
#'
#' @return A numeric vector of the same length as \code{x} containing normal scores
#'   for finite entries and \code{NA_real_} for non-finite entries.
#'
#' @examples
#' normal.scores.vec(c(1, 2, 2, 10, NA, Inf))
#'
#' @export
normal.scores.vec <- function(x) {
    ## ---- validate inputs ----
    if (!is.numeric(x)) {
        suppressWarnings(x <- as.numeric(x))
        if (!is.numeric(x)) stop("`x` must be numeric (or coercible to numeric).")
    }

    ok <- is.finite(x)
    n <- sum(ok)

    if (n <= 1L) return(x * 0)

    r <- rep(NA_real_, length(x))
    r[ok] <- rank(x[ok], ties.method = "average")

    u <- (r[ok] - 0.5) / n
    z <- stats::qnorm(u)

    out <- rep(NA_real_, length(x))
    out[ok] <- z
    out[!ok] <- NA_real_
    out
}

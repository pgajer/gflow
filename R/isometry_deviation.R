.validate.distance.matrix <- function(D, name) {
    if (inherits(D, "dist")) {
        D <- as.matrix(D)
    }
    if (!(is.matrix(D) || is.data.frame(D))) {
        stop(sprintf("'%s' must be a matrix, data frame, or dist object.", name),
             call. = FALSE)
    }
    D <- as.matrix(D)
    if (!is.numeric(D) || nrow(D) != ncol(D)) {
        stop(sprintf("'%s' must be a square numeric distance matrix.", name),
             call. = FALSE)
    }
    if (any(!is.finite(D))) {
        stop(sprintf("'%s' cannot contain NA, NaN, or Inf values.", name),
             call. = FALSE)
    }
    if (!is.double(D)) {
        storage.mode(D) <- "double"
    }
    D
}

.isometry.pair.values <- function(D.estimated, D.true, true.tol = sqrt(.Machine$double.eps)) {
    D.estimated <- .validate.distance.matrix(D.estimated, "D.estimated")
    D.true <- .validate.distance.matrix(D.true, "D.true")
    if (!identical(dim(D.estimated), dim(D.true))) {
        stop("'D.estimated' and 'D.true' must have the same dimensions.", call. = FALSE)
    }
    if (!is.numeric(true.tol) || length(true.tol) != 1L || !is.finite(true.tol) ||
        true.tol < 0) {
        stop("'true.tol' must be a finite non-negative numeric scalar.", call. = FALSE)
    }
    idx <- upper.tri(D.true, diag = FALSE)
    estimated <- as.numeric(D.estimated[idx])
    true <- as.numeric(D.true[idx])
    keep <- is.finite(estimated) & is.finite(true) & true > true.tol
    if (!any(keep)) {
        stop("No usable off-diagonal distance pairs remain after filtering.",
             call. = FALSE)
    }
    list(estimated = estimated[keep], true = true[keep])
}

#' Compute the Optimal Isometry Calibration Scale
#'
#' @description
#' Computes the least-squares scalar \eqn{\alpha} that aligns estimated
#' distances to true distances:
#' \deqn{\alpha^* = \frac{\sum_{i<j} \hat D_{ij}D_{ij}}
#'                  {\sum_{i<j} \hat D_{ij}^2}.}
#'
#' @param D.estimated Estimated distance matrix or `dist` object.
#' @param D.true Ground-truth distance matrix or `dist` object.
#' @param true.tol Non-negative tolerance. Pairs with true distance less than
#'   or equal to this value are excluded from relative metrics.
#'
#' @return Numeric scalar calibration scale.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' isometry.scale(2 * D, D)
#'
#' @export
isometry.scale <- function(D.estimated, D.true, true.tol = sqrt(.Machine$double.eps)) {
    pairs <- .isometry.pair.values(D.estimated, D.true, true.tol)
    denom <- sum(pairs$estimated^2)
    if (denom <= 0) {
        stop("Estimated distances have zero norm; cannot compute calibration scale.",
             call. = FALSE)
    }
    sum(pairs$estimated * pairs$true) / denom
}

#' Compute Relative RMS Isometry Error
#'
#' @description
#' Computes the relative root-mean-square distance-matrix error after optional
#' scalar calibration:
#' \deqn{
#' \left[
#' \frac{\sum_{i<j}(\alpha\hat D_{ij}-D_{ij})^2}
#'      {\sum_{i<j}D_{ij}^2}
#' \right]^{1/2}.
#' }
#'
#' @inheritParams isometry.scale
#' @param scale Logical scalar. If `TRUE`, use `isometry.scale()` to optimally
#'   rescale estimated distances before computing the error.
#'
#' @return Numeric scalar relative RMS error.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' isometry.rel.rms.error(2 * D, D)
#'
#' @export
isometry.rel.rms.error <- function(D.estimated,
                                   D.true,
                                   scale = TRUE,
                                   true.tol = sqrt(.Machine$double.eps)) {
    pairs <- .isometry.pair.values(D.estimated, D.true, true.tol)
    alpha <- if (isTRUE(scale)) {
        isometry.scale(D.estimated, D.true, true.tol)
    } else {
        1
    }
    sqrt(sum((alpha * pairs$estimated - pairs$true)^2) / sum(pairs$true^2))
}

#' Compute Relative Absolute Isometry Errors
#'
#' @description
#' Computes quantiles of pairwise relative absolute errors
#' \deqn{\left|\alpha\hat D_{ij}-D_{ij}\right|/D_{ij}.}
#'
#' @inheritParams isometry.rel.rms.error
#' @param probs Numeric vector of probabilities passed to `stats::quantile()`.
#'
#' @return Named numeric vector of relative absolute error quantiles.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' isometry.rel.abs.error(1.1 * D, D)
#'
#' @export
isometry.rel.abs.error <- function(D.estimated,
                                   D.true,
                                   probs = c(0.5, 0.95),
                                   scale = TRUE,
                                   true.tol = sqrt(.Machine$double.eps)) {
    pairs <- .isometry.pair.values(D.estimated, D.true, true.tol)
    alpha <- if (isTRUE(scale)) {
        isometry.scale(D.estimated, D.true, true.tol)
    } else {
        1
    }
    errors <- abs(alpha * pairs$estimated - pairs$true) / pairs$true
    stats::quantile(errors, probs = probs, names = TRUE, type = 7)
}

#' Compute Multiplicative Distortion Quantiles
#'
#' @description
#' Computes quantiles of the calibrated pairwise distortion
#' \deqn{\rho_{ij}=\alpha\hat D_{ij}/D_{ij}.}
#' Values near 1 indicate approximately isometric distance preservation.
#'
#' @inheritParams isometry.rel.rms.error
#' @param probs Numeric vector of probabilities passed to `stats::quantile()`.
#'
#' @return Named numeric vector of distortion quantiles.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' isometry.distortion.quantiles(1.1 * D, D)
#'
#' @export
isometry.distortion.quantiles <- function(D.estimated,
                                          D.true,
                                          probs = c(0.05, 0.5, 0.95),
                                          scale = TRUE,
                                          true.tol = sqrt(.Machine$double.eps)) {
    pairs <- .isometry.pair.values(D.estimated, D.true, true.tol)
    alpha <- if (isTRUE(scale)) {
        isometry.scale(D.estimated, D.true, true.tol)
    } else {
        1
    }
    distortion <- alpha * pairs$estimated / pairs$true
    stats::quantile(distortion, probs = probs, names = TRUE, type = 7)
}

#' Compute Distance Preservation Correlations
#'
#' @description
#' Computes Pearson and Spearman correlations between estimated and true
#' pairwise distances over the upper-triangular distance entries.
#'
#' @inheritParams isometry.scale
#'
#' @return Named numeric vector with entries `pearson_cor` and `spearman_cor`.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' isometry.distance.correlations(2 * D, D)
#'
#' @export
isometry.distance.correlations <- function(D.estimated,
                                           D.true,
                                           true.tol = sqrt(.Machine$double.eps)) {
    pairs <- .isometry.pair.values(D.estimated, D.true, true.tol)
    if (length(pairs$true) < 2L || stats::sd(pairs$true) == 0 ||
        stats::sd(pairs$estimated) == 0) {
        return(c(pearson_cor = NA_real_, spearman_cor = NA_real_))
    }
    c(
        pearson_cor = stats::cor(pairs$estimated, pairs$true, method = "pearson"),
        spearman_cor = stats::cor(pairs$estimated, pairs$true, method = "spearman")
    )
}

#' Summarize Deviation from Isometry
#'
#' @description
#' Computes the standard benchmark summary used to compare graph geodesic
#' distances with reference geodesic distances: optimal scalar calibration,
#' relative RMS error, relative absolute error quantiles, multiplicative
#' distortion quantiles, and distance correlations.
#'
#' @inheritParams isometry.rel.rms.error
#'
#' @return A one-row data frame.
#'
#' @examples
#' D <- as.matrix(dist(1:4))
#' summarize.isometry.deviation(1.1 * D, D)
#'
#' @export
summarize.isometry.deviation <- function(D.estimated,
                                         D.true,
                                         scale = TRUE,
                                         true.tol = sqrt(.Machine$double.eps)) {
    alpha <- if (isTRUE(scale)) {
        isometry.scale(D.estimated, D.true, true.tol)
    } else {
        1
    }
    rel.abs <- isometry.rel.abs.error(
        D.estimated, D.true, probs = c(0.5, 0.95),
        scale = scale, true.tol = true.tol
    )
    distortion <- isometry.distortion.quantiles(
        D.estimated, D.true, probs = c(0.05, 0.5, 0.95),
        scale = scale, true.tol = true.tol
    )
    cors <- isometry.distance.correlations(D.estimated, D.true, true.tol)
    data.frame(
        scale = alpha,
        rel_rms_error = isometry.rel.rms.error(D.estimated, D.true, scale, true.tol),
        rel_abs_error_median = unname(rel.abs[[1L]]),
        rel_abs_error_q95 = unname(rel.abs[[2L]]),
        distortion_q05 = unname(distortion[[1L]]),
        distortion_median = unname(distortion[[2L]]),
        distortion_q95 = unname(distortion[[3L]]),
        pearson_cor = unname(cors[["pearson_cor"]]),
        spearman_cor = unname(cors[["spearman_cor"]])
    )
}

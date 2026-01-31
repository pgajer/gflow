#' Vertex-wise permutation test for rdgraph regression fitted values
#'
#' @description
#' Performs a vertex-wise permutation test for a fitted rdgraph regression model.
#' The response vector \code{y} is permuted across vertices \code{n.perms} times,
#' the model is refit using \code{refit.rdgraph.regression()}, and a p-value is
#' computed at each vertex by comparing the observed statistic to the permutation
#' null distribution.
#'
#' The default statistic is \eqn{T_v = |\hat y(v) - \bar y|}, where \eqn{\bar y}
#' is the global mean of the original \code{y}. This is appropriate for binary
#' outcomes (incidence) and yields a two-sided test against departures from the
#' global mean incidence at each vertex.
#'
#' @param fitted.model A fitted model object returned by
#'   \code{\link{fit.rdgraph.regression}}.
#'   Must contain \code{fitted.values}.
#'
#' @param y Numeric vector of length n (original response used to produce
#'   \code{fitted.model}). For binary incidence, use 0/1 (numeric or integer).
#'
#' @param n.perms Integer \eqn{\ge 1}. Number of permutations.
#'
#' @param per.column.gcv Logical. Passed to \code{\link{refit.rdgraph.regression}}.
#'   If TRUE, reselect eta by GCV for each permutation column.
#'
#' @param n.cores Integer. Passed to \code{\link{refit.rdgraph.regression}}.
#'
#' @param seed Integer or NULL. If not NULL, sets RNG seed for reproducibility.
#'
#' @param stat Character scalar specifying the global test statistic to compute
#'   in addition to the vertex-wise p-values. Must be one of \code{"l2"},
#'   \code{"log"}, or \code{"logit"}.
#'
#'   The global statistic is computed from the fitted field \eqn{\hat y} and
#'   the global mean incidence \eqn{\bar y}:
#'   \itemize{
#'     \item \code{"l2"}: \eqn{T_2(\hat y) = \sum_{v\in V} (\hat y(v) - \bar y)^2}
#'     \item \code{"log"}: \eqn{T_{2,\log}(\hat y) = \sum_{v\in V} (\log(\hat y_\epsilon(v)) - \log(\bar y_\epsilon))^2}
#'     \item \code{"logit"}: \eqn{T_{2,\mathrm{logit}}(\hat y) = \sum_{v\in V} (\mathrm{logit}(\hat y_\epsilon(v)) - \mathrm{logit}(\bar y_\epsilon))^2}
#'   }
#'
#'   For \code{"log"} and \code{"logit"}, fitted values are clipped to
#'   \eqn{[\epsilon, 1-\epsilon]} with \eqn{\epsilon = \max(10^{-6}, 0.5/n)}
#'   to avoid infinities when \eqn{\hat y(v)} is numerically 0 or 1.
#'
#'   The global permutation p-value is computed as
#'   \eqn{p = (1 + \#\{b: T^{(b)} \ge T^{\mathrm{obs}}\})/(B+1)}.
#'
#' @param two.sided Logical. If TRUE (default), use absolute-value statistic.
#'   If FALSE, uses one-sided statistic \eqn{\hat y(v) - \bar y} (upper-tail).
#'
#' @param return.perm.fits Logical. If TRUE, return the matrix of permuted
#'   fitted values (can be large). Default FALSE.
#'
#' @param verbose Logical. If TRUE, print progress messages. Default TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{p.value}}{Numeric vector of length n, vertex-wise permutation p-values.}
#'     \item{\code{q.value}}{BH-adjusted p-values (FDR), same length as \code{p.value}.}
#'     \item{\code{stat.obs}}{Observed test statistic at each vertex.}
#'     \item{\code{stat.perm}}{Optional matrix of permutation statistics (n x n.perms) if requested.}
#'     \item{\code{y.mean}}{Global mean \eqn{\bar y} used for centering.}
#'     \item{\code{fitted.obs}}{Observed fitted values \eqn{\hat y}.}
#'     \item{\code{perm.refit}}{The refit object returned by \code{refit.rdgraph.regression}.}
#'     \item{\code{perm.fitted.values}}{Optional matrix of permuted fitted values if requested.}
#'     \item{\code{call}}{Matched call.}
#'   }
#'
#' @examples
#' \dontrun{
#' y.fit <- fit.rdgraph.regression(X0.smoothed, y = as.double(sptb.int.im), k = k.opt)
#' perm.res <- permutation.test.rdgraph(
#'   fitted.model = y.fit,
#'   y = as.double(sptb.int.im),
#'   n.perms = 200,
#'   per.column.gcv = TRUE,
#'   n.cores = 14,
#'   seed = 1
#' )
#' summary(perm.res$p.value)
#' }
#'
#' @export
permutation.test.rdgraph <- function(fitted.model,
                                     y,
                                     n.perms,
                                     per.column.gcv = TRUE,
                                     n.cores = 1L,
                                     seed = 12345L,
                                     stat = c("l2","log","logit"),
                                     two.sided = TRUE,
                                     return.perm.fits = FALSE,
                                     verbose = TRUE) {
    ## ---------------------------------------------------------------------
    ## Input validation
    ## ---------------------------------------------------------------------
    if (is.null(fitted.model)) stop("fitted.model cannot be NULL.")
    if (is.null(fitted.model$fitted.values)) {
        stop("fitted.model must contain fitted.values.")
    }

    if (!is.numeric(y)) stop("y must be numeric.")
    y <- as.double(y)

    n <- length(y)
    if (n < 2) stop("y must have length >= 2.")

    fitted.obs <- fitted.model$fitted.values
    if (!is.numeric(fitted.obs)) stop("fitted.model$fitted.values must be numeric.")
    if (length(fitted.obs) != n) {
        stop("Length mismatch: length(y) != length(fitted.model$fitted.values).")
    }

    if (!is.numeric(n.perms) || length(n.perms) != 1 || n.perms < 1 || n.perms != floor(n.perms)) {
        stop("n.perms must be a positive integer.")
    }
    n.perms <- as.integer(n.perms)

    if (!is.logical(per.column.gcv) || length(per.column.gcv) != 1) {
        stop("per.column.gcv must be TRUE/FALSE.")
    }

    if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1 || n.cores != floor(n.cores)) {
        stop("n.cores must be a positive integer.")
    }
    n.cores <- as.integer(n.cores)

    if (!is.logical(two.sided) || length(two.sided) != 1) {
        stop("two.sided must be TRUE/FALSE.")
    }
    if (!is.logical(return.perm.fits) || length(return.perm.fits) != 1) {
        stop("return.perm.fits must be TRUE/FALSE.")
    }
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be TRUE/FALSE.")
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed)) {
            stop("seed must be an integer or NULL.")
        }
        seed <- as.integer(seed)
        set.seed(seed)
    }

    stat <- match.arg(stat)

    ## ---------------------------------------------------------------------
    ## Build permutation matrix (n x n.perms)
    ## ---------------------------------------------------------------------
    if (verbose) {
        message(sprintf("Generating %d permutations for n = %d vertices...", n.perms, n))
    }

    perm.idx <- replicate(n.perms, sample.int(n, size = n, replace = FALSE))
    perm.y <- matrix(y[perm.idx], nrow = n, ncol = n.perms)

    ## ---------------------------------------------------------------------
    ## Refit using cached spectral structure
    ## ---------------------------------------------------------------------
    if (verbose) {
        message("Refitting permuted responses with refit.rdgraph.regression()...")
    }

    perm.refit <- refit.rdgraph.regression(
        fitted.model = fitted.model,
        y.new = perm.y,
        per.column.gcv = per.column.gcv,
        n.cores = n.cores,
        verbose = verbose
    )

    perm.fitted <- perm.refit$fitted.values
    if (is.null(perm.fitted)) stop("refit.rdgraph.regression returned NULL fitted.values.")
    if (!is.matrix(perm.fitted) || nrow(perm.fitted) != n || ncol(perm.fitted) != n.perms) {
        stop("Unexpected shape of permuted fitted.values; expected n x n.perms matrix.")
    }

    ## ---------------------------------------------------------------------
    ## Compute vertex-wise permutation p-values
    ## ---------------------------------------------------------------------
    y.mean <- mean(y)

    ## ---------------------------------------------------------------------
    ## Global permutation test(s): T2 on probability / log / logit scale
    ## ---------------------------------------------------------------------
    ## Note: We compute one global statistic chosen by `stat`, using the same
    ##       permuted fitted fields already available in `perm.fitted`.
    ##       This is independent of the vertex-wise statistic controlled by `two.sided`.
    eps <- max(1e-6, 0.5 / n)

    ## helper: clip to [eps, 1-eps]
    clip01 <- function(x, eps) pmin(1 - eps, pmax(eps, x))

    y.mean.eps <- clip01(y.mean, eps)

    if (stat == "l2") {
        ## T2 = sum (yhat - ybar)^2
        global.stat.obs <- sum((fitted.obs - y.mean)^2)
        ## perm.fitted is n x n.perms: compute column-wise sum of squares
        global.stat.perm <- colSums((perm.fitted - y.mean)^2)
    } else if (stat == "log") {
        ## T2log = sum (log(yhat_eps) - log(ybar_eps))^2
        fitted.obs.eps <- clip01(fitted.obs, eps)
        perm.fitted.eps <- clip01(perm.fitted, eps)

        global.stat.obs <- sum((log(fitted.obs.eps) - log(y.mean.eps))^2)

        ## log(yhat_eps) is n x n.perms, subtract scalar log(ybar_eps)
        global.stat.perm <- colSums((log(perm.fitted.eps) - log(y.mean.eps))^2)
    } else if (stat == "logit") {
        ## T2logit = sum (logit(yhat_eps) - logit(ybar_eps))^2
        logit <- function(p) log(p / (1 - p))

        fitted.obs.eps <- clip01(fitted.obs, eps)
        perm.fitted.eps <- clip01(perm.fitted, eps)

        global.stat.obs <- sum((logit(fitted.obs.eps) - logit(y.mean.eps))^2)

        global.stat.perm <- colSums((logit(perm.fitted.eps) - logit(y.mean.eps))^2)
    }

    ## Upper-tail p-value with +1 correction
    global.p.value <- (sum(global.stat.perm >= global.stat.obs) + 1) / (n.perms + 1)

    if (two.sided) {
        stat.obs <- abs(fitted.obs - y.mean)
        stat.perm <- abs(perm.fitted - y.mean)
        ## upper tail (two-sided via abs)
        ge.count <- rowSums(stat.perm >= stat.obs)
    } else {
        stat.obs <- fitted.obs - y.mean
        stat.perm <- perm.fitted - y.mean
        ## upper tail only
        ge.count <- rowSums(stat.perm >= stat.obs)
    }

    ## Phipson-Smyth style +1 correction for exactness
    p.value <- (ge.count + 1) / (n.perms + 1)

    ## Multiple testing correction (BH/FDR)
    q.value <- stats::p.adjust(p.value, method = "BH")

    out <- list(
        p.value = p.value,
        q.value = q.value,
        stat.obs = stat.obs,
        y.mean = y.mean,
        stat = stat,
        global.stat.obs = global.stat.obs,
        global.p.value = global.p.value,
        fitted.obs = fitted.obs,
        perm.refit = perm.refit,
        eps = eps,
        call = match.call()
    )

    if (return.perm.fits) {
        out$perm.fitted.values <- perm.fitted
        out$stat.perm <- stat.perm
        out$global.stat.perm <- global.stat.perm
    }

    class(out) <- c("rdgraph_permutation_test", "list")
    return(out)
}

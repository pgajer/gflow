#' Box-Cox power transform
#'
#' Applies the Box-Cox transform to a positive numeric vector \eqn{y}.
#' Uses \eqn{\log(y)} when \eqn{\lambda \approx 0}, and
#' \eqn{(y^\lambda - 1)/\lambda} otherwise.
#'
#' @param y A numeric vector of strictly positive and finite values.
#' @param lambda A numeric scalar giving the power parameter \eqn{\lambda}.
#'
#' @details
#' The transform is defined as
#' \deqn{
#'   T_\lambda(y) =
#'   \begin{cases}
#'     \dfrac{y^\lambda - 1}{\lambda}, & \lambda \neq 0, \\
#'     \log(y),                         & \lambda = 0.
#'   \end{cases}
#' }
#' This implementation switches to the log form when \eqn{|\lambda| < 10^{-8}}.
#'
#' @return A numeric vector of the same length as \code{y}, containing
#' the transformed values.
#'
#' @examples
#' y <- rexp(50, rate = 2)          # positive data
#' boxcox.transform(y, lambda = 0)  # log-transform
#' boxcox.transform(y, lambda = 0.5)
#'
#' @seealso \code{\link{boxcox.mle}}
#' @rawNamespace export(boxcox.transform)
boxcox.transform <- function(y, lambda) {
  if (any(y <= 0 | !is.finite(y))) {
    stop("Box-Cox requires y > 0 and finite.")
  }
  if (abs(lambda) < 1e-8) {
    return(log(y))
  } else {
    return((y^lambda - 1) / lambda)
  }
}

#' @title Profile log-likelihood for the Box-Cox transform (internal)
#' @description Computes the profile log-likelihood \eqn{\ell(\lambda)} for a
#'   Gaussian linear model fitted to \eqn{T_\lambda(y)}.
#'
#' @param lambda Numeric scalar \eqn{\lambda}.
#' @param y Numeric response vector with \eqn{y_i > 0}.
#' @param Xmat Model matrix (including intercept if desired).
#'
#' @details
#' Up to an additive constant, the profile log-likelihood is
#' \deqn{
#'   \ell(\lambda) = -\frac{n}{2}\log\!\left(\frac{\mathrm{RSS}(\lambda)}{n}\right)
#'                   + (\lambda - 1)\sum_{i=1}^n \log y_i,
#' }
#' where \eqn{\mathrm{RSS}(\lambda)} is the residual sum of squares from
#' \code{lm.fit(X, T_\lambda(y))}.
#'
#' @return A numeric scalar giving \eqn{\ell(\lambda)} (up to a constant).
#'
#' @keywords internal
#' @noRd
.boxcox.loglik <- function(lambda, y, Xmat) {
  yt <- boxcox.transform(y, lambda)
  fit <- lm.fit(x = Xmat, y = yt)             # fast core of lm()
  rss <- sum(fit$residuals^2)
  n <- length(y)
  # constant terms drop out; this is enough for argmax comparisons
  -(n/2) * log(rss / n) + (lambda - 1) * sum(log(y))
}

#' Box-Cox MLE for the power parameter \eqn{\lambda}
#'
#' Estimates the Box-Cox power parameter \eqn{\lambda} by profiling the
#' Gaussian log-likelihood over a grid with optional local refinement.
#'
#' @param formula A model formula of the form \code{y ~ x1 + x2 + ...}, or
#'   an \code{\link{lm}} object. The response \code{y} must be strictly
#'   positive.
#' @param data Optional data frame for \code{formula}. Ignored if \code{formula}
#'   is an \code{lm} object.
#' @param lambdas Numeric vector of candidate \eqn{\lambda} values used to build
#'   the coarse profile (default \code{seq(-2, 2, by = 0.1)}).
#' @param refine Logical; if \code{TRUE} (default) performs a 1D \code{\link{optimize}}
#'   search near the best grid point to refine \eqn{\lambda}.
#'
#' @details
#' For each \eqn{\lambda} in \code{lambdas}, the response is transformed via
#' \code{\link{boxcox.transform}} and a least-squares fit is computed with
#' \code{\link{lm.fit}}. The profiled log-likelihood
#' \eqn{\ell(\lambda)} is
#' \deqn{
#'   \ell(\lambda) = -\frac{n}{2}\log\!\left(\frac{\mathrm{RSS}(\lambda)}{n}\right)
#'                   + (\lambda - 1)\sum_{i=1}^n \log y_i,
#' }
#' up to an additive constant. A 95\% confidence interval is obtained from the
#' likelihood-ratio cutoff \eqn{2\{\ell(\hat\lambda)-\ell(\lambda)\} \le \chi^2_{1,0.95}}.
#'
#' The implementation covers the standard unweighted Gaussian linear model.
#' (Weighted fits can be added using \code{lm.wfit} and a weighted
#' \eqn{\sum \log y} term.)
#'
#' @return A list with components:
#' \item{lambda}{The MLE \eqn{\hat\lambda}.}
#' \item{loglik}{A data frame with columns \code{lambda} and \code{loglik}
#'   giving the profile evaluated on the coarse and dense grids.}
#' \item{ci95}{Approximate 95\% confidence interval for \eqn{\lambda}
#'   from the LR cutoff, returned as a numeric vector of length 2.}
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- runif(n)
#' y <- exp(1 + 2 * x + rnorm(n, sd = 0.2))  # positive response
#' d <- data.frame(y = y, x = x)
#'
#' fit <- boxcox.mle(y ~ x, data = d)
#' fit$lambda
#' head(fit$loglik)
#' fit$ci95
#'
#' # Transform y with the estimated lambda
#' y.bc <- boxcox.transform(d$y, fit$lambda)
#'
#' @references
#' Box, G. E. P. and Cox, D. R. (1964).
#' An analysis of transformations. \emph{Journal of the Royal Statistical Society. Series B}, \strong{26}(2), 211-252.
#'
#' @rawNamespace export(boxcox.mle)
boxcox.mle <- function(formula,
                       data,
                       lambdas = seq(-2, 2, by = 0.1),
                       refine = TRUE) {
  mf <- model.frame(formula, data = data)
  y  <- model.response(mf)
  if (any(y <= 0 | !is.finite(y))) {
    stop("Box-Cox requires positive, finite response. Consider Yeo-Johnson for y<=0.")
  }
  Xmat <- model.matrix(attr(mf, "terms"), data = mf)

  # Name the vapply args explicitly to avoid mis-matching
  ll <- vapply(X = lambdas,
               FUN = .boxcox.loglik,
               FUN.VALUE = numeric(1),
               y = y, Xmat = Xmat)

  # coarse argmax
  i_max <- which.max(ll)
  lambda_hat <- lambdas[i_max]

  # optional 1D refinement near the best grid point
  if (isTRUE(refine)) {
    step <- if (length(lambdas) > 1) max(diff(lambdas)) else 0.25
    lo <- max(min(lambdas), lambda_hat - step)
    hi <- min(max(lambdas), lambda_hat + step)
    opt <- optimize(function(l) -.boxcox.loglik(l, y, Xmat), interval = c(lo, hi))
    lambda_hat <- opt$minimum
  }

  # approximate 95% CI via LR cutoff
  ll_hat   <- .boxcox.loglik(lambda_hat, y, Xmat)
  cutoff   <- ll_hat - 0.5 * qchisq(0.95, df = 1)
  dense    <- seq(max(min(lambdas), lambda_hat - 2),
                  min(max(lambdas), lambda_hat + 2),
                  length.out = 400)
  ll_dense <- vapply(X = dense,
                     FUN = .boxcox.loglik,
                     FUN.VALUE = numeric(1),
                     y = y, Xmat = Xmat)
  ok <- which(ll_dense >= cutoff)
  ci <- if (length(ok)) range(dense[ok]) else c(NA_real_, NA_real_)

  list(lambda = lambda_hat,
       loglik = data.frame(lambda = c(lambdas, dense),
                           loglik  = c(ll,       ll_dense)),
       ci95 = ci)
}

# Private numerical helper retained only for legacy internal association tests.
# It is deliberately not exported: generic p-value aggregation is not a stable
# gflow estimand.
weighted.p.value <- function(u,
                             mu,
                             sigma,
                             alternative = c("two.sided", "less", "greater")) {
    if (!is.numeric(u) || !length(u) || any(!is.finite(u))) {
        stop("'u' must be a non-empty finite numeric vector")
    }
    if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) {
        stop("'mu' must be a single finite numeric value")
    }
    if (!is.numeric(sigma) || length(sigma) != 1L ||
        !is.finite(sigma) || sigma <= 0) {
        stop("'sigma' must be a single positive finite numeric value")
    }

    alternative <- match.arg(alternative)
    p.values <- switch(
        alternative,
        less = stats::pnorm(u, mean = mu, sd = sigma),
        greater = stats::pnorm(u, mean = mu, sd = sigma, lower.tail = FALSE),
        two.sided = {
            lower <- stats::pnorm(u, mean = mu, sd = sigma)
            upper <- stats::pnorm(u, mean = mu, sd = sigma,
                                  lower.tail = FALSE)
            2 * pmin(lower, upper)
        }
    )

    pmax(0, pmin(1, mean(p.values)))
}

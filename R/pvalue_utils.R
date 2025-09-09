#' Weighted P-value Calculation
#'
#' Computes the weighted p-value for a sample from a distribution quantifying
#' uncertainty of estimation of a value for which the p-value is to be computed,
#' given a normal null distribution with specified mean and standard deviation.
#' The weighted p-value takes into account the uncertainty in the value for
#' which the classical p-value would normally be computed.
#'
#' @param u     A numeric vector representing the sample drawn from a distribution
#'              quantifying uncertainty of estimation of a value for which the
#'              p-value is to be computed.
#' @param mu    The mean of the normal null distribution.
#' @param sigma The standard deviation of the normal null distribution.
#' @param alternative Character string specifying the alternative hypothesis,
#'                    must be one of "two.sided" (default), "greater" or "less".
#'
#' @return A single numeric value representing the weighted p-value.
#'
#' @details
#' The weighted p-value approach is useful when there is uncertainty in the
#' test statistic itself. Instead of computing a single p-value based on a
#' point estimate, this method:
#' \enumerate{
#'   \item Computes p-values for each value in the uncertainty distribution
#'   \item Averages these p-values to obtain a weighted p-value
#' }
#'
#' This approach is particularly valuable in:
#' \itemize{
#'   \item Bayesian contexts where posterior distributions represent uncertainty
#'   \item Measurement error scenarios
#'   \item Bootstrap or simulation-based inference
#' }
#'
#' The weighted p-value can be interpreted as the expected p-value over the
#' uncertainty distribution of the parameter.
#'
#' @references
#' For theoretical foundations of weighted p-values in the context of
#' uncertainty quantification, see relevant literature on Bayesian p-values
#' and posterior predictive checks.
#'
#' @examples
#' # Example 1: Uncertainty represented by a normal distribution
#' u_sample <- rnorm(1000, mean = 1.5, sd = 0.5)
#'
#' # Test against null hypothesis N(0, 1)
#' p_value <- weighted.p.value(u_sample, mu = 0, sigma = 1)
#' cat("Weighted p-value:", p_value, "\n")
#'
#' # Compare with classical p-value using the mean
#' classical_p <- pnorm(mean(u_sample), mean = 0, sd = 1, lower.tail = FALSE)
#' cat("Classical p-value:", classical_p, "\n")
#'
#' # Example 2: Two-sided test
#' u_sample2 <- rnorm(1000, mean = -0.5, sd = 0.3)
#' p_two_sided <- weighted.p.value(u_sample2, mu = 0, sigma = 1,
#'                                  alternative = "two.sided")
#' cat("Two-sided weighted p-value:", p_two_sided, "\n")
#'
#' # Example 3: Using with bootstrap samples
#' # Simulate some data and bootstrap the mean
#' set.seed(123)
#' original_data <- rnorm(30, mean = 1.2, sd = 2)
#' boot_means <- replicate(1000, mean(sample(original_data, replace = TRUE)))
#'
#' # Weighted p-value using bootstrap distribution
#' p_boot <- weighted.p.value(boot_means, mu = 0, sigma = 2/sqrt(30))
#' cat("Bootstrap-based weighted p-value:", p_boot, "\n")
#'
#' @export
weighted.p.value <- function(u, mu, sigma, alternative = c("two.sided", "less", "greater")) {
    # Input validation
    if (!is.numeric(u) || length(u) == 0) {
        stop("'u' must be a non-empty numeric vector")
    }
    if (any(!is.finite(u))) {
        stop("'u' contains non-finite values (NA, NaN, Inf)")
    }
    if (!is.numeric(mu) || length(mu) != 1) {
        stop("'mu' must be a single numeric value")
    }
    if (!is.finite(mu)) {
        stop("'mu' must be finite")
    }
    if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
        stop("'sigma' must be a single positive numeric value")
    }
    if (!is.finite(sigma)) {
        stop("'sigma' must be finite")
    }

    # Match alternative hypothesis
    alternative <- match.arg(alternative)

    # Compute p-values based on alternative hypothesis
    # Using direct vectorization for efficiency
    if (alternative == "less") {
        # H1: true mean < mu
        p.values <- pnorm(u, mean = mu, sd = sigma, lower.tail = TRUE)
    } else if (alternative == "greater") {
        # H1: true mean > mu
        p.values <- pnorm(u, mean = mu, sd = sigma, lower.tail = FALSE)
    } else {  # two.sided
        # H1: true mean != mu
        # For two-sided test, p-value = 2 * min(P(X <= x), P(X >= x))
        p.lower <- pnorm(u, mean = mu, sd = sigma, lower.tail = TRUE)
        p.upper <- pnorm(u, mean = mu, sd = sigma, lower.tail = FALSE)
        p.values <- 2 * pmin(p.lower, p.upper)
    }

    # Compute the weighted p-value as the mean of individual p-values
    weighted.p <- mean(p.values)

    # Ensure p-value is in [0, 1] (numerical precision safeguard)
    weighted.p <- pmax(0, pmin(1, weighted.p))

    return(weighted.p)
}


#' Weighted P-value Summary
#'
#' Provides a comprehensive summary of weighted p-value analysis including
#' comparison with classical p-value and visualization options.
#'
#' @param u     A numeric vector representing the uncertainty distribution.
#' @param mu    The mean of the normal null distribution.
#' @param sigma The standard deviation of the normal null distribution.
#' @param alternative Character string specifying the alternative hypothesis.
#' @param plot  Logical; if TRUE, produces a visualization of the analysis.
#'
#' @return A list containing:
#'   \item{weighted.p.value}{The weighted p-value}
#'   \item{classical.p.value}{Classical p-value based on mean of u}
#'   \item{summary.stats}{Summary statistics of the uncertainty distribution}
#'   \item{interpretation}{Text interpretation of results}
#'
#' @examples
#' u_sample <- rnorm(1000, mean = 1.5, sd = 0.5)
#' summary_results <- weighted.p.value.summary(u_sample, mu = 0, sigma = 1)
#' print(summary_results)
#'
#' @export
weighted.p.value.summary <- function(u, mu, sigma,
                                   alternative = c("two.sided", "less", "greater"),
                                   plot = FALSE) {
    # Match alternative
    alternative <- match.arg(alternative)

    # Calculate weighted p-value
    weighted.p <- weighted.p.value(u, mu, sigma, alternative)

    # Calculate classical p-value using mean of u
    u.mean <- mean(u)
    if (alternative == "less") {
        classical.p <- pnorm(u.mean, mean = mu, sd = sigma, lower.tail = TRUE)
    } else if (alternative == "greater") {
        classical.p <- pnorm(u.mean, mean = mu, sd = sigma, lower.tail = FALSE)
    } else {  # two.sided
        z.score <- abs(u.mean - mu) / sigma
        classical.p <- 2 * pnorm(-z.score)
    }

    # Summary statistics
    summary.stats <- c(
        mean = mean(u),
        sd = sd(u),
        median = median(u),
        q025 = quantile(u, 0.025),
        q975 = quantile(u, 0.975)
    )

    # Interpretation
    sig.level <- 0.05
    weighted.sig <- ifelse(weighted.p < sig.level, "significant", "not significant")
    classical.sig <- ifelse(classical.p < sig.level, "significant", "not significant")

    interpretation <- paste0(
        "Alternative hypothesis: ",
        switch(alternative,
               "two.sided" = "true mean is not equal to ",
               "less" = "true mean is less than ",
               "greater" = "true mean is greater than "),
        mu, "\n",
        "Weighted p-value: ", format(weighted.p, digits = 4),
        " (", weighted.sig, " at alpha = ", sig.level, ")\n",
        "Classical p-value: ", format(classical.p, digits = 4),
        " (", classical.sig, " at alpha = ", sig.level, ")\n"
    )

    # Optional plotting
    if (plot) {
        # Create visualization
        par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

        # Plot 1: Uncertainty distribution with null distribution overlay
        hist(u, breaks = 30, freq = FALSE, col = "lightblue",
             main = "Uncertainty Distribution vs Null Distribution",
             xlab = "Value", ylab = "Density")
        x.range <- range(c(u, mu - 3*sigma, mu + 3*sigma))
        x.seq <- seq(x.range[1], x.range[2], length.out = 200)
        lines(x.seq, dnorm(x.seq, mean = mu, sd = sigma),
              col = "red", lwd = 2)
        abline(v = mu, col = "red", lty = 2)
        abline(v = mean(u), col = "blue", lty = 2)
        legend("topright",
               legend = c("Uncertainty dist.", "Null dist.", "Null mean", "Sample mean"),
               col = c("lightblue", "red", "red", "blue"),
               lty = c(NA, 1, 2, 2), pch = c(15, NA, NA, NA),
               lwd = c(NA, 2, 1, 1))

        # Plot 2: Distribution of p-values
        p.values <- numeric(length(u))
        for (i in seq_along(u)) {
            if (alternative == "less") {
                p.values[i] <- pnorm(u[i], mean = mu, sd = sigma, lower.tail = TRUE)
            } else if (alternative == "greater") {
                p.values[i] <- pnorm(u[i], mean = mu, sd = sigma, lower.tail = FALSE)
            } else {
                p.lower <- pnorm(u[i], mean = mu, sd = sigma, lower.tail = TRUE)
                p.upper <- pnorm(u[i], mean = mu, sd = sigma, lower.tail = FALSE)
                p.values[i] <- 2 * min(p.lower, p.upper)
            }
        }

        hist(p.values, breaks = 30, freq = FALSE, col = "lightgreen",
             main = "Distribution of Individual P-values",
             xlab = "P-value", ylab = "Density", xlim = c(0, 1))
        abline(v = weighted.p, col = "darkgreen", lwd = 2)
        abline(v = classical.p, col = "orange", lwd = 2)
        abline(v = sig.level, col = "red", lty = 3)
        legend("topright",
               legend = c("Weighted p-value", "Classical p-value", "alpha = 0.05"),
               col = c("darkgreen", "orange", "red"),
               lty = c(1, 1, 3), lwd = c(2, 2, 1))

        par(mfrow = c(1, 1))  # Reset plotting parameters
    }

    # Return results
    result <- list(
        weighted.p.value = weighted.p,
        classical.p.value = classical.p,
        summary.stats = summary.stats,
        interpretation = interpretation,
        alternative = alternative,
        null.parameters = c(mu = mu, sigma = sigma)
    )

    class(result) <- "weighted.p.summary"
    return(result)
}


#' Print method for weighted.p.summary objects
#'
#' @param x A weighted.p.summary object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisible copy of x
#' @export
print.weighted.p.summary <- function(x, ...) {
    cat("Weighted P-value Analysis\n")
    cat("========================\n\n")
    cat("Null hypothesis: X ~ N(", x$null.parameters["mu"], ", ",
        x$null.parameters["sigma"], ")\n")
    cat(x$interpretation, "\n")
    cat("Summary of uncertainty distribution:\n")
    print(round(x$summary.stats, 4))
    cat("\n")
    invisible(x)
}

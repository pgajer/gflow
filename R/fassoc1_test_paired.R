## ============================================================================
## First-Order Functional Association Test with Paired Bayesian Bootstrap
## ============================================================================
##
## This implementation uses PAIRED Bayesian bootstrap weights for signal and
## null samples, ensuring that any difference is purely due to the permutation
## destroying the x-y association, not due to different random weights.
##
## Key features:
## - Paired BB weights: same lambda used for signal and null
## - Three test options: paired.t, weighted.pvalue, wilcoxon
## - Automatic Box-Cox transformation with normality testing
##
## @author Pawel Gajer
## @date 2025
## ============================================================================


#' Tests for First-Order Functional Association Between Two Variables
#'
#' Computes and tests the first-order functional association measure between two
#' variables x and y. The function estimates the integral of |dE_x(y)/dx|, where
#' E_x(y) is the conditional mean of y given x. This measure captures the strength
#' of the derivative relationship between variables.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x). Can be binary
#'   or continuous.
#' @param test.type Character string specifying the test for inference. Options are:
#'   \describe{
#'     \item{"paired.t"}{Paired t-test on differences (signal - null). Most powerful
#'       when differences are approximately normal.}
#'     \item{"weighted.pvalue"}{Weighted p-value integrating over signal uncertainty.
#'       Assumes null distribution is approximately normal after Box-Cox.}
#'     \item{"wilcoxon"}{Wilcoxon signed-rank test on differences. Distribution-free,
#'       most robust to non-normality.}
#'   }
#'   Default is "paired.t".
#' @param boxcox Character string controlling Box-Cox transformation:
#'   \describe{
#'     \item{"auto"}{Apply Box-Cox if Shapiro-Wilk test rejects normality (default)}
#'     \item{"always"}{Always apply Box-Cox transformation}
#'     \item{"never"}{Never apply Box-Cox transformation}
#'   }
#' @param boxcox.alpha Significance level for Shapiro-Wilk normality test when
#'   boxcox = "auto". Default is 0.05.
#' @param bw Numeric bandwidth parameter for smoothing. If NULL (default), bandwidth
#'   is automatically selected.
#' @param n.BB Integer specifying the number of paired Bayesian bootstrap samples.
#'   Default is 1000. Each BB sample generates one signal and one null value using
#'   the same weights.
#' @param grid.size Integer specifying the size of evaluation grid. Default is 400.
#' @param degree Integer specifying the polynomial degree (1 or 2). Default is 1.
#' @param min.K Integer specifying minimum number of observations for local
#'   estimation. Default is 5.
#' @param n.cores Integer specifying the number of CPU cores for parallel
#'   computation. Default is 1 (not yet implemented).
#' @param seed Integer random seed for reproducibility. Default is NULL.
#' @param plot.it Logical indicating whether to produce diagnostic plots.
#'   Default is TRUE.
#' @param xlab Character string for x-axis label. Default is "x".
#' @param ylab Character string for y-axis label. Default is "y".
#' @param verbose Logical indicating whether to print progress messages.
#'   Default is TRUE.
#'
#' @return An object of class "assoc1" containing:
#'   \describe{
#'     \item{Delta1}{The total change in conditional mean (E_y(x_max) - E_y(x_min))}
#'     \item{delta1}{The first-order functional association measure (point estimate)}
#'     \item{delta1.z}{Normalized effect size: (delta1 - null.mean) / null.sd}
#'     \item{delta1.robust.z}{Robust normalized effect size: (delta1 - null.median) / null.mad}
#'     \item{p.value}{The p-value for the test of association}
#'     \item{log.p.value}{Natural logarithm of the p-value}
#'     \item{test.type}{The test type used}
#'     \item{signal.delta1}{Vector of signal BB samples}
#'     \item{null.delta1}{Vector of null BB samples (paired with signal)}
#'     \item{diff.delta1}{Vector of differences (signal - null)}
#'     \item{null.mean}{Mean of null distribution}
#'     \item{null.sd}{Standard deviation of null distribution}
#'     \item{null.median}{Median of null distribution}
#'     \item{null.mad}{Median absolute deviation of null distribution (scaled)}
#'     \item{boxcox.applied}{Logical indicating if Box-Cox was applied}
#'     \item{boxcox.lambda}{Box-Cox lambda parameter (if applied)}
#'     \item{shapiro.pvalue.raw}{Shapiro-Wilk p-value on raw differences}
#'     \item{shapiro.pvalue.bc}{Shapiro-Wilk p-value after Box-Cox (if applied)}
#'     \item{signal.Eyg}{Matrix of signal gpredictions (ng x n.BB)}
#'     \item{null.Eyg}{Matrix of null gpredictions (ng x n.BB)}
#'     \item{x}{Input predictor values (sorted)}
#'     \item{y}{Input response values (sorted by x)}
#'     \item{xgrid}{Grid points}
#'     \item{n.BB}{Number of BB samples}
#'   }
#'
#' @details
#' The first-order functional association measure captures the variability in the
#' derivative of the conditional mean function. Unlike the zero-order measure which
#' looks at deviations from the mean, this measure is sensitive to the rate of
#' change in the relationship between x and y.
#'
#' This implementation uses PAIRED Bayesian bootstrap weights: for each bootstrap
#' index b, the same weight vector lambda_b is used for both the signal (original y)
#' and null (permuted y) computations. The only difference between paired samples
#' is the permutation destroying the x-y association. This pairing reduces variance
#' and ensures the test is comparing like with like.
#'
#' The test statistic is delta1 = sum(|diff(Eyg)|), the sum of absolute differences
#' in the conditional mean across the grid. Under H0 (no association), permuting y
#' should not systematically change this measure.
#'
#' @examples
#' \dontrun{
#' # Generate example data with nonlinear relationship
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(4 * pi * x) + rnorm(n, sd = 0.2)
#'
#' # Test for first-order functional association
#' result <- fassoc1.test(x, y, n.BB = 500, seed = 42)
#'
#' # View results
#' print(result)
#'
#' # Plot diagnostics
#' plot(result, type = "diff")
#' }
#'
#' @importFrom stats shapiro.test t.test wilcox.test pnorm sd quantile
#' @importFrom graphics plot points lines abline polygon hist par legend
#' @export
fassoc1.test <- function(x,
                         y,
                         test.type = c("paired.t", "weighted.pvalue", "wilcoxon"),
                         boxcox = c("auto", "always", "never"),
                         boxcox.alpha = 0.05,
                         bw = NULL,
                         n.BB = 1000,
                         grid.size = 400,
                         degree = 1,
                         min.K = 5,
                         n.cores = 1,
                         seed = NULL,
                         plot.it = TRUE,
                         xlab = "x",
                         ylab = "y",
                         verbose = TRUE) {

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    test.type <- match.arg(test.type)
    boxcox <- match.arg(boxcox)

    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Both 'x' and 'y' must be numeric vectors")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (!is.null(bw) && (!is.numeric(bw) || bw <= 0)) {
        stop("'bw' must be NULL or a positive numeric value")
    }

    if (!is.numeric(n.BB) || n.BB < 100) {
        stop("'n.BB' must be at least 100")
    }

    if (!is.numeric(grid.size) || grid.size < 10) {
        stop("'grid.size' must be at least 10")
    }

    if (!is.numeric(degree) || !(degree %in% c(1, 2))) {
        stop("'degree' must be 1 or 2")
    }

    if (!is.numeric(min.K) || min.K < 2) {
        stop("'min.K' must be at least 2")
    }

    if (!is.numeric(boxcox.alpha) || boxcox.alpha <= 0 || boxcox.alpha >= 1) {
        stop("'boxcox.alpha' must be between 0 and 1")
    }

    ## Convert to integers
    n.BB <- as.integer(n.BB)
    grid.size <- as.integer(grid.size)
    degree <- as.integer(degree)
    min.K <- as.integer(min.K)

    ## Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (verbose) {
        routine.ptm <- proc.time()
        message("First-Order Functional Association Test (Paired BB)")
        message(sprintf("  n = %d, n.BB = %d, test.type = %s",
                        length(x), n.BB, test.type))
    }

    nx <- length(x)

    ## Check for binary y
    y.binary <- all(y %in% c(0, 1))

    ## Remove non-finite values
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < length(x)) {
        warning(sprintf("Removed %d non-finite values", length(x) - sum(idx)))
        x <- x[idx]
        y <- y[idx]
        nx <- length(x)
    }

    if (nx <= 2 * min.K) {
        stop(sprintf("Sample size (%d) must be greater than 2 * min.K (%d)",
                     nx, 2 * min.K))
    }

    ## ========================================================================
    ## Generate shared BB weights
    ## ========================================================================

    if (verbose) {
        message("Generating Dirichlet weights...")
        ptm <- proc.time()
    }

    ## Generate weights matrix: nx rows, n.BB columns
    ## Each column is a Dirichlet(1,...,1) sample scaled to sum to nx
    if (exists("generate.dirichlet.weights")) {
        lambda <- generate.dirichlet.weights(nx, n.BB)
    } else {
        ## Fallback R implementation
        lambda <- matrix(0, nrow = nx, ncol = n.BB)
        for (b in seq_len(n.BB)) {
            e <- rexp(nx, rate = 1)
            lambda[, b] <- e / sum(e) * nx
        }
    }

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Compute signal BB samples
    ## ========================================================================

    if (verbose) {
        message("Computing signal BB samples...")
        ptm <- proc.time()
    }

    if (exists("magelo.with.external.BB")) {
        signal.fit <- magelo.with.external.BB(x, y, lambda,
                                              bw = bw,
                                              grid.size = grid.size,
                                              degree = degree,
                                              min.K = min.K)
    } else {
        ## Fallback: use magelo if available, but weights won't be paired
        warning("magelo.with.external.BB not found. Using unpaired BB weights.")

        if (!exists("magelo")) {
            stop("Neither magelo.with.external.BB nor magelo found")
        }

        signal.fit <- magelo(x, y, bw = bw, grid.size = grid.size,
                             degree = degree, min.K = min.K,
                             n.BB = n.BB, get.BB.gpredictions = TRUE)
        signal.fit$BB.gpredictions <- signal.fit$BB.gpredictions
    }

    xgrid <- signal.fit$xgrid
    signal.Eyg <- signal.fit$BB.gpredictions
    bw.used <- signal.fit$opt.bw

    ## Compute delta1 for each BB sample: sum of absolute derivatives
    signal.delta1 <- apply(signal.Eyg, 2, function(z) sum(abs(diff(z))))

    ## Point estimate
    if (!is.null(signal.fit$gpredictions)) {
        Eyg <- signal.fit$gpredictions
        delta1 <- sum(abs(diff(Eyg)))
        Delta1 <- Eyg[length(Eyg)] - Eyg[1]
    } else {
        ## Use median of BB samples
        delta1 <- median(signal.delta1)
        Delta1 <- NA
    }

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Compute null BB samples (permuted y, SAME weights)
    ## ========================================================================

    if (verbose) {
        message("Computing null BB samples with paired weights...")
        ptm <- proc.time()
    }

    ## Generate one permutation for all BB samples
    ## (We could use different permutations for each BB, but using one
    ## permutation keeps the pairing cleaner)
    y.perm <- sample(y)

    if (exists("magelo.with.external.BB")) {
        null.fit <- magelo.with.external.BB(x, y.perm, lambda,
                                            bw = bw.used,  # Use same bw
                                            grid.size = grid.size,
                                            degree = degree,
                                            min.K = min.K)
    } else {
        ## Fallback with unpaired weights
        null.fit <- magelo(x, y.perm, bw = bw.used, grid.size = grid.size,
                           degree = degree, min.K = min.K,
                           n.BB = n.BB, get.BB.gpredictions = TRUE)
    }

    null.Eyg <- null.fit$BB.gpredictions
    null.delta1 <- apply(null.Eyg, 2, function(z) sum(abs(diff(z))))

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Compute paired differences
    ## ========================================================================

    diff.delta1 <- signal.delta1 - null.delta1

    ## ========================================================================
    ## Normality assessment and Box-Cox transformation
    ## ========================================================================

    if (verbose) {
        message("Assessing normality...")
    }

    ## Test normality of differences
    shapiro.raw <- tryCatch(
        shapiro.test(diff.delta1),
        error = function(e) list(p.value = NA)
    )
    shapiro.pvalue.raw <- shapiro.raw$p.value

    ## Decide on Box-Cox
    boxcox.applied <- FALSE
    boxcox.lambda <- NA
    shapiro.pvalue.bc <- NA

    apply.boxcox <- FALSE
    if (boxcox == "always") {
        apply.boxcox <- TRUE
    } else if (boxcox == "auto") {
        if (!is.na(shapiro.pvalue.raw) && shapiro.pvalue.raw < boxcox.alpha) {
            apply.boxcox <- TRUE
        }
    }

    ## For Box-Cox on differences, handle negative values by shifting
    diff.for.test <- diff.delta1

    if (apply.boxcox) {
        diff.for.bc <- diff.delta1
        shift <- 0

        if (min(diff.for.bc) <= 0) {
            shift <- abs(min(diff.for.bc)) + 0.01 * sd(diff.for.bc)
            diff.for.bc <- diff.for.bc + shift
        }

        ## Find optimal lambda via Box-Cox MLE
        if (exists("boxcox.mle")) {
            bc.fit <- tryCatch(
                boxcox.mle(diff.for.bc ~ 1),
                error = function(e) NULL
            )

            if (!is.null(bc.fit)) {
                boxcox.lambda <- bc.fit$lambda
                boxcox.applied <- TRUE

                ## Apply transformation
                if (abs(boxcox.lambda) > .Machine$double.eps) {
                    diff.for.test <- (diff.for.bc^boxcox.lambda - 1) / boxcox.lambda
                } else {
                    diff.for.test <- log(diff.for.bc)
                }

                ## Test normality after transformation
                shapiro.bc <- tryCatch(
                    shapiro.test(diff.for.test),
                    error = function(e) list(p.value = NA)
                )
                shapiro.pvalue.bc <- shapiro.bc$p.value
            }
        } else {
            ## boxcox.mle not available, try log transform
            if (all(diff.for.bc > 0)) {
                diff.for.test <- log(diff.for.bc)
                boxcox.lambda <- 0
                boxcox.applied <- TRUE

                shapiro.bc <- tryCatch(
                    shapiro.test(diff.for.test),
                    error = function(e) list(p.value = NA)
                )
                shapiro.pvalue.bc <- shapiro.bc$p.value
            }
        }
    }

    ## ========================================================================
    ## Perform hypothesis test
    ## ========================================================================

    if (verbose) {
        message(sprintf("Performing %s test...", test.type))
    }

    if (test.type == "paired.t") {
        ## Paired t-test on differences
        ## H0: mean(diff) = 0
        ## H1: mean(diff) > 0 (signal has larger association than null)

        t.result <- t.test(diff.for.test, mu = 0, alternative = "greater")
        p.value <- t.result$p.value

    } else if (test.type == "weighted.pvalue") {
        ## Weighted p-value approach
        ## Fit normal to null distribution, compute weighted p-value over signal

        ## Transform values if Box-Cox was applied
        if (boxcox.applied) {
            null.for.test <- null.delta1
            shift <- 0
            if (min(null.for.test) <= 0) {
                shift <- abs(min(null.for.test)) + 0.01 * sd(null.for.test)
                null.for.test <- null.for.test + shift
            }

            if (abs(boxcox.lambda) > .Machine$double.eps) {
                null.for.test <- (null.for.test^boxcox.lambda - 1) / boxcox.lambda
            } else {
                null.for.test <- log(null.for.test)
            }

            signal.for.test <- signal.delta1
            if (min(signal.for.test) <= 0) {
                signal.for.test <- signal.for.test + shift
            }
            if (abs(boxcox.lambda) > .Machine$double.eps) {
                signal.for.test <- (signal.for.test^boxcox.lambda - 1) / boxcox.lambda
            } else {
                signal.for.test <- log(signal.for.test)
            }
        } else {
            null.for.test <- null.delta1
            signal.for.test <- signal.delta1
        }

        mu <- mean(null.for.test)
        sigma <- sd(null.for.test)

        if (exists("weighted.p.value")) {
            p.value <- weighted.p.value(signal.for.test, mu, sigma,
                                        alternative = "greater")
        } else {
            ## Fallback: compute average p-value manually
            p.values <- pnorm(signal.for.test, mean = mu, sd = sigma,
                              lower.tail = FALSE)
            p.value <- mean(p.values)
        }

    } else if (test.type == "wilcoxon") {
        ## Wilcoxon signed-rank test on differences
        ## H0: median(diff) = 0
        ## H1: median(diff) > 0

        wilcox.result <- wilcox.test(diff.delta1, mu = 0, alternative = "greater")
        p.value <- wilcox.result$p.value
    }

    log.p.value <- log(max(p.value, .Machine$double.eps))

    ## ========================================================================
    ## Compute normalized effect sizes
    ## ========================================================================

    ## Null distribution statistics
    null.mean <- mean(null.delta1)
    null.sd <- sd(null.delta1)
    null.median <- median(null.delta1)
    null.mad <- mad(null.delta1, constant = 1.4826)  # scaled to be consistent with SD

    ## Normalized delta1: z-score relative to null
    ## How many SDs is the point estimate from the null mean?
    if (null.sd > .Machine$double.eps) {
        delta1.z <- (delta1 - null.mean) / null.sd
    } else {
        delta1.z <- NA
    }

    ## Robust normalized delta1: using median and MAD
    ## More robust to outliers in null distribution
    if (null.mad > .Machine$double.eps) {
        delta1.robust.z <- (delta1 - null.median) / null.mad
    } else {
        delta1.robust.z <- NA
    }

    ## ========================================================================
    ## Create diagnostic plots
    ## ========================================================================

    if (plot.it) {
        op <- par(mfrow = c(1, 2), mar = c(3.5, 3.5, 2, 0.5),
                  mgp = c(2, 0.5, 0), tcl = -0.3)
        on.exit(par(op), add = TRUE)

        ## Plot 1: Conditional mean estimates
        signal.Eyg.mean <- rowMeans(signal.Eyg)
        signal.Eyg.q <- apply(signal.Eyg, 1, quantile, probs = c(0.025, 0.975))

        null.Eyg.mean <- rowMeans(null.Eyg)
        null.Eyg.q <- apply(null.Eyg, 1, quantile, probs = c(0.025, 0.975))

        ylim <- range(c(y, signal.Eyg.q, null.Eyg.q), na.rm = TRUE)

        plot(x, y, las = 1, ylim = ylim, xlab = xlab, ylab = ylab,
             main = "Conditional Mean", col = "gray70", pch = 16, cex = 0.5)

        ## Null band
        polygon(c(xgrid, rev(xgrid)),
                c(null.Eyg.q[1, ], rev(null.Eyg.q[2, ])),
                col = adjustcolor("blue", alpha.f = 0.2), border = NA)

        ## Signal band
        polygon(c(xgrid, rev(xgrid)),
                c(signal.Eyg.q[1, ], rev(signal.Eyg.q[2, ])),
                col = adjustcolor("red", alpha.f = 0.2), border = NA)

        lines(xgrid, signal.Eyg.mean, col = "red", lwd = 2)
        lines(xgrid, null.Eyg.mean, col = "blue", lwd = 2, lty = 2)

        legend("topright",
               legend = c("Signal E(y|x)", "Null E(y|x)"),
               col = c("red", "blue"),
               lty = c(1, 2), lwd = 2, cex = 0.8, bty = "n")

        ## Plot 2: Histogram of paired differences
        hist(diff.delta1, breaks = 30, main = "Paired Differences",
             xlab = expression(delta[1]^signal - delta[1]^null),
             col = "lightblue", border = "white", las = 1)
        abline(v = 0, col = "red", lwd = 2, lty = 2)
        abline(v = mean(diff.delta1), col = "darkblue", lwd = 2)

        ## Add p-value text
        text.x <- quantile(diff.delta1, 0.95)
        text.y <- par("usr")[4] * 0.9
        text(text.x, text.y, labels = sprintf("p = %.4f", p.value),
             adj = c(1, 1), cex = 0.9)
    }

    if (verbose) {
        message(sprintf("Done. p-value = %.4g", p.value))
        message(sprintf("Total time: %.2f sec",
                        (proc.time() - routine.ptm)[3]))
    }

    ## ========================================================================
    ## Prepare output
    ## ========================================================================

    output <- list(
        Delta1 = Delta1,
        delta1 = delta1,
        delta1.z = delta1.z,
        delta1.robust.z = delta1.robust.z,
        p.value = p.value,
        log.p.value = log.p.value,
        test.type = test.type,
        signal.delta1 = signal.delta1,
        null.delta1 = null.delta1,
        diff.delta1 = diff.delta1,
        null.mean = null.mean,
        null.sd = null.sd,
        null.median = null.median,
        null.mad = null.mad,
        boxcox.applied = boxcox.applied,
        boxcox.lambda = boxcox.lambda,
        shapiro.pvalue.raw = shapiro.pvalue.raw,
        shapiro.pvalue.bc = shapiro.pvalue.bc,
        signal.Eyg = signal.Eyg,
        null.Eyg = null.Eyg,
        x = x,
        y = y,
        xgrid = xgrid,
        opt.bw = bw.used,
        grid.size = grid.size,
        degree = degree,
        n.BB = n.BB,
        lambda = lambda,
        call = match.call()
    )

    class(output) <- "assoc1"

    return(output)
}


#' @rdname fassoc1.test
#' @export
fofam.test <- fassoc1.test


## ============================================================================
## S3 Methods
## ============================================================================

#' Coef Method for assoc1 Objects
#'
#' Extracts key coefficients from an assoc1 object.
#'
#' @param object An object of class "assoc1".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector containing delta1, normalized effect sizes,
#'   and p-value.
#' @method coef assoc1
#' @export
coef.assoc1 <- function(object, ...) {
    c(delta1 = object$delta1,
      Delta1 = object$Delta1,
      delta1.z = object$delta1.z,
      delta1.robust.z = object$delta1.robust.z,
      p.value = object$p.value)
}


#' Print Method for assoc1 Objects
#'
#' @param x An object of class "assoc1".
#' @param digits Number of significant digits to display.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible x.
#' @method print assoc1
#' @export
print.assoc1 <- function(x, digits = 4, ...) {
    cat("\nFirst-Order Functional Association Test (Paired BB)\n")
    cat("===================================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Test Configuration:\n")
    cat("  Sample size: ", length(x$x), "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n", sep = "")
    cat("  BB samples: ", x$n.BB, "\n", sep = "")
    cat("  Polynomial degree: ", x$degree, "\n", sep = "")
    cat("  Box-Cox applied: ", x$boxcox.applied, "\n", sep = "")
    if (x$boxcox.applied) {
        cat("  Box-Cox lambda: ", signif(x$boxcox.lambda, digits), "\n", sep = "")
    }
    cat("\n")

    cat("Normality Assessment (Shapiro-Wilk):\n")
    cat("  Raw differences: p = ", signif(x$shapiro.pvalue.raw, digits), "\n", sep = "")
    if (x$boxcox.applied) {
        cat("  After Box-Cox: p = ", signif(x$shapiro.pvalue.bc, digits), "\n", sep = "")
    }
    cat("\n")

    cat("Results:\n")
    if (!is.na(x$Delta1)) {
        cat("  Total change (Delta1): ", signif(x$Delta1, digits), "\n", sep = "")
    }
    cat("  Association measure (delta1): ", signif(x$delta1, digits), "\n", sep = "")
    cat("  Normalized effect size (z): ", signif(x$delta1.z, digits), "\n", sep = "")
    cat("  Robust effect size (z.robust): ", signif(x$delta1.robust.z, digits), "\n", sep = "")
    cat("  Mean difference (signal - null): ", signif(mean(x$diff.delta1), digits), "\n", sep = "")
    cat("  p-value: ", signif(x$p.value, digits), "\n", sep = "")
    cat("  log(p-value): ", signif(x$log.p.value, digits), "\n", sep = "")

    invisible(x)
}


#' Summary Method for assoc1 Objects
#'
#' @param object An object of class "assoc1".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class "summary.assoc1" with summary statistics.
#' @method summary assoc1
#' @export
summary.assoc1 <- function(object, ...) {
    out <- list(
        call = object$call,
        n = length(object$x),
        test.type = object$test.type,
        n.BB = object$n.BB,
        degree = object$degree,
        Delta1 = object$Delta1,
        delta1 = object$delta1,
        delta1.z = object$delta1.z,
        delta1.robust.z = object$delta1.robust.z,
        p.value = object$p.value,
        null.mean = object$null.mean,
        null.sd = object$null.sd,
        null.median = object$null.median,
        null.mad = object$null.mad,
        boxcox.applied = object$boxcox.applied,
        boxcox.lambda = object$boxcox.lambda,
        shapiro.raw = object$shapiro.pvalue.raw,
        shapiro.bc = object$shapiro.pvalue.bc,
        signal.summary = summary(object$signal.delta1),
        null.summary = summary(object$null.delta1),
        diff.summary = summary(object$diff.delta1),
        diff.quantiles = quantile(object$diff.delta1,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    )

    class(out) <- "summary.assoc1"
    return(out)
}


#' @method print summary.assoc1
#' @export
print.summary.assoc1 <- function(x, digits = 4, ...) {
    cat("\nFirst-Order Functional Association Test Summary\n")
    cat("===============================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Configuration:\n")
    cat("  Sample size: ", x$n, "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n", sep = "")
    cat("  BB samples: ", x$n.BB, "\n", sep = "")
    cat("  Box-Cox applied: ", x$boxcox.applied, "\n\n", sep = "")

    cat("Null Distribution:\n")
    cat("  Mean: ", signif(x$null.mean, digits), "\n", sep = "")
    cat("  SD: ", signif(x$null.sd, digits), "\n", sep = "")
    cat("  Median: ", signif(x$null.median, digits), "\n", sep = "")
    cat("  MAD: ", signif(x$null.mad, digits), "\n\n", sep = "")

    cat("Results:\n")
    cat("  Association measure (delta1): ", signif(x$delta1, digits), "\n", sep = "")
    cat("  Normalized effect size (z): ", signif(x$delta1.z, digits), "\n", sep = "")
    cat("  Robust effect size (z.robust): ", signif(x$delta1.robust.z, digits), "\n", sep = "")
    cat("  p-value: ", signif(x$p.value, digits), "\n\n", sep = "")

    cat("Signal Distribution:\n")
    print(signif(x$signal.summary, digits))
    cat("\n")

    cat("Null Distribution:\n")
    print(signif(x$null.summary, digits))
    cat("\n")

    cat("Paired Difference (Signal - Null) Quantiles:\n")
    print(signif(x$diff.quantiles, digits))

    invisible(x)
}


#' Plot Method for assoc1 Objects
#'
#' Creates diagnostic plots for objects of class "assoc1".
#'
#' @param x An object of class "assoc1".
#' @param type Character string specifying plot type. Options are:
#'   \describe{
#'     \item{"Exy"}{Conditional mean function with credible intervals (default)}
#'     \item{"dExy"}{Derivative of conditional mean}
#'     \item{"diff"}{Histogram of signal - null differences}
#'     \item{"qq"}{QQ plot of differences}
#'     \item{"comparison"}{Side-by-side signal and null distributions}
#'   }
#' @param ... Additional arguments passed to plotting functions.
#'
#' @method plot assoc1
#' @export
plot.assoc1 <- function(x, type = c("Exy", "dExy", "diff", "qq", "comparison"),
                         ...) {
    type <- match.arg(type)

    if (type == "Exy") {
        ## Conditional mean with CrI
        xgrid <- x$xgrid

        signal.mean <- rowMeans(x$signal.Eyg)
        signal.q <- apply(x$signal.Eyg, 1, quantile, probs = c(0.025, 0.975))

        null.mean <- rowMeans(x$null.Eyg)
        null.q <- apply(x$null.Eyg, 1, quantile, probs = c(0.025, 0.975))

        ylim <- range(c(x$y, signal.q, null.q), na.rm = TRUE)

        plot(x$x, x$y, xlab = "x", ylab = "y", main = "Conditional Mean",
             ylim = ylim, las = 1, col = "gray70", pch = 16, cex = 0.5, ...)

        polygon(c(xgrid, rev(xgrid)), c(null.q[1, ], rev(null.q[2, ])),
                col = adjustcolor("blue", alpha.f = 0.2), border = NA)
        polygon(c(xgrid, rev(xgrid)), c(signal.q[1, ], rev(signal.q[2, ])),
                col = adjustcolor("red", alpha.f = 0.2), border = NA)

        lines(xgrid, signal.mean, col = "red", lwd = 2)
        lines(xgrid, null.mean, col = "blue", lwd = 2, lty = 2)

    } else if (type == "dExy") {
        ## Derivative plot
        xgrid <- x$xgrid
        dxgrid <- xgrid[-1]

        signal.dEyg <- apply(x$signal.Eyg, 2, diff)
        null.dEyg <- apply(x$null.Eyg, 2, diff)

        signal.dEyg.mean <- rowMeans(signal.dEyg)
        signal.dEyg.q <- apply(signal.dEyg, 1, quantile, probs = c(0.025, 0.975))

        null.dEyg.q <- apply(null.dEyg, 1, quantile, probs = c(0.025, 0.975))

        ylim <- range(c(signal.dEyg.q, null.dEyg.q), na.rm = TRUE)

        plot(dxgrid, signal.dEyg.mean, type = "l", ylim = ylim, xlab = "x",
             ylab = "dE(y|x)/dx", main = "Derivative", lwd = 2, col = "red",
             las = 1, ...)

        polygon(c(dxgrid, rev(dxgrid)),
                c(null.dEyg.q[1, ], rev(null.dEyg.q[2, ])),
                col = adjustcolor("blue", alpha.f = 0.2), border = NA)

        abline(h = 0, lty = 2, col = "gray")

    } else if (type == "diff") {
        ## Histogram of differences
        hist(x$diff.delta1, breaks = 30, main = "Paired Differences",
             xlab = expression(delta[1]^signal - delta[1]^null),
             col = "lightblue", border = "white", las = 1, ...)
        abline(v = 0, col = "red", lwd = 2, lty = 2)
        abline(v = mean(x$diff.delta1), col = "darkblue", lwd = 2)

    } else if (type == "qq") {
        ## QQ plot of differences
        qqnorm(x$diff.delta1, main = "QQ Plot of Differences", las = 1, ...)
        qqline(x$diff.delta1, col = "red")

    } else if (type == "comparison") {
        ## Side-by-side comparison
        op <- par(mfrow = c(1, 2))
        on.exit(par(op))

        ## Signal
        hist(x$signal.delta1, breaks = 30, main = "Signal Distribution",
             xlab = expression(delta[1]^signal), col = "lightcoral",
             border = "white", las = 1)
        abline(v = x$delta1, col = "red", lwd = 2)

        ## Null
        hist(x$null.delta1, breaks = 30, main = "Null Distribution",
             xlab = expression(delta[1]^null), col = "lightblue",
             border = "white", las = 1)
        abline(v = mean(x$null.delta1), col = "blue", lwd = 2)
    }

    invisible(NULL)
}


#' Extract Derivative Information from Functional Association Objects
#'
#' Generic function to extract derivative information from functional
#' association test results.
#'
#' @param object An object from a functional association test.
#' @param ... Additional arguments passed to methods.
#'
#' @return Derivative information (format depends on method).
#'
#' @seealso \code{\link{extract.derivatives.assoc1}}
#' @export
extract.derivatives <- function(object, ...) {
    UseMethod("extract.derivatives")
}


#' Extract Derivative Information from assoc1 Object
#'
#' Extracts and summarizes derivative information from the first-order
#' functional association test results.
#'
#' @param object An object of class "assoc1".
#' @param type Character string specifying what to extract:
#'   \describe{
#'     \item{"estimate"}{Point estimate of derivatives (default)}
#'     \item{"credible.interval"}{Credible intervals from BB samples}
#'     \item{"both"}{Both estimate and credible intervals}
#'   }
#' @param probs Numeric vector of probabilities for credible intervals.
#'   Default is c(0.025, 0.975) for 95% intervals.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Depending on type:
#'   \describe{
#'     \item{estimate}{A data.frame with columns x and dEy}
#'     \item{credible.interval}{A data.frame with columns x, lower, upper}
#'     \item{both}{A list containing both estimate and credible.interval data.frames}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- fassoc1.test(x, y)
#'
#' # Get point estimate
#' deriv.est <- extract.derivatives(result)
#'
#' # Get credible intervals
#' deriv.ci <- extract.derivatives(result, type = "credible.interval")
#'
#' # Get both
#' deriv.both <- extract.derivatives(result, type = "both")
#'
#' # Plot derivative with uncertainty
#' plot(deriv.both$estimate$x, deriv.both$estimate$dEy, type = "l")
#' polygon(c(deriv.both$credible.interval$x,
#'           rev(deriv.both$credible.interval$x)),
#'         c(deriv.both$credible.interval$lower,
#'           rev(deriv.both$credible.interval$upper)),
#'         col = "gray90", border = NA)
#' }
#'
#' @method extract.derivatives assoc1
#' @export
extract.derivatives.assoc1 <- function(object,
                                       type = c("estimate", "credible.interval", "both"),
                                       probs = c(0.025, 0.975),
                                       ...) {

    type <- match.arg(type)

    ## Grid points for derivatives (midpoints between grid points)
    xg <- object$xgrid
    ng <- length(xg)
    dxg <- (xg[-1] + xg[-ng]) / 2  # Midpoints

    ## Point estimate of derivatives
    ## Use mean of BB samples for robustness
    signal.dEyg <- apply(object$signal.Eyg, 2, diff)
    dEy <- rowMeans(signal.dEyg)

    if (type == "estimate") {
        return(data.frame(x = dxg, dEy = dEy))
    }

    ## Calculate credible intervals from BB samples
    if (type == "credible.interval" || type == "both") {
        dEy.CrI <- apply(signal.dEyg, 1, stats::quantile, probs = probs)

        ci.df <- data.frame(
            x = dxg,
            lower = dEy.CrI[1, ],
            upper = dEy.CrI[2, ]
        )

        if (type == "credible.interval") {
            return(ci.df)
        }
    }

    ## Return both
    list(
        estimate = data.frame(x = dxg, dEy = dEy),
        credible.interval = ci.df
    )
}

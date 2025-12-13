## ============================================================================
## Zero-Order Functional Association Test with Paired Bayesian Bootstrap
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


#' Tests for Zero-Order Functional Association Between Two Variables
#'
#' Computes and tests the zero-order functional association measure between two
#' variables x and y. The function estimates the mean of |E_x(y) - E(y)|, where
#' E_x(y) is the conditional mean of y given x and E(y) is the marginal mean.
#' This measure captures the deviation of the conditional mean from the overall
#' mean.
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
#' @param degree Integer specifying the polynomial degree (0, 1, or 2). Default is 1.
#'   Note: degree 0 uses weighted means, degree 1 or 2 uses local polynomial regression.
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
#' @return An object of class "assoc0" containing:
#'   \describe{
#'     \item{delta}{The zero-order functional association measure (point estimate)}
#'     \item{delta.z}{Normalized effect size: (delta - null.mean) / null.sd}
#'     \item{delta.robust.z}{Robust normalized effect size: (delta - null.median) / null.mad}
#'     \item{p.value}{The p-value for the test of association}
#'     \item{log.p.value}{Natural logarithm of the p-value}
#'     \item{test.type}{The test type used}
#'     \item{signal.delta}{Vector of signal BB samples}
#'     \item{null.delta}{Vector of null BB samples (paired with signal)}
#'     \item{diff.delta}{Vector of differences (signal - null)}
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
#'     \item{Ey}{Marginal mean of y}
#'     \item{x}{Input predictor values (sorted)}
#'     \item{y}{Input response values (sorted by x)}
#'     \item{xgrid}{Grid points}
#'     \item{n.BB}{Number of BB samples}
#'   }
#'
#' @details
#' The zero-order functional association measure captures the deviation of the
#' conditional mean from the marginal mean. It is defined as:
#'
#' \deqn{\delta_0 = \frac{1}{|X|} \int |E(y|x) - E(y)| dx}
#'
#' where E(y|x) is the conditional mean and E(y) is the marginal mean of y.
#' Under no association, E(y|x) = E(y) for all x, so delta_0 = 0.
#'
#' This implementation uses PAIRED Bayesian bootstrap weights: for each bootstrap
#' index b, the same weight vector lambda_b is used for both the signal (original y)
#' and null (permuted y) computations. The only difference between paired samples
#' is the permutation destroying the x-y association. This pairing reduces variance
#' and ensures the test is comparing like with like.
#'
#' @examples
#' \dontrun{
#' # Generate example data with mean shift relationship
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- ifelse(x > 0.5, 1, 0) + rnorm(n, sd = 0.3)
#'
#' # Test for zero-order functional association
#' result <- fassoc0.test(x, y, n.BB = 500, seed = 42)
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
fassoc0.test <- function(x,
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

    if (!is.numeric(degree) || !(degree %in% c(0, 1, 2))) {
        stop("'degree' must be 0, 1, or 2")
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
        message("Zero-Order Functional Association Test (Paired BB)")
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

    ## Marginal mean of y
    Ey <- mean(y)

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
    }

    xgrid <- signal.fit$xgrid
    signal.Eyg <- signal.fit$BB.gpredictions
    bw.used <- signal.fit$opt.bw

    ## Compute delta0 for each BB sample: mean absolute deviation from marginal mean
    ## Note: We compute the BB-weighted marginal mean for each sample
    signal.delta <- apply(signal.Eyg, 2, function(Eyg.b) {
        mean(abs(Eyg.b - mean(Eyg.b)))
    })

    ## Point estimate
    if (!is.null(signal.fit$gpredictions)) {
        Eyg <- signal.fit$gpredictions
        delta <- mean(abs(Eyg - Ey))
    } else {
        ## Use median of BB samples
        delta <- median(signal.delta)
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

    ## Compute delta0 for null samples
    null.delta <- apply(null.Eyg, 2, function(Eyg.b) {
        mean(abs(Eyg.b - mean(Eyg.b)))
    })

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Compute paired differences
    ## ========================================================================

    diff.delta <- signal.delta - null.delta

    ## ========================================================================
    ## Normality assessment and Box-Cox transformation
    ## ========================================================================

    if (verbose) {
        message("Assessing normality...")
    }

    ## Test normality of differences
    shapiro.raw <- tryCatch(
        shapiro.test(diff.delta),
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
    diff.for.test <- diff.delta

    if (apply.boxcox) {
        diff.for.bc <- diff.delta
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
            null.for.test <- null.delta
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

            signal.for.test <- signal.delta
            if (min(signal.for.test) <= 0) {
                signal.for.test <- signal.for.test + shift
            }
            if (abs(boxcox.lambda) > .Machine$double.eps) {
                signal.for.test <- (signal.for.test^boxcox.lambda - 1) / boxcox.lambda
            } else {
                signal.for.test <- log(signal.for.test)
            }
        } else {
            null.for.test <- null.delta
            signal.for.test <- signal.delta
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

        wilcox.result <- wilcox.test(diff.delta, mu = 0, alternative = "greater")
        p.value <- wilcox.result$p.value
    }

    log.p.value <- log(max(p.value, .Machine$double.eps))

    ## ========================================================================
    ## Compute normalized effect sizes
    ## ========================================================================

    ## Null distribution statistics
    null.mean <- mean(null.delta)
    null.sd <- sd(null.delta)
    null.median <- median(null.delta)
    null.mad <- mad(null.delta, constant = 1.4826)  # scaled to be consistent with SD

    ## Normalized delta: z-score relative to null
    ## How many SDs is the point estimate from the null mean?
    if (null.sd > .Machine$double.eps) {
        delta.z <- (delta - null.mean) / null.sd
    } else {
        delta.z <- NA
    }

    ## Robust normalized delta: using median and MAD
    ## More robust to outliers in null distribution
    if (null.mad > .Machine$double.eps) {
        delta.robust.z <- (delta - null.median) / null.mad
    } else {
        delta.robust.z <- NA
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
        abline(h = Ey, col = "darkgreen", lty = 3, lwd = 1.5)

        legend("topright",
               legend = c("Signal E(y|x)", "Null E(y|x)", "E(y)"),
               col = c("red", "blue", "darkgreen"),
               lty = c(1, 2, 3), lwd = c(2, 2, 1.5), cex = 0.8, bty = "n")

        ## Plot 2: Histogram of paired differences
        hist(diff.delta, breaks = 30, main = "Paired Differences",
             xlab = expression(delta[0]^signal - delta[0]^null),
             col = "lightblue", border = "white", las = 1)
        abline(v = 0, col = "red", lwd = 2, lty = 2)
        abline(v = mean(diff.delta), col = "darkblue", lwd = 2)

        ## Add p-value text
        text.x <- quantile(diff.delta, 0.95)
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
        delta = delta,
        delta.z = delta.z,
        delta.robust.z = delta.robust.z,
        p.value = p.value,
        log.p.value = log.p.value,
        test.type = test.type,
        signal.delta = signal.delta,
        null.delta = null.delta,
        diff.delta = diff.delta,
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
        Ey = Ey,
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

    class(output) <- "assoc0"

    return(output)
}


#' @rdname fassoc0.test
#' @export
zofam.test <- fassoc0.test


## ============================================================================
## S3 Methods
## ============================================================================

#' Coef Method for assoc0 Objects
#'
#' Extracts key coefficients from an assoc0 object.
#'
#' @param object An object of class "assoc0".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector containing delta, normalized effect sizes,
#'   and p-value.
#' @method coef assoc0
#' @export
coef.assoc0 <- function(object, ...) {
    c(delta = object$delta,
      delta.z = object$delta.z,
      delta.robust.z = object$delta.robust.z,
      p.value = object$p.value)
}


#' Print Method for assoc0 Objects
#'
#' @param x An object of class "assoc0".
#' @param digits Number of significant digits to display.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible x.
#' @method print assoc0
#' @export
print.assoc0 <- function(x, digits = 4, ...) {
    cat("\nZero-Order Functional Association Test (Paired BB)\n")
    cat("==================================================\n\n")

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
    cat("  Marginal mean E(y): ", signif(x$Ey, digits), "\n", sep = "")
    cat("  Association measure (delta0): ", signif(x$delta, digits), "\n", sep = "")
    cat("  Normalized effect size (z): ", signif(x$delta.z, digits), "\n", sep = "")
    cat("  Robust effect size (z.robust): ", signif(x$delta.robust.z, digits), "\n", sep = "")
    cat("  Mean difference (signal - null): ", signif(mean(x$diff.delta), digits), "\n", sep = "")
    cat("  p-value: ", signif(x$p.value, digits), "\n", sep = "")
    cat("  log(p-value): ", signif(x$log.p.value, digits), "\n", sep = "")

    invisible(x)
}


#' Summary Method for assoc0 Objects
#'
#' @param object An object of class "assoc0".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class "summary.assoc0" with summary statistics.
#' @method summary assoc0
#' @export
summary.assoc0 <- function(object, ...) {
    out <- list(
        call = object$call,
        n = length(object$x),
        test.type = object$test.type,
        n.BB = object$n.BB,
        degree = object$degree,
        Ey = object$Ey,
        delta = object$delta,
        delta.z = object$delta.z,
        delta.robust.z = object$delta.robust.z,
        p.value = object$p.value,
        null.mean = object$null.mean,
        null.sd = object$null.sd,
        null.median = object$null.median,
        null.mad = object$null.mad,
        boxcox.applied = object$boxcox.applied,
        boxcox.lambda = object$boxcox.lambda,
        shapiro.raw = object$shapiro.pvalue.raw,
        shapiro.bc = object$shapiro.pvalue.bc,
        signal.summary = summary(object$signal.delta),
        null.summary = summary(object$null.delta),
        diff.summary = summary(object$diff.delta),
        diff.quantiles = quantile(object$diff.delta,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    )

    class(out) <- "summary.assoc0"
    return(out)
}


#' @method print summary.assoc0
#' @export
print.summary.assoc0 <- function(x, digits = 4, ...) {
    cat("\nZero-Order Functional Association Test Summary\n")
    cat("==============================================\n\n")

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
    cat("  Marginal mean E(y): ", signif(x$Ey, digits), "\n", sep = "")
    cat("  Association measure (delta0): ", signif(x$delta, digits), "\n", sep = "")
    cat("  Normalized effect size (z): ", signif(x$delta.z, digits), "\n", sep = "")
    cat("  Robust effect size (z.robust): ", signif(x$delta.robust.z, digits), "\n", sep = "")
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


#' Plot Method for assoc0 Objects
#'
#' Creates diagnostic plots for objects of class "assoc0".
#'
#' @param x An object of class "assoc0".
#' @param type Character string specifying plot type. Options are:
#'   \describe{
#'     \item{"Exy"}{Conditional mean function with credible intervals (default)}
#'     \item{"diff"}{Histogram of signal - null differences}
#'     \item{"qq"}{QQ plot of differences}
#'     \item{"comparison"}{Side-by-side signal and null distributions}
#'   }
#' @param ... Additional arguments passed to plotting functions.
#'
#' @method plot assoc0
#' @export
plot.assoc0 <- function(x, type = c("Exy", "diff", "qq", "comparison"),
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
        abline(h = x$Ey, col = "darkgreen", lty = 3, lwd = 1.5)

    } else if (type == "diff") {
        ## Histogram of differences
        hist(x$diff.delta, breaks = 30, main = "Paired Differences",
             xlab = expression(delta[0]^signal - delta[0]^null),
             col = "lightblue", border = "white", las = 1, ...)
        abline(v = 0, col = "red", lwd = 2, lty = 2)
        abline(v = mean(x$diff.delta), col = "darkblue", lwd = 2)

    } else if (type == "qq") {
        ## QQ plot of differences
        qqnorm(x$diff.delta, main = "QQ Plot of Differences", las = 1, ...)
        qqline(x$diff.delta, col = "red")

    } else if (type == "comparison") {
        ## Side-by-side comparison
        op <- par(mfrow = c(1, 2))
        on.exit(par(op))

        ## Signal
        hist(x$signal.delta, breaks = 30, main = "Signal Distribution",
             xlab = expression(delta[0]^signal), col = "lightcoral",
             border = "white", las = 1)
        abline(v = x$delta, col = "red", lwd = 2)

        ## Null
        hist(x$null.delta, breaks = 30, main = "Null Distribution",
             xlab = expression(delta[0]^null), col = "lightblue",
             border = "white", las = 1)
        abline(v = mean(x$null.delta), col = "blue", lwd = 2)
    }

    invisible(NULL)
}


## ============================================================================
## Unified Interface
## ============================================================================

#' Functional Association Test
#'
#' Tests for the existence of non-trivial functional association between two
#' variables x and y using either zero-order or first-order functional
#' association measures. This function serves as a unified interface to both
#' fassoc0.test() and fassoc1.test().
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x).
#' @param order Integer specifying the order of the functional association test.
#'   Can be 0 (zero-order) or 1 (first-order). Default is 1.
#' @param test.type Character string specifying the test type. See
#'   \code{\link{fassoc0.test}} or \code{\link{fassoc1.test}} for options.
#' @param boxcox Character string controlling Box-Cox transformation.
#' @param ... Additional arguments passed to fassoc0.test() or fassoc1.test().
#'
#' @return An object of class "fassoc0" (if order=0) or "fassoc1" (if order=1).
#'
#' @details
#' This function provides a unified interface to test for functional associations
#' between variables. The zero-order measure (order=0) quantifies deviations of
#' the conditional mean from the marginal mean, while the first-order measure
#' (order=1) quantifies the variability in the derivative of the conditional mean.
#'
#' The mathematical notation for these measures:
#' \itemize{
#'   \item Zero-order: \eqn{\delta_0 = \int |E_x(y) - E(y)| dx}
#'   \item First-order: \eqn{\delta_1 = \int |dE_x(y)/dx| dx}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
#'
#' # Test for zero-order functional association
#' result0 <- fassoc.test(x, y, order = 0, n.BB = 500)
#'
#' # Test for first-order functional association
#' result1 <- fassoc.test(x, y, order = 1, n.BB = 500)
#'
#' # Compare p-values
#' print(result0$p.value)
#' print(result1$p.value)
#' }
#'
#' @seealso \code{\link{fassoc0.test}}, \code{\link{fassoc1.test}}
#' @export
fassoc.test <- function(x,
                        y,
                        order = 1,
                        test.type = c("paired.t", "weighted.pvalue", "wilcoxon"),
                        boxcox = c("auto", "always", "never"),
                        ...) {

    test.type <- match.arg(test.type)
    boxcox <- match.arg(boxcox)

    if (!is.numeric(order) || length(order) != 1 || !(order %in% c(0, 1))) {
        stop("'order' must be either 0 or 1")
    }

    if (order == 0) {
        fassoc0.test(x = x, y = y,
                     test.type = test.type,
                     boxcox = boxcox,
                     ...)
    } else {
        fassoc1.test(x = x, y = y,
                     test.type = test.type,
                     boxcox = boxcox,
                     ...)
    }
}


#' @rdname fassoc.test
#' @export
functional.association.test <- fassoc.test

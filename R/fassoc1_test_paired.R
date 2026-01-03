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

#' Tests for First-Order Functional Association Between Two Variables (Paired BB + Permutations; Hybrid Pairing)
#'
#' Computes and tests the first-order functional association measure between two
#' variables x and y. The statistic is based on the total absolute first difference
#' of the conditional mean curve E(y|x) over a uniform grid:
#'   delta1 = sum(|diff(E(y|x))|)
#'
#' This implementation uses paired Bayesian bootstrap (BB) weights for signal and
#' null computations, plus an explicit permutation null with n.perms permutations.
#'
#' Hybrid pairing rule (requested):
#'   - If test.type is diff-based (paired.t / wilcoxon) and n.perms <= n.BB:
#'       use unique BB columns (no cycling).
#'   - If test.type is diff-based and n.perms > n.BB:
#'       * if n.perms <= max.nBB.expand: expand n.BB to n.perms (regenerate lambda and recompute signal)
#'       * else: cycle BB columns and do inference on BB-cluster summaries (valid under cycling-induced dependence)
#'   - If test.type is weighted.pvalue: cycling is allowed (differences are not the inference target).
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values (same length as x). Can be binary or continuous.
#' @param test.type One of "paired.t", "weighted.pvalue", or "wilcoxon".
#' @param boxcox One of "auto", "always", or "never" controlling Box-Cox normalization.
#' @param boxcox.alpha Alpha for Shapiro-Wilk in boxcox="auto".
#' @param bw Bandwidth. If NULL, automatically selected by magelo.
#' @param n.perms Number of permutations for null distribution. Default 1000.
#' @param n.BB Number of BB samples for signal. Default 1000.
#' @param max.nBB.expand Threshold for expanding BB to match permutations in diff-based tests.
#'   If n.perms > n.BB and n.perms <= max.nBB.expand, n.BB is expanded to n.perms.
#'   Default 5000.
#' @param cluster.agg Aggregation for BB-cluster summaries when cycling is used:
#'   "mean", "trimmed.mean", or "median". Default "mean".
#' @param cluster.trim Trim proportion (0 to <0.5) for cluster.agg="trimmed.mean". Default 0.10.
#' @param grid.size Grid size for gpredictions. Default 400.
#' @param degree Local polynomial degree (1 or 2). Default 1.
#' @param min.K Minimum neighbors. Default 5.
#' @param n.cores Cores for permutation loop on Unix-alikes via parallel::mclapply. Default 1.
#' @param seed Optional RNG seed.
#' @param plot.it If TRUE, produces diagnostics plots. Default TRUE.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param verbose If TRUE, prints progress messages. Default TRUE.
#'
#' @return An object of class "assoc1" with signal/null summaries and test results.
#' @export
fassoc1.test <- function(x,
                         y,
                         test.type = c("paired.t", "weighted.pvalue", "wilcoxon"),
                         boxcox = c("auto", "always", "never"),
                         boxcox.alpha = 0.05,
                         bw = NULL,
                         n.perms = 1000,
                         n.BB = 1000,
                         max.nBB.expand = 5000,
                         cluster.agg = c("mean", "trimmed.mean", "median"),
                         cluster.trim = 0.10,
                         grid.size = 400,
                         degree = 1,
                         min.K = 5,
                         n.cores = 1,
                         seed = NULL,
                         plot.it = TRUE,
                         xlab = "x",
                         ylab = "y",
                         verbose = TRUE) {

    ## ------------------------------------------------------------------------
    ## Input validation
    ## ------------------------------------------------------------------------

    test.type <- match.arg(test.type)
    boxcox <- match.arg(boxcox)
    cluster.agg <- match.arg(cluster.agg)

    if (!is.numeric(x) || !is.numeric(y)) stop("Both 'x' and 'y' must be numeric vectors")
    if (length(x) != length(y)) stop("'x' and 'y' must have the same length")

    if (!is.null(bw) && (!is.numeric(bw) || length(bw) != 1 || bw <= 0)) {
        stop("'bw' must be NULL or a positive numeric scalar")
    }

    if (!is.numeric(n.perms) || length(n.perms) != 1 || n.perms < 100) stop("'n.perms' must be at least 100")
    if (!is.numeric(n.BB) || length(n.BB) != 1 || n.BB < 100) stop("'n.BB' must be at least 100")

    if (!is.numeric(max.nBB.expand) || length(max.nBB.expand) != 1 || max.nBB.expand < 100) {
        stop("'max.nBB.expand' must be at least 100")
    }

    if (!is.numeric(cluster.trim) || length(cluster.trim) != 1 || cluster.trim < 0 || cluster.trim >= 0.5) {
        stop("'cluster.trim' must be in [0, 0.5)")
    }

    if (!is.numeric(grid.size) || length(grid.size) != 1 || grid.size < 10) stop("'grid.size' must be at least 10")
    if (!is.numeric(degree) || length(degree) != 1 || !(degree %in% c(1, 2))) stop("'degree' must be 1 or 2")
    if (!is.numeric(min.K) || length(min.K) != 1 || min.K < 2) stop("'min.K' must be at least 2")
    if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1) stop("'n.cores' must be at least 1")

    if (!is.numeric(boxcox.alpha) || length(boxcox.alpha) != 1 || boxcox.alpha <= 0 || boxcox.alpha >= 1) {
        stop("'boxcox.alpha' must be between 0 and 1")
    }

    n.perms <- as.integer(n.perms)
    n.BB <- as.integer(n.BB)
    max.nBB.expand <- as.integer(max.nBB.expand)
    grid.size <- as.integer(grid.size)
    degree <- as.integer(degree)
    min.K <- as.integer(min.K)
    n.cores <- as.integer(n.cores)

    if (!is.null(seed)) {
        if (n.cores > 1) RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
    }

    ## Remove non-finite
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < length(x)) {
        warning(sprintf("Removed %d non-finite values", length(x) - sum(idx)))
        x <- x[idx]
        y <- y[idx]
    }

    nx <- length(x)
    if (nx <= 2 * min.K) stop(sprintf("Sample size (%d) must be greater than 2 * min.K (%d)", nx, 2 * min.K))

    uses.diffs <- test.type %in% c("paired.t", "wilcoxon")

    if (verbose) {
        routine.ptm <- proc.time()
        message("First-Order Functional Association Test (Paired BB + Permutations; Hybrid Pairing)")
        message(sprintf("  n = %d, n.BB = %d, n.perms = %d, test.type = %s", nx, n.BB, n.perms, test.type))
    }

    ## ------------------------------------------------------------------------
    ## Helper: generate Dirichlet weights
    ## ------------------------------------------------------------------------

    generate.lambda <- function(nx, n.BB) {
        if (exists("generate.dirichlet.weights")) {
            generate.dirichlet.weights(nx, n.BB)
        } else {
            ## Fallback: exponential normalization, scaled to sum nx
            lam <- matrix(0, nrow = nx, ncol = n.BB)
            for (b in seq_len(n.BB)) {
                e <- rexp(nx, rate = 1)
                lam[, b] <- e / sum(e) * nx
            }
            lam
        }
    }

    ## ------------------------------------------------------------------------
    ## Hybrid decision: expand or cluster if needed (diff-based tests only)
    ## ------------------------------------------------------------------------

    pairing.strategy <- "standard"
    expanded.nBB <- FALSE
    cycling.used <- FALSE
    do.cluster.inference <- FALSE

    if (uses.diffs && n.perms > n.BB) {
        if (n.perms <= max.nBB.expand) {
            pairing.strategy <- "expand"
            expanded.nBB <- TRUE
            if (verbose) {
                message(sprintf("  Hybrid rule: expanding n.BB from %d to %d (<= max.nBB.expand = %d).",
                                n.BB, n.perms, max.nBB.expand))
            }
            n.BB <- n.perms
        } else {
            pairing.strategy <- "cluster"
            cycling.used <- TRUE
            do.cluster.inference <- TRUE
            if (verbose) {
                message(sprintf("  Hybrid rule: n.perms (%d) > n.BB (%d) and > max.nBB.expand (%d); using BB-cluster inference (%s).",
                                n.perms, n.BB, max.nBB.expand, cluster.agg))
            }
        }
    }

    ## ------------------------------------------------------------------------
    ## Generate paired BB weights (lambda)
    ## ------------------------------------------------------------------------

    if (verbose) {
        message("Generating Dirichlet weights (paired BB)...")
        ptm <- proc.time()
    }

    lambda <- generate.lambda(nx, n.BB)

    if (verbose) message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

    ## ------------------------------------------------------------------------
    ## Signal fit (n.BB BB draws)
    ## ------------------------------------------------------------------------

    if (verbose) {
        message("Computing signal BB samples...")
        ptm <- proc.time()
    }

    if (exists("magelo.with.external.BB")) {
        signal.fit <- magelo.with.external.BB(
            x = x,
            y = y,
            lambda = lambda,
            bw = bw,
            grid.size = grid.size,
            degree = degree,
            min.K = min.K
        )
    } else {
        if (!exists("magelo")) stop("Neither magelo.with.external.BB nor magelo found")
        y.binary <- all(y %in% c(0, 1))
        warning("magelo.with.external.BB not found. Using magelo(); BB pairing is not guaranteed.")
        signal.fit <- magelo(
            x = x,
            y = y,
            bw = bw,
            grid.size = grid.size,
            degree = degree,
            min.K = min.K,
            y.binary = y.binary,
            n.BB = n.BB,
            get.BB.gpredictions = TRUE,
            get.gpredictions.CrI = FALSE,
            get.predictions.CrI = FALSE,
            get.BB.predictions = FALSE
        )
    }

    xgrid <- signal.fit$xgrid
    signal.Eyg <- signal.fit$BB.gpredictions
    bw.used <- signal.fit$opt.bw

    signal.delta1 <- apply(signal.Eyg, 2, function(z) sum(abs(diff(z))))

    if (!is.null(signal.fit$gpredictions)) {
        Eyg <- signal.fit$gpredictions
        delta1 <- sum(abs(diff(Eyg)))
        Delta1 <- Eyg[length(Eyg)] - Eyg[1]
    } else {
        delta1 <- stats::median(signal.delta1)
        Delta1 <- NA_real_
    }

    if (verbose) message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

    ## ------------------------------------------------------------------------
    ## Pairing indices for null permutations
    ## ------------------------------------------------------------------------

    if (uses.diffs) {
        if (!do.cluster.inference) {
            ## No cycling needed: ensure unique BB column for each permutation
            ## Here, n.BB >= n.perms is guaranteed (either original or expanded)
            bb.idx <- seq_len(n.perms)
            cycling.used <- FALSE
        } else {
            ## Cycling (n.perms > n.BB): will do BB-cluster inference
            bb.idx <- ((seq_len(n.perms) - 1L) %% n.BB) + 1L
        }
    } else {
        ## weighted.pvalue does not require independent paired differences
        if (n.perms <= n.BB) {
            bb.idx <- sample.int(n.BB, size = n.perms, replace = FALSE)
        } else {
            bb.idx <- ((seq_len(n.perms) - 1L) %% n.BB) + 1L
        }
    }

    ## ------------------------------------------------------------------------
    ## Compute permutation null
    ## ------------------------------------------------------------------------

    if (verbose) {
        message("Computing permutation null...")
        ptm <- proc.time()
    }

    perm.worker <- function(i) {
        b <- bb.idx[i]
        lambda.i <- lambda[, b, drop = FALSE]
        y.perm <- sample(y)

        if (exists("magelo.with.external.BB")) {
            null.fit.i <- magelo.with.external.BB(
                x = x,
                y = y.perm,
                lambda = lambda.i,
                bw = bw.used,
                grid.size = grid.size,
                degree = degree,
                min.K = min.K
            )
            eyg.i <- null.fit.i$BB.gpredictions[, 1]
        } else {
            ## magelo fallback
            y.binary <- all(y %in% c(0, 1))
            null.fit.i <- magelo(
                x = x,
                y = y.perm,
                bw = bw.used,
                grid.size = grid.size,
                degree = degree,
                min.K = min.K,
                y.binary = y.binary,
                n.BB = 1L,
                get.BB.gpredictions = TRUE,
                get.gpredictions.CrI = FALSE,
                get.predictions.CrI = FALSE,
                get.BB.predictions = FALSE
            )
            eyg.i <- null.fit.i$BB.gpredictions[, 1]
        }

        d1 <- sum(abs(diff(eyg.i)))
        list(null.delta1 = d1, null.Eyg = eyg.i, bb.col = b)
    }

    if (n.cores > 1 &&
        requireNamespace("parallel", quietly = TRUE) &&
        .Platform$OS.type == "unix") {
        res.list <- parallel::mclapply(seq_len(n.perms), perm.worker, mc.cores = n.cores)
    } else {
        if (n.cores > 1 && verbose) {
            message("  Note: parallel permutation requires Unix-alike + parallel::mclapply; running sequentially.")
        }
        res.list <- lapply(seq_len(n.perms), perm.worker)
    }

    null.delta1 <- vapply(res.list, function(z) z$null.delta1, numeric(1))
    null.Eyg <- do.call(cbind, lapply(res.list, function(z) z$null.Eyg))
    null.bb.col <- vapply(res.list, function(z) z$bb.col, integer(1))

    if (verbose) message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

    ## ------------------------------------------------------------------------
    ## Paired differences and (optional) BB-cluster summaries
    ## ------------------------------------------------------------------------

    diff.delta1 <- signal.delta1[null.bb.col] - null.delta1

    inference.vector <- NULL
    inference.unit <- NULL
    diff.delta1.cluster.stat <- NULL
    diff.delta1.cluster.n <- NULL

    if (uses.diffs) {
        if (!do.cluster.inference) {
            inference.vector <- diff.delta1
            inference.unit <- "paired.differences"
        } else {
            ## Cluster summaries by BB column
            agg.fun <- switch(cluster.agg,
                              "mean" = function(z) mean(z),
                              "median" = function(z) stats::median(z),
                              "trimmed.mean" = function(z) mean(z, trim = cluster.trim))

            diff.delta1.cluster.stat <- tapply(diff.delta1, INDEX = null.bb.col, FUN = agg.fun)
            diff.delta1.cluster.stat <- as.numeric(diff.delta1.cluster.stat)

            diff.delta1.cluster.n <- tapply(diff.delta1, INDEX = null.bb.col, FUN = length)
            diff.delta1.cluster.n <- as.integer(diff.delta1.cluster.n)

            inference.vector <- diff.delta1.cluster.stat
            inference.unit <- "bb.cluster.summaries"
        }
    }

    ## ------------------------------------------------------------------------
    ## Testing / transformations
    ## ------------------------------------------------------------------------

    boxcox.applied <- FALSE
    boxcox.lambda <- NA_real_
    shapiro.pvalue.raw <- NA_real_
    shapiro.pvalue.bc <- NA_real_

    p.value <- NA_real_

    if (test.type %in% c("paired.t", "wilcoxon")) {

        sh.raw <- tryCatch(stats::shapiro.test(inference.vector), error = function(e) list(p.value = NA_real_))
        shapiro.pvalue.raw <- sh.raw$p.value

        apply.boxcox <- FALSE
        if (boxcox == "always") apply.boxcox <- TRUE
        if (boxcox == "auto" && !is.na(shapiro.pvalue.raw) && shapiro.pvalue.raw < boxcox.alpha) apply.boxcox <- TRUE

        vec.for.test <- inference.vector

        if (apply.boxcox) {
            if (min(null.for.test) <= 0) {
                s <- stats::sd(null.for.test)
                eps <- ifelse(is.finite(s) && s > 0, 0.01 * s, 1e-8)
                shift <- abs(min(null.for.test)) + eps
                null.for.test <- null.for.test + shift
                signal.for.test <- signal.for.test + shift
            }

            if (exists("boxcox.mle")) {
                bc.fit <- tryCatch(boxcox.mle(null.for.test ~ 1), error = function(e) NULL)
                if (!is.null(bc.fit)) {
                    boxcox.lambda <- bc.fit$lambda
                    boxcox.applied <- TRUE

                    if (abs(boxcox.lambda) > .Machine$double.eps) {
                        null.for.test <- (null.for.test^boxcox.lambda - 1) / boxcox.lambda
                        signal.for.test <- (signal.for.test^boxcox.lambda - 1) / boxcox.lambda
                    } else {
                        null.for.test <- log(null.for.test)
                        signal.for.test <- log(signal.for.test)
                    }

                    sh.bc <- tryCatch(stats::shapiro.test(null.for.test), error = function(e) list(p.value = NA_real_))
                    shapiro.pvalue.bc <- sh.bc$p.value
                }
            } else {
                ## Fallback: log transform if positive
                if (all(null.for.test > 0)) {
                    null.for.test <- log(null.for.test)
                    signal.for.test <- log(signal.for.test)
                    boxcox.lambda <- 0
                    boxcox.applied <- TRUE

                    sh.bc <- tryCatch(stats::shapiro.test(null.for.test), error = function(e) list(p.value = NA_real_))
                    shapiro.pvalue.bc <- sh.bc$p.value
                }
            }
        }

        mu <- mean(null.for.test)
        sigma <- stats::sd(null.for.test)

        if (exists("weighted.p.value")) {
            p.value <- weighted.p.value(signal.for.test, mu, sigma, alternative = "greater")
        } else {
            ## Fallback: average upper-tail normal p-values across signal BB draws
            p.value <- mean(stats::pnorm(signal.for.test, mean = mu, sd = sigma, lower.tail = FALSE))
        }

        inference.unit <- "weighted.pvalue"
    }

    p.value <- max(p.value, .Machine$double.eps)
    log.p.value <- log(p.value)

    ## ------------------------------------------------------------------------
    ## Effect sizes (relative to permutation null)
    ## ------------------------------------------------------------------------

    null.mean <- mean(null.delta1)
    null.sd <- stats::sd(null.delta1)
    null.median <- stats::median(null.delta1)
    null.mad <- stats::mad(null.delta1, constant = 1.4826)

    delta1.z <- if (is.finite(null.sd) && null.sd > .Machine$double.eps) (delta1 - null.mean) / null.sd else NA_real_
    delta1.robust.z <- if (is.finite(null.mad) && null.mad > .Machine$double.eps) (delta1 - null.median) / null.mad else NA_real_

    ## ------------------------------------------------------------------------
    ## Diagnostic plots
    ## ------------------------------------------------------------------------

    if (plot.it) {
        op <- graphics::par(mfrow = c(1, 2), mar = c(3.5, 3.5, 2, 0.5),
                            mgp = c(2, 0.5, 0), tcl = -0.3)
        on.exit(graphics::par(op), add = TRUE)

        ## Plot 1: Conditional mean estimates (signal BB band vs permutation-null band)
        signal.mean <- rowMeans(signal.Eyg)
        signal.q <- apply(signal.Eyg, 1, stats::quantile, probs = c(0.025, 0.975))

        null.mean.curve <- rowMeans(null.Eyg)
        null.q <- apply(null.Eyg, 1, stats::quantile, probs = c(0.025, 0.975))

        ylim <- range(c(y, signal.q, null.q), na.rm = TRUE)

        graphics::plot(x, y, las = 1, ylim = ylim, xlab = xlab, ylab = ylab,
                       main = "Conditional Mean", col = "gray70", pch = 16, cex = 0.5)

        graphics::polygon(c(xgrid, rev(xgrid)),
                          c(null.q[1, ], rev(null.q[2, ])),
                          col = grDevices::adjustcolor("blue", alpha.f = 0.2), border = NA)
        graphics::polygon(c(xgrid, rev(xgrid)),
                          c(signal.q[1, ], rev(signal.q[2, ])),
                          col = grDevices::adjustcolor("red", alpha.f = 0.2), border = NA)

        graphics::lines(xgrid, signal.mean, col = "red", lwd = 2)
        graphics::lines(xgrid, null.mean.curve, col = "blue", lwd = 2, lty = 2)

        graphics::legend("topright",
                         legend = c("Signal E(y|x)", "Null E(y|x)"),
                         col = c("red", "blue"),
                         lty = c(1, 2), lwd = 2, cex = 0.8, bty = "n")

        ## Plot 2: Histogram of inference vector
        if (uses.diffs) {
            vec <- inference.vector
            if (!do.cluster.inference) {
                main.txt <- "Paired Differences"
                xlab.txt <- "signal - null"
            } else {
                main.txt <- sprintf("BB-Cluster Summaries (%s)", cluster.agg)
                xlab.txt <- "cluster(summary(signal - null))"
            }
        } else {
            vec <- null.delta1
            main.txt <- "Permutation Null (delta1)"
            xlab.txt <- "delta1 (null)"
        }

        graphics::hist(vec, breaks = 30, main = main.txt, xlab = xlab.txt,
                       col = "lightblue", border = "white", las = 1)
        graphics::abline(v = 0, col = "red", lwd = 2, lty = 2)
        graphics::abline(v = mean(vec), col = "darkblue", lwd = 2)

        tx <- stats::quantile(vec, 0.95)
        ty <- graphics::par("usr")[4] * 0.9
        graphics::text(tx, ty, labels = sprintf("p = %.4g", p.value),
                       adj = c(1, 1), cex = 0.9)
    }

    if (verbose) {
        message(sprintf("Done. p-value = %.4g", p.value))
        message(sprintf("Total time: %.2f sec", (proc.time() - routine.ptm)[3]))
    }

    ## ------------------------------------------------------------------------
    ## Output
    ## ------------------------------------------------------------------------

    out <- list(
        Delta1 = Delta1,
        delta1 = delta1,
        delta1.z = delta1.z,
        delta1.robust.z = delta1.robust.z,
        p.value = p.value,
        log.p.value = log.p.value,
        test.type = test.type,

        ## Hybrid pairing diagnostics
        pairing.strategy = pairing.strategy,
        expanded.nBB = expanded.nBB,
        cycling.used = cycling.used,
        do.cluster.inference = do.cluster.inference,
        cluster.agg = cluster.agg,
        cluster.trim = cluster.trim,

        ## Signal (BB)
        signal.delta1 = signal.delta1,
        signal.Eyg = signal.Eyg,

        ## Null (permutations)
        null.delta1 = null.delta1,
        null.Eyg = null.Eyg,
        null.bb.col = null.bb.col,

        ## Differences
        diff.delta1 = diff.delta1,
        diff.delta1.cluster.stat = diff.delta1.cluster.stat,
        diff.delta1.cluster.n = diff.delta1.cluster.n,

        ## Inference bookkeeping
        inference.unit = inference.unit,
        inference.n = if (is.null(inference.vector)) NA_integer_ else length(inference.vector),

        ## Null summary
        null.mean = null.mean,
        null.sd = null.sd,
        null.median = null.median,
        null.mad = null.mad,

        ## Box-Cox diagnostics
        boxcox.applied = boxcox.applied,
        boxcox.lambda = boxcox.lambda,
        shapiro.pvalue.raw = shapiro.pvalue.raw,
        shapiro.pvalue.bc = shapiro.pvalue.bc,

        ## Misc
        x = x,
        y = y,
        xgrid = xgrid,
        opt.bw = bw.used,
        grid.size = grid.size,
        degree = degree,
        min.K = min.K,
        n.BB = n.BB,
        n.perms = n.perms,
        max.nBB.expand = max.nBB.expand,
        lambda = lambda,
        call = match.call()
    )

    class(out) <- "assoc1"
    out
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

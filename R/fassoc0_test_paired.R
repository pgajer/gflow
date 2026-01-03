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
## ============================================================================

#' Zero-Order Functional Association Test with Paired Bayesian Bootstrap and Permutations (Hybrid Pairing)
#'
#' Tests for the existence of non-trivial *zero-order* functional association between
#' two variables $x$ and $y$ using the variability of the conditional mean curve
#' $E(y|x)$ relative to the marginal mean $E(y)$.
#'
#' The zero-order association measure is computed on a uniform grid as:
#' $$\delta_0 = \frac{1}{G}\sum_{g=1}^{G}\left|E(y|x_g) - E(y)\right|.$$
#'
#' This implementation uses paired Bayesian bootstrap (BB) weights for signal and
#' permutation-null samples and supports an explicit permutation layer `n.perms`.
#'
#' Hybrid paired-design rule (mirrors the latest paired `fassoc1.test()` logic):
#' - If `test.type` is diff-based ("paired.t" / "wilcoxon") and `n.perms > n.BB`:
#'   - Expand to `n.BB = n.perms` if `n.perms <= max.nBB.expand` (clean one-to-one pairing).
#'   - Otherwise cycle BB columns and perform inference on BB-cluster summaries
#'     (mean/median/trimmed mean) to address dependence induced by cycling.
#' - For `test.type = "weighted.pvalue"`, cycling is permitted (inference is not based on paired diffs).
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x).
#' @param test.type Character string specifying the test for inference:
#'   "paired.t", "weighted.pvalue", or "wilcoxon".
#' @param boxcox Character string controlling Box-Cox transformation:
#'   "auto", "always", or "never".
#' @param boxcox.alpha Alpha level for Shapiro-Wilk normality test when boxcox="auto".
#' @param bw Numeric bandwidth. If NULL (default), bandwidth is automatically selected.
#' @param n.perms Integer specifying the number of permutations for the null distribution. Default 1000.
#' @param n.BB Integer specifying the number of paired BB samples for the signal fit. Default 1000.
#' @param max.nBB.expand Integer threshold for expanding `n.BB` to match `n.perms` when `n.perms > n.BB`
#'   under diff-based inference. Default 5000.
#' @param cluster.agg Aggregation used to summarize BB-clusters when cycling occurs:
#'   "mean", "trimmed.mean", or "median". Default "mean".
#' @param cluster.trim Trim proportion for cluster.agg="trimmed.mean". Must be in [0, 0.5). Default 0.10.
#' @param grid.size Integer specifying the size of evaluation grid. Default 400.
#' @param degree Integer specifying the local polynomial degree (1 or 2). Default 1.
#' @param min.K Integer specifying minimum number of observations for local estimation. Default 5.
#' @param n.cores Integer specifying number of cores for permutation loop on Unix-alikes via mclapply. Default 1.
#' @param seed Optional integer RNG seed for reproducibility.
#' @param plot.it Logical indicating whether to produce diagnostic plots. Default TRUE.
#' @param xlab X-axis label. Default "x".
#' @param ylab Y-axis label. Default "y".
#' @param verbose Logical; if TRUE prints progress messages. Default TRUE.
#'
#' @return An object of class "assoc0" containing effect sizes, p-values,
#'   signal/null distributions, and diagnostic information.
#' @export
fassoc0.test <- function(x,
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

    ## Seed handling (parallel-safe on Unix if n.cores > 1)
    if (!is.null(seed)) {
        if (n.cores > 1) RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
    }

    ## Remove non-finite values
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
        message("Zero-Order Functional Association Test (Paired BB + Permutations; Hybrid Pairing)")
        message(sprintf("  n = %d, n.BB = %d, n.perms = %d, test.type = %s", nx, n.BB, n.perms, test.type))
    }

    ## ------------------------------------------------------------------------
    ## Helpers
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

    cluster.agg.fun <- function(z) {
        if (cluster.agg == "mean") return(mean(z))
        if (cluster.agg == "median") return(stats::median(z))
        mean(z, trim = cluster.trim)
    }

    ## ------------------------------------------------------------------------
    ## Hybrid decision: expand or cluster (diff-based inference only)
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
    ## Shared BB weights (lambda)
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

    ## Marginal mean (unweighted) and BB-weighted marginal means
    Ey <- mean(y)
    Ey.bb <- as.numeric(crossprod(lambda, y) / nx)

    ## Signal BB distribution of delta0
    signal.delta <- vapply(seq_len(n.BB), function(b) {
        mean(abs(signal.Eyg[, b] - Ey.bb[b]))
    }, numeric(1))

    ## Point estimate delta0 (using unweighted Ey)
    if (!is.null(signal.fit$gpredictions)) {
        Eyg <- signal.fit$gpredictions
        delta <- mean(abs(Eyg - Ey))
    } else {
        delta <- stats::median(signal.delta)
    }

    if (verbose) message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

    ## ------------------------------------------------------------------------
    ## Pairing indices for null permutations
    ## ------------------------------------------------------------------------

    if (uses.diffs) {
        if (!do.cluster.inference) {
            ## Ensure unique BB column for each permutation (n.BB >= n.perms by construction here)
            bb.idx <- seq_len(n.perms)
            cycling.used <- FALSE
        } else {
            ## Cycling (n.perms > n.BB): will do BB-cluster inference
            bb.idx <- ((seq_len(n.perms) - 1L) %% n.BB) + 1L
        }
    } else {
        ## weighted.pvalue branch does not require independent paired differences
        if (n.perms <= n.BB) {
            bb.idx <- sample.int(n.BB, size = n.perms, replace = FALSE)
        } else {
            bb.idx <- ((seq_len(n.perms) - 1L) %% n.BB) + 1L
        }
    }

    ## ------------------------------------------------------------------------
    ## Compute permutation null (n.perms permutations, paired via bb.idx)
    ## ------------------------------------------------------------------------

    if (verbose) {
        message("Computing permutation null...")
        ptm <- proc.time()
    }

    perm.worker <- function(i) {
        b <- bb.idx[i]
        lambda.i <- lambda[, b]
        y.perm <- sample(y)

        if (exists("magelo.with.external.BB")) {
            null.fit.i <- magelo.with.external.BB(
                x = x,
                y = y.perm,
                lambda = lambda[, b, drop = FALSE],
                bw = bw.used,
                grid.size = grid.size,
                degree = degree,
                min.K = min.K
            )
            null.Eyg.i <- null.fit.i$BB.gpredictions[, 1]
        } else {
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
            null.Eyg.i <- null.fit.i$BB.gpredictions[, 1]
        }

        Ey.null <- sum(lambda.i * y.perm) / nx
        null.delta.i <- mean(abs(null.Eyg.i - Ey.null))

        list(null.delta = null.delta.i, null.Eyg = null.Eyg.i, bb.col = b, Ey.null = Ey.null)
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

    null.delta <- vapply(res.list, function(z) z$null.delta, numeric(1))
    null.Eyg <- do.call(cbind, lapply(res.list, function(z) z$null.Eyg))
    null.bb.col <- vapply(res.list, function(z) z$bb.col, integer(1))
    null.Ey <- vapply(res.list, function(z) z$Ey.null, numeric(1))

    if (verbose) message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

    ## ------------------------------------------------------------------------
    ## Paired differences and (optional) BB-cluster summaries
    ## ------------------------------------------------------------------------

    diff.delta <- signal.delta[null.bb.col] - null.delta

    inference.vector <- NULL
    inference.unit <- NULL
    diff.delta.cluster.stat <- NULL
    diff.delta.cluster.n <- NULL

    if (uses.diffs) {
        if (!do.cluster.inference) {
            inference.vector <- diff.delta
            inference.unit <- "paired.differences"
        } else {
            diff.delta.cluster.stat <- tapply(diff.delta, INDEX = null.bb.col, FUN = cluster.agg.fun)
            diff.delta.cluster.stat <- as.numeric(diff.delta.cluster.stat)

            diff.delta.cluster.n <- tapply(diff.delta, INDEX = null.bb.col, FUN = length)
            diff.delta.cluster.n <- as.integer(diff.delta.cluster.n)

            inference.vector <- diff.delta.cluster.stat
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
            ## Shift to positive
            vec.bc <- vec.for.test
            if (min(vec.bc) <= 0) {
                s <- stats::sd(vec.bc)
                eps <- ifelse(is.finite(s) && s > 0, 0.01 * s, 1e-8)
                vec.bc <- vec.bc + abs(min(vec.bc)) + eps
            }

            if (exists("boxcox.mle")) {
                bc.fit <- tryCatch(boxcox.mle(vec.bc ~ 1), error = function(e) NULL)
                if (!is.null(bc.fit)) {
                    boxcox.lambda <- bc.fit$lambda
                    boxcox.applied <- TRUE
                    if (abs(boxcox.lambda) > .Machine$double.eps) {
                        vec.for.test <- (vec.bc^boxcox.lambda - 1) / boxcox.lambda
                    } else {
                        vec.for.test <- log(vec.bc)
                    }
                    sh.bc <- tryCatch(stats::shapiro.test(vec.for.test), error = function(e) list(p.value = NA_real_))
                    shapiro.pvalue.bc <- sh.bc$p.value
                }
            } else {
                if (all(vec.bc > 0)) {
                    vec.for.test <- log(vec.bc)
                    boxcox.lambda <- 0
                    boxcox.applied <- TRUE
                    sh.bc <- tryCatch(stats::shapiro.test(vec.for.test), error = function(e) list(p.value = NA_real_))
                    shapiro.pvalue.bc <- sh.bc$p.value
                }
            }
        }

        if (test.type == "paired.t") {
            t.res <- stats::t.test(vec.for.test, mu = 0, alternative = "greater")
            p.value <- t.res$p.value
        } else {
            w.res <- stats::wilcox.test(inference.vector, mu = 0, alternative = "greater")
            p.value <- w.res$p.value
        }

    } else if (test.type == "weighted.pvalue") {

        ## Box-Cox decision is based on the null distribution of delta0
        sh.raw <- tryCatch(stats::shapiro.test(null.delta), error = function(e) list(p.value = NA_real_))
        shapiro.pvalue.raw <- sh.raw$p.value

        apply.boxcox <- FALSE
        if (boxcox == "always") apply.boxcox <- TRUE
        if (boxcox == "auto" && !is.na(shapiro.pvalue.raw) && shapiro.pvalue.raw < boxcox.alpha) apply.boxcox <- TRUE

        null.for.test <- null.delta
        signal.for.test <- signal.delta

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

    null.mean <- mean(null.delta)
    null.sd <- stats::sd(null.delta)
    null.median <- stats::median(null.delta)
    null.mad <- stats::mad(null.delta, constant = 1.4826)

    delta.z <- if (is.finite(null.sd) && null.sd > .Machine$double.eps) (delta - null.mean) / null.sd else NA_real_
    delta.robust.z <- if (is.finite(null.mad) && null.mad > .Machine$double.eps) (delta - null.median) / null.mad else NA_real_

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

        ## Plot 2: Histogram of inference target
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
            vec <- null.delta
            main.txt <- "Permutation Null (delta0)"
            xlab.txt <- "delta0 (null)"
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
    ## Output (keep legacy field names where possible)
    ## ------------------------------------------------------------------------

    output <- list(
        delta = delta,
        delta.z = delta.z,
        delta.robust.z = delta.robust.z,
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
        max.nBB.expand = max.nBB.expand,

        ## Distributions
        signal.delta = signal.delta,
        null.delta = null.delta,
        diff.delta = diff.delta,
        diff.delta.cluster.stat = diff.delta.cluster.stat,
        diff.delta.cluster.n = diff.delta.cluster.n,

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

        ## Curves and means
        signal.Eyg = signal.Eyg,
        null.Eyg = null.Eyg,
        Ey = Ey,
        Ey.bb = Ey.bb,
        null.Ey = null.Ey,

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
        lambda = lambda,
        null.bb.col = null.bb.col,
        inference.unit = inference.unit,
        inference.n = if (is.null(inference.vector)) NA_integer_ else length(inference.vector),
        call = match.call()
    )

    class(output) <- "assoc0"
    output
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

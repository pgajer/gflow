#' Tests for Existence of Non-trivial Functional Association Between Two Variables
#'
#' Computes and tests the zero-order functional association measure between two
#' variables x and y. The function estimates the integral of |E_x(y) - E(y)|, where
#' E_x(y) is the conditional mean of y given x, and E(y) is the marginal mean of y.
#' Statistical significance is assessed using permutation tests with Bayesian bootstrap.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x). Can be binary or continuous.
#' @param two.sample.test Character string specifying the test for comparing distributions.
#'   Options are "Wasserstein", "KS", or "norm". Default is "norm".
#' @param bw Numeric bandwidth parameter for smoothing. If NULL (default), bandwidth is
#'   automatically selected.
#' @param n.perms Integer specifying the number of permutations for null distribution.
#'   Default is 10000.
#' @param n.BB Integer specifying the number of Bayesian bootstrap samples. Default is 1000.
#' @param grid.size Integer specifying the size of evaluation grid. Default is 400.
#' @param min.K Integer specifying minimum number of observations for local estimation.
#'   Default is 5.
#' @param n.cores Integer specifying the number of CPU cores for parallel computation.
#'   Default is 7.
#' @param plot.it Logical indicating whether to produce diagnostic plots. Default is TRUE.
#' @param xlab Character string for x-axis label. Default is "x".
#' @param ylab Character string for y-axis label. Default is "y".
#' @param Eyg.col Color specification for conditional mean curve. Default is "red".
#' @param pt.col Color specification for data points. Default is "black".
#' @param pt.pch Integer or character specifying point type. Default is 1.
#' @param CrI.as.polygon Logical indicating whether to draw credible intervals as
#'   polygons (TRUE) or lines (FALSE). Default is TRUE.
#' @param CrI.polygon.col Color specification for credible interval polygon.
#'   Default is "gray90".
#' @param null.lines.col Color specification for null distribution lines.
#'   Default is "gray95".
#' @param CrI.line.col Color specification for credible interval lines.
#'   Default is "gray".
#' @param CrI.line.lty Line type for credible interval lines. Default is 2.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return An object of class "assoc0" containing:
#'   \item{delta}{The observed functional association measure}
#'   \item{p.value}{The p-value for the test of association}
#'   \item{log.p.value}{Natural logarithm of the p-value}
#'   \item{null.delta}{Vector of null distribution values}
#'   \item{null.delta.boxcox}{Box-Cox transformed null distribution (if applicable)}
#'   \item{BB.delta}{Bayesian bootstrap samples of delta}
#'   \item{BB.delta.boxcox}{Box-Cox transformed bootstrap samples (if applicable)}
#'   \item{delta.boxcox}{Box-Cox transformed observed delta (if applicable)}
#'   \item{Ey.res}{Results from conditional mean estimation}
#'   \item{x}{Input predictor values}
#'   \item{y}{Input response values}
#'   \item{Ey}{Marginal mean of y}
#'   \item{xg}{Grid points for evaluation}
#'   \item{Exyg}{Conditional mean estimates at grid points}
#'
#' @details
#' The function tests for functional association by comparing the observed measure
#' of association against a null distribution generated through permutation. The
#' conditional mean function E_x(y) is estimated using local polynomial methods
#' with Bayesian bootstrap for uncertainty quantification.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
#'
#' # Test for functional association
#' result <- fassoc0.test(x, y, n.cores = 2, n.perms = 1000)
#'
#' # Custom plot
#' plot(result, plot = "Exy")
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom stats ks.test sd
#' @importFrom graphics plot points lines abline polygon hist box text par matlines
#' @export
fassoc0.test <- function(x,
                        y,
                        two.sample.test = c("norm", "Wasserstein", "KS"),
                        bw = NULL,
                        n.perms = 10000,
                        n.BB = 1000,
                        grid.size = 400,
                        min.K = 5,
                        n.cores = 7,
                        plot.it = TRUE,
                        xlab = "x",
                        ylab = "y",
                        Eyg.col = "red",
                        pt.col = "black",
                        pt.pch = 1,
                        CrI.as.polygon = TRUE,
                        CrI.polygon.col = "gray90",
                        null.lines.col = "gray95",
                        CrI.line.col = "gray",
                        CrI.line.lty = 2,
                        verbose = TRUE) {

    # Input validation
    two.sample.test <- match.arg(two.sample.test)

    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Both 'x' and 'y' must be numeric vectors")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    # Additional parameter validation
    if (!is.null(bw) && (!is.numeric(bw) || bw <= 0)) {
        stop("'bw' must be NULL or a positive numeric value")
    }

    if (!is.numeric(n.perms) || n.perms < 100) {
        stop("'n.perms' must be at least 100")
    }

    if (!is.numeric(n.BB) || n.BB < 100) {
        stop("'n.BB' must be at least 100")
    }

    if (!is.numeric(grid.size) || grid.size < 10) {
        stop("'grid.size' must be at least 10")
    }

    if (!is.numeric(min.K) || min.K < 2) {
        stop("'min.K' must be at least 2")
    }

    if (!is.numeric(n.cores) || n.cores < 1) {
        stop("'n.cores' must be at least 1")
    }

    # Convert integers to ensure proper types
    n.perms <- as.integer(n.perms)
    n.BB <- as.integer(n.BB)
    grid.size <- as.integer(grid.size)
    min.K <- as.integer(min.K)
    n.cores <- as.integer(n.cores)

    if (verbose) {
        routine.ptm <- proc.time()
    }

    n <- length(x)

    y.binary <- FALSE
    if (all(y %in% c(0, 1))) {
        y.binary <- TRUE
    }

    # Remove non-finite values
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < length(x)) {
        warning(sprintf("Removed %d non-finite values", length(x) - sum(idx)))
    }
    x <- x[idx]
    y <- y[idx]
    n <- length(x)

    if (n <= 2 * min.K) {
        stop(sprintf("Sample size (%d) must be greater than 2 * min.K (%d)", n, 2 * min.K))
    }

    if (verbose) {
        cat("Estimating the zero-order functional association measure for x and y ... ")
        ptm <- proc.time()
    }

    # Check if required functions exist
    if (!exists("magelo")) {
        stop("Required function 'magelo' not found")
    }

    # Estimate conditional mean function
    Ey.res <- magelo(x, y,
                     bw = bw,
                     grid.size = grid.size,
                     y.binary = y.binary,
                     min.K = min.K,
                     n.BB = n.BB,
                     get.BB.predictions = FALSE,
                     get.predictions.CrI = FALSE,
                     get.BB.gpredictions = TRUE,
                     get.gpredictions.CrI = TRUE)

    xg <- Ey.res$xgrid
    Exyg <- Ey.res$gpredictions
    BB.Exyg <- Ey.res$BB.gpredictions
    bw <- Ey.res$opt.bw

    Ey <- mean(y)
    delta <- mean(abs(Exyg - Ey))
    BB.delta <- apply(BB.Exyg, 2, function(x) mean(abs(x - mean(x))))

    if (verbose) {
        elapsed.time(ptm)
        cat("Generating the null distribution of the zero-order functional association measure ... ")
        ptm <- proc.time()
    }

    # Generate null distribution
    if (n.cores == 1) {
        null.ExYg <- matrix(nrow = grid.size, ncol = n.perms)
        for (i in seq_len(n.perms)) {
            y.rand <- sample(y)
            r <- magelo(x, y.rand, bw = bw, n.BB = 1, get.BB.gpredictions = TRUE)
            null.ExYg[, i] <- r$BB.gpredictions[, 1]
        }
    } else {
        # Parallel computation
        get.BB.gpredictions <- function(i) {
            y.rand <- sample(y)
            r <- magelo(x, y.rand, bw = bw, n.BB = 1, get.BB.gpredictions = TRUE)
            return(r$BB.gpredictions[, 1])
        }

        # Use mclapply for parallel processing
        results <- mclapply(seq_len(n.perms), get.BB.gpredictions, mc.cores = n.cores)

        # Check for errors in parallel results
        failed <- sapply(results, inherits, "try-error")
        if (any(failed)) {
            warning(sprintf("%d permutations failed in parallel computation", sum(failed)))
            results <- results[!failed]
        }

        null.ExYg <- do.call(cbind, results)
    }

    null.delta <- apply(null.ExYg, 2, function(x) mean(abs(x - mean(x))))

    # Initialize Box-Cox variables
    null.delta.boxcox <- NULL
    delta.boxcox <- NULL
    BB.delta.boxcox <- NULL
    lambda <- NULL

    # Perform hypothesis test
    if (two.sample.test == "Wasserstein") {
        if (!exists("wasserstein1d.test")) {
            stop("Function 'wasserstein1d.test' not found")
        }
        dW.res <- wasserstein1d.test(BB.delta, null.delta, n.perms = n.perms, n.cores = n.cores)
        p.value <- dW.res$p.value
        log.p.value <- log(p.value)

    } else if (two.sample.test == "KS") {
        KS.res <- ks.test(BB.delta, null.delta)
        p.value <- KS.res$p.value
        log.p.value <- log(p.value)

    } else {  # "norm" test
        # Box-Cox transformation
        ## Find optimal lambda
        fit <- boxcox.mle(null.delta ~ 1)
        lambda <- fit$lambda

        # Apply Box-Cox transform
        if (abs(lambda) > .Machine$double.eps) {
            null.delta.boxcox <- (null.delta^lambda - 1) / lambda
            delta.boxcox <- (delta^lambda - 1) / lambda
            BB.delta.boxcox <- (BB.delta^lambda - 1) / lambda
        } else {
            null.delta.boxcox <- log(null.delta)
            delta.boxcox <- log(delta)
            BB.delta.boxcox <- log(BB.delta)
        }

        mu <- mean(null.delta.boxcox)
        sigma <- sd(null.delta.boxcox)

        if (!exists("weighted.p.value")) {
            # Simple p-value calculation if weighted.p.value not available
            p.value <- mean(null.delta.boxcox >= mean(BB.delta.boxcox))
            p.value <- max(p.value, 1/n.perms)  # Avoid p-value of 0
        } else {
            p.value <- weighted.p.value(BB.delta.boxcox, mu, sigma)
        }
        log.p.value <- log(p.value)
    }

    # Create diagnostic plots
    if (plot.it) {
        op <- par(mfrow = c(1, 2), mar = c(3.75, 3.75, 0.5, 0.5),
                  mgp = c(2, 0.5, 0), tcl = -0.3)
        on.exit(par(op), add = TRUE)

        # Plot 1: Dependence of y on x
        yu <- Ey.res$Eyg.CrI[1, ]
        yl <- Ey.res$Eyg.CrI[2, ]
        ylim <- range(c(Ey.res$y, Ey.res$yg, yl, yu), na.rm = TRUE)

        plot(Ey.res$x, Ey.res$y, las = 1, ylim = ylim, xlab = xlab,
             ylab = ylab, type = "n")

        if (CrI.as.polygon) {
            polygon(c(Ey.res$xg, rev(Ey.res$xg)), c(yl, rev(yu)),
                   col = CrI.polygon.col, border = NA)
        } else {
            matlines(Ey.res$xg, cbind(yl, yu), col = CrI.line.col,
                    lty = CrI.line.lty)
        }

        # Draw sample of null lines
        n.lines <- min(100, ncol(null.ExYg))
        for (i in seq_len(n.lines)) {
            lines(Ey.res$xg, null.ExYg[, i], col = null.lines.col)
        }

        abline(h = mean(y), col = 'blue')

        points(Ey.res$x, Ey.res$y, col = pt.col, pch = pt.pch)
        lines(Ey.res$xg, Ey.res$Eyg, col = Eyg.col)

        # Plot 2: Histogram of null distribution
        if (exists("hist2")) {
            h <- hist2(null.delta, BB.delta, x1.lab = "null.fassoc0",
                      x2.lab = "BB.fassoc0", n.x1.breaks = 15, n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            hist(null.delta, breaks = 30, main = "", xlab = "null.fassoc0",
                 col = "lightgray", border = "white")
            xlim_used <- par("usr")[1:2]
            ylim_used <- par("usr")[3:4]
        }

        abline(v = delta, col = 'red', lwd = 2)

        if (exists("mode.1D")) {
            abline(v = mode.1D(null.delta), col = "blue")
        }

        # Add p-value text
        text(x = mean(xlim_used), y = 0.7 * diff(ylim_used) + ylim_used[1],
             labels = paste0("p-val: ", signif(p.value, digits = 2)))
        box()
    }

    if (verbose) {
        elapsed.time(ptm)
        txt <- sprintf("Total elapsed time")
        elapsed.time(routine.ptm, txt, with.brackets = FALSE)
    }

    # Prepare output
    output <- list(
        delta = delta,
        log.p.value = log.p.value,
        p.value = p.value,
        null.delta = null.delta,
        null.delta.boxcox = null.delta.boxcox,
        null.Ey = null.ExYg,
        BB.Ey = Ey.res$BB.gpredictions,
        BB.delta = BB.delta,
        BB.delta.boxcox = BB.delta.boxcox,
        delta.boxcox = delta.boxcox,
        lambda = lambda,
        Ey.res = Ey.res,
        x = x,
        y = y,
        Ey = Ey,
        xg = xg,
        Exyg = Exyg,
        two.sample.test = two.sample.test,
        n.perms = n.perms,
        n.BB = n.BB,
        call = match.call()
    )

    class(output) <- "assoc0"

    return(output)
}

#' @rdname fassoc0.test
#' @export
zofam.test <- fassoc0.test


#' Plot Method for assoc0 Objects
#'
#' Creates diagnostic plots for objects of class "assoc0" produced by fassoc0.test().
#'
#' @param x An object of class "assoc0" from fassoc0.test().
#' @param plot Character string specifying plot type. Options are:
#'   \describe{
#'     \item{"Exy"}{Conditional mean function with credible intervals}
#'     \item{"d1hist"}{Histogram of null distribution with observed value}
#'     \item{"bc.d1hist"}{Histogram of Box-Cox transformed null distribution}
#'   }
#' @param xlab Character string for x-axis label. Default is "x".
#' @param ylab Character string for y-axis label. Default is "y".
#' @param Eyg.col Color specification for conditional mean curve. Default is "blue".
#' @param Eyg.lwd Line width for conditional mean curve. Default is 2.
#' @param pt.col Color specification for data points. Default is "black".
#' @param pt.pch Integer or character specifying point type. Default is 1.
#' @param CrI.as.polygon Logical indicating whether to draw credible intervals as
#'   polygons (TRUE) or lines (FALSE). Default is TRUE.
#' @param CrI.polygon.col Color specification for credible interval polygon.
#'   Default is "gray90".
#' @param null.line.col Color specification for null hypothesis line.
#'   Default is "gray70".
#' @param null.CrI.col Color specification for null distribution credible intervals.
#'   Default is "gray90".
#' @param CrI.line.col Color specification for credible interval lines.
#'   Default is "gray".
#' @param CrI.line.lty Line type for credible interval lines. Default is 2.
#' @param dEy.CrI.col Color specification for conditional mean credible intervals.
#'   Default is "cornflowerblue".
#' @param d1hist.x Numeric x-coordinate for p-value label in histogram plots.
#'   If NULL, automatically positioned.
#' @param d1hist.y Numeric y-coordinate for p-value label in histogram plots.
#'   If NULL, automatically positioned.
#' @param title Character string for plot title. Default is "".
#' @param ylim Numeric vector of length 2 giving y-axis limits. If NULL,
#'   automatically determined.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible NULL. Called for side effect of creating plot.
#'
#' @examples
#' \dontrun{
#' # Generate example data and run test
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
#' result <- fassoc0.test(x, y, n.cores = 2, n.perms = 1000, plot.it = FALSE)
#'
#' # Create different plots
#' plot(result, plot = "Exy")
#' plot(result, plot = "d1hist")
#' plot(result, plot = "bc.d1hist")
#' }
#'
#' @method plot assoc0
#' @export
plot.assoc0 <- function(x,
                       plot = c("Exy", "d1hist", "bc.d1hist"),
                       xlab = "x",
                       ylab = "y",
                       Eyg.col = "blue",
                       Eyg.lwd = 2,
                       pt.col = "black",
                       pt.pch = 1,
                       CrI.as.polygon = TRUE,
                       CrI.polygon.col = "gray90",
                       null.line.col = "gray70",
                       null.CrI.col = "gray90",
                       CrI.line.col = "gray",
                       CrI.line.lty = 2,
                       dEy.CrI.col = "cornflowerblue",
                       d1hist.x = NULL,
                       d1hist.y = NULL,
                       title = "",
                       ylim = NULL,
                       ...) {

    # Rename first argument for clarity
    res <- x

    # Validate input
    if (!inherits(res, "assoc0")) {
        stop("'x' must be an object of class 'assoc0'")
    }

    plot <- match.arg(plot)

    # Extract results
    Ey.res <- res$Ey.res
    xg <- Ey.res$xg
    ng <- length(xg)
    Ey <- Ey.res$Eyg

    if (plot == "Exy") {
        # Dependence of y on x
        null.Ey.CrI <- apply(res$null.Ey, 1, quantile, probs = c(0.025, 0.975))
        null.Ey.u <- null.Ey.CrI[2, ]
        null.Ey.d <- null.Ey.CrI[1, ]

        # Check if plot.magelo exists, otherwise use base plot
        if (exists("plot.magelo")) {
            plot(Ey.res, xlab = xlab, ylab = ylab, Eyg.col = Eyg.col,
                 title = title, ylim = ylim)
        } else {
            # Fallback plotting
            if (is.null(ylim)) {
                ylim <- range(c(Ey.res$y, Ey.res$Eyg, Ey.res$Eyg.CrI), na.rm = TRUE)
            }

            plot(Ey.res$x, Ey.res$y, xlab = xlab, ylab = ylab,
                 main = title, ylim = ylim, las = 1)

            if (CrI.as.polygon && !is.null(Ey.res$Eyg.CrI)) {
                yu <- Ey.res$Eyg.CrI[2, ]
                yl <- Ey.res$Eyg.CrI[1, ]
                polygon(c(xg, rev(xg)), c(yl, rev(yu)),
                       col = CrI.polygon.col, border = NA)
            }
        }

        # Add null distribution credible interval
        polygon(c(xg, rev(xg)), c(null.Ey.d, rev(null.Ey.u)),
               col = null.CrI.col, border = NA)

        # Add mean line
        abline(h = mean(res$null.Ey), col = null.line.col)

        # Add conditional mean estimate
        lines(xg, Ey.res$Eyg, col = Eyg.col, lwd = Eyg.lwd)

        # Add credible interval lines if available
        if (!is.null(Ey.res$Eyg.CrI)) {
            lines(xg, Ey.res$Eyg.CrI[1, ], col = dEy.CrI.col)
            lines(xg, Ey.res$Eyg.CrI[2, ], col = dEy.CrI.col)
        }

    } else if (plot == "d1hist") {
        # Histogram of null delta with observed delta
        if (exists("hist2")) {
            h <- hist2(res$null.delta, res$BB.delta,
                      x1.lab = "null.fassoc0",
                      x2.lab = "BB.fassoc0",
                      n.x1.breaks = 15,
                      n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            h <- hist(res$null.delta, breaks = 30, main = title,
                     xlab = "null.fassoc0", las = 1)
            xlim_used <- range(h$breaks)
            ylim_used <- c(0, max(h$counts))
        }

        # Position p-value text
        if (is.null(d1hist.x)) {
            d1hist.x <- mean(xlim_used)
        }
        if (is.null(d1hist.y)) {
            d1hist.y <- 0.7 * diff(ylim_used) + ylim_used[1]
        }

        # Add reference lines
        abline(v = res$delta, col = 'red', lwd = 2)

        if (exists("mode.1D")) {
            abline(v = mode.1D(res$null.delta), col = "blue")
        }

        # Add p-value text
        text(x = d1hist.x, y = d1hist.y,
             labels = paste0("p-val: ", signif(res$p.value, digits = 2)))

    } else if (plot == "bc.d1hist") {
        # Box-Cox transformed histogram
        if (is.null(res$null.delta.boxcox)) {
            stop("Box-Cox transformed values not available. ",
                 "Run fassoc0.test with two.sample.test = 'norm'")
        }

        if (exists("hist2")) {
            h <- hist2(res$null.delta.boxcox, res$BB.delta.boxcox,
                      x1.lab = "BC(null.fassoc0)",
                      x2.lab = "BC(BB.fassoc0)",
                      n.x1.breaks = 15,
                      n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            h <- hist(res$null.delta.boxcox, breaks = 30, main = title,
                     xlab = "BC(null.fassoc0)", las = 1)
            xlim_used <- range(h$breaks)
            ylim_used <- c(0, max(h$counts))
        }

        # Position p-value text
        if (is.null(d1hist.x)) {
            d1hist.x <- mean(xlim_used)
        }
        if (is.null(d1hist.y)) {
            d1hist.y <- 0.7 * diff(ylim_used) + ylim_used[1]
        }

        # Add reference lines
        abline(v = res$delta.boxcox, col = 'red', lwd = 2)

        if (exists("mode.1D")) {
            abline(v = mode.1D(res$null.delta.boxcox), col = "blue")
        }

        # Add p-value text
        text(x = d1hist.x, y = d1hist.y,
             labels = paste0("p-val: ", signif(res$p.value, digits = 2)))
    }

    invisible(NULL)
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
    cat("\nZero-Order Functional Association Test\n")
    cat("=====================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Test Summary:\n")
    cat("  Sample size: ", length(x$x), "\n", sep = "")
    cat("  Test type: ", x$two.sample.test, "\n", sep = "")
    cat("  Permutations: ", x$n.perms, "\n", sep = "")
    cat("  Bootstrap samples: ", x$n.BB, "\n\n", sep = "")

    cat("Results:\n")
    cat("  Association measure (delta): ", signif(x$delta, digits), "\n", sep = "")
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
        test.type = object$two.sample.test,
        n.perms = object$n.perms,
        n.BB = object$n.BB,
        delta = object$delta,
        p.value = object$p.value,
        null.mean = mean(object$null.delta),
        null.sd = stats::sd(object$null.delta),
        null.quantiles = stats::quantile(object$null.delta, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    )

    class(out) <- "summary.assoc0"
    return(out)
}

#' @method print summary.assoc0
#' @export
print.summary.assoc0 <- function(x, digits = 4, ...) {
    cat("\nZero-Order Functional Association Test\n")
    cat("=====================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Test Configuration:\n")
    cat("  Sample size: ", x$n, "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n", sep = "")
    cat("  Permutations: ", x$n.perms, "\n", sep = "")
    cat("  Bootstrap samples: ", x$n.BB, "\n\n", sep = "")

    cat("Results:\n")
    cat("  Association measure (delta): ", signif(x$delta, digits), "\n", sep = "")
    cat("  p-value: ", signif(x$p.value, digits), "\n\n", sep = "")

    cat("Null Distribution Summary:\n")
    cat("  Mean: ", signif(x$null.mean, digits), "\n", sep = "")
    cat("  SD: ", signif(x$null.sd, digits), "\n", sep = "")
    cat("  Quantiles:\n")
    print(signif(x$null.quantiles, digits))

    invisible(x)
}

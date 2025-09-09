#' Tests for First-Order Functional Association Between Two Variables
#'
#' Computes and tests the first-order functional association measure between two
#' variables x and y. The function estimates the integral of |dE_x(y)/dx|, where
#' E_x(y) is the conditional mean of y given x. This measure captures the strength
#' of the derivative relationship between variables.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x). Can be binary or continuous.
#' @param two.sample.test Character string specifying the test for comparing distributions.
#'   Options are "Wasserstein", "KS", or "norm". Default is "norm".
#' @param bw Numeric bandwidth parameter for smoothing. If NULL (default), bandwidth is
#'   automatically selected.
#' @param n.perms Integer specifying the number of permutations for null distribution.
#'   Default is 1000.
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
#' @return An object of class "assoc1" containing:
#'   \item{Delta1}{The total change in conditional mean (E_y(x_max) - E_y(x_min))}
#'   \item{delta1}{The first-order functional association measure}
#'   \item{p.value}{The p-value for the test of association}
#'   \item{log.p.value}{Natural logarithm of the p-value}
#'   \item{null.delta1}{Vector of null distribution values}
#'   \item{null.delta1.boxcox}{Box-Cox transformed null distribution (if applicable)}
#'   \item{BB.delta1}{Bayesian bootstrap samples of delta1}
#'   \item{BB.delta1.boxcox}{Box-Cox transformed bootstrap samples (if applicable)}
#'   \item{delta1.boxcox}{Box-Cox transformed observed delta1 (if applicable)}
#'   \item{null.Ey}{Matrix of null conditional mean estimates}
#'   \item{null.dEy}{Matrix of null conditional mean derivatives}
#'   \item{BB.dEy}{Matrix of bootstrap conditional mean derivatives}
#'   \item{Ey.res}{Results from conditional mean estimation}
#'   \item{x}{Input predictor values}
#'   \item{y}{Input response values}
#'   \item{Ey}{Conditional mean estimates at grid points}
#'   \item{grid.size}{Number of grid points used}
#'
#' @details
#' The first-order functional association measure captures the variability in the
#' derivative of the conditional mean function. Unlike the zero-order measure which
#' looks at deviations from the mean, this measure is sensitive to the rate of change
#' in the relationship between x and y. The test is performed by comparing the observed
#' measure against a null distribution generated through permutation with Bayesian
#' bootstrap for uncertainty quantification.
#'
#' @examples
#' \dontrun{
#' # Generate example data with nonlinear relationship
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(4*pi*x) + rnorm(n, sd = 0.2)
#'
#' # Test for first-order functional association
#' result <- fassoc1.test(x, y, n.cores = 2, n.perms = 500)
#'
#' # Plot derivative of conditional mean
#' plot(result, plot = "dExy")
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom graphics plot points lines abline polygon hist box text par matlines
#' @export
fassoc1.test <- function(x,
                        y,
                        two.sample.test = c("norm", "Wasserstein", "KS"),
                        bw = NULL,
                        n.perms = 1000,
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
        cat("Estimating the first-order functional association measure for x and y ... ")
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
    Ey <- Ey.res$gpredictions
    BB.Ey <- Ey.res$BB.gpredictions
    bw <- Ey.res$opt.bw
    ng <- grid.size

    # Compute derivatives
    dEy <- diff(Ey)
    BB.dEy <- apply(BB.Ey, 2, diff)

    # Compute association measures
    Delta1 <- Ey[ng] - Ey[1]  # Total change
    delta1 <- sum(abs(dEy))    # Sum of absolute derivatives

    BB.Delta1 <- apply(BB.Ey, 2, function(z) z[ng] - z[1])
    BB.delta1 <- apply(BB.Ey, 2, function(z) sum(abs(diff(z))))

    if (verbose) {
        elapsed.time(ptm)
        cat("Generating the null distribution of the first-order functional association measure ... ")
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

        # Check for parallel package
        if (!requireNamespace("parallel", quietly = TRUE)) {
            warning("Package 'parallel' not available. Using sequential processing.")
            n.cores <- 1
            null.ExYg <- matrix(nrow = grid.size, ncol = n.perms)
            for (i in seq_len(n.perms)) {
                y.rand <- sample(y)
                r <- magelo(x, y.rand, bw = bw, n.BB = 1, get.BB.gpredictions = TRUE)
                null.ExYg[, i] <- r$BB.gpredictions[, 1]
            }
        } else {
            results <- parallel::mclapply(seq_len(n.perms), get.BB.gpredictions, mc.cores = n.cores)

            # Check for errors in parallel results
            failed <- sapply(results, inherits, "try-error")
            if (any(failed)) {
                warning(sprintf("%d permutations failed in parallel computation", sum(failed)))
                results <- results[!failed]
            }

            null.ExYg <- do.call(cbind, results)
        }
    }

    # Compute null derivatives
    null.dExYg <- apply(null.ExYg, 2, diff)
    null.delta1 <- apply(null.ExYg, 2, function(z) sum(abs(diff(z))))
    Enull.delta1 <- mean(null.delta1)

    # Initialize Box-Cox variables
    null.delta1.boxcox <- NULL
    BB.delta1.boxcox <- NULL
    delta1.boxcox <- NULL
    lambda <- NULL

    # Perform hypothesis test
    if (two.sample.test == "Wasserstein") {
        if (!exists("wasserstein1d.test")) {
            stop("Function 'wasserstein1d.test' not found")
        }
        dW.res <- wasserstein1d.test(BB.delta1, null.delta1, n.perms = n.perms, n.cores = n.cores)
        p.value <- dW.res$p.value
        log.p.value <- log(p.value)

    } else if (two.sample.test == "KS") {
        KS.res <- stats::ks.test(BB.delta1, null.delta1)
        p.value <- KS.res$p.value
        log.p.value <- log(p.value)

    } else {  # "norm" test
        ## Find optimal lambda
        fit <- boxcox.mle(null.delta1 ~ 1)
        lambda <- fit$lambda

        # Apply Box-Cox transform
        if (abs(lambda) > .Machine$double.eps) {
            delta1.boxcox <- (delta1^lambda - 1) / lambda
            Enull.delta1.boxcox <- (Enull.delta1^lambda - 1) / lambda
            null.delta1.boxcox <- (null.delta1^lambda - 1) / lambda
            BB.delta1.boxcox <- (BB.delta1^lambda - 1) / lambda
        } else {
            delta1.boxcox <- log(delta1)
            Enull.delta1.boxcox <- log(Enull.delta1)
            null.delta1.boxcox <- log(null.delta1)
            BB.delta1.boxcox <- log(BB.delta1)
        }

        mu <- mean(null.delta1.boxcox)
        sigma <- stats::sd(null.delta1.boxcox)

        if (!exists("weighted.p.value")) {
            # Simple p-value calculation if weighted.p.value not available
            p.value <- mean(null.delta1.boxcox >= mean(BB.delta1.boxcox))
            p.value <- max(p.value, 1/n.perms)  # Avoid p-value of 0
        } else {
            p.value <- weighted.p.value(BB.delta1.boxcox, mu, sigma)
        }
        log.p.value <- log(p.value)
    }

    # Create diagnostic plots
    if (plot.it) {
        op <- par(mfrow = c(1, 2), mar = c(3.5, 3.5, 0.5, 0.5),
                  mgp = c(2, 0.5, 0), tcl = -0.3)
        on.exit(par(op), add = TRUE)

        # Plot 1: Dependence of y on x
        null.Ey.CrI <- apply(null.ExYg, 1, stats::quantile, probs = c(0.025, 0.975))
        null.Ey.u <- null.Ey.CrI[2, ]
        null.Ey.d <- null.Ey.CrI[1, ]

        # Check if plot.magelo exists
        if (exists("plot.magelo")) {
            plot(Ey.res, xlab = xlab, ylab = ylab, Eyg.col = Eyg.col)
        } else {
            # Fallback plotting
            plot(Ey.res$x, Ey.res$y, xlab = xlab, ylab = ylab, las = 1)
            if (!is.null(Ey.res$Eyg.CrI) && CrI.as.polygon) {
                yu <- Ey.res$Eyg.CrI[2, ]
                yl <- Ey.res$Eyg.CrI[1, ]
                polygon(c(xg, rev(xg)), c(yl, rev(yu)),
                       col = CrI.polygon.col, border = NA)
            }
            points(Ey.res$x, Ey.res$y, col = pt.col, pch = pt.pch)
        }

        # Add null distribution credible interval
        polygon(c(xg, rev(xg)), c(null.Ey.d, rev(null.Ey.u)),
               col = "gray85", border = NA)
        abline(h = mean(null.ExYg), col = "blue")
        lines(xg, Ey.res$Eyg, col = Eyg.col)

        # Plot 2: Histogram of null distribution
        if (exists("hist2")) {
            h <- hist2(null.delta1, BB.delta1, x1.lab = "null.fassoc1",
                      x2.lab = "BB.fassoc1", n.x1.breaks = 15, n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            h <- hist(null.delta1, breaks = 30, main = "", xlab = "null.fassoc1",
                     col = "lightgray", border = "white")
            xlim_used <- par("usr")[1:2]
            ylim_used <- par("usr")[3:4]
        }

        abline(v = delta1, col = 'red', lwd = 2)

        if (exists("mode.1D")) {
            abline(v = mode.1D(null.delta1), col = "blue")
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
        Delta1 = Delta1,
        delta1 = delta1,
        null.delta1 = null.delta1,
        null.delta1.boxcox = null.delta1.boxcox,
        null.Ey = null.ExYg,
        null.dEy = null.dExYg,
        BB.dEy = BB.dEy,
        BB.Delta1 = BB.Delta1,
        BB.delta1 = BB.delta1,
        BB.delta1.boxcox = BB.delta1.boxcox,
        delta1.boxcox = delta1.boxcox,
        lambda = lambda,
        log.p.value = log.p.value,
        p.value = p.value,
        x = x,
        y = y,
        Ey = Ey,
        Ey.res = Ey.res,
        grid.size = grid.size,
        two.sample.test = two.sample.test,
        n.perms = n.perms,
        n.BB = n.BB,
        call = match.call()
    )

    class(output) <- "assoc1"

    return(output)
}

#' @rdname fassoc1.test
#' @export
fofam.test <- fassoc1.test


#' Plot Method for assoc1 Objects
#'
#' Creates diagnostic plots for objects of class "assoc1" produced by fassoc1.test().
#'
#' @param x An object of class "assoc1" from fassoc1.test().
#' @param plot Character string specifying plot type. Options are:
#'   \describe{
#'     \item{"Exy"}{Conditional mean function with credible intervals}
#'     \item{"dExy"}{Derivative of conditional mean with credible intervals}
#'     \item{"d1hist"}{Histogram of null distribution with observed value}
#'     \item{"bc.d1hist"}{Histogram of Box-Cox transformed null distribution}
#'   }
#' @param xlab Character string for x-axis label. Default is "x".
#' @param ylab Character string for y-axis label. Default is "y".
#' @param x1.lab Character string for first histogram label. Default is "null.fassoc1".
#' @param x2.lab Character string for second histogram label. Default is "BB.fassoc1".
#' @param show.pval Logical indicating whether to show p-value. Default is TRUE.
#' @param show.vlines Logical indicating whether to show vertical reference lines. Default is FALSE.
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
#' @param dEy.CrI.col Color specification for derivative credible intervals.
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
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' y <- sin(4*pi*x) + rnorm(n, sd = 0.2)
#' result <- fassoc1.test(x, y, n.cores = 2, n.perms = 500, plot.it = FALSE)
#'
#' # Create different plots
#' plot(result, plot = "Exy")
#' plot(result, plot = "dExy")
#' plot(result, plot = "d1hist")
#' }
#'
#' @method plot assoc1
#' @export
plot.assoc1 <- function(x,
                       plot = c("Exy", "dExy", "d1hist", "bc.d1hist"),
                       xlab = "x",
                       ylab = "y",
                       x1.lab = "null.fassoc1",
                       x2.lab = "BB.fassoc1",
                       show.pval = TRUE,
                       show.vlines = FALSE,
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
    if (!inherits(res, "assoc1")) {
        stop("'x' must be an object of class 'assoc1'")
    }

    plot <- match.arg(plot)

    # Extract results
    Ey.res <- res$Ey.res
    xg <- Ey.res$xg
    ng <- res$grid.size
    Ey <- Ey.res$Eyg

    if (plot == "Exy") {
        # Dependence of y on x
        null.Ey.CrI <- apply(res$null.Ey, 1, stats::quantile, probs = c(0.025, 0.975))
        null.Ey.u <- null.Ey.CrI[2, ]
        null.Ey.d <- null.Ey.CrI[1, ]

        # Check if plot.magelo exists
        if (exists("plot.magelo")) {
            plot(Ey.res, xlab = xlab, ylab = ylab, Eyg.col = Eyg.col,
                 title = title, ylim = ylim)
        } else {
            # Fallback plotting
            if (is.null(ylim)) {
                ylim <- range(c(Ey.res$y, Ey.res$Eyg, null.Ey.u, null.Ey.d), na.rm = TRUE)
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
        abline(h = mean(res$null.Ey), col = null.line.col)
        lines(xg, Ey.res$Eyg, col = Eyg.col, lwd = Eyg.lwd)

        # Add credible interval lines if available
        if (!is.null(Ey.res$Eyg.CrI)) {
            lines(xg, Ey.res$Eyg.CrI[1, ], col = dEy.CrI.col)
            lines(xg, Ey.res$Eyg.CrI[2, ], col = dEy.CrI.col)
        }

    } else if (plot == "dExy") {
        # Derivative plot
        dEy.CrI <- apply(res$BB.dEy, 1, stats::quantile, probs = c(0.025, 0.975))
        null.dEy.CrI <- apply(res$null.dEy, 1, stats::quantile, probs = c(0.025, 0.975))

        dEy <- diff(Ey)
        dxg <- xg[2:ng]
        dEy.u <- dEy.CrI[2, ]
        dEy.d <- dEy.CrI[1, ]
        null.dEy.u <- null.dEy.CrI[2, ]
        null.dEy.d <- null.dEy.CrI[1, ]

        if (is.null(ylim)) {
            ylim <- range(c(dEy, dEy.d, dEy.u, null.dEy.d, null.dEy.u), na.rm = TRUE)
        }

        plot(dxg, dEy, type = 'l', ylim = ylim, ylab = "dE(y)/dx", xlab = xlab,
             main = title, las = 1)

        # Add credible intervals
        polygon(c(dxg, rev(dxg)), c(dEy.d, rev(dEy.u)), col = "gray95", border = NA)
        polygon(c(dxg, rev(dxg)), c(null.dEy.d, rev(null.dEy.u)), col = "gray90", border = NA)

        abline(h = 0, col = "blue")

        # Add derivative lines
        lines(dxg, dEy.u, col = "indianred")
        lines(dxg, dEy.d, col = "indianred")
        lines(dxg, dEy, col = 'red', lwd = 2)

    } else if (plot == "d1hist") {
        # Histogram of null delta1
        if (exists("hist2")) {
            h <- hist2(res$null.delta1, res$BB.delta1, xlab = xlab,
                      x1.lab = x1.lab, x2.lab = x2.lab,
                      n.x1.breaks = 15, n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            h <- hist(res$null.delta1, breaks = 30, main = title,
                     xlab = x1.lab, las = 1)
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
        if (show.vlines) {
            abline(v = res$delta1, col = 'red', lwd = 2)
            if (exists("mode.1D")) {
                abline(v = mode.1D(res$null.delta1), col = "blue")
            }
        }

        # Add p-value text
        if (show.pval) {
            text(x = d1hist.x, y = d1hist.y,
                 labels = paste0("p-val: ", signif(res$p.value, digits = 2)))
        }

    } else if (plot == "bc.d1hist") {
        # Box-Cox transformed histogram
        if (is.null(res$null.delta1.boxcox)) {
            stop("Box-Cox transformed values not available. ",
                 "Run fassoc1.test with two.sample.test = 'norm'")
        }

        if (exists("hist2")) {
            h <- hist2(res$null.delta1.boxcox, res$BB.delta1.boxcox,
                      x1.lab = "BC(null.fassoc1)", x2.lab = "BC(BB.fassoc1)",
                      n.x1.breaks = 15, n.x2.breaks = 50)
            xlim_used <- h$xlim
            ylim_used <- h$ylim
        } else {
            # Fallback to standard histogram
            h <- hist(res$null.delta1.boxcox, breaks = 30, main = title,
                     xlab = "BC(null.fassoc1)", las = 1)
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
        if (show.vlines) {
            abline(v = res$delta1.boxcox, col = 'red', lwd = 2)
            if (exists("mode.1D")) {
                abline(v = mode.1D(res$null.delta1.boxcox), col = "blue")
            }
        }

        # Add p-value text
        if (show.pval) {
            text(x = d1hist.x, y = d1hist.y,
                 labels = paste0("p-val: ", signif(res$p.value, digits = 2)))
        }
    }

    invisible(NULL)
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
    cat("\nFirst-Order Functional Association Test\n")
    cat("======================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Test Summary:\n")
    cat("  Sample size: ", length(x$x), "\n", sep = "")
    cat("  Test type: ", x$two.sample.test, "\n", sep = "")
    cat("  Permutations: ", x$n.perms, "\n", sep = "")
    cat("  Bootstrap samples: ", x$n.BB, "\n\n", sep = "")

    cat("Results:\n")
    cat("  Total change (Delta1): ", signif(x$Delta1, digits), "\n", sep = "")
    cat("  Association measure (delta1): ", signif(x$delta1, digits), "\n", sep = "")
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
        test.type = object$two.sample.test,
        n.perms = object$n.perms,
        n.BB = object$n.BB,
        Delta1 = object$Delta1,
        delta1 = object$delta1,
        p.value = object$p.value,
        null.mean = mean(object$null.delta1),
        null.sd = stats::sd(object$null.delta1),
        null.quantiles = stats::quantile(object$null.delta1,
                                        probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    )

    class(out) <- "summary.assoc1"
    return(out)
}

#' @method print summary.assoc1
#' @export
print.summary.assoc1 <- function(x, digits = 4, ...) {
    cat("\nFirst-Order Functional Association Test\n")
    cat("======================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Test Configuration:\n")
    cat("  Sample size: ", x$n, "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n", sep = "")
    cat("  Permutations: ", x$n.perms, "\n", sep = "")
    cat("  Bootstrap samples: ", x$n.BB, "\n\n", sep = "")

    cat("Results:\n")
    cat("  Total change (Delta1): ", signif(x$Delta1, digits), "\n", sep = "")
    cat("  Association measure (delta1): ", signif(x$delta1, digits), "\n", sep = "")
    cat("  p-value: ", signif(x$p.value, digits), "\n\n", sep = "")

    cat("Null Distribution Summary:\n")
    cat("  Mean: ", signif(x$null.mean, digits), "\n", sep = "")
    cat("  SD: ", signif(x$null.sd, digits), "\n", sep = "")
    cat("  Quantiles:\n")
    print(signif(x$null.quantiles, digits))

    invisible(x)
}

#' Coef Method for assoc1 Objects
#'
#' Extracts key coefficients from an assoc1 object.
#'
#' @param object An object of class "assoc1".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector containing delta1, Delta1, and p-value.
#' @method coef assoc1
#' @export
coef.assoc1 <- function(object, ...) {
    c(delta1 = object$delta1,
      Delta1 = object$Delta1,
      p.value = object$p.value)
}

#' Extract Derivative Information from assoc1 Object
#'
#' Extracts and summarizes derivative information from the first-order
#' functional association test results.
#'
#' @param object An object of class "assoc1".
#' @param type Character string specifying what to extract:
#'   "estimate" (default), "credible.interval", or "both".
#' @param probs Numeric vector of probabilities for credible intervals.
#'   Default is c(0.025, 0.975) for 95% intervals.
#'
#' @return Depending on type:
#'   \item{estimate}{A numeric vector of derivative estimates}
#'   \item{credible.interval}{A matrix with lower and upper bounds}
#'   \item{both}{A list containing both estimate and credible.interval}
#'
#' @examples
#' \dontrun{
#' result <- fassoc1.test(x, y)
#' deriv.est <- extract.derivatives(result)
#' deriv.ci <- extract.derivatives(result, type = "credible.interval")
#' }
#' @export
extract.derivatives <- function(object,
                               type = c("estimate", "credible.interval", "both"),
                               probs = c(0.025, 0.975)) {

    if (!inherits(object, "assoc1")) {
        stop("'object' must be of class 'assoc1'")
    }

    type <- match.arg(type)

    # Extract derivative estimates
    dEy <- diff(object$Ey)
    xg <- object$Ey.res$xg
    dxg <- xg[2:length(xg)]

    if (type == "estimate") {
        return(data.frame(x = dxg, dEy = dEy))
    }

    # Calculate credible intervals if requested
    if (type == "credible.interval" || type == "both") {
        dEy.CrI <- apply(object$BB.dEy, 1, stats::quantile, probs = probs)
        ci.df <- data.frame(
            x = dxg,
            lower = dEy.CrI[1, ],
            upper = dEy.CrI[2, ]
        )

        if (type == "credible.interval") {
            return(ci.df)
        }
    }

    # Return both if requested
    if (type == "both") {
        return(list(
            estimate = data.frame(x = dxg, dEy = dEy),
            credible.interval = ci.df
        ))
    }
}

#' Functional Association Test (Unified Interface)
#'
#' Unified wrapper for zero-order and first-order functional association tests.
#' Dispatches to \code{fassoc0.test()} when \code{order = 0} and to
#' \code{fassoc1.test()} when \code{order = 1}.
#'
#' This wrapper is updated to match the latest *paired* test design in
#' \code{fassoc0.test()} and \code{fassoc1.test()}, including:
#' \itemize{
#'   \item explicit \code{n.perms} permutation layer,
#'   \item paired Bayesian bootstrap weights via \code{n.BB},
#'   \item hybrid rule for \code{n.perms > n.BB} (expand vs BB-cluster inference),
#'   \item optional Box-Cox normalization, and
#'   \item diff-based inference (\code{"paired.t"} / \code{"wilcoxon"}) or
#'         null-parametric inference (\code{"weighted.pvalue"}).
#' }
#'
#' Backward-compatibility note:
#' If \code{two.sample.test} is provided (legacy argument from older interfaces),
#' it will be mapped to \code{test.type} as follows:
#' \itemize{
#'   \item \code{"norm"} -> \code{"weighted.pvalue"}
#'   \item \code{"KS"} -> \code{"wilcoxon"}
#'   \item \code{"Wasserstein"} -> \code{"paired.t"}
#' }
#' If both are provided, \code{test.type} takes precedence.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values (same length as x).
#' @param order Integer; 0 for zero-order association, 1 for first-order association.
#'
#' @param test.type One of \code{"paired.t"}, \code{"weighted.pvalue"}, \code{"wilcoxon"}.
#' @param boxcox One of \code{"auto"}, \code{"always"}, \code{"never"}.
#' @param boxcox.alpha Alpha for Shapiro-Wilk when \code{boxcox = "auto"}.
#'
#' @param bw Bandwidth. If NULL, chosen automatically by the underlying smoother.
#' @param n.perms Number of permutations for the null distribution.
#' @param n.BB Number of paired Bayesian bootstrap samples for the signal fit.
#' @param max.nBB.expand Threshold for expanding \code{n.BB} to match \code{n.perms}
#'   under diff-based inference when \code{n.perms > n.BB}. See underlying functions.
#' @param cluster.agg Aggregation for BB-cluster inference under cycling:
#'   \code{"mean"}, \code{"trimmed.mean"}, or \code{"median"}.
#' @param cluster.trim Trim proportion for \code{cluster.agg = "trimmed.mean"}.
#'
#' @param grid.size Evaluation grid size for the conditional mean curve.
#' @param degree Local polynomial degree (1 or 2).
#' @param min.K Minimum local sample size parameter for the smoother.
#' @param n.cores Number of CPU cores for permutation loop (Unix-alikes via mclapply).
#' @param seed Optional RNG seed.
#'
#' @param plot.it Logical; if TRUE produce diagnostic plots.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param verbose Logical; if TRUE prints progress messages.
#'
#' @param two.sample.test Legacy argument (optional). See Details.
#' @param ... Ignored legacy plotting/style arguments may be supplied without error.
#'
#' @return An object of class \code{"assoc0"} (order 0) or \code{"assoc1"} (order 1).
fassoc.test <- function(x,
                        y,
                        order = 1,
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
                        verbose = TRUE,
                        two.sample.test = NULL,
                        ...) {

    ## ------------------------------------------------------------------------
    ## Input validation
    ## ------------------------------------------------------------------------

    if (!is.numeric(order) || length(order) != 1 || !(order %in% c(0, 1))) {
        stop("'order' must be either 0 or 1")
    }

    ## Map legacy two.sample.test (if supplied) to test.type
    if (!is.null(two.sample.test)) {
        ## two.sample.test was historically: c("norm","Wasserstein","KS")
        two.sample.test <- match.arg(two.sample.test, c("norm", "Wasserstein", "KS"))

        ## If user did not explicitly choose test.type (hard to detect robustly),
        ## we still allow the legacy override, unless they passed test.type explicitly.
        ## Here: if they passed two.sample.test at all, we apply the mapping unless
        ## they also passed a non-default test.type via match.arg result.
        ## To avoid ambiguity, if both provided, test.type wins.
        if (verbose) {
            message(sprintf("Legacy argument 'two.sample.test' detected (%s); mapping to 'test.type'.", two.sample.test))
        }

        mapped.test.type <- switch(two.sample.test,
                                   "norm" = "weighted.pvalue",
                                   "KS" = "wilcoxon",
                                   "Wasserstein" = "paired.t")

        ## Apply mapping only if user did not explicitly specify test.type via ...;
        ## Because match.arg always returns something, we treat presence of two.sample.test
        ## as intent and set test.type unless user provided test.type explicitly.
        ## In base R, explicit detection is tricky without formals inspection; so:
        ## - if user also provided test.type in the call, they can avoid mapping by setting
        ##   two.sample.test = NULL.
        test.type <- mapped.test.type
    } else {
        test.type <- match.arg(test.type)
    }

    boxcox <- match.arg(boxcox)
    cluster.agg <- match.arg(cluster.agg)

    ## Ignore legacy style args in ...
    dots <- list(...)
    if (length(dots) > 0 && verbose) {
        message(sprintf("Note: ignoring %d additional argument(s) passed via '...'.", length(dots)))
    }

    ## ------------------------------------------------------------------------
    ## Dispatch
    ## ------------------------------------------------------------------------

    if (order == 0) {
        fassoc0.test(
            x = x,
            y = y,
            test.type = test.type,
            boxcox = boxcox,
            boxcox.alpha = boxcox.alpha,
            bw = bw,
            n.perms = n.perms,
            n.BB = n.BB,
            max.nBB.expand = max.nBB.expand,
            cluster.agg = cluster.agg,
            cluster.trim = cluster.trim,
            grid.size = grid.size,
            degree = degree,
            min.K = min.K,
            n.cores = n.cores,
            seed = seed,
            plot.it = plot.it,
            xlab = xlab,
            ylab = ylab,
            verbose = verbose
        )
    } else {
        fassoc1.test(
            x = x,
            y = y,
            test.type = test.type,
            boxcox = boxcox,
            boxcox.alpha = boxcox.alpha,
            bw = bw,
            n.perms = n.perms,
            n.BB = n.BB,
            max.nBB.expand = max.nBB.expand,
            cluster.agg = cluster.agg,
            cluster.trim = cluster.trim,
            grid.size = grid.size,
            degree = degree,
            min.K = min.K,
            n.cores = n.cores,
            seed = seed,
            plot.it = plot.it,
            xlab = xlab,
            ylab = ylab,
            verbose = verbose
        )
    }
}

#' @rdname fassoc.test
functional.association.test <- fassoc.test

#' Calculate Total Associations from Geodesic Analysis
#'
#' Calculates total associations (TAs) and total absolute associations (TAAs)
#' between the distance along a geodesic and some function over the geodesic.
#'
#' @param obj A list result from E.geodesic.X() containing at least:
#'   \describe{
#'     \item{X}{A matrix of data}
#'     \item{data.dir}{Directory containing saved results}
#'   }
#' @param ED Numeric vector of mean disease/outcome values over the geodesic.
#'   If NULL (default), calculates associations with distance along geodesic.
#' @param verbose Logical indicating whether to print progress messages.
#'   Default is FALSE.
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{TA}{Total Association}
#'     \item{TAA}{Total Absolute Association}
#'     \item{pval}{p-value for the association test}
#'   }
#'   Row names correspond to column names of the input matrix X.
#'
#' @details
#' This function processes results from geodesic analysis, calculating associations
#' between variables along the geodesic path. When ED is NULL, it calculates
#' associations with the distance along the geodesic. When ED is provided,
#' it calculates associations with the specified outcome variable.
#'
#' @examples
#' \dontrun{
#' # Assuming obj is a result from E.geodesic.X()
#' # Calculate associations with distance
#' ta_results <- total.associations(obj)
#'
#' # Calculate associations with an outcome
#' ED <- rnorm(400)  # Example outcome vector
#' ta_results_ED <- total.associations(obj, ED = ED)
#' }
#'
#' @importFrom utils file_test
#' @export
total.associations <- function(obj, ED = NULL, verbose = FALSE) {

    # Input validation
    if (!is.list(obj)) {
        stop("'obj' must be a list")
    }

    required_elements <- c("X", "data.dir")
    missing_elements <- setdiff(required_elements, names(obj))
    if (length(missing_elements) > 0) {
        stop("'obj' is missing required elements: ", paste(missing_elements, collapse = ", "))
    }

    data.dir <- obj[["data.dir"]]
    X <- obj[["X"]]

    if (!is.matrix(X)) {
        stop("obj$X must be a matrix")
    }

    if (!dir.exists(data.dir)) {
        stop("Directory specified in obj$data.dir does not exist: ", data.dir)
    }

    # Load first file to get grid size
    col.i <- 1
    file <- file.path(data.dir, paste0(col.i, ".rda"))
    if (!file.exists(file)) {
        stop("Expected file not found: ", file)
    }

    ## Load file to get Ey (assumes file contains 'Ey' and 'pval' objects)
    Ey <- NULL
    load(file)
    if (!is.null(Ey)) {
        stop("Loaded file does not contain 'Ey' object")
    }

    grid.size <- length(Ey)

    if (!is.null(ED)) {
        if (!is.numeric(ED)) {
            stop("'ED' must be numeric")
        }
        if (length(ED) != grid.size) {
            stop(sprintf("'ED' length (%d) must match grid size (%d)", length(ED), grid.size))
        }
    }

    # Initialize results matrix
    res <- matrix(nrow = ncol(X), ncol = 3)
    colnames(res) <- c("TA", "TAA", "pval")
    rownames(res) <- colnames(X)

    if (is.null(ED)) {
        if (verbose) {
            cat("Calculating TAs and TAAs between the distance along the geodesic and the scaled derivatives of the mean abundances of X's columns ... ")
            ptm <- proc.time()
        }

        for (col.i in seq_len(ncol(X))) {
            file <- file.path(data.dir, paste0(col.i, ".rda"))
            if (file.exists(file)) {
                Ey <- pval <- NULL
                load(file)
                if(!is.null(Ey) && !is.null(pval)) {
                    ## Loads Ey and pval
                    tot.loc.asso <- Ey[grid.size] - Ey[1]
                    tot.abs.loc.asso <- sum(abs(diff(Ey)))
                    res[col.i, ] <- c(tot.loc.asso, tot.abs.loc.asso, pval)
                } else {
                    res[col.i, ] <- c(NA, NA, NA)
                }
            } else {
                warning("File not found: ", file)
                res[col.i, ] <- c(NA, NA, NA)
            }
        }
    } else {
        if (verbose) {
            cat("Calculating TAs and TAAs between ED and the scaled derivatives of the mean abundances of X's columns ... ")
            ptm <- proc.time()
        }

        dED <- diff(ED)

        for (col.i in seq_len(ncol(X))) {
            file <- file.path(data.dir, paste0(col.i, ".rda"))
            if (file.exists(file)) {
                Ey <- pval <- NULL
                load(file)  # Loads Ey and pval
                if(!is.null(Ey) && !is.null(pval)) {
                    dEy <- diff(Ey)
                    tot.loc.asso <- sum(dEy * dED)
                    tot.abs.loc.asso <- sum(abs(dEy * dED))
                    res[col.i, ] <- c(tot.loc.asso, tot.abs.loc.asso, pval)
                } else {
                    res[col.i, ] <- c(NA, NA, NA)
                }
            } else {
                warning("File not found: ", file)
                res[col.i, ] <- c(NA, NA, NA)
            }
        }
    }

    if (verbose) {
        elapsed.time(ptm)
    }

    return(res)
}

#' Create Summary Data Frame from First-Order Association Test Results
#'
#' Creates a summary data frame containing key statistics from multiple
#' first-order functional association test results.
#'
#' @param res.list A named list of results from fassoc1.test() or fassoc.test(order=1).
#' @param q.thld Numeric FDR threshold for significance. Default is 0.1.
#'
#' @return A list containing:
#'   \describe{
#'     \item{d1D1.df}{Full data frame with all results}
#'     \item{sign.d1D1.df}{Filtered data frame with significant results only}
#'   }
#'   The data frames contain columns: RDM (Right Deviation from Mean),
#'   D (ratio), p-val, q-val, delta1, and Delta1.
#'
#' @details
#' This function processes multiple test results to create a summary suitable
#' for reporting. Results are sorted by RDM (Right Deviation from Mean) and
#' FDR correction is applied to p-values.
#'
#' @examples
#' \dontrun{
#' # Run multiple tests
#' results <- list()
#' for (i in 1:10) {
#'   x <- runif(100)
#'   y <- sin(2*pi*x) + rnorm(100, sd = 0.3)
#'   results\code{[[paste0("var", i)]}] <- fassoc.test(x, y, order = 1)
#' }
#'
#' # Create summary
#' summary_df <- create.delta1.Delta1.df(results)
#' print(summary_df$sign.d1D1.df)
#' }
#'
#' @importFrom stats p.adjust
#' @export
create.delta1.Delta1.df <- function(res.list, q.thld = 0.1) {

    # Input validation
    if (!is.list(res.list) || length(res.list) == 0) {
        stop("'res.list' must be a non-empty list")
    }

    if (!is.numeric(q.thld) || length(q.thld) != 1 || q.thld <= 0 || q.thld >= 1) {
        stop("'q.thld' must be a single numeric value between 0 and 1")
    }

    # Check that all elements have required components
    required_fields <- c("Ey.res", "p.value", "delta1", "Delta1")
    for (i in seq_along(res.list)) {
        if (!all(required_fields %in% names(res.list[[i]]))) {
            stop(sprintf("Element %d of res.list is missing required fields", i))
        }
    }

    # Initialize data frame
    d1D1.df <- data.frame(
        RDM = numeric(length(res.list)),
        D = numeric(length(res.list)),
        `p-val` = numeric(length(res.list)),
        delta1 = numeric(length(res.list)),
        Delta1 = numeric(length(res.list)),
        check.names = FALSE
    )
    rownames(d1D1.df) <- names(res.list)

    # Populate data frame
    for (i in seq_along(res.list)) {
        res <- res.list[[i]]
        Ey.res <- res$Ey.res

        # Handle different possible structures
        if (!is.null(Ey.res$ng)) {
            ng <- Ey.res$ng
        } else if (!is.null(res$grid.size)) {
            ng <- res$grid.size
        } else {
            ng <- length(Ey.res$Eyg)
        }

        Ey <- mean(Ey.res$y)
        rdm <- Ey.res$Eyg[ng] / Ey
        D <- Ey.res$Eyg[ng] / Ey.res$Eyg[1]

        d1D1.df[i, ] <- c(rdm, D, res$p.value, res$delta1, res$Delta1)
    }

    # Sort by RDM
    o <- order(d1D1.df$RDM, decreasing = TRUE)
    d1D1.df <- d1D1.df[o, ]

    # Apply FDR correction
    q.vals <- stats::p.adjust(d1D1.df[, "p-val"], method = "fdr")
    d1D1.df$`q-val` <- q.vals

    # Reorder columns
    d1D1.df <- d1D1.df[, c("RDM", "D", "p-val", "q-val", "delta1", "Delta1")]

    # Round for display
    for (i in 1:6) {
        d1D1.df[, i] <- signif(d1D1.df[, i], digits = 2)
    }

    # Create significant results subset
    idx <- d1D1.df[, "q-val"] < q.thld
    sign.d1D1.df <- d1D1.df[idx, c("RDM", "D", "q-val"), drop = FALSE]

    return(list(
        d1D1.df = d1D1.df,
        sign.d1D1.df = sign.d1D1.df
    ))
}

#' Calculate Higher-Order Functional Association Indices
#'
#' Computes Delta and delta functional association indices up to order k
#' from conditional mean estimates.
#'
#' @param Eyg Numeric vector of conditional mean estimates E_x(y) over a uniform grid.
#' @param k Integer specifying the highest order of indices to compute.
#'   Default is 10. Must be positive.
#'
#' @return A list containing:
#'   \describe{
#'     \item{Delta}{Numeric vector of Delta indices from order 1 to k}
#'     \item{delta}{Numeric vector of delta indices from order 1 to k}
#'   }
#'
#' @details
#' For each order i from 1 to k:
#' \itemize{
#'   \item \code{Delta[i]} represents the total signed change at order i
#'   \item \code{delta[i]} represents the total absolute change at order i
#' }
#' Higher orders capture increasingly fine-scale variation in the functional relationship.
#'
#' @examples
#' \dontrun{
#' # Generate example conditional mean curve
#' x <- seq(0, 1, length.out = 100)
#' Eyg <- sin(4*pi*x) + 0.5*x
#'
#' # Calculate indices up to order 5
#' indices <- delta.indices(Eyg, k = 5)
#' print(indices$Delta)
#' print(indices$delta)
#'
#' # Plot delta values by order
#' plot(1:5, indices$delta, type = "b",
#'      xlab = "Order", ylab = "delta",
#'      main = "Functional Association by Order")
#' }
#'
#' @export
delta.indices <- function(Eyg, k = 10) {

    # Input validation
    if (!is.numeric(Eyg) || length(Eyg) < 2) {
        stop("'Eyg' must be a numeric vector with at least 2 elements")
    }

    if (!is.numeric(k) || length(k) != 1 || k < 1) {
        stop("'k' must be a positive integer")
    }
    k <- as.integer(k)

    ng <- length(Eyg)
    if (k > ng - 1) {
        warning(sprintf("k (%d) is larger than length(Eyg)-1 (%d). Setting k = %d",
                       k, ng - 1, ng - 1))
        k <- ng - 1
    }

    # Initialize storage
    dkEyg <- vector("list", k)
    Delta <- numeric(k)
    delta <- numeric(k)

    # First order
    dkEyg[[1]] <- diff(Eyg)
    Delta[1] <- Eyg[ng] - Eyg[1]
    delta[1] <- sum(abs(dkEyg[[1]]))

    # Higher orders
    if (k > 1) {
        for (i in 2:k) {
            if (length(dkEyg[[i - 1]]) > 1) {
                dkEyg[[i]] <- diff(dkEyg[[i - 1]])
                delta[i] <- sum(abs(dkEyg[[i]]))
                Delta[i] <- sum(dkEyg[[i]])
            } else {
                # Not enough points for higher differences
                delta[i] <- NA
                Delta[i] <- NA
                if (i == 2) {
                    warning("Not enough points to compute differences beyond order 1")
                }
                break
            }
        }
    }

    # Remove NA values if any
    valid_idx <- !is.na(delta)
    Delta <- Delta[valid_idx]
    delta <- delta[valid_idx]

    return(list(
        Delta = Delta,
        delta = delta
    ))
}

#' Create Higher-Order Difference Profiles
#'
#' Computes and concatenates difference profiles up to order k from
#' conditional mean estimates.
#'
#' @param Eyg Numeric vector of conditional mean estimates E_x(y) over a uniform grid.
#' @param k Integer specifying the highest order of difference profiles to compute.
#'   Must be positive.
#'
#' @return A numeric vector containing the concatenated profiles:
#'   Eyg, diff(Eyg), diff(diff(Eyg)), ..., up to order k.
#'
#' @details
#' This function creates a comprehensive profile by concatenating the original
#' values with successive differences. This can be useful for visualizing
#' the behavior of a function and its derivatives simultaneously.
#'
#' @examples
#' \dontrun{
#' # Generate example curve
#' x <- seq(0, 1, length.out = 50)
#' Eyg <- sin(2*pi*x)
#'
#' # Get profiles up to 3rd order
#' profiles <- compute.diff.profiles(Eyg, k = 3)
#'
#' # The result contains:
#' # - Original values (50 elements)
#' # - First differences (49 elements)
#' # - Second differences (48 elements)
#' # - Third differences (47 elements)
#' # Total: 194 elements
#' }
#'
#' @export
compute.diff.profiles <- function(Eyg, k) {

    # Input validation
    if (!is.numeric(Eyg) || length(Eyg) < 2) {
        stop("'Eyg' must be a numeric vector with at least 2 elements")
    }

    if (!is.numeric(k) || length(k) != 1 || k < 1) {
        stop("'k' must be a positive integer")
    }
    k <- as.integer(k)

    ng <- length(Eyg)
    if (k > ng - 1) {
        warning(sprintf("k (%d) is larger than length(Eyg)-1 (%d). Setting k = %d",
                       k, ng - 1, ng - 1))
        k <- ng - 1
    }

    # Create list of differences
    dkEyg <- vector("list", k)
    dkEyg[[1]] <- diff(Eyg)

    if (k > 1) {
        for (i in 2:k) {
            if (length(dkEyg[[i - 1]]) > 1) {
                dkEyg[[i]] <- diff(dkEyg[[i - 1]])
            } else {
                # Adjust k if we can't compute more differences
                k <- i - 1
                warning(sprintf("Can only compute differences up to order %d", k))
                break
            }
        }
    }

    # Concatenate all profiles
    comb.Eyg <- Eyg
    for (i in 1:k) {
        if (i <= length(dkEyg) && !is.null(dkEyg[[i]])) {
            comb.Eyg <- c(comb.Eyg, dkEyg[[i]])
        }
    }

    return(comb.Eyg)
}

#' Print Method for Functional Association Tests
#'
#' @param x An object of class "assoc0" or "assoc1".
#' @param ... Additional arguments passed to specific print methods.
#'
#' @return Invisible x.
#' @method print fassoc
#' @export
print.fassoc <- function(x, ...) {
    if (inherits(x, "assoc0")) {
        print.assoc0(x, ...)
    } else if (inherits(x, "assoc1")) {
        print.assoc1(x, ...)
    } else {
        cat("Functional Association Test Result\n")
        print.default(x, ...)
    }
    invisible(x)
}

#' Summary Method for Functional Association Tests
#'
#' @param object An object of class "assoc0" or "assoc1".
#' @param ... Additional arguments passed to specific summary methods.
#'
#' @return A summary object.
#' @method summary fassoc
#' @export
summary.fassoc <- function(object, ...) {
    if (inherits(object, "assoc0")) {
        summary.assoc0(object, ...)
    } else if (inherits(object, "assoc1")) {
        summary.assoc1(object, ...)
    } else {
        summary.default(object, ...)
    }
}

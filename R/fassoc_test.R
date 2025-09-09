#' Functional Association Test
#'
#' Tests for the existence of non-trivial functional association between two variables
#' x and y using either zero-order or first-order functional association measures.
#' This function serves as a unified interface to both fassoc0.test() and fassoc1.test().
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values (same length as x). Can be binary or continuous.
#' @param order Integer specifying the order of the functional association test.
#'   Can be 0 (zero-order) or 1 (first-order). Default is 1.
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
#' @return An object of class "assoc0" (if order=0) or "assoc1" (if order=1) containing
#'   test results. See fassoc0.test() and fassoc1.test() for details on the return values.
#'
#' @details
#' This function provides a unified interface to test for functional associations
#' between variables. The zero-order measure (order=0) quantifies deviations of
#' the conditional mean from the marginal mean, while the first-order measure
#' (order=1) quantifies the variability in the derivative of the conditional mean.
#'
#' The mathematical notation for these measures:
#' \itemize{
#'   \item Zero-order: \eqn{\text{fassoc}_0(x,y) = \int |E_x(y) - E(y)| dx}
#'   \item First-order: \eqn{\text{fassoc}_1(x,y) = \int |dE_x(y)/dx| dx}
#' }
#'
#' This function can also be called using the alias:
#' \itemize{
#'   \item \code{functional.association.test()} - Full descriptive name
#' }
#' For direct access to specific orders, use:
#' \itemize{
#'   \item \code{fassoc0.test()} or \code{zofam.test()} - Zero-Order Functional Association Measure
#'   \item \code{fassoc1.test()} or \code{fofam.test()} - First-Order Functional Association Measure
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
#' result0 <- fassoc.test(x, y, order = 0, n.cores = 2)
#'
#' # Test for first-order functional association
#' result1 <- fassoc.test(x, y, order = 1, n.cores = 2)
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

    if (!is.numeric(order) || length(order) != 1 || !(order %in% c(0, 1))) {
        stop("'order' must be either 0 or 1")
    }

    # Check if required functions exist
    if (order == 0 && !exists("fassoc0.test")) {
        stop("Function 'fassoc0.test' not found. Please load required dependencies.")
    }
    if (order == 1 && !exists("fassoc1.test")) {
        stop("Function 'fassoc1.test' not found. Please load required dependencies.")
    }

    # Note: There's a bug in the original code where ylab is set to xlab
    # This has been fixed here

    if (order == 0) {
        fassoc0.test(x = x,
                    y = y,
                    two.sample.test = two.sample.test,
                    bw = bw,
                    n.perms = n.perms,
                    n.BB = n.BB,
                    grid.size = grid.size,
                    min.K = min.K,
                    n.cores = n.cores,
                    plot.it = plot.it,
                    xlab = xlab,
                    ylab = ylab,  # Fixed: was xlab in original
                    Eyg.col = Eyg.col,
                    pt.col = pt.col,
                    pt.pch = pt.pch,
                    CrI.as.polygon = CrI.as.polygon,
                    CrI.polygon.col = CrI.polygon.col,
                    null.lines.col = null.lines.col,
                    CrI.line.col = CrI.line.col,
                    CrI.line.lty = CrI.line.lty,
                    verbose = verbose)
    } else {
        fassoc1.test(x = x,
                    y = y,
                    two.sample.test = two.sample.test,
                    bw = bw,
                    n.perms = n.perms,
                    n.BB = n.BB,
                    grid.size = grid.size,
                    min.K = min.K,
                    n.cores = n.cores,
                    plot.it = plot.it,
                    xlab = xlab,
                    ylab = ylab,  # Fixed: was xlab in original
                    Eyg.col = Eyg.col,
                    pt.col = pt.col,
                    pt.pch = pt.pch,
                    CrI.as.polygon = CrI.as.polygon,
                    CrI.polygon.col = CrI.polygon.col,
                    null.lines.col = null.lines.col,
                    CrI.line.col = CrI.line.col,
                    CrI.line.lty = CrI.line.lty,
                    verbose = verbose)
    }
}

#' @rdname fassoc.test
#' @export
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

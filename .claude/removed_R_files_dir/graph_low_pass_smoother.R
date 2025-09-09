#' Graph Low-Pass Smoother using Laplacian Eigenvectors
#'
#' @description
#' Implements a graph-based low-pass filter by projecting vertex-valued signals onto
#' the subspace spanned by the first k Laplacian eigenvectors. The optimal number of
#' eigenvectors is selected using one of three criteria: eigengap heuristic,
#' Generalized Cross-Validation (GCV), or spectral energy threshold.
#'
#' This method is particularly useful for denoising graph signals while preserving
#' the underlying graph structure. The smoother uses powers of (I - L) where L is
#' the normalized graph Laplacian, allowing for different smoothing strengths.
#'
#' @param adj.list List of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list[[i]]} contains the 1-based indices of vertices
#'   adjacent to vertex i.
#' @param weight.list List of numeric vectors containing edge weights corresponding
#'   to the adjacencies. Must have the same structure as \code{adj.list}, where
#'   \code{weight.list[[i]][j]} is the weight of edge from vertex i to vertex
#'   \code{adj.list[[i]][j]}.
#' @param y Numeric vector of length equal to the number of vertices, containing
#'   the signal value at each vertex.
#' @param n.evectors.to.compute Integer; number of Laplacian eigenvectors to compute.
#'   Defaults to \code{min(10, n.vertices)}.
#' @param min.num.eigenvectors Integer; minimum number of eigenvectors to consider
#'   in the candidate set. Must be at least 1.
#' @param max.num.eigenvectors Integer; maximum number of eigenvectors to consider.
#'   Must be between \code{min.num.eigenvectors} and \code{n.evectors.to.compute}.
#' @param tau.factor Positive real number in (0,1] that scales the kernel bandwidth
#'   as a fraction of the graph diameter. Smaller values result in more localized
#'   smoothing. Default is \code{1/n.evectors.to.compute}.
#' @param radius.factor Positive real number >= 1 that scales the effective radius
#'   of the smoothing kernel. Larger values extend the smoothing influence.
#'   Default is 10.0.
#' @param laplacian.power Positive odd integer between 1 and 10 specifying the power
#'   to which (I - L) is raised. Higher powers apply stronger smoothing by repeatedly
#'   applying the low-pass filter. Must be odd to preserve sign. Default is 3.
#' @param n.candidates Integer >= 1; number of candidate k values to evaluate in
#'   the grid search. Default is 20.
#' @param log.grid Logical; if TRUE, use logarithmically-spaced grid of k values,
#'   otherwise use linearly-spaced grid. Default is FALSE.
#' @param energy.threshold Numeric in (0,1]; fraction of total spectral energy to
#'   retain for the spectral energy selection criterion. Default is 0.9.
#' @param with.k.predictions Logical; if TRUE, return the full matrix of predictions
#'   for each candidate k value. Useful for diagnostics but increases memory usage.
#'   Default is TRUE.
#' @param verbose Logical; if TRUE, print progress messages during computation.
#'   Default is TRUE.
#'
#' @return An object of class \code{"graph_low_pass_smoother"}, which is a list
#' containing:
#' \describe{
#'   \item{evalues}{Numeric vector of computed Laplacian eigenvalues in ascending order.}
#'   \item{evectors}{Matrix with Laplacian eigenvectors as columns, corresponding to evalues.}
#'   \item{candidate.ks}{Integer vector of candidate k values that were evaluated.}
#'   \item{eigengaps}{Numeric vector of eigengaps (evalues[i+1] - evalues[i]).}
#'   \item{gcv.scores}{Numeric vector of GCV scores for each candidate k.}
#'   \item{spectral.energy}{Numeric vector of cumulative spectral energy for each k.}
#'   \item{opt.k.eigengap}{Integer index in candidate.ks corresponding to largest eigengap.}
#'   \item{opt.k.gcv}{Integer index in candidate.ks corresponding to minimum GCV score.}
#'   \item{opt.k.spectral.energy}{Integer index of first k meeting energy threshold.}
#'   \item{used.method}{Character string indicating which selection method was used
#'     ("eigengap", "GCV", or "spectral.energy").}
#'   \item{predictions}{Numeric vector of smoothed signal values at optimal k.}
#'   \item{k.predictions}{If with.k.predictions=TRUE, a matrix where column j contains
#'     predictions using candidate.ks[j] eigenvectors.}
#' }
#'
#' @details
#' The function implements three methods for selecting the optimal number of eigenvectors:
#' \itemize{
#'   \item \strong{Eigengap}: Selects k where the gap between consecutive eigenvalues is largest.
#'   \item \strong{GCV}: Minimizes the Generalized Cross-Validation score.
#'   \item \strong{Spectral Energy}: Selects smallest k retaining the specified energy fraction.
#' }
#'
#' The actual method used is determined automatically based on the data characteristics
#' and is reported in the \code{used.method} field of the output.
#'
#' @seealso \code{\link{print.graph_low_pass_smoother}}, \code{\link{summary.graph_low_pass_smoother}}
#'
#' @examples
#' \dontrun{
#' # Create a simple ring graph with noisy signal
#' n <- 50
#' adj.list <- lapply(1:n, function(i) c(ifelse(i==1, n, i-1), ifelse(i==n, 1, i+1)))
#' weight.list <- lapply(1:n, function(i) c(1, 1))
#'
#' # Generate smooth signal with noise
#' true.signal <- sin(2 * pi * (1:n) / n)
#' y <- true.signal + rnorm(n, sd = 0.3)
#'
#' # Apply low-pass smoother
#' result <- graph.low.pass.smoother(adj.list, weight.list, y, verbose = FALSE)
#'
#' # Plot results
#' plot(y, pch = 20, col = "gray", main = "Graph Low-Pass Smoothing")
#' lines(result$predictions, col = "red", lwd = 2)
#' lines(true.signal, col = "blue", lty = 2)
#' legend("topright", c("Noisy", "Smoothed", "True"),
#'        col = c("gray", "red", "blue"), lty = c(NA, 1, 2), pch = c(20, NA, NA))
#' }
#'
#' @export
graph.low.pass.smoother <- function(
    adj.list,
    weight.list,
    y,
    n.evectors.to.compute = min(10, length(adj.list)),
    min.num.eigenvectors = 1L,
    max.num.eigenvectors = n.evectors.to.compute,
    tau.factor = 1 / n.evectors.to.compute,
    radius.factor = 10.0,
    laplacian.power = 3,
    n.candidates = 20L,
    log.grid = FALSE,
    energy.threshold = 0.9,
    with.k.predictions = TRUE,
    verbose = TRUE
) {
    ## --- Input validation ---
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("'adj.list' and 'weight.list' must both be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }

    n.vertices <- length(adj.list)

    if (n.vertices == 0) {
        stop("Graph must have at least one vertex")
    }

    # Check structure consistency
    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Vertex %d: lengths of adj.list[[i]] and weight.list[[i]] differ", i))
        }

        # Check for valid indices
        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n.vertices)) {
                stop(sprintf("Vertex %d: adjacency list contains invalid vertex indices", i))
            }

            # Check for positive weights
            if (any(weight.list[[i]] <= 0)) {
                stop(sprintf("Vertex %d: all edge weights must be positive", i))
            }
        }
    }

    # Validate signal vector
    if (!is.numeric(y) || length(y) != n.vertices) {
        stop("'y' must be a numeric vector of length equal to number of vertices")
    }

    if (any(!is.finite(y))) {
        stop("'y' must not contain NA, NaN, or infinite values")
    }

    # Validate numeric parameters
    if (!is.numeric(n.evectors.to.compute) || length(n.evectors.to.compute) != 1 ||
        n.evectors.to.compute < 1 || n.evectors.to.compute > n.vertices) {
        stop("'n.evectors.to.compute' must be a single integer in [1, n.vertices]")
    }
    n.evectors.to.compute <- as.integer(n.evectors.to.compute)

    if (!is.numeric(min.num.eigenvectors) || length(min.num.eigenvectors) != 1 ||
        min.num.eigenvectors < 1) {
        stop("'min.num.eigenvectors' must be a single integer >= 1")
    }
    min.num.eigenvectors <- as.integer(min.num.eigenvectors)

    if (!is.numeric(max.num.eigenvectors) || length(max.num.eigenvectors) != 1 ||
        max.num.eigenvectors < min.num.eigenvectors ||
        max.num.eigenvectors > n.evectors.to.compute) {
        stop("'max.num.eigenvectors' must satisfy min <= max <= n.evectors.to.compute")
    }
    max.num.eigenvectors <- as.integer(max.num.eigenvectors)

    if (!is.numeric(tau.factor) || length(tau.factor) != 1 ||
        tau.factor <= 0 || tau.factor > 1) {
        stop("'tau.factor' must be a single number in (0,1]")
    }

    if (!is.numeric(radius.factor) || length(radius.factor) != 1 ||
        radius.factor < 1) {
        stop("'radius.factor' must be a single number >= 1")
    }

    if (!is.numeric(laplacian.power) || length(laplacian.power) != 1 ||
        laplacian.power %% 2 != 1 || laplacian.power < 1 || laplacian.power > 10) {
        stop("'laplacian.power' must be an odd integer in [1, 10]")
    }
    laplacian.power <- as.integer(laplacian.power)

    if (!is.numeric(n.candidates) || length(n.candidates) != 1 || n.candidates < 1) {
        stop("'n.candidates' must be a single integer >= 1")
    }
    n.candidates <- as.integer(n.candidates)

    if (!is.logical(log.grid) || length(log.grid) != 1) {
        stop("'log.grid' must be a single logical value")
    }

    if (!is.numeric(energy.threshold) || length(energy.threshold) != 1 ||
        energy.threshold <= 0 || energy.threshold > 1) {
        stop("'energy.threshold' must be a single number in (0,1]")
    }

    if (!is.logical(with.k.predictions) || length(with.k.predictions) != 1) {
        stop("'with.k.predictions' must be a single logical value")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value")
    }

    ## --- Prepare inputs for C++ (0-based indexing) ---
    adj.list.0based <- lapply(adj.list, function(v) as.integer(v - 1L))

    ## --- Call into C++ ---
    result <- .Call("S_graph_low_pass_smoother",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.integer(n.evectors.to.compute),
                    as.integer(min.num.eigenvectors),
                    as.integer(max.num.eigenvectors),
                    as.numeric(tau.factor),
                    as.numeric(radius.factor),
                    as.integer(laplacian.power),
                    as.integer(n.candidates),
                    as.logical(log.grid),
                    as.numeric(energy.threshold),
                    as.logical(with.k.predictions),
                    as.logical(verbose)
    )

    class(result) <- "graph_low_pass_smoother"
    result
}

#' Print Method for Graph Low-Pass Smoother
#'
#' @description
#' Prints a concise summary of a graph low-pass smoother object.
#'
#' @param x An object of class \code{graph_low_pass_smoother}
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.graph_low_pass_smoother <- function(x, ...) {
    cat("Graph Low-Pass Smoother\n")
    cat("=======================\n")
    cat(sprintf("Number of vertices:           %d\n", length(x$predictions)))
    cat(sprintf("Number of candidates tested:  %d\n", length(x$candidate.ks)))
    cat(sprintf("Method used:                  %s\n", x$used.method))
    cat(sprintf("Optimal k:                    %d\n",
                x$candidate.ks[switch(x$used.method,
                                    "eigengap" = x$opt.k.eigengap,
                                    "GCV" = x$opt.k.gcv,
                                    "spectral.energy" = x$opt.k.spectral.energy)]))
    cat("\nUse summary() for detailed diagnostics.\n")
    invisible(x)
}

#' Summary Method for Graph Low-Pass Smoother
#'
#' @description
#' Provides a detailed summary of the graph low-pass smoother results, including
#' optimization criteria values and selected parameters.
#'
#' @param object An object of class \code{graph_low_pass_smoother}
#' @param ... Additional arguments (ignored)
#'
#' @return A list of class \code{summary.graph_low_pass_smoother} containing:
#' \describe{
#'   \item{n.vertices}{Number of vertices in the graph}
#'   \item{n.candidates}{Number of candidate k values tested}
#'   \item{opt.k.eigengap}{Optimal k according to eigengap criterion}
#'   \item{opt.k.gcv}{Optimal k according to GCV criterion}
#'   \item{opt.k.spectral}{Optimal k according to spectral energy criterion}
#'   \item{used.method}{The method that was actually used}
#'   \item{gcv.range}{Range of GCV scores}
#'   \item{energy.at.optimal}{Spectral energy retained at optimal k}
#' }
#'
#' @export
summary.graph_low_pass_smoother <- function(object, ...) {
    # Get optimal k index based on method used
    opt.idx <- switch(object$used.method,
                      "eigengap" = object$opt.k.eigengap,
                      "GCV" = object$opt.k.gcv,
                      "spectral.energy" = object$opt.k.spectral.energy)

    out <- list(
        n.vertices = length(object$predictions),
        n.candidates = length(object$candidate.ks),
        opt.k.eigengap = object$candidate.ks[object$opt.k.eigengap],
        opt.k.gcv = object$candidate.ks[object$opt.k.gcv],
        opt.k.spectral = object$candidate.ks[object$opt.k.spectral.energy],
        used.method = object$used.method,
        optimal.k = object$candidate.ks[opt.idx],
        gcv.range = range(object$gcv.scores),
        energy.at.optimal = object$spectral.energy[opt.idx]
    )

    class(out) <- "summary.graph_low_pass_smoother"
    out
}

#' Print Summary of Graph Low-Pass Smoother
#'
#' @description
#' Prints a formatted summary of graph low-pass smoother results.
#'
#' @param x An object of class \code{summary.graph_low_pass_smoother}
#' @param digits Number of digits to display for numeric values
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.graph_low_pass_smoother <- function(x, digits = 4, ...) {
    cat("Summary of Graph Low-Pass Smoother\n")
    cat("==================================\n")
    cat(sprintf("Number of vertices:        %d\n", x$n.vertices))
    cat(sprintf("Candidates tested:         %d\n", x$n.candidates))
    cat("\nOptimal k by criterion:\n")
    cat(sprintf("  Eigengap:                %d\n", x$opt.k.eigengap))
    cat(sprintf("  GCV:                     %d\n", x$opt.k.gcv))
    cat(sprintf("  Spectral energy:         %d\n", x$opt.k.spectral))
    cat(sprintf("\nMethod chosen:             %s\n", x$used.method))
    cat(sprintf("Optimal k:                 %d\n", x$optimal.k))
    cat(sprintf("GCV score range:           [%.*f, %.*f]\n",
                digits, x$gcv.range[1], digits, x$gcv.range[2]))
    cat(sprintf("Energy retained:           %.*f%%\n",
                digits - 2, x$energy.at.optimal * 100))
    invisible(x)
}

#' Plot Method for Graph Low-Pass Smoother
#'
#' @description
#' Creates diagnostic plots for the graph low-pass smoother results.
#'
#' @param x An object of class \code{graph_low_pass_smoother}
#' @param which Integer or character specifying which plot to create:
#'   \itemize{
#'     \item 1 or "criteria": Plot all selection criteria
#'     \item 2 or "predictions": Plot original vs smoothed signal
#'     \item 3 or "eigenvalues": Plot eigenvalue spectrum
#'     \item 4 or "eigengaps": Plot eigengaps
#'   }
#' @param ... Additional graphical parameters passed to plot functions
#'
#' @return Invisibly returns the input object.
#'
#' @export
plot.graph_low_pass_smoother <- function(x, which = 1, ...) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    if (is.character(which)) {
        which <- match(which, c("criteria", "predictions", "eigenvalues", "eigengaps"))
        if (is.na(which)) stop("Invalid plot type specified")
    }

    if (which == 1) {
        # Plot all selection criteria
        par(mfrow = c(2, 2))

        # GCV scores
        plot(x$candidate.ks, x$gcv.scores, type = "b", pch = 19,
             xlab = "Number of eigenvectors (k)", ylab = "GCV Score",
             main = "GCV Criterion", ...)
        abline(v = x$candidate.ks[x$opt.k.gcv], col = "red", lty = 2)

        # Spectral energy
        plot(x$candidate.ks, x$spectral.energy, type = "b", pch = 19,
             xlab = "Number of eigenvectors (k)", ylab = "Cumulative Energy",
             main = "Spectral Energy", ...)
        abline(v = x$candidate.ks[x$opt.k.spectral.energy], col = "red", lty = 2)
        abline(h = x$spectral.energy[x$opt.k.spectral.energy], col = "red", lty = 3)

        # Eigengaps
        if (length(x$eigengaps) > 0) {
            plot(seq_along(x$eigengaps), x$eigengaps, type = "b", pch = 19,
                 xlab = "Eigenvalue index", ylab = "Eigengap",
                 main = "Eigengaps", ...)
            max.gap.idx <- which.max(x$eigengaps)
            abline(v = max.gap.idx, col = "red", lty = 2)
        }

        # Show optimal k
        plot.new()
        text(0.5, 0.5, paste("Method used:", x$used.method, "\n",
                             "Optimal k:", x$candidate.ks[switch(x$used.method,
                                                               "eigengap" = x$opt.k.eigengap,
                                                               "GCV" = x$opt.k.gcv,
                                                               "spectral.energy" = x$opt.k.spectral.energy)]),
             cex = 1.5)
    } else if (which == 2) {
        # Plot predictions
        n <- length(x$predictions)
        plot(1:n, x$y, pch = 19, col = "gray",
             xlab = "Vertex index", ylab = "Signal value",
             main = "Original vs Smoothed Signal", ...)
        lines(1:n, x$predictions, col = "red", lwd = 2)
        legend("topright", c("Original", "Smoothed"),
               col = c("gray", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
    } else if (which == 3) {
        # Plot eigenvalues
        plot(seq_along(x$evalues), x$evalues, type = "b", pch = 19,
             xlab = "Eigenvalue index", ylab = "Eigenvalue",
             main = "Laplacian Eigenvalue Spectrum", log = "y", ...)
    } else if (which == 4) {
        # Plot eigengaps
        if (length(x$eigengaps) > 0) {
            plot(seq_along(x$eigengaps), x$eigengaps, type = "b", pch = 19,
                 xlab = "Eigenvalue index", ylab = "Eigengap",
                 main = "Eigengaps (Consecutive Eigenvalue Differences)", ...)
            max.gap.idx <- which.max(x$eigengaps)
            abline(v = max.gap.idx, col = "red", lty = 2)
            text(max.gap.idx, x$eigengaps[max.gap.idx],
                 paste("Max gap at", max.gap.idx), pos = 4, col = "red")
        } else {
            plot.new()
            text(0.5, 0.5, "No eigengaps available", cex = 1.5)
        }
    } else {
        stop("Invalid value for 'which'. Must be 1-4.")
    }

    invisible(x)
}

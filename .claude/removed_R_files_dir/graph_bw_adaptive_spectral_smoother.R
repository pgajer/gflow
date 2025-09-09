#' Graph Bandwidth-Adaptive Spectral Smoother
#'
#' @description
#' Applies a graph-based locally linear smoother using spectral embedding and a
#' globally selected bandwidth (via leave-one-out cross-validation). The smoothing
#' is performed by fitting a locally weighted regression model over neighborhoods
#' defined via bandwidths in a Laplacian eigenvector space.
#'
#' @details
#' This function performs global bandwidth selection across all vertices to
#' optimize the LOOCV error, and then returns either:
#' \itemize{
#'   \item predictions at each vertex (based on optimal bandwidth),
#'   \item optionally, predictions for all bandwidths (for diagnostics),
#'   \item and vertex-specific minimum bandwidth constraints.
#' }
#'
#' The function uses spectral embedding based on the graph Laplacian to define
#' a metric space where local smoothing is performed. The optimal bandwidth is
#' selected by minimizing the leave-one-out cross-validation error across all
#' vertices simultaneously.
#'
#' @param adj.list List of integer vectors. The adjacency list of the graph
#'        (each entry gives neighbor vertex indices, 1-based).
#' @param weight.list List of numeric vectors. Edge weights matching the
#'        structure of \code{adj.list}. All weights must be non-negative.
#' @param y Numeric vector of responses at graph vertices. Missing values
#'        are not allowed.
#' @param n.evectors Integer. Number of Laplacian eigenvectors to use for
#'        spectral embedding. Default is 5.
#' @param n.bws Integer. Number of bandwidths to evaluate. Must be at least 20.
#'        Default is 20.
#' @param log.grid Logical. Use logarithmically spaced bandwidths? Otherwise,
#'        linear spacing. Default is \code{TRUE}.
#' @param min.bw.factor Numeric. Factor multiplied by graph diameter to get the
#'        minimum bandwidth. Must be positive and not greater than 0.15.
#'        Default is 0.025.
#' @param max.bw.factor Numeric. Factor multiplied by graph diameter to get the
#'        maximum bandwidth. Must be greater than \code{min.bw.factor} and not
#'        greater than 1.0. Default is 0.5.
#' @param dist.normalization.factor Numeric. Scale factor used in distance
#'        normalization before computing kernel weights. Must be at least 1.1.
#'        Default is 1.1.
#' @param kernel.type Integer. Code for kernel function to use for weight
#'        calculation:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#' @param precision Numeric. Precision used when generating the bandwidth grid.
#'        Must be positive. Default is 1e-5.
#' @param verbose Logical. If \code{TRUE}, prints progress during computation.
#'        Default is \code{FALSE}.
#' @param use.global.bw.grid Logical. If \code{TRUE}, use global candidate
#'        bandwidths schema. Default is \code{TRUE}.
#' @param with.bw.predictions Logical. If \code{TRUE}, returns prediction matrix
#'        for all bandwidths. Default is \code{TRUE}.
#' @param with.vertex.bw.errors Logical. If \code{TRUE}, returns per-vertex LOOCV
#'        error curves (one per bandwidth). Default is \code{FALSE}.
#'
#' @return A list with class \code{"graph_bw_adaptive_spectral_smoother"}
#'         containing the following components:
#' \describe{
#'   \item{predictions}{Numeric vector of smoothed predictions at each vertex
#'         (based on optimal bandwidth).}
#'   \item{bw.predictions}{Matrix of predictions for each vertex at each
#'         bandwidth (if \code{with.bw.predictions = TRUE}).}
#'   \item{bw.mean.abs.errors}{Mean absolute errors for each bandwidth.}
#'   \item{vertex.min.bws}{Minimum admissible bandwidth for each vertex.}
#'   \item{opt.bw.idx}{The (1-based) index of the globally optimal bandwidth.}
#'   \item{vertex.bw.errors}{Matrix of per-vertex LOOCV errors at each bandwidth
#'         (if \code{with.vertex.bw.errors = TRUE}).}
#' }
#'
#' @references
#' Belkin, M. and Niyogi, P. (2003). Laplacian eigenmaps for dimensionality
#' reduction and data representation. \emph{Neural Computation}, 15(6), 1373-1396.
#'
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. \emph{Journal of the American Statistical Association},
#' 74(368), 829-836.
#'
#' @examples
#' # Create a simple ring graph
#' n <- 20
#' adj.list <- lapply(1:n, function(i) {
#'   c(ifelse(i == 1, n, i - 1), ifelse(i == n, 1, i + 1))
#' })
#' weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
#'
#' # Generate noisy signal
#' y <- sin(2 * pi * (1:n) / n) + rnorm(n, sd = 0.2)
#'
#' # Apply smoother (simplified call for example)
#' ##\dontrun{
#' result <- graph.bw.adaptive.spectral.smoother(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   n.evectors = 3,
#'   n.bws = 20,
#'   verbose = TRUE
#' )
#'
#' # Examine results
#' print(result)
#' summary(result)
#' ##}
#'
#' @export
graph.bw.adaptive.spectral.smoother <- function(
    adj.list,
    weight.list,
    y,
    n.evectors = 5,
    min.bw.factor = 0.025,
    max.bw.factor = 0.5,
    n.bws = 20,
    log.grid = TRUE,
    kernel.type = 7L,
    dist.normalization.factor = 1.1,
    precision = 1e-5,
    use.global.bw.grid = TRUE,
    with.bw.predictions = TRUE,
    with.vertex.bw.errors = FALSE,
    verbose = FALSE
) {
    ## Input validation

    ## Check graph structure
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("'adj.list' and 'weight.list' must both be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }

    n.vertices <- length(adj.list)

    if (n.vertices < 2) {
        stop("Graph must have at least 2 vertices")
    }

    ## Check each vertex's edge/weight consistency and validate indices
    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("At vertex %d: lengths of adjacency and weight lists must be equal", i))
        }

        ## Check for valid vertex indices
        if (length(adj.list[[i]]) > 0) {
            if (!all(adj.list[[i]] %in% 1:n.vertices)) {
                stop(sprintf("At vertex %d: invalid neighbor indices (must be between 1 and %d)",
                           i, n.vertices))
            }

            ## Check for self-loops
            if (i %in% adj.list[[i]]) {
                warning(sprintf("At vertex %d: self-loop detected and will be ignored", i))
            }

            ## Check for non-negative weights
            if (any(weight.list[[i]] < 0)) {
                stop(sprintf("At vertex %d: all weights must be non-negative", i))
            }

            ## Check for NA/NaN/Inf in weights
            if (any(!is.finite(weight.list[[i]]))) {
                stop(sprintf("At vertex %d: weights contain NA, NaN, or infinite values", i))
            }
        }
    }

    ## Response vector checks
    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector")
    }

    if (length(y) != n.vertices) {
        stop(sprintf("Length of 'y' (%d) must match the number of vertices (%d)",
                   length(y), n.vertices))
    }

    if (any(!is.finite(y))) {
        stop("'y' contains NA, NaN, or infinite values")
    }

    ## Bandwidth parameters
    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1) {
        stop("'min.bw.factor' must be a single numeric value")
    }

    if (min.bw.factor <= 0 || min.bw.factor > 0.15) {
        stop("'min.bw.factor' must be positive and not greater than 0.15")
    }

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1) {
        stop("'max.bw.factor' must be a single numeric value")
    }

    if (max.bw.factor <= min.bw.factor || max.bw.factor > 1.0) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor' and not greater than 1.0")
    }

    if (!is.numeric(n.bws) || length(n.bws) != 1) {
        stop("'n.bws' must be a single numeric value")
    }

    n.bws <- as.integer(n.bws)
    if (n.bws < 20) {
        stop("'n.bws' must be at least 20")
    }

    ## Grid type and embedding size
    if (!is.logical(log.grid) || length(log.grid) != 1) {
        stop("'log.grid' must be a single logical value")
    }

    if (!is.numeric(n.evectors) || length(n.evectors) != 1) {
        stop("'n.evectors' must be a single numeric value")
    }

    n.evectors <- as.integer(n.evectors)
    if (n.evectors < 1 || n.evectors > n.vertices) {
        stop(sprintf("'n.evectors' must be between 1 and %d", n.vertices))
    }

    ## Kernel type
    if (!is.numeric(kernel.type) || length(kernel.type) != 1) {
        stop("'kernel.type' must be a single numeric value")
    }

    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel.type' must be one of: 1 (Epanechnikov), 2 (Triangular), ",
             "4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    ## Distance normalization factor
    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1) {
        stop("'dist.normalization.factor' must be a single numeric value")
    }

    if (!is.finite(dist.normalization.factor) || dist.normalization.factor < 1.1) {
        stop("'dist.normalization.factor' must be finite and at least 1.1")
    }

    ## Precision
    if (!is.numeric(precision) || length(precision) != 1) {
        stop("'precision' must be a single numeric value")
    }

    if (!is.finite(precision) || precision <= 0) {
        stop("'precision' must be a finite positive number")
    }

    ## Logical flags
    if (!is.logical(use.global.bw.grid) || length(use.global.bw.grid) != 1) {
        stop("'use.global.bw.grid' must be a single logical value")
    }

    if (!is.logical(with.bw.predictions) || length(with.bw.predictions) != 1) {
        stop("'with.bw.predictions' must be a single logical value")
    }

    if (!is.logical(with.vertex.bw.errors) || length(with.vertex.bw.errors) != 1) {
        stop("'with.vertex.bw.errors' must be a single logical value")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value")
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call(C_graph_bw_adaptive_spectral_smoother,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.integer(n.evectors),
                    as.numeric(min.bw.factor),
                    as.numeric(max.bw.factor),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.integer(kernel.type),
                    as.numeric(dist.normalization.factor),
                    as.numeric(precision),
                    as.logical(use.global.bw.grid),
                    as.logical(with.bw.predictions),
                    as.logical(with.vertex.bw.errors),
                    as.logical(verbose))

    ## Ensure consistent naming in return value
    names(result)[names(result) == "opt_bw_idx"] <- "opt.bw.idx"

    class(result) <- "graph_bw_adaptive_spectral_smoother"

    return(result)
}


#' Print Method for Graph Bandwidth-Adaptive Spectral Smoother
#'
#' @description
#' Prints a concise summary of the output from
#' \code{\link{graph.bw.adaptive.spectral.smoother}}, including the number of
#' vertices, number of bandwidths evaluated, and the index of the optimal
#' bandwidth.
#'
#' @param x An object of class \code{"graph_bw_adaptive_spectral_smoother"} as
#'        returned by \code{graph.bw.adaptive.spectral.smoother}.
#' @param ... Additional arguments passed to or from other methods (currently
#'        ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{graph.bw.adaptive.spectral.smoother}},
#' \code{\link{summary.graph_bw_adaptive_spectral_smoother}}
#'
#' @examples
#' # See examples in graph.bw.adaptive.spectral.smoother
#'
#' @export
#' @method print graph_bw_adaptive_spectral_smoother
print.graph_bw_adaptive_spectral_smoother <- function(x, ...) {
    cat("Graph Bandwidth-Adaptive Spectral Smoother\n")
    cat("==========================================\n")
    cat(sprintf("Number of vertices: %d\n", length(x$predictions)))
    cat(sprintf("Number of bandwidths evaluated: %d\n", length(x$bw.mean.abs.errors)))
    cat(sprintf("Optimal bandwidth index: %d\n", x$opt.bw.idx))
    cat(sprintf("Minimum mean absolute error: %.4f\n",
                x$bw.mean.abs.errors[x$opt.bw.idx]))
    cat("\nCall summary() for detailed diagnostics.\n")
    invisible(x)
}


#' Summary Method for Graph Bandwidth-Adaptive Spectral Smoother
#'
#' @description
#' Computes and returns a summary of the key results from
#' \code{\link{graph.bw.adaptive.spectral.smoother}}, including prediction and
#' error summaries across bandwidths and vertices.
#'
#' @param object An object of class \code{"graph_bw_adaptive_spectral_smoother"}
#'        as returned by \code{graph.bw.adaptive.spectral.smoother}.
#' @param ... Additional arguments passed to or from other methods (currently
#'        ignored).
#'
#' @return
#' An object of class \code{"summary.graph_bw_adaptive_spectral_smoother"}
#' containing:
#' \describe{
#'   \item{n.vertices}{Number of vertices in the graph.}
#'   \item{n.bandwidths}{Number of bandwidths evaluated.}
#'   \item{opt.bw.idx}{1-based index of optimal bandwidth.}
#'   \item{opt.bw.error}{Mean absolute error at optimal bandwidth.}
#'   \item{bw.error.summary}{Summary statistics for the bandwidth-wise mean
#'         absolute errors (min, 1st quartile, median, mean, 3rd quartile, max).}
#'   \item{vertex.min.bw.summary}{Summary statistics for the per-vertex minimum
#'         bandwidths (if available).}
#'   \item{prediction.summary}{Summary statistics for the final predictions.}
#' }
#'
#' @seealso
#' \code{\link{graph.bw.adaptive.spectral.smoother}},
#' \code{\link{print.summary.graph_bw_adaptive_spectral_smoother}}
#'
#' @examples
#' # See examples in graph.bw.adaptive.spectral.smoother
#'
#' @export
#' @method summary graph_bw_adaptive_spectral_smoother
summary.graph_bw_adaptive_spectral_smoother <- function(object, ...) {
    out <- list(
        n.vertices = length(object$predictions),
        n.bandwidths = length(object$bw.mean.abs.errors),
        opt.bw.idx = object$opt.bw.idx,
        opt.bw.error = object$bw.mean.abs.errors[object$opt.bw.idx],
        bw.error.summary = summary(object$bw.mean.abs.errors),
        prediction.summary = summary(object$predictions)
    )

    ## Add vertex minimum bandwidth summary if available
    if (!is.null(object$vertex.min.bws)) {
        out$vertex.min.bw.summary <- summary(object$vertex.min.bws)
    }

    class(out) <- "summary.graph_bw_adaptive_spectral_smoother"
    out
}


#' Print Method for Summary of Graph Bandwidth-Adaptive Spectral Smoother
#'
#' @description
#' Prints a formatted summary of a
#' \code{"summary.graph_bw_adaptive_spectral_smoother"} object, displaying key
#' statistics about the smoothing results.
#'
#' @param x An object of class \code{"summary.graph_bw_adaptive_spectral_smoother"}
#'        as returned by \code{summary.graph_bw_adaptive_spectral_smoother}.
#' @param digits Integer. Number of significant digits to display. Default is 4.
#' @param ... Additional arguments passed to or from other methods (currently
#'        ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{summary.graph_bw_adaptive_spectral_smoother}},
#' \code{\link{graph.bw.adaptive.spectral.smoother}}
#'
#' @examples
#' # See examples in graph.bw.adaptive.spectral.smoother
#'
#' @export
#' @method print summary.graph_bw_adaptive_spectral_smoother
print.summary.graph_bw_adaptive_spectral_smoother <- function(x,
                                                              digits = 4,
                                                              ...) {
    cat("\nSummary of Graph Bandwidth-Adaptive Spectral Smoother\n")
    cat("=====================================================\n\n")

    cat("Graph structure:\n")
    cat(sprintf("  Number of vertices: %d\n", x$n.vertices))
    cat("\n")

    cat("Bandwidth selection:\n")
    cat(sprintf("  Number of bandwidths evaluated: %d\n", x$n.bandwidths))
    cat(sprintf("  Optimal bandwidth index: %d\n", x$opt.bw.idx))
    cat(sprintf("  Mean absolute error (optimal): %.*f\n", digits, x$opt.bw.error))
    cat("\n")

    cat("Bandwidth error distribution:\n")
    print(x$bw.error.summary, digits = digits)
    cat("\n")

    cat("Prediction summary:\n")
    print(x$prediction.summary, digits = digits)

    if (!is.null(x$vertex.min.bw.summary)) {
        cat("\nVertex minimum bandwidth summary:\n")
        print(x$vertex.min.bw.summary, digits = digits)
    }

    invisible(x)
}


#' Plot Method for Graph Bandwidth-Adaptive Spectral Smoother
#'
#' @description
#' Creates diagnostic plots for the results of
#' \code{\link{graph.bw.adaptive.spectral.smoother}}. Can display bandwidth
#' selection curves and optionally per-vertex error curves.
#'
#' @param x An object of class \code{"graph_bw_adaptive_spectral_smoother"} as
#'        returned by \code{graph.bw.adaptive.spectral.smoother}.
#' @param which Integer or character string specifying which plot to create:
#'        \itemize{
#'          \item 1 or \code{"bw.errors"}: Mean absolute error vs bandwidth index
#'          \item 2 or \code{"vertex.errors"}: Per-vertex errors (if available)
#'        }
#'        Default is 1.
#' @param main Character string. Main title for the plot. If \code{NULL}, a
#'        default title is used.
#' @param xlab Character string. Label for x-axis. If \code{NULL}, a default
#'        label is used.
#' @param ylab Character string. Label for y-axis. If \code{NULL}, a default
#'        label is used.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{graph.bw.adaptive.spectral.smoother}}
#'
#' @examples
#' # See examples in graph.bw.adaptive.spectral.smoother
#'
#' @export
#' @method plot graph_bw_adaptive_spectral_smoother
plot.graph_bw_adaptive_spectral_smoother <- function(x,
                                                     which = 1,
                                                     main = NULL,
                                                     xlab = NULL,
                                                     ylab = NULL,
                                                     ...) {
    which <- if (is.character(which)) {
        match(which, c("bw.errors", "vertex.errors"))
    } else {
        as.integer(which)
    }

    if (is.na(which) || !which %in% 1:2) {
        stop("'which' must be 1, 2, 'bw.errors', or 'vertex.errors'")
    }

    if (which == 1) {
        ## Bandwidth error plot
        if (is.null(main)) main <- "Bandwidth Selection Curve"
        if (is.null(xlab)) xlab <- "Bandwidth Index"
        if (is.null(ylab)) ylab <- "Mean Absolute Error"

        n.bws <- length(x$bw.mean.abs.errors)
        plot(1:n.bws, x$bw.mean.abs.errors,
             type = "b", pch = 19, cex = 0.7,
             main = main, xlab = xlab, ylab = ylab, ...)

        ## Highlight optimal bandwidth
        points(x$opt.bw.idx, x$bw.mean.abs.errors[x$opt.bw.idx],
               pch = 19, cex = 1.5, col = "red")
        abline(v = x$opt.bw.idx, lty = 2, col = "red")

        ## Add legend
        legend("topright",
               legend = c("Error curve", "Optimal bandwidth"),
               pch = c(19, 19),
               col = c("black", "red"),
               pt.cex = c(0.7, 1.5))

    } else if (which == 2) {
        ## Per-vertex error plot (if available)
        if (is.null(x$vertex.bw.errors)) {
            stop("Per-vertex errors not available. Run with 'with.vertex.bw.errors = TRUE'")
        }

        if (is.null(main)) main <- "Per-Vertex Error Curves"
        if (is.null(xlab)) xlab <- "Bandwidth Index"
        if (is.null(ylab)) ylab <- "LOOCV Error"

        n.bws <- ncol(x$vertex.bw.errors)
        n.vertices <- nrow(x$vertex.bw.errors)

        ## Plot first vertex
        plot(1:n.bws, x$vertex.bw.errors[1, ],
             type = "l", col = rgb(0, 0, 0, 0.2),
             main = main, xlab = xlab, ylab = ylab,
             ylim = range(x$vertex.bw.errors, na.rm = TRUE), ...)

        ## Add remaining vertices
        for (i in 2:min(n.vertices, 100)) {  # Limit to 100 vertices for clarity
            lines(1:n.bws, x$vertex.bw.errors[i, ], col = rgb(0, 0, 0, 0.2))
        }

        ## Add mean error curve
        lines(1:n.bws, x$bw.mean.abs.errors, col = "red", lwd = 3)

        ## Highlight optimal bandwidth
        abline(v = x$opt.bw.idx, lty = 2, col = "red")

        if (n.vertices > 100) {
            mtext(sprintf("(Showing %d of %d vertices)", 100, n.vertices),
                  side = 3, line = 0, cex = 0.8)
        }
    }

    invisible(x)
}

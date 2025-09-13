#' Graph Kernel Smoother
#'
#' Smooth a scalar signal on a graph using kernel weights along edges.
#' This is a thin user-facing wrapper that forwards to the Rcpp
#' implementation `Rcpp_graph_kernel_smoother()`.
#'
#' @param adj A list of integer vectors; the adjacency list of the graph.
#'   Each element \code{adj[[i]]} contains the (1-based) neighbor indices of vertex \code{i}.
#' @param w A list of numeric vectors of the same shape as \code{adj}, or \code{NULL}.
#'   If provided, \code{w[[i]]} must have length \code{length(adj[[i]])} and
#'   contains edge weights for the neighbors of vertex \code{i}. Use \code{NULL}
#'   for unweighted graphs (defaults to weight 1 per edge).
#' @param y A numeric vector of length \code{length(adj)}; the signal to be smoothed.
#' @param bandwidth Integer bandwidth (method-specific; typically controls kernel
#'   neighborhood/scale). Must be \eqn{\ge 1}.
#' @param with_details Logical; if \code{TRUE}, include additional diagnostics
#'   (e.g., per-vertex bandwidth usage, iteration counts). Defaults to \code{FALSE}.
#'
#' @return A named list with components:
#' \itemize{
#'   \item \code{fitted} — numeric vector of length \code{length(y)} with smoothed values.
#'   \item \code{bandwidth} — integer bandwidth actually used.
#'   \item \code{details} — (optional) a list of diagnostics, present when \code{with_details = TRUE}.
#' }
#'
#' @section Input requirements:
#' \itemize{
#'   \item \code{adj} must be a list the same length as \code{y}; neighbors are 1-based indices.
#'   \item If \code{w} is not \code{NULL}, it must be a list the same length as \code{adj} with
#'         \code{length(w[[i]]) == length(adj[[i]])} for all \code{i}.
#' }
#'
#' @section Notes:
#' This wrapper only orchestrates I/O and delegates computation to the Rcpp routine
#' \code{Rcpp_graph_kernel_smoother()}. For reproducibility, ensure any stochastic
#' components in the underlying method are deterministically seeded upstream.
#'
#' @examples
#' \dontrun{
#' # Tiny chain graph with unit weights
#' n <- 5L
#' adj <- vector("list", n)
#' for (i in seq_len(n)) {
#'   adj[[i]] <- unique(c(if (i > 1) i - 1L else integer(), if (i < n) i + 1L else integer()))
#' }
#' w <- lapply(adj, function(nei) rep(1, length(nei)))
#' y <- as.numeric(1:n)
#'
#' fit <- graph.kernel.smoother(adj, w, y, bandwidth = 2L, with_details = FALSE)
#' str(fit)
#' }
#'
#' @export
graph.kernel.smoother <- function(adj, w, y, bandwidth, with_details = FALSE) {
  # Basic validation at the R layer (fast errors, helpful messages)
  if (!is.list(adj)) stop("`adj` must be a list of integer vectors.")
  n <- length(adj)
  if (!is.numeric(y) || length(y) != n) {
    stop("`y` must be numeric and length(y) must equal length(adj).")
  }
  if (!is.null(w)) {
    if (!is.list(w) || length(w) != n) stop("`w` must be NULL or a list matching `adj` in length.")
    for (i in seq_len(n)) {
      if (length(w[[i]]) && length(w[[i]]) != length(adj[[i]])) {
        stop(sprintf("`w[[%d]]` must have the same length as `adj[[%d]]`.", i, i))
      }
    }
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || bandwidth < 1) {
    stop("`bandwidth` must be a single integer >= 1.")
  }

  # Delegate to the Rcpp implementation (exported via Rcpp::compileAttributes)
  Rcpp_graph_kernel_smoother(adj, w, y, as.integer(bandwidth), isTRUE(with_details))
}


#' Graph-Based Kernel Smoother with Buffer Zone Cross-Validation
#'
#' @description
#' Performs graph-based locally weighted smoothing (LOWESS-like) with adaptive
#' bandwidth selection using spatially-aware cross-validation. The method
#' implements buffer zones around test vertices during cross-validation to
#' prevent spatial autocorrelation from biasing bandwidth selection.
#'
#' @details
#' This function implements a graph-based extension of LOWESS (degree-0, i.e.,
#' locally weighted average) with spatially-stratified cross-validation using
#' buffer zones. The algorithm proceeds as follows:
#'
#' \enumerate{
#'   \item \strong{Spatial fold creation}: Creates a maximal packing of vertices
#'         to serve as fold seed points, then assigns all vertices to the nearest
#'         seed point to form spatially coherent folds.
#'   \item \strong{Buffer zone construction}: For each fold, creates a buffer zone
#'         around test vertices by excluding vertices within a specified graph
#'         distance (hops) from the training set.
#'   \item \strong{Cross-validation}: For each candidate bandwidth, performs
#'         cross-validation across the spatially-separated folds.
#'   \item \strong{Bandwidth selection}: Selects the bandwidth that minimizes
#'         the cross-validation error.
#'   \item \strong{Final fitting}: Fits the final model using all data with the
#'         optimal bandwidth.
#' }
#'
#' The kernel smoother at vertex \eqn{i} computes:
#' \deqn{\hat{y}_i = \sum_{j \in N(i, h)} w_{ij} y_j / \sum_{j \in N(i, h)} w_{ij}}
#' where \eqn{N(i, h)} is the neighborhood of vertex \eqn{i} within bandwidth \eqn{h},
#' and \eqn{w_{ij}} are kernel weights based on graph distance.
#'
#' @param adj.list A list of integer vectors representing the adjacency list of
#'   the graph. Each element \code{adj.list[[i]]} contains the indices of
#'   vertices adjacent to vertex \code{i}. Vertices should be indexed from 1 to n.
#' @param weight.list A list of numeric vectors with edge weights corresponding
#'   to adjacencies. Each element \code{weight.list[[i]][j]} is the weight of
#'   the edge from vertex \code{i} to \code{adj.list[[i]][j]}. Must be positive values.
#' @param y A numeric vector of response values at each vertex. Length must
#'   match the number of vertices in the graph.
#' @param min.bw.factor Minimum bandwidth as a factor of graph diameter. Must
#'   be between 0 and 1. Default is 0.01.
#' @param max.bw.factor Maximum bandwidth as a factor of graph diameter. Must
#'   be between \code{min.bw.factor} and 1. Default is 0.5.
#' @param n.bws Number of bandwidths to test. Must be at least 2. Default is 20.
#' @param log.grid Logical. If \code{TRUE}, use logarithmic spacing for
#'   bandwidth grid; if \code{FALSE}, use linear spacing. Default is \code{TRUE}.
#' @param vertex.hbhd.min.size Integer. Minimum number of vertices required in
#'   any neighborhood. Must be between 1 and \eqn{10\%} of total vertices. Default is 1.
#' @param kernel.type Integer specifying the kernel function:
#'   \itemize{
#'     \item 1: Uniform kernel
#'     \item 2: Triangular kernel
#'     \item 3: Epanechnikov kernel
#'     \item 4: Quartic (biweight) kernel
#'     \item 5: Triweight kernel
#'     \item 6: Gaussian kernel
#'     \item 7: Tricube kernel (default)
#'   }
#' @param dist.normalization.factor Positive factor for normalizing distances
#'   in kernel weight computation. Values > 1 result in wider effective
#'   bandwidths. Default is 1.1.
#' @param use.uniform.weights Logical. If \code{TRUE}, use uniform weights
#'   instead of kernel-based weights. Default is \code{FALSE}.
#' @param buffer.hops Integer. Number of hops for buffer zone around test
#'   vertices during cross-validation. Must be non-negative. Default is 2.
#' @param auto.buffer.hops Logical. If \code{TRUE}, automatically determine
#'   optimal buffer size based on spatial autocorrelation analysis. Default
#'   is \code{TRUE}.
#' @param n.folds Integer. Number of cross-validation folds. Default is 5.
#'   The effective range is \eqn{[2, \lfloor n/2 \rfloor]}, where \eqn{n} is the
#'   number of vertices in the input graph. Values outside this range are
#'   automatically clamped: if \code{n.folds < 2}, it is reset to 2; if
#'   \code{n.folds > floor(n/2)}, it is reset to \eqn{\lfloor n/2 \rfloor}.
#'   A message is issued when values are clamped if \code{verbose = TRUE}.
#' @param with.bw.predictions Logical. If \code{TRUE}, compute and store
#'   predictions for all tested bandwidths. Default is \code{FALSE}.
#' @param precision Numeric. Precision for bandwidth grid computation. Must be
#'   positive. Default is 1e-6.
#' @param verbose Logical. If \code{TRUE}, print progress information during
#'   computation. Default is \code{FALSE}.
#'
#' @return An object of class \code{"graph_kernel_smoother"}, which is a list
#'   containing:
#'   \describe{
#'     \item{\code{predictions}}{Numeric vector of smoothed predictions at each vertex.}
#'     \item{\code{bw_predictions}}{If \code{with.bw.predictions = TRUE}, a list
#'       of prediction vectors for each tested bandwidth; otherwise \code{NULL}.}
#'     \item{\code{bw_errors}}{Numeric vector of cross-validation errors for each
#'       tested bandwidth.}
#'     \item{\code{bws}}{Numeric vector of bandwidth values tested.}
#'     \item{\code{opt_bw}}{Numeric. The optimal bandwidth value selected by
#'       cross-validation.}
#'     \item{\code{opt_bw_idx}}{Integer. Index of the optimal bandwidth in
#'       \code{bws}.}
#'     \item{\code{buffer_hops_used}}{Integer. Number of hops used for the buffer
#'       zone (useful when \code{auto.buffer.hops = TRUE}).}
#'   }
#'
#' @examples
#' \dontrun{
#' n.pts <- 100
#' gm <- generate.1d.gaussian.mixture(
#'     n.points = n.pts,
#'     x.knot = c(0, 10),
#'     y.knot = c(10, 2.5),
#'     sd.knot = 1.5,
#'     x.offset = 3)
#' x <- sort(runif(n.pts, min = min(gm$x), max = max(gm$x)))
#' x.graph <- create.bi.kNN.chain.graph(k = 1, x = x, y = gm$y)
#'
#' y.smooth <- approx(gm$x, gm$y, xout = x)$y
#' sigma <- 1
#' eps <- rnorm(n.pts, 0, sigma)
#' y <- y.smooth + eps
#'
#' g <- ggraph(x.graph$adj.list, x.graph$edge.lengths)
#' plot(g, y.smooth)
#'
#' plot(gm$x, gm$y, type = "l", las = 1, ylim = range(y), col = "red", xlab = "x", ylab = "y")
#' points(x, y)
#' legend("topright", legend = c("y.smooth", "y"), lty = c(1,NA), pch = c(NA,1),
#' col = c("red", "black"), inset = 0.1)
#'
#' gks.res <- graph.kernel.smoother(x.graph$adj.list,
#'                                  x.graph$edge.lengths,
#'                                  y,
#'                                  min.bw.factor = 0.025,
#'                                  max.bw.factor = 0.5,
#'                                  n.bws = 20,
#'                                  log.grid = TRUE,
#'                                  vertex.hbhd.min.size = 3,
#'                                  dist.normalization.factor = 1.1,
#'                                  use.uniform.weights = FALSE,
#'                                  buffer.hops = 1,
#'                                  auto.buffer.hops = FALSE,
#'                                  kernel.type = 7L,
#'                                  n.folds = 10,
#'                                  with.bw.predictions = TRUE,
#'                                  verbose = TRUE)
#'
#' # View results
#' print(gks.res)
#' summary(gks.res)
#'
#' # Computing Mean Absolute and Mean Squared Errors
#' mae <- c()
#' for (i in seq(ncol(gks.res$bw_predictions))) {
#'     mae[i] <- mean(abs(y.smooth - gks.res$bw_predictions[,i]))
#' }
#' plot(mae, las = 1, type = 'b', xlab = "Bandwidth indices", ylab = "Mean Absolute Error")
#' which.min(mae)
#' gks.res$opt_bw_idx
#'
#' plot(gks.res$bw_mean_abs_errors, las = 1, type = "b")
#' abline(v = gks.res$opt_bw_idx, lty = 2)
#' which.min(gks.res$bw_mean_abs_errors)
#'
#' plot(gm$x, gm$y, type = "l", las = 1, ylim = range(y), col = "red")
#' points(x, y)
#' lines(x, gks.res$predictions, col = "blue")
#' }
#'
S.version.graph.kernel.smoother <- function(
    adj.list,
    weight.list,
    y,
    min.bw.factor = 0.01,
    max.bw.factor = 0.5,
    n.bws = 20,
    log.grid = TRUE,
    vertex.hbhd.min.size = 1,
    dist.normalization.factor = 1.1,
    use.uniform.weights = FALSE,
    buffer.hops = 2,
    auto.buffer.hops = TRUE,
    kernel.type = 7L,
    n.folds = 5,
    with.bw.predictions = FALSE,
    precision = 1e-6,
    verbose = FALSE
) {
    ## Input validation with informative error messages
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
    }

    n.vertices <- length(adj.list)
    if (n.vertices < 4) {
        stop("The graph must have at least 4 vertices")
    }

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    if (length(y) != n.vertices) {
        stop("Length of 'y' must match the number of vertices in the graph", call. = FALSE)
    }

    if (anyNA(y)) {
        stop("'y' cannot contain NA values", call. = FALSE)
    }

    ## Validate adjacency and weight lists
    for (i in seq_along(adj.list)) {
        if (!is.numeric(adj.list[[i]]) && !is.integer(adj.list[[i]])) {
            stop(sprintf("adj.list[[%d]] must be numeric or integer", i), call. = FALSE)
        }

        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Length mismatch between adj.list[[%d]] and weight.list[[%d]]", i, i),
                 call. = FALSE)
        }

        if (!is.numeric(weight.list[[i]])) {
            stop(sprintf("weight.list[[%d]] must be numeric", i), call. = FALSE)
        }

        if (any(weight.list[[i]] <= 0)) {
            stop(sprintf("All weights in weight.list[[%d]] must be positive", i), call. = FALSE)
        }

        if (any(adj.list[[i]] < 1 | adj.list[[i]] > n.vertices)) {
            stop(sprintf("Invalid vertex indices in adj.list[[%d]]: indices must be between 1 and %d",
                         i, n.vertices), call. = FALSE)
        }
    }

    ## Validate bandwidth parameters
    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1 ||
        min.bw.factor <= 0 || min.bw.factor > 1) {
        stop("'min.bw.factor' must be a single number between 0 and 1", call. = FALSE)
    }

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1 ||
        max.bw.factor <= min.bw.factor || max.bw.factor > 1) {
        stop("'max.bw.factor' must be a single number between 'min.bw.factor' and 1",
             call. = FALSE)
    }

    if (!is.numeric(n.bws) || length(n.bws) != 1 || n.bws < 2) {
        stop("'n.bws' must be a single integer >= 2", call. = FALSE)
    }
    n.bws <- as.integer(n.bws)

    if (!is.logical(log.grid) || length(log.grid) != 1) {
        stop("'log.grid' must be a single logical value", call. = FALSE)
    }

    ## Validate other parameters
    if (!is.numeric(vertex.hbhd.min.size) || length(vertex.hbhd.min.size) != 1 ||
        vertex.hbhd.min.size < 1 || vertex.hbhd.min.size > max(c(as.integer(0.1 * n.vertices), 1))) {
        stop("'vertex.hbhd.min.size' must be a positive integer not greater than 10% of total vertices",
             call. = FALSE)
    }
    vertex.hbhd.min.size <- as.integer(vertex.hbhd.min.size)

    if (!is.numeric(kernel.type) || length(kernel.type) != 1 ||
        !(kernel.type %in% 1:7)) {
        stop("'kernel.type' must be an integer between 1 and 7", call. = FALSE)
    }
    kernel.type <- as.integer(kernel.type)

    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1 ||
        dist.normalization.factor <= 0) {
        stop("'dist.normalization.factor' must be a positive number", call. = FALSE)
    }

    if (!is.logical(use.uniform.weights) || length(use.uniform.weights) != 1) {
        stop("'use.uniform.weights' must be a single logical value", call. = FALSE)
    }

    if (!is.numeric(buffer.hops) || length(buffer.hops) != 1 || buffer.hops < 0) {
        stop("'buffer.hops' must be a non-negative integer", call. = FALSE)
    }
    buffer.hops <- as.integer(buffer.hops)

    if (!is.logical(auto.buffer.hops) || length(auto.buffer.hops) != 1) {
        stop("'auto.buffer.hops' must be a single logical value", call. = FALSE)
    }


    ## derive allowable max folds for this graph
    n.max.folds <- max(2L, floor(n.vertices / 2L))

    ## sanitize n.folds
    if (length(n.folds) != 1L || is.na(n.folds) || !is.finite(n.folds)) {
        stop("'n.folds' must be a single finite number", call. = FALSE)
    }
    n.folds <- as.integer(n.folds)

    ## auto-clamp to valid range, with a gentle message (or use warning())
    if (n.folds < 2L) {
        if (verbose) message("n.folds < 2; using 2")
        n.folds <- 2L
    } else if (n.folds > n.max.folds) {
        if (verbose) message(sprintf("n.folds > n/2; using %d", n.max.folds))
        n.folds <- n.max.folds
    }

    if (!is.numeric(n.folds) || length(n.folds) != 1 || n.folds < 2 ||
        n.folds > n.vertices / 2) {
        stop("'n.folds' must be an integer between 2 and n/2", call. = FALSE)
    }

    if (!is.logical(with.bw.predictions) || length(with.bw.predictions) != 1) {
        stop("'with.bw.predictions' must be a single logical value", call. = FALSE)
    }

    if (!is.numeric(precision) || length(precision) != 1 || precision <= 0) {
        stop("'precision' must be a positive number", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value", call. = FALSE)
    }

    ## Convert to 0-based indexing for C++ code
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call(S_graph_kernel_smoother,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.double(min.bw.factor),
                    as.double(max.bw.factor),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.integer(vertex.hbhd.min.size),
                    as.integer(kernel.type),
                    as.double(dist.normalization.factor),
                    as.logical(use.uniform.weights),
                    as.integer(buffer.hops),
                    as.logical(auto.buffer.hops),
                    as.integer(n.folds),
                    as.logical(with.bw.predictions),
                    as.double(precision),
                    as.logical(verbose))
        
    ## Add class attribute to result for S3 methods
    class(result) <- c("graph_kernel_smoother", "list")
    
    return(result)
}

#' Estimate Optimal Bandwidth Using Local Extrema Count Elbow Method
#'
#' @description
#' Estimates the optimal bandwidth for kernel smoothing by analyzing the
#' relationship between bandwidth and the number of local extrema. The method
#' identifies an "elbow point" in this relationship using piecewise linear
#' regression and adjusts it based on model residuals.
#'
#' @details
#' The number of local extrema in a kernel-smoothed function typically decreases
#' as bandwidth increases, often exhibiting an elbow-shaped curve. This function:
#'
#' \enumerate{
#'   \item Fits a piecewise linear model to the extrema counts
#'   \item Identifies the breakpoint (elbow) in the relationship
#'   \item Adds a specified multiple of the residual standard deviation to the
#'         breakpoint location
#' }
#'
#' The optimal bandwidth is often slightly beyond the elbow point, as the exact
#' elbow typically represents the transition from undersmoothing to appropriate
#' smoothing.
#'
#' @param extrema.counts A numeric vector of local extrema counts corresponding
#'   to different bandwidths, ordered from smallest to largest bandwidth. Must
#'   contain at least 5 values.
#' @param sd.multiplier The number of residual standard deviations to add to
#'   the elbow point. Must be non-negative. Default is 1.0.
#' @param plot.results Logical. If \code{TRUE}, produces a diagnostic plot
#'   showing the fitted model and selected bandwidth. Default is \code{FALSE}.
#'
#' @return An integer representing the index of the estimated optimal bandwidth
#'   in the input vector.
#'
#' @examples
#' # Simulated extrema counts for increasing bandwidths
#' extrema.counts <- c(150, 142, 128, 95, 72, 58, 45, 38, 35, 33, 32, 31)
#'
#' # Estimate optimal bandwidth
#' opt.idx <- estimate.optimal.bandwidth.from.extrema.elbow(
#'   extrema.counts,
#'   sd.multiplier = 1.5,
#'   plot.results = TRUE
#' )
#'
#' print(paste("Optimal bandwidth index:", opt.idx))
#'
#' @importFrom segmented segmented
#' @importFrom stats lm predict sd
#' @importFrom graphics plot lines abline legend par
#'
#' @seealso \code{\link{graph.kernel.smoother}} for the main smoothing function.
#'
#' @references
#' Muggeo, V. M. (2003). Estimating regression models with unknown break-points.
#' \emph{Statistics in Medicine}, 22(19), 3055-3071.
#'
#' @export
estimate.optimal.bandwidth.from.extrema.elbow <- function(
    extrema.counts,
    sd.multiplier = 1.0,
    plot.results = FALSE
) {
    ## Input validation
    if (!is.numeric(extrema.counts) || !is.vector(extrema.counts)) {
        stop("'extrema.counts' must be a numeric vector", call. = FALSE)
    }

    if (length(extrema.counts) < 5) {
        stop("'extrema.counts' must contain at least 5 values", call. = FALSE)
    }

    if (anyNA(extrema.counts)) {
        stop("'extrema.counts' cannot contain NA values", call. = FALSE)
    }

    if (!is.numeric(sd.multiplier) || length(sd.multiplier) != 1 ||
        sd.multiplier < 0) {
        stop("'sd.multiplier' must be a single non-negative number", call. = FALSE)
    }

    if (!is.logical(plot.results) || length(plot.results) != 1) {
        stop("'plot.results' must be a single logical value", call. = FALSE)
    }

    bandwidth.indices <- seq_along(extrema.counts)

    # Fit initial linear model
    initial.model <- lm(extrema.counts ~ bandwidth.indices)

    # Fit segmented (piecewise) linear model
    tryCatch({
        piecewise.model <- segmented(initial.model, seg.Z = ~bandwidth.indices)
    }, error = function(e) {
        stop("Failed to fit piecewise linear model: ", e$message, call. = FALSE)
    })

    # Extract breakpoint (elbow point)
    breakpoint <- piecewise.model$psi[1, 2]

    # Extract fitted values
    fitted.extrema.counts <- predict(piecewise.model)

    # Calculate residuals and their standard deviation
    model.residuals <- extrema.counts - fitted.extrema.counts
    residual.sd <- sd(model.residuals)

    # Calculate optimal bandwidth index
    optimal.bandwidth.idx <- min(
        round(breakpoint + sd.multiplier * residual.sd),
        length(bandwidth.indices)
    )

    # Ensure index is at least 1
    optimal.bandwidth.idx <- max(optimal.bandwidth.idx, 1)

    # Optionally plot the results
    if (plot.results) {
        # Save current par settings
        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))

        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

        plot(bandwidth.indices, extrema.counts,
             type = "o", pch = 16, cex = 0.8,
             xlab = "Bandwidth index",
             ylab = "Number of local extrema",
             main = "Bandwidth Selection via Local Extrema Elbow Method")

        lines(bandwidth.indices, fitted.extrema.counts,
              col = "red", lwd = 2)
        abline(v = breakpoint, col = "blue", lty = 2)
        abline(v = optimal.bandwidth.idx, col = "purple", lty = 2, lwd = 2)

        legend("topright",
               legend = c("Observed",
                          "Fitted piecewise model",
                          "Breakpoint",
                          paste0("Optimal (BP + ", sd.multiplier, " SD)")),
               col = c("black", "red", "blue", "purple"),
               lty = c(1, 1, 2, 2),
               pch = c(16, NA, NA, NA),
               lwd = c(1, 2, 1, 2),
               cex = 0.9,
               bg = "white",
               box.lty = 1)
    }

    return(optimal.bandwidth.idx)
}

#' Print Method for graph_kernel_smoother Objects
#'
#' @description
#' Prints a concise summary of a \code{graph_kernel_smoother} object.
#'
#' @param x A \code{graph_kernel_smoother} object.
#' @param digits Number of significant digits to display. Default is 4.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
#' @method print graph_kernel_smoother
print.graph_kernel_smoother <- function(x, digits = 4, ...) {
    cat("Graph-Based Kernel Smoother with Buffer Zone Cross-Validation\n")
    cat("==============================================================\n")
    cat("Number of vertices:", length(x$predictions), "\n")
    cat("Optimal bandwidth:", format(x$opt_bw, digits = digits), "\n")
    cat("Buffer hops used:", x$buffer_hops_used, "\n")
    cat("Bandwidths tested:", length(x$bws), "\n")
    cat("Min CV error:", format(min(x$bw_errors), digits = digits), "\n")
    cat("\nUse summary() for detailed statistics.\n")

    invisible(x)
}

#' Summary Method for graph_kernel_smoother Objects
#'
#' @description
#' Computes and displays detailed summary statistics for a
#' \code{graph_kernel_smoother} object.
#'
#' @param object A \code{graph_kernel_smoother} object.
#' @param digits Number of significant digits to display. Default is 4.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns a list containing summary statistics with components:
#'   \describe{
#'     \item{\code{n_vertices}}{Number of vertices in the graph}
#'     \item{\code{n_bandwidths}}{Number of bandwidths tested}
#'     \item{\code{opt_bw}}{Optimal bandwidth value}
#'     \item{\code{opt_bw_idx}}{Index of optimal bandwidth}
#'     \item{\code{buffer_hops_used}}{Buffer hop distance used}
#'     \item{\code{min_error}}{Minimum cross-validation error}
#'     \item{\code{bw_errors}}{Vector of all cross-validation errors}
#'     \item{\code{pred_mean}}{Mean of predictions}
#'     \item{\code{pred_sd}}{Standard deviation of predictions}
#'     \item{\code{pred_range}}{Range of predictions}
#'     \item{\code{pred_quantiles}}{Quantiles of predictions}
#'   }
#'
#' @export
#' @method summary graph_kernel_smoother
summary.graph_kernel_smoother <- function(object, digits = 4, ...) {
    result <- list()

    # Basic information
    result$n_vertices <- length(object$predictions)
    result$n_bandwidths <- length(object$bws)
    result$opt_bw <- object$opt_bw
    result$opt_bw_idx <- object$opt_bw_idx
    result$buffer_hops_used <- object$buffer_hops_used

    # Error statistics
    result$min_error <- min(object$bw_errors)
    result$bw_errors <- object$bw_errors

    # Prediction statistics
    result$pred_mean <- mean(object$predictions)
    result$pred_sd <- sd(object$predictions)
    result$pred_range <- range(object$predictions)
    result$pred_quantiles <- quantile(object$predictions,
                                     probs = c(0.25, 0.5, 0.75))

    # Display summary
    cat("Summary of Graph-Based Kernel Smoother\n")
    cat("======================================\n\n")

    cat("Data:\n")
    cat("  Number of vertices:", result$n_vertices, "\n")
    cat("  Buffer hop distance:", result$buffer_hops_used, "\n\n")

    cat("Bandwidth selection:\n")
    cat("  Optimal bandwidth:", format(result$opt_bw, digits = digits),
        "(index", result$opt_bw_idx, "of", result$n_bandwidths, ")\n")
    cat("  CV error at optimum:", format(result$min_error, digits = digits), "\n")
    cat("  CV error range: [", format(min(result$bw_errors), digits = digits),
        ", ", format(max(result$bw_errors), digits = digits), "]\n\n", sep = "")

    cat("Prediction summary:\n")
    cat("  Mean:", format(result$pred_mean, digits = digits), "\n")
    cat("  SD:", format(result$pred_sd, digits = digits), "\n")
    cat("  Min:", format(result$pred_range[1], digits = digits), "\n")
    cat("  1st Qu:", format(result$pred_quantiles[1], digits = digits), "\n")
    cat("  Median:", format(result$pred_quantiles[2], digits = digits), "\n")
    cat("  3rd Qu:", format(result$pred_quantiles[3], digits = digits), "\n")
    cat("  Max:", format(result$pred_range[2], digits = digits), "\n")

    invisible(result)
}

#' Plot Method for graph_kernel_smoother Objects
#'
#' @description
#' Creates diagnostic plots for \code{graph_kernel_smoother} objects.
#'
#' @param x A \code{graph_kernel_smoother} object.
#' @param which Character string or integer specifying the plot type:
#'   \itemize{
#'     \item \code{"cv"} or 1: Cross-validation error vs bandwidth
#'     \item \code{"predictions"} or 2: Predictions vs vertex index
#'     \item \code{"both"} or 3: Both plots in a 2x1 layout
#'   }
#' @param main Optional main title for the plot(s). If \code{NULL}, default
#'   titles are used.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' # See examples in graph.kernel.smoother()
#'
#' @importFrom graphics plot points abline legend par
#'
#' @export
#' @method plot graph_kernel_smoother
plot.graph_kernel_smoother <- function(
    x,
    which = "cv",
    main = NULL,
    ...
) {
    # Convert 'which' to consistent format
    which <- switch(as.character(which),
                   "1" = "cv",
                   "2" = "predictions",
                   "3" = "both",
                   which)

    if (!which %in% c("cv", "predictions", "both")) {
        stop("'which' must be 'cv', 'predictions', 'both', or 1, 2, 3",
             call. = FALSE)
    }

    # Save current par settings
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    if (which == "both") {
        par(mfrow = c(2, 1))
    }

    if (which %in% c("cv", "both")) {
        # Plot cross-validation errors vs bandwidths
        plot_main <- if (!is.null(main) && which == "cv") {
            main
        } else {
            "Bandwidth Selection via Buffer Zone CV"
        }

        plot(x$bws, x$bw_errors,
             type = "b", log = "x",
             xlab = "Bandwidth",
             ylab = "Cross-validation error",
             main = plot_main,
             pch = 16, cex = 0.8, ...)

        # Highlight optimal bandwidth
        points(x$opt_bw, x$bw_errors[x$opt_bw_idx],
               col = "red", pch = 19, cex = 1.5)
        abline(v = x$opt_bw, col = "red", lty = 2)

        # Add legend
        legend("topright",
               legend = c("CV errors", "Optimal bandwidth"),
               col = c("black", "red"),
               pch = c(16, 19),
               lty = c(1, 2),
               bg = "white",
               box.lty = 1)
    }

    if (which %in% c("predictions", "both")) {
        # Plot predictions vs vertex indices
        plot_main <- if (!is.null(main) && which == "predictions") {
            main
        } else {
            "Smoothed Predictions"
        }

        plot(seq_along(x$predictions), x$predictions,
             type = "p", pch = 16, cex = 0.6,
             xlab = "Vertex index",
             ylab = "Predicted value",
             main = plot_main, ...)

        # Add horizontal line at mean
        abline(h = mean(x$predictions), col = "blue", lty = 2)
    }

    invisible(x)
}

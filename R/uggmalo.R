#' Uniform Grid Graph Model-Averaged Local Linear Regression (UGGMALO)
#'
#' @description
#' \strong{Note: This is an experimental function and may produce errors under
#'     certain conditions. Use with caution and please report issues.}
#'
#' Implements the Uniform Grid Graph Model-Averaged Local linear regression
#' (UGGMALO) algorithm for estimating conditional expectations of functions
#' defined over vertices of a graph using local path linear models with
#' model averaging.
#'
#' The algorithm performs the following main steps:
#' \enumerate{
#'   \item Computes graph diameter and determines bandwidth range
#'   \item Creates a uniform grid representation of the input graph
#'   \item For each candidate bandwidth:
#'   \itemize{
#'     \item Processes paths through grid vertices
#'     \item Fits local logistic models to path data
#'     \item Computes weighted predictions and errors
#'   }
#'   \item Determines optimal bandwidth based on cross-validation errors
#' }
#'
#' The algorithm uses weighted logistic regression on paths through the graph
#' to create local models, which are then combined using weighted averaging.
#' Model evaluation is performed using leave-one-out cross-validation with
#' Brier score errors.
#'
#' @param adj.list A list of integer vectors representing the adjacency list
#'   of the graph. Each element i contains the indices of vertices adjacent
#'   to vertex i (1-based indexing).
#' @param weight.list A list of numeric vectors representing the weights of
#'   the edges. Must have the same structure as \code{adj.list}.
#' @param y A numeric vector of observations at each vertex. Length must
#'   match the number of vertices.
#' @param best.models.coverage.factor Numeric scalar between 0.5 and 1.0
#'   controlling model coverage. Default: 0.9
#' @param min.bw.factor Numeric scalar. Minimum bandwidth factor relative
#'   to graph diameter. Must be positive. Default: 0.05
#' @param max.bw.factor Numeric scalar. Maximum bandwidth factor relative
#'   to graph diameter. Must be greater than \code{min.bw.factor}.
#'   Default: 0.5
#' @param n.bws Integer. Number of bandwidths to test between
#'   \code{min.bw.factor} and \code{max.bw.factor}. Must be positive.
#'   Default: 50
#' @param grid.size Integer. Size of the evaluation grid for predictions.
#'   Must be positive. Default: 100
#' @param start.vertex Integer. Index of the starting vertex (1-based).
#'   Must be between 1 and the number of vertices. Default: 1
#' @param snap.tolerance Numeric scalar. Tolerance for snapping distances
#'   to grid points. Must be positive. Default: 0.1
#' @param dist.normalization.factor Numeric scalar. Factor for normalizing
#'   distances. Must be greater than 1. Default: 1.01
#' @param min.path.size Integer. Minimum path size for distance calculations.
#'   Must be at least 5. Default: 5
#' @param diff.threshold Integer. Threshold for difference in path lengths.
#'   Must be at least 5. Default: 5
#' @param kernel.type Integer between 0 and 7. Type of kernel to use:
#'   \itemize{
#'     \item 0: Uniform kernel
#'     \item 1: Triangular kernel
#'     \item 2: Epanechnikov kernel
#'     \item 3: Quartic kernel
#'     \item 4: Triweight kernel
#'     \item 5: Tricube kernel
#'     \item 6: Gaussian kernel
#'     \item 7: Cosine kernel
#'   }
#'   Default: 7
#' @param fit.quadratic Logical. Whether to fit quadratic terms in the
#'   local models. Default: FALSE
#' @param tolerance Numeric scalar. Convergence tolerance for optimization.
#'   Must be positive. Default: 1e-8
#' @param n.bb Integer. Number of bag bootstrap iterations. Must be
#'   non-negative. Default: 0
#' @param p Numeric scalar between 0 and 1. Probability parameter for
#'   bootstrap. Default: 0.95
#' @param n.perms Integer. Number of permutations for uncertainty
#'   estimation. Must be non-negative. Default: 0
#' @param verbose Logical. Whether to print progress messages.
#'   Default: FALSE
#'
#' @return A list of class "uggmalo" containing:
#' \describe{
#' \item{candidate_bws}{Numeric vector of candidate bandwidths tested}
#' \item{bw_predictions}{Matrix of predictions for each bandwidth (vertices x bandwidths)}
#' \item{mean_errors}{Numeric vector of mean cross-validation errors
#'   for each bandwidth}
#' \item{opt_bw_idx}{Integer. Index of the optimal bandwidth (1-based)}
#' \item{predictions}{Numeric vector of predictions using the optimal
#'   bandwidth}
#' \item{graph_diameter}{Numeric. Diameter of the graph}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Create a simple chain graph
#' set.seed(123)  # For reproducibility
#' n.vertices <- 20
#' adj.list <- lapply(1:n.vertices, function(i) {
#'   # Simple chain graph
#'   if (i == 1) c(2L)
#'   else if (i == n.vertices) c(i - 1L)
#'   else c(i - 1L, i + 1L)
#' })
#'
#' dist.list <- lapply(1:n.vertices, function(i) {
#'   if (i == 1) c(1.0)
#'   else if (i == n.vertices) c(1.0)
#'   else c(1.0, 1.0)
#' })
#'
#' # Create parabolic signal with noise
#' x <- seq(-1, 1, length.out = n.vertices)
#' true.signal <- 3 * x^2 - 2 * x + 1  # Parabola
#' noise <- rnorm(n.vertices, mean = 0, sd = 0.3)
#' y <- true.signal + noise
#'
#' # Run estimation with default parameters
#' result <- uggmalo(adj.list, dist.list, y)
#'
#' # Compare true vs estimated values
#' plot(1:n.vertices, y, pch = 19, col = "gray50",
#'      xlab = "Vertex", ylab = "Value",
#'      main = "UGGMALO: True vs Estimated")
#' lines(1:n.vertices, true.signal, col = "blue", lwd = 2)
#' lines(1:n.vertices, result$predictions, col = "red", lwd = 2)
#' legend("topright", c("Observations", "True signal", "UGGMALO estimate"),
#'        col = c("gray50", "blue", "red"),
#'        pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = 2)
#'
#' # Run with custom parameters for comparison
#' result2 <- uggmalo(adj.list, dist.list, y,
#'                    min.bw.factor = 0.1,
#'                    max.bw.factor = 0.8,
#'                    n.bws = 30,
#'                    verbose = TRUE)
#' }
#'
#' @export
uggmalo <- function(adj.list,
                    weight.list,
                    y,
                    best.models.coverage.factor = 0.9,
                    min.bw.factor = 0.05,
                    max.bw.factor = 0.5,
                    n.bws = 50L,
                    grid.size = 100L,
                    start.vertex = 1L,
                    snap.tolerance = 0.1,
                    dist.normalization.factor = 1.01,
                    min.path.size = 5L,
                    diff.threshold = 5L,
                    kernel.type = 7L,
                    fit.quadratic = FALSE,
                    tolerance = 1e-8,
                    n.bb = 0L,
                    p = 0.95,
                    n.perms = 0L,
                    verbose = FALSE) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!all(vapply(adj.list, is.numeric, logical(1)))) {
        stop("All elements of 'adj.list' must be numeric vectors", call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    if (!all(vapply(weight.list, is.numeric, logical(1)))) {
        stop("All elements of 'weight.list' must be numeric vectors", call. = FALSE)
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
    }

    # Check matching lengths for each vertex
    lens_match <- vapply(seq_along(adj.list), function(i) {
        length(adj.list[[i]]) == length(weight.list[[i]])
    }, logical(1))

    if (!all(lens_match)) {
        stop("Corresponding elements in 'adj.list' and 'weight.list' must have the same length",
             call. = FALSE)
    }

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    if (length(y) != length(adj.list)) {
        stop("Length of 'y' must match the number of vertices", call. = FALSE)
    }

    # Validate numeric parameters
    if (!is.numeric(best.models.coverage.factor) || length(best.models.coverage.factor) != 1) {
        stop("'best.models.coverage.factor' must be a single numeric value", call. = FALSE)
    }

    if (best.models.coverage.factor <= 0.5 || best.models.coverage.factor > 1.0) {
        stop("'best.models.coverage.factor' must be greater than 0.5 and not greater than 1",
             call. = FALSE)
    }

    if (!is.numeric(min.bw.factor) || length(min.bw.factor) != 1 || min.bw.factor <= 0) {
        stop("'min.bw.factor' must be a single positive number", call. = FALSE)
    }

    if (!is.numeric(max.bw.factor) || length(max.bw.factor) != 1) {
        stop("'max.bw.factor' must be a single numeric value", call. = FALSE)
    }

    if (max.bw.factor <= min.bw.factor) {
        stop("'max.bw.factor' must be greater than 'min.bw.factor'", call. = FALSE)
    }

    # Validate integer parameters
    n.bws <- as.integer(n.bws)
    if (length(n.bws) != 1 || is.na(n.bws) || n.bws < 1) {
        stop("'n.bws' must be a positive integer", call. = FALSE)
    }

    grid.size <- as.integer(grid.size)
    if (length(grid.size) != 1 || is.na(grid.size) || grid.size < 1) {
        stop("'grid.size' must be a positive integer", call. = FALSE)
    }

    start.vertex <- as.integer(start.vertex)
    if (length(start.vertex) != 1 || is.na(start.vertex) ||
        start.vertex < 1 || start.vertex > length(adj.list)) {
        stop("'start.vertex' must be an integer between 1 and the number of vertices",
             call. = FALSE)
    }

    if (!is.numeric(snap.tolerance) || length(snap.tolerance) != 1 || snap.tolerance <= 0) {
        stop("'snap.tolerance' must be a single positive number", call. = FALSE)
    }

    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1 ||
        dist.normalization.factor <= 1) {
        stop("'dist.normalization.factor' must be a single number greater than 1",
             call. = FALSE)
    }

    min.path.size <- as.integer(min.path.size)
    if (length(min.path.size) != 1 || is.na(min.path.size) || min.path.size < 5) {
        stop("'min.path.size' must be an integer greater than or equal to 5",
             call. = FALSE)
    }

    diff.threshold <- as.integer(diff.threshold)
    if (length(diff.threshold) != 1 || is.na(diff.threshold) || diff.threshold < 5) {
        stop("'diff.threshold' must be an integer greater than or equal to 5",
             call. = FALSE)
    }

    kernel.type <- as.integer(kernel.type)
    if (length(kernel.type) != 1 || is.na(kernel.type) ||
        kernel.type < 0 || kernel.type > 7) {
        stop("'kernel.type' must be an integer between 0 and 7", call. = FALSE)
    }

    if (!is.logical(fit.quadratic) || length(fit.quadratic) != 1) {
        stop("'fit.quadratic' must be a single logical value", call. = FALSE)
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0) {
        stop("'tolerance' must be a single positive number", call. = FALSE)
    }

    n.bb <- as.integer(n.bb)
    if (length(n.bb) != 1 || is.na(n.bb) || n.bb < 0) {
        stop("'n.bb' must be a non-negative integer", call. = FALSE)
    }

    if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
        stop("'p' must be a single numeric value between 0 and 1", call. = FALSE)
    }

    n.perms <- as.integer(n.perms)
    if (length(n.perms) != 1 || is.na(n.perms) || n.perms < 0) {
        stop("'n.perms' must be a non-negative integer", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value", call. = FALSE)
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))
    start.vertex.0based <- as.integer(start.vertex - 1L)

    # Call the C++ function
    result <- .Call("S_uggmalo",
                   adj.list.0based,
                   weight.list,
                   as.numeric(y),
                   as.numeric(best.models.coverage.factor),
                   as.numeric(min.bw.factor),
                   as.numeric(max.bw.factor),
                   n.bws,
                   grid.size,
                   start.vertex.0based,
                   as.numeric(snap.tolerance),
                   as.numeric(dist.normalization.factor),
                   min.path.size,
                   diff.threshold,
                   kernel.type,
                   fit.quadratic,
                   as.numeric(tolerance),
                   n.bb,
                   as.numeric(p),
                   n.perms,
                   verbose)

    # Add class for S3 methods
    class(result) <- c("uggmalo", class(result))

    return(result)
}

#' Bayesian Bootstrap with Uncertainty for UGGMALO Predictions
#'
#' @description
#' \strong{Note: This is an experimental function and may produce errors under certain conditions.
#' Use with caution and please report issues.}
#'
#' Combines Bayesian bootstrap with permutation-based uncertainty estimation
#' for UGGMALO predictions on a graph. This function accounts for both the
#' uncertainty in p-value estimation from permutation tests and the spatial
#' structure of the data.
#'
#' The function performs these steps:
#' \enumerate{
#'   \item Generates Dirichlet weights for Bayesian bootstrap
#'   \item Resamples from permutation distributions to account for p-value
#'     uncertainty
#'   \item Applies spatial smoothing using UGGMALO to account for graph
#'     structure
#'   \item Computes credible intervals from the bootstrap distribution
#' }
#'
#' @param perm.results Numeric matrix where each row corresponds to a vertex
#'   and each column contains predictions from a permutation run. Dimensions
#'   should be (number of vertices x number of permutations).
#' @param true.predictions Numeric vector of original UGGMALO predictions
#'   for each vertex. Length must match the number of rows in
#'   \code{perm.results}.
#' @param graph List containing graph structure with required components:
#'   \describe{
#'     \item{pruned_adj_list}{Adjacency list representation of the graph
#'       (list of integer vectors)}
#'     \item{pruned_dist_list}{List of numeric vectors containing edge
#'       weights/distances}
#'   }
#' @param n.bootstrap Integer. Number of bootstrap iterations. Must be
#'   positive. Default: 1000
#' @param n.cores Integer. Number of CPU cores to use for parallel
#'   processing. Default: 14
#'
#' @return A list of class "uggmalo_bootstrap" containing:
#' \item{ci.lower}{Numeric vector of lower bounds of \eqn{95\%} credible intervals}
#' \item{ci.upper}{Numeric vector of upper bounds of \eqn{95\%} credible intervals}
#' \item{bootstrap.distribution}{Matrix of bootstrap estimates where each
#'   column represents one bootstrap iteration (vertices x iterations)}
#'
#' @seealso \code{\link{uggmalo}}
#'
#' @section Known Issues:
#' This function currently fails with the error "Reference vertex not found in path vertices"
#' (line 1944 in centered_paths.cpp). This is a known bug that will be addressed in future releases.
#'
#' @examples
#' \dontrun{
#' # WARNING: This is an experimental function with known issues
#' # The following example demonstrates the intended usage but currently fails
#' # with error: "Reference vertex not found in path vertices"
#'
#' if (requireNamespace("foreach", quietly = TRUE) &&
#'     requireNamespace("doParallel", quietly = TRUE)) {
#'
#'   # Create example data
#'   n.vertices <- 20
#'   n.perms <- 100
#'
#'   # Simulated permutation results
#'   perm.results <- matrix(rnorm(n.vertices * n.perms),
#'                          nrow = n.vertices, ncol = n.perms)
#'
#'   # True predictions
#'   true.predictions <- rnorm(n.vertices)
#'
#'   # Simple graph structure
#'   graph <- list(
#'     pruned_adj_list = lapply(1:n.vertices, function(i) {
#'       # Simple chain graph
#'       if (i == 1) c(2L)
#'       else if (i == n.vertices) c(i - 1L)
#'       else c(i - 1L, i + 1L)
#'     }),
#'     pruned_dist_list = lapply(1:n.vertices, function(i) {
#'       if (i == 1) c(1.0)
#'       else if (i == n.vertices) c(1.0)
#'       else c(1.0, 1.0)
#'     })
#'   )
#'
#'   # Run bootstrap (with fewer iterations for example)
#'   # NOTE: This currently fails with a known bug
#'   results <- uggmalo.bayesian.bootstrap.with.uncertainty(
#'     perm.results = perm.results,
#'     true.predictions = true.predictions,
#'     graph = graph,
#'     n.bootstrap = 100,
#'     n.cores = 2
#'   )
#'
#'   # Plot credible intervals (if successful)
#'   plot(1:n.vertices, true.predictions, pch = 19,
#'        ylim = range(c(results$ci.lower, results$ci.upper)),
#'        xlab = "Vertex", ylab = "Prediction",
#'        main = "UGGMALO Predictions with 95% Credible Intervals")
#'   arrows(1:n.vertices, results$ci.lower,
#'          1:n.vertices, results$ci.upper,
#'          length = 0.05, angle = 90, code = 3)
#' }
#' }
#'
#' @importFrom stats median quantile
#' @export
uggmalo.bayesian.bootstrap.with.uncertainty <- function(
    perm.results,
    true.predictions,
    graph,
    n.bootstrap = 1000L,
    n.cores = 14L
    ) {

    if (!isTRUE(getOption("gflow.suppress.experimental.warnings", FALSE))) {
        warning("This function is experimental and may not work reliably in all cases. ",
                "Set options(gflow.suppress.experimental.warnings = TRUE) to suppress this warning.",
                call. = FALSE)
     }

    # Check if required packages are available
    if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("Package 'foreach' is required. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("Package 'doParallel' is required. Please install it.", call. = FALSE)
    }

    # Input validation
    if (!is.matrix(perm.results) || !is.numeric(perm.results)) {
        stop("'perm.results' must be a numeric matrix", call. = FALSE)
    }

    if (!is.numeric(true.predictions) || !is.vector(true.predictions)) {
        stop("'true.predictions' must be a numeric vector", call. = FALSE)
    }

    if (nrow(perm.results) != length(true.predictions)) {
        stop("Number of rows in 'perm.results' must match length of 'true.predictions'",
             call. = FALSE)
    }

    if (!is.list(graph)) {
        stop("'graph' must be a list", call. = FALSE)
    }

    if (!all(c("pruned_adj_list", "pruned_dist_list") %in% names(graph))) {
        stop("'graph' must contain 'pruned_adj_list' and 'pruned_dist_list' components",
             call. = FALSE)
    }

    n.bootstrap <- as.integer(n.bootstrap)
    if (length(n.bootstrap) != 1 || is.na(n.bootstrap) || n.bootstrap < 1) {
        stop("'n.bootstrap' must be a positive integer", call. = FALSE)
    }

    n.cores <- as.integer(n.cores)
    if (length(n.cores) != 1 || is.na(n.cores) || n.cores < 1) {
        stop("'n.cores' must be a positive integer", call. = FALSE)
    }

    # Get dimensions
    n.vertices <- nrow(perm.results)
    n.perms <- ncol(perm.results)

    # Set up parallel processing
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (isTRUE(getOption("uggmalo.verbose", FALSE))) {
        message("Starting Bayesian bootstrap with ", n.bootstrap,
                " iterations on ", n.cores, " cores...")
        ptm <- proc.time()
    }

    # Export required objects to cluster
    parallel::clusterExport(cl, c("perm.results", "true.predictions", "graph",
                                  "n.vertices", "n.perms", "uggmalo"),
                            envir = environment())

    # Perform bootstrap iterations in parallel
    bootstrap.estimates <- foreach::foreach(
        b = seq_len(n.bootstrap),
        .combine = 'cbind',
        .packages = c("stats"),
        .errorhandling = 'stop'
    ) %dopar% {

        # Generate Bayesian bootstrap weights using package function
        weights <- runif.simplex(n.vertices)

        # For each vertex, resample from permutation results
        resampled.p.values <- vapply(seq_len(n.vertices), function(i) {
            # Resample from permutation distribution
            null.dist <- sample(perm.results[i, ], size = n.perms, replace = TRUE)

            # Calculate p-value using resampled distribution
            null.center <- stats::median(null.dist)
            observed.deviation <- abs(true.predictions[i] - null.center)
            p.value <- mean(abs(null.dist - null.center) >= observed.deviation)

            # Avoid log(0)
            p.value <- max(p.value, .Machine$double.eps)

            # Apply bootstrap weight
            -log(p.value) * weights[i]
        }, numeric(1))

        # Fit spatial model to resampled, weighted values
        pvals.res <- uggmalo(
            adj.list = graph$pruned_adj_list,
            weight.list = graph$pruned_dist_list,
            y = resampled.p.values,
            min.bw.factor = 0.1,
            max.bw.factor = 0.99,
            n.bws = 50L,
            grid.size = 100L,
            start.vertex = 1L,
            snap.tolerance = 0.1,
            min.path.size = 5L,
            diff.threshold = 5L,
            verbose = FALSE
        )

        # Return fitted values
        pvals.res$predictions
    }

    if (isTRUE(getOption("uggmalo.verbose", FALSE))) {
        elapsed <- proc.time() - ptm
        message("Bootstrap completed in ", round(elapsed[3], 2), " seconds")
    }

    # Calculate credible intervals from bootstrap distribution
    ci.lower <- apply(bootstrap.estimates, 1, stats::quantile,
                      probs = 0.025, na.rm = TRUE)
    ci.upper <- apply(bootstrap.estimates, 1, stats::quantile,
                      probs = 0.975, na.rm = TRUE)

    # Create result object
    result <- list(
        ci.lower = ci.lower,
        ci.upper = ci.upper,
        bootstrap.distribution = bootstrap.estimates
    )

    class(result) <- c("uggmalo_bootstrap", class(result))

    return(result)
}

#' Print Method for UGGMALO Results
#'
#' @param x An object of class "uggmalo"
#' @param digits Integer. Number of significant digits to print
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
#' @method print uggmalo
print.uggmalo <- function(x, digits = 4, ...) {
    cat("UGGMALO Results\n")
    cat("===============\n\n")

    cat("Graph diameter:", round(x$graph_diameter, digits), "\n")
    cat("Number of vertices:", length(x$predictions), "\n")
    cat("Optimal bandwidth index:", x$opt_bw_idx, "\n")
    cat("Optimal bandwidth value:",
        round(x$candidate_bws[x$opt_bw_idx], digits), "\n\n")

    cat("Cross-validation error summary:\n")
    err_summary <- summary(x$mean_errors)
    print(round(err_summary, digits))

    cat("\nPrediction summary:\n")
    pred_summary <- summary(x$predictions)
    print(round(pred_summary, digits))

    invisible(x)
}


#' Summary Method for UGGMALO Results
#'
#' @param object An object of class "uggmalo"
#' @param ... Additional arguments (currently ignored)
#'
#' @return A list of class "summary.uggmalo" containing summary statistics
#'
#' @export
#' @method summary uggmalo
summary.uggmalo <- function(object, ...) {
    result <- list(
        n.vertices = length(object$predictions),
        graph.diameter = object$graph_diameter,
        n.bandwidths = length(object$candidate_bws),
        bandwidth.range = range(object$candidate_bws),
        optimal.bandwidth = object$candidate_bws[object$opt_bw_idx],
        optimal.bw.index = object$opt_bw_idx,
        cv.error.range = range(object$mean_errors),
        cv.error.optimal = object$mean_errors[object$opt_bw_idx],
        prediction.summary = summary(object$predictions)
    )

    class(result) <- "summary.uggmalo"
    return(result)
}


#' Print Method for UGGMALO Summary
#'
#' @param x An object of class "summary.uggmalo"
#' @param digits Integer. Number of significant digits to print
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
#' @method print summary.uggmalo
print.summary.uggmalo <- function(x, digits = 4, ...) {
    cat("UGGMALO Model Summary\n")
    cat("====================\n\n")

    cat("Graph information:\n")
    cat("  Number of vertices:", x$n.vertices, "\n")
    cat("  Graph diameter:", round(x$graph.diameter, digits), "\n\n")

    cat("Bandwidth selection:\n")
    cat("  Number tested:", x$n.bandwidths, "\n")
    cat("  Range: [", round(x$bandwidth.range[1], digits), ", ",
        round(x$bandwidth.range[2], digits), "]\n", sep = "")
    cat("  Optimal:", round(x$optimal.bandwidth, digits),
        " (index ", x$optimal.bw.index, ")\n\n", sep = "")

    cat("Cross-validation errors:\n")
    cat("  Range: [", round(x$cv.error.range[1], digits), ", ",
        round(x$cv.error.range[2], digits), "]\n", sep = "")
    cat("  At optimal bandwidth:", round(x$cv.error.optimal, digits), "\n\n")

    cat("Predictions:\n")
    print(round(x$prediction.summary, digits))

    invisible(x)
}

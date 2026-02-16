#' Compute the Kernel-weighted Graph Neighbor Mean
#'
#' @description
#' Calculates the kernel-weighted mean value of the neighbors for each vertex
#' in a graph. The graph is represented as an adjacency list, with distances
#' provided as edge lengths, and values associated with each vertex.
#'
#' @param adj.list A list of integer vectors representing the adjacency list of
#'   the graph. Each element contains the indices of neighboring vertices
#'   (1-based indexing).
#' @param edge.lengths A list of numeric vectors containing the edge lengths
#'   (distances) corresponding to the neighbors in \code{adj.list}. Must have
#'   the same structure as \code{adj.list}.
#' @param y A numeric vector of values associated with each vertex in the graph.
#'   Length must equal the number of vertices.
#' @param kernel An integer specifying the kernel function to use:
#'   \itemize{
#'     \item 1 = Epanechnikov kernel
#'     \item 2 = Triangular kernel
#'     \item 3 = Truncated exponential kernel
#'     \item 4 = Normal (Gaussian) kernel
#'   }
#'   Default is 1.
#' @param dist.normalization.factor A positive numeric value (>= 1) used to
#'   normalize distances in the kernel computation. Default is 1.01.
#'
#' @return A numeric vector containing the kernel-weighted neighbor mean for
#'   each vertex.
#'
#' @details
#' The function computes weighted averages where weights are determined by a
#' kernel function applied to the edge distances. Smaller distances result in
#' larger weights. The distance normalization factor controls the bandwidth of
#' the kernel.
#'
#' @examples
#' \dontrun{
#' # Simple triangle graph
#' graph <- list(c(2, 3), c(1, 3), c(1, 2))
#' edge.lengths <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
#' y <- c(1.0, 2.0, 3.0)
#'
#' # Compute kernel-weighted means with Epanechnikov kernel
#' y.kwmean <- graph.kmean(graph, edge.lengths, y, kernel = 1)
#' print(y.kwmean)
#' }
#'
#' @seealso \code{\link{graph.kmean.cv}} for cross-validation
#' @export
graph.kmean <- function(adj.list,
                        edge.lengths,
                        y,
                        kernel = 1,
                        dist.normalization.factor = 1.01) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }
    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list")
    }

    # Check lengths
    n.vertices <- length(y)
    if (length(adj.list) != n.vertices) {
        stop("Length of 'adj.list' (", length(adj.list),
             ") does not match length of 'y' (", n.vertices, ")")
    }
    if (length(edge.lengths) != n.vertices) {
        stop("Length of 'edge.lengths' (", length(edge.lengths),
             ") does not match length of 'y' (", n.vertices, ")")
    }

    # Check structure consistency
    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) != length(edge.lengths[[i]])) {
            stop("Structure mismatch at vertex ", i,
                 ": adj.list has ", length(adj.list[[i]]),
                 " neighbors but edge.lengths has ", length(edge.lengths[[i]]))
        }
    }

    # Validate parameters
    if (!is.numeric(y)) {
        stop("'y' must be numeric")
    }
    if (any(is.na(y))) {
        stop("'y' cannot contain NA values")
    }
    if (!kernel %in% 1:4) {
        stop("'kernel' must be an integer between 1 and 4")
    }
    if (!is.numeric(dist.normalization.factor) ||
        length(dist.normalization.factor) != 1 ||
        dist.normalization.factor < 1) {
        stop("'dist.normalization.factor' must be a single numeric value >= 1")
    }

    # Convert to 0-based indexing for C++ interface
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ implementation
    .Call("S_graph_kmean",
          adj.list.0based,
          edge.lengths,
          as.numeric(y),
          as.integer(kernel),
          as.numeric(dist.normalization.factor))
}

#' Cross-validation for Kernel-weighted Graph Neighbor Mean
#'
#' @description
#' Performs k-fold cross-validation to evaluate the performance of the
#' kernel-weighted neighbor mean algorithm on a graph. This function helps in
#' selecting optimal parameters such as the kernel type and distance
#' normalization factor.
#'
#' @param graph A list of integer vectors representing the adjacency list
#'   (1-based indexing).
#' @param edge.lengths A list of numeric vectors containing edge lengths.
#'   Structure must match \code{graph}.
#' @param y A numeric vector of response values at each vertex.
#' @param kernel Either an integer (1-4) or a character string specifying the
#'   kernel function:
#'   \itemize{
#'     \item 1 or "epanechnikov" = Epanechnikov kernel (default)
#'     \item 2 or "triangular" = Triangular kernel
#'     \item 3 or "truncated_exponential" = Truncated exponential kernel
#'     \item 4 or "normal" = Normal (Gaussian) kernel
#'   }
#' @param dist.normalization.factor A positive numeric value (>= 1) for
#'   distance normalization. Default is 1.01.
#' @param n.CVs Number of cross-validation iterations. Default is 10.
#' @param n.CV.folds Number of folds for each CV iteration. Must be >= 2.
#'   Default is 10.
#' @param seed Integer seed for reproducibility. Default is 0.
#' @param use.weighted.MAD.error Logical; if TRUE, uses weighted Mean Absolute
#'   Deviation for binary classification problems to handle class imbalance.
#'   Default is FALSE.
#'
#' @return A numeric vector of cross-validation errors for each vertex. Values
#'   may be NaN for vertices that were excluded in all CV iterations.
#'
#' @details
#' The function performs repeated k-fold cross-validation. In each iteration,
#' the data is randomly partitioned into k folds. Each fold is held out once
#' while the model is trained on the remaining folds, and predictions are made
#' for the held-out vertices.
#'
#' When \code{use.weighted.MAD.error = TRUE}, the function applies class
#' weights to handle imbalanced binary classification:
#' \itemize{
#'   \item Weight for class 0: 1 / (1 - q)
#'   \item Weight for class 1: 1 / q
#' }
#' where q is the proportion of class 1 samples.
#'
#' @examples
#' # Triangle graph example
#' graph <- list(c(2, 3), c(1, 3), c(1, 2))
#' edge.lengths <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
#' y <- c(1.0, 2.0, 3.0)
#'
#' # Perform cross-validation
#' cv.errors <- graph.kmean.cv(graph, edge.lengths, y,
#'                             kernel = "epanechnikov",
#'                             n.CVs = 5, n.CV.folds = 3)
#' print(cv.errors)
#'
#' @seealso \code{\link{graph.kmean}} for the main function
#' @export
graph.kmean.cv <- function(graph,
                          edge.lengths,
                          y,
                          kernel = 1,
                          dist.normalization.factor = 1.01,
                          n.CVs = 10,
                          n.CV.folds = 10,
                          seed = 0,
                          use.weighted.MAD.error = FALSE) {

    # Input validation
    if (!is.list(graph)) {
        stop("'graph' must be a list")
    }
    if (!is.list(edge.lengths)) {
        stop("'edge.lengths' must be a list")
    }

    # Check dimensions
    n.vertices <- length(y)
    if (length(graph) != n.vertices) {
        stop("Length of 'graph' (", length(graph),
             ") does not match length of 'y' (", n.vertices, ")")
    }
    if (length(edge.lengths) != n.vertices) {
        stop("Length of 'edge.lengths' (", length(edge.lengths),
             ") does not match length of 'y' (", n.vertices, ")")
    }

    # Check structure consistency
    for (i in seq_along(graph)) {
        if (length(graph[[i]]) != length(edge.lengths[[i]])) {
            stop("Structure mismatch at vertex ", i)
        }
    }

    # Handle kernel parameter
    if (is.character(kernel)) {
        kernel <- match.arg(kernel,
                           c("epanechnikov", "triangular",
                             "truncated_exponential", "normal"))
        kernel <- switch(kernel,
                        epanechnikov = 1,
                        triangular = 2,
                        truncated_exponential = 3,
                        normal = 4)
    } else if (!is.numeric(kernel) || !kernel %in% 1:4) {
        stop("'kernel' must be 1-4 or a valid kernel name")
    }

    # Validate other parameters
    if (!is.numeric(dist.normalization.factor) ||
        length(dist.normalization.factor) != 1 ||
        dist.normalization.factor < 1) {
        stop("'dist.normalization.factor' must be a single numeric value >= 1")
    }

    if (!is.numeric(n.CVs) || length(n.CVs) != 1 || n.CVs < 1) {
        stop("'n.CVs' must be a positive integer")
    }
    n.CVs <- as.integer(n.CVs)

    if (!is.numeric(n.CV.folds) || length(n.CV.folds) != 1 || n.CV.folds < 2) {
        stop("'n.CV.folds' must be an integer >= 2")
    }
    n.CV.folds <- as.integer(n.CV.folds)

    if (n.CV.folds > n.vertices) {
        warning("'n.CV.folds' exceeds number of vertices; setting to ",
                n.vertices)
        n.CV.folds <- n.vertices
    }

    if (!is.numeric(seed) || length(seed) != 1 || seed < 0) {
        stop("'seed' must be a non-negative integer")
    }
    seed <- as.integer(seed)

    if (!is.logical(use.weighted.MAD.error) ||
        length(use.weighted.MAD.error) != 1) {
        stop("'use.weighted.MAD.error' must be TRUE or FALSE")
    }

    # Convert to 0-based indexing
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    # Call appropriate C++ function
    if (use.weighted.MAD.error) {
        .Call("S_graph_kmean_wmad_cv",
              graph.0based,
              edge.lengths,
              as.numeric(y),
              as.integer(kernel),
              as.numeric(dist.normalization.factor),
              as.integer(n.CVs),
              as.integer(n.CV.folds),
              as.integer(seed),
              as.logical(use.weighted.MAD.error))
    } else {
        .Call("S_graph_kmean_cv",
              graph.0based,
              edge.lengths,
              as.numeric(y),
              as.integer(kernel),
              as.numeric(dist.normalization.factor),
              as.integer(n.CVs),
              as.integer(n.CV.folds),
              as.integer(seed))
    }
}


#' Adaptive Neighborhood Size Graph K-Means for Univariate Data
#'
#' @description
#' Performs nonparametric regression using graph-based k-means with adaptive
#' neighborhood size selection. The method constructs chain graphs from sorted
#' predictor values and uses cross-validation to select the optimal
#' neighborhood size.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param y.true Optional numeric vector of true response values for performance
#'   evaluation.
#' @param use.median Logical; if TRUE uses median, if FALSE uses mean for
#'   predictions. Default is FALSE.
#' @param h.min Minimum neighborhood size to consider (>= 2). Default is 2.
#' @param h.max Maximum neighborhood size to consider. Default is
#'   min(30, length(x)).
#' @param n.CVs Number of cross-validation iterations. Default is 1000.
#' @param n.CV.folds Number of CV folds (>= 2). Default is 10.
#' @param p Probability level for credible intervals (0 < p < 1). Default is
#'   0.95.
#' @param n.bb Number of Bayesian bootstrap iterations (0 to skip). Default is
#'   500.
#' @param ikernel Integer specifying kernel function (1-6). Default is 1.
#' @param n.cores Number of CPU cores for parallel processing. Default is 1.
#' @param dist.normalization.factor Distance normalization factor (> 1).
#'   Default is 1.01.
#' @param epsilon Small positive value for numerical stability. Default is
#'   1e-15.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return An object of class "ugkmm" containing:
#'   \item{h_values}{Vector of tested neighborhood sizes}
#'   \item{h_cv_errors}{Cross-validation errors for each h}
#'   \item{mean_cv_error}{Mean CV error across all h values}
#'   \item{opt_h}{Optimal neighborhood size}
#'   \item{predictions}{Fitted values}
#'   \item{bb_predictions}{Bootstrap central tendency estimates}
#'   \item{opt_ci_lower}{Lower credible interval bounds}
#'   \item{opt_ci_upper}{Upper credible interval bounds}
#'   \item{x_sorted}{Sorted predictor values}
#'   \item{y_sorted}{Response values (sorted by x)}
#'   \item{y_true_sorted}{True values (sorted by x) if provided}
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Sorts data by predictor values
#'   \item Constructs chain graphs with varying neighborhood sizes
#'   \item Uses cross-validation to select optimal neighborhood size
#'   \item Computes predictions using kernel-weighted means
#'   \item Optionally performs Bayesian bootstrap for uncertainty quantification
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#'
#' # Fit model
#' result <- univariate.gkmm(x, y, h.min = 2, h.max = 20, n.bb = 100)
#'
#' # Plot results
#' plot(result)
#'
#' # Summary
#' summary(result)
#' }
#'
#' @seealso
#' \code{\link{plot.ugkmm}} for plotting methods
#' \code{\link{summary.ugkmm}} for summary statistics
#'
#' @export
univariate.gkmm <- function(x, y, y.true = NULL,
                           use.median = FALSE,
                           h.min = 2,
                           h.max = min(30, length(x)),
                           n.CVs = 1000,
                           n.CV.folds = 10,
                           p = 0.95,
                           n.bb = 500,
                           ikernel = 1,
                           n.cores = 1,
                           dist.normalization.factor = 1.01,
                           epsilon = 1e-15,
                           seed = NULL) {

    # Store the call
    cl <- match.call()

    # Input validation
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("'x' and 'y' must be numeric vectors")
    }
    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }
    if (length(x) < 3) {
        stop("Need at least 3 observations")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop("'x' and 'y' cannot contain NA values")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("'x' and 'y' must contain finite values only")
    }

    # Sort data by x
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (!is.null(y.true)) {
        y.true <- y.true[ord]
    }

    # Validate y.true
    if (!is.null(y.true)) {
        if (!is.numeric(y.true) || length(y.true) != length(x)) {
            stop("'y.true' must be numeric with same length as 'x'")
        }
        if (any(is.na(y.true))) {
            stop("'y.true' cannot contain NA values")
        }
    }

    # Validate other parameters
    if (!is.logical(use.median) || length(use.median) != 1) {
        stop("'use.median' must be TRUE or FALSE")
    }

    h.min <- as.integer(h.min)
    h.max <- as.integer(h.max)
    if (h.min < 2) {
        stop("'h.min' must be at least 2")
    }
    if (h.max <= h.min) {
        stop("'h.max' must be greater than 'h.min'")
    }
    if (h.max > length(x)) {
        h.max <- length(x)
        warning("'h.max' exceeds data size; setting to ", h.max)
    }

    n.CVs <- as.integer(n.CVs)
    if (n.CVs < 1) {
        stop("'n.CVs' must be positive")
    }

    n.CV.folds <- as.integer(n.CV.folds)
    if (n.CV.folds < 2 || n.CV.folds > length(x)) {
        stop("'n.CV.folds' must be between 2 and length(x)")
    }

    if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
        stop("'p' must be a single value between 0 and 1")
    }

    n.bb <- as.integer(n.bb)
    if (n.bb < 0) {
        stop("'n.bb' must be non-negative")
    }

    if (!ikernel %in% 1:6) {
        stop("'ikernel' must be between 1 and 6")
    }

    n.cores <- as.integer(n.cores)
    if (n.cores < 1) {
        stop("'n.cores' must be positive")
    }

    if (!is.numeric(dist.normalization.factor) ||
        dist.normalization.factor <= 1) {
        stop("'dist.normalization.factor' must be > 1")
    }

    if (!is.numeric(epsilon) || epsilon <= 0) {
        stop("'epsilon' must be positive")
    }

    # Handle seed
    if (!is.null(seed)) {
        seed <- as.integer(seed)
        if (seed < 0) {
            stop("'seed' must be non-negative")
        }
        set.seed(seed)
    } else {
        seed <- 0L
    }

    # Call C++ implementation
    result <- .Call("S_univariate_gkmm",
                   as.double(x),
                   as.double(y),
                   if (is.null(y.true)) double() else as.double(y.true),
                   as.logical(use.median),
                   as.integer(h.min),
                   as.integer(h.max),
                   as.integer(n.CVs),
                   as.integer(n.CV.folds),
                   as.double(p),
                   as.integer(n.bb),
                   as.integer(ikernel),
                   as.integer(n.cores),
                   as.double(dist.normalization.factor),
                   as.double(epsilon),
                   seed)
    
    # Add sorted data to result
    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- y.true
    result$call <- cl

    # Set class
    class(result) <- "ugkmm"

    return(result)
}

#' Plot ugkmm Objects
#'
#' @description
#' Visualizes results from adaptive neighborhood size graph k-means regression,
#' including fitted values, credible intervals, residuals, and diagnostics.
#'
#' @param x An object of class "ugkmm".
#' @param type Plot type: "fit" (default), "diagnostic", "residuals", or
#'   "residuals_hist".
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ... Additional graphical parameters.
#'
#' @details
#' Plot types:
#' \itemize{
#'   \item "fit": Shows fitted curve with optional credible intervals
#'   \item "diagnostic": Cross-validation error vs neighborhood size
#'   \item "residuals": Residuals vs fitted values
#'   \item "residuals_hist": Histogram of residuals
#' }
#'
#' @return NULL (invisibly)
#'
#' @examples
#' # See examples in univariate.gkmm
#'
#' @method plot ugkmm
#' @export
plot.ugkmm <- function(x, type = c("fit", "diagnostic", "residuals",
                                   "residuals_hist"),
                       main = "", xlab = "", ylab = "", ...) {

    type <- match.arg(type)

    switch(type,
        fit = {
            # Basic fit plot
            ylim <- range(c(x$y_sorted, x$predictions), na.rm = TRUE)
            if (!is.null(x$opt_ci_lower)) {
                ylim <- range(c(ylim, x$opt_ci_lower, x$opt_ci_upper),
                             na.rm = TRUE)
            }

            plot(x$x_sorted, x$y_sorted,
                 main = main, xlab = xlab, ylab = ylab,
                 ylim = ylim, las = 1, ...)

            # Add credible intervals
            if (!is.null(x$opt_ci_lower)) {
                polygon(c(x$x_sorted, rev(x$x_sorted)),
                       c(x$opt_ci_lower, rev(x$opt_ci_upper)),
                       col = "gray90", border = NA)
            }

            # Add predictions
            lines(x$x_sorted, x$predictions, col = "blue", lwd = 2)

            # Add true values if available
            if (!is.null(x$y_true_sorted)) {
                lines(x$x_sorted, x$y_true_sorted, col = "red", lwd = 2, lty = 2)
                legend("topright", c("Fitted", "True"),
                       col = c("blue", "red"), lty = c(1, 2), lwd = 2)
            }
        },

        diagnostic = {
            if (!is.null(x$h_values) && !is.null(x$h_cv_errors)) {
                plot(x$h_values, x$h_cv_errors, type = 'b',
                     main = if (main == "") "Cross-validation Error" else main,
                     xlab = if (xlab == "") "Neighborhood size (h)" else xlab,
                     ylab = if (ylab == "") "CV Error" else ylab,
                     las = 1, ...)
                abline(v = x$opt_h, col = "red", lty = 2)
                text(x$opt_h, par("usr")[4],
                     paste("h =", x$opt_h), pos = 3, col = "red")
            }
        },

        residuals = {
            residuals <- x$y_sorted - x$predictions
            plot(x$x_sorted, residuals,
                 main = if (main == "") "Residual Plot" else main,
                 xlab = if (xlab == "") "x" else xlab,
                 ylab = if (ylab == "") "Residuals" else ylab,
                 las = 1, ...)
            abline(h = 0, col = "gray", lty = 2)

            # Add loess smoother
            lo <- try(loess(residuals ~ x$x_sorted), silent = TRUE)
            if (!inherits(lo, "try-error")) {
                lines(x$x_sorted, fitted(lo), col = "red", lwd = 2)
            }
        },

        residuals_hist = {
            residuals <- x$y_sorted - x$predictions
            hist(residuals,
                 main = if (main == "") "Residual Distribution" else main,
                 xlab = if (xlab == "") "Residuals" else xlab,
                 col = "lightblue", border = "white", las = 1, ...)

            # Add normal curve
            graphics::curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)) *
                  length(residuals) * diff(hist(residuals, plot = FALSE)$breaks)[1],
                  add = TRUE, col = "red", lwd = 2)
        }
    )

    invisible(NULL)
}

#' Summary Method for ugkmm Objects
#'
#' @description
#' Provides comprehensive summary statistics for adaptive neighborhood size
#' graph k-means regression results.
#'
#' @param object An object of class "ugkmm".
#' @param digits Number of significant digits for output. Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "summary.ugkmm" containing model statistics,
#'   cross-validation results, and fit diagnostics.
#'
#' @method summary ugkmm
#' @export
summary.ugkmm <- function(object, digits = 4, ...) {

    # Validate input
    if (!inherits(object, "ugkmm")) {
        stop("Input must be a 'ugkmm' object")
    }

    # Calculate residuals
    if (is.null(object$predictions) || is.null(object$y_sorted)) {
        stop("Object is missing required fields 'predictions' or 'y_sorted'")
    }
    residuals <- object$y_sorted - object$predictions

    # Model information - with safe access to all fields
    model_info <- list(
        n = length(object$x_sorted),
        optimal_h = if (!is.null(object$opt_h)) object$opt_h else NA,
        x_range = if (!is.null(object$x_sorted)) range(object$x_sorted) else c(NA, NA),
        method = if (!is.null(object$use.median)) {
            if (object$use.median) "Median" else "Mean"
        } else "Mean"  # Default to "Mean" if use.median is missing
    )

    # Fit statistics
    fit_stats <- list(
        mse = mean(residuals^2),
        rmse = sqrt(mean(residuals^2)),
        mae = mean(abs(residuals)),
        median_ae = median(abs(residuals))
    )

    # Calculate R-squared if we have y values
    if (!is.null(object$y_sorted) && length(object$y_sorted) > 0) {
        ss_tot <- sum((object$y_sorted - mean(object$y_sorted))^2)
        if (ss_tot > 0) {
            fit_stats$r_squared <- 1 - sum(residuals^2) / ss_tot
        } else {
            fit_stats$r_squared <- NA
        }
    } else {
        fit_stats$r_squared <- NA
    }

    # True error statistics (if available)
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted) && length(object$y_true_sorted) == length(object$predictions)) {
        true_residuals <- object$y_true_sorted - object$predictions
        true_error_stats <- list(
            true_mse = mean(true_residuals^2),
            true_rmse = sqrt(mean(true_residuals^2)),
            true_mae = mean(abs(true_residuals)),
            true_median_ae = median(abs(true_residuals))
        )
    }

    # Cross-validation statistics
    cv_stats <- NULL
    if (!is.null(object$h_cv_errors) && length(object$h_cv_errors) > 0) {
        cv_stats <- list(
            cv_errors = object$h_cv_errors,
            mean_cv_error = mean(object$h_cv_errors, na.rm = TRUE),
            min_cv_error = min(object$h_cv_errors, na.rm = TRUE)
        )

        # Add h_values if available
        if (!is.null(object$h_values)) {
            cv_stats$h_values <- object$h_values
        }
    }

    # Bootstrap information
    bootstrap_info <- NULL
    if (!is.null(object$bb_predictions) || !is.null(object$opt_ci_lower)) {
        bootstrap_info <- list(
            n_bootstrap = if (!is.null(object$n.bb)) {
                object$n.bb
            } else if (!is.null(object$bb_predictions)) {
                if (is.matrix(object$bb_predictions)) {
                    ncol(object$bb_predictions)
                } else {
                    length(object$bb_predictions)
                }
            } else {
                NA
            },
            has_ci = !is.null(object$opt_ci_lower) && !is.null(object$opt_ci_upper)
        )

        if (bootstrap_info$has_ci &&
            length(object$opt_ci_upper) == length(object$opt_ci_lower)) {
            bootstrap_info$mean_ci_width <-
                mean(object$opt_ci_upper - object$opt_ci_lower, na.rm = TRUE)
        }
    }

    # Residual analysis
    residual_stats <- list(
        mean = mean(residuals),
        sd = sd(residuals),
        min = min(residuals),
        quantiles = quantile(residuals, c(0.25, 0.5, 0.75)),
        max = max(residuals)
    )

    # Add Shapiro-Wilk test only if we have enough observations
    n_resid <- length(residuals)
    if (n_resid >= 3 && n_resid <= 5000) {
        tryCatch({
            residual_stats$shapiro_test <- shapiro.test(residuals)
        }, error = function(e) {
            # If shapiro.test fails, just skip it
            residual_stats$shapiro_test <- NULL
        })
    } else {
        residual_stats$shapiro_test <- NULL
    }

    # Create summary object
    result <- list(
        call = object$call,
        model_info = model_info,
        fit_stats = fit_stats,
        true_error_stats = true_error_stats,
        cv_stats = cv_stats,
        residual_stats = residual_stats,
        bootstrap_info = bootstrap_info
    )

    class(result) <- "summary.ugkmm"
    return(result)
}

#' Print summary.ugkmm Objects
#'
#' @param x An object of class "summary.ugkmm".
#' @param digits Number of significant digits.
#' @param ... Additional arguments (unused).
#'
#' @method print summary.ugkmm
#' @export
print.summary.ugkmm <- function(x, digits = 4, ...) {
    cat("\nAdaptive Neighborhood Size Graph K-Means Regression\n")
    cat(rep("=", 60), "\n", sep = "")

    # Print call if available
    if (!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }

    # Model info
    cat("\nModel Information:\n")
    cat("  Observations:     ", x$model_info$n, "\n")
    if (!is.na(x$model_info$optimal_h)) {
        cat("  Optimal h:        ", x$model_info$optimal_h, "\n")
    }
    if (!any(is.na(x$model_info$x_range))) {
        cat("  X range:          [",
            format(x$model_info$x_range[1], digits = digits), ", ",
            format(x$model_info$x_range[2], digits = digits), "]\n", sep = "")
    }
    cat("  Central tendency: ", x$model_info$method, "\n")

    # Fit statistics
    cat("\nFit Statistics:\n")
    cat("  MSE:              ", format(x$fit_stats$mse, digits = digits), "\n")
    cat("  RMSE:             ", format(x$fit_stats$rmse, digits = digits), "\n")
    cat("  MAE:              ", format(x$fit_stats$mae, digits = digits), "\n")
    cat("  Median AE:        ", format(x$fit_stats$median_ae, digits = digits), "\n")
    if (!is.na(x$fit_stats$r_squared)) {
        cat("  R-squared:        ", format(x$fit_stats$r_squared, digits = digits), "\n")
    }

    # True errors if available
    if (!is.null(x$true_error_stats)) {
        cat("\nTrue Error Statistics:\n")
        cat("  True MSE:         ", format(x$true_error_stats$true_mse, digits = digits), "\n")
        cat("  True RMSE:        ", format(x$true_error_stats$true_rmse, digits = digits), "\n")
        cat("  True MAE:         ", format(x$true_error_stats$true_mae, digits = digits), "\n")
        cat("  True Median AE:   ", format(x$true_error_stats$true_median_ae, digits = digits), "\n")
    }

    # Cross-validation statistics
    if (!is.null(x$cv_stats)) {
        cat("\nCross-validation Statistics:\n")
        if (!is.na(x$cv_stats$mean_cv_error)) {
            cat("  Mean CV Error:    ", format(x$cv_stats$mean_cv_error, digits = digits), "\n")
        }
        if (!is.na(x$cv_stats$min_cv_error)) {
            cat("  Min CV Error:     ", format(x$cv_stats$min_cv_error, digits = digits), "\n")
        }
        if (!is.null(x$cv_stats$h_values) && length(x$cv_stats$h_values) > 0) {
            h_range <- range(x$cv_stats$h_values)
            cat("  h range tested:   [", h_range[1], ", ", h_range[2], "]\n", sep = "")
        }
    }

    # Residual diagnostics
    cat("\nResidual Diagnostics:\n")
    cat("  Mean:             ", format(x$residual_stats$mean, digits = digits), "\n")
    cat("  SD:               ", format(x$residual_stats$sd, digits = digits), "\n")
    cat("  Range:            [",
        format(x$residual_stats$min, digits = digits), ", ",
        format(x$residual_stats$max, digits = digits), "]\n", sep = "")

    # Print quantiles
    cat("  Quantiles:\n")
    q_names <- names(x$residual_stats$quantiles)
    for (i in seq_along(x$residual_stats$quantiles)) {
        cat("    ", format(q_names[i], width = 7), " : ",
            format(x$residual_stats$quantiles[i], digits = digits), "\n")
    }

    # Shapiro-Wilk test (if available)
    if (!is.null(x$residual_stats$shapiro_test)) {
        cat("  Normality test:\n")
        cat("    Shapiro-Wilk W = ",
            format(x$residual_stats$shapiro_test$statistic, digits = digits),
            ", p-value = ",
            format(x$residual_stats$shapiro_test$p.value, digits = digits), "\n")
    }

    # Bootstrap info
    if (!is.null(x$bootstrap_info)) {
        cat("\nBootstrap Information:\n")
        if (!is.na(x$bootstrap_info$n_bootstrap)) {
            cat("  Bootstrap samples:", x$bootstrap_info$n_bootstrap, "\n")
        }
        cat("  Credible intervals: ",
            if (x$bootstrap_info$has_ci) "Available" else "Not available", "\n")
        if (x$bootstrap_info$has_ci && !is.null(x$bootstrap_info$mean_ci_width)) {
            cat("  Mean CI width:    ",
                format(x$bootstrap_info$mean_ci_width, digits = digits), "\n")
        }
    }

    cat("\n")
    invisible(x)
}

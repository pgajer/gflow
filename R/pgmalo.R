#' Path Graph Model Averaging Local Linear Model
#'
#' @description
#' Fits a path graph model averaging local linear model using cross-validation
#' to select optimal neighborhood sizes. The function implements adaptive
#' bandwidth selection for graph-structured data with optional bootstrap
#' confidence intervals.
#'
#' @param neighbors List of integer vectors containing adjacency lists for each
#'   vertex. Each element should contain valid vertex indices (1 to n).
#' @param edge_lengths List of numeric vectors containing positive edge lengths
#'   corresponding to each adjacency. Must have same structure as neighbors.
#' @param y Numeric vector of response variables for each vertex.
#' @param y.true Optional numeric vector of true response values for computing
#'   prediction errors. Must have same length as y.
#' @param use.median Logical; if TRUE uses median instead of mean for bootstrap
#'   interval estimation (default: TRUE).
#' @param h.min Integer; minimum neighborhood size to consider. Must be even
#'   and at least 4 (default: 4).
#' @param h.max Integer; maximum neighborhood size to consider. Must be even
#'   and less than n-2 where n is the number of vertices (default: min(30, n-2)).
#' @param p Numeric; confidence level for bootstrap intervals. Must be between
#'   0 and 1 (default: 0.95).
#' @param n.bb Integer; number of bootstrap iterations. Set to 0 to skip
#'   bootstrap (default: 50).
#' @param bb.max.distance.deviation Integer; maximum distance deviation allowed
#'   in bootstrap samples. Must be >= -1 (default: 1).
#' @param n.CVs Integer; number of cross-validation iterations (default: 100).
#' @param n.CV.folds Integer; number of cross-validation folds. Must be between
#'   2 and n (default: 10).
#' @param seed Integer; random seed for reproducibility (default: 0).
#' @param kernel.type Integer between 1 and 10 specifying the kernel function:
#'   \itemize{
#'     \item 1 = Gaussian
#'     \item 2 = Epanechnikov
#'     \item 3 = Quartic
#'     \item 4 = Triweight
#'     \item 5 = Tricubic
#'     \item 6 = Cosine
#'     \item 7 = Uniform (default)
#'     \item 8-10 = Additional kernels (see documentation)
#'   }
#' @param dist.normalization.factor Numeric; factor for normalizing distances.
#'   Must be greater than 1 (default: 1.1).
#' @param epsilon Numeric; numerical stability threshold. Must be positive
#'   (default: 1e-15).
#' @param verbose Logical; if TRUE prints progress messages (default: FALSE).
#'
#' @return An S3 object of class "pgmalo" containing:
#' \describe{
#'   \item{h_values}{Integer vector of tested neighborhood sizes}
#'   \item{opt_h}{Optimal neighborhood size selected by cross-validation}
#'   \item{opt_h_idx}{Index of optimal h in h_values vector}
#'   \item{h_cv_errors}{Numeric vector of cross-validation errors for each h}
#'   \item{true_errors}{Numeric vector of true prediction errors (if y.true provided)}
#'   \item{predictions}{Numeric vector of predictions using optimal h}
#'   \item{local_predictions}{Numeric vector of locally adaptive predictions}
#'   \item{h_predictions}{List of prediction vectors for each h value}
#'   \item{bb_predictions}{Matrix of bootstrap predictions (n.bb x n)}
#'   \item{ci_lower}{Numeric vector of lower confidence bounds}
#'   \item{ci_upper}{Numeric vector of upper confidence bounds}
#'   \item{has_bootstrap}{Logical indicating if bootstrap was performed}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The Path Graph Model Averaging Local Linear (PGMALO) method performs
#' local linear regression on graph-structured data. For each vertex, it
#' considers neighborhoods of different sizes and uses cross-validation to
#' select the optimal neighborhood size globally or locally.
#'
#' The method accounts for the graph structure through path distances and
#' applies kernel weighting based on these distances. Bootstrap sampling
#' can be used to quantify prediction uncertainty.
#'
#' @examples
#' \donttest{
#' # Create a simple chain graph
#' n <- 50
#' neighbors <- vector("list", n)
#' edge_lengths <- vector("list", n)
#'
#' # Build chain structure
#' for(i in 1:n) {
#'   if(i == 1) {
#'     neighbors[[i]] <- 2L
#'     edge_lengths[[i]] <- 1.0
#'   } else if(i == n) {
#'     neighbors[[i]] <- (n-1L)
#'     edge_lengths[[i]] <- 1.0
#'   } else {
#'     neighbors[[i]] <- c(i-1L, i+1L)
#'     edge_lengths[[i]] <- c(1.0, 1.0)
#'   }
#' }
#'
#' # Generate smooth response with noise
#' x_pos <- seq(0, 1, length.out = n)
#' y <- sin(2 * pi * x_pos) + rnorm(n, 0, 0.1)
#'
#' # Fit model with small h range for speed
#'   fit <- pgmalo(neighbors, edge_lengths, y, h.max = 10, n.CVs = 10)
#'   summary(fit)
#'   print(fit)
#' }
#'
#' @seealso
#' \code{\link{upgmalo}} for univariate version,
#' \code{\link{plot.pgmalo}} for visualization
#'
#' @references
#' Goldstein, L., Chu, B., Nayak, S., and Minsker, S. (2024).
#' "Path Graph Model Averaging Local Linear Estimation."
#' \emph{Journal of Statistical Software} (forthcoming).
#'
#' @export
pgmalo <- function(neighbors,
                   edge_lengths,
                   y,
                   y.true = NULL,
                   use.median = TRUE,
                   h.min = 4L,
                   h.max = min(30L, length(y) - 2L),
                   p = 0.95,
                   n.bb = 50L,
                   bb.max.distance.deviation = 1L,
                   n.CVs = 100L,
                   n.CV.folds = 10L,
                   seed = 0L,
                   kernel.type = 7L,
                   dist.normalization.factor = 1.1,
                   epsilon = 1e-15,
                   verbose = FALSE) {

    # Store the call for the return object
    cl <- match.call()

    # Input validation
    if (!is.list(neighbors)) {
        stop("'neighbors' must be a list", call. = FALSE)
    }
    if (!is.list(edge_lengths)) {
        stop("'edge_lengths' must be a list", call. = FALSE)
    }

    n <- length(y)

    if (length(neighbors) != n) {
        stop("'neighbors' must have same length as 'y'", call. = FALSE)
    }
    if (length(edge_lengths) != n) {
        stop("'edge_lengths' must have same length as 'y'", call. = FALSE)
    }

    # Validate y
    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }
    if (any(is.na(y))) {
        stop("'y' cannot contain NA values", call. = FALSE)
    }
    if (any(!is.finite(y))) {
        stop("'y' cannot contain infinite values", call. = FALSE)
    }

    # Validate neighbors and edge_lengths structure
    for (i in seq_len(n)) {
        if (!is.numeric(neighbors[[i]]) || !is.vector(neighbors[[i]])) {
            stop(sprintf("'neighbors[[%d]]' must be a numeric vector", i), call. = FALSE)
        }
        if (!is.numeric(edge_lengths[[i]]) || !is.vector(edge_lengths[[i]])) {
            stop(sprintf("'edge_lengths[[%d]]' must be a numeric vector", i), call. = FALSE)
        }
        if (length(neighbors[[i]]) != length(edge_lengths[[i]])) {
            stop(sprintf("'neighbors[[%d]]' and 'edge_lengths[[%d]]' must have same length", i, i),
                 call. = FALSE)
        }
        if (length(neighbors[[i]]) > 0) {
            if (any(neighbors[[i]] != as.integer(neighbors[[i]]))) {
                stop(sprintf("'neighbors[[%d]]' must contain integer values", i), call. = FALSE)
            }
            if (any(neighbors[[i]] < 1L | neighbors[[i]] > n)) {
                stop(sprintf("'neighbors[[%d]]' must contain indices between 1 and %d", i, n),
                     call. = FALSE)
            }
            if (any(edge_lengths[[i]] <= 0)) {
                stop(sprintf("'edge_lengths[[%d]]' must contain positive values", i), call. = FALSE)
            }
            if (any(!is.finite(edge_lengths[[i]]))) {
                stop(sprintf("'edge_lengths[[%d]]' cannot contain non-finite values", i), call. = FALSE)
            }
        }
    }

    # Validate y.true if provided
    if (!is.null(y.true)) {
        if (!is.numeric(y.true) || !is.vector(y.true)) {
            stop("'y.true' must be a numeric vector", call. = FALSE)
        }
        if (length(y.true) != n) {
            stop("'y.true' must have same length as 'y'", call. = FALSE)
        }
        if (any(!is.finite(y.true))) {
            stop("'y.true' cannot contain non-finite values", call. = FALSE)
        }
    }

    # Validate scalar parameters
    use.median <- as.logical(use.median)[1]
    if (is.na(use.median)) {
        stop("'use.median' must be TRUE or FALSE", call. = FALSE)
    }

    # Validate h parameters
    h.min <- as.integer(h.min)[1]
    h.max <- as.integer(h.max)[1]
    if (is.na(h.min) || h.min < 4L || h.min %% 2L != 0L) {
        stop("'h.min' must be an even integer >= 4", call. = FALSE)
    }
    if (is.na(h.max) || h.max <= h.min || h.max %% 2L != 0L) {
        stop("'h.max' must be an even integer > h.min", call. = FALSE)
    }
    if (h.max > n - 2L) {
        stop(sprintf("'h.max' cannot exceed %d (n - 2)", n - 2L), call. = FALSE)
    }

    # Validate p
    p <- as.numeric(p)[1]
    if (is.na(p) || p <= 0 || p >= 1) {
        stop("'p' must be a number between 0 and 1", call. = FALSE)
    }

    # Validate bootstrap parameters
    n.bb <- as.integer(n.bb)[1]
    if (is.na(n.bb) || n.bb < 0L) {
        stop("'n.bb' must be a non-negative integer", call. = FALSE)
    }

    bb.max.distance.deviation <- as.integer(bb.max.distance.deviation)[1]
    if (is.na(bb.max.distance.deviation) || bb.max.distance.deviation < -1L) {
        stop("'bb.max.distance.deviation' must be an integer >= -1", call. = FALSE)
    }

    # Validate CV parameters
    n.CVs <- as.integer(n.CVs)[1]
    if (is.na(n.CVs) || n.CVs < 0L) {
        stop("'n.CVs' must be a non-negative integer", call. = FALSE)
    }

    n.CV.folds <- as.integer(n.CV.folds)[1]
    if (is.na(n.CV.folds) || n.CV.folds < 2L || n.CV.folds > n) {
        stop(sprintf("'n.CV.folds' must be an integer between 2 and %d", n), call. = FALSE)
    }

    # Validate other parameters
    seed <- as.integer(seed)[1]
    if (is.na(seed)) {
        stop("'seed' must be an integer", call. = FALSE)
    }

    kernel.type <- as.integer(kernel.type)[1]
    if (is.na(kernel.type) || kernel.type < 1L || kernel.type > 10L) {
        stop("'kernel.type' must be an integer between 1 and 10", call. = FALSE)
    }

    dist.normalization.factor <- as.numeric(dist.normalization.factor)[1]
    if (is.na(dist.normalization.factor) || dist.normalization.factor <= 1) {
        stop("'dist.normalization.factor' must be a number > 1", call. = FALSE)
    }

    epsilon <- as.numeric(epsilon)[1]
    if (is.na(epsilon) || epsilon <= 0) {
        stop("'epsilon' must be a positive number", call. = FALSE)
    }

    verbose <- as.logical(verbose)[1]
    if (is.na(verbose)) {
        stop("'verbose' must be TRUE or FALSE", call. = FALSE)
    }

    # Convert to 0-based indexing for C++
    neighbors.0based <- lapply(neighbors, function(x) {
        if (length(x) > 0) as.integer(x - 1L) else integer(0)
    })

    # Call C++ implementation
    result <- .Call("S_pgmalo",
                    neighbors.0based,
                    edge_lengths,
                    as.double(y),
                    if (is.null(y.true)) double() else as.double(y.true),
                    use.median,
                    h.min,
                    h.max,
                    p,
                    n.bb,
                    bb.max.distance.deviation,
                    n.CVs,
                    n.CV.folds,
                    seed,
                    kernel.type,
                    dist.normalization.factor,
                    epsilon,
                    verbose)

    # Add call to result
    result$call <- cl

    # Set class
    class(result) <- "pgmalo"

    return(result)
}

#' Summary Method for pgmalo Objects
#'
#' @description
#' Produces summary statistics for pgmalo model fits including model parameters,
#' fit quality measures, and cross-validation results.
#'
#' @param object An object of class "pgmalo" from \code{\link{pgmalo}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "summary.pgmalo" containing:
#' \describe{
#'   \item{call}{The matched call}
#'   \item{n_vertices}{Number of vertices in the graph}
#'   \item{optimal_h}{Optimal neighborhood size selected by CV}
#'   \item{cv_info}{Cross-validation statistics including range and errors}
#'   \item{bootstrap_info}{Bootstrap information if available}
#'   \item{true_error_info}{True error statistics if y.true was provided}
#'   \item{prediction_summary}{Summary statistics of predictions}
#' }
#'
#' @details
#' The summary includes basic model information, cross-validation results,
#' and bootstrap information if available. When true values are provided,
#' additional error metrics are computed.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' n <- 50
#' neighbors <- vector("list", n)
#' edge_lengths <- vector("list", n)
#' for(i in 1:n) {
#'   if(i == 1) {
#'     neighbors[[i]] <- 2L
#'     edge_lengths[[i]] <- 1.0
#'   } else if(i == n) {
#'     neighbors[[i]] <- (n-1L)
#'     edge_lengths[[i]] <- 1.0
#'   } else {
#'     neighbors[[i]] <- c(i-1L, i+1L)
#'     edge_lengths[[i]] <- c(1.0, 1.0)
#'   }
#' }
#' y <- rnorm(n)
#'
#' # Fit model and get summary
#' fit <- pgmalo(neighbors, edge_lengths, y)
#' summary(fit)
#' }
#'
#' @seealso
#' \code{\link{pgmalo}} for model fitting,
#' \code{\link{print.summary.pgmalo}} for formatted output
#'
#' @export
summary.pgmalo <- function(object, ...) {

    cat("In summary.pgmalo  \n\n")

    # Input validation
    if (!inherits(object, "pgmalo")) {
        stop("'object' must be a pgmalo object", call. = FALSE)
    }

    # Basic information
    n_vertices <- if (!is.null(object$predictions)) {
        length(object$predictions)
    } else if (!is.null(object$y)) {
        length(object$y)
    } else {
        NA
    }

    # Cross-validation information
    cv_info <- NULL
    if (!is.null(object$h_values) && !is.null(object$h_errors)) {
        opt_idx <- which(object$h_values == object$opt_h)
        cv_info <- list(
            h_tested = length(object$h_values),
            h_range = range(object$h_values),
            h_values = object$h_values,
            cv_errors = object$h_errors,
            min_error = min(object$h_errors),
            opt_error = if (length(opt_idx) > 0) object$h_errors[opt_idx[1]] else NA,
            relative_errors = object$h_errors / min(object$h_errors)
        )
    } else if (!is.null(object$h_values) && !is.null(object$h_cv_errors)) {
        # Handle alternative naming
        opt_idx <- which(object$h_values == object$opt_h)
        cv_info <- list(
            h_tested = length(object$h_values),
            h_range = range(object$h_values),
            h_values = object$h_values,
            cv_errors = object$h_cv_errors,
            min_error = min(object$h_cv_errors),
            opt_error = if (length(opt_idx) > 0) object$h_cv_errors[opt_idx[1]] else NA,
            relative_errors = object$h_cv_errors / min(object$h_cv_errors)
        )
    }

    # Bootstrap information
    has_bootstrap <- FALSE
    if (!is.null(object$has_bootstrap)) {
        has_bootstrap <- object$has_bootstrap
    } else if (!is.null(object$ci_lower) && !is.null(object$ci_upper)) {
        has_bootstrap <- TRUE
    }

    bootstrap_info <- NULL
    if (has_bootstrap && !is.null(object$ci_lower) && !is.null(object$ci_upper)) {
        ci_widths <- object$ci_upper - object$ci_lower
        bootstrap_info <- list(
            available = TRUE,
            n_bootstrap = if (!is.null(object$n_bb)) object$n_bb else NA,
            confidence_level = if (!is.null(object$p) && length(object$p) == 1) object$p else 0.95,
            mean_width = mean(ci_widths, na.rm = TRUE),
            median_width = median(ci_widths, na.rm = TRUE),
            min_width = min(ci_widths, na.rm = TRUE),
            max_width = max(ci_widths, na.rm = TRUE)
        )

        # Coverage if true values available
        if (!is.null(object$y_true) && length(object$y_true) == length(object$ci_lower)) {
            bootstrap_info$empirical_coverage <- mean(
                object$y_true >= object$ci_lower &
                object$y_true <= object$ci_upper,
                na.rm = TRUE
            )
        }
    }

    # True error information
    true_error_info <- NULL
    if (!is.null(object$true_errors) && length(object$true_errors) > 0) {
        true_error_info <- list(
            available = TRUE,
            mean_error = mean(object$true_errors, na.rm = TRUE),
            median_error = median(object$true_errors, na.rm = TRUE),
            sd_error = sd(object$true_errors, na.rm = TRUE),
            mse = mean(object$true_errors^2, na.rm = TRUE),
            rmse = sqrt(mean(object$true_errors^2, na.rm = TRUE)),
            mae = mean(abs(object$true_errors), na.rm = TRUE),
            median_ae = median(abs(object$true_errors), na.rm = TRUE)
        )
    }

    # Prediction summary
    prediction_summary <- NULL
    if (!is.null(object$predictions)) {
        prediction_summary <- list(
            mean = mean(object$predictions, na.rm = TRUE),
            median = median(object$predictions, na.rm = TRUE),
            sd = sd(object$predictions, na.rm = TRUE),
            range = range(object$predictions, na.rm = TRUE),
            quantiles = quantile(object$predictions,
                                 probs = c(0.25, 0.5, 0.75),
                                 na.rm = TRUE)
        )
    }

    # Local predictions information
    local_pred_info <- NULL
    if (!is.null(object$local_predictions) && length(object$local_predictions) > 0) {
        local_pred_info <- list(
            available = TRUE,
            differs_from_global = !isTRUE(all.equal(object$predictions,
                                                     object$local_predictions))
        )
    }

    # Create summary object
    structure(
        list(
            call = object$call,
            n_vertices = n_vertices,
            optimal_h = object$opt_h,
            optimal_h_idx = object$opt_h_idx,
            cv_info = cv_info,
            bootstrap_info = bootstrap_info,
            true_error_info = true_error_info,
            prediction_summary = prediction_summary,
            local_pred_info = local_pred_info
        ),
        class = "summary.pgmalo"
    )
}

#' Print Method for summary.pgmalo Objects
#'
#' @description
#' Displays formatted summary of pgmalo model fits with organized sections
#' for model information, cross-validation results, and diagnostic statistics.
#'
#' @param x An object of class "summary.pgmalo" from \code{\link{summary.pgmalo}}.
#' @param digits Number of significant digits for numeric output (default: 4).
#' @param show_cv_details Logical; if TRUE shows detailed CV results (default: FALSE).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object invisibly.
#'
#' @examples
#' \dontrun{
#' fit <- pgmalo(neighbors, edge_lengths, y)
#' summary(fit)
#'
#' # Show more CV details
#' print(summary(fit), show_cv_details = TRUE)
#' }
#'
#' @export
print.summary.pgmalo <- function(x, digits = 4, show_cv_details = FALSE, ...) {
    # Helper function for separator lines
    hr <- function(width = 65) {
        cat(paste0(rep("-", width), collapse = ""), "\n")
    }

    # Header
    cat("\nPGMALO Model Summary\n")
    cat("====================\n")

    # Call
    if (!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }

    # Model Information section
    cat("\nModel Information:\n")
    hr()
    if (!is.na(x$n_vertices)) {
        cat("  Number of vertices:     ", x$n_vertices, "\n")
    }
    cat("  Optimal h:              ", x$optimal_h, "\n")
    if (!is.null(x$optimal_h_idx)) {
        cat("  Optimal h index:        ", x$optimal_h_idx, "\n")
    }

    # Cross-validation section
    if (!is.null(x$cv_info)) {
        cat("\nCross-validation Results:\n")
        hr()
        cat("  h values tested:        ", x$cv_info$h_tested, "\n")
        cat("  h range:                ",
            sprintf("[%d, %d]", x$cv_info$h_range[1], x$cv_info$h_range[2]), "\n")
        cat("  Minimum CV error:       ",
            format(x$cv_info$min_error, digits = digits), "\n")
        if (!is.na(x$cv_info$opt_error)) {
            cat("  Optimal h CV error:     ",
                format(x$cv_info$opt_error, digits = digits), "\n")
            cat("  Relative CV error:      ",
                format(x$cv_info$opt_error / x$cv_info$min_error, digits = digits), "\n")
        }

        if (show_cv_details && !is.null(x$cv_info$h_values)) {
            cat("\n  Detailed CV errors:\n")
            cv_table <- data.frame(
                h = x$cv_info$h_values,
                CV_Error = format(x$cv_info$cv_errors, digits = digits),
                Relative = format(x$cv_info$relative_errors, digits = digits)
            )
            print(cv_table, row.names = FALSE)
        }
    }

    # Bootstrap information section
    if (!is.null(x$bootstrap_info) && x$bootstrap_info$available) {
        cat("\nBootstrap Confidence Intervals:\n")
        hr()
        cat("  Status:                 Available\n")
        if (!is.na(x$bootstrap_info$confidence_level)) {
            cat("  Confidence level:       ",
                sprintf("%.1f%%", 100 * x$bootstrap_info$confidence_level), "\n")
        }
        if (!is.na(x$bootstrap_info$n_bootstrap)) {
            cat("  Bootstrap samples:      ", x$bootstrap_info$n_bootstrap, "\n")
        }
        cat("  Mean CI width:          ",
            format(x$bootstrap_info$mean_width, digits = digits), "\n")
        cat("  Median CI width:        ",
            format(x$bootstrap_info$median_width, digits = digits), "\n")
        cat("  CI width range:         ",
            sprintf("[%.4f, %.4f]", x$bootstrap_info$min_width, x$bootstrap_info$max_width), "\n")

        if (!is.null(x$bootstrap_info$empirical_coverage)) {
            cat("  Empirical coverage:     ",
                sprintf("%.1f%%", 100 * x$bootstrap_info$empirical_coverage), "\n")
        }
    } else if (!is.null(x$bootstrap_info)) {
        cat("\nBootstrap Confidence Intervals: Not available\n")
    }

    # True error section
    if (!is.null(x$true_error_info) && x$true_error_info$available) {
        cat("\nTrue Error Analysis:\n")
        hr()
        cat("  Mean error:             ",
            format(x$true_error_info$mean_error, digits = digits), "\n")
        cat("  SD of errors:           ",
            format(x$true_error_info$sd_error, digits = digits), "\n")
        cat("  MSE:                    ",
            format(x$true_error_info$mse, digits = digits), "\n")
        cat("  RMSE:                   ",
            format(x$true_error_info$rmse, digits = digits), "\n")
        cat("  MAE:                    ",
            format(x$true_error_info$mae, digits = digits), "\n")
        cat("  Median AE:              ",
            format(x$true_error_info$median_ae, digits = digits), "\n")
    }

    # Prediction summary section
    if (!is.null(x$prediction_summary)) {
        cat("\nPrediction Summary:\n")
        hr()
        cat("  Mean:                   ",
            format(x$prediction_summary$mean, digits = digits), "\n")
        cat("  SD:                     ",
            format(x$prediction_summary$sd, digits = digits), "\n")
        cat("  Range:                  ",
            sprintf("[%.4f, %.4f]",
                    x$prediction_summary$range[1],
                    x$prediction_summary$range[2]), "\n")
        cat("  Quartiles (Q1,Q2,Q3):   ",
            sprintf("%.4f, %.4f, %.4f",
                    x$prediction_summary$quantiles[1],
                    x$prediction_summary$quantiles[2],
                    x$prediction_summary$quantiles[3]), "\n")
    }

    # Local predictions information
    if (!is.null(x$local_pred_info) && x$local_pred_info$available) {
        cat("\nLocal Predictions:        ",
            ifelse(x$local_pred_info$differs_from_global,
                   "Available (differ from global)",
                   "Available (same as global)"), "\n")
    }

    cat("\n")
    invisible(x)
}

#' Univariate Path Graph Model Averaging Local Linear Model
#'
#' @description
#' Implements adaptive neighborhood selection for univariate data using path
#' graph model averaging local linear (UPGMALO) estimation. The function
#' automatically constructs a chain graph from sorted predictor values and
#' applies the PGMALO methodology.
#'
#' @param x Numeric vector of predictor values. Will be automatically sorted
#'   if not already in ascending order.
#' @param y Numeric vector of response values corresponding to x.
#' @param y.true Optional numeric vector of true response values for computing
#'   prediction errors.
#' @param use.median Logical; if TRUE uses median instead of mean for bootstrap
#'   interval estimation (default: TRUE).
#' @param h.min Integer; minimum neighborhood size to consider. Must be even
#'   and at least 4 (default: 4).
#' @param h.max Integer; maximum neighborhood size to consider. Must be even
#'   and at most n-2 (default: min(30, n-2)).
#' @param p Numeric; confidence level for bootstrap intervals (default: 0.95).
#' @param n.bb Integer; number of bootstrap iterations (default: 50).
#' @param bb.max.distance.deviation Integer; maximum distance deviation for
#'   bootstrap samples (default: 1).
#' @param n.CVs Integer; number of cross-validation iterations (default: 100).
#' @param n.CV.folds Integer; number of cross-validation folds (default: 10).
#' @param seed Integer; random seed for reproducibility (default: 0).
#' @param kernel.type Integer between 1 and 10 specifying kernel function
#'   (default: 7 for uniform kernel).
#' @param dist.normalization.factor Numeric; distance normalization factor
#'   (default: 1.01).
#' @param epsilon Numeric; numerical stability parameter (default: 1e-15).
#' @param verbose Logical; if TRUE prints progress messages (default: FALSE).
#'
#' @return An S3 object of class "upgmalo" containing all components from
#' \code{\link{pgmalo}} plus:
#' \describe{
#'   \item{x_sorted}{Sorted predictor values}
#'   \item{y_sorted}{Response values corresponding to sorted x}
#'   \item{y_true_sorted}{True values corresponding to sorted x (if provided)}
#'   \item{h_min}{Minimum h value used}
#'   \item{h_max}{Maximum h value used}
#'   \item{max_h_index}{Index range for h values}
#' }
#'
#' @details
#' This function is a convenience wrapper that constructs a chain graph from
#' univariate data and applies the PGMALO methodology. The chain graph connects
#' each point to its immediate neighbors in the sorted order of x values.
#'
#' @examples
#' \dontrun{
#' # Generate nonlinear data
#' n <- 50
#' x <- seq(0, 2*pi, length.out = n)
#' y_true <- sin(x)
#' y <- y_true + rnorm(n, 0, 0.2)
#'
#' # Fit model
#'   fit <- upgmalo(x, y, y.true = y_true, h.max = 20, n.CVs = 20)
#'
#'   # Plot results
#'   plot(fit)
#' }
#'
#' @seealso
#' \code{\link{pgmalo}} for general graph version,
#' \code{\link{plot.upgmalo}} for visualization,
#' \code{\link{summary.upgmalo}} for model summaries
#'
#' @export
upgmalo <- function(x,
                    y,
                    y.true = NULL,
                    use.median = TRUE,
                    h.min = 4L,
                    h.max = min(30L, length(x) - 2L),
                    p = 0.95,
                    n.bb = 50L,
                    bb.max.distance.deviation = 1L,
                    n.CVs = 100L,
                    n.CV.folds = 10L,
                    seed = 0L,
                    kernel.type = 7L,
                    dist.normalization.factor = 1.01,
                    epsilon = 1e-15,
                    verbose = FALSE) {

    # Store the call
    cl <- match.call()

    # Input validation
    if (!is.numeric(x) || !is.vector(x)) {
        stop("'x' must be a numeric vector", call. = FALSE)
    }
    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    n <- length(x)

    if (length(y) != n) {
        stop("'x' and 'y' must have the same length", call. = FALSE)
    }
    if (n < 4L) {
        stop("Need at least 4 observations", call. = FALSE)
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop("'x' and 'y' cannot contain NA values", call. = FALSE)
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("'x' and 'y' cannot contain infinite values", call. = FALSE)
    }

    # Check if sorting is needed
    if (!identical(x, sort(x))) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
        if (!is.null(y.true)) {
            if (!is.numeric(y.true) || !is.vector(y.true)) {
                stop("'y.true' must be a numeric vector", call. = FALSE)
            }
            if (length(y.true) != n) {
                stop("'y.true' must have same length as 'x'", call. = FALSE)
            }
            if (any(!is.finite(y.true))) {
                stop("'y.true' cannot contain non-finite values", call. = FALSE)
            }
            y.true <- y.true[ord]
        }
    } else if (!is.null(y.true)) {
        # Still validate y.true even if no sorting needed
        if (!is.numeric(y.true) || !is.vector(y.true)) {
            stop("'y.true' must be a numeric vector", call. = FALSE)
        }
        if (length(y.true) != n) {
            stop("'y.true' must have same length as 'x'", call. = FALSE)
        }
        if (any(!is.finite(y.true))) {
            stop("'y.true' cannot contain non-finite values", call. = FALSE)
        }
    }

    # Additional validation specific to upgmalo
    # Ensure h.max is reasonable for univariate case
    h.max <- min(h.max, n - 2L)
    if (h.max %% 2L != 0L) {
        h.max <- h.max - 1L  # Make it even
    }

    # Call C++ implementation for univariate case
    result <- .Call("S_upgmalo",
                    as.double(x),
                    as.double(y),
                    if (is.null(y.true)) double() else as.double(y.true),
                    use.median,
                    as.integer(h.min),
                    as.integer(h.max),
                    p,
                    as.integer(n.bb),
                    as.integer(bb.max.distance.deviation),
                    as.integer(n.CVs),
                    as.integer(n.CV.folds),
                    as.integer(seed),
                    as.integer(kernel.type),
                    dist.normalization.factor,
                    epsilon,
                    verbose)

    # Add additional information
    result$x_sorted <- x
    result$y_sorted <- y
    result$y_true_sorted <- y.true
    result$h_min <- h.min
    result$h_max <- h.max
    result$max_h_index <- (h.max - h.min) %/% 2L
    result$call <- cl

    # Set class
    class(result) <- c("upgmalo", "pgmalo")

    return(result)
}

#' Plot Method for pgmalo Objects
#'
#' @description
#' Creates diagnostic and visualization plots for pgmalo model fits.
#'
#' @param x An object of class "pgmalo".
#' @param type Character string specifying plot type (default: "diagnostic").
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return NULL (invisibly). Called for side effect of creating plots.
#'
#' @details
#' Currently only implements diagnostic plot showing cross-validation errors
#' across different h values. For more detailed plotting options, use
#' \code{\link{plot.upgmalo}} with upgmalo objects.
#'
#' @examples
#' # See examples in pgmalo()
#'
#' @export
plot.pgmalo <- function(x, type = "diagnostic", ...) {
    if (!inherits(x, "pgmalo")) {
        stop("'x' must be a pgmalo object", call. = FALSE)
    }

    if (type == "diagnostic") {
        # Check for either h_cv_errors or h_errors
        h_errors <- if (!is.null(x$h_cv_errors)) x$h_cv_errors else x$h_errors

        if (is.null(x$h_values) || is.null(h_errors)) {
            stop("Diagnostic plot requires h_values and h_errors/h_cv_errors", call. = FALSE)
        }

        plot(x$h_values, h_errors,
             type = 'b',
             las = 1,
             xlab = "Neighborhood Size (h)",
             ylab = "Cross-validation Error",
             main = "Cross-validation Diagnostic Plot",
             ...)

        # Add vertical line at optimal h
        abline(v = x$opt_h, col = "red", lty = 2)

        # Add text annotation
        text(x$opt_h, min(h_errors),
             sprintf("Optimal h = %d", x$opt_h),
             pos = 4, col = "red")
    } else {
        stop("Currently only 'diagnostic' plot type is implemented for pgmalo objects",
             call. = FALSE)
    }

    invisible(NULL)
}

#' Plot Method for upgmalo Objects
#'
#' @description
#' Creates comprehensive visualizations for univariate PGMALO results including
#' fitted values with confidence bands, diagnostic plots, and residual analysis.
#'
#' @param x An object of class "upgmalo" returned by \code{\link{upgmalo}}.
#' @param type Character string specifying the plot type:
#'   \describe{
#'     \item{"fit"}{Fitted values with optional confidence bands (default)}
#'     \item{"diagnostic"}{Cross-validation error by neighborhood size}
#'     \item{"residuals"}{Residual scatter plot}
#'     \item{"residuals_hist"}{Histogram of residuals}
#'   }
#' @param title Main title for the plot (default: "").
#' @param xlab X-axis label (default: "x" for most plots).
#' @param ylab Y-axis label (default: "y" for fit plot, appropriate label for others).
#' @param predictions.type Character string specifying which predictions to plot:
#'   \describe{
#'     \item{"predictions"}{Global optimal h predictions (default)}
#'     \item{"bb.predictions"}{Bootstrap-based predictions}
#'     \item{"local.predictions"}{Locally adaptive predictions}
#'   }
#' @param with.y.true Logical; include true values if available (default: TRUE).
#' @param with.pts Logical; show original data points (default: FALSE).
#' @param with.CrI Logical; show credible/confidence intervals if available (default: TRUE).
#' @param y.true.col Color for true values line (default: "red").
#' @param predictions.col Color for predictions line (default: "blue").
#' @param with.predictions.pts Logical; show prediction points (default: FALSE).
#' @param predictions.pts.col Color for prediction points (default: "blue").
#' @param predictions.pts.pch Point character for predictions (default: 20).
#' @param CrI.as.polygon Logical; show confidence band as polygon (default: TRUE).
#' @param CrI.polygon.col Color for confidence band polygon (default: "gray85").
#' @param CrI.line.col Color for confidence band lines if not polygon (default: "gray10").
#' @param CrI.line.lty Line type for confidence bands if not polygon (default: 2).
#' @param ylim Y-axis limits (default: NULL for automatic).
#' @param ... Additional graphical parameters passed to base plotting functions.
#'
#' @return NULL (invisibly). Called for side effect of creating plots.
#'
#' @details
#' The function provides multiple visualization options:
#' \itemize{
#'   \item Fit plots show the estimated function with optional confidence bands
#'   \item Diagnostic plots help assess the cross-validation process
#'   \item Residual plots aid in model checking
#' }
#'
#' Confidence intervals are only shown for "predictions" and "bb.predictions"
#' types when bootstrap was performed (n.bb > 0 in the original call).
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' n <- 100
#' x <- seq(0, 2*pi, length.out = n)
#' y <- sin(x) + rnorm(n, 0, 0.2)
#'
#' # Fit model
#' fit <- upgmalo(x, y, n.bb = 100)
#'
#' # Create various plots
#' plot(fit, type = "fit", with.pts = TRUE)
#' plot(fit, type = "diagnostic")
#' plot(fit, type = "residuals")
#' }
#'
#' @seealso
#' \code{\link{upgmalo}} for model fitting,
#' \code{\link{summary.upgmalo}} for numerical summaries
#'
#' @export
plot.upgmalo <- function(x,
                         type = c("fit", "diagnostic", "residuals", "residuals_hist"),
                         title = "",
                         xlab = "",
                         ylab = "",
                         predictions.type = c("predictions", "bb.predictions", "local.predictions"),
                         with.y.true = TRUE,
                         with.pts = FALSE,
                         with.CrI = TRUE,
                         y.true.col = "red",
                         predictions.col = "blue",
                         with.predictions.pts = FALSE,
                         predictions.pts.col = "blue",
                         predictions.pts.pch = 20,
                         CrI.as.polygon = TRUE,
                         CrI.polygon.col = "gray85",
                         CrI.line.col = "gray10",
                         CrI.line.lty = 2,
                         ylim = NULL,
                         ...) {

    # Input validation
    if (!inherits(x, "upgmalo")) {
        stop("'x' must be a upgmalo object", call. = FALSE)
    }

    # Match arguments
    type <- match.arg(type)
    predictions.type <- match.arg(predictions.type)

    # Set default labels
    if (xlab == "") {
        xlab <- if (type %in% c("fit", "residuals")) "x" else {
            if (type == "diagnostic") "Neighborhood Size (h)" else "Residuals"
        }
    }

    if (ylab == "") {
        ylab <- switch(type,
                       fit = "y",
                       diagnostic = "Cross-validation Error",
                       residuals = "Residuals",
                       residuals_hist = "Frequency")
    }

    # Call appropriate plotting function
    switch(type,
           fit = .plot.fit.upgmalo(x, title, xlab, ylab, predictions.type,
                                   with.y.true, with.pts, with.CrI,
                                   y.true.col, predictions.col,
                                   with.predictions.pts, predictions.pts.col,
                                   predictions.pts.pch, CrI.as.polygon,
                                   CrI.polygon.col, CrI.line.col,
                                   CrI.line.lty, ylim, ...),
           diagnostic = .plot.diagnostic.upgmalo(x, title, xlab, ylab, ...),
           residuals = .plot.residuals.upgmalo(x, title, xlab, ylab,
                                              predictions.type, ...),
           residuals_hist = .plot.residuals.hist.upgmalo(x, title, xlab, ylab,
                                                         predictions.type, ...))

    invisible(NULL)
}

# Internal plotting functions (not exported)
.plot.fit.upgmalo <- function(res, title, xlab, ylab, predictions.type,
                              with.y.true, with.pts, with.CrI,
                              y.true.col, predictions.col,
                              with.predictions.pts, predictions.pts.col,
                              predictions.pts.pch, CrI.as.polygon,
                              CrI.polygon.col, CrI.line.col,
                              CrI.line.lty, ylim, ...) {

    # Get appropriate predictions
    predictions <- switch(predictions.type,
                          predictions = res$opt_predictions,
                          bb.predictions = {
                              if (is.null(res$opt_bb_predictions)) {
                                  stop("Bootstrap predictions not available (set n.bb > 0)",
                                       call. = FALSE)
                              }
                              res$opt_bb_predictions
                          },
                          local.predictions = {
                              if (is.null(res$opt_local_predictions)) {
                                  stop("Local predictions not available", call. = FALSE)
                              }
                              res$opt_local_predictions
                          })

    # Determine y-axis limits
    if (is.null(ylim)) {
        ylim_vals <- c(res$y_sorted, predictions)
        if (with.CrI && !is.null(res$opt_ci_lower) &&
            predictions.type %in% c("predictions", "bb.predictions")) {
            ylim_vals <- c(ylim_vals, res$opt_ci_lower, res$opt_ci_upper)
        }
        if (with.y.true && !is.null(res$y_true_sorted)) {
            ylim_vals <- c(ylim_vals, res$y_true_sorted)
        }
        ylim <- range(ylim_vals, na.rm = TRUE)
        ylim <- ylim + c(-1, 1) * 0.05 * diff(ylim)  # Add 5% margin
    }

    # Create base plot
    plot(res$x_sorted, res$y_sorted,
         type = "n",
         las = 1,
         ylim = ylim,
         xlab = xlab,
         ylab = ylab,
         main = title,
         ...)

    # Add confidence bands
    if (with.CrI && !is.null(res$opt_ci_lower) &&
        predictions.type %in% c("predictions", "bb.predictions")) {
        if (CrI.as.polygon) {
            polygon(c(res$x_sorted, rev(res$x_sorted)),
                    c(res$opt_ci_lower, rev(res$opt_ci_upper)),
                    col = CrI.polygon.col,
                    border = NA)
        } else {
            lines(res$x_sorted, res$opt_ci_lower,
                  col = CrI.line.col, lty = CrI.line.lty)
            lines(res$x_sorted, res$opt_ci_upper,
                  col = CrI.line.col, lty = CrI.line.lty)
        }
    }

    # Add prediction line
    lines(res$x_sorted, predictions, col = predictions.col, lwd = 2)

    # Add data points
    if (with.pts) {
        points(res$x_sorted, res$y_sorted, pch = 19, cex = 0.5)
    }

    # Add prediction points
    if (with.predictions.pts) {
        points(res$x_sorted, predictions,
               col = predictions.pts.col,
               pch = predictions.pts.pch,
               cex = 0.8)
    }

    # Add true values
    if (with.y.true && !is.null(res$y_true_sorted)) {
        lines(res$x_sorted, res$y_true_sorted,
              col = y.true.col, lwd = 2, lty = 2)
    }

    # Create legend
    legend_items <- character()
    legend_cols <- character()
    legend_ltys <- numeric()
    legend_lwds <- numeric()
    legend_pchs <- numeric()

    # Add items to legend
    pred_label <- switch(predictions.type,
                         predictions = "Global Estimate",
                         bb.predictions = "Bootstrap Estimate",
                         local.predictions = "Local Estimate")

    legend_items <- c(legend_items, pred_label)
    legend_cols <- c(legend_cols, predictions.col)
    legend_ltys <- c(legend_ltys, 1)
    legend_lwds <- c(legend_lwds, 2)
    legend_pchs <- c(legend_pchs, NA)

    if (with.y.true && !is.null(res$y_true_sorted)) {
        legend_items <- c(legend_items, "True Function")
        legend_cols <- c(legend_cols, y.true.col)
        legend_ltys <- c(legend_ltys, 2)
        legend_lwds <- c(legend_lwds, 2)
        legend_pchs <- c(legend_pchs, NA)
    }

    if (with.pts) {
        legend_items <- c(legend_items, "Data")
        legend_cols <- c(legend_cols, "black")
        legend_ltys <- c(legend_ltys, NA)
        legend_lwds <- c(legend_lwds, NA)
        legend_pchs <- c(legend_pchs, 19)
    }

    if (with.CrI && !is.null(res$opt_ci_lower) &&
        predictions.type %in% c("predictions", "bb.predictions")) {
        conf_level <- if (!is.null(res$p)) res$p else 0.95
        legend_items <- c(legend_items, sprintf("%.0f%% CI", 100 * conf_level))
        legend_cols <- c(legend_cols,
                          if (CrI.as.polygon) CrI.polygon.col else CrI.line.col)
        legend_ltys <- c(legend_ltys, if (CrI.as.polygon) 1 else CrI.line.lty)
        legend_lwds <- c(legend_lwds, if (CrI.as.polygon) 10 else 1)
        legend_pchs <- c(legend_pchs, if (CrI.as.polygon) 15 else NA)
    }

    # Place legend
    legend("topright",
           legend = legend_items,
           col = legend_cols,
           lty = legend_ltys,
           lwd = legend_lwds,
           pch = legend_pchs,
           bg = "white",
           box.lty = 1,
           cex = 0.9)
}

.plot.diagnostic.upgmalo <- function(res, title, xlab, ylab, ...) {
    if (is.null(res$h_values) || is.null(res$h_cv_errors)) {
        stop("Diagnostic plot requires h_values and h_cv_errors", call. = FALSE)
    }

    # Create plot
    plot(res$h_values, res$h_cv_errors,
         type = 'b',
         las = 1,
         xlab = xlab,
         ylab = ylab,
         main = if (title == "") "Cross-validation Diagnostic" else title,
         pch = 19,
         ...)

    # Mark optimal h
    opt_idx <- which(res$h_values == res$opt_h)
    if (length(opt_idx) > 0) {
        points(res$opt_h, res$h_cv_errors[opt_idx],
               pch = 19, cex = 1.5, col = "red")
        abline(v = res$opt_h, col = "red", lty = 2)

        # Add annotation
        text(res$opt_h,
             min(res$h_cv_errors) + 0.1 * diff(range(res$h_cv_errors)),
             sprintf("Optimal h = %d", res$opt_h),
             pos = 4, col = "red")
    }

    # Add grid
    grid(col = "gray90")
}

.plot.residuals.upgmalo <- function(res, title, xlab, ylab, predictions.type, ...) {
    # Get predictions
    predictions <- switch(predictions.type,
                          predictions = res$opt_predictions,
                          bb.predictions = res$opt_bb_predictions,
                          local.predictions = res$opt_local_predictions)

    if (is.null(predictions)) {
        stop(sprintf("%s predictions not available", predictions.type), call. = FALSE)
    }

    residuals <- res$y_sorted - predictions

    # Create plot
    plot(res$x_sorted, residuals,
         type = "p",
         las = 1,
         xlab = xlab,
         ylab = ylab,
         main = if (title == "") "Residual Plot" else title,
         pch = 19,
         cex = 0.7,
         ...)

    # Add reference line at zero
    abline(h = 0, col = "red", lty = 2, lwd = 2)

    # Add smooth line
    if (length(residuals) >= 10) {
        lo <- try(loess(residuals ~ res$x_sorted, degree = 2), silent = TRUE)
        if (!inherits(lo, "try-error")) {
            lines(res$x_sorted, fitted(lo), col = "blue", lwd = 2)
        }
    }

    # Add grid
    grid(col = "gray90")
}

.plot.residuals.hist.upgmalo <- function(res, title, xlab, ylab, predictions.type, ...) {
    # Get predictions
    predictions <- switch(predictions.type,
                          predictions = res$opt_predictions,
                          bb.predictions = res$opt_bb_predictions,
                          local.predictions = res$opt_local_predictions)

    if (is.null(predictions)) {
        stop(sprintf("%s predictions not available", predictions.type), call. = FALSE)
    }

    residuals <- res$y_sorted - predictions

    # Create histogram
    hist(residuals,
         breaks = "FD",  # Freedman-Diaconis rule
         col = "lightblue",
         border = "darkblue",
         las = 1,
         xlab = xlab,
         main = if (title == "") "Residual Distribution" else title,
         ...)

    # Add normal curve overlay
    x_seq <- seq(min(residuals), max(residuals), length.out = 100)
    y_norm <- dnorm(x_seq, mean = mean(residuals), sd = sd(residuals))

    # Scale to histogram
    hist_obj <- hist(residuals, breaks = "FD", plot = FALSE)
    y_norm <- y_norm * diff(hist_obj$breaks)[1] * length(residuals)

    lines(x_seq, y_norm, col = "red", lwd = 2)

    # Add vertical line at mean
    abline(v = mean(residuals), col = "darkgreen", lty = 2, lwd = 2)

    # Add legend
    legend("topright",
           legend = c("Normal Curve", "Mean"),
           col = c("red", "darkgreen"),
           lty = c(1, 2),
           lwd = 2,
           bg = "white")

    # Add text with statistics
    mtext(sprintf("Mean = %.3f, SD = %.3f", mean(residuals), sd(residuals)),
          side = 3, line = -2, adj = 0.05, cex = 0.9)
}


#' Summary Method for upgmalo Objects
#'
#' @description
#' Produces comprehensive summary statistics for UPGMALO model fits including
#' fit quality measures, cross-validation results, residual diagnostics, and
#' optional comparison with true values.
#'
#' @param object An object of class "upgmalo" from \code{\link{upgmalo}}.
#' @param quantiles Numeric vector of probabilities for residual quantiles
#'   (default: c(0, 0.25, 0.5, 0.75, 1)).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "summary.upgmalo" containing:
#' \describe{
#'   \item{call}{The matched call}
#'   \item{model_info}{Basic model information}
#'   \item{fit_stats}{Fit quality statistics (MSE, RMSE, MAE)}
#'   \item{cv_stats}{Cross-validation summary}
#'   \item{residual_stats}{Residual distribution summary}
#'   \item{true_error_stats}{True error analysis (if y.true was provided)}
#'   \item{bootstrap_info}{Bootstrap summary (if n.bb > 0)}
#' }
#'
#' @details
#' When true values are available, the summary includes additional diagnostics:
#' \itemize{
#'   \item Equivalent normal standard deviation
#'   \item Shapiro-Wilk test for normality of standardized errors
#'   \item Relative efficiency compared to optimal linear estimator
#' }
#'
#' @examples
#' # See examples in upgmalo()
#'
#' @seealso
#' \code{\link{upgmalo}} for model fitting,
#' \code{\link{print.summary.upgmalo}} for formatted output
#'
#' @export
summary.upgmalo <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    # Input validation
    if (!inherits(object, "upgmalo")) {
        stop("'object' must be a upgmalo object", call. = FALSE)
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("'quantiles' must be probabilities between 0 and 1", call. = FALSE)
    }

    # Check for potential issues
    if (any(is.na(object$y_sorted)) || any(is.na(object$opt_predictions))) {
        warning("Missing values detected in model results", call. = FALSE)
    }

    # Calculate residuals
    residuals <- object$y_sorted - object$opt_predictions

    # Basic model information
    model_info <- list(
        n_observations = length(object$x_sorted),
        optimal_h = object$opt_h,
        h_range = c(object$h_min, object$h_max),
        x_range = range(object$x_sorted),
        n_cv_folds = if (!is.null(object$n_CV_folds)) object$n_CV_folds else NA,
        n_cv_iterations = if (!is.null(object$n_CVs)) object$n_CVs else NA
    )

    # Fit statistics
    fit_stats <- list(
        mse = mean(residuals^2),
        rmse = sqrt(mean(residuals^2)),
        mae = mean(abs(residuals)),
        median_ae = median(abs(residuals)),
        r_squared = 1 - var(residuals) / var(object$y_sorted)
    )

    # Cross-validation statistics
    cv_stats <- NULL
    if (!is.null(object$h_values) && !is.null(object$h_cv_errors)) {
        cv_stats <- list(
            h_tested = length(object$h_values),
            min_cv_error = min(object$h_cv_errors),
            opt_cv_error = object$h_cv_errors[which(object$h_values == object$opt_h)[1]],
            cv_error_range = range(object$h_cv_errors)
        )
    }

    # Residual statistics
    residual_stats <- list(
        mean = mean(residuals),
        sd = sd(residuals),
        skewness = mean(((residuals - mean(residuals)) / sd(residuals))^3),
        kurtosis = mean(((residuals - mean(residuals)) / sd(residuals))^4) - 3,
        quantiles = setNames(quantile(residuals, probs = quantiles),
                             paste0(100 * quantiles, "%"))
    )

    # True error analysis
    true_error_stats <- NULL
    if (!is.null(object$y_true_sorted)) {
        true_errors <- object$y_true_sorted - object$opt_predictions
        true_mse <- mean(true_errors^2)

        # Standardized errors
        std_errors <- true_errors / sqrt(true_mse)

        # Normality test (only if n >= 3)
        shapiro_p <- if (length(std_errors) >= 3 && length(std_errors) <= 5000) {
            shapiro.test(std_errors)$p.value
        } else NA

        true_error_stats <- list(
            true_mse = true_mse,
            true_rmse = sqrt(true_mse),
            true_mae = mean(abs(true_errors)),
            equiv_normal_sd = sqrt(true_mse),
            shapiro_p_value = shapiro_p,
            relative_efficiency = mean(abs(true_errors)) / (sqrt(2/pi) * sqrt(true_mse))
        )
    }

    # Bootstrap information
    bootstrap_info <- NULL
    if (!is.null(object$opt_ci_lower) && !is.null(object$opt_ci_upper)) {
        ci_widths <- object$opt_ci_upper - object$opt_ci_lower
        bootstrap_info <- list(
            has_intervals = TRUE,
            confidence_level = if (!is.null(object$p)) object$p else 0.95,
            mean_ci_width = mean(ci_widths),
            median_ci_width = median(ci_widths),
            coverage = if (!is.null(object$y_true_sorted)) {
                mean(object$y_true_sorted >= object$opt_ci_lower &
                     object$y_true_sorted <= object$opt_ci_upper)
            } else NA
        )
    }

    # Create summary object
    structure(
        list(
            call = object$call,
            model_info = model_info,
            fit_stats = fit_stats,
            cv_stats = cv_stats,
            residual_stats = residual_stats,
            true_error_stats = true_error_stats,
            bootstrap_info = bootstrap_info
        ),
        class = "summary.upgmalo"
    )
}

#' Print Method for summary.upgmalo Objects
#'
#' @description
#' Formats and displays UPGMALO model summaries in a clean, organized format.
#'
#' @param x An object of class "summary.upgmalo".
#' @param digits Number of significant digits for numeric output (default: 4).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object invisibly.
#'
#' @examples
#' # See examples in upgmalo() and summary.upgmalo()
#'
#' @export
print.summary.upgmalo <- function(x, digits = 4, ...) {
    # Helper for horizontal rules
    hr <- function(char = "-", width = 70) {
        cat(rep(char, width), sep = "", "\n")
    }

    # Header
    cat("\n")
    hr("=")
    cat("UPGMALO Model Summary\n")
    hr("=")

    # Call
    cat("\nCall:\n")
    print(x$call)

    # Model Information
    cat("\n")
    hr()
    cat("Model Information:\n")
    hr()
    cat("Number of observations:", x$model_info$n_observations, "\n")
    cat("Optimal h:            ", x$model_info$optimal_h, "\n")
    cat("h range tested:       ", sprintf("[%d, %d]",
        x$model_info$h_range[1], x$model_info$h_range[2]), "\n")
    cat("x range:              ", sprintf("[%.3f, %.3f]",
        x$model_info$x_range[1], x$model_info$x_range[2]), "\n")

    if (!is.na(x$model_info$n_cv_iterations)) {
        cat("CV iterations:        ", x$model_info$n_cv_iterations, "\n")
        cat("CV folds:             ", x$model_info$n_cv_folds, "\n")
    }

    # Fit Statistics
    cat("\n")
    hr()
    cat("Fit Statistics:\n")
    hr()
    cat("MSE:                  ", format(x$fit_stats$mse, digits = digits), "\n")
    cat("RMSE:                 ", format(x$fit_stats$rmse, digits = digits), "\n")
    cat("MAE:                  ", format(x$fit_stats$mae, digits = digits), "\n")
    cat("Median AE:            ", format(x$fit_stats$median_ae, digits = digits), "\n")
    cat("R-squared:            ", format(x$fit_stats$r_squared, digits = digits), "\n")

    # Cross-validation Statistics
    if (!is.null(x$cv_stats)) {
        cat("\n")
        hr()
        cat("Cross-validation Results:\n")
        hr()
        cat("h values tested:      ", x$cv_stats$h_tested, "\n")
        cat("Minimum CV error:     ", format(x$cv_stats$min_cv_error, digits = digits), "\n")
        cat("Optimal h CV error:   ", format(x$cv_stats$opt_cv_error, digits = digits), "\n")
        cat("CV error range:       ", sprintf("[%.4f, %.4f]",
            x$cv_stats$cv_error_range[1], x$cv_stats$cv_error_range[2]), "\n")
    }

    # Residual Statistics
    cat("\n")
    hr()
    cat("Residual Analysis:\n")
    hr()
    cat("Mean:                 ", format(x$residual_stats$mean, digits = digits), "\n")
    cat("Std. deviation:       ", format(x$residual_stats$sd, digits = digits), "\n")
    cat("Skewness:             ", format(x$residual_stats$skewness, digits = digits), "\n")
    cat("Excess kurtosis:      ", format(x$residual_stats$kurtosis, digits = digits), "\n")
    cat("\nQuantiles:\n")
    print(x$residual_stats$quantiles, digits = digits)

    # True Error Statistics
    if (!is.null(x$true_error_stats)) {
        cat("\n")
        hr()
        cat("True Error Analysis:\n")
        hr()
        cat("True MSE:             ", format(x$true_error_stats$true_mse, digits = digits), "\n")
        cat("True RMSE:            ", format(x$true_error_stats$true_rmse, digits = digits), "\n")
        cat("True MAE:             ", format(x$true_error_stats$true_mae, digits = digits), "\n")
        cat("Equivalent Normal SD: ", format(x$true_error_stats$equiv_normal_sd, digits = digits), "\n")

        if (!is.na(x$true_error_stats$shapiro_p_value)) {
            cat("Normality (p-value):  ",
                format(x$true_error_stats$shapiro_p_value, digits = digits), "\n")
        }
        cat("Relative efficiency:  ",
            format(x$true_error_stats$relative_efficiency, digits = digits), "\n")
    }

    # Bootstrap Information
    if (!is.null(x$bootstrap_info) && x$bootstrap_info$has_intervals) {
        cat("\n")
        hr()
        cat("Bootstrap Confidence Intervals:\n")
        hr()
        cat("Confidence level:     ",
            sprintf("%.1f%%", 100 * x$bootstrap_info$confidence_level), "\n")
        cat("Mean CI width:        ",
            format(x$bootstrap_info$mean_ci_width, digits = digits), "\n")
        cat("Median CI width:      ",
            format(x$bootstrap_info$median_ci_width, digits = digits), "\n")

        if (!is.na(x$bootstrap_info$coverage)) {
            cat("Empirical coverage:   ",
                sprintf("%.1f%%", 100 * x$bootstrap_info$coverage), "\n")
        }
    }

    cat("\n")
    invisible(x)
}

#' Print Method for pgmalo Objects
#'
#' @description
#' Displays a concise summary of a pgmalo model fit.
#'
#' @param x An object of class "pgmalo".
#' @param digits Number of significant digits for numeric output (default: 3).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object invisibly.
#'
#' @export
print.pgmalo <- function(x, digits = 3, ...) {
    cat("Path Graph Model Averaging Local Linear (PGMALO) Fit\n")
    cat("----------------------------------------------------\n")

    if (!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }

    cat("\nOptimal neighborhood size (h):", x$opt_h, "\n")

    if (!is.null(x$h_cv_errors) && length(x$h_cv_errors) > 0) {
        cat("Cross-validation error:       ",
            format(x$h_cv_errors[which(x$h_values == x$opt_h)[1]],
                   digits = digits), "\n")
    }

    if (!is.null(x$predictions)) {
        cat("Number of vertices:           ", length(x$predictions), "\n")
    }

    has_bootstrap <- FALSE
    if (!is.null(x$has_bootstrap)) {
        has_bootstrap <- x$has_bootstrap
    } else if (!is.null(x$ci_lower) && !is.null(x$ci_upper)) {
        has_bootstrap <- TRUE
    }

    if (has_bootstrap && !is.null(x$ci_lower)) {
        cat("Bootstrap intervals:           Available\n")
    }

    cat("\n")
    invisible(x)
}

#' Print Method for upgmalo Objects
#'
#' @description
#' Displays a concise summary of a upgmalo model fit.
#'
#' @param x An object of class "upgmalo".
#' @param digits Number of significant digits for numeric output (default: 3).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object invisibly.
#'
#' @export
print.upgmalo <- function(x, digits = 3, ...) {
    cat("Univariate PGMALO Fit\n")
    cat("---------------------\n")

    if (!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }

    cat("\nNumber of observations:      ", length(x$x_sorted), "\n")
    cat("Optimal neighborhood size:   ", x$opt_h, "\n")
    cat("Range of x:                  ",
        sprintf("[%.3f, %.3f]", min(x$x_sorted), max(x$x_sorted)), "\n")

    # Fit quality
    if (!is.null(x$opt_predictions)) {
        residuals <- x$y_sorted - x$opt_predictions
        cat("RMSE:                        ",
            format(sqrt(mean(residuals^2)), digits = digits), "\n")
    }

    # True error if available
    if (!is.null(x$y_true_sorted) && !is.null(x$opt_predictions)) {
        true_errors <- x$y_true_sorted - x$opt_predictions
        cat("True RMSE:                   ",
            format(sqrt(mean(true_errors^2)), digits = digits), "\n")
    }

    if (!is.null(x$opt_ci_lower)) {
        cat("Confidence intervals:         Available\n")
    }

    cat("\n")
    invisible(x)
}


#' Adaptive Model Averaged GEodesic LOcal linear smoothing
#'
#' @description
#' Performs nonparametric smoothing of 1D data using adaptive model averaged local linear regression.
#' AMAGELO implements a grid-based approach with model averaging and automatic bandwidth selection
#' to produce smooth predictions while adapting to local data characteristics. The method also
#' identifies local extrema and provides measures of their significance.
#'
#' @param x Numeric vector of predictor values
#' @param y Numeric vector of response values
#' @param grid.size Integer specifying the number of grid points (default: 100)
#' @param min.bw.factor Numeric minimum bandwidth factor (default: 0.01)
#' @param max.bw.factor Numeric maximum bandwidth factor (default: 0.5)
#' @param n.bws Integer number of bandwidths to evaluate (default: 30)
#' @param use.global.bw.grid Logical: use same bandwidth grid for all points? (default: TRUE)
#' @param with.bw.predictions Logical: return predictions for all bandwidths? (default: FALSE)
#' @param log.grid Logical: use logarithmic bandwidth spacing? (default: TRUE)
#' @param domain.min.size Integer minimum data points per local model (default: 10)
#' @param kernel.type Integer code for kernel function (default: 1 = Gaussian)
#' @param dist.normalization.factor Numeric scale factor for distances (default: 1)
#' @param n.cleveland.iterations Integer robustness iterations (default: 3)
#' @param blending.coef Numeric model blending coefficient (default: 0.5)
#' @param use.linear.blending Logical: use linear blending? (default: FALSE)
#' @param precision Numeric precision threshold for optimization (default: 1e-6)
#' @param small.depth.threshold Numeric threshold for wiggle detection (default: 0.05)
#' @param depth.similarity.tol Numeric tolerance for depth similarity (default: 0.001)
#' @param verbose Logical: print progress information? (default: FALSE)
#'
#' @return A list containing:
#'   \item{x_sorted}{Sorted predictor values}
#'   \item{y_sorted}{Response values ordered by sorted x}
#'   \item{order}{Original indices of sorted data (1-based)}
#'   \item{grid_coords}{x-coordinates of grid points}
#'   \item{predictions}{Fitted values at optimal bandwidth}
#'   \item{bw_predictions}{Matrix of predictions by bandwidth (if requested)}
#'   \item{grid_predictions}{Predictions at grid points}
#'   \item{harmonic_predictions}{Predictions after triplet harmonic smoothing}
#'   \item{local_extrema}{Matrix of detected extrema with columns:
#'     idx (Index in sorted data, 1-based),
#'     x (x-coordinate of extremum),
#'     y (y-value at extremum),
#'     is_max (1 if maximum, 0 if minimum),
#'     depth (Vertical prominence of extremum),
#'     depth_idx (Index where min/max descent terminates, 1-based),
#'     rel_depth (Depth relative to total depth of all extrema),
#'     range_rel_depth (Depth relative to range of predictions)}
#'   \item{monotonic_interval_proportions}{Relative lengths of monotonic intervals}
#'   \item{change_scaled_monotonicity_index}{Weighted signed average of directional changes, quantifying monotonicity strength and directionality; values close to \eqn{+1} or \eqn{-1} indicate strong global monotonic trends.}
#'   \item{bw_errors}{Cross-validation errors for each bandwidth}
#'   \item{opt_bw_idx}{Index of optimal bandwidth (1-based)}
#'   \item{min_bw}{Minimum bandwidth value}
#'   \item{max_bw}{Maximum bandwidth value}
#'   \item{bws}{Vector of evaluated bandwidths}
#'
#' @details
#' AMAGELO constructs a uniform grid over the data domain and fits local linear models
#' at each grid point using a range of bandwidths. Models are then averaged with weights
#' that depend on both spatial proximity and model quality. The optimal bandwidth is
#' selected by cross-validation.
#'
#' The method identifies local extrema (maxima and minima) and calculates their prominence
#' using both absolute and relative depth measures. It also provides measures of overall
#' monotonicity through the TVMI and Simpson index.
#'
#' Triplet harmonic smoothing is applied to remove small wiggles while preserving
#' significant features of the data.
#'
#' @examples
#' \dontrun{
#' # Simulate data with smooth trend and noise
#' x <- seq(0, 10, length.out = 200)
#' y <- sin(x) + 0.2 * rnorm(length(x))
#'
#' # Apply AMAGELO smoothing
#' result <- amagelo(x, y, grid_size = 100)
#'
#' # Plot results
#' plot(x, y, pch = 16, col = "gray")
#' lines(x[result$order], result$predictions, col = "red", lwd = 2)
#'
#' # Examine local extrema
#' extrema <- result$local_extrema
#' points(extrema[,"x"], extrema[,"y"],
#'        pch = ifelse(extrema[,"is_max"] == 1, 24, 25),
#'        bg = ifelse(extrema[,"is_max"] == 1, "red", "blue"),
#'        cex = 2 * extrema[,"range_rel_depth"])
#' }
#'
#' @export
amagelo <- function(
                    x,
                    y,
                    grid.size,
                    min.bw.factor,
                    max.bw.factor,
                    n.bws,
                    use.global.bw.grid = TRUE,
                    with.bw.predictions = FALSE,
                    log.grid = FALSE,
                    domain.min.size = 4,
                    kernel.type = 7L,
                    dist.normalization.factor = 1.1,
                    n.cleveland.iterations = 1,
                    blending.coef = 0,
                    use.linear.blending = TRUE,
                    precision = 1e-6,
                    small.depth.threshold = 0.05,
                    depth.similarity.tol = 0.0001,
                    verbose = FALSE
                    ) {

    ##-- Input validation --------------------------------------------------------
    if (!is.numeric(x))      stop("x must be numeric")
    if (!is.numeric(y))      stop("y must be numeric")
    if (length(x) != length(y))
        stop("x and y must have the same length")
    if (!is.numeric(grid.size) || grid.size <= 0)
        stop("grid.size must be a positive integer")
    if (!is.numeric(min.bw.factor) || min.bw.factor <= 0)
        stop("min.bw.factor must be positive")
    if (!is.numeric(max.bw.factor) || max.bw.factor <= min.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")
    if (!is.numeric(n.bws) || n.bws <= 0)
        stop("n.bws must be a positive integer")
    if (!is.logical(use.global.bw.grid))
        stop("use.global.bw.grid must be logical")
    if (!is.logical(with.bw.predictions))
        stop("with.bw.predictions must be logical")
    if (!is.logical(log.grid))
        stop("log.grid must be logical")
    if (!is.numeric(domain.min.size) || domain.min.size <= 0)
        stop("domain.min.size must be a positive integer")
    if (!is.numeric(kernel.type) || kernel.type < 0)
        stop("kernel.type must be a non-negative integer code")
    if (!is.numeric(dist.normalization.factor) || dist.normalization.factor <= 0)
        stop("dist.normalization.factor must be positive")
    if (!is.numeric(n.cleveland.iterations) || n.cleveland.iterations <= 0)
        stop("n.cleveland.iterations must be a positive integer")
    if (!is.numeric(blending.coef) || blending.coef < 0 || blending.coef > 1)
        stop("blending.coef must be in [0,1]")
    if (!is.logical(use.linear.blending))
        stop("use.linear.blending must be logical")
    if (!is.numeric(precision) || precision <= 0)
        stop("precision must be positive")
    if (!is.numeric(small.depth.threshold) || small.depth.threshold <= 0)
        stop("small.depth.threshold must be positive")
    if (!is.numeric(depth.similarity.tol) || depth.similarity.tol <= 0)
        stop("depth.similarity.tol must be positive")
    if (!is.logical(verbose))
        stop("verbose must be logical")

    ##-- Call C backend ---------------------------------------------------------
    result <- .Call(
        "S_amagelo",
        as.double(x),
        as.double(y),
        as.integer(grid.size),
        as.double(min.bw.factor),
        as.double(max.bw.factor),
        as.integer(n.bws),
        as.logical(use.global.bw.grid),
        as.logical(with.bw.predictions),
        as.logical(log.grid),
        as.integer(domain.min.size),
        as.integer(kernel.type),
        as.double(dist.normalization.factor),
        as.integer(n.cleveland.iterations),
        as.double(blending.coef),
        as.logical(use.linear.blending),
        as.double(precision),
        as.double(small.depth.threshold),
        as.double(depth.similarity.tol),
        as.logical(verbose)
    )

    ## assign class for S3
    class(result) <- "amagelo"
    result
}

#' Print method for amagelo objects
#'
#' @param x An object of class 'amagelo'
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments passed to print
#' @export
print.amagelo <- function(x, digits = 4, ...) {
  cat("AMAGELO: Adaptive Model Averaged GEodesic LOcal linear smoothing\n")
  cat("----------------------------------------------------------\n")

  # Basic information
  cat(sprintf("Number of data points: %d\n", length(x$x_sorted)))
  cat(sprintf("Grid size: %d\n", length(x$grid_coords)))

  # Bandwidth information
  cat("\nBandwidth Information:\n")
  cat(sprintf("  Optimal bandwidth: %.4g (index %d of %d)\n",
              x$bws[x$opt_bw_idx], x$opt_bw_idx, length(x$bws)))
  cat(sprintf("  Bandwidth range: [%.4g, %.4g]\n", x$min_bw, x$max_bw))
  cat(sprintf("  Min error: %.4g (at optimal bandwidth)\n",
              x$bw_errors[x$opt_bw_idx]))

  # Extrema information
  n_maxima <- sum(x$local_extrema[,"is_max"] == 1)
  n_minima <- sum(x$local_extrema[,"is_max"] == 0)
  cat("\nLocal Extrema:\n")
  cat(sprintf("  Number of local maxima: %d\n", n_maxima))
  cat(sprintf("  Number of local minima: %d\n", n_minima))

  # Monotonicity measures
  cat("\nMonotonicity Measures:\n")
  ## cat(sprintf("  TVMI (Total-Variation Monotonicity Index): %.4g\n", x$tvmi))
  ## cat(sprintf("  Simpson index: %.4g\n", x$simpson_index))
  cat(sprintf("  Change-Scaled Monotonicity Index: %.4g\n", x$change_scaled_monotonicity_index))
  cat(sprintf("  Number of monotonic intervals: %d\n",
              length(x$monotonic_interval_proportions)))

  # Data range information
  y_range <- range(x$predictions)
  cat("\nPrediction Range:\n")
  cat(sprintf("  Min: %.4g, Max: %.4g, Range: %.4g\n",
              y_range[1], y_range[2], diff(y_range)))

  cat("\nCall:\n")
  print(attr(x, "call"))

  invisible(x)
}

#' Summary method for amagelo objects
#'
#' @param object An object of class 'amagelo'
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments passed to summary
#' @export
summary.amagelo <- function(object, digits = 4, ...) {
  res <- list()

  # Basic information
  res$n_points <- length(object$x_sorted)
  res$grid_size <- length(object$grid_coords)
  res$x_range <- range(object$x_sorted)
  res$y_range <- range(object$y_sorted)

  # Bandwidth information
  res$opt_bw <- object$bws[object$opt_bw_idx]
  res$bw_range <- c(object$min_bw, object$max_bw)
  res$min_error <- object$bw_errors[object$opt_bw_idx]

  # Prediction statistics
  res$pred_stats <- c(
    Min = min(object$predictions),
    Q1 = quantile(object$predictions, 0.25),
    Median = median(object$predictions),
    Mean = mean(object$predictions),
    Q3 = quantile(object$predictions, 0.75),
    Max = max(object$predictions),
    SD = sd(object$predictions)
  )

  # Extrema information
  res$n_maxima <- sum(object$harmonic_predictions_local_extrema[,"is_max"] == 1)
  res$n_minima <- sum(object$harmonic_predictions_local_extrema[,"is_max"] == 0)

  # Top extrema by prominence
  top_n <- min(5, nrow(object$harmonic_predictions_local_extrema))
  sorted_extrema <- object$harmonic_predictions_local_extrema[order(-object$harmonic_predictions_local_extrema[,"range_rel_depth"]), , drop = FALSE]
  res$top_extrema <- sorted_extrema[1:top_n, c("x", "y", "is_max", "range_rel_depth"), drop = FALSE]

  # Monotonicity measures
  ## res$tvmi <- object$tvmi
  ## res$simpson_index <- object$simpson_index
  res$change_scaled_monotonicity_index <- object$change_scaled_monotonicity_index
  res$monotonic_intervals <- length(object$monotonic_interval_proportions)
  res$monotonic_proportions <- object$monotonic_interval_proportions

  # Calculate Ljung-Box test for autocorrelation in residuals
  if (requireNamespace("stats", quietly = TRUE)) {
    # Calculate residuals from original data
    original_indices <- match(1:res$n_points, object$order)
    residuals <- object$y_sorted[original_indices] - object$predictions
    res$lb_test <- try(stats::Box.test(residuals, lag = min(10, length(residuals)-1),
                                       type = "Ljung-Box"), silent = TRUE)
  }

  class(res) <- "summary.amagelo"
  return(res)
}

#' Print method for summary.amagelo objects
#'
#' @param x An object of class 'summary.amagelo'
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments passed to print
#' @export
print.summary.amagelo <- function(x, digits = 4, ...) {
  cat("AMAGELO: Adaptive Model Averaged GEodesic LOcal linear smoothing\n")
  cat("==============================================================\n\n")

  # Basic information
  cat("Data Summary:\n")
  cat(sprintf("  Number of data points: %d\n", x$n_points))
  cat(sprintf("  Grid size: %d\n", x$grid_size))
  cat(sprintf("  X range: [%.4g, %.4g] (width: %.4g)\n",
              x$x_range[1], x$x_range[2], diff(x$x_range)))
  cat(sprintf("  Y range: [%.4g, %.4g] (width: %.4g)\n",
              x$y_range[1], x$y_range[2], diff(x$y_range)))

  # Bandwidth information
  cat("\nBandwidth Information:\n")
  cat(sprintf("  Optimal bandwidth: %.4g\n", x$opt_bw))
  cat(sprintf("  Bandwidth range: [%.4g, %.4g]\n", x$bw_range[1], x$bw_range[2]))
  cat(sprintf("  Cross-validation error at optimal bandwidth: %.4g\n", x$min_error))

  # Prediction statistics
  cat("\nPrediction Statistics:\n")
  print(round(x$pred_stats, digits))

  # Extrema information
  cat("\nLocal Extrema Summary:\n")
  cat(sprintf("  Number of local maxima: %d\n", x$n_maxima))
  cat(sprintf("  Number of local minima: %d\n", x$n_minima))

  # Top extrema by prominence
  cat("\nTop Extrema by Prominence:\n")
  type_labels <- ifelse(x$top_extrema[,"is_max"] == 1, "Maximum", "Minimum")
  extrema_df <- data.frame(
    Type = type_labels,
    x = round(x$top_extrema[,"x"], digits),
    y = round(x$top_extrema[,"y"], digits),
    RelDepth = round(x$top_extrema[,"range_rel_depth"], digits)
  )
  print(extrema_df, row.names = FALSE)

  # Monotonicity measures
  cat("\nMonotonicity Measures:\n")
  cat(sprintf("  TVMI (Total-Variation Monotonicity Index): %.4g\n", x$tvmi))
  cat(sprintf("  Simpson index: %.4g\n", x$simpson_index))
  cat(sprintf("  Number of monotonic intervals: %d\n", x$monotonic_intervals))

  if (length(x$monotonic_proportions) <= 10) {
    cat("  Monotonic interval proportions: ")
    cat(paste(round(x$monotonic_proportions, digits), collapse = ", "))
    cat("\n")
  } else {
    cat("  Monotonic interval proportions: (first 5) ")
    cat(paste(round(x$monotonic_proportions[1:5], digits), collapse = ", "))
    cat(" ...\n")
  }

  # Ljung-Box test results
  if (!inherits(x$lb_test, "try-error") && !is.null(x$lb_test)) {
    cat("\nLjung-Box test for autocorrelation in residuals:\n")
    cat(sprintf("  Chi-squared = %.4g, df = %d, p-value = %.4g\n",
                x$lb_test$statistic, x$lb_test$parameter, x$lb_test$p.value))
    cat(sprintf("  %s\n", ifelse(x$lb_test$p.value < 0.05,
                              "Residuals may have remaining autocorrelation",
                              "No significant autocorrelation in residuals")))
  }

  invisible(x)
}

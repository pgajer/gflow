#' Local Regression on Graphs Using Spectral Embedding
#'
#' @description Performs local regression on graph data using spectral embeddings
#' with adaptive bandwidth selection.
#'
#' @details This function implements a graph-based extension of LOWESS (Locally
#' Weighted Scatterplot Smoothing) that uses spectral embedding to transform
#' graph distances into a Euclidean space suitable for local linear regression.
#' For each vertex, the function:
#' \enumerate{
#'   \item Finds all vertices within the maximum bandwidth radius
#'   \item Creates a local spectral embedding using graph Laplacian eigenvectors
#'   \item Fits weighted linear models at multiple candidate bandwidths
#'   \item Selects the optimal bandwidth based on leave-one-out cross-validation error
#'   \item Computes smoothed predictions using the optimal model
#' }
#'
#' @param adj.list A list of integer vectors representing the adjacency list of the graph.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of numeric vectors with edge weights corresponding to adjacencies.
#'   Each element \code{weight.list[[i]][j]} is the weight of the edge from vertex i to
#'   \code{adj.list[[i]][j]}.
#' @param y A numeric vector of response values for each vertex in the graph.
#' @param n.evectors Integer specifying the number of eigenvectors to use in the spectral
#'   embedding (default: 5).
#' @param n.bws Integer specifying the number of candidate bandwidths to evaluate (default: 10).
#' @param log.grid Logical indicating whether to use logarithmic spacing for bandwidth
#'   grid (default: TRUE).
#' @param min.bw.factor Numeric value specifying the minimum bandwidth as a fraction of
#'   graph diameter (default: 0.05).
#' @param max.bw.factor Numeric value specifying the maximum bandwidth as a fraction of
#'   graph diameter (default: 0.25).
#' @param dist.normalization.factor Numeric factor for normalizing distances when calculating
#'   kernel weights (default: 1.0).
#' @param kernel.type Integer specifying the kernel function for weighting vertices:
#'        \itemize{
#'          \item 1: Epanechnikov
#'          \item 2: Triangular
#'          \item 4: Laplace
#'          \item 5: Normal
#'          \item 6: Biweight
#'          \item 7: Tricube (default)
#'        }
#'        Default is 7.
#' @param precision Numeric value specifying the precision tolerance for binary search and
#'   optimization algorithms (default: 0.001).
#' @param n.cleveland.iterations Number of Cleveland's robustness iterations (default: 0)
#' @param verbose Logical indicating whether to display progress information (default: FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{predictions}: Numeric vector of smoothed values for each vertex
#'   \item \code{errors}: Numeric vector of leave-one-out cross-validation errors
#'   \item \code{scale}: Numeric vector of optimal bandwidths (local scales) for each vertex
#'   \item \code{graph.diameter}: Numeric scalar with computed graph diameter
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 100 vertices
#' n <- 100
#' set.seed(123)
#'
#' # Create a ring graph
#' adj.list <- vector("list", n)
#' weight.list <- vector("list", n)
#'
#' for (i in 1:n) {
#'   neighbors <- c(i-1, i+1)
#'   # Handle wrap-around for ring structure
#'   neighbors[neighbors == 0] <- n
#'   neighbors[neighbors == n+1] <- 1
#'
#'   adj.list[[i]] <- neighbors
#'   weight.list[[i]] <- rep(1, length(neighbors))
#' }
#'
#' # Generate response values with spatial pattern
#' y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
#'
#' # Apply spectral LOWESS
#' result <- graph.spectral.lowess(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   n.evectors = 5,
#'   verbose = TRUE
#' )
#'
#' # Plot results
#' plot(y, type="l", col="gray", main="Graph Spectral LOWESS")
#' lines(result$predictions, col="red", lwd=2)
#' }
#'
#' @export
graph.spectral.lowess <- function(adj.list,
                                  weight.list,
                                  y,
                                  n.evectors = 5,
                                  n.bws = 20,
                                  log.grid = TRUE,
                                  min.bw.factor = 0.05,
                                  max.bw.factor = 0.5,
                                  dist.normalization.factor = 1.1,
                                  kernel.type = 7L,
                                  precision = 0.001,
                                  n.cleveland.iterations = 1L,
                                  verbose = FALSE) {

    ## Check for connected components (assuming graph.connected.components exists)
    component_ids <- graph.connected.components(adj.list)
    n.ccomp <- length(unique(component_ids))
    if (n.ccomp > 1) {
        stop("The graph has ", n.ccomp, " connected components. ",
             "This function requires a connected graph (single component).")
    }

    ## Basic parameter validation
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices in the graph")
    }

    ## for (i in seq_along(adj.list)) {
    ##   if (length(adj.list[[i]]) != length(weight.list[[i]])) {
    ##     stop(paste0("Adjacency and weight lists have different lengths at vertex ", i))
    ##   }
    ## }

    if (n.evectors < 2) {
        warning("n.evectors should be at least 2; setting to 2")
        n.evectors <- 2
    }

    if (n.bws < 1)
        stop("n.bws must be positive")

    if (!is.numeric(min.bw.factor))
        stop("min.bw.factor must be numeric")

    if (!is.numeric(max.bw.factor) || max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (dist.normalization.factor < 1.05)
        stop("dist.normalization.factor must be greater than or equal to 1.05")

    kernel.type <- as.integer(kernel.type)
    if (!kernel.type %in% c(1L, 2L, 4L, 5L, 6L, 7L)) {
        stop("'kernel.type' must be one of: 1 (Epanechnikov), 2 (Triangular),
             4 (Laplace), 5 (Normal), 6 (Biweight), 7 (Tricube)")
    }

    if (precision <= 0)
        stop("precision must be positive")

    if(!n.cleveland.iterations %in% 0:10) {
        stop("'n.cleveland.iterations' must be an integer between 0 and 10")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ implementation
    result <- .Call("S_graph_spectral_lowess",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.integer(n.evectors),
                    as.integer(n.bws),
                    as.logical(log.grid),
                    as.numeric(min.bw.factor),
                    as.numeric(max.bw.factor),
                    as.numeric(dist.normalization.factor),
                    as.integer(kernel.type),
                    as.numeric(precision),
                    as.integer(n.cleveland.iterations),
                    as.logical(verbose))

    return(result)
}

#' Utility functions for graph.spectral.lowess
#'
#' @description
#' This file contains utility functions for working with graph.spectral.lowess objects,
#' including summary statistics, predictions, residuals, and other methods.

#' Summary Statistics for Graph Spectral LOWESS
#'
#' @description
#' Provides comprehensive summary statistics for graph spectral LOWESS fits,
#' including model information, fit statistics, residuals analysis, and
#' bandwidth selection details.
#'
#' @param object A 'graph.spectral.lowess' object
#' @param quantiles Numeric vector of quantiles for residual statistics (default: c(0, 0.25, 0.5, 0.75, 1))
#' @param ... Additional arguments (currently unused)
#'
#' @return A 'summary.graph.spectral.lowess' object containing:
#' \itemize{
#'   \item model_info: Basic information about the model fit
#'   \item fit_stats: Prediction accuracy metrics
#'   \item bandwidth_stats: Bandwidth selection statistics
#'   \item residual_stats: Residual analysis results
#'   \item component_info: Graph component information
#' }
#'
#' @export
summary.graph.spectral.lowess <- function(object, quantiles = c(0, 0.25, 0.5, 0.75, 1), ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    if (!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
        stop("quantiles must be numeric values between 0 and 1")
    }

    # Calculate residuals
    residuals <- object$y - object$predictions

    # Basic model information
    n_vertices <- length(object$predictions)
    model_info <- list(
        n_vertices = n_vertices,
        n_edges = sum(lengths(object$adj.list)) / 2,  # Assuming undirected
        n_evectors = object$n.evectors,
        graph_diameter = object$graph.diameter,
        kernel_type = object$kernel.type,
        n_bws = object$n.bws
    )

    # Fit statistics
    fit_stats <- list(
        mse = mean(residuals^2),
        rmse = sqrt(mean(residuals^2)),
        mae = mean(abs(residuals)),
        median_ae = median(abs(residuals)),
        r_squared = 1 - (sum(residuals^2) / sum((object$y - mean(object$y))^2))
    )

    # Bandwidth statistics
    bandwidth_stats <- list(
        mean_bandwidth = mean(object$scale),
        sd_bandwidth = sd(object$scale),
        min_bandwidth = min(object$scale),
        max_bandwidth = max(object$scale),
        bandwidth_quantiles = quantile(object$scale, probs = quantiles)
    )

    # Residual statistics
    residual_stats <- list(
        mean = mean(residuals),
        sd = sd(residuals),
        quantiles = quantile(residuals, probs = quantiles),
        shapiro_test = if (n_vertices <= 5000) shapiro.test(residuals) else NULL
    )

    # Component information
    component_info <- list(
        n_components = object$n.components,
        component_sizes = object$component.sizes,
        isolated_vertices = sum(lengths(object$adj.list) == 0)
    )

    # Create summary object
    result <- list(
        call = object$call,
        model_info = model_info,
        fit_stats = fit_stats,
        bandwidth_stats = bandwidth_stats,
        residual_stats = residual_stats,
        component_info = component_info
    )

    class(result) <- "summary.graph.spectral.lowess"
    return(result)
}

#' Print Summary Statistics for Graph Spectral LOWESS
#'
#' @description
#' Formats and displays comprehensive summary statistics for graph spectral LOWESS fits.
#'
#' @param x A 'summary.graph.spectral.lowess' object
#' @param digits Number of significant digits for numerical output (default: 4)
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.summary.graph.spectral.lowess <- function(x, digits = 4, ...) {
    ## Helper function for separator lines
    hr <- function() cat(paste0(rep("\u2500", 80), collapse = ""), "\n")

    # Helper function for section headers
    section_header <- function(text) {
        cat("\n")
        hr()
        cat(text, "\n")
        hr()
    }

    # Main title
    cat("\n")
    section_header("GRAPH SPECTRAL LOWESS SUMMARY")

    # Model Information
    cat("\nModel Information:\n")
    cat(sprintf("Number of vertices:             %d\n", x$model_info$n_vertices))
    cat(sprintf("Number of edges:                %d\n", x$model_info$n_edges))
    cat(sprintf("Number of eigenvectors:         %d\n", x$model_info$n_evectors))
    cat(sprintf("Graph diameter:                 %.3f\n", x$model_info$graph_diameter))
    cat(sprintf("Number of bandwidths tested:    %d\n", x$model_info$n_bws))

    # Component Information
    if (!is.null(x$component_info)) {
        section_header("GRAPH STRUCTURE")
        cat(sprintf("Connected components:           %d\n", x$component_info$n_components))
        if (x$component_info$n_components > 1) {
            cat("Component sizes:                ")
            cat(paste(x$component_info$component_sizes, collapse = ", "), "\n")
        }
        if (x$component_info$isolated_vertices > 0) {
            cat(sprintf("Isolated vertices:              %d\n", x$component_info$isolated_vertices))
        }
    }

    # Fit Statistics
    section_header("FIT STATISTICS")
    cat(sprintf("MSE:                            %.4f\n", x$fit_stats$mse))
    cat(sprintf("RMSE:                           %.4f\n", x$fit_stats$rmse))
    cat(sprintf("MAE:                            %.4f\n", x$fit_stats$mae))
    cat(sprintf("Median AE:                      %.4f\n", x$fit_stats$median_ae))
    cat(sprintf("R-squared:                      %.4f\n", x$fit_stats$r_squared))

    # Bandwidth Statistics
    section_header("BANDWIDTH STATISTICS")
    cat(sprintf("Mean bandwidth:                 %.4f\n", x$bandwidth_stats$mean_bandwidth))
    cat(sprintf("SD bandwidth:                   %.4f\n", x$bandwidth_stats$sd_bandwidth))
    cat(sprintf("Range:                          [%.4f, %.4f]\n",
                x$bandwidth_stats$min_bandwidth, x$bandwidth_stats$max_bandwidth))
    cat("\nBandwidth quantiles:\n")
    print(format(data.frame(
        Quantile = names(x$bandwidth_stats$bandwidth_quantiles),
        Value = x$bandwidth_stats$bandwidth_quantiles
    ), digits = digits), row.names = FALSE)

    # Residual Statistics
    section_header("RESIDUAL STATISTICS")
    cat(sprintf("Mean:                           %.4f\n", x$residual_stats$mean))
    cat(sprintf("Standard Deviation:             %.4f\n", x$residual_stats$sd))
    if (!is.null(x$residual_stats$shapiro_test)) {
        cat(sprintf("Shapiro-Wilk p-value:          %.4f\n",
                    x$residual_stats$shapiro_test$p.value))
    }
    cat("\nResidual quantiles:\n")
    print(format(data.frame(
        Quantile = names(x$residual_stats$quantiles),
        Value = x$residual_stats$quantiles
    ), digits = digits), row.names = FALSE)

    cat("\n")  # Final newline for spacing
    invisible(x)
}

#' Print Graph Spectral LOWESS Object
#'
#' @description
#' Provides a concise summary of a graph spectral LOWESS object
#'
#' @param x A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns x invisibly
#'
#' @export
print.graph.spectral.lowess <- function(x, ...) {
    cat("Graph Spectral LOWESS\n")
    cat("=====================\n")
    cat(sprintf("Number of vertices:    %d\n", length(x$predictions)))
    cat(sprintf("Graph diameter:        %.3f\n", x$graph.diameter))
    cat(sprintf("Eigenvectors used:     %d\n", x$n.evectors))
    cat(sprintf("RMSE:                  %.4f\n", sqrt(mean(x$errors^2))))
    cat(sprintf("Mean bandwidth:        %.4f\n", mean(x$scale)))
    invisible(x)
}

#' Extract Fitted Values from Graph Spectral LOWESS
#'
#' @description
#' Extracts the fitted (smoothed) values from a graph spectral LOWESS object
#'
#' @param object A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector of fitted values
#'
#' @export
fitted.graph.spectral.lowess <- function(object, ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }
    return(object$predictions)
}

#' Extract Residuals from Graph Spectral LOWESS
#'
#' @description
#' Computes and returns residuals from a graph spectral LOWESS fit
#'
#' @param object A 'graph.spectral.lowess' object
#' @param type Character string specifying residual type: "response" (default), "pearson", or "deviance"
#' @param ... Additional arguments (currently unused)
#'
#' @return Numeric vector of residuals
#'
#' @export
residuals.graph.spectral.lowess <- function(object, type = c("response", "pearson", "deviance"), ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    type <- match.arg(type)

    # Basic residuals
    res <- object$y - object$predictions

    switch(type,
           "response" = res,
           "pearson" = {
               # Pearson residuals: residual / sqrt(variance)
               # Estimate variance from local fits
               local_var <- pmax(object$errors^2, .Machine$double.eps)
               res / sqrt(local_var)
           },
           "deviance" = {
               # For normal errors, deviance residuals = pearson residuals
               local_var <- pmax(object$errors^2, .Machine$double.eps)
               sign(res) * sqrt(abs(res)^2 / local_var)
           }
    )
}

#' Extract Coefficients from Graph Spectral LOWESS
#'
#' @description
#' Returns the smoothed values (predictions) as coefficients.
#' Note: Unlike linear models, LOWESS doesn't have traditional coefficients.
#'
#' @param object A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Named numeric vector with smoothed values
#'
#' @export
coef.graph.spectral.lowess <- function(object, ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    coefs <- object$predictions
    names(coefs) <- paste0("vertex_", seq_along(coefs))
    return(coefs)
}

#' Predict Method for Graph Spectral LOWESS
#'
#' @description
#' Makes predictions from a graph spectral LOWESS fit.
#' For new data, this would require extending the graph structure.
#'
#' @param object A 'graph.spectral.lowess' object
#' @param newdata Optional new data (currently not implemented)
#' @param se.fit Logical; should standard errors be returned? (currently not implemented)
#' @param ... Additional arguments (currently unused)
#'
#' @return If newdata is missing, returns fitted values. Otherwise, not yet implemented.
#'
#' @export
predict.graph.spectral.lowess <- function(object, newdata, se.fit = FALSE, ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    if (missing(newdata)) {
        # Return fitted values for original data
        return(object$predictions)
    } else {
        stop("Prediction for new data not yet implemented for graph spectral LOWESS")
    }
}

#' Confidence Intervals for Graph Spectral LOWESS
#'
#' @description
#' Computes confidence intervals for the smoothed values.
#' Note: This requires bootstrap or other uncertainty quantification methods.
#'
#' @param object A 'graph.spectral.lowess' object
#' @param parm Parameter specification (currently unused)
#' @param level Confidence level (default: 0.95)
#' @param ... Additional arguments (currently unused)
#'
#' @return Matrix with lower and upper confidence bounds
#'
#' @export
confint.graph.spectral.lowess <- function(object, parm, level = 0.95, ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    # Check if confidence intervals are available in the object
    if (!is.null(object$ci.lower) && !is.null(object$ci.upper)) {
        ci <- cbind(lower = object$ci.lower, upper = object$ci.upper)
        rownames(ci) <- paste0("vertex_", seq_len(nrow(ci)))
        return(ci)
    } else {
        warning("Confidence intervals not available. Consider using bootstrap.")
        return(NULL)
    }
}

#' Extract Log-Likelihood from Graph Spectral LOWESS
#'
#' @description
#' Computes the log-likelihood assuming normal errors
#'
#' @param object A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Log-likelihood value
#'
#' @export
logLik.graph.spectral.lowess <- function(object, ...) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    n <- length(object$y)
    residuals <- object$y - object$predictions
    sigma2 <- mean(residuals^2)

    # Log-likelihood for normal distribution
    ll <- -n/2 * log(2 * pi) - n/2 * log(sigma2) - sum(residuals^2) / (2 * sigma2)

    # Add attributes
    attr(ll, "df") <- mean(object$scale) * n  # Approximate degrees of freedom
    attr(ll, "nobs") <- n
    class(ll) <- "logLik"

    return(ll)
}

#' AIC for Graph Spectral LOWESS
#'
#' @description
#' Computes the Akaike Information Criterion
#'
#' @param object A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#' @param k Penalty parameter (default: 2)
#'
#' @return AIC value
#'
#' @export
AIC.graph.spectral.lowess <- function(object, ..., k = 2) {
    ll <- logLik(object)
    df <- attr(ll, "df")
    return(-2 * as.numeric(ll) + k * df)
}

#' BIC for Graph Spectral LOWESS
#'
#' @description
#' Computes the Bayesian Information Criterion
#'
#' @param object A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return BIC value
#'
#' @export
BIC.graph.spectral.lowess <- function(object, ...) {
    ll <- logLik(object)
    df <- attr(ll, "df")
    n <- attr(ll, "nobs")
    return(-2 * as.numeric(ll) + log(n) * df)
}

#' Extract Model Frame from Graph Spectral LOWESS
#'
#' @description
#' Constructs a data frame containing the graph structure and response values
#'
#' @param formula A 'graph.spectral.lowess' object
#' @param ... Additional arguments (currently unused)
#'
#' @return Data frame with vertex information
#'
#' @export
model.frame.graph.spectral.lowess <- function(formula, ...) {
    object <- formula  # R convention uses 'formula' as first argument

    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    # Create a data frame with vertex information
    df <- data.frame(
        vertex = seq_along(object$y),
        y = object$y,
        fitted = object$predictions,
        residuals = object$y - object$predictions,
        bandwidth = object$scale,
        degree = lengths(object$adj.list)
    )

    return(df)
}

#' Update Graph Spectral LOWESS Fit
#'
#' @description
#' Updates a graph spectral LOWESS fit with new parameters
#'
#' @param object A 'graph.spectral.lowess' object
#' @param formula Changes to the model (currently unused)
#' @param ... Additional arguments to override in the new call
#' @param evaluate Should the updated call be evaluated? (default: TRUE)
#'
#' @return Updated graph.spectral.lowess object (if evaluate = TRUE) or updated call
#'
#' @export
update.graph.spectral.lowess <- function(object, formula, ..., evaluate = TRUE) {
    if (!inherits(object, "graph.spectral.lowess")) {
        stop("Input must be a 'graph.spectral.lowess' object")
    }

    # Get original call
    call <- object$call
    if (is.null(call)) {
        stop("Need an object with call component")
    }

    # Update call with new arguments
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    if (evaluate) {
        eval(call, parent.frame())
    } else {
        call
    }
}

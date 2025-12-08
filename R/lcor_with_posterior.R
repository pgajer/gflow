## ============================================================================
## lslope with Posterior Uncertainty Propagation
## ============================================================================
##
## This function computes lslope() statistics with full posterior uncertainty
## quantification arising from spectral smoothing. Two modes are supported:
##
## Mode 1 (Original): Pass pre-computed posterior samples in Z.hat.samples
##   - Useful when you already have samples from refit.rdgraph.regression()
##   - Memory cost: O(n × B × p) for storing all samples
##
## Mode 2 (Memory-efficient): Pass fitted.model and Z.abundances with Z.hat.samples = NULL
##   - Computes posterior samples and lslope on-the-fly in C++
##   - Memory cost: O(n × B) per feature, O(n × p) total output
##   - Supports OpenMP parallelization over features
##
## ============================================================================

#' Compute Local Correlation with Posterior Uncertainty Propagation
#'
#' Propagates posterior uncertainty from spectral smoothing into local
#' correlation (lcor) estimates, providing posterior mean, standard deviation,
#' and credible intervals for lcor at each vertex for each feature.
#'
#' @param adj.list Adjacency list (1-based R indexing)
#' @param weight.list Edge weight list
#' @param y.hat Smoothed response values (length n)
#' @param Z.hat.samples Either:
#'   \itemize{
#'     \item A list of length p, where each element is an n x B matrix of
#'       posterior samples for one feature (from \code{refit.rdgraph.regression()}
#'       with \code{return.posterior.samples = TRUE}), OR
#'     \item \code{NULL} to use memory-efficient C++ computation (requires
#'       \code{fitted.model} and \code{Z.abundances})
#'   }
#' @param fitted.model Fitted model from \code{fit.rdgraph.regression()}.
#'   Required when \code{Z.hat.samples = NULL}. Must contain spectral
#'   decomposition (\code{spectral$vectors}, \code{spectral$values}).
#' @param Z.abundances Original (unsmoothed) feature matrix (n x p).
#'   Required when \code{Z.hat.samples = NULL}.
#' @param lcor.type Type of local correlation weighting: "derivative" (default),
#'   "unit", or "sign". See \code{\link{lcor}} for details.
#' @param per.column.gcv Logical. When using C++ mode, select optimal eta
#'   for each feature via GCV (default TRUE). Ignored when Z.hat.samples
#'   is provided.
#' @param n.posterior.samples Number of posterior samples when using C++ mode
#'   (default 500). Ignored when Z.hat.samples is provided.
#' @param credible.level Credible interval level (default 0.95)
#' @param seed Random seed for posterior sampling (default 12345)
#' @param n.cores Number of cores for parallel processing (default 1).
#'   Only used in C++ mode. Requires OpenMP support.
#' @param return.samples Logical. Return individual lcor samples (default FALSE).
#'   Only available in R mode (when Z.hat.samples is provided). Setting TRUE
#'   with many features will use substantial memory.
#' @param verbose Logical. Print progress (default TRUE)
#'
#' @return A list of class "lcor.posterior" containing:
#'   \describe{
#'     \item{mean}{Matrix (p x n) of posterior mean lcor values}
#'     \item{sd}{Matrix (p x n) of posterior standard deviations}
#'     \item{lower}{Matrix (p x n) of lower credible bounds}
#'     \item{upper}{Matrix (p x n) of upper credible bounds}
#'     \item{eta.used}{Vector (length p) of eta values used (C++ mode only)}
#'     \item{effective.df}{Vector (length p) of effective df (C++ mode only)}
#'     \item{samples}{List of sample matrices (only if return.samples = TRUE)}
#'     \item{n.samples}{Number of posterior samples used}
#'     \item{n.features}{Number of features}
#'     \item{n.vertices}{Number of vertices}
#'     \item{lcor.type}{Type of lcor weighting used}
#'     \item{mode}{Either "R" or "C++" indicating computation mode}
#'   }
#'
#' @details
#' This function supports two computational modes:
#'
#' \strong{R Mode} (Z.hat.samples provided): Uses pre-computed posterior samples
#' from \code{refit.rdgraph.regression()}. For each posterior sample of the
#' smoothed features, lcor is computed against the fixed smoothed response.
#' This mode is useful when you want fine control over the smoothing parameters
#' or need to reuse the posterior samples for other analyses.
#'
#' \strong{C++ Mode} (Z.hat.samples = NULL): Computes everything internally in
#' C++ without materializing all posterior samples in memory. For each feature:
#' \enumerate{
#'   \item Select optimal eta via GCV (if per.column.gcv = TRUE)
#'   \item Generate B posterior samples of smoothed feature values
#'   \item Compute lcor for each sample
#'   \item Compute summary statistics (mean, sd, quantiles)
#'   \item Discard samples before processing next feature
#' }
#' This mode is recommended when memory is limited or when processing many
#' features. OpenMP parallelization over features is supported.
#'
#' The posterior sampling model assumes:
#' \deqn{\hat{z} = V F_\eta(\Lambda) V^T z + \epsilon}
#' where the spectral coefficients have posterior distribution:
#' \deqn{\alpha_j | z \sim N(f(\lambda_j) (V^T z)_j, \sigma^2 / (1 + \eta \lambda_j))}
#'
#' @examples
#' \dontrun{
#' ## Fit model
#' fit <- fit.rdgraph.regression(X, y, k = 15)
#'
#' ## ----- Mode 1: Using pre-computed samples -----
#' Z.refit <- refit.rdgraph.regression(
#'     fit, Z[, 1:10],
#'     per.column.gcv = TRUE,
#'     with.posterior = TRUE,
#'     return.posterior.samples = TRUE,
#'     n.posterior.samples = 500
#' )
#'
#' lcor.post <- lcor.with.posterior(
#'     fit$graph$adj.list,
#'     fit$graph$edge.length.list,
#'     fit$fitted.values,
#'     Z.refit$posterior$samples
#' )
#'
#' ## ----- Mode 2: Memory-efficient C++ computation -----
#' lcor.post <- lcor.with.posterior(
#'     fit$graph$adj.list,
#'     fit$graph$edge.length.list,
#'     fit$fitted.values,
#'     Z.hat.samples = NULL,       # Triggers C++ mode
#'     fitted.model = fit,
#'     Z.abundances = Z,
#'     per.column.gcv = TRUE,
#'     n.posterior.samples = 500,
#'     n.cores = 4
#' )
#' }
#'
#' @seealso \code{\link{lcor}} for single-sample local correlation computation
#'
#' @export
lcor.with.posterior <- function(adj.list,
                                 weight.list,
                                 y.hat,
                                 Z.hat.samples = NULL,
                                 fitted.model = NULL,
                                 Z.abundances = NULL,
                                 lcor.type = c("derivative", "unit", "sign"),
                                 per.column.gcv = TRUE,
                                 n.posterior.samples = 500L,
                                 credible.level = 0.95,
                                 seed = 12345L,
                                 n.cores = 1L,
                                 return.samples = FALSE,
                                 verbose = TRUE) {

    lcor.type <- match.arg(lcor.type)
    n <- length(y.hat)

    ## ========================================================================
    ## Determine computation mode
    ## ========================================================================

    use.cpp.mode <- is.null(Z.hat.samples)

    if (use.cpp.mode) {
        ## Validate required arguments for C++ mode
        if (is.null(fitted.model)) {
            stop("fitted.model is required when Z.hat.samples = NULL")
        }
        if (is.null(Z.abundances)) {
            stop("Z.abundances is required when Z.hat.samples = NULL")
        }
        if (!inherits(fitted.model, "knn.riem.fit")) {
            stop("fitted.model must be a 'knn.riem.fit' object")
        }
        if (is.null(fitted.model$spectral)) {
            stop("fitted.model must contain spectral decomposition")
        }

        if (is.vector(Z.abundances)) {
            Z.abundances <- matrix(Z.abundances, ncol = 1)
        }
        if (nrow(Z.abundances) != n) {
            stop("nrow(Z.abundances) must equal length(y.hat)")
        }

        result <- lcor.with.posterior.cpp(
            adj.list = adj.list,
            weight.list = weight.list,
            y.hat = y.hat,
            Z.abundances = Z.abundances,
            fitted.model = fitted.model,
            lcor.type = lcor.type,
            per.column.gcv = per.column.gcv,
            n.posterior.samples = n.posterior.samples,
            credible.level = credible.level,
            seed = seed,
            n.cores = n.cores,
            verbose = verbose
        )

    } else {
        ## R mode: use pre-computed samples
        result <- lcor.with.posterior.R(
            adj.list = adj.list,
            weight.list = weight.list,
            y.hat = y.hat,
            Z.hat.samples = Z.hat.samples,
            lcor.type = lcor.type,
            credible.level = credible.level,
            return.samples = return.samples,
            verbose = verbose
        )
    }

    return(result)
}


## ============================================================================
## R Mode Implementation
## ============================================================================

#' @keywords internal
lcor.with.posterior.R <- function(adj.list,
                                   weight.list,
                                   y.hat,
                                   Z.hat.samples,
                                   lcor.type,
                                   credible.level,
                                   return.samples,
                                   verbose) {

    n <- length(y.hat)

    ## Handle single feature case
    if (is.matrix(Z.hat.samples)) {
        Z.hat.samples <- list(Z.hat.samples)
    }

    ## Validate Z.hat.samples
    if (length(Z.hat.samples) == 0) {
        stop("Z.hat.samples cannot be empty")
    }

    p <- length(Z.hat.samples)

    ## Validate first element has correct dimensions
    if (!is.matrix(Z.hat.samples[[1]])) {
        stop("Each element of Z.hat.samples must be a matrix")
    }
    if (nrow(Z.hat.samples[[1]]) != n) {
        stop(sprintf("nrow(Z.hat.samples[[1]]) = %d must equal length(y.hat) = %d",
                     nrow(Z.hat.samples[[1]]), n))
    }

    n.samples <- ncol(Z.hat.samples[[1]])

    if (verbose) {
        message(sprintf("lcor with posterior (R mode): %d features, %d samples",
                        p, n.samples))
    }

    ## Initialize output storage
    lcor.mean <- matrix(NA_real_, nrow = p, ncol = n)
    lcor.sd <- matrix(NA_real_, nrow = p, ncol = n)
    lcor.lower <- matrix(NA_real_, nrow = p, ncol = n)
    lcor.upper <- matrix(NA_real_, nrow = p, ncol = n)

    if (return.samples) {
        all.samples <- vector("list", p)
    }

    alpha.lower <- (1 - credible.level) / 2
    alpha.upper <- (1 + credible.level) / 2

    if (verbose && p > 5) {
        pb <- txtProgressBar(min = 0, max = p, style = 3)
    }

    for (j in seq_len(p)) {
        ## Z.hat.samples[[j]] is n x n.samples matrix
        samples.j <- Z.hat.samples[[j]]

        ## Compute lcor for each posterior sample
        lcor.samples.j <- matrix(NA_real_, nrow = n, ncol = n.samples)

        for (b in seq_len(n.samples)) {
            z.hat.b <- samples.j[, b]

            ## lcor() with instrumented=FALSE returns a vector directly
            lcor.coeffs <- lcor(
                adj.list, weight.list,
                y.hat, z.hat.b,
                type = lcor.type,
                y.diff.type = "difference",
                z.diff.type = "difference",
                epsilon = 0,
                winsorize.quantile = 0,
                instrumented = FALSE
            )

            lcor.samples.j[, b] <- lcor.coeffs
        }

        ## Compute summary statistics
        lcor.mean[j, ] <- rowMeans(lcor.samples.j, na.rm = TRUE)
        lcor.sd[j, ] <- apply(lcor.samples.j, 1, sd, na.rm = TRUE)
        lcor.lower[j, ] <- apply(lcor.samples.j, 1, quantile,
                                  probs = alpha.lower, na.rm = TRUE)
        lcor.upper[j, ] <- apply(lcor.samples.j, 1, quantile,
                                  probs = alpha.upper, na.rm = TRUE)

        if (return.samples) {
            all.samples[[j]] <- lcor.samples.j
        }

        if (verbose && p > 5) {
            setTxtProgressBar(pb, j)
        }
    }

    if (verbose && p > 5) {
        close(pb)
    }

    ## Build result
    result <- list(
        mean = lcor.mean,
        sd = lcor.sd,
        lower = lcor.lower,
        upper = lcor.upper,
        n.samples = n.samples,
        n.features = p,
        n.vertices = n,
        lcor.type = lcor.type,
        credible.level = credible.level,
        mode = "R"
    )

    if (return.samples) {
        result$samples <- all.samples
    }

    class(result) <- c("lcor.posterior", "list")
    return(result)
}


## ============================================================================
## C++ Mode Implementation (Memory-Efficient)
## ============================================================================

#' @keywords internal
lcor.with.posterior.cpp <- function(adj.list,
                                     weight.list,
                                     y.hat,
                                     Z.abundances,
                                     fitted.model,
                                     lcor.type,
                                     per.column.gcv,
                                     n.posterior.samples,
                                     credible.level,
                                     seed,
                                     n.cores,
                                     verbose) {

    n <- length(y.hat)
    p <- ncol(Z.abundances)

    if (verbose) {
        message(sprintf("lcor with posterior (C++ mode): %d features, %d samples",
                        p, n.posterior.samples))
        message(sprintf("  Memory-efficient: processing features sequentially"))
        if (n.cores > 1) {
            message(sprintf("  Parallel: %d cores", n.cores))
        }
    }

    ## Extract spectral components
    V <- fitted.model$spectral$vectors
    eigenvalues <- fitted.model$spectral$values
    filter.type <- fitted.model$spectral$filter.type
    if (is.null(filter.type)) filter.type <- "heat_kernel"

    ## Get eta (for fixed eta mode)
    eta.fixed <- fitted.model$spectral$eta
    if (is.null(eta.fixed)) eta.fixed <- 1.0

    n.gcv.candidates <- 50L

    ## Convert adjacency list to 0-based indexing for C++
    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1L))

    ## Ensure Z is a matrix with proper storage
    Z.mat <- as.matrix(Z.abundances)
    storage.mode(Z.mat) <- "double"

    ## Call C++ implementation
    cpp.result <- .Call(
        S_lcor_with_posterior_internal,
        adj.list.0,
        weight.list,
        as.double(y.hat),
        Z.mat,
        as.matrix(V),
        as.double(eigenvalues),
        as.double(eta.fixed),
        as.character(lcor.type),
        as.character(filter.type),
        as.logical(per.column.gcv),
        as.integer(n.gcv.candidates),
        as.integer(n.posterior.samples),
        as.double(credible.level),
        as.integer(seed),
        as.integer(n.cores),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## Add row names if available
    if (!is.null(colnames(Z.abundances))) {
        rownames(cpp.result$mean) <- colnames(Z.abundances)
        rownames(cpp.result$sd) <- colnames(Z.abundances)
        rownames(cpp.result$lower) <- colnames(Z.abundances)
        rownames(cpp.result$upper) <- colnames(Z.abundances)
        names(cpp.result$eta.used) <- colnames(Z.abundances)
        names(cpp.result$effective.df) <- colnames(Z.abundances)
    }

    ## Build result object matching R mode output
    result <- list(
        mean = cpp.result$mean,
        sd = cpp.result$sd,
        lower = cpp.result$lower,
        upper = cpp.result$upper,
        eta.used = cpp.result$eta.used,
        effective.df = cpp.result$effective.df,
        n.samples = cpp.result$n.samples,
        n.features = p,
        n.vertices = n,
        lcor.type = lcor.type,
        credible.level = cpp.result$credible.level,
        mode = "C++"
    )

    class(result) <- c("lcor.posterior", "list")
    return(result)
}


#' @export
print.lcor.posterior <- function(x, ...) {
    cat("\nLocal Correlation with Posterior Uncertainty Propagation\n")
    cat("========================================================\n\n")

    cat(sprintf("Mode: %s\n", x$mode))
    cat(sprintf("Features: %d\n", x$n.features))
    cat(sprintf("Vertices: %d\n", x$n.vertices))
    cat(sprintf("Posterior samples: %d\n", x$n.samples))
    cat(sprintf("Credible level: %.0f%%\n", 100 * x$credible.level))
    cat(sprintf("lcor type: %s\n", x$lcor.type))

    cat(sprintf("\nPosterior mean lcor (summary):\n"))
    cat(sprintf("  Global mean: %.4f\n", mean(x$mean, na.rm = TRUE)))
    cat(sprintf("  Global SD: %.4f\n", sd(as.vector(x$mean), na.rm = TRUE)))

    cat(sprintf("\nPosterior SD of lcor (summary):\n"))
    cat(sprintf("  Mean SD: %.4f\n", mean(x$sd, na.rm = TRUE)))
    cat(sprintf("  Median SD: %.4f\n", median(x$sd, na.rm = TRUE)))

    ## Credible interval width
    ci.width <- x$upper - x$lower
    cat(sprintf("\n%.0f%% CI width (summary):\n", 100 * x$credible.level))
    cat(sprintf("  Mean width: %.4f\n", mean(ci.width, na.rm = TRUE)))
    cat(sprintf("  Median width: %.4f\n", median(ci.width, na.rm = TRUE)))

    if (x$mode == "C++" && !is.null(x$eta.used)) {
        cat(sprintf("\nSmoothing parameters (eta):\n"))
        cat(sprintf("  Mean: %.4f\n", mean(x$eta.used, na.rm = TRUE)))
        cat(sprintf("  Range: [%.4f, %.4f]\n",
                    min(x$eta.used, na.rm = TRUE),
                    max(x$eta.used, na.rm = TRUE)))
    }

    invisible(x)
}


#' @export
summary.lcor.posterior <- function(object, threshold = 0.5, ...) {

    ## Identify feature-vertex pairs where CI excludes zero
    ci.excludes.zero <- (object$lower > 0) | (object$upper < 0)

    ## Proportion of posterior mass on same side as mean
    ## (crude measure of "significance")
    ## If CI is entirely positive or negative, this is high
    prop.significant <- mean(ci.excludes.zero, na.rm = TRUE)

    ## Get feature names (or generate them)
    feature.names <- rownames(object$mean)
    if (is.null(feature.names)) {
        feature.names <- paste0("Feature", seq_len(object$n.features))
    }

    ## Per-feature summary
    feature.summary <- data.frame(
        feature = feature.names,
        mean.lcor = rowMeans(object$mean, na.rm = TRUE),
        mean.sd = rowMeans(object$sd, na.rm = TRUE),
        prop.ci.excludes.zero = rowMeans(ci.excludes.zero, na.rm = TRUE),
        max.abs.mean = apply(abs(object$mean), 1, max, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    ## Sort by proportion of significant vertices
    feature.summary <- feature.summary[order(-feature.summary$prop.ci.excludes.zero), ]

    result <- list(
        feature.summary = feature.summary,
        prop.significant.overall = prop.significant,
        n.features = object$n.features,
        n.vertices = object$n.vertices,
        n.samples = object$n.samples,
        credible.level = object$credible.level
    )

    class(result) <- "summary.lcor.posterior"
    return(result)
}


#' @export
print.summary.lcor.posterior <- function(x, n.top = 20, ...) {
    cat("\nSummary: Local Correlation with Posterior Uncertainty\n")
    cat("=====================================================\n\n")

    cat(sprintf("Features: %d, Vertices: %d\n", x$n.features, x$n.vertices))
    cat(sprintf("Posterior samples: %d\n", x$n.samples))
    cat(sprintf("Credible level: %.0f%%\n\n", 100 * x$credible.level))

    cat(sprintf("Overall proportion of (feature, vertex) pairs\n"))
    cat(sprintf("  where %.0f%% CI excludes zero: %.1f%%\n\n",
                100 * x$credible.level,
                100 * x$prop.significant.overall))

    cat("Top features by proportion of vertices with CI excluding zero:\n\n")

    top.features <- head(x$feature.summary, n.top)
    print(top.features, row.names = FALSE, digits = 3)

    invisible(x)
}

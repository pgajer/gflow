#' Refit Riemannian Graph Regression with New Response(s)
#'
#' Apply spectral filtering from a fitted model to new response variable(s).
#' This function uses the cached eigendecomposition from an existing fit,
#' avoiding the computational cost of rebuilding the geometric structure.
#'
#' When \code{per.column.gcv = TRUE}, each column of \code{y.new} receives its
#' own GCV-selected smoothing parameter eta. This is appropriate when columns
#' represent distinct features (e.g., phylotype abundances) with different
#' smoothness characteristics. When \code{FALSE} (default), all columns use the
#' same filtered eigenvalues from the original fit, which is faster but assumes
#' similar smoothness structure across responses.
#'
#' @param fitted.model A fitted model object of class \code{"knn.riem.fit"}
#'   returned by \code{\link{fit.rdgraph.regression}}. Must contain the
#'   \code{spectral} component with cached eigendecomposition.
#'
#' @param y.new New response variable(s). Either a numeric vector of length n,
#'   or a numeric matrix with n rows where each column is a separate response.
#'
#' @param per.column.gcv Logical. If \code{TRUE}, select optimal smoothing
#'   parameter eta independently for each column via GCV. If \code{FALSE}
#'   (default), use the filtered eigenvalues from the original fit for all
#'   columns. Setting \code{TRUE} is recommended when columns represent
#'   features with heterogeneous smoothness (e.g., different phylotypes).
#'
#' @param filter.type Character string specifying the spectral filter type.
#'   Only used when \code{per.column.gcv = TRUE}. Options are:
#'   \describe{
#'     \item{\code{"heat_kernel"}}{exp(-eta * lambda). Default and recommended
#'       for general use.}
#'     \item{\code{"tikhonov"}}{1/(1 + eta * lambda). Gentler attenuation,
#'       better preserves linear trends.}
#'     \item{\code{"cubic_spline"}}{1/(1 + eta * lambda^2). Very smooth results.}
#'     \item{\code{"gaussian"}}{exp(-eta * lambda^2). Aggressive smoothing.}
#'   }
#'   If \code{NULL}, uses the filter type from the original fit.
#'
#' @param n.candidates Integer. Number of candidate eta values for GCV grid
#'   search. Only used when \code{per.column.gcv = TRUE}. Default is 40.
#'
#' @param n.cores Integer. Number of CPU cores for parallel processing.
#'   Default is 1 (sequential). When \code{n.cores > 1} and per-column GCV is
#'   enabled, columns are processed in parallel using \pkg{foreach} and
#'   \pkg{doParallel}. Ignored when \code{per.column.gcv = FALSE}.
#'
#' @param verbose Logical. If \code{TRUE}, print progress information during
#'   per-column GCV selection. Default is \code{FALSE}. Note that progress
#'   reporting is limited in parallel mode.
#'
#' @return A list of class \code{"knn.riem.refit"} containing:
#'   \describe{
#'     \item{\code{fitted.values}}{Smoothed response(s). Vector if single
#'       response, matrix if multiple.}
#'     \item{\code{residuals}}{Residuals (y.new - fitted.values).}
#'     \item{\code{n.responses}}{Number of response columns processed.}
#'     \item{\code{per.column.gcv}}{Logical indicating whether per-column GCV
#'       was used.}
#'     \item{\code{eta.optimal}}{(Only if \code{per.column.gcv = TRUE}) Vector
#'       of optimal eta values, one per column.}
#'     \item{\code{gcv.scores}}{(Only if \code{per.column.gcv = TRUE}) Vector
#'       of minimum GCV scores achieved, one per column.}
#'     \item{\code{effective.df}}{(Only if \code{per.column.gcv = TRUE}) Vector
#'       of effective degrees of freedom (trace of smoother matrix).}
#'     \item{\code{n.cores.used}}{(Only if \code{per.column.gcv = TRUE}) Number
#'       of cores actually used for computation.}
#'     \item{\code{elapsed.time}}{(Only if \code{per.column.gcv = TRUE})
#'       Elapsed wall-clock time in seconds.}
#'   }
#'
#' @details
#' The function implements spectral filtering in the eigenbasis of the Hodge
#' Laplacian:
#' \deqn{\hat{y} = V F_\eta(\Lambda) V^T y}
#' where V contains eigenvectors, Lambda contains eigenvalues, and F_eta is the
#' filter function evaluated at each eigenvalue.
#'
#' For per-column GCV, the optimal eta minimizes the generalized cross-validation
#' score:
#' \deqn{GCV(\eta) = \frac{\|y - \hat{y}_\eta\|^2}{(n - \mathrm{tr}(S_\eta))^2}}
#' where tr(S_eta) = sum of filter weights is the effective degrees of freedom.
#'
#' \strong{Parallel processing:} When \code{n.cores > 1}, the function uses
#' \pkg{foreach} with \pkg{doParallel} backend. Each worker receives a copy of
#' the eigenvector matrix V and filter weights, so memory usage scales with
#' the number of cores. For very large n (>10000), consider limiting n.cores
#' to avoid memory pressure.
#'
#' \strong{Performance considerations:} Per-column GCV adds computational cost
#' proportional to n.candidates times the number of columns. For large matrices
#' (>100 columns), parallel processing with \code{n.cores > 1} is recommended.
#' Typical speedup is near-linear up to 4-8 cores, with diminishing returns
#' beyond that due to memory bandwidth limitations.
#'
#' @examples
#' \dontrun{
#' # Fit initial model on binary outcome
#' fit <- fit.rdgraph.regression(X, y.binary, k = 15)
#'
#' # Smooth phylotype abundances with per-column GCV (parallel)
#' Z.hat <- refit.rdgraph.regression(fit, Z.abundances,
#'                                    per.column.gcv = TRUE,
#'                                    n.cores = 4)
#'
#' # Examine the selected smoothing parameters
#' hist(Z.hat$eta.optimal, main = "Optimal eta by phylotype")
#'
#' # Compare effective degrees of freedom across phylotypes
#' plot(Z.hat$effective.df, Z.hat$gcv.scores,
#'      xlab = "Effective df", ylab = "GCV score")
#'
#' # Check timing
#' cat(sprintf("Processed %d columns in %.1f seconds using %d cores\n",
#'             Z.hat$n.responses, Z.hat$elapsed.time, Z.hat$n.cores.used))
#' }
#'
#' @seealso \code{\link{fit.rdgraph.regression}} for initial model fitting
#'
#' @export
refit.rdgraph.regression <- function(fitted.model,
                                     y.new,
                                     per.column.gcv = TRUE,
                                     filter.type = NULL,
                                     n.candidates = 40L,
                                     n.cores = 1L,
                                     verbose = FALSE) {

    ## ================================================================
    ## INPUT VALIDATION
    ## ================================================================

    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be a 'knn.riem.fit' object from fit.rdgraph.regression()")
    }

    if (is.null(fitted.model$spectral)) {
        stop("Fitted model does not contain spectral component. ",
             "Cannot refit without cached eigendecomposition.")
    }

    ## Validate n.cores
    n.cores <- as.integer(n.cores)
    if (n.cores < 1L) {
        warning("n.cores must be >= 1; setting to 1")
        n.cores <- 1L
    }

    ## Extract spectral components
    V <- fitted.model$spectral$eigenvectors
    eigenvalues <- fitted.model$spectral$filtered.eigenvalues
    n <- nrow(V)
    m <- ncol(V)

    ## Validate and prepare y.new
    is.matrix.input <- is.matrix(y.new)

    if (is.matrix.input) {
        if (nrow(y.new) != n) {
            stop(sprintf("Number of rows in y.new (%d) must match fitted model (%d)",
                         nrow(y.new), n))
        }
        n.responses <- ncol(y.new)
        col.names <- colnames(y.new)
    } else {
        if (length(y.new) != n) {
            stop(sprintf("Length of y.new (%d) must match fitted model (%d)",
                         length(y.new), n))
        }
        y.new <- matrix(y.new, ncol = 1)
        n.responses <- 1L
        col.names <- NULL
    }

    ## Check for non-finite values
    if (any(!is.finite(y.new))) {
        bad.cols <- which(apply(y.new, 2, function(x) any(!is.finite(x))))
        stop(sprintf("y.new contains non-finite values in column(s): %s",
                     paste(bad.cols, collapse = ", ")))
    }

    ## ================================================================
    ## BRANCH: Per-column GCV vs. fixed filter
    ## ================================================================

    if (per.column.gcv) {

        start.time <- Sys.time()

        ## Determine filter type
        if (is.null(filter.type)) {
            filter.type <- fitted.model$spectral$filter.type
            if (is.null(filter.type)) {
                filter.type <- "heat_kernel"
                if (verbose) {
                    message("No filter.type specified and none in fitted model; ",
                            "using 'heat_kernel'")
                }
            }
        }

        filter.type <- match.arg(filter.type,
                                  c("heat_kernel", "tikhonov", "cubic_spline",
                                    "gaussian", "exponential", "butterworth"))

        ## Project all responses onto eigenbasis once: V^T Y (m x p)
        Vt.Y <- crossprod(V, y.new)

        ## Generate eta grid
        eta.grid <- generate.eta.grid(eigenvalues, filter.type, n.candidates)

        ## Pre-compute filter weights for all eta values: m x n.candidates
        filter.weights.matrix <- compute.filter.weights.matrix(
            eigenvalues, eta.grid, filter.type
        )

        ## Pre-compute trace(S) for all eta (used in all columns)
        trace.S.all <- colSums(filter.weights.matrix)

        ## ============================================================
        ## PARALLEL vs SEQUENTIAL PROCESSING
        ## ============================================================

        use.parallel <- (n.cores > 1L) && (n.responses > 1L)

        if (use.parallel) {
            ## Check for required packages
            if (!requireNamespace("foreach", quietly = TRUE)) {
                warning("Package 'foreach' not available; falling back to sequential processing")
                use.parallel <- FALSE
            }
            if (use.parallel && !requireNamespace("doParallel", quietly = TRUE)) {
                warning("Package 'doParallel' not available; falling back to sequential processing")
                use.parallel <- FALSE
            }
        }

        if (use.parallel) {

            ## --------------------------------------------------------
            ## PARALLEL EXECUTION
            ## --------------------------------------------------------

            ## Determine actual number of cores to use
            max.cores <- parallel::detectCores(logical = FALSE)
            n.cores.use <- min(n.cores, max.cores, n.responses)

            if (verbose) {
                message(sprintf("Processing %d columns in parallel using %d cores...",
                                n.responses, n.cores.use))
            }

            ## Register parallel backend
            cl <- parallel::makeCluster(n.cores.use)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            ## Define iterator variable for foreach
            j <- NULL  # Avoid R CMD check NOTE about undefined global

            ## Parallel loop over columns
            results.list <- foreach::foreach(
                j = seq_len(n.responses),
                .combine = "rbind",
                .packages = character(0),
                .export = c("select.eta.gcv.single.fast")
            ) %dopar% {

                y.obs.j <- y.new[, j]
                y.spectral.j <- Vt.Y[, j]

                ## Fast GCV selection (no intermediate y.hat storage)
                gcv.result <- select.eta.gcv.single.fast(
                    y.obs.j, y.spectral.j, V,
                    filter.weights.matrix, eta.grid, trace.S.all
                )

                ## Return as single-row data frame for rbind combining
                data.frame(
                    col.idx = j,
                    eta.optimal = gcv.result$eta.optimal,
                    gcv.min = gcv.result$gcv.min,
                    effective.df = gcv.result$effective.df,
                    best.idx = gcv.result$best.idx
                )
            }

            ## Reconstruct Y.hat using optimal indices
            Y.hat <- matrix(0, nrow = n, ncol = n.responses)
            for (i in seq_len(nrow(results.list))) {
                j <- results.list$col.idx[i]
                best.idx <- results.list$best.idx[i]
                y.spectral.j <- Vt.Y[, j]
                filtered.spectral <- y.spectral.j * filter.weights.matrix[, best.idx]
                Y.hat[, j] <- V %*% filtered.spectral
            }

            eta.optimal <- results.list$eta.optimal
            gcv.scores <- results.list$gcv.min
            effective.df <- results.list$effective.df
            n.cores.used <- n.cores.use

        } else {

            ## --------------------------------------------------------
            ## SEQUENTIAL EXECUTION
            ## --------------------------------------------------------

            n.cores.used <- 1L

            ## Initialize output storage
            Y.hat <- matrix(0, nrow = n, ncol = n.responses)
            eta.optimal <- numeric(n.responses)
            gcv.scores <- numeric(n.responses)
            effective.df <- numeric(n.responses)

            if (verbose && n.responses > 1L) {
                message(sprintf("Selecting eta via GCV for %d response(s)...",
                                n.responses))
                pb <- txtProgressBar(min = 0, max = n.responses, style = 3)
            }

            ## Process each column
            for (j in seq_len(n.responses)) {

                y.spectral <- Vt.Y[, j]

                ## Compute GCV scores for all eta candidates (vectorized)
                gcv.result <- select.eta.gcv.single(
                    y.new[, j], y.spectral, V, filter.weights.matrix, eta.grid
                )

                Y.hat[, j] <- gcv.result$y.hat
                eta.optimal[j] <- gcv.result$eta.optimal
                gcv.scores[j] <- gcv.result$gcv.min
                effective.df[j] <- gcv.result$effective.df

                if (verbose && n.responses > 1L) setTxtProgressBar(pb, j)
            }

            if (verbose && n.responses > 1L) {
                close(pb)
            }
        }

        elapsed.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

        if (verbose) {
            message(sprintf("Done. Processed %d columns in %.1f seconds (%.1f cols/sec)",
                            n.responses, elapsed.time, n.responses / elapsed.time))
            message(sprintf("Eta range: [%.2e, %.2e]",
                            min(eta.optimal), max(eta.optimal)))
        }

        ## Preserve column names if present
        if (!is.null(col.names)) {
            colnames(Y.hat) <- col.names
            names(eta.optimal) <- col.names
            names(gcv.scores) <- col.names
            names(effective.df) <- col.names
        }

        ## Compute residuals
        residuals <- y.new - Y.hat

        ## Build result
        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = TRUE,
            filter.type = filter.type,
            eta.optimal = if (n.responses == 1L) eta.optimal[1] else eta.optimal,
            gcv.scores = if (n.responses == 1L) gcv.scores[1] else gcv.scores,
            effective.df = if (n.responses == 1L) effective.df[1] else effective.df,
            eta.grid = eta.grid,
            n.cores.used = n.cores.used,
            elapsed.time = elapsed.time
        )

    } else {

        ## ============================================================
        ## ORIGINAL BEHAVIOR: Use fixed filtered eigenvalues from fit
        ## ============================================================

        f.lambda <- fitted.model$spectral$filtered.eigenvalues

        if (is.null(f.lambda)) {
            stop("Fitted model does not contain filtered.eigenvalues. ",
                 "Use per.column.gcv = TRUE or refit the original model.")
        }

        ## Apply spectral filter: Y_hat = V diag(f_lambda) V^T Y
        Vt.Y <- crossprod(V, y.new)           # V^T Y: m x p
        filtered <- f.lambda * Vt.Y           # Element-wise broadcast
        Y.hat <- V %*% filtered               # n x p

        ## Compute residuals
        residuals <- y.new - Y.hat

        ## Build result
        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = FALSE
        )
    }

    class(result) <- c("knn.riem.refit", "list")
    return(result)
}


## ============================================================
## HELPER FUNCTIONS FOR PER-COLUMN GCV
## ============================================================

#' Generate logarithmically-spaced eta grid for GCV search
#'
#' @param eigenvalues Vector of eigenvalues from spectral decomposition
#' @param filter.type Character string specifying filter type
#' @param n.candidates Number of grid points
#' @return Numeric vector of eta candidates
#' @keywords internal
generate.eta.grid <- function(eigenvalues, filter.type, n.candidates) {

    eps <- 1e-10
    lambda.max <- max(eigenvalues)

    ## Determine eta bounds based on filter type
    eta.min <- eps

    eta.max <- switch(filter.type,
        "heat_kernel" = ,
        "gaussian" = ,
        "exponential" = if (lambda.max > 0) -log(eps) / lambda.max else 1.0,
        "tikhonov" = ,
        "cubic_spline" = ,
        "butterworth" = 1.0 / eps,
        1.0  # default
    )

    ## Generate log-spaced grid
    t.frac <- seq_len(n.candidates) / n.candidates
    eta.grid <- exp(log(eta.min) + t.frac * (log(eta.max) - log(eta.min)))

    return(eta.grid)
}


#' Compute filter weights for all eigenvalues and all eta candidates
#'
#' @param eigenvalues Vector of eigenvalues (length m)
#' @param eta.grid Vector of eta candidates (length n.candidates)
#' @param filter.type Character string specifying filter type
#' @return Matrix of filter weights (m x n.candidates)
#' @keywords internal
compute.filter.weights.matrix <- function(eigenvalues, eta.grid, filter.type) {

    m <- length(eigenvalues)
    n.eta <- length(eta.grid)

    ## Create outer product structure for vectorized computation
    ## Lambda: m x 1, Eta: 1 x n.eta -> Lambda.Eta: m x n.eta
    lambda.mat <- matrix(eigenvalues, nrow = m, ncol = n.eta)
    eta.mat <- matrix(eta.grid, nrow = m, ncol = n.eta, byrow = TRUE)

    weights <- switch(filter.type,
        "heat_kernel" = exp(-eta.mat * lambda.mat),
        "tikhonov" = 1 / (1 + eta.mat * lambda.mat),
        "cubic_spline" = 1 / (1 + eta.mat * lambda.mat^2),
        "gaussian" = exp(-eta.mat * lambda.mat^2),
        "exponential" = exp(-eta.mat * sqrt(pmax(lambda.mat, 0))),
        "butterworth" = {
            ratio <- lambda.mat / eta.mat
            1 / (1 + ratio^4)
        },
        stop("Unknown filter type: ", filter.type)
    )

    return(weights)
}


#' Select optimal eta for a single response column via GCV
#'
#' @param y.obs Observed response vector (length n)
#' @param y.spectral Spectral coefficients V^T y (length m)
#' @param V Eigenvector matrix (n x m)
#' @param filter.weights.matrix Pre-computed filter weights (m x n.candidates)
#' @param eta.grid Vector of eta candidates
#' @return List with y.hat, eta.optimal, gcv.min, effective.df
#' @keywords internal
select.eta.gcv.single <- function(y.obs, y.spectral, V,
                                   filter.weights.matrix, eta.grid) {

    n <- length(y.obs)
    n.eta <- length(eta.grid)

    ## Compute filtered spectral coefficients for all eta: m x n.eta
    ## y.spectral is m x 1, filter.weights.matrix is m x n.eta
    filtered.spectral <- y.spectral * filter.weights.matrix  # broadcast

    ## Compute predictions for all eta: n x n.eta
    ## V is n x m, filtered.spectral is m x n.eta
    Y.hat.all <- V %*% filtered.spectral

    ## Compute GCV scores for all eta (vectorized)
    ## Residuals: n x n.eta
    residuals.all <- y.obs - Y.hat.all

    ## RSS for each eta
    rss <- colSums(residuals.all^2)

    ## Effective df = trace(S) = sum of filter weights
    trace.S <- colSums(filter.weights.matrix)

    ## GCV score = RSS / (n - trace(S))^2
    denom <- pmax(n - trace.S, 1e-10)
    gcv.scores <- rss / (denom^2)

    ## Find optimal
    best.idx <- which.min(gcv.scores)

    return(list(
        y.hat = Y.hat.all[, best.idx],
        eta.optimal = eta.grid[best.idx],
        gcv.min = gcv.scores[best.idx],
        effective.df = trace.S[best.idx]
    ))
}


#' Fast GCV selection without storing intermediate predictions
#'
#' This version is optimized for parallel execution where we only need
#' the optimal index, not the full prediction matrix. The actual y.hat
#' is computed after parallel collection using the returned best.idx.
#'
#' @param y.obs Observed response vector (length n)
#' @param y.spectral Spectral coefficients V^T y (length m)
#' @param V Eigenvector matrix (n x m)
#' @param filter.weights.matrix Pre-computed filter weights (m x n.candidates)
#' @param eta.grid Vector of eta candidates
#' @param trace.S.all Pre-computed trace(S) for all eta candidates
#' @return List with eta.optimal, gcv.min, effective.df, best.idx
#' @keywords internal
select.eta.gcv.single.fast <- function(y.obs, y.spectral, V,
                                        filter.weights.matrix, eta.grid,
                                        trace.S.all) {

    n <- length(y.obs)
    n.eta <- length(eta.grid)

    ## Compute filtered spectral coefficients for all eta: m x n.eta
    filtered.spectral <- y.spectral * filter.weights.matrix

    ## Compute predictions for all eta: n x n.eta
    Y.hat.all <- V %*% filtered.spectral

    ## Compute RSS for all eta (vectorized)
    residuals.all <- y.obs - Y.hat.all
    rss <- colSums(residuals.all^2)

    ## GCV score = RSS / (n - trace(S))^2
    denom <- pmax(n - trace.S.all, 1e-10)
    gcv.scores <- rss / (denom^2)

    ## Find optimal
    best.idx <- which.min(gcv.scores)

    return(list(
        eta.optimal = eta.grid[best.idx],
        gcv.min = gcv.scores[best.idx],
        effective.df = trace.S.all[best.idx],
        best.idx = best.idx
    ))
}


## ============================================================
## UPDATED PRINT AND SUMMARY METHODS
## ============================================================

#' @export
print.knn.riem.refit <- function(x, ...) {
    cat("\nRefitted Riemannian Graph Regression\n")
    cat("====================================\n\n")

    if (x$n.responses == 1) {
        cat(sprintf("Single response variable (n = %d)\n", length(x$fitted.values)))
        cat(sprintf("RMSE: %.4f\n", sqrt(mean(x$residuals^2))))
        cat(sprintf("MAE:  %.4f\n", mean(abs(x$residuals))))

        if (isTRUE(x$per.column.gcv)) {
            cat(sprintf("\nPer-column GCV (%s filter):\n", x$filter.type))
            cat(sprintf("  Optimal eta:    %.4e\n", x$eta.optimal))
            cat(sprintf("  Effective df:   %.1f\n", x$effective.df))
            cat(sprintf("  GCV score:      %.4e\n", x$gcv.scores))
        }
    } else {
        cat(sprintf("Multiple responses (n = %d, p = %d)\n",
                    nrow(x$fitted.values), x$n.responses))
        rmse.by.col <- apply(x$residuals, 2, function(r) sqrt(mean(r^2)))
        cat(sprintf("Mean RMSE across responses: %.4f\n", mean(rmse.by.col)))
        cat(sprintf("RMSE range: [%.4f, %.4f]\n", min(rmse.by.col), max(rmse.by.col)))

        if (isTRUE(x$per.column.gcv)) {
            cat(sprintf("\nPer-column GCV selection (%s filter):\n", x$filter.type))
            cat(sprintf("  Eta range:          [%.2e, %.2e]\n",
                        min(x$eta.optimal), max(x$eta.optimal)))
            cat(sprintf("  Effective df range: [%.1f, %.1f]\n",
                        min(x$effective.df), max(x$effective.df)))

            if (!is.null(x$n.cores.used) && !is.null(x$elapsed.time)) {
                cat(sprintf("\nComputation: %.1f sec using %d core%s (%.1f cols/sec)\n",
                            x$elapsed.time, x$n.cores.used,
                            if (x$n.cores.used > 1) "s" else "",
                            x$n.responses / x$elapsed.time))
            }
        }
    }

    invisible(x)
}


#' @export
summary.knn.riem.refit <- function(object, ...) {
    n.resp <- object$n.responses
    n.obs <- if (n.resp == 1) length(object$fitted.values) else nrow(object$fitted.values)

    if (n.resp == 1) {
        ## Single response summary
        y.hat <- object$fitted.values
        resid <- object$residuals
        y.obs <- y.hat + resid

        rss <- sum(resid^2)
        tss <- sum((y.obs - mean(y.obs))^2)
        r.squared <- 1 - rss / tss
        rmse <- sqrt(mean(resid^2))
        mae <- mean(abs(resid))

        result <- list(
            n.obs = n.obs,
            n.responses = 1L,
            r.squared = r.squared,
            rmse = rmse,
            mae = mae,
            per.column.gcv = isTRUE(object$per.column.gcv),
            eta.optimal = object$eta.optimal,
            effective.df = object$effective.df,
            filter.type = object$filter.type,
            n.cores.used = object$n.cores.used,
            elapsed.time = object$elapsed.time
        )
    } else {
        ## Multiple response summary
        Y.hat <- object$fitted.values
        residuals <- object$residuals
        Y.obs <- Y.hat + residuals

        ## Per-column metrics
        rmse.by.col <- apply(residuals, 2, function(r) sqrt(mean(r^2)))
        mae.by.col <- apply(residuals, 2, function(r) mean(abs(r)))

        r.squared.by.col <- sapply(seq_len(n.resp), function(j) {
            rss <- sum(residuals[, j]^2)
            tss <- sum((Y.obs[, j] - mean(Y.obs[, j]))^2)
            if (tss > 1e-15) 1 - rss / tss else NA_real_
        })

        result <- list(
            n.obs = n.obs,
            n.responses = n.resp,
            rmse.by.column = rmse.by.col,
            mae.by.column = mae.by.col,
            r.squared.by.column = r.squared.by.col,
            per.column.gcv = isTRUE(object$per.column.gcv),
            eta.optimal = object$eta.optimal,
            effective.df = object$effective.df,
            filter.type = object$filter.type,
            n.cores.used = object$n.cores.used,
            elapsed.time = object$elapsed.time
        )
    }

    class(result) <- "summary.knn.riem.refit"
    return(result)
}


#' @export
print.summary.knn.riem.refit <- function(x, digits = 4, ...) {
    cat("\nSummary: Refitted Riemannian Graph Regression\n")
    cat("==============================================\n\n")

    cat(sprintf("Observations: %d\n", x$n.obs))
    cat(sprintf("Responses: %d\n", x$n.responses))

    if (x$n.responses == 1) {
        cat(sprintf("\nFit quality:\n"))
        cat(sprintf("  R-squared: %.*f\n", digits, x$r.squared))
        cat(sprintf("  RMSE:      %.*f\n", digits, x$rmse))
        cat(sprintf("  MAE:       %.*f\n", digits, x$mae))

        if (x$per.column.gcv) {
            cat(sprintf("\nGCV-selected parameters (%s filter):\n", x$filter.type))
            cat(sprintf("  Optimal eta:    %.*e\n", digits, x$eta.optimal))
            cat(sprintf("  Effective df:   %.1f\n", x$effective.df))
        }
    } else {
        cat(sprintf("\nFit quality summary across %d responses:\n", x$n.responses))

        print.quantile.summary <- function(vals, name) {
            q <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
            cat(sprintf("  %s: min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n",
                        name, q[1], q[2], q[3], q[4], q[5]))
        }

        print.quantile.summary(x$r.squared.by.column, "R-sq")
        print.quantile.summary(x$rmse.by.column, "RMSE")

        if (x$per.column.gcv) {
            cat(sprintf("\nGCV-selected parameters (%s filter):\n", x$filter.type))
            print.quantile.summary(x$eta.optimal, "Eta ")
            print.quantile.summary(x$effective.df, "Eff.df")

            if (!is.null(x$n.cores.used) && !is.null(x$elapsed.time)) {
                cat(sprintf("\nComputation: %.1f sec using %d core%s\n",
                            x$elapsed.time, x$n.cores.used,
                            if (x$n.cores.used > 1) "s" else ""))
            }
        }
    }

    invisible(x)
}

#' Refit Riemannian Graph Regression with New Response(s)
#'
#' Apply spectral filtering from a fitted model to new response variable(s),
#' with optional Bayesian posterior inference for uncertainty quantification.
#'
#' This function uses the cached eigendecomposition from an existing fit,
#' avoiding the computational cost of rebuilding the geometric structure.
#' When posterior sampling is enabled, vertex-wise credible intervals are
#' computed using Monte Carlo sampling from the posterior distribution of
#' the smoothed response.
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
#'   \code{spectral} component with cached eigendecomposition. For per-column
#'   GCV, the model must include raw eigenvalues.
#'
#' @param y.new New response variable(s). Either a numeric vector of length n,
#'   or a numeric matrix with n rows where each column is a separate response.
#'
#' @param per.column.gcv Logical. If \code{TRUE}, select optimal smoothing
#'   parameter eta independently for each column via GCV. If \code{FALSE}
#'   (default), use the filtered eigenvalues from the original fit for all
#'   columns. Setting \code{TRUE} is recommended when columns represent
#'   features with heterogeneous smoothness (e.g., different phylotypes).
#'   Requires fitted model with raw eigenvalues stored.
#'
#' @param n.candidates Integer. Number of candidate eta values for GCV grid
#'   search. Only used when \code{per.column.gcv = TRUE}. Default is 40.
#'
#' @param n.cores Integer. Number of CPU cores for parallel processing.
#'   Default is 1 (sequential). When \code{n.cores > 1}, columns are processed
#'   in parallel using \pkg{foreach} and \pkg{doParallel}.
#'
#' @param verbose Logical. If \code{TRUE}, print progress information during
#'   per-column GCV selection and posterior computation. Default is \code{FALSE}.
#'
#' @param with.posterior Logical. If \code{TRUE}, compute Bayesian posterior
#'   credible intervals for the smoothed response at each vertex. Default is
#'   \code{FALSE}. Posterior computation adds moderate computational cost
#'   (approximately 20-50% overhead per column).
#'
#' @param return.posterior.samples Logical. If \code{TRUE} and
#'   \code{with.posterior = TRUE}, return the full matrix of posterior samples.
#'   Default is \code{FALSE}. For multi-column input with many responses,
#'   storing samples may require substantial memory (n * n.posterior.samples *
#'   n.responses * 8 bytes).
#'
#' @param credible.level Numeric in (0, 1). Coverage probability for credible
#'   intervals. Default is 0.95 for 95% credible intervals.
#'
#' @param n.posterior.samples Integer. Number of Monte Carlo samples for
#'   posterior inference. Larger values provide more accurate quantile estimates.
#'   Default is 1000. Values below 100 are not recommended.
#'
#' @param posterior.seed Integer. Random seed for reproducible posterior samples.
#'   Default is 12345.
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
#'     \item{\code{posterior}}{(Only if \code{with.posterior = TRUE}) A list containing:
#'       \itemize{
#'         \item \code{lower}: Lower credible bounds (vector or matrix).
#'         \item \code{upper}: Upper credible bounds (vector or matrix).
#'         \item \code{sd}: Posterior standard deviations (vector or matrix).
#'         \item \code{sigma}: Estimated residual SD (scalar or vector).
#'         \item \code{credible.level}: Coverage probability used.
#'         \item \code{samples}: (Optional) Posterior samples if requested.
#'       }
#'     }
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
#' \strong{Posterior inference:} When \code{with.posterior = TRUE}, the function
#' computes Bayesian credible intervals by drawing samples from the posterior
#' distribution of the smoothed response. The posterior follows from a Bayesian
#' interpretation of spectral smoothing where coefficients in the eigenbasis
#' have posterior variance inversely proportional to (1 + eta * lambda_j).
#'
#' \strong{Parallel processing:} When \code{n.cores > 1}, the function uses
#' \pkg{foreach} with \pkg{doParallel} backend. Each worker receives a copy of
#' the eigenvector matrix V and filter weights, so memory usage scales with
#' the number of cores. For very large n (>10000), consider limiting n.cores
#' to avoid memory pressure.
#'
#' @examples
#' \dontrun{
#' ## Fit initial model on binary outcome
#' fit <- fit.rdgraph.regression(X, y.binary, k = 15)
#'
#' ## Quick refit using same smoothing (original behavior)
#' Z.hat.quick <- refit.rdgraph.regression(fit, Z.abundances)
#'
#' ## Smooth phylotype abundances with per-column GCV (parallel)
#' Z.hat <- refit.rdgraph.regression(fit, Z.abundances,
#'                                    per.column.gcv = TRUE,
#'                                    n.cores = 4)
#'
#' ## With posterior inference for uncertainty quantification
#' Z.hat.uq <- refit.rdgraph.regression(fit, Z.abundances,
#'                                       per.column.gcv = TRUE,
#'                                       with.posterior = TRUE,
#'                                       credible.level = 0.95,
#'                                       n.posterior.samples = 1000)
#'
#' ## Access credible intervals
#' CI.width <- Z.hat.uq$posterior$upper - Z.hat.uq$posterior$lower
#' }
#'
#' @seealso \code{\link{fit.rdgraph.regression}} for initial model fitting
#'
#' @export
refit.rdgraph.regression <- function(fitted.model,
                                      y.new,
                                      per.column.gcv = FALSE,
                                      n.candidates = 40L,
                                      n.cores = 1L,
                                      verbose = FALSE,
                                      with.posterior = FALSE,
                                      return.posterior.samples = FALSE,
                                      credible.level = 0.95,
                                      n.posterior.samples = 1000L,
                                      posterior.seed = 12345L) {

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

    ## Validate posterior parameters
    if (with.posterior) {
        if (credible.level <= 0 || credible.level >= 1) {
            stop("credible.level must be in (0, 1)")
        }
        n.posterior.samples <- as.integer(n.posterior.samples)
        if (n.posterior.samples < 100L) {
            warning("n.posterior.samples < 100 may produce unreliable quantiles")
        }
        posterior.seed <- as.integer(posterior.seed)
    }

    ## Extract spectral components
    V <- fitted.model$spectral$eigenvectors
    f.lambda <- fitted.model$spectral$filtered.eigenvalues
    n <- nrow(V)
    m <- ncol(V)

    ## Check for raw eigenvalues (required for per-column GCV and posterior)
    has.raw.eigenvalues <- !is.null(fitted.model$spectral$eigenvalues)

    if (per.column.gcv && !has.raw.eigenvalues) {
        stop("Per-column GCV requires raw eigenvalues in the fitted model.\n",
             "Your fitted model only contains filtered.eigenvalues.\n",
             "Please refit the model with an updated version of gflow that\n",
             "exports raw eigenvalues, or use per.column.gcv = FALSE.")
    }

    if (with.posterior && !has.raw.eigenvalues) {
        stop("Posterior inference requires raw eigenvalues in the fitted model.\n",
             "Your fitted model only contains filtered.eigenvalues.\n",
             "Please refit the model with an updated version of gflow.")
    }

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

    ## Warn about memory for posterior samples with many columns
    if (with.posterior && return.posterior.samples && n.responses > 10L) {
        mem.mb <- n * n.posterior.samples * n.responses * 8 / 1e6
        warning(sprintf(
            "Storing posterior samples for %d columns will require ~%.0f MB.\n",
            n.responses, mem.mb),
            "Consider setting return.posterior.samples = FALSE.",
            call. = FALSE)
    }

    ## ================================================================
    ## BRANCH: Per-column GCV vs. fixed filter
    ## ================================================================

    if (per.column.gcv) {

        start.time <- Sys.time()

        ## Extract raw eigenvalues and filter type from fitted model
        eigenvalues <- fitted.model$spectral$eigenvalues
        filter.type <- fitted.model$spectral$filter.type

        ## Default to heat_kernel if filter.type not stored (older models)
        if (is.null(filter.type)) {
            filter.type <- "heat_kernel"
            if (verbose) {
                message("Filter type not stored in model; using 'heat_kernel'")
            }
        }

        if (verbose) {
            message(sprintf("Using '%s' filter for per-column GCV", filter.type))
        }

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
            best.idx.vec <- results.list$best.idx
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
            best.idx.vec <- integer(n.responses)

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
                best.idx.vec[j] <- gcv.result$best.idx

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

        ## ============================================================
        ## POSTERIOR INFERENCE (PER-COLUMN GCV CASE)
        ## ============================================================

        posterior <- NULL

        if (with.posterior) {

            if (verbose) {
                message("Computing posterior credible intervals...")
            }

            ## Initialize storage for posterior summaries
            posterior.lower <- matrix(0, nrow = n, ncol = n.responses)
            posterior.upper <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sd <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sigma <- numeric(n.responses)

            ## Optionally store samples
            if (return.posterior.samples) {
                posterior.samples <- vector("list", n.responses)
            }

            if (verbose && n.responses > 1L) {
                pb.post <- txtProgressBar(min = 0, max = n.responses, style = 3)
            }

            for (j in seq_len(n.responses)) {

                ## Get the filtered eigenvalues for this column's optimal eta
                filtered.eigenvalues.j <- filter.weights.matrix[, best.idx.vec[j]]

                ## Call C++ posterior computation
                post.j <- .Call(
                    S_compute_posterior_summary,
                    V,
                    eigenvalues,
                    filtered.eigenvalues.j,
                    y.new[, j],
                    Y.hat[, j],
                    eta.optimal[j],
                    credible.level,
                    n.posterior.samples,
                    posterior.seed + j - 1L,  # Different seed per column
                    return.posterior.samples,
                    PACKAGE = "gflow"
                )

                posterior.lower[, j] <- post.j$lower
                posterior.upper[, j] <- post.j$upper
                posterior.sd[, j] <- post.j$sd
                posterior.sigma[j] <- post.j$sigma

                if (return.posterior.samples) {
                    posterior.samples[[j]] <- post.j$samples
                }

                if (verbose && n.responses > 1L) {
                    setTxtProgressBar(pb.post, j)
                }
            }

            if (verbose && n.responses > 1L) {
                close(pb.post)
            }

            ## Preserve column names
            if (!is.null(col.names)) {
                colnames(posterior.lower) <- col.names
                colnames(posterior.upper) <- col.names
                colnames(posterior.sd) <- col.names
                names(posterior.sigma) <- col.names
                if (return.posterior.samples) {
                    names(posterior.samples) <- col.names
                }
            }

            ## Assemble posterior component
            posterior <- list(
                lower = if (n.responses == 1L) as.vector(posterior.lower) else posterior.lower,
                upper = if (n.responses == 1L) as.vector(posterior.upper) else posterior.upper,
                sd = if (n.responses == 1L) as.vector(posterior.sd) else posterior.sd,
                sigma = if (n.responses == 1L) posterior.sigma[1] else posterior.sigma,
                credible.level = credible.level,
                n.samples = n.posterior.samples
            )

            if (return.posterior.samples) {
                posterior$samples <- if (n.responses == 1L) {
                    posterior.samples[[1]]
                } else {
                    posterior.samples
                }
            }

            if (verbose) {
                message("Posterior computation complete.")
            }
        }

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

        if (with.posterior) {
            result$posterior <- posterior
        }

    } else {

        ## ============================================================
        ## ORIGINAL BEHAVIOR: Use fixed filtered eigenvalues from fit
        ## ============================================================

        if (is.null(f.lambda)) {
            stop("Fitted model does not contain filtered.eigenvalues. ",
                 "Cannot refit without spectral filter information.")
        }

        ## Apply spectral filter: Y_hat = V diag(f_lambda) V^T Y
        Vt.Y <- crossprod(V, y.new)           # V^T Y: m x p
        filtered <- f.lambda * Vt.Y           # Element-wise broadcast
        Y.hat <- V %*% filtered               # n x p

        ## Compute residuals
        residuals <- y.new - Y.hat

        ## ============================================================
        ## POSTERIOR INFERENCE (FIXED FILTER CASE)
        ## ============================================================

        posterior <- NULL

        if (with.posterior) {

            ## For fixed filter case, we need the original eta
            ## Try to extract from gcv component
            eta.original <- NULL

            if (!is.null(fitted.model$gcv) && !is.null(fitted.model$gcv$eta.optimal)) {
                ## Use the last iteration's optimal eta
                eta.vec <- fitted.model$gcv$eta.optimal
                eta.original <- eta.vec[length(eta.vec)]
            } else if (!is.null(fitted.model$iteration) &&
                       !is.null(fitted.model$iteration$eta.optimal)) {
                eta.vec <- fitted.model$iteration$eta.optimal
                eta.original <- eta.vec[length(eta.vec)]
            }

            if (is.null(eta.original)) {
                stop("Cannot perform posterior inference without eta value.\n",
                     "The fitted model does not contain gcv$eta.optimal.\n",
                     "Please use per.column.gcv = TRUE for posterior inference,\n",
                     "or refit with a newer version of gflow.")
            }

            eigenvalues <- fitted.model$spectral$eigenvalues

            if (verbose) {
                message(sprintf("Computing posterior with eta = %.4e from original fit",
                                eta.original))
            }

            ## Initialize storage for posterior summaries
            posterior.lower <- matrix(0, nrow = n, ncol = n.responses)
            posterior.upper <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sd <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sigma <- numeric(n.responses)

            if (return.posterior.samples) {
                posterior.samples <- vector("list", n.responses)
            }

            if (verbose && n.responses > 1L) {
                message(sprintf("Computing posterior for %d response(s)...", n.responses))
                pb.post <- txtProgressBar(min = 0, max = n.responses, style = 3)
            }

            for (j in seq_len(n.responses)) {

                ## Call C++ posterior computation
                ## All columns use the same eta and filtered eigenvalues
                post.j <- .Call(
                    S_compute_posterior_summary,
                    V,
                    eigenvalues,
                    f.lambda,
                    y.new[, j],
                    Y.hat[, j],
                    eta.original,
                    credible.level,
                    n.posterior.samples,
                    posterior.seed + j - 1L,
                    return.posterior.samples,
                    PACKAGE = "gflow"
                )

                posterior.lower[, j] <- post.j$lower
                posterior.upper[, j] <- post.j$upper
                posterior.sd[, j] <- post.j$sd
                posterior.sigma[j] <- post.j$sigma

                if (return.posterior.samples) {
                    posterior.samples[[j]] <- post.j$samples
                }

                if (verbose && n.responses > 1L) {
                    setTxtProgressBar(pb.post, j)
                }
            }

            if (verbose && n.responses > 1L) {
                close(pb.post)
            }

            ## Preserve column names
            if (!is.null(col.names)) {
                colnames(posterior.lower) <- col.names
                colnames(posterior.upper) <- col.names
                colnames(posterior.sd) <- col.names
                names(posterior.sigma) <- col.names
                if (return.posterior.samples) {
                    names(posterior.samples) <- col.names
                }
            }

            ## Assemble posterior component
            posterior <- list(
                lower = if (n.responses == 1L) as.vector(posterior.lower) else posterior.lower,
                upper = if (n.responses == 1L) as.vector(posterior.upper) else posterior.upper,
                sd = if (n.responses == 1L) as.vector(posterior.sd) else posterior.sd,
                sigma = if (n.responses == 1L) posterior.sigma[1] else posterior.sigma,
                credible.level = credible.level,
                n.samples = n.posterior.samples,
                eta.used = eta.original
            )

            if (return.posterior.samples) {
                posterior$samples <- if (n.responses == 1L) {
                    posterior.samples[[1]]
                } else {
                    posterior.samples
                }
            }

            if (verbose) {
                message("Posterior computation complete.")
            }
        }

        ## Build result
        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = FALSE
        )

        if (with.posterior) {
            result$posterior <- posterior
        }
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
#' @param eigenvalues Vector of eigenvalues
#' @param eta.grid Vector of eta values
#' @param filter.type Filter type string
#' @return Matrix of filter weights (m x length(eta.grid))
#' @keywords internal
compute.filter.weights.matrix <- function(eigenvalues, eta.grid, filter.type) {

    m <- length(eigenvalues)
    n.eta <- length(eta.grid)

    ## Pre-allocate matrix (eigenvalues in rows, eta in columns)
    weights <- matrix(0, nrow = m, ncol = n.eta)

    for (k in seq_len(n.eta)) {
        eta <- eta.grid[k]

        weights[, k] <- switch(filter.type,
            "heat_kernel" = exp(-eta * eigenvalues),
            "tikhonov" = 1.0 / (1.0 + eta * eigenvalues),
            "cubic_spline" = 1.0 / (1.0 + eta * eigenvalues^2),
            "gaussian" = exp(-eta * eigenvalues^2),
            "exponential" = exp(-eta * sqrt(pmax(eigenvalues, 0))),
            "butterworth" = {
                ## Order-2 Butterworth: 1/(1 + (lambda/eta)^4)
                x <- eigenvalues / eta
                1.0 / (1.0 + x^4)
            },
            ## Default to heat kernel
            exp(-eta * eigenvalues)
        )
    }

    return(weights)
}


#' Single-column GCV selection with full fitted values
#'
#' @param y.obs Observed response vector
#' @param y.spectral Spectral coefficients V^T y
#' @param V Eigenvector matrix
#' @param filter.weights.matrix Pre-computed filter weights (m x n.eta)
#' @param eta.grid Vector of eta candidates
#' @return List with y.hat, eta.optimal, gcv.min, effective.df
#' @keywords internal
select.eta.gcv.single <- function(y.obs, y.spectral, V,
                                   filter.weights.matrix, eta.grid) {

    n <- length(y.obs)
    m <- length(y.spectral)
    n.eta <- length(eta.grid)

    ## Compute filtered spectral coefficients for all eta: m x n.eta
    filtered.spectral <- y.spectral * filter.weights.matrix

    ## Compute predictions for all eta: n x n.eta
    Y.hat.all <- V %*% filtered.spectral

    ## Compute RSS and effective df for all eta
    residuals.all <- y.obs - Y.hat.all
    rss <- colSums(residuals.all^2)
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
        effective.df = trace.S[best.idx],
        best.idx = best.idx
    ))
}


#' Fast single-column GCV selection (returns only scalars)
#'
#' Optimized version that doesn't return full y.hat vector.
#' Used in parallel processing where we reconstruct y.hat afterward.
#'
#' @param y.obs Observed response vector
#' @param y.spectral Spectral coefficients V^T y
#' @param V Eigenvector matrix
#' @param filter.weights.matrix Pre-computed filter weights (m x n.eta)
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
## PRINT AND SUMMARY METHODS
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

        if (!is.null(x$posterior)) {
            cat(sprintf("\nPosterior inference (%.0f%% credible intervals):\n",
                        x$posterior$credible.level * 100))
            cat(sprintf("  Residual sigma: %.4f\n", x$posterior$sigma))
            ci.width <- x$posterior$upper - x$posterior$lower
            cat(sprintf("  Mean CI width:  %.4f\n", mean(ci.width)))
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

        if (!is.null(x$posterior)) {
            cat(sprintf("\nPosterior inference (%.0f%% credible intervals):\n",
                        x$posterior$credible.level * 100))
            cat(sprintf("  Sigma range:    [%.4f, %.4f]\n",
                        min(x$posterior$sigma), max(x$posterior$sigma)))
            ci.width <- x$posterior$upper - x$posterior$lower
            mean.ci.width <- colMeans(ci.width)
            cat(sprintf("  Mean CI width:  [%.4f, %.4f]\n",
                        min(mean.ci.width), max(mean.ci.width)))
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

        if (!is.null(object$posterior)) {
            result$posterior.summary <- list(
                credible.level = object$posterior$credible.level,
                sigma = object$posterior$sigma,
                mean.ci.width = mean(object$posterior$upper - object$posterior$lower)
            )
        }
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

        if (!is.null(object$posterior)) {
            ci.width <- object$posterior$upper - object$posterior$lower
            result$posterior.summary <- list(
                credible.level = object$posterior$credible.level,
                sigma.range = range(object$posterior$sigma),
                mean.ci.width.by.col = colMeans(ci.width)
            )
        }
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

        if (!is.null(x$posterior.summary)) {
            cat(sprintf("\nPosterior inference (%.0f%% CI):\n",
                        x$posterior.summary$credible.level * 100))
            cat(sprintf("  Residual sigma: %.*f\n", digits, x$posterior.summary$sigma))
            cat(sprintf("  Mean CI width:  %.*f\n", digits, x$posterior.summary$mean.ci.width))
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

        if (!is.null(x$posterior.summary)) {
            cat(sprintf("\nPosterior inference (%.0f%% CI):\n",
                        x$posterior.summary$credible.level * 100))
            cat(sprintf("  Sigma range: [%.4f, %.4f]\n",
                        x$posterior.summary$sigma.range[1],
                        x$posterior.summary$sigma.range[2]))
            print.quantile.summary(x$posterior.summary$mean.ci.width.by.col, "CI width")
        }
    }

    invisible(x)
}

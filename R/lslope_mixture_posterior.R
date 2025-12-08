## ============================================================================
## Phase 2: Posterior Uncertainty Propagation for lslope
## ============================================================================
##
## This module extends the basic mixture model framework to propagate
## uncertainty from the spectral smoothing step into the lslope estimates.
## Rather than using point estimates of z.hat, we use posterior samples
## to obtain a distribution of lslope values at each vertex.
##
## The key insight is that the smoothed values z.hat have uncertainty
## characterized by the posterior distribution from spectral filtering.
## By computing lslope() for each posterior sample, we obtain a distribution
## of lslope values that incorporates this uncertainty.
##
## This provides:
## 1. More robust lslope estimates (posterior mean instead of point estimate)
## 2. Uncertainty quantification for lslope values
## 3. Better calibrated mixture model inference
## ============================================================================


#' Compute lslope with Posterior Uncertainty Propagation
#'
#' Propagates uncertainty from spectral smoothing into lslope estimates by
#' computing lslope for each posterior sample of the smoothed feature values.
#'
#' @param adj.list Adjacency list (1-based R indexing)
#' @param weight.list Edge weight list
#' @param y.hat Smoothed response values (length n). This is treated as fixed
#'   since it comes from the original fit.
#' @param Z.hat.samples List of posterior sample matrices for smoothed features.
#'   Each element is an n x B matrix of posterior samples for one feature,
#'   or a single n x B matrix if only one feature.
#' @param lslope.type Type of lslope measure: "normalized", "slope", or "sign"
#' @param ascending Logical. Use ascending gradient direction
#' @param return.samples Logical. If TRUE, return full sample matrices
#' @param verbose Logical. Print progress
#'
#' @return A list of class "lslope.posterior" containing:
#'   \describe{
#'     \item{mean}{Matrix (p x n) of posterior mean lslope values}
#'     \item{sd}{Matrix (p x n) of posterior standard deviations}
#'     \item{lower}{Matrix (p x n) of lower credible bounds (2.5%)}
#'     \item{upper}{Matrix (p x n) of upper credible bounds (97.5%)}
#'     \item{samples}{(Optional) List of sample matrices}
#'     \item{n.samples}{Number of posterior samples used}
#'   }
#'
#' @details
#' For each posterior sample b of the smoothed feature values z.hat^(b),
#' this function computes lslope(y.hat, z.hat^(b)). The resulting collection
#' of lslope values at each vertex characterizes the uncertainty in the
#' local association measure due to smoothing uncertainty.
#'
#' The posterior samples should come from \code{refit.rdgraph.regression()}
#' with \code{with.posterior = TRUE} and \code{return.posterior.samples = TRUE}.
#'
#' @examples
#' \dontrun{
#' ## Fit original model
#' fit <- fit.rdgraph.regression(X, y, k = 15, with.posterior = TRUE)
#'
#' ## Refit with posterior samples for features
#' Z.refit <- refit.rdgraph.regression(
#'     fit, Z[, 1:10],  # First 10 features
#'     per.column.gcv = TRUE,
#'     with.posterior = TRUE,
#'     return.posterior.samples = TRUE,
#'     n.posterior.samples = 500
#' )
#'
#' ## Compute lslope with uncertainty propagation
#' lslope.post <- lslope.with.posterior(
#'     fit$graph$adj.list,
#'     fit$graph$edge.length.list,
#'     fit$fitted.values,
#'     Z.refit$posterior$samples  # List of sample matrices
#' )
#' }
#'
#' @export
lslope.with.posterior <- function(adj.list,
                                   weight.list,
                                   y.hat,
                                   Z.hat.samples,
                                   lslope.type = c("normalized", "slope", "sign"),
                                   ascending = TRUE,
                                   return.samples = FALSE,
                                   verbose = TRUE) {

    lslope.type <- match.arg(lslope.type)

    n <- length(y.hat)

    ## Handle input format
    if (is.matrix(Z.hat.samples)) {
        ## Single feature: n x B matrix
        Z.hat.samples <- list(Z.hat.samples)
    }

    if (!is.list(Z.hat.samples)) {
        stop("Z.hat.samples must be a list of matrices or a single matrix")
    }

    p <- length(Z.hat.samples)
    B <- ncol(Z.hat.samples[[1]])

    ## Validate dimensions
    for (j in seq_len(p)) {
        if (nrow(Z.hat.samples[[j]]) != n) {
            stop(sprintf("Feature %d: nrow(samples) = %d, but n = %d",
                         j, nrow(Z.hat.samples[[j]]), n))
        }
        if (ncol(Z.hat.samples[[j]]) != B) {
            stop(sprintf("Feature %d has %d samples, but feature 1 has %d",
                         j, ncol(Z.hat.samples[[j]]), B))
        }
    }

    if (verbose) {
        message(sprintf("Computing lslope with posterior propagation:"))
        message(sprintf("  Features: %d, Vertices: %d, Samples: %d", p, n, B))
    }

    ## Storage for results
    lslope.mean <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.sd <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.lower <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.upper <- matrix(NA_real_, nrow = p, ncol = n)

    if (return.samples) {
        lslope.samples <- vector("list", p)
    }

    if (verbose) {
        pb <- txtProgressBar(min = 0, max = p, style = 3)
    }

    for (j in seq_len(p)) {

        ## Storage for this feature's samples: n x B
        samples.j <- matrix(NA_real_, nrow = n, ncol = B)

        ## Compute lslope for each posterior sample
        for (b in seq_len(B)) {
            z.hat.b <- Z.hat.samples[[j]][, b]

            ## lslope.gradient with instrumented=FALSE returns a vector directly
            lslope.coeffs <- lslope.gradient(
                adj.list, weight.list,
                y.hat, z.hat.b,
                type = lslope.type,
                ascending = ascending,
                instrumented = FALSE
            )

            samples.j[, b] <- lslope.coeffs
        }

        ## Compute summary statistics at each vertex
        lslope.mean[j, ] <- apply(samples.j, 1, mean, na.rm = TRUE)
        lslope.sd[j, ] <- apply(samples.j, 1, sd, na.rm = TRUE)
        lslope.lower[j, ] <- apply(samples.j, 1, quantile,
                                    probs = 0.025, na.rm = TRUE)
        lslope.upper[j, ] <- apply(samples.j, 1, quantile,
                                    probs = 0.975, na.rm = TRUE)

        if (return.samples) {
            lslope.samples[[j]] <- samples.j
        }

        if (verbose) {
            setTxtProgressBar(pb, j)
        }
    }

    if (verbose) {
        close(pb)
    }

    result <- list(
        mean = lslope.mean,
        sd = lslope.sd,
        lower = lslope.lower,
        upper = lslope.upper,
        n.samples = B,
        n.features = p,
        n.vertices = n,
        lslope.type = lslope.type
    )

    if (return.samples) {
        result$samples <- lslope.samples
    }

    class(result) <- c("lslope.posterior", "list")
    return(result)
}


#' @export
print.lslope.posterior <- function(x, ...) {
    cat("\nlslope with Posterior Uncertainty Propagation\n")
    cat("=============================================\n\n")

    cat(sprintf("Features: %d\n", x$n.features))
    cat(sprintf("Vertices: %d\n", x$n.vertices))
    cat(sprintf("Posterior samples: %d\n", x$n.samples))
    cat(sprintf("lslope type: %s\n", x$lslope.type))

    cat(sprintf("\nPosterior mean lslope (summary):\n"))
    cat(sprintf("  Global mean: %.4f\n", mean(x$mean, na.rm = TRUE)))
    cat(sprintf("  Global SD: %.4f\n", sd(as.vector(x$mean), na.rm = TRUE)))

    cat(sprintf("\nPosterior SD of lslope (summary):\n"))
    cat(sprintf("  Mean SD: %.4f\n", mean(x$sd, na.rm = TRUE)))
    cat(sprintf("  Median SD: %.4f\n", median(x$sd, na.rm = TRUE)))

    ## 95% CI width
    ci.width <- x$upper - x$lower
    cat(sprintf("\n95%% CI width (summary):\n"))
    cat(sprintf("  Mean width: %.4f\n", mean(ci.width, na.rm = TRUE)))
    cat(sprintf("  Median width: %.4f\n", median(ci.width, na.rm = TRUE)))

    invisible(x)
}


#' Fit Mixture Model with Posterior-Averaged lslope
#'
#' Extends the mixture model fitting to use posterior mean lslope values
#' and incorporate the posterior uncertainty into the inference.
#'
#' @param lslope.post Result from \code{lslope.with.posterior()}
#' @param null.dist Null distribution from \code{lslope.null.distribution()}
#' @param use.posterior.sd Logical. If TRUE, inflate null variance by
#'   posterior SD to account for smoothing uncertainty
#' @param pi0.method Method for estimating pi0
#' @param verbose Logical. Print progress
#'
#' @return A list of class "lslope.mixture" with posterior-adjusted results
#'
#' @export
fit.lslope.mixture.posterior <- function(lslope.post,
                                          null.dist,
                                          use.posterior.sd = TRUE,
                                          pi0.method = "convex",
                                          verbose = TRUE) {

    if (!inherits(lslope.post, "lslope.posterior")) {
        stop("lslope.post must be from lslope.with.posterior()")
    }

    ## Use posterior mean as the observed values
    obs.lslope <- lslope.post$mean

    ## If requested, inflate null SD by posterior uncertainty
    if (use.posterior.sd) {
        if (verbose) {
            message("Adjusting null distribution for posterior uncertainty...")
        }

        ## Adjust null SD: sqrt(null.sd^2 + posterior.sd^2)
        ## This accounts for the additional uncertainty from smoothing
        adjusted.null.params <- null.dist$null.params

        ## Mean posterior SD across features at each vertex
        mean.post.sd <- apply(lslope.post$sd, 2, mean, na.rm = TRUE)

        adjusted.null.params$sd <- sqrt(
            null.dist$null.params$sd^2 + mean.post.sd^2
        )

        ## Create adjusted null.dist object
        null.dist.adj <- null.dist
        null.dist.adj$null.params <- adjusted.null.params
        null.dist.adj$adjusted.for.posterior <- TRUE

    } else {
        null.dist.adj <- null.dist
    }

    ## Fit mixture model using adjusted null
    mixture <- fit.lslope.mixture(
        obs.lslope, null.dist.adj,
        pi0.method = pi0.method,
        null.model = "normal",
        verbose = verbose
    )

    ## Add posterior information
    mixture$lslope.posterior <- list(
        sd = lslope.post$sd,
        lower = lslope.post$lower,
        upper = lslope.post$upper,
        n.samples = lslope.post$n.samples
    )
    mixture$adjusted.for.posterior <- use.posterior.sd

    return(mixture)
}


## ============================================================================
## INTEGRATION WITH refit.rdgraph.regression()
## ============================================================================

#' Full Pipeline with Posterior Propagation
#'
#' Complete pipeline that:
#' 1. Computes smoothed features with posterior samples
#' 2. Propagates uncertainty into lslope estimates
#' 3. Generates null distribution (permuting original y, then smoothing)
#' 4. Fits mixture model with posterior adjustment
#'
#' @param fitted.model Fitted model from \code{fit.rdgraph.regression()}
#' @param y Original (unsmoothed) response vector. This is permuted to
#'   generate the null distribution.
#' @param Z Matrix of features to test (n x p)
#' @param n.posterior.samples Number of posterior samples for smoothing
#' @param n.perm Number of permutations for null distribution
#' @param lslope.type Type of lslope measure
#' @param batch.size Number of features to process per batch (for memory)
#' @param seed Random seed
#' @param verbose Print progress
#' @param n.cores Number of cores for parallel processing
#'
#' @return A comprehensive result object with posterior-propagated inference
#'
#' @details
#' This function orchestrates the full posterior-propagated inference pipeline.
#' Due to memory constraints (storing posterior samples for many features),
#' features are processed in batches when p is large.
#'
#' IMPORTANT: The null distribution is generated by permuting the ORIGINAL
#' response y (not the smoothed y.hat), then smoothing the permuted response.
#' This ensures the null properly captures what lslope values would look like
#' when there's no true signal but the geometric smoothing is still applied.
#'
#' The recommended workflow is:
#' 1. Fit the response model with \code{fit.rdgraph.regression()}
#' 2. Call this function with the original y and feature matrix Z
#' 3. Examine the results for significant context-dependent associations
#'
#' @examples
#' \dontrun{
#' ## Fit model
#' fit <- fit.rdgraph.regression(X, y, k = 15)
#'
#' ## Run full posterior-propagated pipeline
#' ## Note: pass original y, not fit$fitted.values!
#' result <- lslope.test.with.posterior(
#'     fitted.model = fit,
#'     y = y,  # Original response
#'     Z = Z,
#'     n.posterior.samples = 500,
#'     n.perm = 200,
#'     n.cores = 4
#' )
#'
#' ## Examine results
#' print(result$mixture)
#' summary(result$mixture, threshold = 0.9)
#' }
#'
#' @export
lslope.test.with.posterior <- function(fitted.model,
                                        y,
                                        Z,
                                        n.posterior.samples = 500L,
                                        n.perm = 200L,
                                        lslope.type = "normalized",
                                        batch.size = 50L,
                                        seed = 12345L,
                                        verbose = TRUE,
                                        n.cores = 1L) {

    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be from fit.rdgraph.regression()")
    }

    ## Extract components
    adj.list <- fitted.model$graph$adj.list
    weight.list <- fitted.model$graph$edge.length.list
    y.hat <- fitted.model$fitted.values
    n <- length(y.hat)

    ## Validate original y
    if (length(y) != n) {
        stop("length(y) must equal length(fitted.model$fitted.values)")
    }

    if (is.vector(Z)) {
        Z <- matrix(Z, ncol = 1)
    }
    if (nrow(Z) != n) {
        stop("nrow(Z) must equal length(y)")
    }
    p <- ncol(Z)

    if (verbose) {
        message("===========================================")
        message("lslope Association Test with Posterior Propagation")
        message("===========================================\n")
        message(sprintf("Features: %d, Vertices: %d", p, n))
        message(sprintf("Posterior samples: %d, Permutations: %d\n",
                        n.posterior.samples, n.perm))
    }

    ## ========================================================================
    ## Step 1: Smooth features
    ## ========================================================================

    if (verbose) message("\n=== Step 1: Computing smoothed features ===\n")

    Z.hat.point <- refit.rdgraph.regression(
        fitted.model, Z,
        per.column.gcv = TRUE,
        verbose = verbose
    )

    ## ========================================================================
    ## Step 2: Generate null distribution (permuting original y!)
    ## ========================================================================

    if (verbose) message("\n=== Step 2: Generating null distribution ===\n")

    null.dist <- lslope.null.distribution(
        fitted.model = fitted.model,
        y = y,  # Original y, not y.hat!
        Z.hat = Z.hat.point$fitted.values,
        n.perm = n.perm,
        lslope.type = lslope.type,
        per.column.gcv = FALSE,  # Use fixed eta for efficiency
        seed = seed,
        verbose = verbose,
        n.cores = n.cores
    )

    ## ========================================================================
    ## Step 2: Compute lslope with posterior propagation (in batches)
    ## ========================================================================

    if (verbose) message("\n=== Step 3: Posterior-propagated lslope ===\n")

    ## Initialize storage for aggregated results
    lslope.mean.all <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.sd.all <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.lower.all <- matrix(NA_real_, nrow = p, ncol = n)
    lslope.upper.all <- matrix(NA_real_, nrow = p, ncol = n)

    if (!is.null(colnames(Z))) {
        rownames(lslope.mean.all) <- colnames(Z)
        rownames(lslope.sd.all) <- colnames(Z)
        rownames(lslope.lower.all) <- colnames(Z)
        rownames(lslope.upper.all) <- colnames(Z)
    }

    ## Process in batches to manage memory
    n.batches <- ceiling(p / batch.size)

    for (batch in seq_len(n.batches)) {
        batch.start <- (batch - 1) * batch.size + 1
        batch.end <- min(batch * batch.size, p)
        batch.idx <- batch.start:batch.end

        if (verbose) {
            message(sprintf("  Batch %d/%d: features %d-%d",
                            batch, n.batches, batch.start, batch.end))
        }

        ## Get posterior samples for this batch
        Z.batch <- Z[, batch.idx, drop = FALSE]

        Z.hat.post <- refit.rdgraph.regression(
            fitted.model, Z.batch,
            per.column.gcv = TRUE,
            with.posterior = TRUE,
            return.posterior.samples = TRUE,
            n.posterior.samples = n.posterior.samples,
            posterior.seed = seed + batch,
            verbose = FALSE
        )

        ## Compute lslope with posterior propagation
        lslope.post.batch <- lslope.with.posterior(
            adj.list, weight.list,
            y.hat, Z.hat.post$posterior$samples,
            lslope.type = lslope.type,
            return.samples = FALSE,
            verbose = FALSE
        )

        ## Store results
        lslope.mean.all[batch.idx, ] <- lslope.post.batch$mean
        lslope.sd.all[batch.idx, ] <- lslope.post.batch$sd
        lslope.lower.all[batch.idx, ] <- lslope.post.batch$lower
        lslope.upper.all[batch.idx, ] <- lslope.post.batch$upper
    }

    ## Create lslope.posterior object
    lslope.post <- list(
        mean = lslope.mean.all,
        sd = lslope.sd.all,
        lower = lslope.lower.all,
        upper = lslope.upper.all,
        n.samples = n.posterior.samples,
        n.features = p,
        n.vertices = n,
        lslope.type = lslope.type
    )
    class(lslope.post) <- c("lslope.posterior", "list")

    ## ========================================================================
    ## Step 3: Fit mixture model with posterior adjustment
    ## ========================================================================

    if (verbose) message("\n=== Step 4: Fitting mixture model ===\n")

    mixture <- fit.lslope.mixture.posterior(
        lslope.post, null.dist,
        use.posterior.sd = TRUE,
        pi0.method = "convex",
        verbose = verbose
    )

    ## ========================================================================
    ## Assemble results
    ## ========================================================================

    result <- list(
        lslope.posterior = lslope.post,
        null.dist = null.dist,
        mixture = mixture,
        Z.hat.point = Z.hat.point$fitted.values,
        parameters = list(
            n.posterior.samples = n.posterior.samples,
            n.perm = n.perm,
            lslope.type = lslope.type,
            batch.size = batch.size,
            seed = seed
        )
    )

    class(result) <- c("lslope.posterior.test", "list")

    if (verbose) {
        message("\n=== Pipeline complete ===\n")
        print(mixture)
    }

    return(result)
}


#' @export
print.lslope.posterior.test <- function(x, ...) {
    cat("\n================================================\n")
    cat("Bayesian lslope Test with Posterior Propagation\n")
    cat("================================================\n\n")

    cat(sprintf("Features: %d\n", x$lslope.posterior$n.features))
    cat(sprintf("Vertices: %d\n", x$lslope.posterior$n.vertices))
    cat(sprintf("Posterior samples: %d\n", x$parameters$n.posterior.samples))
    cat(sprintf("Permutations: %d\n", x$parameters$n.perm))
    cat(sprintf("lslope type: %s\n\n", x$parameters$lslope.type))

    print(x$mixture)

    invisible(x)
}


## ============================================================================
## VISUALIZATION HELPERS
## ============================================================================

#' Plot Posterior Distribution of lslope at a Vertex
#'
#' Visualizes the posterior distribution of lslope values for a specific
#' feature at a specific vertex, showing the null distribution and
#' posterior probability of association.
#'
#' @param lslope.samples Matrix (n x B) of posterior samples for one feature
#' @param vertex Vertex index
#' @param null.params Null distribution parameters from lslope.null.distribution()
#' @param prob.alt Posterior probability of alternative (from mixture model)
#' @param feature.name Optional feature name for plot title
#' @param ... Additional arguments passed to plot
#'
#' @export
plot.lslope.posterior.vertex <- function(lslope.samples,
                                          vertex,
                                          null.params,
                                          prob.alt = NULL,
                                          feature.name = NULL,
                                          ...) {

    samples.v <- lslope.samples[vertex, ]
    samples.v <- samples.v[is.finite(samples.v)]

    if (length(samples.v) < 10) {
        warning("Too few valid samples at this vertex")
        return(invisible(NULL))
    }

    ## Compute kernel density for posterior
    post.dens <- density(samples.v, bw = "SJ")

    ## Null distribution
    null.mean <- null.params$mean[vertex]
    null.sd <- null.params$sd[vertex]
    null.x <- seq(min(post.dens$x), max(post.dens$x), length.out = 200)
    null.y <- dnorm(null.x, mean = null.mean, sd = null.sd)

    ## Scale null to match posterior density
    scale.factor <- max(post.dens$y) / max(null.y)
    null.y.scaled <- null.y * scale.factor * 0.8

    ## Plot
    plot(post.dens,
         main = if (!is.null(feature.name)) {
             sprintf("Posterior lslope: %s at vertex %d", feature.name, vertex)
         } else {
             sprintf("Posterior lslope at vertex %d", vertex)
         },
         xlab = "lslope",
         ylab = "Density",
         col = "steelblue",
         lwd = 2,
         ...)

    ## Add null distribution
    lines(null.x, null.y.scaled, col = "gray50", lty = 2, lwd = 1.5)

    ## Add credible interval
    ci <- quantile(samples.v, probs = c(0.025, 0.975))
    abline(v = ci, col = "steelblue", lty = 3)

    ## Add posterior mean
    abline(v = mean(samples.v), col = "steelblue", lwd = 2)

    ## Add legend
    legend.text <- c("Posterior", "Null",
                     sprintf("95%% CI: [%.3f, %.3f]", ci[1], ci[2]))
    legend.col <- c("steelblue", "gray50", "steelblue")
    legend.lty <- c(1, 2, 3)

    if (!is.null(prob.alt)) {
        legend.text <- c(legend.text, sprintf("P(assoc) = %.3f", prob.alt))
        legend.col <- c(legend.col, "black")
        legend.lty <- c(legend.lty, NA)
    }

    legend("topright",
           legend = legend.text,
           col = legend.col,
           lty = legend.lty,
           lwd = c(2, 1.5, 1, NA),
           bty = "n")

    invisible(list(
        posterior.density = post.dens,
        null.x = null.x,
        null.y = null.y,
        ci = ci
    ))
}

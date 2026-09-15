## ============================================================================
## lcor with Posterior Uncertainty Propagation
## ============================================================================
##
## This function computes lcor() statistics with full posterior uncertainty
## quantification represented by supplied samples of vertex fields. Sampling
## and smoothing belong to the model that produced those fields; this file
## contains only the model-independent local-association summary.
##
## ============================================================================

#' Compute Local Correlation with Posterior Uncertainty Propagation
#'
#' Summarizes local-correlation estimates over supplied draws of vertex fields,
#' providing a mean, standard deviation, and interval at each vertex for each
#' feature. The function does not fit or refit a smoother.
#'
#' @param adj.list Adjacency list (1-based R indexing)
#' @param weight.list Edge weight list
#' @param y.hat Finite numeric response vector of positive length n, in graph
#'   vertex order. Optional names identify vertices; see Details.
#' @param Z.hat.samples A list of length p whose elements are n by B matrices
#'   of finite numeric sampled vertex fields, or one n by B matrix for a single
#'   feature. Every feature must have the same positive number B of draws.
#'   Optional list names identify features and matrix row names identify vertices.
#' @param lcor.type Type of local correlation weighting: "derivative" (default),
#'   "unit", or "sign". See \code{\link{lcor}} for details.
#' @param credible.level Finite numeric interval level strictly between 0 and 1
#'   (default 0.95).
#' @param return.samples Logical. Return individual lcor samples (default FALSE).
#'   Setting TRUE with many features will use substantial memory.
#' @param verbose Logical. Print progress (default TRUE)
#'
#' @section Memory:
#' The four feature-by-vertex summary matrices require about 32*p*n bytes.
#' Processing one feature also materializes an n-by-B double matrix (8*n*B
#' bytes), even when `return.samples = FALSE`. Returning samples retains another
#' 8*p*n*B bytes across features, in addition to input draws and temporaries.
#' Process subsets of features separately when needed; keep the same response,
#' graph, settings and complete draw count, and preserve feature identities when
#' combining summaries. This does not reduce the per-feature draw allocation.
#'
#' @return A list of class "lcor.posterior" containing:
#'   \describe{
#'     \item{mean}{Matrix (p x n) of posterior mean lcor values}
#'     \item{sd}{Matrix (p x n) of posterior standard deviations}
#'     \item{lower}{Matrix (p x n) of lower credible bounds}
#'     \item{upper}{Matrix (p x n) of upper credible bounds}
#'     \item{samples}{List of sample matrices (only if return.samples = TRUE)}
#'     \item{n.samples}{Number of posterior samples used}
#'     \item{n.features}{Number of features}
#'     \item{n.vertices}{Number of vertices}
#'     \item{lcor.type}{Type of lcor weighting used}
#'   }
#'
#' @details
#' Each supplied draw is passed to \code{\link{lcor}} against the fixed
#' response field. This makes the estimand independent of how the draws were
#' generated. The graph and response remain fixed. Quantile intervals are
#' pointwise, not simultaneous; their posterior interpretation depends on the
#' supplied draws and excludes uncertainty not represented in those draws.
#' Archived graph-regression objects can be adapted in \code{gflowx}.
#'
#' All features are validated before any local correlation is computed; unequal
#' draw counts are rejected rather than truncated. With one draw, the mean and
#' both interval endpoints equal its local correlation, and the standard
#' deviation is \code{NA}.
#'
#' Supplied feature names and vertex names must be unique, nonmissing, and
#' nonempty. Vertex names in \code{y.hat} and row names in any supplied matrices
#' must agree in order; inputs are never reordered by name. Names do not verify
#' graph alignment: rows must still follow the adjacency list's vertex order.
#' Feature names become row names of the summary matrices and names of the
#' optional samples list. Vertex names become column names of the summaries and
#' row names of sample matrices; draw column names are preserved separately for
#' each feature. Unnamed dimensions remain unnamed, and \code{summary()} labels
#' unnamed features \code{Feature1}, \code{Feature2}, and so on.
#'
#' @examples
#' adj <- list(c(2L), c(1L, 3L), c(2L))
#' weights <- lapply(adj, function(x) rep(1, length(x)))
#' draws <- cbind(c(0, 1, 2), c(0.1, 0.9, 2.1), c(-0.1, 1.1, 1.9))
#' post <- lcor.with.posterior(adj, weights, c(0, 1, 2), draws,
#'                             verbose = FALSE)
#'
#' @seealso \code{\link{lcor}} for single-sample local correlation computation
#'
#' @export
lcor.with.posterior <- function(adj.list,
                                 weight.list,
                                 y.hat,
                                 Z.hat.samples,
                                 lcor.type = c("derivative", "unit", "sign"),
                                 credible.level = 0.95,
                                 return.samples = FALSE,
                                 verbose = TRUE) {

    lcor.type <- match.arg(lcor.type)
    lcor.with.posterior.R(
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

    if (!is.numeric(y.hat) || !is.null(dim(y.hat)) ||
        length(y.hat) == 0L || any(!is.finite(y.hat))) {
        stop("y.hat must be a nonempty finite numeric vector.", call. = FALSE)
    }
    n <- length(y.hat)
    if (!is.numeric(credible.level) || length(credible.level) != 1L ||
        !is.finite(credible.level) || credible.level <= 0 || credible.level >= 1) {
        stop("credible.level must be a finite number strictly between 0 and 1.",
             call. = FALSE)
    }
    for (control in c("return.samples", "verbose")) {
        value <- get(control)
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop(control, " must be TRUE or FALSE.", call. = FALSE)
        }
    }

    ## Handle single feature case
    if (is.matrix(Z.hat.samples)) {
        Z.hat.samples <- list(Z.hat.samples)
    }

    if (!is.list(Z.hat.samples) || is.data.frame(Z.hat.samples) ||
        !is.null(dim(Z.hat.samples)) || length(Z.hat.samples) == 0L) {
        stop("Z.hat.samples must be a numeric matrix or a nonempty list of numeric matrices.",
             call. = FALSE)
    }
    p <- length(Z.hat.samples)
    validate.names <- function(value, label) {
        if (!is.null(value) &&
            (anyNA(value) || any(!nzchar(value)) || anyDuplicated(value))) {
            stop(label, " must be unique, nonmissing, and nonempty.", call. = FALSE)
        }
    }
    feature.names <- names(Z.hat.samples)
    validate.names(feature.names, "names(Z.hat.samples)")
    vertex.ids <- names(y.hat)
    validate.names(vertex.ids, "names(y.hat)")
    n.samples <- NULL
    ## Validate every feature before allocating results or calling lcor().
    for (j in seq_len(p)) {
        label <- sprintf("Z.hat.samples[[%d]]", j)
        if (!is.null(feature.names)) {
            label <- paste0(label, " ('", feature.names[[j]], "')")
        }
        samples.j <- Z.hat.samples[[j]]
        if (!is.matrix(samples.j) || !is.numeric(samples.j)) {
            stop(label, " must be a numeric matrix.", call. = FALSE)
        }
        if (nrow(samples.j) != n) {
            stop(sprintf("%s has %d rows; expected %d to match length(y.hat).",
                         label, nrow(samples.j), n), call. = FALSE)
        }
        if (ncol(samples.j) == 0L) {
            stop(label, " must contain at least one draw (column).", call. = FALSE)
        }
        if (is.null(n.samples)) {
            n.samples <- ncol(samples.j)
        } else if (ncol(samples.j) != n.samples) {
            stop(sprintf("%s has %d columns; expected %d draws as in Z.hat.samples[[1]].",
                         label, ncol(samples.j), n.samples), call. = FALSE)
        }
        if (any(!is.finite(samples.j))) {
            stop(label, " must contain only finite values.", call. = FALSE)
        }
        ids <- rownames(samples.j)
        validate.names(ids, paste0("rownames(", label, ")"))
        if (!is.null(ids)) {
            if (is.null(vertex.ids)) {
                vertex.ids <- ids
            } else if (!identical(ids, vertex.ids)) {
                stop(label, " must use the same vertex names in the same order as y.hat or the preceding named feature.",
                     call. = FALSE)
            }
        }
    }

    if (verbose) {
        message(sprintf("lcor with posterior (R mode): %d features, %d samples",
                        p, n.samples))
    }

    ## Initialize output storage
    summary.names <- if (is.null(feature.names) && is.null(vertex.ids)) NULL else
        list(feature.names, vertex.ids)
    lcor.mean <- matrix(NA_real_, nrow = p, ncol = n, dimnames = summary.names)
    lcor.sd <- matrix(NA_real_, nrow = p, ncol = n, dimnames = summary.names)
    lcor.lower <- matrix(NA_real_, nrow = p, ncol = n, dimnames = summary.names)
    lcor.upper <- matrix(NA_real_, nrow = p, ncol = n, dimnames = summary.names)

    if (return.samples) {
        all.samples <- vector("list", p)
        names(all.samples) <- feature.names
    }

    alpha.lower <- (1 - credible.level) / 2
    alpha.upper <- (1 + credible.level) / 2

    if (verbose && p > 5) {
        pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
    }

    for (j in seq_len(p)) {
        ## Z.hat.samples[[j]] is n x n.samples matrix
        samples.j <- Z.hat.samples[[j]]

        ## Compute lcor for each posterior sample
        lcor.samples.j <- matrix(NA_real_, nrow = n, ncol = n.samples)
        if (!is.null(vertex.ids) || !is.null(colnames(samples.j))) {
            dimnames(lcor.samples.j) <- list(vertex.ids, colnames(samples.j))
        }

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
            utils::setTxtProgressBar(pb, j)
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


#' Print a Local-Correlation Posterior Result
#'
#' @param x An object returned by \code{\link{lcor.with.posterior}}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.lcor.posterior <- function(x, ...) {
    .validate.lcor.posterior(x)
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


#' Summarize a Local-Correlation Posterior Result
#'
#' @param object An object returned by \code{\link{lcor.with.posterior}}.
#' @param threshold Reserved for compatibility with earlier summary calls.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A \code{summary.lcor.posterior} object containing per-feature and
#'   overall credible-interval summaries.
#'
#' @export
summary.lcor.posterior <- function(object, threshold = 0.5, ...) {
    .validate.lcor.posterior(object)

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


#' Print a Local-Correlation Posterior Summary
#'
#' @param x An object returned by \code{summary.lcor.posterior()}.
#' @param n.top Maximum number of features to print.
#' @param ... Additional arguments, currently ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.summary.lcor.posterior <- function(x, n.top = 20, ...) {
    .validate.summary.lcor.posterior(x)
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

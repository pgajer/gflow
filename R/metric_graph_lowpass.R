#' Construct a Metric-Conductance Graph Low-Pass Operator
#'
#' Builds a weighted graph Laplacian by transforming metric edge lengths into
#' conductances. This operator is a direct metric-conductance comparator to
#' \code{\link{fit.rdgraph.regression}}, not a replacement for the
#' Riemannian-complex/overlap-density smoother.
#'
#' @param adj.list List of integer neighbor vectors using 1-based vertex indices.
#' @param weight.list List of positive metric edge lengths parallel to
#'   \code{adj.list}. These are interpreted as edge lengths, not conductances.
#' @param conductance.rule Character scalar. One of
#'   \code{"inverse.length.power"}, \code{"exp.length"},
#'   \code{"exp.length.squared"}, or \code{"self.tuned.gaussian"}.
#' @param conductance.epsilon Positive numeric regularizer used by inverse-power
#'   conductances and as a local-scale floor.
#' @param conductance.alpha Positive numeric exponent for
#'   \code{"inverse.length.power"}.
#' @param conductance.sigma Optional positive global scale for exponential rules.
#'   If \code{NULL}, it is selected by \code{conductance.sigma.rule}.
#' @param conductance.sigma.rule Rule for selecting a global scale when needed.
#' @param conductance.sigma.quantile Quantile used when
#'   \code{conductance.sigma.rule = "edge.quantile"}.
#' @param conductance.local.k Positive integer local incident-edge order
#'   statistic for \code{"self.tuned.gaussian"}.
#' @param laplacian.type Phase-1 supports only \code{"unnormalized"}, the
#'   weighted graph Laplacian \eqn{L = D - C}.
#' @param return.sparse Logical. If \code{TRUE}, attach a \pkg{Matrix} sparse
#'   Laplacian when \pkg{Matrix} is available.
#' @param verbose Logical. Reserved for future diagnostic messages.
#'
#' @details
#' The current \code{fit.rdgraph.regression()} precomputed-graph path uses
#' supplied \code{weight.list} values as edge lengths for neighborhood ordering.
#' The spectral conductance in that smoother is overlap-density based:
#' \deqn{c_e^\rho = 1 / \max(\rho_1(e), 10^{-10}),}
#' where \eqn{\rho_1(e)} is computed by the Riemannian-complex density machinery.
#'
#' This function instead constructs conductances directly from metric lengths:
#' \deqn{c_{ij} = \phi(\ell_{ij}).}
#' Supported phase-1 transforms are
#' \deqn{c_{ij}=(\ell_{ij}+\epsilon)^{-\alpha},}
#' \deqn{c_{ij}=\exp(-\ell_{ij}/\sigma),}
#' \deqn{c_{ij}=\exp(-\ell_{ij}^{2}/\sigma^{2}),}
#' and
#' \deqn{c_{ij}=\exp(-\ell_{ij}^{2}/(\sigma_i\sigma_j)).}
#'
#' @return A list of class \code{"metric.graph.lowpass.operator"} containing
#'   the edge table, conductances, degree vector, Laplacian triplets, summaries,
#'   and optionally a sparse Laplacian matrix.
#' @export
metric.graph.lowpass.operator <- function(
    adj.list,
    weight.list,
    conductance.rule = c("inverse.length.power", "exp.length",
                         "exp.length.squared", "self.tuned.gaussian"),
    conductance.epsilon = 1e-8,
    conductance.alpha = 1,
    conductance.sigma = NULL,
    conductance.sigma.rule = c("edge.quantile", "median", "local.k"),
    conductance.sigma.quantile = 0.75,
    conductance.local.k = 5L,
    laplacian.type = c("unnormalized"),
    return.sparse = TRUE,
    verbose = FALSE) {

    graph <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    args <- .validate.metric.graph.lowpass.operator.args(
        conductance.rule = conductance.rule,
        conductance.epsilon = conductance.epsilon,
        conductance.alpha = conductance.alpha,
        conductance.sigma = conductance.sigma,
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = conductance.sigma.quantile,
        conductance.local.k = conductance.local.k,
        laplacian.type = laplacian.type,
        verbose = verbose
    )

    out <- .Call(
        "S_metric_graph_lowpass_operator",
        graph$adj.list.0based,
        graph$weight.list.cpp,
        args$conductance.rule,
        args$conductance.epsilon,
        args$conductance.alpha,
        args$conductance.sigma,
        args$conductance.sigma.rule,
        args$conductance.sigma.quantile,
        args$conductance.local.k,
        args$laplacian.type,
        PACKAGE = "gflow"
    )

    out$graph$adj.list <- graph$adj.list
    out$graph$weight.list <- graph$weight.list
    out <- .attach.metric.graph.lowpass.laplacian(out, return.sparse)
    class(out) <- c("metric.graph.lowpass.operator", "list")
    out
}

#' Fit Metric-Conductance Graph Low-Pass Regression
#'
#' Fits graph-spectral low-pass regression on a supplied graph by transforming
#' metric edge lengths into conductances and smoothing the response in the
#' eigenbasis of the weighted graph Laplacian.
#'
#' @inheritParams metric.graph.lowpass.operator
#' @param y Numeric response vector of length \code{length(adj.list)}.
#' @param n.eigenpairs Positive integer number of eigenpairs to compute.
#' @param filter.type Spectral low-pass filter family.
#' @param eta.grid Optional positive numeric eta grid. If \code{NULL}, the
#'   existing package helper \code{generate.eta.grid()} is used.
#' @param n.candidates Number of eta candidates when \code{eta.grid = NULL}.
#' @param eigen.solver \code{"auto"}, \code{"sparse"}, or \code{"dense"}.
#'   \code{"auto"} uses dense decomposition only for
#'   \code{n <= dense.eigen.threshold}, then sparse-first.
#' @param dense.eigen.threshold Exact dense threshold for auto mode. Default
#'   \code{200L}, intended for small reference/testing problems.
#' @param dense.fallback.threshold Maximum graph size for emergency dense
#'   fallback when sparse decomposition fails and fallback is allowed.
#' @param dense.fallback \code{"auto"}, \code{"never"}, or \code{"always"}.
#'
#' @return A list of class \code{"metric.graph.lowpass.fit"}.
#' @export
fit.metric.graph.lowpass <- function(
    adj.list,
    weight.list,
    y,
    conductance.rule = c("inverse.length.power", "exp.length",
                         "exp.length.squared", "self.tuned.gaussian"),
    conductance.epsilon = 1e-8,
    conductance.alpha = 1,
    conductance.sigma = NULL,
    conductance.sigma.rule = c("edge.quantile", "median", "local.k"),
    conductance.sigma.quantile = 0.75,
    conductance.local.k = 5L,
    laplacian.type = c("unnormalized"),
    n.eigenpairs = 50L,
    filter.type = c("heat_kernel", "tikhonov", "cubic_spline",
                    "gaussian", "exponential", "butterworth"),
    eta.grid = NULL,
    n.candidates = 40L,
    eigen.solver = c("auto", "sparse", "dense"),
    dense.eigen.threshold = 200L,
    dense.fallback.threshold = 5000L,
    dense.fallback = c("auto", "never", "always"),
    verbose = FALSE) {

    graph <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    n <- length(graph$adj.list)
    y <- .validate.metric.graph.lowpass.response(y, n, "y")

    args <- .validate.metric.graph.lowpass.operator.args(
        conductance.rule = conductance.rule,
        conductance.epsilon = conductance.epsilon,
        conductance.alpha = conductance.alpha,
        conductance.sigma = conductance.sigma,
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = conductance.sigma.quantile,
        conductance.local.k = conductance.local.k,
        laplacian.type = laplacian.type,
        verbose = verbose
    )
    solver <- .validate.metric.graph.lowpass.solver.args(
        n.eigenpairs = n.eigenpairs,
        n = n,
        eigen.solver = eigen.solver,
        dense.eigen.threshold = dense.eigen.threshold,
        dense.fallback.threshold = dense.fallback.threshold,
        dense.fallback = dense.fallback
    )
    filter.type <- match.arg(filter.type)
    n.candidates <- .validate.positive.integer.scalar(n.candidates, "n.candidates")

    raw <- .Call(
        "S_metric_graph_lowpass_spectrum",
        graph$adj.list.0based,
        graph$weight.list.cpp,
        args$conductance.rule,
        args$conductance.epsilon,
        args$conductance.alpha,
        args$conductance.sigma,
        args$conductance.sigma.rule,
        args$conductance.sigma.quantile,
        args$conductance.local.k,
        args$laplacian.type,
        solver$n.eigenpairs,
        solver$eigen.solver,
        solver$dense.eigen.threshold,
        solver$dense.fallback.threshold,
        solver$dense.fallback,
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    operator <- raw$operator
    operator$graph$adj.list <- graph$adj.list
    operator$graph$weight.list <- graph$weight.list
    operator <- .attach.metric.graph.lowpass.laplacian(operator, TRUE)
    class(operator) <- c("metric.graph.lowpass.operator", "list")

    spectral <- raw$spectral
    eigenvalues <- as.double(spectral$eigenvalues)
    V <- as.matrix(spectral$eigenvectors)

    eta.grid <- .prepare.metric.graph.lowpass.eta.grid(
        eta.grid = eta.grid,
        eigenvalues = eigenvalues,
        filter.type = filter.type,
        n.candidates = n.candidates
    )
    filter.weights.matrix <- compute.filter.weights.matrix(eigenvalues, eta.grid, filter.type)
    y.spectral <- as.vector(crossprod(V, y))
    gcv.result <- select.eta.gcv.single(y, y.spectral, V, filter.weights.matrix, eta.grid)
    best.idx <- gcv.result$best.idx

    spectral$filtered.eigenvalues <- filter.weights.matrix[, best.idx]
    spectral$eta.optimal <- gcv.result$eta.optimal
    spectral$filter.type <- filter.type
    spectral$n.eigenpairs <- length(eigenvalues)

    fitted.grid <- V %*% (y.spectral * filter.weights.matrix)
    rss.grid <- colSums((y - fitted.grid)^2)
    df.grid <- colSums(filter.weights.matrix)
    gcv.scores <- rss.grid / pmax(length(y) - df.grid, 1e-10)^2

    residuals <- y - gcv.result$y.hat
    result <- list(
        fitted.values = as.vector(gcv.result$y.hat),
        residuals = as.vector(residuals),
        y = y,
        graph = operator$graph,
        operator = operator,
        conductance = operator$conductance,
        laplacian = operator$laplacian,
        laplacian.type = operator$laplacian.type,
        spectral = spectral,
        gcv = list(
            eta.grid = eta.grid,
            gcv.scores = as.vector(gcv.scores),
            eta.optimal = gcv.result$eta.optimal,
            gcv.optimal = gcv.result$gcv.min,
            effective.df = gcv.result$effective.df,
            best.idx = best.idx
        ),
        parameters = c(args, solver, list(filter.type = filter.type)),
        timing = NULL
    )
    attr(result, "call") <- match.call()
    class(result) <- c("metric.graph.lowpass.fit", "list")
    result
}

#' Refit Metric-Conductance Graph Low-Pass Regression
#'
#' Reuses a fitted metric graph low-pass eigensystem to smooth new responses.
#'
#' @param fitted.model A \code{"metric.graph.lowpass.fit"} object.
#' @param y.new Numeric vector or matrix with one row per graph vertex.
#' @param per.column.gcv Logical. If \code{TRUE}, select eta independently for
#'   each response column using the cached eigenbasis.
#' @param eta.grid Optional positive numeric eta grid for per-column GCV.
#' @param n.candidates Number of eta candidates when \code{eta.grid = NULL}.
#' @param n.cores Number of cores for per-column GCV. Phase 1 uses sequential
#'   processing if optional parallel packages are unavailable.
#' @param block.size Optional block size for fixed-eta multi-column refits.
#' @param verbose Logical progress flag.
#'
#' @return A list of class \code{"metric.graph.lowpass.refit"}.
#' @export
refit.metric.graph.lowpass <- function(fitted.model,
                                       y.new,
                                       per.column.gcv = FALSE,
                                       eta.grid = NULL,
                                       n.candidates = 40L,
                                       n.cores = 1L,
                                       block.size = NULL,
                                       verbose = FALSE) {
    if (!inherits(fitted.model, "metric.graph.lowpass.fit")) {
        stop("fitted.model must be a 'metric.graph.lowpass.fit' object.")
    }
    spectral <- fitted.model$spectral
    V <- spectral$eigenvectors
    eigenvalues <- spectral$eigenvalues
    if (is.null(V) || !is.matrix(V)) stop("fitted.model$spectral$eigenvectors must be a matrix.")
    if (is.null(eigenvalues) || !is.numeric(eigenvalues)) {
        stop("fitted.model$spectral$eigenvalues must be numeric.")
    }

    n <- nrow(V)
    y.info <- .prepare.metric.graph.lowpass.response.matrix(y.new, n)
    Y <- y.info$Y
    n.responses <- ncol(Y)
    col.names <- y.info$col.names

    n.candidates <- .validate.positive.integer.scalar(n.candidates, "n.candidates")
    n.cores <- .validate.positive.integer.scalar(n.cores, "n.cores")
    block.size <- .validate.optional.block.size(block.size)
    filter.type <- spectral$filter.type
    if (is.null(filter.type)) filter.type <- "heat_kernel"

    if (isTRUE(per.column.gcv)) {
        eta.grid <- .prepare.metric.graph.lowpass.eta.grid(
            eta.grid = eta.grid,
            eigenvalues = eigenvalues,
            filter.type = filter.type,
            n.candidates = n.candidates
        )
        filter.weights.matrix <- compute.filter.weights.matrix(eigenvalues, eta.grid, filter.type)
        trace.S.all <- colSums(filter.weights.matrix)
        Vt.Y <- crossprod(V, Y)

        Y.hat <- matrix(0, nrow = n, ncol = n.responses)
        eta.optimal <- numeric(n.responses)
        gcv.scores <- numeric(n.responses)
        effective.df <- numeric(n.responses)
        best.idx <- integer(n.responses)

        if (verbose && n.responses > 1L) {
            message(sprintf("Selecting eta via GCV for %d response(s).", n.responses))
        }
        for (j in seq_len(n.responses)) {
            gcv.result <- select.eta.gcv.single(
                Y[, j], Vt.Y[, j], V, filter.weights.matrix, eta.grid
            )
            Y.hat[, j] <- gcv.result$y.hat
            eta.optimal[j] <- gcv.result$eta.optimal
            gcv.scores[j] <- gcv.result$gcv.min
            effective.df[j] <- gcv.result$effective.df
            best.idx[j] <- gcv.result$best.idx
        }
        if (!is.null(col.names)) {
            colnames(Y.hat) <- col.names
            names(eta.optimal) <- col.names
            names(gcv.scores) <- col.names
            names(effective.df) <- col.names
        }
        residuals <- Y - Y.hat
        out <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = TRUE,
            filter.type = filter.type,
            eta.optimal = if (n.responses == 1L) eta.optimal[1] else eta.optimal,
            gcv.scores = if (n.responses == 1L) gcv.scores[1] else gcv.scores,
            effective.df = if (n.responses == 1L) effective.df[1] else effective.df,
            eta.grid = eta.grid,
            best.idx = if (n.responses == 1L) best.idx[1] else best.idx,
            n.cores.used = 1L
        )
        class(out) <- c("metric.graph.lowpass.refit", "list")
        return(out)
    }

    f.lambda <- spectral$filtered.eigenvalues
    eta.used <- spectral$eta.optimal
    if (is.null(f.lambda) || length(f.lambda) != ncol(V)) {
        stop("fitted.model$spectral$filtered.eigenvalues is missing or has the wrong length.")
    }
    block.index <- .make.metric.graph.lowpass.block.index(n.responses, block.size)
    Y.hat <- matrix(0, nrow = n, ncol = n.responses)
    residuals <- matrix(0, nrow = n, ncol = n.responses)

    for (cols in block.index) {
        y.block <- Y[, cols, drop = FALSE]
        Vt.Y.block <- crossprod(V, y.block)
        Y.hat.block <- V %*% (f.lambda * Vt.Y.block)
        Y.hat[, cols] <- Y.hat.block
        residuals[, cols] <- y.block - Y.hat.block
    }
    if (!is.null(col.names)) {
        colnames(Y.hat) <- col.names
        colnames(residuals) <- col.names
    }
    out <- list(
        fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
        residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
        n.responses = n.responses,
        per.column.gcv = FALSE,
        eta.used = eta.used,
        block.size.used = if (length(block.index) > 1L) block.size else NA_integer_
    )
    class(out) <- c("metric.graph.lowpass.refit", "list")
    out
}

.validate.metric.graph.lowpass.graph <- function(adj.list, weight.list) {
    if (!is.list(adj.list)) stop("adj.list must be a list.")
    if (!is.list(weight.list)) stop("weight.list must be a list.")
    n <- length(adj.list)
    if (n < 2L) stop("adj.list must contain at least two vertices.")
    if (length(weight.list) != n) stop("adj.list and weight.list must have the same length.")

    adj.norm <- vector("list", n)
    weight.norm <- vector("list", n)
    tol <- 1e-10
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        wts <- weight.list[[i]]
        if (!is.numeric(nbrs) && !is.integer(nbrs)) stop(sprintf("adj.list[[%d]] must be numeric/integer.", i))
        if (!is.numeric(wts)) stop(sprintf("weight.list[[%d]] must be numeric.", i))
        nbrs <- as.integer(nbrs)
        wts <- as.double(wts)
        if (length(nbrs) != length(wts)) {
            stop(sprintf("Length mismatch at vertex %d: adj.list=%d, weight.list=%d.",
                         i, length(nbrs), length(wts)))
        }
        if (anyNA(nbrs)) stop(sprintf("adj.list[[%d]] contains NA.", i))
        if (any(nbrs < 1L | nbrs > n)) stop(sprintf("adj.list[[%d]] has indices outside 1..n.", i))
        if (any(nbrs == i)) stop(sprintf("adj.list[[%d]] contains self-loops.", i))
        if (anyDuplicated(nbrs)) stop(sprintf("adj.list[[%d]] contains duplicate neighbors.", i))
        if (any(!is.finite(wts))) stop(sprintf("weight.list[[%d]] contains non-finite values.", i))
        if (any(wts <= 0)) stop(sprintf("weight.list[[%d]] contains non-positive values.", i))
        adj.norm[[i]] <- nbrs
        weight.norm[[i]] <- wts
    }
    for (i in seq_len(n)) {
        for (idx in seq_along(adj.norm[[i]])) {
            j <- adj.norm[[i]][idx]
            rev.idx <- match(i, adj.norm[[j]])
            if (is.na(rev.idx)) {
                stop(sprintf("Graph must be undirected: edge %d -> %d has no reciprocal entry.", i, j))
            }
            w.ij <- weight.norm[[i]][idx]
            w.ji <- weight.norm[[j]][rev.idx]
            if (abs(w.ij - w.ji) > tol * max(1, abs(w.ij), abs(w.ji))) {
                stop(sprintf("Reciprocal edge weights mismatch for (%d, %d): %.12g vs %.12g.",
                             i, j, w.ij, w.ji))
            }
        }
    }
    list(
        adj.list = adj.norm,
        weight.list = weight.norm,
        adj.list.0based = lapply(adj.norm, function(v) as.integer(v - 1L)),
        weight.list.cpp = lapply(weight.norm, as.double)
    )
}

.validate.metric.graph.lowpass.operator.args <- function(conductance.rule,
                                                         conductance.epsilon,
                                                         conductance.alpha,
                                                         conductance.sigma,
                                                         conductance.sigma.rule,
                                                         conductance.sigma.quantile,
                                                         conductance.local.k,
                                                         laplacian.type,
                                                         verbose) {
    conductance.rule <- match.arg(
        conductance.rule,
        choices = c("inverse.length.power", "exp.length",
                    "exp.length.squared", "self.tuned.gaussian")
    )
    conductance.sigma.rule <- match.arg(
        conductance.sigma.rule,
        choices = c("edge.quantile", "median", "local.k")
    )
    laplacian.type <- match.arg(laplacian.type, choices = "unnormalized")
    if (!is.numeric(conductance.epsilon) || length(conductance.epsilon) != 1L ||
        !is.finite(conductance.epsilon) || conductance.epsilon <= 0) {
        stop("conductance.epsilon must be a finite positive numeric scalar.")
    }
    if (!is.numeric(conductance.alpha) || length(conductance.alpha) != 1L ||
        !is.finite(conductance.alpha) || conductance.alpha <= 0) {
        stop("conductance.alpha must be a finite positive numeric scalar.")
    }
    if (is.null(conductance.sigma)) {
        conductance.sigma <- NaN
    } else if (!is.numeric(conductance.sigma) || length(conductance.sigma) != 1L ||
               !is.finite(conductance.sigma) || conductance.sigma <= 0) {
        stop("conductance.sigma must be NULL or a finite positive numeric scalar.")
    }
    if (!is.numeric(conductance.sigma.quantile) || length(conductance.sigma.quantile) != 1L ||
        !is.finite(conductance.sigma.quantile) ||
        conductance.sigma.quantile < 0 || conductance.sigma.quantile > 1) {
        stop("conductance.sigma.quantile must be in [0, 1].")
    }
    conductance.local.k <- .validate.positive.integer.scalar(conductance.local.k, "conductance.local.k")
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.")
    list(
        conductance.rule = conductance.rule,
        conductance.epsilon = as.double(conductance.epsilon),
        conductance.alpha = as.double(conductance.alpha),
        conductance.sigma = as.double(conductance.sigma),
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = as.double(conductance.sigma.quantile),
        conductance.local.k = as.integer(conductance.local.k),
        laplacian.type = laplacian.type,
        verbose = as.logical(verbose)
    )
}

.validate.metric.graph.lowpass.solver.args <- function(n.eigenpairs,
                                                       n,
                                                       eigen.solver,
                                                       dense.eigen.threshold,
                                                       dense.fallback.threshold,
                                                       dense.fallback) {
    n.eigenpairs <- .validate.positive.integer.scalar(n.eigenpairs, "n.eigenpairs")
    if (n.eigenpairs > n) {
        warning("n.eigenpairs exceeds number of vertices; using n.", call. = FALSE)
        n.eigenpairs <- n
    }
    eigen.solver <- match.arg(eigen.solver, choices = c("auto", "sparse", "dense"))
    dense.fallback <- match.arg(dense.fallback, choices = c("auto", "never", "always"))
    dense.eigen.threshold <- .validate.positive.integer.scalar(dense.eigen.threshold, "dense.eigen.threshold")
    dense.fallback.threshold <- .validate.positive.integer.scalar(dense.fallback.threshold, "dense.fallback.threshold")
    list(
        n.eigenpairs = as.integer(n.eigenpairs),
        eigen.solver = eigen.solver,
        dense.eigen.threshold = as.integer(dense.eigen.threshold),
        dense.fallback.threshold = as.integer(dense.fallback.threshold),
        dense.fallback = dense.fallback
    )
}

.validate.positive.integer.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x != floor(x)) {
        stop(sprintf("%s must be a positive integer scalar.", name))
    }
    as.integer(x)
}

.validate.metric.graph.lowpass.response <- function(y, n, name) {
    if (!is.numeric(y) || length(y) != n) {
        stop(sprintf("%s must be a numeric vector of length %d.", name, n))
    }
    y <- as.double(y)
    if (any(!is.finite(y))) stop(sprintf("%s cannot contain NA/NaN/Inf.", name))
    y
}

.prepare.metric.graph.lowpass.response.matrix <- function(y.new, n) {
    is.matrix.input <- is.matrix(y.new) || inherits(y.new, "Matrix")
    if (is.matrix.input) {
        if (nrow(y.new) != n) stop(sprintf("nrow(y.new) must be %d.", n))
        col.names <- colnames(y.new)
        Y <- if (inherits(y.new, "Matrix")) as.matrix(y.new) else y.new
        if (!is.numeric(Y)) stop("y.new must be numeric.")
    } else {
        Y <- matrix(.validate.metric.graph.lowpass.response(y.new, n, "y.new"), ncol = 1)
        col.names <- NULL
    }
    storage.mode(Y) <- "double"
    if (any(!is.finite(Y))) stop("y.new cannot contain NA/NaN/Inf.")
    list(Y = Y, col.names = col.names)
}

.prepare.metric.graph.lowpass.eta.grid <- function(eta.grid, eigenvalues, filter.type, n.candidates) {
    if (is.null(eta.grid)) {
        eta.grid <- generate.eta.grid(eigenvalues, filter.type, n.candidates)
    } else {
        if (!is.numeric(eta.grid) || length(eta.grid) < 1L ||
            any(!is.finite(eta.grid)) || any(eta.grid <= 0)) {
            stop("eta.grid must be NULL or a positive finite numeric vector.")
        }
        eta.grid <- as.double(eta.grid)
    }
    eta.grid
}

.attach.metric.graph.lowpass.laplacian <- function(out, return.sparse) {
    if (isTRUE(return.sparse)) {
        if (requireNamespace("Matrix", quietly = TRUE)) {
            trip <- out$laplacian
            out$laplacian$matrix <- Matrix::sparseMatrix(
                i = trip$i + 1L,
                j = trip$j + 1L,
                x = trip$x,
                dims = trip$dim,
                giveCsparse = TRUE
            )
        } else {
            warning("Matrix is unavailable; returning Laplacian triplets only.", call. = FALSE)
        }
    }
    out
}

.validate.optional.block.size <- function(block.size) {
    if (is.null(block.size)) return(NA_integer_)
    .validate.positive.integer.scalar(block.size, "block.size")
}

.make.metric.graph.lowpass.block.index <- function(p, block.size) {
    if (is.na(block.size) || p <= 1L) return(list(seq_len(p)))
    split(seq_len(p), ceiling(seq_len(p) / block.size))
}

#' @export
print.metric.graph.lowpass.fit <- function(x, ...) {
    cat("\nMetric-Conductance Graph Low-Pass Fit\n")
    cat("====================================\n\n")
    cat(sprintf("Vertices: %d; edges: %d\n", x$graph$n.vertices, x$graph$n.edges))
    cat(sprintf("Conductance rule: %s\n", x$conductance$rule))
    cat(sprintf("Laplacian: %s\n", x$laplacian.type))
    cat(sprintf("Filter: %s; eta: %.4e; GCV: %.4e; effective df: %.2f\n",
                x$spectral$filter.type, x$spectral$eta.optimal,
                x$gcv$gcv.optimal, x$gcv$effective.df))
    invisible(x)
}

#' @export
print.metric.graph.lowpass.refit <- function(x, ...) {
    cat("\nMetric-Conductance Graph Low-Pass Refit\n")
    cat("======================================\n\n")
    cat(sprintf("Responses: %d\n", x$n.responses))
    if (isTRUE(x$per.column.gcv)) {
        cat(sprintf("Per-column GCV: yes; filter: %s\n", x$filter.type))
    } else {
        cat(sprintf("Per-column GCV: no; eta used: %.4e\n", x$eta.used))
    }
    invisible(x)
}

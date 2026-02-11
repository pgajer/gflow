#' PHATE Core Computation (Kernel -> Diffusion Operator -> Potential Distances)
#'
#' @description
#' Minimal native PHATE-style core pipeline for gflow.
#'
#' Given one of \code{X}, \code{D}, \code{K}, or \code{P}, the function:
#' \enumerate{
#'   \item Builds an adaptive alpha-decay kernel \code{K} (if needed),
#'   \item Row-normalizes to a Markov diffusion operator \code{P} (if needed),
#'   \item Chooses diffusion time \code{t} (optional automatic selection),
#'   \item Computes \eqn{P^t},
#'   \item Builds diffusion-potential coordinates \code{U.pot = -log(P^t)},
#'   \item Optionally computes pairwise diffusion-potential distances \code{D.pot}.
#' }
#'
#' This function focuses on robust, reproducible core quantities that can feed
#' downstream graph and conditional-expectation workflows in gflow.
#'
#' @param X Optional numeric data matrix (rows are observations).
#' @param D Optional numeric distance matrix (or \code{dist} object).
#' @param K Optional numeric affinity/kernel matrix.
#' @param P Optional numeric Markov diffusion operator.
#' @param k Integer k-NN neighborhood size for kernel construction from
#'   \code{X} or \code{D}.
#' @param alpha.decay Positive kernel decay exponent.
#' @param knn.use Logical; if \code{TRUE}, build sparse kNN kernel.
#' @param knn.symmetrize Character: \code{"mean"} or \code{"max"} for
#'   symmetrizing directed kNN affinities.
#' @param t Either \code{"auto"} or a positive integer diffusion time.
#' @param t.max Positive integer maximum diffusion time when \code{t="auto"}.
#' @param pca.dim Optional integer. If provided and \code{ncol(X) > pca.dim},
#'   run PCA before distance computation and keep up to \code{pca.dim} PCs.
#' @param pca.center Logical; passed to \code{stats::prcomp}.
#' @param pca.scale Logical; passed to \code{stats::prcomp}.
#' @param kernel.eps Small positive floor for kernel construction.
#' @param potential.eps Small positive floor used in \code{-log(P^t)}.
#' @param compute.D.pot Logical; if \code{TRUE}, compute \code{D.pot}.
#' @param verbose Logical; print progress diagnostics.
#'
#' @return A list of class \code{"phate_core"} with:
#' \describe{
#'   \item{input_type}{One of \code{"X"}, \code{"D"}, \code{"K"}, \code{"P"}.}
#'   \item{K}{Affinity/kernel matrix used for diffusion (if available).}
#'   \item{P}{Row-stochastic diffusion operator.}
#'   \item{Pt}{Diffusion operator raised to \code{t}.}
#'   \item{U.pot}{Diffusion-potential coordinates \eqn{-\log(P^t)}.}
#'   \item{D.pot}{Pairwise Euclidean distances between rows of \code{U.pot} (optional).}
#'   \item{t}{Diffusion time used.}
#'   \item{diagnostics}{List with entropy curve (for auto-\code{t}),
#'     sparsity summaries, and basic dimensions.}
#'   \item{call}{Matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(120), ncol = 4)
#'
#' fit <- phate.core(
#'   X = X,
#'   k = 10,
#'   alpha.decay = 40,
#'   t = "auto",
#'   t.max = 30
#' )
#'
#' fit$t
#' fit$diagnostics$t_auto
#' fit$D.pot[1:5, 1:5]
#' }
#'
#' @export
phate.core <- function(X = NULL,
                       D = NULL,
                       K = NULL,
                       P = NULL,
                       k = 15L,
                       alpha.decay = 40,
                       knn.use = TRUE,
                       knn.symmetrize = c("mean", "max"),
                       t = "auto",
                       t.max = 50L,
                       pca.dim = NULL,
                       pca.center = TRUE,
                       pca.scale = FALSE,
                       kernel.eps = 1e-12,
                       potential.eps = 1e-12,
                       compute.D.pot = TRUE,
                       verbose = FALSE) {

    knn.symmetrize <- match.arg(knn.symmetrize)

    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
        k < 1 || k != floor(k)) {
        stop("k must be a positive integer.")
    }
    k <- as.integer(k)

    if (!is.numeric(alpha.decay) || length(alpha.decay) != 1L ||
        !is.finite(alpha.decay) || alpha.decay <= 0) {
        stop("alpha.decay must be a single positive finite number.")
    }

    if (!is.logical(knn.use) || length(knn.use) != 1L) {
        stop("knn.use must be TRUE/FALSE.")
    }

    if (!is.numeric(t.max) || length(t.max) != 1L || !is.finite(t.max) ||
        t.max < 1 || t.max != floor(t.max)) {
        stop("t.max must be a positive integer.")
    }
    t.max <- as.integer(t.max)

    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1L ||
            !is.finite(pca.dim) || pca.dim < 1 || pca.dim != floor(pca.dim)) {
            stop("pca.dim must be NULL or a positive integer.")
        }
        pca.dim <- as.integer(pca.dim)
    }

    if (!is.numeric(kernel.eps) || length(kernel.eps) != 1L ||
        !is.finite(kernel.eps) || kernel.eps <= 0) {
        stop("kernel.eps must be a single positive finite number.")
    }

    if (!is.numeric(potential.eps) || length(potential.eps) != 1L ||
        !is.finite(potential.eps) || potential.eps <= 0) {
        stop("potential.eps must be a single positive finite number.")
    }

    if (!is.logical(compute.D.pot) || length(compute.D.pot) != 1L) {
        stop("compute.D.pot must be TRUE/FALSE.")
    }

    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("verbose must be TRUE/FALSE.")
    }

    provided <- c(!is.null(X), !is.null(D), !is.null(K), !is.null(P))
    if (sum(provided) != 1L) {
        stop("Provide exactly one of X, D, K, or P.")
    }

    input_type <- c("X", "D", "K", "P")[which(provided)]

    if (isTRUE(verbose)) {
        message(sprintf("PHATE core input: %s", input_type))
    }

    pca.info <- NULL

    if (!is.null(X)) {
        X <- .phc_as_numeric_matrix(X, "X")
        if (!is.null(pca.dim) && ncol(X) > pca.dim) {
            if (isTRUE(verbose)) {
                message(sprintf("Applying PCA to %d dimensions.", pca.dim))
            }
            pca.res <- stats::prcomp(X, center = pca.center, scale. = pca.scale)
            X <- pca.res$x[, seq_len(pca.dim), drop = FALSE]
            pca.info <- list(
                original_dim = ncol(pca.res$rotation),
                projected_dim = pca.dim,
                variance_explained = sum((pca.res$sdev[seq_len(pca.dim)]^2)) /
                    sum(pca.res$sdev^2)
            )
        }
        D <- as.matrix(stats::dist(X))
    }

    if (!is.null(D)) {
        D <- .phc_as_distance_matrix(D, "D")
    }

    n <- if (!is.null(P)) {
        nrow(.phc_as_square_numeric_matrix(P, "P"))
    } else if (!is.null(K)) {
        nrow(.phc_as_square_numeric_matrix(K, "K"))
    } else {
        nrow(D)
    }

    if (n < 2L) {
        stop("PHATE requires at least 2 observations.")
    }

    k.use <- min(k, n - 1L)

    K.used <- NULL
    P.used <- NULL

    if (!is.null(P)) {
        P.used <- .phc_as_square_numeric_matrix(P, "P")
        if (any(P.used < 0)) {
            stop("P must be non-negative.")
        }
        P.used <- .phc_row_normalize(P.used, eps = kernel.eps)
    } else if (!is.null(K)) {
        K.used <- .phc_as_square_numeric_matrix(K, "K")
        if (any(K.used < 0)) {
            stop("K must be non-negative.")
        }
        K.used <- 0.5 * (K.used + t(K.used))
        diag(K.used) <- 0
        P.used <- .phc_row_normalize(K.used, eps = kernel.eps)
    } else {
        K.used <- .phc_build_alpha_kernel(
            D = D,
            k = k.use,
            alpha.decay = alpha.decay,
            knn.use = knn.use,
            knn.symmetrize = knn.symmetrize,
            eps = kernel.eps
        )
        P.used <- .phc_row_normalize(K.used, eps = kernel.eps)
    }

    t.auto <- FALSE
    t.grid <- NULL
    vne <- NULL

    if (is.character(t)) {
        if (length(t) != 1L || t != "auto") {
            stop("t must be a positive integer or the string 'auto'.")
        }
        t.auto <- TRUE
    } else if (is.numeric(t)) {
        if (length(t) != 1L || !is.finite(t) || t < 1 || t != floor(t)) {
            stop("t must be a positive integer or the string 'auto'.")
        }
        t <- as.integer(t)
    } else {
        stop("t must be a positive integer or the string 'auto'.")
    }

    if (isTRUE(t.auto)) {
        if (isTRUE(verbose)) {
            message("Selecting diffusion time t by entropy-knee heuristic.")
        }
        vne.info <- .phc_vne_curve(P.used, t.max = t.max, eps = kernel.eps)
        t <- vne.info$t.opt
        t.grid <- vne.info$t.grid
        vne <- vne.info$entropy
    }

    if (isTRUE(verbose)) {
        message(sprintf("Using diffusion time t = %d", t))
    }

    Pt <- .phc_matrix_power(P.used, t)
    U.pot <- -log(pmax(Pt, potential.eps))

    D.pot <- NULL
    if (isTRUE(compute.D.pot)) {
        if (n > 4000L) {
            warning("Computing D.pot is O(n^2) memory/time and may be expensive for n > 4000.")
        }
        D.pot <- as.matrix(stats::dist(U.pot))
    }

    diagnostics <- list(
        n_vertices = n,
        k_used = k.use,
        alpha_decay = alpha.decay,
        t_auto = t.auto,
        t_grid = t.grid,
        vne = vne,
        kernel_sparsity = if (!is.null(K.used)) mean(K.used == 0) else NA_real_,
        operator_sparsity = mean(P.used == 0),
        pca = pca.info
    )

    res <- list(
        input_type = input_type,
        K = K.used,
        P = P.used,
        Pt = Pt,
        U.pot = U.pot,
        D.pot = D.pot,
        t = as.integer(t),
        diagnostics = diagnostics,
        call = match.call()
    )

    class(res) <- c("phate_core", "list")
    return(res)
}


.phc_as_numeric_matrix <- function(X, name = "matrix") {
    M <- as.matrix(X)
    if (!is.numeric(M) || length(dim(M)) != 2L) {
        stop(sprintf("%s must be a numeric matrix.", name))
    }
    if (any(!is.finite(M))) {
        stop(sprintf("%s contains non-finite values.", name))
    }
    return(M)
}


.phc_as_square_numeric_matrix <- function(X, name = "matrix") {
    M <- .phc_as_numeric_matrix(X, name = name)
    if (nrow(M) != ncol(M)) {
        stop(sprintf("%s must be square.", name))
    }
    return(M)
}


.phc_as_distance_matrix <- function(D, name = "D") {
    M <- if (inherits(D, "dist")) as.matrix(D) else as.matrix(D)
    M <- .phc_as_square_numeric_matrix(M, name = name)
    if (any(M < 0)) {
        stop(sprintf("%s must be non-negative.", name))
    }
    M <- 0.5 * (M + t(M))
    diag(M) <- 0
    return(M)
}


.phc_row_normalize <- function(M, eps = 1e-12) {
    n <- nrow(M)
    rs <- rowSums(M)
    P <- matrix(0, nrow = n, ncol = n)
    ok <- is.finite(rs) & rs > eps

    if (any(ok)) {
        P[ok, ] <- M[ok, , drop = FALSE] / rs[ok]
    }

    if (any(!ok)) {
        idx <- which(!ok)
        P[cbind(idx, idx)] <- 1
    }
    return(P)
}


.phc_matrix_power <- function(M, t) {
    if (t <= 1L) {
        return(M)
    }
    out <- M
    for (iter in 2:t) {
        out <- out %*% M
    }
    return(out)
}


.phc_build_alpha_kernel <- function(D,
                                    k = 15L,
                                    alpha.decay = 40,
                                    knn.use = TRUE,
                                    knn.symmetrize = c("mean", "max"),
                                    eps = 1e-12) {
    knn.symmetrize <- match.arg(knn.symmetrize)
    D <- .phc_as_distance_matrix(D, "D")
    n <- nrow(D)
    k <- min(as.integer(k), n - 1L)

    sigma <- numeric(n)
    for (i in seq_len(n)) {
        row.i <- D[i, ]
        row.i[i] <- Inf
        ord <- order(row.i, decreasing = FALSE)
        kth <- row.i[ord[k]]
        if (!is.finite(kth) || kth <= 0) {
            finite.row <- row.i[is.finite(row.i) & row.i > 0]
            if (length(finite.row) == 0L) {
                kth <- 1
            } else {
                kth <- median(finite.row)
            }
        }
        sigma[i] <- max(kth, eps)
    }

    K.dir <- matrix(0, nrow = n, ncol = n)

    if (isTRUE(knn.use)) {
        for (i in seq_len(n)) {
            row.i <- D[i, ]
            row.i[i] <- Inf
            ord <- order(row.i, decreasing = FALSE)
            nbr <- ord[seq_len(k)]
            vals <- exp(- (row.i[nbr] / sigma[i])^alpha.decay)
            vals[!is.finite(vals)] <- 0
            K.dir[i, nbr] <- vals
        }
    } else {
        for (i in seq_len(n)) {
            vals <- exp(- (D[i, ] / sigma[i])^alpha.decay)
            vals[!is.finite(vals)] <- 0
            vals[i] <- 0
            K.dir[i, ] <- vals
        }
    }

    K <- if (knn.symmetrize == "max") {
        pmax(K.dir, t(K.dir))
    } else {
        0.5 * (K.dir + t(K.dir))
    }

    diag(K) <- 0
    K[K < eps] <- 0
    return(K)
}


.phc_vne_curve <- function(P, t.max = 50L, eps = 1e-12) {
    P <- .phc_as_square_numeric_matrix(P, "P")
    Ps <- 0.5 * (P + t(P))
    eig <- suppressWarnings(
        tryCatch(
            eigen(Ps, symmetric = TRUE, only.values = TRUE)$values,
            error = function(e) NULL
        )
    )

    if (is.null(eig)) {
        eig <- svd(P, nu = 0, nv = 0)$d
    }

    lam <- as.numeric(Re(eig))
    lam[!is.finite(lam)] <- 0
    lam[lam < 0] <- 0

    t.grid <- seq_len(as.integer(t.max))
    entropy <- numeric(length(t.grid))

    for (idx in seq_along(t.grid)) {
        tt <- t.grid[idx]
        wt <- lam^tt
        s <- sum(wt)
        if (!is.finite(s) || s <= eps) {
            entropy[idx] <- 0
        } else {
            p <- wt / s
            p <- p[p > eps]
            entropy[idx] <- -sum(p * log(p))
        }
    }

    t.opt <- .phc_knee_point(entropy)
    if (!is.finite(t.opt) || t.opt < 1L) {
        t.opt <- 1L
    }

    list(
        t.grid = t.grid,
        entropy = entropy,
        t.opt = as.integer(t.opt)
    )
}


.phc_knee_point <- function(y) {
    y <- as.numeric(y)
    n <- length(y)
    if (n <= 2L) {
        return(1L)
    }

    x <- seq_len(n)
    x1 <- x[1L]
    x2 <- x[n]
    y1 <- y[1L]
    y2 <- y[n]

    den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    if (!is.finite(den) || den <= 0) {
        return(which.max(y))
    }

    num <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    idx <- which.max(num)
    if (length(idx) == 0L || !is.finite(idx)) {
        return(1L)
    }
    return(as.integer(idx[1L]))
}


#' Embed Diffusion-Potential Distances into 2D/3D Coordinates
#'
#' @description
#' Computes a low-dimensional embedding from PHATE diffusion-potential geometry
#' using either metric MDS (\code{cmdscale}) or non-metric MDS
#' (\code{MASS::isoMDS}).
#'
#' @param core Optional object returned by \code{\link{phate.core}}.
#' @param D.pot Optional diffusion-potential distance matrix (or \code{dist}).
#' @param U.pot Optional diffusion-potential coordinate matrix. Used when
#'   \code{D.pot} is not provided.
#' @param ndim Embedding dimension (typically 2 or 3).
#' @param method Character: \code{"metric_mds"} or \code{"nonmetric_mds"}.
#' @param maxit Positive integer maximum iterations for non-metric MDS.
#' @param tol Positive convergence tolerance for non-metric MDS.
#' @param cmdscale.add Logical; passed to \code{stats::cmdscale(add = ...)}.
#' @param seed Optional integer seed (used for random fallback initialization).
#' @param verbose Logical; print progress messages.
#'
#' @return A list of class \code{"phate_embedding"} with:
#' \describe{
#'   \item{embedding}{Numeric matrix (\code{n x ndim}) of embedded coordinates.}
#'   \item{method}{Embedding method used.}
#'   \item{stress}{Normalized distance reconstruction stress.}
#'   \item{stress_raw}{Raw non-metric stress from \code{isoMDS} (if applicable).}
#'   \item{cor_spearman}{Spearman correlation between input and embedded distances.}
#'   \item{diagnostics}{List with method-specific diagnostics.}
#'   \item{call}{Matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(120), ncol = 3)
#' core <- phate.core(X = X, t = "auto", t.max = 20)
#' emb <- phate.embed(core = core, ndim = 2, method = "metric_mds")
#' emb$stress
#' }
#'
#' @export
phate.embed <- function(core = NULL,
                        D.pot = NULL,
                        U.pot = NULL,
                        ndim = 2L,
                        method = c("metric_mds", "nonmetric_mds"),
                        maxit = 200L,
                        tol = 1e-3,
                        cmdscale.add = TRUE,
                        seed = NULL,
                        verbose = FALSE) {

    method <- match.arg(method)

    if (!is.numeric(ndim) || length(ndim) != 1L || !is.finite(ndim) ||
        ndim < 1 || ndim != floor(ndim)) {
        stop("ndim must be a positive integer.")
    }
    ndim <- as.integer(ndim)

    if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
        maxit < 1 || maxit != floor(maxit)) {
        stop("maxit must be a positive integer.")
    }
    maxit <- as.integer(maxit)

    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("tol must be a single positive finite number.")
    }

    if (!is.logical(cmdscale.add) || length(cmdscale.add) != 1L) {
        stop("cmdscale.add must be TRUE/FALSE.")
    }

    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("verbose must be TRUE/FALSE.")
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
            stop("seed must be NULL or a single finite numeric value.")
        }
        set.seed(as.integer(seed))
    }

    d.info <- .phe_resolve_distance_input(core = core, D.pot = D.pot, U.pot = U.pot)
    D.in <- d.info$D
    n <- nrow(D.in)

    if (ndim >= n) {
        stop("ndim must be strictly smaller than the number of observations.")
    }

    if (isTRUE(verbose)) {
        message(sprintf("PHATE embedding: method=%s ndim=%d n=%d", method, ndim, n))
    }

    stress.raw <- NA_real_
    method.info <- list()

    if (method == "metric_mds") {
        fit <- stats::cmdscale(stats::as.dist(D.in),
                               k = ndim,
                               eig = TRUE,
                               add = cmdscale.add)
        Y <- as.matrix(fit$points)
        colnames(Y) <- paste0("PHATE", seq_len(ncol(Y)))
        method.info <- list(
            eig = fit$eig,
            ac = if (!is.null(fit$ac)) fit$ac else NA_real_,
            gof = if (!is.null(fit$GOF)) fit$GOF else NA_real_
        )
    } else {
        if (!requireNamespace("MASS", quietly = TRUE)) {
            stop("method='nonmetric_mds' requires package 'MASS'.")
        }

        init <- tryCatch(
            {
                init0 <- stats::cmdscale(stats::as.dist(D.in), k = ndim, add = cmdscale.add)
                if (is.list(init0) && !is.null(init0$points)) {
                    as.matrix(init0$points)
                } else {
                    as.matrix(init0)
                }
            },
            error = function(e) matrix(stats::rnorm(n * ndim), nrow = n, ncol = ndim)
        )

        fit <- MASS::isoMDS(
            d = stats::as.dist(D.in),
            y = init,
            k = ndim,
            maxit = maxit,
            trace = isTRUE(verbose),
            tol = tol
        )

        Y <- as.matrix(fit$points)
        colnames(Y) <- paste0("PHATE", seq_len(ncol(Y)))
        stress.raw <- as.numeric(fit$stress)
        method.info <- list(
            maxit = maxit,
            tol = tol
        )
    }

    D.emb <- as.matrix(stats::dist(Y))
    up <- upper.tri(D.in, diag = FALSE)
    d.ref <- D.in[up]
    d.emb <- D.emb[up]

    denom <- sum(d.ref^2)
    if (!is.finite(denom) || denom <= 0) {
        stress <- sqrt(mean((d.emb - d.ref)^2))
    } else {
        stress <- sqrt(sum((d.emb - d.ref)^2) / denom)
    }

    rho <- suppressWarnings(stats::cor(d.ref, d.emb, method = "spearman"))
    if (!is.finite(rho) || is.na(rho)) {
        rho <- 0
    }

    out <- list(
        embedding = Y,
        method = method,
        stress = stress,
        stress_raw = stress.raw,
        cor_spearman = rho,
        diagnostics = c(
            list(
                n_vertices = n,
                ndim = ndim,
                distance_source = d.info$source
            ),
            method.info
        ),
        call = match.call()
    )

    class(out) <- c("phate_embedding", "list")
    return(out)
}


#' End-to-End PHATE (Core + Embedding)
#'
#' @description
#' One-call wrapper that runs \code{\link{phate.core}} followed by
#' \code{\link{phate.embed}}.
#'
#' @param X,D,K,P Inputs passed to \code{\link{phate.core}} (provide exactly one).
#' @param k,alpha.decay,knn.use,knn.symmetrize,t,t.max,pca.dim,pca.center,pca.scale,
#'   kernel.eps,potential.eps,compute.D.pot Passed to \code{\link{phate.core}}.
#' @param ndim Embedding dimension passed to \code{\link{phate.embed}}.
#' @param embed.method Embedding method passed to \code{\link{phate.embed}}.
#' @param embed.maxit,embed.tol,cmdscale.add,seed Passed to \code{\link{phate.embed}}.
#' @param verbose Logical; print progress.
#'
#' @return A list of class \code{"phate"} with components:
#' \describe{
#'   \item{embedding}{Final low-dimensional coordinates.}
#'   \item{core}{Object from \code{\link{phate.core}}.}
#'   \item{embed}{Object from \code{\link{phate.embed}}.}
#'   \item{t}{Diffusion time used in the core step.}
#'   \item{call}{Matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(150), ncol = 5)
#' fit <- phate(X = X, k = 10, t = "auto", ndim = 2)
#' fit$embedding[1:5, ]
#' fit$t
#' }
#'
#' @export
phate <- function(X = NULL,
                  D = NULL,
                  K = NULL,
                  P = NULL,
                  k = 15L,
                  alpha.decay = 40,
                  knn.use = TRUE,
                  knn.symmetrize = c("mean", "max"),
                  t = "auto",
                  t.max = 50L,
                  pca.dim = NULL,
                  pca.center = TRUE,
                  pca.scale = FALSE,
                  kernel.eps = 1e-12,
                  potential.eps = 1e-12,
                  compute.D.pot = FALSE,
                  ndim = 2L,
                  embed.method = c("metric_mds", "nonmetric_mds"),
                  embed.maxit = 200L,
                  embed.tol = 1e-3,
                  cmdscale.add = TRUE,
                  seed = NULL,
                  verbose = FALSE) {

    embed.method <- match.arg(embed.method)

    core <- phate.core(
        X = X,
        D = D,
        K = K,
        P = P,
        k = k,
        alpha.decay = alpha.decay,
        knn.use = knn.use,
        knn.symmetrize = knn.symmetrize,
        t = t,
        t.max = t.max,
        pca.dim = pca.dim,
        pca.center = pca.center,
        pca.scale = pca.scale,
        kernel.eps = kernel.eps,
        potential.eps = potential.eps,
        compute.D.pot = compute.D.pot,
        verbose = verbose
    )

    emb <- phate.embed(
        core = core,
        ndim = ndim,
        method = embed.method,
        maxit = embed.maxit,
        tol = embed.tol,
        cmdscale.add = cmdscale.add,
        seed = seed,
        verbose = verbose
    )

    out <- list(
        embedding = emb$embedding,
        core = core,
        embed = emb,
        t = core$t,
        call = match.call()
    )
    class(out) <- c("phate", "list")
    return(out)
}


.phe_resolve_distance_input <- function(core = NULL, D.pot = NULL, U.pot = NULL) {
    if (!is.null(D.pot)) {
        D <- .phc_as_distance_matrix(D.pot, "D.pot")
        return(list(D = D, source = "D.pot"))
    }

    if (!is.null(U.pot)) {
        U <- .phc_as_numeric_matrix(U.pot, "U.pot")
        D <- as.matrix(stats::dist(U))
        return(list(D = D, source = "U.pot"))
    }

    if (!is.null(core)) {
        if (!is.list(core)) {
            stop("core must be NULL or a list-like object from phate.core.")
        }
        if (!is.null(core$D.pot)) {
            D <- .phc_as_distance_matrix(core$D.pot, "core$D.pot")
            return(list(D = D, source = "core$D.pot"))
        }
        if (!is.null(core$U.pot)) {
            U <- .phc_as_numeric_matrix(core$U.pot, "core$U.pot")
            D <- as.matrix(stats::dist(U))
            return(list(D = D, source = "core$U.pot"))
        }
    }

    stop("Provide one of: D.pot, U.pot, or core with D.pot/U.pot.")
}

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
#' @param kernel.mode Character: \code{"gflow"} preserves the original gflow
#'   PHATE-style kernel; \code{"phate"} uses the clean-room original
#'   PHATE-compatible kernel path (diagonal-one self affinity, thresholded
#'   alpha-decay affinities, and additive symmetrization). The older
#'   \code{"python"} spelling is accepted as a deprecated alias for
#'   \code{"phate"}.
#' @param kernel.threshold Optional non-negative affinity threshold. If
#'   \code{NULL}, defaults to \code{kernel.eps} in \code{"gflow"} mode and
#'   \code{1e-4} in \code{"phate"} mode.
#' @param bandwidth.scale Positive multiplier for adaptive kernel bandwidths in
#'   \code{"phate"} mode.
#' @param t Either \code{"auto"} or a positive integer diffusion time.
#' @param t.max Positive integer maximum diffusion time when \code{t="auto"}.
#' @param vne.method Character: \code{"auto"} selects \code{"phate"} when
#'   \code{kernel.mode="phate"} and \code{"gflow"} otherwise. \code{"phate"}
#'   uses singular values of the diffusion operator and the PHATE two-line knee
#'   heuristic. \code{"gflow"} preserves the original symmetrized-eigenvalue
#'   heuristic. The older \code{"python"} spelling is accepted as a deprecated
#'   alias for \code{"phate"}.
#' @param pca.dim Optional integer. If provided and \code{ncol(X) > pca.dim},
#'   run PCA before distance computation and keep up to \code{pca.dim} PCs.
#' @param pca.center Logical; passed to \code{stats::prcomp}.
#' @param pca.scale Logical; passed to \code{stats::prcomp}.
#' @param kernel.eps Small positive floor for kernel construction.
#' @param gamma Numeric scalar in \code{[-1, 1]} controlling the PHATE
#'   diffusion-potential transform. \code{gamma=1} is the log potential,
#'   \code{gamma=0} is the square-root potential, and \code{gamma=-1} returns
#'   raw diffused probabilities.
#' @param potential.mode Character: \code{"auto"} selects \code{"phate"} when
#'   \code{kernel.mode="phate"} and \code{"gflow"} otherwise. \code{"phate"}
#'   uses PHATE's \code{-log(P^t + potential.eps)} log transform; \code{"gflow"}
#'   preserves the original floored \code{-log(pmax(P^t, potential.eps))}
#'   behavior when \code{gamma=1}. The older \code{"python"} spelling is
#'   accepted as a deprecated alias for \code{"phate"}.
#' @param potential.eps Small positive floor/offset used in the log potential.
#'   If \code{NULL}, defaults to \code{1e-7} in original PHATE-compatible
#'   potential mode and \code{1e-12} otherwise.
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
	                       kernel.mode = c("gflow", "phate"),
	                       kernel.threshold = NULL,
	                       bandwidth.scale = 1,
	                       t = "auto",
                       t.max = 50L,
                       vne.method = c("auto", "gflow", "phate"),
                       pca.dim = NULL,
                       pca.center = TRUE,
                       pca.scale = FALSE,
                       kernel.eps = 1e-12,
                       gamma = 1,
                       potential.mode = c("auto", "gflow", "phate"),
                       potential.eps = NULL,
                       compute.D.pot = TRUE,
                       verbose = FALSE) {

    knn.symmetrize <- match.arg(knn.symmetrize)
    kernel.mode <- .phc_match_phate_mode(kernel.mode, c("gflow", "phate"), "kernel.mode")
    vne.method <- .phc_match_phate_mode(vne.method, c("auto", "gflow", "phate"), "vne.method")
    potential.mode <- .phc_match_phate_mode(potential.mode, c("auto", "gflow", "phate"),
                                            "potential.mode")

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

    if (!is.null(kernel.threshold)) {
        if (!is.numeric(kernel.threshold) || length(kernel.threshold) != 1L ||
            !is.finite(kernel.threshold) || kernel.threshold < 0 ||
            kernel.threshold >= 1) {
            stop("kernel.threshold must be NULL or a single finite number in [0, 1).")
        }
    }

    if (!is.numeric(bandwidth.scale) || length(bandwidth.scale) != 1L ||
        !is.finite(bandwidth.scale) || bandwidth.scale <= 0) {
        stop("bandwidth.scale must be a single positive finite number.")
    }

    if (!is.numeric(t.max) || length(t.max) != 1L || !is.finite(t.max) ||
        t.max < 1 || t.max != floor(t.max)) {
        stop("t.max must be a positive integer.")
    }
    t.max <- as.integer(t.max)

    vne.method.used <- if (vne.method == "auto") {
        if (kernel.mode == "phate") "phate" else "gflow"
    } else {
        vne.method
    }

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

    kernel.threshold.used <- if (is.null(kernel.threshold)) {
        if (kernel.mode == "phate") 1e-4 else kernel.eps
    } else {
        as.numeric(kernel.threshold)
    }

    if (!is.numeric(gamma) || length(gamma) != 1L ||
        !is.finite(gamma) || gamma < -1 || gamma > 1) {
        stop("gamma must be a single finite number in [-1, 1].")
    }

    potential.mode.used <- if (potential.mode == "auto") {
        if (kernel.mode == "phate") "phate" else "gflow"
    } else {
        potential.mode
    }

    potential.eps.used <- if (is.null(potential.eps)) {
        if (potential.mode.used == "phate") {
            1e-7
        } else {
            1e-12
        }
    } else {
        if (!is.numeric(potential.eps) || length(potential.eps) != 1L ||
            !is.finite(potential.eps) || potential.eps <= 0) {
            stop("potential.eps must be NULL or a single positive finite number.")
        }
        as.numeric(potential.eps)
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
    kernel.info <- NULL

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
        if (kernel.mode == "phate") {
            if (!isTRUE(knn.use)) {
                warning("kernel.mode='phate' uses thresholded adaptive affinities; knn.use=FALSE is ignored.")
            }
            kernel.info <- .phc_build_phate_kernel_cpp(
                D = D,
                k = k.use,
                alpha.decay = alpha.decay,
                knn.symmetrize = knn.symmetrize,
                threshold = kernel.threshold.used,
                bandwidth.scale = bandwidth.scale
            )
            K.used <- kernel.info$K
            P.used <- kernel.info$P
        } else {
            K.used <- .phc_build_alpha_kernel(
                D = D,
                k = k.use,
                alpha.decay = alpha.decay,
                knn.use = knn.use,
                knn.symmetrize = knn.symmetrize,
                eps = kernel.eps,
                threshold = kernel.threshold.used
            )
            P.used <- .phc_row_normalize(K.used, eps = kernel.eps)
        }
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
            message(sprintf("Selecting diffusion time t by %s entropy-knee heuristic.",
                            vne.method.used))
        }
        vne.info <- if (vne.method.used == "phate") {
            .phc_vne_curve_phate(P.used, t.max = t.max)
        } else {
            .phc_vne_curve_gflow(P.used, t.max = t.max, eps = kernel.eps)
        }
        t <- vne.info$t.opt
        t.grid <- vne.info$t.grid
        vne <- vne.info$entropy
    }

    if (isTRUE(verbose)) {
        message(sprintf("Using diffusion time t = %d", t))
    }

    Pt <- .phc_matrix_power(P.used, t)
    U.pot <- .phc_diffusion_potential(
        Pt = Pt,
        gamma = gamma,
        potential.eps = potential.eps.used,
        potential.mode = potential.mode.used
    )

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
        vne_method = vne.method.used,
        kernel_mode = kernel.mode,
        kernel_threshold = kernel.threshold.used,
        bandwidth_scale = bandwidth.scale,
        gamma = gamma,
        potential_mode = potential.mode.used,
        potential_eps = potential.eps.used,
        bandwidth = if (!is.null(kernel.info)) kernel.info$bandwidth else NULL,
        n_retained_directed = if (!is.null(kernel.info)) kernel.info$n_retained_directed else NULL,
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


.phc_match_phate_mode <- function(x, choices, arg) {
    if (missing(x) || is.null(x)) {
        x <- choices[1L]
    }
    if (!is.character(x) || anyNA(x)) {
        stop(sprintf("%s must be one of: %s.", arg, paste(shQuote(choices), collapse = ", ")))
    }
    choices.with.alias <- if ("phate" %in% choices) c(choices, "python") else choices
    if (length(x) > 1L) {
        if (!all(x %in% choices.with.alias)) {
            stop(sprintf("%s must be one of: %s.", arg, paste(shQuote(choices), collapse = ", ")))
        }
        x <- x[1L]
    }
    x <- match.arg(x, choices.with.alias)
    if (identical(x, "python") && "phate" %in% choices) {
        warning(sprintf("%s='python' is deprecated; use %s='phate' for the original PHATE-compatible algorithm.",
                        arg, arg),
                call. = FALSE)
        x <- "phate"
    }
    x
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
    M <- .phc_as_square_numeric_matrix(M, "M")
    if (!is.numeric(t) || length(t) != 1L || !is.finite(t) ||
        t < 0 || t != floor(t)) {
        stop("t must be a non-negative integer for matrix power.")
    }
    t <- as.integer(t)

    if (t == 0L) {
        return(diag(nrow(M)))
    }
    if (t == 1L) {
        return(M)
    }

    out <- diag(nrow(M))
    base <- M
    exp <- t
    while (exp > 0L) {
        if (exp %% 2L == 1L) {
            out <- out %*% base
        }
        exp <- exp %/% 2L
        if (exp > 0L) {
            base <- base %*% base
        }
    }
    return(out)
}


.phc_diffusion_potential <- function(Pt,
                                     gamma = 1,
                                     potential.eps = 1e-12,
                                     potential.mode = c("gflow", "phate")) {
    potential.mode <- .phc_match_phate_mode(potential.mode, c("gflow", "phate"),
                                            "potential.mode")
    Pt <- .phc_as_square_numeric_matrix(Pt, "Pt")
    if (any(Pt < 0)) {
        stop("Pt must be non-negative.")
    }

    if (isTRUE(gamma == 1)) {
        if (potential.mode == "phate") {
            return(-log(Pt + potential.eps))
        }
        return(-log(pmax(Pt, potential.eps)))
    }

    if (isTRUE(gamma == -1)) {
        return(Pt)
    }

    c.gamma <- (1 - gamma) / 2
    Pt.safe <- pmax(Pt, 0)
    return((Pt.safe^c.gamma) / c.gamma)
}


.phc_build_alpha_kernel <- function(D,
                                    k = 15L,
                                    alpha.decay = 40,
                                    knn.use = TRUE,
                                    knn.symmetrize = c("mean", "max"),
                                    eps = 1e-12,
                                    threshold = eps) {
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
    K[K < threshold] <- 0
    return(K)
}


.phc_build_phate_kernel_cpp <- function(D,
                                        k = 15L,
                                        alpha.decay = 40,
                                        knn.symmetrize = c("mean", "max"),
                                        threshold = 1e-4,
                                        bandwidth.scale = 1) {
    knn.symmetrize <- match.arg(knn.symmetrize)
    D <- .phc_as_distance_matrix(D, "D")

    symm.code <- switch(knn.symmetrize,
                        "mean" = 0L,
                        "max" = 1L,
                        stop("Unsupported knn.symmetrize."))

    res <- .Call(
        "S_phate_build_kernel",
        D,
        as.integer(k),
        as.numeric(alpha.decay),
        as.numeric(threshold),
        as.numeric(bandwidth.scale),
        as.integer(symm.code),
        TRUE,
        PACKAGE = "gflow"
    )

    res$K <- .phc_as_square_numeric_matrix(res$K, "K")
    res$P <- .phc_as_square_numeric_matrix(res$P, "P")
    res$bandwidth <- as.numeric(res$bandwidth)
    res$n_retained_directed <- as.integer(res$n_retained_directed)
    return(res)
}


.phc_vne_curve_gflow <- function(P, t.max = 50L, eps = 1e-12) {
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

    t.opt <- .phc_knee_point_gflow(entropy)
    if (!is.finite(t.opt) || t.opt < 1L) {
        t.opt <- 1L
    }

    list(
        t.grid = t.grid,
        entropy = entropy,
        t.opt = as.integer(t.opt)
    )
}


.phc_vne_curve_phate <- function(P, t.max = 50L) {
    P <- .phc_as_square_numeric_matrix(P, "P")
    sv <- svd(P, nu = 0, nv = 0)$d
    sv <- as.numeric(sv)
    sv[!is.finite(sv)] <- 0
    sv[sv < 0] <- 0

    t.grid <- seq.int(0L, as.integer(t.max) - 1L)
    entropy <- numeric(length(t.grid))
    sv.t <- sv

    for (idx in seq_along(t.grid)) {
        s <- sum(sv.t)
        if (!is.finite(s) || s <= 0) {
            entropy[idx] <- 0
        } else {
            prob <- sv.t / s
            prob <- prob + .Machine$double.eps
            entropy[idx] <- -sum(prob * log(prob))
        }
        sv.t <- sv.t * sv
    }

    t.opt <- .phc_knee_point_phate(entropy, x = t.grid)
    if (!is.finite(t.opt) || t.opt < 0L) {
        t.opt <- 0L
    }

    list(
        t.grid = t.grid,
        entropy = entropy,
        t.opt = as.integer(t.opt)
    )
}


.phc_knee_point_gflow <- function(y) {
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


.phc_knee_point_phate <- function(y, x = NULL) {
    y <- as.numeric(y)
    if (length(y) < 3L) {
        stop("Cannot find knee point on vector of length less than 3.")
    }
    if (any(!is.finite(y))) {
        stop("y must contain only finite values.")
    }

    if (is.null(x)) {
        x <- seq_along(y) - 1L
    } else {
        x <- as.numeric(x)
        if (length(x) != length(y)) {
            stop("x and y must have the same length.")
        }
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
    }

    n <- as.numeric(seq.int(2L, length(y)))

    sigma.xy <- cumsum(x * y)[-1L]
    sigma.x <- cumsum(x)[-1L]
    sigma.y <- cumsum(y)[-1L]
    sigma.xx <- cumsum(x * x)[-1L]
    det <- n * sigma.xx - sigma.x * sigma.x
    mfwd <- (n * sigma.xy - sigma.x * sigma.y) / det
    bfwd <- -(sigma.x * sigma.xy - sigma.xx * sigma.y) / det

    xr <- rev(x)
    yr <- rev(y)
    sigma.xy <- cumsum(xr * yr)[-1L]
    sigma.x <- cumsum(xr)[-1L]
    sigma.y <- cumsum(yr)[-1L]
    sigma.xx <- cumsum(xr * xr)[-1L]
    det <- n * sigma.xx - sigma.x * sigma.x
    mbck <- rev((n * sigma.xy - sigma.x * sigma.y) / det)
    bbck <- rev(-(sigma.x * sigma.xy - sigma.xx * sigma.y) / det)

    error.curve <- rep(NA_real_, length(y))
    for (breakpt in seq.int(1L, length(y) - 2L)) {
        left.idx <- seq_len(breakpt + 1L)
        right.idx <- seq.int(breakpt + 1L, length(y))
        dels.fwd <- (mfwd[breakpt] * x[left.idx] + bfwd[breakpt]) - y[left.idx]
        dels.bck <- (mbck[breakpt] * x[right.idx] + bbck[breakpt]) - y[right.idx]
        error.curve[breakpt + 1L] <- sum(abs(dels.fwd)) + sum(abs(dels.bck))
    }

    loc <- which.min(error.curve[seq.int(2L, length(y) - 1L)]) + 1L
    knee.point <- x[loc]
    return(as.integer(knee.point))
}


#' Embed Diffusion-Potential Distances into 2D/3D Coordinates
#'
#' @description
#' Computes a low-dimensional embedding from PHATE diffusion-potential geometry
#' using classical MDS (\code{stats::cmdscale}), metric SMACOF
#' (\code{smacof::smacofSym(type = "ratio")}), or nonmetric SMACOF
#' (\code{smacof::smacofSym(type = "ordinal")}).
#'
#' @param core Optional object returned by \code{\link{phate.core}}.
#' @param D.pot Optional diffusion-potential distance matrix (or \code{dist}).
#' @param U.pot Optional diffusion-potential coordinate matrix. Used when
#'   \code{D.pot} is not provided.
#' @param ndim Embedding dimension (typically 2 or 3).
#' @param method Character: \code{"classic"}, \code{"metric"}, or
#'   \code{"nonmetric"}. Deprecated aliases \code{"metric_mds"} and
#'   \code{"nonmetric_mds"} are accepted for compatibility and map to
#'   \code{"metric"} and \code{"nonmetric"}, respectively.
#' @param maxit Positive integer maximum iterations for SMACOF MDS.
#' @param tol Positive convergence tolerance for SMACOF MDS.
#' @param cmdscale.add Logical; passed to \code{stats::cmdscale(add = ...)}.
#' @param seed Optional integer seed (used for random fallback initialization).
#' @param verbose Logical; print progress messages.
#'
#' @return A list of class \code{"phate_embedding"} with:
#' \describe{
#'   \item{embedding}{Numeric matrix (\code{n x ndim}) of embedded coordinates.}
#'   \item{method}{Canonical embedding method used.}
#'   \item{method_requested}{Embedding method name requested by the caller.}
#'   \item{stress}{Normalized distance reconstruction stress.}
#'   \item{stress_raw}{Raw backend stress from SMACOF (if applicable).}
#'   \item{cor_spearman}{Spearman correlation between input and embedded distances.}
#'   \item{diagnostics}{List with embedding-distance source and method-specific diagnostics.}
#'   \item{call}{Matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(120), ncol = 3)
#' core <- phate.core(X = X, t = "auto", t.max = 20)
#' emb <- phate.embed(core = core, ndim = 2, method = "classic")
#' emb$stress
#' }
#'
#' @export
phate.embed <- function(core = NULL,
                        D.pot = NULL,
                        U.pot = NULL,
                        ndim = 2L,
                        method = c("classic", "metric", "nonmetric",
                                   "metric_mds", "nonmetric_mds"),
                        maxit = 200L,
                        tol = 1e-3,
                        cmdscale.add = TRUE,
                        seed = NULL,
                        verbose = FALSE) {

    method.requested <- if (length(method) > 1L) method[1L] else method
    method <- .phe_match_embedding_method(method)

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

    d.info <- .phe_build_embedding_distance(core = core, D.pot = D.pot, U.pot = U.pot)
    D.in <- d.info$D
    n <- nrow(D.in)

    if (ndim >= n) {
        stop("ndim must be strictly smaller than the number of observations.")
    }

    if (isTRUE(verbose)) {
        message(sprintf("PHATE embedding: method=%s ndim=%d n=%d", method, ndim, n))
    }

    backend.fit <- .phe_run_mds_backend(
        D.in = D.in,
        method = method,
        ndim = ndim,
        maxit = maxit,
        tol = tol,
        cmdscale.add = cmdscale.add,
        verbose = verbose
    )
    Y <- backend.fit$points
    stress.raw <- backend.fit$stress_raw
    method.info <- backend.fit$diagnostics

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
        method_requested = method.requested,
        stress = stress,
        stress_raw = stress.raw,
        cor_spearman = rho,
        diagnostics = c(
            list(
                n_vertices = n,
                ndim = ndim,
                distance_source = d.info$source,
                embedding_distance_source = d.info$embedding_distance_source,
                embedding_distance_input = d.info$embedding_distance_input,
                embedding_distance_backend = d.info$embedding_distance_backend,
                mds_backend = method.info$backend
            ),
            method.info
        ),
        call = match.call()
    )

    class(out) <- c("phate_embedding", "list")
    return(out)
}


.phe_build_embedding_distance <- function(core = NULL, D.pot = NULL, U.pot = NULL) {
    d.info <- .phe_resolve_distance_input(core = core, D.pot = D.pot, U.pot = U.pot)
    d.info$embedding_distance_source <- "potential"
    d.info$embedding_distance_input <- d.info$source
    d.info$embedding_distance_backend <- if (identical(d.info$source, "D.pot") ||
                                            identical(d.info$source, "core$D.pot")) {
        "provided_distance"
    } else {
        "euclidean_rows"
    }
    d.info
}


.phe_run_mds_backend <- function(D.in,
                                 method = c("classic", "metric", "nonmetric"),
                                 ndim = 2L,
                                 maxit = 200L,
                                 tol = 1e-3,
                                 cmdscale.add = TRUE,
                                 verbose = FALSE) {
    method <- match.arg(method)

    if (method == "classic") {
        fit <- .phe_cmdscale_fit(D.in = D.in,
                                 ndim = ndim,
                                 cmdscale.add = cmdscale.add)
        return(list(
            points = fit$points,
            stress_raw = NA_real_,
            diagnostics = c(
                list(backend = "stats::cmdscale"),
                fit$diagnostics
            )
        ))
    }

    if (!requireNamespace("smacof", quietly = TRUE)) {
        stop("method='", method, "' requires package 'smacof'. Install smacof or use method='classic'.")
    }

    init.fit <- .phe_cmdscale_fit(D.in = D.in,
                                  ndim = ndim,
                                  cmdscale.add = cmdscale.add)
    smacof.type <- if (method == "metric") "ratio" else "ordinal"
    fit <- smacof::smacofSym(
        delta = stats::as.dist(D.in),
        ndim = ndim,
        type = smacof.type,
        init = init.fit$points,
        itmax = maxit,
        eps = tol,
        verbose = isTRUE(verbose)
    )

    Y <- as.matrix(fit$conf)
    colnames(Y) <- paste0("PHATE", seq_len(ncol(Y)))
    stress.raw <- as.numeric(fit$stress)

    list(
        points = Y,
        stress_raw = stress.raw,
        diagnostics = list(
            backend = "smacof::smacofSym",
            smacof_type = smacof.type,
            smacof_stress = stress.raw,
            niter = fit$niter,
            init = "classic",
            init_diagnostics = init.fit$diagnostics,
            maxit = maxit,
            tol = tol
        )
    )
}


.phe_match_embedding_method <- function(method) {
    choices <- c("classic", "metric", "nonmetric", "metric_mds", "nonmetric_mds")
    if (!is.character(method) || anyNA(method)) {
        stop("method must be one of: 'classic', 'metric', 'nonmetric'.")
    }
    if (length(method) > 1L) {
        if (!all(method %in% choices)) {
            stop("method must be one of: 'classic', 'metric', 'nonmetric'.")
        }
        method <- method[1L]
    }
    method <- match.arg(method, choices)
    if (method == "metric_mds") {
        warning("method='metric_mds' is deprecated; use method='metric' for metric SMACOF MDS.",
                call. = FALSE)
        return("metric")
    }
    if (method == "nonmetric_mds") {
        warning("method='nonmetric_mds' is deprecated; use method='nonmetric' for nonmetric SMACOF MDS.",
                call. = FALSE)
        return("nonmetric")
    }
    method
}


.phe_cmdscale_fit <- function(D.in, ndim = 2L, cmdscale.add = TRUE) {
    fit <- stats::cmdscale(stats::as.dist(D.in),
                           k = ndim,
                           eig = TRUE,
                           add = cmdscale.add)
    Y <- as.matrix(fit$points)
    if (ncol(Y) < ndim) {
        Y <- cbind(Y, matrix(0, nrow = nrow(Y), ncol = ndim - ncol(Y)))
    }
    if (ncol(Y) > ndim) {
        Y <- Y[, seq_len(ndim), drop = FALSE]
    }
    colnames(Y) <- paste0("PHATE", seq_len(ncol(Y)))
    list(
        points = Y,
        diagnostics = list(
            eig = fit$eig,
            ac = if (!is.null(fit$ac)) fit$ac else NA_real_,
            gof = if (!is.null(fit$GOF)) fit$GOF else NA_real_
        )
    )
}


#' End-to-End PHATE (Core + Embedding)
#'
#' @description
#' One-call wrapper that runs \code{\link{phate.core}} followed by
#' \code{\link{phate.embed}}.
#'
#' @param X Optional numeric data matrix passed to \code{\link{phate.core}}.
#' @param D Optional numeric distance matrix (or \code{dist}) passed to
#'   \code{\link{phate.core}}.
#' @param K Optional affinity/kernel matrix passed to \code{\link{phate.core}}.
#' @param P Optional Markov diffusion operator passed to \code{\link{phate.core}}.
#' @param k Integer k-NN neighborhood size passed to \code{\link{phate.core}}.
#' @param alpha.decay Positive kernel decay exponent passed to
#'   \code{\link{phate.core}}.
#' @param knn.use Logical; if \code{TRUE}, build a sparse kNN kernel in
#'   \code{\link{phate.core}}.
#' @param knn.symmetrize Character in \code{c("mean", "max")} passed to
#'   \code{\link{phate.core}}.
#' @param kernel.mode,kernel.threshold,bandwidth.scale Passed to
#'   \code{\link{phate.core}}.
#' @param t Either \code{"auto"} or a positive integer diffusion time passed to
#'   \code{\link{phate.core}}.
#' @param t.max Positive integer upper bound for automatic \code{t} selection in
#'   \code{\link{phate.core}}.
#' @param vne.method Automatic diffusion-time method passed to
#'   \code{\link{phate.core}}.
#' @param pca.dim Optional integer PCA dimension passed to
#'   \code{\link{phate.core}}.
#' @param pca.center Logical PCA centering flag passed to
#'   \code{\link{phate.core}}.
#' @param pca.scale Logical PCA scaling flag passed to \code{\link{phate.core}}.
#' @param kernel.eps Small positive floor used in kernel normalization by
#'   \code{\link{phate.core}}.
#' @param gamma,potential.mode,potential.eps Diffusion-potential transform
#'   controls passed to \code{\link{phate.core}}.
#' @param compute.D.pot Logical; if \code{TRUE}, compute diffusion-potential
#'   distances in \code{\link{phate.core}}.
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
	                  kernel.mode = c("gflow", "phate"),
	                  kernel.threshold = NULL,
	                  bandwidth.scale = 1,
	                  t = "auto",
                  t.max = 50L,
                  vne.method = c("auto", "gflow", "phate"),
                  pca.dim = NULL,
                  pca.center = TRUE,
                  pca.scale = FALSE,
                  kernel.eps = 1e-12,
                  gamma = 1,
                  potential.mode = c("auto", "gflow", "phate"),
                  potential.eps = NULL,
                  compute.D.pot = FALSE,
                  ndim = 2L,
                  embed.method = c("classic", "metric", "nonmetric",
                                   "metric_mds", "nonmetric_mds"),
                  embed.maxit = 200L,
                  embed.tol = 1e-3,
                  cmdscale.add = TRUE,
                  seed = NULL,
                  verbose = FALSE) {

    embed.method <- .phe_match_embedding_method(embed.method)
    kernel.mode <- .phc_match_phate_mode(kernel.mode, c("gflow", "phate"), "kernel.mode")
    vne.method <- .phc_match_phate_mode(vne.method, c("auto", "gflow", "phate"), "vne.method")
    potential.mode <- .phc_match_phate_mode(potential.mode, c("auto", "gflow", "phate"),
                                            "potential.mode")

    core <- phate.core(
        X = X,
        D = D,
        K = K,
        P = P,
        k = k,
        alpha.decay = alpha.decay,
        knn.use = knn.use,
        knn.symmetrize = knn.symmetrize,
        kernel.mode = kernel.mode,
        kernel.threshold = kernel.threshold,
        bandwidth.scale = bandwidth.scale,
        t = t,
        t.max = t.max,
        vne.method = vne.method,
        pca.dim = pca.dim,
        pca.center = pca.center,
        pca.scale = pca.scale,
        kernel.eps = kernel.eps,
        gamma = gamma,
        potential.mode = potential.mode,
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

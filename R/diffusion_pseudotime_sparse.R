#' Build Sparse Row-Stochastic Transition Matrix from Graph Lists
#'
#' @description
#' Builds a sparse row-stochastic transition operator \eqn{P} from adjacency and
#' edge-weight lists. This is the sparse backend used by diffusion/potential
#' pseudotime routines.
#'
#' @param adj.list Graph adjacency list.
#' @param weight.list Optional edge-weight list aligned with \code{adj.list}. If
#'   \code{NULL}, unit weights are used.
#' @param weight.mode How to convert edge values into transition affinities:
#'   \code{"inverse"} (default), \code{"identity"}, or \code{"exp_neg"}.
#' @param weight.param Positive scalar used by \code{weight.mode}.
#'   For \code{"inverse"} it is the denominator floor; for \code{"exp_neg"} it
#'   is the exponential scale.
#' @param lazy Lazy-walk mixing in \eqn{[0,1]}. \code{1} means full random walk
#'   with no extra self-loop blending.
#'
#' @return A list with:
#' \describe{
#'   \item{P}{Sparse \code{dgCMatrix} transition operator.}
#'   \item{strength}{Numeric vector of transformed outgoing strengths before row normalization.}
#'   \item{indexing}{Detected indexing convention in input adjacency (\code{"zero_based"} or \code{"one_based"}).}
#' }
#'
#' @export
build.sparse.transition <- function(adj.list,
                                    weight.list = NULL,
                                    weight.mode = c("inverse", "identity", "exp_neg"),
                                    weight.param = 1e-6,
                                    lazy = 1) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.")
    }

    weight.mode <- match.arg(weight.mode)
    prep <- .dpt_prepare_graph_inputs(adj.list = adj.list, weight.list = weight.list)
    wm <- .dpt_weight_mode_code(weight.mode)

    res <- .Call(
        "S_build_sparse_transition",
        prep$adj.0based,
        prep$weight.list,
        as.integer(wm),
        as.numeric(weight.param),
        as.numeric(lazy),
        PACKAGE = "gflow"
    )

    list(
        P = .dpt_transition_slots_to_dgCMatrix(res),
        strength = as.numeric(res$strength),
        indexing = prep$indexing
    )
}


#' Compute Rooted Diffusion Pseudotime on Sparse Graphs
#'
#' @description
#' Computes rooted diffusion pseudotime without materializing dense matrix powers.
#' The method propagates a multi-root diffusion score \eqn{q_t = (P^T)^t q_0}
#' with sparse operator applications and defines distance as
#' \eqn{\max(q_t) - q_t}.
#'
#' @param adj.list Graph adjacency list.
#' @param weight.list Optional edge-weight list aligned with \code{adj.list}.
#' @param root.vertices Integer vector of root vertices (1-based or 0-based).
#' @param root.weights Optional non-negative root weights. If \code{NULL},
#'   roots are weighted uniformly.
#' @param t.steps Positive integer diffusion time.
#' @param n.probes Reserved for future stochastic approximations (currently not
#'   used).
#' @param seed Reserved RNG seed argument for forward compatibility (currently
#'   not used).
#' @param weight.mode How to convert edge values into transition affinities:
#'   \code{"inverse"} (default), \code{"identity"}, or \code{"exp_neg"}.
#' @param weight.param Positive scalar used by \code{weight.mode}.
#' @param lazy Lazy-walk mixing in \eqn{[0,1]}.
#' @param normalize Logical; if \code{TRUE}, min-max normalize pseudotime to
#'   \eqn{[0,1]}.
#' @param return.transition Logical; if \code{TRUE}, include sparse transition
#'   matrix in output.
#' @param verbose Logical; print a short summary.
#'
#' @return A list with:
#' \describe{
#'   \item{pseudotime}{Numeric pseudotime vector (length n).}
#'   \item{diffusion.distance}{Rooted diffusion distance proxy (length n).}
#'   \item{root.score}{Root-to-all diffusion score \eqn{q_t}.}
#'   \item{root.vertices}{Root vertices in 1-based indexing.}
#'   \item{method}{\code{"diffusion_sparse"}.}
#'   \item{params}{Parameter list used for the run.}
#'   \item{transition}{Optional sparse \code{dgCMatrix} transition operator.}
#' }
#'
#' @export
compute.diffusion.pseudotime.sparse <- function(adj.list,
                                                weight.list = NULL,
                                                root.vertices,
                                                root.weights = NULL,
                                                t.steps = 8L,
                                                n.probes = 32L,
                                                seed = 1L,
                                                weight.mode = c("inverse", "identity", "exp_neg"),
                                                weight.param = 1e-6,
                                                lazy = 1,
                                                normalize = TRUE,
                                                return.transition = FALSE,
                                                verbose = FALSE) {
    weight.mode <- match.arg(weight.mode)
    prep <- .dpt_prepare_graph_inputs(adj.list = adj.list, weight.list = weight.list)
    n <- length(prep$adj.0based)
    roots.0 <- .dpt_vertices_to_0based(root.vertices = root.vertices, n.vertices = n)
    wm <- .dpt_weight_mode_code(weight.mode)

    res <- .Call(
        "S_compute_diffusion_pseudotime_sparse",
        prep$adj.0based,
        prep$weight.list,
        as.integer(roots.0),
        if (is.null(root.weights)) NULL else as.numeric(root.weights),
        as.integer(t.steps),
        as.integer(n.probes),
        as.integer(seed),
        as.integer(wm),
        as.numeric(weight.param),
        as.numeric(lazy),
        as.logical(normalize),
        as.logical(return.transition),
        PACKAGE = "gflow"
    )

    out <- list(
        pseudotime = as.numeric(res$pseudotime),
        diffusion.distance = as.numeric(res$diffusion.distance),
        root.score = as.numeric(res$root.score),
        root.vertices = as.integer(roots.0 + 1L),
        method = "diffusion_sparse",
        params = list(
            t.steps = as.integer(t.steps),
            n.probes = as.integer(n.probes),
            seed = as.integer(seed),
            weight.mode = weight.mode,
            weight.param = as.numeric(weight.param),
            lazy = as.numeric(lazy),
            normalize = isTRUE(normalize)
        )
    )

    if (isTRUE(return.transition) && !is.null(res$transition)) {
        out$transition <- .dpt_transition_slots_to_dgCMatrix(res$transition)
    }

    if (isTRUE(verbose)) {
        rng <- range(out$pseudotime, na.rm = TRUE)
        cat(sprintf("diffusion pseudotime range: [%.4f, %.4f]\n", rng[1], rng[2]))
    }

    class(out) <- c("diffusion_pseudotime_sparse", "list")
    out
}


#' Compute Rooted Diffusion-Potential Pseudotime on Sparse Graphs
#'
#' @description
#' Computes rooted potential pseudotime using sparse diffusion:
#' \eqn{\tau_i = -\log((P^t r)_i + \epsilon)}.
#' Optionally computes landmark potential coordinates for approximate pairwise
#' potential geometry.
#'
#' @param adj.list Graph adjacency list.
#' @param weight.list Optional edge-weight list aligned with \code{adj.list}.
#' @param root.vertices Integer vector of root vertices (1-based or 0-based).
#' @param root.weights Optional non-negative root weights. If \code{NULL},
#'   roots are weighted uniformly.
#' @param t.steps Positive integer diffusion time.
#' @param potential.eps Positive floor used in \code{-log(pmax(score, eps))}.
#' @param landmark.vertices Optional landmark vertices (1-based or 0-based) for
#'   approximate potential geometry.
#' @param return.landmark.distances Logical; if \code{TRUE} and landmarks are
#'   provided, compute pairwise Euclidean distances in landmark-potential space.
#' @param max.landmark.distance.n Maximum \code{n} for dense landmark-distance
#'   materialization.
#' @param weight.mode How to convert edge values into transition affinities:
#'   \code{"inverse"} (default), \code{"identity"}, or \code{"exp_neg"}.
#' @param weight.param Positive scalar used by \code{weight.mode}.
#' @param lazy Lazy-walk mixing in \eqn{[0,1]}.
#' @param normalize Logical; if \code{TRUE}, min-max normalize pseudotime to
#'   \eqn{[0,1]}.
#' @param return.transition Logical; if \code{TRUE}, include sparse transition
#'   matrix in output.
#' @param verbose Logical; print a short summary.
#'
#' @return A list with:
#' \describe{
#'   \item{pseudotime}{Numeric pseudotime vector (length n).}
#'   \item{root.score}{Root-to-all diffusion score \eqn{P^t r}.}
#'   \item{root.potential}{Unnormalized \eqn{-\log(P^t r + \epsilon)}.}
#'   \item{landmark.potential}{Optional \eqn{n \times L} landmark potential matrix.}
#'   \item{landmark.distances}{Optional dense pairwise distances in landmark space.}
#'   \item{root.vertices}{Root vertices in 1-based indexing.}
#'   \item{method}{\code{"potential_sparse"}.}
#'   \item{params}{Parameter list used for the run.}
#'   \item{transition}{Optional sparse \code{dgCMatrix} transition operator.}
#' }
#'
#' @export
compute.potential.pseudotime.sparse <- function(adj.list,
                                                weight.list = NULL,
                                                root.vertices,
                                                root.weights = NULL,
                                                t.steps = 8L,
                                                potential.eps = 1e-12,
                                                landmark.vertices = NULL,
                                                return.landmark.distances = FALSE,
                                                max.landmark.distance.n = 5000L,
                                                weight.mode = c("inverse", "identity", "exp_neg"),
                                                weight.param = 1e-6,
                                                lazy = 1,
                                                normalize = TRUE,
                                                return.transition = FALSE,
                                                verbose = FALSE) {
    weight.mode <- match.arg(weight.mode)
    prep <- .dpt_prepare_graph_inputs(adj.list = adj.list, weight.list = weight.list)
    n <- length(prep$adj.0based)
    roots.0 <- .dpt_vertices_to_0based(root.vertices = root.vertices, n.vertices = n)
    landmarks.0 <- if (is.null(landmark.vertices)) integer(0) else {
        .dpt_vertices_to_0based(root.vertices = landmark.vertices, n.vertices = n)
    }
    wm <- .dpt_weight_mode_code(weight.mode)

    res <- .Call(
        "S_compute_potential_pseudotime_sparse",
        prep$adj.0based,
        prep$weight.list,
        as.integer(roots.0),
        if (is.null(root.weights)) NULL else as.numeric(root.weights),
        as.integer(t.steps),
        as.numeric(potential.eps),
        as.integer(landmarks.0),
        as.integer(wm),
        as.numeric(weight.param),
        as.numeric(lazy),
        as.logical(normalize),
        as.logical(return.transition),
        PACKAGE = "gflow"
    )

    out <- list(
        pseudotime = as.numeric(res$pseudotime),
        root.score = as.numeric(res$root.score),
        root.potential = as.numeric(res$root.potential),
        root.vertices = as.integer(roots.0 + 1L),
        method = "potential_sparse",
        params = list(
            t.steps = as.integer(t.steps),
            potential.eps = as.numeric(potential.eps),
            landmarks = as.integer(landmarks.0 + 1L),
            weight.mode = weight.mode,
            weight.param = as.numeric(weight.param),
            lazy = as.numeric(lazy),
            normalize = isTRUE(normalize)
        )
    )

    if (!is.null(res$landmark.potential)) {
        out$landmark.potential <- as.matrix(res$landmark.potential)
    }

    if (isTRUE(return.transition) && !is.null(res$transition)) {
        out$transition <- .dpt_transition_slots_to_dgCMatrix(res$transition)
    }

    if (isTRUE(return.landmark.distances) && !is.null(out$landmark.potential)) {
        if (n <= as.integer(max.landmark.distance.n)) {
            out$landmark.distances <- as.matrix(stats::dist(out$landmark.potential))
        } else {
            warning(sprintf(
                paste0("Skipping landmark distance matrix (n=%d > max.landmark.distance.n=%d). ",
                       "Increase max.landmark.distance.n to force computation."),
                n, as.integer(max.landmark.distance.n)
            ))
        }
    }

    if (isTRUE(verbose)) {
        rng <- range(out$pseudotime, na.rm = TRUE)
        cat(sprintf("potential pseudotime range: [%.4f, %.4f]\n", rng[1], rng[2]))
    }

    class(out) <- c("potential_pseudotime_sparse", "list")
    out
}


.dpt_weight_mode_code <- function(weight.mode) {
    switch(
        weight.mode,
        identity = 0L,
        inverse = 1L,
        exp_neg = 2L,
        stop("Unknown weight.mode: ", weight.mode)
    )
}


.dpt_transition_slots_to_dgCMatrix <- function(x) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.")
    }
    Matrix::sparseMatrix(
        i = as.integer(x$i),
        p = as.integer(x$p),
        x = as.numeric(x$x),
        dims = as.integer(x$Dim),
        index1 = FALSE
    )
}


.dpt_vertices_to_0based <- function(root.vertices, n.vertices) {
    if (!is.numeric(root.vertices)) {
        stop("root.vertices must be numeric/integer.")
    }
    rv <- as.integer(root.vertices)
    if (length(rv) < 1L || anyNA(rv)) {
        stop("root.vertices must contain at least one non-NA index.")
    }

    if (all(rv >= 1L & rv <= n.vertices)) {
        return(rv - 1L)
    }
    if (all(rv >= 0L & rv <= (n.vertices - 1L))) {
        return(rv)
    }
    stop("root.vertices must be either 1-based (1..n) or 0-based (0..n-1).")
}


.dpt_prepare_graph_inputs <- function(adj.list, weight.list = NULL) {
    if (!is.list(adj.list)) {
        stop("adj.list must be a list.")
    }
    n <- length(adj.list)
    if (n < 2L) {
        stop("adj.list must have length >= 2.")
    }

    if (is.null(weight.list)) {
        weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
    }
    if (!is.list(weight.list) || length(weight.list) != n) {
        stop("weight.list must be NULL or a list with same length as adj.list.")
    }

    adj.int <- lapply(adj.list, function(x) as.integer(x))
    w.num <- vector("list", n)
    for (i in seq_len(n)) {
        w.num[[i]] <- as.numeric(weight.list[[i]])
        if (length(w.num[[i]]) != length(adj.int[[i]])) {
            stop(sprintf("weight.list[[%d]] length must match adj.list[[%d]] length.", i, i))
        }
    }

    all.nbr <- unlist(adj.int, use.names = FALSE)
    if (length(all.nbr) < 1L) {
        stop("adj.list has no edges.")
    }

    if (all(all.nbr >= 1L & all.nbr <= n)) {
        indexing <- "one_based"
        adj.0 <- lapply(adj.int, function(x) x - 1L)
    } else if (all(all.nbr >= 0L & all.nbr <= (n - 1L))) {
        indexing <- "zero_based"
        adj.0 <- adj.int
    } else {
        stop("adj.list indices must be either fully 1-based (1..n) or 0-based (0..n-1).")
    }

    list(
        adj.0based = adj.0,
        weight.list = w.num,
        indexing = indexing
    )
}

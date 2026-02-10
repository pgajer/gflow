#' Construct a Sparse Graph Metric Aligned to Diffusion-Potential Distances
#'
#' @description
#' Fits positive edge lengths on a fixed sparse undirected topology so that
#' graph shortest-path distances approximate target diffusion-potential distances.
#' The optimization focuses on mid/long-range vertex pairs while anchoring local
#' edge lengths to direct targets on existing edges.
#'
#' @details
#' The function keeps the input edge topology fixed and iterates:
#' \enumerate{
#'   \item Build shortest paths for selected vertex pairs under current edge lengths.
#'   \item Linearize pair distances via path-incidence matrix \eqn{A}.
#'   \item Solve a regularized weighted least-squares problem for updated edge lengths.
#'   \item Recompute shortest paths and repeat until convergence.
#' }
#'
#' For fixed paths, the subproblem is:
#' \deqn{
#' \min_{w \ge \epsilon}
#' \sum_{(i,j)\in\mathcal{P}} \omega_{ij} \left((Aw)_{ij} - D_{ij}\right)^2
#' + \lambda_{\text{local}}\|w - d_E\|_2^2
#' + \lambda_{\text{ridge}}\|w\|_2^2
#' }
#' where \eqn{D_{ij}} are target diffusion-potential distances and \eqn{d_E} are
#' edge-local targets from \code{D.pot} on existing edges.
#'
#' @param adj.list List of integer neighbor vectors (1-based indices), defining
#'   the fixed sparse topology.
#' @param weight.list Optional list of numeric edge lengths matching
#'   \code{adj.list}. If \code{NULL}, initialization uses edge targets from
#'   \code{D.pot}.
#' @param D.pot Optional numeric matrix (or \code{dist}) of target
#'   diffusion-potential distances.
#' @param U.pot Optional numeric matrix of potential coordinates. Used only when
#'   \code{D.pot} is \code{NULL}; then \code{stats::dist(U.pot)} defines targets.
#' @param pair.quantiles Length-2 numeric vector in \eqn{[0,1]} giving lower and
#'   upper quantiles used to select target pairs (default mid/long range).
#' @param pair.mode Character: \code{"landmark"} or \code{"all"}.
#' @param n.landmarks Integer number of landmarks used when
#'   \code{pair.mode = "landmark"}.
#' @param max.pairs Maximum number of selected pairs used in fitting.
#' @param long.range.power Non-negative exponent for pair weights:
#'   \eqn{\omega_{ij} \propto D_{ij}^{\beta}}.
#' @param eta.blend Blend weight in \eqn{[0,1]} for initialization:
#'   \eqn{w_0 = \eta d_E + (1-\eta) w_{\text{orig}}}.
#' @param lambda.local Non-negative regularization weight for local edge anchor.
#' @param lambda.ridge Non-negative ridge regularization on edge lengths.
#' @param max.iter Maximum number of outer iterations.
#' @param tol Relative objective-improvement tolerance for convergence.
#' @param min.edge.length Strictly positive lower bound for edge lengths.
#' @param seed Optional integer seed used only when random subsampling is needed.
#' @param verbose Logical; print iteration diagnostics.
#'
#' @return A list with components:
#' \describe{
#'   \item{adj_list}{Adjacency list of fitted graph (same topology).}
#'   \item{weight_list}{Fitted edge lengths aligned with \code{adj_list}.}
#'   \item{edge_matrix}{Two-column undirected edge matrix (1-based).}
#'   \item{edge_weights}{Numeric vector of fitted edge lengths per row in \code{edge_matrix}.}
#'   \item{pairs}{Two-column matrix of fitted pair indices.}
#'   \item{pair_targets}{Target distances for \code{pairs}.}
#'   \item{pair_predictions}{Final fitted shortest-path distances for \code{pairs}.}
#'   \item{pair_weights}{Pair weights used in objective.}
#'   \item{diagnostics}{List with objective trace, stress, correlation, convergence info.}
#' }
#'
#' @examples
#' \dontrun{
#' # Chain graph on 8 vertices
#' n <- 8
#' adj <- vector("list", n)
#' w <- vector("list", n)
#' for (i in seq_len(n)) {
#'   nbr <- integer(0)
#'   if (i > 1) nbr <- c(nbr, i - 1L)
#'   if (i < n) nbr <- c(nbr, i + 1L)
#'   adj[[i]] <- nbr
#'   w[[i]] <- rep(1, length(nbr))
#' }
#'
#' # Potential coordinates and corresponding target distances
#' U <- cbind(seq_len(n))
#'
#' fit <- construct.potential.metric.graph(
#'   adj.list = adj,
#'   weight.list = w,
#'   U.pot = U,
#'   pair.mode = "all",
#'   max.iter = 5,
#'   verbose = TRUE
#' )
#'
#' fit$diagnostics$stress_weighted
#' }
#'
#' @export
construct.potential.metric.graph <- function(adj.list,
                                             weight.list = NULL,
                                             D.pot = NULL,
                                             U.pot = NULL,
                                             pair.quantiles = c(0.4, 0.95),
                                             pair.mode = c("landmark", "all"),
                                             n.landmarks = 128L,
                                             max.pairs = 5000L,
                                             long.range.power = 1.5,
                                             eta.blend = 0.7,
                                             lambda.local = 1.0,
                                             lambda.ridge = 1e-3,
                                             max.iter = 15L,
                                             tol = 1e-3,
                                             min.edge.length = 1e-8,
                                             seed = NULL,
                                             verbose = FALSE) {

    .pmg_validate_graph_inputs(adj.list, weight.list)
    n.vertices <- length(adj.list)

    pair.mode <- match.arg(pair.mode)

    if (!is.numeric(pair.quantiles) || length(pair.quantiles) != 2L ||
        any(!is.finite(pair.quantiles)) || any(pair.quantiles < 0) ||
        any(pair.quantiles > 1) || pair.quantiles[1] >= pair.quantiles[2]) {
        stop("pair.quantiles must be length-2 numeric in [0,1] with q1 < q2.")
    }

    if (!is.numeric(long.range.power) || length(long.range.power) != 1L ||
        !is.finite(long.range.power) || long.range.power < 0) {
        stop("long.range.power must be a single non-negative finite number.")
    }

    if (!is.numeric(eta.blend) || length(eta.blend) != 1L ||
        !is.finite(eta.blend) || eta.blend < 0 || eta.blend > 1) {
        stop("eta.blend must be a single numeric in [0,1].")
    }

    if (!is.numeric(lambda.local) || length(lambda.local) != 1L ||
        !is.finite(lambda.local) || lambda.local < 0) {
        stop("lambda.local must be a single non-negative finite number.")
    }

    if (!is.numeric(lambda.ridge) || length(lambda.ridge) != 1L ||
        !is.finite(lambda.ridge) || lambda.ridge < 0) {
        stop("lambda.ridge must be a single non-negative finite number.")
    }

    if (!is.numeric(max.iter) || length(max.iter) != 1L ||
        !is.finite(max.iter) || max.iter < 1 || max.iter != floor(max.iter)) {
        stop("max.iter must be a positive integer.")
    }
    max.iter <- as.integer(max.iter)

    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("tol must be a single positive finite number.")
    }

    if (!is.numeric(min.edge.length) || length(min.edge.length) != 1L ||
        !is.finite(min.edge.length) || min.edge.length <= 0) {
        stop("min.edge.length must be a single positive finite number.")
    }

    if (!is.numeric(n.landmarks) || length(n.landmarks) != 1L ||
        !is.finite(n.landmarks) || n.landmarks < 2 ||
        n.landmarks != floor(n.landmarks)) {
        stop("n.landmarks must be an integer >= 2.")
    }
    n.landmarks <- as.integer(n.landmarks)

    if (!is.numeric(max.pairs) || length(max.pairs) != 1L ||
        !is.finite(max.pairs) || max.pairs < 1 ||
        max.pairs != floor(max.pairs)) {
        stop("max.pairs must be a positive integer.")
    }
    max.pairs <- as.integer(max.pairs)

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
            stop("seed must be NULL or a single finite numeric value.")
        }
        set.seed(as.integer(seed))
    }

    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("verbose must be TRUE/FALSE.")
    }

    D <- .pmg_resolve_target_distance_matrix(D.pot, U.pot, n.vertices)

    edge.info <- .pmg_extract_undirected_edges(adj.list, weight.list)
    edge.mat <- edge.info$edge_matrix
    w.orig <- edge.info$edge_weights
    n.edges <- nrow(edge.mat)

    if (n.edges < 1L) {
        stop("Input graph has no undirected edges.")
    }

    g0 <- .pmg_graph_from_edges(edge.mat, rep(1, n.edges), n.vertices)
    if (igraph::components(g0)$no > 1L) {
        stop("Input topology must be connected for shortest-path fitting.")
    }

    d.edge.target <- D[cbind(edge.mat[, 1L], edge.mat[, 2L])]
    d.edge.target[!is.finite(d.edge.target) | d.edge.target <= 0] <- min.edge.length

    if (is.null(w.orig)) {
        w.orig <- d.edge.target
    } else {
        w.orig <- as.numeric(w.orig)
        w.orig[!is.finite(w.orig) | w.orig <= 0] <- d.edge.target[!is.finite(w.orig) | w.orig <= 0]
    }

    w.current <- eta.blend * d.edge.target + (1 - eta.blend) * w.orig
    w.current <- pmax(min.edge.length, w.current)

    pair.sel <- .pmg_select_pairs(D = D,
                                  pair.quantiles = pair.quantiles,
                                  pair.mode = pair.mode,
                                  n.landmarks = n.landmarks,
                                  max.pairs = max.pairs)

    pairs <- pair.sel$pairs
    n.pairs.candidate <- pair.sel$n_pairs_candidate
    d.target <- D[cbind(pairs[, 1L], pairs[, 2L])]
    pair.weights <- d.target^long.range.power
    if (!all(is.finite(pair.weights)) || any(pair.weights <= 0)) {
        pair.weights <- rep(1.0, length(d.target))
    }

    eval.current <- .pmg_evaluate_weights(edge.mat = edge.mat,
                                          edge.weights = w.current,
                                          n.vertices = n.vertices,
                                          pairs = pairs,
                                          d.target = d.target,
                                          pair.weights = pair.weights,
                                          d.edge.target = d.edge.target,
                                          lambda.local = lambda.local,
                                          lambda.ridge = lambda.ridge)

    objective.trace <- c(eval.current$objective)
    stress.trace <- c(eval.current$stress_weighted)
    corr.trace <- c(eval.current$corr_spearman)
    rel.improvement.trace <- numeric(0)

    converged <- FALSE
    n.iters <- 0L

    if (isTRUE(verbose)) {
        message(sprintf("PMGraph init: obj=%.6g stress=%.6g rho=%.4f n_pairs=%d n_edges=%d",
                        eval.current$objective,
                        eval.current$stress_weighted,
                        eval.current$corr_spearman,
                        nrow(pairs),
                        n.edges))
    }

    for (iter in seq_len(max.iter)) {
        n.iters <- iter

        w.proposed <- .pmg_solve_fixed_paths(
            A = eval.current$incidence,
            d.target = d.target,
            pair.weights = pair.weights,
            d.edge.target = d.edge.target,
            lambda.local = lambda.local,
            lambda.ridge = lambda.ridge,
            min.edge.length = min.edge.length
        )

        eval.proposed <- .pmg_evaluate_weights(edge.mat = edge.mat,
                                               edge.weights = w.proposed,
                                               n.vertices = n.vertices,
                                               pairs = pairs,
                                               d.target = d.target,
                                               pair.weights = pair.weights,
                                               d.edge.target = d.edge.target,
                                               lambda.local = lambda.local,
                                               lambda.ridge = lambda.ridge)

        best.eval <- eval.current
        best.w <- w.current

        if (eval.proposed$objective < best.eval$objective) {
            best.eval <- eval.proposed
            best.w <- w.proposed
        }

        if (best.eval$objective > eval.current$objective) {
            alpha <- 0.5
            while (alpha >= (1 / 64)) {
                w.try <- pmax(min.edge.length, alpha * w.proposed + (1 - alpha) * w.current)
                eval.try <- .pmg_evaluate_weights(edge.mat = edge.mat,
                                                  edge.weights = w.try,
                                                  n.vertices = n.vertices,
                                                  pairs = pairs,
                                                  d.target = d.target,
                                                  pair.weights = pair.weights,
                                                  d.edge.target = d.edge.target,
                                                  lambda.local = lambda.local,
                                                  lambda.ridge = lambda.ridge)
                if (eval.try$objective < best.eval$objective) {
                    best.eval <- eval.try
                    best.w <- w.try
                }
                alpha <- alpha * 0.5
            }
        }

        rel.improvement <- (eval.current$objective - best.eval$objective) /
            max(1, abs(eval.current$objective))
        rel.improvement.trace <- c(rel.improvement.trace, rel.improvement)

        w.current <- best.w
        eval.current <- best.eval

        objective.trace <- c(objective.trace, eval.current$objective)
        stress.trace <- c(stress.trace, eval.current$stress_weighted)
        corr.trace <- c(corr.trace, eval.current$corr_spearman)

        if (isTRUE(verbose)) {
            message(sprintf("PMGraph iter=%d obj=%.6g rel_impr=%.3e stress=%.6g rho=%.4f",
                            iter,
                            eval.current$objective,
                            rel.improvement,
                            eval.current$stress_weighted,
                            eval.current$corr_spearman))
        }

        if (rel.improvement >= 0 && rel.improvement < tol) {
            converged <- TRUE
            break
        }
    }

    fitted.graph <- .pmg_edges_to_lists(edge.mat,
                                        w.current,
                                        n.vertices,
                                        names(adj.list))

    diagnostics <- list(
        objective_trace = objective.trace,
        relative_improvement_trace = rel.improvement.trace,
        stress_trace = stress.trace,
        corr_spearman_trace = corr.trace,
        stress_weighted = eval.current$stress_weighted,
        corr_spearman_midlong = eval.current$corr_spearman,
        local_edge_rmse = sqrt(mean((w.current - d.edge.target)^2)),
        pair_coverage = if (n.pairs.candidate > 0) nrow(pairs) / n.pairs.candidate else 1.0,
        n_pairs = nrow(pairs),
        n_pairs_candidate = n.pairs.candidate,
        n_edges = n.edges,
        n_iters = n.iters,
        converged = converged,
        pair_quantile_range = pair.sel$distance_quantile_range
    )

    result <- list(
        adj_list = fitted.graph$adj_list,
        weight_list = fitted.graph$weight_list,
        edge_matrix = edge.mat,
        edge_weights = w.current,
        pairs = pairs,
        pair_targets = d.target,
        pair_predictions = eval.current$d.pred,
        pair_weights = pair.weights,
        diagnostics = diagnostics
    )

    class(result) <- c("potential_metric_graph", "list")
    return(result)
}


#' Construct Potential-Metric Graph from PHATE-Style Inputs
#'
#' @description
#' Convenience wrapper for \code{\link{construct.potential.metric.graph}} that
#' accepts PHATE-style objects directly: affinity/kernel matrix \code{K},
#' diffusion operator \code{P}, potential coordinates \code{U.pot}, or target
#' potential-distance matrix \code{D.pot}.
#'
#' The wrapper can automatically construct a sparse undirected topology from:
#' \itemize{
#'   \item nonzero operator entries (\code{topology.method = "operator_nonzero"}),
#'   \item operator kNN neighborhoods (\code{topology.method = "operator_knn"}),
#'   \item kNN in potential space (\code{topology.method = "potential_knn"}).
#' }
#'
#' @param K Optional numeric square affinity/kernel matrix.
#' @param P Optional numeric square diffusion/transition operator matrix.
#' @param U.pot Optional numeric matrix of potential coordinates.
#' @param D.pot Optional numeric matrix (or \code{dist}) of target
#'   diffusion-potential distances.
#' @param adj.list Optional fixed adjacency list. If provided, automatic topology
#'   construction is skipped.
#' @param weight.list Optional edge-length list parallel to \code{adj.list}.
#' @param topology.method Character: \code{"auto"}, \code{"operator_nonzero"},
#'   \code{"operator_knn"}, or \code{"potential_knn"}.
#' @param k.topology Integer k used when a kNN topology mode is selected.
#' @param knn.symmetrize Character: \code{"union"} or \code{"mutual"} for
#'   converting directed kNN neighborhoods into an undirected graph.
#' @param min.operator.weight Non-negative threshold for retaining operator edges.
#' @param diffusion.time Positive integer \(t\) used to form
#'   diffusion potential from \code{P^t} when only \code{K}/\code{P} are provided.
#' @param operator.as.markov Logical. If \code{TRUE}, row-normalize operator
#'   before diffusion-potential construction.
#' @param potential.eps Small positive floor used in \code{-log(pmax(P^t, eps))}.
#' @param edge.length.transform Character: \code{"neglog"} or
#'   \code{"reciprocal"} mapping operator similarities to edge lengths.
#' @param edge.length.eps Positive floor used in similarity-to-length transform.
#' @param pair.quantiles,pair.mode,n.landmarks,max.pairs,long.range.power,eta.blend,
#'   lambda.local,lambda.ridge,max.iter,tol,min.edge.length,seed,verbose Passed to
#'   \code{\link{construct.potential.metric.graph}}.
#'
#' @return The same object as \code{construct.potential.metric.graph}, with an
#' additional \code{$wrapper_info} list describing how topology/targets were derived.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(40), ncol = 2)
#' D <- as.matrix(stats::dist(X))
#' K <- exp(-D^2)
#'
#' fit <- construct.potential.metric.graph.from.phate(
#'   K = K,
#'   topology.method = "operator_knn",
#'   k.topology = 5,
#'   diffusion.time = 1,
#'   pair.mode = "landmark"
#' )
#'
#' fit$wrapper_info
#' fit$diagnostics$stress_weighted
#' }
#'
#' @export
construct.potential.metric.graph.from.phate <- function(K = NULL,
                                                        P = NULL,
                                                        U.pot = NULL,
                                                        D.pot = NULL,
                                                        adj.list = NULL,
                                                        weight.list = NULL,
                                                        topology.method = c("auto",
                                                                            "operator_nonzero",
                                                                            "operator_knn",
                                                                            "potential_knn"),
                                                        k.topology = 15L,
                                                        knn.symmetrize = c("union", "mutual"),
                                                        min.operator.weight = 0,
                                                        diffusion.time = 1L,
                                                        operator.as.markov = TRUE,
                                                        potential.eps = 1e-12,
                                                        edge.length.transform = c("neglog", "reciprocal"),
                                                        edge.length.eps = 1e-12,
                                                        pair.quantiles = c(0.4, 0.95),
                                                        pair.mode = c("landmark", "all"),
                                                        n.landmarks = 128L,
                                                        max.pairs = 5000L,
                                                        long.range.power = 1.5,
                                                        eta.blend = 0.7,
                                                        lambda.local = 1.0,
                                                        lambda.ridge = 1e-3,
                                                        max.iter = 15L,
                                                        tol = 1e-3,
                                                        min.edge.length = 1e-8,
                                                        seed = NULL,
                                                        verbose = FALSE) {

    topology.method <- match.arg(topology.method)
    knn.symmetrize <- match.arg(knn.symmetrize)
    edge.length.transform <- match.arg(edge.length.transform)

    if (!is.numeric(k.topology) || length(k.topology) != 1L ||
        !is.finite(k.topology) || k.topology < 1 ||
        k.topology != floor(k.topology)) {
        stop("k.topology must be a positive integer.")
    }
    k.topology <- as.integer(k.topology)

    if (!is.numeric(min.operator.weight) || length(min.operator.weight) != 1L ||
        !is.finite(min.operator.weight) || min.operator.weight < 0) {
        stop("min.operator.weight must be a single non-negative finite number.")
    }

    if (!is.numeric(diffusion.time) || length(diffusion.time) != 1L ||
        !is.finite(diffusion.time) || diffusion.time < 1 ||
        diffusion.time != floor(diffusion.time)) {
        stop("diffusion.time must be a positive integer.")
    }
    diffusion.time <- as.integer(diffusion.time)

    if (!is.logical(operator.as.markov) || length(operator.as.markov) != 1L) {
        stop("operator.as.markov must be TRUE/FALSE.")
    }

    if (!is.numeric(potential.eps) || length(potential.eps) != 1L ||
        !is.finite(potential.eps) || potential.eps <= 0) {
        stop("potential.eps must be a single positive finite number.")
    }

    if (!is.numeric(edge.length.eps) || length(edge.length.eps) != 1L ||
        !is.finite(edge.length.eps) || edge.length.eps <= 0) {
        stop("edge.length.eps must be a single positive finite number.")
    }

    if (!is.null(K) && !is.null(P)) {
        stop("Provide at most one of K or P.")
    }

    operator.info <- .pmg_prepare_operator_inputs(K = K,
                                                  P = P,
                                                  operator.as.markov = operator.as.markov,
                                                  potential.eps = potential.eps)

    D.target <- NULL
    if (!is.null(D.pot)) {
        D.target <- .pmg_as_square_numeric_matrix(D.pot, "D.pot")
    }

    U.target <- NULL
    if (!is.null(U.pot)) {
        U.target <- .pmg_as_numeric_matrix(U.pot, "U.pot")
    }

    inferred.n <- c()
    if (!is.null(adj.list)) inferred.n <- c(inferred.n, length(adj.list))
    if (!is.null(D.target)) inferred.n <- c(inferred.n, nrow(D.target))
    if (!is.null(U.target)) inferred.n <- c(inferred.n, nrow(U.target))
    if (!is.null(operator.info$operator)) inferred.n <- c(inferred.n, nrow(operator.info$operator))

    if (length(inferred.n) == 0L) {
        stop("Provide at least one of: adj.list, D.pot, U.pot, K, or P.")
    }

    inferred.n <- unique(inferred.n)
    if (length(inferred.n) != 1L) {
        stop("Input dimensions are inconsistent across adj.list / D.pot / U.pot / K / P.")
    }
    n.vertices <- inferred.n[1L]

    if (!is.null(U.target) && nrow(U.target) != n.vertices) {
        stop("nrow(U.pot) must match graph size.")
    }
    if (!is.null(D.target) && any(dim(D.target) != c(n.vertices, n.vertices))) {
        stop("D.pot must be an n x n matrix matching graph size.")
    }

    derived.potential.from.operator <- FALSE
    if (is.null(U.target) && is.null(D.target) && !is.null(operator.info$P.diffusion)) {
        Pt <- .pmg_matrix_power_int(operator.info$P.diffusion, diffusion.time)
        U.target <- -log(pmax(Pt, potential.eps))
        derived.potential.from.operator <- TRUE
    }

    topology.method.used <- topology.method
    topology.info <- NULL

    if (is.null(adj.list)) {
        if (topology.method.used == "auto") {
            if (!is.null(operator.info$operator)) {
                topology.method.used <- "operator_nonzero"
            } else {
                topology.method.used <- "potential_knn"
            }
        }

        if (topology.method.used == "operator_nonzero") {
            if (is.null(operator.info$operator)) {
                stop("topology.method = 'operator_nonzero' requires K or P.")
            }
            topology.info <- .pmg_build_topology_from_operator_nonzero(
                operator = operator.info$operator,
                min.operator.weight = min.operator.weight,
                edge.length.transform = edge.length.transform,
                edge.length.eps = edge.length.eps
            )
        } else if (topology.method.used == "operator_knn") {
            if (is.null(operator.info$operator)) {
                stop("topology.method = 'operator_knn' requires K or P.")
            }
            topology.info <- .pmg_build_topology_from_operator_knn(
                operator = operator.info$operator,
                k.topology = k.topology,
                knn.symmetrize = knn.symmetrize,
                min.operator.weight = min.operator.weight,
                edge.length.transform = edge.length.transform,
                edge.length.eps = edge.length.eps
            )
        } else if (topology.method.used == "potential_knn") {
            if (is.null(U.target) && is.null(D.target)) {
                stop("topology.method = 'potential_knn' requires U.pot, D.pot, K, or P.")
            }
            topology.info <- .pmg_build_topology_from_potential_knn(
                U.pot = U.target,
                D.pot = D.target,
                k.topology = k.topology,
                knn.symmetrize = knn.symmetrize,
                edge.length.eps = edge.length.eps
            )
        } else {
            stop("Unsupported topology.method.")
        }

        adj.list <- topology.info$adj.list
        if (is.null(weight.list)) {
            weight.list <- topology.info$weight.list
        }
    } else {
        topology.method.used <- "provided_adj"
        .pmg_validate_graph_inputs(adj.list, weight.list)
    }

    fit <- construct.potential.metric.graph(
        adj.list = adj.list,
        weight.list = weight.list,
        D.pot = D.target,
        U.pot = U.target,
        pair.quantiles = pair.quantiles,
        pair.mode = pair.mode,
        n.landmarks = n.landmarks,
        max.pairs = max.pairs,
        long.range.power = long.range.power,
        eta.blend = eta.blend,
        lambda.local = lambda.local,
        lambda.ridge = lambda.ridge,
        max.iter = max.iter,
        tol = tol,
        min.edge.length = min.edge.length,
        seed = seed,
        verbose = verbose
    )

    fit$wrapper_info <- list(
        topology_method_requested = topology.method,
        topology_method_used = topology.method.used,
        operator_source = operator.info$source,
        diffusion_time = diffusion.time,
        derived_potential_from_operator = derived.potential.from.operator,
        knn_symmetrize = knn.symmetrize,
        k_topology = k.topology,
        n_vertices = n.vertices,
        n_edges_topology = if (!is.null(topology.info)) topology.info$n.edges else fit$diagnostics$n_edges
    )

    class(fit) <- c("potential_metric_graph_phate", class(fit))
    return(fit)
}


.pmg_as_square_numeric_matrix <- function(X, name = "matrix") {
    if (inherits(X, "dist")) {
        M <- as.matrix(X)
    } else {
        M <- as.matrix(X)
    }
    if (!is.numeric(M) || length(dim(M)) != 2L || nrow(M) != ncol(M)) {
        stop(sprintf("%s must be a numeric square matrix (or dist).", name))
    }
    if (any(!is.finite(M))) {
        stop(sprintf("%s contains non-finite values.", name))
    }
    return(M)
}


.pmg_as_numeric_matrix <- function(X, name = "matrix") {
    M <- as.matrix(X)
    if (!is.numeric(M) || length(dim(M)) != 2L) {
        stop(sprintf("%s must be a numeric matrix.", name))
    }
    if (any(!is.finite(M))) {
        stop(sprintf("%s contains non-finite values.", name))
    }
    return(M)
}


.pmg_prepare_operator_inputs <- function(K,
                                         P,
                                         operator.as.markov = TRUE,
                                         potential.eps = 1e-12) {
    if (is.null(K) && is.null(P)) {
        return(list(operator = NULL, P.diffusion = NULL, source = NULL))
    }

    if (!is.null(K)) {
        operator <- .pmg_as_square_numeric_matrix(K, "K")
        source <- "K"
    } else {
        operator <- .pmg_as_square_numeric_matrix(P, "P")
        source <- "P"
    }

    if (any(operator < 0)) {
        stop(sprintf("%s must be non-negative.", source))
    }

    P.diffusion <- operator
    if (isTRUE(operator.as.markov)) {
        P.diffusion <- .pmg_row_normalize_operator(P.diffusion,
                                                   eps = potential.eps)
    }

    return(list(operator = operator, P.diffusion = P.diffusion, source = source))
}


.pmg_row_normalize_operator <- function(M, eps = 1e-12) {
    n <- nrow(M)
    rs <- rowSums(M)
    out <- matrix(0, nrow = n, ncol = n)
    good <- is.finite(rs) & rs > eps

    if (any(good)) {
        out[good, ] <- M[good, , drop = FALSE] / rs[good]
    }

    if (any(!good)) {
        idx <- which(!good)
        out[cbind(idx, idx)] <- 1
    }

    return(out)
}


.pmg_matrix_power_int <- function(M, t) {
    if (t <= 1L) {
        return(M)
    }

    out <- M
    for (iter in 2:t) {
        out <- out %*% M
    }
    return(out)
}


.pmg_similarity_to_edge_length <- function(similarity,
                                           transform = c("neglog", "reciprocal"),
                                           eps = 1e-12) {
    transform <- match.arg(transform)
    s <- pmax(as.numeric(similarity), eps)

    if (transform == "neglog") {
        return(-log(s))
    }
    return(1 / s)
}


.pmg_build_topology_from_operator_nonzero <- function(operator,
                                                      min.operator.weight = 0,
                                                      edge.length.transform = "neglog",
                                                      edge.length.eps = 1e-12) {
    n <- nrow(operator)
    S <- 0.5 * (operator + t(operator))
    diag(S) <- 0

    edge.idx <- which(upper.tri(S) & S > min.operator.weight, arr.ind = TRUE)
    if (nrow(edge.idx) < 1L) {
        stop("No edges retained from operator_nonzero topology. Lower min.operator.weight or use a kNN topology.")
    }

    edge.mat <- cbind(edge.idx[, 1L], edge.idx[, 2L])
    sim <- S[cbind(edge.mat[, 1L], edge.mat[, 2L])]
    edge.length <- .pmg_similarity_to_edge_length(similarity = sim,
                                                  transform = edge.length.transform,
                                                  eps = edge.length.eps)

    graph.lists <- .pmg_edges_to_lists(edge.mat = edge.mat,
                                       edge.weights = edge.length,
                                       n.vertices = n)

    return(list(
        adj.list = graph.lists$adj_list,
        weight.list = graph.lists$weight_list,
        n.edges = nrow(edge.mat)
    ))
}


.pmg_build_topology_from_operator_knn <- function(operator,
                                                  k.topology = 15L,
                                                  knn.symmetrize = c("union", "mutual"),
                                                  min.operator.weight = 0,
                                                  edge.length.transform = "neglog",
                                                  edge.length.eps = 1e-12) {
    knn.symmetrize <- match.arg(knn.symmetrize)
    n <- nrow(operator)
    k.use <- min(as.integer(k.topology), n - 1L)
    if (k.use < 1L) {
        stop("k.topology must be at least 1 for operator_knn.")
    }

    dir.mat <- matrix(FALSE, nrow = n, ncol = n)
    for (i in seq_len(n)) {
        row.i <- as.numeric(operator[i, ])
        row.i[i] <- -Inf
        ord <- order(row.i, decreasing = TRUE)
        ord <- ord[is.finite(row.i[ord]) & row.i[ord] > min.operator.weight]
        if (length(ord) > 0L) {
            keep <- ord[seq_len(min(k.use, length(ord)))]
            dir.mat[i, keep] <- TRUE
        }
    }

    undir <- if (knn.symmetrize == "mutual") {
        dir.mat & t(dir.mat)
    } else {
        dir.mat | t(dir.mat)
    }

    edge.idx <- which(upper.tri(undir) & undir, arr.ind = TRUE)
    if (nrow(edge.idx) < 1L) {
        stop("No edges retained from operator_knn topology. Increase k.topology or lower min.operator.weight.")
    }

    edge.mat <- cbind(edge.idx[, 1L], edge.idx[, 2L])
    sim.uv <- operator[cbind(edge.mat[, 1L], edge.mat[, 2L])]
    sim.vu <- operator[cbind(edge.mat[, 2L], edge.mat[, 1L])]
    if (knn.symmetrize == "mutual") {
        sim <- 0.5 * (sim.uv + sim.vu)
    } else {
        sim <- pmax(sim.uv, sim.vu)
    }
    sim <- pmax(sim, edge.length.eps)

    edge.length <- .pmg_similarity_to_edge_length(similarity = sim,
                                                  transform = edge.length.transform,
                                                  eps = edge.length.eps)

    graph.lists <- .pmg_edges_to_lists(edge.mat = edge.mat,
                                       edge.weights = edge.length,
                                       n.vertices = n)

    return(list(
        adj.list = graph.lists$adj_list,
        weight.list = graph.lists$weight_list,
        n.edges = nrow(edge.mat)
    ))
}


.pmg_build_topology_from_potential_knn <- function(U.pot = NULL,
                                                   D.pot = NULL,
                                                   k.topology = 15L,
                                                   knn.symmetrize = c("union", "mutual"),
                                                   edge.length.eps = 1e-12) {
    knn.symmetrize <- match.arg(knn.symmetrize)

    if (is.null(U.pot) && is.null(D.pot)) {
        stop("Provide U.pot or D.pot for potential_knn topology.")
    }

    if (!is.null(U.pot)) {
        U <- .pmg_as_numeric_matrix(U.pot, "U.pot")
        n <- nrow(U)
    } else {
        D <- .pmg_as_square_numeric_matrix(D.pot, "D.pot")
        n <- nrow(D)
        U <- NULL
    }

    if (is.null(D.pot)) {
        D <- as.matrix(stats::dist(U))
    } else {
        D <- .pmg_as_square_numeric_matrix(D.pot, "D.pot")
    }
    D <- 0.5 * (D + t(D))
    diag(D) <- 0

    k.use <- min(as.integer(k.topology), n - 1L)
    if (k.use < 1L) {
        stop("k.topology must be at least 1 for potential_knn.")
    }

    dir.mat <- matrix(FALSE, nrow = n, ncol = n)
    if (!is.null(U)) {
        knn <- FNN::get.knn(U, k = k.use)
        for (i in seq_len(n)) {
            dir.mat[i, knn$nn.index[i, ]] <- TRUE
        }
    } else {
        for (i in seq_len(n)) {
            row.i <- D[i, ]
            row.i[i] <- Inf
            ord <- order(row.i, decreasing = FALSE)
            keep <- ord[seq_len(k.use)]
            dir.mat[i, keep] <- TRUE
        }
    }

    undir <- if (knn.symmetrize == "mutual") {
        dir.mat & t(dir.mat)
    } else {
        dir.mat | t(dir.mat)
    }

    edge.idx <- which(upper.tri(undir) & undir, arr.ind = TRUE)
    if (nrow(edge.idx) < 1L) {
        stop("No edges retained from potential_knn topology. Increase k.topology.")
    }

    edge.mat <- cbind(edge.idx[, 1L], edge.idx[, 2L])
    edge.length <- D[cbind(edge.mat[, 1L], edge.mat[, 2L])]
    edge.length <- pmax(edge.length, edge.length.eps)

    graph.lists <- .pmg_edges_to_lists(edge.mat = edge.mat,
                                       edge.weights = edge.length,
                                       n.vertices = n)

    return(list(
        adj.list = graph.lists$adj_list,
        weight.list = graph.lists$weight_list,
        n.edges = nrow(edge.mat)
    ))
}


.pmg_validate_graph_inputs <- function(adj.list, weight.list = NULL) {
    if (!is.list(adj.list) || length(adj.list) < 2L) {
        stop("adj.list must be a list with at least 2 vertices.")
    }

    n <- length(adj.list)
    for (i in seq_len(n)) {
        nbr <- adj.list[[i]]
        if (!is.numeric(nbr)) {
            stop(sprintf("adj.list[[%d]] must be numeric/integer.", i))
        }
        if (length(nbr) > 0L) {
            if (any(!is.finite(nbr)) || any(nbr != floor(nbr))) {
                stop(sprintf("adj.list[[%d]] must contain finite integer indices.", i))
            }
            if (any(nbr < 1L) || any(nbr > n)) {
                stop(sprintf("adj.list[[%d]] has indices outside [1, %d].", i, n))
            }
        }
    }

    if (!is.null(weight.list)) {
        if (!is.list(weight.list) || length(weight.list) != n) {
            stop("weight.list must be NULL or a list with the same length as adj.list.")
        }
        for (i in seq_len(n)) {
            if (length(weight.list[[i]]) != length(adj.list[[i]])) {
                stop(sprintf("weight.list[[%d]] must match length(adj.list[[%d]]).", i, i))
            }
            if (length(weight.list[[i]]) > 0L) {
                w <- as.numeric(weight.list[[i]])
                if (any(!is.finite(w)) || any(w <= 0)) {
                    stop(sprintf("weight.list[[%d]] must contain finite positive values.", i))
                }
            }
        }
    }
}


.pmg_resolve_target_distance_matrix <- function(D.pot, U.pot, n.vertices) {
    if (is.null(D.pot)) {
        if (is.null(U.pot)) {
            stop("Provide either D.pot or U.pot.")
        }
        U <- as.matrix(U.pot)
        if (!is.numeric(U) || nrow(U) != n.vertices) {
            stop("U.pot must be numeric with nrow(U.pot) == length(adj.list).")
        }
        if (any(!is.finite(U))) {
            stop("U.pot contains non-finite values.")
        }
        D <- as.matrix(stats::dist(U))
    } else {
        if (inherits(D.pot, "dist")) {
            D <- as.matrix(D.pot)
        } else {
            D <- as.matrix(D.pot)
        }
        if (!is.numeric(D) || any(dim(D) != c(n.vertices, n.vertices))) {
            stop("D.pot must be a numeric n x n matrix (or dist object).")
        }
        if (any(!is.finite(D))) {
            stop("D.pot contains non-finite values.")
        }
    }

    D <- 0.5 * (D + t(D))
    diag(D) <- 0

    if (max(abs(D - t(D))) > 1e-8) {
        stop("D.pot must be symmetric.")
    }
    if (any(D < 0)) {
        stop("D.pot must be non-negative.")
    }

    return(D)
}


.pmg_extract_undirected_edges <- function(adj.list, weight.list = NULL) {
    n <- length(adj.list)
    key_env <- new.env(hash = TRUE, parent = emptyenv())
    edge_u <- integer(0)
    edge_v <- integer(0)
    edge_w <- numeric(0)

    for (i in seq_len(n)) {
        nbr <- as.integer(adj.list[[i]])
        has_weights <- !is.null(weight.list)
        if (has_weights) {
            w <- as.numeric(weight.list[[i]])
        } else {
            w <- NULL
        }

        if (length(nbr) == 0L) next
        for (k in seq_along(nbr)) {
            j <- nbr[k]
            if (i == j) next
            u <- min(i, j)
            v <- max(i, j)
            key <- paste0(u, "-", v)

            if (!exists(key, envir = key_env, inherits = FALSE)) {
                idx <- length(edge_u) + 1L
                assign(key, idx, envir = key_env)
                edge_u[idx] <- u
                edge_v[idx] <- v
                if (has_weights) {
                    edge_w[idx] <- w[k]
                }
            } else if (has_weights) {
                idx <- get(key, envir = key_env, inherits = FALSE)
                edge_w[idx] <- 0.5 * (edge_w[idx] + w[k])
            }
        }
    }

    if (length(edge_u) == 0L) {
        return(list(edge_matrix = matrix(integer(0), ncol = 2), edge_weights = NULL))
    }

    ord <- order(edge_u, edge_v)
    edge.mat <- cbind(edge_u[ord], edge_v[ord])

    if (!is.null(weight.list)) {
        edge.w <- edge_w[ord]
    } else {
        edge.w <- NULL
    }

    return(list(edge_matrix = edge.mat, edge_weights = edge.w))
}


.pmg_select_pairs <- function(D,
                              pair.quantiles,
                              pair.mode = c("landmark", "all"),
                              n.landmarks = 128L,
                              max.pairs = 5000L) {
    pair.mode <- match.arg(pair.mode)
    n <- nrow(D)

    upper.mask <- upper.tri(D, diag = FALSE)
    vals <- D[upper.mask]
    valid <- is.finite(vals) & vals > 0
    if (!any(valid)) {
        stop("No positive finite pairwise target distances were found.")
    }

    qvals <- as.numeric(stats::quantile(vals[valid],
                                        probs = pair.quantiles,
                                        na.rm = TRUE,
                                        names = FALSE))

    cand <- which(upper.mask & D >= qvals[1] & D <= qvals[2], arr.ind = TRUE)
    if (nrow(cand) == 0L) {
        cand <- which(upper.mask & D > 0, arr.ind = TRUE)
    }

    if (pair.mode == "landmark") {
        nL <- min(as.integer(n.landmarks), n)
        landmarks <- unique(as.integer(round(seq(1, n, length.out = nL))))
        keep <- cand[, 1L] %in% landmarks | cand[, 2L] %in% landmarks
        cand <- cand[keep, , drop = FALSE]
        if (nrow(cand) == 0L) {
            cand <- which(upper.mask & D > 0, arr.ind = TRUE)
        }
    }

    n.pairs.candidate <- nrow(cand)
    if (n.pairs.candidate < 1L) {
        stop("No pair candidates available for fitting.")
    }

    if (nrow(cand) > max.pairs) {
        d.cand <- D[cbind(cand[, 1L], cand[, 2L])]
        prob <- d.cand / sum(d.cand)
        idx <- sample.int(nrow(cand), size = max.pairs, replace = FALSE, prob = prob)
        cand <- cand[idx, , drop = FALSE]
    }

    storage.mode(cand) <- "integer"

    return(list(
        pairs = cand,
        n_pairs_candidate = n.pairs.candidate,
        distance_quantile_range = qvals
    ))
}


.pmg_graph_from_edges <- function(edge.mat, edge.weights, n.vertices) {
    g <- igraph::make_empty_graph(n = n.vertices, directed = FALSE)
    if (nrow(edge.mat) > 0L) {
        g <- igraph::add_edges(g, as.vector(t(edge.mat)))
        igraph::E(g)$weight <- edge.weights
    }
    return(g)
}


.pmg_build_incidence <- function(edge.mat, edge.weights, n.vertices, pairs) {
    n.edges <- nrow(edge.mat)
    n.pairs <- nrow(pairs)

    g <- .pmg_graph_from_edges(edge.mat, edge.weights, n.vertices)

    edge.keys <- paste(edge.mat[, 1L], edge.mat[, 2L], sep = "-")
    edge.index <- setNames(seq_len(n.edges), edge.keys)

    A <- matrix(0, nrow = n.pairs, ncol = n.edges)
    d.pred <- numeric(n.pairs)

    for (r in seq_len(n.pairs)) {
        i <- pairs[r, 1L]
        j <- pairs[r, 2L]

        sp <- igraph::shortest_paths(g,
                                     from = i,
                                     to = j,
                                     weights = igraph::E(g)$weight,
                                     output = "vpath")
        vpath <- sp$vpath[[1L]]

        if (length(vpath) == 0L) {
            d.pred[r] <- Inf
            next
        }

        verts <- as.integer(igraph::as_ids(vpath))
        if (length(verts) < 2L) {
            d.pred[r] <- 0
            next
        }

        u <- pmin(verts[-length(verts)], verts[-1L])
        v <- pmax(verts[-length(verts)], verts[-1L])
        keys <- paste(u, v, sep = "-")
        idx <- unname(edge.index[keys])

        if (any(is.na(idx))) {
            stop("Internal error: failed to map shortest-path edges to edge index.")
        }

        A[r, idx] <- A[r, idx] + 1
        d.pred[r] <- sum(edge.weights[idx])
    }

    if (any(!is.finite(d.pred))) {
        stop("Graph became disconnected under current edge weights.")
    }

    return(list(incidence = A, d.pred = d.pred))
}


.pmg_compute_objective <- function(d.pred,
                                   d.target,
                                   pair.weights,
                                   edge.weights,
                                   d.edge.target,
                                   lambda.local,
                                   lambda.ridge) {
    residual <- d.pred - d.target
    fit.term <- sum(pair.weights * residual^2)
    local.term <- lambda.local * sum((edge.weights - d.edge.target)^2)
    ridge.term <- lambda.ridge * sum(edge.weights^2)
    objective <- fit.term + local.term + ridge.term

    denom <- sum(pair.weights * d.target^2)
    if (denom <= 0) {
        stress <- sqrt(mean(residual^2))
    } else {
        stress <- sqrt(fit.term / denom)
    }

    rho <- suppressWarnings(stats::cor(d.pred, d.target, method = "spearman"))
    if (is.na(rho) || !is.finite(rho)) {
        rho <- 0
    }

    return(list(
        objective = objective,
        fit_term = fit.term,
        local_term = local.term,
        ridge_term = ridge.term,
        stress_weighted = stress,
        corr_spearman = rho
    ))
}


.pmg_solve_fixed_paths <- function(A,
                                   d.target,
                                   pair.weights,
                                   d.edge.target,
                                   lambda.local,
                                   lambda.ridge,
                                   min.edge.length) {
    ws <- sqrt(pair.weights)
    A.weighted <- A * ws
    b.weighted <- d.target * ws

    H <- crossprod(A.weighted)
    reg <- lambda.local + lambda.ridge
    if (reg <= 0) reg <- 1e-8
    H <- H + diag(reg, ncol(A))

    rhs <- crossprod(A.weighted, b.weighted) + lambda.local * d.edge.target

    sol <- tryCatch(
        as.numeric(solve(H, rhs)),
        error = function(e) as.numeric(qr.solve(H, rhs))
    )

    sol[!is.finite(sol)] <- d.edge.target[!is.finite(sol)]
    sol <- pmax(min.edge.length, sol)
    return(sol)
}


.pmg_evaluate_weights <- function(edge.mat,
                                  edge.weights,
                                  n.vertices,
                                  pairs,
                                  d.target,
                                  pair.weights,
                                  d.edge.target,
                                  lambda.local,
                                  lambda.ridge) {
    inc <- .pmg_build_incidence(edge.mat = edge.mat,
                                edge.weights = edge.weights,
                                n.vertices = n.vertices,
                                pairs = pairs)

    obj <- .pmg_compute_objective(d.pred = inc$d.pred,
                                  d.target = d.target,
                                  pair.weights = pair.weights,
                                  edge.weights = edge.weights,
                                  d.edge.target = d.edge.target,
                                  lambda.local = lambda.local,
                                  lambda.ridge = lambda.ridge)

    return(c(inc, obj))
}


.pmg_edges_to_lists <- function(edge.mat, edge.weights, n.vertices, vertex.names = NULL) {
    adj.list <- vector("list", n.vertices)
    weight.list <- vector("list", n.vertices)

    for (e in seq_len(nrow(edge.mat))) {
        i <- edge.mat[e, 1L]
        j <- edge.mat[e, 2L]
        w <- edge.weights[e]

        adj.list[[i]] <- c(adj.list[[i]], j)
        weight.list[[i]] <- c(weight.list[[i]], w)

        adj.list[[j]] <- c(adj.list[[j]], i)
        weight.list[[j]] <- c(weight.list[[j]], w)
    }

    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) > 1L) {
            ord <- order(adj.list[[i]])
            adj.list[[i]] <- as.integer(adj.list[[i]][ord])
            weight.list[[i]] <- as.numeric(weight.list[[i]][ord])
        } else {
            adj.list[[i]] <- as.integer(adj.list[[i]])
            weight.list[[i]] <- as.numeric(weight.list[[i]])
        }
    }

    if (!is.null(vertex.names) && length(vertex.names) == n.vertices) {
        names(adj.list) <- vertex.names
        names(weight.list) <- vertex.names
    }

    return(list(adj_list = adj.list, weight_list = weight.list))
}

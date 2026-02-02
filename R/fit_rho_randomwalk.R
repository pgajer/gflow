#' Fit a Lazy Random-Walk Smoother to a Target Vertex Density
#'
#' @description
#' Fits a simple surrogate mass/density evolution model on a weighted graph by
#' applying a lazy random-walk smoothing operator \eqn{P^m} to an initial vertex
#' density \eqn{\rho_0} and choosing smoothing parameters to best match a target
#' density \eqn{\rho_T}. This is intended as a robust, eigen-free preconditioner
#' or diagnostic for graph-based density evolution pipelines.
#'
#' The graph is supplied as an adjacency/weight-list pair, where each vertex
#' \eqn{i} has neighbors \code{adj.list[[i]]} with corresponding nonnegative
#' weights \code{weight.list[[i]]}. The transition is a lazy random walk
#' \eqn{P = (1-\eta)I + \eta D^{-1}W}, where \eqn{D} is the weighted degree matrix.
#'
#' @details
#' The routine performs a grid search over the number of random-walk steps
#' \code{m} and laziness \code{eta}. For each parameter pair, the algorithm:
#' \enumerate{
#'   \item Takes \eqn{\rho_0 =} \code{rho.vertex[[1]]} and \eqn{\rho_T =} \code{rho.vertex[[11]]}.
#'   \item Applies \eqn{P^m} to \eqn{\rho_0} using repeated sparse neighborhood aggregation
#'         without forming an explicit \eqn{n\times n} matrix.
#'   \item Optionally floors densities at \code{min.rho} and renormalizes to sum \eqn{n}.
#'   \item Computes a loss (RMSE) to \eqn{\rho_T}.
#' }
#'
#' As a simple dimensionless smoothing scale, the function can estimate a robust
#' distance scale using unweighted BFS shortest paths and report \code{m/scale}.
#' The default scale is the 95th percentile of sampled shortest-path distances.
#'
#' @param rho.vertex List of numeric vectors, each of length \eqn{n}, containing
#'   vertex density vectors at successive iterations. The function uses
#'   \code{rho.vertex[[1]]} as the initial density and \code{rho.vertex[[11]]} as
#'   the target density.
#' @param adj.list List of integer vectors giving neighbor indices for each
#'   vertex. Indices must be in \code{1:n}. Length must equal \eqn{n}.
#' @param weight.list List of numeric vectors of the same structure as
#'   \code{adj.list}, giving nonnegative edge weights for each neighbor.
#' @param m.grid Integer vector of candidate random-walk step counts \eqn{m \ge 0}.
#' @param eta.grid Numeric vector of candidate laziness parameters
#'   \eqn{\eta \in (0, 1]}. Larger values produce stronger smoothing per step.
#' @param distance.scale Character scalar specifying the distance scale used to
#'   report a dimensionless smoothing ratio. One of \code{"p95_geodesic"} (default),
#'   \code{"diameter"}, or \code{"none"}. Distances are estimated using unweighted
#'   BFS sampling.
#' @param sample.pairs Integer scalar. Number of sampled source--target pairs
#'   used to estimate the distance scale when \code{distance.scale != "none"}.
#' @param seed Integer scalar random seed used for distance-scale sampling.
#' @param normalize.sum.to.n Logical; if \code{TRUE} (default), renormalize each
#'   smoothed density so that \code{sum(rho) == n}.
#' @param min.rho Numeric scalar floor applied to densities before renormalization.
#' @param verbose Logical; if \code{TRUE}, print progress and the estimated scale.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{rho0}}{Initial density vector \code{rho.vertex[[1]]}, normalized if requested.}
#'   \item{\code{rhoT}}{Target density vector \code{rho.vertex[[11]]}, normalized if requested.}
#'   \item{\code{scale}}{List with \code{type} and \code{value} for the distance scale estimate.}
#'   \item{\code{grid}}{Data frame of all tested parameter pairs with RMSE and correlation,
#'     sorted by RMSE (best first).}
#'   \item{\code{best}}{List with best \code{m}, \code{eta}, \code{rmse}, \code{corr}, and
#'     \code{m_over_scale}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## rho.vertex is typically taken from a fitted object:
#' ## rho.vertex <- fit$density$rho.vertex
#'
#' res <- fit.rho.randomwalk(
#'   rho.vertex = rho.vertex,
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   m.grid = 0:40,
#'   eta.grid = c(0.1, 0.2, 0.3, 0.5, 0.75, 1),
#'   distance.scale = "p95_geodesic",
#'   sample.pairs = 5000L
#' )
#'
#' res$best
#' head(res$grid, 10)
#' }
#'
#' @seealso \code{\link[igraph]{cluster_louvain}} for community detection and
#'   \code{\link[stats]{cor}} for correlation.
#'
#' @export
fit.rho.randomwalk <- function(rho.vertex,
                               adj.list,
                               weight.list,
                               m.grid = 0:30,
                               eta.grid = c(0.1, 0.2, 0.3, 0.5, 0.75, 1.0),
                               distance.scale = c("p95_geodesic", "diameter", "none"),
                               sample.pairs = 5000L,
                               seed = 1L,
                               normalize.sum.to.n = TRUE,
                               min.rho = 0,
                               verbose = TRUE) {
    ## rho.vertex: list of numeric vectors (length n), e.g. J.fit$density$rho.vertex
    ## adj.list, weight.list: adjacency + weights per vertex (same length n)
    ## m.grid: integer steps for P^m
    ## eta.grid: laziness parameter eta in (0,1]; P = (1-eta)I + eta D^{-1}W
    ## distance.scale: how to compute dimensionless scale
    ## sample.pairs: number of random pairs to estimate p95 shortest path
    ## normalize.sum.to.n: enforce sum(rho)=n after each smoothing
    ## min.rho: floor (0 recommended; you can also use 1e-8)
    ## Returns: list with best params, diagnostics, and trajectories

    distance.scale <- match.arg(distance.scale)

    stopifnot(is.list(rho.vertex), length(rho.vertex) >= 11L)
    rho0 <- rho.vertex[[1L]]
    rhoT <- rho.vertex[[11L]]

    n <- length(rho0)
    stopifnot(length(rhoT) == n)
    stopifnot(length(adj.list) == n)
    stopifnot(length(weight.list) == n)

    ## ------------------------------
    ## Basic validation + cleaning
    ## ------------------------------
    if (any(!is.finite(rho0)) || any(!is.finite(rhoT))) stop("rho0/rhoT must be finite.")
    if (any(rho0 < 0) || any(rhoT < 0)) stop("rho0/rhoT must be nonnegative.")

    ## normalize rho0 and rhoT for fair comparison
    normalize.rho <- function(rho) {
        rho <- pmax(rho, min.rho)
        if (normalize.sum.to.n) {
            s <- sum(rho)
            if (!(s > 0)) stop("rho has nonpositive sum after flooring.")
            rho <- rho * (n / s)
        }
        rho
    }

    rho0 <- normalize.rho(rho0)
    rhoT <- normalize.rho(rhoT)

    ## ------------------------------
    ## Build row-stochastic transition application: v <- P v
    ## (lazy RW): P = (1-eta)I + eta D^{-1}W
    ## ------------------------------
    ## We avoid forming an n x n matrix; we apply it as an operator for speed.
    deg.w <- numeric(n)
    for (i in seq_len(n)) {
        wi <- weight.list[[i]]
        if (length(wi) != length(adj.list[[i]])) stop("weight.list[[i]] length mismatch.")
        if (any(!is.finite(wi)) || any(wi < 0)) stop("weights must be finite and nonnegative.")
        deg.w[i] <- sum(wi)
    }
    if (any(deg.w <= 0)) stop("Some vertices have zero weighted degree; RW undefined.")

    apply.P <- function(v, eta) {
        ## returns P v
        out <- numeric(n)
        for (i in seq_len(n)) {
            nb <- adj.list[[i]]
            w <- weight.list[[i]]
            ## weighted neighbor average
            out[i] <- sum(w * v[nb]) / deg.w[i]
        }
        (1 - eta) * v + eta * out
    }

    ## ------------------------------
    ## Distance scale for dimensionless bandwidth
    ## Use unweighted shortest paths estimated via BFS sampling
    ## ------------------------------
    estimate.scale <- function() {
        if (distance.scale == "none") return(NA_real_)

        set.seed(seed)
        s.idx <- sample.int(n, size = min(n, max(10L, as.integer(sqrt(sample.pairs)))), replace = FALSE)

        ## BFS distances from one source (unweighted)
        bfs.dist <- function(src) {
            dist <- rep.int(NA_integer_, n)
            dist[src] <- 0L
            q <- integer(n)
            head <- 1L; tail <- 1L
            q[tail] <- src

            while (head <= tail) {
                v <- q[head]; head <- head + 1L
                dv <- dist[v]
                nb <- adj.list[[v]]
                for (u in nb) {
                    if (is.na(dist[u])) {
                        dist[u] <- dv + 1L
                        tail <- tail + 1L
                        q[tail] <- u
                    }
                }
            }
            dist
        }

        ## sample pair distances by BFS from random sources
        d.samples <- integer(0L)
        for (src in s.idx) {
            d <- bfs.dist(src)
            ok <- which(!is.na(d))
            if (length(ok) <= 1L) next
            ## sample targets
            t.idx <- sample(ok, size = min(length(ok), as.integer(ceiling(sample.pairs / length(s.idx)))), replace = TRUE)
            d.samples <- c(d.samples, d[t.idx])
        }

        d.samples <- d.samples[d.samples > 0L]
        if (length(d.samples) == 0L) return(NA_real_)

        if (distance.scale == "diameter") {
            return(as.double(max(d.samples)))
        }
        ## p95_geodesic
        as.double(stats::quantile(d.samples, probs = 0.95, type = 7, names = FALSE))
    }

    scale.val <- estimate.scale()

    ## ------------------------------
    ## Loss function
    ## ------------------------------
    rmse <- function(x, y) sqrt(mean((x - y)^2))

    ## ------------------------------
    ## Grid search
    ## ------------------------------
    results <- expand.grid(m = as.integer(m.grid),
                           eta = as.double(eta.grid),
                           stringsAsFactors = FALSE)

    results$rmse <- NA_real_
    results$corr <- NA_real_

    if (verbose) {
        message("Grid size: ", nrow(results), " (m x eta)")
        if (!is.na(scale.val)) message("Distance scale (", distance.scale, "): ", signif(scale.val, 4))
    }

    for (k in seq_len(nrow(results))) {
        m <- results$m[k]
        eta <- results$eta[k]

        v <- rho0
        if (m > 0L) {
            for (step in seq_len(m)) {
                v <- apply.P(v, eta)
                v <- normalize.rho(v)
            }
        }

        results$rmse[k] <- rmse(v, rhoT)
        results$corr[k] <- suppressWarnings(stats::cor(v, rhoT))
    }

    best.idx <- which.min(results$rmse)
    best <- results[best.idx, , drop = FALSE]

    ## Dimensionless scale suggestions
    best$m_over_scale <- if (!is.na(scale.val) && scale.val > 0) best$m / scale.val else NA_real_

    ## also report m/diameter if p95 used
    if (distance.scale == "p95_geodesic") {
        ## estimate a rough diameter too, from same samples
        ## (cheap: just use max sample as a proxy)
        ## If you want exact diameter, that’s expensive.
        ## Here we reuse scale.val for p95; diameter proxy omitted unless requested.
        best$m_over_diam_proxy <- NA_real_
    } else if (distance.scale == "diameter") {
        best$m_over_diam_proxy <- best$m_over_scale
    } else {
        best$m_over_diam_proxy <- NA_real_
    }

    list(
        rho0 = rho0,
        rhoT = rhoT,
        scale = list(type = distance.scale, value = scale.val),
        grid = results[order(results$rmse), ],
        best = list(
            m = best$m,
            eta = best$eta,
            rmse = best$rmse,
            corr = best$corr,
            m_over_scale = best$m_over_scale
        )
    )
}

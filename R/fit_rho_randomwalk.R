
#' Fit a Lazy Random-Walk Smoother to a Target Vertex Density
#'
#' @description
#' Fits a simple surrogate density evolution model on a weighted graph by
#' applying a lazy random-walk smoothing operator to an initial vertex density
#' \eqn{\rho_0} and choosing parameters to best match a target density
#' \eqn{\rho_T}. Optionally applies a post-smoothing power transform
#' \eqn{\rho \leftarrow \rho^{\kappa}} (followed by renormalization) to better
#' emulate nonlinear compression in downstream pipelines.
#'
#' The graph is supplied as an adjacency/weight-list pair, where each vertex
#' \eqn{i} has neighbors \code{adj.list[[i]]} with corresponding nonnegative
#' weights \code{weight.list[[i]]}. The transition is a lazy random walk
#' \eqn{P = (1-\eta)I + \eta D^{-1}W}, where \eqn{D} is the weighted degree matrix.
#'
#' @details
#' The function performs a grid search over:
#' \itemize{
#'   \item \code{m.grid}: number of random-walk steps \eqn{m}
#'   \item \code{eta.grid}: laziness parameter \eqn{\eta \in (0,1]}
#'   \item \code{kappa.grid}: post-smoothing power \eqn{\kappa > 0}
#' }
#'
#' After each random-walk step, the density is optionally floored by \code{min.rho}
#' and renormalized so that \code{sum(rho) = n} (where \eqn{n} is the number of
#' vertices). After \eqn{m} steps, the optional power transform is applied and the
#' density is renormalized again.
#'
#' The optimization objective can be selected via \code{objective}:
#' \itemize{
#'   \item \code{"rmse"}: minimize RMSE on the density scale
#'   \item \code{"corr"}: maximize Pearson correlation with the target
#'   \item \code{"rel_rmse"}: minimize relative RMSE using \eqn{(x-y)/(y+\varepsilon)}
#'   \item \code{"log_rmse"}: minimize RMSE on \eqn{\log(\rho+\varepsilon)}
#'   \item \code{"hybrid"}: maximize \code{corr - lambda * rmse}
#' }
#'
#' @param rho.vertex List of numeric vectors of length \eqn{n}. The function uses
#'   \code{rho.vertex[[1]]} as \eqn{\rho_0} and \code{rho.vertex[[11]]} as \eqn{\rho_T}.
#' @param adj.list List of integer vectors giving neighbors (indices in \code{1:n}).
#' @param weight.list List of numeric vectors of nonnegative edge weights aligned to \code{adj.list}.
#' @param m.grid Integer vector of candidate step counts \eqn{m \ge 0}.
#' @param eta.grid Numeric vector of candidate laziness parameters \eqn{\eta \in (0,1]}.
#' @param kappa.grid Numeric vector of candidate power-transform exponents \eqn{\kappa > 0}.
#'   Use \code{1} to disable the transform.
#' @param objective Character scalar specifying the optimization objective.
#' @param lambda Numeric scalar penalty used only when \code{objective = "hybrid"}.
#' @param eps.rel Small positive scalar used in relative RMSE denominators.
#' @param eps.log Small positive scalar used inside logarithms for \code{"log_rmse"}.
#' @param distance.scale Character scalar controlling how a graph distance scale is estimated.
#' @param sample.pairs Integer number of random source-target samples used to estimate \code{p95_geodesic}.
#' @param seed Integer seed used for distance-scale sampling.
#' @param normalize.sum.to.n Logical; if \code{TRUE}, enforce \code{sum(rho) = n}.
#' @param min.rho Nonnegative floor applied to densities before renormalization.
#' @param return.best.rho Logical; if \code{TRUE}, also return the best fitted density.
#' @param verbose Logical; if \code{TRUE}, print basic progress information.
#'
#' @return A list with elements:
#' \describe{
#'   \item{rho0}{Normalized initial density.}
#'   \item{rhoT}{Normalized target density.}
#'   \item{scale}{List with the distance scale type and estimated value.}
#'   \item{grid}{Data frame of grid results, sorted by the chosen objective.}
#'   \item{best}{List with best parameters and metrics under the chosen objective.}
#'   \item{rho.best}{(Optional) Best fitted density under the chosen objective.}
#' }
#'
#' @examples
#' ## Not run: requires an adjacency/weight list pair
#' ## res <- fit.rho.randomwalk(rho.vertex, adj.list, weight.list,
#' ##                           m.grid = 0:300,
#' ##                           eta.grid = c(0.2, 0.5, 1.0),
#' ##                           kappa.grid = seq(0.5, 2, by = 0.1),
#' ##                           objective = "corr")
#' ## End(Not run)
#'
#' @export
fit.rho.randomwalk <- function(rho.vertex,
                               adj.list,
                               weight.list,
                               m.grid = 0:30,
                               eta.grid = c(0.1, 0.2, 0.3, 0.5, 0.75, 1.0),
                               kappa.grid = 1,
                               objective = c("rmse", "corr", "rel_rmse", "log_rmse", "hybrid"),
                               lambda = 1,
                               eps.rel = 1e-8,
                               eps.log = 1e-8,
                               distance.scale = c("p95_geodesic", "diameter", "none"),
                               sample.pairs = 5000L,
                               seed = 1L,
                               normalize.sum.to.n = TRUE,
                               min.rho = 0,
                               return.best.rho = TRUE,
                               verbose = TRUE) {

    objective <- match.arg(objective)
    distance.scale <- match.arg(distance.scale)

    if (!is.list(rho.vertex) || length(rho.vertex) < 11L) {
        stop("rho.vertex must be a list with at least 11 entries (uses [[1]] and [[11]]).")
    }
    rho0 <- rho.vertex[[1L]]
    rhoT <- rho.vertex[[11L]]

    if (!is.numeric(rho0) || !is.numeric(rhoT)) stop("rho0/rhoT must be numeric vectors.")
    n <- length(rho0)
    if (n < 2L) stop("rho0 must have length >= 2.")
    if (length(rhoT) != n) stop("rhoT must have same length as rho0.")
    if (!is.list(adj.list) || length(adj.list) != n) stop("adj.list must be a list of length n.")
    if (!is.list(weight.list) || length(weight.list) != n) stop("weight.list must be a list of length n.")

    if (any(!is.finite(rho0)) || any(!is.finite(rhoT))) stop("rho0/rhoT must be finite.")
    if (any(rho0 < 0) || any(rhoT < 0)) stop("rho0/rhoT must be nonnegative.")

    if (!is.numeric(m.grid) || any(!is.finite(m.grid)) || any(m.grid < 0) || any(m.grid != as.integer(m.grid))) {
        stop("m.grid must be a nonnegative integer vector.")
    }
    m.grid <- as.integer(unique(m.grid))

    if (!is.numeric(eta.grid) || any(!is.finite(eta.grid)) || any(eta.grid <= 0) || any(eta.grid > 1)) {
        stop("eta.grid must be numeric in (0, 1].")
    }
    eta.grid <- as.double(unique(eta.grid))

    if (!is.numeric(kappa.grid) || any(!is.finite(kappa.grid)) || any(kappa.grid <= 0)) {
        stop("kappa.grid must be positive numeric.")
    }
    kappa.grid <- as.double(unique(kappa.grid))

    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop("lambda must be a single nonnegative finite number.")
    }
    if (!is.numeric(eps.rel) || length(eps.rel) != 1L || !is.finite(eps.rel) || eps.rel <= 0) {
        stop("eps.rel must be a single positive finite number.")
    }
    if (!is.numeric(eps.log) || length(eps.log) != 1L || !is.finite(eps.log) || eps.log <= 0) {
        stop("eps.log must be a single positive finite number.")
    }
    if (!is.numeric(min.rho) || length(min.rho) != 1L || !is.finite(min.rho) || min.rho < 0) {
        stop("min.rho must be a single nonnegative finite number.")
    }

    ## ------------------------------
    ## Normalization helper
    ## ------------------------------
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
    ## Build transition application: v <- P v
    ## (lazy RW): P = (1-eta)I + eta D^{-1}W
    ## ------------------------------
    deg.w <- numeric(n)
    for (i in seq_len(n)) {
        nb <- adj.list[[i]]
        w <- weight.list[[i]]

        if (!is.integer(nb)) nb <- as.integer(nb)
        if (!is.numeric(w)) w <- as.double(w)

        if (length(nb) != length(w)) stop("weight.list[[i]] length mismatch with adj.list[[i]].")
        if (length(nb) == 0L) stop("Some vertices have empty neighbor sets; RW undefined.")
        if (any(nb < 1L) || any(nb > n)) stop("adj.list indices must be in 1..n.")
        if (any(!is.finite(w)) || any(w < 0)) stop("weights must be finite and nonnegative.")

        deg.w[i] <- sum(w)
    }
    if (any(deg.w <= 0)) stop("Some vertices have zero weighted degree; RW undefined.")

    apply.P <- function(v, eta) {
        out <- numeric(n)
        for (i in seq_len(n)) {
            nb <- adj.list[[i]]
            w <- weight.list[[i]]
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
        s.size <- min(n, max(10L, as.integer(sqrt(sample.pairs))))
        s.idx <- sample.int(n, size = s.size, replace = FALSE)

        bfs.dist <- function(src) {
            dist <- rep.int(NA_integer_, n)
            dist[src] <- 0L
            q <- integer(n)
            head <- 1L
            tail <- 1L
            q[tail] <- src

            while (head <= tail) {
                v <- q[head]
                head <- head + 1L
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

        d.samples <- integer(0L)
        per.src <- as.integer(ceiling(sample.pairs / length(s.idx)))
        for (src in s.idx) {
            d <- bfs.dist(src)
            ok <- which(!is.na(d))
            if (length(ok) <= 1L) next
            t.idx <- sample(ok, size = min(length(ok), per.src), replace = TRUE)
            d.samples <- c(d.samples, d[t.idx])
        }

        d.samples <- d.samples[d.samples > 0L]
        if (length(d.samples) == 0L) return(NA_real_)

        if (distance.scale == "diameter") return(as.double(max(d.samples)))

        ## p95_geodesic
        as.double(stats::quantile(d.samples, probs = 0.95, type = 7, names = FALSE))
    }

    scale.val <- estimate.scale()

    ## ------------------------------
    ## Metrics
    ## ------------------------------
    rmse <- function(x, y) sqrt(mean((x - y)^2))

    rel.rmse <- function(x, y, eps) {
        denom <- pmax(y, eps)
        sqrt(mean(((x - y) / denom)^2))
    }

    log.rmse <- function(x, y, eps) {
        lx <- log(pmax(x, 0) + eps)
        ly <- log(pmax(y, 0) + eps)
        sqrt(mean((lx - ly)^2))
    }

    ## ------------------------------
    ## Grid search
    ## ------------------------------
    results <- expand.grid(m = m.grid,
                           eta = eta.grid,
                           kappa = kappa.grid,
                           stringsAsFactors = FALSE)

    results$rmse <- NA_real_
    results$corr <- NA_real_
    results$rel_rmse <- NA_real_
    results$log_rmse <- NA_real_
    results$obj <- NA_real_

    if (verbose) {
        message("Grid size: ", nrow(results), " (m x eta x kappa)")
        if (!is.na(scale.val)) message("Distance scale (", distance.scale, "): ", signif(scale.val, 4))
    }

    for (k in seq_len(nrow(results))) {
        m <- as.integer(results$m[k])
        eta <- as.double(results$eta[k])
        kappa <- as.double(results$kappa[k])

        v <- rho0
        if (m > 0L) {
            for (step in seq_len(m)) {
                v <- apply.P(v, eta)
                v <- normalize.rho(v)
            }
        }

        ## Post-smoothing power transform (Step C)
        if (!isTRUE(all.equal(kappa, 1))) {
            v <- v^kappa
            v <- normalize.rho(v)
        }

        ## Metrics
        results$rmse[k] <- rmse(v, rhoT)
        results$corr[k] <- suppressWarnings(stats::cor(v, rhoT))
        results$rel_rmse[k] <- rel.rmse(v, rhoT, eps = max(eps.rel, min.rho))
        results$log_rmse[k] <- log.rmse(v, rhoT, eps = max(eps.log, min.rho))

        ## Objective (lower is better)
        if (objective == "rmse") {
            results$obj[k] <- results$rmse[k]
        } else if (objective == "corr") {
            results$obj[k] <- -results$corr[k]
        } else if (objective == "rel_rmse") {
            results$obj[k] <- results$rel_rmse[k]
        } else if (objective == "log_rmse") {
            results$obj[k] <- results$log_rmse[k]
        } else if (objective == "hybrid") {
            results$obj[k] <- -(results$corr[k] - lambda * results$rmse[k])
        } else {
            stop("Unknown objective.")
        }
    }

    ## Select best
    best.idx <- which.min(results$obj)
    best.row <- results[best.idx, , drop = FALSE]

    best.list <- list(
        m = as.integer(best.row$m),
        eta = as.double(best.row$eta),
        kappa = as.double(best.row$kappa),
        rmse = as.double(best.row$rmse),
        corr = as.double(best.row$corr),
        rel_rmse = as.double(best.row$rel_rmse),
        log_rmse = as.double(best.row$log_rmse),
        objective = objective,
        lambda = lambda,
        m_over_scale = if (!is.na(scale.val) && scale.val > 0) as.double(best.row$m) / scale.val else NA_real_
    )

    ## Compute rho.best if requested
    rho.best <- NULL
    if (isTRUE(return.best.rho)) {
        v <- rho0
        if (best.list$m > 0L) {
            for (step in seq_len(best.list$m)) {
                v <- apply.P(v, best.list$eta)
                v <- normalize.rho(v)
            }
        }
        if (!isTRUE(all.equal(best.list$kappa, 1))) {
            v <- v^best.list$kappa
            v <- normalize.rho(v)
        }
        rho.best <- v
    }

    ## Sort grid by objective (ascending)
    results <- results[order(results$obj), ]

    list(
        rho0 = rho0,
        rhoT = rhoT,
        scale = list(type = distance.scale, value = scale.val),
        grid = results,
        best = best.list,
        rho.best = rho.best
    )
}

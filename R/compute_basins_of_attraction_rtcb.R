# RTCB basin construction (internal helper used by compute.gfc(modulation = "RTCB"))
#
# Implements relaxed trajectory-constrained basin construction with:
# - Stage A: n_min-limited candidate region from Dijkstra
# - Stage B: constrained multi-label search with Option 1/2
# - Stage C: basin extraction from accepted boundary sinks
compute.basins.of.attraction.rtcb <- function(adj.list,
                                              weight.list,
                                              y,
                                              edge.length.quantile.thld = 0.9,
                                              with.trajectories = FALSE,
                                              n.min = NULL,
                                              m.min = 0.7,
                                              q.min = NULL,
                                              run.max = 2L,
                                              tau0 = 0.1,
                                              kappa = 1.5,
                                              k.max = 24L,
                                              h.max = 1000L,
                                              d.path.max = Inf,
                                              eta.step = 0,
                                              epsilon.d = 1e-12,
                                              eps.num = 1e-12,
                                              sink.prune = FALSE,
                                              sink.prune.max.iter = 5L,
                                              max.paths.per.sink = 64L) {
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }
    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop("y must be a numeric vector with length equal to length(adj.list)")
    }
    if (!is.numeric(edge.length.quantile.thld) ||
        length(edge.length.quantile.thld) != 1L ||
        is.na(edge.length.quantile.thld) ||
        edge.length.quantile.thld <= 0 ||
        edge.length.quantile.thld > 1) {
        stop("edge.length.quantile.thld must be a single numeric value in (0, 1]")
    }

    n <- length(adj.list)
    if (is.null(n.min)) {
        n.min <- max(20L, as.integer(ceiling(sqrt(n))))
    }

    if (!is.null(q.min)) {
        if (!is.numeric(q.min) || length(q.min) != 1L || is.na(q.min) || q.min <= 0 || q.min > 1) {
            stop("q.min must be NULL or a single numeric value in (0, 1]")
        }
    }
    if (!is.numeric(m.min) || length(m.min) != 1L || is.na(m.min) || m.min <= -1 || m.min >= 1) {
        stop("m.min must be a single numeric value in (-1, 1)")
    }

    params <- list(
        n_min = as.integer(max(1L, n.min)),
        m_min = as.numeric(m.min),
        q_min = if (is.null(q.min)) as.numeric(-1) else as.numeric(q.min),
        run_max = as.integer(max(0L, run.max)),
        tau0 = as.numeric(max(0, tau0)),
        kappa = as.numeric(max(1e-9, kappa)),
        k_max = as.integer(max(1L, k.max)),
        h_max = as.integer(max(1L, h.max)),
        d_path_max = as.numeric(d.path.max),
        eta_step = as.numeric(max(0, eta.step)),
        epsilon_d = as.numeric(max(1e-18, epsilon.d)),
        eps_num = as.numeric(max(1e-18, eps.num)),
        sink_prune = as.logical(sink.prune),
        sink_prune_max_iter = as.integer(max(1L, sink.prune.max.iter)),
        max_paths_per_sink = as.integer(max(1L, max.paths.per.sink))
    )

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))

    result <- .Call(
        "S_compute_basins_of_attraction_rtcb",
        adj.list.0based,
        weight.list,
        as.numeric(y),
        as.numeric(edge.length.quantile.thld),
        as.logical(with.trajectories),
        params,
        PACKAGE = "gflow"
    )

    result$n_vertices <- n
    result$y <- y
    result$rtcb.params <- params
    class(result) <- c("basins_of_attraction", "list")
    result
}


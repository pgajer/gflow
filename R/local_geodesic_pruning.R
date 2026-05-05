.normalize.local.prune.method <- function(prune.method, supported = c("none", "local.geodesic")) {
    prune.method <- match.arg(prune.method, supported)
    prune.method
}

.normalize.local.prune.controls <- function(n,
                                            k,
                                            prune.tau,
                                            prune.local.k,
                                            with.pruned.edge.stats) {
    if (is.null(prune.tau)) {
        prune.tau <- 1.05
    }
    if (!is.numeric(prune.tau) || length(prune.tau) != 1L ||
        !is.finite(prune.tau) || prune.tau <= 1) {
        stop("'prune.tau' must be a finite numeric scalar greater than 1.", call. = FALSE)
    }
    if (is.null(prune.local.k)) {
        prune.local.k <- k
    }
    if (!is.numeric(prune.local.k) || length(prune.local.k) != 1L ||
        !is.finite(prune.local.k) || prune.local.k != floor(prune.local.k) ||
        prune.local.k < 1L || prune.local.k >= n) {
        stop("'prune.local.k' must be a positive integer smaller than nrow(X).",
             call. = FALSE)
    }
    if (!is.logical(with.pruned.edge.stats) || length(with.pruned.edge.stats) != 1L ||
        is.na(with.pruned.edge.stats)) {
        stop("'with.pruned.edge.stats' must be TRUE or FALSE.", call. = FALSE)
    }
    list(
        prune.tau = as.numeric(prune.tau),
        prune.local.k = as.integer(prune.local.k),
        with.pruned.edge.stats = isTRUE(with.pruned.edge.stats)
    )
}

.empty.pruned.edge.stats <- function() {
    data.frame(
        u = integer(),
        v = integer(),
        edge_length = numeric(),
        alt_path_length = numeric(),
        path_edge_ratio = numeric()
    )
}

.exact.knn.index <- function(X, k) {
    n <- nrow(X)
    out <- matrix(NA_integer_, nrow = n, ncol = k)
    for (i in seq_len(n)) {
        d <- rowSums((t(t(X) - X[i, ]))^2)
        d[[i]] <- Inf
        out[i, ] <- order(d, seq_len(n))[seq_len(k)]
    }
    out
}

.local.dijkstra.distance <- function(adj.list,
                                     weight.list,
                                     local.vertices,
                                     source,
                                     target,
                                     excluded.u,
                                     excluded.v,
                                     cutoff) {
    n <- length(adj.list)
    in.local <- rep(FALSE, n)
    in.local[local.vertices] <- TRUE
    in.local[[source]] <- TRUE
    in.local[[target]] <- TRUE
    dist <- rep(Inf, n)
    visited <- rep(FALSE, n)
    dist[[source]] <- 0
    repeat {
        candidates <- which(in.local & !visited & is.finite(dist))
        if (!length(candidates)) {
            return(Inf)
        }
        u <- candidates[[which.min(dist[candidates])]]
        du <- dist[[u]]
        if (du > cutoff) {
            return(Inf)
        }
        if (u == target) {
            return(du)
        }
        visited[[u]] <- TRUE
        nbrs <- adj.list[[u]]
        weights <- weight.list[[u]]
        if (!length(nbrs)) {
            next
        }
        for (pos in seq_along(nbrs)) {
            v <- nbrs[[pos]]
            if (!in.local[[v]] || visited[[v]]) {
                next
            }
            if ((u == excluded.u && v == excluded.v) ||
                (u == excluded.v && v == excluded.u)) {
                next
            }
            nd <- du + weights[[pos]]
            if (nd < dist[[v]] && nd <= cutoff) {
                dist[[v]] <- nd
            }
        }
    }
}

.remove.undirected.edge <- function(adj.list, weight.list, u, v) {
    remove.one <- function(a, w, from, to) {
        pos <- which(a[[from]] == to)
        if (length(pos)) {
            keep <- seq_along(a[[from]])[-pos[[1]]]
            a[[from]] <- as.integer(a[[from]][keep])
            w[[from]] <- as.numeric(w[[from]][keep])
        }
        list(adj = a, weight = w)
    }
    first <- remove.one(adj.list, weight.list, u, v)
    second <- remove.one(first$adj, first$weight, v, u)
    list(adj_list = second$adj, weight_list = second$weight)
}

.prune.graph.local.geodesic <- function(X,
                                        adj.list,
                                        weight.list,
                                        k,
                                        prune.tau = 1.05,
                                        prune.local.k = NULL,
                                        with.pruned.edge.stats = FALSE) {
    n <- nrow(X)
    controls <- .normalize.local.prune.controls(
        n, k, prune.tau, prune.local.k, with.pruned.edge.stats
    )
    edges <- .graph.edge.table(adj.list, weight.list)
    n.edges.before <- nrow(edges)
    if (!n.edges.before) {
        return(list(
            adj_list = adj.list,
            weight_list = weight.list,
            n_edges_before_pruning = 0L,
            n_edges_after_pruning = 0L,
            n_pruned_edges = 0L,
            pruned_edge_stats = .empty.pruned.edge.stats(),
            prune_tau = controls$prune.tau,
            prune_local_k = controls$prune.local.k,
            with_pruned_edge_stats = controls$with.pruned.edge.stats
        ))
    }

    local.index <- .exact.knn.index(X, controls$prune.local.k)
    edges <- edges[order(-edges$weight, edges$from, edges$to), , drop = FALSE]
    stats <- vector("list", n.edges.before)
    stats.cursor <- 0L
    tol <- 1e-12
    for (r in seq_len(nrow(edges))) {
        u <- edges$from[[r]]
        v <- edges$to[[r]]
        current <- .graph.edge.table(adj.list, weight.list)
        key <- current$from == u & current$to == v
        if (!any(key)) {
            next
        }
        edge.length <- current$weight[which(key)[[1]]]
        local.vertices <- sort(unique(c(u, v, local.index[u, ], local.index[v, ])))
        cutoff <- controls$prune.tau * edge.length
        alt <- .local.dijkstra.distance(
            adj.list, weight.list, local.vertices, u, v, u, v, cutoff + tol
        )
        if (is.finite(alt) && alt <= cutoff + tol) {
            removed <- .remove.undirected.edge(adj.list, weight.list, u, v)
            adj.list <- removed$adj_list
            weight.list <- removed$weight_list
            if (controls$with.pruned.edge.stats) {
                stats.cursor <- stats.cursor + 1L
                stats[[stats.cursor]] <- data.frame(
                    u = u,
                    v = v,
                    edge_length = edge.length,
                    alt_path_length = alt,
                    path_edge_ratio = alt / edge.length
                )
            }
        }
    }
    stats.df <- if (controls$with.pruned.edge.stats && stats.cursor > 0L) {
        do.call(rbind, stats[seq_len(stats.cursor)])
    } else {
        .empty.pruned.edge.stats()
    }
    n.edges.after <- nrow(.graph.edge.table(adj.list, weight.list))
    list(
        adj_list = adj.list,
        weight_list = weight.list,
        n_edges_before_pruning = n.edges.before,
        n_edges_after_pruning = n.edges.after,
        n_pruned_edges = n.edges.before - n.edges.after,
        pruned_edge_stats = stats.df,
        prune_tau = controls$prune.tau,
        prune_local_k = controls$prune.local.k,
        with_pruned_edge_stats = controls$with.pruned.edge.stats
    )
}

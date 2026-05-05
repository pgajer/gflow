.graph.edge.table <- function(adj.list, weight.list) {
    n <- length(adj.list)
    rows <- vector("list", sum(lengths(adj.list)))
    cursor <- 0L
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        weights <- weight.list[[i]]
        if (!length(nbrs)) {
            next
        }
        for (pos in seq_along(nbrs)) {
            j <- as.integer(nbrs[[pos]])
            if (i < j) {
                cursor <- cursor + 1L
                rows[[cursor]] <- data.frame(
                    from = i,
                    to = j,
                    weight = as.numeric(weights[[pos]])
                )
            }
        }
    }
    if (cursor == 0L) {
        return(data.frame(from = integer(), to = integer(), weight = numeric()))
    }
    do.call(rbind, rows[seq_len(cursor)])
}

.graph.from.edge.table <- function(n, edges) {
    adj <- vector("list", n)
    weights <- vector("list", n)
    if (is.null(edges) || !nrow(edges)) {
        for (i in seq_len(n)) {
            adj[[i]] <- integer(0)
            weights[[i]] <- numeric(0)
        }
        return(list(adj_list = adj, weight_list = weights))
    }
    edges <- edges[order(edges$from, edges$to), , drop = FALSE]
    for (r in seq_len(nrow(edges))) {
        u <- as.integer(edges$from[[r]])
        v <- as.integer(edges$to[[r]])
        w <- as.numeric(edges$weight[[r]])
        adj[[u]] <- c(adj[[u]], v)
        weights[[u]] <- c(weights[[u]], w)
        adj[[v]] <- c(adj[[v]], u)
        weights[[v]] <- c(weights[[v]], w)
    }
    for (i in seq_len(n)) {
        if (length(adj[[i]]) > 1L) {
            ord <- order(adj[[i]])
            adj[[i]] <- as.integer(adj[[i]][ord])
            weights[[i]] <- as.numeric(weights[[i]][ord])
        } else {
            adj[[i]] <- as.integer(adj[[i]])
            weights[[i]] <- as.numeric(weights[[i]])
        }
    }
    list(adj_list = adj, weight_list = weights)
}

.graph.components <- function(adj.list) {
    n <- length(adj.list)
    comp <- rep.int(NA_integer_, n)
    comp.id <- 0L
    for (start in seq_len(n)) {
        if (!is.na(comp[[start]])) {
            next
        }
        comp.id <- comp.id + 1L
        comp[[start]] <- comp.id
        stack <- start
        while (length(stack)) {
            u <- stack[[length(stack)]]
            stack <- stack[-length(stack)]
            for (v in adj.list[[u]]) {
                if (is.na(comp[[v]])) {
                    comp[[v]] <- comp.id
                    stack <- c(stack, v)
                }
            }
        }
    }
    list(component_id = comp, n_components = comp.id)
}

.euclidean.distance <- function(X, i, j) {
    sqrt(sum((X[i, , drop = TRUE] - X[j, , drop = TRUE])^2))
}

.component.kruskal <- function(candidates, n.components) {
    if (n.components <= 1L) {
        return(data.frame(from = integer(), to = integer(), weight = numeric()))
    }
    if (is.null(candidates) || !nrow(candidates)) {
        return(NULL)
    }
    candidates <- candidates[order(candidates$weight, candidates$from, candidates$to), , drop = FALSE]
    parent <- seq_len(n.components)
    rank <- integer(n.components)
    find <- function(x) {
        while (parent[[x]] != x) {
            parent[[x]] <<- parent[[parent[[x]]]]
            x <- parent[[x]]
        }
        x
    }
    unite <- function(a, b) {
        ra <- find(a)
        rb <- find(b)
        if (ra == rb) {
            return(FALSE)
        }
        if (rank[[ra]] < rank[[rb]]) {
            tmp <- ra
            ra <- rb
            rb <- tmp
        }
        parent[[rb]] <<- ra
        if (rank[[ra]] == rank[[rb]]) {
            rank[[ra]] <<- rank[[ra]] + 1L
        }
        TRUE
    }
    keep <- vector("list", n.components - 1L)
    cursor <- 0L
    for (r in seq_len(nrow(candidates))) {
        if (unite(candidates$component_from[[r]], candidates$component_to[[r]])) {
            cursor <- cursor + 1L
            keep[[cursor]] <- candidates[r, c("from", "to", "weight"), drop = FALSE]
            if (cursor == n.components - 1L) {
                break
            }
        }
    }
    if (cursor != n.components - 1L) {
        return(NULL)
    }
    do.call(rbind, keep[seq_len(cursor)])
}

.best.component.bridges <- function(X, comp, bridge.index = NULL, bridge.k = NULL) {
    n <- nrow(X)
    n.components <- max(comp)
    best <- new.env(parent = emptyenv())
    update.best <- function(i, j) {
        ci <- comp[[i]]
        cj <- comp[[j]]
        if (ci == cj) {
            return()
        }
        ca <- min(ci, cj)
        cb <- max(ci, cj)
        u <- min(i, j)
        v <- max(i, j)
        key <- paste(ca, cb, sep = "-")
        d <- .euclidean.distance(X, u, v)
        old <- best[[key]]
        if (is.null(old) || d < old$weight ||
            (d == old$weight && (u < old$from || (u == old$from && v < old$to)))) {
            best[[key]] <- list(
                component_from = ca,
                component_to = cb,
                from = u,
                to = v,
                weight = d
            )
        }
    }
    if (is.null(bridge.index)) {
        for (i in seq_len(n - 1L)) {
            for (j in (i + 1L):n) {
                update.best(i, j)
            }
        }
    } else {
        for (i in seq_len(n)) {
            for (j in bridge.index[i, seq_len(bridge.k)]) {
                update.best(i, as.integer(j))
            }
        }
    }
    keys <- ls(best, all.names = TRUE)
    if (!length(keys)) {
        return(data.frame(
            component_from = integer(), component_to = integer(),
            from = integer(), to = integer(), weight = numeric()
        ))
    }
    rows <- lapply(keys, function(key) as.data.frame(best[[key]]))
    out <- do.call(rbind, rows)
    out[order(out$component_from, out$component_to), , drop = FALSE]
}

.full.euclidean.mst.edges <- function(X) {
    n <- nrow(X)
    candidates <- vector("list", n * (n - 1L) / 2L)
    cursor <- 0L
    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            cursor <- cursor + 1L
            candidates[[cursor]] <- data.frame(
                component_from = i,
                component_to = j,
                from = i,
                to = j,
                weight = .euclidean.distance(X, i, j)
            )
        }
    }
    .component.kruskal(do.call(rbind, candidates), n)
}

.normalize.bridge.controls <- function(n, k, bridge.k, bridge.k.max, bridge.growth) {
    if (!is.numeric(bridge.growth) || length(bridge.growth) != 1L ||
        !is.finite(bridge.growth) || bridge.growth <= 1) {
        stop("'bridge.growth' must be a finite numeric scalar greater than 1.", call. = FALSE)
    }
    if (is.null(bridge.k)) {
        bridge.k <- min(n - 1L, max(as.integer(k), 10L))
    }
    if (!is.numeric(bridge.k) || length(bridge.k) != 1L || !is.finite(bridge.k) ||
        bridge.k != floor(bridge.k) || bridge.k < 1L || bridge.k >= n) {
        stop("'bridge.k' must be a positive integer smaller than nrow(X).", call. = FALSE)
    }
    if (is.null(bridge.k.max)) {
        bridge.k.max <- min(n - 1L, max(50L, 8L * as.integer(k), as.integer(bridge.k)))
    }
    if (!is.numeric(bridge.k.max) || length(bridge.k.max) != 1L ||
        !is.finite(bridge.k.max) || bridge.k.max != floor(bridge.k.max) ||
        bridge.k.max < bridge.k || bridge.k.max >= n) {
        stop("'bridge.k.max' must be an integer between bridge.k and nrow(X) - 1.",
             call. = FALSE)
    }
    list(
        bridge.k = as.integer(bridge.k),
        bridge.k.max = as.integer(bridge.k.max),
        bridge.growth = as.numeric(bridge.growth)
    )
}

.augment.graph.with.component.mst <- function(X,
                                              adj.list,
                                              weight.list,
                                              k,
                                              connect.components = FALSE,
                                              connect.method = c("component.mst", "component.mst.ann", "global.mst"),
                                              bridge.k = NULL,
                                              bridge.k.max = NULL,
                                              bridge.growth = 2) {
    connect.method <- match.arg(connect.method)
    n <- nrow(X)
    controls <- .normalize.bridge.controls(n, k, bridge.k, bridge.k.max, bridge.growth)
    edges <- .graph.edge.table(adj.list, weight.list)
    comps.before <- .graph.components(adj.list)
    added <- data.frame(from = integer(), to = integer(), weight = numeric())
    bridge.method <- "none"
    bridge.k.used <- NA_integer_
    exact.fallback <- FALSE

    if (isTRUE(connect.components) && comps.before$n_components > 1L) {
        if (identical(connect.method, "component.mst")) {
            candidates <- .best.component.bridges(X, comps.before$component_id)
            added <- .component.kruskal(candidates, comps.before$n_components)
            bridge.method <- "exact"
        } else if (identical(connect.method, "global.mst")) {
            mst <- .full.euclidean.mst.edges(X)
            existing <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "-")
            mst.key <- paste(mst$from, mst$to, sep = "-")
            added <- mst[!(mst.key %in% existing), c("from", "to", "weight"), drop = FALSE]
            bridge.method <- "global.mst"
        } else {
            raw.bridge <- .Call(
                "S_kNN",
                X,
                as.integer(controls$bridge.k.max + 1L),
                PACKAGE = "gflow"
            )$indices
            bridge.index <- matrix(NA_integer_, nrow = n, ncol = controls$bridge.k.max)
            for (i in seq_len(n)) {
                row <- raw.bridge[i, ]
                non.self <- row[row != (i - 1L)]
                if (length(non.self) < controls$bridge.k.max) {
                    stop("ANN bridge search returned too few non-self neighbors.", call. = FALSE)
                }
                bridge.index[i, ] <- non.self[seq_len(controls$bridge.k.max)] + 1L
            }
            current.k <- controls$bridge.k
            repeat {
                candidates <- .best.component.bridges(
                    X, comps.before$component_id,
                    bridge.index = bridge.index,
                    bridge.k = current.k
                )
                proposed <- .component.kruskal(candidates, comps.before$n_components)
                if (!is.null(proposed)) {
                    added <- proposed
                    bridge.method <- "ann"
                    bridge.k.used <- current.k
                    break
                }
                if (current.k >= controls$bridge.k.max) {
                    candidates <- .best.component.bridges(X, comps.before$component_id)
                    added <- .component.kruskal(candidates, comps.before$n_components)
                    bridge.method <- "ann_then_exact"
                    exact.fallback <- TRUE
                    bridge.k.used <- current.k
                    break
                }
                current.k <- min(
                    controls$bridge.k.max,
                    max(current.k + 1L, as.integer(ceiling(current.k * controls$bridge.growth)))
                )
            }
        }
        if (is.null(added)) {
            stop("Internal error: component MST bridge construction failed.", call. = FALSE)
        }
        if (nrow(added)) {
            edges <- rbind(edges, added)
            edges <- edges[!duplicated(paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "-")),
                           , drop = FALSE]
        }
    }

    graph <- .graph.from.edge.table(n, edges)
    comps.after <- .graph.components(graph$adj_list)
    list(
        adj_list = graph$adj_list,
        weight_list = graph$weight_list,
        n_components_before = comps.before$n_components,
        n_components_after = comps.after$n_components,
        component_id_before = comps.before$component_id,
        component_id_after = comps.after$component_id,
        mst_edge_matrix = as.matrix(added[, c("from", "to"), drop = FALSE]),
        mst_edge_weight = as.numeric(added$weight),
        n_mst_edges_added = nrow(added),
        connect_components = isTRUE(connect.components),
        connect_method = connect.method,
        bridge_method = bridge.method,
        bridge_k = controls$bridge.k,
        bridge_k_max = controls$bridge.k.max,
        bridge_growth = controls$bridge.growth,
        bridge_k_used = bridge.k.used,
        bridge_exact_fallback_used = exact.fallback
    )
}

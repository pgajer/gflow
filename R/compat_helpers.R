## Compatibility/internal helpers used by optional workflows.
## These definitions keep package checks clean and provide sensible fallbacks.

.postprocess.gfc.flow.v2 <- function(result) {
    .postprocess.gfc.flow(result)
}

.vertex.cell.gfc.flow.v1 <- function(gfc.flow, vertex) {
    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0L) {
        return(data.frame(
            vertex = integer(0),
            min.basin.idx = integer(0),
            max.basin.idx = integer(0),
            min.label = character(0),
            max.label = character(0),
            min.is.spurious = logical(0),
            max.is.spurious = logical(0),
            stringsAsFactors = FALSE
        ))
    }

    summary.df <- gfc.flow$summary.all
    if (is.null(summary.df) || nrow(summary.df) == 0L) {
        summary.df <- gfc.flow$summary
    }

    min.rows <- if (!is.null(summary.df)) summary.df[summary.df$type == "min", , drop = FALSE] else NULL
    max.rows <- if (!is.null(summary.df)) summary.df[summary.df$type == "max", , drop = FALSE] else NULL
    min.vertex.to.idx <- if (!is.null(min.rows)) setNames(seq_len(nrow(min.rows)), as.character(min.rows$vertex)) else integer(0)
    max.vertex.to.idx <- if (!is.null(max.rows)) setNames(seq_len(nrow(max.rows)), as.character(max.rows$vertex)) else integer(0)
    min.vertex.to.label <- if (!is.null(min.rows)) setNames(as.character(min.rows$label), as.character(min.rows$vertex)) else character(0)
    max.vertex.to.label <- if (!is.null(max.rows)) setNames(as.character(max.rows$label), as.character(max.rows$vertex)) else character(0)
    min.vertex.to.spurious <- if (!is.null(min.rows) && "is.spurious" %in% names(min.rows)) setNames(as.logical(min.rows$is.spurious), as.character(min.rows$vertex)) else logical(0)
    max.vertex.to.spurious <- if (!is.null(max.rows) && "is.spurious" %in% names(max.rows)) setNames(as.logical(max.rows$is.spurious), as.character(max.rows$vertex)) else logical(0)

    hits <- list()
    for (tr in gfc.flow$trajectories) {
        v <- as.integer(tr$vertices)
        if (length(v) == 0L || !(vertex %in% v)) {
            next
        }
        min.v <- as.integer(tr$start.vertex)
        max.v <- as.integer(tr$end.vertex)
        hits[[length(hits) + 1L]] <- data.frame(
            vertex = as.integer(vertex),
            min.basin.idx = if (as.character(min.v) %in% names(min.vertex.to.idx)) as.integer(min.vertex.to.idx[[as.character(min.v)]]) else NA_integer_,
            max.basin.idx = if (as.character(max.v) %in% names(max.vertex.to.idx)) as.integer(max.vertex.to.idx[[as.character(max.v)]]) else NA_integer_,
            min.label = if (as.character(min.v) %in% names(min.vertex.to.label)) min.vertex.to.label[[as.character(min.v)]] else paste0("v", min.v),
            max.label = if (as.character(max.v) %in% names(max.vertex.to.label)) max.vertex.to.label[[as.character(max.v)]] else paste0("v", max.v),
            min.is.spurious = if (as.character(min.v) %in% names(min.vertex.to.spurious)) isTRUE(min.vertex.to.spurious[[as.character(min.v)]]) else FALSE,
            max.is.spurious = if (as.character(max.v) %in% names(max.vertex.to.spurious)) isTRUE(max.vertex.to.spurious[[as.character(max.v)]]) else FALSE,
            stringsAsFactors = FALSE
        )
    }

    if (length(hits) == 0L) {
        return(data.frame(
            vertex = integer(0),
            min.basin.idx = integer(0),
            max.basin.idx = integer(0),
            min.label = character(0),
            max.label = character(0),
            min.is.spurious = logical(0),
            max.is.spurious = logical(0),
            stringsAsFactors = FALSE
        ))
    }

    unique(do.call(rbind, hits))
}

apply.spectral.smoothing <- function(y, V, filter.weights) {
    y <- as.numeric(y)
    V <- as.matrix(V)
    filter.weights <- as.numeric(filter.weights)
    if (nrow(V) != length(y)) {
        stop("nrow(V) must equal length(y).")
    }
    if (ncol(V) != length(filter.weights)) {
        stop("length(filter.weights) must equal ncol(V).")
    }
    coeff <- crossprod(V, y)
    as.numeric(V %*% (filter.weights * coeff))
}

generate.eta.grid <- function(eigenvalues, filter.type, n.candidates = 40L) {
    n.candidates <- max(5L, as.integer(n.candidates))
    ev <- sort(as.numeric(eigenvalues))
    ev <- ev[is.finite(ev) & ev > .Machine$double.eps]
    if (length(ev) == 0L) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }

    lo <- 1 / max(ev)
    hi <- 1 / min(ev)
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= lo) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }

    if (identical(filter.type, "butterworth")) {
        return(exp(seq(log(max(lo, 1e-6)), log(hi), length.out = n.candidates)))
    }
    exp(seq(log(max(lo, 1e-6)), log(hi), length.out = n.candidates))
}

prepare.basin.da.inputs <- function(X.rel, y, basin.vertices, counts,
                                    min.prevalence = 0.1, min.count = 10L) {
    X.rel <- as.matrix(X.rel)
    counts <- as.matrix(counts)
    y <- as.integer(y)
    idx <- as.integer(intersect(as.integer(basin.vertices), seq_len(nrow(X.rel))))
    if (length(idx) < 2L) {
        stop("Need at least 2 basin vertices after mapping.")
    }

    y.C <- y[idx]
    if (ncol(counts) == nrow(X.rel)) {
        counts.C <- counts[, idx, drop = FALSE]
    } else if (nrow(counts) == nrow(X.rel)) {
        counts.C <- t(counts[idx, , drop = FALSE])
    } else {
        stop("counts dimensions are not compatible with X.rel.")
    }

    keep <- rowMeans(counts.C > as.integer(min.count), na.rm = TRUE) >= as.numeric(min.prevalence)
    counts.C <- counts.C[keep, , drop = FALSE]
    list(counts.C = counts.C, y.C = y.C)
}

run.da.edger <- function(counts.C, y.C) {
    counts.C <- as.matrix(counts.C)
    y.C <- as.integer(y.C)
    grp1 <- which(y.C == 1L)
    grp0 <- which(y.C == 0L)
    if (length(grp1) < 2L || length(grp0) < 2L) {
        stop("Need at least 2 samples per class for DA fallback.")
    }

    n.features <- nrow(counts.C)
    pvals <- rep(NA_real_, n.features)
    logfc <- rep(NA_real_, n.features)

    for (i in seq_len(n.features)) {
        a <- log1p(counts.C[i, grp1])
        b <- log1p(counts.C[i, grp0])
        logfc[i] <- mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)
        pvals[i] <- tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
    }

    out <- data.frame(
        logFC = logfc,
        PValue = pvals,
        FDR = stats::p.adjust(pvals, method = "BH"),
        stringsAsFactors = FALSE
    )
    rownames(out) <- rownames(counts.C)
    out
}

find.path.in.subgraph <- function(adj.list, start, end, max.path.length = 10L) {
    start <- as.integer(start)
    end <- as.integer(end)
    max.path.length <- as.integer(max.path.length)
    n <- length(adj.list)
    if (start < 1L || start > n || end < 1L || end > n) return(NULL)
    if (start == end) return(start)

    q <- list(start)
    visited <- rep(FALSE, n)
    visited[start] <- TRUE
    prev <- rep(NA_integer_, n)
    depth <- rep(Inf, n)
    depth[start] <- 0L

    while (length(q) > 0L) {
        v <- q[[1L]]
        q <- q[-1L]
        if (depth[v] >= max.path.length) next
        for (u in as.integer(adj.list[[v]])) {
            if (u < 1L || u > n || visited[u]) next
            visited[u] <- TRUE
            prev[u] <- v
            depth[u] <- depth[v] + 1L
            if (u == end) {
                path <- end
                cur <- end
                while (!is.na(prev[cur])) {
                    cur <- prev[cur]
                    path <- c(cur, path)
                }
                return(path)
            }
            q[[length(q) + 1L]] <- u
        }
    }
    NULL
}

summarize.clinical <- function(df) {
    if (is.null(df) || nrow(df) == 0L) return(list())
    out <- lapply(df, function(col) {
        if (is.numeric(col)) {
            c(mean = mean(col, na.rm = TRUE), sd = stats::sd(col, na.rm = TRUE))
        } else {
            sort(table(col, useNA = "ifany"), decreasing = TRUE)
        }
    })
    out
}

compute.subgraph.density <- function(adj.list) {
    n <- length(adj.list)
    if (n <= 1L) return(0)
    edges <- 0L
    for (i in seq_len(n)) {
        nbr <- as.integer(adj.list[[i]])
        edges <- edges + sum(nbr > i & nbr <= n, na.rm = TRUE)
    }
    edges / (n * (n - 1) / 2)
}

compute.diameter <- function(adj.list) {
    n <- length(adj.list)
    if (n <= 1L) return(0L)
    bfs_dist <- function(src) {
        d <- rep(Inf, n)
        d[src] <- 0L
        q <- src
        head <- 1L
        while (head <= length(q)) {
            v <- q[head]
            head <- head + 1L
            for (u in as.integer(adj.list[[v]])) {
                if (u >= 1L && u <= n && !is.finite(d[u])) {
                    d[u] <- d[v] + 1L
                    q <- c(q, u)
                }
            }
        }
        d
    }
    diam <- 0L
    for (i in seq_len(n)) {
        di <- bfs_dist(i)
        finite <- di[is.finite(di)]
        if (length(finite) > 0L) diam <- max(diam, max(finite))
    }
    as.integer(diam)
}

determine.optimal.k <- function(similarity.matrix) {
    n <- nrow(similarity.matrix)
    if (is.null(n) || n < 2L) return(1L)
    max(2L, min(n - 1L, as.integer(round(sqrt(n)))))
}

utils::globalVariables(c("k.val", "n.abundance", "n.total", "graph.3d", "p", "feature"))

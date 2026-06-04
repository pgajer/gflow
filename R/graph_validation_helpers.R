## Private graph/response validators shared by graph-domain methods.

.validate.positive.integer.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x != floor(x)) {
        stop(sprintf("%s must be a positive integer scalar.", name))
    }
    as.integer(x)
}

.validate.metric.graph.lowpass.response <- function(y, n, name) {
    if (!is.numeric(y) || length(y) != n) {
        stop(sprintf("%s must be a numeric vector of length %d.", name, n))
    }
    y <- as.double(y)
    if (any(!is.finite(y))) stop(sprintf("%s cannot contain NA/NaN/Inf.", name))
    y
}

.validate.metric.graph.lowpass.graph <- function(adj.list, weight.list) {
    if (!is.list(adj.list)) stop("adj.list must be a list.")
    if (!is.list(weight.list)) stop("weight.list must be a list.")
    n <- length(adj.list)
    if (n < 2L) stop("adj.list must contain at least two vertices.")
    if (length(weight.list) != n) {
        stop("adj.list and weight.list must have the same length.")
    }

    adj.norm <- vector("list", n)
    weight.norm <- vector("list", n)
    tol <- 1e-10
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        wts <- weight.list[[i]]
        if (!is.numeric(nbrs) && !is.integer(nbrs)) {
            stop(sprintf("adj.list[[%d]] must be numeric/integer.", i))
        }
        if (!is.numeric(wts)) {
            stop(sprintf("weight.list[[%d]] must be numeric.", i))
        }
        nbrs <- as.integer(nbrs)
        wts <- as.double(wts)
        if (length(nbrs) != length(wts)) {
            stop(sprintf(
                "Length mismatch at vertex %d: adj.list=%d, weight.list=%d.",
                i, length(nbrs), length(wts)
            ))
        }
        if (anyNA(nbrs)) stop(sprintf("adj.list[[%d]] contains NA.", i))
        if (any(nbrs < 1L | nbrs > n)) {
            stop(sprintf("adj.list[[%d]] has indices outside 1..n.", i))
        }
        if (any(nbrs == i)) stop(sprintf("adj.list[[%d]] contains self-loops.", i))
        if (anyDuplicated(nbrs)) {
            stop(sprintf("adj.list[[%d]] contains duplicate neighbors.", i))
        }
        if (any(!is.finite(wts))) {
            stop(sprintf("weight.list[[%d]] contains non-finite values.", i))
        }
        if (any(wts <= 0)) {
            stop(sprintf("weight.list[[%d]] contains non-positive values.", i))
        }
        adj.norm[[i]] <- nbrs
        weight.norm[[i]] <- wts
    }
    for (i in seq_len(n)) {
        for (idx in seq_along(adj.norm[[i]])) {
            j <- adj.norm[[i]][idx]
            rev.idx <- match(i, adj.norm[[j]])
            if (is.na(rev.idx)) {
                stop(sprintf(
                    "Graph must be undirected: edge %d -> %d has no reciprocal entry.",
                    i, j
                ))
            }
            w.ij <- weight.norm[[i]][idx]
            w.ji <- weight.norm[[j]][rev.idx]
            if (abs(w.ij - w.ji) > tol * max(1, abs(w.ij), abs(w.ji))) {
                stop(sprintf(
                    "Reciprocal edge weights mismatch for (%d, %d): %.12g vs %.12g.",
                    i, j, w.ij, w.ji
                ))
            }
        }
    }
    list(
        adj.list = adj.norm,
        weight.list = weight.norm,
        adj.list.0based = lapply(adj.norm, function(v) as.integer(v - 1L)),
        weight.list.cpp = lapply(weight.norm, as.double)
    )
}

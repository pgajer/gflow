# Shared helpers for MST-completion benchmark reports.

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

html.escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
}

html.table <- function(x, digits = 4) {
    if (is.null(x) || !nrow(x)) return("<p><em>No rows.</em></p>")
    y <- x
    for (j in seq_along(y)) {
        if (is.numeric(y[[j]])) {
            y[[j]] <- formatC(y[[j]], digits = digits, format = "fg", flag = "#")
        }
    }
    header <- paste(sprintf("<th>%s</th>", html.escape(names(y))), collapse = "")
    rows <- apply(y, 1, function(row) {
        paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(row)), collapse = ""), "</tr>")
    })
    paste0("<table><thead><tr>", header, "</tr></thead><tbody>",
           paste(rows, collapse = "\n"), "</tbody></table>")
}

save.png <- function(path, expr, width = 1400, height = 1000, res = 150) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(path, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr)
    invisible(path)
}

circle.dist.matrix <- function(theta, radius = 1) {
    d <- abs(outer(theta, theta, "-"))
    radius * pmin(d, 2 * pi - d)
}

euclidean.dist.matrix <- function(X) {
    as.matrix(stats::dist(X))
}

weighted.circle.graph <- function(theta, radius = 1) {
    theta <- as.double(theta) %% (2 * pi)
    n <- length(theta)
    ord <- order(theta)
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    add.edge <- function(i, j, w) {
        adj.list[[i]] <<- c(adj.list[[i]], as.integer(j))
        weight.list[[i]] <<- c(weight.list[[i]], as.double(w))
        adj.list[[j]] <<- c(adj.list[[j]], as.integer(i))
        weight.list[[j]] <<- c(weight.list[[j]], as.double(w))
    }
    for (pos in seq_len(n)) {
        i <- ord[pos]
        j <- ord[if (pos == n) 1L else pos + 1L]
        gap <- if (pos == n) theta[ord[1L]] + 2 * pi - theta[i] else theta[j] - theta[i]
        add.edge(i, j, radius * gap)
    }
    list(adj_list = adj.list, weight_list = weight.list)
}

weighted.angular.path.graph <- function(theta, radius = 1) {
    theta <- as.double(theta) %% (2 * pi)
    n <- length(theta)
    ord <- order(theta)
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    add.edge <- function(i, j, w) {
        adj.list[[i]] <<- c(adj.list[[i]], as.integer(j))
        weight.list[[i]] <<- c(weight.list[[i]], as.double(w))
        adj.list[[j]] <<- c(adj.list[[j]], as.integer(i))
        weight.list[[j]] <<- c(weight.list[[j]], as.double(w))
    }
    for (pos in seq_len(n - 1L)) {
        i <- ord[pos]
        j <- ord[pos + 1L]
        add.edge(i, j, radius * (theta[j] - theta[i]))
    }
    list(adj_list = adj.list, weight_list = weight.list)
}

edge.table <- function(adj.list, weight.list = NULL) {
    n <- length(adj.list)
    out <- vector("list", n)
    ptr <- 0L
    for (i in seq_len(n)) {
        nbrs <- as.integer(adj.list[[i]])
        if (!length(nbrs)) next
        wts <- if (is.null(weight.list)) rep(1, length(nbrs)) else as.double(weight.list[[i]])
        keep <- nbrs > i
        if (!any(keep)) next
        ptr <- ptr + 1L
        out[[ptr]] <- data.frame(from = i, to = nbrs[keep], weight = wts[keep])
    }
    if (!ptr) return(data.frame(from = integer(), to = integer(), weight = numeric()))
    do.call(rbind, out[seq_len(ptr)])
}

graph.from.edges <- function(n, edges) {
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)
    if (nrow(edges)) {
        for (r in seq_len(nrow(edges))) {
            i <- as.integer(edges$from[[r]])
            j <- as.integer(edges$to[[r]])
            w <- as.double(edges$weight[[r]])
            if (!is.finite(w) || w < 0) next
            adj.list[[i]] <- c(adj.list[[i]], j)
            weight.list[[i]] <- c(weight.list[[i]], w)
            adj.list[[j]] <- c(adj.list[[j]], i)
            weight.list[[j]] <- c(weight.list[[j]], w)
        }
    }
    list(adj_list = adj.list, weight_list = weight.list)
}

edge.key <- function(i, j) {
    paste(pmin(i, j), pmax(i, j), sep = "-")
}

union.graphs <- function(...) {
    graphs <- list(...)
    n <- length(graphs[[1]]$adj_list)
    edge.map <- new.env(parent = emptyenv())
    for (g in graphs) {
        edges <- edge.table(g$adj_list, g$weight_list)
        if (!nrow(edges)) next
        for (r in seq_len(nrow(edges))) {
            key <- edge.key(edges$from[[r]], edges$to[[r]])
            w <- edges$weight[[r]]
            if (!is.finite(w) || w < 0) next
            if (!exists(key, edge.map, inherits = FALSE) || get(key, edge.map) > w) {
                assign(key, w, edge.map)
            }
        }
    }
    keys <- ls(edge.map)
    if (!length(keys)) return(graph.from.edges(n, data.frame(from = integer(), to = integer(), weight = numeric())))
    ij <- do.call(rbind, strsplit(keys, "-", fixed = TRUE))
    edges <- data.frame(
        from = as.integer(ij[, 1]),
        to = as.integer(ij[, 2]),
        weight = vapply(keys, function(k) get(k, edge.map), numeric(1))
    )
    graph.from.edges(n, edges)
}

igraph.from.graph <- function(adj.list, weight.list = NULL) {
    edges <- edge.table(adj.list, weight.list)
    if (nrow(edges)) {
        edges <- edges[is.finite(edges$weight) & edges$weight >= 0, , drop = FALSE]
    }
    g <- igraph::make_empty_graph(n = length(adj.list), directed = FALSE)
    if (nrow(edges)) {
        g <- igraph::add_edges(g, as.vector(t(as.matrix(edges[, c("from", "to")]))))
        igraph::E(g)$weight <- edges$weight
    }
    g
}

graph.distances <- function(adj.list, weight.list) {
    g <- igraph.from.graph(adj.list, weight.list)
    as.matrix(igraph::distances(g, weights = igraph::E(g)$weight))
}

component.count <- function(adj.list) {
    igraph::components(igraph.from.graph(adj.list))$no
}

cycle.rank <- function(adj.list) {
    n.edges <- nrow(edge.table(adj.list))
    n.edges - length(adj.list) + component.count(adj.list)
}

mst.graph.from.cmst <- function(X) {
    g <- create.cmst.graph(X, q.thld = 0.9, pca.dim = NULL, verbose = FALSE)
    list(adj_list = g$mst_adj_list, weight_list = g$mst_weight_list)
}

sanitize.graph.weights <- function(graph, X = NULL, reweight = FALSE) {
    edges <- edge.table(graph$adj_list, graph$weight_list)
    if (!nrow(edges)) return(graph)
    if (reweight) {
        if (is.null(X)) stop("X is required when reweight = TRUE.")
        D <- euclidean.dist.matrix(X)
        edges$weight <- D[cbind(edges$from, edges$to)]
        edges <- edges[is.finite(edges$weight) & edges$weight >= 0, , drop = FALSE]
        return(graph.from.edges(length(graph$adj_list), edges))
    }
    bad <- !is.finite(edges$weight) | edges$weight < 0
    if (!any(bad)) return(graph)
    if (is.null(X)) {
        edges <- edges[!bad, , drop = FALSE]
    } else {
        D <- euclidean.dist.matrix(X)
        edges$weight[bad] <- D[cbind(edges$from[bad], edges$to[bad])]
        edges <- edges[is.finite(edges$weight) & edges$weight >= 0, , drop = FALSE]
    }
    graph.from.edges(length(graph$adj_list), edges)
}

local.scales <- function(X, mst.graph = NULL, source = c("knn", "mst"), k = 5L) {
    source <- match.arg(source)
    n <- nrow(X)
    D <- euclidean.dist.matrix(X)
    if (identical(source, "knn")) {
        k <- min(as.integer(k), n - 1L)
        apply(D, 1, function(x) sort(x[x > 0])[k])
    } else {
        if (is.null(mst.graph)) mst.graph <- mst.graph.from.cmst(X)
        s <- vapply(mst.graph$weight_list, function(w) if (length(w)) max(w) else NA_real_, numeric(1))
        s[is.na(s)] <- stats::median(s, na.rm = TRUE)
        s
    }
}

adaptive.threshold.graph <- function(X, mst.graph, scale.source = c("knn", "mst"),
                                     multiplier = 2, k.scale = 5L) {
    scale.source <- match.arg(scale.source)
    D <- euclidean.dist.matrix(X)
    n <- nrow(X)
    s <- local.scales(X, mst.graph, source = scale.source, k = k.scale)
    keep <- which(upper.tri(D) & D <= multiplier * pmax(outer(s, s, pmax)), arr.ind = TRUE)
    edges <- data.frame(from = keep[, 1], to = keep[, 2], weight = D[keep])
    union.graphs(mst.graph, graph.from.edges(n, edges))
}

iknn.graph <- function(X, k) {
    res <- create.iknn.graphs(
        X, kmin = k, kmax = k,
        max.path.edge.ratio.deviation.thld = 0,
        path.edge.ratio.percentile = 0.5,
        threshold.percentile = 0,
        compute.full = TRUE,
        with.isize.pruning = FALSE,
        pca.dim = NULL,
        n.cores = 1L,
        verbose = FALSE
    )
    sanitize.graph.weights(res$geom_pruned_graphs[[1L]], X = X, reweight = TRUE)
}

mknn.graph <- function(X, k) {
    g <- create.mknn.graph(X, k = k)
    sanitize.graph.weights(list(adj_list = g$adj_list, weight_list = g$weight_list), X = X, reweight = TRUE)
}

prune.redundant.edges <- function(graph, alt.path.ratio = 1.15) {
    edges <- edge.table(graph$adj_list, graph$weight_list)
    if (!nrow(edges)) return(graph)
    edges <- edges[order(edges$weight, decreasing = TRUE), ]
    kept <- edges
    for (r in seq_len(nrow(edges))) {
        key <- edge.key(edges$from[[r]], edges$to[[r]])
        kept.keys <- edge.key(kept$from, kept$to)
        idx <- match(key, kept.keys)
        if (is.na(idx) || nrow(kept) <= length(graph$adj_list) - 1L) next
        candidate <- kept[-idx, , drop = FALSE]
        g.cand <- graph.from.edges(length(graph$adj_list), candidate)
        if (component.count(g.cand$adj_list) > 1L) next
        d <- graph.distances(g.cand$adj_list, g.cand$weight_list)[edges$from[[r]], edges$to[[r]]]
        if (is.finite(d) && d <= alt.path.ratio * edges$weight[[r]]) {
            kept <- candidate
        }
    }
    graph.from.edges(length(graph$adj_list), kept)
}

cycle.repair.graph <- function(base.graph, X, n.add = 1L, local.multiplier = 2,
                               k.scale = 5L, min.ratio = 10) {
    n <- nrow(X)
    D <- euclidean.dist.matrix(X)
    base.dist <- graph.distances(base.graph$adj_list, base.graph$weight_list)
    s <- local.scales(X, base.graph, source = "knn", k = k.scale)
    existing <- edge.key(edge.table(base.graph$adj_list)$from, edge.table(base.graph$adj_list)$to)
    keep <- which(upper.tri(D), arr.ind = TRUE)
    cand <- data.frame(from = keep[, 1], to = keep[, 2], dE = D[keep])
    cand$key <- edge.key(cand$from, cand$to)
    cand <- cand[!(cand$key %in% existing), , drop = FALSE]
    cand$scale <- pmax(s[cand$from], s[cand$to])
    cand$dG <- base.dist[cbind(cand$from, cand$to)]
    cand$ratio <- cand$dG / pmax(cand$dE, .Machine$double.eps)
    cand <- cand[is.finite(cand$ratio) &
                     cand$dE <= local.multiplier * cand$scale &
                     cand$ratio >= min.ratio, , drop = FALSE]
    if (!nrow(cand)) return(base.graph)
    cand <- cand[order(cand$ratio, decreasing = TRUE), ]
    cand <- cand[seq_len(min(n.add, nrow(cand))), ]
    add.graph <- graph.from.edges(n, data.frame(from = cand$from, to = cand$to, weight = cand$dE))
    union.graphs(base.graph, add.graph)
}

oracle.edge.retention <- function(adj.list, oracle.adj.list, theta, radius = 1) {
    graph.edges <- edge.table(adj.list)
    oracle.edges <- edge.table(oracle.adj.list)
    graph.keys <- if (nrow(graph.edges)) edge.key(graph.edges$from, graph.edges$to) else character()
    oracle.keys <- edge.key(oracle.edges$from, oracle.edges$to)
    retained <- oracle.keys %in% graph.keys
    oracle.arc <- circle.dist.matrix(theta, radius = radius)[cbind(oracle.edges$from, oracle.edges$to)]
    data.frame(
        oracle_edge_count = nrow(oracle.edges),
        retained_oracle_edges = sum(retained),
        missing_oracle_edges = sum(!retained),
        max_missing_oracle_arc = if (any(!retained)) max(oracle.arc[!retained]) else 0,
        median_missing_oracle_arc = if (any(!retained)) stats::median(oracle.arc[!retained]) else 0
    )
}

geodesic.metrics <- function(adj.list, weight.list, true.dist, k.neighborhood = 10L) {
    graph.dist <- graph.distances(adj.list, weight.list)
    n <- nrow(true.dist)
    keep <- upper.tri(true.dist)
    finite <- keep & is.finite(graph.dist)
    finite.fraction <- sum(finite) / sum(keep)
    if (!any(finite)) {
        return(data.frame(
            finite_pair_fraction = finite.fraction,
            pearson = NA_real_, spearman = NA_real_, scale = NA_real_,
            median_relative_error = NA_real_, q90_relative_error = NA_real_,
            mean_neighbor_overlap = NA_real_
        ))
    }
    gd <- graph.dist[finite]
    td <- true.dist[finite]
    scale <- sum(gd * td) / sum(gd^2)
    rel.err <- abs(scale * gd - td) / pmax(td, .Machine$double.eps)
    rel.err <- rel.err[is.finite(rel.err)]
    neighbor.overlap <- vapply(seq_len(n), function(i) {
        true.order <- order(true.dist[i, ], seq_len(n))
        graph.order <- order(graph.dist[i, ], seq_len(n), na.last = NA)
        true.nn <- setdiff(true.order, i)[seq_len(min(k.neighborhood, n - 1L))]
        graph.nn <- setdiff(graph.order, i)[seq_len(min(k.neighborhood, length(graph.order) - 1L))]
        if (!length(graph.nn) || anyNA(graph.nn)) return(NA_real_)
        length(intersect(true.nn, graph.nn)) / length(true.nn)
    }, numeric(1))
    data.frame(
        finite_pair_fraction = finite.fraction,
        pearson = suppressWarnings(stats::cor(gd, td, method = "pearson")),
        spearman = suppressWarnings(stats::cor(gd, td, method = "spearman")),
        scale = scale,
        median_relative_error = if (length(rel.err)) stats::median(rel.err) else NA_real_,
        q90_relative_error = if (length(rel.err)) as.numeric(stats::quantile(rel.err, 0.9, names = FALSE, na.rm = TRUE)) else NA_real_,
        mean_neighbor_overlap = mean(neighbor.overlap, na.rm = TRUE)
    )
}

graph.summary.metrics <- function(adj.list, weight.list, true.dist = NULL,
                                  shortcut.threshold = pi / 4) {
    deg <- lengths(adj.list)
    edges <- edge.table(adj.list, weight.list)
    if (nrow(edges) && !is.null(true.dist)) {
        edge.true <- true.dist[cbind(edges$from, edges$to)]
        false.shortcut.rate <- mean(edge.true > shortcut.threshold)
        max.edge.true.geodesic <- max(edge.true)
    } else {
        false.shortcut.rate <- NA_real_
        max.edge.true.geodesic <- NA_real_
    }
    data.frame(
        n_edges = nrow(edges),
        n_components = component.count(adj.list),
        cycle_rank = cycle.rank(adj.list),
        mean_degree = mean(deg),
        min_degree = min(deg),
        max_degree = max(deg),
        edge_length_median = if (nrow(edges)) stats::median(edges$weight) else NA_real_,
        edge_length_max = if (nrow(edges)) max(edges$weight) else NA_real_,
        false_shortcut_rate = false.shortcut.rate,
        max_edge_true_geodesic = max.edge.true.geodesic
    )
}

score.graph.metrics <- function(metrics) {
    bad <- !is.finite(metrics$median_relative_error) | metrics$finite_pair_fraction < 1
    score <- metrics$median_relative_error +
        0.5 * metrics$q90_relative_error +
        2.0 * metrics$false_shortcut_rate +
        ifelse(metrics$n_components > 1, 10, 0) +
        ifelse(metrics$cycle_rank < 1, 2, 0) -
        0.25 * metrics$mean_neighbor_overlap
    score[bad] <- score[bad] + 10
    score
}

draw.graph.overlay <- function(X, values, adj.list, weight.list, main,
                               shortcut.fun = NULL, palette.name = "Viridis") {
    pal <- grDevices::hcl.colors(128, palette.name)
    z01 <- (values - min(values)) / max(diff(range(values)), .Machine$double.eps)
    cols <- pal[pmax(1L, pmin(128L, floor(z01 * 127) + 1L))]
    plot(X[, 1], X[, 2], asp = 1, pch = 19, col = cols, xlab = "x", ylab = "y",
         main = main, cex = 0.62)
    edges <- edge.table(adj.list, weight.list)
    if (nrow(edges)) {
        edge.cols <- rep("gray35", nrow(edges))
        edge.lwd <- rep(0.8, nrow(edges))
        if (!is.null(shortcut.fun)) {
            shortcut <- shortcut.fun(edges)
            edge.cols[shortcut] <- "red3"
            edge.lwd[shortcut] <- 2.1
        }
        segments(X[edges$from, 1], X[edges$from, 2], X[edges$to, 1], X[edges$to, 2],
                 col = edge.cols, lwd = edge.lwd)
        points(X[, 1], X[, 2], pch = 19, col = cols, cex = 0.56)
    }
}

draw.distance.scatter <- function(adj.list, weight.list, true.dist, main) {
    gd <- graph.distances(adj.list, weight.list)
    keep <- upper.tri(true.dist) & is.finite(gd)
    plot(true.dist[keep], gd[keep], pch = 16, cex = 0.35,
         col = grDevices::adjustcolor("black", alpha.f = 0.25),
         xlab = "true circular geodesic distance",
         ylab = "graph shortest-path distance",
         main = main)
    if (any(keep)) abline(stats::lm(gd[keep] ~ true.dist[keep]), col = "steelblue", lwd = 2)
}

tube.coords <- function(theta, u, radius = 1) {
    r <- radius + u
    cbind(x = r * cos(theta), y = r * sin(theta))
}

tube.density <- function(u, sigma, type = c("normal", "laplace"),
                         density.floor = 1e-4) {
    type <- match.arg(type)
    if (sigma <= 0) stop("sigma must be positive.")
    rho <- if (identical(type, "normal")) {
        exp(-(u / sigma)^2)
    } else {
        exp(-abs(u) / sigma)
    }
    pmax(rho, density.floor)
}

latent.tube.oracle <- function(theta.sample, u.sample, radius = 1,
                               sigma = 0.08,
                               density.type = c("normal", "laplace"),
                               alpha = 1,
                               n.theta = 360L,
                               u.max = 4 * sigma,
                               n.u = 41L,
                               density.floor = 1e-4,
                               include.diagonal = TRUE,
                               attach.extra = 1L) {
    density.type <- match.arg(density.type)
    theta.sample <- as.double(theta.sample) %% (2 * pi)
    u.sample <- as.double(u.sample)
    if (length(theta.sample) != length(u.sample)) {
        stop("theta.sample and u.sample must have the same length.")
    }
    n.theta <- as.integer(n.theta)
    n.u <- as.integer(n.u)
    attach.extra <- as.integer(attach.extra)
    if (n.theta < 8L || n.u < 3L) stop("oracle grid is too coarse.")

    theta.grid <- seq(0, 2 * pi, length.out = n.theta + 1L)[- (n.theta + 1L)]
    u.grid <- seq(-u.max, u.max, length.out = n.u)
    grid <- expand.grid(theta = theta.grid, u = u.grid, KEEP.OUT.ATTRS = FALSE)
    grid$kind <- "grid"
    n.grid <- nrow(grid)
    sample <- data.frame(theta = theta.sample, u = u.sample, kind = "sample")
    vertices <- rbind(grid, sample)
    coords <- tube.coords(vertices$theta, vertices$u, radius = radius)
    vertices$x <- coords[, 1]
    vertices$y <- coords[, 2]
    vertices$density <- tube.density(
        vertices$u,
        sigma = sigma,
        type = density.type,
        density.floor = density.floor
    )

    grid.id <- matrix(seq_len(n.grid), nrow = n.theta, ncol = n.u)
    add.edge <- function(from, to) {
        from <- as.integer(from)
        to <- as.integer(to)
        ok <- from != to
        from <- from[ok]
        to <- to[ok]
        if (!length(from)) {
            return(data.frame(from = integer(), to = integer()))
        }
        data.frame(from = pmin(from, to), to = pmax(from, to))
    }

    pieces <- list()
    ptr <- 0L
    for (a in seq_len(n.theta)) {
        b <- if (a == n.theta) 1L else a + 1L
        ptr <- ptr + 1L
        pieces[[ptr]] <- add.edge(grid.id[a, ], grid.id[b, ])
    }
    if (n.u > 1L) {
        for (v in seq_len(n.u - 1L)) {
            ptr <- ptr + 1L
            pieces[[ptr]] <- add.edge(grid.id[, v], grid.id[, v + 1L])
        }
    }
    if (include.diagonal && n.u > 1L) {
        for (a in seq_len(n.theta)) {
            b <- if (a == n.theta) 1L else a + 1L
            for (v in seq_len(n.u - 1L)) {
                ptr <- ptr + 1L
                pieces[[ptr]] <- add.edge(grid.id[a, v], grid.id[b, v + 1L])
                ptr <- ptr + 1L
                pieces[[ptr]] <- add.edge(grid.id[a, v + 1L], grid.id[b, v])
            }
        }
    }

    for (i in seq_along(theta.sample)) {
        sample.id <- n.grid + i
        theta.pos <- theta.sample[[i]] / (2 * pi) * n.theta
        theta.near <- ((floor(theta.pos) + (-attach.extra):(attach.extra + 1L)) %% n.theta) + 1L
        u.order <- order(abs(u.grid - u.sample[[i]]))
        u.near <- u.order[seq_len(min(length(u.order), 2L * attach.extra + 2L))]
        targets <- as.vector(grid.id[theta.near, u.near, drop = FALSE])
        ptr <- ptr + 1L
        pieces[[ptr]] <- add.edge(rep(sample.id, length(targets)), targets)
    }

    edges <- unique(do.call(rbind, pieces))
    mid.u <- (vertices$u[edges$from] + vertices$u[edges$to]) / 2
    rho <- tube.density(mid.u, sigma = sigma, type = density.type,
                        density.floor = density.floor)
    dx <- vertices$x[edges$from] - vertices$x[edges$to]
    dy <- vertices$y[edges$from] - vertices$y[edges$to]
    euclidean <- sqrt(dx^2 + dy^2)
    edges$weight <- euclidean / (rho^alpha)

    graph <- igraph::make_empty_graph(n = nrow(vertices), directed = FALSE)
    graph <- igraph::add_edges(graph, as.vector(t(as.matrix(edges[, c("from", "to")]))))
    igraph::E(graph)$weight <- edges$weight
    list(
        graph = graph,
        vertices = vertices,
        edges = edges,
        sample.vertices = n.grid + seq_along(theta.sample),
        params = list(
            radius = radius,
            sigma = sigma,
            density.type = density.type,
            alpha = alpha,
            n.theta = n.theta,
            n.u = n.u,
            u.max = u.max,
            density.floor = density.floor,
            include.diagonal = include.diagonal,
            attach.extra = attach.extra
        )
    )
}

latent.tube.oracle.distances <- function(oracle) {
    as.matrix(igraph::distances(
        oracle$graph,
        v = oracle$sample.vertices,
        to = oracle$sample.vertices,
        weights = igraph::E(oracle$graph)$weight
    ))
}

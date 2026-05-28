#' Construct PTTF Local Geometry
#'
#' Builds the geometry layer for the experimental parallel-transport trend
#' filtering (PTTF) route. The returned object contains graph supports, local
#' PCA tangent frames, edge directions, orthogonal edge transports, sampling
#' density metadata, and diagnostics. It does not assemble transported
#' difference operators and it does not fit a response.
#'
#' @param X Numeric data matrix with one row per vertex.
#' @param adj.list Optional supplied graph adjacency list using 1-based vertex
#'   indices.
#' @param weight.list Optional positive edge-length list parallel to
#'   \code{adj.list}.
#' @param graph Graph source. \code{"rknn"} builds an adaptive-radius graph from
#'   \code{X}; \code{"supplied"} uses \code{adj.list} and \code{weight.list}.
#' @param tangent.dim Fixed global tangent dimension for phase-1 geometry.
#' @param local.support Local support construction rule. \code{"graph.disk"}
#'   uses graph-hop disks; \code{"neighbors"} starts from the closed one-hop
#'   neighborhood and tops up by graph-hop expansion.
#' @param min.support Minimum support size. Defaults to
#'   \eqn{\max(2m+2,m+3)}, where \eqn{m=\code{tangent.dim}}.
#' @param max.hops Maximum hop radius used when topping up local supports.
#' @param support.k Optional cap on local support size after top-up. The center
#'   vertex is always retained and remaining vertices are selected by graph
#'   metric distance, with vertex-index tie breaking.
#' @param graph.k.scale,graph.radius.factor,graph.radius.rule Adaptive-radius
#'   graph controls passed to \code{\link{create.rknn.graph}} when
#'   \code{graph = "rknn"}.
#' @param transport.rule Edge transport rule. Phase 1 implements
#'   \code{"procrustes"}.
#' @param synchronize.orientation Logical. If \code{TRUE}, return deterministic
#'   spanning-tree orientation-synchronization metadata. Phase 1 does not mutate
#'   the returned frames.
#' @param density.method Sampling-density metadata rule. Phase 1 implements
#'   \code{"support.radius"}.
#' @param alpha Numeric exponent used only to report proposed future
#'   density-normalization weights.
#' @param diagnostics Logical. If \code{TRUE}, compute cycle diagnostics.
#'
#' @details
#' The stored transport convention is column-vector based. For an oriented edge
#' \eqn{i \leftarrow j}, \code{Oij} maps tangent coordinates from frame
#' \eqn{j} into frame \eqn{i}:
#' \deqn{v_j \mapsto O_{ij} v_j.}
#'
#' For shared-support Procrustes, row-coordinate matrices \eqn{Z_i} and
#' \eqn{Z_j} are formed on the same vertex set. The row problem
#' \eqn{\min_Q \|Z_j Q - Z_i\|_F^2} gives \eqn{Q = U V^\top} when
#' \eqn{Z_j^\top Z_i = U\Sigma V^\top}. The stored column map is
#' \eqn{O_{ij}=Q^\top=VU^\top}; the reported residual is
#' \eqn{\|Z_j O_{ij}^\top - Z_i\|_F / \max(\|Z_i\|_F,\epsilon)}.
#'
#' @return A list of class \code{"pttf_geometry"}.
#' @export
pttf.geometry <- function(
    X,
    adj.list = NULL,
    weight.list = NULL,
    graph = c("rknn", "supplied"),
    tangent.dim,
    local.support = c("graph.disk", "neighbors"),
    min.support = NULL,
    max.hops = 3L,
    support.k = NULL,
    graph.k.scale = 1L,
    graph.radius.factor = 1,
    graph.radius.rule = "geomean",
    transport.rule = c("procrustes"),
    synchronize.orientation = TRUE,
    density.method = c("support.radius"),
    alpha = 1,
    diagnostics = TRUE) {

    X <- .validate.numeric.data.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    graph <- match.arg(graph)
    local.support <- match.arg(local.support)
    transport.rule <- match.arg(transport.rule)
    density.method <- match.arg(density.method)

    tangent.dim <- .validate.positive.integer.scalar(tangent.dim, "tangent.dim")
    if (tangent.dim > p) {
        stop("'tangent.dim' cannot exceed the ambient dimension of 'X'.",
             call. = FALSE)
    }
    max.hops <- .validate.positive.integer.scalar(max.hops, "max.hops")
    if (is.null(min.support)) {
        min.support <- max(2L * tangent.dim + 2L, tangent.dim + 3L)
    } else {
        min.support <- .validate.positive.integer.scalar(min.support, "min.support")
    }
    if (min.support < tangent.dim + 2L) {
        stop("'min.support' must be at least tangent.dim + 2.", call. = FALSE)
    }
    if (!is.null(support.k)) {
        support.k <- .validate.positive.integer.scalar(support.k, "support.k")
        if (support.k < min.support) {
            stop("'support.k' must be at least 'min.support'.", call. = FALSE)
        }
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)) {
        stop("'alpha' must be a finite numeric scalar.", call. = FALSE)
    }

    g <- .pttf.geometry.prepare.graph(
        X = X,
        adj.list = adj.list,
        weight.list = weight.list,
        graph = graph,
        graph.k.scale = graph.k.scale,
        graph.radius.factor = graph.radius.factor,
        graph.radius.rule = graph.radius.rule
    )
    adj <- g$adj.list
    weights <- g$weight.list
    edge.table <- .pttf.geometry.edge.table(adj, weights)

    metric.distances <- .pttf.geometry.all.source.distances(adj, weights)
    supports <- .pttf.geometry.supports(
        adj = adj,
        metric.distances = metric.distances,
        local.support = local.support,
        min.support = min.support,
        max.hops = max.hops,
        support.k = support.k
    )
    frames <- .pttf.geometry.frames(
        X = X,
        supports = supports,
        tangent.dim = tangent.dim
    )
    density <- .pttf.geometry.density(
        supports = supports,
        tangent.dim = tangent.dim,
        alpha = alpha
    )
    transports <- .pttf.geometry.transports(
        X = X,
        adj = adj,
        weights = weights,
        edge.table = edge.table,
        supports = supports,
        frames = frames,
        tangent.dim = tangent.dim
    )
    sync <- .pttf.geometry.orientation.sync(
        adj = adj,
        transports = transports,
        tangent.dim = tangent.dim,
        enabled = isTRUE(synchronize.orientation)
    )
    cycle.diag <- if (isTRUE(diagnostics)) {
        .pttf.geometry.cycle.diagnostics(adj, transports, tangent.dim)
    } else {
        list(triangles = data.frame(), fundamental.cycles = data.frame())
    }

    out <- list(
        graph = c(
            list(
                adj.list = adj,
                weight.list = weights,
                edge.table = edge.table,
                constructor = g$constructor,
                parameters = g$parameters
            ),
            g$diagnostics
        ),
        frames = list(
            tangent.dim = tangent.dim,
            centers = X,
            frame.list = frames$frame.list,
            pca.spectra = frames$pca.spectra,
            support.index.list = supports$index.list,
            support.diagnostics = frames$support.diagnostics
        ),
        edges = transports$edges,
        transport = list(
            convention = "Oij maps frame j coordinates to frame i coordinates",
            oriented.edge.table = transports$oriented.edge.table,
            Oij = transports$Oij,
            uij = transports$uij,
            normalized.uij = transports$normalized.uij
        ),
        density = c(
            list(
                method = density.method,
                alpha = alpha,
                used.in.operator = FALSE
            ),
            density
        ),
        diagnostics = list(
            orientation.synchronization = sync,
            cycles = cycle.diag,
            status.summary = list(
                support = sort(table(frames$support.diagnostics$status)),
                edge = sort(table(transports$edges$status))
            )
        ),
        parameters = list(
            tangent.dim = tangent.dim,
            local.support = local.support,
            min.support = min.support,
            max.hops = max.hops,
            support.k = support.k,
            transport.rule = transport.rule,
            synchronize.orientation = isTRUE(synchronize.orientation),
            density.method = density.method
        )
    )
    class(out) <- c("pttf_geometry", "list")
    out
}

.pttf.geometry.prepare.graph <- function(X, adj.list, weight.list, graph,
                                         graph.k.scale, graph.radius.factor,
                                         graph.radius.rule) {
    if (identical(graph, "supplied") || !is.null(adj.list) || !is.null(weight.list)) {
        if (is.null(adj.list) || is.null(weight.list)) {
            stop("Supplied PTTF geometry graphs require both 'adj.list' and 'weight.list'.",
                 call. = FALSE)
        }
        validated <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
        return(list(
            adj.list = validated$adj.list,
            weight.list = validated$weight.list,
            constructor = "supplied",
            parameters = list(),
            diagnostics = list(
                n.vertices = length(validated$adj.list),
                n.edges = nrow(.pttf.geometry.edge.table(validated$adj.list,
                                                         validated$weight.list))
            )
        ))
    }

    graph.radius.rule <- match.arg(graph.radius.rule, c("max", "min", "geomean"))
    built <- create.rknn.graph(
        X,
        type = "adaptive.radius",
        k.scale = graph.k.scale,
        radius.factor = graph.radius.factor,
        radius.rule = graph.radius.rule,
        prune.method = "none",
        connect.components = TRUE,
        connect.method = "component.mst",
        graph.detail = "full"
    )
    validated <- .validate.metric.graph.lowpass.graph(built$adj_list,
                                                      built$weight_list)
    list(
        adj.list = validated$adj.list,
        weight.list = validated$weight.list,
        constructor = "create.rknn.graph",
        parameters = list(
            type = "adaptive.radius",
            k.scale = graph.k.scale,
            radius.factor = graph.radius.factor,
            radius.rule = graph.radius.rule,
            prune.method = "none",
            connect.components = TRUE,
            connect.method = "component.mst"
        ),
        diagnostics = list(
            n.vertices = built$n_vertices,
            n.edges = built$n_edges,
            n.components.before.mst = built$n_components_before,
            n.components.after.mst = built$n_components_after,
            n.mst.edges.added = built$n_mst_edges_added,
            n.edges.before.mst = built$n_edges_before_mst,
            n.edges.after.mst = built$n_edges_after_mst
        )
    )
}

.pttf.geometry.edge.table <- function(adj.list, weight.list) {
    rows <- vector("list", sum(lengths(adj.list)))
    cursor <- 0L
    for (i in seq_along(adj.list)) {
        nbrs <- adj.list[[i]]
        w <- weight.list[[i]]
        keep <- nbrs > i
        if (any(keep)) {
            for (k in which(keep)) {
                cursor <- cursor + 1L
                rows[[cursor]] <- data.frame(
                    from = as.integer(i),
                    to = as.integer(nbrs[k]),
                    length = as.numeric(w[k])
                )
            }
        }
    }
    if (!cursor) {
        return(data.frame(from = integer(), to = integer(), length = numeric()))
    }
    out <- do.call(rbind, rows[seq_len(cursor)])
    rownames(out) <- NULL
    out
}

.pttf.geometry.hop.distances <- function(adj, start, max.hops) {
    n <- length(adj)
    dist <- rep.int(Inf, n)
    dist[start] <- 0
    frontier <- start
    for (h in seq_len(max.hops)) {
        if (!length(frontier)) break
        next.nodes <- sort(unique(unlist(adj[frontier], use.names = FALSE)))
        next.nodes <- next.nodes[is.infinite(dist[next.nodes])]
        if (!length(next.nodes)) {
            frontier <- integer()
            next
        }
        dist[next.nodes] <- h
        frontier <- next.nodes
    }
    dist
}

.pttf.geometry.single.source.distances <- function(adj, weights, start) {
    n <- length(adj)
    dist <- rep.int(Inf, n)
    visited <- rep.int(FALSE, n)
    dist[start] <- 0
    for (step in seq_len(n)) {
        candidates <- which(!visited & is.finite(dist))
        if (!length(candidates)) break
        u <- candidates[which.min(dist[candidates])]
        visited[u] <- TRUE
        if (!length(adj[[u]])) next
        for (k in seq_along(adj[[u]])) {
            v <- adj[[u]][k]
            alt <- dist[u] + weights[[u]][k]
            if (alt < dist[v]) dist[v] <- alt
        }
    }
    dist
}

.pttf.geometry.all.source.distances <- function(adj, weights) {
    n <- length(adj)
    out <- matrix(Inf, nrow = n, ncol = n)
    for (i in seq_len(n)) {
        out[i, ] <- .pttf.geometry.single.source.distances(adj, weights, i)
    }
    out
}

.pttf.geometry.supports <- function(adj, metric.distances, local.support,
                                    min.support, max.hops, support.k) {
    n <- length(adj)
    index.list <- vector("list", n)
    hop.radius <- integer(n)
    metric.radius <- numeric(n)
    initial.size <- integer(n)
    final.size <- integer(n)
    raw.status <- character(n)

    for (i in seq_len(n)) {
        hdist <- .pttf.geometry.hop.distances(adj, i, max.hops)
        if (identical(local.support, "neighbors")) {
            initial <- sort(unique(c(i, adj[[i]])))
            initial.size[i] <- length(initial)
            if (length(initial) >= min.support) {
                ids <- initial
                h.used <- 1L
            } else {
                ids <- initial
                h.used <- 1L
                for (h in seq.int(2L, max.hops)) {
                    ids <- which(hdist <= h)
                    h.used <- h
                    if (length(ids) >= min.support) break
                }
            }
        } else {
            ids <- integer()
            initial.size[i] <- sum(hdist <= 1L)
            h.used <- 1L
            for (h in seq_len(max.hops)) {
                ids <- which(hdist <= h)
                h.used <- h
                if (length(ids) >= min.support) break
            }
        }

        if (!is.null(support.k) && length(ids) > support.k) {
            candidates <- setdiff(ids, i)
            ord <- order(metric.distances[i, candidates], candidates)
            ids <- sort(c(i, candidates[ord][seq_len(support.k - 1L)]))
        } else {
            ids <- sort(ids)
        }
        index.list[[i]] <- ids
        hop.radius[i] <- h.used
        final.size[i] <- length(ids)
        finite.dist <- metric.distances[i, ids]
        metric.radius[i] <- if (length(finite.dist) && any(is.finite(finite.dist))) {
            max(finite.dist[is.finite(finite.dist)])
        } else {
            0
        }
        raw.status[i] <- if (!length(ids)) {
            "empty_support"
        } else if (length(ids) < min.support) {
            "too_few_support_points"
        } else {
            "pending"
        }
    }

    list(
        index.list = index.list,
        diagnostics = data.frame(
            vertex = seq_len(n),
            initial.size = initial.size,
            support.size = final.size,
            hop.radius = hop.radius,
            metric.radius = metric.radius,
            raw.status = raw.status,
            stringsAsFactors = FALSE
        )
    )
}

.pttf.geometry.frames <- function(X, supports, tangent.dim) {
    n <- nrow(X)
    p <- ncol(X)
    frame.list <- vector("list", n)
    spectra <- matrix(NA_real_, nrow = n, ncol = p)
    colnames(spectra) <- paste0("sv", seq_len(p))
    diag <- supports$diagnostics
    diag$rank <- integer(n)
    diag$condition <- rep(Inf, n)
    diag$status <- diag$raw.status

    for (i in seq_len(n)) {
        ids <- supports$index.list[[i]]
        frame.list[[i]] <- matrix(NA_real_, nrow = p, ncol = tangent.dim)
        if (identical(diag$raw.status[i], "empty_support") ||
            identical(diag$raw.status[i], "too_few_support_points")) {
            next
        }
        centered <- sweep(X[ids, , drop = FALSE], 2L, X[i, ], "-")
        sv <- tryCatch(svd(centered, nu = 0L, nv = tangent.dim),
                       error = function(e) NULL)
        if (is.null(sv) || !length(sv$d) || !all(is.finite(sv$d))) {
            diag$status[i] <- "rank_deficient"
            next
        }
        spectra[i, seq_along(sv$d)] <- sv$d
        s1 <- max(sv$d[1L], .Machine$double.eps)
        rank <- sum(sv$d > 1e-10 * s1)
        diag$rank[i] <- rank
        if (rank < tangent.dim || ncol(sv$v) < tangent.dim) {
            diag$status[i] <- "rank_deficient"
            next
        }
        cond <- sv$d[1L] / max(sv$d[tangent.dim], .Machine$double.eps)
        diag$condition[i] <- cond
        if (cond > 1e8) {
            diag$status[i] <- "ill_conditioned"
            next
        }
        frame.list[[i]] <- sv$v[, seq_len(tangent.dim), drop = FALSE]
        diag$status[i] <- "ok"
    }

    list(
        frame.list = frame.list,
        pca.spectra = as.data.frame(spectra),
        support.diagnostics = diag[, c("vertex", "support.size", "initial.size",
                                       "hop.radius", "metric.radius", "rank",
                                       "condition", "status")]
    )
}

.pttf.geometry.density <- function(supports, tangent.dim, alpha) {
    diag <- supports$diagnostics
    radius <- pmax(diag$metric.radius, 1e-12)
    rho <- diag$support.size / (radius^tangent.dim)
    rho <- pmax(rho, 1e-12)
    node.mass <- 1 / rho
    list(
        rho.hat = rho,
        node.mass = node.mass,
        proposed.alpha.weight = node.mass^alpha
    )
}

.pttf.geometry.valid.frame <- function(frame) {
    is.matrix(frame) && all(is.finite(frame))
}

.pttf.geometry.orthogonal.factor <- function(M) {
    sv <- svd(M)
    sv$u %*% t(sv$v)
}

.pttf.geometry.matrix.condition <- function(M) {
    sv <- svd(M, nu = 0L, nv = 0L)$d
    if (!length(sv)) return(Inf)
    max(sv) / max(min(sv), .Machine$double.eps)
}

.pttf.geometry.procrustes.edge <- function(X, i, j, supports, frames,
                                           tangent.dim) {
    Fi <- frames$frame.list[[i]]
    Fj <- frames$frame.list[[j]]
    I <- diag(1, tangent.dim)
    if (!.pttf.geometry.valid.frame(Fi) || !.pttf.geometry.valid.frame(Fj)) {
        return(list(O = I, residual = NA_real_, shared = 0L,
                    frame.residual = NA_real_,
                    frame.inverse.residual = NA_real_,
                    status = "invalid_frame"))
    }

    frame.map <- .pttf.geometry.orthogonal.factor(t(Fi) %*% Fj)
    frame.inverse.residual <- norm(
        frame.map - t(.pttf.geometry.orthogonal.factor(t(Fj) %*% Fi)),
        type = "F"
    )

    shared <- sort(intersect(supports$index.list[[i]], supports$index.list[[j]]))
    if (length(shared) < tangent.dim + 1L) {
        return(list(O = frame.map, residual = NA_real_, shared = length(shared),
                    frame.residual = 0,
                    frame.inverse.residual = frame.inverse.residual,
                    status = "too_few_shared_points"))
    }

    Zi <- sweep(X[shared, , drop = FALSE], 2L, X[i, ], "-") %*% Fi
    Zj <- sweep(X[shared, , drop = FALSE], 2L, X[j, ], "-") %*% Fj
    Zi.centered <- scale(Zi, center = TRUE, scale = FALSE)
    Zj.centered <- scale(Zj, center = TRUE, scale = FALSE)
    sv.i <- svd(Zi.centered, nu = 0L, nv = 0L)$d
    sv.j <- svd(Zj.centered, nu = 0L, nv = 0L)$d
    rank.i <- sum(sv.i > 1e-10 * max(sv.i, 0))
    rank.j <- sum(sv.j > 1e-10 * max(sv.j, 0))
    if (min(rank.i, rank.j) < tangent.dim) {
        return(list(O = frame.map, residual = NA_real_, shared = length(shared),
                    frame.residual = 0,
                    frame.inverse.residual = frame.inverse.residual,
                    status = "shared_rank_deficient"))
    }

    cross <- crossprod(Zj, Zi)
    if (.pttf.geometry.matrix.condition(cross) > 1e8) {
        return(list(O = frame.map, residual = NA_real_, shared = length(shared),
                    frame.residual = 0,
                    frame.inverse.residual = frame.inverse.residual,
                    status = "shared_ill_conditioned"))
    }
    sv <- svd(cross)
    Q <- sv$u %*% t(sv$v)
    O <- t(Q)
    residual <- norm(Zj %*% t(O) - Zi, type = "F") /
        max(norm(Zi, type = "F"), .Machine$double.eps)
    frame.residual <- norm(frame.map - O, type = "F") / sqrt(tangent.dim)
    list(O = O, residual = residual, shared = length(shared),
         frame.residual = frame.residual,
         frame.inverse.residual = frame.inverse.residual,
         status = "ok")
}

.pttf.geometry.transports <- function(X, adj, weights, edge.table, supports,
                                      frames, tangent.dim) {
    O.list <- list()
    u.list <- list()
    un.list <- list()
    oriented <- vector("list", 2L * nrow(edge.table))
    edge.rows <- vector("list", nrow(edge.table))
    cursor <- 0L

    for (e in seq_len(nrow(edge.table))) {
        i <- edge.table$from[e]
        j <- edge.table$to[e]
        pr <- .pttf.geometry.procrustes.edge(X, i, j, supports, frames,
                                             tangent.dim)
        Oij <- pr$O
        Oji <- t(Oij)
        inv.res <- norm(Oji - t(Oij), type = "F")

        for (pair in list(c(i, j), c(j, i))) {
            a <- pair[1L]
            b <- pair[2L]
            cursor <- cursor + 1L
            key <- paste(a, b, sep = "->")
            O <- if (a == i && b == j) Oij else Oji
            O.list[[key]] <- O
            oriented[[cursor]] <- data.frame(from = a, to = b, key = key,
                                             edge.index = e,
                                             stringsAsFactors = FALSE)
            Fa <- frames$frame.list[[a]]
            vec <- X[b, ] - X[a, ]
            if (.pttf.geometry.valid.frame(Fa)) {
                u <- as.vector(crossprod(Fa, vec))
            } else {
                u <- rep(NA_real_, tangent.dim)
            }
            u.list[[key]] <- u
            un.list[[key]] <- if (all(is.finite(u)) && sqrt(sum(u^2)) > 0) {
                u / sqrt(sum(u^2))
            } else {
                rep(NA_real_, tangent.dim)
            }
        }

        edge.rows[[e]] <- data.frame(
            from = i,
            to = j,
            length = edge.table$length[e],
            support.shared = pr$shared,
            procrustes.residual = pr$residual,
            frame.transport.residual = pr$frame.residual,
            frame.inverse.residual = pr$frame.inverse.residual,
            inverse.residual = inv.res,
            det.transport = det(Oij),
            status = pr$status,
            stringsAsFactors = FALSE
        )
    }

    list(
        edges = if (length(edge.rows)) do.call(rbind, edge.rows) else data.frame(),
        oriented.edge.table = if (length(oriented)) do.call(rbind, oriented)
                              else data.frame(),
        Oij = O.list,
        uij = u.list,
        normalized.uij = un.list
    )
}

.pttf.geometry.orientation.sync <- function(adj, transports, tangent.dim, enabled) {
    n <- length(adj)
    if (!enabled) {
        return(list(enabled = FALSE, method = "none"))
    }
    parent <- rep.int(NA_integer_, n)
    order <- integer()
    queue <- 1L
    parent[1L] <- 0L
    gauge <- vector("list", n)
    gauge[[1L]] <- diag(1, tangent.dim)
    while (length(queue)) {
        u <- queue[1L]
        queue <- queue[-1L]
        order <- c(order, u)
        for (v in sort(adj[[u]])) {
            if (is.na(parent[v])) {
                parent[v] <- u
                key <- paste(u, v, sep = "->")
                gauge[[v]] <- gauge[[u]] %*% transports$Oij[[key]]
                queue <- c(queue, v)
            }
        }
    }
    n.reached <- length(order)
    non.tree <- .pttf.geometry.non.tree.sync.residuals(
        adj = adj,
        parent = parent,
        gauge = gauge,
        transports = transports,
        tangent.dim = tangent.dim
    )
    non.tree.summary <- .pttf.geometry.non.tree.summary(non.tree)

    if (tangent.dim == 1L) {
        sign <- rep.int(NA_real_, n)
        sign[order] <- vapply(gauge[order], function(G) {
            val <- G[1L, 1L]
            if (is.finite(val) && val < 0) -1 else 1
        }, numeric(1L))
        return(list(enabled = TRUE, method = "bfs.sign", root = 1L,
                    metadata.only = TRUE,
                    n.reached = n.reached,
                    n.flips = sum(sign[order] < 0),
                    parent = parent,
                    sign = sign,
                    non.tree.edge.residuals = non.tree,
                    non.tree.edge.summary = non.tree.summary))
    }
    list(enabled = TRUE, method = "bfs.metadata.only", root = 1L,
         metadata.only = TRUE,
         n.reached = n.reached,
         n.flips = NA_integer_,
         parent = parent,
         non.tree.edge.residuals = non.tree,
         non.tree.edge.summary = non.tree.summary)
}

.pttf.geometry.non.tree.sync.residuals <- function(adj, parent, gauge, transports,
                                                  tangent.dim) {
    rows <- list()
    cursor <- 0L
    for (i in seq_along(adj)) {
        for (j in adj[[i]][adj[[i]] > i]) {
            is.tree <- identical(parent[j], i) || identical(parent[i], j)
            if (is.tree || is.na(parent[i]) || is.na(parent[j])) next
            key <- paste(i, j, sep = "->")
            stored <- transports$Oij[[key]]
            if (is.null(stored) || is.null(gauge[[i]]) || is.null(gauge[[j]])) next
            predicted <- t(gauge[[i]]) %*% gauge[[j]]
            cursor <- cursor + 1L
            rows[[cursor]] <- data.frame(
                from = i,
                to = j,
                residual = norm(stored - predicted, type = "F") / sqrt(tangent.dim)
            )
        }
    }
    if (length(rows)) do.call(rbind, rows) else data.frame(
        from = integer(), to = integer(), residual = numeric()
    )
}

.pttf.geometry.non.tree.summary <- function(non.tree) {
    if (!nrow(non.tree)) {
        return(list(n.edges = 0L, median.residual = NA_real_,
                    max.residual = NA_real_))
    }
    list(
        n.edges = nrow(non.tree),
        median.residual = stats::median(non.tree$residual, na.rm = TRUE),
        max.residual = max(non.tree$residual, na.rm = TRUE)
    )
}

.pttf.geometry.get.O <- function(transports, from, to, tangent.dim) {
    key <- paste(to, from, sep = "->")
    O <- transports$Oij[[key]]
    if (is.null(O)) diag(1, tangent.dim) else O
}

.pttf.geometry.cycle.product.residual <- function(path, transports, tangent.dim) {
    prod <- diag(1, tangent.dim)
    for (k in seq_len(length(path) - 1L)) {
        prod <- .pttf.geometry.get.O(transports, path[k], path[k + 1L],
                                     tangent.dim) %*% prod
    }
    norm(prod - diag(1, tangent.dim), type = "F") / sqrt(tangent.dim)
}

.pttf.geometry.cycle.diagnostics <- function(adj, transports, tangent.dim,
                                             max.cycles = 500L) {
    n <- length(adj)
    triangles <- list()
    cursor <- 0L
    for (i in seq_len(n)) {
        for (j in adj[[i]][adj[[i]] > i]) {
            common <- intersect(adj[[i]], adj[[j]])
            for (k in common[common > j]) {
                cursor <- cursor + 1L
                if (cursor > max.cycles) break
                path <- c(i, j, k, i)
                triangles[[cursor]] <- data.frame(
                    v1 = i, v2 = j, v3 = k,
                    residual = .pttf.geometry.cycle.product.residual(
                        path, transports, tangent.dim
                    )
                )
            }
            if (cursor > max.cycles) break
        }
        if (cursor > max.cycles) break
    }
    tri <- if (length(triangles)) do.call(rbind, triangles) else data.frame()

    fund <- .pttf.geometry.fundamental.cycles(adj, transports, tangent.dim,
                                             max.cycles = max.cycles)
    list(triangles = tri, fundamental.cycles = fund)
}

.pttf.geometry.fundamental.cycles <- function(adj, transports, tangent.dim,
                                             max.cycles = 500L) {
    n <- length(adj)
    parent <- rep.int(NA_integer_, n)
    depth <- rep.int(Inf, n)
    parent[1L] <- 0L
    depth[1L] <- 0L
    queue <- 1L
    tree.edges <- character()
    while (length(queue)) {
        u <- queue[1L]
        queue <- queue[-1L]
        for (v in sort(adj[[u]])) {
            if (is.na(parent[v])) {
                parent[v] <- u
                depth[v] <- depth[u] + 1L
                tree.edges <- c(tree.edges, paste(min(u, v), max(u, v), sep = "-"))
                queue <- c(queue, v)
            }
        }
    }
    rows <- list()
    cursor <- 0L
    for (i in seq_len(n)) {
        for (j in adj[[i]][adj[[i]] > i]) {
            ekey <- paste(i, j, sep = "-")
            if (ekey %in% tree.edges) next
            path <- .pttf.geometry.tree.path(i, j, parent, depth)
            if (length(path) < 2L) next
            cycle.path <- c(path, i)
            cursor <- cursor + 1L
            if (cursor > max.cycles) break
            rows[[cursor]] <- data.frame(
                from = i,
                to = j,
                length = length(cycle.path) - 1L,
                residual = .pttf.geometry.cycle.product.residual(
                    cycle.path, transports, tangent.dim
                )
            )
        }
        if (cursor > max.cycles) break
    }
    if (length(rows)) do.call(rbind, rows) else data.frame()
}

.pttf.geometry.tree.path <- function(i, j, parent, depth) {
    if (!is.finite(depth[i]) || !is.finite(depth[j])) return(integer())
    left <- i
    right <- j
    a <- i
    b <- j
    while (depth[a] > depth[b]) {
        a <- parent[a]
        left <- c(left, a)
    }
    while (depth[b] > depth[a]) {
        b <- parent[b]
        right <- c(b, right)
    }
    while (a != b) {
        a <- parent[a]
        b <- parent[b]
        left <- c(left, a)
        right <- c(b, right)
    }
    c(left, right[-1L])
}

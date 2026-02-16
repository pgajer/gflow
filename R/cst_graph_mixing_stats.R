#' CST mixing statistics on a graph with optional weight handling and permutation null
#'
#' @description
#' Computes label mixing metrics for a graph:
#'   - weighted edge homophily (same-label edge fraction),
#'   - weighted nominal assortativity (categorical; mixing-matrix formulation),
#'   - per-label conductance (weighted cut / min(vol, vol_complement)),
#'   - permutation null distributions with effect sizes and z-scores.
#'
#' Graph can be provided as an igraph object or as (adj.list, weight.list).
#'
#' @param igraph.obj igraph graph object (preferred). If NULL, adj.list is used.
#' @param adj.list Optional adjacency list: list of integer neighbor indices (1..n).
#' @param weight.list Optional weight list parallel to adj.list (same structure).
#' @param labels Vector of labels aligned to vertices (length vcount(g)); NA allowed.
#' @param n.perm Integer number of permutations for null distribution.
#' @param perm.blocks Optional vector/factor (length vcount(g)) to permute labels within blocks.
#' @param use.weights Logical; if TRUE, uses E(g)$weight when present.
#' @param weights.are.edge.lengths Logical; if TRUE, interpret weights as edge lengths and convert to affinities.
#' @param affinity.method Character: "exp" for exp(-(d/sigma)^2) or "inv" for 1/(d+eps).
#' @param sigma Optional positive numeric; length scale for exp affinity. If NULL, median of positive lengths.
#' @param affinity.eps Small positive scalar for inverse affinity.
#' @param simplify.multiple Logical; if TRUE, removes loops and merges multi-edges.
#' @param seed Integer RNG seed.
#'
#' @return Object of class "cst_graph_mixing_stats".
cst.graph.mixing.stats <- function(igraph.obj = NULL,
                                  adj.list = NULL,
                                  weight.list = NULL,
                                  labels,
                                  n.perm = 200L,
                                  perm.blocks = NULL,
                                  use.weights = TRUE,
                                  weights.are.edge.lengths = FALSE,
                                  affinity.method = c("exp", "inv"),
                                  sigma = NULL,
                                  affinity.eps = 1e-8,
                                  simplify.multiple = TRUE,
                                  seed = 1L) {

    ## ---- build / validate graph ----
    affinity.method <- match.arg(affinity.method)

    adjlist.weightlist.to.igraph <- function(adj.list, weight.list = NULL) {
        if (is.null(adj.list) || !is.list(adj.list)) stop("`adj.list` must be a list.")
        n <- length(adj.list)
        if (n < 2L) stop("`adj.list` must have length >= 2.")

        has.w <- !is.null(weight.list)
        if (has.w) {
            if (!is.list(weight.list) || length(weight.list) != n) stop("`weight.list` must be a list of same length as adj.list.")
        }

        ## collect edges i<j to avoid duplicates
        e1 <- integer(0); e2 <- integer(0); ew <- numeric(0)

        for (i in seq_len(n)) {
            nb <- adj.list[[i]]
            if (length(nb) == 0L) next
            nb <- as.integer(nb)
            if (any(nb < 1L | nb > n)) stop("`adj.list` has out-of-range neighbor indices.")
            if (!has.w) {
                jj <- nb[nb > i]
                if (length(jj) > 0L) {
                    e1 <- c(e1, rep.int(i, length(jj)))
                    e2 <- c(e2, jj)
                    ew <- c(ew, rep.int(1.0, length(jj)))
                }
            } else {
                wv <- weight.list[[i]]
                if (length(wv) != length(nb)) stop("weight.list[[i]] length must match adj.list[[i]].")
                wv <- as.double(wv)
                keep <- (nb > i)
                if (any(keep)) {
                    e1 <- c(e1, rep.int(i, sum(keep)))
                    e2 <- c(e2, nb[keep])
                    ew <- c(ew, wv[keep])
                }
            }
        }

        if (length(e1) == 0L) {
            g <- igraph::make_empty_graph(n = n, directed = FALSE)
            return(g)
        }

        em <- cbind(e1, e2)
        g <- igraph::make_empty_graph(n = n, directed = FALSE)
        g <- igraph::add_edges(g, as.vector(t(em)))
        igraph::E(g)$weight <- ew
        g
    }

    if (!is.null(igraph.obj)) {
        if (!inherits(igraph.obj, "igraph")) stop("`igraph.obj` must be an igraph object.")
        g <- igraph.obj
    } else {
        g <- adjlist.weightlist.to.igraph(adj.list, weight.list)
    }

    n <- igraph::vcount(g)
    m <- igraph::ecount(g)
    if (missing(labels) || is.null(labels)) stop("`labels` must be provided.")
    if (length(labels) != n) stop("`labels` must have length vcount(g).")

    labels <- as.character(labels)

    ## ---- handle multi-edges / loops ----
    w.raw <- NULL
    if (isTRUE(use.weights) && "weight" %in% igraph::edge_attr_names(g)) {
        w.raw <- as.double(igraph::E(g)$weight)
        if (length(w.raw) != m) w.raw <- NULL
    }
    if (is.null(w.raw)) w.raw <- rep(1.0, m)

    ## simplify: decide how to combine weights
    if (isTRUE(simplify.multiple)) {
        comb.fun <- if (isTRUE(weights.are.edge.lengths)) "min" else "max"
        ## if unweighted, "max" is fine (all 1)
        g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                              edge.attr.comb = list(weight = comb.fun, "ignore"))
        ## refresh after simplify
        n <- igraph::vcount(g)
        m <- igraph::ecount(g)
        w.raw <- if ("weight" %in% igraph::edge_attr_names(g)) as.double(igraph::E(g)$weight) else rep(1.0, m)
        if (length(w.raw) != m) w.raw <- rep(1.0, m)
    }

    ## ---- convert weights if needed (lengths -> affinities) ----
    w.used <- w.raw
    w.aff <- NULL

    if (isTRUE(use.weights)) {
        if (isTRUE(weights.are.edge.lengths)) {
            d <- w.raw
            d <- d[is.finite(d) & d > 0]
            if (length(d) == 0L) {
                w.aff <- rep(1.0, m)
            } else {
                sig <- sigma
                if (is.null(sig)) sig <- stats::median(d)
                if (!is.finite(sig) || sig <= 0) sig <- stats::median(d)
                if (!is.finite(sig) || sig <= 0) sig <- 1.0

                if (affinity.method == "exp") {
                    w.aff <- exp(- (w.raw / sig)^2)
                } else {
                    w.aff <- 1 / (w.raw + affinity.eps)
                }
                w.aff[!is.finite(w.aff)] <- 0
                w.aff[w.aff < 0] <- 0
            }
            w.used <- w.aff
        } else {
            ## treat as affinity/similarity already
            w.used <- w.raw
            w.used[!is.finite(w.used)] <- 0
            w.used[w.used < 0] <- 0
        }
    } else {
        w.used <- rep(1.0, m)
    }

    ## ---- core edge endpoint arrays ----
    ends <- igraph::ends(g, igraph::E(g), names = FALSE)
    u <- ends[, 1]; v <- ends[, 2]
    lab.u <- labels[u]; lab.v <- labels[v]
    ok.e <- !is.na(lab.u) & !is.na(lab.v)

    ## labeled vertex set
    ok.v <- !is.na(labels)
    idx.v <- which(ok.v)

    ## ---- weighted homophily ----
    w.ok <- w.used[ok.e]
    same <- (lab.u[ok.e] == lab.v[ok.e])
    homophily <- if (sum(w.ok) > 0) sum(w.ok[same]) / sum(w.ok) else NA_real_

    ## ---- weighted mixing matrix + assortativity (nominal; weighted) ----
    assort <- NA_real_
    mix.mat <- NULL
    lev <- sort(unique(labels[ok.v]))
    if (length(lev) >= 2L && sum(ok.e) > 0L) {
        map <- setNames(seq_along(lev), lev)
        a <- map[lab.u[ok.e]]
        b <- map[lab.v[ok.e]]

        M <- matrix(0.0, nrow = length(lev), ncol = length(lev),
                    dimnames = list(lev, lev))

        ## for undirected edges, add weight to both (a,b) and (b,a)
        for (i in seq_along(a)) {
            w0 <- w.used[which(ok.e)[i]]
            M[a[i], b[i]] <- M[a[i], b[i]] + w0
            M[b[i], a[i]] <- M[b[i], a[i]] + w0
        }
        sM <- sum(M)
        if (sM > 0) {
            e <- M / sM
            aa <- rowSums(e)
            tr <- sum(diag(e))
            denom <- 1 - sum(aa^2)
            assort <- if (denom > 0) (tr - sum(aa^2)) / denom else NA_real_
        }
        mix.mat <- M
    }

    ## ---- strength-weighted homophily null baseline and adjusted homophily ----
    strength <- igraph::strength(g, weights = w.used)
    strength[!is.finite(strength)] <- 0
    strength.ok <- strength[ok.v]
    labels.ok <- labels[ok.v]

    h.null <- NA_real_
    if (sum(strength.ok) > 0 && length(lev) >= 1L) {
        p.l <- tapply(strength.ok, labels.ok, sum)
        p.l <- p.l / sum(strength.ok)
        h.null <- sum(as.numeric(p.l)^2)
    }
    homophily.adjusted <- if (is.finite(h.null) && h.null < 1 && is.finite(homophily)) (homophily - h.null) / (1 - h.null) else NA_real_

    ## ---- conductance per label (weighted) ----
    conductance.by.label <- NULL
    conductance.summary <- NULL
    if (length(lev) >= 2L && sum(ok.e) > 0L) {
        ## compute vol and cut with weights
        vol.total <- sum(strength.ok)
        vol.l <- tapply(strength.ok, labels.ok, sum)
        vol.l <- vol.l[lev]

        cut.l <- setNames(rep(0.0, length(lev)), lev)
        for (i in which(ok.e)) {
            a0 <- labels[u[i]]
            b0 <- labels[v[i]]
            if (!is.na(a0) && !is.na(b0) && a0 != b0) {
                w0 <- w.used[i]
                cut.l[a0] <- cut.l[a0] + w0
                cut.l[b0] <- cut.l[b0] + w0
            }
        }

        cond <- rep(NA_real_, length(lev)); names(cond) <- lev
        for (l in lev) {
            va <- as.numeric(vol.l[l])
            vb <- vol.total - va
            denom <- min(va, vb)
            cond[l] <- if (is.finite(denom) && denom > 0) as.numeric(cut.l[l]) / denom else NA_real_
        }

        conductance.by.label <- data.frame(
            cst = lev,
            vol = as.numeric(vol.l),
            cut = as.numeric(cut.l[lev]),
            conductance = as.numeric(cond),
            stringsAsFactors = FALSE
        )

        wv <- conductance.by.label$vol
        ok.c <- is.finite(conductance.by.label$conductance) & wv > 0
        conductance.summary <- list(
            conductance.median = stats::median(conductance.by.label$conductance[ok.c], na.rm = TRUE),
            conductance.vol.weighted.mean = sum(conductance.by.label$conductance[ok.c] * wv[ok.c]) / sum(wv[ok.c])
        )
    }

    ## ---- permutation nulls ----
    n.perm <- as.integer(n.perm)
    perm <- NULL
    if (n.perm > 0L && sum(ok.v) >= 10L && sum(ok.e) >= 10L) {

        permute.labels <- function(lbl, blocks = NULL) {
            lbl2 <- lbl
            idx <- which(!is.na(lbl2))
            if (is.null(blocks)) {
                lbl2[idx] <- sample(lbl2[idx], replace = FALSE)
            } else {
                if (length(blocks) != length(lbl2)) stop("perm.blocks must have length vcount(g).")
                b <- blocks
                for (bb in unique(b[idx])) {
                    ii <- idx[b[idx] == bb]
                    if (length(ii) >= 2L) lbl2[ii] <- sample(lbl2[ii], replace = FALSE)
                }
            }
            lbl2
        }

        metric.from.labels <- function(lbl) {
            ## recompute homophily + assortativity with fixed weights and edges
            lab.u2 <- lbl[u]; lab.v2 <- lbl[v]
            ok.e2 <- !is.na(lab.u2) & !is.na(lab.v2)
            if (sum(ok.e2) == 0L) return(c(h = NA_real_, r = NA_real_))

            w.ok2 <- w.used[ok.e2]
            same2 <- (lab.u2[ok.e2] == lab.v2[ok.e2])
            h2 <- if (sum(w.ok2) > 0) sum(w.ok2[same2]) / sum(w.ok2) else NA_real_

            ## weighted assortativity
            ok.v2 <- !is.na(lbl)
            lev2 <- sort(unique(lbl[ok.v2]))
            r2 <- NA_real_
            if (length(lev2) >= 2L) {
                map2 <- setNames(seq_along(lev2), lev2)
                a2 <- map2[lab.u2[ok.e2]]
                b2 <- map2[lab.v2[ok.e2]]
                M2 <- matrix(0.0, nrow = length(lev2), ncol = length(lev2))
                idx.e2 <- which(ok.e2)
                for (ii in seq_along(a2)) {
                    w0 <- w.used[idx.e2[ii]]
                    M2[a2[ii], b2[ii]] <- M2[a2[ii], b2[ii]] + w0
                    M2[b2[ii], a2[ii]] <- M2[b2[ii], a2[ii]] + w0
                }
                sM2 <- sum(M2)
                if (sM2 > 0) {
                    e2 <- M2 / sM2
                    aa2 <- rowSums(e2)
                    tr2 <- sum(diag(e2))
                    denom2 <- 1 - sum(aa2^2)
                    r2 <- if (denom2 > 0) (tr2 - sum(aa2^2)) / denom2 else NA_real_
                }
            }
            c(h = h2, r = r2)
        }

        set.seed(seed)
        h.null.vec <- rep(NA_real_, n.perm)
        r.null.vec <- rep(NA_real_, n.perm)

        for (b in seq_len(n.perm)) {
            lbl.p <- permute.labels(labels, blocks = perm.blocks)
            mm <- metric.from.labels(lbl.p)
            h.null.vec[b] <- mm["h"]
            r.null.vec[b] <- mm["r"]
        }

        summarize.null <- function(obs, null.vec) {
            mu <- mean(null.vec, na.rm = TRUE)
            sd0 <- stats::sd(null.vec, na.rm = TRUE)
            eff <- obs - mu
            z <- if (is.finite(sd0) && sd0 > 0) eff / sd0 else NA_real_
            null.finite <- null.vec[is.finite(null.vec)]
            p.upper <- if (length(null.finite) >= 10L) (1 + sum(null.finite >= obs)) / (1 + length(null.finite)) else NA_real_
            list(obs = obs, mu = mu, sd = sd0, effect = eff, z = z, p = p.upper)
        }

        perm <- list(
            n.perm = n.perm,
            homophily.null = h.null.vec,
            assortativity.null = r.null.vec,
            homophily = summarize.null(homophily, h.null.vec),
            assortativity = summarize.null(assort, r.null.vec)
        )
    }

    out <- list(
        n.vertices = n,
        n.edges = m,
        labels.levels = lev,
        weights.used = isTRUE(use.weights),
        weights.are.edge.lengths = isTRUE(weights.are.edge.lengths),
        affinity.method = affinity.method,
        sigma = sigma,
        homophily = homophily,
        homophily.null = h.null,
        homophily.adjusted = homophily.adjusted,
        assortativity = assort,
        mixing.matrix = mix.mat,
        conductance.by.label = conductance.by.label,
        conductance.summary = conductance.summary,
        permutation = perm
    )
    class(out) <- "cst_graph_mixing_stats"
    out
}

#' Plot method for cst_graph_mixing_stats
#'
#' @param x Object from cst.graph.mixing.stats().
#' @param which Character vector specifying panels to plot.
#'   Options: "null.homophily", "null.assortativity", "mixing.matrix", "conductance".
#' @param ... Base graphics args.
#' @export
plot.cst_graph_mixing_stats <- function(x,
                                       which = c("null.homophily", "null.assortativity", "mixing.matrix", "conductance"),
                                       ...) {
    if (!inherits(x, "cst_graph_mixing_stats")) stop("x must be class 'cst_graph_mixing_stats'.")

    which <- unique(which)
    np <- length(which)
    if (np == 0L) return(invisible(NULL))

    ## layout
    nr <- ceiling(np / 2)
    nc <- if (np == 1L) 1 else 2
    op <- graphics::par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2.0, 0.8), mgp = c(2.0, 0.6, 0))
    on.exit(graphics::par(op), add = TRUE)

    for (w in which) {
        if (w == "mixing.matrix") {
            M <- x$mixing.matrix
            if (is.null(M)) {
                graphics::plot.new(); graphics::title("mixing.matrix (none)")
            } else {
                ## simple heatmap via image
                z <- log1p(M)
                graphics::image(t(z[nrow(z):1, , drop = FALSE]),
                                axes = FALSE, main = "log1p(weighted mixing matrix)")
                graphics::axis(1, at = seq(0, 1, length.out = ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.6)
                graphics::axis(2, at = seq(0, 1, length.out = nrow(z)), labels = rev(rownames(z)), las = 2, cex.axis = 0.6)
            }
        } else if (w == "conductance") {
            df <- x$conductance.by.label
            if (is.null(df)) {
                graphics::plot.new(); graphics::title("conductance (none)")
            } else {
                o <- order(df$conductance, decreasing = FALSE, na.last = NA)
                graphics::barplot(df$conductance[o], names.arg = df$cst[o], las = 2,
                                  ylab = "conductance", main = "Per-CST conductance", cex.names = 0.6)
            }
        } else if (w == "null.homophily") {
            if (is.null(x$permutation)) {
                graphics::plot.new(); graphics::title("homophily null (none)")
            } else {
                h0 <- x$permutation$homophily.null
                graphics::hist(h0, breaks = 30, col = "gray", border = "white",
                               main = "Homophily null", xlab = "homophily")
                graphics::abline(v = x$homophily, lwd = 2)
                graphics::mtext(sprintf("obs=%.3f  z=%.2f  p=%.3g",
                                        x$permutation$homophily$obs,
                                        x$permutation$homophily$z,
                                        x$permutation$homophily$p),
                                side = 3, line = 0.1, cex = 0.8)
            }
        } else if (w == "null.assortativity") {
            if (is.null(x$permutation)) {
                graphics::plot.new(); graphics::title("assortativity null (none)")
            } else {
                r0 <- x$permutation$assortativity.null
                graphics::hist(r0, breaks = 30, col = "gray", border = "white",
                               main = "Assortativity null", xlab = "assortativity")
                graphics::abline(v = x$assortativity, lwd = 2)
                graphics::mtext(sprintf("obs=%.3f  z=%.2f  p=%.3g",
                                        x$permutation$assortativity$obs,
                                        x$permutation$assortativity$z,
                                        x$permutation$assortativity$p),
                                side = 3, line = 0.1, cex = 0.8)
            }
        } else {
            graphics::plot.new()
            graphics::title(paste0("Unknown panel: ", w))
        }
    }

    invisible(NULL)
}

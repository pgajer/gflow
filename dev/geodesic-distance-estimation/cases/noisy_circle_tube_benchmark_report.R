# Density-aware noisy-circle / latent-tube benchmark report.

options(stringsAsFactors = FALSE)

args <- commandArgs(FALSE)
file.arg <- args[grepl("^--file=", args)]
script.path <- if (length(file.arg)) sub("^--file=", "", file.arg[[1L]]) else ""
root.dir <- if (nzchar(script.path)) {
    normalizePath(file.path(dirname(script.path), "../../.."), mustWork = TRUE)
} else {
    normalizePath(getwd(), mustWork = TRUE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to run this development report.")
}
pkgload::load_all(root.dir, quiet = TRUE)
source(file.path(root.dir, "dev/geodesic-distance-estimation/R/graph_construction.R"))

out.root <- file.path(root.dir, "dev/geodesic-distance-estimation")
notes.dir <- file.path(out.root, "notes")
fig.dir <- file.path(out.root, "figures/noisy-circle-tube-benchmark")
cache.dir <- file.path(out.root, "cache/noisy-circle-tube-benchmark")
dir.create(notes.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
    radius = 1,
    n = 160L,
    seed = 4101L,
    sigma = 0.08,
    noise_type = "normal",
    alpha_grid = c(0, 0.25, 0.5, 0.75, 1),
    primary_alpha = 0.5,
    stress_alpha = 1,
    density_floor_grid = c(1e-4, 1e-3, 1e-2),
    penalty_cap_grid = c(Inf, 20, 10),
    attach_extra_grid = c(0L, 1L, 2L),
    primary_floor = 1e-4,
    primary_penalty_cap = Inf,
    primary_attach_extra = 1L,
    oracle_n_theta = 360L,
    oracle_n_u = 41L,
    oracle_u_max_multiplier = 4,
    graph_k = 8L,
    iknn_pruning_delta = 0.10,
    diffusion_t = 8L,
    diffusion_k = 10L
)

md.path <- file.path(notes.dir, "density_aware_geodesic_distance_noisy_circle_tube_benchmark_report.md")

md.escape <- function(x) {
    x <- as.character(x)
    x <- gsub("\\|", "\\\\|", x)
    x <- gsub("\n", " ", x, fixed = TRUE)
    x
}

fmt <- function(x, digits = 4) {
    ifelse(
        is.na(x), "NA",
        ifelse(is.infinite(x), as.character(x),
               formatC(x, digits = digits, format = "fg", flag = "#"))
    )
}

md.table <- function(x, digits = 4) {
    if (is.null(x) || !nrow(x)) return("_No rows._\n")
    y <- x
    for (j in seq_along(y)) {
        if (is.numeric(y[[j]])) y[[j]] <- fmt(y[[j]], digits = digits)
        y[[j]] <- md.escape(y[[j]])
    }
    header <- paste(names(y), collapse = " | ")
    sep <- paste(rep("---", ncol(y)), collapse = " | ")
    rows <- apply(y, 1, paste, collapse = " | ")
    paste(c(header, sep, rows), collapse = "\n")
}

short.table <- function(x, names) {
    names(x) <- names
    x
}

fig.link <- function(path) {
    file.path("../figures/noisy-circle-tube-benchmark", basename(path))
}

make.case <- function(params) {
    X.df <- generate.circle.data(
        n = params$n,
        radius = params$radius,
        noise = params$sigma,
        type = "random",
        noise.type = params$noise_type,
        seed = params$seed
    )
    X <- as.matrix(X.df[, c("x", "y")])
    theta <- as.double(X.df$angles) %% (2 * pi)
    u <- sqrt(rowSums(X^2)) - params$radius
    list(
        X.df = X.df,
        X = X,
        theta = theta,
        u = u,
        arc.dist = circle.dist.matrix(theta, radius = params$radius),
        euclidean.dist = euclidean.dist.matrix(X)
    )
}

reweight.oracle <- function(oracle, alpha, penalty.cap = Inf) {
    edges <- oracle$edges
    vertices <- oracle$vertices
    mid.u <- (vertices$u[edges$from] + vertices$u[edges$to]) / 2
    rho <- tube.density(
        mid.u,
        sigma = oracle$params$sigma,
        type = oracle$params$density.type,
        density.floor = oracle$params$density.floor
    )
    dx <- vertices$x[edges$from] - vertices$x[edges$to]
    dy <- vertices$y[edges$from] - vertices$y[edges$to]
    euclidean <- sqrt(dx^2 + dy^2)
    penalty <- rho^(-alpha)
    if (is.finite(penalty.cap)) penalty <- pmin(penalty, penalty.cap)
    edges$base_weight <- euclidean
    edges$density <- rho
    edges$penalty <- penalty
    edges$weight <- euclidean * penalty
    igraph::E(oracle$graph)$weight <- edges$weight
    oracle$edges <- edges
    oracle$params$alpha <- alpha
    oracle$params$penalty.cap <- penalty.cap
    oracle
}

build.oracle <- function(case, alpha, density.floor, penalty.cap, attach.extra, params) {
    oracle <- latent.tube.oracle(
        theta.sample = case$theta,
        u.sample = case$u,
        radius = params$radius,
        sigma = params$sigma,
        density.type = "normal",
        alpha = alpha,
        n.theta = params$oracle_n_theta,
        u.max = params$oracle_u_max_multiplier * params$sigma,
        n.u = params$oracle_n_u,
        density.floor = density.floor,
        include.diagonal = TRUE,
        attach.extra = attach.extra
    )
    reweight.oracle(oracle, alpha = alpha, penalty.cap = penalty.cap)
}

oracle.dist <- function(oracle) {
    latent.tube.oracle.distances(oracle)
}

distance.summary <- function(name, D, case, truth = NULL) {
    keep <- upper.tri(D) & is.finite(D)
    d <- D[keep]
    arc <- case$arc.dist[keep]
    eu <- case$euclidean.dist[keep]
    if (is.null(truth)) truth <- case$arc.dist
    td <- truth[keep]
    scale <- sum(d * td) / sum(d^2)
    rel.err <- abs(scale * d - td) / pmax(td, .Machine$double.eps)
    data.frame(
        distance = name,
        finite_pair_fraction = sum(keep) / sum(upper.tri(D)),
        median = stats::median(d),
        q90 = as.numeric(stats::quantile(d, 0.9, names = FALSE)),
        max = max(d),
        spearman_vs_arc = suppressWarnings(stats::cor(d, arc, method = "spearman")),
        spearman_vs_euclidean = suppressWarnings(stats::cor(d, eu, method = "spearman")),
        spearman_vs_primary_truth = suppressWarnings(stats::cor(d, td, method = "spearman")),
        median_relative_error_to_primary = stats::median(rel.err),
        q90_relative_error_to_primary = as.numeric(stats::quantile(rel.err, 0.9, names = FALSE))
    )
}

base.cmst.graph <- function(X, q.thld = 0.5) {
    g <- create.cmst.graph(X, q.thld = q.thld, pca.dim = NULL, verbose = FALSE)
    sanitize.graph.weights(
        list(adj_list = g$cmst_adj_list, weight_list = g$cmst_weight_list),
        X = X,
        reweight = TRUE
    )
}

base.iknn.graph <- function(X, k, delta) {
    g <- create.iknn.graphs(
        X,
        kmin = k,
        kmax = k,
        max.path.edge.ratio.deviation.thld = delta,
        path.edge.ratio.percentile = 0.5,
        threshold.percentile = 0,
        compute.full = TRUE,
        with.isize.pruning = FALSE,
        pca.dim = NULL,
        n.cores = 1L,
        verbose = FALSE
    )$geom_pruned_graphs[[1L]]
    sanitize.graph.weights(g, X = X, reweight = TRUE)
}

safe.mknn.graph <- function(X, k) {
    g <- create.mknn.graph(X, k = k)
    sanitize.graph.weights(list(adj_list = g$adj_list, weight_list = g$weight_list), X = X, reweight = TRUE)
}

density.weight.sample.graph <- function(graph, case, alpha, density.floor, penalty.cap = Inf) {
    edges <- edge.table(graph$adj_list, graph$weight_list)
    mid.u <- (case$u[edges$from] + case$u[edges$to]) / 2
    rho <- tube.density(mid.u, sigma = params$sigma, type = "normal", density.floor = density.floor)
    penalty <- rho^(-alpha)
    if (is.finite(penalty.cap)) penalty <- pmin(penalty, penalty.cap)
    edges$base_weight <- edges$weight
    edges$density <- rho
    edges$penalty <- penalty
    edges$weight <- edges$base_weight * penalty
    out <- graph.from.edges(length(graph$adj_list), edges[, c("from", "to", "weight")])
    out$edge_density <- edges
    out
}

fermat.graph.distances <- function(graph, p, rooted = TRUE) {
    edges <- edge.table(graph$adj_list, graph$weight_list)
    edges$weight <- edges$weight^p
    g <- graph.from.edges(length(graph$adj_list), edges)
    D <- graph.distances(g$adj_list, g$weight_list)
    if (rooted) D <- D^(1 / p)
    D
}

kernel.transition <- function(D, k = 10L, t = 8L) {
    n <- nrow(D)
    diag(D) <- Inf
    sigma <- apply(D, 1, function(x) sort(x)[min(k, length(x) - 1L)])
    eps <- stats::median(sigma)^2
    K <- exp(-(D^2) / pmax(eps, .Machine$double.eps))
    diag(K) <- 0
    keep <- t(apply(D, 1, function(x) rank(x, ties.method = "first") <= k))
    K[!keep] <- 0
    K <- pmax(K, t(K))
    rs <- rowSums(K)
    rs[rs == 0] <- 1
    P <- K / rs
    Pt <- diag(n)
    for (i in seq_len(t)) Pt <- Pt %*% P
    Pt
}

row.euclidean.dist <- function(M) {
    as.matrix(stats::dist(M))
}

path.info <- function(graph, from, to, case) {
    g <- igraph.from.graph(graph$adj_list, graph$weight_list)
    sp <- igraph::shortest_paths(g, from = from, to = to, weights = igraph::E(g)$weight)$vpath[[1L]]
    v <- as.integer(sp)
    if (length(v) < 2L) {
        return(list(vertices = v, edges = data.frame(), summary = data.frame(
            n_edges = 0L, base_length = NA_real_, weighted_length = NA_real_,
            min_edge_density = NA_real_, max_edge_penalty = NA_real_
        )))
    }
    e <- data.frame(from = head(v, -1L), to = tail(v, -1L))
    e$key <- edge.key(e$from, e$to)
    all.edges <- edge.table(graph$adj_list, graph$weight_list)
    all.edges$key <- edge.key(all.edges$from, all.edges$to)
    idx <- match(e$key, all.edges$key)
    e$weight <- all.edges$weight[idx]
    e$base_weight <- if ("edge_density" %in% names(graph)) {
        de <- graph$edge_density
        de$key <- edge.key(de$from, de$to)
        de$base_weight[match(e$key, de$key)]
    } else {
        case$euclidean.dist[cbind(e$from, e$to)]
    }
    if ("edge_density" %in% names(graph)) {
        de <- graph$edge_density
        de$key <- edge.key(de$from, de$to)
        e$density <- de$density[match(e$key, de$key)]
        e$penalty <- de$penalty[match(e$key, de$key)]
    } else {
        e$density <- NA_real_
        e$penalty <- NA_real_
    }
    finite.min <- function(x) {
        x <- x[is.finite(x)]
        if (length(x)) min(x) else NA_real_
    }
    finite.max <- function(x) {
        x <- x[is.finite(x)]
        if (length(x)) max(x) else NA_real_
    }
    list(
        vertices = v,
        edges = e,
        summary = data.frame(
            n_edges = nrow(e),
            base_length = sum(e$base_weight, na.rm = TRUE),
            weighted_length = sum(e$weight, na.rm = TRUE),
            min_edge_density = finite.min(e$density),
            max_edge_penalty = finite.max(e$penalty)
        )
    )
}

select.pairs <- function(case, primary.truth, graph.dist) {
    keep <- which(upper.tri(primary.truth), arr.ind = TRUE)
    arc <- case$arc.dist[keep]
    eu <- case$euclidean.dist[keep]
    truth <- primary.truth[keep]
    ratio <- truth / pmax(eu, .Machine$double.eps)
    short.eu <- eu <= stats::quantile(eu, 0.15)
    radial <- arc <= stats::quantile(arc, 0.15)
    gd <- graph.dist[keep]
    under <- truth / pmax(gd, .Machine$double.eps)
    idx <- c(
        which.max(arc),
        which.max(ifelse(short.eu, ratio, -Inf)),
        which.max(ifelse(radial, abs(case$u[keep[, 1]] - case$u[keep[, 2]]), -Inf)),
        which.max(ifelse(is.finite(under), under, -Inf))
    )
    labels <- c("opposite-angle", "euclidean-close/high-truth-ratio",
                "near-angle/radial-separation", "largest-graph-underestimate")
    data.frame(
        label = labels,
        from = keep[idx, 1],
        to = keep[idx, 2],
        arc = arc[idx],
        euclidean = eu[idx],
        primary_truth = truth[idx],
        base_graph_distance = gd[idx],
        truth_to_euclidean = ratio[idx],
        truth_to_graph = under[idx],
        u_from = case$u[keep[idx, 1]],
        u_to = case$u[keep[idx, 2]]
    )
}

case <- make.case(params)

primary.oracle <- build.oracle(
    case,
    alpha = params$primary_alpha,
    density.floor = params$primary_floor,
    penalty.cap = params$primary_penalty_cap,
    attach.extra = params$primary_attach_extra,
    params = params
)
primary.truth <- oracle.dist(primary.oracle)

alpha.rows <- lapply(params$alpha_grid, function(alpha) {
    o <- build.oracle(case, alpha, params$primary_floor, Inf, params$primary_attach_extra, params)
    D <- oracle.dist(o)
    cbind(data.frame(alpha = alpha, density_floor = params$primary_floor,
                     penalty_cap = Inf, attach_extra = params$primary_attach_extra),
          distance.summary(sprintf("tube alpha=%s", alpha), D, case, primary.truth))
})
alpha.summary <- do.call(rbind, alpha.rows)

floor.cap.grid <- expand.grid(
    alpha = c(params$primary_alpha, params$stress_alpha),
    density_floor = params$density_floor_grid,
    penalty_cap = params$penalty_cap_grid,
    KEEP.OUT.ATTRS = FALSE
)
floor.cap.rows <- lapply(seq_len(nrow(floor.cap.grid)), function(i) {
    row <- floor.cap.grid[i, ]
    o <- build.oracle(case, row$alpha, row$density_floor, row$penalty_cap,
                      params$primary_attach_extra, params)
    D <- oracle.dist(o)
    cbind(row, attach_extra = params$primary_attach_extra,
          distance.summary("oracle variant", D, case, primary.truth))
})
floor.cap.summary <- do.call(rbind, floor.cap.rows)

attach.rows <- lapply(params$attach_extra_grid, function(ae) {
    o <- build.oracle(case, params$primary_alpha, params$primary_floor, Inf, ae, params)
    D <- oracle.dist(o)
    cbind(data.frame(alpha = params$primary_alpha, density_floor = params$primary_floor,
                     penalty_cap = Inf, attach_extra = ae),
          distance.summary("attachment variant", D, case, primary.truth))
})
attachment.summary <- do.call(rbind, attach.rows)

mst.graph <- mst.graph.from.cmst(case$X)
cmst.graph <- base.cmst.graph(case$X, q.thld = 0.5)
iknn.graph.primary <- base.iknn.graph(case$X, k = params$graph_k, delta = params$iknn_pruning_delta)
mknn.graph.primary <- safe.mknn.graph(case$X, k = params$graph_k)
density.iknn <- density.weight.sample.graph(
    iknn.graph.primary,
    case,
    alpha = params$primary_alpha,
    density.floor = params$primary_floor,
    penalty.cap = Inf
)
density.iknn.capped <- density.weight.sample.graph(
    iknn.graph.primary,
    case,
    alpha = params$stress_alpha,
    density.floor = params$primary_floor,
    penalty.cap = 20
)

graph.distances.list <- list(
    Euclidean = case$euclidean.dist,
    `latent arc` = case$arc.dist,
    MST = graph.distances(mst.graph$adj_list, mst.graph$weight_list),
    CMST = graph.distances(cmst.graph$adj_list, cmst.graph$weight_list),
    `iKNN base` = graph.distances(iknn.graph.primary$adj_list, iknn.graph.primary$weight_list),
    `mKNN base` = graph.distances(mknn.graph.primary$adj_list, mknn.graph.primary$weight_list),
    `iKNN density alpha=0.5` = graph.distances(density.iknn$adj_list, density.iknn$weight_list),
    `iKNN density alpha=1 capped` = graph.distances(density.iknn.capped$adj_list, density.iknn.capped$weight_list),
    `Fermat rooted p=2 on iKNN` = fermat.graph.distances(iknn.graph.primary, p = 2, rooted = TRUE)
)

Pt <- kernel.transition(case$euclidean.dist, k = params$diffusion_k, t = params$diffusion_t)
graph.distances.list$`diffusion row distance` <- row.euclidean.dist(Pt)
potential <- -log(pmax(Pt, 1e-12))
graph.distances.list$`PHATE-like log potential` <- row.euclidean.dist(potential)

method.summary <- do.call(rbind, lapply(names(graph.distances.list), function(name) {
    distance.summary(name, graph.distances.list[[name]], case, primary.truth)
}))
method.summary <- method.summary[order(method.summary$median_relative_error_to_primary), ]

edge.attr <- density.iknn$edge_density
edge.attr$arc <- case$arc.dist[cbind(edge.attr$from, edge.attr$to)]
edge.attr$primary_truth <- primary.truth[cbind(edge.attr$from, edge.attr$to)]
edge.attr$truth_to_euclidean <- edge.attr$primary_truth / pmax(edge.attr$base_weight, .Machine$double.eps)
edge.attr$inflation <- edge.attr$weight / pmax(edge.attr$base_weight, .Machine$double.eps)
top.edge.attr <- edge.attr[order(edge.attr$truth_to_euclidean, decreasing = TRUE), ]
top.edge.attr <- top.edge.attr[seq_len(min(12L, nrow(top.edge.attr))), ]

pair.table <- select.pairs(case, primary.truth, graph.distances.list$`iKNN base`)
path.rows <- list()
for (r in seq_len(nrow(pair.table))) {
    for (graph.name in c("iKNN base", "iKNN density alpha=0.5")) {
        graph <- if (identical(graph.name, "iKNN base")) iknn.graph.primary else density.iknn
        info <- path.info(graph, pair.table$from[[r]], pair.table$to[[r]], case)
        row <- cbind(pair.table[r, c("label", "from", "to", "primary_truth")],
                     graph = graph.name,
                     info$summary)
        path.rows[[length(path.rows) + 1L]] <- row
    }
}
path.summary <- do.call(rbind, path.rows)

fig.files <- list(
    geometry = file.path(fig.dir, "01_geometry_density.png"),
    oracle.sensitivity = file.path(fig.dir, "02_oracle_sensitivity.png"),
    method.comparison = file.path(fig.dir, "03_method_comparison.png"),
    edge.attribution = file.path(fig.dir, "04_edge_attribution.png"),
    path.attribution = file.path(fig.dir, "05_path_attribution.png")
)

save.png(fig.files$geometry, {
    rho <- tube.density(case$u, sigma = params$sigma, type = "normal",
                        density.floor = params$primary_floor)
    pal <- grDevices::hcl.colors(128, "YlGnBu", rev = TRUE)
    z <- (rho - min(rho)) / max(diff(range(rho)), .Machine$double.eps)
    cols <- pal[pmax(1L, pmin(128L, floor(z * 127) + 1L))]
    oldpar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    plot(case$X[, 1], case$X[, 2], asp = 1, pch = 19, col = cols,
         xlab = "x", ylab = "y", main = "Observed noisy circle, colored by density")
    th <- seq(0, 2 * pi, length.out = 400)
    lines(cos(th), sin(th), col = "gray25", lwd = 1.5)
    plot(case$theta, case$u, pch = 19, col = cols,
         xlab = "theta", ylab = "u = ||x|| - radius",
         main = "Latent tube coordinates")
    abline(h = 0, col = "gray50", lty = 2)
}, width = 1600, height = 800)

save.png(fig.files$oracle.sensitivity, {
    oldpar <- par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    plot(alpha.summary$alpha, alpha.summary$median, type = "b", pch = 19,
         xlab = "alpha", ylab = "median oracle distance",
         main = "Oracle distance grows with alpha")
    plot(alpha.summary$alpha, alpha.summary$max, type = "b", pch = 19,
         xlab = "alpha", ylab = "max oracle distance",
         main = "Tail sensitivity")
    stress <- floor.cap.summary[floor.cap.summary$alpha == params$stress_alpha, ]
    stress$cap.label <- ifelse(is.infinite(stress$penalty_cap), "Inf", as.character(stress$penalty_cap))
    cols <- setNames(c("#333333", "#1f77b4", "#d95f02"), unique(stress$cap.label))
    plot(NA, xlim = range(log10(params$density_floor_grid)), ylim = range(stress$max),
         xlab = "log10 density floor", ylab = "max distance",
         main = "Alpha=1 floor/cap control")
    for (cap in unique(stress$cap.label)) {
        x <- stress[stress$cap.label == cap, ]
        x <- x[order(x$density_floor), ]
        lines(log10(x$density_floor), x$max, type = "b", pch = 19, col = cols[[cap]])
    }
    legend("topright", legend = paste("cap", names(cols)), col = cols, lty = 1, pch = 19, bty = "n")
    plot(attachment.summary$attach_extra, attachment.summary$spearman_vs_primary_truth,
         type = "b", pch = 19, xlab = "attach_extra", ylab = "Spearman vs primary truth",
         main = "Attachment sensitivity")
}, width = 1600, height = 1200)

save.png(fig.files$method.comparison, {
    oldpar <- par(mfrow = c(1, 2), mar = c(8.0, 4.8, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    top <- method.summary
    top$distance <- factor(top$distance, levels = top$distance)
    barplot(top$median_relative_error_to_primary, names.arg = top$distance,
            las = 2, cex.names = 0.72, ylab = "median relative error",
            main = "Distance recovery vs tube oracle alpha=0.5",
            col = "#7aa6c2")
    barplot(top$spearman_vs_primary_truth, names.arg = top$distance,
            las = 2, cex.names = 0.72, ylab = "Spearman correlation",
            main = "Rank agreement with tube oracle alpha=0.5",
            col = "#9ccf8b")
}, width = 1700, height = 950)

save.png(fig.files$edge.attribution, {
    oldpar <- par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    plot(edge.attr$base_weight, edge.attr$primary_truth, pch = 16, cex = 0.55,
         col = grDevices::adjustcolor("black", alpha.f = 0.45),
         xlab = "Euclidean edge length", ylab = "tube-oracle edge distance",
         main = "Edge truth vs Euclidean length")
    abline(0, 1, col = "gray55", lty = 2)
    plot(edge.attr$density, edge.attr$inflation, pch = 16, cex = 0.55,
         col = grDevices::adjustcolor("firebrick", alpha.f = 0.45),
         xlab = "mid-edge tube density", ylab = "edge inflation ratio",
         main = "Density penalty on iKNN edges")
}, width = 1500, height = 750)

save.png(fig.files$path.attribution, {
    oldpar <- par(mfrow = c(nrow(pair.table), 2), mar = c(2.6, 2.6, 2.6, 0.8))
    on.exit(par(oldpar), add = TRUE)
    for (r in seq_len(nrow(pair.table))) {
        for (graph.name in c("iKNN base", "iKNN density alpha=0.5")) {
            graph <- if (identical(graph.name, "iKNN base")) iknn.graph.primary else density.iknn
            info <- path.info(graph, pair.table$from[[r]], pair.table$to[[r]], case)
            plot(case$X[, 1], case$X[, 2], asp = 1, pch = 19, cex = 0.35,
                 col = "gray75", xlab = "", ylab = "",
                 main = sprintf("%s: %s", pair.table$label[[r]], graph.name))
            if (length(info$vertices) > 1L) {
                lines(case$X[info$vertices, 1], case$X[info$vertices, 2], col = "red3", lwd = 2.2)
            }
            points(case$X[c(pair.table$from[[r]], pair.table$to[[r]]), 1],
                   case$X[c(pair.table$from[[r]], pair.table$to[[r]]), 2],
                   pch = 21, bg = c("dodgerblue3", "goldenrod2"), col = "black", cex = 1.0)
        }
    }
}, width = 1500, height = 2100)

saveRDS(
    list(
        params = params,
        case = case,
        alpha.summary = alpha.summary,
        floor.cap.summary = floor.cap.summary,
        attachment.summary = attachment.summary,
        method.summary = method.summary,
        edge.attr = edge.attr,
        pair.table = pair.table,
        path.summary = path.summary
    ),
    file.path(cache.dir, "noisy_circle_tube_benchmark_results.rds")
)

utils::write.csv(alpha.summary, file.path(cache.dir, "oracle_alpha_summary.csv"), row.names = FALSE)
utils::write.csv(floor.cap.summary, file.path(cache.dir, "oracle_floor_cap_summary.csv"), row.names = FALSE)
utils::write.csv(attachment.summary, file.path(cache.dir, "oracle_attachment_summary.csv"), row.names = FALSE)
utils::write.csv(method.summary, file.path(cache.dir, "method_summary.csv"), row.names = FALSE)
utils::write.csv(edge.attr, file.path(cache.dir, "edge_attribution.csv"), row.names = FALSE)
utils::write.csv(path.summary, file.path(cache.dir, "path_summary.csv"), row.names = FALSE)

best.method <- method.summary[1L, ]
alpha05 <- alpha.summary[alpha.summary$alpha == params$primary_alpha, ]
alpha1 <- alpha.summary[alpha.summary$alpha == params$stress_alpha, ]

report <- paste0(
    "# Density-Aware Geodesic Distance Noisy-Circle / Tube Benchmark Report\n\n",
    "<style>\n",
    "table { font-size: 10px; width: 100%; }\n",
    "th, td { padding: 2px 4px; }\n",
    "img { max-width: 100%; }\n",
    "code { white-space: pre-wrap; }\n",
    "</style>\n\n",
    "Date: ", format(Sys.time(), "%Y-%m-%d"), "\n\n",
    "This report is Deliverable 4 for the density-aware geodesic-distance research track. ",
    "It turns the noisy-circle latent-tube oracle from a sanity check into a benchmark ",
    "with density sensitivity, floor/cap/attachment sweeps, distance-family comparisons, ",
    "and edge/path attribution.\n\n",
    "## 1. Benchmark Setup\n\n",
    "The observed data are generated as a noisy circle with radius ", params$radius,
    ", radial noise scale \\(\\sigma=", params$sigma, "\\), sample size \\(n=",
    params$n, "\\), and seed `", params$seed, "`. Each sample has latent tube ",
    "coordinates \\((\\theta_i,u_i)\\), where\n\n",
    "$$\n",
    "u_i = \\|x_i\\| - r.\n",
    "$$\n\n",
    "The density-aware oracle lives on a structured periodic grid in \\((\\theta,u)\\). ",
    "For an oracle edge \\((a,b)\\), the baseline weight is\n\n",
    "$$\n",
    "w_{ab}=\\|z_a-z_b\\|\\,\\min\\{\\rho(m_{ab})^{-\\alpha}, c_{\\max}\\},\n",
    "$$\n\n",
    "where \\(m_{ab}\\) is the edge midpoint in radial coordinate, ",
    "\\(\\rho(u)=\\max\\{\\exp[-(u/\\sigma)^2],\\rho_{\\min}\\}\\), and ",
    "\\(c_{\\max}=\\infty\\) unless a penalty cap is specified.\n\n",
    "Primary truth for method comparison is the tube oracle with ",
    "\\(\\alpha=", params$primary_alpha, "\\), \\(\\rho_{\\min}=",
    params$primary_floor, "\\), no penalty cap, and `attach_extra = ",
    params$primary_attach_extra, "`.\n\n",
    "![Observed geometry and latent tube coordinates](", fig.link(fig.files$geometry), ")\n\n",
    "## 2. Oracle Sensitivity\n\n",
    "The first question is whether the oracle itself is stable enough to serve as truth. ",
    "The table below varies \\(\\alpha\\) while holding the density floor and attachment rule fixed.\n\n",
    md.table(short.table(
        alpha.summary[, c("alpha", "median", "q90", "max", "spearman_vs_arc",
                          "spearman_vs_euclidean")],
        c("alpha", "median", "q90", "max", "rho.arc", "rho.euc")
    )), "\n\n",
    "At \\(\\alpha=0.5\\), median oracle distance is ", fmt(alpha05$median),
    " and Spearman correlation with latent arc distance is ", fmt(alpha05$spearman_vs_arc),
    ". At \\(\\alpha=1\\), maximum oracle distance rises to ", fmt(alpha1$max),
    " and rank agreement with latent arc drops to ", fmt(alpha1$spearman_vs_arc),
    ". This repeats the earlier warning: \\(\\alpha=1\\) is tail dominated unless floors, caps, ",
    "or attachment rules are handled carefully.\n\n",
    "![Oracle sensitivity to alpha, floors/caps, and attachment](", fig.link(fig.files$oracle.sensitivity), ")\n\n",
    "### Floor And Cap Sweep\n\n",
    "The floor/cap sweep below focuses on \\(\\alpha=0.5\\) and \\(\\alpha=1\\). ",
    "Rows are scored against the primary \\(\\alpha=0.5\\) oracle, not against their own ",
    "variant, so high relative error for \\(\\alpha=1\\) means that the strong-density truth is ",
    "a different object.\n\n",
    md.table(short.table(
        floor.cap.summary[, c("alpha", "density_floor", "penalty_cap", "median",
                              "q90", "max", "spearman_vs_primary_truth",
                              "median_relative_error_to_primary")],
        c("alpha", "floor", "cap", "median", "q90", "max", "rho.primary", "med.err")
    )), "\n\n",
    "### Attachment Sweep\n\n",
    md.table(short.table(
        attachment.summary[, c("attach_extra", "median", "q90", "max",
                               "spearman_vs_primary_truth",
                               "median_relative_error_to_primary")],
        c("attach", "median", "q90", "max", "rho.primary", "med.err")
    )), "\n\n",
    "Attachment changes are modest for \\(\\alpha=0.5\\) in this case, but attachment remains ",
    "a required diagnostic because tail samples can otherwise inherit an artificial radial cost.\n\n",
    "## 3. Distance-Family Comparison\n\n",
    "The comparison includes Euclidean distance, latent arc distance, MST, CMST, iKNN, mKNN, ",
    "density-weighted iKNN, Fermat/PWSPD on iKNN, diffusion row distance, and a PHATE-like ",
    "log-potential distance. Diffusion and PHATE-like distances are included as support-aware ",
    "comparators, not as geodesic distances.\n\n",
    "![Distance-family comparison](", fig.link(fig.files$method.comparison), ")\n\n",
    md.table(short.table(
        method.summary[, c("distance", "finite_pair_fraction", "spearman_vs_primary_truth",
                           "median_relative_error_to_primary",
                           "q90_relative_error_to_primary",
                           "spearman_vs_arc", "spearman_vs_euclidean")],
        c("distance", "finite", "rho.primary", "med.err", "q90.err", "rho.arc", "rho.euc")
    )), "\n\n",
    "Best median relative error against the primary tube oracle is achieved by `",
    best.method$distance, "` with median error ", fmt(best.method$median_relative_error_to_primary),
    " and Spearman correlation ", fmt(best.method$spearman_vs_primary_truth), ". ",
    "This table should not be read as a final method ranking. It is a calibration check: the ",
    "primary tube oracle is still close enough to latent circle geometry that ordinary graph ",
    "geodesics remain competitive, while density-weighting changes the geometry in the intended ",
    "direction but needs better graph-topology and tail diagnostics.\n\n",
    "## 4. Edge Attribution\n\n",
    "For the primary density-weighted iKNN graph, each edge is attributed by Euclidean length, ",
    "mid-edge tube density, density penalty, density-weighted length, arc distance, primary ",
    "tube truth, and truth/Euclidean ratio.\n\n",
    "![Edge attribution](", fig.link(fig.files$edge.attribution), ")\n\n",
    "Top edges by tube-truth/Euclidean ratio:\n\n",
    md.table(short.table(
        top.edge.attr[, c("from", "to", "base_weight", "density", "penalty",
                         "weight", "arc", "primary_truth", "truth_to_euclidean")],
        c("from", "to", "ell", "rho", "penalty", "w", "arc", "truth", "truth/ell")
    )), "\n\n",
    "These edges are the first audit targets for false-shortcut and tail-dominance analysis. ",
    "Large truth/Euclidean ratios mean the observed endpoints are geometrically close but the ",
    "oracle regards travel between them as expensive, usually because the edge crosses radial ",
    "low-density support.\n\n",
    "## 5. Path Attribution\n\n",
    "Four diagnostic pairs were selected: an opposite-angle pair, a Euclidean-close pair with ",
    "high truth/Euclidean ratio, a near-angle pair with large radial separation, and the largest ",
    "base-graph underestimate relative to primary tube truth.\n\n",
    md.table(short.table(
        pair.table[, c("label", "from", "to", "arc", "euclidean", "primary_truth",
                       "base_graph_distance", "truth_to_euclidean", "truth_to_graph",
                       "u_from", "u_to")],
        c("label", "from", "to", "arc", "euc", "truth", "base.graph",
          "truth/euc", "truth/base", "u.from", "u.to")
    )), "\n\n",
    "![Path attribution](", fig.link(fig.files$path.attribution), ")\n\n",
    "Path summaries:\n\n",
    md.table(short.table(
        path.summary[, c("label", "from", "to", "primary_truth", "graph", "n_edges",
                         "base_length", "weighted_length", "min_edge_density",
                         "max_edge_penalty")],
        c("label", "from", "to", "truth", "graph", "edges", "base.len",
          "weighted.len", "min.rho", "max.penalty")
    )), "\n\n",
    "The path panels show where density-aware edge weights alter the cheapest route through the ",
    "sample graph. In the present graph, topology still matters more than edge weighting: if an ",
    "unsupported local edge is admitted, density weighting can penalize it only when the midpoint ",
    "density summary sees the low-density crossing.\n\n",
    "## 6. Interpretation\n\n",
    "1. \\(\\alpha=0.5\\) is a usable first tube-oracle truth. It remains strongly related to latent ",
    "arc geometry while adding a real radial-support penalty.\n\n",
    "2. \\(\\alpha=1\\) should be treated as a stress test. Without robust caps or stronger ",
    "attachment diagnostics, it is dominated by tail samples and produces extremely large maximum ",
    "distances.\n\n",
    "3. Density floors and penalty caps are not cosmetic. They define the tail behavior of the ",
    "metric. Any future report using density-aware distance must show floor/cap hit rates and ",
    "edge inflation distributions.\n\n",
    "4. Fermat/PWSPD is a necessary comparator. In this benchmark, rooted \\(p=2\\) Fermat on the ",
    "iKNN topology is a direct implicit-density alternative to explicit \\(\\rho^{-\\alpha}\\) ",
    "edge weighting.\n\n",
    "5. Diffusion and PHATE-like distances behave as support-aware comparators but should remain ",
    "separate from graph geodesic distances. Their low or high agreement with tube truth answers ",
    "a different question: whether points diffuse to similar neighborhoods, not whether the ",
    "cheapest supported path is short.\n\n",
    "## 7. Next Implementation Steps\n\n",
    "1. Add reusable density-transform and edge-attribution helpers so this report does not carry ",
    "one-off implementations.\n\n",
    "2. Extend `latent.tube.oracle()` to return edge density, base length, penalty, floor/cap hits, ",
    "and attachment diagnostics directly.\n\n",
    "3. Run a larger size/seed sweep with \\(\\alpha\\in\\{0,0.25,0.5,0.75,1\\}\\), but treat ",
    "\\(\\alpha=1\\) as a stress setting unless robust caps are enabled.\n\n",
    "4. Add nonuniform angular density, gaps, nearby non-touching arcs, and outliers before moving ",
    "to compositional microbiome-like data.\n"
)

writeLines(report, md.path)
cat("Wrote Markdown report: ", md.path, "\n", sep = "")
cat("Wrote figures under: ", fig.dir, "\n", sep = "")
cat("Wrote cache under: ", cache.dir, "\n", sep = "")

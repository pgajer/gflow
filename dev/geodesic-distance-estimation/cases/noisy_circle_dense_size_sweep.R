# Dense Noisy Circle Size Sweep: iKNN, mKNN, Pruning, Components

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
report.dir <- file.path(out.root, "reports/noisy-circle-dense-size-sweep")
fig.dir <- file.path(out.root, "figures/noisy-circle-dense-size-sweep")
cache.dir <- file.path(out.root, "cache/noisy-circle-dense-size-sweep")
run.cache.dir <- file.path(cache.dir, "runs")
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run.cache.dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
    radius = 1,
    geometry_noise = 0.08,
    geometry_noise_type = "normal",
    n_grid = seq(50L, 500L, by = 10L),
    seeds = 4101:4103,
    k_grid = 3L:10L,
    pruning_delta_grid = c(0, 0.05, 0.10),
    path_edge_ratio_percentile = 0.5,
    large_component_min_fraction = 0.05,
    large_component_min_size = 8L,
    shortcut_arc_threshold = pi / 4
)

make.noisy.circle.case <- function(n, seed, params) {
    X.df <- generate.circle.data(
        n = n,
        radius = params$radius,
        noise = params$geometry_noise,
        type = "random",
        noise.type = params$geometry_noise_type,
        seed = seed
    )
    X <- as.matrix(X.df[, c("x", "y")])
    theta <- as.double(X.df$angles)
    list(
        X = X,
        theta = theta,
        true.dist = circle.dist.matrix(theta, radius = params$radius),
        oracle = weighted.circle.graph(theta, radius = params$radius),
        mst = mst.graph.from.cmst(X)
    )
}

pruned.iknn.graph <- function(X, k, delta, path.percentile) {
    g <- create.iknn.graphs(
        X,
        kmin = k,
        kmax = k,
        max.path.edge.ratio.deviation.thld = delta,
        path.edge.ratio.percentile = path.percentile,
        threshold.percentile = 0,
        compute.full = TRUE,
        with.isize.pruning = FALSE,
        pca.dim = NULL,
        n.cores = 1L,
        verbose = FALSE
    )$geom_pruned_graphs[[1L]]
    sanitize.graph.weights(g, X = X, reweight = TRUE)
}

local.mknn.graph <- function(X, k) {
    n <- nrow(X)
    k <- min(as.integer(k), n - 1L)
    D <- euclidean.dist.matrix(X)
    diag(D) <- Inf
    knn <- t(apply(D, 1, function(x) order(x, seq_along(x))[seq_len(k)]))
    keep <- matrix(FALSE, n, n)
    for (i in seq_len(n)) {
        keep[i, knn[i, ]] <- TRUE
    }
    mutual <- keep & t(keep)
    ij <- which(upper.tri(mutual) & mutual, arr.ind = TRUE)
    if (!nrow(ij)) {
        return(graph.from.edges(n, data.frame(from = integer(), to = integer(), weight = numeric())))
    }
    graph.from.edges(n, data.frame(from = ij[, 1], to = ij[, 2], weight = D[ij]))
}

pruned.mknn.graph <- function(X, k, delta, path.percentile) {
    if (delta <= 0) {
        return(local.mknn.graph(X, k = k))
    }
    g <- create.mknn.graphs(
        X,
        kmin = k,
        kmax = k,
        max.path.edge.ratio.thld = delta,
        path.edge.ratio.percentile = path.percentile,
        compute.full = TRUE,
        pca.dim = NULL,
        verbose = FALSE
    )$pruned_graphs[[1L]]
    sanitize.graph.weights(g, X = X, reweight = TRUE)
}

component.info <- function(graph) {
    comp <- igraph::components(igraph.from.graph(graph$adj_list, graph$weight_list))
    sizes <- as.integer(comp$csize)
    list(
        membership = as.integer(comp$membership),
        sizes = sizes,
        largest_id = which.max(sizes),
        largest_size = max(sizes),
        n_components = comp$no
    )
}

induced.graph <- function(graph, vertices) {
    vertices <- sort(unique(as.integer(vertices)))
    n <- length(vertices)
    if (!n) return(graph.from.edges(0L, data.frame(from = integer(), to = integer(), weight = numeric())))
    edges <- edge.table(graph$adj_list, graph$weight_list)
    if (!nrow(edges)) return(graph.from.edges(n, data.frame(from = integer(), to = integer(), weight = numeric())))
    keep <- edges$from %in% vertices & edges$to %in% vertices
    edges <- edges[keep, , drop = FALSE]
    if (!nrow(edges)) return(graph.from.edges(n, data.frame(from = integer(), to = integer(), weight = numeric())))
    map <- seq_along(vertices)
    names(map) <- as.character(vertices)
    edges$from <- unname(map[as.character(edges$from)])
    edges$to <- unname(map[as.character(edges$to)])
    graph.from.edges(n, edges)
}

dsu.make <- function(keys) {
    parent <- seq_along(keys)
    names(parent) <- as.character(keys)
    find <- function(x) {
        idx <- match(as.character(x), names(parent))
        while (parent[[idx]] != idx) {
            parent[[idx]] <<- parent[[parent[[idx]]]]
            idx <- parent[[idx]]
        }
        idx
    }
    union <- function(a, b) {
        ra <- find(a)
        rb <- find(b)
        if (ra == rb) return(FALSE)
        parent[[rb]] <<- ra
        TRUE
    }
    list(find = find, union = union)
}

join.large.components <- function(base.graph, source.graph, min.size) {
    n <- length(base.graph$adj_list)
    comp <- component.info(base.graph)
    large.ids <- which(comp$sizes >= min.size)
    large.vertices <- which(comp$membership %in% large.ids)
    if (length(large.ids) <= 1L) {
        return(list(graph = base.graph, added_edges = 0L, large_vertices = large.vertices))
    }

    source.edges <- edge.table(source.graph$adj_list, source.graph$weight_list)
    if (!nrow(source.edges)) {
        return(list(graph = base.graph, added_edges = 0L, large_vertices = large.vertices))
    }
    source.edges$from_comp <- comp$membership[source.edges$from]
    source.edges$to_comp <- comp$membership[source.edges$to]
    keep <- source.edges$from_comp != source.edges$to_comp &
        source.edges$from_comp %in% large.ids &
        source.edges$to_comp %in% large.ids
    source.edges <- source.edges[keep, , drop = FALSE]
    if (!nrow(source.edges)) {
        return(list(graph = base.graph, added_edges = 0L, large_vertices = large.vertices))
    }
    source.edges <- source.edges[order(source.edges$weight), ]

    dsu <- dsu.make(large.ids)
    add <- vector("list", length(large.ids) - 1L)
    ptr <- 0L
    for (r in seq_len(nrow(source.edges))) {
        if (dsu$union(source.edges$from_comp[[r]], source.edges$to_comp[[r]])) {
            ptr <- ptr + 1L
            add[[ptr]] <- source.edges[r, c("from", "to", "weight")]
            if (ptr == length(large.ids) - 1L) break
        }
    }
    if (!ptr) {
        return(list(graph = base.graph, added_edges = 0L, large_vertices = large.vertices))
    }
    add.edges <- do.call(rbind, add[seq_len(ptr)])
    repaired <- union.graphs(base.graph, graph.from.edges(n, add.edges))
    list(graph = repaired, added_edges = ptr, large_vertices = large.vertices)
}

candidate.metrics <- function(graph, case, eval.vertices, params,
                              base.graph.type, component.policy,
                              k, pruning.delta, bridge.source = "none",
                              added.edges = 0L) {
    full.comp <- component.info(graph)
    eval.vertices <- sort(unique(as.integer(eval.vertices)))
    eval.graph <- induced.graph(graph, eval.vertices)
    eval.true <- case$true.dist[eval.vertices, eval.vertices, drop = FALSE]

    if (length(eval.vertices) < 2L) {
        geo <- data.frame(
            finite_pair_fraction = NA_real_, pearson = NA_real_, spearman = NA_real_,
            scale = NA_real_, median_relative_error = NA_real_, q90_relative_error = NA_real_,
            mean_neighbor_overlap = NA_real_
        )
        eval.summary <- data.frame(
            n_edges = 0L, n_components = NA_integer_, cycle_rank = NA_real_,
            mean_degree = NA_real_, min_degree = NA_real_, max_degree = NA_real_,
            edge_length_median = NA_real_, edge_length_max = NA_real_,
            false_shortcut_rate = NA_real_, max_edge_true_geodesic = NA_real_
        )
    } else {
        geo <- geodesic.metrics(eval.graph$adj_list, eval.graph$weight_list, eval.true)
        eval.summary <- graph.summary.metrics(
            eval.graph$adj_list,
            eval.graph$weight_list,
            true.dist = eval.true,
            shortcut.threshold = params$shortcut_arc_threshold
        )
    }

    row <- cbind(
        data.frame(
            base_graph = base.graph.type,
            component_policy = component.policy,
            bridge_source = bridge.source,
            k = k,
            pruning_delta = pruning.delta,
            added_bridge_edges = added.edges,
            n_eval = length(eval.vertices),
            eval_fraction = length(eval.vertices) / nrow(case$X),
            n_components_full = full.comp$n_components,
            largest_component_fraction_full = full.comp$largest_size / nrow(case$X)
        ),
        eval.summary,
        geo
    )
    row$score <- score.graph.metrics(row)
    row
}

evaluate.base.graph <- function(base.graph, case, params, base.graph.type, k,
                                pruning.delta, iknn.source = NULL) {
    n <- nrow(case$X)
    min.large <- max(params$large_component_min_size,
                     ceiling(params$large_component_min_fraction * n))
    comp <- component.info(base.graph)
    largest.vertices <- which(comp$membership == comp$largest_id)
    large.vertices <- which(comp$membership %in% which(comp$sizes >= min.large))

    rows <- list()
    rows[[length(rows) + 1L]] <- candidate.metrics(
        base.graph, case, seq_len(n), params, base.graph.type,
        component.policy = "all_components", k = k, pruning.delta = pruning.delta
    )
    rows[[length(rows) + 1L]] <- candidate.metrics(
        base.graph, case, largest.vertices, params, base.graph.type,
        component.policy = "major_component_only", k = k, pruning.delta = pruning.delta
    )

    if (identical(base.graph.type, "mKNN")) {
        mst.joined <- join.large.components(base.graph, case$mst, min.size = min.large)
        rows[[length(rows) + 1L]] <- candidate.metrics(
            mst.joined$graph, case, mst.joined$large_vertices, params,
            base.graph.type, component.policy = "join_large_components",
            k = k, pruning.delta = pruning.delta, bridge.source = "MST",
            added.edges = mst.joined$added_edges
        )
        if (!is.null(iknn.source)) {
            iknn.joined <- join.large.components(base.graph, iknn.source, min.size = min.large)
            rows[[length(rows) + 1L]] <- candidate.metrics(
                iknn.joined$graph, case, iknn.joined$large_vertices, params,
                base.graph.type, component.policy = "join_large_components",
                k = k, pruning.delta = pruning.delta, bridge.source = "iKNN",
                added.edges = iknn.joined$added_edges
            )
        }
    }

    do.call(rbind, rows)
}

safe.mean <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe.sd <- function(x) {
    if (sum(!is.na(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
}

summarize.metrics <- function(metrics) {
    keys <- unique(metrics[, c("n", "base_graph", "component_policy", "bridge_source", "k", "pruning_delta")])
    rows <- lapply(seq_len(nrow(keys)), function(i) {
        key <- keys[i, ]
        x <- metrics[
            metrics$n == key$n &
                metrics$base_graph == key$base_graph &
                metrics$component_policy == key$component_policy &
                metrics$bridge_source == key$bridge_source &
                metrics$k == key$k &
                metrics$pruning_delta == key$pruning_delta,
            , drop = FALSE
        ]
        data.frame(
            n = key$n,
            base_graph = key$base_graph,
            component_policy = key$component_policy,
            bridge_source = key$bridge_source,
            k = key$k,
            pruning_delta = key$pruning_delta,
            runs = nrow(x),
            score_mean = safe.mean(x$score),
            score_sd = safe.sd(x$score),
            median_relative_error_mean = safe.mean(x$median_relative_error),
            median_relative_error_sd = safe.sd(x$median_relative_error),
            q90_relative_error_mean = safe.mean(x$q90_relative_error),
            spearman_mean = safe.mean(x$spearman),
            finite_pair_fraction_mean = safe.mean(x$finite_pair_fraction),
            mean_neighbor_overlap_mean = safe.mean(x$mean_neighbor_overlap),
            eval_fraction_mean = safe.mean(x$eval_fraction),
            n_components_full_mean = safe.mean(x$n_components_full),
            largest_component_fraction_full_mean = safe.mean(x$largest_component_fraction_full),
            n_edges_mean = safe.mean(x$n_edges),
            added_bridge_edges_mean = safe.mean(x$added_bridge_edges)
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$n, out$score_mean, out$median_relative_error_mean), ]
}

method.label <- function(x) {
    paste(x$base_graph, x$component_policy, x$bridge_source, sep = " / ")
}

short.method.label <- function(base.graph, policy, bridge) {
    if (identical(policy, "all_components")) return(base.graph)
    if (identical(policy, "major_component_only")) return(paste(base.graph, "major", sep = "+"))
    paste(base.graph, "join", bridge, sep = "+")
}

grid <- expand.grid(
    n = params$n_grid,
    seed = params$seeds,
    KEEP.OUT.ATTRS = FALSE
)

run.cache.path <- function(n, seed) {
    file.path(run.cache.dir, sprintf("noisy_circle_dense_size_sweep_n%04d_seed%s.rds", n, seed))
}

run.single.case <- function(n.run, seed.run, params) {
    message(sprintf("Dense noisy circle sweep n=%s seed=%s", n.run, seed.run))
    case <- make.noisy.circle.case(n.run, seed.run, params)
    metric.rows <- list()
    row.ptr <- 0L
    for (k in params$k_grid) {
        for (delta in params$pruning_delta_grid) {
            iknn <- pruned.iknn.graph(case$X, k, delta, params$path_edge_ratio_percentile)
            mknn <- pruned.mknn.graph(case$X, k, delta, params$path_edge_ratio_percentile)

            iknn.metrics <- evaluate.base.graph(
                iknn, case, params, base.graph.type = "iKNN",
                k = k, pruning.delta = delta
            )
            mknn.metrics <- evaluate.base.graph(
                mknn, case, params, base.graph.type = "mKNN",
                k = k, pruning.delta = delta, iknn.source = iknn
            )

            for (block in list(iknn.metrics, mknn.metrics)) {
                block$n <- n.run
                block$seed <- seed.run
                block$large_component_min_size <- max(
                    params$large_component_min_size,
                    ceiling(params$large_component_min_fraction * n.run)
                )
                row.ptr <- row.ptr + 1L
                metric.rows[[row.ptr]] <- block
            }
        }
    }
    do.call(rbind, metric.rows)
}

metric.rows <- list()
for (g in seq_len(nrow(grid))) {
    n.run <- grid$n[[g]]
    seed.run <- grid$seed[[g]]
    path <- run.cache.path(n.run, seed.run)
    if (file.exists(path)) {
        message(sprintf("Reading cached dense sweep n=%s seed=%s", n.run, seed.run))
        metric.rows[[g]] <- readRDS(path)
    } else {
        metric.rows[[g]] <- run.single.case(n.run, seed.run, params)
        saveRDS(metric.rows[[g]], path)
    }
}

metrics <- do.call(rbind, metric.rows)
metrics$method <- mapply(
    short.method.label,
    metrics$base_graph,
    metrics$component_policy,
    metrics$bridge_source,
    USE.NAMES = FALSE
)
summary.metrics <- summarize.metrics(metrics)
summary.metrics$method <- mapply(
    short.method.label,
    summary.metrics$base_graph,
    summary.metrics$component_policy,
    summary.metrics$bridge_source,
    USE.NAMES = FALSE
)

best.by.n <- do.call(rbind, lapply(split(summary.metrics, summary.metrics$n), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[seq_len(min(12L, nrow(x))), , drop = FALSE]
}))

best.by.method.n <- do.call(rbind, lapply(split(summary.metrics, list(summary.metrics$n, summary.metrics$method), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
best.by.method.n <- best.by.method.n[order(best.by.method.n$n, best.by.method.n$score_mean), ]
best.by.n <- best.by.method.n

winner.by.n <- do.call(rbind, lapply(split(best.by.method.n, best.by.method.n$n), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
winner.by.n <- winner.by.n[order(winner.by.n$n), ]

raw.all.components <- summary.metrics[summary.metrics$component_policy == "all_components", ]
raw.base.best <- do.call(rbind, lapply(split(raw.all.components, list(raw.all.components$n, raw.all.components$base_graph), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
raw.base.best <- raw.base.best[order(raw.base.best$n, raw.base.best$base_graph), ]

best.pruning.by.base <- do.call(rbind, lapply(split(summary.metrics, list(summary.metrics$n, summary.metrics$base_graph, summary.metrics$pruning_delta), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
best.pruning.by.base <- best.pruning.by.base[order(best.pruning.by.base$n, best.pruning.by.base$base_graph, best.pruning.by.base$pruning_delta), ]

disconnected.mknn <- raw.all.components[
    raw.all.components$base_graph == "mKNN" &
        raw.all.components$n_components_full_mean > 1,
    , drop = FALSE
]
worst.disconnected.mknn <- disconnected.mknn[order(-disconnected.mknn$n_components_full_mean), ]
worst.disconnected.mknn <- worst.disconnected.mknn[seq_len(min(8L, nrow(worst.disconnected.mknn))), , drop = FALSE]

best.delta.curves <- do.call(rbind, lapply(split(summary.metrics, list(summary.metrics$n, summary.metrics$method, summary.metrics$pruning_delta), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
best.delta.curves <- best.delta.curves[order(best.delta.curves$n, best.delta.curves$method, best.delta.curves$pruning_delta), ]

method.win.counts <- as.data.frame(table(winner.by.n$method), stringsAsFactors = FALSE)
names(method.win.counts) <- c("method", "n_size_wins")
method.win.counts <- method.win.counts[order(-method.win.counts$n_size_wins, method.win.counts$method), ]

checkpoint.n <- params$n_grid[params$n_grid %% 50L == 0L]
winner.checkpoints <- winner.by.n[winner.by.n$n %in% checkpoint.n, ]
best.method.checkpoints <- best.by.method.n[best.by.method.n$n %in% checkpoint.n, ]

raw.by.n.k.base <- do.call(rbind, lapply(split(raw.all.components, list(raw.all.components$n, raw.all.components$base_graph, raw.all.components$k), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
raw.by.n.k.base <- raw.by.n.k.base[order(raw.by.n.k.base$n, raw.by.n.k.base$base_graph, raw.by.n.k.base$k), ]

raw.unpruned <- raw.all.components[raw.all.components$pruning_delta == 0, ]
raw.unpruned.best <- do.call(rbind, lapply(split(raw.unpruned, list(raw.unpruned$n, raw.unpruned$base_graph), drop = TRUE), function(x) {
    x <- x[order(x$score_mean, x$median_relative_error_mean), ]
    x[1L, , drop = FALSE]
}))
raw.pruned.best <- raw.base.best
raw.pruning.benefit <- merge(
    raw.pruned.best[, c("n", "base_graph", "median_relative_error_mean", "q90_relative_error_mean", "k", "pruning_delta")],
    raw.unpruned.best[, c("n", "base_graph", "median_relative_error_mean", "q90_relative_error_mean", "k")],
    by = c("n", "base_graph"),
    suffixes = c("_best", "_unpruned")
)
raw.pruning.benefit$median_error_delta <- raw.pruning.benefit$median_relative_error_mean_best -
    raw.pruning.benefit$median_relative_error_mean_unpruned
raw.pruning.benefit$q90_error_delta <- raw.pruning.benefit$q90_relative_error_mean_best -
    raw.pruning.benefit$q90_relative_error_mean_unpruned
raw.pruning.benefit <- raw.pruning.benefit[order(raw.pruning.benefit$n, raw.pruning.benefit$base_graph), ]

transition.summary <- do.call(rbind, lapply(split(raw.all.components, list(raw.all.components$n, raw.all.components$base_graph, raw.all.components$k), drop = TRUE), function(x) {
    x0 <- x[x$pruning_delta == 0, , drop = FALSE]
    if (!nrow(x0)) x0 <- x
    x0 <- x0[order(x0$score_mean, x0$median_relative_error_mean), ][1L, , drop = FALSE]
    x0
}))
transition.summary <- transition.summary[order(transition.summary$n, transition.summary$base_graph, transition.summary$k), ]

utils::write.csv(metrics, file.path(cache.dir, "noisy_circle_dense_size_sweep_metrics.csv"), row.names = FALSE)
utils::write.csv(summary.metrics, file.path(cache.dir, "noisy_circle_dense_size_sweep_summary.csv"), row.names = FALSE)

saveRDS(
    list(
        params = params,
        metrics = metrics,
        summary = summary.metrics,
        best.by.n = best.by.n,
        best.by.method.n = best.by.method.n,
        best.delta.curves = best.delta.curves,
        winner.by.n = winner.by.n,
        method.win.counts = method.win.counts,
        winner.checkpoints = winner.checkpoints,
        best.method.checkpoints = best.method.checkpoints,
        raw.base.best = raw.base.best,
        raw.by.n.k.base = raw.by.n.k.base,
        best.pruning.by.base = best.pruning.by.base,
        raw.pruning.benefit = raw.pruning.benefit,
        transition.summary = transition.summary,
        worst.disconnected.mknn = worst.disconnected.mknn
    ),
    file.path(cache.dir, "noisy_circle_dense_size_sweep_results.rds")
)

fig.files <- list()
method.levels <- unique(best.by.method.n$method)
method.cols <- setNames(grDevices::hcl.colors(length(method.levels), "Dark 3"), method.levels)

dense.matrix <- function(x, value.column, base.graph) {
    x <- x[x$base_graph == base.graph, ]
    out <- matrix(NA_real_, nrow = length(params$k_grid), ncol = length(params$n_grid),
                  dimnames = list(paste0("k=", params$k_grid), params$n_grid))
    for (i in seq_len(nrow(x))) {
        out[match(x$k[[i]], params$k_grid), match(x$n[[i]], params$n_grid)] <- x[[value.column]][[i]]
    }
    out
}

draw.heatmap <- function(mat, main, zlim = range(mat, na.rm = TRUE), palette.name = "YlOrRd") {
    cols <- grDevices::hcl.colors(128, palette.name, rev = FALSE)
    image(
        x = seq_len(ncol(mat)),
        y = seq_len(nrow(mat)),
        z = t(mat),
        col = cols,
        zlim = zlim,
        axes = FALSE,
        xlab = "sample size n",
        ylab = "k",
        main = main
    )
    axis(1, at = seq(1, ncol(mat), by = 5), labels = colnames(mat)[seq(1, ncol(mat), by = 5)], las = 2, cex.axis = 0.75)
    axis(2, at = seq_len(nrow(mat)), labels = params$k_grid, las = 1)
    box()
}

fig.files$best.error <- file.path(fig.dir, "01_best_error_by_size.png")
save.png(fig.files$best.error, {
    layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1.1))
    oldpar <- par(mar = c(4.8, 4.8, 3.2, 1.2), xpd = FALSE)
    on.exit(par(oldpar), add = TRUE)
    plot(NA, xlim = range(params$n_grid), ylim = range(best.by.method.n$median_relative_error_mean, na.rm = TRUE),
         xlab = "sample size n", ylab = "best mean median relative error",
         main = "Best GDA by Method Family")
    for (method in method.levels) {
        x <- best.by.method.n[best.by.method.n$method == method, ]
        lines(x$n, x$median_relative_error_mean, type = "b", pch = 19, col = method.cols[[method]])
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = method.levels,
           col = method.cols[method.levels], lty = 1, pch = 19,
           bty = "n", cex = 0.9, ncol = 3)
}, width = 1700, height = 950)

fig.files$heatmaps <- file.path(fig.dir, "02_raw_error_heatmaps.png")
save.png(fig.files$heatmaps, {
    oldpar <- par(mfrow = c(2, 2), mar = c(5.8, 4.5, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    iknn.err <- dense.matrix(raw.by.n.k.base, "median_relative_error_mean", "iKNN")
    mknn.err <- dense.matrix(raw.by.n.k.base, "median_relative_error_mean", "mKNN")
    iknn.q90 <- dense.matrix(raw.by.n.k.base, "q90_relative_error_mean", "iKNN")
    mknn.q90 <- dense.matrix(raw.by.n.k.base, "q90_relative_error_mean", "mKNN")
    err.zlim <- range(iknn.err, mknn.err, na.rm = TRUE)
    q90.zlim <- range(iknn.q90, mknn.q90, na.rm = TRUE)
    draw.heatmap(iknn.err, "iKNN: best median relative error by n,k", zlim = err.zlim)
    draw.heatmap(mknn.err, "mKNN: best median relative error by n,k", zlim = err.zlim)
    draw.heatmap(iknn.q90, "iKNN: best q90 relative error by n,k", zlim = q90.zlim)
    draw.heatmap(mknn.q90, "mKNN: best q90 relative error by n,k", zlim = q90.zlim)
}, width = 1800, height = 1300)

fig.files$pruning <- file.path(fig.dir, "03_pruning_benefit_by_size.png")
save.png(fig.files$pruning, {
    oldpar <- par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.2, 1.2))
    on.exit(par(oldpar), add = TRUE)
    base.cols <- c(iKNN = "#1f77b4", mKNN = "#d95f02")
    ylim <- range(raw.pruning.benefit$median_error_delta, raw.pruning.benefit$q90_error_delta, na.rm = TRUE)
    for (metric in c("median_error_delta", "q90_error_delta")) {
        plot(NA, xlim = range(params$n_grid), ylim = ylim,
             xlab = "sample size n",
             ylab = "best pruned error - best unpruned error",
             main = if (metric == "median_error_delta") "Median error pruning benefit" else "Q90 error pruning benefit")
        abline(h = 0, col = "gray55", lty = 2)
        for (base in names(base.cols)) {
            x <- raw.pruning.benefit[raw.pruning.benefit$base_graph == base, ]
            lines(x$n, x[[metric]], type = "b", pch = 19, cex = 0.65, col = base.cols[[base]])
        }
        legend("topright", legend = names(base.cols), col = base.cols, lty = 1, pch = 19, bty = "n")
    }
}, width = 1700, height = 800)

fig.files$components <- file.path(fig.dir, "04_component_policy_diagnostics.png")
save.png(fig.files$components, {
    oldpar <- par(mfrow = c(2, 2), mar = c(5.8, 4.8, 3.2, 1.2))
    on.exit(par(oldpar), add = TRUE)
    plot(best.by.method.n$eval_fraction_mean, best.by.method.n$median_relative_error_mean,
         pch = 19, col = method.cols[best.by.method.n$method],
         xlab = "mean evaluated vertex fraction", ylab = "best mean median relative error",
         main = "CEEP filtering tradeoff")
    plot(best.by.method.n$n_components_full_mean, best.by.method.n$finite_pair_fraction_mean,
         pch = 19, col = method.cols[best.by.method.n$method],
         xlab = "mean full graph component count", ylab = "mean finite pair fraction",
         main = "Connectivity vs GDA support")
    legend("bottomleft", legend = method.levels, col = method.cols[method.levels],
           pch = 19, bty = "n", cex = 0.72)
    mknn.comp <- dense.matrix(transition.summary, "n_components_full_mean", "mKNN")
    mknn.finite <- dense.matrix(transition.summary, "finite_pair_fraction_mean", "mKNN")
    draw.heatmap(mknn.comp, "unpruned mKNN: mean component count", palette.name = "Reds")
    draw.heatmap(mknn.finite, "unpruned mKNN: finite-pair fraction", zlim = c(0, 1), palette.name = "YlGnBu")
}, width = 1800, height = 1300)

fig.files$winner <- file.path(fig.dir, "05_winner_map.png")
save.png(fig.files$winner, {
    oldpar <- par(mfrow = c(3, 1), mar = c(4.8, 4.8, 3.0, 1.2))
    on.exit(par(oldpar), add = TRUE)
    y <- match(winner.by.n$method, method.levels)
    plot(winner.by.n$n, y, pch = 15, cex = 1.5, col = method.cols[winner.by.n$method],
         yaxt = "n", ylim = c(0.5, length(method.levels) + 0.5),
         xlab = "sample size n", ylab = "winning method",
         main = "Winning method by sample size")
    axis(2, at = seq_along(method.levels), labels = method.levels, las = 1, cex.axis = 0.75)
    grid(col = "gray90")
    plot(winner.by.n$n, winner.by.n$k, type = "b", pch = 19,
         xlab = "sample size n", ylab = "selected k",
         main = "Winner selected k")
    plot(winner.by.n$n, winner.by.n$pruning_delta, type = "b", pch = 19,
         xlab = "sample size n", ylab = "selected pruning delta",
         main = "Winner selected pruning delta")
}, width = 1600, height = 1300)

rel.fig <- function(path) {
    file.path("../../figures/noisy-circle-dense-size-sweep", basename(path))
}

commands <- paste(
    "X.df <- gflow::generate.circle.data(",
    "  n = n, radius = 1, noise = 0.08,",
    "  type = \"random\", noise.type = \"normal\", seed = seed",
    ")",
    "X <- as.matrix(X.df[, c(\"x\", \"y\")])",
    "iknn <- gflow::create.iknn.graphs(",
    "  X, kmin = k, kmax = k,",
    "  max.path.edge.ratio.deviation.thld = delta,",
    "  path.edge.ratio.percentile = 0.5,",
    "  threshold.percentile = 0, compute.full = TRUE,",
    "  pca.dim = NULL, n.cores = 1L, verbose = FALSE",
    ")$geom_pruned_graphs[[1L]]",
    "mknn <- if (delta == 0) {",
    "  local.mknn.graph(X, k = k)",
    "} else {",
    "  gflow::create.mknn.graphs(",
    "    X, kmin = k, kmax = k,",
    "    max.path.edge.ratio.thld = delta,",
    "    path.edge.ratio.percentile = 0.5,",
    "    compute.full = TRUE, pca.dim = NULL, verbose = FALSE",
    "  )$pruned_graphs[[1L]]",
    "}",
    sep = "\n"
)

dense.winner.text <- paste(sprintf(
    "%s wins %s of %s sample sizes",
    method.win.counts$method,
    method.win.counts$n_size_wins,
    length(params$n_grid)
), collapse = "; ")

best.overall.method <- method.win.counts$method[[1L]]
best.overall.count <- method.win.counts$n_size_wins[[1L]]

report.path <- file.path(report.dir, "noisy_circle_dense_size_sweep.html")
html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Noisy Circle Dense Size Sweep</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;max-width:1240px;}h1,h2,h3{line-height:1.15;}code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f6f8fb;}figure{margin:24px 0;}img{max-width:100%;border:1px solid #d9dee7;}figcaption{font-size:13px;color:#52616f;margin-top:6px;}details{border:1px solid #d9dee7;padding:10px 14px;margin:18px 0;background:#fbfcfe;}summary{font-weight:600;cursor:pointer;}</style>",
    "</head><body>",
    "<h1>Noisy Circle Dense Size Sweep</h1>",
    sprintf("<p>Generated at %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "<h2>Objective</h2>",
    "<p>This dense sweep asks how iKNN, mKNN, geometric pruning, and CEEP-style component policies behave as sample size increases in the simplest closed-manifold benchmark: a noisy circle with known circular geodesic truth.</p>",
    "<p>The true geodesic distance between samples is the shorter arc length on the latent circle, computed from the generating latent angle before observational noise is added to the embedded coordinates. Graph shortest-path distances are compared to this truth after an optimal scalar rescaling.</p>",
    "<h2>Grid</h2>",
    sprintf("<p><code>n</code>: %s to %s by %s; seeds: %s; <code>k</code>: %s; pruning deltas: %s; large component threshold: <code>max(%s, ceiling(%s * n))</code>.</p>",
            min(params$n_grid), max(params$n_grid), unique(diff(params$n_grid)),
            paste(params$seeds, collapse = ", "),
            paste(params$k_grid, collapse = ", "), paste(params$pruning_delta_grid, collapse = ", "),
            params$large_component_min_size, params$large_component_min_fraction),
    "<p>The complete run evaluates ", format(nrow(metrics), big.mark = ","), " run-level candidate rows and ",
    format(nrow(summary.metrics), big.mark = ","), " averaged candidate rows. Complete tables are written to the cache directory as CSV and RDS files.</p>",
    "<h2>Construction Commands</h2><pre><code>", html.escape(commands), "</code></pre>",
    "<p><strong>Implementation note:</strong> for unpruned mKNN (<code>delta = 0</code>), this dense sweep uses the local deterministic mutual-kNN constructor defined in the script. It implements the same mutual-neighbor rule from the Euclidean distance matrix and avoids a low-level <code>create.mknn.graph()</code> segfault encountered during the long sweep.</p>",
    "<h2>Candidate Families</h2>",
    "<ul><li><strong>iKNN</strong>: inclusive k-nearest-neighbor graph, optionally geometrically pruned.</li><li><strong>iKNN+major</strong>: iKNN evaluated only on the largest connected component.</li><li><strong>mKNN</strong>: mutual k-nearest-neighbor graph, optionally geometrically pruned.</li><li><strong>mKNN+major</strong>: mKNN evaluated only on the largest connected component.</li><li><strong>mKNN+join+MST</strong>: large mKNN components joined using sparse MST-supported bridge edges.</li><li><strong>mKNN+join+iKNN</strong>: large mKNN components joined using sparse iKNN-supported bridge edges.</li></ul>",
    "<h2>Main Figures</h2>",
    sprintf("<figure><img src=\"%s\"><figcaption>Best mean median relative geodesic error by sample size, after optimizing over k and pruning delta inside each method family.</figcaption></figure>", rel.fig(fig.files$best.error)),
    sprintf("<figure><img src=\"%s\"><figcaption>Raw iKNN and mKNN heatmaps over sample size and k. Each cell uses the best pruning delta for that n, k, and base graph.</figcaption></figure>", rel.fig(fig.files$heatmaps)),
    sprintf("<figure><img src=\"%s\"><figcaption>Pruning benefit relative to the best unpruned raw graph. Values below zero mean pruning improved geodesic recovery.</figcaption></figure>", rel.fig(fig.files$pruning)),
    sprintf("<figure><img src=\"%s\"><figcaption>Component-policy diagnostics, including the unpruned mKNN connectivity transition over n and k.</figcaption></figure>", rel.fig(fig.files$components)),
    sprintf("<figure><img src=\"%s\"><figcaption>Winning method, selected k, and selected pruning delta as functions of sample size.</figcaption></figure>", rel.fig(fig.files$winner)),
    "<h2>Initial Automated Analysis</h2>",
    "<h3>Winner Counts</h3>",
    "<p>The current dense sweep identifies <strong>", html.escape(best.overall.method), "</strong> as the most frequent winner, with ",
    best.overall.count, " wins across ", length(params$n_grid), " sample sizes. The full count summary is: ",
    html.escape(dense.winner.text), ".</p>",
    html.table(method.win.counts),
    "<h3>What This Sweep Is Designed To Reveal</h3>",
    "<p>The important question is not only which method wins at a few hand-picked sizes. The dense grid should reveal whether winners change smoothly with sample size, whether optimal k drifts upward with n, whether pruning has a stable useful range, and whether mKNN failures are mostly connectivity failures or genuine geodesic distortions after connectivity is repaired.</p>",
    "<p>Because only three random samplings are used at each sample size, narrow local reversals should be treated as signals for follow-up rather than final claims. Stable bands across adjacent sample sizes are more meaningful than isolated one-size wins.</p>",
    "<h3>How To Read The Heatmaps</h3>",
    "<p>The raw iKNN/mKNN heatmaps show the base geometric behavior before CEEP component filtering or repair. A good region is not merely dark at one sample size; it should form a coherent band over neighboring values of n and k. A jagged isolated optimum is more likely to be seed noise or scoring sensitivity.</p>",
    "<p>The mKNN connectivity heatmap is especially important. If mKNN performs poorly where component counts are high and finite-pair fractions are low, then its main issue is disconnectedness. If it performs poorly even after connectivity has stabilized, then the mutual-neighbor rule is removing edges that are geometrically important for circular geodesic reconstruction.</p>",
    "<h3>Current Reporting Policy</h3>",
    "<p>This HTML report intentionally avoids printing the complete candidate table. The full table is too large for human inspection and is stored instead as <code>noisy_circle_dense_size_sweep_summary.csv</code> and <code>noisy_circle_dense_size_sweep_metrics.csv</code> in the cache directory. The tables below are human-facing summaries.</p>",
    "<h2>Winner By Sample Size</h2>",
    html.table(winner.by.n[, c("n", "method", "base_graph", "component_policy", "bridge_source", "k",
                               "pruning_delta", "runs", "score_mean", "median_relative_error_mean",
                               "q90_relative_error_mean", "spearman_mean", "finite_pair_fraction_mean",
                               "eval_fraction_mean", "n_components_full_mean", "largest_component_fraction_full_mean",
                               "n_edges_mean", "added_bridge_edges_mean")]),
    "<h2>Checkpoint Detail</h2>",
    "<p>These checkpoint rows show the best candidate for each method family every 50 samples. They are meant for sanity checking the trends visible in the figures.</p>",
    html.table(best.method.checkpoints[, c("n", "method", "base_graph", "component_policy", "bridge_source", "k",
                                           "pruning_delta", "runs", "score_mean", "median_relative_error_mean",
                                           "q90_relative_error_mean", "spearman_mean", "finite_pair_fraction_mean",
                                           "eval_fraction_mean", "n_components_full_mean", "n_edges_mean")]),
    "<h2>Raw Base-Graph Checkpoints</h2>",
    "<p>These rows compare the best raw all-component iKNN and mKNN candidates at the same 50-sample checkpoints. This is the cleanest view of the base graph constructors before component filtering or joining policies enter.</p>",
    html.table(raw.base.best[raw.base.best$n %in% checkpoint.n,
                             c("n", "base_graph", "k", "pruning_delta", "score_mean",
                               "median_relative_error_mean", "q90_relative_error_mean", "spearman_mean",
                               "finite_pair_fraction_mean", "n_components_full_mean",
                               "largest_component_fraction_full_mean", "n_edges_mean")]),
    "<h2>Pruning Benefit Checkpoints</h2>",
    "<p>Negative deltas mean that the best pruned raw graph improved over the best unpruned raw graph for the same base constructor.</p>",
    html.table(raw.pruning.benefit[raw.pruning.benefit$n %in% checkpoint.n,
                                   c("n", "base_graph", "k_best", "pruning_delta",
                                     "median_relative_error_mean_best", "median_relative_error_mean_unpruned",
                                     "median_error_delta", "q90_error_delta")]),
    "<details><summary>Worst disconnected raw mKNN settings</summary>",
    html.table(worst.disconnected.mknn[, c("n", "k", "pruning_delta", "median_relative_error_mean",
                                           "finite_pair_fraction_mean", "n_components_full_mean",
                                           "largest_component_fraction_full_mean", "n_edges_mean")]),
    "</details>",
    "</body></html>"
)
writeLines(html, report.path)

cat("Wrote report:\n", report.path, "\n", sep = "")
cat("Wrote figures under:\n", fig.dir, "\n", sep = "")
cat("Wrote cache under:\n", cache.dir, "\n", sep = "")

if (FALSE) {

winner.text <- paste(vapply(seq_len(nrow(winner.by.n)), function(i) {
    sprintf(
        "n = %s: %s, k = %s, pruning delta = %s, mean median relative error = %.4f, mean q90 relative error = %.4f",
        winner.by.n$n[[i]],
        winner.by.n$method[[i]],
        winner.by.n$k[[i]],
        winner.by.n$pruning_delta[[i]],
        winner.by.n$median_relative_error_mean[[i]],
        winner.by.n$q90_relative_error_mean[[i]]
    )
}, character(1)), collapse = "; ")

report.path <- file.path(report.dir, "noisy_circle_gda_ceep_iknn_mknn.html")
html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Noisy Circle GDA/CEEP iKNN-mKNN Benchmark</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;max-width:1220px;}h1,h2,h3{line-height:1.15;}code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f6f8fb;}figure{margin:24px 0;}img{max-width:100%;border:1px solid #d9dee7;}figcaption{font-size:13px;color:#52616f;margin-top:6px;}details{border:1px solid #d9dee7;padding:10px 14px;margin:18px 0;background:#fbfcfe;}summary{font-weight:600;cursor:pointer;}</style>",
    "</head><body>",
    "<h1>Noisy Circle GDA/CEEP iKNN-mKNN Benchmark</h1>",
    sprintf("<p>Generated at %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "<h2>Objective</h2>",
    "<p>This report treats geodesic distance approximation (GDA) as the primary metric-recovery problem, while allowing CEEP-motivated component policies such as major-component filtering or large-component repair.</p>",
    "<h2>Grid</h2>",
    sprintf("<p><code>n</code>: %s; seeds: %s; <code>k</code>: %s; pruning deltas: %s; large component threshold: <code>max(%s, ceiling(%s * n))</code>.</p>",
            paste(params$n_grid, collapse = ", "), paste(params$seeds, collapse = ", "),
            paste(params$k_grid, collapse = ", "), paste(params$pruning_delta_grid, collapse = ", "),
            params$large_component_min_size, params$large_component_min_fraction),
    "<h2>Construction Commands</h2><pre><code>", html.escape(commands), "</code></pre>",
    "<h2>Component Policies</h2>",
    "<ul><li><strong>all_components</strong>: evaluate the graph as constructed.</li><li><strong>major_component_only</strong>: evaluate only the largest connected component.</li><li><strong>join_large_components / MST</strong>: add minimal MST-supported bridges among large mKNN components.</li><li><strong>join_large_components / iKNN</strong>: add minimal iKNN-supported bridges among large mKNN components.</li></ul>",
    "<h2>Main Results</h2>",
    sprintf("<figure><img src=\"%s\"><figcaption>Best mean median relative geodesic error by sample size, after optimizing over k and pruning delta inside each method family.</figcaption></figure>", rel.fig(fig.files$best.error)),
    sprintf("<figure><img src=\"%s\"><figcaption>Effect of geometric pruning delta after optimizing over k for each method family.</figcaption></figure>", rel.fig(fig.files$pruning)),
    sprintf("<figure><img src=\"%s\"><figcaption>Component policy diagnostics: filtering/repair tradeoff and connectivity support.</figcaption></figure>", rel.fig(fig.files$components)),
    sprintf("<figure><img src=\"%s\"><figcaption>Run-to-run error distributions for each method family's best candidate at each sample size.</figcaption></figure>", rel.fig(fig.files$box)),
    "<h2>Analysis Of Results</h2>",
    "<h3>Can We Identify A Winner?</h3>",
    "<p>Within this pilot grid, there is no single method that cleanly wins at every sample size. The size-specific winners are: ",
    html.escape(winner.text), ".</p>",
    "<p>The most defensible current conclusion is that <strong>mildly pruned iKNN is the best scalable strategy in this noisy-circle family</strong>. It wins at <code>n = 160</code> and <code>n = 320</code>, and at <code>n = 80</code> it is only slightly behind mKNN. The small-sample mKNN win is real in this pilot, but it does not persist as sample size increases.</p>",
    "<p>If we force a single strategy recommendation from this data, it would be <strong>iKNN with mild geometric pruning</strong>, using <code>k</code> in the upper part of the tested grid and a positive pruning delta. In this run, the selected iKNN settings are <code>k = 8, delta = 0.05</code> for <code>n = 160</code> and <code>k = 10, delta = 0.10</code> for <code>n = 320</code>. This should be treated as a pilot conclusion, not a universal parameter rule.</p>",
    "<h3>iKNN Versus mKNN</h3>",
    "<p>Raw iKNN and raw mKNN are both strong local-geometry estimators on the noisy circle when <code>k</code> is sufficiently large. At <code>n = 80</code>, mKNN with <code>k = 10</code> and no pruning has the best mean median relative error. At <code>n = 160</code> and <code>n = 320</code>, iKNN is better. The mKNN result at <code>n = 320</code> is notably worse even after optimizing over the tested pruning grid, which suggests that mKNN's conservative edge rule can under-represent circular geometry at this sample size/grid range.</p>",
    html.table(raw.base.best[, c("n", "base_graph", "method", "k", "pruning_delta",
                                 "score_mean", "median_relative_error_mean",
                                 "q90_relative_error_mean", "spearman_mean",
                                 "finite_pair_fraction_mean", "n_components_full_mean",
                                 "largest_component_fraction_full_mean", "n_edges_mean")]),
    "<h3>Effect Of Geometric Pruning</h3>",
    "<p>Geometric pruning helps iKNN in this pilot. The best iKNN settings use positive pruning deltas at every tested size: <code>0.05</code> at <code>n = 80</code> and <code>n = 160</code>, and <code>0.10</code> at <code>n = 320</code>. The improvement is not enormous, but it is consistent enough to matter.</p>",
    "<p>For mKNN, pruning is not clearly beneficial in the same way. At <code>n = 80</code> and <code>n = 160</code>, the best mKNN-family setting is unpruned. At <code>n = 320</code>, the best mKNN-family setting uses <code>delta = 0.10</code>, but its error remains much worse than iKNN. This suggests that pruning alone does not solve the mKNN geometry issue in the larger noisy-circle case.</p>",
    html.table(best.pruning.by.base[, c("n", "base_graph", "method", "k", "pruning_delta",
                                        "score_mean", "median_relative_error_mean",
                                        "q90_relative_error_mean", "spearman_mean",
                                        "finite_pair_fraction_mean", "n_components_full_mean")]),
    "<h3>Component Policies And CEEP Perspective</h3>",
    "<p>The component-policy results are useful, but mostly as diagnostics. In the best mKNN settings, the graph is either connected or has only tiny extra components. Therefore <code>mKNN</code>, <code>mKNN+major</code>, <code>mKNN+join+MST</code>, and <code>mKNN+join+iKNN</code> often have identical geodesic-error summaries. In those cases the joining rules do not add meaningful bridge edges; they mainly change the evaluation scope or finite-pair coverage bookkeeping.</p>",
    "<p>At low <code>k</code>, mKNN can fragment severely. For example, the worst disconnected mKNN rows include cases with many components and a small largest component fraction. Major-component filtering can make CEEP-style evaluation possible, but it may discard too much of the circle to be a satisfactory GDA solution. Joining large components improves finite-pair coverage, but in the worst low-<code>k</code> cases it does not make the method competitive with tuned iKNN.</p>",
    html.table(worst.disconnected.mknn[, c("n", "k", "pruning_delta",
                                           "median_relative_error_mean",
                                           "finite_pair_fraction_mean",
                                           "n_components_full_mean",
                                           "largest_component_fraction_full_mean",
                                           "n_edges_mean")]),
    "<h3>Reading The Figures</h3>",
    "<p>The best-error figure shows the main conclusion visually: mKNN is slightly better at <code>n = 80</code>, but iKNN dominates by <code>n = 160</code> and remains strong at <code>n = 320</code>. The pruning curves show that positive pruning is helpful for iKNN, while mKNN is more sensitive and less reliable. The component-policy diagnostic shows that component filtering is mainly a CEEP tool, not a complete GDA fix. The boxplots show that the apparent iKNN advantage at larger sizes is not driven by a single seed; it is stable across the three random samplings in this pilot.</p>",
    "<h3>Current Recommendation</h3>",
    "<p>For the next round of noisy-circle GDA experiments, promote <strong>iKNN with mild geometric pruning</strong> as the primary baseline to beat. Keep mKNN in the benchmark because it can win at smaller sample size and because its conservative topology is scientifically attractive, but focus mKNN development on avoiding fragmentation without losing local geometry. The component-repair layer should remain in the framework, but this pilot does not show it as the winning mechanism for noisy circles.</p>",
    "<h2>Best Candidates By Sample Size</h2>",
    html.table(best.by.n[, c("n", "base_graph", "component_policy", "bridge_source", "k", "pruning_delta", "runs",
                             "score_mean", "median_relative_error_mean", "q90_relative_error_mean",
                             "spearman_mean", "finite_pair_fraction_mean", "eval_fraction_mean",
                             "n_components_full_mean", "largest_component_fraction_full_mean",
                             "n_edges_mean", "added_bridge_edges_mean")]),
    "<h2>Interpretation</h2>",
    "<p>The benchmark separates raw GDA from CEEP-relevant evaluation. A graph may be poor as a global metric if it is disconnected, but useful for CEEP if the major component has strong geodesic recovery and the discarded vertices are small low-support components.</p>",
    "<p>The mKNN repair policies test whether mKNN's weakness is mostly disconnectedness. MST bridges provide guaranteed sparse repair; iKNN bridges ask whether local-neighbor evidence can repair mKNN components more geometrically.</p>",
    "<details><summary>Best candidate per method family and sample size</summary>",
    html.table(best.by.method.n[, c("n", "method", "base_graph", "component_policy", "bridge_source", "k",
                                    "pruning_delta", "runs", "score_mean", "median_relative_error_mean",
                                    "q90_relative_error_mean", "spearman_mean", "finite_pair_fraction_mean",
                                    "eval_fraction_mean", "n_components_full_mean", "n_edges_mean")]),
    "</details>",
    "<details><summary>Complete averaged candidate summary</summary>",
    html.table(summary.metrics[, c("n", "method", "base_graph", "component_policy", "bridge_source", "k",
                                   "pruning_delta", "runs", "score_mean", "score_sd",
                                   "median_relative_error_mean", "median_relative_error_sd",
                                   "q90_relative_error_mean", "spearman_mean", "finite_pair_fraction_mean",
                                   "eval_fraction_mean", "n_components_full_mean",
                                   "largest_component_fraction_full_mean", "n_edges_mean",
                                   "added_bridge_edges_mean")]),
    "</details>",
    "</body></html>"
)
writeLines(html, report.path)

cat("Wrote report:\n", report.path, "\n", sep = "")
cat("Wrote figures under:\n", fig.dir, "\n", sep = "")
cat("Wrote cache under:\n", cache.dir, "\n", sep = "")
}

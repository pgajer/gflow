# Noisy Circle Strategy Comparison for Geodesic Distance Estimation
#
# Focus: graph geodesic recovery only. The response variable is intentionally
# ignored in this report.

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
report.dir <- file.path(out.root, "reports/noisy-circle-strategy-comparison")
fig.dir <- file.path(out.root, "figures/noisy-circle-strategy-comparison")
cache.dir <- file.path(out.root, "cache/noisy-circle-strategy-comparison")
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
    seed = 1002L,
    n = 160L,
    radius = 1,
    geometry_noise = 0.08,
    geometry_noise_type = "normal",
    q_grid = c(0.9, 0.95, 0.99, 0.9999),
    prune_alt_path_ratio = c(1.05, 1.15, 1.3),
    adaptive_scale_source = c("mst", "knn"),
    adaptive_multiplier = c(1.25, 1.5, 2, 2.5),
    adaptive_k_scale = c(3L, 5L, 8L),
    graph_k_grid = 2:10,
    cycle_base = c("mst", "cmst_q0.5", "cmst_q0.7"),
    cycle_n_add = c(1L, 3L, 5L, 10L),
    cycle_local_multiplier = c(1.5, 2, 2.5),
    cycle_min_ratio = c(5, 10, 20),
    shortcut_arc_threshold = pi / 4,
    replicate_n = c(80L, 160L, 320L),
    replicate_seeds = 3001:3005
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
        oracle.path = weighted.angular.path.graph(theta, radius = params$radius),
        mst = mst.graph.from.cmst(X)
    )
}

new.candidate.list <- function() {
    out <- list()
    add <- function(family, strategy, parameter, graph) {
        graph <- sanitize.graph.weights(graph)
        out[[length(out) + 1L]] <<- list(
            family = family,
            strategy = strategy,
            parameter = parameter,
            graph = graph
        )
        invisible(out)
    }
    env <- environment(add)
    list(add = add, get = function() env$out)
}

make.cmst.graph <- function(X, q.thld) {
    g <- create.cmst.graph(X, q.thld = q.thld, pca.dim = NULL, verbose = FALSE)
    list(adj_list = g$cmst_adj_list, weight_list = g$cmst_weight_list)
}

make.cycle.base <- function(name, case) {
    if (identical(name, "mst")) return(case$mst)
    if (identical(name, "cmst_q0.5")) return(make.cmst.graph(case$X, q.thld = 0.5))
    if (identical(name, "cmst_q0.7")) return(make.cmst.graph(case$X, q.thld = 0.7))
    stop("Unknown cycle base: ", name)
}

build.exhaustive.candidates <- function(case, params) {
    candidates <- new.candidate.list()
    candidates$add("reference", "MST only", NA_character_, case$mst)
    candidates$add("reference", "oracle weighted path", NA_character_, case$oracle.path)
    candidates$add("reference", "oracle weighted circle", NA_character_, case$oracle)

    for (q in params$q_grid) {
        base.graph <- make.cmst.graph(case$X, q.thld = q)
        candidates$add("cmst", "CMST", sprintf("q=%s", q), base.graph)
        for (ratio in params$prune_alt_path_ratio) {
            candidates$add(
                "high_q_pruned",
                "high-q CMST + shortcut pruning",
                sprintf("q=%s; alt_ratio=%s", q, ratio),
                prune.redundant.edges(base.graph, alt.path.ratio = ratio)
            )
        }
    }

    for (source.name in params$adaptive_scale_source) {
        for (mult in params$adaptive_multiplier) {
            for (k.scale in params$adaptive_k_scale) {
                candidates$add(
                    "local_adaptive",
                    "MST + local adaptive threshold",
                    sprintf("scale=%s; c=%s; k_scale=%s", source.name, mult, k.scale),
                    adaptive.threshold.graph(
                        case$X,
                        mst.graph = case$mst,
                        scale.source = source.name,
                        multiplier = mult,
                        k.scale = k.scale
                    )
                )
            }
        }
    }

    for (k in params$graph_k_grid) {
        iknn <- iknn.graph(case$X, k)
        candidates$add("iknn", "iKNN", sprintf("k=%s", k), iknn)
        candidates$add("mst_union_iknn", "MST union iKNN", sprintf("k=%s", k), union.graphs(case$mst, iknn))

        mknn <- mknn.graph(case$X, k)
        candidates$add("mknn", "mKNN", sprintf("k=%s", k), mknn)
        candidates$add("mst_union_mknn", "MST union mKNN", sprintf("k=%s", k), union.graphs(case$mst, mknn))
    }

    for (base.name in params$cycle_base) {
        base.graph <- make.cycle.base(base.name, case)
        for (n.add in params$cycle_n_add) {
            for (mult in params$cycle_local_multiplier) {
                for (ratio in params$cycle_min_ratio) {
                    candidates$add(
                        "cycle_repair",
                        "geodesic-discrepancy cycle repair",
                        sprintf("base=%s; n_add=%s; c=%s; min_ratio=%s", base.name, n.add, mult, ratio),
                        cycle.repair.graph(
                            base.graph,
                            case$X,
                            n.add = n.add,
                            local.multiplier = mult,
                            k.scale = 5L,
                            min.ratio = ratio
                        )
                    )
                }
            }
        }
    }

    candidates$get()
}

build.representative.candidates <- function(case) {
    candidates <- new.candidate.list()
    candidates$add("reference", "MST only", NA_character_, case$mst)
    candidates$add("reference", "oracle weighted circle", NA_character_, case$oracle)
    candidates$add("cmst", "CMST", "q=0.9", make.cmst.graph(case$X, q.thld = 0.9))
    candidates$add("cmst", "CMST", "q=0.99", make.cmst.graph(case$X, q.thld = 0.99))

    high.q <- make.cmst.graph(case$X, q.thld = 0.9999)
    candidates$add(
        "high_q_pruned",
        "high-q CMST + shortcut pruning",
        "q=0.9999; alt_ratio=1.05",
        prune.redundant.edges(high.q, alt.path.ratio = 1.05)
    )

    candidates$add(
        "local_adaptive",
        "MST + local adaptive threshold",
        "scale=knn; c=1.25; k_scale=5",
        adaptive.threshold.graph(
            case$X,
            mst.graph = case$mst,
            scale.source = "knn",
            multiplier = 1.25,
            k.scale = 5L
        )
    )

    iknn <- iknn.graph(case$X, 5L)
    candidates$add("iknn", "iKNN", "k=5", iknn)
    candidates$add("mst_union_iknn", "MST union iKNN", "k=5", union.graphs(case$mst, iknn))

    mknn <- mknn.graph(case$X, 9L)
    candidates$add("mknn", "mKNN", "k=9", mknn)
    candidates$add("mst_union_mknn", "MST union mKNN", "k=9", union.graphs(case$mst, mknn))

    candidates$add(
        "cycle_repair",
        "geodesic-discrepancy cycle repair",
        "base=cmst_q0.7; n_add=10; c=1.5; min_ratio=5",
        cycle.repair.graph(
            make.cmst.graph(case$X, q.thld = 0.7),
            case$X,
            n.add = 10L,
            local.multiplier = 1.5,
            k.scale = 5L,
            min.ratio = 5
        )
    )

    candidates$get()
}

score.candidates <- function(candidate.graphs, case, params) {
    metric.rows <- lapply(seq_along(candidate.graphs), function(i) {
        cand <- candidate.graphs[[i]]
        g <- cand$graph
        cbind(
            data.frame(
                id = i,
                family = cand$family,
                strategy = cand$strategy,
                parameter = cand$parameter
            ),
            graph.summary.metrics(
                g$adj_list,
                g$weight_list,
                true.dist = case$true.dist,
                shortcut.threshold = params$shortcut_arc_threshold
            ),
            geodesic.metrics(g$adj_list, g$weight_list, case$true.dist),
            oracle.edge.retention(g$adj_list, case$oracle$adj_list, case$theta, radius = params$radius)
        )
    })
    metrics <- do.call(rbind, metric.rows)
    metrics$score <- score.graph.metrics(metrics)
    metrics[order(metrics$score, metrics$median_relative_error), ]
}

top.by.family <- function(metrics, n = 5L) {
    rows <- lapply(split(metrics, metrics$family), function(x) {
        x <- x[order(x$score, x$median_relative_error), ]
        x[seq_len(min(n, nrow(x))), , drop = FALSE]
    })
    out <- do.call(rbind, rows)
    out[order(out$family, out$score), ]
}

summarize.replicates <- function(metrics) {
    group.keys <- unique(metrics[, c("n", "family", "strategy", "parameter")])
    param.key <- function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "<NA>"
        x
    }
    metric.param <- param.key(metrics$parameter)
    rows <- lapply(seq_len(nrow(group.keys)), function(i) {
        key <- group.keys[i, ]
        x <- metrics[
            metrics$n == key$n &
                metrics$family == key$family &
                metrics$strategy == key$strategy &
                metric.param == param.key(key$parameter),
        ]
        data.frame(
            n = key$n,
            family = key$family,
            strategy = key$strategy,
            parameter = key$parameter,
            runs = nrow(x),
            score_mean = mean(x$score, na.rm = TRUE),
            score_sd = stats::sd(x$score, na.rm = TRUE),
            median_relative_error_mean = mean(x$median_relative_error, na.rm = TRUE),
            median_relative_error_sd = stats::sd(x$median_relative_error, na.rm = TRUE),
            q90_relative_error_mean = mean(x$q90_relative_error, na.rm = TRUE),
            spearman_mean = mean(x$spearman, na.rm = TRUE),
            finite_pair_fraction_mean = mean(x$finite_pair_fraction, na.rm = TRUE),
            n_components_mean = mean(x$n_components, na.rm = TRUE),
            n_edges_mean = mean(x$n_edges, na.rm = TRUE),
            false_shortcut_rate_mean = mean(x$false_shortcut_rate, na.rm = TRUE),
            retained_oracle_edges_mean = mean(x$retained_oracle_edges, na.rm = TRUE)
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$n, out$score_mean, out$median_relative_error_mean), ]
}

strategy.docs <- data.frame(
    strategy = c(
        "CMST",
        "High-q CMST + shortcut pruning",
        "MST + local adaptive threshold",
        "iKNN and MST union iKNN",
        "mKNN and MST union mKNN",
        "Geodesic-discrepancy cycle repair"
    ),
    implementation = c(
        "Calls create.cmst.graph(). The graph starts with the MST and adds Euclidean edges whose lengths are no larger than the selected MST-edge quantile.",
        "Calls create.cmst.graph() with a permissive q.thld, then removes an edge when the graph remains connected and an alternate path between its endpoints is not much longer than the edge.",
        "Computes a local radius for each point and adds edges satisfying d(i,j) <= c * max(scale_i, scale_j), then unions those edges with the MST.",
        "Uses create.iknn.graphs() to create an intersection-kNN topology; the development wrapper reweights retained edges by Euclidean distance. The MST-union variant adds the MST as a connectedness scaffold.",
        "Uses create.mknn.graph() to create a mutual-kNN topology; the development wrapper reweights retained edges by Euclidean distance. The MST-union variant adds the MST as a connectedness scaffold.",
        "Starts from MST or conservative CMST. Candidate repair edges must be locally short in Euclidean distance but have a large current graph-distance / Euclidean-distance ratio."
    ),
    parameters = c(
        "q.thld controls the MST-edge quantile used as the completion distance threshold; pca.dim = NULL keeps the original two-dimensional noisy circle coordinates.",
        "q.thld controls how permissive the initial graph is; alt_ratio controls how short an alternate path must be before an edge is considered redundant.",
        "scale chooses the local scale source: knn uses the kth nearest-neighbor distance, mst uses local MST incident lengths. c is the multiplier on that scale. k_scale is the neighbor rank used when scale = knn.",
        "k is the neighbor rank. In this report, k = 5 is the representative recipe and k = 2:10 is swept in the diagnostic run.",
        "k is the mutual-neighbor rank. In this report, k = 9 is the representative recipe and k = 2:10 is swept in the diagnostic run.",
        "base chooses the starting graph; n_add is the number of repair edges; c is the local scale multiplier; min_ratio is the minimum graph-distance / Euclidean-distance discrepancy."
    ),
    command = c(
        "cmst <- gflow::create.cmst.graph(X, q.thld = 0.9, pca.dim = NULL, verbose = FALSE)\ng <- list(adj_list = cmst$cmst_adj_list, weight_list = cmst$cmst_weight_list)",
        "base <- gflow::create.cmst.graph(X, q.thld = 0.9999, pca.dim = NULL, verbose = FALSE)\nbase.graph <- list(adj_list = base$cmst_adj_list, weight_list = base$cmst_weight_list)\ng <- prune.redundant.edges(base.graph, alt.path.ratio = 1.05)",
        "mst <- gflow::create.cmst.graph(X, q.thld = 0.9, pca.dim = NULL, verbose = FALSE)\nmst.graph <- list(adj_list = mst$mst_adj_list, weight_list = mst$mst_weight_list)\ng <- adaptive.threshold.graph(X, mst.graph, scale.source = \"knn\", multiplier = 1.25, k.scale = 5L)",
        "iknn <- iknn.graph(X, k = 5L)\ng.local <- iknn\ng.union <- union.graphs(mst.graph, iknn)",
        "mknn <- mknn.graph(X, k = 9L)\ng.local <- mknn\ng.union <- union.graphs(mst.graph, mknn)",
        "base <- gflow::create.cmst.graph(X, q.thld = 0.7, pca.dim = NULL, verbose = FALSE)\nbase.graph <- list(adj_list = base$cmst_adj_list, weight_list = base$cmst_weight_list)\ng <- cycle.repair.graph(base.graph, X, n.add = 10L, local.multiplier = 1.5, k.scale = 5L, min.ratio = 5)"
    )
)

html.strategy.docs <- function(docs) {
    paste(vapply(seq_len(nrow(docs)), function(i) {
        paste0(
            "<section class=\"strategy-doc\"><h3>", html.escape(docs$strategy[[i]]), "</h3>",
            "<p>", html.escape(docs$implementation[[i]]), "</p>",
            "<p><strong>Parameters:</strong> ", html.escape(docs$parameters[[i]]), "</p>",
            "<pre><code>", html.escape(docs$command[[i]]), "</code></pre></section>"
        )
    }, character(1)), collapse = "\n")
}

single.case <- make.noisy.circle.case(params$n, params$seed, params)
exhaustive.graphs <- build.exhaustive.candidates(single.case, params)
metrics <- score.candidates(exhaustive.graphs, single.case, params)

best.by.family <- do.call(rbind, lapply(split(metrics, metrics$family), function(x) {
    x[which.min(x$score), , drop = FALSE]
}))
best.by.family <- best.by.family[order(best.by.family$score), ]
diagnostic.top <- top.by.family(metrics, n = 5L)

union.comparison <- metrics[metrics$family %in% c("iknn", "mst_union_iknn", "mknn", "mst_union_mknn"), ]
union.comparison$k <- as.integer(sub("^k=", "", union.comparison$parameter))
union.comparison <- union.comparison[order(union.comparison$strategy, union.comparison$k), ]

union.delta <- do.call(rbind, lapply(list(
    c(base = "iKNN", union = "MST union iKNN"),
    c(base = "mKNN", union = "MST union mKNN")
), function(pair) {
    base <- union.comparison[union.comparison$strategy == pair[["base"]], ]
    union <- union.comparison[union.comparison$strategy == pair[["union"]], ]
    base <- base[order(base$k), ]
    union <- union[match(base$k, union$k), ]
    data.frame(
        graph = pair[["base"]],
        k = base$k,
        base_components = base$n_components,
        union_components = union$n_components,
        delta_finite_pair_fraction = union$finite_pair_fraction - base$finite_pair_fraction,
        delta_median_relative_error = union$median_relative_error - base$median_relative_error,
        delta_q90_relative_error = union$q90_relative_error - base$q90_relative_error,
        extra_edges = union$n_edges - base$n_edges
    )
}))

selected.ids <- unique(c(
    best.by.family$id[best.by.family$family %in% c("high_q_pruned", "local_adaptive",
                                                   "iknn", "mst_union_iknn",
                                                   "mknn", "mst_union_mknn",
                                                   "cycle_repair")],
    metrics$id[metrics$family == "reference" & metrics$strategy == "oracle weighted circle"]
))
selected <- metrics[match(selected.ids, metrics$id), ]
selected <- selected[!is.na(selected$id), ]

replicate.grid <- expand.grid(
    n = params$replicate_n,
    seed = params$replicate_seeds,
    KEEP.OUT.ATTRS = FALSE
)
replicate.rows <- vector("list", nrow(replicate.grid))
for (r in seq_len(nrow(replicate.grid))) {
    n.run <- replicate.grid$n[[r]]
    seed.run <- replicate.grid$seed[[r]]
    message(sprintf("Replicate benchmark n=%s seed=%s", n.run, seed.run))
    case.run <- make.noisy.circle.case(n.run, seed.run, params)
    graphs.run <- build.representative.candidates(case.run)
    m.run <- score.candidates(graphs.run, case.run, params)
    m.run$n <- n.run
    m.run$seed <- seed.run
    replicate.rows[[r]] <- m.run
}
replicate.metrics <- do.call(rbind, replicate.rows)
replicate.summary <- summarize.replicates(replicate.metrics)
replicate.best <- do.call(rbind, lapply(split(replicate.summary, replicate.summary$n), function(x) {
    x <- x[x$strategy != "oracle weighted circle", , drop = FALSE]
    x[seq_len(min(8L, nrow(x))), , drop = FALSE]
}))

saveRDS(
    list(
        params = params,
        X = single.case$X,
        theta = single.case$theta,
        true.dist = single.case$true.dist,
        metrics = metrics,
        best.by.family = best.by.family,
        diagnostic.top = diagnostic.top,
        union.comparison = union.comparison,
        union.delta = union.delta,
        selected = selected,
        candidate.graphs = exhaustive.graphs,
        replicate.metrics = replicate.metrics,
        replicate.summary = replicate.summary,
        replicate.best = replicate.best,
        strategy.docs = strategy.docs
    ),
    file.path(cache.dir, "noisy_circle_strategy_comparison_results.rds")
)

graph.by.id <- function(id) exhaustive.graphs[[id]]$graph
label.by.row <- function(row) sprintf("%s\n%s", row$strategy, row$parameter)
shortcut.fun <- function(edges) single.case$true.dist[cbind(edges$from, edges$to)] > params$shortcut_arc_threshold

fig.files <- list()

fig.files$best.overlays <- file.path(fig.dir, "01_best_strategy_overlays.png")
save.png(fig.files$best.overlays, {
    oldpar <- par(mfrow = c(2, 4), mar = c(3.6, 3.6, 3.2, 0.8))
    on.exit(par(oldpar), add = TRUE)
    for (r in seq_len(min(8L, nrow(selected)))) {
        g <- graph.by.id(selected$id[[r]])
        draw.graph.overlay(
            single.case$X, single.case$theta, g$adj_list, g$weight_list,
            main = label.by.row(selected[r, ]),
            shortcut.fun = shortcut.fun
        )
    }
}, width = 2100, height = 1100)

fig.files$best.distances <- file.path(fig.dir, "02_best_strategy_distances.png")
save.png(fig.files$best.distances, {
    oldpar <- par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    for (r in seq_len(min(8L, nrow(selected)))) {
        g <- graph.by.id(selected$id[[r]])
        draw.distance.scatter(g$adj_list, g$weight_list, single.case$true.dist, main = selected$strategy[[r]])
    }
}, width = 2100, height = 1200)

fig.files$pareto <- file.path(fig.dir, "03_strategy_pareto.png")
save.png(fig.files$pareto, {
    fam <- factor(metrics$family)
    cols <- grDevices::hcl.colors(nlevels(fam), "Dark 3")
    oldpar <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    plot(metrics$n_edges, metrics$median_relative_error, pch = 19, col = cols[fam],
         xlab = "edge count", ylab = "median relative geodesic error",
         main = "Error vs graph size")
    plot(metrics$false_shortcut_rate, metrics$median_relative_error, pch = 19, col = cols[fam],
         xlab = "false shortcut rate", ylab = "median relative geodesic error",
         main = "Error vs shortcuts")
    plot(metrics$cycle_rank, metrics$spearman, pch = 19, col = cols[fam],
         xlab = "cycle rank", ylab = "Spearman(graph, truth)",
         main = "Topology vs recovery")
    legend("bottomright", legend = levels(fam), col = cols, pch = 19, cex = 0.72, bty = "n")
}, width = 1900, height = 700)

fig.files$union <- file.path(fig.dir, "04_mst_union_effects.png")
save.png(fig.files$union, {
    oldpar <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    strategy.cols <- c(
        "iKNN" = "#1f77b4",
        "MST union iKNN" = "#1f77b4",
        "mKNN" = "#d95f02",
        "MST union mKNN" = "#d95f02"
    )
    strategy.lty <- c(
        "iKNN" = 2,
        "MST union iKNN" = 1,
        "mKNN" = 2,
        "MST union mKNN" = 1
    )
    plot(NA, xlim = range(params$graph_k_grid), ylim = range(union.comparison$median_relative_error, na.rm = TRUE),
         xlab = "k", ylab = "median relative geodesic error",
         main = "Does MST union improve error?")
    for (strategy in unique(union.comparison$strategy)) {
        x <- union.comparison[union.comparison$strategy == strategy, ]
        lines(x$k, x$median_relative_error, type = "b", pch = 19,
              col = strategy.cols[[strategy]], lty = strategy.lty[[strategy]])
    }
    legend("topright", legend = unique(union.comparison$strategy),
           col = strategy.cols[unique(union.comparison$strategy)],
           lty = strategy.lty[unique(union.comparison$strategy)],
           pch = 19, cex = 0.75, bty = "n")
    plot(NA, xlim = range(params$graph_k_grid), ylim = range(union.comparison$finite_pair_fraction, na.rm = TRUE),
         xlab = "k", ylab = "finite pair fraction",
         main = "Does MST union fix disconnectedness?")
    for (strategy in unique(union.comparison$strategy)) {
        x <- union.comparison[union.comparison$strategy == strategy, ]
        lines(x$k, x$finite_pair_fraction, type = "b", pch = 19,
              col = strategy.cols[[strategy]], lty = strategy.lty[[strategy]])
    }
    plot(NA, xlim = range(params$graph_k_grid), ylim = range(union.comparison$n_edges, na.rm = TRUE),
         xlab = "k", ylab = "edge count",
         main = "MST union cost")
    for (strategy in unique(union.comparison$strategy)) {
        x <- union.comparison[union.comparison$strategy == strategy, ]
        lines(x$k, x$n_edges, type = "b", pch = 19,
              col = strategy.cols[[strategy]], lty = strategy.lty[[strategy]])
    }
}, width = 1900, height = 700)

strategy.levels <- unique(replicate.summary$strategy)
strategy.cols <- setNames(grDevices::hcl.colors(length(strategy.levels), "Dark 3"), strategy.levels)
strategy.labels <- c(
    "MST only" = "MST",
    "CMST" = "CMST",
    "high-q CMST + shortcut pruning" = "CMST+prune",
    "MST + local adaptive threshold" = "local adaptive",
    "iKNN" = "iKNN",
    "MST union iKNN" = "MST+iKNN",
    "mKNN" = "mKNN",
    "MST union mKNN" = "MST+mKNN",
    "geodesic-discrepancy cycle repair" = "cycle repair",
    "oracle weighted circle" = "oracle"
)
label.strategy <- function(x) {
    out <- unname(strategy.labels[as.character(x)])
    out[is.na(out)] <- as.character(x)[is.na(out)]
    out
}

fig.files$rep.error <- file.path(fig.dir, "05_replicated_error_by_size.png")
save.png(fig.files$rep.error, {
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1.15))
    oldpar <- par(mar = c(4.8, 4.8, 3.2, 1))
    on.exit(par(oldpar), add = TRUE)
    human <- replicate.summary[replicate.summary$strategy != "oracle weighted circle", ]
    plot(NA, xlim = range(params$replicate_n), ylim = range(human$median_relative_error_mean, na.rm = TRUE),
         xlab = "sample size n", ylab = "mean median relative geodesic error",
         main = "Replicated geodesic error")
    for (strategy in unique(human$strategy)) {
        x <- human[human$strategy == strategy, ]
        lines(x$n, x$median_relative_error_mean, type = "b", pch = 19, col = strategy.cols[[strategy]])
    }
    plot(NA, xlim = range(params$replicate_n), ylim = range(human$q90_relative_error_mean, na.rm = TRUE),
         xlab = "sample size n", ylab = "mean q90 relative geodesic error",
         main = "Upper-tail distortion")
    for (strategy in unique(human$strategy)) {
        x <- human[human$strategy == strategy, ]
        lines(x$n, x$q90_relative_error_mean, type = "b", pch = 19, col = strategy.cols[[strategy]])
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = label.strategy(unique(human$strategy)),
           col = strategy.cols[unique(human$strategy)], lty = 1, pch = 19,
           cex = 0.9, bty = "n", ncol = 5)
}, width = 1800, height = 900)

fig.files$rep.connectivity <- file.path(fig.dir, "06_replicated_connectivity_by_size.png")
save.png(fig.files$rep.connectivity, {
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1.15))
    oldpar <- par(mar = c(4.8, 4.8, 3.2, 1))
    on.exit(par(oldpar), add = TRUE)
    human <- replicate.summary[replicate.summary$strategy != "oracle weighted circle", ]
    plot(NA, xlim = range(params$replicate_n), ylim = c(0, 1),
         xlab = "sample size n", ylab = "mean finite pair fraction",
         main = "Graph-distance coverage")
    for (strategy in unique(human$strategy)) {
        x <- human[human$strategy == strategy, ]
        lines(x$n, x$finite_pair_fraction_mean, type = "b", pch = 19, col = strategy.cols[[strategy]])
    }
    plot(NA, xlim = range(params$replicate_n), ylim = range(human$n_components_mean, na.rm = TRUE),
         xlab = "sample size n", ylab = "mean component count",
         main = "Disconnectedness")
    for (strategy in unique(human$strategy)) {
        x <- human[human$strategy == strategy, ]
        lines(x$n, x$n_components_mean, type = "b", pch = 19, col = strategy.cols[[strategy]])
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = label.strategy(unique(human$strategy)),
           col = strategy.cols[unique(human$strategy)], lty = 1, pch = 19,
           cex = 0.9, bty = "n", ncol = 5)
}, width = 1800, height = 900)

fig.files$rep.box <- file.path(fig.dir, "07_replicated_error_distributions.png")
save.png(fig.files$rep.box, {
    oldpar <- par(mfrow = c(1, length(params$replicate_n)), mar = c(6.5, 4, 3, 0.8))
    on.exit(par(oldpar), add = TRUE)
    human <- replicate.metrics[replicate.metrics$strategy != "oracle weighted circle", ]
    human$strategy_label <- factor(label.strategy(human$strategy), levels = label.strategy(unique(human$strategy)))
    for (n.value in params$replicate_n) {
        x <- human[human$n == n.value, ]
        boxplot(
            median_relative_error ~ strategy_label,
            data = x,
            las = 2,
            cex.axis = 0.76,
            ylab = if (n.value == params$replicate_n[[1]]) "median relative error" else "",
            xlab = "strategy",
            main = sprintf("n = %s", n.value),
            col = "gray90"
        )
    }
}, width = 2100, height = 900)

rel.fig <- function(path) {
    file.path("../../figures/noisy-circle-strategy-comparison", basename(path))
}

commands <- paste(
    "set.seed(1002L)",
    "X.df <- gflow::generate.circle.data(",
    "  n = 160L, radius = 1, noise = 0.08,",
    "  type = \"random\", noise.type = \"normal\", seed = 1002L",
    ")",
    "X <- as.matrix(X.df[, c(\"x\", \"y\")])",
    "theta <- as.double(X.df$angles)",
    "true.dist <- pmin(abs(outer(theta, theta, \"-\")),",
    "                   2 * pi - abs(outer(theta, theta, \"-\")))",
    "mst <- gflow::create.cmst.graph(X, q.thld = 0.9, pca.dim = NULL, verbose = FALSE)",
    "## The replicated benchmark uses n = 80, 160, 320 and seeds = 3001:3005.",
    sep = "\n"
)

report.path <- file.path(report.dir, "noisy_circle_strategy_comparison.html")
html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Noisy Circle Strategy Comparison</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;max-width:1220px;}h1,h2,h3{line-height:1.15;}code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f6f8fb;}figure{margin:24px 0;}img{max-width:100%;border:1px solid #d9dee7;}figcaption{font-size:13px;color:#52616f;margin-top:6px;}.strategy-doc{border-top:1px solid #d9dee7;padding-top:10px;margin-top:14px;}details{border:1px solid #d9dee7;padding:10px 14px;margin:18px 0;background:#fbfcfe;}summary{font-weight:600;cursor:pointer;}</style>",
    "</head><body>",
    "<h1>Noisy Circle Strategy Comparison</h1>",
    sprintf("<p>Generated at %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "<h2>Objective</h2>",
    "<p>This report tests MST-completion strategy families on noisy circle examples. The response variable is ignored. The only target is recovery of the true circular geodesic distance.</p>",
    "<pre><code>d_true(i, j) = min(|theta_i - theta_j|, 2*pi - |theta_i - theta_j|)</code></pre>",
    "<h2>Benchmark Design</h2>",
    sprintf("<p>The diagnostic sweep uses one medium noisy circle with <code>n = %s</code> and <code>seed = %s</code>. The human-facing replicated benchmark uses sample sizes <code>%s</code> and seeds <code>%s</code>. This separates parameter exploration from stability assessment across random samplings.</p>",
            params$n, params$seed, paste(params$replicate_n, collapse = ", "), paste(params$replicate_seeds, collapse = ", ")),
    "<h2>Reproducible Data Commands</h2><pre><code>", html.escape(commands), "</code></pre>",
    "<h2>Strategy Recipes</h2>",
    "<p>Each recipe below names the construction, the implementation used in this development report, the parameters that matter, and the R command skeleton. Helper functions such as <code>adaptive.threshold.graph()</code>, <code>prune.redundant.edges()</code>, and <code>cycle.repair.graph()</code> live in <code>dev/geodesic-distance-estimation/R/graph_construction.R</code>.</p>",
    html.strategy.docs(strategy.docs),
    "<h2>Replicated Size Benchmark</h2>",
    "<p>The figures below are the main human-facing comparison. They average representative strategy recipes across several random noisy-circle samplings at each size. The goal is to see whether a strategy is stable, not merely best in one draw.</p>",
    sprintf("<figure><img src=\"%s\"><figcaption>Mean geodesic distance error and upper-tail distortion across sample sizes and random seeds.</figcaption></figure>", rel.fig(fig.files$rep.error)),
    sprintf("<figure><img src=\"%s\"><figcaption>Connectivity behavior across sample sizes. MST-union variants should reach finite-pair coverage by construction; local-only graphs may not.</figcaption></figure>", rel.fig(fig.files$rep.connectivity)),
    sprintf("<figure><img src=\"%s\"><figcaption>Run-to-run distribution of median relative geodesic error. This is more useful than a single huge candidate table for judging robustness.</figcaption></figure>", rel.fig(fig.files$rep.box)),
    "<h3>Best Replicated Recipes By Size</h3>",
    html.table(replicate.best[, c("n", "family", "strategy", "parameter", "runs", "score_mean",
                                  "median_relative_error_mean", "q90_relative_error_mean",
                                  "spearman_mean", "finite_pair_fraction_mean", "n_edges_mean")]),
    "<details><summary>Replicated summary table for completeness</summary>",
    html.table(replicate.summary[, c("n", "family", "strategy", "parameter", "runs", "score_mean", "score_sd",
                                     "median_relative_error_mean", "median_relative_error_sd",
                                     "q90_relative_error_mean", "spearman_mean",
                                     "finite_pair_fraction_mean", "n_components_mean", "n_edges_mean")]),
    "</details>",
    "<h2>Single-Run Diagnostic Sweep</h2>",
    "<p>The single-run sweep is still useful for understanding parameter sensitivity and for selecting recipes to promote into replicated comparisons. Its full table is long, so the report first shows best candidates by family and top candidates within each family; the complete table is folded below for completeness.</p>",
    html.table(best.by.family[, c("family", "strategy", "parameter", "score", "n_edges", "n_components",
                                  "cycle_rank", "spearman", "median_relative_error", "q90_relative_error",
                                  "mean_neighbor_overlap", "false_shortcut_rate", "retained_oracle_edges",
                                  "missing_oracle_edges")]),
    sprintf("<figure><img src=\"%s\"><figcaption>Best-scoring representative graph from each strategy family for the diagnostic run. Red edges exceed the report's pi/4 true-arc shortcut diagnostic.</figcaption></figure>", rel.fig(fig.files$best.overlays)),
    sprintf("<figure><img src=\"%s\"><figcaption>Graph shortest-path distances versus true circular geodesic distances for best-scoring diagnostic representatives.</figcaption></figure>", rel.fig(fig.files$best.distances)),
    "<h3>MST Union Diagnostic</h3>",
    "<p>The pairwise delta table uses <code>union - local</code>, so negative error deltas are improvements and positive finite-pair deltas mean the MST scaffold fixed disconnected graph distances.</p>",
    html.table(union.delta),
    sprintf("<figure><img src=\"%s\"><figcaption>Effect of MST union on iKNN and mKNN across k in the diagnostic run.</figcaption></figure>", rel.fig(fig.files$union)),
    "<h3>Top Diagnostic Candidates Within Each Family</h3>",
    html.table(diagnostic.top[, c("family", "strategy", "parameter", "score", "n_edges", "n_components",
                                  "cycle_rank", "finite_pair_fraction", "spearman",
                                  "median_relative_error", "q90_relative_error",
                                  "mean_neighbor_overlap", "retained_oracle_edges")]),
    sprintf("<figure><img src=\"%s\"><figcaption>Pareto-style diagnostics across all single-run candidates.</figcaption></figure>", rel.fig(fig.files$pareto)),
    "<details><summary>Complete single-run candidate table for completeness</summary>",
    html.table(metrics[, c("family", "strategy", "parameter", "score", "n_edges", "n_components",
                           "cycle_rank", "finite_pair_fraction", "spearman", "median_relative_error",
                           "q90_relative_error", "mean_neighbor_overlap", "false_shortcut_rate",
                           "retained_oracle_edges", "missing_oracle_edges")]),
    "</details>",
    "<h2>Initial Interpretation</h2>",
    "<p>For this noisy circle, a successful graph needs at least one cycle and should recover circular shortest-path distances without adding chord shortcuts. The replicated benchmark is the more reliable comparison; the single-run sweep is a parameter-discovery tool.</p>",
    "<p><strong>MST union iKNN/mKNN:</strong> MST union is most useful when the local graph is disconnected. It improves finite-pair coverage by construction. When iKNN or mKNN is already connected, the union should be judged empirically because extra MST edges can shorten paths in ways that help sparse regions or inject MST artifacts.</p>",
    "<h2>Next Checks</h2>",
    "<ul><li>Repeat the replicated strategy comparison on a nonuniform noisy circle.</li><li>Run the same representative candidates on a crescent or spiral to see which circle-winning repairs create false shortcuts.</li><li>Promote only stable strategy families into package-facing prototypes.</li></ul>",
    "</body></html>"
)
writeLines(html, report.path)

cat("Wrote report:\n", report.path, "\n", sep = "")
cat("Wrote figures under:\n", fig.dir, "\n", sep = "")
cat("Wrote cache under:\n", cache.dir, "\n", sep = "")

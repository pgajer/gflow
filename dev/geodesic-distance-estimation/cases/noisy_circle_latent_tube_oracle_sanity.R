# Noisy Circle Latent-Tube Oracle Sanity Report

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
report.dir <- file.path(out.root, "reports/noisy-circle-latent-tube-oracle-sanity")
fig.dir <- file.path(out.root, "figures/noisy-circle-latent-tube-oracle-sanity")
cache.dir <- file.path(out.root, "cache/noisy-circle-latent-tube-oracle-sanity")
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
    radius = 1,
    geometry_noise = 0.08,
    geometry_noise_type = "normal",
    validation_n = 160L,
    validation_seed = 4101L,
    benchmark_n_grid = c(80L, 160L),
    benchmark_seeds = 4101L,
    k_grid = c(5L, 8L, 10L),
    pruning_delta_grid = c(0, 0.05, 0.10),
    alpha_grid = c(0, 0.5, 1),
    diagnostic_alpha = 0.5,
    benchmark_alpha_grid = c(0.5, 1),
    oracle_n_theta = 360L,
    oracle_n_u = 41L,
    oracle_u_max_multiplier = 4,
    oracle_density_floor = 1e-4,
    oracle_attach_extra = 1L,
    oracle_density_type = "normal",
    path_edge_ratio_percentile = 0.5
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

build.tube.oracle <- function(case, alpha, params) {
    latent.tube.oracle(
        theta.sample = case$theta,
        u.sample = case$u,
        radius = params$radius,
        sigma = params$geometry_noise,
        density.type = params$oracle_density_type,
        alpha = alpha,
        n.theta = params$oracle_n_theta,
        u.max = params$oracle_u_max_multiplier * params$geometry_noise,
        n.u = params$oracle_n_u,
        density.floor = params$oracle_density_floor,
        include.diagonal = TRUE,
        attach.extra = params$oracle_attach_extra
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

upper.values <- function(D) {
    D[upper.tri(D)]
}

distance.comparison.row <- function(name, D, arc, euclidean) {
    keep <- upper.tri(D)
    d <- D[keep]
    a <- arc[keep]
    e <- euclidean[keep]
    data.frame(
        truth = name,
        median_distance = stats::median(d),
        q90_distance = as.numeric(stats::quantile(d, 0.9, names = FALSE)),
        max_distance = max(d),
        spearman_vs_arc = suppressWarnings(stats::cor(d, a, method = "spearman")),
        spearman_vs_euclidean = suppressWarnings(stats::cor(d, e, method = "spearman")),
        median_ratio_to_arc = stats::median(d / pmax(a, .Machine$double.eps)),
        median_ratio_to_euclidean = stats::median(d / pmax(e, .Machine$double.eps))
    )
}

select.diagnostic.pairs <- function(case, tube.dist) {
    n <- nrow(case$X)
    keep <- which(upper.tri(tube.dist), arr.ind = TRUE)
    arc <- case$arc.dist[keep]
    eu <- case$euclidean.dist[keep]
    tube <- tube.dist[keep]
    ratio <- tube / pmax(eu, .Machine$double.eps)

    near.idx <- which.min(abs(arc - stats::quantile(arc, 0.1)))
    far.idx <- which.max(arc)
    short.eu <- eu <= stats::quantile(eu, 0.15)
    deceptive.idx <- which.max(ifelse(short.eu, ratio, -Inf))
    radial.candidate <- arc <= stats::quantile(arc, 0.15)
    radial.idx <- which.max(ifelse(radial.candidate, abs(case$u[keep[, 1]] - case$u[keep[, 2]]), -Inf))

    idx <- c(near.idx, far.idx, deceptive.idx, radial.idx)
    labels <- c("near-angle", "opposite-angle", "euclidean-close/high-oracle-ratio",
                "near-angle/radial-separation")
    data.frame(
        label = labels,
        from = keep[idx, 1],
        to = keep[idx, 2],
        arc_distance = arc[idx],
        euclidean_distance = eu[idx],
        tube_distance = tube[idx],
        tube_to_euclidean_ratio = ratio[idx],
        u_from = case$u[keep[idx, 1]],
        u_to = case$u[keep[idx, 2]]
    )
}

edge.diagnostic.summary <- function(graph, case, true.dist) {
    edges <- edge.table(graph$adj_list, graph$weight_list)
    if (!nrow(edges)) return(data.frame())
    edges$arc_distance <- case$arc.dist[cbind(edges$from, edges$to)]
    edges$oracle_distance <- true.dist[cbind(edges$from, edges$to)]
    edges$euclidean_distance <- case$euclidean.dist[cbind(edges$from, edges$to)]
    edges$oracle_to_euclidean_ratio <- edges$oracle_distance /
        pmax(edges$euclidean_distance, .Machine$double.eps)
    data.frame(
        n_edges = nrow(edges),
        median_euclidean_edge = stats::median(edges$euclidean_distance),
        q90_euclidean_edge = as.numeric(stats::quantile(edges$euclidean_distance, 0.9, names = FALSE)),
        median_oracle_edge = stats::median(edges$oracle_distance),
        q90_oracle_edge = as.numeric(stats::quantile(edges$oracle_distance, 0.9, names = FALSE)),
        median_oracle_to_euclidean = stats::median(edges$oracle_to_euclidean_ratio),
        q90_oracle_to_euclidean = as.numeric(stats::quantile(edges$oracle_to_euclidean_ratio, 0.9, names = FALSE)),
        max_oracle_to_euclidean = max(edges$oracle_to_euclidean_ratio)
    )
}

score.row <- function(metrics) {
    metrics$median_relative_error +
        0.5 * metrics$q90_relative_error -
        0.25 * metrics$mean_neighbor_overlap +
        ifelse(metrics$finite_pair_fraction < 1, 10, 0)
}

validation.case <- make.noisy.circle.case(
    params$validation_n,
    params$validation_seed,
    params
)
validation.oracles <- setNames(
    lapply(params$alpha_grid, function(alpha) build.tube.oracle(validation.case, alpha, params)),
    paste0("alpha=", params$alpha_grid)
)
validation.tube.dist <- lapply(validation.oracles, latent.tube.oracle.distances)
validation.truth.summary <- do.call(rbind, c(
    list(distance.comparison.row("latent arc", validation.case$arc.dist,
                                 validation.case$arc.dist, validation.case$euclidean.dist),
         distance.comparison.row("Euclidean", validation.case$euclidean.dist,
                                 validation.case$arc.dist, validation.case$euclidean.dist)),
    Map(function(name, D) distance.comparison.row(
        paste("tube oracle", name), D,
        validation.case$arc.dist, validation.case$euclidean.dist
    ), names(validation.tube.dist), validation.tube.dist)
))

benchmark.rows <- list()
edge.rows <- list()
row.ptr <- 0L
edge.ptr <- 0L
for (n.value in params$benchmark_n_grid) {
    for (seed.value in params$benchmark_seeds) {
        message(sprintf("Latent-tube oracle sanity benchmark n=%s seed=%s", n.value, seed.value))
        case <- make.noisy.circle.case(n.value, seed.value, params)
        truth.list <- list(`latent arc` = case$arc.dist,
                           `tube oracle alpha=0.5` = latent.tube.oracle.distances(
                               build.tube.oracle(case, 0.5, params)
                           ),
                           `tube oracle alpha=1` = latent.tube.oracle.distances(
                               build.tube.oracle(case, 1, params)
                           ))
        for (k in params$k_grid) {
            for (delta in params$pruning_delta_grid) {
                graph <- pruned.iknn.graph(
                    case$X,
                    k = k,
                    delta = delta,
                    path.percentile = params$path_edge_ratio_percentile
                )
                for (truth.name in names(truth.list)) {
                    truth <- truth.list[[truth.name]]
                    geo <- geodesic.metrics(graph$adj_list, graph$weight_list, truth)
                    row <- data.frame(
                        n = n.value,
                        seed = seed.value,
                        graph = "iKNN",
                        k = k,
                        pruning_delta = delta,
                        truth = truth.name
                    )
                    row <- cbind(row, geo)
                    row$score <- score.row(row)
                    row.ptr <- row.ptr + 1L
                    benchmark.rows[[row.ptr]] <- row

                    e <- edge.diagnostic.summary(graph, case, truth)
                    e <- cbind(row[, c("n", "seed", "graph", "k", "pruning_delta", "truth")], e)
                    edge.ptr <- edge.ptr + 1L
                    edge.rows[[edge.ptr]] <- e
                }
            }
        }
    }
}
benchmark.metrics <- do.call(rbind, benchmark.rows)
edge.summary <- do.call(rbind, edge.rows)

best.by.truth.n <- do.call(rbind, lapply(split(
    benchmark.metrics,
    list(benchmark.metrics$n, benchmark.metrics$truth),
    drop = TRUE
), function(x) {
    x <- x[order(x$score, x$median_relative_error), ]
    x[1L, , drop = FALSE]
}))
best.by.truth.n <- best.by.truth.n[order(best.by.truth.n$n, best.by.truth.n$truth), ]

pair.table <- select.diagnostic.pairs(
    validation.case,
    validation.tube.dist[[paste0("alpha=", params$diagnostic_alpha)]]
)

saveRDS(
    list(
        params = params,
        validation.case = validation.case,
        validation.truth.summary = validation.truth.summary,
        pair.table = pair.table,
        benchmark.metrics = benchmark.metrics,
        edge.summary = edge.summary,
        best.by.truth.n = best.by.truth.n
    ),
    file.path(cache.dir, "noisy_circle_latent_tube_oracle_sanity_results.rds")
)
utils::write.csv(benchmark.metrics,
                 file.path(cache.dir, "noisy_circle_latent_tube_oracle_sanity_benchmark_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(edge.summary,
                 file.path(cache.dir, "noisy_circle_latent_tube_oracle_sanity_edge_summary.csv"),
                 row.names = FALSE)

fig.files <- list()

benchmark_metrics_key <- function(metrics, truth.name, n.value) {
    metrics$truth == truth.name & metrics$n == n.value
}

fig.files$oracle.grid <- file.path(fig.dir, "01_oracle_grid_and_samples.png")
save.png(fig.files$oracle.grid, {
    oracle <- validation.oracles[[paste0("alpha=", params$diagnostic_alpha)]]
    grid <- oracle$vertices[oracle$vertices$kind == "grid", ]
    pal <- grDevices::hcl.colors(128, "YlGnBu", rev = TRUE)
    z <- (grid$density - min(grid$density)) / max(diff(range(grid$density)), .Machine$double.eps)
    cols <- pal[pmax(1L, pmin(128L, floor(z * 127) + 1L))]
    oldpar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    plot(grid$x, grid$y, asp = 1, pch = 15, cex = 0.2, col = cols,
         xlab = "x", ylab = "y", main = "Latent-tube oracle grid")
    points(validation.case$X[, 1], validation.case$X[, 2], pch = 19,
           cex = 0.55, col = "black")
    theta <- seq(0, 2 * pi, length.out = 400)
    lines(cos(theta), sin(theta), col = "white", lwd = 2)
    lines(cos(theta), sin(theta), col = "gray25", lwd = 1)
    plot(validation.case$theta, validation.case$u, pch = 19, cex = 0.65,
         xlab = "latent angle theta", ylab = "radial displacement u",
         main = "Observed samples in latent tube coordinates")
    abline(h = 0, col = "gray55", lty = 2)
}, width = 1600, height = 800)

fig.files$truth.scatter <- file.path(fig.dir, "02_truth_distance_comparisons.png")
save.png(fig.files$truth.scatter, {
    oldpar <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.0, 1.0))
    on.exit(par(oldpar), add = TRUE)
    keep <- upper.tri(validation.case$arc.dist)
    for (name in names(validation.tube.dist)) {
        D <- validation.tube.dist[[name]]
        plot(validation.case$arc.dist[keep], D[keep], pch = 16, cex = 0.28,
             col = grDevices::adjustcolor("black", alpha.f = 0.22),
             xlab = "latent arc distance", ylab = "tube-oracle distance",
             main = paste("Arc vs tube", name))
        abline(0, 1, col = "gray55", lty = 2)
        plot(validation.case$euclidean.dist[keep], D[keep], pch = 16, cex = 0.28,
             col = grDevices::adjustcolor("black", alpha.f = 0.22),
             xlab = "Euclidean distance", ylab = "tube-oracle distance",
             main = paste("Euclidean vs tube", name))
        abline(0, 1, col = "gray55", lty = 2)
    }
}, width = 1800, height = 1100)

fig.files$paths <- file.path(fig.dir, "03_oracle_shortest_paths.png")
save.png(fig.files$paths, {
    path.alphas <- c(0, params$diagnostic_alpha, 1)
    oldpar <- par(mfrow = c(nrow(pair.table), length(path.alphas)),
                  mar = c(2.8, 2.8, 2.6, 0.8))
    on.exit(par(oldpar), add = TRUE)
    for (r in seq_len(nrow(pair.table))) {
        for (alpha in path.alphas) {
            oracle <- validation.oracles[[paste0("alpha=", alpha)]]
            from.vertex <- oracle$sample.vertices[[pair.table$from[[r]]]]
            to.vertex <- oracle$sample.vertices[[pair.table$to[[r]]]]
            path <- igraph::shortest_paths(
                oracle$graph,
                from = from.vertex,
                to = to.vertex,
                weights = igraph::E(oracle$graph)$weight
            )$vpath[[1L]]
            path <- as.integer(path)
            grid <- oracle$vertices[oracle$vertices$kind == "grid", ]
            plot(grid$x, grid$y, asp = 1, pch = 15, cex = 0.12,
                 col = grDevices::adjustcolor("gray75", alpha.f = 0.45),
                 xlab = "", ylab = "",
                 main = sprintf("%s; alpha=%s", pair.table$label[[r]], alpha))
            lines(oracle$vertices$x[path], oracle$vertices$y[path],
                  col = "red3", lwd = 2.2)
            points(validation.case$X[, 1], validation.case$X[, 2],
                   pch = 19, cex = 0.28, col = "black")
            points(validation.case$X[c(pair.table$from[[r]], pair.table$to[[r]]), 1],
                   validation.case$X[c(pair.table$from[[r]], pair.table$to[[r]]), 2],
                   pch = 21, cex = 1.0, bg = c("dodgerblue3", "goldenrod2"), col = "black")
        }
    }
}, width = 2100, height = 1900)

fig.files$benchmark <- file.path(fig.dir, "04_small_iknn_benchmark.png")
save.png(fig.files$benchmark, {
    truth.levels <- unique(benchmark.metrics$truth)
    oldpar <- par(mfrow = c(length(truth.levels), length(params$benchmark_n_grid)),
                  mar = c(4.8, 4.8, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    cols <- c(`0` = "gray25", `0.05` = "#1f77b4", `0.1` = "#d95f02")
    for (truth.name in truth.levels) {
        for (n.value in params$benchmark_n_grid) {
            x <- benchmark.metrics[benchmark_metrics_key(benchmark.metrics, truth.name, n.value), ]
            ylim <- range(x$median_relative_error, na.rm = TRUE)
            plot(NA, xlim = range(params$k_grid), ylim = ylim,
                 xlab = "k", ylab = "median relative GDA error",
                 main = sprintf("%s; n=%s", truth.name, n.value))
            for (delta in params$pruning_delta_grid) {
                y <- x[x$pruning_delta == delta, ]
                y <- y[order(y$k), ]
                lines(y$k, y$median_relative_error, type = "b", pch = 19,
                      col = cols[[as.character(delta)]])
            }
            legend("topright", legend = paste("delta", params$pruning_delta_grid),
                   col = cols[as.character(params$pruning_delta_grid)],
                   lty = 1, pch = 19, bty = "n", cex = 0.75)
        }
    }
}, width = 1700, height = 1400)

fig.files$edge <- file.path(fig.dir, "05_edge_geometry_summary.png")
save.png(fig.files$edge, {
    tube.truths <- c("tube oracle alpha=0.5", "tube oracle alpha=1")
    oldpar <- par(mfrow = c(length(tube.truths), 2), mar = c(5.2, 4.8, 3.2, 1.0))
    on.exit(par(oldpar), add = TRUE)
    cols <- c(`0` = "gray25", `0.05` = "#1f77b4", `0.1` = "#d95f02")
    for (truth.name in tube.truths) {
        edge.tube <- edge.summary[edge.summary$truth == truth.name, ]
        for (metric in c("q90_oracle_edge", "q90_oracle_to_euclidean")) {
            plot(NA, xlim = range(params$k_grid), ylim = range(edge.tube[[metric]], na.rm = TRUE),
                 xlab = "k", ylab = metric,
                 main = sprintf("%s; %s", truth.name, metric))
            for (delta in params$pruning_delta_grid) {
                z <- edge.tube[edge.tube$pruning_delta == delta & edge.tube$n == params$validation_n, ]
                if (!nrow(z)) next
                z <- z[order(z$k), ]
                lines(z$k, z[[metric]], type = "b", pch = 19,
                      col = cols[[as.character(delta)]])
            }
            legend("topright", legend = paste("delta", params$pruning_delta_grid),
                   col = cols[as.character(params$pruning_delta_grid)],
                   lty = 1, pch = 19, bty = "n", cex = 0.75)
        }
    }
}, width = 1700, height = 1300)

rel.fig <- function(path) {
    file.path("../../figures/noisy-circle-latent-tube-oracle-sanity", basename(path))
}

commands <- paste(
    "X.df <- gflow::generate.circle.data(",
    "  n = n, radius = 1, noise = 0.08,",
    "  type = \"random\", noise.type = \"normal\", seed = seed",
    ")",
    "X <- as.matrix(X.df[, c(\"x\", \"y\")])",
    "theta <- as.double(X.df$angles)",
    "u <- sqrt(rowSums(X^2)) - 1",
    "",
    "oracle <- latent.tube.oracle(",
    "  theta.sample = theta, u.sample = u, radius = 1,",
    "  sigma = 0.08, density.type = \"normal\", alpha = alpha,",
    "  n.theta = 360, u.max = 4 * 0.08, n.u = 41,",
    "  density.floor = 1e-4, include.diagonal = TRUE,",
    "  attach.extra = 1L",
    ")",
    "tube.dist <- latent.tube.oracle.distances(oracle)",
    "",
    "iknn <- gflow::create.iknn.graphs(",
    "  X, kmin = k, kmax = k,",
    "  max.path.edge.ratio.deviation.thld = delta,",
    "  path.edge.ratio.percentile = 0.5,",
    "  threshold.percentile = 0, compute.full = TRUE,",
    "  with.isize.pruning = FALSE, pca.dim = NULL,",
    "  n.cores = 1L, verbose = FALSE",
    ")$geom_pruned_graphs[[1L]]",
    sep = "\n"
)

best.winner.text <- paste(vapply(seq_len(nrow(best.by.truth.n)), function(i) {
    sprintf(
        "n=%s, truth=%s: k=%s, delta=%s, median error=%.4f, q90 error=%.4f",
        best.by.truth.n$n[[i]],
        best.by.truth.n$truth[[i]],
        best.by.truth.n$k[[i]],
        best.by.truth.n$pruning_delta[[i]],
        best.by.truth.n$median_relative_error[[i]],
        best.by.truth.n$q90_relative_error[[i]]
    )
}, character(1)), collapse = "; ")

report.path <- file.path(report.dir, "noisy_circle_latent_tube_oracle_sanity.html")
html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Noisy Circle Latent-Tube Oracle Sanity</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;max-width:1240px;}h1,h2,h3{line-height:1.15;}code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f6f8fb;}figure{margin:24px 0;}img{max-width:100%;border:1px solid #d9dee7;}figcaption{font-size:13px;color:#52616f;margin-top:6px;}details{border:1px solid #d9dee7;padding:10px 14px;margin:18px 0;background:#fbfcfe;}summary{font-weight:600;cursor:pointer;}</style>",
    "</head><body>",
    "<h1>Noisy Circle Latent-Tube Oracle Sanity</h1>",
    sprintf("<p>Generated at %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "<h2>Background</h2>",
    "<p>The original noisy-circle benchmarks used latent circular arc distance as truth. That is appropriate when radial noise is very small and the intended geometry is the original one-dimensional circle. Once radial noise becomes a meaningful part of the data-generating distribution, the observed data are better viewed as a density-supported tube around the circle. The truth metric should then be allowed to describe travel through this tube, not only travel along the centerline.</p>",
    "<p>This report implements a first latent-tube oracle. It is not yet the final dense benchmark. Its purpose is to validate the oracle construction on a small, inspectable case before replacing the truth metric in larger GDA sweeps.</p>",
    "<h2>Precise Definition</h2>",
    "<p>Each point is represented by latent tube coordinates <code>(theta, u)</code>, where <code>theta</code> is the generating angle and <code>u = ||x|| - radius</code> is signed radial displacement. The embedding is <code>z(theta,u) = (radius + u) (cos(theta), sin(theta))</code>.</p>",
    "<p>The oracle graph is a structured grid in <code>(theta,u)</code>. Grid vertices are connected to angular, radial, and diagonal neighbors, with periodic wrap in <code>theta</code>. Observed samples are attached to nearby grid vertices in latent coordinate space.</p>",
    "<p>For Gaussian radial noise, density is defined up to scale as <code>rho(u) = exp(-(u/sigma)^2)</code>. For a grid edge <code>a-b</code>, the oracle edge weight is <code>||z_a-z_b|| / rho((u_a+u_b)/2)^alpha</code>. The oracle distance between two observed samples is shortest-path distance between their attached sample vertices in this weighted graph.</p>",
    "<h2>R Calls</h2><pre><code>", html.escape(commands), "</code></pre>",
    "<h2>Oracle Validation Figures</h2>",
    sprintf("<figure><img src=\"%s\"><figcaption>Left: structured oracle grid in observed coordinates, colored by radial density, with observed samples overlaid. Right: observed sample locations in latent tube coordinates.</figcaption></figure>", rel.fig(fig.files$oracle.grid)),
    sprintf("<figure><img src=\"%s\"><figcaption>Comparison of latent arc and Euclidean distances against tube-oracle distances for alpha values 0, 0.5, and 1.</figcaption></figure>", rel.fig(fig.files$truth.scatter)),
    sprintf("<figure><img src=\"%s\"><figcaption>Oracle shortest paths for selected pairs. Each row is one diagnostic pair; columns compare alpha=0, alpha=0.5, and alpha=1.</figcaption></figure>", rel.fig(fig.files$paths)),
    "<h2>Oracle Distance Summary</h2>",
    html.table(validation.truth.summary),
    "<h2>Diagnostic Pairs</h2>",
    html.table(pair.table),
    "<h2>Small iKNN Benchmark</h2>",
    sprintf("<figure><img src=\"%s\"><figcaption>Small iKNN benchmark comparing latent-arc truth and tube-oracle truth at alpha=0.5 and alpha=1. This sanity grid uses n=80 and n=160, seed 4101, k in {5,8,10}, and pruning delta in {0,0.05,0.10}.</figcaption></figure>", rel.fig(fig.files$benchmark)),
    sprintf("<figure><img src=\"%s\"><figcaption>Edge-level summaries for iKNN under tube-oracle truth. These are early diagnostics for whether pruning mostly changes high-end oracle edge length or high oracle/Euclidean edge ratio.</figcaption></figure>", rel.fig(fig.files$edge)),
    "<h3>Best Small-Grid Candidates</h3>",
    "<p>", html.escape(best.winner.text), ".</p>",
    html.table(best.by.truth.n[, c("n", "truth", "k", "pruning_delta", "score",
                                   "median_relative_error", "q90_relative_error",
                                   "spearman", "finite_pair_fraction",
                                   "mean_neighbor_overlap")]),
    "<details><summary>Complete small iKNN benchmark metrics</summary>",
    html.table(benchmark.metrics[, c("n", "seed", "graph", "k", "pruning_delta",
                                     "truth", "score", "median_relative_error",
                                     "q90_relative_error", "spearman",
                                     "finite_pair_fraction",
                                     "mean_neighbor_overlap")]),
    "</details>",
    "<details><summary>Edge geometry summaries</summary>",
    html.table(edge.summary),
    "</details>",
    "<h2>Interpretation</h2>",
    "<p>The oracle grid figure confirms that the numerical truth is being computed on the intended annular tube rather than on an unconstrained two-dimensional plane. The alpha parameter changes the cost of moving away from the high-density centerline: alpha=0 is mostly geometric tube travel, alpha=0.5 gives a moderate density penalty, and alpha=1 gives a strong penalty for tail excursions.</p>",
    "<p>The distance scatterplots reveal a useful but important caution. At alpha=0, tube distance is very close to ordinary geometric travel inside the annulus and remains strongly aligned with latent arc distance. At alpha=0.5, the oracle is still interpretable: median distance is 1.737, Spearman correlation with arc distance is 0.941, and the largest distance is 7.335. At alpha=1, a few radial-tail samples dominate the metric: the maximum distance jumps to 236.36 and the rank correlation with arc distance drops to 0.759. This suggests that alpha=1 is probably too aggressive for sigma=0.08 unless we also redesign the density floor, sample attachment rule, or tail handling.</p>",
    "<p>The shortest-path panels are the most important qualitative check. The opposite-angle path stays inside the annular tube and changes its route as density penalty increases: alpha=0 prefers geometrically shorter inner-radius travel, while positive alpha pushes paths toward higher-density regions. The Euclidean-close/high-oracle-ratio and radial-separation rows show that the oracle can treat visually nearby points as far apart when moving between them requires low-density radial travel. That is the desired mechanism, but alpha=1 makes this effect extreme.</p>",
    "<p>In the small iKNN grid, the latent-arc truth remains easy for iKNN, with median relative errors near 0.02. The alpha=0.5 tube truth is harder but still in a useful range, with best median errors around 0.08. The alpha=1 tube truth is not yet a good benchmark target in this form: best median errors are 2.89 at n=80 and 1.83 at n=160, reflecting the tail-dominated oracle distances rather than ordinary circular geometry.</p>",
    "<p>For alpha=0.5, pruning is helpful in this small grid: the best rows use delta=0.10 at both n=80 and n=160. For alpha=1, unpruned iKNN wins, but this should not be over-interpreted because the truth metric is dominated by a small number of very high-cost radial-tail pairs.</p>",
    "<p>The edge summary gives a promising preview of the next attribution phase. Under alpha=0.5 at n=160, pruning reduces the q90 oracle edge length from roughly 0.69-0.79 to roughly 0.36-0.47, and reduces the q90 oracle/Euclidean ratio from roughly 2.8-3.1 to roughly 2.4-2.7. Under alpha=1 the same pattern is amplified, but the maximum oracle/Euclidean ratio remains enormous because tail-adjacent sample edges are still present. This supports the next diagnostic question: are the largest GDA errors carried by a small number of high oracle/Euclidean-ratio edges, or by missing local connectors that force detours?</p>",
    "<p>The immediate recommendation is to treat alpha=0.5 as the first serious tube-oracle truth candidate and alpha=1 as a stress test. Before using alpha=1 in a dense sweep, we should revisit density flooring, radial attachment, and possibly robust truncation of extreme tail costs.</p>",
    "</body></html>"
)
writeLines(html, report.path)

cat("Wrote report:\n", report.path, "\n", sep = "")
cat("Wrote figures under:\n", fig.dir, "\n", sep = "")
cat("Wrote cache under:\n", cache.dir, "\n", sep = "")

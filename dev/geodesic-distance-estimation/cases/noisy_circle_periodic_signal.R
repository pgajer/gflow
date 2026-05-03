# Noisy Circle: Periodic Signal Recovery CMST Benchmark
#
# This script generates the first MST-completion benchmark report. It is kept
# mostly self-contained on purpose; once a few reports exist, repeated pieces can
# be moved into dev/geodesic-distance-estimation/R/.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

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
        gap <- if (pos == n) {
            theta[ord[1L]] + 2 * pi - theta[i]
        } else {
            theta[j] - theta[i]
        }
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

igraph.from.graph <- function(adj.list, weight.list = NULL) {
    edges <- edge.table(adj.list, weight.list)
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
    g <- igraph.from.graph(adj.list)
    igraph::components(g)$no
}

graph.summary.row <- function(name, adj.list, weight.list, parameter = NA_real_,
                              true.dist = NULL, shortcut.threshold = pi / 4) {
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
        graph = name,
        parameter = parameter,
        n_edges = nrow(edges),
        n_components = component.count(adj.list),
        mean_degree = mean(deg),
        min_degree = min(deg),
        max_degree = max(deg),
        edge_length_median = if (nrow(edges)) stats::median(edges$weight) else NA_real_,
        edge_length_max = if (nrow(edges)) max(edges$weight) else NA_real_,
        false_shortcut_rate = false.shortcut.rate,
        max_edge_true_geodesic = max.edge.true.geodesic
    )
}

edge.key <- function(i, j) {
    paste(pmin(i, j), pmax(i, j), sep = "-")
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

geodesic.metrics <- function(name, adj.list, weight.list, true.dist, parameter = NA_real_,
                             k.neighborhood = 10L) {
    graph.dist <- graph.distances(adj.list, weight.list)
    n <- nrow(true.dist)
    keep <- upper.tri(true.dist)
    finite <- keep & is.finite(graph.dist)
    finite.fraction <- sum(finite) / sum(keep)

    if (!any(finite)) {
        return(data.frame(
            graph = name, parameter = parameter, finite_pair_fraction = finite.fraction,
            pearson = NA_real_, spearman = NA_real_, scale = NA_real_,
            median_relative_error = NA_real_, q90_relative_error = NA_real_,
            mean_neighbor_overlap = NA_real_
        ))
    }

    gd <- graph.dist[finite]
    td <- true.dist[finite]
    scale <- sum(gd * td) / sum(gd^2)
    rel.err <- abs(scale * gd - td) / pmax(td, .Machine$double.eps)

    neighbor.overlap <- vapply(seq_len(n), function(i) {
        true.order <- order(true.dist[i, ], seq_len(n))
        graph.order <- order(graph.dist[i, ], seq_len(n), na.last = NA)
        true.nn <- setdiff(true.order, i)[seq_len(min(k.neighborhood, n - 1L))]
        graph.nn <- setdiff(graph.order, i)[seq_len(min(k.neighborhood, length(graph.order) - 1L))]
        if (!length(graph.nn) || anyNA(graph.nn)) return(NA_real_)
        length(intersect(true.nn, graph.nn)) / length(true.nn)
    }, numeric(1))

    data.frame(
        graph = name,
        parameter = parameter,
        finite_pair_fraction = finite.fraction,
        pearson = suppressWarnings(stats::cor(gd, td, method = "pearson")),
        spearman = suppressWarnings(stats::cor(gd, td, method = "spearman")),
        scale = scale,
        median_relative_error = stats::median(rel.err),
        q90_relative_error = as.numeric(stats::quantile(rel.err, 0.9, names = FALSE)),
        mean_neighbor_overlap = mean(neighbor.overlap, na.rm = TRUE)
    )
}

draw.graph.overlay <- function(X, values, adj.list, weight.list, main, palette.name = "Viridis",
                               max.edges = Inf, highlight.shortcut = NULL) {
    pal <- grDevices::hcl.colors(128, palette.name)
    z <- values
    z01 <- (z - min(z)) / max(diff(range(z)), .Machine$double.eps)
    cols <- pal[pmax(1L, pmin(128L, floor(z01 * 127) + 1L))]
    plot(X[, 1], X[, 2], asp = 1, pch = 19, col = cols, xlab = "x", ylab = "y",
         main = main, cex = 0.65)
    edges <- edge.table(adj.list, weight.list)
    if (nrow(edges)) {
        if (is.finite(max.edges) && nrow(edges) > max.edges) {
            edges <- edges[order(edges$weight, decreasing = TRUE), ][seq_len(max.edges), ]
        }
        edge.cols <- rep("gray35", nrow(edges))
        edge.lwd <- rep(0.9, nrow(edges))
        if (!is.null(highlight.shortcut)) {
            shortcut <- highlight.shortcut(edges)
            edge.cols[shortcut] <- "red3"
            edge.lwd[shortcut] <- 2.2
        }
        segments(X[edges$from, 1], X[edges$from, 2], X[edges$to, 1], X[edges$to, 2],
                 col = edge.cols, lwd = edge.lwd)
        points(X[, 1], X[, 2], pch = 19, col = cols, cex = 0.58)
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
    if (any(keep)) {
        abline(stats::lm(gd[keep] ~ true.dist[keep]), col = "steelblue", lwd = 2)
    }
}

out.root <- file.path(root.dir, "dev/geodesic-distance-estimation")
report.dir <- file.path(out.root, "reports/noisy-circle-periodic-signal")
fig.dir <- file.path(out.root, "figures/noisy-circle-periodic-signal")
cache.dir <- file.path(out.root, "cache/noisy-circle-periodic-signal")
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
    seed = 1002L,
    n = 160L,
    radius = 1,
    geometry_noise = 0.08,
    geometry_noise_type = "normal",
    response_noise_sd = 0.12,
    q_grid = c(0.5, 0.7, 0.9, 0.95, 0.99),
    iknn_k_grid = 3:10,
    shortcut_arc_threshold = pi / 4
)

set.seed(params$seed)
X.df <- generate.circle.data(
    n = params$n,
    radius = params$radius,
    noise = params$geometry_noise,
    type = "random",
    noise.type = params$geometry_noise_type,
    seed = params$seed
)
X <- as.matrix(X.df[, c("x", "y")])
theta <- as.double(X.df$angles)
y.true <- circular.synthetic.mixture.of.gaussians(
    x = theta,
    x.knot = c(0.45, 2.35, 5.35),
    y.knot = c(1.25, 0.85, 1.7),
    sd.knot = 0.32
)
y <- as.double(y.true + stats::rnorm(params$n, sd = params$response_noise_sd))
true.dist <- circle.dist.matrix(theta, radius = params$radius)

oracle <- weighted.circle.graph(theta, radius = params$radius)
oracle.path <- weighted.angular.path.graph(theta, radius = params$radius)

cmst.graphs <- lapply(params$q_grid, function(q) {
    create.cmst.graph(X, q.thld = q, pca.dim = NULL, verbose = FALSE)
})
names(cmst.graphs) <- paste0("q=", params$q_grid)

iknn.res <- create.iknn.graphs(
    X,
    kmin = min(params$iknn_k_grid),
    kmax = max(params$iknn_k_grid),
    max.path.edge.ratio.deviation.thld = 0,
    path.edge.ratio.percentile = 0.5,
    threshold.percentile = 0,
    compute.full = TRUE,
    with.isize.pruning = FALSE,
    pca.dim = NULL,
    n.cores = 1L,
    verbose = FALSE
)
iknn.graphs <- iknn.res$geom_pruned_graphs
names(iknn.graphs) <- paste0("k=", params$iknn_k_grid)

mst.graph <- list(
    adj_list = cmst.graphs[["q=0.9"]]$mst_adj_list,
    weight_list = cmst.graphs[["q=0.9"]]$mst_weight_list
)

graph.entries <- list(
    list(name = "MST only", parameter = NA_real_, adj = mst.graph$adj_list, wt = mst.graph$weight_list),
    list(name = "Oracle weighted path", parameter = NA_real_, adj = oracle.path$adj_list, wt = oracle.path$weight_list),
    list(name = "Oracle weighted circle", parameter = NA_real_, adj = oracle$adj_list, wt = oracle$weight_list)
)
graph.entries <- c(graph.entries, lapply(seq_along(cmst.graphs), function(i) {
    g <- cmst.graphs[[i]]
    list(name = "CMST", parameter = params$q_grid[[i]], adj = g$cmst_adj_list, wt = g$cmst_weight_list)
}))
graph.entries <- c(graph.entries, lapply(seq_along(iknn.graphs), function(i) {
    g <- iknn.graphs[[i]]
    list(name = "iKNN", parameter = params$iknn_k_grid[[i]], adj = g$adj_list, wt = g$weight_list)
}))

graph.summary <- do.call(rbind, lapply(graph.entries, function(g) {
    graph.summary.row(
        g$name, g$adj, g$wt, parameter = g$parameter,
        true.dist = true.dist,
        shortcut.threshold = params$shortcut_arc_threshold
    )
}))

geodesic.summary <- do.call(rbind, lapply(graph.entries, function(g) {
    geodesic.metrics(g$name, g$adj, g$wt, true.dist, parameter = g$parameter)
}))

cmst.thresholds <- data.frame(
    q_thld = params$q_grid,
    cmst_distance_threshold = vapply(cmst.graphs, function(g) g$cmst_distance_threshold, numeric(1)),
    n_edges = vapply(cmst.graphs, function(g) nrow(edge.table(g$cmst_adj_list, g$cmst_weight_list)), integer(1)),
    false_shortcut_rate = graph.summary$false_shortcut_rate[graph.summary$graph == "CMST"],
    spearman_geodesic = geodesic.summary$spearman[geodesic.summary$graph == "CMST"],
    median_relative_error = geodesic.summary$median_relative_error[geodesic.summary$graph == "CMST"]
)

iknn.sweep <- data.frame(
    k = params$iknn_k_grid,
    n_edges = vapply(iknn.graphs, function(g) nrow(edge.table(g$adj_list, g$weight_list)), integer(1)),
    n_components = vapply(iknn.graphs, function(g) component.count(g$adj_list), numeric(1)),
    finite_pair_fraction = geodesic.summary$finite_pair_fraction[geodesic.summary$graph == "iKNN"],
    spearman_geodesic = geodesic.summary$spearman[geodesic.summary$graph == "iKNN"],
    median_relative_error = geodesic.summary$median_relative_error[geodesic.summary$graph == "iKNN"]
)

oracle.retention <- do.call(rbind, lapply(graph.entries, function(g) {
    cbind(
        data.frame(graph = g$name, parameter = g$parameter),
        oracle.edge.retention(g$adj, oracle$adj_list, theta, radius = params$radius)
    )
}))

saveRDS(
    list(params = params, X.df = X.df, X = X, theta = theta, y.true = y.true, y = y,
         true.dist = true.dist, graph.summary = graph.summary, geodesic.summary = geodesic.summary,
         cmst.thresholds = cmst.thresholds, iknn.sweep = iknn.sweep,
         oracle.retention = oracle.retention),
    file.path(cache.dir, "noisy_circle_periodic_signal_results.rds")
)

fig.files <- list()

fig.files$data_geometry <- file.path(fig.dir, "01_data_geometry.png")
save.png(fig.files$data_geometry, {
    oldpar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    draw.graph.overlay(X, theta, oracle$adj_list, oracle$weight_list,
                       "Latent angle with oracle circle")
    draw.graph.overlay(X, y.true, oracle$adj_list, oracle$weight_list,
                       "Periodic truth with oracle circle", palette.name = "Plasma")
    hist(theta, breaks = 24, col = "gray85", border = "white",
         main = "Random angular sampling", xlab = "theta")
}, width = 1800, height = 650)

fig.files$graph_overlays <- file.path(fig.dir, "02_graph_overlays.png")
cmst.q09 <- cmst.graphs[["q=0.9"]]
cmst.q99 <- cmst.graphs[["q=0.99"]]
iknn.k5 <- iknn.graphs[[which(params$iknn_k_grid == 5L)]]
save.png(fig.files$graph_overlays, {
    oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    shortcut.fun <- function(edges) true.dist[cbind(edges$from, edges$to)] > params$shortcut_arc_threshold
    draw.graph.overlay(X, theta, mst.graph$adj_list, mst.graph$weight_list,
                       "MST only", highlight.shortcut = shortcut.fun)
    draw.graph.overlay(X, theta, cmst.q09$cmst_adj_list, cmst.q09$cmst_weight_list,
                       "CMST q=0.90", highlight.shortcut = shortcut.fun)
    draw.graph.overlay(X, theta, cmst.q99$cmst_adj_list, cmst.q99$cmst_weight_list,
                       "CMST q=0.99", highlight.shortcut = shortcut.fun)
    draw.graph.overlay(X, theta, iknn.k5$adj_list, iknn.k5$weight_list,
                       "iKNN k=5", highlight.shortcut = shortcut.fun)
}, width = 1500, height = 1400)

fig.files$distance_scatter <- file.path(fig.dir, "03_distance_scatter.png")
save.png(fig.files$distance_scatter, {
    oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    draw.distance.scatter(mst.graph$adj_list, mst.graph$weight_list, true.dist, "MST only")
    draw.distance.scatter(cmst.q09$cmst_adj_list, cmst.q09$cmst_weight_list, true.dist, "CMST q=0.90")
    draw.distance.scatter(iknn.k5$adj_list, iknn.k5$weight_list, true.dist, "iKNN k=5")
    draw.distance.scatter(oracle$adj_list, oracle$weight_list, true.dist, "Oracle weighted circle")
}, width = 1500, height = 1400)

fig.files$sweeps <- file.path(fig.dir, "04_parameter_sweeps.png")
save.png(fig.files$sweeps, {
    oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    plot(cmst.thresholds$q_thld, cmst.thresholds$spearman_geodesic, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "CMST q.thld", ylab = "Spearman(graph, truth)",
         main = "CMST geodesic recovery")
    plot(cmst.thresholds$q_thld, cmst.thresholds$false_shortcut_rate, type = "b", pch = 19,
         xlab = "CMST q.thld", ylab = "false shortcut rate",
         main = "CMST shortcut diagnostics")
    plot(iknn.sweep$k, iknn.sweep$spearman_geodesic, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "iKNN k", ylab = "Spearman(graph, truth)",
         main = "iKNN geodesic recovery")
    plot(iknn.sweep$k, iknn.sweep$n_components, type = "b", pch = 19,
         xlab = "iKNN k", ylab = "component count",
         main = "iKNN connectedness")
}, width = 1500, height = 1200)

fig.files$response <- file.path(fig.dir, "05_periodic_response.png")
save.png(fig.files$response, {
    ord <- order(theta)
    plot(theta[ord], y[ord], pch = 16, cex = 0.55,
         col = grDevices::adjustcolor("gray30", alpha.f = 0.55),
         xlab = "theta", ylab = "response",
         main = "Periodic response over latent angle")
    lines(theta[ord], y.true[ord], col = "firebrick", lwd = 2)
    points(theta[ord], y[ord], pch = 16, cex = 0.45,
           col = grDevices::adjustcolor("black", alpha.f = 0.35))
    legend("topright", legend = c("truth", "noisy observations"),
           col = c("firebrick", "gray30"), lwd = c(2, NA), pch = c(NA, 16),
           bty = "n")
}, width = 1200, height = 800)

rel.fig <- function(path) {
    file.path("../../figures/noisy-circle-periodic-signal", basename(path))
}

commands <- paste(
    "set.seed(1002L)",
    "X.df <- gflow::generate.circle.data(",
    "  n = 160L, radius = 1, noise = 0.08,",
    "  type = \"random\", noise.type = \"normal\", seed = 1002L",
    ")",
    "X <- as.matrix(X.df[, c(\"x\", \"y\")])",
    "theta <- as.double(X.df$angles)",
    "y.true <- gflow::circular.synthetic.mixture.of.gaussians(",
    "  x = theta, x.knot = c(0.45, 2.35, 5.35),",
    "  y.knot = c(1.25, 0.85, 1.7), sd.knot = 0.32",
    ")",
    "y <- y.true + rnorm(length(theta), sd = 0.12)",
    "true.dist <- radius * pmin(abs(outer(theta, theta, \"-\")),",
    "                           2 * pi - abs(outer(theta, theta, \"-\")))",
    "oracle.circle <- weighted circle graph from observed sorted angles",
    "oracle.path <- weighted angular path graph from observed sorted angles, without cyclic closure",
    "cmst <- lapply(c(0.5, 0.7, 0.9, 0.95, 0.99), function(q) {",
    "  gflow::create.cmst.graph(X, q.thld = q, pca.dim = NULL, verbose = FALSE)",
    "})",
    "iknn <- gflow::create.iknn.graphs(",
    "  X, kmin = 3L, kmax = 10L, compute.full = TRUE,",
    "  max.path.edge.ratio.deviation.thld = 0, threshold.percentile = 0,",
    "  pca.dim = NULL, n.cores = 1L, verbose = FALSE",
    ")",
    sep = "\n"
)

report.path <- file.path(report.dir, "noisy_circle_periodic_signal.html")
html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Noisy Circle CMST Benchmark</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;max-width:1180px;}h1,h2,h3{line-height:1.15;}code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f6f8fb;}figure{margin:24px 0;}img{max-width:100%;border:1px solid #d9dee7;}figcaption{font-size:13px;color:#52616f;margin-top:6px;}.note{background:#fff8e6;border-left:4px solid #e6a700;padding:10px 14px;}</style>",
    "</head><body>",
    "<h1>Noisy Circle: Periodic Signal Recovery</h1>",
    sprintf("<p>Generated at %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "<h2>Why This Example</h2>",
    "<p>This example is the first CMST benchmark because a noisy circle is the smallest setting where connectedness, topology, periodic boundary behavior, local edge completion, and harmful ambient-space shortcuts are all visible. A graph can be connected and still fail by adding chord shortcuts that collapse circular geodesic distances.</p>",
    "<p>The observed data are noisy radial perturbations of points on a unit circle. The response is a smooth periodic mixture of circular Gaussian bumps. The oracle circle graph connects adjacent samples in latent angular order and includes the cyclic closure edge; the oracle path graph uses the same angular order but deliberately omits that closure edge.</p>",
    "<h2>True Geodesic Definition</h2>",
    "<p>Each sample has a latent angle <code>theta_i</code> on a circle of radius <code>r = 1</code>. The true geodesic distance is the circular arc length</p>",
    "<pre><code>d_true(i, j) = r * min(|theta_i - theta_j|, 2*pi - |theta_i - theta_j|)</code></pre>",
    "<p>Radial noise changes the observed coordinates but not the latent geodesic truth. This lets us ask whether each graph reconstructs the circular metric from noisy Euclidean observations.</p>",
    "<h2>Expected Challenges</h2>",
    "<ul><li>Low-k local graphs may be disconnected around sparse angular regions.</li><li>High-threshold graphs may add chord shortcuts across the circle.</li><li>The MST guarantees connectedness but can be too sparse and path-like.</li><li>CMST may fail by preserving a connected path while missing the cyclic closure edge.</li><li>The <code>0 / 2pi</code> boundary should be treated as adjacent, not distant.</li></ul>",
    "<h2>Reproducible Commands</h2>",
    "<pre><code>", html.escape(commands), "</code></pre>",
    "<h2>Data Geometry</h2>",
    sprintf("<figure><img src=\"%s\"><figcaption>Observed noisy circle. Left: points colored by latent angle with oracle cyclic edges. Middle: points colored by periodic truth. Right: angular sampling histogram.</figcaption></figure>", rel.fig(fig.files$data_geometry)),
    sprintf("<figure><img src=\"%s\"><figcaption>Periodic response over latent angle. The true periodic signal is shown in red, with noisy observations overlaid.</figcaption></figure>", rel.fig(fig.files$response)),
    "<h2>Graph Visualizations</h2>",
    "<p>Edges shown in red have true circular arc distance greater than pi/4 and are treated as suspected chord shortcuts for this report. This threshold is a diagnostic, not a mathematical definition of failure.</p>",
    sprintf("<figure><img src=\"%s\"><figcaption>Graph overlays for MST-only, CMST at q=0.90 and q=0.99, and iKNN at k=5. Red edges are suspected shortcuts by true circular distance.</figcaption></figure>", rel.fig(fig.files$graph_overlays)),
    "<h2>Graph Structural Summary</h2>",
    html.table(graph.summary),
    "<h2>Oracle Circle Edge Retention</h2>",
    "<p>This table asks a deliberately simple topological question: how many adjacent-in-angle oracle circle edges are present in each graph? A graph can have many short local edges and still fail to close the circle if it misses enough oracle cyclic adjacencies.</p>",
    html.table(oracle.retention),
    "<h2>Geodesic Recovery</h2>",
    "<p>Graph shortest-path distances are compared to the true circular arc-length distance. Relative errors are computed after a single least-squares global scaling of graph distances to true distances.</p>",
    html.table(geodesic.summary),
    sprintf("<figure><img src=\"%s\"><figcaption>Graph shortest-path distances versus true circular geodesic distances for representative graphs.</figcaption></figure>", rel.fig(fig.files$distance_scatter)),
    "<h2>Parameter Sweeps</h2>",
    "<h3>CMST q.thld Sweep</h3>",
    html.table(cmst.thresholds),
    "<h3>iKNN k Sweep</h3>",
    html.table(iknn.sweep),
    sprintf("<figure><img src=\"%s\"><figcaption>Parameter sweep diagnostics for CMST and iKNN.</figcaption></figure>", rel.fig(fig.files$sweeps)),
    "<h2>Initial Interpretation</h2>",
    "<p>The oracle weighted circle graph is the reference target: it has exactly one cyclic edge per sample and its graph shortest-path metric should closely match the true circular geodesic. The oracle weighted path graph is included as a useful foil: it has the same latent ordering but omits cyclic closure, so it exposes the cost of turning a circle into a path.</p>",
    "<p>In this first run, CMST is connected and adds local support, but its geodesic recovery remains close to MST-only and the oracle path rather than the oracle circle. This is an important arising issue: the completion threshold can add many short local edges without restoring the missing circular closure. The problem here is less about obvious long chord shortcuts and more about topological under-completion.</p>",
    "<p>iKNN behaves differently: low <code>k</code> values can be disconnected, but once connected it recovers the circular geodesic metric very well. This makes noisy circle a strong benchmark for future CMST variants, because any MST-based method needs a way to recognize or repair cycle closure when the true object is periodic.</p>",
    "<h2>Open Issues Raised By This Case</h2>",
    "<ul><li>Should CMST include a diagnostic for topological under-completion, such as missing local cycles or unusually path-like global structure?</li><li>Should completion use a global MST-edge quantile or a local/adaptive threshold on nonuniform circles?</li><li>Should long MST bridges and missing closure edges be flagged separately from completion edges?</li><li>Can a cycle-repair step be added without creating harmful chord shortcuts?</li><li>How should downstream regression balance a graph that has excellent local neighborhoods but biased long-range geodesics?</li></ul>",
    "<h2>Next Benchmark</h2>",
    "<p>The next natural case is a nonuniform noisy circle, because it keeps the same true geodesic definition while stressing the global threshold assumption more directly.</p>",
    "</body></html>"
)
writeLines(html, report.path)

cat("Wrote report:\n", report.path, "\n", sep = "")
cat("Wrote figures under:\n", fig.dir, "\n", sep = "")
cat("Wrote cache under:\n", cache.dir, "\n", sep = "")

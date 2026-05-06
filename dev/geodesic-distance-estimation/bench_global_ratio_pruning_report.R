#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
})

out.dir <- file.path(
  "dev", "geodesic-distance-estimation", "reports",
  "global_ratio_pruning_benchmark"
)
fig.dir <- file.path(out.dir, "figures")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260506)

html.escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

median.kth.distance <- function(X, k) {
  D <- as.matrix(stats::dist(X))
  diag(D) <- Inf
  stats::median(apply(D, 1L, function(row) sort(row, partial = k)[[k]]))
}

make.circle <- function(n, noise = 0.035) {
  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  cbind(cos(theta), sin(theta)) + matrix(stats::rnorm(n * 2L, sd = noise), ncol = 2L)
}

make.swiss <- function(n, noise = 0.02) {
  t <- stats::runif(n, 1.5 * pi, 4.5 * pi)
  h <- stats::runif(n, -0.45, 0.45)
  X <- cbind(t * cos(t), h * 6, t * sin(t)) / (4.5 * pi)
  X + matrix(stats::rnorm(n * 3L, sd = noise), ncol = 3L)
}

make.branches <- function(n, noise = 0.025) {
  arm <- sample.int(3L, n, replace = TRUE)
  r <- stats::runif(n, 0, 1)
  angles <- c(pi / 2, 7 * pi / 6, 11 * pi / 6)
  X <- cbind(r * cos(angles[arm]), r * sin(angles[arm]))
  X + matrix(stats::rnorm(n * 2L, sd = noise), ncol = 2L)
}

cases <- list(
  circle.n500 = list(kind = "noisy circle", n = 500L, k = 12L, X = make.circle(500L)),
  circle.n900 = list(kind = "noisy circle", n = 900L, k = 16L, X = make.circle(900L)),
  swiss.n600 = list(kind = "swiss roll", n = 600L, k = 16L, X = make.swiss(600L)),
  branches.n700 = list(kind = "three branches", n = 700L, k = 14L, X = make.branches(700L))
)

prune.methods <- c("none", "local.geodesic", "global.geodesic.ratio")
families <- c("sknn", "mknn", "radius", "adaptive.radius", "iknn")
bench.reps <- 2L

graph.edge.count <- function(adj.list) {
  sum(lengths(adj.list)) / 2
}

build.raw.graph <- function(case, family) {
  X <- case$X
  k <- case$k
  radius <- median.kth.distance(X, k) * 1.12
  start <- proc.time()[["elapsed"]]
  graph <- switch(
    family,
    sknn = create.sknn.graph(X, k = k, prune.method = "none",
                             connect.components = FALSE),
    mknn = create.mknn.graph(X, k = max(2L, k), prune.method = "none",
                             connect.components = FALSE),
    radius = create.radius.graph(X, radius = radius, prune.method = "none",
                                 connect.components = FALSE),
    adaptive.radius = create.adaptive.radius.graph(
      X, k.scale = k, radius.factor = 1.12, radius.rule = "max",
      prune.method = "none", connect.components = FALSE
    ),
    iknn = create.single.iknn.graph(
      X, k = k, prune.method = "none", connect.components = FALSE,
      with.lifecycle.branches = FALSE, verbose = FALSE
    ),
    stop("unknown family")
  )
  elapsed <- proc.time()[["elapsed"]] - start
  raw.adj <- if (!is.null(graph$raw_adj_list)) graph$raw_adj_list else graph$adj_list
  raw.weight <- if (!is.null(graph$raw_weight_list)) graph$raw_weight_list else graph$weight_list
  prune.k <- if (family %in% c("radius", "adaptive.radius")) {
    .default.radius.prune.k(raw.adj)
  } else {
    k
  }
  list(
    X = X,
    family = family,
    adj_list = raw.adj,
    weight_list = raw.weight,
    build_elapsed_sec = elapsed,
    n_edges = graph.edge.count(raw.adj),
    prune_k = prune.k
  )
}

time.prune <- function(raw.graph, method) {
  gc(verbose = FALSE)
  start <- proc.time()[["elapsed"]]
  out <- .prune.graph.by.method(
    X = raw.graph$X,
    adj.list = raw.graph$adj_list,
    weight.list = raw.graph$weight_list,
    k = raw.graph$prune_k,
    prune.method = method,
    max.path.edge.ratio.deviation.thld = 0.1,
    path.edge.ratio.percentile = 0.5,
    threshold.percentile = 0,
    prune.tau = 1.05,
    prune.local.k = raw.graph$prune_k,
    with.pruned.edge.stats = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - start
  data.frame(
    elapsed_sec = elapsed,
    edges_before = out$n_edges_before_pruning,
    edges_after = out$n_edges_after_pruning,
    pruned_edges = out$n_pruned_edges,
    stringsAsFactors = FALSE
  )
}

message("Building raw graphs and timing pruning methods...")
raw.rows <- list()
bench.rows <- list()
cursor <- 0L
raw.cursor <- 0L
for (case.name in names(cases)) {
  case <- cases[[case.name]]
  for (family in families) {
    message(sprintf("  raw graph: %s / %s", case.name, family))
    raw.graph <- build.raw.graph(case, family)
    raw.cursor <- raw.cursor + 1L
    raw.rows[[raw.cursor]] <- data.frame(
      case = case.name,
      kind = case$kind,
      n = case$n,
      k = case$k,
      family = family,
      raw_edges = raw.graph$n_edges,
      prune_k = raw.graph$prune_k,
      build_elapsed_sec = raw.graph$build_elapsed_sec,
      stringsAsFactors = FALSE
    )
    for (method in prune.methods) {
      for (rep in seq_len(bench.reps)) {
        message(sprintf("    prune: %s rep %d", method, rep))
        timed <- time.prune(raw.graph, method)
        cursor <- cursor + 1L
        bench.rows[[cursor]] <- cbind(
          data.frame(
            case = case.name,
            kind = case$kind,
            n = case$n,
            k = case$k,
            family = family,
            method = method,
            rep = rep,
            stringsAsFactors = FALSE
          ),
          timed
        )
      }
    }
  }
}

raw.graphs <- do.call(rbind, raw.rows)
benchmark <- do.call(rbind, bench.rows)
benchmark$pruned_fraction <- with(benchmark, ifelse(edges_before > 0, pruned_edges / edges_before, 0))
benchmark$edges_per_sec <- with(benchmark, ifelse(elapsed_sec > 0, edges_before / elapsed_sec, NA_real_))
benchmark$pruned_edges_per_sec <- with(benchmark, ifelse(elapsed_sec > 0, pruned_edges / elapsed_sec, NA_real_))

summarize.benchmark <- function(x) {
  split.keys <- interaction(x$case, x$kind, x$n, x$k, x$family, x$method, drop = TRUE)
  pieces <- lapply(split(x, split.keys), function(df) {
    data.frame(
      case = df$case[[1L]],
      kind = df$kind[[1L]],
      n = df$n[[1L]],
      k = df$k[[1L]],
      family = df$family[[1L]],
      method = df$method[[1L]],
      median_elapsed_sec = stats::median(df$elapsed_sec),
      min_elapsed_sec = min(df$elapsed_sec),
      max_elapsed_sec = max(df$elapsed_sec),
      edges_before = df$edges_before[[1L]],
      median_edges_after = stats::median(df$edges_after),
      median_pruned_edges = stats::median(df$pruned_edges),
      median_pruned_fraction = stats::median(df$pruned_fraction),
      median_edges_per_sec = stats::median(df$edges_per_sec, na.rm = TRUE),
      median_pruned_edges_per_sec = stats::median(df$pruned_edges_per_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(out$case, out$family, match(out$method, prune.methods)), , drop = FALSE]
}

benchmark.summary <- summarize.benchmark(benchmark)

run.parity.audit <- function() {
  rows <- list()
  cursor <- 0L
  scenarios <- expand.grid(
    seed = c(31L, 37L, 41L, 43L),
    percentile = c(0, 0.25, 0.5, 0.75, 1),
    ratio = c(1.03, 1.08, 1.12),
    with.stats = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE
  )
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    set.seed(scenario$seed)
    X <- matrix(stats::rnorm(12L * 2L), ncol = 2L)
    graph <- create.sknn.graph(X, k = 5L, prune.method = "none")
    r.result <- .prune.graph.global.geodesic(
      graph$raw_adj_list, graph$raw_weight_list,
      max.ratio.threshold = scenario$ratio,
      path.edge.ratio.percentile = scenario$percentile,
      with.pruned.edge.stats = scenario$with.stats
    )
    cpp.result <- .prune.graph.global.geodesic.ratio(
      graph$raw_adj_list, graph$raw_weight_list,
      max.ratio.threshold = scenario$ratio,
      path.edge.ratio.percentile = scenario$percentile,
      with.pruned.edge.stats = scenario$with.stats
    )
    edge.equal <- isTRUE(all.equal(
      .graph.edge.table(r.result$adj_list, r.result$weight_list),
      .graph.edge.table(cpp.result$adj_list, cpp.result$weight_list),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
    count.equal <- identical(
      c(r.result$n_edges_before_pruning, r.result$n_edges_after_pruning, r.result$n_pruned_edges),
      c(cpp.result$n_edges_before_pruning, cpp.result$n_edges_after_pruning, cpp.result$n_pruned_edges)
    )
    stats.equal <- if (isTRUE(scenario$with.stats)) {
      isTRUE(all.equal(r.result$pruned_edge_stats, cpp.result$pruned_edge_stats,
                       tolerance = 1e-12, check.attributes = FALSE))
    } else {
      nrow(r.result$pruned_edge_stats) == 0L &&
        nrow(cpp.result$pruned_edge_stats) == 0L
    }
    cursor <- cursor + 1L
    rows[[cursor]] <- data.frame(
      seed = scenario$seed,
      percentile = scenario$percentile,
      ratio = scenario$ratio,
      with_stats = scenario$with.stats,
      edge_table_equal = edge.equal,
      count_equal = count.equal,
      stats_equal = stats.equal,
      all_equal = edge.equal && count.equal && stats.equal,
      r_pruned_edges = r.result$n_pruned_edges,
      cpp_pruned_edges = cpp.result$n_pruned_edges,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

message("Running independent parity audit...")
parity <- run.parity.audit()

write.csv(raw.graphs, file.path(out.dir, "raw_graphs.csv"), row.names = FALSE)
write.csv(benchmark, file.path(out.dir, "benchmark_replicates.csv"), row.names = FALSE)
write.csv(benchmark.summary, file.path(out.dir, "benchmark_summary.csv"), row.names = FALSE)
write.csv(parity, file.path(out.dir, "parity_audit.csv"), row.names = FALSE)

png(file.path(fig.dir, "synthetic_data.png"), width = 1200, height = 900, res = 130)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (case.name in names(cases)) {
  X <- cases[[case.name]]$X
  plot(X[, 1], X[, 2], pch = 16, cex = 0.45, col = "#2f6f8f",
       xlab = "coordinate 1", ylab = "coordinate 2",
       main = paste0(case.name, " (n=", cases[[case.name]]$n, ")"))
}
par(op)
dev.off()

png(file.path(fig.dir, "elapsed_by_method.png"), width = 1400, height = 900, res = 130)
plot.df <- benchmark.summary
plot.df$label <- paste(plot.df$case, plot.df$family, sep = "\n")
plot.df <- plot.df[plot.df$method != "none", , drop = FALSE]
labels <- unique(plot.df$label)
ymax <- max(plot.df$median_elapsed_sec, na.rm = TRUE) * 1.18
plot(seq_along(labels), rep(NA_real_, length(labels)), ylim = c(0, ymax),
     xaxt = "n", xlab = "", ylab = "median pruning elapsed time (seconds)",
     main = "Pruning elapsed time by graph and method")
axis(1, at = seq_along(labels), labels = labels, las = 2, cex.axis = 0.62)
cols <- c(local.geodesic = "#2f6f8f", global.geodesic.ratio = "#c4512f")
offsets <- c(local.geodesic = -0.14, global.geodesic.ratio = 0.14)
for (method in names(cols)) {
  y <- plot.df$median_elapsed_sec[plot.df$method == method]
  x <- match(plot.df$label[plot.df$method == method], labels) + offsets[[method]]
  points(x, y, pch = 19, col = cols[[method]], cex = 1.15)
  segments(x, 0, x, y, col = cols[[method]], lwd = 2)
}
legend("topleft", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()

png(file.path(fig.dir, "edge_reduction_by_method.png"), width = 1400, height = 900, res = 130)
plot.df <- benchmark.summary[benchmark.summary$method != "none", , drop = FALSE]
plot.df$label <- paste(plot.df$case, plot.df$family, sep = "\n")
labels <- unique(plot.df$label)
plot(seq_along(labels), rep(NA_real_, length(labels)), ylim = c(0, 1),
     xaxt = "n", xlab = "", ylab = "median fraction of raw edges pruned",
     main = "Edge reduction by graph and method")
axis(1, at = seq_along(labels), labels = labels, las = 2, cex.axis = 0.62)
for (method in names(cols)) {
  y <- plot.df$median_pruned_fraction[plot.df$method == method]
  x <- match(plot.df$label[plot.df$method == method], labels) + offsets[[method]]
  points(x, y, pch = 19, col = cols[[method]], cex = 1.15)
  segments(x, 0, x, y, col = cols[[method]], lwd = 2)
}
legend("topright", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()

png(file.path(fig.dir, "pruning_rate.png"), width = 1200, height = 850, res = 130)
plot.df <- benchmark.summary[benchmark.summary$method != "none", , drop = FALSE]
rate <- pmax(plot.df$median_pruned_edges_per_sec, 0)
plot(plot.df$median_elapsed_sec, rate,
     log = "xy", pch = 19,
     col = cols[plot.df$method],
     xlab = "median elapsed time (seconds, log scale)",
     ylab = "median pruned edges per second (log scale)",
     main = "Runtime and pruning throughput")
text(plot.df$median_elapsed_sec, rate, labels = plot.df$family,
     pos = 4, cex = 0.55, col = "#333333")
legend("topright", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()

png(file.path(fig.dir, "raw_edge_counts.png"), width = 1200, height = 850, res = 130)
raw.graphs$label <- paste(raw.graphs$case, raw.graphs$family, sep = "\n")
barplot(raw.graphs$raw_edges, names.arg = raw.graphs$label, las = 2,
        cex.names = 0.62, col = "#6f8f3a", border = NA,
        ylab = "raw undirected edges",
        main = "Raw graph density before pruning")
dev.off()

png(file.path(fig.dir, "parity_audit.png"), width = 1000, height = 650, res = 130)
parity.counts <- c(
  "edge tables" = sum(parity$edge_table_equal),
  "edge counts" = sum(parity$count_equal),
  "stats" = sum(parity$stats_equal),
  "all checks" = sum(parity$all_equal)
)
barplot(parity.counts, ylim = c(0, nrow(parity)), col = "#2f6f8f", border = NA,
        ylab = sprintf("passing scenarios out of %d", nrow(parity)),
        main = "R-vs-C++ global-ratio parity audit")
abline(h = nrow(parity), lty = 2, col = "#555555")
dev.off()

table.html <- function(df, n = nrow(df)) {
  df <- head(df, n)
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", html.escape(names(df))), collapse = ""), "</tr>")
  rows <- apply(df, 1L, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(as.character(row))), collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste(rows, collapse = "\n"), "</table>")
}

summary.for.html <- benchmark.summary
summary.for.html$median_elapsed_sec <- fmt(summary.for.html$median_elapsed_sec, 4)
summary.for.html$min_elapsed_sec <- fmt(summary.for.html$min_elapsed_sec, 4)
summary.for.html$max_elapsed_sec <- fmt(summary.for.html$max_elapsed_sec, 4)
summary.for.html$median_pruned_fraction <- fmt(summary.for.html$median_pruned_fraction, 3)
summary.for.html$median_edges_per_sec <- fmt(summary.for.html$median_edges_per_sec, 0)
summary.for.html$median_pruned_edges_per_sec <- fmt(summary.for.html$median_pruned_edges_per_sec, 0)

parity.pass <- sum(parity$all_equal)
parity.total <- nrow(parity)
local.rows <- benchmark.summary[benchmark.summary$method == "local.geodesic", ]
global.rows <- benchmark.summary[benchmark.summary$method == "global.geodesic.ratio", ]
paired <- merge(
  local.rows[, c("case", "family", "median_elapsed_sec", "median_pruned_fraction")],
  global.rows[, c("case", "family", "median_elapsed_sec", "median_pruned_fraction")],
  by = c("case", "family"),
  suffixes = c("_local", "_global")
)
median.global.to.local.time <- stats::median(
  paired$median_elapsed_sec_global / paired$median_elapsed_sec_local,
  na.rm = TRUE
)
median.global.minus.local.prune <- stats::median(
  paired$median_pruned_fraction_global - paired$median_pruned_fraction_local,
  na.rm = TRUE
)

html <- c(
  "<!doctype html>",
  "<html lang=\"en\">",
  "<head>",
  "<meta charset=\"utf-8\">",
  "<title>Global Geodesic Ratio Pruning Benchmark and Parity Report</title>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;line-height:1.45;margin:32px;color:#1f2933;max-width:1180px}",
  "h1,h2{color:#16202a}",
  ".note{background:#f3f6f8;border-left:4px solid #2f6f8f;padding:12px 16px;margin:16px 0}",
  "img{max-width:100%;border:1px solid #d7dee5;margin:8px 0 24px 0}",
  "table{border-collapse:collapse;font-size:13px;margin:12px 0 24px 0}",
  "th,td{border:1px solid #d7dee5;padding:6px 8px;text-align:right}",
  "th:first-child,td:first-child,td:nth-child(2),td:nth-child(5),td:nth-child(6){text-align:left}",
  "th{background:#eef3f7}",
  "code{background:#f2f4f6;padding:1px 4px;border-radius:3px}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>Global Geodesic Ratio Pruning Benchmark and Parity Report</h1>",
  sprintf("<p>Generated on %s from the local <code>gflow</code> repository.</p>", Sys.time()),
  "<div class=\"note\">",
  "<strong>Scope.</strong> This report checks the separately named <code>global.geodesic.ratio</code> pruning method added after the local C++ pruning work. It does two things: an R-vs-C++ parity audit against the previous R whole-graph implementation, and a medium synthetic benchmark comparing <code>none</code>, <code>local.geodesic</code>, and <code>global.geodesic.ratio</code> on the same raw graph inputs.",
  "</div>",
  "<h2>What Was Tested</h2>",
  "<p>The parity audit generated small deterministic sKNN graphs and compared the old R helper <code>.prune.graph.global.geodesic()</code> to the new C++ helper behind <code>.prune.graph.global.geodesic.ratio()</code>. The audit varied random seed, ratio threshold, candidate edge percentile, and stats collection. For each scenario it compared final edge tables, pruning counts, and optional per-edge pruning statistics.</p>",
  sprintf("<p><strong>Parity outcome:</strong> %d of %d scenarios passed all checks.</p>", parity.pass, parity.total),
  "<img src=\"figures/parity_audit.png\" alt=\"Parity audit results\">",
  "<h2>Benchmark Design</h2>",
  "<p>The benchmark built raw graphs once per data case and graph family, then timed only the pruning dispatcher on those fixed adjacency and weight lists. This separates pruning cost from graph construction cost. Each pruning method was run twice per raw graph and summarized by the median elapsed time. The global-ratio settings were <code>max.path.edge.ratio.deviation.thld = 0.1</code> and <code>path.edge.ratio.percentile = 0.5</code>; local pruning used <code>prune.tau = 1.05</code>.</p>",
  "<p>The synthetic data cases were noisy circles, a simple Swiss roll projection, and a three-branch structure:</p>",
  "<img src=\"figures/synthetic_data.png\" alt=\"Synthetic benchmark data sets\">",
  "<h2>Raw Graph Density</h2>",
  "<p>Raw edge count is the main context for interpreting pruning runtime. Denser graphs give the whole-graph ratio method more candidate alternatives to inspect and can also give local pruning more local paths to certify.</p>",
  "<img src=\"figures/raw_edge_counts.png\" alt=\"Raw graph edge counts\">",
  "<h2>Runtime Results</h2>",
  sprintf("<p>Across paired local/global cases, the median global-to-local runtime ratio was %s. Values vary by graph family and density, so the plot below is more informative than a single average.</p>", fmt(median.global.to.local.time, 2)),
  "<img src=\"figures/elapsed_by_method.png\" alt=\"Elapsed time by method\">",
  "<h2>Pruning Strength</h2>",
  sprintf("<p>The median global-minus-local edge reduction difference was %s. Positive values mean global-ratio removed a larger fraction of raw edges than local pruning on the same graph.</p>", fmt(median.global.minus.local.prune, 3)),
  "<img src=\"figures/edge_reduction_by_method.png\" alt=\"Edge reduction by method\">",
  "<h2>Throughput View</h2>",
  "<p>This figure combines elapsed time and pruning yield. Points farther upward remove more edges per second; points farther left finish faster. It is useful for seeing whether a method is merely more aggressive or actually efficient for a graph family.</p>",
  "<img src=\"figures/pruning_rate.png\" alt=\"Pruning throughput\">",
  "<h2>Benchmark Summary Table</h2>",
  table.html(summary.for.html),
  "<h2>Generated Artifacts</h2>",
  "<ul>",
  "<li><code>raw_graphs.csv</code>: raw graph sizes and build timings.</li>",
  "<li><code>benchmark_replicates.csv</code>: all timing replicates.</li>",
  "<li><code>benchmark_summary.csv</code>: median benchmark summaries used in this report.</li>",
  "<li><code>parity_audit.csv</code>: R-vs-C++ parity scenarios and outcomes.</li>",
  "</ul>",
  "<h2>Interpretation Notes</h2>",
  "<p>This is a medium synthetic benchmark, not a claim about all production data. The useful outcome is comparative: the C++ global-ratio path has parity with the old R reference on the tested cases, and the report shows where it buys stronger pruning versus where it costs more runtime. If this method becomes a serious default candidate, the next benchmark should add downstream graph-geodesic fidelity metrics, not just pruning runtime and edge counts.</p>",
  "</body>",
  "</html>"
)
writeLines(html, file.path(out.dir, "index.html"))

message("Wrote report: ", normalizePath(file.path(out.dir, "index.html")))

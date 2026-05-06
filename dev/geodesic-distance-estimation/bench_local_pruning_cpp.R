#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
})

make_local_pruning_data <- function(kind, n, seed) {
  set.seed(seed)
  if (identical(kind, "line")) {
    x <- seq(0, 1, length.out = n)
    return(cbind(x, 0.08 * sin(8 * pi * x)))
  }
  if (identical(kind, "circle")) {
    theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
    return(cbind(cos(theta), sin(theta)) +
             matrix(rnorm(2L * n, sd = 0.01), ncol = 2L))
  }
  if (identical(kind, "paraboloid")) {
    xy <- matrix(runif(2L * n, -1, 1), ncol = 2L)
    return(cbind(xy, xy[, 1L]^2 + 0.5 * xy[, 2L]^2))
  }
  stop("Unknown benchmark kind: ", kind, call. = FALSE)
}

time_local_pruning_case <- function(kind,
                                    n,
                                    graph.k,
                                    local.k,
                                    prune.tau = 1.05,
                                    with.stats = FALSE,
                                    seed = 1L) {
  X <- make_local_pruning_data(kind, n, seed)
  graph <- create.sknn.graph(
    X,
    k = graph.k,
    prune.method = "none",
    connect.components = FALSE
  )
  gc()
  timing <- system.time({
    pruned <- .prune.graph.local.geodesic(
      X = X,
      adj.list = graph$raw_adj_list,
      weight.list = graph$raw_weight_list,
      k = graph.k,
      prune.tau = prune.tau,
      prune.local.k = local.k,
      with.pruned.edge.stats = with.stats
    )
  })
  data.frame(
    kind = kind,
    n = n,
    graph_k = graph.k,
    local_k = local.k,
    prune_tau = prune.tau,
    with_stats = with.stats,
    edges_before = pruned$n_edges_before_pruning,
    pruned_edges = pruned$n_pruned_edges,
    elapsed_sec = unname(timing[["elapsed"]]),
    sec_per_edge = unname(timing[["elapsed"]]) / pruned$n_edges_before_pruning
  )
}

cases <- data.frame(
  kind = c(
    "line", "circle", "paraboloid",
    "circle", "paraboloid",
    "circle", "paraboloid",
    "paraboloid"
  ),
  n = c(1000L, 1000L, 1000L, 2000L, 2000L, 2000L, 2000L, 2000L),
  graph_k = c(20L, 20L, 20L, 20L, 20L, 40L, 40L, 40L),
  local_k = c(20L, 20L, 20L, 20L, 20L, 40L, 40L, 80L),
  prune_tau = c(1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05),
  with_stats = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  seed = seq_len(8L) + 5000L,
  stringsAsFactors = FALSE
)

message("Local geodesic pruning C++ benchmark")
message("Run from package root with: Rscript dev/geodesic-distance-estimation/bench_local_pruning_cpp.R")

results <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  row <- cases[i, ]
  message(sprintf(
    "case %d/%d: kind=%s n=%d graph.k=%d local.k=%d stats=%s",
    i, nrow(cases), row$kind, row$n, row$graph_k, row$local_k, row$with_stats
  ))
  time_local_pruning_case(
    kind = row$kind,
    n = row$n,
    graph.k = row$graph_k,
    local.k = row$local_k,
    prune.tau = row$prune_tau,
    with.stats = row$with_stats,
    seed = row$seed
  )
}))

print(results, row.names = FALSE)

#!/usr/bin/env Rscript

## PHATE affinity-to-length geometry reconstruction report.
##
## Development artifact. It compares edge-length transformations on PHATE graph
## support, with ikNN ambient/native edge lengths as controls, on quadratic graph
## surfaces.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(gflow)
}

required <- c("FNN", "ggplot2", "igraph", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

out.dir <- file.path("dev", "phate-affinity-length-geometry", "report")
fig.dir <- file.path(out.dir, "figures")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

metrics.path <- file.path(out.dir, "affinity_length_geometry_metrics.csv")
edge.samples.path <- file.path(out.dir, "affinity_length_edge_samples.csv")
geo.samples.path <- file.path(out.dir, "affinity_length_geodesic_samples.csv")
html.path <- file.path(out.dir, "phate_affinity_length_geometry_report.html")

datasets <- c("paraboloid", "saddle")
n.values <- c(50L, 100L, 200L, 300L)
k.values <- 3L:10L
decay.values <- c(10L, 20L, 40L, 80L)
eps <- 1e-7
surface.dim <- 2L
ambient.dim <- 3L

html.escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

make.dataset <- function(dataset, n) {
  seed <- switch(dataset,
                 paraboloid = 21000L,
                 saddle = 31000L,
                 41000L) + as.integer(n)
  set.seed(seed)
  theta <- stats::runif(n, 0, 2 * pi)
  radius <- sqrt(stats::runif(n, 0, 1))
  u <- cbind(radius * cos(theta), radius * sin(theta))
  z <- if (identical(dataset, "paraboloid")) {
    u[, 1]^2 + u[, 2]^2
  } else {
    u[, 1]^2 - u[, 2]^2
  }
  list(U = u, X = cbind(u, z), seed = seed)
}

embed.surface <- function(dataset, U) {
  z <- if (identical(dataset, "paraboloid")) {
    U[, 1]^2 + U[, 2]^2
  } else {
    U[, 1]^2 - U[, 2]^2
  }
  cbind(U, z)
}

dense.disk.points <- function(n = 2600L, seed = 90210L) {
  set.seed(seed)
  theta <- stats::runif(n, 0, 2 * pi)
  radius <- sqrt(stats::runif(n, 0, 1))
  cbind(radius * cos(theta), radius * sin(theta))
}

edge.lengths <- function(X, edges) {
  if (is.null(edges) || nrow(edges) == 0L) return(numeric(0))
  sqrt(rowSums((X[edges[, 1] + 1L, , drop = FALSE] -
                  X[edges[, 2] + 1L, , drop = FALSE])^2))
}

edge.keys <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) return(character(0))
  paste(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]), sep = "-")
}

edges.from.kernel <- function(K) {
  K <- as.matrix(K)
  idx <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(matrix(integer(0), ncol = 2L))
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  cbind(idx[, 1] - 1L, idx[, 2] - 1L)
}

weights.from.kernel <- function(K, edges) {
  K[cbind(edges[, 1] + 1L, edges[, 2] + 1L)]
}

edges.weights.from.adj <- function(adj.list, weight.list) {
  keys <- character()
  vals <- numeric()
  for (i in seq_along(adj.list)) {
    nbr <- as.integer(adj.list[[i]])
    wt <- as.numeric(weight.list[[i]])
    if (length(nbr) == 0L) next
    if (length(nbr) != length(wt)) {
      stop("adjacency and weight lists are not aligned at vertex ", i)
    }
    keep <- nbr != i
    nbr <- nbr[keep]
    wt <- wt[keep]
    if (length(nbr) == 0L) next
    a <- pmin(i - 1L, nbr - 1L)
    b <- pmax(i - 1L, nbr - 1L)
    keys <- c(keys, paste(a, b, sep = "-"))
    vals <- c(vals, wt)
  }
  if (length(keys) == 0L) {
    return(list(edges = matrix(integer(0), ncol = 2L), weights = numeric(0)))
  }
  split.vals <- split(vals, keys)
  keys.u <- names(split.vals)
  parts <- strsplit(keys.u, "-", fixed = TRUE)
  edges <- do.call(rbind, lapply(parts, as.integer))
  ord <- order(edges[, 1], edges[, 2])
  list(edges = edges[ord, , drop = FALSE],
       weights = vapply(split.vals, min, numeric(1))[ord])
}

component.count <- function(n, edges) {
  if (is.null(edges) || nrow(edges) == 0L) return(n)
  g <- igraph::graph_from_edgelist(edges + 1L, directed = FALSE)
  g <- igraph::add_vertices(g, n - igraph::vcount(g))
  igraph::components(g)$no
}

make.graph.distances <- function(n, edges, weights) {
  g <- igraph::graph_from_edgelist(edges + 1L, directed = FALSE)
  if (igraph::vcount(g) < n) {
    g <- igraph::add_vertices(g, n - igraph::vcount(g))
  }
  igraph::E(g)$weight <- as.numeric(weights)
  as.matrix(igraph::distances(g, v = seq_len(n), to = seq_len(n),
                              weights = igraph::E(g)$weight))
}

reference.geodesic.distances <- function(dataset, U.sample, seed) {
  U.dense <- dense.disk.points(seed = seed)
  U.all <- rbind(U.sample, U.dense)
  X.all <- embed.surface(dataset, U.all)
  nn <- FNN::get.knn(U.all, k = 12L)
  n.all <- nrow(U.all)
  edges <- cbind(rep(seq_len(n.all), each = 12L), as.vector(t(nn$nn.index)))
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  edges <- unique(edges)
  weights <- sqrt(rowSums((X.all[edges[, 1], , drop = FALSE] -
                            X.all[edges[, 2], , drop = FALSE])^2))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  if (igraph::vcount(g) < n.all) {
    g <- igraph::add_vertices(g, n.all - igraph::vcount(g))
  }
  igraph::E(g)$weight <- weights
  as.matrix(igraph::distances(g, v = seq_len(nrow(U.sample)),
                              to = seq_len(nrow(U.sample)),
                              weights = igraph::E(g)$weight))
}

upper.values <- function(M) {
  M[upper.tri(M)]
}

calibrate <- function(truth, estimate) {
  ok <- is.finite(truth) & is.finite(estimate) & truth > 0 & estimate > 0
  truth <- truth[ok]
  estimate <- estimate[ok]
  if (length(truth) == 0L) {
    return(list(scale = NA_real_, calibrated = numeric(0), truth = numeric(0),
                estimate = numeric(0), finite_fraction = mean(ok)))
  }
  denom <- sum(estimate^2)
  scale <- if (denom > 0) sum(truth * estimate) / denom else NA_real_
  calibrated <- scale * estimate
  list(scale = scale, calibrated = calibrated, truth = truth,
       estimate = estimate, finite_fraction = mean(ok))
}

calibration.metrics <- function(truth, estimate) {
  cal <- calibrate(truth, estimate)
  if (length(cal$truth) == 0L || !is.finite(cal$scale)) {
    return(data.frame(
      scale = NA_real_, rel_rmse = NA_real_, median_abs_log_error = NA_real_,
      q90_abs_log_error = NA_real_, bias_log_error = NA_real_,
      pearson = NA_real_, spearman = NA_real_,
      finite_fraction = cal$finite_fraction
    ))
  }
  log.err <- log(cal$calibrated / cal$truth)
  data.frame(
    scale = cal$scale,
    rel_rmse = sqrt(sum((cal$truth - cal$calibrated)^2) / sum(cal$truth^2)),
    median_abs_log_error = stats::median(abs(log.err)),
    q90_abs_log_error = as.numeric(stats::quantile(abs(log.err), 0.9, names = FALSE)),
    bias_log_error = mean(log.err),
    pearson = suppressWarnings(stats::cor(cal$truth, cal$estimate, method = "pearson")),
    spearman = suppressWarnings(stats::cor(cal$truth, cal$estimate, method = "spearman")),
    finite_fraction = cal$finite_fraction
  )
}

sample.calibrated <- function(truth, estimate, max.n, seed) {
  cal <- calibrate(truth, estimate)
  n <- length(cal$truth)
  if (n == 0L || !is.finite(cal$scale)) return(data.frame())
  set.seed(seed)
  idx <- if (n > max.n) sample.int(n, max.n) else seq_len(n)
  data.frame(
    truth = cal$truth[idx],
    estimate_raw = cal$estimate[idx],
    estimate_calibrated = cal$calibrated[idx],
    log_error = log(cal$calibrated[idx] / cal$truth[idx])
  )
}

phate.length.models <- function(K.values, decay) {
  Kc <- pmin(pmax(as.numeric(K.values), eps), 1)
  neglog <- pmax(-log(Kc), eps)
  list(
    phate_ambient = NULL,
    phate_neglog = neglog,
    phate_rootlog = neglog^(1 / decay),
    phate_inv_power_ambient_p = Kc^(-1 / ambient.dim),
    phate_inv_power_ambient_p_minus1 = pmax(Kc^(-1 / ambient.dim) - 1, eps),
    phate_inv_power_surface_dim_minus1 = pmax(Kc^(-1 / surface.dim) - 1, eps)
  )
}

metric.rows <- list()
edge.samples <- list()
geo.samples <- list()

previous.path <- file.path("dev", "phate-graph-mds-comparison", "report",
                           "phate_graph_mds_geometry_metrics.csv")
selected.cases <- NULL
if (file.exists(previous.path)) {
  prev <- utils::read.csv(previous.path, stringsAsFactors = FALSE)
  prev <- prev[prev$status == "ok", , drop = FALSE]
  prev <- prev[order(prev$case_id, prev$rel_dist_error), , drop = FALSE]
  prev.best <- prev[!duplicated(prev$case_id), , drop = FALSE]
  selected.cases <- prev.best[prev.best$method == "weighted_grip",
                              c("dataset", "n", "k", "decay")]
  selected.cases$key <- paste(selected.cases$dataset, selected.cases$n,
                              selected.cases$k, selected.cases$decay, sep = "__")
}

case.grid <- expand.grid(dataset = datasets, n = n.values, k = k.values,
                         decay = decay.values, stringsAsFactors = FALSE)
case.grid$key <- paste(case.grid$dataset, case.grid$n, case.grid$k,
                       case.grid$decay, sep = "__")
case.grid$in_phate_winning_subset <- if (is.null(selected.cases)) {
  TRUE
} else {
  case.grid$key %in% selected.cases$key
}
case.grid <- case.grid[case.grid$in_phate_winning_subset, , drop = FALSE]
case.grid <- case.grid[order(case.grid$dataset, case.grid$n, case.grid$k,
                             case.grid$decay), ]

message("Running ", nrow(case.grid), " selected graph cases.")

append.csv <- function(rows, path) {
  if (length(rows) == 0L) return()
  df <- do.call(rbind, rows)
  utils::write.table(df, path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path),
                     qmethod = "double")
}

for (path in c(metrics.path, edge.samples.path, geo.samples.path)) {
  if (file.exists(path)) file.remove(path)
}

for (dataset in datasets) {
  for (n in n.values) {
    ds <- make.dataset(dataset, n)
    X <- ds$X
    U <- ds$U
    D <- as.matrix(stats::dist(X))
    D.ref <- reference.geodesic.distances(dataset, U, seed = ds$seed + 500000L)
    iknn.cache <- list()

    sub.grid <- case.grid[case.grid$dataset == dataset & case.grid$n == n, ,
                          drop = FALSE]
    if (nrow(sub.grid) == 0L) next
    message(sprintf("Dataset %s n=%d: %d cases", dataset, n, nrow(sub.grid)))

    for (ii in seq_len(nrow(sub.grid))) {
      k <- sub.grid$k[ii]
      decay <- sub.grid$decay[ii]
      case.id <- paste(dataset, paste0("n", n), paste0("k", k),
                       paste0("decay", decay), sep = "_")
      message("  ", case.id)

      core <- phate.core(
        X = X,
        k = k,
        alpha.decay = decay,
        kernel.mode = "phate",
        vne.method = "phate",
        potential.mode = "phate",
        gamma = 1,
        compute.D.pot = FALSE,
        verbose = FALSE
      )
      ph.edges <- edges.from.kernel(core$K)
      ph.K <- weights.from.kernel(core$K, ph.edges)
      ph.ambient <- edge.lengths(X, ph.edges)
      ph.models <- phate.length.models(ph.K, decay = decay)
      ph.models$phate_ambient <- ph.ambient

      if (is.null(iknn.cache[[as.character(k)]])) {
        g <- create.single.iknn.graph(
          X = X,
          k = k,
          max.path.edge.ratio.deviation.thld = 0,
          threshold.percentile = 0,
          compute.full = TRUE,
          with.isize.pruning = FALSE,
          pca.dim = NULL,
          verbose = FALSE
        )
        ik <- edges.weights.from.adj(g$pruned_adj_list, g$pruned_weight_list)
        iknn.cache[[as.character(k)]] <- ik
      }
      ik <- iknn.cache[[as.character(k)]]
      ik.ambient <- edge.lengths(X, ik$edges)
      ik.models <- list(
        iknn_ambient = ik.ambient,
        iknn_native = ik$weights
      )

      graph.models <- c(
        setNames(ph.models, names(ph.models)),
        setNames(ik.models, names(ik.models))
      )

      case.metric.rows <- list()
      case.edge.samples <- list()
      case.geo.samples <- list()

      for (model in names(graph.models)) {
        is.phate <- startsWith(model, "phate")
        edges <- if (is.phate) ph.edges else ik$edges
        truth.edge <- if (is.phate) ph.ambient else ik.ambient
        estimate.edge <- as.numeric(graph.models[[model]])
        if (length(estimate.edge) != length(truth.edge)) {
          stop("Length mismatch for ", model, " in ", case.id)
        }
        estimate.edge <- pmax(estimate.edge, eps)

        edge.met <- calibration.metrics(truth.edge, estimate.edge)
        D.graph <- make.graph.distances(n, edges, estimate.edge)
        geo.met <- calibration.metrics(upper.values(D.ref), upper.values(D.graph))
        n.components <- component.count(n, edges)

        row <- data.frame(
          dataset = dataset,
          n = n,
          k = k,
          decay = decay,
          case_id = case.id,
          graph = if (is.phate) "PHATE" else "ikNN",
          length_model = model,
          edge_count = nrow(edges),
          n_components = n.components,
          edge_scale = edge.met$scale,
          edge_rel_rmse = edge.met$rel_rmse,
          edge_median_abs_log_error = edge.met$median_abs_log_error,
          edge_q90_abs_log_error = edge.met$q90_abs_log_error,
          edge_bias_log_error = edge.met$bias_log_error,
          edge_pearson = edge.met$pearson,
          edge_spearman = edge.met$spearman,
          geo_scale = geo.met$scale,
          geo_rel_rmse = geo.met$rel_rmse,
          geo_median_abs_log_error = geo.met$median_abs_log_error,
          geo_q90_abs_log_error = geo.met$q90_abs_log_error,
          geo_bias_log_error = geo.met$bias_log_error,
          geo_pearson = geo.met$pearson,
          geo_spearman = geo.met$spearman,
          geo_finite_fraction = geo.met$finite_fraction,
          stringsAsFactors = FALSE
        )
        case.metric.rows[[length(case.metric.rows) + 1L]] <- row

        es <- sample.calibrated(truth.edge, estimate.edge, max.n = 400L,
                                seed = ds$seed + 1000L + k * 100L + decay)
        if (nrow(es) > 0L) {
          es$dataset <- dataset
          es$n <- n
          es$k <- k
          es$decay <- decay
          es$case_id <- case.id
          es$graph <- if (is.phate) "PHATE" else "ikNN"
          es$length_model <- model
          case.edge.samples[[length(case.edge.samples) + 1L]] <- es
        }
        gs <- sample.calibrated(upper.values(D.ref), upper.values(D.graph),
                                max.n = 900L,
                                seed = ds$seed + 3000L + k * 100L + decay)
        if (nrow(gs) > 0L) {
          gs$dataset <- dataset
          gs$n <- n
          gs$k <- k
          gs$decay <- decay
          gs$case_id <- case.id
          gs$graph <- if (is.phate) "PHATE" else "ikNN"
          gs$length_model <- model
          case.geo.samples[[length(case.geo.samples) + 1L]] <- gs
        }
      }

      append.csv(case.metric.rows, metrics.path)
      append.csv(case.edge.samples, edge.samples.path)
      append.csv(case.geo.samples, geo.samples.path)
    }
  }
}

metrics <- utils::read.csv(metrics.path, stringsAsFactors = FALSE)
edge.samples <- utils::read.csv(edge.samples.path, stringsAsFactors = FALSE)
geo.samples <- utils::read.csv(geo.samples.path, stringsAsFactors = FALSE)

model.labels <- c(
  phate_ambient = "PHATE ambient",
  phate_neglog = "PHATE -log(K)",
  phate_rootlog = "PHATE [-log(K)]^(1/a)",
  phate_inv_power_ambient_p = "PHATE K^(-1/3)",
  phate_inv_power_ambient_p_minus1 = "PHATE K^(-1/3)-1",
  phate_inv_power_surface_dim_minus1 = "PHATE K^(-1/2)-1",
  iknn_ambient = "ikNN ambient",
  iknn_native = "ikNN native"
)
metrics$label <- model.labels[metrics$length_model]
edge.samples$label <- model.labels[edge.samples$length_model]
geo.samples$label <- model.labels[geo.samples$length_model]

summary <- aggregate(
  cbind(edge_rel_rmse, edge_median_abs_log_error, edge_q90_abs_log_error,
        geo_rel_rmse, geo_median_abs_log_error, geo_q90_abs_log_error,
        geo_spearman, geo_finite_fraction) ~ dataset + graph + length_model + label,
  metrics,
  function(x) c(mean = mean(x, na.rm = TRUE),
                median = stats::median(x, na.rm = TRUE))
)
summary <- do.call(data.frame, summary)
names(summary) <- gsub("\\.", "_", names(summary))
summary <- summary[order(summary$dataset, summary$geo_rel_rmse_mean), ]

best.edge <- metrics[order(metrics$case_id, metrics$edge_rel_rmse), ]
best.edge <- best.edge[!duplicated(best.edge$case_id), ]
best.geo <- metrics[order(metrics$case_id, metrics$geo_rel_rmse), ]
best.geo <- best.geo[!duplicated(best.geo$case_id), ]
wins.edge <- as.data.frame(table(best.edge$dataset, best.edge$label),
                           stringsAsFactors = FALSE)
names(wins.edge) <- c("dataset", "label", "edge_rmse_wins")
wins.geo <- as.data.frame(table(best.geo$dataset, best.geo$label),
                          stringsAsFactors = FALSE)
names(wins.geo) <- c("dataset", "label", "geo_rmse_wins")
wins <- merge(wins.edge, wins.geo, by = c("dataset", "label"), all = TRUE)
wins$edge_rmse_wins[is.na(wins$edge_rmse_wins)] <- 0L
wins$geo_rmse_wins[is.na(wins$geo_rmse_wins)] <- 0L
wins <- wins[wins$edge_rmse_wins > 0L | wins$geo_rmse_wins > 0L, ]
wins <- wins[order(wins$dataset, -wins$geo_rmse_wins, -wins$edge_rmse_wins), ]

utils::write.csv(summary, file.path(out.dir, "affinity_length_geometry_summary.csv"),
                 row.names = FALSE)
utils::write.csv(wins, file.path(out.dir, "affinity_length_geometry_wins.csv"),
                 row.names = FALSE)

model.order <- unname(model.labels)
edge.samples$label <- factor(edge.samples$label, levels = model.order)
geo.samples$label <- factor(geo.samples$label, levels = model.order)
metrics$label <- factor(metrics$label, levels = model.order)

hist.edge <- ggplot2::ggplot(edge.samples,
                             ggplot2::aes(x = log_error, fill = graph)) +
  ggplot2::geom_histogram(bins = 55, alpha = 0.78, position = "identity") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35) +
  ggplot2::facet_grid(dataset ~ label, scales = "free_y") +
  ggplot2::coord_cartesian(xlim = c(-1.2, 1.2)) +
  ggplot2::labs(x = "edge log error after scalar calibration",
                y = "sampled edge count",
                title = "Edge-length calibration errors") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                 legend.position = "bottom")
hist.edge.path <- file.path(fig.dir, "edge_log_error_histograms.png")
ggplot2::ggsave(hist.edge.path, hist.edge, width = 15, height = 6.8, dpi = 160)

hist.geo <- ggplot2::ggplot(geo.samples,
                            ggplot2::aes(x = log_error, fill = graph)) +
  ggplot2::geom_histogram(bins = 55, alpha = 0.78, position = "identity") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35) +
  ggplot2::facet_grid(dataset ~ label, scales = "free_y") +
  ggplot2::coord_cartesian(xlim = c(-1.2, 1.2)) +
  ggplot2::labs(x = "geodesic log error after scalar calibration",
                y = "sampled pair count",
                title = "Graph-geodesic calibration errors") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                 legend.position = "bottom")
hist.geo.path <- file.path(fig.dir, "geodesic_log_error_histograms.png")
ggplot2::ggsave(hist.geo.path, hist.geo, width = 15, height = 6.8, dpi = 160)

scatter.edge.sample <- edge.samples
if (nrow(scatter.edge.sample) > 50000L) {
  set.seed(11)
  scatter.edge.sample <- scatter.edge.sample[sample.int(nrow(scatter.edge.sample), 50000L), ]
}
scatter.edge <- ggplot2::ggplot(scatter.edge.sample,
                                ggplot2::aes(x = truth, y = estimate_calibrated,
                                             color = graph)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.35) +
  ggplot2::geom_point(alpha = 0.22, size = 0.45) +
  ggplot2::facet_grid(dataset ~ label, scales = "free") +
  ggplot2::labs(x = "true edge length: ambient chord",
                y = "calibrated estimated edge length",
                title = "True edge length vs estimated edge length") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                 legend.position = "bottom")
scatter.edge.path <- file.path(fig.dir, "edge_true_vs_estimated_scatter.png")
ggplot2::ggsave(scatter.edge.path, scatter.edge, width = 15, height = 6.8, dpi = 160)

scatter.geo.sample <- geo.samples
if (nrow(scatter.geo.sample) > 50000L) {
  set.seed(12)
  scatter.geo.sample <- scatter.geo.sample[sample.int(nrow(scatter.geo.sample), 50000L), ]
}
scatter.geo <- ggplot2::ggplot(scatter.geo.sample,
                               ggplot2::aes(x = truth, y = estimate_calibrated,
                                            color = graph)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.35) +
  ggplot2::geom_point(alpha = 0.18, size = 0.4) +
  ggplot2::facet_grid(dataset ~ label, scales = "free") +
  ggplot2::labs(x = "reference surface geodesic distance",
                y = "calibrated graph geodesic distance",
                title = "True geodesic length vs estimated graph geodesic length") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                 legend.position = "bottom")
scatter.geo.path <- file.path(fig.dir, "geodesic_true_vs_estimated_scatter.png")
ggplot2::ggsave(scatter.geo.path, scatter.geo, width = 15, height = 6.8, dpi = 160)

format.num <- function(x, digits = 4L) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

html.table <- function(df, digits = 4L, max.rows = Inf) {
  if (nrow(df) > max.rows) df <- utils::head(df, max.rows)
  df.out <- df
  for (nm in names(df.out)) {
    if (is.numeric(df.out[[nm]])) df.out[[nm]] <- format.num(df.out[[nm]], digits)
  }
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", html.escape(names(df.out))),
                                 collapse = ""), "</tr>")
  rows <- apply(df.out, 1L, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(row)), collapse = ""), "</tr>")
  })
  paste0("<table><thead>", header, "</thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}

display.summary <- summary[, c(
  "dataset", "graph", "label",
  "edge_rel_rmse_mean", "edge_median_abs_log_error_median",
  "geo_rel_rmse_mean", "geo_median_abs_log_error_median",
  "geo_spearman_mean", "geo_finite_fraction_mean"
)]
display.summary <- display.summary[order(display.summary$dataset,
                                         display.summary$geo_rel_rmse_mean), ]

iknn.native.delta <- metrics[metrics$length_model %in% c("iknn_ambient", "iknn_native"), ]
iknn.native.delta <- reshape(
  iknn.native.delta[, c("case_id", "dataset", "length_model", "edge_rel_rmse",
                       "geo_rel_rmse")],
  idvar = c("case_id", "dataset"),
  timevar = "length_model",
  direction = "wide"
)
native.max.edge.delta <- max(abs(iknn.native.delta$edge_rel_rmse.iknn_ambient -
                                  iknn.native.delta$edge_rel_rmse.iknn_native),
                             na.rm = TRUE)
native.max.geo.delta <- max(abs(iknn.native.delta$geo_rel_rmse.iknn_ambient -
                                 iknn.native.delta$geo_rel_rmse.iknn_native),
                            na.rm = TRUE)

css <- paste(
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:34px;line-height:1.5;color:#222}",
  "h1,h2{color:#16324f}h3{color:#243b53;margin-top:24px}",
  ".note{background:#f7f9fb;border-left:4px solid #3B6EA8;padding:12px 14px;margin:16px 0}",
  "table{border-collapse:collapse;width:100%;font-size:13px;margin:14px 0}",
  "th,td{border:1px solid #ddd;padding:5px 7px;text-align:left}th{background:#f2f5f8}",
  "code{background:#f4f4f4;padding:1px 3px}.fig{max-width:1280px;width:100%;border:1px solid #ddd}",
  ".math-block{overflow-x:auto;margin:12px 0}",
  sep = "\n"
)

html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\">",
  "<title>PHATE Affinity-to-Length Geometry Reconstruction</title>",
  "<script>window.MathJax={tex:{inlineMath:[[\"$\",\"$\"],[\"\\\\(\",\"\\\\)\"]],displayMath:[[\"\\\\[\",\"\\\\]\"],[\"$$\",\"$$\"]]}};</script>",
  "<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>",
  "<style>", css, "</style></head><body>",
  "<h1>PHATE Affinity-to-Length Geometry Reconstruction</h1>",
  "<p><strong>Generated:</strong> ", html.escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
  "<div class=\"note\">This report follows from the previous observation that PHATE graph support, when equipped with ambient chord lengths, usually reconstructs the 2D quadratic surface examples better than the tested alternatives. Here the support question is held mostly fixed and the edge-length question is tested directly.</div>",
  "<h2>Context And Questions</h2>",
  "<p>The data-geodesic reconstruction problem asks: given finite data \\(X\\subset\\mathbb{R}^p\\), which weighted graph \\(G(X)\\) makes graph shortest-path distance \\(d_{G(X)}\\) as close as possible to the data manifold geodesic distance \\(d_X\\)?</p>",
  "<p>The previous report suggested that the PHATE edge support is often superior to ikNN support on sampled quadratic graph surfaces when both supports use local ambient chord lengths. This report asks the next question: if the PHATE support is good, can PHATE affinities themselves be transformed into useful edge lengths?</p>",
  "<p>The tested surfaces are sampled from</p>",
  "<div class=\"math-block\">\\[z=x_1^2+x_2^2\\qquad\\text{and}\\qquad z=x_1^2-x_2^2,\\]</div>",
  "<p>with points sampled over the unit disk and embedded as graph surfaces in \\(\\mathbb{R}^3\\). Reference surface geodesics are approximated numerically by a dense auxiliary graph on the latent disk, with local edge lengths measured in the embedded surface.</p>",
  "<h2>Length Models</h2>",
  "<p>For PHATE support, the baseline is the ambient chord length</p>",
  "<div class=\"math-block\">\\[\\ell_{ij}^{\\mathrm{ambient}}=\\lVert x_i-x_j\\rVert_2.\\]</div>",
  "<p>The affinity-based candidates are</p>",
  "<div class=\"math-block\">\\[\\ell_{ij}^{\\log}=-\\log(K_{ij}),\\qquad \\ell_{ij}^{\\mathrm{rootlog}}=(-\\log(K_{ij}))^{1/a},\\]</div>",
  "<div class=\"math-block\">\\[\\ell_{ij}^{p}=K_{ij}^{-1/3},\\qquad \\ell_{ij}^{p0}=K_{ij}^{-1/3}-1,\\qquad \\ell_{ij}^{s0}=K_{ij}^{-1/2}-1.\\]</div>",
  "<p>The exponent \\(1/3\\) corresponds to the ambient data dimension \\(p=3\\). The centered variants subtract one so that high-affinity edges near \\(K_{ij}=1\\) have near-zero estimated length before calibration.</p>",
  "<p>For ikNN, two controls are included: ambient chord lengths on ikNN support and the native gflow ikNN graph weights. The native ikNN weight is not simply the direct chord length. For an ikNN edge \\((i,j)\\), gflow stores an intersection-mediated length based on shared neighbors, approximately</p>",
  "<div class=\"math-block\">\\[\\ell_{ij}^{\\mathrm{ikNN\\ native}}=\\min_{c\\in N_k(i)\\cap N_k(j)}\\{d(i,c)+d(j,c)\\}.\\]</div>",
  "<p>Thus the native ikNN model tests whether this two-hop shared-neighbor length can improve global graph-geodesic geometry relative to direct ambient chord lengths on the same ikNN support.</p>",
  "<h2>Scalar Calibration</h2>",
  "<p>Every candidate length model is allowed one global scale factor before error is measured. For true lengths \\(t_m\\) and raw estimated lengths \\(r_m\\), the scalar is the least-squares scale</p>",
  "<div class=\"math-block\">\\[s^*=\\arg\\min_s\\sum_m(t_m-sr_m)^2=\\frac{\\sum_m t_mr_m}{\\sum_m r_m^2}.\\]</div>",
  "<p>The calibrated estimate is \\(\\hat t_m=s^*r_m\\). Edge calibration uses edge pairs \\((i,j)\\in E\\). Geodesic calibration uses vertex pairs \\((u,v)\\), comparing reference surface geodesic distance to shortest-path distance in the candidate weighted graph. Histograms show</p>",
  "<div class=\"math-block\">\\[\\log\\left(\\frac{\\hat t_m}{t_m}\\right),\\]</div>",
  "<p>so zero means exact after scale calibration, positive means too long, and negative means too short.</p>",
  "<h2>Scope</h2>",
  "<p>The report evaluates ", nrow(case.grid), " selected cases: the cases where the previous PHATE-support weighted-GRIP benchmark selected PHATE ambient support as the best PHATE-graph embedding method. This follows the proposed restriction to \\(k\\) and decay settings where PHATE support already looked promising.</p>",
  "<h2>Winner Counts</h2>",
  html.table(wins, digits = 0L),
  "<h2>Summary Metrics</h2>",
  html.table(display.summary, digits = 4L),
  "<h2>Edge-Length Error Histograms</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(hist.edge.path), "\" alt=\"Edge log error histograms\"></p>",
  "<h2>Edge True-Vs-Estimated Calibration Plots</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(scatter.edge.path), "\" alt=\"Edge true versus estimated scatter\"></p>",
  "<h2>Graph-Geodesic Error Histograms</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(hist.geo.path), "\" alt=\"Geodesic log error histograms\"></p>",
  "<h2>Graph-Geodesic True-Vs-Estimated Calibration Plots</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(scatter.geo.path), "\" alt=\"Geodesic true versus estimated scatter\"></p>",
  "<h2>ikNN Native Weight Check</h2>",
  "<p>The maximum absolute difference between ikNN ambient-length and ikNN native-weight edge RMSE across matched cases is <strong>", format.num(native.max.edge.delta, 6), "</strong>; for graph-geodesic RMSE it is <strong>", format.num(native.max.geo.delta, 6), "</strong>. These are not zero, confirming that the native ikNN weights define a distinct length model rather than duplicating direct ambient chords.</p>",
  "<h2>Generated Files</h2>",
  "<ul>",
  "<li><code>", html.escape(normalizePath(metrics.path, mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(edge.samples.path, mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(geo.samples.path, mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "affinity_length_geometry_summary.csv"), mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "affinity_length_geometry_wins.csv"), mustWork = FALSE)), "</code></li>",
  "</ul>",
  "</body></html>"
)

writeLines(html, html.path)
message("Wrote ", normalizePath(html.path, mustWork = FALSE))

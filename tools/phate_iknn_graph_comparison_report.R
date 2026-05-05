#!/usr/bin/env Rscript

## Interactive PHATE-vs-ikNN graph comparison report.
## This is a development/report artifact; it does not modify package APIs.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(gflow)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required to generate this report.")
}
if (!requireNamespace("grip", quietly = TRUE)) {
  stop("Package 'grip' is required to generate weighted GRIP embeddings.")
}

out.dir <- file.path("dev", "phate-graph-comparison", "report")
lib.dir <- file.path(out.dir, "lib")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lib.dir, recursive = TRUE, showWarnings = FALSE)

plotly.bundle <- system.file(
  "htmlwidgets", "lib", "plotlyjs", "plotly-latest.min.js",
  package = "plotly"
)
if (!nzchar(plotly.bundle) || !file.exists(plotly.bundle)) {
  stop("Could not locate the installed plotly JavaScript bundle.")
}
file.copy(plotly.bundle, file.path(lib.dir, "plotly-latest.min.js"), overwrite = TRUE)

datasets <- c("noisy_circle", "paraboloid", "saddle")
n.values <- c(50L, 100L, 200L, 300L)
k.values <- 3L:10L
decay.values <- c(10L, 20L, 40L, 80L)
iknn.candidate.k <- 1L:30L

html.escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

dataset.label <- function(dataset) {
  switch(dataset,
         noisy_circle = "Noisy Circle",
         paraboloid = "2D Paraboloid Graph",
         saddle = "2D Saddle Graph",
         dataset)
}

dataset.description <- function(dataset) {
  switch(dataset,
         noisy_circle = "A noisy circle embedded in three dimensions.",
         paraboloid = "A quadratic-form graph z = x^2 + y^2 sampled over a disk.",
         saddle = "A quadratic-form graph z = x^2 - y^2 sampled over a disk.",
         "")
}

make.dataset <- function(dataset, n) {
  seed <- switch(dataset,
                 noisy_circle = 11000L,
                 paraboloid = 21000L,
                 saddle = 31000L,
                 41000L) + as.integer(n)
  set.seed(seed)

  if (identical(dataset, "noisy_circle")) {
    theta <- sort(stats::runif(n, 0, 2 * pi))
    radius <- 1 + stats::rnorm(n, sd = 0.055)
    X <- cbind(
      radius * cos(theta) + stats::rnorm(n, sd = 0.025),
      radius * sin(theta) + stats::rnorm(n, sd = 0.025),
      stats::rnorm(n, sd = 0.075)
    )
    color.value <- theta
  } else {
    theta <- stats::runif(n, 0, 2 * pi)
    radius <- sqrt(stats::runif(n, 0, 1))
    x <- radius * cos(theta)
    y <- radius * sin(theta)
    z <- if (identical(dataset, "paraboloid")) x^2 + y^2 else x^2 - y^2
    X <- cbind(x, y, z)
    color.value <- z
  }

  list(X = X, color.value = color.value, seed = seed)
}

scale.coords <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X)
  X <- sweep(X, 2L, center, "-")
  scale <- max(sqrt(rowSums(X^2)))
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  X / scale
}

point.colors <- function(x) {
  r <- rank(x, ties.method = "first")
  grDevices::hcl.colors(length(x), "Viridis")[r]
}

edges.from.kernel <- function(K) {
  K <- as.matrix(K)
  idx <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(list(edges = matrix(integer(0), ncol = 2), weights = numeric(0)))
  }
  ord <- order(idx[, 1], idx[, 2])
  idx <- idx[ord, , drop = FALSE]
  weights <- K[idx]
  list(edges = cbind(idx[, 1] - 1L, idx[, 2] - 1L), weights = as.numeric(weights))
}

edges.from.adj <- function(adj.list) {
  out <- list()
  m <- 0L
  for (i in seq_along(adj.list)) {
    nbr <- as.integer(adj.list[[i]])
    if (length(nbr) == 0L) next
    for (j in nbr) {
      if (!is.na(j) && j > i) {
        m <- m + 1L
        out[[m]] <- c(i - 1L, j - 1L)
      }
    }
  }
  if (m == 0L) {
    matrix(integer(0), ncol = 2)
  } else {
    do.call(rbind, out)
  }
}

edge.keys <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(character(0))
  }
  paste(edges[, 1], edges[, 2], sep = "-")
}

edge.pairs.for.json <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(list())
  }
  lapply(seq_len(nrow(edges)), function(i) as.integer(edges[i, ]))
}

component.stats <- function(n, edges) {
  parent <- seq_len(n)
  size <- rep.int(1L, n)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  unite <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return()
    if (size[ra] < size[rb]) {
      tmp <- ra
      ra <- rb
      rb <- tmp
    }
    parent[rb] <<- ra
    size[ra] <<- size[ra] + size[rb]
  }
  if (!is.null(edges) && nrow(edges) > 0L) {
    for (i in seq_len(nrow(edges))) {
      unite(edges[i, 1] + 1L, edges[i, 2] + 1L)
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  tab <- table(roots)
  list(
    n_components = length(tab),
    largest_component_fraction = as.numeric(max(tab)) / n
  )
}

edge.lengths <- function(X, edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(numeric(0))
  }
  sqrt(rowSums((X[edges[, 1] + 1L, , drop = FALSE] -
                  X[edges[, 2] + 1L, , drop = FALSE])^2))
}

graph.metrics <- function(X, edges, weights = NULL, weighted.degree = NULL) {
  n <- nrow(X)
  deg <- integer(n)
  if (!is.null(edges) && nrow(edges) > 0L) {
    deg <- tabulate(c(edges[, 1], edges[, 2]) + 1L, nbins = n)
  }
  comps <- component.stats(n, edges)
  lengths <- edge.lengths(X, edges)
  q <- if (length(lengths) > 0L) {
    as.numeric(stats::quantile(lengths, probs = c(0.1, 0.5, 0.9), names = FALSE))
  } else {
    c(NA_real_, NA_real_, NA_real_)
  }
  out <- list(
    edge_count = if (is.null(edges)) 0L else nrow(edges),
    mean_degree = mean(deg),
    median_degree = stats::median(deg),
    max_degree = max(deg),
    n_components = comps$n_components,
    largest_component_fraction = comps$largest_component_fraction,
    edge_length_q10 = q[1],
    edge_length_q50 = q[2],
    edge_length_q90 = q[3]
  )
  if (!is.null(weighted.degree)) {
    out$weighted_degree_mean <- mean(weighted.degree)
    out$weighted_degree_q10 <- as.numeric(stats::quantile(weighted.degree, 0.1, names = FALSE))
    out$weighted_degree_q50 <- as.numeric(stats::quantile(weighted.degree, 0.5, names = FALSE))
    out$weighted_degree_q90 <- as.numeric(stats::quantile(weighted.degree, 0.9, names = FALSE))
  }
  out
}

compare.edge.sets <- function(phate.edges, iknn.edges) {
  p <- edge.keys(phate.edges)
  i <- edge.keys(iknn.edges)
  inter <- intersect(p, i)
  union <- union(p, i)
  list(
    shared = length(inter),
    phate_only = length(setdiff(p, i)),
    iknn_only = length(setdiff(i, p)),
    jaccard = if (length(union) == 0L) 1 else length(inter) / length(union),
    phate_only_fraction = if (length(p) == 0L) 0 else length(setdiff(p, i)) / length(p),
    iknn_only_fraction = if (length(i) == 0L) 0 else length(setdiff(i, p)) / length(i)
  )
}

coords.record <- function(X, color) {
  Xs <- scale.coords(X)
  lapply(seq_len(nrow(Xs)), function(i) {
    list(
      x = unname(Xs[i, 1]),
      y = unname(Xs[i, 2]),
      z = unname(Xs[i, 3]),
      color = unname(color[i])
    )
  })
}

coords.record.from.scaled <- function(Xs, color) {
  lapply(seq_len(nrow(Xs)), function(i) {
    list(
      x = unname(Xs[i, 1]),
      y = unname(Xs[i, 2]),
      z = unname(Xs[i, 3]),
      color = unname(color[i])
    )
  })
}

embedding.record <- function(Y, color, label, kind, diagnostics = list()) {
  list(
    label = label,
    kind = kind,
    points = coords.record(Y, color),
    diagnostics = diagnostics
  )
}

procrustes.overlay.record <- function(target, source, target.label, source.label) {
  target <- as.matrix(target)
  source <- as.matrix(source)
  target.centered <- sweep(target, 2L, colMeans(target), "-")
  source.centered <- sweep(source, 2L, colMeans(source), "-")
  sv <- svd(t(source.centered) %*% target.centered)
  R <- sv$u %*% t(sv$v)
  source.rot <- source.centered %*% R
  denom <- sum(source.rot^2)
  scale <- if (is.finite(denom) && denom > 0) {
    sum(target.centered * source.rot) / denom
  } else {
    1
  }
  source.aligned <- scale * source.rot
  residual <- sqrt(sum((target.centered - source.aligned)^2) / sum(target.centered^2))
  combined <- scale.coords(rbind(target.centered, source.aligned))
  n <- nrow(target)
  list(
    target_label = target.label,
    source_label = source.label,
    residual = residual,
    target = coords.record.from.scaled(combined[seq_len(n), , drop = FALSE],
                                       rep("#111111", n)),
    source = coords.record.from.scaled(combined[n + seq_len(n), , drop = FALSE],
                                       rep("#D7191C", n))
  )
}

embedding.edge.weights <- function(X, edges) {
  edge.lengths(X, edges)
}

weighted.grip.embedding <- function(X, edges, seed) {
  n <- nrow(X)
  if (is.null(edges) || nrow(edges) == 0L) {
    stop("Cannot compute weighted GRIP embedding for graph with no edges.")
  }
  weights <- embedding.edge.weights(X, edges)
  grip::grip.layout.weighted(
    edges = edges + 1L,
    edge_weights = weights,
    n = n,
    dim = 3L,
    rounds = 24L,
    final_rounds = 32L,
    num_init = min(8L, n),
    num_nbrs = min(16L, n - 1L),
    coarse_repulsion_factor = 0.8,
    coarse_repulsion_sample = min(12L, n),
    coarse_repulsion_exact_below = 48L,
    seed = seed,
    disconnected = "components"
  )
}

phate.support.kernel <- function(D, edges, k, decay, threshold = 1e-4) {
  D <- as.matrix(D)
  n <- nrow(D)
  k <- min(as.integer(k), n - 1L)
  sigma <- numeric(n)
  for (i in seq_len(n)) {
    row <- D[i, ]
    row[i] <- Inf
    sigma[i] <- max(sort(row, partial = k)[k], .Machine$double.eps)
  }

  A <- matrix(0, nrow = n, ncol = n)
  diag(A) <- 1
  if (!is.null(edges) && nrow(edges) > 0L) {
    for (r in seq_len(nrow(edges))) {
      i <- edges[r, 1] + 1L
      j <- edges[r, 2] + 1L
      vij <- exp(-((D[i, j] / sigma[i])^decay))
      vji <- exp(-((D[j, i] / sigma[j])^decay))
      if (is.finite(vij) && vij >= threshold) A[i, j] <- vij
      if (is.finite(vji) && vji >= threshold) A[j, i] <- vji
    }
  }
  K <- 0.5 * (A + t(A))
  rs <- rowSums(K)
  P <- K / rs
  P[!is.finite(P)] <- 0
  list(K = K, P = P, sigma = sigma)
}

phate.embedding.from.P <- function(P, seed, t.max = 30L) {
  core <- phate.core(
    P = P,
    t = "auto",
    t.max = t.max,
    vne.method = "phate",
    potential.mode = "phate",
    gamma = 1,
    compute.D.pot = TRUE,
    verbose = FALSE
  )
  emb <- phate.embed(
    core = core,
    ndim = 3L,
    method = "classic",
    seed = seed,
    verbose = FALSE
  )$embedding
  list(embedding = emb, t = core$t)
}

graph.records <- list()
point.records <- list()
embedding.records <- list()
embedding.matrices <- list()
overlay.records <- list()
case.records <- list()
summary.rows <- list()

message("Generating PHATE-vs-ikNN graph comparison grid...")

for (dataset in datasets) {
  for (n in n.values) {
    ds <- make.dataset(dataset, n)
    X <- ds$X
    D <- as.matrix(stats::dist(X))
    color <- point.colors(ds$color.value)
    point.id <- paste(dataset, n, sep = "__")
    point.records[[point.id]] <- list(
      id = point.id,
      dataset = dataset,
      dataset_label = dataset.label(dataset),
      description = dataset.description(dataset),
      n = n,
      seed = ds$seed,
      points = coords.record(X, color)
    )

    message(sprintf("  Dataset %s, n=%d: ikNN candidates", dataset, n))
    iknn.by.k <- list()
    max.candidate <- min(max(iknn.candidate.k), n - 1L)
    for (ik in seq_len(max.candidate)) {
      g <- create.single.iknn.graph(
        X = X,
        k = ik,
        max.path.edge.ratio.deviation.thld = 0,
        threshold.percentile = 0,
        compute.full = FALSE,
        with.isize.pruning = FALSE,
        pca.dim = NULL,
        verbose = FALSE
      )
      edges <- edges.from.adj(g$pruned_adj_list)
      graph.id <- paste("iknn", dataset, n, ik, sep = "__")
      metrics <- graph.metrics(X, edges)
      graph.records[[graph.id]] <- list(
        id = graph.id,
        type = "iknn",
        point_id = point.id,
        dataset = dataset,
        n = n,
        k = ik,
        decay = NA,
        edges = edge.pairs.for.json(edges),
        metrics = metrics
      )
      grip.id <- paste("grip", graph.id, sep = "__")
      grip.Y <- weighted.grip.embedding(X, edges, seed = ds$seed + 5000L + ik)
      embedding.records[[grip.id]] <- embedding.record(
        grip.Y,
        color = color,
        label = sprintf("Weighted GRIP: ikNN k=%d", ik),
        kind = "weighted_grip",
        diagnostics = list(edge_length = "ambient_euclidean")
      )
      embedding.matrices[[grip.id]] <- grip.Y
      graph.records[[graph.id]]$grip_embedding <- grip.id
      iknn.by.k[[as.character(ik)]] <- list(id = graph.id, edges = edges, metrics = metrics)
    }

    for (k in k.values) {
      for (decay in decay.values) {
        message(sprintf("  Dataset %s, n=%d: PHATE k=%d decay=%d", dataset, n, k, decay))
        fit <- phate.core(
          X = X,
          k = k,
          alpha.decay = decay,
          t = 1,
          kernel.mode = "phate",
          vne.method = "phate",
          potential.mode = "phate",
          compute.D.pot = FALSE,
          verbose = FALSE
        )
        ph <- edges.from.kernel(fit$K)
        graph.id <- paste("phate", dataset, n, k, decay, sep = "__")
        metrics <- graph.metrics(
          X = X,
          edges = ph$edges,
          weights = ph$weights,
          weighted.degree = rowSums(fit$K)
        )
        graph.records[[graph.id]] <- list(
          id = graph.id,
          type = "phate",
          point_id = point.id,
          dataset = dataset,
          n = n,
          k = k,
          decay = decay,
          edges = edge.pairs.for.json(ph$edges),
          metrics = metrics
        )
        grip.phate.id <- paste("grip", graph.id, sep = "__")
        grip.phate.Y <- weighted.grip.embedding(X, ph$edges, seed = ds$seed + 7000L + 100L * k + decay)
        embedding.records[[grip.phate.id]] <- embedding.record(
          grip.phate.Y,
          color = color,
          label = sprintf("Weighted GRIP: PHATE k=%d decay=%d", k, decay),
          kind = "weighted_grip",
          diagnostics = list(edge_length = "ambient_euclidean")
        )
        embedding.matrices[[grip.phate.id]] <- grip.phate.Y
        graph.records[[graph.id]]$grip_embedding <- grip.phate.id

        diffusion.phate.id <- paste("diffusion", graph.id, sep = "__")
        diffusion.phate <- phate.embedding.from.P(
          P = fit$P,
          seed = ds$seed + 9000L + 100L * k + decay
        )
        embedding.records[[diffusion.phate.id]] <- embedding.record(
          diffusion.phate$embedding,
          color = color,
          label = sprintf("PHATE embedding: PHATE graph k=%d decay=%d", k, decay),
          kind = "phate_diffusion_embedding",
          diagnostics = list(t = diffusion.phate$t, method = "classic")
        )
        embedding.matrices[[diffusion.phate.id]] <- diffusion.phate$embedding
        graph.records[[graph.id]]$diffusion_embedding <- diffusion.phate.id

        same <- iknn.by.k[[as.character(k)]]
        edge.counts <- vapply(iknn.by.k, function(x) x$metrics$edge_count, numeric(1))
        candidate.ks <- as.integer(names(iknn.by.k))
        diff.edge.count <- abs(edge.counts - metrics$edge_count)
        best.idx <- order(diff.edge.count, candidate.ks)[1]
        matched <- iknn.by.k[[best.idx]]
        matched.k <- candidate.ks[best.idx]

        cmp.same <- compare.edge.sets(ph$edges, same$edges)
        cmp.matched <- compare.edge.sets(ph$edges, matched$edges)

        case.id <- paste(dataset, n, k, decay, sep = "__")
        diffusion.iknn.same.id <- paste("diffusion__iknn_same", case.id, sep = "__")
        iknn.same.kernel <- phate.support.kernel(D, same$edges, k = k, decay = decay)
        diffusion.iknn.same <- phate.embedding.from.P(
          P = iknn.same.kernel$P,
          seed = ds$seed + 11000L + 100L * k + decay
        )
        embedding.records[[diffusion.iknn.same.id]] <- embedding.record(
          diffusion.iknn.same$embedding,
          color = color,
          label = sprintf("PHATE embedding: ikNN same-k support, k=%d decay=%d", k, decay),
          kind = "phate_diffusion_embedding",
          diagnostics = list(t = diffusion.iknn.same$t, method = "classic",
                             support = "iknn_same_k")
        )
        embedding.matrices[[diffusion.iknn.same.id]] <- diffusion.iknn.same$embedding

        diffusion.iknn.matched.id <- paste("diffusion__iknn_matched", case.id, sep = "__")
        iknn.matched.kernel <- phate.support.kernel(D, matched$edges, k = k, decay = decay)
        diffusion.iknn.matched <- phate.embedding.from.P(
          P = iknn.matched.kernel$P,
          seed = ds$seed + 13000L + 100L * k + decay
        )
        embedding.records[[diffusion.iknn.matched.id]] <- embedding.record(
          diffusion.iknn.matched$embedding,
          color = color,
          label = sprintf("PHATE embedding: ikNN density-matched support, k=%d decay=%d", matched.k, decay),
          kind = "phate_diffusion_embedding",
          diagnostics = list(t = diffusion.iknn.matched$t, method = "classic",
                             support = "iknn_density_matched")
        )
        embedding.matrices[[diffusion.iknn.matched.id]] <- diffusion.iknn.matched$embedding

        grip.overlay.same.id <- paste("overlay__grip_same", case.id, sep = "__")
        overlay.records[[grip.overlay.same.id]] <- procrustes.overlay.record(
          target = embedding.matrices[[grip.phate.id]],
          source = embedding.matrices[[graph.records[[same$id]]$grip_embedding]],
          target.label = "PHATE graph weighted GRIP",
          source.label = "ikNN same-k weighted GRIP"
        )
        grip.overlay.matched.id <- paste("overlay__grip_matched", case.id, sep = "__")
        overlay.records[[grip.overlay.matched.id]] <- procrustes.overlay.record(
          target = embedding.matrices[[grip.phate.id]],
          source = embedding.matrices[[graph.records[[matched$id]]$grip_embedding]],
          target.label = "PHATE graph weighted GRIP",
          source.label = "ikNN density-matched weighted GRIP"
        )
        diffusion.overlay.same.id <- paste("overlay__diffusion_same", case.id, sep = "__")
        overlay.records[[diffusion.overlay.same.id]] <- procrustes.overlay.record(
          target = embedding.matrices[[diffusion.phate.id]],
          source = embedding.matrices[[diffusion.iknn.same.id]],
          target.label = "PHATE graph PHATE embedding",
          source.label = "ikNN same-k PHATE embedding"
        )
        diffusion.overlay.matched.id <- paste("overlay__diffusion_matched", case.id, sep = "__")
        overlay.records[[diffusion.overlay.matched.id]] <- procrustes.overlay.record(
          target = embedding.matrices[[diffusion.phate.id]],
          source = embedding.matrices[[diffusion.iknn.matched.id]],
          target.label = "PHATE graph PHATE embedding",
          source.label = "ikNN density-matched PHATE embedding"
        )

        case.records[[case.id]] <- list(
          id = case.id,
          point_id = point.id,
          dataset = dataset,
          n = n,
          k = k,
          decay = decay,
          phate_graph = graph.id,
          iknn_same_graph = same$id,
          iknn_matched_graph = matched$id,
          grip_phate_embedding = grip.phate.id,
          grip_iknn_same_embedding = graph.records[[same$id]]$grip_embedding,
          grip_iknn_matched_embedding = graph.records[[matched$id]]$grip_embedding,
          diffusion_phate_embedding = diffusion.phate.id,
          diffusion_iknn_same_embedding = diffusion.iknn.same.id,
          diffusion_iknn_matched_embedding = diffusion.iknn.matched.id,
          grip_same_overlay = grip.overlay.same.id,
          grip_matched_overlay = grip.overlay.matched.id,
          diffusion_same_overlay = diffusion.overlay.same.id,
          diffusion_matched_overlay = diffusion.overlay.matched.id,
          iknn_matched_k = matched.k,
          same = cmp.same,
          matched = cmp.matched
        )
        summary.rows[[length(summary.rows) + 1L]] <- data.frame(
          dataset = dataset,
          n = n,
          k = k,
          decay = decay,
          phate_edges = metrics$edge_count,
          iknn_same_k = k,
          iknn_same_edges = same$metrics$edge_count,
          iknn_same_jaccard = cmp.same$jaccard,
          iknn_matched_k = matched.k,
          iknn_matched_edges = matched$metrics$edge_count,
          iknn_matched_jaccard = cmp.matched$jaccard,
          phate_components = metrics$n_components,
          iknn_same_components = same$metrics$n_components,
          iknn_matched_components = matched$metrics$n_components,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

summary.df <- do.call(rbind, summary.rows)
write.csv(summary.df, file.path(out.dir, "phate_iknn_graph_comparison_summary.csv"),
          row.names = FALSE)

payload <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  datasets = point.records,
  graphs = graph.records,
  embeddings = embedding.records,
  overlays = overlay.records,
  cases = case.records,
  controls = list(
    dataset = datasets,
    n = n.values,
    k = k.values,
    decay = decay.values,
    iknn_mode = c("same", "matched"),
    edge_mode = c("none", "own", "both", "difference")
  )
)

payload.json <- jsonlite::toJSON(
  payload,
  auto_unbox = TRUE,
  null = "null",
  digits = 7,
  dataframe = "rows"
)
writeLines(payload.json, file.path(out.dir, "phate_iknn_graph_payload.json"))

summary.table <- utils::head(
  summary.df[order(summary.df$iknn_same_jaccard), ],
  12L
)
summary.table.html <- paste0(
  "<table><thead><tr>",
  paste(sprintf("<th>%s</th>", html.escape(names(summary.table))), collapse = ""),
  "</tr></thead><tbody>",
  paste(apply(summary.table, 1, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(row)), collapse = ""), "</tr>")
  }), collapse = "\n"),
  "</tbody></table>"
)

js <- '
const DATA = __PAYLOAD__;
const datasetSel = document.getElementById("dataset");
const nSel = document.getElementById("n");
const kSel = document.getElementById("k");
const decaySel = document.getElementById("decay");
const iknnModeSel = document.getElementById("iknnMode");
const edgeModeSel = document.getElementById("edgeMode");
const metricsBox = document.getElementById("metricsBox");

function optionList(select, values, labels) {
  select.innerHTML = "";
  values.forEach((v) => {
    const opt = document.createElement("option");
    opt.value = v;
    opt.textContent = labels && labels[v] ? labels[v] : v;
    select.appendChild(opt);
  });
}

function fmt(x, digits=3) {
  if (x === null || x === undefined || Number.isNaN(Number(x))) return "NA";
  return Number(x).toFixed(digits);
}

function currentCase() {
  const id = [datasetSel.value, nSel.value, kSel.value, decaySel.value].join("__");
  return DATA.cases[id];
}

function pointArrays(pointId) {
  const pts = DATA.datasets[pointId].points;
  return recordArrays(pts);
}

function recordArrays(pts) {
  return {
    x: pts.map(p => p.x),
    y: pts.map(p => p.y),
    z: pts.map(p => p.z),
    color: pts.map(p => p.color)
  };
}

function embeddingArrays(embeddingId) {
  return recordArrays(DATA.embeddings[embeddingId].points);
}

function overlayArrays(overlayId) {
  const overlay = DATA.overlays[overlayId];
  return {
    target: recordArrays(overlay.target),
    source: recordArrays(overlay.source),
    residual: overlay.residual,
    targetLabel: overlay.target_label,
    sourceLabel: overlay.source_label
  };
}

function edgeKey(edge) {
  return edge[0] + "-" + edge[1];
}

function edgeSet(edges) {
  const out = new Set();
  edges.forEach(e => out.add(edgeKey(e)));
  return out;
}

function edgeTrace(points, edges, color, name, width=2, opacity=0.7) {
  const x = [], y = [], z = [];
  edges.forEach((e) => {
    const i = e[0], j = e[1];
    x.push(points.x[i], points.x[j], null);
    y.push(points.y[i], points.y[j], null);
    z.push(points.z[i], points.z[j], null);
  });
  return {
    type: "scatter3d",
    mode: "lines",
    x, y, z,
    name,
    hoverinfo: "skip",
    line: {color, width},
    opacity
  };
}

function pointTrace(points, name) {
  return {
    type: "scatter3d",
    mode: "markers",
    x: points.x,
    y: points.y,
    z: points.z,
    name,
    hoverinfo: "skip",
    marker: {
      size: 4.5,
      color: points.color,
      opacity: 0.95,
      symbol: "circle",
      line: {color: "rgba(40,40,40,0.25)", width: 0.4}
    }
  };
}

function layout(title) {
  return {
    title: {text: title, font: {size: 14}},
    margin: {l: 0, r: 0, t: 34, b: 0},
    showlegend: false,
    paper_bgcolor: "white",
    scene: {
      aspectmode: "cube",
      xaxis: {visible: false},
      yaxis: {visible: false},
      zaxis: {visible: false},
      camera: {eye: {x: 1.45, y: 1.45, z: 1.05}}
    }
  };
}

function edgeDifference(phateEdges, iknnEdges) {
  const pSet = edgeSet(phateEdges);
  const iSet = edgeSet(iknnEdges);
  const shared = [], phateOnly = [], iknnOnly = [];
  phateEdges.forEach(e => {
    if (iSet.has(edgeKey(e))) shared.push(e);
    else phateOnly.push(e);
  });
  iknnEdges.forEach(e => {
    if (!pSet.has(edgeKey(e))) iknnOnly.push(e);
  });
  return {shared, phateOnly, iknnOnly};
}

function renderPanel(divId, title, points, traces) {
  Plotly.react(divId, traces, layout(title), {displayModeBar: true, responsive: true});
}

function render() {
  const c = currentCase();
  const mode = iknnModeSel.value;
  const edgeMode = edgeModeSel.value;
  const phateGraph = DATA.graphs[c.phate_graph];
  const iknnGraph = DATA.graphs[mode === "same" ? c.iknn_same_graph : c.iknn_matched_graph];
  const pts = pointArrays(c.point_id);
  const gripPhate = embeddingArrays(c.grip_phate_embedding);
  const gripIknn = embeddingArrays(mode === "same" ? c.grip_iknn_same_embedding : c.grip_iknn_matched_embedding);
  const diffusionPhate = embeddingArrays(c.diffusion_phate_embedding);
  const diffusionIknn = embeddingArrays(mode === "same" ? c.diffusion_iknn_same_embedding : c.diffusion_iknn_matched_embedding);
  const gripOverlay = overlayArrays(mode === "same" ? c.grip_same_overlay : c.grip_matched_overlay);
  const diffusionOverlay = overlayArrays(mode === "same" ? c.diffusion_same_overlay : c.diffusion_matched_overlay);

  const iknnTitle = mode === "same" ?
    `ikNN graph, k=${iknnGraph.k}` :
    `ikNN graph, density-matched k=${iknnGraph.k}`;
  const phateTitle = `PHATE graph, k=${phateGraph.k}, decay=${phateGraph.decay}`;

  const iknnTraces = [pointTrace(pts, "points")];
  const phateTraces = [pointTrace(pts, "points")];
  if (edgeMode === "own" || edgeMode === "both" || edgeMode === "difference") {
    iknnTraces.unshift(edgeTrace(pts, iknnGraph.edges, "rgba(40,91,168,0.45)", "ikNN edges", 1.6, 0.65));
    phateTraces.unshift(edgeTrace(pts, phateGraph.edges, "rgba(202,45,38,0.42)", "PHATE edges", 1.6, 0.65));
  }

  const overlayTraces = [pointTrace(pts, "points")];
  if (edgeMode === "both") {
    overlayTraces.unshift(edgeTrace(pts, iknnGraph.edges, "rgba(35,84,164,0.42)", "ikNN", 1.5, 0.55));
    overlayTraces.unshift(edgeTrace(pts, phateGraph.edges, "rgba(215,48,39,0.42)", "PHATE", 1.5, 0.55));
  } else if (edgeMode === "difference") {
    const diff = edgeDifference(phateGraph.edges, iknnGraph.edges);
    overlayTraces.unshift(edgeTrace(pts, diff.iknnOnly, "rgba(35,84,164,0.62)", "ikNN only", 1.7, 0.7));
    overlayTraces.unshift(edgeTrace(pts, diff.phateOnly, "rgba(215,48,39,0.62)", "PHATE only", 1.7, 0.7));
    overlayTraces.unshift(edgeTrace(pts, diff.shared, "rgba(35,35,35,0.45)", "shared", 1.4, 0.55));
  }

  renderPanel("plotIknn", iknnTitle, pts, iknnTraces);
  renderPanel("plotPhate", phateTitle, pts, phateTraces);
  renderPanel("plotOverlay", edgeMode === "difference" ? "Overlay: edge differences" : "Overlay", pts, overlayTraces);

  renderPanel("plotGripIknn", `Weighted GRIP embedding: ${iknnTitle}`, gripIknn, [pointTrace(gripIknn, "ikNN weighted GRIP")]);
  renderPanel("plotGripPhate", `Weighted GRIP embedding: ${phateTitle}`, gripPhate, [pointTrace(gripPhate, "PHATE weighted GRIP")]);
  renderPanel(
    "plotGripOverlay",
    `Weighted GRIP overlay, residual ${fmt(gripOverlay.residual, 3)}`,
    gripOverlay.target,
    [pointTrace(gripOverlay.target, "PHATE graph"), pointTrace(gripOverlay.source, "ikNN aligned")]
  );

  renderPanel("plotDiffusionIknn", `PHATE embedding: ${iknnTitle}`, diffusionIknn, [pointTrace(diffusionIknn, "ikNN-support PHATE")]);
  renderPanel("plotDiffusionPhate", `PHATE embedding: ${phateTitle}`, diffusionPhate, [pointTrace(diffusionPhate, "native PHATE")]);
  renderPanel(
    "plotDiffusionOverlay",
    `PHATE embedding overlay, residual ${fmt(diffusionOverlay.residual, 3)}`,
    diffusionOverlay.target,
    [pointTrace(diffusionOverlay.target, "PHATE graph"), pointTrace(diffusionOverlay.source, "ikNN aligned")]
  );

  const cmp = mode === "same" ? c.same : c.matched;
  const im = iknnGraph.metrics;
  const pm = phateGraph.metrics;
  const ds = DATA.datasets[c.point_id];
  metricsBox.innerHTML = `
    <div><strong>${ds.dataset_label}</strong>, n=${c.n}, nominal k=${c.k}, PHATE decay=${c.decay}, ikNN mode=${mode === "same" ? "same k" : "density matched"}</div>
    <div class="metric-grid">
      <span>|E_PHATE| <b>${pm.edge_count}</b></span>
      <span>|E_ikNN| <b>${im.edge_count}</b></span>
      <span>Jaccard <b>${fmt(cmp.jaccard, 3)}</b></span>
      <span>Shared <b>${cmp.shared}</b></span>
      <span>PHATE-only <b>${fmt(cmp.phate_only_fraction, 3)}</b></span>
      <span>ikNN-only <b>${fmt(cmp.iknn_only_fraction, 3)}</b></span>
      <span>PHATE comps <b>${pm.n_components}</b></span>
      <span>ikNN comps <b>${im.n_components}</b></span>
      <span>PHATE LCC <b>${fmt(pm.largest_component_fraction, 3)}</b></span>
      <span>ikNN LCC <b>${fmt(im.largest_component_fraction, 3)}</b></span>
      <span>Matched ikNN k <b>${c.iknn_matched_k}</b></span>
      <span>PHATE median degree <b>${fmt(pm.median_degree, 1)}</b></span>
      <span>GRIP overlay residual <b>${fmt(gripOverlay.residual, 3)}</b></span>
      <span>PHATE embedding residual <b>${fmt(diffusionOverlay.residual, 3)}</b></span>
    </div>`;
}

function initialize() {
  const labels = {};
  DATA.controls.dataset.forEach(d => labels[d] = DATA.datasets[d + "__50"].dataset_label);
  optionList(datasetSel, DATA.controls.dataset, labels);
  optionList(nSel, DATA.controls.n);
  optionList(kSel, DATA.controls.k);
  optionList(decaySel, DATA.controls.decay);
  optionList(iknnModeSel, ["same", "matched"], {"same": "same k", "matched": "density matched"});
  optionList(edgeModeSel, ["difference", "own", "both", "none"], {
    "difference": "difference overlay",
    "own": "own panels",
    "both": "both overlay",
    "none": "none"
  });
  datasetSel.value = "noisy_circle";
  nSel.value = "100";
  kSel.value = "5";
  decaySel.value = "40";
  iknnModeSel.value = "same";
  edgeModeSel.value = "difference";
  [datasetSel, nSel, kSel, decaySel, iknnModeSel, edgeModeSel].forEach(el => el.addEventListener("change", render));
  render();
}

window.addEventListener("load", initialize);
'

js <- sub("__PAYLOAD__", payload.json, js, fixed = TRUE)

html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\">",
  "<title>PHATE and ikNN Graph Construction Comparison</title>",
  "<script>window.MathJax={tex:{inlineMath:[[\"$\",\"$\"],[\"\\\\(\",\"\\\\)\"]],displayMath:[[\"\\\\[\",\"\\\\]\"],[\"$$\",\"$$\"]]}};</script>",
  "<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>",
  "<script src=\"lib/plotly-latest.min.js\"></script>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:34px;line-height:1.48;color:#222;background:#fff}",
  "h1,h2{color:#16324f}h3{margin-top:26px;color:#243b53}.note{background:#f7f9fb;border-left:4px solid #3B6EA8;padding:12px 14px;margin:16px 0}",
  ".controls{display:grid;grid-template-columns:repeat(6,minmax(120px,1fr));gap:12px;align-items:end;margin:18px 0;padding:14px;background:#f6f8fa;border:1px solid #d9e2ec}",
  "label{font-size:13px;font-weight:600;color:#334e68}select{width:100%;margin-top:4px;padding:7px 8px;border:1px solid #bcccdc;background:white;border-radius:4px}",
  ".plots{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-top:14px}.plot{height:430px;border:1px solid #d9e2ec;background:white}",
  ".metrics{border:1px solid #d9e2ec;padding:12px;background:#fff;margin:12px 0 14px}.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(130px,1fr));gap:6px 14px;margin-top:8px;font-size:13px}",
  "table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0}th,td{border:1px solid #ddd;padding:5px 7px;text-align:left}th{background:#f2f5f8}",
  "code{background:#f4f4f4;padding:1px 3px}.math-block{overflow-x:auto;margin:12px 0}",
  "@media(max-width:1100px){.controls{grid-template-columns:repeat(2,1fr)}.plots{grid-template-columns:1fr}.plot{height:460px}.metric-grid{grid-template-columns:repeat(2,1fr)}}",
  "</style></head><body>",
  "<h1>PHATE and ikNN Graph Construction Comparison</h1>",
  sprintf("<p><strong>Generated:</strong> %s</p>", html.escape(payload$generated)),
  "<div class=\"note\">This report compares graph supports on the same observed 3D point clouds. PHATE uses an adaptive affinity kernel and row-normalized diffusion operator; ikNN uses shared-neighborhood graph construction. mKNN is explained conceptually but is intentionally omitted from the interactive example panels.</div>",
  "<h2>PHATE Graph Construction</h2>",
  "<p>For observations \\(x_1,\\ldots,x_n\\), PHATE first computes pairwise distances \\(d_{ij}=d(x_i,x_j)\\), usually Euclidean distances after optional PCA. For each point \\(i\\), the local bandwidth \\(\\sigma_i\\) is based on the distance to its \\(k\\)-th non-self nearest neighbor:</p>",
  "<div class=\"math-block\">\\[\\sigma_i = b\\,d_{i,(k)},\\]</div>",
  "<p>where \\(b\\) is the bandwidth scale. PHATE then forms a directed adaptive alpha-decay affinity</p>",
  "<div class=\"math-block\">\\[A_{ij}=\\exp\\left[-\\left(\\frac{d_{ij}}{\\sigma_i}\\right)^a\\right],\\qquad A_{ii}=1,\\]</div>",
  "<p>where \\(a\\) is the PHATE decay parameter. Affinities below the threshold \\(\\tau=10^{-4}\\) are set to zero:</p>",
  "<div class=\"math-block\">\\[A_{ij}=0\\quad\\text{if}\\quad A_{ij}\\lt \\tau.\\]</div>",
  "<p>The default PHATE/graphtools symmetrization is additive:</p>",
  "<div class=\"math-block\">\\[K_{ij}=\\frac{A_{ij}+A_{ji}}{2}.\\]</div>",
  "<p>The resulting kernel is finally row-normalized into a Markov diffusion operator:</p>",
  "<div class=\"math-block\">\\[P_{ij}=\\frac{K_{ij}}{\\sum_{\\ell=1}^{n}K_{i\\ell}}.\\]</div>",
  "<p>The critical point is that PHATE's \\(k\\) defines local bandwidths; it does not simply keep exactly \\(k\\) graph edges per point. The decay and threshold control the effective radius of retained affinities. Smaller decay values broaden the graph, while larger decay values make the support closer to a hard local neighborhood.</p>",
  "<h2>How This Differs From ikNN and mKNN</h2>",
  "<p><strong>ikNN</strong> connects points using overlap between their local neighbor sets. In gflow's intersection kNN graph, two points are connected when their \\(k\\)-nearest-neighbor sets share at least one common point. The graph is therefore based on shared local context rather than an adaptive exponential affinity.</p>",
  "<p><strong>mKNN</strong> is stricter: points \\(i\\) and \\(j\\) are connected only when each is in the other's \\(k\\)-nearest-neighbor list. This report does not visualize mKNN, but it is conceptually useful as a contrast: mKNN asks for reciprocity, ikNN asks for shared neighbors, and PHATE asks for non-negligible adaptive affinity after local bandwidth scaling.</p>",
  "<h2>Embedding Rows Added To The Graph View</h2>",
  "<p>The first interactive row still draws graph edges on the original sampled coordinates. This is the right view for graph-support comparison, but it is not a graph layout or manifold embedding. To ask how each graph behaves as a geometry, the report adds two embedding rows.</p>",
  "<h3>Weighted GRIP Embeddings</h3>",
  "<p>For weighted GRIP, both graph supports use the same edge-length convention: the Euclidean length of an edge in the original observed coordinates,</p>",
  "<div class=\"math-block\">\\[\\ell_{ij}=\\lVert x_i-x_j\\rVert_2,\\qquad (i,j)\\in E.\\]</div>",
  "<p>This intentionally isolates the effect of graph support. PHATE affinities are not used as GRIP lengths in this row, because an affinity is a similarity rather than a distance. The weighted GRIP row therefore asks: if the selected edge support is treated as a weighted geometric graph with ambient edge lengths, what 3D layout does it imply?</p>",
  "<h3>PHATE Embeddings From Each Graph Support</h3>",
  "<p>The native PHATE embedding row uses the PHATE kernel and Markov operator described above. For the ikNN version, the ikNN edge set supplies the support, but PHATE-style adaptive affinities are used on that support:</p>",
  "<div class=\"math-block\">\\[K^{I,\\mathrm{dir}}_{ij}=\\mathbf{1}\\{(i,j)\\in E_{\\mathrm{ikNN}}\\}\\exp\\left[-\\left(\\frac{d_{ij}}{\\sigma_i}\\right)^a\\right],\\qquad K^{I,\\mathrm{dir}}_{ii}=1.\\]</div>",
  "<p>The directed support-restricted affinity is then symmetrized and row-normalized:</p>",
  "<div class=\"math-block\">\\[K^I_{ij}=\\frac{K^{I,\\mathrm{dir}}_{ij}+K^{I,\\mathrm{dir}}_{ji}}{2},\\qquad P^I_{ij}=\\frac{K^I_{ij}}{\\sum_{\\ell=1}^{n}K^I_{i\\ell}}.\\]</div>",
  "<p>Both the native PHATE graph and the ikNN-supported PHATE graph are embedded through PHATE diffusion-potential distances. In this interactive grid, the final 3D coordinates use classical MDS for speed and determinism; the purpose is comparative geometry over the full parameter grid, not final publication-quality MDS optimization.</p>",
  "<h2>Interactive Comparison</h2>",
  "<p>The controls select a dataset, sample size, nominal \\(k\\), PHATE decay, and whether the ikNN panel uses the same \\(k\\) or a density-matched ikNN graph. Density matching chooses</p>",
  "<div class=\"math-block\">\\[k^*_{\\mathrm{ikNN}}=\\arg\\min_{k'}\\left|\\,|E_{\\mathrm{ikNN}}(k')|-|E_{\\mathrm{PHATE}}(k,a)|\\,\\right|.\\]</div>",
  "<p>The overlay difference mode colors shared edges gray, PHATE-only edges red, and ikNN-only edges blue. In the two embedding overlay panels, the PHATE graph embedding is black and the ikNN-derived embedding is red after Procrustes centering, rotation/reflection, and scalar alignment.</p>",
  "<div class=\"controls\">",
  "<div><label for=\"dataset\">Dataset</label><select id=\"dataset\"></select></div>",
  "<div><label for=\"n\">n</label><select id=\"n\"></select></div>",
  "<div><label for=\"k\">k</label><select id=\"k\"></select></div>",
  "<div><label for=\"decay\">PHATE decay</label><select id=\"decay\"></select></div>",
  "<div><label for=\"iknnMode\">ikNN mode</label><select id=\"iknnMode\"></select></div>",
  "<div><label for=\"edgeMode\">Edge display</label><select id=\"edgeMode\"></select></div>",
  "</div>",
  "<div id=\"metricsBox\" class=\"metrics\"></div>",
  "<h3>Row 1: Original-Coordinate Graph View</h3>",
  "<p>Vertices are plotted at the original sampled coordinates; only the edge supports differ.</p>",
  "<div class=\"plots\"><div id=\"plotIknn\" class=\"plot\"></div><div id=\"plotPhate\" class=\"plot\"></div><div id=\"plotOverlay\" class=\"plot\"></div></div>",
  "<h3>Row 2: Weighted GRIP Embeddings</h3>",
  "<p>Each graph support is embedded by weighted GRIP using Euclidean edge lengths on the selected support.</p>",
  "<div class=\"plots\"><div id=\"plotGripIknn\" class=\"plot\"></div><div id=\"plotGripPhate\" class=\"plot\"></div><div id=\"plotGripOverlay\" class=\"plot\"></div></div>",
  "<h3>Row 3: PHATE Diffusion Embeddings</h3>",
  "<p>The PHATE graph uses its native affinity support. The ikNN version uses ikNN support with PHATE-style adaptive affinities and row-normalized diffusion.</p>",
  "<div class=\"plots\"><div id=\"plotDiffusionIknn\" class=\"plot\"></div><div id=\"plotDiffusionPhate\" class=\"plot\"></div><div id=\"plotDiffusionOverlay\" class=\"plot\"></div></div>",
  "<h2>Graph Support Metrics</h2>",
  "<p>For each same- or matched-density comparison, the report computes the edge-support Jaccard index</p>",
  "<div class=\"math-block\">\\[J(E_P,E_I)=\\frac{|E_P\\cap E_I|}{|E_P\\cup E_I|},\\]</div>",
  "<p>as well as PHATE-only and ikNN-only fractions. Low Jaccard values mean the two methods are selecting substantially different local connections, even when the edge counts are similar.</p>",
  "<h3>Lowest Same-k Jaccard Cases</h3>",
  summary.table.html,
  "<script>", js, "</script>",
  "</body></html>"
)

html.path <- file.path(out.dir, "phate_iknn_graph_comparison_report.html")
writeLines(html, html.path)

message(sprintf("Wrote %s", normalizePath(html.path, mustWork = FALSE)))
message(sprintf("Wrote %s", normalizePath(file.path(out.dir, "phate_iknn_graph_payload.json"), mustWork = FALSE)))
message(sprintf("Wrote %s", normalizePath(file.path(out.dir, "phate_iknn_graph_comparison_summary.csv"), mustWork = FALSE)))

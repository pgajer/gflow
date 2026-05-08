#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg.value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) {
    return(default)
  }
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

as.logical.arg <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "yes", "y")
}

repo.dir <- normalizePath(getwd(), mustWork = TRUE)
mode <- arg.value("mode", "smoke")
if (!mode %in% c("smoke", "full")) {
  stop("--mode must be smoke or full.")
}
output.dir <- arg.value(
  "output.dir",
  file.path(repo.dir, "dev", "data-geodesic-reconstruction", "quadform-3d-tier2-benchmark", "runs", mode)
)
n.ref <- as.integer(arg.value("n.ref", "5000"))
k.min <- as.integer(arg.value("k.min", "3"))
k.max <- as.integer(arg.value("k.max", "25"))
max.widgets <- as.integer(arg.value("max.widgets", if (identical(mode, "smoke")) "4" else "8"))
batch.index <- as.integer(arg.value("batch.index", "1"))
batch.count <- as.integer(arg.value("batch.count", "1"))
dataset.ids.arg <- arg.value("dataset.ids", NULL)
report.only <- as.logical.arg(arg.value("report.only"), FALSE)

if (!is.finite(k.min) || !is.finite(k.max) || k.min < 1L || k.max < k.min) {
  stop("Invalid k range: expected 1 <= k.min <= k.max.")
}
if (!is.finite(batch.index) || !is.finite(batch.count) || batch.count < 1L || batch.index < 1L || batch.index > batch.count) {
  stop("Invalid batch arguments: expected 1 <= batch.index <= batch.count.")
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The pkgload package is required.")
}
pkgload::load_all(repo.dir, quiet = TRUE)

required.packages <- "ggplot2"
if (is.finite(max.widgets) && max.widgets > 0L) {
  required.packages <- c(required.packages, "plotly", "htmlwidgets", "grip")
}
missing.packages <- required.packages[!vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing.packages)) {
  stop("Missing required packages: ", paste(missing.packages, collapse = ", "))
}

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
results.dir <- file.path(output.dir, "results")
report.dir <- file.path(output.dir, "report")
figure.dir <- file.path(report.dir, "figures")
widget.dir <- file.path(report.dir, "widgets")
for (path in c(results.dir, report.dir, figure.dir, widget.dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

metrics.file <- file.path(results.dir, "metrics.csv")
datasets.file <- file.path(results.dir, "datasets.csv")
settings.file <- file.path(results.dir, "settings.csv")
diagnostics.file <- file.path(results.dir, "graph_diagnostics.csv")
best.file <- file.path(results.dir, "best_settings.csv")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

sample.ball <- function(n, dim, radius = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Z <- matrix(stats::rnorm(n * dim), nrow = n, ncol = dim)
  norms <- sqrt(rowSums(Z^2))
  Z <- Z / pmax(norms, .Machine$double.eps)
  radii <- radius * stats::runif(n)^(1 / dim)
  Z * radii
}

sample.cube <- function(n, dim, radius = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  matrix(stats::runif(n * dim, min = -radius, max = radius), nrow = n, ncol = dim)
}

sample.domain <- function(n, dim, shape, radius = 1, seed = NULL) {
  if (identical(shape, "ball")) {
    return(sample.ball(n, dim, radius = radius, seed = seed))
  }
  if (identical(shape, "cube")) {
    return(sample.cube(n, dim, radius = radius, seed = seed))
  }
  stop("Unsupported domain shape: ", shape)
}

surface.label <- function(index.k, coefficients) {
  paste0(
    "index k = ", index.k,
    ", coefficients = (", paste(coefficients, collapse = ", "), ")"
  )
}

make.datasets <- function(mode) {
  if (identical(mode, "smoke")) {
    datasets <- data.frame(
      dataset_id = c("ball_k3_c111_n80", "ball_k3_c111_n120", "ball_k1_c124_n80", "ball_k1_c124_n120"),
      surface = c("positive_curvature", "positive_curvature", "mixed_anisotropic", "mixed_anisotropic"),
      n = c(80L, 120L, 80L, 120L),
      index.k = c(3L, 3L, 1L, 1L),
      coeff_label = c("c111", "c111", "c124", "c124"),
      coefficients = I(list(c(1, 1, 1), c(1, 1, 1), c(1, 2, 4), c(1, 2, 4))),
      seed = c(3101L, 3102L, 4101L, 4102L),
      domain.shape = "ball",
      domain.radius = 1,
      stringsAsFactors = FALSE
    )
    return(datasets)
  }

  surface.spec <- data.frame(
    surface = c(
      "positive_c111", "positive_c124", "positive_c144",
      "mixed_c111", "mixed_c124", "mixed_c144"
    ),
    index.k = c(3L, 3L, 3L, 1L, 1L, 1L),
    coeff_label = c("c111", "c124", "c144", "c111", "c124", "c144"),
    stringsAsFactors = FALSE
  )
  coeff.map <- list(c111 = c(1, 1, 1), c124 = c(1, 2, 4), c144 = c(1, 4, 4))
  sizes <- c(80L, 120L, 200L)
  rows <- list()
  seed.base <- 7000L
  idx <- 0L
  for (si in seq_len(nrow(surface.spec))) {
    for (n in sizes) {
      idx <- idx + 1L
      spec <- surface.spec[si, ]
      rows[[idx]] <- data.frame(
        dataset_id = paste("ball", paste0("k", spec$index.k), spec$coeff_label, paste0("n", n), sep = "_"),
        surface = spec$surface,
        n = n,
        index.k = spec$index.k,
        coeff_label = spec$coeff_label,
        coefficients = I(list(coeff.map[[spec$coeff_label]])),
        seed = seed.base + idx,
        domain.shape = "ball",
        domain.radius = 1,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

datasets <- make.datasets(mode)
if (!is.null(dataset.ids.arg) && nzchar(dataset.ids.arg)) {
  requested.ids <- strsplit(dataset.ids.arg, ",", fixed = TRUE)[[1]]
  requested.ids <- trimws(requested.ids[nzchar(requested.ids)])
  missing.ids <- setdiff(requested.ids, datasets$dataset_id)
  if (length(missing.ids)) {
    stop("Unknown dataset ids: ", paste(missing.ids, collapse = ", "))
  }
  datasets <- datasets[datasets$dataset_id %in% requested.ids, , drop = FALSE]
}
if (batch.count > 1L) {
  batch.slot <- ((seq_len(nrow(datasets)) - 1L) %% batch.count) + 1L
  datasets <- datasets[batch.slot == batch.index, , drop = FALSE]
}
if (!nrow(datasets)) {
  stop("No datasets selected for this run.")
}
datasets$description <- vapply(
  seq_len(nrow(datasets)),
  function(i) surface.label(datasets$index.k[i], datasets$coefficients[[i]]),
  character(1)
)

make.settings <- function() {
  rows <- list()
  add <- function(row) rows[[length(rows) + 1L]] <<- row
  for (k in k.min:k.max) {
    add(data.frame(family = "sknn", k = k, k_scale = NA_integer_, radius.rule = NA_character_, radius.factor = NA_real_, delta = NA_real_))
    add(data.frame(family = "iknn", k = k, k_scale = NA_integer_, radius.rule = NA_character_, radius.factor = NA_real_, delta = NA_real_))
  }
  for (k.scale in k.min:k.max) {
    for (radius.rule in c("max", "geomean")) {
      for (radius.factor in c(1, 1.25, 1.5)) {
        add(data.frame(family = "adaptive_radius", k = NA_integer_, k_scale = k.scale, radius.rule = radius.rule, radius.factor = radius.factor, delta = NA_real_))
      }
    }
    for (delta in c(1, 1.25, 1.5)) {
      add(data.frame(family = "cknn", k = NA_integer_, k_scale = k.scale, radius.rule = "geomean", radius.factor = NA_real_, delta = delta))
    }
  }
  settings <- do.call(rbind, rows)
  settings$setting_id <- sprintf("s%03d", seq_len(nrow(settings)))
  settings[, c("setting_id", "family", "k", "k_scale", "radius.rule", "radius.factor", "delta")]
}

settings <- make.settings()

build.graph <- function(X, setting) {
  family <- as.character(setting$family[[1]])
  if (identical(family, "sknn")) {
    return(create.sknn.graph(
      X,
      k = as.integer(setting$k[[1]]),
      prune.method = "none",
      connect.components = TRUE,
      connect.method = "component.mst"
    ))
  }
  if (identical(family, "iknn")) {
    return(create.single.iknn.graph(
      X,
      k = as.integer(setting$k[[1]]),
      prune.method = "none",
      threshold.percentile = 0,
      connect.components = TRUE,
      connect.method = "component.mst",
      with.lifecycle.branches = TRUE,
      pca.dim = NULL,
      verbose = FALSE
    ))
  }
  if (identical(family, "adaptive_radius")) {
    return(create.adaptive.radius.graph(
      X,
      k.scale = as.integer(setting$k_scale[[1]]),
      radius.rule = as.character(setting$radius.rule[[1]]),
      radius.factor = as.numeric(setting$radius.factor[[1]]),
      prune.method = "none",
      connect.components = TRUE,
      connect.method = "component.mst"
    ))
  }
  if (identical(family, "cknn")) {
    return(create.cknn.graph(
      X,
      k.scale = as.integer(setting$k_scale[[1]]),
      delta = as.numeric(setting$delta[[1]]),
      prune.method = "none",
      connect.components = TRUE,
      connect.method = "component.mst"
    ))
  }
  stop("Unsupported graph family: ", family)
}

adj.edges <- function(adj.list) {
  rows <- do.call(
    rbind,
    lapply(seq_along(adj.list), function(i) {
      nbrs <- adj.list[[i]]
      if (!length(nbrs)) {
        return(NULL)
      }
      j <- nbrs[nbrs > i]
      if (!length(j)) {
        return(NULL)
      }
      cbind(i = i, j = j)
    })
  )
  if (is.null(rows)) {
    matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j")))
  } else {
    rows
  }
}

graph.edge.count <- function(graph) {
  adj <- graph$adj_list %||% graph$final_adj_list %||% graph$raw_adj_list
  if (is.null(adj)) {
    return(NA_integer_)
  }
  sum(lengths(adj)) / 2
}

graph.components <- function(adj.list) {
  n <- length(adj.list)
  seen <- rep(FALSE, n)
  count <- 0L
  for (start in seq_len(n)) {
    if (seen[start]) {
      next
    }
    count <- count + 1L
    queue <- start
    seen[start] <- TRUE
    while (length(queue)) {
      v <- queue[[1]]
      queue <- queue[-1]
      nbrs <- adj.list[[v]]
      unseen <- nbrs[!seen[nbrs]]
      if (length(unseen)) {
        seen[unseen] <- TRUE
        queue <- c(queue, unseen)
      }
    }
  }
  count
}

make.plot <- function(data, filename, width = 8, height = 5) {
  ggplot2::ggsave(
    filename = file.path(figure.dir, filename),
    plot = data,
    width = width,
    height = height,
    dpi = 150
  )
  file.path("figures", filename)
}

html.escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

format.num <- function(x, digits = 4) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "fg"))
}

table.html <- function(df, digits = 4) {
  if (!nrow(df)) {
    return("<p>No rows available.</p>")
  }
  df2 <- df
  for (j in seq_along(df2)) {
    if (is.numeric(df2[[j]])) {
      df2[[j]] <- format.num(df2[[j]], digits = digits)
    } else {
      df2[[j]] <- html.escape(as.character(df2[[j]]))
    }
  }
  header <- paste0("<th>", html.escape(names(df2)), "</th>", collapse = "")
  rows <- apply(df2, 1, function(row) {
    paste0("<tr>", paste0("<td>", row, "</td>", collapse = ""), "</tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

plotly.graph.view <- function(coords, edges, title, point.color = "#1f77b4", edge.color = "rgba(60,60,60,0.28)") {
  coords <- as.matrix(coords)
  p <- plotly::plot_ly(type = "scatter3d", mode = "markers") |>
    plotly::add_markers(
      x = coords[, 1], y = coords[, 2], z = coords[, 3],
      marker = list(size = 3, color = point.color),
      hoverinfo = "none",
      showlegend = FALSE
    )
  if (nrow(edges)) {
    edge.limit <- min(nrow(edges), 1200L)
    edges <- edges[seq_len(edge.limit), , drop = FALSE]
    xe <- as.vector(rbind(coords[edges[, 1], 1], coords[edges[, 2], 1], NA))
    ye <- as.vector(rbind(coords[edges[, 1], 2], coords[edges[, 2], 2], NA))
    ze <- as.vector(rbind(coords[edges[, 1], 3], coords[edges[, 2], 3], NA))
    p <- p |>
      plotly::add_trace(
        x = xe, y = ye, z = ze,
        type = "scatter3d",
        mode = "lines",
        line = list(width = 1, color = edge.color),
        hoverinfo = "none",
        showlegend = FALSE
      )
  }
  p |>
    plotly::layout(
      title = list(text = title, font = list(size = 12)),
      scene = list(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        zaxis = list(title = ""),
        aspectmode = "data"
      ),
      margin = list(l = 0, r = 0, b = 0, t = 30)
    )
}

save.widget.pair <- function(dataset.row, graph, setting, widget.name) {
  adj <- graph$adj_list
  weight <- graph$weight_list
  edges <- adj.edges(adj)
  original.coords <- dataset.row$X_param[[1]]
  layout.coords <- tryCatch(
    grip::grip.layout.weighted(
      adj_list = adj,
      weight_list = weight,
      dim = 3,
      rounds = 60,
      final_rounds = 96,
      seed = dataset.row$seed[[1]]
    ),
    error = function(e) {
      warning("GRIP layout failed for ", widget.name, ": ", conditionMessage(e))
      stats::cmdscale(graph.geodesic.distances(graph, stage = "final"), k = 3)
    }
  )
  layout.coords <- as.matrix(layout.coords)
  p1 <- plotly.graph.view(original.coords, edges, "Parameter-domain coordinates", "#555555", "rgba(80,80,80,0.20)")
  p2 <- plotly.graph.view(layout.coords, edges, "Weighted GRIP layout", "#1f77b4", "rgba(20,80,160,0.25)")
  widget <- plotly::subplot(p1, p2, nrows = 1, shareX = FALSE, shareY = FALSE, titleX = FALSE, titleY = FALSE) |>
    plotly::layout(
      title = list(
        text = paste0(
          dataset.row$dataset_id[[1]], " / ",
          setting$family, " / setting ", setting$setting_id
        ),
        font = list(size = 14)
      )
    )
  path <- file.path(widget.dir, paste0(widget.name, ".html"))
  suppressWarnings(htmlwidgets::saveWidget(widget, path, selfcontained = TRUE))
  file.path("widgets", basename(path))
}

if (!report.only) {
  dataset.rows <- list()
  metric.rows <- list()
  diagnostic.rows <- list()
  graph.cache <- new.env(parent = emptyenv())

  utils::write.csv(datasets, datasets.file, row.names = FALSE)
  utils::write.csv(settings, settings.file, row.names = FALSE)

  for (di in seq_len(nrow(datasets))) {
    drow <- datasets[di, ]
    message("[dataset] ", drow$dataset_id, " (n = ", drow$n, ", ", drow$description, ")")
    X.param <- sample.domain(drow$n, 3, drow$domain.shape, radius = drow$domain.radius, seed = drow$seed)
    X.surface <- quadform.embed(X.param, index.k = drow$index.k, coefficients = drow$coefficients[[1]])
    ref <- quadform.delaunay.geodesic.distances(
      X.param,
      index.k = drow$index.k,
      coefficients = drow$coefficients[[1]],
      domain.radius = drow$domain.radius,
      domain.shape = drow$domain.shape,
      n.ref = n.ref,
      seed = drow$seed + 9000L,
      candidate.multiplier = 6,
      boundary.fraction = 0.25,
      edge.length.factor = 4,
      delaunay.backend = "cpp"
    )
    dataset.rows[[length(dataset.rows) + 1L]] <- data.frame(
      dataset_id = drow$dataset_id,
      surface = drow$surface,
      n = drow$n,
      index.k = drow$index.k,
      coefficients = paste(drow$coefficients[[1]], collapse = ";"),
      seed = drow$seed,
      n.ref = n.ref,
      n.ref.actual = ref$n_reference_vertices %||% NA_integer_,
      n.delaunay.edges = ref$n_delaunay_edges %||% NA_integer_,
      n.retained.edges = ref$n_edges %||% NA_integer_,
      edge.length.factor = ref$filter_factor_used %||% 4,
      stringsAsFactors = FALSE
    )

    for (si in seq_len(nrow(settings))) {
      setting <- settings[si, ]
      key <- paste(drow$dataset_id, setting$setting_id, sep = "__")
      t0 <- proc.time()[["elapsed"]]
      result <- tryCatch({
        graph <- build.graph(X.surface, setting)
        Dg <- graph.geodesic.distances(graph, stage = "final")
        summary <- summarize.isometry.deviation(Dg, ref$distances, scale = TRUE)
        elapsed <- proc.time()[["elapsed"]] - t0
        n.edges <- graph.edge.count(graph)
        metric.rows[[length(metric.rows) + 1L]] <- cbind(
          data.frame(
            dataset_id = drow$dataset_id,
            surface = drow$surface,
            n = drow$n,
            setting_id = setting$setting_id,
            family = setting$family,
            k = setting$k,
            k_scale = setting$k_scale,
            radius.rule = setting$radius.rule,
            radius.factor = setting$radius.factor,
            delta = setting$delta,
            n_edges = n.edges,
            elapsed_sec = elapsed,
            stringsAsFactors = FALSE
          ),
          summary
        )
        diagnostic.rows[[length(diagnostic.rows) + 1L]] <- data.frame(
          dataset_id = drow$dataset_id,
          setting_id = setting$setting_id,
          family = setting$family,
          n_edges = n.edges,
          n_components_final = graph.components(graph$adj_list),
          n_mst_edges_added = graph$n_mst_edges_added %||% graph$mst_edges_added %||% NA_integer_,
          stringsAsFactors = FALSE
        )
        assign(key, list(dataset = drow, X_param = X.param, graph = graph, setting = setting), envir = graph.cache)
        TRUE
      }, error = function(e) {
        message("    setting failed: ", drow$dataset_id, " / ", setting$setting_id, " / ", setting$family, ": ", conditionMessage(e))
        elapsed <- proc.time()[["elapsed"]] - t0
        metric.rows <<- c(metric.rows, list(data.frame(
          dataset_id = drow$dataset_id,
          surface = drow$surface,
          n = drow$n,
          setting_id = setting$setting_id,
          family = setting$family,
          k = setting$k,
          k_scale = setting$k_scale,
          radius.rule = setting$radius.rule,
          radius.factor = setting$radius.factor,
          delta = setting$delta,
          n_edges = NA_real_,
          elapsed_sec = elapsed,
          scale = NA_real_,
          rel_rms_error = NA_real_,
          rel_abs_error_median = NA_real_,
          rel_abs_error_q95 = NA_real_,
          distortion_q05 = NA_real_,
          distortion_median = NA_real_,
          distortion_q95 = NA_real_,
          pearson_cor = NA_real_,
          spearman_cor = NA_real_,
          rel_geodesic_stress = NA_real_,
          signed_bias = NA_real_,
          shortcut_fraction = NA_real_,
          q50_rel_abs_residual = NA_real_,
          q90_rel_abs_residual = NA_real_,
          q95_rel_abs_residual = NA_real_,
          short_band_bias = NA_real_,
          mid_band_bias = NA_real_,
          long_band_bias = NA_real_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )))
        FALSE
      })
      if (si %% 20L == 0L || identical(si, nrow(settings))) {
        message("  settings complete: ", si, "/", nrow(settings), "; last ok = ", result)
      }
    }
  }
  metrics <- do.call(rbind, metric.rows)
  if (!"error" %in% names(metrics)) {
    metrics$error <- NA_character_
  }
  utils::write.csv(do.call(rbind, dataset.rows), datasets.file, row.names = FALSE)
  utils::write.csv(metrics, metrics.file, row.names = FALSE)
  utils::write.csv(do.call(rbind, diagnostic.rows), diagnostics.file, row.names = FALSE)

  ok <- metrics[is.na(metrics$error) & is.finite(metrics$rel_rms_error), ]
  best <- do.call(rbind, lapply(split(ok, list(ok$dataset_id, ok$family), drop = TRUE), function(df) {
    df[which.min(df$rel_rms_error), ]
  }))
  best <- best[order(best$dataset_id, best$rel_rms_error), ]
  utils::write.csv(best, best.file, row.names = FALSE)

  best.overall <- do.call(rbind, lapply(split(ok, ok$dataset_id), function(df) df[which.min(df$rel_rms_error), ]))
  if (is.finite(max.widgets) && nrow(best.overall) > max.widgets) {
    best.overall <- best.overall[seq_len(max.widgets), ]
  }
  for (i in seq_len(nrow(best.overall))) {
    row <- best.overall[i, ]
    key <- paste(row$dataset_id, row$setting_id, sep = "__")
    payload <- get(key, envir = graph.cache)
    payload$dataset$X_param <- list(payload$X_param)
    rel <- save.widget.pair(
      payload$dataset,
      payload$graph,
      payload$setting,
      paste0(row$dataset_id, "_", row$setting_id)
    )
    best.overall$widget[i] <- rel
  }
  utils::write.csv(best.overall, file.path(results.dir, "best_overall_widgets.csv"), row.names = FALSE)
}

metrics <- utils::read.csv(metrics.file, stringsAsFactors = FALSE)
datasets.out <- utils::read.csv(datasets.file, stringsAsFactors = FALSE)
best <- utils::read.csv(best.file, stringsAsFactors = FALSE)
best.widgets.file <- file.path(results.dir, "best_overall_widgets.csv")
best.widgets <- if (file.exists(best.widgets.file)) utils::read.csv(best.widgets.file, stringsAsFactors = FALSE) else data.frame()

ok <- metrics[is.na(metrics$error) & is.finite(metrics$rel_rms_error), ]

best.plot.data <- best
best.plot.data$surface_n <- paste0(best.plot.data$surface, "\nN = ", best.plot.data$n)
p.best <- ggplot2::ggplot(best.plot.data, ggplot2::aes(x = family, y = rel_rms_error, fill = family)) +
  ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
  ggplot2::facet_wrap(~ surface_n, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Best relative RMS error", title = "Best surface-target performance by graph family") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
best.plot <- make.plot(p.best, "family_best_surface_error.png", width = 9, height = 5)

knn.data <- ok[ok$family %in% c("sknn", "iknn"), ]
p.knn <- ggplot2::ggplot(knn.data, ggplot2::aes(x = k, y = rel_rms_error, color = family)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::facet_grid(surface ~ n, labeller = ggplot2::label_both) +
  ggplot2::labs(x = "k", y = "Relative RMS error", title = "sKNN and iKNN k sensitivity") +
  ggplot2::theme_bw(base_size = 11)
knn.plot <- make.plot(p.knn, "knn_k_sensitivity.png", width = 9, height = 5)

adaptive.data <- ok[ok$family == "adaptive_radius", ]
adaptive.data$radius.factor <- factor(adaptive.data$radius.factor)
p.adaptive <- ggplot2::ggplot(adaptive.data, ggplot2::aes(x = k_scale, y = rel_rms_error, color = radius.factor, linetype = radius.rule)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.1) +
  ggplot2::facet_grid(surface ~ n, labeller = ggplot2::label_both) +
  ggplot2::labs(x = "k.scale", y = "Relative RMS error", color = "radius factor", linetype = "radius rule", title = "Adaptive-radius parameter sensitivity") +
  ggplot2::theme_bw(base_size = 11)
adaptive.plot <- make.plot(p.adaptive, "adaptive_radius_parameter_sensitivity.png", width = 10, height = 5.5)

cknn.data <- ok[ok$family == "cknn", ]
cknn.data$delta <- factor(cknn.data$delta)
p.cknn <- ggplot2::ggplot(cknn.data, ggplot2::aes(x = k_scale, y = rel_rms_error, color = delta)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::facet_grid(surface ~ n, labeller = ggplot2::label_both) +
  ggplot2::labs(x = "k.scale", y = "Relative RMS error", color = "delta", title = "Continuous-kNN parameter sensitivity") +
  ggplot2::theme_bw(base_size = 11)
cknn.plot <- make.plot(p.cknn, "cknn_parameter_sensitivity.png", width = 9, height = 5)

p.edges <- ggplot2::ggplot(ok, ggplot2::aes(x = n_edges, y = rel_rms_error, color = family)) +
  ggplot2::geom_point(alpha = 0.7, size = 1.6) +
  ggplot2::facet_grid(surface ~ n, labeller = ggplot2::label_both) +
  ggplot2::labs(x = "Number of graph edges", y = "Relative RMS error", title = "Graph size versus surface-target error") +
  ggplot2::theme_bw(base_size = 11)
edges.plot <- make.plot(p.edges, "graph_size_vs_error.png", width = 9, height = 5)

best.summary <- best[, c("dataset_id", "surface", "n", "family", "setting_id", "k", "k_scale", "radius.rule", "radius.factor", "delta", "n_edges", "rel_rms_error", "rel_geodesic_stress", "signed_bias", "shortcut_fraction")]
best.summary <- best.summary[order(best.summary$dataset_id, best.summary$rel_rms_error), ]

family.summary <- aggregate(
  cbind(rel_rms_error, rel_geodesic_stress, n_edges) ~ surface + n + family,
  data = ok,
  FUN = function(x) median(x, na.rm = TRUE)
)
family.summary <- family.summary[order(family.summary$surface, family.summary$n, family.summary$rel_rms_error), ]

widget.blocks <- ""
if (nrow(best.widgets)) {
  widget.blocks <- paste(vapply(seq_len(nrow(best.widgets)), function(i) {
    row <- best.widgets[i, ]
    paste0(
      "<section class=\"widget-card\"><h3>", html.escape(row$dataset_id), " / ",
      html.escape(row$family), " / ", html.escape(row$setting_id), "</h3>",
      "<iframe src=\"", html.escape(row$widget), "\" loading=\"lazy\"></iframe></section>"
    )
  }, character(1)), collapse = "\n")
}

report.path <- file.path(report.dir, "quadform_3d_tier2_benchmark_report.html")
html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>3D Quadratic Hypersurface Tier 2 Benchmark</title>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;line-height:1.45;margin:28px;color:#222;}",
  "h1,h2{line-height:1.15;} img{max-width:100%;border:1px solid #ddd;} table{border-collapse:collapse;font-size:13px;margin:16px 0;width:100%;}",
  "th,td{border:1px solid #ddd;padding:5px 7px;text-align:left;} th{background:#f4f4f4;} .note{color:#555;max-width:980px;}",
  ".widget-card{margin:24px 0;} iframe{width:100%;height:640px;border:1px solid #ccc;border-radius:4px;}",
  "code{background:#f7f7f7;padding:1px 4px;border-radius:3px;}",
  "</style></head><body>",
  "<h1>3D Quadratic Hypersurface Tier 2 Benchmark</h1>",
  "<p class=\"note\">Build timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "</p>",
  "<p>Run mode: <code>", html.escape(mode), "</code>. Parameter range: <code>k, k.scale = ", k.min, ":", k.max, "</code>. Reference target size: <code>", n.ref, "</code>. ",
  "Batch: <code>", batch.index, "/", batch.count, "</code>.</p>",
  "<p>This benchmark checks data-to-graph reconstruction on three-dimensional quadratic hypersurface domains. ",
  "It uses the C++ Delaunay reference oracle in <code>quadform.delaunay.geodesic.distances()</code> with <code>edge.length.factor = 4</code>. ",
  "Only the surface geodesic target is evaluated; sample-oracle targets and pruning are intentionally excluded.</p>",
  "<p>The graph families are <code>sKNN</code>, <code>iKNN</code>, adaptive-radius graphs, and continuous-kNN/geometric-mean adaptive-radius graphs. ",
  "All graphs are repaired with the common component-MST method before graph-geodesic distances are compared to the surface reference.</p>",
  "<h2>Datasets</h2>",
  table.html(datasets.out[, c("dataset_id", "surface", "n", "index.k", "coefficients", "n.ref", "n.ref.actual", "n.delaunay.edges", "n.retained.edges")]),
  "<h2>Main Results</h2>",
  "<p>The first figure shows the best setting found within each graph family for each dataset. These are the settings a later full run should study more densely.</p>",
  "<img src=\"", best.plot, "\" alt=\"Best family performance\">",
  "<p>The next plots expand the parameter sensitivity for the individual graph families. They are deliberately shown before summary tables because the trend shape matters more than a single winning row in this smoke run.</p>",
  "<img src=\"", knn.plot, "\" alt=\"KNN sensitivity\">",
  "<img src=\"", adaptive.plot, "\" alt=\"Adaptive radius sensitivity\">",
  "<img src=\"", cknn.plot, "\" alt=\"CKNN sensitivity\">",
  "<p>The graph-size plot is retained as a diagnostic: low error with a very dense graph is less interesting than low error with a compact support.</p>",
  "<img src=\"", edges.plot, "\" alt=\"Graph size versus error\">",
  "<h2>Best Settings</h2>",
  table.html(best.summary),
  "<h2>Family-Level Median Performance</h2>",
  table.html(family.summary),
  "<h2>Interactive 3D Diagnostics</h2>",
  "<p>The left panel uses the three-dimensional parameter-domain coordinates of the sampled hypersurface graph. ",
  "The right panel is a weighted GRIP layout of the selected data-derived graph. These views are meant as visual diagnostics, not as the error metric itself.</p>",
  widget.blocks,
  "<h2>Artifacts</h2>",
  "<p>CSV outputs are written under <code>", html.escape(results.dir), "</code>. The report assets live under <code>", html.escape(report.dir), "</code>.</p>",
  "</body></html>"
)
writeLines(html, report.path)
message("Report written to: ", report.path)

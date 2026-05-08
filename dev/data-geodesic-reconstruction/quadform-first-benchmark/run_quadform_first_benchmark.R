#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(RGL_USE_NULL = "TRUE")

parse_args <- function(args) {
  out <- list(
    mode = "smoke",
    grid.size = NA_integer_,
    output.dir = "dev/data-geodesic-reconstruction/quadform-first-benchmark",
    resume = TRUE,
    layouts = FALSE,
    max.layout.rows = 6L,
    setting.timeout.sec = NA_integer_,
    workers = NA_integer_,
    rerun.errors = FALSE,
    report.only = FALSE,
    layout.graph.key = NA_character_,
    save.graph.assets = TRUE,
    save.layout.assets = TRUE,
    asset.stages = "all"
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = FALSE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    if (key %in% names(out)) {
      out[[key]] <- value
    }
  }
  out$mode <- match.arg(out$mode, c("smoke", "full"))
  out$resume <- as.logical(out$resume)
  out$layouts <- as.logical(out$layouts)
  out$grid.size <- suppressWarnings(as.integer(out$grid.size))
  out$max.layout.rows <- suppressWarnings(as.integer(out$max.layout.rows))
  out$setting.timeout.sec <- suppressWarnings(as.integer(out$setting.timeout.sec))
  out$workers <- suppressWarnings(as.integer(out$workers))
  out$rerun.errors <- as.logical(out$rerun.errors)
  out$report.only <- as.logical(out$report.only)
  out$save.graph.assets <- as.logical(out$save.graph.assets)
  out$save.layout.assets <- as.logical(out$save.layout.assets)
  out$asset.stages <- match.arg(
    as.character(out$asset.stages),
    c("benchmark", "available", "all")
  )
  out$layout.graph.key <- as.character(out$layout.graph.key)
  if (!length(out$layout.graph.key) || is.na(out$layout.graph.key) ||
      !nzchar(out$layout.graph.key)) {
    out$layout.graph.key <- NA_character_
  }
  if (is.na(out$grid.size)) {
    out$grid.size <- if (identical(out$mode, "smoke")) 101L else 501L
  }
  if (is.na(out$max.layout.rows) || out$max.layout.rows < 0L) {
    out$max.layout.rows <- 0L
  }
  if (is.na(out$setting.timeout.sec)) {
    out$setting.timeout.sec <- if (identical(out$mode, "smoke")) 120L else 180L
  }
  if (is.na(out$workers) || out$workers < 1L) {
    cores <- parallel::detectCores(logical = TRUE)
    out$workers <- max(1L, min(4L, cores - 1L))
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for this dev benchmark.", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE)

required.pkgs <- c("jsonlite", "ggplot2", "htmltools")
missing.pkgs <- required.pkgs[!vapply(required.pkgs, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing.pkgs)) {
  stop("Missing required packages: ", paste(missing.pkgs, collapse = ", "), call. = FALSE)
}

base.dir <- normalizePath(args$output.dir, mustWork = FALSE)
run.dir <- file.path(base.dir, "runs", args$mode)
cache.dir <- file.path(run.dir, "cache")
asset.dir <- file.path(run.dir, "assets")
dataset.asset.dir <- file.path(asset.dir, "datasets")
graph.asset.dir <- file.path(asset.dir, "graphs")
layout.asset.dir <- file.path(asset.dir, "layouts")
report.dir <- file.path(run.dir, "report")
fig.dir <- file.path(report.dir, "figures")
widget.dir <- file.path(report.dir, "widgets")
widget.cache.dir <- file.path(report.dir, "widget_cache")
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dataset.asset.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graph.asset.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(layout.asset.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(widget.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(widget.cache.dir, recursive = TRUE, showWarnings = FALSE)

log.path <- file.path(run.dir, "progress.log")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log.path, append = TRUE)
}

surface_specs <- function() {
  data.frame(
    surface = c("paraboloid", "saddle"),
    index.k = c(2L, 1L),
    formula = c("q(u,v) = u^2 + v^2", "q(u,v) = u^2 - v^2"),
    stringsAsFactors = FALSE
  )
}

mode_config <- function(mode, grid.size) {
  if (identical(mode, "smoke")) {
    list(
      surfaces = surface_specs(),
      n.values = 50L,
      seeds = 1L,
      k.values = 3:4,
      radius.rank = 3:4,
      k.scale = 3:4,
      radius.rule = c("min", "max"),
      radius.factor = c(1, 1.25),
      cknn.delta = c(1, 1.25),
      prune.method = c("none", "local.geodesic", "global.geodesic.ratio"),
      local.prune.tau = 1.05,
      global.max.path.edge.ratio.deviation.thld = 0.1,
      global.path.edge.ratio.percentile = 0.5,
      grid.size = grid.size,
      sample.connection.k = 8L,
      oracle.tube.k = 8L,
      domain.radius = 1
    )
  } else {
    list(
      surfaces = surface_specs(),
      n.values = c(50L, 100L, 200L),
      seeds = 1:5,
      k.values = 3:10,
      radius.rank = 3:10,
      k.scale = 3:10,
      radius.rule = c("min", "max"),
      radius.factor = c(1, 1.25, 1.5),
      cknn.delta = c(1, 1.25, 1.5),
      prune.method = c("none", "local.geodesic", "global.geodesic.ratio"),
      local.prune.tau = 1.05,
      global.max.path.edge.ratio.deviation.thld = 0.1,
      global.path.edge.ratio.percentile = 0.5,
      grid.size = grid.size,
      sample.connection.k = 8L,
      oracle.tube.k = 8L,
      domain.radius = 1
    )
  }
}

config <- mode_config(args$mode, args$grid.size)
config$mode <- args$mode
config$created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
config$output_dir <- run.dir
config$layouts <- args$layouts
config$save_graph_assets <- args$save.graph.assets
config$save_layout_assets <- args$save.layout.assets
config$asset_stages <- args$asset.stages
config$setting_timeout_sec <- args$setting.timeout.sec
config$workers <- args$workers
run.config.path <- file.path(run.dir, "run_config.json")
if (isTRUE(args$report.only) && file.exists(run.config.path)) {
  existing.config <- tryCatch(
    jsonlite::fromJSON(run.config.path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.list(existing.config)) {
    config <- existing.config
    config$layouts <- args$layouts
    config$save_graph_assets <- args$save.graph.assets
    config$save_layout_assets <- args$save.layout.assets
    config$asset_stages <- args$asset.stages
  }
} else {
  jsonlite::write_json(config, run.config.path, pretty = TRUE,
                       auto_unbox = TRUE)
}

escape_html <- function(x) {
  htmltools::htmlEscape(as.character(x))
}

table_html <- function(df, digits = 4L, max.rows = Inf) {
  if (is.null(df) || !nrow(df)) {
    return("<p>No rows.</p>")
  }
  if (is.finite(max.rows) && nrow(df) > max.rows) {
    df <- df[seq_len(max.rows), , drop = FALSE]
  }
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      ifelse(is.na(col), NA_character_, format(round(col, digits), nsmall = 0,
                                               scientific = FALSE, trim = TRUE))
    } else {
      as.character(col)
    }
  })
  head <- paste0("<tr>", paste(sprintf("<th>%s</th>", escape_html(names(df))),
                               collapse = ""), "</tr>")
  rows <- apply(df, 1L, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", escape_html(row)), collapse = ""), "</tr>")
  })
  paste0("<table>", head, paste(rows, collapse = "\n"), "</table>")
}

upper_tri <- function(D) {
  D[upper.tri(D)]
}

safe_summary <- function(D.estimated, D.true) {
  if (any(!is.finite(upper_tri(D.estimated)))) {
    stop("Graph geodesic distance matrix contains non-finite off-diagonal values.",
         call. = FALSE)
  }
  if (any(!is.finite(upper_tri(D.true)))) {
    stop("Reference distance matrix contains non-finite off-diagonal values.",
         call. = FALSE)
  }
  summarize.isometry.deviation(D.estimated, D.true)
}

with_elapsed_timeout <- function(expr, timeout.sec) {
  if (!is.finite(timeout.sec) || timeout.sec <= 0L) {
    return(force(expr))
  }
  old <- setTimeLimit(cpu = Inf, elapsed = timeout.sec, transient = TRUE)
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)
  force(expr)
}

pairwise_nonzero_distances <- function(X) {
  as.numeric(stats::dist(X))
}

build_settings <- function(cfg, n) {
  prune <- cfg$prune.method
  settings <- list()
  add <- function(df) {
    settings[[length(settings) + 1L]] <<- df
  }
  add(expand.grid(
    graph_family = "sknn", k = cfg$k.values, prune_method = prune,
    stringsAsFactors = FALSE
  ))
  add(expand.grid(
    graph_family = "mknn", k = cfg$k.values, prune_method = prune,
    stringsAsFactors = FALSE
  ))
  add(expand.grid(
    graph_family = "iknn", k = cfg$k.values, prune_method = prune,
    stringsAsFactors = FALSE
  ))
  add(expand.grid(
    graph_family = "fixed_radius", radius_rank = cfg$radius.rank,
    prune_method = prune, stringsAsFactors = FALSE
  ))
  add(expand.grid(
    graph_family = "adaptive_radius", k_scale = cfg$k.scale,
    radius_rule = cfg$radius.rule, radius_factor = cfg$radius.factor,
    prune_method = prune, stringsAsFactors = FALSE
  ))
  cknn.settings <- expand.grid(
    graph_family = "cknn", k_scale = cfg$k.scale, delta = cfg$cknn.delta,
    prune_method = prune, stringsAsFactors = FALSE
  )
  cknn.settings$radius_rule <- "geomean"
  cknn.settings$radius_factor <- cknn.settings$delta
  add(cknn.settings)
  out <- rbind.fill.simple(settings)
  out$stage <- ifelse(out$prune_method == "none", "raw.repaired", "repaired.pruned")
  out$setting_id <- sprintf("g%04d", seq_len(nrow(out)))
  out
}

rbind.fill.simple <- function(x) {
  cols <- unique(unlist(lapply(x, names)))
  rows <- lapply(x, function(df) {
    missing <- setdiff(cols, names(df))
    for (col in missing) {
      df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

make_dataset_id <- function(surface, n, seed) {
  sprintf("%s_n%d_seed%03d", surface, n, seed)
}

dataset_cache_path <- function(dataset.id) {
  file.path(cache.dir, paste0("dataset_", dataset.id, ".rds"))
}

reference_cache_path <- function(dataset.id) {
  file.path(cache.dir, paste0("reference_", dataset.id, "_grid", config$grid.size, ".rds"))
}

result_cache_path <- function(dataset.id, setting.id) {
  file.path(cache.dir, "graph_results", dataset.id, paste0(setting.id, ".rds"))
}

safe_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
}

dataset_asset_path <- function(dataset.id) {
  file.path(dataset.asset.dir, paste0(safe_stem(dataset.id), ".rds"))
}

graph_stage_asset_path <- function(dataset.id, setting.id, stage) {
  file.path(
    graph.asset.dir,
    safe_stem(dataset.id),
    safe_stem(setting.id),
    paste0(safe_stem(stage), ".rds")
  )
}

layout_stage_asset_path <- function(dataset.id, setting.id, stage) {
  file.path(
    layout.asset.dir,
    safe_stem(dataset.id),
    safe_stem(setting.id),
    paste0(safe_stem(stage), "_weighted_grip_3d.rds")
  )
}

all_graph_stages <- function() {
  c("raw", "raw.repaired", "pruned", "pruned.repaired", "repaired.pruned", "final")
}

asset_stages_for_setting <- function(setting) {
  if (identical(args$asset.stages, "benchmark")) {
    unique(as.character(setting$stage))
  } else {
    all_graph_stages()
  }
}

stage_payload <- function(g, stage) {
  switch(stage,
         raw = list(adj = g$raw_adj_list, weight = g$raw_weight_list),
         raw.repaired = list(adj = g$raw_repaired_adj_list,
                             weight = g$raw_repaired_weight_list),
         pruned = list(adj = g$pruned_adj_list, weight = g$pruned_weight_list),
         pruned.repaired = list(adj = g$pruned_repaired_adj_list,
                                weight = g$pruned_repaired_weight_list),
         repaired.pruned = list(adj = g$repaired_pruned_adj_list,
                                weight = g$repaired_pruned_weight_list),
         final = list(adj = g$adj_list, weight = g$weight_list),
         stop("Unknown stage: ", stage, call. = FALSE))
}

compute_grip_layout <- function(adj, weight, seed) {
  if (!requireNamespace("grip", quietly = TRUE)) return(NULL)
  if (!exists("grip.layout.weighted", envir = asNamespace("grip"), inherits = FALSE)) {
    return(NULL)
  }
  fun <- get("grip.layout.weighted", envir = asNamespace("grip"), inherits = FALSE)
  tryCatch(
    fun(adj_list = adj, weight_list = weight, dim = 3, rounds = 8,
        final_rounds = 12, seed = seed),
    error = function(e) NULL
  )
}

save_dataset_asset <- function(dataset.id, ds, meta) {
  path <- dataset_asset_path(dataset.id)
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      list(
        dataset_id = dataset.id,
        surface = as.character(meta$surface),
        index_k = as.integer(meta$index.k),
        n = as.integer(meta$n),
        seed = as.integer(meta$seed),
        X = ds$X_embed,
        X_embed = ds$X_embed,
        X_param = ds$X_param,
        source = "quadform.sample.dataset",
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      ),
      path
    )
  }
  path
}

save_graph_and_layout_assets <- function(g, dataset.id, setting, meta) {
  stages <- asset_stages_for_setting(setting)
  graph.rows <- list()
  layout.rows <- list()

  for (stage in stages) {
    payload <- tryCatch(stage_payload(g, stage), error = function(e) NULL)
    if (is.null(payload) || is.null(payload$adj) || is.null(payload$weight)) {
      next
    }

    graph.path <- graph_stage_asset_path(dataset.id, setting$setting_id, stage)
    if (isTRUE(args$save.graph.assets)) {
      dir.create(dirname(graph.path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(
        list(
          dataset_id = dataset.id,
          setting_id = as.character(setting$setting_id),
          stage = stage,
          graph_family = as.character(setting$graph_family),
          parameters = as.list(setting),
          adj_list = payload$adj,
          weight_list = payload$weight,
          n_vertices = length(payload$adj),
          n_edges = count_edges(payload$adj),
          n_components = count_components(payload$adj),
          diagnostics = extract_graph_diagnostics(g, meta),
          created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
        ),
        graph.path
      )
    }

    graph.rows[[length(graph.rows) + 1L]] <- data.frame(
      dataset_id = dataset.id,
      setting_id = as.character(setting$setting_id),
      stage = stage,
      graph_family = as.character(setting$graph_family),
      graph_asset_file = normalizePath(graph.path, mustWork = FALSE),
      n_vertices = length(payload$adj),
      n_edges = count_edges(payload$adj),
      n_components = count_components(payload$adj),
      stringsAsFactors = FALSE
    )

    if (isTRUE(args$save.layout.assets)) {
      layout.path <- layout_stage_asset_path(dataset.id, setting$setting_id, stage)
      layout <- if (file.exists(layout.path)) {
        existing.layout <- tryCatch(readRDS(layout.path), error = function(e) NULL)
        if (is.list(existing.layout) && !is.null(existing.layout$coords)) {
          existing.layout$coords
        } else {
          existing.layout
        }
      } else {
        compute_grip_layout(
          payload$adj,
          payload$weight,
          seed = 100000L + as.integer(meta$seed[[1L]]) * 1000L +
            suppressWarnings(as.integer(sub("^g", "", setting$setting_id))) +
            match(stage, all_graph_stages(), nomatch = 0L)
        )
      }
      if (!is.null(layout)) {
        dir.create(dirname(layout.path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(
          list(
            dataset_id = dataset.id,
            setting_id = as.character(setting$setting_id),
            stage = stage,
            method = "weighted_grip",
            coords = layout,
            params = list(dim = 3, rounds = 8, final_rounds = 12),
            graph_asset_file = normalizePath(graph.path, mustWork = FALSE),
            created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
          ),
          layout.path
        )
        layout.rows[[length(layout.rows) + 1L]] <- data.frame(
          dataset_id = dataset.id,
          setting_id = as.character(setting$setting_id),
          stage = stage,
          method = "weighted_grip",
          layout_asset_file = normalizePath(layout.path, mustWork = FALSE),
          n_vertices = nrow(as.matrix(layout)),
          stringsAsFactors = FALSE
        )
      } else {
        layout.rows[[length(layout.rows) + 1L]] <- data.frame(
          dataset_id = dataset.id,
          setting_id = as.character(setting$setting_id),
          stage = stage,
          method = "weighted_grip",
          layout_asset_file = NA_character_,
          n_vertices = NA_integer_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    graph_assets = if (length(graph.rows)) rbind.fill.simple(graph.rows) else data.frame(),
    layout_assets = if (length(layout.rows)) rbind.fill.simple(layout.rows) else data.frame()
  )
}

result_has_requested_assets <- function(res) {
  if (!isTRUE(args$save.graph.assets) && !isTRUE(args$save.layout.assets)) {
    return(TRUE)
  }
  if (!is.list(res) || !is.data.frame(res$metrics) ||
      !all(res$metrics$status == "ok")) {
    return(TRUE)
  }
  if (isTRUE(args$save.graph.assets) &&
      (!is.data.frame(res$graph_assets) || nrow(res$graph_assets) < 1L)) {
    return(FALSE)
  }
  if (isTRUE(args$save.layout.assets) &&
      (!is.data.frame(res$layout_assets) || nrow(res$layout_assets) < 1L ||
       !any(!is.na(res$layout_assets$layout_asset_file)))) {
    return(FALSE)
  }
  TRUE
}

count_edges <- function(adj.list) {
  if (is.null(adj.list)) {
    return(NA_real_)
  }
  sum(lengths(adj.list)) / 2
}

count_components <- function(adj.list) {
  if (is.null(adj.list)) {
    return(NA_integer_)
  }
  n <- length(adj.list)
  seen <- rep(FALSE, n)
  n.comp <- 0L
  for (i in seq_len(n)) {
    if (seen[[i]]) {
      next
    }
    n.comp <- n.comp + 1L
    queue <- i
    seen[[i]] <- TRUE
    cursor <- 1L
    while (cursor <= length(queue)) {
      v <- queue[[cursor]]
      cursor <- cursor + 1L
      nb <- adj.list[[v]]
      nb <- nb[nb >= 1L & nb <= n]
      new <- nb[!seen[nb]]
      if (length(new)) {
        seen[new] <- TRUE
        queue <- c(queue, new)
      }
    }
  }
  n.comp
}

extract_graph_diagnostics <- function(g, meta) {
  stages <- c("raw", "raw.repaired", "pruned", "pruned.repaired", "repaired.pruned", "final")
  edge.names <- paste0("n_edges_", gsub("\\.", "_", stages))
  comp.names <- paste0("n_components_", gsub("\\.", "_", stages))
  edges <- comps <- vector("list", length(stages))
  for (i in seq_along(stages)) {
    payload <- tryCatch(stage_payload(g, stages[[i]]), error = function(e) NULL)
    edges[[i]] <- if (is.null(payload)) NA_real_ else count_edges(payload$adj)
    comps[[i]] <- if (is.null(payload)) NA_integer_ else count_components(payload$adj)
  }
  names(edges) <- edge.names
  names(comps) <- comp.names
  data.frame(
    meta,
    n_vertices = length(g$adj_list %||% g$raw_adj_list),
    as.data.frame(edges, check.names = FALSE),
    as.data.frame(comps, check.names = FALSE),
    n_components_before = scalar_or_na(g$n_components_before %||% g$n_components_before_mst),
    n_components_after = scalar_or_na(g$n_components_after %||% g$n_components_after_mst),
    n_mst_edges_added = scalar_or_na(g$n_mst_edges_added),
    bridge_method = scalar_or_na(g$bridge_method),
    bridge_k = scalar_or_na(g$bridge_k),
    bridge_k_max = scalar_or_na(g$bridge_k_max),
    bridge_k_used = scalar_or_na(g$bridge_k_used),
    bridge_exact_fallback_used = scalar_or_na(g$bridge_exact_fallback_used),
    n_edges_before_pruning = scalar_or_na(g$n_edges_before_pruning),
    n_edges_after_pruning = scalar_or_na(g$n_edges_after_pruning),
    n_pruned_edges = scalar_or_na(g$n_pruned_edges),
    prune_tau = scalar_or_na(g$prune_tau),
    prune_local_k = scalar_or_na(g$prune_local_k),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x)) {
    return(NA)
  }
  x[[1L]]
}

build_graph <- function(X, setting, radius.lookup) {
  family <- setting$graph_family
  prune.method <- setting$prune_method
  if (identical(family, "sknn")) {
    create.sknn.graph(
      X, k = setting$k, prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$k,
      connect.components = TRUE,
      connect.method = "component.mst"
    )
  } else if (identical(family, "mknn")) {
    create.mknn.graph(
      X, k = setting$k, prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$k,
      connect.components = TRUE,
      connect.method = "component.mst"
    )
  } else if (identical(family, "iknn")) {
    create.single.iknn.graph(
      X, k = setting$k, prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$k,
      threshold.percentile = 0,
      connect.components = TRUE, connect.method = "component.mst",
      with.lifecycle.branches = TRUE, pca.dim = NULL, verbose = FALSE
    )
  } else if (identical(family, "fixed_radius")) {
    radius <- radius.lookup[[as.character(setting$radius_rank)]]
    create.radius.graph(
      X, radius = radius, prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$radius_rank,
      connect.components = TRUE,
      connect.method = "component.mst"
    )
  } else if (identical(family, "adaptive_radius")) {
    create.adaptive.radius.graph(
      X, k.scale = setting$k_scale, radius.rule = setting$radius_rule,
      radius.factor = setting$radius_factor, prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$k_scale,
      connect.components = TRUE, connect.method = "component.mst"
    )
  } else if (identical(family, "cknn")) {
    create.cknn.graph(
      X, k.scale = setting$k_scale, delta = setting$delta,
      prune.method = prune.method,
      max.path.edge.ratio.deviation.thld = config$global.max.path.edge.ratio.deviation.thld,
      path.edge.ratio.percentile = config$global.path.edge.ratio.percentile,
      prune.tau = config$local.prune.tau, prune.local.k = setting$k_scale,
      connect.components = TRUE, connect.method = "component.mst"
    )
  } else {
    stop("Unknown graph family: ", family, call. = FALSE)
  }
}

save_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
}

dataset.manifest.path <- file.path(run.dir, "dataset_manifest.csv")
if (isTRUE(args$report.only) && file.exists(dataset.manifest.path)) {
  dataset_manifest <- utils::read.csv(dataset.manifest.path, stringsAsFactors = FALSE)
} else {
  dataset_manifest <- do.call(rbind, lapply(seq_len(nrow(config$surfaces)), function(i) {
    surface <- config$surfaces$surface[[i]]
    index.k <- config$surfaces$index.k[[i]]
    expand.grid(surface = surface, index.k = index.k, n = config$n.values,
                seed = config$seeds, stringsAsFactors = FALSE)
  }))
  dataset_manifest$dataset_id <- mapply(make_dataset_id, dataset_manifest$surface,
                                        dataset_manifest$n, dataset_manifest$seed)
  save_csv(dataset_manifest, dataset.manifest.path)
}

run_one_dataset <- function(row) {
  dataset.id <- row$dataset_id
  log_msg("Dataset ", dataset.id)

  ds.path <- dataset_cache_path(dataset.id)
  if (args$resume && file.exists(ds.path)) {
    ds <- readRDS(ds.path)
  } else {
    ds <- quadform.sample.dataset(
      n = row$n,
      index.k = row$index.k,
      domain.radius = config$domain.radius,
      sample.method = "uniform.parameter.disk",
      seed = row$seed,
      grid.size = 5,
      sample.connection.k = min(4L, row$n - 1L)
    )
    saveRDS(ds, ds.path)
  }
  ds.asset.path <- save_dataset_asset(dataset.id, ds, row)

  ref.path <- reference_cache_path(dataset.id)
  if (args$resume && file.exists(ref.path)) {
    ref <- readRDS(ref.path)
  } else {
    ref <- quadform.grid.geodesic.distances(
      X = ds$X_param,
      index.k = row$index.k,
      domain.radius = config$domain.radius,
      grid.size = config$grid.size,
      sample.connection.k = config$sample.connection.k,
      oracle = "sample.path",
      oracle.tube.k = config$oracle.tube.k,
      oracle.tube.radius = NULL,
      return.oracle.paths = FALSE
    )
    saveRDS(ref, ref.path)
  }

  settings <- build_settings(config, row$n)
  d <- pairwise_nonzero_distances(ds$X_embed)
  radius.values <- stats::quantile(
    d,
    probs = config$radius.rank / (row$n - 1),
    names = FALSE,
    type = 7
  )
  radius.lookup <- stats::setNames(as.numeric(radius.values), as.character(config$radius.rank))

  out.dir <- dirname(result_cache_path(dataset.id, settings$setting_id[[1L]]))
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

  run_one_setting <- function(idx) {
    setting <- settings[idx, , drop = FALSE]
    result.path <- result_cache_path(dataset.id, setting$setting_id)
    if (args$resume && file.exists(result.path)) {
      res <- readRDS(result.path)
      if ((!isTRUE(args$rerun.errors) || all(res$metrics$status == "ok")) &&
          result_has_requested_assets(res)) {
        return(res)
      }
    }

    meta <- data.frame(
      dataset_id = dataset.id,
      surface = row$surface,
      index_k = row$index.k,
      n = row$n,
      seed = row$seed,
      setting,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    res <- tryCatch({
      with_elapsed_timeout({
        g <- build_graph(ds$X_embed, setting, radius.lookup)
        D.graph <- graph.geodesic.distances(g, stage = setting$stage)
        surface.summary <- safe_summary(D.graph, ref$distances)
        oracle.summary <- safe_summary(D.graph, ref$oracle_distances)
        surface.summary$target <- "surface"
        oracle.summary$target <- "sample_oracle"
        metric <- rbind(surface.summary, oracle.summary)
        metric <- cbind(meta, status = "ok", error = NA_character_, metric)
        diag <- extract_graph_diagnostics(g, meta)
        diag$status <- "ok"
        diag$error <- NA_character_
        assets <- save_graph_and_layout_assets(g, dataset.id, setting, meta)
        list(
          metrics = metric,
          diagnostics = diag,
          dataset_asset = data.frame(
            dataset_id = dataset.id,
            dataset_asset_file = normalizePath(ds.asset.path, mustWork = FALSE),
            stringsAsFactors = FALSE
          ),
          graph_assets = assets$graph_assets,
          layout_assets = assets$layout_assets
        )
      }, args$setting.timeout.sec)
    }, error = function(e) {
      metric <- cbind(
        meta,
        status = "error",
        error = conditionMessage(e),
        target = c("surface", "sample_oracle"),
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
        long_band_bias = NA_real_
      )
      diag <- cbind(meta, status = "error", error = conditionMessage(e))
      list(
        metrics = metric,
        diagnostics = diag,
        dataset_asset = data.frame(
          dataset_id = dataset.id,
          dataset_asset_file = normalizePath(ds.asset.path, mustWork = FALSE),
          stringsAsFactors = FALSE
        ),
        graph_assets = data.frame(),
        layout_assets = data.frame()
      )
    })
    saveRDS(res, result.path)
    res
  }

  log_msg("  graph settings: ", nrow(settings), " (workers=", args$workers,
          ", timeout=", args$setting.timeout.sec, " sec)")
  if (args$workers > 1L) {
    setting.results <- parallel::mclapply(
      seq_len(nrow(settings)),
      run_one_setting,
      mc.cores = args$workers,
      mc.preschedule = FALSE
    )
  } else {
    setting.results <- lapply(seq_len(nrow(settings)), run_one_setting)
  }
  ok.count <- sum(vapply(setting.results, function(x) {
    all(x$metrics$status == "ok")
  }, logical(1L)))
  log_msg("  graph settings complete: ", ok.count, " ok, ",
          length(setting.results) - ok.count, " with errors/timeouts")
  list(metrics = rbind.fill.simple(lapply(setting.results, `[[`, "metrics")),
       diagnostics = rbind.fill.simple(lapply(setting.results, `[[`, "diagnostics")),
       dataset_assets = unique(rbind.fill.simple(lapply(setting.results, `[[`, "dataset_asset"))),
       graph_assets = rbind.fill.simple(lapply(setting.results, `[[`, "graph_assets")),
       layout_assets = rbind.fill.simple(lapply(setting.results, `[[`, "layout_assets")))
}

if (isTRUE(args$report.only)) {
  if (!file.exists(file.path(run.dir, "results.rds"))) {
    stop("Cannot use --report.only=true because results.rds does not exist.",
         call. = FALSE)
  }
  log_msg("Report-only mode: reading existing benchmark artifacts.")
} else {
  all.results <- vector("list", nrow(dataset_manifest))
  for (i in seq_len(nrow(dataset_manifest))) {
    all.results[[i]] <- run_one_dataset(dataset_manifest[i, , drop = FALSE])
    metrics.tmp <- rbind.fill.simple(lapply(all.results[seq_len(i)], `[[`, "metrics"))
    diagnostics.tmp <- rbind.fill.simple(lapply(all.results[seq_len(i)], `[[`, "diagnostics"))
    dataset.assets.tmp <- unique(rbind.fill.simple(
      lapply(all.results[seq_len(i)], `[[`, "dataset_assets")
    ))
    graph.assets.tmp <- rbind.fill.simple(
      lapply(all.results[seq_len(i)], `[[`, "graph_assets")
    )
    layout.assets.tmp <- rbind.fill.simple(
      lapply(all.results[seq_len(i)], `[[`, "layout_assets")
    )
    saveRDS(
      list(
        metrics = metrics.tmp,
        diagnostics = diagnostics.tmp,
        dataset_assets = dataset.assets.tmp,
        graph_assets = graph.assets.tmp,
        layout_assets = layout.assets.tmp,
        config = config
      ),
      file.path(run.dir, "results.rds")
    )
    save_csv(metrics.tmp, file.path(run.dir, "metrics.csv"))
    save_csv(diagnostics.tmp, file.path(run.dir, "graph_diagnostics.csv"))
    save_csv(dataset.assets.tmp, file.path(run.dir, "dataset_assets.csv"))
    save_csv(graph.assets.tmp, file.path(run.dir, "graph_assets.csv"))
    save_csv(layout.assets.tmp, file.path(run.dir, "layout_assets.csv"))
  }
}

results <- readRDS(file.path(run.dir, "results.rds"))
metrics <- results$metrics
diagnostics <- results$diagnostics
dataset.assets <- results$dataset_assets %||% data.frame()
graph.assets <- results$graph_assets %||% data.frame()
layout.assets <- results$layout_assets %||% data.frame()

geodesic.metric.cols <- c(
  "rel_geodesic_stress",
  "signed_bias",
  "shortcut_fraction",
  "q50_rel_abs_residual",
  "q90_rel_abs_residual",
  "q95_rel_abs_residual",
  "short_band_bias",
  "mid_band_bias",
  "long_band_bias"
)
for (col in geodesic.metric.cols) {
  if (!col %in% names(metrics)) {
    metrics[[col]] <- NA_real_
  }
}
if ("rel_geodesic_stress" %in% names(metrics) &&
    all(is.na(metrics$rel_geodesic_stress)) &&
    "rel_rms_error" %in% names(metrics)) {
  metrics$rel_geodesic_stress <- metrics$rel_rms_error
}
save_csv(metrics, file.path(run.dir, "metrics.csv"))

write_asset_manifests <- function() {
  if (nrow(dataset.assets)) {
    save_csv(dataset.assets, file.path(run.dir, "dataset_assets.csv"))
  }
  if (nrow(graph.assets)) {
    save_csv(graph.assets, file.path(run.dir, "graph_assets.csv"))
  }
  if (nrow(layout.assets)) {
    save_csv(layout.assets, file.path(run.dir, "layout_assets.csv"))
  }

  setting.cols <- intersect(
    c("dataset_id", "setting_id", "surface", "index_k", "n", "seed",
      "graph_family", "k", "radius_rank", "k_scale", "radius_rule",
      "radius_factor", "delta", "prune_method", "stage"),
    names(metrics)
  )
  settings <- unique(metrics[, setting.cols, drop = FALSE])
  manifest <- list(
    version = "1",
    project = "quadform_first_benchmark",
    mode = args$mode,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    run_dir = normalizePath(run.dir, mustWork = FALSE),
    asset_dir = normalizePath(asset.dir, mustWork = FALSE),
    report_file = normalizePath(
      file.path(report.dir, "quadform_first_benchmark_report.html"),
      mustWork = FALSE
    ),
    config = config,
    dataset_manifest_file = normalizePath(file.path(run.dir, "dataset_manifest.csv"), mustWork = FALSE),
    metrics_file = normalizePath(file.path(run.dir, "metrics.csv"), mustWork = FALSE),
    graph_diagnostics_file = normalizePath(file.path(run.dir, "graph_diagnostics.csv"), mustWork = FALSE),
    dataset_assets_file = normalizePath(file.path(run.dir, "dataset_assets.csv"), mustWork = FALSE),
    graph_assets_file = normalizePath(file.path(run.dir, "graph_assets.csv"), mustWork = FALSE),
    layout_assets_file = normalizePath(file.path(run.dir, "layout_assets.csv"), mustWork = FALSE),
    dataset_assets = dataset.assets,
    graph_settings = settings,
    graph_assets = graph.assets,
    layout_assets = layout.assets
  )
  saveRDS(manifest, file.path(run.dir, "quadform_benchmark_manifest.rds"))
  jsonlite::write_json(
    manifest,
    file.path(run.dir, "quadform_benchmark_manifest.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null"
  )
  invisible(manifest)
}

asset.manifest <- write_asset_manifests()

plot_family_performance <- function(metrics, path, target) {
  ok <- metrics[metrics$status == "ok" & metrics$target == target, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  agg <- stats::aggregate(rel_rms_error ~ surface + n + graph_family + prune_method,
                          data = ok, FUN = median)
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = n, y = rel_rms_error, color = graph_family,
                 linetype = prune_method, group = interaction(graph_family, prune_method))
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ surface, scales = "free_y") +
    ggplot2::labs(
      x = "Sample size", y = "Median relative RMS error",
      color = "Graph family", linetype = "Pruning",
      title = paste("Method performance against", gsub("_", " ", target), "target")
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(path, p, width = 11, height = 6, dpi = 150)
  path
}

plot_rank_distribution <- function(metrics, path, target) {
  ok <- metrics[metrics$status == "ok" & metrics$target == target, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  p <- ggplot2::ggplot(ok, ggplot2::aes(x = graph_family, y = rel_rms_error,
                                        fill = prune_method)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::facet_wrap(~ surface, scales = "free_y") +
    ggplot2::labs(x = "Graph family", y = "Relative RMS error",
                  fill = "Pruning",
                  title = paste("Distribution of graph-distance errors:",
                                gsub("_", " ", target))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(path, p, width = 11, height = 6, dpi = 150)
  path
}

plot_pruning_effect <- function(metrics, path, target) {
  ok <- metrics[metrics$status == "ok" & metrics$target == target, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  agg <- stats::aggregate(rel_rms_error ~ surface + graph_family + prune_method,
                          data = ok, FUN = median)
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = graph_family, y = rel_rms_error,
                                         fill = prune_method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ surface, scales = "free_y") +
    ggplot2::labs(x = "Graph family", y = "Median relative RMS error",
                  fill = "Pruning",
      title = paste("Effect of geometric pruning method:",
                    gsub("_", " ", target))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(path, p, width = 11, height = 6, dpi = 150)
  path
}

plot_diagnostics <- function(diagnostics, path) {
  ok <- diagnostics[diagnostics$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  for (nm in c("n_mst_edges_added", "n_pruned_edges")) {
    if (!nm %in% names(ok)) ok[[nm]] <- NA_real_
  }
  agg <- stats::aggregate(cbind(n_mst_edges_added, n_pruned_edges) ~
                            surface + n + graph_family + prune_method,
                          data = ok,
                          FUN = function(x) median(as.numeric(x), na.rm = TRUE))
  long <- rbind(
    data.frame(agg[, c("surface", "n", "graph_family", "prune_method")],
               diagnostic = "MST bridge edges added", value = agg$n_mst_edges_added),
    data.frame(agg[, c("surface", "n", "graph_family", "prune_method")],
               diagnostic = "Edges pruned", value = agg$n_pruned_edges)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(x = graph_family, y = value,
                                          fill = prune_method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(diagnostic ~ surface + n, scales = "free_y") +
    ggplot2::labs(x = "Graph family", y = "Median count", fill = "Pruning",
                  title = "Connectivity repair and pruning diagnostics") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 14, height = 7, dpi = 150)
  path
}

edge_table_from_adj <- function(adj.list, weight.list = NULL, max.edges = 800L) {
  rows <- list()
  cursor <- 0L
  for (i in seq_along(adj.list)) {
    nb <- adj.list[[i]]
    if (!length(nb)) next
    keep <- nb > i
    nb <- nb[keep]
    if (!length(nb)) next
    for (j in nb) {
      cursor <- cursor + 1L
      rows[[cursor]] <- c(i, j)
    }
  }
  if (!cursor) return(matrix(integer(), ncol = 2L))
  edges <- do.call(rbind, rows)
  if (nrow(edges) > max.edges) {
    set.seed(17)
    edges <- edges[sort(sample(seq_len(nrow(edges)), max.edges)), , drop = FALSE]
  }
  edges
}

sphere_radius <- function(coords) {
  rng <- apply(coords, 2L, range, finite = TRUE)
  span <- max(rng[2L, ] - rng[1L, ])
  if (!is.finite(span) || span <= 0) return(0.025)
  span / 90
}

save_rgl_scene <- function(coords, file, title, edges = NULL, color = NULL,
                           edge.color = "#8a8f98", point.radius = NULL) {
  if (!requireNamespace("rgl", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(NULL)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::bg3d("white")
  if (!is.null(edges) && nrow(edges)) {
    seg <- matrix(NA_real_, nrow = 2L * nrow(edges), ncol = 3L)
    seg[seq(1L, nrow(seg), by = 2L), ] <- coords[edges[, 1L], , drop = FALSE]
    seg[seq(2L, nrow(seg), by = 2L), ] <- coords[edges[, 2L], , drop = FALSE]
    rgl::segments3d(seg, color = edge.color, alpha = 0.45)
  }
  if (is.null(color)) color <- grDevices::hcl.colors(nrow(coords), "viridis")
  if (is.null(point.radius)) point.radius <- sphere_radius(coords)
  rgl::spheres3d(coords[, 1], coords[, 2], coords[, 3],
                 radius = point.radius, color = color)
  rgl::axes3d(color = "#777777")
  rgl::title3d(main = title, color = "#222222")
  widget <- rgl::rglwidget()
  htmlwidgets::saveWidget(widget, file, selfcontained = FALSE)
  file
}

safe_stem <- function(x) {
  x <- paste(as.character(x), collapse = "_")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

stage_col <- function(prefix, stage) {
  paste0(prefix, "_", gsub("\\.", "_", stage))
}

available_graph_stages <- function(row) {
  stages <- c("raw", "raw.repaired", "pruned", "pruned.repaired",
              "repaired.pruned", "final")
  ok <- vapply(stages, function(stage) {
    edge.name <- stage_col("n_edges", stage)
    edge.name %in% names(row) && is.finite(suppressWarnings(as.numeric(row[[edge.name]])))
  }, logical(1L))
  stages[ok]
}

graph_setting_key <- function(row) {
  parts <- c(
    row$dataset_id,
    row$graph_family,
    if (!is.na(row$k)) paste0("k", row$k),
    if (!is.na(row$radius_rank)) paste0("rr", row$radius_rank),
    if (!is.na(row$k_scale)) paste0("ks", row$k_scale),
    if (!is.na(row$radius_rule)) paste0("rule", row$radius_rule),
    if (!is.na(row$radius_factor)) paste0("rf", row$radius_factor),
    if ("delta" %in% names(row) && !is.na(row$delta)) paste0("delta", row$delta),
    paste0("prune", row$prune_method)
  )
  safe_stem(parts)
}

layout_stage_key <- function(graph.key, stage) {
  paste0(graph.key, "__stage_", safe_stem(stage))
}

relative_widget_path <- function(path) {
  if (is.null(path) || !length(path) || is.na(path)) return(NA_character_)
  sub(paste0("^", report.dir, "/?"), "", normalizePath(path, mustWork = FALSE))
}

write_placeholder_widget <- function(file, title, message) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
    "margin:0;background:#fbfcfe;color:#26323f;height:100vh;display:flex;",
    "align-items:center;justify-content:center}.box{max-width:520px;padding:24px;",
    "border:1px solid #d8dde6;background:white}h2{font-size:18px;margin:0 0 10px}",
    "code{background:#f2f4f7;padding:2px 4px}</style></head><body><div class=\"box\"><h2>",
    escape_html(title), "</h2><p>", message, "</p></div></body></html>"
  )
  writeLines(html, file)
  file
}

ensure_widget_placeholders <- function() {
  list(
    original = write_placeholder_widget(
      file.path(widget.dir, "missing_original.html"),
      "Original data widget not generated",
      "Run the report with <code>--layouts=true</code> to precompute 3D widgets."
    ),
    grip = write_placeholder_widget(
      file.path(widget.dir, "missing_grip.html"),
      "Weighted GRIP layout not precomputed",
      "Use the graph key shown in the metrics box with <code>--layout.graph.key=&lt;graph_key&gt;</code> to cache this layout and regenerate the report."
    ),
    stage = write_placeholder_widget(
      file.path(widget.dir, "stage_unavailable.html"),
      "Selected graph stage unavailable",
      "The selected lifecycle stage is not present for this graph object."
    )
  )
}

metric_for_target <- function(metrics, target) {
  keep <- metrics$status == "ok" & metrics$target == target
  cols <- c("dataset_id", "setting_id", "rel_rms_error", "rel_abs_error_median",
            "rel_abs_error_q95", "distortion_q05", "distortion_median",
            "distortion_q95", "pearson_cor", "spearman_cor",
            geodesic.metric.cols)
  out <- metrics[keep, intersect(cols, names(metrics)), drop = FALSE]
  metric.names <- setdiff(names(out), c("dataset_id", "setting_id"))
  names(out)[match(metric.names, names(out))] <- paste0(target, "_", metric.names)
  out
}

build_widget_index <- function(metrics, diagnostics) {
  if (!nrow(diagnostics)) return(data.frame())
  key.cols <- c("dataset_id", "surface", "index_k", "n", "seed", "setting_id",
                "graph_family", "k", "radius_rank", "k_scale", "radius_rule",
                "radius_factor", "delta", "prune_method", "stage")
  stage.cols <- grep("^n_(edges|components)_", names(diagnostics), value = TRUE)
  diag.cols <- c("n_vertices", stage.cols, "n_mst_edges_added", "n_pruned_edges",
                 "n_edges_before_pruning", "n_edges_after_pruning",
                 "n_components_before", "n_components_after")
  cols <- intersect(c(key.cols, diag.cols), names(diagnostics))
  idx <- unique(diagnostics[, cols, drop = FALSE])
  ok.keys <- paste(metrics$dataset_id[metrics$status == "ok"],
                   metrics$setting_id[metrics$status == "ok"], sep = "|")
  idx <- idx[paste(idx$dataset_id, idx$setting_id, sep = "|") %in% ok.keys, ,
             drop = FALSE]
  if (!nrow(idx)) return(data.frame())

  surface.metric <- metric_for_target(metrics, "surface")
  sample.metric <- metric_for_target(metrics, "sample_oracle")
  idx <- merge(idx, surface.metric, by = c("dataset_id", "setting_id"),
               all.x = TRUE, sort = FALSE)
  idx <- merge(idx, sample.metric, by = c("dataset_id", "setting_id"),
               all.x = TRUE, sort = FALSE)

  idx$default_stage <- idx$stage
  idx$graph_key <- vapply(seq_len(nrow(idx)), function(i) {
    graph_setting_key(idx[i, , drop = FALSE])
  }, character(1L))
  idx$available_stages <- vapply(seq_len(nrow(idx)), function(i) {
    paste(available_graph_stages(idx[i, , drop = FALSE]), collapse = ",")
  }, character(1L))
  idx$n_edges_stage <- vapply(seq_len(nrow(idx)), function(i) {
    edge.name <- stage_col("n_edges", idx$default_stage[[i]])
    if (edge.name %in% names(idx)) as.numeric(idx[[edge.name]][[i]]) else NA_real_
  }, numeric(1L))
  idx$n_components_stage <- vapply(seq_len(nrow(idx)), function(i) {
    comp.name <- stage_col("n_components", idx$default_stage[[i]])
    if (comp.name %in% names(idx)) as.numeric(idx[[comp.name]][[i]]) else NA_real_
  }, numeric(1L))
  idx$original_widget <- NA_character_
  idx$grip_widget <- NA_character_
  idx$widget_stage <- NA_character_
  idx$grip_cached <- FALSE
  idx[order(idx$surface, idx$n, idx$seed, idx$graph_family, idx$setting_id), ,
      drop = FALSE]
}

preset_row_indices <- function(widget.index) {
  ok <- widget.index[is.finite(widget.index$surface_rel_rms_error) |
                       is.finite(widget.index$sample_oracle_rel_rms_error), ,
                     drop = FALSE]
  if (!nrow(ok)) return(integer())
  groups <- split(seq_len(nrow(ok)), paste(ok$surface, ok$n, ok$seed, sep = "|"))
  selected <- integer()
  for (idx in groups) {
    sub <- ok[idx, , drop = FALSE]
    for (target.col in c("surface_rel_rms_error", "sample_oracle_rel_rms_error")) {
      values <- suppressWarnings(as.numeric(sub[[target.col]]))
      finite <- which(is.finite(values))
      if (!length(finite)) next
      ordered <- finite[order(values[finite])]
      picked <- c(ordered[[1L]], ordered[[ceiling(length(ordered) / 2)]],
                  ordered[[length(ordered)]])
      selected <- c(selected, idx[picked])
    }
  }
  unique(selected)
}

select_layout_rows <- function(widget.index, max.rows, requested.key) {
  if (!nrow(widget.index)) return(integer())
  selected <- integer()
  if (!is.na(requested.key)) {
    exact <- which(widget.index$graph_key == requested.key |
                     vapply(seq_len(nrow(widget.index)), function(i) {
                       any(vapply(c("raw", "raw.repaired", "pruned",
                                    "pruned.repaired", "repaired.pruned", "final"),
                                  function(stage) {
                                    layout_stage_key(widget.index$graph_key[[i]],
                                                     stage) == requested.key
                                  }, logical(1L)))
                     }, logical(1L)))
    selected <- c(selected, exact)
  }
  min.n <- min(widget.index$n, na.rm = TRUE)
  selected <- c(selected, which(widget.index$n == min.n &
                                  widget.index$prune_method == "none"))
  selected <- c(selected, preset_row_indices(widget.index))
  selected <- c(selected, which(widget.index$n == min.n))
  selected <- unique(selected)
  if (!length(selected)) return(integer())
  limit <- max.rows
  if (!is.na(requested.key)) {
    limit <- max(limit, length(selected[widget.index$graph_key[selected] == requested.key]))
  }
  if (!isTRUE(args$layouts)) {
    return(integer())
  }
  if (limit <= 0L) return(integer())
  selected[seq_len(min(length(selected), limit))]
}

requested_render_stage <- function(graph.key, default.stage, requested.key) {
  if (is.na(requested.key)) return(default.stage)
  stages <- c("raw", "raw.repaired", "pruned", "pruned.repaired",
              "repaired.pruned", "final")
  for (stage in stages) {
    if (identical(layout_stage_key(graph.key, stage), requested.key)) {
      return(stage)
    }
  }
  default.stage
}

asset_row_for_stage <- function(tbl, dataset.id, setting.id, stage) {
  if (!is.data.frame(tbl) || !nrow(tbl)) {
    return(NULL)
  }
  keep <- tbl$dataset_id == dataset.id &
    tbl$setting_id == setting.id &
    tbl$stage == stage
  idx <- which(keep)
  if (length(idx) < 1L) {
    return(NULL)
  }
  tbl[idx[[1L]], , drop = FALSE]
}

load_graph_stage_asset <- function(dataset.id, setting.id, stage) {
  row <- asset_row_for_stage(graph.assets, dataset.id, setting.id, stage)
  if (is.null(row) || !("graph_asset_file" %in% names(row))) {
    return(NULL)
  }
  path <- as.character(row$graph_asset_file[[1L]])
  if (!nzchar(path) || is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  g <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(g) || is.null(g$adj_list) || is.null(g$weight_list)) {
    return(NULL)
  }
  list(adj = g$adj_list, weight = g$weight_list, file = path)
}

load_layout_stage_asset <- function(dataset.id, setting.id, stage) {
  row <- asset_row_for_stage(layout.assets, dataset.id, setting.id, stage)
  if (is.null(row) || !("layout_asset_file" %in% names(row))) {
    return(NULL)
  }
  path <- as.character(row$layout_asset_file[[1L]])
  if (!nzchar(path) || is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  coords <- if (is.list(x) && !is.null(x$coords)) x$coords else x
  if (is.null(coords)) {
    return(NULL)
  }
  coords <- as.matrix(coords)
  if (nrow(coords) < 1L || ncol(coords) < 2L) {
    return(NULL)
  }
  list(coords = coords, file = path)
}

setting_from_index_row <- function(row) {
  data.frame(
    graph_family = row$graph_family,
    k = suppressWarnings(as.integer(row$k)),
    radius_rank = suppressWarnings(as.integer(row$radius_rank)),
    k_scale = suppressWarnings(as.integer(row$k_scale)),
    radius_rule = as.character(row$radius_rule),
    radius_factor = suppressWarnings(as.numeric(row$radius_factor)),
    delta = suppressWarnings(as.numeric(scalar_or_na(row$delta))),
    prune_method = row$prune_method,
    stage = row$default_stage,
    setting_id = row$setting_id,
    stringsAsFactors = FALSE
  )
}

make_radius_lookup <- function(ds, n) {
  radius.values <- stats::quantile(
    pairwise_nonzero_distances(ds$X_embed),
    probs = config$radius.rank / (n - 1),
    names = FALSE,
    type = 7
  )
  stats::setNames(as.numeric(radius.values), as.character(config$radius.rank))
}

precompute_widgets <- function(widget.index, max.rows, requested.key) {
  unlink(Sys.glob(file.path(widget.dir, "view_*.html")), recursive = TRUE)
  unlink(Sys.glob(file.path(widget.dir, "view_*_files")), recursive = TRUE)
  placeholders <- ensure_widget_placeholders()
  if (!nrow(widget.index)) return(widget.index)
  widget.index$original_widget <- relative_widget_path(placeholders$original)
  widget.index$grip_widget <- relative_widget_path(placeholders$grip)
  widget.index$widget_stage <- widget.index$default_stage
  widget.index$stage_unavailable_widget <- relative_widget_path(placeholders$stage)
  widget.index$missing_grip_widget <- relative_widget_path(placeholders$grip)

  if (!requireNamespace("rgl", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE) ||
      !isTRUE(args$layouts)) {
    return(widget.index)
  }

  dataset.ids <- unique(widget.index$dataset_id)
  for (dataset.id in dataset.ids) {
    ds.path <- dataset_cache_path(dataset.id)
    if (!file.exists(ds.path)) next
    ds <- readRDS(ds.path)
    out.file <- file.path(widget.dir, paste0("original_", safe_stem(dataset.id), ".html"))
    if (!file.exists(out.file)) {
      save_rgl_scene(
        ds$X_embed, out.file,
        paste(dataset.id, "original data"),
        color = grDevices::hcl.colors(nrow(ds$X_embed), "viridis")
      )
    }
    widget.index$original_widget[widget.index$dataset_id == dataset.id] <-
      relative_widget_path(out.file)
  }

  rows.to.render <- select_layout_rows(widget.index, max.rows, requested.key)
  if (!length(rows.to.render)) return(widget.index)
  for (ii in seq_along(rows.to.render)) {
    i <- rows.to.render[[ii]]
    row <- widget.index[i, , drop = FALSE]
    stage.to.render <- requested_render_stage(row$graph_key, row$default_stage,
                                              requested.key)
    payload <- load_graph_stage_asset(row$dataset_id, row$setting_id, stage.to.render)
    if (is.null(payload)) {
      ds.path <- dataset_cache_path(row$dataset_id)
      if (!file.exists(ds.path)) next
      ds <- readRDS(ds.path)
      setting <- setting_from_index_row(row)
      radius.lookup <- make_radius_lookup(ds, row$n)
      g <- tryCatch(build_graph(ds$X_embed, setting, radius.lookup),
                    error = function(e) NULL)
      if (is.null(g)) next
      payload <- tryCatch(stage_payload(g, stage.to.render), error = function(e) NULL)
    }
    if (is.null(payload) || is.null(payload$adj) || is.null(payload$weight)) next
    edges <- edge_table_from_adj(payload$adj, payload$weight, max.edges = 900L)
    key <- layout_stage_key(row$graph_key, stage.to.render)
    grip.file <- file.path(widget.dir, paste0(key, "_weighted_grip.html"))
    cache.file <- file.path(widget.cache.dir, paste0(key, ".rds"))
    layout <- if (file.exists(cache.file)) {
      readRDS(cache.file)
    } else {
      asset.layout <- load_layout_stage_asset(row$dataset_id, row$setting_id,
                                              stage.to.render)
      if (!is.null(asset.layout)) {
        asset.layout$coords
      } else {
        compute_grip_layout(payload$adj, payload$weight, seed = 2000L + ii)
      }
    }
    if (!is.null(layout)) {
      saveRDS(layout, cache.file)
      if (!file.exists(grip.file)) {
        save_rgl_scene(
          layout, grip.file,
          paste(row$graph_family, stage.to.render, "weighted GRIP"),
          edges = edges,
          color = rep("#111111", nrow(layout)),
          edge.color = "#d14b3f"
        )
      }
      widget.index$grip_cached[[i]] <- TRUE
      widget.index$grip_widget[[i]] <- relative_widget_path(grip.file)
    }
    widget.index$widget_stage[[i]] <- stage.to.render
  }
  widget.index
}

ok.metrics <- metrics[metrics$status == "ok", , drop = FALSE]
error.metrics <- metrics[metrics$status != "ok", , drop = FALSE]
best.by.group <- if (nrow(ok.metrics)) {
  split.idx <- split(seq_len(nrow(ok.metrics)),
                     paste(ok.metrics$surface, ok.metrics$n, ok.metrics$target, sep = "|"))
  do.call(rbind, lapply(split.idx, function(idx) {
    sub <- ok.metrics[idx, , drop = FALSE]
    sub[which.min(sub$rel_rms_error), , drop = FALSE]
  }))
} else data.frame()

family.summary <- if (nrow(ok.metrics)) {
  stats::aggregate(
    cbind(rel_rms_error, rel_abs_error_q95, pearson_cor,
          rel_geodesic_stress, signed_bias, shortcut_fraction,
          q50_rel_abs_residual, q90_rel_abs_residual, q95_rel_abs_residual,
          short_band_bias, mid_band_bias, long_band_bias) ~ target + surface + graph_family + prune_method,
    data = ok.metrics,
    FUN = function(x) stats::median(as.numeric(x), na.rm = TRUE),
    na.action = stats::na.pass
  )
} else data.frame()

diag.summary <- if (nrow(diagnostics) && "status" %in% names(diagnostics)) {
  ok.diag <- diagnostics[diagnostics$status == "ok", , drop = FALSE]
  if (nrow(ok.diag) && "n_mst_edges_added" %in% names(ok.diag)) {
    stats::aggregate(cbind(n_mst_edges_added, n_pruned_edges) ~ surface + graph_family + prune_method,
                     data = ok.diag,
                     FUN = function(x) median(as.numeric(x), na.rm = TRUE))
  } else data.frame()
} else data.frame()

figs <- list(
  sample_perf = plot_family_performance(metrics, file.path(fig.dir, "sample_oracle_method_performance.png"), "sample_oracle"),
  surface_perf = plot_family_performance(metrics, file.path(fig.dir, "surface_method_performance.png"), "surface"),
  sample_dist = plot_rank_distribution(metrics, file.path(fig.dir, "sample_oracle_error_distribution.png"), "sample_oracle"),
  surface_dist = plot_rank_distribution(metrics, file.path(fig.dir, "surface_error_distribution.png"), "surface"),
  sample_prune = plot_pruning_effect(metrics, file.path(fig.dir, "sample_oracle_pruning_effect.png"), "sample_oracle"),
  surface_prune = plot_pruning_effect(metrics, file.path(fig.dir, "surface_pruning_effect.png"), "surface"),
  diagnostics = plot_diagnostics(diagnostics, file.path(fig.dir, "connectivity_pruning_diagnostics.png"))
)

widget.index <- build_widget_index(metrics, diagnostics)
widget.index <- precompute_widgets(widget.index, args$max.layout.rows, args$layout.graph.key)
jsonlite::write_json(widget.index, file.path(report.dir, "widget_index.json"),
                     dataframe = "rows", pretty = TRUE, auto_unbox = TRUE,
                     na = "null")

rel_path <- function(path) {
  if (is.null(path) || !length(path)) return(character())
  if (startsWith(path, report.dir)) {
    return(htmltools::htmlEscape(sub(paste0("^", report.dir, "/?"), "", path)))
  }
  htmltools::htmlEscape(path)
}

img_tag <- function(path, alt) {
  if (is.null(path) || !file.exists(path)) return("")
  sprintf("<figure><img src=\"%s\" alt=\"%s\" /><figcaption>%s</figcaption></figure>",
          rel_path(path), escape_html(alt), escape_html(alt))
}

cols.best <- intersect(c("surface", "n", "graph_family", "k", "radius_rank", "k_scale",
                         "radius_rule", "radius_factor", "prune_method", "stage",
                         "rel_rms_error", "rel_abs_error_q95", "pearson_cor",
                         geodesic.metric.cols),
                       names(best.by.group))
cols.family <- intersect(c("surface", "graph_family", "prune_method", "rel_rms_error",
                           "rel_abs_error_q95", "pearson_cor", geodesic.metric.cols),
                         names(family.summary))

best_sample <- best.by.group[best.by.group$target == "sample_oracle", cols.best, drop = FALSE]
best_surface <- best.by.group[best.by.group$target == "surface", cols.best, drop = FALSE]
fam_sample <- family.summary[family.summary$target == "sample_oracle", cols.family, drop = FALSE]
fam_surface <- family.summary[family.summary$target == "surface", cols.family, drop = FALSE]
if (nrow(fam_sample)) fam_sample <- fam_sample[order(fam_sample$surface, fam_sample$rel_rms_error), ]
if (nrow(fam_surface)) fam_surface <- fam_surface[order(fam_surface$surface, fam_surface$rel_rms_error), ]

widget_controls <- if (nrow(widget.index)) {
  data.js <- jsonlite::toJSON(widget.index, dataframe = "rows", auto_unbox = TRUE,
                              na = "null")
  paste0(
    "<div class=\"controls\" id=\"graphControls\">",
    "<div><label for=\"presetSel\">Preset</label><select id=\"presetSel\"></select></div>",
    "<div><label for=\"surfaceSel\">Surface</label><select id=\"surfaceSel\"></select></div>",
    "<div><label for=\"nSel\">n</label><select id=\"nSel\"></select></div>",
    "<div><label for=\"seedSel\">seed</label><select id=\"seedSel\"></select></div>",
    "<div><label for=\"familySel\">Graph family</label><select id=\"familySel\"></select></div>",
    "<div><label for=\"pruneSel\">Pruning</label><select id=\"pruneSel\"></select></div>",
    "<div><label for=\"stageSel\">Stage</label><select id=\"stageSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"sknn,mknn,iknn\"><label for=\"kSel\">k</label><select id=\"kSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"fixed_radius\"><label for=\"radiusRankSel\">radius rank</label><select id=\"radiusRankSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"adaptive_radius\"><label for=\"kScaleSel\">k scale</label><select id=\"kScaleSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"adaptive_radius\"><label for=\"radiusRuleSel\">radius rule</label><select id=\"radiusRuleSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"adaptive_radius\"><label for=\"radiusFactorSel\">radius factor</label><select id=\"radiusFactorSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"cknn\"><label for=\"cknnKScaleSel\">k scale</label><select id=\"cknnKScaleSel\"></select></div>",
    "<div class=\"param-group\" data-family=\"cknn\"><label for=\"deltaSel\">delta</label><select id=\"deltaSel\"></select></div>",
    "</div>",
    "<div id=\"widgetMeta\" class=\"note\"></div>",
    "<div class=\"widget-grid\"><figure><figcaption>Original data</figcaption><iframe id=\"origFrame\" class=\"widget\"></iframe></figure><figure><figcaption>Weighted GRIP layout</figcaption><iframe id=\"gripFrame\" class=\"widget\"></iframe></figure></div>",
    "<script>const GRAPH_INDEX=", data.js, ";
const PRESETS=['manual','best_by_surface','best_by_sample_oracle','median_by_surface','median_by_sample_oracle','worst_by_surface','worst_by_sample_oracle'];
const STAGES=['raw','raw.repaired','pruned','pruned.repaired','repaired.pruned','final'];
const LABELS={manual:'manual',best_by_surface:'best by surface',best_by_sample_oracle:'best by sample oracle',median_by_surface:'median by surface',median_by_sample_oracle:'median by sample oracle',worst_by_surface:'worst by surface',worst_by_sample_oracle:'worst by sample oracle',adaptive_radius:'adaptive radius',cknn:'continuous kNN',fixed_radius:'fixed radius',sknn:'sKNN',mknn:'mKNN',iknn:'iKNN','global.geodesic.ratio':'global geodesic ratio','local.geodesic':'local geodesic','none':'none'};
const ids=['surfaceSel','nSel','seedSel','familySel','pruneSel','stageSel','kSel','radiusRankSel','kScaleSel','radiusRuleSel','radiusFactorSel','cknnKScaleSel','deltaSel'];
function el(id){return document.getElementById(id);}
function present(x){return x!==null && x!==undefined && x!=='' && String(x)!=='NA';}
function num(x){return present(x) && isFinite(Number(x)) ? Number(x).toFixed(5) : 'NA';}
function fmt(x){return present(x) ? String(x) : 'NA';}
function uniq(xs){return Array.from(new Set(xs.filter(present).map(String))).sort((a,b)=>a.localeCompare(b,undefined,{numeric:true}));}
function setOptions(sel, values, labels){
  const old=sel.value; sel.innerHTML='';
  values.forEach(v=>{const o=document.createElement('option'); o.value=String(v); o.textContent=(labels&&labels[v])?labels[v]:String(v); sel.appendChild(o);});
  if(values.map(String).includes(old)) sel.value=old;
}
function familyParams(row){
  if(row.graph_family==='sknn'||row.graph_family==='mknn'||row.graph_family==='iknn') return ['k'];
  if(row.graph_family==='fixed_radius') return ['radius_rank'];
  if(row.graph_family==='adaptive_radius') return ['k_scale','radius_rule','radius_factor'];
  if(row.graph_family==='cknn') return ['k_scale','delta'];
  return [];
}
function wanted(){
  const fam=el('familySel').value;
  return {surface:el('surfaceSel').value,n:el('nSel').value,seed:el('seedSel').value,graph_family:fam,prune_method:el('pruneSel').value,k:el('kSel').value,radius_rank:el('radiusRankSel').value,k_scale:fam==='cknn'?el('cknnKScaleSel').value:el('kScaleSel').value,radius_rule:el('radiusRuleSel').value,radius_factor:el('radiusFactorSel').value,delta:el('deltaSel').value};
}
function rowMatches(row, ignore){
  const w=wanted();
  for(const key of ['surface','n','seed','graph_family','prune_method']){
    if(key===ignore) continue;
    if(present(w[key]) && String(row[key])!==String(w[key])) return false;
  }
  for(const key of familyParams(row)){
    if(key===ignore) continue;
    if(present(w[key]) && String(row[key])!==String(w[key])) return false;
  }
  return true;
}
function candidates(ignore){return GRAPH_INDEX.filter(r=>rowMatches(r,ignore));}
function refreshOptions(){
  setOptions(el('surfaceSel'), uniq(candidates('surface').map(r=>r.surface)));
  setOptions(el('nSel'), uniq(candidates('n').map(r=>r.n)));
  setOptions(el('seedSel'), uniq(candidates('seed').map(r=>r.seed)));
  setOptions(el('familySel'), uniq(candidates('graph_family').map(r=>r.graph_family)), LABELS);
  setOptions(el('pruneSel'), uniq(candidates('prune_method').map(r=>r.prune_method)), LABELS);
  setOptions(el('stageSel'), STAGES);
  setOptions(el('kSel'), uniq(candidates('k').map(r=>r.k)));
  setOptions(el('radiusRankSel'), uniq(candidates('radius_rank').map(r=>r.radius_rank)));
  setOptions(el('kScaleSel'), uniq(candidates('k_scale').map(r=>r.k_scale)));
  setOptions(el('radiusRuleSel'), uniq(candidates('radius_rule').map(r=>r.radius_rule)));
  setOptions(el('radiusFactorSel'), uniq(candidates('radius_factor').map(r=>r.radius_factor)));
  setOptions(el('cknnKScaleSel'), uniq(candidates('k_scale').map(r=>r.k_scale)));
  setOptions(el('deltaSel'), uniq(candidates('delta').map(r=>r.delta)));
  updateParamVisibility();
}
function updateParamVisibility(){
  const fam=el('familySel').value;
  document.querySelectorAll('.param-group').forEach(g=>{
    const families=g.getAttribute('data-family').split(',');
    const show=families.includes(fam);
    g.style.display=show?'block':'none';
    g.querySelectorAll('select').forEach(s=>s.disabled=!show);
  });
}
function setControlsFromRow(row){
  el('surfaceSel').value=String(row.surface); el('nSel').value=String(row.n);
  el('seedSel').value=String(row.seed); el('familySel').value=String(row.graph_family);
  el('pruneSel').value=String(row.prune_method);
  refreshOptions();
  if(present(row.k)) el('kSel').value=String(row.k);
  if(present(row.radius_rank)) el('radiusRankSel').value=String(row.radius_rank);
  if(present(row.k_scale)) el('kScaleSel').value=String(row.k_scale);
  if(present(row.k_scale)) el('cknnKScaleSel').value=String(row.k_scale);
  if(present(row.radius_rule)) el('radiusRuleSel').value=String(row.radius_rule);
  if(present(row.radius_factor)) el('radiusFactorSel').value=String(row.radius_factor);
  if(present(row.delta)) el('deltaSel').value=String(row.delta);
  el('stageSel').value=String(row.default_stage);
}
function applyPreset(){
  const preset=el('presetSel').value;
  if(preset==='manual') return;
  const target=preset.includes('sample_oracle')?'sample_oracle':'surface';
  const metric=target+'_rel_rms_error';
  let rows=GRAPH_INDEX.filter(r=>String(r.surface)===String(el('surfaceSel').value)&&String(r.n)===String(el('nSel').value)&&String(r.seed)===String(el('seedSel').value)&&present(r[metric]));
  if(!rows.length) rows=GRAPH_INDEX.filter(r=>present(r[metric]));
  rows.sort((a,b)=>Number(a[metric])-Number(b[metric]));
  let row=rows[0];
  if(preset.startsWith('worst')) row=rows[rows.length-1];
  if(preset.startsWith('median')) row=rows[Math.floor((rows.length-1)/2)];
  if(row) setControlsFromRow(row);
}
function exactRow(){
  const rows=GRAPH_INDEX.filter(r=>rowMatches(r,''));
  return rows.length ? rows[0] : null;
}
function stageName(stage){return stage.replaceAll('.','_');}
function stageValue(row,prefix,stage){const key=prefix+'_'+stageName(stage); return present(row[key]) ? row[key] : null;}
function renderWidget(){
  refreshOptions();
  const row=exactRow();
  if(!row){
    el('widgetMeta').innerHTML='No graph setting matches the current controls.';
    return;
  }
  const stage=el('stageSel').value || row.default_stage;
  const available=String(row.available_stages||'').split(',').filter(Boolean);
  const stageOK=available.includes(stage);
  const stageCached=String(stage)===String(row.widget_stage);
  el('origFrame').src=row.original_widget;
  el('gripFrame').src=stageOK && stageCached && row.grip_cached ? row.grip_widget : (stageOK ? row.missing_grip_widget : row.stage_unavailable_widget);
  const edges=stageValue(row,'n_edges',stage);
  const comps=stageValue(row,'n_components',stage);
  const layoutKey=row.graph_key+'__stage_'+stage.replace(/[^A-Za-z0-9]+/g,'_');
  let note = '<strong>'+LABELS[row.graph_family]+'</strong> on '+row.surface+', n='+row.n+', seed='+row.seed+
    '<br><strong>selected parameters:</strong> k='+fmt(row.k)+', radius_rank='+fmt(row.radius_rank)+', k_scale='+fmt(row.k_scale)+', radius_rule='+fmt(row.radius_rule)+', radius_factor='+fmt(row.radius_factor)+', delta='+fmt(row.delta)+
    '<br><strong>pruning:</strong> '+LABELS[row.prune_method]+'; <strong>stage:</strong> '+stage+
    '<br><strong>surface rel_rms_error:</strong> '+num(row.surface_rel_rms_error)+'; <strong>sample-oracle rel_rms_error:</strong> '+num(row.sample_oracle_rel_rms_error)+
    '<br><strong>surface signed bias:</strong> '+num(row.surface_signed_bias)+'; <strong>surface shortcut fraction:</strong> '+num(row.surface_shortcut_fraction)+
    '<br><strong>n_edges:</strong> '+fmt(edges)+'; <strong>n_components:</strong> '+fmt(comps)+'; <strong>n_mst_edges_added:</strong> '+fmt(row.n_mst_edges_added)+'; <strong>n_pruned_edges:</strong> '+fmt(row.n_pruned_edges)+
    '<br><strong>graph_key:</strong> <code>'+row.graph_key+'</code>; <strong>layout key:</strong> <code>'+layoutKey+'</code>';
  if(!stageOK) note += '<br><span class=\"warn-inline\">The selected lifecycle stage is unavailable for this graph object.</span>';
  else if(!row.grip_cached || !stageCached) note += '<br><span class=\"warn-inline\">Weighted GRIP layout is not cached for this selected graph/stage. Run the on-demand layout command shown below.</span>';
  el('widgetMeta').innerHTML=note;
}
function manualChange(){el('presetSel').value='manual'; renderWidget();}
function initWidgets(){
  setOptions(el('presetSel'), PRESETS, LABELS);
  setOptions(el('surfaceSel'), uniq(GRAPH_INDEX.map(r=>r.surface)));
  if(uniq(GRAPH_INDEX.map(r=>r.surface)).includes('paraboloid')) el('surfaceSel').value='paraboloid';
  setOptions(el('nSel'), uniq(GRAPH_INDEX.filter(r=>String(r.surface)===String(el('surfaceSel').value)).map(r=>r.n)));
  setOptions(el('seedSel'), uniq(GRAPH_INDEX.filter(r=>String(r.surface)===String(el('surfaceSel').value)&&String(r.n)===String(el('nSel').value)).map(r=>r.seed)));
  setOptions(el('familySel'), uniq(GRAPH_INDEX.map(r=>r.graph_family)), LABELS);
  if(uniq(GRAPH_INDEX.map(r=>r.graph_family)).includes('adaptive_radius')) el('familySel').value='adaptive_radius';
  setOptions(el('pruneSel'), uniq(GRAPH_INDEX.map(r=>r.prune_method)), LABELS);
  setOptions(el('stageSel'), STAGES);
  refreshOptions();
  el('presetSel').value='best_by_surface';
  applyPreset();
  el('presetSel').addEventListener('change',()=>{applyPreset();renderWidget();});
  ids.forEach(id=>el(id).addEventListener('change',manualChange));
  renderWidget();
}
window.addEventListener('load', initWidgets);
</script>")
} else {
  "<p>No rglwidget views were generated for this run.</p>"
}

timeout.text <- if (args$setting.timeout.sec <= 0L || !is.finite(args$setting.timeout.sec)) {
  "No elapsed-time watchdog was used."
} else {
  paste0("A per-setting elapsed-time watchdog of ", args$setting.timeout.sec,
         " seconds was used; rows exceeding it are marked as errors.")
}

run.command <- paste(
  "Rscript dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R",
  paste0("--mode=", args$mode),
  paste0("--output.dir=", base.dir),
  paste0("--grid.size=", config$grid.size),
  paste0("--workers=", args$workers),
  paste0("--setting.timeout.sec=", args$setting.timeout.sec),
  paste0("--layouts=", tolower(as.character(args$layouts))),
  paste0("--save.graph.assets=", tolower(as.character(args$save.graph.assets))),
  paste0("--save.layout.assets=", tolower(as.character(args$save.layout.assets))),
  paste0("--asset.stages=", args$asset.stages)
)

prune.methods.in.results <- sort(unique(as.character(ok.metrics$prune_method)))
prune.methods.in.results <- prune.methods.in.results[nzchar(prune.methods.in.results)]
prune.method.html <- paste(sprintf("<code>%s</code>", escape_html(prune.methods.in.results)),
                           collapse = ", ")
is.pruning.comparison.report <- "global.geodesic.ratio" %in% prune.methods.in.results ||
  grepl("pruning-method-comparison", normalizePath(base.dir, mustWork = FALSE))
report.title <- if (isTRUE(is.pruning.comparison.report)) {
  "Quadratic-Surface Pruning-Method Comparison Benchmark"
} else {
  "Quadratic-Surface Data-to-Graph Geodesic Reconstruction Benchmark"
}
context.paragraph <- if (isTRUE(is.pruning.comparison.report)) {
  paste0(
    "<p>This pruning-method comparison reuses the quadratic-surface data-to-graph benchmark from ",
    "<code>dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md</code>. ",
    "It asks how ", prune.method.html,
    " pruning choices affect graph shortest-path fidelity to the intrinsic geodesic metric of the sampled surface. ",
    "Runtime, connectivity, and edge sparsity are reported as supporting diagnostics, while relative graph-geodesic error is the primary comparison.</p>"
  )
} else {
  paste0(
    "<p>This report follows <code>dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md</code>. ",
    "The data geodesic geometric reconstruction problem asks how to turn a finite sample <code>X</code> into a weighted graph whose shortest-path metric is close to the intrinsic geodesic metric of the sampled space. ",
    "This matters when the ambient Euclidean geometry is misleading, for example on curved, folded, branched, or nearly self-intersecting data. ",
    "The pruning choices represented in this completed run are ", prune.method.html, ".</p>"
  )
}
diagnostic.paragraph <- paste0(
  "<p>The benchmark uses MST component repair before evaluation. ",
  "The pruning choices represented in this run are ", prune.method.html,
  ". The diagnostic figure shows how many MST bridge edges and pruned edges are typical for each family.</p>"
)

html <- paste0("<!doctype html><html><head><meta charset=\"utf-8\"><title>", escape_html(report.title), "</title>",
"<script>window.MathJax={tex:{inlineMath:[[\"\\\\(\",\"\\\\)\"],[\"$\",\"$\"]],displayMath:[[\"\\\\[\",\"\\\\]\"],[\"$$\",\"$$\"]]}};</script><script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>",
"<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:34px;line-height:1.52;color:#20242a;background:#fff}h1,h2{color:#16324f}h3{color:#243b53;margin-top:24px}table{border-collapse:collapse;width:100%;margin:14px 0 30px;font-size:13px}th,td{border:1px solid #d8dde6;padding:6px 8px;text-align:left;vertical-align:top}th{background:#eef2f6}code{background:#f2f4f7;padding:1px 4px;border-radius:3px}.math-block{overflow-x:auto;background:#fbfcfe;border:1px solid #e5e9f0;padding:10px 14px;margin:12px 0}.note{background:#f7f9fc;border-left:4px solid #4f77aa;padding:12px 14px;margin:18px 0}.warn{background:#fff7ed;border-left:4px solid #d97706;padding:12px 14px;margin:18px 0}.warn-inline{color:#a55307;font-weight:600}img{max-width:100%;border:1px solid #d8dde6;margin:8px 0 6px}.controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;align-items:end;margin:16px 0;padding:12px;background:#f6f8fa;border:1px solid #d9e2ec}.controls label{display:block;font-weight:600;font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:#52616f}.controls select{box-sizing:border-box;width:100%;padding:7px;border:1px solid #bcccdc;background:white}.widget-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}.widget-grid.three{grid-template-columns:1fr 1fr 1fr}.widget-grid figure{margin:0}.widget{width:100%;height:520px;border:1px solid #d8dde6}.metric{font-weight:600}.small{color:#5b6470;font-size:13px}figure{margin:14px 0 26px}figcaption{font-size:13px;color:#4b5563}@media(max-width:1300px){.widget-grid.three{grid-template-columns:1fr}.widget{height:520px}}</style></head><body>",
"<h1>", escape_html(report.title), "</h1>",
"<p><strong>Generated:</strong> ", escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "<br><strong>Mode:</strong> ", escape_html(args$mode), "<br><strong>Grid size:</strong> ", config$grid.size, "<br><strong>Workers:</strong> ", args$workers, "<br><strong>Run directory:</strong> <code>", escape_html(run.dir), "</code></p>",
"<h2>Context and Question</h2>", context.paragraph,
"<div class=\"math-block\">\\[G(X)=(V,E,\\ell),\\qquad V=X,\\qquad d_G(x_i,x_j)=\\min_{\\gamma:i\\leadsto j}\\sum_{[u,v]\\in\\gamma}\\ell_{uv}.\\]</div>",
"<p>The benchmark asks which graph construction gives <code>d_G</code> closest to the reference geodesic geometry for two quadratic graph surfaces:</p><div class=\"math-block\">\\[q(u,v)=u^2+v^2\\qquad\\text{and}\\qquad q(u,v)=u^2-v^2.\\]</div>",
"<p>All graph edges in this benchmark use ambient Euclidean lengths in the embedded sample, <span class=\"math-inline\">\\(\\ell_{ij}=\\|x_i-x_j\\|_2\\)</span>. The graph family and parameter controls therefore change the edge support, not the edge-length convention.</p>",
"<h2>Evaluation Targets</h2><p>Each graph is compared with two targets. The <strong>surface</strong> target is the numerical continuum geodesic distance on the quadratic surface:</p><div class=\"math-block\">\\[D_{\\mathrm{surface}}=\\bigl(d_M(x_i,x_j)\\bigr)_{ij}.\\]</div><p>The <strong>sample-oracle</strong> target follows the estimated continuum geodesic and measures the best route available through nearby observed sample points:</p><div class=\"math-block\">\\[D_{\\mathrm{oracle}}=\\bigl(d_X^{\\mathrm{oracle}}(x_i,x_j)\\bigr)_{ij}.\\]</div><p>These targets affect ranking and preset selection only. They do not change graph construction, graph pruning, graph repair, edge lengths, or weighted-GRIP visualization.</p>",
"<h2>Metric Definitions</h2><p>The main score is a scale-calibrated relative error. Given an estimated graph distance matrix <span class=\"math-inline\">\\(D_G\\)</span> and a reference matrix <span class=\"math-inline\">\\(D_{\\mathrm{ref}}\\)</span>, the scalar calibration is</p><div class=\"math-block\">\\[\\alpha^*=\\arg\\min_{\\alpha>0}\\sum_{i<j}\\bigl(\\alpha D_G(i,j)-D_{\\mathrm{ref}}(i,j)\\bigr)^2.\\]</div><p>The relative RMS error is</p><div class=\"math-block\">\\[\\operatorname{rel\\_rms\\_error}=\\left(\\frac{\\sum_{i<j}\\bigl(\\alpha^*D_G(i,j)-D_{\\mathrm{ref}}(i,j)\\bigr)^2}{\\sum_{i<j}D_{\\mathrm{ref}}(i,j)^2}\\right)^{1/2}.\\]</div><p>The tables also report relative absolute error quantiles, computed from <span class=\"math-inline\">\\(|\\alpha^*D_G(i,j)-D_{\\mathrm{ref}}(i,j)|/D_{\\mathrm{ref}}(i,j)</span>, and distortion quantiles, computed from <span class=\"math-inline\">\\(\\alpha^*D_G(i,j)/D_{\\mathrm{ref}}(i,j)</span>. Lower relative error is better; distortion values closer to one are better.</p>",
"<h2>Execution Summary</h2><p>Configured datasets: ", nrow(dataset_manifest), ". Datasets represented in current results: ", length(unique(metrics$dataset_id)), ". Metric rows: ", nrow(metrics), ". Successful metric rows: ", nrow(ok.metrics), ". Error metric rows: ", nrow(error.metrics), ".</p>",
"<p>A <span class=\"metric\">successful metric row</span> is one row in <code>metrics.csv</code> with <code>status='ok'</code>: graph construction completed, the selected lifecycle stage was available, graph shortest paths were computed, and the resulting distance matrix was compared with one evaluation target. An <span class=\"metric\">error metric row</span> has <code>status='error'</code> and is retained for bookkeeping with the graph setting, target, and error message. ", escape_html(timeout.text), "</p>",
if (nrow(error.metrics)) paste0("<div class=\"warn\"><p>Some settings did not complete. The first distinct errors are shown below.</p>", table_html(unique(error.metrics[, c("dataset_id", "graph_family", "setting_id", "error")]), max.rows = 12), "</div>") else "",
"<h2>Sample-Oracle Target: Method Performance</h2><p>This target asks whether the graph metric matches the best finite-sample path suggested by the true surface geodesic. Lower relative RMS error is better.</p>",
img_tag(figs$sample_perf, "Median method performance against the sample-oracle target"), img_tag(figs$sample_dist, "Distribution of sample-oracle relative RMS errors"),
"<h3>Best Settings: Sample-Oracle Target</h3>", table_html(best_sample, digits = 5, max.rows = 80),
"<h3>Family-Level Median Performance: Sample-Oracle Target</h3>", table_html(fam_sample, digits = 5, max.rows = 120),
"<h2>Surface Target: Method Performance</h2><p>This target asks whether the graph metric recovers the continuum surface geodesic distance. It is the cleanest measure of data-to-graph geodesic reconstruction when the underlying surface is known.</p>",
img_tag(figs$surface_perf, "Median method performance against the surface target"), img_tag(figs$surface_dist, "Distribution of surface-target relative RMS errors"),
"<h3>Best Settings: Surface Target</h3>", table_html(best_surface, digits = 5, max.rows = 80),
"<h3>Family-Level Median Performance: Surface Target</h3>", table_html(fam_surface, digits = 5, max.rows = 120),
"<h2>Geodesic-Isometry Diagnostics</h2>",
"<p>These diagnostics summarize calibrated pairwise graph-distance residuals. Signed bias is the mean signed residual divided by the mean reference distance; negative values indicate systematic shortcuts, positive values indicate systematic detours. Shortcut fraction is the fraction of calibrated graph distances below the reference. The short, mid, and long band biases are median signed relative residuals over reference-distance tertiles.</p>",
"<p>For rows inherited from older cached runs, <code>rel_geodesic_stress</code> is backfilled from <code>rel_rms_error</code>. Directional residual diagnostics require pairwise residuals and appear only for settings recomputed with the current package helper.</p>",
"<h3>Sample-Oracle Diagnostics</h3>",
table_html(fam_sample[, intersect(c("surface", "graph_family", "prune_method",
                                    geodesic.metric.cols),
                                  names(fam_sample)), drop = FALSE],
           digits = 5, max.rows = 120),
"<h3>Surface Diagnostics</h3>",
table_html(fam_surface[, intersect(c("surface", "graph_family", "prune_method",
                                     geodesic.metric.cols),
                                   names(fam_surface)), drop = FALSE],
           digits = 5, max.rows = 120),
"<h2>Pruning and Connectivity Diagnostics</h2>", diagnostic.paragraph,
img_tag(figs$sample_prune, "Pruning effect on the sample-oracle target"), img_tag(figs$surface_prune, "Pruning effect on the surface target"), img_tag(figs$diagnostics, "Connectivity repair and pruning diagnostics"), table_html(diag.summary, digits = 3, max.rows = 100),
"<h2>Interactive 3D Graph Diagnostics</h2><p>The controls below expose the graph-construction layer directly. The <strong>preset</strong> control jumps to useful settings, but <strong>manual</strong> mode can select every graph setting represented in the benchmark results. The <strong>target</strong> is deliberately absent from these controls because evaluation targets rank graph settings; they are not graph-construction parameters.</p><p>The single 3D row has two panels: original data and a weighted-GRIP layout of the selected graph stage. Points are rendered as spheres. Edges are sampled when necessary to keep widgets responsive. If a selected weighted-GRIP layout is missing, the report shows the graph key needed to cache it.</p>", widget_controls,
"<div class=\"note\"><p>On-demand layout refresh command:</p><pre><code>", escape_html(paste0("Rscript dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R --mode=", args$mode, " --output.dir=", base.dir, " --report.only=true --layouts=true --layout.graph.key=<graph_key>")), "</code></pre></div>",
"<h2>Artifacts and Reproducibility</h2><ul><li><code>results.rds</code>: metrics, graph diagnostics, asset manifests, and run config.</li><li><code>metrics.csv</code>: one row per graph setting and target.</li><li><code>graph_diagnostics.csv</code>: edge/component/bridge/pruning diagnostics.</li><li><code>dataset_manifest.csv</code>: generated datasets.</li><li><code>dataset_assets.csv</code>: dataset asset files used by interactive explorers.</li><li><code>graph_assets.csv</code>: graph-stage <code>adj_list</code>/<code>weight_list</code> asset files.</li><li><code>layout_assets.csv</code>: weighted-GRIP 3D layout asset files.</li><li><code>quadform_benchmark_manifest.rds</code> and <code>quadform_benchmark_manifest.json</code>: consolidated benchmark asset manifest for downstream tools such as <code>gflowui</code>.</li><li><code>run_config.json</code>: parameter grid and run settings.</li><li><code>report/widget_index.json</code>: compact JSON index used by the static 3D controls.</li><li><code>report/widgets/*.html</code>: generated rglwidget and placeholder panels.</li><li><code>report/widget_cache/*.rds</code>: static-report widget cache, separate from benchmark layout assets.</li></ul><p><strong>Command shape:</strong></p><pre><code>", escape_html(run.command), "</code></pre>",
"</body></html>")

report.path <- file.path(report.dir, "quadform_first_benchmark_report.html")
writeLines(html, report.path)
index.path <- file.path(base.dir, paste0("quadform_first_benchmark_", args$mode, "_report.html"))
index.html <- paste0("<!doctype html><html><head><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"0; url=", file.path("runs", args$mode, "report", "quadform_first_benchmark_report.html"), "\"><title>", escape_html(report.title), "</title></head><body><p>Open <a href=\"", file.path("runs", args$mode, "report", "quadform_first_benchmark_report.html"), "\">the benchmark report</a>.</p></body></html>")
writeLines(index.html, index.path)
log_msg("Report written to ", report.path)
log_msg("Index written to ", index.path)

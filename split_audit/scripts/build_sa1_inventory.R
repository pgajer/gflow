#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run from the gflow repository root.")
}

path_exists <- function(path) file.exists(file.path(repo_root, path))
pseudo_asset <- function(path) grepl("\\*$|/$", path)

list_assets <- function(dir, pattern = NULL, max_depth = 1L) {
  root <- file.path(repo_root, dir)
  if (!dir.exists(root)) {
    return(character())
  }
  files <- list.files(root, recursive = max_depth > 1L, full.names = FALSE, all.files = FALSE)
  if (max_depth > 1L) {
    depth <- lengths(gregexpr("/", files, fixed = TRUE))
    depth[files == ""] <- 0L
    files <- files[depth < max_depth]
  }
  files <- file.path(dir, files)
  files <- files[file.info(file.path(repo_root, files))$isdir %in% FALSE]
  if (!is.null(pattern)) {
    files <- files[grepl(pattern, files)]
  }
  sort(files)
}

read_namespace_exports <- function() {
  ns <- file.path(repo_root, "NAMESPACE")
  if (!file.exists(ns)) {
    return(character())
  }
  lines <- readLines(ns, warn = FALSE)
  exported <- sub("^export\\((.*)\\)$", "\\1", lines[grepl("^export\\(", lines)])
  sort(exported)
}

exports <- read_namespace_exports()

asset_exports <- function(path) {
  full <- file.path(repo_root, path)
  if (!file.exists(full) || !grepl("^R/.*\\.R$", path)) {
    return("")
  }
  txt <- readLines(full, warn = FALSE)
  roxy <- grep("^#'\\s*@export", txt)
  if (!length(roxy)) {
    return("")
  }
  candidates <- character()
  for (i in roxy) {
    window <- txt[seq.int(i + 1L, min(length(txt), i + 20L))]
    hit <- grep("^[A-Za-z.][A-Za-z0-9._]*\\s*(<-|=)\\s*function\\s*\\(", window, value = TRUE)
    if (length(hit)) {
      candidates <- c(candidates, sub("^([A-Za-z.][A-Za-z0-9._]*)\\s*(<-|=).*", "\\1", hit[1]))
    }
  }
  paste(sort(unique(candidates)), collapse = ";")
}

classify_asset <- function(path) {
  p <- tolower(path)
  category <- "undecided"
  rationale <- "Needs manual review."
  confidence <- "low"

  set <- function(cat, why, conf = "medium") {
    category <<- cat
    rationale <<- why
    confidence <<- conf
  }

  if (grepl("(^|/)rcppexports\\.|(^|/)init\\.c$|makevars|compat|globals|utils|stats_utils|cpp_utils|memory_utils|progress_utils|io_utils|kernel_utils|kernels|local_pca|knn_cache|knn_metric|row_wise|random_sampling|synthetic_data|spline_replacements|inst/include|man/\\*\\.rd", p)) {
    set("shared", "Cross-cutting infrastructure, generated registration, utilities, synthetic/test helpers, or shared numerical support.", "medium")
  }
  if (grepl("graph|geodesic|iknn|rknn|mknn|sknn|cknn|radius|mst|shortest_path|dijkstra|isometry|packing|nerve|phate|diffusion|endpoint|spectrum|embedding|pruning|distance_graph", p)) {
    set("gflow", "Graph, geodesic, geometry, dimensionality, or gradient-flow infrastructure.", "medium")
  }
  if (grepl("gfc|gradient_flow|basin|morse|trajectory|extrema|separatrice|cell_|gflow_cx|gflow_graph", p)) {
    set("gflow", "Gradient-flow cells, Morse-Smale objects, basins, extrema, or trajectory infrastructure.", "high")
  }
  if (grepl("kernel_local_polynomial|lpl_tf|slpl_tf|malps|ssrhe|graph_trend_filtering|metric_graph_lowpass|harmonic_smoother|harmonic_extension|mean_shift_smoother|data_smoother|riem_dcx|rdgraph|pttf|lowess|regression|smoother|smooth|conditional|ulogit", p)) {
    set("geosmooth", "Conditional-expectation, smoothing, regression, trend-filtering, or fitted-response model.", "high")
  }
  if (grepl("assoc|lcor|lslope|concordance|permutation|wasserstein|clustering|cluster|cst|phylotype|trajectory_hurdle|metabolon|visualize_mgcp|isa_|bicluster|two_factor|compositional|comono|outlier|dbscan|consensus", p)) {
    set("gflowx", "Analysis workflow, association testing, clustering, visualization, microbiome-specific, or experimental application layer.", "medium")
  }
  if (grepl("agemalo|pgmalo|pgmalog|uggmalo|malo|klaps|adaptive_grid|ray_agemalo", p)) {
    set("gflowx", "Older or experimental conditional-expectation method family; candidate for GitHub-only archive unless still central.", "high")
  }
  if (grepl("description|namespace|license|copyright|wordlist|readme|news|vignette|examples", p)) {
    set("shared", "Package-level metadata, examples, or documentation that must be reassigned after package boundaries are finalized.", "medium")
  }
  if (grepl("^tests/testthat/test-(kernel-local-polynomial|lpl-tf|slpl-tf|malps|ssrhe|graph-trend-filtering|metric-graph-lowpass|pttf)", p)) {
    set("geosmooth", "Test coverage for smoothing/regression/trend-filtering method.", "high")
  }
  if (grepl("^tests/testthat/test-(geodesic|graph|radius|mknn|sknn|mst|nerve|phate|isometry|distance-graph|component-mst|quadform|local-geodesic)", p)) {
    set("gflow", "Test coverage for graph/geodesic/geometry infrastructure.", "high")
  }
  if (grepl("^tests/testthat/helper-", p)) {
    set("shared", "Test helper; final home depends on tests retained by each split package.", "medium")
  }

  list(category = category, rationale = rationale, confidence = confidence)
}

line_count <- function(path) {
  if (pseudo_asset(path)) {
    return(NA_integer_)
  }
  full <- file.path(repo_root, path)
  if (!file.exists(full)) {
    return(NA_integer_)
  }
  length(readLines(full, warn = FALSE))
}

assets <- c(
  "DESCRIPTION",
  "NAMESPACE",
  "man/*.Rd",
  "inst/include/Eigen/",
  "inst/include/libqhull/",
  list_assets("R", "\\.R$"),
  list_assets("src", "\\.(c|cpp|h|hpp|win|local|Makevars)$"),
  list_assets("tests/testthat", "\\.R$"),
  list_assets("vignettes", "\\.(Rmd|html)$"),
  list_assets("inst", NULL, max_depth = 2L)
)
assets <- sort(unique(assets[path_exists(assets) | pseudo_asset(assets)]))

rows <- lapply(assets, function(path) {
  cinfo <- classify_asset(path)
  data.frame(
    asset = path,
    asset_type = if (path == "DESCRIPTION") "metadata" else if (path == "NAMESPACE") "namespace" else sub("/.*", "", path),
    proposed_package = cinfo$category,
    confidence = cinfo$confidence,
    lines = line_count(path),
    exported_functions_detected = asset_exports(path),
    rationale = cinfo$rationale,
    stringsAsFactors = FALSE
  )
})

inventory <- do.call(rbind, rows)

out_dir <- file.path(repo_root, "split_audit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(inventory, file.path(out_dir, "gflow_split_inventory.csv"), row.names = FALSE)

summary <- aggregate(asset ~ proposed_package, inventory, length)
names(summary)[2] <- "asset_count"
summary <- summary[order(summary$proposed_package), ]
write.csv(summary, file.path(out_dir, "gflow_split_inventory_summary.csv"), row.names = FALSE)

cat("Wrote ", nrow(inventory), " inventory rows\n", sep = "")
print(summary, row.names = FALSE)

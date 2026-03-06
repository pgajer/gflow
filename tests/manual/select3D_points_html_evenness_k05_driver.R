#!/usr/bin/env Rscript

`%or%` <- function(x, y) if (is.null(x)) y else x

timestamp.now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

parse.bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(isTRUE(default))
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  isTRUE(default)
}

auto.sphere.radius <- function(X) {
  rng <- apply(X, 2, function(v) diff(range(v, na.rm = TRUE)))
  max(1e-08, 0.01 * mean(rng))
}

detect.repo.root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  script.file <- NULL
  if (length(file.arg) > 0L) {
    script.file <- normalizePath(sub("^--file=", "", file.arg[1L]), mustWork = FALSE)
  }
  if (is.null(script.file)) {
    frm <- try(sys.frames()[[1L]], silent = TRUE)
    if (!inherits(frm, "try-error") && !is.null(frm$ofile)) {
      script.file <- normalizePath(frm$ofile, mustWork = FALSE)
    }
  }
  if (!is.null(script.file) && nzchar(script.file)) {
    return(normalizePath(file.path(dirname(script.file), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

ensure.gflow <- function(repo.root) {
  has.gflow <- requireNamespace("gflow", quietly = TRUE)
  if (has.gflow) {
    if ("select3D.points.html" %in% getNamespaceExports("gflow")) return(invisible(TRUE))
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop(
      "Need either an updated installed gflow package or pkgload to load local sources.\n",
      "Install pkgload or reinstall gflow from this repo."
    )
  }
  pkgload::load_all(path = repo.root, quiet = TRUE)
  if (!("select3D.points.html" %in% getNamespaceExports("gflow"))) {
    stop("Local gflow load succeeded, but select3D.points.html() is still unavailable.")
  }
  invisible(TRUE)
}

parse.args <- function(args) {
  out <- list(
    bundle.file = "/Users/pgajer/current_projects/symptoms/results/asv_full_graph_evenness_endpoints_k05/evenness.endpoints.k05.bundle.rds",
    run.meta.file = "/Users/pgajer/current_projects/symptoms/results/asv_full_graph_evenness_endpoints_k05/run.metadata.rds",
    layout.file = NULL,
    out.dir = "/Users/pgajer/current_projects/symptoms/results/asv_full_graph_evenness_endpoints_k05/manual_selection",
    selector = "plotly",
    hide.gray = TRUE,
    launch.browser = TRUE
  )

  if (length(args) < 1L) return(out)

  for (aa in args) {
    if (grepl("^--bundle-file=", aa)) out$bundle.file <- sub("^--bundle-file=", "", aa)
    if (grepl("^--run-meta-file=", aa)) out$run.meta.file <- sub("^--run-meta-file=", "", aa)
    if (grepl("^--layout-file=", aa)) out$layout.file <- sub("^--layout-file=", "", aa)
    if (grepl("^--out-dir=", aa)) out$out.dir <- sub("^--out-dir=", "", aa)
    if (grepl("^--selector=", aa)) out$selector <- tolower(trimws(sub("^--selector=", "", aa)))
    if (grepl("^--hide-gray=", aa)) out$hide.gray <- parse.bool(sub("^--hide-gray=", "", aa), default = out$hide.gray)
    if (grepl("^--launch-browser=", aa)) out$launch.browser <- parse.bool(sub("^--launch-browser=", "", aa), default = out$launch.browser)
  }

  if (!(out$selector %in% c("plotly", "rglhtml"))) {
    stop("--selector must be one of: plotly, rglhtml")
  }

  out
}

main <- function() {
  repo.root <- detect.repo.root()
  ensure.gflow(repo.root)

  opts <- parse.args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(opts$bundle.file)) stop("Bundle file not found: ", opts$bundle.file)
  bundle <- readRDS(opts$bundle.file)

  if (is.null(opts$layout.file) && file.exists(opts$run.meta.file)) {
    meta <- readRDS(opts$run.meta.file)
    opts$layout.file <- meta$opts$layout.file %or% opts$layout.file
  }
  if (is.null(opts$layout.file) || !file.exists(opts$layout.file)) {
    stop("Layout file not found. Provide --layout-file=... or a valid --run-meta-file.")
  }

  layout.3d <- as.matrix(readRDS(opts$layout.file))
  if (!is.matrix(layout.3d) || ncol(layout.3d) != 3L) {
    stop("Layout must be an n x 3 matrix: ", opts$layout.file)
  }

  y.full <- rep(NA_real_, nrow(layout.3d))
  lcc.idx <- as.integer(bundle$lcc$lcc.index.global %or% integer(0))
  yhat <- as.numeric(bundle$fitted.evenness.lcc %or% numeric(0))
  n.assign <- min(length(lcc.idx), length(yhat))
  if (n.assign > 0L) {
    y.full[lcc.idx[seq_len(n.assign)]] <- yhat[seq_len(n.assign)]
  }
  subset.fit <- is.finite(y.full)

  end.labels.global <- setNames(
    as.character(bundle$end.labels %or% character(0)),
    as.character(bundle$end.vertices.global %or% integer(0))
  )

  plot.X <- layout.3d
  plot.y <- y.full
  plot.end.labels <- end.labels.global
  plot.args <- list(
    subset = subset.fit,
    non.highlight.type = "point",
    highlight.type = "sphere",
    non.highlight.color = "gray78",
    point.size = 2.5,
    radius = auto.sphere.radius(layout.3d) * 0.8,
    quantize.method = "quantile",
    n.levels = 11,
    legend.title = "Smoothed evenness (k=5)",
    end.labels = plot.end.labels,
    end.labels.col = "black",
    end.labels.cex = 1.6,
    end.labels.radius = 0.2,
    background.color = "white"
  )

  if (isTRUE(opts$hide.gray)) {
    keep <- which(subset.fit)
    if (length(keep) < 1L) stop("No finite fitted values in y.full; cannot hide gray component.")
    plot.X <- layout.3d[keep, , drop = FALSE]
    plot.y <- y.full[keep]
    loc <- match(as.integer(names(end.labels.global)), keep)
    ok <- is.finite(loc)
    plot.end.labels <- setNames(as.character(end.labels.global[ok]), as.character(loc[ok]))
    plot.args <- list(
      subset = NULL,
      highlight.type = "sphere",
      point.size = 2.5,
      radius = auto.sphere.radius(plot.X) * 0.8,
      quantize.method = "quantile",
      n.levels = 11,
      legend.title = "Smoothed evenness (k=5, gray hidden)",
      end.labels = plot.end.labels,
      end.labels.col = "black",
      end.labels.cex = 1.6,
      end.labels.radius = 0.2,
      background.color = "white"
    )
  }

  cat(sprintf(
    "Launching %s selector on %d points (hide.gray=%s)\n",
    opts$selector, nrow(plot.X), opts$hide.gray
  ))

  if (identical(opts$selector, "rglhtml")) {
    sel <- gflow::select3D.points.html(
      X = plot.X,
      plot.type = "cont",
      y = plot.y,
      plot.args = plot.args,
      title = "Manual endpoint selection for symptoms k=5 graph",
      widget.width = 1700L,
      widget.height = 1000L,
      launch.browser = isTRUE(opts$launch.browser)
    )
  } else {
    sel <- gflow::select3D.points.plotly(
      X = plot.X,
      plot.type = "cont",
      y = plot.y,
      subset = NULL,
      point.size = 3,
      title = "Manual endpoint selection for symptoms k=5 graph (Plotly)",
      end.labels = plot.end.labels,
      launch.browser = isTRUE(opts$launch.browser),
      save.dir = opts$out.dir,
      save.prefix = "manual_selected_local",
      save.on.done = TRUE,
      allow.edit.save.dir = TRUE
    )
  }

  sel.local <- as.integer(sel$idx %or% integer(0))
  if (isTRUE(opts$hide.gray)) {
    keep <- which(subset.fit)
    sel.global <- keep[sel.local]
  } else {
    sel.global <- sel.local
  }
  sel.global <- as.integer(sel.global)
  sel.global <- sel.global[is.finite(sel.global) & sel.global >= 1L & sel.global <= nrow(layout.3d)]
  sel.global <- sort(unique(sel.global))

  dir.create(opts$out.dir, recursive = TRUE, showWarnings = FALSE)
  tag <- timestamp.now()

  out.tbl <- data.frame(
    vertex.global = sel.global,
    label = as.character(end.labels.global[as.character(sel.global)]),
    stringsAsFactors = FALSE
  )
  out.tbl$label[is.na(out.tbl$label)] <- ""

  csv.file <- file.path(opts$out.dir, paste0("manual_selected_vertices_", tag, ".csv"))
  rds.file <- file.path(opts$out.dir, paste0("manual_selected_vertices_", tag, ".rds"))

  write.csv(out.tbl, csv.file, row.names = FALSE)
  saveRDS(list(
    selected = out.tbl,
    raw.selection = sel,
    options = opts,
    bundle.file = normalizePath(opts$bundle.file, mustWork = FALSE),
    layout.file = normalizePath(opts$layout.file, mustWork = FALSE),
    generated.at = as.character(Sys.time())
  ), rds.file)

  cat(sprintf("Selected %d global vertices.\n", nrow(out.tbl)))
  cat(sprintf("CSV: %s\n", csv.file))
  cat(sprintf("RDS: %s\n", rds.file))
}

if (sys.nframe() == 0L) {
  main()
}

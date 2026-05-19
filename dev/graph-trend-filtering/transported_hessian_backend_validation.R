`%||%` <- function(a, b) if (is.null(a)) b else a

transported.hessian.load.packages <- function(repo.dir,
                                              use.local.grip = TRUE,
                                              local.grip.dir = "/Users/pgajer/current_projects/grip") {
  if (isTRUE(use.local.grip) &&
      dir.exists(local.grip.dir) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(local.grip.dir, quiet = TRUE)
  }
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo.dir, quiet = TRUE)
  } else {
    library(gflow)
  }
}

transported.hessian.edge.kk.available <- function() {
  if (!requireNamespace("grip", quietly = TRUE)) return(FALSE)
  ns <- asNamespace("grip")
  exists("grip.optimize.edge.kk.layout", envir = ns, inherits = FALSE) ||
    exists("grip.optimize.edge.isometric.layout", envir = ns, inherits = FALSE)
}

transported.hessian.edge.kk.function.name <- function() {
  if (!requireNamespace("grip", quietly = TRUE)) return(NA_character_)
  ns <- asNamespace("grip")
  if (exists("grip.optimize.edge.kk.layout", envir = ns, inherits = FALSE)) {
    return("grip.optimize.edge.kk.layout")
  }
  if (exists("grip.optimize.edge.isometric.layout", envir = ns, inherits = FALSE)) {
    return("grip.optimize.edge.isometric.layout")
  }
  NA_character_
}

transported.hessian.path.graph <- function(n = 12L) {
  x <- seq(0, 1, length.out = n)
  adj <- vector("list", n)
  w <- vector("list", n)
  for (i in seq_len(n)) {
    nbrs <- integer()
    lens <- numeric()
    if (i > 1L) {
      nbrs <- c(nbrs, i - 1L)
      lens <- c(lens, x[i] - x[i - 1L])
    }
    if (i < n) {
      nbrs <- c(nbrs, i + 1L)
      lens <- c(lens, x[i + 1L] - x[i])
    }
    adj[[i]] <- nbrs
    w[[i]] <- lens
  }
  list(
    case.id = "path_1d",
    graph.type = "path",
    expected.dimension = 1L,
    adj.list = adj,
    weight.list = w,
    coordinates = matrix(x, ncol = 1)
  )
}

transported.hessian.grid.graph <- function(nx = 4L, ny = 4L, jitter = 0,
                                           seed = 1L, case.id = "grid_2d") {
  set.seed(seed)
  coords <- expand.grid(x = seq(0, 1, length.out = nx),
                        y = seq(0, 1, length.out = ny))
  coords <- as.matrix(coords)
  if (jitter > 0) {
    coords <- coords + matrix(stats::rnorm(length(coords), sd = jitter),
                              nrow = nrow(coords))
  }
  index <- function(ix, iy) (iy - 1L) * nx + ix
  n <- nx * ny
  adj <- vector("list", n)
  w <- vector("list", n)
  for (iy in seq_len(ny)) {
    for (ix in seq_len(nx)) {
      i <- index(ix, iy)
      nbrs <- integer()
      if (ix > 1L) nbrs <- c(nbrs, index(ix - 1L, iy))
      if (ix < nx) nbrs <- c(nbrs, index(ix + 1L, iy))
      if (iy > 1L) nbrs <- c(nbrs, index(ix, iy - 1L))
      if (iy < ny) nbrs <- c(nbrs, index(ix, iy + 1L))
      adj[[i]] <- nbrs
      w[[i]] <- sqrt(rowSums((coords[nbrs, , drop = FALSE] -
                              matrix(coords[i, ], nrow = length(nbrs),
                                     ncol = 2L, byrow = TRUE))^2))
    }
  }
  list(
    case.id = case.id,
    graph.type = if (jitter > 0) "jittered_grid" else "grid",
    expected.dimension = 2L,
    adj.list = adj,
    weight.list = w,
    coordinates = coords
  )
}

transported.hessian.sknn.graph <- function(n = 24L, k = 5L, seed = 4L,
                                           case.id = sprintf("sampled_sknn_2d_n%d", n)) {
  set.seed(seed)
  coords <- cbind(stats::runif(n), stats::runif(n))
  graph <- create.sknn.graph(coords, k = k, connect.components = TRUE)
  list(
    case.id = case.id,
    graph.type = "connected_sknn",
    expected.dimension = 2L,
    adj.list = graph$adj_list,
    weight.list = graph$weight_list,
    coordinates = coords,
    graph.metadata = list(
      k = k,
      n.components.before = graph$n_components_before,
      n.components.after = graph$n_components_after,
      n.mst.edges.added = graph$n_mst_edges_added
    )
  )
}

transported.hessian.polynomial.probes <- function(coords) {
  coords <- as.matrix(coords)
  if (ncol(coords) == 1L) {
    x <- coords[, 1L]
    return(cbind(
      constant = rep(1, length(x)),
      x = x,
      x2 = x^2
    ))
  }
  x <- coords[, 1L]
  y <- coords[, 2L]
  cbind(
    constant = rep(1, length(x)),
    x = x,
    y = y,
    x2 = x^2,
    xy = x * y,
    y2 = y^2
  )
}

transported.hessian.residual.summary <- function(op) {
  tab <- op$diagnostics$polynomial.residuals$per.column
  named <- stats::setNames(tab$residual, tab$probe)
  linear.names <- intersect(c("x", "y", "z"), names(named))
  quadratic.names <- intersect(c("x2", "xy", "y2", "xz", "yz", "z2"), names(named))
  data.frame(
    residual.constant = unname(named["constant"] %||% NA_real_),
    residual.linear.mean = if (length(linear.names)) {
      mean(named[linear.names], na.rm = TRUE)
    } else {
      NA_real_
    },
    residual.quadratic.mean = if (length(quadratic.names)) {
      mean(named[quadratic.names], na.rm = TRUE)
    } else {
      NA_real_
    },
    residual.overall = op$diagnostics$polynomial.residuals$overall %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.summarize.operator <- function(op, case, method,
                                                   elapsed.sec, status,
                                                   filter.name = "none",
                                                   message = NA_character_) {
  emb <- op$transport$embedding.table %||% data.frame()
  row.tab <- op$row.table %||% data.frame()
  unique.rows <- if (nrow(row.tab) && "row" %in% names(row.tab)) {
    row.tab[!duplicated(row.tab$row), , drop = FALSE]
  } else {
    data.frame()
  }
  residuals <- transported.hessian.residual.summary(op)
  cbind(
    data.frame(
      case.id = case$case.id,
      graph.type = case$graph.type,
      n.vertices = length(case$adj.list),
      expected.dimension = case$expected.dimension,
      local.embedding.method = method,
      transport.filter = filter.name,
      status = status,
      runtime.sec = elapsed.sec,
      nullity.estimate = op$diagnostics$nullity.estimate,
      n.candidate.rows = op$diagnostics$summary$n.candidate.rows,
      n.rows = op$diagnostics$summary$n.rows,
      n.row.records = op$diagnostics$summary$n.row.records,
      n.dropped = op$diagnostics$summary$n.dropped,
      row.retention = if (op$diagnostics$summary$n.candidate.rows > 0) {
        op$diagnostics$summary$n.rows / op$diagnostics$summary$n.candidate.rows
      } else {
        NA_real_
      },
      median.edge.stress = if (nrow(emb)) stats::median(emb$edge.stress, na.rm = TRUE) else NA_real_,
      median.graph.distance.stress = if (nrow(emb)) stats::median(emb$graph.distance.stress, na.rm = TRUE) else NA_real_,
      mean.transport.entropy = if (nrow(unique.rows)) mean(unique.rows$transport.entropy, na.rm = TRUE) else NA_real_,
      mean.effective.matches = if (nrow(unique.rows)) mean(unique.rows$effective.matches, na.rm = TRUE) else NA_real_,
      mean.effective.match.fraction = if (nrow(unique.rows) &&
                                           "effective.match.fraction" %in% names(unique.rows)) {
        mean(unique.rows$effective.match.fraction, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.match.score = if (nrow(row.tab)) mean(row.tab$match.score, na.rm = TRUE) else NA_real_,
      mean.best.score.robust.z = if (nrow(unique.rows) &&
                                      "best.score.robust.z" %in% names(unique.rows)) {
        mean(unique.rows$best.score.robust.z, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.margin.robust.z = if (nrow(unique.rows) &&
                                  "margin.robust.z" %in% names(unique.rows)) {
        mean(unique.rows$margin.robust.z, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.best.score.local.quantile = if (nrow(unique.rows) &&
                                            "best.score.local.quantile" %in% names(unique.rows)) {
        mean(unique.rows$best.score.local.quantile, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.margin.local.quantile = if (nrow(unique.rows) &&
                                        "margin.local.quantile" %in% names(unique.rows)) {
        mean(unique.rows$margin.local.quantile, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.best.match.margin = if (nrow(unique.rows) && "best.match.margin" %in% names(unique.rows)) {
        mean(unique.rows$best.match.margin, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.best.match.angle = if (nrow(unique.rows) && "best.match.angle" %in% names(unique.rows)) {
        mean(unique.rows$best.match.angle, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean.best.length.relative.error = if (nrow(unique.rows) &&
                                            "best.length.relative.error" %in% names(unique.rows)) {
        mean(unique.rows$best.length.relative.error, na.rm = TRUE)
      } else {
        NA_real_
      },
      backend.used = if (nrow(emb)) paste(sort(unique(emb$backend.used)), collapse = "; ") else NA_character_,
      fallback.reason = if (nrow(emb)) {
        vals <- unique(emb$fallback.reason[!is.na(emb$fallback.reason)])
        if (length(vals)) paste(vals, collapse = "; ") else NA_character_
      } else {
        NA_character_
      },
      message = message,
      stringsAsFactors = FALSE
    ),
    residuals
  )
}

transported.hessian.run.one <- function(case, method, local.disk.hops = 1L,
                                        local.max.vertices = 18L,
                                        soft.bandwidth = 0.18,
                                        filter.name = "none",
                                        filter.args = list()) {
  probes <- transported.hessian.polynomial.probes(case$coordinates)
  start <- proc.time()[["elapsed"]]
  op.args <- c(
    list(
      case$adj.list,
      case$weight.list,
      transport.rule = "local.embedding.soft",
      polynomial.probes = probes,
      local.embedding.method = method,
      local.embedding.dim = max(2L, case$expected.dimension),
      local.disk.hops = local.disk.hops,
      local.max.vertices = local.max.vertices,
      soft.bandwidth = soft.bandwidth
    ),
    filter.args
  )
  fit <- try(
    do.call(transported.graph.hessian.operator, op.args),
    silent = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(fit, "try-error")) {
    summary <- data.frame(
      case.id = case$case.id,
      graph.type = case$graph.type,
      n.vertices = length(case$adj.list),
      expected.dimension = case$expected.dimension,
      local.embedding.method = method,
      transport.filter = filter.name,
      status = "error",
      runtime.sec = elapsed,
      nullity.estimate = NA_integer_,
      n.candidate.rows = NA_integer_,
      n.rows = NA_integer_,
      n.row.records = NA_integer_,
      n.dropped = NA_integer_,
      row.retention = NA_real_,
      median.edge.stress = NA_real_,
      median.graph.distance.stress = NA_real_,
      mean.transport.entropy = NA_real_,
      mean.effective.matches = NA_real_,
      mean.effective.match.fraction = NA_real_,
      mean.match.score = NA_real_,
      mean.best.score.robust.z = NA_real_,
      mean.margin.robust.z = NA_real_,
      mean.best.score.local.quantile = NA_real_,
      mean.margin.local.quantile = NA_real_,
      mean.best.match.margin = NA_real_,
      mean.best.match.angle = NA_real_,
      mean.best.length.relative.error = NA_real_,
      backend.used = NA_character_,
      fallback.reason = NA_character_,
      message = as.character(fit),
      residual.constant = NA_real_,
      residual.linear.mean = NA_real_,
      residual.quadratic.mean = NA_real_,
      residual.overall = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, embedding = data.frame(), rows = data.frame(),
                dropped = data.frame()))
  }
  summary <- transported.hessian.summarize.operator(
    fit, case, method, elapsed, "ok",
    filter.name = filter.name
  )
  embedding <- fit$transport$embedding.table %||% data.frame()
  if (nrow(embedding)) {
    embedding$case.id <- case$case.id
    embedding$graph.type <- case$graph.type
    embedding$local.embedding.method <- method
    embedding$transport.filter <- filter.name
  }
  rows <- fit$row.table %||% data.frame()
  if (nrow(rows)) {
    rows <- rows[!duplicated(rows$row), , drop = FALSE]
    rows$case.id <- case$case.id
    rows$graph.type <- case$graph.type
    rows$local.embedding.method <- method
    rows$transport.filter <- filter.name
  }
  dropped <- fit$dropped.table %||% data.frame()
  if (nrow(dropped)) {
    dropped$case.id <- case$case.id
    dropped$graph.type <- case$graph.type
    dropped$local.embedding.method <- method
    dropped$transport.filter <- filter.name
  }
  list(summary = summary, embedding = embedding, rows = rows, dropped = dropped)
}

transported.hessian.backend.validation.cases <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  cases <- list(
    transported.hessian.path.graph(12L),
    transported.hessian.grid.graph(4L, 4L, jitter = 0, case.id = "grid_4x4"),
    transported.hessian.grid.graph(4L, 4L, jitter = 0.025, seed = 3L,
                                   case.id = "jittered_grid_4x4"),
    transported.hessian.sknn.graph(24L, k = 5L, seed = 5L,
                                   case.id = "sampled_sknn_2d_n24")
  )
  if (identical(mode, "full")) {
    cases <- c(cases, list(
      transported.hessian.grid.graph(5L, 5L, jitter = 0, case.id = "grid_5x5"),
      transported.hessian.grid.graph(5L, 5L, jitter = 0.02, seed = 7L,
                                     case.id = "jittered_grid_5x5"),
      transported.hessian.sknn.graph(32L, k = 6L, seed = 11L,
                                     case.id = "sampled_sknn_2d_n32")
    ))
  }
  cases
}

run.transported.hessian.backend.validation <- function(mode = c("smoke", "full"),
                                                       repo.dir = "/Users/pgajer/current_projects/gflow",
                                                       output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
                                                       use.local.grip = TRUE) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.load.packages(repo.dir, use.local.grip = use.local.grip)

  methods <- c("cmdscale", "mds.edge.kk", "grip.edge.kk")
  cases <- transported.hessian.backend.validation.cases(mode)
  runs <- list()
  counter <- 0L
  for (case in cases) {
    for (method in methods) {
      counter <- counter + 1L
      runs[[counter]] <- transported.hessian.run.one(case, method)
    }
  }

  summary <- do.call(rbind, lapply(runs, `[[`, "summary"))
  embedding <- do.call(rbind, lapply(runs, `[[`, "embedding"))
  rows <- do.call(rbind, lapply(runs, `[[`, "rows"))

  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    use.local.grip = use.local.grip,
    grip.available = requireNamespace("grip", quietly = TRUE),
    edge.kk.available = transported.hessian.edge.kk.available(),
    edge.kk.function = transported.hessian.edge.kk.function.name(),
    gflow.repo = repo.dir
  )
  out <- list(
    metadata = metadata,
    summary = summary,
    embedding = embedding,
    rows = rows
  )

  prefix <- file.path(output.dir, paste0("transported_hessian_backend_validation_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(embedding, paste0(prefix, "_embedding.csv"), row.names = FALSE)
  utils::write.csv(rows, paste0(prefix, "_rows.csv"), row.names = FALSE)
  out
}

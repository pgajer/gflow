repo.dir <- "/Users/pgajer/current_projects/gflow"

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
} else {
  library(gflow)
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required.")
}

add_undirected_edge <- function(adj, weights, i, j, weight = 1) {
  adj[[i]] <- c(adj[[i]], j)
  weights[[i]] <- c(weights[[i]], weight)
  adj[[j]] <- c(adj[[j]], i)
  weights[[j]] <- c(weights[[j]], weight)
  list(adj = adj, weights = weights)
}

make_graph_from_edges <- function(n, edges, weights = rep(1, nrow(edges))) {
  adj <- vector("list", n)
  w <- vector("list", n)
  for (idx in seq_len(nrow(edges))) {
    out <- add_undirected_edge(adj, w, edges[idx, 1], edges[idx, 2], weights[idx])
    adj <- out$adj
    w <- out$weights
  }
  list(adj.list = lapply(adj, as.integer), weight.list = lapply(w, as.double))
}

make_path_graph <- function(n, weights = rep(1, n - 1L)) {
  make_graph_from_edges(n, cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L), weights)
}

make_three_arm_star <- function(arm.length = 2L, weights = NULL) {
  n <- 1L + 3L * arm.length
  edges <- matrix(NA_integer_, nrow = 3L * arm.length, ncol = 2L)
  row <- 1L
  for (arm in 0:2) {
    first <- 2L + arm * arm.length
    edges[row, ] <- c(1L, first)
    row <- row + 1L
    if (arm.length > 1L) {
      for (step in seq_len(arm.length - 1L)) {
        edges[row, ] <- c(first + step - 1L, first + step)
        row <- row + 1L
      }
    }
  }
  if (is.null(weights)) weights <- rep(1, nrow(edges))
  make_graph_from_edges(n, edges, weights)
}

make_cycle_graph <- function(n) {
  edges <- cbind(seq_len(n), c(seq_len(n - 1L) + 1L, 1L))
  make_graph_from_edges(n, edges, rep(1, n))
}

make_grid_graph <- function(nrow = 3L, ncol = 3L) {
  idx <- function(i, j) (i - 1L) * ncol + j
  edges <- matrix(integer(), ncol = 2L)
  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      if (j < ncol) edges <- rbind(edges, c(idx(i, j), idx(i, j + 1L)))
      if (i < nrow) edges <- rbind(edges, c(idx(i, j), idx(i + 1L, j)))
    }
  }
  make_graph_from_edges(nrow * ncol, edges, rep(1, nrow(edges)))
}

make_symmetric_knn_1d_graph <- function(x, k = 3L) {
  n <- length(x)
  dist.mat <- abs(outer(x, x, "-"))
  edge.map <- new.env(parent = emptyenv())
  for (i in seq_len(n)) {
    ord <- order(dist.mat[i, ], seq_len(n))
    nbrs <- setdiff(ord, i)[seq_len(k)]
    for (j in nbrs) {
      key <- paste(sort(c(i, j)), collapse = ":")
      edge.map[[key]] <- TRUE
    }
  }
  keys <- ls(edge.map)
  edges <- do.call(rbind, strsplit(keys, ":", fixed = TRUE))
  edges <- matrix(as.integer(edges), ncol = 2L)
  weights <- abs(x[edges[, 1]] - x[edges[, 2]])
  make_graph_from_edges(n, edges, weights)
}

path_vertices <- function(path.string) {
  as.integer(strsplit(path.string, "-", fixed = TRUE)[[1]])
}

path_edges <- function(vertices) {
  if (length(vertices) < 2L) return(character())
  vapply(seq_len(length(vertices) - 1L), function(i) {
    paste(sort(vertices[i + 0:1]), collapse = ":")
  }, character(1))
}

path_overlap_matrix <- function(path.list, metric = c("vertices", "edges")) {
  metric <- match.arg(metric)
  sets <- lapply(path.list, function(v) {
    if (metric == "vertices") as.character(v) else path_edges(v)
  })
  n <- length(sets)
  mat <- matrix(0L, nrow = n, ncol = n)
  if (n < 1L) return(mat)
  for (i in seq_len(n)) {
    for (j in i:n) {
      val <- length(intersect(sets[[i]], sets[[j]]))
      mat[i, j] <- val
      mat[j, i] <- val
    }
  }
  diag(mat) <- 0L
  mat
}

coverage_counts <- function(path.list, n.vertices) {
  vertex.count <- integer(n.vertices)
  edge.env <- new.env(parent = emptyenv())
  for (path in path.list) {
    vertex.count[path] <- vertex.count[path] + 1L
    for (edge in path_edges(path)) edge.env[[edge]] <- (edge.env[[edge]] %||% 0L) + 1L
  }
  list(vertex = vertex.count, edge = unlist(as.list(edge.env), use.names = TRUE))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

select_overlap_thinned <- function(path.list, max.vertex.overlap = 1L,
                                   min.vertex.coverage = 1L) {
  n <- length(path.list)
  selected <- integer()
  if (n < 1L) return(selected)
  overlap <- path_overlap_matrix(path.list, "vertices")
  path.degree <- rowSums(overlap > 0)
  order.idx <- order(path.degree, vapply(path.list, `[`, integer(1), 1L), seq_len(n))
  for (idx in order.idx) {
    if (length(selected) < 1L || all(overlap[idx, selected] <= max.vertex.overlap)) {
      selected <- c(selected, idx)
    }
  }
  n.vertices <- max(unlist(path.list, use.names = FALSE))
  covered <- coverage_counts(path.list[selected], n.vertices)$vertex
  while (any(covered < min.vertex.coverage)) {
    candidates <- setdiff(seq_len(n), selected)
    if (length(candidates) < 1L) break
    gains <- vapply(candidates, function(idx) {
      sum(unique(path.list[[idx]])[covered[unique(path.list[[idx]])] < min.vertex.coverage] > 0L)
    }, integer(1))
    if (max(gains) <= 0L) break
    selected <- c(selected, candidates[which.max(gains)])
    covered <- coverage_counts(path.list[selected], n.vertices)$vertex
  }
  sort(selected)
}

select_coverage_greedy <- function(path.list, min.vertex.coverage = 1L,
                                   min.edge.coverage = 1L) {
  n <- length(path.list)
  selected <- integer()
  if (n < 1L) return(selected)
  n.vertices <- max(unlist(path.list, use.names = FALSE))
  all.edges <- sort(unique(unlist(lapply(path.list, path_edges), use.names = FALSE)))
  vertex.cover <- integer(n.vertices)
  edge.cover <- setNames(integer(length(all.edges)), all.edges)
  while (any(vertex.cover < min.vertex.coverage) ||
         any(edge.cover < min.edge.coverage)) {
    candidates <- setdiff(seq_len(n), selected)
    if (length(candidates) < 1L) break
    gains <- vapply(candidates, function(idx) {
      verts <- unique(path.list[[idx]])
      edges <- path_edges(path.list[[idx]])
      sum(vertex.cover[verts] < min.vertex.coverage) +
        sum(edge.cover[edges] < min.edge.coverage)
    }, integer(1))
    if (max(gains) <= 0L) break
    selected <- c(selected, candidates[which.max(gains)])
    path <- path.list[[tail(selected, 1L)]]
    vertex.cover[unique(path)] <- vertex.cover[unique(path)] + 1L
    edges <- path_edges(path)
    edge.cover[edges] <- edge.cover[edges] + 1L
  }
  sort(selected)
}

operator_from_selected_rows <- function(base.operator, selected.rows) {
  mat <- base.operator$penalty$matrix[selected.rows, , drop = FALSE]
  base.operator$penalty$matrix <- mat
  base.operator$penalty$i <- NULL
  base.operator$penalty$j <- NULL
  base.operator$penalty$x <- NULL
  trip <- methods::as(mat, "TsparseMatrix")
  base.operator$penalty$i <- as.integer(trip@i + 1L)
  base.operator$penalty$j <- as.integer(trip@j + 1L)
  base.operator$penalty$x <- as.double(trip@x)
  base.operator$penalty$dim <- dim(mat)
  base.operator$path.table <- base.operator$path.table[selected.rows, , drop = FALSE]
  base.operator$path.table$row <- seq_len(nrow(base.operator$path.table))
  base.operator$nullity.estimate <- if (ncol(mat) <= 200L) {
    as.integer(ncol(mat) - Matrix::rankMatrix(as.matrix(mat))[1])
  } else {
    NA_integer_
  }
  base.operator
}

diagnose_path_selection <- function(graph.name, graph, order = 1L,
                                    selection = c("all", "branch.continuation",
                                                  "overlap.thinned",
                                                  "coverage.greedy")) {
  selection <- match.arg(selection)
  path.family <- if (selection == "branch.continuation") "branch.continuation" else "all.simple"
  op <- tryCatch(
    graph.trend.filtering.operator(
      graph$adj.list,
      graph$weight.list,
      order = order,
      operator.family = "path.divided.difference",
      path.family = path.family,
      path.weighting = "unit",
      return.sparse = TRUE
    ),
    error = function(e) e
  )
  if (inherits(op, "error")) {
    n.vertices <- length(graph$adj.list)
    n.edges <- sum(lengths(graph$adj.list)) / 2
    empty.summary <- data.frame(
      graph = graph.name,
      selection = selection,
      n.vertices = n.vertices,
      n.edges = n.edges,
      order = order,
      difference.order = order + 1L,
      candidate.rows = 0L,
      selected.rows = 0L,
      nullity = NA_integer_,
      min.vertex.coverage = 0L,
      median.vertex.coverage = 0,
      max.vertex.coverage = 0L,
      mean.path.overlap.degree = 0,
      max.path.overlap.degree = 0,
      status = "error",
      message = conditionMessage(op),
      stringsAsFactors = FALSE
    )
    return(list(
      graph.name = graph.name,
      selection = selection,
      order = order,
      difference.order = order + 1L,
      operator = NULL,
      path.list = list(),
      overlap = matrix(0L, nrow = 0L, ncol = 0L),
      summary = empty.summary,
      vertex.coverage = data.frame(
        graph = graph.name,
        selection = selection,
        vertex = seq_len(n.vertices),
        coverage = integer(n.vertices)
      ),
      path.table = data.frame(
        row = integer(),
        start = integer(),
        end = integer(),
        vertices = character(),
        metric.length = numeric(),
        row.weight = numeric(),
        path.family = character(),
        path.weighting = character(),
        duplicates.removed = integer(),
        graph = character(),
        selection = character()
      )
    ))
  }
  paths <- lapply(op$path.table$vertices, path_vertices)
  selected <- seq_along(paths)
  if (selection == "overlap.thinned") {
    selected <- select_overlap_thinned(paths)
    op <- operator_from_selected_rows(op, selected)
    paths <- paths[selected]
  } else if (selection == "coverage.greedy") {
    selected <- select_coverage_greedy(paths)
    op <- operator_from_selected_rows(op, selected)
    paths <- paths[selected]
  }

  overlap <- path_overlap_matrix(paths, "vertices")
  cov <- coverage_counts(paths, op$graph$n.vertices)
  row.count <- nrow(op$penalty$matrix)
  total.rows <- if (selection %in% c("all", "branch.continuation")) row.count else nrow(
    graph.trend.filtering.operator(
      graph$adj.list, graph$weight.list,
      order = order,
      operator.family = "path.divided.difference",
      path.family = "all.simple",
      path.weighting = "unit"
    )$penalty$matrix
  )
  out <- list(
    graph.name = graph.name,
    selection = selection,
    order = order,
    difference.order = order + 1L,
    operator = op,
    path.list = paths,
    overlap = overlap,
    summary = data.frame(
      graph = graph.name,
      selection = selection,
      n.vertices = op$graph$n.vertices,
      n.edges = op$graph$n.edges,
      order = order,
      difference.order = order + 1L,
      candidate.rows = total.rows,
      selected.rows = row.count,
      nullity = op$nullity.estimate,
      min.vertex.coverage = min(cov$vertex),
      median.vertex.coverage = stats::median(cov$vertex),
      max.vertex.coverage = max(cov$vertex),
      mean.path.overlap.degree = if (row.count > 1L) mean(rowSums(overlap > 0)) else 0,
      max.path.overlap.degree = if (row.count > 1L) max(rowSums(overlap > 0)) else 0,
      status = "ok",
      message = "",
      stringsAsFactors = FALSE
    ),
    vertex.coverage = data.frame(
      graph = graph.name,
      selection = selection,
      vertex = seq_along(cov$vertex),
      coverage = cov$vertex
    ),
    path.table = transform(op$path.table, graph = graph.name, selection = selection)
  )
  out
}

build_phase4_path_selection_diagnostics <- function() {
  set.seed(20260516)
  graphs <- list(
    path7 = make_path_graph(7),
    star3_arm2 = make_three_arm_star(2),
    star3_arm4 = make_three_arm_star(4),
    cycle8 = make_cycle_graph(8),
    grid3x3 = make_grid_graph(3, 3),
    sknn1d12 = make_symmetric_knn_1d_graph(sort(stats::runif(12)), k = 3)
  )
  selections <- c("all", "branch.continuation", "overlap.thinned", "coverage.greedy")
  diagnostics <- list()
  idx <- 1L
  for (graph.name in names(graphs)) {
    for (selection in selections) {
      diagnostics[[idx]] <- diagnose_path_selection(
        graph.name = graph.name,
        graph = graphs[[graph.name]],
        order = 1L,
        selection = selection
      )
      idx <- idx + 1L
    }
  }
  summaries <- do.call(rbind, lapply(diagnostics, `[[`, "summary"))
  vertex.coverage <- do.call(rbind, lapply(diagnostics, `[[`, "vertex.coverage"))
  path.tables <- do.call(rbind, lapply(diagnostics, `[[`, "path.table"))
  list(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    selections = selections,
    graph.names = names(graphs),
    summaries = summaries,
    vertex.coverage = vertex.coverage,
    path.tables = path.tables,
    diagnostics = diagnostics
  )
}

if (sys.nframe() == 0L) {
  out.dir <- file.path(repo.dir, "dev/graph-trend-filtering/reports")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  out <- build_phase4_path_selection_diagnostics()
  saveRDS(out, file.path(out.dir, "path_selection_diagnostics.rds"))
  utils::write.csv(out$summaries,
                   file.path(out.dir, "path_selection_diagnostics_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(out$vertex.coverage,
                   file.path(out.dir, "path_selection_vertex_coverage.csv"),
                   row.names = FALSE)
  utils::write.csv(out$path.tables,
                   file.path(out.dir, "path_selection_path_table.csv"),
                   row.names = FALSE)
}

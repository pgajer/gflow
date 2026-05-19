repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")

source(file.path(dev.dir, "local-dimension-validation", "compute_local_dimension_validation.R"))

monomial.design.exact.nd <- function(X, degree) {
  X <- as.matrix(X)
  powers <- monomial.powers.nd(ncol(X), degree)
  powers <- powers[rowSums(powers) == degree, , drop = FALSE]
  out <- matrix(1, nrow = nrow(X), ncol = nrow(powers))
  for (col in seq_len(nrow(powers))) {
    for (j in seq_len(ncol(X))) {
      if (powers[col, j] > 0L) {
        out[, col] <- out[, col] * X[, j]^powers[col, j]
      }
    }
  }
  colnames(out) <- apply(powers, 1L, function(z) {
    paste0("x", seq_along(z), "^", z, collapse = " ")
  })
  out
}

make.rank.calibration.case <- function(case.id,
                                       family,
                                       true.dimension,
                                       graph,
                                       intrinsic.coords,
                                       graph.type) {
  list(
    case.id = case.id,
    family = family,
    true.dimension = as.integer(true.dimension),
    graph = graph,
    intrinsic.coords = as.matrix(intrinsic.coords),
    graph.type = graph.type
  )
}

make.rank.calibration.cases <- function(mode = c("smoke", "full"),
                                        seed.base = 20260517L) {
  mode <- match.arg(mode)
  path.graph <- make.path.validation.graph(if (mode == "smoke") 40L else 70L,
                                           closed = FALSE)
  cycle.graph <- make.path.validation.graph(if (mode == "smoke") 48L else 80L,
                                            closed = TRUE)
  grid.graph <- make.grid.validation.graph(
    if (mode == "smoke") 6L else 8L,
    if (mode == "smoke") 6L else 8L,
    "none"
  )
  grid.warp.graph <- make.grid.validation.graph(
    if (mode == "smoke") 6L else 8L,
    if (mode == "smoke") 6L else 8L,
    "nonuniform"
  )
  set.seed(seed.base + 101L)
  U <- sample.square.points(if (mode == "smoke") 60L else 90L,
                            seed = seed.base + 102L)
  square.graph <- make.experimental.graph(U, "delaunay.1skeleton",
                                          metric.coords = U)
  square.graph$param <- U
  X3 <- embed.quadric(U, "paraboloid")
  quad.graph <- make.experimental.graph(U, "delaunay.1skeleton",
                                        metric.coords = X3)
  quad.graph$param <- U
  cube.graph <- make.cube.lattice.validation.graph(
    if (mode == "smoke") 4L else 5L,
    if (mode == "smoke") 4L else 5L,
    if (mode == "smoke") 4L else 5L
  )

  cases <- list(
    make.rank.calibration.case(
      "path", "path", 1L, path.graph,
      intrinsic.coords = path.graph$coords[, 1L, drop = FALSE],
      graph.type = "path"
    ),
    make.rank.calibration.case(
      "cycle", "cycle", 1L, cycle.graph,
      intrinsic.coords = matrix(seq_len(length(cycle.graph$adj.list)),
                                ncol = 1L),
      graph.type = "cycle"
    ),
    make.rank.calibration.case(
      "grid", "grid", 2L, grid.graph,
      intrinsic.coords = grid.graph$param,
      graph.type = "grid"
    ),
    make.rank.calibration.case(
      "grid_nonuniform", "grid", 2L, grid.warp.graph,
      intrinsic.coords = grid.warp.graph$param,
      graph.type = "grid"
    ),
    make.rank.calibration.case(
      "sampled_square_delaunay", "sampled square", 2L, square.graph,
      intrinsic.coords = square.graph$param,
      graph.type = "delaunay.1skeleton"
    ),
    make.rank.calibration.case(
      "paraboloid_delaunay", "quadric surface", 2L, quad.graph,
      intrinsic.coords = quad.graph$param,
      graph.type = "delaunay.1skeleton"
    ),
    make.rank.calibration.case(
      "cube_lattice", "cube lattice", 3L, cube.graph,
      intrinsic.coords = cube.graph$coords,
      graph.type = "cube.lattice"
    )
  )
  if (mode == "smoke") {
    cases[c(1L, 3L, 6L, 7L)]
  } else {
    cases
  }
}

score.rank.calibration.operator <- function(A, coords, degree) {
  intrinsic <- monomial.design.nd(coords, degree)
  control <- monomial.design.exact.nd(coords, degree + 1L)
  data.frame(
    intrinsic.residual = frobenius.residual.ratio(A, intrinsic),
    control.residual = frobenius.residual.ratio(A, control),
    stringsAsFactors = FALSE
  )
}

rbind.fill.base <- function(tables) {
  tables <- tables[vapply(tables, function(x) {
    is.data.frame(x) && nrow(x) > 0L
  }, logical(1))]
  if (length(tables) < 1L) return(data.frame())
  all.names <- unique(unlist(lapply(tables, names), use.names = FALSE))
  normalized <- lapply(tables, function(tab) {
    missing <- setdiff(all.names, names(tab))
    for (nm in missing) tab[[nm]] <- NA
    tab[, all.names, drop = FALSE]
  })
  do.call(rbind, normalized)
}

path.string.to.integer <- function(x) {
  as.integer(strsplit(as.character(x), "-", fixed = TRUE)[[1L]])
}

sparse.from.row.objects <- function(row.objects, n) {
  if (length(row.objects) < 1L) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                dims = c(0L, n)))
  }
  i.idx <- integer()
  j.idx <- integer()
  x.val <- numeric()
  for (row in seq_along(row.objects)) {
    i.idx <- c(i.idx, rep(row, length(row.objects[[row]]$j)))
    j.idx <- c(j.idx, row.objects[[row]]$j)
    x.val <- c(x.val, row.objects[[row]]$x)
  }
  Matrix::sparseMatrix(i = i.idx, j = j.idx, x = x.val,
                       dims = c(length(row.objects), n))
}

candidate.block.for.path <- function(path,
                                     adj.list,
                                     weight.list,
                                     n,
                                     degree,
                                     row.normalization = "unit.l2") {
  path.object <- list(
    vertices = as.integer(path),
    role = "global.rank.connector",
    pass = NA_integer_
  )
  rows <- build.skeleton.row.objects(
    paths = list(path.object),
    adj.list = adj.list,
    weight.list = weight.list,
    n = n,
    degree = degree,
    row.normalization = row.normalization
  )
  sparse.from.row.objects(rows, n)
}

matrix.rank.dense <- function(A, tol = 1e-8) {
  dense <- as.matrix(A)
  if (nrow(dense) < 1L || ncol(dense) < 1L) return(0L)
  sv <- svd(dense, nu = 0, nv = 0)$d
  if (length(sv) < 1L) return(0L)
  as.integer(sum(sv > tol * max(1, max(sv))))
}

candidate.pool.from.local.rank <- function(local.rank.op) {
  candidates <- local.rank.op$skeleton.connector.candidate.table %||% data.frame()
  if (nrow(candidates) < 1L) return(data.frame())
  candidates$path.key <- vapply(
    candidates$path,
    function(z) path.key(path.string.to.integer(z)),
    character(1)
  )
  candidates <- candidates[order(
    -candidates$rank.gain,
    -candidates$coverage.gain,
    -candidates$score,
    candidates$path.length,
    candidates$path.key
  ), , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$path.key), , drop = FALSE]
  row.names(candidates) <- NULL
  candidates
}

candidate.pool.for.row.family <- function(row.family) {
  if (identical(row.family, "global.rank.rich")) {
    "rich.local"
  } else {
    "skeleton.membership"
  }
}

fit.global.rank.calibration.config <- function(case,
                                               cfg,
                                               seed = 1L,
                                               alpha = 0.25,
                                               beta = 0.05,
                                               gamma = 0.01,
                                               eta = 1.00) {
  graph <- case$graph
  n <- length(graph$adj.list)
  expected.nullity <- choose(cfg$degree + case$true.dimension, case$true.dimension)
  base.op <- build.geodesic.annihilator.operator(
    graph$adj.list,
    graph$weight.list,
    degree = cfg$degree,
    path.family = "geodesic.skeleton",
    disk.k = cfg$disk.k,
    max.skeleton.paths = cfg$skeleton.max.paths,
    skeleton.min.coverage = 1L,
    skeleton.max.passes = 2L,
    skeleton.max.cross.paths = 0L,
    row.normalization = "unit.l2",
    vertex.normalization = "mean",
    seed = seed
  )
  local.rank.op <- build.geodesic.annihilator.operator(
    graph$adj.list,
    graph$weight.list,
    degree = cfg$degree,
    path.family = "geodesic.skeleton.cross.enriched",
    disk.k = cfg$disk.k,
    max.skeleton.paths = cfg$skeleton.max.paths,
    skeleton.min.coverage = 1L,
    skeleton.max.passes = 2L,
    skeleton.max.cross.paths = cfg$skeleton.max.cross.paths,
    skeleton.cross.rule = "local.dimension.rank",
    skeleton.cross.embedding.method = "cmdscale",
    skeleton.cross.candidate.pool = cfg$skeleton.cross.candidate.pool,
    skeleton.cross.local.dims = 1:5,
    skeleton.cross.rank.tol = cfg$skeleton.cross.rank.tol,
    skeleton.cross.max.candidates.per.center = 50L,
    skeleton.cross.length.penalty = cfg$skeleton.cross.length.penalty,
    skeleton.cross.coverage.reward = cfg$skeleton.cross.coverage.reward,
    skeleton.cross.iterations = cfg$skeleton.cross.iterations,
    row.normalization = "unit.l2",
    vertex.normalization = "mean",
    seed = seed
  )
  pool <- candidate.pool.from.local.rank(local.rank.op)
  A.current <- base.op$A
  rank.current <- matrix.rank.dense(A.current)
  nullity.current <- n - rank.current
  P <- monomial.design.nd(case$intrinsic.coords, cfg$degree)
  residual.max <- if (case$true.dimension == 1L) 1e-6 else 1e-2
  selected.records <- list()
  candidate.records <- list()
  selected.keys <- character()
  max.steps <- min(as.integer(cfg$skeleton.max.cross.paths), nrow(pool))
  if (max.steps > 0L && nrow(pool) > 0L) {
    for (step in seq_len(max.steps)) {
      if (nullity.current <= expected.nullity) break
      remaining <- pool[!pool$path.key %in% selected.keys, , drop = FALSE]
      if (nrow(remaining) < 1L) break
      scored <- vector("list", nrow(remaining))
      for (idx in seq_len(nrow(remaining))) {
        path <- path.string.to.integer(remaining$path[idx])
        A.block <- candidate.block.for.path(
          path,
          graph$adj.list,
          graph$weight.list,
          n = n,
          degree = cfg$degree,
          row.normalization = "unit.l2"
        )
        if (nrow(A.block) < 1L) next
        polynomial.residual <- frobenius.residual.ratio(A.block, P)
        if (!is.finite(polynomial.residual)) polynomial.residual <- Inf
        rank.after <- matrix.rank.dense(rbind(A.current, A.block))
        global.rank.gain <- max(0L, rank.after - rank.current)
        score <- global.rank.gain +
          alpha * remaining$rank.gain[idx] +
          beta * remaining$coverage.gain[idx] -
          gamma * remaining$path.length[idx] -
          eta * polynomial.residual
        scored[[idx]] <- data.frame(
          step = step,
          candidate.pool.index = idx,
          center = remaining$center[idx],
          candidate.type = if ("candidate.type" %in% names(remaining)) {
            remaining$candidate.type[idx]
          } else {
            NA_character_
          },
          path = remaining$path[idx],
          path.key = remaining$path.key[idx],
          path.length = remaining$path.length[idx],
          n.block.rows = nrow(A.block),
          rank.before = rank.current,
          rank.after = rank.after,
          global.rank.gain = global.rank.gain,
          local.rank.gain = remaining$rank.gain[idx],
          coverage.gain = remaining$coverage.gain[idx],
          polynomial.residual = polynomial.residual,
          score = score,
          selected = FALSE,
          rejected.residual = polynomial.residual > residual.max,
          stringsAsFactors = FALSE
        )
      }
      step.table <- rbind.fill.base(scored)
      if (nrow(step.table) < 1L) break
      candidate.records[[length(candidate.records) + 1L]] <- step.table
      eligible <- subset(step.table, !rejected.residual & global.rank.gain > 0)
      if (nrow(eligible) < 1L) break
      eligible <- eligible[order(
        -eligible$global.rank.gain,
        -eligible$local.rank.gain,
        -eligible$coverage.gain,
        -eligible$score,
        eligible$path.length,
        eligible$path.key
      ), , drop = FALSE]
      chosen <- eligible[1L, , drop = FALSE]
      chosen$selected <- TRUE
      selected.records[[length(selected.records) + 1L]] <- chosen
      selected.keys <- c(selected.keys, chosen$path.key)
      A.block <- candidate.block.for.path(
        path.string.to.integer(chosen$path),
        graph$adj.list,
        graph$weight.list,
        n = n,
        degree = cfg$degree,
        row.normalization = "unit.l2"
      )
      A.current <- rbind(A.current, A.block)
      rank.current <- chosen$rank.after
      nullity.current <- n - rank.current
    }
  }
  nullity <- estimate.nullity(A.current)
  score <- score.rank.calibration.operator(A.current, case$intrinsic.coords, cfg$degree)
  direction <- local.rank.op$skeleton.direction.table %||% data.frame()
  local.candidates <- local.rank.op$skeleton.connector.candidate.table %||% data.frame()
  skeleton.table <- base.op$skeleton.table %||% data.frame()
  coverage <- base.op$skeleton.coverage %||% data.frame()
  selection.table <- rbind.fill.base(selected.records)
  global.candidates <- rbind.fill.base(candidate.records)
  if (nrow(global.candidates) > 0L && nrow(selection.table) > 0L) {
    global.candidates$selected <- global.candidates$path.key %in% selection.table$path.key
  }
  summary <- data.frame(
    case.id = case$case.id,
    family = case$family,
    graph.type = case$graph.type,
    n = n,
    true.dimension = case$true.dimension,
    expected.nullity = expected.nullity,
    row.family = cfg$row.family,
    path.family = "geodesic.skeleton.cross.enriched",
    skeleton.cross.rule = "global.rank",
    skeleton.cross.candidate.pool = cfg$skeleton.cross.candidate.pool,
    degree = cfg$degree,
    skeleton.max.paths = cfg$skeleton.max.paths,
    disk.k = cfg$disk.k,
    skeleton.max.cross.paths = cfg$skeleton.max.cross.paths,
    skeleton.cross.iterations = cfg$skeleton.cross.iterations,
    skeleton.cross.rank.tol = cfg$skeleton.cross.rank.tol,
    skeleton.cross.length.penalty = cfg$skeleton.cross.length.penalty,
    skeleton.cross.coverage.reward = cfg$skeleton.cross.coverage.reward,
    status = "ok",
    message = "",
    n.operator.rows = nrow(A.current),
    nullity.1e.6 = nullity$table$nullity[nullity$table$rel.tol == 1e-6],
    intrinsic.residual = score$intrinsic.residual,
    control.residual = score$control.residual,
    n.skeleton.paths = nrow(skeleton.table),
    n.connectors = nrow(selection.table),
    n.direction.disks = nrow(direction),
    mean.deficit.before = if (nrow(direction)) {
      mean(direction$direction.deficit.before, na.rm = TRUE)
    } else {
      NA_real_
    },
    mean.deficit.after = if (nrow(direction)) {
      mean(direction$direction.deficit.after, na.rm = TRUE)
    } else {
      NA_real_
    },
    max.deficit.before = if (nrow(direction)) {
      max(direction$direction.deficit.before, na.rm = TRUE)
    } else {
      NA_real_
    },
    max.deficit.after = if (nrow(direction)) {
      max(direction$direction.deficit.after, na.rm = TRUE)
    } else {
      NA_real_
    },
    n.rank.connectors.added = nrow(selection.table),
    n.connector.candidates = nrow(pool),
    n.global.candidates.scored = nrow(global.candidates),
    global.rank.initial = matrix.rank.dense(base.op$A),
    global.rank.final = rank.current,
    global.rank.gain = rank.current - matrix.rank.dense(base.op$A),
    global.nullity.initial = n - matrix.rank.dense(base.op$A),
    global.nullity.final = n - rank.current,
    median.skeleton.coverage = if (nrow(coverage)) {
      stats::median(coverage$skeleton.coverage, na.rm = TRUE)
    } else {
      NA_real_
    },
    min.skeleton.coverage = if (nrow(coverage)) {
      min(coverage$skeleton.coverage, na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  tag <- summary[, c(
    "case.id", "family", "graph.type", "n", "true.dimension", "row.family",
    "path.family", "skeleton.cross.rule", "skeleton.cross.candidate.pool",
    "degree", "skeleton.max.paths",
    "disk.k", "skeleton.max.cross.paths", "skeleton.cross.iterations",
    "skeleton.cross.rank.tol", "skeleton.cross.length.penalty",
    "skeleton.cross.coverage.reward"
  ), drop = FALSE]
  add.tag <- function(tab) {
    if (is.null(tab) || nrow(tab) < 1L) return(data.frame())
    cbind(tag[rep(1L, nrow(tab)), , drop = FALSE], tab)
  }
  list(
    summary = summary,
    direction = add.tag(direction),
    connector.candidates = add.tag(local.candidates),
    skeleton = add.tag(skeleton.table),
    coverage = add.tag(coverage),
    global.candidates = add.tag(global.candidates),
    global.selection = add.tag(selection.table)
  )
}

rank.calibration.config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  rank.grid <- if (mode == "smoke") {
    expand.grid(
      row.family = "rank.cross",
      degree = 1:2,
      skeleton.max.paths = c(4L, 8L),
      disk.k = c(8L, 14L),
      skeleton.max.cross.paths = c(10L, 20L, 40L, 80L),
      skeleton.cross.iterations = c(1L, 2L, 4L),
      skeleton.cross.rank.tol = 1e-6,
      skeleton.cross.length.penalty = c(0.00, 0.05),
      skeleton.cross.coverage.reward = 0.10,
      skeleton.cross.candidate.pool = "skeleton.membership",
      stringsAsFactors = FALSE
    )
  } else {
    expand.grid(
      row.family = "rank.cross",
      degree = 1:3,
      skeleton.max.paths = c(4L, 8L, 16L, 32L),
      disk.k = c(8L, 14L, 22L),
      skeleton.max.cross.paths = c(10L, 20L, 40L, 80L, 160L),
      skeleton.cross.iterations = c(1L, 2L, 4L),
      skeleton.cross.rank.tol = c(1e-7, 1e-6, 1e-4),
      skeleton.cross.length.penalty = c(0.00, 0.05, 0.10),
      skeleton.cross.coverage.reward = c(0.05, 0.10),
      skeleton.cross.candidate.pool = "skeleton.membership",
      stringsAsFactors = FALSE
    )
  }
  baseline.grid <- expand.grid(
    row.family = c("skeleton.only", "coverage.cross"),
    degree = sort(unique(rank.grid$degree)),
    skeleton.max.paths = if (mode == "smoke") c(8L) else c(16L),
    disk.k = if (mode == "smoke") c(14L) else c(14L, 22L),
    skeleton.max.cross.paths = if (mode == "smoke") 20L else 40L,
    skeleton.cross.iterations = 1L,
    skeleton.cross.rank.tol = 1e-6,
    skeleton.cross.length.penalty = 0.05,
    skeleton.cross.coverage.reward = 0.10,
    skeleton.cross.candidate.pool = "skeleton.membership",
    stringsAsFactors = FALSE
  )
  global.grid <- if (mode == "smoke") {
    expand.grid(
      row.family = c("global.rank.cross", "global.rank.rich"),
      degree = 1:2,
      skeleton.max.paths = c(4L, 8L),
      disk.k = c(8L, 14L),
      skeleton.max.cross.paths = c(10L, 20L),
      skeleton.cross.iterations = 2L,
      skeleton.cross.rank.tol = 1e-6,
      skeleton.cross.length.penalty = 0.00,
      skeleton.cross.coverage.reward = 0.10,
      skeleton.cross.candidate.pool = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    expand.grid(
      row.family = c("global.rank.cross", "global.rank.rich"),
      degree = 1:3,
      skeleton.max.paths = c(4L, 8L, 16L),
      disk.k = c(8L, 14L, 22L),
      skeleton.max.cross.paths = c(10L, 20L, 40L, 80L),
      skeleton.cross.iterations = c(2L, 4L),
      skeleton.cross.rank.tol = 1e-6,
      skeleton.cross.length.penalty = c(0.00, 0.05),
      skeleton.cross.coverage.reward = 0.10,
      skeleton.cross.candidate.pool = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  global.grid$skeleton.cross.candidate.pool <- vapply(
    global.grid$row.family,
    candidate.pool.for.row.family,
    character(1)
  )
  unique(rbind(rank.grid, global.grid, baseline.grid))
}

fit.rank.calibration.config <- function(case,
                                        cfg,
                                        seed = 1L) {
  graph <- case$graph
  if (cfg$row.family %in% c("global.rank.cross", "global.rank.rich")) {
    return(fit.global.rank.calibration.config(case, cfg, seed = seed))
  }
  path.family <- if (cfg$row.family == "skeleton.only") {
    "geodesic.skeleton"
  } else {
    "geodesic.skeleton.cross.enriched"
  }
  cross.rule <- if (cfg$row.family == "rank.cross") {
    "local.dimension.rank"
  } else {
    "coverage"
  }
  min.coverage <- if (cfg$row.family == "rank.cross") 1L else 2L
  op <- build.geodesic.annihilator.operator(
    graph$adj.list,
    graph$weight.list,
    degree = cfg$degree,
    path.family = path.family,
    disk.k = cfg$disk.k,
    max.skeleton.paths = cfg$skeleton.max.paths,
    skeleton.min.coverage = min.coverage,
    skeleton.max.passes = 2L,
    skeleton.max.cross.paths = cfg$skeleton.max.cross.paths,
    skeleton.cross.rule = cross.rule,
    skeleton.cross.embedding.method = "cmdscale",
    skeleton.cross.candidate.pool = cfg$skeleton.cross.candidate.pool,
    skeleton.cross.local.dims = 1:5,
    skeleton.cross.rank.tol = cfg$skeleton.cross.rank.tol,
    skeleton.cross.max.candidates.per.center = 50L,
    skeleton.cross.length.penalty = cfg$skeleton.cross.length.penalty,
    skeleton.cross.coverage.reward = cfg$skeleton.cross.coverage.reward,
    skeleton.cross.iterations = cfg$skeleton.cross.iterations,
    row.normalization = "unit.l2",
    vertex.normalization = "mean",
    seed = seed
  )
  nullity <- estimate.nullity(op$A)
  score <- score.rank.calibration.operator(op$A, case$intrinsic.coords, cfg$degree)
  expected.nullity <- choose(cfg$degree + case$true.dimension, case$true.dimension)
  direction <- op$skeleton.direction.table %||% data.frame()
  connector.candidates <- op$skeleton.connector.candidate.table %||% data.frame()
  skeleton.table <- op$skeleton.table %||% data.frame()
  coverage <- op$skeleton.coverage %||% data.frame()
  summary <- data.frame(
    case.id = case$case.id,
    family = case$family,
    graph.type = case$graph.type,
    n = length(graph$adj.list),
    true.dimension = case$true.dimension,
    expected.nullity = expected.nullity,
    row.family = cfg$row.family,
    path.family = path.family,
    skeleton.cross.rule = cross.rule,
    skeleton.cross.candidate.pool = cfg$skeleton.cross.candidate.pool,
    degree = cfg$degree,
    skeleton.max.paths = cfg$skeleton.max.paths,
    disk.k = cfg$disk.k,
    skeleton.max.cross.paths = cfg$skeleton.max.cross.paths,
    skeleton.cross.iterations = cfg$skeleton.cross.iterations,
    skeleton.cross.rank.tol = cfg$skeleton.cross.rank.tol,
    skeleton.cross.length.penalty = cfg$skeleton.cross.length.penalty,
    skeleton.cross.coverage.reward = cfg$skeleton.cross.coverage.reward,
    status = "ok",
    message = "",
    n.operator.rows = nrow(op$A),
    nullity.1e.6 = nullity$table$nullity[nullity$table$rel.tol == 1e-6],
    intrinsic.residual = score$intrinsic.residual,
    control.residual = score$control.residual,
    n.skeleton.paths = nrow(skeleton.table),
    n.connectors = if (nrow(skeleton.table)) {
      sum(grepl("connector", skeleton.table$skeleton.role))
    } else {
      0L
    },
    n.direction.disks = nrow(direction),
    mean.deficit.before = if (nrow(direction)) {
      mean(direction$direction.deficit.before, na.rm = TRUE)
    } else {
      NA_real_
    },
    mean.deficit.after = if (nrow(direction)) {
      mean(direction$direction.deficit.after, na.rm = TRUE)
    } else {
      NA_real_
    },
    max.deficit.before = if (nrow(direction)) {
      max(direction$direction.deficit.before, na.rm = TRUE)
    } else {
      NA_real_
    },
    max.deficit.after = if (nrow(direction)) {
      max(direction$direction.deficit.after, na.rm = TRUE)
    } else {
      NA_real_
    },
    n.rank.connectors.added = if (nrow(direction)) {
      sum(direction$connector.added, na.rm = TRUE)
    } else {
      0L
    },
    n.connector.candidates = nrow(connector.candidates),
    median.skeleton.coverage = if (nrow(coverage)) {
      stats::median(coverage$skeleton.coverage, na.rm = TRUE)
    } else {
      NA_real_
    },
    min.skeleton.coverage = if (nrow(coverage)) {
      min(coverage$skeleton.coverage, na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  tag <- summary[, c(
    "case.id", "family", "graph.type", "n", "true.dimension", "row.family",
    "path.family", "skeleton.cross.rule", "skeleton.cross.candidate.pool",
    "degree", "skeleton.max.paths",
    "disk.k", "skeleton.max.cross.paths", "skeleton.cross.iterations",
    "skeleton.cross.rank.tol", "skeleton.cross.length.penalty",
    "skeleton.cross.coverage.reward"
  ), drop = FALSE]
  add.tag <- function(tab) {
    if (is.null(tab) || nrow(tab) < 1L) return(data.frame())
    cbind(tag[rep(1L, nrow(tab)), , drop = FALSE], tab)
  }
  list(
    summary = summary,
    direction = add.tag(direction),
    connector.candidates = add.tag(connector.candidates),
    skeleton = add.tag(skeleton.table),
    coverage = add.tag(coverage)
  )
}

run.geodesic.skeleton.rank.calibration <- function(mode = c("smoke", "full"),
                                                   output.dir = file.path(
                                                     dev.dir,
                                                     "reports"
                                                   ),
                                                   seed.base = 20260518L) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  cases <- make.rank.calibration.cases(mode = mode, seed.base = seed.base)
  config <- rank.calibration.config(mode)
  records <- list()
  idx <- 0L
  total <- length(cases) * nrow(config)
  for (case.idx in seq_along(cases)) {
    for (cfg.idx in seq_len(nrow(config))) {
      idx <- idx + 1L
      if (idx == 1L || idx %% 25L == 0L || idx == total) {
        message(sprintf("rank calibration %s run %d/%d", mode, idx, total))
      }
      cfg <- config[cfg.idx, , drop = FALSE]
      records[[idx]] <- tryCatch(
        fit.rank.calibration.config(
          cases[[case.idx]], cfg, seed = seed.base + idx
        ),
        error = function(e) {
          list(
            summary = data.frame(
              case.id = cases[[case.idx]]$case.id,
              family = cases[[case.idx]]$family,
              graph.type = cases[[case.idx]]$graph.type,
              n = length(cases[[case.idx]]$graph$adj.list),
              true.dimension = cases[[case.idx]]$true.dimension,
              expected.nullity = choose(cfg$degree + cases[[case.idx]]$true.dimension,
                                        cases[[case.idx]]$true.dimension),
              row.family = cfg$row.family,
              path.family = NA_character_,
              skeleton.cross.rule = NA_character_,
              skeleton.cross.candidate.pool = cfg$skeleton.cross.candidate.pool,
              degree = cfg$degree,
              skeleton.max.paths = cfg$skeleton.max.paths,
              disk.k = cfg$disk.k,
              skeleton.max.cross.paths = cfg$skeleton.max.cross.paths,
              skeleton.cross.iterations = cfg$skeleton.cross.iterations,
              skeleton.cross.rank.tol = cfg$skeleton.cross.rank.tol,
              skeleton.cross.length.penalty = cfg$skeleton.cross.length.penalty,
              skeleton.cross.coverage.reward = cfg$skeleton.cross.coverage.reward,
              status = "error",
              message = conditionMessage(e),
              n.operator.rows = NA_integer_,
              nullity.1e.6 = NA_integer_,
              intrinsic.residual = NA_real_,
              control.residual = NA_real_,
              n.skeleton.paths = NA_integer_,
              n.connectors = NA_integer_,
              n.direction.disks = NA_integer_,
              mean.deficit.before = NA_real_,
              mean.deficit.after = NA_real_,
              max.deficit.before = NA_real_,
              max.deficit.after = NA_real_,
              n.rank.connectors.added = NA_integer_,
              n.connector.candidates = NA_integer_,
              n.global.candidates.scored = NA_integer_,
              global.rank.initial = NA_integer_,
              global.rank.final = NA_integer_,
              global.rank.gain = NA_integer_,
              global.nullity.initial = NA_integer_,
              global.nullity.final = NA_integer_,
              median.skeleton.coverage = NA_real_,
              min.skeleton.coverage = NA_real_,
              stringsAsFactors = FALSE
            ),
            direction = data.frame(),
            connector.candidates = data.frame(),
            skeleton = data.frame(),
            coverage = data.frame(),
            global.candidates = data.frame(),
            global.selection = data.frame()
          )
        }
      )
    }
  }
  summary <- rbind.fill.base(lapply(records, `[[`, "summary"))
  if (!"status" %in% names(summary)) summary$status <- "ok"
  if (!"message" %in% names(summary)) summary$message <- ""
  direction <- rbind.fill.base(lapply(records, `[[`, "direction"))
  connector.candidates <- rbind.fill.base(lapply(records, `[[`, "connector.candidates"))
  skeleton <- rbind.fill.base(lapply(records, `[[`, "skeleton"))
  coverage <- rbind.fill.base(lapply(records, `[[`, "coverage"))
  global.candidates <- rbind.fill.base(lapply(records, `[[`, "global.candidates"))
  global.selection <- rbind.fill.base(lapply(records, `[[`, "global.selection"))
  diagnostics <- list(
    mode = mode,
    config = config,
    summary = summary,
    direction = direction,
    connector.candidates = connector.candidates,
    skeleton = skeleton,
    coverage = coverage,
    global.candidates = global.candidates,
    global.selection = global.selection,
    created = Sys.time()
  )
  saveRDS(diagnostics, file.path(output.dir,
                                 paste0("geodesic_skeleton_rank_calibration_", mode, ".rds")))
  utils::write.csv(summary, file.path(output.dir,
                                      paste0("geodesic_skeleton_rank_calibration_", mode, "_summary.csv")),
                   row.names = FALSE)
  utils::write.csv(direction, file.path(output.dir,
                                        paste0("geodesic_skeleton_rank_calibration_", mode, "_direction.csv")),
                   row.names = FALSE)
  utils::write.csv(connector.candidates, file.path(output.dir,
                                                   paste0("geodesic_skeleton_rank_calibration_", mode, "_connector_candidates.csv")),
                   row.names = FALSE)
  utils::write.csv(skeleton, file.path(output.dir,
                                       paste0("geodesic_skeleton_rank_calibration_", mode, "_skeleton.csv")),
                   row.names = FALSE)
  utils::write.csv(coverage, file.path(output.dir,
                                       paste0("geodesic_skeleton_rank_calibration_", mode, "_coverage.csv")),
                   row.names = FALSE)
  utils::write.csv(global.candidates, file.path(output.dir,
                                                paste0("geodesic_skeleton_rank_calibration_", mode, "_global_candidates.csv")),
                   row.names = FALSE)
  utils::write.csv(global.selection, file.path(output.dir,
                                               paste0("geodesic_skeleton_rank_calibration_", mode, "_global_selection.csv")),
                   row.names = FALSE)
  diagnostics
}

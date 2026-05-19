repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")

source(file.path(dev.dir, "geodesic_annihilator_experiments.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

make.path.validation.graph <- function(n = 60L, closed = FALSE) {
  x <- seq(0, 1, length.out = n)
  coords <- cbind(x, 0)
  edges <- cbind(seq_len(n - 1L), seq.int(2L, n))
  if (closed) {
    theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
    coords <- cbind(cos(theta), sin(theta))
    edges <- rbind(edges, c(n, 1L))
  }
  graph <- make.graph.from.edges(n, edges, edge.lengths.from.coords(edges, coords))
  graph$coords <- coords
  graph$boundary.score <- if (closed) rep(0, n) else pmin(x, 1 - x)
  graph
}

make.noisy.circle.validation.graph <- function(n = 80L,
                                               sigma = 0.02,
                                               graph.type = c("cycle", "symmetric.knn"),
                                               k = 8L,
                                               seed = 1L) {
  graph.type <- match.arg(graph.type)
  set.seed(seed)
  theta <- sort(stats::runif(n, 0, 2 * pi))
  clean <- cbind(cos(theta), sin(theta))
  coords <- clean + matrix(stats::rnorm(2L * n, sd = sigma), ncol = 2L)
  if (graph.type == "cycle") {
    edges <- cbind(seq_len(n), c(seq.int(2L, n), 1L))
    graph <- make.graph.from.edges(n, edges, edge.lengths.from.coords(edges, coords))
  } else {
    graph <- make.symmetric.knn.graph(coords, k = k, metric.coords = coords)
  }
  graph$coords <- coords
  graph$latent.theta <- theta
  graph$boundary.score <- rep(0, n)
  graph
}

make.grid.validation.graph <- function(nx = 8L, ny = 8L,
                                       warp = c("none", "smooth", "nonuniform")) {
  warp <- match.arg(warp)
  gx <- seq(0, 1, length.out = nx)
  gy <- seq(0, 1, length.out = ny)
  if (warp == "nonuniform") {
    gx <- gx^1.45
    gy <- 1 - (1 - gy)^1.35
  }
  U <- as.matrix(expand.grid(gx, gy))
  colnames(U) <- c("x", "y")
  coords <- U
  if (warp == "smooth") {
    coords <- cbind(
      U[, 1L] + 0.08 * sin(2 * pi * U[, 2L]),
      U[, 2L] + 0.06 * sin(2 * pi * U[, 1L])
    )
  }
  idx <- function(i, j) (j - 1L) * nx + i
  edges <- matrix(NA_integer_, nrow = 0L, ncol = 2L)
  for (j in seq_len(ny)) {
    for (i in seq_len(nx)) {
      if (i < nx) edges <- rbind(edges, c(idx(i, j), idx(i + 1L, j)))
      if (j < ny) edges <- rbind(edges, c(idx(i, j), idx(i, j + 1L)))
    }
  }
  graph <- make.graph.from.edges(nrow(U), edges, edge.lengths.from.coords(edges, coords))
  graph$coords <- coords
  graph$param <- U
  graph$boundary.score <- pmin(U[, 1L], 1 - U[, 1L], U[, 2L], 1 - U[, 2L])
  graph
}

make.cube.lattice.validation.graph <- function(nx = 5L, ny = 5L, nz = 5L) {
  gx <- seq(0, 1, length.out = nx)
  gy <- seq(0, 1, length.out = ny)
  gz <- seq(0, 1, length.out = nz)
  X <- as.matrix(expand.grid(gx, gy, gz))
  colnames(X) <- c("x", "y", "z")
  idx <- function(i, j, k) (k - 1L) * nx * ny + (j - 1L) * nx + i
  edges <- matrix(NA_integer_, nrow = 0L, ncol = 2L)
  for (k in seq_len(nz)) {
    for (j in seq_len(ny)) {
      for (i in seq_len(nx)) {
        if (i < nx) edges <- rbind(edges, c(idx(i, j, k), idx(i + 1L, j, k)))
        if (j < ny) edges <- rbind(edges, c(idx(i, j, k), idx(i, j + 1L, k)))
        if (k < nz) edges <- rbind(edges, c(idx(i, j, k), idx(i, j, k + 1L)))
      }
    }
  }
  graph <- make.graph.from.edges(nrow(X), edges, edge.lengths.from.coords(edges, X))
  graph$coords <- X
  graph$boundary.score <- do.call(pmin, as.data.frame(cbind(X, 1 - X)))
  graph
}

make.sampled.cube.validation.graph <- function(n = 100L, k = 10L, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::runif(3L * n), ncol = 3L)
  graph <- make.symmetric.knn.graph(X, k = k, metric.coords = X)
  graph$coords <- X
  graph$boundary.score <- do.call(pmin, as.data.frame(cbind(X, 1 - X)))
  graph
}

sample.ball.validation <- function(n, radius = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z <- matrix(stats::rnorm(3L * n), ncol = 3L)
  Z <- Z / pmax(sqrt(rowSums(Z^2)), .Machine$double.eps)
  radii <- radius * stats::runif(n)^(1 / 3)
  Z * radii
}

sample.cube.validation <- function(n, radius = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(stats::runif(3L * n, min = -radius, max = radius), ncol = 3L)
}

make.quadform3d.validation.graph <- function(n = 80L,
                                             index.k = 3L,
                                             coefficients = c(1, 1, 1),
                                             domain.shape = c("ball", "cube"),
                                             domain.radius = 1,
                                             seed = 1L,
                                             graph.type = c("delaunay.1skeleton",
                                                            "symmetric.knn"),
                                             k = 10L) {
  domain.shape <- match.arg(domain.shape)
  graph.type <- match.arg(graph.type)
  X.param <- if (domain.shape == "ball") {
    sample.ball.validation(n, radius = domain.radius, seed = seed)
  } else {
    sample.cube.validation(n, radius = domain.radius, seed = seed)
  }
  X.embed <- quadform.embed(X.param, index.k = index.k, coefficients = coefficients)
  if (graph.type == "delaunay.1skeleton") {
    edges <- delaunay.edges(X.param)
    weights <- quadform.edge.lengths(
      X.param[edges[, 1L], , drop = FALSE],
      X.param[edges[, 2L], , drop = FALSE],
      index.k = index.k,
      coefficients = coefficients
    )
    graph <- make.graph.from.edges(nrow(X.param), edges, weights)
  } else {
    graph <- make.symmetric.knn.graph(X.embed, k = k, metric.coords = X.embed)
  }
  graph$coords <- X.embed
  graph$param <- X.param
  graph$quadform.index.k <- index.k
  graph$quadform.coefficients <- coefficients
  graph$domain.shape <- domain.shape
  graph$boundary.score <- if (domain.shape == "ball") {
    domain.radius - sqrt(rowSums(X.param^2))
  } else {
    domain.radius - do.call(pmax, as.data.frame(abs(X.param)))
  }
  graph
}

make.local.dimension.case <- function(case.id,
                                      family,
                                      true.dimension,
                                      expected.class = as.character(true.dimension),
                                      graph,
                                      dataset = family,
                                      graph.type = "custom",
                                      transition = FALSE,
                                      noise.sigma = NA_real_,
                                      seed = NA_integer_) {
  list(
    case.id = case.id,
    family = family,
    dataset = dataset,
    graph.type = graph.type,
    true.dimension = true.dimension,
    expected.class = expected.class,
    transition = transition,
    noise.sigma = noise.sigma,
    seed = seed,
    graph = graph
  )
}

make.local.dimension.validation.cases <- function(mode = c("smoke", "full"),
                                                  seed.base = 20260517L) {
  mode <- match.arg(mode)
  cases <- list(
    make.local.dimension.case(
      "path_60", "path", 1L, graph = make.path.validation.graph(60L, closed = FALSE),
      graph.type = "path"
    ),
    make.local.dimension.case(
      "cycle_72", "cycle", 1L, graph = make.path.validation.graph(72L, closed = TRUE),
      graph.type = "cycle"
    ),
    make.local.dimension.case(
      "grid_8x8", "grid", 2L, graph = make.grid.validation.graph(8L, 8L, "none"),
      graph.type = "grid"
    ),
    make.local.dimension.case(
      "warped_grid_8x8", "warped grid", 2L,
      graph = make.grid.validation.graph(8L, 8L, "smooth"),
      graph.type = "grid"
    ),
    make.local.dimension.case(
      "cube_lattice_5x5x5", "cube lattice", 3L,
      graph = make.cube.lattice.validation.graph(5L, 5L, 5L),
      graph.type = "grid"
    )
  )
  noisy.sigmas <- if (mode == "smoke") c(0, 0.04, 0.10) else c(0, 0.015, 0.03, 0.06, 0.10, 0.16)
  for (sigma in noisy.sigmas) {
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      sprintf("noisy_circle_sknn_sigma_%0.3f", sigma),
      "noisy circle",
      1L,
      expected.class = "1-to-2 transition",
      graph = make.noisy.circle.validation.graph(
        n = if (mode == "smoke") 72L else 100L,
        sigma = sigma,
        graph.type = "symmetric.knn",
        k = if (mode == "smoke") 6L else 8L,
        seed = seed.base + round(1000 * sigma)
      ),
      graph.type = "symmetric.knn",
      transition = TRUE,
      noise.sigma = sigma,
      seed = seed.base + round(1000 * sigma)
    )
  }
  if (mode == "smoke") {
    U <- sample.square.points(70L, seed.base + 101L)
    graph <- make.experimental.graph(U, "delaunay.1skeleton", metric.coords = U)
    graph$coords <- U
    graph$boundary.score <- pmin(U[, 1L], 1 - U[, 1L], U[, 2L], 1 - U[, 2L])
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "sampled_square_delaunay_seed_1",
      "sampled square", 2L, graph = graph,
      graph.type = "delaunay.1skeleton", seed = 1L
    )
    X3 <- embed.quadric(U, "paraboloid")
    qgraph <- make.experimental.graph(U, "delaunay.1skeleton", metric.coords = X3)
    qgraph$coords <- X3
    qgraph$param <- U
    qgraph$boundary.score <- graph$boundary.score
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "paraboloid_delaunay_seed_1",
      "paraboloid", 2L, graph = qgraph,
      graph.type = "delaunay.1skeleton", seed = 1L
    )
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "quadform3d_ball_k3_c111_seed_1",
      "3D quadform positive", 3L,
      graph = make.quadform3d.validation.graph(
        n = 80L,
        index.k = 3L,
        coefficients = c(1, 1, 1),
        domain.shape = "ball",
        seed = seed.base + 501L,
        graph.type = "delaunay.1skeleton"
      ),
      graph.type = "delaunay.1skeleton",
      seed = 1L
    )
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "quadform3d_ball_k1_c124_seed_1",
      "3D quadform mixed", 3L,
      graph = make.quadform3d.validation.graph(
        n = 80L,
        index.k = 1L,
        coefficients = c(1, 2, 4),
        domain.shape = "ball",
        seed = seed.base + 502L,
        graph.type = "delaunay.1skeleton"
      ),
      graph.type = "delaunay.1skeleton",
      seed = 1L
    )
  }
  if (mode == "full") {
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "nonuniform_grid_9x9", "nonuniform grid", 2L,
      graph = make.grid.validation.graph(9L, 9L, "nonuniform"),
      graph.type = "grid"
    )
    for (seed in 1:2) {
      U <- sample.square.points(90L, seed.base + seed)
      graph <- make.experimental.graph(U, "delaunay.1skeleton", metric.coords = U)
      graph$coords <- U
      graph$boundary.score <- pmin(U[, 1L], 1 - U[, 1L], U[, 2L], 1 - U[, 2L])
      cases[[length(cases) + 1L]] <- make.local.dimension.case(
        sprintf("sampled_square_delaunay_seed_%d", seed),
        "sampled square", 2L, graph = graph,
        graph.type = "delaunay.1skeleton", seed = seed
      )
      for (surface in c("paraboloid", "saddle")) {
        X3 <- embed.quadric(U, surface)
        qgraph <- make.experimental.graph(U, "delaunay.1skeleton", metric.coords = X3)
        qgraph$coords <- X3
        qgraph$param <- U
        qgraph$boundary.score <- graph$boundary.score
        cases[[length(cases) + 1L]] <- make.local.dimension.case(
          sprintf("%s_delaunay_seed_%d", surface, seed),
          surface, 2L, graph = qgraph,
          graph.type = "delaunay.1skeleton", seed = seed
        )
      }
    }
    cases[[length(cases) + 1L]] <- make.local.dimension.case(
      "sampled_cube_sknn_seed_1", "sampled cube", 3L,
      graph = make.sampled.cube.validation.graph(110L, k = 10L, seed = seed.base + 301L),
      graph.type = "symmetric.knn", seed = 1L
    )
    quad3.specs <- data.frame(
      case.prefix = c("quadform3d_ball_k3_c111", "quadform3d_ball_k1_c124",
                      "quadform3d_cube_k3_c124", "quadform3d_cube_k1_c144"),
      family = c("3D quadform positive", "3D quadform mixed",
                 "3D quadform positive", "3D quadform mixed"),
      domain.shape = c("ball", "ball", "cube", "cube"),
      index.k = c(3L, 1L, 3L, 1L),
      coeff.label = c("c111", "c124", "c124", "c144"),
      stringsAsFactors = FALSE
    )
    coeff.map.3d <- list(c111 = c(1, 1, 1), c124 = c(1, 2, 4), c144 = c(1, 4, 4))
    for (ii in seq_len(nrow(quad3.specs))) {
      for (seed in 1:2) {
        spec <- quad3.specs[ii, ]
        cases[[length(cases) + 1L]] <- make.local.dimension.case(
          sprintf("%s_seed_%d", spec$case.prefix, seed),
          spec$family, 3L,
          graph = make.quadform3d.validation.graph(
            n = 110L,
            index.k = spec$index.k,
            coefficients = coeff.map.3d[[spec$coeff.label]],
            domain.shape = spec$domain.shape,
            seed = seed.base + 600L + 10L * ii + seed,
            graph.type = "delaunay.1skeleton"
          ),
          graph.type = "delaunay.1skeleton",
          seed = seed
        )
      }
    }
  }
  cases
}

disk.vertices.by.k <- function(D, center, disk.k) {
  d <- D[center, ]
  finite <- which(is.finite(d))
  finite <- finite[order(d[finite])]
  finite[seq_len(min(length(finite), disk.k + 1L))]
}

run.case.local.dimension.validation <- function(case,
                                                method = "cmdscale",
                                                disk.k = 18L,
                                                centers.max = 40L,
                                                dims = 1:5,
                                                seed = 1L) {
  graph <- case$graph
  n <- length(graph$adj.list)
  D <- shortest.path(graph$adj.list, graph$weight.list, seq_len(n))
  set.seed(seed)
  centers <- seq_len(n)
  if (is.finite(centers.max) && centers.max > 0L && centers.max < n) {
    boundary.score <- graph$boundary.score %||% rep(1, n)
    interior <- centers[boundary.score > stats::median(boundary.score, na.rm = TRUE)]
    boundary <- centers[boundary.score <= stats::median(boundary.score, na.rm = TRUE)]
    n.interior <- min(length(interior), ceiling(centers.max / 2))
    n.boundary <- min(length(boundary), centers.max - n.interior)
    centers <- sort(c(
      if (n.interior) sample(interior, n.interior) else integer(),
      if (n.boundary) sample(boundary, n.boundary) else integer()
    ))
  }
  records <- vector("list", length(centers))
  for (idx in seq_along(centers)) {
    center <- centers[idx]
    disk.vertices <- disk.vertices.by.k(D, center, disk.k)
    out <- local.dimension.diagnostics(
      adj.list = graph$adj.list,
      weight.list = graph$weight.list,
      D = D,
      center = center,
      disk.vertices = disk.vertices,
      method = method,
      dims = dims,
      seed = seed + idx
    )
    out$case.id <- case$case.id
    out$family <- case$family
    out$dataset <- case$dataset
    out$graph.type <- case$graph.type
    out$true.dimension <- case$true.dimension
    out$expected.class <- case$expected.class
    out$transition <- case$transition
    out$noise.sigma <- case$noise.sigma
    out$case.seed <- case$seed
    out$disk.k <- disk.k
    out$center.boundary.score <- (graph$boundary.score %||% rep(NA_real_, n))[center]
    out$center.region <- if (is.na(out$center.boundary.score[1L])) {
      "unknown"
    } else if (out$center.boundary.score[1L] <= stats::quantile(graph$boundary.score, 0.35, na.rm = TRUE)) {
      "boundary"
    } else {
      "interior"
    }
    records[[idx]] <- out
  }
  do.call(rbind, records)
}

summarize.local.dimension.validation <- function(local.table) {
  center.table <- unique(local.table[
    , c("case.id", "family", "dataset", "graph.type", "true.dimension",
        "expected.class", "transition", "noise.sigma", "case.seed",
        "center", "center.region", "local.embedding.method",
        "selected.local.dim", "best.local.dim", "dimension.cap.hit",
        "dimension.selection.rule", "mds.spectrum.local.dim",
        "graph.spectrum.local.dim", "consensus.local.dim")
  ])
  center.table$correct.strict <- with(
    center.table,
    !transition & !is.na(selected.local.dim) & selected.local.dim == true.dimension
  )
  center.table$noise.sigma.summary <- ifelse(
    is.na(center.table$noise.sigma),
    -1,
    center.table$noise.sigma
  )
  summary <- aggregate(
    cbind(correct.strict, dimension.cap.hit) ~
      case.id + family + dataset + graph.type + true.dimension + expected.class +
      transition + noise.sigma.summary + local.embedding.method,
    data = center.table,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  names(summary)[names(summary) == "noise.sigma.summary"] <- "noise.sigma"
  summary$noise.sigma[summary$noise.sigma < 0] <- NA_real_
  count <- aggregate(
    selected.local.dim ~ case.id + local.embedding.method,
    data = center.table,
    FUN = length
  )
  names(count)[names(count) == "selected.local.dim"] <- "n.centers"
  summary <- merge(summary, count, by = c("case.id", "local.embedding.method"),
                   all.x = TRUE, sort = FALSE)
  names(summary)[names(summary) == "correct.strict"] <- "strict.accuracy"
  names(summary)[names(summary) == "dimension.cap.hit"] <- "cap.hit.rate"
  list(center.table = center.table, summary = summary)
}

run.local.dimension.validation <- function(mode = c("smoke", "full"),
                                           output.dir = file.path(
                                             dev.dir,
                                             "local-dimension-validation",
                                             "reports"
                                           ),
                                           method = c("cmdscale",
                                                      "metric.mds.edge.kk",
                                                      "hybrid.edge.kk"),
                                           seed.base = 20260517L) {
  mode <- match.arg(mode)
  method <- match.arg(method)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  cases <- make.local.dimension.validation.cases(mode, seed.base)
  disk.k <- if (mode == "smoke") 16L else 22L
  centers.max <- if (mode == "smoke") 24L else 48L
  records <- vector("list", length(cases))
  for (idx in seq_along(cases)) {
    message(sprintf(
      "local dimension validation %s case %d/%d: %s",
      mode, idx, length(cases), cases[[idx]]$case.id
    ))
    records[[idx]] <- tryCatch(
        run.case.local.dimension.validation(
        cases[[idx]],
        method = method,
        disk.k = disk.k,
        centers.max = centers.max,
        dims = 1:5,
        seed = seed.base + idx
      ),
      error = function(e) {
        data.frame(
          case.id = cases[[idx]]$case.id,
          family = cases[[idx]]$family,
          dataset = cases[[idx]]$dataset,
          graph.type = cases[[idx]]$graph.type,
          true.dimension = cases[[idx]]$true.dimension,
          expected.class = cases[[idx]]$expected.class,
          transition = cases[[idx]]$transition,
          noise.sigma = cases[[idx]]$noise.sigma,
          local.embedding.method = method,
          embedding.dim = 1:5,
          selected.local.dim = NA_integer_,
          best.local.dim = NA_integer_,
          dimension.cap.hit = NA,
          status = "case.error",
          message = conditionMessage(e),
          edge.rrmse = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    )
  }
  local.table <- do.call(rbind, records)
  summaries <- summarize.local.dimension.validation(local.table)
  diagnostics <- list(
    mode = mode,
    local.dimension = local.table,
    center.selection = summaries$center.table,
    summary = summaries$summary,
    created = Sys.time(),
    control = list(method = method, disk.k = disk.k, centers.max = centers.max)
  )
  saveRDS(diagnostics, file.path(output.dir,
                                 paste0("local_dimension_validation_", mode, ".rds")))
  utils::write.csv(local.table, file.path(output.dir,
                                          paste0("local_dimension_validation_", mode, "_local_dimension.csv")),
                   row.names = FALSE)
  utils::write.csv(summaries$center.table, file.path(output.dir,
                                                     paste0("local_dimension_validation_", mode, "_center_selection.csv")),
                   row.names = FALSE)
  utils::write.csv(summaries$summary, file.path(output.dir,
                                                paste0("local_dimension_validation_", mode, "_summary.csv")),
                   row.names = FALSE)
  diagnostics
}

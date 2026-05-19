repo.dir <- "/Users/pgajer/current_projects/gflow"

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
} else {
  library(gflow)
}

required.packages <- c("Matrix", "FNN", "geometry", "ggplot2")
missing.packages <- required.packages[
  !vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing.packages) > 0L) {
  stop("Missing required packages: ", paste(missing.packages, collapse = ", "))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

add.undirected.edge <- function(adj, weights, i, j, weight = 1) {
  adj[[i]] <- c(adj[[i]], j)
  weights[[i]] <- c(weights[[i]], weight)
  adj[[j]] <- c(adj[[j]], i)
  weights[[j]] <- c(weights[[j]], weight)
  list(adj = adj, weights = weights)
}

make.graph.from.edges <- function(n, edges, weights = rep(1, nrow(edges))) {
  if (nrow(edges) < 1L) {
    stop("At least one edge is required.")
  }
  edges <- as.matrix(edges)
  edge.keys <- apply(edges, 1L, function(z) paste(sort(z), collapse = ":"))
  keep <- !duplicated(edge.keys) & edges[, 1L] != edges[, 2L]
  edges <- edges[keep, , drop = FALSE]
  weights <- weights[keep]
  adj <- vector("list", n)
  w <- vector("list", n)
  for (idx in seq_len(nrow(edges))) {
    out <- add.undirected.edge(adj, w, edges[idx, 1L], edges[idx, 2L], weights[idx])
    adj <- out$adj
    w <- out$weights
  }
  list(adj.list = lapply(adj, as.integer),
       weight.list = lapply(w, as.double),
       edges = edges,
       edge.weights = as.double(weights))
}

edge.lengths.from.coords <- function(edges, coords) {
  sqrt(rowSums((coords[edges[, 1L], , drop = FALSE] -
                  coords[edges[, 2L], , drop = FALSE])^2))
}

delaunay.edges <- function(U) {
  tri <- geometry::delaunayn(U, options = "QJ")
  edge.mat <- do.call(rbind, lapply(seq_len(nrow(tri)), function(i) {
    cmb <- utils::combn(tri[i, ], 2L)
    t(apply(cmb, 2L, sort))
  }))
  edge.mat <- unique(edge.mat)
  edge.mat[order(edge.mat[, 1L], edge.mat[, 2L]), , drop = FALSE]
}

make.delaunay.graph <- function(U, metric.coords = U) {
  edges <- delaunay.edges(U)
  make.graph.from.edges(nrow(U), edges, edge.lengths.from.coords(edges, metric.coords))
}

make.symmetric.knn.graph <- function(U, k = 8L, metric.coords = U) {
  knn <- FNN::get.knn(U, k = k)
  edges <- do.call(rbind, lapply(seq_len(nrow(U)), function(i) {
    cbind(i, knn$nn.index[i, ])
  }))
  edges <- t(apply(edges, 1L, sort))
  edges <- unique(edges)
  make.graph.from.edges(nrow(U), edges, edge.lengths.from.coords(edges, metric.coords))
}

make.radius.graph <- function(U, radius, metric.coords = U) {
  D <- as.matrix(stats::dist(U))
  edges <- which(upper.tri(D) & D <= radius, arr.ind = TRUE)
  if (nrow(edges) < 1L) {
    stop("Radius graph has no edges.")
  }
  make.graph.from.edges(nrow(U), edges, edge.lengths.from.coords(edges, metric.coords))
}

make.cknn.graph <- function(U, k = 8L, delta = 1.35, metric.coords = U) {
  D <- as.matrix(stats::dist(U))
  kth <- apply(D, 1L, function(z) sort(z)[k + 1L])
  thresh <- delta * sqrt(outer(kth, kth))
  edges <- which(upper.tri(D) & D <= thresh, arr.ind = TRUE)
  if (nrow(edges) < 1L) {
    stop("cKNN graph has no edges.")
  }
  make.graph.from.edges(nrow(U), edges, edge.lengths.from.coords(edges, metric.coords))
}

make.experimental.graph <- function(U,
                                    graph.type = c("delaunay.1skeleton",
                                                   "symmetric.knn",
                                                   "radius",
                                                   "cknn"),
                                    metric.coords = U,
                                    k = 8L,
                                    radius.multiplier = 2.25,
                                    cknn.delta = 1.35) {
  graph.type <- match.arg(graph.type)
  if (graph.type == "delaunay.1skeleton") {
    graph <- make.delaunay.graph(U, metric.coords)
  } else if (graph.type == "symmetric.knn") {
    graph <- make.symmetric.knn.graph(U, k = k, metric.coords = metric.coords)
  } else if (graph.type == "radius") {
    radius <- radius.multiplier * sqrt(log(nrow(U)) / nrow(U))
    graph <- make.radius.graph(U, radius = radius, metric.coords = metric.coords)
  } else {
    graph <- make.cknn.graph(U, k = k, delta = cknn.delta, metric.coords = metric.coords)
  }
  graph$graph.type <- graph.type
  graph
}

path.edge.length <- function(adj.list, weight.list, i, j) {
  pos <- match(j, adj.list[[i]])
  if (is.na(pos)) {
    stop("Path uses a missing graph edge: ", i, "-", j)
  }
  weight.list[[i]][pos]
}

path.cumulative.distances <- function(vertices, adj.list, weight.list) {
  if (length(vertices) == 1L) return(0)
  increments <- vapply(seq_len(length(vertices) - 1L), function(idx) {
    path.edge.length(adj.list, weight.list, vertices[idx], vertices[idx + 1L])
  }, numeric(1))
  c(0, cumsum(increments))
}

deterministic.shortest.path <- function(adj.list,
                                        weight.list,
                                        source,
                                        target,
                                        tol = 1e-12) {
  n <- length(adj.list)
  source <- as.integer(source)
  target <- as.integer(target)
  if (source == target) return(source)
  dist <- rep(Inf, n)
  prev <- rep(NA_integer_, n)
  visited <- rep(FALSE, n)
  dist[source] <- 0
  for (step in seq_len(n)) {
    available <- which(!visited & is.finite(dist))
    if (length(available) < 1L) break
    u <- available[order(dist[available], available)][1L]
    if (u == target) break
    visited[u] <- TRUE
    nbr.order <- order(adj.list[[u]])
    nbrs <- adj.list[[u]][nbr.order]
    weights <- weight.list[[u]][nbr.order]
    for (idx in seq_along(nbrs)) {
      v <- nbrs[idx]
      if (visited[v]) next
      alt <- dist[u] + weights[idx]
      if (alt < dist[v] - tol ||
          (abs(alt - dist[v]) <= tol && (is.na(prev[v]) || u < prev[v]))) {
        dist[v] <- alt
        prev[v] <- u
      }
    }
  }
  if (!is.finite(dist[target])) {
    return(integer())
  }
  path <- target
  current <- target
  guard <- 0L
  while (current != source && guard <= n) {
    current <- prev[current]
    if (is.na(current)) return(integer())
    path <- c(current, path)
    guard <- guard + 1L
  }
  as.integer(path)
}

make.composite.path <- function(path.i, path.j, adj.list, weight.list) {
  xi <- path.cumulative.distances(path.i, adj.list, weight.list)
  xj <- path.cumulative.distances(path.j, adj.list, weight.list)
  vertices <- c(rev(path.i), path.j[-1L])
  x <- c(-rev(xi), xj[-1L])
  list(vertices = as.integer(vertices), x = as.double(x))
}

path.key <- function(vertices) {
  forward <- paste(vertices, collapse = "-")
  reverse <- paste(rev(vertices), collapse = "-")
  if (forward <= reverse) forward else reverse
}

poly.design.1d <- function(x, degree) {
  outer(as.double(x), seq_len(degree + 1L) - 1L, `^`)
}

annihilator.rows.for.path <- function(vertices,
                                      x,
                                      n.vertices,
                                      degree,
                                      row.normalization = c("unit.l2", "none")) {
  row.normalization <- match.arg(row.normalization)
  keep <- !duplicated(vertices)
  vertices <- vertices[keep]
  x <- x[keep]
  if (length(vertices) < degree + 2L) {
    return(NULL)
  }
  P <- poly.design.1d(x, degree)
  q <- qr(P, LAPACK = TRUE)
  rank <- q$rank
  if (rank >= length(vertices)) {
    return(NULL)
  }
  Q <- qr.Q(q, complete = TRUE)
  basis <- Q[, seq.int(rank + 1L, length(vertices)), drop = FALSE]
  rows <- vector("list", ncol(basis))
  for (idx in seq_len(ncol(basis))) {
    coefs <- basis[, idx]
    if (row.normalization == "unit.l2") {
      norm <- sqrt(sum(coefs^2))
      if (is.finite(norm) && norm > 0) coefs <- coefs / norm
    }
    rows[[idx]] <- list(j = as.integer(vertices), x = as.double(coefs))
  }
  rows
}

monomial.design.2d <- function(U, degree, exact.degree = FALSE) {
  powers <- do.call(rbind, lapply(0:degree, function(total) {
    cbind(0:total, total - 0:total)
  }))
  if (exact.degree) {
    powers <- powers[rowSums(powers) == degree, , drop = FALSE]
  }
  X <- sapply(seq_len(nrow(powers)), function(idx) {
    U[, 1L]^powers[idx, 1L] * U[, 2L]^powers[idx, 2L]
  })
  colnames(X) <- sprintf("u^%d v^%d", powers[, 1L], powers[, 2L])
  X
}

monomial.powers.nd <- function(dim, degree) {
  if (dim < 1L) return(matrix(0L, nrow = 1L, ncol = 0L))
  grid <- expand.grid(rep(list(0:degree), dim))
  powers <- as.matrix(grid[rowSums(grid) <= degree, , drop = FALSE])
  powers <- powers[order(rowSums(powers), do.call(order, as.data.frame(powers))), , drop = FALSE]
  storage.mode(powers) <- "integer"
  powers
}

monomial.design.nd <- function(X, degree) {
  X <- as.matrix(X)
  powers <- monomial.powers.nd(ncol(X), degree)
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

frobenius.residual.ratio <- function(A, X) {
  AX <- as.matrix(A %*% X)
  denom <- sqrt(sum(X^2))
  if (denom == 0) return(NA_real_)
  sqrt(sum(AX^2)) / denom
}

estimate.nullity <- function(A, rel.tol = c(1e-8, 1e-7, 1e-6, 1e-5)) {
  dense <- as.matrix(A)
  sv <- if (nrow(dense) == 0L || ncol(dense) == 0L) numeric() else svd(dense, nu = 0, nv = 0)$d
  if (length(sv) == 0L) {
    return(list(
      singular.values = sv,
      table = data.frame(rel.tol = rel.tol, rank = 0L,
                         nullity = ncol(A), stringsAsFactors = FALSE)
    ))
  }
  smax <- max(sv)
  tab <- data.frame(
    rel.tol = rel.tol,
    rank = vapply(rel.tol, function(tol) sum(sv > tol * smax), integer(1)),
    stringsAsFactors = FALSE
  )
  tab$nullity <- ncol(A) - tab$rank
  list(singular.values = sv, table = tab)
}

induced.edge.data <- function(adj.list, weight.list, vertices) {
  vertices <- sort(unique(as.integer(vertices)))
  local.index <- match(seq_along(adj.list), vertices)
  in.disk <- rep(FALSE, length(adj.list))
  in.disk[vertices] <- TRUE
  edges <- matrix(NA_integer_, nrow = 0L, ncol = 2L)
  weights <- numeric()
  for (global.i in vertices) {
    nbrs <- adj.list[[global.i]]
    w <- weight.list[[global.i]]
    keep <- in.disk[nbrs] & global.i < nbrs
    if (!any(keep)) next
    local.i <- local.index[global.i]
    local.j <- local.index[nbrs[keep]]
    edges <- rbind(edges, cbind(local.i, local.j))
    weights <- c(weights, w[keep])
  }
  list(edges = edges, edge.weights = as.double(weights), vertices = vertices)
}

edge.stress.diagnostics <- function(coords, edges, edge.weights) {
  if (is.null(coords) || nrow(edges) < 1L) {
    return(data.frame(
      edge.rrmse = NA_real_,
      edge.ratio.q95 = NA_real_,
      edge.ratio.gt.1p10 = NA_real_,
      spread.trace = NA_real_
    ))
  }
  d <- sqrt(rowSums((coords[edges[, 1L], , drop = FALSE] -
                       coords[edges[, 2L], , drop = FALSE])^2))
  ok <- is.finite(d) & is.finite(edge.weights) & edge.weights > 0
  if (!any(ok)) {
    return(data.frame(
      edge.rrmse = NA_real_,
      edge.ratio.q95 = NA_real_,
      edge.ratio.gt.1p10 = NA_real_,
      spread.trace = NA_real_
    ))
  }
  ratio <- d[ok] / edge.weights[ok]
  data.frame(
    edge.rrmse = sqrt(sum((d[ok] - edge.weights[ok])^2) / sum(edge.weights[ok]^2)),
    edge.ratio.q95 = as.numeric(stats::quantile(ratio, 0.95, names = FALSE)),
    edge.ratio.gt.1p10 = mean(ratio > 1.10),
    spread.trace = if (ncol(coords) > 0L && nrow(coords) > 1L) {
      sum(diag(stats::cov(coords)))
    } else {
      0
    }
  )
}

scale.coords.to.edges <- function(coords, edges, edge.weights) {
  if (is.null(coords) || nrow(edges) < 1L) return(coords)
  d <- sqrt(rowSums((coords[edges[, 1L], , drop = FALSE] -
                       coords[edges[, 2L], , drop = FALSE])^2))
  ok <- is.finite(d) & d > 0 & is.finite(edge.weights) & edge.weights > 0
  if (!any(ok)) return(coords)
  scale <- sum(d[ok] * edge.weights[ok]) / sum(d[ok]^2)
  if (is.finite(scale) && scale > 0) coords * scale else coords
}

grip.function.or.null <- function(name) {
  if (!requireNamespace("grip", quietly = TRUE)) return(NULL)
  ns <- asNamespace("grip")
  if (exists(name, envir = ns, inherits = FALSE)) {
    get(name, envir = ns, inherits = FALSE)
  } else {
    NULL
  }
}

local.cmdscale.coords <- function(D.sub, dim) {
  if (nrow(D.sub) <= 1L || any(!is.finite(D.sub))) return(NULL)
  k <- min(as.integer(dim), max(1L, nrow(D.sub) - 1L))
  coords <- tryCatch(
    stats::cmdscale(D.sub, k = k, eig = FALSE, add = FALSE),
    error = function(e) NULL
  )
  if (is.null(coords)) return(NULL)
  coords <- as.matrix(coords)
  if (ncol(coords) < dim) {
    coords <- cbind(coords, matrix(0, nrow(coords), dim - ncol(coords)))
  }
  coords
}

local.edge.kk.function <- function() {
  edge.kk.fn <- grip.function.or.null("grip.optimize.edge.kk.layout")
  if (!is.null(edge.kk.fn)) return(edge.kk.fn)
  grip.function.or.null("grip.optimize.edge.isometric.layout")
}

local.metric.mds.edge.kk.coords <- function(edges,
                                            edge.weights,
                                            n,
                                            dim,
                                            D.sub = NULL,
                                            edge.kk.max.iter = 30L) {
  if (n <= 1L || nrow(edges) < 1L) return(NULL)
  metric.mds.fn <- grip.function.or.null("grip.metric.mds.layout")
  if (!is.null(metric.mds.fn)) {
    fit <- metric.mds.fn(
      edges = edges,
      edge_weights = edge.weights,
      n = n,
      dim = dim,
      diagnostics = FALSE
    )
    coords <- fit$coords
    note <- "metric.mds"
  } else {
    if (is.null(D.sub)) {
      stop("grip.metric.mds.layout() is unavailable and no local distance matrix was supplied.")
    }
    coords <- local.cmdscale.coords(D.sub, dim)
    note <- "cmdscale.metric.mds.fallback"
  }
  coords <- scale.coords.to.edges(as.matrix(coords), edges, edge.weights)
  edge.kk.fn <- local.edge.kk.function()
  if (is.null(edge.kk.fn)) {
    attr(coords, "embedding.note") <- paste(note, "edge.kk.unavailable", sep = ".")
    return(coords)
  }
  fit <- edge.kk.fn(
    coords = coords,
    edges = edges,
    edge_weights = edge.weights,
    n = n,
    dim = dim,
    stiffness_method = "density",
    stiffness_transform = "sqrt",
    density_mix_schedule = c(0, 0.5, 1),
    scale_mode = "identity",
    max_iter = edge.kk.max.iter,
    initial_step = 0.2,
    return_trace = FALSE,
    diagnostics = FALSE,
    engine = "cpp"
  )
  out <- fit$coords
  attr(out, "embedding.note") <- paste(note, "edge.kk", sep = ".")
  out
}

local.weighted.grip.edge.kk.coords <- function(edges,
                                               edge.weights,
                                               n,
                                               dim,
                                               seed,
                                               grip.rounds = 40L,
                                               grip.final.rounds = 80L,
                                               edge.kk.max.iter = 30L) {
  if (!requireNamespace("grip", quietly = TRUE)) {
    stop("Package 'grip' is required for local weighted-GRIP -> edge-KK diagnostics.")
  }
  if (!dim %in% c(2L, 3L)) {
    stop("weighted GRIP local diagnostics currently support dim = 2 or 3")
  }
  if (n <= 1L || nrow(edges) < 1L) return(NULL)
  init <- grip::grip.layout.weighted(
    edges = edges,
    edge_weights = edge.weights,
    n = n,
    dim = dim,
    preset = "irregular",
    rounds = grip.rounds,
    final_rounds = grip.final.rounds,
    seed = seed,
    disconnected = "error"
  )
  init <- scale.coords.to.edges(as.matrix(init), edges, edge.weights)
  edge.kk.fn <- local.edge.kk.function()
  if (is.null(edge.kk.fn)) {
    attr(init, "embedding.note") <- "weighted.grip.only.edge.kk.unavailable"
    return(init)
  }
  fit <- edge.kk.fn(
    coords = init,
    edges = edges,
    edge_weights = edge.weights,
    n = n,
    dim = dim,
    stiffness_method = "density",
    stiffness_transform = "sqrt",
    density_mix_schedule = c(0, 0.5, 1),
    scale_mode = "identity",
    max_iter = edge.kk.max.iter,
    initial_step = 0.2,
    return_trace = FALSE,
    diagnostics = FALSE,
    engine = "cpp"
  )
  coords <- fit$coords
  attr(coords, "embedding.note") <- "weighted.grip.edge.kk"
  coords
}

local.hybrid.edge.kk.coords <- function(edges,
                                        edge.weights,
                                        n,
                                        dim,
                                        seed,
                                        D.sub = NULL) {
  if (dim %in% c(2L, 3L) &&
      !is.null(grip.function.or.null("grip.layout.weighted"))) {
    return(local.weighted.grip.edge.kk.coords(
      edges = edges,
      edge.weights = edge.weights,
      n = n,
      dim = dim,
      seed = seed
    ))
  }
  local.metric.mds.edge.kk.coords(
    edges = edges,
    edge.weights = edge.weights,
    n = n,
    dim = dim,
    D.sub = D.sub
  )
}

normalize.local.dimension.method.label <- function(method) {
  method <- gsub("weighted.grip.edge.kk", "weighted GRIP -> edge-KK", method)
  method <- gsub("metric.mds.edge.kk", "metric MDS -> edge-KK", method)
  method <- gsub("hybrid.edge.kk", "weighted GRIP / metric MDS -> edge-KK", method)
  gsub("cmdscale", "cmdscale", method)
}

select.local.dimension <- function(dim.table,
                                   improvement.threshold = 0.15,
                                   cap.stress.threshold = 0.15) {
  ok <- dim.table[is.finite(dim.table$edge.rrmse), , drop = FALSE]
  if (nrow(ok) < 1L) {
    return(list(selected.dim = NA_integer_, best.dim = NA_integer_,
                dimension.cap.hit = NA, selection.rule = "no.finite.fit"))
  }
  ok <- ok[order(ok$embedding.dim), , drop = FALSE]
  best.dim <- ok$embedding.dim[which.min(ok$edge.rrmse)]
  selected.dim <- ok$embedding.dim[1L]
  if (nrow(ok) > 1L) {
    for (idx in 2:nrow(ok)) {
      previous <- ok$edge.rrmse[idx - 1L]
      current <- ok$edge.rrmse[idx]
      improvement <- if (is.finite(previous) && previous > 0) {
        (previous - current) / previous
      } else {
        0
      }
      if (is.finite(improvement) && improvement >= improvement.threshold) {
        selected.dim <- ok$embedding.dim[idx]
      } else {
        break
      }
    }
  }
  max.dim <- max(ok$embedding.dim)
  final.stress <- ok$edge.rrmse[ok$embedding.dim == max.dim][1L]
  list(
    selected.dim = as.integer(selected.dim),
    best.dim = as.integer(best.dim),
    dimension.cap.hit = isTRUE(best.dim == max.dim && final.stress > cap.stress.threshold),
    selection.rule = sprintf("stress.elbow.relative.%0.2f", improvement.threshold)
  )
}

elbow.dimension.from.variances <- function(values,
                                           dims = 1:5,
                                           improvement.threshold = 0.15) {
  values <- as.double(values)
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 1L) {
    return(list(
      dimension = NA_integer_,
      confidence = NA_real_,
      value = NA_real_,
      rule = "no.positive.spectrum"
    ))
  }
  max.dim <- min(max(dims), length(values))
  dims <- dims[dims <= max.dim]
  selected <- min(dims)
  improvements <- numeric()
  for (q in dims[dims > min(dims)]) {
    remaining.before <- sum(values[q:length(values)])
    improvement <- if (q <= length(values) && remaining.before > 0) {
      values[q] / remaining.before
    } else {
      0
    }
    improvements <- c(improvements, improvement)
    if (is.finite(improvement) && improvement >= improvement.threshold) {
      selected <- q
    } else {
      break
    }
  }
  list(
    dimension = as.integer(selected),
    confidence = if (length(improvements)) max(improvements, na.rm = TRUE) else NA_real_,
    value = if (selected <= length(values)) values[selected] / sum(values) else NA_real_,
    rule = sprintf("spectrum.elbow.relative.%0.2f", improvement.threshold)
  )
}

local.mds.spectrum.vote <- function(D.sub,
                                    dims = 1:5,
                                    variance.threshold = 0.90) {
  out <- tryCatch({
    fit <- stats::cmdscale(
      D.sub,
      k = min(max(dims), max(1L, nrow(D.sub) - 1L)),
      eig = TRUE,
      add = FALSE
    )
    eig <- fit$eig
    eig <- eig[is.finite(eig) & eig > 0]
    if (length(eig) < 1L) {
      stop("No positive metric-MDS eigenvalues.")
    }
    max.dim <- min(max(dims), length(eig))
    dims.use <- dims[dims <= max.dim]
    explained <- cumsum(eig) / sum(eig)
    selected <- dims.use[
      which(explained[dims.use] >= variance.threshold)[1L]
    ]
    if (is.na(selected)) selected <- max(dims.use)
    data.frame(
      dimension.method = "mds.spectrum.elbow",
      estimated.dimension = as.integer(selected),
      confidence = explained[selected],
      diagnostic.value = eig[selected] / sum(eig),
      rule = sprintf("positive.eigen.cumulative.%0.2f", variance.threshold),
      status = "ok",
      message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      dimension.method = "mds.spectrum.elbow",
      estimated.dimension = NA_integer_,
      confidence = NA_real_,
      diagnostic.value = NA_real_,
      rule = "spectrum.elbow",
      status = "error",
      message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  out
}

local.graph.laplacian.vote <- function(edges,
                                       edge.weights,
                                       n,
                                       dims = 1:5) {
  if (n < 3L || nrow(edges) < 1L) {
    return(data.frame(
      dimension.method = "graph.laplacian.weyl",
      estimated.dimension = NA_integer_,
      confidence = NA_real_,
      diagnostic.value = NA_real_,
      rule = "insufficient.local.graph",
      status = "insufficient.local.graph",
      message = "",
      stringsAsFactors = FALSE
    ))
  }
  out <- tryCatch({
    C <- Matrix::sparseMatrix(
      i = c(edges[, 1L], edges[, 2L]),
      j = c(edges[, 2L], edges[, 1L]),
      x = c(edge.weights, edge.weights),
      dims = c(n, n)
    )
    deg <- Matrix::rowSums(C)
    L <- Matrix::Diagonal(x = as.numeric(deg)) - C
    eig <- eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values
    eig <- sort(eig[is.finite(eig) & eig > sqrt(.Machine$double.eps)])
    if (length(eig) < 4L) {
      vote <- elbow.dimension.from.variances(1 / eig, dims)
      return(data.frame(
        dimension.method = "graph.laplacian.weyl",
        estimated.dimension = vote$dimension,
        confidence = vote$confidence,
        diagnostic.value = vote$value,
        rule = "inverse.laplacian.spectrum.fallback",
        status = "ok",
        message = "",
        stringsAsFactors = FALSE
      ))
    }
    k <- min(length(eig), max(8L, min(16L, floor(length(eig) / 2L))))
    lambda <- eig[seq_len(k)]
    counts <- seq_along(lambda)
    keep <- lambda > 0 & counts > 1L
    lambda <- lambda[keep]
    counts <- counts[keep]
    if (length(lambda) < 3L) {
      stop("Too few positive Laplacian eigenvalues for Weyl slope.")
    }
    fit <- stats::lm(log(counts) ~ log(lambda))
    slope <- unname(stats::coef(fit)[2L])
    estimated <- 2 * slope
    dim.est <- as.integer(round(estimated))
    dim.est <- max(min(dims), min(max(dims), dim.est))
    r2 <- summary(fit)$r.squared
    data.frame(
      dimension.method = "graph.laplacian.weyl",
      estimated.dimension = dim.est,
      confidence = r2,
      diagnostic.value = estimated,
      rule = "weyl.logN.loglambda.slope",
      status = "ok",
      message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      dimension.method = "graph.laplacian.weyl",
      estimated.dimension = NA_integer_,
      confidence = NA_real_,
      diagnostic.value = NA_real_,
      rule = "weyl.logN.loglambda.slope",
      status = "error",
      message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  out
}

local.dimension.vote.table <- function(dim.table,
                                       D.sub,
                                       edge.data,
                                       dims = 1:5) {
  selected <- select.local.dimension(dim.table)
  stress.row <- data.frame(
    dimension.method = "edge.stress.elbow",
    estimated.dimension = selected$selected.dim,
    confidence = {
      ok <- dim.table[is.finite(dim.table$edge.rrmse), , drop = FALSE]
      ok <- ok[order(ok$embedding.dim), , drop = FALSE]
      if (nrow(ok) > 1L && selected$selected.dim > min(ok$embedding.dim)) {
        idx <- match(selected$selected.dim, ok$embedding.dim)
        if (!is.na(idx) && idx > 1L && ok$edge.rrmse[idx - 1L] > 0) {
          (ok$edge.rrmse[idx - 1L] - ok$edge.rrmse[idx]) / ok$edge.rrmse[idx - 1L]
        } else {
          NA_real_
        }
      } else {
        NA_real_
      }
    },
    diagnostic.value = {
      ok <- dim.table[is.finite(dim.table$edge.rrmse), , drop = FALSE]
      if (nrow(ok) && !is.na(selected$selected.dim)) {
        ok$edge.rrmse[match(selected$selected.dim, ok$embedding.dim)]
      } else {
        NA_real_
      }
    },
    rule = selected$selection.rule,
    status = if (is.na(selected$selected.dim)) "error" else "ok",
    message = "",
    stringsAsFactors = FALSE
  )
  rbind(
    stress.row,
    local.mds.spectrum.vote(D.sub, dims = dims),
    local.graph.laplacian.vote(
      edges = edge.data$edges,
      edge.weights = edge.data$edge.weights,
      n = length(edge.data$vertices),
      dims = dims
    )
  )
}

consensus.dimension.from.votes <- function(vote.table) {
  ok <- vote.table[
    vote.table$status == "ok" & !is.na(vote.table$estimated.dimension),
    ,
    drop = FALSE
  ]
  if (nrow(ok) < 1L) return(NA_integer_)
  tab <- sort(table(ok$estimated.dimension), decreasing = TRUE)
  as.integer(names(tab)[1L])
}

order.centers.by.rule <- function(adj.list,
                                  center.order = c("natural",
                                                   "increasing.degree",
                                                   "decreasing.degree")) {
  center.order <- match.arg(center.order)
  centers <- seq_along(adj.list)
  degree <- lengths(adj.list)
  if (center.order == "increasing.degree") {
    centers[order(degree, centers)]
  } else if (center.order == "decreasing.degree") {
    centers[order(-degree, centers)]
  } else {
    centers
  }
}

candidate.length <- function(candidate) {
  if (is.null(candidate$x) || length(candidate$x) < 1L) return(Inf)
  max(candidate$x, na.rm = TRUE) - min(candidate$x, na.rm = TRUE)
}

candidate.canonical.key <- function(candidate) {
  if (is.null(candidate$vertices) || length(candidate$vertices) < 1L) return("")
  path.key(candidate$vertices)
}

sort.candidates.deterministic <- function(candidates) {
  if (length(candidates) < 2L) return(candidates)
  tab <- data.frame(
    index = seq_along(candidates),
    source = vapply(candidates, function(z) z$source %||% "", character(1)),
    path.length = vapply(candidates, candidate.length, numeric(1)),
    n.vertices = vapply(candidates, function(z) length(z$vertices %||% integer()), integer(1)),
    key = vapply(candidates, candidate.canonical.key, character(1)),
    stringsAsFactors = FALSE
  )
  candidates[tab$index[order(tab$source, tab$path.length, tab$n.vertices, tab$key, tab$index)]]
}

skeleton.path.length <- function(path, adj.list, weight.list) {
  if (length(path) < 2L) return(0)
  tail(path.cumulative.distances(path, adj.list, weight.list), 1L)
}

finite.diameter.pair <- function(D) {
  upper <- upper.tri(D) & is.finite(D)
  if (!any(upper)) {
    return(c(1L, 1L))
  }
  candidates <- which(upper & D == max(D[upper]), arr.ind = TRUE)
  candidates <- candidates[order(candidates[, 1L], candidates[, 2L]), , drop = FALSE]
  as.integer(candidates[1L, ])
}

coverage.distance.to_skeleton <- function(D, skeleton.vertices) {
  if (length(skeleton.vertices) < 1L) {
    return(rep(Inf, nrow(D)))
  }
  apply(D[, skeleton.vertices, drop = FALSE], 1L, function(z) {
    finite <- z[is.finite(z)]
    if (length(finite) < 1L) Inf else min(finite)
  })
}

add.skeleton.path <- function(paths,
                              path,
                              pass,
                              role,
                              adj.list,
                              weight.list,
                              path.keys) {
  path <- as.integer(path)
  if (length(path) < 2L) {
    return(list(paths = paths, path.keys = path.keys, added = FALSE))
  }
  key <- path.key(path)
  if (key %in% path.keys) {
    return(list(paths = paths, path.keys = path.keys, added = FALSE))
  }
  paths[[length(paths) + 1L]] <- list(
    vertices = path,
    x = path.cumulative.distances(path, adj.list, weight.list),
    source = role,
    pass = as.integer(pass),
    role = role,
    path.length = skeleton.path.length(path, adj.list, weight.list)
  )
  list(paths = paths, path.keys = c(path.keys, key), added = TRUE)
}

path.vertex.membership <- function(paths, n) {
  membership <- vector("list", n)
  for (idx in seq_along(paths)) {
    for (v in unique(paths[[idx]]$vertices)) {
      membership[[v]] <- c(membership[[v]], idx)
    }
  }
  membership
}

rank.from.directions <- function(directions, tol = 1e-6) {
  directions <- as.matrix(directions)
  if (nrow(directions) < 1L || ncol(directions) < 1L) return(0L)
  norms <- sqrt(rowSums(directions^2))
  keep <- is.finite(norms) & norms > sqrt(.Machine$double.eps)
  if (!any(keep)) return(0L)
  directions <- directions[keep, , drop = FALSE] / norms[keep]
  sv <- svd(t(directions), nu = 0, nv = 0)$d
  if (length(sv) < 1L) return(0L)
  as.integer(sum(sv > tol * max(1, max(sv))))
}

path.local.direction <- function(path.vertices,
                                 disk.vertices,
                                 coords,
                                 local.index,
                                 center) {
  hits <- path.vertices[path.vertices %in% disk.vertices]
  if (length(hits) < 2L) return(NULL)
  if (center %in% hits) {
    pos <- which(hits == center)[1L]
    if (pos > 1L && pos < length(hits)) {
      endpoints <- c(hits[pos - 1L], hits[pos + 1L])
    } else if (pos == 1L) {
      endpoints <- c(hits[1L], hits[min(length(hits), 2L)])
    } else {
      endpoints <- c(hits[max(1L, length(hits) - 1L)], hits[length(hits)])
    }
  } else {
    endpoints <- c(hits[1L], tail(hits, 1L))
  }
  idx <- local.index[endpoints]
  if (any(is.na(idx)) || idx[1L] == idx[2L]) return(NULL)
  direction <- coords[idx[2L], , drop = TRUE] - coords[idx[1L], , drop = TRUE]
  if (!any(is.finite(direction)) || sqrt(sum(direction^2)) <= 0) return(NULL)
  as.double(direction)
}

local.skeleton.direction.rank <- function(paths,
                                          disk.vertices,
                                          coords,
                                          center,
                                          rank.tol = 1e-6) {
  local.index <- integer(max(disk.vertices))
  local.index[disk.vertices] <- seq_along(disk.vertices)
  dirs <- lapply(paths, function(path) {
    path.local.direction(path$vertices, disk.vertices, coords, local.index, center)
  })
  dirs <- dirs[!vapply(dirs, is.null, logical(1))]
  if (length(dirs) < 1L) {
    return(list(rank = 0L, n.directions = 0L, directions = matrix(0, nrow = 0L, ncol = ncol(coords))))
  }
  mat <- do.call(rbind, dirs)
  list(rank = rank.from.directions(mat, rank.tol),
       n.directions = nrow(mat),
       directions = mat)
}

select.spread.local.vertices <- function(vertices,
                                         coords,
                                         local.index,
                                         center,
                                         max.vertices = 12L) {
  vertices <- sort(unique(as.integer(vertices)))
  vertices <- vertices[vertices != center]
  max.vertices <- max(1L, as.integer(max.vertices))
  if (length(vertices) <= max.vertices) return(vertices)
  idx <- local.index[vertices]
  keep <- !is.na(idx) & idx > 0L
  vertices <- vertices[keep]
  idx <- idx[keep]
  if (length(vertices) <= max.vertices) return(vertices)
  X <- as.matrix(coords[idx, , drop = FALSE])
  norms <- sqrt(rowSums(X^2))
  selected <- integer()
  first <- order(-norms, vertices)[1L]
  selected <- c(selected, first)
  while (length(selected) < max.vertices) {
    remaining <- setdiff(seq_along(vertices), selected)
    if (length(remaining) < 1L) break
    min.dist <- vapply(remaining, function(i) {
      sqrt.min <- sqrt(rowSums((t(t(X[selected, , drop = FALSE]) - X[i, ]))^2))
      min(sqrt.min, na.rm = TRUE)
    }, numeric(1))
    next.idx <- remaining[order(-min.dist, -norms[remaining], vertices[remaining])][1L]
    selected <- c(selected, next.idx)
  }
  sort(vertices[selected])
}

build.geodesic.skeleton.paths <- function(adj.list,
                                          weight.list,
                                          D,
                                          min.coverage = 2L,
                                          max.paths = NULL,
                                          max.passes = 2L,
                                          geodesic.tol = 1e-8) {
  n <- length(adj.list)
  max.paths <- as.integer(max.paths %||% max(4L, ceiling(n / 2L)))
  min.coverage <- max(1L, as.integer(min.coverage))
  max.passes <- max(1L, as.integer(max.passes))
  paths <- list()
  path.keys <- character()
  diameter <- finite.diameter.pair(D)
  first.path <- deterministic.shortest.path(
    adj.list, weight.list, diameter[1L], diameter[2L], tol = geodesic.tol
  )
  added <- add.skeleton.path(
    paths, first.path, pass = 0L, role = "diameter",
    adj.list = adj.list, weight.list = weight.list, path.keys = path.keys
  )
  paths <- added$paths
  path.keys <- added$path.keys
  coverage <- integer(n)
  if (length(paths) > 0L) {
    coverage[unique(paths[[1L]]$vertices)] <- coverage[unique(paths[[1L]]$vertices)] + 1L
  }
  pass.records <- list()
  pass.records[[1L]] <- data.frame(
    skeleton.path.id = 1L,
    skeleton.pass = 0L,
    skeleton.role = "diameter",
    source = diameter[1L],
    target = diameter[2L],
    attach = NA_integer_,
    n.path.vertices = length(first.path),
    path.length = if (length(paths)) paths[[1L]]$path.length else NA_real_,
    new.vertices = length(unique(first.path)),
    stringsAsFactors = FALSE
  )

  for (pass in seq_len(max.passes)) {
    target.coverage <- min(min.coverage, pass + 1L)
    blocked.targets <- integer()
    guard <- 0L
    repeat {
      guard <- guard + 1L
      if (length(paths) >= max.paths || guard > n * 2L) break
      deficit <- setdiff(which(coverage < target.coverage), blocked.targets)
      if (length(deficit) < 1L) break
      skeleton.vertices <- sort(unique(unlist(lapply(paths, `[[`, "vertices"), use.names = FALSE)))
      d.to.skeleton <- coverage.distance.to_skeleton(D, skeleton.vertices)
      deficit.tab <- data.frame(
        vertex = deficit,
        coverage = coverage[deficit],
        dist.to.skeleton = d.to.skeleton[deficit],
        eccentricity = apply(D[deficit, , drop = FALSE], 1L, function(z) {
          finite <- z[is.finite(z)]
          if (length(finite) < 1L) -Inf else max(finite)
        })
      )
      deficit.tab <- deficit.tab[order(
        deficit.tab$coverage,
        -deficit.tab$dist.to.skeleton,
        -deficit.tab$eccentricity,
        deficit.tab$vertex
      ), , drop = FALSE]
      target <- deficit.tab$vertex[1L]
      finite.to.skeleton <- D[target, skeleton.vertices]
      finite.to.skeleton[!is.finite(finite.to.skeleton)] <- Inf
      if (length(skeleton.vertices) > 0L &&
          any(is.finite(finite.to.skeleton)) &&
          min(finite.to.skeleton) > geodesic.tol) {
        attach <- skeleton.vertices[order(finite.to.skeleton, skeleton.vertices)][1L]
        path <- deterministic.shortest.path(adj.list, weight.list, attach, target,
                                            tol = geodesic.tol)
        role <- "farthest.to.skeleton"
      } else {
        distances <- D[target, ]
        distances[!is.finite(distances)] <- -Inf
        distances[target] <- -Inf
        endpoint <- which(distances == max(distances))[1L]
        attach <- target
        path <- deterministic.shortest.path(adj.list, weight.list, target, endpoint,
                                            tol = geodesic.tol)
        role <- "coverage.pass"
      }
      old.coverage <- coverage
      added <- add.skeleton.path(
        paths, path, pass = pass, role = role,
        adj.list = adj.list, weight.list = weight.list, path.keys = path.keys
      )
      if (!added$added) {
        blocked.targets <- union(blocked.targets, target)
        next
      }
      paths <- added$paths
      path.keys <- added$path.keys
      path.vertices <- unique(path)
      coverage[path.vertices] <- coverage[path.vertices] + 1L
      pass.records[[length(pass.records) + 1L]] <- data.frame(
        skeleton.path.id = length(paths),
        skeleton.pass = pass,
        skeleton.role = role,
        source = path[1L],
        target = tail(path, 1L),
        attach = attach,
        n.path.vertices = length(path),
        path.length = tail(path.cumulative.distances(path, adj.list, weight.list), 1L),
        new.vertices = sum(old.coverage[path.vertices] == 0L),
        stringsAsFactors = FALSE
      )
    }
  }
  list(paths = paths, coverage = coverage,
       pass.table = do.call(rbind, pass.records))
}

add.local.cross.skeleton.paths <- function(skeleton,
                                           adj.list,
                                           weight.list,
                                           D,
                                           disk.k = 12L,
                                           max.cross.paths = 30L,
                                           min.coverage = 2L,
                                           geodesic.tol = 1e-8) {
  paths <- skeleton$paths
  if (length(paths) < 2L || max.cross.paths < 1L) return(skeleton)
  n <- length(adj.list)
  path.keys <- vapply(paths, function(p) path.key(p$vertices), character(1))
  coverage <- skeleton$coverage
  membership <- path.vertex.membership(paths, n)
  centers <- seq_len(n)
  centers <- centers[order(coverage, -lengths(adj.list), centers)]
  cross.records <- list()
  for (center in centers) {
    if (length(cross.records) >= max.cross.paths) break
    finite.dist <- D[center, is.finite(D[center, ]) & D[center, ] > 0]
    if (length(finite.dist) < 2L) next
    radius <- sort(finite.dist)[min(as.integer(disk.k), length(finite.dist))]
    disk.vertices <- which(is.finite(D[center, ]) & D[center, ] <= radius)
    disk.vertices <- disk.vertices[lengths(membership[disk.vertices]) > 0L]
    if (length(disk.vertices) < 2L) next
    best <- NULL
    for (a in seq_len(length(disk.vertices) - 1L)) {
      va <- disk.vertices[a]
      for (b in seq.int(a + 1L, length(disk.vertices))) {
        vb <- disk.vertices[b]
        if (!is.finite(D[va, vb]) || D[va, vb] <= geodesic.tol) next
        if (length(intersect(membership[[va]], membership[[vb]])) > 0L) next
        score <- D[va, center] + D[vb, center] + D[va, vb]
        candidate <- data.frame(va = va, vb = vb, score = score)
        if (is.null(best) ||
            candidate$score > best$score ||
            (candidate$score == best$score && paste(va, vb) < paste(best$va, best$vb))) {
          best <- candidate
        }
      }
    }
    if (is.null(best)) next
    path <- deterministic.shortest.path(adj.list, weight.list, best$va, best$vb,
                                        tol = geodesic.tol)
    added <- add.skeleton.path(
      paths, path, pass = max(skeleton$pass.table$skeleton.pass) + 1L,
      role = "local.cross.connector", adj.list = adj.list,
      weight.list = weight.list, path.keys = path.keys
    )
    if (!added$added) next
    paths <- added$paths
    path.keys <- added$path.keys
    path.vertices <- unique(path)
    old.coverage <- coverage
    coverage[path.vertices] <- coverage[path.vertices] + 1L
    cross.records[[length(cross.records) + 1L]] <- data.frame(
      skeleton.path.id = length(paths),
      skeleton.pass = max(skeleton$pass.table$skeleton.pass) + 1L,
      skeleton.role = "local.cross.connector",
      source = path[1L],
      target = tail(path, 1L),
      attach = center,
      n.path.vertices = length(path),
      path.length = tail(path.cumulative.distances(path, adj.list, weight.list), 1L),
      new.vertices = sum(old.coverage[path.vertices] == 0L),
      stringsAsFactors = FALSE
    )
    membership <- path.vertex.membership(paths, n)
  }
  skeleton$paths <- paths
  skeleton$coverage <- coverage
  if (length(cross.records) > 0L) {
    skeleton$pass.table <- rbind(skeleton$pass.table, do.call(rbind, cross.records))
  }
  skeleton$direction.table <- data.frame()
  skeleton$connector.candidate.table <- data.frame()
  skeleton$cross.rule <- "coverage"
  skeleton
}

add.dimension.rank.cross.skeleton.paths <- function(skeleton,
                                                    adj.list,
                                                    weight.list,
                                                    D,
                                                    disk.k = 12L,
                                                    max.cross.paths = 30L,
                                                    local.dims = 1:5,
                                                    embedding.method = c("cmdscale", "hybrid.edge.kk"),
                                                    candidate.pool = c("skeleton.membership",
                                                                       "rich.local"),
                                                    rank.tol = 1e-6,
                                                    max.candidates.per.center = 50L,
                                                    length.penalty = 0.05,
                                                    coverage.reward = 0.10,
                                                    iterations = 1L,
                                                    geodesic.tol = 1e-8,
                                                    seed = 1L) {
  embedding.method <- match.arg(embedding.method)
  candidate.pool <- match.arg(candidate.pool)
  paths <- skeleton$paths
  if (length(paths) < 2L || max.cross.paths < 1L) {
    skeleton$direction.table <- data.frame()
    skeleton$connector.candidate.table <- data.frame()
    skeleton$cross.rule <- "local.dimension.rank"
    return(skeleton)
  }
  n <- length(adj.list)
  path.keys <- vapply(paths, function(p) path.key(p$vertices), character(1))
  coverage <- skeleton$coverage
  if (!"skeleton.cross.iteration" %in% names(skeleton$pass.table)) {
    skeleton$pass.table$skeleton.cross.iteration <- NA_integer_
  }
  direction.records <- list()
  candidate.records <- list()
  added.cross <- 0L
  iterations <- max(1L, as.integer(iterations))
  for (iteration in seq_len(iterations)) {
    centers <- seq_len(n)
    centers <- centers[order(coverage, -lengths(adj.list), centers)]
    for (center in centers) {
    if (added.cross >= max.cross.paths) break
    finite.dist <- D[center, is.finite(D[center, ]) & D[center, ] > 0]
    if (length(finite.dist) < 2L) next
    radius <- sort(finite.dist)[min(as.integer(disk.k), length(finite.dist))]
    disk.vertices <- which(is.finite(D[center, ]) & D[center, ] <= radius)
    disk.vertices <- sort(unique(as.integer(disk.vertices)))
    if (length(disk.vertices) < 3L) next
    D.sub <- D[disk.vertices, disk.vertices, drop = FALSE]
    vote <- local.mds.spectrum.vote(D.sub, dims = local.dims)
    local.dim <- vote$estimated.dimension[1L]
    if (!is.finite(local.dim) || is.na(local.dim) || local.dim < 1L) local.dim <- 1L
    local.dim <- min(as.integer(local.dim), max(local.dims), max(1L, length(disk.vertices) - 1L))
    embed.dim <- min(max(local.dims), max(local.dim, 2L), max(1L, length(disk.vertices) - 1L))
    edge.data <- induced.edge.data(adj.list, weight.list, disk.vertices)
    coords <- tryCatch({
      if (embedding.method == "hybrid.edge.kk") {
        local.hybrid.edge.kk.coords(
          edges = edge.data$edges,
          edge.weights = edge.data$edge.weights,
          n = length(disk.vertices),
          dim = embed.dim,
          seed = seed + 1000L * center,
          D.sub = D.sub
        )
      } else {
        local.cmdscale.coords(D.sub, embed.dim)
      }
    }, error = function(e) NULL)
    if (is.null(coords)) next
    coords <- as.matrix(coords)
    if (ncol(coords) < embed.dim) {
      coords <- cbind(coords, matrix(0, nrow(coords), embed.dim - ncol(coords)))
    }
    membership <- path.vertex.membership(paths, n)
    local.path.ids <- sort(unique(unlist(membership[disk.vertices], use.names = FALSE)))
    local.path.ids <- local.path.ids[!is.na(local.path.ids)]
    if (length(local.path.ids) < 2L) {
      direction.records[[length(direction.records) + 1L]] <- data.frame(
        center = center,
        skeleton.cross.iteration = iteration,
        disk.radius = radius,
        n.disk.vertices = length(disk.vertices),
        local.dimension = local.dim,
        embedding.dim = embed.dim,
        direction.rank.before = 0L,
        direction.rank.after = 0L,
        direction.deficit.before = local.dim,
        direction.deficit.after = local.dim,
        n.local.skeleton.paths = length(local.path.ids),
        n.local.connector.candidates = 0L,
        connector.added = FALSE,
        connector.path.id = NA_integer_,
        connector.path = NA_character_,
        connector.length = NA_real_,
        rank.gain = 0L,
        coverage.gain = 0L,
        score = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    rank.before <- local.skeleton.direction.rank(
      paths[local.path.ids], disk.vertices, coords, center, rank.tol
    )$rank
    deficit.before <- max(0L, local.dim - rank.before)
    if (deficit.before < 1L) {
      direction.records[[length(direction.records) + 1L]] <- data.frame(
        center = center,
        skeleton.cross.iteration = iteration,
        disk.radius = radius,
        n.disk.vertices = length(disk.vertices),
        local.dimension = local.dim,
        embedding.dim = embed.dim,
        direction.rank.before = rank.before,
        direction.rank.after = rank.before,
        direction.deficit.before = deficit.before,
        direction.deficit.after = deficit.before,
        n.local.skeleton.paths = length(local.path.ids),
        n.local.connector.candidates = 0L,
        connector.added = FALSE,
        connector.path.id = NA_integer_,
        connector.path = NA_character_,
        connector.length = NA_real_,
        rank.gain = 0L,
        coverage.gain = 0L,
        score = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    candidate.vertices <- disk.vertices[lengths(membership[disk.vertices]) > 0L]
    candidate.list <- list()
    candidate.keys <- character()
    local.index <- integer(max(disk.vertices))
    local.index[disk.vertices] <- seq_along(disk.vertices)
    current.dirs <- local.skeleton.direction.rank(
      paths[local.path.ids], disk.vertices, coords, center, rank.tol
    )$directions
    add.candidate <- function(va, vb, candidate.type) {
      if (!is.finite(D[va, vb]) || D[va, vb] <= geodesic.tol) return(NULL)
      path <- deterministic.shortest.path(adj.list, weight.list, va, vb,
                                          tol = geodesic.tol)
      if (length(path) < 2L) return(NULL)
      key <- path.key(path)
      if (key %in% path.keys || key %in% candidate.keys) return(NULL)
      path.local <- path[path %in% disk.vertices]
      if (length(path.local) < 2L) return(NULL)
      new.direction <- path.local.direction(path, disk.vertices, coords,
                                            local.index, center)
      if (is.null(new.direction)) return(NULL)
      rank.after <- rank.from.directions(rbind(current.dirs, new.direction), rank.tol)
      rank.gain <- max(0L, rank.after - rank.before)
      path.length <- skeleton.path.length(path, adj.list, weight.list)
      coverage.gain <- sum(coverage[unique(path)] == 0L) +
        0.25 * sum(coverage[unique(path)] < local.dim)
      score <- rank.gain + coverage.reward * coverage.gain -
        length.penalty * path.length
      candidate.keys <<- c(candidate.keys, key)
      candidate.list[[length(candidate.list) + 1L]] <<- list(
        path = path,
        key = key,
        va = va,
        vb = vb,
        candidate.type = candidate.type,
        path.length = path.length,
        rank.after = rank.after,
        rank.gain = rank.gain,
        coverage.gain = coverage.gain,
        score = score
      )
      NULL
    }
    if (length(candidate.vertices) >= 2L) {
      for (a in seq_len(length(candidate.vertices) - 1L)) {
        va <- candidate.vertices[a]
        for (b in seq.int(a + 1L, length(candidate.vertices))) {
          vb <- candidate.vertices[b]
          if (length(intersect(membership[[va]], membership[[vb]])) > 0L) next
          add.candidate(va, vb, "skeleton.membership")
        }
      }
    }
    if (candidate.pool == "rich.local") {
      radial.dist <- D[center, disk.vertices]
      boundary.vertices <- disk.vertices[
        is.finite(radial.dist) & radial.dist >= 0.70 * radius
      ]
      annulus.vertices <- disk.vertices[
        is.finite(radial.dist) & radial.dist >= 0.35 * radius &
          radial.dist < 0.70 * radius
      ]
      boundary.vertices <- select.spread.local.vertices(
        boundary.vertices, coords, local.index, center, max.vertices = 12L
      )
      annulus.vertices <- select.spread.local.vertices(
        annulus.vertices, coords, local.index, center, max.vertices = 12L
      )
      if (length(boundary.vertices) > 0L && center %in% disk.vertices) {
        for (vb in boundary.vertices) {
          add.candidate(center, vb, "radial.center.to.boundary")
        }
      }
      if (length(boundary.vertices) >= 2L) {
        for (a in seq_len(length(boundary.vertices) - 1L)) {
          va <- boundary.vertices[a]
          for (b in seq.int(a + 1L, length(boundary.vertices))) {
            vb <- boundary.vertices[b]
            add.candidate(va, vb, "boundary.to.boundary")
          }
        }
      }
      if (length(annulus.vertices) > 0L && length(boundary.vertices) > 0L) {
        for (va in annulus.vertices) {
          for (vb in boundary.vertices) {
            if (va == vb) next
            add.candidate(va, vb, "annulus.to.boundary")
          }
        }
      }
    }
    if (length(candidate.list) > max.candidates.per.center) {
      tab <- data.frame(
        idx = seq_along(candidate.list),
        rank.gain = vapply(candidate.list, `[[`, numeric(1), "rank.gain"),
        coverage.gain = vapply(candidate.list, `[[`, numeric(1), "coverage.gain"),
        score = vapply(candidate.list, `[[`, numeric(1), "score"),
        path.length = vapply(candidate.list, `[[`, numeric(1), "path.length"),
        key = vapply(candidate.list, `[[`, character(1), "key"),
        stringsAsFactors = FALSE
      )
      keep <- head(tab$idx[order(-tab$rank.gain, -tab$coverage.gain, -tab$score,
                                 tab$path.length, tab$key)], max.candidates.per.center)
      candidate.list <- candidate.list[keep]
    }
    n.candidates <- length(candidate.list)
    candidate.table <- if (n.candidates > 0L) {
      do.call(rbind, lapply(seq_along(candidate.list), function(idx) {
        z <- candidate.list[[idx]]
        data.frame(
          center = center,
          skeleton.cross.iteration = iteration,
          candidate.index = idx,
          source = z$va,
          target = z$vb,
          candidate.type = z$candidate.type %||% "skeleton.membership",
          path = paste(z$path, collapse = "-"),
          path.length = z$path.length,
          direction.rank.before = rank.before,
          direction.rank.after = z$rank.after,
          rank.gain = z$rank.gain,
          coverage.gain = z$coverage.gain,
          score = z$score,
          selected = FALSE,
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame()
    }
    selected <- NULL
    if (n.candidates > 0L) {
      tab <- data.frame(
        idx = seq_along(candidate.list),
        rank.gain = vapply(candidate.list, `[[`, numeric(1), "rank.gain"),
        coverage.gain = vapply(candidate.list, `[[`, numeric(1), "coverage.gain"),
        score = vapply(candidate.list, `[[`, numeric(1), "score"),
        path.length = vapply(candidate.list, `[[`, numeric(1), "path.length"),
        key = vapply(candidate.list, `[[`, character(1), "key"),
        stringsAsFactors = FALSE
      )
      tab <- tab[order(-tab$rank.gain, -tab$coverage.gain, -tab$score,
                       tab$path.length, tab$key), , drop = FALSE]
      if (tab$rank.gain[1L] > 0L || tab$coverage.gain[1L] > 0L) {
        selected <- candidate.list[[tab$idx[1L]]]
        if (nrow(candidate.table) > 0L) {
          candidate.table$selected[candidate.table$candidate.index == tab$idx[1L]] <- TRUE
        }
      }
    }
    if (nrow(candidate.table) > 0L) {
      candidate.records[[length(candidate.records) + 1L]] <- candidate.table
    }
    if (is.null(selected)) {
      direction.records[[length(direction.records) + 1L]] <- data.frame(
       center = center,
        skeleton.cross.iteration = iteration,
        disk.radius = radius,
        n.disk.vertices = length(disk.vertices),
        local.dimension = local.dim,
        embedding.dim = embed.dim,
        direction.rank.before = rank.before,
        direction.rank.after = rank.before,
        direction.deficit.before = deficit.before,
        direction.deficit.after = deficit.before,
        n.local.skeleton.paths = length(local.path.ids),
        n.local.connector.candidates = n.candidates,
        connector.added = FALSE,
        connector.path.id = NA_integer_,
        connector.path = NA_character_,
        connector.length = NA_real_,
        rank.gain = 0L,
        coverage.gain = 0L,
        score = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    old.coverage <- coverage
    added <- add.skeleton.path(
      paths, selected$path,
      pass = max(skeleton$pass.table$skeleton.pass) + 1L,
      role = "local.dimension.rank.connector",
      adj.list = adj.list,
      weight.list = weight.list,
      path.keys = path.keys
    )
    if (!added$added) next
    paths <- added$paths
    path.keys <- added$path.keys
    path.vertices <- unique(selected$path)
    coverage[path.vertices] <- coverage[path.vertices] + 1L
    added.cross <- added.cross + 1L
    rank.after <- selected$rank.after
    deficit.after <- max(0L, local.dim - rank.after)
    skeleton$pass.table <- rbind(skeleton$pass.table, data.frame(
      skeleton.path.id = length(paths),
      skeleton.pass = max(skeleton$pass.table$skeleton.pass) + 1L,
      skeleton.cross.iteration = iteration,
      skeleton.role = "local.dimension.rank.connector",
      source = selected$path[1L],
      target = tail(selected$path, 1L),
      attach = center,
      n.path.vertices = length(selected$path),
      path.length = selected$path.length,
      new.vertices = sum(old.coverage[path.vertices] == 0L),
      stringsAsFactors = FALSE
    ))
    direction.records[[length(direction.records) + 1L]] <- data.frame(
      center = center,
      skeleton.cross.iteration = iteration,
      disk.radius = radius,
      n.disk.vertices = length(disk.vertices),
      local.dimension = local.dim,
      embedding.dim = embed.dim,
      direction.rank.before = rank.before,
      direction.rank.after = rank.after,
      direction.deficit.before = deficit.before,
      direction.deficit.after = deficit.after,
      n.local.skeleton.paths = length(local.path.ids),
      n.local.connector.candidates = n.candidates,
      connector.added = TRUE,
      connector.path.id = length(paths),
      connector.path = paste(selected$path, collapse = "-"),
      connector.length = selected$path.length,
      rank.gain = selected$rank.gain,
      coverage.gain = selected$coverage.gain,
      score = selected$score,
      stringsAsFactors = FALSE
    )
    }
  }
  skeleton$paths <- paths
  skeleton$coverage <- coverage
  skeleton$direction.table <- if (length(direction.records) > 0L) {
    do.call(rbind, direction.records)
  } else {
    data.frame()
  }
  skeleton$connector.candidate.table <- if (length(candidate.records) > 0L) {
    do.call(rbind, candidate.records)
  } else {
    data.frame()
  }
  skeleton$cross.rule <- "local.dimension.rank"
  skeleton
}

build.skeleton.row.objects <- function(paths,
                                       adj.list,
                                       weight.list,
                                       n,
                                       degree,
                                       row.normalization) {
  row.objects <- list()
  window.size <- degree + 2L
  for (path.id in seq_along(paths)) {
    vertices <- paths[[path.id]]$vertices
    if (length(vertices) < window.size) next
    for (window.start in seq_len(length(vertices) - window.size + 1L)) {
      window.vertices <- vertices[window.start:(window.start + window.size - 1L)]
      window.x <- path.cumulative.distances(window.vertices, adj.list, weight.list)
      rows <- annihilator.rows.for.path(
        window.vertices, window.x, n.vertices = n,
        degree = degree, row.normalization = row.normalization
      )
      if (length(rows) < 1L) next
      for (basis.idx in seq_along(rows)) {
        row.objects[[length(row.objects) + 1L]] <- list(
          j = rows[[basis.idx]]$j,
          x = rows[[basis.idx]]$x,
          candidate.idx = path.id,
          basis.idx = basis.idx,
          source = paste0("geodesic.skeleton.", paths[[path.id]]$role),
          n.path.vertices = length(window.vertices),
          path.length = tail(window.x, 1L),
          path.vertices = paste(window.vertices, collapse = "-"),
          skeleton.path.id = path.id,
          skeleton.pass = paths[[path.id]]$pass,
          skeleton.role = paths[[path.id]]$role,
          window.start = window.start
        )
      }
    }
  }
  row.objects
}

build.geodesic.skeleton.annihilator.operator <- function(adj.list,
                                                         weight.list,
                                                         degree = 2L,
                                                         path.family = c("geodesic.skeleton",
                                                                         "geodesic.skeleton.cross.enriched"),
                                                         disk.k = 12L,
                                                         max.skeleton.paths = NULL,
                                                         skeleton.min.coverage = 2L,
                                                         skeleton.max.passes = 2L,
                                                         skeleton.max.cross.paths = 30L,
                                                         skeleton.cross.rule = c("coverage",
                                                                                 "local.dimension.rank"),
                                                         skeleton.cross.embedding.method = c("cmdscale",
                                                                                             "hybrid.edge.kk"),
                                                         skeleton.cross.candidate.pool = c("skeleton.membership",
                                                                                           "rich.local"),
                                                         skeleton.cross.local.dims = 1:5,
                                                         skeleton.cross.rank.tol = 1e-6,
                                                         skeleton.cross.max.candidates.per.center = 50L,
                                                         skeleton.cross.length.penalty = 0.05,
                                                         skeleton.cross.coverage.reward = 0.10,
                                                         skeleton.cross.iterations = 1L,
                                                         row.normalization = c("unit.l2", "none"),
                                                         vertex.normalization = c("mean", "none"),
                                                         seed = 1L,
                                                         geodesic.tol = 1e-8) {
  path.family <- match.arg(path.family)
  skeleton.cross.rule <- match.arg(skeleton.cross.rule)
  skeleton.cross.embedding.method <- match.arg(skeleton.cross.embedding.method)
  skeleton.cross.candidate.pool <- match.arg(skeleton.cross.candidate.pool)
  row.normalization <- match.arg(row.normalization)
  vertex.normalization <- match.arg(vertex.normalization)
  n <- length(adj.list)
  D <- shortest.path(adj.list, weight.list, seq_len(n))
  skeleton <- build.geodesic.skeleton.paths(
    adj.list = adj.list,
    weight.list = weight.list,
    D = D,
    min.coverage = skeleton.min.coverage,
    max.paths = max.skeleton.paths,
    max.passes = skeleton.max.passes,
    geodesic.tol = geodesic.tol
  )
  if (path.family == "geodesic.skeleton.cross.enriched") {
    if (skeleton.cross.rule == "local.dimension.rank") {
      skeleton <- add.dimension.rank.cross.skeleton.paths(
        skeleton = skeleton,
        adj.list = adj.list,
        weight.list = weight.list,
        D = D,
        disk.k = disk.k,
        max.cross.paths = skeleton.max.cross.paths,
        local.dims = skeleton.cross.local.dims,
        embedding.method = skeleton.cross.embedding.method,
        candidate.pool = skeleton.cross.candidate.pool,
        rank.tol = skeleton.cross.rank.tol,
        max.candidates.per.center = skeleton.cross.max.candidates.per.center,
        length.penalty = skeleton.cross.length.penalty,
        coverage.reward = skeleton.cross.coverage.reward,
        iterations = skeleton.cross.iterations,
        geodesic.tol = geodesic.tol,
        seed = seed
      )
    } else {
      skeleton <- add.local.cross.skeleton.paths(
        skeleton = skeleton,
        adj.list = adj.list,
        weight.list = weight.list,
        D = D,
        disk.k = disk.k,
        max.cross.paths = skeleton.max.cross.paths,
        min.coverage = skeleton.min.coverage,
        geodesic.tol = geodesic.tol
      )
    }
  } else {
    skeleton$direction.table <- data.frame()
    skeleton$connector.candidate.table <- data.frame()
    skeleton$cross.rule <- "none"
  }
  row.objects <- build.skeleton.row.objects(
    paths = skeleton$paths,
    adj.list = adj.list,
    weight.list = weight.list,
    n = n,
    degree = degree,
    row.normalization = row.normalization
  )
  i.idx <- integer()
  j.idx <- integer()
  x.val <- numeric()
  path.records <- list()
  if (length(row.objects) > 0L) {
    for (row in seq_along(row.objects)) {
      row.obj <- row.objects[[row]]
      row.x <- row.obj$x
      if (vertex.normalization == "mean") {
        row.x <- row.x / length(row.objects)
      }
      i.idx <- c(i.idx, rep(row, length(row.obj$j)))
      j.idx <- c(j.idx, row.obj$j)
      x.val <- c(x.val, row.x)
      path.records[[row]] <- data.frame(
        row = row,
        center = NA_integer_,
        candidate = row.obj$candidate.idx,
        basis.row = row.obj$basis.idx,
        degree = degree,
        source = row.obj$source,
        n.path.vertices = row.obj$n.path.vertices,
        path.length = row.obj$path.length,
        path.vertices = row.obj$path.vertices,
        row.selection = "none",
        center.order = "skeleton",
        center.processing.rank = NA_integer_,
        candidate.cap.mode = "deterministic",
        path.family = path.family,
        skeleton.cross.rule = skeleton$cross.rule %||% "none",
        skeleton.path.id = row.obj$skeleton.path.id,
        skeleton.pass = row.obj$skeleton.pass,
        skeleton.role = row.obj$skeleton.role,
        window.start = row.obj$window.start,
        polynomial.residual = NA_real_,
        selection.rank = row,
        stringsAsFactors = FALSE
      )
    }
  }
  A <- Matrix::sparseMatrix(
    i = i.idx, j = j.idx, x = x.val, dims = c(length(row.objects), n)
  )
  row.coverage <- tabulate(unlist(lapply(row.objects, `[[`, "j"), use.names = FALSE),
                           nbins = n)
  skeleton.vertices <- sort(unique(unlist(lapply(skeleton$paths, `[[`, "vertices"),
                                          use.names = FALSE)))
  d.to.skeleton <- coverage.distance.to_skeleton(D, skeleton.vertices)
  vertex.table <- data.frame(
    center = seq_len(n),
    center.order = "skeleton",
    center.processing.rank = seq_len(n),
    degree = degree,
    disk.radius.rule = "geodesic.skeleton",
    disk.k = disk.k,
    disk.k.start = NA_integer_,
    disk.k.step = NA_integer_,
    disk.k.max = NA_integer_,
    selected.disk.k = NA_integer_,
    attempted.disk.k = NA_character_,
    radius.status = "geodesic.skeleton",
    disk.radius = d.to.skeleton,
    n.disk.vertices = NA_integer_,
    n.rays = NA_integer_,
    n.composite.available = 0L,
    n.paths.used = skeleton$coverage,
    n.operator.rows = row.coverage,
    row.selection = "none",
    candidate.cap.mode = "deterministic",
    path.family = path.family,
    skeleton.cross.rule = skeleton$cross.rule %||% "none",
    skeleton.coverage = skeleton$coverage,
    skeleton.distance = d.to.skeleton,
    candidate.operator.rows = length(row.objects),
    selected.operator.rows = length(row.objects),
    thinning.local.dim = NA_integer_,
    thinning.polynomial.dim = NA_integer_,
    thinning.rank = NA_integer_,
    thinning.coverage = sum(skeleton$coverage > 0L),
    thinning.residual.min = NA_real_,
    thinning.residual.median = NA_real_,
    thinning.residual.max = NA_real_,
    used.mirror.proxy = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    A = A,
    path.table = if (length(path.records)) do.call(rbind, path.records) else data.frame(),
    vertex.table = vertex.table,
    control = list(
      degree = degree,
      path.family = path.family,
      disk.k = disk.k,
      skeleton.min.coverage = skeleton.min.coverage,
      skeleton.max.passes = skeleton.max.passes,
      skeleton.max.cross.paths = skeleton.max.cross.paths,
      skeleton.cross.rule = skeleton$cross.rule %||% "none",
      skeleton.cross.embedding.method = skeleton.cross.embedding.method,
      skeleton.cross.candidate.pool = skeleton.cross.candidate.pool,
      skeleton.cross.local.dims = paste(skeleton.cross.local.dims, collapse = ","),
      skeleton.cross.rank.tol = skeleton.cross.rank.tol,
      skeleton.cross.max.candidates.per.center = skeleton.cross.max.candidates.per.center,
      skeleton.cross.length.penalty = skeleton.cross.length.penalty,
      skeleton.cross.coverage.reward = skeleton.cross.coverage.reward,
      skeleton.cross.iterations = skeleton.cross.iterations,
      max.skeleton.paths = max.skeleton.paths,
      row.normalization = row.normalization,
      vertex.normalization = vertex.normalization
    ),
    local.dimension.table = data.frame(),
    local.dimension.vote.table = data.frame(),
    skeleton.table = skeleton$pass.table,
    skeleton.direction.table = skeleton$direction.table %||% data.frame(),
    skeleton.connector.candidate.table = skeleton$connector.candidate.table %||% data.frame(),
    skeleton.coverage = data.frame(
      vertex = seq_len(n),
      skeleton.coverage = skeleton$coverage,
      skeleton.distance = d.to.skeleton,
      row.coverage = row.coverage,
      stringsAsFactors = FALSE
    )
  )
}

row.object.local.vector <- function(row.object, disk.vertices) {
  idx <- match(row.object$j, disk.vertices)
  keep <- !is.na(idx)
  out <- numeric(length(disk.vertices))
  if (any(keep)) {
    out[idx[keep]] <- row.object$x[keep]
  }
  out
}

score.row.polynomial.residual <- function(row.object, disk.vertices, P) {
  idx <- match(row.object$j, disk.vertices)
  keep <- !is.na(idx)
  if (!any(keep)) return(Inf)
  value <- drop(row.object$x[keep] %*% P[idx[keep], , drop = FALSE])
  denom <- sqrt(sum(P[idx[keep], , drop = FALSE]^2))
  if (!is.finite(denom) || denom <= 0) denom <- 1
  sqrt(sum(value^2)) / denom
}

thin.center.rows.polynomial.residual <- function(row.objects,
                                                 disk.vertices,
                                                 D.sub,
                                                 degree,
                                                 dims = 1:5,
                                                 max.rows = NULL,
                                                 min.rows = 1L) {
  n.candidates <- length(row.objects)
  if (n.candidates < 1L) {
    return(list(
      rows = row.objects,
      diagnostics = data.frame(
        row.selection = "polynomial.residual.greedy",
        candidate.rows = 0L,
        selected.rows = 0L,
        thinning.local.dim = NA_integer_,
        thinning.polynomial.dim = NA_integer_,
        thinning.rank = 0L,
        thinning.coverage = 0L,
        residual.min = NA_real_,
        residual.median = NA_real_,
        residual.max = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }
  vote <- local.mds.spectrum.vote(D.sub, dims = dims)
  local.dim <- vote$estimated.dimension[1L]
  if (is.na(local.dim) || local.dim < 1L) local.dim <- 1L
  local.dim <- min(local.dim, max(1L, nrow(D.sub) - 1L), max(dims))
  coords <- local.cmdscale.coords(D.sub, local.dim)
  if (is.null(coords)) {
    return(list(
      rows = row.objects,
      diagnostics = data.frame(
        row.selection = "polynomial.residual.greedy.fallback.all",
        candidate.rows = n.candidates,
        selected.rows = n.candidates,
        thinning.local.dim = local.dim,
        thinning.polynomial.dim = NA_integer_,
        thinning.rank = NA_integer_,
        thinning.coverage = length(unique(unlist(lapply(row.objects, `[[`, "j"), use.names = FALSE))),
        residual.min = NA_real_,
        residual.median = NA_real_,
        residual.max = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }
  P <- monomial.design.nd(coords, degree)
  p.dim <- ncol(P)
  target.rows <- max(as.integer(min.rows), length(disk.vertices) - p.dim)
  if (!is.null(max.rows) && is.finite(max.rows)) {
    target.rows <- min(target.rows, as.integer(max.rows))
  }
  target.rows <- max(1L, min(n.candidates, target.rows))
  residuals <- vapply(
    row.objects,
    score.row.polynomial.residual,
    numeric(1),
    disk.vertices = disk.vertices,
    P = P
  )
  finite.res <- residuals[is.finite(residuals)]
  residual.threshold <- if (length(finite.res) > 0L) {
    stats::quantile(finite.res, 0.75, names = FALSE, na.rm = TRUE)
  } else {
    Inf
  }
  path.lengths <- vapply(row.objects, function(z) z$path.length %||% Inf, numeric(1))
  n.path.vertices <- vapply(row.objects, function(z) z$n.path.vertices %||% 0L, integer(1))
  path.keys <- vapply(row.objects, function(z) z$path.vertices %||% "", character(1))
  basis.idx <- vapply(row.objects, function(z) z$basis.idx %||% 0L, integer(1))
  order.idx <- order(residuals, path.lengths, n.path.vertices, path.keys, basis.idx,
                     na.last = TRUE)
  selected <- integer()
  covered <- integer()
  selected.mat <- matrix(0, nrow = 0L, ncol = length(disk.vertices))
  current.rank <- 0L
  for (idx in order.idx) {
    row.vec <- row.object.local.vector(row.objects[[idx]], disk.vertices)
    if (!any(row.vec != 0)) next
    trial.mat <- rbind(selected.mat, row.vec)
    trial.rank <- qr(trial.mat, LAPACK = TRUE)$rank
    new.coverage <- setdiff(row.objects[[idx]]$j, covered)
    keeps.rank <- trial.rank > current.rank
    keeps.coverage <- length(new.coverage) > 0L && length(selected) < target.rows
    keeps.residual <- is.finite(residuals[idx]) && residuals[idx] <= residual.threshold
    if (keeps.residual && (keeps.rank || keeps.coverage || length(selected) < min.rows)) {
      selected <- c(selected, idx)
      selected.mat <- trial.mat
      current.rank <- trial.rank
      covered <- union(covered, row.objects[[idx]]$j)
    }
    if (length(selected) >= target.rows) break
  }
  if (length(selected) < min.rows) {
    selected <- unique(c(selected, head(order.idx, min.rows)))
  }
  if (length(selected) < 1L) selected <- order.idx[seq_len(min(1L, length(order.idx)))]
  for (idx in seq_along(row.objects)) {
    row.objects[[idx]]$polynomial.residual <- residuals[idx]
  }
  selected.objects <- row.objects[selected]
  for (idx in seq_along(selected.objects)) {
    selected.objects[[idx]]$selection.rank <- idx
  }
  list(
    rows = selected.objects,
    diagnostics = data.frame(
      row.selection = "polynomial.residual.greedy",
      candidate.rows = n.candidates,
      selected.rows = length(selected.objects),
      thinning.local.dim = local.dim,
      thinning.polynomial.dim = p.dim,
      thinning.rank = if (nrow(selected.mat)) qr(selected.mat, LAPACK = TRUE)$rank else 0L,
      thinning.coverage = length(unique(unlist(lapply(selected.objects, `[[`, "j"), use.names = FALSE))),
      residual.min = if (length(finite.res)) min(finite.res) else NA_real_,
      residual.median = if (length(finite.res)) stats::median(finite.res) else NA_real_,
      residual.max = if (length(finite.res)) max(finite.res) else NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

select.center.row.objects <- function(row.objects,
                                      row.selection = c("none", "polynomial.residual.greedy"),
                                      disk.vertices,
                                      D.sub,
                                      degree,
                                      dims = 1:5,
                                      max.rows = NULL,
                                      min.rows = 1L) {
  row.selection <- match.arg(row.selection)
  if (row.selection == "none") {
    for (idx in seq_along(row.objects)) {
      row.objects[[idx]]$polynomial.residual <- NA_real_
      row.objects[[idx]]$selection.rank <- idx
    }
    return(list(
      rows = row.objects,
      diagnostics = data.frame(
        row.selection = "none",
        candidate.rows = length(row.objects),
        selected.rows = length(row.objects),
        thinning.local.dim = NA_integer_,
        thinning.polynomial.dim = NA_integer_,
        thinning.rank = NA_integer_,
        thinning.coverage = length(unique(unlist(lapply(row.objects, `[[`, "j"), use.names = FALSE))),
        residual.min = NA_real_,
        residual.median = NA_real_,
        residual.max = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }
  thin.center.rows.polynomial.residual(
    row.objects = row.objects,
    disk.vertices = disk.vertices,
    D.sub = D.sub,
    degree = degree,
    dims = dims,
    max.rows = max.rows,
    min.rows = min.rows
  )
}

local.dimension.diagnostics <- function(adj.list,
                                        weight.list,
                                        D,
                                        center,
                                        disk.vertices,
                                        method = c("cmdscale",
                                                   "weighted.grip.edge.kk",
                                                   "metric.mds.edge.kk",
                                                   "hybrid.edge.kk"),
                                        dims = 1:5,
                                        seed = 1L) {
  method <- match.arg(method)
  disk.vertices <- sort(unique(as.integer(disk.vertices)))
  edge.data <- induced.edge.data(adj.list, weight.list, disk.vertices)
  if (length(disk.vertices) < 2L || nrow(edge.data$edges) < 1L) {
    out <- data.frame(
      center = center,
      local.embedding.method = method,
      embedding.dim = dims,
      selected.local.dim = NA_integer_,
      best.local.dim = NA_integer_,
      dimension.cap.hit = NA,
      dimension.selection.rule = "insufficient.local.graph",
      n.local.vertices = length(disk.vertices),
      n.local.edges = nrow(edge.data$edges),
      edge.rrmse = NA_real_,
      edge.ratio.q95 = NA_real_,
      edge.ratio.gt.1p10 = NA_real_,
      spread.trace = NA_real_,
      status = "insufficient.local.graph",
      message = "",
      stringsAsFactors = FALSE
    )
    attr(out, "vote.table") <- data.frame(
      center = center,
      local.embedding.method = method,
      dimension.method = c("edge.stress.elbow", "mds.spectrum.elbow", "graph.laplacian.weyl"),
      estimated.dimension = NA_integer_,
      confidence = NA_real_,
      diagnostic.value = NA_real_,
      rule = "insufficient.local.graph",
      status = "insufficient.local.graph",
      message = "",
      stringsAsFactors = FALSE
    )
    return(out)
  }
  D.sub <- D[disk.vertices, disk.vertices, drop = FALSE]
  records <- lapply(seq_along(dims), function(idx) {
    dim <- as.integer(dims[idx])
    if (method == "weighted.grip.edge.kk" && !dim %in% c(2L, 3L)) {
      return(data.frame(
        center = center,
        local.embedding.method = method,
        embedding.dim = dim,
        n.local.vertices = length(disk.vertices),
        n.local.edges = nrow(edge.data$edges),
        status = "unsupported.dimension",
        message = "weighted GRIP local diagnostics currently support dim = 2 or 3",
        edge.rrmse = NA_real_,
        edge.ratio.q95 = NA_real_,
        edge.ratio.gt.1p10 = NA_real_,
        spread.trace = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    fit <- tryCatch({
      coords <- if (method == "cmdscale") {
        local.cmdscale.coords(D.sub, dim)
      } else if (method == "metric.mds.edge.kk") {
        local.metric.mds.edge.kk.coords(
          edges = edge.data$edges,
          edge.weights = edge.data$edge.weights,
          n = length(disk.vertices),
          dim = dim,
          D.sub = D.sub
        )
      } else if (method == "hybrid.edge.kk") {
        local.hybrid.edge.kk.coords(
          edges = edge.data$edges,
          edge.weights = edge.data$edge.weights,
          n = length(disk.vertices),
          dim = dim,
          seed = seed + 1000L * center + dim,
          D.sub = D.sub
        )
      } else {
        local.weighted.grip.edge.kk.coords(
          edges = edge.data$edges,
          edge.weights = edge.data$edge.weights,
          n = length(disk.vertices),
          dim = dim,
          seed = seed + 1000L * center + dim
        )
      }
      note <- attr(coords, "embedding.note") %||% ""
      coords <- scale.coords.to.edges(coords, edge.data$edges, edge.data$edge.weights)
      diag <- edge.stress.diagnostics(coords, edge.data$edges, edge.data$edge.weights)
      cbind(data.frame(
        center = center,
        local.embedding.method = method,
        embedding.dim = dim,
        n.local.vertices = length(disk.vertices),
        n.local.edges = nrow(edge.data$edges),
        status = if (is.null(coords)) "embedding.failed" else "ok",
        message = note,
        stringsAsFactors = FALSE
      ), diag)
    }, error = function(e) {
      data.frame(
        center = center,
        local.embedding.method = method,
        embedding.dim = dim,
        n.local.vertices = length(disk.vertices),
        n.local.edges = nrow(edge.data$edges),
        status = "error",
        message = conditionMessage(e),
        edge.rrmse = NA_real_,
        edge.ratio.q95 = NA_real_,
        edge.ratio.gt.1p10 = NA_real_,
        spread.trace = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    fit
  })
  out <- do.call(rbind, records)
  selected <- select.local.dimension(out)
  out$selected.local.dim <- selected$selected.dim
  out$best.local.dim <- selected$best.dim
  out$dimension.cap.hit <- selected$dimension.cap.hit
  out$dimension.selection.rule <- selected$selection.rule
  vote.table <- local.dimension.vote.table(
    dim.table = out,
    D.sub = D.sub,
    edge.data = edge.data,
    dims = dims
  )
  vote.table$center <- center
  vote.table$local.embedding.method <- method
  consensus <- consensus.dimension.from.votes(vote.table)
  out$mds.spectrum.local.dim <- vote.table$estimated.dimension[
    vote.table$dimension.method == "mds.spectrum.elbow"
  ][1L]
  out$graph.spectrum.local.dim <- vote.table$estimated.dimension[
    vote.table$dimension.method == "graph.laplacian.weyl"
  ][1L]
  out$consensus.local.dim <- consensus
  attr(out, "vote.table") <- vote.table
  out
}

candidate.geodesics.for.center <- function(adj.list,
                                           weight.list,
                                           D,
                                           center,
                                           disk.k,
                                           geodesic.tol = 1e-8) {
  finite.dist <- D[center, is.finite(D[center, ]) & D[center, ] > 0]
  if (length(finite.dist) < 1L) {
    return(list(
      disk.radius = NA_real_,
      disk.k.actual = NA_integer_,
      ray.paths = list(),
      disk.vertices = center,
      candidates = list()
    ))
  }
  disk.k.actual <- min(as.integer(disk.k), length(finite.dist))
  disk.radius <- sort(finite.dist)[disk.k.actual]
  sp <- find.shortest.paths.within.radius(adj.list, weight.list, center, disk.radius)
  ray.paths <- sp$paths
  ray.paths <- ray.paths[vapply(ray.paths, length, integer(1)) >= 2L]
  endpoints <- vapply(ray.paths, function(p) tail(p, 1L), integer(1))
  keep <- endpoints != center
  ray.paths <- ray.paths[keep]
  endpoints <- endpoints[keep]
  disk.vertices <- unique(unlist(ray.paths, use.names = FALSE))
  if (!center %in% disk.vertices) disk.vertices <- c(center, disk.vertices)

  candidates <- list()
  candidate.keys <- character()
  if (length(ray.paths) >= 2L) {
    for (a in seq_len(length(ray.paths) - 1L)) {
      for (b in seq.int(a + 1L, length(ray.paths))) {
        shared <- intersect(ray.paths[[a]][-1L], ray.paths[[b]][-1L])
        if (length(shared) > 0L) next
        endpoint.a <- endpoints[a]
        endpoint.b <- endpoints[b]
        composite.length <- D[center, endpoint.a] + D[center, endpoint.b]
        if (!is.finite(D[endpoint.a, endpoint.b])) next
        if (abs(D[endpoint.a, endpoint.b] - composite.length) >
            geodesic.tol * max(1, composite.length)) next
        cp <- make.composite.path(ray.paths[[a]], ray.paths[[b]], adj.list, weight.list)
        key <- path.key(cp$vertices)
        if (key %in% candidate.keys) next
        candidate.keys <- c(candidate.keys, key)
        candidates[[length(candidates) + 1L]] <- list(
          vertices = cp$vertices,
          x = cp$x,
          source = "composite.geodesic"
        )
      }
    }
  }
  list(
    disk.radius = disk.radius,
    disk.k.actual = disk.k.actual,
    ray.paths = ray.paths,
    disk.vertices = disk.vertices,
    candidates = candidates
  )
}

build.geodesic.annihilator.operator <- function(adj.list,
                                                weight.list,
                                                degree = 2L,
                                                path.family = c("local.composite.geodesics",
                                                                "geodesic.skeleton",
                                                                "geodesic.skeleton.cross.enriched"),
                                                disk.k = 25L,
                                                disk.radius.rule = c("fixed.k",
                                                                     "min.composite.geodesics"),
                                                disk.k.start = disk.k,
                                                disk.k.step = 4L,
                                                disk.k.max = max(disk.k, 60L),
                                                max.composite.paths.per.vertex = 50L,
                                                min.composite.paths.per.vertex = 2L,
                                                mirror = c("auto", "never", "always"),
                                                row.normalization = c("unit.l2", "none"),
                                                vertex.normalization = c("mean", "none"),
                                                row.selection = c("none",
                                                                  "polynomial.residual.greedy"),
                                                row.selection.max.rows.per.vertex = NULL,
                                                row.selection.min.rows.per.vertex = 1L,
                                                center.order = c("increasing.degree",
                                                                 "natural",
                                                                 "decreasing.degree"),
                                                candidate.cap.mode = c("deterministic",
                                                                       "random"),
                                                max.skeleton.paths = NULL,
                                                skeleton.min.coverage = 2L,
                                                skeleton.max.passes = 2L,
                                                skeleton.max.cross.paths = 30L,
                                                skeleton.cross.rule = c("coverage",
                                                                        "local.dimension.rank"),
                                                skeleton.cross.embedding.method = c("cmdscale",
                                                                                    "hybrid.edge.kk"),
                                                skeleton.cross.candidate.pool = c("skeleton.membership",
                                                                                  "rich.local"),
                                                skeleton.cross.local.dims = 1:5,
                                                skeleton.cross.rank.tol = 1e-6,
                                                skeleton.cross.max.candidates.per.center = 50L,
                                                skeleton.cross.length.penalty = 0.05,
                                                skeleton.cross.coverage.reward = 0.10,
                                                skeleton.cross.iterations = 1L,
                                                local.dimension.method = c("none",
                                                                           "cmdscale",
                                                                           "weighted.grip.edge.kk",
                                                                           "metric.mds.edge.kk",
                                                                           "hybrid.edge.kk"),
                                                local.dimension.max.centers = 0L,
                                                local.dimension.dims = 1:5,
                                                seed = 1L,
                                                geodesic.tol = 1e-8) {
  path.family <- match.arg(path.family)
  disk.radius.rule <- match.arg(disk.radius.rule)
  mirror <- match.arg(mirror)
  row.normalization <- match.arg(row.normalization)
  vertex.normalization <- match.arg(vertex.normalization)
  row.selection <- match.arg(row.selection)
  center.order <- match.arg(center.order)
  candidate.cap.mode <- match.arg(candidate.cap.mode)
  skeleton.cross.rule <- match.arg(skeleton.cross.rule)
  skeleton.cross.embedding.method <- match.arg(skeleton.cross.embedding.method)
  skeleton.cross.candidate.pool <- match.arg(skeleton.cross.candidate.pool)
  local.dimension.method <- match.arg(local.dimension.method)
  if (path.family != "local.composite.geodesics") {
    return(build.geodesic.skeleton.annihilator.operator(
      adj.list = adj.list,
      weight.list = weight.list,
      degree = degree,
      path.family = path.family,
      disk.k = disk.k,
      max.skeleton.paths = max.skeleton.paths,
      skeleton.min.coverage = skeleton.min.coverage,
      skeleton.max.passes = skeleton.max.passes,
      skeleton.max.cross.paths = skeleton.max.cross.paths,
      skeleton.cross.rule = skeleton.cross.rule,
      skeleton.cross.embedding.method = skeleton.cross.embedding.method,
      skeleton.cross.candidate.pool = skeleton.cross.candidate.pool,
      skeleton.cross.local.dims = skeleton.cross.local.dims,
      skeleton.cross.rank.tol = skeleton.cross.rank.tol,
      skeleton.cross.max.candidates.per.center = skeleton.cross.max.candidates.per.center,
      skeleton.cross.length.penalty = skeleton.cross.length.penalty,
      skeleton.cross.coverage.reward = skeleton.cross.coverage.reward,
      skeleton.cross.iterations = skeleton.cross.iterations,
      row.normalization = row.normalization,
      vertex.normalization = vertex.normalization,
      seed = seed,
      geodesic.tol = geodesic.tol
    ))
  }
  n <- length(adj.list)
  set.seed(seed)
  local.dimension.max.centers <- as.integer(local.dimension.max.centers)
  ordered.centers <- order.centers.by.rule(adj.list, center.order)
  diagnostic.centers <- integer()
  if (local.dimension.method != "none" && local.dimension.max.centers != 0L) {
    if (is.infinite(local.dimension.max.centers) || local.dimension.max.centers >= n) {
      diagnostic.centers <- ordered.centers
    } else if (local.dimension.max.centers > 0L) {
      diagnostic.centers <- head(ordered.centers, local.dimension.max.centers)
    }
  }
  D <- shortest.path(adj.list, weight.list, seq_len(n))
  i.idx <- integer()
  j.idx <- integer()
  x.val <- numeric()
  row <- 0L
  path.records <- list()
  vertex.records <- vector("list", n)
  local.dimension.records <- list()
  local.dimension.vote.records <- list()

  for (center.position in seq_along(ordered.centers)) {
    center <- ordered.centers[center.position]
    if (disk.radius.rule == "fixed.k") {
      selected <- candidate.geodesics.for.center(
        adj.list, weight.list, D, center, disk.k, geodesic.tol
      )
      attempted.disk.k <- selected$disk.k.actual
      radius.status <- "fixed.k"
    } else {
      attempted <- seq.int(as.integer(disk.k.start), as.integer(disk.k.max),
                           by = as.integer(disk.k.step))
      if (!as.integer(disk.k.max) %in% attempted) {
        attempted <- c(attempted, as.integer(disk.k.max))
      }
      selected <- NULL
      for (candidate.disk.k in attempted) {
        selected <- candidate.geodesics.for.center(
          adj.list, weight.list, D, center, candidate.disk.k, geodesic.tol
        )
        if (length(selected$candidates) >= min.composite.paths.per.vertex) break
      }
      attempted.disk.k <- paste(attempted, collapse = ",")
      radius.status <- if (length(selected$candidates) >= min.composite.paths.per.vertex) {
        "met.min.composite.geodesics"
      } else {
        "hit.disk.k.max"
      }
    }

    disk.radius <- selected$disk.radius
    selected.disk.k <- selected$disk.k.actual
    ray.paths <- selected$ray.paths
    disk.vertices <- selected$disk.vertices
    candidates <- selected$candidates
    n.composite.available <- length(candidates)
    if (center %in% diagnostic.centers) {
      ldt <- local.dimension.diagnostics(
        adj.list = adj.list,
        weight.list = weight.list,
        D = D,
        center = center,
        disk.vertices = disk.vertices,
        method = local.dimension.method,
        dims = local.dimension.dims,
        seed = seed
      )
      local.dimension.records[[length(local.dimension.records) + 1L]] <- ldt
      vote.table <- attr(ldt, "vote.table")
      if (!is.null(vote.table) && nrow(vote.table) > 0L) {
        local.dimension.vote.records[[length(local.dimension.vote.records) + 1L]] <- vote.table
      }
    }
    if (length(candidates) > max.composite.paths.per.vertex) {
      if (candidate.cap.mode == "deterministic") {
        candidates <- head(sort.candidates.deterministic(candidates),
                           max.composite.paths.per.vertex)
      } else {
        candidates <- candidates[sample(seq_along(candidates), max.composite.paths.per.vertex)]
      }
    } else {
      candidates <- sort.candidates.deterministic(candidates)
    }

    use.one.sided <- mirror == "always" ||
      (mirror == "auto" && length(candidates) < min.composite.paths.per.vertex)
    if (use.one.sided && length(ray.paths) > 0L) {
      one.sided <- lapply(ray.paths, function(p) {
        list(vertices = as.integer(p),
             x = path.cumulative.distances(p, adj.list, weight.list),
             source = "one.sided.mirror.proxy")
      })
      enough <- vapply(one.sided, function(z) length(z$vertices) >= degree + 2L, logical(1))
      one.sided <- one.sided[enough]
      slots <- max.composite.paths.per.vertex - length(candidates)
      if (length(one.sided) > slots) {
        if (candidate.cap.mode == "deterministic") {
          one.sided <- head(sort.candidates.deterministic(one.sided), slots)
        } else {
          one.sided <- one.sided[sample(seq_along(one.sided), slots)]
        }
      } else {
        one.sided <- sort.candidates.deterministic(one.sided)
      }
      candidates <- c(candidates, one.sided)
    }

    center.row.objects <- list()
    for (candidate.idx in seq_along(candidates)) {
      rows <- annihilator.rows.for.path(
        candidates[[candidate.idx]]$vertices,
        candidates[[candidate.idx]]$x,
        n.vertices = n,
        degree = degree,
        row.normalization = row.normalization
      )
      if (length(rows) < 1L) next
      for (basis.idx in seq_along(rows)) {
        center.row.objects[[length(center.row.objects) + 1L]] <- list(
          j = rows[[basis.idx]]$j,
          x = rows[[basis.idx]]$x,
          candidate.idx = candidate.idx,
          basis.idx = basis.idx,
          source = candidates[[candidate.idx]]$source,
          n.path.vertices = length(candidates[[candidate.idx]]$vertices),
          path.length = max(candidates[[candidate.idx]]$x) -
            min(candidates[[candidate.idx]]$x),
          path.vertices = paste(candidates[[candidate.idx]]$vertices, collapse = "-")
        )
      }
    }
    D.sub <- D[disk.vertices, disk.vertices, drop = FALSE]
    selected.center <- select.center.row.objects(
      center.row.objects,
      row.selection = row.selection,
      disk.vertices = disk.vertices,
      D.sub = D.sub,
      degree = degree,
      dims = local.dimension.dims,
      max.rows = row.selection.max.rows.per.vertex,
      min.rows = row.selection.min.rows.per.vertex
    )
    selected.row.objects <- selected.center$rows
    center.start.row <- row + 1L
    for (row.obj.idx in seq_along(selected.row.objects)) {
      row.obj <- selected.row.objects[[row.obj.idx]]
      row <- row + 1L
      row.x <- row.obj$x
      if (vertex.normalization == "mean" && length(selected.row.objects) > 0L) {
        row.x <- row.x / length(selected.row.objects)
      }
      i.idx <- c(i.idx, rep(row, length(row.obj$j)))
      j.idx <- c(j.idx, row.obj$j)
      x.val <- c(x.val, row.x)
      path.records[[row]] <- data.frame(
          row = row,
          center = center,
          candidate = row.obj$candidate.idx,
          basis.row = row.obj$basis.idx,
          degree = degree,
          source = row.obj$source,
          n.path.vertices = row.obj$n.path.vertices,
          path.length = row.obj$path.length,
          path.vertices = row.obj$path.vertices,
          row.selection = row.selection,
          center.order = center.order,
          center.processing.rank = center.position,
          candidate.cap.mode = candidate.cap.mode,
          path.family = path.family,
          skeleton.path.id = NA_integer_,
          skeleton.pass = NA_integer_,
          skeleton.role = NA_character_,
          window.start = NA_integer_,
          polynomial.residual = row.obj$polynomial.residual %||% NA_real_,
          selection.rank = row.obj$selection.rank %||% NA_integer_,
          stringsAsFactors = FALSE
        )
    }
    center.end.row <- row
    thinning.diag <- selected.center$diagnostics
    vertex.records[[center]] <- data.frame(
      center = center,
      center.order = center.order,
      center.processing.rank = center.position,
      degree = degree,
      disk.radius.rule = disk.radius.rule,
      disk.k = disk.k,
      disk.k.start = disk.k.start,
      disk.k.step = disk.k.step,
      disk.k.max = disk.k.max,
      selected.disk.k = selected.disk.k,
      attempted.disk.k = attempted.disk.k,
      radius.status = radius.status,
      disk.radius = disk.radius,
      n.disk.vertices = length(disk.vertices),
      n.rays = length(ray.paths),
      n.composite.available = n.composite.available,
      n.paths.used = length(candidates),
      n.operator.rows = max(0L, center.end.row - center.start.row + 1L),
      row.selection = row.selection,
      candidate.cap.mode = candidate.cap.mode,
      path.family = path.family,
      skeleton.coverage = NA_integer_,
      skeleton.distance = NA_real_,
      candidate.operator.rows = thinning.diag$candidate.rows,
      selected.operator.rows = thinning.diag$selected.rows,
      thinning.local.dim = thinning.diag$thinning.local.dim,
      thinning.polynomial.dim = thinning.diag$thinning.polynomial.dim,
      thinning.rank = thinning.diag$thinning.rank,
      thinning.coverage = thinning.diag$thinning.coverage,
      thinning.residual.min = thinning.diag$residual.min,
      thinning.residual.median = thinning.diag$residual.median,
      thinning.residual.max = thinning.diag$residual.max,
      used.mirror.proxy = any(vapply(candidates, function(z) {
        identical(z$source, "one.sided.mirror.proxy")
      }, logical(1))),
      stringsAsFactors = FALSE
    )
  }

  if (row == 0L) {
    A <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                              dims = c(0L, n))
  } else {
    A <- Matrix::sparseMatrix(i = i.idx, j = j.idx, x = x.val, dims = c(row, n))
  }
  list(
    A = A,
    path.table = if (length(path.records)) do.call(rbind, path.records) else data.frame(),
    vertex.table = do.call(rbind, vertex.records),
    control = list(
      degree = degree,
      path.family = path.family,
      disk.k = disk.k,
      disk.radius.rule = disk.radius.rule,
      disk.k.start = disk.k.start,
      disk.k.step = disk.k.step,
      disk.k.max = disk.k.max,
      max.composite.paths.per.vertex = max.composite.paths.per.vertex,
      min.composite.paths.per.vertex = min.composite.paths.per.vertex,
      mirror = mirror,
      row.normalization = row.normalization,
      vertex.normalization = vertex.normalization,
      row.selection = row.selection,
      row.selection.max.rows.per.vertex = row.selection.max.rows.per.vertex,
      row.selection.min.rows.per.vertex = row.selection.min.rows.per.vertex,
      center.order = center.order,
      candidate.cap.mode = candidate.cap.mode,
      max.skeleton.paths = max.skeleton.paths,
      skeleton.min.coverage = skeleton.min.coverage,
      skeleton.max.passes = skeleton.max.passes,
      skeleton.max.cross.paths = skeleton.max.cross.paths,
      local.dimension.method = local.dimension.method,
      local.dimension.max.centers = local.dimension.max.centers,
      local.dimension.dims = paste(local.dimension.dims, collapse = ","),
      seed = seed
    ),
    local.dimension.table = if (length(local.dimension.records)) {
      do.call(rbind, local.dimension.records)
    } else {
      data.frame()
    },
    local.dimension.vote.table = if (length(local.dimension.vote.records)) {
      do.call(rbind, local.dimension.vote.records)
    } else {
      data.frame()
    }
  )
}

sample.square.points <- function(n, seed) {
  set.seed(seed)
  cbind(stats::runif(n), stats::runif(n))
}

embed.quadric <- function(U, surface = c("paraboloid", "saddle", "elliptic")) {
  surface <- match.arg(surface)
  u <- U[, 1L]
  v <- U[, 2L]
  z <- switch(surface,
              paraboloid = 0.7 * ((u - 0.5)^2 + (v - 0.5)^2),
              saddle = 0.7 * ((u - 0.5)^2 - (v - 0.5)^2),
              elliptic = 0.5 * (2 * (u - 0.5)^2 + 0.5 * (v - 0.5)^2))
  cbind(u, v, z)
}

score.operator.polynomials <- function(A, U, X3, degree) {
  intrinsic <- monomial.design.2d(U, degree)
  intrinsic.control <- monomial.design.2d(U, degree + 1L, exact.degree = TRUE)
  ambient.linear <- cbind(1, X3)
  colnames(ambient.linear) <- c("1", "X", "Y", "Z")
  ambient.quadratic <- cbind(
    ambient.linear,
    X3[, 1L]^2,
    X3[, 1L] * X3[, 2L],
    X3[, 2L]^2,
    X3[, 1L] * X3[, 3L],
    X3[, 2L] * X3[, 3L],
    X3[, 3L]^2
  )
  colnames(ambient.quadratic) <- c("1", "X", "Y", "Z", "X2", "XY", "Y2",
                                   "XZ", "YZ", "Z2")
  data.frame(
    target = c("intrinsic.degree.leq.d",
               "intrinsic.exact.degree.d.plus.1",
               "ambient.linear",
               "ambient.quadratic"),
    residual.ratio = c(
      frobenius.residual.ratio(A, intrinsic),
      frobenius.residual.ratio(A, intrinsic.control),
      frobenius.residual.ratio(A, ambient.linear),
      frobenius.residual.ratio(A, ambient.quadratic)
    ),
    stringsAsFactors = FALSE
  )
}

run.geodesic.annihilator.experiment <- function(mode = c("smoke", "full"),
                                                output.dir = file.path(repo.dir,
                                                  "dev/graph-trend-filtering/reports"),
                                                seed.base = 20260517L) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  if (mode == "smoke") {
    config <- expand.grid(
      dataset = c("square", "paraboloid"),
      n = 60L,
      seed = 1L,
      graph.type = c("delaunay.1skeleton", "symmetric.knn"),
      degree = 1:2,
      path.family = "local.composite.geodesics",
      disk.radius.rule = c("fixed.k", "min.composite.geodesics"),
      disk.k = 15L,
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 48L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = c("none", "polynomial.residual.greedy"),
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "cmdscale",
      local.dimension.max.centers = 8L,
      stringsAsFactors = FALSE
    )
    config$skeleton.cross.rule <- "coverage"
    config$skeleton.min.coverage <- 2L
    config$skeleton.max.paths <- NA_integer_
    skeleton.config <- expand.grid(
      dataset = c("square", "paraboloid"),
      n = 60L,
      seed = 1L,
      graph.type = c("delaunay.1skeleton", "symmetric.knn"),
      degree = 1:2,
      path.family = c("geodesic.skeleton", "geodesic.skeleton.cross.enriched"),
      disk.radius.rule = "fixed.k",
      disk.k = 15L,
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 48L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = "none",
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "none",
      local.dimension.max.centers = 0L,
      stringsAsFactors = FALSE
    )
    skeleton.config$skeleton.cross.rule <- "coverage"
    skeleton.config$skeleton.min.coverage <- 2L
    skeleton.config$skeleton.max.paths <- 60L
    skeleton.rank.config <- skeleton.config[
      skeleton.config$path.family == "geodesic.skeleton.cross.enriched",
      ,
      drop = FALSE
    ]
    skeleton.rank.config$skeleton.cross.rule <- "local.dimension.rank"
    skeleton.rank.config$skeleton.min.coverage <- 1L
    skeleton.rank.config$skeleton.max.paths <- 8L
    skeleton.config <- rbind(skeleton.config, skeleton.rank.config)
    weighted.config <- expand.grid(
      dataset = "square",
      n = 60L,
      seed = 1L,
      graph.type = c("delaunay.1skeleton", "symmetric.knn"),
      degree = 2L,
      path.family = "local.composite.geodesics",
      disk.radius.rule = "min.composite.geodesics",
      disk.k = 15L,
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 40L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = "none",
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "hybrid.edge.kk",
      local.dimension.max.centers = 3L,
      stringsAsFactors = FALSE
    )
    weighted.config$skeleton.cross.rule <- "coverage"
    weighted.config$skeleton.min.coverage <- 2L
    weighted.config$skeleton.max.paths <- NA_integer_
    config <- unique(rbind(config, skeleton.config, weighted.config))
  } else {
    fixed.config <- expand.grid(
      dataset = c("square", "paraboloid", "saddle"),
      n = 80L,
      seed = 1:2,
      graph.type = c("delaunay.1skeleton", "symmetric.knn", "radius", "cknn"),
      degree = 1:3,
      path.family = "local.composite.geodesics",
      disk.k = c(18L, 30L),
      disk.radius.rule = "fixed.k",
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 60L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = "none",
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "cmdscale",
      local.dimension.max.centers = 8L,
      stringsAsFactors = FALSE
    )
    fixed.config$skeleton.cross.rule <- "coverage"
    fixed.config$skeleton.min.coverage <- 2L
    fixed.config$skeleton.max.paths <- NA_integer_
    adaptive.config <- expand.grid(
      dataset = c("square", "paraboloid", "saddle"),
      n = 80L,
      seed = 1:2,
      graph.type = c("delaunay.1skeleton", "symmetric.knn", "radius", "cknn"),
      degree = 1:3,
      path.family = "local.composite.geodesics",
      disk.k = 18L,
      disk.radius.rule = "min.composite.geodesics",
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 48L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = c("none", "polynomial.residual.greedy"),
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "cmdscale",
      local.dimension.max.centers = 8L,
      stringsAsFactors = FALSE
    )
    adaptive.config$skeleton.cross.rule <- "coverage"
    adaptive.config$skeleton.min.coverage <- 2L
    adaptive.config$skeleton.max.paths <- NA_integer_
    base.config <- rbind(fixed.config, adaptive.config)
    base.config$analysis.group <- "primary"
    skeleton.config <- expand.grid(
      dataset = c("square", "paraboloid", "saddle"),
      n = 80L,
      seed = 1:2,
      graph.type = c("delaunay.1skeleton", "symmetric.knn", "cknn"),
      degree = 1:3,
      path.family = c("geodesic.skeleton", "geodesic.skeleton.cross.enriched"),
      disk.radius.rule = "fixed.k",
      disk.k = 18L,
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 48L,
      min.composite.geodesics = 10L,
      row.normalization = "unit.l2",
      vertex.normalization = "mean",
      row.selection = "none",
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "none",
      local.dimension.max.centers = 0L,
      stringsAsFactors = FALSE
    )
    skeleton.config$skeleton.cross.rule <- "coverage"
    skeleton.config$skeleton.min.coverage <- 2L
    skeleton.config$skeleton.max.paths <- 100L
    skeleton.rank.config <- skeleton.config[
      skeleton.config$path.family == "geodesic.skeleton.cross.enriched",
      ,
      drop = FALSE
    ]
    skeleton.rank.config$skeleton.cross.rule <- "local.dimension.rank"
    skeleton.rank.config$skeleton.min.coverage <- 1L
    skeleton.rank.config$skeleton.max.paths <- 16L
    skeleton.config <- rbind(skeleton.config, skeleton.rank.config)
    skeleton.config$analysis.group <- "geodesic.skeleton"
    norm.config <- expand.grid(
      dataset = c("square", "paraboloid"),
      n = 80L,
      seed = 1L,
      graph.type = c("delaunay.1skeleton", "symmetric.knn"),
      degree = 1:3,
      path.family = "local.composite.geodesics",
      disk.radius.rule = c("fixed.k", "min.composite.geodesics"),
      disk.k = 18L,
      disk.k.start = 12L,
      disk.k.step = 4L,
      disk.k.max = 48L,
      min.composite.geodesics = 10L,
      row.normalization = c("unit.l2", "none"),
      vertex.normalization = c("mean", "none"),
      row.selection = "none",
      center.order = "increasing.degree",
      candidate.cap.mode = "deterministic",
      local.dimension.method = "cmdscale",
      local.dimension.max.centers = 8L,
      stringsAsFactors = FALSE
    )
    norm.config$skeleton.cross.rule <- "coverage"
    norm.config$skeleton.min.coverage <- 2L
    norm.config$skeleton.max.paths <- NA_integer_
    norm.config$analysis.group <- "normalization.sensitivity"
    config <- unique(rbind(base.config, skeleton.config, norm.config))
  }
  config$run.id <- seq_len(nrow(config))

  result.records <- list()
  nullity.records <- list()
  singular.records <- list()
  vertex.records <- list()
  path.records <- list()
  skeleton.records <- list()
  skeleton.coverage.records <- list()
  skeleton.direction.records <- list()
  skeleton.connector.candidate.records <- list()
  local.dimension.records <- list()
  local.dimension.vote.records <- list()
  for (idx in seq_len(nrow(config))) {
    if (idx == 1L || idx %% 10L == 0L || idx == nrow(config)) {
      message(sprintf("geodesic annihilator %s run %d/%d", mode, idx, nrow(config)))
    }
    cfg <- config[idx, ]
    row.names(cfg) <- NULL
    U <- sample.square.points(cfg$n, seed.base + cfg$seed)
    X3 <- if (cfg$dataset == "square") cbind(U, 0) else embed.quadric(U, cfg$dataset)
    metric.coords <- if (cfg$dataset == "square") U else X3
    graph <- tryCatch(
      make.experimental.graph(
        U = U,
        graph.type = cfg$graph.type,
        metric.coords = metric.coords,
        k = 8L
      ),
      error = function(e) e
    )
    if (inherits(graph, "error")) {
      result.records[[idx]] <- cbind(cfg, data.frame(
        status = "graph.error",
        message = conditionMessage(graph),
        n.edges = NA_real_,
        n.operator.rows = NA_integer_,
        nullity.1e.6 = NA_integer_,
        expected.nullity = choose(cfg$degree + 2L, 2L),
        intrinsic.residual = NA_real_,
        control.residual = NA_real_,
        stringsAsFactors = FALSE
      ))
      next
    }
    op <- tryCatch(
      build.geodesic.annihilator.operator(
        graph$adj.list,
        graph$weight.list,
        degree = cfg$degree,
        path.family = cfg$path.family %||% "local.composite.geodesics",
        disk.k = cfg$disk.k,
        disk.radius.rule = cfg$disk.radius.rule,
        disk.k.start = cfg$disk.k.start,
        disk.k.step = cfg$disk.k.step,
        disk.k.max = cfg$disk.k.max,
        max.composite.paths.per.vertex = if (mode == "smoke") 30L else 50L,
        min.composite.paths.per.vertex = cfg$min.composite.geodesics,
        mirror = "auto",
        row.normalization = cfg$row.normalization,
        vertex.normalization = cfg$vertex.normalization,
        row.selection = cfg$row.selection,
        row.selection.max.rows.per.vertex = NULL,
        row.selection.min.rows.per.vertex = 1L,
        center.order = cfg$center.order %||% "increasing.degree",
        candidate.cap.mode = cfg$candidate.cap.mode %||% "deterministic",
        max.skeleton.paths = if (!is.null(cfg$skeleton.max.paths) &&
            !is.na(cfg$skeleton.max.paths)) {
          as.integer(cfg$skeleton.max.paths)
        } else if (mode == "smoke") {
          60L
        } else {
          100L
        },
        skeleton.min.coverage = cfg$skeleton.min.coverage %||% 2L,
        skeleton.max.passes = 2L,
        skeleton.max.cross.paths = if (mode == "smoke") 12L else 30L,
        skeleton.cross.rule = cfg$skeleton.cross.rule %||% "coverage",
        skeleton.cross.embedding.method = "cmdscale",
        skeleton.cross.local.dims = 1:5,
        skeleton.cross.rank.tol = 1e-6,
        skeleton.cross.max.candidates.per.center = if (mode == "smoke") 30L else 50L,
        skeleton.cross.length.penalty = 0.05,
        skeleton.cross.coverage.reward = 0.10,
        local.dimension.method = cfg$local.dimension.method,
        local.dimension.max.centers = cfg$local.dimension.max.centers,
        local.dimension.dims = cfg$local.dimension.dims %||% 1:5,
        seed = seed.base + idx
      ),
      error = function(e) e
    )
    if (inherits(op, "error")) {
      result.records[[idx]] <- cbind(cfg, data.frame(
        status = "operator.error",
        message = conditionMessage(op),
        n.edges = length(graph$edge.weights),
        n.operator.rows = NA_integer_,
        nullity.1e.6 = NA_integer_,
        expected.nullity = choose(cfg$degree + 2L, 2L),
        intrinsic.residual = NA_real_,
        control.residual = NA_real_,
        stringsAsFactors = FALSE
      ))
      next
    }
    nullity <- estimate.nullity(op$A)
    poly.scores <- score.operator.polynomials(op$A, U, X3, cfg$degree)
    nullity.records[[idx]] <- cbind(cfg, nullity$table)
    singular.records[[idx]] <- if (length(nullity$singular.values) > 0L) {
      data.frame(
        run.id = cfg$run.id,
        singular.index = seq_along(nullity$singular.values),
        singular.value = nullity$singular.values
      )
    } else {
      data.frame(
        run.id = integer(),
        singular.index = integer(),
        singular.value = numeric()
      )
    }
    vt <- op$vertex.table
    vt$run.id <- cfg$run.id
    vt$dataset <- cfg$dataset
    vt$n <- cfg$n
    vt$seed <- cfg$seed
      vt$graph.type <- cfg$graph.type
      vt$path.family <- cfg$path.family %||% "local.composite.geodesics"
      vt$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      vt$disk.radius.rule <- cfg$disk.radius.rule
      vt$row.normalization <- cfg$row.normalization
      vt$vertex.normalization <- cfg$vertex.normalization
      vt$row.selection <- cfg$row.selection
      vt$center.order <- cfg$center.order %||% "increasing.degree"
      vt$candidate.cap.mode <- cfg$candidate.cap.mode %||% "deterministic"
      vertex.records[[idx]] <- vt
    if (nrow(op$path.table) > 0L) {
      pt <- op$path.table
      pt$run.id <- cfg$run.id
      pt$dataset <- cfg$dataset
      pt$n <- cfg$n
      pt$seed <- cfg$seed
      pt$graph.type <- cfg$graph.type
      pt$path.family <- cfg$path.family %||% "local.composite.geodesics"
      pt$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      pt$disk.radius.rule <- cfg$disk.radius.rule
      pt$row.normalization <- cfg$row.normalization
      pt$vertex.normalization <- cfg$vertex.normalization
      pt$row.selection <- cfg$row.selection
      pt$center.order <- cfg$center.order %||% "increasing.degree"
      pt$candidate.cap.mode <- cfg$candidate.cap.mode %||% "deterministic"
      path.records[[idx]] <- pt
    }
    if (!is.null(op$skeleton.table) && nrow(op$skeleton.table) > 0L) {
      st <- op$skeleton.table
      st$run.id <- cfg$run.id
      st$dataset <- cfg$dataset
      st$n <- cfg$n
      st$seed <- cfg$seed
      st$graph.type <- cfg$graph.type
      st$path.family <- cfg$path.family %||% "local.composite.geodesics"
      st$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      st$degree <- cfg$degree
      st$disk.k <- cfg$disk.k
      st$row.normalization <- cfg$row.normalization
      st$vertex.normalization <- cfg$vertex.normalization
      skeleton.records[[idx]] <- st
    }
    if (!is.null(op$skeleton.direction.table) &&
        nrow(op$skeleton.direction.table) > 0L) {
      sdt <- op$skeleton.direction.table
      sdt$run.id <- cfg$run.id
      sdt$dataset <- cfg$dataset
      sdt$n <- cfg$n
      sdt$seed <- cfg$seed
      sdt$graph.type <- cfg$graph.type
      sdt$path.family <- cfg$path.family %||% "local.composite.geodesics"
      sdt$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      sdt$degree <- cfg$degree
      sdt$disk.k <- cfg$disk.k
      sdt$row.normalization <- cfg$row.normalization
      sdt$vertex.normalization <- cfg$vertex.normalization
      skeleton.direction.records[[idx]] <- sdt
    }
    if (!is.null(op$skeleton.connector.candidate.table) &&
        nrow(op$skeleton.connector.candidate.table) > 0L) {
      sct <- op$skeleton.connector.candidate.table
      sct$run.id <- cfg$run.id
      sct$dataset <- cfg$dataset
      sct$n <- cfg$n
      sct$seed <- cfg$seed
      sct$graph.type <- cfg$graph.type
      sct$path.family <- cfg$path.family %||% "local.composite.geodesics"
      sct$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      sct$degree <- cfg$degree
      sct$disk.k <- cfg$disk.k
      sct$row.normalization <- cfg$row.normalization
      sct$vertex.normalization <- cfg$vertex.normalization
      skeleton.connector.candidate.records[[idx]] <- sct
    }
    if (!is.null(op$skeleton.coverage) && nrow(op$skeleton.coverage) > 0L) {
      sc <- op$skeleton.coverage
      sc$run.id <- cfg$run.id
      sc$dataset <- cfg$dataset
      sc$n <- cfg$n
      sc$seed <- cfg$seed
      sc$graph.type <- cfg$graph.type
      sc$path.family <- cfg$path.family %||% "local.composite.geodesics"
      sc$skeleton.cross.rule <- cfg$skeleton.cross.rule %||% "coverage"
      sc$degree <- cfg$degree
      sc$disk.k <- cfg$disk.k
      sc$row.normalization <- cfg$row.normalization
      sc$vertex.normalization <- cfg$vertex.normalization
      skeleton.coverage.records[[idx]] <- sc
    }
    if (nrow(op$local.dimension.table) > 0L) {
      ldt <- op$local.dimension.table
      ldt$run.id <- cfg$run.id
      ldt$dataset <- cfg$dataset
      ldt$n <- cfg$n
      ldt$seed <- cfg$seed
      ldt$graph.type <- cfg$graph.type
      ldt$path.family <- cfg$path.family %||% "local.composite.geodesics"
      ldt$degree <- cfg$degree
      ldt$disk.radius.rule <- cfg$disk.radius.rule
      ldt$row.normalization <- cfg$row.normalization
      ldt$vertex.normalization <- cfg$vertex.normalization
      ldt$row.selection <- cfg$row.selection
      ldt$center.order <- cfg$center.order %||% "increasing.degree"
      ldt$candidate.cap.mode <- cfg$candidate.cap.mode %||% "deterministic"
      local.dimension.records[[idx]] <- ldt
    }
    if (!is.null(op$local.dimension.vote.table) &&
        nrow(op$local.dimension.vote.table) > 0L) {
      ldv <- op$local.dimension.vote.table
      ldv$run.id <- cfg$run.id
      ldv$dataset <- cfg$dataset
      ldv$n <- cfg$n
      ldv$seed <- cfg$seed
      ldv$graph.type <- cfg$graph.type
      ldv$path.family <- cfg$path.family %||% "local.composite.geodesics"
      ldv$degree <- cfg$degree
      ldv$disk.radius.rule <- cfg$disk.radius.rule
      ldv$row.normalization <- cfg$row.normalization
      ldv$vertex.normalization <- cfg$vertex.normalization
      ldv$row.selection <- cfg$row.selection
      ldv$center.order <- cfg$center.order %||% "increasing.degree"
      ldv$candidate.cap.mode <- cfg$candidate.cap.mode %||% "deterministic"
      local.dimension.vote.records[[idx]] <- ldv
    }
    result.records[[idx]] <- cbind(cfg, data.frame(
      status = "ok",
      message = "",
      n.edges = length(graph$edge.weights),
      n.operator.rows = nrow(op$A),
      nullity.1e.6 = nullity$table$nullity[nullity$table$rel.tol == 1e-6],
      expected.nullity = choose(cfg$degree + 2L, 2L),
      intrinsic.residual = poly.scores$residual.ratio[
        poly.scores$target == "intrinsic.degree.leq.d"
      ],
      control.residual = poly.scores$residual.ratio[
        poly.scores$target == "intrinsic.exact.degree.d.plus.1"
      ],
      stringsAsFactors = FALSE
    ))
    result.records[[idx]]$score.table <- list(poly.scores)
  }

  results <- do.call(rbind, result.records)
  nullity.table <- if (length(nullity.records)) do.call(rbind, nullity.records) else data.frame()
  singular.table <- if (length(singular.records)) do.call(rbind, singular.records) else data.frame()
  vertex.table <- if (length(vertex.records)) do.call(rbind, vertex.records) else data.frame()
  path.table <- if (length(path.records)) do.call(rbind, path.records) else data.frame()
  skeleton.table <- if (length(skeleton.records)) do.call(rbind, skeleton.records) else data.frame()
  skeleton.coverage.table <- if (length(skeleton.coverage.records)) {
    do.call(rbind, skeleton.coverage.records)
  } else {
    data.frame()
  }
  skeleton.direction.table <- if (length(skeleton.direction.records)) {
    do.call(rbind, skeleton.direction.records)
  } else {
    data.frame()
  }
  skeleton.connector.candidate.table <- if (length(skeleton.connector.candidate.records)) {
    do.call(rbind, skeleton.connector.candidate.records)
  } else {
    data.frame()
  }
  local.dimension.table <- if (length(local.dimension.records)) {
    do.call(rbind, local.dimension.records)
  } else {
    data.frame()
  }
  local.dimension.vote.table <- if (length(local.dimension.vote.records)) {
    do.call(rbind, local.dimension.vote.records)
  } else {
    data.frame()
  }
  score.table <- do.call(rbind, lapply(seq_len(nrow(results)), function(i) {
    score <- results$score.table[[i]]
    if (is.null(score) || !is.data.frame(score)) return(NULL)
    base <- results[i, setdiff(names(results), "score.table"), drop = FALSE]
    row.names(base) <- NULL
    row.names(score) <- NULL
    cbind(base, score)
  }))
  results$score.table <- NULL

  diagnostics <- list(
    mode = mode,
    config = config,
    results = results,
    nullity = nullity.table,
    singular.values = singular.table,
    vertex.table = vertex.table,
    path.table = path.table,
    skeleton.table = skeleton.table,
    skeleton.coverage.table = skeleton.coverage.table,
    skeleton.direction.table = skeleton.direction.table,
    skeleton.connector.candidate.table = skeleton.connector.candidate.table,
    local.dimension.table = local.dimension.table,
    local.dimension.vote.table = local.dimension.vote.table,
    score.table = score.table,
    created = Sys.time()
  )
  saveRDS(diagnostics, file.path(output.dir,
                                 paste0("geodesic_annihilator_", mode, ".rds")))
  utils::write.csv(results, file.path(output.dir,
                                      paste0("geodesic_annihilator_", mode, "_summary.csv")),
                   row.names = FALSE)
  utils::write.csv(nullity.table, file.path(output.dir,
                                            paste0("geodesic_annihilator_", mode, "_nullity.csv")),
                   row.names = FALSE)
  utils::write.csv(vertex.table, file.path(output.dir,
                                           paste0("geodesic_annihilator_", mode, "_vertex_diagnostics.csv")),
                   row.names = FALSE)
  utils::write.csv(path.table, file.path(output.dir,
                                         paste0("geodesic_annihilator_", mode, "_path_table.csv")),
                   row.names = FALSE)
  utils::write.csv(skeleton.table, file.path(output.dir,
                                             paste0("geodesic_annihilator_", mode, "_skeleton_table.csv")),
                   row.names = FALSE)
  utils::write.csv(skeleton.coverage.table, file.path(output.dir,
                                                      paste0("geodesic_annihilator_", mode, "_skeleton_coverage.csv")),
                   row.names = FALSE)
  utils::write.csv(skeleton.direction.table, file.path(output.dir,
                                                       paste0("geodesic_annihilator_", mode, "_skeleton_direction_diagnostics.csv")),
                   row.names = FALSE)
  utils::write.csv(skeleton.connector.candidate.table, file.path(output.dir,
                                                                paste0("geodesic_annihilator_", mode, "_skeleton_connector_candidates.csv")),
                   row.names = FALSE)
  utils::write.csv(local.dimension.table, file.path(output.dir,
                                                    paste0("geodesic_annihilator_", mode, "_local_dimension.csv")),
                   row.names = FALSE)
  utils::write.csv(local.dimension.vote.table, file.path(output.dir,
                                                         paste0("geodesic_annihilator_", mode, "_local_dimension_votes.csv")),
                   row.names = FALSE)
  utils::write.csv(score.table, file.path(output.dir,
                                          paste0("geodesic_annihilator_", mode, "_polynomial_scores.csv")),
                   row.names = FALSE)
  diagnostics
}

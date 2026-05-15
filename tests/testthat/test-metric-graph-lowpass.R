make_path_graph_lengths <- function(lengths) {
  n <- length(lengths) + 1L
  adj <- vector("list", n)
  w <- vector("list", n)
  add_edge <- function(i, j, ell) {
    adj[[i]] <<- c(adj[[i]], j)
    w[[i]] <<- c(w[[i]], ell)
    adj[[j]] <<- c(adj[[j]], i)
    w[[j]] <<- c(w[[j]], ell)
  }
  for (i in seq_along(lengths)) add_edge(i, i + 1L, lengths[i])
  list(adj.list = lapply(adj, as.integer), weight.list = lapply(w, as.double))
}

laplacian_dense <- function(op) {
  as.matrix(op$laplacian$matrix)
}

test_that("metric operator validates undirected positive length graphs", {
  graph <- make_path_graph_lengths(c(1, 2))
  graph.bad <- graph
  graph.bad$adj.list[[2]] <- as.integer(1L)
  graph.bad$weight.list[[2]] <- 1

  expect_error(
    metric.graph.lowpass.operator(graph.bad$adj.list, graph.bad$weight.list),
    "no reciprocal entry"
  )

  graph.bad <- graph
  graph.bad$weight.list[[1]][1] <- -1
  expect_error(
    metric.graph.lowpass.operator(graph.bad$adj.list, graph.bad$weight.list),
    "non-positive"
  )
})

test_that("metric conductance transforms match hand calculations", {
  graph <- make_path_graph_lengths(c(2, 4))
  eps <- 1e-8

  op.inv <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "inverse.length.power",
    conductance.epsilon = eps,
    conductance.alpha = 2
  )
  expect_equal(op.inv$edge.table$conductance, (c(2, 4) + eps)^(-2), tolerance = 1e-12)

  op.exp <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "exp.length",
    conductance.sigma = 2
  )
  expect_equal(op.exp$edge.table$conductance, exp(-c(2, 4) / 2), tolerance = 1e-12)

  op.exp2 <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "exp.length.squared",
    conductance.sigma = 2
  )
  expect_equal(op.exp2$edge.table$conductance, exp(-(c(2, 4)^2) / 4), tolerance = 1e-12)
})

test_that("self-tuned Gaussian uses local incident edge scales", {
  graph <- make_path_graph_lengths(c(2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "self.tuned.gaussian",
    conductance.local.k = 1L
  )
  expect_equal(op$conductance$local.scales, c(2, 2, 4), tolerance = 1e-12)
  expect_equal(
    op$edge.table$conductance,
    c(exp(-(2^2) / (2 * 2)), exp(-(4^2) / (2 * 4))),
    tolerance = 1e-12
  )
})

test_that("unnormalized weighted Laplacian assembly matches R reference", {
  graph <- make_path_graph_lengths(c(2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "inverse.length.power",
    conductance.epsilon = 1e-8
  )
  c12 <- 1 / (2 + 1e-8)
  c23 <- 1 / (4 + 1e-8)
  L.ref <- matrix(c(
     c12, -c12,     0,
    -c12, c12+c23, -c23,
        0, -c23,    c23
  ), nrow = 3, byrow = TRUE)
  L <- laplacian_dense(op)
  expect_equal(L, L.ref, tolerance = 1e-12)
  expect_equal(rowSums(L), rep(0, 3), tolerance = 1e-12)
  expect_true(isTRUE(all.equal(L, t(L), tolerance = 1e-12)))
  expect_true(min(eigen(L, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
})

test_that("C++ dense spectral path matches R dense eigen reference", {
  graph <- make_path_graph_lengths(c(1, 2, 3, 4))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list,
    y = sin(seq_len(5)),
    n.eigenpairs = 5L,
    eigen.solver = "dense",
    conductance.rule = "inverse.length.power",
    conductance.epsilon = 1e-8
  )
  L <- laplacian_dense(fit$operator)
  ev <- eigen(L, symmetric = TRUE)
  expect_equal(fit$spectral$eigenvalues, sort(ev$values), tolerance = 1e-9)

  U.ref <- ev$vectors[, order(ev$values), drop = FALSE]
  cors <- abs(diag(crossprod(fit$spectral$eigenvectors, U.ref)))
  expect_true(all(cors > 1 - 1e-7))
})

test_that("sparse spectral path handles inverse-square metric conductances", {
  lengths <- exp(seq(log(1e-4), log(2e-2), length.out = 79))
  graph <- make_path_graph_lengths(lengths)
  y <- sin(seq_len(80) / 7)
  eta.grid <- c(1e-8, 1e-6, 1e-4)

  fit.sparse <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "inverse.length.power",
    conductance.alpha = 2,
    n.eigenpairs = 20L,
    eigen.solver = "sparse",
    dense.fallback = "never",
    eta.grid = eta.grid
  )
  fit.dense <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "inverse.length.power",
    conductance.alpha = 2,
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )

  expect_equal(fit.sparse$spectral$backend, "sparse.shift")
  expect_equal(fit.sparse$spectral$eigenvalues, fit.dense$spectral$eigenvalues, tolerance = 1e-6)
  expect_equal(fit.sparse$fitted.values, fit.dense$fitted.values, tolerance = 1e-8)
})

test_that("metric low-pass preserves constant responses and refit fixed eta", {
  graph <- make_path_graph_lengths(rep(1, 4))
  y <- rep(3, 5)
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    n.eigenpairs = 5L,
    eigen.solver = "dense",
    eta.grid = c(0.1, 1, 10)
  )
  expect_equal(fit$fitted.values, y, tolerance = 1e-10)

  y2 <- seq_len(5)
  refit <- refit.metric.graph.lowpass(fit, y2)
  V <- fit$spectral$eigenvectors
  f <- fit$spectral$filtered.eigenvalues
  y2.hat <- as.vector(V %*% (f * as.vector(crossprod(V, y2))))
  expect_equal(refit$fitted.values, y2.hat, tolerance = 1e-10)
})

test_that("effective degrees of freedom decrease with eta for heat kernel", {
  graph <- make_path_graph_lengths(rep(1, 5))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list,
    y = sin(seq_len(6)),
    n.eigenpairs = 6L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 0.1, 1, 10),
    filter.type = "heat_kernel"
  )
  weights <- compute.filter.weights.matrix(
    fit$spectral$eigenvalues,
    fit$gcv$eta.grid,
    fit$spectral$filter.type
  )
  expect_true(all(diff(colSums(weights)) <= 1e-12))
})

test_that("per-column GCV refit matches independent fits", {
  graph <- make_path_graph_lengths(rep(1, 7))
  Y <- cbind(
    sin(seq_len(8) / 2),
    cos(seq_len(8) / 3)
  )
  eta.grid <- c(0.01, 0.1, 1, 10)
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 1],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )
  refit <- refit.metric.graph.lowpass(
    fit, Y,
    per.column.gcv = TRUE,
    eta.grid = eta.grid
  )

  fit1 <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 1],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )
  fit2 <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 2],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )

  expect_equal(refit$eta.optimal, c(fit1$gcv$eta.optimal, fit2$gcv$eta.optimal), tolerance = 1e-12)
  expect_equal(refit$fitted.values[, 1], fit1$fitted.values, tolerance = 1e-10)
  expect_equal(refit$fitted.values[, 2], fit2$fitted.values, tolerance = 1e-10)
})

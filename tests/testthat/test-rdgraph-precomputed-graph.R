as_edge_table_undirected <- function(adj.list, weight.list) {
  expect_true(is.list(adj.list))
  expect_true(is.list(weight.list))
  expect_identical(length(adj.list), length(weight.list))

  u <- integer(0)
  v <- integer(0)
  w <- numeric(0)

  for (i in seq_along(adj.list)) {
    nbrs <- as.integer(adj.list[[i]])
    wts <- as.double(weight.list[[i]])
    expect_identical(length(nbrs), length(wts))

    keep <- nbrs > i
    if (!any(keep)) next

    u <- c(u, rep.int(i, sum(keep)))
    v <- c(v, nbrs[keep])
    w <- c(w, wts[keep])
  }

  if (!length(u)) {
    return(data.frame(u = integer(0), v = integer(0), w = numeric(0)))
  }

  ord <- order(u, v)
  data.frame(u = u[ord], v = v[ord], w = w[ord])
}

test_that("fit accepts precomputed graph and matches standard fit-path result", {
  set.seed(2611)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))
  k <- 9L

  iknn <- create.single.iknn.graph(
    X,
    k = k,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    compute.full = FALSE,
    verbose = FALSE
  )

  fit.standard <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  fit.injected <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    adj.list = iknn$pruned_adj_list,
    weight.list = iknn$pruned_weight_list,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  expect_identical(fit.injected$parameters$graph.source, "precomputed")
  expect_identical(fit.standard$graph$n.edges, fit.injected$graph$n.edges)

  edges.standard <- as_edge_table_undirected(
    fit.standard$graph$adj.list,
    fit.standard$graph$edge.length.list
  )
  edges.injected <- as_edge_table_undirected(
    fit.injected$graph$adj.list,
    fit.injected$graph$edge.length.list
  )

  expect_equal(edges.standard[, c("u", "v")], edges.injected[, c("u", "v")])
  expect_equal(edges.standard$w, edges.injected$w, tolerance = 1e-12)
  expect_identical(length(fit.standard$fitted.values), length(fit.injected$fitted.values))
  expect_true(all(is.finite(fit.injected$fitted.values)))
  expect_true(all(is.finite(fit.injected$residuals)))
})

test_that("fit rejects malformed precomputed graph inputs", {
  set.seed(2612)
  X <- matrix(rnorm(80 * 6), nrow = 80, ncol = 6)
  y <- rnorm(nrow(X))
  k <- 8L

  iknn <- create.single.iknn.graph(
    X,
    k = k,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    compute.full = FALSE,
    verbose = FALSE
  )

  expect_error(
    fit.rdgraph.regression(
      X,
      y,
      k = k,
      adj.list = iknn$pruned_adj_list,
      weight.list = NULL,
      max.iterations = 1L,
      n.eigenpairs = 40L,
      pca.dim = NULL,
      verbose.level = 0L
    ),
    "provided together"
  )

  adj.bad <- iknn$pruned_adj_list
  w.bad <- iknn$pruned_weight_list
  i <- which(lengths(adj.bad) > 0L)[1L]
  j <- adj.bad[[i]][1L]
  rev.idx <- match(i, adj.bad[[j]])
  adj.bad[[j]] <- adj.bad[[j]][-rev.idx]
  w.bad[[j]] <- w.bad[[j]][-rev.idx]

  expect_error(
    fit.rdgraph.regression(
      X,
      y,
      k = k,
      adj.list = adj.bad,
      weight.list = w.bad,
      max.iterations = 1L,
      n.eigenpairs = 40L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      max.ratio.threshold = 0,
      threshold.percentile = 0,
      dense.fallback = "never",
      verbose.level = 0L
    ),
    "undirected|reciprocal"
  )
})

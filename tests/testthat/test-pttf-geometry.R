make_pttf_path_graph <- function(n) {
  adj <- vector("list", n)
  w <- vector("list", n)
  x <- seq(0, 1, length.out = n)
  for (i in seq_len(n - 1L)) {
    ell <- x[i + 1L] - x[i]
    adj[[i]] <- c(adj[[i]], i + 1L)
    w[[i]] <- c(w[[i]], ell)
    adj[[i + 1L]] <- c(adj[[i + 1L]], i)
    w[[i + 1L]] <- c(w[[i + 1L]], ell)
  }
  list(
    X = matrix(x, ncol = 1L),
    adj.list = lapply(adj, as.integer),
    weight.list = lapply(w, as.double)
  )
}

test_that("pttf.geometry validates fixed tangent dimension and supplied graphs", {
  graph <- make_pttf_path_graph(12)

  expect_error(
    pttf.geometry(graph$X, graph = "supplied", tangent.dim = 2L,
                  adj.list = graph$adj.list, weight.list = graph$weight.list),
    "tangent.dim"
  )
  expect_error(
    pttf.geometry(graph$X, graph = "supplied", tangent.dim = 1L,
                  adj.list = graph$adj.list),
    "both 'adj.list' and 'weight.list'"
  )
  expect_error(
    pttf.geometry(graph$X, graph = "supplied", tangent.dim = 1L,
                  adj.list = graph$adj.list, weight.list = graph$weight.list,
                  min.support = 2L),
    "min.support"
  )

  geom <- pttf.geometry(
    graph$X,
    graph = "supplied",
    adj.list = graph$adj.list,
    weight.list = graph$weight.list,
    tangent.dim = 1L
  )
  expect_s3_class(geom, "pttf_geometry")
  expect_equal(geom$graph$constructor, "supplied")
  expect_equal(geom$frames$tangent.dim, 1L)
  expect_equal(geom$parameters$min.support, 4L)
  expect_true(all(geom$frames$support.diagnostics$status == "ok"))
  expect_equal(geom$transport$convention,
               "Oij maps frame j coordinates to frame i coordinates")
  expect_false(geom$density$used.in.operator)
})

test_that("pttf.geometry returns oriented transports with transpose convention", {
  set.seed(20260526)
  X <- matrix(runif(120), ncol = 2L)
  geom <- pttf.geometry(
    X,
    tangent.dim = 2L,
    graph.k.scale = 4L,
    max.hops = 3L,
    min.support = 8L
  )

  expect_true(any(geom$edges$status == "ok"))
  ok.edge <- geom$edges[geom$edges$status == "ok", ][1L, ]
  key.ij <- paste(ok.edge$from, ok.edge$to, sep = "->")
  key.ji <- paste(ok.edge$to, ok.edge$from, sep = "->")
  Oij <- geom$transport$Oij[[key.ij]]
  Oji <- geom$transport$Oij[[key.ji]]
  expect_equal(Oji, t(Oij), tolerance = 1e-12)
  expect_equal(as.matrix(crossprod(Oij)), diag(2), tolerance = 1e-10)
  expect_true(is.finite(ok.edge$procrustes.residual))

  shared <- intersect(
    geom$frames$support.index.list[[ok.edge$from]],
    geom$frames$support.index.list[[ok.edge$to]]
  )
  Fi <- geom$frames$frame.list[[ok.edge$from]]
  Fj <- geom$frames$frame.list[[ok.edge$to]]
  Zi <- sweep(X[shared, , drop = FALSE], 2L, X[ok.edge$from, ], "-") %*% Fi
  Zj <- sweep(X[shared, , drop = FALSE], 2L, X[ok.edge$to, ], "-") %*% Fj
  residual <- norm(Zj %*% t(Oij) - Zi, type = "F") /
    max(norm(Zi, type = "F"), .Machine$double.eps)
  expect_equal(ok.edge$procrustes.residual, residual, tolerance = 1e-12)

  frame.svd <- svd(t(Fi) %*% Fj)
  frame.map <- frame.svd$u %*% t(frame.svd$v)
  expect_equal(
    ok.edge$frame.transport.residual,
    norm(frame.map - Oij, type = "F") / sqrt(2),
    tolerance = 1e-12
  )
  expect_true("frame.inverse.residual" %in% names(geom$edges))
})

test_that("pttf.geometry rknn mode records graph, density, and cycle diagnostics", {
  set.seed(20260527)
  theta <- seq(0, 2 * pi, length.out = 40L)
  X <- cbind(cos(theta), sin(theta))

  geom <- pttf.geometry(
    X,
    tangent.dim = 1L,
    graph = "rknn",
    graph.k.scale = 3L,
    max.hops = 3L,
    diagnostics = TRUE
  )

  expect_equal(geom$graph$constructor, "create.rknn.graph")
  expect_true(nrow(geom$graph$edge.table) > 0L)
  expect_true(all(is.finite(geom$density$rho.hat)))
  expect_true(all(is.finite(geom$density$node.mass)))
  expect_true(all(is.finite(geom$density$proposed.alpha.weight)))
  expect_true("support" %in% names(geom$diagnostics$status.summary))
  expect_true("edge" %in% names(geom$diagnostics$status.summary))
  expect_true(is.list(geom$diagnostics$cycles))
  expect_equal(geom$diagnostics$orientation.synchronization$method, "bfs.sign")
  expect_true(isTRUE(geom$diagnostics$orientation.synchronization$metadata.only))
  expect_equal(geom$diagnostics$orientation.synchronization$n.reached, nrow(X))
  expect_true(is.numeric(geom$diagnostics$orientation.synchronization$n.flips))
  expect_true(
    all(c("n.edges", "median.residual", "max.residual") %in%
          names(geom$diagnostics$orientation.synchronization$non.tree.edge.summary))
  )
})

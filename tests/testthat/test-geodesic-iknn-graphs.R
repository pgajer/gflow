canonicalize_geodesic_edges <- function(graph) {
  u <- integer(0)
  v <- integer(0)
  w <- numeric(0)
  isize <- integer(0)

  for (i in seq_along(graph$adj_list)) {
    nbrs <- as.integer(graph$adj_list[[i]])
    weights <- as.double(graph$weight_list[[i]])
    isizes <- if (!is.null(graph$isize_list)) {
      as.integer(graph$isize_list[[i]])
    } else {
      rep.int(NA_integer_, length(nbrs))
    }

    keep <- nbrs > i
    if (!any(keep)) {
      next
    }

    u <- c(u, rep.int(i, sum(keep)))
    v <- c(v, nbrs[keep])
    w <- c(w, weights[keep])
    isize <- c(isize, isizes[keep])
  }

  if (!length(u)) {
    return(data.frame(u = integer(0), v = integer(0), w = numeric(0), isize = integer(0)))
  }

  ord <- order(u, v)
  data.frame(u = u[ord], v = v[ord], w = w[ord], isize = isize[ord])
}

test_that("create.geodesic.iknn.graph applies the nerve rule with graph distances", {
  path.graph <- list(
    adj_list = list(2L, c(1L, 3L), 2L),
    weight_list = list(1, c(1, 1), 1)
  )

  g1 <- create.geodesic.iknn.graph(path.graph, k = 1L)
  edges <- canonicalize_geodesic_edges(g1)

  expect_s3_class(g1, "geodesic_iknn_graph")
  expect_equal(
    edges[, c("u", "v", "w")],
    data.frame(
      u = c(1L, 1L, 2L),
      v = c(2L, 3L, 3L),
      w = c(1, 2, 1)
    ),
    tolerance = 1e-12
  )
  expect_equal(edges$isize, c(2L, 1L, 1L))
  expect_equal(g1$n_edges, 3L)
})

test_that("create.iterated.iknn.graphs starts from create.iknn.graphs output", {
  set.seed(20260423)
  X <- matrix(rnorm(30), ncol = 2)

  iterated <- create.iterated.iknn.graphs(
    X,
    kmin = 2L,
    kmax = 2L,
    n.iterations = 2L,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    pca.dim = NULL,
    n.cores = 1L,
    verbose = FALSE
  )

  g0.reference <- create.iknn.graphs(
    X,
    kmin = 2L,
    kmax = 2L,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    compute.full = TRUE,
    pca.dim = NULL,
    n.cores = 1L,
    verbose = FALSE
  )$geom_pruned_graphs[[1L]]

  expect_s3_class(iterated, "iterated_iknn_graphs")
  expect_named(iterated$graphs, c("G0", "G1", "G2"))
  expect_equal(
    canonicalize_geodesic_edges(iterated$graphs$G0[["2"]])[, c("u", "v", "w")],
    canonicalize_geodesic_edges(g0.reference)[, c("u", "v", "w")],
    tolerance = 1e-12
  )
  expect_equal(nrow(iterated$summary), 3L)
  expect_equal(iterated$summary$graph, c("G0", "G1", "G2"))
})

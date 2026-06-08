.edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}

.circle_dist <- function(n) {
  theta <- 2 * pi * (seq_len(n) - 1L) / n
  X <- cbind(cos(theta), sin(theta))
  as.matrix(stats::dist(X))
}

.adj_edge_keys <- function(adj.list) {
  edges <- do.call(rbind, lapply(seq_along(adj.list), function(i) {
    nbrs <- adj.list[[i]]
    nbrs <- nbrs[nbrs > i]
    if (length(nbrs) == 0L) {
      return(NULL)
    }
    cbind(i, nbrs)
  }))
  .edge_keys(edges)
}


test_that("distance-matrix mKNN helper matches mutual neighbor definition", {
  D <- matrix(c(
    0, 1, 2, 4,
    1, 0, 1, 3,
    2, 1, 0, 1,
    4, 3, 1, 0
  ), nrow = 4, byrow = TRUE)

  g <- gflow:::.dgh_mknn_graph_from_dist(D, k = 2)

  expect_equal(.edge_keys(g$edge.matrix), c("1-2", "2-3", "3-4"))
  expect_equal(g$edge.weight, c(1, 1, 1))
  expect_equal(g$knn.sets[[1]], c(1L, 2L, 3L))
  expect_equal(g$adj.list[[2]], c(1L, 3L))
  expect_equal(g$weight.list[[2]], c(1, 1))
  expect_equal(g$weight.type, "distance")
})


test_that("open and closed mKNN neighborhoods give the same edge predicate", {
  D <- matrix(c(
    0, 1, 1.2, 3.5, 4.0,
    1, 0, 1.1, 2.8, 4.2,
    1.2, 1.1, 0, 1.3, 3.9,
    3.5, 2.8, 1.3, 0, 1.0,
    4.0, 4.2, 3.9, 1.0, 0
  ), nrow = 5, byrow = TRUE)

  open.nn <- gflow:::.dgh_knn_sets_from_dist(D, k = 2)
  closed.nn <- gflow:::.dgh_closed_knn_sets_from_dist(D, k = 2)

  open.edges <- matrix(integer(0), ncol = 2L)
  closed.edges <- matrix(integer(0), ncol = 2L)
  for (i in 1:4) {
    for (j in (i + 1):5) {
      if (j %in% open.nn[[i]] && i %in% open.nn[[j]]) {
        open.edges <- rbind(open.edges, c(i, j))
      }
      if (j %in% closed.nn[[i]] && i %in% closed.nn[[j]]) {
        closed.edges <- rbind(closed.edges, c(i, j))
      }
    }
  }

  g <- gflow:::.dgh_mknn_graph_from_dist(D, k = 2)
  expect_equal(.edge_keys(open.edges), .edge_keys(closed.edges))
  expect_equal(.edge_keys(g$edge.matrix), .edge_keys(closed.edges))
})


test_that("distance-matrix ikNN helper matches closed shared-neighbor definition and detour weights", {
  D <- matrix(c(
    0, 1, 1, 3,
    1, 0, 2, 2,
    1, 2, 0, 2,
    3, 2, 2, 0
  ), nrow = 4, byrow = TRUE)

  g <- gflow:::.dgh_iknn_graph_from_dist(D, k = 2)

  expect_equal(.edge_keys(g$edge.matrix), c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4"))
  expect_equal(g$intersection.size, c(3L, 3L, 2L, 3L, 2L, 2L))
  expect_equal(g$direct.distance, c(1, 1, 3, 2, 2, 2))
  expect_equal(g$edge.weight, c(1, 1, 3, 2, 2, 2))
  expect_equal(g$weight.type, "closed_shared_neighbor_detour")
})


test_that("distance-matrix PHATE support helper matches adaptive OR radius", {
  D <- matrix(c(
    0, 1.0, 1.8, 5.0,
    1.0, 0, 1.2, 2.1,
    1.8, 1.2, 0, 1.0,
    5.0, 2.1, 1.0, 0
  ), nrow = 4, byrow = TRUE)

  g <- gflow:::.dgh_phate_support_graph_from_dist(
    D,
    k = 1,
    decay = 40,
    threshold = 1e-4
  )

  r.tau <- (-log(1e-4))^(1 / 40)
  expected <- matrix(integer(0), ncol = 2L)
  for (i in 1:3) {
    for (j in (i + 1):4) {
      if (D[i, j] <= r.tau * max(g$sigma[i], g$sigma[j])) {
        expected <- rbind(expected, c(i, j))
      }
    }
  }

  expect_equal(.edge_keys(g$edge.matrix), .edge_keys(expected))
  expect_true(all(g$K[upper.tri(g$K)][g$K[upper.tri(g$K)] > 0] > 0))
  expect_equal(rowSums(g$P), rep(1, 4), tolerance = 1e-12)
  expect_equal(diag(g$A), rep(1, 4), tolerance = 1e-12)
  expect_equal(g$weight.type, "symmetrized_affinity")
})


test_that("root-of-unity examples match the graph-construction definitions", {
  D3 <- .circle_dist(3)
  ph3 <- gflow:::.dgh_phate_support_graph_from_dist(D3, k = 2)
  mk3 <- gflow:::.dgh_mknn_graph_from_dist(D3, k = 2)
  ik3 <- gflow:::.dgh_iknn_graph_from_dist(D3, k = 2)

  complete3 <- c("1-2", "1-3", "2-3")
  expect_equal(.edge_keys(ph3$edge.matrix), complete3)
  expect_equal(.edge_keys(mk3$edge.matrix), complete3)
  expect_equal(.edge_keys(ik3$edge.matrix), complete3)

  D4 <- .circle_dist(4)
  ph4 <- gflow:::.dgh_phate_support_graph_from_dist(D4, k = 2)
  mk4 <- gflow:::.dgh_mknn_graph_from_dist(D4, k = 2)
  ik4 <- gflow:::.dgh_iknn_graph_from_dist(D4, k = 2)

  cycle4 <- c("1-2", "1-4", "2-3", "3-4")
  complete4 <- c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4")
  expect_equal(.edge_keys(ph4$edge.matrix), cycle4)
  expect_equal(.edge_keys(mk4$edge.matrix), cycle4)
  expect_equal(.edge_keys(ik4$edge.matrix), complete4)

  D5 <- .circle_dist(5)
  ph5 <- gflow:::.dgh_phate_support_graph_from_dist(D5, k = 2)
  mk5 <- gflow:::.dgh_mknn_graph_from_dist(D5, k = 2)
  ik5 <- gflow:::.dgh_iknn_graph_from_dist(D5, k = 2)

  cycle5 <- c("1-2", "1-5", "2-3", "3-4", "4-5")
  complete5 <- c("1-2", "1-3", "1-4", "1-5", "2-3",
                 "2-4", "2-5", "3-4", "3-5", "4-5")

  expect_equal(.edge_keys(ph5$edge.matrix), cycle5)
  expect_equal(.edge_keys(mk5$edge.matrix), cycle5)
  expect_equal(.edge_keys(ik5$edge.matrix), complete5)
})

test_that("distance-matrix ikNN helper matches native gflow ikNN on root examples", {
  roots <- function(n) {
    theta <- 2 * pi * (seq_len(n) - 1L) / n
    cbind(cos(theta), sin(theta))
  }

  for (n in c(4L, 5L)) {
    X <- roots(n)
    helper <- gflow:::.dgh_iknn_graph_from_dist(as.matrix(stats::dist(X)), k = 2)
    native <- dgraphs::create.single.iknn.graph(
      X,
      k = 2,
      max.path.edge.ratio.deviation.thld = 0,
      threshold.percentile = 0,
      compute.full = TRUE,
      pca.dim = NULL,
      verbose = FALSE
    )

    expect_equal(.edge_keys(helper$edge.matrix), .adj_edge_keys(native$adj_list))
  }
})


test_that("mKNN support is contained in PHATE support under standard assumptions", {
  set.seed(42)
  X <- matrix(stats::rnorm(24), ncol = 3)
  D <- as.matrix(stats::dist(X))

  ph <- gflow:::.dgh_phate_support_graph_from_dist(D, k = 3, decay = 40,
                                                   threshold = 1e-4)
  mk <- gflow:::.dgh_mknn_graph_from_dist(D, k = 3)

  expect_true(all(.edge_keys(mk$edge.matrix) %in% .edge_keys(ph$edge.matrix)))
})

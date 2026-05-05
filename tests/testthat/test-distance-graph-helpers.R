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
  expect_equal(g$adj.list[[2]], c(1L, 3L))
  expect_equal(g$weight.list[[2]], c(1, 1))
  expect_equal(g$weight.type, "distance")
})


test_that("distance-matrix ikNN helper matches shared-neighbor definition and detour weights", {
  D <- matrix(c(
    0, 1, 1, 3,
    1, 0, 2, 2,
    1, 2, 0, 2,
    3, 2, 2, 0
  ), nrow = 4, byrow = TRUE)

  g <- gflow:::.dgh_iknn_graph_from_dist(D, k = 2)

  expect_equal(.edge_keys(g$edge.matrix), c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4"))
  expect_equal(g$intersection.size, c(1L, 1L, 2L, 1L, 1L, 1L))
  expect_equal(g$direct.distance, c(1, 1, 3, 2, 2, 2))
  expect_equal(g$edge.weight, c(3, 3, 3, 2, 4, 4))
  expect_equal(g$weight.type, "shared_neighbor_detour")
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

  D5 <- .circle_dist(5)
  ph5 <- gflow:::.dgh_phate_support_graph_from_dist(D5, k = 2)
  mk5 <- gflow:::.dgh_mknn_graph_from_dist(D5, k = 2)
  ik5 <- gflow:::.dgh_iknn_graph_from_dist(D5, k = 2)

  cycle5 <- c("1-2", "1-5", "2-3", "3-4", "4-5")
  pentagram5 <- c("1-3", "1-4", "2-4", "2-5", "3-5")

  expect_equal(.edge_keys(ph5$edge.matrix), cycle5)
  expect_equal(.edge_keys(mk5$edge.matrix), cycle5)
  expect_equal(.edge_keys(ik5$edge.matrix), pentagram5)
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

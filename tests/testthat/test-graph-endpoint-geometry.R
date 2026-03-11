test_that("compute.graph.endpoint.scores assigns highest score to chain tips", {
  n <- 21L
  adj.list <- lapply(seq_len(n), function(i) {
    nbrs <- integer(0)
    if (i > 1L) nbrs <- c(nbrs, i - 1L)
    if (i < n) nbrs <- c(nbrs, i + 1L)
    nbrs
  })
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  layout.3d <- cbind(seq_len(n), 0, 0)

  res <- compute.graph.endpoint.scores(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    k = c(5L, 7L, 9L),
    min.neighborhood.size = 3L
  )

  expect_equal(res$score[c(1, n)], c(1, 1), tolerance = 1e-8)
  expect_true(all(abs(res$score[2:(n - 1)]) < 1e-8))
})

test_that("C++ endpoint scorer matches the reference R implementation", {
  n <- 17L
  adj.list <- lapply(seq_len(n), function(i) integer(0))
  for (i in 1:(n - 1L)) {
    adj.list[[i]] <- c(adj.list[[i]], i + 1L)
    adj.list[[i + 1L]] <- c(adj.list[[i + 1L]], i)
  }
  adj.list[[5]] <- c(adj.list[[5]], 12L)
  adj.list[[12]] <- c(adj.list[[12]], 5L)

  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  layout.3d <- cbind(seq_len(n), sin(seq_len(n) / 2), cos(seq_len(n) / 3))

  ref <- gflow:::.compute.graph.endpoint.scores.reference(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    scales = c(4, 6),
    neighborhood = "geodesic_k",
    q = 0.1,
    neighbor.weighting = "uniform",
    gaussian.sigma = NULL,
    min.neighborhood.size = 3L,
    prefer.cpp = FALSE
  )

  cpp <- gflow:::.compute.graph.endpoint.scores.reference(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    scales = c(4, 6),
    neighborhood = "geodesic_k",
    q = 0.1,
    neighbor.weighting = "uniform",
    gaussian.sigma = NULL,
    min.neighborhood.size = 3L,
    prefer.cpp = TRUE
  )

  expect_equal(cpp$s.min, ref$s.min, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(cpp$s.q, ref$s.q, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(cpp$m, ref$m, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(cpp$score, ref$score, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(cpp$neighborhood.size, ref$neighborhood.size, ignore_attr = TRUE)
  expect_equal(cpp$distance.scale, ref$distance.scale, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("C++ scale-stability support matches the igraph reference", {
  adj.list <- list(
    c(2L),
    c(1L, 3L, 6L),
    c(2L, 4L),
    c(3L, 5L),
    c(4L),
    c(2L, 7L),
    c(6L, 8L),
    c(7L)
  )
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  local.max.by.scale <- matrix(FALSE, nrow = length(adj.list), ncol = 3L)
  colnames(local.max.by.scale) <- c("s1", "s2", "s3")
  local.max.by.scale[c(1L, 5L), 1L] <- TRUE
  local.max.by.scale[c(5L, 8L), 2L] <- TRUE
  local.max.by.scale[c(1L, 8L), 3L] <- TRUE

  ref <- gflow:::.compute.graph.endpoint.support.by.scale(
    adj.list = adj.list,
    weight.list = weight.list,
    local.max.by.scale = local.max.by.scale,
    radius = 2,
    prefer.cpp = FALSE
  )
  cpp <- gflow:::.compute.graph.endpoint.support.by.scale(
    adj.list = adj.list,
    weight.list = weight.list,
    local.max.by.scale = local.max.by.scale,
    radius = 2,
    prefer.cpp = TRUE
  )

  expect_equal(cpp, ref, ignore_attr = TRUE)
})

test_that("C++ greedy maxima suppression matches the reference fallback", {
  adj.list <- list(
    c(2L),
    c(1L, 3L, 6L),
    c(2L, 4L),
    c(3L, 5L),
    c(4L),
    c(2L, 7L),
    c(6L, 8L),
    c(7L)
  )
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  local.max.by.scale <- matrix(FALSE, nrow = length(adj.list), ncol = 2L)
  score.by.scale <- matrix(0, nrow = length(adj.list), ncol = 2L)
  colnames(local.max.by.scale) <- c("s1", "s2")
  colnames(score.by.scale) <- c("s1", "s2")

  local.max.by.scale[c(1L, 2L, 5L, 8L), 1L] <- TRUE
  score.by.scale[c(1L, 2L, 5L, 8L), 1L] <- c(0.7, 0.9, 0.8, 0.95)
  local.max.by.scale[c(1L, 5L, 6L), 2L] <- TRUE
  score.by.scale[c(1L, 5L, 6L), 2L] <- c(0.6, 0.75, 0.85)

  ref <- gflow:::.suppress.graph.endpoint.maxima.by.scale(
    adj.list = adj.list,
    weight.list = weight.list,
    local.max.by.scale = local.max.by.scale,
    score.by.scale = score.by.scale,
    radius = 1,
    prefer.cpp = FALSE
  )
  cpp <- gflow:::.suppress.graph.endpoint.maxima.by.scale(
    adj.list = adj.list,
    weight.list = weight.list,
    local.max.by.scale = local.max.by.scale,
    score.by.scale = score.by.scale,
    radius = 1,
    prefer.cpp = TRUE
  )

  expect_equal(cpp, ref, ignore_attr = TRUE)
  expect_setequal(which(cpp[, 1L]), c(2L, 5L, 8L))
  expect_setequal(which(cpp[, 2L]), c(1L, 5L, 6L))
})

test_that("compute.graph.endpoint.scores supports robust quantile scoring", {
  n <- 31L
  adj.list <- lapply(seq_len(n), function(i) integer(0))
  for (i in 1:29) {
    adj.list[[i]] <- c(adj.list[[i]], i + 1L)
    adj.list[[i + 1L]] <- c(adj.list[[i + 1L]], i)
  }
  adj.list[[2]] <- c(adj.list[[2]], 31L)
  adj.list[[31]] <- c(2L)
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  layout.3d <- matrix(0, nrow = n, ncol = 3)
  layout.3d[1:30, ] <- cbind(seq(0, 29), 0, 0)
  layout.3d[31, ] <- c(1, 8, 0)

  res <- compute.graph.endpoint.scores(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    k = 21L,
    q = 0.1,
    min.neighborhood.size = 3L
  )

  expect_true(res$s.q[1] > res$s.min[1])
  expect_true(res$score[1] > 0.9)
})

test_that("detect.graph.endpoints finds major arm tips on a Y graph", {
  adj.list <- list(
    c(2L),
    c(1L, 3L),
    c(2L, 4L),
    c(3L, 5L),
    c(4L, 6L, 9L),
    c(5L, 7L),
    c(6L, 8L),
    c(7L),
    c(5L, 10L),
    c(9L, 11L),
    c(10L)
  )
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  layout.3d <- matrix(0, nrow = 11, ncol = 3)
  layout.3d[1:5, ] <- cbind(seq(-4, 0, length.out = 5), 0, 0)
  layout.3d[6:8, ] <- cbind(c(1, 2, 3), c(1, 2, 3), 0)
  layout.3d[9:11, ] <- cbind(c(1, 2, 3), c(-1, -2, -3), 0)

  res <- detect.graph.endpoints(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    k = c(4L, 6L, 8L),
    min.neighborhood.size = 3L,
    detect.min.neighborhood.size = 2L
  )

  expect_setequal(res$endpoints, c(1L, 8L, 11L))
  expect_true(all(res$summary$scale.stability[res$endpoints] >= 1))
})

test_that("detect.graph.endpoints supports radius neighborhoods and optional smoothing", {
  n <- 21L
  adj.list <- lapply(seq_len(n), function(i) {
    nbrs <- integer(0)
    if (i > 1L) nbrs <- c(nbrs, i - 1L)
    if (i < n) nbrs <- c(nbrs, i + 1L)
    nbrs
  })
  weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
  layout.3d <- cbind(seq_len(n), 0, 0)

  radius.res <- detect.graph.endpoints(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    neighborhood = "geodesic_radius",
    radius = c(4, 6, 8),
    min.neighborhood.size = 3L,
    detect.min.neighborhood.size = 2L
  )

  expect_setequal(radius.res$endpoints, c(1L, n))

  smooth.res <- detect.graph.endpoints(
    adj.list = adj.list,
    weight.list = weight.list,
    layout.3d = layout.3d,
    k = c(5L, 7L, 9L),
    min.neighborhood.size = 3L,
    detect.min.neighborhood.size = 2L,
    smooth = TRUE,
    smooth.fit.args = list(k = 2L, verbose.level = 0L),
    smooth.refit.args = list(per.column.gcv = FALSE, verbose = FALSE)
  )

  expect_equal(sort(smooth.res$endpoints), c(2L, n - 1L))
  expect_true(!is.null(smooth.res$smoothed.scores))
  expect_true(all(smooth.res$scale.stability[smooth.res$endpoints] >= 1))
})

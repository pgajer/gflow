test_that("compute.tube.lens.corridor builds base and excess variants from shared searches", {
  adj.list <- list(
    c(2, 7),    # 1
    c(1, 3, 5, 6), # 2
    c(2, 4, 5),    # 3
    c(3),          # 4
    c(2, 3),       # 5
    c(2),          # 6
    c(1)           # 7
  )
  weight.list <- lapply(adj.list, function(x) rep(1, length(x)))

  res.base <- compute.tube.lens.corridor(
    adj.list = adj.list,
    weight.list = weight.list,
    start.vertex = 1,
    end.vertex = 4,
    path.relative.radius = 0.4,
    mode = "base"
  )

  expect_equal(res.base$path.vertices, c(1L, 2L, 3L, 4L))
  expect_equal(res.base$path.arc.length, c(0, 1/3, 2/3, 1))
  expect_equal(res.base$path.length, 3)
  expect_equal(res.base$tube.radius, 1.2)
  expect_equal(res.base$corridor.vertices, c(1L, 2L, 3L, 4L, 5L, 6L))
  expect_equal(res.base$t.balance, c(0, 1/3, 2/3, 1, 1/2, 1/3))
  expect_equal(res.base$harmonic.t, c(0, 1/3, 2/3, 1, 1/2, 1/3))
  expect_equal(res.base$distance.to.path, c(0, 0, 0, 0, 1, 1))
  expect_equal(res.base$excess, c(0, 0, 0, 0, 1, 2))
  expect_false(7L %in% res.base$corridor.vertices)
  expect_true(is.null(res.base$excess.vertices))
  expect_equal(res.base$selected.vertices, res.base$corridor.vertices)

  res.both <- compute.tube.lens.corridor(
    adj.list = adj.list,
    weight.list = weight.list,
    start.vertex = 1,
    end.vertex = 4,
    path.relative.radius = 0.4,
    excess.tol = 0.5,
    mode = "both"
  )

  expect_equal(res.both$corridor.vertices, c(1L, 2L, 3L, 4L, 5L, 6L))
  expect_equal(res.both$excess.vertices, c(1L, 2L, 3L, 4L))
  expect_equal(res.both$selected.vertices, res.both$corridor.vertices)
  expect_equal(res.both$excess.tolerance, 0.5)

  res.excess.default <- compute.tube.lens.corridor(
    adj.list = adj.list,
    weight.list = weight.list,
    start.vertex = 1,
    end.vertex = 4,
    path.relative.radius = 0.4,
    mode = "excess"
  )

  expect_equal(res.excess.default$excess.tolerance, res.excess.default$tube.radius)
  expect_equal(res.excess.default$selected.vertices, res.excess.default$excess.vertices)
})

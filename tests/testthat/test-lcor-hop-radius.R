test_that("lcor supports hop.radius > 1 for unit and derivative weighting", {
  adj.list <- list(
    c(2L),      # 1 -- 2
    c(1L, 3L),  # 2 -- 1,3
    c(2L)       # 3 -- 2
  )
  weight.list <- list(
    c(1.0),     # length 1-2
    c(1.0, 2.0),# lengths 2-1 and 2-3
    c(2.0)      # length 3-2
  )

  y <- c(1.0, 2.0, 4.0)
  z <- c(2.0, 3.0, 1.0)

  # hop.radius = 2 expands neighborhood to include all vertices
  unit.res <- lcor(
    adj.list, weight.list,
    y, z,
    type = "unit",
    hop.radius = 2
  )
  deriv.res <- lcor(
    adj.list, weight.list,
    y, z,
    type = "derivative",
    hop.radius = 2
  )

  expect_equal(
    as.numeric(unit.res),
    c(-0.4472136, -0.6, -0.8682431),
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(deriv.res),
    c(0.4472136, 0.0, -0.8944272),
    tolerance = 1e-6
  )
})

test_that("graph.cltr.imputation imputes using k nearest labeled neighbors", {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    skip("igraph not installed")
  }

  adj.list <- list(
    c(2L),        # 1
    c(1L, 3L),    # 2
    c(2L, 4L),    # 3
    c(3L),        # 4
    c(6L),        # 5 (disconnected component)
    c(5L)         # 6
  )
  weight.list <- list(
    c(1),
    c(1, 1),
    c(1, 1),
    c(1),
    c(1),
    c(1)
  )

  y <- c("A", NA, "B", NA, NA, NA)
  names(y) <- paste0("v", seq_along(y))

  res <- graph.cltr.imputation(adj.list, weight.list, y, k = 2)

  expect_equal(unname(res[1]), "A")
  expect_equal(unname(res[2]), "A")  # tie broken by nearest-neighbor order
  expect_equal(unname(res[3]), "B")
  expect_equal(unname(res[4]), "B")
  expect_true(is.na(res[5]))  # no labeled vertices reachable
  expect_true(is.na(res[6]))
  expect_equal(names(res), names(y))
})

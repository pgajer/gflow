.radius_edge_keys <- function(edge.matrix) {
  if (is.null(edge.matrix) || nrow(edge.matrix) == 0L) {
    return(character(0))
  }
  sort(paste(edge.matrix[, 1L], edge.matrix[, 2L], sep = "-"))
}


test_that("create.radius.graph matches the fixed-radius definition", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  g <- create.radius.graph(X, radius = 1.1)

  expect_s3_class(g, "radius_graph")
  expect_equal(g$graph_rule, "fixed.radius")
  expect_equal(g$radius, 1.1)
  expect_equal(.radius_edge_keys(g$edge_matrix), "1-2")
  expect_equal(g$edge_weight, 1, tolerance = 1e-12)
  expect_equal(g$n_components_before, 2L)
  expect_equal(g$n_components_after, 2L)
})


test_that("create.radius.graph can repair components by component MST", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  g <- create.radius.graph(X, radius = 1.1, connect.components = TRUE)

  expect_equal(g$n_components_before, 2L)
  expect_equal(g$n_components_after, 1L)
  expect_equal(g$n_mst_edges_added, 1L)
  expect_equal(.radius_edge_keys(g$mst_edge_matrix), "2-3")
  expect_equal(.radius_edge_keys(g$edge_matrix), c("1-2", "2-3"))
})


test_that("create.adaptive.radius.graph distinguishes max and min rules", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  max.g <- create.adaptive.radius.graph(X, k.scale = 1, radius.rule = "max")
  min.g <- create.adaptive.radius.graph(X, k.scale = 1, radius.rule = "min")

  expect_s3_class(max.g, "adaptive_radius_graph")
  expect_equal(max.g$sigma, c(1, 1, 2), tolerance = 1e-12)
  expect_equal(min.g$sigma, c(1, 1, 2), tolerance = 1e-12)
  expect_equal(.radius_edge_keys(max.g$edge_matrix), c("1-2", "2-3"))
  expect_equal(.radius_edge_keys(min.g$edge_matrix), "1-2")
  expect_equal(max.g$radius_rule, "max")
  expect_equal(min.g$radius_rule, "min")
})


test_that("radius graph constructors validate inputs", {
  X <- matrix(c(0, 1, 3), ncol = 1)

  expect_error(create.radius.graph(X, radius = 0), "radius")
  expect_error(create.radius.graph(X, radius = Inf), "radius")
  expect_error(create.radius.graph(c("a", "b"), radius = 1), "matrix or data frame")
  expect_error(create.adaptive.radius.graph(X, k.scale = 0), "k.scale")
  expect_error(create.adaptive.radius.graph(X, k.scale = 3), "k.scale")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1, radius.factor = 0),
               "radius.factor")
  expect_error(create.adaptive.radius.graph(X, k.scale = 1, radius.rule = "mean"),
               "'arg' should be one of")
})

test_that("quadform embedding, gradient, and metric match definitions", {
  X <- rbind(c(1, 2), c(2, 1))

  expect_equal(quadform.embed(X, index.k = 1)[, "q"], c(-3, 3))
  expect_equal(quadform.embed(X, index.k = 2)[, "q"], c(5, 5))

  grad <- quadform.gradient(matrix(c(1, 2), nrow = 1), index.k = 1)
  expect_equal(grad, matrix(c(2, -4), nrow = 1))

  metric <- quadform.metric(matrix(c(1, 2), nrow = 1), index.k = 1)
  expect_equal(metric, diag(2) + tcrossprod(c(2, -4)))
})


test_that("quadform.edge.length matches numerical quadrature", {
  u <- c(0.2, -0.4)
  v <- c(0.9, 0.3)
  h <- v - u
  signs <- c(1, -1)
  integrand <- function(t) {
    vapply(t, function(tt) {
      x <- u + tt * h
      dq <- 2 * sum(x * signs * h)
      sqrt(sum(h^2) + dq^2)
    }, numeric(1))
  }

  numeric.length <- stats::integrate(integrand, lower = 0, upper = 1,
                                     rel.tol = 1e-12)$value

  expect_equal(quadform.edge.length(u, v, index.k = 1), numeric.length,
               tolerance = 1e-10)
  expect_equal(quadform.edge.length(u, u, index.k = 1), 0)
})


test_that("quadform.reference.geodesics returns finite symmetric sample distances", {
  X <- rbind(
    c(0, 0),
    c(0.5, 0),
    c(0, 0.5)
  )

  ref <- quadform.reference.geodesics(
    X,
    index.k = 2,
    domain.radius = 1,
    grid.size = 9,
    sample.connection.k = 4
  )

  expect_equal(dim(ref$distances), c(3L, 3L))
  expect_true(all(is.finite(ref$distances)))
  expect_equal(ref$distances, t(ref$distances), tolerance = 1e-12)
  expect_equal(diag(ref$distances), rep(0, 3), tolerance = 1e-12)
  expect_equal(ref$n_sample_vertices, 3L)
  expect_true(ref$n_reference_vertices > ref$n_sample_vertices)
  expect_true(ref$n_edges > 0)
})


test_that("quadform.sample.dataset samples parameter disk data with reference distances", {
  ds <- quadform.sample.dataset(
    n = 12,
    index.k = 1,
    domain.radius = 1.25,
    grid.size = 11,
    sample.connection.k = 4,
    seed = 10
  )

  expect_s3_class(ds, "quadform_sample_dataset")
  expect_equal(dim(ds$X_param), c(12L, 2L))
  expect_equal(dim(ds$X_embed), c(12L, 3L))
  expect_equal(ds$X_embed, quadform.embed(ds$X_param, index.k = 1))
  expect_equal(ds$q, as.numeric(ds$X_embed[, "q"]))
  expect_true(all(sqrt(rowSums(ds$X_param^2)) <= 1.25 * (1 + 1e-12)))
  expect_equal(dim(ds$D_geodesic), c(12L, 12L))
  expect_equal(ds$distances, ds$D_geodesic)
  expect_true(all(is.finite(ds$D_geodesic)))
  expect_equal(ds$D_geodesic, t(ds$D_geodesic), tolerance = 1e-12)
  expect_equal(diag(ds$D_geodesic), rep(0, 12), tolerance = 1e-12)

  expect_equal(ds$metadata$sample_method, "uniform.parameter.disk")
  expect_equal(ds$metadata$index_k, 1L)
  expect_equal(ds$metadata$dim, 2L)
  expect_equal(ncol(ds$reference$grid_param), 2L)
  expect_equal(ncol(ds$reference$grid_embed), 3L)
  expect_equal(ds$reference$vertices_embed,
               quadform.embed(ds$reference$vertices_param, index.k = 1))
  expect_equal(ds$reference$sample_vertex, seq_len(12L))
  expect_true(nrow(ds$reference$edge_matrix) > 0)
})


test_that("quadform.sample.dataset seed is reproducible and local", {
  set.seed(99)
  before <- runif(3)
  state <- .Random.seed
  first <- quadform.sample.dataset(n = 5, index.k = 2, grid.size = 7, seed = 123)
  expect_equal(.Random.seed, state)
  second <- quadform.sample.dataset(n = 5, index.k = 2, grid.size = 7, seed = 123)
  expect_equal(first$X_param, second$X_param)
  expect_equal(first$D_geodesic, second$D_geodesic)

  set.seed(99)
  expect_equal(before, runif(3))
})


test_that("quadform utilities validate inputs", {
  expect_error(quadform.embed(matrix(1:4, ncol = 2), index.k = 3), "index.k")
  expect_error(quadform.edge.length(c(1, 2), c(1, 2, 3), index.k = 1),
               "equal length")
  expect_error(quadform.reference.geodesics(matrix(1:9, ncol = 3), index.k = 1),
               "2D")
  expect_error(
    quadform.reference.geodesics(matrix(c(2, 0, 0, 0), ncol = 2),
                                 index.k = 1, domain.radius = 1),
    "inside"
  )
  expect_error(quadform.sample.dataset(n = 0, index.k = 1), "positive integer")
  expect_error(quadform.sample.dataset(n = 3, index.k = 3), "index.k")
  expect_error(quadform.sample.dataset(n = 3, index.k = 1, domain.radius = 0),
               "domain.radius")
  expect_error(quadform.sample.dataset(n = 3, index.k = 1, seed = Inf), "seed")
})

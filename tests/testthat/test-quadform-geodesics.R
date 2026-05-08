test_that("quadform embedding, gradient, and metric match definitions", {
  X <- rbind(c(1, 2), c(2, 1))

  expect_equal(quadform.embed(X, index.k = 1)[, "q"], c(-3, 3))
  expect_equal(quadform.embed(X, index.k = 2)[, "q"], c(5, 5))
  expect_equal(quadform.embed(X, index.k = 2, coefficients = c(2, 4))[, "q"],
               c(18, 12))

  grad <- quadform.gradient(matrix(c(1, 2), nrow = 1), index.k = 1)
  expect_equal(grad, matrix(c(2, -4), nrow = 1))
  grad.coeff <- quadform.gradient(matrix(c(1, 2), nrow = 1), index.k = 2,
                                  coefficients = c(2, 4))
  expect_equal(grad.coeff, matrix(c(4, 16), nrow = 1))

  metric <- quadform.metric(matrix(c(1, 2), nrow = 1), index.k = 1)
  expect_equal(metric, diag(2) + tcrossprod(c(2, -4)))
  metric.coeff <- quadform.metric(matrix(c(1, 2), nrow = 1), index.k = 2,
                                  coefficients = c(2, 4))
  expect_equal(metric.coeff, diag(2) + tcrossprod(c(4, 16)))
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

  coefficients <- c(2, 4)
  coeff.integrand <- function(t) {
    vapply(t, function(tt) {
      x <- u + tt * h
      dq <- 2 * sum(x * signs * coefficients * h)
      sqrt(sum(h^2) + dq^2)
    }, numeric(1))
  }
  coeff.numeric.length <- stats::integrate(coeff.integrand, lower = 0,
                                           upper = 1, rel.tol = 1e-12)$value
  expect_equal(quadform.edge.length(u, v, index.k = 1,
                                    coefficients = coefficients),
               coeff.numeric.length, tolerance = 1e-10)
})


test_that("quadform.edge.lengths matches scalar edge lengths in 2D and 3D", {
  U <- rbind(
    c(0.2, -0.4),
    c(-0.5, 0.3),
    c(0, 0)
  )
  V <- rbind(
    c(0.9, 0.3),
    c(0.1, -0.2),
    c(0, 0)
  )
  expected <- vapply(seq_len(nrow(U)), function(i) {
    quadform.edge.length(U[i, ], V[i, ], index.k = 1,
                         coefficients = c(2, 4))
  }, numeric(1))
  expect_equal(quadform.edge.lengths(U, V, index.k = 1,
                                     coefficients = c(2, 4)),
               expected, tolerance = 1e-12)

  U3 <- rbind(
    c(0.1, -0.2, 0.3),
    c(-0.4, 0.2, -0.1)
  )
  V3 <- rbind(
    c(0.6, 0.1, -0.2),
    c(0.2, -0.3, 0.4)
  )
  expected3 <- vapply(seq_len(nrow(U3)), function(i) {
    quadform.edge.length(U3[i, ], V3[i, ], index.k = 2,
                         coefficients = c(1, 2, 3))
  }, numeric(1))
  expect_equal(quadform.edge.lengths(U3, V3, index.k = 2,
                                     coefficients = c(1, 2, 3)),
               expected3, tolerance = 1e-12)
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


test_that("quadform.grid.geodesic.distances matches R reference on small grids", {
  X <- rbind(
    c(0, 0),
    c(0.5, 0),
    c(0, 0.5)
  )

  r.ref <- quadform.reference.geodesics(
    X,
    index.k = 2,
    domain.radius = 1,
    grid.size = 11,
    sample.connection.k = 4
  )
  cpp.ref <- quadform.grid.geodesic.distances(
    X,
    index.k = 2,
    domain.radius = 1,
    grid.size = 11,
    sample.connection.k = 4
  )

  expect_equal(cpp.ref$distances, r.ref$distances, tolerance = 1e-12)
  expect_equal(cpp.ref$n_sample_vertices, 3L)
  expect_equal(cpp.ref$grid_size, 11L)
  expect_equal(cpp.ref$sample_connection_k, 4L)
  expect_equal(cpp.ref$index_k, 2L)

  r.coeff <- quadform.reference.geodesics(
    X,
    index.k = 2,
    coefficients = c(2, 4),
    domain.radius = 1,
    grid.size = 11,
    sample.connection.k = 4
  )
  cpp.coeff <- quadform.grid.geodesic.distances(
    X,
    index.k = 2,
    coefficients = c(2, 4),
    domain.radius = 1,
    grid.size = 11,
    sample.connection.k = 4,
    oracle = "none"
  )

  expect_equal(cpp.coeff$distances, r.coeff$distances, tolerance = 1e-12)
  expect_equal(cpp.coeff$coefficients, c(2, 4))
})


test_that("quadform grid geodesics support square domains", {
  X <- rbind(
    c(-1, -1),
    c(1, 1),
    c(0, 0),
    c(0.5, -0.75)
  )

  r.ref <- quadform.reference.geodesics(
    X,
    index.k = 2,
    coefficients = c(2, 4),
    domain.radius = 1,
    domain.shape = "square",
    grid.size = 11,
    sample.connection.k = 4
  )
  cpp.ref <- quadform.grid.geodesic.distances(
    X,
    index.k = 2,
    coefficients = c(2, 4),
    domain.radius = 1,
    domain.shape = "square",
    grid.size = 11,
    sample.connection.k = 4,
    oracle = "none"
  )

  expect_equal(cpp.ref$distances, r.ref$distances, tolerance = 1e-12)
  expect_equal(cpp.ref$domain_shape, "square")
  expect_equal(cpp.ref$n_reference_vertices, 11L * 11L)
  expect_true(all(is.finite(cpp.ref$distances)))
  expect_equal(cpp.ref$distances, t(cpp.ref$distances), tolerance = 1e-12)
})


test_that("quadform.grid.geodesic.distances can return sample-path oracle distances", {
  X <- cbind(x1 = c(0, 0.5, 1), x2 = c(0, 0, 0))

  ref <- quadform.grid.geodesic.distances(
    X,
    index.k = 2,
    domain.radius = 1,
    grid.size = 21,
    sample.connection.k = 4,
    oracle = "sample.path",
    oracle.tube.radius = 1e-6,
    return.oracle.paths = TRUE
  )

  X.embed <- quadform.embed(X, index.k = 2)
  chord.length <- function(i, j) {
    sqrt(sum((X.embed[i, ] - X.embed[j, ])^2))
  }
  expected.13 <- chord.length(1, 2) + chord.length(2, 3)

  expect_equal(ref$oracle_distances[1, 3], expected.13, tolerance = 1e-12)
  expect_equal(ref$oracle_distances, t(ref$oracle_distances), tolerance = 1e-12)
  expect_equal(diag(ref$oracle_distances), rep(0, 3), tolerance = 1e-12)
  expect_equal(ref$oracle_n_points[1, 3], 3L)
  expect_equal(ref$oracle_status[1, 3], "ok")
  expect_equal(ref$oracle_paths[[3]], 1:3)
  expect_equal(ref$oracle_method, "sample.path")
  expect_equal(ref$oracle_tube_radius, 1e-6)
})


test_that("quadform.delaunay.geodesic.distances returns finite 3D sample distances", {
  skip_if_not_installed("geometry")

  X <- rbind(
    c(0, 0, 0),
    c(0.5, 0, 0),
    c(0, 0.5, 0),
    c(0, 0, 0.5),
    c(-0.35, -0.2, 0.25)
  )

  ref <- quadform.delaunay.geodesic.distances(
    X,
    index.k = 2,
    coefficients = c(1, 2, 3),
    domain.radius = 1,
    domain.shape = "cube",
    n.ref = 80,
    seed = 13,
    candidate.multiplier = 3,
    boundary.fraction = 0.25,
    edge.length.factor = 3
  )

  expect_s3_class(ref, "quadform_delaunay_geodesics")
  expect_equal(dim(ref$distances), c(nrow(X), nrow(X)))
  expect_true(all(is.finite(ref$distances)))
  expect_equal(ref$distances, t(ref$distances), tolerance = 1e-12)
  expect_equal(diag(ref$distances), rep(0, nrow(X)), tolerance = 1e-12)
  expect_equal(ref$sample_vertex, seq_len(nrow(X)))
  expect_true(ref$n_reference_vertices >= nrow(X))
  expect_true(ref$n_edges > 0)
  expect_equal(ref$n_components, 1L)
  expect_equal(ncol(ref$vertices_param), 3L)
  expect_equal(ncol(ref$vertices_embed), 4L)
  expect_equal(ref$n_edges_unfiltered, ref$n_delaunay_edges)
  expect_true(is.numeric(ref$retained_edge_fraction))
  expect_true(ref$retained_edge_fraction > 0 && ref$retained_edge_fraction <= 1)
  expect_true(is.data.frame(ref$filter_attempts))
  expect_true(all(c("filter_factor", "n_edges", "n_components", "connected") %in%
                    names(ref$filter_attempts)))
  expect_equal(tail(ref$filter_attempts$n_components, 1L), ref$n_components)
})


test_that("C++ Delaunay edge extractor matches R/Qhull edge sets", {
  skip_if_not_installed("geometry")

  set.seed(42)
  X <- matrix(rnorm(75), ncol = 3)
  X[, 3] <- X[, 3] + seq_len(nrow(X)) * 1e-4

  r.edges <- gflow:::.quadform.delaunay.edges.3d(X, backend = "geometry")
  cpp <- gflow:::rcpp_quadform_delaunay_edges_3d(X, "Qt Qbb Qc")

  edge.keys <- function(edge.matrix) {
    apply(edge.matrix, 1L, function(edge) paste(sort(edge), collapse = "-"))
  }

  expect_true(is.matrix(cpp$edge_matrix))
  expect_true(is.matrix(cpp$tetrahedra))
  expect_equal(ncol(cpp$edge_matrix), 2L)
  expect_equal(ncol(cpp$tetrahedra), 4L)
  expect_equal(cpp$n_edges, nrow(cpp$edge_matrix))
  expect_equal(cpp$n_tetrahedra, nrow(cpp$tetrahedra))
  expect_equal(cpp$exit_code, 0L)
  expect_equal(sort(edge.keys(cpp$edge_matrix)),
               sort(edge.keys(r.edges)))
})


test_that("quadform.delaunay.geodesic.distances C++ backend matches geometry backend", {
  skip_if_not_installed("geometry")

  X <- rbind(
    c(0, 0, 0),
    c(0.5, 0, 0),
    c(0, 0.5, 0),
    c(0, 0, 0.5),
    c(-0.35, -0.2, 0.25)
  )
  args <- list(
    X = X,
    index.k = 2,
    coefficients = c(1, 2, 3),
    domain.radius = 1,
    domain.shape = "cube",
    n.ref = 90,
    seed = 13,
    candidate.multiplier = 3,
    boundary.fraction = 0.25,
    edge.length.factor = 4
  )

  cpp <- do.call(quadform.delaunay.geodesic.distances,
                 c(args, list(delaunay.backend = "cpp")))
  geom <- do.call(quadform.delaunay.geodesic.distances,
                  c(args, list(delaunay.backend = "geometry")))

  edge.keys <- function(edge.matrix) {
    apply(edge.matrix, 1L, function(edge) paste(sort(edge), collapse = "-"))
  }

  expect_equal(cpp$delaunay_backend, "cpp")
  expect_equal(geom$delaunay_backend, "geometry")
  expect_equal(cpp$vertices_param, geom$vertices_param, tolerance = 1e-14)
  expect_equal(cpp$epsilon, geom$epsilon, tolerance = 1e-14)
  expect_equal(cpp$n_delaunay_edges, geom$n_delaunay_edges)
  expect_equal(cpp$n_edges, geom$n_edges)
  expect_equal(sort(edge.keys(cpp$edge_matrix)),
               sort(edge.keys(geom$edge_matrix)))
  expect_equal(cpp$filter_attempts$n_edges, geom$filter_attempts$n_edges)
  expect_equal(cpp$filter_attempts$n_components,
               geom$filter_attempts$n_components)
  expect_equal(cpp$filter_attempts$connected, geom$filter_attempts$connected)
  expect_equal(cpp$filter_factor_used, geom$filter_factor_used)
  expect_equal(cpp$n_components, geom$n_components)
  expect_equal(cpp$distances, geom$distances, tolerance = 1e-10)
})


test_that("quadform.delaunay.geodesic.distances supports disabling edge filtering", {
  skip_if_not_installed("geometry")

  X <- rbind(
    c(0, 0, 0),
    c(0.5, 0, 0),
    c(0, 0.5, 0),
    c(0, 0, 0.5),
    c(-0.35, -0.2, 0.25)
  )

  ref <- quadform.delaunay.geodesic.distances(
    X,
    index.k = 2,
    coefficients = c(1, 2, 3),
    domain.radius = 1,
    domain.shape = "cube",
    n.ref = 80,
    seed = 13,
    candidate.multiplier = 3,
    boundary.fraction = 0.25,
    edge.length.factor = Inf
  )

  expect_s3_class(ref, "quadform_delaunay_geodesics")
  expect_true(all(is.finite(ref$distances)))
  expect_true(is.infinite(ref$filter_factor_requested))
  expect_true(is.infinite(ref$filter_factor_used))
  expect_equal(ref$n_edges, ref$n_delaunay_edges)
  expect_equal(ref$n_edges_unfiltered, ref$n_delaunay_edges)
  expect_equal(ref$retained_edge_fraction, 1)
  expect_false(ref$filter_relaxation_happened)
  expect_false(ref$relaxed_to_inf)
  expect_equal(nrow(ref$filter_attempts), 1L)
  expect_true(is.infinite(ref$filter_attempts$filter_factor))
  expect_equal(ref$filter_attempts$retained_edge_fraction, 1)
})


test_that("quadform.delaunay.geodesic.distances defaults to factor 4 diagnostics", {
  skip_if_not_installed("geometry")

  X <- rbind(
    c(0, 0, 0),
    c(0.5, 0, 0),
    c(0, 0.5, 0),
    c(0, 0, 0.5),
    c(-0.35, -0.2, 0.25)
  )

  ref <- quadform.delaunay.geodesic.distances(
    X,
    index.k = 2,
    coefficients = c(1, 2, 3),
    domain.radius = 1,
    domain.shape = "cube",
    n.ref = 80,
    seed = 13,
    candidate.multiplier = 3,
    boundary.fraction = 0.25
  )

  expect_equal(ref$filter_factor_requested, 4)
  expect_true(ref$filter_factor_used >= 4 || is.infinite(ref$filter_factor_used))
  expect_true(all(ref$filter_attempts$n_edges <= ref$n_delaunay_edges))
  expect_equal(tail(ref$filter_attempts$n_edges, 1L), ref$n_edges)
})


test_that("quadform.grid.geodesic.calibration summarizes grid convergence", {
  cal <- quadform.grid.geodesic.calibration(
    index.k = 1,
    grid.sizes = c(11, 21),
    reference.grid.size = 31,
    n.sources = 4,
    targets.per.source = 3,
    seed = 11
  )

  expect_s3_class(cal, "quadform_grid_geodesic_calibration")
  expect_equal(cal$metadata$n_pairs, 12L)
  expect_equal(sort(cal$summary$grid_size), c(11L, 21L, 31L))
  expect_true(all(c("frob_rel_error", "q95_rel_error", "runtime_sec") %in%
                    names(cal$summary)))
  ref.row <- cal$summary[cal$summary$grid_size == 31L, ]
  expect_equal(ref.row$frob_rel_error, 0, tolerance = 1e-12)
  expect_equal(ref.row$max_rel_error, 0, tolerance = 1e-12)
  expect_equal(nrow(cal$pair_points), 12L)
  expect_equal(ncol(cal$pair_points), 4L)

  cal2 <- quadform.grid.geodesic.calibration(
    index.k = 1,
    grid.sizes = c(11, 21),
    reference.grid.size = 31,
    n.sources = 4,
    targets.per.source = 3,
    seed = 11
  )
  expect_equal(cal$pair_points, cal2$pair_points)
  expect_equal(cal$summary$grid_size, cal2$summary$grid_size)
  expect_equal(cal$summary$frob_rel_error, cal2$summary$frob_rel_error)
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

  ds.coeff <- quadform.sample.dataset(
    n = 6,
    index.k = 2,
    coefficients = c(2, 4),
    grid.size = 9,
    seed = 2
  )
  expect_equal(ds.coeff$X_embed,
               quadform.embed(ds.coeff$X_param, index.k = 2,
                              coefficients = c(2, 4)))
  expect_equal(ds.coeff$metadata$coefficients, c(2, 4))

  ds.square <- quadform.sample.dataset(
    n = 10,
    index.k = 2,
    coefficients = c(2, 4),
    domain.radius = 1.5,
    sample.method = "uniform.parameter.square",
    grid.size = 9,
    seed = 3
  )
  expect_true(all(abs(ds.square$X_param) <= 1.5))
  expect_equal(ds.square$metadata$domain_shape, "square")
  expect_equal(ds.square$reference$domain_shape, "square")
  expect_equal(ds.square$X_embed,
               quadform.embed(ds.square$X_param, index.k = 2,
                              coefficients = c(2, 4)))
})


test_that("quadform.sample.dataset supports radial density variants", {
  methods <- c(
    "radial.center.parameter.disk",
    "radial.boundary.parameter.disk",
    "radial.center.parameter.square",
    "radial.boundary.parameter.square"
  )
  for (method in methods) {
    ds <- quadform.sample.dataset(
      n = 12,
      index.k = 1,
      sample.method = method,
      grid.size = 7,
      seed = 4
    )
    if (grepl("square$", method)) {
      expect_true(all(abs(ds$X_param) <= 1 + 1e-12))
      expect_equal(ds$metadata$domain_shape, "square")
    } else {
      expect_true(all(sqrt(rowSums(ds$X_param^2)) <= 1 + 1e-12))
      expect_equal(ds$metadata$domain_shape, "disk")
    }
    expect_true(all(is.finite(ds$D_geodesic)))
  }
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
  expect_error(quadform.grid.geodesic.distances(matrix(1:9, ncol = 3), index.k = 1),
               "2D")
  expect_error(
    quadform.grid.geodesic.distances(matrix(c(2, 0, 0, 0), ncol = 2),
                                     index.k = 1, domain.radius = 1),
    "inside"
  )
  expect_silent(
    quadform.grid.geodesic.distances(matrix(c(1, 1, 0, 0), ncol = 2),
                                     index.k = 1, domain.radius = 1,
                                     domain.shape = "square")
  )
  expect_error(
    quadform.grid.geodesic.distances(matrix(c(1.1, 0, 0, 0), ncol = 2),
                                     index.k = 1, domain.radius = 1,
                                     domain.shape = "square"),
    "inside"
  )
  expect_error(quadform.grid.geodesic.calibration(index.k = 3), "index.k")
  expect_error(quadform.grid.geodesic.calibration(index.k = 1, n.sources = 0),
               "n.sources")
  expect_error(quadform.grid.geodesic.distances(matrix(c(0, 0, 1, 0), ncol = 2),
                                                index.k = 1, oracle = "sample.path",
                                                oracle.tube.radius = -1),
               "oracle.tube.radius")
  expect_error(quadform.grid.geodesic.distances(matrix(c(0, 0, 1, 0), ncol = 2),
                                                index.k = 1, oracle.tube.k = 0),
               "oracle.tube.k")
  expect_error(quadform.sample.dataset(n = 0, index.k = 1), "positive integer")
  expect_error(quadform.sample.dataset(n = 3, index.k = 3), "index.k")
  expect_error(quadform.sample.dataset(n = 3, index.k = 1, domain.radius = 0),
               "domain.radius")
  expect_error(quadform.embed(matrix(1:4, ncol = 2), index.k = 1,
                              coefficients = c(1, 0)), "coefficients")
  expect_error(quadform.grid.geodesic.distances(matrix(c(0, 0, 1, 0), ncol = 2),
                                                index.k = 1,
                                                coefficients = c(1, 2, 3)),
               "coefficients")
  expect_error(quadform.sample.dataset(n = 3, index.k = 1, seed = Inf), "seed")
})

make_path_graph_weights <- function(weights) {
  n <- length(weights) + 1L
  adj <- vector("list", n)
  w <- vector("list", n)
  add_edge <- function(i, j, weight) {
    adj[[i]] <<- c(adj[[i]], j)
    w[[i]] <<- c(w[[i]], weight)
    adj[[j]] <<- c(adj[[j]], i)
    w[[j]] <<- c(w[[j]], weight)
  }
  for (i in seq_along(weights)) add_edge(i, i + 1L, weights[i])
  list(adj.list = lapply(adj, as.integer), weight.list = lapply(w, as.double))
}

make_three_arm_star_weights <- function(weights = rep(1, 6)) {
  adj <- vector("list", 7)
  w <- vector("list", 7)
  add_edge <- function(i, j, weight) {
    adj[[i]] <<- c(adj[[i]], j)
    w[[i]] <<- c(w[[i]], weight)
    adj[[j]] <<- c(adj[[j]], i)
    w[[j]] <<- c(w[[j]], weight)
  }
  add_edge(1, 2, weights[1])
  add_edge(2, 3, weights[2])
  add_edge(1, 4, weights[3])
  add_edge(4, 5, weights[4])
  add_edge(1, 6, weights[5])
  add_edge(6, 7, weights[6])
  list(adj.list = lapply(adj, as.integer), weight.list = lapply(w, as.double))
}

make_grid_graph_weights <- function(nx, ny) {
  idx <- function(x, y) (y - 1L) * nx + x
  adj <- vector("list", nx * ny)
  w <- vector("list", nx * ny)
  add_edge <- function(i, j, weight = 1) {
    adj[[i]] <<- c(adj[[i]], j)
    w[[i]] <<- c(w[[i]], weight)
    adj[[j]] <<- c(adj[[j]], i)
    w[[j]] <<- c(w[[j]], weight)
  }
  for (x in seq_len(nx)) {
    for (y in seq_len(ny)) {
      if (x < nx) add_edge(idx(x, y), idx(x + 1L, y))
      if (y < ny) add_edge(idx(x, y), idx(x, y + 1L))
    }
  }
  coords <- do.call(rbind, lapply(seq_len(ny), function(y) {
    cbind(x = seq_len(nx) - 1, y = y - 1)
  }))
  list(
    adj.list = lapply(adj, as.integer),
    weight.list = lapply(w, as.double),
    coordinates = coords
  )
}

load_ssrhe_all_labeled_validation_definitions <- function() {
  runner <- "/Users/pgajer/current_projects/trend_filtering/development/ssrhe_hessian_energy/ssrhe_all_labeled_comparator_validation.R"
  skip_if_not(file.exists(runner),
              "SSRHE all-labeled validation runner is not available.")
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    skip("Package 'pkgload' is required for SSRHE validation regression cases.")
  }
  env <- new.env(parent = globalenv())
  exprs <- parse(runner, keep.source = FALSE)
  for (expr in exprs) {
    txt <- paste(deparse(expr), collapse = "\n")
    if (grepl("^manifest <- build\\.dataset\\.manifest", txt)) break
    if (grepl("pkgload::load_all\\(gflow.dir", txt, fixed = FALSE)) next
    eval(expr, envir = env)
  }
  env
}

fit_ssrhe_graph_trend_filtering_case <- function(dataset.id, order, variant) {
  env <- load_ssrhe_all_labeled_validation_definitions()
  manifest <- env$build.dataset.manifest()
  row <- manifest[manifest$dataset.id == dataset.id, , drop = FALSE]
  expect_equal(nrow(row), 1L)
  ds <- env$materialize.dataset(row)
  graph.info <- env$make.graph.payload(ds)
  weights <- if (identical(variant, "unit")) {
    graph.info$metric.weight.list
  } else {
    graph.info$conductance.weight.list
  }
  fit.graph.trend.filtering(
    adj.list = graph.info$metric.adj.list,
    weight.list = weights,
    y = ds$y,
    order = order,
    lambda.selection = "cv",
    weight.rule = variant,
    n.lambda = 20L,
    nfolds = 4L,
    maxsteps = 500L
  )
}

ssrhe_like_truth_function <- function(U) {
  y <- sin(pi * U[, 1])
  if (ncol(U) >= 2L) {
    y <- y + 0.7 * cos(pi * U[, 2]) +
      0.45 * U[, 1]^2 - 0.25 * U[, 1] * U[, 2]
  }
  as.numeric(scale(y, center = TRUE, scale = FALSE))
}

adj_weight_from_edge_matrix <- function(n, edge.matrix, edge.weight) {
  adj <- vector("list", n)
  weights <- vector("list", n)
  for (i in seq_len(n)) {
    adj[[i]] <- integer()
    weights[[i]] <- numeric()
  }
  for (r in seq_len(nrow(edge.matrix))) {
    i <- as.integer(edge.matrix[r, 1L])
    j <- as.integer(edge.matrix[r, 2L])
    w <- as.numeric(edge.weight[r])
    adj[[i]] <- c(adj[[i]], j)
    weights[[i]] <- c(weights[[i]], w)
    adj[[j]] <- c(adj[[j]], i)
    weights[[j]] <- c(weights[[j]], w)
  }
  list(adj.list = adj, weight.list = weights)
}

make_ssrhe_like_graph_trend_case <- function(kind = c("flat", "quadform")) {
  kind <- match.arg(kind)
  n <- 60L
  dim <- 2L
  if (identical(kind, "flat")) {
    seed <- 13200L
    set.seed(seed)
    U <- matrix(stats::runif(n * dim, -1, 1), ncol = dim)
    X <- U
    index.k <- NA_integer_
    coefficients <- rep(NA_real_, dim)
    order <- 0L
    variant <- "sqrt.conductance"
  } else {
    seed <- 16104L
    set.seed(seed)
    U <- matrix(stats::runif(n * dim, -1, 1), ncol = dim)
    index.k <- 1L
    coefficients <- 0.35 * c(1, 1)
    X <- quadform.embed(U, index.k = index.k, coefficients = coefficients)
    order <- 2L
    variant <- "conductance"
  }
  truth <- ssrhe_like_truth_function(U)
  set.seed(seed + 900000L)
  y <- truth + stats::rnorm(length(truth), sd = 0.08)
  graph <- create.adaptive.radius.graph(
    X,
    k.scale = 12L,
    radius.factor = 1.25,
    radius.rule = "geomean",
    prune.method = "none",
    connect.components = TRUE,
    connect.method = "component.mst"
  )
  edges <- graph$edge_matrix
  metric.length <- if (identical(kind, "flat")) {
    sqrt(rowSums((U[edges[, 1L], , drop = FALSE] -
                    U[edges[, 2L], , drop = FALSE])^2))
  } else {
    quadform.edge.lengths(
      U[edges[, 1L], , drop = FALSE],
      U[edges[, 2L], , drop = FALSE],
      index.k = index.k,
      coefficients = coefficients
    )
  }
  metric.graph <- adj_weight_from_edge_matrix(nrow(X), edges, metric.length)
  conductance.graph <- adj_weight_from_edge_matrix(
    nrow(X), edges, 1 / pmax(metric.length, 1e-8)
  )
  list(
    adj.list = metric.graph$adj.list,
    weight.list = if (identical(variant, "unit")) {
      metric.graph$weight.list
    } else {
      conductance.graph$weight.list
    },
    y = y,
    order = order,
    variant = variant,
    graph = graph
  )
}

test_that("graph trend-filtering operator validates and orients weighted graphs", {
  graph <- make_path_graph_weights(c(4, 9))
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    weight.rule = "sqrt.conductance"
  )

  expect_s3_class(op, "graph.trend.filtering.operator")
  expect_equal(op$edge.table$from, c(1L, 2L))
  expect_equal(op$edge.table$to, c(2L, 3L))
  expect_equal(op$edge.table$trend.weight, c(2, 3), tolerance = 1e-12)
  expect_equal(as.matrix(op$incidence$matrix), matrix(c(
    -2,  2,  0,
     0, -3,  3
  ), nrow = 2, byrow = TRUE), tolerance = 1e-12)
})

test_that("graph trend-filtering weight rules are explicit", {
  graph <- make_path_graph_weights(c(4, 9))

  op.conductance <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    weight.rule = "conductance"
  )
  op.sqrt <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    weight.rule = "sqrt.conductance"
  )
  op.unit <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    weight.rule = "unit"
  )

  expect_equal(op.conductance$edge.table$trend.weight, c(4, 9), tolerance = 1e-12)
  expect_equal(op.sqrt$edge.table$trend.weight, c(2, 3), tolerance = 1e-12)
  expect_equal(op.unit$edge.table$trend.weight, c(1, 1), tolerance = 1e-12)
})

test_that("graph trend filtering rejects unsupported orders above phase 2 scope", {
  graph <- make_path_graph_weights(c(1, 1))
  expect_error(
    graph.trend.filtering.operator(graph$adj.list, graph$weight.list, order = 3L),
    "supports only order = 0L, 1L, or 2L"
  )
})

test_that("higher-order unit path operators match dense matrix references", {
  graph <- make_path_graph_weights(rep(1, 3))
  op0 <- graph.trend.filtering.operator(
    graph$adj.list, graph$weight.list,
    order = 0L, weight.rule = "unit"
  )
  D.ref <- matrix(c(
    -1,  1,  0,  0,
     0, -1,  1,  0,
     0,  0, -1,  1
  ), nrow = 3, byrow = TRUE)
  L.ref <- t(D.ref) %*% D.ref

  expect_equal(as.matrix(op0$penalty$matrix), D.ref, tolerance = 1e-12)
  expect_equal(op0$penalty$kind, "incidence")
  expect_equal(op0$nullity.estimate, 1L)

  op1 <- graph.trend.filtering.operator(
    graph$adj.list, graph$weight.list,
    order = 1L, weight.rule = "unit"
  )
  expect_equal(as.matrix(op1$penalty$matrix), L.ref, tolerance = 1e-12)
  expect_equal(op1$penalty$kind, "laplacian")
  expect_equal(op1$nullity.estimate, 1L)

  op2 <- graph.trend.filtering.operator(
    graph$adj.list, graph$weight.list,
    order = 2L, weight.rule = "unit"
  )
  expect_equal(as.matrix(op2$penalty$matrix), D.ref %*% L.ref, tolerance = 1e-12)
  expect_equal(op2$penalty$kind, "incidence_laplacian")
  expect_equal(op2$nullity.estimate, 1L)
})

test_that("higher-order weighted path operators use selected weight convention", {
  graph <- make_path_graph_weights(c(4, 9, 16))
  op <- graph.trend.filtering.operator(
    graph$adj.list, graph$weight.list,
    order = 1L, weight.rule = "sqrt.conductance"
  )
  D.ref <- matrix(c(
    -2,  2,  0,  0,
     0, -3,  3,  0,
     0,  0, -4,  4
  ), nrow = 3, byrow = TRUE)
  expect_equal(as.matrix(op$incidence$matrix), D.ref, tolerance = 1e-12)
  expect_equal(as.matrix(op$laplacian$matrix), t(D.ref) %*% D.ref, tolerance = 1e-12)
  expect_equal(as.matrix(op$penalty$matrix), t(D.ref) %*% D.ref, tolerance = 1e-12)
})

test_that("path divided-difference path graph operators match genlasso references", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  y <- seq_len(6)
  graph <- make_path_graph_weights(rep(1, length(y) - 1L))

  for (ord in 0:2) {
    op <- graph.trend.filtering.operator(
      graph$adj.list,
      graph$weight.list,
      order = ord,
      operator.family = "path.divided.difference",
      path.weighting = "unit"
    )
    ref <- genlasso::trendfilter(y, ord = ord)
    expect_equal(as.matrix(op$penalty$matrix),
                 as.matrix(ref$pathobjs$D0),
                 tolerance = 1e-12)
    expect_equal(op$operator.family, "path.divided.difference")
    expect_equal(op$penalty$kind, "path_divided_difference")
    expect_equal(op$penalty$difference.order, ord + 1L)
  }
})

test_that("path divided-difference operators use metric path lengths", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  pos <- c(0, 1, 3, 6, 10, 15)
  graph <- make_path_graph_weights(diff(pos))
  y <- seq_along(pos)
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    order = 2L,
    operator.family = "path.divided.difference",
    path.weighting = "unit"
  )
  ref <- genlasso::trendfilter(y, pos = pos, ord = 2)

  expect_equal(as.matrix(op$penalty$matrix),
               as.matrix(ref$pathobjs$D0),
               tolerance = 1e-12)
  expect_equal(op$path.table$vertices, c("1-2-3-4", "2-3-4-5", "3-4-5-6"))
  expect_equal(op$path.table$metric.length, c(6, 9, 12), tolerance = 1e-12)
})

test_that("path divided-difference canonicalization removes reversed path rows", {
  graph <- make_path_graph_weights(rep(1, 3))
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    order = 1L,
    operator.family = "path.divided.difference",
    path.family = "all.simple",
    path.weighting = "unit"
  )

  expect_equal(nrow(op$penalty$matrix), 2L)
  expect_equal(op$path.table$vertices, c("1-2-3", "2-3-4"))
  expect_equal(op$path.table$duplicates.removed[1], 2L)
})

test_that("branch-continuation star graph operator preserves arm slopes", {
  graph <- make_three_arm_star_weights()
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    order = 1L,
    operator.family = "path.divided.difference",
    path.family = "branch.continuation",
    path.weighting = "unit"
  )

  expect_equal(op$path.table$vertices, c("1-2-3", "1-4-5", "1-6-7"))
  expect_equal(as.matrix(op$penalty$matrix), matrix(c(
     1, -2,  1,  0,  0,  0,  0,
     1,  0,  0, -2,  1,  0,  0,
     1,  0,  0,  0,  0, -2,  1
  ), nrow = 3, byrow = TRUE), tolerance = 1e-12)
  expect_equal(op$nullity.estimate, 4L)

  arm.linear <- c(10, 11, 12, 8, 6, 13, 16)
  expect_equal(as.vector(op$penalty$matrix %*% arm.linear),
               rep(0, 3), tolerance = 1e-12)
})

test_that("all-simple star graph operator includes cross-arm paths", {
  graph <- make_three_arm_star_weights()
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    order = 1L,
    operator.family = "path.divided.difference",
    path.family = "all.simple",
    path.weighting = "unit"
  )

  expect_true(all(c("1-2-3", "1-4-5", "1-6-7",
                    "2-1-4", "2-1-6", "4-1-6") %in%
                    op$path.table$vertices))
  expect_equal(op$nullity.estimate, 1L)
})

test_that("anchor-mean path weighting balances canonical start vertices", {
  graph <- make_three_arm_star_weights()
  op <- graph.trend.filtering.operator(
    graph$adj.list,
    graph$weight.list,
    order = 1L,
    operator.family = "path.divided.difference",
    path.family = "all.simple",
    path.weighting = "anchor.mean"
  )
  starts <- op$path.table$start
  expected <- as.double(1 / table(starts)[as.character(starts)])
  expect_equal(op$path.table$row.weight, expected, tolerance = 1e-12)
})

test_that("fixed lambda fit agrees with direct genlasso reference", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  graph <- make_path_graph_weights(c(1, 2, 1))
  y <- c(0, 0.2, 2.0, 2.1)
  lambda <- 0.35
  fit <- fit.graph.trend.filtering(
    graph$adj.list,
    graph$weight.list,
    y,
    lambda.grid = lambda,
    lambda.selection = "fixed",
    weight.rule = "conductance"
  )

  D <- fit$operator$incidence$matrix
  ref <- genlasso::genlasso(y = y, D = D)
  beta.ref <- as.vector(stats::coef(ref, lambda = lambda)$beta)
  expect_equal(fit$fitted.values, beta.ref, tolerance = 1e-9)
  expect_equal(fit$lambda, lambda, tolerance = 1e-12)
})

test_that("package-local SSRHE-like graph trend-filtering cases return finite fits", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  for (kind in c("flat", "quadform")) {
    case <- make_ssrhe_like_graph_trend_case(kind)
    set.seed(910000L + match(kind, c("flat", "quadform")))
    fit <- fit.graph.trend.filtering(
      adj.list = case$adj.list,
      weight.list = case$weight.list,
      y = case$y,
      order = case$order,
      lambda.selection = "cv",
      weight.rule = case$variant,
      n.lambda = 20L,
      nfolds = 4L,
      maxsteps = 500L
    )

    expect_s3_class(fit, "graph.trend.filtering.fit")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(is.finite(fit$lambda))
    expect_true(fit$path$svd)
    expect_equal(length(fit$fitted.values), length(case$y))
  }
})

test_that("SSRHE affected graph trend-filtering cases return finite fitted values", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  cases <- list(
    list(
      dataset.id = "flat_d2_rep1",
      order = 0L,
      variant = "sqrt.conductance"
    ),
    list(
      dataset.id = "quadform_d2_idx1_curv035_rep1",
      order = 2L,
      variant = "conductance"
    )
  )

  for (case in cases) {
    fit <- fit_ssrhe_graph_trend_filtering_case(
      dataset.id = case$dataset.id,
      order = case$order,
      variant = case$variant
    )
    expect_s3_class(fit, "graph.trend.filtering.fit")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(is.finite(fit$lambda))
    expect_equal(length(fit$fitted.values), fit$graph$n.vertices)
  }
})

test_that("path divided-difference fit matches genlasso trendfilter on a path graph", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  pos <- c(0, 0.5, 1.5, 2.5, 4, 7)
  graph <- make_path_graph_weights(diff(pos))
  y <- c(0, 0.2, 1.1, 0.9, 2.2, 2.8)
  lambda <- 0.08

  fit <- fit.graph.trend.filtering(
    graph$adj.list,
    graph$weight.list,
    y,
    order = 2L,
    operator.family = "path.divided.difference",
    path.weighting = "unit",
    lambda.grid = lambda,
    lambda.selection = "fixed"
  )
  ref <- genlasso::trendfilter(y, pos = pos, ord = 2)
  beta.ref <- as.vector(stats::coef(ref, lambda = lambda)$beta)

  expect_equal(fit$fitted.values, beta.ref, tolerance = 1e-8)
  expect_equal(fit$operator.family, "path.divided.difference")
  expect_equal(fit$path.weighting, "unit")
})

test_that("higher-order fixed lambda fits agree with direct genlasso references", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  graph <- make_path_graph_weights(rep(1, 5))
  y <- sin(seq_len(6))
  lambda <- 0.1

  for (ord in 1:2) {
    fit <- fit.graph.trend.filtering(
      graph$adj.list,
      graph$weight.list,
      y,
      order = ord,
      lambda.grid = lambda,
      lambda.selection = "fixed",
      weight.rule = "unit"
    )

    D <- fit$operator$penalty$matrix
    if (ord == 1L) {
      ref <- genlasso::genlasso(y = y, D = as.matrix(D), svd = TRUE)
    } else {
      ref <- genlasso::genlasso(y = y, D = D)
    }
    beta.ref <- as.vector(stats::coef(ref, lambda = lambda)$beta)
    expect_equal(fit$fitted.values, beta.ref, tolerance = 1e-9)
  }
})

test_that("unweighted path graph agrees with genlasso fusedlasso1d", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  y <- c(0, 0.1, 1.8, 2.0, 2.1)
  lambda <- 0.2
  graph <- make_path_graph_weights(rep(1, length(y) - 1L))

  fit <- fit.graph.trend.filtering(
    graph$adj.list,
    graph$weight.list,
    y,
    lambda.grid = lambda,
    lambda.selection = "fixed",
    weight.rule = "unit"
  )
  ref <- genlasso::fusedlasso1d(y)
  beta.ref <- as.vector(stats::coef(ref, lambda = lambda)$beta)

  expect_equal(fit$fitted.values, beta.ref, tolerance = 1e-8)
})

test_that("large lambda fuses connected graph to the response mean", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  graph <- make_path_graph_weights(rep(1, 4))
  y <- c(-1, 0, 1, 3, 7)
  fit <- fit.graph.trend.filtering(
    graph$adj.list,
    graph$weight.list,
    y,
    lambda.grid = 1e6,
    lambda.selection = "fixed",
    weight.rule = "unit"
  )

  expect_equal(fit$fitted.values, rep(mean(y), length(y)), tolerance = 1e-6)
})

test_that("small vertex CV returns selected lambda and diagnostics", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  graph <- make_path_graph_weights(rep(1, 5))
  y <- c(0, 0.1, 0.2, 2.0, 2.1, 2.2)
  lambda.grid <- c(1, 0.1, 0)
  fit <- fit.graph.trend.filtering(
    graph$adj.list,
    graph$weight.list,
    y,
    lambda.grid = lambda.grid,
    lambda.selection = "cv",
    weight.rule = "unit",
    foldid = c(1, 2, 3, 1, 2, 3),
    maxsteps = 50L
  )

  expect_s3_class(fit, "graph.trend.filtering.fit")
  expect_true(fit$lambda %in% lambda.grid)
  expect_equal(dim(fit$cv$fold.errors), c(3, 3))
  expect_equal(length(fit$fitted.values), length(y))
})

test_that("transported Hessian exact path operator annihilates affine probes", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- 0:5
  probes <- cbind(
    constant = rep(1, length(x)),
    affine = x,
    quadratic = x^2
  )
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    coordinates = matrix(x, ncol = 1),
    polynomial.probes = probes
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$diagnostics$nullity.estimate, 2L)
  expect_equal(op$diagnostics$summary$n.direction.labels, 2L)
  expect_gt(op$diagnostics$summary$n.rows, 0L)
  expect_gt(op$diagnostics$summary$n.dropped, 0L)
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0, tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0, tolerance = 1e-12)
  expect_gt(residuals$residual[residuals$probe == "quadratic"], 0)
})

test_that("transported Hessian exact grid operator matches only parallel directions", {
  graph <- make_grid_graph_weights(3, 3)
  coords <- graph$coordinates
  probes <- cbind(
    constant = rep(1, nrow(coords)),
    x = coords[, 1],
    y = coords[, 2],
    xy = coords[, 1] * coords[, 2]
  )
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    coordinates = coords,
    polynomial.probes = probes
  )

  expect_equal(op$diagnostics$summary$n.direction.labels, 4L)
  expect_true(all(op$row.table$direction.label == op$row.table$matched.direction.label))
  expect_false(any(op$row.table$direction.label == "axis1+" &
                   op$row.table$matched.direction.label == "axis2+"))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0, tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "x"], 0, tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "y"], 0, tolerance = 1e-12)
  expect_gt(residuals$residual[residuals$probe == "xy"], 0)
})

test_that("transported Hessian exposes row metadata and accepts explicit labels", {
  graph <- make_grid_graph_weights(3, 2)
  labels <- lapply(seq_along(graph$adj.list), function(i) {
    nbrs <- graph$adj.list[[i]]
    apply(graph$coordinates[nbrs, , drop = FALSE], 1, function(z) {
      delta <- z - graph$coordinates[i, ]
      if (abs(delta[1]) > 0) {
        paste0("horizontal", if (delta[1] > 0) "+" else "-")
      } else {
        paste0("vertical", if (delta[2] > 0) "+" else "-")
      }
    })
  })
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    direction.labels = labels
  )

  expect_true(all(c("row", "base.dart", "direction.dart", "matched.dart",
                    "direction.label", "matched.direction.label") %in%
                    names(op$row.table)))
  expect_true(all(op$row.table$direction.label == op$row.table$matched.direction.label))
  expect_true(all(grepl("horizontal|vertical", op$darts$direction.label)))
  expect_equal(ncol(op$hessian$matrix), length(graph$adj.list))
})

test_that("soft transported Hessian normalizes local embedding transport weights", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = graph$coordinates,
    soft.bandwidth = 0.1
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$rule, "local.embedding.soft")
  expect_equal(nrow(op$hessian$matrix), op$diagnostics$summary$n.rows)
  row.weight.sum <- rowsum(op$row.table$transport.weight, op$row.table$row)
  expect_equal(as.vector(row.weight.sum), rep(1, nrow(row.weight.sum)),
               tolerance = 1e-12)
  expect_equal(op$diagnostics$summary$n.dropped, 0L)
  expect_true(is.finite(op$diagnostics$summary$mean.transport.entropy))
  expect_true(is.finite(op$transport$scales$angle.scale))
})

test_that("soft transported Hessian strongly favors parallel grid directions", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = graph$coordinates,
    soft.bandwidth = 0.1
  )
  center.up.from.right <- subset(
    op$row.table,
    base.from == 5L & base.to == 6L &
      direction.from == 5L & direction.to == 8L
  )

  expect_gt(nrow(center.up.from.right), 1L)
  best <- center.up.from.right[which.max(center.up.from.right$transport.weight), ]
  expect_equal(best$matched.from, 6L)
  expect_equal(best$matched.to, 9L)
  expect_gt(best$transport.weight, 0.999)
  expect_lt(best$angle, 1e-10)
  expect_true(all(c("best.match.margin", "best.match.angle",
                    "best.length.relative.error", "best.transport.weight",
                    "reciprocal.best.match") %in% names(op$row.table)))
})

test_that("soft transported Hessian reports polynomial probe diagnostics", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- 0:5
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = matrix(x, ncol = 1),
    polynomial.probes = cbind(constant = rep(1, length(x)), affine = x),
    soft.bandwidth = 0.1
  )

  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_true(is.finite(residuals$residual[residuals$probe == "affine"]))
  expect_gte(residuals$residual[residuals$probe == "affine"], 0)
})

test_that("edge-angle hard transported Hessian preserves path affine probes", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "edge.angle.hard",
    coordinates = x,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1],
                              quadratic = x[, 1]^2),
    edge.angle.max.angle.difference = 1e-8,
    edge.angle.max.length.relative.error = 1e-8
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$rule, "edge.angle.hard")
  expect_true(all(op$row.table$match.status == "edge_angle_hard_matched"))
  expect_true(all(c("source.angle", "target.angle", "length.relative.error") %in%
                    names(op$row.table)))
  expect_equal(op$row.table$transport.weight, rep(1, nrow(op$row.table)))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-12)
  expect_gt(residuals$residual[residuals$probe == "quadratic"], 0)
})

test_that("edge-angle soft transported Hessian records mirror ambiguity on a grid", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "edge.angle.soft",
    coordinates = graph$coordinates,
    edge.angle.bandwidth = 0.1
  )
  center.up.from.right <- subset(
    op$row.table,
    base.from == 5L & base.to == 6L &
      direction.from == 5L & direction.to == 8L
  )

  expect_gt(nrow(center.up.from.right), 1L)
  row.weight.sum <- rowsum(op$row.table$transport.weight, op$row.table$row)
  expect_equal(as.vector(row.weight.sum), rep(1, nrow(row.weight.sum)),
               tolerance = 1e-12)
  best <- center.up.from.right[which.max(center.up.from.right$transport.weight), ]
  expect_equal(best$matched.from, 6L)
  expect_true(best$matched.to %in% c(3L, 9L))
  mirror <- subset(center.up.from.right, matched.to %in% c(3L, 9L))
  expect_equal(sum(mirror$transport.weight), 1, tolerance = 1e-8)
  expect_true(all(mirror$transport.weight > 0.49))
  expect_lt(best$angle, 1e-10)
})

test_that("soft transported Hessian confidence gates drop low-margin rows", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = graph$coordinates,
    soft.bandwidth = 0.2,
    min.match.margin = 0.5
  )

  expect_gt(nrow(op$dropped.table), 0L)
  expect_true("match_margin_below_threshold" %in% op$dropped.table$drop.reason)
  expect_true(all(op$row.table$match.accepted))
  expect_equal(op$diagnostics$summary$n.dropped, nrow(op$dropped.table))
  row.weight.sum <- rowsum(op$row.table$transport.weight, op$row.table$row)
  expect_equal(as.vector(row.weight.sum), rep(1, nrow(row.weight.sum)),
               tolerance = 1e-12)
})

test_that("soft transported Hessian supports adaptive local confidence gates", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = graph$coordinates,
    soft.bandwidth = 0.2,
    match.threshold.rule = "local.robust.z",
    min.best.score.z = 0,
    min.margin.z = 0,
    max.effective.match.fraction = 0.75
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_true(all(c("n.match.candidates", "best.score.local.quantile",
                    "best.score.robust.z", "margin.local.quantile",
                    "margin.robust.z", "effective.match.fraction") %in%
                    names(op$row.table)))
  expect_true(all(is.finite(op$row.table$effective.match.fraction)))
  expect_true(all(op$row.table$effective.match.fraction > 0))
  expect_lte(max(op$row.table$effective.match.fraction), 0.75)

  strict <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    coordinates = graph$coordinates,
    soft.bandwidth = 0.2,
    match.threshold.rule = "local.quantile",
    match.score.quantile = 0.25,
    match.margin.quantile = 1
  )
  expect_gt(nrow(strict$dropped.table), 0L)
  expect_true(any(strict$dropped.table$drop.reason %in%
                    c("best_score_quantile_above_threshold",
                      "margin_quantile_below_threshold")))
})

test_that("regression-gradient transported Hessian annihilates affine path probes", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1],
                              quadratic = x[, 1]^2),
    gradient.ridge = 0
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$rule, "regression.gradient")
  expect_equal(nrow(op$hessian$matrix), op$diagnostics$summary$n.rows)
  expect_true(all(c("component", "coefficient.vertex", "coefficient.value",
                    "local.design.rank.from", "local.design.rank.to") %in%
                    names(op$row.table)))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-12)
  expect_gt(residuals$residual[residuals$probe == "quadratic"], 0)
})

test_that("regression-gradient transported Hessian annihilates affine grid probes", {
  graph <- make_grid_graph_weights(3, 3)
  coords <- graph$coordinates
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = coords,
    polynomial.probes = cbind(constant = rep(1, nrow(coords)),
                              x = coords[, 1],
                              y = coords[, 2],
                              xy = coords[, 1] * coords[, 2]),
    gradient.ridge = 0
  )

  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "x"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "y"], 0,
               tolerance = 1e-12)
  expect_gt(residuals$residual[residuals$probe == "xy"], 0)
  expect_equal(ncol(op$hessian$matrix), length(graph$adj.list))
  expect_gt(nrow(op$transport$gradient.diagnostics), 0L)
})

test_that("regression-quadratic transported operator annihilates path quadratics", {
  graph <- make_path_graph_weights(rep(1, 7))
  x <- matrix(0:7, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.order = 3L,
    transport.rule = "regression.gradient",
    coordinates = x,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1],
                              quadratic = x[, 1]^2,
                              cubic = x[, 1]^3),
    gradient.ridge = 0,
    gradient.quadratic.disk.hops = 2L
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport.order, 3L)
  expect_equal(op$transport$scales$representation,
               "local_weighted_quadratic_least_squares")
  expect_equal(op$penalty$order, 2L)
  expect_equal(op$penalty$derivative.order, 3L)
  expect_equal(op$diagnostics$nullity.estimate, 3L)
  expect_true(all(c("hessian.component", "component.a", "component.b",
                    "gradient.quadratic.disk.hops") %in% names(op$row.table)))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-10)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-10)
  expect_equal(residuals$residual[residuals$probe == "quadratic"], 0,
               tolerance = 1e-10)
  expect_gt(residuals$residual[residuals$probe == "cubic"], 1e-3)
})

test_that("regression-quadratic transported operator annihilates grid quadratics", {
  graph <- make_grid_graph_weights(5, 5)
  coords <- graph$coordinates
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.order = 3L,
    transport.rule = "regression.gradient",
    coordinates = coords,
    polynomial.probes = cbind(constant = rep(1, nrow(coords)),
                              x = coords[, 1],
                              y = coords[, 2],
                              x2 = coords[, 1]^2,
                              xy = coords[, 1] * coords[, 2],
                              y2 = coords[, 2]^2,
                              x3 = coords[, 1]^3),
    gradient.ridge = 0,
    gradient.quadratic.disk.hops = 2L,
    gradient.quadratic.max.vertices = 25L
  )

  residuals <- op$diagnostics$polynomial.residuals$per.column
  for (probe in c("constant", "x", "y", "x2", "xy", "y2")) {
    expect_equal(residuals$residual[residuals$probe == probe], 0,
                 tolerance = 1e-9)
  }
  expect_gt(residuals$residual[residuals$probe == "x3"], 1e-3)
  expect_true(all(op$transport$gradient.diagnostics$rank == 5L))
  expect_equal(unique(op$row.table$hessian.component),
               c("H11", "H12", "H22"))
})

test_that("regression-quadratic transported operator is scoped to supplied coordinates", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  expect_error(
    transported.graph.hessian.operator(
      graph$adj.list,
      graph$weight.list,
      transport.order = 3L,
      transport.rule = "regression.gradient",
      coordinates = x,
      gradient.coordinate.method = "local.embedding"
    ),
    "supports only"
  )
  expect_error(
    transported.graph.hessian.operator(
      graph$adj.list,
      graph$weight.list,
      transport.order = 3L,
      transport.rule = "local.embedding.soft",
      coordinates = x
    ),
    "supports only"
  )
})

test_that("regression-gradient transported Hessian supports graph-derived local embeddings", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    gradient.coordinate.method = "local.embedding",
    gradient.embedding.method = "cmdscale",
    gradient.embedding.dim = 1L,
    gradient.disk.hops = 1L,
    gradient.max.vertices = 6L,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1],
                              quadratic = x[, 1]^2),
    gradient.ridge = 0
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$scales$gradient.coordinate.method, "local.embedding")
  expect_true(all(c("requested.method", "backend.used", "disk.size",
                    "edge.stress", "ambient.affine.residual",
                    "ambient.linear.gradient.residual.mean") %in%
                    names(op$transport$embedding.table)))
  expect_true(all(op$transport$embedding.table$requested.method == "cmdscale"))
  expect_equal(max(op$transport$embedding.table$ambient.affine.residual,
                   na.rm = TRUE), 0, tolerance = 1e-8)
  expect_equal(max(op$transport$embedding.table$ambient.linear.gradient.residual.mean,
                   na.rm = TRUE), 0, tolerance = 1e-8)
  expect_true(all(op$row.table$gradient.coordinate.method == "local.embedding"))
  expect_true(all(c("base.dart", "endpoint.role", "embedding.backend") %in%
                    names(op$transport$gradient.diagnostics)))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-10)
  expect_gt(residuals$residual[residuals$probe == "quadratic"], 0)
})

test_that("regression-gradient local embedding does not require supplied coordinates", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    gradient.coordinate.method = "local.embedding",
    gradient.embedding.method = "cmdscale",
    gradient.embedding.dim = 2L,
    gradient.disk.hops = 1L,
    gradient.max.vertices = 9L,
    gradient.ridge = 1e-8
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_gt(nrow(op$row.table), 0L)
  expect_gt(nrow(op$transport$embedding.table), 0L)
  expect_equal(ncol(op$hessian$matrix), length(graph$adj.list))
  expect_true(all(is.finite(op$row.table$coefficient.value)))
  expect_true(all(is.na(op$transport$embedding.table$ambient.affine.residual)))
  expect_true(all(is.na(op$transport$embedding.table$ambient.linear.gradient.residual.mean)))
})

test_that("regression-gradient local embedding supports adaptive chart selection", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    gradient.coordinate.method = "local.embedding",
    gradient.chart.selection = "adaptive",
    gradient.embedding.candidates = c("cmdscale", "mds.edge.kk"),
    gradient.disk.rule.candidates = "hops",
    gradient.disk.hops.candidates = 1:5,
    gradient.embedding.dim = 1L,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1],
                              quadratic = x[, 1]^2)
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$scales$gradient.chart.selection, "adaptive")
  expect_true(all(c("chart.selected", "chart.oracle", "chart.selection.score",
                    "disk.fraction", "disk.weighted.radius") %in%
                    names(op$transport$embedding.table)))
  selected <- subset(op$transport$embedding.table, chart.selected)
  expect_equal(nrow(selected), op$graph$n.darts)
  expect_true(all(selected$gradient.disk.hops %in% 1:5))
  expect_true(all(is.finite(selected$chart.selection.score)))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-8)
  expect_gt(residuals$residual[residuals$probe == "quadratic"], 0)
})

test_that("regression-gradient local embedding supports metric diameter fraction disks", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    gradient.coordinate.method = "local.embedding",
    gradient.embedding.method = "cmdscale",
    gradient.embedding.dim = 1L,
    gradient.disk.rule = "metric.diameter.fraction",
    gradient.disk.radius.fraction = 0.30,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1])
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_true(all(op$transport$embedding.table$gradient.disk.rule ==
                    "metric.diameter.fraction"))
  expect_true(all(op$transport$embedding.table$gradient.disk.radius.fraction ==
                    0.30))
  expect_true(all(is.finite(op$transport$embedding.table$gradient.disk.metric.radius)))
  expect_true(all(op$transport$embedding.table$disk.size >= 2L))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-8)
})

test_that("regression-gradient adaptive chart selection can select metric disks", {
  graph <- make_path_graph_weights(rep(1, 5))
  x <- matrix(0:5, ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    gradient.coordinate.method = "local.embedding",
    gradient.chart.selection = "adaptive",
    gradient.embedding.candidates = "cmdscale",
    gradient.disk.rule.candidates = "metric.diameter.fraction",
    gradient.disk.radius.fraction.candidates = c(0.20, 0.40),
    gradient.embedding.dim = 1L
  )

  selected <- subset(op$transport$embedding.table, chart.selected)
  expect_equal(nrow(selected), op$graph$n.darts)
  expect_true(all(selected$gradient.disk.rule == "metric.diameter.fraction"))
  expect_true(all(selected$gradient.disk.radius.fraction %in% c(0.20, 0.40)))
  expect_true(all(is.na(selected$gradient.disk.hops)))
})

test_that("regression-gradient local embedding supports metric local-scale disks", {
  graph <- make_path_graph_weights(c(1, 1, 2, 2, 3))
  x <- matrix(cumsum(c(0, c(1, 1, 2, 2, 3))), ncol = 1)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "regression.gradient",
    coordinates = x,
    gradient.coordinate.method = "local.embedding",
    gradient.embedding.method = "cmdscale",
    gradient.embedding.dim = 1L,
    gradient.disk.rule = "metric.local.scale",
    gradient.disk.local.scale.method = "median.incident.length",
    gradient.disk.local.scale.multiplier = 2,
    gradient.disk.min.vertices = 4L,
    polynomial.probes = cbind(constant = rep(1, nrow(x)),
                              affine = x[, 1])
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_true(all(op$transport$embedding.table$gradient.disk.rule ==
                    "metric.local.scale"))
  expect_true(all(op$transport$embedding.table$gradient.disk.local.scale.multiplier ==
                    2))
  expect_true(all(is.finite(op$transport$embedding.table$gradient.disk.local.scale)))
  expect_true(all(is.finite(op$transport$embedding.table$gradient.disk.metric.radius)))
  expect_true(all(op$transport$embedding.table$disk.size >= 4L))
  residuals <- op$diagnostics$polynomial.residuals$per.column
  expect_equal(residuals$residual[residuals$probe == "constant"], 0,
               tolerance = 1e-12)
  expect_equal(residuals$residual[residuals$probe == "affine"], 0,
               tolerance = 1e-8)
})

test_that("soft transported Hessian can use local cmdscale graph-disk embeddings", {
  graph <- make_grid_graph_weights(3, 3)
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    local.embedding.method = "cmdscale",
    local.embedding.dim = 2L,
    local.disk.hops = 1L,
    local.max.vertices = 9L,
    soft.bandwidth = 0.1
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_equal(op$transport$rule, "local.embedding.soft")
  expect_true(all(c("backend.used", "disk.size", "edge.stress") %in%
                    names(op$transport$embedding.table)))
  expect_true(all(op$transport$embedding.table$backend.used == "cmdscale"))
  expect_true(all(op$transport$embedding.table$disk.size <= 9L))
  row.weight.sum <- rowsum(op$row.table$transport.weight, op$row.table$row)
  expect_equal(as.vector(row.weight.sum), rep(1, nrow(row.weight.sum)),
               tolerance = 1e-12)
})

test_that("soft transported Hessian auto backend records graph-derived embeddings", {
  graph <- make_path_graph_weights(rep(1, 4))
  op <- transported.graph.hessian.operator(
    graph$adj.list,
    graph$weight.list,
    transport.rule = "local.embedding.soft",
    local.embedding.method = "auto",
    local.embedding.dim = 2L,
    local.disk.hops = 1L,
    local.max.vertices = 5L,
    soft.bandwidth = 0.2
  )

  expect_s3_class(op, "transported.graph.hessian.operator")
  expect_gt(nrow(op$transport$embedding.table), 0L)
  expect_true(all(is.finite(op$transport$embedding.table$edge.stress) |
                    is.na(op$transport$embedding.table$edge.stress)))
  expect_true(any(op$transport$embedding.table$backend.used %in%
                    c("grip.optimize.edge.kk.layout",
                      "grip.layout.weighted+grip.optimize.edge.kk.layout",
                      "cmdscale+grip.optimize.edge.kk.layout",
                      "grip.layout.weighted+grip.optimize.edge.isometric.layout",
                      "cmdscale+grip.optimize.edge.isometric.layout",
                      "cmdscale")))
  row.weight.sum <- rowsum(op$row.table$transport.weight, op$row.table$row)
  expect_equal(as.vector(row.weight.sum), rep(1, nrow(row.weight.sum)),
               tolerance = 1e-12)
})

test_that("higher-order small vertex CV returns selected lambda and diagnostics", {
  skip_if_not_installed("genlasso")
  skip_if_not_installed("Matrix")

  graph <- make_path_graph_weights(rep(1, 5))
  y <- sin(seq_len(6))
  lambda.grid <- c(1, 0.1, 0.01)

  for (ord in 1:2) {
    fit <- fit.graph.trend.filtering(
      graph$adj.list,
      graph$weight.list,
      y,
      order = ord,
      lambda.grid = lambda.grid,
      lambda.selection = "cv",
      weight.rule = "unit",
      foldid = c(1, 2, 3, 1, 2, 3),
      maxsteps = 100L
    )
    expect_s3_class(fit, "graph.trend.filtering.fit")
    expect_true(fit$lambda %in% lambda.grid)
    expect_equal(dim(fit$cv$fold.errors), c(3, 3))
    expect_equal(length(fit$fitted.values), length(y))
    expect_equal(fit$order, ord)
  }
})

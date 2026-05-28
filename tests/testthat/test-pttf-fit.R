make_pttf_fit_path_graph <- function(n) {
  adj <- vector("list", n)
  w <- vector("list", n)
  x <- seq(0, 1, length.out = n)
  for (i in seq_len(n - 1L)) {
    ell <- x[i + 1L] - x[i]
    adj[[i]] <- c(adj[[i]], i + 1L)
    w[[i]] <- c(w[[i]], ell)
    adj[[i + 1L]] <- c(adj[[i + 1L]], i)
    w[[i + 1L]] <- c(w[[i + 1L]], ell)
  }
  list(X = matrix(x, ncol = 1L),
       adj.list = lapply(adj, as.integer),
       weight.list = lapply(w, as.double))
}

make_pttf_fit_operator <- function(n = 24L) {
  graph <- make_pttf_fit_path_graph(n)
  geom <- pttf.geometry(
    graph$X,
    graph = "supplied",
    adj.list = graph$adj.list,
    weight.list = graph$weight.list,
    tangent.dim = 1L
  )
  op <- pttf.operator(
    geom,
    derivative.order = 3L,
    row.mass.rule = "none",
    row.normalize = "none"
  )
  list(graph = graph, geom = geom, op = op)
}

reconstruct_pttf_triplet <- function(triplet) {
  Matrix::sparseMatrix(
    i = triplet$i,
    j = triplet$j,
    x = triplet$x,
    dims = triplet$dim,
    giveCsparse = TRUE
  )
}

test_that("pttf.operator.filter.rows remaps triplets to filtered rows", {
  obj <- make_pttf_fit_operator(24L)
  keep <- which(obj$op$row.table$vertex > 3L &
                  obj$op$row.table$vertex <= 21L)
  filtered <- pttf.operator.filter.rows(obj$op, keep, reason = "test")

  expect_s3_class(filtered, "pttf_operator")
  expect_equal(dim(filtered$A), c(length(keep), ncol(obj$op$A)))
  expect_equal(as.matrix(filtered$B), as.matrix(Matrix::crossprod(filtered$A)),
               tolerance = 1e-12)
  expect_true(all(filtered$A.triplet$i %in% seq_len(nrow(filtered$A))))
  expect_equal(filtered$A.triplet$dim, dim(filtered$A))
  expect_equal(as.matrix(reconstruct_pttf_triplet(filtered$A.triplet)),
               as.matrix(filtered$A), tolerance = 1e-12)
  expect_equal(filtered$row.table$compact.row, keep)
  expect_equal(filtered$row.table$fit.row, seq_len(length(keep)))
  expect_equal(filtered$row.filter$dropped.nrow, nrow(obj$op$A) - length(keep))
})

test_that("PTTF L1 fixed-lambda fit respects drop-line-boundary policy", {
  skip_if_not_installed("genlasso")
  obj <- make_pttf_fit_operator(24L)
  x <- obj$graph$X[, 1L]
  y <- sin(2 * pi * x)

  fit <- fit.pttf.trend.filtering(
    operator = obj$op,
    y = y,
    penalty = "l1",
    lambda.selection = "fixed",
    lambda.grid = 0.01,
    operator.row.policy = "drop.line.boundary",
    line.order = seq_along(y)
  )

  expect_s3_class(fit, "pttf.trend.filtering.fit")
  expect_length(fit$fitted.values, length(y))
  expect_true(all(is.finite(fit$fitted.values)))
  expect_equal(fit$row.policy$name, "drop.line.boundary")
  expect_equal(nrow(fit$fit.operator$A), nrow(obj$op$A) - 6L)
  expect_true(all(fit$fit.operator$A.triplet$i %in%
                    seq_len(nrow(fit$fit.operator$A))))
  expect_equal(fit$lambda.boundary, "none")
})

test_that("PTTF CV uses supplied foldid deterministically and records metadata", {
  skip_if_not_installed("genlasso")
  obj <- make_pttf_fit_operator(18L)
  x <- obj$graph$X[, 1L]
  y <- cos(2 * pi * x)
  foldid <- rep(1:3, length.out = length(y))

  fit1 <- fit.pttf.trend.filtering(
    operator = obj$op,
    y = y,
    penalty = "l1",
    lambda.selection = "cv",
    lambda.grid = c(1e-2, 1e-3),
    foldid = foldid,
    cv.loss = "mse",
    selection = "min",
    operator.row.policy = "all",
    maxsteps = 500L
  )
  fit2 <- fit.pttf.trend.filtering(
    operator = obj$op,
    y = y,
    penalty = "l1",
    lambda.selection = "cv",
    lambda.grid = c(1e-2, 1e-3),
    foldid = foldid,
    cv.loss = "mse",
    selection = "min",
    operator.row.policy = "all",
    maxsteps = 500L
  )

  expect_equal(fit1$fitted.values, fit2$fitted.values, tolerance = 1e-12)
  expect_equal(fit1$lambda, fit2$lambda)
  expect_equal(fit1$foldid, foldid)
  expect_equal(fit1$fold.source, "supplied")
  expect_equal(fit1$n.lambda, 80L)
  expect_equal(fit1$cv.loss, "mse")
  expect_equal(fit1$selection, "min")
  expect_true(all(c("mean.error", "se", "selected.idx") %in%
                    names(fit1$cv)))
})

test_that("PTTF L2 fit supports supplied folds and finite predictions", {
  obj <- make_pttf_fit_operator(18L)
  x <- obj$graph$X[, 1L]
  y <- x^2 + 0.1 * sin(2 * pi * x)
  foldid <- rep(1:3, length.out = length(y))

  fit <- fit.pttf.trend.filtering(
    operator = obj$op,
    y = y,
    penalty = "l2",
    lambda.selection = "cv",
    lambda.grid = c(1, 0.1, 0.01),
    foldid = foldid,
    operator.row.policy = "all"
  )

  expect_s3_class(fit, "pttf.trend.filtering.fit")
  expect_true(all(is.finite(fit$fitted.values)))
  expect_equal(fit$foldid, foldid)
  expect_equal(fit$fold.source, "supplied")
  expect_equal(fit$cv$loss, "mse")
})

test_that("PTTF fit validates row-policy and response contracts", {
  obj <- make_pttf_fit_operator(12L)
  y <- seq_len(12L)
  expect_error(
    fit.pttf.trend.filtering(
      operator = obj$op,
      y = y,
      lambda.selection = "fixed",
      lambda.grid = 0.01,
      operator.row.policy = "drop.line.boundary"
    ),
    "line.order"
  )
  expect_error(
    fit.pttf.trend.filtering(
      operator = obj$op,
      y = cbind(y, y),
      lambda.selection = "fixed",
      lambda.grid = 0.01
    ),
    "one response"
  )
})

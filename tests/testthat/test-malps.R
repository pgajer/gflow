test_that("fit.malps reproduces constant, affine, and quadratic signals", {
  x1 <- seq(-1, 1, length.out = 7L)
  x2 <- seq(-0.6, 0.6, length.out = 6L)
  X <- as.matrix(expand.grid(x1, x2))

  y0 <- rep(3.5, nrow(X))
  fit0 <- fit.malps(
    X, y0,
    degree = 0L,
    support.type = "knn",
    support.size = nrow(X),
    min.support = 1L,
    kernel = "gaussian"
  )
  expect_s3_class(fit0, "malps")
  expect_equal(fit0$fitted.values, y0, tolerance = 1e-10)

  y1 <- 1 + 2 * X[, 1L] - 0.5 * X[, 2L]
  fit1 <- fit.malps(
    X, y1,
    degree = 1L,
    support.type = "knn",
    support.size = nrow(X),
    min.support = 4L,
    kernel = "gaussian"
  )
  expect_equal(fit1$fitted.values, y1, tolerance = 1e-9)

  y2 <- 0.7 + 1.2 * X[, 1L] - 0.8 * X[, 2L] +
    0.3 * X[, 1L]^2 - 0.4 * X[, 1L] * X[, 2L] + 0.9 * X[, 2L]^2
  fit2 <- fit.malps(
    X, y2,
    degree = 2L,
    support.type = "knn",
    support.size = nrow(X),
    min.support = 8L,
    kernel = "gaussian"
  )
  expect_equal(fit2$fitted.values, y2, tolerance = 1e-8)
  expect_true(all(is.finite(fit2$diagnostics$condition.number)))
  expect_false(any(fit2$diagnostics$fallback.used))
})

test_that("fit.malps records observed anchors and kNN self-inclusion", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 4L,
    min.support = 4L,
    kernel = "epanechnikov"
  )

  expect_equal(fit$anchor.index, seq_len(nrow(X)))
  expect_equal(fit$anchors, X)
  expect_true(all(vapply(seq_along(fit$supports), function(i) {
    fit$anchor.index[i] %in% fit$supports[[i]]$index
  }, logical(1L))))
  expect_true(all(fit$diagnostics$support.size == 4L))
  expect_true(all(fit$diagnostics$positive.weight.support.size == 4L))
  expect_true(all(rowSums(fit$averaging.weights) > 0))
  expect_equal(rowSums(fit$averaging.weights), rep(1, nrow(X)),
               tolerance = 1e-12)
})

test_that("fit.malps supports model-quality averaging weights", {
  X <- matrix(c(0, 0.1, 0.2, 1), ncol = 1L)
  y <- c(0, 1, 2, 10)

  default <- fit.malps(
    X, y,
    anchor.index = c(1L, 4L),
    degree = 0L,
    support.type = "fixed.radius",
    radius = 0.25,
    min.support = 1L,
    kernel = "triangular"
  )
  weighted <- fit.malps(
    X, y,
    anchor.index = c(1L, 4L),
    degree = 0L,
    support.type = "fixed.radius",
    radius = 0.25,
    min.support = 1L,
    kernel = "triangular",
    model.weight.rule = "support"
  )

  expect_equal(default$model.weight.rule, "none")
  expect_equal(default$model.weights, rep(1, 2L))
  expect_equal(weighted$model.weight.rule, "support")
  expect_equal(weighted$diagnostics$model.weight.rule, "support")
  expect_gt(weighted$model.weights[1L], weighted$model.weights[2L])
  expect_equal(weighted$diagnostics$model.weight, weighted$model.weights)
  expect_true(all(c("model.weight.raw", "model.weight.floor",
                    "model.weight.condition", "model.weight.support",
                    "model.weight.boundary",
                    "model.weight.source.condition.number",
                    "model.weight.fixed.from.original") %in%
                    names(weighted$diagnostics)))
  expect_false(weighted$diagnostics$model.weight.fixed.from.original)
  expect_equal(rowSums(weighted$averaging.weights), rep(1, nrow(X)),
               tolerance = 1e-12)

  refit <- refit.malps(weighted, y = y + 1)
  expect_equal(refit$model.weights, weighted$model.weights, tolerance = 1e-12)
  expect_equal(refit$averaging.weights, weighted$averaging.weights,
               tolerance = 1e-12)
  expect_true(refit$diagnostics$model.weight.fixed.from.original)
  expect_equal(refit$diagnostics$model.weight.source.condition.number,
               weighted$diagnostics$model.weight.source.condition.number,
               tolerance = 1e-12)
})

test_that("fit.malps quality weighting records component multipliers", {
  X <- matrix(c(0, 0.05, 0.1, 0.15, 0.6, 0.9, 0.95, 1), ncol = 1L)
  y <- as.numeric(X[, 1L]^2)

  fit <- fit.malps(
    X, y,
    anchor.index = c(1L, 4L, 5L, 7L),
    degree = 1L,
    support.type = "fixed.radius",
    radius = 0.35,
    min.support = 2L,
    kernel = "triangular",
    model.weight.rule = "quality"
  )

  condition.weight <- 1 / pmax(fit$diagnostics$model.condition.number, 1)
  condition.weight[!is.finite(condition.weight) | condition.weight < 0] <- 0
  support.weight <- fit$diagnostics$positive.weight.support.size /
      max(fit$diagnostics$positive.weight.support.size)
  boundary.weight <- 1 / (1 + pmax(fit$diagnostics$boundary.score, 0))
  raw.weight <- condition.weight * support.weight * boundary.weight
  expected.weight <- raw.weight / stats::median(raw.weight[raw.weight > 0])
  expected.weight[expected.weight <= 0 | !is.finite(expected.weight)] <-
      .Machine$double.eps

  expect_equal(fit$model.weight.rule, "quality")
  expect_equal(fit$diagnostics$model.weight.condition, condition.weight,
               tolerance = 1e-12)
  expect_equal(fit$diagnostics$model.weight.support, support.weight,
               tolerance = 1e-12)
  expect_equal(fit$diagnostics$model.weight.boundary, boundary.weight,
               tolerance = 1e-12)
  expect_equal(fit$diagnostics$model.weight.raw, raw.weight,
               tolerance = 1e-12)
  expect_equal(fit$model.weights, expected.weight, tolerance = 1e-12)
  expect_equal(fit$diagnostics$model.weight, fit$model.weights,
               tolerance = 1e-12)
  expect_equal(rowSums(fit$averaging.weights), rep(1, nrow(X)),
               tolerance = 1e-12)
})

test_that("model-quality zero raw weights preserve prediction support membership", {
  X <- matrix(c(0, 0.1, 0.2, 0.8, 0.9, 1), ncol = 1L)
  y <- as.numeric(X[, 1L]^2)

  base <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "fixed.radius",
    radius = 0.16,
    min.support = 2L,
    kernel = "triangular",
    model.weight.rule = "none"
  )
  quality <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "fixed.radius",
    radius = 0.16,
    min.support = 2L,
    kernel = "triangular",
    model.weight.rule = "quality"
  )

  expect_true(any(!is.finite(quality$diagnostics$model.condition.number)))
  expect_true(any(quality$diagnostics$model.weight.raw == 0))
  expect_true(all(quality$model.weights > 0))
  expect_equal(
    lapply(quality$prediction.supports, `[[`, "index"),
    lapply(base$prediction.supports, `[[`, "index")
  )
  expect_equal(quality$diagnostics$coverage.count,
               base$diagnostics$coverage.count)
  expect_equal(quality$diagnostics$model.weight.floor,
               .Machine$double.eps)
})

test_that("fit.malps kNN self-inclusion survives duplicate coordinates", {
  X <- matrix(c(0, 0, 1, 2), ncol = 1L)
  y <- c(1, 2, 3, 4)
  fit <- fit.malps(
    X, y,
    degree = 0L,
    support.type = "knn",
    support.size = 1L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_equal(fit$supports[[2L]]$index, 2L)
  expect_true(all(fit$diagnostics$support.size == 1L))
  expect_true(fit$diagnostics$has.duplicate.rows)
  expect_equal(fit$diagnostics$n.duplicate.rows, 2L)
  expect_equal(fit$diagnostics$duplicate.groups, list(c(1L, 2L)))
  expect_equal(fit$duplicate.action, "keep")
})

test_that("fit.malps can reject duplicate coordinate rows explicitly", {
  X <- matrix(c(0, 0, 1, 2), ncol = 1L)
  y <- c(1, 2, 3, 4)

  expect_error(
    fit.malps(
      X, y,
      degree = 0L,
      support.type = "knn",
      support.size = 1L,
      min.support = 1L,
      duplicate.action = "error"
    ),
    "duplicate coordinate row"
  )
})

test_that("fit.malps supports restricted observed anchors and detects no coverage", {
  X <- matrix(seq(0, 1, length.out = 9L), ncol = 1L)
  y <- X[, 1L]^2

  ok <- fit.malps(
    X, y,
    anchor.index = c(1L, 5L, 9L),
    degree = 1L,
    support.type = "fixed.radius",
    radius = 0.51,
    min.support = 2L,
    kernel = "triangular"
  )
  expect_equal(ok$anchor.index, c(1L, 5L, 9L))
  expect_true(all(ok$diagnostics$coverage.count > 0L))
  expect_equal(ok$diagnostics$n.no.coverage, 0L)

  expect_error(
    fit.malps(
      X, y,
      anchor.index = 1L,
      degree = 1L,
      support.type = "fixed.radius",
      radius = 0.2,
      min.support = 2L,
      kernel = "triangular"
    ),
    "uncovered.index"
  )
})

test_that("fit.malps handles adaptive-radius support and solver fallback", {
  X <- cbind(
    rep(seq(0, 1, length.out = 6L), each = 2L),
    rep(c(0, 0.1), times = 6L)
  )
  y <- 1 + X[, 1L] + X[, 2L]

  fit <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "adaptive.radius",
    min.support = 6L,
    support.buffer = 0L,
    kernel = "tricube",
    local.solver = "auto"
  )

  expect_s3_class(fit, "malps")
  expect_true(all(is.finite(fit$fitted.values)))
  expect_true(all(fit$diagnostics$positive.weight.support.size >= 6L))
  expect_true(any(fit$diagnostics$solver.used == "svd") ||
                any(fit$diagnostics$fallback.used))
  expect_true(any(fit$diagnostics$rank.deficient))
  expect_true(any(is.infinite(fit$diagnostics$condition.number)))
})

test_that("predict.malps returns training fits and predicts new coordinates", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- 1 + 2 * X[, 1L]
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "fixed.radius",
    radius = 2,
    min.support = 2L,
    kernel = "gaussian"
  )
  newdata <- matrix(c(0.1, 0.35, 0.9), ncol = 1L)

  expect_equal(predict(fit), fit$fitted.values)
  expect_equal(as.numeric(predict(fit, newdata = newdata)),
               1 + 2 * newdata[, 1L],
               tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, newdata = newdata[, 1L])),
               1 + 2 * newdata[, 1L],
               tolerance = 1e-10)
  expect_error(predict(fit, newdata = matrix(0, ncol = 2L)),
               "same number of columns")
  expect_error(predict.malps(list()), "malps")
})

test_that("predict.malps reports uncovered new coordinates", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]
  fit <- fit.malps(
    X, y,
    degree = 0L,
    support.type = "fixed.radius",
    radius = 0.25,
    min.support = 1L,
    kernel = "epanechnikov"
  )
  newdata <- matrix(c(0.5, 3), ncol = 1L)

  expect_error(predict(fit, newdata = newdata), "no positive averaging coverage")
  pred <- predict(fit, newdata = newdata, allow.incomplete = TRUE)
  expect_true(is.finite(pred[1L]))
  expect_true(is.na(pred[2L]))
})

test_that("fit.malps supports local PCA chart coordinates", {
  t <- seq(-1, 1, length.out = 9L)
  X <- cbind(t, 2 * t)
  y <- 1 + 3 * t
  fit <- fit.malps(
    X, y,
    degree = 1L,
    coordinate.method = "local.pca",
    chart.dim = 1L,
    support.type = "knn",
    support.size = nrow(X),
    min.support = 2L,
    kernel = "gaussian"
  )
  newdata <- cbind(c(-0.5, 0.25), 2 * c(-0.5, 0.25))

  expect_equal(fit$chart.dim, 1L)
  expect_equal(fit$diagnostics$required.support, 2L)
  expect_equal(fit$fitted.values, y, tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, newdata = newdata)),
               1 + 3 * newdata[, 1L],
               tolerance = 1e-10)
  expect_equal(ncol(fit$charts[[1L]]$basis), 1L)
  expect_equal(fit$coordinate.method, "local.pca")
})

test_that("fit.malps auto solver handles rank-deficient local PCA charts", {
  t <- seq(-1, 1, length.out = 11L)
  X <- cbind(t, 2 * t, 0.25 * t)
  y <- 1 - 2 * t

  fit <- fit.malps(
    X, y,
    degree = 1L,
    coordinate.method = "local.pca",
    support.type = "knn",
    support.size = 7L,
    min.support = 4L,
    kernel = "gaussian"
  )

  expect_equal(fit$chart.dim, ncol(X))
  expect_true(all(is.finite(fit$fitted.values)))
  expect_true(any(fit$diagnostics$solver.used == "svd"))
  expect_true(any(fit$diagnostics$fallback.used))
  expect_true(any(!is.finite(fit$diagnostics$condition.number)))
  expect_equal(fit$fitted.values, y, tolerance = 1e-9)
})

test_that("fit.malps local PCA charts work with CV and refit", {
  t <- seq(-1, 1, length.out = 15L)
  X <- cbind(t, 0.5 * t)
  y <- 2 - t
  fit <- fit.malps(
    X, y,
    coordinate.method = "local.pca",
    chart.dim = 1L,
    support.selection = "cv",
    support.type = "knn",
    degree.grid = 1L,
    support.grid = c(5L, 15L),
    min.support.grid = 2L,
    kernel.grid = "gaussian",
    cv.folds = 3L,
    cv.seed = 17L
  )
  refit <- refit.malps(fit, y = y)

  expect_equal(fit$selection$selected.params$chart.dim, 1L)
  expect_true("chart.dim" %in% names(fit$selection$candidate.table))
  expect_true(all(fit$selection$candidate.table$chart.dim == 1L))
  expect_equal(fit$chart.dim, 1L)
  expect_equal(refit$fitted.values, fit$fitted.values, tolerance = 1e-12)
  expect_equal(predict(refit, newdata = X[1:3, , drop = FALSE]),
               predict(fit, newdata = X[1:3, , drop = FALSE]),
               tolerance = 1e-12)
})

test_that("fit.malps uses supplied graph geodesic support distances", {
  X <- matrix(c(0, 100, 1, 2), ncol = 1L)
  y <- c(0, 1, 2, 3)
  adj <- list(
    2L,
    c(1L, 3L),
    c(2L, 4L),
    3L
  )
  weights <- list(
    1,
    c(1, 1),
    c(1, 1),
    1
  )

  fit <- fit.malps(
    X, y,
    adj.list = adj,
    weight.list = weights,
    support.metric = "graph.geodesic",
    degree = 0L,
    support.type = "knn",
    support.size = 2L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_equal(fit$support.metric, "graph.geodesic")
  expect_equal(fit$graph$source, "supplied")
  expect_equal(fit$supports[[1L]]$index, c(1L, 2L))
  expect_equal(fit$supports[[1L]]$distance, c(0, 1), tolerance = 1e-12)
  expect_true(all(is.finite(fit$fitted.values)))

  coordinate.fit <- fit.malps(
    X, y,
    support.metric = "coordinates",
    degree = 0L,
    support.type = "knn",
    support.size = 2L,
    min.support = 1L,
    kernel = "gaussian"
  )
  expect_equal(coordinate.fit$supports[[1L]]$index, c(1L, 3L))
})

test_that("fit.malps auto support metric uses graph input when supplied", {
  X <- matrix(c(0, 100, 1, 2), ncol = 1L)
  y <- c(0, 1, 2, 3)
  adj <- list(
    2L,
    c(1L, 3L),
    c(2L, 4L),
    3L
  )
  weights <- list(
    1,
    c(1, 1),
    c(1, 1),
    1
  )

  fit <- fit.malps(
    X, y,
    adj.list = adj,
    weight.list = weights,
    degree = 0L,
    support.type = "knn",
    support.size = 2L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_equal(fit$support.metric, "graph.geodesic")
  expect_equal(fit$supports[[1L]]$index, c(1L, 2L))
  expect_error(
    predict(fit, newdata = matrix(c(0, 1), ncol = 1L)),
    "graph.geodesic"
  )
})

test_that("fit.malps accepts supported gflow graph objects for geodesic supports", {
  X <- matrix(c(0, 1, 2, 3), ncol = 1L)
  y <- c(0, 1, 2, 3)
  graph <- create.rknn.graph(
    X,
    type = "fixed",
    radius = 1.1,
    connect.components = TRUE
  )

  fit <- fit.malps(
    X, y,
    graph = graph,
    degree = 0L,
    support.type = "knn",
    support.size = 2L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_equal(fit$support.metric, "graph.geodesic")
  expect_equal(fit$graph$source, "graph.final")
  expect_equal(fit$graph$stage, "final")
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.malps coordinate metric ignores unused invalid graph payloads", {
  X <- matrix(c(0, 1, 2), ncol = 1L)
  y <- c(0, 1, 4)
  adj <- list(2L, c(1L, 3L), 2L)

  fit <- fit.malps(
    X, y,
    adj.list = adj,
    support.metric = "coordinates",
    degree = 0L,
    support.type = "knn",
    support.size = 2L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_equal(fit$support.metric, "coordinates")
  expect_null(fit$graph)
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.malps graph geodesic support works with CV diagnostics", {
  X <- matrix(seq_len(6), ncol = 1L)
  y <- as.numeric(seq_len(6))
  adj <- vector("list", 6L)
  weights <- vector("list", 6L)
  for (i in seq_len(6L)) {
    nbrs <- c(i - 1L, i + 1L)
    nbrs <- nbrs[nbrs >= 1L & nbrs <= 6L]
    adj[[i]] <- nbrs
    weights[[i]] <- rep(1, length(nbrs))
  }

  fit <- fit.malps(
    X, y,
    adj.list = adj,
    weight.list = weights,
    support.metric = "graph.geodesic",
    support.selection = "cv",
    degree.grid = 0L,
    support.type = "knn",
    support.grid = c(2L, 3L),
    min.support.grid = 1L,
    kernel.grid = "gaussian",
    cv.folds = 3L,
    cv.seed = 42L
  )

  expect_equal(fit$selection$selected.params$support.metric, "graph.geodesic")
  expect_true("support.metric" %in% names(fit$selection$candidate.table))
  expect_true(all(fit$selection$candidate.table$support.metric ==
                    "graph.geodesic"))
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit.malps CV records model-quality averaging rule", {
  X <- matrix(seq(0, 1, length.out = 18L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])

  fit <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "adaptive.radius",
    degree.grid = 0L,
    min.support.grid = c(3L, 5L),
    kernel.grid = "gaussian",
    model.weight.rule = "boundary",
    cv.folds = 3L,
    cv.seed = 124L
  )

  expect_equal(fit$model.weight.rule, "boundary")
  expect_equal(fit$selection$selected.params$model.weight.rule, "boundary")
  expect_true("model.weight.rule" %in% names(fit$selection$candidate.table))
  expect_true(all(fit$selection$candidate.table$model.weight.rule == "boundary"))
  expect_true(all(is.finite(fit$model.weights)))
})

test_that("fit.malps validates graph geodesic support inputs", {
  X <- matrix(seq_len(3), ncol = 1L)
  y <- as.numeric(seq_len(3))
  adj <- list(2L, c(1L, 3L), 2L)
  weights <- list(1, c(1, 1), 1)

  expect_error(
    fit.malps(
      X, y,
      support.metric = "graph.geodesic",
      degree = 0L,
      support.type = "knn",
      support.size = 2L,
      min.support = 1L
    ),
    "requires graph input"
  )
  expect_error(
    fit.malps(
      X, y,
      adj.list = adj,
      support.metric = "graph.geodesic",
      degree = 0L,
      support.type = "knn",
      support.size = 2L,
      min.support = 1L
    ),
    "Both adj.list and weight.list"
  )
  bad.weights <- weights
  bad.weights[[1L]] <- 2
  expect_error(
    fit.malps(
      X, y,
      adj.list = adj,
      weight.list = bad.weights,
      support.metric = "graph.geodesic",
      degree = 0L,
      support.type = "knn",
      support.size = 2L,
      min.support = 1L
    ),
    "Reciprocal edge weights mismatch"
  )
})

test_that("fit.malps supports off-sample coordinate anchors", {
  X <- matrix(seq(0, 1, length.out = 9L), ncol = 1L)
  y <- 1 + 2 * X[, 1L]
  anchor.coordinates <- matrix(c(-0.25, 0.5, 1.25), ncol = 1L)

  fit <- fit.malps(
    X, y,
    anchor.coordinates = anchor.coordinates,
    degree = 1L,
    support.type = "knn",
    support.size = nrow(X),
    min.support = 2L,
    kernel = "gaussian"
  )

  expect_true(all(is.na(fit$anchor.index)))
  expect_equal(fit$anchors, anchor.coordinates)
  expect_equal(fit$anchor.labels, paste0("offsample.", 1:3))
  expect_equal(colnames(fit$averaging.weights),
               paste0("anchor.offsample.", 1:3))
  expect_equal(fit$fitted.values, y, tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, newdata = matrix(c(0.2, 0.8), ncol = 1L))),
               c(1.4, 2.6), tolerance = 1e-10)
})

test_that("fit.malps off-sample coordinate anchors work with CV and refit", {
  X <- matrix(seq(0, 1, length.out = 15L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  anchors <- matrix(seq(0.1, 0.9, length.out = 4L), ncol = 1L)

  fit <- fit.malps(
    X, y,
    anchor.coordinates = anchors,
    support.selection = "cv",
    support.type = "knn",
    degree.grid = 1L,
    support.grid = c(5L, 9L),
    min.support.grid = 2L,
    kernel.grid = "gaussian",
    cv.folds = 3L,
    cv.seed = 99L
  )
  refit <- refit.malps(fit, y = y)

  expect_equal(fit$selection$selected.params$support.metric, "coordinates")
  expect_true(all(is.na(fit$anchor.index)))
  expect_equal(refit$anchors, fit$anchors)
  expect_equal(refit$anchor.labels, fit$anchor.labels)
  expect_equal(refit$fitted.values, fit$fitted.values, tolerance = 1e-12)
})

test_that("fit.malps rejects unsupported off-sample graph-geodesic anchors", {
  X <- matrix(seq_len(3), ncol = 1L)
  y <- as.numeric(seq_len(3))
  adj <- list(2L, c(1L, 3L), 2L)
  weights <- list(1, c(1, 1), 1)

  expect_error(
    fit.malps(
      X, y,
      adj.list = adj,
      weight.list = weights,
      anchor.coordinates = matrix(1.5, ncol = 1L),
      support.metric = "graph.geodesic",
      degree = 0L,
      support.type = "knn",
      support.size = 2L,
      min.support = 1L
    ),
    "anchor.coordinates is not supported"
  )
})

test_that("refit.malps reuses the fixed support profile", {
  X <- matrix(seq(0, 1, length.out = 12L), ncol = 1L)
  y1 <- sin(2 * pi * X[, 1L])
  y2 <- 1 + X[, 1L] + 0.25 * X[, 1L]^2
  fit <- fit.malps(
    X, y1,
    degree = 2L,
    support.type = "knn",
    support.size = 7L,
    min.support = 5L,
    kernel = "gaussian"
  )

  refit <- refit.malps(fit, y = y2)
  direct <- fit.malps(
    X, y2,
    degree = fit$degree,
    support.type = fit$support.type,
    support.size = fit$selection$selected.params$support.size,
    min.support = fit$selection$selected.params$min.support,
    radius = fit$selection$selected.params$radius,
    kernel = fit$kernel,
    local.solver = fit$local.solver,
    normal.equations.max.condition = fit$normal.equations.max.condition
  )

  expect_s3_class(refit, "malps")
  expect_equal(refit$fitted.values, direct$fitted.values, tolerance = 1e-12)
  expect_equal(refit$averaging.weights, fit$averaging.weights,
               tolerance = 1e-12)
  expect_equal(refit$supports, fit$supports)
  expect_equal(refit$selection, fit$selection)
  expect_equal(refit$refit$reuse.selection, TRUE)
})

test_that("refit.malps reuses selected CV parameters without retuning", {
  X <- matrix(seq(0, 1, length.out = 24L), ncol = 1L)
  y1 <- sin(2 * pi * X[, 1L])
  y2 <- cos(2 * pi * X[, 1L])
  fit <- fit.malps(
    X, y1,
    support.selection = "cv",
    support.type = "adaptive.radius",
    degree.grid = c(0L, 1L),
    min.support.grid = c(5L, 8L),
    kernel.grid = "gaussian",
    cv.folds = 3L,
    cv.seed = 11L
  )

  refit <- refit.malps(fit, y = y2)
  params <- fit$selection$selected.params
  direct <- fit.malps(
    X, y2,
    anchor.index = fit$anchor.index,
    degree = params$degree,
    support.type = params$support.type,
    support.size = params$support.size,
    radius = params$radius,
    min.support = params$min.support,
    kernel = params$kernel,
    support.selection = "fixed",
    local.solver = fit$local.solver,
    normal.equations.max.condition = fit$normal.equations.max.condition
  )

  expect_equal(refit$fitted.values, direct$fitted.values, tolerance = 1e-12)
  expect_equal(refit$selection$method, "cv")
  expect_equal(refit$selection$selected.params, params)
  expect_equal(refit$cv, fit$cv)
})

test_that("refit.malps supports weighted coefficient refits", {
  X <- matrix(c(0, 1, 2), ncol = 1L)
  y <- c(0, 10, 10)
  fit <- fit.malps(
    X, y,
    degree = 0L,
    support.type = "knn",
    support.size = 3L,
    min.support = 1L,
    kernel = "gaussian"
  )

  weighted <- refit.malps(fit, weights = c(1, 0, 0))
  expect_equal(weighted$fitted.values, rep(0, nrow(X)), tolerance = 1e-12)
  expect_equal(weighted$refit$weights, c(1, 0, 0))
  expect_equal(weighted$case.weights, c(1, 0, 0))
  expect_true(all(weighted$diagnostics$effective.fit.support.size == 1L))
  expect_equal(weighted$averaging.weights, fit$averaging.weights,
               tolerance = 1e-12)
  expect_true(all(vapply(weighted$fit.weights, function(w) sum(w > 0) == 1L,
                         logical(1L))))

  unweighted <- refit.malps(fit, weights = rep(1, nrow(X)))
  expect_equal(unweighted$fitted.values, fit$fitted.values, tolerance = 1e-12)
})

test_that("refit.malps rejects globally nonzero weights that zero a local support", {
  X <- matrix(seq(0, 1, length.out = 4L), ncol = 1L)
  y <- X[, 1L]
  fit <- fit.malps(
    X, y,
    degree = 0L,
    support.type = "knn",
    support.size = 1L,
    min.support = 1L,
    kernel = "gaussian"
  )

  expect_error(
    refit.malps(fit, weights = c(1, 0, 0, 0)),
    "anchor 2 has no positive effective fit weights"
  )
})

test_that("refit.malps validates current contracts", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 4L,
    min.support = 4L
  )

  expect_equal(refit.malps(fit)$fitted.values, fit$fitted.values,
               tolerance = 1e-12)
  expect_error(refit.malps(list()), "malps")
  expect_error(refit.malps(fit, y = y[-1L]), "length nrow")
  expect_error(refit.malps(fit, weights = rep(1, nrow(X) - 1L)),
               "length nrow")
  expect_error(refit.malps(fit, weights = rep(NA_real_, nrow(X))),
               "finite")
  expect_error(refit.malps(fit, weights = rep(-1, nrow(X))),
               "nonnegative")
  expect_error(refit.malps(fit, weights = rep(0, nrow(X))),
               "positive")
  expect_error(refit.malps(fit, reuse.selection = FALSE), "reuse.selection")
  expect_error(
    refit.malps(fit, refit.local.coefficients = FALSE),
    "refit.local.coefficients"
  )
})

test_that("fit.malps validates current coordinate contracts", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]

  expect_error(fit.malps(X, y[-1L]), "length nrow")
  expect_error(fit.malps(X, y, degree = 3L), "degree")
  expect_error(fit.malps(X, y, anchor.index = 0L), "anchor.index")
  expect_error(fit.malps(X, y, chart.dim = 1L), "chart.dim")
  expect_error(
    fit.malps(cbind(X, X), y, coordinate.method = "local.pca", chart.dim = 3L),
    "chart.dim"
  )
  expect_error(
    fit.malps(X, y, support.type = "knn", support.size = 2L,
              min.support = 4L),
    "support.size"
  )
  expect_error(
    fit.malps(X, y, support.type = "fixed.radius"),
    "radius"
  )
})

test_that("fit.malps CV uses supplied foldid deterministically", {
  X <- matrix(seq(0, 1, length.out = 24L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L]) + 0.1 * X[, 1L]
  foldid <- rep(1:4, length.out = nrow(X))

  fit1 <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "knn",
    degree.grid = c(0L, 1L),
    support.grid = c(6L, 9L),
    min.support.grid = 4L,
    kernel.grid = "gaussian",
    foldid = foldid,
    cv.repeats = 3L,
    cv.loss = "mse"
  )
  fit2 <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "knn",
    degree.grid = c(0L, 1L),
    support.grid = c(6L, 9L),
    min.support.grid = 4L,
    kernel.grid = "gaussian",
    foldid = foldid,
    cv.repeats = 3L,
    cv.loss = "mse"
  )

  expect_s3_class(fit1, "malps")
  expect_equal(fit1$fitted.values, fit2$fitted.values, tolerance = 1e-12)
  expect_equal(fit1$selection$method, "cv")
  expect_equal(fit1$selection$foldid.source, "supplied")
  expect_equal(fit1$selection$cv.folds, 4L)
  expect_equal(fit1$selection$cv.folds.requested, 5L)
  expect_equal(fit1$selection$cv.folds.effective, 4L)
  expect_equal(fit1$selection$cv.repeats.requested, 3L)
  expect_equal(fit1$selection$cv.repeats.effective, 1L)
  expect_match(fit1$selection$fold.message, "forced to 1")
  expect_equal(nrow(fit1$selection$candidate.table), 4L)
  expect_true(all(is.finite(fit1$fitted.values)))
  expect_true(is.finite(fit1$selection$candidate.table$loss.mean[
    fit1$selection$selected.index
  ]))
})

test_that("fit.malps CV seed makes generated folds deterministic", {
  X <- matrix(seq(0, 1, length.out = 18L), ncol = 1L)
  y <- cos(2 * pi * X[, 1L])

  fit1 <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "adaptive.radius",
    degree.grid = 1L,
    min.support.grid = c(5L, 7L),
    kernel.grid = c("epanechnikov", "gaussian"),
    cv.folds = 3L,
    cv.repeats = 2L,
    cv.seed = 42L
  )
  fit2 <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "adaptive.radius",
    degree.grid = 1L,
    min.support.grid = c(5L, 7L),
    kernel.grid = c("epanechnikov", "gaussian"),
    cv.folds = 3L,
    cv.repeats = 2L,
    cv.seed = 42L
  )

  expect_equal(fit1$fitted.values, fit2$fitted.values, tolerance = 1e-12)
  expect_equal(fit1$selection$foldid.source, "generated_seeded")
  expect_equal(fit1$selection$cv.folds, 3L)
  expect_equal(fit1$selection$cv.folds.requested, 3L)
  expect_equal(fit1$selection$cv.folds.effective, 3L)
  expect_equal(fit1$selection$cv.repeats.effective, 2L)
  expect_equal(fit1$selection$foldid.list, fit2$selection$foldid.list)
  expect_equal(nrow(fit1$selection$candidate.table), 4L)
})

test_that("fit.malps default adaptive-radius CV tunes support", {
  X <- matrix(seq(0, 1, length.out = 30L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])

  fit <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "adaptive.radius",
    cv.folds = 3L,
    cv.seed = 7L
  )

  expect_gt(nrow(fit$selection$candidate.table), 1L)
  expect_gt(length(unique(fit$selection$candidate.table$min.support)), 1L)
  expect_true(all(fit$selection$candidate.table$min.support <= nrow(X)))
})

test_that("fit.malps generated CV warns when unseeded and validates fold count", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]

  expect_warning(
    fit <- fit.malps(
      X, y,
      support.selection = "cv",
      support.type = "adaptive.radius",
      degree.grid = 0L,
      min.support.grid = 2L,
      cv.folds = 2L
    ),
    "not reproducible"
  )
  expect_equal(fit$selection$foldid.source, "generated_unseeded")

  expect_error(
    fit.malps(
      X, y,
      support.selection = "cv",
      support.type = "adaptive.radius",
      cv.folds = nrow(X) + 1L
    ),
    "cv.folds cannot exceed nrow"
  )
})

test_that("fit.malps CV records held-out no-coverage diagnostics", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]^2

  fit <- fit.malps(
    X, y,
    anchor.index = c(1L, 8L),
    support.selection = "cv",
    support.type = "fixed.radius",
    degree.grid = 0L,
    radius.grid = c(0.2, 1.1),
    min.support.grid = 1L,
    kernel.grid = "epanechnikov",
    foldid = rep(1:2, length.out = nrow(X))
  )

  candidate <- fit$selection$candidate.table
  small <- candidate[candidate$radius == 0.2, , drop = FALSE]
  large <- candidate[candidate$radius == 1.1, , drop = FALSE]
  expect_gt(small$n.no.coverage, 0L)
  expect_gt(small$n.fit.failures, 0L)
  expect_false(is.finite(small$loss.mean))
  expect_equal(large$n.no.coverage, 0L)
  expect_true(is.finite(large$loss.mean))
  expect_match(
    paste(fit$selection$fold.table$message, collapse = "\n"),
    "no positive averaging coverage"
  )
})

test_that("fit.malps CV fails explicitly when no candidate is valid", {
  X <- matrix(seq(0, 1, length.out = 8L), ncol = 1L)
  y <- X[, 1L]

  expect_error(
    fit.malps(
      X, y,
      support.selection = "cv",
      support.type = "fixed.radius",
      degree.grid = 1L,
      radius.grid = 1e-6,
      min.support.grid = 3L,
      foldid = rep(1:2, length.out = nrow(X))
    ),
    "All MALPS CV candidates failed"
  )
})

test_that("fit.malps GCV selects the minimum exact dense GCV candidate", {
  X <- matrix(seq(0, 1, length.out = 18L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L]) + 0.15 * X[, 1L]

  fit <- fit.malps(
    X, y,
    support.selection = "gcv",
    support.type = "adaptive.radius",
    degree.grid = 0:1,
    min.support.grid = c(4L, 7L),
    kernel.grid = c("epanechnikov", "gaussian")
  )
  candidates <- fit$selection$candidate.table
  selected <- fit$selection$selected.index

  expect_s3_class(fit, "malps")
  expect_equal(fit$selection$method, "gcv")
  expect_equal(fit$selection$loss.name, "gcv")
  expect_equal(nrow(candidates), 8L)
  expect_true(all(candidates$valid))
  expect_equal(selected, candidates$candidate.id[which.min(candidates$gcv)])
  expect_equal(fit$degree, candidates$degree[selected])
  expect_equal(fit$kernel, candidates$kernel[selected])
  expect_equal(fit$diagnostics$min.support, candidates$min.support[selected])
  expect_true(is.finite(candidates$edf[selected]))
  expect_true(is.finite(candidates$mse[selected]))
})

test_that("fit.malps GCV validates exact dense linear-smoother contracts", {
  X <- matrix(seq(0, 1, length.out = 10L), ncol = 1L)
  y <- cos(2 * pi * X[, 1L])

  expect_error(
    fit.malps(
      X, y,
      support.selection = "gcv",
      degree.grid = 1L,
      min.support.grid = 4L,
      robust.iterations = 1L
    ),
    "requires robust.iterations = 0"
  )
  expect_error(
    fit.malps(
      X, y,
      support.selection = "gcv",
      degree.grid = 1L,
      min.support.grid = 4L,
      gcv.exact.max.n = 5L
    ),
    "exceeds gcv.exact.max.n"
  )
})

test_that("malps.smoother.matrix reproduces fixed-support fits and refits", {
  X <- matrix(seq(0, 1, length.out = 14L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])

  fit <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "knn",
    support.size = 8L,
    min.support = 4L,
    kernel = "gaussian"
  )
  S <- malps.smoother.matrix(fit)

  expect_equal(dim(S), c(nrow(X), nrow(X)))
  expect_equal(as.numeric(S %*% y), fit$fitted.values, tolerance = 1e-10)
  expect_false(isTRUE(attr(S, "conditional.on.selection")))
  expect_false(isTRUE(attr(S, "robust.fixed.weight.linearization")))

  y2 <- cos(2 * pi * X[, 1L])
  refit <- refit.malps(fit, y = y2)
  expect_equal(as.numeric(S %*% y2), refit$fitted.values, tolerance = 1e-10)
})

test_that("malps.smoother.matrix handles weighted refits conditionally", {
  X <- matrix(seq(0, 1, length.out = 12L), ncol = 1L)
  y <- 1 + X[, 1L] + X[, 1L]^2
  fit <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "knn",
    support.size = 9L,
    min.support = 4L,
    kernel = "gaussian"
  )
  weights <- rep(1, nrow(X))
  weights[c(3L, 9L)] <- 0
  y2 <- sin(X[, 1L])
  refit <- refit.malps(fit, y = y2, weights = weights)
  S <- malps.smoother.matrix(refit)

  expect_equal(as.numeric(S %*% y2), refit$fitted.values, tolerance = 1e-10)
  expect_false(any(abs(S[, c(3L, 9L)]) > 0))
})

test_that("malps.gcv reports EDF, GCV, and analytic LOOCV diagnostics", {
  X <- matrix(seq(0, 1, length.out = 16L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 7L,
    min.support = 3L,
    kernel = "gaussian"
  )
  S <- malps.smoother.matrix(fit)
  gcv <- malps.gcv(fit, smoother.matrix = S)
  residuals <- y - as.numeric(S %*% y)
  edf <- sum(base::diag(S))
  expected.gcv <- mean(residuals^2) / (1 - edf / length(y))^2

  expect_s3_class(gcv, "malps_gcv")
  expect_equal(gcv$fitted.values, fit$fitted.values, tolerance = 1e-10)
  expect_equal(gcv$edf, edf, tolerance = 1e-12)
  expect_equal(gcv$gcv, expected.gcv, tolerance = 1e-12)
  expect_equal(gcv$loocv.residuals,
               residuals / (1 - base::diag(S)),
               tolerance = 1e-12)
  expect_true(is.finite(gcv$loocv.mse))
})

test_that("malps.smoother.matrix validates robust and size contracts", {
  X <- matrix(seq(0, 1, length.out = 10L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 6L,
    min.support = 2L,
    kernel = "gaussian",
    robust.iterations = 1L
  )

  expect_error(malps.smoother.matrix(fit), "response-dependent")
  S <- malps.smoother.matrix(fit, allow.robust = TRUE)
  expect_equal(as.numeric(S %*% y), fit$fitted.values, tolerance = 1e-10)
  expect_true(isTRUE(attr(S, "robust.fixed.weight.linearization")))
  expect_error(malps.smoother.matrix(fit, max.n = 5L), "exceeds max.n")
  expect_error(malps.gcv(fit), "response-dependent")
})

test_that("fit.malps robust iterations preserve exact local polynomial fits", {
  X <- matrix(seq(0, 1, length.out = 12L), ncol = 1L)
  y <- 1 + 2 * X[, 1L] - 0.5 * X[, 1L]^2

  plain <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "knn",
    support.size = 8L,
    min.support = 6L,
    kernel = "gaussian",
    robust.iterations = 0L
  )
  robust <- fit.malps(
    X, y,
    degree = 2L,
    support.type = "knn",
    support.size = 8L,
    min.support = 6L,
    kernel = "gaussian",
    robust.iterations = 2L
  )

  expect_equal(robust$fitted.values, plain$fitted.values, tolerance = 1e-10)
  expect_equal(robust$diagnostics$robust.iterations, 2L)
  expect_true(robust$diagnostics$robust.used)
  expect_true(all(unlist(robust$robust.weights) == 1))
})

test_that("fit.malps robust iterations downweight local outliers", {
  X <- matrix(seq(0, 1, length.out = 21L), ncol = 1L)
  truth <- 1 + X[, 1L]
  y <- truth
  y[11L] <- y[11L] + 20

  plain <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 21L,
    min.support = 2L,
    kernel = "gaussian",
    robust.iterations = 0L
  )
  robust <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 21L,
    min.support = 2L,
    kernel = "gaussian",
    robust.iterations = 4L
  )

  plain.err <- sqrt(mean((plain$fitted.values - truth)^2))
  robust.err <- sqrt(mean((robust$fitted.values - truth)^2))
  expect_lt(robust.err, plain.err)
  expect_true(any(unlist(robust$robust.weights) < 0.5))
  expect_true(any(robust$diagnostics$robust.downweighted.count > 0L))
  expect_true(all(robust$diagnostics$robust.weight.min >= 0))
  expect_true(all(robust$diagnostics$robust.weight.max <= 1))
})

test_that("robust fitting does not change model-quality averaging weights", {
  X <- matrix(seq(0, 1, length.out = 21L), ncol = 1L)
  truth <- 1 + X[, 1L]
  y <- truth
  y[11L] <- y[11L] + 20

  for (rule in c("condition", "quality")) {
    plain <- fit.malps(
      X, y,
      degree = 1L,
      support.type = "knn",
      support.size = 15L,
      min.support = 2L,
      kernel = "gaussian",
      model.weight.rule = rule,
      robust.iterations = 0L
    )
    robust <- fit.malps(
      X, y,
      degree = 1L,
      support.type = "knn",
      support.size = 15L,
      min.support = 2L,
      kernel = "gaussian",
      model.weight.rule = rule,
      robust.iterations = 4L
    )

    expect_equal(robust$model.weights, plain$model.weights, tolerance = 1e-12)
    expect_equal(robust$averaging.weights, plain$averaging.weights,
                 tolerance = 1e-12)
    expect_true(any(unlist(robust$robust.weights) < 0.5))
    expect_equal(robust$diagnostics$model.condition.number,
                 plain$diagnostics$model.condition.number,
                 tolerance = 1e-12)
    expect_equal(robust$diagnostics$model.weight.raw,
                 plain$diagnostics$model.weight.raw,
                 tolerance = 1e-12)
  }
})

test_that("refit.malps preserves robust controls and recomputes robust weights", {
  X <- matrix(seq(0, 1, length.out = 18L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  y[9L] <- y[9L] + 8

  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 12L,
    min.support = 2L,
    kernel = "gaussian",
    robust.iterations = 2L,
    robust.tuning.constant = 5
  )
  refit <- refit.malps(fit, y = y + 0.1)

  expect_equal(refit$robust.iterations, 2L)
  expect_equal(refit$robust.tuning.constant, 5)
  expect_equal(refit$diagnostics$robust.iterations, 2L)
  expect_equal(length(refit$robust.weights), length(fit$supports))
  expect_true(all(is.finite(refit$fitted.values)))
})

test_that("weighted robust refits ignore zero-case-weight rows in robust scale", {
  X <- matrix(seq(0, 1, length.out = 12L), ncol = 1L)
  truth <- 1 + X[, 1L]
  fit <- fit.malps(
    X, truth,
    degree = 1L,
    support.type = "knn",
    support.size = 12L,
    min.support = 2L,
    kernel = "gaussian",
    robust.iterations = 3L
  )

  y <- truth
  y[1:6] <- 1e6
  y[9L] <- y[9L] + 20
  case.weights <- c(rep(0, 6L), rep(1, 6L))
  refit <- refit.malps(fit, y = y, weights = case.weights)

  zero.row.weights <- unlist(Map(function(s, rw) {
    rw[s$index <= 6L]
  }, refit$supports, refit$robust.weights))
  outlier.row.weights <- unlist(Map(function(s, rw) {
    rw[s$index == 9L]
  }, refit$supports, refit$robust.weights))

  expect_equal(zero.row.weights, rep(1, length(zero.row.weights)))
  expect_true(any(outlier.row.weights < 0.5))
  expect_true(all(refit$diagnostics$robust.downweighted.count <= 6L))
  expect_true(all(is.finite(refit$fitted.values)))
})

test_that("fit.malps CV records robust controls on candidates", {
  X <- matrix(seq(0, 1, length.out = 18L), ncol = 1L)
  y <- cos(2 * pi * X[, 1L])

  fit <- fit.malps(
    X, y,
    support.selection = "cv",
    support.type = "adaptive.radius",
    degree.grid = 1L,
    min.support.grid = c(5L, 7L),
    kernel.grid = "gaussian",
    cv.folds = 3L,
    cv.seed = 22L,
    robust.iterations = 1L,
    robust.tuning.constant = 4
  )

  expect_equal(fit$diagnostics$robust.iterations, 1L)
  expect_equal(fit$diagnostics$robust.tuning.constant, 4)
  expect_true(all(fit$selection$candidate.table$robust.iterations == 1L))
  expect_true(all(fit$selection$candidate.table$robust.tuning.constant == 4))
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("bootstrap.malps returns deterministic Bayesian uncertainty summaries", {
  X <- matrix(seq(0, 1, length.out = 14L), ncol = 1L)
  y <- sin(2 * pi * X[, 1L])
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 8L,
    min.support = 3L,
    kernel = "gaussian"
  )

  boot1 <- bootstrap.malps(fit, B = 6L, seed = 101L, keep.weights = TRUE)
  boot2 <- bootstrap.malps(fit, B = 6L, seed = 101L, keep.weights = TRUE)

  expect_s3_class(boot1, "malps_bootstrap")
  expect_equal(boot1$B.completed, 6L)
  expect_equal(boot1$n.failures, 0L)
  expect_equal(dim(boot1$replicate.fitted.values), c(nrow(X), 6L))
  expect_equal(dim(boot1$replicate.weights), c(nrow(X), 6L))
  expect_equal(boot1$weights, boot1$replicate.weights)
  expect_s3_class(boot1$failures, "data.frame")
  expect_named(boot1$failures, c("attempt", "failure", "message"))
  expect_equal(nrow(boot1$failures), 0L)
  expect_equal(boot1$replicate.fitted.values,
               boot2$replicate.fitted.values,
               tolerance = 1e-12)
  expect_equal(boot1$replicate.weights, boot2$replicate.weights,
               tolerance = 1e-12)
  expect_equal(boot1$summary$index, seq_len(nrow(X)))
  expect_equal(boot1$summary$fitted, fit$fitted.values, tolerance = 1e-12)
  expect_true(all(boot1$summary$lower <= boot1$summary$upper))

  refit1 <- refit.malps(fit, weights = boot1$replicate.weights[, 1L])
  expect_equal(boot1$replicate.fitted.values[, 1L],
               refit1$fitted.values,
               tolerance = 1e-12)
})

test_that("bootstrap.malps supports multinomial weights and validates inputs", {
  X <- matrix(seq(0, 1, length.out = 10L), ncol = 1L)
  y <- 1 + X[, 1L]
  fit <- fit.malps(
    X, y,
    degree = 1L,
    support.type = "knn",
    support.size = 10L,
    min.support = 2L,
    kernel = "gaussian"
  )

  boot <- bootstrap.malps(
    fit,
    B = 4L,
    weight.type = "multinomial",
    seed = 22L,
    keep.weights = TRUE
  )
  expect_equal(boot$B.completed, 4L)
  expect_true(all(colSums(boot$replicate.weights) == nrow(X)))
  expect_true(all(is.finite(boot$replicate.fitted.values)))

  expect_error(bootstrap.malps(fit, B = 0L), "B must be a positive integer")
  expect_error(bootstrap.malps(fit, conf.level = 1), "conf.level")
  expect_error(
    bootstrap.malps(fit, weight.type = "bayesian", probs = rep(1, nrow(X))),
    "probs"
  )
  expect_error(
    bootstrap.malps(fit, weight.type = "multinomial",
                    probs = c(rep(1, nrow(X) - 1L), NA_real_)),
    "probs"
  )
})

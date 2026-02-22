test_that("safe spline predict handles single-point input", {
  res <- gflow:::.gflow.safe.spline.predict(
    x = 2,
    y = 5,
    xout = c(1, 2, 3)
  )

  expect_identical(res$method, "constant")
  expect_equal(res$yhat.in, 5)
  expect_equal(res$yhat.out, c(5, 5, 5))
})

test_that("prepare spline data collapses duplicated x values with weights", {
  dat <- gflow:::.gflow.prepare.spline.data(
    x = c(1, 1, 2, 3, 3),
    y = c(2, 4, 6, 10, 14),
    w = c(1, 3, 2, 1, 1)
  )

  expect_equal(dat$x, c(1, 2, 3))
  expect_equal(unname(dat$y[1]), (2 * 1 + 4 * 3) / 4)
  expect_equal(unname(dat$y[2]), 6)
  expect_equal(unname(dat$y[3]), 12)
  expect_equal(unname(dat$w), c(4, 2, 2))
})

test_that("gcv smoother returns expected structure and dimensions", {
  set.seed(11)
  x <- seq(0, 1, length.out = 40)
  y <- sin(2 * pi * x) + rnorm(length(x), sd = 0.05)

  fit <- gflow:::.gflow.fit.gcv.smoother(
    x = x,
    y = y,
    params = list(grid.size = 75)
  )

  expect_s3_class(fit, "gflow_spline_fit")
  expect_equal(length(fit$xgrid), 75)
  expect_equal(length(fit$gpredictions), 75)
  expect_equal(dim(fit$gpredictions.CrI), c(2, 75))
  expect_true(all(is.finite(fit$gpredictions)))
})

test_that("spline wrapper supports robust repeated-CV selection with 1SE rule", {
  set.seed(44)
  x <- sort(runif(80, -1, 1))
  y <- sin(3 * x) + rnorm(length(x), sd = 0.12)

  fit <- gflow::gflow.smooth.spline(
    x = x,
    y = y,
    use.gcv = FALSE,
    cv.folds = 4,
    cv.repeats = 2,
    cv.one.se = TRUE,
    cv.seed = 123
  )

  expect_true(!is.null(fit))
  expect_true(is.list(fit$gflow.selection))
  expect_true(grepl("^cv_", fit$gflow.selection$method))
  expect_true(is.finite(fit$gflow.selection$selected.spar))
  expect_true(nrow(fit$gflow.selection$cv.table) > 0)
})

test_that("spline wrapper enforces df.max cap", {
  set.seed(45)
  x <- sort(runif(100))
  y <- sin(8 * x) + rnorm(length(x), sd = 0.08)

  fit <- gflow::gflow.smooth.spline(
    x = x,
    y = y,
    use.gcv = FALSE,
    df.max = 6,
    cv.folds = 5,
    cv.repeats = 2,
    cv.seed = 321
  )

  expect_true(!is.null(fit))
  expect_true(is.finite(fit$df))
  expect_true(fit$df <= 6 + 1e-8)
})

test_that("external BB spline fit returns expected prediction matrix shape", {
  set.seed(22)
  x <- seq(0, 1, length.out = 30)
  y <- cos(2 * pi * x) + rnorm(length(x), sd = 0.02)

  lambda <- replicate(
    6,
    {
      w <- rexp(length(x))
      w / sum(w) * length(x)
    }
  )

  fit <- gflow:::.gflow.with.external.BB.spline(
    x = x,
    y = y,
    lambda = lambda,
    grid.size = 64
  )

  expect_equal(length(fit$xgrid), 64)
  expect_equal(length(fit$gpredictions), 64)
  expect_equal(dim(fit$BB.gpredictions), c(64, 6))
  expect_true(all(is.finite(fit$BB.gpredictions)))
})

test_that("external BB spline reuses a single selected spar across replicates", {
  set.seed(46)
  x <- sort(runif(70, -1, 1))
  y <- cos(4 * x) + rnorm(length(x), sd = 0.1)

  lambda <- replicate(
    5,
    {
      w <- rexp(length(x))
      w / sum(w) * length(x)
    }
  )

  fit <- gflow:::.gflow.with.external.BB.spline(
    x = x,
    y = y,
    lambda = lambda,
    use.gcv = FALSE,
    cv.folds = 4,
    cv.repeats = 2,
    cv.seed = 456
  )

  expect_true(is.finite(fit$opt.bw))
  expect_equal(dim(fit$BB.gpredictions), c(400, 5))
  expect_true(grepl("smooth.spline", fit$fit.method))
})

test_that("iknn trend helper returns a breakpoint for adequate input", {
  x <- 1:12
  y <- c(6, 5, 4, 3, 2.5, 2.2, 2.1, 2.2, 2.6, 3.3, 4.1, 5.2)

  trend <- gflow:::.gflow.fit.iknn.trend(x, y)

  expect_s3_class(trend, "gflow_trend_fit")
  expect_true(is.finite(trend$breakpoint))
  expect_true(trend$breakpoint >= min(x))
  expect_true(trend$breakpoint <= max(x))
})

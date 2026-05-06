test_that("isometry helpers identify scaled exact distances", {
  D.true <- as.matrix(stats::dist(matrix(1:5, ncol = 1)))
  D.estimated <- 2 * D.true

  expect_equal(isometry.scale(D.estimated, D.true), 0.5, tolerance = 1e-12)
  expect_equal(isometry.rel.rms.error(D.estimated, D.true), 0, tolerance = 1e-12)
  expect_equal(isometry.rel.abs.error(D.estimated, D.true), c(`50%` = 0, `95%` = 0),
               tolerance = 1e-12)
  expect_equal(
    isometry.distortion.quantiles(D.estimated, D.true),
    c(`5%` = 1, `50%` = 1, `95%` = 1),
    tolerance = 1e-12
  )
  expect_equal(isometry.distance.correlations(D.estimated, D.true),
               c(pearson_cor = 1, spearman_cor = 1),
               tolerance = 1e-12)
})


test_that("summarize.isometry.deviation returns benchmark columns", {
  D.true <- as.matrix(stats::dist(matrix(1:4, ncol = 1)))
  D.estimated <- D.true
  D.estimated[1, 4] <- D.estimated[4, 1] <- 4

  out <- summarize.isometry.deviation(D.estimated, D.true, scale = FALSE)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_true(all(c(
    "scale", "rel_rms_error", "rel_abs_error_median", "rel_abs_error_q95",
    "distortion_q05", "distortion_median", "distortion_q95",
    "pearson_cor", "spearman_cor"
  ) %in% names(out)))
  expect_equal(out$scale, 1)
  expect_gt(out$rel_rms_error, 0)
  expect_gt(out$rel_abs_error_q95, 0)
})


test_that("isometry helpers validate distance matrices", {
  D <- diag(3)

  expect_error(isometry.scale(D[1:2, ], D), "square")
  expect_error(isometry.scale(D, D[1:2, 1:2]), "same dimensions")
  D.bad <- D
  D.bad[1, 2] <- NA
  expect_error(isometry.scale(D.bad, D), "cannot contain")

  D.true <- matrix(1, 3, 3)
  diag(D.true) <- 0
  D.zero <- matrix(0, 3, 3)
  expect_error(isometry.scale(D.zero, D.true), "zero norm")
})

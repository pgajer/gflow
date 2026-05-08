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
  expect_equal(
    isometry.geodesic.diagnostics(D.estimated, D.true),
    c(
      rel_geodesic_stress = 0,
      signed_bias = 0,
      shortcut_fraction = 0,
      q50_rel_abs_residual = 0,
      q90_rel_abs_residual = 0,
      q95_rel_abs_residual = 0,
      short_band_bias = 0,
      mid_band_bias = 0,
      long_band_bias = 0
    ),
    tolerance = 1e-12
  )
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
  expect_equal(out$rel_geodesic_stress, out$rel_rms_error)
  expect_true(all(c(
    "signed_bias", "shortcut_fraction",
    "q50_rel_abs_residual", "q90_rel_abs_residual", "q95_rel_abs_residual",
    "short_band_bias", "mid_band_bias", "long_band_bias"
  ) %in% names(out)))
})


test_that("geodesic diagnostics report signed shortcut and detour behavior", {
  D.true <- as.matrix(stats::dist(matrix(1:4, ncol = 1)))
  D.short <- D.true
  D.short[1, 4] <- D.short[4, 1] <- 1

  short <- isometry.geodesic.diagnostics(D.short, D.true, scale = FALSE)

  expect_gt(short[["rel_geodesic_stress"]], 0)
  expect_lt(short[["signed_bias"]], 0)
  expect_gt(short[["shortcut_fraction"]], 0)
  expect_gt(short[["q95_rel_abs_residual"]], 0)

  D.long <- D.true
  D.long[1, 4] <- D.long[4, 1] <- 5
  detour <- isometry.geodesic.diagnostics(D.long, D.true, scale = FALSE)

  expect_gt(detour[["signed_bias"]], 0)
  expect_lt(detour[["shortcut_fraction"]], 1)
  expect_gt(detour[["long_band_bias"]], 0)
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
  expect_error(isometry.geodesic.diagnostics(D, D, band.probs = c(0.8, 0.2)),
               "band.probs")
})

test_that("geosmooth public APIs are not exposed or defined by gflow", {
  migrated <- c(
    "bootstrap.malps",
    "fit.graph.trend.filtering",
    "fit.lpl.tf",
    "fit.lps",
    "fit.malps",
    "fit.metric.graph.lowpass",
    "fit.ps.lps",
    "fit.pttf.trend.filtering",
    "fit.slpl.tf",
    "fit.ssrhe.hessian.l1.regression",
    "fit.ssrhe.hessian.regression",
    "fit.ssrhe.hessian.regression.cv",
    "fit.ssrhe.hessian.regression.gcv",
    "get.region.boundary",
    "graph.low.pass.filter",
    "graph.trend.filtering.operator",
    "harmonic.smoother",
    "lpl.tf.operator",
    "lps.backend.diagnostics",
    "malps.gcv",
    "malps.smoother.matrix",
    "metric.graph.lowpass.operator",
    "perform.harmonic.smoothing",
    "pttf.geometry",
    "pttf.operator",
    "pttf.operator.filter.rows",
    "refit.lpl.tf",
    "refit.malps",
    "refit.metric.graph.lowpass",
    "refit.slpl.tf",
    "refit.ssrhe.hessian.l1.regression",
    "refit.ssrhe.hessian.regression",
    "slpl.tf.operator",
    "ssrhe.hessian.operator",
    "ssrhe.support.grid",
    "transported.graph.hessian.operator"
  )

  expect_false(any(migrated %in% getNamespaceExports("gflow")))

  still.defined <- vapply(
    migrated,
    exists,
    logical(1),
    envir = asNamespace("gflow"),
    inherits = FALSE
  )
  expect_false(any(still.defined))
})

test_that("geosmooth local PCA chart-dimension helpers are not defined by gflow", {
  migrated.helpers <- c(
    ".local.pca.auto.chart.dim",
    ".local.pca.auto.chart.dim.result",
    ".local.pca.auto.chart.dim.row",
    ".local.pca.auto.chart.dim.with.metric",
    ".local.pca.max.chart.dim.for.support"
  )

  still.defined <- vapply(
    migrated.helpers,
    exists,
    logical(1),
    envir = asNamespace("gflow"),
    inherits = FALSE
  )
  expect_false(any(still.defined))
})

test_that("catch-all support files expose only owned widget workflows", {
  source.files <- c(
    "stats_utils.R",
    "plot_utils.R",
    "synthetic_data_utils.R",
    "grids.R",
    "random_sampling.R",
    "preprocess_matrix.R",
    "hist_utils.R",
    "boxcox_utils.R",
    "divergences.R",
    "wasserstein_dist.R"
  )
  source.files <- testthat::test_path("..", "..", "R", source.files)

  top.level.functions <- function(path) {
    expressions <- parse(path)
    out <- character()
    for (expr in expressions) {
      is.assignment <- is.call(expr) &&
        is.symbol(expr[[1L]]) &&
        as.character(expr[[1L]]) %in% c("<-", "=")
      if (is.assignment &&
          is.symbol(expr[[2L]]) &&
          is.call(expr[[3L]]) &&
          identical(expr[[3L]][[1L]], as.name("function"))) {
        out <- c(out, as.character(expr[[2L]]))
      }
    }
    out
  }

  defined <- unique(unlist(lapply(source.files, top.level.functions)))
  public <- intersect(defined, getNamespaceExports("gflow"))
  expect_setequal(
    public,
    c(
      "plot3D.plain.widget",
      "plot3D.cont.widget",
      "plot3D.cltrs.widget"
    )
  )
})

test_that("matrix preprocessing helpers have one top-level definition", {
  source.files <- list.files(
    testthat::test_path("..", "..", "R"),
    pattern = "\\.R$",
    full.names = TRUE
  )
  assignments <- unlist(lapply(source.files, function(path) {
    expressions <- parse(path)
    vapply(expressions, function(expr) {
      if (is.call(expr) &&
          is.symbol(expr[[1L]]) &&
          as.character(expr[[1L]]) %in% c("<-", "=") &&
          is.symbol(expr[[2L]]) &&
          is.call(expr[[3L]]) &&
          identical(expr[[3L]][[1L]], as.name("function"))) {
        return(as.character(expr[[2L]]))
      }
      NA_character_
    }, character(1L))
  }))
  assignments <- assignments[!is.na(assignments)]

  expect_identical(sum(assignments == "winsorize.vec"), 1L)
  expect_identical(sum(assignments == "robust.zscore.vec"), 1L)
})

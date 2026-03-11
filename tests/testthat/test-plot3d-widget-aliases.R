test_that("plot3D widget helpers return htmlwidgets and html aliases remain available", {
  skip_if_not_installed("rgl")
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("htmltools")

  X <- matrix(
    c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ),
    ncol = 3,
    byrow = TRUE
  )

  w_plain <- plot3D.plain.widget(X, size = 2)
  expect_s3_class(w_plain, "htmlwidget")
  expect_false(is.null(attr(w_plain, "ids")))

  w_plain_alias <- plot3D.plain.html(X, size = 2)
  expect_s3_class(w_plain_alias, "htmlwidget")

  y <- c(0.1, 0.3, 0.7, 0.9)
  w_cont <- plot3D.cont.widget(X, y = y, legend.show = FALSE)
  expect_s3_class(w_cont, "htmlwidget")
  expect_false(is.null(attr(w_cont, "y.col.tbl")))

  w_cont_alias <- plot3D.cont.html(X, y = y, legend.show = FALSE)
  expect_s3_class(w_cont_alias, "htmlwidget")

  cltr <- c("A", "A", "B", "B")
  w_cltr <- plot3D.cltrs.widget(X, cltr = cltr, show.legend = FALSE, show.cltr.labels = FALSE)
  expect_s3_class(w_cltr, "htmlwidget")
  expect_false(is.null(attr(w_cltr, "cltr.col.tbl")))

  w_cltr_alias <- plot3D.cltrs.html(X, cltr = cltr, show.legend = FALSE, show.cltr.labels = FALSE)
  expect_s3_class(w_cltr_alias, "htmlwidget")
})

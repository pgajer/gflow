test_that("generic browser implementations have moved out of gflow", {
    retired <- c("plot3D.plain.widget", "plot3D.cont.widget", "plot3D.cltrs.widget",
                 "plot3D.plain.html", "plot3D.cont.html", "plot3D.cltrs.html",
                 ".run_plot3d_html_layers", "quantize.for.legend",
                 "plot3D.plain", "plot3D.cont", "plot3D.cltrs",
                 "plot3D.tree", "plot3D.path")
    ns <- asNamespace("gflow")
    expect_false(any(vapply(retired, exists, logical(1), envir = ns, inherits = FALSE)))
    expect_true(exists("quantize.cont.var", envir = ns, inherits = FALSE))
    expect_true(exists(".compute.igraph.layout", envir = ns, inherits = FALSE))
})

test_that("gflow and ivue coexist without changing analysis or native helpers", {
    skip_if_not_installed("ivue")
    skip_if_not_installed("rgl")
    X <- matrix(seq_len(30), ncol = 3)
    w <- ivue::plot3D.cont(X, values = seq_len(nrow(X)))
    expect_s3_class(w, "htmlwidget")
    expect_equal(attr(w, "ivue")$X, X)
    expect_equal(length(rgl::rgl.dev.list()), 0L)
    expect_true(exists("plot3D.graph", envir = asNamespace("gflow"), inherits = FALSE))
})

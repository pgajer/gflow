test_that("retired generic names are not called by remaining adapters", {
    for (name in c("show.3d.cl", "plot3D.cl", "map.S.to.X", "plot3D.diskEmbdg")) {
        fun <- get(name, envir = asNamespace("gflow"))
        expect_false(grepl(
            "(?<![[:alnum:]_.:])plot3D[.](plain|cont|cltrs|tree|path)\\s*\\(",
            paste(deparse(body(fun)), collapse = "\n"), perl = TRUE))
    }
})

test_that("cluster and selection adapters return complete ivue scenes", {
    skip_if_not_installed("ivue")
    skip_if_not_installed("rgl")
    .require.ivue.plotting()
    old <- options(rgl.useNULL = TRUE)
    on.exit(options(old), add = TRUE)
    X <- matrix(seq_len(18) / 10, ncol = 3,
                dimnames = list(letters[1:6], c("x", "y", "z")))
    groups <- c("A", "A", "B", "B", "C", "C")
    before <- rgl::rgl.dev.list()
    w <- show.3d.cl("B", groups, X, cl.radius = 0.08, show.labels = TRUE,
                    show.ref.cltr = FALSE, show.cltr.labels = FALSE)
    info <- attr(w, "ivue")
    expect_s3_class(w, "htmlwidget")
    expect_equal(info$X, X)
    expect_equal(which(info$highlight), c(3L, 4L))
    labels <- Filter(function(x) x$type == "text" && any(nzchar(x$texts)), info$scene$objects)
    expect_setequal(unlist(lapply(labels, `[[`, "texts")), c("c", "d"))
    spheres <- Filter(function(x) x$type == "spheres", info$scene$objects)
    expect_true(any(vapply(spheres, function(x) all(abs(x$radii - 0.08) < 1e-7), logical(1))))
    points <- Filter(function(x) x$type == "points", info$scene$objects)
    expect_true(all(vapply(points, function(x) all(x$colors[, 4] == 0), logical(1))))

    w <- plot3D.cl(c("A", "C"), groups, as.data.frame(X), legend.show = FALSE)
    expect_s3_class(w, "htmlwidget")
    labels <- Filter(function(x) x$type == "text" && any(nzchar(x$texts)), attr(w, "ivue")$scene$objects)
    expect_setequal(unlist(lapply(labels, `[[`, "texts")), c("A", "B", "C"))
    for (label in labels) {
        expected <- apply(X[groups == as.character(label$texts), , drop = FALSE], 2, median)
        expect_equal(as.numeric(label$vertices[1, ]), unname(expected), tolerance = 1e-6)
    }
    expect_s3_class(plot3D.cl("A", rep("A", 6), X), "htmlwidget")
    expect_s3_class(plot3D.cl("A", factor(groups, levels = c("A", "B", "C", "unused")), X), "htmlwidget")
    expect_s3_class(show.3d.cl(factor("A"), groups, X), "htmlwidget")
    expect_error(plot3D.cl("absent", groups, X), "existing clusters")

    w <- map.S.to.X(c("b", "e", "absent"), X, radius = 0.09)
    expect_s3_class(w, "htmlwidget")
    expect_equal(which(attr(w, "ivue")$highlight), c(2L, 5L))
    expect_warning(w <- map.S.to.X("absent", X), "No samples")
    expect_false(any(attr(w, "ivue")$highlight))
    expect_identical(rgl::rgl.dev.list(), before)
})

test_that("disk decorations are captured without changing the caller device", {
    skip_if_not_installed("ivue")
    skip_if_not_installed("rgl")
    .require.ivue.plotting()
    old <- options(rgl.useNULL = TRUE)
    on.exit(options(old), add = TRUE)
    device <- rgl::open3d(useNULL = TRUE, silent = TRUE)
    on.exit(rgl::close3d(device), add = TRUE)
    rgl::points3d(0, 0, 0)
    before <- rgl::scene3d(minimal = FALSE)
    object <- list(axis.pos = diag(3), X.ebdg = matrix(seq_len(18) / 20, ncol = 3),
                   X = matrix(1, 6, 3, dimnames = list(NULL, c("a", "b", "c"))), n = 3L)
    w <- plot3D.diskEmbdg(object, edge.col = "magenta", adj.df = matrix(0.5, 3, 2))
    expect_s3_class(w, "htmlwidget")
    objects <- attr(w, "ivue")$scene$objects
    edges <- Filter(function(x) x$type == "lines" &&
        isTRUE(all.equal(as.numeric(x$colors[1, 1:3]), c(1, 0, 1))), objects)
    expect_equal(length(edges), 3L)
    labels <- Filter(function(x) x$type == "text" && any(nzchar(x$texts)), objects)
    expect_setequal(unlist(lapply(labels, `[[`, "texts")), c("a", "b", "c"))
    expect_equal(rgl::cur3d(), device)
    expect_equal(rgl::scene3d(minimal = FALSE), before)
    expect_error(show.3d.cl("A", rep("A", 6), object$X.ebdg,
        layers = list(ivue::layer3D.callback(function(ctx) stop("layer failure")))), "layer failure")
    expect_equal(rgl::scene3d(minimal = FALSE), before)
})

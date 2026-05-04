test_that("linf.simplex kNN uses unfolded distances across faces", {
    X <- rbind(
        c(1.0, 0.2, 0.4),
        c(1.0, 0.7, 0.4),
        c(0.3, 1.0, 0.4),
        c(1.0, 0.2, 0.9)
    )

    res <- .Call("S_linf_simplex_knn",
                 X,
                 as.integer(4L),
                 as.double(sqrt(.Machine$double.eps)),
                 PACKAGE = "gflow")

    expect_equal(res$indices[1, ], c(0L, 1L, 3L, 2L))
    expect_equal(res$distances[1, ], c(0, 0.5, 0.5, 1.5), tolerance = 1e-12)
})

test_that("linf.simplex kNN matches Euclidean kNN within one face", {
    X <- rbind(
        c(1.0, 0.1, 0.1),
        c(1.0, 0.2, 0.1),
        c(1.0, 0.4, 0.1),
        c(1.0, 0.8, 0.1)
    )

    linf <- .Call("S_linf_simplex_knn",
                  X,
                  as.integer(4L),
                  as.double(sqrt(.Machine$double.eps)),
                  PACKAGE = "gflow")
    euclidean <- .Call("S_kNN", X, as.integer(4L), PACKAGE = "gflow")

    expect_equal(linf$indices, euclidean$indices)
    expect_equal(linf$distances, euclidean$distances, tolerance = 1e-12)
})

test_that("linf.simplex public API validates normalized simplex input", {
    X <- rbind(
        c(1.0, 0.2, 0.4),
        c(0.3, 1.0, 0.4),
        c(1.0, 0.2, 0.9)
    )

    graph <- create.single.iknn.graph(
        X,
        k = 1,
        max.path.edge.ratio.deviation.thld = 0,
        threshold.percentile = 0,
        compute.full = FALSE,
        pca.dim = NULL,
        knn.metric = "linf.simplex",
        verbose = FALSE
    )
    expect_equal(attr(graph, "knn.metric"), "linf.simplex")

    expect_error(
        create.single.iknn.graph(
            X,
            k = 1,
            pca.dim = 2,
            knn.metric = "linf.simplex",
            verbose = FALSE
        ),
        "pca.dim must be NULL"
    )

    X.bad <- X
    X.bad[1, ] <- c(0.8, 0.2, 0.4)
    expect_error(
        create.single.iknn.graph(
            X.bad,
            k = 1,
            pca.dim = NULL,
            knn.metric = "linf.simplex",
            verbose = FALSE
        ),
        "L-infinity normalized"
    )
})

test_that("kNN cache distinguishes euclidean and linf.simplex metrics", {
    X <- rbind(
        c(1.0, 0.2, 0.4),
        c(1.0, 0.7, 0.4),
        c(0.3, 1.0, 0.4),
        c(1.0, 0.2, 0.9)
    )
    cache.path <- tempfile(fileext = ".bin")

    create.single.iknn.graph(
        X,
        k = 1,
        max.path.edge.ratio.deviation.thld = 0,
        threshold.percentile = 0,
        compute.full = FALSE,
        pca.dim = NULL,
        knn.metric = "euclidean",
        knn.cache.path = cache.path,
        knn.cache.mode = "write",
        verbose = FALSE
    )

    expect_error(
        create.single.iknn.graph(
            X,
            k = 1,
            max.path.edge.ratio.deviation.thld = 0,
            threshold.percentile = 0,
            compute.full = FALSE,
            pca.dim = NULL,
            knn.metric = "linf.simplex",
            knn.cache.path = cache.path,
            knn.cache.mode = "read",
            verbose = FALSE
        ),
        "metric mismatch"
    )
})

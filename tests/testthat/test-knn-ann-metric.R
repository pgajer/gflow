test_that("low-level ANN kNN uses Euclidean distance", {
    X <- rbind(
        c(0, 0),
        c(3, 4),
        c(10, 0)
    )

    knn <- .Call("S_kNN", X, as.integer(2L), PACKAGE = "gflow")

    expect_equal(knn$indices[1, ], c(0L, 1L))
    expect_equal(knn$distances[1, ], c(0, 5), tolerance = 1e-12)
    expect_equal(knn$indices[3, ], c(2L, 1L))
    expect_equal(knn$distances[3, ], c(0, sqrt(65)), tolerance = 1e-12)
})

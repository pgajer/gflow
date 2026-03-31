test_that("adjusted.rand.index matches expected partition behavior", {
    expect_equal(adjusted.rand.index(c(1, 1, 2, 2), c("a", "a", "b", "b")), 1)
    expect_equal(adjusted.rand.index(c(1, 1, 1, 1), c(2, 2, 2, 2)), 1)
    expect_error(adjusted.rand.index(c(1), c(1)), "at least 2 observations")
    expect_error(adjusted.rand.index(c(1, NA), c(1, 2)), "must not contain NA")
})

test_that("congruence.with.labels computes ARI with native helper", {
    sample.ids <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10")
    cluster.labels <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
    labels <- c(
        s1 = "I", s2 = "I", s3 = "I", s4 = "I", s5 = "I",
        s6 = "III", s7 = "III", s8 = "III", s9 = "III", s10 = "III"
    )

    res <- congruence.with.labels(sample.ids, cluster.labels, labels)
    expect_equal(res$n, 10L)
    expect_equal(unname(res$ari), 1)
    expect_equal(dim(res$cst.table), c(2L, 2L))
})

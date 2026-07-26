test_that("every in-scope S3 method has an enforced class contract", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    source(file.path(root, "tools/cleanup_ledger_lib.R"), local = TRUE)
    source(file.path(root, "tools/s3_namespace_lib.R"), local = TRUE)

    expect_identical(s3.validate.namespace(root), character())
})

test_that("all retained validators reject structurally empty objects", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    ownership <- read.csv(
        file.path(root, "split_audit/cleanup/s3-class-ownership.csv"),
        stringsAsFactors = FALSE
    )
    namespace <- asNamespace("gflow")

    for (i in seq_len(nrow(ownership))) {
        validator <- get(ownership$validator[[i]], envir = namespace)
        empty <- structure(list(), class = ownership$class[[i]])
        expect_error(
            validator(empty),
            info = paste("validator for", ownership$class[[i]])
        )
    }
})

test_that("distance quantile class validator accepts constructor output", {
    result <- distance.quantile.bin.analysis(
        x = seq_len(12),
        d = seq_len(12),
        n.bins = 3L,
        min.per.bin = 2L
    )
    validator <- get(
        ".validate.distance.quantile.bins",
        envir = asNamespace("gflow")
    )

    expect_identical(validator(result), result)
})

test_that("association constructors satisfy their retained S3 contracts", {
    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- lapply(adj.list, function(neighbors) rep(1, length(neighbors)))
    y <- c(1, 2, 4)
    z <- c(2, 5, 3)
    Y <- cbind(y1 = y, y2 = c(4, 1, 3))
    Z <- cbind(z1 = z, z2 = c(3, 1, 5))

    results <- list(
        lcor.vector.matrix = lcor(adj.list, weight.list, y, Z),
        lcor.matrix.matrix = lcor(adj.list, weight.list, Y, Z),
        lslope.gradient = lslope(
            adj.list, weight.list, y, z,
            y.diff.type = "difference",
            z.diff.type = "difference",
            instrumented = TRUE
        ),
        lslope.vector.matrix = lslope(
            adj.list, weight.list, y, Z,
            y.diff.type = "difference",
            z.diff.type = "difference"
        ),
        lslope.matrix.matrix = lslope(
            adj.list, weight.list, Y, Z,
            y.diff.type = "difference",
            z.diff.type = "difference"
        ),
        lslope.neighborhood = lslope.neighborhood(
            adj.list, weight.list, y, z,
            y.diff.type = "difference",
            z.diff.type = "difference"
        )
    )

    for (result in results) {
        expect_silent(capture.output(print(result)))
    }

    expect_silent(capture.output(summary(results$lcor.vector.matrix)))
    expect_silent(capture.output(summary(results$lcor.matrix.matrix)))
    expect_silent(capture.output(summary(results$lslope.vector.matrix)))
    expect_silent(capture.output(summary(results$lslope.matrix.matrix)))
})

test_that("analysis constructors satisfy their retained S3 contracts", {
    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- lapply(adj.list, function(neighbors) rep(1, length(neighbors)))

    modulation <- compute.gfc.modulation(
        adj.list,
        weight.list,
        modulation = "density",
        density = c(1, 2, 3),
        verbose = FALSE
    )
    p.summary <- weighted.p.value.summary(
        c(-1, 0, 1),
        mu = 0,
        sigma = 1
    )

    expect_silent(capture.output(print(modulation)))
    expect_silent(capture.output(plot(modulation)))
    expect_silent(capture.output(print(p.summary)))
})

phase6.registered.s3 <- function() {
    registrations <- getNamespaceInfo(asNamespace("gflow"), "S3methods")
    paste(registrations[, 1L], registrations[, 2L], sep = ":")
}

test_that("orphan and relocated S3 registrations are absent", {
    retired <- c(
        "identify:bottlenecks",
        "plot:chain.with.path",
        "plot:comono.diagnostics",
        "plot:component",
        "plot:components.multi",
        "plot:ggraph",
        "plot:kh.matrix",
        "plot:model.errors",
        "plot:prediction.errors",
        "plot:star_object"
    )

    expect_length(intersect(retired, phase6.registered.s3()), 0L)
    expect_true("plot:graphMScx" %in% phase6.registered.s3())
})

test_that("retired example constructors are private and MADAG API is explicit", {
    exports <- getNamespaceExports("gflow")

    expect_false(any(c("generate.star.dataset", "ggraph") %in% exports))
    expect_true("madag.bottlenecks" %in% exports)
    expect_false("identify.bottlenecks" %in% exports)
})

test_that("distance quantile analysis constructs its registered S3 class", {
    result <- distance.quantile.bin.analysis(
        x = seq_len(12),
        d = seq_len(12),
        n.bins = 3L,
        min.per.bin = 2L
    )

    expect_s3_class(result, "distance.quantile.bins")
    expect_s3_class(result$bins, "data.frame")
    expect_true("plot:distance.quantile.bins" %in% phase6.registered.s3())
})

test_that("retired private helpers are not mistaken for S3 methods", {
    namespace <- asNamespace("gflow")
    old.names <- c(
        "coef.assoc0", "coef.assoc1",
        "print.assoc0", "print.assoc1",
        "summary.assoc0", "summary.assoc1",
        "plot.assoc0", "plot.assoc1",
        "plot.gaussian_mixture", "print.gaussian_mixture",
        "summary.gaussian_mixture", "plot.gaussian_mixture_data",
        "plot.gridX", "plot.box.tiling",
        "transform.logit", "compare.representative.methods"
    )

    expect_false(any(vapply(
        old.names,
        exists,
        logical(1L),
        envir = namespace,
        inherits = FALSE
    )))
})

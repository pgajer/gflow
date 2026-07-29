test_that("breaking release exposes only the maintained non-protected API", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(
            root, "split_audit/cleanup/api-ownership.csv"
        )),
        "source-only release ledger is not shipped in the package tarball"
    )
    ledger <- read.csv(
        file.path(root, "split_audit/cleanup/api-ownership.csv"),
        stringsAsFactors = FALSE
    )
    maintained <- c(
        "analyze.harmonic.extensions",
        "apply.harmonic.extension",
        "build.se.tree",
        "compare.harmonic.methods",
        "compute.gfc.modulation",
        "compute.harmonic.extension",
        "compute.tube.lens.corridor",
        "construct.gflow.graph",
        "construct.madag",
        "extremality.summary",
        "gfassoc.deviation",
        "gfassoc.membership",
        "gfassoc.overlap",
        "gfassoc.polarity",
        "gfcor",
        "label.extremality.3d",
        "lcor",
        "lcor.with.posterior",
        "lslope",
        "lslope.neighborhood",
        "madag.bottlenecks",
        "permutation.test.lcor",
        "select.max.density.trajectory"
    )
    actual <- ledger$item_name[
        ledger$item_type == "export" & ledger$ownership != "PROTECTED"
    ]

    expect_setequal(actual, maintained)
    expect_true(all(
        ledger$ownership[ledger$item_type == "export"] %in%
            c("PROTECTED", "CORE-ANALYSIS", "CORE-ASSOCIATION")
    ))
})

test_that("retired Phase 8 families are absent from the namespace", {
    retired <- c(
        "cluster.graph.louvain",
        "congruence.with.labels",
        "detect.graph.endpoints",
        "find.shortest.paths.within.radius",
        "compute.diffusion.pseudotime.sparse",
        "compute.potential.pseudotime.sparse",
        "phate",
        "phate.core",
        "phate.embed",
        "quadform.metric",
        "fit.rho.randomwalk",
        "distance.quantile.bin.analysis",
        "weighted.p.value",
        "weighted.p.value.summary",
        "select3D.points",
        "plot3D.plain.widget",
        "plot3D.cont.widget",
        "plot3D.cltrs.widget"
    )

    expect_length(intersect(retired, getNamespaceExports("gflow")), 0L)
})

test_that("retained in-scope S3 methods have explicit release ownership", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(
            root, "split_audit/cleanup/s3-class-ownership.csv"
        )),
        "source-only release ledger is not shipped in the package tarball"
    )
    ownership <- read.csv(
        file.path(root, "split_audit/cleanup/s3-class-ownership.csv"),
        stringsAsFactors = FALSE
    )
    methods <- c(
        "as.data.frame.harmonic_extension",
        "plot.gfc_modulation",
        "plot.gflow_graph",
        "print.gfc_modulation",
        "print.gfcor",
        "print.gflow_graph",
        "print.harmonic_extension",
        "print.lcor_matrix_matrix_result",
        "print.lcor_vector_matrix_result",
        "print.lcor.posterior",
        "print.lslope_gradient_result",
        "print.lslope_matrix_matrix_result",
        "print.lslope_neighborhood_result",
        "print.lslope_vector_matrix_result",
        "print.madag",
        "print.se_tree",
        "print.summary.gfcor",
        "print.summary.gflow_graph",
        "print.summary.lcor.posterior",
        "summary.gfcor",
        "summary.gflow_graph",
        "summary.lcor_matrix_matrix_result",
        "summary.lcor_vector_matrix_result",
        "summary.lcor.posterior",
        "summary.lslope_matrix_matrix_result",
        "summary.lslope_vector_matrix_result",
        "summary.madag"
    )
    owned <- unlist(strsplit(ownership$methods, ";", fixed = TRUE))

    expect_setequal(owned, methods)
    expect_true(all(
        ownership$ownership %in% c("CORE-ANALYSIS", "CORE-ASSOCIATION")
    ))
})

test_that("public package story and migration guide are present", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(root, "README.md")),
        "source-only package story is not shipped in the package tarball"
    )
    expect_true(all(file.exists(file.path(
        root,
        c(
            "README.md",
            "REFERENCE.md",
            "split_audit/cleanup/breaking-release-migration.md"
        )
    ))))
    expect_identical(
        unname(utils::packageDescription("gflow", fields = "Version")),
        "0.2.0"
    )
})

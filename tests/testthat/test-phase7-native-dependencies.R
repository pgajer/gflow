test_that("dependency and native ownership ledgers cover the live package", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    command <- file.path(root, "tools/check_phase7_ownership.R")
    skip_if_not(
        file.exists(command),
        "source-only ownership audit is not shipped in the package tarball"
    )
    output <- system2(
        file.path(R.home("bin"), "Rscript"),
        command,
        stdout = TRUE,
        stderr = TRUE
    )

    expect_identical(attr(output, "status"), NULL)
    expect_match(
        paste(output, collapse = "\n"),
        "Phase 7 ownership passes"
    )
})

test_that("the optional OpenMP diagnostic has a registered tested caller", {
    info <- .Call("S_gflow_openmp_diag", PACKAGE = "gflow")

    expect_type(info, "list")
    expect_named(
        info,
        c(
            "openmp_compiled", "_OPENMP",
            "omp_get_max_threads", "omp_get_num_procs"
        )
    )
    expect_type(info$openmp_compiled, "logical")
    expect_length(info$openmp_compiled, 1L)
    expect_true(info$omp_get_max_threads >= 1L)
    expect_true(info$omp_get_num_procs >= 1L)
})

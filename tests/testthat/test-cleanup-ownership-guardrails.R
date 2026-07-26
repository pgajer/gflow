test_that("every export and S3 registration has a complete ownership row", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(root, "tools/cleanup_ledger_lib.R")),
        "source-only ownership audit is not shipped in the package tarball"
    )
    source(file.path(root, "tools/cleanup_ledger_lib.R"), local = TRUE)

    errors <- cleanup.validate.ledger(root)
    expect_identical(errors, character())
})

test_that("protected basin and extrema construction surface is unchanged", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(root, "tools/cleanup_ledger_lib.R")),
        "source-only ownership audit is not shipped in the package tarball"
    )
    source(file.path(root, "tools/cleanup_ledger_lib.R"), local = TRUE)

    errors <- cleanup.check.protected.surface(root)
    expect_identical(errors, character())
})

test_that("retired exports do not reappear", {
    root <- normalizePath(file.path(testthat::test_path(), "../.."))
    skip_if_not(
        file.exists(file.path(root, "tools/cleanup_ledger_lib.R")),
        "source-only ownership audit is not shipped in the package tarball"
    )
    source(file.path(root, "tools/cleanup_ledger_lib.R"), local = TRUE)

    retired <- readLines(
        file.path(root, "split_audit/cleanup/retired-exports.txt"),
        warn = FALSE
    )
    retired <- retired[
        nzchar(trimws(retired)) & !grepl("^#", trimws(retired))
    ]
    inventory <- cleanup.namespace.inventory(root)
    exports <- inventory$item_name[inventory$item_type == "export"]

    expect_length(intersect(retired, exports), 0L)
})

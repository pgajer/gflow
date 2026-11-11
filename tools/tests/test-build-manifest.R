#!/usr/bin/env Rscript
# Check manifest stability across collation and sensitivity to source changes.
run_test <- function() {
    builder <- normalizePath("tools/update_gflow_build_manifest.R")
    identity <- normalizePath("R/basin_identity.R")
    original <- getwd()
    collation <- Sys.getlocale("LC_COLLATE")
    fixture <- tempfile("gflow-manifest-test-")
    dir.create(file.path(fixture, "R"), recursive = TRUE)
    on.exit({
        setwd(original)
        Sys.setlocale("LC_COLLATE", collation)
        unlink(fixture, recursive = TRUE)
    })
    file.copy(identity, file.path(fixture, "R", "basin_identity.R"))
    setwd(fixture)
    writeLines(c("Package: example", "Version: 0.0.0"), "DESCRIPTION")
    writeLines("increment <- function(x) x + 1", "R/fixture.R")
    update <- function() suppressMessages(sys.source(builder, envir = new.env()))
    path <- file.path("inst", "extdata", "gflow-code-manifest.tsv")
    update()
    rows <- utils::read.delim(path, stringsAsFactors = FALSE)
    rows <- rows[rev(seq_len(nrow(rows))), ]
    utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
    digest <- unname(tools::md5sum(path))
    Sys.setlocale("LC_COLLATE", "C")
    update()
    stopifnot(identical(unname(tools::md5sum(path)), digest))
    writeLines("increment <- function(x) x + 2", "R/fixture.R")
    update()
    changed <- utils::read.delim(path, stringsAsFactors = FALSE)
    stopifnot(!identical(unname(tools::md5sum(path)), digest),
              identical(changed$md5[changed$path == "R/fixture.R"],
                        unname(tools::md5sum("R/fixture.R"))))
    writeLines("square <- function(x) x^2", "R/additional.R")
    update()
    stopifnot("R/additional.R" %in% utils::read.delim(path)$path)
    unlink("R/additional.R")
    update()
    stopifnot(!"R/additional.R" %in% utils::read.delim(path)$path)
}
run_test()
cat("Manifest checks passed: stable row order; changed, added and removed files detected.\n")

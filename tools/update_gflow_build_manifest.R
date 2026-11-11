source(file.path("R", "basin_identity.R"), local = environment())

path <- file.path("inst", "extdata", "gflow-code-manifest.tsv")
current <- .gflow.compute.code.manifest(".")
recorded <- if (file.exists(path))
    utils::read.delim(path, stringsAsFactors = FALSE) else NULL
# Preserve the stored row order when only the system collation differs.
# A changed, added or removed source file still produces a new manifest.
same <- !is.null(recorded) &&
    identical(sort(names(recorded)), sort(names(current))) &&
    !anyDuplicated(recorded$path) &&
    setequal(recorded$path, current$path) &&
    identical(recorded$md5[match(current$path, recorded$path)], current$md5)
if (same) {
    message("Source manifest matches the current files.")
} else {
    path <- .gflow.write.code.manifest(".")
    message("Updated ", path)
}

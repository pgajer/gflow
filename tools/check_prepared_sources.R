#!/usr/bin/env Rscript
# Verify the prepared source and vignette results without regenerating files.
source(file.path("R", "basin_identity.R"), local = environment())
manifest_file <- file.path("inst", "extdata", "gflow-code-manifest.tsv")
recorded <- utils::read.delim(manifest_file, stringsAsFactors = FALSE)
current <- .gflow.compute.code.manifest(".")
changed <- union(setdiff(recorded$path, current$path),
                 setdiff(current$path, recorded$path))
shared <- intersect(recorded$path, current$path)
changed <- union(changed, shared[
    recorded$md5[match(shared, recorded$path)] !=
        current$md5[match(shared, current$path)]])
if (anyDuplicated(recorded$path) || length(changed)) {
    stop("Source manifest does not match the prepared sources: ",
         paste(changed, collapse = ", "),
         ". Regenerate the manifest and vignette results before building.")
}
results <- readRDS(file.path("inst", "extdata", "vignette-results", "knn-basins.rds"))
if (!identical(readLines(manifest_file), results$provenance$code.manifest) ||
    !identical(unname(tools::md5sum("vignettes/knn-basins.Rmd")),
               results$provenance$vignette.md5)) {
    stop("Vignette results are stale; regenerate with vignettes/precompute-knn-basins.R.")
}
message("Prepared sources and vignette results match.")

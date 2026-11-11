#!/usr/bin/env Rscript
# Render the matrix-input vignette against optimized local package sources.
# Build products, cached fits, and numerical verification stay under build/.
root <- normalizePath(".")
out <- file.path(root, "build", "knn-basins")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
           VECLIB_MAXIMUM_THREADS = 1, MAKEFLAGS = "-j1")
dependency <- file.path(dirname(root), "dgraphs")
pkgbuild::compile_dll(dependency, debug = FALSE)
pkgload::load_all(dependency, recompile = FALSE, quiet = TRUE)
pkgbuild::compile_dll(root, debug = FALSE)
pkgload::load_all(root, recompile = FALSE, quiet = TRUE)
identity <- gflow::get.gflow.build.identity(refresh = TRUE)
saveRDS(identity, file.path(out, "build-identity.rds"))
environment <- new.env(parent = globalenv())
environment$knn_artifact_dir <- out
environment$knn_build_identity <- identity
# This source preview recomputes all comparisons against the loaded implementation.
options(gflow.recompute.vignettes = TRUE)
result <- rmarkdown::render(
    file.path(root, "vignettes", "knn-basins.Rmd"),
    output_dir = out, intermediates_dir = file.path(out, "intermediates"),
    envir = environment, quiet = FALSE, clean = FALSE
)
if (!isTRUE(environment$numerical_checks_passed)) {
    stop("Vignette numerical verification did not finish.")
}
write.csv(environment$stage_table, file.path(out, "stage-counts.csv"), row.names = FALSE)
write.csv(environment$sensitivity_table, file.path(out, "sensitivity.csv"), row.names = FALSE)
write.csv(environment$terrain_mode_table, file.path(out, "terrain-mode-comparison.csv"), row.names = FALSE)
write.csv(environment$terrain_interior_table, file.path(out, "terrain-interior-comparison.csv"), row.names = FALSE)
write.csv(environment$synthetic_results, file.path(out, "synthetic-mode-results.csv"), row.names = FALSE)
write.csv(environment$synthetic_summary, file.path(out, "synthetic-mode-summary.csv"), row.names = FALSE)
write.csv(environment$synthetic_fragments, file.path(out, "synthetic-fragment-example.csv"), row.names = FALSE)
write.csv(environment$synthetic_reference_checks, file.path(out, "synthetic-reference-checks.csv"), row.names = FALSE)
write.csv(environment$verification, file.path(out, "numerical-verification.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
saveRDS(list(one = environment$one, noisy = environment$noisy,
             two = environment$two, terrain = environment$terrain,
             terrain.asymmetric = environment$terrain_asymmetric,
             synthetic.comparison = environment$synthetic_cases,
             density = environment$density_control,
             sensitivity = environment$sensitivity_fits),
        file.path(out, "vignette-results.rds"))
cat("\nRendered and numerically verified:", result, "\n")

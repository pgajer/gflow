#!/usr/bin/env Rscript
# Run from the package source directory after installing the current sources.
input <- normalizePath("vignettes/knn-basins.Rmd", mustWork = TRUE)
output <- file.path("inst", "extdata", "vignette-results", "knn-basins.rds")
installed_manifest <- system.file("extdata", "gflow-code-manifest.tsv", package = "gflow")
source_manifest <- file.path("inst", "extdata", "gflow-code-manifest.tsv")
if (!identical(readLines(source_manifest), readLines(installed_manifest)))
    stop("Install the current package sources before regenerating vignette results.")
previous <- options(gflow.recompute.vignettes = TRUE)
work <- tempfile("knn-basins-render-")
dir.create(work)
evaluation <- new.env(parent = globalenv())
evaluation$knn_artifact_dir <- file.path(work, "artifacts")
compact_fit <- function(fit) {
    # Retain the values used by figures, tables and independent consistency checks.
    list(max = fit$max, min = fit$min, y = fit$y, parameters = fit$parameters,
         neighborhood = list(n.components = fit$neighborhood$n.components))
}
tryCatch({
    rmarkdown::render(input, output_dir = work, intermediates_dir = work,
                      quiet = TRUE, envir = evaluation)
    stopifnot(isTRUE(evaluation$numerical_checks_passed))
    cases <- lapply(evaluation$synthetic_cases, function(z) {
        z$fit <- compact_fit(z$fit)
        z
    })
    sensitivity <- lapply(evaluation$sensitivity_fits, compact_fit)
    provenance <- list(
        vignette.md5 = unname(tools::md5sum(input)),
        code.manifest = readLines(installed_manifest),
        package.versions = vapply(c("gflow", "dgraphs", "Rcpp", "knitr", "rmarkdown"),
            function(package) as.character(utils::packageVersion(package)), character(1)),
        R.version = R.version.string, platform = R.version$platform,
        seeds = evaluation$synthetic_seeds, observations.per.synthetic.case = 2500L,
        synthetic.settings = list(k = 50L, metric = "L2", direction = "max",
            min.persistence = .15, min.basin.size = 100L,
            surfaces = c("Two peaks", "Four peaks"), outcomes = c("Exact", "Rounded"),
            modes = c("Symmetric", "Asymmetric")),
        sensitivity.settings = lapply(evaluation$sensitivity_inputs,
            function(x) list(observations = nrow(x$X), persistence = x$p,
                             metrics = x$metrics, k = c(25L, 50L, 100L),
                             min.basin.size = 100L)),
        paired.inputs.checked = 32L, sign.reversal.checked = 64L)
    payload <- list(format.version = 1L, provenance = provenance,
        seeds = evaluation$synthetic_seeds, cases = cases,
        equalities = evaluation$synthetic_equalities, sensitivity = sensitivity)
    stopifnot(length(cases) == 64L, nrow(payload$equalities) == 32L,
              length(sensitivity) == 15L)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("knn-basins-", tmpdir = dirname(output), fileext = ".rds")
    saveRDS(payload, temporary, version = 3, compress = "xz")
    if (!file.rename(temporary, output)) stop("Could not replace the stored results.")
    message("Saved ", output, " (", file.info(output)$size, " bytes)")
}, finally = {
    options(previous)
    unlink(work, recursive = TRUE)
})

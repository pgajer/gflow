#!/usr/bin/env Rscript
# Run from a source checkout, using the same R libraries as `make install`.
metadata <- read.dcf("DESCRIPTION")[1L, ]
fields <- metadata[intersect(c("Depends", "Imports", "LinkingTo"), names(metadata))]
required <- unique(trimws(sub("\\s*\\(.*", "", unlist(strsplit(paste(fields, collapse = ","), ",")))))
required <- setdiff(required, "R")
required <- unique(c(required, "roxygen2", "knitr", "rmarkdown"))
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (any(!vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Some build dependencies are unavailable. Check installation messages and R library permissions.")
}
if (!rmarkdown::pandoc_available()) {
    stop("Pandoc is required to build installed guides. Install Pandoc or run with RStudio's Pandoc on PATH.")
}
if (packageVersion("dgraphs") < "0.2.0" || length(unclass(packageVersion("dgraphs"))[[1L]]) > 3L) {
    stop("Use published dgraphs 0.2.0 or a compatible release for this build; the development graph API is not supported. Install the CRAN release in the selected R library.")
}
cat("Build dependencies and Pandoc are available. Next: make install R_ENV=\"env -u R_HOME\"\n")

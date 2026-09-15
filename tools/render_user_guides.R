#!/usr/bin/env Rscript
# Build static, self-contained previews from an installed package.
root <- normalizePath(".")
out <- file.path(root, "build", "vignettes")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (source in list.files("vignettes", "[.]Rmd$", full.names = TRUE)) {
    rmarkdown::render(source, output_dir = out,
                      envir = new.env(parent = globalenv()), quiet = TRUE)
}
cat("Rendered installed-package workflows into build/vignettes.\n")

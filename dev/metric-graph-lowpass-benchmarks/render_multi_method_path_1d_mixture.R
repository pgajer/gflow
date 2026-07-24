script.arg <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
script.dir <- if (length(script.arg)) dirname(script.arg) else getwd()
repo.dir <- normalizePath(file.path(script.dir, "..", ".."), mustWork = TRUE)
geosmooth.dir <- normalizePath(file.path(repo.dir, "..", "geosmooth"), mustWork = TRUE)
report.dir <- file.path(geosmooth.dir, "dev", "metric-graph-lowpass-benchmarks")
out.dir <- file.path(report.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(report.dir)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
}

rmarkdown::render(
  input = file.path(report.dir, "multi_method_path_1d_mixture.Rmd"),
  output_file = file.path(out.dir, "multi_method_path_1d_mixture.html"),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

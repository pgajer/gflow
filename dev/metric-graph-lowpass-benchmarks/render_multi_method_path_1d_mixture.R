repo.dir <- "/Users/pgajer/current_projects/gflow"
report.dir <- file.path(repo.dir, "dev/metric-graph-lowpass-benchmarks")
out.dir <- file.path(report.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(report.dir)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
}

rmarkdown::render(
  input = "multi_method_path_1d_mixture.Rmd",
  output_file = file.path(out.dir, "multi_method_path_1d_mixture.html"),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

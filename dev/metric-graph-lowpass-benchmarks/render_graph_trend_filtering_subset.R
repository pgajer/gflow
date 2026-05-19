repo.dir <- "/Users/pgajer/current_projects/gflow"
report.dir <- file.path(repo.dir, "dev/metric-graph-lowpass-benchmarks")
out.dir <- file.path(report.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(report.dir, "graph_trend_filtering_subset.Rmd"),
  output_file = file.path(out.dir, "graph_trend_filtering_subset.html"),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

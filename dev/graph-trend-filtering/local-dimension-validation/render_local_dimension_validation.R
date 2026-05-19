repo.dir <- "/Users/pgajer/current_projects/gflow"
validation.dir <- file.path(repo.dir, "dev/graph-trend-filtering/local-dimension-validation")
out.dir <- file.path(validation.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(validation.dir, "compute_local_dimension_validation.R"))

mode <- Sys.getenv("GFLOW_LOCAL_DIM_VALIDATION_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))
method <- Sys.getenv("GFLOW_LOCAL_DIM_VALIDATION_METHOD", unset = "cmdscale")
method <- match.arg(method, c("cmdscale", "metric.mds.edge.kk", "hybrid.edge.kk"))

diagnostics <- run.local.dimension.validation(mode = mode, output.dir = out.dir,
                                              method = method)
message("Computed ", mode, " local dimension validation with ",
        nrow(diagnostics$center.selection), " sampled centers.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(validation.dir, "local_dimension_validation.Rmd"),
  output_file = file.path(out.dir, paste0("local_dimension_validation_", mode, ".html")),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

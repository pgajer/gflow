repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "transported_hessian_quadric_coordinate_validation.R"))

mode <- Sys.getenv("GFLOW_QUADRIC_COORDINATE_VALIDATION_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))

diagnostics <- run.transported.hessian.quadric.coordinate.validation(
  mode = mode,
  repo.dir = repo.dir,
  output.dir = out.dir
)
message("Computed transported Hessian quadric coordinate validation with ",
        nrow(diagnostics$summary), " graph/policy runs.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir, "transported_hessian_quadric_coordinate_validation.Rmd"),
  output_file = file.path(out.dir, paste0("transported_hessian_quadric_coordinate_validation_", mode, ".html")),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "transported_hessian_quadric_multimethod_validation.R"))

mode <- Sys.getenv("GFLOW_QUADRIC_MULTIMETHOD_VALIDATION_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))

diagnostics <- run.transported.hessian.quadric.multimethod.validation(
  mode = mode,
  repo.dir = repo.dir,
  output.dir = out.dir,
  graph.weights = "ambient"
)
message("Computed transported Hessian quadric multi-method validation with ",
        nrow(diagnostics$summary), " graph/method runs.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir,
                    "transported_hessian_quadric_multimethod_validation.Rmd"),
  output_file = file.path(
    out.dir,
    paste0("transported_hessian_quadric_multimethod_validation_", mode, ".html")
  ),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

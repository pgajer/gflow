repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "geodesic_annihilator_experiments.R"))

mode <- Sys.getenv("GFLOW_GEO_ANNIHILATOR_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))

diagnostics <- run.geodesic.annihilator.experiment(mode = mode, output.dir = out.dir)
message("Computed ", mode, " geodesic annihilator experiment with ",
        nrow(diagnostics$results), " runs.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir, "geodesic_annihilator_nullspace_experiments.Rmd"),
  output_file = file.path(out.dir, paste0("geodesic_annihilator_", mode, ".html")),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

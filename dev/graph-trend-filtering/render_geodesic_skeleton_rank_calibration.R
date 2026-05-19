repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "geodesic_skeleton_rank_calibration.R"))

mode <- Sys.getenv("GFLOW_GEO_SKELETON_RANK_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))

diagnostics <- run.geodesic.skeleton.rank.calibration(
  mode = mode,
  output.dir = out.dir
)
message("Computed ", mode, " geodesic skeleton rank calibration with ",
        nrow(diagnostics$summary), " runs.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir, "geodesic_skeleton_rank_calibration.Rmd"),
  output_file = file.path(out.dir, paste0("geodesic_skeleton_rank_calibration_", mode, ".html")),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)


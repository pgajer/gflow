repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "transported_hessian_adaptive_threshold_sweep.R"))

mode <- Sys.getenv("GFLOW_ADAPTIVE_THRESHOLD_SWEEP_MODE", unset = "smoke")
mode <- match.arg(mode, c("smoke", "full"))
use.local.grip <- tolower(Sys.getenv("GFLOW_USE_LOCAL_GRIP", unset = "true")) %in%
  c("true", "1", "yes", "y")

diagnostics <- run.transported.hessian.adaptive.threshold.sweep(
  mode = mode,
  repo.dir = repo.dir,
  output.dir = out.dir,
  use.local.grip = use.local.grip
)
message("Computed transported Hessian adaptive-threshold sweep with ",
        nrow(diagnostics$summary), " method/case/filter runs.")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir, "transported_hessian_adaptive_threshold_sweep.Rmd"),
  output_file = file.path(out.dir, paste0("transported_hessian_adaptive_threshold_sweep_", mode, ".html")),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

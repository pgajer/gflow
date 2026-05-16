repo.dir <- "/Users/pgajer/current_projects/gflow"
report.dir <- file.path(repo.dir, "dev/metric-graph-lowpass-benchmarks")
input <- file.path(report.dir, "single_scenario_smoke.Rmd")
output.dir <- file.path(report.dir, "reports")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("The rmarkdown package is required to render this report.", call. = FALSE)
}

rmarkdown::render(
    input = input,
    output_file = "single_scenario_smoke.html",
    output_dir = output.dir,
    knit_root_dir = repo.dir,
    envir = new.env(parent = globalenv()),
    quiet = FALSE
)

repo.dir <- "/Users/pgajer/current_projects/gflow"
report.dir <- file.path(repo.dir, "dev/metric-graph-lowpass-benchmarks")
out.dir <- file.path(report.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
} else {
  library(gflow)
}

source(file.path(report.dir, "multi_scenario_1d_mixture_helpers.R"))

results <- run_multi_method_benchmark(
  repo.dir = repo.dir,
  graph.type = "path",
  replicates = 1L,
  shape.ids = c("S01", "S05"),
  variant.ids = "V1",
  progress = TRUE
)

csv.file <- file.path(out.dir, "multi_method_path_1d_mixture_smoke_results.csv")
rds.file <- file.path(out.dir, "multi_method_path_1d_mixture_smoke_results.rds")
write.csv(results, csv.file, row.names = FALSE)
saveRDS(results, rds.file)

message(sprintf("Wrote %s", csv.file))
message(sprintf("Wrote %s", rds.file))
message(sprintf("Rows: %d; datasets: %d; errors: %d",
                nrow(results), length(unique(results$dataset.id)),
                sum(results$status != "ok")))

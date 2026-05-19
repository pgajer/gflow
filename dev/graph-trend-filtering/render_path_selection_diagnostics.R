repo.dir <- "/Users/pgajer/current_projects/gflow"
dev.dir <- file.path(repo.dir, "dev/graph-trend-filtering")
out.dir <- file.path(dev.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(dev.dir, "path_selection_diagnostics.R"))
diagnostics <- build_phase4_path_selection_diagnostics()
saveRDS(diagnostics, file.path(out.dir, "path_selection_diagnostics.rds"))
utils::write.csv(diagnostics$summaries,
                 file.path(out.dir, "path_selection_diagnostics_summary.csv"),
                 row.names = FALSE)
utils::write.csv(diagnostics$vertex.coverage,
                 file.path(out.dir, "path_selection_vertex_coverage.csv"),
                 row.names = FALSE)
utils::write.csv(diagnostics$path.tables,
                 file.path(out.dir, "path_selection_path_table.csv"),
                 row.names = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render this report.")
}

rmarkdown::render(
  input = file.path(dev.dir, "path_selection_diagnostics.Rmd"),
  output_file = file.path(out.dir, "path_selection_diagnostics.html"),
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)

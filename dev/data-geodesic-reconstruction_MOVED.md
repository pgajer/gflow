# Data Geodesic Reconstruction Assets Moved

The project-layer data geodesic reconstruction assets that were previously
tracked under `gflow/dev/data-geodesic-reconstruction/` now live in:

`/Users/pgajer/current_projects/dgraphs_manuscripts/geodesic_data_geometry/experiments/quadform-benchmarks/`

Curated reports and manuscript-style records are kept in:

`/Users/pgajer/current_projects/dgraphs_manuscripts/geodesic_data_geometry/reports/curated/`

The reusable implementation is now owned by the `dgraphs` package: graph
constructors, pruning methods, reference geodesic utilities, and isometry
diagnostic helpers.

Large raw run outputs, `.rds` graph/layout assets, caches, and widget payloads
are derived artifacts and should not be version controlled in `gflow`.

Older report exports recovered on September 12, 2026 are preserved separately
under `dgraphs_manuscripts/geodesic_data_geometry/archive/gflow-20260912/`.
Unmoved raw run directories remain local; see the manuscript relocation record.

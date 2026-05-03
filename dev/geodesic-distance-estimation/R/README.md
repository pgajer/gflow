# Benchmark Helpers

Place reusable R helper source here. Suggested modules:

- `graph_construction.R`: wrappers for CMST, iKNN, oracle graphs, and parameter
  sweeps.
- `metrics.R`: geodesic distortion, false-shortcut rate, neighborhood
  preservation, connectivity, and degree summaries.
- `preprocessing.R`: CLR/Aitchison-style preprocessing for compositional
  examples.
- `plotting.R`: geometry panels, graph overlays, distance scatterplots, and
  benchmark summaries.

Keep helpers source-only. Generated outputs should be written under
`../reports/`, `../figures/`, or `../cache/`.

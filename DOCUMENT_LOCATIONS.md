# Document locations

Scientific documents are kept in sibling manuscript workspaces, not in this
R package checkout:

- `../gflow_manuscripts/`: graph-Laplacian notes and research presentations.
- `../dgraphs_manuscripts/geodesic_data_geometry/`: graph-geometry studies.
  Older recovered exports are in `archive/gflow-20260912/`; they do not replace
  the current manuscript sources.
- `../geosmooth_manuscripts/reports/metric-graph-lowpass/from-gflow/`:
  historical smoothing comparisons and their scripts.

These are local sibling workspace paths, not package runtime dependencies.
User help, examples, test fixtures, and public API ownership metadata remain
here. Cleanup tables and the migration documents consumed by package guardrail
checks also remain here; those checks do not require private documents.
Old run directories were not swept into this document move.

New scientific writing belongs in the appropriate manuscript workspace;
agent correspondence, prompts, execution checklists, and internal reviews
belong in the user's private project storage. No compatibility symlinks to
external documents are required to build or test the package.

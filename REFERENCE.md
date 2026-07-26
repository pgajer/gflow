# gflow public API reference

The supported `gflow` API has four user-facing families. Basin and extrema
construction internals are intentionally outside the cleanup that produced
this index; the canonical basin-complex API remains the starting point.

## Basin and flow objects

- `create.basin.complex()` constructs the canonical `basin_complex` object.
- `as.basin.complex()` converts supported archived construction results.
- `get.basin.table()`, `get.basin.membership()`, and
  `get.basin.assignment()` expose stable basin tables.
- `get.basin.merge.tree()`, `get.basin.trajectory.forest()`, and
  `get.basin.cells()` expose method-specific structures without requiring
  users to inspect internal lists.
- Standard `print()`, `summary()`, `plot()`, and `as.data.frame()` methods
  describe canonical objects.

## Complex and trajectory exploration

- `construct.gflow.graph()` summarizes basin intersections as a flow graph.
- `construct.madag()` and `madag.bottlenecks()` analyze directed basin/cell
  structure.
- `compute.harmonic.extension()`, `apply.harmonic.extension()`,
  `analyze.harmonic.extensions()`, and `compare.harmonic.methods()` extend and
  compare coordinates around trajectories.
- `select.max.density.trajectory()` selects a representative trajectory.
- `compute.gfc.modulation()` applies a named modulation to a gradient-flow
  complex.
- `extremality.summary()` and `label.extremality.3d()` support extrema-focused
  interpretation.

## Local association

- `lcor()` estimates symmetric graph-local association.
- `lslope()` estimates directed graph-local response.
- `lslope.neighborhood()` reports neighborhood-level directed response.
- `lcor.with.posterior()` summarizes supplied posterior field draws; it does
  not fit a conditional-expectation model.
- `permutation.test.lcor()` provides permutation inference.

## Flow-aware association

- `gfcor()` computes global and basin-aware flow association.
- `gfassoc.membership()`, `gfassoc.polarity()`, `gfassoc.overlap()`, and
  `gfassoc.deviation()` expose the component flow-aware summaries.

## Explicitly outside this package

- Graph construction, conversion, paths, components, spectral routines,
  endpoint diagnostics, graph selection, and generic graph plotting:
  `dgraphs`.
- Graph regression, response smoothing, conditional-expectation estimation,
  PHATE, diffusion/potential pseudotime, and quadratic-form geodesic
  experiments: not part of the supported `gflow` API. Archived estimators live
  in `gflowx`; experimental families without a verified successor are recorded
  as removed in the migration guide.
- Interactive selection widgets and domain-specific analysis pipelines:
  specialist UI or analysis repositories, not `gflow`.

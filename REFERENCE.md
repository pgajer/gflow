# gflow public API reference

Start with [Finding your way around gflow](https://pgajer.github.io/gflow/articles/function-guide.html) for
an exhaustive, status-labeled catalog and [Example graphs and scalar
fields](https://pgajer.github.io/gflow/articles/example-graphs-and-fields.html) for runnable inputs.
Both are installed vignettes. The canonical basin-complex API is the starting
point for new basin analyses; the families below also contain advanced and
archived-object interfaces.

## Basin and flow objects

- `detect.adaptive.extrema()` detects graph extrema in adaptive neighborhoods;
  `summary()`, `plot()`, and `vertices()` operate on its `gflow_local_extrema`
  results. `dgraphs::detect.local.extrema()` has different, fixed-radius semantics.
- `create.basin.complex()` constructs the canonical `basin_complex` object.
- `as.basin.complex()` converts supported archived construction results.
- `get.basin.table()`, `get.basin.membership()`, and
  `get.basin.assignment()` expose stable basin tables.
- `get.basin.merge.tree()` returns a complete canonical merge-tree object;
  `plot.basin.merge.tree()`, `cut.basin.merge.tree()`, and
  `as.dendrogram.basin.merge.tree()` visualize, query, and coerce its exact
  graph level-set hierarchy.
- `get.basin.trajectory.forest()` and `get.basin.cells()` expose other
  method-specific structures without requiring users to inspect internal
  lists.
- Standard `print()`, `summary()`, `plot()`, and `as.data.frame()` methods
  describe canonical objects.

## Complex and trajectory exploration

- `construct.gflow.graph()` summarizes archived basin intersections as a flow graph.
- `construct.madag()` and `madag.bottlenecks()` analyze directed basin/cell
  structure.
- `compute.harmonic.extension()`, `apply.harmonic.extension()`,
  `analyze.harmonic.extensions()`, and `compare.harmonic.methods()` extend and
  compare coordinates around trajectories.
- `select.max.density.trajectory()` selects a representative trajectory.
- `compute.gfc.modulation()` computes edge modulation weights from graph
  lengths and optional density.
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

- `gfcor()` computes global and basin-aware flow association from two fields
  and two archived `basins_of_attraction` objects. It currently rejects
  `basin_complex` objects; no public canonical adapter exists.
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

## Migration

Use `help("gflow-migration", package = "gflow")` or the migration section of
the function guide for installed-package advice. Exported retirement stubs
raise migration errors and must not be used as estimators.

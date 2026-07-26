# Phase 2 smoothing and conditional-expectation boundary

## Ownership decisions

| Family | Owner after Phase 2 | Action |
|---|---|---|
| `lslope.test()` and its S3/extractor family | `gflowx` | Move. The API requires a `knn.riem.fit` and performs spectral refits. |
| `subject.neighborhood.stats()` and rdgraph weight helpers | `gflowx` | Move. The input contract is an archived rdgraph fit. |
| fitted-model posterior `lcor` adapter | `gflowx` | Replace the native mixed-mode path with `lcor.with.rdgraph.posterior()`, which obtains draws through `refit.rdgraph.regression()` and delegates association summaries to `gflow`. |
| `lcor.with.posterior()` | `gflow` | Retain only the model-independent estimand over supplied vertex-field draws. Remove fitted-model and smoothing parameters from its public contract. |
| `permutation.test.lcor()` | `gflow` | Retain as model-independent local-association inference; remove the archived-fit example. |
| `gflow.smooth.spline()` | internal `gflow` helper | De-export. It remains private only for trajectory and diagnostic plotting code. |

## Protected references

The remaining optional `gflowx` hooks in extrema and basin/complex code are
outside this cleanup by explicit project scope. In particular,
`graph_endpoint_geometry.R`, `extremality_utils.R`,
`compute_pextrema_nbhds.R`, `compute_gfc_modulation.R`, and
`vertex_density.R` are marked protected for the later basin/extrema merge.
They justify retaining `gflowx` in `Suggests` for now.

## Boundary invariant

No in-scope exported `gflow` function accepts or requires a
`knn.riem.fit`, `riem.dcx`, `fit.rdgraph.regression()`, or
`refit.rdgraph.regression()` object. Sampling and smoothing are performed by
their owning package; `gflow` receives graph fields or already-generated
draws for local and flow-aware analysis.

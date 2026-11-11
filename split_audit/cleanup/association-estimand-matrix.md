# Association estimand matrix

This table distinguishes the association estimands and their public interfaces. “Local”
means that the estimand is evaluated at graph vertices from a graph
neighborhood; “flow-aware” means that basin, cell, polarity, or trajectory
structure is part of the estimand.

| Family | Estimand | Locality | Input contract | Canonical public API | Result contract | Availability |
|---|---|---|---|---|---|---|
| `lcor` | Symmetric alignment of edge differences | vertex-local, optionally multi-hop | adjacency/weight lists plus vertex field(s) | `lcor()` | numeric vector or classed matrix/list/array according to input shape | public; dispatch helpers are internal |
| `lslope` | Directed response slope along the local gradient of `y` | vertex-local and direction-aware | adjacency/weight lists plus directing and response fields | `lslope()` | numeric vector or classed matrix/array; instrumented vector mode is `lslope_gradient_result` | public |
| neighborhood slope | Regression slope over all local incident edges | vertex-local | adjacency/weight lists plus two fields | `lslope.neighborhood()` | `lslope_neighborhood_result` | public, with a distinct regression estimand |
| posterior `lcor` | Distribution of `lcor` over supplied field draws | vertex-local uncertainty | graph fields, fixed response, sampled feature fields | `lcor.with.posterior()` | `lcor.posterior` | public and model-independent; fitted-model adapter is in `gflowx` |
| permutation `lcor` | Feature-level permutation significance of a chosen aggregate of `lcor` | global inference over vertex-local scores | graph fields, response/features, permutation controls | `permutation.test.lcor()` | `lcor_permutation_test` | public |
| paired `lslope` | Paired-bootstrap inference after rdgraph smoothing | vertex-local but smoother-dependent | `knn.riem.fit`, raw response/feature | `gflowx::lslope.test()` | `lslope.test` | archived in `gflowx` |
| `gfcor` | Polarity concordance between two basin decompositions | vertex-, basin-, and global flow-aware | two fields plus their basin structures | `gfcor()` | `gfcor` | public |
| `gfassoc.*` components | Membership, polarity, overlap, and deviation diagnostics | basin/cell flow-aware | basin structures or derived membership objects | `gfassoc.membership()`, `gfassoc.polarity()`, `gfassoc.overlap()`, `gfassoc.deviation()` | named classed structures or matrices | public advanced decomposition APIs, each with a distinct intermediate estimand |
| trajectory hurdle display | Detection-aware association display along selected paths | trajectory-aware | response, trajectory coordinate, paths | `plot.trajectory.hurdle.association()` | invisible long-form data plus optional fit metadata | public trajectory analysis display, not a second numerical estimator |
| `fassoc`, `fassoc0`, `fassoc1` | Generic one-dimensional conditional-mean/slope association | not graph-local or flow-aware | paired numeric vectors plus smoother controls | none in `gflow` | internal `assoc0`/`assoc1` objects | archived in `gflowx`; removed from `gflow` |

## Distinctions between interfaces

- `lcor.vector.matrix()` and `lcor.matrix.matrix()` are dispatch
  implementations, not separate estimands. They are private.
- `lslope()` and `lslope.neighborhood()` are not aliases: one follows the
  selected gradient edge and the other fits over the full neighborhood.
- `gfcor()` is the normal user entry point. The `gfassoc.*` functions expose
  interpretable intermediate flow quantities for advanced analysis and testing.
- Posterior and permutation functions are inference layers around `lcor()`;
  they do not redefine the coefficient.
- The generic `fassoc*` family is not part of the package’s local/flow-aware
  scope and is now archived in `gflowx`.

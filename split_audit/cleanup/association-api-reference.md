# Association API reference

Use `lcor()` for symmetric graph-local association and `lslope()` when one
field directs the local gradient and the other is the response. Both dispatch
over vector and matrix inputs; callers should use these front doors rather than
shape-specific helpers.

Use `lcor.with.posterior()` when draws of a vertex field already exist, and
`permutation.test.lcor()` for feature-level permutation inference. Archived
rdgraph fits are adapted by `gflowx::lcor.with.rdgraph.posterior()`; they are
not accepted by the core functions.

Use `gfcor()` for a complete flow-aware association analysis between two
basin decompositions. The `gfassoc.*` functions are advanced decomposition
tools for membership, polarity, overlap, and deviation from independence.

Use `lslope.neighborhood()` only when the full-neighborhood regression
estimand is intended; it is not an alias for gradient-restricted `lslope()`.

`fassoc.test()`, `fassoc0.test()`, and `fassoc1.test()` are legacy generic
conditional-mean methods and are internal. Fit-dependent paired local-slope
testing is archived as `gflowx::lslope.test()`.

## Stable arguments

- Graph-local methods start with `adj.list`, `weight.list`.
- Symmetric association uses `y`, `z`; directed slope treats `y` as directing
  and `z` as response.
- Difference controls use `y.diff.type`, `z.diff.type`.
- Matrix inputs retain one row per graph vertex.
- Inference functions use `seed` and an explicit sample/permutation count.

## Stable results

- `lcor()` and `lslope()` return one vertex-level coefficient per requested
  field pair, with classed matrix/array results for multi-field inputs.
- `lcor.with.posterior()` returns `mean`, `sd`, `lower`, and `upper` matrices
  in a `lcor.posterior` object.
- `permutation.test.lcor()` returns a ranked feature table plus named observed,
  p-value, and q-value vectors in a `lcor_permutation_test` object.
- `gfcor()` returns `global`, `vertex`, `basin_character`, `overlap`, and
  `membership` components in a `gfcor` object.

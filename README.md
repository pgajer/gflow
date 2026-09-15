# gflow

`gflow` constructs, explores, and analyzes basin and gradient-flow complexes
on structured high-dimensional data. Its supported public API is deliberately
narrow: basin/flow objects, their summaries and trajectories, graph-local
association, and flow-aware association.

Generic graph construction and graph algorithms belong to
[`dgraphs`](https://github.com/pgajer/dgraphs). Retired response-smoothing and
conditional-expectation estimators are archived in
[`gflowx`](https://github.com/pgajer/gflowx); they are not dependencies of
`gflow`.

## Core workflow

```r
library(gflow)
adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
edge_lengths <- lapply(adjacency, function(v) rep(1, length(v)))
field <- c(0, 3, 1, 2, 0)
bc <- create.basin.complex(
  adjacency, edge_lengths, field,
  method = "superlevel_merge_tree", direction = "max"
)
summary(bc)
get.basin.table(bc)
get.basin.membership(bc)
plot(get.basin.merge.tree(bc))
```

For data-derived graphs, choose a documented constructor in `dgraphs` and
supply its adjacency and aligned edge-length lists.

## User guides

- [Finding your way around gflow](vignettes/function-guide.Rmd): task map,
  complete export catalog, method availability, and installed migration advice.
- [Example graphs and scalar fields](vignettes/example-graphs-and-fields.Rmd):
  paths, a grid, and disconnected components with reproducible comparisons.
- [Canonical basin workflow](vignettes/basin_complex_workflow_vignette.Rmd).
- [Noisy-circle workflow](vignettes/noisy_circle_core_workflow_vignette.Rmd).

After installing a build with vignettes, use
`vignette("function-guide", package = "gflow")` or
`vignette("example-graphs-and-fields", package = "gflow")` for rendered guides.
`help("gflow-migration", package = "gflow")` provides installed migration help.

## Supported API map

| Purpose | Canonical entry points |
|---|---|
| Construct/convert basin complexes | `create.basin.complex()`, `as.basin.complex()` |
| Inspect basin objects | `summary()`, `plot()`, `get.basin.table()`, `get.basin.membership()`, `get.basin.assignment()` |
| Explore trajectories and cells | `get.basin.trajectory.forest()`, `get.basin.cells()`, `compute.harmonic.extension()`, `construct.madag()` |
| Local association | `lcor()`, `lslope()`, `lslope.neighborhood()`, `permutation.test.lcor()` |
| Flow-aware association on archived basin/membership objects | `gfcor()`, `gfassoc.membership()`, `gfassoc.polarity()`, `gfassoc.overlap()`, `gfassoc.deviation()` |

See [REFERENCE.md](REFERENCE.md) for the maintained public families and
[the public migration guidance](vignettes/function-guide.Rmd#migration-and-package-boundaries)
for removed names and verified package boundaries. `gfcor()` and
`gfassoc.membership()` currently reject canonical basin objects.

## Installation

For this source checkout, install the published `dgraphs` dependency first,
then build `gflow` with its required dependencies. The current workflows were
checked with dgraphs 0.2.0; development versions may change input contracts.

```bash
R -q -e 'install.packages("dgraphs", repos="https://cloud.r-project.org")'
R -q -e 'remotes::install_local(".", dependencies=c("Depends","Imports","LinkingTo"), upgrade="never")'
```

The default `cran-safe` profile is intended to support serial builds and uses
OpenMP when the toolchain supplies it. Portability must be checked for each
release candidate; the profile name alone is not evidence of a successful
non-OpenMP build. A performance-oriented build can require OpenMP:

```bash
R -q -e 'Sys.setenv(GFLOW_BUILD_PROFILE="dev"); remotes::install_local(".", dependencies=c("Depends","Imports","LinkingTo"), upgrade="never")'
```

Detailed toolchain instructions are in [INSTALL.md](INSTALL.md).

## Development QA

```bash
make document
make audit-api-guide
make check-fast
make check
make audit-final-acceptance
```

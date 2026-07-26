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
library(dgraphs)
library(gflow)

# Build the graph with dgraphs, then construct a basin complex in gflow.
# graph <- dgraphs::<graph-constructor>(...)
bc <- create.basin.complex(
  adj.list = graph$adj.list,
  weight.list = graph$weight.list,
  f = field,
  method = "merge_tree"
)

summary(bc)
get.basin.table(bc)
get.basin.membership(bc)
plot(bc)
```

The graph construction line is intentionally schematic because the appropriate
`dgraphs` constructor depends on the data and graph model. `gflow` consumes
adjacency and weight lists; it does not duplicate generic graph infrastructure.

## Supported API map

| Purpose | Canonical entry points |
|---|---|
| Construct/convert basin complexes | `create.basin.complex()`, `as.basin.complex()` |
| Inspect basin objects | `summary()`, `plot()`, `get.basin.table()`, `get.basin.membership()`, `get.basin.assignment()` |
| Explore trajectories and cells | `get.basin.trajectory.forest()`, `get.basin.cells()`, `compute.harmonic.extension()`, `construct.madag()` |
| Local association | `lcor()`, `lslope()`, `lslope.neighborhood()`, `permutation.test.lcor()` |
| Flow-aware association | `gfcor()`, `gfassoc.membership()`, `gfassoc.polarity()`, `gfassoc.overlap()`, `gfassoc.deviation()` |

See [REFERENCE.md](REFERENCE.md) for the maintained public families and
[the 0.2.0 migration guide](split_audit/cleanup/breaking-release-migration.md)
for removed names and successor packages.

## Installation

Install `dgraphs` first, then install `gflow` with its required dependencies:

```bash
R -q -e 'remotes::install_local("../dgraphs", dependencies=c("Depends","Imports","LinkingTo"), upgrade="never")'
R -q -e 'remotes::install_local(".", dependencies=c("Depends","Imports","LinkingTo"), upgrade="never")'
```

The default `cran-safe` build is portable and uses OpenMP when the active R
toolchain supplies it. A performance-oriented build can require OpenMP:

```bash
R -q -e 'Sys.setenv(GFLOW_BUILD_PROFILE="dev"); remotes::install_local(".", dependencies=c("Depends","Imports","LinkingTo"), upgrade="never")'
```

Detailed toolchain instructions are in [INSTALL.md](INSTALL.md).

## Development QA

```bash
make document
make check-fast
make check
make audit-final-acceptance
```

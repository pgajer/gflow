# gflow Split Audit SA0: Boundary Contract

Generated: 2026-06-03

This document defines the proposed package boundaries for splitting the current
`gflow` package into a focused graph/gradient-flow core, a geometric smoothing
package, and a GitHub-only archive/lab package.  This is an audit contract only:
no source files are moved as part of SA0/SA1.

## Proposed Package Names

- `gflow`: keep the current package name for graph construction, graph/geodesic
  infrastructure, gradient-flow objects, and Morse-Smale style graph analysis.
- `geosmooth`: new CRAN-facing package for conditional expectation,
  smoothing, regression, and local-polynomial or trend-filtering models over
  geometric supports.
- `gflowx`: public GitHub-only package for experimental, retired, highly
  application-specific, or unstable methods that should not burden the CRAN
  package surface.

R package names may contain dots, but the proposed split avoids dots in new
package names for readability and to avoid confusion with S3 method notation.

## Package Responsibilities

### `gflow`

`gflow` should contain durable infrastructure whose central object is a graph,
geodesic structure, gradient flow, or geometric decomposition of a point cloud.

Belongs in `gflow`:

- kNN/radius/adaptive-radius graph construction and repair;
- graph shortest paths, graph geodesic distances, graph spectra, and graph
  embeddings;
- graph/geometry diagnostics, isometry diagnostics, graph pruning, local
  geodesic pruning, and endpoint/core geometry;
- gradient-flow cells, basins, Morse-Smale complexes, extrema, trajectories,
  and related plotting/summary methods;
- PHATE/diffusion/geodesic geometry primitives when used as graph or geometry
  infrastructure;
- low-level graph utilities that multiple downstream packages can import.

Does not belong in `gflow` long-term:

- fitted-response smoothers whose main user-facing purpose is estimating
  `E[Y | X]` or predicting `y`;
- method families that only use graphs as one ingredient in a larger smoother
  or regression estimator;
- microbiome/application-specific analysis workflows unless they are general
  graph primitives.

### `geosmooth`

`geosmooth` should contain fitted-response methods that estimate conditional
expectations, smooth observed responses, or compare smoothing/regression
strategies over geometric supports.

Belongs in `geosmooth`:

- `kernel.local.polynomial.cv`;
- LPL-TF and SLPLiFT/S-LPL-TF operators, fitters, refitters, and predictors;
- MALPS;
- graph trend filtering and metric-graph low-pass smoothers when used as
  response smoothers rather than graph infrastructure;
- SSRHE response smoothers and Hessian-energy smoothers;
- harmonic extension/smoothing used as regression or imputation baselines;
- shared local chart, support, weighted least-squares, CV, and selector
  infrastructure needed by these methods.

`geosmooth` may import graph builders and graph distance functions from
`gflow`.  It should not duplicate durable graph-construction logic.

### `gflowx`

`gflowx` should collect assets that are valuable to preserve publicly but are
not appropriate for the stable CRAN-facing core.

Belongs in `gflowx`:

- retired or superseded conditional-expectation experiments;
- exploratory MALO/GFA method families unless actively maintained;
- highly application-specific microbiome, trajectory, concordance, and
  clustering workflows;
- visualization helpers tied to one-off analysis workflows;
- association-test pipelines that are not required by the core graph package.

`gflowx` does not need to be CRAN-submitted, but it should still use clean
R-package conventions so code remains installable and testable.

## Shared or Undecided Assets

Some assets should remain undecided until SA2/SA3 dependency analysis:

- generated `RcppExports` and native registration files;
- vendored third-party headers and licenses;
- basic kernels, local PCA chart helpers, kNN cache helpers, random/synthetic
  data helpers, plotting helpers, and p-value utilities;
- general grids, density/divergence helpers, low-level linear model solvers,
  and OpenMP diagnostics;
- generated `man/*.Rd` pages, which should follow their source files after the
  split rather than being classified independently.

For shared infrastructure, the default rule is:

1. if it is graph/geodesic infrastructure, keep it in `gflow`;
2. if it is response-model infrastructure used only by smoothers, move it to
   `geosmooth`;
3. if it is needed by both, keep one implementation in the lowest-level package
   that avoids cyclic dependencies;
4. if it is obsolete or application-specific, move it to `gflowx`.

## Migration Principles

- Do not move files until the SA1 inventory and SA2 dependency graph agree.
- Keep public API wrappers in `gflow` temporarily for moved CRAN-facing
  functions, with deprecation warnings only after downstream scripts are
  migrated.
- Move one pilot method first.  The recommended pilot is
  `kernel.local.polynomial.cv`, because it is important, relatively self
  contained, and exercises shared support-size, kernel, local-chart, and CV
  infrastructure.
- Do not split compiled code by copying large helper blocks.  SA3 should decide
  whether each C++ helper stays in `gflow`, moves to `geosmooth`, or is exposed
  through a small shared API.
- Keep generated files generated.  Do not hand-edit `NAMESPACE`, `man/*.Rd`, or
  Rcpp generated files during the extraction.

## SA0 Exit Criteria

SA0 is complete when:

- the three package roles are defined;
- belongs/does-not-belong rules are explicit;
- shared/undecided rules are recorded;
- the first pilot extraction candidate is identified;
- no source files have been moved.

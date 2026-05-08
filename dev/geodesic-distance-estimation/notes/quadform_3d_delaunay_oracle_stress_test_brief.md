# 3D Quadform Delaunay Oracle Stress-Test Brief

This brief covers background, source assets, and proposed tasks for stress-testing the current R/Qhull-backed Delaunay reference-geodesic oracle for 3D quadratic hypersurfaces.

## Background

The current gflow benchmark program is focused on the data geodesic geometric reconstruction problem:

\[
X=\{x_1,\ldots,x_n\}\subset \Gamma,
\qquad
D_G(i,j)\approx D_\Gamma(i,j).
\]

Here \(X\) is a finite sample from a known quadratic graph hypersurface \(\Gamma\), \(D_G\) is the shortest-path metric of a data-derived weighted graph on \(X\), and \(D_\Gamma\) is a numerical reference approximation to the true surface geodesic distance. The 2D benchmark currently uses regular parameter-domain grids. For parameter dimension 3, a regular tensor grid becomes too large and too axis-aligned, so we started a new reference-oracle path based on an approximate epsilon-net plus the Delaunay one-skeleton.

For a 3D parameter point \(x=(x_1,x_2,x_3)\), the quadratic graph hypersurface is

\[
\Gamma=\{(x,q(x)):x\in\Omega\subset\mathbb{R}^3\},
\]

with

\[
q(x)=\sum_{i=1}^{k}c_i x_i^2-\sum_{i=k+1}^{3}c_i x_i^2.
\]

The current experimental oracle:

1. Builds an approximate epsilon-net in the 3D parameter domain \(\Omega\).
2. Forcibly includes the sample points \(X\) as reference vertices.
3. Adds boundary candidates for ball/cube boundary coverage.
4. Computes a Delaunay tessellation in parameter space via `geometry::delaunayn()` (Qhull).
5. Extracts the Delaunay one-skeleton.
6. Optionally filters long Delaunay edges by parameter-space length.
7. Weights retained edges by exact quadratic-hypersurface segment length.
8. Computes shortest-path distances between the sample vertices.

The goal of this stress-test lane is to decide whether that oracle design is stable enough to become the high-dimensional reference geometry used in the next graph-construction benchmarks. It is not yet a request to port Delaunay to C++.

## Current Implementation Assets

Package code:

- `R/quadform_geodesics.R`
  - `quadform.edge.length()`
  - `quadform.edge.lengths()`
  - `quadform.delaunay.geodesic.distances()`
  - Internal helpers for 3D candidate sampling, boundary sampling, epsilon-net construction, Delaunay edge extraction, edge filtering, adjacency construction, and component checks.
- `src/quadform_grid_geodesics.cpp`
  - Generic C++ exact edge-length kernel used by `quadform.edge.lengths()`.
  - Existing 2D grid geodesic C++ code remains separate.
- `src/RcppExports.cpp`, `R/RcppExports.R`, `src/init.c`
  - Rcpp wrapper and native routine registration for `rcpp_quadform_edge_lengths`.
- `DESCRIPTION`
  - `geometry` is listed in `Suggests` for the current R/Qhull Delaunay backend.

Tests:

- `tests/testthat/test-quadform-geodesics.R`
  - Edge-length parity tests between scalar R and vectorized C++ implementations.
  - Small Delaunay oracle smoke test, skipped if `geometry` is unavailable.

Current report runner:

- `dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/run_quadform_3d_delaunay_reference_report.R`

Current generated report from the first convergence run:

- `dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/report/quadform_3d_delaunay_reference_report.html`

Current generated CSV summaries:

- `dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/results/delaunay_reference_run_summary.csv`
- `dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/results/delaunay_reference_error_summary.csv`

Related benchmark/report assets:

- `dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- `dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`
- `dev/geodesic-distance-estimation/notes/geodesic_isometry_diagnostics_report_brief.md`
- `dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R`

## Current First-Pass Findings

The first readiness report ran target reference sizes

\[
N_{\mathrm{ref}}\in\{2500,5000,10000,20000\}
\]

on three 3D cases:

- ball domain, \(k=3\), elliptic paraboloid-type hypersurface;
- ball domain, \(k=1\), mixed-sign quadratic hypersurface;
- cube domain, \(k=2\), mixed-sign quadratic hypersurface.

Approximate results from the first run:

- 20k-target cases produced about 29k actual reference vertices and about 215k retained Delaunay edges.
- Runtime was roughly 17 seconds per 20k-target case on the local machine.
- Relative RMS error versus the 20k reference was around \(0.016\) to \(0.018\) at the 10k-target level.

These results are promising, but they are not enough to declare the oracle ready. We need a more systematic stress-test of stability, boundary behavior, filtering behavior, curvature sensitivity, and failure modes.

## Stress-Test Questions

The stress-test should answer the following questions.

### 1. Epsilon-Net Stability

Does the approximate epsilon-net produce enough quasi-uniform coverage for the oracle to stabilize as \(N_{\mathrm{ref}}\) grows?

Diagnostics should include:

- target reference size versus actual reference vertex count;
- estimated epsilon;
- nearest-neighbor distance distribution among non-forced reference vertices;
- sample-to-reference inclusion check;
- boundary candidate retention rate;
- runtime and memory where feasible.

### 2. Boundary Coverage

Does the oracle behave differently near the boundary of a ball or cube domain?

Stress cases should include:

- samples concentrated in the interior;
- samples concentrated near the boundary;
- mixed interior/boundary samples;
- ball and cube domains separately.

Useful diagnostics:

- pairwise error stratified by whether endpoints are near boundary;
- shortest paths that visibly hug or cut across boundary regions;
- convergence of boundary-heavy samples versus interior samples.

### 3. Delaunay Edge Filtering

The current oracle filters Delaunay edges by parameter-space edge length:

\[
\|u-v\|_2 \le \lambda\epsilon,
\]

where \(\lambda=\texttt{edge.length.factor}\). If filtering disconnects the graph, the implementation progressively relaxes \(\lambda\), falling back to the unfiltered one-skeleton if necessary.

The stress-test should examine:

- \(\lambda\in\{2, 2.5, 3, 4, 6, \infty\}\);
- number/proportion of edges removed;
- whether filtering disconnects the graph before relaxation;
- final `filter_factor_used`;
- impact on geodesic-distance convergence;
- evidence of long-edge shortcuts when filtering is too permissive.

### 4. Curvature And Index Sensitivity

The oracle should be tested across quadratic coefficients and index values:

\[
(c_1,c_2,c_3)\in\{(1,1,1),(1,1,2),(1,2,4),(1,4,4)\}
\]

and

\[
k\in\{0,1,2,3\}.
\]

The goal is not to run a huge benchmark immediately, but to see whether convergence or edge filtering becomes unstable when curvature is anisotropic or the signature changes.

### 5. Reference-Size Convergence

Use a finest available reference as a numerical target and compare coarser references against it:

\[
N_{\mathrm{ref}}\in\{2500,5000,10000,20000\}
\]

with an optional larger sentinel, for example \(40000\), on one or two cases if runtime permits.

Compare distance matrices using the existing diagnostics:

- `rel_rms_error`;
- `rel_geodesic_stress`;
- `signed_bias`;
- `shortcut_fraction`;
- `q50_rel_abs_residual`;
- `q90_rel_abs_residual`;
- `q95_rel_abs_residual`;
- distance correlation.

### 6. Failure Modes

The stress-test should explicitly report failures or warning-worthy behavior:

- Delaunay failures from Qhull;
- disconnected filtered graphs;
- excessive relaxation to `Inf`;
- strong shortcut bias;
- large boundary-specific error;
- non-monotone convergence;
- runtime or memory blowups.

## Recommended Output

Create a new stress-test directory:

```text
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/
```

Recommended assets:

- `run_quadform_3d_delaunay_oracle_stress_test.R`
- `report/quadform_3d_delaunay_oracle_stress_test_report.html`
- `results/oracle_stress_run_summary.csv`
- `results/oracle_stress_error_summary.csv`
- `results/oracle_stress_edge_filter_summary.csv`
- `results/oracle_stress_boundary_summary.csv`
- optional cached RDS payloads under `runs/` or `cache/` if needed.

The report should read like a research-facing validation report, not only a table dump. Figures should precede detailed tables.

Recommended report sections:

1. Motivation and oracle definition.
2. Experimental design.
3. Reference-size convergence.
4. Edge-filtering sensitivity.
5. Boundary stress tests.
6. Curvature and index sensitivity.
7. Failure modes and warnings.
8. Recommendation: whether the current R/Qhull oracle is stable enough to serve as behavioral truth for a future C++ Delaunay backend.

## Testing And QA Expectations

Before modifying package code, inspect the existing implementation and determine whether the stress-test can be built as a dev-report-only script. Prefer not to change package code unless the stress-test reveals a correctness bug or a small missing diagnostic that belongs in the package API.

If package code is changed:

```bash
make document
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-quadform-geodesics.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
```

If only dev-report scripts are added:

```bash
Rscript dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R --mode=smoke
```

Then run the fuller validation only after the smoke report is credible.

## Important Constraints

- Do not port Delaunay to C++ in this lane.
- Treat `geometry::delaunayn()` as the current behavioral implementation under test.
- Do not overwrite the existing first-pass report unless explicitly asked.
- Keep generated heavy report outputs local/ignored unless asked to track them.
- Preserve existing 2D grid-oracle behavior.
- If failures are found, document them plainly and propose fixes rather than hiding them behind new defaults.

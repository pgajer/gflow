# Prompt: Stress-Test the 3D Quadform Delaunay Reference Oracle

You are working in the `gflow` repository:

```text
/Users/pgajer/current_projects/gflow
```

Your task is to stress-test the current R/Qhull-backed Delaunay reference-geodesic oracle for 3D quadratic hypersurfaces.

First read this brief carefully:

```text
dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_brief.md
```

Then inspect the relevant implementation and test assets:

```text
R/quadform_geodesics.R
src/quadform_grid_geodesics.cpp
tests/testthat/test-quadform-geodesics.R
dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/run_quadform_3d_delaunay_reference_report.R
dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/report/quadform_3d_delaunay_reference_report.html
dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/results/delaunay_reference_run_summary.csv
dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/results/delaunay_reference_error_summary.csv
```

## Goal

Design and execute a stress-test report that decides whether the current 3D oracle is stable enough to serve as behavioral truth for a future C++ Delaunay backend.

The oracle under test is:

```r
quadform.delaunay.geodesic.distances()
```

It currently uses:

- approximate epsilon-net reference vertices;
- forced inclusion of sample points;
- boundary candidates;
- `geometry::delaunayn()` as the Qhull-backed Delaunay implementation;
- exact C++ quadratic edge lengths via `quadform.edge.lengths()`;
- optional long-edge filtering controlled by `edge.length.factor`.

Do not port Delaunay to C++ in this task. This lane is about stress-testing the current R/Qhull oracle.

## Required Workflow

1. Read the brief and inspect the implementation.
2. Generate a concrete action plan before implementation.
3. Share the plan with the user for approval before executing any large run or package-code change.
4. After approval, implement the stress-test script/report.
5. Run a smoke test first.
6. Inspect the smoke report and only then run the fuller stress test.
7. Summarize findings, including failures and recommended defaults.

## Recommended Output Directory

Create a new directory:

```text
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/
```

Recommended output files:

```text
run_quadform_3d_delaunay_oracle_stress_test.R
report/quadform_3d_delaunay_oracle_stress_test_report.html
results/oracle_stress_run_summary.csv
results/oracle_stress_error_summary.csv
results/oracle_stress_edge_filter_summary.csv
results/oracle_stress_boundary_summary.csv
```

Keep heavy/generated outputs local unless the user explicitly asks to track them.

## Stress-Test Dimensions

The stress-test should cover:

- reference-size convergence:
  \[
  N_{\mathrm{ref}}\in\{2500,5000,10000,20000\}
  \]
  with optional \(40000\) sentinel on a small number of cases if runtime permits;
- ball and cube domains;
- interior, boundary-heavy, and mixed sample placement;
- index values:
  \[
  k\in\{0,1,2,3\};
  \]
- curvature coefficients:
  \[
  (c_1,c_2,c_3)\in\{(1,1,1),(1,1,2),(1,2,4),(1,4,4)\};
  \]
- edge filter values:
  \[
  \texttt{edge.length.factor}\in\{2,2.5,3,4,6,\infty\}.
  \]

Do not blindly run the full Cartesian product if it is too large. Use a staged design: smoke, focused full grid, then optional sentinel cases.

## Required Diagnostics

At minimum, compute and report:

- target versus actual reference vertex count;
- epsilon;
- number of Delaunay edges before filtering;
- number/proportion of retained edges after filtering;
- final `filter_factor_used`;
- connectedness and any filter relaxation;
- runtime;
- scale-calibrated relative RMS error versus the finest reference;
- `rel_geodesic_stress`;
- `signed_bias`;
- `shortcut_fraction`;
- `q50_rel_abs_residual`;
- `q90_rel_abs_residual`;
- `q95_rel_abs_residual`;
- boundary-stratified errors where possible.

Use existing helper diagnostics such as:

```r
summarize.isometry.deviation()
```

## Report Expectations

The HTML report should read like a research validation report. It should include explanatory prose, formulas where helpful, figures before detailed tables, and a clear final recommendation.

Required sections:

1. Motivation and oracle definition.
2. Experimental design.
3. Reference-size convergence.
4. Edge-filtering sensitivity.
5. Boundary stress tests.
6. Curvature and index sensitivity.
7. Failure modes and warnings.
8. Recommendation.

The recommendation should answer:

- Is the current R/Qhull oracle stable enough to use as behavioral truth?
- What default `edge.length.factor` is recommended?
- Which cases remain concerning?
- What should be fixed before a C++ Delaunay backend is ported?

## QA

If you only add dev-report scripts, run at least:

```bash
Rscript dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R --mode=smoke
```

If you modify package code, also run:

```bash
make document
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-quadform-geodesics.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
```

Report all warnings, failures, and skipped checks explicitly.

# Follow-Up Prompt: Finalize 3D Quadform Delaunay Oracle Stress Test

You are working in the `gflow` repository:

```text
/Users/pgajer/current_projects/gflow
```

Read the original stress-test brief and prompt first:

```text
dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_brief.md
dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_agent_prompt.md
```

Then inspect your current stress-test runner and generated outputs:

```text
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/report/quadform_3d_delaunay_oracle_stress_test_report.html
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/results/
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test-sentinel/report/quadform_3d_delaunay_oracle_stress_test_report.html
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test-sentinel/results/
```

## Audit Feedback

Your controlled rerun fixed the main methodological issue in the first edge-filter sweep. The original sweep was confounded because the oracle seed depended on `edge.length.factor`; changing the seed also changed the epsilon-net/Delaunay graph. The updated sweep now uses the same sample, same oracle seed, and same target reference size across edge factors, so it isolates filtering.

The sanity check now behaves correctly:

- when `edge.length.factor = 6` retained all Delaunay edges, it matched `Inf` exactly;
- `Inf` correctly represents the no-filter path;
- no oracle runs failed;
- no warnings, disconnected graphs, or filter relaxations were observed.

The updated edge-filter conclusion is more convincing:

- `2` is too aggressive;
- `2.5` is still visibly different from `Inf`;
- `3` is usually small but still non-negligible in the controlled sweep;
- `4` is nearly indistinguishable from `Inf`;
- `6` and `Inf` were identical in the controlled sweep.

Therefore the current empirical evidence favors `edge.length.factor = 4` as the provisional reference-oracle default if filtering remains enabled. Factor `3` should be described as a locality-preserving alternative, not as the default for generating behavioral-truth reference distances.

The 40k sentinel results are reassuring, but the sentinel mode appears to have used factor `3`. Those results show that the hard cases converge acceptably, but they do not directly validate factor `4` as the hard-case reference policy. One final sentinel pass is needed.

## Your Remaining Task

Finish the empirical validation and recommendation layer. Do not start the C++ Delaunay port. Do not add new package-level diagnostic APIs unless explicitly asked later.

## Required Next Step 1: Final Factor-4 Sentinel Pass

Run a final sentinel validation on the two hard cases:

```text
ball_mixed_k1_c124
cube_mixed_k0_c144
```

Use:

```text
N_ref = 10000, 20000, 40000
edge.length.factor = 4
```

Also include the `Inf` no-filter path at least at `N_ref = 40000`, and preferably at both `20000` and `40000` if runtime is acceptable.

The purpose is to answer:

1. Does factor `4` remain effectively indistinguishable from `Inf` in the high-anisotropy/mixed hard cases?
2. Do boundary-pair q95 residuals continue to improve with resolution?
3. Are there any warnings, Qhull failures, disconnections, or filter relaxations?
4. Should factor `4` become the provisional default for reference-truth generation?

Recommended output directory:

```text
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test-factor4-sentinel/
```

It should contain:

```text
report/quadform_3d_delaunay_oracle_stress_test_report.html
results/oracle_stress_run_summary.csv
results/oracle_stress_error_summary.csv
results/oracle_stress_boundary_summary.csv
```

If you need to modify the runner, keep the modification small and clear. The runner should support this sentinel mode without disrupting the already generated controlled full report.

## Required Next Step 2: Update Recommendation Text

After the final factor-4 sentinel pass, update the stress-test report text or create a short follow-up markdown note summarizing:

- whether factor `4` is confirmed as the provisional reference-oracle default;
- whether factor `Inf` should remain available as the no-filter behavioral baseline;
- how factor `3` should be described;
- what tolerance ranges are reasonable for high-anisotropy/mixed boundary cases;
- which cases remain the hardest.

The language should be careful:

- “provisional behavioral truth” is appropriate;
- “mathematically exact truth” is not appropriate;
- “C++ backend behavioral harness” is appropriate once the final sentinel pass supports it.

## Required Next Step 3: Prepare a Regression-Fixture Table

Create a small markdown or CSV asset listing proposed regression fixtures for the future C++ Delaunay backend. This is not the C++ implementation; it is the fixture design.

Suggested file:

```text
dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_cpp_regression_fixtures.md
```

The fixture table should include:

- fixture id;
- domain shape;
- placement;
- index \(k\);
- coefficients;
- sample seed;
- sample size;
- `N_ref`;
- `edge.length.factor`;
- reason for inclusion;
- expected comparison target.

At minimum include:

1. Easy ball interior case:
   - `ball_interior_k3_c111`
   - factor `4`
   - optionally `Inf`

2. High-anisotropy ball mixed case:
   - `ball_mixed_k1_c124`
   - factor `4`
   - `Inf`

3. High-anisotropy cube mixed case:
   - `cube_mixed_k0_c144`
   - factor `4`
   - `Inf`

The expected comparison target should be the R/Qhull oracle output, not a hand-derived exact distance.

## QA Expectations

If you only modify the dev stress-test runner and notes, run the relevant stress-test commands and report exact outputs.

If you modify package code, run:

```bash
make document
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-quadform-geodesics.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
```

If no package code is modified, still run:

```bash
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-quadform-geodesics.R")'
```

## Final Response Requirements

Report:

- exact commands run;
- whether all factor-4 sentinel runs succeeded;
- factor `4` versus `Inf` comparison in the hard cases;
- boundary-pair tail behavior;
- final recommendation for `edge.length.factor`;
- path to the updated/follow-up report;
- path to the regression-fixture table;
- any residual concerns before C++ porting.

Do not commit or push unless explicitly asked.

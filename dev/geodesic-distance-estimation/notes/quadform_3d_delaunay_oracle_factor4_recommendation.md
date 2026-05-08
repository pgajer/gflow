# 3D Quadform Delaunay Oracle: Factor-4 Follow-Up Recommendation

Generated after the controlled edge-filter rerun and the final factor-4
sentinel pass. This note updates the recommendation for the current
R/Qhull-backed oracle:

```r
quadform.delaunay.geodesic.distances()
```

The oracle remains a numerical reference oracle. It is appropriate to call it
provisional behavioral truth for the C++ backend harness, not mathematically
exact truth.

## Final Factor-4 Sentinel Design

The final sentinel pass used the two hardest mixed/anisotropic cases:

- `ball_mixed_k1_c124`
- `cube_mixed_k0_c144`

It ran:

- `edge.length.factor = 4` at `N_ref = 10000, 20000, 40000`;
- `edge.length.factor = Inf` at `N_ref = 20000, 40000`.

Output directory:

```text
dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test-factor4-sentinel/
```

All 10 oracle runs succeeded. There were no Qhull warnings, disconnected
graphs, or filter relaxations.

## Factor 4 Versus Inf

At matching reference sizes, factor `4` remained effectively indistinguishable
from the no-filter path.

| case_id | N_ref | rel_rms_error vs Inf | q95_rel_abs_residual vs Inf |
|---|---:|---:|---:|
| `ball_mixed_k1_c124` | 20000 | 0.001042 | 0.000126 |
| `ball_mixed_k1_c124` | 40000 | 0.000000 | 0.000000 |
| `cube_mixed_k0_c144` | 20000 | 0.002865 | 0.000576 |
| `cube_mixed_k0_c144` | 40000 | 0.001644 | 0.000236 |

The cube case remains the harder filter case, but the residuals are small
relative to the discretization error of the finite reference itself.

## Resolution And Boundary Tails

Against the 40k factor-4 reference:

| case_id | N_ref | rel_rms_error | q95_rel_abs_residual |
|---|---:|---:|---:|
| `ball_mixed_k1_c124` | 10000 | 0.022519 | 0.074470 |
| `ball_mixed_k1_c124` | 20000 | 0.019531 | 0.058881 |
| `cube_mixed_k0_c144` | 10000 | 0.021086 | 0.075865 |
| `cube_mixed_k0_c144` | 20000 | 0.018005 | 0.064376 |

Boundary-pair q95 residuals improved with resolution. The largest boundary or
interior stratum q95 residuals in the factor-4 sentinel were:

- cube mixed, 10k: about `0.102` in the interior-interior stratum and `0.090`
  in the boundary-boundary stratum;
- cube mixed, 20k: about `0.086` in the interior-interior stratum and `0.070`
  in the boundary-boundary stratum;
- ball mixed, 20k: about `0.068` in both boundary-boundary and
  interior-interior strata.

Reasonable tolerance bands for these hard high-anisotropy/mixed cases are:

- 20k versus 40k reference-size convergence: expect relative RMS around
  `0.018-0.020`, with q95 residual around `0.06-0.07`;
- boundary-stratified 20k versus 40k: allow q95 residual up to about `0.09`
  for the hardest cube mixed strata;
- factor `4` versus `Inf` at matched resolution: expect relative RMS below
  `0.003` and q95 residual below `0.001` in these sentinel cases.

These are tolerance ranges for validating the reference-oracle lane, not
scientific claims about exact geodesic error.

## Recommendation

Use `edge.length.factor = 4` as the provisional reference-oracle default when
generating behavioral-truth distances for the future C++ Delaunay backend.

Keep `edge.length.factor = Inf` available as the no-filter behavioral baseline.
It is essential for debugging and for confirming that finite filter factors are
not changing sample geodesic distances in sensitive cases.

Describe `edge.length.factor = 3` as a locality-preserving alternative. In the
controlled sweep it was usually small, but still detectably different from
`Inf`, especially in the high-anisotropy mixed cases. It should not be the
default for behavioral-truth reference generation.

The hardest cases remain:

- `cube_mixed_k0_c144`, especially boundary and near-boundary pair strata;
- `ball_mixed_k1_c124`, mainly as the representative high-anisotropy ball
  mixed case.

With the final sentinel pass, the R/Qhull oracle is credible as the first C++
backend behavioral harness, provided the future comparison tolerances are
explicit and the hard mixed/boundary cases are included in regression fixtures.

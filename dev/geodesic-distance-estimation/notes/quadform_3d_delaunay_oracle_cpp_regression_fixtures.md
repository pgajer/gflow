# Proposed C++ Regression Fixtures For The 3D Quadform Delaunay Oracle

These fixtures are for a future C++ Delaunay backend. They do not imply a C++
implementation in the current lane. The expected comparison target is the
current R/Qhull oracle output from:

```r
quadform.delaunay.geodesic.distances()
```

The target is behavioral parity with the R/Qhull oracle, not a hand-derived
closed-form geodesic distance.

Recommended default fixture policy:

- use `edge.length.factor = 4` for reference-truth generation;
- include `edge.length.factor = Inf` no-filter fixtures as debugging and
  behavioral-baseline controls;
- use fixed sample seeds, oracle seeds, sample size, reference size, candidate
  multiplier, and boundary fraction when comparing R/Qhull and C++ outputs.

Unless a fixture says otherwise, use:

- `domain.radius = 1`;
- `sample_n = 36`;
- `candidate.multiplier = 4`;
- `boundary.fraction = 0.25`;
- `oracle_seed = sample_seed + N_ref`;
- expected target: R/Qhull oracle distance matrix and graph diagnostics.

| fixture_id | domain_shape | placement | index_k | coefficients | sample_seed | sample_n | N_ref | edge.length.factor | reason for inclusion | expected comparison target |
|---|---|---|---:|---|---:|---:|---:|---:|---|---|
| `cpp_ball_interior_k3_c111_f4_10k` | ball | interior | 3 | `(1,1,1)` | 101 | 36 | 10000 | 4 | Easy ball interior baseline; low anisotropy; checks standard filtered path. | R/Qhull oracle distances and diagnostics at same seed/settings. |
| `cpp_ball_interior_k3_c111_inf_10k` | ball | interior | 3 | `(1,1,1)` | 101 | 36 | 10000 | `Inf` | No-filter control for the easy baseline. | R/Qhull oracle no-filter distances and diagnostics. |
| `cpp_ball_mixed_k1_c124_f4_20k` | ball | mixed | 1 | `(1,2,4)` | 103 | 36 | 20000 | 4 | High-anisotropy ball mixed case; representative hard but stable reference-truth fixture. | R/Qhull oracle distances and diagnostics at same seed/settings. |
| `cpp_ball_mixed_k1_c124_inf_20k` | ball | mixed | 1 | `(1,2,4)` | 103 | 36 | 20000 | `Inf` | No-filter baseline for high-anisotropy ball mixed case; verifies filter impact remains small. | R/Qhull oracle no-filter distances and diagnostics. |
| `cpp_ball_mixed_k1_c124_f4_40k_sentinel` | ball | mixed | 1 | `(1,2,4)` | 103 | 36 | 40000 | 4 | Sentinel-resolution target for the high-anisotropy ball mixed case. | R/Qhull oracle distances and diagnostics at same seed/settings. |
| `cpp_ball_mixed_k1_c124_inf_40k_sentinel` | ball | mixed | 1 | `(1,2,4)` | 103 | 36 | 40000 | `Inf` | Sentinel no-filter baseline; factor 4 matched this exactly in the final pass. | R/Qhull oracle no-filter distances and diagnostics. |
| `cpp_cube_mixed_k0_c144_f4_20k` | cube | mixed | 0 | `(1,4,4)` | 106 | 36 | 20000 | 4 | Hardest mixed/anisotropic cube case; stresses cube boundary behavior and filtering. | R/Qhull oracle distances and diagnostics at same seed/settings. |
| `cpp_cube_mixed_k0_c144_inf_20k` | cube | mixed | 0 | `(1,4,4)` | 106 | 36 | 20000 | `Inf` | No-filter baseline for the hardest cube case. | R/Qhull oracle no-filter distances and diagnostics. |
| `cpp_cube_mixed_k0_c144_f4_40k_sentinel` | cube | mixed | 0 | `(1,4,4)` | 106 | 36 | 40000 | 4 | Sentinel-resolution target for the hardest cube mixed case; includes boundary-pair tails. | R/Qhull oracle distances and diagnostics at same seed/settings. |
| `cpp_cube_mixed_k0_c144_inf_40k_sentinel` | cube | mixed | 0 | `(1,4,4)` | 106 | 36 | 40000 | `Inf` | Sentinel no-filter baseline; checks whether finite filtering changes hard-case distances. | R/Qhull oracle no-filter distances and diagnostics. |

## Suggested Comparison Metrics

For same-settings R/Qhull versus C++ comparisons, compute:

- `summarize.isometry.deviation(D_cpp, D_r_qhull, scale = TRUE)`;
- max absolute pairwise distance difference;
- max relative pairwise distance difference for nonzero R/Qhull distances;
- graph diagnostics parity or near-parity, including actual reference vertex
  count, epsilon, unfiltered edge count, retained edge count, retained edge
  fraction, requested/used filter factor, component count, and relaxation state.

For exact same Delaunay triangulations, distance matrices should be nearly
identical up to floating-point tolerance. If the C++ triangulation has valid
but non-identical tie-breaking, evaluate with the geodesic-deviation summaries
and inspect the hard cube mixed fixtures first.

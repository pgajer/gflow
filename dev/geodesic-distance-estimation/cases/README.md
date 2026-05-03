# Benchmark Cases

Place source scripts for individual synthetic benchmark cases here. Each case
should define:

- a deterministic random seed,
- latent coordinates or known intrinsic geodesic distances,
- observed feature matrix construction,
- optional response field for downstream regression checks,
- expected failure modes for graph construction.

Generated outputs from these cases should go to `../reports/`, `../figures/`,
or `../cache/`, not this directory.

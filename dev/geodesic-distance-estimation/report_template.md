# Geodesic Distance Estimation HTML Report Template

This is a living template for benchmark reports in the geodesic distance estimation
development project. It is guidance, not a rigid specification. Each benchmark
report may depart from this structure when the example demands a different
presentation. Update this file when repeated departures reveal a better pattern.

The purpose of each report is to understand one benchmark example deeply: why it
was chosen, what geodesic-distance-estimation challenge it poses, how candidate
graph constructions and truth definitions behave, and what this example teaches
us about the next iteration of the method family.

## Report-Level Requirements

Each report should be reproducible, visual, and interpretive.

- Use deterministic random seeds.
- Show explicit parameter values for data generation and graph construction.
- Include commands or code blocks sufficient to rerun the case.
- Prioritize figures and tables that reveal graph behavior, not just aggregate
  scores.
- Compare relevant candidate methods: CMST/MST, iKNN, mKNN, pruning variants,
  component repair policies, and oracle graphs when available.
- Interpret results in text, including failures and ambiguous behavior.
- Record open questions raised by the example.

Generated HTML reports should be written under `dev/geodesic-distance-estimation/reports/`.
Generated figures should be written under `dev/geodesic-distance-estimation/figures/`.
Temporary data or model objects should be written under `dev/geodesic-distance-estimation/cache/`.
These output directories are ignored by git.

## Suggested Section Order

### 1. Title and Metadata

Include:

- benchmark case name,
- generation timestamp,
- git commit or working tree note if available,
- random seed,
- sample size,
- main graph parameters,
- preprocessing mode.

### 2. Why This Example

Explain the scientific and algorithmic reason for the example.

Prompts:

- What latent geometry is being simulated?
- Why is this geometry important for geodesic distance estimation?
- What does a successful estimator need to preserve?
- What failure modes are expected?
- How does this case relate to downstream gflow use cases?

### 3. Expected Challenges

List the anticipated graph-construction issues before showing results.

Possible issues:

- nonuniform sampling density,
- ambient-space shortcuts,
- disconnected local graphs,
- forced MST bridges,
- noisy off-manifold points,
- false cross-branch edges,
- periodic boundary continuity,
- compositional zeros or library-size variation,
- high-dimensional nuisance variation,
- sensitivity to `q.thld` or `k`.

### 4. Reproducible Commands

Show the exact commands used to generate the data and graphs.

Include:

- `set.seed(...)`,
- data-generating function call,
- response-generating function call if applicable,
- preprocessing steps,
- CMST construction call,
- iKNN construction calls,
- mKNN construction calls when applicable,
- oracle graph construction call when applicable,
- parameter sweep values.

The goal is not to print every helper implementation, but to make the case
reconstructible from source files in `dev/geodesic-distance-estimation/`.

### 5. Data-Generating Mechanism

Describe the synthetic truth.

Include:

- latent coordinate definition,
- observed feature matrix definition,
- noise model,
- response field if present,
- known intrinsic geodesic distance,
- known forbidden shortcuts or true bridges.

For compositional examples, include:

- latent ecological coordinates,
- taxon abundance functions,
- count model,
- library-size distribution,
- zero-inflation mechanism,
- preprocessing used before graph construction.

### 6. Data Geometry Figures

Recommended figures:

- observed data colored by latent coordinate,
- observed data colored by sampling density or region,
- observed data colored by true response if present,
- oracle path, circle, tree, or known adjacency overlaid on geometry,
- high-dimensional projection used for visualization if the data are not 2D.

For circular examples:

- angle histogram,
- geometry colored by angle,
- geometry colored by periodic response,
- boundary region near `0 / 2pi` highlighted.

### 7. Graph Constructions

Describe each graph included in the comparison.

Recommended graph set:

- **MST-only**: the connected sparse backbone.
- **CMST**: the method under development, with `q.thld` and
  `cmst_distance_threshold` reported.
- **iKNN**: local graph competitor, ideally over a `k` sweep.
- **Oracle graph**: weighted path, circle, tree, or known latent graph when
  available.

For each graph, report:

- edge count,
- component count,
- mean, min, and max degree,
- edge-length summary,
- long-edge or suspicious-edge counts,
- parameter values.

### 8. Graph Visualization Figures

Recommended figures:

- MST-only edge overlay,
- CMST edge overlay,
- iKNN edge overlay for selected `k`,
- oracle edge overlay,
- CMST edges colored by edge length,
- vertices colored by degree,
- long MST bridges highlighted,
- suspected shortcut edges highlighted.

When a parameter sweep is included, show:

- CMST graph overlays for representative `q.thld` values,
- edge count versus `q.thld`,
- false-shortcut rate versus `q.thld`,
- geodesic recovery metric versus `q.thld`.

### 9. Truth Comparison

Compare graph shortest-path distances to known latent geodesic distances.

Recommended tables:

- distance correlation by graph,
- median and upper-quantile distortion by graph,
- false-shortcut rate by graph,
- neighborhood preservation by graph,
- connectedness and graph complexity by graph.

Recommended figures:

- scatterplot of graph distances versus latent geodesic distances,
- residual/distortion versus latent distance,
- heatmap of graph distance matrix ordered by latent coordinate,
- heatmap of latent geodesic distance matrix,
- difference heatmap between graph and latent distances.

For circular examples:

- boundary-distance errors near `0 / 2pi`,
- chord-shortcut diagnostics,
- graph distance versus circular arc distance.

### 10. Downstream Response Check

Include this section when the example has a known response field.

Recommended analyses:

- graph regression or smoothing on CMST,
- same analysis on iKNN,
- same analysis on oracle graph,
- fitted response versus true response,
- residuals versus latent coordinate,
- boundary RMSE for periodic examples,
- RMSE and correlation summary table.

This section should not replace geodesic-distance evaluation. It is a downstream
sanity check: does graph geometry support the kind of analysis gflow wants to
perform?

### 11. Parameter Stability

Summarize how results change across key parameters.

For CMST:

- `q.thld`,
- PCA dimensionality or variance threshold,
- preprocessing mode,
- optional sample size and noise sweeps.

For iKNN:

- `k`,
- graph connectivity,
- edge count,
- recovery metrics.

The central question is the stable parameter region: not only which value wins,
but which range behaves acceptably.

### 12. Interpretation

Write a direct, case-specific interpretation.

Prompts:

- Did CMST recover the intended geodesic metric?
- Where did it outperform MST-only or iKNN?
- Where did it fail or behave ambiguously?
- Did the actual distance threshold make sense?
- Were long MST edges legitimate bridges or warning signs?
- Did completion add useful local edges or harmful shortcuts?
- What parameter range looks defensible?
- What does this example teach us about the algorithm?

### 13. Open Issues

Record unresolved questions raised by the example.

Examples:

- Should long MST bridges be flagged or pruned after connectivity is achieved?
- Should completion use a global threshold or local/adaptive thresholds?
- Should edge weights be raw Euclidean distances or transformed weights?
- Should duplicate points and zero-distance edges receive special treatment?
- Should compositional data use a domain-specific distance before CMST?
- What diagnostics should become package-level outputs?

### 14. Next Benchmark

End with a short note explaining what the next benchmark should test, based on
what this example revealed.

## Figure Checklist

Use this checklist as a starting point:

- observed geometry colored by latent coordinate,
- observed geometry colored by response truth or density,
- MST-only overlay,
- CMST overlay,
- selected iKNN overlay,
- oracle graph overlay,
- degree map,
- edge-length distribution,
- graph distance versus truth scatterplot,
- distance distortion heatmap,
- parameter sweep summary,
- response fit and residual plot if applicable.

Not every report needs every figure. Prefer fewer figures with clear captions
over many repetitive panels.

## Table Checklist

Useful tables include:

- data-generation parameters,
- graph-construction parameters,
- graph structural summary,
- edge-length summary,
- geodesic recovery summary,
- false-shortcut summary,
- parameter sweep summary,
- response recovery summary.

## Style Notes

- Write captions that explain what the reader should notice.
- Keep code blocks reproducible but not overwhelming.
- Put interpretation next to the relevant figure or table when possible.
- Name uncertainty explicitly.
- Preserve failures; they are often the most useful part of the benchmark.
- Revise this template as the benchmark ladder teaches us what matters.

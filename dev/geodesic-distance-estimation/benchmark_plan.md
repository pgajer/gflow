# Geodesic Distance Estimation Benchmark Plan

## Project Objective

The objective of this project is to develop and evaluate methods for estimating
useful intrinsic geodesic distances from finite, noisy, unevenly sampled, and
eventually high-dimensional biological data.

The project began with the Minimal Spanning Tree Completion Graph (CMST), as
constructed by `create.cmst.graph()`. CMST remains one important candidate
because it starts with an MST to keep all samples connected, then adds local
Euclidean edges using a data-adaptive threshold. The broader scientific question
is now method-agnostic: when does a graph construction, oracle truth definition,
or pruning/repair strategy produce shortest-path distances that behave like the
intended intrinsic geometry?

A good geodesic distance estimator should preserve local neighborhoods, avoid
unsupported ambient-space shortcuts, bridge sampling variability only when
appropriate, handle disconnected local-graph components in a principled way,
and support downstream tasks such as trajectory reconstruction, basin detection,
conditional expectation estimation, and graph-based regression.

The benchmark ladder should begin with simple non-symmetric synthetic examples
where the correct geodesic structure is known exactly. It should then move
quickly toward data shapes that resemble real high-dimensional biological data:
nonuniform sampling, branching, loops, gaps, high-dimensional nuisance
variation, compositional sparsity, count noise, batch effects, and community
state transitions. For noisy examples, the project should explicitly examine
which truth metric is appropriate, including latent arc/path distances and
density-weighted oracle distances on known generative manifolds or tubes.

## Core Evaluation Metrics

- **Latent geodesic recovery**: Spearman and Pearson correlation between graph
  shortest-path distances and known latent geodesic distances.
- **Distance distortion**: relative error of graph distances after optimal
  global scaling, summarized by median, upper quantiles, and worst cases.
- **False-shortcut rate**: fraction of edges or shortest paths that cross known
  gaps, branches, holes, or non-adjacent latent regions.
- **Connectivity**: component count, largest component fraction, and whether
  the MST backbone prevents disconnected failures seen in local kNN graphs.
- **Local neighborhood preservation**: overlap between graph neighborhoods and
  latent-neighborhood sets at multiple radii.
- **Graph complexity**: edge count, degree distribution, high-degree hubs, and
  sensitivity to `q.thld`.
- **Downstream stability**: stability of graph regression, basin calls,
  endpoint calls, and trajectory summaries across seeds and parameter sweeps.
- **Boundary behavior**: for circular or periodic examples, error near the
  `0 / 2pi` boundary and preservation of periodic continuity.

## Benchmark Ladder

### 1. Uneven Noisy Line

Generate a one-dimensional curve with asymmetric density: one dense end, one
sparse end, and mild off-manifold noise. The latent geodesic distance is the
absolute difference in the one-dimensional coordinate.

Purpose: establish the simplest nonuniform-density case and test whether a
global MST-edge quantile over-connects dense regions or under-connects sparse
ones.

Expected comparisons: CMST, iKNN, oracle weighted path graph.

### 2. Crescent or Open Arc

Generate a curved open arc with nonuniform angular sampling and radial noise.
The endpoints may be close in ambient space for tighter crescents.

Purpose: test whether CMST follows curved geodesics rather than adding
ambient-space chords across the arc.

Expected failure mode: high `q.thld` may introduce shortcut edges across the
inside of the crescent.

### 3. Two Nearby Non-Touching Arcs

Generate two arcs separated by a small ambient gap but distinct latent
components or branches.

Purpose: test false bridge formation. This case forces a distinction between
sampling gaps within one object and true separation between objects.

Expected failure mode: CMST must connect all points through the MST, so the
long forced bridge should be identifiable as a large MST edge. Completion edges
should not add many additional cross-gap shortcuts.

### 4. Noisy Circle: Periodic Signal Recovery

Use `generate.circle.data()` with random angles and radial noise, then define a
periodic response with `circular.synthetic.mixture.of.gaussians()`. This is the
case represented in
`tests/manual/reports/rdgraph-regression-correctness/rdgraph_regression_correctness.html`
as `circular_gaussian_mixture_recovers_periodic_signal`.

Purpose: test recovery of circular geodesic distance and preservation of
periodic continuity at the `0 / 2pi` boundary. This should be a central CMST
benchmark because it exposes topology, connectedness, local edge completion,
and false chord shortcuts in one small example.

Expected comparisons: CMST, iKNN across `k`, oracle weighted circle graph.

Key metrics: circular distance correlation, boundary RMSE, chord-shortcut rate,
edge count, component count, and regression RMSE against known periodic truth.

### 5. Noisy Circle With Nonuniform Sampling

Generate a circle with dense angular sectors and sparse angular sectors. Keep
radial noise moderate.

Purpose: evaluate whether the global MST-edge quantile is too permissive in
dense sectors or too conservative across sparse sectors.

Expected failure mode: low `q.thld` may leave the graph close to an MST in
sparse regions, while high `q.thld` may add radial or chord shortcuts in dense
regions.

### 6. Thick Noisy Circle or Annulus

Generate samples from a tubular neighborhood of a circle or from a thin
annulus. Optionally add radial structure so points have both angular and radial
variation.

Purpose: test whether CMST reconstructs the intended circular skeleton or
instead treats radial variation as equally important. This is a useful bridge
between pure manifolds and real data with biological dispersion around a
trajectory.

Expected failure mode: cross-thickness shortcuts can reduce angular geodesic
distances and blur circular topology.

### 7. Noisy Circle With a Gap

Generate a circle with an intentionally missing angular sector. Run two
interpretations:

- **Sampling-gap interpretation**: the missing sector is unobserved but the
  latent process is still continuous.
- **True-gap interpretation**: the missing sector represents a real forbidden
  or absent state.

Purpose: make the scientific prior explicit. CMST always returns a connected
graph, so the question becomes whether the forced MST bridge should be treated
as a legitimate geodesic connection or flagged as a long uncertain bridge.

Key metrics: bridge length relative to MST-edge distribution, number of
completion edges spanning the gap, and downstream sensitivity to removing or
downweighting the bridge.

### 8. High-Dimensional Embedded Noisy Circle

Embed a noisy circle in many dimensions with nuisance coordinates, anisotropic
noise, and optional low-rank batch effects. Evaluate with and without PCA.

Purpose: test the default PCA path in `create.cmst.graph()` and determine when
ambient nuisance variation destroys local Euclidean geometry.

Expected comparisons: CMST on raw features, CMST after PCA, CMST after known
latent two-dimensional projection, and oracle weighted circle graph.

### 9. Spiral and S-Curve

Generate open nonlinear manifolds where nearby ambient points can be far apart
along latent geodesics.

Purpose: stress-test chord avoidance and quantify shortest-path distortion on
classic manifold-learning examples.

Expected failure mode: overly high completion thresholds create edges between
adjacent turns of the spiral or folds of the S-curve.

### 10. Branching Y and Asymmetric Tree

Generate a Y-shaped graph with unequal branch lengths, unequal branch
densities, and off-branch noise. Add variants with one short dense branch and
one long sparse branch.

Purpose: test branch-point recovery, endpoint distances, and whether completion
creates shortcuts between branches.

Key metrics: branch assignment of edges, endpoint-to-endpoint distance
distortion, and false cross-branch edge rate.

### 11. Loops With Self-Approach: Figure-Eight

Generate a figure-eight or lemniscate with noise. The crossing can be either a
true shared vertex or two nearby but separate strands.

Purpose: distinguish topological crossing from ambient proximity.

Expected failure mode: CMST may connect nearby strands at the crossing even
when the latent structure treats them as separate.

### 12. Two Basins Connected by a Narrow Transition

Generate two dense clusters connected by a narrow curved bridge. Vary bridge
density and noise.

Purpose: model transitional biological states. CMST should preserve the bridge
without flooding the graph with direct cluster-to-cluster shortcuts.

Key metrics: bridge edge retention, cluster shortcut rate, and geodesic
distance from each cluster center through the transition.

### 13. Outliers and Contamination

Add isolated points, small off-manifold clusters, or contamination samples to
otherwise clean manifolds.

Purpose: determine how forced MST connectivity affects geodesic distances in
the presence of outliers.

Expected output: identify long MST edges and quantify whether completion
creates additional support for outlier connections.

### 14. Compositional Cyclic Gradient

Simulate microbiome-like count data from a latent circular ecological gradient.
Let latent angle drive taxa through periodic abundance waves. For taxon `j`,
define an angular optimum, amplitude, width, and baseline; transform these
periodic abundance functions into sample-specific compositions. Generate
observed sample counts from either a Dirichlet-multinomial model or a
logistic-normal multinomial model, with library-size variation and zero
inflation.

Evaluate CMST after compositional preprocessing:

- raw relative abundance,
- pseudocount plus CLR transform,
- Aitchison distance or Euclidean distance after CLR,
- optional filtering of rare taxa,
- optional PCA on transformed features.

Purpose: connect the toy noisy-circle examples to 16S rRNA and metagenomic
settings where cyclic or recurrent community-state transitions may be hidden in
sparse compositional count data.

Key metrics: recovery of latent circular geodesic distance, robustness to
library size and zeros, sensitivity to pseudocount choice, false shortcuts
between distinct ecological states, and stability of graph regression on
periodic response fields.

### 15. Compositional Branching Gradient

Simulate a latent Y-shaped ecological gradient where taxa respond to trunk and
branch coordinates. Include branch-specific taxa, shared core taxa, rare taxa,
and overdispersed counts.

Purpose: approximate developmental, disease-progression, or environmental
transition data where community composition follows branching trajectories.

Key metrics: branch geodesic recovery, cross-branch shortcut rate, endpoint
recovery, and basin stability.

### 16. Realistic 16S or Metagenomic Community-State Mixtures

Simulate enterotype-like basins or community state types connected by
transitional samples. Include batch effects, subject effects, sequencing depth
variation, taxon block swaps, and rare-feature sparsity.

Purpose: evaluate whether CMST can serve as a robust geometric backbone for
real microbiome data analysis, not just clean manifold examples.

Key metrics: basin separation, transitional path recovery, robustness to batch
perturbation, and agreement between graph geodesics and known latent ecological
gradients.

## Parameter Sweeps

For each case, sweep at least:

- `q.thld`: e.g. `0.5`, `0.7`, `0.9`, `0.95`, `0.99`.
- preprocessing: raw, scaled, PCA, and domain-specific transforms.
- noise level: low, medium, high.
- density imbalance: uniform, moderate imbalance, strong imbalance.
- sample size: small, medium, large.

The output should report not just the best setting, but the stability region:
the range of parameters for which CMST gives similar geodesic recovery without
false shortcuts.

## Source and Output Policy

Source assets live in this directory and should be version controlled. Generated
outputs should not be committed:

- reports go to `dev/geodesic-distance-estimation/reports/`,
- figures go to `dev/geodesic-distance-estimation/figures/`,
- temporary simulation outputs go to `dev/geodesic-distance-estimation/cache/`.

Small invariant tests for the package implementation belong in
`tests/testthat/`; large visual reports and exploratory benchmark sweeps belong
in this development workspace.

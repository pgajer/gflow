# Research Brief: Density-Aware Geodesic Distance Estimation

## Purpose

This brief defines the density-aware geodesic distance research track inside the
`gflow` geodesic distance estimation project.

The purpose of this track is to pursue density-aware distances for their own
sake: as a
mathematical object, as a practical graph/geodesic construction, as a bridge to
diffusion and PHATE-like geometry, and as a possible component of the broader
gradient-flow analysis framework supported by `gflow`.

Any agent or developer working from this brief should proceed slowly and
critically. Do not assume the idea is novel. Do not reduce the task to
implementation before clarifying definitions, existing literature, failure
modes, and biological use cases.

## Immediate Context

The current development workspace is:

```text
dev/geodesic-distance-estimation/
```

This workspace began as a Minimal Spanning Tree Completion Graph project for
`create.cmst.graph()`, but it has grown into a broader study of geodesic
distance estimation from sampled data. CMST, MST, iKNN, mKNN, geometric pruning,
component repair, oracle graph construction, and edge/path attribution are now
candidate tools rather than the whole objective.

The current density-aware idea emerged from noisy-circle work. The original
truth metric was latent circular arc distance. That is appropriate when radial
noise is tiny and the intended geometry is the original one-dimensional circle.
For thicker noisy circles, however, the sample is better viewed as a
density-supported tube around a circle. This led to a latent-tube oracle with
edge weights of the form:

```text
w(a,b) = ||z_a - z_b|| / rho((u_a + u_b)/2)^alpha
```

where `u` is signed radial displacement, `rho(u)` is the radial density, and
`alpha` controls the strength of density weighting.

The first sanity report found that `alpha = 0.5` behaves like a plausible
moderate density-aware truth, while `alpha = 1` can become tail-dominated unless
density floors, sample attachment, or robust truncation are handled carefully.

## Existing Project Assets

Read these files first:

```text
dev/geodesic-distance-estimation/README.md
dev/geodesic-distance-estimation/benchmark_plan.md
dev/geodesic-distance-estimation/report_template.md
dev/geodesic-distance-estimation/R/graph_construction.R
```

Important case scripts:

```text
dev/geodesic-distance-estimation/cases/noisy_circle_periodic_signal.R
dev/geodesic-distance-estimation/cases/noisy_circle_strategy_comparison.R
dev/geodesic-distance-estimation/cases/noisy_circle_gda_ceep_iknn_mknn.R
dev/geodesic-distance-estimation/cases/noisy_circle_dense_size_sweep.R
dev/geodesic-distance-estimation/cases/noisy_circle_latent_tube_oracle_sanity.R
```

Important notes:

```text
dev/geodesic-distance-estimation/notes/iknn_mknn_geodesic_geometry_reflections.md
dev/geodesic-distance-estimation/notes/noisy_circle_gda_ceep_testing_strategy.md
```

Generated reports and figures are ignored by git, but they are locally useful:

```text
dev/geodesic-distance-estimation/reports/noisy-circle-latent-tube-oracle-sanity/noisy_circle_latent_tube_oracle_sanity.html
dev/geodesic-distance-estimation/reports/noisy-circle-dense-size-sweep/noisy_circle_dense_size_sweep.html
dev/geodesic-distance-estimation/figures/noisy-circle-latent-tube-oracle-sanity/
dev/geodesic-distance-estimation/figures/noisy-circle-dense-size-sweep/
```

Core helper functions already implemented in `R/graph_construction.R`:

```r
tube.coords()
tube.density()
latent.tube.oracle()
latent.tube.oracle.distances()
geodesic.metrics()
graph.summary.metrics()
score.graph.metrics()
```

## Scientific Motivation

The central question is:

> Can we define and estimate geodesic distances from data in a way that rewards
> paths through high-density, well-supported regions and penalizes paths through
> low-density, weakly supported, or biologically implausible regions?

This matters because real data, especially microbiome and metagenomic data, are
not clean samples from noiseless manifolds. They are finite, uneven, noisy,
compositional, sparse, batch-affected, and shaped by biological population
structure. A purely Euclidean graph geodesic can add shortcuts across ambient
space. A purely latent manifold metric may ignore meaningful biological
dispersion. A density-aware distance may provide a controlled interpolation
between geometry and support.

## Publication Position To Evaluate

The idea should not be presented as wholly new. Density-driven and
density-weighted distances already exist under names such as Fermat distance,
power-weighted shortest path distance, density-based metrics, density-weighted
geodesics, and diffusion-geometry normalizations.

Current working view:

- A standalone publication claiming invention of density-aware distance is
  probably not appropriate.
- A standalone methods/application paper may become appropriate if `gflow`
  contributes a practical, validated framework for real biological data:
  density estimation, graph construction, parameter calibration, diagnostics,
  edge/path attribution, and downstream gradient-flow analysis.
- Otherwise the concept should be promoted as a principled module inside the
  larger gradient-flow analysis framework.

Your task is to make this judgment evidence-based.

## Known External Anchors

Investigate at least the following literature and connect it explicitly to the
`gflow` direction.

### Fermat Distances

Fermat distances are density-driven metrics on sampled manifolds. A sample
version often has the form:

```text
d_p(x,y) = inf over paths sum ||x_{i+1} - x_i||^p
```

for `p > 1`. Larger `p` favors paths made of many short edges, which tends to
keep geodesics in high-density regions. Continuum limits can be interpreted as
density-distorted Riemannian metrics.

Primary starting points:

- `Fermat Distances: Metric Approximation, Spectral Convergence, and Clustering Algorithms`
  JMLR 2024, https://jmlr.org/papers/v25/23-0939.html
- `Intrinsic Persistent Homology via Density-based Metric Learning`
  JMLR 2023, https://www.jmlr.org/papers/volume24/21-1044/21-1044.pdf

Questions:

- How does the Fermat exponent `p` translate to an explicit density exponent
  like `rho^-alpha`?
- What dimensionality assumptions enter this relationship?
- Can a `gflow` density-aware distance be parameterized so that Fermat distance
  is a special case?
- Does the noisy-circle latent-tube oracle reproduce known Fermat behavior?

### Power-Weighted Shortest Path Distances

Power-weighted shortest path distances explicitly study the tradeoff between
geometry and density in high-dimensional data.

Primary starting point:

- `Balancing Geometry and Density: Path Distances on High-Dimensional Data`,
  SIAM Journal on Mathematics of Data Science 2022,
  https://epubs.siam.org/doi/10.1137/20M1386657

Questions:

- What parameter regimes are recommended?
- What finite-sample guarantees exist for kNN graph approximations?
- Are there computational shortcuts we can reuse?
- How do PWSPDs differ from explicit density-estimated edge weights?

### Diffusion Maps, Diffusion Distance, and Density Normalization

Diffusion maps do not define a shortest-path metric, but their transition
operator is density-sensitive unless normalized. The diffusion maps
`alpha`-normalization controls the influence of sampling density. Diffusion
distance is small when many high-probability paths connect two points.

Starting points:

- Coifman and Lafon-style diffusion maps and anisotropic normalization.
- pydiffmap documentation for `alpha` and density functions:
  https://pydiffmap.readthedocs.io/en/master/reference/diffusion_map.html
- PHATE documentation:
  https://phate.readthedocs.io/en/stable/

Questions:

- In what precise sense does diffusion distance prefer high-density paths?
- Is this preference shortest-path-like, random-walk-like, or potential-based?
- How does PHATE potential distance transform diffusion probabilities?
- Can PHATE potential distance be treated as an implicit density-aware distance?
- When would a density-aware shortest-path distance and a diffusion/potential
  distance agree or disagree?

### PHATE and Potential Distance

PHATE builds an affinity graph, diffuses it for time `t`, transforms diffusion
probabilities into potentials, and embeds using distances between potentials.
The PHATE documentation describes `t` as the diffusion level and `gamma` as the
informational distance constant, with `gamma = 1` corresponding to log
potential and `gamma = 0` to square-root potential.

Questions:

- Does PHATE potential distance effectively measure support through
  high-density paths?
- Is the effect mediated mostly by the kernel graph, the diffusion operator, the
  diffusion time, or the potential transform?
- Can we define a PHATE-derived geodesic distance or a density-aware graph
  distance that aligns with PHATE potential distances?
- How should `gflow` expose PHATE-like density-aware options without confusing
  shortest-path geometry with diffusion geometry?

### Current Internal References

The user specifically noted:

- arXiv:2511.03817, `Adaptive Geometric Regression for High-Dimensional
  Structured Data`, where data density is estimated.
- arXiv:2508.02080, `The Geometry of Machine Learning Models`, by Pawel Gajer
  and Jacques Ravel, https://arxiv.org/abs/2508.02080

The quoted internal text states:

```text
Drawing inspiration from optimal transport theory and information geometry,
we define a density-weighted metric:
d_rho(e) = ell(e) / rho(e)^alpha
where ell(e) is the original Riemannian length of edge e and alpha in (0, 1]
controls the strength of density weighting. This modification makes paths
through high-density regions effectively shorter than those through sparse
regions.
```

Treat this as an internal conceptual anchor, not as a novelty claim. Determine
where it overlaps with known Fermat/PWSPD/diffusion geometry and where the
`gflow` implementation can still be distinctive.

## Formalization Tasks

Develop a taxonomy of density-aware distances. At minimum include these
families.

### 1. Continuous Density-Weighted Path Metrics

General form:

```text
L_{rho,alpha}(gamma) = integral_gamma rho(x)^(-alpha) ds
d_{rho,alpha}(x,y) = inf_gamma L_{rho,alpha}(gamma)
```

Extensions to consider:

- density floor: `rho_floor`
- robust density transform: `phi(rho)` instead of raw `rho`
- capped penalty: `min(rho^-alpha, penalty_max)`
- local uncertainty penalty
- anisotropic local metric tensor
- density on a manifold, tube, graph, or ambient feature space

Questions:

- When is this a true metric?
- What happens when density is zero or badly estimated?
- How should `alpha` scale with intrinsic dimension?
- Should the path length be symmetric, or can directed/asymmetric variants be
  useful for trajectory or gradient-flow analysis?

### 2. Graph Edge Reweighting

Given a graph with base edge length `ell_ij`, define:

```text
w_ij = ell_ij / phi(rho_i, rho_j, rho_ij)^alpha
```

Candidate edge density summaries:

- arithmetic mean of endpoint densities
- geometric mean of endpoint densities
- harmonic mean
- minimum endpoint density
- midpoint density, if ambient coordinates permit
- edge-neighborhood density estimated from samples near the edge
- random-walk stationary mass or personalized mass

Questions:

- Which summary is most robust?
- Which one best penalizes edges crossing low-density gaps?
- Should long edges and short edges use different density summaries?
- How should floors/caps be chosen?

### 3. Density-Aware Neighbor Selection

Modify graph construction before shortest paths:

- density-aware kNN
- radius graphs with density-dependent radius
- mutual kNN with density correction
- local scaling using kNN radii or density estimates
- bridge admission rules that require both geometric and density support
- density-aware geometric pruning thresholds

Questions:

- Is it better to change the graph topology first or only reweight edges?
- Can density-aware kNN reduce shortcut edges before pruning?
- How does this interact with mKNN fragmentation?
- Can MST/CMST bridges be flagged or weighted by density support?

### 4. Fermat / Power-Weighted Sample Distances

Study sample path metrics:

```text
d_p(x,y) = inf_path sum ||x_{i+1} - x_i||^p
```

Questions:

- When should `gflow` implement this directly?
- Is it computationally feasible on kNN/iKNN/mKNN graphs?
- How should distances be rescaled across `n`, `p`, and intrinsic dimension?
- Can this serve as a density-aware truth for benchmarks?

### 5. Diffusion and PHATE-Derived Distances

Study distances derived from transition probabilities:

- diffusion distance
- PHATE potential distance
- commute time / resistance distance
- heat kernel distance
- personalized PageRank distance
- random-walk hitting or commute distances

Questions:

- Which are genuinely density-aware?
- Which are robust to sampling density versus intentionally density-sensitive?
- Which are suitable as geodesic distances for downstream `gflow` regression?
- Which are only embedding distances and should not be called geodesic?

### 6. Density-Aware Oracle Truth

For synthetic benchmarks with known generative density, define oracle truth via
structured grids, dense samples, or analytic approximations.

Existing implementation:

```r
latent.tube.oracle()
latent.tube.oracle.distances()
```

Questions:

- How should `alpha`, density floor, and support truncation be chosen?
- How should observed samples attach to the oracle graph?
- How sensitive are oracle distances to grid resolution?
- Should oracle distance be compared to latent arc, Euclidean, Fermat, and
  diffusion distances?

## Practical `gflow` Design Questions

Investigate possible API directions, but do not implement public APIs before
the definitions stabilize.

Possible future options:

```r
estimate.density(...)
create.density.aware.graph(...)
density.weight.graph(...)
geodesic.distance(..., metric = "density_weighted")
create.iknn.graphs(..., density.weighting = ...)
create.cmst.graph(..., density.weighting = ...)
phate.core(..., density.aware = ...)
```

Questions:

- Should density-aware distance be a graph-construction option, a distance
  post-processing option, or both?
- How should density estimates be represented and reused?
- Should density be estimated in raw feature space, PCA space, CLR/Aitchison
  space, graph space, or latent/oracle space?
- What diagnostics must be reported whenever density-aware distance is used?
- How do we protect users from extreme tail behavior?

## Density Estimation Ideas To Explore

Pre-graph density estimates:

- kNN radius density
- kernel density in reduced space
- local PCA dimension-corrected density
- adaptive bandwidth density
- compositional density after CLR/Aitchison preprocessing
- model-based density, including normalizing-flow or mixture models where
  appropriate

Graph-based density estimates:

- degree or weighted degree
- stationary distribution of a random walk
- local volume growth
- heat kernel mass at scale `t`
- personalized random-walk mass from a center vertex
- return probability or diffusion entropy
- density estimated from graph balls at multiple geodesic radii

Random-walk centered idea from the user:

For a random walk initiated at vertex `i`, the distribution after time `t`
defines a local cloud with center `i`, scale, entropy, and effective support.
This can be interpreted as a local density/volume probe around `i`.

Questions:

- Can the scale of this random-walk distribution define a vertex density?
- Can it define density along an edge?
- Does it connect to diffusion distance, PHATE potential distance, or heat
  kernel signatures?
- Does it reveal local bottlenecks or sparse bridges better than kNN density?

## Benchmark Program

Start from existing noisy-circle assets, then expand.

### Phase A: Noisy Circle / Tube

Use:

```text
cases/noisy_circle_latent_tube_oracle_sanity.R
cases/noisy_circle_dense_size_sweep.R
```

Tasks:

- Compare latent arc, Euclidean, tube oracle, Fermat, and diffusion/PHATE
  distances.
- Stress-test `alpha = 0, 0.25, 0.5, 0.75, 1`.
- Vary density floors and robust caps.
- Diagnose whether extreme tail behavior is an oracle artifact or desired
  behavior.
- Add edge/path attribution:
  - high oracle/Euclidean edge ratio,
  - high true-distance endpoint edges,
  - missing local connectors,
  - graph geodesics that underestimate truth,
  - graph geodesics that overestimate truth.

### Phase B: Nonuniform Curves and Gaps

Use synthetic cases where density variation is intentional:

- uneven noisy line
- crescent/open arc
- two nearby non-touching arcs
- circle with gap
- circle with dense and sparse angular sectors

Tasks:

- Ask whether density-aware distances distinguish sampling gaps from true
  forbidden gaps.
- Compare density-aware MST/CMST/iKNN/mKNN strategies.
- Determine when density-aware weighting helps and when it falsely separates
  valid sparse regions.

### Phase C: Compositional and Microbiome-Like Data

Use:

- compositional cyclic gradient
- Dirichlet-multinomial or logistic-normal multinomial counts
- CLR/Aitchison-like preprocessing
- 16S rRNA and metagenomic-inspired synthetic examples

Tasks:

- Decide where density should be estimated: counts, relative abundance,
  transformed space, graph space, or model latent space.
- Test sensitivity to library size, zero inflation, and batch effects.
- Ask whether density-aware distances improve conditional expectation
  estimation or trajectory/regression stability.

### Phase D: Real Data

Only after the synthetic program is stable:

- identify real datasets with plausible gradients, trajectories, or community
  states;
- compute density-aware and density-unaware distances;
- compare biological interpretability, stability, and downstream model fit.

## Diagnostics To Require

Every density-aware distance report should include:

- density estimate distribution;
- density floor/cap diagnostics;
- edge length versus density-weighted edge length;
- edge density summaries;
- path examples through high-density and low-density regions;
- sensitivity to `alpha`;
- sensitivity to density estimator;
- correlation with Euclidean, arc/oracle, Fermat, diffusion, and PHATE
  distances when available;
- underestimation and overestimation path attribution;
- tail/outlier influence diagnostics;
- stability across seeds and sample sizes.

## Deliverables

Produce these deliverables in order:

1. **Literature map**
   - Fermat distances, PWSPD, diffusion maps, PHATE, density-based persistent
     homology, conformal metrics, weighted geodesics, random-walk distances.
   - Include precise equations and parameter correspondences.

2. **Definition memo**
   - Formal family of density-aware distances.
   - Continuous, graph, sample-Fermat, diffusion, and PHATE variants.
   - Conditions under which each is a metric or pseudometric.

3. **Implementation design memo**
   - Proposed internal R helpers.
   - Proposed future public API.
   - Required diagnostics and defaults.

4. **Noisy-circle/tube benchmark report**
   - Refine the existing latent-tube oracle.
   - Compare `alpha`, floors, caps, and attachment rules.
   - Include edge/path attribution.

5. **Publication recommendation**
   - Decide whether this is:
     - a section inside the larger gradient-flow analysis framework,
     - a standalone methods note,
     - an application paper focused on microbiome/metagenomic data,
     - or not yet publication-ready.

6. **Risk register**
   - Novelty risk.
   - Overfitting/parameter-tuning risk.
   - Density estimation instability.
   - Tail domination.
   - Mislabeling diffusion distances as geodesics.
   - Biological interpretability risk.

## Working Hypotheses

Treat these as hypotheses, not conclusions:

1. Density-aware distance is not new, but it may be underused in practical
   microbiome/metagenomic graph construction.
2. `rho^-alpha` edge weighting is a useful unifying language, but Fermat and
   diffusion distances may provide better-studied alternatives.
3. Moderate density weighting, such as `alpha` near `0.5`, may be more useful
   than aggressive `alpha = 1` in noisy finite samples.
4. Density-aware weighting should probably be paired with explicit tail
   handling: floors, caps, robust density transforms, or outlier filtering.
5. Diffusion/PHATE distances are density-aware in a broad support/connectivity
   sense, but they are not the same object as shortest-path geodesic distances.
6. The practical contribution may be not the distance formula alone, but the
   full workflow: density estimation, graph construction, diagnostics,
   benchmarks, and downstream `gflow` integration.

## Style Of Work

- Be explicit about what is known, what is borrowed, and what is new.
- Prefer primary sources.
- Keep mathematical definitions separate from implementation choices.
- Run small sanity examples before dense sweeps.
- Generate HTML reports with figures and interpretation.
- Do not promote a public API until diagnostics are convincing.
- Always ask whether the distance improves the scientific task, not only an
  abstract benchmark score.

## First Task

Begin by producing Deliverable 1: a literature map of density-aware and
density-adjacent distances.

This literature map should be mathematically precise and historically grounded.
For each concept or construct, include:

- the original motivation and problem setting;
- a short historical note, including key papers and terminology;
- the formal definition, with equations and clearly stated objects;
- whether the construct is a metric, pseudometric, graph distance, embedding
  distance, random-walk distance, or sampling-limit object;
- the role played by data density, sampling density, graph degree, transition
  probability, or kernel normalization;
- parameter meanings and correspondences across frameworks, where known;
- computational construction from finite data;
- known asymptotic or finite-sample behavior;
- practical strengths, failure modes, and diagnostic requirements;
- relevance to `gflow`, geodesic distance estimation, noisy-tube benchmarks,
  and microbiome/metagenomic data.

In addition to mapping the theory, perform an explicit practical-adoption
review. Search for evidence that each concept has been used in real data
analysis, not only proposed mathematically. If the agent environment supports
delegation, consider spawning focused sub-agents for separate application
areas. At minimum, investigate:

- bioinformatics broadly defined;
- microbiome, 16S rRNA, metagenomics, metatranscriptomics, and compositional
  count data;
- single-cell genomics, spatial transcriptomics, developmental trajectories,
  and cell-state manifolds;
- ecological gradients and community-composition data;
- medical imaging, morphology, or other biological shape/phenotype data;
- non-biological applied domains where practical lessons may transfer.

For each application area, report:

- whether density-aware, Fermat-like, power-weighted, diffusion/potential, or
  density-normalized graph distances are actually used;
- the exact name under which the method appears in that literature;
- what data type and scientific question motivated its use;
- whether the distance was used for visualization, clustering, trajectory
  inference, regression, classification, topology, denoising, or downstream
  modeling;
- whether the paper provides software, defaults, diagnostics, or sensitivity
  analyses;
- whether density weighting was central to the method or only an implicit
  property of a graph/diffusion construction;
- gaps where the concept appears theoretically mature but practically
  underused.

At minimum, cover:

- density-weighted Riemannian or conformal path metrics;
- Fermat distances;
- power-weighted shortest path distances;
- density-based metric learning for topology and persistent homology;
- kNN, mKNN, iKNN, and locally scaled graph geodesics;
- diffusion distance and diffusion-map normalization;
- PHATE potential distance;
- random-walk, commute-time, resistance, heat-kernel, and personalized
  PageRank-style distances;
- density-aware neighbor selection and edge reweighting;
- density-aware oracle constructions for synthetic benchmarks.

End the literature map with a synthesis table that separates:

- constructs that are explicitly density-weighted geodesic distances;
- constructs that implicitly favor high-density paths;
- constructs that correct for sampling density rather than exploit it;
- constructs that are useful analogies but should not be called geodesic
  distances.

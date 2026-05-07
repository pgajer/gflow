# Data Graph Constructions Report Brief

This brief covers background information, assets, and tasks related to the
LaTeX report

- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/phate_knn_graph_constructions.tex`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/build/phate_knn_graph_constructions.pdf`

The existing report began as an explanation of the relationship between PHATE's
adaptive affinity support and k-nearest-neighbor graph variants. The report now
needs to become a broader, self-contained document about data-to-graph
constructions used for geodesic distance approximation.

## Recommended Report Rename

The current filename and title are too narrow:

```tex
\title{PHATE, Symmetric kNN, Mutual kNN, and Intersection kNN Graph Constructions}
```

Recommended title:

```tex
\title{Data Graph Constructions for Geodesic Distance Approximation}
```

Recommended output name:

- `data_graph_constructions.pdf`, or
- `data_to_graph_geodesic_constructions.pdf`

The shorter `data_graph_constructions.pdf` is probably best. The report should
still contain PHATE, but PHATE should be positioned as an adaptive affinity and
diffusion construction whose support is closely related to symmetric kNN and
adaptive-radius graphs. PHATE is not primarily a graph-geodesic distance
construction in the same sense as sKNN, mKNN, iKNN, fixed-radius, and
adaptive-radius graphs with edge lengths.

Also replace the current date command

```tex
\date{\today}
```

with a build timestamp line, matching other Codex project reports. For example:

```tex
\date{Build timestamp: \today}
```

If the build script can inject a full timestamp, prefer that over `\today`.

## Scientific Context

The motivating problem is the data geodesic geometric reconstruction problem.
Given a finite sample

$$
X = \{x_1,\ldots,x_N\}\subset \mathbb{R}^p
$$

from a structured geometric object, we want to construct a weighted graph

$$
G(X) = (V,E,\ell), \qquad V = X,
$$

whose shortest-path metric approximates the intrinsic geodesic geometry of the
sampled object. If \(d_M\) denotes the geodesic distance on the underlying
space \(M\), and \(d_{G(X)}\) denotes graph shortest-path distance, the guiding
goal is

$$
d_{G(X)}(x_i,x_j) \approx d_M(x_i,x_j)
$$

for all sampled pairs \(x_i,x_j\). In the finite-sample setting, we also compare
against a sample-oracle geodesic distance \(d_X^{\mathrm{oracle}}\), which asks
how well distances can be reconstructed using paths through the observed sample
points near a reference geodesic.

This problem matters because Euclidean distances in the ambient space can be
misleading on nonlinear or branched data. Two points can be close in
\(\mathbb{R}^p\) while far along the intrinsic geometry, for example across a
fold, near a self-approach, or across two crossing branches. Conversely, graph
shortest paths can recover manifold-like or branch-like geometry if the graph
support and edge lengths are chosen well.

The current benchmark focus is deliberately narrow:

$$
y = x_0^2 + x_1^2
$$

and

$$
y = x_0^2 - x_1^2,
$$

sampled over a parameter disk. These two quadratic graph surfaces are useful
because they are simple enough to visualize and to compute reference geodesics,
but still curved enough to expose distortion in graph construction.

## How PHATE Led To This Report

The investigation started with an audit of `gflow`'s PHATE implementation. A key
observation was that PHATE graph construction is not an ordinary kNN graph. It
first builds a directed adaptive affinity

$$
A_{ij}
= \exp\left[-\left(\frac{d_{ij}}{\sigma_i}\right)^a\right],
$$

where \(\sigma_i\) is typically the distance from \(x_i\) to its \(k\)-th
nearest neighbor. After thresholding small affinities at \(\tau\), a directed
edge survives when

$$
d_{ij} \le r_\tau \sigma_i,
\qquad
r_\tau = (-\log\tau)^{1/a}.
$$

With additive symmetrization, PHATE has undirected support

$$
d_{ij} \le r_\tau \max(\sigma_i,\sigma_j).
$$

For PHATE's common values \(a=40\) and \(\tau=10^{-4}\),

$$
r_\tau \approx 1.057.
$$

Thus PHATE's support is very close to a slightly inflated symmetric kNN support.
This observation led to the broader question: if the PHATE support geometry is
often useful, should we study sKNN, mKNN, iKNN, fixed-radius, and adaptive-radius
graphs directly as data-to-graph geodesic approximation methods?

The answer is yes. The new report should present these graph constructions as
competing ways to approximate intrinsic geodesic geometry.

## Current Report Assets

Primary report files:

- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/phate_knn_graph_constructions.tex`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/build/phate_knn_graph_constructions.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figure_generation.log`

Existing figures:

- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/roots_n3_k2.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/roots_n4_k2.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/roots_n5_k2.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/roots_n100_k10.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n50_k8.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n50_k8_diagnostics.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n100_k8.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n100_k8_diagnostics.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n150_k8.pdf`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/figures/noisy_circle_fermat_n150_k8_diagnostics.pdf`

Related PHATE graph comparison report:

- `/Users/pgajer/current_projects/gflow/dev/phate-graph-comparison/report/phate_iknn_graph_comparison_report.html`
- `/Users/pgajer/current_projects/gflow/dev/phate-graph-comparison/phate_iknn_graph_comparison_action_plan.md`

Current geodesic benchmark design notes:

- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/geodesic_distance_literature_map.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/density_aware_geodesic_distance_literature_map.md`

Current benchmark runner area:

- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/`
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R`

Earlier noisy-circle cases that may contain useful examples and diagnostics:

- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_gda_ceep_iknn_mknn.R`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_tube_benchmark_report.R`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_latent_tube_oracle_sanity.R`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_periodic_signal.R`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_strategy_comparison.R`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/cases/noisy_circle_dense_size_sweep.R`

PHATE clean-room implementation planning assets:

- `/Users/pgajer/current_projects/gflow/dev/notes/phate_clean_room_cpp_plan_2026-05-03.md`
- `/Users/pgajer/current_projects/gflow/dev/phate-phase6-validation/phase6_action_plan.md`

## Current gflow Implementation Assets

The report should point readers to the package functions used in current
geodesic distance approximation explorations.

Graph constructors:

- `create.sknn.graph()` in `/Users/pgajer/current_projects/gflow/R/sknn_graphs.R`
- `create.mknn.graph()` in `/Users/pgajer/current_projects/gflow/R/mknn_graphs.R`
- `create.single.iknn.graph()` in `/Users/pgajer/current_projects/gflow/R/iknn_graph.R`
- `create.radius.graph()` in `/Users/pgajer/current_projects/gflow/R/radius_graphs.R`
- `create.adaptive.radius.graph()` in `/Users/pgajer/current_projects/gflow/R/radius_graphs.R`

Core C++ graph implementations:

- `/Users/pgajer/current_projects/gflow/src/sknn_graphs.cpp`
- `/Users/pgajer/current_projects/gflow/src/mknn_graphs.cpp`
- `/Users/pgajer/current_projects/gflow/src/iknn_graphs.cpp`

Connectivity repair:

- `/Users/pgajer/current_projects/gflow/R/component_mst_connectivity.R`

Geometric pruning:

- `/Users/pgajer/current_projects/gflow/R/local_geodesic_pruning.R`
- C++ local pruning support in `/Users/pgajer/current_projects/gflow/src/sknn_graphs.cpp`
- iKNN-specific pruning support in `/Users/pgajer/current_projects/gflow/src/iknn_graphs.cpp`

Graph shortest-path wrapper:

- `graph.geodesic.distances()` in `/Users/pgajer/current_projects/gflow/R/graph_geodesic_distances.R`

Quadratic-surface reference geodesic utilities:

- `quadform.embed()`
- `quadform.edge.length()`
- `quadform.reference.geodesics()`
- `quadform.grid.geodesic.distances()`
- `quadform.grid.geodesic.calibration()`
- `quadform.sample.dataset()`

These are in `/Users/pgajer/current_projects/gflow/R/quadform_geodesics.R`, with
C++ support in `/Users/pgajer/current_projects/gflow/src/quadform_grid_geodesics.cpp`.

Deviation-from-isometry helpers:

- `isometry.scale()`
- `isometry.rel.rms.error()`
- `isometry.rel.abs.error()`
- `isometry.distortion.quantiles()`
- `isometry.distance.correlations()`
- `summarize.isometry.deviation()`

These are in `/Users/pgajer/current_projects/gflow/R/isometry_deviation.R`.

## Graph Constructors To Describe

The report should contain a precise algorithmic description of every graph
constructor currently used in the geodesic distance approximation benchmark.

### Directed kNN

Given distances \(d_{ij}\), define \(N_k(i)\) as the \(k\) nearest non-self
neighbors of \(i\). The directed kNN graph has edge

$$
i\to j
\quad\Longleftrightarrow\quad
j\in N_k(i).
$$

This is a search primitive rather than the final graph used for graph geodesics.

### Symmetric kNN / union kNN / sKNN

The symmetric kNN graph has edge

$$
\{i,j\}\in E_{\mathrm{sKNN}}
\quad\Longleftrightarrow\quad
j\in N_k(i)\ \text{or}\ i\in N_k(j).
$$

When there are no ties, if \(\sigma_i=d_{i,(k)}\), this is equivalent to

$$
d_{ij}\le \max(\sigma_i,\sigma_j).
$$

Use `sKNN` as the main abbreviation. Mention that `uKNN` can be useful when
emphasizing the union/OR symmetrization, but `sKNN` is more standard.

### Mutual kNN / mKNN

The mutual kNN graph has edge

$$
\{i,j\}\in E_{\mathrm{mKNN}}
\quad\Longleftrightarrow\quad
j\in N_k(i)\ \text{and}\ i\in N_k(j).
$$

Without ties this becomes

$$
d_{ij}\le \min(\sigma_i,\sigma_j).
$$

This graph is more conservative than sKNN and is more prone to disconnected
components at small \(k\), but it can suppress cross-density and false-shortcut
edges.

### Intersection kNN / iKNN

The iKNN graph uses closed neighborhoods

$$
\widehat{N}_k(i)=\{i\}\cup N_k(i)
$$

and has edge

$$
\{i,j\}\in E_{\mathrm{iKNN}}
\quad\Longleftrightarrow\quad
\widehat{N}_k(i)\cap \widehat{N}_k(j)\ne\varnothing.
$$

Unlike sKNN and mKNN, iKNN is not equivalent to a simple radius predicate in
\(d_{ij}\). In `gflow`, the native iKNN edge length is an intersection-mediated
detour length

$$
\ell_{ij}^{\mathrm{iKNN}}
=
\min_{c\in \widehat{N}_k(i)\cap\widehat{N}_k(j)}
\left(d_{ic}+d_{jc}\right).
$$

For the current geodesic reconstruction benchmark, be explicit about whether
ambient Euclidean edge lengths or native iKNN edge lengths are being used.

### Fixed-Radius Graph

The fixed-radius graph has edge

$$
\{i,j\}\in E_\varepsilon
\quad\Longleftrightarrow\quad
d_{ij}\le \varepsilon.
$$

In the benchmark, \(\varepsilon\) should be chosen from interpretable quantiles
or neighbor-rank distances so that fixed-radius graphs can be compared fairly
against kNN-style graphs. Edge lengths are normally

$$
\ell_{ij}=\|x_i-x_j\|_2.
$$

This section still needs to be added to the LaTeX report.

### Adaptive-Radius Graph

Define a local scale

$$
\sigma_i=d_{i,(k_{\mathrm{scale}})}.
$$

The adaptive-radius graph uses edge rules such as

$$
d_{ij}\le c\max(\sigma_i,\sigma_j)
$$

or

$$
d_{ij}\le c\min(\sigma_i,\sigma_j),
$$

where \(c\) is `radius.factor`. The `max` rule is sKNN-like; the `min` rule is
mKNN-like. With \(c=1\) and no ties, these match the sKNN and mKNN support
conditions above.

This section also needs to be added to the LaTeX report. It should explicitly
connect PHATE's thresholded support to the max-rule adaptive-radius graph:

$$
d_{ij}\le r_\tau \max(\sigma_i,\sigma_j).
$$

### PHATE Adaptive Affinity Support

PHATE should be included as a bridge between dimensionality reduction and
data-to-graph construction. Its edge support is adaptive-radius-like, but its
weights are affinities, not edge lengths. The row-normalized version is a Markov
diffusion operator, not a graph-geodesic metric by itself.

The report should warn against applying graph-geodesic methods directly to
PHATE affinities without converting similarities to lengths and justifying the
conversion.

## Pruning Methods To Add

The report should have a standalone pruning section with figures. This is
currently missing or underdeveloped.

### Local Geodesic Pruning

Local geometric pruning is an opt-in sparsification stage. For a candidate edge
\(\{u,v\}\) with length \(\ell_{uv}\), define a local vertex set
\(L(u,v)\), usually from the union of local neighborhoods around the endpoints.
Remove the edge temporarily and compute the shortest alternative path length

$$
d_{L(u,v)\setminus\{uv\}}(u,v).
$$

For tolerance \(\tau>1\), prune the edge if

$$
d_{L(u,v)\setminus\{uv\}}(u,v)
\le
\tau\,\ell_{uv}.
$$

Interpretation: if an edge can be replaced by a local path that is only a small
fraction longer, the direct edge is redundant for local geodesic geometry.

This is a local spanner-like pruning rule. It is related in spirit to geometric
spanners, but the current implementation is a pragmatic local data-graph
procedure rather than a formal spanner construction with proven global stretch.

Important implementation note: local pruning is currently much slower than
expected in benchmark runs, even on moderate examples. The HTML benchmark work
has exposed this as an implementation issue to investigate separately.

Suggested figures:

1. A small point cloud with a long chord and a short multi-edge path.
2. The same graph after the chord is removed.
3. A local neighborhood view showing that the alternative path is restricted to
   nearby vertices.

### iKNN Global Geometric Pruning

`create.single.iknn.graph()` also has older iKNN-specific pruning behavior. This
is faster in current tests. The report should distinguish it from the local
geodesic pruning now available across graph families.

The key point to communicate: global iKNN pruning and local geodesic pruning are
not the same algorithm. Do not describe them as interchangeable.

### Pruning Order

The package now keeps graph lifecycle fields so that benchmark code can compare
different stages:

- `raw`
- `raw.repaired`
- `pruned`
- `pruned.repaired`
- `repaired.pruned`
- `final`

For the first quadratic-surface benchmark, the preferred stage is:

```r
stage <- if (prune.method == "none") "raw.repaired" else "repaired.pruned"
```

This means we first repair connectivity, then prune locally. This avoids
comparing finite graph-geodesic matrices against disconnected graphs with
infinite distances.

## MST Repair To Add

The report should include a detailed MST repair section and an illustration.

If a graph has connected components

$$
C_1,\ldots,C_m,
$$

MST repair adds bridge edges between components so the final graph is connected.
The exact component-MST construction treats components as supernodes and defines
an intercomponent bridge length

$$
w(C_a,C_b)
=
\min_{i\in C_a,\ j\in C_b} \|x_i-x_j\|_2.
$$

Then compute an MST on the complete component graph and add the selected
shortest bridge edges back to the original graph.

This gives \(m-1\) bridge edges when there are \(m\) components. It does not
connect every pair of components; it connects the component graph minimally.

Current implementation supports:

- `connect.method = "component.mst"`
- `connect.method = "component.mst.ann"`
- `connect.method = "global.mst"`

The current benchmark uses `connect.components = TRUE` and
`connect.method = "component.mst"`.

Suggested figures:

1. Disconnected graph with components colored separately.
2. Complete component graph shown schematically with candidate bridge edges.
3. Repaired graph with only the MST bridge edges added.

## Extra 2D Examples To Add

The current root-of-unity and noisy-circle examples are useful, but the expanded
report should add a few more 2D examples that expose different failure modes.

Recommended additions:

1. Figure eight / lemniscate.
   - Purpose: show false shortcuts near self-approach or crossing-like geometry.
   - Useful for comparing sKNN, mKNN, iKNN, fixed radius, adaptive radius, and
     PHATE support.

2. Two line segments crossing or nearly crossing.
   - Purpose: show the ambiguity between ambient proximity and intrinsic branch
     membership.
   - Include both exact crossing and near crossing if possible.

3. Tree or branching structure.
   - Purpose: connects the report to PHATE-style biological trajectory examples.
   - Shows whether graph construction preserves branches or creates shortcuts
     between nearby arms.

4. Quadratic graph surfaces projected to 2D parameter and 3D embedded views.
   - Paraboloid: \(y=x_0^2+x_1^2\).
   - Saddle: \(y=x_0^2-x_1^2\).
   - These tie the graph-construction report to the current geodesic distance
     approximation benchmark.

5. Nonuniform noisy circle.
   - Purpose: show sensitivity to sampling density and why adaptive local scale
     constructions are important.

Optional additions:

- Swiss roll or S-curve for a familiar manifold-learning example.
- Two moons for clustering-oriented graph behavior.
- Spiral curve to show local-vs-global shortcut behavior.

## Historical Perspective And References

The expanded report should include a historical perspective section. The
following are good reference targets for the next agent to verify and add to the
LaTeX bibliography.

### Nearest-Neighbor Graphs And Manifold Learning

Isomap is a central motivation for graph-geodesic reconstruction. It estimates
manifold geodesic distances by shortest paths on a neighborhood graph and then
applies MDS. It explicitly uses kNN or fixed-radius neighborhood graphs.

Candidate citation:

- Tenenbaum, de Silva, and Langford, "A global geometric framework for nonlinear
  dimensionality reduction", Science, 2000.

### Fixed-Radius Graphs And Random Geometric Graphs

Fixed-radius graphs are the deterministic data-analysis analogue of random
geometric graphs. The classical random geometric graph origin is Gilbert's
random plane network model.

Candidate citations:

- Gilbert, "Random Plane Networks", SIAM Journal on Applied Mathematics, 1961.
- Penrose, "Random Geometric Graphs", 2003.

### Adaptive Local Scaling

Adaptive-radius and PHATE-like support rules are closely related to local-scale
kernel ideas in spectral clustering.

Candidate citation:

- Zelnik-Manor and Perona, "Self-Tuning Spectral Clustering", NeurIPS, 2004.

### Mutual kNN

Mutual kNN graphs are widely used to make kNN graph construction more
conservative and less sensitive to asymmetric neighbor relations. They are
common in clustering, spectral graph construction, and semi-supervised learning.

The next agent should verify the best historical citation. Useful search terms:

- "mutual k nearest neighbor graph clustering"
- "mutual kNN graph spectral clustering"
- "optimal construction k nearest neighbor graphs noisy clusters mutual symmetric"

### Shared Nearest Neighbor And Single-Cell Context

Shared nearest-neighbor graphs are not the same as iKNN. SNN graphs typically
weight edges by neighborhood overlap, often Jaccard similarity. They are
important in single-cell workflows, for example Seurat's `FindNeighbors()`.

The report should avoid confusion:

- Raw KNN in Seurat means each cell's nearest-neighbor list in a chosen
  embedding or feature space.
- SNN is a derived graph whose weights come from overlap of those neighbor
  lists.
- iKNN in `gflow` uses nonempty intersection of closed neighborhoods as an edge
  predicate and stores an intersection-mediated length.

Candidate references:

- Jarvis and Patrick, "Clustering Using a Similarity Measure Based on Shared
  Near Neighbors", IEEE Transactions on Computers, 1973.
- Seurat `FindNeighbors()` documentation for modern single-cell usage.

### Geometric Pruning And Spanners

Local geodesic pruning is related to spanner ideas: keep a sparse graph while
preserving shortest-path distances up to a multiplicative stretch.

Candidate references:

- Althoefer et al., "On sparse spanners of weighted graphs", 1993.
- Narasimhan and Smid, "Geometric Spanner Networks", 2007.
- Recent survey/tutorial references on graph spanners.

### MST Repair

The MST repair stage is a pragmatic graph-construction step rather than a
standard named manifold-learning graph. Cite classical MST algorithms only if
the report discusses algorithms:

- Kruskal, 1956.
- Prim, 1957.

## Other Additions I Recommend

### Terminology Table

Add a table separating:

- graph support,
- edge weights,
- whether the weights are lengths or affinities,
- whether row-normalization produces a Markov operator,
- whether the graph is intended for graph-geodesic distances.

This table will prevent the common confusion between PHATE affinities and graph
edge lengths.

### Edge Length Conventions

The report should state clearly that in the geodesic reconstruction benchmark,
most graphs use ambient Euclidean edge lengths:

$$
\ell_{ij}=\|x_i-x_j\|_2.
$$

iKNN can also use its native intersection-mediated length, and PHATE affinities
can be transformed into lengths only after choosing an explicit transformation,
such as

$$
\ell_{ij}=-\log(K_{ij}+\epsilon)
$$

or

$$
\ell_{ij}=K_{ij}^{-1/p}.
$$

Those transformations are research questions, not default graph-geodesic
conventions.

### Nesting And Non-Nesting Diagram

Add a small diagram or table showing no-tie support relations:

$$
E_{\mathrm{mKNN}}
\subseteq
E_{\mathrm{sKNN}}
$$

and, for PHATE with additive symmetrization,

$$
E_{\mathrm{sKNN}}
\subseteq
E_{\mathrm{PHATE}}
$$

when \(r_\tau\ge 1\). iKNN is not generally nested with PHATE because it uses a
shared-closed-neighborhood predicate, not a direct radius condition.

### Computational Complexity Section

Add a short practical section:

- exact pairwise distances are simple and deterministic but \(O(n^2)\);
- ANN neighborhood search is faster for kNN support;
- exact component MST bridge search can be expensive if many large components
  exist;
- local geodesic pruning is currently the main computational bottleneck in the
  benchmark.

### Evaluation Metrics

Briefly connect the graph-construction report to the benchmark metrics:

$$
\alpha^*
=
\arg\min_{\alpha>0}
\sum_{i<j}
\left(\alpha\,d_G(i,j)-d_{\mathrm{ref}}(i,j)\right)^2.
$$

Then the relative RMS error is

$$
\left(
\frac{
\sum_{i<j}(\alpha^* d_G(i,j)-d_{\mathrm{ref}}(i,j))^2
}{
\sum_{i<j}d_{\mathrm{ref}}(i,j)^2
}
\right)^{1/2}.
$$

This makes the report self-contained for readers who then look at the HTML
benchmark.

## Suggested Expansion Task List For The Next Agent

1. Rename the report source/output or add a clear plan for the rename.
2. Replace the date stamp under the title with `Build timestamp`.
3. Rewrite the introduction around geodesic distance approximation.
4. Add a formal problem statement for data-to-graph geodesic reconstruction.
5. Reframe PHATE as an affinity/diffusion method whose support motivates sKNN
   and adaptive-radius comparisons.
6. Add fixed-radius graph construction with formulas, algorithm, use cases, and
   references.
7. Add adaptive-radius graph construction with max/min rules, formulas,
   algorithm, use cases, and references.
8. Add detailed pruning algorithms, especially local geodesic pruning, with at
   least one explanatory figure.
9. Add detailed MST repair algorithm with at least one explanatory figure.
10. Add examples: figure eight, crossing line segments, tree/branching data,
    nonuniform noisy circle, and quadratic graph surfaces.
11. Add a historical perspective section with verified references.
12. Add a terminology table separating supports, affinities, lengths, Markov
    operators, and graph-geodesic metrics.
13. Add a computational considerations section, including the current local
    pruning performance issue.
14. Rebuild the LaTeX report and check that every figure included in the TeX has
    a current source script or clear generation path.

## Current Open Implementation Issues Related To The Report

These do not all need to be solved in the LaTeX report, but the report should
not hide them.

1. Local geodesic pruning is slow in current benchmark runs. It should be
   profiled and optimized before large parameter sweeps are trusted.
2. The fastest iKNN global pruning and the cross-family local geodesic pruning
   are different algorithms; comparisons must label them explicitly.
3. PHATE edge weights are affinities, not lengths. Any PHATE-as-graph-geodesic
   experiment must state the affinity-to-length transformation.
4. Fixed-radius and adaptive-radius graphs now have pruning and MST repair
   support in `gflow`, so the report should include them as first-class graph
   constructors.
5. The current benchmark runner/report work is still evolving. Treat the
   benchmark HTML as an active exploration, not a finalized result.

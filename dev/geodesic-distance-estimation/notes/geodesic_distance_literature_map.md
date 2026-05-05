# Geodesic Distance Literature Map

Date: 2026-05-04

This memo maps the broader, mostly density-unaware literature on estimating
intrinsic or geodesic distance from finite data. It is intentionally separate
from the density-aware lane. Density enters here mainly through sampling,
neighbor radii, local scales, and graph connectivity; it is not the primary
edge-cost target unless stated explicitly.

The central lesson for `gflow` is simple but important: graph topology and edge
weights are different design decisions. A shortest-path distance can only be as
good as the graph on which paths are allowed to travel. Reweighting a bad
shortcut may help if the weight detects that it is unsupported, but a cleaner
workflow first asks whether the neighbor graph should have admitted the edge at
all.

## Executive Map

The classical graph-geodesic estimator is the Isomap distance:

$$
d_G(i,j)
= \min_{\pi:i\to j}
\sum_{(u,v)\in\pi} \|x_u-x_v\|.
$$

The graph `G` is usually a `k`-nearest-neighbor graph or an epsilon/radius
graph, and the edge weights are local ambient distances. If the data densely
sample a smooth manifold, and the neighborhood scale shrinks slowly enough,
then these graph shortest paths can converge to the manifold geodesic distance.
This is the borrowed backbone behind many practical data-geodesic ideas,
starting with Isomap ([Tenenbaum, de Silva, and Langford 2000](https://doi.org/10.1126/science.290.5500.2319))
and its theoretical analysis ([Bernstein et al. 2000/2001](https://www.researchgate.net/publication/2466340_Graph_Approximations_to_Geodesics_on_Embedded_Manifolds)).

For `gflow`, the competing families should be evaluated fairly:

- ordinary directed/symmetrized kNN shortest paths;
- mutual kNN and shared-neighbor graphs;
- adaptive-k, local-scale, radius, and hybrid graphs;
- pruning by geometric empty-region rules, local alternative paths, or topology;
- connectivity repair by MST, component bridges, landmarks, or scaffold edges;
- diffusion, resistance, PHATE, UMAP, and pseudotime distances as comparators,
  but not automatically as geodesic distances.

No method family is the answer by definition. CMST/MST completion is best viewed
as one connected scaffold strategy. It can be useful when finite graph distances,
edge attribution, and bridge diagnostics matter, but it can also impose tree
artifacts or false bridges. A well-tuned kNN or mKNN graph may be the cleaner
geodesic estimator in some regimes.

## Basic Objects and Vocabulary

Let

$$
X_n = \{x_1,\ldots,x_n\}\subset \mathbb{R}^D
$$

be sampled near a lower-dimensional object `M`, often imagined as a compact
`d`-dimensional Riemannian manifold, stratified manifold, tube, curve network,
or branching biological state space.

A finite weighted graph is

$$
G=(V,E,w), \qquad V=\{1,\ldots,n\}, \qquad w_{ij}\ge 0.
$$

The graph shortest-path metric is

$$
d_G(i,j)
= \min_{i=v_0,\ldots,v_m=j}
\sum_{r=0}^{m-1}w_{v_rv_{r+1}},
$$

with `d_G(i,j)=infinity` if `i` and `j` are disconnected. On each connected
component with nonnegative weights and no zero-weight identifications, this is a
metric. With zero weights or directed/asymmetric reachability, it may become a
pseudometric or a directed quasi-distance.

The continuous target, when it exists, is the Riemannian geodesic distance

$$
d_M(x,y)
= \inf_{\gamma:x\to y}
\int_0^1
\sqrt{g_{\gamma(t)}(\gamma'(t),\gamma'(t))}\,dt.
$$

The data problem is not just "estimate this equation." Biological data add
finite-sample gaps, high-dimensional nuisance variation, batch effects, count
noise, compositional transforms, branches, loops, boundaries, outliers, and
scientific ambiguity about whether a missing region is an unobserved transition
or a true forbidden state.

## Historical Spine

1. **Computational geometry and graph shortest paths.** Dijkstra-style graph
   shortest paths are older than manifold learning, but proximity graphs such as
   Delaunay triangulations, Gabriel graphs, relative neighborhood graphs, beta
   skeletons, and geometric spanners supplied a vocabulary for sparse graphs
   that preserve local geometry or Euclidean stretch.

2. **Isomap made graph geodesics central to manifold learning.** Isomap builds a
   kNN or epsilon graph, weights edges by local Euclidean distances, computes
   all-pairs shortest paths, and then uses classical MDS on those approximate
   geodesic distances ([Science 2000](https://doi.org/10.1126/science.290.5500.2319)).

3. **Theory clarified the asymptotic promise.** Bernstein, de Silva, Langford,
   and Tenenbaum analyzed when graph distances approximate manifold geodesics
   under dense enough sampling, bounded curvature, and appropriate neighborhood
   scale ([graph approximation report](https://www.researchgate.net/publication/2466340_Graph_Approximations_to_Geodesics_on_Embedded_Manifolds)).
   Later random geometric graph theory made connectivity and shortest-path
   behavior more precise ([Penrose 2003](https://academic.oup.com/book/9064)).

4. **Failure modes were recognized almost immediately.** Too small a
   neighborhood disconnects the graph or produces tree-like detours. Too large a
   neighborhood creates "short circuits" across folds, holes, branches, gaps,
   or noisy off-manifold points. Balasubramanian and Schwartz highlighted
   topological instability in early Isomap discussions; modern software still
   exposes the same parameter tradeoff.

5. **Scalable variants replaced all-pairs geodesics with landmarks.** Landmark
   Isomap / Landmark MDS embeds a subset of landmarks and places the remaining
   points from landmark distances, reducing memory and time at the cost of
   landmark-selection sensitivity ([de Silva and Tenenbaum 2004](https://graphics.stanford.edu/courses/cs468-05-winter/Papers/Landmarks/Silva_landmarks5.pdf)).

6. **Out-of-sample extension became a separate issue.** A new point must be
   attached to the training graph, distances from it to training points must be
   estimated, and the embedding or downstream distance must be extended without
   rebuilding everything. Bengio et al. framed Isomap, MDS, LLE, eigenmaps, and
   spectral clustering as eigenfunction extensions ([NeurIPS 2003](https://papers.nips.cc/paper_files/paper/2003/file/cf05968255451bdefe3c5bc64d550517-Paper.pdf)).

7. **Applied biology mostly adopted graphs, not pure Isomap distances.**
   Single-cell workflows widely use kNN graphs, MSTs, principal graphs,
   diffusion maps, PAGA, and pseudotime graph distances. Microbiome and
   metagenomics workflows more often use beta-diversity, Aitchison/UniFrac,
   ordination, PHATE, UMAP, and network/TDA methods than explicit kNN
   shortest-path distance matrices.

## Construct 1: Classical Isomap and kNN/Epsilon Graph Geodesics

### Motivation and Problem Setting

Isomap was introduced for nonlinear dimensionality reduction when pairwise
ambient distances are misleading because the data lie on or near a curved
low-dimensional manifold. The problem setting is geometric rather than
probabilistic: preserve intrinsic distances on the manifold, then find a
low-dimensional Euclidean representation that best matches those distances.

### Definition

Choose a neighborhood graph:

$$
E_k^{\mathrm{dir}}
= \{(i,j): x_j \in N_k(x_i)\},
$$

or

$$
E_{\varepsilon}
= \{(i,j): \|x_i-x_j\|\le \varepsilon\}.
$$

For undirected Isomap, directed kNN is usually symmetrized by union:

$$
E_k^{\cup}
= \{\{i,j\}: j\in N_k(i)\ \mathrm{or}\ i\in N_k(j)\}.
$$

Weights are usually

$$
w_{ij}=\|x_i-x_j\|.
$$

The graph-geodesic estimate is

$$
\widehat d_M(x_i,x_j)=d_G(i,j).
$$

Isomap then embeds with classical scaling:

$$
B = -\frac{1}{2} H D_G^{\circ 2} H,
\qquad
H=I-\frac{1}{n}\mathbf{1}\mathbf{1}^T,
$$

and takes the leading eigenvectors of `B`.

Metric type: weighted graph metric on connected components; continuous
geodesic estimator under assumptions; embedding distance only after MDS.

### Parameters and Practical Defaults

The main parameter is neighborhood scale. `k` controls degree; `epsilon`
controls a global radius. There is no universal default. Current
scikit-learn `Isomap` exposes `n_neighbors=5`, `radius=None`,
`path_method='auto'`, and `n_components=2` by default
([scikit-learn Isomap docs](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html)).
Single-cell software often uses larger graph neighborhoods for other purposes:
Scanpy defaults to `n_neighbors=15` for `scanpy.pp.neighbors`
([Scanpy docs](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html)),
while Seurat `FindNeighbors` defaults to `k.param=20` and often computes a
shared-nearest-neighbor graph ([Seurat docs](https://satijalab.org/seurat/reference/findneighbors)).
Those are adoption signals, not geodesic-optimality guarantees.

Common selection heuristics:

- sweep `k` or `epsilon`;
- require one large component but inspect bridges;
- use residual variance or embedding stress for Isomap;
- check neighborhood overlap against known latent coordinates in simulations;
- check false-shortcut edges across folds, holes, branches, or gaps;
- prefer stability across a range rather than a single lucky parameter.

### Asymptotic and Finite-Data Behavior

For a smooth compact manifold with enough samples and an appropriate shrinking
radius, local Euclidean edge lengths approximate local manifold arc lengths, and
shortest paths on the graph approximate continuous geodesics. The intuition is:
each true geodesic can be covered by small sample-to-sample steps, and each graph
edge is short enough that its ambient chord is close to its intrinsic local
distance. The theory requires sampling density bounded away from zero on the
support, a scale larger than the covering radius, and local geometry controlled
by curvature/reach.

Finite data violate these assumptions in exactly the ways `gflow` cares about:
uneven sampling, high-dimensional noise, off-manifold thickness, true gaps,
branch points, and boundaries.

### Strengths

- Directly estimates an interpretable graph geodesic.
- Simple to implement and inspect.
- Edge/path attribution is straightforward.
- Works well on clean S-curves, spirals, circles, and tubes when `k` is tuned.
- Separates local metric choice from downstream regression or embedding.

### Failure Modes

- **Shortcuts:** an edge across a fold, crescent, annulus, branch, gap, or
  figure-eight self-approach can shorten many paths.
- **Disconnectedness:** small `k` or small `epsilon` fragments the graph,
  yielding infinite distances.
- **Tree artifacts:** too sparse a graph approximates a traversal tree rather
  than the manifold.
- **Boundary effects:** boundary points have fewer neighbors, can be pulled into
  long detours, or can be connected across periodic boundaries incorrectly.
- **Holes and voids:** a missing sector can be treated as either a real hole or
  a sample gap; the graph does not know which prior is intended.
- **Noise:** off-manifold radial or nuisance variation can create ambient chords
  that are not intrinsic adjacencies.
- **Nonuniform sampling:** fixed `k` uses larger radii in sparse regions and
  smaller radii in dense regions; fixed radius does the reverse with degree.

### gflow Diagnostics

For every graph-geodesic estimator, `gflow` should report:

- component count and largest component fraction;
- edge-length distribution by graph family;
- degree distribution and high-degree hubs;
- cycle rank and whether a circle-like structure is path-like or loop-like;
- false-shortcut rate on synthetic truth;
- long bridge edges and their path usage;
- sensitivity over `k`, `epsilon`, preprocessing, and PCA dimension;
- shortest-path attribution for influential downstream regression pairs.

Relevance: the classical estimator is the baseline to beat. On noisy circles,
tubes, gaps, branches, and compositional gradients, it makes the right question
visible: which topology plus which weights recover the latent distances useful
for gflow regression and gradient-flow summaries?

## Construct 2: Landmark Isomap, Scalable Isomap, and Out-of-Sample Extension

### Landmark Isomap

Landmark Isomap chooses a subset `L` of `m << n` landmark points, computes
shortest-path distances from landmarks to all samples, embeds the landmarks by
MDS, and places non-landmarks from their landmark-distance profiles. Landmark
MDS uses the fact that distances to enough well-spread landmarks can locate a
point in the embedding.

Let `D_LL` be landmark-to-landmark graph distances and `D_LX` be
landmark-to-all graph distances. Classical MDS is applied to `D_LL`; nonlandmark
points are projected by a Nyström-like / triangulation formula from `D_LX`.

Metric type: the underlying distances are graph geodesics; the scalable output
is an approximate embedding and approximate coordinate extension.

Strengths:

- avoids storing or eigendecomposing the full `n x n` matrix;
- can support large biological data sets;
- gives natural landmark/witness diagnostics.

Failure modes:

- poor landmark coverage misses branches, rare states, holes, or sparse
  transitions;
- shortest-path errors from a bad graph are inherited;
- landmark distances can smooth over local defects.

gflow relevance: landmarks are useful for large `n`, but rare biological states
and compositional tails should be actively sampled as landmarks rather than
left to random chance.

### Out-of-Sample Extension

For a new sample `x_*`, one common Isomap extension attaches it to the training
graph through its nearest training neighbors:

$$
\widehat d_G(*,j)
= \min_{i\in N_k(*)}
\|x_* - x_i\| + d_G(i,j).
$$

The new squared-distance vector is then centered and projected into the learned
MDS/kernel-PCA coordinates. This is the practical logic exposed in
scikit-learn's `Isomap.transform`, which links query points to the training
geodesic graph and projects the resulting kernel onto the training embedding
([scikit-learn Isomap docs](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html)).
Bengio et al. give the broader eigenfunction-extension framing
([NeurIPS 2003](https://papers.nips.cc/paper_files/paper/2003/file/cf05968255451bdefe3c5bc64d550517-Paper.pdf)).

For `gflow`, out-of-sample distance extension matters more than embedding
extension. New samples should report which training edges or bridge components
control their geodesic attachment; otherwise, a single off-manifold query can
silently inherit a misleading downstream gradient.

## Construct 3: kNN Graph Variants and Local Topology Choices

The neighborhood graph is the main finite-data object. Edge weights can be
changed later, but an absent edge cannot be used and an admitted shortcut can
dominate all-pairs paths.

### Directed kNN

Definition:

$$
(i,j)\in E_k^{\mathrm{dir}}
\quad \Longleftrightarrow \quad
x_j \in N_k(x_i).
$$

Metric type: directed graph distance if paths respect edge direction. This is
not symmetric:

$$
d_G(i,j)\neq d_G(j,i)
$$

in general.

Use case: directed kNN is natural when modeling flow, velocity, or asymmetric
reachability. For ordinary geodesic distance, it is usually an intermediate
construction before symmetrization.

Failure mode: directed reachability can create one-way finite distances and
break ordinary metric assumptions in regression unless explicitly intended.

### Symmetrized kNN

Union symmetrization:

$$
\{i,j\}\in E_k^{\cup}
\Longleftrightarrow
j\in N_k(i)\ \mathrm{or}\ i\in N_k(j).
$$

This keeps edges if either endpoint selected the other. It is more connected
and can preserve sparse-region reachability, but admits asymmetric hub edges.

Intersection / mutual kNN:

$$
\{i,j\}\in E_k^{\cap}
\Longleftrightarrow
j\in N_k(i)\ \mathrm{and}\ i\in N_k(j).
$$

Mutual kNN is conservative. It often reduces hubness and false long-radius edges
from sparse points into dense regions, but it fragments easily.

Metric type: graph construction; becomes a graph metric after weights and
shortest paths are defined.

gflow relevance: this is a direct candidate family. The noisy-circle note
already observed a plausible division: iKNN/symmetrized kNN can recover loops
but may become dense; mKNN can be geometrically cleaner but disconnected.

### Shared Nearest Neighbor Graphs

Shared-neighbor methods weight an edge by overlap between neighbor sets:

$$
s_{ij}
= \frac{|N_k(i)\cap N_k(j)|}{|N_k(i)\cup N_k(j)|}
$$

or related counts. Seurat's SNN graph uses neighborhood overlap and prunes low
Jaccard-weight edges by default ([Seurat docs](https://satijalab.org/seurat/reference/findneighbors)).

Metric type: similarity graph; not a geodesic distance unless transformed into
edge lengths and used in shortest paths.

Strength: robust to local high-dimensional noise and hubness.

Failure: shared-neighbor weights can emphasize cluster membership over
intrinsic arclength, which may be good for clustering and bad for distance
estimation.

### Radius / Epsilon Graphs

Definition:

$$
\{i,j\}\in E_{\epsilon}
\Longleftrightarrow
\|x_i-x_j\|\le \epsilon.
$$

Metric type: graph construction; weighted shortest paths estimate geodesics
under classical random geometric graph assumptions.

Strength: simple and theoretically clean for bounded-density samples.

Failure: a global radius gives high degree in dense regions and isolated points
in sparse regions. In biological data with uneven sampling, this can be worse
than fixed `k`.

### Hybrid kNN/Radius Graphs

Hybrid admission rules combine local degree and maximum edge length, for
example:

$$
\{i,j\}\in E
\Longleftrightarrow
\big(j\in N_k(i)\ \mathrm{or}\ i\in N_k(j)\big)
\ \mathrm{and}\
\|x_i-x_j\|\le r_{\max}.
$$

or

$$
\|x_i-x_j\|\le c\max(r_i,r_j),
\qquad
r_i = \|x_i-x_{i(k)}\|.
$$

Metric type: graph construction. The second rule is local-scale or adaptive
radius rather than global radius.

Strength: prevents some long sparse-region kNN shortcuts while preserving
minimum degree.

Failure: local scale must be capped or robustified around outliers; otherwise
outlier neighborhoods can authorize exactly the bridges one wants to flag.

### Locally Scaled Graphs

Self-tuning spectral clustering defines a local scale `sigma_i`, often the
distance to the `k`th neighbor, and uses affinities like

$$
A_{ij}
= \exp\left(
-\frac{\|x_i-x_j\|^2}{\sigma_i\sigma_j}
\right).
$$

The original motivation was multi-scale clustering and clutter robustness
([Zelnik-Manor and Perona 2004](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering)).

Metric type: kernel/spectral graph. It is not a shortest-path geodesic unless
one deliberately converts affinities into path costs, for example
`w_ij=-log A_ij` or a locally scaled edge length.

gflow relevance: local scales are natural diagnostics and edge-admission tools.
They should be kept conceptually separate from density-aware `ell/rho^alpha`
edge costs, because local scaling can either exploit density, normalize it away,
or merely adapt the graph to uneven sampling.

### Adaptive Neighborhoods

Adaptive-neighborhood methods try to infer the right local graph degree from
the data rather than fixing `k`. IAN starts from a conservative Gabriel graph,
then iteratively sparsifies/adapts neighborhoods using a global optimization and
outlier checks, with applications to nonlinear dimension reduction, geodesic
computation, and dimension estimation
([Dyballa and Zucker 2023](https://direct.mit.edu/neco/article/35/3/453/114733/IAN-Iterated-Adaptive-Neighborhoods-for-Manifold),
[arXiv](https://arxiv.org/abs/2208.09123)).

Metric type: graph construction; final graph distances may be geodesic
estimators.

gflow relevance: adaptive `k` is likely important for uneven noisy circles,
sparse transitions, branches, and compositional data where fixed-neighborhood
rules confound sampling density with biology.

## Construct 4: Connectivity Repair and MST Scaffolds

### MST and CMST as Scaffolds

The Euclidean MST on `X_n` is

$$
T = \arg\min_{\text{spanning trees }T}
\sum_{\{i,j\}\in T}\|x_i-x_j\|.
$$

It guarantees connectedness with `n-1` edges. CMST-style completion starts from
an MST and adds local edges under a threshold or rule:

$$
E_{\mathrm{CMST}}
= E_{\mathrm{MST}}\cup E_{\mathrm{local}}.
$$

Metric type: connected graph shortest-path metric after weights are assigned.

Strengths:

- finite distances for all pairs;
- long forced bridges are inspectable;
- sparse and interpretable;
- useful when downstream methods require one connected graph.

Failures:

- an MST is not a manifold skeleton in general;
- it can break loops into paths;
- it can force false bridges between true components;
- adding MST edges to an already good kNN/mKNN graph can shorten paths through
  tree artifacts.

gflow position: MST/CMST is a candidate scaffold and diagnostic, not a privileged
source of intrinsic geometry.

### Component Repair for Conservative Graphs

A conservative graph such as mKNN may fragment. Instead of unioning with the
full MST, one can connect components by selected bridges:

$$
e^*_{ab}
= \arg\min_{i\in C_a,\ j\in C_b}
\|x_i-x_j\|,
$$

possibly subject to local-scale or geometric plausibility constraints:

$$
\|x_i-x_j\|\le c\max(r_i,r_j),
\qquad
\mathrm{or}
\qquad
\frac{\|x_i-x_j\|}{\mathrm{median\ local\ edge}}\le c.
$$

Repair can also use a component-level MST, where vertices are components and
candidate edges are nearest inter-component bridges.

Metric type: graph repair before shortest paths.

Strength: preserves the low-shortcut behavior of mKNN while making distances
finite.

Failure: if the component split represents a true biological separation,
repairing it erases meaningful disconnectedness.

Diagnostics: number of repaired components, bridge length quantiles, bridge
path usage, downstream sensitivity to removing or inflating bridge weights, and
whether bridges cross known latent gaps in simulations.

## Construct 5: Pruning and Shortcut Control

Pruning asks whether an admitted edge should remain. This is separate from
edge reweighting.

### Local Alternative-Path Pruning

An edge `e={i,j}` can be removed if the graph without it still has a short
alternative path:

$$
d_{G\setminus e}(i,j)
\le \tau w_{ij}.
$$

For `tau` near 1, this removes redundant edges while preserving local path
lengths. A greedy spanner uses the dual rule: add an edge only if no existing
path is short enough.

Metric type: topology modification before graph shortest paths.

Strength: reduces dense kNN/iKNN shortcuts and can preserve cycle geometry with
fewer edges.

Failure: if the alternate path itself contains a shortcut, the rule ratifies
bad topology. If `tau` is too strict, it leaves a near-tree; if too loose, it
may overprune legitimate local supports.

gflow diagnostics: removed-edge attribution, stretch distribution, component
check after each pruning level, and false-shortcut rate on known latent cases.

### Relative Neighborhood Graph

The relative neighborhood graph connects `i` and `j` when there is no third
point closer to both endpoints than they are to each other:

$$
\{i,j\}\in E_{\mathrm{RNG}}
\Longleftrightarrow
\nexists k:
\max\{\|x_i-x_k\|,\|x_j-x_k\|\}<\|x_i-x_j\|.
$$

Toussaint introduced RNGs as a proximity structure for finite planar sets
([Toussaint 1980 record](https://cir.nii.ac.jp/crid/1361981468897626368)).
RNGs are subgraphs of Delaunay triangulations in the plane and contain the
Euclidean MST, hence are connected for planar Euclidean point sets in general
position ([overview](https://en.wikipedia.org/wiki/Relative_neighborhood_graph)).

Metric type: sparse proximity graph; shortest paths on it are graph geodesics
but not automatically manifold-geodesic estimators.

Strength: removes many edges with obvious local alternatives.

Failure: can be too sparse and has no direct adaptation to high-dimensional
noise, compositional metrics, or manifold folds.

### Gabriel Graph

The Gabriel graph connects `i` and `j` if the closed ball with diameter
`x_i x_j` contains no other sample point:

$$
\{i,j\}\in E_{\mathrm{Gabriel}}
\Longleftrightarrow
\nexists k:
\left\|x_k-\frac{x_i+x_j}{2}\right\|
< \frac{1}{2}\|x_i-x_j\|.
$$

Gabriel and Sokal introduced the graph in geographic variation analysis; it is
also a beta skeleton and a subgraph of the Delaunay triangulation in the plane
([overview](https://en.wikipedia.org/wiki/Gabriel_graph)).

Metric type: empty-ball proximity graph.

Strength: scale-free local criterion; useful conservative starting graph.

Failure: empty-ball tests in high dimensions are fragile because points are
sparse and distances concentrate; a Gabriel edge can still be a manifold
shortcut across a fold if the ambient ball is empty.

### Delaunay-Based Graphs

Delaunay triangulation connects points whose Voronoi cells touch. In low
dimensions it is a rich proximity supergraph. Many proximity graphs, including
Gabriel and RNG, are subgraphs of Delaunay in the plane.

Metric type: proximity graph; weighted shortest paths can be used as graph
geodesics.

Strength: good computational geometry foundation and spanner properties in
Euclidean planar settings.

Failure: high-dimensional Delaunay is often infeasible and may be too dense;
ambient Delaunay does not know the manifold topology.

### Beta Skeletons

Beta skeletons use an empty lens or circle region controlled by `beta`; Gabriel
is the `beta=1` case and the relative neighborhood graph is a common `beta=2`
lune-based case. Increasing `beta` usually makes the graph more conservative
([beta skeleton overview](https://en.wikipedia.org/wiki/Beta_skeleton)).

Metric type: geometric pruning family.

Strength: tunable shortcut control.

Failure: planar Euclidean intuition may not transfer to high-dimensional
biological feature spaces without preprocessing.

### Geometric Spanners

A weighted graph is a `t`-spanner of a metric space if

$$
d_G(i,j)\le t\,d(i,j)
$$

for all pairs. Greedy spanners add an edge only when the current graph lacks a
path of length at most `t` times the direct distance. Theta/Yao graphs use cones
around each point to keep a bounded number of directional neighbors. Geometric
spanners were developed for sparse networks that preserve Euclidean distances
([geometric spanner overview](https://en.wikipedia.org/wiki/Geometric_spanner)).

Metric type: sparse graph with bounded stretch relative to a chosen base
metric, usually Euclidean.

Important caution: Euclidean spanner guarantees are not manifold-geodesic
guarantees. A graph that well preserves ambient Euclidean distances may be
exactly the graph that preserves harmful chords.

gflow relevance: spanner logic is valuable when applied to a local graph or
latent/oracle metric, not blindly to ambient Euclidean all-pairs distances.

### Witness and Landmark Graphs

Witness graphs and witness complexes use a smaller landmark set `L` and allow
nonlandmark points to witness relationships among landmarks. A simple graph
version connects landmarks `a,b` when some witness `x` has both among its near
landmarks:

$$
\{a,b\}\in E
\Longleftrightarrow
\exists x:
a,b\in N_m^L(x).
$$

Metric type: sparse topology surrogate; graph geodesic only after edge weights
are defined.

Strength: scalable and useful for rare-state coverage if landmark selection is
designed.

Failure: compresses cell/sample-level paths; branch tips and small holes can be
lost if no landmark/witness supports them.

### Topology-Aware Pruning

Topology-aware pruning uses homology, cycle persistence, branch consistency, or
cluster abstraction to decide whether edges preserve or destroy intended
structure. In single-cell analysis, PAGA estimates connectivity of manifold
partitions from a symmetrized kNN-like graph and provides a topology-preserving
abstracted graph ([PAGA Genome Biology 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x)).

Metric type: usually graph abstraction or topology diagnostic; not necessarily
a sample-level geodesic metric.

gflow relevance: topology-aware pruning could flag false circle chords,
figure-eight crossovers, branch-to-branch shortcuts, and gap bridges before
distance estimation.

## Construct 6: Theoretical Foundations

### Convergence of Graph Shortest Paths

The classical convergence story requires:

- samples lie on or close to a compact smooth manifold;
- sampling density is bounded above and away from zero on the support;
- the neighborhood radius `epsilon_n` tends to zero but remains above the
  sample covering radius;
- local chords approximate local manifold distances;
- curvature, reach, and boundary do not destroy local Euclidean approximations.

In this regime, graph shortest paths over epsilon or local graphs converge to
manifold geodesics. Bernstein et al. provided an early proof for Isomap-style
claims. Later work generalized shortest-path approximation on sampled subsets,
including curvature-constrained variants ([Arias-Castro and Le Gouic 2017](https://arxiv.org/abs/1706.09441)).

### Random Geometric Graph Connectivity

A random geometric graph connects points within radius `r_n`. Connectivity
thresholds scale like

$$
r_n \asymp \left(\frac{\log n}{n}\right)^{1/d}
$$

up to constants and boundary/density effects. Penrose's monograph is the
standard reference for rigorous random geometric graph theory
([Penrose 2003](https://academic.oup.com/book/9064)).

For kNN graphs, the right order for connectivity is also logarithmic in `n` in
classical uniform settings. Balister, Bollobas, Sarkar, Walters, and related
work sharpen constants in planar random kNN graphs; a recent survey-style
introduction notes that `k ~ log n` is the right order for connectivity
([PMC review/introduction](https://pmc.ncbi.nlm.nih.gov/articles/PMC7480951/)).

Practical implication: fixed small `k` is not an asymptotic connectivity rule.
It is a finite-data modeling choice. In single-cell workflows, fixed defaults
like 15 or 20 are pragmatic, not theoretical guarantees.

### Density and Edge Weight Regimes

Alamgir and von Luxburg analyzed shortest-path limits in random kNN graphs and
showed that unweighted kNN shortest paths can induce an undesirable
density-dependent geometry, preferring low-density regions in a way that is
opposite to many machine-learning intuitions
([ICML 2012 PDF](https://www.tml.cs.uni-tuebingen.de/team/luxburg/publications/AlamgirLuxburg2012.pdf),
[arXiv](https://arxiv.org/abs/1206.6381)).

This is directly relevant to `gflow`: hop count is not a neutral geodesic
distance. Euclidean edge weights are usually safer as a default for geometric
distance, while density-aware or power-weighted alternatives should be named as
separate metric choices.

Power-weighted shortest paths

$$
d_{p}(i,j)
= \min_{\pi:i\to j}
\sum_{(u,v)\in\pi}\|x_u-x_v\|^p
$$

have a density-weighted continuum limit for `p>1` and are covered in the
density-aware report. Here they matter mainly as a warning: changing edge
weights changes the target metric, not just numerical stability.

### Intrinsic Dimension, Reach, Curvature, Boundary

Intrinsic dimension controls covering rates and connectivity thresholds.
Reach controls how far one can move before local Euclidean neighborhoods become
ambiguous across folds or self-approaches. Curvature controls chord-vs-arc
error:

$$
\mathrm{arc}(i,j) - \|x_i-x_j\|
= O(\kappa^2\|x_i-x_j\|^3)
$$

for small edges on a smooth curve with curvature `kappa`. Boundaries reduce
neighbor support and make symmetric local neighborhoods false. Branch points and
stratified intersections violate the smooth-manifold assumptions entirely.

Practical implication: a method can be asymptotically consistent on smooth
manifolds and still fail on the exact biological structures `gflow` wants:
branches, basins, loops, gaps, tubes, and compositional sparse count clouds.

### kNN Spanner Results

Spanner results say when a sparse graph approximately preserves shortest-path
distances of another graph or metric. Little, McKenzie, and Murphy showed that
kNN graphs can approximate complete-graph power-weighted shortest-path
distances under suitable sampling and `k` scaling
([SIAM MDS 2022](https://epubs.siam.org/doi/10.1137/20M1386657)).

For density-unaware `p=1` Euclidean-weighted graph geodesics, spanner intuition
supports using sparse local graphs for computation, but only after the target
metric and admissible neighborhood topology are specified.

## Practical Adoption Review

### Bioinformatics and Computational Biology

General bioinformatics has used Isomap and manifold learning for visualization,
feature extraction, image and microarray analysis, and disease prediction, but
explicit all-pairs graph-geodesic distances are less commonly treated as the
final scientific object. For example, microarray feature-selection work describes
Isomap's kNN graph and Dijkstra/Floyd-Warshall geodesic matrix as part of a
manifold-learning pipeline ([PMC example](https://pmc.ncbi.nlm.nih.gov/articles/PMC3940899/)).
Metagenomic deep-learning work has used Isomap as a 2D representation method,
again defining geodesic distances by shortest paths on a kNN graph
([Cell Reports Methods 2022](https://www.sciencedirect.com/science/article/pii/S2666389922002987)).

Transferable lesson: practical users often consume the embedding, not the
distance matrix. `gflow` is more distance-forward, so it needs stronger
diagnostics than typical visualization papers provide.

### Microbiome, Metagenomics, and Compositional Count Data

The dominant microbiome distance ecosystem is not Isomap-style graph geodesics.
It is beta-diversity and compositional geometry: Bray-Curtis, Jaccard,
Jensen-Shannon, weighted/unweighted UniFrac, generalized UniFrac, robust
Aitchison, DEICODE, ordination, and PERMANOVA-style workflows.

Evidence for routine kNN shortest-path geodesics in microbiome analysis is
thin. Isomap appears as a representation option in some metagenomic machine
learning workflows, and PHATE/tmap/diffusion/TDA approaches show that graph
geometry is acceptable in the domain, but explicit geodesic distance matrices
are not standard microbiome outputs.

gflow implications:

- use Aitchison/CLR/robust Aitchison or other compositional preprocessing before
  Euclidean neighbor search;
- compare to UniFrac and beta-diversity baselines rather than only Euclidean
  kNN;
- treat zeros, rare taxa, library size, and batch effects as topology risks;
- inspect whether graph paths follow ecological gradients or merely sequencing
  depth and compositional sparsity.

### Single-Cell Genomics and Trajectory Inference

Single-cell analysis is the strongest adoption area for graph paths, but the
term "geodesic" often refers to distance along an inferred trajectory or
principal graph rather than all-pairs Isomap distance on the original kNN graph.

Examples:

- Monocle 2 learns a principal tree and calculates pseudotime as geodesic
  distance from a root after projecting cells to the principal graph
  ([Monocle 2 / DDRTree PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5764547/)).
- Reviews of trajectory methods summarize many MST, principal-graph, and kNN
  graph workflows where pseudotime is a distance along a learned graph or curve
  ([scTEP BMC Bioinformatics 2023](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05179-2),
  [Slingshot comparison table](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0/tables/1)).
- Wanderlust uses ensembles of kNN graphs and shortest-path-like temporal
  distances; Wishbone and DPT use diffusion geometry; PAGA abstracts a kNN-like
  graph into a topology-preserving partition graph; Palantir uses diffusion maps
  and shortest paths in diffusion-component space in parts of its workflow
  ([VIA review/background](https://pmc.ncbi.nlm.nih.gov/articles/PMC8452770/)).

This is a strong adoption signal for graph-based geometry but a caution for
terminology. A root-to-cell pseudotime on a principal tree is a graph geodesic
on that tree; it is not the same as all-pairs manifold geodesic distance on the
sample cloud.

gflow implications:

- principal-graph distances should be comparators for branches and basins;
- all-pairs graph geodesics may capture more geometry than a fitted tree, but
  can also overfit noise;
- pseudotime methods provide practical diagnostics: root sensitivity, branch
  assignment, projection residuals, and lineage stability.

### Spatial Transcriptomics

Spatial transcriptomics often uses graph distances in physical space or
segmentation graphs. A recent point-cloud segmentation framework computes
geodesic distance as graph shortest path between RNA nodes and cell centroids
using Dijkstra's algorithm ([PMC spatial transcriptomics example](https://pmc.ncbi.nlm.nih.gov/articles/PMC11227573/)).

For cell-state manifolds in expression space, adoption appears more diffusion,
neighbor graph, and graph neural network oriented than explicit Isomap-style
geodesic-distance oriented.

gflow implications: distinguish physical-space geodesics, expression-space
geodesics, and tissue-neighborhood graph distances. They answer different
questions.

### Ecology and Community Composition

Ecology uses distance along environmental gradients, resistance surfaces, and
ordination heavily. Graph shortest paths are natural in landscape connectivity,
but the graph usually represents geography, dispersal cost, or habitat
resistance rather than a sample cloud in high-dimensional composition space.

Transferable lesson: when a path distance is scientifically useful, the cost
surface is interpretable. For `gflow` microbiome gradients, the analogous cost
surface might be compositional dissimilarity, ecological transition support, or
sample-density-aware plausibility; it should not be left implicit.

### Medical Imaging, Morphology, Shape, and Phenotyping

Shape analysis uses geodesics on meshes, heat kernels, diffusion descriptors,
and intrinsic surface distances. Mesh geodesic algorithms differ from sample
kNN graph geodesics because the surface connectivity is often observed or
reconstructed. Graph-based discrete geodesic methods on meshes can carry
accuracy guarantees when vertices are well distributed
([discrete geodesic graph example](https://www.sciencedirect.com/science/article/abs/pii/S0167839617300390)).

Transferable lesson: when topology is known from a mesh, shortest paths are
much safer. In `gflow`, topology is inferred and should be diagnosed as the
main uncertain object.

### Non-Biological Domains

Isomap-style geodesics are widely used in computer vision, image manifolds,
robotics/path planning, remote sensing, speech/audio manifolds, and general
nonlinear dimensionality reduction. Non-biological lessons transfer cleanly:

- nearest-neighbor scale selection controls success;
- holes and folds punish shortcuts;
- landmarks are necessary for scale but can miss rare modes;
- high-dimensional approximate nearest-neighbor errors can change graph
  topology;
- embedding quality does not prove distance quality.

## Related but Distinct Distances

### Diffusion Distance

Diffusion maps build a kernel graph, normalize it into a Markov transition
operator `P`, and compare transition distributions after `t` steps:

$$
D_t^2(i,j)
= \sum_k
\frac{(P^t(i,k)-P^t(j,k))^2}{\pi_k}.
$$

Coifman and Lafon introduced diffusion maps as a framework relating Markov
processes and geometric descriptions of data
([Coifman and Lafon 2006](https://doi.org/10.1016/j.acha.2006.04.006)).

Metric type: random-walk/operator distance, often represented by diffusion
coordinates. It integrates over many paths; it is not a graph shortest path.

When to compare: diffusion distance is a useful robust comparator for noisy
tubes and branches. It may be better for basin structure and worse for exact
arc-length recovery.

### Commute and Resistance Distance

For a connected weighted graph with Laplacian `L`, effective resistance is

$$
R_{ij}
= (e_i-e_j)^TL^+(e_i-e_j).
$$

Commute time is proportional to resistance:

$$
C_{ij}=\operatorname{vol}(G)R_{ij}.
$$

Metric type: electrical/random-walk distance. It depends on all parallel paths,
not the shortest path. In large random neighborhood graphs, commute times can be
dominated by degree and volume effects rather than global geometry
([von Luxburg, Radl, and Hein 2014](https://jmlr.csail.mit.edu/beta/papers/v15/vonluxburg14a.html)).

gflow use: comparator and graph redundancy diagnostic, not a geodesic metric.

### Heat-Kernel Distance

The manifold heat kernel satisfies Varadhan's small-time relation:

$$
d_M(x,y)^2
= -\lim_{t\to 0}4t\log H_t(x,y).
$$

On graphs, heat kernels use

$$
H_t = \exp(-tL).
$$

Metric type: heat/operator distance or descriptor; only the small-time limit
has a direct geodesic interpretation.

gflow use: multiscale stability diagnostic, especially for local neighborhoods
and graph smoothness.

### PHATE Potential Distance

PHATE builds an adaptive affinity graph, diffuses transition probabilities,
transforms them into potentials, and embeds potential distances by MDS
([Moon et al. 2019](https://www.nature.com/articles/s41587-019-0336-3),
[PHATE docs](https://phate.readthedocs.io/en/stable/)).

Sketch:

$$
P_t=P^t,
\qquad
D_{\mathrm{PHATE}}(i,j)
= \|U(P_t(i,\cdot))-U(P_t(j,\cdot))\|_2.
$$

Metric type: diffusion/potential distance and embedding distance. It is
support-aware, but not a shortest-path geodesic.

gflow use: important comparator for biological transitions and microbiome
branches, especially because it has practical adoption.

### UMAP Graph Geometry

UMAP constructs a weighted kNN fuzzy simplicial set. It estimates local
connectivity using `rho_i` and local scales `sigma_i`, then combines asymmetric
local fuzzy graphs by fuzzy set union. The UMAP docs describe the construction
as locally approximating geodesic distance at each point and combining local
fuzzy simplicial sets into a global graph
([UMAP API docs](https://umap-learn.readthedocs.io/en/latest/api.html?highlight=precomputed),
[UMAP overview](https://umap-learn.readthedocs.io/en/latest/index.html)).

Metric type: fuzzy topological graph plus embedding objective. The final UMAP
coordinates are embedding distances, not graph shortest-path distances. The
high-dimensional fuzzy graph can be used as a weighted graph, but then one must
define a path cost from membership strengths.

gflow use: comparator graph topology and local-scale diagnostic; avoid calling
UMAP embedding distances geodesic distances.

### Trajectory and Pseudotime Graph Distances

Pseudotime often uses distance from a root along an inferred graph, principal
curve, principal tree, cluster MST, or diffusion graph. This can be a true graph
geodesic on the inferred trajectory:

$$
\mathrm{ptime}(i)
= d_{G_{\mathrm{traj}}}(r, \operatorname{proj}(i)).
$$

But it is root-specific, often one-dimensional per lineage, and depends on the
trajectory abstraction.

gflow use: strong downstream comparator for branches and gradients, but not a
replacement for all-pairs intrinsic distance unless the intended geometry is a
tree/trajectory.

## Relevance to gflow Benchmarks

### Noisy Circle and Thick Tube

Best stressors:

- does the graph recover a loop rather than a path?
- are radial/chord shortcuts admitted?
- does mKNN fragment at sparse angular sectors?
- does MST/CMST force a long closure bridge or add false chords?
- does pruning remove redundant chords while preserving the fundamental cycle?

Recommended candidates:

- symmetrized kNN over `k` sweep;
- mKNN over `k` sweep;
- mKNN with component repair;
- iKNN/symmetrized kNN with local alternative-path pruning;
- CMST/MST completion as connected scaffold;
- UMAP/PHATE/diffusion as comparators, not targets.

### Gaps and Holes

The graph cannot decide whether a missing region is a sampling gap or a true
hole. `gflow` should therefore report two truth interpretations when possible:
continuous latent manifold with missing samples, and disconnected/forbidden
support with the gap treated as real.

Diagnostics: long bridge edges, number of completion edges crossing the gap,
path dependence on bridge edges, and downstream regression sensitivity when
bridges are removed or inflated.

### Branches and Figure-Eights

Branch points and self-approaches punish ambient-nearest-neighbor assumptions.
Good diagnostics are edge branch labels, endpoint-to-endpoint distortion,
false cross-branch edge rate, and whether pruning destroys or preserves true
branch-point connectivity.

### High-Dimensional Nuisance Variation

Neighbor search should be evaluated in raw features, PCA space, known latent
space for simulations, and compositional transforms. Approximate nearest-neighbor
recall is also a topology variable: a few wrong neighbors can become global
shortcuts.

### Compositional Data

For microbiome and metagenomics, Euclidean distances on raw relative abundance
are rarely the right local metric. Candidate spaces should include CLR with
pseudocount sensitivity, robust CLR/DEICODE-like representations, Aitchison
distance, filtering variants, and possibly UniFrac-derived neighbor graphs.

Downstream gflow regression should track whether distance improvements in
latent benchmarks translate into smoother conditional expectation estimates,
more stable gradient-flow calls, and interpretable edge/path attribution.

## Synthesis Table

| Family | Primary role | Metric type | Should gflow call it geodesic? | Main parameters | Main caution |
|---|---|---|---|---|---|
| Isomap kNN/epsilon graph | Explicit graph-geodesic estimator | Weighted graph metric; continuous geodesic estimator under assumptions | Yes, graph geodesic | `k`, `epsilon`, preprocessing, edge weights | Shortcuts vs disconnectedness |
| Directed kNN | Local directed topology | Directed graph distance | Only if asymmetric reachability intended | `k`, direction rule | Not a symmetric metric |
| Symmetrized kNN / iKNN | Local topology before shortest paths | Graph construction | Yes after weights and shortest paths | `k`, union rule | Hubs and dense-region shortcuts |
| Mutual kNN | Conservative topology | Graph construction | Yes after repair/weights if connected | `k` | Fragmentation |
| Shared-neighbor graph | Robust clustering topology | Similarity graph | Only after explicit path-cost definition | `k`, overlap, pruning threshold | Cluster similarity not arclength |
| Radius graph | Classical random geometric graph | Graph construction | Yes after weights | `epsilon` | Poor under nonuniform density |
| Hybrid kNN/radius | Topology control | Graph construction | Yes after weights | `k`, radius cap/local scale | Local scale can authorize outliers |
| Local scaling | Multi-scale affinity | Kernel/spectral graph | No, unless converted to path costs | scale neighbor, kernel form | May normalize rather than estimate geometry |
| Adaptive neighborhoods / IAN | Data-adaptive topology | Graph construction | Yes after weights | initial graph, optimization/pruning rules | Complexity and tuning |
| MST/CMST | Connectivity scaffold and diagnostic | Connected graph metric | Yes, graph geodesic on scaffold | MST metric, completion rule | Tree artifacts and forced bridges |
| Component repair | Connectivity repair | Graph modification | Yes after weights | bridge rule, local plausibility | May erase true disconnection |
| Alternative-path pruning | Shortcut/redundancy control | Graph modification | Yes after weights | stretch `tau` | Can ratify existing shortcuts |
| RNG/Gabriel/beta skeleton | Empty-region geometric pruning | Proximity graph | Yes after weights, cautiously | beta/empty-region rule | Ambient planar assumptions |
| Delaunay graph | Proximity supergraph | Proximity graph | Yes after weights, cautiously | triangulation metric | High-dimensional infeasibility/density |
| Geometric spanner | Sparse stretch control | Graph metric preserving chosen base metric | Only if base metric is appropriate | stretch `t`, cone count | Euclidean spanner can preserve shortcuts |
| Witness/landmark graph | Scalability/topology surrogate | Sparse graph/complex | Yes after weights, approximate | landmark selection, witness rule | Misses rare states |
| PAGA/topology abstraction | Topology summary | Cluster/partition graph | Usually no for sample distances | clustering, connectivity model | Abstract graph, not all-pairs metric |
| Diffusion distance | Random-walk comparator | Operator/embedding distance | No | kernel scale, alpha, diffusion time | Many-path basin effect |
| Commute/resistance | Electrical/random-walk comparator | Resistance metric / commute time | No | graph weights, Laplacian | Degree/volume domination |
| Heat-kernel distance | Multiscale operator geometry | Heat/operator distance | Only small-time limit | time `t`, Laplacian | Scale changes interpretation |
| PHATE potential distance | Biological transition comparator | Diffusion/potential/embedding distance | No | `knn`, decay, `t`, `gamma` | Not shortest path |
| UMAP fuzzy graph/embedding | Fuzzy topology and visualization | Fuzzy graph; embedding distance | No for embedding; graph only after path-cost definition | `n_neighbors`, `min_dist`, metric | Coordinates are not geodesic distances |
| Pseudotime graph distance | Rooted trajectory distance | Graph/curve distance on inferred trajectory | Yes on trajectory graph only | root, trajectory graph, projection | Root and topology dependent |

## What Is Known, Borrowed, and Distinctive for gflow

Known and borrowed:

- kNN/epsilon shortest paths are the classical finite-data geodesic estimator.
- asymptotic consistency needs shrinking local neighborhoods, sufficient
  sampling, and controlled smooth geometry.
- proximity graphs and spanners provide pruning language but not automatic
  manifold truth.
- single-cell analysis already accepts graph distances, MSTs, principal graphs,
  and diffusion comparators.

Unsettled:

- which graph topology best serves finite, noisy, uneven biological data;
- when forced connectedness is a feature versus a false bridge;
- how to choose among kNN, mKNN, adaptive graphs, pruning, and repair without
  overfitting synthetic cases;
- how to define truth when a gap may be biological absence or sampling absence.

Potentially distinctive for `gflow`:

- treating graph topology, edge weighting, component repair, and path
  attribution as separate inspectable modules;
- evaluating distances by downstream regression and gradient-flow stability,
  not only embedding appearance;
- compositional-data-aware neighbor construction for microbiome/metagenomics;
- explicit bridge and shortcut diagnostics on noisy circles, tubes, branches,
  gaps, and high-dimensional nuisance benchmarks;
- allowing ordinary kNN/mKNN/iKNN to beat CMST when the evidence says they do.

## Recommendations for the Benchmark Program

1. Keep the first benchmark axis topology-only: compare kNN, mKNN, radius,
   hybrid, adaptive, pruned, component-repaired, MST, and CMST graphs with the
   same Euclidean or compositional edge weights.
2. Separate topology diagnostics from edge-weight diagnostics.
3. Require parameter sweeps, not single-method exemplars.
4. For connectedness repair, report bridge edge identity and shortest-path usage.
5. For pruning, report stretch and false-shortcut effects, not only edge count.
6. Include diffusion, PHATE, UMAP, resistance, and pseudotime as comparators
   under their correct metric names.
7. In applied microbiome work, benchmark against Aitchison/UniFrac/Bray-Curtis
   baselines before claiming biological usefulness.
8. In single-cell-style branches, compare all-pairs graph distances to
   root-based pseudotime and principal-graph distances.
9. For every downstream gflow regression, retain path attribution so distance
   failures can be traced to specific graph edges.

## Source Index

- Classical graph geodesics and Isomap: [Tenenbaum et al. 2000 Science](https://doi.org/10.1126/science.290.5500.2319), [Bernstein et al. graph approximation report](https://www.researchgate.net/publication/2466340_Graph_Approximations_to_Geodesics_on_Embedded_Manifolds), [scikit-learn Isomap docs](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html), [scikit-learn manifold guide](https://scikit-learn.org/dev/modules/manifold.html).
- Landmark and out-of-sample Isomap: [de Silva and Tenenbaum Landmark MDS](https://graphics.stanford.edu/courses/cs468-05-winter/Papers/Landmarks/Silva_landmarks5.pdf), [Bengio et al. NeurIPS 2003](https://papers.nips.cc/paper_files/paper/2003/file/cf05968255451bdefe3c5bc64d550517-Paper.pdf).
- Local scaling and adaptive neighborhoods: [Self-tuning spectral clustering](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering), [IAN Neural Computation](https://direct.mit.edu/neco/article/35/3/453/114733/IAN-Iterated-Adaptive-Neighborhoods-for-Manifold), [IAN arXiv](https://arxiv.org/abs/2208.09123).
- Random geometric graph and kNN graph theory: [Penrose Random Geometric Graphs](https://academic.oup.com/book/9064), [randomized near-neighbor graph review](https://pmc.ncbi.nlm.nih.gov/articles/PMC7480951/), [Alamgir and von Luxburg ICML 2012 PDF](https://www.tml.cs.uni-tuebingen.de/team/luxburg/publications/AlamgirLuxburg2012.pdf), [Arias-Castro and Le Gouic 2017](https://arxiv.org/abs/1706.09441).
- Proximity graphs and pruning: [Relative neighborhood graph overview](https://en.wikipedia.org/wiki/Relative_neighborhood_graph), [Toussaint 1980 record](https://cir.nii.ac.jp/crid/1361981468897626368), [Gabriel graph overview](https://en.wikipedia.org/wiki/Gabriel_graph), [beta skeleton overview](https://en.wikipedia.org/wiki/Beta_skeleton), [Delaunay/RNG/beta skeleton construction](https://www.sciencedirect.com/science/article/pii/0925772194900183), [geometric spanner overview](https://en.wikipedia.org/wiki/Geometric_spanner).
- Single-cell graph adoption: [Monocle 2 / DDRTree](https://pmc.ncbi.nlm.nih.gov/articles/PMC5764547/), [scTEP trajectory review/background](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05179-2), [Slingshot methods table](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0/tables/1), [PAGA Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x), [VIA background](https://pmc.ncbi.nlm.nih.gov/articles/PMC8452770/), [DensityPath](https://academic.oup.com/bioinformatics/article/35/15/2593/5233001).
- Microbiome and bioinformatics adoption: [microarray Isomap example](https://pmc.ncbi.nlm.nih.gov/articles/PMC3940899/), [metagenomic 2D representation example](https://www.sciencedirect.com/science/article/pii/S2666389922002987), [Scanpy neighbors docs](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html), [Seurat FindNeighbors docs](https://satijalab.org/seurat/reference/findneighbors).
- Spatial/shape examples: [spatial transcriptomics point-cloud graph geodesics](https://pmc.ncbi.nlm.nih.gov/articles/PMC11227573/), [discrete geodesic graph on meshes](https://www.sciencedirect.com/science/article/abs/pii/S0167839617300390).
- Distinct comparator distances: [Coifman and Lafon diffusion maps](https://doi.org/10.1016/j.acha.2006.04.006), [von Luxburg et al. commute/hitting times](https://jmlr.csail.mit.edu/beta/papers/v15/vonluxburg14a.html), [PHATE Nature Biotechnology](https://www.nature.com/articles/s41587-019-0336-3), [PHATE docs](https://phate.readthedocs.io/en/stable/), [UMAP docs](https://umap-learn.readthedocs.io/en/latest/index.html), [UMAP API fuzzy simplicial set](https://umap-learn.readthedocs.io/en/latest/api.html?highlight=precomputed).

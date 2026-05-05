# Density-Aware Geodesic Distance Literature Map

Date: 2026-05-04

This memo is Deliverable 1 for the density-aware geodesic-distance research
track. It maps density-aware and density-adjacent distance constructions, then
reviews practical adoption in biology-facing data analysis. The main conclusion
is conservative: density-aware geodesic distances are not new, but explicit
Fermat/PWSPD-style geodesic distances appear underused in microbiome and
metagenomic workflows compared with diffusion, PHATE, UMAP, Aitchison/UniFrac,
and other standard beta-diversity approaches.

## Executive Map

The closest mathematical match to the current `gflow` proposal

$$
w_{ij} = \frac{\ell_{ij}}{\rho_{ij}^{\alpha}}.
$$

is the continuum Fermat / power-weighted shortest-path metric

$$
d_{\rho,\alpha}(x,y)
= \inf_{\gamma:x\to y}
\int_{\gamma} \rho(\gamma(s))^{-\alpha}\,ds.
$$

If the data are sampled from an intrinsic `d`-dimensional manifold with density
`rho`, and one uses the unrooted Fermat convention

$$
d_{F,p}(x,y)
= \inf_{\gamma:x\to y}
\int_{\gamma} \rho(\gamma(s))^{-(p-1)/d}\,ds,
$$

then

$$
\alpha = \frac{p-1}{d},
\qquad
p = 1 + \alpha d.
$$

Equivalently, this is the geodesic distance of the conformally changed
Riemannian metric tensor

$$
g_p = \rho^{-2(p-1)/d} g.
$$

Some papers write the tensor exponent as `2(p - 1) / d`; the path-length
integrand exponent, which corresponds to `ell / rho^alpha`, is half of that.

There is a convention trap. Some Fermat/PWSPD papers define the discrete
distance as

$$
\ell_p(x,y)
= \min_{\text{path}}
\left(\sum_i \|x_i - x_{i+1}\|^p\right)^{1/p},
$$

while persistent-homology papers often use the unrooted version

$$
d_{X,p}(x,y)
= \min_{\text{path}} \sum_i \|x_i - x_{i+1}\|^p.
$$

The continuum density-weighted geodesic that matches `gflow` edge weighting is
the unrooted path metric, or `ell_p^p` under the rooted convention. This should
be stated explicitly in any future design memo.

## Historical Spine

1. Classical geodesic manifold learning starts from local Euclidean graphs:
   Isomap estimates manifold geodesic distances by shortest paths on a kNN or
   epsilon graph before MDS ([Tenenbaum et al. 2000](https://doi.org/10.1126/science.290.5500.2319)).

2. Density-based distances entered semi-supervised learning as a way to make
   points close when they are connected through high-density unlabeled data.
   Sajama and Orlitsky's ICML 2005 work is an early explicit source for
   density-based Riemannian path metrics ([DBLP record](https://dblp.org/rec/conf/icml/SajamaO05)).
   Bijral, Ratliff, and Srebro later proposed graph-based density-based
   distances for semi-supervised learning ([arXiv:1202.3702](https://arxiv.org/abs/1202.3702)).

3. Fermat distances / power-weighted shortest paths gave a density-sensitive
   graph-distance estimator that avoids explicit density estimation by raising
   edge lengths to a power. The theory connects to nonhomogeneous Euclidean
   first-passage percolation and Fermat's principle in optics
   ([Groisman et al.](https://arxiv.org/abs/1810.09398)).

4. Hwang, Damelin, and Hero analyzed shortest paths through random points and
   proved convergence of power-weighted sample paths to a density-weighted
   continuum distance ([arXiv:1202.0045](https://arxiv.org/abs/1202.0045)).
   This is a key bridge from density-based distances to the modern Fermat
   literature.

5. Recent theory made the continuum limit and spectral behavior precise:
   Garcia Trillos, Little, McKenzie, and Murphy prove convergence and spectral
   results for Fermat distances ([JMLR 2024](https://jmlr.org/papers/v25/23-0939.html)).
   Fernandez, Borghini, Mindlin, and Groisman use Fermat distance for
   intrinsic persistent homology ([JMLR 2023](https://jmlr.org/papers/v24/21-1044.html)).
   Little, McKenzie, and Murphy analyze PWSPDs and kNN graph spanners
   ([SIAM MDS 2022](https://epubs.siam.org/doi/10.1137/20M1386657),
   [arXiv](https://arxiv.org/abs/2012.09385)).

6. Diffusion geometry follows a parallel lineage. Diffusion maps define a
   Markov operator on a kernel graph and compare random-walk transition
   distributions; the alpha normalization controls whether sampling density is
   exploited or removed ([Coifman and Lafon 2006](https://doi.org/10.1016/j.acha.2006.04.006)).

7. PHATE uses diffusion probabilities and then transforms them into potentials
   before embedding. It is support-aware and density-adjacent, but it is not a
   shortest-path geodesic metric ([Moon et al. 2019](https://www.nature.com/articles/s41587-019-0336-3),
   [PHATE docs](https://phate.readthedocs.io/en/stable/)).

## Construct 1: Continuous Density-Weighted Path Metrics

Objects:

$$
\begin{aligned}
(M,g) &: \text{smooth connected Riemannian manifold},\\
\rho &: M \to \mathbb{R}_{+} \quad \text{sampling/support density},\\
\alpha &\ge 0 \quad \text{density penalty exponent},\\
\gamma &: \text{piecewise } C^1 \text{ curve from } x \text{ to } y,\\
ds_g &: \text{Riemannian arclength}.
\end{aligned}
$$

Definition:

$$
L_{\rho,\alpha}(\gamma)
= \int_0^1
\rho(\gamma(t))^{-\alpha}
\sqrt{g_{\gamma(t)}(\gamma'(t),\gamma'(t))}\,dt,
$$

$$
d_{\rho,\alpha}(x,y)
= \inf_{\gamma:x\to y} L_{\rho,\alpha}(\gamma).
$$

If `rho` is positive and finite on each connected component, this is a true
path metric on that component. It is the Riemannian geodesic distance for the
conformal tensor

$$
g_{\alpha} = \rho^{-2\alpha} g.
$$

If `rho` reaches zero and no floor is imposed, distances crossing zero-density
regions can become infinite. If `rho` is estimated and clipped by a floor, the
result is a metric for the clipped density but may encode the floor more than
the biology in tails.

Density role: high-density regions reduce path cost. Low-density regions act as
barriers. This is an exploitation of density as support or plausibility, not a
normalization that removes density.

Parameter correspondence:

$$
\begin{aligned}
\text{Fermat/PWSPD path integrand} &: \rho^{-(p-1)/d}\,ds,\\
\text{\texttt{gflow} }\alpha &: \alpha = \frac{p-1}{d},\\
\text{Riemannian tensor exponent} &: 2\alpha = \frac{2(p-1)}{d}.
\end{aligned}
$$

Finite-data construction:

1. Estimate `rho_i` or `rho_ij` in an appropriate space.
2. Build a graph, usually kNN, mutual kNN, iKNN/mKNN, radius, MST/CMST, or a
   repaired local graph.
3. Define edge weights such as

   $$
   w_{ij}
   = \frac{\ell_{ij}}{\phi(\rho_i,\rho_j,\rho_{ij})^{\alpha}}.
   $$

4. Compute Dijkstra or all-pairs shortest paths.

Strengths: clear mathematical object; direct control of density sensitivity;
compatible with edge/path attribution; naturally matches noisy-tube oracle
truths where support is part of the geometry.

Failure modes: density-estimation instability; dimensional scaling errors;
tail domination for `alpha` near 1 in low-density regions; conflating true
biological rarity with sampling failure; bridges can be made too long or too
short depending on density summary.

Diagnostic requirements: density distribution, floor/cap fraction, edge
`ell_ij` versus `w_ij`, path examples, alpha sensitivity, estimator sensitivity,
tail/outlier influence, and comparison to Euclidean, latent, Fermat, and
diffusion distances when available.

## Construct 2: Fermat Distances and PWSPD

Discrete rooted PWSPD:

$$
\ell_p(x,y;X)
= \min_{x=x_0,\ldots,x_m=y}
\left(\sum_{r=0}^{m-1}\|x_{r+1}-x_r\|^p\right)^{1/p},
\qquad p \ge 1.
$$

Unrooted Fermat distance:

$$
d_{X,p}(x,y)
= \min_{\text{path}}
\sum_{r=0}^{m-1}\|x_{r+1}-x_r\|^p.
$$

For `p = 1`, the complete-graph rooted distance collapses to Euclidean
distance; on a restricted local graph, it is the usual Isomap-style graph
geodesic. For `p > 1`, many short hops are favored, so paths tend to stay in
dense regions. As `p -> infinity`, the optimization approaches a bottleneck /
longest-edge path distance.

Continuum limit:

$$
d_{F,p}(x,y)
= \inf_{\gamma:x\to y}
\int_{\gamma} \rho(\gamma(s))^{-(p-1)/d}\,ds.
$$

With appropriate normalization, sample Fermat distances converge to the
population density-weighted geodesic on an intrinsic `d`-dimensional manifold.
The normalization often involves `n^((p - 1) / d)` for the unrooted convention
or `n^((p - 1) / (p d))` for the rooted convention, up to a percolation
constant depending on `p` and `d`.

Hwang, Damelin, and Hero's earlier notation makes the same correspondence
explicit:

$$
n^{(p-1)/d} d_{X_n,p}(x,y)
\to C(d,p)d_{F,p}(x,y),
$$

$$
g_p = \rho^{2(1-p)/d} g.
$$

The sign in the tensor exponent is the same as `rho^(-2(p - 1) / d)`.

Asymptotic and finite-sample behavior:

- JMLR 2024 establishes local convergence of discrete Fermat distances to
  continuum Fermat distances and spectral convergence of Fermat graph
  Laplacians ([Garcia Trillos et al.](https://jmlr.org/papers/v25/23-0939.html)).
- SIAM MDS 2022 proves that PWSPD on a kNN graph can match the complete-graph
  PWSPD with high probability under k of order `log n`, with constants depending
  on intrinsic dimension, density bounds, `p`, and manifold geometry
  ([Little et al.](https://arxiv.org/abs/2012.09385)).
- Parameter choice is not benign: Chazal et al. emphasize a tradeoff between
  large enough `alpha`/`p` to follow geometry and small enough values to avoid
  noise-amplified distance estimation ([arXiv:2311.18663](https://arxiv.org/abs/2311.18663)).

Practical strengths:

- No explicit density estimate is required.
- kNN graph computation is feasible for large `n`.
- Directly implements the geometry-density interpolation wanted by `gflow`.

Practical failures:

- Scaling across `n`, `p`, and `d` is easy to mishandle.
- Large `p` can become bottleneck-like and unstable.
- Complete-graph intuition and kNN-graph behavior differ if k is too small or
  if local neighborhoods are distorted by high-dimensional noise.

Relevance to `gflow`: this family should be a named comparator and possibly a
special case of density-aware graph distance. The noisy-tube oracle can be
parameterized to reproduce its continuum limit by setting
`alpha = (p - 1) / d`.

## Construct 3: Density-Based Persistent Homology

Motivation: persistent homology of a point cloud in ambient Euclidean distance
is sensitive to embedding distortion and outliers. Fermat distances provide an
intrinsic metric that accounts for both manifold geometry and sampling density.

Definition: compute Vietoris-Rips persistence on `(X_n, d_{X_n,p})`, where
`d_{X_n,p}` is the sample Fermat distance. In the continuum limit, this
approximates persistent homology of `(M, d_{F,p})`.

Metric type: metric-space input for persistent homology; the Fermat component
is a density-weighted geodesic sampling-limit object.

Key source: Fernandez et al. prove convergence of persistence diagrams and
robustness advantages for positive-degree homology, then apply the method to
time-series pattern recognition and anomaly detection, including ECG signals
and canary song ([JMLR 2023](https://jmlr.org/papers/v24/21-1044.html)).

Adoption status: strong as a theory-plus-application paper for time series;
not yet a common microbiome or metagenomics tool.

Relevance to `gflow`: provides a topology-facing benchmark for density-aware
distances. For noisy circles, annuli, figure-eights, and gaps, persistence under
Euclidean, latent, Fermat, and tube-oracle distances could diagnose whether
density weighting preserves or destroys intended topology.

## Construct 4: Graph Geodesics, Neighbor Selection, and Local Scaling

### Isomap-style graph geodesics

Definition:

$$
d_G(i,j)
= \min_{\text{paths } i\to j}
\sum_{(u,v)\in\text{path}} \|x_u-x_v\|.
$$

This is a graph metric on each connected component. Density enters implicitly
through graph topology: dense regions have many short edges, sparse regions
have fewer or longer edges. It is not explicitly density weighted.

Important caution: unweighted shortest path distances on random kNN graphs can
prefer low-density regions because a fixed number of hops travels farther where
points are sparse. This is a useful warning for any `gflow` graph that replaces
length weights with hop counts or over-normalizes edge lengths
([Alamgir and von Luxburg 2012](https://icml.cc/2012/papers/545.pdf)).

### kNN, mKNN, iKNN, and adaptive k

kNN graphs are computationally attractive but may over-connect dense regions
and under-connect sparse ones. Mutual kNN is more conservative and can fragment.
iKNN can repair connectivity but can become too dense. Adaptive k and local
repair should be evaluated as graph-topology choices before applying
density-weighted edge lengths.

Metric type: graph distance after a neighbor-selection rule.

Density role: implicit through neighbor radii and degree; explicit only if k or
edge admission depends on density estimates.

Relevance to `gflow`: density-aware weighting should not be the only lever.
Graph topology decides whether unsupported shortcuts exist at all. Weighting a
bad shortcut can help only if the density summary recognizes the unsupported
gap.

### Local scaling and density-aware kernels

Self-tuning spectral clustering uses a local scale, often the distance to the
`k`th neighbor:

$$
A_{ij}
= \exp\left(
-\frac{\|x_i-x_j\|^2}{\sigma_i\sigma_j}
\right).
$$

This adapts to multi-scale data without explicit density estimation
([Zelnik-Manor and Perona 2004](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering)).
In biology, Spectrum develops an adaptive density-aware kernel for clustering
single and multi-omic data, using common-neighbor information to strengthen
within-density connections and providing R software
([Spectrum, Bioinformatics](https://academic.oup.com/bioinformatics/article/36/4/1159/5566508)).

Metric type: kernel/spectral geometry, not a shortest-path geodesic unless
converted into graph distances.

Relevance to `gflow`: local scales are natural density surrogates for edge
admission and for robust floors/caps.

### Mutual reachability distances

HDBSCAN transforms pairwise distances by local core distances:

$$
d_{\mathrm{mreach}}(i,j)
= \max\{\mathrm{core}_k(i),\mathrm{core}_k(j),d(i,j)\}.
$$

This spreads sparse points apart and stabilizes density-based clustering
([HDBSCAN docs](https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html)).

Metric type: graph/clustering dissimilarity; commonly used through an MST and
cluster hierarchy.

Relevance to `gflow`: mutual reachability is a robust alternative to
`ell / rho^alpha` for penalizing low local support. It is especially relevant
to bridge admission and outlier handling.

## Construct 5: Diffusion Maps and Diffusion Distance

Objects:

$$
\begin{aligned}
K_{\varepsilon}(i,j)
&= \exp\left(-\frac{\|x_i-x_j\|^2}{\varepsilon}\right),\\
q_i &= \sum_j K_{\varepsilon}(i,j),\\
K_a(i,j) &= \frac{K_{\varepsilon}(i,j)}{q_i^a q_j^a},\\
d_i &= \sum_j K_a(i,j),\\
P_a(i,j) &= \frac{K_a(i,j)}{d_i}.
\end{aligned}
$$

Here `a` is the diffusion-map alpha normalization, not the `gflow` density
exponent. The diffusion distance at time `t` can be written as a weighted
distance between transition distributions:

$$
D_t^2(i,j)
= \sum_k
\frac{\left(P_a^t(i,k)-P_a^t(j,k)\right)^2}{\pi_k}.
$$

where `pi` is the stationary distribution or an associated reference measure.
Spectrally, this is represented by diffusion coordinates with eigenvalues
raised to time `t`.

Density role:

- `a = 0`: no kernel density renormalization; the limiting operator includes a
  density drift and random walks are strongly density influenced.
- `a = 1`: approximates the Laplace-Beltrami operator and removes sampling
  density from the continuum geometry.
- intermediate values interpolate between density-sensitive and
  density-corrected diffusion. Many implementations default near `a = 0.5`.

pydiffmap describes `alpha` as the exponent used in left normalization and
supports density functions and variable bandwidths
([pydiffmap docs](https://pydiffmap.readthedocs.io/en/master/reference/diffusion_map.html)).
datafold explicitly describes `alpha = 1` as correcting sampling density
([datafold docs](https://datafold-dev.gitlab.io/datafold/api/datafold.dynfold.DiffusionMaps.html)).

Metric type: random-walk / embedding distance, not a shortest-path geodesic.
It is a metric or pseudometric depending on retained eigenvectors, time, graph
connectivity, and whether identical transition profiles occur.

Practical strengths:

- Integrates many paths rather than choosing one.
- Robust to small local graph noise.
- Very mature in single-cell trajectory work.

Failure modes:

- Diffusion time changes scale and topology perception.
- Density normalization may remove biologically meaningful population mass.
- A random-walk distance can make two points close because they diffuse to the
  same basin, even if the shortest supported path between them is long.

Relevance to `gflow`: diffusion distances are essential comparators and
possibly useful diagnostics, but should not be labeled geodesic distances
unless a heat-kernel small-time/geodesic approximation is explicitly intended.

## Construct 6: PHATE Potential Distance

PHATE builds an adaptive affinity graph, constructs a Markov diffusion
operator, chooses a diffusion time `t`, transforms transition probabilities
into potentials, and embeds by MDS.

Current software defaults include:

```text
knn = 5
decay = 40
n_landmark = 2000
t = "auto"
gamma = 1
n_pca = 100
knn_dist = "euclidean"
```

`t` is the power of the diffusion operator and controls diffusion scale.
`gamma = 1` gives log potential; `gamma = 0` gives square-root potential
([PHATE docs](https://phate.readthedocs.io/en/stable/)).

Formal sketch:

$$
P_t = P^t.
$$

$$
U_{\gamma}(P_t(i,\cdot)) =
\begin{cases}
\log P_t(i,\cdot), & \gamma = 1,\\
\text{power or square-root-like potential}, & \gamma \approx 0,
\end{cases}
$$

$$
D_{\mathrm{PHATE}}(i,j) =
\left\|U_{\gamma}(P_t(i,\cdot)) -
U_{\gamma}(P_t(j,\cdot))\right\|_2.
$$

Metric type: potential/embedding distance derived from diffusion probabilities.
It is support-aware and density-adjacent, but not a shortest-path geodesic.

Density role: mediated by kNN graph construction, kernel tails, Markov
normalization, diffusion time, and potential transform. High-density paths
increase transition probabilities, but row normalization and diffusion can also
smooth or erase density contrast.

Practical adoption: strong in single-cell, mass cytometry, Hi-C, and some
microbiome work. The PHATE paper reports gut microbiome data among applications
([Moon et al. 2019](https://www.nature.com/articles/s41587-019-0336-3)).

Relevance to `gflow`: PHATE should be exposed as a diffusion/potential
comparator, not as a density-weighted geodesic. A future `gflow` PHATE-like
option should name the output as potential or diffusion distance.

## Construct 7: Random-Walk, Commute, Resistance, Heat-Kernel, and PPR Distances

### Commute time and resistance distance

For a connected weighted graph with Laplacian `L`, effective resistance is

$$
R_{ij}
= (e_i-e_j)^{T}L^{+}(e_i-e_j).
$$

where `L^+` is the Moore-Penrose pseudoinverse. Commute time is proportional
to resistance:

$$
C_{ij} = \operatorname{vol}(G) R_{ij}.
$$

Metric type: graph metric for resistance; random-walk expected-time distance
for commute.

Density role: high-conductance, many-path regions have low resistance. This is
not a shortest path; parallel paths matter.

Failure mode: in large random neighborhood graphs, commute and hitting times
can become dominated by vertex degrees rather than meaningful global geometry
([von Luxburg et al. 2014](https://jmlr.csail.mit.edu/beta/papers/v15/vonluxburg14a.html)).

Biology adoption: landscape ecology uses resistance and circuit theory to model
movement and gene flow through environmental resistance surfaces; this is
explicit environmental resistance, not sample-density weighting. R packages
such as `gdistance` implement commute/resistance distances
([gdistance docs](https://search.r-project.org/CRAN/refmans/gdistance/html/commuteDistance-methods.html)).

### Heat-kernel distances

On a manifold, the heat kernel `H_t(x, y)` is linked to geodesic distance by
Varadhan's small-time formula:

$$
d_g(x,y)^2
= -\lim_{t\to 0} 4t \log H_t(x,y).
$$

On graphs, one often uses

$$
H_t = \exp(-tL).
$$

Distances can be formed from rows of `H_t`, from heat-kernel signatures
`H_t(x, x)`, or from small-time geodesic approximations.

Metric type: diffusion/heat distance or shape descriptor; only in a small-time
limit does it recover geodesic distance.

Adoption: heat-kernel signatures are well established for shape analysis and
medical/biological morphology, but usually as descriptors for matching,
retrieval, segmentation, or classification rather than density-aware sample
geodesics ([Sun et al. HKS page](https://geometry.stanford.edu/paper.php?id=sog-hks-09)).

### Personalized PageRank and random-walk-with-restart distances

PPR solves

$$
r_i = \beta e_i + (1-\beta)r_iP.
$$

Distances or similarities are built from personalized stationary vectors.
They emphasize local graph basins and are useful for ranking, clustering, and
semi-supervised propagation.

Heat-kernel PageRank replaces the geometric mixture of walk lengths with a
Poisson mixture:

$$
h_{t,i} = e_i\exp(-t(I-P)).
$$

Metric type: usually asymmetric proximity; symmetric variants may be
pseudometrics after transformation.

Biology adoption: random-walk-with-restart and local random walks are common in
network biology and single-cell fate tools, but their density role is generally
implicit through graph construction.

Relevance to `gflow`: random-walk clouds from a vertex can be used as local
volume/density probes, but that should be treated as a diagnostic or density
estimator unless a distance is explicitly defined and tested.

## Construct 8: Density-Aware Oracle Truth for Synthetic Benchmarks

For a known generative density on a manifold, tube, grid, or graph, define

$$
d_{\mathrm{oracle}}(x,y)
= \inf_{\gamma:x\to y}
\int_{\gamma} \rho(\gamma(s))^{-\alpha}\,ds,
$$

or build a dense oracle graph with edge weights

$$
w_{ab}
= \frac{\|z_a-z_b\|}
{\rho(\operatorname{midpoint}(a,b))^{\alpha}}.
$$

Metric type: synthetic truth metric, not a data estimator.

Conditions:

- positive density floor or explicit support boundary;
- convergence checks over oracle grid resolution;
- transparent sample-to-oracle attachment;
- sensitivity over `alpha`, density floor, cap, and attachment rule.

Relevance to `gflow`: the noisy-tube oracle is a correct and useful internal
object if presented as an oracle benchmark for a chosen density-weighted truth,
not as a new distance theory.

## Practical Adoption Review

### Microbiome, 16S rRNA, Metagenomics, and Compositional Counts

Observed adoption:

- Standard microbiome workflows overwhelmingly use beta-diversity distances:
  Bray-Curtis, Jaccard, Jensen-Shannon, weighted/unweighted UniFrac,
  generalized UniFrac, Aitchison/robust Aitchison, and related ordinations.
  These are not sample-density geodesics.
- PHATE has concrete microbiome use. A large gut microbiome analysis used PHATE
  to map local ecological states into global branches and then performed
  pseudotime analysis along those branches
  ([Nature Communications 2023](https://www.nature.com/articles/s41467-023-38558-7)).
  The PHATE effect is implicit diffusion/potential geometry, not an explicit
  `rho^-alpha` geodesic.
- tmap uses topological data analysis and network representation for
  population-scale microbiome stratification
  ([Genome Biology 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1871-4)).
  It is important adoption context for graph/TDA microbiome workflows, but not
  a Fermat geodesic method.
- Diffusion maps appear in ecological/microbiome-adjacent trait-space work,
  including bacterial metabolic niche maps based on genome-scale metabolic
  traits combined with Earth Microbiome Project 16S community data
  ([Fahimipour and Gross 2020](https://www.nature.com/articles/s41467-020-18695-z)).
- Compositional methods such as DEICODE / robust Aitchison PCA address sparse
  microbiome counts directly, but are not graph geodesic density metrics
  ([Martino et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6372836/)).

Evidence gap:

- I found little evidence that Fermat distance, PWSPD, or explicit
  `ell / rho^alpha` graph geodesics are standard in microbiome, 16S,
  metagenomic, or metatranscriptomic analysis.
- This gap is a plausible `gflow` opportunity, but only if density estimation
  is compositional-data-aware and benchmarked against Aitchison/UniFrac/PHATE.

Scientific questions motivating use:

- gradients and branches of community states;
- stratification of host/environment associations;
- visualization and ordination;
- pseudotime-like ecological transitions.

Software/defaults/diagnostics:

- PHATE and phateR provide defaults and usable software.
- Microbiome-specific density-aware geodesic diagnostics appear largely absent.

### Single-Cell Genomics and Developmental Trajectories

Observed adoption:

- Diffusion maps and diffusion pseudotime are deeply adopted. `destiny`
  computes cell-to-cell transition probabilities, applies density normalization,
  and uses eigendecomposition for diffusion components
  ([destiny paper](https://academic.oup.com/bioinformatics/article/32/8/1241/1744143),
  [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/destiny/inst/doc/Diffusion-Maps.html)).
- DPT reconstructs lineage branching from diffusion geometry
  ([Haghverdi et al. 2016](https://www.nature.com/articles/nmeth.3971)).
- PHATE is heavily adopted for visualizing high-dimensional biological
  transitions, including scRNA-seq and mass cytometry
  ([Moon et al. 2019](https://www.nature.com/articles/s41587-019-0336-3)).
- MAGIC uses graph diffusion for denoising/imputation rather than reporting a
  distance matrix ([MAGIC paper](https://pubmed.ncbi.nlm.nih.gov/29961576/)).
- Palantir, Wishbone, MARGARET, and related tools use diffusion maps,
  Markov chains, waypoints, random walks, and fate probabilities rather than
  explicit density-weighted shortest paths.
- scPMP is the clearest biology adoption of PWSPD. It uses power-weighted path
  metrics for scRNA-seq clustering and visualization, computes shortest paths
  on kNN graphs with edge weights `||x_i - x_j||^p`, and provides software and
  defaults ([PLOS Computational Biology 2024](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012014)).
- DensityPath builds a density landscape after dimensional reduction, extracts
  high-density representative cell states, and connects them by geodesic MST on
  the density surface ([Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/15/2593/5233001)).

Density interpretation:

- Single-cell diffusion methods often correct away sampling density because
  rare and abundant populations should both be represented.
- scPMP and DensityPath are the strongest precedents for exploiting density as
  biological support.

Relevance to `gflow`:

- Single-cell gives both caution and precedent. It shows adoption is possible,
  but also that the field distinguishes density correction from density
  exploitation depending on biological question.

### Spatial Transcriptomics

Observed adoption:

- Spatial methods commonly use spatial neighbor graphs, graph convolution,
  diffusion propagation, or physical-distance kernels.
- I found less evidence for Fermat/PWSPD-style density-weighted geodesic
  distances over spatial transcriptomic cell-state manifolds.

Likely transferable lesson:

- Separate density in measurement space, spatial density of spots/cells, and
  biological density of states. A single `rho` may conflate them.

### Ecology and Community Composition

Observed adoption:

- Ecology uses ordination, beta-diversity, distance decay, environmental
  gradients, and resistance surfaces extensively.
- Diffusion-map dissimilarity has been proposed for ecological community
  composition when communities have high species turnover and few shared
  species ([Gault et al. 2023](https://epic.awi.de/id/eprint/58032/)).
- Circuit/resistance distances are common for movement, dispersal, and gene
  flow through landscapes, but their "density" is environmental conductance or
  resistance, not sample density.

Relevance to `gflow`:

- Ecology is conceptually close to microbiome state-space work. It reinforces
  that support-aware paths are meaningful, but also that the resistance surface
  must be scientifically interpretable.

### Medical Imaging, Morphology, and Biological Shape

Observed adoption:

- Heat-kernel and diffusion descriptors are used for shape matching,
  segmentation, and morphology because they are intrinsic and multiscale.
- These methods generally use geometry of an observed surface or mesh rather
  than sampling-density-weighted data geodesics.

Relevance to `gflow`:

- Heat-kernel diagnostics may be useful for graph stability and multiscale
  structure, but should be kept separate from shortest-path density-weighted
  distances.

### Non-Biological Applied Domains

Observed adoption:

- Fermat/PWSPD has applications in clustering, semi-supervised learning, image
  segmentation, and time-series pattern recognition.
- Density-based graph distances are also central to DBSCAN, OPTICS, HDBSCAN,
  and density-based cluster validation.

Transferable lessons:

- Explicit density-aware path metrics work best when the density is meaningful
  and stable.
- Large density exponents are powerful but fragile.
- kNN graph spanner results justify computational shortcuts, but constants and
  preprocessing matter in finite biological data.

## Relevance to gflow Design

Known/borrowed:

- `ell / rho^alpha` is a standard conformal path metric.
- Fermat/PWSPD supplies a density-free estimator of a related continuum metric.
- Diffusion/PHATE supplies support-aware random-walk comparators.

Potentially distinctive:

- compositional and microbiome-aware density estimation;
- diagnostics for floors, caps, tails, and bridge attribution;
- integration with CMST/iKNN/mKNN graph construction and repair;
- noisy-tube oracle benchmarks;
- downstream gradient-flow and conditional expectation analysis.

Immediate recommendations for later deliverables:

1. Parameterize any explicit density-weighted graph metric by `alpha`, and
   report the implied Fermat `p = 1 + alpha d_hat` when an intrinsic dimension
   estimate is available.
2. Implement Fermat/PWSPD as a comparator, not only explicit density weighting.
3. Keep diffusion/PHATE distances in a separate family named "diffusion" or
   "potential", not "geodesic".
4. Require density floors/caps and tail diagnostics before any density-aware
   distance is used downstream.
5. Benchmark against microbiome-native comparators: Aitchison/robust
   Aitchison, Bray-Curtis, UniFrac, PHATE, UMAP/densMAP, and graph-geodesic
   baselines.

## Synthesis Table

| Family | Density role | Metric type | Should gflow call it geodesic? | Adoption signal | Main caution |
|---|---|---|---|---|---|
| Continuous `integral rho^-alpha ds` | Explicitly exploits density as support | Path/geodesic metric if `rho > 0` | Yes | Theory, semi-supervised learning | Density estimation and tails |
| Explicit edge reweighting `ell / phi(rho)^alpha` | Explicitly exploits estimated density | Graph shortest-path metric | Yes, graph geodesic | Limited biology evidence | Bad graph topology can still shortcut |
| Fermat/PWSPD | Implicitly exploits density through edge powers | Sample graph metric with continuum geodesic limit | Yes, with convention stated | Strong theory; scPMP in scRNA-seq | Scaling with `n,p,d`; large `p` noise |
| Density-based persistent homology | Uses Fermat metric for topology | Metric-space/TDA object | The distance yes; PH output no | Time series applications | Less used in microbiome |
| Isomap/kNN geodesic | Density affects topology implicitly | Graph shortest-path metric | Yes, graph geodesic | Broad | kNN shortcuts/fragments; hop-count pathology |
| mKNN/iKNN/adaptive k | Density affects topology/degree | Graph construction | Only after weights define paths | Local graph practice | Can fragment or overconnect |
| Local scaling/self-tuning kernels | Usually corrects/adapts to local scale | Kernel/spectral geometry | No, unless converted to paths | Spectrum in omics | Kernel effects can oppose density weighting |
| HDBSCAN mutual reachability | Penalizes low local density/core support | Clustering dissimilarity/MST | Not usually | Broad clustering | Clustering, not distance-estimation goal |
| Diffusion maps | Can exploit or correct density via alpha | Random-walk/embedding distance | No | Very strong single-cell | Alpha meaning differs from gflow alpha |
| PHATE potential distance | Implicit support through diffusion probabilities | Potential/embedding distance | No | Strong single-cell; some microbiome | Not shortest-path; depends on `t`, `gamma` |
| Commute/resistance | Many high-conductance paths reduce distance | Random-walk/electrical graph metric | No | Ecology/network biology | Degree domination in large graphs |
| Heat-kernel distance/signature | Diffusion over geometry | Heat/diffusion descriptor | Only in small-time limit | Shape/morphology | Time scale changes interpretation |
| Personalized PageRank/RWR | Local basin probability | Proximity, often asymmetric | No | Network biology | Not automatically metric |
| Density-aware oracle | Known generative density defines truth | Synthetic benchmark metric | Yes, if path metric | Internal benchmark | Oracle floor/attachment artifacts |

## Source Index

- Fermat/PWSPD theory: [Hwang et al. arXiv:1202.0045](https://arxiv.org/abs/1202.0045), [JMLR 2024](https://jmlr.org/papers/v25/23-0939.html), [SIAM MDS 2022](https://epubs.siam.org/doi/10.1137/20M1386657), [arXiv:2012.09385](https://arxiv.org/abs/2012.09385), [arXiv:1810.09398](https://arxiv.org/abs/1810.09398).
- Fermat persistent homology: [JMLR 2023](https://jmlr.org/papers/v24/21-1044.html).
- Fermat parameter selection: [arXiv:2311.18663](https://arxiv.org/abs/2311.18663).
- Semi-supervised density-based distances: [arXiv:1202.3702](https://arxiv.org/abs/1202.3702), [ICML 2005 DBLP](https://dblp.org/rec/conf/icml/SajamaO05).
- Diffusion maps: [Coifman and Lafon DOI](https://doi.org/10.1016/j.acha.2006.04.006), [pydiffmap docs](https://pydiffmap.readthedocs.io/en/master/reference/diffusion_map.html), [datafold docs](https://datafold-dev.gitlab.io/datafold/api/datafold.dynfold.DiffusionMaps.html).
- PHATE: [Nature Biotechnology 2019](https://www.nature.com/articles/s41587-019-0336-3), [PHATE docs](https://phate.readthedocs.io/en/stable/).
- Single-cell diffusion adoption: [destiny Bioinformatics](https://academic.oup.com/bioinformatics/article/32/8/1241/1744143), [DPT Nature Methods](https://www.nature.com/articles/nmeth.3971), [MAGIC Cell](https://pubmed.ncbi.nlm.nih.gov/29961576/).
- Single-cell density/path adoption: [scPMP PLOS Computational Biology 2024](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012014), [DensityPath Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/15/2593/5233001).
- Microbiome adoption: [Gut microbiome branches with PHATE](https://www.nature.com/articles/s41467-023-38558-7), [tmap Genome Biology 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1871-4), [DEICODE/robust Aitchison](https://pmc.ncbi.nlm.nih.gov/articles/PMC6372836/).
- Local scaling and density-aware kernels: [Self-tuning spectral clustering](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering), [Spectrum Bioinformatics](https://academic.oup.com/bioinformatics/article/36/4/1159/5566508).
- Random-walk and resistance cautions: [von Luxburg et al. JMLR](https://jmlr.csail.mit.edu/beta/papers/v15/vonluxburg14a.html), [gdistance commuteDistance](https://search.r-project.org/CRAN/refmans/gdistance/html/commuteDistance-methods.html).

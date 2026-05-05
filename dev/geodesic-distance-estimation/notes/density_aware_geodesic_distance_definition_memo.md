# Density-Aware Geodesic Distance Definition Memo

Date: 2026-05-04

This memo is Deliverable 2 for the density-aware geodesic-distance research
track. It fixes the working `gflow` vocabulary for density-aware distances. The
companion literature map covers history, sources, and adoption; this memo is
more axiomatic. Its goal is to say what objects `gflow` should mean when it
says density-aware geodesic, graph distance, Fermat distance, diffusion
distance, PHATE potential distance, or oracle truth.

## 1. Notation

Let

$$
(M,g)
$$

be a connected smooth Riemannian manifold, possibly with boundary, with
arclength element \(ds_g\). Let

$$
\rho: M \to [0,\infty)
$$

denote a sampling, support, or biological-state density. In finite data, let

$$
X = \{x_1,\ldots,x_n\}\subset \mathbb{R}^D
$$

or more generally \(X\subset M\). A graph on \(X\) is

$$
G=(V,E),\qquad V=\{1,\ldots,n\}.
$$

For an edge \(e=(i,j)\), let \(\ell_{ij}\ge 0\) be its base geometric length,
usually

$$
\ell_{ij} = \|x_i-x_j\|
$$

after the chosen preprocessing. Let \(\widehat{\rho}_i\) be an estimated vertex
density and \(\widehat{\rho}_{ij}\) an edge density summary.

Throughout this memo, \(\alpha\ge 0\) is the `gflow` density penalty exponent:
larger \(\alpha\) makes low-density regions more expensive. This is not the
same parameter as diffusion-map alpha normalization. When a diffusion-map
normalization exponent is needed, write it as \(a\).

## 2. Basic Distance Vocabulary

A function \(d:X\times X\to[0,\infty]\) is an extended metric if it is
nonnegative, symmetric, zero exactly on the diagonal, satisfies the triangle
inequality, and may take value \(\infty\). It is a pseudometric if distinct
points may have distance zero. It is a directed distance or quasi-metric if the
triangle inequality holds but symmetry may fail. It is a graph metric when it
is the shortest-path distance induced by nonnegative edge weights on an
undirected connected graph.

This memo reserves the phrase **geodesic distance** for:

1. a continuous path metric obtained by minimizing a path length over curves;
2. a graph shortest-path metric intended to estimate such a continuous path
   metric;
3. an oracle path metric on a known synthetic support.

Diffusion, PHATE, commute, resistance, heat-kernel, and PageRank distances are
support-aware or random-walk distances, but should not be called geodesic
distances unless a specific small-time or shortest-path limit is being invoked.

## 3. Continuous Density-Weighted Geodesic

### Definition

For \(\alpha\ge 0\), define the density-weighted length of a piecewise \(C^1\)
curve \(\gamma:[0,1]\to M\) by

$$
L_{\rho,\alpha}(\gamma)
= \int_0^1
\rho(\gamma(t))^{-\alpha}
\sqrt{g_{\gamma(t)}(\gamma'(t),\gamma'(t))}\,dt.
$$

The associated density-weighted distance is

$$
d_{\rho,\alpha}(x,y)
= \inf_{\gamma(0)=x,\ \gamma(1)=y} L_{\rho,\alpha}(\gamma).
$$

When \(\alpha=0\), this reduces to the ordinary Riemannian geodesic distance
on \((M,g)\).

### Conformal Metric Interpretation

If \(\rho\) is positive and smooth, then \(d_{\rho,\alpha}\) is the geodesic
distance of the conformally transformed Riemannian metric

$$
g_{\rho,\alpha} = \rho^{-2\alpha}g.
$$

The factor of two appears because path length uses the square root of the
metric tensor:

$$
\sqrt{g_{\rho,\alpha}(\gamma',\gamma')}
= \rho^{-\alpha}\sqrt{g(\gamma',\gamma')}.
$$

### Metric Status

If \(M\) is connected and \(\rho\) is positive and finite everywhere, then
\(d_{\rho,\alpha}\) is a metric on \(M\). If \(\rho\) vanishes on barriers,
then distances crossing those barriers may be infinite; this is an extended
metric on each reachable support component. If \(\rho\) is unbounded, distances
between distinct points can collapse toward zero in degenerate cases; `gflow`
should avoid this by working with bounded or transformed densities.

### Practical Robust Version

For finite samples, use a transformed density

$$
\rho^{\star}(x)=T(\widehat{\rho}(x)),
$$

where \(T\) is a monotone positive transform. Important examples are:

$$
\rho^{\star}(x) = \max\{\widehat{\rho}(x),\rho_{\min}\},
$$

$$
\rho^{\star}(x)
= \min\{\max\{\widehat{\rho}(x),\rho_{\min}\},\rho_{\max}\},
$$

and a penalty cap

$$
\rho^{\star}(x)^{-\alpha}
= \min\{\widehat{\rho}(x)^{-\alpha},c_{\max}\}.
$$

These produce a metric for the transformed density. Diagnostics must report
how much of the data or edge set is controlled by the floor or cap.

## 4. Graph Density-Weighted Geodesic

### Definition

Let \(G=(V,E)\) be an undirected graph. For each edge \((i,j)\in E\), define

$$
w_{ij}
= \frac{\ell_{ij}}{\phi_{ij}^{\alpha}},
\qquad
\phi_{ij}>0,
$$

where \(\phi_{ij}\) summarizes density along or near the edge. The graph
density-weighted distance is

$$
d_{G,\phi,\alpha}(i,j)
= \min_{\pi:i\to j}\sum_{(u,v)\in\pi} w_{uv}.
$$

### Edge Density Summaries

Candidate summaries include

$$
\phi_{ij}^{\mathrm{arith}}
= \frac{\widehat{\rho}_i+\widehat{\rho}_j}{2},
$$

$$
\phi_{ij}^{\mathrm{geom}}
= \sqrt{\widehat{\rho}_i\widehat{\rho}_j},
$$

$$
\phi_{ij}^{\mathrm{harm}}
= \frac{2}{\widehat{\rho}_i^{-1}+\widehat{\rho}_j^{-1}},
$$

$$
\phi_{ij}^{\min}
= \min\{\widehat{\rho}_i,\widehat{\rho}_j\},
$$

and, when coordinates and an estimator are available,

$$
\phi_{ij}^{\mathrm{mid}}
= \widehat{\rho}\!\left(\frac{x_i+x_j}{2}\right).
$$

The minimum and harmonic summaries penalize low-density endpoints more
strongly. Midpoint or edge-neighborhood summaries are better for detecting
edges that cross gaps between two dense endpoint regions, but they require a
meaningful ambient or latent coordinate system.

### Metric Status

If \(G\) is undirected, connected, and all \(w_{ij}>0\), then
\(d_{G,\phi,\alpha}\) is a metric on vertices. If \(G\) is disconnected, it is
an extended metric with infinite distances between components. If zero-weight
edges are allowed, the result is a pseudometric. If \(w_{ij}\ne w_{ji}\) or
the graph is directed, the shortest-path distance is a directed distance, not a
metric.

### Relationship to the Continuous Metric

The graph distance estimates \(d_{\rho,\alpha}\) only when:

1. the graph topology is locally consistent with the intended support;
2. \(\ell_{ij}\) estimates local arclength at the edge scale;
3. \(\phi_{ij}\) estimates the density encountered along the edge;
4. graph resolution increases while edge lengths shrink.

Edge reweighting cannot fully fix a graph that contains unsupported shortcuts
unless the edge density summary sees the low-density gap crossed by the edge.

## 5. Fermat / Power-Weighted Sample Distance

### Rooted Convention

For \(p\ge 1\), define

$$
\ell_p(x,y;X)
= \min_{x=x_0,\ldots,x_m=y}
\left(
\sum_{r=0}^{m-1}\|x_{r+1}-x_r\|^p
\right)^{1/p}.
$$

### Unrooted Convention

Many topology papers use

$$
d_{X,p}(x,y)
= \min_{x=x_0,\ldots,x_m=y}
\sum_{r=0}^{m-1}\|x_{r+1}-x_r\|^p.
$$

These are related by \(d_{X,p}=\ell_p^p\) when the same admissible graph or
complete point cloud is used.

### Continuum Correspondence

For samples from an intrinsic \(d\)-dimensional density \(\rho\), the
continuum Fermat geodesic has integrand

$$
\rho^{-(p-1)/d}\,ds_g.
$$

Thus the `gflow` correspondence is

$$
\alpha = \frac{p-1}{d},
\qquad
p = 1+\alpha d.
$$

The corresponding conformal metric tensor is

$$
g_p = \rho^{-2(p-1)/d} g.
$$

### Metric Status

On a connected graph with positive Euclidean edge lengths, \(\ell_p\) is a
metric for \(p\ge 1\). The unrooted \(d_{X,p}\) is also a path metric because
it is a shortest-path distance with edge weights \(\|x_i-x_j\|^p\). In the
complete graph and \(p=1\), \(\ell_1\) collapses to ordinary Euclidean
distance. For \(p>1\), paths through dense regions are favored because they can
be decomposed into many short edges.

### `gflow` Interpretation

Fermat/PWSPD should be treated as both:

1. a well-studied comparator for density-aware graph geodesics;
2. a special implicit-density case of \(d_{G,\phi,\alpha}\), where no explicit
   density estimate is used and density enters through local sample spacing.

When reporting Fermat results, `gflow` should state the convention, \(p\), the
graph used for computation, any rescaling, and the implied
\(\alpha=(p-1)/\widehat{d}\) if an intrinsic dimension estimate is available.

## 6. Density-Aware Neighbor Selection

Density-aware distance can enter before edge weighting by changing graph
topology. Let \(r_k(i)\) be the distance from \(x_i\) to its \(k\)th nearest
neighbor. A kNN-radius density surrogate is

$$
\widehat{\rho}_i
\propto \frac{k}{n\,r_k(i)^{\widehat{d}}},
$$

where \(\widehat{d}\) is an intrinsic dimension estimate.

Possible topology rules include:

1. choose smaller effective \(k_i\) in dense regions and larger \(k_i\) in
   sparse regions, subject to a maximum edge length;
2. admit an edge only if it satisfies both a geometric length rule and a
   density-support rule;
3. connect mutual-kNN components using bridges ranked by geometric length and
   density support;
4. prune edges whose midpoint or edge-neighborhood density is too low.

These rules do not themselves define a distance. They define \(G\), after
which a graph metric is obtained from edge weights.

## 7. Diffusion Distance

### Definition

Given a kernel \(K_\varepsilon\), define a density estimate

$$
q_i = \sum_j K_\varepsilon(i,j).
$$

Diffusion-map alpha normalization uses

$$
K_a(i,j)
= \frac{K_\varepsilon(i,j)}{q_i^a q_j^a},
\qquad
P_a(i,j)
= \frac{K_a(i,j)}{\sum_k K_a(i,k)}.
$$

The diffusion distance at time \(t\) is

$$
D_t^2(i,j)
= \sum_k
\frac{(P_a^t(i,k)-P_a^t(j,k))^2}{\pi_k},
$$

where \(\pi\) is the stationary or reference measure used by the diffusion
construction.

### Density Role

The parameter \(a\) controls sampling-density normalization:

$$
a=0 \quad \text{keeps strong density drift},
$$

$$
a=1 \quad \text{removes sampling-density drift in the continuum limit}.
$$

This \(a\) is unrelated to the `gflow` density exponent \(\alpha\).

### Metric Status

Diffusion distance is a random-walk distributional distance. It is a metric on
states with distinct transition profiles under the chosen time and weighting;
it becomes a pseudometric if two points have identical diffusion profiles or
if the embedding truncates discriminating eigenvectors. It is not a
shortest-path geodesic.

### `gflow` Interpretation

Diffusion distance should be a comparator or diagnostic family. It answers:
"Do two points diffuse to similar neighborhoods?" It does not answer exactly:
"What is the cheapest supported path between two points?"

## 8. PHATE Potential Distance

PHATE constructs a Markov operator \(P\), chooses diffusion time \(t\), forms
\(P_t=P^t\), transforms rows of \(P_t\) into potentials, then compares those
potentials.

A simplified definition is

$$
U_{\gamma}(P_t(i,\cdot)) =
\begin{cases}
\log P_t(i,\cdot), & \gamma = 1,\\
\text{power or square-root-like potential}, & \gamma \approx 0,
\end{cases}
$$

and

$$
D_{\mathrm{PHATE}}(i,j)
= \left\|
U_{\gamma}(P_t(i,\cdot))
- U_{\gamma}(P_t(j,\cdot))
\right\|_2.
$$

PHATE is density-adjacent because high-density paths affect transition
probabilities, but the effect is mediated by graph construction, kernel
normalization, diffusion time, and the potential transform.

Metric status: PHATE potential distance is an embedding/potential distance,
usually a pseudometric if potential vectors coincide or if approximations
collapse points. It should not be labeled a geodesic distance.

## 9. Resistance, Commute, Heat-Kernel, and PageRank Distances

### Resistance and Commute

For graph Laplacian \(L\), effective resistance is

$$
R_{ij}
= (e_i-e_j)^T L^+(e_i-e_j).
$$

Commute time satisfies

$$
C_{ij} = \operatorname{vol}(G)R_{ij}.
$$

Resistance is a graph metric on connected undirected weighted graphs with
positive conductances. It rewards many parallel high-conductance paths. It is
not a shortest-path metric.

### Heat Kernel

The graph heat kernel is

$$
H_t = \exp(-tL).
$$

Heat-kernel distances compare rows of \(H_t\). On a smooth manifold, the
small-time heat kernel recovers geodesic distance through Varadhan's formula:

$$
d_g(x,y)^2
= -\lim_{t\to 0}4t\log H_t(x,y).
$$

Away from this limit, heat distances are diffusion-scale distances rather than
geodesic distances.

### Personalized PageRank

Personalized PageRank from source \(i\) solves

$$
r_i = \beta e_i + (1-\beta)r_iP.
$$

It is generally an asymmetric proximity. Symmetrized transformations may be
useful, but they are not automatically metrics.

## 10. Directed and Gradient-Flow Variants

For trajectory or gradient-flow analysis, `gflow` may eventually need
asymmetric path costs. A generic directed edge weight is

$$
w_{ij}^{\rightarrow}
= \frac{\ell_{ij}}{\phi_{ij}^{\alpha}}
+ \lambda\,B_{ij},
$$

where \(B_{ij}\) penalizes movement against a vector field, potential, time
arrow, or biological gradient. The resulting shortest-path distance

$$
d^{\rightarrow}(i,j)
= \min_{\pi:i\to j}\sum_{(u,v)\in\pi} w_{uv}^{\rightarrow}
$$

is generally not symmetric. It should be called a directed distance or
quasi-metric, not a metric.

This is out of scope for the first density-aware implementation but should be
kept conceptually separate from symmetric density-aware geodesics.

## 11. Density-Aware Oracle Truth

For synthetic data with known density \(\rho_{\mathrm{true}}\), define oracle
truth by

$$
d_{\mathrm{oracle},\alpha}(x,y)
= \inf_{\gamma:x\to y}
\int_{\gamma}
\rho_{\mathrm{true}}(\gamma(s))^{-\alpha}\,ds.
$$

In a dense oracle graph with oracle nodes \(z_a\), use

$$
w_{ab}^{\mathrm{oracle}}
= \frac{\|z_a-z_b\|}
{\rho_{\mathrm{true}}(m_{ab})^{\alpha}},
$$

where \(m_{ab}\) is a midpoint or edge-density quadrature point.

Oracle distances are not estimators; they define the benchmark target for a
chosen scientific interpretation of support. Reports must state \(\alpha\),
density floor, support truncation, grid resolution, and sample-attachment rule.

## 12. Recommended `gflow` Definition Set

For the next implementation and benchmark phase, use these names:

1. **density-weighted continuous geodesic**:
   \(d_{\rho,\alpha}\), defined by \(\int \rho^{-\alpha}\,ds\).
2. **density-weighted graph geodesic**:
   \(d_{G,\phi,\alpha}\), shortest paths with
   \(w_{ij}=\ell_{ij}/\phi_{ij}^{\alpha}\).
3. **Fermat/PWSPD distance**:
   shortest paths with edge powers \(\ell_{ij}^p\), reported with convention
   and implied \(\alpha=(p-1)/\widehat{d}\).
4. **diffusion distance**:
   row distance between \(P_a^t\) transition distributions.
5. **PHATE potential distance**:
   distance between transformed diffusion-potential rows.
6. **resistance/commute/heat/PageRank distances**:
   random-walk or operator distances, useful analogies and comparators but not
   density-weighted geodesics.
7. **density-aware oracle distance**:
   synthetic truth metric using known density and explicit support rules.

## 13. Minimal Diagnostics Required

Every report using \(d_{G,\phi,\alpha}\), Fermat/PWSPD, or an oracle
density-weighted distance should include:

1. density estimate distribution and transform;
2. floor/cap counts and fractions;
3. edge length versus density-weighted edge length;
4. density summary used for each edge;
5. sensitivity over \(\alpha\) or \(p\);
6. sample paths through dense and sparse regions;
7. underestimation and overestimation path attribution;
8. comparison to Euclidean, unweighted graph geodesic, Fermat/PWSPD,
   diffusion, and PHATE distances when available;
9. tail and outlier influence diagnostics;
10. stability across seeds and sample sizes.

## 14. Boundary Between Definition and Implementation

The definitions above do not decide:

1. whether density is estimated in raw, PCA, CLR/Aitchison, latent, graph, or
   oracle space;
2. which graph constructor is best;
3. which edge density summary is best;
4. which default \(\alpha\), floor, or cap should be used;
5. whether density should be exploited or corrected away for a given biology
   problem.

Those questions belong to the implementation design memo and benchmark reports.
This memo only fixes the mathematical language that those reports should use.

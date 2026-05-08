# Geodesic-Isometry Diagnostics Report Brief

This brief covers background, source assets, and proposed report tasks for a LaTeX report on geodesic-isometry diagnostics for data-derived graphs on quadratic surfaces. The report should serve as the diagnostic companion to the graph-construction reports: graph-construction reports explain how a graph is built from data, while this report explains how we judge whether the resulting graph metric recovers the surface geodesic geometry.

## Background

The current gflow development thread is focused on the data geodesic geometric reconstruction problem. Given a finite sample

\[
X=\{x_1,\ldots,x_n\}\subset \Gamma\subset \mathbb{R}^p,
\]

where \(\Gamma\) is a known quadratic surface or hypersurface, we construct a weighted graph \(G(X)\) on \(X\). The graph shortest-path metric is denoted

\[
D_G(i,j)=d_{G(X)}(x_i,x_j),
\]

and the numerical surface geodesic reference metric is denoted

\[
D_\Gamma(i,j)=d_\Gamma(x_i,x_j).
\]

The practical question is:

\[
D_G \approx D_\Gamma?
\]

That is, which data-to-graph construction makes \(X\), equipped with graph geodesic distances, as close as possible to the sampled surface geometry?

The current benchmark suite has evaluated graph families such as sKNN, iKNN, mKNN, fixed-radius, adaptive-radius, and continuous-kNN/cKNN on quadratic graph surfaces. Results so far suggest that adaptive-radius, cKNN, sKNN, and iKNN should be treated as first-class competitors in the next tests, while mKNN and fixed-radius can be kept as sparse sentinel baselines.

The next diagnostic report should not primarily introduce new graph constructors. Instead, it should define a principled diagnostic layer for graph geodesic recovery, with formulas, interpretation, and concrete examples.

## Connection to Geodesic MDS

The diagnostics are motivated by the Geodesic MDS project. The relevant conceptual distinction is between:

- graph geodesic distances in the input graph;
- the fixed shortest-path family used to realize those distances;
- embedded path lengths along those same fixed graph geodesics;
- reference manifold or surface geodesic distances.

In Geodesic MDS, for a fixed graph \(G=(V,E,w)\) and a fixed family of shortest paths

\[
\Gamma_G=\{\gamma^G_{ij}: i<j\},
\]

an embedding \(Z=(z_1,\ldots,z_n)\) is evaluated by embedded path lengths

\[
d_{\mathrm{emb},\Gamma_G}(i,j;Z)
=
\sum_{(u,v)\in\gamma^G_{ij}}
\|z_u-z_v\|_2.
\]

The Geodesic MDS stress compares these embedded path lengths to the input graph distances \(D_G(i,j)\). In our gflow benchmark setting, however, graph edges usually have ambient Euclidean edge lengths:

\[
w_{uv}=\|x_u-x_v\|_2.
\]

At the original surface coordinates \(Z=X\), the embedded length of each graph shortest path equals the graph shortest-path distance:

\[
d_{\mathrm{emb},\Gamma_G}(i,j;X)=D_G(i,j).
\]

Therefore ordinary Geodesic MDS stress against \(D_G\) is zero at the original coordinates and cannot diagnose whether \(G\) recovers the surface metric. The useful adaptation is to compare the graph path geometry to the surface reference:

\[
D_G(i,j)\quad\text{versus}\quad D_\Gamma(i,j).
\]

Equivalently, the current surface-target relative RMS metric can be interpreted as a cross-geodesic stress comparing the graph path metric with the surface geodesic metric.

## Core Diagnostic Definitions

The report should define the following diagnostics precisely.

### 1. Scale-Calibrated Relative Geodesic Stress

Let \(D_G\) be the graph geodesic distance matrix and \(D_\Gamma\) the reference surface geodesic distance matrix. Define the optimal scale

\[
\alpha^\ast
=
\arg\min_{\alpha>0}
\sum_{i<j}
\left(\alpha D_G(i,j)-D_\Gamma(i,j)\right)^2.
\]

The closed-form solution is

\[
\alpha^\ast
=
\frac{\sum_{i<j}D_G(i,j)D_\Gamma(i,j)}
       {\sum_{i<j}D_G(i,j)^2},
\]

assuming \(D_G\) is not identically zero on off-diagonal pairs.

The primary diagnostic is

\[
\sigma_\Gamma(G)
=
\left(
\frac{
\sum_{i<j}
\left(
\alpha^\ast D_G(i,j)-D_\Gamma(i,j)
\right)^2
}{
\sum_{i<j}D_\Gamma(i,j)^2
}
\right)^{1/2}.
\]

Interpretation: lower is better. Because the denominator is fixed by the reference surface metric, this is a relative RMS distance-matrix error and is comparable across datasets, with the usual caveat that some geometries are intrinsically harder.

### 2. Signed Residuals and Bias

Define the calibrated residual

\[
r_{ij}
=
\alpha^\ast D_G(i,j)-D_\Gamma(i,j).
\]

The signed bias is

\[
\operatorname{bias}(G)
=
\frac{
\frac{2}{n(n-1)}\sum_{i<j} r_{ij}
}{
\frac{2}{n(n-1)}\sum_{i<j} D_\Gamma(i,j)
}.
\]

Interpretation:

- \(\operatorname{bias}(G)>0\): calibrated graph geodesics are systematically too long; this often indicates a graph that is too sparse or path-constrained.
- \(\operatorname{bias}(G)<0\): calibrated graph geodesics are systematically too short; this often indicates shortcut edges.
- \(\operatorname{bias}(G)\approx0\): no strong signed bias, though error variance and tail failures may still be large.

Also define a shortcut fraction

\[
p_-(G)=
\frac{2}{n(n-1)}
\#\{(i,j): i<j,\ r_{ij}<0\}.
\]

This is the fraction of pairs whose calibrated graph geodesic distance is shorter than the surface reference.

### 3. Residual Tail Diagnostics

Define the relative absolute residual

\[
e_{ij}
=
\frac{|r_{ij}|}{D_\Gamma(i,j)}.
\]

The report should define and recommend reporting:

\[
Q_{50}(e),\quad Q_{90}(e),\quad Q_{95}(e),\quad Q_{99}(e).
\]

These tail metrics identify graph settings that look acceptable in RMS average but fail badly for a subset of pairs.

### 4. Scale-Regime Diagnostics

The Geodesic MDS manuscript treats \(k\) as a geometric scale parameter: small neighborhoods produce sparse, often over-long graph geodesics; large neighborhoods create shortcuts and approach the ambient Euclidean regime.

For a graph construction setting, define pairwise distortion ratios

\[
\eta_{ij}
=
\frac{\alpha^\ast D_G(i,j)}{D_\Gamma(i,j)}.
\]

Then bin pairs by \(D_\Gamma(i,j)\), for example:

- local pairs: lower third or lower quartile of \(D_\Gamma\);
- middle pairs: middle third or middle quantiles;
- long-range pairs: upper third or upper quartile.

For each band, report median signed residual or median distortion:

\[
\operatorname{median}_{(i,j)\in B}
\left(
\frac{r_{ij}}{D_\Gamma(i,j)}
\right)
\quad\text{or}\quad
\operatorname{median}_{(i,j)\in B}\eta_{ij}.
\]

Interpretation:

- local-band failure suggests local edge support or local edge-length geometry is wrong;
- long-range negative residuals suggest shortcuts;
- long-range positive residuals suggest excessive sparsity or bottlenecks;
- good local behavior with bad long-range behavior suggests the support graph has correct local geometry but wrong global topology.

### 5. Path-Level Diagnostics

Distance-matrix diagnostics do not reveal where graph shortest paths travel. For selected pairs \((i,j)\), especially worst residual pairs, compare the graph shortest path

\[
\gamma^G_{ij}=(i=v_0,v_1,\ldots,v_L=j)
\]

with a reference surface geodesic path \(\gamma^\Gamma_{ij}\) when available from the surface reference oracle.

Recommended diagnostics:

#### Tube Deviation

\[
\tau_{ij}
=
\frac{1}{L+1}
\sum_{s=0}^{L}
\operatorname{dist}_\Gamma(v_s,\gamma^\Gamma_{ij}).
\]

This measures whether the graph shortest path stays near the reference surface geodesic route.

#### Maximum Path Deviation

\[
\tau^{\max}_{ij}
=
\max_{0\le s\le L}
\operatorname{dist}_\Gamma(v_s,\gamma^\Gamma_{ij}).
\]

This catches graph paths that briefly wander far away.

#### Worst-Pair Path Examples

For each selected dataset and graph family, show:

- the largest positive residual pairs;
- the largest negative residual pairs;
- the graph shortest path on the surface;
- the reference geodesic path if available;
- the path-level diagnostic values.

These examples should be visual, not color-only. Use line thickness, point symbols, direct labels, and separate panels when possible.

### 6. Normal-Displacement Visualization

For 2D quadratic surfaces \(\Gamma\subset\mathbb{R}^3\), the report should define an experimental visualization based on normal displacement.

Let \(x_i\in\Gamma\) be the original sample point and let \(n_i\) be a consistently oriented unit normal. Define

\[
y_i(t)=x_i+t_i n_i.
\]

The proposed visualization optimizes scalar displacements \(t=(t_1,\ldots,t_n)\), while keeping points attached to their surface normals.

Using a selected pair set \(\mathcal P\), define

\[
L_{\mathrm{normal}}(t)
=
\sum_{(i,j)\in\mathcal P}
\rho_{ij}
\left(
\alpha
\sum_{(u,v)\in\gamma^G_{ij}}
\|y_u(t)-y_v(t)\|_2
-
D_\Gamma(i,j)
\right)^2
+
\lambda\sum_{i=1}^n t_i^2.
\]

Here:

- \(\gamma^G_{ij}\) is the fixed graph shortest path in the input graph;
- \(D_\Gamma(i,j)\) is the surface reference distance;
- \(\rho_{ij}\) is a pair weight, for example

\[
\rho_{ij}=\frac{1}{D_\Gamma(i,j)^2+\epsilon};
\]

- \(\lambda\) controls how strongly the visualization resists moving points away from the original surface.

The displacement magnitudes \(|t_i|\) form a geometric distortion field:

- small displacements suggest the graph metric is already compatible with the surface geometry;
- large structured displacement suggests systematic metric distortion;
- localized spikes suggest local shortcut or bottleneck failures.

This visualization is explicitly diagnostic and supervised by the known surface geometry. It should not be presented as an unsupervised layout algorithm.

## Proposed Example Structure

The report should follow the style of the existing graph-construction reports: start with a self-contained mathematical explanation, then illustrate the concepts on small and medium examples.

Recommended examples:

1. **Toy graph examples**
   - A path graph with correct edge lengths.
   - A path graph with one shortcut edge.
   - A sparse graph that forces detours.
   These examples can be schematic and do not need full quadratic-surface machinery.

2. **2D paraboloid**

\[
z=x^2+y^2.
\]

Use one or two graph settings:
   - a near-optimal adaptive-radius or cKNN graph;
   - an intentionally sparse sKNN/iKNN graph;
   - optionally a shortcut-heavy high-parameter graph.

3. **2D saddle**

\[
z=x^2-y^2.
\]

Repeat the same diagnostic layout to show whether residual signs and distance-band behavior differ from the paraboloid.

4. **Worst-pair path panels**
   For at least one paraboloid and one saddle dataset, show paths for:
   - largest positive residual;
   - largest negative residual;
   - a representative low-error pair.

5. **Normal-displacement visualization**
   Include this as a design/prototype section. If implementation is not yet complete, include formulas and a mock/schematic figure. If implementation is feasible, include one small example.

## Existing Assets

### gflow Package Assets

- Main repo: `/Users/pgajer/current_projects/gflow`
- Quadratic geodesic utilities: `/Users/pgajer/current_projects/gflow/R/quadform_geodesics.R`
- Deviation helpers: `/Users/pgajer/current_projects/gflow/R/isometry_deviation.R` if present, otherwise search `R/` for `summarize.isometry.deviation`
- Graph geodesic wrapper: `/Users/pgajer/current_projects/gflow/R/graph_geodesic_distances.R`
- Graph constructors:
  - `/Users/pgajer/current_projects/gflow/R/sknn_graphs.R`
  - `/Users/pgajer/current_projects/gflow/R/mknn_graphs.R`
  - `/Users/pgajer/current_projects/gflow/R/iknn_graphs.R`
  - `/Users/pgajer/current_projects/gflow/R/radius_graphs.R`

### Existing Reports and Notes

- Graph construction report source:
  `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/phate_knn_graph_constructions.tex`
- Graph construction report PDF:
  `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/build/phate_knn_graph_constructions.pdf`
- First quadform benchmark design:
  `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- Benchmark HTML report spec:
  `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`
- Current Tier 2 benchmark report:
  `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-tier2-domain-sampling-benchmark/runs/full/report/quadform_tier2_domain_sampling_benchmark_report.html`
- Surface-only unpruned benchmark report:
  `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-surface-unpruned-benchmark/report/quadform_surface_unpruned_benchmark_report.html`

### Geodesic MDS Assets

- Main geodesic MDS manuscript:
  `/Users/pgajer/current_projects/geodesic_MDS/manuscript/geodesic_mds_v3.tex`
- Built PDF:
  `/Users/pgajer/current_projects/geodesic_MDS/build/geodesic_mds_v3.pdf`
- Diagnostic GMDS report:
  `/Users/pgajer/current_projects/geodesic_MDS/gmds_diagnostic_report.tex`
- Diagnostic GMDS PDF:
  `/Users/pgajer/current_projects/geodesic_MDS/gmds_diagnostic_report.pdf`
- MISF note:
  `/Users/pgajer/current_projects/geodesic_MDS/notes/misf_for_GMDS.md`

Key sections in the geodesic MDS manuscript:

- Introduction and output-side inconsistency.
- Self-consistency principle.
- Triangle inequality bias.
- Geodesic stress definition.
- \(k\)-regime / geometric scale parameter.
- Current implementation and diagnostic lessons.

## Recommended Deliverables

Create a LaTeX report under a new report directory, for example:

```text
dev/data-geodesic-reconstruction/geodesic-isometry-diagnostics-report/
```

Suggested files:

```text
geodesic_isometry_diagnostics.tex
geodesic_isometry_diagnostics.bib
build/geodesic_isometry_diagnostics.pdf
figures/
scripts/
```

The first milestone should produce a mathematically complete report even if the examples are schematic. A second milestone can add generated examples and figures.

## Suggested Report Title

Possible title:

```text
Geodesic-Isometry Diagnostics for Data-Derived Graphs on Quadratic Surfaces
```

Alternative shorter title:

```text
Graph Geodesic Recovery Diagnostics
```

## Open Design Questions

The report agent should explicitly address these rather than silently choosing:

1. Should the primary diagnostic be called `rel_geodesic_stress`, `surface_geodesic_stress`, or keep the current `rel_rms_error` name while explaining the interpretation?
2. Should pair-distance bands use tertiles, quartiles, or fixed distance thresholds?
3. Should residual quantiles be computed on \(|r_{ij}|/D_\Gamma(i,j)\) or on \(|r_{ij}|/\operatorname{median}(D_\Gamma)\) to reduce instability for tiny distances?
4. For path-level diagnostics, what reference path payload is currently available from the gflow surface oracle?
5. Should the first report include implemented normal-displacement examples, or define the method and leave implementation to a later agent?

## Recommended Implementation Order

1. Read this brief and the listed Geodesic MDS sections.
2. Inspect the current gflow quadform benchmark utilities and reports.
3. Produce an action plan before editing.
4. Build the LaTeX report skeleton with all mathematical definitions.
5. Add schematic toy examples if fast.
6. Add quadratic-surface examples if the data/assets are readily reusable.
7. Build the PDF and fix LaTeX issues.
8. Summarize remaining implementation gaps, especially normal-displacement visualization.


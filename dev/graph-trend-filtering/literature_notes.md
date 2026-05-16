# Literature Notes: Graph Trend Filtering

These notes summarize the initial paper set for implementing graph trend
filtering in `gflow`. The goal is not to reproduce each paper, but to extract
the operator definitions, solver implications, and design questions that matter
for an R/C++ package implementation.

PDFs are stored locally in `papers/` for convenience and are ignored by git.

## 1. Wang, Sharpnack, Smola, Tibshirani: Trend Filtering on Graphs

Local PDF:

- `papers/wang_sharpnack_smola_tibshirani_2016_trend_filtering_on_graphs.pdf`

Links:

- arXiv: <https://arxiv.org/abs/1410.7690>
- JMLR PDF: <https://jmlr.csail.mit.edu/papers/volume17/15-147/15-147.pdf>

### Core Estimator

For graph \(G = (V, E)\), observations \(y \in \mathbb R^n\), and graph signal
\(\beta \in \mathbb R^n\), graph trend filtering solves

```math
\hat\beta =
\arg\min_{\beta \in \mathbb R^n}
\frac{1}{2}\|y-\beta\|_2^2
+
\lambda \|\Delta^{(k+1)}\beta\|_1.
```

The important shift from Laplacian smoothing is replacing a quadratic
roughness penalty with an \(\ell_1\) penalty on graph differences.

### Difference Operators

Let \(D\) be an oriented incidence matrix and \(L = D^\top D\) be the graph
Laplacian. In the unweighted construction,

```math
\Delta^{(1)} = D.
```

For higher orders, the paper defines graph difference operators recursively by
alternating incidence and Laplacian-style operations. One convenient expression
is:

```math
\Delta^{(k+1)}
=
\begin{cases}
L^{(k+1)/2}, & k \text{ odd},\\
D L^{k/2}, & k \text{ even}.
\end{cases}
```

Interpretation:

- \(k=0\): graph fused lasso / graph total variation; piecewise constant.
- \(k=1\): \(\|\Delta^{(2)}\beta\|_1 = \|L\beta\|_1\); graph analogue of
  piecewise linearity.
- \(k=2\): \(\|D L\beta\|_1\); edge differences of graph Laplacian values,
  graph analogue of piecewise quadratic behavior.

### Implementation Implications for gflow

This is the primary blueprint for `fit.graph.trend.filtering()`.

Likely phase order:

1. Implement \(k=0\): weighted graph fused lasso.
2. Add operator construction diagnostics for \(D_w\), \(L_w\), and
   \(\Delta_w^{(k+1)}\).
3. Add \(k=1\) and \(k=2\) only after the weighting convention and solver path
   are stable.

The paper emphasizes local adaptivity relative to Laplacian smoothing, which
fits exactly with `gflow`'s comparator goal: graph trend filtering should be
presented as an \(\ell_1\)-adaptive comparator to metric graph low-pass
regression, not as a replacement.

### Open Questions

- What should the weighted incidence convention be?
  - \(D_w[e, i] = -w_e,\; D_w[e, j] = w_e\)
  - \(D_w[e, i] = -\sqrt{c_e},\; D_w[e, j] = \sqrt{c_e}\)
  - another length-scaled convention?
- Should weights represent conductances directly, or should metric lengths be
  transformed into conductances first, as in `fit.metric.graph.lowpass()`?
- Should disconnected graphs be supported by fitting each component separately,
  or rejected initially?

## 2. Madrid Padilla, Sharpnack, Scott, Tibshirani: The DFS Fused Lasso

Local PDF:

- `papers/madrid_padilla_sharpnack_scott_tibshirani_2018_dfs_fused_lasso.pdf`

Link:

- JMLR page: <https://www.jmlr.org/beta/papers/v18/16-532.html>

### Core Focus

This paper focuses on the graph fused lasso / graph total variation problem:

```math
\hat\beta =
\arg\min_{\beta \in \mathbb R^n}
\frac{1}{2}\|y-\beta\|_2^2
+
\lambda \sum_{(i,j)\in E} |\beta_i-\beta_j|.
```

This is the \(k=0\) graph trend filtering case.

### Main Computational Idea

The paper shows that a depth-first-search traversal induces a chain graph whose
total variation is controlled by the original graph total variation. This
allows a simple 1D fused lasso computation on the DFS chain to achieve strong
statistical guarantees for graph denoising.

### Implementation Implications for gflow

This paper suggests two possible implementation paths:

1. Exact graph fused lasso solver, for correctness and smaller graphs.
2. DFS-chain approximation, for large graphs or benchmark comparison.

The exact solver should be the first target if feasible, because it gives a
clear mathematical operator. The DFS method could later be exposed as an
approximate/fast mode.

### Open Questions

- Should `gflow` expose a DFS approximation as a solver mode?
- If yes, should it be framed as an approximation to graph fused lasso or as a
  separate estimator?
- How should edge weights affect DFS traversal and chain weights?

## 3. Pastukhov: Fused \(\ell_1\) Trend Filtering on Graphs

Local PDF:

- `papers/pastukhov_2024_fused_l1_trend_filtering_on_graphs.pdf`

Link:

- arXiv: <https://arxiv.org/abs/2401.05250>

### Core Focus

This paper studies fused trend filtering on general graphs, combining graph
trend filtering with additional fusion regularizers such as anisotropic total
variation and nearly-isotonic restrictions. It also discusses computational
approaches with per-iteration complexity linear in the number of edges.

### Implementation Implications for gflow

This paper is useful after the base graph trend filtering API is working. It
points toward richer operators:

```math
\frac{1}{2}\|y-\beta\|_2^2
+
\lambda_{\mathrm{tf}}\|\Delta^{(k+1)}\beta\|_1
+
\lambda_{\mathrm{fuse}}\|D\beta\|_1.
```

This may be attractive for SIMODS-style experiments if we want both
higher-order graph smoothness and explicit edge-wise fusion.

### Open Questions

- Should fusion and trend-filtering penalties share one lambda or use two
  separate tuning parameters?
- Is this too broad for the first `gflow` implementation?
- Should nearly-isotonic constraints be kept out of scope initially?

## Proposed Near-Term gflow Scope

Start with a deliberately narrow API:

```r
fit.graph.trend.filtering(
  adj.list,
  weight.list = NULL,
  y,
  order = 0L,
  lambda.grid = NULL,
  lambda.selection = c("cv", "fixed"),
  weight.rule = c("conductance", "sqrt.conductance", "unit"),
  ...
)
```

The first landing patch should focus on:

- graph canonicalization;
- weighted oriented incidence construction;
- \(k=0\) graph fused lasso objective;
- tiny exact reference tests;
- path-graph agreement checks against known 1D fused-lasso/trend-filtering
  behavior where possible;
- clear return object structure;
- documentation that contrasts \(\ell_1\)-adaptive graph trend filtering with
  \(\ell_2\) metric graph low-pass smoothing.

Higher-order \(k=1,2\) graph trend filtering should be a second phase unless
the solver design makes them almost free once \(k=0\) is implemented.

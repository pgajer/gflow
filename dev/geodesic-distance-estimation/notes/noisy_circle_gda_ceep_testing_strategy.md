# Noisy Circle GDA/CEEP Testing Strategy

## Purpose

This document describes a benchmark strategy for studying geodesic distance
approximation (GDA) in noisy circle examples while keeping the downstream
conditional expectation estimation problem (CEEP) in view.

The central question is:

```text
Which graph construction, pruning level, and component policy gives the best
approximation to the true latent circular geodesic distance?
```

The benchmark should allow iKNN, mKNN, repaired local-neighbor graphs, or
MST/CMST-derived graphs to win. The target is not method identity. The target
is recovery of the true data geometry.

## GDA and CEEP

For GDA, disconnectedness is a geometric outcome and should be measured. A
graph with many components does not define finite graph distances for all pairs.
That matters if the goal is a global geodesic metric.

For CEEP, however, the graph is often used to estimate or smooth a response on
well-supported regions of the data. Very small connected components may be
better treated as outliers or low-support regions. In that context, we can
reasonably:

- filter tiny components,
- perform CEE on the major connected component,
- or join only large components when there is local-neighbor or MST evidence
  that they should belong to the same geometric object.

This creates two layers of decisions.

## Decision Layer 1: Base Local Geometry

The first layer asks whether iKNN or mKNN gives a better raw approximation to
the true latent geodesic distance.

### iKNN

iKNN tends to be more connected and often gives good circular geodesic
recovery at moderate `k`, but it can become too dense. As `k` grows it may add
edges that shorten graph paths more than the true manifold geometry supports.

### mKNN

mKNN is more conservative. This can reduce spurious edges, but it often has
more connected components. From a pure GDA perspective this is a liability.
From a CEEP perspective it may be acceptable if small components are filtered
or if only large components are repaired.

## Decision Layer 2: Graph Modification

The second layer asks how to modify the base graph.

### Component Policy

The benchmark should compare at least:

1. **All components**: evaluate the graph as constructed.
2. **Major component only**: filter to the largest connected component.
3. **Join large components using MST edges**: add only enough MST-supported
   edges to connect large components, leaving small components as outliers.
4. **Join large components using iKNN edges**: add only enough iKNN-supported
   edges to connect large mKNN components when iKNN provides local-neighbor
   evidence.

The threshold for a large component should be explicit, for example:

```text
large component size >= max(8, ceiling(0.05 * n))
```

This threshold should be treated as a benchmark parameter, not a hidden
constant.

### Geometric Pruning

The benchmark should sweep geometric pruning parameters for both iKNN and
mKNN.

For iKNN, the relevant gflow controls are:

```r
create.iknn.graphs(
  X,
  kmin = k,
  kmax = k,
  max.path.edge.ratio.deviation.thld = delta,
  path.edge.ratio.percentile = p,
  threshold.percentile = q,
  compute.full = TRUE,
  pca.dim = NULL
)
```

`max.path.edge.ratio.deviation.thld` is a deviation threshold. Internally,
the path-to-edge ratio threshold is `1 + delta`. Larger values retain more
near-redundant edges; smaller positive values prune more aggressively. Setting
`delta = 0` disables geometric pruning in this benchmark.

For mKNN, the corresponding gflow controls are:

```r
create.mknn.graphs(
  X,
  kmin = k,
  kmax = k,
  max.path.edge.ratio.thld = delta,
  path.edge.ratio.percentile = p,
  compute.full = TRUE,
  pca.dim = NULL
)
```

In this development benchmark we treat `max.path.edge.ratio.thld` as the
public pruning amount and use `0` to disable pruning.

## Noisy Circle Grid

The initial executable benchmark should be a pilot that is large enough to
show patterns but small enough to iterate.

Recommended pilot grid:

```text
n in {80, 160, 320}
seed in {4101, 4102, 4103}
k in {4, 6, 8, 10}
geometric pruning delta in {0, 0.05, 0.10}
path.edge.ratio.percentile = 0.5
large component threshold = max(8, ceiling(0.05 * n))
```

The grid can be expanded once the report is stable.

## Metrics

For each candidate graph, report:

- sample size `n`,
- seed,
- base graph type: iKNN or mKNN,
- `k`,
- pruning `delta`,
- component policy,
- number of full-graph components,
- largest component fraction,
- evaluation scope: all vertices, largest component, or large components,
- number and fraction of evaluated vertices,
- graph-distance finite-pair fraction in the evaluation scope,
- Spearman and Pearson correlation with true circular geodesics,
- median and 90th percentile relative geodesic error after optimal global
  scaling,
- mean true-neighbor overlap,
- edge count,
- false shortcut rate.

The primary GDA metric should be median relative error, with q90 relative error
and finite-pair fraction as important safeguards.

## Primary Comparisons

The report should answer:

1. Does iKNN or mKNN produce better raw GDA at the same `k` and pruning level?
2. Does geometric pruning improve iKNN or mKNN GDA?
3. Does mKNN mainly fail because of disconnectedness, or because of poor
   within-component geodesic geometry?
4. Does repairing large mKNN components with MST edges improve GDA?
5. Does repairing large mKNN components with iKNN edges improve GDA more than
   MST repair?
6. Does filtering to the major component yield better CEEP-relevant geometry
   than forcing global connectedness?

## Reporting Style

The report should be figure-first:

- best mean GDA error by sample size,
- pruning curves averaged across seeds,
- component-policy comparison,
- connectivity and evaluated-vertex fraction,
- small best-candidate tables by sample size.

Complete candidate tables should be included only for reproducibility or as
collapsed detail.

## Interpretation Principle

The benchmark should not assume that connectedness is always good.

For GDA, disconnectedness limits global metric recovery.

For CEEP, disconnectedness can be acceptable if small components are treated as
outliers and the major component supports stable conditional expectation
estimation.

The correct graph is therefore the one that best serves the stated objective:
global geodesic reconstruction, major-component CEE, or sparse connected
geometry for downstream modeling.

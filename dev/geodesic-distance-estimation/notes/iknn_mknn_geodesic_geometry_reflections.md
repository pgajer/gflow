# Reflections on iKNN, mKNN, MST, and True Geodesic Recovery

## Context

The noisy circle benchmark has clarified an important point: the goal is not to
make an MST-completion graph look better than iKNN or mKNN as a matter of
method identity. The goal is to recover the latent geodesic geometry of the
data. In this setting, iKNN and mKNN are not merely competitors to CMST; they
are useful local-geometry estimators that may already contain better evidence
about the manifold than the MST backbone.

For the noisy circle, plain MST and conservative CMST tend to remain too
path-like. They preserve connectedness, but they do not necessarily recover the
circle topology. In contrast, tuned iKNN and mKNN graphs often recover circular
geodesic distances much better, although they introduce their own problems:
iKNN can become too dense, and mKNN can fragment into too many connected
components.

This suggests a reframing.

## Reframed Question

Instead of asking:

```text
Can geodesic distance estimation better approximate true geodesic geometry than iKNN?
```

we should ask:

```text
Given iKNN/mKNN as imperfect local-geometry estimators,
can their graph structure or graph metric be modified to better approximate
the true latent geodesic geometry?
```

Under this framing, MST/CMST is one possible scaffold or diagnostic device, not
the privileged source of geometry.

## What Does MST Add to iKNN or mKNN?

MST does not usually add better local geometry to iKNN or mKNN. If iKNN or mKNN
already gives a connected graph whose shortest-path metric approximates the
latent geodesic metric, simply adding MST edges can be neutral or harmful. It
can shorten paths through tree artifacts that are not true manifold
adjacencies.

MST is useful mainly in three roles:

1. **Connectedness scaffold**: it guarantees finite graph distances when a
   conservative local graph is disconnected.
2. **Sparse baseline**: it supplies a minimal connected graph that can be
   repaired selectively rather than replaced by a dense local-neighbor graph.
3. **Failure diagnostic**: large discrepancies between an MST path distance and
   an alternative local-neighbor graph distance can identify where the tree has
   stretched the manifold geometry too much.

Thus, an MST+iKNN strategy should not be interpreted as "MST improves iKNN" by
default. A better interpretation is:

```text
Can iKNN evidence repair a sparse connected MST/CMST graph enough to recover
geodesic structure?
```

This is a valid goal when sparsity, connectedness, or interpretability matters.
If those constraints are not important, a well-tuned iKNN or mKNN graph may be
the cleaner estimator.

## Detour-Based Repair

One proposed repair idea compares MST geodesic distance to a local comparator
distance.

The useful signal is not:

```text
d_MST(i, j) <= c * d_comparator(i, j)
```

because that condition identifies pairs whose MST path is already reasonably
short and therefore mostly redundant.

The more useful signal is a large detour:

```text
d_MST(i, j) >= c * d_comparator(i, j)
```

or equivalently:

```text
d_MST(i, j) / d_comparator(i, j) >= c
```

For a Euclidean comparator, this can detect missing circle-closure edges, but
it is dangerous on crescents, spirals, figure-eights, and nearby non-touching
arcs, where Euclidean-near points can be far apart in the true geodesic metric.

For an iKNN comparator, the idea is more defensible:

```text
i and j are in the same iKNN connected component
and
d_MST(i, j) / d_iKNN(i, j) >= c
```

This says that a local-neighbor graph already supports a relatively short
within-graph route, while the MST forces a long detour. That is evidence that
the MST backbone may be missing a topological closure or alternate path.

However, this still does not prove that MST+iKNN is better than iKNN. It only
supports a sparse-repair strategy:

```text
Use MST/CMST as a sparse connected scaffold.
Use iKNN/mKNN as evidence for where the scaffold needs repair.
Evaluate the final shortest-path metric against known truth.
```

## Directly Improving iKNN and mKNN

The noisy circle results suggest that the next development step should include
direct modifications of iKNN and mKNN, not only MST-completion variants.

### 1. Connected mKNN by Local Bridges

mKNN is conservative and often geometrically cleaner than iKNN, but it can have
too many connected components. Instead of unioning with the full MST, we can
connect mKNN components using only the shortest locally plausible
inter-component bridges.

Questions to test:

- Does component repair preserve the low-shortcut behavior of mKNN?
- Does it outperform full `MST union mKNN`?
- Are the added bridges concentrated in true sampling gaps or do they create
  false shortcuts?

### 2. iKNN or mKNN Shortcut Pruning

iKNN can become too dense as `k` increases. Some edges may be locally short in
Euclidean distance but redundant or harmful for graph geodesics.

A pruning rule can remove an edge when:

```text
the graph remains connected
and
there is an alternate path between the endpoints whose length is not much
larger than the edge length
```

This may preserve circular topology while reducing chord-like shortcuts.

### 3. Adaptive k

A fixed `k` is sensitive to nonuniform sampling. Dense regions may be
over-connected, while sparse regions may remain under-connected.

An adaptive-k graph could use:

- smaller `k` in dense regions,
- larger `k` in sparse regions,
- local scale estimates from nearest-neighbor distances or MST incident edge
  lengths,
- constraints that prevent high-k sparse-region repairs from becoming long
  ambient shortcuts.

### 4. Geodesic-Calibrated Edge Weights

Even if the topology is good, raw Euclidean edge weights may underestimate true
arc length, especially on noisy curved manifolds. We should test whether the
same iKNN/mKNN topology improves when edge weights are locally calibrated.

Possible weight modifications:

- Euclidean edge weights with local scale correction,
- weights inflated by local curvature or neighbor-spacing asymmetry,
- weights derived from short local graph paths rather than direct chords,
- post-hoc global scaling plus local residual diagnostics.

### 5. Cycle-Aware Closure

For circle-like data, a graph can be connected but topologically path-like. A
cycle-aware repair should add only a small number of closure edges supported by
local-neighbor evidence.

This should be tested cautiously. A rule that works on a circle can fail badly
on a crescent or spiral. The benchmark ladder must therefore include examples
where the correct topology is open, folded, branching, or separated.

## Decision Rule for Final Graph Selection

A practical benchmark should allow iKNN/mKNN to win. The decision rule should
be geometry-first:

```text
If iKNN or mKNN is connected, stable across seeds, sparse enough for the
application, and has lower geodesic distortion than MST-repair variants,
prefer iKNN/mKNN.

If iKNN or mKNN is disconnected, too dense, unstable, or creates shortcuts,
use it as evidence to repair a sparser connected graph.
```

This prevents the evaluation from becoming biased toward CMST as the intended
answer. CMST is one construction in a broader search for graph metrics that
recover latent geodesic structure.

## Proposed Next Report

The next report should be explicitly titled around modified iKNN/mKNN geometry,
for example:

```text
Noisy Circle: Can Modified iKNN/mKNN Improve True Geodesic Recovery?
```

It should compare:

- plain iKNN across `k`,
- plain mKNN across `k`,
- `MST union iKNN(k)`,
- `MST union mKNN(k)`,
- component-repaired mKNN,
- shortcut-pruned iKNN,
- shortcut-pruned mKNN,
- adaptive-k iKNN/mKNN,
- iKNN/mKNN with calibrated edge weights,
- CMST and local-adaptive CMST as references,
- oracle weighted circle graph as the target reference.

The primary outcome should be true geodesic metric recovery, summarized across
sample sizes and random seeds. The main communication should be figures and
small best-method tables, with exhaustive candidate tables included only as
collapsible or cached completeness material.

## Working Hypothesis

For noisy circles, the strongest candidates may be:

- mKNN with careful component repair,
- iKNN with shortcut pruning,
- local-adaptive MST/CMST in moderate sample sizes,
- high-q CMST with pruning at larger sample sizes.

The hypothesis should remain provisional until tested on examples that punish
false Euclidean shortcuts, especially crescents, spirals, nearby non-touching
arcs, figure-eights, and compositional cyclic gradients.

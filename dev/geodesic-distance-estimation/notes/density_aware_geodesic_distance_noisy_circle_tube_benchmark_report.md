# Density-Aware Geodesic Distance Noisy-Circle / Tube Benchmark Report

<style>
table { font-size: 10px; width: 100%; }
th, td { padding: 2px 4px; }
img { max-width: 100%; }
code { white-space: pre-wrap; }
</style>

Date: 2026-05-04

This report is Deliverable 4 for the density-aware geodesic-distance research track. It turns the noisy-circle latent-tube oracle from a sanity check into a benchmark with density sensitivity, floor/cap/attachment sweeps, distance-family comparisons, and edge/path attribution.

## 1. Benchmark Setup

The observed data are generated as a noisy circle with radius 1, radial noise scale \(\sigma=0.08\), sample size \(n=160\), and seed `4101`. Each sample has latent tube coordinates \((\theta_i,u_i)\), where

$$
u_i = \|x_i\| - r.
$$

The density-aware oracle lives on a structured periodic grid in \((\theta,u)\). For an oracle edge \((a,b)\), the baseline weight is

$$
w_{ab}=\|z_a-z_b\|\,\min\{\rho(m_{ab})^{-\alpha}, c_{\max}\},
$$

where \(m_{ab}\) is the edge midpoint in radial coordinate, \(\rho(u)=\max\{\exp[-(u/\sigma)^2],\rho_{\min}\}\), and \(c_{\max}=\infty\) unless a penalty cap is specified.

Primary truth for method comparison is the tube oracle with \(\alpha=0.5\), \(\rho_{\min}=1e-04\), no penalty cap, and `attach_extra = 1`.

![Observed geometry and latent tube coordinates](../figures/noisy-circle-tube-benchmark/01_geometry_density.png)

## 2. Oracle Sensitivity

The first question is whether the oracle itself is stable enough to serve as truth. The table below varies \(\alpha\) while holding the density floor and attachment rule fixed.

alpha | median | q90 | max | rho.arc | rho.euc
--- | --- | --- | --- | --- | ---
0 | 1.430 | 2.315 | 2.749 | 0.9943 | 0.9927
0.2500 | 1.633 | 2.903 | 3.784 | 0.9916 | 0.9712
0.5000 | 1.737 | 3.032 | 7.335 | 0.9407 | 0.9209
0.7500 | 1.858 | 3.230 | 33.62 | 0.8562 | 0.8400
1.000 | 1.959 | 4.721 | 236.4 | 0.7589 | 0.7467

At \(\alpha=0.5\), median oracle distance is 1.737 and Spearman correlation with latent arc distance is 0.9407. At \(\alpha=1\), maximum oracle distance rises to 236.4 and rank agreement with latent arc drops to 0.7589. This repeats the earlier warning: \(\alpha=1\) is tail dominated unless floors, caps, or attachment rules are handled carefully.

![Oracle sensitivity to alpha, floors/caps, and attachment](../figures/noisy-circle-tube-benchmark/02_oracle_sensitivity.png)

### Floor And Cap Sweep

The floor/cap sweep below focuses on \(\alpha=0.5\) and \(\alpha=1\). Rows are scored against the primary \(\alpha=0.5\) oracle, not against their own variant, so high relative error for \(\alpha=1\) means that the strong-density truth is a different object.

alpha | floor | cap | median | q90 | max | rho.primary | med.err
--- | --- | --- | --- | --- | --- | --- | ---
0.5000 | 0.0001000 | Inf | 1.737 | 3.032 | 7.335 | 1.000 | 0
1.000 | 0.0001000 | Inf | 1.959 | 4.721 | 236.4 | 0.9029 | 0.9587
0.5000 | 0.001000 | Inf | 1.736 | 3.027 | 6.017 | 0.9995 | 0.01859
1.000 | 0.001000 | Inf | 1.959 | 4.721 | 86.32 | 0.9029 | 0.8422
0.5000 | 0.01000 | Inf | 1.717 | 2.993 | 4.752 | 0.9919 | 0.03823
1.000 | 0.01000 | Inf | 1.959 | 4.433 | 17.99 | 0.9051 | 0.3868
0.5000 | 0.0001000 | 20.00 | 1.734 | 3.015 | 5.492 | 0.9979 | 0.02719
1.000 | 0.0001000 | 20.00 | 1.919 | 3.290 | 7.358 | 0.9662 | 0.09119
0.5000 | 0.001000 | 20.00 | 1.734 | 3.015 | 5.492 | 0.9979 | 0.02719
1.000 | 0.001000 | 20.00 | 1.919 | 3.290 | 7.358 | 0.9662 | 0.09119
0.5000 | 0.01000 | 20.00 | 1.717 | 2.993 | 4.752 | 0.9919 | 0.03823
1.000 | 0.01000 | 20.00 | 1.919 | 3.290 | 7.358 | 0.9662 | 0.09119
0.5000 | 0.0001000 | 10.00 | 1.717 | 2.993 | 4.752 | 0.9919 | 0.03823
1.000 | 0.0001000 | 10.00 | 1.850 | 3.132 | 5.524 | 0.9864 | 0.03024
0.5000 | 0.001000 | 10.00 | 1.717 | 2.993 | 4.752 | 0.9919 | 0.03823
1.000 | 0.001000 | 10.00 | 1.850 | 3.132 | 5.524 | 0.9864 | 0.03024
0.5000 | 0.01000 | 10.00 | 1.717 | 2.993 | 4.752 | 0.9919 | 0.03823
1.000 | 0.01000 | 10.00 | 1.850 | 3.132 | 5.524 | 0.9864 | 0.03024

### Attachment Sweep

attach | median | q90 | max | rho.primary | med.err
--- | --- | --- | --- | --- | ---
0 | 1.745 | 3.042 | 7.373 | 0.9999 | 0.002894
1.000 | 1.737 | 3.032 | 7.335 | 1.000 | 0
2.000 | 1.730 | 3.025 | 7.137 | 1.000 | 0.004558

Attachment changes are modest for \(\alpha=0.5\) in this case, but attachment remains a required diagnostic because tail samples can otherwise inherit an artificial radial cost.

## 3. Distance-Family Comparison

The comparison includes Euclidean distance, latent arc distance, MST, CMST, iKNN, mKNN, density-weighted iKNN, Fermat/PWSPD on iKNN, diffusion row distance, and a PHATE-like log-potential distance. Diffusion and PHATE-like distances are included as support-aware comparators, not as geodesic distances.

![Distance-family comparison](../figures/noisy-circle-tube-benchmark/03_method_comparison.png)

distance | finite | rho.primary | med.err | q90.err | rho.arc | rho.euc
--- | --- | --- | --- | --- | --- | ---
iKNN density alpha=0.5 | 1.000 | 0.9805 | 0.02985 | 0.1673 | 0.9772 | 0.9576
iKNN density alpha=1 capped | 1.000 | 0.9861 | 0.04925 | 0.1723 | 0.9398 | 0.9227
mKNN base | 1.000 | 0.9507 | 0.07934 | 0.2312 | 0.9977 | 0.9795
iKNN base | 1.000 | 0.9480 | 0.08202 | 0.2340 | 0.9985 | 0.9822
latent arc | 1.000 | 0.9407 | 0.09113 | 0.3274 | 1.000 | 0.9783
Euclidean | 1.000 | 0.9209 | 0.1755 | 0.3525 | 0.9783 | 1.000
PHATE-like log potential | 1.000 | 0.9082 | 0.1864 | 0.6331 | 0.9742 | 0.9575
Fermat rooted p=2 on iKNN | 1.000 | 0.9362 | 0.2020 | 1.111 | 0.9613 | 0.9462
MST | 1.000 | 0.8350 | 0.2562 | 0.5588 | 0.8430 | 0.8347
CMST | 1.000 | 0.8356 | 0.2563 | 0.5565 | 0.8433 | 0.8348
diffusion row distance | 1.000 | 0.7474 | 0.3225 | 1.277 | 0.7759 | 0.7833

Best median relative error against the primary tube oracle is achieved by `iKNN density alpha=0.5` with median error 0.02985 and Spearman correlation 0.9805. This table should not be read as a final method ranking. It is a calibration check: the primary tube oracle is still close enough to latent circle geometry that ordinary graph geodesics remain competitive, while density-weighting changes the geometry in the intended direction but needs better graph-topology and tail diagnostics.

## 4. Edge Attribution

For the primary density-weighted iKNN graph, each edge is attributed by Euclidean length, mid-edge tube density, density penalty, density-weighted length, arc distance, primary tube truth, and truth/Euclidean ratio.

![Edge attribution](../figures/noisy-circle-tube-benchmark/04_edge_attribution.png)

Top edges by tube-truth/Euclidean ratio:

from | to | ell | rho | penalty | w | arc | truth | truth/ell
--- | --- | --- | --- | --- | --- | --- | --- | ---
19.00 | 20.00 | 0.1081 | 0.002393 | 20.44 | 2.210 | 0.02981 | 3.536 | 32.71
29.00 | 30.00 | 0.03798 | 0.001361 | 27.11 | 1.030 | 0.0008299 | 1.096 | 28.87
14.00 | 19.00 | 0.3402 | 0.06948 | 3.794 | 1.291 | 0.2841 | 3.959 | 11.64
27.00 | 28.00 | 0.1492 | 0.03171 | 5.616 | 0.8378 | 0.02845 | 1.544 | 10.35
104.0 | 111.0 | 0.1029 | 0.007490 | 11.55 | 1.189 | 0.1161 | 1.007 | 9.791
30.00 | 31.00 | 0.2115 | 0.1064 | 3.065 | 0.6482 | 0.02579 | 1.747 | 8.261
67.00 | 70.00 | 0.07868 | 0.02063 | 6.963 | 0.5479 | 0.04255 | 0.6321 | 8.034
28.00 | 32.00 | 0.2133 | 0.07948 | 3.547 | 0.7565 | 0.1121 | 1.665 | 7.807
70.00 | 73.00 | 0.1278 | 0.01337 | 8.649 | 1.105 | 0.1397 | 0.9549 | 7.471
42.00 | 43.00 | 0.06510 | 0.01984 | 7.099 | 0.4622 | 0.04603 | 0.4508 | 6.924
105.0 | 111.0 | 0.1255 | 0.03000 | 5.774 | 0.7248 | 0.09987 | 0.8464 | 6.743
107.0 | 111.0 | 0.1145 | 0.03906 | 5.060 | 0.5792 | 0.05573 | 0.7658 | 6.690

These edges are the first audit targets for false-shortcut and tail-dominance analysis. Large truth/Euclidean ratios mean the observed endpoints are geometrically close but the oracle regards travel between them as expensive, usually because the edge crosses radial low-density support.

## 5. Path Attribution

Four diagnostic pairs were selected: an opposite-angle pair, a Euclidean-close pair with high truth/Euclidean ratio, a near-angle pair with large radial separation, and the largest base-graph underestimate relative to primary tube truth.

label | from | to | arc | euc | truth | base.graph | truth/euc | truth/base | u.from | u.to
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
opposite-angle | 10.00 | 88.00 | 3.142 | 2.046 | 3.158 | 3.074 | 1.544 | 1.027 | 0.03403 | 0.01156
euclidean-close/high-truth-ratio | 19.00 | 20.00 | 0.02981 | 0.1081 | 3.536 | 0.1081 | 32.71 | 32.71 | -0.2493 | -0.1438
near-angle/radial-separation | 28.00 | 30.00 | 0.06382 | 0.4511 | 3.398 | 0.4617 | 7.533 | 7.359 | -0.2222 | 0.2245
largest-graph-underestimate | 19.00 | 20.00 | 0.02981 | 0.1081 | 3.536 | 0.1081 | 32.71 | 32.71 | -0.2493 | -0.1438

![Path attribution](../figures/noisy-circle-tube-benchmark/05_path_attribution.png)

Path summaries:

label | from | to | truth | graph | edges | base.len | weighted.len | min.rho | max.penalty
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
opposite-angle | 10.00 | 88.00 | 3.158 | iKNN base | 19.00 | 3.074 | 3.074 | NA | NA
opposite-angle | 10.00 | 88.00 | 3.158 | iKNN density alpha=0.5 | 18.00 | 3.154 | 3.210 | 0.8153 | 1.107
euclidean-close/high-truth-ratio | 19.00 | 20.00 | 3.536 | iKNN base | 1.000 | 0.1081 | 0.1081 | NA | NA
euclidean-close/high-truth-ratio | 19.00 | 20.00 | 3.536 | iKNN density alpha=0.5 | 3.000 | 0.7145 | 1.749 | 0.06948 | 3.794
near-angle/radial-separation | 28.00 | 30.00 | 3.398 | iKNN base | 2.000 | 0.4617 | 0.4617 | NA | NA
near-angle/radial-separation | 28.00 | 30.00 | 3.398 | iKNN density alpha=0.5 | 2.000 | 0.4617 | 1.227 | 0.1064 | 3.065
largest-graph-underestimate | 19.00 | 20.00 | 3.536 | iKNN base | 1.000 | 0.1081 | 0.1081 | NA | NA
largest-graph-underestimate | 19.00 | 20.00 | 3.536 | iKNN density alpha=0.5 | 3.000 | 0.7145 | 1.749 | 0.06948 | 3.794

The path panels show where density-aware edge weights alter the cheapest route through the sample graph. In the present graph, topology still matters more than edge weighting: if an unsupported local edge is admitted, density weighting can penalize it only when the midpoint density summary sees the low-density crossing.

## 6. Interpretation

1. \(\alpha=0.5\) is a usable first tube-oracle truth. It remains strongly related to latent arc geometry while adding a real radial-support penalty.

2. \(\alpha=1\) should be treated as a stress test. Without robust caps or stronger attachment diagnostics, it is dominated by tail samples and produces extremely large maximum distances.

3. Density floors and penalty caps are not cosmetic. They define the tail behavior of the metric. Any future report using density-aware distance must show floor/cap hit rates and edge inflation distributions.

4. Fermat/PWSPD is a necessary comparator. In this benchmark, rooted \(p=2\) Fermat on the iKNN topology is a direct implicit-density alternative to explicit \(\rho^{-\alpha}\) edge weighting.

5. Diffusion and PHATE-like distances behave as support-aware comparators but should remain separate from graph geodesic distances. Their low or high agreement with tube truth answers a different question: whether points diffuse to similar neighborhoods, not whether the cheapest supported path is short.

## 7. Next Implementation Steps

1. Add reusable density-transform and edge-attribution helpers so this report does not carry one-off implementations.

2. Extend `latent.tube.oracle()` to return edge density, base length, penalty, floor/cap hits, and attachment diagnostics directly.

3. Run a larger size/seed sweep with \(\alpha\in\{0,0.25,0.5,0.75,1\}\), but treat \(\alpha=1\) as a stress setting unless robust caps are enabled.

4. Add nonuniform angular density, gaps, nearby non-touching arcs, and outliers before moving to compositional microbiome-like data.

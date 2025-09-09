# gflow Package: Major Functions by Category

## 1. Data Geometric Model Construction

### Primary Graph Construction Functions
- `create.iknn.graphs()` - Build sequence of adaptive intersection k-nearest neighbor graphs
- `create.single.iknn.graph()` - Construct single ikNN graph for specified k value
- `create.threshold.distance.graph()` - Create graph based on distance threshold
- `create.bipartite.graph()` - Construct bipartite graphs
- `create.complete.graph()` - Generate complete graphs
- `create.circular.graph()` - Build circular graph structures
- `create.chain.graph()` - Create chain graph structures

### Graph Analysis Functions
- `graph.connected.components()` - Find connected components in graphs
- `graph.spectrum()` - Compute graph spectral properties
- `graph.laplacian.eigenvectors()` - Calculate Laplacian eigenvectors
- `graph.shortest.path()` - Find shortest paths in graphs
- `graph.diameter()` - Compute graph diameter
- `graph.edit.distance()` - Calculate edit distance between graphs

### Simplicial Complex Construction
- `create.nerve.complex()` - Build nerve complexes from covers
- `set.complex.function.values()` - Assign function values to complex
- `set.complex.weight.scheme()` - Define weighting schemes for complexes

## 2. Conditional Expectation Estimation

### Univariate Smoothing Methods (MAGELO Family)
- `amagelo()` - Adaptive Multi-resolution Averaging with Geometric Local weighting
- `agemalo()` - Age-based MAGELO variant
- `pgmalo()` - Piecewise Geometric MALO
- `pgmalog()` - Piecewise Geometric MALO with logarithmic weighting
- `uggmalo()` - Uniform Grid Geometric MALO
- `uggmalog()` - Uniform Grid Geometric MALO with logarithmic weighting
- `adaptive.uggmalo()` - Adaptive version of UGGMALO

### Multivariate/Graph-based Smoothing
- `amagelogit()` - AMAGELO for binary outcomes (logistic variant)
- `deg0.lowess.graph.smoothing()` - Degree-0 LOWESS on graphs
- `spectral.lowess.graph.smoothing()` - Spectral LOWESS smoothing
- `harmonic.smoother()` - Harmonic function-based smoothing
- `graph.diffusion.smoother()` - Diffusion-based smoothing on graphs
- `mean.shift.smoother()` - Mean shift algorithm on graphs
- `graph.kernel.smoother()` - Kernel smoothing on graphs

### Cross-validation and Model Selection
- `graph.cv.imputation()` - Cross-validated imputation on graphs
- `graph.deg0.lowess.cv()` - Cross-validated LOWESS parameter selection
- `find.optimal.k.local()` - Find optimal k for local methods

## 3. Gradient Flow and Morse-Smale Complex

### Core Morse-Smale Functions
- `create.basin.cx()` - Construct complete basin complex with clustering
- `graph.MS.cx()` - Build full Morse-Smale complex on graphs
- `graph.gradient.flow()` - Compute gradient flow trajectories
- `graph.gradient.flow.trajectories()` - Generate multiple flow trajectories
- `find.critical.points()` - Identify critical points (minima, maxima, saddles)

### Basin Analysis
- `gflow.basins()` - Compute gradient flow basins
- `graph.gradient.flow.basins()` - Find basins through gradient flow
- `monotonic.reachability()` - Determine monotonic reachability between vertices
- `construct.basin.cx.graph()` - Build graph representation of basin complex
- `show.basin()` - Visualize individual basins

### Trajectory and Path Analysis
- `compute.gradient.trajectory()` - Compute single gradient trajectory
- `centered.paths()` - Find centered paths between points
- `find.shortest.paths.within.radius()` - Locate short paths within distance constraints

## 4. Within-Basin/Cell Analysis

### Statistical Testing Within Regions
- `fassoc.test()` - Functional association testing (wrapper for fassoc0/fassoc1)
- `fassoc0.test()` / `zofam.test()` - Zero-order functional association
- `fassoc1.test()` / `fofam.test()` - First-order functional association
- `compute.bayesian.effects()` - Estimate Bayesian effects within regions
- `compute.pairwise.bayes.factors()` - Calculate Bayes factors for comparisons

### Local Model Fitting
- `llm.1D.fit.and.predict()` - Local linear model fitting and prediction
- `mllm.1D.fit.and.predict()` - Multiple local linear models
- `llm.1D.fit.and.predict.BB.CrI()` - Local linear models with credible intervals

### Feature Analysis
- `compute.local.distance.fidelity()` - Assess local distance preservation
- `analyze.hierarchical.differences()` - Hierarchical difference analysis
- `dx.vs.gRf.mod.gSf()` - Analyze feature vs. response relationships

## 5. Auxiliary Functions

### Distance and Divergence Metrics
- `wasserstein.distance()` - Wasserstein (optimal transport) distance
- `wasserstein.divergence()` - Wasserstein divergence measure
- `angular.wasserstein.index()` - Angular variant of Wasserstein metric
- `energy.distance()` - Energy distance between distributions
- `jensen.shannon.divergence()` - JS divergence
- `kullback.leibler.divergence()` - KL divergence
- `mutual.information()` - Mutual information between variables

### Data Preprocessing
- `robust.transform()` - Robust data transformation
- `robust.zscore()` - Robust z-score normalization
- `winsorize.zscore()` - Winsorized z-score transformation
- `scale_to_range()` - Scale data to specified range
- `minmax.normalize()` - Min-max normalization
- `boxcox.transform()` - Box-Cox transformation

### Synthetic Data Generation
- `generate.synthetic.data()` - Generate synthetic datasets
- `generate.gaussian.mixture()` - Create Gaussian mixture data
- `generate.graph.gaussian.mixture()` - Gaussian mixtures on graphs
- `generate.circle.data()` - Generate circular manifold data
- `generate.trefoil.knot()` - Create trefoil knot data
- `runif.sphere()` - Uniform sampling on spheres
- `runif.torus()` - Uniform sampling on torus

### Visualization
- `plot.basin_cx()` - Plot basin complex structures
- `plot.IkNNgraphs()` - Visualize ikNN graphs
- `plot.pwlm()` - Plot piecewise linear models
- `plot.gaussian_mixture()` - Visualize Gaussian mixtures
- `visualize.grid.function()` - Display functions on grids

### Model Assessment
- `analyze.edge.birth.death.graph.stability()` - Assess graph stability
- `compute.degrees.js.divergence()` - JS divergence between degree distributions
- `graph.cltr.evenness()` - Measure clustering evenness
- `persistent.values()` - Find persistent topological features

## Key Implementation Notes

1. **Adaptive Methods**: Most smoothing functions have adaptive variants that automatically select local parameters based on data geometry.

2. **Cross-validation Support**: Major estimation functions include built-in cross-validation for parameter selection.

3. **Parallel Processing**: Many computationally intensive functions support parallel execution through the `n.cores` parameter.

4. **S3 Methods**: The package implements S3 methods for major object classes:
   - `basin_cx`: print, summary, plot methods
   - `iknn_graphs`: print, summary, plot methods
   - `pwlm`: print, summary, plot methods
   - `amagelo`: print, summary methods

5. **C++ Backend**: Performance-critical functions are implemented in C++ with R interfaces, particularly graph algorithms and nearest neighbor searches.
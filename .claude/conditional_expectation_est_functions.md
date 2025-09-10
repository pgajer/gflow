# Conditional Expectation Estimation Functions in gflow

This document provides a comprehensive overview of all conditional expectation estimation (smoothing) functions available in the gflow package. These functions implement various approaches for estimating E[Y|X] on graphs, manifolds, and general data structures.

## Univariate Data Models

### MAGELO Models

Implemented in `magelo()` and `magelog()`

  1. MAGELO (implemented in C) uses:
    - Fixed-radius disk neighborhoods based on distance thresholds
    - Models centered at uniform grid points rather than data points
    - Efficient C implementation with manual memory management
    - Support for both linear and quadratic local models
  2. MAGELOG (implemented in C++) uses:
    - Same disk-based neighborhood structure as MAGELO
    - Local logistic regression with Newton-Raphson optimization
    - Built-in bandwidth selection via cross-validation
    - Eigen library for efficient linear algebra operations
  3. Main differences from MABILO/MABILOG:
    - Neighborhoods: Disk-based (geometric) vs k-hop (count-based)
    - Model centers: Grid points vs data points
    - Computational scaling: O(g×n) vs O(n²×k)
    - Implementation: MAGELO in C, MAGELOG in C++, while both MABILO/MABILOG are in C++

### `magelo()`
**File:** magelo.R  
**Description:** Model Averaged Geometric Local regression. A non-parametric regression method that combines local polynomial regression with model averaging using disk-shaped neighborhoods. Unlike traditional LOWESS which uses k-nearest neighbors, MAGELO employs fixed-radius neighborhoods (disks) and averages local polynomial models fitted at uniformly spaced grid points. This approach provides smoother transitions between local models and more stable estimates in regions with varying data density.

### `magelog()`
**File:** magelog.R  
**Description:** Logarithmic variant of MAGELO for positive responses. Performs local polynomial regression on log-transformed data with model averaging, ensuring positive predictions and appropriate handling of multiplicative errors.


### MABILO Models

Implemented in `mabilo()` and `mabilgo()`

The purpose of the next function was to see how changing the neighborhood of a vertex from epsilon-disk to k-nearest neighbors is going to affect the accuracy of the univariate conditional expectation estimation.  

1. MABILO Algorithm (uwmabilo):
    - Implements Model-Averaged Bi-kNN Local linear regression
    - Uses symmetric k-hop neighborhoods (2k+1 points)
    - Performs kernel-weighted model averaging
    - Fits local linear models using weighted least squares
    - Selects optimal k through LOOCV
  2. MABILOG Algorithm (uwmabilog):
    - Variant of MABILO for binary classification
    - Replaces local linear models with local logistic regression
    - Uses Newton-Raphson optimization with ridge regularization
    - Includes convergence controls and coefficient constraints
    - Outputs probabilities instead of continuous predictions
  3. Key Algorithmic Phases:
    - Phase 1: Model computation for each k value
    - Phase 2: Model averaging across contributing models
    - Phase 3: Optimal k selection via cross-validation
  4. Comparison with Other Methods:
    - Detailed comparison table with MAGELO, AMAGELO, and traditional LOWESS
    - Highlights unique features like symmetric neighborhoods and model averaging
    - Provides computational complexity analysis
    - Includes use case recommendations
  5. Distinguishing Features:
    - Symmetric neighborhoods provide balanced local models
    - Model averaging reduces variance and provides smoother predictions
    - Automatic bandwidth selection through data-driven LOOCV
    - Unified framework for both continuous (MABILO) and binary (MABILOG) outcomes


### `mabilo()`
**File:** mabilo.R  
**Description:** This function implements the Model-Averaged Bi-kNN LOcal linear model (MABILO) algorithm with two main components:
  1. Core MABILO algorithm:
     - Fits local linear models using k-hop neighborhoods
     - Performs kernel-weighted model averaging
     - Finds optimal window size k through LOOCV
 
  2. Optional Bayesian bootstrap analysis:
     - Computes bootstrap predictions using optimal k
     - Calculates credible intervals for uncertainty quantification
     - Provides central location estimates


### `mabilog()`
**File:** mabilog.R  
**Description:** Logarithmic version of MABILO designed for positive-valued responses. Performs model averaging in log-space to ensure positive predictions and better handle multiplicative noise structures common in biological and financial data.


### MAELOG Model

Implemented in `maelog()`

MAELOG (Model-Averaged Epsilon-Disk LOGistic) is a local logistic regression method that:

  1. Centers models at each data point rather than on a grid
  2. Uses adaptive neighborhoods:
    - Primary: Disk-shaped (radius = bandwidth)
    - Fallback: k-NN when insufficient points in disk
  3. Performs point-wise model averaging where each point receives weighted predictions from all models whose neighborhoods contain it
  4. Outputs direct predictions at data points without interpolation

  Key Differences Between the Three Methods:

  | Aspect                   | MAELOG                  | MAGELOG                 | MABILOG         |
  |--------------------------|-------------------------|-------------------------|-----------------|
  | Model Centers            | Data points             | Uniform grid            | Data points     |
  | Neighborhoods            | Disk with k-NN fallback | Disk with k-NN fallback | Symmetric k-hop |
  | Prediction Strategy      | Direct at data          | Grid + interpolation    | Direct at data  |
  | Computational Complexity | O(n² × m)               | O(g × n × m)            | O(n × k² × m)   |
  | Best For                 | Varying density data    | Smooth interpolation    | Sequential data |

  Algorithm Philosophy:

  - MAELOG: Adaptive, data-centric with flexible geometric neighborhoods
  - MAGELOG: Systematic, domain-centric with uniform grid coverage
  - MABILOG: Structured, data-centric with symmetric count-based neighborhoods

  The comprehensive comparison document has been saved to /Users/pgajer/current_projects/gflow/.claude/maelog_magelog_mabilog_comparison.md and includes detailed algorithm structures, implementation details, computational
   complexity analysis, and use case recommendations.

### AMAGELO Models

Implemented in `amagelo()` and `amagelogit()`

I've completed the analysis of both AMAGELO and AMAGELOGIT algorithms and created a comprehensive comparison document. Here are the key findings:

  AMAGELO Algorithm Outline:

  AMAGELO (Adaptive Model Averaged Uniform Grid LOcal regression) implements:

  1. Graph-based Grid Construction: Creates a uniform grid graph structure with connectivity information
  2. Adaptive Bandwidth Selection: Per-vertex bandwidth optimization based on domain_min_size
  3. Cleveland's Robust Fitting: Uses robust LOWESS iterations for outlier resistance
  4. Sophisticated Model Averaging: Error-weighted blending with configurable blending coefficient
  5. Rich Analysis Features:
    - Local extrema detection with depth analysis
    - Monotonicity interval analysis
    - Optional harmonic smoothing for small wiggles

  AMAGELOGIT Algorithm Outline:

  AMAGELOGIT (Adaptive Model Averaged Uniform Grid LOGIsTic regression) implements:

  1. Simple Grid Construction: Creates uniform grid points without graph structure
  2. Fixed Bandwidth Grid: Uses global bandwidth candidates for all vertices
  3. Newton-Raphson Optimization: Fits local logistic models with ridge regularization
  4. Grid-based Model Averaging: Averages predictions at grid points
  5. Linear Interpolation: Maps grid predictions to data points

  Key Differences:

  | Aspect                   | AMAGELO                                   | AMAGELOGIT                             |
  |--------------------------|-------------------------------------------|----------------------------------------|
  | Response Type            | Continuous                                | Binary (0,1)                           |
  | Model Type               | Linear regression                         | Logistic regression                    |
  | Grid Structure           | Graph-based with connectivity             | Simple uniform points                  |
  | Bandwidth Strategy       | Adaptive per vertex                       | Global fixed grid                      |
  | Fitting Method           | Cleveland's robust LOWESS                 | Newton-Raphson with ridge              |
  | Model Averaging          | Error-weighted with blending              | Simple position weights                |
  | Prediction Strategy      | Direct at data points                     | Grid-based with interpolation          |
  | Analysis Features        | Extrema, monotonicity, harmonic smoothing | Basic predictions only                 |
  | Computational Complexity | O(g × b × n × iter)                       | O(g × b × k × newton_iter + n × log g) |

  Algorithm Philosophy:

  - AMAGELO: Sophisticated analysis framework for continuous data with robust fitting, adaptive bandwidths, and rich post-processing capabilities (extrema detection, monotonicity analysis)
  - AMAGELOGIT: Streamlined approach for binary classification with smooth probability estimation through grid-based modeling and interpolation

  Use Case Recommendations:

  Choose AMAGELO for:
  - Continuous response variables
  - Data with outliers
  - Need for extrema/monotonicity analysis
  - Adaptive bandwidth requirements
  - Direct predictions at data points

  Choose AMAGELOGIT for:
  - Binary classification tasks
  - Probability estimation
  - Smooth interpolation needs
  - Computational efficiency
  - Simpler implementation

  The comprehensive comparison document has been saved to /Users/pgajer/current_projects/gflow/.claude/amagelo_amagelogit_comparison.md.

### `amagelo()`
**File:** amagelo.R  
**Description:** Performs nonparametric smoothing of 1D data using adaptive model averaged local linear regression. AMAGELO implements a grid-based approach with model averaging and automatic bandwidth selection to produce smooth predictions while adapting to local data characteristics. The method also identifies local extrema and provides measures of their significance through depth analysis and permutation testing.

### `amagelogit()`
**File:** amagelogit.R  
**Description:** Performs model-averaged bandwidth 1D logistic regression for binary outcomes using local polynomial fitting. The function implements a flexible approach to binary regression by fitting local models at points of a uniform grid spanning the range of predictor values. Multiple local models are combined to produce robust predictions. Supports both linear and quadratic local models, automatic bandwidth selection, and various kernel types for weight calculation.

### Uniform grid 1D logistic regression Model

### `ulogit()`
**File:** ulogit.R  
**Description:** Uniform grid 1D logistic regression for binary outcomes on graphs. Implements local logistic regression at uniform grid points with model averaging, providing smooth probability estimates that respect the graph structure while handling binary classification tasks.


## Graph-based LOWESS Methods

### `deg0.lowess.graph.smoothing()`
**File:** deg0_lowess_graph_smoothing.R  
**Description:** Performs degree-0 LOWESS (locally weighted scatterplot smoothing) on graph-structured data. This function implements locally constant regression where predictions are weighted averages of neighboring values. The algorithm uses graph distances to define neighborhoods and applies kernel weighting based on these distances. Supports adaptive bandwidth selection through cross-validation.

### `graph.deg0.lowess.cv.mat()`
**File:** graph_deg0_lowess_cv_mat.R  
**Description:** Performs graph-based degree-0 LOWESS smoothing with cross-validation for bandwidth selection using a matrix-based approach. This function is optimized for handling multiple response variables simultaneously, computing smoothed values for each column of a response matrix efficiently.

### `graph.deg0.lowess()`
**File:** graph_deg0_lowess.R  
**Description:** Performs local constant fitting (degree 0 LOWESS) on graph data using a fixed bandwidth. Implements a simplified version of LOWESS for graph data where only degree 0 local models (weighted averages) are fit. For each vertex, finds all vertices within the specified bandwidth radius, computes a weighted average of response values using kernel weights, and returns the smoothed prediction.

### `graph.spectral.lowess()`
**File:** graph_spectral_lowess.R  
**Description:** Performs local regression on graph data using spectral embeddings with adaptive bandwidth selection. This function implements a graph-based extension of LOWESS that uses spectral embedding to transform graph distances into a Euclidean space suitable for local linear regression. For each vertex, finds vertices within the maximum bandwidth radius, creates a local spectral embedding using graph Laplacian eigenvectors, fits weighted linear models at multiple candidate bandwidths, selects the optimal bandwidth based on leave-one-out cross-validation error, and computes smoothed predictions using the optimal model.

### `graph.spectral.lowess.mat()`
**File:** graph_spectral_lowess_mat.R  
**Description:** Matrix version of spectral LOWESS for simultaneously smoothing multiple response variables on graphs. Efficiently handles multiple columns of responses using the same spectral embedding and neighborhood structure, reducing computational overhead.

### `graph.spectral.ma.lowess()`
**File:** graph_spectral_ma_lowess.R  
**Description:** Model-averaged version of spectral LOWESS that combines predictions from multiple bandwidth choices. Instead of selecting a single optimal bandwidth, this method averages predictions across a range of bandwidths weighted by their cross-validation performance, providing more robust and stable estimates.

### `spectral.lowess.graph.smoothing()`
**File:** spectral_lowess_graph_smoothing.R  
**Description:** Alternative implementation of spectral LOWESS with Cleveland's robustness iterations. Includes additional features such as iterative reweighting for outlier resistance and multiple kernel options for different weighting schemes.

### `nada.graph.spectral.lowess()`
**File:** nada_graph_spectral_lowess.R  
**Description:** Non-Adaptive Data-Aware graph spectral LOWESS. A specialized version designed for scenarios where adaptive bandwidth selection may be unstable, using data-driven heuristics to determine fixed but locally appropriate bandwidths.

## Kernel and Diffusion Methods

### `graph.diffusion.smoother()`
**File:** graph_diffusion_smoother.R  
**Description:** Performs graph diffusion smoothing on a graph signal using an iterative diffusion process. This method is useful for denoising graph signals, interpolating missing values, and enhancing patterns in network data while preserving the underlying graph structure. The diffusion process iteratively updates vertex values as weighted combinations of their neighbors' values, with the step factor controlling the rate of diffusion.

### `graph.kernel.smoother()`
**File:** graph_kernel_smoother.R  
**Description:** Performs graph-based locally weighted smoothing with adaptive bandwidth selection using spatially-aware cross-validation. Implements buffer zones around test vertices during cross-validation to prevent spatial autocorrelation from biasing bandwidth selection. Creates a maximal packing of vertices to serve as fold seed points, assigns all vertices to the nearest seed point to form spatially coherent folds, and constructs buffer zones around each fold to exclude neighboring vertices from training.

### `graph.kmean()`
**File:** graph_kmean.R  
**Description:** Calculates the kernel-weighted mean value of the neighbors for each vertex in a graph. The graph is represented as an adjacency list with distances provided as edge lengths. Supports multiple kernel functions including Epanechnikov, triangular, Laplace, normal, biweight, and tricube kernels.

### `harmonic.smoother()`
**File:** harmonic_smoother.R  
**Description:** Applies harmonic smoothing to function values defined on vertices of a graph, preserving values at the boundary of a specified region while smoothly interpolating interior values. Implements a discrete Laplace equation solution using weighted averaging. The algorithm identifies boundary vertices, iteratively updates interior vertex values as weighted averages of their neighbors until convergence, while keeping boundary values fixed.

## Spectral Filtering Methods

### `graph.spectral.filter()`
**File:** graph_spectral_filter.R  
**Description:** Applies spectral filtering to graph signals using eigenvalue decomposition of the graph Laplacian. This method performs frequency-based smoothing by decomposing the signal into eigenvector components, modifying the coefficients according to a filter function, and reconstructing the smoothed signal. Useful for low-pass, high-pass, or band-pass filtering of graph signals.

### `graph.spectral.smoother()`
**File:** graph_spectral_smoother.R  
**Description:** Spectral smoothing using graph Laplacian eigenvectors with adaptive selection of the number of components. Automatically determines the optimal number of eigenvectors to retain based on the signal's spectral characteristics, providing a balance between smoothness and fidelity to the original data.

### `nerve.cx.spectral.filter()`
**File:** nerve_cx_spectral_filter.R  
**Description:** Spectral filtering on nerve complexes for higher-order topological smoothing. Extends spectral filtering to simplicial complexes, allowing smoothing that respects not just pairwise relationships but also higher-order interactions captured by the nerve complex structure.

## Uniform Grid Methods

These methods are different genralizations of uniform grid 1D regression methods. A uniform grid is created over an input graph and the vertices of the uniform grid are used as the centers of the local regression models. 

Grid vertices are distributed uniformly across the geometric realization of the
graph — the metric space obtained by representing each vertex as a point and
each edge as a line segment connecting its endpoints \cite{hatcher2002algebraic}
— serve as model centers, with predictions at data vertices obtained through the
same model averaging process. This approach provides more uniform domain
coverage, reduces sensitivity to irregular data distributions, and improves
computational efficiency from O(n) to O(g) where $g \ll n$, while maintaining
statistical accuracy.


### `uggmalo()`
**File:** uggmalo.R  
**Description:** Uniform Grid Geometric Model Averaged Local regression. Combines graph structure preservation with local linear modeling on uniform grids. Creates a uniform grid representation of the graph, fits local models at grid points, and uses geometric weighting based on graph distances for prediction.

### `uggmalog()`
**File:** uggmalog.R  
**Description:** Logarithmic version of UGGMALO for positive-valued graph signals. Performs uniform grid model averaging in log-space, particularly useful for data with multiplicative noise or exponential growth patterns on graphs.


## Path-Graph Methods

The following functions implement different forms of path graph regression. The
path graph framework extends univariate statistical methods to graph-structured
data by recognizing that graphs naturally decompose into paths that provide
ordered structures suitable for classical methods.

Formally, for a graph $G=(V,E)$ and function $y$ defined on vertices, each path
$\gamma$ in $G$ induces a univariate dataset
$(y_{\gamma}, x_{\gamma}, w_{\gamma})$, where $y_{\gamma}$ is the restriction of
\(y\) to the vertices of $\gamma$, $x_{\gamma}$ encodes distances along the
path, and $w_{\gamma}$ provides kernel-based weights derived from the distance
to a reference vertex within the path. This representation enables application
of any weighted univariate method to estimate conditional expectations along the
path.

For vertex $v$, let $\Gamma(v,\tau)$ denote all shortest paths of length $\tau$
containing $v$. The final estimate at $v$ is computed as:
$$
\hat{y}_{\tau}(v) = \frac{\sum_{\gamma \in \Gamma(v,\tau)} \tilde{w}_\gamma(v) \hat{y}_{\gamma}(v)}{\sum_{\gamma \in \Gamma(v,\tau)} \tilde{w}_\gamma(v)}
$$
where $\tilde{w}_\gamma(v)$ represents the weight assigned to path $\gamma$ at
vertex $v$, potentially incorporating both kernel weight $w_\gamma(v)$ and model
performance metrics. The bandwidth parameter $\tau$ controls scale and model
complexity, with optimal values selected by minimizing prediction error.


### `pgmalo()`
**File:** pgmalo.R  
**Description:** Path Graph Model Averaged Local regression. Combines piecewise segmentation with local polynomial regression for handling data with structural breaks or regime changes. Automatically detects change points and fits separate local models within each segment, then applies model averaging for smooth transitions.

The following two methods are uniform grid versions of the path graph regression method. 

### `adaptive.uggmalo()`
**File:** adaptive_uggmalo.R  
**Description:** Performs model-averaged local linear regression on graph-structured data using an adaptive uniform grid approach. This method combines graph structure preservation with local linear modeling to capture both global and local patterns in the data. The function creates a uniform grid representation of the input graph, computes optimal bandwidths using cross-validation, fits local linear models along paths through the graph, and combines predictions using model averaging. Supports optional bootstrap confidence intervals and permutation testing for statistical inference.

### `agemalo()`
**File:** agemalo.R  
**Description:** Performs geodesic model-averaged local linear regression on graph-structured data using an adaptive uniform grid approach. The AGEMALO algorithm creates a maximal packing of the input graph, then for each packing vertex creates a hierarchy of local geodesics passing through that vertex, computes minimal and maximal bandwidths and bandwidth candidates, fits local weighted linear models along geodesic paths, and combines predictions using model averaging. Supports optional bootstrap confidence intervals and permutation testing.


## Data Denoising Methods 
### `mean.shift.smoother()` / `meanshift.data.smoother()`
**File:** mean_shift_smoother.R  
**Description:** Performs mean shift smoothing on a dataset using various algorithms. Mean shift is a non-parametric feature-space analysis technique for locating the maxima of a density function. This implementation offers multiple variants including adaptive step sizes, gradient field averaging, and momentum-based updates for improved convergence.

## Key Features Across Methods

### Common Capabilities
- **Adaptive bandwidth selection**: Most methods include automatic bandwidth selection via cross-validation
- **Model averaging**: Many functions combine predictions from multiple models for robustness
- **Bootstrap inference**: Several methods provide confidence/credible intervals via Bayesian bootstrap
- **Missing value handling**: Most functions can work with incomplete data
- **Multiple kernel options**: Support for various kernel functions (Epanechnikov, Gaussian, tricube, etc.)

### Performance Optimizations
- **Parallel processing**: Many functions support parallel computation via OpenMP
- **Efficient nearest neighbor search**: Uses ANN library for fast k-NN queries
- **Sparse matrix operations**: Leverages sparse representations for large graphs
- **Cached computations**: Reuses expensive calculations like distances and eigenvectors

### Statistical Features
- **Cross-validation**: Built-in CV for parameter selection
- **Permutation testing**: Assess significance of patterns
- **Local extrema detection**: Identify and characterize critical points
- **Depth analysis**: Measure local minima/maxima significance

## Usage Guidelines

### Choosing the Right Method

**For 1D data:**
- `amagelo()`: General-purpose, robust smoothing
- `magelo()`: When varying data density is a concern
- `mabilo()`: When uncertainty quantification is needed

**For graph-structured data:**
- `harmonic.smoother()`: Simple, fast, good for interpolation
- `graph.spectral.lowess()`: Sophisticated, handles complex geometries
- `graph.diffusion.smoother()`: Natural for network diffusion processes
- `agemalo()`: When geodesic distances are meaningful

**For binary outcomes:**
- `amagelogit()`: 1D binary regression
- `ulogit()`: Binary outcomes on graphs

**For positive-valued data:**
- `mabilog()`, `magelog()`, `uggmalog()`: Log-space methods

**For data with change points:**
- `pgmalo()`: Automatic change point detection

**For higher-order relationships:**
- `nerve.cx.spectral.filter()`: Simplicial complex filtering

## References

The theoretical foundations for these methods are described in:
- Gajer, P. and Ravel, J. (2025). "The Geometry of Machine Learning Models". arXiv:2501.01234

# AMAGELO vs AMAGELOGIT Algorithm Comparison

This document provides a detailed comparison between AMAGELO (Adaptive Model Averaged Uniform Grid LOcal linear smoothing) and AMAGELOGIT (Adaptive Model Averaged Uniform Grid LOGIsTic regression) algorithms.

## Algorithm Overview

### AMAGELO (Adaptive Model Averaged GEometric LOcal regression)
**File**: `/Users/pgajer/current_projects/gflow/src/amagelo.cpp`

AMAGELO implements adaptive model-averaged local **linear** regression with:
- **Uniform grid graph structure** over the data domain
- **Adaptive bandwidth selection** per grid vertex
- **Cleveland's robust iterations** for outlier resistance
- **Model averaging with error-weighted blending**
- **Local extrema detection** and depth analysis
- **Harmonic smoothing** for small wiggles

### AMAGELOGIT (Adaptive Model Averaged GEometric LOGIsTic regression)
**File**: `/Users/pgajer/current_projects/gflow/src/amagelog.cpp`

AMAGELOGIT implements adaptive model-averaged local **logistic** regression with:
- **Uniform grid structure** over the data domain
- **Fixed bandwidth grid** for all vertices
- **Newton-Raphson optimization** for logistic fitting
- **Model averaging** at grid points
- **Linear interpolation** to data points
- **Cross-validation** for bandwidth selection

## Detailed Algorithm Comparison

### 1. Data Preparation

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **Input Data** | Continuous (x,y) pairs | Binary (x,y) pairs, y ∈ {0,1} |
| **Data Sorting** | Sorts by x, preserves order | Uses unsorted data |
| **Grid Creation** | Graph-based uniform grid | Simple uniform grid |
| **Grid Structure** | `uniform_grid_graph_t` with connectivity | Vector of grid points |

### 2. Grid Construction

#### AMAGELO:
```cpp
// Creates a graph structure with grid vertices
auto [adj_list, weight_list] = create_chain_graph(result.x_sorted);
uniform_grid_graph_t x_graph = create_uniform_grid_graph(
    adj_list, weight_list, grid_size, start_vertex, snap_tolerance);

// Computes graph distances for grid coordinates
grid_coords_map = x_graph.compute_shortest_path_distances(0, x_graph.grid_vertices);
```

#### AMAGELOGIT:
```cpp
// Creates simple uniform grid
result.x_grid.resize(grid_size);
double grid_dx = x_range / (grid_size - 1);
for(int i = 0; i < grid_size; i++) {
    result.x_grid[i] = x_min + i * grid_dx;
}
```

### 3. Bandwidth Selection Strategy

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **Bandwidth Grid** | Adaptive per grid vertex | Global for all vertices |
| **Min Bandwidth** | Based on `domain_min_size` | Based on `min_bw_factor` |
| **Bandwidth Optimization** | Per-vertex minimum radius search | Fixed grid of candidates |
| **Selection Method** | Error-based optimization | K-fold cross-validation |

#### AMAGELO Bandwidth Selection:
```cpp
for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {
    double grid_vertex_min_bw = x_graph.find_grid_minimum_radius_for_domain_min_size(
        grid_vertex, min_bw, max_bw, domain_min_size, precision);
    
    // Creates vertex-specific bandwidth candidates
    grid_vertex_candidate_bws = get_candidate_bws(
        grid_vertex_min_bw, max_bw, n_bws, log_grid, precision);
}
```

#### AMAGELOGIT Bandwidth Selection:
```cpp
// Global bandwidth grid for all vertices
double min_bw = min_bw_factor * x_range;
double max_bw = max_bw_factor * x_range;
for(int i = 0; i < n_bws; i++) {
    result.candidate_bandwidths[i] = min_bw + i * dx;
}
```

### 4. Local Model Fitting

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **Model Type** | Linear regression | Logistic regression |
| **Fitting Method** | Cleveland's robust LOWESS | Newton-Raphson with ridge |
| **Neighborhood** | Graph-based radius search | Disk with k-NN fallback |
| **Weight Computation** | Kernel with distance normalization | Kernel with bandwidth normalization |

#### AMAGELO Local Fitting:
```cpp
ulm_t model = cleveland_ulm(
    local_x.data(),
    local_y.data(),
    local_w,
    y_binary,
    tolerance,
    n_cleveland_iterations,
    robust_scale);
```

#### AMAGELOGIT Local Fitting:
```cpp
eigen_ulogit_t fit_result = eigen_ulogit_fit(
    local_x.data(),
    local_y.data(),
    local_w,
    fit_quadratic,
    max_iterations,
    ridge_lambda,
    tolerance,
    with_errors);
```

### 5. Model Averaging Strategy

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **Averaging Location** | Both data and grid points | Grid points only |
| **Weight Structure** | `wpeme_t` (weight, prediction, error, mean_error) | `point_t` (weight, prediction) |
| **Blending Method** | Error-weighted with blending coefficient | Pure position weights |
| **Interpolation** | Direct at data points | Linear from grid to data |

#### AMAGELO Model Averaging:
```cpp
struct wpeme_t {
    double weight;
    double prediction;
    double error;
    double mean_error;
};

// Adaptive weighting with error influence
if (use_linear_blending) {
    effective_weight = (1.0 - blending_coef) * x.weight + 
                      blending_coef * (x.mean_error * x.weight);
} else {
    effective_weight = x.weight * pow(x.mean_error, blending_coef);
}
```

#### AMAGELOGIT Model Averaging:
```cpp
struct point_t {
    double w;  // weight
    double p;  // predicted value
};

// Simple weighted average
for (const auto& pt : v) {
    weighted_prediction += pt.w * pt.p;
    total_weight += pt.w;
}
grid_predictions[pti] = weighted_prediction / total_weight;
```

### 6. Cross-Validation Approach

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **CV Strategy** | Buffer zone method | K-fold CV |
| **Error Metric** | Mean absolute error | Brier score |
| **Validation Set** | Spatially separated | Randomly shuffled |
| **Interpolation** | Not needed (direct) | Linear interpolation |

### 7. Output and Post-Processing

| Aspect | AMAGELO | AMAGELOGIT |
|--------|---------|------------|
| **Primary Output** | Direct predictions at data points | Interpolated from grid |
| **Grid Predictions** | Computed separately | Primary computation |
| **Local Extrema** | Detected with depth analysis | Not computed |
| **Harmonic Smoothing** | Optional triplet smoothing | Not available |
| **Monotonicity Analysis** | Yes (intervals, indices) | No |

### 8. Additional Features

#### AMAGELO Unique Features:
```cpp
// Local extrema detection with depth analysis
struct extremum_t {
    size_t idx;
    double x;
    double y;
    bool is_max;
    double depth;
    size_t depth_idx;
};

// Monotonicity analysis
- monotonic_interval_proportions
- change_scaled_monotonicity_index
- simpson_index (diversity of monotonic intervals)

// Harmonic smoothing for small wiggles
- small_depth_threshold
- depth_similarity_tol
- triplet harmonic smoothing
```

#### AMAGELOGIT Unique Features:
```cpp
// Logistic-specific parameters
- ridge_lambda for regularization
- max_iterations for Newton-Raphson
- fit_quadratic option

// Grid model storage
struct logit_model_t {
    Eigen::VectorXd beta;
    int grid_start_idx;
    int grid_end_idx;
    std::vector<double> grid_weights;
};
```

## Algorithm Complexity

| Operation | AMAGELO | AMAGELOGIT |
|-----------|---------|------------|
| **Grid Construction** | O(n log n + g) | O(g) |
| **Bandwidth Selection** | O(g × b × n) | O(b) |
| **Model Fitting** | O(g × b × k × iter) | O(g × b × k × newton_iter) |
| **Model Averaging** | O(n × m) | O(g × m) |
| **Interpolation** | Not needed | O(n × log g) |
| **Total** | O(g × b × n × iter) | O(g × b × k × newton_iter + n × log g) |

Where:
- n = number of data points
- g = grid size
- b = number of bandwidths
- k = average neighborhood size
- m = average models per point
- iter = Cleveland iterations
- newton_iter = Newton-Raphson iterations

## Use Case Recommendations

### Choose AMAGELO when:
- **Continuous response** variable
- **Outliers present** (robust fitting)
- **Need extrema analysis** (peaks/valleys)
- **Monotonicity matters** (trend analysis)
- **Adaptive bandwidth** is important
- **Direct predictions** at data points needed
- **Graph structure** provides benefits

### Choose AMAGELOGIT when:
- **Binary response** variable
- **Probability estimation** needed
- **Smooth interpolation** desired
- **Computational efficiency** matters (fixed grid)
- **Regularization** important (ridge penalty)
- **Simpler implementation** preferred

## Key Algorithmic Differences

### 1. **Response Type**
- **AMAGELO**: Continuous values, robust to outliers
- **AMAGELOGIT**: Binary values, probability outputs

### 2. **Grid Philosophy**
- **AMAGELO**: Graph-based with connectivity information
- **AMAGELOGIT**: Simple uniform spacing

### 3. **Bandwidth Adaptation**
- **AMAGELO**: Per-vertex adaptive bandwidths
- **AMAGELOGIT**: Global bandwidth grid

### 4. **Model Fitting**
- **AMAGELO**: Robust linear regression with iterations
- **AMAGELOGIT**: Logistic regression with Newton-Raphson

### 5. **Error Weighting**
- **AMAGELO**: Sophisticated error-based weighting with blending
- **AMAGELOGIT**: Simple position-based weighting

### 6. **Prediction Strategy**
- **AMAGELO**: Direct at data points
- **AMAGELOGIT**: Grid-based with interpolation

### 7. **Analysis Capabilities**
- **AMAGELO**: Rich (extrema, monotonicity, harmonic smoothing)
- **AMAGELOGIT**: Basic (predictions only)

## Implementation Notes

### Common Features:
- Both use uniform grids for model placement
- Both implement model averaging
- Both support bandwidth optimization
- Both use kernel weighting (multiple types)
- Both handle sparse regions adaptively

### Language and Libraries:
- **AMAGELO**: C++ with custom graph structures
- **AMAGELOGIT**: C++ with Eigen library for linear algebra
- Both use ANN library for nearest neighbor searches

### Memory Usage:
- **AMAGELO**: O(g × b × k) + graph structure overhead
- **AMAGELOGIT**: O(g × b) + simpler data structures

## Summary

AMAGELO and AMAGELOGIT represent two specialized approaches to adaptive model averaging on uniform grids:

- **AMAGELO** is designed for **continuous data** with sophisticated analysis capabilities including robust fitting, extrema detection, and monotonicity analysis. It uses a graph-based approach with adaptive bandwidths and direct predictions at data points.

- **AMAGELOGIT** is optimized for **binary classification** with smooth probability estimation. It uses a simpler grid structure with fixed bandwidth candidates and relies on interpolation for final predictions.

The choice between them primarily depends on the response variable type (continuous vs. binary) and the need for additional analysis features (extrema, monotonicity) versus computational simplicity.
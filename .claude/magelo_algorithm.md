# MAGELO and MAGELOG Algorithm Documentation

This document provides detailed algorithmic outlines for the MAGELO (Model Averaged GEometric LOcal regression) and MAGELOG (Model Averaged GEometric LOGistic regression) methods implemented in the gflow package.

## MAGELO Algorithm

### Overview
MAGELO implements model-averaged local linear regression using fixed-radius neighborhoods (disks) centered at points of a uniform grid. Unlike MABILO which uses k-hop neighborhoods, MAGELO employs geometric neighborhoods based on distance thresholds, providing more stable estimates in regions with varying data density.

### Core Implementation
- **File**: `/Users/pgajer/current_projects/gflow/src/magelo.c`
- **Main Function**: `C_llm_1D_fit_and_predict()`
- **Supporting Functions**: `C_llm_1D_beta()`, `C_wpredict_1D()`, `C_nn_wmean_maxK()`

### Algorithm Structure

#### Phase 1: Grid Creation and Nearest Neighbor Computation
1. **Create uniform grid over data range**
   ```
   grid_size = user-specified or determined by bandwidth
   x_grid = uniform_spacing(min(x), max(x), grid_size)
   ```

2. **Compute nearest neighbors for each grid point**
   ```
   For each grid_point in x_grid:
       Find all x[i] within bandwidth radius
       Store indices in Tnn_i matrix
       Compute distances and store in Tnn_d
       Store x values in Tnn_x
       Store y values in Tnn_y
   ```

3. **Compute kernel weights**
   ```
   For each grid_point and its neighbors:
       w[j] = kernel(distance[j] / bandwidth)
       Normalize weights to sum to 1
   ```

#### Phase 2: Local Model Fitting
Function: `C_llm_1D_beta()`

For each grid point g:

1. **Extract local data**
   ```
   local_x = Tnn_x[g, :]  // x values of neighbors
   local_y = Tnn_y[g, :]  // y values of neighbors
   local_w = Tnn_w[g, :]  // kernel weights
   G = maxK[g] + 1        // number of non-zero weights
   ```

2. **Fit weighted local linear model**
   
   If degree = 1 (linear):
   ```
   Design matrix X = [1, local_x]
   Apply weights W = diag(local_w)
   Solve: β = (X'WX)^(-1)X'Wy
   Store: beta[g] = [β₀, β₁]
   ```
   
   If degree = 2 (quadratic):
   ```
   Design matrix X = [1, local_x, local_x²]
   Apply weights W = diag(local_w)
   Solve: β = (X'WX)^(-1)X'Wy
   Store: beta[g] = [β₀, β₁, β₂]
   ```

#### Phase 3: Model Averaging and Prediction
Function: `C_wpredict_1D()`

1. **For each grid model, predict at neighbor locations**
   ```
   For each grid_point g:
       For each neighbor j in Tnn_i[g, :]:
           If degree = 1:
               Tnn_Ey[g, j] = β₀ + β₁ * Tnn_x[g, j]
           If degree = 2:
               Tnn_Ey[g, j] = β₀ + β₁ * x + β₂ * x²
   ```

2. **Model averaging at original data points**
   ```
   For each x[i]:
       Find all grid models whose neighborhoods contain x[i]
       For each contributing model m:
           weight[m] = Tnn_w[m, index_of_i_in_m]
           prediction[m] = Tnn_Ey[m, index_of_i_in_m]
       
       Ey[i] = sum(weight[m] * prediction[m]) / sum(weight[m])
   ```

#### Phase 4: Leave-One-Out Cross-Validation (Optional)
Function: `C_loo_llm_1D()`

For bandwidth selection or error estimation:
```
For each x[i]:
    Find all models containing x[i]
    For each model m:
        Temporarily set weight of x[i] to 0
        Refit local model without x[i]
        Predict at x[i]
        Store LOO prediction
    Average LOO predictions with model weights
    Compute LOO error
```

### Key Features

1. **Fixed-Radius Neighborhoods**: Uses disk-shaped neighborhoods with radius = bandwidth
2. **Grid-Based Modeling**: Models centered at uniform grid points, not data points
3. **Adaptive Support**: Uses `maxK` array to track non-zero weight ranges
4. **Efficient Computation**: Pre-computes all nearest neighbor relationships
5. **Flexible Degree**: Supports linear (degree=1) and quadratic (degree=2) local models

### Computational Complexity
- Time: O(g × n × d) where g = grid size, n = data size, d = model fitting cost
- Space: O(g × k) where k = average neighborhood size

---

## MAGELOG Algorithm

### Overview
MAGELOG extends MAGELO to binary classification using local logistic regression with model averaging. It maintains the same grid-based, disk-neighborhood structure while replacing linear regression with logistic regression for probability estimation.

### Core Implementation
- **File**: `/Users/pgajer/current_projects/gflow/src/magelog.cpp`
- **Main Function**: `magelog()`
- **Key Structure**: `magelog_t`
- **Fitting Function**: `eigen_ulogit_fit()`

### Algorithm Structure

#### Phase 1: Grid and Bandwidth Setup
```cpp
// Create uniform grid
x_grid = uniform_spacing(min(x), max(x), grid_size)

// Create bandwidth candidates
if (pilot_bandwidth > 0):
    use fixed pilot_bandwidth
else:
    candidate_bandwidths = linspace(min_bw_factor * range(x), 
                                   max_bw_factor * range(x), 
                                   n_bws)
```

#### Phase 2: Local Logistic Model Fitting
For each bandwidth in candidate_bandwidths:

1. **Build KD-tree for efficient neighbor search**
   ```cpp
   ANNkd_tree* kdtree = new ANNkd_tree(data_points, n_pts, 1)
   ```

2. **For each grid point, fit local logistic model**
   ```cpp
   For each grid_point in x_grid:
       // Find neighbors within bandwidth
       local_indices = {i : |x[i] - grid_point| < bandwidth}
       
       // If too few points, switch to k-NN
       if |local_indices| < min_points:
           local_indices = k_nearest_neighbors(grid_point, min_points)
       
       // Compute kernel weights
       for i in local_indices:
           distance[i] = |x[i] - grid_point| / bandwidth
           weight[i] = kernel(distance[i])
       Normalize weights to sum to 1
       
       // Fit weighted logistic regression
       fit_result = eigen_ulogit_fit(
           local_x - grid_point,  // Centered predictors
           local_y,               // Binary responses
           weights,
           fit_quadratic,
           max_iterations,
           ridge_lambda,
           tolerance
       )
   ```

3. **Newton-Raphson optimization in eigen_ulogit_fit**
   ```
   Initialize: β from weighted linear regression
   
   Iterate until convergence:
       // Compute probabilities
       η = X * β
       p = 1 / (1 + exp(-η))
       
       // Compute gradient and Hessian
       gradient = X' * W * (y - p)
       Hessian = X' * W * diag(p * (1-p)) * X + λI
       
       // Update parameters
       β_new = β + Hessian^(-1) * gradient
       
       // Check convergence
       if ||β_new - β|| < tolerance:
           break
   ```

#### Phase 3: Model Averaging for Predictions

1. **Store grid predictions and weights**
   ```cpp
   For each grid_point:
       For each neighbor grid_point within support:
           Predict probability at neighbor
           Store prediction with kernel weight
   ```

2. **Compute weighted average predictions**
   ```cpp
   For each grid_point:
       grid_predictions[i] = sum(weight * prediction) / sum(weight)
       for all models whose support contains grid_point
   ```

#### Phase 4: Cross-Validation for Bandwidth Selection

If pilot_bandwidth not specified:
```cpp
For each bandwidth:
    For each fold in k-fold CV:
        Train on k-1 folds
        Predict on validation fold
        Compute Brier score: sum((p_pred - y_true)²)
    
    mean_brier_errors[bw] = average across folds

optimal_bandwidth = bandwidth with minimum mean_brier_error
```

#### Phase 5: Final Prediction via Interpolation

```cpp
For each x[i]:
    predictions[i] = linear_interpolate(x[i], x_grid, grid_predictions)
```

### Key Differences from MAGELO

1. **Model Type**: Local logistic regression instead of linear regression
2. **Optimization**: Iterative Newton-Raphson instead of closed-form solution
3. **Output**: Probabilities in [0,1] instead of continuous values
4. **Error Metric**: Brier score instead of squared error
5. **Regularization**: Ridge penalty to prevent coefficient explosion
6. **Implementation**: C++ with Eigen library vs pure C

### Additional Features

1. **Adaptive Neighborhood**: Switches from disk to k-NN when insufficient points
2. **Grid Interpolation**: Predictions at data points via linear interpolation from grid
3. **Multiple Error Metrics**: Supports Brier score for model selection
4. **Efficient Search**: Uses ANN library for fast nearest neighbor queries

### Parameters

- `grid_size`: Number of uniform grid points (default varies)
- `fit_quadratic`: Include quadratic term (true/false)
- `pilot_bandwidth`: Fixed bandwidth if > 0
- `kernel_type`: 0=Gaussian, 1=Epanechnikov
- `min_points`: Minimum points for local fitting
- `cv_folds`: Number of CV folds for bandwidth selection
- `n_bws`: Number of bandwidth candidates to test
- `max_iterations`: Newton-Raphson iteration limit
- `ridge_lambda`: Ridge regularization parameter
- `tolerance`: Convergence threshold

---

## Comparison: MAGELO/MAGELOG vs MABILO/MABILOG

### Neighborhood Structure

| Aspect | MAGELO/MAGELOG | MABILO/MABILOG |
|--------|----------------|----------------|
| **Type** | Fixed-radius disk | Symmetric k-hop |
| **Definition** | All points within distance r | k points on each side |
| **Adaptivity** | Radius-based | Count-based |
| **Boundary behavior** | Natural shrinkage | Adjusted windows |
| **Data density sensitivity** | Less sensitive | More sensitive |

### Model Placement

| Aspect | MAGELO/MAGELOG | MABILO/MABILOG |
|--------|----------------|----------------|
| **Model centers** | Uniform grid points | Data points |
| **Number of models** | Fixed (grid_size) | Varies with n |
| **Coverage** | Uniform over domain | Concentrated at data |
| **Extrapolation** | Limited | Not applicable |

### Computational Characteristics

| Aspect | MAGELO/MAGELOG | MABILO/MABILOG |
|--------|----------------|----------------|
| **Time complexity** | O(g × n × d) | O(n × k² × d) |
| **Space complexity** | O(g × k_avg) | O(n × k) |
| **Parallelization** | Grid points independent | Data points independent |
| **Memory pattern** | Predictable (grid-based) | Data-dependent |

### Statistical Properties

| Aspect | MAGELO/MAGELOG | MABILO/MABILOG |
|--------|----------------|----------------|
| **Bias** | Lower in sparse regions | Higher in sparse regions |
| **Variance** | Higher (fewer models) | Lower (more models) |
| **Smoothness** | Very smooth | Moderately smooth |
| **Edge effects** | Handled naturally | Requires adjustment |

### Use Case Recommendations

**Choose MAGELO/MAGELOG when:**
- Data has varying density
- Smooth interpolation is important
- Grid-based representation is natural
- Computational efficiency matters for large n
- Need predictions at many new points

**Choose MABILO/MABILOG when:**
- Data is relatively uniform
- Need high accuracy at data points
- Sample size is moderate
- Bootstrap inference is important
- Working with 1D sequential data

### Implementation Languages

- **MAGELO**: Pure C with manual memory management
- **MAGELOG**: C++ with Eigen library and STL
- **MABILO/MABILOG**: C++ with custom optimization

### Special Features Comparison

| Feature | MAGELO | MAGELOG | MABILO | MABILOG |
|---------|--------|---------|--------|---------|
| **Bootstrap CI** | ✓ | ✗ | ✓ | ✓ |
| **LOOCV** | ✓ | ✓ | ✓ | ✓ |
| **Permutation testing** | ✓ | ✗ | ✗ | ✗ |
| **Auto bandwidth selection** | Via parent R function | ✓ Built-in | ✓ | ✓ |
| **Ridge regularization** | ✗ | ✓ | ✗ | ✓ |
| **Quadratic models** | ✓ | ✓ | ✗ | ✗ |

## Summary

Both algorithm families implement sophisticated model averaging approaches for nonparametric regression:

- **MAGELO/MAGELOG**: Emphasizes geometric neighborhoods and uniform domain coverage through grid-based modeling
- **MABILO/MABILOG**: Focuses on symmetric local neighborhoods and data-centric modeling

The choice between them depends on data characteristics, computational resources, and the specific requirements of the analysis task. MAGELO/MAGELOG excels in scenarios with irregular data distributions and when smooth interpolation is crucial, while MABILO/MABILOG provides superior local accuracy and uncertainty quantification capabilities.
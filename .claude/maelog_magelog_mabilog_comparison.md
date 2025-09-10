# Comprehensive Comparison: MAELOG vs MAGELOG vs MABILOG

This document provides a detailed comparison of three logistic regression methods for binary classification with local model averaging: MAELOG, MAGELOG, and MABILOG.

## Algorithm Overview

### MAELOG (Model-Averaged Exponential LOGistic)
**File**: `/Users/pgajer/current_projects/gflow/src/maelog.cpp`

MAELOG implements model-averaged local logistic regression where:
- Models are centered at **each data point** (not on a grid)
- Uses **adaptive neighborhoods** (disk-based or k-NN fallback)
- Performs **point-wise model averaging**
- Focuses on **direct prediction** at data locations

### MAGELOG (Model-Averaged GEometric LOGistic)
**File**: `/Users/pgajer/current_projects/gflow/src/magelog.cpp`

MAGELOG implements grid-based local logistic regression where:
- Models are centered at **uniform grid points**
- Uses **fixed-radius disk neighborhoods**
- Performs **grid-based model averaging**
- Requires **interpolation** for predictions at data points

### MABILOG (Model-Averaged Bi-kNN LOGistic)
**File**: `/Users/pgajer/current_projects/gflow/src/mabilog.cpp`

MABILOG implements symmetric k-hop local logistic regression where:
- Models are centered at **each data point**
- Uses **symmetric k-hop neighborhoods** (2k+1 points)
- Performs **overlapping model averaging**
- Provides **direct predictions** without interpolation

## Detailed Algorithm Comparison

### 1. Model Placement Strategy

| Method | Model Centers | Number of Models | Coverage |
|--------|--------------|------------------|----------|
| **MAELOG** | Data points | n (number of data points) | Concentrated at data |
| **MAGELOG** | Uniform grid | g (grid_size parameter) | Uniform over domain |
| **MABILOG** | Data points | n (number of data points) | Concentrated at data |

### 2. Neighborhood Definition

| Method | Primary Neighborhood | Fallback Strategy | Adaptivity |
|--------|---------------------|-------------------|------------|
| **MAELOG** | Disk (radius = bandwidth) | k-NN when < min_points | Hybrid adaptive |
| **MAGELOG** | Disk (radius = bandwidth) | k-NN when < min_points | Hybrid adaptive |
| **MABILOG** | Symmetric k-hop | Boundary adjustment | Count-based |

### 3. Algorithm Structure

#### MAELOG Algorithm:
```
For each data point i:
    1. Define neighborhood:
       - Try disk of radius bandwidth around x[i]
       - If < min_points, use k-NN fallback
    2. Fit local logistic model:
       - Center predictors at x[i]
       - Apply kernel weights
       - Newton-Raphson optimization with ridge
    3. Store predictions for all neighbors
    4. Model averaging:
       - Each point receives predictions from multiple models
       - Weight by kernel and average
```

#### MAGELOG Algorithm:
```
Create uniform grid over data range:
For each grid point g:
    1. Define neighborhood:
       - Disk of radius bandwidth around grid[g]
       - If < min_points, use k-NN fallback
    2. Fit local logistic model:
       - Center predictors at grid[g]
       - Apply kernel weights
       - Newton-Raphson optimization with ridge
    3. Store predictions at grid points
    4. Model averaging at grid:
       - Grid points receive predictions from overlapping models
       - Weight by kernel and average
    5. Interpolate to data points:
       - Linear interpolation from grid predictions
```

#### MABILOG Algorithm:
```
For each data point i:
    1. Define symmetric window:
       - [i-k, i+k] for interior points
       - Adjusted for boundaries
    2. Fit local logistic model:
       - Use all 2k+1 points
       - Apply kernel weights
       - Newton-Raphson optimization with ridge
    3. Store model and predictions
    4. Model averaging:
       - Each point covered by ~2k models
       - Weight by kernel position and average
```

### 4. Key Implementation Details

#### MAELOG:
```cpp
// Neighborhood selection (adaptive)
for (int i = 0; i < n_pts; ++i) {
    double dist = std::abs(x[i] - center) / bandwidth;
    if (dist < 1.0) {
        local_indices.push_back(i);
    }
}
if (local_indices.size() < min_points) {
    // Switch to k-NN
    kdtree->annkSearch(query_point, min_points, nn_idx, nn_dists);
}

// Model averaging
for (const auto& pt : v) {
    weighted_prediction += pt.w * pt.p;
    weighted_beta1 += pt.w * pt.beta1;
    total_weight += pt.w;
}
predictions[pti] = weighted_prediction / total_weight;
```

#### MAGELOG:
```cpp
// Grid-based fitting
for (int pti = 0; pti < grid_size; pti++) {
    double center = result.x_grid[pti];
    // Fit local model at grid point
    // Store predictions at grid points
}

// Interpolation to data points
for(int i = 0; i < n_pts; i++) {
    result.predictions[i] = interpolate_grid(x[i], result.x_grid, grid_preds);
}
```

#### MABILOG:
```cpp
// Symmetric window
if (i in interior):
    window = [i-k, i+k]
else if (i near left boundary):
    window = [0, 2k]
else if (i near right boundary):
    window = [n-2k-1, n-1]

// Direct model averaging
For each point i:
    contributing_models = models whose windows contain i
    prediction[i] = weighted_average(contributing_models)
```

### 5. Computational Complexity

| Method | Time Complexity | Space Complexity | Memory Pattern |
|--------|----------------|------------------|----------------|
| **MAELOG** | O(n² × m) | O(n × k_avg) | Data-dependent |
| **MAGELOG** | O(g × n × m) | O(g × k_avg) | Predictable |
| **MABILOG** | O(n × k² × m) | O(n × k) | Predictable |

Where:
- n = number of data points
- g = grid size
- k = neighborhood size
- m = cost of fitting one local model

### 6. Bandwidth Selection

| Method | Selection Method | Error Metric | Implementation |
|--------|-----------------|--------------|----------------|
| **MAELOG** | CV or LOOCV | Brier score | Built-in |
| **MAGELOG** | CV or LOOCV | Brier score | Built-in |
| **MABILOG** | LOOCV | Brier score | Built-in |

### 7. Output and Predictions

| Method | Primary Output | Coefficient Storage | Interpolation Required |
|--------|---------------|--------------------|-----------------------|
| **MAELOG** | Predictions at data points | Beta1, Beta2 arrays | No |
| **MAGELOG** | Grid predictions | Not stored | Yes (linear) |
| **MABILOG** | Predictions at data points | Not stored | No |

## Statistical Properties

### Bias-Variance Trade-off

| Method | Bias | Variance | Smoothness |
|--------|------|----------|------------|
| **MAELOG** | Medium | Low-Medium | Good |
| **MAGELOG** | Low | Higher | Very smooth |
| **MABILOG** | Low-Medium | Low | Moderate |

### Edge Effects

| Method | Boundary Handling | Extrapolation |
|--------|------------------|---------------|
| **MAELOG** | Natural (disk shrinks) | Limited |
| **MAGELOG** | Natural (grid extends) | Better |
| **MABILOG** | Explicit adjustment | None |

## Use Case Recommendations

### Choose MAELOG when:
- You need predictions exactly at data points
- Data density varies significantly
- You want a balance between accuracy and smoothness
- Memory usage needs to be data-adaptive

### Choose MAGELOG when:
- You need predictions at many new points
- Smooth interpolation is important
- Data covers the domain relatively uniformly
- You can afford the computational cost of grid-based fitting

### Choose MABILOG when:
- Working with ordered/sequential data
- Need symmetric neighborhoods
- Want predictable memory usage
- Bootstrap inference is important

## Key Differences Summary

### 1. **Model Center Philosophy**
- **MAELOG & MABILOG**: Data-centric (models at data points)
- **MAGELOG**: Domain-centric (models on uniform grid)

### 2. **Neighborhood Type**
- **MAELOG & MAGELOG**: Geometric (disk-based with k-NN fallback)
- **MABILOG**: Topological (k-hop based)

### 3. **Prediction Strategy**
- **MAELOG & MABILOG**: Direct prediction at data points
- **MAGELOG**: Grid prediction + interpolation

### 4. **Implementation Language**
- **MAELOG & MAGELOG**: C++ with Eigen library
- **MABILOG**: C++ with custom optimization

### 5. **Memory Scaling**
- **MAELOG**: O(n × average_neighborhood_size)
- **MAGELOG**: O(grid_size × average_neighborhood_size)
- **MABILOG**: O(n × k)

## Performance Characteristics

### Speed (fastest to slowest):
1. **MAGELOG** (when g << n)
2. **MABILOG** (predictable k)
3. **MAELOG** (variable neighborhoods)

### Memory Efficiency:
1. **MAGELOG** (controlled by grid_size)
2. **MABILOG** (predictable)
3. **MAELOG** (data-dependent)

### Smoothness:
1. **MAGELOG** (grid interpolation)
2. **MAELOG** (disk averaging)
3. **MABILOG** (k-hop averaging)

### Accuracy at Data Points:
1. **MABILOG** (direct fitting)
2. **MAELOG** (direct fitting)
3. **MAGELOG** (interpolation error)

## Implementation Notes

### Common Features:
- All use Newton-Raphson optimization for logistic regression
- All support ridge regularization
- All implement kernel weighting (multiple kernel types)
- All have automatic bandwidth selection via CV
- All handle sparse regions with k-NN fallback (MAELOG, MAGELOG) or window adjustment (MABILOG)

### Unique Features:
- **MAELOG**: Stores beta coefficients for each point
- **MAGELOG**: Grid-based structure enables efficient new predictions
- **MABILOG**: Symmetric neighborhoods provide balanced local models

## Conclusion

The three methods represent different philosophies in local logistic regression:

- **MAELOG**: Adaptive, data-centric approach with flexible neighborhoods
- **MAGELOG**: Systematic, grid-based approach with uniform coverage
- **MABILOG**: Structured, symmetric approach with predictable behavior

The choice depends on:
1. Whether predictions are needed at data points only (MAELOG, MABILOG) or throughout the domain (MAGELOG)
2. Whether data is sequential/ordered (MABILOG) or scattered (MAELOG, MAGELOG)
3. Computational resources and scalability requirements
4. Importance of smoothness vs. accuracy at data points
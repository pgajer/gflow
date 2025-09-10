# MABILO and MABILOG Algorithm Documentation

This document provides detailed algorithmic outlines for the MABILO (Model-Averaged Bi-kNN LOcal linear) and MABILOG (Model-Averaged Bi-kNN LOcal logistic) regression methods implemented in the gflow package.

## MABILO Algorithm (uwmabilo)

### Overview
MABILO extends traditional LOWESS by incorporating model averaging over symmetric k-hop neighborhoods with kernel-weighted contributions. The algorithm fits local linear models and combines predictions through weighted averaging.

### Algorithm Structure

#### Phase 1: Model Computation and Storage
For each k value in [k_min, k_max]:

1. **Initialize data structures**
   - Pre-allocate model storage for each point: `pt_models[n_points][2k+1]`
   - Pre-allocate error vectors for LOOCV and true errors

2. **For each data point i:**
   
   a. **Determine symmetric window**
   ```
   if (i in interior):
       window = [i-k, i+k]  // Symmetric k-hop neighborhood
   else if (i near left boundary):
       window = [0, 2k]     // Adjusted for boundary
   else if (i near right boundary):
       window = [n-2k-1, n-1]  // Adjusted for boundary
   ```
   
   b. **Compute kernel weights for window**
   ```
   - Calculate distances: d[j] = |x[j] - x[i]| for j in window
   - Normalize by max distance: d[j] = d[j] / (max_dist * normalization_factor)
   - Apply kernel function: w[j] = kernel(d[j])
   - Normalize weights: w[j] = w[j] / sum(w)
   - Apply sample weights (if provided): w[j] = w[j] * sample_weight[j]
   ```
   
   c. **Fit weighted local linear model**
   ```
   - Build design matrix X = [1, x_local]
   - Apply weights to observations
   - Solve weighted least squares: β = (X'WX)^(-1)X'Wy
   - Compute predictions for all points in window
   - Calculate LOOCV error at point i
   ```
   
   d. **Store model for later use**
   ```
   pt_models[i].push_back(model)
   ```

#### Phase 2: Model Averaging
For each k value:

1. **For each prediction point i:**
   
   a. **Collect contributing models**
   ```
   - Find all models whose windows contain point i
   - For each model m with window [start, end]:
       if (start <= i <= end):
           contributing_models.add(m)
   ```
   
   b. **Compute model weights**
   ```
   - Extract kernel weight from each model for point i
   - Normalize weights across contributing models
   ```
   
   c. **Compute weighted average prediction**
   ```
   prediction[i] = sum(model_weight[m] * model_prediction[m][i])
                   for all contributing models m
   ```

2. **Calculate mean errors for this k**
   ```
   - LOOCV error: mean of squared differences between observed and predicted
   - True error (if available): mean absolute error
   ```

#### Phase 3: Optimal k Selection
1. **Find k with minimum mean LOOCV error**
   ```
   opt_k = argmin(k_mean_errors)
   ```

2. **Return results**
   ```
   - Optimal k value
   - Predictions using optimal k
   - All k predictions for analysis
   - Error metrics for each k
   ```

### Key Algorithmic Features

1. **Symmetric Neighborhoods**: Uses k points on each side (2k+1 total) rather than k nearest neighbors
2. **Boundary Handling**: Adjusts windows near data boundaries to maintain model stability
3. **Kernel Weighting**: Distance-based weights using tricube, Epanechnikov, or exponential kernels
4. **Model Averaging**: Each point's prediction is a weighted average of all models whose windows contain it
5. **LOOCV Error**: Leave-one-out cross-validation for unbiased error estimation

### Computational Complexity
- Time: O(n × k_max² × m) where m is the cost of fitting a local linear model
- Space: O(n × k_max) for storing models and predictions

---

## MABILOG Algorithm (uwmabilog)

### Overview
MABILOG adapts the MABILO framework for binary classification by replacing local linear models with local logistic regression. The algorithm maintains the same model averaging structure while handling binary outcomes.

### Algorithm Structure

#### Phase 1: Model Computation and Storage
For each k value in [k_min, k_max]:

1. **Initialize data structures** (same as MABILO)

2. **For each data point i:**
   
   a. **Determine symmetric window** (same as MABILO)
   
   b. **Compute kernel weights** (same as MABILO)
   
   c. **Fit weighted local logistic model**
   ```
   - Initialize β = [β₀, β₁] using weighted linear regression
   - Iterate until convergence:
       * Compute probabilities: p[j] = 1/(1 + exp(-X[j]β))
       * Compute gradient: g = X'W(y - p)
       * Compute Hessian: H = X'W diag(p(1-p)) WX + λI (ridge)
       * Update: β = β + H^(-1)g
       * Check convergence: ||β_new - β_old|| < tolerance
   - Apply max_beta constraint if needed
   - Compute predictions: p[j] = 1/(1 + exp(-X[j]β))
   - Calculate LOOCV error using log-loss
   ```
   
   d. **Store model** (same as MABILO)

#### Phase 2: Model Averaging
For each k value:

1. **For each prediction point i:**
   
   a. **Collect contributing models** (same as MABILO)
   
   b. **Compute model weights** (same as MABILO)
   
   c. **Compute weighted average probability**
   ```
   probability[i] = sum(model_weight[m] * model_probability[m][i])
                    for all contributing models m
   ```

2. **Calculate mean errors**
   ```
   - LOOCV error: mean log-loss or classification error
   - True error (if available): classification accuracy
   ```

#### Phase 3: Optimal k Selection (same as MABILO)

### Key Differences from MABILO

1. **Model Type**: Local logistic regression instead of linear regression
2. **Optimization**: Iterative Newton-Raphson instead of closed-form solution
3. **Regularization**: Ridge penalty (λ) to prevent coefficient explosion
4. **Constraints**: Maximum beta constraint to ensure numerical stability
5. **Error Metric**: Log-loss or classification error instead of squared error
6. **Predictions**: Probabilities in [0,1] instead of continuous values

### Additional Parameters

- `max_iterations`: Maximum Newton-Raphson iterations (default: 25)
- `ridge_lambda`: Ridge regularization parameter (default: 1e-6)
- `max_beta`: Maximum absolute value for coefficients (default: 10)
- `tolerance`: Convergence threshold (default: 1e-8)

### Computational Complexity
- Time: O(n × k_max² × m × iter) where iter is the average number of Newton-Raphson iterations
- Space: O(n × k_max) for storing models and predictions

---

## Comparison with Other Univariate Methods

### Structural Comparison

| Method | Neighborhood | Model Type | Averaging | Bandwidth |
|--------|-------------|------------|-----------|-----------|
| **MABILO** | Symmetric k-hop | Local linear | Yes, kernel-weighted | Adaptive via CV |
| **MABILOG** | Symmetric k-hop | Local logistic | Yes, kernel-weighted | Adaptive via CV |
| **MAGELO** | Fixed-radius disks | Local polynomial | Yes, grid-based | Fixed radius |
| **AMAGELO** | Adaptive radius | Local linear | Yes, grid-based | Adaptive |
| **Traditional LOWESS** | k-nearest | Local linear/quadratic | No | Fixed k |

### Algorithmic Features

| Feature | MABILO/MABILOG | MAGELO | AMAGELO | LOWESS |
|---------|----------------|--------|---------|---------|
| **Model centers** | Data points | Grid points | Grid points | Data points |
| **Window symmetry** | Symmetric | Symmetric | Adaptive | Asymmetric possible |
| **Model averaging** | Yes | Yes | Yes | No |
| **LOOCV selection** | Yes | Yes | Yes | Optional |
| **Boundary handling** | Adjusted windows | Natural | Natural | Adjusted |
| **Computational scaling** | O(n²k) | O(ng) | O(ng) | O(nk) |

### Advantages of MABILO/MABILOG

1. **Symmetric Neighborhoods**: More balanced local models in 1D
2. **Model Averaging**: Smoother predictions with reduced variance
3. **Automatic Bandwidth Selection**: Data-driven optimal k via LOOCV
4. **Robustness**: Multiple models contribute to each prediction
5. **Uncertainty Quantification**: Natural framework for bootstrap inference

### Use Case Recommendations

- **MABILO**: General nonparametric regression with continuous outcomes
- **MABILOG**: Binary classification with smooth probability estimates
- **MAGELO**: When data density varies significantly across domain
- **AMAGELO**: When adaptive bandwidth is crucial for local features
- **Traditional LOWESS**: Simple, fast smoothing without model averaging

### Implementation Notes

Both algorithms share:
- Common kernel weight computation
- Similar window determination logic
- Parallel bootstrap capability
- Efficient memory management with pre-allocation
- Progress tracking for long computations

The main implementation difference is in the model fitting step:
- MABILO uses closed-form weighted least squares
- MABILOG uses iterative weighted logistic regression with Newton-Raphson
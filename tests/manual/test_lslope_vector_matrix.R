## test_lslope_vector_matrix.R
##
## Test script to compare the R loop implementation vs C++ implementation
## of lslope.grad.vector.matrix
##
## Run this after integrating the C++ code into the gflow package.

## library(gflow)

## Create a test graph (random geometric graph on 200 vertices)
set.seed(42)
n <- 200
coords <- matrix(runif(n * 3), ncol = 3)  # 3D for variety

## Build k-NN graph
k <- 10
dist.mat <- as.matrix(dist(coords))

adj.list <- vector("list", n)
weight.list <- vector("list", n)

for (i in 1:n) {
    dists <- dist.mat[i, ]
    dists[i] <- Inf  # exclude self
    neighbors <- order(dists)[1:k]
    adj.list[[i]] <- as.integer(neighbors)
    weight.list[[i]] <- dists[neighbors]
}

## Make symmetric
for (i in 1:n) {
    for (j in seq_along(adj.list[[i]])) {
        neighbor <- adj.list[[i]][j]
        weight <- weight.list[[i]][j]
        
        if (!(i %in% adj.list[[neighbor]])) {
            adj.list[[neighbor]] <- c(adj.list[[neighbor]], as.integer(i))
            weight.list[[neighbor]] <- c(weight.list[[neighbor]], weight)
        }
    }
}

## Create test data
## y: directing function (smooth gradient + noise)
y <- coords[, 1] + 0.2 * coords[, 2] + rnorm(n, sd = 0.1)

## Z: response matrix with various patterns
n.cols <- 20
Z <- matrix(NA, nrow = n, ncol = n.cols)
colnames(Z) <- paste0("Z", 1:n.cols)

for (j in 1:n.cols) {
    if (j <= 5) {
        ## Positively correlated with y
        Z[, j] <- y + rnorm(n, sd = 0.3)
    } else if (j <= 10) {
        ## Negatively correlated with y
        Z[, j] <- -y + rnorm(n, sd = 0.3)
    } else if (j <= 15) {
        ## Uncorrelated
        Z[, j] <- rnorm(n)
    } else {
        ## Compositional (relative abundances)
        Z[, j] <- exp(y * runif(1, -1, 1)) + runif(n, 0.01, 0.1)
        Z[, j] <- Z[, j] / sum(Z[, j])  # normalize to sum to 1
    }
}

cat("Test graph: n =", n, "vertices, k =", k, "neighbors\n")
cat("Response matrix Z:", n, "x", n.cols, "\n\n")

## ============================================================================
## Test 1: Compare R loop vs C++ for type = "normalized"
## ============================================================================

cat("Test 1: type = 'normalized', diff.type = 'difference'\n")
cat("=" , rep("-", 50), "=\n", sep = "")

## R loop implementation (from lslope_extended.R)
## This calls S_lslope_gradient for each column
t1.r <- system.time({
    result.r <- lslope.grad.vector.matrix( ## this is C++ implementation and R implementation is lslope.grad.vector.matrix.R()
        adj.list, weight.list, y, Z,
        type = "normalized",
        y.diff.type = "difference",
        z.diff.type = "difference",
        epsilon = 0,
        sigmoid.alpha = 0,
        ascending = TRUE
    )
})

cat("R loop time:", t1.r["elapsed"], "seconds\n")
cat("Dimensions:", dim(result.r), "\n")
cat("Sigmoid alpha:", attr(result.r, "sigmoid.alpha"), "\n")
cat("Local extrema:", attr(result.r, "n.local.maxima") + attr(result.r, "n.local.minima"), "\n\n")

## C++ implementation (once integrated)
## Uncomment when S_lslope_vector_matrix is available:
#
# t1.cpp <- system.time({
#     result.cpp <- lslope.grad.vector.matrix.cpp(
#         adj.list, weight.list, y, Z,
#         type = "normalized",
#         y.diff.type = "difference",
#         z.diff.type = "difference",
#         epsilon = 0,
#         sigmoid.alpha = 0,
#         ascending = TRUE,
#         n.threads = 4L
#     )
# })
# 
# cat("C++ time:", t1.cpp["elapsed"], "seconds\n")
# cat("Speedup:", t1.r["elapsed"] / t1.cpp["elapsed"], "x\n\n")
#
# ## Compare results
# max.diff <- max(abs(result.r - result.cpp))
# cat("Maximum difference between R and C++:", max.diff, "\n")
# 
# if (max.diff < 1e-10) {
#     cat("✓ PASSED: Results match within numerical tolerance\n\n")
# } else {
#     cat("✗ FAILED: Results differ significantly\n")
#     cat("  Correlation:", cor(as.vector(result.r), as.vector(result.cpp)), "\n\n")
# }

## ============================================================================
## Test 2: type = "slope" (raw, unbounded)
## ============================================================================

cat("Test 2: type = 'slope', diff.type = 'difference'\n")
cat("=" , rep("-", 50), "=\n", sep = "")

result.slope <- lslope.grad.vector.matrix(
    adj.list, weight.list, y, Z,
    type = "slope",
    y.diff.type = "difference",
    z.diff.type = "difference"
)

cat("Range of coefficients:", range(result.slope), "\n")
cat("Column means (first 5):", round(colMeans(result.slope)[1:5], 4), "\n\n")

## ============================================================================
## Test 3: type = "sign" 
## ============================================================================

cat("Test 3: type = 'sign', diff.type = 'difference'\n")
cat("=" , rep("-", 50), "=\n", sep = "")

result.sign <- lslope.grad.vector.matrix(
    adj.list, weight.list, y, Z,
    type = "sign",
    y.diff.type = "difference",
    z.diff.type = "difference"
)

cat("Unique values:", sort(unique(as.vector(result.sign))), "\n")
cat("Proportion positive (col 1):", mean(result.sign[, 1] > 0), "\n")
cat("Proportion positive (col 6):", mean(result.sign[, 6] > 0), "\n\n")

## ============================================================================
## Test 4: Compositional data with logratio
## ============================================================================

cat("Test 4: Compositional data with logratio\n")
cat("=" , rep("-", 50), "=\n", sep = "")

result.comp <- lslope.grad.vector.matrix(
    adj.list, weight.list, y, Z[, 16:20],
    type = "normalized",
    y.diff.type = "difference",
    z.diff.type = "logratio",
    epsilon = 0
)

cat("Compositional columns (16-20) with logratio\n")
cat("Column means:", round(colMeans(result.comp), 4), "\n\n")

## ============================================================================
## Test 5: Verify consistency with vector-vector lslope.grad
## ============================================================================

cat("Test 5: Verify consistency with vector-vector lslope.grad\n")
cat("=" , rep("-", 50), "=\n", sep = "")

## Compute for first column using vector-vector
result.vv <- lslope.grad(
    adj.list, weight.list, y, Z[, 1],
    type = "normalized",
    y.diff.type = "difference",
    z.diff.type = "difference",
    sigmoid.alpha = attr(result.r, "sigmoid.alpha"),  # Use same alpha
    ascending = TRUE
)

## Compare
diff.col1 <- max(abs(result.r[, 1] - result.vv))
cat("Max diff between vector-matrix[,1] and vector-vector:", diff.col1, "\n")

if (diff.col1 < 1e-10) {
    cat("✓ PASSED: Vector-matrix is consistent with vector-vector\n\n")
} else {
    cat("✗ FAILED: Inconsistency detected\n\n")
}

## ============================================================================
## Test 6: Matrix-matrix
## ============================================================================

cat("Test 6: Matrix-matrix lslope.grad\n")
cat("=" , rep("-", 50), "=\n", sep = "")

## Small test: 2 directing functions, 5 response functions
Y.small <- cbind(y, -y + rnorm(n, sd = 0.1))
colnames(Y.small) <- c("Y1", "Y2")

result.mm <- lslope.grad.matrix.matrix(
    adj.list, weight.list,
    Y.small, Z[, 1:5],
    type = "normalized",
    mc.cores = 1
)

cat("Dimensions:", dim(result.mm), "\n")
cat("Array structure: [vertices, Y-cols, Z-cols]\n")
cat("Mean coefficients for Y1 vs Z1-5:", round(colMeans(result.mm[, 1, ]), 4), "\n")
cat("Mean coefficients for Y2 vs Z1-5:", round(colMeans(result.mm[, 2, ]), 4), "\n\n")

## ============================================================================
## Summary
## ============================================================================

cat("=" , rep("=", 50), "=\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 50), "=\n", sep = "")
cat("All basic tests completed.\n")
cat("To compare R vs C++ implementation:\n")
cat("1. Integrate lslope_vector_matrix.cpp into gflow\n")
cat("2. Add lslope_vector_matrix_r.cpp and .h\n")
cat("3. Register S_lslope_vector_matrix in init.c\n")
cat("4. Rebuild package\n")
cat("5. Uncomment the C++ test section above\n")

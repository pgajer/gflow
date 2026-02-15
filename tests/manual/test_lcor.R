# =============================================================================
# LCOR VALIDATION TEST SUITE (R Version)
# =============================================================================
# 
# Comprehensive validation of lcor() implementation by comparing all output
# fields against hand-calculated expected values.
#
# Run with: source("tests/manual/test_lcor_validation.R")
#
# Each test includes:
# - Expected vertex.coefficients
# - Expected vertex.delta.y, vertex.delta.z, vertex.weights
# - Expected all.delta.y, all.delta.z
# - Expected winsorization bounds (when applicable)
# =============================================================================

## Assumes package is loaded with devtools::load_all()

# Helper functions ============================================================

# Tolerance for floating-point comparisons
TOLERANCE <- 1e-10

# Check if two numbers are approximately equal
approx_equal <- function(a, b, tol = TOLERANCE) {
  abs(a - b) < tol
}

# Check if two vectors are approximately equal
vectors_equal <- function(a, b, tol = TOLERANCE) {
  if (length(a) != length(b)) return(FALSE)
  all(abs(a - b) < tol)
}

# Print test header
print_test_header <- function(test_name) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("TEST:", test_name, "\n")
  cat(rep("=", 70), "\n", sep = "")
}

# Print pass/fail
print_pass <- function() {
  cat("\u2713 PASS\n", sep = "")
}

print_fail <- function(message) {
  cat("\u2717 FAIL:", message, "\n", sep = "")
}

# Test result tracking
test_results <- list()

# =============================================================================
# TEST 1: Simple 3-vertex path with unit weights, standard differences
# =============================================================================
#
# Graph: 0 -- 1 -- 2 (using 1-based indexing in R)
#        1 -- 2 -- 3
# Edges: [1,2] with length 1.0
#        [2,3] with length 1.0
#
# Functions: y = c(1, 2, 3), z = c(2, 4, 6)
# Type: unit (w_e = 1)
# Diff types: both difference
#
# Manual calculations (using 0-based internal indexing for clarity):
#
# Vertex 0 (R index 1):
#   Neighbors: [1] (R: [2])
#   Edge [0,1]: Δy = y[1] - y[0] = 2 - 1 = 1
#               Δz = z[1] - z[0] = 4 - 2 = 2
#   weight = 1
#   numerator = 1 * 1 * 2 = 2
#   sum_y² = 1 * 1² = 1
#   sum_z² = 1 * 2² = 4
#   coeff[0] = 2 / (√1 * √4) = 2/2 = 1.0
#
# Vertex 1 (R index 2):
#   Neighbors: [0, 2] (R: [1, 3])
#   Edge [1,0]: Δy = y[0] - y[1] = 1 - 2 = -1
#               Δz = z[0] - z[1] = 2 - 4 = -2
#   Edge [1,2]: Δy = y[2] - y[1] = 3 - 2 = 1
#               Δz = z[2] - z[1] = 6 - 4 = 2
#   weights = [1, 1]
#   numerator = 1*(-1)*(-2) + 1*1*2 = 2 + 2 = 4
#   sum_y² = 1*1 + 1*1 = 2
#   sum_z² = 1*4 + 1*4 = 8
#   coeff[1] = 4 / (√2 * √8) = 4/4 = 1.0
#
# Vertex 2 (R index 3):
#   Neighbors: [1] (R: [2])
#   Edge [2,1]: Δy = y[1] - y[2] = 2 - 3 = -1
#               Δz = z[1] - z[2] = 4 - 6 = -2
#   weight = 1
#   numerator = 1 * (-1) * (-2) = 2
#   sum_y² = 1 * 1 = 1
#   sum_z² = 1 * 4 = 4
#   coeff[2] = 2 / (√1 * √4) = 1.0
# =============================================================================

test_simple_path_unit_weights <- function() {
  print_test_header("Simple 3-vertex path, unit weights, standard differences")
  
  # Build graph: 1 -- 2 -- 3 (R 1-based indexing)
  adj.list <- list(
    c(2),      # vertex 1 neighbors
    c(1, 3),   # vertex 2 neighbors
    c(2)       # vertex 3 neighbors
  )
  
  weight.list <- list(
    c(1.0),       # edges from vertex 1
    c(1.0, 1.0),  # edges from vertex 2
    c(1.0)        # edges from vertex 3
  )
  
  # Define functions
  y <- c(1.0, 2.0, 3.0)
  z <- c(2.0, 4.0, 6.0)
  
  # Compute lcor
  result <- lcor(
    adj.list, weight.list, y, z,
    type = "unit",
    y.diff.type = "difference",
    z.diff.type = "difference",
    epsilon = 0,
    winsorize.quantile = 0
  )
  
  # Expected values
  expected_coeffs <- c(1.0, 1.0, 1.0)
  
  expected_delta_y <- list(
    c(1.0),         # vertex 1
    c(-1.0, 1.0),   # vertex 2
    c(-1.0)         # vertex 3
  )
  
  expected_delta_z <- list(
    c(2.0),         # vertex 1
    c(-2.0, 2.0),   # vertex 2
    c(-2.0)         # vertex 3
  )
  
  expected_weights <- list(
    c(1.0),         # vertex 1
    c(1.0, 1.0),    # vertex 2
    c(1.0)          # vertex 3
  )
  
  # Note: all.delta.y/z order depends on iteration order (vertex 1, then 2, then 3)
  expected_all_delta_y <- c(1.0, -1.0, 1.0, -1.0)
  expected_all_delta_z <- c(2.0, -2.0, 2.0, -2.0)
  
  # Verify vertex coefficients
  coeffs_ok <- TRUE
  for (i in 1:3) {
    if (!approx_equal(result$vertex.coefficients[i], expected_coeffs[i])) {
      print_fail(sprintf("vertex.coefficients[%d] = %.6f, expected %.6f",
                        i, result$vertex.coefficients[i], expected_coeffs[i]))
      coeffs_ok <- FALSE
    }
  }
  if (coeffs_ok) {
    cat("\u2713 vertex.coefficients correct\n")
  }
  
  # Verify vertex.delta.y
  delta_y_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.delta.y[[v]], expected_delta_y[[v]])) {
      print_fail(sprintf("vertex.delta.y[[%d]] mismatch", v))
      delta_y_ok <- FALSE
    }
  }
  if (delta_y_ok) {
    cat("\u2713 vertex.delta.y correct\n")
  }
  
  # Verify vertex.delta.z
  delta_z_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.delta.z[[v]], expected_delta_z[[v]])) {
      print_fail(sprintf("vertex.delta.z[[%d]] mismatch", v))
      delta_z_ok <- FALSE
    }
  }
  if (delta_z_ok) {
    cat("\u2713 vertex.delta.z correct\n")
  }
  
  # Verify vertex.weights
  weights_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.weights[[v]], expected_weights[[v]])) {
      print_fail(sprintf("vertex.weights[[%d]] mismatch", v))
      weights_ok <- FALSE
    }
  }
  if (weights_ok) {
    cat("\u2713 vertex.weights correct\n")
  }
  
  # Verify all.delta.y
  all_delta_y_ok <- TRUE
  if (length(result$all.delta.y) != length(expected_all_delta_y)) {
    print_fail(sprintf("all.delta.y size mismatch: got %d, expected %d",
                      length(result$all.delta.y), length(expected_all_delta_y)))
    all_delta_y_ok <- FALSE
  } else if (!vectors_equal(result$all.delta.y, expected_all_delta_y)) {
    print_fail("all.delta.y values mismatch")
    all_delta_y_ok <- FALSE
  }
  if (all_delta_y_ok) {
    cat("\u2713 all.delta.y correct\n")
  }
  
  # Verify all.delta.z
  all_delta_z_ok <- TRUE
  if (length(result$all.delta.z) != length(expected_all_delta_z)) {
    print_fail("all.delta.z size mismatch")
    all_delta_z_ok <- FALSE
  } else if (!vectors_equal(result$all.delta.z, expected_all_delta_z)) {
    print_fail("all.delta.z values mismatch")
    all_delta_z_ok <- FALSE
  }
  if (all_delta_z_ok) {
    cat("\u2713 all.delta.z correct\n")
  }
  
  all_pass <- coeffs_ok && delta_y_ok && delta_z_ok && 
              weights_ok && all_delta_y_ok && all_delta_z_ok
  
  if (all_pass) {
    print_pass()
  }
  
  return(all_pass)
}

# =============================================================================
# TEST 2: Path with derivative weights and log-ratios
# =============================================================================
#
# Graph: Same as Test 1 but with varying edge lengths
# Edge lengths: [1,2] = 1.0, [2,3] = 2.0
#
# Functions: y = c(1, 2, 4), z = c(0.01, 0.10, 0.50)
# Type: derivative (w_e = 1/ℓ²)
# y.diff.type: difference
# z.diff.type: logratio with adaptive epsilon
#
# Adaptive epsilon for z: min_nonzero = 0.01, ε_z = 1e-6 * 0.01 = 1e-8
# (essentially negligible for the ratios involved)
#
# Manual calculations:
#
# Vertex 0 (R index 1):
#   Edge [0,1]: length = 1.0, weight = 1/1² = 1.0
#               Δy = 2 - 1 = 1
#               Δz = log((0.10)/(0.01)) = log(10) ≈ 2.302585
#   numerator = 1.0 * 1 * 2.302585 = 2.302585
#   sum_y² = 1.0 * 1² = 1.0
#   sum_z² = 1.0 * (2.302585)² ≈ 5.301903
#   coeff[0] = 2.302585 / (1.0 * √5.301903) = 2.302585 / 2.302585 = 1.0
#
# Vertex 1 (R index 2):
#   Edge [1,0]: length = 1.0, weight = 1.0
#               Δy = 1 - 2 = -1
#               Δz = log((0.01)/(0.10)) = log(0.1) ≈ -2.302585
#   Edge [1,2]: length = 2.0, weight = 1/4 = 0.25
#               Δy = 4 - 2 = 2
#               Δz = log((0.50)/(0.10)) = log(5) ≈ 1.609438
#   numerator = 1.0*(-1)*(-2.302585) + 0.25*2*1.609438
#             = 2.302585 + 0.804719 = 3.107304
#   sum_y² = 1.0*1 + 0.25*4 = 2.0
#   sum_z² = 1.0*(2.302585)² + 0.25*(1.609438)²
#          = 5.301903 + 0.647580 = 5.949483
#   coeff[1] = 3.107304 / (√2.0 * √5.949483)
#            = 3.107304 / (1.414214 * 2.439292)
#            = 3.107304 / 3.450327 ≈ 0.900631
#
# Vertex 2 (R index 3):
#   Edge [2,1]: length = 2.0, weight = 0.25
#               Δy = 2 - 4 = -2
#               Δz = log((0.10)/(0.50)) = log(0.2) ≈ -1.609438
#   numerator = 0.25 * (-2) * (-1.609438) = 0.804719
#   sum_y² = 0.25 * 4 = 1.0
#   sum_z² = 0.25 * (1.609438)² = 0.647580
#   coeff[2] = 0.804719 / (√1.0 * √0.647580) = 0.804719 / 0.804719 = 1.0
# =============================================================================

test_derivative_weights_logratio <- function() {
  print_test_header("Path with derivative weights and log-ratios")
  
  # Build graph with different edge lengths
  adj.list <- list(
    c(2),
    c(1, 3),
    c(2)
  )
  
  weight.list <- list(
    c(1.0),
    c(1.0, 2.0),  # Edge to vertex 3 has length 2.0
    c(2.0)
  )
  
  y <- c(1.0, 2.0, 4.0)
  z <- c(0.01, 0.10, 0.50)
  
  result <- lcor(
    adj.list, weight.list, y, z,
    type = "derivative",
    y.diff.type = "difference",
    z.diff.type = "logratio",
    epsilon = 0,  # adaptive
    winsorize.quantile = 0
  )
  
  # Expected values
  expected_coeffs <- c(1.0, 0.900803, 1.0)

  log10 <- log(10.0)
  log5 <- log(5.0)
  log_one_fifth <- log(0.2)
  
  expected_delta_y <- list(
    c(1.0),
    c(-1.0, 2.0),
    c(-2.0)
  )
  
  expected_delta_z <- list(
    c(log10),
    c(-log10, log5),
    c(log_one_fifth)
  )
  
  expected_weights <- list(
    c(1.0),
    c(1.0, 0.25),
    c(0.25)
  )
  
  # Verify coefficients (with slightly larger tolerance for accumulated rounding)
  coeffs_ok <- TRUE
  for (i in 1:3) {
    if (!approx_equal(result$vertex.coefficients[i], expected_coeffs[i], 1e-6)) {
      print_fail(sprintf("vertex.coefficients[%d] = %.6f, expected %.6f",
                        i, result$vertex.coefficients[i], expected_coeffs[i]))
      coeffs_ok <- FALSE
    }
  }
  if (coeffs_ok) {
    cat("\u2713 vertex.coefficients correct\n")
  }
  
  # Verify delta_y (should be exact)
  delta_y_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.delta.y[[v]], expected_delta_y[[v]])) {
      print_fail(sprintf("vertex.delta.y[[%d]] mismatch", v))
      delta_y_ok <- FALSE
    }
  }
  if (delta_y_ok) {
    cat("\u2713 vertex.delta.y correct\n")
  }
  
  # Verify delta_z (log-ratios, small tolerance)
  delta_z_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.delta.z[[v]], expected_delta_z[[v]], 1e-6)) {
      print_fail(sprintf("vertex.delta.z[[%d]] mismatch", v))
      delta_z_ok <- FALSE
    }
  }
  if (delta_z_ok) {
    cat("\u2713 vertex.delta.z correct\n")
  }
  
  # Verify weights
  weights_ok <- TRUE
  for (v in 1:3) {
    if (!vectors_equal(result$vertex.weights[[v]], expected_weights[[v]])) {
      print_fail(sprintf("vertex.weights[[%d]] mismatch", v))
      weights_ok <- FALSE
    }
  }
  if (weights_ok) {
    cat("\u2713 vertex.weights correct\n")
  }
  
  all_pass <- coeffs_ok && delta_y_ok && delta_z_ok && weights_ok
  
  if (all_pass) {
    print_pass()
  }
  
  return(all_pass)
}

# =============================================================================
# TEST 3: Winsorization test with outliers
# =============================================================================
#
# Graph: 1 -- 2 -- 3 -- 4 (4-vertex path)
# Edge lengths: all 1.0
#
# Functions: y = c(1, 2, 100, 3) (outlier at vertex 3)
#            z = c(1, 2, 3, 4)
# Type: unit
# Diff types: both difference
# winsorize.quantile: 0.05
#
# Edge differences (before winsorization):
# From vertex 1: Δy = 2-1 = 1, Δz = 2-1 = 1
# From vertex 2: Δy = 1-2 = -1, Δy = 100-2 = 98, Δz = 1-2 = -1, Δz = 3-2 = 1
# From vertex 3: Δy = 2-100 = -98, Δy = 3-100 = -97, Δz = 2-3 = -1, Δz = 4-3 = 1
# From vertex 4: Δy = 100-3 = 97, Δz = 3-4 = -1
#
# all.delta.y = c(1, -1, 98, -98, -97, 97)
# all.delta.z = c(1, -1, 1, -1, 1, -1)
#
# Sorted all.delta.y = c(-98, -97, -1, 1, 97, 98)
# n = 6
# lower_idx = floor(6 * 0.05) = floor(0.3) = 0 (R: 1) → y.lower = -98
# upper_idx = ceil(6 * 0.95) = ceil(5.7) = 6 (R: 6) - 1 + 1 = 6 → y.upper = 98
#
# Sorted all.delta.z = c(-1, -1, -1, 1, 1, 1)
# z.lower = -1, z.upper = 1
#
# After winsorization (clipping to bounds):
# Vertex 1: Δy = c(1)           (unchanged)
# Vertex 2: Δy = c(-1, 98)      (98 at upper bound)
# Vertex 3: Δy = c(-98, -98)    (-97 clipped to -98)
# Vertex 4: Δy = c(97)          (97 within bounds)
#
# Final coefficients (after winsorization):
#
# Vertex 1:
#   Δy = c(1), Δz = c(1), w = c(1)
#   numerator = 1*1*1 = 1
#   sum_y² = 1*1 = 1, sum_z² = 1*1 = 1
#   coeff[1] = 1 / (1*1) = 1.0
#
# Vertex 2:
#   Δy = c(-1, 98), Δz = c(-1, 1), w = c(1, 1)
#   numerator = 1*(-1)*(-1) + 1*98*1 = 1 + 98 = 99
#   sum_y² = 1*1 + 1*98² = 1 + 9604 = 9605
#   sum_z² = 1*1 + 1*1 = 2
#   coeff[2] = 99 / (√9605 * √2) = 99 / (98.009 * 1.414) ≈ 0.714286
#
# Vertex 3:
#   Δy = c(-98, -98), Δz = c(-1, 1), w = c(1, 1)
#   numerator = 1*(-98)*(-1) + 1*(-98)*1 = 98 - 98 = 0
#   coeff[3] = 0
#
# Vertex 4:
#   Δy = c(97), Δz = c(-1), w = c(1)
#   numerator = 1*97*(-1) = -97
#   sum_y² = 1*97² = 9409, sum_z² = 1*1 = 1
#   coeff[4] = -97 / (√9409 * 1) = -97 / 97 = -1.0
# =============================================================================

test_winsorization <- function() {
  print_test_header("Winsorization with outliers")
  
  # Build 4-vertex path graph
  adj.list <- list(
    c(2),
    c(1, 3),
    c(2, 4),
    c(3)
  )
  
  weight.list <- list(
    c(1.0),
    c(1.0, 1.0),
    c(1.0, 1.0),
    c(1.0)
  )
  
  y <- c(1.0, 2.0, 100.0, 3.0)  # Outlier at vertex 3
  z <- c(1.0, 2.0, 3.0, 4.0)
  
  result <- lcor(
    adj.list, weight.list, y, z,
    type = "unit",
    y.diff.type = "difference",
    z.diff.type = "difference",
    epsilon = 0,
    winsorize.quantile = 0.05  # 5% winsorization
  )
  
  # Expected values
  expected_coeffs <- c(1.0, 0.714286, 0.005128, -1.0)

  # Expected bounds
  expected_y_lower <- -98.0
  expected_y_upper <- 98.0
  expected_z_lower <- -1.0
  expected_z_upper <- 1.0
  
  # Expected all.delta.y (BEFORE winsorization)
  expected_all_delta_y <- c(1.0, -1.0, 98.0, -98.0, -97.0, 97.0)
  
  # Expected all.delta.z (BEFORE winsorization)
  expected_all_delta_z <- c(1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
  
  # Expected vertex.delta.y (AFTER winsorization - clipped)
  expected_vertex_delta_y <- list(
    c(1.0),
    c(-1.0, 98.0),
    c(-98.0, -97.0),  # ✓ -97 is NOT clipped (within bounds [-98, 98])
    c(97.0)
 )

  expected_vertex_delta_z <- list(
    c(1.0),
    c(-1.0, 1.0),
    c(-1.0, 1.0),
    c(-1.0)
  )
  
  # Verify coefficients
  coeffs_ok <- TRUE
  for (i in 1:4) {
    if (!approx_equal(result$vertex.coefficients[i], expected_coeffs[i], 1e-5)) {
      print_fail(sprintf("vertex.coefficients[%d] = %.6f, expected %.6f",
                        i, result$vertex.coefficients[i], expected_coeffs[i]))
      coeffs_ok <- FALSE
    }
  }
  if (coeffs_ok) {
    cat("\u2713 vertex.coefficients correct\n")
  }
  
  # Verify winsorization bounds
  bounds_ok <- TRUE
  if (!approx_equal(result$y.lower, expected_y_lower)) {
    print_fail(sprintf("y.lower = %.6f, expected %.6f",
                      result$y.lower, expected_y_lower))
    bounds_ok <- FALSE
  }
  if (!approx_equal(result$y.upper, expected_y_upper)) {
    print_fail(sprintf("y.upper = %.6f, expected %.6f",
                      result$y.upper, expected_y_upper))
    bounds_ok <- FALSE
  }
  if (!approx_equal(result$z.lower, expected_z_lower)) {
    print_fail(sprintf("z.lower = %.6f, expected %.6f",
                      result$z.lower, expected_z_lower))
    bounds_ok <- FALSE
  }
  if (!approx_equal(result$z.upper, expected_z_upper)) {
    print_fail(sprintf("z.upper = %.6f, expected %.6f",
                      result$z.upper, expected_z_upper))
    bounds_ok <- FALSE
  }
  if (bounds_ok) {
    cat("\u2713 Winsorization bounds correct\n")
    cat(sprintf("  y.lower = %.1f, y.upper = %.1f\n", 
               result$y.lower, result$y.upper))
    cat(sprintf("  z.lower = %.1f, z.upper = %.1f\n",
               result$z.lower, result$z.upper))
  }
  
  # Verify all.delta.y (pre-winsorization)
  all_delta_y_ok <- TRUE
  if (length(result$all.delta.y) != length(expected_all_delta_y)) {
    print_fail("all.delta.y size mismatch")
    all_delta_y_ok <- FALSE
  } else {
    # Sort both for comparison (order might differ)
    sorted_result <- sort(result$all.delta.y)
    sorted_expected <- sort(expected_all_delta_y)
    
    if (!vectors_equal(sorted_result, sorted_expected)) {
      print_fail("all.delta.y values mismatch")
      all_delta_y_ok <- FALSE
    }
  }
  if (all_delta_y_ok) {
    cat("\u2713 all.delta.y correct (pre-winsorization)\n")
  }
  
  # Verify all.delta.z (pre-winsorization)
  all_delta_z_ok <- TRUE
  if (length(result$all.delta.z) != length(expected_all_delta_z)) {
    print_fail("all.delta.z size mismatch")
    all_delta_z_ok <- FALSE
  } else {
    sorted_result <- sort(result$all.delta.z)
    sorted_expected <- sort(expected_all_delta_z)
    
    if (!vectors_equal(sorted_result, sorted_expected)) {
      print_fail("all.delta.z values mismatch")
      all_delta_z_ok <- FALSE
    }
  }
  if (all_delta_z_ok) {
    cat("\u2713 all.delta.z correct (pre-winsorization)\n")
  }
  
  # Verify vertex.delta.y (post-winsorization, clipped values)
  vertex_delta_y_ok <- TRUE
  for (v in 1:4) {
    if (!vectors_equal(result$vertex.delta.y[[v]], expected_vertex_delta_y[[v]])) {
      print_fail(sprintf("vertex.delta.y[[%d]] mismatch (post-winsorization)", v))
      vertex_delta_y_ok <- FALSE
    }
  }
  if (vertex_delta_y_ok) {
    cat("\u2713 vertex.delta.y correct (post-winsorization)\n")
  }
  
  # Verify vertex.delta.z (should be unchanged by winsorization)
  vertex_delta_z_ok <- TRUE
  for (v in 1:4) {
    if (!vectors_equal(result$vertex.delta.z[[v]], expected_vertex_delta_z[[v]])) {
      print_fail(sprintf("vertex.delta.z[[%d]] mismatch", v))
      vertex_delta_z_ok <- FALSE
    }
  }
  if (vertex_delta_z_ok) {
    cat("\u2713 vertex.delta.z correct\n")
  }
  
  all_pass <- coeffs_ok && bounds_ok && all_delta_y_ok && 
              all_delta_z_ok && vertex_delta_y_ok && vertex_delta_z_ok
  
  if (all_pass) {
    print_pass()
  }
  
  return(all_pass)
}

# =============================================================================
# TEST 4: Negative correlation test
# =============================================================================
#
# Graph: Triangle (1 -- 2 -- 3 -- 1)
#
# Functions: y = c(1, 2, 3) (increasing)
#            z = c(3, 2, 1) (decreasing)
# Type: unit
# Diff types: both difference
#
# Expected: Negative correlations (functions change in opposite directions)
#
# Vertex 0 (R index 1):
#   Neighbors: [1, 2] (R: [2, 3])
#   Edge [0,1]: Δy = 2-1 = 1, Δz = 2-3 = -1
#   Edge [0,2]: Δy = 3-1 = 2, Δz = 1-3 = -2
#   numerator = 1*1*(-1) + 1*2*(-2) = -1 - 4 = -5
#   sum_y² = 1*1 + 1*4 = 5
#   sum_z² = 1*1 + 1*4 = 5
#   coeff[0] = -5 / (√5 * √5) = -5/5 = -1.0
#
# By symmetry, all vertices should have coeff = -1.0
# =============================================================================

test_negative_correlation <- function() {
  print_test_header("Negative correlation (inverse pattern)")
  
  # Build triangle graph
  adj.list <- list(
    c(2, 3),
    c(1, 3),
    c(1, 2)
  )
  
  weight.list <- list(
    c(1.0, 1.0),
    c(1.0, 1.0),
    c(1.0, 1.0)
  )
  
  y <- c(1.0, 2.0, 3.0)  # Increasing
  z <- c(3.0, 2.0, 1.0)  # Decreasing (inverse)
  
  result <- lcor(
    adj.list, weight.list, y, z,
    type = "unit",
    y.diff.type = "difference",
    z.diff.type = "difference",
    epsilon = 0,
    winsorize.quantile = 0
  )
  
  # All vertices should have correlation = -1.0
  expected_coeffs <- c(-1.0, -1.0, -1.0)
  
  coeffs_ok <- TRUE
  for (i in 1:3) {
    if (!approx_equal(result$vertex.coefficients[i], expected_coeffs[i])) {
      print_fail(sprintf("vertex.coefficients[%d] = %.6f, expected -1.0",
                        i, result$vertex.coefficients[i]))
      coeffs_ok <- FALSE
    }
  }
  
  if (coeffs_ok) {
    cat("\u2713 vertex.coefficients all equal to -1.0 (perfect negative correlation)\n")
    print_pass()
  }
  
  return(coeffs_ok)
}

# =============================================================================
# Main test runner
# =============================================================================

run_all_tests <- function() {
  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("  LCOR VALIDATION TEST SUITE (R Version)\n")
  cat("  Testing all output fields against hand-calculated values\n")
  cat(rep("=", 70), "\n", sep = "")
  
  passed <- 0
  total <- 4
  
  tests <- list(
    test_simple_path_unit_weights,
    test_derivative_weights_logratio,
    test_winsorization,
    test_negative_correlation
  )
  
  for (test_fn in tests) {
    if (test_fn()) {
      passed <- passed + 1
    }
  }
  
  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("  TEST SUMMARY\n")
  cat(rep("=", 70), "\n", sep = "")
  cat(sprintf("\nPassed: %d/%d tests\n", passed, total))
  
  if (passed == total) {
    cat("\n\u2713\u2713\u2713 ALL TESTS PASSED \u2713\u2713\u2713\n")
    invisible(TRUE)
  } else {
    cat("\n\u2717\u2717\u2717 SOME TESTS FAILED \u2717\u2717\u2717\n")
    invisible(FALSE)
  }
}

# Run tests if sourced directly
##if (sys.nframe() == 0) {
  run_all_tests()
##}

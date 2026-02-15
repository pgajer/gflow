# =============================================================================
# LSLOPE VALIDATION TEST SUITE (R Version)
# =============================================================================
#
# Comprehensive validation of lslope.gradient() and lslope.neighborhood()
# implementations by comparing all output fields against hand-calculated
# expected values.
#
# Run with: source("tests/manual/test_lslope_validation.R")
# =============================================================================

## Assumes package is loaded with devtools::load_all()

# Helper functions ============================================================

TOLERANCE <- 1e-10

approx_equal <- function(a, b, tol = TOLERANCE) {
    abs(a - b) < tol
}

vectors_equal <- function(a, b, tol = TOLERANCE) {
    if (length(a) != length(b)) return(FALSE)
    all(abs(a - b) < tol)
}

print_test_header <- function(test_name) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("TEST:", test_name, "\n")
    cat(rep("=", 70), "\n", sep = "")
}

print_pass <- function() {
    cat("\u2713 PASS\n", sep = "")
}

print_fail <- function(message) {
    cat("\u2717 FAIL:", message, "\n", sep = "")
}

# =============================================================================
# TEST 1: Simple path with gradient slope (ascending)
# =============================================================================
# Graph: 1 -- 2 -- 3
# y = c(1, 2, 3), z = c(2, 5, 6)
# Ascending gradient: vertex 1 -> 2, vertex 2 -> 3, vertex 3 is local max
# Expected slopes: Δz/Δy = 3/1=3.0, 1/1=1.0, 0 (extremum)
# =============================================================================

test_gradient_slope_ascending <- function() {
    print_test_header("Gradient slope (ascending), type='slope'")

    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(2.0, 5.0, 6.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        y.diff.type = "difference",
        z.diff.type = "difference",
        ascending = TRUE,
        instrumented = TRUE
    )

    expected_coeffs <- c(3.0, 1.0, 0.0)
    expected_neighbors <- c(2L, 3L, NA)
    expected_is_extremum <- c(FALSE, FALSE, TRUE)

    coeffs_ok <- vectors_equal(result$vertex.coefficients, expected_coeffs)
    if (coeffs_ok) cat("\u2713 vertex.coefficients correct: c(3.0, 1.0, 0.0)\n")
    else print_fail("vertex.coefficients mismatch")

    neighbors_ok <- is.na(result$gradient.neighbors[3]) &&
                    result$gradient.neighbors[1] == 2L &&
                    result$gradient.neighbors[2] == 3L
    if (neighbors_ok) cat("\u2713 gradient.neighbors correct: c(2, 3, NA)\n")
    else print_fail("gradient.neighbors mismatch")

    extremum_ok <- all(result$is.local.extremum == expected_is_extremum)
    if (extremum_ok) cat("\u2713 is.local.extremum correct\n")
    else print_fail("is.local.extremum mismatch")

    counts_ok <- result$n.local.maxima == 1 && result$n.local.minima == 0
    if (counts_ok) cat("\u2713 n.local.maxima=1, n.local.minima=0\n")
    else print_fail("extrema counts mismatch")

    all_pass <- coeffs_ok && neighbors_ok && extremum_ok && counts_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 2: Gradient slope (descending) - local minimum detection
# =============================================================================
# Same graph, ascending = FALSE
# Descending gradient: vertex 3 -> 2, vertex 2 -> 1, vertex 1 is local min
# Expected slopes: 0 (extremum), 3.0, 1.0
# =============================================================================

test_gradient_slope_descending <- function() {
    print_test_header("Gradient slope (descending), type='slope'")

    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(2.0, 5.0, 6.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        y.diff.type = "difference",
        z.diff.type = "difference",
        ascending = FALSE,
        instrumented = TRUE
    )

    expected_coeffs <- c(0.0, 3.0, 1.0)
    expected_is_extremum <- c(TRUE, FALSE, FALSE)

    coeffs_ok <- vectors_equal(result$vertex.coefficients, expected_coeffs)
    if (coeffs_ok) cat("\u2713 vertex.coefficients correct: c(0.0, 3.0, 1.0)\n")
    else print_fail("vertex.coefficients mismatch")

    extremum_ok <- all(result$is.local.extremum == expected_is_extremum)
    if (extremum_ok) cat("\u2713 Vertex 1 is local minimum\n")
    else print_fail("is.local.extremum mismatch")

    counts_ok <- result$n.local.maxima == 0 && result$n.local.minima == 1
    if (counts_ok) cat("\u2713 n.local.maxima=0, n.local.minima=1\n")
    else print_fail("extrema counts mismatch")

    all_pass <- coeffs_ok && extremum_ok && counts_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 3: Gradient sign (type = "sign")
# =============================================================================
# Graph: 1 -- 2 -- 3 -- 4
# y = c(1, 2, 3, 4) (monotonic), z = c(5, 3, 8, 2) (non-monotonic)
# Gradient signs: sign(3-5)=-1, sign(8-3)=+1, sign(2-8)=-1, 0 (extremum)
# =============================================================================

test_gradient_sign <- function() {
    print_test_header("Gradient sign, type='sign'")

    adj.list <- list(c(2L), c(1L, 3L), c(2L, 4L), c(3L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0, 1.0), c(1.0))
    y <- c(1.0, 2.0, 3.0, 4.0)
    z <- c(5.0, 3.0, 8.0, 2.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "sign",
        y.diff.type = "difference",
        z.diff.type = "difference",
        ascending = TRUE,
        instrumented = TRUE
    )

    expected_coeffs <- c(-1.0, 1.0, -1.0, 0.0)

    coeffs_ok <- vectors_equal(result$vertex.coefficients, expected_coeffs)
    if (coeffs_ok) cat("\u2713 vertex.coefficients correct: c(-1, +1, -1, 0)\n")
    else print_fail(sprintf("Got c(%.0f, %.0f, %.0f, %.0f)",
                            result$vertex.coefficients[1], result$vertex.coefficients[2],
                            result$vertex.coefficients[3], result$vertex.coefficients[4]))

    extremum_ok <- result$is.local.extremum[4] == TRUE &&
                   all(result$is.local.extremum[1:3] == FALSE)
    if (extremum_ok) cat("\u2713 Vertex 4 is local maximum\n")
    else print_fail("Extremum detection mismatch")

    all_pass <- coeffs_ok && extremum_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 4: Gradient slope normalized (bounded to [-1, 1])
# =============================================================================
# Graph: Triangle, y = c(1, 3, 2), z = c(0, 10, 4)
# Vertex 2 is local max (ascending), raw slopes at v1 = 5.0, v3 = 6.0
# Normalized with tanh should be bounded
# =============================================================================

test_gradient_slope_normalized <- function() {
    print_test_header("Gradient slope normalized, type='normalized'")

    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
    y <- c(1.0, 3.0, 2.0)
    z <- c(0.0, 10.0, 4.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "normalized",
        y.diff.type = "difference",
        z.diff.type = "difference",
        sigmoid.alpha = 0,
        ascending = TRUE,
        instrumented = TRUE
    )

    bounded_ok <- all(result$vertex.coefficients >= -1 & result$vertex.coefficients <= 1)
    if (bounded_ok) cat("\u2713 All coefficients bounded in [-1, 1]\n")
    else print_fail("Coefficients not bounded")

    extremum_ok <- result$is.local.extremum[2] == TRUE
    if (extremum_ok) cat("\u2713 Vertex 2 is local maximum\n")
    else print_fail("Local extremum detection incorrect")

    extremum_coef_ok <- approx_equal(result$vertex.coefficients[2], 0.0)
    if (extremum_coef_ok) cat("\u2713 Extremum coefficient is 0\n")
    else print_fail("Extremum coefficient should be 0")

    positive_ok <- result$vertex.coefficients[1] > 0 && result$vertex.coefficients[3] > 0
    if (positive_ok) cat("\u2713 Non-extremum coefficients are positive\n")
    else print_fail("Non-extremum coefficients should be positive")

    alpha_ok <- result$sigmoid.alpha > 0
    if (alpha_ok) cat(sprintf("\u2713 sigmoid.alpha = %.6f\n", result$sigmoid.alpha))
    else print_fail("sigmoid.alpha not computed")

    all_pass <- bounded_ok && extremum_ok && extremum_coef_ok && positive_ok && alpha_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 5: Neighborhood local regression coefficient (unit weights)
# =============================================================================
# Triangle graph, y = c(1, 2, 3), z = c(2, 4, 6)
# Perfect positive linear relationship: z = 2*y
# Expected β = 2.0 everywhere, lcor = 1.0 everywhere
# =============================================================================

test_neighborhood_slope_unit <- function() {
    print_test_header("Neighborhood regression coefficient, unit weights")

    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(2.0, 4.0, 6.0)

    result <- lslope.neighborhood(
        adj.list, weight.list, y, z,
        weight.type = "unit",
        y.diff.type = "difference",
        z.diff.type = "difference"
    )

    expected_coeffs <- c(2.0, 2.0, 2.0)
    expected_lcor <- c(1.0, 1.0, 1.0)

    coeffs_ok <- vectors_equal(result$vertex.coefficients, expected_coeffs)
    if (coeffs_ok) cat("\u2713 vertex.coefficients correct: c(2.0, 2.0, 2.0)\n")
    else print_fail(sprintf("Got c(%.4f, %.4f, %.4f)",
                            result$vertex.coefficients[1], result$vertex.coefficients[2],
                            result$vertex.coefficients[3]))

    lcor_ok <- vectors_equal(result$lcor, expected_lcor)
    if (lcor_ok) cat("\u2713 lcor correct: c(1.0, 1.0, 1.0)\n")
    else print_fail("lcor mismatch")

    ## Verify β = lcor × sd_z / sd_y
    relation_ok <- TRUE
    for (i in 1:3) {
        if (result$sd.y[i] > 1e-10) {
            expected_beta <- result$lcor[i] * result$sd.z[i] / result$sd.y[i]
            if (!approx_equal(result$vertex.coefficients[i], expected_beta, tol = 1e-6)) {
                relation_ok <- FALSE
            }
        }
    }
    if (relation_ok) cat("\u2713 Relationship beta = lcor * sd_z/sd_y verified\n")
    else print_fail("Relationship not satisfied")

    all_pass <- coeffs_ok && lcor_ok && relation_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 6: Neighborhood slope with derivative weights
# =============================================================================
# Path 1 -- 2 -- 3, edge lengths 1.0 and 2.0
# y = c(1, 2, 4), z = c(0, 2, 10)
# Derivative weights: w = 1/length^2
# =============================================================================

test_neighborhood_slope_derivative <- function() {
    print_test_header("Neighborhood regression coefficient, derivative weights")

    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 2.0), c(2.0))
    y <- c(1.0, 2.0, 4.0)
    z <- c(0.0, 2.0, 10.0)

    result <- lslope.neighborhood(
        adj.list, weight.list, y, z,
        weight.type = "derivative",
        y.diff.type = "difference",
        z.diff.type = "difference"
    )

    ## Manual calculations:
    ## Vertex 1: w=1, Dy=1, Dz=2, beta = 1*1*2 / (1*1) = 2.0
    ## Vertex 2: w1=1, Dy1=-1, Dz1=-2; w2=0.25, Dy2=2, Dz2=8
    ##           numer = 1*1*2 + 0.25*2*8 = 2+4 = 6
    ##           denom = 1*1 + 0.25*4 = 1+1 = 2
    ##           beta = 6/2 = 3.0
    ## Vertex 3: w=0.25, Dy=-2, Dz=-8, beta = 0.25*4*8 / (0.25*4) = 8/1 = 4.0
    ## Wait, let me recalculate...
    ## Vertex 3: Dy = 2-4 = -2, Dz = 2-10 = -8
    ##           numer = 0.25 * (-2) * (-8) = 4
    ##           denom = 0.25 * 4 = 1
    ##           beta = 4/1 = 4.0

    expected_coeffs <- c(2.0, 3.0, 4.0)

    coeffs_ok <- vectors_equal(result$vertex.coefficients, expected_coeffs)
    if (coeffs_ok) {
        cat("\u2713 vertex.coefficients correct: c(2.0, 3.0, 4.0)\n")
    } else {
        print_fail(sprintf("Got c(%.4f, %.4f, %.4f)",
                           result$vertex.coefficients[1], result$vertex.coefficients[2],
                           result$vertex.coefficients[3]))
    }

    if (coeffs_ok) print_pass()
    return(coeffs_ok)
}

# =============================================================================
# TEST 7: Negative slopes (inverse relationship)
# =============================================================================
# Triangle, y increasing, z decreasing
# Expected: negative slopes everywhere
# =============================================================================

test_negative_slopes <- function() {
    print_test_header("Negative slopes (inverse relationship)")

    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(6.0, 4.0, 2.0)

    ## Gradient slope (ascending)
    result_grad <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        ascending = TRUE,
        instrumented = TRUE
    )

    ## Neighborhood slope
    result_nbhd <- lslope.neighborhood(
        adj.list, weight.list, y, z,
        weight.type = "unit"
    )

    ## z = -2*y + 8, so slope should be -2 everywhere
    expected_nbhd <- c(-2.0, -2.0, -2.0)

    ## For gradient slope:
    ## v1: grad neighbor = 2 or 3, Dy to 2 = 1, Dz to 2 = -2, slope = -2
    ## v2: grad neighbor = 3, Dy = 1, Dz = -2, slope = -2
    ## v3: local max, slope = 0

    grad_ok <- all(result_grad$vertex.coefficients[1:2] < 0) &&
               approx_equal(result_grad$vertex.coefficients[3], 0.0)
    if (grad_ok) cat("\u2713 Gradient slopes are negative (except extremum)\n")
    else print_fail("Gradient slopes not negative")

    nbhd_ok <- vectors_equal(result_nbhd$vertex.coefficients, expected_nbhd)
    if (nbhd_ok) cat("\u2713 Neighborhood slopes correct: c(-2.0, -2.0, -2.0)\n")
    else print_fail(sprintf("Got c(%.4f, %.4f, %.4f)",
                            result_nbhd$vertex.coefficients[1],
                            result_nbhd$vertex.coefficients[2],
                            result_nbhd$vertex.coefficients[3]))

    lcor_ok <- all(result_nbhd$lcor < -0.99)
    if (lcor_ok) cat("\u2713 Local correlations are -1.0\n")
    else print_fail("lcor should be -1.0")

    all_pass <- grad_ok && nbhd_ok && lcor_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 8: Log-ratio differences for compositional data
# =============================================================================
# Path graph, y continuous, z compositional
# z = c(0.01, 0.10, 1.00) representing 10-fold increases
# =============================================================================

test_logratio_differences <- function() {
    print_test_header("Log-ratio differences for compositional data")

    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(0.01, 0.10, 1.00)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        y.diff.type = "difference",
        z.diff.type = "logratio",
        ascending = TRUE,
        instrumented = TRUE
    )

    ## For log-ratios: Δz = log(z[neighbor]/z[v])
    ## v1 -> v2: Dy = 1, Dz_log = log(0.10/0.01) = log(10) ≈ 2.303
    ##           slope = 2.303 / 1 = 2.303
    ## v2 -> v3: Dy = 1, Dz_log = log(1.00/0.10) = log(10) ≈ 2.303
    ##           slope = 2.303 / 1 = 2.303
    ## v3: local max, slope = 0

    expected_slope <- log(10)  # ≈ 2.302585

    slope_v1_ok <- approx_equal(result$vertex.coefficients[1], expected_slope, tol = 0.01)
    if (slope_v1_ok) cat(sprintf("\u2713 Vertex 1 slope = %.4f (log(10))\n", result$vertex.coefficients[1]))
    else print_fail(sprintf("Vertex 1 slope = %.4f, expected %.4f", result$vertex.coefficients[1], expected_slope))

    slope_v2_ok <- approx_equal(result$vertex.coefficients[2], expected_slope, tol = 0.01)
    if (slope_v2_ok) cat(sprintf("\u2713 Vertex 2 slope = %.4f (log(10))\n", result$vertex.coefficients[2]))
    else print_fail(sprintf("Vertex 2 slope = %.4f", result$vertex.coefficients[2]))

    extremum_ok <- result$is.local.extremum[3] == TRUE
    if (extremum_ok) cat("\u2713 Vertex 3 is local maximum\n")
    else print_fail("Vertex 3 should be local maximum")

    all_pass <- slope_v1_ok && slope_v2_ok && extremum_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 9: Multiple local extrema (star graph)
# =============================================================================
# Star graph: vertex 1 connected to 2, 3, 4, 5
# y = c(5, 1, 2, 3, 4) - vertex 1 is local max
# =============================================================================

test_star_graph_extrema <- function() {
    print_test_header("Star graph with central local maximum")

    ## Star graph: center vertex 1 connected to leaves 2,3,4,5
    adj.list <- list(
        c(2L, 3L, 4L, 5L),  # vertex 1 (center)
        c(1L),              # vertex 2
        c(1L),              # vertex 3
        c(1L),              # vertex 4
        c(1L)               # vertex 5
    )
    weight.list <- list(
        c(1.0, 1.0, 1.0, 1.0),
        c(1.0),
        c(1.0),
        c(1.0),
        c(1.0)
    )

    y <- c(5.0, 1.0, 2.0, 3.0, 4.0)  # center is max
    z <- c(10.0, 2.0, 4.0, 6.0, 8.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        ascending = TRUE,
        instrumented = TRUE
    )

    ## Vertex 1: all neighbors have lower y, so it's local max
    ## Vertices 2-5: only neighbor is vertex 1 which has higher y, so gradient -> 1

    extremum_v1_ok <- result$is.local.extremum[1] == TRUE
    if (extremum_v1_ok) cat("\u2713 Vertex 1 (center) is local maximum\n")
    else print_fail("Vertex 1 should be local maximum")

    leaves_not_extremum <- all(result$is.local.extremum[2:5] == FALSE)
    if (leaves_not_extremum) cat("\u2713 Leaf vertices are not extrema\n")
    else print_fail("Leaf vertices should not be extrema")

    ## All leaves should point to center
    all_point_to_center <- all(result$gradient.neighbors[2:5] == 1L)
    if (all_point_to_center) cat("\u2713 All leaf gradient edges point to center\n")
    else print_fail("Leaf gradient edges should point to center")

    ## Check slopes for leaves
    ## v2: Dy = 5-1 = 4, Dz = 10-2 = 8, slope = 8/4 = 2
    ## v3: Dy = 5-2 = 3, Dz = 10-4 = 6, slope = 6/3 = 2
    ## etc.
    expected_slopes <- c(0.0, 2.0, 2.0, 2.0, 2.0)
    slopes_ok <- vectors_equal(result$vertex.coefficients, expected_slopes)
    if (slopes_ok) cat("\u2713 All slopes correct: c(0, 2, 2, 2, 2)\n")
    else print_fail(sprintf("Got: c(%.2f, %.2f, %.2f, %.2f, %.2f)",
                            result$vertex.coefficients[1], result$vertex.coefficients[2],
                            result$vertex.coefficients[3], result$vertex.coefficients[4],
                            result$vertex.coefficients[5]))

    all_pass <- extremum_v1_ok && leaves_not_extremum && all_point_to_center && slopes_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 10: Verify gradient delta values
# =============================================================================
# Check that gradient.delta.y and gradient.delta.z are computed correctly
# =============================================================================

test_gradient_deltas <- function() {
    print_test_header("Gradient delta values verification")

    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
    y <- c(1.0, 5.0, 3.0)  # v2 is max
    z <- c(10.0, 20.0, 15.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        ascending = TRUE,
        instrumented = TRUE
    )

    ## Vertex 1: neighbors 2,3; y[2]=5 > y[3]=3, grad neighbor = 2
    ##           Dy = 5-1 = 4, Dz = 20-10 = 10
    ## Vertex 2: local max (both neighbors have lower y)
    ##           Dy = 0, Dz = 0
    ## Vertex 3: neighbors 1,2; y[2]=5 > y[1]=1, grad neighbor = 2
    ##           Dy = 5-3 = 2, Dz = 20-15 = 5

    expected_delta_y <- c(4.0, 0.0, 2.0)
    expected_delta_z <- c(10.0, 0.0, 5.0)

    delta_y_ok <- vectors_equal(result$gradient.delta.y, expected_delta_y)
    if (delta_y_ok) cat("\u2713 gradient.delta.y correct: c(4, 0, 2)\n")
    else print_fail(sprintf("Got: c(%.1f, %.1f, %.1f)",
                            result$gradient.delta.y[1], result$gradient.delta.y[2],
                            result$gradient.delta.y[3]))

    delta_z_ok <- vectors_equal(result$gradient.delta.z, expected_delta_z)
    if (delta_z_ok) cat("\u2713 gradient.delta.z correct: c(10, 0, 5)\n")
    else print_fail(sprintf("Got: c(%.1f, %.1f, %.1f)",
                            result$gradient.delta.z[1], result$gradient.delta.z[2],
                            result$gradient.delta.z[3]))

    ## Verify slopes = delta_z / delta_y
    expected_slopes <- c(10.0/4.0, 0.0, 5.0/2.0)  # c(2.5, 0, 2.5)
    slopes_ok <- vectors_equal(result$vertex.coefficients, expected_slopes)
    if (slopes_ok) cat("\u2713 Slopes = delta_z/delta_y verified\n")
    else print_fail("Slope calculation mismatch")

    all_pass <- delta_y_ok && delta_z_ok && slopes_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 11: Production version (non-instrumented) returns vector
# =============================================================================

test_production_version <- function() {
    print_test_header("Production version (instrumented=FALSE)")

    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0))
    y <- c(1.0, 2.0, 3.0)
    z <- c(2.0, 4.0, 6.0)

    result <- lslope.gradient(
        adj.list, weight.list, y, z,
        type = "slope",
        ascending = TRUE,
        instrumented = FALSE
    )

    ## Should return a numeric vector, not a list
    is_vector <- is.numeric(result) && !is.list(result)
    if (is_vector) cat("\u2713 Returns numeric vector (not list)\n")
    else print_fail("Should return numeric vector")

    length_ok <- length(result) == 3
    if (length_ok) cat("\u2713 Length equals number of vertices\n")
    else print_fail("Length mismatch")

    ## Check values
    expected <- c(2.0, 2.0, 0.0)  # z = 2*y, slope = 2 (except extremum)
    values_ok <- vectors_equal(result, expected)
    if (values_ok) cat("\u2713 Values correct: c(2.0, 2.0, 0.0)\n")
    else print_fail(sprintf("Got: c(%.2f, %.2f, %.2f)", result[1], result[2], result[3]))

    all_pass <- is_vector && length_ok && values_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# Main test runner
# =============================================================================

run_all_lslope_tests <- function() {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("  LSLOPE VALIDATION TEST SUITE (R Version)\n")
    cat("  Testing gradient-restricted and neighborhood local slopes\n")
    cat(rep("=", 70), "\n", sep = "")

    tests <- list(
        test_gradient_slope_ascending,
        test_gradient_slope_descending,
        test_gradient_sign,
        test_gradient_slope_normalized,
        test_neighborhood_slope_unit,
        test_neighborhood_slope_derivative,
        test_negative_slopes,
        test_logratio_differences,
        test_star_graph_extrema,
        test_gradient_deltas,
        test_production_version
    )

    passed <- 0
    total <- length(tests)

    for (test_fn in tests) {
        tryCatch({
            if (test_fn()) {
                passed <- passed + 1
            }
        }, error = function(e) {
            cat("\u2717 ERROR:", conditionMessage(e), "\n")
        })
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

# Run tests
run_all_lslope_tests()

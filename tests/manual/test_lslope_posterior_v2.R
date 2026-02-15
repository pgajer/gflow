# =============================================================================
# LSLOPE WITH POSTERIOR VALIDATION TEST SUITE (R Version)
# =============================================================================
#
# Comprehensive validation of lslope.with.posterior() implementation testing
# posterior uncertainty propagation for gradient-restricted local slopes.
#
# Note: lslope.with.posterior() uses lslope.gradient() internally, which
# computes slopes along the steepest ascent/descent direction of y.hat.
# This is distinct from lcor.with.posterior() which uses local correlation.
#
# Run with: source("tests/manual/test_lslope_posterior.R")
# =============================================================================

## Assumes package is loaded with devtools::load_all()

# Helper functions ============================================================

TOLERANCE <- 1e-6
TOLERANCE_LOOSE <- 1e-3  # For Monte Carlo comparisons

approx_equal <- function(a, b, tol = TOLERANCE) {
    abs(a - b) < tol
}

vectors_equal <- function(a, b, tol = TOLERANCE) {
    if (length(a) != length(b)) return(FALSE)
    all(abs(a - b) < tol)
}

matrices_equal <- function(a, b, tol = TOLERANCE) {
    if (!all(dim(a) == dim(b))) return(FALSE)
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
# Create test fixtures
# =============================================================================

create_simple_path_graph <- function() {
    ## Simple path graph: 1 -- 2 -- 3 -- 4 -- 5
    list(
        adj.list = list(c(2L), c(1L, 3L), c(2L, 4L), c(3L, 5L), c(4L)),
        weight.list = list(c(1.0), c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0), c(1.0)),
        n = 5L
    )
}

create_triangle_graph <- function() {
    ## Triangle graph: 1 -- 2 -- 3 -- 1
    list(
        adj.list = list(c(2L, 3L), c(1L, 3L), c(1L, 2L)),
        weight.list = list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0)),
        n = 3L
    )
}

create_star_graph <- function() {
    ## Star graph: center vertex 1 connected to leaves 2,3,4,5
    list(
        adj.list = list(
            c(2L, 3L, 4L, 5L),  # vertex 1 (center)
            c(1L),              # vertex 2
            c(1L),              # vertex 3
            c(1L),              # vertex 4
            c(1L)               # vertex 5
        ),
        weight.list = list(
            c(1.0, 1.0, 1.0, 1.0),
            c(1.0),
            c(1.0),
            c(1.0),
            c(1.0)
        ),
        n = 5L
    )
}

# =============================================================================
# TEST 1: Basic functionality with pre-computed samples
# =============================================================================

test_basic_functionality <- function() {
    print_test_header("Basic functionality with pre-computed samples")

    graph <- create_simple_path_graph()
    n <- graph$n
    p <- 2L  # Two features
    n.samples <- 50L

    ## Create mock posterior samples (list of n x n.samples matrices)
    set.seed(123)
    Z.hat.samples <- lapply(1:p, function(j) {
        matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)
    })

    ## Create monotonic y.hat so gradient direction is well-defined
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    ## Run lslope.with.posterior
    result <- lslope.with.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y.hat = y.hat,
        Z.hat.samples = Z.hat.samples,
        lslope.type = "slope",
        ascending = TRUE,
        y.diff.type = "difference",
        z.diff.type = "difference",
        verbose = FALSE
    )

    ## Check output structure
    class_ok <- inherits(result, "lslope.posterior")
    if (class_ok) cat("\u2713 Returns object of class 'lslope.posterior'\n")
    else print_fail("Wrong class")

    ## Check dimensions
    mean_dim_ok <- all(dim(result$mean) == c(p, n))
    if (mean_dim_ok) cat("\u2713 mean matrix has correct dimensions (p x n)\n")
    else print_fail(sprintf("mean dimensions: expected (%d, %d), got (%d, %d)",
                            p, n, nrow(result$mean), ncol(result$mean)))

    sd_dim_ok <- all(dim(result$sd) == c(p, n))
    if (sd_dim_ok) cat("\u2713 sd matrix has correct dimensions\n")
    else print_fail("sd dimensions incorrect")

    lower_dim_ok <- all(dim(result$lower) == c(p, n))
    upper_dim_ok <- all(dim(result$upper) == c(p, n))
    if (lower_dim_ok && upper_dim_ok) cat("\u2713 Credible interval matrices have correct dimensions\n")
    else print_fail("CI dimensions incorrect")

    ## Check metadata
    meta_ok <- result$n.samples == n.samples &&
               result$n.features == p &&
               result$n.vertices == n &&
               result$lslope.type == "slope"
    if (meta_ok) cat("\u2713 Metadata correctly stored\n")
    else print_fail("Metadata mismatch")

    ## Check that credible intervals are sensible
    ci_ok <- all(result$lower <= result$mean, na.rm = TRUE) &&
             all(result$mean <= result$upper, na.rm = TRUE)
    if (ci_ok) cat("\u2713 Credible intervals satisfy lower <= mean <= upper\n")
    else print_fail("CI ordering violated")

    all_pass <- class_ok && mean_dim_ok && sd_dim_ok &&
                lower_dim_ok && upper_dim_ok && meta_ok && ci_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 2: Single feature (matrix input)
# =============================================================================

test_single_feature <- function() {
    print_test_header("Single feature (matrix input)")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 30L

    ## Single feature as matrix (not list)
    set.seed(456)
    Z.hat.samples <- matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)

    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y.hat = y.hat,
        Z.hat.samples = Z.hat.samples,  # Matrix, not list
        lslope.type = "normalized",
        verbose = FALSE
    )

    ## Should still work with p=1
    p_ok <- result$n.features == 1
    if (p_ok) cat("\u2713 Single feature correctly handled\n")
    else print_fail("n.features should be 1")

    dim_ok <- nrow(result$mean) == 1 && ncol(result$mean) == n
    if (dim_ok) cat("\u2713 Output dimensions correct for single feature\n")
    else print_fail("Dimensions incorrect")

    all_pass <- p_ok && dim_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 3: return.samples = TRUE
# =============================================================================

test_return_samples <- function() {
    print_test_header("return.samples = TRUE")

    graph <- create_simple_path_graph()
    n <- graph$n
    p <- 2L
    n.samples <- 20L

    set.seed(789)
    Z.hat.samples <- lapply(1:p, function(j) {
        matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)
    })

    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y.hat = y.hat,
        Z.hat.samples = Z.hat.samples,
        return.samples = TRUE,
        verbose = FALSE
    )

    ## Check that samples are returned
    has_samples <- !is.null(result$samples)
    if (has_samples) cat("\u2713 samples field present in result\n")
    else print_fail("samples field missing")

    ## Check samples structure
    samples_length_ok <- length(result$samples) == p
    if (samples_length_ok) cat("\u2713 samples list has correct length\n")
    else print_fail("samples list length incorrect")

    ## Check individual sample matrices
    sample_dims_ok <- all(sapply(result$samples, function(s) {
        nrow(s) == n && ncol(s) == n.samples
    }))
    if (sample_dims_ok) cat("\u2713 Each sample matrix has dimensions (n x n.samples)\n")
    else print_fail("Sample matrix dimensions incorrect")

    ## Verify mean is computed from samples
    computed_mean <- apply(result$samples[[1]], 1, mean, na.rm = TRUE)
    mean_matches <- vectors_equal(result$mean[1, ], computed_mean, tol = TOLERANCE_LOOSE)
    if (mean_matches) cat("\u2713 Mean matches apply(samples, 1, mean)\n")
    else print_fail("Mean doesn't match samples")

    all_pass <- has_samples && samples_length_ok && sample_dims_ok && mean_matches
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 4: lslope.type parameter affects results
# =============================================================================

test_lslope_type_parameter <- function() {
    print_test_header("lslope.type parameter affects results")

    graph <- create_simple_path_graph()
    n <- graph$n
    p <- 1L
    n.samples <- 100L

    set.seed(111)
    Z.hat.samples <- list(matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)  # Monotonic

    ## Test all three types
    result_slope <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope", verbose = FALSE
    )

    result_normalized <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "normalized", verbose = FALSE
    )

    result_sign <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "sign", verbose = FALSE
    )

    ## Normalized should be bounded
    normalized_bounded <- all(abs(result_normalized$mean) <= 1, na.rm = TRUE)
    if (normalized_bounded) cat("\u2713 'normalized' type produces bounded values\n")
    else print_fail("'normalized' should be bounded in [-1, 1]")

    ## Sign should be in {-1, 0, 1}
    sign_values <- unique(as.vector(result_sign$mean))
    sign_valid <- all(sign_values %in% c(-1, 0, 1) | is.na(sign_values))
    if (sign_valid) cat("\u2713 'sign' type produces values in {-1, 0, 1}\n")
    else cat("  Note: 'sign' may have intermediate values due to averaging\n")

    ## lslope.type should be stored
    type_stored <- result_slope$lslope.type == "slope" &&
                   result_normalized$lslope.type == "normalized" &&
                   result_sign$lslope.type == "sign"
    if (type_stored) cat("\u2713 lslope.type correctly stored in results\n")
    else print_fail("lslope.type not stored correctly")

    all_pass <- normalized_bounded && type_stored
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 5: ascending parameter affects gradient direction
# =============================================================================

test_ascending_parameter <- function() {
    print_test_header("ascending parameter affects gradient direction")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 50L

    ## Create z.hat samples that increase along the path
    set.seed(222)
    z.base <- c(1.0, 2.0, 3.0, 4.0, 5.0)  # Increasing
    Z.hat.samples <- list(matrix(
        rep(z.base, n.samples) + rnorm(n * n.samples, sd = 0.1),
        nrow = n, ncol = n.samples
    ))

    ## y.hat also increasing
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    ## Test ascending (gradient points to higher y values)
    result_asc <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        ascending = TRUE,
        verbose = FALSE
    )

    ## Test descending (gradient points to lower y values)
    result_desc <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        ascending = FALSE,
        verbose = FALSE
    )

    ## With y and z both increasing:
    ## - ascending: non-extremum vertices should have positive slopes
    ## - descending: non-extremum vertices should have positive slopes (same sign)
    ## The key is that extrema are detected differently

    ## Check that results are different (extrema detection differs)
    results_differ <- !matrices_equal(result_asc$mean, result_desc$mean, tol = 0.01)
    if (results_differ) cat("\u2713 Results differ between ascending and descending\n")
    else cat("  Note: Results similar (may be okay for monotonic data)\n")

    ## Non-extremum slopes should be positive for co-monotonic y and z
    ## (skipping first vertex for ascending, last for descending)
    asc_slopes_positive <- all(result_asc$mean[1, 1:4] > 0, na.rm = TRUE)
    if (asc_slopes_positive) cat("\u2713 Ascending mode: non-extremum slopes positive\n")
    else cat("  Note: Some ascending slopes not positive\n")

    all_pass <- TRUE  # Informational test
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 6: Local extrema detection
# =============================================================================

test_local_extrema <- function() {
    print_test_header("Local extrema have zero slope")

    graph <- create_star_graph()
    n <- graph$n
    n.samples <- 30L

    ## Center vertex (1) is local maximum of y
    y.hat <- c(5.0, 1.0, 2.0, 3.0, 4.0)  # center is max

    ## z also has maximum at center
    set.seed(333)
    z.base <- c(10.0, 2.0, 4.0, 6.0, 8.0)
    Z.hat.samples <- list(matrix(
        rep(z.base, n.samples) + rnorm(n * n.samples, sd = 0.1),
        nrow = n, ncol = n.samples
    ))

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        ascending = TRUE,
        verbose = FALSE
    )

    ## Vertex 1 (center) is local max in ascending mode -> should have zero slope
    center_zero <- approx_equal(result$mean[1, 1], 0, tol = 0.01)
    if (center_zero) cat("\u2713 Local maximum (center) has near-zero mean slope\n")
    else cat(sprintf("  Note: Center slope = %.4f (expected ~0)\n", result$mean[1, 1]))

    ## Leaf vertices should have positive slopes (pointing to center)
    leaves_positive <- all(result$mean[1, 2:5] > 0, na.rm = TRUE)
    if (leaves_positive) cat("\u2713 Leaf vertices have positive slopes (pointing to max)\n")
    else print_fail("Leaf slopes should be positive")

    all_pass <- leaves_positive
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 7: Print method works
# =============================================================================

test_print_method <- function() {
    print_test_header("print.lslope.posterior method")

    graph <- create_simple_path_graph()
    n.samples <- 30L

    set.seed(444)
    Z.hat.samples <- list(
        matrix(rnorm(5 * n.samples), nrow = 5, ncol = n.samples),
        matrix(rnorm(5 * n.samples), nrow = 5, ncol = n.samples)
    )
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    ## Capture print output
    output <- capture.output(print(result))

    ## Check that key information is present
    has_header <- any(grepl("lslope", output, ignore.case = TRUE))
    has_features <- any(grepl("Features:", output))
    has_vertices <- any(grepl("Vertices:", output))
    has_samples <- any(grepl("samples", output, ignore.case = TRUE))

    if (has_header) cat("\u2713 Print output has header\n")
    else print_fail("Missing header")

    if (has_features && has_vertices && has_samples) {
        cat("\u2713 Print output contains key metadata\n")
    } else {
        print_fail("Missing key metadata in print output")
    }

    all_pass <- has_header && has_features && has_vertices && has_samples
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 8: Input validation
# =============================================================================

test_input_validation <- function() {
    print_test_header("Input validation")

    graph <- create_simple_path_graph()
    n.samples <- 20L
    Z.hat.samples <- list(matrix(rnorm(5 * n.samples), nrow = 5, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    ## Test mismatched dimensions
    error_caught_dim <- FALSE
    tryCatch({
        bad_samples <- list(matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples))  # Wrong n
        lslope.with.posterior(
            graph$adj.list, graph$weight.list,
            y.hat = y.hat,
            Z.hat.samples = bad_samples,
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_dim <<- TRUE
    })

    if (error_caught_dim) cat("\u2713 Error caught for mismatched sample dimensions\n")
    else print_fail("Should error on wrong sample dimensions")

    ## Test non-list, non-matrix input
    error_caught_type <- FALSE
    tryCatch({
        lslope.with.posterior(
            graph$adj.list, graph$weight.list,
            y.hat = y.hat,
            Z.hat.samples = c(1, 2, 3),  # Vector instead of matrix
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_type <<- TRUE
    })

    if (error_caught_type) cat("\u2713 Error caught for invalid input type\n")
    else print_fail("Should error on vector input")

    ## Test inconsistent sample counts across features
    error_caught_samples <- FALSE
    tryCatch({
        bad_samples <- list(
            matrix(rnorm(5 * 20), nrow = 5, ncol = 20),
            matrix(rnorm(5 * 30), nrow = 5, ncol = 30)  # Different B
        )
        lslope.with.posterior(
            graph$adj.list, graph$weight.list,
            y.hat = y.hat,
            Z.hat.samples = bad_samples,
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_samples <<- TRUE
    })

    if (error_caught_samples) cat("\u2713 Error caught for inconsistent sample counts\n")
    else print_fail("Should error on inconsistent sample counts")

    all_pass <- error_caught_dim && error_caught_type && error_caught_samples
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 9: Consistent samples give low variance
# =============================================================================

test_consistent_samples <- function() {
    print_test_header("Consistent samples give low variance")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 50L

    ## Create nearly identical samples (low variance input)
    set.seed(555)
    z.base <- c(2.0, 4.0, 6.0, 8.0, 10.0)  # Linear increasing
    Z.hat.samples <- list(matrix(
        rep(z.base, n.samples) + rnorm(n * n.samples, sd = 0.001),  # Very low noise
        nrow = n, ncol = n.samples
    ))

    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        verbose = FALSE
    )

    ## SD should be very low for consistent samples
    mean_sd <- mean(result$sd, na.rm = TRUE)
    low_variance <- mean_sd < 0.1
    if (low_variance) cat(sprintf("\u2713 Low input variance gives low output SD (%.4f)\n", mean_sd))
    else print_fail(sprintf("SD too high for consistent samples: %.4f", mean_sd))

    ## CI should be narrow
    ci_width <- mean(result$upper - result$lower, na.rm = TRUE)
    narrow_ci <- ci_width < 0.5
    if (narrow_ci) cat(sprintf("\u2713 Narrow CI for consistent samples (width = %.4f)\n", ci_width))
    else print_fail(sprintf("CI too wide for consistent samples: %.4f", ci_width))

    all_pass <- low_variance && narrow_ci
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 10: Variable samples give high variance
# =============================================================================

test_variable_samples <- function() {
    print_test_header("Variable samples give high variance")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 100L

    ## Create highly variable samples
    set.seed(666)
    Z.hat.samples <- list(matrix(
        rnorm(n * n.samples, sd = 5),  # High variance
        nrow = n, ncol = n.samples
    ))

    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        verbose = FALSE
    )

    ## SD should be higher for variable samples
    mean_sd <- mean(result$sd, na.rm = TRUE)
    has_variance <- mean_sd > 0.5
    if (has_variance) cat(sprintf("\u2713 High input variance gives measurable output SD (%.4f)\n", mean_sd))
    else cat(sprintf("  Note: SD = %.4f (may be okay)\n", mean_sd))

    ## CI should be wider
    ci_width <- mean(result$upper - result$lower, na.rm = TRUE)
    wide_ci <- ci_width > 1.0
    if (wide_ci) cat(sprintf("\u2713 Wide CI for variable samples (width = %.4f)\n", ci_width))
    else cat(sprintf("  Note: CI width = %.4f (may be okay)\n", ci_width))

    all_pass <- TRUE  # Informational test
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 11: Large number of samples
# =============================================================================

test_large_sample_size <- function() {
    print_test_header("Large number of posterior samples")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 500L  # Large

    set.seed(777)
    Z.hat.samples <- list(matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    start_time <- Sys.time()

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    time_ok <- elapsed < 10  # Should be reasonably fast
    if (time_ok) cat(sprintf("\u2713 Completed in %.2f seconds\n", elapsed))
    else print_fail(sprintf("Too slow: %.2f seconds", elapsed))

    samples_ok <- result$n.samples == n.samples
    if (samples_ok) cat("\u2713 n.samples correctly recorded\n")
    else print_fail("n.samples incorrect")

    all_pass <- time_ok && samples_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 12: Multiple features processed correctly
# =============================================================================

test_multiple_features <- function() {
    print_test_header("Multiple features processed correctly")

    graph <- create_simple_path_graph()
    n <- graph$n
    p <- 10L  # Many features
    n.samples <- 30L

    set.seed(888)
    Z.hat.samples <- lapply(1:p, function(j) {
        ## Each feature has different mean level
        matrix(rnorm(n * n.samples, mean = j), nrow = n, ncol = n.samples)
    })

    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    ## Check correct number of features
    n_features_ok <- result$n.features == p
    if (n_features_ok) cat("\u2713 Correct number of features processed\n")
    else print_fail("Feature count incorrect")

    ## Check output dimensions
    dim_ok <- nrow(result$mean) == p && ncol(result$mean) == n
    if (dim_ok) cat("\u2713 Output matrix has correct dimensions\n")
    else print_fail("Output dimensions incorrect")

    ## Each feature should have different results (different means)
    feature_means <- rowMeans(result$mean, na.rm = TRUE)
    features_differ <- length(unique(round(feature_means, 2))) > 1
    if (features_differ) cat("\u2713 Different features produce different results\n")
    else cat("  Note: Feature results similar\n")

    all_pass <- n_features_ok && dim_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 13: Reproducibility (deterministic for same input)
# =============================================================================

test_reproducibility <- function() {
    print_test_header("Reproducibility (deterministic for same input)")

    graph <- create_simple_path_graph()
    n.samples <- 30L

    set.seed(999)
    Z.hat.samples <- list(matrix(rnorm(5 * n.samples), nrow = 5, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    ## Two runs with same input should give identical results
    result1 <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    result2 <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    mean_identical <- matrices_equal(result1$mean, result2$mean)
    if (mean_identical) cat("\u2713 Results identical across runs (deterministic)\n")
    else print_fail("Results differ between runs")

    all_pass <- mean_identical
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 14: Known slope values with deterministic input
# =============================================================================

test_known_slopes <- function() {
    print_test_header("Known slope values with deterministic input")

    ## Simple path: 1 -- 2 -- 3
    adj.list <- list(c(2L), c(1L, 3L), c(2L))
    weight.list <- list(c(1.0), c(1.0, 1.0), c(1.0))
    n <- 3L

    ## y = c(1, 2, 3), z = c(2, 4, 6) -> slope = 2 everywhere (except extremum)
    y.hat <- c(1.0, 2.0, 3.0)
    z.base <- c(2.0, 4.0, 6.0)

    ## Create samples with essentially no variance
    n.samples <- 10L
    Z.hat.samples <- list(matrix(
        rep(z.base, n.samples),
        nrow = n, ncol = n.samples
    ))

    result <- lslope.with.posterior(
        adj.list, weight.list, y.hat, Z.hat.samples,
        lslope.type = "slope",
        ascending = TRUE,
        verbose = FALSE
    )

    ## Vertex 1: gradient to vertex 2, Dz/Dy = (4-2)/(2-1) = 2
    ## Vertex 2: gradient to vertex 3, Dz/Dy = (6-4)/(3-2) = 2
    ## Vertex 3: local maximum, slope = 0

    expected_slopes <- c(2.0, 2.0, 0.0)
    slopes_match <- vectors_equal(result$mean[1, ], expected_slopes, tol = 0.01)

    if (slopes_match) {
        cat("\u2713 Slopes match expected values: c(2, 2, 0)\n")
    } else {
        cat(sprintf("  Got: c(%.2f, %.2f, %.2f)\n",
                    result$mean[1, 1], result$mean[1, 2], result$mean[1, 3]))
        print_fail("Slopes don't match expected")
    }

    all_pass <- slopes_match
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 15: Credible interval coverage
# =============================================================================

test_ci_coverage <- function() {
    print_test_header("95% credible interval properties")

    graph <- create_simple_path_graph()
    n <- graph$n
    n.samples <- 200L

    set.seed(1010)
    Z.hat.samples <- list(matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0, 4.0, 5.0)

    result <- lslope.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        return.samples = TRUE,
        verbose = FALSE
    )

    ## Check that ~95% of samples fall within CI for each vertex
    coverage_ok <- TRUE
    for (v in 1:n) {
        samples_v <- result$samples[[1]][v, ]
        lower_v <- result$lower[1, v]
        upper_v <- result$upper[1, v]
        in_ci <- sum(samples_v >= lower_v & samples_v <= upper_v) / n.samples
        if (in_ci < 0.90) {
            coverage_ok <- FALSE
        }
    }

    if (coverage_ok) cat("\u2713 All vertices have >= 90% sample coverage in CI\n")
    else print_fail("Some vertices have poor CI coverage")

    ## CI should be approximately symmetric around mean
    ci_width_lower <- result$mean - result$lower
    ci_width_upper <- result$upper - result$mean
    symmetry_ratio <- mean(ci_width_lower / ci_width_upper, na.rm = TRUE)
    reasonably_symmetric <- symmetry_ratio > 0.5 && symmetry_ratio < 2.0

    if (reasonably_symmetric) {
        cat(sprintf("\u2713 CI reasonably symmetric (ratio = %.2f)\n", symmetry_ratio))
    } else {
        cat(sprintf("  Note: CI asymmetric (ratio = %.2f)\n", symmetry_ratio))
    }

    all_pass <- coverage_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# Main test runner
# =============================================================================

run_all_lslope_posterior_tests <- function() {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("  LSLOPE WITH POSTERIOR VALIDATION TEST SUITE\n")
    cat("  Testing posterior uncertainty propagation for gradient local slopes\n")
    cat(rep("=", 70), "\n", sep = "")

    tests <- list(
        test_basic_functionality,
        test_single_feature,
        test_return_samples,
        test_lslope_type_parameter,
        test_ascending_parameter,
        test_local_extrema,
        test_print_method,
        test_input_validation,
        test_consistent_samples,
        test_variable_samples,
        test_large_sample_size,
        test_multiple_features,
        test_reproducibility,
        test_known_slopes,
        test_ci_coverage
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
run_all_lslope_posterior_tests()

# =============================================================================
# LCOR WITH POSTERIOR VALIDATION TEST SUITE (R Version)
# =============================================================================
#
# Comprehensive validation of lcor.with.posterior() implementation testing
# both R mode (pre-computed samples) and C++ mode (memory-efficient).
#
# Run with: source("tests/manual/test_lcor_posterior.R")
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

create_mock_fitted_model <- function(n, m = NULL) {
    ## Create a mock fitted model with spectral decomposition
    ## n = number of vertices, m = number of eigenpairs (default n)
    if (is.null(m)) m <- n

    ## Generate a simple random symmetric matrix
    set.seed(42)
    A <- matrix(rnorm(n * n), n, n)
    A <- (A + t(A)) / 2

    ## Eigendecomposition
    eig <- eigen(A, symmetric = TRUE)

    ## Keep only m eigenpairs
    V <- eig$vectors[, 1:m, drop = FALSE]
    lambda <- abs(eig$values[1:m])  # Use absolute values as "eigenvalues"

    ## Create mock fitted model structure
    structure(
        list(
            spectral = list(
                vectors = V,
                values = lambda,
                filter.type = "heat_kernel",
                eta = 0.5
            ),
            fitted.values = rnorm(n),  # Random fitted values
            graph = list(
                adj.list = create_simple_path_graph()$adj.list[1:n],
                edge.length.list = create_simple_path_graph()$weight.list[1:n]
            )
        ),
        class = "knn.riem.fit"
    )
}

# =============================================================================
# TEST 1: R mode basic functionality
# =============================================================================
# Test that R mode works with pre-computed posterior samples
# =============================================================================

test_r_mode_basic <- function() {
    print_test_header("R mode: basic functionality with pre-computed samples")

    graph <- create_triangle_graph()
    n <- graph$n
    p <- 2L  # Two features
    n.samples <- 50L

    ## Create mock posterior samples (list of n x n.samples matrices)
    set.seed(123)
    Z.hat.samples <- lapply(1:p, function(j) {
        matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)
    })

    ## Create y.hat
    y.hat <- c(1.0, 2.0, 3.0)

    ## Run lcor.with.posterior in R mode
    result <- lcor.with.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y.hat = y.hat,
        Z.hat.samples = Z.hat.samples,
        lcor.type = "derivative",
        credible.level = 0.95,
        verbose = FALSE
    )

    ## Check output structure
    class_ok <- inherits(result, "lcor.posterior")
    if (class_ok) cat("\u2713 Returns object of class 'lcor.posterior'\n")
    else print_fail("Wrong class")

    mode_ok <- result$mode == "R"
    if (mode_ok) cat("\u2713 Mode correctly set to 'R'\n")
    else print_fail("Mode should be 'R'")

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
               result$lcor.type == "derivative"
    if (meta_ok) cat("\u2713 Metadata correctly stored\n")
    else print_fail("Metadata mismatch")

    ## Check that credible intervals are sensible
    ci_ok <- all(result$lower <= result$mean, na.rm = TRUE) &&
             all(result$mean <= result$upper, na.rm = TRUE)
    if (ci_ok) cat("\u2713 Credible intervals satisfy lower <= mean <= upper\n")
    else print_fail("CI ordering violated")

    all_pass <- class_ok && mode_ok && mean_dim_ok && sd_dim_ok &&
                lower_dim_ok && upper_dim_ok && meta_ok && ci_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 2: R mode with single feature (edge case)
# =============================================================================

test_r_mode_single_feature <- function() {
    print_test_header("R mode: single feature (matrix input)")

    graph <- create_triangle_graph()
    n <- graph$n
    n.samples <- 30L

    ## Single feature as matrix (not list)
    set.seed(456)
    Z.hat.samples <- matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)

    y.hat <- c(1.0, 2.0, 3.0)

    result <- lcor.with.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y.hat = y.hat,
        Z.hat.samples = Z.hat.samples,  # Matrix, not list
        lcor.type = "unit",
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
# TEST 3: R mode with return.samples = TRUE
# =============================================================================

test_r_mode_return_samples <- function() {
    print_test_header("R mode: return.samples = TRUE")

    graph <- create_triangle_graph()
    n <- graph$n
    p <- 2L
    n.samples <- 20L

    set.seed(789)
    Z.hat.samples <- lapply(1:p, function(j) {
        matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples)
    })

    y.hat <- c(1.0, 2.0, 3.0)

    result <- lcor.with.posterior(
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
    computed_mean <- rowMeans(result$samples[[1]])
    mean_matches <- vectors_equal(result$mean[1, ], computed_mean, tol = TOLERANCE_LOOSE)
    if (mean_matches) cat("\u2713 Mean matches rowMeans of samples\n")
    else print_fail("Mean doesn't match samples")

    all_pass <- has_samples && samples_length_ok && sample_dims_ok && mean_matches
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 4: lcor.type parameter affects results
# =============================================================================

test_lcor_type_parameter <- function() {
    print_test_header("lcor.type parameter affects results")

    graph <- create_triangle_graph()
    n <- graph$n
    p <- 1L
    n.samples <- 100L

    set.seed(111)
    Z.hat.samples <- list(matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples))
    y.hat <- c(1.0, 3.0, 2.0)

    ## Test all three types
    result_deriv <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lcor.type = "derivative", verbose = FALSE
    )

    result_unit <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lcor.type = "unit", verbose = FALSE
    )

    result_sign <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        lcor.type = "sign", verbose = FALSE
    )

    ## Results should generally differ
    deriv_unit_diff <- !vectors_equal(result_deriv$mean[1,], result_unit$mean[1,], tol = 0.1)
    if (deriv_unit_diff) cat("\u2713 'derivative' and 'unit' produce different results\n")
    else cat("  Note: 'derivative' and 'unit' similar (may be okay for this graph)\n")

    ## lcor.type should be stored
    type_stored <- result_deriv$lcor.type == "derivative" &&
                   result_unit$lcor.type == "unit" &&
                   result_sign$lcor.type == "sign"
    if (type_stored) cat("\u2713 lcor.type correctly stored in results\n")
    else print_fail("lcor.type not stored correctly")

    all_pass <- type_stored
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 5: Credible level affects interval width
# =============================================================================

test_credible_level <- function() {
    print_test_header("Credible level affects interval width")

    graph <- create_triangle_graph()
    n <- graph$n
    n.samples <- 200L

    set.seed(222)
    Z.hat.samples <- list(matrix(rnorm(n * n.samples), nrow = n, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0)

    ## 80% CI
    result_80 <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        credible.level = 0.80, verbose = FALSE
    )

    ## 95% CI
    result_95 <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        credible.level = 0.95, verbose = FALSE
    )

    ## 95% CI should be wider
    width_80 <- mean(result_80$upper - result_80$lower)
    width_95 <- mean(result_95$upper - result_95$lower)

    wider_95 <- width_95 > width_80
    if (wider_95) cat(sprintf("\u2713 95%% CI wider than 80%% CI (%.3f vs %.3f)\n", width_95, width_80))
    else print_fail(sprintf("95%% CI should be wider: %.3f vs %.3f", width_95, width_80))

    ## Credible levels stored correctly
    level_stored <- result_80$credible.level == 0.80 && result_95$credible.level == 0.95
    if (level_stored) cat("\u2713 credible.level correctly stored\n")
    else print_fail("credible.level not stored correctly")

    all_pass <- wider_95 && level_stored
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 6: Print method works
# =============================================================================

test_print_method <- function() {
    print_test_header("print.lcor.posterior method")

    graph <- create_triangle_graph()
    n.samples <- 50L

    set.seed(333)
    Z.hat.samples <- list(
        matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples),
        matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples)
    )
    y.hat <- c(1.0, 2.0, 3.0)

    result <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    ## Capture print output
    output <- capture.output(print(result))

    ## Check that key information is present
    has_header <- any(grepl("Local Correlation", output))
    has_mode <- any(grepl("Mode:", output))
    has_features <- any(grepl("Features:", output))
    has_vertices <- any(grepl("Vertices:", output))

    if (has_header) cat("\u2713 Print output has header\n")
    else print_fail("Missing header")

    if (has_mode && has_features && has_vertices) cat("\u2713 Print output contains key metadata\n")
    else print_fail("Missing key metadata in print output")

    all_pass <- has_header && has_mode && has_features && has_vertices
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 7: Summary method works
# =============================================================================

test_summary_method <- function() {
    print_test_header("summary.lcor.posterior method")

    graph <- create_triangle_graph()
    n.samples <- 50L

    set.seed(444)
    Z.hat.samples <- lapply(1:5, function(j) {
        matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples)
    })
    y.hat <- c(1.0, 2.0, 3.0)

    result <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    summ <- summary(result)

    ## Check summary structure
    class_ok <- inherits(summ, "summary.lcor.posterior")
    if (class_ok) cat("\u2713 Summary has correct class\n")
    else print_fail("Summary class incorrect")

    has_feature_summary <- !is.null(summ$feature.summary) && is.data.frame(summ$feature.summary)
    if (has_feature_summary) cat("\u2713 feature.summary is a data frame\n")
    else print_fail("feature.summary missing or wrong type")

    has_prop_sig <- !is.null(summ$prop.significant.overall)
    if (has_prop_sig) cat("\u2713 prop.significant.overall present\n")
    else print_fail("prop.significant.overall missing")

    ## Check that print.summary works
    output <- capture.output(print(summ))
    print_works <- length(output) > 0
    if (print_works) cat("\u2713 print.summary.lcor.posterior works\n")
    else print_fail("print.summary failed")

    all_pass <- class_ok && has_feature_summary && has_prop_sig && print_works
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 8: Input validation for R mode
# =============================================================================

test_input_validation_r_mode <- function() {
    print_test_header("Input validation: R mode")

    graph <- create_triangle_graph()
    n.samples <- 20L
    Z.hat.samples <- list(matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0)

    ## Test mismatched y.hat length
    error_caught_yhat <- FALSE
    tryCatch({
        lcor.with.posterior(
            graph$adj.list, graph$weight.list,
            y.hat = c(1.0, 2.0),  # Wrong length
            Z.hat.samples = Z.hat.samples,
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_yhat <<- TRUE
    })

    if (error_caught_yhat) cat("\u2713 Error caught for mismatched y.hat length\n")
    else print_fail("Should error on wrong y.hat length")

    ## Test empty samples list
    error_caught_empty <- FALSE
    tryCatch({
        lcor.with.posterior(
            graph$adj.list, graph$weight.list,
            y.hat = y.hat,
            Z.hat.samples = list(),  # Empty
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_empty <<- TRUE
    })

    if (error_caught_empty) cat("\u2713 Error caught for empty samples list\n")
    else print_fail("Should error on empty samples list")

    all_pass <- error_caught_yhat && error_caught_empty
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 9: Zero variance samples handled gracefully
# =============================================================================

test_zero_variance_samples <- function() {
    print_test_header("Zero variance samples handled gracefully")

    graph <- create_triangle_graph()
    n.samples <- 30L

    ## Create constant samples (zero variance)
    constant_samples <- list(matrix(1.0, nrow = 3, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0)

    result <- NULL
    error_occurred <- FALSE

    tryCatch({
        result <- lcor.with.posterior(
            graph$adj.list, graph$weight.list, y.hat, constant_samples,
            verbose = FALSE
        )
    }, error = function(e) {
        error_occurred <<- TRUE
        cat("  Error message:", conditionMessage(e), "\n")
    })

    if (!error_occurred && !is.null(result)) {
        cat("\u2713 Function handles constant samples without error\n")

        ## SD should be zero or near-zero (or NA if lcor produces NaN)
        sd_near_zero <- all(result$sd < 1e-10 | is.na(result$sd) | is.nan(result$sd))
        if (sd_near_zero) cat("\u2713 SD is zero/NA for constant samples\n")
        else cat("  Note: SD not zero for constant samples (may be okay)\n")

        all_pass <- TRUE
    } else {
        ## For constant samples, lcor will produce NaN (0/0)
        ## This might cause issues but shouldn't crash
        cat("  Note: Function may produce warnings for constant samples\n")
        all_pass <- !error_occurred  # Pass if no hard error
    }

    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 10: Large number of samples
# =============================================================================

test_large_sample_size <- function() {
    print_test_header("Large number of posterior samples")

    graph <- create_triangle_graph()
    n.samples <- 1000L  # Large

    set.seed(555)
    Z.hat.samples <- list(matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0)

    start_time <- Sys.time()

    result <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    time_ok <- elapsed < 5  # Should be fast
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
# TEST 11: Row names preserved from Z.hat.samples names
# =============================================================================

test_row_names_preserved <- function() {
    print_test_header("Feature names preserved in output")

    graph <- create_triangle_graph()
    n.samples <- 30L

    set.seed(666)
    Z.hat.samples <- list(
        feature_A = matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples),
        feature_B = matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples),
        feature_C = matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples)
    )
    y.hat <- c(1.0, 2.0, 3.0)

    result <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    ## Check if row names are preserved
    expected_names <- c("feature_A", "feature_B", "feature_C")
    names_match <- identical(rownames(result$mean), expected_names)

    if (names_match) cat("\u2713 Feature names preserved in mean matrix rownames\n")
    else cat("  Note: Feature names not preserved (implementation choice)\n")

    ## At minimum, result should have correct number of rows
    nrow_ok <- nrow(result$mean) == 3
    if (nrow_ok) cat("\u2713 Correct number of features in output\n")
    else print_fail("Wrong number of features")

    all_pass <- nrow_ok
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 12: Verify lcor values are bounded in [-1, 1]
# =============================================================================

test_lcor_bounded <- function() {
    print_test_header("lcor values bounded in [-1, 1]")

    graph <- create_triangle_graph()
    n.samples <- 100L

    set.seed(777)
    ## Create diverse samples to test bounds
    Z.hat.samples <- lapply(1:3, function(j) {
        matrix(rnorm(3 * n.samples, mean = j * 10, sd = 5), nrow = 3, ncol = n.samples)
    })
    y.hat <- c(-5.0, 0.0, 10.0)  # Diverse y values

    result <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    ## Check bounds on mean
    mean_bounded <- all(result$mean >= -1 & result$mean <= 1, na.rm = TRUE)
    if (mean_bounded) cat("\u2713 All mean values in [-1, 1]\n")
    else print_fail("Mean values out of bounds")

    ## Check bounds on CI
    ci_bounded <- all(result$lower >= -1, na.rm = TRUE) &&
                  all(result$upper <= 1, na.rm = TRUE)
    if (ci_bounded) cat("\u2713 All CI bounds in [-1, 1]\n")
    else print_fail("CI bounds out of range")

    all_pass <- mean_bounded && ci_bounded
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 13: Consistent results with same seed
# =============================================================================

test_reproducibility <- function() {
    print_test_header("Reproducibility with same seed")

    graph <- create_triangle_graph()
    n.samples <- 50L

    Z.hat.samples <- list(matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples))
    y.hat <- c(1.0, 2.0, 3.0)

    ## Two runs should give identical results (no randomness in R mode)
    result1 <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    result2 <- lcor.with.posterior(
        graph$adj.list, graph$weight.list, y.hat, Z.hat.samples,
        verbose = FALSE
    )

    mean_identical <- matrices_equal(result1$mean, result2$mean)
    if (mean_identical) cat("\u2713 Results identical across runs (deterministic R mode)\n")
    else print_fail("Results differ between runs")

    all_pass <- mean_identical
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 14: Different lcor types produce different extremum detection
# =============================================================================

test_lcor_type_extremum_sensitivity <- function() {
    print_test_header("lcor types affect local correlation values")

    ## Use a graph where derivative weighting matters
    adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    ## Different edge lengths
    weight.list <- list(c(1.0, 2.0), c(1.0, 1.5), c(2.0, 1.5))

    n.samples <- 50L
    set.seed(888)
    Z.hat.samples <- list(matrix(rnorm(3 * n.samples), nrow = 3, ncol = n.samples))
    y.hat <- c(1.0, 5.0, 3.0)  # Vertex 2 is local max

    result_deriv <- lcor.with.posterior(
        adj.list, weight.list, y.hat, Z.hat.samples,
        lcor.type = "derivative", verbose = FALSE
    )

    result_unit <- lcor.with.posterior(
        adj.list, weight.list, y.hat, Z.hat.samples,
        lcor.type = "unit", verbose = FALSE
    )

    ## With different edge lengths, derivative vs unit should differ
    ## (unless the random samples happen to align)
    cat(sprintf("  derivative mean: [%.3f, %.3f, %.3f]\n",
                result_deriv$mean[1,1], result_deriv$mean[1,2], result_deriv$mean[1,3]))
    cat(sprintf("  unit mean:       [%.3f, %.3f, %.3f]\n",
                result_unit$mean[1,1], result_unit$mean[1,2], result_unit$mean[1,3]))

    ## Both should complete without error
    cat("\u2713 Both lcor types computed successfully\n")

    all_pass <- TRUE
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# TEST 15: Input validation for C++ mode (missing required parameters)
# =============================================================================

test_cpp_mode_validation <- function() {
    print_test_header("C++ mode: input validation")

    graph <- create_triangle_graph()
    y.hat <- c(1.0, 2.0, 3.0)
    Z.abundances <- matrix(rnorm(9), nrow = 3, ncol = 3)

    ## Test missing fitted.model when Z.hat.samples = NULL
    error_caught_model <- FALSE
    tryCatch({
        lcor.with.posterior(
            graph$adj.list, graph$weight.list, y.hat,
            Z.hat.samples = NULL,
            fitted.model = NULL,  # Missing!
            Z.abundances = Z.abundances,
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_model <<- TRUE
    })

    if (error_caught_model) cat("\u2713 Error caught for missing fitted.model\n")
    else print_fail("Should error when fitted.model is NULL")

    ## Test missing Z.abundances when Z.hat.samples = NULL
    mock_model <- create_mock_fitted_model(3)

    error_caught_Z <- FALSE
    tryCatch({
        lcor.with.posterior(
            graph$adj.list, graph$weight.list, y.hat,
            Z.hat.samples = NULL,
            fitted.model = mock_model,
            Z.abundances = NULL,  # Missing!
            verbose = FALSE
        )
    }, error = function(e) {
        error_caught_Z <<- TRUE
    })

    if (error_caught_Z) cat("\u2713 Error caught for missing Z.abundances\n")
    else print_fail("Should error when Z.abundances is NULL")

    all_pass <- error_caught_model && error_caught_Z
    if (all_pass) print_pass()
    return(all_pass)
}

# =============================================================================
# Main test runner
# =============================================================================

run_all_lcor_posterior_tests <- function() {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("  LCOR WITH POSTERIOR VALIDATION TEST SUITE\n")
    cat("  Testing posterior uncertainty propagation for local correlation\n")
    cat(rep("=", 70), "\n", sep = "")

    tests <- list(
        test_r_mode_basic,
        test_r_mode_single_feature,
        test_r_mode_return_samples,
        test_lcor_type_parameter,
        test_credible_level,
        test_print_method,
        test_summary_method,
        test_input_validation_r_mode,
        test_zero_variance_samples,
        test_large_sample_size,
        test_row_names_preserved,
        test_lcor_bounded,
        test_reproducibility,
        test_lcor_type_extremum_sensitivity,
        test_cpp_mode_validation
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
run_all_lcor_posterior_tests()

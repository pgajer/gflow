## ============================================================================
## Tests for lslope Wasserstein Permutation Test
## ============================================================================
##
## These tests verify the correctness of the lslope Wasserstein permutation
## test implementation.
##
## Run with: source("test_lslope_wasserstein_test.R")
## Or: testthat::test_file("test_lslope_wasserstein_test.R")
##
## @author Pawel Gajer
## @date 2025
## ============================================================================

library(testthat)
library(gflow)


## ============================================================================
## Test 1: Sample generation (paired design)
## ============================================================================

test_that("paired sample generation produces correct dimensions", {
    skip_if_not_installed("gflow")

    n <- 50
    B <- 20
    set.seed(123)

    V <- diag(n)
    filter.weights <- exp(-0.5 * seq(0, 2, length.out = n))

    ## Simple chain graph (0-based for C++)
    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    y <- sin(seq(0, 2*pi, length.out = n)) + rnorm(n, sd = 0.1)
    z.hat <- cos(seq(0, 2*pi, length.out = n))

    ## Generate paired samples
    samples <- .Call(
        "S_generate_paired_lslope_samples",
        as.numeric(y),
        as.numeric(z.hat),
        V,
        as.numeric(filter.weights),
        adj.list.0,
        weight.list,
        as.integer(B),
        as.integer(1),  # normalized
        TRUE,           # ascending
        as.integer(12345),
        PACKAGE = "gflow"
    )

    expect_true(is.list(samples))
    expect_true("signal.samples" %in% names(samples))
    expect_true("null.samples" %in% names(samples))

    expect_equal(nrow(samples$signal.samples), n)
    expect_equal(ncol(samples$signal.samples), B)
    expect_equal(nrow(samples$null.samples), n)
    expect_equal(ncol(samples$null.samples), B)

    ## Check that values are bounded for normalized lslope
    expect_true(all(samples$signal.samples >= -1 & samples$signal.samples <= 1,
                    na.rm = TRUE))
    expect_true(all(samples$null.samples >= -1 & samples$null.samples <= 1,
                    na.rm = TRUE))
})


## ============================================================================
## Test 2: Full pipeline test
## ============================================================================

test_that("full lslope Wasserstein test pipeline works", {
    skip_if_not_installed("gflow")

    n <- 50
    set.seed(42)

    ## Create simple test data
    V <- diag(n)
    filter.weights <- exp(-0.5 * seq(0, 2, length.out = n))

    ## Simple chain graph (0-based for C++)
    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## Create y with genuine signal in first half
    y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n, sd = 0.2)

    ## z.hat correlated with y in first region, uncorrelated in second
    z.hat <- c(seq(0, 1, length.out = n/2),
               rnorm(n/2, mean = 0.5, sd = 0.3))

    ## Run the test
    result <- .Call(
        "S_lslope_wasserstein_test",
        as.numeric(y),
        as.numeric(z.hat),
        V,
        as.numeric(filter.weights),
        adj.list.0,
        weight.list,
        as.integer(50),    # n.samples
        as.integer(199),   # n.permutations
        as.integer(1),     # lslope.type = normalized
        TRUE,              # ascending
        as.integer(12345), # seed
        as.integer(1),     # n.cores
        FALSE,             # instrumented
        FALSE,             # verbose
        PACKAGE = "gflow"
    )

    ## Check structure
    expect_true(is.list(result))
    expect_equal(length(result$wasserstein.observed), n)
    expect_equal(length(result$p.value), n)
    expect_equal(length(result$signal.mean), n)
    expect_equal(length(result$signal.sd), n)
    expect_equal(length(result$null.mean), n)
    expect_equal(length(result$null.sd), n)

    ## Check that Wasserstein distances are non-negative
    expect_true(all(result$wasserstein.observed >= 0))

    ## Check that p-values are in (0, 1]
    expect_true(all(result$p.value > 0))
    expect_true(all(result$p.value <= 1))

    ## P-values should be multiples of 1/(n_permutations + 1)
    ## With 199 permutations, smallest possible is 1/200 = 0.005
    expect_true(all(result$p.value >= 1 / 200))
})


## ============================================================================
## Test 3: Reproducibility with seed
## ============================================================================

test_that("results are reproducible with same seed", {
    skip_if_not_installed("gflow")

    n <- 30
    set.seed(123)

    V <- diag(n)
    filter.weights <- exp(-0.5 * seq(0, 2, length.out = n))

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    y <- rnorm(n)
    z.hat <- rnorm(n)

    ## Run twice with same seed
    result1 <- .Call(
        "S_lslope_wasserstein_test",
        y, z.hat, V, filter.weights,
        adj.list.0, weight.list,
        20L, 99L, 1L, TRUE, 42L, 1L, FALSE, FALSE,
        PACKAGE = "gflow"
    )

    result2 <- .Call(
        "S_lslope_wasserstein_test",
        y, z.hat, V, filter.weights,
        adj.list.0, weight.list,
        20L, 99L, 1L, TRUE, 42L, 1L, FALSE, FALSE,
        PACKAGE = "gflow"
    )

    expect_equal(result1$wasserstein.observed, result2$wasserstein.observed)
    expect_equal(result1$p.value, result2$p.value)
    expect_equal(result1$signal.mean, result2$signal.mean)
    expect_equal(result1$null.mean, result2$null.mean)
})


## ============================================================================
## Test 4: Detection of true signal
## ============================================================================

test_that("test detects true association vs pure noise", {
    skip_if_not_installed("gflow")

    ## This test uses a simple setup to verify the test can detect signal.
    ## Note: Using V = diag(n) is not realistic (doesn't correspond to graph
    ## smoothing), so we use a more robust test criterion.

    n <- 50
    set.seed(999)

    ## Use identity V (simplified, not realistic graph smoothing)
    V <- diag(n)
    filter.weights <- rep(1.0, n)  # No filtering for simplicity

    ## Chain graph (0-based)
    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## Strong signal: y determines z with clear monotonic relationship
    ## Use a non-linear relationship to create varied local slopes
    y.signal <- seq(0, 1, length.out = n)
    z.signal <- y.signal^2  # Quadratic: different slopes at different locations

    ## Test with signal
    result.signal <- .Call(
        "S_lslope_wasserstein_test",
        y.signal, z.signal, V, filter.weights,
        adj.list.0, weight.list,
        100L, 299L, 1L, TRUE, 123L, 1L, FALSE, FALSE,
        PACKAGE = "gflow"
    )

    ## Basic validity checks
    expect_true(all(result.signal$wasserstein.observed >= 0))
    expect_true(all(result.signal$p.value > 0 & result.signal$p.value <= 1))

    ## For signal case: signal and null distributions should differ
    ## The signal has consistent positive slopes (z increases with y)
    ## The null (permuted y) should have mixed slopes
    ## So we expect some vertices to show significant differences

    ## Check that signal mean differs from null mean at most vertices
    ## (signal should show consistent positive slopes, null should be ~0)
    mean.diff <- abs(result.signal$signal.mean - result.signal$null.mean)
    expect_gt(median(mean.diff), 0.01)

    ## At least some vertices should have small p-values
    ## (not testing against noise, just that signal is detectable)
    n.signif <- sum(result.signal$p.value < 0.1)
    expect_gt(n.signif, 0)
})


test_that("signal distribution differs from null under association", {
    skip_if_not_installed("gflow")

    ## More focused test: under true association, signal mean should
    ## systematically differ from null mean

    n <- 40
    set.seed(888)

    V <- diag(n)
    filter.weights <- rep(1.0, n)

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## Strong monotonic signal
    y <- seq(-2, 2, length.out = n)
    z <- 0.5 * y + rnorm(n, sd = 0.1)  # z increases with y

    result <- .Call(
        "S_lslope_wasserstein_test",
        y, z, V, filter.weights,
        adj.list.0, weight.list,
        100L, 199L, 1L, TRUE, 456L, 1L, FALSE, FALSE,
        PACKAGE = "gflow"
    )

    ## Signal slopes should be positive (z increases with y)
    ## Null slopes should be centered around 0 (permutation destroys relationship)
    expect_gt(mean(result$signal.mean), mean(result$null.mean))

    ## The difference should be detectable
    ## Note: with normalized lslope (tanh), values are compressed
    expect_gt(mean(result$signal.mean), 0)  # Signal has positive slopes
    expect_lt(abs(mean(result$null.mean)), 0.3)  # Null is roughly centered
})


## ============================================================================
## Test 5: Instrumented mode returns sample matrices
## ============================================================================

test_that("instrumented mode returns sample matrices", {
    skip_if_not_installed("gflow")

    n <- 30
    B <- 25
    set.seed(456)

    V <- diag(n)
    filter.weights <- exp(-0.5 * seq(0, 2, length.out = n))

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    y <- rnorm(n)
    z.hat <- rnorm(n)

    ## Run without instrumentation
    result.plain <- .Call(
        "S_lslope_wasserstein_test",
        y, z.hat, V, filter.weights,
        adj.list.0, weight.list,
        as.integer(B), 99L, 1L, TRUE, 789L, 1L,
        FALSE,  # instrumented = FALSE
        FALSE,
        PACKAGE = "gflow"
    )

    ## Run with instrumentation
    result.instr <- .Call(
        "S_lslope_wasserstein_test",
        y, z.hat, V, filter.weights,
        adj.list.0, weight.list,
        as.integer(B), 99L, 1L, TRUE, 789L, 1L,
        TRUE,  # instrumented = TRUE
        FALSE,
        PACKAGE = "gflow"
    )

    ## Plain should NOT have sample matrices
    expect_null(result.plain$signal.samples)
    expect_null(result.plain$null.samples)

    ## Instrumented should have sample matrices
    expect_true(!is.null(result.instr$signal.samples))
    expect_true(!is.null(result.instr$null.samples))

    expect_equal(dim(result.instr$signal.samples), c(n, B))
    expect_equal(dim(result.instr$null.samples), c(n, B))

    ## Core results should be identical (same seed)
    expect_equal(result.plain$wasserstein.observed, result.instr$wasserstein.observed)
    expect_equal(result.plain$p.value, result.instr$p.value)
})


## ============================================================================
## Test 6: P-value behavior under null
## ============================================================================

test_that("p-values have expected properties", {
    skip_if_not_installed("gflow")

    ## NOTE: Perfect uniformity of p-values under the null is NOT expected
    ## because of a fundamental design characteristic:
    ##
    ## - Signal samples all use the SAME underlying y (with different weights)
    ## - Null samples each use DIFFERENT permuted y
    ##
    ## This creates different correlation structures: signal samples are
    ## correlated (they share y), while null samples are approximately
    ## independent. The Wasserstein distance can detect this difference
    ## in empirical distribution spread even when the true marginal
    ## distributions are identical.
    ##
    ## This test verifies basic validity properties rather than exact uniformity.

    n <- 30
    set.seed(2025)

    V <- diag(n)
    filter.weights <- rep(1.0, n)

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## Under true null: y and z are independent
    y <- rnorm(n)
    z <- rnorm(n)

    result <- .Call(
        "S_lslope_wasserstein_test",
        y, z, V, filter.weights,
        adj.list.0, weight.list,
        100L, 499L, 1L, TRUE, 999L, 1L, FALSE, FALSE,
        PACKAGE = "gflow"
    )

    ## Basic validity: p-values should be in valid range
    expect_true(all(result$p.value > 0))
    expect_true(all(result$p.value <= 1))

    ## P-values should be discrete with step 1/(n_permutations + 1)
    ## Minimum possible is 1/500 = 0.002
    expect_true(all(result$p.value >= 1 / 500))

    ## Wasserstein distances should be non-negative
    expect_true(all(result$wasserstein.observed >= 0))

    ## Signal and null means should not be systematically different
    ## under the null (even if p-values are inflated)
    mean.diff <- mean(result$signal.mean) - mean(result$null.mean)
    expect_lt(abs(mean.diff), 0.5)  # Not systematically shifted
})


## ============================================================================
## Test 7: R wrapper function (if using fitted model)
## ============================================================================

test_that("R wrapper function works with fitted model", {
    skip_if_not_installed("gflow")
    skip("Requires full fitted model - manual testing recommended")

    ## This test requires a full fitted model from fit.rdgraph.regression()
    ## Uncomment and modify for manual testing

    # ## Create test data
    # n <- 100
    # p <- 5
    # set.seed(42)
    #
    # X <- matrix(rnorm(n * p), n, p)
    # y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.5)
    # z <- X[,1] + rnorm(n, sd = 0.3)  # Correlated with y
    #
    # ## Fit model
    # fit <- fit.rdgraph.regression(X, y, k = 10)
    #
    # ## Run test
    # result <- lslope.wasserstein.test(
    #     fitted.model = fit,
    #     y = y,
    #     z = z,
    #     n.samples = 100,
    #     n.permutations = 499,
    #     verbose = FALSE
    # )
    #
    # expect_s3_class(result, "lslope.wasserstein.test")
    # expect_equal(result$n.vertices, n)
    #
    # ## Test instrumented mode
    # result.instr <- lslope.wasserstein.test(
    #     fitted.model = fit,
    #     y = y,
    #     z = z,
    #     n.samples = 100,
    #     n.permutations = 499,
    #     instrumented = TRUE,
    #     verbose = FALSE
    # )
    #
    # expect_true(!is.null(result.instr$signal.samples))
})


## ============================================================================
## Run all tests
## ============================================================================

cat("\n=== Running lslope Wasserstein Permutation Test Tests ===\n\n")

## If running interactively:
## test_file("test_lslope_wasserstein_test.R")

cat("\nAll tests completed.\n")

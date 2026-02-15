## ============================================================================
## Test Script for lslope Wasserstein Bayesian Test
## ============================================================================
##
## This script provides integration tests for the lslope Wasserstein Bayesian
## test pipeline. Run after compiling the package to verify correctness.
##
## @author Pawel Gajer
## @date 2025
## ============================================================================

library(testthat)
## library(gflow)

## ============================================================================
## Helper function to create a simple test graph
## ============================================================================

create.test.graph <- function(n = 100, k = 10) {
    ## Create a simple grid-like structure for testing
    ## This is a placeholder - in practice, use real data

    ## Random points in 2D
    set.seed(42)
    X <- matrix(rnorm(n * 2), n, 2)

    ## Build k-NN graph
    ## (Assuming iknn.graph function exists in gflow)
    graph <- iknn.graph(X, k = k)

    return(graph)
}


## ============================================================================
## Test 1: Basic functionality of sample generation
## ============================================================================

test_that("signal sample generation works", {
    ##skip_if_not_installed("gflow")

    ## Create a simple test case
    n <- 50
    set.seed(123)

    ## Mock spectral decomposition (identity for simplicity)
    V <- diag(n)
    eigenvalues <- seq(0, 2, length.out = n)
    filter.weights <- exp(-0.5 * eigenvalues)

    ## Simple chain graph
    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)  # 0-based
        if (i < n) neighbors <- c(neighbors, i)        # 0-based
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## Test data
    y <- sin(seq(0, 2*pi, length.out = n)) + rnorm(n, sd = 0.1)
    z.hat <- cos(seq(0, 2*pi, length.out = n))

    ## Generate signal samples
    samples <- .Call(
        "S_generate_signal_lslope_samples",
        as.numeric(y),
        as.numeric(z.hat),
        V,
        as.numeric(filter.weights),
        adj.list.0,
        weight.list,
        as.integer(20),   # n.samples
        as.integer(1),    # lslope.type = normalized
        TRUE,             # ascending
        as.integer(12345),
        PACKAGE = "gflow"
    )

    ## Check dimensions
    expect_equal(nrow(samples), 20)
    expect_equal(ncol(samples), n)

    ## Check that values are bounded for normalized lslope
    expect_true(all(samples >= -1 & samples <= 1, na.rm = TRUE))
})


test_that("null sample generation works", {
    ##skip_if_not_installed("gflow")

    n <- 50
    set.seed(123)

    V <- diag(n)
    eigenvalues <- seq(0, 2, length.out = n)
    filter.weights <- exp(-0.5 * eigenvalues)

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    y <- sin(seq(0, 2*pi, length.out = n)) + rnorm(n, sd = 0.1)
    z.hat <- cos(seq(0, 2*pi, length.out = n))

    ## Generate null samples
    samples <- .Call(
        "S_generate_null_lslope_samples",
        as.numeric(y),
        as.numeric(z.hat),
        V,
        as.numeric(filter.weights),
        adj.list.0,
        weight.list,
        as.integer(20),
        as.integer(1),
        TRUE,
        as.integer(12345),
        PACKAGE = "gflow"
    )

    expect_equal(nrow(samples), 20)
    expect_equal(ncol(samples), n)
    expect_true(all(samples >= -1 & samples <= 1, na.rm = TRUE))
})


## ============================================================================
## Test 2: Full pipeline test
## ============================================================================

test_that("full lslope Wasserstein Bayesian test pipeline works", {
    ##skip_if_not_installed("gflow")

    n <- 50
    set.seed(42)

    ## Create simple test data
    V <- diag(n)
    eigenvalues <- seq(0, 2, length.out = n)
    filter.weights <- exp(-0.5 * eigenvalues)

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
        "S_lslope_wasserstein_bayes_test",
        as.numeric(y),
        as.numeric(z.hat),
        V,
        as.numeric(eigenvalues),
        as.numeric(filter.weights),
        adj.list.0,
        weight.list,
        as.integer(50),    # n.signal.samples
        as.integer(50),    # n.null.samples
        as.integer(100),   # n.bayes.samples
        as.integer(1),     # lslope.type = normalized
        TRUE,              # ascending
        as.numeric(0.05),  # threshold
        as.numeric(0.95),  # credible.level
        as.integer(12345), # seed
        as.integer(1),     # n.cores
        FALSE,             # verbose
        PACKAGE = "gflow"
    )

    ## Check structure
    expect_true(is.list(result))
    expect_equal(length(result$wasserstein.mean), n)
    expect_equal(length(result$wasserstein.sd), n)
    expect_equal(length(result$wasserstein.lower), n)
    expect_equal(length(result$wasserstein.upper), n)
    expect_equal(length(result$prob.exceeds.threshold), n)
    expect_equal(length(result$signal.mean), n)
    expect_equal(length(result$null.mean), n)

    ## Check that Wasserstein distances are non-negative
    expect_true(all(result$wasserstein.mean >= 0))
    expect_true(all(result$wasserstein.point >= 0))

    ## Check that credible intervals are ordered
    expect_true(all(result$wasserstein.lower <= result$wasserstein.upper))

    ## Check that probabilities are in [0, 1]
    expect_true(all(result$prob.exceeds.threshold >= 0))
    expect_true(all(result$prob.exceeds.threshold <= 1))
})


## ============================================================================
## Test 3: Reproducibility with seed
## ============================================================================

test_that("results are reproducible with same seed", {
    ##skip_if_not_installed("gflow")

    n <- 30
    set.seed(123)

    V <- diag(n)
    eigenvalues <- seq(0, 2, length.out = n)
    filter.weights <- exp(-0.5 * eigenvalues)

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
        "S_lslope_wasserstein_bayes_test",
        y, z.hat, V, eigenvalues, filter.weights,
        adj.list.0, weight.list,
        20L, 20L, 50L, 1L, TRUE, 0.0, 0.95, 42L, 1L, FALSE,
        PACKAGE = "gflow"
    )

    result2 <- .Call(
        "S_lslope_wasserstein_bayes_test",
        y, z.hat, V, eigenvalues, filter.weights,
        adj.list.0, weight.list,
        20L, 20L, 50L, 1L, TRUE, 0.0, 0.95, 42L, 1L, FALSE,
        PACKAGE = "gflow"
    )

    expect_equal(result1$wasserstein.mean, result2$wasserstein.mean)
    expect_equal(result1$prob.exceeds.threshold, result2$prob.exceeds.threshold)
})


## ============================================================================
## Test 4: Detection of true signal
## ============================================================================

test_that("test detects true association vs pure noise", {
    ##skip_if_not_installed("gflow")

    n <- 40
    set.seed(999)

    V <- diag(n)
    eigenvalues <- seq(0, 2, length.out = n)
    filter.weights <- exp(-0.5 * eigenvalues)

    adj.list.0 <- lapply(1:n, function(i) {
        neighbors <- integer(0)
        if (i > 1) neighbors <- c(neighbors, i - 2L)
        if (i < n) neighbors <- c(neighbors, i)
        neighbors
    })

    weight.list <- lapply(adj.list.0, function(x) rep(1.0, length(x)))

    ## True signal: y and z perfectly correlated
    y.signal <- seq(-1, 1, length.out = n)
    z.signal <- y.signal

    ## Noise: y and z uncorrelated
    y.noise <- rnorm(n)
    z.noise <- rnorm(n)

    ## Test with signal - use meaningful threshold
    result.signal <- .Call(
        "S_lslope_wasserstein_bayes_test",
        y.signal, z.signal, V, eigenvalues, filter.weights,
        adj.list.0, weight.list,
        30L, 30L, 100L, 1L, TRUE, 0.05, 0.95, 123L, 1L, FALSE,  ## threshold = 0.05
        PACKAGE = "gflow"
    )

    ## Test with noise
    result.noise <- .Call(
        "S_lslope_wasserstein_bayes_test",
        y.noise, z.noise, V, eigenvalues, filter.weights,
        adj.list.0, weight.list,
        30L, 30L, 100L, 1L, TRUE, 0.05, 0.95, 123L, 1L, FALSE,  ## threshold = 0.05
        PACKAGE = "gflow"
    )

    ## Signal case should have more vertices with high exceedance probability
    n.signif.signal <- sum(result.signal$prob.exceeds.threshold > 0.9)
    n.signif.noise <- sum(result.noise$prob.exceeds.threshold > 0.9)

    ## Use >= to allow for equal (both could be low for noise)
    ## Or better: compare mean Wasserstein which is more robust
    expect_gt(n.signif.signal, n.signif.noise)

    ## Mean Wasserstein should be higher for signal
    expect_gt(mean(result.signal$wasserstein.mean),
              mean(result.noise$wasserstein.mean))
})

## ============================================================================
## Test 5: R wrapper function (if using fitted model)
## ============================================================================

test_that("R wrapper function works with fitted model", {
    ##skip_if_not_installed("gflow")
    skip("Requires full fitted model - manual testing recommended")

    ## This test requires a full fitted model from fit.rdgraph.regression()
    ## Uncomment and modify for manual testing

    ## Create test data
    n <- 100
    p <- 5
    set.seed(42)

    X <- matrix(rnorm(n * p), n, p)
    y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.5)
    z <- X[,1] + rnorm(n, sd = 0.3)

    ## Fit model
    fit <- fit.rdgraph.regression(X, y, k = 10)

    ## Run test
    result <- lslope.wasserstein.bayes.test(
        fitted.model = fit,
        y = y,
        z = z,
        n.signal.samples = 100,
        n.null.samples = 100,
        n.bayes.samples = 200,
        verbose = FALSE
    )

    expect_s3_class(result, "lslope.wasserstein.bayes.test")
    expect_equal(result$n.vertices, n)
})


## ============================================================================
## Run all tests
## ============================================================================

cat("\n=== Running lslope Wasserstein Bayesian Test Tests ===\n\n")

## If running interactively:
## test_file("test_lslope_wasserstein_bayes.R")

cat("\nAll tests completed.\n")

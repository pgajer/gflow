## ============================================================================
## Test Script for Weighted Wasserstein Distance Functions
## ============================================================================
##
## This script provides unit tests for the weighted Wasserstein distance
## implementation. Run after compiling the package to verify correctness.
##
## @author Pawel Gajer
## @date 2025
## ============================================================================

library(testthat)

## ============================================================================
## Test 1: Equal-weight Wasserstein distance (basic)
## ============================================================================

test_that("wasserstein.distance computes correct value for simple case", {
    ## Two identical distributions should have distance 0
    x <- c(1, 2, 3, 4, 5)
    expect_equal(wasserstein.distance(x, x), 0.0)

    ## Shifted distributions
    y <- x + 1
    ## W_1 = (1/5) * sum(|x_i - y_i|) = (1/5) * 5 = 1
    expect_equal(wasserstein.distance(x, y), 1.0)

    ## Different distributions
    x <- c(1, 2, 3)
    y <- c(2, 3, 4)
    ## After sorting: x = (1,2,3), y = (2,3,4)
    ## W_1 = (1/3) * (|1-2| + |2-3| + |3-4|) = (1/3) * 3 = 1
    expect_equal(wasserstein.distance(x, y), 1.0)
})


## ============================================================================
## Test 2: Weighted Wasserstein distance (manual calculation)
## ============================================================================

test_that("wasserstein.distance.weighted computes correct value", {
    ## Example from algorithm derivation:
    ## x = [1, 3], weights_x = [0.4, 0.6]
    ## y = [2, 4], weights_y = [0.5, 0.5]
    ##
    ## Intervals:
    ## [0, 0.4]: |1-2| * 0.4 = 0.4
    ## [0.4, 0.5]: |3-2| * 0.1 = 0.1
    ## [0.5, 1.0]: |3-4| * 0.5 = 0.5
    ## Total: 1.0

    x <- c(1, 3)
    w.x <- c(0.4, 0.6)
    y <- c(2, 4)
    w.y <- c(0.5, 0.5)

    expect_equal(wasserstein.distance.weighted(x, w.x, y, w.y), 1.0)
})

test_that("wasserstein.distance.weighted equals unweighted for uniform weights", {
    set.seed(123)
    n <- 50
    x <- rnorm(n)
    y <- rnorm(n, mean = 0.5)

    ## Uniform weights
    w <- rep(1/n, n)

    W.weighted <- wasserstein.distance.weighted(x, w, y, w)
    W.unweighted <- wasserstein.distance(x, y)

    expect_equal(W.weighted, W.unweighted, tolerance = 1e-10)
})


## ============================================================================
## Test 3: Bayesian bootstrap posterior
## ============================================================================

test_that("bayesian.bootstrap.wasserstein returns correct length", {
    set.seed(42)
    signal <- rnorm(100, mean = 0.5)
    null <- rnorm(100, mean = 0)

    posterior <- bayesian.bootstrap.wasserstein(signal, null, n.bootstrap = 500)

    expect_length(posterior, 500)
    expect_true(all(posterior >= 0))  # Wasserstein distance is non-negative
})

test_that("bayesian.bootstrap.wasserstein is reproducible with seed", {
    signal <- rnorm(50, mean = 0.3)
    null <- rnorm(50, mean = 0)

    p1 <- bayesian.bootstrap.wasserstein(signal, null, n.bootstrap = 100, seed = 123)
    p2 <- bayesian.bootstrap.wasserstein(signal, null, n.bootstrap = 100, seed = 123)

    expect_equal(p1, p2)
})

test_that("bayesian.bootstrap.wasserstein has sensible posterior for identical distributions", {
    set.seed(42)
    x <- rnorm(200)

    ## Same distribution - Wasserstein should be near 0
    posterior <- bayesian.bootstrap.wasserstein(x, x, n.bootstrap = 500)

    ## Posterior mean should be near 0
    expect_lt(mean(posterior), 0.1)
})


## ============================================================================
## Test 4: Vertex-wise Bayesian bootstrap
## ============================================================================

test_that("vertex.bayesian.bootstrap.wasserstein returns correct structure", {
    set.seed(42)
    n.samples <- 50
    n.vertices <- 10

    signal.matrix <- matrix(rnorm(n.samples * n.vertices), n.samples, n.vertices)
    null.matrix <- matrix(rnorm(n.samples * n.vertices), n.samples, n.vertices)

    result <- vertex.bayesian.bootstrap.wasserstein(
        signal.matrix, null.matrix,
        n.bootstrap = 100,
        threshold = 0.1,
        credible.level = 0.95
    )

    expect_s3_class(result, "vertex.wasserstein.bayes.test")
    expect_length(result$wasserstein.mean, n.vertices)
    expect_length(result$wasserstein.sd, n.vertices)
    expect_length(result$wasserstein.lower, n.vertices)
    expect_length(result$wasserstein.upper, n.vertices)
    expect_length(result$prob.exceeds.threshold, n.vertices)

    ## All distances should be non-negative
    expect_true(all(result$wasserstein.mean >= 0))

    ## Credible intervals should be ordered
    expect_true(all(result$wasserstein.lower <= result$wasserstein.upper))

    ## Probabilities should be in [0, 1]
    expect_true(all(result$prob.exceeds.threshold >= 0))
    expect_true(all(result$prob.exceeds.threshold <= 1))
})

test_that("vertex.bayesian.bootstrap.wasserstein detects shifted distributions", {
    set.seed(42)
    n.samples <- 100
    n.vertices <- 20

    ## Null distribution at all vertices
    null.matrix <- matrix(rnorm(n.samples * n.vertices), n.samples, n.vertices)

    ## Signal: shifted at first 5 vertices, same as null at rest
    signal.matrix <- null.matrix
    signal.matrix[, 1:5] <- signal.matrix[, 1:5] + 1.0  # Clear shift

    result <- vertex.bayesian.bootstrap.wasserstein(
        signal.matrix, null.matrix,
        n.bootstrap = 200,
        threshold = 0.3,
        credible.level = 0.95
    )

    ## Shifted vertices should have higher P(W > threshold)
    shifted.probs <- result$prob.exceeds.threshold[1:5]
    unshifted.probs <- result$prob.exceeds.threshold[6:20]

    expect_true(mean(shifted.probs) > mean(unshifted.probs))
    expect_true(all(shifted.probs > 0.8))  # Should confidently detect shift
})


## ============================================================================
## Test 5: Edge cases
## ============================================================================

test_that("functions handle edge cases gracefully", {
    ## Empty vectors
    expect_equal(wasserstein.distance(numeric(0), numeric(0)), 0.0)
    expect_equal(wasserstein.distance.weighted(numeric(0), numeric(0),
                                                numeric(0), numeric(0)), 0.0)

    ## Single element
    expect_equal(wasserstein.distance(1, 2), 1.0)
    expect_equal(wasserstein.distance.weighted(1, 1, 2, 1), 1.0)
})


## ============================================================================
## Test 6: Comparison with R implementation (validation)
## ============================================================================

test_that("weighted Wasserstein matches pure R implementation", {
    ## Pure R implementation for validation
    wasserstein.r <- function(x, w.x, y, w.y) {
        ## Sort both distributions
        ord.x <- order(x)
        ord.y <- order(y)

        x.sorted <- x[ord.x]
        w.x.sorted <- w.x[ord.x]
        y.sorted <- y[ord.y]
        w.y.sorted <- w.y[ord.y]

        ## Cumulative weights
        cum.x <- c(0, cumsum(w.x.sorted))
        cum.y <- c(0, cumsum(w.y.sorted))

        ## Two-pointer merge
        W <- 0
        i <- 1
        j <- 1
        n <- length(x)
        m <- length(y)

        while (i <= n && j <= m) {
            t.left <- max(cum.x[i], cum.y[j])
            t.right <- min(cum.x[i + 1], cum.y[j + 1])

            if (t.right > t.left) {
                W <- W + abs(x.sorted[i] - y.sorted[j]) * (t.right - t.left)
            }

            if (cum.x[i + 1] < cum.y[j + 1] - 1e-12) {
                i <- i + 1
            } else if (cum.y[j + 1] < cum.x[i + 1] - 1e-12) {
                j <- j + 1
            } else {
                i <- i + 1
                j <- j + 1
            }
        }

        return(W)
    }

    ## Test with random data
    set.seed(456)
    for (trial in 1:5) {
        n <- sample(10:50, 1)
        m <- sample(10:50, 1)

        x <- rnorm(n)
        y <- rnorm(m, mean = runif(1, -1, 1))

        ## Random Dirichlet weights (simulated)
        w.x <- rgamma(n, 1, 1)
        w.x <- w.x / sum(w.x)
        w.y <- rgamma(m, 1, 1)
        w.y <- w.y / sum(w.y)

        W.cpp <- wasserstein.distance.weighted(x, w.x, y, w.y)
        W.r <- wasserstein.r(x, w.x, y, w.y)

        expect_equal(W.cpp, W.r, tolerance = 1e-10,
                     label = sprintf("Trial %d (n=%d, m=%d)", trial, n, m))
    }
})


## ============================================================================
## Run all tests
## ============================================================================

cat("\n=== Running Weighted Wasserstein Distance Tests ===\n\n")

## If running interactively, uncomment:
## test_file("test_wasserstein_weighted.R")

cat("\nAll tests completed.\n")

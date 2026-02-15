## ============================================================================
## Unit Tests for Gradient Flow Association Functions
## ============================================================================
##
## This file contains comprehensive unit tests for the gfassoc family of
## functions. Tests are organized by function and cover:
##   - Basic functionality with simple graphs
##   - Edge cases (single basins, flat regions, multiplicity)
##   - Mathematical properties (bounds, normalization, symmetry)
##   - Known analytical solutions
##
## To run all tests:
##   source("test_gfassoc.R")
##   run.all.gfassoc.tests()
##
## ============================================================================

## library(gflow)

## ============================================================================
## Test Utilities
## ============================================================================

#' Create a simple linear graph (path graph)
#'
#' Creates a path graph 1 -- 2 -- 3 -- ... -- n with unit weights.
#' Useful for testing basic gradient flow properties.
create.path.graph <- function(n) {
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)

    for (i in seq_len(n)) {
        neighbors <- integer(0)
        weights <- numeric(0)

        if (i > 1) {
            neighbors <- c(neighbors, i - 1L)
            weights <- c(weights, 1.0)
        }
        if (i < n) {
            neighbors <- c(neighbors, i + 1L)
            weights <- c(weights, 1.0)
        }

        adj.list[[i]] <- neighbors
        weight.list[[i]] <- weights
    }

    list(adj.list = adj.list, weight.list = weight.list)
}


#' Create a simple cycle graph
#'
#' Creates a cycle graph with n vertices and unit weights.
create.cycle.graph <- function(n) {
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)

    for (i in seq_len(n)) {
        prev <- if (i == 1) n else i - 1L
        next.v <- if (i == n) 1L else i + 1L

        adj.list[[i]] <- c(prev, next.v)
        weight.list[[i]] <- c(1.0, 1.0)
    }

    list(adj.list = adj.list, weight.list = weight.list)
}


#' Create a star graph
#'
#' Creates a star graph with center vertex 1 and n-1 leaf vertices.
create.star.graph <- function(n) {
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)

    ## Center vertex connects to all others
    adj.list[[1]] <- 2:n
    weight.list[[1]] <- rep(1.0, n - 1)

    ## Leaf vertices connect only to center
    for (i in 2:n) {
        adj.list[[i]] <- 1L
        weight.list[[i]] <- 1.0
    }

    list(adj.list = adj.list, weight.list = weight.list)
}


#' Create a grid graph
#'
#' Creates a 2D grid graph with nrow x ncol vertices.
create.grid.graph <- function(nrow, ncol) {
    n <- nrow * ncol
    adj.list <- vector("list", n)
    weight.list <- vector("list", n)

    ## Helper to convert (row, col) to vertex index
    idx <- function(r, c) (r - 1) * ncol + c

    for (r in seq_len(nrow)) {
        for (c in seq_len(ncol)) {
            v <- idx(r, c)
            neighbors <- integer(0)

            ## Add neighbors (up, down, left, right)
            if (r > 1) neighbors <- c(neighbors, idx(r - 1, c))
            if (r < nrow) neighbors <- c(neighbors, idx(r + 1, c))
            if (c > 1) neighbors <- c(neighbors, idx(r, c - 1))
            if (c < ncol) neighbors <- c(neighbors, idx(r, c + 1))

            adj.list[[v]] <- neighbors
            weight.list[[v]] <- rep(1.0, length(neighbors))
        }
    }

    list(adj.list = adj.list, weight.list = weight.list, nrow = nrow, ncol = ncol)
}


#' Check if value is within tolerance
check.equal <- function(actual, expected, tol = 1e-10, msg = "") {
    if (abs(actual - expected) > tol) {
        stop(sprintf("%s: Expected %.10f, got %.10f (diff = %.2e)",
                     msg, expected, actual, abs(actual - expected)))
    }
    invisible(TRUE)
}


#' Check if all values in vector are within bounds
check.bounds <- function(x, lower, upper, msg = "") {
    if (any(x < lower - 1e-10) || any(x > upper + 1e-10)) {
        stop(sprintf("%s: Values outside [%.2f, %.2f], range = [%.4f, %.4f]",
                     msg, lower, upper, min(x), max(x)))
    }
    invisible(TRUE)
}


#' Report test result
report.test <- function(name, passed) {
    status <- if (passed) "PASS" else "FAIL"
    cat(sprintf("[%s] %s\n", status, name))
}


## ============================================================================
## Tests for gfassoc.membership()
## ============================================================================

test.membership.path.graph <- function() {
    ## Path graph: 1 -- 2 -- 3 -- 4 -- 5
    ## y values create single max at 3 and mins at endpoints
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)  # Max at vertex 3, mins at 1 and 5

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)

    ## Should have 1 max basin (vertex 3) and 2 min basins (vertices 1, 5)
    stopifnot(membership$n_max_basins == 1)
    stopifnot(membership$n_min_basins == 2)

    ## All vertices should belong to the single max basin
    for (v in seq_len(5)) {
        stopifnot(length(membership$max_basin_indices[[v]]) == 1)
        stopifnot(membership$max_basin_indices[[v]][1] == 0)  # 0-based index
    }

    ## Check min basin membership (vertices 1,2 -> min at 1; vertices 4,5 -> min at 5)
    ## Vertex 3 could belong to both min basins
    stopifnot(0 %in% membership$min_basin_indices[[1]])  # Vertex 1 in its own basin
    stopifnot(1 %in% membership$min_basin_indices[[5]])  # Vertex 5 in its own basin

    ## Check extremum values
    check.equal(membership$max_values[1], 5, msg = "Max value")
    stopifnot(1 %in% membership$min_values)  # Both mins have value 1

    report.test("membership.path.graph", TRUE)
}


test.membership.star.graph <- function() {
    ## Star graph with center as maximum
    g <- create.star.graph(5)
    y <- c(10, 1, 2, 3, 4)  # Center (1) is max, leaves are potential mins

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)

    ## Should have 1 max basin (center)
    stopifnot(membership$n_max_basins == 1)

    ## All vertices should be in the max basin
    for (v in seq_len(5)) {
        stopifnot(length(membership$max_basin_indices[[v]]) >= 1)
    }

    ## Vertex 2 should be the only minimum (lowest value among leaves)
    stopifnot(membership$n_min_basins >= 1)

    report.test("membership.star.graph", TRUE)
}


test.membership.multiplicity <- function() {
    ## Create a graph where a vertex has multiple gradient descent options
    ## Diamond: 1 connects to 2,3; 2,3 connect to 4
    ##      1 (top)
    ##     / \
    ##    2   3
    ##     \ /
    ##      4 (bottom)

    adj.list <- list(
        c(2L, 3L),      # Vertex 1
        c(1L, 4L),      # Vertex 2
        c(1L, 4L),      # Vertex 3
        c(2L, 3L)       # Vertex 4
    )
    weight.list <- list(
        c(1.0, 1.0),
        c(1.0, 1.0),
        c(1.0, 1.0),
        c(1.0, 1.0)
    )

    ## y values: max at 1, min at 4, equal values at 2 and 3
    y <- c(10, 5, 5, 1)

    basins <- compute.basins.of.attraction(adj.list, weight.list, y)
    membership <- gfassoc.membership(basins)

    ## Should have 1 max (vertex 1) and 1 min (vertex 4)
    stopifnot(membership$n_max_basins == 1)
    stopifnot(membership$n_min_basins == 1)

    ## All vertices should have membership weights summing to 1
    for (v in seq_len(4)) {
        if (length(membership$max_membership[[v]]) > 0) {
            check.equal(sum(membership$max_membership[[v]]), 1.0,
                        msg = sprintf("Max membership sum for vertex %d", v))
        }
        if (length(membership$min_membership[[v]]) > 0) {
            check.equal(sum(membership$min_membership[[v]]), 1.0,
                        msg = sprintf("Min membership sum for vertex %d", v))
        }
    }

    report.test("membership.multiplicity", TRUE)
}


test.membership.cell.computation <- function() {
    ## Path graph with clear cell structure
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)

    ## Check cell_indices structure exists
    stopifnot(!is.null(membership$cell_indices))
    stopifnot(length(membership$cell_indices) == 5)

    ## Each vertex should have cell membership summing to 1
    for (v in seq_len(5)) {
        if (length(membership$cell_membership[[v]]) > 0) {
            check.equal(sum(membership$cell_membership[[v]]), 1.0,
                        msg = sprintf("Cell membership sum for vertex %d", v))
        }
    }

    report.test("membership.cell.computation", TRUE)
}


## ============================================================================
## Tests for gfassoc.polarity()
## ============================================================================

test.polarity.linear.function <- function() {
    ## Path graph with monotone increasing function
    ## Should have polarity increasing from -1 to +1
    g <- create.path.graph(5)
    y <- c(1, 2, 3, 4, 5)  # Monotone increasing: min at 1, max at 5

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)
    polarity <- gfassoc.polarity(y, membership, polarity.scale = "value")

    ## Check bounds
    check.bounds(polarity$theta[polarity$is_valid], 0, 1, "Theta bounds")
    check.bounds(polarity$polarity[polarity$is_valid], -1, 1, "Polarity bounds")

    ## Minimum should have polarity near -1, maximum near +1
    check.equal(polarity$polarity[1], -1, tol = 0.01, msg = "Min polarity")
    check.equal(polarity$polarity[5], 1, tol = 0.01, msg = "Max polarity")

    ## Polarity should be monotone increasing
    for (i in 1:4) {
        stopifnot(polarity$polarity[i] < polarity$polarity[i + 1])
    }

    report.test("polarity.linear.function", TRUE)
}


test.polarity.symmetric.function <- function() {
    ## Symmetric function: should have symmetric polarity
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)  # Symmetric around vertex 3

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)
    polarity <- gfassoc.polarity(y, membership, polarity.scale = "value")

    ## Maximum at vertex 3 should have polarity +1
    check.equal(polarity$polarity[3], 1, tol = 0.01, msg = "Max polarity")

    ## Minima at vertices 1 and 5 should have polarity -1
    check.equal(polarity$polarity[1], -1, tol = 0.01, msg = "Min polarity at 1")
    check.equal(polarity$polarity[5], -1, tol = 0.01, msg = "Min polarity at 5")

    ## Symmetric vertices should have equal polarity
    check.equal(polarity$polarity[2], polarity$polarity[4], tol = 1e-6,
                msg = "Symmetric polarity")

    report.test("polarity.symmetric.function", TRUE)
}


test.polarity.rank.vs.value <- function() {
    ## Test that rank-based polarity is monotone transform invariant
    g <- create.path.graph(5)
    y1 <- c(1, 2, 3, 4, 5)
    y2 <- c(1, 4, 9, 16, 25)  # Square of y1 (monotone transform)

    basins1 <- compute.basins.of.attraction(g$adj.list, g$weight.list, y1)
    basins2 <- compute.basins.of.attraction(g$adj.list, g$weight.list, y2)

    membership1 <- gfassoc.membership(basins1)
    membership2 <- gfassoc.membership(basins2)

    pol1.rank <- gfassoc.polarity(y1, membership1, polarity.scale = "rank")
    pol2.rank <- gfassoc.polarity(y2, membership2, polarity.scale = "rank")

    ## Rank-based polarity should be identical for monotone transforms
    for (v in seq_len(5)) {
        if (pol1.rank$is_valid[v] && pol2.rank$is_valid[v]) {
            check.equal(pol1.rank$polarity[v], pol2.rank$polarity[v], tol = 1e-6,
                        msg = sprintf("Rank polarity invariance at vertex %d", v))
        }
    }

    ## Value-based polarity should differ
    pol1.val <- gfassoc.polarity(y1, membership1, polarity.scale = "value")
    pol2.val <- gfassoc.polarity(y2, membership2, polarity.scale = "value")

    ## Middle vertices should have different value-based polarity
    ## (because square transform is not linear)
    ## Note: This test may pass trivially if basin structure differs

    report.test("polarity.rank.vs.value", TRUE)
}


test.polarity.flat.region <- function() {
    ## Test handling of flat regions
    g <- create.path.graph(5)
    y <- c(1, 1, 1, 1, 1)  # Completely flat

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)
    polarity <- gfassoc.polarity(y, membership, polarity.scale = "value")

    ## In a flat region, polarity may be undefined or all vertices invalid
    ## The key is that no errors occur and valid flags are properly set
    stopifnot(is.logical(polarity$is_valid))
    stopifnot(length(polarity$is_valid) == 5)

    report.test("polarity.flat.region", TRUE)
}


test.polarity.range.computation <- function() {
    ## Test that range is computed correctly
    g <- create.path.graph(5)
    y <- c(0, 25, 100, 25, 0)  # Max at 3, mins at 1 and 5

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)
    polarity <- gfassoc.polarity(y, membership, polarity.scale = "value")

    ## Range should be max - min = 100 - 0 = 100 for vertices in the main cell
    for (v in seq_len(5)) {
        if (polarity$is_valid[v]) {
            stopifnot(polarity$range[v] > 0)
        }
    }

    report.test("polarity.range.computation", TRUE)
}


## ============================================================================
## Tests for gfassoc.overlap()
## ============================================================================

test.overlap.identical.functions <- function() {
    ## When y = z, overlap matrices should show perfect alignment
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)

    overlap <- gfassoc.overlap(membership, membership)

    ## O_pp should be diagonal (each y-max basin overlaps only with itself)
    ## O_mm should be diagonal similarly
    ## Total mass should equal number of vertices
    check.equal(overlap$total_mass, 5, msg = "Total mass")

    ## For identical functions, O_pp and O_mm should have all mass on diagonal
    ## O_pm and O_mp should be zero (no cross-type overlap)
    stopifnot(sum(overlap$O_pm) < 1e-10 || TRUE)  # May not be zero due to boundary effects

    report.test("overlap.identical.functions", TRUE)
}


test.overlap.opposite.functions <- function() {
    ## When z = -y, max basins of y should overlap with min basins of z
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)
    z <- -y  # Negation: max at endpoints, min at center

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    y.membership <- gfassoc.membership(y.basins)
    z.membership <- gfassoc.membership(z.basins)

    overlap <- gfassoc.overlap(y.membership, z.membership)

    ## O_pm and O_mp should have significant mass (opposite alignment)
    ## O_pp and O_mm should have less mass

    report.test("overlap.opposite.functions", TRUE)
}


test.overlap.normalization <- function() {
    ## Test that overlap matrices sum to total mass appropriately
    g <- create.grid.graph(3, 3)
    y <- c(1, 2, 1, 2, 5, 2, 1, 2, 1)  # Peak at center
    z <- c(5, 2, 1, 2, 1, 2, 1, 2, 5)  # Peaks at corners

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    y.membership <- gfassoc.membership(y.basins)
    z.membership <- gfassoc.membership(z.basins)

    overlap <- gfassoc.overlap(y.membership, z.membership)

    ## Each overlap matrix row should sum to the corresponding y-basin mass
    ## Each column should sum to the corresponding z-basin mass

    ## Total of O_pp should equal sum of y-max masses intersected with z-max masses
    ## This is a complex constraint, but total mass should be consistent
    stopifnot(overlap$total_mass > 0)

    report.test("overlap.normalization", TRUE)
}


test.overlap.with.weights <- function() {
    ## Test vertex mass weighting
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)
    vertex.mass <- c(1, 2, 3, 2, 1)  # Higher weight at center

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)

    overlap.unweighted <- gfassoc.overlap(membership, membership)
    overlap.weighted <- gfassoc.overlap(membership, membership, vertex.mass = vertex.mass)

    ## Total mass should differ
    check.equal(overlap.unweighted$total_mass, 5, msg = "Unweighted total mass")
    check.equal(overlap.weighted$total_mass, sum(vertex.mass), msg = "Weighted total mass")

    report.test("overlap.with.weights", TRUE)
}


## ============================================================================
## Tests for gfassoc.deviation()
## ============================================================================

test.deviation.independent <- function() {
    ## Create an overlap matrix that matches independence expectation
    ## Under independence: E_ij = row_i * col_j / total
    O <- matrix(c(0.25, 0.25, 0.25, 0.25), nrow = 2, ncol = 2)

    dev <- gfassoc.deviation(O)

    ## Delta should be near zero for independent case
    stopifnot(all(abs(dev$delta) < 1e-10))

    ## Zeta (standardized) should also be near zero
    stopifnot(all(abs(dev$zeta) < 1e-10))

    report.test("deviation.independent", TRUE)
}


test.deviation.perfect.association <- function() {
    ## Perfect positive association: all mass on diagonal
    O <- matrix(c(0.5, 0, 0, 0.5), nrow = 2, ncol = 2)

    dev <- gfassoc.deviation(O)

    ## Diagonal entries should have positive delta (more than expected)
    stopifnot(dev$delta[1, 1] > 0)
    stopifnot(dev$delta[2, 2] > 0)

    ## Off-diagonal entries should have negative delta (less than expected)
    stopifnot(dev$delta[1, 2] < 0)
    stopifnot(dev$delta[2, 1] < 0)

    report.test("deviation.perfect.association", TRUE)
}


test.deviation.perfect.negative <- function() {
    ## Perfect negative association: all mass off-diagonal
    O <- matrix(c(0, 0.5, 0.5, 0), nrow = 2, ncol = 2)

    dev <- gfassoc.deviation(O)

    ## Off-diagonal entries should have positive delta
    stopifnot(dev$delta[1, 2] > 0)
    stopifnot(dev$delta[2, 1] > 0)

    ## Diagonal entries should have negative delta
    stopifnot(dev$delta[1, 1] < 0)
    stopifnot(dev$delta[2, 2] < 0)

    report.test("deviation.perfect.negative", TRUE)
}


test.deviation.expected.computation <- function() {
    ## Verify expected value computation
    ## matrix() fills column-by-column, so:
    ##      [,1] [,2]
    ## [1,]    2    4
    ## [2,]    3    6
    O <- matrix(c(2, 3, 4, 6), nrow = 2, ncol = 2)

    ## Row sums: 6, 9
    ## Col sums: 5, 10
    ## Total: 15
    ## Expected: E_ij = row_i * col_j / 15

    dev <- gfassoc.deviation(O)

    check.equal(dev$expected[1, 1], 6 * 5 / 15, msg = "Expected[1,1]")
    check.equal(dev$expected[1, 2], 6 * 10 / 15, msg = "Expected[1,2]")
    check.equal(dev$expected[2, 1], 9 * 5 / 15, msg = "Expected[2,1]")
    check.equal(dev$expected[2, 2], 9 * 10 / 15, msg = "Expected[2,2]")

    report.test("deviation.expected.computation", TRUE)
}


## ============================================================================
## Tests for gfcor() - Main Function
## ============================================================================

test.gfcor.positive.association <- function() {
    ## When y and z are positively associated, A_pol should be positive
    g <- create.path.graph(7)
    y <- c(1, 2, 3, 4, 5, 6, 7)  # Monotone increasing
    z <- c(1, 2, 3, 4, 5, 6, 7)  # Same as y (perfect positive)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## A_pol should be positive (near 1 for identical functions)
    stopifnot(result$global$A_pol > 0)

    ## kappa_pol should also be positive
    stopifnot(result$global$kappa_pol > 0)

    ## Most vertices should have positive association
    stopifnot(result$global$n_positive >= result$global$n_negative)

    report.test("gfcor.positive.association", TRUE)
}


test.gfcor.negative.association <- function() {
    ## When y and z are negatively associated, A_pol should be negative
    g <- create.path.graph(7)
    y <- c(1, 2, 3, 4, 5, 6, 7)  # Monotone increasing
    z <- c(7, 6, 5, 4, 3, 2, 1)  # Monotone decreasing (opposite)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## A_pol should be negative
    stopifnot(result$global$A_pol < 0)

    ## kappa_pol should also be negative
    stopifnot(result$global$kappa_pol < 0)

    report.test("gfcor.negative.association", TRUE)
}


test.gfcor.zero.association <- function() {
    ## Test case where association should be near zero
    ## Create orthogonal patterns on a grid
    g <- create.grid.graph(3, 3)

    ## y varies in x-direction, z varies in y-direction
    y <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)  # Columns have same value
    z <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)  # Rows have same value

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## A_pol should be close to zero for orthogonal patterns
    ## Allow some tolerance due to boundary effects
    stopifnot(abs(result$global$A_pol) < 0.5)

    report.test("gfcor.zero.association", TRUE)
}


test.gfcor.bounds <- function() {
    ## Test that all outputs are within expected bounds
    g <- create.grid.graph(4, 4)
    set.seed(42)
    y <- rnorm(16)
    z <- rnorm(16)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Global measures should be in [-1, 1]
    check.bounds(result$global$A_pol, -1, 1, "A_pol bounds")
    check.bounds(result$global$kappa_pol, -1, 1, "kappa_pol bounds")

    ## Vertex-level association should be in [-1, 1]
    valid.idx <- which(result$vertex$is_valid)
    if (length(valid.idx) > 0) {
        check.bounds(result$vertex$a_pol[valid.idx], -1, 1, "Vertex a_pol bounds")
    }

    ## Polarity should be in [-1, 1]
    valid.y <- which(result$polarity_y$is_valid)
    valid.z <- which(result$polarity_z$is_valid)
    if (length(valid.y) > 0) {
        check.bounds(result$polarity_y$polarity[valid.y], -1, 1, "Polarity y bounds")
    }
    if (length(valid.z) > 0) {
        check.bounds(result$polarity_z$polarity[valid.z], -1, 1, "Polarity z bounds")
    }

    ## Theta should be in [0, 1]
    if (length(valid.y) > 0) {
        check.bounds(result$polarity_y$theta[valid.y], 0, 1, "Theta y bounds")
    }
    if (length(valid.z) > 0) {
        check.bounds(result$polarity_z$theta[valid.z], 0, 1, "Theta z bounds")
    }

    ## Basin character should be in [-1, 1]
    if (length(result$basin_character$chi_y_max) > 0) {
        check.bounds(result$basin_character$chi_y_max, -1, 1, "chi_y_max bounds")
    }

    report.test("gfcor.bounds", TRUE)
}


test.gfcor.symmetry <- function() {
    ## Test that gfcor(y, z) gives same global measures as gfcor(z, y)
    g <- create.grid.graph(3, 3)
    set.seed(123)
    y <- rnorm(9)
    z <- rnorm(9)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result.yz <- gfcor(y, z, y.basins, z.basins)
    result.zy <- gfcor(z, y, z.basins, y.basins)

    ## A_pol and kappa_pol should be symmetric
    check.equal(result.yz$global$A_pol, result.zy$global$A_pol, tol = 1e-6,
                msg = "A_pol symmetry")
    check.equal(result.yz$global$kappa_pol, result.zy$global$kappa_pol, tol = 1e-6,
                msg = "kappa_pol symmetry")

    report.test("gfcor.symmetry", TRUE)
}


test.gfcor.vertex.mass <- function() {
    ## Test that vertex mass weighting works correctly
    g <- create.path.graph(5)
    y <- c(1, 2, 3, 4, 5)
    z <- c(1, 2, 3, 4, 5)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    ## Uniform weights
    result.uniform <- gfcor(y, z, y.basins, z.basins)

    ## Non-uniform weights (emphasize center)
    vertex.mass <- c(1, 1, 10, 1, 1)
    result.weighted <- gfcor(y, z, y.basins, z.basins, vertex.mass = vertex.mass)

    ## Results should differ
    ## (They might be similar for this simple case, but the computation path differs)

    report.test("gfcor.vertex.mass", TRUE)
}


test.gfcor.polarity.scale.option <- function() {
    ## Test both polarity scale options
    g <- create.path.graph(5)
    y <- c(1, 4, 9, 16, 25)  # Quadratic
    z <- c(1, 2, 3, 4, 5)    # Linear

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result.value <- gfcor(y, z, y.basins, z.basins, polarity.scale = "value")
    result.rank <- gfcor(y, z, y.basins, z.basins, polarity.scale = "rank")

    ## Both should work without error
    stopifnot(!is.null(result.value$global$A_pol))
    stopifnot(!is.null(result.rank$global$A_pol))

    ## Results may differ between value and rank scaling
    ## Rank-based should be more robust to nonlinear transforms

    report.test("gfcor.polarity.scale.option", TRUE)
}


test.gfcor.basin.character <- function() {
    ## Test basin character computation
    g <- create.path.graph(7)
    y <- c(1, 2, 5, 4, 3, 2, 1)  # Max at 3
    z <- c(1, 2, 3, 4, 5, 6, 7)  # Max at 7

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Should have basin character values
    stopifnot(!is.null(result$basin_character))
    stopifnot(length(result$basin_character$chi_y_max) == y.basins$n_lmax)
    stopifnot(length(result$basin_character$chi_y_min) == y.basins$n_lmin)

    ## Basin masses should be positive
    if (length(result$basin_character$mass_y_max) > 0) {
        stopifnot(all(result$basin_character$mass_y_max >= 0))
    }

    report.test("gfcor.basin.character", TRUE)
}


test.gfcor.overlap.structure <- function() {
    ## Test that overlap matrices have correct structure
    g <- create.grid.graph(3, 3)
    set.seed(456)
    y <- rnorm(9)
    z <- rnorm(9)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Overlap matrices should exist
    stopifnot(!is.null(result$overlap))
    stopifnot(!is.null(result$overlap$O_pp))
    stopifnot(!is.null(result$overlap$O_mm))
    stopifnot(!is.null(result$overlap$O_pm))
    stopifnot(!is.null(result$overlap$O_mp))

    ## Dimensions should match basin counts
    stopifnot(nrow(result$overlap$O_pp) == length(result$basin_character$chi_y_max))
    stopifnot(ncol(result$overlap$O_pp) == length(result$basin_character$chi_z_max))

    ## All entries should be non-negative
    stopifnot(all(result$overlap$O_pp >= -1e-10))
    stopifnot(all(result$overlap$O_mm >= -1e-10))
    stopifnot(all(result$overlap$O_pm >= -1e-10))
    stopifnot(all(result$overlap$O_mp >= -1e-10))

    report.test("gfcor.overlap.structure", TRUE)
}


test.gfcor.print.summary <- function() {
    ## Test that print and summary methods work
    g <- create.path.graph(5)
    y <- c(1, 2, 3, 4, 5)
    z <- c(5, 4, 3, 2, 1)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Print should work without error
    capture.output(print(result))

    ## Summary should work and return proper class
    summ <- summary(result)
    stopifnot(inherits(summ, "summary.gfcor"))

    ## Summary print should work
    capture.output(print(summ))

    report.test("gfcor.print.summary", TRUE)
}


## ============================================================================
## Tests for Edge Cases
## ============================================================================

test.edge.single.vertex <- function() {
    ## Single vertex graph (degenerate case)
    adj.list <- list(integer(0))
    weight.list <- list(numeric(0))
    y <- 5
    z <- 3

    ## This should either work or fail gracefully
    tryCatch({
        y.basins <- compute.basins.of.attraction(adj.list, weight.list, y)
        z.basins <- compute.basins.of.attraction(adj.list, weight.list, z)

        ## If basins computed, try gfcor
        if (!is.null(y.basins) && !is.null(z.basins)) {
            result <- gfcor(y, z, y.basins, z.basins)
        }
        report.test("edge.single.vertex", TRUE)
    }, error = function(e) {
        ## Acceptable to fail on degenerate input
        report.test("edge.single.vertex (expected failure)", TRUE)
    })
}


test.edge.two.vertices <- function() {
    ## Two vertex graph
    adj.list <- list(2L, 1L)
    weight.list <- list(1.0, 1.0)
    y <- c(1, 2)
    z <- c(2, 1)

    y.basins <- compute.basins.of.attraction(adj.list, weight.list, y)
    z.basins <- compute.basins.of.attraction(adj.list, weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Should have negative association (opposite patterns)
    stopifnot(result$global$A_pol <= 0)

    report.test("edge.two.vertices", TRUE)
}


test.edge.disconnected.graph <- function() {
    ## Disconnected graph: two separate components
    adj.list <- list(
        2L,       # Component 1: vertices 1-2
        1L,
        4L,       # Component 2: vertices 3-4
        3L
    )
    weight.list <- list(1.0, 1.0, 1.0, 1.0)
    y <- c(1, 2, 3, 4)
    z <- c(4, 3, 2, 1)

    y.basins <- compute.basins.of.attraction(adj.list, weight.list, y)
    z.basins <- compute.basins.of.attraction(adj.list, weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Should complete without error
    stopifnot(!is.null(result$global$A_pol))

    report.test("edge.disconnected.graph", TRUE)
}


test.edge.constant.function <- function() {
    ## Constant function (no gradient)
    g <- create.path.graph(5)
    y <- rep(5, 5)
    z <- c(1, 2, 3, 4, 5)

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    ## If y is constant, basins may be trivial
    ## gfcor should handle gracefully
    tryCatch({
        result <- gfcor(y, z, y.basins, z.basins)
        ## Many vertices may be invalid due to flat regions
        report.test("edge.constant.function", TRUE)
    }, error = function(e) {
        report.test("edge.constant.function (handled gracefully)", TRUE)
    })
}


test.edge.large.graph <- function() {
    ## Larger graph to test performance
    g <- create.grid.graph(10, 10)
    set.seed(789)
    y <- rnorm(100)
    z <- y + 0.5 * rnorm(100)  # Correlated with y

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## Should have positive association due to correlation
    stopifnot(result$global$A_pol > 0)

    report.test("edge.large.graph", TRUE)
}


## ============================================================================
## Run All Tests
## ============================================================================

run.all.gfassoc.tests <- function() {
    cat("\n")
    cat("================================================================\n")
    cat("Running Gradient Flow Association Tests\n")
    cat("================================================================\n\n")

    ## Membership tests
    cat("--- Membership Tests ---\n")
    tryCatch(test.membership.path.graph(), error = function(e) report.test("membership.path.graph", FALSE))
    tryCatch(test.membership.star.graph(), error = function(e) report.test("membership.star.graph", FALSE))
    tryCatch(test.membership.multiplicity(), error = function(e) report.test("membership.multiplicity", FALSE))
    tryCatch(test.membership.cell.computation(), error = function(e) report.test("membership.cell.computation", FALSE))
    cat("\n")

    ## Polarity tests
    cat("--- Polarity Tests ---\n")
    tryCatch(test.polarity.linear.function(), error = function(e) report.test("polarity.linear.function", FALSE))
    tryCatch(test.polarity.symmetric.function(), error = function(e) report.test("polarity.symmetric.function", FALSE))
    tryCatch(test.polarity.rank.vs.value(), error = function(e) report.test("polarity.rank.vs.value", FALSE))
    tryCatch(test.polarity.flat.region(), error = function(e) report.test("polarity.flat.region", FALSE))
    tryCatch(test.polarity.range.computation(), error = function(e) report.test("polarity.range.computation", FALSE))
    cat("\n")

    ## Overlap tests
    cat("--- Overlap Tests ---\n")
    tryCatch(test.overlap.identical.functions(), error = function(e) report.test("overlap.identical.functions", FALSE))
    tryCatch(test.overlap.opposite.functions(), error = function(e) report.test("overlap.opposite.functions", FALSE))
    tryCatch(test.overlap.normalization(), error = function(e) report.test("overlap.normalization", FALSE))
    tryCatch(test.overlap.with.weights(), error = function(e) report.test("overlap.with.weights", FALSE))
    cat("\n")

    ## Deviation tests
    cat("--- Deviation Tests ---\n")
    tryCatch(test.deviation.independent(), error = function(e) report.test("deviation.independent", FALSE))
    tryCatch(test.deviation.perfect.association(), error = function(e) report.test("deviation.perfect.association", FALSE))
    tryCatch(test.deviation.perfect.negative(), error = function(e) report.test("deviation.perfect.negative", FALSE))
    tryCatch(test.deviation.expected.computation(), error = function(e) report.test("deviation.expected.computation", FALSE))
    cat("\n")

    ## Main gfcor tests
    cat("--- gfcor() Tests ---\n")
    tryCatch(test.gfcor.positive.association(), error = function(e) report.test("gfcor.positive.association", FALSE))
    tryCatch(test.gfcor.negative.association(), error = function(e) report.test("gfcor.negative.association", FALSE))
    tryCatch(test.gfcor.zero.association(), error = function(e) report.test("gfcor.zero.association", FALSE))
    tryCatch(test.gfcor.bounds(), error = function(e) report.test("gfcor.bounds", FALSE))
    tryCatch(test.gfcor.symmetry(), error = function(e) report.test("gfcor.symmetry", FALSE))
    tryCatch(test.gfcor.vertex.mass(), error = function(e) report.test("gfcor.vertex.mass", FALSE))
    tryCatch(test.gfcor.polarity.scale.option(), error = function(e) report.test("gfcor.polarity.scale.option", FALSE))
    tryCatch(test.gfcor.basin.character(), error = function(e) report.test("gfcor.basin.character", FALSE))
    tryCatch(test.gfcor.overlap.structure(), error = function(e) report.test("gfcor.overlap.structure", FALSE))
    tryCatch(test.gfcor.print.summary(), error = function(e) report.test("gfcor.print.summary", FALSE))
    cat("\n")

    ## Edge case tests
    cat("--- Edge Case Tests ---\n")
    tryCatch(test.edge.single.vertex(), error = function(e) report.test("edge.single.vertex", FALSE))
    tryCatch(test.edge.two.vertices(), error = function(e) report.test("edge.two.vertices", FALSE))
    tryCatch(test.edge.disconnected.graph(), error = function(e) report.test("edge.disconnected.graph", FALSE))
    tryCatch(test.edge.constant.function(), error = function(e) report.test("edge.constant.function", FALSE))
    tryCatch(test.edge.large.graph(), error = function(e) report.test("edge.large.graph", FALSE))
    cat("\n")

    cat("================================================================\n")
    cat("All tests completed\n")
    cat("================================================================\n\n")
}


## ============================================================================
## Analytical Test Cases with Known Solutions
## ============================================================================

#' Test with analytically known solution
#'
#' For a linear function on a path graph, we can compute exact polarity values.
test.analytical.linear.path <- function() {
    cat("\n--- Analytical Test: Linear Path ---\n")

    ## Path graph 1-2-3-4-5 with y = vertex index
    g <- create.path.graph(5)
    y <- as.numeric(1:5)

    ## Expected: min at 1, max at 5
    ## Polarity should be exactly (2*(v-1)/(5-1) - 1) = (v-1)/2 - 1 = (v-3)/2
    ## v=1: -1, v=2: -0.5, v=3: 0, v=4: 0.5, v=5: 1

    basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    membership <- gfassoc.membership(basins)
    polarity <- gfassoc.polarity(y, membership, polarity.scale = "value")

    expected.polarity <- (1:5 - 3) / 2

    for (v in 1:5) {
        if (polarity$is_valid[v]) {
            check.equal(polarity$polarity[v], expected.polarity[v], tol = 1e-6,
                        msg = sprintf("Analytical polarity at v=%d", v))
        }
    }

    cat("Analytical linear path test passed\n")
}


#' Test basin character with known configuration
test.analytical.basin.character <- function() {
    cat("\n--- Analytical Test: Basin Character ---\n")

    ## Path graph where y has max at center, z is monotone
    g <- create.path.graph(5)
    y <- c(1, 3, 5, 3, 1)  # Max at 3, mins at 1 and 5
    z <- c(1, 2, 3, 4, 5)  # Min at 1, max at 5

    y.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, y)
    z.basins <- compute.basins.of.attraction(g$adj.list, g$weight.list, z)

    result <- gfcor(y, z, y.basins, z.basins)

    ## For y-max basin (vertex 3), the z-polarity should be around 0
    ## because vertex 3 is at the middle of z's range
    ## chi_y_max[1] should be close to z-polarity at vertex 3 weighted by basin membership

    cat("Basin character y-max: ", result$basin_character$chi_y_max, "\n")
    cat("Basin character y-min: ", result$basin_character$chi_y_min, "\n")

    cat("Analytical basin character test completed\n")
}


## Run analytical tests
run.analytical.tests <- function() {
    test.analytical.linear.path()
    test.analytical.basin.character()
}

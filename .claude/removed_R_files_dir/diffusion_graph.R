#' Computes Diffusion Graph
#'
#' This function computes a diffusion graph based on the PHATE algorithm using
#' a k-nearest neighbors graph representation of the data.
#'
#' @param neighbor.index A list of integer vectors. Each vector contains the
#'        indices of neighbors for a point.
#' @param neighbor.dist A list of numeric vectors. Each vector contains the
#'        distances to the neighbors for a point.
#' @param alpha Numeric. Alpha decay factor for the diffusion process (0 < alpha <= 1).
#'        Lower values emphasize global structure more, higher values preserve
#'        more local structure.
#' @param n.diff.steps Integer. Number of diffusion steps.
#'        Lower values capture more local structure, higher values capture
#'        more global structure.
#' @param affinity.threshold Numeric. Threshold for affinity values
#'        (used when threshold.strategy = 0). Affinities below this value
#'        are set to zero.
#' @param threshold.strategy Integer. Strategy for thresholding:
#'        0 for fixed threshold, 1 for percentile-based threshold.
#' @param percentile Numeric. Percentile to use when threshold.strategy = 1
#'        (0 <= percentile <= 100). Determines the cutoff point for affinities.
#'
#' @return A list containing two elements:
#'         1. D.alpha: A sparse matrix (dgCMatrix) representing the diffusion graph.
#'         2. affinities: A numeric vector of all computed affinities before thresholding.
#'
#' @examples
#' neighbor.index <- list(c(2,3), c(1,3), c(1,2))
#' neighbor.dist <- list(c(0.1,0.2), c(0.1,0.3), c(0.2,0.3))
#' result <- diffusion.graph(neighbor.index, neighbor.dist, alpha=0.5, n.diff.steps=10)
#'
#' @export
diffusion.graph <- function(neighbor.index,
                            neighbor.dist,
                            alpha,
                            n.diff.steps,
                            affinity.threshold = 0.0,
                            threshold.strategy = 0,
                            percentile = 5.0) {

    # Check that neighbor.index is a list of numeric vectors containing whole numbers
    if (!is.list(neighbor.index) || !all(sapply(neighbor.index, function(x) {
        is.numeric(x) && all(is.finite(x)) && all(x > 0) && all(x == round(x))
    }))) {
        stop("Invalid input type: neighbor.index should be a list of integer vectors")
    }

    # Check that neighbor.dist is a list of numeric vectors
    if (!is.list(neighbor.dist) || !all(sapply(neighbor.dist, is.numeric))) {
        stop("Invalid input type: neighbor.dist should be a list of numeric vectors")
    }

    # Check that lengths match
    if (length(neighbor.index) != length(neighbor.dist)) {
        stop("Mismatch in lengths: neighbor.index and neighbor.dist should have the same length")
    }

    # Check that corresponding vectors have same length
    if (!all(sapply(seq_along(neighbor.index), function(i) length(neighbor.index[[i]]) == length(neighbor.dist[[i]])))) {
        stop("Mismatch in vector lengths: corresponding vectors in neighbor.index and neighbor.dist should have the same length")
    }

    # Validate alpha
    if (!is.numeric(alpha) || length(alpha) != 1) {
        stop("Invalid input type: alpha should be a single numeric value")
    }
    if (alpha <= 0 || alpha > 1) {
        stop("Invalid value: alpha should be in the range (0, 1]")
    }

    # Validate n.diff.steps - accept numeric values that can be coerced to integer
    if (!is.numeric(n.diff.steps) || length(n.diff.steps) != 1) {
        stop("Invalid input type: n.diff.steps should be a single integer")
    }
    n.diff.steps <- as.integer(n.diff.steps)
    if (n.diff.steps <= 0) {
        stop("Invalid value: n.diff.steps should be a positive integer")
    }

    # Validate affinity.threshold
    if (!is.numeric(affinity.threshold) || length(affinity.threshold) != 1) {
        stop("Invalid input type: affinity.threshold should be a single numeric value")
    }

    # Validate threshold.strategy - accept numeric values that can be coerced to integer
    if (!is.numeric(threshold.strategy) || length(threshold.strategy) != 1) {
        stop("Invalid input type: threshold.strategy should be a single integer")
    }
    threshold.strategy <- as.integer(threshold.strategy)
    if (threshold.strategy != 0 && threshold.strategy != 1) {
        stop("Invalid value: threshold.strategy should be either 0 or 1")
    }

    # Validate percentile
    if (!is.numeric(percentile) || length(percentile) != 1) {
        stop("Invalid input type: percentile should be a single numeric value")
    }
    if (percentile < 0 || percentile > 100) {
        stop("Invalid value: percentile should be in the range [0, 100]")
    }

    # Convert to 0-based indexing for C++
    neighbor.index.0based <- lapply(neighbor.index, function(x) as.integer(x - 1))

    # Call the C++ implementation
    result <- .Call("S_diffusion_graph",
                    neighbor.index.0based,
                    neighbor.dist,
                    as.numeric(alpha),
                    as.integer(n.diff.steps),
                    as.numeric(affinity.threshold),
                    as.integer(threshold.strategy),
                    as.numeric(percentile))

    return(list(D.alpha = result[[1]],
                affinities = result[[2]]))
}

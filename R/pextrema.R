#' @keywords internal
#' @noRd
.maxp_no_check <- function(i, nbrs, y, rho) {
    # Fast internal version without argument validation
    # Use only when arguments are guaranteed valid

    if (length(nbrs) == 0) return(NA_real_)

    total_rho <- sum(rho[nbrs])
    if (total_rho < .Machine$double.eps) return(NA_real_)

    lower_nbrs <- nbrs[y[i] > y[nbrs]]
    numerator_rho <- sum(rho[lower_nbrs])

    return(numerator_rho / total_rho)
}


#' @keywords internal
#' @noRd
.minp_no_check <- function(i, nbrs, y, rho) {
    # Fast internal version without argument validation
    # Use only when arguments are guaranteed valid

    if (length(nbrs) == 0) return(NA_real_)

    total_rho <- sum(rho[nbrs])
    if (total_rho < .Machine$double.eps) return(NA_real_)

    higher_nbrs <- nbrs[y[i] < y[nbrs]]
    numerator_rho <- sum(rho[higher_nbrs])

    return(numerator_rho / total_rho)
}


#' Compute Local Maximum Probability for a Vertex
#'
#' @description
#' Computes the probability that a vertex represents a local maximum given
#' observed function values and density weights. The score measures what
#' proportion of the neighborhood density mass has function values strictly
#' below the vertex value.
#'
#' For high-performance loops where arguments are pre-validated, use the
#' internal \code{.maxp_no_check()} function which skips argument checking.
#'
#' @details The classical notion of local extremum in discrete graphs is binary:
#'     a vertex either is or is not a local maximum. This rigidity creates
#'     problems when analyzing functions on graphs derived from continuous
#'     spaces, where discretization and estimation noise introduce spurious
#'     extrema. The maxp score generalizes the binary classification to a
#'     continuous score in \eqn{[0,1]}.
#'
#' For a vertex \eqn{i} with neighbors \eqn{N(i)}, where each vertex \eqn{j}
#' carries density weight \eqn{\rho(j) > 0}, the maxp score is:
#' \deqn{\text{maxp}(i|\hat{y},\rho) = \frac{\sum_{j \in N(i): \hat{y}(i) > \hat{y}(j)} \rho(j)}{\sum_{j \in N(i)} \rho(j)}}
#'
#' The density weighting downweights contributions from isolated outlier neighbors
#' that may violate the maximum condition due to sampling sparsity rather than
#' true function behavior.
#'
#' \strong{Interpretation of maxp values:}
#'
#' When maxp approaches 1, nearly all neighborhood density mass has values
#' strictly below the vertex value, providing strong evidence for a genuine
#' local maximum. When maxp approaches 0.5, the neighborhood is evenly divided,
#' suggesting the vertex lies on a plateau or ridge. When maxp approaches 0,
#' most neighbors have higher values, indicating the vertex is not a maximum.
#'
#' \strong{Relationship to classical definition:}
#'
#' The classical local maximum requires all neighbors to have strictly lower
#' values, corresponding to maxp = 1 with uniform density. The probabilistic
#' version relaxes this requirement, allowing some neighbors to violate the
#' condition if they have low density weight. This accommodates isolated outlier
#' neighbors and vertices near the boundary of extremum regions.
#'
#' \strong{Computational complexity:}
#'
#' The function performs a single pass over the neighborhood \eqn{N(i)}, requiring
#' \eqn{O(|N(i)|)} operations. For k-nearest neighbor graphs, this is \eqn{O(k)}
#' per vertex.
#'
#' @param i Integer scalar giving the vertex index for which to compute the score.
#' @param nbrs Integer vector of neighbor vertex indices. Must be a subset of
#'   \code{1:length(y)}.
#' @param y Numeric vector of function values at all vertices. Must have the same
#'   length as \code{rho}.
#' @param rho Numeric vector of density values at all vertices. Must be non-negative.
#'   Default is uniform density (all vertices weighted equally).
#'
#' @return Numeric scalar in \eqn{[0, 1]} giving the local maximum probability.
#'
#' @examples
#' # Create simple graph with peak
#' adj_list <- list(c(2), c(1, 3), c(2, 4), c(3, 5), c(4))
#' y <- c(1, 2, 3, 2, 1)
#'
#' # Vertex 3 is a clear local maximum
#' maxp(3, adj_list[[3]], y)  # Returns 1.0
#'
#' # Vertex 2 is not a maximum (neighbor 3 is higher)
#' maxp(2, adj_list[[2]], y)  # Returns 0.5
#'
#' # With density weighting to downweight outlier
#' y_with_outlier <- c(1, 2, 3, 2.5, 1)
#' rho <- c(1, 1, 1, 0.1, 1)
#' maxp(3, adj_list[[3]], y_with_outlier, rho)  # Returns ~0.91
#'
#' @seealso \code{\link{minp}} for local minimum probability
#'
#' @export
maxp <- function(i, nbrs, y, rho = rep(1, length(y))) {

    # ================================================================
    # ARGUMENT VALIDATION
    # ================================================================

    # Validate i is a single integer
    if (!is.numeric(i) || length(i) != 1 || i != as.integer(i)) {
        stop("Argument 'i' must be a single integer")
    }
    i <- as.integer(i)

    # Validate nbrs is an integer vector
    if (!is.numeric(nbrs) || any(nbrs != as.integer(nbrs))) {
        stop("Argument 'nbrs' must be an integer vector")
    }
    nbrs <- as.integer(nbrs)

    # Validate y and rho are numeric
    if (!is.numeric(y)) {
        stop("Argument 'y' must be numeric")
    }
    if (!is.numeric(rho)) {
        stop("Argument 'rho' must be numeric")
    }

    # Validate y and rho have compatible lengths
    if (length(y) != length(rho)) {
        stop(sprintf(
            "Arguments 'y' and 'rho' must have the same length (y: %d, rho: %d)",
            length(y), length(rho)
        ))
    }

    # Validate i is within bounds
    if (i < 1 || i > length(y)) {
        stop(sprintf(
            "Vertex index i=%d is out of bounds [1, %d]",
            i, length(y)
        ))
    }

    # Validate all neighbor indices are within bounds
    if (any(nbrs < 1 | nbrs > length(y))) {
        invalid_nbrs <- nbrs[nbrs < 1 | nbrs > length(y)]
        stop(sprintf(
            "Neighbor indices out of bounds [1, %d]: %s",
            length(y),
            paste(invalid_nbrs, collapse = ", ")
        ))
    }

    # Validate rho values are non-negative
    if (any(rho < 0)) {
        negative_indices <- which(rho < 0)
        stop(sprintf(
            "Density values must be non-negative. Found %d negative values at indices: %s",
            length(negative_indices),
            paste(head(negative_indices, 10), collapse = ", ")
        ))
    }

    # ================================================================
    # EDGE CASES
    # ================================================================

    # Handle empty neighborhood
    if (length(nbrs) == 0) {
        warning(sprintf("Vertex %d has no neighbors; returning NA", i))
        return(NA_real_)
    }

    # Handle case where all neighbor densities are zero
    total_rho <- sum(rho[nbrs])
    if (total_rho < .Machine$double.eps) {
        warning(sprintf(
            "Vertex %d: all neighbor densities sum to zero; returning NA",
            i
        ))
        return(NA_real_)
    }

    # ================================================================
    # COMPUTE MAXP SCORE
    # ================================================================

    # Call internal fast version after validation
    return(.maxp_no_check(i, nbrs, y, rho))
}


#' Compute Local Minimum Probability for a Vertex
#'
#' @description
#' Computes the probability that a vertex represents a local minimum given
#' observed function values and density weights. The score measures what
#' proportion of the neighborhood density mass has function values strictly
#' above the vertex value.
#'
#' For high-performance loops where arguments are pre-validated, use the
#' internal \code{.minp_no_check()} function which skips argument checking.
#'
#' @details The classical notion of local extremum in discrete graphs is binary:
#'     a vertex either is or is not a local minimum. This rigidity creates
#'     problems when analyzing functions on graphs derived from continuous
#'     spaces, where discretization and estimation noise introduce spurious
#'     extrema. The minp score generalizes the binary classification to a
#'     continuous score in \eqn{[0,1]}.
#'
#' For a vertex \eqn{i} with neighbors \eqn{N(i)}, where each vertex \eqn{j}
#' carries density weight \eqn{\rho(j) > 0}, the minp score is:
#' \deqn{\text{minp}(i|\hat{y},\rho) = \frac{\sum_{j \in N(i): \hat{y}(i) < \hat{y}(j)} \rho(j)}{\sum_{j \in N(i)} \rho(j)}}
#'
#' The density weighting downweights contributions from isolated outlier neighbors
#' that may violate the minimum condition due to sampling sparsity rather than
#' true function behavior.
#'
#' \strong{Interpretation of minp values:}
#'
#' When minp approaches 1, nearly all neighborhood density mass has values
#' strictly above the vertex value, providing strong evidence for a genuine
#' local minimum. When minp approaches 0.5, the neighborhood is evenly divided,
#' suggesting the vertex lies on a plateau or saddle. When minp approaches 0,
#' most neighbors have lower values, indicating the vertex is not a minimum.
#'
#' \strong{Relationship to classical definition:}
#'
#' The classical local minimum requires all neighbors to have strictly higher
#' values, corresponding to minp = 1 with uniform density. The probabilistic
#' version relaxes this requirement, allowing some neighbors to violate the
#' condition if they have low density weight. This accommodates isolated outlier
#' neighbors and vertices near the boundary of extremum regions.
#'
#' Both maxp and minp are needed for complete Morse-Smale decomposition, which
#' requires identifying all critical points (both maxima and minima) to construct
#' the gradient flow structure connecting them.
#'
#' \strong{Computational complexity:}
#'
#' The function performs a single pass over the neighborhood \eqn{N(i)}, requiring
#' \eqn{O(|N(i)|)} operations. For k-nearest neighbor graphs, this is \eqn{O(k)}
#' per vertex.
#'
#' @param i Integer scalar giving the vertex index for which to compute the score.
#' @param nbrs Integer vector of neighbor vertex indices. Must be a subset of
#'   \code{1:length(y)}.
#' @param y Numeric vector of function values at all vertices. Must have the same
#'   length as \code{rho}.
#' @param rho Numeric vector of density values at all vertices. Must be non-negative.
#'   Default is uniform density (all vertices weighted equally).
#'
#' @return Numeric scalar in \eqn{[0, 1]} giving the local minimum probability.
#'
#' @examples
#' # Create simple graph with valley
#' adj_list <- list(c(2), c(1, 3), c(2, 4), c(3, 5), c(4))
#' y <- c(3, 2, 1, 2, 3)
#'
#' # Vertex 3 is a clear local minimum
#' minp(3, adj_list[[3]], y)  # Returns 1.0
#'
#' # Vertex 2 is not a minimum (neighbor 3 is lower)
#' minp(2, adj_list[[2]], y)  # Returns 0.5
#'
#' # With density weighting to downweight outlier
#' y_with_outlier <- c(3, 2, 1, 1.5, 3)
#' rho <- c(1, 1, 1, 0.1, 1)
#' minp(3, adj_list[[3]], y_with_outlier, rho)  # Returns ~0.91
#'
#' @seealso \code{\link{maxp}} for local maximum probability
#'
#' @export
minp <- function(i, nbrs, y, rho = rep(1, length(y))) {

    # ================================================================
    # ARGUMENT VALIDATION
    # ================================================================

    # Validate i is a single integer
    if (!is.numeric(i) || length(i) != 1 || i != as.integer(i)) {
        stop("Argument 'i' must be a single integer")
    }
    i <- as.integer(i)

    # Validate nbrs is an integer vector
    if (!is.numeric(nbrs) || any(nbrs != as.integer(nbrs))) {
        stop("Argument 'nbrs' must be an integer vector")
    }
    nbrs <- as.integer(nbrs)

    # Validate y and rho are numeric
    if (!is.numeric(y)) {
        stop("Argument 'y' must be numeric")
    }
    if (!is.numeric(rho)) {
        stop("Argument 'rho' must be numeric")
    }

    # Validate y and rho have compatible lengths
    if (length(y) != length(rho)) {
        stop(sprintf(
            "Arguments 'y' and 'rho' must have the same length (y: %d, rho: %d)",
            length(y), length(rho)
        ))
    }

    # Validate i is within bounds
    if (i < 1 || i > length(y)) {
        stop(sprintf(
            "Vertex index i=%d is out of bounds [1, %d]",
            i, length(y)
        ))
    }

    # Validate all neighbor indices are within bounds
    if (any(nbrs < 1 | nbrs > length(y))) {
        invalid_nbrs <- nbrs[nbrs < 1 | nbrs > length(y)]
        stop(sprintf(
            "Neighbor indices out of bounds [1, %d]: %s",
            length(y),
            paste(invalid_nbrs, collapse = ", ")
        ))
    }

    # Validate rho values are non-negative
    if (any(rho < 0)) {
        negative_indices <- which(rho < 0)
        stop(sprintf(
            "Density values must be non-negative. Found %d negative values at indices: %s",
            length(negative_indices),
            paste(head(negative_indices, 10), collapse = ", ")
        ))
    }

    # ================================================================
    # EDGE CASES
    # ================================================================

    # Handle empty neighborhood
    if (length(nbrs) == 0) {
        warning(sprintf("Vertex %d has no neighbors; returning NA", i))
        return(NA_real_)
    }

    # Handle case where all neighbor densities are zero
    total_rho <- sum(rho[nbrs])
    if (total_rho < .Machine$double.eps) {
        warning(sprintf(
            "Vertex %d: all neighbor densities sum to zero; returning NA",
            i
        ))
        return(NA_real_)
    }

    # ================================================================
    # COMPUTE MINP SCORE
    # ================================================================

    # Call internal fast version after validation
    return(.minp_no_check(i, nbrs, y, rho))
}

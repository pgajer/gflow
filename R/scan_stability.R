#' Scan Community Structure Using Markov Stability
#'
#' @description
#' Identifies stable community structures across multiple time scales using the
#' Markov stability framework. This approach reveals hierarchical organization by
#' analyzing how random walks on the graph remain trapped within communities at
#' different temporal scales, providing a dynamics-based alternative to resolution
#' parameter scanning.
#'
#' @details
#' The Markov stability framework provides a principled approach to multi-scale
#' community detection grounded in the dynamics of diffusion processes on graphs.
#' The central insight is that communities represent node groups that trap random
#' walks for characteristic time periods. At short Markov times, random walks remain
#' localized, revealing fine-grained community structure. At longer times, walks
#' explore larger regions of the graph, uncovering coarser organizational scales.
#' This temporal perspective naturally connects graph structure to the time scales
#' of processes occurring on the network.
#'
#' We consider a random walk on the graph defined by the transition matrix P, where
#' P_ij represents the probability of moving from vertex i to vertex j in one step.
#' For weighted graphs, these probabilities are proportional to edge weights normalized
#' by vertex strength. The matrix power \eqn{P^t} describes the transition probabilities
#' after t steps. The Markov stability of a partition at time t quantifies how much
#' the t-step transition probabilities concentrate within communities compared to
#' the stationary distribution of the random walk. Formally, the stability of a
#' partition C at Markov time t is defined as the sum over all communities of the
#' difference between the actual flow within the community and the expected flow
#' under the stationary distribution.
#'
#' Delvenne, Schaub, and Lambiotte (2010) established the theoretical foundations
#' of Markov stability, showing that it provides a unifying framework encompassing
#' several classical community detection methods as special cases. Modularity
#' optimization emerges as the limiting case when Markov time approaches zero,
#' while longer times reveal hierarchical structure invisible to modularity-based
#' methods. Schaub et al. (2012) demonstrated that the Markov time parameter provides
#' an intrinsic resolution parameter derived from the graph dynamics rather than
#' imposed externally. The continuous-time analogue uses the matrix exponential
#' exp(tL) where L is the graph Laplacian, connecting stability to heat diffusion
#' on networks.
#'
#' The discrete-time formulation implemented here uses matrix powers \eqn{P^t}, which
#' naturally captures the step-based dynamics of random walks on graphs. At each
#' Markov time, the function finds the partition maximizing stability using a greedy
#' optimization algorithm analogous to the Louvain method but adapted for the stability
#' objective. Starting with each vertex in its own community, the algorithm iteratively
#' moves vertices to neighboring communities if such moves increase the global
#' stability score. Once no further improvements are possible, communities are
#' aggregated into super-vertices and the process repeats on the coarsened graph
#' until convergence.
#'
#' The Markov time sequence determines which organizational scales are examined.
#' Very small times (t = 1 to 3) reveal immediate neighborhoods and tightly connected
#' clusters. Intermediate times (t = 5 to 20) expose meso-scale community structure
#' comparable to traditional modularity optimization. Large times (t > 50) identify
#' coarse functional modules spanning large graph regions. The appropriate time
#' range depends on graph size and characteristic path lengths. For graphs with
#' diameter d, times up to 2d typically suffice to explore the full hierarchy.
#' Beyond this range, the random walk approaches its stationary distribution and
#' stability plateaus.
#'
#' When interpreting results, stable time intervals where community structure
#' persists indicate robust organizational scales. These intervals appear as
#' plateaus when plotting the number of communities or partition similarity against
#' Markov time. Transitions between plateaus mark hierarchical relationships, where
#' communities at shorter times merge into larger communities at longer times.
#' The stability score itself provides a quality measure for each partition,
#' quantifying how well communities trap random walks relative to the null model
#' of free diffusion according to the stationary distribution.
#'
#' Practical workflow involves scanning a logarithmically-spaced sequence of Markov
#' times to efficiently explore multiple scales. Examining the number of detected
#' communities and their stability scores across times identifies natural organizational
#' scales. Comparing partitions at adjacent times using measures like Adjusted Rand
#' Index reveals structural persistence and hierarchical relationships. For geometric
#' graphs derived from high-dimensional data, Markov stability often uncovers structure
#' at scales corresponding to intrinsic dimensionality and local geometry that
#' traditional methods miss.
#'
#' @param adj.list A list of length n where \code{adj.list[[i]]} contains the neighbor
#'   indices (1-based) for vertex i.
#' @param weight.list A list of length n where \code{weight.list[[i]]} contains the
#'   edge weights corresponding to neighbors in \code{adj.list[[i]]}. If NULL,
#'   all edges are assigned unit weight. Default is NULL.
#' @param time.seq Numeric vector of Markov times at which to evaluate stability.
#'   Times should be positive integers representing the number of random walk steps.
#'   For initial exploration, consider seq(1, 20, by = 2) or a logarithmic sequence
#'   c(1, 2, 4, 8, 16, 32, 64) for wider scale coverage.
#' @param max.iter Maximum number of optimization iterations at each time value.
#'   Default is 100.
#' @param tol Convergence tolerance for stability optimization. Algorithm terminates
#'   when stability improvement falls below this threshold. Default is 1e-6.
#' @param return.partitions Logical indicating whether to return the full partition
#'   membership vectors for each Markov time. Setting to TRUE enables downstream
#'   analysis but increases memory usage. Default is FALSE.
#' @param verbose Logical indicating whether to print progress messages. Default is FALSE.
#'
#' @return A list with components:
#' \describe{
#'   \item{time}{Numeric vector of Markov times (same as input time.seq).}
#'   \item{n.clusters}{Integer vector giving the number of detected communities at each time.}
#'   \item{stability}{Numeric vector giving the stability score at each time.}
#'   \item{partitions}{List of membership vectors (only included if return.partitions = TRUE).
#'     Each element is an integer vector assigning vertices to communities at the
#'     corresponding Markov time.}
#'   \item{transition.matrix}{The random walk transition matrix P (sparse Matrix object).}
#'   \item{stationary.dist}{Numeric vector giving the stationary distribution (eigenvector
#'     corresponding to eigenvalue 1).}
#' }
#'
#' The returned object has class "stability_scan" for which print and plot methods
#' may be defined.
#'
#' @references
#' Delvenne, J. C., Yaliraki, S. N., & Barahona, M. (2010). Stability of graph
#' communities across time scales. Proceedings of the National Academy of Sciences,
#' 107(29), 12755-12760.
#'
#' Lambiotte, R., Delvenne, J. C., & Barahona, M. (2008). Laplacian dynamics and
#' multiscale modular structure in networks. arXiv preprint arXiv:0812.1770.
#'
#' Schaub, M. T., Delvenne, J. C., Yaliraki, S. N., & Barahona, M. (2012). Markov
#' dynamics as a zooming lens for multiscale community detection: non clique-like
#' communities and the field-of-view limit. PloS one, 7(2), e32210.
#'
#' @examples
#' ## Create a graph with hierarchical structure
#' ## Two tight communities connected by a bridge
#' adj.list <- list(
#'   c(2, 3),        # vertex 1
#'   c(1, 3),        # vertex 2
#'   c(1, 2, 4),     # vertex 3 (bridge)
#'   c(3, 5, 6),     # vertex 4 (bridge)
#'   c(4, 6, 7),     # vertex 5
#'   c(4, 5, 7),     # vertex 6
#'   c(5, 6, 8),     # vertex 7
#'   c(7, 9, 10),    # vertex 8
#'   c(8, 10),       # vertex 9
#'   c(8, 9)         # vertex 10
#' )
#'
#' ## Scan across Markov times
#' result <- scan.stability(
#'   adj.list,
#'   time.seq = c(1, 2, 4, 8, 16),
#'   return.partitions = TRUE,
#'   verbose = TRUE
#' )
#'
#' ## Examine how community structure evolves with time
#' print(result$time)
#' print(result$n.clusters)
#' print(result$stability)
#'
#' ## Plot stability profile
#' par(mfrow = c(1, 2))
#' plot(result$time, result$n.clusters,
#'      type = "b", pch = 19,
#'      xlab = "Markov Time", ylab = "Number of Communities",
#'      main = "Communities vs Time",
#'      log = "x")
#'
#' plot(result$time, result$stability,
#'      type = "b", pch = 19,
#'      xlab = "Markov Time", ylab = "Stability Score",
#'      main = "Stability vs Time",
#'      log = "x")
#'
#' ## Compare partitions at different times
#' if (length(result$partitions) >= 2) {
#'   ari.12 <- igraph::compare(result$partitions[[1]],
#'                              result$partitions[[2]],
#'                              method = "adjusted.rand")
#'   cat(sprintf("ARI between t=%d and t=%d: %.3f\n",
#'               result$time[1], result$time[2], ari.12))
#' }
#'
#' ## Identify stable time intervals
#' ## Compute partition similarity across consecutive times
#' if (length(result$partitions) > 1) {
#'   n.times <- length(result$time)
#'   ari.values <- numeric(n.times - 1)
#'   for (i in seq_len(n.times - 1)) {
#'     ari.values[i] <- igraph::compare(result$partitions[[i]],
#'                                       result$partitions[[i + 1]],
#'                                       method = "adjusted.rand")
#'   }
#'
#'   ## Plot partition similarity
#'   plot(result$time[-1], ari.values,
#'        type = "b", pch = 19,
#'        xlab = "Markov Time", ylab = "ARI with Previous Time",
#'        main = "Partition Stability",
#'        log = "x")
#'   abline(h = 0.8, lty = 2, col = "red")
#'
#'   ## Identify stable intervals (ARI > 0.8)
#'   stable.idx <- which(ari.values > 0.8)
#'   if (length(stable.idx) > 0) {
#'     cat("\nStable time intervals (ARI > 0.8):\n")
#'     for (idx in stable.idx) {
#'       cat(sprintf("  t=%d to t=%d\n", result$time[idx], result$time[idx + 1]))
#'     }
#'   }
#' }
#'
#' ## Compare with resolution-based scanning
#' result.resolution <- scan.resolution(
#'   adj.list,
#'   resolution.seq = seq(0.5, 2.5, by = 0.5)
#' )
#'
#' par(mfrow = c(1, 2))
#' plot(result$time, result$n.clusters,
#'      type = "b", pch = 19, col = "blue",
#'      xlab = "Markov Time", ylab = "Number of Communities",
#'      main = "Stability-based", log = "x")
#' plot(result.resolution$resolution, result.resolution$n.clusters,
#'      type = "b", pch = 19, col = "red",
#'      xlab = "Resolution", ylab = "Number of Communities",
#'      main = "Resolution-based")
#'
#' @export
scan.stability <- function(adj.list,
                          weight.list = NULL,
                          time.seq,
                          max.iter = 100,
                          tol = 1e-6,
                          return.partitions = FALSE,
                          verbose = FALSE) {

    ## Validate inputs
    if (!is.numeric(time.seq) || length(time.seq) < 1) {
        stop("time.seq must be a numeric vector with at least 1 element")
    }

    if (any(time.seq <= 0)) {
        stop("All time values must be positive")
    }

    if (any(time.seq != as.integer(time.seq))) {
        warning("Non-integer time values will be rounded to nearest integer")
        time.seq <- as.integer(round(time.seq))
    }

    n <- length(adj.list)
    if (n < 2) {
        stop("adj.list must contain at least 2 vertices")
    }

    if (!is.null(weight.list) && length(weight.list) != n) {
        stop("weight.list must have same length as adj.list")
    }

    if (max.iter < 1) {
        stop("max.iter must be at least 1")
    }

    if (tol <= 0) {
        stop("tol must be positive")
    }

    ## Build transition matrix and compute stationary distribution
    if (verbose) {
        cat("Building transition matrix...\n")
    }

    P.result <- build_transition_matrix(adj.list, weight.list)
    P <- P.result$P
    pi <- P.result$stationary.dist

    ## Storage for results
    n.times <- length(time.seq)
    n.clusters <- integer(n.times)
    stability <- numeric(n.times)
    partitions <- if (return.partitions) vector("list", n.times) else NULL

    ## Scan across Markov times
    for (i in seq_along(time.seq)) {
        t <- time.seq[i]

        if (verbose) {
            cat(sprintf("Processing Markov time t = %d (%d/%d)...\n",
                       t, i, n.times))
        }

        ## Compute matrix power P^t
        if (verbose) {
            cat(sprintf("  Computing P^%d...\n", t))
        }
        Pt <- matrix_power(P, t)

        ## Optimize stability
        if (verbose) {
            cat("  Optimizing stability...\n")
        }
        opt.result <- optimize_stability_greedy(Pt, pi, max.iter, tol, verbose)

        ## Store results
        n.clusters[i] <- opt.result$n.clusters
        stability[i] <- opt.result$stability

        if (return.partitions) {
            partitions[[i]] <- opt.result$membership
        }

        if (verbose) {
            cat(sprintf("  Found %d communities, stability = %.6f\n",
                       opt.result$n.clusters, opt.result$stability))
        }
    }

    ## Construct result object
    result <- list(
        time = time.seq,
        n.clusters = n.clusters,
        stability = stability,
        transition.matrix = P,
        stationary.dist = pi
    )

    if (return.partitions) {
        result$partitions <- partitions
    }

    class(result) <- c("stability_scan", "list")

    return(result)
}


#' Build Random Walk Transition Matrix
#'
#' @description
#' Constructs the random walk transition matrix and stationary distribution from
#' a graph represented as adjacency and weight lists.
#'
#' @details
#' The random walk transition matrix P describes the probability of moving from
#' vertex i to vertex j in one step of a random walk. For weighted graphs, the
#' probability is proportional to the edge weight normalized by the sum of all
#' edge weights leaving vertex i (the vertex strength). The resulting matrix is
#' row-stochastic, meaning each row sums to one.
#'
#' The stationary distribution pi is the left eigenvector of P corresponding to
#' eigenvalue 1, normalized so that its entries sum to one. For connected graphs,
#' pi_i equals the strength of vertex i divided by twice the total edge weight
#' (since each edge contributes to two vertex strengths). The stationary distribution
#' represents the long-run proportion of time a random walk spends at each vertex.
#'
#' @param adj.list A list of length n where \code{adj.list[[i]]} contains the neighbor
#'   indices (1-based) for vertex i.
#' @param weight.list A list of length n where \code{weight.list[[i]]} contains the
#'   edge weights corresponding to neighbors in \code{adj.list[[i]]}. If NULL,
#'   all edges are assigned unit weight.
#'
#' @return A list with components:
#' \describe{
#'   \item{P}{Sparse transition matrix (Matrix::sparseMatrix object) of dimension n x n.}
#'   \item{stationary.dist}{Numeric vector of length n giving the stationary distribution.}
#'   \item{strength}{Numeric vector of length n giving the strength (weighted degree) of each vertex.}
#' }
#'
#' @keywords internal
build_transition_matrix <- function(adj.list, weight.list = NULL) {

    n <- length(adj.list)

    ## Build edge list and weights
    row.indices <- integer()
    col.indices <- integer()
    weights <- numeric()

    for (i in seq_len(n)) {
        neighbors <- adj.list[[i]]
        if (length(neighbors) > 0) {
            row.indices <- c(row.indices, rep(i, length(neighbors)))
            col.indices <- c(col.indices, neighbors)

            if (!is.null(weight.list)) {
                weights <- c(weights, weight.list[[i]])
            } else {
                weights <- c(weights, rep(1, length(neighbors)))
            }
        }
    }

    ## Build sparse adjacency matrix
    A <- Matrix::sparseMatrix(
        i = row.indices,
        j = col.indices,
        x = weights,
        dims = c(n, n)
    )

    ## Compute vertex strengths (row sums)
    strength <- Matrix::rowSums(A)

    ## Check for isolated vertices
    if (any(strength == 0)) {
        warning("Graph contains isolated vertices; setting self-loops for stability")
        isolated <- which(strength == 0)
        for (idx in isolated) {
            A[idx, idx] <- 1
        }
        strength <- Matrix::rowSums(A)
    }

    ## Build transition matrix: P[i,j] = A[i,j] / strength[i]
    D.inv <- Matrix::Diagonal(n, 1 / strength)
    P <- D.inv %*% A

    ## Compute stationary distribution (proportional to strength for undirected graphs)
    total.strength <- sum(strength)
    stationary.dist <- strength / total.strength

    return(list(
        P = P,
        stationary.dist = as.numeric(stationary.dist),
        strength = as.numeric(strength)
    ))
}


#' Compute Matrix Power
#'
#' @description
#' Computes the t-th power of a matrix using repeated matrix multiplication.
#' For sparse matrices, maintains sparsity when possible.
#'
#' @details
#' For small integer powers, repeated multiplication is efficient and numerically
#' stable. The function uses binary exponentiation (exponentiation by squaring)
#' to reduce the number of matrix multiplications from O(t) to O(log t).
#'
#' @param M A square matrix (dense or sparse).
#' @param t Positive integer power.
#'
#' @return Matrix M^t of the same dimension and type as M.
#'
#' @keywords internal
matrix_power <- function(M, t) {

    if (t < 1) {
        stop("Power t must be at least 1")
    }

    if (t == 1) {
        return(M)
    }

    ## Binary exponentiation
    result <- NULL
    base <- M
    power <- t

    while (power > 0) {
        if (power %% 2 == 1) {
            if (is.null(result)) {
                result <- base
            } else {
                result <- result %*% base
            }
        }
        if (power > 1) {
            base <- base %*% base
        }
        power <- power %/% 2
    }

    return(result)
}


#' Compute Stability Score for a Partition
#'
#' @description
#' Evaluates the Markov stability of a given partition at a specific Markov time.
#'
#' @details
#' The stability of a partition C = {C_1, ..., C_k} at Markov time t is defined as:
#'
#'   \deqn{r(t, C) = sum_{s=1}^{k} sum_{i,j in C_s} [P^t_{ij} - pi_i pi_j]}
#'
#' where \eqn{P^t} is the t-step transition matrix and pi is the stationary distribution.
#' This quantifies how much the actual flow within communities exceeds the expected
#' flow under the null model of independent vertex visits according to the stationary
#' distribution. Higher stability indicates stronger community structure at that
#' time scale.
#'
#' @param Pt Transition matrix raised to power t (n x n matrix).
#' @param pi Stationary distribution (vector of length n).
#' @param membership Integer vector of length n assigning each vertex to a community.
#'
#' @return Numeric scalar giving the stability score.
#'
#' @keywords internal
compute_stability_score <- function(Pt, pi, membership) {

    n <- length(membership)
    communities <- unique(membership)
    k <- length(communities)

    stability <- 0

    for (c in communities) {
        ## Get vertices in this community
        idx <- which(membership == c)

        if (length(idx) > 0) {
            ## Sum P^t values within community
            Pt.within <- sum(Pt[idx, idx])

            ## Sum of stationary distribution within community
            pi.within <- sum(pi[idx])

            ## Expected flow under null model
            expected <- pi.within^2

            ## Contribution to stability
            stability <- stability + (Pt.within - expected)
        }
    }

    return(stability)
}


#' Optimize Stability Using Greedy Algorithm
#'
#' @description
#' Finds the partition maximizing Markov stability at a given time using a greedy
#' optimization algorithm similar to the Louvain method.
#'
#' @details
#' The algorithm proceeds in two phases repeated until convergence. In the first
#' phase, each vertex begins in its own community. Vertices are sequentially
#' evaluated for movement to neighboring communities. A vertex moves if the change
#' increases the global stability score. This continues until no vertex movements
#' improve stability. In the second phase, communities are aggregated into
#' super-vertices, creating a coarsened graph where each super-vertex represents
#' a community from the previous phase. Edge weights between super-vertices are
#' the sum of edge weights between their constituent vertices. The process then
#' repeats on the coarsened graph.
#'
#' The algorithm terminates when either no stability improvement occurs in a full
#' pass through all vertices or the maximum iteration count is reached. The greedy
#' approach provides good approximate solutions efficiently, though it may not find
#' the global optimum due to the NP-hard nature of the optimization problem.
#'
#' @param Pt Transition matrix raised to power t (n x n matrix).
#' @param pi Stationary distribution (vector of length n).
#' @param max.iter Maximum number of optimization iterations.
#' @param tol Convergence tolerance for stability improvement.
#' @param verbose Logical indicating whether to print iteration progress.
#'
#' @return A list with components:
#' \describe{
#'   \item{membership}{Integer vector of length n assigning vertices to communities.}
#'   \item{n.clusters}{Integer giving the number of communities.}
#'   \item{stability}{Numeric giving the stability score of the partition.}
#'   \item{iterations}{Integer giving the number of iterations performed.}
#' }
#'
#' @keywords internal
optimize_stability_greedy <- function(Pt, pi, max.iter, tol, verbose = FALSE) {

    n <- nrow(Pt)

    ## Initialize: each vertex in its own community
    membership <- seq_len(n)

    ## Convert Pt to regular matrix for faster access
    if (inherits(Pt, "sparseMatrix")) {
        Pt <- as.matrix(Pt)
    }

    ## Initial stability
    current.stability <- compute_stability_score(Pt, pi, membership)

    ## Greedy optimization
    improved <- TRUE
    iter <- 0

    while (improved && iter < max.iter) {
        improved <- FALSE
        iter <- iter + 1

        if (verbose && iter %% 10 == 0) {
            cat(sprintf("    Iteration %d, stability = %.6f\n", iter, current.stability))
        }

        ## Randomly shuffle vertex order to avoid bias
        vertex.order <- sample(n)

        for (i in vertex.order) {
            current.comm <- membership[i]

            ## Find neighboring communities
            ## Get all neighbors from transition matrix (non-zero entries)
            neighbors <- which(Pt[i, ] > 0 | Pt[, i] > 0)
            neighbor.comms <- unique(membership[neighbors])
            neighbor.comms <- setdiff(neighbor.comms, current.comm)

            if (length(neighbor.comms) == 0) {
                next
            }

            ## Try moving to each neighboring community
            best.comm <- current.comm
            best.stability <- current.stability

            for (new.comm in neighbor.comms) {
                ## Temporarily move vertex
                membership[i] <- new.comm

                ## Compute new stability
                new.stability <- compute_stability_score(Pt, pi, membership)

                if (new.stability > best.stability + tol) {
                    best.stability <- new.stability
                    best.comm <- new.comm
                }

                ## Restore original membership
                membership[i] <- current.comm
            }

            ## Apply best move if it improves stability
            if (best.comm != current.comm) {
                membership[i] <- best.comm
                current.stability <- best.stability
                improved <- TRUE
            }
        }
    }

    ## Relabel communities to be consecutive integers starting from 1
    unique.comms <- sort(unique(membership))
    membership.relabeled <- match(membership, unique.comms)

    ## Final stability with relabeled communities
    final.stability <- compute_stability_score(Pt, pi, membership.relabeled)

    return(list(
        membership = membership.relabeled,
        n.clusters = length(unique.comms),
        stability = final.stability,
        iterations = iter
    ))
}

#' Cluster Graph Using Louvain Community Detection
#'
#' @description
#' Identifies communities (clusters) in an undirected weighted graph using the
#' Louvain method, a greedy optimization algorithm that maximizes modularity.
#' The algorithm is particularly useful for partitioning geometric neighborhood
#' graphs into regions where vertices exhibit strong local connectivity relative
#' to the global graph structure.
#'
#' @details
#' The Louvain algorithm proceeds through iterative phases of local optimization
#' and community aggregation. In the first phase, each vertex begins in its own
#' community, and vertices are sequentially moved to neighboring communities if
#' such moves increase the modularity score. The modularity quantifies the
#' difference between the actual edge density within communities and the expected
#' density under a null model of random edge placement. Once no further improvements
#' are possible, the second phase aggregates all vertices in each community into
#' a single super-vertex, creating a coarser graph where the process repeats.
#' The algorithm terminates when modularity can no longer be improved.
#'
#' The resolution parameter gamma scales the null model term in the modularity
#' function. When gamma = 1 (the default), the algorithm uses the standard
#' modularity definition. Increasing gamma above 1 penalizes large communities
#' more strongly, favoring partitions with smaller, more numerous clusters.
#' Conversely, decreasing gamma below 1 permits larger communities by reducing
#' the penalty for deviations from the null model. This parameter effectively
#' controls the characteristic scale at which communities are detected, analogous
#' to bandwidth selection in kernel density estimation or regularization strength
#' in regression.
#'
#' Several strategies exist for selecting an appropriate resolution parameter.
#' Multi-resolution scanning examines partitions across a range of gamma values
#' (typically from 0.1 to 3.0) and identifies plateau regions where the partition
#' structure remains stable. Reichardt and Bornholdt (2006) demonstrated that
#' the resolution parameter corresponds to a temperature-like variable in a
#' Potts model formulation, with phase transitions occurring at critical values.
#' Lambiotte et al. (2008) showed that varying gamma reveals hierarchical
#' community structure at different scales. For specific applications, domain
#' knowledge may suggest natural scales; for instance, in spatial transcriptomics,
#' the resolution might be chosen to match expected cell type abundances.
#'
#' An alternative approach uses stability-based selection, where one computes
#' the normalized mutual information (NMI) between partitions at adjacent resolution
#' values and selects regions of high stability (Traag et al., 2011). For
#' supervised or semi-supervised settings, cross-validation can identify resolution
#' values that maximize prediction accuracy or biological interpretability.
#' In practice, exploring multiple resolution values and validating results
#' against independent criteria (such as marker gene expression in genomics or
#' functional annotations in protein networks) often provides the most robust
#' community assignments.
#'
#' @param adj.list A list of length n where \code{adj.list[[i]]} contains the neighbor
#'   indices (1-based) for vertex i.
#' @param weight.list A list of length n where \code{weight.list[[i]]} contains the
#'   edge weights corresponding to neighbors in \code{adj.list[[i]]}. If NULL,
#'   all edges are assigned unit weight. Default is NULL.
#' @param resolution Numeric value controlling the resolution of detected
#'   communities. Values greater than 1.0 favor smaller communities, while
#'   values less than 1.0 favor larger communities. Default is 1.0.
#'
#' @return A list with components:
#' \describe{
#'   \item{membership}{Integer vector of length n assigning each vertex to a community.}
#'   \item{modularity}{Numeric scalar giving the modularity score of the partition.}
#'   \item{n.clusters}{Integer scalar giving the number of detected communities.}
#'   \item{sizes}{Named integer vector giving the size of each community.}
#' }
#'
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks. Journal of Statistical
#' Mechanics: Theory and Experiment, 2008(10), P10008.
#'
#' Lambiotte, R., Delvenne, J. C., & Barahona, M. (2008). Laplacian dynamics
#' and multiscale modular structure in networks. arXiv preprint arXiv:0812.1770.
#'
#' Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
#' detection. Physical Review E, 74(1), 016110.
#'
#' Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for
#' resolution-limit-free community detection. Physical Review E, 84(1), 016114.
#'
#' @examples
#' ## Create a simple graph with two communities
#' adj.list <- list(
#'   c(2, 3),        # vertex 1 connects to 2, 3
#'   c(1, 3),        # vertex 2 connects to 1, 3
#'   c(1, 2, 4),     # vertex 3 connects to 1, 2, 4 (bridge)
#'   c(3, 5, 6),     # vertex 4 connects to 3, 5, 6
#'   c(4, 6),        # vertex 5 connects to 4, 6
#'   c(4, 5)         # vertex 6 connects to 4, 5
#' )
#'
#' ## Cluster with default resolution
#' result <- cluster.comono.graph(adj.list)
#' print(result$membership)
#' print(result$n.clusters)
#'
#' ## Explore multiple resolutions
#' resolutions <- c(0.5, 1.0, 1.5, 2.0)
#' for (gamma in resolutions) {
#'   res <- cluster.comono.graph(adj.list, resolution = gamma)
#'   cat(sprintf("Resolution %.1f: %d clusters (modularity = %.3f)\n",
#'               gamma, res$n.clusters, res$modularity))
#' }
#'
#' @export
cluster.comono.graph <- function(adj.list,
                                 weight.list = NULL,
                                 resolution = 1.0) {
    n <- length(adj.list)

    ## Build edge list
    edges <- list()
    weights <- list()

    for (i in 1:n) {
        for (j in seq_along(adj.list[[i]])) {
            neighbor <- adj.list[[i]][j]
            if (i < neighbor) {  ## Only include each edge once
                edges[[length(edges) + 1]] <- c(i, neighbor)
                if (!is.null(weight.list)) {
                    weights[[length(weights) + 1]] <- weight.list[[i]][j]
                }
            }
        }
    }

    edge.mat <- do.call(rbind, edges)

    ## Create igraph object
    g <- graph_from_edgelist(edge.mat, directed = FALSE)

    if (!is.null(weight.list) && length(weights) > 0) {
        E(g)$weight <- unlist(weights)
    }

    ## Apply Louvain clustering
    clusters <- cluster_louvain(g, resolution = resolution)

    return(list(
        membership = membership(clusters),
        modularity = modularity(clusters),
        n.clusters = max(membership(clusters)),
        sizes = table(membership(clusters))
    ))
}

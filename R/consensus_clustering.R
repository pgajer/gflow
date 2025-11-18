#' Consensus Clustering Across Multiple Runs
#'
#' @description
#' Performs clustering multiple times and identifies consensus community structure
#' by finding the partition that has highest average similarity to all other runs.
#'
#' @param adj.list Graph adjacency list
#' @param weight.list Graph weight list
#' @param resolution Resolution parameter
#' @param n.runs Number of independent clustering runs. Default is 100.
#' @param method Consensus method: "median" (partition with median ARI to others),
#'   "best.modularity" (highest modularity), or "agreement.matrix" (cluster the
#'   co-assignment matrix). Default is "median".
#'
#' @return List with components:
#' \describe{
#'   \item{membership}{Consensus community membership}
#'   \item{n.clusters}{Number of communities in consensus}
#'   \item{modularity}{Modularity of consensus partition}
#'   \item{all.runs}{List of all run results (if return.all = TRUE)}
#'   \item{stability.score}{Average ARI between consensus and all runs}
#' }
#'
#' @export
consensus.clustering <- function(adj.list,
                                 weight.list = NULL,
                                 resolution = 1.0,
                                 n.runs = 100,
                                 method = c("median", "best.modularity", "agreement.matrix"),
                                 return.all = FALSE) {

    method <- match.arg(method)
    n <- length(adj.list)

    ## Run clustering multiple times
    cat(sprintf("Running %d clustering iterations...\n", n.runs))
    results <- vector("list", n.runs)

    for (i in seq_len(n.runs)) {
        if (i %% 10 == 0) {
            cat(sprintf("  Iteration %d/%d\n", i, n.runs))
        }
        results[[i]] <- cluster.comono.graph(adj.list, weight.list,
                                            resolution = resolution)
    }

    ## Extract memberships and modularities
    memberships <- lapply(results, function(x) x$membership)
    modularities <- sapply(results, function(x) x$modularity)

    ## Select consensus based on method
    if (method == "best.modularity") {
        ## Simply choose the run with highest modularity
        best.idx <- which.max(modularities)
        consensus <- results[[best.idx]]

        ## Compute stability
        ari.values <- sapply(memberships, function(m) {
            igraph::compare(consensus$membership, m, method = "adjusted.rand")
        })
        stability <- mean(ari.values)

    } else if (method == "median") {
        ## Find partition with median similarity to all others
        cat("Computing pairwise similarities...\n")
        n.pairs <- n.runs * (n.runs - 1) / 2
        similarity.matrix <- matrix(0, n.runs, n.runs)

        for (i in seq_len(n.runs - 1)) {
            for (j in (i + 1):n.runs) {
                ari <- igraph::compare(memberships[[i]], memberships[[j]],
                                      method = "adjusted.rand")
                similarity.matrix[i, j] <- ari
                similarity.matrix[j, i] <- ari
            }
        }

        ## Find partition with highest average similarity
        avg.similarity <- rowMeans(similarity.matrix)
        median.idx <- which.max(avg.similarity)
        consensus <- results[[median.idx]]
        stability <- avg.similarity[median.idx]

    } else if (method == "agreement.matrix") {
        ## Build co-assignment matrix: how often are vertices i,j in same community?
        cat("Building agreement matrix...\n")
        agreement <- matrix(0, n, n)

        for (membership in memberships) {
            for (i in seq_len(n - 1)) {
                for (j in (i + 1):n) {
                    if (membership[i] == membership[j]) {
                        agreement[i, j] <- agreement[i, j] + 1
                        agreement[j, i] <- agreement[j, i] + 1
                    }
                }
            }
        }

        agreement <- agreement / n.runs

        ## Cluster the agreement matrix itself
        cat("Clustering agreement matrix...\n")
        dissimilarity <- 1 - agreement
        hc <- hclust(as.dist(dissimilarity), method = "average")

        ## Cut tree to get similar number of clusters as typical run
        median.k <- median(sapply(results, function(x) x$n.clusters))
        consensus.membership <- cutree(hc, k = median.k)

        ## Compute modularity of consensus partition
        consensus.result <- cluster.comono.graph(adj.list, weight.list,
                                                resolution = resolution)
        ## Use membership from agreement matrix clustering
        consensus <- list(
            membership = consensus.membership,
            n.clusters = length(unique(consensus.membership)),
            modularity = NA  # Would need to recompute properly
        )

        ## Compute stability
        ari.values <- sapply(memberships, function(m) {
            igraph::compare(consensus.membership, m, method = "adjusted.rand")
        })
        stability <- mean(ari.values)
    }

    result <- list(
        membership = consensus$membership,
        n.clusters = consensus$n.clusters,
        modularity = consensus$modularity,
        stability.score = stability,
        modularity.range = range(modularities),
        n.clusters.range = range(sapply(results, function(x) x$n.clusters))
    )

    if (return.all) {
        result$all.runs <- results
    }

    class(result) <- c("consensus_clustering", "list")

    return(result)
}

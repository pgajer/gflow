#' Scan Community Structure Across Resolution Parameters
#'
#' @description
#' Identifies stable regions in the resolution parameter space of graph clustering
#' by comparing partition similarity across a sequence of resolution values. This
#' approach reveals characteristic scales at which community structure persists and
#' helps select appropriate resolution parameters for downstream analysis.
#'
#' @details
#' The resolution parameter in modularity-based clustering controls the characteristic
#' scale at which communities are detected. However, selecting an appropriate value
#' often requires exploring how partition structure varies across the resolution
#' spectrum. This function systematically evaluates clustering solutions at multiple
#' resolution values and quantifies structural stability by measuring similarity
#' between partitions at consecutive resolutions.
#'
#' A stable resolution region, or plateau, occurs when partition structure remains
#' largely unchanged across a range of resolution values. Such plateaus indicate
#' natural scales of community organization in the graph. Conversely, rapid changes
#' in partition similarity suggest transitions between different organizational scales,
#' analogous to phase transitions in statistical physics. Reichardt and Bornholdt (2006)
#' demonstrated that these transitions correspond to changes in the dominant community
#' size through their mapping of modularity optimization to Potts model ground states.
#' The existence of stable plateaus suggests that certain organizational scales are
#' robust to parameter perturbations, while sharp transitions between plateaus indicate
#' hierarchical relationships between community structures at different scales.
#'
#' Three partition similarity metrics are available, each with distinct properties
#' and interpretations. The Adjusted Rand Index (ARI) measures agreement between
#' partitions while correcting for chance, ranging from -1 (less agreement than
#' expected by chance) to 1 (perfect agreement). Values above 0.8 typically indicate
#' high structural similarity, while values below 0.3 suggest substantially different
#' partitions. The ARI has the advantage of adjustment for chance agreement but can
#' be sensitive to partition granularity, particularly when comparing partitions with
#' very different numbers of clusters. Hubert and Arabie (1985) established the
#' theoretical foundations of the Rand Index and its adjusted variant, showing that
#' it properly accounts for the combinatorial structure of partition comparison.
#'
#' Normalized Mutual Information (NMI) quantifies the shared information between
#' partitions normalized by their entropies, ranging from 0 (independent partitions)
#' to 1 (identical partitions). The NMI remains interpretable across partitions of
#' different sizes and provides an information-theoretic measure of structural
#' similarity. Values above 0.9 indicate nearly identical community structure, while
#' values below 0.5 suggest fundamentally different organizational patterns. The
#' normalization makes NMI particularly robust when the number of communities changes
#' substantially across the resolution range. Strehl and Ghosh (2002) demonstrated
#' that NMI provides a principled basis for ensemble clustering methods where multiple
#' partitions must be compared and reconciled.
#'
#' Variation of Information (VI) measures the information lost and gained when
#' transitioning from one partition to another, providing a true metric on the
#' space of clusterings. Unlike ARI and NMI, smaller VI values indicate greater
#' similarity, with 0 representing identical partitions and larger values indicating
#' increasing dissimilarity. The VI satisfies the triangle inequality and provides
#' distance-based interpretation, making it particularly suitable for tracking
#' gradual structural evolution across the resolution spectrum. Typical plateau
#' regions show VI values below 0.5, while structural transitions exhibit VI values
#' above 1.5. The metric property also enables sophisticated hierarchical analyses
#' where one seeks to quantify the total structural change across extended resolution
#' intervals. Meilă (2007) proved that VI is the unique metric on partitions satisfying
#' natural additivity properties under refinement and coarsening operations.
#'
#' The practical workflow involves examining the similarity profile across resolutions
#' to identify plateau regions. A stable plateau is characterized by consistently
#' high similarity (ARI > 0.8, NMI > 0.9, or VI < 0.5) across multiple consecutive
#' resolution values. The choice among metrics depends on the application. The ARI
#' is widely used and intuitive for comparisons where partition sizes remain similar,
#' providing a direct measure of clustering agreement corrected for chance. The NMI
#' is preferred when comparing partitions of very different sizes or when seeking
#' information-theoretic interpretation, as it remains bounded and interpretable
#' regardless of the number of clusters. The VI provides metric properties useful
#' for hierarchical analysis or when measuring cumulative structural change across
#' extended parameter ranges, as distances can be meaningfully summed and compared.
#'
#' When n.consensus.runs is greater than 1, the function performs consensus
#' clustering at each resolution value to account for the stochastic nature of
#' the Louvain algorithm. At each resolution, the clustering is repeated
#' n.consensus.runs times, and a representative partition is selected according
#' to consensus.method. The "best.modularity" method selects the run with highest
#' modularity score, while "median" selects the partition with highest average
#' similarity to all other runs at that resolution. The median method provides
#' more robust results but requires computing pairwise partition similarities.
#'
#' Consensus clustering substantially increases computational cost (factor of
#' n.consensus.runs) but provides more stable and reproducible parameter scanning.
#' For exploratory analysis, n.consensus.runs = 1 (default) is sufficient. For
#' publication-quality results, n.consensus.runs = 50 to 100 is recommended.
#' When consensus is used, additional fields are returned including modularity.sd
#' (standard deviation of modularity across runs at each resolution) and
#' n.clusters.sd (standard deviation of cluster count), which quantify the
#' stability of clustering at each resolution value.
#'
#' @param adj.list A list of length n where \code{adj.list[[i]]} contains the neighbor
#'   indices (1-based) for vertex i.
#' @param weight.list A list of length n where \code{weight.list[[i]]} contains the
#'   edge weights corresponding to neighbors in \code{adj.list[[i]]}. If NULL,
#'   all edges are assigned unit weight. Default is NULL.
#' @param resolution.seq Numeric vector of resolution values to scan. Should be
#'   in increasing order for interpretable similarity profiles. Common ranges
#'   are seq(0.1, 3.0, by = 0.1) for fine-grained scanning or seq(0.5, 2.5, by = 0.25)
#'   for coarser exploration.
#' @param method Character string specifying the partition similarity metric.
#'   Options are "ari" (Adjusted Rand Index), "nmi" (Normalized Mutual Information),
#'   or "vi" (Variation of Information). Default is "ari".
#' @param n.consensus.runs Integer giving the number of independent clustering runs
#'   at each resolution for consensus clustering. When set to 1 (default), performs
#'   a single clustering run at each resolution. Values of 50-100 are recommended
#'   for robust results. Higher values increase computational cost linearly.
#' @param consensus.method Character string specifying how to select the consensus
#'   partition from multiple runs. Options are "max.modularity" (select run with
#'   highest modularity) or "max.avg.similarity" (select partition with highest average similarity
#'   to other runs). Only used when n.consensus.runs > 1. Default is "median".
#' @param return.partitions Logical indicating whether to return the full partition
#'   membership vectors for each resolution. Setting to TRUE enables downstream
#'   analysis but increases memory usage. Default is TRUE.
#'
#' @return A list with components:
#' \describe{
#'   \item{resolution}{Numeric vector of resolution values (same as input resolution.seq).}
#'   \item{n.clusters}{Integer vector giving the number of detected communities at each resolution.}
#'   \item{modularity}{Numeric vector giving the modularity score at each resolution.}
#'   \item{similarity}{Numeric vector of length length(resolution.seq) - 1 giving the
#'     similarity between partition at \code{resolution[i]} and \code{resolution[i+1]}. For VI,
#'     smaller values indicate greater similarity; for ARI and NMI, larger values
#'     indicate greater similarity.}
#'   \item{method}{Character string indicating the similarity metric used.}
#'   \item{partitions}{List of membership vectors (only included if return.partitions = TRUE).
#'     Each element is an integer vector assigning vertices to communities at the
#'     corresponding resolution.}
#' }
#'
#' The returned object has class "resolution_scan" for which print and plot methods
#' may be defined.
#'
#' @references
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. Journal of Classification,
#' 2(1), 193-218.
#'
#' Meilă, M. (2007). Comparing clusterings: an information based distance. Journal of
#' Multivariate Analysis, 98(5), 873-895.
#'
#' Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community detection.
#' Physical Review E, 74(1), 016110.
#'
#' Strehl, A., & Ghosh, J. (2002). Cluster ensembles: a knowledge reuse framework for
#' combining multiple partitions. Journal of Machine Learning Research, 3, 583-617.
#'
#' @examples
#' ## Create a graph with potential multi-scale structure
#' adj.list <- list(
#'   c(2, 3, 4),
#'   c(1, 3, 4),
#'   c(1, 2, 4),
#'   c(1, 2, 3, 5),
#'   c(4, 6, 7),
#'   c(5, 7, 8),
#'   c(5, 6, 8),
#'   c(6, 7)
#' )
#'
#' ## Scan resolutions using ARI
#' result.ari <- scan.resolution(
#'   adj.list,
#'   resolution.seq = seq(0.5, 2.5, by = 0.25),
#'   method = "ari"
#' )
#'
#' ## Identify plateau regions (ARI > 0.8 indicates stability)
#' stable.indices <- which(result.ari$similarity > 0.8)
#' if (length(stable.indices) > 0) {
#'   stable.resolutions <- result.ari$resolution[stable.indices]
#'   cat("Stable resolution range:",
#'       range(stable.resolutions), "\n")
#' }
#'
#' ## Compare metrics
#' result.nmi <- scan.resolution(adj.list,
#'                               resolution.seq = seq(0.5, 2.5, by = 0.25),
#'                               method = "nmi")
#' result.vi <- scan.resolution(adj.list,
#'                              resolution.seq = seq(0.5, 2.5, by = 0.25),
#'                              method = "vi")
#'
#' ## Plot stability profile
#' par(mfrow = c(2, 2))
#'
#' ## ARI profile
#' plot(result.ari$resolution[-1], result.ari$similarity,
#'      type = "b", pch = 19,
#'      xlab = "Resolution", ylab = "ARI",
#'      main = "Adjusted Rand Index")
#' abline(h = 0.8, lty = 2, col = "red")
#'
#' ## NMI profile
#' plot(result.nmi$resolution[-1], result.nmi$similarity,
#'      type = "b", pch = 19,
#'      xlab = "Resolution", ylab = "NMI",
#'      main = "Normalized Mutual Information")
#' abline(h = 0.9, lty = 2, col = "red")
#'
#' ## VI profile (note: lower is more similar)
#' plot(result.vi$resolution[-1], result.vi$similarity,
#'      type = "b", pch = 19,
#'      xlab = "Resolution", ylab = "VI",
#'      main = "Variation of Information")
#' abline(h = 0.5, lty = 2, col = "red")
#'
#' ## Number of clusters vs resolution
#' plot(result.ari$resolution, result.ari$n.clusters,
#'      type = "b", pch = 19,
#'      xlab = "Resolution", ylab = "Number of Clusters",
#'      main = "Community Count")
#'
#' ## Select optimal resolution from largest plateau
#' run.lengths <- rle(result.ari$similarity > 0.8)
#' if (any(run.lengths$values)) {
#'   longest.run.idx <- which.max(run.lengths$lengths * run.lengths$values)
#'   run.start <- sum(run.lengths$lengths[seq_len(longest.run.idx - 1)]) + 1
#'   run.end <- run.start + run.lengths$lengths[longest.run.idx] - 1
#'   optimal.idx <- round((run.start + run.end) / 2)
#'   optimal.resolution <- result.ari$resolution[optimal.idx]
#'   cat(sprintf("\nSuggested resolution from largest plateau: %.2f\n",
#'               optimal.resolution))
#'   cat(sprintf("  Plateau spans resolutions %.2f to %.2f\n",
#'               result.ari$resolution[run.start],
#'               result.ari$resolution[run.end]))
#'   cat(sprintf("  Produces %d clusters with modularity %.3f\n",
#'               result.ari$n.clusters[optimal.idx],
#'               result.ari$modularity[optimal.idx]))
#' }
#'
#' @export
scan.resolution <- function(adj.list,
                           weight.list = NULL,
                           resolution.seq = seq(0.1, 2.5, by = 0.1),
                           method = c("ari", "nmi", "vi"),
                           n.consensus.runs = 1,
                           consensus.method = "max.avg.similarity",
                           return.partitions = TRUE) {

    method <- match.arg(method)
    consensus.method <- match.arg(consensus.method)

    if (!is.numeric(resolution.seq) || length(resolution.seq) < 2) {
        stop("resolution.seq must be a numeric vector with at least 2 elements")
    }

    if (any(resolution.seq <= 0)) {
        stop("All resolution values must be positive")
    }

    if (is.unsorted(resolution.seq)) {
        warning("resolution.seq is not sorted; sorting for interpretable output")
        resolution.seq <- sort(resolution.seq)
    }

    n <- length(adj.list)
    if (n < 2) {
        stop("adj.list must contain at least 2 vertices")
    }

    if (!is.null(weight.list) && length(weight.list) != n) {
        stop("weight.list must have same length as adj.list")
    }

    ## Map method names to igraph::compare method strings
    igraph.method <- switch(method,
                            ari = "adjusted.rand",
                            nmi = "nmi",
                            vi = "vi")

    ## Storage for results
    n.res <- length(resolution.seq)
    n.clusters <- integer(n.res)
    modularity <- numeric(n.res)
    similarity <- numeric(n.res - 1)

    ## NEW: Storage for consensus statistics
    if (n.consensus.runs > 1) {
        modularity.sd <- numeric(n.res)
        n.clusters.sd <- numeric(n.res)
        stability.score <- numeric(n.res)
    }

    partitions <- if (return.partitions) vector("list", n.res) else NULL

    ## Initialize previous membership
    prev.membership <- NULL

    ## Scan across resolutions
    for (i in seq_along(resolution.seq)) {
        gamma <- resolution.seq[i]

        if (n.consensus.runs == 1) {
            ## Original behavior: single run
            result <- cluster.comono.graph(adj.list, weight.list, resolution = gamma)
            curr.membership <- result$membership

        } else {
            ## NEW: Consensus clustering at this resolution
            cat(sprintf("Resolution %.2f: Running %d consensus iterations...\n",
                       gamma, n.consensus.runs))

            ## Run multiple times
            runs <- vector("list", n.consensus.runs)
            for (j in seq_len(n.consensus.runs)) {
                runs[[j]] <- cluster.comono.graph(adj.list, weight.list,
                                                 resolution = gamma)
            }

            ## Extract metrics
            all.modularities <- sapply(runs, function(x) x$modularity)
            all.n.clusters <- sapply(runs, function(x) x$n.clusters)
            all.memberships <- lapply(runs, function(x) x$membership)

            ## Select consensus partition
            if (consensus.method == "max.modularity") {
                best.idx <- which.max(all.modularities)
                result <- runs[[best.idx]]

            } else if (consensus.method == "max.avg.similarity") {
                ## Find partition with highest average similarity to others
                similarity.matrix <- matrix(0, n.consensus.runs, n.consensus.runs)
                for (ii in seq_len(n.consensus.runs - 1)) {
                    for (jj in (ii + 1):n.consensus.runs) {
                        ari <- igraph::compare(all.memberships[[ii]],
                                              all.memberships[[jj]],
                                              method = "adjusted.rand")
                        similarity.matrix[ii, jj] <- ari
                        similarity.matrix[jj, ii] <- ari
                    }
                }
                avg.similarity <- rowMeans(similarity.matrix)
                median.idx <- which.max(avg.similarity)
                result <- runs[[median.idx]]

                ## Store stability score
                stability.score[i] <- avg.similarity[median.idx]
            }

            curr.membership <- result$membership

            ## Store consensus statistics
            modularity.sd[i] <- sd(all.modularities)
            n.clusters.sd[i] <- sd(all.n.clusters)
        }

        ## Store results (same for both paths)
        n.clusters[i] <- result$n.clusters
        modularity[i] <- result$modularity

        if (return.partitions) {
            partitions[[i]] <- curr.membership
        }

        ## Compute similarity with previous partition
        if (i > 1) {
            similarity[i - 1] <- igraph::compare(prev.membership,
                                                 curr.membership,
                                                 method = igraph.method)
        }

        prev.membership <- curr.membership
    }

    ## Construct result object
    result <- list(
        resolution = resolution.seq,
        n.clusters = n.clusters,
        modularity = modularity,
        similarity = similarity,
        method = method,
        n.consensus.runs = n.consensus.runs
    )

    ## Add consensus-specific fields
    if (n.consensus.runs > 1) {
        result$modularity.sd <- modularity.sd
        result$n.clusters.sd <- n.clusters.sd
        result$consensus.method <- consensus.method
        if (consensus.method == "median") {
            result$stability.score <- stability.score
        }
    }

    if (return.partitions) {
        result$partitions <- partitions
    }

    class(result) <- c("resolution_scan", "list")

    return(result)
}

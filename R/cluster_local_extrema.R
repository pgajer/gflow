#' Cluster Local Extrema Based on Basin Overlap
#'
#' @description
#' Identifies clusters of local extrema (either maxima or minima) whose basins
#' of attraction exhibit substantial overlap or containment relationships. The
#' function uses basin labels from \code{\link{summary.basins_of_attraction}}
#' to ensure consistency with other analyses and returns both numeric cluster
#' assignments and a convenient cluster membership list.
#'
#' @details
#' The clustering problem arises naturally in Morse-Smale analysis when multiple
#' local extrema are positioned close together in the function landscape, creating
#' basins that either overlap substantially or exhibit containment relationships
#' where a smaller basin lies entirely within a larger one. Traditional similarity
#' measures like the Jaccard index may fail to detect such containment patterns
#' because they normalize by the union size rather than the smaller set size.
#'
#' We address this by employing the overlap coefficient, also known as the
#' Szymkiewicz-Simpson index, which measures basin similarity through the ratio
#' of intersection size to the minimum basin size. For two basins \eqn{A} and \eqn{B},
#' the overlap coefficient is defined as
#' \deqn{\text{OC}(A,B) = \frac{|A \cap B|}{\min(|A|, |B|)}}
#' This measure equals one when either basin is completely contained in the other,
#' making it ideal for detecting nested or highly overlapping basin structures.
#'
#' The overlap distance, given by \eqn{d(A,B) = 1 - \text{OC}(A,B)}, converts this
#' similarity into a proper metric. We construct a threshold graph where basins
#' are connected by edges whenever their overlap distance falls below the
#' specified threshold. Connected components of this graph define clusters of
#' similar basins. Within each cluster, basins likely represent the same or
#' closely related features in the underlying function landscape.
#'
#' The function internally calls \code{\link{summary.basins_of_attraction}} to
#' obtain consistent basin labels. These labels follow the convention that minima
#' are labeled m1, m2, ... in order of increasing function value (so m1 is the
#' global minimum), while maxima are labeled M1, M2, ... in order of decreasing
#' function value (so M1 is the global maximum). This labeling scheme ensures
#' that clustering results align with other basin analyses.
#'
#' The function operates exclusively on basins of a single extremum type to
#' maintain the interpretability of clusters. Mixing ascending and descending
#' basins would conflate fundamentally different topological features.
#'
#' The overlap graph returned by the function can be visualized using standard
#' graph visualization tools to examine the structure of basin relationships.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}}. This object contains the
#'   complete basin structure for both local minima and local maxima.
#' @param adj.list Adjacency list of the graph.
#' @param edgelen.list A list of numeric vectors containing edge lengths.
#'   Element \code{i} contains the lengths of edges incident to vertex \code{i}.
#'   Must have length equal to the number of vertices in the graph. This parameter
#'   is passed to \code{\link{summary.basins_of_attraction}} to generate
#'   consistent basin labels and compute basin characteristics.
#' @param extrema.type Character string specifying which type of extrema to
#'   cluster. Must be either \code{"max"} for local maxima (descending basins)
#'   or \code{"min"} for local minima (ascending basins). Default is \code{"max"}.
#' @param overlap.threshold Numeric value in the interval \eqn{[0, 1]} specifying the
#'   threshold for the overlap distance below which basins are considered similar
#'   and connected in the clustering graph. Smaller values produce more conservative
#'   clustering with tighter similarity requirements. For example, a threshold of
#'   0.15 requires that basins share at least 85% overlap relative to the smaller
#'   basin. Default is 0.1, corresponding to 90% overlap.
#'
#' @return A list containing the clustering results with the following components:
#'   \describe{
#'     \item{cluster.assignments}{Named integer vector where names are basin
#'       labels from \code{summary.basins_of_attraction} (\code{M1}, \code{M2}, ...
#'       for maxima or \code{m1}, \code{m2}, ... for minima) and values are cluster
#'       identifiers. Basins with the same cluster identifier belong to the same
#'       cluster.}
#'     \item{clusters}{Named list where names are cluster identifiers (as character
#'       strings) and values are character vectors of basin labels belonging to each
#'       cluster. This provides a convenient inverse mapping from clusters to their
#'       member basins.}
#'     \item{overlap.distances}{Symmetric matrix of pairwise overlap distances
#'       between all basins. Rows and columns are labeled with basin identifiers
#'       from \code{summary.basins_of_attraction}. Diagonal entries are zero.}
#'     \item{basin.vertices}{Named list where each element contains the integer
#'       vector of vertex indices belonging to the corresponding basin. Names
#'       match basin labels from \code{summary.basins_of_attraction}.}
#'     \item{overlap.graph}{List with components \code{adj_list} and
#'       \code{weight_list} representing the overlap graph structure. In this graph,
#'       vertices correspond to basins and edges connect basins whose overlap
#'       distance is below the threshold. Edge weights are the overlap distances.
#'       Vertex names in the adjacency list match basin labels.}
#'     \item{basin.summary}{Data frame from \code{summary.basins_of_attraction}
#'       filtered to include only the extrema of the specified type. Contains
#'       detailed information about each basin including label, vertex, value,
#'       basin size, and density metrics.}
#'     \item{n.clusters}{Integer scalar giving the total number of distinct
#'       clusters identified.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins of attraction for a function on a graph
#' basins <- compute.basins.of.attraction(adj.list, adj.list, weight.list, y)
#'
#' # Generate basin summary with consistent labels
#' basin.df <- summary(basins, adj.list, edgelen.list)
#' print(basin.df)
#'
#' # Cluster local maxima using the same labels
#' max.clusters <- cluster.local.extrema(basins,
#'                                       adj.list,
#'                                       edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#'
#' # Report clustering results
#' cat("Found", max.clusters$n.clusters, "clusters among",
#'     nrow(max.clusters$basin.summary), "maxima\n")
#'
#' # Examine cluster assignments
#' print(max.clusters$cluster.assignments)
#'
#' # View clusters and their members
#' print(max.clusters$clusters)
#'
#' # Identify multi-basin clusters
#' multi.basin.clusters <- Filter(function(x) length(x) > 1,
#'                                max.clusters$clusters)
#' if (length(multi.basin.clusters) > 0) {
#'   cat("\nClusters with multiple basins:\n")
#'   for (cluster.id in names(multi.basin.clusters)) {
#'     cat("Cluster", cluster.id, ":",
#'         paste(multi.basin.clusters[[cluster.id]], collapse = ", "), "\n")
#'   }
#' }
#'
#' # Inspect basin characteristics for clustered maxima
#' clustered.labels <- unlist(multi.basin.clusters)
#' clustered.basins <- max.clusters$basin.summary[
#'   max.clusters$basin.summary$label %in% clustered.labels,
#' ]
#' print(clustered.basins[, c("label", "vertex", "value", "basin.size")])
#'
#' # Check overlap distances
#' print(round(max.clusters$overlap.distances, 3))
#'
#' # Visualize the overlap graph
#' library(igraph)
#' n.basins <- length(max.clusters$overlap.graph$adj_list)
#' edge.list <- matrix(nrow = 0, ncol = 2)
#' edge.weights <- numeric(0)
#' for (i in seq_len(n.basins)) {
#'   neighbors <- max.clusters$overlap.graph$adj_list[[i]]
#'   if (length(neighbors) > 0) {
#'     for (j in seq_along(neighbors)) {
#'       neighbor <- neighbors[j]
#'       if (i < neighbor) {  # Avoid duplicate edges
#'         edge.list <- rbind(edge.list, c(i, neighbor))
#'         edge.weights <- c(edge.weights,
#'                          max.clusters$overlap.graph$weight_list[[i]][j])
#'       }
#'     }
#'   }
#' }
#' if (nrow(edge.list) > 0) {
#'   g <- graph_from_edgelist(edge.list, directed = FALSE)
#'   V(g)$name <- names(max.clusters$overlap.graph$adj_list)
#'   E(g)$weight <- edge.weights
#'   plot(g, vertex.label = V(g)$name,
#'        edge.label = round(E(g)$weight, 2),
#'        main = "Basin Overlap Graph")
#' }
#'
#' # Similarly cluster local minima
#' min.clusters <- cluster.local.extrema(basins,
#'                                       adj.list,
#'                                       edgelen.list,
#'                                       extrema.type = "min",
#'                                       overlap.threshold = 0.15)
#' print(min.clusters$clusters)
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing basins of attraction,
#' \code{\link{summary.basins_of_attraction}} for generating basin summaries with
#' consistent labels,
#' \code{\link{create.basin.cx}} for comprehensive basin complex analysis with
#' automatic clustering and simplification
#'
#' @export
cluster.local.extrema <- function(basins.obj,
                                  adj.list,
                                  edgelen.list,
                                  extrema.type = "max",
                                  overlap.threshold = 0.1) {
    ## Input validation
    if (!inherits(basins.obj, "basins_of_attraction")) {
        stop("basins.obj must be of class 'basins_of_attraction'")
    }

    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }

    if (!is.numeric(overlap.threshold) || length(overlap.threshold) != 1) {
        stop("overlap.threshold must be a single numeric value")
    }

    if (overlap.threshold < 0 || overlap.threshold > 1) {
        stop("overlap.threshold must be in the interval [0, 1]")
    }

    ## Generate basin summary to get consistent labels
    basin.summary <- summary(basins.obj, adj.list, edgelen.list)

    ## Filter to selected extrema type
    if (extrema.type == "max") {
        basin.summary <- basin.summary[basin.summary$type == "max", ]
    } else {
        basin.summary <- basin.summary[basin.summary$type == "min", ]
    }

    ## Handle case with no extrema of this type
    if (nrow(basin.summary) == 0) {
        warning(paste("No", extrema.type, "basins found in basins.obj"))
        return(list(
            cluster.assignments = integer(0),
            clusters = list(),
            overlap.distances = matrix(numeric(0), 0, 0),
            basin.vertices = list(),
            overlap.graph = list(adj_list = list(), weight_list = list()),
            basin.summary = basin.summary,
            n.clusters = 0
        ))
    }

    ## Extract the appropriate basins using the summary order
    if (extrema.type == "max") {
        basin.list <- basins.obj$lmax_basins
    } else {
        basin.list <- basins.obj$lmin_basins
    }

    ## Extract vertices for each basin in summary order
    basin.vertices.list <- list()
    for (i in seq_len(nrow(basin.summary))) {
        label <- basin.summary$label[i]
        vertex <- basin.summary$vertex[i]

        ## Find the basin with this vertex
        for (basin in basin.list) {
            if (basin$vertex == vertex) {
                basin.vertices.list[[label]] <- basin$basin_df[, 1]
                break
            }
        }
    }

    ## Compute overlap distance matrix
    overlap.dists <- compute.overlap.distance.matrix(basin.vertices.list)

    ## Create threshold graph
    overlap.graph <- .create.threshold.distance.graph(
        overlap.dists,
        threshold = overlap.threshold
    )

    ## Find connected components (clusters)
    cluster.assignments <- graph.connected.components(overlap.graph$adj_list)
    names(cluster.assignments) <- rownames(overlap.dists)

    ## Create inverse mapping: clusters -> basin labels
    cluster.ids <- sort(unique(cluster.assignments))
    clusters.list <- list()
    for (cluster.id in cluster.ids) {
        basin.labels <- names(cluster.assignments)[cluster.assignments == cluster.id]
        clusters.list[[as.character(cluster.id)]] <- basin.labels
    }

    ## Return clustering information
    list(
        cluster.assignments = cluster.assignments,
        clusters = clusters.list,
        overlap.distances = overlap.dists,
        basin.vertices = basin.vertices.list,
        overlap.graph = overlap.graph,
        basin.summary = basin.summary,
        n.clusters = length(cluster.ids)
    )
}

#' Compute Pairwise Overlap Distance Matrix
#'
#' @param basin_vertices_list A list of integer vectors.
#' @return A symmetric matrix where entry \eqn{(i, j) = 1 - (|A \cap B| / \min(|A|, |B|)),}
#'         measuring the overlap distance between vectors i and j.
compute.overlap.distance.matrix <- function(basin_vertices_list) {
    n_basins <- length(basin_vertices_list)
    if (n_basins == 0) return(matrix(0, 0, 0))

    basin_ids <- names(basin_vertices_list)
    result <- matrix(0, nrow = n_basins, ncol = n_basins)
    rownames(result) <- basin_ids
    colnames(result) <- basin_ids

    ## Calculate overlap distances (1 - Jaccard index)
    for (i in 1:n_basins) {
        vertices_i <- basin_vertices_list[[i]]
        size_i <- length(vertices_i)
        result[i, i] <- 0  ## Self distance = 0

        for (j in (i+1):n_basins) {
            if (j <= n_basins) {  ## Ensure we don't go out of bounds
                vertices_j <- basin_vertices_list[[j]]
                size_j <- length(vertices_j)

                ## Calculate overlap distance
                intersection_size <- length(intersect(vertices_i, vertices_j))

                min_size <- min(size_i, size_j)

                if (min_size > 0) {
                    dist <- 1 - (intersection_size / min_size)
                } else {
                    dist <- 1  ## If both sets are empty, set distance to 1
                }

                ## Store the overlap distance in both directions
                result[i, j] <- dist
                result[j, i] <- dist
            }
        }
    }

    return(result)
}

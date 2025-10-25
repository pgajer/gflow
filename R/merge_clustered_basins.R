#' Merge Clustered Basins of Attraction with Trajectory Preservation
#'
#' @description
#' Merges basins of attraction for local extrema that belong to the same cluster,
#' creating a hierarchical basin structure where each cluster is represented by a
#' single basin with its highest-valued (for maxima) or lowest-valued (for minima)
#' extremum as the representative. The function preserves complete gradient flow
#' information including predecessors and terminal extrema for all absorbed basins,
#' enabling complete trajectory reconstruction in downstream analyses.
#'
#' @details
#' When multiple local extrema have highly overlapping or nested basins of attraction,
#' they often represent the same geometric feature fragmented by noise or sampling
#' artifacts. This function addresses the problem by consolidating such extrema into
#' unified basin structures while maintaining complete trajectory information.
#'
#' We begin by identifying cluster representatives according to a natural ordering
#' principle. For maxima, we select the extremum with the largest function value as
#' the representative, since it corresponds to the highest point in the cluster. For
#' minima, we select the extremum with the smallest function value, corresponding to
#' the lowest point. This choice ensures that the representative captures the most
#' extreme manifestation of the geometric feature.
#'
#' The merged basin is constructed as the set-theoretic union of all basins within
#' the cluster. This union represents the complete region of influence for the
#' geometric feature. We preserve the basin structure of the representative extremum
#' but extend its vertex set to include all vertices from absorbed basins. Critically,
#' we maintain the complete gradient flow structure by storing absorbed basins as
#' hierarchical substructures.
#'
#' For each absorbed extremum, we preserve its complete original basin structure
#' including the predecessor map and terminal extrema. This enables reconstruction of
#' gradient trajectories that terminate at absorbed extrema. Additionally, we compute
#' absorption metadata including the geodesic path connecting the absorbed extremum
#' to the representative and the height of the monotonicity barrier crossed during
#' absorption. This information proves essential for gradient flow cell analysis,
#' where trajectories may pass through or terminate at absorbed extrema.
#'
#' The hierarchical structure supports flexible downstream analyses. For strict
#' monotonicity testing, trajectories can be truncated at absorbed extrema. For
#' permissive analyses, trajectories can be extended through absorbed extrema using
#' the stored geodesic paths, with appropriate flagging of monotonicity violations.
#' The absorption metadata enables assessment of spuriousness, as extrema absorbed
#' over small barriers are more likely to be noise artifacts.
#'
#' The output object retains the structure of a standard basins of attraction result
#' but includes additional hierarchical metadata. This design ensures compatibility
#' with existing analysis pipelines while providing the complete information needed
#' for advanced gradient flow analyses.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}} with gradient basin structures
#'   including predecessors and terminal extrema. This object contains the original
#'   basin structure before merging.
#' @param clustering.result An object returned by \code{\link{cluster.local.extrema}}
#'   containing cluster assignments and basin metadata for the extrema type to be
#'   merged.
#' @param extrema.type Character string specifying which type of extrema were clustered.
#'   Must be either \code{"max"} for local maxima or \code{"min"} for local minima.
#'   This should match the \code{extrema.type} used in the clustering step.
#' @param edgelen.list Optional list of numeric vectors containing edge lengths,
#'   with the same structure as the adjacency list. Used to compute geodesic paths
#'   between absorbed and representative extrema. If NULL, hop distances are used.
#'
#' @return An object of class \code{"basins_of_attraction"} with the same structure
#'   as the input \code{basins.obj} but with clustered basins merged. The object
#'   contains:
#'   \describe{
#'     \item{lmin_basins}{List of basin structures for local minima (merged if
#'       \code{extrema.type = "min"}). Each basin now includes:
#'       \itemize{
#'         \item Standard basin components (vertex, value, hop_idx, basin_df, etc.)
#'         \item \code{predecessors}: Integer vector for trajectory reconstruction
#'         \item \code{terminal_extrema}: Vector of terminal extrema indices
#'         \item \code{absorbed_extrema}: List of absorbed basin structures (if any)
#'       }}
#'     \item{lmax_basins}{List of basin structures for local maxima (merged if
#'       \code{extrema.type = "max"}), with the same structure as lmin_basins.}
#'     \item{n_vertices}{Total number of vertices in the graph.}
#'     \item{y}{Copy of the input function values.}
#'     \item{merge.info}{List containing detailed information about the merging
#'       operation with components:
#'       \describe{
#'         \item{extrema.type}{Type of extrema that were merged (\code{"max"} or
#'           \code{"min"}).}
#'         \item{n.clusters}{Number of clusters identified.}
#'         \item{n.merged}{Number of clusters that involved merging (size > 1).}
#'         \item{cluster.representatives}{Named character vector mapping cluster IDs
#'           to representative basin labels.}
#'         \item{absorbed.extrema}{Data frame with one row per absorbed extremum
#'           containing columns:
#'           \describe{
#'             \item{absorbed.label}{Original label of the absorbed extremum.}
#'             \item{absorbed.vertex}{Vertex index of the absorbed extremum.}
#'             \item{absorbed.value}{Function value at the absorbed extremum.}
#'             \item{representative.label}{Label of the representative basin that
#'               absorbed this extremum.}
#'             \item{representative.vertex}{Vertex index of the representative.}
#'             \item{representative.value}{Function value at the representative.}
#'             \item{y.barrier}{Height of monotonicity barrier crossed.}
#'             \item{cluster.id}{Cluster identifier.}
#'           }}
#'         \item{basin.size.changes}{Data frame showing how basin sizes changed due
#'           to merging, with columns \code{label}, \code{original.size},
#'           \code{merged.size}, \code{size.increase}.}
#'       }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute gradient basins of attraction
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Generate basin summary
#' basin.df <- summary(basins, edgelen.list)
#'
#' # Cluster local maxima
#' max.clusters <- cluster.local.extrema(basins,
#'                                       edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#'
#' # Merge clustered basins with trajectory preservation
#' merged.basins <- merge.clustered.basins(basins,
#'                                         max.clusters,
#'                                         extrema.type = "max",
#'                                         edgelen.list = edgelen.list)
#'
#' # Examine merge information
#' print(merged.basins$merge.info$absorbed.extrema)
#'
#' # Access absorbed basin structures for trajectory reconstruction
#' rep.basin <- merged.basins$lmax_basins[[1]]
#' if (length(rep.basin$absorbed_extrema) > 0) {
#'   absorbed <- rep.basin$absorbed_extrema[[1]]
#'   cat("Absorbed extremum at vertex:", absorbed$vertex, "\n")
#'   cat("Original basin size:", nrow(absorbed$original_basin_df), "\n")
#'   cat("Y barrier crossed:", absorbed$absorption_info$y_barrier, "\n")
#'
#'   # Can still reconstruct trajectories in absorbed basin
#'   v <- absorbed$original_basin_df[1, 1]  # Pick a vertex
#'   traj <- extract.gradient.trajectory.from.structure(absorbed, v)
#' }
#'
#' # For gradient flow cell analysis
#' # Trajectories terminating at absorbed extrema can be identified
#' absorbed.vertices <- sapply(rep.basin$absorbed_extrema, function(a) a$vertex)
#' # Check if any terminal extrema are absorbed
#' absorbed.terminals <- intersect(rep.basin$terminal_extrema, absorbed.vertices)
#' }
#'
#' @seealso
#' \code{\link{cluster.local.extrema}} for clustering basins,
#' \code{\link{compute.basins.of.attraction}} for computing basins,
#' \code{\link{extract.gradient.trajectory}} for trajectory reconstruction
#'
#' @export
merge.clustered.basins <- function(basins.obj,
                                   clustering.result,
                                   extrema.type,
                                   edgelen.list = NULL) {
    ## Input validation
    if (!inherits(basins.obj, "basins_of_attraction")) {
        stop("basins.obj must be of class 'basins_of_attraction'")
    }

    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }

    ## Extract components
    basin.summary <- clustering.result$basin.summary
    clusters <- clustering.result$clusters
    cluster.assignments <- clustering.result$cluster.assignments
    basin.vertices <- clustering.result$basin.vertices

    ## Get original basin list
    if (extrema.type == "max") {
        original.basins <- basins.obj$lmax_basins
    } else {
        original.basins <- basins.obj$lmin_basins
    }

    ## Get y values for barrier computation
    y <- basins.obj$y

    ## Initialize merged basin list and tracking structures
    merged.basins <- list()
    cluster.representatives <- character(0)
    absorbed.extrema.df <- data.frame(
        absorbed.label = character(),
        absorbed.vertex = integer(),
        absorbed.value = numeric(),
        representative.label = character(),
        representative.vertex = integer(),
        representative.value = numeric(),
        y.barrier = numeric(),
        cluster.id = character(),
        stringsAsFactors = FALSE
    )

    basin.size.changes <- data.frame(
        label = character(),
        original.size = integer(),
        merged.size = integer(),
        size.increase = integer(),
        stringsAsFactors = FALSE
    )

    ## Helper function to find basin by vertex
    find.basin.by.vertex <- function(vertex.id, basin.list) {
        for (basin in basin.list) {
            if (basin$vertex == vertex.id) {
                return(basin)
            }
        }
        return(NULL)
    }

    ## Helper function to compute geodesic path between two vertices
    ## Returns list(path = c(from, ..., to), distance = numeric)
    compute.geodesic.path <- function(from.vertex, to.vertex,
                                      adj.list, edge.lengths = NULL) {
        if (is.null(edge.lengths)) {
            ## Use hop distance (all edges weight 1)
            edge.lengths <- lapply(adj.list, function(x) rep(1, length(x)))
        }

        n <- length(adj.list)
        dist <- rep(Inf, n)
        prev <- rep(NA, n)
        dist[from.vertex] <- 0

        ## Simple Dijkstra implementation
        visited <- rep(FALSE, n)

        for (i in seq_len(n)) {
            ## Find unvisited vertex with minimum distance
            min.dist <- Inf
            u <- NA
            for (v in seq_len(n)) {
                if (!visited[v] && dist[v] < min.dist) {
                    min.dist <- dist[v]
                    u <- v
                }
            }

            if (is.na(u) || u == to.vertex) break

            visited[u] <- TRUE

            ## Update distances to neighbors
            if (length(adj.list[[u]]) > 0) {
                for (idx in seq_along(adj.list[[u]])) {
                    v <- adj.list[[u]][idx]
                    edge.len <- edge.lengths[[u]][idx]
                    alt <- dist[u] + edge.len

                    if (alt < dist[v]) {
                        dist[v] <- alt
                        prev[v] <- u
                    }
                }
            }
        }

        ## Reconstruct path
        if (is.infinite(dist[to.vertex])) {
            return(list(path = integer(), distance = Inf))
        }

        path <- integer()
        current <- to.vertex
        while (!is.na(current)) {
            path <- c(current, path)
            current <- prev[current]
        }

        list(path = path, distance = dist[to.vertex])
    }

    ## Helper function to compute y barrier between two extrema
    compute.y.barrier <- function(from.vertex, to.vertex, y.values,
                                  path.vertices = NULL) {
        if (is.null(path.vertices) || length(path.vertices) == 0) {
            ## Direct barrier
            return(abs(y.values[from.vertex] - y.values[to.vertex]))
        }

        ## Maximum deviation from monotonicity along path
        if (extrema.type == "max") {
            ## For maxima: barrier is max deviation above the lower value
            min.y <- min(y.values[from.vertex], y.values[to.vertex])
            max.y.along.path <- max(y.values[path.vertices])
            return(max(0, max.y.along.path - min.y))
        } else {
            ## For minima: barrier is max deviation below the higher value
            max.y <- max(y.values[from.vertex], y.values[to.vertex])
            min.y.along.path <- min(y.values[path.vertices])
            return(max(0, max.y - min.y.along.path))
        }
    }

    ## Process each cluster
    for (cluster.id in names(clusters)) {
        basin.labels <- clusters[[cluster.id]]

        if (length(basin.labels) == 1) {
            ## Single-basin cluster: just copy the basin
            label <- basin.labels[1]
            vertex.idx <- basin.summary$vertex[basin.summary$label == label]

            basin <- find.basin.by.vertex(vertex.idx, original.basins)
            if (!is.null(basin)) {
                ## Initialize absorbed_extrema as empty list
                basin$absorbed_extrema <- list()
                merged.basins[[label]] <- basin
            }

            cluster.representatives[cluster.id] <- label

        } else {
            ## Multi-basin cluster: merge required

            ## Find representative (max value for maxima, min value for minima)
            cluster.basin.info <- basin.summary[basin.summary$label %in% basin.labels, ]

            if (extrema.type == "max") {
                rep.idx <- which.max(cluster.basin.info$value)
            } else {
                rep.idx <- which.min(cluster.basin.info$value)
            }

            rep.label <- cluster.basin.info$label[rep.idx]
            rep.vertex <- cluster.basin.info$vertex[rep.idx]
            rep.value <- cluster.basin.info$value[rep.idx]

            cluster.representatives[cluster.id] <- rep.label

            ## Get representative's original basin
            rep.basin <- find.basin.by.vertex(rep.vertex, original.basins)
            if (is.null(rep.basin)) {
                warning("Representative basin not found for cluster ", cluster.id)
                next
            }

            ## Compute union of all basin vertices in cluster
            union.vertices <- unique(unlist(basin.vertices[basin.labels]))
            union.vertices <- sort(union.vertices)

            ## Create merged basin structure (copy representative)
            merged.basin <- rep.basin

            ## Update basin_df to include all union vertices
            original.vertices <- rep.basin$basin_df[, 1]
            max.hop <- max(rep.basin$basin_df[, 2])

            new.vertices <- setdiff(union.vertices, original.vertices)
            if (length(new.vertices) > 0) {
                ## Add new vertices with hop distance = max.hop + 1
                new.rows <- cbind(new.vertices, rep(max.hop + 1, length(new.vertices)))
                merged.basin$basin_df <- rbind(rep.basin$basin_df, new.rows)

                ## Sort by vertex index
                merged.basin$basin_df <- merged.basin$basin_df[
                    order(merged.basin$basin_df[, 1]),
                ]
            }

            ## Update terminal_extrema: remove absorbed extrema from terminal list
            absorbed.labels <- setdiff(basin.labels, rep.label)
            absorbed.vertices.in.cluster <- sapply(absorbed.labels, function(lbl) {
                cluster.basin.info$vertex[cluster.basin.info$label == lbl]
            })

            ## Remove absorbed extrema from terminal_extrema list
            merged.basin$terminal_extrema <- setdiff(
                merged.basin$terminal_extrema,
                absorbed.vertices.in.cluster
            )

            ## Initialize absorbed_extrema list
            merged.basin$absorbed_extrema <- list()

            ## Process each absorbed extremum
            for (absorbed.label in absorbed.labels) {
                absorbed.info <- cluster.basin.info[
                    cluster.basin.info$label == absorbed.label,
                ]
                absorbed.vertex <- absorbed.info$vertex
                absorbed.value <- absorbed.info$value

                ## Get original basin structure for absorbed extremum
                absorbed.basin <- find.basin.by.vertex(absorbed.vertex,
                                                      original.basins)

                if (is.null(absorbed.basin)) {
                    warning("Absorbed basin not found for ", absorbed.label)
                    next
                }

                ## Compute geodesic path from absorbed to representative
                geodesic.result <- compute.geodesic.path(
                    absorbed.vertex,
                    rep.vertex,
                    basins.obj$adj.list,
                    edgelen.list
                )

                ## Compute y barrier
                y.barrier <- compute.y.barrier(
                    absorbed.vertex,
                    rep.vertex,
                    y,
                    geodesic.result$path
                )

                ## Create absorbed extremum structure with complete basin info
                absorbed.structure <- list(
                    vertex = absorbed.vertex,
                    value = absorbed.value,
                    label = absorbed.label,

                    ## Original basin structure (for trajectory reconstruction)
                    original_basin_df = absorbed.basin$basin_df,
                    original_basin_bd_df = absorbed.basin$basin_bd_df,
                    original_predecessors = absorbed.basin$predecessors,
                    original_terminal_extrema = absorbed.basin$terminal_extrema,
                    original_hop_idx = absorbed.basin$hop_idx,

                    ## Absorption metadata
                    absorption_info = list(
                        representative_vertex = rep.vertex,
                        representative_value = rep.value,
                        representative_label = rep.label,
                        geodesic_path = geodesic.result$path,
                        geodesic_distance = geodesic.result$distance,
                        y_barrier = y.barrier
                    )
                )

                ## Add to absorbed_extrema list
                merged.basin$absorbed_extrema[[absorbed.label]] <- absorbed.structure

                ## Record in tracking dataframe
                absorbed.extrema.df <- rbind(
                    absorbed.extrema.df,
                    data.frame(
                        absorbed.label = absorbed.label,
                        absorbed.vertex = absorbed.vertex,
                        absorbed.value = absorbed.value,
                        representative.label = rep.label,
                        representative.vertex = rep.vertex,
                        representative.value = rep.value,
                        y.barrier = y.barrier,
                        cluster.id = cluster.id,
                        stringsAsFactors = FALSE
                    )
                )
            }

            ## Store merged basin with representative's label
            merged.basins[[rep.label]] <- merged.basin

            ## Record basin size change
            original.size <- nrow(rep.basin$basin_df)
            merged.size <- nrow(merged.basin$basin_df)

            basin.size.changes <- rbind(
                basin.size.changes,
                data.frame(
                    label = rep.label,
                    original.size = original.size,
                    merged.size = merged.size,
                    size.increase = merged.size - original.size,
                    stringsAsFactors = FALSE
                )
            )
        }
    }

    ## Create result object
    result <- basins.obj

    if (extrema.type == "max") {
        result$lmax_basins <- merged.basins
    } else {
        result$lmin_basins <- merged.basins
    }

    ## Add merge information
    result$merge.info <- list(
        extrema.type = extrema.type,
        n.clusters = length(clusters),
        n.merged = sum(sapply(clusters, length) > 1),
        cluster.representatives = cluster.representatives,
        absorbed.extrema = absorbed.extrema.df,
        basin.size.changes = basin.size.changes
    )

    return(result)
}


#' Extract Gradient Trajectory from Absorbed Basin Structure
#'
#' @description
#' Helper function to reconstruct gradient trajectories from absorbed basin
#' structures stored within merged basins.
#'
#' @param absorbed.structure Absorbed extremum structure from
#'   merged.basin$absorbed_extrema
#' @param vertex.id Vertex identifier (1-based)
#'
#' @return Integer vector of trajectory or NULL if vertex not in basin
#'
#' @keywords internal
extract.gradient.trajectory.from.structure <- function(absorbed.structure,
                                                       vertex.id) {
    ## Check if vertex is in original basin
    if (!(vertex.id %in% absorbed.structure$original_basin_df[, 1])) {
        return(NULL)
    }

    trajectory <- integer()
    current <- vertex.id
    predecessors <- absorbed.structure$original_predecessors

    ## Backtrack through predecessors
    while (current != 0) {
        trajectory <- c(current, trajectory)
        current <- predecessors[current]
    }

    return(trajectory)
}


#' Get All Absorbed Extrema for a Merged Basin
#'
#' @description
#' Extracts information about all extrema absorbed by a given basin during merging.
#'
#' @param merged.basin A basin structure from a merged basins object
#'
#' @return Data frame with columns: vertex, value, label, y.barrier,
#'   geodesic.distance
#'
#' @export
get.absorbed.extrema.info <- function(merged.basin) {
    if (length(merged.basin$absorbed_extrema) == 0) {
        return(data.frame(
            vertex = integer(),
            value = numeric(),
            label = character(),
            y.barrier = numeric(),
            geodesic.distance = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    info.list <- lapply(merged.basin$absorbed_extrema, function(abs) {
        data.frame(
            vertex = abs$vertex,
            value = abs$value,
            label = abs$label,
            y.barrier = abs$absorption_info$y_barrier,
            geodesic.distance = abs$absorption_info$geodesic_distance,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, info.list)
}

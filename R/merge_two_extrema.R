#' Merge Two Specified Local Extrema with Trajectory Preservation
#'
#' @description
#' Merges two specific local extrema of the same type by absorbing one basin into
#' another, creating a unified hierarchical basin structure. This function enables
#' manual expert-based simplification while preserving complete gradient flow
#' information including predecessors and terminal extrema for downstream analysis.
#'
#' @details
#' The merging problem arises when visual examination of the response surface or
#' consideration of the physical process suggests that two nearby local extrema
#' should be treated as a single feature. Perhaps the separation between them results
#' from noise in the data, or perhaps the analyst judges based on domain knowledge
#' that the distinction lacks practical significance. In such cases, we wish to
#' consolidate the two basins while preserving complete information about the
#' operation for downstream analysis.
#'
#' We begin by identifying the two extrema to merge through their basin labels,
#' which provide unambiguous references to specific basins. The basin lists
#' (\code{lmax_basins} and \code{lmin_basins}) are named with these labels,
#' enabling direct access without requiring auxiliary lookup structures. The user
#' designates one extremum as the winner, which will represent the merged feature,
#' while the other becomes the loser whose basin will be absorbed. This choice may
#' reflect which extremum has the more extreme function value, which lies in a more
#' central location, or which better represents the feature according to domain expertise.
#'
#' The merged basin is constructed as the set-theoretic union of the two original
#' basins. All vertices from the loser's basin are incorporated into the winner's
#' basin structure, extending the region of influence. The winner's extremum
#' vertex and function value remain as the representative for the merged basin,
#' while the basin vertex set expands to include all vertices from both basins.
#'
#' Critically, this implementation preserves the complete gradient flow structure
#' of the absorbed basin. The loser's original basin structure, including its
#' predecessor map and terminal extrema, is stored hierarchically within the
#' winner's basin. This enables reconstruction of gradient trajectories that
#' terminate at the absorbed extremum, which proves essential for gradient flow
#' cell analysis and monotonicity testing. Additionally, we compute absorption
#' metadata including the geodesic path connecting the loser to the winner and
#' the height of the monotonicity barrier crossed. This information quantifies
#' the "spuriousness" of the absorbed extremum and supports flexible downstream
#' analysis strategies.
#'
#' The winner's terminal extrema are updated to remove the loser if it was
#' previously a terminal extremum of the winner's basin. This maintains consistency
#' in the gradient flow cell structure, as the loser no longer exists as an
#' independent extremum but is preserved as an absorbed substructure.
#'
#' We preserve the hop distance structure from the winner's basin for its original
#' vertices. For vertices absorbed from the loser's basin, we assign a hop distance
#' one greater than the maximum hop distance in the winner's original basin. This
#' convention maintains the breadth-first search structure while clearly indicating
#' which vertices came from the absorbed basin.
#'
#' The function returns a modified basins of attraction object maintaining the
#' same structure as the input but with hierarchical absorption tracking. Complete
#' metadata about the operation appears in a merge.info component, documenting
#' which extremum was absorbed, both extrema's vertices and function values, the
#' y-barrier crossed, and how the basin size changed. This hierarchical structure
#' supports flexible analysis strategies ranging from strict monotonicity testing
#' (respecting absorbed extrema as endpoints) to permissive approaches (extending
#' trajectories through absorbed extrema).
#'
#' The user may apply this function sequentially to merge multiple pairs of extrema,
#' building up a progressively simplified basin structure through a series of
#' pairwise operations. Each merge operation preserves information about previous
#' merges in both the basin-level absorbed_extrema lists and the global merge.info
#' component, allowing complete reconstruction of the simplification history.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}} with gradient basin structures
#'   including predecessors and terminal extrema. The basin lists (\code{lmax_basins}
#'   and \code{lmin_basins}) must be named with extremum labels.
#' @param winner.label Character string or numeric value specifying the label of
#'   the winning basin that will absorb the other basin. This label must appear
#'   as a name in the appropriate basin list (e.g., in \code{names(basins.obj$lmin_basins)}
#'   for minima).
#' @param loser.label Character string or numeric value specifying the label of
#'   the losing basin to be absorbed. This basin will be removed from the basin list,
#'   its vertices incorporated into the winner's basin, and its complete structure
#'   preserved in the winner's absorbed_extrema list.
#'
#' @return An object of class \code{"basins_of_attraction"} with the same structure
#'   as the input \code{basins.obj} but with the specified basins merged. The object
#'   contains:
#'   \describe{
#'     \item{lmin_basins}{Named list of basin structures for local minima (modified
#'       if \code{extrema.type = "min"}). The winner basin now includes an
#'       \code{absorbed_extrema} list with the loser's complete structure. The loser's
#'       entry is removed from the list.}
#'     \item{lmax_basins}{Named list of basin structures for local maxima (modified
#'       if \code{extrema.type = "max"}), with the same hierarchical structure.}
#'     \item{n_vertices}{Total number of vertices in the graph.}
#'     \item{y}{Copy of the input function values.}
#'     \item{merge.info}{List containing information about merge operations. If
#'       previous merges were performed, this list accumulates information from
#'       all operations. For the current merge, the following components are added
#'       or updated:
#'       \describe{
#'         \item{operations}{Data frame with one row per merge operation.}
#'         \item{absorbed.extrema}{Data frame tracking all absorbed extrema across
#'           all merge operations.}
#'       }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute refined basins with named lists
#' result <- compute.refined.basins(adj.list, edgelen.list, fitted.values)
#' basins <- result$basins
#'
#' # Examine the summary to identify basins to merge
#' print(result$summary[result$summary$type == "min", ])
#'
#' # Merge two minima with trajectory preservation
#' # Labels are now directly accessible as names
#' merged.basins <- merge.two.extrema(basins,
#'                                    winner.label = "42",
#'                                    loser.label = "37")
#'
#' # Examine merge information
#' print(merged.basins$merge.info$operations)
#'
#' # Access absorbed basin structure
#' winner.basin <- merged.basins$lmin_basins[["42"]]
#' if (length(winner.basin$absorbed_extrema) > 0) {
#'   absorbed <- winner.basin$absorbed_extrema[["37"]]
#'   cat("Absorbed extremum:", absorbed$vertex, "\n")
#'   cat("Y barrier:", absorbed$absorption_info$y_barrier, "\n")
#' }
#' }
#'
#' @export
merge.two.extrema <- function(basins.obj,
                              winner.label,
                              loser.label,
                              verbose = FALSE) {
    
    ## Validate inputs
    if (!is.list(basins.obj)) {
        stop("basins.obj must be a list")
    }

    required.components <- c("lmax_basins", "lmin_basins", "y")
    if (!all(required.components %in% names(basins.obj$basins))) {
        stop("basins.obj$basins must contain lmax_basins, lmin_basins, and y components")
    }

    required.components <- c("lmax_basins", "lmin_basins", "y", "adj.list", "edge.length.list")
    missing.components <- setdiff(required.components, names(basins.obj$basins))
    if (length(missing.components) > 0) {
        stop(sprintf("basins.obj must contain: %s. Missing: %s",
                     paste(required.components, collapse = ", "),
                     paste(missing.components, collapse = ", ")))
    }

    adj.list <- basins.obj$basins$adj.list
    edgelen.list <- basins.obj$basins$edge.length.list

    ## Convert labels to character for consistent indexing
    winner.label.char <- as.character(winner.label)
    loser.label.char <- as.character(loser.label)

    ## Automatically determine extrema type by checking which basin list contains both labels
    if (verbose) cat("  Detecting extrema type from labels...\n")

    winner.loser.labels <- c(winner.label.char, loser.label.char)

    ## Check if both labels are in lmax_basins
    ## NOTE: Use intersect() not intersection()
    common.in.max <- intersect(winner.loser.labels, names(basins.obj$basins$lmax_basins))
    extrema.type <- NA
    if (length(common.in.max) == 2) {
        ## Both labels found in maxima
        basin.list <- basins.obj$basins$lmax_basins
        list.name <- "lmax_basins"
        extrema.type <- "max"
        if (verbose) {
            cat(sprintf("    Detected: maximum (both labels in lmax_basins)\n"))
        }
    } else {
        ## Check if both labels are in lmin_basins
        common.in.min <- intersect(winner.loser.labels, names(basins.obj$basins$lmin_basins))
        if (length(common.in.min) == 2) {
            ## Both labels found in minima
            basin.list <- basins.obj$basins$lmin_basins
            list.name <- "lmin_basins"
            extrema.type <- "min"
            if (verbose) {
                cat(sprintf("    Detected: minimum (both labels in lmin_basins)\n"))
            }
        } else {
            ## ERROR CASE: Labels not both in same list
            ## Determine what went wrong and provide helpful error message

            ## Check where each label was found
            winner.in.max <- winner.label.char %in% names(basins.obj$basins$lmax_basins)
            winner.in.min <- winner.label.char %in% names(basins.obj$basins$lmin_basins)
            loser.in.max <- loser.label.char %in% names(basins.obj$basins$lmax_basins)
            loser.in.min <- loser.label.char %in% names(basins.obj$basins$lmin_basins)

            ## Build informative error message based on what was found
            error.msg <- sprintf("Cannot merge extrema '%s' and '%s': ",
                                 winner.label.char, loser.label.char)

            if (winner.in.max && loser.in.min) {
                ## Winner is max, loser is min
                error.msg <- paste0(error.msg,
                                    sprintf("'%s' is a maximum but '%s' is a minimum. ",
                                            winner.label.char, loser.label.char),
                                    "Cannot merge extrema of different types.")

            } else if (winner.in.min && loser.in.max) {
                ## Winner is min, loser is max
                error.msg <- paste0(error.msg,
                                    sprintf("'%s' is a minimum but '%s' is a maximum. ",
                                            winner.label.char, loser.label.char),
                                    "Cannot merge extrema of different types.")

            } else if (!winner.in.max && !winner.in.min && !loser.in.max && !loser.in.min) {
                ## Neither label found anywhere
                error.msg <- paste0(error.msg,
                                    "Neither label found in basin lists.\n",
                                    sprintf("  Available maxima: %s\n",
                                            paste(names(basins.obj$basins$lmax_basins), collapse = ", ")),
                                    sprintf("  Available minima: %s",
                                            paste(names(basins.obj$basins$lmin_basins), collapse = ", ")))

            } else if (!winner.in.max && !winner.in.min) {
                ## Only winner not found
                error.msg <- paste0(error.msg,
                                    sprintf("Winner label '%s' not found in basin lists.\n",
                                            winner.label.char),
                                    sprintf("  Available maxima: %s\n",
                                            paste(names(basins.obj$basins$lmax_basins), collapse = ", ")),
                                    sprintf("  Available minima: %s",
                                            paste(names(basins.obj$basins$lmin_basins), collapse = ", ")))

            } else if (!loser.in.max && !loser.in.min) {
                ## Only loser not found
                error.msg <- paste0(error.msg,
                                    sprintf("Loser label '%s' not found in basin lists.\n",
                                            loser.label.char),
                                    sprintf("  Available maxima: %s\n",
                                            paste(names(basins.obj$basins$lmax_basins), collapse = ", ")),
                                    sprintf("  Available minima: %s",
                                            paste(names(basins.obj$basins$lmin_basins), collapse = ", ")))

            } else {
                ## Unexpected case (e.g., label in both lists somehow)
                error.msg <- paste0(error.msg,
                                    "Labels found in unexpected configuration.\n",
                                    sprintf("  Winner '%s': in_maxima=%s, in_minima=%s\n",
                                            winner.label.char, winner.in.max, winner.in.min),
                                    sprintf("  Loser '%s': in_maxima=%s, in_minima=%s",
                                            loser.label.char, loser.in.max, loser.in.min))
            }

            stop(error.msg)
        }
    }

    ## Validate that basin lists are named
    if (is.null(names(basin.list))) {
        stop(sprintf("%s must be a named list with extremum labels as names", list.name))
    }

    ## Extract winner and loser basins using label names
    winner.basin <- basin.list[[winner.label.char]]
    loser.basin <- basin.list[[loser.label.char]]
    
    ## Extract basin information
    winner.vertex <- winner.basin$vertex
    loser.vertex <- loser.basin$vertex
    
    y <- basins.obj$basins$y
    winner.value <- y[winner.vertex]
    loser.value <- y[loser.vertex]
    
    winner.vertices <- winner.basin$basin_df[, 1]
    loser.vertices <- loser.basin$basin_df[, 1]
    
    winner.original.size <- length(winner.vertices)
    loser.size <- length(loser.vertices)
    
    ## Helper function to compute geodesic path using Dijkstra's algorithm
    compute.geodesic.path <- function(from.vertex, to.vertex, adj.list, edge.lengths) {
        n.vertices <- length(adj.list)
        dist <- rep(Inf, n.vertices)
        prev <- rep(NA, n.vertices)
        visited <- rep(FALSE, n.vertices)
        
        dist[from.vertex] <- 0
        
        while (TRUE) {
            unvisited.dist <- dist
            unvisited.dist[visited] <- Inf
            
            if (all(is.infinite(unvisited.dist))) {
                break
            }
            
            current <- which.min(unvisited.dist)
            
            if (current == to.vertex) {
                break
            }
            
            visited[current] <- TRUE
            
            neighbors <- adj.list[[current]]
            if (is.null(edge.lengths)) {
                edge.lens <- rep(1, length(neighbors))
            } else {
                edge.lens <- edge.lengths[[current]]
            }
            
            for (i in seq_along(neighbors)) {
                neighbor <- neighbors[i]
                if (!visited[neighbor]) {
                    alt <- dist[current] + edge.lens[i]
                    if (alt < dist[neighbor]) {
                        dist[neighbor] <- alt
                        prev[neighbor] <- current
                    }
                }
            }
        }
        
        ## Reconstruct path
        if (is.infinite(dist[to.vertex])) {
            return(list(path = integer(0), distance = Inf))
        }
        
        path <- c()
        current <- to.vertex
        while (!is.na(prev[current])) {
            path <- c(prev[current], path)
            current <- prev[current]
        }
        path <- c(path, to.vertex)
        
        list(path = path, distance = dist[to.vertex])
    }
    
    ## Helper function to compute y barrier
    compute.y.barrier <- function(from.vertex, to.vertex, y.values,
                                  path.vertices = NULL) {
        if (is.null(path.vertices) || length(path.vertices) == 0) {
            return(abs(y.values[from.vertex] - y.values[to.vertex]))
        }
        
        if (extrema.type == "max") {
            min.y <- min(y.values[from.vertex], y.values[to.vertex])
            max.y.along.path <- max(y.values[path.vertices])
            return(max(0, max.y.along.path - min.y))
        } else {
            max.y <- max(y.values[from.vertex], y.values[to.vertex])
            min.y.along.path <- min(y.values[path.vertices])
            return(max(0, max.y - min.y.along.path))
        }
    }
    
    ## Compute geodesic path and y barrier
    geodesic.result <- compute.geodesic.path(
        loser.vertex,
        winner.vertex,
        adj.list,
        edgelen.list
    )
    
    y.barrier <- compute.y.barrier(
        loser.vertex,
        winner.vertex,
        y,
        geodesic.result$path
    )
    
    ## Create merged basin
    merged.basin <- winner.basin
    
    ## Initialize absorbed_extrema if it doesn't exist
    if (is.null(merged.basin$absorbed_extrema)) {
        merged.basin$absorbed_extrema <- list()
    }
    
    ## Determine hop distance for absorbed vertices
    max.winner.hop <- max(winner.basin$basin_df[, 2])
    absorbed.hop <- max.winner.hop + 1
    
    ## Add loser vertices that are not in winner
    new.vertices <- setdiff(loser.vertices, winner.vertices)
    if (length(new.vertices) > 0) {
        new.rows <- cbind(new.vertices, rep(absorbed.hop, length(new.vertices)))
        merged.basin$basin_df <- rbind(winner.basin$basin_df, new.rows)
        
        ## Sort by vertex index
        merged.basin$basin_df <- merged.basin$basin_df[
            order(merged.basin$basin_df[, 1]),
        ]
    }
    
    ## Update hop_idx
    merged.basin$hop_idx <- max(merged.basin$basin_df[, 2])
    
    ## Update terminal_extrema: remove loser if present
    if (!is.null(merged.basin$terminal_extrema)) {
        merged.basin$terminal_extrema <- setdiff(
            merged.basin$terminal_extrema,
            loser.vertex
        )
    }
    
    ## Create absorbed extremum structure with complete basin info
    absorbed.structure <- list(
        vertex = loser.vertex,
        value = loser.value,
        label = loser.label.char,
        
        ## Original basin structure (for trajectory reconstruction)
        original_basin_df = loser.basin$basin_df,
        original_basin_bd_df = loser.basin$basin_bd_df,
        original_predecessors = loser.basin$predecessors,
        original_terminal_extrema = loser.basin$terminal_extrema,
        original_hop_idx = loser.basin$hop_idx,
        
        ## Absorption metadata
        absorption_info = list(
            representative_vertex = winner.vertex,
            representative_value = winner.value,
            representative_label = winner.label.char,
            geodesic_path = geodesic.result$path,
            geodesic_distance = geodesic.result$distance,
            y_barrier = y.barrier
        )
    )
    
    ## Add to absorbed_extrema list
    merged.basin$absorbed_extrema[[loser.label.char]] <- absorbed.structure
    
    ## Merged size
    winner.merged.size <- nrow(merged.basin$basin_df)
    
    ## Create merge operation record
    operation.record <- data.frame(
        extrema.type = extrema.type,
        winner.label = winner.label.char,
        winner.vertex = winner.vertex,
        winner.value = winner.value,
        loser.label = loser.label.char,
        loser.vertex = loser.vertex,
        loser.value = loser.value,
        y.barrier = y.barrier,
        winner.original.size = winner.original.size,
        winner.merged.size = winner.merged.size,
        loser.size = loser.size,
        stringsAsFactors = FALSE
    )
    
    ## Update basin list using named access
    basin.list[[winner.label.char]] <- merged.basin
    basin.list[[loser.label.char]] <- NULL
    
    ## Update result object
    result <- basins.obj

    if (extrema.type == "max") {
        result$basins$lmax_basins <- basin.list
    } else {
        result$basins$lmin_basins <- basin.list
    }

    ## Update summary
    result$summary <- summary(result$basins, adj.list, edge.length.list, hop.k = basins.obj$basins$hop.k)

    ## Update or create merge.info
    if (is.null(basins.obj$basins$merge.info)) {
        ## First merge operation
        operation.record$operation.id <- 1
        
        result$basins$merge.info <- list(
            operations = operation.record,
            absorbed.extrema = data.frame(
                absorbed.label = loser.label.char,
                absorbed.vertex = loser.vertex,
                absorbed.value = loser.value,
                representative.label = winner.label.char,
                representative.vertex = winner.vertex,
                representative.value = winner.value,
                y.barrier = y.barrier,
                stringsAsFactors = FALSE
            )
        )
    } else {
        ## Append to existing merge information
        n.operations <- nrow(basins.obj$basins$merge.info$operations)
        operation.record$operation.id <- n.operations + 1
        
        result$basins$merge.info$operations <- rbind(
            basins.obj$basins$merge.info$operations,
            operation.record
        )
        
        result$basins$merge.info$absorbed.extrema <- rbind(
            basins.obj$basins$merge.info$absorbed.extrema,
            data.frame(
                absorbed.label = loser.label.char,
                absorbed.vertex = loser.vertex,
                absorbed.value = loser.value,
                representative.label = winner.label.char,
                representative.vertex = winner.vertex,
                representative.value = winner.value,
                y.barrier = y.barrier,
                stringsAsFactors = FALSE
            )
        )
    }

    return(result)
}

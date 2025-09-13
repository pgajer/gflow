#' Graph Gradient Flow Analysis
#'
#' @description
#' Analyzes the structure of a function defined on a weighted graph by computing
#' gradient trajectories, basins, and Morse-Smale pro-cells. This function identifies
#' local extrema (minima and maxima) and constructs trajectories that follow the gradient
#' of the function across the graph.
#'
#' The gradient flow computation relies on a scale parameter that defines the local
#' neighborhood size for each vertex, within which the gradient direction is determined.
#' Trajectories are constructed by following these gradient directions to connect
#' local minima to local maxima.
#'
#' @param adj.list List of integer vectors where \code{adj.list[[i]]} contains the indices
#'   of vertices adjacent to vertex i (using 1-based indexing).
#' @param weight.list List of numeric vectors where \code{weight.list[[i]][j]} is the weight
#'   of the edge from vertex i to \code{adj.list[[i]][j]}.
#' @param y Numeric vector of function values at each graph vertex.
#' @param scale Numeric vector of positive scale values for each vertex, controlling the local
#'   neighborhood size for gradient computation and extrema detection. If a single value,
#'   it is replicated for all vertices.
#' @param quantile.scale.thld Numeric scalar between 0 and 1. Scale values are truncated to
#'   the interval \code{(0, quantile(scale, prob = quantile.scale.thld))}. Default is 0.025.
#' @param with.trajectories Logical; if \code{TRUE}, trajectory information is included in
#'   the output. Default is \code{FALSE}.
#'
#' @return An object of class \code{"ggflow"} containing:
#' \describe{
#'   \item{local_extrema}{Data frame with columns:
#'     \describe{
#'       \item{vertex_index}{Integer; vertex index (1-based)}
#'       \item{is_maximum}{Integer; 1 for maxima, 0 for minima}
#'       \item{label}{Character; automatically generated label (e.g., "M1", "m2")}
#'       \item{fn_value}{Numeric; function value at the extremum}
#'     }
#'   }
#'   \item{trajectories}{List of trajectories (NULL if \code{with.trajectories = FALSE}).
#'     Each trajectory contains:
#'     \describe{
#'       \item{vertices}{Integer vector of vertex indices forming the trajectory}
#'       \item{type}{Character string: "LMIN_LMAX", "LMIN_ONLY", "LMAX_ONLY", or "UNKNOWN"}
#'       \item{label}{Character; trajectory label based on endpoints}
#'     }
#'   }
#'   \item{basins}{List with two components:
#'     \describe{
#'       \item{ascending}{List of basins around local minima}
#'       \item{descending}{List of basins around local maxima}
#'     }
#'     Each basin contains:
#'     \describe{
#'       \item{local_min/local_max}{Integer; index of extremum vertex}
#'       \item{vertices}{Integer vector; vertices in the basin}
#'       \item{label}{Character; basin label}
#'     }
#'   }
#'   \item{cells}{List of gradient flow pro-cells. Each cell contains:
#'     \describe{
#'       \item{local_min}{Integer; index of local minimum vertex}
#'       \item{local_max}{Integer; index of local maximum vertex}
#'       \item{vertices}{Integer vector; vertices in the pro-cell}
#'       \item{label}{Character; cell label}
#'     }
#'   }
#' }
#'
#' @details
#' The gradient flow computation uses the scale parameter to determine the local neighborhood
#' for each vertex. For a vertex \code{v}, the algorithm finds all shortest paths of length
#' \code{<= scale[v]} starting at \code{v}, and selects the path with the steepest rate of
#' change in function value. This approach ensures that the gradient direction is determined
#' at an appropriate scale for each part of the graph.
#'
#' The returned pro-cells partition the graph into regions, where each region contains
#' trajectories flowing from a specific local minimum to a specific local maximum.
#' This decomposition provides insight into the topological structure of the function
#' defined on the graph.
#'
#' @examples
#' \dontrun{
#' # Create a simple grid graph
#' n <- 10
#' adj.list <- list()
#' weight.list <- list()
#'
#' # Build grid adjacency (simplified example)
#' for(i in 1:(n*n)) {
#'   neighbors <- integer()
#'   weights <- numeric()
#'
#'   # Add horizontal neighbors
#'   if(i %% n != 1) {
#'     neighbors <- c(neighbors, i-1)
#'     weights <- c(weights, 1.0)
#'   }
#'   if(i %% n != 0) {
#'     neighbors <- c(neighbors, i+1)
#'     weights <- c(weights, 1.0)
#'   }
#'
#'   adj.list[[i]] <- neighbors
#'   weight.list[[i]] <- weights
#' }
#'
#' # Define a function on the graph (e.g., a Gaussian peak)
#' coords <- expand.grid(x = 1:n/n, y = 1:n/n)
#' z <- exp(-10 * ((coords$x - 0.5)^2 + (coords$y - 0.5)^2))
#'
#' # Compute gradient flow
#' flow <- construct.graph.gradient.flow(
#'   adj.list,
#'   weight.list,
#'   z,
#'   scale = 2.0,
#'   with.trajectories = TRUE
#' )
#'
#' # Examine results
#' print(flow$local_extrema)
#' summary(flow)
#' }
#'
#' @seealso \code{\link{summary.ggflow}}
#'   \code{\link{basins.merge}}, \code{\link{replace.basin.label}}
#'
#' @export
#' @importFrom stats quantile setNames
construct.graph.gradient.flow <- function(adj.list,
                                          weight.list,
                                          y,
                                          scale,
                                          quantile.scale.thld = 0.025,
                                          with.trajectories = FALSE) {

    # Input validation with informative error messages
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!all(vapply(adj.list, function(x) is.numeric(x) || is.integer(x),
                    logical(1L)))) {
        stop("All elements of 'adj.list' must be numeric or integer vectors",
             call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    if (!all(vapply(weight.list, is.numeric, logical(1L)))) {
        stop("All elements of 'weight.list' must be numeric vectors",
             call. = FALSE)
    }

    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length",
             call. = FALSE)
    }

    # Check matching lengths of corresponding elements
    adj.lengths <- lengths(adj.list)
    weight.lengths <- lengths(weight.list)
    if (!identical(adj.lengths, weight.lengths)) {
        mismatches <- which(adj.lengths != weight.lengths)
        stop(sprintf(
            "Mismatched lengths between adj.list and weight.list at indices: %s",
            paste(mismatches, collapse = ", ")
        ), call. = FALSE)
    }

    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    n.vertices <- length(adj.list)
    if (length(y) != n.vertices) {
        stop(sprintf(
            "'y' length (%d) must match number of vertices in 'adj.list' (%d)",
            length(y), n.vertices
        ), call. = FALSE)
    }

    # Handle scale parameter
    if (!is.numeric(scale)) {
        stop("'scale' must be numeric", call. = FALSE)
    }

    if (length(scale) == 1L) {
        scale <- rep(scale, n.vertices)
    } else if (length(scale) != n.vertices) {
        stop(sprintf(
            "'scale' must be a single value or have length %d (same as 'y')",
            n.vertices
        ), call. = FALSE)
    }

    if (any(scale <= 0)) {
        stop("All values in 'scale' must be positive", call. = FALSE)
    }

    # Validate quantile threshold
    if (!is.numeric(quantile.scale.thld) || length(quantile.scale.thld) != 1L) {
        stop("'quantile.scale.thld' must be a single numeric value", call. = FALSE)
    }

    if (quantile.scale.thld <= 0 || quantile.scale.thld > 1) {
        stop("'quantile.scale.thld' must be in the interval (0, 1]", call. = FALSE)
    }

    if (!is.logical(with.trajectories) || length(with.trajectories) != 1L) {
        stop("'with.trajectories' must be a single logical value", call. = FALSE)
    }

    # Validate vertex indices
    all.indices <- unlist(adj.list)
    if (length(all.indices) > 0) {
        if (!all(all.indices >= 1L & all.indices <= n.vertices)) {
            invalid <- unique(all.indices[all.indices < 1L | all.indices > n.vertices])
            stop(sprintf(
                "Invalid vertex indices in 'adj.list': %s. Must be between 1 and %d.",
                paste(invalid, collapse = ", "), n.vertices
            ), call. = FALSE)
        }
    }

    # Check for non-negative weights
    all.weights <- unlist(weight.list)
    if (any(all.weights < 0)) {
        stop("All edge weights must be non-negative", call. = FALSE)
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))

    # Call C++ implementation
    result <- .Call(S_construct_graph_gradient_flow,
                    adj.list.0based,
                    weight.list,
                    y,
                    scale,
                    as.double(quantile.scale.thld),
                    as.logical(with.trajectories))

    # Process messages if present
    if (!is.null(result$messages) && nchar(result$messages) > 0L) {
        attr(result, "messages") <- result$messages
        result$messages <- NULL
    }

    # Enhanced labeling of extrema
    if (nrow(result$local_extrema) > 0L) {
        extrema_df <- .label_extrema(result$local_extrema, y)
        result$local_extrema <- extrema_df

        # Update labels in other components
        vertex_to_label <- setNames(extrema_df$label,
                                    as.character(extrema_df$vertex_index))

        result <- .update_component_labels(result, vertex_to_label, with.trajectories)
    }

    class(result) <- "ggflow"
    result
}


#' @keywords internal
.label_extrema <- function(extrema_matrix, y) {
    n_extrema <- nrow(extrema_matrix)
    extrema_indices <- extrema_matrix[, 1L]
    is_maximum <- extrema_matrix[, 2L] == 1L
    extrema_values <- y[extrema_indices]

    # Initialize labels
    extrema_labels <- character(n_extrema)

    # Label maxima (M1, M2, ...) in descending order
    if (any(is_maximum)) {
        max_idx <- which(is_maximum)
        max_order <- order(extrema_values[max_idx], decreasing = TRUE)
        extrema_labels[max_idx[max_order]] <- paste0("M", seq_along(max_order))
    }

    # Label minima (m1, m2, ...) in ascending order
    if (any(!is_maximum)) {
        min_idx <- which(!is_maximum)
        min_order <- order(extrema_values[min_idx])
        extrema_labels[min_idx[min_order]] <- paste0("m", seq_along(min_order))
    }

    # Create data frame
    data.frame(
        vertex_index = as.integer(extrema_indices),
        is_maximum = as.integer(is_maximum),
        label = extrema_labels,
        fn_value = extrema_values,
        stringsAsFactors = FALSE
    )
}


#' @keywords internal
.update_component_labels <- function(result, vertex_to_label, with.trajectories) {
    # Update basin labels
    if (length(result$basins$ascending) > 0L) {
        for (i in seq_along(result$basins$ascending)) {
            v <- as.character(result$basins$ascending[[i]]$local_min)
            result$basins$ascending[[i]]$label <- vertex_to_label[v]
        }
    }

    if (length(result$basins$descending) > 0L) {
        for (i in seq_along(result$basins$descending)) {
            v <- as.character(result$basins$descending[[i]]$local_max)
            result$basins$descending[[i]]$label <- vertex_to_label[v]
        }
    }

    # Update cell labels
    if (length(result$cells) > 0L) {
        for (i in seq_along(result$cells)) {
            min_v <- as.character(result$cells[[i]]$local_min)
            max_v <- as.character(result$cells[[i]]$local_max)
            min_label <- vertex_to_label[min_v]
            max_label <- vertex_to_label[max_v]

            result$cells[[i]]$min_label <- min_label
            result$cells[[i]]$max_label <- max_label
            result$cells[[i]]$label <- paste0(min_label, "-", max_label)
        }
    }

    # Update trajectory labels if present
    if (with.trajectories && length(result$trajectories) > 0L) {
        result <- .update_trajectory_labels(result, vertex_to_label)
    }

    result
}


#' @keywords internal
.update_trajectory_labels <- function(result, vertex_to_label) {
    for (i in seq_along(result$trajectories)) {
        traj <- result$trajectories[[i]]

        if (traj$type == "LMIN_LMAX") {
            first_v <- as.character(traj$vertices[1L])
            last_v <- as.character(traj$vertices[length(traj$vertices)])

            result$trajectories[[i]]$start_label <- vertex_to_label[first_v]
            result$trajectories[[i]]$end_label <- vertex_to_label[last_v]
            result$trajectories[[i]]$label <- paste0(
                vertex_to_label[first_v], "-", vertex_to_label[last_v]
            )
        } else if (traj$type == "LMIN_ONLY") {
            first_v <- as.character(traj$vertices[1L])
            result$trajectories[[i]]$start_label <- vertex_to_label[first_v]
            result$trajectories[[i]]$label <- paste0(vertex_to_label[first_v], "-*")
        } else if (traj$type == "LMAX_ONLY") {
            last_v <- as.character(traj$vertices[length(traj$vertices)])
            result$trajectories[[i]]$end_label <- vertex_to_label[last_v]
            result$trajectories[[i]]$label <- paste0("*-", vertex_to_label[last_v])
        }
    }
    result
}


#' Merge Two Basins in a Gradient Flow Complex
#'
#' @description
#' Merges two basins in a gradient flow complex by having one basin absorb another.
#' The function updates the basin vertices, recomputes the cells based on the
#' new basin configuration, and updates the local extrema table.
#'
#' @param flow An object of class \code{"ggflow"} returned by
#'   \code{\link{construct.graph.gradient.flow}}.
#' @param absorbing.label Character string specifying the label of the basin that
#'   will absorb the other basin.
#' @param absorbed.label Character string specifying the label of the basin that
#'   will be absorbed.
#'
#' @return A modified gradient flow object with merged basins and updated cells.
#'
#' @details
#' The function identifies the type of basins (ascending or descending) from the
#' label prefixes: "m" for minima (ascending basins) and "M" for maxima
#' (descending basins). The absorbing basin incorporates all vertices from the
#' absorbed basin, and the gradient flow structure is updated accordingly.
#'
#' @examples
#' \dontrun{
#' # Merge descending basins M2 and M3, with M2 absorbing M3
#' merged_flow <- basins.merge(flow, "M2", "M3")
#'
#' # Check the updated structure
#' summary(merged_flow)
#' }
#'
#' @seealso \code{\link{construct.graph.gradient.flow}}, \code{\link{summary.ggflow}}
#'
#' @export
basins.merge <- function(flow, absorbing.label, absorbed.label) {
    if (!inherits(flow, "ggflow")) {
        stop("'flow' must be an object of class 'ggflow'", call. = FALSE)
    }

    if (!is.character(absorbing.label) || length(absorbing.label) != 1L) {
        stop("'absorbing.label' must be a single character string", call. = FALSE)
    }

    if (!is.character(absorbed.label) || length(absorbed.label) != 1L) {
        stop("'absorbed.label' must be a single character string", call. = FALSE)
    }

    # Determine basin type
    basin_type <- .determine_basin_type(flow, absorbing.label, absorbed.label)
    basins <- flow$basins[[basin_type]]

    # Find basin indices
    absorbing_idx <- which(vapply(basins, function(b) b$label == absorbing.label,
                                  logical(1L)))
    absorbed_idx <- which(vapply(basins, function(b) b$label == absorbed.label,
                                 logical(1L)))

    if (length(absorbing_idx) == 0L) {
        stop(sprintf("Basin with label '%s' not found", absorbing.label),
             call. = FALSE)
    }

    if (length(absorbed_idx) == 0L) {
        stop(sprintf("Basin with label '%s' not found", absorbed.label),
             call. = FALSE)
    }

    # Get extremum information
    if (basin_type == "ascending") {
        absorbing_extremum <- basins[[absorbing_idx]]$local_min
        absorbed_extremum <- basins[[absorbed_idx]]$local_min
        extremum_type <- "minimum"
    } else {
        absorbing_extremum <- basins[[absorbing_idx]]$local_max
        absorbed_extremum <- basins[[absorbed_idx]]$local_max
        extremum_type <- "maximum"
    }

    # Update extrema table
    flow <- .update_extrema_for_merge(flow, absorbing_extremum, absorbed_extremum,
                                      absorbing.label, absorbed.label, extremum_type)

    # Merge basin vertices
    merged_vertices <- unique(c(basins[[absorbing_idx]]$vertices,
                                basins[[absorbed_idx]]$vertices))

    # Update the absorbing basin
    flow$basins[[basin_type]][[absorbing_idx]]$vertices <- merged_vertices

    # Add merge history
    if (is.null(flow$basins[[basin_type]][[absorbing_idx]]$merge_history)) {
        flow$basins[[basin_type]][[absorbing_idx]]$merge_history <- character()
    }

    flow$basins[[basin_type]][[absorbing_idx]]$merge_history <- c(
        flow$basins[[basin_type]][[absorbing_idx]]$merge_history,
        sprintf("Absorbed basin %s at vertex %d", absorbed.label, absorbed_extremum)
    )

    # Remove the absorbed basin
    flow$basins[[basin_type]] <- flow$basins[[basin_type]][-absorbed_idx]

    # Recompute cells
    flow$cells <- .recompute_cells(flow$basins$ascending, flow$basins$descending)

    # Record merger
    if (is.null(flow$basin_mergers)) {
        flow$basin_mergers <- list()
    }

    flow$basin_mergers <- c(flow$basin_mergers, list(list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        basin_type = basin_type,
        absorbing.label = absorbing.label,
        absorbing_vertex = absorbing_extremum,
        absorbed.label = absorbed.label,
        absorbed_vertex = absorbed_extremum
    )))

    flow
}


#' @keywords internal
.determine_basin_type <- function(flow, absorbing.label, absorbed.label) {
    # Try prefix-based detection
    absorbing_prefix <- substring(absorbing.label, 1L, 1L)
    absorbed_prefix <- substring(absorbed.label, 1L, 1L)

    if (absorbing_prefix == "m" && absorbed_prefix == "m") {
        return("ascending")
    } else if (absorbing_prefix == "M" && absorbed_prefix == "M") {
        return("descending")
    }

    # Fall back to checking actual basins
    asc_labels <- vapply(flow$basins$ascending, function(b) b$label, character(1L))
    desc_labels <- vapply(flow$basins$descending, function(b) b$label, character(1L))

    absorbing_in_asc <- absorbing.label %in% asc_labels
    absorbed_in_asc <- absorbed.label %in% asc_labels
    absorbing_in_desc <- absorbing.label %in% desc_labels
    absorbed_in_desc <- absorbed.label %in% desc_labels

    if (absorbing_in_asc && absorbed_in_asc) {
        return("ascending")
    } else if (absorbing_in_desc && absorbed_in_desc) {
        return("descending")
    } else {
        stop("Basin labels must both be in the same basin type", call. = FALSE)
    }
}


#' @keywords internal
.update_extrema_for_merge <- function(flow, absorbing_extremum, absorbed_extremum,
                                      absorbing.label, absorbed.label, extremum_type) {
    extrema_absorbed_idx <- which(flow$local_extrema$vertex_index == absorbed_extremum)
    extrema_absorbing_idx <- which(flow$local_extrema$vertex_index == absorbing_extremum)

    if (length(extrema_absorbed_idx) == 0L || length(extrema_absorbing_idx) == 0L) {
        warning("Could not find one or both extrema in the local_extrema table")
        return(flow)
    }

    # Add merge tracking columns if needed
    if (!"merged_into" %in% names(flow$local_extrema)) {
        flow$local_extrema$merged_into <- NA_integer_
    }

    if (!"is_merged" %in% names(flow$local_extrema)) {
        flow$local_extrema$is_merged <- FALSE
    }

    if (!"notes" %in% names(flow$local_extrema)) {
        flow$local_extrema$notes <- ""
    }

    # Update the absorbed extremum
    flow$local_extrema$merged_into[extrema_absorbed_idx] <- absorbing_extremum
    flow$local_extrema$is_merged[extrema_absorbed_idx] <- TRUE
    flow$local_extrema$label[extrema_absorbed_idx] <-
        flow$local_extrema$label[extrema_absorbing_idx]

    flow$local_extrema$notes[extrema_absorbed_idx] <- sprintf(
        "This %s was merged: basin %s was absorbed into basin %s",
        extremum_type, absorbed.label, absorbing.label
    )

    flow
}


#' @keywords internal
.recompute_cells <- function(ascending_basins, descending_basins) {
    cells <- list()

    for (asc_basin in ascending_basins) {
        for (desc_basin in descending_basins) {
            intersection <- intersect(asc_basin$vertices, desc_basin$vertices)

            if (length(intersection) > 0L) {
                cell <- list(
                    local_min = asc_basin$local_min,
                    local_max = desc_basin$local_max,
                    vertices = intersection,
                    min_label = asc_basin$label,
                    max_label = desc_basin$label,
                    label = paste0(asc_basin$label, "-", desc_basin$label)
                )
                cells <- c(cells, list(cell))
            }
        }
    }

    cells
}


#' Summarize Gradient Flow Structure
#'
#' @description
#' Provides a comprehensive summary of a gradient flow structure including
#' information about local extrema, basins, trajectories, cells, and
#' Jaccard indices between basins.
#'
#' @param object An object of class \code{"ggflow"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{"summary.ggflow"} containing:
#' \describe{
#'   \item{n_extrema}{Total number of extrema}
#'   \item{n_maxima}{Number of local maxima}
#'   \item{n_minima}{Number of local minima}
#'   \item{maxima}{Data frame with information about local maxima}
#'   \item{minima}{Data frame with information about local minima}
#'   \item{n_trajectories}{Number of trajectories (if available)}
#'   \item{trajectory_types}{Table of trajectory type counts}
#'   \item{avg_trajectory_length}{Average trajectory length}
#'   \item{connectivity}{Matrix showing min-max connections}
#'   \item{n_cells}{Number of gradient flow cells}
#'   \item{cell_sizes}{Summary statistics of cell sizes}
#'   \item{n_ascending_basins}{Number of ascending basins}
#'   \item{n_descending_basins}{Number of descending basins}
#'   \item{ascending_basins_jaccard_index}{Jaccard similarity matrix for ascending basins}
#'   \item{descending_basins_jaccard_index}{Jaccard similarity matrix for descending basins}
#' }
#'
#' @examples
#' \dontrun{
#' flow <- construct.graph.gradient.flow(adj.list, weight.list, y, scale)
#' summary(flow)
#' }
#'
#' @export
summary.ggflow <- function(object, ...) {
    if (!inherits(object, "ggflow")) {
        stop("'object' must be of class 'ggflow'", call. = FALSE)
    }

    result <- list()
    extrema <- object$local_extrema

    # Basic counts
    result$n_extrema <- nrow(extrema)
    result$n_maxima <- sum(extrema$is_maximum == 1L)
    result$n_minima <- result$n_extrema - result$n_maxima

    # Separate maxima and minima info
    if (result$n_extrema > 0L) {
        result$maxima <- extrema[extrema$is_maximum == 1L, ]
        result$minima <- extrema[extrema$is_maximum == 0L, ]
    } else {
        result$maxima <- result$minima <- data.frame()
    }

    # Trajectory analysis
    if (!is.null(object$trajectories)) {
        result <- .analyze_trajectories(result, object$trajectories,
                                        result$minima, result$maxima)
    } else {
        result$n_trajectories <- 0L
    }

    # Cell analysis
    result <- .analyze_cells(result, object$cells)

    # Basin analysis
    if (!is.null(object$basins)) {
        result <- .analyze_basins(result, object$basins, extrema)
    }

    result$messages <- attr(object, "messages")
    class(result) <- "summary.ggflow"
    result
}


#' @keywords internal
.analyze_trajectories <- function(result, trajectories, minima, maxima) {
    result$n_trajectories <- length(trajectories)

    # Trajectory types
    traj_types <- vapply(trajectories, function(t) t$type, character(1L))
    result$trajectory_types <- table(traj_types)

    # Length statistics
    traj_lengths <- vapply(trajectories, function(t) length(t$vertices), integer(1L))
    result$avg_trajectory_length <- mean(traj_lengths)
    result$max_trajectory_length <- max(traj_lengths)
    result$min_trajectory_length <- min(traj_lengths)

    # Connectivity matrix
    if (nrow(minima) > 0L && nrow(maxima) > 0L) {
        result$connectivity <- .compute_connectivity_matrix(
            trajectories, minima$label, maxima$label
        )
    }

    result
}


#' @keywords internal
.compute_connectivity_matrix <- function(trajectories, min_labels, max_labels) {
    conn <- matrix(0L, nrow = length(min_labels), ncol = length(max_labels))
    rownames(conn) <- min_labels
    colnames(conn) <- max_labels

    for (traj in trajectories) {
        if (traj$type == "LMIN_LMAX" &&
            !is.null(traj$start_label) &&
            !is.null(traj$end_label)) {

            min_idx <- which(min_labels == traj$start_label)
            max_idx <- which(max_labels == traj$end_label)

            if (length(min_idx) > 0L && length(max_idx) > 0L) {
                conn[min_idx, max_idx] <- conn[min_idx, max_idx] + 1L
            }
        }
    }

    conn
}


#' @keywords internal
.analyze_cells <- function(result, cells) {
    result$n_cells <- length(cells)

    if (result$n_cells > 0L) {
        cell_sizes <- vapply(cells, function(p) length(p$vertices), integer(1L))
        result$cell_sizes <- summary(cell_sizes)

        if ("label" %in% names(cells[[1L]])) {
            result$cell_connections <- data.frame(
                min_label = vapply(cells, function(p) p$min_label, character(1L)),
                max_label = vapply(cells, function(p) p$max_label, character(1L)),
                size = cell_sizes,
                stringsAsFactors = FALSE
            )
        }
    }

    result
}


#' @keywords internal
.analyze_basins <- function(result, basins, extrema) {
    result$n_ascending_basins <- length(basins$ascending)
    result$n_descending_basins <- length(basins$descending)

    # Get extrema vertices by type
    max_vertices <- extrema$vertex_index[extrema$is_maximum == 1L]
    min_vertices <- extrema$vertex_index[extrema$is_maximum == 0L]

    # Analyze ascending basins
    if (result$n_ascending_basins > 0L) {
        result <- .analyze_basin_set(result, basins$ascending, min_vertices,
                                     "ascending", "minima")
    }

    # Analyze descending basins
    if (result$n_descending_basins > 0L) {
        result <- .analyze_basin_set(result, basins$descending, max_vertices,
                                     "descending", "maxima")
    }

    result
}


#' @keywords internal
.analyze_basin_set <- function(result, basin_set, valid_vertices, basin_type,
                               extrema_type) {
    basin_sizes <- vapply(basin_set, function(b) length(b$vertices), integer(1L))
    result[[paste0(basin_type, "_basin_sizes")]] <- summary(basin_sizes)

    if ("label" %in% names(basin_set[[1L]])) {
        # Get basin info
        vertices <- vapply(basin_set,
                           function(b) b[[if(basin_type == "ascending") "local_min" else "local_max"]],
                           integer(1L))

        basin_info <- data.frame(
            vertex = vertices,
            label = vapply(basin_set, function(b) b$label, character(1L)),
            size = basin_sizes,
            stringsAsFactors = FALSE
        )

        # Filter valid basins
        is_valid <- vertices %in% valid_vertices
        basin_info <- basin_info[is_valid, ]

        # Add function values
        if (nrow(result[[extrema_type]]) > 0L) {
            vertex_to_value <- setNames(
                result[[extrema_type]]$fn_value,
                as.character(result[[extrema_type]]$vertex_index)
            )
            basin_info$function_value <- vertex_to_value[as.character(basin_info$vertex)]
        }

        result[[paste0(extrema_type, "_info")]] <- basin_info

        # Compute Jaccard indices
        valid_basins <- basin_set[is_valid]
        if (length(valid_basins) > 1L) {
            result[[paste0(basin_type, "_basins_jaccard_index")]] <-
                .compute_jaccard_matrix(valid_basins)
        }
    }

    result
}


#' @keywords internal
.compute_jaccard_matrix <- function(basins) {
    n <- length(basins)
    jaccard <- matrix(0, nrow = n, ncol = n)

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (i != j) {
                set1 <- basins[[i]]$vertices
                set2 <- basins[[j]]$vertices
                intersection_size <- length(intersect(set1, set2))
                union_size <- length(union(set1, set2))
                if (union_size > 0L) {
                    jaccard[i, j] <- intersection_size / union_size
                }
            }
        }
    }

    if ("label" %in% names(basins[[1L]])) {
        labels <- vapply(basins, function(b) b$label, character(1L))
        rownames(jaccard) <- colnames(jaccard) <- labels
    }

    jaccard
}


#' Print Gradient Flow Structure Summary
#'
#' @description
#' Print method for objects of class \code{"summary.ggflow"}. Displays a
#' formatted summary of a gradient flow structure.
#'
#' @param x An object of class \code{"summary.ggflow"}.
#' @param digits Number of significant digits for numeric output.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' flow <- construct.graph.gradient.flow(adj.list, weight.list, y, scale)
#' summary(flow)
#' }
#'
#' @export
print.summary.ggflow <- function(x, digits = 3L, ...) {
    cat("Gradient Flow Structure Summary\n")
    cat(strrep("=", 30), "\n\n")

    # Extrema summary
    cat("Extrema:\n")
    cat("  Total:", x$n_extrema, "\n")
    cat("  Local maxima:", x$n_maxima, "\n")
    cat("  Local minima:", x$n_minima, "\n\n")

    # Detailed extrema info
    if (!is.null(x$maxima_info) && nrow(x$maxima_info) > 0L) {
        cat("Local Maxima Information:\n")
        print(x$maxima_info[order(x$maxima_info$function_value, decreasing = TRUE), ],
              digits = digits)
        cat("\n")
    }

    if (!is.null(x$minima_info) && nrow(x$minima_info) > 0L) {
        cat("Local Minima Information:\n")
        print(x$minima_info[order(x$minima_info$function_value), ],
              digits = digits)
        cat("\n")
    }

    # Trajectory info
    if (!is.null(x$n_trajectories) && x$n_trajectories > 0L) {
        cat("Trajectories:\n")
        cat("  Total:", x$n_trajectories, "\n")

        if (!is.null(x$trajectory_types)) {
            cat("  Types:\n")
            for (type in names(x$trajectory_types)) {
                cat("    ", type, ": ", x$trajectory_types[type], "\n", sep = "")
            }
        }

        cat("  Length statistics:\n")
        cat("    Min:", x$min_trajectory_length, "\n")
        cat("    Avg:", round(x$avg_trajectory_length, digits), "\n")
        cat("    Max:", x$max_trajectory_length, "\n\n")
    }

    # Connectivity
    if (!is.null(x$connectivity)) {
        cat("Connectivity (minima to maxima):\n")
        print(x$connectivity)
        cat("\n")
    }

    # Cells
    cat("Cells:\n")
    cat("  Total:", x$n_cells, "\n")

    if (!is.null(x$cell_sizes)) {
        cat("  Size statistics:\n")
        print(x$cell_sizes)
        cat("\n")
    }

    # Basins
    cat("Basins:\n")
    cat("  Ascending basins:", x$n_ascending_basins, "\n")
    cat("  Descending basins:", x$n_descending_basins, "\n\n")

    # Jaccard indices
    .print_jaccard_matrix(x, "ascending", "minima", digits)
    .print_jaccard_matrix(x, "descending", "maxima", digits)

    # Messages
    if (!is.null(x$messages) && nchar(x$messages) > 0L) {
        cat("Messages:\n")
        cat("  ", x$messages, "\n\n", sep = "")
    }

    invisible(x)
}


#' @keywords internal
.print_jaccard_matrix <- function(x, basin_type, extrema_type, digits) {
    jaccard_name <- paste0(basin_type, "_basins_jaccard_index")

    if (!is.null(x[[jaccard_name]])) {
        cat(tools::toTitleCase(basin_type), "Basins Jaccard Index:\n")

        # Order by function values if available
        info_name <- paste0(extrema_type, "_info")
        if (!is.null(x[[info_name]]) && nrow(x[[info_name]]) > 0L) {
            order_idx <- if (extrema_type == "maxima") {
                order(x[[info_name]]$function_value, decreasing = TRUE)
            } else {
                order(x[[info_name]]$function_value)
            }

            ordered_labels <- x[[info_name]]$label[order_idx]
            jaccard_matrix <- x[[jaccard_name]][ordered_labels, ordered_labels]
            print(round(jaccard_matrix, digits))
        } else {
            print(round(x[[jaccard_name]], digits))
        }
        cat("\n")

        # Summary statistics
        jac_values <- x[[jaccard_name]][lower.tri(x[[jaccard_name]])]
        cat(tools::toTitleCase(basin_type), "Basins Jaccard Index Statistics:\n")
        print(summary(jac_values))
        cat("\n")
    }
}


#' Replace a Basin Label in Gradient Flow Object
#'
#' @description
#' Changes a specific basin label throughout a gradient flow object, ensuring
#' consistency across extrema table, basins, cells, and trajectories.
#'
#' @param flow An object of class \code{"ggflow"}.
#' @param old.label Character string specifying the label to be replaced.
#' @param new.label Character string specifying the new label.
#'
#' @return A modified gradient flow object with the updated label.
#'
#' @details
#' This function is useful for customizing basin labels after computation,
#' for example, to use more meaningful names based on domain knowledge.
#'
#' @examples
#' \dontrun{
#' # Replace a generic label with a meaningful one
#' flow_updated <- replace.basin.label(flow, "M1", "HighState")
#' }
#'
#' @export
replace.basin.label <- function(flow, old.label, new.label) {
    if (!inherits(flow, "ggflow")) {
        stop("'flow' must be an object of class 'ggflow'", call. = FALSE)
    }

    if (!is.character(old.label) || length(old.label) != 1L) {
        stop("'old.label' must be a single character string", call. = FALSE)
    }

    if (!is.character(new.label) || length(new.label) != 1L) {
        stop("'new.label' must be a single character string", call. = FALSE)
    }

    # Find the extremum with old.label
    extrema_idx <- which(flow$local_extrema$label == old.label)

    if (length(extrema_idx) == 0L) {
        stop(sprintf("Label '%s' not found in local extrema", old.label),
             call. = FALSE)
    }

    vertex_idx <- flow$local_extrema$vertex_index[extrema_idx]
    is_max <- flow$local_extrema$is_maximum[extrema_idx] == 1L
    basin_type <- if (is_max) "descending" else "ascending"

    # Update extrema table
    flow$local_extrema$label[extrema_idx] <- new.label

    # Update basins
    flow <- .update_basin_labels(flow, basin_type, vertex_idx, new.label)

    # Update cells
    flow <- .update_cell_labels(flow, old.label, new.label, is_max)

    # Update trajectories if present
    if (!is.null(flow$trajectories)) {
        flow <- .update_trajectory_labels_for_replacement(flow, old.label, new.label, is_max)
    }

    flow
}


#' @keywords internal
.update_basin_labels <- function(flow, basin_type, vertex_idx, new.label) {
    extremum_field <- if (basin_type == "ascending") "local_min" else "local_max"

    for (i in seq_along(flow$basins[[basin_type]])) {
        if (flow$basins[[basin_type]][[i]][[extremum_field]] == vertex_idx) {
            flow$basins[[basin_type]][[i]]$label <- new.label
            break
        }
    }

    flow
}


#' @keywords internal
.update_cell_labels <- function(flow, old.label, new.label, is_max) {
    for (i in seq_along(flow$cells)) {
        cell <- flow$cells[[i]]

        if (!is_max && !is.null(cell$min_label) && cell$min_label == old.label) {
            flow$cells[[i]]$min_label <- new.label
            flow$cells[[i]]$label <- paste0(new.label, "-", cell$max_label)
        } else if (is_max && !is.null(cell$max_label) && cell$max_label == old.label) {
            flow$cells[[i]]$max_label <- new.label
            flow$cells[[i]]$label <- paste0(cell$min_label, "-", new.label)
        }
    }

    flow
}


#' @keywords internal
.update_trajectory_labels_for_replacement <- function(flow, old.label, new.label, is_max) {
    for (i in seq_along(flow$trajectories)) {
        traj <- flow$trajectories[[i]]

        if (!is_max && !is.null(traj$start_label) && traj$start_label == old.label) {
            flow$trajectories[[i]]$start_label <- new.label

            if (!is.null(traj$label)) {
                if (grepl("-\\*$", traj$label)) {
                    flow$trajectories[[i]]$label <- paste0(new.label, "-*")
                } else if (!is.null(traj$end_label)) {
                    flow$trajectories[[i]]$label <- paste0(new.label, "-", traj$end_label)
                }
            }
        } else if (is_max && !is.null(traj$end_label) && traj$end_label == old.label) {
            flow$trajectories[[i]]$end_label <- new.label

            if (!is.null(traj$label)) {
                if (grepl("^\\*-", traj$label)) {
                    flow$trajectories[[i]]$label <- paste0("*-", new.label)
                } else if (!is.null(traj$start_label)) {
                    flow$trajectories[[i]]$label <- paste0(traj$start_label, "-", new.label)
                }
            }
        }
    }

    flow
}

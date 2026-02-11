# R Utility Functions for Gradient Basin Analysis
# ================================================

#' Extract Gradient Trajectory from Basin
#'
#' @description
#' Reconstructs the gradient trajectory from a vertex to the extremum
#' by following predecessor links backwards.
#'
#' @param basin A single basin structure from basins.of.attraction object
#' @param vertex.id Integer vertex identifier (1-based indexing)
#'
#' @return Integer vector representing the path from extremum to vertex.id,
#'   or NULL if vertex.id is not in the basin. The first element is the
#'   extremum, the last element is vertex.id.
#'
#' @examples
#' \dontrun{
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' # Extract trajectory from minimum at vertex 5 to vertex 23
#' traj <- extract.gradient.trajectory(basins$lmin_basins[[1]], 23)
#' # traj might be: c(5, 7, 12, 18, 23)
#' }
#'
#' @export
extract.gradient.trajectory <- function(basin, vertex.id) {
    # Check if vertex is in basin
    if (!(vertex.id %in% basin$basin_df[, 1])) {
        return(NULL)
    }
    
    trajectory <- integer()
    current <- vertex.id
    
    # Backtrack through predecessors
    while (current != 0) {
        trajectory <- c(current, trajectory)
        current <- basin$predecessors[current]
    }
    
    return(trajectory)
}


#' Extract All Trajectories Terminating at Specific Extremum
#'
#' @description
#' For a given basin and terminal extremum, extracts all basin vertices
#' that have trajectories terminating at that extremum.
#'
#' @param basin Basin structure (ascending or descending)
#' @param terminal.vertex Terminal extremum vertex (1-based)
#'
#' @return Integer vector of vertices whose trajectories terminate at terminal.vertex
#'
#' @export
extract.vertices.terminating.at <- function(basin, terminal.vertex) {
    if (!(terminal.vertex %in% basin$terminal_extrema)) {
        return(integer(0))
    }
    
    # Get all basin vertices
    basin.vertices <- basin$basin_df[, 1]
    
    # Check which ones have trajectories ending at terminal.vertex
    terminating.vertices <- integer()
    
    for (v in basin.vertices) {
        traj <- extract.gradient.trajectory(basin, v)
        if (!is.null(traj) && traj[length(traj)] == terminal.vertex) {
            terminating.vertices <- c(terminating.vertices, v)
        }
    }
    
    return(terminating.vertices)
}


#' Identify All Gradient Flow Cells
#'
#' @description
#' Constructs all valid gradient flow cells from basin structures. A cell
#' (m_i, M_j) exists if and only if M_j is reachable from m_i via ascending
#' trajectories and m_i is reachable from M_j via descending trajectories.
#'
#' @param x Object of class "basins_of_attraction" with terminal extrema
#' @param ... Additional arguments (currently ignored).
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item{min.vertex}: Local minimum vertex (1-based)
#'     \item{max.vertex}: Local maximum vertex (1-based)
#'     \item{min.value}: Function value at minimum
#'     \item{max.value}: Function value at maximum
#'     \item{cell.height}: max.value - min.value (monotonicity span)
#'     \item{min.basin.size}: Number of vertices in ascending basin
#'     \item{max.basin.size}: Number of vertices in descending basin
#'     \item{cell.size}: Number of vertices in cell intersection (computed separately)
#'   }
#'
#' @examples
#' \dontrun{
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' cells <- identify.gradient.flow.cells(basins)
#' 
#' # Find cells with largest monotonicity span
#' cells[order(-cells$cell.height), ]
#' 
#' # Find cells containing specific vertex
#' vertex.of.interest <- 42
#' cells[cells$min.vertex == vertex.of.interest | 
#'       cells$max.vertex == vertex.of.interest, ]
#' }
#'
#' @export
identify.gradient.flow.cells <- function(x, ...) {
    basins.obj <- x
    cells.list <- list()
    
    # Iterate through all minima
    for (i in seq_along(basins.obj$basins$lmin_basins)) {
        min.basin <- basins.obj$basins$lmin_basins[[i]]
        m.i <- min.basin$vertex
        min.basin.label <- names(basins.obj$basins$lmin_basins)[i]
        
        if (length(min.basin$terminal_extrema) == 0) next
        
        # For each terminal maximum of this minimum
        for (M.j in min.basin$terminal_extrema) {
            # Find the corresponding maximum basin
            max.basin.idx <- which(sapply(basins.obj$basins$lmax_basins,
                                         function(b) b$vertex == M.j))
            
            if (length(max.basin.idx) == 0) next
            
            max.basin <- basins.obj$basins$lmax_basins[[max.basin.idx]]
            max.basin.label <- names(basins.obj$basins$lmax_basins)[max.basin.idx]
            
            # Verify reciprocal reachability
            if (m.i %in% max.basin$terminal_extrema) {
                # Valid cell found
                cells.list[[length(cells.list) + 1]] <- data.frame(
                    min.label = min.basin.label,
                    max.label = max.basin.label,
                    min.vertex = m.i,
                    max.vertex = M.j,
                    min.value = min.basin$value,
                    max.value = max.basin$value,
                    cell.height = max.basin$value - min.basin$value,
                    min.basin.size = nrow(min.basin$basin_df),
                    max.basin.size = nrow(max.basin$basin_df),
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    
    if (length(cells.list) == 0) {
        return(data.frame(
            min.label = "",
            max.label = "",
            min.vertex = integer(),
            max.vertex = integer(),
            min.value = numeric(),
            max.value = numeric(),
            cell.height = numeric(),
            min.basin.size = integer(),
            max.basin.size = integer()
        ))
    }
    
    do.call(rbind, cells.list)
}


#' Extract Complete Gradient Trajectory Through Cell
#'
#' @description
#' Reconstructs a complete gradient trajectory through a cell (m_i, M_j)
#' passing through a specified interior vertex. The trajectory consists of
#' an ascending segment from m_i to the interior vertex, followed by a
#' descending segment from the interior vertex to M_j.
#'
#' @param basins.obj Basins object with terminal extrema
#' @param min.vertex Minimum vertex m_i (1-based)
#' @param max.vertex Maximum vertex M_j (1-based)
#' @param interior.vertex Vertex in cell (m_i, M_j) (1-based)
#'
#' @return List with components:
#'   \itemize{
#'     \item{ascending.segment}: Trajectory from m_i to interior.vertex
#'     \item{descending.segment}: Trajectory from interior.vertex to M_j
#'     \item{complete.trajectory}: Full path from m_i through interior.vertex to M_j
#'     \item{y.values}: Function values along complete trajectory
#'     \item{is.monotone}: Logical indicating if trajectory is monotone
#'   }
#'
#' @examples
#' \dontrun{
#' # Extract trajectory through cell (5, 10) passing through vertex 15
#' traj.info <- extract.cell.trajectory(basins, 
#'                                       min.vertex = 5, 
#'                                       max.vertex = 10, 
#'                                       interior.vertex = 15)
#' 
#' # Check monotonicity
#' if (!traj.info$is.monotone) {
#'   warning("Trajectory is not monotone - may contain spurious extrema")
#' }
#' 
#' # Plot trajectory
#' plot(seq_along(traj.info$y.values), traj.info$y.values, 
#'      type = "b", xlab = "Position", ylab = "Function value")
#' }
#'
#' @export
extract.cell.trajectory <- function(basins.obj, min.vertex, max.vertex, 
                                    interior.vertex) {
    # Find relevant basins
    min.basin.idx <- which(sapply(basins.obj$basins$lmin_basins,
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$basins$lmax_basins,
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        stop("Specified extrema not found in basins object")
    }
    
    min.basin <- basins.obj$basins$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$basins$lmax_basins[[max.basin.idx]]
    
    # Extract ascending trajectory: min -> interior
    asc.traj <- extract.gradient.trajectory(min.basin, interior.vertex)
    if (is.null(asc.traj)) {
        stop("Interior vertex not reachable from minimum")
    }
    
    # Extract descending trajectory: interior -> max
    desc.traj <- extract.gradient.trajectory(max.basin, interior.vertex)
    if (is.null(desc.traj)) {
        stop("Interior vertex not reachable from maximum")
    }
    
    # Reverse descending trajectory so it goes interior -> max
    desc.traj <- rev(desc.traj)
    
    # Combine trajectories (remove duplicate interior vertex)
    complete.traj <- c(asc.traj, desc.traj[-1])
    
    # Extract y values
    y.vals <- basins.obj$basins$y[complete.traj]
    
    # Check monotonicity
    # Should increase along ascending segment, decrease along descending segment
    asc.monotone <- all(diff(y.vals[seq_len(length(asc.traj))]) > 0)
    desc.monotone <- all(diff(y.vals[length(asc.traj):length(complete.traj)]) < 0)
    is.monotone <- asc.monotone && desc.monotone
    
    list(
        ascending.segment = asc.traj,
        descending.segment = desc.traj,
        complete.trajectory = complete.traj,
        y.values = y.vals,
        is.monotone = is.monotone
    )
}

#' Compute Cell Vertices (Intersection of Ascending and Descending Basins)
#'
#' @description
#' For a valid gradient flow cell (m_i, M_j), computes all vertices that
#' belong to both the ascending basin of m_i and the descending basin of M_j.
#'
#' @param basins.obj Basins object
#' @param min.vertex Minimum vertex m_i
#' @param max.vertex Maximum vertex M_j
#'
#' @return Integer vector of vertices in cell (m_i, M_j)
#'
#' @export
compute.cell.vertices <- function(basins.obj, min.vertex, max.vertex) {
    min.basin.idx <- which(sapply(basins.obj$basins$lmin_basins,
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$basins$lmax_basins,
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        return(integer(0))
    }
    
    min.basin <- basins.obj$basins$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$basins$lmax_basins[[max.basin.idx]]
    
    # Get basin vertex sets
    min.vertices <- min.basin$basin_df[, 1]
    max.vertices <- max.basin$basin_df[, 1]
    
    # Intersection
    intersect(min.vertices, max.vertices)
}


#' Validate Gradient Flow Cell
#'
#' @description
#' Checks if a pair (m_i, M_j) forms a valid gradient flow cell by verifying
#' bidirectional reachability.
#'
#' @param basins.obj Basins object
#' @param min.vertex Minimum vertex m_i
#' @param max.vertex Maximum vertex M_j
#'
#' @return Logical indicating if (m_i, M_j) is a valid cell
#'
#' @export
is.valid.cell <- function(basins.obj, min.vertex, max.vertex) {
    min.basin.idx <- which(sapply(basins.obj$basins$lmin_basins,
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$basins$lmax_basins,
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        return(FALSE)
    }
    
    min.basin <- basins.obj$basins$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$basins$lmax_basins[[max.basin.idx]]
    
    # Check bidirectional reachability
    max.in.min <- max.vertex %in% min.basin$terminal_extrema
    min.in.max <- min.vertex %in% max.basin$terminal_extrema
    
    max.in.min && min.in.max
}

#' Extract Gradient Flow Cells for a Local Maximum Basin
#'
#' Identifies all gradient flow cells within a local maximum basin, where each
#' cell corresponds to a terminal extremum and contains all vertices whose
#' gradient trajectories pass through that extremum on their way to the maximum.
#'
#' @param basins.obj Basin object containing the basin structure and
#'   predecessors vector from \code{compute.basins.of.attraction}
#' @param max.basin.label Character string identifying the local maximum basin
#'   (e.g., "M1", "M2")
#' @param verbose Logical indicating whether to print diagnostic information
#'   about the cell extraction process. Default is FALSE.
#'
#' @return A list with components:
#'   \item{max.basin.label}{The label of the maximum basin}
#'   \item{max.vertex}{The vertex index of the local maximum}
#'   \item{cells}{A named list where each element corresponds to a terminal
#'     extremum}
#'
#' @details
#' For a descending basin from a local maximum, the predecessor chain points
#' backward toward the maximum (i.e., \code{predecessors[v]} is the next vertex on
#' the path from v to the maximum). Terminal extrema are local minima where
#' gradient flow "starts" before flowing upward to the maximum.
#'
#' To construct gradient flow cells, we build a reverse graph (successors) where
#' \code{successors[u]} contains all vertices v such that \code{predecessors[v] = u}. Then,
#' for each terminal extremum, we perform an upstream search to find all vertices
#' whose path to the maximum passes through that terminal.
#'
#' @export
extract.gradient.flow.cells <- function(basins.obj,
                                        max.basin.label,
                                        verbose = FALSE) {

    if (!max.basin.label %in% names(basins.obj$basins$lmax_basins)) {
        stop("Basin label '", max.basin.label, "' not found in lmax_basins")
    }

    max.basin.idx <- which(names(basins.obj$basins$lmax_basins) == max.basin.label)
    max.basin <- basins.obj$basins$lmax_basins[[max.basin.idx]]

    if (is.null(max.basin$trajectory_sets) || length(max.basin$trajectory_sets) == 0) {
        stop("Trajectories not found in basin object.\n",
             "Please recompute basins with with.trajectories = TRUE")
    }

    max.vertex <- max.basin$vertex
    terminal.extrema <- max.basin$terminal_extrema
    basin.vertices <- as.integer(max.basin$basin_df[, 1])
    basin.size <- length(basin.vertices)

    absorbed.extrema <- if (!is.null(max.basin$absorbed_extrema)) {
                            max.basin$absorbed_extrema
                        } else {
                            integer(0)
                        }

    if (verbose) {
        cat("\n===== Extracting Gradient Flow Cells =====\n")
        cat("Basin:", max.basin.label, "\n")
        cat("Maximum vertex:", max.vertex, "\n")
        cat("Basin size:", basin.size, "vertices\n")
        cat("Number of terminal extrema:", length(terminal.extrema), "\n")
        ## cat("Terminal extrema:", paste(terminal.extrema, collapse = ", "), "\n")
        ## if (length(absorbed.extrema) > 0) {
        ##   cat("Absorbed extrema:", paste(absorbed.extrema, collapse = ", "), "\n")
        ##   cat("Note: Basin size may include vertices from absorbed extrema\n")
        ## }
        ## cat("\n")
    }

    ## Extract cells from pre-computed trajectory sets
    ## if (verbose) {
    ##     cat("Extracting cells from pre-computed trajectories...\n")
    ## }

    cells <- vector("list", length(max.basin$trajectory_sets))
    names(cells) <- paste0("cell_", sapply(max.basin$trajectory_sets,
                                           function(ts) ts$terminal_vertex))

    for (i in seq_along(max.basin$trajectory_sets)) {
        traj_set <- max.basin$trajectory_sets[[i]]
        term.vertex <- traj_set$terminal_vertex
        trajectories <- traj_set$trajectories

        ## All trajectories in this set start from the same terminal vertex
        ## Each trajectory is a path: [terminal, ..., maximum]
        ## Extract all unique vertices that appear in any trajectory
        cell.vertices <- unique(unlist(trajectories))
        cell.vertices <- sort(cell.vertices)

        ## Keep trajectories as a simple numbered list
        ## since they all start from the same terminal
        trajectories.list <- trajectories

        ## Build a mapping: which trajectories pass through each vertex?
        vertex.to.trajectories <- vector("list", max(basin.vertices))
        for (j in seq_along(trajectories)) {
            traj <- trajectories[[j]]
            for (v in traj) {
                vertex.to.trajectories[[v]] <- c(vertex.to.trajectories[[v]], j)
            }
        }

        cells[[i]] <- list(
            terminal.vertex = term.vertex,
            cell.vertices = cell.vertices,
            cell.size = length(cell.vertices),
            trajectories = trajectories.list,  # Simple list, all starting from terminal
            n.trajectories = length(trajectories),
            vertex.to.trajectories = vertex.to.trajectories[cell.vertices]  # For querying
        )

        ## if (verbose && i <= 5) {
        ##   cat(sprintf("  Cell %d (terminal %d): %d vertices, %d trajectories\n",
        ##               i, term.vertex, length(cell.vertices), length(trajectories)))
        ##   # Show trajectory length distribution
        ##   traj.lengths <- sapply(trajectories, length)
        ##   cat(sprintf("    Trajectory lengths: min=%d, median=%.0f, max=%d\n",
        ##               min(traj.lengths), median(traj.lengths), max(traj.lengths)))
        ## }
    }

    ## if (verbose && length(max.basin$trajectory_sets) > 5) {
    ##     cat("  ... (", length(max.basin$trajectory_sets) - 5,
    ##         " more cells processed)\n", sep = "")
    ## }

    ## Validation: check coverage and overlaps
    all.cell.vertices <- unique(unlist(lapply(cells, function(x) x$cell.vertices)))
    unassigned <- setdiff(basin.vertices, all.cell.vertices)

    ## Check for vertices appearing in multiple cells
    vertex.to.cells <- vector("list", max(basin.vertices))
    for (i in seq_along(cells)) {
        for (v in cells[[i]]$cell.vertices) {
            vertex.to.cells[[v]] <- c(vertex.to.cells[[v]], i)
        }
    }

    overlapping.vertices <- integer(0)
    for (v in basin.vertices) {
        if (length(vertex.to.cells[[v]]) > 1) {
            overlapping.vertices <- c(overlapping.vertices, v)
        }
    }

    if (verbose) {
        cat("\n===== Validation Summary =====\n")
        cat("Basin size:              ", basin.size, "\n")
        cat("Total assigned vertices: ", length(all.cell.vertices), "\n")
    }

    ## Build diagnostics
    diagnostics <- list(
        basin.size = basin.size,
        total.assigned = length(all.cell.vertices),
        size.match = (length(all.cell.vertices) == basin.size),
        n.unassigned = length(unassigned),
        n.overlapping = length(overlapping.vertices),
        has.absorbed.extrema = length(absorbed.extrema) > 0,
        n.absorbed.extrema = length(absorbed.extrema),
        total.trajectories = sum(sapply(cells, function(x) x$n.trajectories))
    )

    result <- list(
        basin.label = max.basin.label,
        max.vertex = max.vertex,
        cells = cells,
        unassigned.vertices = unassigned,
        overlapping.vertices = overlapping.vertices,
        absorbed.extrema = absorbed.extrema,
        diagnostics = diagnostics,
        n.cells = length(cells),
        total.assigned.vertices = length(all.cell.vertices)
    )

    class(result) <- c("gradient_flow_cells", "list")
    return(result)
}

#' Print Method for Gradient Flow Cells
#'
#' @param x An object of class "gradient_flow_cells"
#' @param ... Additional arguments (not used)
#'
#' @export
print.gradient_flow_cells <- function(x, ...) {
    cat("Gradient Flow Cells for Basin:", x$basin.label, "\n")
    cat("Maximum vertex:", x$max.vertex, "\n")
    cat("Number of cells:", x$n.cells, "\n")
    cat("Total assigned vertices:", x$total.assigned.vertices, "\n")

    if (x$diagnostics$n.unassigned > 0) {
        cat("WARNING: Unassigned vertices:", x$diagnostics$n.unassigned, "\n")
    }

    if (x$diagnostics$n.overlapping > 0) {
        cat("WARNING: Overlapping vertices:", x$diagnostics$n.overlapping, "\n")
    }

    ## if (x$diagnostics$has.absorbed.extrema) {
    ##     cat("Absorbed extrema:", x$diagnostics$n.absorbed.extrema, "\n")
    ## }

    ## if (!x$diagnostics$size.match) {
    ##     cat("WARNING: Cell sizes do not sum to basin size\n")
    ##     cat("  Basin size:", x$diagnostics$basin.size, "\n")
    ##     cat("  Sum of cells:", x$total.assigned.vertices, "\n")
    ## }

    cat("\nCell Summary:\n")
    for (i in seq_along(x$cells)) {
        cell <- x$cells[[i]]
        cat(sprintf("  %s: terminal vertex = %d, size = %d, trajectories = %d\n",
                    names(x$cells)[i], cell$terminal.vertex, cell$cell.size,
                    cell$n.trajectories))
    }

    cat("\nTotal trajectories:", x$diagnostics$total.trajectories, "\n")

    invisible(x)
}

#' Summary Method for Gradient Flow Cells
#'
#' @param object An object of class "gradient_flow_cells"
#' @param ... Additional arguments (not used)
#'
#' @export
summary.gradient_flow_cells <- function(object, ...) {
    cat("===== Gradient Flow Cells Summary =====\n\n")
    cat("Basin:", object$basin.label, "\n")
    cat("Maximum vertex:", object$max.vertex, "\n\n")

    cat("Structure:\n")
    cat("  Number of cells:", object$n.cells, "\n")
    cat("  Basin size:", object$diagnostics$basin.size, "vertices\n")
    cat("  Assigned vertices:", object$total.assigned.vertices, "\n")
    cat("  Unassigned vertices:", object$diagnostics$n.unassigned, "\n")
    cat("  Overlapping vertices:", object$diagnostics$n.overlapping, "\n")
    cat("  Total trajectories:", object$diagnostics$total.trajectories, "\n\n")

    ## if (object$diagnostics$has.absorbed.extrema) {
    ##     cat("Merged structure:\n")
    ##     cat("  Absorbed extrema:", object$diagnostics$n.absorbed.extrema, "\n")
    ##     if (length(object$absorbed.extrema) > 0) {
    ##         cat("  Absorbed vertex indices:",
    ##             paste(head(object$absorbed.extrema, 10), collapse = ", "))
    ##         if (length(object$absorbed.extrema) > 10) cat(", ...")
    ##         cat("\n")
    ##     }
    ##     cat("\n")
    ## }

    cat("Cell details:\n")
    for (i in seq_along(object$cells)) {
        cell <- object$cells[[i]]
        cat(sprintf("  Cell %d: terminal = %d\n", i, cell$terminal.vertex))
        cat(sprintf("    Vertices: %d\n", cell$cell.size))
        cat(sprintf("    Trajectories: %d\n", cell$n.trajectories))

        if (length(cell$trajectories) > 0) {
            traj.lengths <- sapply(cell$trajectories, length)
            traj.lengths <- traj.lengths[traj.lengths > 0]
            if (length(traj.lengths) > 0) {
                cat(sprintf("    Trajectory lengths: min = %d, median = %.1f, max = %d\n",
                            min(traj.lengths), median(traj.lengths), max(traj.lengths)))
            }
        }
        cat("\n")
    }

    if (!object$diagnostics$size.match) {
        cat("Validation:\n")
        cat("  WARNING: Size mismatch detected\n")
        cat("  Basin size:", object$diagnostics$basin.size, "\n")
        cat("  Sum of cells:", object$total.assigned.vertices, "\n")
        cat("  Difference:",
            object$diagnostics$basin.size - object$total.assigned.vertices, "\n\n")
    } else {
        cat("Validation: PASSED (perfect partition)\n\n")
    }

    cat("======================================\n")
    invisible(object)
}

#' Filter Gradient Flow Cells by Vertex Count
#'
#' Filters gradient flow cells to retain only those with at least a specified
#' number of vertices. This is useful for focusing analysis on cells with
#' sufficient statistical power or biological significance.
#'
#' @param flow.cells A gradient flow cells object returned by extract.gradient.flow.cells()
#' @param min.vertices Minimum number of vertices required for a cell to be retained
#' @param verbose If TRUE, print filtering summary
#'
#' @return A filtered gradient flow cells object with the same structure as the input,
#'         containing only cells meeting the vertex count threshold. The object retains:
#'         - basin.label: Basin identifier
#'         - max.vertex: Reference extremum vertex
#'         - basin.size: Original basin size (unchanged)
#'         - cells: List of filtered cells
#'         - n.cells: Number of cells after filtering
#'         - assigned.vertices: Set of vertices in retained cells
#'         - n.assigned: Count of assigned vertices after filtering
#'         - diagnostics: Updated diagnostic information
#'
#' @details The function preserves all trajectory information for retained cells
#'          and recalculates coverage statistics based on the filtered cell set.
#'          Cells are removed if they contain fewer than min.vertices unique vertices.
#'
#' @examples
#' \dontrun{
#' ## Extract cells from a basin
#' flow.cells <- extract.gradient.flow.cells(basins.obj, max.basin.label = "M1")
#'
#' ## Keep only cells with at least 50 vertices
#' large.cells <- filter.cells.by.vertex.count(flow.cells, min.vertices = 50)
#'
#' ## Keep only cells with at least 100 vertices, with verbose output
#' very.large.cells <- filter.cells.by.vertex.count(
#'   flow.cells,
#'   min.vertices = 100,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
filter.cells.by.vertex.count <- function(flow.cells,
                                         min.vertices = 50,
                                         verbose = FALSE) {
    ## Input validation
    if (!inherits(flow.cells, "gradient_flow_cells")) {
        stop("flow.cells must be a gradient_flow_cells object from extract.gradient.flow.cells()")
    }

    if (!is.numeric(min.vertices) || min.vertices < 1) {
        stop("min.vertices must be a positive number")
    }

    ## Store original counts for reporting
    original.n.cells <- flow.cells$n.cells
    original.n.assigned <- flow.cells$n.assigned

    ## Filter cells by vertex count
    filtered.cells <- list()
    retained.terminal.vertices <- c()
    basin.size <- 0
    for (i in seq_along(flow.cells$cells)) {
        cell <- flow.cells$cells[[i]]
        basin.size <- basin.size + cell$cell.size
        if (cell$cell.size >= min.vertices) {
            filtered.cells[[length(filtered.cells) + 1]] <- cell
            retained.terminal.vertices <- c(retained.terminal.vertices, cell$terminal.vertex)
        }
    }

    ## Update cell count
    n.filtered.cells <- length(filtered.cells)

    ## Recalculate assigned vertices from retained cells
    filtered.assigned.vertices <- integer(0)
    for (cell in filtered.cells) {
        filtered.assigned.vertices <- union(filtered.assigned.vertices, cell$cell.vertices)
    }

    n.filtered.assigned <- length(filtered.assigned.vertices)

    ## Calculate unassigned vertices
    n.unassigned <- basin.size - n.filtered.assigned

    ## Calculate overlapping vertices (vertices in multiple cells)
    total.vertex.occurrences <- sum(sapply(filtered.cells, function(x) x$cell.size))
    n.overlapping <- total.vertex.occurrences - n.filtered.assigned

    ## Update diagnostics
    diagnostics <- list(
        n.assigned = n.filtered.assigned,
        n.unassigned = n.unassigned,
        n.overlapping = n.overlapping,
        is.partition = (n.unassigned == 0 && n.overlapping == 0)
    )

    ## Create filtered result object
    result <- list(
        basin.label = flow.cells$basin.label,
        max.vertex = flow.cells$max.vertex,
        basin.size = flow.cells$basin.size,
        cells = filtered.cells,
        n.cells = n.filtered.cells,
        assigned.vertices = filtered.assigned.vertices,
        n.assigned = n.filtered.assigned,
        diagnostics = diagnostics
    )

    class(result) <- c("gradient_flow_cells", "list")

    ## Report filtering results
    if (verbose) {
        cat("\n===== Cell Filtering Summary =====\n")
        cat(sprintf("Basin: %s\n", flow.cells$basin.label))
        cat(sprintf("Minimum vertices threshold: %d\n", min.vertices))
        cat(sprintf("\nOriginal structure:\n"))
        cat(sprintf("  Total cells: %d\n", original.n.cells))
        cat(sprintf("  Assigned vertices: %d/%d (%.1f%%)\n",
                    original.n.assigned,
                    flow.cells$basin.size,
                    100.0 * original.n.assigned / flow.cells$basin.size))
        cat(sprintf("\nFiltered structure:\n"))
        cat(sprintf("  Retained cells: %d (%.1f%%)\n",
                    n.filtered.cells,
                    100.0 * n.filtered.cells / original.n.cells))
        cat(sprintf("  Removed cells: %d\n", original.n.cells - n.filtered.cells))
        cat(sprintf("  Assigned vertices: %d/%d (%.1f%%)\n",
                    n.filtered.assigned,
                    flow.cells$basin.size,
                    100.0 * n.filtered.assigned / flow.cells$basin.size))
        cat(sprintf("  Unassigned vertices: %d\n", n.unassigned))
        cat(sprintf("  Vertex reduction: %d vertices (%.1f%%)\n",
                    original.n.assigned - n.filtered.assigned,
                    100.0 * (original.n.assigned - n.filtered.assigned) / original.n.assigned))

        if (n.filtered.cells > 0) {
            cell.sizes <- sapply(filtered.cells, function(x) x$cell.size)
            cat(sprintf("\nRetained cell sizes:\n"))
            cat(sprintf("  Min: %d vertices\n", min(cell.sizes)))
            cat(sprintf("  Median: %.0f vertices\n", median(cell.sizes)))
            cat(sprintf("  Max: %d vertices\n", max(cell.sizes)))
        }
        cat("==================================\n\n")
    }

    return(result)
}


#' Filter Gradient Flow Cells by Multiple Criteria
#'
#' Advanced filtering function that allows multiple criteria including
#' minimum/maximum vertex count, minimum trajectory count, and vertex count percentiles.
#'
#' @param flow.cells A gradient flow cells object
#' @param min.vertices Minimum number of vertices (default: NULL, no minimum)
#' @param max.vertices Maximum number of vertices (default: NULL, no maximum)
#' @param min.trajectories Minimum number of trajectories (default: NULL, no minimum)
#' @param min.percentile Keep cells above this percentile of cell sizes (default: NULL)
#' @param max.percentile Keep cells below this percentile of cell sizes (default: NULL)
#' @param verbose If TRUE, print filtering summary
#'
#' @return Filtered gradient flow cells object
#'
#' @examples
#' \dontrun{
#' ## Keep cells between 50 and 200 vertices
#' medium.cells <- filter.cells.by.criteria(
#'   flow.cells,
#'   min.vertices = 50,
#'   max.vertices = 200
#' )
#'
#' ## Keep top 25% largest cells
#' top.cells <- filter.cells.by.criteria(
#'   flow.cells,
#'   min.percentile = 0.75,
#'   verbose = TRUE
#' )
#'
#' ## Keep cells with at least 50 vertices and 10 trajectories
#' substantial.cells <- filter.cells.by.criteria(
#'   flow.cells,
#'   min.vertices = 50,
#'   min.trajectories = 10
#' )
#' }
#'
#' @export
filter.cells.by.criteria <- function(flow.cells,
                                     min.vertices = NULL,
                                     max.vertices = NULL,
                                     min.trajectories = NULL,
                                     min.percentile = NULL,
                                     max.percentile = NULL,
                                     verbose = FALSE) {
    ## Input validation
    if (!inherits(flow.cells, "gradient_flow_cells")) {
        stop("flow.cells must be a gradient_flow_cells object")
    }

    ## Calculate cell sizes and trajectory counts
    cell.sizes <- sapply(flow.cells$cells, function(x) x$cell.size)
    n.trajectories <- sapply(flow.cells$cells, function(x) x$n.trajectories)

    ## Initialize filter mask (all TRUE)
    keep.mask <- rep(TRUE, length(flow.cells$cells))

    ## Apply filters
    if (!is.null(min.vertices)) {
        keep.mask <- keep.mask & (cell.sizes >= min.vertices)
    }

    if (!is.null(max.vertices)) {
        keep.mask <- keep.mask & (cell.sizes <= max.vertices)
    }

    if (!is.null(min.trajectories)) {
        keep.mask <- keep.mask & (n.trajectories >= min.trajectories)
    }

    if (!is.null(min.percentile)) {
        threshold <- quantile(cell.sizes, min.percentile)
        keep.mask <- keep.mask & (cell.sizes >= threshold)
    }

    if (!is.null(max.percentile)) {
        threshold <- quantile(cell.sizes, max.percentile)
        keep.mask <- keep.mask & (cell.sizes <= threshold)
    }

    ## Filter cells
    filtered.cells <- flow.cells$cells[keep.mask]

    ## Recalculate assigned vertices
    filtered.assigned.vertices <- integer(0)
    for (cell in filtered.cells) {
        filtered.assigned.vertices <- union(filtered.assigned.vertices, cell$cell.vertices)
    }

    n.filtered.assigned <- length(filtered.assigned.vertices)

    ## Calculate unassigned and overlapping vertices
    n.unassigned <- flow.cells$basin.size - n.filtered.assigned
    total.vertex.occurrences <- sum(sapply(filtered.cells, function(x) x$cell.size))
    n.overlapping <- total.vertex.occurrences - n.filtered.assigned

    ## Update diagnostics
    diagnostics <- list(
        n.assigned = n.filtered.assigned,
        n.unassigned = n.unassigned,
        n.overlapping = n.overlapping,
        is.partition = (n.unassigned == 0 && n.overlapping == 0)
    )

    ## Create result
    result <- list(
        basin.label = flow.cells$basin.label,
        max.vertex = flow.cells$max.vertex,
        basin.size = flow.cells$basin.size,
        cells = filtered.cells,
        n.cells = length(filtered.cells),
        assigned.vertices = filtered.assigned.vertices,
        n.assigned = n.filtered.assigned,
        diagnostics = diagnostics
    )

    class(result) <- c("gradient_flow_cells", "list")

    ## Verbose reporting
    if (verbose) {
        cat("\n===== Multi-Criteria Cell Filtering =====\n")
        cat(sprintf("Basin: %s\n", flow.cells$basin.label))
        cat("\nCriteria applied:\n")
        if (!is.null(min.vertices)) cat(sprintf("  Minimum vertices: %d\n", min.vertices))
        if (!is.null(max.vertices)) cat(sprintf("  Maximum vertices: %d\n", max.vertices))
        if (!is.null(min.trajectories)) cat(sprintf("  Minimum trajectories: %d\n", min.trajectories))
        if (!is.null(min.percentile)) cat(sprintf("  Minimum percentile: %.1f%%\n", 100 * min.percentile))
        if (!is.null(max.percentile)) cat(sprintf("  Maximum percentile: %.1f%%\n", 100 * max.percentile))

        cat(sprintf("\nResults:\n"))
        cat(sprintf("  Original cells: %d\n", flow.cells$n.cells))
        cat(sprintf("  Retained cells: %d (%.1f%%)\n",
                    result$n.cells,
                    100.0 * result$n.cells / flow.cells$n.cells))
        cat(sprintf("  Removed cells: %d\n", flow.cells$n.cells - result$n.cells))

        if (result$n.cells > 0) {
            retained.sizes <- sapply(result$cells, function(x) x$cell.size)
            cat(sprintf("\nRetained cell statistics:\n"))
            cat(sprintf("  Size range: %d - %d vertices\n", min(retained.sizes), max(retained.sizes)))
            cat(sprintf("  Median size: %.0f vertices\n", median(retained.sizes)))
            cat(sprintf("  Total vertices: %d (%.1f%% of basin)\n",
                        result$n.assigned,
                        100.0 * result$n.assigned / flow.cells$basin.size))
        }
        cat("==========================================\n\n")
    }

    return(result)
}

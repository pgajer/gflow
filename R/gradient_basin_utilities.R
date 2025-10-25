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
#' @param basins.obj Object of class "basins_of_attraction" with terminal extrema
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
identify.gradient.flow.cells <- function(basins.obj) {
    cells.list <- list()
    
    # Iterate through all minima
    for (i in seq_along(basins.obj$lmin_basins)) {
        min.basin <- basins.obj$lmin_basins[[i]]
        m.i <- min.basin$vertex
        
        if (length(min.basin$terminal_extrema) == 0) next
        
        # For each terminal maximum of this minimum
        for (M.j in min.basin$terminal_extrema) {
            # Find the corresponding maximum basin
            max.basin.idx <- which(sapply(basins.obj$lmax_basins, 
                                         function(b) b$vertex == M.j))
            
            if (length(max.basin.idx) == 0) next
            
            max.basin <- basins.obj$lmax_basins[[max.basin.idx]]
            
            # Verify reciprocal reachability
            if (m.i %in% max.basin$terminal_extrema) {
                # Valid cell found
                cells.list[[length(cells.list) + 1]] <- data.frame(
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
    min.basin.idx <- which(sapply(basins.obj$lmin_basins, 
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$lmax_basins, 
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        stop("Specified extrema not found in basins object")
    }
    
    min.basin <- basins.obj$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$lmax_basins[[max.basin.idx]]
    
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
    y.vals <- basins.obj$y[complete.traj]
    
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
    min.basin.idx <- which(sapply(basins.obj$lmin_basins, 
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$lmax_basins, 
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        return(integer(0))
    }
    
    min.basin <- basins.obj$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$lmax_basins[[max.basin.idx]]
    
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
    min.basin.idx <- which(sapply(basins.obj$lmin_basins, 
                                  function(b) b$vertex == min.vertex))
    max.basin.idx <- which(sapply(basins.obj$lmax_basins, 
                                  function(b) b$vertex == max.vertex))
    
    if (length(min.basin.idx) == 0 || length(max.basin.idx) == 0) {
        return(FALSE)
    }
    
    min.basin <- basins.obj$lmin_basins[[min.basin.idx]]
    max.basin <- basins.obj$lmax_basins[[max.basin.idx]]
    
    # Check bidirectional reachability
    max.in.min <- max.vertex %in% min.basin$terminal_extrema
    min.in.max <- min.vertex %in% max.basin$terminal_extrema
    
    max.in.min && min.in.max
}

#' Compute Harmonic Extension of Trajectory Coordinates
#'
#' Extends arc-length coordinates from a geodesic trajectory to a tubular
#' neighborhood via the discrete harmonic (Laplace) equation. The resulting
#' coordinates provide a smooth parameterization of the neighborhood that
#' respects the graph structure.
#'
#' @param adj.list Graph adjacency list (1-based indexing).
#' @param weight.list Edge weight list (edge lengths).
#' @param trajectory Integer vector of trajectory vertices in order from
#'   minimum to maximum (1-based indexing).
#' @param tube.radius Integer; radius of tubular neighborhood in hops.
#' @param use.edge.weights Logical; if \code{TRUE} (default), use inverse
#'   edge length as Laplacian weights. If \code{FALSE}, use unit weights.
#' @param max.iterations Integer; maximum iterations for Gauss-Seidel solver.
#' @param tolerance Numeric; convergence tolerance for solver.
#' @param basin.restriction Optional integer vector of vertices to restrict
#'   the tubular neighborhood. If \code{NULL}, no restriction is applied.
#' @param verbose Logical; print progress information.
#'
#' @return An object of class \code{"harmonic_extension"} containing:
#'   \item{trajectory}{Trajectory vertices (1-based).}
#'   \item{trajectory.coords}{Arc-length coordinates for trajectory vertices.}
#'   \item{trajectory.length}{Total geodesic length of trajectory.}
#'   \item{tubular.vertices}{All vertices in tubular neighborhood (1-based).}
#'   \item{hop.distances}{Hop distance from trajectory for each tubular vertex.}
#'   \item{nearest.traj.idx}{Index (1-based) into trajectory of the nearest
#'     trajectory vertex for each tubular vertex.}
#'   \item{extended.coords}{Harmonic extension coordinates for all tubular
#'     vertices, in same order as tubular.vertices.}
#'   \item{n.iterations}{Number of solver iterations.}
#'   \item{final.max.change}{Final maximum coordinate change (convergence).}
#'   \item{tube.radius}{Tube radius used.}
#'
#' @details
#' The harmonic extension is computed by solving the discrete Laplace equation
#' with Dirichlet boundary conditions. Trajectory vertices have fixed
#' coordinates determined by their arc-length position along the path.
#' Non-trajectory vertices in the tubular neighborhood receive coordinates
#' via the unique harmonic function that matches the boundary conditions.
#'
#' The harmonic extension minimizes the Dirichlet energy (sum of squared
#' differences across edges) and satisfies the maximum principle: all
#' extended coordinates lie within \eqn{[0, 1]}.
#'
#' @examples
#' \dontrun{
#' # Extract trajectory for a cell
#' cell.traj <- cell.trajectories(gfc, "m4", "M1")
#'
#' # Select trajectory with maximal mean density
#' best.idx <- select.max.density.trajectory(
#'     cell.traj$trajectories,
#'     smoothed.rho
#' )
#' best.traj <- cell.traj$trajectories[[best.idx]]
#'
#' # Compute harmonic extension
#' hext <- compute.harmonic.extension(
#'     adj.list, weight.list,
#'     trajectory = best.traj,
#'     tube.radius = 3,
#'     verbose = TRUE
#' )
#'
#' # Access extended coordinates
#' coords <- data.frame(
#'     vertex = hext$tubular.vertices,
#'     coord = hext$extended.coords,
#'     hop = hext$hop.distances
#' )
#' }
#'
#' @seealso \code{\link{select.max.density.trajectory}},
#'   \code{\link{cell.trajectories}}
#'
#' @export
compute.harmonic.extension <- function(adj.list,
                                        weight.list,
                                        trajectory,
                                        tube.radius = 2L,
                                        use.edge.weights = TRUE,
                                        max.iterations = 1000L,
                                        tolerance = 1e-8,
                                        basin.restriction = NULL,
                                        verbose = TRUE) {

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }

    if (!is.numeric(trajectory) || length(trajectory) < 2) {
        stop("trajectory must be a numeric vector with at least 2 vertices")
    }

    trajectory <- as.integer(trajectory)

    if (any(trajectory < 1) || any(trajectory > n.vertices)) {
        stop("trajectory vertices must be between 1 and n.vertices")
    }

    tube.radius <- as.integer(tube.radius)
    if (tube.radius < 1) {
        stop("tube.radius must be at least 1")
    }

    ## ========================================================================
    ## Convert to 0-based indexing for C++
    ## ========================================================================

    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            integer(0)
        } else {
            as.integer(x - 1L)
        }
    })

    trajectory.0based <- trajectory - 1L

    basin.restriction.0based <- NULL
    if (!is.null(basin.restriction)) {
        basin.restriction.0based <- as.integer(basin.restriction - 1L)
    }

    ## ========================================================================
    ## Call C++
    ## ========================================================================

    result <- .Call(
        S_compute_harmonic_extension,
        adj.list.0based,
        weight.list,
        trajectory.0based,
        tube.radius,
        as.logical(use.edge.weights),
        as.integer(max.iterations),
        as.numeric(tolerance),
        basin.restriction.0based,
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    return(result)
}


#' Select Trajectory with Maximal Mean Density
#'
#' Given a collection of trajectories (e.g., all paths in a gradient flow cell),
#' selects the one with the highest mean vertex density. This trajectory passes
#' through the densest region of the cell and is suitable for harmonic extension.
#'
#' @param trajectories List of integer vectors, each representing a trajectory
#'   as a sequence of vertex indices (1-based).
#' @param density Numeric vector of density values for all vertices.
#'
#' @return Integer index (1-based) of the trajectory with maximal mean density.
#'
#' @examples
#' \dontrun{
#' # Get cell trajectories
#' cell.traj <- cell.trajectories(gfc, "m4", "M1")
#'
#' # Compute smoothed density
#' smoothed.rho <- compute.smoothed.density(fit)
#'
#' # Select best trajectory
#' best.idx <- select.max.density.trajectory(
#'     cell.traj$trajectories,
#'     smoothed.rho
#' )
#'
#' # Extract the selected trajectory
#' best.traj <- cell.traj$trajectories[[best.idx]]
#' }
#'
#' @export
select.max.density.trajectory <- function(trajectories, density) {

    if (!is.list(trajectories)) {
        stop("trajectories must be a list")
    }

    if (length(trajectories) == 0) {
        return(NA_integer_)
    }

    if (!is.numeric(density)) {
        stop("density must be a numeric vector")
    }

    ## Convert trajectories to 0-based for C++
    trajectories.0based <- lapply(trajectories, function(traj) {
        as.integer(traj - 1L)
    })

    result <- .Call(
        S_select_max_density_trajectory,
        trajectories.0based,
        as.numeric(density),
        PACKAGE = "gflow"
    )

    return(result)
}


#' Print Method for harmonic_extension Objects
#'
#' @param x A \code{harmonic_extension} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.harmonic_extension <- function(x, ...) {

    cat("Harmonic Extension of Trajectory Coordinates\n")
    cat("=============================================\n")
    cat(sprintf("Trajectory: %d vertices, length %.4f\n",
                length(x$trajectory), x$trajectory.length))
    cat(sprintf("Tubular neighborhood: %d vertices (radius %d hops)\n",
                length(x$tubular.vertices), x$tube.radius))
    cat(sprintf("Solver: %d iterations, final change %.2e\n",
                x$n.iterations, x$final.max.change))

    ## Coordinate summary
    cat(sprintf("\nExtended coordinates:\n"))
    cat(sprintf("  Range: [%.4f, %.4f]\n",
                min(x$extended.coords), max(x$extended.coords)))
    cat(sprintf("  Mean: %.4f, SD: %.4f\n",
                mean(x$extended.coords), sd(x$extended.coords)))

    ## Count by hop distance
    hop.table <- table(x$hop.distances)
    cat(sprintf("\nVertices by hop distance:\n"))
    for (h in names(hop.table)) {
        cat(sprintf("  Hop %s: %d vertices\n", h, hop.table[h]))
    }

    invisible(x)
}


#' Get Coordinate for a Specific Vertex
#'
#' @param hext A \code{harmonic_extension} object.
#' @param vertex Vertex index (1-based).
#'
#' @return The extended coordinate for the vertex, or \code{NA} if the vertex
#'   is not in the tubular neighborhood.
#'
#' @export
get.extended.coord <- function(hext, vertex) {

    if (!inherits(hext, "harmonic_extension")) {
        stop("hext must be a harmonic_extension object")
    }

    idx <- match(vertex, hext$tubular.vertices)

    if (is.na(idx)) {
        return(NA_real_)
    }

    return(hext$extended.coords[idx])
}


#' Create Data Frame of Extended Coordinates
#'
#' Extracts the harmonic extension results as a data frame suitable for
#' downstream analysis.
#'
#' @param hext A \code{harmonic_extension} object.
#' @param y Optional numeric vector of response values (e.g., fitted values).
#' @param density Optional numeric vector of density values.
#'
#' @return A data frame with columns:
#'   \item{vertex}{Vertex index (1-based)}
#'   \item{coord}{Extended coordinate in \eqn{[0, 1]}}
#'   \item{hop}{Hop distance from trajectory}
#'   \item{on.trajectory}{Logical; TRUE if vertex is on the trajectory}
#'   \item{y}{Response value (if provided)}
#'   \item{density}{Density value (if provided)}
#'
#' @export
as.data.frame.harmonic_extension <- function(hext, y = NULL, density = NULL) {

    df <- data.frame(
        vertex = hext$tubular.vertices,
        coord = hext$extended.coords,
        hop = hext$hop.distances,
        on.trajectory = hext$tubular.vertices %in% hext$trajectory
    )

    if (!is.null(y)) {
        df$y <- y[hext$tubular.vertices]
    }

    if (!is.null(density)) {
        df$density <- density[hext$tubular.vertices]
    }

    ## Sort by coordinate
    df <- df[order(df$coord), ]
    rownames(df) <- NULL

    return(df)
}

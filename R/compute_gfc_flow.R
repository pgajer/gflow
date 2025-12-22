#' Compute Gradient Flow Complex Using Trajectory-First Approach
#'
#' Computes basins of attraction by tracing gradient flow trajectories rather
#' than growing basins outward from extrema. This approach treats gradient flow
#' lines as the fundamental geometric primitives, with basins emerging as the
#' collection of vertices along trajectories sharing the same endpoints.
#'
#' The algorithm proceeds as follows. First, all local minima are identified
#' using neighborhood comparison. From each local minimum, an ascending
#' trajectory is traced following the steepest gradient (with optional
#' modulation) until reaching a local maximum. All trajectory vertices are
#' assigned to both the starting min-basin and ending max-basin. For any
#' unvisited vertices, both ascending and descending trajectories are traced
#' and joined at that vertex. Finally, the filtering and merging pipeline is
#' applied, including relative value filtering, overlap-based clustering, and
#' geometric filtering.
#'
#' The trajectory-first approach naturally covers all vertices without
#' requiring a separate expansion step.
#'
#' @section Modulation Options:
#' The \code{modulation} parameter controls how the steepest direction is
#' determined. With \code{"NONE"}, standard gradient flow uses raw function
#' differences. With \code{"DENSITY"}, the flow is density-modulated to prefer
#' higher-density regions. With \code{"EDGELEN"}, the flow is edge-length
#' normalized to prevent basin jumping through long edges. With
#' \code{"DENSITY_EDGELEN"}, both modulations are combined.
#'
#' @param adj.list List of integer vectors. Each element \code{adj.list[[i]]}
#'   contains the 1-based indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors. Each element
#'   \code{weight.list[[i]]} contains the edge weights (distances) corresponding
#'   to \code{adj.list[[i]]}.
#' @param y Numeric vector of function values at each vertex.
#' @param density Optional numeric vector of density values at each vertex.
#'   Required if \code{modulation} is \code{"DENSITY"} or
#'   \code{"DENSITY_EDGELEN"}. If not provided, these modulations fall back
#'   to \code{"NONE"} or \code{"EDGELEN"} respectively.
#' @param modulation Character string specifying gradient modulation strategy.
#'   One of \code{"NONE"}, \code{"DENSITY"}, \code{"EDGELEN"}, or
#'   \code{"DENSITY_EDGELEN"}. Default is \code{"NONE"}.
#' @param edge.length.quantile.thld Numeric in (0,1]. Quantile threshold for
#'   edge length filtering during trajectory computation. Edges longer than
#'   this quantile are not traversed. Default is 0.9.
#' @param apply.relvalue.filter Logical. Whether to filter extrema by relative
#'   value. Default is \code{TRUE}.
#' @param min.rel.value.max Numeric. Minimum relative value (value/median) for
#'   retaining maxima. Default is 1.1.
#' @param max.rel.value.min Numeric. Maximum relative value for retaining
#'   minima. Default is 0.9.
#' @param apply.maxima.clustering Logical. Whether to cluster and merge maxima
#'   based on basin overlap. Default is \code{TRUE}.
#' @param apply.minima.clustering Logical. Whether to cluster and merge minima.
#'   Default is \code{TRUE}.
#' @param max.overlap.threshold Numeric in \eqn{[0,1]}. Overlap distance threshold
#'   for merging maxima basins. Default is 0.15.
#' @param min.overlap.threshold Numeric in \eqn{[0,1]}. Overlap distance threshold
#'   for merging minima basins. Default is 0.15.
#' @param apply.geometric.filter Logical. Whether to filter by geometric
#'   characteristics. Default is \code{TRUE}.
#' @param p.mean.nbrs.dist.threshold Numeric in \eqn{[0,1]}. Percentile threshold
#'   for mean neighbor distance (maxima only). Default is 0.9.
#' @param p.mean.hopk.dist.threshold Numeric in \eqn{[0,1]}. Percentile threshold
#'   for hop-k distance. Default is 0.9.
#' @param p.deg.threshold Numeric in \eqn{[0,1]}. Degree percentile threshold.
#'   Default is 0.9.
#' @param min.basin.size Integer. Minimum basin size to retain. Default is 10.
#' @param hop.k Integer. Hop distance for computing summary statistics.
#'   Default is 2.
#' @param store.trajectories Logical. Whether to store full trajectory
#'   information in the result. Default is \code{TRUE}.
#' @param max.trajectory.length Integer. Maximum trajectory length (vertices)
#'   before stopping. Prevents infinite loops in degenerate cases.
#'   Default is 10000.
#' @param verbose Logical. Whether to print progress messages. Default is
#'   \code{FALSE}.
#'
#' @return A list of class \code{"gfc.flow"} containing:
#'   \item{max.basins}{List of maximum basin structures, each containing
#'     \code{extremum.vertex}, \code{extremum.value}, \code{is.maximum},
#'     \code{vertices}, \code{hop.distances}, \code{max.hop.distance}.}
#'   \item{min.basins}{List of minimum basin structures with the same fields.}
#'   \item{max.assignment}{Integer vector where \code{max.assignment[v]} is the
#'     1-based index of the max basin containing vertex \code{v}, or \code{NA}.}
#'   \item{min.assignment}{Integer vector for min basin assignments.}
#'   \item{trajectories}{List of trajectory structures (if
#'     \code{store.trajectories = TRUE}), each containing \code{vertices},
#'     \code{start.vertex}, \code{end.vertex}, \code{starts.at.lmin},
#'     \code{ends.at.lmax}, \code{total.change}, \code{trajectory.id}.}
#'   \item{vertex.trajectory}{Integer vector mapping vertices to trajectories.}
#'   \item{n.lmin.trajectories}{Number of trajectories started from local minima.}
#'   \item{n.join.trajectories}{Number of trajectories from joining operations.}
#'   \item{n.vertices}{Total number of vertices.}
#'   \item{y.median}{Median of function values.}
#'   \item{stage.history}{Data frame tracking extrema counts through pipeline.}
#'   \item{modulation}{Modulation strategy used.}
#'
#' @seealso \code{\link{compute.gfc}} for the extrema-first approach
#'
#' @examples
#' \donttest{
#' ## Create a simple graph
#' adj.list <- list(
#'   c(2L, 3L),
#'   c(1L, 3L, 4L),
#'   c(1L, 2L, 4L),
#'   c(2L, 3L)
#' )
#' weight.list <- list(
#'   c(1.0, 1.5),
#'   c(1.0, 1.2, 1.0),
#'   c(1.5, 1.2, 1.0),
#'   c(1.0, 1.0)
#' )
#' y <- c(0.5, 1.2, 0.8, 2.0)
#'
#' ## Compute GFC using trajectory approach
#' result <- compute.gfc.flow(adj.list, weight.list, y)
#'
#' ## With edge-length modulation
#' result2 <- compute.gfc.flow(adj.list, weight.list, y,
#'                             modulation = "EDGELEN")
#' }
#'
#' @export
compute.gfc.flow <- function(
    adj.list,
    weight.list,
    y,
    density = NULL,
    modulation = c("NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN"),
    edge.length.quantile.thld = 0.9,
    apply.relvalue.filter = TRUE,
    min.rel.value.max = 1.1,
    max.rel.value.min = 0.9,
    apply.maxima.clustering = TRUE,
    apply.minima.clustering = TRUE,
    max.overlap.threshold = 0.15,
    min.overlap.threshold = 0.15,
    apply.geometric.filter = TRUE,
    p.mean.nbrs.dist.threshold = 0.9,
    p.mean.hopk.dist.threshold = 0.9,
    p.deg.threshold = 0.9,
    min.basin.size = 10L,
    hop.k = 2L,
    store.trajectories = TRUE,
    max.trajectory.length = 10000L,
    verbose = FALSE
) {
    ## -------------------------------------------------------------------------
    ## Input validation
    ## -------------------------------------------------------------------------

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
    if (length(y) != n.vertices) {
        stop("y must have the same length as adj.list")
    }

    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }

    ## Validate density if provided
    if (!is.null(density)) {
        if (!is.numeric(density)) {
            stop("density must be a numeric vector")
        }
        if (length(density) != n.vertices) {
            stop("density must have the same length as adj.list")
        }
    }

    ## Match modulation argument
    modulation <- match.arg(modulation)

    ## Validate numeric parameters
    if (edge.length.quantile.thld <= 0 || edge.length.quantile.thld > 1) {
        stop("edge.length.quantile.thld must be in (0, 1]")
    }

    ## -------------------------------------------------------------------------
    ## Convert adjacency list to 0-based indexing for C++
    ## -------------------------------------------------------------------------

    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            return(integer(0))
        }
        as.integer(x - 1L)
    })

    ## -------------------------------------------------------------------------
    ## Build parameter list
    ## -------------------------------------------------------------------------

    params <- list(
        edge_length_quantile_thld = as.double(edge.length.quantile.thld),
        apply_relvalue_filter = as.logical(apply.relvalue.filter),
        min_rel_value_max = as.double(min.rel.value.max),
        max_rel_value_min = as.double(max.rel.value.min),
        apply_maxima_clustering = as.logical(apply.maxima.clustering),
        apply_minima_clustering = as.logical(apply.minima.clustering),
        max_overlap_threshold = as.double(max.overlap.threshold),
        min_overlap_threshold = as.double(min.overlap.threshold),
        apply_geometric_filter = as.logical(apply.geometric.filter),
        p_mean_nbrs_dist_threshold = as.double(p.mean.nbrs.dist.threshold),
        p_mean_hopk_dist_threshold = as.double(p.mean.hopk.dist.threshold),
        p_deg_threshold = as.double(p.deg.threshold),
        min_basin_size = as.integer(min.basin.size),
        hop_k = as.integer(hop.k),
        modulation = as.character(modulation),
        store_trajectories = as.logical(store.trajectories),
        max_trajectory_length = as.integer(max.trajectory.length)
    )

    ## -------------------------------------------------------------------------
    ## Call C++ implementation
    ## -------------------------------------------------------------------------

    result <- .Call(
        S_compute_gfc_flow,
        adj.list.0based,
        weight.list,
        as.double(y),
        if (is.null(density)) double(0) else as.double(density),
        params,
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## -------------------------------------------------------------------------
    ## Convert names from snake_case to dot.snake for R convention
    ## -------------------------------------------------------------------------

    names(result) <- gsub("_", ".", names(result))

    ## Rename basin list elements
    result$max.basins <- lapply(result$max.basins, function(basin) {
        names(basin) <- gsub("_", ".", names(basin))
        basin
    })
    result$min.basins <- lapply(result$min.basins, function(basin) {
        names(basin) <- gsub("_", ".", names(basin))
        basin
    })

    ## Rename trajectory list elements if present
    if (store.trajectories && length(result$trajectories) > 0) {
        result$trajectories <- lapply(result$trajectories, function(traj) {
            names(traj) <- gsub("_", ".", names(traj))
            traj
        })
    }

    ## Rename stage.history columns
    if (!is.null(result$stage.history)) {
        names(result$stage.history) <- gsub("_", ".", names(result$stage.history))
    }

    ## -------------------------------------------------------------------------
    ## Add metadata and class
    ## -------------------------------------------------------------------------

    result$call <- match.call()
    class(result) <- c("gfc.flow", "list")

    return(result)
}


#' Compute GFC Flow for Multiple Functions
#'
#' @description
#' Efficiently computes trajectory-based GFC for each column of a matrix,
#' reusing graph structure across all computations. Supports parallel
#' execution via OpenMP.
#'
#' @param adj.list List of integer vectors (adjacency list, 1-based).
#' @param weight.list List of numeric vectors (edge weights).
#' @param Y Numeric matrix of function values (n.vertices x n.functions).
#' @param density Optional numeric vector of density values.
#' @param modulation Character string specifying gradient modulation.
#' @param n.cores Integer number of OpenMP threads. Default is 1 (sequential).
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments passed to \code{compute.gfc.flow()}.
#'
#' @return List of \code{gfc.flow} results, one per column of \code{Y}.
#'
#' @seealso \code{\link{compute.gfc.flow}}
#'
#' @export
compute.gfc.flow.matrix <- function(
    adj.list,
    weight.list,
    Y,
    density = NULL,
    modulation = c("NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN"),
    n.cores = 1L,
    verbose = FALSE,
    ...
) {
    ## -------------------------------------------------------------------------
    ## Input validation
    ## -------------------------------------------------------------------------

    if (!is.matrix(Y)) {
        stop("Y must be a matrix")
    }

    n.vertices <- length(adj.list)
    if (nrow(Y) != n.vertices) {
        stop("Number of rows in Y must equal length of adj.list")
    }

    modulation <- match.arg(modulation)

    ## -------------------------------------------------------------------------
    ## Convert adjacency list to 0-based
    ## -------------------------------------------------------------------------

    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            return(integer(0))
        }
        as.integer(x - 1L)
    })

    ## -------------------------------------------------------------------------
    ## Build parameter list from ...
    ## -------------------------------------------------------------------------

    dots <- list(...)

    params <- list(
        edge_length_quantile_thld = as.double(dots$edge.length.quantile.thld %||% 0.9),
        apply_relvalue_filter = as.logical(dots$apply.relvalue.filter %||% TRUE),
        min_rel_value_max = as.double(dots$min.rel.value.max %||% 1.1),
        max_rel_value_min = as.double(dots$max.rel.value.min %||% 0.9),
        apply_maxima_clustering = as.logical(dots$apply.maxima.clustering %||% TRUE),
        apply_minima_clustering = as.logical(dots$apply.minima.clustering %||% TRUE),
        max_overlap_threshold = as.double(dots$max.overlap.threshold %||% 0.15),
        min_overlap_threshold = as.double(dots$min.overlap.threshold %||% 0.15),
        apply_geometric_filter = as.logical(dots$apply.geometric.filter %||% TRUE),
        p_mean_nbrs_dist_threshold = as.double(dots$p.mean.nbrs.dist.threshold %||% 0.9),
        p_mean_hopk_dist_threshold = as.double(dots$p.mean.hopk.dist.threshold %||% 0.9),
        p_deg_threshold = as.double(dots$p.deg.threshold %||% 0.9),
        min_basin_size = as.integer(dots$min.basin.size %||% 10L),
        hop_k = as.integer(dots$hop.k %||% 2L),
        modulation = as.character(modulation),
        store_trajectories = as.logical(dots$store.trajectories %||% TRUE),
        max_trajectory_length = as.integer(dots$max.trajectory.length %||% 10000L)
    )

    ## -------------------------------------------------------------------------
    ## Call C++ implementation
    ## -------------------------------------------------------------------------

    results <- .Call(
        S_compute_gfc_flow_matrix,
        adj.list.0based,
        weight.list,
        Y,
        if (is.null(density)) double(0) else as.double(density),
        params,
        as.integer(n.cores),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## -------------------------------------------------------------------------
    ## Post-process each result
    ## -------------------------------------------------------------------------

    results <- lapply(results, function(result) {
        names(result) <- gsub("_", ".", names(result))

        result$max.basins <- lapply(result$max.basins, function(basin) {
            names(basin) <- gsub("_", ".", names(basin))
            basin
        })
        result$min.basins <- lapply(result$min.basins, function(basin) {
            names(basin) <- gsub("_", ".", names(basin))
            basin
        })

        if (length(result$trajectories) > 0) {
            result$trajectories <- lapply(result$trajectories, function(traj) {
                names(traj) <- gsub("_", ".", names(traj))
                traj
            })
        }

        if (!is.null(result$stage.history)) {
            names(result$stage.history) <- gsub("_", ".", names(result$stage.history))
        }

        class(result) <- c("gfc.flow", "list")
        result
    })

    return(results)
}


#' Print Method for gfc.flow Objects
#'
#' @param x A \code{gfc.flow} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.gfc.flow <- function(x, ...) {
    cat("Gradient Flow Complex (Trajectory-First)\n")
    cat("=========================================\n")
    cat(sprintf("Vertices: %d\n", x$n.vertices))
    cat(sprintf("Modulation: %s\n", x$modulation))
    cat(sprintf("Maximum basins: %d\n", length(x$max.basins)))
    cat(sprintf("Minimum basins: %d\n", length(x$min.basins)))
    cat(sprintf("Trajectories: %d (%d from lmin, %d from joins)\n",
                x$n.lmin.trajectories + x$n.join.trajectories,
                x$n.lmin.trajectories,
                x$n.join.trajectories))
    cat(sprintf("y median: %.4f\n", x$y.median))

    if (!is.null(x$stage.history) && nrow(x$stage.history) > 0) {
        cat("\nRefinement pipeline:\n")
        print(x$stage.history, row.names = FALSE)
    }

    invisible(x)
}


#' Summary Method for gfc.flow Objects
#'
#' @param object A \code{gfc.flow} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.gfc.flow <- function(object, ...) {
    cat("Gradient Flow Complex Summary\n")
    cat("=============================\n\n")

    cat("Graph:\n")
    cat(sprintf("  Vertices: %d\n", object$n.vertices))
    cat(sprintf("  y range: [%.4f, %.4f], median: %.4f\n",
                min(sapply(object$min.basins, function(b) b$extremum.value)),
                max(sapply(object$max.basins, function(b) b$extremum.value)),
                object$y.median))

    cat("\nTrajectory computation:\n")
    cat(sprintf("  Modulation: %s\n", object$modulation))
    cat(sprintf("  Total trajectories: %d\n",
                object$n.lmin.trajectories + object$n.join.trajectories))
    cat(sprintf("    From local minima: %d\n", object$n.lmin.trajectories))
    cat(sprintf("    From joining: %d\n", object$n.join.trajectories))

    cat("\nBasins after refinement:\n")
    cat(sprintf("  Maximum basins: %d\n", length(object$max.basins)))
    cat(sprintf("  Minimum basins: %d\n", length(object$min.basins)))

    if (length(object$max.basins) > 0) {
        max.sizes <- sapply(object$max.basins, function(b) length(b$vertices))
        cat(sprintf("  Max basin sizes: min=%d, median=%d, max=%d\n",
                    min(max.sizes), median(max.sizes), max(max.sizes)))
    }

    if (length(object$min.basins) > 0) {
        min.sizes <- sapply(object$min.basins, function(b) length(b$vertices))
        cat(sprintf("  Min basin sizes: min=%d, median=%d, max=%d\n",
                    min(min.sizes), median(min.sizes), max(min.sizes)))
    }

    ## Coverage statistics
    max.covered <- sum(!is.na(object$max.assignment))
    min.covered <- sum(!is.na(object$min.assignment))
    cat(sprintf("\nCoverage:\n"))
    cat(sprintf("  Max basins: %d/%d vertices (%.1f%%)\n",
                max.covered, object$n.vertices, 100 * max.covered / object$n.vertices))
    cat(sprintf("  Min basins: %d/%d vertices (%.1f%%)\n",
                min.covered, object$n.vertices, 100 * min.covered / object$n.vertices))

    invisible(object)
}


## Helper: null-coalescing operator
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

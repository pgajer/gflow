#' Compute Gradient Flow Complex Using Trajectory-First Approach
#'
#' @description
#' Computes basins of attraction by tracing gradient flow trajectories rather
#' than growing basins outward from extrema. This approach treats gradient flow
#' lines as the fundamental geometric primitives, with basins emerging as the
#' collection of vertices along trajectories sharing the same endpoints.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Identify all local minima using neighborhood comparison
#'   \item From each local minimum, trace an ascending trajectory following
#'         the selected modulation strategy until reaching a local maximum
#'   \item Assign all trajectory vertices to both the starting min-basin
#'         and ending max-basin
#'   \item For any unvisited vertices, trace both ascending and descending
#'         trajectories, joining them at that vertex
#'   \item Apply the filtering/merging pipeline (relative value filtering,
#'         overlap-based clustering, geometric filtering)
#' }
#'
#' Retained extrema are numbered by descending basin size.
#'
#' @section Modulation Options:
#' The \code{modulation} parameter controls how the next vertex is selected
#' at each step of the trajectory:
#'
#' \describe{
#'   \item{\code{"NONE"}}{Standard steepest ascent/descent. Selects the neighbor
#'     maximizing the function difference \eqn{|y(u) - y(v)|}. This is the
#'     classical gradient flow approach but may cause basin jumping near
#'     separatrices when adjacent basins have very different gradient magnitudes.}
#'
#'   \item{\code{"DENSITY"}}{Density-modulated flow. Multiplies the gradient by
#'     vertex density to prefer higher-density regions: \eqn{\rho(u) \cdot \Delta y}.}
#'
#'   \item{\code{"EDGELEN"}}{Edge-length modulated flow. Multiplies the gradient
#'     by a KDE-based edge length weight that down-weights atypically long edges:
#'     \eqn{w(d(v,u)) \cdot \Delta y}. This prevents basin jumping through outlier
#'     edges but still uses a multiplicative score that can be dominated by
#'     extreme gradient contrasts.}
#'
#'   \item{\code{"DENSITY_EDGELEN"}}{Combined density and edge-length modulation.}
#'
#'   \item{\code{"CLOSEST"}}{Lexicographic closest ascending neighbor rule. This
#'     implements a two-level selection:
#'     \enumerate{
#'       \item Filter to ascending neighbors: \eqn{A(v) = \{u \in N(v) : y(u) > y(v)\}}
#'       \item Among \eqn{A(v)}, select the closest: \eqn{u^* = \arg\min_{u \in A(v)} d(v,u)}
#'     }
#'     This approach minimizes basin-jumping errors by taking the smallest step
#'     that makes progress. Unlike multiplicative scores, the lexicographic rule
#'     cannot be overridden by extreme gradient contrasts between adjacent basins.
#'     Requires no tuning parameters beyond the standard edge length threshold.
#'
#'     The theoretical justification is that any monotonically ascending discrete
#'     path will reach the correct local maximum provided consecutive steps remain
#'     within the same basin. By selecting the closest ascending neighbor, we
#'     maximize the probability of staying within the basin of the current vertex.}
#' }
#'
#' @param adj.list List of integer vectors. Each element \code{adj.list[[i]]}
#'   contains the 1-based indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors. Each element
#'   \code{weight.list[[i]]} contains the edge weights (distances) corresponding
#'   to \code{adj.list[[i]]}.
#' @param y Numeric vector of function values at each vertex.
#' @param density Optional numeric vector of density values at each vertex.
#'   Required if \code{modulation} is \code{"DENSITY"} or
#'   \code{"DENSITY_EDGELEN"}. Ignored for \code{"CLOSEST"}.
#' @param modulation Character string specifying gradient modulation strategy.
#'   One of \code{"NONE"}, \code{"DENSITY"}, \code{"EDGELEN"},
#'   \code{"DENSITY_EDGELEN"}, or \code{"CLOSEST"}. Default is \code{"CLOSEST"}.
#'   See Details for descriptions of each option.
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
#' @param min.n.trajectories Integer. Minimum number of basin trajectories to retain the basin. Default is 0.
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
#' @return A list of class \code{"gfc.flow"} containing basin and trajectory
#'   information. See the original function documentation for full details.
#'
#' @examples
#' \dontrun{
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
#' ## Standard steepest ascent (may have basin jumping issues)
#' result.steep <- compute.gfc.flow(adj.list, weight.list, y,
#'                                  modulation = "NONE")
#'
#' ## Edge-length modulated (reduces but doesn't eliminate jumping)
#' result.edgelen <- compute.gfc.flow(adj.list, weight.list, y,
#'                                    modulation = "EDGELEN")
#'
#' ## Closest ascending neighbor (lexicographic, no tuning needed)
#' result.closest <- compute.gfc.flow(adj.list, weight.list, y,
#'                                    modulation = "CLOSEST")
#' }
#'
#' @seealso \code{\link{compute.gfc}} for the extrema-first approach
#'
#' @export
compute.gfc.flow <- function(
    adj.list,
    weight.list,
    y,
    density = NULL,
    modulation = c("CLOSEST", "NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN"),
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
    min.n.trajectories = 0L,
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

    if (!is.numeric(min.n.trajectories) || length(min.n.trajectories) != 1L) {
        stop("min.n.trajectories must be a single integer-like value")
    }
    if (min.n.trajectories < 0) {
        stop("min.n.trajectories must be >= 0")
    }

    ## -------------------------------------------------------
    ## Breaking ties (if any)
    ## -------------------------------------------------------

    y <- break.ties(y,
                    noise.scale = 1e-15,
                    min.abs.noise = 1e-16,
                    preserve.bounds = TRUE,
                    seed = NULL,
                    verbose = FALSE)

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
        min_n_trajectories = as.integer(min.n.trajectories),
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
    ## Post-process result
    ## -------------------------------------------------------------------------

    result <- .postprocess.gfc.flow(result)

    result$y <- y

    result$call <- match.call()
    class(result) <- c("gfc.flow", "list")

    return(result)
}

#' Post-process GFC Flow Result
#'
#' Builds combined summary tables from C++ output. Since C++ now returns
#' dot.snake names directly, no name conversion is needed.
#'
#' @param result Raw result from .Call(1S_compute_gfc_flow, ...)
#' @return Processed gfc.flow object
#'
#' @keywords internal
.postprocess.gfc.flow <- function(result) {

    ## -------------------------------------------------------------------------
    ## Build combined summary table (all extrema) with sorting
    ## -------------------------------------------------------------------------

    result$summary.all <- .build.summary.all(result)

    ## -------------------------------------------------------------------------
    ## Build retained-only summary (inherits sorting from summary.all)
    ## -------------------------------------------------------------------------

    result$summary <- .build.summary.retained(result)

    ## -------------------------------------------------------------------------
    ## Set class
    ## -------------------------------------------------------------------------

    class(result) <- c("gfc.flow", "list")

    return(result)
}


#' Build Summary Table for All Extrema
#'
#' Combines min and max summaries from C++, adds type column, reorders columns,
#' and sorts by: (1) is.spurious, (2) type, (3) label number.
#'
#' @param result Raw gfc.flow result from C++
#' @return Sorted data frame with all extrema (retained + spurious)
#'
#' @keywords internal
.build.summary.all <- function(result) {

    summary.table <- NULL

    ## Add minima
    if (!is.null(result$min.summaries.all) && nrow(result$min.summaries.all) > 0) {
        min.summary <- result$min.summaries.all
        min.summary$type <- "min"
        summary.table <- min.summary
    }

    ## Add maxima
    if (!is.null(result$max.summaries.all) && nrow(result$max.summaries.all) > 0) {
        max.summary <- result$max.summaries.all
        max.summary$type <- "max"
        if (is.null(summary.table)) {
            summary.table <- max.summary
        } else {
            summary.table <- rbind(summary.table, max.summary)
        }
    }

    if (is.null(summary.table) || nrow(summary.table) == 0) {
        return(NULL)
    }

    ## -------------------------------------------------------------------------
    ## Reorder columns for clarity
    ## -------------------------------------------------------------------------

    col.order <- c("label", "vertex", "value", "rel.value", "type",
                   "is.spurious", "filter.stage", "merged.into",
                   "basin.size", "hop.index", "p.mean.nbrs.dist",
                   "p.mean.hopk.dist", "degree", "deg.percentile")
    col.order <- col.order[col.order %in% names(summary.table)]
    summary.table <- summary.table[, col.order, drop = FALSE]

    ## -------------------------------------------------------------------------
    ## Sort by: is.spurious (FALSE first), type (min first), label number
    ## -------------------------------------------------------------------------

    summary.table <- .sort.summary.df(summary.table)

    return(summary.table)
}


#' Build Summary Table for Retained Extrema Only
#'
#' Filters summary.all to retained extrema and removes spurious-specific columns.
#' Sorting is inherited from summary.all.
#'
#' @param result Processed gfc.flow result (must have summary.all already built)
#' @return Data frame with retained extrema only
#'
#' @keywords internal
.build.summary.retained <- function(result) {

    summary.all <- result$summary.all

    if (is.null(summary.all) || nrow(summary.all) == 0) {
        return(NULL)
    }

    ## Filter to retained only (already sorted from summary.all)
    summary.retained <- summary.all[!summary.all$is.spurious, , drop = FALSE]

    if (nrow(summary.retained) == 0) {
        return(NULL)
    }

    ## Remove spurious-specific columns for cleaner output
    summary.retained$is.spurious <- NULL
    summary.retained$filter.stage <- NULL
    summary.retained$merged.into <- NULL

    rownames(summary.retained) <- NULL

    return(summary.retained)
}


#' Sort Summary Data Frame
#'
#' Sorts a summary data frame by:
#'   1. is.spurious (FALSE first, TRUE second)
#'   2. type (min first, max second)
#'   3. numeric part of label
#'
#' @param df Data frame with label, is.spurious, type columns
#' @return Sorted data frame with reset row names
#'
#' @keywords internal
.sort.summary.df <- function(df) {

    if (is.null(df) || nrow(df) == 0) {
        return(df)
    }

    ## Extract numeric suffix from label: m1 -> 1, M12 -> 12, sm3 -> 3, sM7 -> 7
    label.number <- as.integer(gsub("^[smSM]+", "", df$label))

    ## Create sort keys
    sort.key.spurious <- as.integer(df$is.spurious)
    sort.key.type <- ifelse(df$type == "min", 0L, 1L)

    ## Get sort order
    ord <- order(sort.key.spurious, sort.key.type, label.number)

    ## Return sorted data frame with reset row names
    sorted.df <- df[ord, , drop = FALSE]
    rownames(sorted.df) <- NULL

    return(sorted.df)
}

#' Print Method for gfc.flow Objects
#'
#' @param x A gfc.flow object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gfc.flow <- function(x, ...) {
    cat("Gradient Flow Complex (trajectory-first approach)\n")
    cat("==================================================\n")
    cat(sprintf("Vertices: %d\n", x$n.vertices))
    cat(sprintf("Median y: %.4f\n", x$y.median))
    cat(sprintf("Modulation: %s\n", x$modulation))
    cat("\n")

    cat("EXTREMA:\n")
    cat(sprintf("  Retained: %d minima (m1-%s), %d maxima (M1-%s)\n",
                x$n.min.retained,
                if (x$n.min.retained > 0) paste0("m", x$n.min.retained) else "none",
                x$n.max.retained,
                if (x$n.max.retained > 0) paste0("M", x$n.max.retained) else "none"))
    cat(sprintf("  Spurious: %d minima (sm1-%s), %d maxima (sM1-%s)\n",
                x$n.min.spurious,
                if (x$n.min.spurious > 0) paste0("sm", x$n.min.spurious) else "none",
                x$n.max.spurious,
                if (x$n.max.spurious > 0) paste0("sM", x$n.max.spurious) else "none"))
    cat(sprintf("  Total: %d minima, %d maxima\n",
                x$n.min.all, x$n.max.all))
    cat("\n")

    cat("TRAJECTORIES:\n")
    cat(sprintf("  Total: %d (%d from minima, %d from joins)\n",
                length(x$trajectories),
                x$n.lmin.trajectories,
                x$n.join.trajectories))

    ## Count trajectories by endpoint status
    if (length(x$trajectories) > 0) {
        n.both.retained <- sum(sapply(x$trajectories, function(t) {
            !t$start.is.spurious && !t$end.is.spurious
        }))
        n.start.spurious <- sum(sapply(x$trajectories, function(t) {
            t$start.is.spurious && !t$end.is.spurious
        }))
        n.end.spurious <- sum(sapply(x$trajectories, function(t) {
            !t$start.is.spurious && t$end.is.spurious
        }))
        n.both.spurious <- sum(sapply(x$trajectories, function(t) {
            t$start.is.spurious && t$end.is.spurious
        }))

        cat(sprintf("  Both endpoints retained: %d\n", n.both.retained))
        cat(sprintf("  Start spurious only: %d\n", n.start.spurious))
        cat(sprintf("  End spurious only: %d\n", n.end.spurious))
        cat(sprintf("  Both endpoints spurious: %d\n", n.both.spurious))
    }
    cat("\n")

    ## -------------------------------------------------------------------------
    ## Coverage information
    ## -------------------------------------------------------------------------

    uncov <- uncovered.vertices.gfc.flow(x)

    cat("COVERAGE:\n")
    cat(sprintf("  Covered: %d/%d vertices (%.1f%%)\n",
                uncov$n.total - uncov$n.uncovered,
                uncov$n.total,
                uncov$coverage * 100))

    if (uncov$n.uncovered > 0) {
        cat(sprintf("  Uncovered: %d (isolated by edge length threshold)\n",
                    uncov$n.uncovered))
        if (uncov$n.local.minima > 0 || uncov$n.local.maxima > 0) {
            cat(sprintf("    Including: %d local minima, %d local maxima without valid basins\n",
                        uncov$n.local.minima, uncov$n.local.maxima))
        }
    }

    cat("\n")
    invisible(x)
}

#' Summary Method for gfc.flow Objects
#'
#' @param object A gfc.flow object
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.gfc.flow <- function(object, ...) {
    print(object)
    cat("SUMMARY TABLE (all extrema):\n")
    if (!is.null(object$summary.all) && nrow(object$summary.all) > 0) {
        print(object$summary.all)
    } else {
        cat("  (none)\n")
    }
    invisible(object)
}


#' Get Basin by Label
#'
#' Retrieves a basin structure by its label (m1, M2, sm3, sM4, etc.)
#'
#' @param gfc.flow A gfc.flow object
#' @param label Character label of the basin
#'
#' @return The basin structure, or NULL if not found
#'
#' @examples
#' \donttest{
#' ## Get retained maximum M1
#' basin.M1 <- get.basin(gfc.flow.res, "M1")
#'
#' ## Get spurious minimum sm2
#' basin.sm2 <- get.basin(gfc.flow.res, "sm2")
#' }
#'
#' @export
get.basin <- function(gfc.flow, label) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    ## Determine if min or max from label prefix
    is.min <- grepl("^s?m[0-9]", label)

    basins <- if (is.min) gfc.flow$min.basins.all else gfc.flow$max.basins.all

    if (is.null(basins)) return(NULL)

    for (basin in basins) {
        if (!is.null(basin$label) && basin$label == label) {
            return(basin)
        }
    }

    return(NULL)
}


#' List Retained Extrema
#'
#' @param gfc.flow A gfc.flow object
#' @param type "min", "max", or "all"
#'
#' @return Data frame of retained extrema
#'
#' @export
list.retained <- function(gfc.flow, type = c("all", "min", "max")) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    type <- match.arg(type)

    if (is.null(gfc.flow$summary.all)) return(NULL)

    retained <- gfc.flow$summary.all[!gfc.flow$summary.all$is.spurious, , drop = FALSE]

    if (type == "min") {
        retained <- retained[retained$type == "min", , drop = FALSE]
    } else if (type == "max") {
        retained <- retained[retained$type == "max", , drop = FALSE]
    }

    rownames(retained) <- NULL
    return(retained)
}


#' List Spurious Extrema
#'
#' @param gfc.flow A gfc.flow object
#' @param type "min", "max", or "all"
#' @param filter.stage Optional: filter by stage ("relvalue", "cluster.merge", "geometric")
#'
#' @return Data frame of spurious extrema
#'
#' @export
list.spurious <- function(gfc.flow, type = c("all", "min", "max"),
                          filter.stage = NULL) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    type <- match.arg(type)

    if (is.null(gfc.flow$summary.all)) return(NULL)

    spurious <- gfc.flow$summary.all[gfc.flow$summary.all$is.spurious, , drop = FALSE]

    if (type == "min") {
        spurious <- spurious[spurious$type == "min", , drop = FALSE]
    } else if (type == "max") {
        spurious <- spurious[spurious$type == "max", , drop = FALSE]
    }

    if (!is.null(filter.stage)) {
        spurious <- spurious[spurious$filter.stage == filter.stage, , drop = FALSE]
    }

    rownames(spurious) <- NULL
    return(spurious)
}


#' Get Vertices for Harmonic Repair
#'
#' Returns the vertices of a spurious basin that need harmonic repair,
#' along with the boundary vertices and their values.
#'
#' @param gfc.flow A gfc.flow object
#' @param label Label of a spurious extremum (e.g., "sm1", "sM2")
#' @param y Original function values (used to set boundary values)
#' @param adj.list Adjacency list (for finding boundary)
#'
#' @return A list with:
#'   \item{interior.vertices}{Vertices to be repaired}
#'   \item{boundary.vertices}{Vertices on the boundary}
#'   \item{boundary.values}{Current y values at boundary}
#'   \item{extremum.vertex}{The spurious extremum vertex}
#'
#' @export
get.harmonic.repair.vertices <- function(gfc.flow, label, y, adj.list) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    basin <- get.basin(gfc.flow, label)
    if (is.null(basin)) {
        stop(sprintf("Basin '%s' not found", label))
    }

    extremum.vertex <- as.integer(basin$extremum.vertex)
    repair.vertices <- unique(as.integer(basin$vertices))

    compute.boundary.interior <- function(vertices) {
        v.set <- as.integer(vertices)
        boundary.vertices <- integer(0)

        for (v in v.set) {
            nbrs <- adj.list[[v]]
            if (length(nbrs) == 0) next
            if (any(!(nbrs %in% v.set))) {
                boundary.vertices <- c(boundary.vertices, v)
            }
        }

        boundary.vertices <- unique(as.integer(boundary.vertices))
        interior.vertices <- setdiff(v.set, boundary.vertices)

        list(boundary.vertices = boundary.vertices,
             interior.vertices = as.integer(interior.vertices))
    }

    bi <- compute.boundary.interior(repair.vertices)

    ## Fallback: ensure extremum.vertex is interior (so it can be modified)
    if (!extremum.vertex %in% bi$interior.vertices) {
        repair.vertices <- unique(c(repair.vertices, as.integer(adj.list[[extremum.vertex]])))
        bi <- compute.boundary.interior(repair.vertices)
    }

    boundary.values <- y[bi$boundary.vertices]
    names(boundary.values) <- as.character(bi$boundary.vertices)

    list(
        interior.vertices = bi$interior.vertices,
        boundary.vertices = bi$boundary.vertices,
        boundary.values = boundary.values,
        extremum.vertex = extremum.vertex,
        label = label,
        repair.vertices = repair.vertices
    )
}

#' Get Trajectories by Endpoint Status
#'
#' @param gfc.flow A gfc.flow object
#' @param start.retained Logical: require start endpoint retained? (NA = any)
#' @param end.retained Logical: require end endpoint retained? (NA = any)
#'
#' @return List of matching trajectories
#'
#' @export
get.trajectories.by.status <- function(gfc.flow,
                                       start.retained = NA,
                                       end.retained = NA) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        return(list())
    }

    result <- list()
    for (traj in gfc.flow$trajectories) {
        match <- TRUE

        if (!is.na(start.retained)) {
            if (start.retained && traj$start.is.spurious) match <- FALSE
            if (!start.retained && !traj$start.is.spurious) match <- FALSE
        }

        if (!is.na(end.retained)) {
            if (end.retained && traj$end.is.spurious) match <- FALSE
            if (!end.retained && !traj$end.is.spurious) match <- FALSE
        }

        if (match) {
            result[[length(result) + 1]] <- traj
        }
    }

    return(result)
}


#' Count Cell Membership Statistics
#'
#' @param gfc.flow A gfc.flow object
#'
#' @return List with membership statistics
#'
#' @export
count.cell.memberships <- function(gfc.flow) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    n <- gfc.flow$n.vertices

    ## Count memberships using ALL basins
    n.max.all <- sapply(gfc.flow$max.membership.all, length)
    n.min.all <- sapply(gfc.flow$min.membership.all, length)

    n.cells.all <- ifelse(
        n.max.all > 0 & n.min.all > 0,
        n.max.all * n.min.all,
        0
    )

    ## Count using RETAINED basins only
    n.max.ret <- sapply(gfc.flow$max.membership.retained, length)
    n.min.ret <- sapply(gfc.flow$min.membership.retained, length)

    n.cells.retained <- ifelse(
        n.max.ret > 0 & n.min.ret > 0,
        n.max.ret * n.min.ret,
        0
    )

    list(
        ## Using ALL basins (should be 100% coverage)
        n.with.any.membership = sum(n.cells.all > 0),
        n.no.any.membership = sum(n.cells.all == 0),
        n.multi.membership.all = sum(n.cells.all > 1),
        max.memberships.all = max(n.cells.all),

        ## Using RETAINED basins only
        n.with.retained.membership = sum(n.cells.retained > 0),
        n.no.retained.membership = sum(n.cells.retained == 0),
        n.multi.membership.retained = sum(n.cells.retained > 1),
        max.memberships.retained = max(n.cells.retained)
    )
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
#' @param ... Additional arguments passed to internal C++ function.
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
    max.trajectory.length = 10000L
) {
    ## Input validation
    if (!is.matrix(Y)) {
        stop("Y must be a matrix")
    }

    n.vertices <- length(adj.list)
    if (nrow(Y) != n.vertices) {
        stop("Number of rows in Y must equal length of adj.list")
    }

    modulation <- match.arg(modulation)

    ## Convert adjacency list to 0-based indexing
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) return(integer(0))
        as.integer(x - 1L)
    })

    ## Build parameters
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

    ## Call C++ implementation
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

    ## Post-process each result
    results <- lapply(results, .postprocess.gfc.flow.v2)

    return(results)
}


#' Extract Cell Trajectories from GFC Flow Result
#'
#' Extracts gradient flow trajectories for a specific cell (min-max pair)
#' from a gfc.flow result object. A cell is defined by a pair of extrema
#' (one minimum, one maximum) that are connected by gradient flow trajectories.
#'
#' @param gfc.flow A gfc.flow object from \code{compute.gfc.flow()} with
#'   trajectories computed (i.e., \code{store.trajectories = TRUE}).
#' @param min.id Minimum extremum identifier: either a label (e.g., "m4", "sm1")
#'   or a vertex index (integer, 1-based).
#' @param max.id Maximum extremum identifier: either a label (e.g., "M1", "sM2")
#'   or a vertex index (integer, 1-based).
#' @param map Optional integer vector mapping graph indices to subgraph
#'   vertices. If provided, trajectory vertices are converted to subgraph
#'   indices.
#'
#' @return A list of class \code{"gfc_cell_trajectories"} containing:
#'   \item{min.vertex}{Minimum vertex index (in output coordinate system).}
#'   \item{max.vertex}{Maximum vertex index (in output coordinate system).}
#'   \item{min.label}{Minimum label (e.g., "m4", "sm1").}
#'   \item{max.label}{Maximum label (e.g., "M1", "sM2").}
#'   \item{min.value}{Function value at minimum.}
#'   \item{max.value}{Function value at maximum.}
#'   \item{min.is.spurious}{TRUE if minimum is spurious.}
#'   \item{max.is.spurious}{TRUE if maximum is spurious.}
#'   \item{trajectories}{List of trajectory vertex vectors.}
#'   \item{n.trajectories}{Number of trajectories.}
#'   \item{mapped}{Logical; TRUE if vertices were mapped to subgraph indices.}
#'
#' @export
cell.trajectories.gfc.flow <- function(gfc.flow,
                                       min.id,
                                       max.id,
                                       map = NULL) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object from compute.gfc.flow()")
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        stop("gfc.flow does not contain trajectory information. ",
             "Use store.trajectories = TRUE in compute.gfc.flow()")
    }

    ## Use summary.all for v2 structure (includes spurious extrema)
    summary.df <- gfc.flow$summary.all %||% gfc.flow$summary
    if (is.null(summary.df) || nrow(summary.df) == 0) {
        stop("gfc.flow does not contain a summary table")
    }

    ## =========================================================================
    ## Resolve min.id to vertex and label
    ## =========================================================================

    if (is.character(min.id)) {
        min.label <- min.id
        idx <- match(min.label, summary.df$label)
        if (is.na(idx)) {
            stop(sprintf("Minimum label '%s' not found in summary", min.label))
        }
        min.vertex <- summary.df$vertex[idx]
        min.value <- summary.df$value[idx]
        min.type <- summary.df$type[idx]
        min.is.spurious <- if ("is.spurious" %in% names(summary.df)) {
            summary.df$is.spurious[idx]
        } else {
            FALSE
        }
    } else if (is.numeric(min.id)) {
        min.vertex <- as.integer(min.id)
        idx <- match(min.vertex, summary.df$vertex)
        if (is.na(idx)) {
            stop(sprintf("Minimum vertex %d not found in summary", min.vertex))
        }
        min.label <- summary.df$label[idx]
        min.value <- summary.df$value[idx]
        min.type <- summary.df$type[idx]
        min.is.spurious <- if ("is.spurious" %in% names(summary.df)) {
            summary.df$is.spurious[idx]
        } else {
            FALSE
        }
    } else {
        stop("min.id must be a character label or numeric vertex index")
    }

    if (min.type != "min") {
        stop(sprintf("'%s' (vertex %d) is not a minimum, it is a %s",
                     min.label, min.vertex, min.type))
    }

    ## =========================================================================
    ## Resolve max.id to vertex and label
    ## =========================================================================

    if (is.character(max.id)) {
        max.label <- max.id
        idx <- match(max.label, summary.df$label)
        if (is.na(idx)) {
            stop(sprintf("Maximum label '%s' not found in summary", max.label))
        }
        max.vertex <- summary.df$vertex[idx]
        max.value <- summary.df$value[idx]
        max.type <- summary.df$type[idx]
        max.is.spurious <- if ("is.spurious" %in% names(summary.df)) {
            summary.df$is.spurious[idx]
        } else {
            FALSE
        }
    } else if (is.numeric(max.id)) {
        max.vertex <- as.integer(max.id)
        idx <- match(max.vertex, summary.df$vertex)
        if (is.na(idx)) {
            stop(sprintf("Maximum vertex %d not found in summary", max.vertex))
        }
        max.label <- summary.df$label[idx]
        max.value <- summary.df$value[idx]
        max.type <- summary.df$type[idx]
        max.is.spurious <- if ("is.spurious" %in% names(summary.df)) {
            summary.df$is.spurious[idx]
        } else {
            FALSE
        }
    } else {
        stop("max.id must be a character label or numeric vertex index")
    }

    if (max.type != "max") {
        stop(sprintf("'%s' (vertex %d) is not a maximum, it is a %s",
                     max.label, max.vertex, max.type))
    }

    ## =========================================================================
    ## Extract trajectories for this cell
    ## =========================================================================

    trajectories <- list()

    for (traj in gfc.flow$trajectories) {
        if (traj$start.vertex == min.vertex && traj$end.vertex == max.vertex) {
            trajectories[[length(trajectories) + 1]] <- traj$vertices
        }
    }

    if (length(trajectories) == 0) {
        warning(sprintf("No trajectories found for cell %s -> %s",
                        min.label, max.label))
    }

    ## =========================================================================
    ## Apply vertex mapping if provided
    ## =========================================================================

    original.min.vertex <- min.vertex
    original.max.vertex <- max.vertex
    mapped <- FALSE

    if (!is.null(map)) {
        mapped <- TRUE
        map <- as.integer(map)

        min.vertex.mapped <- match(min.vertex, map)
        max.vertex.mapped <- match(max.vertex, map)

        if (is.na(min.vertex.mapped)) {
            warning(sprintf("Minimum vertex %d not found in map", min.vertex))
        }
        if (is.na(max.vertex.mapped)) {
            warning(sprintf("Maximum vertex %d not found in map", max.vertex))
        }

        min.vertex <- min.vertex.mapped
        max.vertex <- max.vertex.mapped

        n.unmapped <- 0
        trajectories <- lapply(trajectories, function(path) {
            mapped.path <- match(path, map)
            n.na <- sum(is.na(mapped.path))
            if (n.na > 0) {
                n.unmapped <<- n.unmapped + n.na
            }
            return(mapped.path)
        })

        if (n.unmapped > 0) {
            warning(sprintf("%d trajectory vertices not found in map (returned as NA)",
                            n.unmapped))
        }
    }

    ## =========================================================================
    ## Build result
    ## =========================================================================

    result <- list(
        min.vertex = min.vertex,
        max.vertex = max.vertex,
        min.label = min.label,
        max.label = max.label,
        min.value = min.value,
        max.value = max.value,
        min.is.spurious = min.is.spurious,
        max.is.spurious = max.is.spurious,
        trajectories = trajectories,
        n.trajectories = length(trajectories),
        mapped = mapped,
        original.min.vertex = original.min.vertex,
        original.max.vertex = original.max.vertex
    )

    class(result) <- "gfc_cell_trajectories"

    return(result)
}

#' Print Method for gfc_cell_trajectories Objects
#'
#' @description
#' Prints a compact summary of a \code{gfc_cell_trajectories} object, including
#' the cell endpoints (min/max), function values, number of trajectories, and
#' trajectory length diagnostics. Optionally prints the first few trajectories
#' with truncation for readability.
#'
#' @param x A \code{gfc_cell_trajectories} object.
#' @param max.print Integer scalar. Maximum number of trajectories to print in
#'   detail. Default is 5. Use 0 to suppress printing trajectories.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.gfc_cell_trajectories <- function(x, max.print = 5L, ...) {

    ## ------------------------------------------------------------------------
    ## Input validation (lightweight; do not over-fail in print methods)
    ## ------------------------------------------------------------------------

    if (!is.list(x)) {
        stop("print.gfc_cell_trajectories: x must be a list-like object")
    }

    if (!is.numeric(max.print) || length(max.print) != 1 || is.na(max.print) || max.print < 0) {
        stop("print.gfc_cell_trajectories: max.print must be a single non-negative number")
    }
    max.print <- as.integer(max.print)

    ## ------------------------------------------------------------------------
    ## Header
    ## ------------------------------------------------------------------------

    cat("GFC Cell Trajectories\n")
    cat("=====================\n")

    ## ------------------------------------------------------------------------
    ## Endpoint reporting (prefer original vertices if mapped)
    ## ------------------------------------------------------------------------

    mapped <- isTRUE(x$mapped)

    min.vertex.orig <- if (mapped) x$original.min.vertex else x$min.vertex
    max.vertex.orig <- if (mapped) x$original.max.vertex else x$max.vertex

    min.label <- if (!is.null(x$min.label)) as.character(x$min.label) else "min"
    max.label <- if (!is.null(x$max.label)) as.character(x$max.label) else "max"

    min.spurious <- isTRUE(x$min.is.spurious)
    max.spurious <- isTRUE(x$max.is.spurious)

    min.status <- if (min.spurious) " [SPURIOUS]" else ""
    max.status <- if (max.spurious) " [SPURIOUS]" else ""

    cat(sprintf(
        "Cell: %s%s (vertex %s) <-> %s%s (vertex %s)\n",
        min.label, min.status,
        ifelse(is.null(min.vertex.orig) || is.na(min.vertex.orig), "NA", as.character(min.vertex.orig)),
        max.label, max.status,
        ifelse(is.null(max.vertex.orig) || is.na(max.vertex.orig), "NA", as.character(max.vertex.orig))
    ))

    ## Values (only if present)
    if (!is.null(x$min.value) && !is.null(x$max.value) &&
        is.finite(x$min.value) && is.finite(x$max.value)) {
        delta.value <- x$max.value - x$min.value
        cat(sprintf(
            "Values: %.4f (min) to %.4f (max), delta = %.4f\n",
            x$min.value, x$max.value, delta.value
        ))
    }

    n.trajectories <- if (!is.null(x$n.trajectories)) as.integer(x$n.trajectories) else NA_integer_
    if (!is.na(n.trajectories)) {
        cat(sprintf("Trajectories: %d\n", n.trajectories))
    }

    ## ------------------------------------------------------------------------
    ## Mapped indices reporting (show mapped endpoints when available)
    ## ------------------------------------------------------------------------

    if (mapped) {
        cat(sprintf(
            "Vertices mapped to subgraph indices (min -> %s, max -> %s)\n",
            ifelse(is.null(x$min.vertex) || is.na(x$min.vertex), "NA", as.character(x$min.vertex)),
            ifelse(is.null(x$max.vertex) || is.na(x$max.vertex), "NA", as.character(x$max.vertex))
        ))
    }

    ## ------------------------------------------------------------------------
    ## Trajectory diagnostics
    ## ------------------------------------------------------------------------

    traj.list <- x$trajectories
    has.traj <- is.list(traj.list) && length(traj.list) > 0L

    if (has.traj) {
        lengths <- vapply(traj.list, length, integer(1))

        cat(sprintf(
            "Trajectory lengths: Min: %d, Max: %d, Mean: %.1f\n",
            min(lengths), max(lengths), mean(lengths)
        ))

        ## Unique vertices count (diagnostic)
        all.vertices <- unlist(traj.list, use.names = FALSE)
        n.unique.vertices <- length(unique(all.vertices))
        cat(sprintf("Number of vertices (unique across all trajectories): %d\n", n.unique.vertices))

        ## Optional: show first max.print trajectories
        if (max.print > 0L) {
            n.show <- min(max.print, length(traj.list))
            cat(sprintf("\nFirst %d trajectories:\n", n.show))

            for (i in seq_len(n.show)) {
                path <- traj.list[[i]]

                if (length(path) <= 10L) {
                    path.str <- paste(path, collapse = " -> ")
                } else {
                    path.str <- paste(
                        c(
                            paste(head(path, 4L), collapse = " -> "),
                            "...",
                            paste(tail(path, 3L), collapse = " -> ")
                        ),
                        collapse = " -> "
                    )
                }

                cat(sprintf("  [%d] (%d vertices): %s\n", i, length(path), path.str))
            }

            if (length(traj.list) > n.show) {
                cat(sprintf("  ... and %d more trajectories\n", length(traj.list) - n.show))
            }
        }
    }

    invisible(x)
}


#' Get All Trajectory Vertices as a Single Vector
#'
#' Extracts all unique vertices from cell trajectories.
#'
#' @param cell.traj A gfc_cell_trajectories object.
#' @param unique Logical; if TRUE (default), return unique vertices only.
#'
#' @return Integer vector of vertex indices.
#'
#' @export
trajectory.vertices.gfc.flow <- function(cell.traj, unique = TRUE) {

    if (!inherits(cell.traj, "gfc_cell_trajectories")) {
        stop("cell.traj must be a gfc_cell_trajectories object")
    }

    all.vertices <- unlist(cell.traj$trajectories)

    if (unique) {
        all.vertices <- unique(all.vertices)
    }

    return(all.vertices)
}


#' Convert Trajectory to Edge List
#'
#' Converts trajectory paths to an edge list suitable for graph operations
#' or visualization.
#'
#' @param cell.traj A gfc_cell_trajectories object.
#' @param weighted Logical; if TRUE, include edge counts as weights.
#'
#' @return A data frame with columns \code{from}, \code{to}, and optionally
#'   \code{weight}.
#'
#' @export
trajectory.edges.gfc.flow <- function(cell.traj, weighted = FALSE) {

    if (!inherits(cell.traj, "gfc_cell_trajectories")) {
        stop("cell.traj must be a gfc_cell_trajectories object")
    }

    if (cell.traj$n.trajectories == 0) {
        return(data.frame(from = integer(0), to = integer(0)))
    }

    all.edges <- list()

    for (path in cell.traj$trajectories) {
        if (length(path) < 2) next

        for (i in seq_len(length(path) - 1)) {
            from <- path[i]
            to <- path[i + 1]
            key <- paste(from, to, sep = "-")
            if (is.null(all.edges[[key]])) {
                all.edges[[key]] <- list(from = from, to = to, count = 0)
            }
            all.edges[[key]]$count <- all.edges[[key]]$count + 1
        }
    }

    if (length(all.edges) == 0) {
        return(data.frame(from = integer(0), to = integer(0)))
    }

    result <- data.frame(
        from = sapply(all.edges, function(e) e$from),
        to = sapply(all.edges, function(e) e$to)
    )

    if (weighted) {
        result$weight <- sapply(all.edges, function(e) e$count)
    }

    rownames(result) <- NULL
    return(result)
}

#' List Cells
#'
#' S3 generic for listing cells (min-max pairs) in a gradient-flow style object.
#'
#' @param x An object for which cells can be listed.
#' @param ... Passed to methods.
#'
#' @export
list.cells <- function(x, ...) {
    UseMethod("list.cells")
}

#' @export
list.cells.default <- function(x, ...) {
    stop("No list.cells() method for class: ", paste(class(x), collapse = ", "))
}


#' List All Cells in Gradient Flow Result
#'
#' Enumerates all cells (min-max pairs) in a gfc.flow result, including cells
#' with spurious endpoints (unless excluded).
#'
#' @param gfc.flow A gfc.flow object from \code{compute.gfc.flow()}.
#' @param include.spurious Logical; if TRUE (default), include cells with
#'   spurious endpoints.
#' @param loc.min Optional. If non-NULL, filter to cells with this minimum
#'   endpoint. Can be a vertex index (numeric/integer) or a label (character).
#' @param loc.max Optional. If non-NULL, filter to cells with this maximum
#'   endpoint. Can be a vertex index (numeric/integer) or a label (character).
#' @param include.spurious.flags Logical; if TRUE, include logical columns
#'   \code{min.is.spurious} and \code{max.is.spurious} in the returned data
#'   frame. Default FALSE.
#' @param ... Extra parameters.
#'
#' @return A data frame with columns:
#'   \item{min.label}{Minimum label.}
#'   \item{max.label}{Maximum label.}
#'   \item{min.vertex}{Minimum vertex index.}
#'   \item{max.vertex}{Maximum vertex index.}
#'   \item{n.trajectories}{Number of trajectories connecting this pair.}
#'   \item{unique.vertices}{Number of unique vertices in trajectories.}
#'   If \code{include.spurious.flags = TRUE}, also includes:
#'   \item{min.is.spurious}{TRUE if the minimum endpoint is spurious.}
#'   \item{max.is.spurious}{TRUE if the maximum endpoint is spurious.}
#'
#' @export
list.cells.gfc.flow <- function(x,
                                include.spurious = TRUE,
                                loc.min = NULL,
                                loc.max = NULL,
                                include.spurious.flags = FALSE,
                                ...) {

    gfc.flow <- x

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object from compute.gfc.flow()")
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        result <- data.frame(
            min.label = character(0),
            max.label = character(0),
            min.vertex = integer(0),
            max.vertex = integer(0),
            min.is.spurious = logical(0),
            max.is.spurious = logical(0),
            n.trajectories = integer(0),
            unique.vertices = integer(0),
            stringsAsFactors = FALSE
        )

        if (!include.spurious.flags) {
            result$min.is.spurious <- NULL
            result$max.is.spurious <- NULL
        }

        return(result)
    }

    summary.df <- gfc.flow$summary.all %||% gfc.flow$summary

    ## ------------------------------------------------------------------------
    ## Helper: resolve loc spec (label or vertex) to a vertex index
    ## ------------------------------------------------------------------------
    resolve.loc.to.vertex <- function(loc, labels, vertices) {
        if (is.null(loc)) return(NULL)

        ## Numeric/integer: treat as vertex index
        if (is.numeric(loc) && length(loc) == 1) {
            return(as.integer(loc))
        }

        ## Character: could be label or numeric string
        if (is.character(loc) && length(loc) == 1) {
            if (grepl("^[0-9]+$", loc)) {
                return(as.integer(loc))
            }

            ## Try resolve via summary.df label -> vertex
            if (!is.null(summary.df) && all(c("label", "vertex") %in% names(summary.df))) {
                idx <- match(loc, summary.df$label)
                if (!is.na(idx)) {
                    return(as.integer(summary.df$vertex[idx]))
                }
            }

            ## Fall back: allow filtering by label directly
            if (loc %in% labels) {
                return(NA_integer_)
            }

            stop("Could not resolve loc specification '", loc,
                 "' to a vertex. Provide a valid label or vertex index.")
        }

        stop("loc must be NULL, a single vertex index, or a single label string.")
    }

    ## ------------------------------------------------------------------------
    ## Collect cell information
    ## ------------------------------------------------------------------------
    cells <- list()

    for (traj in gfc.flow$trajectories) {
        min.v <- traj$start.vertex
        max.v <- traj$end.vertex
        key <- paste(min.v, max.v, sep = "-")

        if (is.null(cells[[key]])) {
            cells[[key]] <- list(
                min.vertex = min.v,
                max.vertex = max.v,
                min.is.spurious = traj$start.is.spurious %||% FALSE,
                max.is.spurious = traj$end.is.spurious %||% FALSE,
                n.trajectories = 0L,
                all.vertices = integer(0)
            )
        } else {
            ## If spurious flags disagree across trajectories, keep TRUE
            cells[[key]]$min.is.spurious <- cells[[key]]$min.is.spurious |
                (traj$start.is.spurious %||% FALSE)
            cells[[key]]$max.is.spurious <- cells[[key]]$max.is.spurious |
                (traj$end.is.spurious %||% FALSE)
        }

        cells[[key]]$n.trajectories <- cells[[key]]$n.trajectories + 1L
        cells[[key]]$all.vertices <- c(cells[[key]]$all.vertices, traj$vertices)
    }

    ## ------------------------------------------------------------------------
    ## Build result data frame
    ## ------------------------------------------------------------------------
    result <- data.frame(
        min.label = character(length(cells)),
        max.label = character(length(cells)),
        min.vertex = integer(length(cells)),
        max.vertex = integer(length(cells)),
        min.is.spurious = logical(length(cells)),
        max.is.spurious = logical(length(cells)),
        n.trajectories = integer(length(cells)),
        unique.vertices = integer(length(cells)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(cells)) {
        cell <- cells[[i]]

        ## Look up labels from summary
        min.idx <- match(cell$min.vertex, summary.df$vertex)
        max.idx <- match(cell$max.vertex, summary.df$vertex)

        min.label <- if (!is.na(min.idx)) summary.df$label[min.idx] else paste0("v", cell$min.vertex)
        max.label <- if (!is.na(max.idx)) summary.df$label[max.idx] else paste0("v", cell$max.vertex)

        result$min.label[i] <- min.label
        result$max.label[i] <- max.label
        result$min.vertex[i] <- as.integer(cell$min.vertex)
        result$max.vertex[i] <- as.integer(cell$max.vertex)
        result$min.is.spurious[i] <- isTRUE(cell$min.is.spurious)
        result$max.is.spurious[i] <- isTRUE(cell$max.is.spurious)
        result$n.trajectories[i] <- as.integer(cell$n.trajectories)
        result$unique.vertices[i] <- length(unique(cell$all.vertices))
    }

    ## ------------------------------------------------------------------------
    ## Optionally filter out spurious cells
    ## ------------------------------------------------------------------------
    if (!include.spurious) {
        keep <- !result$min.is.spurious & !result$max.is.spurious
        result <- result[keep, , drop = FALSE]
    }

    ## ------------------------------------------------------------------------
    ## Filter by loc.min / loc.max (vertex index or label)
    ## ------------------------------------------------------------------------
    if (!is.null(loc.min) && nrow(result) > 0) {
        v.min <- resolve.loc.to.vertex(loc.min, result$min.label, result$min.vertex)
        if (is.na(v.min)) {
            result <- result[result$min.label == loc.min, , drop = FALSE]
        } else {
            result <- result[result$min.vertex == v.min, , drop = FALSE]
        }
    }

    if (!is.null(loc.max) && nrow(result) > 0) {
        v.max <- resolve.loc.to.vertex(loc.max, result$max.label, result$max.vertex)
        if (is.na(v.max)) {
            result <- result[result$max.label == loc.max, , drop = FALSE]
        } else {
            result <- result[result$max.vertex == v.max, , drop = FALSE]
        }
    }

    ## ------------------------------------------------------------------------
    ## Sort by number of trajectories (descending)
    ## ------------------------------------------------------------------------
    if (nrow(result) > 0) {
        result <- result[order(-result$n.trajectories), , drop = FALSE]
    }
    rownames(result) <- NULL

    ## ------------------------------------------------------------------------
    ## Optionally drop spurious flag columns from the output
    ## ------------------------------------------------------------------------
    if (!include.spurious.flags) {
        result$min.is.spurious <- NULL
        result$max.is.spurious <- NULL
    }

    return(result)
}

#' List Basins
#'
#' Enumerates basin-level summaries for minima or maxima in a \code{gfc.flow}
#' result. Basin sizes are taken from \code{gfc.flow$summary} /
#' \code{gfc.flow$summary.all}. The number of incident cells is computed from
#' \code{list.cells()}.
#'
#' @param x An object for which basins can be listed.
#' @param ... Passed to methods.
#'
#' @export
list.basins <- function(x, ...) {
    UseMethod("list.basins")
}

#' @export
list.basins.default <- function(x, ...) {
    stop("No list.basins() method for class: ", paste(class(x), collapse = ", "))
}

#' List Basins for gfc.flow
#'
#' @param x A \code{gfc.flow} object from \code{compute.gfc.flow()}.
#' @param type One of \code{"min"} or \code{"max"}.
#' @param include.spurious Logical; if TRUE (default), include basins whose
#'     defining extremum is spurious.
#' @param with.rel.value Numeric. Logical; if TRUE (default), includes relative
#'     value (value/median) for local extrema.
#' @param include.spurious.flags Logical; if TRUE, include spurious-related
#'     columns from \code{summary.all} (e.g., \code{is.spurious},
#'     \code{filter.stage}, \code{merged.into}). Default FALSE.
#' @param include.absorbed Logial; if FALSE, absorbed extrema are not shown.
#'     Default FALSE.
#' @param group.by.spurious Logical; if TRUE, lists retained (refined) extrema
#'     first and then spurious extrema, each block sorted by basin size (and
#'     then ties).
#' @param loc Optional. If non-NULL, filter to a single basin extremum (vertex
#'     index or label).
#'
#' @return A data frame with basin-level columns including:
#'   \item{label}{Extremum label.}
#'   \item{vertex}{Extremum vertex index.}
#'   \item{type}{\code{"min"} or \code{"max"}.}
#'   \item{basin.size}{Number of vertices in the basin (from summary/all).}
#'   \item{n.cells}{Number of incident cells (min-max pairs).}
#'   \item{n.trajectories}{Total number of trajectories across incident cells.}
#'
#' @export
list.basins.gfc.flow <- function(x,
                                 type = c("min", "max"),
                                 include.spurious = TRUE,
                                 with.rel.value = TRUE,
                                 include.spurious.flags = FALSE,
                                 include.absorbed = FALSE,
                                 group.by.spurious = TRUE,
                                 loc = NULL,
                                 ...) {

    type <- match.arg(type)

    if (!inherits(x, "gfc.flow")) {
        stop("x must be a gfc.flow object from compute.gfc.flow()")
    }

    summary.all.df <- x$summary.all %||% x$summary
    if (is.null(summary.all.df) || nrow(summary.all.df) == 0) {
        out <- data.frame(
            label = character(0),
            vertex = integer(0),
            type = character(0),
            basin.size = integer(0),
            n.cells = integer(0),
            n.trajectories = integer(0),
            stringsAsFactors = FALSE
        )
        return(out)
    }

    ## ------------------------------------------------------------------------
    ## Build baseline basin table from summary.all
    ## ------------------------------------------------------------------------
    basins.df <- summary.all.df[summary.all.df$type == type, , drop = FALSE]

    if (isTRUE(group.by.spurious) && !("is.spurious" %in% names(basins.df))) {

        ## Fallback: infer spuriousness from label prefix (sM*, sm*)
        if ("label" %in% names(basins.df)) {
            lab <- as.character(basins.df$label)
            basins.df$is.spurious <- grepl("^(sM|sm)", lab)
        } else {
            basins.df$is.spurious <- FALSE
        }
    }

    ## Ensure required fields exist
    if (!("basin.size" %in% names(basins.df))) {
        stop("summary/all is missing 'basin.size' column")
    }

    ## Endpoint-only spurious filtering (uses summary.all$is.spurious if present)
    if (!include.spurious && ("is.spurious" %in% names(basins.df))) {
        basins.df <- basins.df[!basins.df$is.spurious, , drop = FALSE]
    }

    ## ------------------------------------------------------------------------
    ## Compute n.cells and n.trajectories from list.cells()
    ## ------------------------------------------------------------------------
    cells.df <- list.cells(x, include.spurious = TRUE, include.spurious.flags = FALSE)

    if (nrow(cells.df) > 0 && nrow(basins.df) > 0) {
        if (type == "min") {
            key.v <- cells.df$min.vertex
        } else {
            key.v <- cells.df$max.vertex
        }

        n.cells.by.v <- as.integer(table(key.v))
        v.levels <- as.integer(names(table(key.v)))

        ## Sum trajectories per extremum
        ## (table is used for counts; for sums we do a small loop for clarity)
        n.traj.by.v <- integer(length(v.levels))
        for (i in seq_along(v.levels)) {
            v <- v.levels[i]
            idx <- which(key.v == v)
            n.traj.by.v[i] <- sum(cells.df$n.trajectories[idx])
        }

        agg.df <- data.frame(
            vertex = v.levels,
            n.cells = n.cells.by.v,
            n.trajectories = n.traj.by.v,
            stringsAsFactors = FALSE
        )

        basins.df <- merge(basins.df, agg.df, by = "vertex", all.x = TRUE, sort = FALSE)
        basins.df$n.cells[is.na(basins.df$n.cells)] <- 0L
        basins.df$n.trajectories[is.na(basins.df$n.trajectories)] <- 0L
    } else {
        basins.df$n.cells <- integer(nrow(basins.df))
        basins.df$n.trajectories <- integer(nrow(basins.df))
    }

    ## ------------------------------------------------------------------------
    ## Filter by loc (vertex index or label)
    ## ------------------------------------------------------------------------
    if (!is.null(loc) && nrow(basins.df) > 0) {
        if (is.numeric(loc) && length(loc) == 1) {
            basins.df <- basins.df[basins.df$vertex == as.integer(loc), , drop = FALSE]
        } else if (is.character(loc) && length(loc) == 1) {
            basins.df <- basins.df[basins.df$label == loc, , drop = FALSE]
        } else {
            stop("loc must be NULL, a single vertex index, or a single label string.")
        }
    }

    ## ------------------------------------------------------------------------
    ## Output shaping
    ## ------------------------------------------------------------------------
    ## Keep a clean, stable column set; optionally keep spurious-related fields.
    keep.cols <- c("label", "basin.size", "n.cells", "n.trajectories")

    if (isTRUE(with.rel.value) && ("rel.value" %in% names(basins.df))) {
        keep.cols <- c(keep.cols, "rel.value")
    }

    ## Drop absorbed (size 0) if requested
    if (!isTRUE(include.absorbed) && nrow(basins.df) > 0) {
        basins.df <- basins.df[basins.df$basin.size > 0L, , drop = FALSE]
    }

    ## Sort: either grouped (retained first, spurious second) or global
    if (nrow(basins.df) > 0) {

        if (isTRUE(group.by.spurious)) {

            ## Ensure is.spurious exists (defensive)
            if (!("is.spurious" %in% names(basins.df))) {
                basins.df$is.spurious <- FALSE
            }

            ## Retained first, then spurious; each block sorted by basin.size etc.
            ord <- order(basins.df$is.spurious,
                         -basins.df$basin.size,
                         -basins.df$n.cells,
                         -basins.df$n.trajectories)

            basins.df <- basins.df[ord, , drop = FALSE]

        } else {

            ## Global sort (current behavior)
            basins.df <- basins.df[order(-basins.df$basin.size,
                                         -basins.df$n.cells,
                                         -basins.df$n.trajectories), , drop = FALSE]
        }
    }

    if (include.spurious.flags) {
        extra.cols <- intersect(c("is.spurious", "filter.stage", "merged.into"), names(basins.df))
        keep.cols <- c(keep.cols, extra.cols)
    }

    keep.cols <- intersect(keep.cols, names(basins.df))
    basins.df <- basins.df[, keep.cols, drop = FALSE]

    rownames(basins.df) <- NULL
    basins.df
}

#' Relabel Basins in a Trajectory-First Gradient Flow Complex (gfc.flow)
#'
#' @description
#' Reassigns extremum labels (e.g., \code{M1}, \code{M2}, \code{sM1}, \code{m1}, \code{sm1})
#' in a \code{gfc.flow} object so that numbering reflects the *current* basin ranking
#' (by default, descending basin size). This is particularly useful after one or more
#' rounds of basin merging with \code{merge.basins.gfc.flow()}, where some extrema may
#' be absorbed and end up with \code{basin.size == 0} and/or stale labels.
#'
#' The function updates labels consistently across:
#' \itemize{
#'   \item \code{x$max.summaries.all}, \code{x$min.summaries.all}
#'   \item \code{x$max.basins.all}, \code{x$min.basins.all}
#'   \item \code{x$summary.all} and \code{x$summary} (if present)
#'   \item overlap matrices dimnames (\code{x$max.overlap.dist}, \code{x$min.overlap.dist}) (if present)
#'   \item \code{x$merge.report} label columns (if present)
#' }
#'
#' By default, any extremum with \code{basin.size == 0} is treated as absorbed and
#' forced to be spurious (to avoid retained labels like \code{M8} having empty basins).
#'
#' @param x A \code{gfc.flow} object.
#' @param rank.by Character string specifying ranking criterion for labeling.
#'   Currently supported: \code{"basin.size"} (default). (Reserved: \code{"n.trajectories"}, \code{"value"}.)
#' @param tie.break Character vector specifying tie-break order. Supported values:
#'   \code{"n.trajectories"}, \code{"value"}, \code{"vertex"}.
#'   Default is \code{c("n.trajectories","value","vertex")}.
#' @param compute.n.trajectories Logical. If TRUE and \code{x$trajectories} exists,
#'   computes number of trajectories per extremum and uses it as a tie-breaker (and
#'   stores it internally for mapping decisions; it does not add a new column to summaries).
#' @param enforce.absorbed.spurious Logical. If TRUE (default), any extremum with
#'   \code{basin.size == 0} is forced to \code{is.spurious = TRUE}. If \code{filter.stage}
#'   is present and empty/\code{"none"}, it is set to \code{absorbed.filter.stage}.
#' @param absorbed.filter.stage Character string used for \code{filter.stage} when
#'   forcing absorbed extrema to spurious. Default is \code{"absorbed"}.
#' @param relabel.spurious Logical. If TRUE (default), spurious extrema are relabeled
#'   \code{sM1, sM2, ...} / \code{sm1, sm2, ...} according to the chosen ranking.
#'   If FALSE, spurious labels are left unchanged but retained labels are still relabeled.
#' @param store.label.map Logical. If TRUE (default), stores old->new label maps in
#'   \code{x$label.map.history}.
#' @param verbose Logical.
#'
#' @return A modified \code{gfc.flow} object with relabeled basins/summaries.
#'
#' @export
relabel.basins.gfc.flow <- function(x,
                                    rank.by = c("basin.size"),
                                    tie.break = c("n.trajectories", "value", "vertex"),
                                    compute.n.trajectories = TRUE,
                                    enforce.absorbed.spurious = TRUE,
                                    absorbed.filter.stage = "absorbed",
                                    relabel.spurious = TRUE,
                                    store.label.map = TRUE,
                                    verbose = FALSE) {

    ## ---------------------------------------------------------------------
    ## Basic validation and helpers
    ## ---------------------------------------------------------------------

    if (!inherits(x, "gfc.flow")) {
        stop("x must be a gfc.flow object from compute.gfc.flow()")
    }

    rank.by <- match.arg(rank.by)

    `%||%` <- function(a, b) if (!is.null(a)) a else b

    .is.empty.stage <- function(s) {
        is.null(s) || is.na(s) || (is.character(s) && (s == "" || s == "none"))
    }

    .safe.as.char <- function(v) {
        if (is.null(v)) return(NULL)
        if (is.character(v)) return(v)
        if (is.factor(v)) return(as.character(v))
        as.character(v)
    }

    .update.merged.into <- function(v, label.map) {
        if (is.null(v)) return(v)

        ## In your objects merged.into sometimes is integer NA; sometimes character label.
        if (is.numeric(v) || is.integer(v)) {
            return(v)
        }

        v.ch <- .safe.as.char(v)
        if (is.null(v.ch)) return(v)

        repl <- label.map[v.ch]
        idx <- which(!is.na(repl))
        if (length(idx) > 0L) v.ch[idx] <- repl[idx]
        v.ch
    }

    .order.extrema <- function(df,
                               is.maximum,
                               basin.size,
                               n.trajectories,
                               value,
                               vertex,
                               rank.by,
                               tie.break) {

        ## Primary key (currently only basin.size supported)
        ord.keys <- list()

        if (rank.by == "basin.size") {
            ord.keys <- c(ord.keys, list(-basin.size))
        } else if (rank.by == "n.trajectories") {
            ord.keys <- c(ord.keys, list(-n.trajectories))
        } else if (rank.by == "value") {
            ## maxima: high first, minima: low first
            ord.keys <- c(ord.keys, list(if (is.maximum) -value else value))
        } else {
            stop("Unsupported rank.by='", rank.by, "'.")
        }

        ## Tie-breaks
        for (tb in tie.break) {
            if (tb == "n.trajectories") {
                ord.keys <- c(ord.keys, list(-n.trajectories))
            } else if (tb == "value") {
                ord.keys <- c(ord.keys, list(if (is.maximum) -value else value))
            } else if (tb == "vertex") {
                ord.keys <- c(ord.keys, list(vertex))
            } else {
                stop("Unsupported tie.break element: '", tb, "'.")
            }
        }

        do.call(order, c(ord.keys, list(na.last = TRUE)))
    }

    .compute.ntraj <- function(trajs, type) {
        ## Returns integer named vector: names are vertex indices (as character), values are counts
        if (is.null(trajs) || length(trajs) == 0L) return(integer(0))

        if (type == "max") {
            end.v <- vapply(trajs, function(tr) as.integer(tr$end.vertex), integer(1))
            ok <- vapply(trajs, function(tr) isTRUE(tr$ends.at.lmax), logical(1))
            end.v <- end.v[ok]
            if (length(end.v) == 0L) return(integer(0))
            tab <- table(end.v)
            out <- as.integer(tab)
            names(out) <- names(tab)
            return(out)
        }

        if (type == "min") {
            start.v <- vapply(trajs, function(tr) as.integer(tr$start.vertex), integer(1))
            ok <- vapply(trajs, function(tr) isTRUE(tr$starts.at.lmin), logical(1))
            start.v <- start.v[ok]
            if (length(start.v) == 0L) return(integer(0))
            tab <- table(start.v)
            out <- as.integer(tab)
            names(out) <- names(tab)
            return(out)
        }

        stop("type must be 'max' or 'min'")
    }

    ## ---------------------------------------------------------------------
    ## Identify source summary/basin tables
    ## ---------------------------------------------------------------------

    if (is.null(x$max.summaries.all) || is.null(x$min.summaries.all) ||
        is.null(x$max.basins.all) || is.null(x$min.basins.all)) {
        stop("x must contain max/min summaries and basins: max.summaries.all, min.summaries.all, max.basins.all, min.basins.all.")
    }

    max.sum <- x$max.summaries.all
    min.sum <- x$min.summaries.all

    if (!all(c("vertex", "value", "basin.size", "label") %in% names(max.sum))) {
        stop("x$max.summaries.all is missing required columns (vertex, value, basin.size, label).")
    }
    if (!all(c("vertex", "value", "basin.size", "label") %in% names(min.sum))) {
        stop("x$min.summaries.all is missing required columns (vertex, value, basin.size, label).")
    }

    if (length(x$max.basins.all) != nrow(max.sum) || length(x$min.basins.all) != nrow(min.sum)) {
        stop("Length of basins.all does not match number of rows in summaries.all (max or min).")
    }

    ## ---------------------------------------------------------------------
    ## Compute trajectory counts (optional, for tie-breaking)
    ## ---------------------------------------------------------------------

    ntraj.max <- integer(0)
    ntraj.min <- integer(0)

    if (isTRUE(compute.n.trajectories)) {
        if (is.null(x$trajectories) || length(x$trajectories) == 0L) {
            if (verbose) {
                cat("relabel.basins.gfc.flow(): compute.n.trajectories=TRUE but x$trajectories is missing/empty; using 0 for all.\n")
            }
        } else {
            ntraj.max <- .compute.ntraj(x$trajectories, "max")
            ntraj.min <- .compute.ntraj(x$trajectories, "min")
        }
    }

    ## ---------------------------------------------------------------------
    ## Enforce absorbed->spurious if requested (basin.size == 0)
    ## ---------------------------------------------------------------------

    .enforce.absorbed <- function(sum.df, basins.all) {
        if (!isTRUE(enforce.absorbed.spurious)) return(list(sum.df = sum.df, basins.all = basins.all))

        if (!("is.spurious" %in% names(sum.df))) {
            sum.df$is.spurious <- FALSE
        }
        if (!("filter.stage" %in% names(sum.df))) {
            sum.df$filter.stage <- rep.int("none", nrow(sum.df))
        }

        for (i in seq_len(nrow(sum.df))) {
            if (isTRUE(sum.df$basin.size[i] == 0L)) {
                ## Force spurious at both levels
                sum.df$is.spurious[i] <- TRUE
                if (.is.empty.stage(sum.df$filter.stage[i])) {
                    sum.df$filter.stage[i] <- absorbed.filter.stage
                }

                if (!is.null(basins.all[[i]])) {
                    basins.all[[i]]$is.spurious <- TRUE
                    if (!is.null(basins.all[[i]]$filter.stage) && .is.empty.stage(basins.all[[i]]$filter.stage)) {
                        basins.all[[i]]$filter.stage <- absorbed.filter.stage
                    }
                }
            }
        }

        list(sum.df = sum.df, basins.all = basins.all)
    }

    tmp <- .enforce.absorbed(max.sum, x$max.basins.all)
    max.sum <- tmp$sum.df
    x$max.basins.all <- tmp$basins.all

    tmp <- .enforce.absorbed(min.sum, x$min.basins.all)
    min.sum <- tmp$sum.df
    x$min.basins.all <- tmp$basins.all

    ## Ensure is.spurious exists (after enforcement)
    if (!("is.spurious" %in% names(max.sum))) max.sum$is.spurious <- FALSE
    if (!("is.spurious" %in% names(min.sum))) min.sum$is.spurious <- FALSE

    ## ---------------------------------------------------------------------
    ## Build new labels and maps
    ## ---------------------------------------------------------------------

    .relabel.one <- function(sum.df, basins.all, type, ntraj.map) {
        is.maximum <- (type == "max")

        vertex <- as.integer(sum.df$vertex)
        value <- as.double(sum.df$value)
        basin.size <- as.integer(sum.df$basin.size)

        n.trajectories <- integer(length(vertex))
        if (length(ntraj.map) > 0L) {
            idx <- match(as.character(vertex), names(ntraj.map))
            ok <- which(!is.na(idx))
            if (length(ok) > 0L) n.trajectories[ok] <- as.integer(ntraj.map[idx[ok]])
        }

        is.spurious <- as.logical(sum.df$is.spurious)

        retained.idx <- which(!is.spurious)
        spurious.idx <- which(is.spurious)

        ## Order retained and spurious separately
        ord.ret <- retained.idx
        if (length(ord.ret) > 0L) {
            ord <- .order.extrema(sum.df[ord.ret, , drop = FALSE],
                                  is.maximum = is.maximum,
                                  basin.size = basin.size[ord.ret],
                                  n.trajectories = n.trajectories[ord.ret],
                                  value = value[ord.ret],
                                  vertex = vertex[ord.ret],
                                  rank.by = rank.by,
                                  tie.break = tie.break)
            ord.ret <- ord.ret[ord]
        }

        ord.spu <- spurious.idx
        if (length(ord.spu) > 0L) {
            ord <- .order.extrema(sum.df[ord.spu, , drop = FALSE],
                                  is.maximum = is.maximum,
                                  basin.size = basin.size[ord.spu],
                                  n.trajectories = n.trajectories[ord.spu],
                                  value = value[ord.spu],
                                  vertex = vertex[ord.spu],
                                  rank.by = rank.by,
                                  tie.break = tie.break)
            ord.spu <- ord.spu[ord]
        }

        ## Assign labels
        old.labels <- .safe.as.char(sum.df$label)
        new.labels <- old.labels

        if (is.maximum) {
            retained.prefix <- "M"
            spurious.prefix <- "sM"
        } else {
            retained.prefix <- "m"
            spurious.prefix <- "sm"
        }

        ## Retained renumbering
        if (length(ord.ret) > 0L) {
            new.labels[ord.ret] <- hooking <- paste0(retained.prefix, seq_along(ord.ret))
        }

        ## Spurious renumbering (optional)
        if (isTRUE(relabel.spurious) && length(ord.spu) > 0L) {
            new.labels[ord.spu] <- paste0(spurious.prefix, seq_along(ord.spu))
        }

        ## Build maps
        label.map <- setNames(new.labels, old.labels)
        vertex.map <- setNames(new.labels, as.character(vertex))

        ## Update summaries and basins
        sum.df$label <- as.character(new.labels)
        for (i in seq_along(basins.all)) {
            if (!is.null(basins.all[[i]])) {
                basins.all[[i]]$label <- as.character(new.labels[i])
            }
        }

        ## Update merged.into if it is a label string
        if ("merged.into" %in% names(sum.df)) {
            sum.df$merged.into <- .update.merged.into(sum.df$merged.into, label.map)
        }
        for (i in seq_along(basins.all)) {
            if (!is.null(basins.all[[i]]) && !is.null(basins.all[[i]]$merged.into)) {
                basins.all[[i]]$merged.into <- .update.merged.into(basins.all[[i]]$merged.into, label.map)
            }
        }

        list(sum.df = sum.df,
             basins.all = basins.all,
             label.map = label.map,
             vertex.map = vertex.map,
             n.trajectories = n.trajectories)
    }

    max.res <- .relabel.one(max.sum, x$max.basins.all, "max", ntraj.max)
    min.res <- .relabel.one(min.sum, x$min.basins.all, "min", ntraj.min)

    x$max.summaries.all <- max.res$sum.df
    x$min.summaries.all <- min.res$sum.df
    x$max.basins.all <- max.res$basins.all
    x$min.basins.all <- min.res$basins.all

    ## ---------------------------------------------------------------------
    ## Update summary.all and summary (if present) using vertex->label maps
    ## ---------------------------------------------------------------------

    .update.summary.df <- function(df) {
        if (is.null(df) || nrow(df) == 0L) return(df)
        if (!all(c("vertex", "type", "label") %in% names(df))) return(df)

        for (i in seq_len(nrow(df))) {
            v <- as.character(df$vertex[i])
            if (df$type[i] == "max") {
                lab <- max.res$vertex.map[v]
            } else if (df$type[i] == "min") {
                lab <- min.res$vertex.map[v]
            } else {
                lab <- NA_character_
            }

            if (!is.na(lab)) df$label[i] <- lab
        }

        if ("merged.into" %in% names(df)) {
            ## merged.into could refer to labels; update using both maps
            label.map.all <- c(max.res$label.map, min.res$label.map)
            df$merged.into <- .update.merged.into(df$merged.into, label.map.all)
        }

        df
    }

    x$summary.all <- .update.summary.df(x$summary.all)
    x$summary <- .update.summary.df(x$summary)

    ## ---------------------------------------------------------------------
    ## Update overlap matrices dimnames (if present)
    ## ---------------------------------------------------------------------

    .update.overlap.dimnames <- function(mat, label.map) {
        if (is.null(mat)) return(mat)
        dn <- dimnames(mat)
        if (is.null(dn) || length(dn) != 2L) return(mat)
        rn <- dn[[1]]
        cn <- dn[[2]]
        if (!is.null(rn)) {
            repl <- label.map[rn]
            idx <- which(!is.na(repl))
            if (length(idx) > 0L) rn[idx] <- repl[idx]
        }
        if (!is.null(cn)) {
            repl <- label.map[cn]
            idx <- which(!is.na(repl))
            if (length(idx) > 0L) cn[idx] <- repl[idx]
        }
        dimnames(mat) <- list(rn, cn)
        mat
    }

    if (!is.null(x$max.overlap.dist)) {
        x$max.overlap.dist <- .update.overlap.dimnames(x$max.overlap.dist, max.res$label.map)
    }
    if (!is.null(x$min.overlap.dist)) {
        x$min.overlap.dist <- .update.overlap.dimnames(x$min.overlap.dist, min.res$label.map)
    }

    ## ---------------------------------------------------------------------
    ## Update merge.report label columns (if present)
    ## ---------------------------------------------------------------------

    if (!is.null(x$merge.report) && is.data.frame(x$merge.report)) {
        mr <- x$merge.report

        if ("type" %in% names(mr)) {
            if ("loser.label" %in% names(mr)) {
                idx.max <- which(mr$type == "max")
                idx.min <- which(mr$type == "min")
                if (length(idx.max) > 0L) {
                    repl <- max.res$label.map[as.character(mr$loser.label[idx.max])]
                    ok <- which(!is.na(repl))
                    if (length(ok) > 0L) mr$loser.label[idx.max[ok]] <- repl[ok]
                }
                if (length(idx.min) > 0L) {
                    repl <- min.res$label.map[as.character(mr$loser.label[idx.min])]
                    ok <- which(!is.na(repl))
                    if (length(ok) > 0L) mr$loser.label[idx.min[ok]] <- repl[ok]
                }
            }

            if ("winner.label" %in% names(mr)) {
                idx.max <- which(mr$type == "max")
                idx.min <- which(mr$type == "min")
                if (length(idx.max) > 0L) {
                    repl <- max.res$label.map[as.character(mr$winner.label[idx.max])]
                    ok <- which(!is.na(repl))
                    if (length(ok) > 0L) mr$winner.label[idx.max[ok]] <- repl[ok]
                }
                if (length(idx.min) > 0L) {
                    repl <- min.res$label.map[as.character(mr$winner.label[idx.min])]
                    ok <- which(!is.na(repl))
                    if (length(ok) > 0L) mr$winner.label[idx.min[ok]] <- repl[ok]
                }
            }
        } else {
            ## If type missing, best effort: update both via combined map
            label.map.all <- c(max.res$label.map, min.res$label.map)
            if ("loser.label" %in% names(mr)) mr$loser.label <- .update.merged.into(mr$loser.label, label.map.all)
            if ("winner.label" %in% names(mr)) mr$winner.label <- .update.merged.into(mr$winner.label, label.map.all)
        }

        x$merge.report <- mr
    }

    ## ---------------------------------------------------------------------
    ## Store label maps (optional)
    ## ---------------------------------------------------------------------

    if (isTRUE(store.label.map)) {
        hist <- x$label.map.history %||% list()
        hist[[length(hist) + 1L]] <- list(
            call = match.call(),
            rank.by = rank.by,
            tie.break = tie.break,
            relabel.spurious = relabel.spurious,
            enforce.absorbed.spurious = enforce.absorbed.spurious,
            absorbed.filter.stage = absorbed.filter.stage,
            max.label.map = max.res$label.map,
            min.label.map = min.res$label.map
        )
        x$label.map.history <- hist
    }

    if (verbose) {
        n.max.ret <- sum(!as.logical(x$max.summaries.all$is.spurious))
        n.max.spu <- sum(as.logical(x$max.summaries.all$is.spurious))
        n.min.ret <- sum(!as.logical(x$min.summaries.all$is.spurious))
        n.min.spu <- sum(as.logical(x$min.summaries.all$is.spurious))

        cat("relabel.basins.gfc.flow(): relabeled maxima (retained=", n.max.ret, ", spurious=", n.max.spu, ")\n", sep = "")
        cat("relabel.basins.gfc.flow(): relabeled minima (retained=", n.min.ret, ", spurious=", n.min.spu, ")\n", sep = "")
    }

    class(x) <- unique(c("gfc.flow", class(x)))
    x
}


#' Identify Cell Membership for a Single Vertex
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertex Single integer vertex index (1-based).
#'
#' @return A data frame with all cell memberships for this vertex.
#'
#' @export
vertex.cell.gfc.flow <- function(gfc.flow, vertex) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    vertex <- as.integer(vertex)
    if (length(vertex) != 1) {
        stop("vertex must be a single integer")
    }

    if (vertex < 1 || vertex > gfc.flow$n.vertices) {
        stop(sprintf("vertex must be in range [1, %d]", gfc.flow$n.vertices))
    }

    ## Use v2 structure if available
    if (!is.null(gfc.flow$max.membership.all) && !is.null(gfc.flow$min.membership.all)) {
        max.basin.indices <- gfc.flow$max.membership.all[[vertex]]
        min.basin.indices <- gfc.flow$min.membership.all[[vertex]]

        if (length(max.basin.indices) == 0 || length(min.basin.indices) == 0) {
            return(data.frame(
                vertex = integer(0),
                min.basin.idx = integer(0),
                max.basin.idx = integer(0),
                min.label = character(0),
                max.label = character(0),
                min.is.spurious = logical(0),
                max.is.spurious = logical(0),
                stringsAsFactors = FALSE
            ))
        }

        ## Build all combinations
        combos <- expand.grid(
            min.basin.idx = min.basin.indices,
            max.basin.idx = max.basin.indices
        )

        n.cells <- nrow(combos)

        result <- data.frame(
            vertex = rep(vertex, n.cells),
            min.basin.idx = combos$min.basin.idx,
            max.basin.idx = combos$max.basin.idx,
            min.label = character(n.cells),
            max.label = character(n.cells),
            min.is.spurious = logical(n.cells),
            max.is.spurious = logical(n.cells),
            stringsAsFactors = FALSE
        )

        for (i in seq_len(n.cells)) {
            min.idx <- combos$min.basin.idx[i]
            max.idx <- combos$max.basin.idx[i]

            result$min.label[i] <- gfc.flow$min.basins.all[[min.idx]]$label %||% paste0("m", min.idx)
            result$max.label[i] <- gfc.flow$max.basins.all[[max.idx]]$label %||% paste0("M", max.idx)
            result$min.is.spurious[i] <- gfc.flow$min.basins.all[[min.idx]]$is.spurious %||% FALSE
            result$max.is.spurious[i] <- gfc.flow$max.basins.all[[max.idx]]$is.spurious %||% FALSE
        }

        return(result)
    }

    ## Fall back to trajectory-based lookup (v1 structure)
    return(.vertex.cell.gfc.flow.v1(gfc.flow, vertex))
}

#' Identify Cell Membership for Multiple Vertices
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertices Integer vector of vertex indices (1-based).
#'
#' @return A long-format data frame with all cell memberships.
#'
#' @export
vertex.cells.gfc.flow <- function(gfc.flow, vertices) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    vertices <- as.integer(vertices)

    results <- lapply(vertices, function(v) {
        vertex.cell.gfc.flow(gfc.flow, v)
    })

    do.call(rbind, results)
}


#' Get Primary Cell Membership (Single-Valued)
#'
#' Returns the first cell membership for a vertex (for backward compatibility).
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertex Single integer vertex index.
#'
#' @return A single-row data frame with primary cell membership and count.
#'
#' @export
vertex.primary.cell.gfc.flow <- function(gfc.flow, vertex) {

    all.cells <- vertex.cell.gfc.flow(gfc.flow, vertex)

    if (nrow(all.cells) == 0) {
        return(data.frame(
            vertex = vertex,
            min.label = NA_character_,
            max.label = NA_character_,
            n.total.cells = 0L,
            stringsAsFactors = FALSE
        ))
    }

    primary <- all.cells[1, , drop = FALSE]
    primary$n.total.cells <- nrow(all.cells)

    return(primary)
}

#' Count Multi-Membership Statistics
#'
#' @param gfc.flow A gfc.flow object.
#'
#' @return List with membership statistics (legacy interface).
#'
#' @export
count.multi.memberships.gfc.flow <- function(gfc.flow) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    n <- gfc.flow$n.vertices

    ## Use v2 structure if available
    if (!is.null(gfc.flow$max.membership.all) && !is.null(gfc.flow$min.membership.all)) {
        n.max.all <- sapply(gfc.flow$max.membership.all, length)
        n.min.all <- sapply(gfc.flow$min.membership.all, length)

        n.cells.all <- ifelse(
            n.max.all > 0 & n.min.all > 0,
            n.max.all * n.min.all,
            0
        )

        ## Use retained membership if available
        if (!is.null(gfc.flow$max.membership.retained) && !is.null(gfc.flow$min.membership.retained)) {
            n.max.ret <- sapply(gfc.flow$max.membership.retained, length)
            n.min.ret <- sapply(gfc.flow$min.membership.retained, length)

            n.cells.retained <- ifelse(
                n.max.ret > 0 & n.min.ret > 0,
                n.max.ret * n.min.ret,
                0
            )
        } else {
            n.cells.retained <- n.cells.all
        }

        return(list(
            n.multi.membership = sum(n.cells.retained > 1),
            n.single.membership = sum(n.cells.retained == 1),
            n.no.membership = sum(n.cells.retained == 0),
            max.memberships = max(n.cells.retained),
            multi.membership.vertices = which(n.cells.retained > 1)
        ))
    }

    ## Fall back to v1 trajectory-based counting
    n.max.memberships <- sapply(seq_len(n), function(v) {
        if (is.null(gfc.flow$max.membership)) return(0)
        length(gfc.flow$max.membership[[v]])
    })
    n.min.memberships <- sapply(seq_len(n), function(v) {
        if (is.null(gfc.flow$min.membership)) return(0)
        length(gfc.flow$min.membership[[v]])
    })

    n.cells <- ifelse(
        n.max.memberships > 0 & n.min.memberships > 0,
        n.max.memberships * n.min.memberships,
        0
    )

    list(
        n.multi.membership = sum(n.cells > 1),
        n.single.membership = sum(n.cells == 1),
        n.no.membership = sum(n.cells == 0),
        max.memberships = max(n.cells),
        multi.membership.vertices = which(n.cells > 1)
    )
}

#' Get All Trajectories Passing Through a Vertex
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertex Single integer vertex index.
#'
#' @return List of trajectory structures passing through the vertex.
#'
#' @export
vertex.all.trajectories.gfc.flow <- function(gfc.flow, vertex) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    vertex <- as.integer(vertex)
    if (length(vertex) != 1) {
        stop("vertex must be a single integer")
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        return(list())
    }

    result <- list()
    for (traj in gfc.flow$trajectories) {
        if (vertex %in% traj$vertices) {
            result[[length(result) + 1]] <- traj
        }
    }

    return(result)
}

#' Find Vertices in the Same Cell
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertex Single integer vertex index.
#' @param include.query If TRUE, include the query vertex in result.
#'
#' @return Integer vector of vertex indices in the same cell.
#'
#' @export
cell.vertices.gfc.flow <- function(gfc.flow, vertex, include.query = TRUE) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    vertex <- as.integer(vertex)
    if (length(vertex) != 1) {
        stop("vertex must be a single integer")
    }

    if (vertex < 1 || vertex > gfc.flow$n.vertices) {
        stop(sprintf("vertex must be in range [1, %d]", gfc.flow$n.vertices))
    }

    cell.info <- vertex.cell.gfc.flow(gfc.flow, vertex)

    if (nrow(cell.info) == 0) {
        return(NULL)
    }

    ## Use first cell (primary)
    min.label <- cell.info$min.label[1]
    max.label <- cell.info$max.label[1]

    ## Get vertices from summary
    summary.df <- gfc.flow$summary.all %||% gfc.flow$summary
    min.v <- summary.df$vertex[summary.df$label == min.label]
    max.v <- summary.df$vertex[summary.df$label == max.label]

    if (length(min.v) == 0 || length(max.v) == 0) {
        return(NULL)
    }

    ## Find all vertices in trajectories with this min-max pair
    all.vertices <- c()

    for (traj in gfc.flow$trajectories) {
        if (traj$start.vertex == min.v && traj$end.vertex == max.v) {
            all.vertices <- c(all.vertices, traj$vertices)
        }
    }

    all.vertices <- unique(all.vertices)

    if (!include.query) {
        all.vertices <- setdiff(all.vertices, vertex)
    }

    return(sort(all.vertices))
}


#' Get Summary Row by Label
#'
#' @param gfc.flow A gfc.flow object.
#' @param label Character label (e.g., "m1", "M2", "sm3").
#'
#' @return Data frame row from summary.all, or NULL.
#'
#' @export
get.summary <- function(gfc.flow, label) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    summary.df <- gfc.flow$summary.all %||% gfc.flow$summary

    if (is.null(summary.df)) return(NULL)

    idx <- which(summary.df$label == label)
    if (length(idx) == 0) return(NULL)

    summary.df[idx, , drop = FALSE]
}


#' Get All Cells for a Vertex (alias for vertex.cell.gfc.flow)
#'
#' @param gfc.flow A gfc.flow object.
#' @param vertex Single integer vertex index.
#'
#' @return Data frame with all cell memberships.
#'
#' @export
vertex.cell.all <- function(gfc.flow, vertex) {
    vertex.cell.gfc.flow(gfc.flow, vertex)
}

#' Get Trajectory for a Vertex Within a Specific Cell
#'
#' Returns the trajectory passing through a vertex that belongs to a specific
#' cell defined by its minimum and maximum extrema. A cell is uniquely
#' identified by its (min, max) endpoint pair.
#'
#' @param gfc.flow A gfc.flow object from \code{compute.gfc.flow()} with
#'   trajectories computed (i.e., \code{store.trajectories = TRUE}).
#' @param vertex Single integer vertex index (1-based).
#' @param min.id Minimum extremum identifier: either a label (e.g., "m4", "sm1")
#'   or a vertex index (integer, 1-based).
#' @param max.id Maximum extremum identifier: either a label (e.g., "M1", "sM2")
#'   or a vertex index (integer, 1-based).
#' @param map Optional integer vector mapping graph indices to subgraph
#'   vertices. If provided, trajectory vertices are converted to subgraph
#'   indices.
#'
#' @return A trajectory structure (list) if found, or NULL if no trajectory
#'   passes through the vertex within the specified cell. If \code{map} is
#'   provided, the returned trajectory has its vertex indices remapped.
#'
#' @details
#' This function is useful when a vertex belongs to multiple cells (due to
#' overlapping basins) and you need to retrieve the specific trajectory
#' associated with a particular cell. For vertices with unique cell membership,
#' \code{vertex.trajectory.gfc.flow()} may be simpler to use.
#'
#' @examples
#' \dontrun{
#' ## Get trajectory for vertex 100 in cell (m2, M3)
#' traj <- vertex.cell.trajectory.gfc.flow(gfc.flow.res, 100, "m2", "M3")
#'
#' ## Using vertex indices instead of labels
#' traj <- vertex.cell.trajectory.gfc.flow(gfc.flow.res, 100, 432, 718)
#'
#' ## With subgraph mapping
#' traj <- vertex.cell.trajectory.gfc.flow(gfc.flow.res, 100, "m2", "M3",
#'                                         map = subgraph.map)
#' }
#'
#' @seealso \code{\link{vertex.trajectory.gfc.flow}},
#'   \code{\link{vertex.cell.gfc.flow}}, \code{\link{cell.trajectories.gfc.flow}}
#'
#' @export
vertex.cell.trajectory.gfc.flow <- function(gfc.flow, vertex, min.id, max.id,
                                            map = NULL) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    vertex <- as.integer(vertex)
    if (length(vertex) != 1) {
        stop("vertex must be a single integer")
    }

    if (vertex < 1 || vertex > gfc.flow$n.vertices) {
        stop(sprintf("vertex must be in range [1, %d]", gfc.flow$n.vertices))
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        warning("gfc.flow does not contain trajectory data")
        return(NULL)
    }

    ## -------------------------------------------------------------------------
    ## Resolve min.id to vertex index
    ## -------------------------------------------------------------------------

    if (is.character(min.id)) {
        ## Label provided - look up in summary
        summary.df <- gfc.flow$summary.all %||% gfc.flow$summary
        if (is.null(summary.df)) {
            stop("No summary table available to resolve label")
        }
        min.row <- summary.df[summary.df$label == min.id, ]
        if (nrow(min.row) == 0) {
            stop(sprintf("Minimum label '%s' not found in summary", min.id))
        }
        min.vertex <- min.row$vertex[1]
    } else {
        ## Vertex index provided
        min.vertex <- as.integer(min.id)
        if (min.vertex < 1 || min.vertex > gfc.flow$n.vertices) {
            stop(sprintf("min.id vertex must be in range [1, %d]",
                         gfc.flow$n.vertices))
        }
    }

    ## -------------------------------------------------------------------------
    ## Resolve max.id to vertex index
    ## -------------------------------------------------------------------------

    if (is.character(max.id)) {
        ## Label provided - look up in summary
        summary.df <- gfc.flow$summary.all %||% gfc.flow$summary
        if (is.null(summary.df)) {
            stop("No summary table available to resolve label")
        }
        max.row <- summary.df[summary.df$label == max.id, ]
        if (nrow(max.row) == 0) {
            stop(sprintf("Maximum label '%s' not found in summary", max.id))
        }
        max.vertex <- max.row$vertex[1]
    } else {
        ## Vertex index provided
        max.vertex <- as.integer(max.id)
        if (max.vertex < 1 || max.vertex > gfc.flow$n.vertices) {
            stop(sprintf("max.id vertex must be in range [1, %d]",
                         gfc.flow$n.vertices))
        }
    }

    ## -------------------------------------------------------------------------
    ## Find trajectory matching cell and containing vertex
    ## -------------------------------------------------------------------------

    for (traj in gfc.flow$trajectories) {
        if (traj$start.vertex == min.vertex &&
            traj$end.vertex == max.vertex &&
            vertex %in% traj$vertices) {

            ## Found matching trajectory
            if (!is.null(map)) {
                ## Remap vertices to subgraph indices
                traj$vertices <- map[traj$vertices]
                traj$start.vertex <- map[traj$start.vertex]
                traj$end.vertex <- map[traj$end.vertex]
            }

            return(traj)
        }
    }

    ## No matching trajectory found
    return(NULL)
}

#' Get Uncovered Vertices from GFC Flow Result
#'
#' Returns information about vertices that are not assigned to any basin,
#' typically because they are isolated by the edge length threshold. This
#' includes identification of any local extrema among the uncovered vertices.
#'
#' @param gfc.flow A gfc.flow object from \code{compute.gfc.flow()}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{vertices}{Integer vector of uncovered vertex indices (1-based).}
#'     \item{n.uncovered}{Total number of uncovered vertices.}
#'     \item{n.total}{Total number of vertices in the graph.}
#'     \item{coverage}{Fraction of vertices covered (assigned to at least
#'       one basin).}
#'     \item{local.minima}{Integer vector of uncovered vertices that are
#'       local minima.}
#'     \item{local.maxima}{Integer vector of uncovered vertices that are
#'       local maxima.}
#'     \item{n.local.minima}{Number of uncovered local minima.}
#'     \item{n.local.maxima}{Number of uncovered local maxima.}
#'     \item{by.type}{Data frame with columns: vertex, type ("regular",
#'       "lmin", "lmax").}
#'   }
#'
#' @details
#' Uncovered vertices arise when the edge length threshold prevents
#' trajectories from reaching valid endpoints. These vertices are typically
#' isolated outliers connected to the main graph only through atypically
#' long edges.
#'
#' Local extrema among uncovered vertices represent minima or maxima that
#' were identified in the initial extrema detection but could not form
#' valid basins because no trajectories with valid endpoints passed
#' through them.
#'
#' @examples
#' \dontrun{
#' result <- compute.gfc.flow(adj.list, weight.list, y,
#'                            edge.length.quantile.thld = 0.9)
#' uncov <- uncovered.vertices.gfc.flow(result)
#' cat("Uncovered:", uncov$n.uncovered, "vertices\n")
#' cat("Including:", uncov$n.local.minima, "local minima,",
#'     uncov$n.local.maxima, "local maxima\n")
#' }
#'
#' @seealso \code{\link{compute.gfc.flow}}, \code{\link{summary.gfc.flow}}
#'
#' @export
uncovered.vertices.gfc.flow <- function(gfc.flow) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object")
    }

    n <- gfc.flow$n.vertices

    ## -------------------------------------------------------------------------
    ## Find all vertices in any basin
    ## -------------------------------------------------------------------------

    covered.vertices <- integer(0)

    if (!is.null(gfc.flow$min.basins.all)) {
        for (basin in gfc.flow$min.basins.all) {
            covered.vertices <- c(covered.vertices, basin$vertices)
        }
    }

    if (!is.null(gfc.flow$max.basins.all)) {
        for (basin in gfc.flow$max.basins.all) {
            covered.vertices <- c(covered.vertices, basin$vertices)
        }
    }

    covered.vertices <- unique(covered.vertices)
    uncovered.vertices <- setdiff(seq_len(n), covered.vertices)

    ## -------------------------------------------------------------------------
    ## Identify local extrema from trajectories
    ## -------------------------------------------------------------------------

    all.lmin <- integer(0)
    all.lmax <- integer(0)

    if (!is.null(gfc.flow$trajectories)) {
        for (traj in gfc.flow$trajectories) {
            if (isTRUE(traj$starts.at.lmin)) {
                all.lmin <- c(all.lmin, traj$start.vertex)
            }
            if (isTRUE(traj$ends.at.lmax)) {
                all.lmax <- c(all.lmax, traj$end.vertex)
            }
        }
        all.lmin <- unique(all.lmin)
        all.lmax <- unique(all.lmax)
    }

    ## -------------------------------------------------------------------------
    ## Find uncovered local extrema
    ## -------------------------------------------------------------------------

    uncovered.lmin <- intersect(uncovered.vertices, all.lmin)
    uncovered.lmax <- intersect(uncovered.vertices, all.lmax)

    ## -------------------------------------------------------------------------
    ## Build type classification
    ## -------------------------------------------------------------------------

    if (length(uncovered.vertices) > 0) {
        types <- rep("regular", length(uncovered.vertices))
        types[uncovered.vertices %in% uncovered.lmin] <- "lmin"
        types[uncovered.vertices %in% uncovered.lmax] <- "lmax"

        by.type <- data.frame(
            vertex = uncovered.vertices,
            type = types,
            stringsAsFactors = FALSE
        )
    } else {
        by.type <- data.frame(
            vertex = integer(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }

    ## -------------------------------------------------------------------------
    ## Return result
    ## -------------------------------------------------------------------------

    list(
        vertices = uncovered.vertices,
        n.uncovered = length(uncovered.vertices),
        n.total = n,
        coverage = 1 - length(uncovered.vertices) / n,
        local.minima = uncovered.lmin,
        local.maxima = uncovered.lmax,
        n.local.minima = length(uncovered.lmin),
        n.local.maxima = length(uncovered.lmax),
        by.type = by.type
    )
}

#' Compute Maximum Edge Length for Each Cell Trajectory
#'
#' For each trajectory in a cell.trajectories object, computes the maximum
#' edge length along the trajectory path. This is useful for identifying
#' trajectories that traverse unusually long edges, which may indicate
#' paths through isolated or outlier vertices.
#'
#' @param adj.list List of integer vectors. Each element \code{adj.list[[i]]}
#'   contains the 1-based indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors. Each element
#'   \code{weight.list[[i]]} contains the edge weights (distances) corresponding
#'   to \code{adj.list[[i]]}.
#' @param cell.trajectories Output from \code{cell.trajectories.gfc.flow()},
#'   a list of trajectory objects each containing a \code{vertices} field.
#'
#' @return Numeric vector of length equal to \code{length(cell.trajectories)},
#'   where the i-th element is the maximum edge length along the i-th
#'   trajectory. Returns \code{NA} for trajectories with fewer than 2 vertices
#'   or if an edge is not found in the adjacency structure.
#'
#' @examples
#' \dontrun{
#' cell.trajs <- cell.trajectories.gfc.flow(gfc.flow.res, "m2", "M7")
#' max.lens <- cell.trajectories.max.edge.length(adj.list, weight.list, cell.trajs)
#'
#' ## Find trajectories with edges exceeding the threshold
#' threshold <- gfc.flow.res$edge.length.thld
#' which(max.lens > threshold)
#' }
#'
#' @seealso \code{\link{cell.trajectories.gfc.flow}}
#'
#' @export
cell.trajectories.max.edge.length <- function(adj.list,
                                              weight.list,
                                              cell.trajectories) {

    n.traj <- length(cell.trajectories$trajectories)

    if (n.traj == 0) {
        return(numeric(0))
    }

    max.edge.lengths <- numeric(n.traj)

    for (i in seq_len(n.traj)) {
        vertices <- cell.trajectories$trajectories[[i]]
        n.v <- length(vertices)

        if (n.v < 2) {
            max.edge.lengths[i] <- NA_real_
            next
        }

        max.len <- 0
        edge.not.found <- FALSE

        for (j in seq_len(n.v - 1)) {
            v1 <- vertices[j]
            v2 <- vertices[j + 1]

            ## Find edge weight from v1 to v2
            nbrs <- adj.list[[v1]]
            idx <- which(nbrs == v2)

            if (length(idx) == 0) {
                warning(sprintf("Edge (%d, %d) not found in trajectory %d",
                                v1, v2, i))
                edge.not.found <- TRUE
                break
            }

            edge.len <- weight.list[[v1]][idx[1]]
            if (edge.len > max.len) {
                max.len <- edge.len
            }
        }

        max.edge.lengths[i] <- if (edge.not.found) NA_real_ else max.len
    }

    return(max.edge.lengths)
}

#' Compute Monotonicity Statistics for Cell Trajectories
#'
#' For each trajectory in a cell.trajectories object, evaluates whether the
#' graph distance from the starting vertex increases monotonically along the
#' trajectory. This diagnostic helps identify trajectories that "wander" away
#' from and back toward the starting extremum rather than progressing smoothly.
#'
#' @param adj.list List of integer vectors. Each element \code{adj.list[[i]]}
#'   contains the 1-based indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors. Each element
#'   \code{weight.list[[i]]} contains the edge weights (distances) corresponding
#'   to \code{adj.list[[i]]}.
#' @param cell.trajectories Output from \code{cell.trajectories.gfc.flow()},
#'   containing a \code{trajectories} element which is a list of trajectory
#'   objects each with a \code{vertices} field.
#'
#' @return A data frame with one row per trajectory and three columns:
#'   \describe{
#'     \item{n.mono}{Number of consecutive vertex pairs where the distance
#'       from the start vertex increases, i.e., where
#'       \eqn{d(x_0, v_{j+1}) > d(x_0, v_j)}.}
#'     \item{n.vert}{Number of consecutive vertex pairs in the trajectory
#'       (equal to the number of vertices minus one).}
#'     \item{p.mono}{Proportion of monotonic transitions: \code{n.mono / n.vert}.
#'       A value of 1.0 indicates perfect distance monotonicity.}
#'   }
#'
#' @details
#' The function computes shortest path distances from each trajectory's
#' starting vertex using igraph's optimized shortest path algorithms.
#' For trajectories sharing the same start vertex (common within a cell),
#' distances are computed once and cached for efficiency.
#'
#' A trajectory with \code{p.mono = 1.0} is fully monotonic: every step
#' along the trajectory moves further from the starting extremum in terms
#' of graph distance. Lower values indicate trajectories that backtrack
#' or take detours, which may suggest paths through geometrically isolated
#' vertices.
#'
#' @examples
#' \dontrun{
#' cell.trajs <- cell.trajectories.gfc.flow(gfc.flow.res, "m2", "M7")
#' mono.stats <- cell.trajectories.monotonicity(adj.list, weight.list, cell.trajs)
#'
#' ## Find non-monotonic trajectories
#' which(mono.stats$p.mono < 1.0)
#'
#' ## Summary of monotonicity across all trajectories
#' summary(mono.stats$p.mono)
#' }
#'
#' @seealso \code{\link{cell.trajectories.gfc.flow}},
#'   \code{\link{cell.trajectories.max.edge.length}}
#'
#' @export
cell.trajectories.monotonicity <- function(adj.list,
                                           weight.list,
                                           cell.trajectories) {

    trajs <- cell.trajectories$trajectories
    n.traj <- length(trajs)

    if (n.traj == 0) {
        return(data.frame(
            n.mono = integer(0),
            n.vert = integer(0),
            p.mono = numeric(0)
        ))
    }

    ## Build igraph object
    graph.obj <- convert.adjacency.to.edge.matrix(adj.list, weight.list)
    g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
    igraph::E(g)$weight <- graph.obj$weights

    ## Cache distances by start vertex to avoid recomputation
    dist.cache <- list()

    n.mono <- integer(n.traj)
    n.vert <- integer(n.traj)

    for (i in seq_len(n.traj)) {
        vertices <- trajs[[i]]
        n.v <- length(vertices)

        if (n.v < 2) {
            n.mono[i] <- 0L
            n.vert[i] <- 0L
            next
        }

        start <- vertices[1]
        start.key <- as.character(start)

        ## Get or compute distances from start vertex
        if (is.null(dist.cache[[start.key]])) {
            dist.cache[[start.key]] <- as.numeric(
                igraph::distances(g, v = start, mode = "all")
            )
        }
        dist.from.start <- dist.cache[[start.key]]

        ## Count monotonic transitions
        n.vert[i] <- n.v - 1L
        mono.count <- 0L

        for (j in seq_len(n.v - 1)) {
            d.current <- dist.from.start[vertices[j]]
            d.next <- dist.from.start[vertices[j + 1]]
            if (d.next > d.current) {
                mono.count <- mono.count + 1L
            }
        }

        n.mono[i] <- mono.count
    }

    p.mono <- ifelse(n.vert > 0, n.mono / n.vert, NA_real_)

    data.frame(
        n.mono = n.mono,
        n.vert = n.vert,
        p.mono = p.mono
    )
}


#' Construct a Connected Graph Path Through Waypoint Vertices
#'
#' Given a set of waypoint vertices in a specified order, constructs a connected
#' path in the graph by computing shortest paths between consecutive waypoints
#' and concatenating the results. This is useful when you have a set of vertices
#' that should be traversed in a particular order but are not necessarily
#' adjacent in the graph.
#'
#' The function iteratively computes shortest paths between consecutive
#' waypoints using edge weights if available. The resulting path visits all
#' waypoints in the specified order, with intermediate vertices inserted as
#' needed to maintain graph connectivity.
#'
#' @param igraph.obj An igraph graph object. Should be connected, at least
#'   among the waypoint vertices. Edge weights, if present, are used for
#'   shortest path computation.
#' @param waypoints Integer vector of vertex indices specifying the waypoints
#'   to connect, in the desired traversal order. Must contain at least two
#'   vertices. Vertex indices should be valid for the given graph (1-based).
#'
#' @return Integer vector of vertex indices forming a connected path through
#'   all waypoints. The path starts at \code{waypoints[1]} and ends at
#'   \code{waypoints[length(waypoints)]}, visiting all intermediate waypoints
#'   in order. The returned path may contain vertices not in the original
#'   waypoints if shortest paths require intermediate steps.
#'
#' @details
#' When waypoints are derived from distance-ordering (e.g., sorting vertices
#' by graph distance from a reference point), consecutive waypoints in the
#' sorted order are typically close in the graph but not necessarily adjacent.
#' This function resolves that gap by inserting the shortest path between
#' each consecutive pair.
#'
#' The function assumes the graph is connected between all waypoint pairs.
#' If any two consecutive waypoints are in disconnected components, the
#' shortest path computation will fail.
#'
#' @seealso \code{\link{compute.harmonic.extension}} which requires a connected
#'   path as input
#'
#' @examples
#' \dontrun{
#' ## Suppose we have vertices sorted by distance from a reference point
#' ## but they don't form a connected path
#' ordered.vertices <- c(10, 25, 18, 42, 37)
#'
#' ## Construct a connected path through these waypoints
#' connected.path <- construct.path.through.waypoints(gr, ordered.vertices)
#'
#' ## Now the path can be used for harmonic extension
#' hext <- compute.harmonic.extension(
#'     adj.list = adj.list,
#'     weight.list = weight.list,
#'     trajectory = connected.path,
#'     tube.radius = 1
#' )
#' }
#'
#' @export
construct.path.through.waypoints <- function(igraph.obj, waypoints) {
    ## Start with the first waypoint
    connected.path <- waypoints[1]

    for (i in seq_len(length(waypoints) - 1)) {
        from.v <- waypoints[i]
        to.v <- waypoints[i + 1]

        ## Compute shortest path between consecutive waypoints
        sp <- igraph::shortest_paths(
            igraph.obj,
            from = from.v,
            to = to.v,
            weights = igraph::E(igraph.obj)$weight,
            output = "vpath"
        )

        ## Extract path vertices (excluding the starting vertex to avoid duplicates)
        path.vertices <- as.integer(sp$vpath[[1]])
        if (length(path.vertices) > 1) {
            connected.path <- c(connected.path, path.vertices[-1])
        }
    }

    return(connected.path)
}

#' Compute Trajectory Distances to a Target Vertex
#'
#' For each trajectory in a list, computes the minimum graph distance from any
#' vertex in the trajectory to a specified target vertex. This is useful for
#' identifying trajectories that pass near a landmark vertex (e.g., a community
#' state type centroid) without necessarily including it.
#'
#' @param trajectories A list of integer vectors, where each vector contains
#'   vertex indices forming a trajectory. Vertex indices should be 1-based.
#' @param target.vertex Integer specifying the target vertex index (1-based).
#' @param adj.list List of integer vectors specifying the adjacency structure.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices
#'   adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors specifying edge weights
#'   (lengths). Each element \code{weight.list[[i]]} contains the weights
#'   of edges from vertex \code{i} to its neighbors in \code{adj.list[[i]]}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{min.distances}{Numeric vector of length equal to the number of
#'       trajectories. Each element is the minimum graph distance from any
#'       vertex in the corresponding trajectory to the target vertex.}
#'     \item{nearest.vertices}{Integer vector of the trajectory vertex that
#'       achieves the minimum distance to the target for each trajectory.}
#'     \item{nearest.positions}{Integer vector of the position within each
#'       trajectory where the nearest vertex occurs (1-based).}
#'     \item{target.vertex}{The target vertex index (for reference).}
#'   }
#'
#' @details
#' The function constructs an igraph object from the adjacency and weight lists,
#' then computes shortest path distances from the target vertex to all vertices
#' in the graph. For each trajectory, it identifies the vertex with the smallest
#' distance to the target and records both the distance and the vertex identity.
#'
#' Trajectories that pass through the target vertex will have a minimum distance
#' of zero.
#'
#' @seealso \code{\link{vertex.all.trajectories.gfc.flow}} for extracting
#'   trajectories from a GFC flow result, \code{\link{cell.trajectories}} for
#'   extracting trajectories within a specific cell
#'
#' @examples
#' \dontrun{
#' ## Get all trajectories descending to minimum m2
#' m2.trajs <- vertex.all.trajectories.gfc.flow(gfc.res, m2.vertex)
#'
#' ## Find which trajectories pass closest to the Li landmark
#' traj.dists <- trajectory.distances.to.vertex(
#'     trajectories = m2.trajs,
#'     target.vertex = Li.vertex,
#'     adj.list = adj.list,
#'     weight.list = weight.list
#' )
#'
#' ## Select trajectories within distance 0.5 of Li
#' near.Li.idx <- which(traj.dists$min.distances < 0.5)
#' }
#'
#' @export
trajectory.distances.to.vertex <- function(trajectories,
                                           target.vertex,
                                           adj.list,
                                           weight.list) {

    ## Input validation
    if (!is.list(trajectories)) {
        stop("trajectories must be a list")
    }

    if (length(trajectories) == 0) {
        return(list(
            min.distances = numeric(0),
            nearest.vertices = integer(0),
            nearest.positions = integer(0),
            target.vertex = target.vertex
        ))
    }

    target.vertex <- as.integer(target.vertex)
    if (length(target.vertex) != 1) {
        stop("target.vertex must be a single integer")
    }

    n.vertices <- length(adj.list)
    if (target.vertex < 1 || target.vertex > n.vertices) {
        stop("target.vertex must be between 1 and the number of vertices")
    }

    ## Build adjacency matrix and igraph object
    adj.mat <- convert.adjacency.list.to.adjacency.matrix(
        adj.list,
        weight.list
    )

    gr <- igraph::graph_from_adjacency_matrix(
        adj.mat,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )

    ## Compute distances from target vertex to all vertices in graph
    all.dists <- igraph::distances(
        gr,
        v = target.vertex,
        to = igraph::V(gr),
        weights = igraph::E(gr)$weight
    )[1, ]

    ## For each trajectory, find minimum distance to target
    n.trajs <- length(trajectories)
    min.distances <- numeric(n.trajs)
    nearest.vertices <- integer(n.trajs)
    nearest.positions <- integer(n.trajs)

    for (i in seq_len(n.trajs)) {
        traj <- trajectories[[i]]
        traj.dists <- all.dists[traj$vertices]
        min.pos <- which.min(traj.dists)

        min.distances[i] <- traj.dists[min.pos]
        nearest.vertices[i] <- traj$vertices[min.pos]
        nearest.positions[i] <- min.pos
    }

    return(list(
        min.distances = min.distances,
        nearest.vertices = nearest.vertices,
        nearest.positions = nearest.positions,
        target.vertex = target.vertex
    ))
}

#' Select Best Trajectory by Distance and Length
#'
#' Among trajectories closest to a target vertex, selects the one with the
#' most vertices. This is useful when multiple trajectories are equally close
#' to a landmark and a tiebreaker is needed.
#'
#' @param trajectories A list of integer vectors, each representing a trajectory.
#' @param traj.dists Output from \code{\link{trajectory.distances.to.vertex}}.
#' @param distance.tolerance Numeric tolerance for considering distances equal.
#'   Trajectories within this tolerance of the minimum distance are considered
#'   tied. Default is 1e-9.
#'
#' @return A list with components:
#'   \describe{
#'     \item{best.idx}{Index of the selected trajectory in the input list.}
#'     \item{trajectory}{The selected trajectory (integer vector of vertices).}
#'     \item{min.distance}{Distance from the selected trajectory to the target.}
#'     \item{n.vertices}{Number of vertices in the selected trajectory.}
#'     \item{n.tied}{Number of trajectories tied at the minimum distance.}
#'     \item{tied.idx}{Indices of all tied trajectories.}
#'   }
#'
#' @export
select.closest.longest.trajectory <- function(trajectories,
                                              traj.dists,
                                              distance.tolerance = 1e-9) {

    min.dist <- min(traj.dists$min.distances)
    tied.idx <- which(abs(traj.dists$min.distances - min.dist) < distance.tolerance)

    ## Among tied trajectories, find the longest
    tied.lengths <- sapply(trajectories[tied.idx], length)
    best.local.idx <- which.max(tied.lengths)
    best.idx <- tied.idx[best.local.idx]

    return(list(
        best.idx = best.idx,
        trajectory = trajectories[[best.idx]],
        min.distance = traj.dists$min.distances[best.idx],
        n.vertices = length(trajectories[[best.idx]]$vertices),
        n.tied = length(tied.idx),
        tied.idx = tied.idx
    ))
}

#' Compute Vertex Distances to a Trajectory
#'
#' For each vertex in a given set, computes the minimum graph distance to any
#' vertex in a trajectory. This is useful for characterizing how far vertices
#' in a tubular neighborhood are from the central path.
#'
#' @param vertices Integer vector of vertex indices for which to compute
#'   distances. Vertex indices should be 1-based.
#' @param trajectory.vertices Integer vector of vertex indices forming the
#'   trajectory (or any reference vertex set). Vertex indices should be 1-based.
#' @param adj.list List of integer vectors specifying the adjacency structure.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices
#'   adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors specifying edge weights
#'   (lengths). Each element \code{weight.list[[i]]} contains the weights
#'   of edges from vertex \code{i} to its neighbors in \code{adj.list[[i]]}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{min.distances}{Numeric vector of length equal to the number of
#'       input vertices. Each element is the minimum graph distance from the
#'       corresponding vertex to any vertex in the trajectory.}
#'     \item{nearest.traj.vertices}{Integer vector of the trajectory vertex
#'       that is closest to each input vertex.}
#'     \item{nearest.traj.positions}{Integer vector of the position within
#'       the trajectory where the nearest vertex occurs (1-based). Useful
#'       when trajectory order is meaningful.}
#'     \item{vertices}{The input vertices (for reference).}
#'     \item{trajectory.vertices}{The trajectory vertices (for reference).}
#'   }
#'
#' @details
#' The function constructs an igraph object from the adjacency and weight lists,
#' then computes shortest path distances from all trajectory vertices to all
#' input vertices. For each input vertex, it identifies the trajectory vertex
#' with the smallest distance.
#'
#' Vertices that are part of the trajectory will have a minimum distance of zero.
#'
#' @seealso \code{\link{trajectory.distances.to.vertex}} for the inverse
#'   operation (distance from trajectories to a single vertex),
#'   \code{\link{compute.harmonic.extension}} which produces tubular
#'   neighborhoods that can be analyzed with this function
#'
#' @examples
#' \dontrun{
#' ## Compute harmonic extension to get tubular neighborhood
#' hext <- compute.harmonic.extension(adj.list, weight.list, trajectory,
#'                                    tube.radius = 2)
#'
#' ## Get distances from tubular vertices to the trajectory
#' tube.dists <- vertex.distances.to.trajectory(
#'     vertices = hext$tubular.vertices,
#'     trajectory.vertices = trajectory,
#'     adj.list = adj.list,
#'     weight.list = weight.list
#' )
#'
#' ## Vertices on the trajectory have distance 0
#' sum(tube.dists$min.distances == 0)  # equals length(trajectory)
#' }
#'
#' @export
vertex.distances.to.trajectory <- function(vertices,
                                           trajectory.vertices,
                                           adj.list,
                                           weight.list) {

    ## Input validation
    if (length(vertices) == 0) {
        return(list(
            min.distances = numeric(0),
            nearest.traj.vertices = integer(0),
            nearest.traj.positions = integer(0),
            vertices = vertices,
            trajectory.vertices = trajectory.vertices
        ))
    }

    if (length(trajectory.vertices) == 0) {
        stop("trajectory.vertices must contain at least one vertex")
    }

    vertices <- as.integer(vertices)
    trajectory.vertices <- as.integer(trajectory.vertices)

    n.graph.vertices <- length(adj.list)

    if (any(vertices < 1) || any(vertices > n.graph.vertices)) {
        stop("all vertices must be between 1 and the number of graph vertices")
    }

    if (any(trajectory.vertices < 1) || any(trajectory.vertices > n.graph.vertices)) {
        stop("all trajectory.vertices must be between 1 and the number of graph vertices")
    }

    ## Build adjacency matrix and igraph object
    adj.mat <- convert.adjacency.list.to.adjacency.matrix(adj.list, weight.list)

    gr <- igraph::graph_from_adjacency_matrix(
        adj.mat,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )

    ## Compute distances from all trajectory vertices to all input vertices
    ## Result is a matrix: rows = trajectory vertices, cols = input vertices
    dist.mat <- igraph::distances(
        gr,
        v = trajectory.vertices,
        to = vertices,
        weights = igraph::E(gr)$weight
    )

    ## For each input vertex (column), find the trajectory vertex (row) with min distance
    n.vertices <- length(vertices)
    min.distances <- numeric(n.vertices)
    nearest.traj.vertices <- integer(n.vertices)
    nearest.traj.positions <- integer(n.vertices)

    for (j in seq_len(n.vertices)) {
        min.row.idx <- which.min(dist.mat[, j])
        min.distances[j] <- dist.mat[min.row.idx, j]
        nearest.traj.vertices[j] <- trajectory.vertices[min.row.idx]
        nearest.traj.positions[j] <- min.row.idx
    }

    return(list(
        min.distances = min.distances,
        nearest.traj.vertices = nearest.traj.vertices,
        nearest.traj.positions = nearest.traj.positions,
        vertices = vertices,
        trajectory.vertices = trajectory.vertices
    ))
}

#' Merge Basins in a Trajectory-First Gradient Flow Complex (gfc.flow)
#'
#' @description
#' Merges (absorbs) a set of "loser" extrema basins into corresponding "winner"
#' extrema basins by modifying trajectories: each trajectory that ends at a loser
#' extremum is extended by a connector path from the loser to the winner.
#'
#' The connector is computed on the domain graph. Initially, only the weighted
#' shortest-path connector is implemented (via \pkg{igraph}). The public API is
#' designed so that \code{"high.saddle"} and \code{"high.saddle.shortest"} can be
#' added later without refactoring.
#'
#' @details
#' For each merge pair \eqn{lM \to wM} (or \eqn{lm \to wm} for minima), and for each
#' trajectory \eqn{T=(u_0,\ldots,u_m=lM)} ending at the loser, the function replaces it
#' with \eqn{T'=(u_0,\ldots,u_m=lM,\ldots,wM)} by appending a connector path
#' \eqn{P(lM,wM)} (excluding the first vertex to avoid duplication).
#'
#' After trajectory rewrites, the function rebuilds:
#' \itemize{
#'   \item \code{max.membership.all}, \code{min.membership.all} from trajectories
#'   \item \code{max.basins.all}, \code{min.basins.all} vertices and hop distances
#'   \item \code{basin.size} in \code{min.summaries.all} / \code{max.summaries.all}
#'   \item \code{summary.all} and \code{summary} via \code{.postprocess.gfc.flow()}
#' }
#'
#' Note: fields such as \code{max.overlap.dist} / \code{min.overlap.dist} are not
#' recomputed here and should be treated as stale after merging.
#'
#' @param x A \code{gfc.flow} object from \code{compute.gfc.flow()}.
#' @param adj.list Adjacency list of the domain graph (1-based vertex indices).
#' @param weight.list Parallel list of edge weights; \code{weight.list[[i]]}
#'     must align with \code{adj.list[[i]]}.
#' @param merge.map Named character/integer vector specifying merges. Names are
#'     losers, values are winners. Each element can be a label (e.g.,
#'     \code{"M4"}) or a vertex index (e.g., \code{1814}).
#' @param losers Optional vector of losers (labels or vertex indices). Used if
#'     \code{merge.map} is not supplied.
#' @param winners Optional vector of winners (labels or vertex indices). If
#'     length 1 and \code{losers} has length > 1, the winner is recycled.
#' @param type Extremum type to merge: \code{"max"} or \code{"min"}.
#' @param connector.method Connector policy. Currently only \code{"closest"} and
#'     \code{"shortest"} methods are implemented. \code{"high.saddle"} and
#'     \code{"high.saddle.shortest"} are reserved for future implementations.
#' @param mark.loser.spurious Logical. If TRUE, marks loser extrema as spurious
#'     and sets \code{merged.into} to the winner label in summaries and basin
#'     records.
#' @param merge.filter.stage Character string to write to \code{filter.stage}
#'     for merged losers (only used when \code{mark.loser.spurious=TRUE}).
#' @param verbose Logical.
#'
#' @return A modified \code{gfc.flow} object with updated trajectories and derived
#'   basin/membership/summary fields.
#'
#' @export
merge.basins.gfc.flow <- function(x,
                                  adj.list,
                                  weight.list,
                                  merge.map = NULL,
                                  losers = NULL,
                                  winners = NULL,
                                  type = c("max", "min"),
                                  connector.method = c("closest",
                                                       "shortest",
                                                       "high.saddle",
                                                       "high.saddle.shortest"),
                                  mark.loser.spurious = TRUE,
                                  merge.filter.stage = "MERGE",
                                  verbose = TRUE) {

    ## ------------------------------------------------------------------------
    ## Basic validation
    ## ------------------------------------------------------------------------

    type <- match.arg(type)
    connector.method <- match.arg(connector.method)

    if (!inherits(x, "gfc.flow")) {
        stop("x must be a gfc.flow object from compute.gfc.flow()")
    }

    if (is.null(adj.list) || is.null(weight.list)) {
        stop("adj.list and weight.list must be supplied (domain graph definition).")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length.")
    }

    if (is.null(x$trajectories) || length(x$trajectories) == 0) {
        stop("x$trajectories is empty; merging requires stored trajectories.")
    }

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for connector path finding.")
    }

    ## ------------------------------------------------------------------------
    ## Merge specification (merge.map, or losers+winners)
    ## ------------------------------------------------------------------------

    if (is.null(merge.map)) {
        if (is.null(losers) || is.null(winners)) {
            stop("Provide merge.map, or provide both losers and winners.")
        }
        if (length(winners) == 1L && length(losers) > 1L) {
            winners <- rep(winners, length(losers))
        }
        if (length(losers) != length(winners)) {
            stop("losers and winners must have the same length (or winners length 1).")
        }
        merge.map <- winners
        names(merge.map) <- losers
    }

    if (is.null(names(merge.map)) || any(names(merge.map) == "")) {
        stop("merge.map must be a named vector: names are losers, values are winners.")
    }

    ## ------------------------------------------------------------------------
    ## Resolve labels/vertex indices to vertices for the requested extremum type
    ## ------------------------------------------------------------------------

    summary.df <- x$summary.all
    if (is.null(summary.df) || nrow(summary.df) == 0) {
        summary.df <- x$summary
    }
    if (is.null(summary.df) || nrow(summary.df) == 0) {
        stop("x has no summary tables (summary.all / summary).")
    }

    summary.df.type <- summary.df[summary.df$type == type, , drop = FALSE]
    if (nrow(summary.df.type) == 0) {
        stop("No extrema of requested type='", type, "' found in x$summary(.all).")
    }

    label.to.vertex <- setNames(as.integer(summary.df.type$vertex), summary.df.type$label)
    vertex.to.label <- setNames(as.character(summary.df.type$label), as.character(summary.df.type$vertex))

    resolve.to.vertex <- function(loc) {
        if (is.numeric(loc) && length(loc) == 1L) {
            return(as.integer(loc))
        }
        if (is.character(loc) && length(loc) == 1L) {
            if (grepl("^[0-9]+$", loc)) {
                return(as.integer(loc))
            }
            if (!is.na(label.to.vertex[loc])) {
                return(as.integer(label.to.vertex[loc]))
            }
            stop("Could not resolve extremum label '", loc, "' for type='", type, "'.")
        }
        stop("Each loser/winner must be a single label or a single vertex index.")
    }

    loser.loc <- names(merge.map)
    winner.loc <- as.vector(merge.map)

    loser.vertex <- vapply(loser.loc, resolve.to.vertex, integer(1))
    winner.vertex <- vapply(winner.loc, resolve.to.vertex, integer(1))

    if (any(loser.vertex == winner.vertex)) {
        stop("Some merges map a loser to itself (loser == winner).")
    }

    ## ------------------------------------------------------------------------
    ## Collapse transitive merges (A->B, B->C => A->C), and detect cycles
    ## ------------------------------------------------------------------------

    map.vertex <- setNames(as.integer(winner.vertex), as.character(loser.vertex))

    resolve.winner.vertex <- function(v, map) {
        seen <- integer(0)
        cur <- as.integer(v)
        while (!is.na(map[as.character(cur)])) {
            if (cur %in% seen) {
                stop("Cycle detected in merge mapping involving vertex ", cur, ".")
            }
            seen <- c(seen, cur)
            cur <- as.integer(map[as.character(cur)])
        }
        cur
    }

    final.winner.vertex <- vapply(loser.vertex, resolve.winner.vertex, integer(1), map = map.vertex)
    map.vertex <- setNames(as.integer(final.winner.vertex), as.character(loser.vertex))

    ## Keep unique losers after collapsing
    keep.idx <- !duplicated(names(map.vertex))
    map.vertex <- map.vertex[keep.idx]

    loser.vertex <- as.integer(names(map.vertex))
    winner.vertex <- as.integer(unname(map.vertex))

    loser.label <- ifelse(!is.na(vertex.to.label[as.character(loser.vertex)]),
                          vertex.to.label[as.character(loser.vertex)],
                          paste0("v", loser.vertex))
    winner.label <- ifelse(!is.na(vertex.to.label[as.character(winner.vertex)]),
                           vertex.to.label[as.character(winner.vertex)],
                           paste0("v", winner.vertex))

    if (verbose) {
        cat("merge.basins.gfc.flow(): merging ", length(loser.vertex), " ", type, " basin(s)\n", sep = "")
        for (i in seq_along(loser.vertex)) {
            cat("  ", loser.label[i], " (v", loser.vertex[i], ") -> ",
                winner.label[i], " (v", winner.vertex[i], ")\n", sep = "")
        }
    }

    ## ------------------------------------------------------------------------
    ## Build igraph and compute connector paths
    ## ------------------------------------------------------------------------

    if (!connector.method %in% c("shortest", "closest")) {
        stop("connector.method='", connector.method, "' is reserved but not implemented yet. ",
             "Use connector.method='shortest' or 'closest' for now.")
    }

    n.vertices <- x$n.vertices %||% length(adj.list)

    ## Build directed edge list, then collapse to undirected with min weight
    edge.from <- integer(0)
    edge.to <- integer(0)
    edge.weight <- double(0)

    for (i in seq_len(length(adj.list))) {
        nbrs <- adj.list[[i]]
        if (length(nbrs) == 0) next

        wts <- weight.list[[i]]
        if (is.null(wts) || length(wts) != length(nbrs)) {
            stop("weight.list[[", i, "]] must align with adj.list[[", i, "]].")
        }

        edge.from <- c(edge.from, rep.int(i, length(nbrs)))
        edge.to <- c(edge.to, as.integer(nbrs))
        edge.weight <- c(edge.weight, as.double(wts))
    }

    g.df <- data.frame(
        from = as.character(edge.from),
        to = as.character(edge.to),
        weight = as.double(edge.weight),
        stringsAsFactors = FALSE
    )

    v.df <- data.frame(name = as.character(seq_len(n.vertices)), stringsAsFactors = FALSE)

    g <- igraph::graph_from_data_frame(g.df, directed = TRUE, vertices = v.df)
    g <- igraph::as.undirected(g, mode = "collapse",
                               edge.attr.comb = list(weight = "min"))
    g <- igraph::simplify(g,
                          remove.multiple = TRUE,
                          remove.loops = TRUE,
                          edge.attr.comb = list(weight = "min"))

    e.weight <- igraph::E(g)$weight
    if (is.null(e.weight)) {
        stop("Failed to assign edge weights in igraph object.")
    }

    ## ------------------------------------------------------------------------
    ## Precompute winner basin vertex sets (needed for connector.method='closest')
    ## ------------------------------------------------------------------------
    connector.by.loser <- vector("list", length(loser.vertex))
    names(connector.by.loser) <- as.character(loser.vertex)

    if (connector.method == "closest") {

        basins.all <- if (type == "max") x$max.basins.all else x$min.basins.all
        if (is.null(basins.all) || length(basins.all) == 0L) {
            stop("connector.method='closest' requires ", type, ".basins.all in x.")
        }

        basin.vertices.by.ext <- lapply(basins.all, function(b) as.integer(b$vertices))
        names(basin.vertices.by.ext) <- vapply(basins.all, function(b) as.character(b$extremum.vertex), character(1))

        ## Prefer next.up if present (fast O(path)), fallback to scanning trajectories
        .asc.path.x.to.w <- function(x.vertex,
                                     w.vertex,
                                     trajectories,
                                     next.up = NULL,
                                     max.steps = 100000L) {

            x.vertex <- as.integer(x.vertex)
            w.vertex <- as.integer(w.vertex)

            ## ------------------------------------------------------------
            ## Primary: follow next.up if available
            ## ------------------------------------------------------------
            if (!is.null(next.up)) {
                cur <- x.vertex
                path <- cur
                seen <- integer(0)
                n.steps <- 0L

                while (cur != w.vertex) {
                    n.steps <- n.steps + 1L
                    if (n.steps > max.steps) {
                        break
                    }
                    if (cur %in% seen) {
                        break
                    }
                    seen <- c(seen, cur)

                    nxt <- next.up[cur]
                    if (is.na(nxt)) {
                        break
                    }
                    cur <- as.integer(nxt)
                    path <- c(path, cur)
                }

                if (length(path) >= 1L && tail(path, 1L) == w.vertex) {
                    return(path)
                }
                ## else: fall through to trajectory-suffix fallback
            }

            ## ------------------------------------------------------------
            ## Fallback: find a stored trajectory ending at w.vertex that contains x.vertex
            ## ------------------------------------------------------------
            for (tr in trajectories) {
                if (as.integer(tr$end.vertex) != w.vertex) next
                vv <- as.integer(tr$vertices)
                pos <- match(x.vertex, vv)
                if (!is.na(pos)) {
                    return(vv[pos:length(vv)])
                }
            }

            stop("Could not recover a path from x=", x.vertex, " to winner=", w.vertex,
                 " via next.up nor via stored trajectories.")
        }
    }

    ## ------------------------------------------------------------------------
    ## Compute connector paths
    ## ------------------------------------------------------------------------

    closest.x.by.loser <- list()

    for (i in seq_along(loser.vertex)) {

        l.v <- as.integer(loser.vertex[i])
        w.v <- as.integer(winner.vertex[i])

        from.v <- as.character(l.v)
        to.v <- as.character(w.v)

        if (connector.method == "closest") {

            w.basin <- basin.vertices.by.ext[[as.character(w.v)]]
            if (is.null(w.basin) || length(w.basin) == 0L) {
                stop("Winner vertex ", w.v, " has empty or missing basin vertices.")
            }

            ## Dijkstra distances from loser to all vertices in winner basin
            d <- igraph::distances(
                             g,
                             v = from.v,
                             to = as.character(w.basin),
                             weights = igraph::E(g)$weight,
                             mode = "all"
                         )

            d <- as.numeric(d[1, ])
            if (all(!is.finite(d))) {
                stop("No finite path from loser vertex ", l.v, " to any vertex in winner basin of ", w.v, ".")
            }

            ## Choose closest basin vertex x (deterministic tie-break by vertex id)
            ok <- which(is.finite(d))
            w.basin.ok <- w.basin[ok]
            d.ok <- d[ok]
            ord <- order(d.ok, w.basin.ok)
            x.v <- as.integer(w.basin.ok[ord[1]])

            closest.x.by.loser[[as.character(l.v)]] <- as.integer(x.v)

            ## P(lM, x): weighted shortest path lM -> x
            sp.lx <- igraph::shortest_paths(
                                 g,
                                 from = from.v,
                                 to = as.character(x.v),
                                 weights = igraph::E(g)$weight,
                                 output = "vpath"
                             )

            vpath.lx <- sp.lx$vpath[[1]]
            if (length(vpath.lx) == 0) {
                stop("No path found between loser vertex ", l.v, " and closest basin vertex x=", x.v, ".")
            }

            path.lx <- as.integer(igraph::as_ids(vpath.lx))

            ## P(x, wM): follow ascent map inside winner basin (next.up preferred), otherwise a stored trajectory suffix
            next.up <- x$next.up
            if (!is.null(next.up)) {
                next.up <- as.integer(next.up)  ## already 1-based in R, NA for maxima
            }

            path.xw <- .asc.path.x.to.w(
                x.vertex = x.v,
                w.vertex = w.v,
                trajectories = x$trajectories,
                next.up = next.up
            )

            ## Combine: [lM ... x] + [x ... wM] without duplicating x
            if (length(path.xw) < 2L) {
                stop("Internal error: path.xw too short for x=", x.v, " -> wM=", w.v, ".")
            }

            connector.by.loser[[as.character(l.v)]] <- c(path.lx, path.xw[-1L])

        } else if (connector.method == "shortest") {

            sp <- igraph::shortest_paths(
                              g,
                              from = from.v,
                              to = to.v,
                              weights = igraph::E(g)$weight,
                              output = "vpath"
                          )

            vpath <- sp$vpath[[1]]
            if (length(vpath) == 0) {
                stop("No path found between loser vertex ", l.v,
                     " and winner vertex ", w.v, ".")
            }

            connector.by.loser[[as.character(l.v)]] <- as.integer(igraph::as_ids(vpath))

        }
    }

    ## ------------------------------------------------------------------------
    ## Modify trajectories by appending connectors
    ## ------------------------------------------------------------------------

    map.winner.by.loser <- setNames(as.integer(winner.vertex), as.character(loser.vertex))

    n.modified <- 0L
    for (k in seq_along(x$trajectories)) {
        traj <- x$trajectories[[k]]
        l.v <- as.integer(traj$end.vertex)

        w.v <- map.winner.by.loser[as.character(l.v)]
        if (is.na(w.v)) next

        connector <- connector.by.loser[[as.character(l.v)]]
        if (is.null(connector) || length(connector) < 2L) {
            stop("Internal error: invalid connector for loser vertex ", l.v, ".")
        }

        ## Append connector excluding the first vertex (loser extremum)
        traj$vertices <- c(as.integer(traj$vertices), as.integer(connector[-1L]))
        traj$end.vertex <- as.integer(w.v)

        ## ----------------------------------------------------------------
        ## Attach traj$y.hat.modified (monotone repair along the merged path)
        ## ----------------------------------------------------------------

        y.hat.global <- x$y
        if (is.null(y.hat.global)) {
            stop("merge.basins.gfc.flow(): cannot compute traj$y.hat.modified because x$y.hat (or x$y) is missing.")
        }

        ## Ensure a clean state if this trajectory was previously processed
        traj$y.hat.modified <- NULL

        if (connector.method == "closest") {

            ## Original trajectory length before appending the connector
            orig.len <- length(traj$vertices) - (length(connector) - 1L)

            ## Closest basin vertex x for this loser
            x.v <- closest.x.by.loser[[as.character(l.v)]]
            if (is.null(x.v) || length(x.v) != 1L || is.na(x.v)) {
                stop("merge.basins.gfc.flow(): missing closest basin vertex x for loser vertex ", l.v, ".")
            }

            ## Position of x within the connector (connector includes lM as first element)
            pos.x.in.connector <- match(as.integer(x.v), as.integer(connector))
            if (is.na(pos.x.in.connector)) {
                stop("merge.basins.gfc.flow(): internal error: x not found in connector for loser vertex ", l.v, ".")
            }

            ## Position of x in the new trajectory after appending connector[-1]
            ## - original trajectory ends at lM (connector[1])
            ## - appended part is connector[2..]
            idx.x <- orig.len + (pos.x.in.connector - 1L)

            new.vertices <- as.integer(traj$vertices)
            if (idx.x < 1L || idx.x > length(new.vertices) || new.vertices[idx.x] != as.integer(x.v)) {
                stop("merge.basins.gfc.flow(): internal error locating x within the modified trajectory for loser vertex ", l.v, ".")
            }

            y.path <- as.double(y.hat.global[new.vertices])

            ## Find x' on the ORIGINAL trajectory T (indices 1..orig.len):
            ## last vertex with y.hat(x') < y.hat(x)
            cand <- which(y.path[seq_len(orig.len)] < y.path[idx.x])

            if (length(cand) > 0L) {
                idx.xprime <- max(cand)

                ## Linear interpolation (harmonic extension on a path)
                seg.idx <- idx.xprime:idx.x
                seg.len <- length(seg.idx)

                if (seg.len >= 2L) {
                    y0 <- y.path[idx.xprime]
                    y1 <- y.path[idx.x]
                    t <- (0:(seg.len - 1L)) / (seg.len - 1L)
                    y.path[seg.idx] <- y0 + t * (y1 - y0)
                }
                ## seg.len == 1L implies no change

                traj$y.hat.modified <- y.path
            } else {
                ## No eligible x' found under the strict rule; leave NULL
                ## (You may optionally relax to <= if you want more repairs.)
                traj$y.hat.modified <- NULL
            }
        }

        x$trajectories[[k]] <- traj
        n.modified <- n.modified + 1L
    }

    if (verbose) {
        cat("  Modified trajectories: ", n.modified, "\n", sep = "")
    }

    ## ------------------------------------------------------------------------
    ## Update summaries: mark losers merged (optional), and update basin.size later
    ## ------------------------------------------------------------------------

    if (type == "max" && !is.null(x$max.summaries.all)) {
        max.sum <- x$max.summaries.all

        for (i in seq_along(loser.vertex)) {
            idx.l <- match(loser.vertex[i], max.sum$vertex)
            idx.w <- match(winner.vertex[i], max.sum$vertex)
            if (!is.na(idx.l)) {
                if (mark.loser.spurious) {
                    max.sum$is.spurious[idx.l] <- TRUE
                    if ("filter.stage" %in% names(max.sum)) {
                        max.sum$filter.stage[idx.l] <- merge.filter.stage
                    }
                }
                if ("merged.into" %in% names(max.sum)) {
                    max.sum$merged.into[idx.l] <- if (!is.na(idx.w)) {
                        as.character(max.sum$label[idx.w])
                    } else {
                        winner.label[i]
                    }
                }
            }
        }

        x$max.summaries.all <- max.sum
    }

    if (type == "min" && !is.null(x$min.summaries.all)) {
        min.sum <- x$min.summaries.all

        for (i in seq_along(loser.vertex)) {
            idx.l <- match(loser.vertex[i], min.sum$vertex)
            idx.w <- match(winner.vertex[i], min.sum$vertex)
            if (!is.na(idx.l)) {
                if (mark.loser.spurious) {
                    min.sum$is.spurious[idx.l] <- TRUE
                    if ("filter.stage" %in% names(min.sum)) {
                        min.sum$filter.stage[idx.l] <- merge.filter.stage
                    }
                }
                if ("merged.into" %in% names(min.sum)) {
                    min.sum$merged.into[idx.l] <- if (!is.na(idx.w)) {
                        as.character(min.sum$label[idx.w])
                    } else {
                        winner.label[i]
                    }
                }
            }
        }

        x$min.summaries.all <- min.sum
    }

    ## ------------------------------------------------------------------------
    ## Rebuild membership lists from trajectories
    ## ------------------------------------------------------------------------

    rebuild.membership <- function(trajs, extrema.vertices, which.end = c("start", "end")) {
        which.end <- match.arg(which.end)

        mem <- vector("list", length(extrema.vertices))
        names(mem) <- as.character(extrema.vertices)
        for (i in seq_along(mem)) mem[[i]] <- integer(0)

        for (traj in trajs) {
            ext.v <- if (which.end == "start") as.integer(traj$start.vertex) else as.integer(traj$end.vertex)
            key <- as.character(ext.v)
            if (is.na(key) || is.null(mem[[key]])) next
            mem[[key]] <- c(mem[[key]], as.integer(traj$vertices))
        }

        mem <- lapply(mem, function(v) unique(as.integer(v)))
        mem
    }

    ## Use existing basin lists as the authoritative set of extrema vertices (all)
    max.vertices.all <- if (!is.null(x$max.basins.all) && length(x$max.basins.all) > 0) {
        vapply(x$max.basins.all, function(b) as.integer(b$extremum.vertex), integer(1))
    } else {
        summary.df$vertex[summary.df$type == "max"]
    }

    min.vertices.all <- if (!is.null(x$min.basins.all) && length(x$min.basins.all) > 0) {
        vapply(x$min.basins.all, function(b) as.integer(b$extremum.vertex), integer(1))
    } else {
        summary.df$vertex[summary.df$type == "min"]
    }

    max.mem.all <- rebuild.membership(x$trajectories, max.vertices.all, which.end = "end")
    min.mem.all <- rebuild.membership(x$trajectories, min.vertices.all, which.end = "start")

    x$max.membership.all <- max.mem.all
    x$min.membership.all <- min.mem.all

    ## ------------------------------------------------------------------------
    ## Update basin vertices + hop distances in basins.all
    ## ------------------------------------------------------------------------

    compute.hop.distances <- function(ext.v, verts) {
        if (length(verts) == 0) {
            return(list(hop.distances = integer(0), max.hop.distance = 0L))
        }

        d.mat <- igraph::distances(
            g,
            v = as.character(ext.v),
            to = as.character(verts),
            weights = NULL,
            mode = "all"
        )

        hop <- as.integer(d.mat[1, ])
        if (any(!is.finite(hop))) {
            ## Should not happen in typical connected basins, but guard anyway
            hop[!is.finite(hop)] <- NA_integer_
        }

        list(
            hop.distances = hop,
            max.hop.distance = if (length(hop) > 0) max(hop, na.rm = TRUE) else 0L
        )
    }

    update.basins.all <- function(basins.all, mem.all, type) {
        if (is.null(basins.all) || length(basins.all) == 0) return(basins.all)

        for (i in seq_along(basins.all)) {
            b <- basins.all[[i]]
            ext.v <- as.integer(b$extremum.vertex)

            verts <- mem.all[[as.character(ext.v)]]
            if (is.null(verts)) verts <- integer(0)

            b$vertices <- as.integer(verts)

            hop.info <- compute.hop.distances(ext.v, verts)
            b$hop.distances <- hop.info$hop.distances
            b$max.hop.distance <- as.integer(hop.info$max.hop.distance)

            ## If requested, mark merged losers as spurious in basin records too
            if (mark.loser.spurious && ext.v %in% loser.vertex && type == type) {
                b$is.spurious <- TRUE
                if (!is.null(b$filter.stage)) b$filter.stage <- merge.filter.stage
                if (!is.null(b$merged.into)) {
                    j <- match(ext.v, loser.vertex)
                    b$merged.into <- winner.label[j]
                }
            }

            basins.all[[i]] <- b
        }

        basins.all
    }

    x$max.basins.all <- update.basins.all(x$max.basins.all, x$max.membership.all, "max")
    x$min.basins.all <- update.basins.all(x$min.basins.all, x$min.membership.all, "min")

    ## ------------------------------------------------------------------------
    ## Update basin.size in min/max summaries (ALL), then rebuild summary tables
    ## ------------------------------------------------------------------------

    update.basin.size <- function(sum.df, mem.all) {
        if (is.null(sum.df) || nrow(sum.df) == 0) return(sum.df)
        if (!("basin.size" %in% names(sum.df))) return(sum.df)

        bs <- integer(nrow(sum.df))
        for (i in seq_len(nrow(sum.df))) {
            v <- as.integer(sum.df$vertex[i])
            vv <- mem.all[[as.character(v)]]
            bs[i] <- if (is.null(vv)) 0L else length(vv)
        }
        sum.df$basin.size <- as.integer(bs)
        sum.df
    }

    if (!is.null(x$max.summaries.all)) {
        x$max.summaries.all <- update.basin.size(x$max.summaries.all, x$max.membership.all)
    }
    if (!is.null(x$min.summaries.all)) {
        x$min.summaries.all <- update.basin.size(x$min.summaries.all, x$min.membership.all)
    }

    ## Recompute retained/spurious indices and counts from basins.all flags
    recompute.indices <- function(basins.all) {
        if (is.null(basins.all) || length(basins.all) == 0) {
            return(list(retained = integer(0), spurious = integer(0)))
        }
        is.spur <- vapply(basins.all, function(b) isTRUE(b$is.spurious), logical(1))
        list(
            retained = which(!is.spur),
            spurious = which(is.spur)
        )
    }

    idx.max <- recompute.indices(x$max.basins.all)
    idx.min <- recompute.indices(x$min.basins.all)

    x$retained.max.indices <- as.integer(idx.max$retained)
    x$spurious.max.indices <- as.integer(idx.max$spurious)
    x$retained.min.indices <- as.integer(idx.min$retained)
    x$spurious.min.indices <- as.integer(idx.min$spurious)

    x$n.max.retained <- length(x$retained.max.indices)
    x$n.min.retained <- length(x$retained.min.indices)
    x$n.max.spurious <- length(x$spurious.max.indices)
    x$n.min.spurious <- length(x$spurious.min.indices)

    ## Rebuild summary/all summary using existing internal postprocessor
    x <- .postprocess.gfc.flow(x)

    ## ------------------------------------------------------------------------
    ## Attach merge report
    ## ------------------------------------------------------------------------

    merge.report <- data.frame(
        type = type,
        loser.label = loser.label,
        winner.label = winner.label,
        loser.vertex = loser.vertex,
        winner.vertex = winner.vertex,
        connector.length = vapply(connector.by.loser[as.character(loser.vertex)], length, integer(1)),
        stringsAsFactors = FALSE
    )

    x$merge.report <- merge.report

    class(x) <- unique(c("gfc.flow", class(x)))
    x
}

#' Get Cell Vertices for a (min, max) Pair
#'
#' @description
#' Returns the set of domain-graph vertices belonging to the $\hat{y}$-cell
#' determined by a minimum extremum and a maximum extremum. In a trajectory-first
#' \code{gfc.flow} object, a cell is represented by the union of vertices from all
#' trajectories that start at the given minimum and end at the given maximum.
#'
#' @param gfc.flow A \code{gfc.flow} object from \code{compute.gfc.flow()} with
#'   trajectories computed (i.e., \code{store.trajectories = TRUE}).
#' @param min.id Minimum extremum identifier: either a label (e.g., \code{"m4"},
#'   \code{"sm1"}) or a vertex index (integer, 1-based).
#' @param max.id Maximum extremum identifier: either a label (e.g., \code{"M1"},
#'   \code{"sM2"}) or a vertex index (integer, 1-based).
#'
#' @return An integer vector of unique vertex indices in the requested cell.
#'
#' @export
cell.vertices <- function(gfc.flow, min.id, max.id) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object.")
    }
    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        stop("gfc.flow has no trajectories. Recompute with store.trajectories = TRUE.")
    }

    ## ------------------------------------------------------------------------
    ## Resolve labels / indices to vertices
    ## ------------------------------------------------------------------------

    summary.df <- gfc.flow$summary.all
    if (is.null(summary.df) || nrow(summary.df) == 0) {
        summary.df <- gfc.flow$summary
    }
    if (is.null(summary.df) || nrow(summary.df) == 0) {
        stop("gfc.flow does not contain summary tables (summary.all / summary).")
    }

    min.df <- summary.df[summary.df$type == "min", , drop = FALSE]
    max.df <- summary.df[summary.df$type == "max", , drop = FALSE]

    if (nrow(min.df) == 0) stop("No minima found in gfc.flow summary tables.")
    if (nrow(max.df) == 0) stop("No maxima found in gfc.flow summary tables.")

    min.label.to.vertex <- setNames(as.integer(min.df$vertex), as.character(min.df$label))
    max.label.to.vertex <- setNames(as.integer(max.df$vertex), as.character(max.df$label))

    resolve.extremum.vertex <- function(id, label.to.vertex, what) {

        if (length(id) != 1L) {
            stop(what, ".id must be a single label or a single vertex index.")
        }

        if (is.numeric(id)) {
            v <- as.integer(id)
            if (is.na(v) || v < 1L) stop("Invalid ", what, " vertex index.")
            return(v)
        }

        if (!is.character(id)) {
            stop(what, ".id must be a character label or a numeric vertex index.")
        }

        if (grepl("^[0-9]+$", id)) {
            v <- as.integer(id)
            if (is.na(v) || v < 1L) stop("Invalid ", what, " vertex index string.")
            return(v)
        }

        v <- label.to.vertex[id]
        if (is.na(v)) {
            stop("Could not resolve ", what, " label '", id, "'.")
        }

        as.integer(v)
    }

    min.vertex <- resolve.extremum.vertex(min.id, min.label.to.vertex, "min")
    max.vertex <- resolve.extremum.vertex(max.id, max.label.to.vertex, "max")

    ## ------------------------------------------------------------------------
    ## Collect cell vertices from trajectories
    ## ------------------------------------------------------------------------

    verts <- integer(0)

    for (traj in gfc.flow$trajectories) {

        sv <- as.integer(traj$start.vertex)
        ev <- as.integer(traj$end.vertex)

        if (is.na(sv) || is.na(ev)) next

        if (sv == min.vertex && ev == max.vertex) {
            verts <- c(verts, as.integer(traj$vertices))
        }
    }

    unique(as.integer(verts))
}

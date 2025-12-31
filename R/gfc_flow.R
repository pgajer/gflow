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

    if (!basin$is.spurious) {
        warning(sprintf("Basin '%s' is not spurious", label))
    }

    basin.vertices <- basin$vertices
    basin.set <- as.integer(basin.vertices)

    ## Find boundary: vertices in basin with neighbors outside basin
    boundary.vertices <- c()
    for (v in basin.vertices) {
        nbrs <- adj.list[[v]]
        if (any(!(nbrs %in% basin.set))) {
            boundary.vertices <- c(boundary.vertices, v)
        }
    }

    ## Interior = basin - boundary
    interior.vertices <- setdiff(basin.vertices, boundary.vertices)

    ## Get boundary values
    boundary.values <- y[boundary.vertices]
    names(boundary.values) <- boundary.vertices

    list(
        interior.vertices = interior.vertices,
        boundary.vertices = boundary.vertices,
        boundary.values = boundary.values,
        extremum.vertex = basin$extremum.vertex,
        label = label
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

## Helper: null-coalescing operator
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
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
#' @param x A gfc_cell_trajectories object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.gfc_cell_trajectories <- function(x, ...) {
    cat("GFC Cell Trajectories\n")
    cat("=====================\n")
    ## Show spurious status if available
    min.status <- if (!is.null(x$min.is.spurious) && x$min.is.spurious) " [SPURIOUS]" else ""
    max.status <- if (!is.null(x$max.is.spurious) && x$max.is.spurious) " [SPURIOUS]" else ""
    cat(sprintf("Cell: %s%s (vertex %d) <-> %s%s (vertex %d)\n",
                x$min.label, min.status,
                if (x$mapped) x$original.min.vertex else x$min.vertex,
                x$max.label, max.status,
                if (x$mapped) x$original.max.vertex else x$max.vertex))
    cat(sprintf("Values: %.4f (min) to %.4f (max), delta = %.4f\n",
                x$min.value, x$max.value, x$max.value - x$min.value))
    cat(sprintf("Trajectories: %d\n", x$n.trajectories))
    if (x$mapped) {
        cat(sprintf("Vertices mapped to subgraph indices (min -> %s, max -> %s)\n",
                    ifelse(is.na(x$min.vertex), "NA", as.character(x$min.vertex)),
                    ifelse(is.na(x$max.vertex), "NA", as.character(x$max.vertex))))
    }
    if (x$n.trajectories > 0) {
        lengths <- sapply(x$trajectories, length)
        cat(sprintf("Trajectory lengths: Min: %d, Max: %d, Mean: %.1f\n",
                    min(lengths), max(lengths), mean(lengths)))
        n.unique.vertices <- length(unique(unlist(x$trajectories)))
        cat(sprintf("Number of vertices: %d\n", n.unique.vertices))
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


#' List All Cells in GFC Flow Result
#'
#' Enumerates all cells (min-max pairs) in a gfc.flow result, including cells
#' with spurious endpoints.
#'
#' @param gfc.flow A gfc.flow object from \code{compute.gfc.flow()}.
#' @param include.spurious Logical; if TRUE (default), include cells with
#'   spurious endpoints.
#'
#' @return A data frame with columns:
#'   \item{min.label}{Minimum label.}
#'   \item{max.label}{Maximum label.}
#'   \item{min.vertex}{Minimum vertex index.}
#'   \item{max.vertex}{Maximum vertex index.}
#'   \item{min.is.spurious}{TRUE if minimum is spurious.}
#'   \item{max.is.spurious}{TRUE if maximum is spurious.}
#'   \item{n.trajectories}{Number of trajectories connecting this pair.}
#'   \item{unique.vertices}{Number of unique vertices in trajectories.}
#'
#' @export
list.cells.gfc.flow <- function(gfc.flow, include.spurious = TRUE) {

    if (!inherits(gfc.flow, "gfc.flow")) {
        stop("gfc.flow must be a gfc.flow object from compute.gfc.flow()")
    }

    if (is.null(gfc.flow$trajectories) || length(gfc.flow$trajectories) == 0) {
        return(data.frame(
            min.label = character(0),
            max.label = character(0),
            min.vertex = integer(0),
            max.vertex = integer(0),
            min.is.spurious = logical(0),
            max.is.spurious = logical(0),
            n.trajectories = integer(0),
            unique.vertices = integer(0)
        ))
    }

    summary.df <- gfc.flow$summary.all %||% gfc.flow$summary

    ## Collect cell information
    cells <- list()

    for (traj in gfc.flow$trajectories) {
        min.v <- traj$start.vertex
        max.v <- traj$end.vertex
        key <- paste(min.v, max.v, sep = "-")

        if (is.null(cells[[key]])) {
            cells[[key]] <- list(
                min.vertex = min.v,
                max.vertex = max.v,
                start.is.spurious = traj$start.is.spurious %||% FALSE,
                end.is.spurious = traj$end.is.spurious %||% FALSE,
                n.trajectories = 0,
                all.vertices = c()
            )
        }

        cells[[key]]$n.trajectories <- cells[[key]]$n.trajectories + 1
        cells[[key]]$all.vertices <- c(cells[[key]]$all.vertices, traj$vertices)
    }

    ## Build result data frame
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
        result$min.vertex[i] <- cell$min.vertex
        result$max.vertex[i] <- cell$max.vertex
        result$min.is.spurious[i] <- cell$start.is.spurious
        result$max.is.spurious[i] <- cell$end.is.spurious
        result$n.trajectories[i] <- cell$n.trajectories
        result$unique.vertices[i] <- length(unique(cell$all.vertices))
    }

    ## Optionally filter out spurious cells
    if (!include.spurious) {
        keep <- !result$min.is.spurious & !result$max.is.spurious
        result <- result[keep, , drop = FALSE]
    }

    ## Sort by number of trajectories (descending)
    result <- result[order(-result$n.trajectories), ]
    rownames(result) <- NULL

    return(result)
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

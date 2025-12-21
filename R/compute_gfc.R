#' Compute Refined Gradient Flow Complex
#'
#' Computes the gradient flow complex (GFC) for a scalar function defined on
#' a weighted graph, applying a sequence of refinement operations to produce
#' a robust basin structure suitable for association analysis. This function
#' mirrors the behavior of \code{compute.refined.basins()} for compatibility.
#'
#' @param adj.list A list of integer vectors representing the graph's adjacency
#'   structure. Element \code{i} contains the 1-based indices of vertices
#'   adjacent to vertex \code{i}.
#' @param edge.length.list A list of numeric vectors containing edge lengths.
#'   Element \code{i} contains lengths corresponding to the edges in
#'   \code{adj.list[[i]]}. Must have the same structure as \code{adj.list}.
#' @param fitted.values Numeric vector of function values at each vertex.
#'   The length must equal the number of vertices.
#' @param edge.length.quantile.thld Numeric value in (0, 1] specifying the
#'   quantile of the graph's edge length distribution to use as a threshold
#'   for basin construction. Default is 0.9.
#' @param min.rel.value.max Minimum relative value threshold for retaining
#'   maxima. Maxima with \code{value / median(y)} below this threshold are
#'   removed. Default is 1.1.
#' @param max.rel.value.min Maximum relative value threshold for retaining
#'   minima. Minima with \code{value / median(y)} above this threshold are
#'   removed. Default is 0.9.
#' @param max.overlap.threshold Overlap distance threshold for clustering
#'   maxima. Basins with overlap distance below this threshold are merged.
#'   Default is 0.15.
#' @param min.overlap.threshold Overlap distance threshold for clustering
#'   minima. Default is 0.15.
#' @param p.mean.nbrs.dist.threshold Percentile threshold for mean neighbor
#'   distance in geometric filtering. Extrema at vertices exceeding this
#'   percentile are removed. This filter is applied symmetrically to both
#'   maxima and minima, consistent with \code{compute.refined.basins()}.
#'   Default is 0.9.
#' @param p.mean.hopk.dist.threshold Percentile threshold for mean hop-k
#'   distance in geometric filtering. Extrema exceeding this percentile are
#'   removed. Default is 0.9.
#' @param p.deg.threshold Percentile threshold for vertex degree in geometric
#'   filtering. Extrema at vertices exceeding this percentile are removed.
#'   Default is 0.9.
#' @param min.basin.size Minimum number of vertices required for a basin to
#'   be retained. Default is 10.
#' @param expand.basins Logical indicating whether to expand basins to cover
#'   all graph vertices after filtering. Uncovered vertices are assigned to
#'   the nearest retained basin. Default is TRUE.
#' @param apply.relvalue.filter Logical indicating whether to apply relative
#'   value filtering. Default is TRUE.
#' @param apply.maxima.clustering Logical indicating whether to cluster and
#'   merge similar maximum basins. Default is TRUE.
#' @param apply.minima.clustering Logical indicating whether to cluster and
#'   merge similar minimum basins. Default is TRUE.
#' @param apply.geometric.filter Logical indicating whether to apply geometric
#'   filtering based on hop distance and degree. Default is TRUE.
#' @param hop.k Parameter for hop-k distance calculation in summary statistics.
#'   Default is 2.
#' @param with.trajectories Logical indicating whether to return gradient
#'   trajectories for each basin. When TRUE, the result includes complete
#'   trajectory information enabling reconstruction of all gradient flow paths.
#'   This is computationally more expensive and produces larger output objects.
#'   Default is FALSE.
#' @param max.chain.depth Integer specifying the maximum recursion depth when
#'   extending trajectories through spurious extrema. Higher values allow
#'   trajectories to be joined across more intermediate spurious extrema.
#'   Default is 5.
#' @param verbose Logical indicating whether to print progress messages.
#'   Default is FALSE.
#'
#' @return A list of class \code{"gfc"} containing:
#'   \describe{
#'     \item{max.basins}{Named list of maximum basins. Each basin contains
#'       \code{vertex} (1-based extremum location), \code{value} (function
#'       value at extremum), \code{vertices} (1-based basin member indices),
#'       \code{hop.distances} (hop distance from extremum for each vertex),
#'       and \code{max.hop.distance}.}
#'     \item{min.basins}{Named list of minimum basins with the same structure.}
#'     \item{summary}{Data frame with one row per retained extremum containing
#'       label, vertex, value, rel.value, type, hop.idx, basin.size, and
#'       geometric measures. Minima are labeled m1, m2, ... in order of
#'       increasing value; maxima are labeled M1, M2, ... in order of
#'       decreasing value.}
#'     \item{max.membership}{List of length n.vertices where element v contains
#'       the 1-based indices of maximum basins containing vertex v.}
#'     \item{min.membership}{List of length n.vertices for minimum basins.}
#'     \item{expanded.max.assignment}{Integer vector of length n.vertices
#'       giving the 1-based index of the assigned maximum basin for each
#'       vertex (NA if unassigned). Only populated if expand.basins = TRUE.}
#'     \item{expanded.min.assignment}{Integer vector for minimum basins.}
#'     \item{stage.history}{Data frame recording the number of extrema before
#'       and after each refinement stage.}
#'     \item{n.vertices}{Total number of vertices in the graph.}
#'     \item{y.median}{Median of the function values.}
#'     \item{parameters}{Named list of all parameter values used.}
#'     \item{joined.trajectories}{(If \code{with.trajectories = TRUE}) List of
#'       joined trajectories, each containing \code{min.vertex}, \code{max.vertex},
#'       \code{path}, \code{intermediate.extrema}, \code{total.change}, and
#'       \code{path.length}. These are complete gradient trajectories from
#'       non-spurious minima to non-spurious maxima.}
#'     \item{cell.map}{(If \code{with.trajectories = TRUE}) Named list mapping
#'       "minVertex-maxVertex" pairs to indices in \code{joined.trajectories}.
#'       This organizes trajectories by their Morse-Smale cell.}
#'   }
#'
#' @details
#' The GFC computation proceeds through several stages, each of which can be
#' enabled or disabled via parameters.
#'
#' The initial computation uses \code{compute_geodesic_basin()} to perform
#' monotone BFS with edge length filtering. Starting from each local extremum,
#' the algorithm expands to neighbors with monotonically changing function
#' values (descending for maxima, ascending for minima). Crucially, edges
#' longer than the specified quantile threshold (\code{edge.length.quantile.thld})
#' are skipped during BFS, preventing "basin jumping" through long edges. This
#' matches the behavior of R's \code{compute.basins.of.attraction()} function.
#'
#' The relative value filtering stage removes extrema whose values are too
#' close to the median. This focuses the analysis on prominent features of
#' the function landscape rather than minor fluctuations.
#'
#' The clustering stage identifies groups of extrema whose basins exhibit
#' substantial overlap, indicating that they likely represent the same
#' underlying feature. The overlap coefficient (Szymkiewicz-Simpson index)
#' measures similarity through the ratio of intersection size to minimum
#' basin size. Connected components of the threshold graph define clusters,
#' and basins within each cluster are merged.
#'
#' The geometric filtering stage removes extrema whose basins exhibit unusual
#' structural characteristics. All geometric filters (mean neighbor distance,
#' mean hop-k distance, and degree percentile) are applied symmetrically to
#' both maxima and minima, consistent with \code{compute.refined.basins()}.
#'
#' When \code{with.trajectories = TRUE}, the function additionally computes
#' joined trajectories that connect non-spurious minima to non-spurious maxima.
#' Trajectories that terminate at spurious extrema are extended by recursively
#' exploring the basins of those spurious extrema to find paths to non-spurious
#' targets. The \code{max.chain.depth} parameter limits the recursion depth.
#'
#' @examples
#' \dontrun{
#' ## Compute GFC with default parameters
#' gfc <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                    verbose = TRUE)
#'
#' ## Examine the result
#' print(gfc)
#' print(gfc$summary)
#'
#' ## Use parameters matching compute.refined.basins() defaults
#' gfc <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                    p.mean.nbrs.dist.threshold = 0.9,
#'                    p.mean.hopk.dist.threshold = 0.9,
#'                    p.deg.threshold = 0.9,
#'                    verbose = TRUE)
#'
#' ## With trajectories for downstream analysis
#' gfc <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                    with.trajectories = TRUE,
#'                    verbose = TRUE)
#'
#' ## Access joined trajectories
#' length(gfc$joined.trajectories)  # Total number of trajectories
#' names(gfc$cell.map)              # Cell identifiers
#' }
#'
#' @seealso \code{\link{compute.gfc.matrix}} for computing GFC for multiple
#'   functions, \code{\link{compute.refined.basins}} for the R implementation
#'
#' @export
compute.gfc <- function(adj.list,
                        edge.length.list,
                        fitted.values,
                        edge.length.quantile.thld = 0.9,
                        min.rel.value.max = 1.1,
                        max.rel.value.min = 0.9,
                        max.overlap.threshold = 0.15,
                        min.overlap.threshold = 0.15,
                        p.mean.nbrs.dist.threshold = 0.9,
                        p.mean.hopk.dist.threshold = 0.9,
                        p.deg.threshold = 0.9,
                        min.basin.size = 10L,
                        expand.basins = TRUE,
                        apply.relvalue.filter = TRUE,
                        apply.maxima.clustering = TRUE,
                        apply.minima.clustering = TRUE,
                        apply.geometric.filter = TRUE,
                        hop.k = 2L,
                        with.trajectories = FALSE,
                        max.chain.depth = 5L,
                        verbose = FALSE) {

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(edge.length.list)) {
        stop("edge.length.list must be a list")
    }
    if (!is.numeric(fitted.values)) {
        stop("fitted.values must be numeric")
    }
    if (length(adj.list) != length(edge.length.list)) {
        stop("adj.list and edge.length.list must have the same length")
    }
    if (length(fitted.values) != length(adj.list)) {
        stop("fitted.values must have the same length as adj.list")
    }

    ## Validate numeric parameters
    if (!is.numeric(edge.length.quantile.thld) ||
        edge.length.quantile.thld <= 0 ||
        edge.length.quantile.thld > 1) {
        stop("edge.length.quantile.thld must be in (0, 1]")
    }

    if (!is.numeric(max.overlap.threshold) ||
        max.overlap.threshold < 0 ||
        max.overlap.threshold > 1) {
        stop("max.overlap.threshold must be in [0, 1]")
    }

    if (!is.numeric(min.overlap.threshold) ||
        min.overlap.threshold < 0 ||
        min.overlap.threshold > 1) {
        stop("min.overlap.threshold must be in [0, 1]")
    }

    ## Validate trajectory parameters
    if (!is.logical(with.trajectories) || length(with.trajectories) != 1) {
        stop("with.trajectories must be a single logical value")
    }

    if (!is.numeric(max.chain.depth) || max.chain.depth < 0) {
        stop("max.chain.depth must be a non-negative integer")
    }

    ## ========================================================================
    ## Convert to 0-based indexing for C++
    ## ========================================================================

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))

    ## ========================================================================
    ## Call C++ implementation
    ## ========================================================================

    result <- .Call(
        S_compute_gfc,
        adj.list.0based,
        edge.length.list,
        as.double(fitted.values),
        as.double(edge.length.quantile.thld),
        as.double(min.rel.value.max),
        as.double(max.rel.value.min),
        as.double(max.overlap.threshold),
        as.double(min.overlap.threshold),
        as.double(p.mean.nbrs.dist.threshold),
        as.double(p.mean.hopk.dist.threshold),
        as.double(p.deg.threshold),
        as.integer(min.basin.size),
        as.logical(expand.basins),
        as.logical(apply.relvalue.filter),
        as.logical(apply.maxima.clustering),
        as.logical(apply.minima.clustering),
        as.logical(apply.geometric.filter),
        as.integer(hop.k),
        as.logical(with.trajectories),
        as.integer(max.chain.depth),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## ========================================================================
    ## Add names to basins based on summary labels (match by vertex)
    ## ========================================================================

    if (nrow(result$summary) > 0) {
        ## Build vertex-to-label lookup from summary
        vertex.to.label <- result$summary$label
        names(vertex.to.label) <- as.character(result$summary$vertex)

        ## Name max.basins by matching vertices
        max.basin.names <- sapply(result$max.basins, function(b) {
            vertex.to.label[as.character(b$vertex)]
        })
        names(result$max.basins) <- max.basin.names

        ## Name min.basins by matching vertices
        min.basin.names <- sapply(result$min.basins, function(b) {
            vertex.to.label[as.character(b$vertex)]
        })
        names(result$min.basins) <- min.basin.names
    }

    ## ========================================================================
    ## Store parameters for reproducibility
    ## ========================================================================

    result$parameters <- list(
        edge.length.quantile.thld = edge.length.quantile.thld,
        min.rel.value.max = min.rel.value.max,
        max.rel.value.min = max.rel.value.min,
        max.overlap.threshold = max.overlap.threshold,
        min.overlap.threshold = min.overlap.threshold,
        p.mean.nbrs.dist.threshold = p.mean.nbrs.dist.threshold,
        p.mean.hopk.dist.threshold = p.mean.hopk.dist.threshold,
        p.deg.threshold = p.deg.threshold,
        min.basin.size = min.basin.size,
        expand.basins = expand.basins,
        apply.relvalue.filter = apply.relvalue.filter,
        apply.maxima.clustering = apply.maxima.clustering,
        apply.minima.clustering = apply.minima.clustering,
        apply.geometric.filter = apply.geometric.filter,
        hop.k = hop.k,
        with.trajectories = with.trajectories,
        max.chain.depth = max.chain.depth
    )

    ## ========================================================================
    ## Set class and return
    ## ========================================================================

    class(result) <- c("gfc", "list")

    return(result)
}


#' Print Method for Gradient Flow Complex
#'
#' @param x A gfc object from compute.gfc()
#' @param ... Additional arguments (unused)
#'
#' @export
print.gfc <- function(x, ...) {
    cat("Gradient Flow Complex (GFC)\n")
    cat("===========================\n\n")

    n.max <- length(x$max.basins)
    n.min <- length(x$min.basins)

    cat(sprintf("Vertices: %d\n", x$n.vertices))
    cat(sprintf("Maxima: %d\n", n.max))
    cat(sprintf("Minima: %d\n", n.min))
    cat(sprintf("Median y: %.4f\n\n", x$y.median))

    ## Report coverage
    if (!is.null(x$expanded.max.assignment)) {
        n.covered <- sum(!is.na(x$expanded.max.assignment))
        cat(sprintf("Max basin coverage: %d/%d (%.1f%%)\n",
                    n.covered, x$n.vertices, 100 * n.covered / x$n.vertices))
    }
    if (!is.null(x$expanded.min.assignment)) {
        n.covered <- sum(!is.na(x$expanded.min.assignment))
        cat(sprintf("Min basin coverage: %d/%d (%.1f%%)\n",
                    n.covered, x$n.vertices, 100 * n.covered / x$n.vertices))
    }

    ## Show summary (first few rows)
    if (nrow(x$summary) > 0) {
        cat("\nBasin summary:\n")
        print(head(x$summary, 10))
        if (nrow(x$summary) > 10) {
            cat(sprintf("\n... and %d more rows\n", nrow(x$summary) - 10))
        }
    }

    invisible(x)
}


#' Compute GFC for Multiple Functions
#'
#' Efficiently computes the gradient flow complex for each column of a matrix
#' of function values over the same graph structure. Supports OpenMP
#' parallelization for improved performance with many functions.
#'
#' @param adj.list Adjacency list (1-based R indexing).
#' @param edge.length.list Edge length list.
#' @param Y Matrix of function values with n vertices (rows) and p functions
#'   (columns).
#' @param ... Additional parameters passed to \code{\link{compute.gfc}}.
#' @param n.cores Number of OpenMP threads for parallel processing. Default
#'   is 1 (sequential).
#' @param verbose Logical indicating whether to print progress. Default is
#'   FALSE.
#'
#' @return A list of length p containing GFC results for each function. If Y
#'   has column names, they are used to name the result list elements.
#'
#' @examples
#' \dontrun{
#' ## Compute GFC for 100 features
#' gfc.list <- compute.gfc.matrix(adj.list, edge.length.list, Z.hat,
#'                                n.cores = 4, verbose = TRUE)
#'
#' ## Access results for specific feature
#' gfc.feature1 <- gfc.list[[1]]
#' print(gfc.feature1$summary)
#' }
#'
#' @seealso \code{\link{compute.gfc}} for single function computation
#'
#' @export
compute.gfc.matrix <- function(adj.list,
                               edge.length.list,
                               Y,
                               ...,
                               n.cores = 1L,
                               verbose = FALSE) {

    ## Validate matrix input
    if (!is.matrix(Y)) {
        if (is.numeric(Y)) {
            Y <- matrix(Y, ncol = 1)
        } else {
            stop("Y must be a numeric matrix")
        }
    }

    if (nrow(Y) != length(adj.list)) {
        stop("nrow(Y) must equal length(adj.list)")
    }

    ## Extract optional parameters
    dots <- list(...)

    ## Build parameter set with defaults
    params <- list(
        edge.length.quantile.thld = dots$edge.length.quantile.thld %||% 0.9,
        min.rel.value.max = dots$min.rel.value.max %||% 1.1,
        max.rel.value.min = dots$max.rel.value.min %||% 0.9,
        max.overlap.threshold = dots$max.overlap.threshold %||% 0.15,
        min.overlap.threshold = dots$min.overlap.threshold %||% 0.15,
        p.mean.nbrs.dist.threshold = dots$p.mean.nbrs.dist.threshold %||% 0.9,
        p.mean.hopk.dist.threshold = dots$p.mean.hopk.dist.threshold %||% 0.9,
        p.deg.threshold = dots$p.deg.threshold %||% 0.9,
        min.basin.size = dots$min.basin.size %||% 10L,
        expand.basins = dots$expand.basins %||% TRUE,
        apply.relvalue.filter = dots$apply.relvalue.filter %||% TRUE,
        apply.maxima.clustering = dots$apply.maxima.clustering %||% TRUE,
        apply.minima.clustering = dots$apply.minima.clustering %||% TRUE,
        apply.geometric.filter = dots$apply.geometric.filter %||% TRUE,
        hop.k = dots$hop.k %||% 2L
    )

    ## Convert to 0-based indexing
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))

    ## Ensure Y has proper storage
    storage.mode(Y) <- "double"

    ## Call C++ implementation
    result <- .Call(
        S_compute_gfc_matrix,
        adj.list.0based,
        edge.length.list,
        Y,
        as.double(params$edge.length.quantile.thld),
        as.double(params$min.rel.value.max),
        as.double(params$max.rel.value.min),
        as.double(params$max.overlap.threshold),
        as.double(params$min.overlap.threshold),
        as.double(params$p.mean.nbrs.dist.threshold),
        as.double(params$p.mean.hopk.dist.threshold),
        as.double(params$p.deg.threshold),
        as.integer(params$min.basin.size),
        as.logical(params$expand.basins),
        as.logical(params$apply.relvalue.filter),
        as.logical(params$apply.maxima.clustering),
        as.logical(params$apply.minima.clustering),
        as.logical(params$apply.geometric.filter),
        as.integer(params$hop.k),
        as.integer(n.cores),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## Add class to each result
    for (i in seq_along(result)) {
        class(result[[i]]) <- c("gfc", "list")
        result[[i]]$parameters <- params
    }

    ## Add names from column names if available
    if (!is.null(colnames(Y))) {
        names(result) <- colnames(Y)
    }

    return(result)
}


## Null-coalescing operator (if not already defined)
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}


# =============================================================================
# Helper Functions for Basin Lookup
# =============================================================================

#' Resolve Basin ID to Vertex Index
#'
#' @param object A gfc object.
#' @param basin.id Either a character label (e.g., "M1") or numeric vertex index.
#' @return Integer vertex index (1-based).
#' @keywords internal
resolve.basin.vertex <- function(object, basin.id) {

    if (is.character(basin.id)) {
        ## Label-based lookup
        if (grepl("^M[0-9]+$", basin.id)) {
            idx <- which(names(object$max.basins) == basin.id)
            if (length(idx) == 1) {
                return(object$max.basins[[idx]]$vertex)
            }
        } else if (grepl("^m[0-9]+$", basin.id)) {
            idx <- which(names(object$min.basins) == basin.id)
            if (length(idx) == 1) {
                return(object$min.basins[[idx]]$vertex)
            }
        }
        stop(sprintf("Basin '%s' not found. Available basins: %s",
                     basin.id,
                     paste(c(names(object$max.basins),
                             names(object$min.basins)),
                           collapse = ", ")))

    } else if (is.numeric(basin.id)) {
        vertex.id <- as.integer(basin.id)

        ## Verify it exists
        for (b in object$max.basins) {
            if (b$vertex == vertex.id) return(vertex.id)
        }
        for (b in object$min.basins) {
            if (b$vertex == vertex.id) return(vertex.id)
        }
        stop(sprintf("Vertex %d is not a basin extremum in this GFC", vertex.id))

    } else {
        stop("basin.id must be either a character label or numeric vertex index")
    }
}


#' Get Basin Information
#'
#' @param object A gfc object.
#' @param basin.id Either a character label or numeric vertex index.
#' @return List with vertex, value, type, label, and basin object.
#' @keywords internal
get.basin.info <- function(object, basin.id) {

    basin <- NULL
    basin.label <- NULL
    basin.type <- NULL

    if (is.character(basin.id)) {
        basin.label <- basin.id

        if (grepl("^M[0-9]+$", basin.id)) {
            basin.type <- "max"
            idx <- which(names(object$max.basins) == basin.id)
            if (length(idx) == 1) {
                basin <- object$max.basins[[idx]]
            }
        } else if (grepl("^m[0-9]+$", basin.id)) {
            basin.type <- "min"
            idx <- which(names(object$min.basins) == basin.id)
            if (length(idx) == 1) {
                basin <- object$min.basins[[idx]]
            }
        }

        if (is.null(basin)) {
            stop(sprintf("Basin '%s' not found. Available basins: %s",
                         basin.id,
                         paste(c(names(object$max.basins),
                                 names(object$min.basins)),
                               collapse = ", ")))
        }

    } else if (is.numeric(basin.id)) {
        vertex.id <- as.integer(basin.id)

        for (i in seq_along(object$max.basins)) {
            if (object$max.basins[[i]]$vertex == vertex.id) {
                basin <- object$max.basins[[i]]
                basin.label <- names(object$max.basins)[i]
                basin.type <- "max"
                break
            }
        }

        if (is.null(basin)) {
            for (i in seq_along(object$min.basins)) {
                if (object$min.basins[[i]]$vertex == vertex.id) {
                    basin <- object$min.basins[[i]]
                    basin.label <- names(object$min.basins)[i]
                    basin.type <- "min"
                    break
                }
            }
        }

        if (is.null(basin)) {
            stop(sprintf("Vertex %d is not a basin extremum in this GFC", vertex.id))
        }

    } else {
        stop("basin.id must be either a character label or numeric vertex index")
    }

    list(
        vertex = basin$vertex,
        value = basin$value,
        type = basin.type,
        label = basin.label,
        basin = basin
    )
}


# =============================================================================
# Trajectories Generic and Methods
# =============================================================================

#' Extract Gradient Trajectories from GFC Result
#'
#' @description
#' Generic function for extracting gradient trajectories from a gradient flow
#' complex object.
#'
#' @param object An object containing gradient flow information.
#' @param ... Additional arguments passed to methods.
#'
#' @return Trajectory information, format depends on the method.
#'
#' @seealso \code{\link{trajectories.gfc}}
#'
#' @export
trajectories <- function(object, ...) {
    UseMethod("trajectories")
}


#' Extract Gradient Trajectories from GFC Result
#'
#' @description
#' Extracts gradient trajectories for a specified basin from a GFC result object.
#' Can summarize by raw trajectory terminals or by joined trajectory cells.
#'
#' @param object A gfc object from \code{compute.gfc()} computed with
#'   \code{with.trajectories = TRUE}.
#' @param basin.id Either a character label (e.g., "M1", "m3") or an integer
#'   vertex index identifying the basin.
#' @param summarize.by How to summarize trajectories:
#'   \itemize{
#'     \item \code{"terminal"}: Group by raw trajectory terminals (default).
#'       Shows all terminals including spurious extrema.
#'     \item \code{"cell"}: Group by joined trajectory cells. Shows only
#'       connections to non-spurious extrema after extension.
#'   }
#' @param include.paths Logical indicating whether to include the full
#'   trajectory paths in the output. Default is TRUE. Set to FALSE for
#'   summary information only.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with class \code{"gfc_trajectories"} containing trajectory
#'   information. The structure depends on \code{summarize.by}.
#'
#' @examples
#' \dontrun{
#' gfc <- compute.gfc(adj.list, edge.length.list, fitted.values,
#'                    with.trajectories = TRUE)
#'
#' ## Raw terminal view (shows spurious extrema)
#' traj.M1 <- trajectories(gfc, "M1", summarize.by = "terminal")
#' print(traj.M1)
#'
#' ## Cell view (shows joined trajectories to non-spurious extrema)
#' traj.M1.cells <- trajectories(gfc, "M1", summarize.by = "cell")
#' print(traj.M1.cells)
#' }
#'
#' @seealso \code{\link{compute.gfc}} for computing GFC with trajectories
#'
#' @export
trajectories.gfc <- function(object, basin.id,
                             summarize.by = c("terminal", "cell"),
                             include.paths = TRUE, ...) {

    if (!inherits(object, "gfc")) {
        stop("object must be of class 'gfc'")
    }

    summarize.by <- match.arg(summarize.by)

    if (summarize.by == "cell") {
        return(trajectories.by.cell(object, basin.id, include.paths))
    } else {
        return(trajectories.by.terminal(object, basin.id, include.paths))
    }
}


#' Summarize Trajectories by Raw Terminal
#' @keywords internal
trajectories.by.terminal <- function(object, basin.id, include.paths) {

    ## Get basin info
    info <- get.basin.info(object, basin.id)
    basin <- info$basin

    ## Check that trajectories were computed
    if (is.null(basin$trajectory.sets)) {
        stop("GFC object does not contain trajectory data. ",
             "Recompute with with.trajectories = TRUE")
    }

    ## Build set of non-spurious extrema for terminal classification
    non.spurious.vertices <- object$summary$vertex
    non.spurious.labels <- object$summary$label
    names(non.spurious.labels) <- as.character(non.spurious.vertices)

    ## Process trajectory sets
    terminal.groups <- list()
    total.trajectories <- 0

    for (traj.set in basin$trajectory.sets) {

        terminal.vertex <- traj.set$terminal.vertex
        terminal.vertex.char <- as.character(terminal.vertex)

        ## Classify terminal
        if (terminal.vertex %in% non.spurious.vertices) {
            terminal.label <- non.spurious.labels[terminal.vertex.char]
            terminal.type <- "non.spurious"
            group.name <- terminal.label
        } else {
            terminal.label <- NA_character_
            terminal.type <- "spurious"
            group.name <- paste0("spurious.", terminal.vertex)
        }

        ## Compute trajectory statistics
        n.traj <- length(traj.set$trajectories)
        total.trajectories <- total.trajectories + n.traj

        traj.lengths <- sapply(traj.set$trajectories, length)
        mean.length <- if (n.traj > 0) mean(traj.lengths) else NA_real_

        ## Build group entry
        group <- list(
            terminal.vertex = terminal.vertex,
            terminal.label = terminal.label,
            terminal.type = terminal.type,
            n.trajectories = n.traj,
            mean.length = mean.length
        )

        if (include.paths) {
            group$trajectories <- traj.set$trajectories
        }

        terminal.groups[[group.name]] <- group
    }

    ## Build result
    n.terminals <- length(terminal.groups)
    n.spurious <- sum(sapply(terminal.groups, function(g) g$terminal.type == "spurious"))

    result <- list(
        extremum = list(
            vertex = info$vertex,
            value = info$value,
            type = info$type,
            label = info$label
        ),
        terminal.groups = terminal.groups,
        n.terminals = n.terminals,
        n.spurious.terminals = n.spurious,
        n.trajectories = total.trajectories,
        summary.type = "terminal"
    )

    class(result) <- c("gfc_trajectories", "list")
    return(result)
}


#' Summarize Trajectories by Joined Cell
#' @keywords internal
trajectories.by.cell <- function(object, basin.id, include.paths) {

    ## Get basin info
    info <- get.basin.info(object, basin.id)

    ## Check that joined trajectories exist
    if (is.null(object$joined.trajectories) ||
        length(object$joined.trajectories) == 0) {
        stop("No joined trajectories available. ",
             "Recompute with with.trajectories = TRUE")
    }

    ## Filter joined trajectories for this basin
    if (info$type == "max") {
        cell.trajs <- Filter(function(jt) jt$max.vertex == info$vertex,
                             object$joined.trajectories)
    } else {
        cell.trajs <- Filter(function(jt) jt$min.vertex == info$vertex,
                             object$joined.trajectories)
    }

    ## Build lookup for opposite endpoint labels
    non.spurious.vertices <- object$summary$vertex
    non.spurious.labels <- object$summary$label
    names(non.spurious.labels) <- as.character(non.spurious.vertices)

    ## Group by opposite endpoint
    cell.groups <- list()

    for (jt in cell.trajs) {
        if (info$type == "max") {
            opposite.vertex <- jt$min.vertex
        } else {
            opposite.vertex <- jt$max.vertex
        }

        opposite.vertex.char <- as.character(opposite.vertex)

        ## Find label for opposite vertex
        if (opposite.vertex %in% non.spurious.vertices) {
            opposite.label <- non.spurious.labels[opposite.vertex.char]
        } else {
            opposite.label <- paste0("v", opposite.vertex)
        }

        group <- list(
            opposite.vertex = opposite.vertex,
            opposite.label = opposite.label,
            n.vertices = length(jt$path),
            path.length = jt$path.length,
            total.change = jt$total.change,
            n.intermediates = length(jt$intermediate.extrema),
            intermediate.extrema = jt$intermediate.extrema
        )

        if (include.paths) {
            group$path <- jt$path
        }

        cell.groups[[opposite.label]] <- group
    }

    ## Build result
    result <- list(
        extremum = list(
            vertex = info$vertex,
            value = info$value,
            type = info$type,
            label = info$label
        ),
        cell.groups = cell.groups,
        n.cells = length(cell.groups),
        summary.type = "cell"
    )

    class(result) <- c("gfc_trajectories", "list")
    return(result)
}


#' Print Method for GFC Trajectories
#'
#' @param x A gfc_trajectories object from trajectories.gfc()
#' @param ... Additional arguments (unused)
#'
#' @export
print.gfc_trajectories <- function(x, ...) {

    cat("GFC Trajectories\n")
    cat("================\n\n")

    cat(sprintf("Basin: %s (vertex %d, %s, value = %.4f)\n",
                x$extremum$label,
                x$extremum$vertex,
                x$extremum$type,
                x$extremum$value))

    if (x$summary.type == "terminal") {
        ## Terminal-based summary
        cat(sprintf("\nTerminal vertices: %d (%d non-spurious, %d spurious)\n",
                    x$n.terminals,
                    x$n.terminals - x$n.spurious.terminals,
                    x$n.spurious.terminals))

        cat(sprintf("Total trajectories: %d\n\n", x$n.trajectories))

        if (length(x$terminal.groups) > 0) {
            cat("Terminal groups:\n")

            for (name in names(x$terminal.groups)) {
                g <- x$terminal.groups[[name]]
                status <- if (g$terminal.type == "spurious") "[spurious]" else ""
                cat(sprintf("  %s (v%d) %s: %d trajectories, mean length %.1f\n",
                            name, g$terminal.vertex, status,
                            g$n.trajectories, g$mean.length))
            }
        }

    } else {
        ## Cell-based summary
        cat(sprintf("\nCells: %d (connections to non-spurious extrema)\n\n",
                    x$n.cells))

        if (length(x$cell.groups) > 0) {
            cat("Cells:\n")

            for (name in names(x$cell.groups)) {
                g <- x$cell.groups[[name]]
                intermediates <- if (g$n.intermediates > 0) {
                    sprintf(" via %d spurious", g$n.intermediates)
                } else {
                    " (direct)"
                }
                cat(sprintf("  %s -> %s%s: n.vertices = %d, path.length = %.3f, change = %.4f\n",
                            x$extremum$label, name, intermediates,
                            g$n.vertices, g$path.length, g$total.change))
            }
        }
    }

    invisible(x)
}


#' Extract Cell Trajectories from GFC Result
#'
#' @description
#' Extracts gradient flow trajectories for a specific cell (min-max pair)
#' from a GFC result object. A cell is defined by a pair of non-spurious
#' extrema (one minimum, one maximum).
#'
#' @param gfc A gfc object from compute.gfc() with trajectories computed.
#' @param min.id Minimum extremum identifier: either a label (e.g., "m4") or a
#'     vertex index (integer, 1-based).
#' @param max.id Maximum extremum identifier: either a label (e.g., "M1") or a
#'     vertex index (integer, 1-based).
#' @param map Optional integer vector mapping subgraph indices to original graph
#'     vertices. If provided, trajectory vertices are converted to subgraph
#'     indices. Specifically, \code{map[i]} gives the original graph vertex
#'     corresponding to subgraph vertex i. Vertices not found in map are
#'     returned as NA with a warning.
#'
#' @return A list of class "gfc_cell_trajectories" containing:
#'   \describe{
#'     \item{min.vertex}{Minimum vertex index (in output coordinate system)}
#'     \item{max.vertex}{Maximum vertex index (in output coordinate system)}
#'     \item{min.label}{Minimum label (e.g., "m4")}
#'     \item{max.label}{Maximum label (e.g., "M1")}
#'     \item{min.value}{Function value at minimum}
#'     \item{max.value}{Function value at maximum}
#'     \item{trajectories}{List of trajectory vertex vectors}
#'     \item{n.trajectories}{Number of trajectories}
#'     \item{mapped}{Logical; TRUE if vertices were mapped to subgraph indices}
#'     \item{original.min.vertex}{Original minimum vertex (before mapping)}
#'     \item{original.max.vertex}{Original maximum vertex (before mapping)}
#'   }
#'
#' @details
#' The function first resolves the min.id and max.id to vertex indices and
#' labels using the gfc summary. It then extracts all trajectories that
#' connect these two extrema.
#'
#' When working with a subgraph (e.g., the extended basin of a maximum),
#' the map parameter allows converting trajectory vertices to subgraph
#' indices. If map = M1.vertices where \code{M1.vertices[i]} is the original
#' graph vertex corresponding to subgraph vertex i, then the returned
#' trajectories will use subgraph indices 1, 2, ..., length(M1.vertices).
#'
#' @examples
#' \dontrun{
#' gfc <- compute.gfc(adj.list, weight.list, fitted.values,
#'                     with.trajectories = TRUE)
#'
#' # Extract by labels
#' cell.traj <- cell.trajectories.gfc(gfc, min.id = "m4", max.id = "M1")
#'
#' # Extract by vertex indices
#' cell.traj <- cell.trajectories.gfc(gfc, min.id = 509, max.id = 1147)
#'
#' # Extract and map to subgraph indices
#' M1.vertices <- gfc$max.basins[["M1"]]$vertices
#' cell.traj <- cell.trajectories.gfc(gfc, "m4", "M1", map = M1.vertices)
#' }
#'
#' @seealso \code{\link{compute.gfc}}, \code{\link{trajectories.gfc}}
#'
#' @export
cell.trajectories.gfc <- function(gfc,
                                   min.id,
                                   max.id,
                                   map = NULL) {

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    if (!inherits(gfc, "gfc")) {
        stop("gfc must be a gfc object from compute.gfc()")
    }

    if (is.null(gfc$joined.trajectories) && is.null(gfc$max.basins)) {
        stop("gfc does not contain trajectory information. ",
             "Use with.trajectories = TRUE in compute.gfc()")
    }

    summary.df <- gfc$summary

    ## ========================================================================
    ## Resolve min.id to vertex and label
    ## ========================================================================

    if (is.character(min.id)) {
        ## Label provided
        min.label <- min.id
        idx <- match(min.label, summary.df$label)
        if (is.na(idx)) {
            stop(sprintf("Minimum label '%s' not found in summary", min.label))
        }
        min.vertex <- summary.df$vertex[idx]
        min.value <- summary.df$value[idx]
        min.type <- summary.df$type[idx]
    } else if (is.numeric(min.id)) {
        ## Vertex index provided
        min.vertex <- as.integer(min.id)
        idx <- match(min.vertex, summary.df$vertex)
        if (is.na(idx)) {
            stop(sprintf("Minimum vertex %d not found in summary", min.vertex))
        }
        min.label <- summary.df$label[idx]
        min.value <- summary.df$value[idx]
        min.type <- summary.df$type[idx]
    } else {
        stop("min.id must be a character label or numeric vertex index")
    }

    ## Verify it's a minimum
    if (min.type != "min") {
        stop(sprintf("'%s' (vertex %d) is not a minimum, it is a %s",
                     min.label, min.vertex, min.type))
    }

    ## ========================================================================
    ## Resolve max.id to vertex and label
    ## ========================================================================

    if (is.character(max.id)) {
        ## Label provided
        max.label <- max.id
        idx <- match(max.label, summary.df$label)
        if (is.na(idx)) {
            stop(sprintf("Maximum label '%s' not found in summary", max.label))
        }
        max.vertex <- summary.df$vertex[idx]
        max.value <- summary.df$value[idx]
        max.type <- summary.df$type[idx]
    } else if (is.numeric(max.id)) {
        ## Vertex index provided
        max.vertex <- as.integer(max.id)
        idx <- match(max.vertex, summary.df$vertex)
        if (is.na(idx)) {
            stop(sprintf("Maximum vertex %d not found in summary", max.vertex))
        }
        max.label <- summary.df$label[idx]
        max.value <- summary.df$value[idx]
        max.type <- summary.df$type[idx]
    } else {
        stop("max.id must be a character label or numeric vertex index")
    }

    ## Verify it's a maximum
    if (max.type != "max") {
        stop(sprintf("'%s' (vertex %d) is not a maximum, it is a %s",
                     max.label, max.vertex, max.type))
    }

    ## ========================================================================
    ## Extract trajectories for this cell
    ## ========================================================================

    trajectories <- list()

    ## Check if joined trajectories are available
    if (!is.null(gfc$joined.trajectories)) {
        ## Joined trajectories: look for this specific cell
        for (traj in gfc$joined.trajectories) {
            if (traj$min.vertex == min.vertex && traj$max.vertex == max.vertex) {
                trajectories[[length(trajectories) + 1]] <- traj$path
            }
        }
    } else if (!is.null(gfc$max.basins[[max.label]])) {
        ## Fall back to max.basins trajectory sets
        max.basin <- gfc$max.basins[[max.label]]

        if (!is.null(max.basin$trajectory.sets)) {
            for (traj.set in max.basin$trajectory.sets) {
                ## Check if any trajectory in this set ends at min.vertex
                terminal <- traj.set$terminal.vertex
                if (terminal == min.vertex) {
                    ## Add all paths in this trajectory set
                    for (path in traj.set$trajectories) {
                        trajectories[[length(trajectories) + 1]] <- path
                    }
                }
            }
        }
    }

    if (length(trajectories) == 0) {
        warning(sprintf("No trajectories found for cell %s -> %s",
                        max.label, min.label))
    }

    ## ========================================================================
    ## Apply vertex mapping if provided
    ## ========================================================================

    original.min.vertex <- min.vertex
    original.max.vertex <- max.vertex
    mapped <- FALSE

    if (!is.null(map)) {
        mapped <- TRUE
        map <- as.integer(map)

        ## Map extrema vertices
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

        ## Map trajectory vertices
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

    ## ========================================================================
    ## Build result
    ## ========================================================================

    result <- list(
        min.vertex = min.vertex,
        max.vertex = max.vertex,
        min.label = min.label,
        max.label = max.label,
        min.value = min.value,
        max.value = max.value,
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
#' @param x A gfc_cell_trajectories object
#' @param max.print Maximum number of trajectories to print details for
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gfc_cell_trajectories <- function(x, max.print = 5, ...) {

    cat("GFC Cell Trajectories\n")
    cat("=====================\n")
    cat(sprintf("Cell: %s (vertex %d) -> %s (vertex %d)\n",
                x$max.label,
                if (x$mapped) x$original.max.vertex else x$max.vertex,
                x$min.label,
                if (x$mapped) x$original.min.vertex else x$min.vertex))
    cat(sprintf("Value change: %.4f -> %.4f (delta = %.4f)\n",
                x$max.value, x$min.value, x$max.value - x$min.value))
    cat(sprintf("Trajectories: %d\n", x$n.trajectories))

    if (x$mapped) {
        cat(sprintf("Vertices mapped to subgraph indices (min -> %s, max -> %s)\n",
                    ifelse(is.na(x$min.vertex), "NA", as.character(x$min.vertex)),
                    ifelse(is.na(x$max.vertex), "NA", as.character(x$max.vertex))))
    }

    if (x$n.trajectories > 0) {
        cat("\nTrajectory lengths:\n")
        lengths <- sapply(x$trajectories, length)
        cat(sprintf("  Min: %d, Max: %d, Mean: %.1f\n",
                    min(lengths), max(lengths), mean(lengths)))

        n.show <- min(max.print, x$n.trajectories)
        cat(sprintf("\nFirst %d trajectories:\n", n.show))
        for (i in seq_len(n.show)) {
            path <- x$trajectories[[i]]
            if (length(path) <= 10) {
                path.str <- paste(path, collapse = " -> ")
            } else {
                path.str <- paste(
                    c(paste(head(path, 4), collapse = " -> "),
                      "...",
                      paste(tail(path, 3), collapse = " -> ")),
                    collapse = " -> "
                )
            }
            cat(sprintf("  [%d] (%d vertices): %s\n", i, length(path), path.str))
        }

        if (x$n.trajectories > max.print) {
            cat(sprintf("  ... and %d more trajectories\n",
                        x$n.trajectories - max.print))
        }
    }

    invisible(x)
}


#' Get All Trajectory Vertices as a Single Vector
#'
#' @description
#' Extracts all unique vertices from cell trajectories.
#'
#' @param cell.traj A gfc_cell_trajectories object
#' @param unique Logical; if TRUE (default), return unique vertices only
#'
#' @return Integer vector of vertex indices
#'
#' @export
trajectory.vertices <- function(cell.traj, unique = TRUE) {

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
#' @description
#' Converts trajectory paths to an edge list suitable for graph operations
#' or visualization.
#'
#' @param cell.traj A gfc_cell_trajectories object
#' @param weighted Logical; if TRUE, include edge counts as weights
#'
#' @return A data frame with columns 'from', 'to', and optionally 'weight'
#'
#' @export
trajectory.edges <- function(cell.traj, weighted = FALSE) {

    if (!inherits(cell.traj, "gfc_cell_trajectories")) {
        stop("cell.traj must be a gfc_cell_trajectories object")
    }

    if (cell.traj$n.trajectories == 0) {
        return(data.frame(from = integer(0), to = integer(0)))
    }

    ## Collect all edges
    edges <- list()
    for (path in cell.traj$trajectories) {
        if (length(path) < 2) next
        for (i in seq_len(length(path) - 1)) {
            edges[[length(edges) + 1]] <- c(path[i], path[i + 1])
        }
    }

    if (length(edges) == 0) {
        return(data.frame(from = integer(0), to = integer(0)))
    }

    edge.df <- data.frame(
        from = sapply(edges, `[`, 1),
        to = sapply(edges, `[`, 2)
    )

    if (weighted) {
        ## Count edge occurrences
        edge.key <- paste(edge.df$from, edge.df$to, sep = "-")
        edge.counts <- table(edge.key)

        ## Get unique edges with counts
        unique.edges <- unique(edge.df)
        unique.key <- paste(unique.edges$from, unique.edges$to, sep = "-")
        unique.edges$weight <- as.integer(edge.counts[unique.key])

        return(unique.edges)
    }

    return(unique(edge.df))
}


#' Draw a Single Cell Trajectory in 3D Space
#'
#' Visualizes a single trajectory from a cell object in 3D using rgl graphics.
#' The function draws spheres at trajectory vertices, connecting segments, and
#' optionally directional arrows along the path.
#'
#' @param graph.3d A matrix of 3D coordinates for graph vertices
#' @param i Integer index specifying which trajectory to draw from the cell's
#'   trajectories list
#' @param cell List object containing trajectory data with components:
#'   \itemize{
#'     \item \code{trajectories}: List of trajectory vertex sequences
#'     \item \code{terminal.vertex}: Terminal vertex identifier
#'   }
#' @param with.arrows Logical; if TRUE, draws directional arrows along trajectory
#'   segments. Default is FALSE
#' @param arrow.size Numeric scaling factor for arrow size. Default is 1
#' @param col Character string specifying color for trajectory spheres.
#'   Default is "cyan"
#' @param terminal.radius Numeric radius for terminal vertex spheres.
#'   Default is 0.3
#' @param radius Numeric radius for trajectory vertex spheres. Default is 0.15
#' @param terminal.vertex.adj Numeric vector of length 2 for text label
#'   adjustment. Default is c(0,0)
#' @param terminal.vertex.cex Numeric character expansion factor for terminal
#'   vertex labels. Default is 3
#'
#' @details
#' Trajectory vertices are drawn as cyan spheres connected by gray segments.
#' When \code{with.arrows = TRUE}, red arrows are drawn at 40-60% along each
#' segment to indicate trajectory direction.
#'
#' @return None. Function is called for its side effect of adding 3D graphics
#'   to the current rgl device.
#'
#' @examples
#' \dontrun{
#' # Assumes graph.3d and cell object are defined
#' draw.cell.trajectory(graph.3d, 1, my.cell, with.arrows = TRUE, arrow.size = 1.5)
#' }
#'
#' @export
draw.cell.trajectory <- function(graph.3d,
                                 i,
                                 cell,
                                 with.arrows = FALSE,
                                 arrow.size = 1,
                                 col = "cyan",
                                 terminal.radius = 0.3,
                                 radius = 0.15,
                                 terminal.vertex.adj = c(0,0),
                                 terminal.vertex.cex = 3) {
    traj <- cell$trajectories[[i]]
    n <- length(traj)

    rgl::spheres3d(graph.3d[cell$min.vertex,], radius = terminal.radius, col = col)
    rgl::texts3d(graph.3d[cell$min.vertex,], texts = cell$min.vertex,
                 adj = terminal.vertex.adj, cex = terminal.vertex.cex)

    rgl::spheres3d(graph.3d[cell$max.vertex,], radius = terminal.radius, col = col)
    rgl::texts3d(graph.3d[cell$max.vertex,], texts = cell$max.vertex,
                 adj = terminal.vertex.adj, cex = terminal.vertex.cex)

    rgl::spheres3d(graph.3d[traj,], radius = radius, col = "cyan")

    ## Create pairs of consecutive points
    segment.indices <- c(rbind(traj[-n], traj[-1]))
    M <- graph.3d[segment.indices, ]
    rgl::segments3d(M, col = "gray", lwd = 5)

    if (with.arrows) {
        for (j in 1:(n - 1)) {
            start <- graph.3d[traj[j], ]
            end <- graph.3d[traj[j + 1], ]
            mid.start <- start * 0.6 + end * 0.4  # start arrow at 40% along edge
            mid.end <- start * 0.4 + end * 0.6    # end arrow at 60% along edge
            rgl::arrow3d(mid.start, mid.end, type = "flat", col = "red",
                        width = 0.5, s = arrow.size)
        }
    }
}

#' Draw All Trajectories for a Cell in 3D Space
#'
#' Visualizes all trajectories from a cell object in 3D using rgl graphics.
#' Iterates through all trajectories in the cell, drawing spheres at vertices,
#' connecting segments, and optionally directional arrows.
#'
#' @param cell List object containing trajectory data with components:
#'   \itemize{
#'     \item \code{trajectories}: List of trajectory vertex sequences
#'     \item \code{terminal.vertex}: Terminal vertex identifier
#'   }
#' @param with.edges Logical; if TRUE, draws line segments connecting trajectory
#'   vertices. Default is TRUE
#' @param with.arrows Logical; if TRUE, draws directional arrows along trajectory
#'   segments. Default is TRUE
#' @param arrow.size Numeric scaling factor for arrow size. Default is 1
#' @param arrow.width Numeric line width for connecting segments. Default is 1
#' @param col Character string specifying color for trajectory spheres.
#'   Default is "cyan"
#' @param terminal.radius Numeric radius for terminal vertex spheres.
#'   Default is 0.3
#' @param radius Numeric radius for trajectory vertex spheres. Default is 0.15
#' @param terminal.vertex.adj Numeric vector of length 2 for text label
#'   adjustment. Default is c(0,0)
#' @param terminal.vertex.cex Numeric character expansion factor for terminal
#'   vertex labels. Default is 2
#'
#' @details
#' The function expects a global \code{graph.3d} object containing 3D coordinates
#' for graph vertices. For each trajectory in the cell, vertices are drawn as
#' spheres in the specified color, with optional connecting segments and
#' directional arrows. Arrows are positioned at 40-60% along each segment.
#'
#' @return None. Function is called for its side effect of adding 3D graphics
#'   to the current rgl device.
#'
#' @examples
#' \dontrun{
#' # Assumes graph.3d and cell object are defined
#' draw.cell.trajectories(my.cell, with.arrows = FALSE, with.edges = TRUE)
#' }
#'
#' @export
draw.cell.trajectories <- function(cell,
                                   with.edges = TRUE,
                                   with.arrows = TRUE,
                                   arrow.size = 1,
                                   arrow.width = 1,
                                   col = "cyan",
                                   terminal.radius = 0.3,
                                   radius = 0.15,
                                   terminal.vertex.adj = c(0,0),
                                   terminal.vertex.cex = 2) {
    for (i in seq_along(cell$trajectories)) {
        traj <- cell$trajectories[[i]]
        n <- length(traj)

        # NOTE: Remove or fix this line - M1.cells appears to be hardcoded
        # spheres3d(graph.3d[M1.cells$max.vertex,], radius = terminal.radius, col = col)

        spheres3d(graph.3d[cell$terminal.vertex,], radius = terminal.radius, col = col)
        texts3d(graph.3d[cell$terminal.vertex,], texts = cell$terminal.vertex,
                adj = terminal.vertex.adj, cex = terminal.vertex.cex)
        spheres3d(graph.3d[traj,], radius = radius, col = col)

        ## Create edges
        if (with.edges) {
            segment.indices <- c(rbind(traj[-n], traj[-1]))
            M <- graph.3d[segment.indices, ]
            rgl::segments3d(M, col = "gray", lwd = arrow.width)
        }

        if (with.arrows) {
            for (j in 1:(n - 1)) {
                start <- graph.3d[traj[j], ]
                end <- graph.3d[traj[j + 1], ]
                mid.start <- start * 0.6 + end * 0.4  # start arrow at 40% along edge
                mid.end <- start * 0.4 + end * 0.6    # end arrow at 60% along edge
                rgl::arrow3d(mid.start, mid.end, type = "flat", col = "red",
                            width = 0.5, s = arrow.size)
            }
        }
    }
}

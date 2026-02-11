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

    ## -------------------------------------------------------
    ## Breaking ties (if any)
    ## -------------------------------------------------------
    fitted.values <- break.ties(fitted.values,
                                noise.scale = 1e-15,
                                min.abs.noise = 1e-16,
                                preserve.bounds = TRUE,
                                seed = NULL,
                                verbose = FALSE)

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
    ## Fix assignment index mapping to match summary label order
    ## ========================================================================

    if (nrow(result$summary) > 0) {
        ## The expanded assignments use indices based on basin vector position,
        ## but summary labels (M1, M2, ...) are assigned by sorted value.
        ## We need to create a mapping from vector position to label index.

        ## Get maxima in summary order (this is how labels were assigned)
        max.summary <- result$summary[result$summary$type == "max", ]
        min.summary <- result$summary[result$summary$type == "min", ]

        ## Build mapping: for each basin in max.basins (by position),
        ## find which label index (1, 2, ...) it should have
        if (length(result$max.basins) > 0 && length(result$expanded.max.assignment) > 0) {
            ## Get extremum vertices in basin vector order
            basin.vertices.order <- sapply(result$max.basins, function(b) b$vertex)

            ## For each basin position, find its row in the sorted summary
            ## The summary row index IS the label number (M1 = row 1 among maxima)
            position.to.label.idx <- match(basin.vertices.order, max.summary$vertex)

            ## Now remap the assignment vector
            ## Old: assignment[v] = position in max.basins (1-based)
            ## New: assignment[v] = label index (1 = M1, 2 = M2, ...)
            old.assignment <- result$expanded.max.assignment
            new.assignment <- rep(NA_integer_, length(old.assignment))

            for (v in seq_along(old.assignment)) {
                old.idx <- old.assignment[v]
                if (!is.na(old.idx) && old.idx >= 1 && old.idx <= length(position.to.label.idx)) {
                    new.assignment[v] <- position.to.label.idx[old.idx]
                }
            }
            result$expanded.max.assignment <- new.assignment
        }

        ## Same fix for minima
        if (length(result$min.basins) > 0 && length(result$expanded.min.assignment) > 0) {
            basin.vertices.order <- sapply(result$min.basins, function(b) b$vertex)
            position.to.label.idx <- match(basin.vertices.order, min.summary$vertex)

            old.assignment <- result$expanded.min.assignment
            new.assignment <- rep(NA_integer_, length(old.assignment))

            for (v in seq_along(old.assignment)) {
                old.idx <- old.assignment[v]
                if (!is.na(old.idx) && old.idx >= 1 && old.idx <= length(position.to.label.idx)) {
                    new.assignment[v] <- position.to.label.idx[old.idx]
                }
            }
            result$expanded.min.assignment <- new.assignment
        }

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
        return(trajectories_by_cell(object, basin.id, include.paths))
    } else {
        return(trajectories_by_terminal(object, basin.id, include.paths))
    }
}


#' Summarize Trajectories by Raw Terminal
#' @keywords internal
trajectories_by_terminal <- function(object, basin.id, include.paths) {

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
trajectories_by_cell <- function(object, basin.id, include.paths) {

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

#' Extract Cell Trajectories from GFC Objects
#'
#' Generic function to extract gradient flow trajectories for a specific
#' cell (min-max pair) from GFC-related objects.
#'
#' @param x A GFC object.
#' @param min.id Minimum extremum identifier.
#' @param max.id Maximum extremum identifier.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object containing cell trajectory information.
#'
#' @seealso \code{\link{cell.trajectories.gfc}}
#'
#' @export
cell.trajectories <- function(x, min.id, max.id, ...) {
    UseMethod("cell.trajectories")
}

#' Extract Cell Trajectories from GFC Result
#'
#' @description
#' Extracts gradient flow trajectories for a specific cell (min-max pair)
#' from a GFC result object. A cell is defined by a pair of non-spurious
#' extrema (one minimum, one maximum).
#'
#' @param x A gfc object from compute.gfc() with trajectories computed.
#' @param gfc Backward-compatible alias for \code{x}.
#' @param min.id Minimum extremum identifier: either a label (e.g., "m4") or a
#'     vertex index (integer, 1-based).
#' @param max.id Maximum extremum identifier: either a label (e.g., "M1") or a
#'     vertex index (integer, 1-based).
#' @param map Optional integer vector mapping subgraph indices to original graph
#'     vertices. If provided, trajectory vertices are converted to subgraph
#'     indices. Specifically, \code{map[i]} gives the original graph vertex
#'     corresponding to subgraph vertex i. Vertices not found in map are
#'     returned as NA with a warning.
#' @param ... Additional arguments (currently ignored).
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
cell.trajectories.gfc <- function(x,
                                  min.id,
                                  max.id,
                                  map = NULL,
                                  gfc = x,
                                  ...) {
    x <- gfc
    gfc <- x

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

#' Draw a Trajectory in 3D Space
#'
#' Visualizes a trajectory as a sequence of vertices in 3D using rgl graphics.
#' The function draws spheres at trajectory vertices with highlighted terminal
#' points, connecting segments between consecutive vertices, and optionally
#' directional arrows indicating flow direction along the path.
#'
#' @param graph.3d A matrix of 3D coordinates for graph vertices where each row
#'   corresponds to a vertex and columns represent x, y, z coordinates.
#' @param vertices Integer vector specifying the ordered sequence of vertex
#'   indices defining the trajectory. The first element is treated as the
#'   source and the last as the destination.
#' @param with.arrows Logical; if TRUE, draws directional arrows along trajectory
#'   segments to indicate flow direction. Default is FALSE.
#' @param arrow.size Numeric scaling factor for arrow size when with.arrows is
#'   TRUE. Default is 1.
#' @param col Character string specifying color for trajectory spheres.
#'   Default is "cyan".
#' @param terminal.col Character string specifying color for terminal vertex
#'   spheres. If NULL, uses the value of col. Default is NULL.
#' @param segment.col Character string specifying color for connecting segments.
#'   Default is "gray".
#' @param arrow.col Character string specifying color for directional arrows.
#'   Default is "red".
#' @param terminal.radius Numeric radius for terminal vertex spheres.
#'   Default is 0.3.
#' @param radius Numeric radius for interior trajectory vertex spheres.
#'   Default is 0.15.
#' @param segment.lwd Numeric line width for connecting segments. Default is 5.
#' @param with.labels Logical; if TRUE, draws text labels at terminal vertices
#'   showing vertex indices. Default is TRUE.
#' @param terminal.vertex.adj Numeric vector of length 2 for text label
#'   adjustment relative to terminal vertices. Default is c(0,0).
#' @param terminal.vertex.cex Numeric character expansion factor for terminal
#'   vertex labels. Default is 3.
#'
#' @details
#' The trajectory is rendered as a sequence of spheres connected by line
#' segments. Terminal vertices at the start and end of the trajectory are
#' drawn with larger spheres to distinguish them from interior vertices.
#' When with.arrows is TRUE, flat arrows are drawn at the midpoint of each
#' segment spanning from 40 percent to 60 percent along the edge to indicate
#' the direction of flow from source to destination.
#'
#' @return None. Function is called for its side effect of adding 3D graphics
#'   to the current rgl device.
#'
#' @examples
#' \dontrun{
#' ## Draw a simple trajectory through vertices 1, 5, 12, 8
#' draw.trajectory(graph.3d, vertices = c(1, 5, 12, 8))
#'
#' ## Draw with directional arrows and custom colors
#' draw.trajectory(graph.3d, vertices = c(1, 5, 12, 8),
#'                 with.arrows = TRUE, col = "blue", arrow.col = "orange")
#' }
#'
#' @export
draw.trajectory <- function(graph.3d,
                            vertices,
                            with.arrows = FALSE,
                            arrow.size = 1,
                            col = "cyan",
                            terminal.col = NULL,
                            segment.col = "gray",
                            arrow.col = "red",
                            terminal.radius = 0.3,
                            radius = 0.15,
                            segment.lwd = 5,
                            with.labels = TRUE,
                            terminal.vertex.adj = c(0, 0),
                            terminal.vertex.cex = 3) {

    n <- length(vertices)
    if (n < 2) {
        warning("Trajectory must contain at least 2 vertices")
        return(invisible(NULL))
    }

    if (is.null(terminal.col)) {
        terminal.col <- col
    }

    source.vertex <- vertices[1]
    dest.vertex <- vertices[n]

    ## Draw terminal vertices with larger spheres
    rgl::spheres3d(graph.3d[source.vertex, ], radius = terminal.radius,
                   col = terminal.col)
    rgl::spheres3d(graph.3d[dest.vertex, ], radius = terminal.radius,
                   col = terminal.col)

    ## Draw labels at terminal vertices
    if (with.labels) {
        rgl::texts3d(graph.3d[source.vertex, ], texts = source.vertex,
                     adj = terminal.vertex.adj, cex = terminal.vertex.cex)
        rgl::texts3d(graph.3d[dest.vertex, ], texts = dest.vertex,
                     adj = terminal.vertex.adj, cex = terminal.vertex.cex)
    }

    ## Draw interior vertices with smaller spheres
    if (n > 2) {
        interior.vertices <- vertices[2:(n - 1)]
        rgl::spheres3d(graph.3d[interior.vertices, ], radius = radius, col = col)
    }

    ## Draw connecting segments
    segment.indices <- c(rbind(vertices[-n], vertices[-1]))
    M <- graph.3d[segment.indices, ]
    rgl::segments3d(M, col = segment.col, lwd = segment.lwd)

    ## Draw directional arrows if requested
    if (with.arrows) {
        for (j in 1:(n - 1)) {
            start <- graph.3d[vertices[j], ]
            end <- graph.3d[vertices[j + 1], ]
            mid.start <- start * 0.6 + end * 0.4
            mid.end <- start * 0.4 + end * 0.6
            rgl::arrow3d(mid.start, mid.end, type = "flat", col = arrow.col,
                         width = 0.5, s = arrow.size)
        }
    }

    invisible(NULL)
}

#' Draw a Single Cell Trajectory in 3D Space
#'
#' Visualizes a single trajectory from a cell object in 3D using rgl graphics.
#' This is a convenience wrapper around \code{\link{draw.trajectory}} that
#' extracts the trajectory vertex sequence from a cell's trajectory list.
#'
#' @param graph.3d A matrix of 3D coordinates for graph vertices where each row
#'   corresponds to a vertex and columns represent x, y, z coordinates.
#' @param i Integer index specifying which trajectory to draw from the cell's
#'   trajectories list.
#' @param cell List object containing trajectory data with component
#'   \code{trajectories}, a list of trajectory vertex sequences.
#' @inheritDotParams draw.trajectory -vertices
#'
#' @return None. Function is called for its side effect of adding 3D graphics
#'   to the current rgl device.
#'
#' @seealso \code{\link{draw.trajectory}} for the underlying drawing function
#'   and full parameter documentation.
#'
#' @examples
#' \dontrun{
#' draw.cell.trajectory(graph.3d, 1, my.cell, with.arrows = TRUE)
#' }
#'
#' @export
draw.cell.trajectory <- function(graph.3d, i, cell, ...) {
    draw.trajectory(graph.3d, vertices = cell$trajectories[[i]], ...)
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
        # rgl::spheres3d(graph.3d[M1.cells$max.vertex,], radius = terminal.radius, col = col)

        rgl::spheres3d(graph.3d[cell$terminal.vertex,], radius = terminal.radius, col = col)
        rgl::texts3d(graph.3d[cell$terminal.vertex,], texts = cell$terminal.vertex,
                adj = terminal.vertex.adj, cex = terminal.vertex.cex)
        rgl::spheres3d(graph.3d[traj,], radius = radius, col = col)

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

#' Diagnose Missing Trajectory Vertices
#'
#' Identifies vertices in a basin that are not covered by any trajectory.
#'
#' @param gfc A gfc object with trajectories.
#' @param max.id Maximum basin identifier (label or vertex).
#' @param debug Logical; if TRUE, print detailed intermediate diagnostics.
#'
#' @return List with diagnostic information.
#'
#' @export
diagnose.trajectory.coverage <- function(gfc, max.id, debug = TRUE) {

    ## Resolve max.id
    if (is.character(max.id)) {
        max.label <- max.id
        idx <- match(max.label, gfc$summary$label)
        max.vertex <- gfc$summary$vertex[idx]
    } else {
        max.vertex <- as.integer(max.id)
        idx <- match(max.vertex, gfc$summary$vertex)
        max.label <- gfc$summary$label[idx]
    }

    ## Debug: understand the assignment encoding
    if (debug) {
        cat("=== Debug Information ===\n")
        cat(sprintf("max.label: %s, max.vertex: %d\n", max.label, max.vertex))

        ## Check what values are in expanded.max.assignment
        cat("\nexpanded.max.assignment unique values:\n")
        print(table(gfc$expanded.max.assignment, useNA = "ifany"))

        ## Check the summary to understand maxima ordering
        max.summary <- gfc$summary[gfc$summary$type == "max", ]
        cat("\nMaxima in summary (order matters):\n")
        print(max.summary[, c("label", "vertex")])

        ## What index is M1 in the maxima list?
        max.labels <- gfc$summary$label[gfc$summary$type == "max"]
        max.position <- match(max.label, max.labels)
        cat(sprintf("\n%s is at position %d in maxima list\n", max.label, max.position))
    }

    ## Get basin object
    basin <- gfc$max.basins[[max.label]]
    traj.vertices <- basin$vertices

    if (debug) {
        cat(sprintf("\nbasin$vertices: n=%d, range=[%d, %d]\n",
                    length(traj.vertices), min(traj.vertices), max(traj.vertices)))
        cat(sprintf("First 10: %s\n", paste(head(traj.vertices, 10), collapse = ", ")))

        ## Check if max.vertex is in basin$vertices
        cat(sprintf("max.vertex %d in basin$vertices: %s\n",
                    max.vertex, max.vertex %in% traj.vertices))
    }

    ## Try different interpretations of the assignment
    max.index.from.label <- as.integer(sub("M", "", max.label))
    max.position.in.list <- match(max.label,
                                   gfc$summary$label[gfc$summary$type == "max"])

    if (debug) {
        cat(sprintf("\nTrying different assignment interpretations:\n"))
        cat(sprintf("  Using label number (M1 -> 1): %d vertices\n",
                    sum(gfc$expanded.max.assignment == max.index.from.label, na.rm = TRUE)))
        cat(sprintf("  Using position in max list: %d vertices\n",
                    sum(gfc$expanded.max.assignment == max.position.in.list, na.rm = TRUE)))
        cat(sprintf("  Using vertex index directly: %d vertices\n",
                    sum(gfc$expanded.max.assignment == max.vertex, na.rm = TRUE)))
    }

    ## Determine correct interpretation
    ## Check which interpretation puts max.vertex in its own basin
    for (interp in c("label_number", "position", "vertex")) {
        test.idx <- switch(interp,
                           label_number = max.index.from.label,
                           position = max.position.in.list,
                           vertex = max.vertex)
        test.basin <- which(gfc$expanded.max.assignment == test.idx)
        if (max.vertex %in% test.basin) {
            if (debug) {
                cat(sprintf("\n>>> Correct interpretation: '%s' (value=%d)\n",
                            interp, test.idx))
            }
            basin.vertices <- test.basin
            break
        }
    }

    if (!exists("basin.vertices")) {
        warning("Could not determine correct assignment interpretation")
        basin.vertices <- which(gfc$expanded.max.assignment == max.index.from.label)
    }

    ## Now compare
    missing.vertices <- setdiff(basin.vertices, traj.vertices)
    extra.vertices <- setdiff(traj.vertices, basin.vertices)
    overlap.vertices <- intersect(basin.vertices, traj.vertices)

    if (debug) {
        cat(sprintf("\n=== Coverage Analysis ===\n"))
        cat(sprintf("Expanded basin vertices: %d\n", length(basin.vertices)))
        cat(sprintf("Trajectory vertices: %d\n", length(traj.vertices)))
        cat(sprintf("Overlap: %d\n", length(overlap.vertices)))
        cat(sprintf("Missing (in basin, not in traj): %d\n", length(missing.vertices)))
        cat(sprintf("Extra (in traj, not in basin): %d\n", length(extra.vertices)))
    }

    result <- list(
        max.label = max.label,
        max.vertex = max.vertex,
        basin.vertices = basin.vertices,
        traj.vertices = traj.vertices,
        n.basin.vertices = length(basin.vertices),
        n.traj.vertices = length(traj.vertices),
        n.overlap = length(overlap.vertices),
        n.missing = length(missing.vertices),
        n.extra = length(extra.vertices),
        missing.vertices = sort(missing.vertices),
        extra.vertices = sort(extra.vertices),
        overlap.vertices = sort(overlap.vertices)
    )

    return(invisible(result))
}

old.diagnose.trajectory.coverage <- function(gfc, max.id) {

    ## Resolve max.id
    if (is.character(max.id)) {
        max.label <- max.id
        idx <- match(max.label, gfc$summary$label)
        max.vertex <- gfc$summary$vertex[idx]
    } else {
        max.vertex <- as.integer(max.id)
        idx <- match(max.vertex, gfc$summary$vertex)
        max.label <- gfc$summary$label[idx]
    }

    max.index <- as.integer(sub("M","",max.label))

    ## Get basin vertices
    basin <- gfc$max.basins[[max.label]]
    traj.vertices <- basin$vertices
    basin.vertices <- which(gfc$expanded.max.assignment == max.index)

    ## Find vertices in basin but not in any trajectory
    missing.vertices <- setdiff(basin.vertices, traj.vertices)

    ## Find vertices in trajectories but not in basin (shouldn't happen)
    extra.vertices <- setdiff(traj.vertices, basin.vertices)

    result <- list(
        max.label = max.label,
        max.vertex = max.vertex,
        n.basin.vertices = length(basin.vertices),
        n.traj.vertices = length(traj.vertices),
        n.missing = length(missing.vertices),
        n.extra = length(extra.vertices),
        missing.vertices = sort(missing.vertices),
        extra.vertices = sort(extra.vertices),
        coverage = length(traj.vertices) / length(basin.vertices)
    )

    cat(sprintf("Trajectory Coverage Diagnosis for %s (vertex %d)\n",
                max.label, max.vertex))
    cat(sprintf("  Basin vertices: %d\n", result$n.basin.vertices))
    cat(sprintf("  Trajectory vertices: %d\n", result$n.traj.vertices))
    cat(sprintf("  Missing from trajectories: %d (%.1f%%)\n",
                result$n.missing, 100 * (1 - result$coverage)))
    cat(sprintf("  Extra (not in basin): %d\n", result$n.extra))

    if (result$n.missing > 0 && result$n.missing <= 20) {
        cat(sprintf("  Missing vertices: %s\n",
                    paste(result$missing.vertices, collapse = ", ")))
    }

    return(invisible(result))
}

#' Trace Gradient Flow from a Single Vertex
#'
#' Manually traces the ascending gradient flow from a vertex to diagnose
#' why it may not be part of expected trajectories.
#'
#' @param adj.list Graph adjacency list.
#' @param weight.list Edge weight list (can be modulated).
#' @param y Function values.
#' @param start.vertex Starting vertex (1-based).
#' @param max.steps Maximum steps to prevent infinite loops.
#' @param verbose Print each step.
#'
#' @return List with trajectory and diagnostic info.
#'
#' @export
trace.gradient.flow <- function(adj.list,
                                 weight.list,
                                 y,
                                 start.vertex,
                                 max.steps = 100,
                                 verbose = TRUE) {

    trajectory <- start.vertex
    current <- start.vertex
    step <- 0

    if (verbose) {
        cat(sprintf("Tracing ascending gradient from vertex %d (y = %.4f)\n",
                    start.vertex, y[start.vertex]))
    }

    while (step < max.steps) {
        nbrs <- adj.list[[current]]
        weights <- weight.list[[current]]

        if (length(nbrs) == 0) {
            if (verbose) cat("  No neighbors - isolated vertex\n")
            break
        }

        ## Find neighbor with maximum y value (ascending gradient)
        y.nbrs <- y[nbrs]
        best.idx <- which.max(y.nbrs)
        best.nbr <- nbrs[best.idx]
        best.y <- y.nbrs[best.idx]

        if (verbose) {
            cat(sprintf("  Step %d: v=%d (y=%.4f) -> ", step + 1, current, y[current]))
        }

        ## Check if we can ascend
        if (best.y <= y[current]) {
            if (verbose) {
                cat(sprintf("LOCAL MAX (no higher neighbor)\n"))
                cat(sprintf("    Neighbors: %s\n", paste(nbrs, collapse = ", ")))
                cat(sprintf("    Neighbor y: %s\n",
                            paste(round(y.nbrs, 4), collapse = ", ")))
            }
            break
        }

        if (verbose) {
            cat(sprintf("v=%d (y=%.4f, delta=%.4f)\n",
                        best.nbr, best.y, best.y - y[current]))
        }

        ## Check for cycle
        if (best.nbr %in% trajectory) {
            if (verbose) cat("  CYCLE DETECTED!\n")
            break
        }

        trajectory <- c(trajectory, best.nbr)
        current <- best.nbr
        step <- step + 1
    }

    if (step >= max.steps) {
        if (verbose) cat("  Max steps reached\n")
    }

    terminal <- trajectory[length(trajectory)]

    result <- list(
        start.vertex = start.vertex,
        terminal.vertex = terminal,
        trajectory = trajectory,
        n.steps = length(trajectory) - 1,
        start.y = y[start.vertex],
        terminal.y = y[terminal],
        delta.y = y[terminal] - y[start.vertex]
    )

    if (verbose) {
        cat(sprintf("\nSummary: %d -> %d in %d steps, y: %.4f -> %.4f\n",
                    start.vertex, terminal, result$n.steps,
                    result$start.y, result$terminal.y))
    }

    return(invisible(result))
}

#' Diagnose Gradient Flow for Missing Vertices
#'
#' Traces gradient flow for vertices that are missing from trajectories
#' to understand where they actually flow.
#'
#' @param gfc A gfc object.
#' @param adj.list Graph adjacency list.
#' @param weight.list Edge weight list used in compute.gfc().
#' @param y Function values.
#' @param max.id Maximum basin to diagnose.
#' @param max.diagnose Maximum number of missing vertices to diagnose.
#'
#' @return Data frame with flow destinations for missing vertices.
#'
#' @export
diagnose.missing.flows <- function(gfc,
                                    adj.list,
                                    weight.list,
                                    y,
                                    max.id,
                                    max.diagnose = 50) {

    ## Get missing vertices
    coverage <- diagnose.trajectory.coverage(gfc, max.id)
    missing <- coverage$missing.vertices

    if (length(missing) == 0) {
        cat("No missing vertices to diagnose.\n")
        return(NULL)
    }

    n.diagnose <- min(length(missing), max.diagnose)
    cat(sprintf("\nDiagnosing gradient flow for %d missing vertices...\n\n",
                n.diagnose))

    ## Get expected maximum vertex
    expected.max <- coverage$max.vertex

    ## Trace each missing vertex
    results <- data.frame(
        vertex = integer(n.diagnose),
        terminal = integer(n.diagnose),
        reaches.expected = logical(n.diagnose),
        n.steps = integer(n.diagnose),
        y.start = numeric(n.diagnose),
        y.terminal = numeric(n.diagnose)
    )

    terminal.counts <- list()

    for (i in seq_len(n.diagnose)) {
        v <- missing[i]
        flow <- trace.gradient.flow(adj.list, weight.list, y, v, verbose = FALSE)

        results$vertex[i] <- v
        results$terminal[i] <- flow$terminal.vertex
        results$reaches.expected[i] <- (flow$terminal.vertex == expected.max)
        results$n.steps[i] <- flow$n.steps
        results$y.start[i] <- flow$start.y
        results$y.terminal[i] <- flow$terminal.y

        ## Count terminals
        term.key <- as.character(flow$terminal.vertex)
        if (is.null(terminal.counts[[term.key]])) {
            terminal.counts[[term.key]] <- 0
        }
        terminal.counts[[term.key]] <- terminal.counts[[term.key]] + 1
    }

    ## Summary
    cat("Terminal destination summary:\n")
    for (term in names(terminal.counts)) {
        term.v <- as.integer(term)
        ## Find label if it's a known extremum
        label.idx <- match(term.v, gfc$summary$vertex)
        if (!is.na(label.idx)) {
            label <- gfc$summary$label[label.idx]
        } else {
            label <- "(not in summary)"
        }
        cat(sprintf("  Vertex %s %s: %d vertices\n",
                    term, label, terminal.counts[[term]]))
    }

    cat(sprintf("\nReach expected maximum (%d): %d / %d (%.1f%%)\n",
                expected.max,
                sum(results$reaches.expected),
                n.diagnose,
                100 * mean(results$reaches.expected)))

    return(results)
}

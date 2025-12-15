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
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## ========================================================================
    ## Add names to basins based on summary labels
    ## ========================================================================

    if (nrow(result$summary) > 0) {
        ## Name max.basins
        max.labels <- result$summary$label[result$summary$type == "max"]
        if (length(max.labels) == length(result$max.basins)) {
            names(result$max.basins) <- max.labels
        }

        ## Name min.basins
        min.labels <- result$summary$label[result$summary$type == "min"]
        if (length(min.labels) == length(result$min.basins)) {
            names(result$min.basins) <- min.labels
        }
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
        hop.k = hop.k
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

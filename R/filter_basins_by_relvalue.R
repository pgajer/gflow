#' Filter Basins by Relative Value Threshold
#'
#' @description
#' Filters basins of attraction by removing extrema whose relative function
#' values fall below specified thresholds. This operation proves essential when
#' analyzing biological or physical data where only extrema of sufficient
#' magnitude represent meaningful features rather than noise or artifacts.
#'
#' @details
#' Consider the problem of identifying biologically significant features in a
#' function landscape. Many local extrema emerge from measurement noise, minor
#' perturbations, or fine-scale variations that, while mathematically valid as
#' critical points, lack substantive meaning for downstream interpretation. The
#' relative value of an extremum, computed as its function value divided by the
#' mean across all vertices, provides a natural scale-invariant measure of
#' prominence.
#'
#' We begin by examining which extrema satisfy the threshold criteria. For local
#' maxima, we retain only those whose relative values exceed a specified minimum,
#' typically chosen to represent peaks that rise substantially above the average
#' landscape level. For local minima, we retain those falling below a maximum
#' relative value, capturing valleys that descend significantly beneath the mean.
#' These criteria filter the topological structure to emphasize features of
#' genuine interest.
#'
#' The filtering operation must preserve the integrity of the basin structure.
#' Simply removing extrema from the extrema list would leave orphaned basins in
#' the basin lists, creating inconsistencies that corrupt downstream analyses.
#' Therefore, we simultaneously filter both the extrema vertices and their
#' corresponding basin structures, ensuring that every basin in the output
#' corresponds to a retained extremum and vice versa.
#'
#' After filtering, the basin indices no longer match the original vertex indices.
#' The function reconstructs the complete basin structure with consistent indexing,
#' preserving all basin information (vertex sets, hop distances, boundary vertices)
#' for the retained extrema. This filtered structure integrates seamlessly with
#' clustering and merging operations, as the subsequent analyses depend only on
#' the relative relationships between basins rather than their absolute positions
#' in the original extrema list.
#'
#' The operation proves particularly valuable in biological applications where
#' the signal-to-noise ratio varies across features. Microbiome data, for instance,
#' often exhibits numerous small local maxima in composition space that reflect
#' sampling variation rather than distinct ecological states. Filtering by relative
#' value focuses the analysis on robust features that represent stable configurations
#' worthy of biological interpretation.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}} containing the complete basin
#'   structure before filtering.
#'
#' @param min.rel.value.max Numeric value specifying the minimum relative value
#'   threshold for local maxima. Only maxima with relative values (value divided
#'   by mean of y) greater than or equal to this threshold will be retained.
#'   Default is \code{NULL}, meaning no filtering is applied to maxima. Set to
#'   \code{1.0} to retain only maxima above the mean, or higher values for more
#'   stringent filtering.
#'
#' @param max.rel.value.min Numeric value specifying the maximum relative value
#'   threshold for local minima. Only minima with relative values less than or
#'   equal to this threshold will be retained. Default is \code{NULL}, meaning
#'   no filtering is applied to minima. Set to \code{1.0} to retain only minima
#'   below the mean, or lower values for more stringent filtering.
#'
#' @return An object of class \code{"basins_of_attraction"} with the same structure
#'   as the input but containing only basins for extrema that satisfy the relative
#'   value criteria. The object contains:
#'   \describe{
#'     \item{lmin_basins}{List of basin structures for retained local minima,
#'       with each basin preserving its complete vertex set, hop distances, and
#'       boundary information.}
#'     \item{lmax_basins}{List of basin structures for retained local maxima,
#'       similarly preserving all basin properties.}
#'     \item{n_vertices}{Total number of vertices in the graph (unchanged).}
#'     \item{y}{Copy of the input function values (unchanged).}
#'     \item{filter.info}{List documenting the filtering operation with components:
#'       \describe{
#'         \item{min.rel.value.max}{Threshold used for maxima filtering, or
#'           \code{NULL} if not applied.}
#'         \item{max.rel.value.min}{Threshold used for minima filtering, or
#'           \code{NULL} if not applied.}
#'         \item{n.lmax.original}{Number of maxima before filtering.}
#'         \item{n.lmax.retained}{Number of maxima after filtering.}
#'         \item{n.lmin.original}{Number of minima before filtering.}
#'         \item{n.lmin.retained}{Number of minima after filtering.}
#'         \item{removed.lmax.vertices}{Vertex indices of removed maxima.}
#'         \item{removed.lmin.vertices}{Vertex indices of removed minima.}
#'       }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins of attraction
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Generate summary to examine relative values
#' edgelen.list <- compute.edge.lengths(adj.list, weight.list)
#' basin.summary <- summary(basins, edgelen.list)
#'
#' # Examine distribution of relative values
#' max.basins <- basin.summary[basin.summary$type == "max", ]
#' cat("Maxima relative values:\n")
#' print(sort(max.basins$rel.value, decreasing = TRUE))
#'
#' # Filter to retain only biologically significant maxima (rel.value >= 1.0)
#' filtered.basins <- filter.basins.by.relvalue(basins,
#'                                              min.rel.value.max = 1.0)
#'
#' cat("Original maxima:", length(basins$lmax_basins), "\n")
#' cat("Retained maxima:", length(filtered.basins$lmax_basins), "\n")
#'
#' # Examine what was removed
#' print(filtered.basins$filter.info)
#'
#' # Filter both maxima and minima
#' filtered.both <- filter.basins.by.relvalue(basins,
#'                                           min.rel.value.max = 1.2,
#'                                           max.rel.value.min = 0.8)
#'
#' # Now cluster only the retained extrema
#' max.clusters <- cluster.local.extrema(filtered.basins,
#'                                       edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#'
#' # Continue with merging
#' merged.max <- merge.clustered.basins(filtered.basins,
#'                                     max.clusters,
#'                                     extrema.type = "max")
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing basins,
#' \code{\link{summary.basins_of_attraction}} for computing relative values,
#' \code{\link{cluster.local.extrema}} for clustering filtered basins,
#' \code{\link{merge.clustered.basins}} for merging clustered basins
#'
#' @export
filter.basins.by.relvalue <- function(basins.obj,
                                      min.rel.value.max = NULL,
                                      max.rel.value.min = NULL) {
  
  ## Input validation
  if (!inherits(basins.obj, "basins_of_attraction")) {
    stop("basins.obj must be of class 'basins_of_attraction'")
  }
  
  if (!is.null(min.rel.value.max) && !is.numeric(min.rel.value.max)) {
    stop("min.rel.value.max must be numeric or NULL")
  }
  
  if (!is.null(max.rel.value.min) && !is.numeric(max.rel.value.min)) {
    stop("max.rel.value.min must be numeric or NULL")
  }
  
  ## Extract function values to compute relative values
  y <- basins.obj$y
  y.mean <- mean(y)
  
  ## Store original counts
  n.lmax.original <- length(basins.obj$lmax_basins)
  n.lmin.original <- length(basins.obj$lmin_basins)
  
  ## Initialize lists for filtered basins
  filtered.lmax.basins <- list()
  filtered.lmin.basins <- list()
  
  ## Track removed vertices
  removed.lmax.vertices <- integer(0)
  removed.lmin.vertices <- integer(0)
  
  ## Filter local maxima
  if (length(basins.obj$lmax_basins) > 0) {
    for (i in seq_along(basins.obj$lmax_basins)) {
      basin <- basins.obj$lmax_basins[[i]]
      rel.value <- basin$value / y.mean
      
      ## Check if this maximum should be retained
      retain <- TRUE
      if (!is.null(min.rel.value.max)) {
        retain <- rel.value >= min.rel.value.max
      }
      
      if (retain) {
        ## Add to filtered list
        filtered.lmax.basins <- c(filtered.lmax.basins, list(basin))
      } else {
        ## Track removed vertex
        removed.lmax.vertices <- c(removed.lmax.vertices, basin$vertex)
      }
    }
  }
  
  ## Filter local minima
  if (length(basins.obj$lmin_basins) > 0) {
    for (i in seq_along(basins.obj$lmin_basins)) {
      basin <- basins.obj$lmin_basins[[i]]
      rel.value <- basin$value / y.mean
      
      ## Check if this minimum should be retained
      retain <- TRUE
      if (!is.null(max.rel.value.min)) {
        retain <- rel.value <= max.rel.value.min
      }
      
      if (retain) {
        ## Add to filtered list
        filtered.lmin.basins <- c(filtered.lmin.basins, list(basin))
      } else {
        ## Track removed vertex
        removed.lmin.vertices <- c(removed.lmin.vertices, basin$vertex)
      }
    }
  }
  
  ## Count retained basins
  n.lmax.retained <- length(filtered.lmax.basins)
  n.lmin.retained <- length(filtered.lmin.basins)
  
  ## Create result object preserving original structure
  result <- list(
    lmin_basins = filtered.lmin.basins,
    lmax_basins = filtered.lmax.basins,
    n_vertices = basins.obj$n_vertices,
    y = basins.obj$y
  )
  
  ## Add filtering information
  result$filter.info <- list(
    min.rel.value.max = min.rel.value.max,
    max.rel.value.min = max.rel.value.min,
    n.lmax.original = n.lmax.original,
    n.lmax.retained = n.lmax.retained,
    n.lmax.removed = n.lmax.original - n.lmax.retained,
    n.lmin.original = n.lmin.original,
    n.lmin.retained = n.lmin.retained,
    n.lmin.removed = n.lmin.original - n.lmin.retained,
    removed.lmax.vertices = removed.lmax.vertices,
    removed.lmin.vertices = removed.lmin.vertices
  )
  
  ## Preserve class
  class(result) <- "basins_of_attraction"
  
  return(result)
}


#' Print Method for Filtered Basins Filter Info
#'
#' @description
#' Prints a summary of basin filtering information.
#'
#' @param filter.info List containing filter information from
#'   \code{\link{filter.basins.by.relvalue}}.
#'
#' @return Invisibly returns the filter.info object.
#'
#' @keywords internal
print.basin.filter.info <- function(filter.info) {
  cat("Basin Filtering Summary\n")
  cat("=======================\n\n")
  
  if (!is.null(filter.info$min.rel.value.max)) {
    cat("Local Maxima Filtering:\n")
    cat(sprintf("  Threshold: rel.value >= %.3f\n",
                filter.info$min.rel.value.max))
    cat(sprintf("  Original count: %d\n", filter.info$n.lmax.original))
    cat(sprintf("  Retained count: %d\n", filter.info$n.lmax.retained))
    cat(sprintf("  Removed count:  %d\n", filter.info$n.lmax.removed))
    
    if (filter.info$n.lmax.removed > 0) {
      cat(sprintf("  Removed vertices: %s\n",
                 paste(filter.info$removed.lmax.vertices, collapse = ", ")))
    }
    cat("\n")
  } else {
    cat("Local Maxima: No filtering applied\n\n")
  }
  
  if (!is.null(filter.info$max.rel.value.min)) {
    cat("Local Minima Filtering:\n")
    cat(sprintf("  Threshold: rel.value <= %.3f\n",
                filter.info$max.rel.value.min))
    cat(sprintf("  Original count: %d\n", filter.info$n.lmin.original))
    cat(sprintf("  Retained count: %d\n", filter.info$n.lmin.retained))
    cat(sprintf("  Removed count:  %d\n", filter.info$n.lmin.removed))
    
    if (filter.info$n.lmin.removed > 0) {
      cat(sprintf("  Removed vertices: %s\n",
                 paste(filter.info$removed.lmin.vertices, collapse = ", ")))
    }
    cat("\n")
  } else {
    cat("Local Minima: No filtering applied\n\n")
  }
  
  invisible(filter.info)
}

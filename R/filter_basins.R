#' Filter Basins by Retaining Only Specified Extrema
#'
#' @description
#' Filters a basin structure to retain only basins corresponding to a specified
#' set of extrema vertices. This function proves essential when extrema have been
#' filtered based on various criteria (such as isolation measures, basin size,
#' or other quality metrics) and the basin structure must be updated to reflect
#' only the retained features.
#'
#' @details
#' The filtering problem arises naturally in topological data analysis when
#' initial basin computation identifies numerous local extrema, many of which
#' fail subsequent quality checks. We begin with a complete basin structure
#' containing all detected extrema, then apply domain-specific filters to
#' identify which extrema represent genuine features worthy of further analysis.
#'
#' Consider the case of filtering isolated extrema using density-based measures.
#' The `p.mean.nbrs.dist` statistic quantifies how densely surrounded an extremum
#' is by other data points, with high values indicating isolation in sparse regions
#' of the feature space. By examining the distribution of this statistic across
#' extrema and identifying natural gaps, we can establish thresholds that separate
#' well-supported features from artifacts of sparse sampling. However, this filtering
#' operates at the summary level, leaving us with a vector of vertices corresponding
#' to extrema that passed our criteria.
#'
#' The function addresses the subsequent step of updating the basin structure itself.
#' Simply discarding rows from a summary table proves insufficient because the
#' basins object maintains complex nested structures linking extrema to their
#' associated vertex sets, hop distances, and boundary information. Each basin
#' in the `lmin_basins` and `lmax_basins` lists must be examined to determine
#' whether its extremum vertex appears in the retention set, with only matching
#' basins preserved in the output.
#'
#' This approach offers maximum flexibility regarding filtering criteria. Whether
#' extrema are selected based on isolation measures, basin sizes, relative values,
#' spatial location, or any combination of quality metrics, the function operates
#' uniformly on the resulting vertex set. The retention vector may contain only
#' minima, only maxima, or a mixed set, with the function automatically determining
#' the type of each extremum from the basin structure and filtering accordingly.
#'
#' The filtering operation preserves complete basin information for retained extrema.
#' All basin properties including vertex membership, hop distance maps, and boundary
#' vertices carry forward unchanged. The function also records comprehensive metadata
#' about the filtering operation, documenting which extrema were removed and providing
#' counts by extremum type. This transparency proves essential for understanding
#' how filtering decisions impact the subsequent analysis pipeline.
#'
#' After filtering, the resulting basin structure integrates seamlessly with
#' downstream operations. The filtered basins can proceed directly to gradient
#' flow graph construction, or they can undergo additional clustering and merging
#' if further consolidation seems warranted. The filtering step thus serves as
#' a quality control gate, ensuring that only extrema meeting domain-specific
#' criteria propagate through the analysis workflow.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}} or \code{\link{merge.clustered.basins}}.
#'   This object contains the complete basin structure before filtering.
#'
#' @param good.extrema Integer vector of vertex indices (1-based) for extrema to
#'   retain. These are the extrema that survived external filtering based on
#'   quality metrics, isolation measures, or other criteria. The vector may
#'   contain only minima, only maxima, or both types mixed together. The function
#'   automatically determines which type each vertex represents by examining the
#'   basin structure.
#'
#' @return An object of class \code{"basins_of_attraction"} with the same structure
#'   as the input but containing only basins for extrema whose vertices appear in
#'   \code{good.extrema}. The object contains:
#'   \describe{
#'     \item{lmin_basins}{List of basin structures for retained local minima,
#'       preserving all basin properties (vertex sets, hop distances, boundaries).}
#'     \item{lmax_basins}{List of basin structures for retained local maxima,
#'       similarly preserving complete basin information.}
#'     \item{n_vertices}{Total number of vertices in the graph (unchanged).}
#'     \item{y}{Copy of the input function values (unchanged).}
#'     \item{filter.info}{List documenting the filtering operation with components:
#'       \describe{
#'         \item{n.extrema.input}{Total number of extrema in input (minima + maxima).}
#'         \item{n.extrema.retained}{Total number of extrema after filtering.}
#'         \item{n.extrema.removed}{Number of extrema removed by filtering.}
#'         \item{n.lmin.original}{Number of minima before filtering.}
#'         \item{n.lmin.retained}{Number of minima after filtering.}
#'         \item{n.lmin.removed}{Number of minima removed.}
#'         \item{n.lmax.original}{Number of maxima before filtering.}
#'         \item{n.lmax.retained}{Number of maxima after filtering.}
#'         \item{n.lmax.removed}{Number of maxima removed.}
#'         \item{retained.lmin.vertices}{Vertex indices of retained minima.}
#'         \item{retained.lmax.vertices}{Vertex indices of retained maxima.}
#'         \item{removed.lmin.vertices}{Vertex indices of removed minima.}
#'         \item{removed.lmax.vertices}{Vertex indices of removed maxima.}
#'       }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute and merge basins
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' edgelen.list <- compute.edge.lengths(adj.list, weight.list)
#'
#' # Cluster and merge extrema
#' max.clusters <- cluster.local.extrema(basins, edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#' merged.max <- merge.clustered.basins(basins, max.clusters, extrema.type = "max")
#'
#' min.clusters <- cluster.local.extrema(merged.max, edgelen.list,
#'                                       extrema.type = "min",
#'                                       overlap.threshold = 0.15)
#' merged.basins <- merge.clustered.basins(merged.max, min.clusters,
#'                                        extrema.type = "min")
#'
#' # Generate summary
#' merged.basin.df <- summary(merged.basins, edgelen.list)
#'
#' # Filter maxima by isolation measure
#' max.basin.df <- merged.basin.df[merged.basin.df$type == "max", ]
#' max.basin.df <- max.basin.df[max.basin.df$p.mean.nbrs.dist < 0.8, ]
#' good.max.vertices <- max.basin.df$vertex
#'
#' # Filter minima by isolation measure
#' min.basin.df <- merged.basin.df[merged.basin.df$type == "min", ]
#' min.basin.df <- min.basin.df[min.basin.df$p.mean.nbrs.dist < 0.75, ]
#' good.min.vertices <- min.basin.df$vertex
#'
#' # Combine filtered extrema
#' good.extrema <- c(good.min.vertices, good.max.vertices)
#'
#' # Filter basin structure
#' filtered.basins <- filter.basins(merged.basins, good.extrema)
#'
#' cat("Original:", length(merged.basins$lmin_basins), "minima,",
#'     length(merged.basins$lmax_basins), "maxima\n")
#' cat("Filtered:", length(filtered.basins$lmin_basins), "minima,",
#'     length(filtered.basins$lmax_basins), "maxima\n")
#'
#' # Examine filtering details
#' print(filtered.basins$filter.info)
#'
#' # Continue with gradient flow graph construction
#' gflow.graph <- construct.gflow.graph(filtered.basins)
#' plot(gflow.graph)
#'
#' # Alternative: Filter only maxima
#' filtered.max.only <- filter.basins(merged.basins, good.max.vertices)
#'
#' # Alternative: Filter by basin size
#' large.basins <- merged.basin.df[merged.basin.df$basin.size >= 100, ]
#' filtered.by.size <- filter.basins(merged.basins, large.basins$vertex)
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing basins,
#' \code{\link{merge.clustered.basins}} for merging basins,
#' \code{\link{summary.basins_of_attraction}} for generating basin summaries,
#' \code{\link{filter.basins.by.relvalue}} for filtering by relative value thresholds,
#' \code{\link{construct.gflow.graph}} for constructing gradient flow graphs
#'
#' @export
filter.basins <- function(basins.obj, good.extrema) {
  
  ## Input validation
  if (!inherits(basins.obj, "basins_of_attraction")) {
    stop("basins.obj must be of class 'basins_of_attraction'")
  }
  
  if (!is.numeric(good.extrema) && !is.integer(good.extrema)) {
    stop("good.extrema must be a numeric or integer vector of vertex indices")
  }
  
  if (length(good.extrema) == 0) {
    stop("good.extrema must contain at least one vertex index")
  }
  
  ## Convert to integer and ensure uniqueness
  good.extrema <- as.integer(unique(good.extrema))
  
  ## Store original counts
  n.lmin.original <- length(basins.obj$lmin_basins)
  n.lmax.original <- length(basins.obj$lmax_basins)
  n.extrema.original <- n.lmin.original + n.lmax.original
  
  ## Initialize lists for filtered basins
  filtered.lmin.basins <- list()
  filtered.lmax.basins <- list()
  
  ## Track retained and removed vertices
  retained.lmin.vertices <- integer(0)
  retained.lmax.vertices <- integer(0)
  removed.lmin.vertices <- integer(0)
  removed.lmax.vertices <- integer(0)
  
  ## Filter local minima
  if (n.lmin.original > 0) {
    for (i in seq_along(basins.obj$lmin_basins)) {
      basin <- basins.obj$lmin_basins[[i]]
      
      ## Check if this minimum's vertex is in the retention set
      if (basin$vertex %in% good.extrema) {
        ## Retain this basin
        filtered.lmin.basins <- c(filtered.lmin.basins, list(basin))
        retained.lmin.vertices <- c(retained.lmin.vertices, basin$vertex)
      } else {
        ## Track removed vertex
        removed.lmin.vertices <- c(removed.lmin.vertices, basin$vertex)
      }
    }
  }
  
  ## Filter local maxima
  if (n.lmax.original > 0) {
    for (i in seq_along(basins.obj$lmax_basins)) {
      basin <- basins.obj$lmax_basins[[i]]
      
      ## Check if this maximum's vertex is in the retention set
      if (basin$vertex %in% good.extrema) {
        ## Retain this basin
        filtered.lmax.basins <- c(filtered.lmax.basins, list(basin))
        retained.lmax.vertices <- c(retained.lmax.vertices, basin$vertex)
      } else {
        ## Track removed vertex
        removed.lmax.vertices <- c(removed.lmax.vertices, basin$vertex)
      }
    }
  }
  
  ## Count retained basins
  n.lmin.retained <- length(filtered.lmin.basins)
  n.lmax.retained <- length(filtered.lmax.basins)
  n.extrema.retained <- n.lmin.retained + n.lmax.retained
  
  ## Check that we retained at least some extrema
  if (n.extrema.retained == 0) {
    warning("No extrema retained after filtering. All input vertices may be invalid.")
  }
  
  ## Check for vertices in good.extrema that were not found
  found.vertices <- c(retained.lmin.vertices, retained.lmax.vertices)
  missing.vertices <- setdiff(good.extrema, found.vertices)
  if (length(missing.vertices) > 0) {
    warning(sprintf(
      "The following vertices in good.extrema were not found in basin structure: %s",
      paste(missing.vertices, collapse = ", ")
    ))
  }
  
  ## Create result object preserving original structure
  result <- list(
    lmin_basins = filtered.lmin.basins,
    lmax_basins = filtered.lmax.basins,
    n_vertices = basins.obj$n_vertices,
    y = basins.obj$y
  )
  
  ## Add filtering information
  result$filter.info <- list(
    method = "vertex_selection",
    n.extrema.input = n.extrema.original,
    n.extrema.retained = n.extrema.retained,
    n.extrema.removed = n.extrema.original - n.extrema.retained,
    n.lmin.original = n.lmin.original,
    n.lmin.retained = n.lmin.retained,
    n.lmin.removed = n.lmin.original - n.lmin.retained,
    n.lmax.original = n.lmax.original,
    n.lmax.retained = n.lmax.retained,
    n.lmax.removed = n.lmax.original - n.lmax.retained,
    retained.lmin.vertices = retained.lmin.vertices,
    retained.lmax.vertices = retained.lmax.vertices,
    removed.lmin.vertices = removed.lmin.vertices,
    removed.lmax.vertices = removed.lmax.vertices,
    n.good.extrema.input = length(good.extrema),
    n.vertices.not.found = length(missing.vertices),
    vertices.not.found = missing.vertices
  )
  
  ## Preserve class
  class(result) <- "basins_of_attraction"
  
  return(result)
}


#' Print Method for Basin Filter Info (Enhanced)
#'
#' @description
#' Prints a summary of basin filtering information, supporting both
#' relative value filtering and vertex selection filtering.
#'
#' @param filter.info List containing filter information.
#'
#' @return Invisibly returns the filter.info object.
#'
#' @keywords internal
print_basin_filter_info <- function(filter.info) {
  cat("Basin Filtering Summary\n")
  cat("=======================\n\n")
  
  ## Check filtering method
  method <- filter.info$method
  if (is.null(method)) {
    method <- "relative_value"  # Backward compatibility
  }
  
  if (method == "vertex_selection") {
    ## Vertex selection filtering
    cat("Filtering Method: Vertex Selection\n")
    cat(sprintf("  Input vertices to retain: %d\n",
                filter.info$n.good.extrema.input))
    
    if (filter.info$n.vertices.not.found > 0) {
      cat(sprintf("  Warning: %d vertices not found in basin structure\n",
                 filter.info$n.vertices.not.found))
    }
    
    cat("\nOverall Filtering Results:\n")
    cat(sprintf("  Total extrema (input):    %d\n", filter.info$n.extrema.input))
    cat(sprintf("  Total extrema (retained): %d\n", filter.info$n.extrema.retained))
    cat(sprintf("  Total extrema (removed):  %d (%.1f%%)\n",
               filter.info$n.extrema.removed,
               100 * filter.info$n.extrema.removed / filter.info$n.extrema.input))
    
    cat("\nLocal Minima:\n")
    cat(sprintf("  Original count: %d\n", filter.info$n.lmin.original))
    cat(sprintf("  Retained count: %d\n", filter.info$n.lmin.retained))
    cat(sprintf("  Removed count:  %d\n", filter.info$n.lmin.removed))
    
    if (filter.info$n.lmin.retained > 0) {
      cat(sprintf("  Retained vertices: %s\n",
                 paste(head(filter.info$retained.lmin.vertices, 10), collapse = ", ")))
      if (filter.info$n.lmin.retained > 10) {
        cat(sprintf("    ... and %d more\n", filter.info$n.lmin.retained - 10))
      }
    }
    
    if (filter.info$n.lmin.removed > 0) {
      cat(sprintf("  Removed vertices: %s\n",
                 paste(head(filter.info$removed.lmin.vertices, 10), collapse = ", ")))
      if (filter.info$n.lmin.removed > 10) {
        cat(sprintf("    ... and %d more\n", filter.info$n.lmin.removed - 10))
      }
    }
    
    cat("\nLocal Maxima:\n")
    cat(sprintf("  Original count: %d\n", filter.info$n.lmax.original))
    cat(sprintf("  Retained count: %d\n", filter.info$n.lmax.retained))
    cat(sprintf("  Removed count:  %d\n", filter.info$n.lmax.removed))
    
    if (filter.info$n.lmax.retained > 0) {
      cat(sprintf("  Retained vertices: %s\n",
                 paste(head(filter.info$retained.lmax.vertices, 10), collapse = ", ")))
      if (filter.info$n.lmax.retained > 10) {
        cat(sprintf("    ... and %d more\n", filter.info$n.lmax.retained - 10))
      }
    }
    
    if (filter.info$n.lmax.removed > 0) {
      cat(sprintf("  Removed vertices: %s\n",
                 paste(head(filter.info$removed.lmax.vertices, 10), collapse = ", ")))
      if (filter.info$n.lmax.removed > 10) {
        cat(sprintf("    ... and %d more\n", filter.info$n.lmax.removed - 10))
      }
    }
    
  } else {
    ## Relative value filtering (backward compatibility)
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
  }
  
  invisible(filter.info)
}

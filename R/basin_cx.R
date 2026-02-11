
#' Create a Gradient Flow Basin Complex on a Weighted Graph
#'
#' @description
#' Constructs a comprehensive topological analysis of a scalar function defined on a graph
#' through identification, clustering, and simplification of gradient flow basins. The function
#' identifies local extrema (minima and maxima), constructs monotonic basins around them,
#' clusters similar basins based on their overlap, and creates a cell complex representing
#' their relationships.
#'
#' @details
#' This function performs the following key operations:
#' \enumerate{
#'   \item Detects local extrema (minima and maxima) on the graph
#'   \item Constructs monotonic basins around each extremum
#'   \item Clusters similar basins based on their vertex overlap
#'   \item Merges basins within clusters to simplify the analysis
#'   \item Constructs a cell complex from intersections between basins
#' }
#'
#' Each basin represents a region of the graph where the function values change
#' monotonically (either increasing or decreasing) from the local extremum. Basins
#' are labeled as "mx" for minima and "Mx" for maxima, where x is a sequential number.
#'
#' The cell complex represents intersections between basins and provides insights
#' into the topological structure of the function on the graph. Three types of cells
#' are created:
#' \itemize{
#'   \item Ascending-descending cells: Intersections between minimum and maximum basins
#'   \item Ascending-ascending cells: Intersections between two minimum basins
#'   \item Descending-descending cells: Intersections between two maximum basins
#' }
#'
#' @param adj.list List of integer vectors. Each vector contains indices of vertices
#'        adjacent to the corresponding vertex. Indices should be 1-based.
#' @param weight.list List of numeric vectors. Each vector contains weights of edges
#'        corresponding to adjacencies in adj.list.
#' @param y Numeric vector of length \eqn{n}, the response values at each graph vertex.
#' @param basin.merge.overlap.thld Numeric value between 0 and 1 specifying the threshold
#'        for basin overlap distance below which basins should be merged (default: 0.1).
#'        Lower values result in more conservative merging.
#' @param min.asc.desc.cell.size.thld Minimum size threshold for ascending-descending cells
#'        (default: 1). Cells with fewer vertices are excluded.
#' @param min.asc.asc.cell.size.thld Minimum size threshold for ascending-ascending cells
#'        (default: 1). Cells with fewer vertices are excluded.
#' @param min.desc.desc.cell.size.thld Minimum size threshold for descending-descending cells
#'        (default: 1). Cells with fewer vertices are excluded.
#' @param graph.params List of graph parameters.
#'
#' @return An object of class "basin_cx" containing:
#' \describe{
#'   \item{basins}{A list with elements \code{ascending} (minima) and \code{descending} (maxima),
#'        each containing the merged basins following clustering.}
#'   \item{local_extrema}{A data frame with vertex indices, types (maximum/minimum), labels,
#'        and function values for all extrema.}
#'   \item{original_y}{The original input function values.}
#'   \item{harmonic_predictions}{Numeric vector of repaired function values after simplification.}
#'   \item{cluster_assignments}{A named vector showing which cluster each original basin was assigned to.}
#'   \item{cluster_mappings}{A data frame showing the mapping between original basin labels and
#'        merged basin labels.}
#'   \item{initial_basin_cx}{The original basin complex before clustering and merging.}
#'   \item{cells}{A cell complex structure with ascending-descending, ascending-ascending,
#'        and descending-descending cells representing basin intersections.}
#'   \item{basins_df}{A data frame summarizing all basins with properties like extremum vertex,
#'        function value, relative span, and size.}
#'   \item{lmin_basins_df, lmax_basins_df}{Data frames containing only minima or maxima basins.}
#'   \item{*_dist_mat, *_overlap_dists}{Various distance and similarity matrices between basins.}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a graph with adjacency list and weights
#' adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
#' weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
#'
#' # Define a function on vertices
#' y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
#'
#' # Create basin complex
#' basin_cx <- create.basin.cx(adj_list, weight_list, y, basin.merge.overlap.thld = 0.15)
#'
#' # Examine results
#' summary(basin_cx)
#'
#' # Visualize original vs. simplified function
#' plot(basin_cx, type = "comparison")
#'
#' # Extract vertices from a specific basin
#' m1_vertices <- get_basin_vertices(basin_cx, "m1")
#' }
#'
#' @seealso \code{\link{plot.basin_cx}}
#'
#' @export
create.basin.cx <- function(adj.list,
                           weight.list,
                           y,
                           basin.merge.overlap.thld = 0.1,
                           min.asc.desc.cell.size.thld = 1,
                           min.asc.asc.cell.size.thld = 1,
                           min.desc.desc.cell.size.thld = 1,
                           graph.params = list()
                           ) {
    ## Input validation
    if (!is.list(adj.list) || !is.list(weight.list))
        stop("adj.list and weight.list must be lists")

    if (length(adj.list) != length(weight.list))
        stop("adj.list and weight.list must have the same length")

    if (length(y) != length(adj.list))
        stop("y must have the same length as adj.list")

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector")
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function (initial basin detection)
    initial_basin_cx <- .Call("S_create_basin_cx",
                              adj.list.0based,
                              weight.list,
                              as.numeric(y))

    ## Store original y values for comparison
    initial_basin_cx$original_y <- y

    ## Build unified local_extrema table and label basins
    ## Extract minima and maxima
    mins <- lapply(initial_basin_cx$lmin_basins, function(b) list(vertex=b$vertex, value=b$value))
    maxs <- lapply(initial_basin_cx$lmax_basins, function(b) list(vertex=b$vertex, value=b$value))
    n_min <- length(mins)
    n_max <- length(maxs)

    ## Create extrema data frame
    local_extrema <- create.extrema.df(mins, maxs, n_min, n_max)
    initial_basin_cx$local_extrema <- local_extrema

    ## Attach labels to basin lists
    initial_basin_cx <- attach.basin.labels(initial_basin_cx, local_extrema)

    ## Create basins data frame
    initial_basin_cx <- create.basins.df(initial_basin_cx, y)

    ## Create basins structure
    initial_basin_cx$basins <- list(
        ascending  = initial_basin_cx$lmin_basins,
        descending = initial_basin_cx$lmax_basins
    )
    initial_basin_cx$lmin_basins <- initial_basin_cx$lmax_basins <- NULL

    ## Compute the Jaccard index matrix between all basins
    ## initial_basin_cx$basins_jaccard_mat <- compute_basins_jaccard_matrix(initial_basin_cx$basins, initial_basin_cx$local_extrema)

    ## Relabel distance matrices to match basin labels
    initial_basin_cx <- relabel_distance_matrices(initial_basin_cx, local_extrema)

    ## Get basin IDs in order for consistent presentation
    lmin.basin.ids <- sort(initial_basin_cx$basins_df$label[initial_basin_cx$basins_df$is_max == 0])
    lmax.basin.ids <- sort(initial_basin_cx$basins_df$label[initial_basin_cx$basins_df$is_max == 1])

    ## Reorder basins for consistency
    initial_basin_cx <- basin_dataframes_reorder(initial_basin_cx, lmin.basin.ids, lmax.basin.ids)

    ## Extract vertices for all basins for further analysis
    lmin.basins.vertices.list <- extract_basins_vertices(initial_basin_cx, lmin.basin.ids)
    lmax.basins.vertices.list <- extract_basins_vertices(initial_basin_cx, lmax.basin.ids)

    ## Calculate various distance matrices for basin analysis
    initial_basin_cx <- calculate_basin_distances(
        initial_basin_cx,
        lmin.basins.vertices.list,
        lmax.basins.vertices.list,
        lmin.basin.ids,
        lmax.basin.ids
    )

    ## SIMPLIFICATION: Cluster basins based on overlap distances
    cluster_assignments <- cluster_basins_by_overlap(
        initial_basin_cx,
        lmin.basin.ids,
        lmax.basin.ids,
        basin.merge.overlap.thld
    )

    ## Merge basins within clusters
    basin_cx <- merge_basin_clusters(initial_basin_cx, cluster_assignments)

    ## Add original basin complex for reference
    basin_cx$initial_basin_cx <- initial_basin_cx

    ## Create cell complex for the merged basins
    basin_cx$cells <- create_basin_cell_complex(
        basin_cx,
        min.asc.desc.cell.size.thld,
        min.asc.asc.cell.size.thld,
        min.desc.desc.cell.size.thld
    )

    ## Use default parameters or override with user-provided ones
    default_params <- list(
        min_intersection = 1,
        ms_edges_only = FALSE,
        weight_type = "dice"
    )

    ## Merge with user-provided parameters (user params take precedence)
    params <- utils::modifyList(default_params, graph.params)

    class(basin_cx) <- "basin_cx"

    ## Create the graph with the specified parameters
    basin_cx$graph <- do.call(
        construct.basin.cx.graph,
        c(list(basin_cx = basin_cx), params)
    )

    return(basin_cx)
}

#' Print Method for Basin Complex Objects
#'
#' @description
#' Prints a concise summary of a gradient flow basin complex object, displaying
#' key statistics about merged minima/maxima, function value ranges, and cell
#' complex structure.
#'
#' @param x An object of class "basin_cx" as created by \code{\link{create.basin.cx}}.
#' @param ... Additional arguments passed to print methods (currently unused).
#'
#' @return
#' Invisibly returns the input object \code{x}. This function is called for its
#' side effect of printing a summary to the console.
#'
#' @details
#' The print method displays:
#' \itemize{
#'   \item Number of merged minima (ascending basins)
#'   \item Number of merged maxima (descending basins)
#'   \item Range of original function values
#'   \item Cell complex summary including total number of cells (if available)
#' }
#'
#' For more detailed information about the basin complex, use \code{summary()}.
#'
#' @examples
#' \dontrun{
#' # Create a basin complex
#' adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
#' weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
#' y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
#' basin_cx <- create.basin.cx(adj_list, weight_list, y)
#'
#' # Print basic information
#' print(basin_cx)
#' # or simply:
#' basin_cx
#' }
#'
#' @seealso
#' \code{\link{create.basin.cx}} for creating basin complex objects,
#' \code{\link{summary.basin_cx}} for detailed basin complex information
#'
#' @export
print.basin_cx <- function(x, ...) {
    cat("Gradient Flow Basin Complex\n")
    cat("---------------------------\n")
    cat(sprintf("Number of merged minima: %d\n", length(x$basins$ascending)))
    cat(sprintf("Number of merged maxima: %d\n", length(x$basins$descending)))
    cat(sprintf("Original function values range: [%.4f, %.4f]\n",
                min(x$original_y), max(x$original_y)))

    ## Print cell complex summary if available
    if (!is.null(x$cells)) {
        cat("\nCell Complex Summary:\n")
        cat(sprintf("  Total cells: %d\n", sum(
            length(x$cells$asc_desc_cells),
            length(x$cells$asc_asc_cells),
            length(x$cells$desc_desc_cells)
        )))
    }

    cat("\nUse summary() for more detailed information.\n")
    invisible(x)
}

#' Summary Method for Basin Complex Objects
#'
#' @description
#' Generates a comprehensive summary of a gradient flow basin complex, including
#' detailed information about extrema, basins, basin overlaps, and cell complex
#' structure.
#'
#' @param object An object of class "basin_cx" as created by \code{\link{create.basin.cx}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' An object of class "summary.basin_cx" containing:
#' \describe{
#'   \item{n_extrema}{Total number of local extrema}
#'   \item{n_minima}{Number of local minima}
#'   \item{n_maxima}{Number of local maxima}
#'   \item{minima}{Data frame of local minima information}
#'   \item{maxima}{Data frame of local maxima information}
#'   \item{n_ascending_basins}{Number of merged ascending basins (minima)}
#'   \item{n_descending_basins}{Number of merged descending basins (maxima)}
#'   \item{ascending_basin_sizes}{Summary statistics of ascending basin sizes}
#'   \item{descending_basin_sizes}{Summary statistics of descending basin sizes}
#'   \item{minima_info}{Detailed information about each minimum basin including
#'         vertex, label, size, and function value}
#'   \item{maxima_info}{Detailed information about each maximum basin including
#'         vertex, label, size, and function value}
#'   \item{ascending_basins_jaccard_index}{Matrix of Jaccard indices between
#'         pairs of ascending basins}
#'   \item{descending_basins_jaccard_index}{Matrix of Jaccard indices between
#'         pairs of descending basins}
#'   \item{cells_summary}{Summary of cell complex structure}
#'   \item{messages}{Any messages generated during basin complex creation}
#' }
#'
#' @details
#' The summary method computes several key statistics and relationships:
#'
#' \strong{Extrema Analysis:} Counts and details of all local minima and maxima
#' in the original function.
#'
#' \strong{Basin Statistics:} For both ascending (minimum) and descending (maximum)
#' basins, provides counts, size distributions, and detailed information including
#' the representative vertex, label, size (number of vertices), and function value.
#'
#' \strong{Basin Overlap:} Computes Jaccard indices between all pairs of basins
#' of the same type. The Jaccard index measures the overlap between two basins
#' as the ratio of their intersection size to their union size, ranging from 0
#' (no overlap) to 1 (complete overlap).
#'
#' \strong{Cell Complex:} Summarizes the structure of cells representing
#' intersections between basins.
#'
#' @examples
#' \dontrun{
#' # Create a basin complex
#' adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
#' weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
#' y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
#' basin_cx <- create.basin.cx(adj_list, weight_list, y)
#'
#' # Get detailed summary
#' basin_summary <- summary(basin_cx)
#'
#' # Access specific components
#' basin_summary$minima_info
#' basin_summary$ascending_basins_jaccard_index
#' }
#'
#' @seealso
#' \code{\link{create.basin.cx}} for creating basin complex objects,
#' \code{\link{print.basin_cx}} for basic basin complex information,
#' \code{\link{print.summary.basin_cx}} for printing summary objects
#'
#' @export
summary.basin_cx <- function(object, ...) {
  if (!inherits(object, "basin_cx")) {
    stop("Object must be of class 'basin_cx'")
  }

  extrema <- object$local_extrema
  n_extrema <- nrow(extrema)
  n_maxima <- sum(extrema$is_maximum)
  n_minima <- n_extrema - n_maxima

  # Maxima and minima tables
  max_indices <- which(extrema$is_maximum == 1)
  min_indices <- which(extrema$is_maximum == 0)

  maxima_df <- if (length(max_indices) > 0) extrema[max_indices, ] else data.frame()
  minima_df <- if (length(min_indices) > 0) extrema[min_indices, ] else data.frame()

  # Basin sizes
  asc_sizes <- vapply(object$basins$ascending, function(b) length(b$vertices), integer(1))
  desc_sizes <- vapply(object$basins$descending, function(b) length(b$vertices), integer(1))

  # Minima/maxima info tables with function values and sizes
  minima_info <- if (length(asc_sizes) > 0) {
    data.frame(
      vertex = vapply(object$basins$ascending, function(b) b$vertex, integer(1)),
      label = names(object$basins$ascending),
      size = asc_sizes,
      fn_value = vapply(object$basins$ascending, function(b) b$value, numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  maxima_info <- if (length(desc_sizes) > 0) {
    data.frame(
      vertex = vapply(object$basins$descending, function(b) b$vertex, integer(1)),
      label = names(object$basins$descending),
      size = desc_sizes,
      fn_value = vapply(object$basins$descending, function(b) b$value, numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  # Jaccard index matrices
  jaccard_index <- function(basins) {
    n <- length(basins)
    if (n <= 1) return(NULL)
    mat <- matrix(0, n, n)
    colnames(mat) <- rownames(mat) <- names(basins)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        vi <- basins[[i]]$vertices
        vj <- basins[[j]]$vertices
        int_len <- length(intersect(vi, vj))
        union_len <- length(union(vi, vj))
        jaccard <- if (union_len > 0) int_len / union_len else 0
        mat[i, j] <- mat[j, i] <- jaccard
      }
    }
    mat
  }

  asc_jaccard <- jaccard_index(object$basins$ascending)
  desc_jaccard <- jaccard_index(object$basins$descending)

  result <- list(
    n_extrema = n_extrema,
    n_minima = n_minima,
    n_maxima = n_maxima,
    minima = minima_df,
    maxima = maxima_df,
    n_ascending_basins = length(object$basins$ascending),
    n_descending_basins = length(object$basins$descending),
    ascending_basin_sizes = if (length(asc_sizes) > 0) summary(asc_sizes) else NULL,
    descending_basin_sizes = if (length(desc_sizes) > 0) summary(desc_sizes) else NULL,
    minima_info = minima_info,
    maxima_info = maxima_info,
    ascending_basins_jaccard_index = asc_jaccard,
    descending_basins_jaccard_index = desc_jaccard,
    cells_summary = summarize_cells(object$cells),
    messages = if (!is.null(object$messages)) object$messages else NULL
  )

  class(result) <- "summary.basin_cx"
  return(result)
}


#' Print Method for Basin Complex Summary Objects
#'
#' @description
#' Prints a formatted, comprehensive summary of a gradient flow basin complex,
#' including detailed statistics about extrema, basins, basin overlaps, and
#' cell complex structure.
#'
#' @param x An object of class "summary.basin_cx" as created by
#'        \code{\link{summary.basin_cx}}.
#' @param ... Additional arguments passed to print methods (currently unused).
#'
#' @return
#' Invisibly returns the input summary object \code{x}. This function is called
#' for its side effect of printing detailed information to the console.
#'
#' @details
#' The print method displays the following sections:
#'
#' \strong{Extrema:} Total counts of local extrema, broken down by minima and maxima.
#'
#' \strong{Local Minima/Maxima Information:} Tables showing each extremum with its
#' vertex index, label, basin size, and function value. Minima are sorted by
#' increasing function value, maxima by decreasing function value.
#'
#' \strong{Basins:} Counts of merged ascending (minimum) and descending (maximum)
#' basins, along with size statistics.
#'
#' \strong{Basin Size Statistics:} Summary statistics (min, quartiles, mean, max)
#' for the number of vertices in each type of basin.
#'
#' \strong{Cell Complex Summary:} Total number of cells and breakdown by type
#' (ascending-descending, ascending-ascending, descending-descending).
#'
#' \strong{Jaccard Index Matrices:} Matrices showing pairwise Jaccard indices
#' (overlap measures) between basins of the same type, rounded to 3 decimal places.
#'
#' \strong{Messages:} Any informational messages generated during basin complex
#' creation.
#'
#' @examples
#' \dontrun{
#' # Create a basin complex
#' adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
#' weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
#' y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
#' basin_cx <- create.basin.cx(adj_list, weight_list, y)
#'
#' # Get and print detailed summary
#' basin_summary <- summary(basin_cx)
#' print(basin_summary)
#' # or simply:
#' basin_summary
#' }
#'
#' @seealso
#' \code{\link{summary.basin_cx}} for creating summary objects,
#' \code{\link{create.basin.cx}} for creating basin complex objects,
#' \code{\link{print.basin_cx}} for basic basin complex printing
#'
#' @export
print.summary.basin_cx <- function(x, ...) {
  cat("Gradient Flow Basin Complex Summary\n")
  cat("==================================\n\n")

  # Extrema
  cat("Extrema:\n")
  cat(sprintf("  Total: %d\n", x$n_extrema))
  cat(sprintf("  Local maxima: %d\n", x$n_maxima))
  cat(sprintf("  Local minima: %d\n\n", x$n_minima))

  # Minima info
  if (!is.null(x$minima_info) && nrow(x$minima_info) > 0) {
    cat("Local Minima Information:\n")
    print(x$minima_info[order(x$minima_info$fn_value), ])
    cat("\n")
  }

  # Maxima info
  if (!is.null(x$maxima_info) && nrow(x$maxima_info) > 0) {
    cat("Local Maxima Information:\n")
    print(x$maxima_info[order(x$maxima_info$fn_value, decreasing = TRUE), ])
    cat("\n")
  }

  # Basins
  cat("Basins:\n")
  cat(sprintf("  Ascending basins (minima): %d\n", x$n_ascending_basins))
  cat(sprintf("  Descending basins (maxima): %d\n\n", x$n_descending_basins))

  if (!is.null(x$ascending_basin_sizes)) {
    cat("Ascending Basin Size Statistics:\n")
    print(x$ascending_basin_sizes)
    cat("\n")
  }

  if (!is.null(x$descending_basin_sizes)) {
    cat("Descending Basin Size Statistics:\n")
    print(x$descending_basin_sizes)
    cat("\n")
  }

  # Cell Complex
  if (!is.null(x$cells_summary)) {
    cat("Cell Complex Summary:\n")
    cat(sprintf("  Total cells: %d\n", x$cells_summary$total_cells))
    cat(sprintf("  Ascending-Descending cells: %d\n", x$cells_summary$n_asc_desc_cells))
    cat(sprintf("  Ascending-Ascending cells: %d\n", x$cells_summary$n_asc_asc_cells))
    cat(sprintf("  Descending-Descending cells: %d\n\n", x$cells_summary$n_desc_desc_cells))
  }

  # Jaccard Index for Ascending
  if (!is.null(x$ascending_basins_jaccard_index)) {
    cat("Ascending Basins Jaccard Index:\n")
    print(round(x$ascending_basins_jaccard_index, 3))
    cat("\n")
  }

  # Jaccard Index for Descending
  if (!is.null(x$descending_basins_jaccard_index)) {
    cat("Descending Basins Jaccard Index:\n")
    print(round(x$descending_basins_jaccard_index, 3))
    cat("\n")
  }

  # Messages
  if (!is.null(x$messages)) {
    cat("Messages:\n")
    cat(x$messages, "\n")
  }

  invisible(x)
}

## Get cell vertices
get.cell.vertices <- function(x, cell_id) {
    if (!inherits(x, "basin_cx") || is.null(x$cells)) {
        stop("'x' must be a basin_cx object with cells")
    }

    ## Check each cell type
    cell_types <- c("asc_desc_cells", "asc_asc_cells", "desc_desc_cells")
    for (type in cell_types) {
        if (cell_id %in% names(x$cells[[type]])) {
            return(sort(x$cells[[type]][[cell_id]]$vertices))
        }
    }

    warning(paste("Cell with ID", cell_id, "not found"))
    return(integer(0))
}

#' Create Data Frame of Local Extrema
#'
#' @description
#' Internal function that creates a data frame containing information about local
#' extrema (minima and maxima) from lists of minimum and maximum basin data.
#' The function assigns standardized labels to extrema based on their function
#' values, with minima labeled in ascending order of value and maxima labeled
#' in descending order.
#'
#' @param mins List of minima, where each element contains at least 'vertex'
#'        (vertex index) and 'value' (function value) components.
#' @param maxs List of maxima, where each element contains at least 'vertex'
#'        (vertex index) and 'value' (function value) components.
#' @param n_min Integer specifying the number of minima in the mins list.
#' @param n_max Integer specifying the number of maxima in the maxs list.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{vertex_index}{Integer indices of extrema vertices}
#'   \item{is_maximum}{Integer flag (0 for minimum, 1 for maximum)}
#'   \item{label}{Character labels assigned to extrema ("m1", "m2", ... for minima;
#'         "M1", "M2", ... for maxima)}
#'   \item{fn_value}{Numeric function values at the extrema}
#' }
#'
#' @details
#' Labels are assigned based on function value ordering:
#' \itemize{
#'   \item Minima are labeled "m1", "m2", etc., where "m1" has the lowest function value
#'   \item Maxima are labeled "M1", "M2", etc., where "M1" has the highest function value
#' }
#'
#' If no extrema are found (n_min + n_max = 0), an empty data frame with the
#' appropriate column structure is returned.
#'
#' @keywords internal
#' @noRd
create.extrema.df <- function(mins, maxs, n_min, n_max) {
  ## Initialize with empty data frame if no extrema found
  if (n_min + n_max == 0) {
    return(data.frame(
      vertex_index = integer(),
      is_maximum = integer(),
      label = character(),
      fn_value = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  ## Extract data for data frame construction
  e_idx <- c(sapply(mins, `[[`, "vertex"), sapply(maxs, `[[`, "vertex"))
  is_max <- c(rep(FALSE, n_min), rep(TRUE, n_max))
  fn_val <- c(sapply(mins, `[[`, "value"), sapply(maxs, `[[`, "value"))

  ## Create labels for extrema
  labels <- character(length(e_idx))

  ## Label maxima in descending order of function value
  if (n_max > 0) {
    ord <- order(fn_val[(n_min+1):(n_min+n_max)], decreasing = TRUE)
    labels[(n_min+1):(n_min+n_max)][ord] <- paste0("M", seq_len(n_max))
  }

  ## Label minima in ascending order of function value
  if (n_min > 0) {
    ord <- order(fn_val[1:n_min], decreasing = FALSE)
    labels[1:n_min][ord] <- paste0("m", seq_len(n_min))
  }

  ## Create and return the data frame
  data.frame(
    vertex_index = e_idx,
    is_maximum = as.integer(is_max),
    label = labels,
    fn_value = fn_val,
    stringsAsFactors = FALSE
  )
}

#' Attach Labels to Basin Lists
#'
#' @description
#' Internal function that assigns labels from the local extrema data frame to
#' corresponding basins in the basin complex. Labels are matched based on the
#' extremum vertex index.
#'
#' @param basin_cx A basin complex object containing at least 'lmin_basins' and
#'        'lmax_basins' lists. Each basin in these lists must have a 'vertex'
#'        component indicating the extremum vertex index.
#' @param local_extrema A data frame containing extrema information with at least
#'        'vertex_index' and 'label' columns, as created by create_extrema_df().
#'
#' @return The modified basin_cx object with labels attached to each basin.
#'         Each basin in lmin_basins and lmax_basins will have a 'label'
#'         component added if a matching vertex is found in local_extrema.
#'
#' @details
#' The function iterates through all minima basins (lmin_basins) and maxima
#' basins (lmax_basins), finding the corresponding label from the local_extrema
#' data frame based on vertex index matching. If a matching vertex is found,
#' the label is added to the basin structure. If no match is found, the basin
#' remains unchanged.
#'
#' This function modifies the basin_cx object in place by adding label components
#' to individual basins, establishing the connection between the extrema table
#' and the basin structures.
#'
#' @keywords internal
#' @noRd
attach.basin.labels <- function(basin_cx, local_extrema) {
  ## Attach labels to minima basins
  for (i in seq_along(basin_cx$lmin_basins)) {
    vertex <- basin_cx$lmin_basins[[i]]$vertex
    label <- local_extrema$label[local_extrema$vertex_index == vertex]
    if (length(label) > 0) {
      basin_cx$lmin_basins[[i]]$label <- label
    }
  }

  ## Attach labels to maxima basins
  for (i in seq_along(basin_cx$lmax_basins)) {
    vertex <- basin_cx$lmax_basins[[i]]$vertex
    label <- local_extrema$label[local_extrema$vertex_index == vertex]
    if (length(label) > 0) {
      basin_cx$lmax_basins[[i]]$label <- label
    }
  }

  return(basin_cx)
}

#' Create Data Frame Representation of Basins
#'
#' @description
#' Internal function that converts basin matrix data into structured data frames,
#' creating separate representations for all basins, minima basins, and maxima basins.
#' The function processes raw basin matrix data and adds appropriate labels, types,
#' and formatting.
#'
#' @param basin_cx A basin complex object containing 'lmin_basins', 'lmax_basins',
#'        'basins_matrix', and 'local_extrema' components.
#' @param y Numeric vector of function values (currently unused but kept for
#'        potential future use or backwards compatibility).
#'
#' @return The modified basin_cx object with three data frame components added:
#' \describe{
#'   \item{basins_df}{Complete data frame of all basins with columns: evertex,
#'         hop_idx, value, is_max, min_rspan, max_rspan, Drspan, size, rsize, label}
#'   \item{lmin_basins_df}{Subset containing only minima basins, with labels as row names}
#'   \item{lmax_basins_df}{Subset containing only maxima basins, with labels as row names}
#' }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Checks if basins exist; returns empty structures if none found
#'   \item Extracts basin information from the basins_matrix
#'   \item Matches basins to their labels using local_extrema data
#'   \item Converts matrix to data frame with proper column types
#'   \item Rounds numeric values to 2 decimal places for readability
#'   \item Creates separate filtered data frames for minima and maxima
#' }
#'
#' Column definitions in the output data frames:
#' \itemize{
#'   \item evertex: Extremum vertex index
#'   \item hop_idx: Hop index from extremum
#'   \item value: Function value at extremum
#'   \item is_max: Binary indicator (0 = minimum, 1 = maximum)
#'   \item min_rspan: Minimum relative span
#'   \item max_rspan: Maximum relative span
#'   \item Drspan: Delta relative span
#'   \item size: Number of vertices in basin
#'   \item rsize: Relative size of basin
#'   \item label: Basin label (e.g., "m1", "M2")
#' }
#'
#' @keywords internal
#' @noRd
create.basins.df <- function(basin_cx, y) {
    ## Count total number of basins (minima + maxima)
    total_basins <- length(basin_cx$lmin_basins) + length(basin_cx$lmax_basins)

    if (total_basins == 0) {
        ## Return empty data frame if no basins
        basin_cx$basins_df <- data.frame(
            evertex = integer(),
            hop_idx = integer(),
            value = numeric(),
            is_max = integer(),
            min_rspan = numeric(),
            max_rspan = numeric(),
            Drspan = numeric(),
            size = integer(),
            rsize = numeric(),
            label = character(),
            stringsAsFactors = FALSE
        )
        return(basin_cx)
    }

    ## Create basin matrix
    basin_matrix <- basin_cx$basins_matrix

    ## Add labels as a new column if not already present
    if (!("label" %in% colnames(basin_matrix))) {
        basin_labels <- character(nrow(basin_matrix))

        ## Match each basin to its label in local_extrema
        for (i in 1:nrow(basin_matrix)) {
            vertex_idx <- basin_matrix[i, "extremum_vertex"]
            is_max <- as.logical(basin_matrix[i, "is_maximum"])

            ## Find matching row in local_extrema
            matching_row <- which(basin_cx$local_extrema$vertex_index == vertex_idx &
                                  basin_cx$local_extrema$is_maximum == is_max)

            if (length(matching_row) > 0) {
                basin_labels[i] <- basin_cx$local_extrema$label[matching_row]
            } else {
                basin_labels[i] <- NA_character_
            }
        }

        ## Add labels column to basin_matrix
        basin_matrix <- cbind(basin_matrix, label = basin_labels)
    }

    ## Convert matrix to data frame with appropriate column names and types
    basin_cx$basins_df <- as.data.frame(basin_matrix)
    basin_cx$basins_df$extremum_vertex <- as.integer(basin_cx$basins_df$extremum_vertex)
    basin_cx$basins_df$hop_idx <- as.integer(basin_cx$basins_df$hop_idx)
    basin_cx$basins_df$is_maximum <- as.integer(basin_cx$basins_df$is_maximum)
    basin_cx$basins_df$size <- as.integer(basin_cx$basins_df$size)
    basin_cx$basins_df$value <- round(as.numeric(basin_cx$basins_df$value), digits = 2)
    basin_cx$basins_df$rel_min_span <- round(as.numeric(basin_cx$basins_df$rel_min_span), digits = 2)
    basin_cx$basins_df$rel_max_span <- round(as.numeric(basin_cx$basins_df$rel_max_span), digits = 2)
    basin_cx$basins_df$delta_rel_span <- round(as.numeric(basin_cx$basins_df$delta_rel_span), digits = 2)
    basin_cx$basins_df$rel_size <- round(as.numeric(basin_cx$basins_df$rel_size), digits = 2)

    ## Rename columns for consistency
    colnames(basin_cx$basins_df) <- c("evertex", "hop_idx", "value", "is_max", "min_rspan",
                                      "max_rspan", "Drspan", "size", "rsize", "label")

    ## Create filtered dataframes for minima and maxima
    basin_cx$lmin_basins_df <- basin_cx$basins_df[basin_cx$basins_df$is_max == 0, ]
    rownames(basin_cx$lmin_basins_df) <- basin_cx$lmin_basins_df$label
    basin_cx$lmin_basins_df$label <- NULL

    basin_cx$lmax_basins_df <- basin_cx$basins_df[basin_cx$basins_df$is_max == 1, ]
    rownames(basin_cx$lmax_basins_df) <- basin_cx$lmax_basins_df$label
    basin_cx$lmax_basins_df$label <- NULL

    return(basin_cx)
}


#' Reorder Basin Data Frames for Consistent Presentation
#'
#' @description
#' Internal function that reorders the rows of minima and maxima basin data frames
#' to match a specified ordering of basin IDs. This ensures consistent ordering
#' across different representations of basin data.
#'
#' @param basin_cx A basin complex object containing 'lmin_basins_df' and
#'        'lmax_basins_df' data frames with basin labels as row names.
#' @param lmin_basin_ids Character vector specifying the desired order of minimum
#'        basin labels (e.g., c("m1", "m2", "m3")).
#' @param lmax_basin_ids Character vector specifying the desired order of maximum
#'        basin labels (e.g., c("M1", "M2", "M3")).
#'
#' @return The modified basin_cx object with reordered basin data frames.
#'         The row order of lmin_basins_df and lmax_basins_df will match the
#'         order specified in the ID vectors.
#'
#' @details
#' The function performs safe reordering by checking that:
#' \itemize{
#'   \item The ID vectors are not empty
#'   \item The corresponding data frames have rows to reorder
#' }
#'
#' If either condition is not met for a basin type, that data frame remains
#' unchanged. The function uses row name indexing to reorder the data frames,
#' so the basin labels must exist as row names in the respective data frames.
#'
#' This function is typically called after basin labeling to ensure that basins
#' are presented in a consistent order throughout the analysis pipeline.
#'
#' @keywords internal
#' @noRd
basin_dataframes_reorder <- function(basin_cx, lmin_basin_ids, lmax_basin_ids) {
  ## Reorder minima basins dataframe
  if (length(lmin_basin_ids) > 0 && nrow(basin_cx$lmin_basins_df) > 0) {
    basin_cx$lmin_basins_df <- basin_cx$lmin_basins_df[lmin_basin_ids, ]
  }

  ## Reorder maxima basins dataframe
  if (length(lmax_basin_ids) > 0 && nrow(basin_cx$lmax_basins_df) > 0) {
    basin_cx$lmax_basins_df <- basin_cx$lmax_basins_df[lmax_basin_ids, ]
  }

  return(basin_cx)
}

## 5. Gets vertices for each basin
extract_basins_vertices <- function(basin_cx, basin_ids) {
  result_list <- list()

  for (id in basin_ids) {
    ## Determine if minimum or maximum basin from prefix
    is_max <- substr(id, 1, 1) == "M"
    basin_list <- if (is_max) basin_cx$basins$descending else basin_cx$basins$ascending

    ## Find basin with matching label
    basin_idx <- NULL
    for (i in seq_along(basin_list)) {
      if (!is.null(basin_list[[i]]$label) && basin_list[[i]]$label == id) {
        basin_idx <- i
        break
      }
    }

    if (!is.null(basin_idx)) {
      basin <- basin_list[[basin_idx]]
      ## Extract vertices (first column of basin matrix)
      result_list[[id]] <- basin$basin[, 1]
    } else {
      result_list[[id]] <- integer(0)
    }
  }

  return(result_list)
}

## 6. Calculates various distance matrices
calculate_basin_distances <- function(basin_cx, lmin_basins_vertices_list, lmax_basins_vertices_list, lmin_basin_ids, lmax_basin_ids) {
  ## Calculate intersection sizes
  basin_cx$lmin_intersection_dists <- intersection_size_matrix(lmin_basins_vertices_list)
  basin_cx$lmax_intersection_dists <- intersection_size_matrix(lmax_basins_vertices_list)

  ## Calculate overlap distances
  basin_cx$lmin_overlap_dists <- compute.overlap.distance.matrix(lmin_basins_vertices_list)
  basin_cx$lmax_overlap_dists <- compute.overlap.distance.matrix(lmax_basins_vertices_list)

  ## Calculate relative distance matrices
  if (!is.null(basin_cx$lmin_dist_mat) && nrow(basin_cx$lmin_dist_mat) > 0) {
    basin_cx$lmin_rel_dist_mat <- basin_cx$lmin_dist_mat / basin_cx$graph_diameter
  }

  if (!is.null(basin_cx$lmax_dist_mat) && nrow(basin_cx$lmax_dist_mat) > 0) {
    basin_cx$lmax_rel_dist_mat <- basin_cx$lmax_dist_mat / basin_cx$graph_diameter
  }

  return(basin_cx)
}

## 7. Creates clusters based on overlap
cluster_basins_by_overlap <- function(basin_cx, lmin_basin_ids, lmax_basin_ids, basin_merge_overlap_thld) {
  ## Create overlap distance graphs
  lmax_overlap_dists_graph <- .create.threshold.distance.graph(
    basin_cx$lmax_overlap_dists,
    threshold = basin_merge_overlap_thld
  )

  lmin_overlap_dists_graph <- .create.threshold.distance.graph(
    basin_cx$lmin_overlap_dists,
    threshold = basin_merge_overlap_thld
  )

  ## Find connected components
  lmin_overlap_dists_graph_cc <- graph.connected.components(lmin_overlap_dists_graph$adj_list)
  names(lmin_overlap_dists_graph_cc) <- lmin_basin_ids

  lmax_overlap_dists_graph_cc <- graph.connected.components(lmax_overlap_dists_graph$adj_list)
  names(lmax_overlap_dists_graph_cc) <- lmax_basin_ids

  ## Combine minima and maxima clusters
  cluster_assignments <- c(lmin_overlap_dists_graph_cc, lmax_overlap_dists_graph_cc)

  return(cluster_assignments)
}

## 8. Merges basins within clusters
merge_basin_clusters <- function(basin_cx, cluster_assignments) {
  ## Separate min and max basin labels based on first letter
  min_basin_labels <- names(cluster_assignments)[substr(names(cluster_assignments), 1, 1) == "m"]
  max_basin_labels <- names(cluster_assignments)[substr(names(cluster_assignments), 1, 1) == "M"]

  ## Get cluster IDs for minima and maxima
  min_cluster_ids <- unique(cluster_assignments[min_basin_labels])
  max_cluster_ids <- unique(cluster_assignments[max_basin_labels])

  ## Initialize lists for merged basins (using temporary IDs for now)
  merged_lmin_basins_tmp <- list()
  merged_lmax_basins_tmp <- list()

  ## Get original function values
  original_y <- basin_cx$original_y

  ## Tracking basin mapping (original to merged)
  basin_mapping <- data.frame(
    original_label = character(),
    cluster_id = integer(),
    merged_label = character(),
    extremum_vertex = integer(),
    stringsAsFactors = FALSE
  )

  ## Process minima clusters
  min_extrema_info <- data.frame(
    cluster_id = integer(),
    extremum_vertex = integer(),
    value = numeric(),
    n_vertices = integer(),
    stringsAsFactors = FALSE
  )

  ## Process each minimum cluster
  for (cluster_id in min_cluster_ids) {
    ## Process minima cluster
    basin_labels <- min_basin_labels[cluster_assignments[min_basin_labels] == cluster_id]
    if (length(basin_labels) == 0) next

    ## Extract basin objects for this cluster
    cluster_basins <- list()
    for (label in basin_labels) {
      basin <- get_basin_by_label(basin_cx, label)
      if (!is.null(basin)) {
        cluster_basins[[label]] <- basin
      }
    }

    if (length(cluster_basins) == 0) next

    ## Find extremum with lowest value (for minima)
    extremum_values <- sapply(cluster_basins, `[[`, "value")
    extremum_label <- names(extremum_values)[which.min(extremum_values)]
    extremum_vertex <- cluster_basins[[extremum_label]]$vertex
    extremum_value <- min(extremum_values)

    ## Collect all vertices from all basins in the cluster
    all_vertices <- integer(0)
    for (basin in cluster_basins) {
      basin_vertices <- basin$basin[, 1]  ## First column has vertex indices
      all_vertices <- union(all_vertices, basin_vertices)
    }

    ## Create merged basin
    merged_basin <- list(
      vertex = extremum_vertex,
      value = extremum_value,
      vertices = sort(all_vertices)
    )

    ## Add to merged basins list (with temporary ID)
    merged_lmin_basins_tmp[[as.character(cluster_id)]] <- merged_basin

    ## Add to extrema info
    min_extrema_info <- rbind(min_extrema_info, data.frame(
      cluster_id = cluster_id,
      extremum_vertex = extremum_vertex,
      value = extremum_value,
      n_vertices = length(all_vertices),
      stringsAsFactors = FALSE
    ))

    ## Update basin mapping
    for (label in basin_labels) {
      basin_mapping <- rbind(basin_mapping, data.frame(
        original_label = label,
        cluster_id = cluster_id,
        merged_label = NA, ## Will be filled later
        extremum_vertex = extremum_vertex,
        stringsAsFactors = FALSE
      ))
    }
  }

  ## Process maxima clusters
  max_extrema_info <- data.frame(
    cluster_id = integer(),
    extremum_vertex = integer(),
    value = numeric(),
    n_vertices = integer(),
    stringsAsFactors = FALSE
  )

  ## Process each maximum cluster
  for (cluster_id in max_cluster_ids) {
    basin_labels <- max_basin_labels[cluster_assignments[max_basin_labels] == cluster_id]
    if (length(basin_labels) == 0) next

    ## Extract basin objects for this cluster
    cluster_basins <- list()
    for (label in basin_labels) {
      basin <- get_basin_by_label(basin_cx, label)
      if (!is.null(basin)) {
        cluster_basins[[label]] <- basin
      }
    }

    if (length(cluster_basins) == 0) next

    ## Find extremum with highest value (for maxima)
    extremum_values <- sapply(cluster_basins, `[[`, "value")
    extremum_label <- names(extremum_values)[which.max(extremum_values)]
    extremum_vertex <- cluster_basins[[extremum_label]]$vertex
    extremum_value <- max(extremum_values)

    ## Collect all vertices from all basins in the cluster
    all_vertices <- integer(0)
    for (basin in cluster_basins) {
      basin_vertices <- basin$basin[, 1]  ## First column has vertex indices
      all_vertices <- union(all_vertices, basin_vertices)
    }

    ## Create merged basin
    merged_basin <- list(
      vertex = extremum_vertex,
      value = extremum_value,
      vertices = sort(all_vertices)
    )

    ## Add to merged basins list (with temporary ID)
    merged_lmax_basins_tmp[[as.character(cluster_id)]] <- merged_basin

    ## Add to extrema info
    max_extrema_info <- rbind(max_extrema_info, data.frame(
      cluster_id = cluster_id,
      extremum_vertex = extremum_vertex,
      value = extremum_value,
      n_vertices = length(all_vertices),
      stringsAsFactors = FALSE
    ))

    ## Update basin mapping
    for (label in basin_labels) {
      basin_mapping <- rbind(basin_mapping, data.frame(
        original_label = label,
        cluster_id = cluster_id,
        merged_label = NA, ## Will be filled later
        extremum_vertex = extremum_vertex,
        stringsAsFactors = FALSE
      ))
    }
  }

  ## Create final basins with proper labels
  merged_lmin_basins <- list()
  merged_lmax_basins <- list()

  ## Sort minima by increasing value
  if (nrow(min_extrema_info) > 0) {
    min_extrema_info <- min_extrema_info[order(min_extrema_info$value), ]

    ## Create properly labeled minima basins
    for (i in 1:nrow(min_extrema_info)) {
      cluster_id <- min_extrema_info$cluster_id[i]
      new_label <- paste0("m", i)

      ## Copy the basin with the new label
      merged_lmin_basins[[new_label]] <- merged_lmin_basins_tmp[[as.character(cluster_id)]]
      merged_lmin_basins[[new_label]]$label <- new_label

      ## Update basin mapping
      basin_mapping$merged_label[basin_mapping$cluster_id == cluster_id &
                               substr(basin_mapping$original_label, 1, 1) == "m"] <- new_label
    }
  }

  ## Sort maxima by decreasing value
  if (nrow(max_extrema_info) > 0) {
    max_extrema_info <- max_extrema_info[order(max_extrema_info$value, decreasing = TRUE), ]

    ## Create properly labeled maxima basins
    for (i in 1:nrow(max_extrema_info)) {
      cluster_id <- max_extrema_info$cluster_id[i]
      new_label <- paste0("M", i)

      ## Copy the basin with the new label
      merged_lmax_basins[[new_label]] <- merged_lmax_basins_tmp[[as.character(cluster_id)]]
      merged_lmax_basins[[new_label]]$label <- new_label

      ## Update basin mapping
      basin_mapping$merged_label[basin_mapping$cluster_id == cluster_id &
                               substr(basin_mapping$original_label, 1, 1) == "M"] <- new_label
    }
  }

  ## Create local_extrema data frame for the merged basin complex
  local_extrema <- data.frame(
    vertex_index = integer(),
    is_maximum = integer(),
    label = character(),
    fn_value = numeric(),
    stringsAsFactors = FALSE
  )

  ## Add minima to local_extrema
  for (i in seq_along(merged_lmin_basins)) {
    basin <- merged_lmin_basins[[i]]
    local_extrema <- rbind(local_extrema, data.frame(
      vertex_index = basin$vertex,
      is_maximum = 0,
      label = names(merged_lmin_basins)[i],
      fn_value = basin$value,
      stringsAsFactors = FALSE
    ))
  }

  ## Add maxima to local_extrema
  for (i in seq_along(merged_lmax_basins)) {
    basin <- merged_lmax_basins[[i]]
    local_extrema <- rbind(local_extrema, data.frame(
      vertex_index = basin$vertex,
      is_maximum = 1,
      label = names(merged_lmax_basins)[i],
      fn_value = basin$value,
      stringsAsFactors = FALSE
    ))
  }

  ## Create the merged basin complex
  result <- list(
    basins = list(
      ascending = merged_lmin_basins,
      descending = merged_lmax_basins
    ),
    local_extrema = local_extrema,
    cluster_assignments = cluster_assignments,
    original_y = original_y,
    cluster_mappings = basin_mapping,
    harmonic_predictions = basin_cx$harmonic_predictions
  )

  return(result)
}



## Utility function to get a basin by its label
get_basin_by_label <- function(basin_cx, label) {
  ## Determine if minimum or maximum basin from the label prefix
  is_max <- substr(label, 1, 1) == "M"
  basin_list <- if (is_max) basin_cx$basins$descending else basin_cx$basins$ascending

  ## Find basin with matching label
  for (i in seq_along(basin_list)) {
    if (!is.null(basin_list[[i]]$label) && basin_list[[i]]$label == label) {
      return(basin_list[[i]])
    }
  }

  return(NULL)
}

## Calculate intersection size matrix between all pairs of basin vertex sets
intersection_size_matrix <- function(basin_vertices_list) {
  n_basins <- length(basin_vertices_list)
  if (n_basins == 0) return(matrix(0, 0, 0))

  basin_ids <- names(basin_vertices_list)
  result <- matrix(0, nrow = n_basins, ncol = n_basins)
  rownames(result) <- basin_ids
  colnames(result) <- basin_ids

  ## Calculate intersection sizes
  for (i in 1:n_basins) {
    vertices_i <- basin_vertices_list[[i]]
    result[i, i] <- length(vertices_i)  ## Self-intersection = size

    for (j in (i+1):n_basins) {
      if (j <= n_basins) {  ## Ensure we don't go out of bounds
        vertices_j <- basin_vertices_list[[j]]
        intersection_size <- length(intersect(vertices_i, vertices_j))

        ## Store the intersection size in both directions
        result[i, j] <- intersection_size
        result[j, i] <- intersection_size
      }
    }
  }

  return(result)
}


## Create a graph from a distance matrix with edges connecting vertices under a threshold
.create.threshold.distance.graph <- function(distance_matrix, threshold, include.names = TRUE) {

    n_vertices <- nrow(distance_matrix)
    if (n_vertices == 0) {
        return(list(
            adj_list = list(),
            weight_list = list()
        ))
    }

    ## Get vertex names if available
    vertex.names <- rownames(distance_matrix)
    if (is.null(vertex.names)) {
        vertex.names <- 1:n_vertices
    }

    adj_list <- vector("list", n_vertices)
    weight_list <- vector("list", n_vertices)

    ## Initialize empty lists
    ## for (i in 1:n_vertices) {
    ##     adj_list[[i]] <- integer(0)
    ##     weight_list[[i]] <- numeric(0)
    ## }

    ## Fill adjacency and weight lists
    ## for (i in 1:n_vertices) {
    ##   for (j in 1:n_vertices) {
    ##     if (i != j && distance_matrix[i, j] <= threshold) {
    ##       ## 1-based indices for R
    ##       adj_list[[i]] <- c(adj_list[[i]], j)
    ##       weight_list[[i]] <- c(weight_list[[i]], distance_matrix[i, j])
    ##     }
    ##   }
    ## }

    ## Populate the adjacency and weight lists
    for (i in 1:n_vertices) {
        ## Find all vertices j where dist(i,j) < threshold and i != j
        neighbors <- which(distance_matrix[i, ] < threshold & (1:n_vertices != i))

        ## Add to adjacency list
        adj_list[[i]] <- neighbors

        ## Add corresponding weights
        weight_list[[i]] <- distance_matrix[i, neighbors]
    }

    ## Add names if requested
    if (include.names && !is.null(vertex.names)) {
        names(adj_list) <- vertex.names
        names(weight_list) <- vertex.names
    }

    return(list(
        adj_list = adj_list,
        weight_list = weight_list
    ))
}

## Function to summarize cells for summary.basin_cx
summarize_cells <- function(cells) {
  if (is.null(cells)) return(NULL)

  n_asc_desc_cells <- length(cells$asc_desc_cells)
  n_asc_asc_cells <- length(cells$asc_asc_cells)
  n_desc_desc_cells <- length(cells$desc_desc_cells)
  total_cells <- n_asc_desc_cells + n_asc_asc_cells + n_desc_desc_cells

  ## Calculate cell size statistics
  asc_desc_sizes <- if (n_asc_desc_cells > 0) {
    sapply(cells$asc_desc_cells, function(cell) length(cell$vertices))
  } else {
    numeric(0)
  }

  asc_asc_sizes <- if (n_asc_asc_cells > 0) {
    sapply(cells$asc_asc_cells, function(cell) length(cell$vertices))
  } else {
    numeric(0)
  }

  desc_desc_sizes <- if (n_desc_desc_cells > 0) {
    sapply(cells$desc_desc_cells, function(cell) length(cell$vertices))
  } else {
    numeric(0)
  }

  ## Create size summaries
  asc_desc_summary <- if (length(asc_desc_sizes) > 0) summary(asc_desc_sizes) else NULL
  asc_asc_summary <- if (length(asc_asc_sizes) > 0) summary(asc_asc_sizes) else NULL
  desc_desc_summary <- if (length(desc_desc_sizes) > 0) summary(desc_desc_sizes) else NULL

  return(list(
    n_asc_desc_cells = n_asc_desc_cells,
    n_asc_asc_cells = n_asc_asc_cells,
    n_desc_desc_cells = n_desc_desc_cells,
    total_cells = total_cells,
    asc_desc_size_summary = asc_desc_summary,
    asc_asc_size_summary = asc_asc_summary,
    desc_desc_size_summary = desc_desc_summary
  ))
}

## Relabel distance matrices to match basin labels
relabel_distance_matrices <- function(basin_cx, local_extrema) {
  ## Relabel minimum distance matrix
  if (!is.null(basin_cx$lmin_dist_mat) && nrow(basin_cx$lmin_dist_mat) > 0) {
    ## Get the mapping of vertex indices to labels for minima
    min_vertex_to_label <- setNames(
      local_extrema$label[local_extrema$is_maximum == 0],
      as.character(local_extrema$vertex_index[local_extrema$is_maximum == 0])
    )

    ## Extract the vertex indices from row names
    rownames_vertices <- as.character(as.numeric(substring(rownames(basin_cx$lmin_dist_mat), 2)))
    colnames_vertices <- as.character(as.numeric(substring(colnames(basin_cx$lmin_dist_mat), 2)))

    ## Replace with corresponding labels
    rownames(basin_cx$lmin_dist_mat) <- min_vertex_to_label[rownames_vertices]
    colnames(basin_cx$lmin_dist_mat) <- min_vertex_to_label[colnames_vertices]
  }

  ## Relabel maximum distance matrix
  if (!is.null(basin_cx$lmax_dist_mat) && nrow(basin_cx$lmax_dist_mat) > 0) {
    ## Get the mapping of vertex indices to labels for maxima
    max_vertex_to_label <- setNames(
      local_extrema$label[local_extrema$is_maximum == 1],
      as.character(local_extrema$vertex_index[local_extrema$is_maximum == 1])
    )

    ## Extract the vertex indices from row names
    rownames_vertices <- as.character(as.numeric(substring(rownames(basin_cx$lmax_dist_mat), 2)))
    colnames_vertices <- as.character(as.numeric(substring(colnames(basin_cx$lmax_dist_mat), 2)))

    ## Replace with corresponding labels
    rownames(basin_cx$lmax_dist_mat) <- max_vertex_to_label[rownames_vertices]
    colnames(basin_cx$lmax_dist_mat) <- max_vertex_to_label[colnames_vertices]
  }

  return(basin_cx)
}



## Create basin cell complex
create_basin_cell_complex <- function(basin_cx,
                                     min.asc.desc.cell.size.thld = 1,
                                     min.asc.asc.cell.size.thld = 1,
                                     min.desc.desc.cell.size.thld = 1) {
  ## Extract basins
  asc_basins <- basin_cx$basins$ascending  ## minima
  desc_basins <- basin_cx$basins$descending  ## maxima

  ## Create empty lists for different types of cells
  asc_desc_cells <- list()
  asc_asc_cells <- list()
  desc_desc_cells <- list()

  ## 1. Create ascending-descending cells (classical Morse-Smale)
  if (length(asc_basins) > 0 && length(desc_basins) > 0) {
    cell_idx <- 1

    for (min_label in names(asc_basins)) {
      min_basin <- asc_basins[[min_label]]
      lmin_vertex <- min_basin$vertex
      lmin_value <- min_basin$value
      min_vertices <- min_basin$vertices

      for (max_label in names(desc_basins)) {
        max_basin <- desc_basins[[max_label]]
        lmax_vertex <- max_basin$vertex
        lmax_value <- max_basin$value
        max_vertices <- max_basin$vertices

        ## Find intersection
        cell_vertices <- intersect(min_vertices, max_vertices)

        ## Check if intersection meets size threshold
        if (length(cell_vertices) >= min.asc.desc.cell.size.thld) {
          cell <- list(
            lmin_vertex = lmin_vertex,
            lmax_vertex = lmax_vertex,
            lmin_value = lmin_value,
            lmax_value = lmax_value,
            lmin_label = min_label,
            lmax_label = max_label,
            vertices = sort(cell_vertices)
          )

          asc_desc_cells[[paste0("cell_", cell_idx)]] <- cell
          cell_idx <- cell_idx + 1
        }
      }
    }
  }

  ## 2. Create ascending-ascending cells (min-min intersections)
  if (length(asc_basins) >= 2) {
    cell_idx <- 1
    asc_labels <- names(asc_basins)

    for (i in 1:(length(asc_labels) - 1)) {
      min1_label <- asc_labels[i]
      min1_basin <- asc_basins[[min1_label]]
      lmin1_vertex <- min1_basin$vertex
      lmin1_value <- min1_basin$value
      min1_vertices <- min1_basin$vertices

      for (j in (i+1):length(asc_labels)) {
        min2_label <- asc_labels[j]
        min2_basin <- asc_basins[[min2_label]]
        lmin2_vertex <- min2_basin$vertex
        lmin2_value <- min2_basin$value
        min2_vertices <- min2_basin$vertices

        ## Find intersection
        cell_vertices <- intersect(min1_vertices, min2_vertices)

        ## Check if intersection meets size threshold
        if (length(cell_vertices) >= min.asc.asc.cell.size.thld) {
          cell <- list(
            lmin1_vertex = lmin1_vertex,
            lmin2_vertex = lmin2_vertex,
            lmin1_value = lmin1_value,
            lmin2_value = lmin2_value,
            lmin1_label = min1_label,
            lmin2_label = min2_label,
            vertices = sort(cell_vertices)
          )

          asc_asc_cells[[paste0("cell_", cell_idx)]] <- cell
          cell_idx <- cell_idx + 1
        }
      }
    }
  }

  ## 3. Create descending-descending cells (max-max intersections)
  if (length(desc_basins) >= 2) {
    cell_idx <- 1
    desc_labels <- names(desc_basins)

    for (i in 1:(length(desc_labels) - 1)) {
      max1_label <- desc_labels[i]
      max1_basin <- desc_basins[[max1_label]]
      lmax1_vertex <- max1_basin$vertex
      lmax1_value <- max1_basin$value
      max1_vertices <- max1_basin$vertices

      for (j in (i+1):length(desc_labels)) {
        max2_label <- desc_labels[j]
        max2_basin <- desc_basins[[max2_label]]
        lmax2_vertex <- max2_basin$vertex
        lmax2_value <- max2_basin$value
        max2_vertices <- max2_basin$vertices

        ## Find intersection
        cell_vertices <- intersect(max1_vertices, max2_vertices)

        ## Check if intersection meets size threshold
        if (length(cell_vertices) >= min.desc.desc.cell.size.thld) {
          cell <- list(
            lmax1_vertex = lmax1_vertex,
            lmax2_vertex = lmax2_vertex,
            lmax1_value = lmax1_value,
            lmax2_value = lmax2_value,
            lmax1_label = max1_label,
            lmax2_label = max2_label,
            vertices = sort(cell_vertices)
          )

          desc_desc_cells[[paste0("cell_", cell_idx)]] <- cell
          cell_idx <- cell_idx + 1
        }
      }
    }
  }

  ## Create result
  result <- list(
    asc_desc_cells = asc_desc_cells,
    asc_asc_cells = asc_asc_cells,
    desc_desc_cells = desc_desc_cells
  )

  ## Add summary information
  attr(result, "n_asc_desc_cells") <- length(asc_desc_cells)
  attr(result, "n_asc_asc_cells") <- length(asc_asc_cells)
  attr(result, "n_desc_desc_cells") <- length(desc_desc_cells)
  attr(result, "total_cells") <- length(asc_desc_cells) + length(asc_asc_cells) + length(desc_desc_cells)

  return(result)
}

#' Show Basin Details from a Basin Complex
#'
#' @description
#' Retrieves and returns a specific basin from a basin complex object by its label.
#' The function searches through both ascending (minima) and descending (maxima) basins.
#'
#' @param x A basin complex object returned by create.basin.cx().
#' @param ... Additional arguments passed to methods.
#' @export
show.basin <- function(x, ...) {
  UseMethod("show.basin")
}

#' Show Basin Details from a Basin Complex
#'
#' @description
#' Retrieves and returns a specific basin from a basin complex object by its label.
#' The function searches through both ascending (minima) and descending (maxima) basins.
#'
#' @param x A basin complex object returned by create.basin.cx().
#' @param basin.label Character string specifying the basin label to retrieve (e.g., "M1" or "m2").
#' @param ... Additional arguments (not used).
#'
#' @return A basin list corresponding to the requested label. If no basin with the specified
#'        label is found, the function returns NULL and prints a message.
#'
#' @details
#' Basin labels follow the naming convention of the basin complex: "Mx" for maxima (descending)
#' and "mx" for minima (ascending), where x is a sequential number assigned during the basin
#' complex construction.
#'
#' The returned basin list typically contains:
#' - vertex: 1-based index of the extremum vertex
#' - value: Function value at the extremum
#' - label: The basin label matching the requested label
#' - basin: Matrix with vertex indices and distances from the extremum
#' - basin_bd: Matrix of boundary vertices and their monotonicity spans
#'
#' @export
show.basin.basin_cx <- function(x, basin.label, ...) {
    ## Check if input is a basin_cx object
    if (!inherits(x, "basin_cx")) {
        stop("The input 'x' must be a basin_cx object.")
    }

    ## Check if basin.label is provided
    if (missing(basin.label)) {
        stop("basin.label must be specified")
    }

    ## Determine if we're looking for a minimum or maximum basin
    is_max <- substr(basin.label, 1, 1) == "M"

    ## Select the appropriate basin list to search
    basin_list <- if (is_max) x$basins$descending else x$basins$ascending

    ## Search for the basin with the matching label
    matching_basin <- NULL
    for (basin in basin_list) {
        if (!is.null(basin$label) && basin$label == basin.label) {
            matching_basin <- basin
            break
        }
    }

    ## Return the result with appropriate message
    if (is.null(matching_basin)) {
        message(paste0("Basin with label '", basin.label, "' not found."))
        return(NULL)
    } else {
        return(matching_basin)
    }
}

#' Union of Gradient Flow Basins
#'
#' Get the union of vertices from specified basins in a basin complex.
#'
#' @param object An object of class "basin_cx" returned by create.basin.cx().
#' @param x Backward-compatible alias for \code{object}.
#' @param ids Character vector of basin labels (e.g., c("M1", "m2", "M3")).
#' @param ... Additional arguments (ignored).
#'
#' @return An integer vector of vertex indices contained in the union of the specified basins.
#'
#' @export
basins.union <- function(object, ids, x = object, ...) {
    object <- x
    x <- object
    ## Check if x is a basin_cx object
    if (!inherits(x, "basin_cx")) {
        stop("'x' must be an object of class 'basin_cx'")
    }

    ## Check if ids is a character vector
    if (!is.character(ids)) {
        stop("'ids' must be a character vector of basin labels")
    }

    ## Initialize an empty set for the union
    union_vertices <- integer(0)

    ## Process each basin label
    for (id in ids) {
        ## Determine if this is a minimum (m) or maximum (M) basin
        is_max <- substr(id, 1, 1) == "M"

        ## Get the list of basins (either ascending for minima or descending for maxima)
        basins_list <- if (is_max) x$basins$descending else x$basins$ascending

        ## Find the basin with the matching label
        basin_idx <- which(sapply(basins_list, function(b) b$label == id))

        if (length(basin_idx) > 0) {
            ## Extract vertices from the basin matrix (first column contains vertex indices)
            basin_vertices <- basins_list[[basin_idx]]$basin[, 1]

            ## Add to the union
            union_vertices <- union(union_vertices, basin_vertices)
        } else {
            warning(paste("Basin with label", id, "not found"))
        }
    }

    ## Sort the vertex indices for consistency
    union_vertices <- sort(union_vertices)

    return(union_vertices)
}


#' Extract vertices belonging to a specific basin from a basin complex
#'
#' @description
#' This function extracts all vertex indices that belong to a specified basin from
#' a gradient flow basin complex object. The basin is identified by its label (e.g., "M1" or "m2").
#'
#' @param x An object of class 'basin_cx', typically returned by \code{create.basin.cx()}
#' @param id A character string specifying the basin label to extract (e.g., "M1" for first maximum, "m2" for second minimum)
#'
#' @return An integer vector containing the sorted vertex indices belonging to the specified basin
#'
#' @details
#' Basin labels follow the naming convention of the basin complex: "Mx" for maxima
#' basins and "mx" for minima basins, where x is a sequential number assigned based
#' on function value ordering.
#'
#' @seealso \code{\link{create.basin.cx}}
#' @export
get.basin.vertices <- function(x, id) {
    ## Check if x is a basin_cx object
    if (!inherits(x, "basin_cx")) {
        stop("'x' must be an object of class 'basin_cx'")
    }

    ## Check if ids is a character vector
    if (!is.character(id)) {
        stop("'id' must be a character basin label")
    }

    ## Initialize an empty set for the union
    vertices <- integer(0)

    ## Determine if this is a minimum (m) or maximum (M) basin
    is_max <- substr(id, 1, 1) == "M"

    ## Get the list of basins (either ascending for minima or descending for maxima)
    basins_list <- if (is_max) x$basins$descending else x$basins$ascending

    ## Find the basin with the matching label
    basin_idx <- which(sapply(basins_list, function(b) b$label == id))

    if (length(basin_idx) > 0) {
        ## Extract vertices from the basin matrix (first column contains vertex indices)
        vertices <- basins_list[[basin_idx]]$basin[, 1]
    } else {
        warning(paste("Basin with label", id, "not found"))
    }

    ## Sort the vertex indices for consistency
    vertices <- sort(vertices)

    return(vertices)
}


#' Construct Basin Complex Graph
#'
#' @description
#' Constructs a graph representation of a basin complex, where nodes represent basins
#' (local maxima and minima) and edges represent their relationships. Edges are weighted
#' based on basin overlap, function value differences, or other basin metrics.
#'
#' @param basin_cx A basin complex object of class "basin_cx" returned by create.basin.cx()
#' @param min_intersection Minimum number of vertices that must be shared between basins
#'        for an edge to be created (default: 1)
#' @param ms_edges_only Logical; if TRUE, only include edges between ascending and descending basins,
#'        ignoring ascending-ascending and descending-descending connections (default: FALSE)
#' @param weight_type Character string specifying the edge weight metric to use. Options are:
#'        "dice" (Dice-Sorensen similarity), "jaccard" (Jaccard index), "overlap" (overlap size),
#'        or "y_diff" (function value difference) (default: "dice")
#'
#' @details
#' This function creates a graph where nodes represent basins (either local minima or maxima)
#' and edges represent relationships between them. The resulting graph provides insight into
#' the topological structure of the function on the original graph.
#'
#' Edge weights can be based on different metrics:
#' \itemize{
#'   \item dice: Dice-Sorensen similarity index \eqn{(2|A \cap B|/(|A| + |B|))}
#'   \item jaccard: Jaccard similarity index \eqn{(|A \cap B|/|A \cup B|)}
#'   \item overlap: Raw number of vertices in the intersection
#'   \item y_diff: Absolute difference in function value between extrema
#' }
#'
#' When ms_edges_only = TRUE, only edges between an ascending basin (minimum) and a
#' descending basin (maximum) are included, producing a bipartite graph.
#'
#' @return A list with class "basin_cx_graph" containing:
#' \describe{
#'   \item{adjacency_list}{List of adjacency lists for the basin complex graph}
#'   \item{weights_list}{List of edge weights}
#'   \item{intersection_matrix}{Matrix where \code{[i,j]} is the number of vertices shared between basins i and j}
#'   \item{similarity_matrix}{Matrix of similarity values (based on weight_type) between basins}
#'   \item{y_diff_matrix}{Matrix of function value differences between extrema}
#'   \item{basin_metadata}{Data frame with information about each basin (type, extremum index, etc.)}
#'   \item{n_ascending}{Number of ascending basins}
#'   \item{n_descending}{Number of descending basins}
#'   \item{weight_type}{The type of weight metric used}
#' }
#'
#' @seealso \code{\link{create.basin.cx}}
#'
#' @export
construct.basin.cx.graph <- function(basin_cx,
                                     min_intersection = 1,
                                     ms_edges_only = FALSE,
                                     weight_type = c("dice", "jaccard", "overlap", "y_diff")) {

  ## Check input
  if (!inherits(basin_cx, "basin_cx")) {
    stop("Input must be an object of class 'basin_cx'")
  }

  ## Match weight_type argument
  weight_type <- match.arg(weight_type)

  ## Extract basins
  ascending_basins <- basin_cx$basins$ascending
  descending_basins <- basin_cx$basins$descending

  ## Extract vertices from each basin
  ascending_vertices <- lapply(ascending_basins, function(b) b$vertices)
  descending_vertices <- lapply(descending_basins, function(b) b$vertices)

  ## Create a list of all basin vertices
  all_basin_vertices <- c(ascending_vertices, descending_vertices)

  ## Get basin sizes
  basin_sizes <- sapply(all_basin_vertices, length)

  ## Get number of basins
  n_ascending <- length(ascending_basins)
  n_descending <- length(descending_basins)
  n_total <- n_ascending + n_descending

  if (n_total == 0) {
    stop("No basins found in the basin complex")
  }

  ## Create basin labels
  ascending_labels <- names(ascending_basins)
  descending_labels <- names(descending_basins)
  all_labels <- c(ascending_labels, descending_labels)

  ## Get the extrema indices and values
  min_indices <- sapply(ascending_basins, function(b) b$vertex)
  max_indices <- sapply(descending_basins, function(b) b$vertex)
  min_values <- sapply(ascending_basins, function(b) b$value)
  max_values <- sapply(descending_basins, function(b) b$value)

  ## Compute intersection matrix
  intersection_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  rownames(intersection_matrix) <- all_labels
  colnames(intersection_matrix) <- all_labels

  ## Calculate intersection sizes
  for (i in 1:n_total) {
    for (j in i:n_total) {
      if (i == j) {
        ## Self-intersection = basin size
        intersection_matrix[i, i] <- basin_sizes[i]
      } else {
        ## Calculate intersection size
        intersection_size <- length(intersect(all_basin_vertices[[i]], all_basin_vertices[[j]]))

        ## Store the result in both directions
        if (intersection_size >= min_intersection) {
          intersection_matrix[i, j] <- intersection_size
          intersection_matrix[j, i] <- intersection_size
        }
      }
    }
  }

  ## If ms_edges_only is TRUE, zero out non-ascending-descending connections
  if (ms_edges_only) {
    ## Zero out ascending-ascending connections (upper-left block)
    if (n_ascending > 1) {
      intersection_matrix[1:n_ascending, 1:n_ascending] <- 0
    }

    ## Zero out descending-descending connections (lower-right block)
    if (n_descending > 1) {
      descending_start <- n_ascending + 1
      descending_end <- n_total
      intersection_matrix[descending_start:descending_end, descending_start:descending_end] <- 0
    }
  }

  ## Calculate similarity matrix based on weight_type
  similarity_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  rownames(similarity_matrix) <- all_labels
  colnames(similarity_matrix) <- all_labels

  for (i in 1:n_total) {
    for (j in 1:n_total) {
      if (i == j) {
        ## Self-similarity = 1
        similarity_matrix[i, i] <- 1
      } else if (intersection_matrix[i, j] > 0) {
        ## Calculate similarity based on selected metric
        if (weight_type == "dice") {
          ## Dice-Sorensen index = 2|intersection| / (|A| + |B|)
          similarity <- 2 * intersection_matrix[i, j] / (basin_sizes[i] + basin_sizes[j])
        } else if (weight_type == "jaccard") {
          ## Jaccard index = |intersection| / |union|
          union_size <- basin_sizes[i] + basin_sizes[j] - intersection_matrix[i, j]
          similarity <- intersection_matrix[i, j] / union_size
        } else if (weight_type == "overlap") {
          ## Overlap index
          similarity <- intersection_matrix[i, j] / min(c(basin_sizes[i], basin_sizes[j]))
        } else {
          ## Default to Dice-Sorensen index
          similarity <- 2 * intersection_matrix[i, j] / (basin_sizes[i] + basin_sizes[j])
        }

        similarity_matrix[i, j] <- similarity
      }
    }
  }

  ## Calculate y-value differences matrix
  y_diff_matrix <- matrix(NA, nrow = n_total, ncol = n_total)
  rownames(y_diff_matrix) <- all_labels
  colnames(y_diff_matrix) <- all_labels

  for (i in 1:n_ascending) {
    for (j in (n_ascending+1):n_total) {
      if (intersection_matrix[i, j] > 0) {
        ## Get the function values
        min_value <- min_values[i]
        max_value <- max_values[j - n_ascending]

        ## Calculate difference
        y_diff <- max_value - min_value
        y_diff_matrix[i, j] <- y_diff
        y_diff_matrix[j, i] <- y_diff
      }
    }
  }

  ## If weight_type is "y_diff", use y_diff_matrix for the similarity matrix where defined
  if (weight_type == "y_diff") {
    for (i in 1:n_total) {
      for (j in 1:n_total) {
        if (!is.na(y_diff_matrix[i, j])) {
          ## Invert the difference so larger values = more similar
          max_diff <- max(y_diff_matrix, na.rm = TRUE)
          if (max_diff > 0) {
            similarity_matrix[i, j] <- 1 - (y_diff_matrix[i, j] / max_diff)
          }
        }
      }
    }
  }

  ## Convert similarity matrix to adjacency list and weights list
  adjacency_list <- vector("list", n_total)
  weights_list <- vector("list", n_total)
  names(adjacency_list) <- all_labels
  names(weights_list) <- all_labels

  for (i in 1:n_total) {
    adjacency_list[[i]] <- integer(0)
    weights_list[[i]] <- numeric(0)

    for (j in 1:n_total) {
      if (i != j && similarity_matrix[i, j] > 0) {
        ## Add to adjacency and weights lists
        adjacency_list[[i]] <- c(adjacency_list[[i]], j)
        weights_list[[i]] <- c(weights_list[[i]], similarity_matrix[i, j])
      }
    }
  }

  ## Create basin metadata
  basin_metadata <- data.frame(
    index = 1:n_total,
    label = all_labels,
    type = c(rep("ascending", n_ascending), rep("descending", n_descending)),
    original_index = c(1:n_ascending, 1:n_descending),
    size = basin_sizes,
    extremum_index = c(min_indices, max_indices),
    y_value = c(min_values, max_values),
    stringsAsFactors = FALSE
  )

  ## Create and return result
  result <- list(
    adjacency_list = adjacency_list,
    weights_list = weights_list,
    intersection_matrix = intersection_matrix,
    similarity_matrix = similarity_matrix,
    y_diff_matrix = y_diff_matrix,
    basin_metadata = basin_metadata,
    n_ascending = n_ascending,
    n_descending = n_descending,
    weight_type = weight_type
  )

  class(result) <- "basin_cx_graph"
  return(result)
}

#' Plot Basin Complex Graph
#'
#' @description
#' Creates a visualization of a basin complex graph showing connections between
#' basins of attraction for local extrema. The plot shows the basins as nodes
#' and their relationships as edges, with customizable styling options.
#'
#' @param basin_cx_graph An object of class "basin_cx_graph" returned by
#'        construct.basin.cx.graph(), or a basin complex object of class "basin_cx"
#' @param coords Matrix of coordinates for nodes, or NULL to compute them automatically
#' @param method Embedding method for node layout: "fr" for Fruchterman-Reingold
#'        or "kk" for Kamada-Kawai (default: "fr")
#' @param use_weights Whether to use edge weights in the layout calculation (default: TRUE)
#' @param invert_weights Whether to invert weights for embedding (default: TRUE)
#' @param edge_color Color for edges between ascending and descending basins (default: "gray50")
#' @param edge_width_factor Factor to scale edge width by similarity (default: 1)
#' @param edge_lmin_lmin_color Color for edges between ascending basins (default: "gray20")
#' @param edge_lmax_lmax_color Color for edges between descending basins (default: "gray20")
#' @param max_color Color for local maximum nodes (default: "red")
#' @param min_color Color for local minimum nodes (default: "blue")
#' @param label_max_color Color for maximum labels (default: "black")
#' @param label_min_color Color for minimum labels (default: "black")
#' @param node_size Base size for nodes (default: 3)
#' @param node_size_factor Factor to scale node size by basin size (default: 0.01)
#' @param arrow_length Length of arrows (default: 0.15)
#' @param arrow_gap Gap between arrow and edge (default: 0.05)
#' @param text_size Base text size for labels (default: 0.8)
#' @param show_y_values Whether to show y-values at nodes (default: TRUE)
#' @param show_basin_size Whether to show basin sizes (default: TRUE)
#' @param show_cell_size Whether to show intersection sizes on edges (default: TRUE)
#' @param add Whether to add to existing plot (default: FALSE)
#' @param add_box Whether to add a box around the plot (default: FALSE)
#' @param x_lim X-axis limits, or NULL for automatic
#' @param y_lim Y-axis limits, or NULL for automatic
#' @param x_lim_pad_factor Padding factor for x limits (default: 0.15)
#' @param y_lim_pad_factor Padding factor for y limits (default: 0.15)
#' @param main Plot title (default: "")
#' @param xlab X-axis label (default: "")
#' @param ylab Y-axis label (default: "")
#' @param label_adj_max List of adj values for maximum labels, with names matching max labels
#' @param label_adj_min List of adj values for minimum labels, with names matching min labels
#' @param info_adj_max List of adj values for maximum info text, with names matching max labels
#' @param info_adj_min List of adj values for minimum info text, with names matching min labels
#' @param edge_label_adj Named list of adj values for edge labels
#' @param ... Additional parameters passed to plot
#'
#' @details
#' This function creates a visual representation of the basin complex graph, which shows
#' the relationships between basins of local extrema. Different colors and shapes are
#' used to distinguish between ascending (minimum) and descending (maximum) basins.
#'
#' Local maxima are shown as upward-pointing triangles, and local minima as downward-pointing
#' triangles. Edge styling differs based on the basin types being connected:
#' \itemize{
#'   \item Ascending-descending edges: Solid lines with a single arrow in the middle
#'   \item Ascending-ascending edges: Dashed lines with arrows pointing inward
#'   \item Descending-descending edges: Dashed lines with arrows pointing outward
#' }
#'
#' Node positions are calculated using force-directed graph layout algorithms:
#' \itemize{
#'   \item "fr": Fruchterman-Reingold algorithm, which tends to distribute nodes evenly
#'   \item "kk": Kamada-Kawai algorithm, which attempts to preserve graph-theoretic distances
#' }
#'
#' The size of nodes reflects the size of the corresponding basin, and edge widths
#' represent the strength of the relationship (based on similarity or intersection size).
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{coords}{The coordinates of all nodes}
#'   \item{adjustments}{All adjustment parameters used for labels}
#'   \item{edge_label_info}{Detailed information about edge labels including positions}
#'   \item{basin_meta}{Basin metadata}
#'   \item{n_ascending}{Number of ascending basins}
#'   \item{n_descending}{Number of descending basins}
#' }
#'
#' @seealso \code{\link{construct.basin.cx.graph}}, \code{\link{create.basin.cx}}
#'
#' @keywords internal
#' @noRd
plot_basin_cx_graph <- function(basin_cx_graph,
                              coords = NULL,
                              method = c("fr", "kk"),
                              use_weights = TRUE,
                              invert_weights = TRUE,
                              edge_color = "gray50",
                              edge_width_factor = 1,
                              edge_lmin_lmin_color = "gray20",
                              edge_lmax_lmax_color = "gray20",
                              max_color = "red",
                              min_color = "blue",
                              label_max_color = "black",
                              label_min_color = "black",
                              node_size = 3,
                              node_size_factor = 0.01,
                              arrow_length = 0.15,
                              arrow_gap = 0.05,
                              text_size = 0.8,
                              show_y_values = TRUE,
                              show_basin_size = TRUE,
                              show_cell_size = TRUE,
                              add = FALSE,
                              add_box = FALSE,
                              x_lim = NULL,
                              y_lim = NULL,
                              x_lim_pad_factor = 0.15,
                              y_lim_pad_factor = 0.15,
                              main = "",
                              xlab = "",
                              ylab = "",
                              label_adj_max = NULL,
                              label_adj_min = NULL,
                              info_adj_max = NULL,
                              info_adj_min = NULL,
                              edge_label_adj = NULL,
                              ...) {

  ## Match method argument
  method <- match.arg(method)

  ## Check if input is a basin_cx_graph or a basin_cx and convert if needed
  if (inherits(basin_cx_graph, "basin_cx")) {
    basin_cx_graph <- construct.basin.cx.graph(basin_cx_graph)
  } else if (!inherits(basin_cx_graph, "basin_cx_graph")) {
    stop("Input must be of class 'basin_cx_graph' or 'basin_cx'")
  }

  ## Extract graph properties
  n_ascending <- basin_cx_graph$n_ascending
  n_descending <- basin_cx_graph$n_descending
  n_total <- n_ascending + n_descending

  ## Check if the graph is empty
  if (n_total == 0) {
    warning("No basins found in the basin complex graph")
    plot(0, 0, type = "n", main = main, xlab = xlab, ylab = ylab,
         xlim = c(-1, 1), ylim = c(-1, 1), ...)
    text(0, 0, "Empty basin complex", cex = 1.2)
    return(invisible(NULL))
  }

  ## Prepare weights list for embedding
  weights_list <- if (use_weights) basin_cx_graph$weights_list else NULL

  ## Create 2D embedding of the graph if not provided
  if (is.null(coords)) {
    ## Load the graph.embedding function or implement it
    ## This is a placeholder - we'd need to implement or import a graph embedding function
    coords <- graph.embedding(
      adj.list = basin_cx_graph$adjacency_list,
      weights.list = weights_list,
      invert.weights = invert_weights,
      dim = 2,
      method = method,
      verbose = FALSE
    )
  }

  ## Extract basin metadata
  basin_meta <- basin_cx_graph$basin_metadata

  ## Get basin labels
  min_labels <- basin_meta$label[1:n_ascending]
  max_labels <- basin_meta$label[(n_ascending+1):n_total]

  ## Initialize storage for all adj values used in the plot
  all_adj_values <- list(
    label_adj_max = list(),
    label_adj_min = list(),
    info_adj_max = list(),
    info_adj_min = list(),
    edge_label_adj = list()
  )

  ## Set up default adjustment values for labels if not provided
  if (is.null(label_adj_max)) {
    label_adj_max <- lapply(1:n_descending, function(i) c(0.3, -0.8))
    names(label_adj_max) <- max_labels
  }

  if (is.null(label_adj_min)) {
    label_adj_min <- lapply(1:n_ascending, function(i) c(0.3, 1.5))
    names(label_adj_min) <- min_labels
  }

  if (is.null(info_adj_max)) {
    info_adj_max <- lapply(1:n_descending, function(i) c(1.2, 0.5))
    names(info_adj_max) <- max_labels
  }

  if (is.null(info_adj_min)) {
    info_adj_min <- lapply(1:n_ascending, function(i) c(1.2, 0.5))
    names(info_adj_min) <- min_labels
  }

  ## Store the adj values
  all_adj_values$label_adj_max <- label_adj_max
  all_adj_values$label_adj_min <- label_adj_min
  all_adj_values$info_adj_max <- info_adj_max
  all_adj_values$info_adj_min <- info_adj_min

  ## Set up plot if not adding to existing one
  if (!add) {
    ## Calculate plot bounds with some padding
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    x_pad <- diff(x_range) * x_lim_pad_factor
    y_pad <- diff(y_range) * y_lim_pad_factor

    if (is.null(x_lim)) {
      x_lim <- c(min(x_range) - x_pad, max(x_range) + x_pad)
    }

    if (is.null(y_lim)) {
      y_lim <- c(min(y_range) - y_pad, max(y_range) + y_pad)
    }

    plot(0, 0, type = "n",
         xlim = x_lim,
         ylim = y_lim,
         xlab = xlab, ylab = ylab, main = main,
         axes = FALSE,
         ...)

    if (add_box) {
      box()
    }
  }

  ## Helper function for drawing an arrow in the middle of an edge
  draw_mid_arrow <- function(x0, y0, x1, y1, eps = arrow_gap, length = arrow_length, col = edge_color) {
    ## Calculate midpoint
    mid_x <- (x0 + x1) / 2
    mid_y <- (y0 + y1) / 2

    ## Calculate direction vector
    dir_x <- x1 - x0
    dir_y <- y1 - y0
    dir_len <- sqrt(dir_x^2 + dir_y^2)

    ## Normalize direction vector
    unit_dir_x <- dir_x / dir_len
    unit_dir_y <- dir_y / dir_len

    ## Calculate points before and after midpoint
    before_x <- mid_x - eps * unit_dir_x
    before_y <- mid_y - eps * unit_dir_y
    after_x <- mid_x + eps * unit_dir_x
    after_y <- mid_y + eps * unit_dir_y

    ## Draw arrow
    graphics::arrows(before_x, before_y, after_x, after_y, length = length, col = col)
  }

  ## Helper function for drawing an arrow at a specific position along an edge
  draw_p_mid_arrow <- function(x0, y0, x1, y1, p = 0.5, eps = arrow_gap,
                             length = arrow_length, col = edge_color) {
    ## Calculate point at position p along the line
    p_x <- x0 + p * (x1 - x0)
    p_y <- y0 + p * (y1 - y0)

    ## Calculate direction vector
    dir_x <- x1 - x0
    dir_y <- y1 - y0
    dir_len <- sqrt(dir_x^2 + dir_y^2)

    ## Normalize direction vector
    unit_dir_x <- dir_x / dir_len
    unit_dir_y <- dir_y / dir_len

    ## Calculate points before and after midpoint
    before_x <- p_x - eps * unit_dir_x
    before_y <- p_y - eps * unit_dir_y
    after_x <- p_x + eps * unit_dir_x
    after_y <- p_y + eps * unit_dir_y

    ## Draw arrow
    graphics::arrows(before_x, before_y, after_x, after_y, length = length, col = col)
  }

  ## Initialize storage for edge label adjustments
  edge_label_info <- list()

  ## Draw edges based on basin types
  adj_list <- basin_cx_graph$adjacency_list
  intersection_sizes <- basin_cx_graph$intersection_matrix
  y_diffs <- basin_cx_graph$y_diff_matrix

  ## Process each vertex to draw edges
  for (i in 1:n_total) {
    if (length(adj_list[[i]]) == 0) next

    for (j in adj_list[[i]]) {
      ## Only process each edge once (when i < j)
      if (i < j) {
        x0 <- coords[i, 1]
        y0 <- coords[i, 2]
        x1 <- coords[j, 1]
        y1 <- coords[j, 2]

        ## Get intersection size and basin sizes for scaling
        intersection_size <- intersection_sizes[i, j]
        basin_i_size <- basin_meta$size[i]
        basin_j_size <- basin_meta$size[j]

        ## Calculate edge thickness based on intersection size
        rel_size <- intersection_size / max(basin_i_size, basin_j_size)
        edge_width <- 1 + edge_width_factor * rel_size

        ## Get actual basin labels for the edge
        from_label <- basin_meta$label[i]
        to_label <- basin_meta$label[j]

        ## Different edge styling based on basin types
        if (i <= n_ascending && j > n_ascending) {
          ## lmin-lmax connection (ascending to descending)

          ## Create edge key using actual basin labels
          edge_key <- paste0(from_label, "-", to_label)

          ## Draw the line first
          segments(x0, y0, x1, y1, col = edge_color, lwd = edge_width)

          ## Draw arrow in the middle
          draw_mid_arrow(x0, y0, x1, y1, col = edge_color)

          ## Calculate y-value difference
          y_diff <- y_diffs[i, j]

          ## Add text for y-diff and cell size if requested
          if (!is.na(y_diff) && (show_cell_size || show_y_values)) {
            ## Position text near the middle of the edge
            text_x <- (x0 + x1) / 2
            text_y <- (y0 + y1) / 2

            ## Calculate offset perpendicular to the edge
            edge_dir_x <- x1 - x0
            edge_dir_y <- y1 - y0
            edge_len <- sqrt(edge_dir_x^2 + edge_dir_y^2)
            perp_x <- -edge_dir_y / edge_len
            perp_y <- edge_dir_x / edge_len

            ## Position text near the middle of the edge with offset
            text_x <- (x0 + x1) / 2 + perp_x * edge_len * 0.07
            text_y <- (y0 + y1) / 2 + perp_y * edge_len * 0.07

            ## Get adjustment for this specific edge if available
            edge_adj <- NULL
            if (!is.null(edge_label_adj) && edge_key %in% names(edge_label_adj)) {
              edge_adj <- edge_label_adj[[edge_key]]
            } else {
              ## Default: offset perpendicular to edge
              edge_adj <- c(0.5, 0.5)
              offset_factor <- 0.07
              text_x <- text_x + perp_x * edge_len * offset_factor
              text_y <- text_y + perp_y * edge_len * offset_factor
            }

            ## Store edge adjustment for return
            all_adj_values$edge_label_adj[[edge_key]] <- edge_adj
            edge_label_info[[edge_key]] <- list(
              x = text_x,
              y = text_y,
              adj = edge_adj,
              from_label = from_label,
              to_label = to_label,
              y_diff = y_diff,
              intersection_size = intersection_size
            )

            ## Use separate text elements for y-diff and intersection size
            if (show_y_values) {
              ## Create LaTeX expression for delta y
              delta_text <- bquote(Delta*y == .(format(round(y_diff, 2), nsmall = 2)))
              text(text_x, text_y - 0.03, delta_text, cex = text_size, adj = edge_adj)
            }

            if (show_cell_size) {
              ## Create LaTeX expression for intersection size
              intersect_text <- bquote("|B" %intersection% "B'|" == .(intersection_size))
              text(text_x, text_y + (if(show_y_values) 0.03 else 0),
                   intersect_text, cex = text_size, adj = edge_adj)
            }
          }

        } else if (i <= n_ascending && j <= n_ascending) {
          ## lmin-lmin connection (between ascending basins)

          ## Create edge key using actual basin labels
          edge_key <- paste0(from_label, "-", to_label)

          ## Draw dashed line
          segments(x0, y0, x1, y1, col = edge_lmin_lmin_color, lty = 2, lwd = edge_width)

          ## Draw arrows pointing inward (25% from each end)
          draw_p_mid_arrow(x0, y0, x1, y1, p = 0.25, col = edge_lmin_lmin_color,
                         length = arrow_length, eps = arrow_gap)
          draw_p_mid_arrow(x1, y1, x0, y0, p = 0.25, col = edge_lmin_lmin_color,
                         length = arrow_length, eps = arrow_gap)

        } else if (i > n_ascending && j > n_ascending) {
          ## lmax-lmax connection (between descending basins)

          ## Create edge key using actual basin labels
          edge_key <- paste0(from_label, "-", to_label)

          ## Draw dashed line
          segments(x0, y0, x1, y1, col = edge_lmax_lmax_color, lty = 2, lwd = edge_width)

          ## Draw arrows pointing outward (25% from each end)
          draw_p_mid_arrow(x1, y1, x0, y0, p = 0.75, col = edge_lmax_lmax_color,
                         length = arrow_length, eps = arrow_gap)
          draw_p_mid_arrow(x0, y0, x1, y1, p = 0.75, col = edge_lmax_lmax_color,
                         length = arrow_length, eps = arrow_gap)
        }
      }
    }
  }

  ## Draw nodes and apply customized label positioning
  for (i in 1:n_total) {
    x <- coords[i, 1]
    y <- coords[i, 2]

    ## Get basin info
    is_max <- i > n_ascending
    basin_size <- basin_meta$size[i]
    extremum_index <- basin_meta$extremum_index[i]
    y_value <- basin_meta$y_value[i]
    basin_label <- basin_meta$label[i]

    ## Scale node size by basin size
    point_size <- node_size + node_size_factor * sqrt(basin_size)

    ## Draw node and apply appropriate label adjustments
    if (is_max) {
      ## Up-facing triangle for maxima
      points(x, y, pch = 24, col = max_color, bg = max_color, cex = point_size)

      ## Use actual basin label as key
      label_key <- basin_label

      ## Get adjustment for this specific node label
      node_adj <- if (!is.null(label_adj_max) && label_key %in% names(label_adj_max)) {
        label_adj_max[[label_key]]
      } else {
        c(0.3, -0.8)  ## Default if not found
      }

      ## Add label
      text(x, y, basin_label, adj = node_adj, cex = text_size * 1.2, col = label_max_color)

      ## Add information text with custom adjustment
      if (show_y_values || show_basin_size) {
        info_text <- ""
        if (show_y_values) {
          info_text <- paste0(info_text, sprintf("y=%.2f", y_value))
        }
        if (show_basin_size) {
          if (info_text != "") info_text <- paste0(info_text, "\n")
          info_text <- paste0(info_text, sprintf("|B|=%d", basin_size))
        }

        ## Get adjustment for this specific info label
        info_adj <- if (!is.null(info_adj_max) && label_key %in% names(info_adj_max)) {
          info_adj_max[[label_key]]
        } else {
          c(1.2, 0.5)  ## Default if not found
        }

        text(x, y, info_text, adj = info_adj, cex = text_size)
      }
    } else {
      ## Down-facing triangle for minima
      points(x, y, pch = 25, col = min_color, bg = min_color, cex = point_size)

      ## Use actual basin label as key
      label_key <- basin_label

      ## Get adjustment for this specific node label
      node_adj <- if (!is.null(label_adj_min) && label_key %in% names(label_adj_min)) {
        label_adj_min[[label_key]]
      } else {
        c(0.3, 1.5)  ## Default if not found
      }

      ## Add label
      text(x, y, basin_label, adj = node_adj, cex = text_size * 1.2, col = label_min_color)

      ## Add information text with custom adjustment
      if (show_y_values || show_basin_size) {
        info_text <- ""
        if (show_y_values) {
          info_text <- paste0(info_text, sprintf("y=%.2f", y_value))
        }
        if (show_basin_size) {
          if (info_text != "") info_text <- paste0(info_text, "\n")
          info_text <- paste0(info_text, sprintf("|B|=%d", basin_size))
        }

        ## Get adjustment for this specific info label
        info_adj <- if (!is.null(info_adj_min) && label_key %in% names(info_adj_min)) {
          info_adj_min[[label_key]]
        } else {
          c(1.2, 0.5)  ## Default if not found
        }

        text(x, y, info_text, adj = info_adj, cex = text_size)
      }
    }
  }

  ## Return comprehensive information about the plot
  return_value <- list(
    coords = coords,
    adjustments = all_adj_values,
    edge_label_info = edge_label_info,
    basin_meta = basin_meta,
    n_ascending = n_ascending,
    n_descending = n_descending
  )

  ## Return invisibly
  invisible(return_value)
}

#' Plot method for basin_cx objects
#'
#' @param x A basin_cx object
#' @param y Not used
#' @param type Type of plot: "graph" (default), "basins", "cells", or "comparison"
#' @param ... Additional arguments passed to specific plotting functions
#'
#'@export
plot.basin_cx <- function(x, y, ..., type = c("graph", "basins", "cells", "comparison")) {
  type <- match.arg(type)

  switch(type,
         comparison = {
           ## Plot original vs simplified function values
           if (missing(y)) {
             y <- seq_along(x$original_y)
           }

           plot(y, x$original_y, type = "l", col = "gray",
                main = "Original vs Simplified Function",
                xlab = "Index", ylab = "Value", ...)
           lines(y, x$harmonic_predictions, col = "red")
           legend("topright", legend = c("Original", "Simplified"),
                  col = c("gray", "red"), lty = 1)
         },
         basins = {
           ## Plot individual basins
           ## Implementation depends on the specific requirements
           warning("Basins plot type not yet implemented")
         },
         cells = {
           ## Plot cell complex
           ## Implementation depends on the specific requirements
           warning("Cells plot type not yet implemented")
         },
         graph = {
           ## Plot basin complex graph
           plot_basin_cx_graph(x, ...)
         })

  invisible(x)
}

#' Merge Two Basins in a Basin Complex
#'
#' @description
#' Merges one basin into another within a \code{basin_cx} object by updating
#' vertices, labels, and the corresponding cell complex.
#'
#' @param object An object of class \code{basin_cx}, as returned by \code{create.basin.cx()}.
#' @param basin_cx Backward-compatible alias for \code{object}.
#' @param absorbing.label Label of the basin that will absorb another (e.g., "m1", "M2").
#' @param absorbed.label Label of the basin to be absorbed (must be same type as \code{absorbing.label}).
#' @param ... Additional arguments (ignored).
#'
#' @return A new \code{basin_cx} object with updated basins and cell complex.
#' @export
basin.cx.merge <- function(object, absorbing.label, absorbed.label, basin_cx = object, ...) {
    object <- basin_cx
    basin_cx <- object
    stopifnot(inherits(basin_cx, "basin_cx"))

    # Declare variables to avoid R CMD check warnings
    label <- NULL
    is_max <- NULL

    ## Determine basin type
    if (substr(absorbing.label, 1, 1) == "m" && substr(absorbed.label, 1, 1) == "m") {
        basin_type <- "ascending"
    } else if (substr(absorbing.label, 1, 1) == "M" && substr(absorbed.label, 1, 1) == "M") {
        basin_type <- "descending"
    } else {
        stop("Both basin labels must start with 'm' or both with 'M'")
    }
    basins <- basin_cx$basins[[basin_type]]
    absorb_idx <- which(names(basins) == absorbing.label)
    absorbed_idx <- which(names(basins) == absorbed.label)
    if (length(absorb_idx) == 0 || length(absorbed_idx) == 0) {
        stop("One or both basin labels not found")
    }
    ## Merge vertices
    merged_vertices <- sort(unique(c(
        basins[[absorb_idx]]$vertices,
        basins[[absorbed_idx]]$vertices
    )))
    basins[[absorb_idx]]$vertices <- merged_vertices
    ## Optional: record merge history
    if (is.null(basins[[absorb_idx]]$merge_history))
        basins[[absorb_idx]]$merge_history <- character(0)
    basins[[absorb_idx]]$merge_history <- c(
        basins[[absorb_idx]]$merge_history,
        paste("Absorbed", absorbed.label)
    )
    ## Remove absorbed basin
    basins[[absorbed_idx]] <- NULL
    basin_cx$basins[[basin_type]] <- basins
    ## Remove absorbed from local_extrema
    absorbed_vertex <- basin_cx$local_extrema$vertex_index[
                                                  basin_cx$local_extrema$label == absorbed.label
                                              ]
    basin_cx$local_extrema <- subset(
        basin_cx$local_extrema,
        basin_cx$local_extrema$label != absorbed.label
    )
    ## Update basins_df and *_basins_df
    basin_cx$basins_df <- do.call(rbind, lapply(unlist(basin_cx$basins, recursive = FALSE), function(b) {
        data.frame(
            evertex = b$vertex,
            value = b$value,
            is_max = as.integer(substr(b$label, 1, 1) == "M"),
            min_rspan = NA,
            max_rspan = NA,
            Drspan = NA,
            size = length(b$vertices),
            rsize = length(b$vertices) / length(basin_cx$original_y),
            label = b$label
        )
    }))
    basin_cx$lmin_basins_df <- subset(basin_cx$basins_df, basin_cx$basins_df$is_max == 0)
    rownames(basin_cx$lmin_basins_df) <- basin_cx$lmin_basins_df$label
    basin_cx$lmin_basins_df$label <- NULL
    basin_cx$lmax_basins_df <- subset(basin_cx$basins_df, basin_cx$basins_df$is_max == 1)
    rownames(basin_cx$lmax_basins_df) <- basin_cx$lmax_basins_df$label
    basin_cx$lmax_basins_df$label <- NULL
    ## Recompute cell complex
    basin_cx$cells <- create_basin_cell_complex(
        basin_cx,
        min.asc.desc.cell.size.thld = 1,
        min.asc.asc.cell.size.thld = 1,
        min.desc.desc.cell.size.thld = 1
    )
    ## Recompute basin graph
    basin_cx$graph <- construct.basin.cx.graph(basin_cx)
    return(basin_cx)
}

#' Compute Ascending-Descending Cell Intersection Matrix
#'
#' @param basin_cx An object of class 'basin_cx'
#' @return A square matrix where (i,j) is the number of shared vertices between cell i and cell j,
#'         with labels as lmin_label:lmax_label. Final row and column contain cell sizes.
compute_asc_desc_cell_intersection_matrix <- function(basin_cx) {
  if (!inherits(basin_cx, "basin_cx")) {
    stop("Input must be of class 'basin_cx'")
  }

  cells <- basin_cx$cells$asc_desc_cells
  if (is.null(cells) || length(cells) == 0) {
    return(matrix(0, 0, 0))
  }

  labels <- vapply(cells, function(cell) paste0(cell$lmin_label, ":", cell$lmax_label), character(1))
  n <- length(labels)

  mat <- matrix(0, n + 1, n + 1)
  rownames(mat) <- colnames(mat) <- c(labels, "cell_size")

  vertices_list <- lapply(cells, function(cell) cell$vertices)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- length(intersect(vertices_list[[i]], vertices_list[[j]]))
    }
    mat[n + 1, i] <- length(vertices_list[[i]]) # last row
    mat[i, n + 1] <- length(vertices_list[[i]]) # last column
  }

  mat[n + 1, n + 1] <- NA  # corner cell
  return(mat)
}


extract.basins.vertices <- function(basin_cx, basin_ids) {
    result_list <- list()

    for (id in basin_ids) {
        ## Determine if minimum or maximum basin from prefix
        is_max <- substr(id, 1, 1) == "M"
        basin_list <- if (is_max) basin_cx$basins$descending else basin_cx$basins$ascending

        ## Find basin with matching label
        basin_idx <- NULL
        for (i in seq_along(basin_list)) {
            if (!is.null(basin_list[[i]]$label) && basin_list[[i]]$label == id) {
                basin_idx <- i
                break
            }
        }

        if (!is.null(basin_idx)) {
            basin <- basin_list[[basin_idx]]
            ## Extract vertices (first column of basin matrix)
            result_list[[id]] <- basin$vertices
        } else {
            result_list[[id]] <- integer(0)
        }
    }

    return(result_list)
}

#' Create Taxonomic Labels for Vertices
#'
#' @description
#' Creates unique shortcut labels for vertices based on their taxonomic profiles.
#' The function generates two-letter shortcuts from species names and combines them
#' to create concise labels that represent the most abundant taxa in each vertex.
#'
#' @param vertices A numeric vector of vertex indices to be labeled.
#' @param state.space A matrix or data frame containing the state space data,
#'        where rows correspond to different states/samples.
#' @param taxonomy A data structure containing taxonomic information used by the
#'        prof.fn function to generate taxonomic profiles.
#' @param min.relAb.thld Numeric value (default: 0.05). Minimum relative abundance
#'        threshold. Taxa with abundance below this threshold will not be included
#'        in the label.
#'
#' @return A list with two components:
#'   \item{labels}{A named character vector containing the generated labels, with
#'         names corresponding to vertex indices.}
#'   \item{profiles}{A list of data frames containing the taxonomic profiles for
#'         each vertex.}
#'
#' @details The function first creates taxonomic profiles for each vertex using
#'     the prof.fn function. It then generates labels based on the most abundant
#'     taxa above the specified threshold. If multiple vertices get the same
#'     label, the function ensures uniqueness by adding additional taxa
#'     information. The labels consist of two-letter shortcuts for each taxon
#'     name, created by taking the first letter of the first two parts of the
#'     underscore-separated name.
#'
#' @export
create.taxonomic.labels <- function(vertices,
                                    state.space,
                                    taxonomy = NULL,
                                    min.relAb.thld = 0.05
                                    ) {

    # Helper function to create two-letter shortcuts from species names
    create.shortcut <- function(sp.name) {
        parts <- strsplit(sp.name, "_")[[1]]
        if (length(parts) >= 2) {
            paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))
        } else {
            paste0(substr(parts[1], 1, 1), "")
        }
    }

    labels <- character(length(vertices))
    names(labels) <- as.character(vertices)
    profiles <- list()

    # Process vertices
    for (i in seq_along(vertices)) {
        vertex <- vertices[i]
        id <- rownames(state.space)[vertex]

        # Get profile for this vertex
        profile <- prof.fn(id, state.space, taxonomy)
        # Store profile
        profiles[[as.character(vertex)]] <- profile

        abundances <- as.numeric(profile[, 1])
        above.threshold <- which(abundances >= min.relAb.thld)

        if (length(above.threshold) > 0) {
            taxa.names <- rownames(profile)[above.threshold][1]
            shortcuts <- sapply(taxa.names, create.shortcut)
            labels[i] <- paste(shortcuts, collapse = "")
        } else {
            labels[i] <- create.shortcut(rownames(profile)[1])
        }
    }

    # Ensure uniqueness of labels
    while (any(duplicated(labels))) {
        dups <- labels[duplicated(labels)]
        for (dup.label in unique(dups)) {
            dup.indices <- which(labels == dup.label)
            for (idx in dup.indices) {
                vertex <- vertices[idx]
                profile <- profiles[[as.character(vertex)]]

                current.taxa <- ceiling(nchar(labels[idx])/2)
                if (current.taxa < nrow(profile)) {
                    new.taxon <- create.shortcut(rownames(profile)[current.taxa + 1][1])
                    labels[idx] <- paste0(labels[idx], new.taxon)
                }
            }
        }
    }

    return(list(
        labels = labels,
        profiles = profiles
    ))
}


#' Replace a Basin Label in a Basin Complex Object
#'
#' @description
#' Replaces a specified label for a local extremum in a basin_cx object.
#' Updates all references to this label in extrema, basins, and cells.
#'
#' @param bcx An object of class "basin_cx"
#' @param old.label The existing label to be replaced
#' @param new.label The new label to replace with
#'
#' @return A modified object of class "basin_cx" with updated labels
#'
#' @export
replace.basin.label.basin_cx <- function(bcx, old.label, new.label) {
  if (!inherits(bcx, "basin_cx")) {
    stop("Input must be a 'basin_cx' object.")
  }

  # Update in local_extrema
  idx <- which(bcx$local_extrema$label == old.label)
  if (length(idx) == 0) {
    stop(sprintf("Label '%s' not found in local_extrema.", old.label))
  }
  vertex <- bcx$local_extrema$vertex[idx]
  is_max <- bcx$local_extrema$is_maximum[idx]
  bcx$local_extrema$label[idx] <- new.label

  # Update in lmin_basins or lmax_basins
  if (!is.null(bcx$lmin_basins) && !is.null(rownames(bcx$lmin_basins)) && !is_max) {
    if (old.label %in% rownames(bcx$lmin_basins)) {
      rownames(bcx$lmin_basins)[rownames(bcx$lmin_basins) == old.label] <- new.label
    }
  }
  if (!is.null(bcx$lmax_basins) && !is.null(rownames(bcx$lmax_basins)) && is_max) {
    if (old.label %in% rownames(bcx$lmax_basins)) {
      rownames(bcx$lmax_basins)[rownames(bcx$lmax_basins) == old.label] <- new.label
    }
  }

  # Update in cells
  for (i in seq_along(bcx$cells$asc_desc_cells)) {
    if (!is_max && bcx$cells$asc_desc_cells[[i]]$lmin_label == old.label) {
      bcx$cells$asc_desc_cells[[i]]$min_label <- new.label
    }
    if (is_max && bcx$cells$asc_desc_cells[[i]]$lmax_label == old.label) {
      bcx$cells$asc_desc_cells[[i]]$max_label <- new.label
    }

    # Update full label string
    bcx$cells$asc_desc_cells[[i]]$label <- paste0(bcx$cells$asc_desc_cells[[i]]$lmin_label, "-", bcx$cells$asc_desc_cells[[i]]$lmax_label)
  }

  return(bcx)
}

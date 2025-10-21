#' Relabel extrema with sequential labels after filtering
#'
#' @param extrema.df A data frame of extrema (output from extrema.summary or
#'   a filtered subset) containing columns: label, vertex, value, type, etc.
#' @param sort.ascending Logical. If TRUE, sort by value ascending (for minima).
#'   If FALSE, sort by value descending (for maxima). Default is NULL, which
#'   automatically uses FALSE for maxima and TRUE for minima based on the 'type'
#'   column.
#'
#' @return A data frame with the same structure as input but with relabeled
#'   extrema. Maxima are labeled M1, M2, ... in decreasing order of value.
#'   Minima are labeled m1, m2, ... in increasing order of value.
#'
#' @details
#' After filtering extrema based on criteria such as hop.idx, extremality
#' radius, or other properties, the labels (M1, M2, ..., m1, m2, ...) become
#' non-sequential. This function reassigns labels to maintain sequential
#' numbering in the appropriate order.
#'
#' If the input contains both maxima and minima, they are relabeled separately
#' within their respective types.
#'
#' @examples
#' # Filter maxima by hop index
#' filtered.max <- extr.summary[extr.summary$type == "max" &
#'                               extr.summary$hop.idx > 2, ]
#' # Relabel
#' relabeled.max <- relabel.extrema(filtered.max)
#'
#' # Works with minima too
#' filtered.min <- extr.summary[extr.summary$type == "min" &
#'                               extr.summary$hop.idx > 2, ]
#' relabeled.min <- relabel.extrema(filtered.min)
#'
#' # Or relabel both at once
#' filtered.both <- extr.summary[extr.summary$hop.idx > 2, ]
#' relabeled.both <- relabel.extrema(filtered.both)
#'
#' @export
relabel.extrema <- function(extrema.df, sort.ascending = NULL) {

    # Input validation
    if (!is.data.frame(extrema.df)) {
        stop("extrema.df must be a data frame")
    }

    required.cols <- c("label", "vertex", "value", "type")
    missing.cols <- setdiff(required.cols, names(extrema.df))
    if (length(missing.cols) > 0) {
        stop(paste("Missing required columns:", paste(missing.cols, collapse = ", ")))
    }

    if (nrow(extrema.df) == 0) {
        warning("Input data frame is empty; returning as-is")
        return(extrema.df)
    }

    # Check if we have mixed types
    types.present <- unique(extrema.df$type)

    if (length(types.present) > 1) {
        # Mixed types: process separately and combine
        result <- data.frame()

        if ("max" %in% types.present) {
            max.df <- extrema.df[extrema.df$type == "max", ]
            max.relabeled <- relabel.extrema(max.df, sort.ascending = FALSE)
            result <- rbind(result, max.relabeled)
        }

        if ("min" %in% types.present) {
            min.df <- extrema.df[extrema.df$type == "min", ]
            min.relabeled <- relabel.extrema(min.df, sort.ascending = TRUE)
            result <- rbind(result, min.relabeled)
        }

        rownames(result) <- NULL
        return(result)
    }

    # Single type: proceed with relabeling
    type <- types.present[1]

    # Determine sort order
    if (is.null(sort.ascending)) {
        sort.ascending <- (type == "min")
    }

    # Sort by value
    if (sort.ascending) {
        sorted.df <- extrema.df[order(extrema.df$value), ]
        label.prefix <- "m"
    } else {
        sorted.df <- extrema.df[order(-extrema.df$value), ]
        label.prefix <- "M"
    }

    # Assign new sequential labels
    sorted.df$label <- paste0(label.prefix, seq_len(nrow(sorted.df)))

    # Reset row names
    rownames(sorted.df) <- NULL

    return(sorted.df)
}

#' Computes Local Extrema Hop Neighborhoods in a Graph Function
#'
#' @description
#' Identifies local extrema (minima and maxima) of a function defined on a graph and computes
#' their extremum hop neighborhoods. For each extremum, this function determines the maximum hop distance
#' at which it maintains its extremum property and identifies boundary vertices where this
#' property is first violated.
#'
#' @details
#' A local extremum is a vertex whose function value is more extreme (higher for maxima, lower for minima)
#' than all its adjacent neighbors. The hop neighborhood analysis extends this concept by examining
#' how far the extremum property persists across the graph structure.
#'
#' For each extremum, the function computes:
#' \itemize{
#'   \item Its maximum hop radius (hop_idx) where it remains a strict extremum
#'   \item All vertices within this hop radius and their distances
#'   \item The boundary vertices at hop_idx + 1 where the extremum property is first violated
#' }
#'
#' This information is useful for identification of spurious local extrema.
#'
#' @param adj.list A list where each element contains the indices of vertices adjacent to vertex i.
#'                The list should be 1-indexed (as is standard in R).
#' @param weight.list A list of the same structure as adj.list, where each element contains the
#'                   weights of edges connecting vertex i to its adjacent vertices.
#' @param y A numeric vector containing function values at each vertex of the graph.
#'
#' @return A list with two components:
#' \describe{
#'   \item{lmin.hop.nbhds}{A list of hop neighborhoods for local minima. Each neighborhood
#'                        is a list containing:
#'     \itemize{
#'       \item{vertex: The vertex index (1-indexed) that is a local minimum}
#'       \item{hop.idx: The maximum hop distance at which the minimum property holds}
#'       \item{nbhd.df: A matrix with two columns - vertex indices and their hop distances from the minimum}
#'       \item{nbhd.bd.df: A matrix with two columns - boundary vertex indices and their function values}
#'     }
#'   }
#'   \item{lmax.hop.nbhds}{A list with the same structure as lmin.hop.nbhds, but for local maxima}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' adj.list <- list(c(2,3), c(1,3,4), c(1,2,4), c(2,3))
#' weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1,1))
#' y <- c(0.5, 0.2, 0.8, 0.3)
#'
#' # Calculate extrema hop neighborhoods
#' result <- compute.extrema.hop.nbhds(adj.list, weight.list, y)
#'
#' # Examine minima
#' minima <- result$lmin.hop.nbhds
#' # Examine maxima
#' maxima <- result$lmax.hop.nbhds
#'
#' # Plot the graph and highlight extrema and their hop neighborhoods
#' # (requires igraph package)
#' library(igraph)
#' g <- graph.from.adj.list(adj.list)
#' E(g)$weight <- unlist(weight.list)
#' vertex.colors <- rep("grey", length(y))
#'
#' # Highlight minima
#' for(min.nbhd in minima) {
#'   min.vertex <- min.nbhd$vertex
#'   vertex.colors[min.vertex] <- "blue"
#'   # Highlight vertices in neighborhood with decreasing intensity
#'   if(nrow(min.nbhd$nbhd.df) > 0) {
#'     for(i in 1:nrow(min.nbhd$nbhd.df)) {
#'       v <- min.nbhd$nbhd.df[i,1]
#'       if(v != min.vertex) {
#'         vertex.colors[v] <- "lightblue"
#'       }
#'     }
#'   }
#' }
#'
#' # Highlight maxima
#' for(max.nbhd in maxima) {
#'   max.vertex <- max.nbhd$vertex
#'   vertex.colors[max.vertex] <- "red"
#'   # Highlight vertices in neighborhood with decreasing intensity
#'   if(nrow(max.nbhd$nbhd.df) > 0) {
#'     for(i in 1:nrow(max.nbhd$nbhd.df)) {
#'       v <- max.nbhd$nbhd.df[i,1]
#'       if(v != max.vertex && vertex.colors[v] == "grey") {
#'         vertex.colors[v] <- "pink"
#'       }
#'     }
#'   }
#' }
#'
#' plot(g, vertex.color = vertex.colors, vertex.label = y)
#' }
#'
#' @seealso
#' \code{\link{compute.graph.basin.complex}} for computing the basin complex of extrema,
#' \code{\link{graph.watershed}} for watershed segmentation based on extrema.
#'
#' @export
compute.extrema.hop.nbhds <- function(adj.list,
                                      weight.list,
                                      y) {
    ## Input validation
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices (length of adj.list)")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call(S_compute_extrema_hop_nbhds,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    PACKAGE = "gflow")

    ## Create extrema data frame
    result$extrema_df <- create.hop.nbhd.extrema.df(result)

    ## Return the processed results
    return(result)
}

#' Create a summary data frame of extrema information
#'
#' @param res An object from compute.extrema.hop.nbhds() containing extrema_df.
#' @param edgelen.list A list where \code{edgelen.list[[i]]} contains edge lengths
#'   from vertex i to its neighbors.
#' @param y Numeric vector of outcome values at all vertices.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item label: Label for the extremum (e.g., "m1" for minima, "M1" for maxima).
#'     \item vertex: Vertex index.
#'     \item value: Function value at the extremum.
#'     \item rel.value: Relative function value (value / mean(y)).
#'     \item type: "min" for minima, "max" for maxima.
#'     \item hop.idx: Hop index of the extremum.
#'     \item p.d1: Percentile of nearest neighbor distance. Higher values (approaching 1)
#'           indicate the vertex is in a sparser region where neighbors are farther away.
#'     \item p.mean.nbrs.dist: Percentile of mean neighbor distance. Higher values
#'           indicate sparser neighborhoods.
#'   }
#'
#'   Maxima are sorted by value (descending), minima by value (ascending).
#'
#' @details
#' The inverse density measures (p.d1, p.mean.nbrs.dist) quantify local
#' sparsity. A value of 0.95 means 95% of vertices have distances less than
#' or equal to this vertex's distance, indicating this vertex is in a sparse
#' region. Isolated vertices (no neighbors) are assigned Inf distances and
#' will have inv measures of 1.0.
#'
#' @export
extrema.summary <- function(res, edgelen.list, y) {

    # Input validation
    if (!"extrema_df" %in% names(res)) {
        stop("Input 'res' must have the 'extrema_df' component")
    }

    if (length(edgelen.list) != length(y)) {
        stop("Length of edgelen.list must equal length of y")
    }

    extrema.df <- res$extrema_df

    # Initialize results data frame
    results <- data.frame(
        label = character(),
        vertex = integer(),
        value = numeric(),
        rel.value = numeric(),
        type = character(),
        hop.idx = numeric(),
        p.d1 = numeric(),
        p.mean.nbrs.dist = numeric(),
        stringsAsFactors = FALSE
    )

    # Compute distances to nearest neighbor and mean neighbor distances
    n <- length(edgelen.list)
    d1 <- numeric(n)
    mean.nbrs.dist <- numeric(n)

    for (i in seq_len(n)) {
        if (length(edgelen.list[[i]]) > 0) {
            d1[i] <- min(edgelen.list[[i]])
            mean.nbrs.dist[i] <- mean(edgelen.list[[i]])
        } else {
            # Isolated vertex: assign Inf to indicate no neighbors
            d1[i] <- Inf
            mean.nbrs.dist[i] <- Inf
        }
    }

    # Compute inverse density measures (percentiles)
    # Higher values indicate sparser regions
    p.d1 <- numeric(n)
    p.mean.nbrs.dist <- numeric(n)

    for (i in seq_len(n)) {
        # Proportion of vertices with distance <= current vertex's distance
        p.d1[i] <- sum(d1 <= d1[i]) / n
        p.mean.nbrs.dist[i] <- sum(mean.nbrs.dist <= mean.nbrs.dist[i]) / n
    }

    # Separate maxima and minima candidates
    # Get row indices, not vertex values
    max_rows <- which(extrema.df$is_max == 1)
    min_rows <- which(extrema.df$is_max == 0)

    # Process maxima
    if (length(max_rows) > 0) {
        max.vertices <- extrema.df$vertex[max_rows]

        max.df <- data.frame(
            label = "",  # Will be filled after sorting
            vertex = max.vertices,
            value = extrema.df$value[max_rows],
            rel.value = extrema.df$value[max_rows] / mean(y),
            type = "max",
            hop.idx = extrema.df$hop_idx[max_rows],
            p.d1 = p.d1[max.vertices],
            p.mean.nbrs.dist = p.mean.nbrs.dist[max.vertices],
            stringsAsFactors = FALSE
        )

        # Sort by value (descending) and assign labels
        max.df <- max.df[order(-max.df$value), ]
        max.df$label <- paste0("M", seq_len(nrow(max.df)))

        results <- rbind(results, max.df)
    }

    # Process minima
    if (length(min_rows) > 0) {
        min.vertices <- extrema.df$vertex[min_rows]

        min.df <- data.frame(
            label = "",  # Will be filled after sorting
            vertex = min.vertices,
            value = extrema.df$value[min_rows],
            rel.value = extrema.df$value[min_rows] / mean(y),
            type = "min",
            hop.idx = extrema.df$hop_idx[min_rows],
            p.d1 = p.d1[min.vertices],
            p.mean.nbrs.dist = p.mean.nbrs.dist[min.vertices],
            stringsAsFactors = FALSE
        )

        # Sort by value (ascending for minima) and assign labels
        min.df <- min.df[order(min.df$value), ]
        min.df$label <- paste0("m", seq_len(nrow(min.df)))  # lowercase 'm' for minima

        results <- rbind(results, min.df)
    }

    # Reset row names
    rownames(results) <- NULL

    return(results)
}

#' Create a Data Frame of Extrema Information
#'
#' Creates a data frame containing information about local extrema from either
#' compute.extrema.hop.nbhds() or create.gflow.cx() results.
#'
#' @param result An object from compute.extrema.hop.nbhds() or create.gflow.cx()
#' @param include_spurious Logical; whether to include spurious extrema (default TRUE)
#' @param threshold Optional threshold for hop_idx; extrema with hop_idx <= threshold
#'   are considered spurious. If NULL, uses the threshold from the result object
#'   or defaults to 2.
#' @param sort_by_value Logical; if TRUE, orders extrema by function value
#'   rather than sequentially (default FALSE)
#' @param vertex_column_name Character; name of the vertex column in the output
#'   (default "vertex")
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item vertex (or evertex): Vertex index
#'     \item hop_idx: Hop index of the extremum
#'     \item is_max: 0 for minima, 1 for maxima
#'     \item label: Label for the extremum (e.g., "min1", "max1")
#'     \item value: Function value at the extremum (if available)
#'     \item spurious: Logical; whether the extremum is considered spurious
#'   }
#'
#' @examples
#' \dontrun{
#' # With compute.extrema.hop.nbhds() result
#' hop_nbhds_result <- compute.extrema.hop.nbhds(graph$adj_list, graph$weight_list, y)
#' extrema_df <- create_extrema_df(hop_nbhds_result)
#'
#' # With create.gflow.cx() result
#' gflow_result <- create.gflow.cx(graph$adj_list, graph$weight_list, y)
#' extrema_df <- create.hop.nbhd.extrema.df(gflow_result, include_spurious = FALSE)
#' }
#'
#' @export
create.hop.nbhd.extrema.df <- function(result,
                              include_spurious = TRUE,
                              threshold = NULL,
                              sort_by_value = FALSE,
                              vertex_column_name = "vertex") {

  # Check if the input has required components
  if (!all(c("lmin_hop_nbhds", "lmax_hop_nbhds") %in% names(result))) {
    stop("Input must have 'lmin_hop_nbhds' and 'lmax_hop_nbhds' components")
  }

  # Determine if we're dealing with gflow_cx object
  is_gflow <- inherits(result, "gflow_cx")

  # Determine threshold for spurious extrema
  if (is.null(threshold)) {
    if (is_gflow) {
      threshold <- attr(result, "hop_idx_threshold")
    }
    if (is.null(threshold)) {
      threshold <- 2  # Default if not found
    }
  }

  # Extract data from minima
  lmin_data <- lapply(result$lmin_hop_nbhds, function(x) {
    list(
      vertex = x$vertex,
      hop_idx = x$hop_idx,
      value = if ("value" %in% names(x)) x$value else NA
    )
  })

  # Extract data from maxima
  lmax_data <- lapply(result$lmax_hop_nbhds, function(x) {
    list(
      vertex = x$vertex,
      hop_idx = x$hop_idx,
      value = if ("value" %in% names(x)) x$value else NA
    )
  })

  n_min <- length(lmin_data)
  n_max <- length(lmax_data)

  # Handle empty case
  if (n_min + n_max == 0) {
    return(data.frame(
      vertex = integer(),
      hop_idx = integer(),
      is_max = integer(),
      label = character(),
      value = numeric(),
      spurious = logical(),
      stringsAsFactors = FALSE
    ))
  }

  # Extract data
  vertices <- c(sapply(lmin_data, `[[`, "vertex"), sapply(lmax_data, `[[`, "vertex"))
  hop_indices <- c(sapply(lmin_data, `[[`, "hop_idx"), sapply(lmax_data, `[[`, "hop_idx"))
  is_max <- c(rep(0, n_min), rep(1, n_max))

  # Get function values if available
  values <- c(sapply(lmin_data, `[[`, "value"), sapply(lmax_data, `[[`, "value"))
  has_values <- !all(is.na(values))

  # Determine spurious extrema
  spurious <- hop_indices <= threshold

  # Filter spurious extrema if requested
  if (!include_spurious) {
    keep <- !spurious
    vertices <- vertices[keep]
    hop_indices <- hop_indices[keep]
    is_max <- is_max[keep]
    values <- values[keep]
    spurious <- spurious[keep]

    # Adjust counts after filtering
    n_min <- sum(is_max == 0)
    n_max <- sum(is_max == 1)
  }

  # Create labels
  if (sort_by_value && has_values) {
    # Label based on function values
    labels <- character(length(vertices))

    # Label maxima in descending order of function value
    if (n_max > 0) {
      max_indices <- which(is_max == 1)
      max_values <- values[max_indices]
      ord <- order(max_values, decreasing = TRUE)
      labels[max_indices[ord]] <- paste0("max", seq_len(n_max))
    }

    # Label minima in ascending order of function value
    if (n_min > 0) {
      min_indices <- which(is_max == 0)
      min_values <- values[min_indices]
      ord <- order(min_values, decreasing = FALSE)
      labels[min_indices[ord]] <- paste0("min", seq_len(n_min))
    }
  } else {
    # Sequential labeling
    labels <- c(
      paste0("min", seq_len(n_min)),
      paste0("max", seq_len(n_max))
    )
  }

  # Create data frame with appropriate column name for vertex
  df <- data.frame(
    hop_idx = hop_indices,
    is_max = is_max,
    label = labels,
    spurious = spurious,
    stringsAsFactors = FALSE
  )

  # Add function values if available
  if (has_values) {
    df$value <- values
  }

  # Add vertex column with requested name
  df[[vertex_column_name]] <- vertices

  # Reorder columns to put vertex first
  col_order <- c(vertex_column_name, "hop_idx", "is_max", "label", "spurious")
  if (has_values) {
    col_order <- c(col_order, "value")
  }

  df <- df[, col_order]

  return(df)
}

#' Label Local Extrema in 3D Plot
#'
#' Adds labeled line segments to a 3D plot at positions of local extrema
#'
#' @param graph.3d A matrix with 3 columns representing 3D coordinates
#' @param extrema.df Data frame with columns: vertex, label, and either
#'   is_max (0/1) or type ("min"/"max")
#' @param extrema.type Character string specifying which extrema to plot:
#'   "both" (default), "maxima", or "minima".
#' @param offset Numeric vector of length 3 specifying the offset for line
#'   segments. Default c(0, 0, 0.25).
#' @param with.labels Logical. Whether to show labels. Default TRUE.
#' @param lab.cex Character expansion factor for labels. Default 1.5.
#' @param lab.adj Adjustment parameter for label positioning. Default c(0, 0).
#' @param C Scaling factor for label position relative to stick center. Default -1.
#' @param pwd Line width for segments. Default 5.
#' @param separate.colors Logical. If TRUE and extrema.type="both", plot maxima
#'   and minima with different colors. Default TRUE.
#' @param col.max Color for maxima. Default "red".
#' @param col.min Color for minima. Default "blue".
#' @param col Color for all extrema when separate.colors=FALSE. Default "black".
#' @param ... Additional arguments passed to bin.segments3d().
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' # Plot both maxima and minima with different colors
#' label.extrema.3d(graph.3d, extrema.df)
#'
#' # Plot only maxima
#' label.extrema.3d(graph.3d, extrema.df, extrema.type = "maxima")
#'
#' # Plot only minima
#' label.extrema.3d(graph.3d, extrema.df, extrema.type = "minima")
#'
#' # Plot both with same color
#' label.extrema.3d(graph.3d, extrema.df, separate.colors = FALSE, col = "purple")
#' }
#'
#' @export
label.extrema.3d <- function(graph.3d,
                             extrema.df,
                             extrema.type = c("both", "maxima", "minima"),
                             offset = c(0, 0, 0.25),
                             with.labels = TRUE,
                             lab.cex = 1.5,
                             lab.adj = c(0, 0),
                             C = -1,
                             pwd = 5,
                             separate.colors = TRUE,
                             col.max = "red",
                             col.min = "blue",
                             col = "black",
                             ...) {

    # Match extrema.type argument
    extrema.type <- match.arg(extrema.type)

    # Validate graph.3d
    S <- graph.3d
    if (!is.matrix(S) && !is.data.frame(S)) {
        stop("graph.3d must be a matrix or data.frame")
    }
    if (ncol(S) != 3) {
        stop("graph.3d must have 3 columns")
    }

    # Ensure rownames exist
    if (is.null(rownames(S))) {
        rownames(S) <- as.character(1:nrow(S))
    }

    # Check which column format is present and standardize to is_max
    has_is_max <- "is_max" %in% names(extrema.df)
    has_type <- "type" %in% names(extrema.df)

    if (!has_is_max && !has_type) {
        stop("extrema.df must have either 'is_max' column (0/1) or 'type' column ('min'/'max')")
    }

    # Standardize to is_max format internally
    if (has_type && !has_is_max) {
        # Convert type to is_max
        extrema.df$is_max <- ifelse(extrema.df$type == "max", 1, 0)
    } else if (has_is_max && has_type) {
        # Both present - verify consistency and prefer is_max
        type_derived <- ifelse(extrema.df$type == "max", 1, 0)
        if (!all(extrema.df$is_max == type_derived)) {
            warning("Both 'is_max' and 'type' columns present with inconsistent values; using 'is_max'")
        }
    }
    # If only has_is_max, nothing to do

    # Filter extrema.df based on extrema.type
    if (extrema.type == "maxima") {
        extrema.df <- extrema.df[extrema.df$is_max == 1, ]
    } else if (extrema.type == "minima") {
        extrema.df <- extrema.df[extrema.df$is_max == 0, ]
    }

    # Check if we have any extrema to plot
    if (nrow(extrema.df) == 0) {
        warning("No extrema to plot based on extrema.type='", extrema.type, "'")
        return(invisible(NULL))
    }

    # Create binary vector indicating which vertices are extrema
    y <- rep(0, nrow(S))
    names(y) <- rownames(S)

    # Mark extrema vertices as 1
    extrema_vertices <- as.character(extrema.df$vertex)
    y[extrema_vertices] <- 1

    # Create label table mapping vertex IDs to labels
    lab.tbl <- extrema.df$label
    names(lab.tbl) <- as.character(extrema.df$vertex)

    # Plot based on settings
    if (extrema.type == "both" && separate.colors) {
        # Plot maxima and minima separately with different colors

        # Maxima
        is_maxima <- extrema.df$is_max == 1
        if (any(is_maxima)) {
            y_max <- rep(0, nrow(S))
            names(y_max) <- rownames(S)
            max_vertices <- as.character(extrema.df$vertex[is_maxima])
            y_max[max_vertices] <- 1

            lab.tbl.max <- extrema.df$label[is_maxima]
            names(lab.tbl.max) <- max_vertices

            bin.segments3d(S, y_max,
                          offset = offset,
                          with.labels = with.labels,
                          lab.tbl = lab.tbl.max,
                          lab.cex = lab.cex,
                          lab.adj = lab.adj,
                          C = C,
                          pwd = pwd,
                          col = col.max,
                          ...)
        }

        # Minima
        is_minima <- extrema.df$is_max == 0
        if (any(is_minima)) {
            y_min <- rep(0, nrow(S))
            names(y_min) <- rownames(S)
            min_vertices <- as.character(extrema.df$vertex[is_minima])
            y_min[min_vertices] <- 1

            lab.tbl.min <- extrema.df$label[is_minima]
            names(lab.tbl.min) <- min_vertices

            bin.segments3d(S, y_min,
                          offset = offset,
                          with.labels = with.labels,
                          lab.tbl = lab.tbl.min,
                          lab.cex = lab.cex,
                          lab.adj = lab.adj,
                          C = C,
                          pwd = pwd,
                          col = col.min,
                          ...)
        }

    } else {
        # Plot all selected extrema with the same style
        # Use appropriate color based on extrema.type
        plot.col <- col
        if (extrema.type == "maxima") {
            plot.col <- col.max
        } else if (extrema.type == "minima") {
            plot.col <- col.min
        }

        bin.segments3d(S, y,
                      offset = offset,
                      with.labels = with.labels,
                      lab.tbl = lab.tbl,
                      lab.cex = lab.cex,
                      lab.adj = lab.adj,
                      C = C,
                      pwd = pwd,
                      col = plot.col,
                      ...)
    }

    invisible(NULL)
}

#' Extract Neighborhood Component by Label ID
#'
#' Retrieves the neighborhood component from an extrema result object based on
#' the specified label identifier. The function automatically determines whether
#' to extract from local maxima or local minima neighborhoods based on the
#' extremum type.
#'
#' @param extr.obj A list object containing extrema analysis results with the
#'   following components:
#'   \itemize{
#'     \item \code{extrema_df}: A data frame with columns 'vertex', 'label',
#'       and 'is_max' (among others)
#'     \item \code{lmax_hop_nbhds}: A list of neighborhood objects for local maxima
#'     \item \code{lmin_hop_nbhds}: A list of neighborhood objects for local minima
#'   }
#' @param id A character string specifying the label of the extremum to extract
#'   (e.g., "M1", "M12", "m11", etc.). Labels starting with "M" typically indicate
#'   maxima, while "m" indicates minima.
#'
#' @return A list containing the neighborhood information for the specified extremum
#'   with the following components:
#'   \itemize{
#'     \item \code{vertex}: Integer vertex index
#'     \item \code{value}: Numeric value at the vertex
#'     \item \code{hop_idx}: Numeric hop index
#'     \item \code{nbhd_df}: Numeric matrix of neighborhood vertices
#'     \item \code{nbhd_bd_df}: Numeric matrix of neighborhood boundary vertices
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Looks up the label in \code{extrema_df} to find the corresponding vertex
#'   \item Determines whether the extremum is a maximum (is_max = 1) or minimum
#'     (is_max = 0)
#'   \item Selects the appropriate neighborhood list (\code{lmax_hop_nbhds} or
#'     \code{lmin_hop_nbhds})
#'   \item Extracts the component with the matching vertex index
#' }
#'
#' @examples
#' # Extract neighborhood for maximum M1
#' M1_nbhd <- get.nbrs(Z.extr.res, "M1")
#'
#' # Extract neighborhood for minimum m11
#' m11_nbhd <- get.nbrs(Z.extr.res, "m11")
#'
#' # Access components of the result
#' M1_nbhd$vertex      # Get vertex index
#' M1_nbhd$value       # Get value at vertex
#' M1_nbhd$nbhd_df     # Get neighborhood dataframe
#'
#' @seealso
#' Related extrema analysis functions (if applicable)
#'
#' @export
get.nbrs <- function(extr.obj, id) {
  # Find the row in extrema_df with the matching label
  row_idx <- which(extr.obj$extrema_df$label == id)

  # Check if label exists
  if (length(row_idx) == 0) {
    stop(paste("Label", id, "not found in extrema_df"))
  }

  # Get the vertex and is_max flag from extrema_df
  vertex_val <- extr.obj$extrema_df$vertex[row_idx]
  is_max <- extr.obj$extrema_df$is_max[row_idx]

  # Choose the appropriate list based on is_max
  nbhd_list <- if (is_max == 1) {
    extr.obj$lmax_hop_nbhds
  } else {
    extr.obj$lmin_hop_nbhds
  }

  # Find the component with matching vertex
  comp_idx <- which(sapply(nbhd_list, function(x) x$vertex == vertex_val))

  if (length(comp_idx) == 0) {
    stop(paste("Vertex", vertex_val, "not found in",
               ifelse(is_max == 1, "lmax_hop_nbhds", "lmin_hop_nbhds")))
  }

  return(nbhd_list[[comp_idx]])
}

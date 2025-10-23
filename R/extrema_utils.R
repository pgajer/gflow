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

#' Compute Basins of Attraction for Local Extrema
#'
#' @description
#' This function identifies all local extrema in a weighted graph and computes
#' their basins of attraction. A basin of attraction for a local extremum
#' consists of all vertices reachable via monotone paths: ascending paths for
#' minima and descending paths for maxima.
#'
#' @details
#' The basin of attraction captures the region of influence for each local
#' extremum. For a local minimum, the basin includes all vertices from which
#' the function value strictly increases along every edge of the path leading
#' back to the minimum. Similarly, for a local maximum, the basin includes
#' vertices from which the function value strictly decreases along paths to
#' the maximum.
#'
#' The algorithm performs a breadth-first search from each local extremum,
#' exploring only those edges that maintain the monotonicity property. A vertex
#' belongs to the basin if there exists at least one monotone path connecting
#' it to the extremum. The search continues until no new vertices satisfying
#' the monotonicity constraint can be reached.
#'
#' Each basin structure returned by the function contains the extremum vertex
#' index, its function value, the maximum hop distance within the basin, a
#' complete mapping of all basin vertices to their hop distances from the
#' extremum, and a boundary map identifying vertices adjacent to but outside
#' the basin.
#'
#' The function uses 0-based indexing internally for C++ compatibility but
#' returns results with 1-based indexing following R conventions.
#'
#' @param adj.list A list of integer vectors representing the graph's adjacency
#'   structure. Element \code{i} contains the indices of vertices adjacent to
#'   vertex \code{i}. Indices should be 1-based.
#' @param weight.list A list of numeric vectors containing edge weights.
#'   Element \code{i} contains weights corresponding to the edges in
#'   \code{adj.list[[i]]}. Must have the same structure as \code{adj.list}.
#' @param y A numeric vector of function values at each vertex. The length
#'   must equal the number of vertices (i.e., \code{length(adj.list)}).
#'
#' @return An object of class \code{"basins_of_attraction"} containing:
#'   \item{lmin_basins}{A list of basin structures for local minima. Each
#'     structure contains:
#'     \itemize{
#'       \item \code{vertex}: The vertex index (1-based)
#'       \item \code{value}: Function value at the vertex
#'       \item \code{hop_idx}: Maximum hop distance within the basin
#'       \item \code{basin_df}: Matrix with columns (vertex, hop_distance)
#'         for all basin members
#'       \item \code{basin_bd_df}: Matrix with columns (vertex, y_value)
#'         for boundary vertices
#'     }}
#'   \item{lmax_basins}{A list of basin structures for local maxima with
#'     the same structure as \code{lmin_basins}}
#'   \item{n_vertices}{Total number of vertices in the graph}
#'   \item{y}{Copy of the input function values}
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' adj.list <- list(
#'   c(2, 3),      # vertex 1 connects to 2, 3
#'   c(1, 3, 4),   # vertex 2 connects to 1, 3, 4
#'   c(1, 2, 4),   # vertex 3 connects to 1, 2, 4
#'   c(2, 3)       # vertex 4 connects to 2, 3
#' )
#'
#' weight.list <- list(
#'   c(1.0, 1.0),
#'   c(1.0, 1.0, 1.0),
#'   c(1.0, 1.0, 1.0),
#'   c(1.0, 1.0)
#' )
#'
#' # Function values with one minimum and one maximum
#' y <- c(5.0, 2.0, 3.0, 6.0)
#'
#' # Compute basins
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Examine structure
#' print(basins)
#'
#' # Access individual basins
#' if (length(basins$lmin_basins) > 0) {
#'   min_basin <- basins$lmin_basins[[1]]
#'   cat("Minimum at vertex:", min_basin$vertex, "\n")
#'   cat("Basin size:", nrow(min_basin$basin_df), "\n")
#' }
#' }
#'
#' @seealso \code{\link{summary.basins_of_attraction}} for generating
#'   summary statistics
#'
#' @export
compute.basins.of.attraction <- function(adj.list, weight.list, y) {
    # Input validation
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }

    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != length(adj.list)) {
        stop("Length of y must match the number of vertices")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call C++ function
    result <- .Call(S_compute_basins_of_attraction,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    PACKAGE = "gflow")

    # Add metadata
    result$n_vertices <- length(y)
    result$y <- y

    # Set class for S3 method dispatch
    class(result) <- c("basins_of_attraction", "list")

    return(result)
}

#' Summarize Basins of Attraction
#'
#' @description
#' Generates a comprehensive summary data frame of basins of attraction,
#' including basin characteristics and local density metrics for each extremum.
#'
#' @details
#' This function processes the raw basin structures returned by
#' \code{compute.basins.of.attraction} and produces a summary table suitable
#' for analysis and reporting. The summary includes both basin-specific
#' properties and local geometric characteristics of the graph at each extremum.
#'
#' The density metrics provide information about the local sparsity of the
#' graph around each extremum. The \code{d1} measure represents the distance
#' to the nearest neighbor, while \code{mean.nbrs.dist} captures the average
#' distance to all neighbors. These are converted to percentile ranks
#' (\code{p.d1} and \code{p.mean.nbrs.dist}) to facilitate comparison across
#' vertices. Higher percentile values indicate that an extremum lies in a
#' sparser region of the graph.
#'
#' The hop-k distance measure provides complementary information about isolation
#' at a specified graph distance. By examining vertices at exactly \code{hop_k}
#' steps from each extremum, we assess whether the extremum lies in a locally
#' sparse or dense region at that specific scale. This measure is particularly
#' informative when the immediate neighborhood (captured by \code{mean.nbrs.dist})
#' does not fully characterize the local geometry. The percentile rank
#' \code{p.mean.hopk.dist} enables comparison of this extended neighborhood
#' isolation across all vertices.
#'
#' The basin size quantifies the extent of the extremum's region of influence.
#' Larger basins suggest more prominent features in the function landscape.
#' The hop index indicates the maximum graph distance from the extremum to
#' any vertex within its basin, providing a measure of basin elongation.
#'
#' Extrema are labeled systematically: minima receive lowercase labels
#' (\code{m1}, \code{m2}, ...) in order of increasing function value, while
#' maxima receive uppercase labels (\code{M1}, \code{M2}, ...) in order of
#' decreasing function value. This labeling convention facilitates identification
#' of the most significant extrema.
#'
#' @param object An object of class \code{"basins_of_attraction"} returned by
#'   \code{compute.basins.of.attraction}.
#' @param adj.list A list of integer vectors representing the graph adjacency structure.
#'   Element \code{i} contains the vertex indices of neighbors of vertex \code{i}.
#'   Must have length equal to the number of vertices in the graph.
#' @param edgelen.list A list of numeric vectors containing edge lengths.
#'   Element \code{i} contains the lengths of edges incident to vertex \code{i},
#'   parallel to \code{adj.list[[i]]}. Must have length equal to the number of
#'   vertices in the graph.
#' @param hop_k Integer specifying the hop distance for computing extended
#'   neighborhood isolation. Default is 2. Must be a positive integer.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"basin_summary"} with one row per
#'   local extremum and the following columns:
#'   \item{label}{Character label for the extremum (\code{m1}, \code{m2}, ...
#'     for minima; \code{M1}, \code{M2}, ... for maxima)}
#'   \item{vertex}{Vertex index (1-based)}
#'   \item{value}{Function value at the extremum}
#'   \item{rel.value}{Function value relative to the mean (\code{value / mean(y)})}
#'   \item{type}{Extremum type: \code{"min"} or \code{"max"}}
#'   \item{hop.idx}{Maximum hop distance within the basin}
#'   \item{basin.size}{Number of vertices in the basin}
#'   \item{p.d1}{Percentile rank of nearest neighbor distance (density measure)}
#'   \item{p.mean.nbrs.dist}{Percentile rank of mean neighbor distance (density measure)}
#'   \item{p.mean.hopk.dist}{Percentile rank of mean hop-k distance (extended isolation measure)}
#'   \item{deg}{Degree of the vertex}
#'
#' @examples
#' \dontrun{
#' # Compute basins (see compute.basins.of.attraction examples)
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Create edge length list (Euclidean distances for example)
#' edgelen.list <- lapply(adj.list, function(nbrs) {
#'   rep(1.0, length(nbrs))  # uniform edge lengths
#' })
#'
#' # Generate summary with default hop_k = 2
#' basin_summary <- summary(basins, adj.list, edgelen.list)
#' print(basin_summary)
#'
#' # Generate summary with custom hop distance
#' basin_summary_h3 <- summary(basins, adj.list, edgelen.list, hop_k = 3)
#'
#' # Identify most prominent maximum
#' max_basins <- subset(basin_summary, type == "max")
#' most_prominent <- max_basins[which.max(max_basins$basin.size), ]
#' cat("Most prominent maximum:", most_prominent$label, "\n")
#' cat("Basin size:", most_prominent$basin.size, "\n")
#'
#' # Find extrema in sparse regions (both immediate and extended neighborhoods)
#' sparse_extrema <- subset(basin_summary,
#'                          p.mean.nbrs.dist > 0.75 & p.mean.hopk.dist > 0.75)
#' print(sparse_extrema[, c("label", "type", "value",
#'                          "p.mean.nbrs.dist", "p.mean.hopk.dist")])
#' }
#'
#' @seealso \code{\link{compute.basins.of.attraction}} for computing the basins
#'
#' @export
summary.basins_of_attraction <- function(object, adj.list, edgelen.list, hop_k = 2, ...) {

    if (!inherits(object, "basins_of_attraction")) {
        stop("Input must be of class 'basins_of_attraction'")
    }

    if (length(adj.list) != object$n_vertices) {
        stop("Length of adj.list must equal number of vertices")
    }

    if (length(edgelen.list) != object$n_vertices) {
        stop("Length of edgelen.list must equal number of vertices")
    }

    # Validate that adj.list and edgelen.list are parallel structures
    for (i in seq_len(object$n_vertices)) {
        if (length(adj.list[[i]]) != length(edgelen.list[[i]])) {
            stop(sprintf("Mismatch at vertex %d: adj.list has %d neighbors but edgelen.list has %d lengths",
                        i, length(adj.list[[i]]), length(edgelen.list[[i]])))
        }
    }

    if (!is.numeric(hop_k) || length(hop_k) != 1 || hop_k < 1 || hop_k != as.integer(hop_k)) {
        stop("hop_k must be a positive integer")
    }
    hop_k <- as.integer(hop_k)

    y <- object$y
    n <- object$n_vertices

    # Compute distance metrics
    d1 <- numeric(n)
    mean.nbrs.dist <- numeric(n)
    deg <- numeric(n)

    for (i in seq_len(n)) {
        if (length(edgelen.list[[i]]) > 0) {
            d1[i] <- min(edgelen.list[[i]])
            mean.nbrs.dist[i] <- mean(edgelen.list[[i]])
            deg[i] <- length(edgelen.list[[i]])
        } else {
            d1[i] <- Inf
            mean.nbrs.dist[i] <- Inf
            deg[i] <- 0
        }
    }

    # Compute hop-k distance metric
    mean.hopk.dist <- numeric(n)

    # Compute hop-k neighborhoods for each vertex using BFS
    for (i in seq_len(n)) {
        mean.hopk.dist[i] <- compute_mean_hopk_distance(i, adj.list, edgelen.list, hop_k)
    }

    # Compute density percentiles
    p.d1 <- sapply(d1, function(x) sum(d1 <= x) / n)
    p.mean.nbrs.dist <- sapply(mean.nbrs.dist, function(x)
                                sum(mean.nbrs.dist <= x) / n)
    p.mean.hopk.dist <- sapply(mean.hopk.dist, function(x)
                                sum(mean.hopk.dist <= x) / n)

    # Initialize results
    results <- data.frame(
        label = character(),
        vertex = integer(),
        value = numeric(),
        rel.value = numeric(),
        type = character(),
        hop.idx = numeric(),
        basin.size = integer(),
        p.d1 = numeric(),
        p.mean.nbrs.dist = numeric(),
        p.mean.hopk.dist = numeric(),
        deg = numeric(),
        stringsAsFactors = FALSE
    )

    # Process minima basins
    if (length(object$lmin_basins) > 0) {
        min.df <- data.frame(
            label = character(length(object$lmin_basins)),
            vertex = sapply(object$lmin_basins, function(b) b$vertex),
            value = sapply(object$lmin_basins, function(b) b$value),
            rel.value = sapply(object$lmin_basins, function(b) b$value / mean(y)),
            type = "min",
            hop.idx = sapply(object$lmin_basins, function(b) b$hop_idx),
            basin.size = sapply(object$lmin_basins, function(b) nrow(b$basin_df)),
            stringsAsFactors = FALSE
        )

        # Add density metrics
        min.df$p.d1 <- p.d1[min.df$vertex]
        min.df$p.mean.nbrs.dist <- p.mean.nbrs.dist[min.df$vertex]
        min.df$p.mean.hopk.dist <- p.mean.hopk.dist[min.df$vertex]
        min.df$deg = deg[min.df$vertex]

        # Sort by value and assign labels
        min.df <- min.df[order(min.df$value), ]
        min.df$label <- paste0("m", seq_len(nrow(min.df)))

        results <- rbind(results, min.df)
    }

    # Process maxima basins
    if (length(object$lmax_basins) > 0) {
        max.df <- data.frame(
            label = character(length(object$lmax_basins)),
            vertex = sapply(object$lmax_basins, function(b) b$vertex),
            value = sapply(object$lmax_basins, function(b) b$value),
            rel.value = sapply(object$lmax_basins, function(b) b$value / mean(y)),
            type = "max",
            hop.idx = sapply(object$lmax_basins, function(b) b$hop_idx),
            basin.size = sapply(object$lmax_basins, function(b) nrow(b$basin_df)),
            stringsAsFactors = FALSE
        )

        # Add density metrics
        max.df$p.d1 <- p.d1[max.df$vertex]
        max.df$p.mean.nbrs.dist <- p.mean.nbrs.dist[max.df$vertex]
        max.df$p.mean.hopk.dist <- p.mean.hopk.dist[max.df$vertex]
        max.df$deg = deg[max.df$vertex]

        # Sort by value and assign labels
        max.df <- max.df[order(-max.df$value), ]
        max.df$label <- paste0("M", seq_len(nrow(max.df)))

        results <- rbind(results, max.df)
    }

    rownames(results) <- NULL

    # Return with class for potential print/plot methods
    class(results) <- c("basin_summary", "data.frame")
    return(results)
}

#' Compute mean hop-k distance for a vertex
#'
#' @description
#' Helper function to compute the mean distance from a vertex to all vertices
#' at exactly hop-k distance in the graph.
#'
#' @details
#' We perform breadth-first search to identify all vertices at exactly \code{hop_k}
#' graph distance from the given vertex. For each such vertex found, we compute
#' the geodesic distance by summing edge lengths along the shortest path. The
#' function returns the mean of these distances. If no vertices exist at hop-k
#' distance, the function returns infinity to indicate isolation.
#'
#' @param vertex Integer vertex index (1-based)
#' @param adj.list List of integer vectors representing graph adjacency
#' @param edgelen.list List of edge length vectors parallel to adjacency list
#' @param hop_k Integer hop distance to query
#'
#' @return Numeric value representing mean distance to hop-k neighbors,
#'   or Inf if no vertices exist at that hop distance
#'
#' @keywords internal
compute_mean_hopk_distance <- function(vertex, adj.list, edgelen.list, hop_k) {
    n <- length(adj.list)

    # BFS initialization
    visited <- rep(FALSE, n)
    hop_distance <- rep(-1, n)
    predecessor <- rep(NA_integer_, n)

    # Queue for BFS: stores vertex indices
    queue <- integer(0)

    # Initialize with starting vertex
    visited[vertex] <- TRUE
    hop_distance[vertex] <- 0
    queue <- c(queue, vertex)

    # Perform BFS
    while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]

        # Stop if we've gone beyond hop_k
        if (hop_distance[current] >= hop_k) {
            next
        }

        # Explore neighbors
        neighbors <- adj.list[[current]]
        if (length(neighbors) > 0) {
            for (neighbor in neighbors) {
                if (!visited[neighbor]) {
                    visited[neighbor] <- TRUE
                    hop_distance[neighbor] <- hop_distance[current] + 1
                    predecessor[neighbor] <- current
                    queue <- c(queue, neighbor)
                }
            }
        }
    }

    # Find all vertices at exactly hop_k distance
    hopk_vertices <- which(hop_distance == hop_k)

    if (length(hopk_vertices) == 0) {
        return(Inf)
    }

    # Compute geodesic distances to hop-k vertices by tracing back paths
    distances <- numeric(length(hopk_vertices))

    for (i in seq_along(hopk_vertices)) {
        v <- hopk_vertices[i]
        total_dist <- 0
        current <- v

        # Trace back to source vertex
        while (!is.na(predecessor[current])) {
            prev <- predecessor[current]

            # Find edge length from prev to current
            # Need to find which neighbor index current is in prev's adjacency list
            neighbor_idx <- which(adj.list[[prev]] == current)
            if (length(neighbor_idx) > 0) {
                total_dist <- total_dist + edgelen.list[[prev]][neighbor_idx[1]]
            }

            current <- prev
        }

        distances[i] <- total_dist
    }

    return(mean(distances))
}

#' @export
print.basins_of_attraction <- function(x, ...) {
    cat("Basins of Attraction\n")
    cat("====================\n")
    cat(sprintf("Number of vertices: %d\n", x$n_vertices))
    cat(sprintf("Local minima: %d\n", length(x$lmin_basins)))
    cat(sprintf("Local maxima: %d\n", length(x$lmax_basins)))
    cat("\nUse summary() to generate detailed basin statistics\n")
    invisible(x)
}


#' Cluster Local Extrema Based on Basin Overlap
#'
#' @description
#' Identifies clusters of local extrema (either maxima or minima) whose basins
#' of attraction exhibit substantial overlap or containment relationships. The
#' function uses basin labels from \code{\link{summary.basins_of_attraction}}
#' to ensure consistency with other analyses and returns both numeric cluster
#' assignments and a convenient cluster membership list.
#'
#' @details
#' The clustering problem arises naturally in Morse-Smale analysis when multiple
#' local extrema are positioned close together in the function landscape, creating
#' basins that either overlap substantially or exhibit containment relationships
#' where a smaller basin lies entirely within a larger one. Traditional similarity
#' measures like the Jaccard index may fail to detect such containment patterns
#' because they normalize by the union size rather than the smaller set size.
#'
#' We address this by employing the overlap coefficient, also known as the
#' Szymkiewicz-Simpson index, which measures basin similarity through the ratio
#' of intersection size to the minimum basin size. For two basins $A$ and $B$,
#' the overlap coefficient is defined as
#' $$
#' \text{OC}(A,B) = \frac{|A \cap B|}{\min(|A|, |B|)}
#' $$
#' This measure equals one when either basin is completely contained in the other,
#' making it ideal for detecting nested or highly overlapping basin structures.
#'
#' The overlap distance, given by $d(A,B) = 1 - \text{OC}(A,B)$, converts this
#' similarity into a proper metric. We construct a threshold graph where basins
#' are connected by edges whenever their overlap distance falls below the
#' specified threshold. Connected components of this graph define clusters of
#' similar basins. Within each cluster, basins likely represent the same or
#' closely related features in the underlying function landscape.
#'
#' The function internally calls \code{\link{summary.basins_of_attraction}} to
#' obtain consistent basin labels. These labels follow the convention that minima
#' are labeled m1, m2, ... in order of increasing function value (so m1 is the
#' global minimum), while maxima are labeled M1, M2, ... in order of decreasing
#' function value (so M1 is the global maximum). This labeling scheme ensures
#' that clustering results align with other basin analyses.
#'
#' The function operates exclusively on basins of a single extremum type to
#' maintain the interpretability of clusters. Mixing ascending and descending
#' basins would conflate fundamentally different topological features.
#'
#' The overlap graph returned by the function can be visualized using standard
#' graph visualization tools to examine the structure of basin relationships.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}}. This object contains the
#'   complete basin structure for both local minima and local maxima.
#' @param edgelen.list A list of numeric vectors containing edge lengths.
#'   Element \code{i} contains the lengths of edges incident to vertex \code{i}.
#'   Must have length equal to the number of vertices in the graph. This parameter
#'   is passed to \code{\link{summary.basins_of_attraction}} to generate
#'   consistent basin labels and compute basin characteristics.
#' @param extrema.type Character string specifying which type of extrema to
#'   cluster. Must be either \code{"max"} for local maxima (descending basins)
#'   or \code{"min"} for local minima (ascending basins). Default is \code{"max"}.
#' @param overlap.threshold Numeric value in the interval \eqn{[0, 1]} specifying the
#'   threshold for the overlap distance below which basins are considered similar
#'   and connected in the clustering graph. Smaller values produce more conservative
#'   clustering with tighter similarity requirements. For example, a threshold of
#'   0.15 requires that basins share at least 85% overlap relative to the smaller
#'   basin. Default is 0.1, corresponding to 90% overlap.
#'
#' @return A list containing the clustering results with the following components:
#'   \describe{
#'     \item{cluster.assignments}{Named integer vector where names are basin
#'       labels from \code{summary.basins_of_attraction} (\code{M1}, \code{M2}, ...
#'       for maxima or \code{m1}, \code{m2}, ... for minima) and values are cluster
#'       identifiers. Basins with the same cluster identifier belong to the same
#'       cluster.}
#'     \item{clusters}{Named list where names are cluster identifiers (as character
#'       strings) and values are character vectors of basin labels belonging to each
#'       cluster. This provides a convenient inverse mapping from clusters to their
#'       member basins.}
#'     \item{overlap.distances}{Symmetric matrix of pairwise overlap distances
#'       between all basins. Rows and columns are labeled with basin identifiers
#'       from \code{summary.basins_of_attraction}. Diagonal entries are zero.}
#'     \item{basin.vertices}{Named list where each element contains the integer
#'       vector of vertex indices belonging to the corresponding basin. Names
#'       match basin labels from \code{summary.basins_of_attraction}.}
#'     \item{overlap.graph}{List with components \code{adj_list} and
#'       \code{weight_list} representing the overlap graph structure. In this graph,
#'       vertices correspond to basins and edges connect basins whose overlap
#'       distance is below the threshold. Edge weights are the overlap distances.
#'       Vertex names in the adjacency list match basin labels.}
#'     \item{basin.summary}{Data frame from \code{summary.basins_of_attraction}
#'       filtered to include only the extrema of the specified type. Contains
#'       detailed information about each basin including label, vertex, value,
#'       basin size, and density metrics.}
#'     \item{n.clusters}{Integer scalar giving the total number of distinct
#'       clusters identified.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins of attraction for a function on a graph
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Generate basin summary with consistent labels
#' basin.df <- summary(basins, edgelen.list)
#' print(basin.df)
#'
#' # Cluster local maxima using the same labels
#' max.clusters <- cluster.local.extrema(basins,
#'                                       edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#'
#' # Report clustering results
#' cat("Found", max.clusters$n.clusters, "clusters among",
#'     nrow(max.clusters$basin.summary), "maxima\n")
#'
#' # Examine cluster assignments
#' print(max.clusters$cluster.assignments)
#'
#' # View clusters and their members
#' print(max.clusters$clusters)
#'
#' # Identify multi-basin clusters
#' multi.basin.clusters <- Filter(function(x) length(x) > 1,
#'                                max.clusters$clusters)
#' if (length(multi.basin.clusters) > 0) {
#'   cat("\nClusters with multiple basins:\n")
#'   for (cluster.id in names(multi.basin.clusters)) {
#'     cat("Cluster", cluster.id, ":",
#'         paste(multi.basin.clusters[[cluster.id]], collapse = ", "), "\n")
#'   }
#' }
#'
#' # Inspect basin characteristics for clustered maxima
#' clustered.labels <- unlist(multi.basin.clusters)
#' clustered.basins <- max.clusters$basin.summary[
#'   max.clusters$basin.summary$label %in% clustered.labels,
#' ]
#' print(clustered.basins[, c("label", "vertex", "value", "basin.size")])
#'
#' # Check overlap distances
#' print(round(max.clusters$overlap.distances, 3))
#'
#' # Visualize the overlap graph
#' library(igraph)
#' n.basins <- length(max.clusters$overlap.graph$adj_list)
#' edge.list <- matrix(nrow = 0, ncol = 2)
#' edge.weights <- numeric(0)
#' for (i in seq_len(n.basins)) {
#'   neighbors <- max.clusters$overlap.graph$adj_list[[i]]
#'   if (length(neighbors) > 0) {
#'     for (j in seq_along(neighbors)) {
#'       neighbor <- neighbors[j]
#'       if (i < neighbor) {  # Avoid duplicate edges
#'         edge.list <- rbind(edge.list, c(i, neighbor))
#'         edge.weights <- c(edge.weights,
#'                          max.clusters$overlap.graph$weight_list[[i]][j])
#'       }
#'     }
#'   }
#' }
#' if (nrow(edge.list) > 0) {
#'   g <- graph_from_edgelist(edge.list, directed = FALSE)
#'   V(g)$name <- names(max.clusters$overlap.graph$adj_list)
#'   E(g)$weight <- edge.weights
#'   plot(g, vertex.label = V(g)$name,
#'        edge.label = round(E(g)$weight, 2),
#'        main = "Basin Overlap Graph")
#' }
#'
#' # Similarly cluster local minima
#' min.clusters <- cluster.local.extrema(basins,
#'                                       edgelen.list,
#'                                       extrema.type = "min",
#'                                       overlap.threshold = 0.15)
#' print(min.clusters$clusters)
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing basins of attraction,
#' \code{\link{summary.basins_of_attraction}} for generating basin summaries with
#' consistent labels,
#' \code{\link{create.basin.cx}} for comprehensive basin complex analysis with
#' automatic clustering and simplification
#'
#' @export
cluster.local.extrema <- function(basins.obj,
                                  edgelen.list,
                                  extrema.type = "max",
                                  overlap.threshold = 0.1) {
    ## Input validation
    if (!inherits(basins.obj, "basins_of_attraction")) {
        stop("basins.obj must be of class 'basins_of_attraction'")
    }

    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }

    if (!is.numeric(overlap.threshold) || length(overlap.threshold) != 1) {
        stop("overlap.threshold must be a single numeric value")
    }

    if (overlap.threshold < 0 || overlap.threshold > 1) {
        stop("overlap.threshold must be in the interval [0, 1]")
    }

    ## Generate basin summary to get consistent labels
    basin.summary <- summary(basins.obj, edgelen.list)

    ## Filter to selected extrema type
    if (extrema.type == "max") {
        basin.summary <- basin.summary[basin.summary$type == "max", ]
    } else {
        basin.summary <- basin.summary[basin.summary$type == "min", ]
    }

    ## Handle case with no extrema of this type
    if (nrow(basin.summary) == 0) {
        warning(paste("No", extrema.type, "basins found in basins.obj"))
        return(list(
            cluster.assignments = integer(0),
            clusters = list(),
            overlap.distances = matrix(numeric(0), 0, 0),
            basin.vertices = list(),
            overlap.graph = list(adj_list = list(), weight_list = list()),
            basin.summary = basin.summary,
            n.clusters = 0
        ))
    }

    ## Extract the appropriate basins using the summary order
    if (extrema.type == "max") {
        basin.list <- basins.obj$lmax_basins
    } else {
        basin.list <- basins.obj$lmin_basins
    }

    ## Extract vertices for each basin in summary order
    basin.vertices.list <- list()
    for (i in seq_len(nrow(basin.summary))) {
        label <- basin.summary$label[i]
        vertex <- basin.summary$vertex[i]

        ## Find the basin with this vertex
        for (basin in basin.list) {
            if (basin$vertex == vertex) {
                basin.vertices.list[[label]] <- basin$basin_df[, 1]
                break
            }
        }
    }

    ## Compute overlap distance matrix
    overlap.dists <- overlap_distance_matrix(basin.vertices.list)

    ## Create threshold graph
    overlap.graph <- create_threshold_distance_graph(
        overlap.dists,
        threshold = overlap.threshold
    )

    ## Find connected components (clusters)
    cluster.assignments <- graph.connected.components(overlap.graph$adj_list)
    names(cluster.assignments) <- rownames(overlap.dists)

    ## Create inverse mapping: clusters -> basin labels
    cluster.ids <- sort(unique(cluster.assignments))
    clusters.list <- list()
    for (cluster.id in cluster.ids) {
        basin.labels <- names(cluster.assignments)[cluster.assignments == cluster.id]
        clusters.list[[as.character(cluster.id)]] <- basin.labels
    }

    ## Return clustering information
    list(
        cluster.assignments = cluster.assignments,
        clusters = clusters.list,
        overlap.distances = overlap.dists,
        basin.vertices = basin.vertices.list,
        overlap.graph = overlap.graph,
        basin.summary = basin.summary,
        n.clusters = length(cluster.ids)
    )
}

#' Merge Clustered Basins of Attraction
#'
#' @description
#' Merges basins of attraction for local extrema that belong to the same cluster,
#' creating a simplified basin structure where each cluster is represented by a
#' single basin with its highest-valued (for maxima) or lowest-valued (for minima)
#' extremum as the representative. The function preserves complete information about
#' which extrema were merged for downstream analysis and gradient trajectory construction.
#'
#' @details
#' When multiple local extrema have highly overlapping or nested basins of attraction,
#' they often represent the same topological feature fragmented by noise or sampling
#' artifacts. This function addresses the problem by consolidating such extrema into
#' unified basin structures.
#'
#' We begin by identifying cluster representatives according to a natural ordering
#' principle. For maxima, we select the extremum with the largest function value as
#' the representative, since it corresponds to the highest point in the cluster. For
#' minima, we select the extremum with the smallest function value, corresponding to
#' the lowest point. This choice ensures that the representative captures the most
#' extreme manifestation of the topological feature.
#'
#' The merged basin is constructed as the set-theoretic union of all basins within
#' the cluster. This union represents the complete region of influence for the
#' topological feature. We preserve the basin structure of the representative extremum
#' but extend its vertex set to include all vertices from absorbed basins.
#'
#' The function maintains complete transparency about the merging operation through
#' detailed metadata. For each merged basin, we record which extrema were absorbed,
#' providing their vertex indices, labels, and function values. This information proves
#' essential for downstream analyses, particularly gradient trajectory construction,
#' where paths may terminate at absorbed extrema that no longer exist as separate
#' basins but must still be interpretable.
#'
#' The output object retains the structure of a standard basins of attraction result
#' but includes additional metadata fields documenting the simplification. This design
#' ensures compatibility with existing analysis pipelines while providing the
#' information needed to handle edge cases involving absorbed extrema.
#'
#' @param basins.obj An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{compute.basins.of.attraction}}. This object contains the original
#'   basin structure before merging.
#' @param clustering.result An object returned by \code{\link{cluster.local.extrema}}
#'   containing cluster assignments and basin metadata for the extrema type to be
#'   merged.
#' @param extrema.type Character string specifying which type of extrema were clustered.
#'   Must be either \code{"max"} for local maxima or \code{"min"} for local minima.
#'   This should match the \code{extrema.type} used in the clustering step.
#'
#' @return An object of class \code{"basins_of_attraction"} with the same structure
#'   as the input \code{basins.obj} but with clustered basins merged. The object
#'   contains:
#'   \describe{
#'     \item{lmin_basins}{List of basin structures for local minima (merged if
#'       \code{extrema.type = "min"}).}
#'     \item{lmax_basins}{List of basin structures for local maxima (merged if
#'       \code{extrema.type = "max"}).}
#'     \item{n_vertices}{Total number of vertices in the graph.}
#'     \item{y}{Copy of the input function values.}
#'     \item{merge.info}{List containing detailed information about the merging
#'       operation with components:
#'       \describe{
#'         \item{extrema.type}{Type of extrema that were merged (\code{"max"} or
#'           \code{"min"}).}
#'         \item{n.clusters}{Number of clusters identified.}
#'         \item{n.merged}{Number of clusters that involved merging (size > 1).}
#'         \item{cluster.representatives}{Named character vector mapping cluster IDs
#'           to representative basin labels.}
#'         \item{absorbed.extrema}{Data frame with one row per absorbed extremum
#'           containing columns:
#'           \describe{
#'             \item{absorbed.label}{Original label of the absorbed extremum.}
#'             \item{absorbed.vertex}{Vertex index of the absorbed extremum.}
#'             \item{absorbed.value}{Function value at the absorbed extremum.}
#'             \item{representative.label}{Label of the representative basin that
#'               absorbed this extremum.}
#'             \item{representative.vertex}{Vertex index of the representative.}
#'             \item{representative.value}{Function value at the representative.}
#'             \item{cluster.id}{Cluster identifier.}
#'           }}
#'         \item{basin.size.changes}{Data frame showing how basin sizes changed due
#'           to merging, with columns \code{label}, \code{original.size},
#'           \code{merged.size}, \code{size.increase}.}
#'       }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins of attraction
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' # Generate basin summary
#' basin.df <- summary(basins, edgelen.list)
#'
#' # Cluster local maxima
#' max.clusters <- cluster.local.extrema(basins,
#'                                       edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#'
#' # Merge clustered basins
#' merged.basins <- merge.clustered.basins(basins,
#'                                         max.clusters,
#'                                         extrema.type = "max")
#'
#' # Examine merge information
#' print(merged.basins$merge.info$absorbed.extrema)
#'
#' # See which basins grew
#' print(merged.basins$merge.info$basin.size.changes)
#'
#' # Verify the number of maxima decreased
#' cat("Original maxima:", length(basins$lmax_basins), "\n")
#' cat("Merged maxima:", length(merged.basins$lmax_basins), "\n")
#'
#' # Use merged basins for downstream analysis
#' # Note: gradient trajectories computed with original y values may
#' # terminate at absorbed extrema vertices, which can be mapped back
#' # to representatives using merge.info$absorbed.extrema
#' }
#'
#' @seealso
#' \code{\link{cluster.local.extrema}} for clustering basins,
#' \code{\link{compute.basins.of.attraction}} for computing basins,
#' \code{\link{summary.basins_of_attraction}} for basin summaries
#'
#' @export
merge.clustered.basins <- function(basins.obj,
                                   clustering.result,
                                   extrema.type) {
    ## Input validation
    if (!inherits(basins.obj, "basins_of_attraction")) {
        stop("basins.obj must be of class 'basins_of_attraction'")
    }

    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }

    ## Extract components
    basin.summary <- clustering.result$basin.summary
    clusters <- clustering.result$clusters
    cluster.assignments <- clustering.result$cluster.assignments
    basin.vertices <- clustering.result$basin.vertices

    ## Get original basin list
    if (extrema.type == "max") {
        original.basins <- basins.obj$lmax_basins
    } else {
        original.basins <- basins.obj$lmin_basins
    }

    ## Initialize merged basin list and tracking structures
    merged.basins <- list()
    cluster.representatives <- character(0)
    absorbed.extrema.df <- data.frame(
        absorbed.label = character(),
        absorbed.vertex = integer(),
        absorbed.value = numeric(),
        representative.label = character(),
        representative.vertex = integer(),
        representative.value = numeric(),
        cluster.id = character(),
        stringsAsFactors = FALSE
    )

    basin.size.changes <- data.frame(
        label = character(),
        original.size = integer(),
        merged.size = integer(),
        size.increase = integer(),
        stringsAsFactors = FALSE
    )

    ## Process each cluster
    for (cluster.id in names(clusters)) {
        basin.labels <- clusters[[cluster.id]]

        if (length(basin.labels) == 1) {
            ## Single-basin cluster: just copy the basin
            label <- basin.labels[1]
            vertex.idx <- basin.summary$vertex[basin.summary$label == label]

            ## Find original basin
            for (basin in original.basins) {
                if (basin$vertex == vertex.idx) {
                    merged.basins[[label]] <- basin
                    break
                }
            }

            cluster.representatives[cluster.id] <- label

        } else {
            ## Multi-basin cluster: merge required

            ## Find representative (max value for maxima, min value for minima)
            cluster.basin.info <- basin.summary[basin.summary$label %in% basin.labels, ]

            if (extrema.type == "max") {
                rep.idx <- which.max(cluster.basin.info$value)
            } else {
                rep.idx <- which.min(cluster.basin.info$value)
            }

            rep.label <- cluster.basin.info$label[rep.idx]
            rep.vertex <- cluster.basin.info$vertex[rep.idx]
            rep.value <- cluster.basin.info$value[rep.idx]

            cluster.representatives[cluster.id] <- rep.label

            ## Get representative's original basin
            rep.basin <- NULL
            for (basin in original.basins) {
                if (basin$vertex == rep.vertex) {
                    rep.basin <- basin
                    break
                }
            }

            ## Compute union of all basin vertices in cluster
            union.vertices <- unique(unlist(basin.vertices[basin.labels]))
            union.vertices <- sort(union.vertices)

            ## Create merged basin structure
            merged.basin <- rep.basin

            ## Update basin_df to include all union vertices
            ## Keep hop distances from representative; set others to NA or max hop + 1
            original.vertices <- rep.basin$basin_df[, 1]
            max.hop <- max(rep.basin$basin_df[, 2])

            new.vertices <- setdiff(union.vertices, original.vertices)
            if (length(new.vertices) > 0) {
                ## Add new vertices with hop distance = max.hop + 1
                new.rows <- cbind(new.vertices, rep(max.hop + 1, length(new.vertices)))
                merged.basin$basin_df <- rbind(rep.basin$basin_df, new.rows)

                ## Sort by vertex index
                merged.basin$basin_df <- merged.basin$basin_df[order(merged.basin$basin_df[, 1]), ]
            }

            ## Store merged basin with representative's label
            merged.basins[[rep.label]] <- merged.basin

            ## Record absorbed extrema
            absorbed.labels <- setdiff(basin.labels, rep.label)
            for (absorbed.label in absorbed.labels) {
                absorbed.info <- cluster.basin.info[cluster.basin.info$label == absorbed.label, ]

                absorbed.extrema.df <- rbind(absorbed.extrema.df, data.frame(
                                                                      absorbed.label = absorbed.label,
                                                                      absorbed.vertex = absorbed.info$vertex,
                                                                      absorbed.value = absorbed.info$value,
                                                                      representative.label = rep.label,
                                                                      representative.vertex = rep.vertex,
                                                                      representative.value = rep.value,
                                                                      cluster.id = cluster.id,
                                                                      stringsAsFactors = FALSE
                                                                  ))
            }

            ## Record basin size change
            original.size <- nrow(rep.basin$basin_df)
            merged.size <- nrow(merged.basin$basin_df)

            basin.size.changes <- rbind(basin.size.changes, data.frame(
                                                                label = rep.label,
                                                                original.size = original.size,
                                                                merged.size = merged.size,
                                                                size.increase = merged.size - original.size,
                                                                stringsAsFactors = FALSE
                                                            ))
        }
    }

    ## Create result object
    result <- basins.obj

    if (extrema.type == "max") {
        result$lmax_basins <- merged.basins
    } else {
        result$lmin_basins <- merged.basins
    }

    ## Add merge information
    result$merge.info <- list(
        extrema.type = extrema.type,
        n.clusters = length(clusters),
        n.merged = sum(sapply(clusters, length) > 1),
        cluster.representatives = cluster.representatives,
        absorbed.extrema = absorbed.extrema.df,
        basin.size.changes = basin.size.changes
    )

    return(result)
}

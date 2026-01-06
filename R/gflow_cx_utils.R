#' Extract Local Subgraphs Around Extrema
#'
#' @description
#' Creates a list of subgraphs centered around each local extremum in the input data.
#' Each subgraph includes vertices within the hop radius of the extremum plus an
#' additional offset.
#'
#' @param adj.list A list where each element contains the neighbors of the corresponding vertex.
#'   Must be a valid adjacency list representation of an undirected graph.
#' @param weight.list A list of numeric vectors where each element contains the edge weights
#'   corresponding to the neighbors in \code{adj.list}. Must have the same structure as
#'   \code{adj.list}.
#' @param y Numeric vector of function values on vertices. Length must equal the number of
#'   vertices in the graph.
#' @param predictions Numeric vector of smoothed predictions (e.g., from spectral filtering).
#'   Must have the same length as \code{y}.
#' @param extrema.df Data frame containing extrema information with required columns:
#'   \itemize{
#'     \item \code{vertex}: Integer vertex index
#'     \item \code{hop_idx}: Integer hop radius of the extremum
#'     \item \code{is_max}: Logical indicating if extremum is a maximum
#'     \item \code{label}: Character label for the extremum
#'     \item \code{value}: Numeric function value at the extremum
#'   }
#' @param hop.offset Integer specifying additional hop distance to include beyond the
#'   extremum's hop radius (default: 2). Must be non-negative.
#' @param min.hop.idx Integer specifying minimum hop index for extrema to process
#'   (default: 0). Must be non-negative.
#' @param max.hop.idx Integer or \code{Inf} specifying maximum hop index for extrema to
#'   process (default: \code{Inf}).
#'
#' @return A list where each element corresponds to a processed extremum and contains:
#' \describe{
#'   \item{adj_list}{Adjacency list of the subgraph (1-indexed)}
#'   \item{weight_list}{Weight list of the subgraph}
#'   \item{y}{Original function values for subgraph vertices}
#'   \item{predictions}{Predicted function values for subgraph vertices}
#'   \item{extremum_info}{Data frame row containing information about the central extremum,
#'     with vertex index mapped to subgraph numbering}
#'   \item{region_vertices}{Integer vector of subgraph indices for vertices within the
#'     extremum's hop radius}
#'   \item{boundary_vertices}{Integer vector of subgraph indices for vertices at exactly
#'     hop_idx+1 distance}
#'   \item{boundary_values}{Named numeric vector of y values at boundary vertices,
#'     names are subgraph indices}
#'   \item{vertex_mapping}{List with two components:
#'     \itemize{
#'       \item \code{orig.to.new}: Named integer vector mapping original to subgraph indices
#'       \item \code{new.to.orig}: Integer vector mapping subgraph to original indices
#'     }}
#'   \item{hop_distances}{Numeric vector of hop distances from center for all subgraph vertices}
#'   \item{extrema_df}{Data frame of all extrema within the subgraph with vertex indices
#'     mapped to subgraph numbering}
#'   \item{orig_region_vertices}{Integer vector of original indices for region vertices}
#'   \item{orig_boundary_vertices}{Integer vector of original indices for boundary vertices}
#'   \item{orig_boundary_values}{Named numeric vector of y values with original vertex indices}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' adj_list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
#' weight_list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
#' y <- c(1, 2, 3, 1)
#' predictions <- c(1.1, 2.1, 2.9, 1.1)
#'
#' # Create extrema data frame
#' extrema_df <- data.frame(
#'   vertex = 3,
#'   hop_idx = 1,
#'   is_max = TRUE,
#'   label = "Max1",
#'   fn_value = 3,
#'   stringsAsFactors = FALSE
#' )
#'
#' # Extract subgraph
#' result <- extract.extrema.subgraphs(
#'   adj_list, weight_list, y, predictions, extrema_df, hop.offset = 1
#' )
#' }
#' @export
extract.extrema.subgraphs <- function(adj.list,
                                     weight.list,
                                     y,
                                     predictions,
                                     extrema.df,
                                     hop.offset = 2,
                                     min.hop.idx = 0,
                                     max.hop.idx = Inf) {

  # Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }

  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }

  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have the same length")
  }

  n.vertices <- length(adj.list)

  if (!is.numeric(y) || length(y) != n.vertices) {
    stop("y must be a numeric vector with length equal to the number of vertices")
  }

  if (!is.numeric(predictions) || length(predictions) != n.vertices) {
    stop("predictions must be a numeric vector with length equal to the number of vertices")
  }

  if (!is.data.frame(extrema.df)) {
    stop("extrema.df must be a data frame")
  }

  required.cols <- c("vertex", "hop_idx", "is_max", "label", "value")
  missing.cols <- setdiff(required.cols, colnames(extrema.df))
  if (length(missing.cols) > 0) {
    stop("extrema.df is missing required columns: ", paste(missing.cols, collapse = ", "))
  }

  if (!is.numeric(hop.offset) || length(hop.offset) != 1 || hop.offset < 0) {
    stop("hop.offset must be a non-negative scalar")
  }

  if (!is.numeric(min.hop.idx) || length(min.hop.idx) != 1 || min.hop.idx < 0) {
    stop("min.hop.idx must be a non-negative scalar")
  }

  if (!is.numeric(max.hop.idx) || length(max.hop.idx) != 1) {
    stop("max.hop.idx must be a scalar numeric value or Inf")
  }

  # Ensure integer values where needed
  hop.offset <- as.integer(hop.offset)
  min.hop.idx <- as.integer(min.hop.idx)

  # Filter extrema by hop_idx if requested
  extrema.to.process <- extrema.df[extrema.df$hop_idx >= min.hop.idx &
                                  extrema.df$hop_idx <= max.hop.idx, ]

  if (nrow(extrema.to.process) == 0) {
    warning("No extrema found within the specified hop_idx range")
    return(list())
  }

  results <- vector("list", nrow(extrema.to.process))

  for (i in seq_len(nrow(extrema.to.process))) {
    extremum.info <- extrema.to.process[i, ]
    center.vertex <- as.integer(extremum.info$vertex)

    # Validate vertex index
    if (center.vertex < 1 || center.vertex > n.vertices) {
      warning(sprintf("Extremum at row %d has invalid vertex index %d, skipping",
                      i, center.vertex))
      next
    }

    hop.radius <- extremum.info$hop_idx
    if (is.infinite(hop.radius)) {
      # For global extrema, use a reasonable radius
      hop.radius <- 5L  # Arbitrary but reasonable value
    }
    hop.radius <- as.integer(hop.radius)

    total.radius <- hop.radius + hop.offset

    # Find all vertices within total.radius hops
    included.vertices <- find.vertices.within.hops(adj.list, center.vertex, total.radius)

    # Create mapping between original and new indices
    vertex.mapping <- list(
      orig.to.new = setNames(seq_along(included.vertices), as.character(included.vertices)),
      new.to.orig = included.vertices
    )

    # Create subgraph adjacency and weight lists
    n.subgraph <- length(included.vertices)
    subgraph.adj.list <- vector("list", n.subgraph)
    subgraph.weight.list <- vector("list", n.subgraph)

    for (v_idx in seq_along(included.vertices)) {
      orig.v <- included.vertices[v_idx]
      neighbors <- adj.list[[orig.v]]
      weights <- weight.list[[orig.v]]

      # Validate structure
      if (length(neighbors) != length(weights)) {
        stop(sprintf("Mismatch between neighbors and weights for vertex %d", orig.v))
      }

      # Only include edges to vertices that are in our subgraph
      valid.edges <- neighbors %in% included.vertices
      if (any(valid.edges)) {
        subgraph.neighbors <- vertex.mapping$orig.to.new[as.character(neighbors[valid.edges])]
        subgraph.weights <- weights[valid.edges]

        subgraph.adj.list[[v_idx]] <- as.integer(subgraph.neighbors)
        subgraph.weight.list[[v_idx]] <- as.numeric(subgraph.weights)
      } else {
        subgraph.adj.list[[v_idx]] <- integer(0)
        subgraph.weight.list[[v_idx]] <- numeric(0)
      }
    }

    # Find region vertices (within hop.radius) and boundary vertices (at hop.radius + 1)
    center.new.idx <- vertex.mapping$orig.to.new[as.character(center.vertex)]
    hop.distances <- find.hop.distances(subgraph.adj.list, center.new.idx, total.radius)

    region.indices <- which(hop.distances <= hop.radius)
    boundary.indices <- which(hop.distances == hop.radius + 1)

    region.vertices <- vertex.mapping$new.to.orig[region.indices]
    boundary.vertices <- vertex.mapping$new.to.orig[boundary.indices]

    # Handle empty boundary case
    if (length(boundary.vertices) > 0) {
      boundary.values <- setNames(y[boundary.vertices], as.character(boundary.vertices))
    } else {
      boundary.values <- numeric(0)
    }

    # Extract y and predictions for the subgraph vertices
    subgraph.y <- y[included.vertices]
    subgraph.predictions <- predictions[included.vertices]

    # Find extrema contained in this subgraph
    extrema.in.subgraph <- extrema.df[extrema.df$vertex %in% included.vertices, ]

    # Map vertex indices to subgraph numbering
    if (nrow(extrema.in.subgraph) > 0) {
      # Clone the dataframe to avoid modifying the original
      subgraph.extrema.df <- extrema.in.subgraph

      # Map vertex column to subgraph indices
      vertex.col.idx <- which(colnames(subgraph.extrema.df) == "vertex")
      if (length(vertex.col.idx) > 0) {
        subgraph.extrema.df$vertex <- sapply(
          subgraph.extrema.df$vertex,
          function(v) {
            new.idx <- vertex.mapping$orig.to.new[as.character(v)]
            if (is.null(new.idx)) {
              warning(sprintf("Vertex %d in extrema.df not found in subgraph", v))
              return(NA_integer_)
            }
            return(as.integer(new.idx))
          }
        )
      }
    } else {
      subgraph.extrema.df <- extrema.in.subgraph  # Empty dataframe
    }

    # Also map the vertices in extremum.info to subgraph indices
    subgraph.extremum.info <- extremum.info
    subgraph.extremum.info$vertex <- as.integer(
      vertex.mapping$orig.to.new[as.character(extremum.info$vertex)]
    )

    # Prepare boundary values with subgraph indices as names
    if (length(boundary.values) > 0) {
      boundary.values.subgraph <- setNames(
        boundary.values,
        sapply(names(boundary.values), function(v) {
          as.character(vertex.mapping$orig.to.new[v])
        })
      )
    } else {
      boundary.values.subgraph <- numeric(0)
    }

    # Store results
    results[[i]] <- list(
      adj_list = subgraph.adj.list,
      weight_list = subgraph.weight.list,
      y = subgraph.y,
      predictions = subgraph.predictions,
      extremum_info = subgraph.extremum.info,
      region_vertices = as.integer(sapply(region.vertices, function(v) {
        vertex.mapping$orig.to.new[as.character(v)]
      })),
      boundary_vertices = as.integer(sapply(boundary.vertices, function(v) {
        vertex.mapping$orig.to.new[as.character(v)]
      })),
      boundary_values = boundary.values.subgraph,
      vertex_mapping = vertex.mapping,
      hop_distances = hop.distances,
      extrema_df = subgraph.extrema.df,
      orig_region_vertices = region.vertices,
      orig_boundary_vertices = boundary.vertices,
      orig_boundary_values = boundary.values
    )
  }

  # Remove any NULL results (from skipped vertices)
  results <- results[!sapply(results, is.null)]

  return(results)
}


#' Find Vertices Within Specified Hop Distance
#'
#' @description
#' Uses breadth-first search to find all vertices within a specified number of hops
#' from a starting vertex.
#'
#' @param adj.list Adjacency list representation of the graph
#' @param start.vertex Integer index of the starting vertex (1-indexed)
#' @param max.hops Integer maximum number of hops to search
#'
#' @return Integer vector of vertex indices that are within max.hops of start.vertex
#'
#' @keywords internal
find.vertices.within.hops <- function(adj.list, start.vertex, max.hops) {
  n.vertices <- length(adj.list)

  # Input validation
  if (!is.numeric(start.vertex) || length(start.vertex) != 1 ||
      start.vertex < 1 || start.vertex > n.vertices) {
    stop("start.vertex must be a valid vertex index")
  }

  if (!is.numeric(max.hops) || length(max.hops) != 1 || max.hops < 0) {
    stop("max.hops must be a non-negative integer")
  }

  start.vertex <- as.integer(start.vertex)
  max.hops <- as.integer(max.hops)

  visited <- logical(n.vertices)
  visited[start.vertex] <- TRUE

  if (max.hops == 0) {
    return(start.vertex)
  }

  current.level <- start.vertex

  for (hop in seq_len(max.hops)) {
    next.level <- integer()

    for (v in current.level) {
      neighbors <- adj.list[[v]]
      if (length(neighbors) > 0) {
        unvisited <- neighbors[!visited[neighbors]]
        if (length(unvisited) > 0) {
          visited[unvisited] <- TRUE
          next.level <- c(next.level, unvisited)
        }
      }
    }

    if (length(next.level) == 0) break
    current.level <- next.level
  }

  return(which(visited))
}


#' Compute Hop Distances from Starting Vertex
#'
#' @description
#' Computes the hop distance from a starting vertex to all other vertices in the graph
#' using breadth-first search.
#'
#' @param adj.list Adjacency list representation of the graph
#' @param start.vertex Integer index of the starting vertex (1-indexed)
#' @param max.hops Integer or Inf specifying maximum hops to compute (default: Inf)
#'
#' @return Numeric vector of hop distances, with Inf for unreachable vertices
#'
#' @keywords internal
find.hop.distances <- function(adj.list, start.vertex, max.hops = Inf) {
  n.vertices <- length(adj.list)

  # Input validation
  if (!is.numeric(start.vertex) || length(start.vertex) != 1 ||
      start.vertex < 1 || start.vertex > n.vertices) {
    stop("start.vertex must be a valid vertex index")
  }

  if (!is.numeric(max.hops) || length(max.hops) != 1 || max.hops < 0) {
    stop("max.hops must be a non-negative numeric value or Inf")
  }

  start.vertex <- as.integer(start.vertex)

  distances <- rep(Inf, n.vertices)
  distances[start.vertex] <- 0

  if (max.hops == 0) {
    return(distances)
  }

  current.level <- start.vertex
  current.hop <- 0

  while (length(current.level) > 0 && current.hop < max.hops) {
    current.hop <- current.hop + 1
    next.level <- integer()

    for (v in current.level) {
      neighbors <- adj.list[[v]]
      if (length(neighbors) > 0) {
        unvisited <- neighbors[is.infinite(distances[neighbors])]
        if (length(unvisited) > 0) {
          distances[unvisited] <- current.hop
          next.level <- c(next.level, unvisited)
        }
      }
    }

    current.level <- next.level
  }

  return(distances)
}


#' Apply Harmonic Extension to Local Subgraph
#'
#' @description
#' Applies a specified harmonic extension method to a local subgraph around an extremum.
#' This function interfaces with C++ implementations of various harmonic extension algorithms.
#'
#' @param subgraph A subgraph object created by \code{extract.extrema.subgraphs}
#' @param method Character string specifying the smoothing method:
#' \describe{
#'   \item{"weighted_mean"}{Weighted mean hop disk extension}
#'   \item{"harmonic_iterative"}{Iterative harmonic extension}
#'   \item{"harmonic_eigen"}{Eigen-based harmonic extension}
#'   \item{"biharmonic_harmonic"}{Hybrid biharmonic-harmonic extension}
#'   \item{"boundary_smoothed"}{Boundary smoothed harmonic extension}
#' }
#' @param max_iterations Integer maximum number of iterations for iterative methods (default: 100)
#' @param tolerance Numeric convergence tolerance (default: 1e-6)
#' @param sigma Numeric parameter for weighted mean method (default: 1.0)
#' @param record_iterations Logical whether to record intermediate states (default: TRUE)
#' @param verbose Logical whether to print progress information (default: FALSE)
#'
#' @return A list containing:
#' \describe{
#'   \item{original}{Original function values}
#'   \item{smoothed}{Smoothed function values after harmonic extension}
#'   \item{iterations}{List of function values at each iteration (if \code{record_iterations = TRUE})}
#'   \item{method}{Method used for smoothing}
#'   \item{extremum_info}{Information about the extremum}
#'   \item{convergence_info}{List with convergence information:
#'     \itemize{
#'       \item \code{iterations_performed}: Number of iterations
#'       \item \code{final_change}: Maximum change in final iteration
#'       \item \code{converged}: Whether convergence was achieved
#'     }}
#' }
#'
#' @examples
#' \dontrun{
#' # Assumes subgraph was created by extract.extrema.subgraphs
#' result <- apply.harmonic.extension(subgraph, method = "harmonic_eigen")
#' }
#'
#' @export
apply.harmonic.extension <- function(subgraph,
                                     method = c("weighted_mean", "harmonic_iterative",
                                                "harmonic_eigen", "biharmonic_harmonic",
                                                "boundary_smoothed"),
                                     max_iterations = 100,
                                     tolerance = 1e-6,
                                     sigma = 1.0,
                                     record_iterations = TRUE,
                                     verbose = FALSE) {

  method <- match.arg(method)

  # Input validation
  if (!is.list(subgraph)) {
    stop("subgraph must be a list object created by extract.extrema.subgraphs")
  }

  required.components <- c("adj_list", "weight_list", "y", "region_vertices",
                          "boundary_values", "extremum_info")
  missing.components <- setdiff(required.components, names(subgraph))
  if (length(missing.components) > 0) {
    stop("subgraph is missing required components: ",
         paste(missing.components, collapse = ", "))
  }

  if (!is.numeric(max_iterations) || length(max_iterations) != 1 || max_iterations < 1) {
    stop("max_iterations must be a positive integer")
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0) {
    stop("tolerance must be a positive numeric value")
  }

  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("sigma must be a positive numeric value")
  }

  if (!is.logical(record_iterations) || length(record_iterations) != 1) {
    stop("record_iterations must be a single logical value")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  max_iterations <- as.integer(max_iterations)

  # Note: The actual C++ implementation would be called here
  # This is a placeholder showing the expected interface
  warning("C++ implementation not available. Returning mock results for demonstration.")

  # Mock implementation for demonstration
  n_vertices <- length(subgraph$y)
  smoothed_values <- subgraph$y

  # Simple averaging for mock demonstration
  for (i in seq_len(min(10, max_iterations))) {
    old_values <- smoothed_values
    for (v in seq_len(n_vertices)) {
      if (length(subgraph$adj_list[[v]]) > 0) {
        neighbors <- subgraph$adj_list[[v]]
        weights <- subgraph$weight_list[[v]]
        smoothed_values[v] <- sum(weights * old_values[neighbors]) / sum(weights)
      }
    }
  }

  convergence_info <- list(
    iterations_performed = 10L,
    final_change = 0.001,
    converged = TRUE
  )

  return(list(
    original = subgraph$y,
    smoothed = smoothed_values,
    iterations = if (record_iterations) list(smoothed_values) else list(),
    method = method,
    extremum_info = subgraph$extremum_info,
    convergence_info = convergence_info
  ))
}


#' Compare Multiple Harmonic Extension Methods
#'
#' @description
#' Applies multiple harmonic extension methods to a subgraph around an extremum
#' and compares their results using various metrics.
#'
#' @param subgraph A subgraph object created by \code{extract.extrema.subgraphs}
#' @param methods Character vector of method names to test (default: all available methods)
#' @param max_iterations Maximum number of iterations for iterative methods (default: 100)
#' @param tolerance Convergence tolerance (default: 1e-6)
#' @param sigma Parameter for weighted mean method (default: 1.0)
#' @param plot_results Logical whether to visualize the results (default: TRUE)
#' @param plot_type Character string specifying visualization type if \code{plot_results = TRUE}:
#'   "2d" or "3d" (default: "3d")
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{Named list of results from each method}
#'   \item{metrics}{Data frame with comparison metrics for each method:
#'     \itemize{
#'       \item \code{method}: Method name
#'       \item \code{iterations}: Number of iterations performed
#'       \item \code{final_change}: Final iteration change
#'       \item \code{converged}: Whether method converged
#'       \item \code{total_change}: Sum of absolute changes from original
#'       \item \code{extremum_removal}: Change at the extremum vertex
#'     }}
#'   \item{subgraph}{The input subgraph object}
#' }
#'
#' @examples
#' \dontrun{
#' # Assumes subgraph was created by extract.extrema.subgraphs
#' comparison <- compare.harmonic.methods(
#'   subgraph,
#'   methods = c("harmonic_iterative", "harmonic_eigen")
#' )
#' }
#'
#' @export compare.harmonic.methods
compare.harmonic.methods <- function(subgraph,
                                   methods = c("weighted_mean", "harmonic_iterative",
                                              "harmonic_eigen", "biharmonic_harmonic",
                                              "boundary_smoothed"),
                                   max_iterations = 100,
                                   tolerance = 1e-6,
                                   sigma = 1.0,
                                   plot_results = TRUE,
                                   plot_type = c("2d", "3d")) {

  # Validate inputs
  valid_methods <- c("weighted_mean", "harmonic_iterative", "harmonic_eigen",
                    "biharmonic_harmonic", "boundary_smoothed")
  methods <- match.arg(methods, valid_methods, several.ok = TRUE)
  plot_type <- match.arg(plot_type)

  if (!is.list(subgraph)) {
    stop("subgraph must be a list object created by extract.extrema.subgraphs")
  }

  # Apply each method
  results <- list()

  for (method in methods) {
    results[[method]] <- apply.harmonic.extension(
      subgraph,
      method = method,
      max_iterations = max_iterations,
      tolerance = tolerance,
      sigma = sigma,
      record_iterations = TRUE,
      verbose = FALSE
    )
  }

  # Compute comparison metrics
  metrics_list <- lapply(methods, function(method) {
    r <- results[[method]]

    # Find extremum vertex in subgraph
    extremum_idx <- which(subgraph$vertex_mapping$new.to.orig ==
                         subgraph$extremum_info$vertex)

    if (length(extremum_idx) == 0) {
      extremum_idx <- subgraph$extremum_info$vertex  # Already mapped
    }

    data.frame(
      method = method,
      iterations = ifelse(!is.null(r$convergence_info$iterations_performed),
                         r$convergence_info$iterations_performed, NA_integer_),
      final_change = ifelse(!is.null(r$convergence_info$final_change),
                           r$convergence_info$final_change, NA_real_),
      converged = ifelse(!is.null(r$convergence_info$converged),
                        r$convergence_info$converged, NA),
      total_change = sum(abs(r$smoothed - subgraph$y)),
      extremum_removal = abs(r$smoothed[extremum_idx] - subgraph$y[extremum_idx]),
      stringsAsFactors = FALSE
    )
  })

  metrics <- do.call(rbind, metrics_list)
  rownames(metrics) <- NULL

  # Plot results if requested
  if (plot_results && interactive()) {
    message("Visualization not implemented in this version")
  }

  # Return comparison results
  return(list(
    results = results,
    metrics = metrics,
    subgraph = subgraph
  ))
}


#' Analyze Harmonic Extensions for Local Extrema
#'
#' @description
#' Complete workflow for analyzing how different harmonic extension methods perform
#' on local extrema in a graph function. This function extracts subgraphs around
#' each extremum and applies various harmonic extension methods for comparison.
#'
#' @param adj_list Original graph adjacency list
#' @param weight_list Original graph weight list
#' @param y Original function values on vertices
#' @param predictions Smoothed predictions (e.g., from spectral filtering)
#' @param extrema_df Data frame containing extrema information
#' @param hop_offset Additional hop distance beyond extremum's hop radius (default: 2)
#' @param min_hop_idx Minimum hop index for extrema to process (default: 0)
#' @param max_hop_idx Maximum hop index for extrema to process (default: Inf)
#' @param methods Vector of harmonic extension methods to test
#' @param max_iterations Maximum iterations for iterative methods (default: 100)
#' @param tolerance Convergence tolerance (default: 1e-6)
#' @param sigma Parameter for weighted mean method (default: 1.0)
#' @param visualize Whether to visualize results (default: TRUE)
#' @param save_results Whether to save detailed results (default: FALSE)
#' @param output_dir Directory for saving results (default: "harmonic_extension_analysis")
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \describe{
#'   \item{subgraphs}{List of extracted subgraphs for each extremum}
#'   \item{results}{List of results from applying harmonic extension methods}
#'   \item{metrics}{Data frame with detailed metrics for all extrema and methods}
#'   \item{summary}{Data frame with summary statistics across methods}
#'   \item{extrema_info}{Data frame with information about processed extrema}
#' }
#'
#' @examples
#' \dontrun{
#' result <- analyze.harmonic.extensions(
#'   graph$adj_list,
#'   graph$weight_list,
#'   y,
#'   gsf_res$predictions,
#'   gsf_b_cx$extrema_df,
#'   hop_offset = 2,
#'   min_hop_idx = 1,
#'   max_hop_idx = 5
#' )
#' }
#'
#' @export
analyze.harmonic.extensions <- function(adj_list,
                                       weight_list,
                                       y,
                                       predictions,
                                       extrema_df,
                                       hop_offset = 2,
                                       min_hop_idx = 0,
                                       max_hop_idx = Inf,
                                       methods = c("weighted_mean", "harmonic_iterative",
                                                  "harmonic_eigen", "biharmonic_harmonic",
                                                  "boundary_smoothed"),
                                       max_iterations = 100,
                                       tolerance = 1e-6,
                                       sigma = 1.0,
                                       visualize = TRUE,
                                       save_results = FALSE,
                                       output_dir = "harmonic_extension_analysis",
                                       verbose = TRUE) {

  # Create output directory if needed
  if (save_results && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Extract subgraphs around each extremum
  if (verbose) {
    message("Extracting subgraphs for extrema...")
  }

  extrema_subgraphs <- extract.extrema.subgraphs(
    adj_list,
    weight_list,
    y,
    predictions,
    extrema_df,
    hop.offset = hop_offset,
    min.hop.idx = min_hop_idx,
    max.hop.idx = max_hop_idx
  )

  n_extrema <- length(extrema_subgraphs)

  if (n_extrema == 0) {
    warning("No extrema found to process")
    return(list(
      subgraphs = list(),
      results = list(),
      metrics = data.frame(),
      summary = data.frame(),
      extrema_info = data.frame()
    ))
  }

  if (verbose) {
    message("Found ", n_extrema, " extrema to process")
  }

  # Apply harmonic extension methods to each subgraph
  all_results <- vector("list", n_extrema)
  comparison_metrics <- vector("list", n_extrema)

  for (i in seq_len(n_extrema)) {
    subgraph <- extrema_subgraphs[[i]]
    extremum_info <- subgraph$extremum_info

    if (verbose) {
      message(sprintf("\nProcessing extremum %d of %d - %s at vertex %d with hop_idx = %d",
                     i, n_extrema, extremum_info$label,
                     subgraph$orig_region_vertices[1],  # Use original vertex index
                     extremum_info$hop_idx))
    }

    # Apply methods and compare
    comparison <- compare.harmonic.methods(
      subgraph,
      methods = methods,
      max_iterations = max_iterations,
      tolerance = tolerance,
      sigma = sigma,
      plot_results = FALSE  # Disable plotting in batch processing
    )

    all_results[[i]] <- comparison$results

    # Add extremum index to metrics
    comparison$metrics$extremum_idx <- i
    comparison_metrics[[i]] <- comparison$metrics

    # Save detailed results if requested
    if (save_results) {
      extremum_dir <- file.path(output_dir,
                               sprintf("extremum_%d", subgraph$orig_region_vertices[1]))
      dir.create(extremum_dir, recursive = TRUE, showWarnings = FALSE)

      # Save metrics
      write.csv(comparison$metrics,
               file.path(extremum_dir, "method_comparison.csv"),
               row.names = FALSE)

      # Save original and smoothed values
      values_df <- data.frame(
        vertex_original = subgraph$vertex_mapping$new.to.orig,
        vertex_subgraph = seq_along(subgraph$y),
        original = subgraph$y
      )

      for (method in names(comparison$results)) {
        values_df[[method]] <- comparison$results[[method]]$smoothed
      }

      write.csv(values_df,
               file.path(extremum_dir, "function_values.csv"),
               row.names = FALSE)
    }
  }

  # Combine comparison metrics
  all_metrics <- do.call(rbind, comparison_metrics)
  rownames(all_metrics) <- NULL

  # Create extremum information data frame
  extrema_info_list <- lapply(seq_len(n_extrema), function(i) {
    extremum <- extrema_subgraphs[[i]]$extremum_info
    orig_vertex <- extrema_subgraphs[[i]]$orig_region_vertices[1]
    data.frame(
      extremum_idx = i,
      vertex = orig_vertex,
      label = extremum$label,
      hop_idx = extremum$hop_idx,
      is_max = extremum$is_max,
      value = extremum$fn_value,
      stringsAsFactors = FALSE
    )
  })
  extrema_info <- do.call(rbind, extrema_info_list)

  # Create summary statistics
  summary_stats <- aggregate(
    cbind(iterations, final_change, converged, total_change, extremum_removal) ~ method,
    data = all_metrics,
    FUN = function(x) {
      if (is.logical(x)) {
        mean(x, na.rm = TRUE)
      } else {
        mean(x, na.rm = TRUE)
      }
    }
  )

  # Rename converged column to pct_converged
  names(summary_stats)[names(summary_stats) == "converged"] <- "pct_converged"
  summary_stats$pct_converged <- summary_stats$pct_converged * 100

  # Save overall comparison if requested
  if (save_results) {
    write.csv(summary_stats,
             file.path(output_dir, "overall_comparison.csv"),
             row.names = FALSE)

    write.csv(all_metrics,
             file.path(output_dir, "detailed_metrics.csv"),
             row.names = FALSE)

    write.csv(extrema_info,
             file.path(output_dir, "extrema_info.csv"),
             row.names = FALSE)
  }

  # Return results
  return(list(
    subgraphs = extrema_subgraphs,
    results = all_results,
    metrics = all_metrics,
    summary = summary_stats,
    extrema_info = extrema_info
  ))
}

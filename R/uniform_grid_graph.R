#' Create a Refined Graph with Uniformly Spaced Grid Vertices
#'
#' @description
#' Creates a refined version of an input graph by adding grid vertices (points) along
#' its edges. The grid vertices are placed to maintain approximately uniform spacing
#' throughout the graph structure. This function is particularly useful for tasks
#' that require a denser sampling of points along the graph edges, such as
#' graph-based interpolation or spatial analysis.
#'
#' @param adj.list A list where each element \code{i} is an integer vector containing the
#'   indices of vertices adjacent to vertex \code{i}. Vertex indices must be 1-based
#'   (following R's convention). The graph structure must be undirected, meaning if
#'   vertex \code{j} appears in \code{adj.list[[i]]}, then vertex \code{i} must appear
#'   in \code{adj.list[[j]]}.
#' @param weight.list A list matching the structure of \code{adj.list}, where each element
#'   contains the corresponding edge weights (typically distances or lengths).
#'   \code{weight.list[[i]][j]} should contain the weight of the edge between vertex
#'   \code{i} and vertex \code{adj.list[[i]][j]}. Weights must be positive numbers.
#' @param grid.size A positive integer specifying the desired number of grid vertices
#'   to add. Must be at least 2. The actual number of added vertices may differ
#'   slightly from this target due to the distribution of edge lengths in the graph.
#' @param start.vertex An integer specifying the starting vertex for graph traversal.
#'   Must be between 1 and the number of vertices in the input graph. Defaults to 1.
#' @param snap.tolerance A numeric value between 0 and 0.5 controlling the snapping
#'   behavior when placing grid vertices near existing vertices. When a grid vertex
#'   would be placed within this fraction of an edge length from an existing vertex,
#'   it is merged with that vertex instead. Defaults to 0.1.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{adj_list}}{A list representing the adjacency structure of the refined
#'     graph, including both original and grid vertices. The structure follows the same
#'     format as the input \code{adj.list}.}
#'   \item{\code{weight_list}}{A list containing edge weights for the refined graph,
#'     structured to match the new \code{adj.list}. Edge weights are adjusted to
#'     reflect the new distances between connected vertices.}
#'   \item{\code{grid_vertices}}{An integer vector containing the 1-based indices of
#'     the newly added grid vertices in the refined graph.}
#' }
#'
#' @details
#' The function performs a breadth-first traversal starting from \code{start.vertex}
#' to determine the order in which edges are processed. Grid vertices are distributed
#' across edges based on their relative lengths, with longer edges receiving more grid
#' vertices. The algorithm ensures that the spacing between consecutive vertices
#' (original or grid) along any edge is as uniform as possible.
#'
#' The \code{snap.tolerance} parameter helps prevent the creation of vertices that are
#' too close to existing ones, which can cause numerical issues in downstream analyses.
#'
#' @examples
#' # Create a simple path graph with 3 vertices
#' adj <- list(c(2L), c(1L, 3L), c(2L))
#' weights <- list(c(1.0), c(1.0, 2.0), c(2.0))
#'
#' # Add approximately 5 grid vertices
#' result <- create.grid.graph(adj, weights, grid.size = 5)
#'
#' # Examine the results
#' cat("Original graph had", length(adj), "vertices\\n")
#' cat("Refined graph has", length(result$adj.list), "vertices\\n")
#' cat("Grid vertices added:", length(result$grid.vertices), "\\n")
#' cat("Indices of grid vertices:", result$grid.vertices, "\\n")
#'
#' \dontrun{
#' # Create a more complex graph (a cycle with varying edge weights)
#' n <- 6
#' adj <- lapply(1:n, function(i) c(ifelse(i == 1, n, i - 1),
#'                                  ifelse(i == n, 1, i + 1)))
#' weights <- lapply(1:n, function(i) runif(2, 0.5, 2.0))
#'
#' # Refine with more grid vertices
#' result <- create.grid.graph(adj, weights, grid.size = 20)
#'
#' g <- igraph::graph_from_adj_list(result$adj.list, mode = "undirected")
#' plot(g, vertex.color = ifelse(1:length(result$adj.list) %in%
#'                                 result$grid.vertices, "red", "blue"))
#' }
#'
#' @seealso
#' \code{\link[igraph]{graph_from_adj_list}} for converting adjacency lists to igraph objects
#'
#' @export
create.grid.graph <- function(adj.list,
                              weight.list,
                              grid.size,
                              start.vertex = 1L,
                              snap.tolerance = 0.1) {

  # Input validation with informative error messages

  # Check that inputs are lists
  if (!is.list(adj.list)) {
    stop("'adj.list' must be a list", call. = FALSE)
  }

  if (!is.list(weight.list)) {
    stop("'weight.list' must be a list", call. = FALSE)
  }

  # Get number of vertices and check minimum requirement
  n.vertices <- length(adj.list)

  if (n.vertices < 2L) {
    stop("Graph must have at least 2 vertices (adj.list length >= 2)",
         call. = FALSE)
  }

  if (length(weight.list) != n.vertices) {
    stop("'adj.list' and 'weight.list' must have the same length",
         call. = FALSE)
  }

  # Validate grid.size
  if (!is.numeric(grid.size) || length(grid.size) != 1L) {
    stop("'grid.size' must be a single numeric value", call. = FALSE)
  }

  grid.size <- as.integer(grid.size)
  if (is.na(grid.size) || grid.size < 2L) {
    stop("'grid.size' must be an integer >= 2", call. = FALSE)
  }

  # Validate start.vertex
  if (!is.numeric(start.vertex) || length(start.vertex) != 1L) {
    stop("'start.vertex' must be a single numeric value", call. = FALSE)
  }

  start.vertex <- as.integer(start.vertex)
  if (is.na(start.vertex) || start.vertex < 1L || start.vertex > n.vertices) {
    stop("'start.vertex' must be an integer between 1 and ", n.vertices,
         call. = FALSE)
  }

  # Validate snap.tolerance
  if (!is.numeric(snap.tolerance) || length(snap.tolerance) != 1L) {
    stop("'snap.tolerance' must be a single numeric value", call. = FALSE)
  }

  if (is.na(snap.tolerance) || snap.tolerance < 0 || snap.tolerance > 0.5) {
    stop("'snap.tolerance' must be between 0 and 0.5", call. = FALSE)
  }

  # Validate structure of adj.list and weight.list elements
  for (i in seq_len(n.vertices)) {
    # Check adj.list[[i]]
    if (!is.numeric(adj.list[[i]])) {
      stop("adj.list[[", i, "]] must be a numeric vector", call. = FALSE)
    }

    # Convert to integer and check for valid vertex indices
    adj.list[[i]] <- as.integer(adj.list[[i]])
    if (any(is.na(adj.list[[i]]))) {
      stop("adj.list[[", i, "]] contains non-integer values", call. = FALSE)
    }

    if (length(adj.list[[i]]) > 0L) {
      if (any(adj.list[[i]] < 1L | adj.list[[i]] > n.vertices)) {
        stop("adj.list[[", i, "]] contains vertex indices outside the valid range [1, ",
             n.vertices, "]", call. = FALSE)
      }
    }

    # Check weight.list[[i]]
    if (!is.numeric(weight.list[[i]])) {
      stop("weight.list[[", i, "]] must be a numeric vector", call. = FALSE)
    }

    if (any(is.na(weight.list[[i]]))) {
      stop("weight.list[[", i, "]] contains NA values", call. = FALSE)
    }

    if (any(weight.list[[i]] <= 0)) {
      stop("weight.list[[", i, "]] contains non-positive values; all edge weights must be positive",
           call. = FALSE)
    }

    # Check matching lengths
    if (length(adj.list[[i]]) != length(weight.list[[i]])) {
      stop("Length mismatch at vertex ", i, ": adj.list[[", i, "]] has length ",
           length(adj.list[[i]]), " but weight.list[[", i, "]] has length ",
           length(weight.list[[i]]), call. = FALSE)
    }
  }

  # Optional: Check for graph symmetry (undirected property)
  # This is commented out in the original code, but could be enabled
  if (FALSE) {  # Set to TRUE to enable symmetry checking
    for (i in seq_len(n.vertices)) {
      for (j in seq_along(adj.list[[i]])) {
        neighbor <- adj.list[[i]][j]
        # Find if vertex i appears in neighbor's adjacency list
        idx <- which(adj.list[[neighbor]] == i)
        if (length(idx) == 0L) {
          stop("Graph is not undirected: edge (", i, " -> ", neighbor,
               ") exists but edge (", neighbor, " -> ", i, ") does not",
               call. = FALSE)
        }
        # Also check if weights match (for undirected graphs)
        neighbor.weight.to.i <- weight.list[[neighbor]][idx[1L]]
        if (abs(weight.list[[i]][j] - neighbor.weight.to.i) > .Machine$double.eps^0.5) {
          warning("Edge weight mismatch for undirected edge between vertices ",
                  i, " and ", neighbor, ": ", weight.list[[i]][j], " != ",
                  neighbor.weight.to.i, call. = FALSE)
        }
      }
    }
  }

  # Convert to 0-based indexing for C++ (R uses 1-based indexing)
  adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1L))

  # Call the C++ implementation
  # The C++ function expects 0-based indices
  result <- .Call(S_create_uniform_grid_graph,
                  adj.list.0based,
                  weight.list,
                  as.integer(grid.size),
                  as.integer(start.vertex - 1L),
                  as.numeric(snap.tolerance))

  return(result)
}

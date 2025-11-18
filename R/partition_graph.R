#' Compute Partition Graph from Vertex Partition
#'
#' Given a graph and a partition of its vertices, this function constructs an
#' aggregated graph where vertices represent partition cells and edges capture
#' connectivity patterns between cells. Two cells are connected if at least one
#' vertex in the first cell has a neighbor in the second cell.
#'
#' The function operates on the geometric structure of the partition, treating
#' each cell as a macro-vertex and aggregating edge information across the
#' partition boundary. This operation is fundamental in hierarchical graph
#' analysis and multi-scale geometric data analysis.
#'
#' @param adj.list List of integer vectors representing the adjacency list of
#'   the graph. Each element \code{adj.list[[i]]} contains the neighbor indices
#'   of vertex \code{i} (1-based indexing).
#' @param weight.list List of numeric vectors containing edge weights. Each
#'   element \code{weight.list[[i]]} contains weights corresponding to edges in
#'   \code{adj.list[[i]]}.
#' @param partition Integer vector of length equal to the number of vertices,
#'   where \code{partition[i]} indicates the cell label to which vertex \code{i}
#'   belongs. Cell labels can be any integers and need not be contiguous.
#' @param weight.type Character string specifying how edge weights should be
#'   computed in the partition graph. Options are:
#'   \itemize{
#'     \item \code{"count"}: Raw count of edges connecting the two cells
#'     (default).
#'     \item \code{"normalized"}: Edge count normalized by the geometric mean of
#'     the two cell sizes, providing a scale-invariant measure.
#'     \item \code{"jaccard"}: Jaccard index based on boundary vertices,
#'     measuring the overlap of cells' neighborhoods relative to their union.
#'   }
#'
#' @return A list with two components:
#'   \describe{
#'     \item{adj.list}{Adjacency list for the partition graph, where each cell
#'       is represented by an index from 1 to the number of unique partition
#'       labels.}
#'     \item{weight.list}{Edge weights for the partition graph according to the
#'       specified \code{weight.type}.}
#'   }
#'
#' @details
#' The partition graph provides a coarse-grained view of connectivity patterns
#' across a partition. For each pair of cells C and C', the function counts all
#' edges with one endpoint in C and the other in C'. The raw edge count can be
#' normalized in several ways to account for cell sizes and provide
#' interpretable edge weights.
#'
#' The normalization options serve different purposes. The count option
#' preserves the total edge weight information, which may be useful when cells
#' have similar sizes. The normalized option divides by the geometric mean of
#' cell sizes, making edge weights comparable across cell pairs of different
#' sizes. The Jaccard option provides a similarity measure based on the overlap
#' of boundary vertices, which can be interpreted as the fraction of shared
#' connectivity between cells.
#'
#' The resulting partition graph can be used for hierarchical analysis,
#' community detection at the cell level, or as input to further geometric
#' analysis methods. Cell indices in the output are assigned as contiguous
#' integers from 1 to the number of cells, ordered by the sorted unique values
#' of the input partition labels.
#'
#' @examples
#' \dontrun{
#' ## Create a simple graph with 6 vertices
#' adj.list <- list(
#'   c(2, 3),    # vertex 1
#'   c(1, 3),    # vertex 2
#'   c(1, 2, 4), # vertex 3
#'   c(3, 5, 6), # vertex 4
#'   c(4, 6),    # vertex 5
#'   c(4, 5)     # vertex 6
#' )
#' weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
#'
#' ## Define a partition: cells {1,2,3} and {4,5,6}
#' partition <- c(1, 1, 1, 2, 2, 2)
#'
#' ## Compute partition graph with different weight types
#' pg.count <- partition.graph(adj.list, weight.list, partition, "count")
#' pg.norm <- partition.graph(adj.list, weight.list, partition, "normalized")
#' pg.jacc <- partition.graph(adj.list, weight.list, partition, "jaccard")
#' }
#'
#' @export
partition.graph <- function(adj.list, weight.list, partition,
                           weight.type = c("count", "normalized", "jaccard")) {

  ## Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }

  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }

  n.vertices <- length(adj.list)

  if (length(weight.list) != n.vertices) {
    stop("adj.list and weight.list must have the same length")
  }

  if (!is.numeric(partition) && !is.integer(partition)) {
    stop("partition must be a numeric or integer vector")
  }

  if (length(partition) != n.vertices) {
    stop("partition must have length equal to number of vertices")
  }

  if (any(is.na(partition))) {
    stop("partition contains NA values")
  }

  ## Match weight type argument
  weight.type <- match.arg(weight.type)

  ## Convert partition to integer
  partition <- as.integer(partition)

  ## Convert adjacency list to 0-based indexing for C++
  adj.list.0 <- lapply(adj.list, function(x) {
    if (length(x) == 0) return(integer(0))
    as.integer(x - 1)
  })

  ## Ensure weight list is numeric
  weight.list.numeric <- lapply(weight.list, function(x) {
    if (length(x) == 0) return(numeric(0))
    as.numeric(x)
  })

  ## Validate matching lengths
  for (i in seq_along(adj.list.0)) {
    if (length(adj.list.0[[i]]) != length(weight.list.numeric[[i]])) {
      stop(sprintf("Mismatch in lengths for vertex %d: adj.list has %d neighbors, weight.list has %d weights",
                   i, length(adj.list.0[[i]]), length(weight.list.numeric[[i]])))
    }
  }

  ## Call C++ function
  result <- .Call(S_partition_graph,
                  adj.list.0,
                  weight.list.numeric,
                  partition,
                  weight.type,
                  PACKAGE = "gflow")

  ## Add metadata
  n.cells <- length(unique(partition))
  attr(result, "n_vertices") <- n.vertices
  attr(result, "n_cells") <- n.cells
  attr(result, "weight_type") <- weight.type

  ## Assign class for method dispatch
  class(result) <- c("partition_graph", "wgraph", "list")

  return(result)
}

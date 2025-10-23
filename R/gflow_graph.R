#' Construct Gradient Flow Graph from Merged Basins
#'
#' @description
#' Constructs the gradient flow graph, which is the 1-skeleton of the nerve
#' simplicial complex associated with a collection of merged basins of attraction.
#' This graph captures the topological structure of how basins cover the point set
#' through their vertex intersections.
#'
#' @details
#' Consider a collection of basins, each representing a region of the graph where
#' the function exhibits monotonic behavior toward a local extremum. These basins
#' form a covering of the vertex set, though not necessarily a disjoint one. Two
#' basins may share vertices, particularly near boundaries where the gradient flow
#' exhibits ambiguous behavior or where noise in the data creates overlapping
#' influence regions.
#'
#' We begin with a fundamental construction from algebraic topology. Given any
#' covering of a space by open sets, the nerve complex records which collections
#' of sets have nonempty intersection. A vertex in the nerve corresponds to each
#' set in the covering, and a simplex appears whenever the corresponding sets
#' intersect. This construction preserves essential topological information about
#' the covering and, under suitable conditions, about the underlying space itself.
#'
#' For our gradient flow basins, the nerve complex takes a particularly natural
#' form. Each basin becomes a vertex in the complex. An edge connects two vertices
#' when the corresponding basins share at least one point from the original graph.
#' Higher dimensional simplices would appear when three or more basins have common
#' vertices, but for the gradient flow graph we focus only on the 1-skeleton, which
#' captures the essential connectivity structure.
#'
#' The motivation for studying this 1-skeleton becomes clear when we consider the
#' problem of understanding function structure on discrete graphs. Local extrema and
#' their basins fragment the space into regions of monotonic behavior. The gradient
#' flow graph reveals how these fragments connect, which basins border each other,
#' and ultimately the large-scale structure of the function landscape. An edge in
#' this graph indicates that gradient flow from different extrema can reach a common
#' region, suggesting either a saddle point connection or a region where multiple
#' basins compete for influence.
#'
#' The construction proceeds in three stages. First, we extract the vertex sets from
#' all basins, both ascending and descending. Second, we compute pairwise intersections
#' between all basin vertex sets to identify which pairs share vertices. Third, we
#' build the graph adjacency structure, creating an edge between basins whenever their
#' intersection is nonempty. The resulting graph encodes the combinatorial structure
#' of the basin covering, with edge weights reflecting the size of shared vertex sets.
#'
#' This graph serves multiple purposes in data analysis. It provides a discrete
#' representation of the Morse-Smale complex and enables downstream regression
#' and classification tasks that respect the structure of the feature space. The
#' vertex set of the gradient flow graph consists of both ascending and
#' descending basins, while edges may connect any pair of basins regardless of
#' type, capturing the full structure of how different flow regions interact.
#'
#' @param merged.basins An object of class \code{"basins_of_attraction"} returned by
#'   \code{\link{merge.clustered.basins}} containing merged basin structures. This
#'   object should have both ascending and descending basins merged according to their
#'   clustering. The function works with the basin structure after both maximum and
#'   minimum clustering and merging have been performed.
#'
#' @param min.intersection Integer specifying the minimum number of shared vertices
#'   required for an edge to appear in the gradient flow graph. Default is 1, meaning
#'   any nonempty intersection creates an edge. Larger values filter out weak
#'   connections and focus on more substantial basin overlaps.
#'
#' @param edge.type Character string specifying which types of edges to include in
#'   the graph. Options are:
#'   \describe{
#'     \item{\code{"all"}}{Include all edges between any pair of basins with sufficient
#'       intersection (default).}
#'     \item{\code{"ms_only"}}{Include only Morse-Smale edges connecting ascending
#'       (minimum) basins to descending (maximum) basins, creating a bipartite graph
#'       that captures saddle connections.}
#'     \item{\code{"same_type"}}{Include only edges between basins of the same type,
#'       either both ascending or both descending, revealing competition within extrema
#'       of the same orientation.}
#'   }
#'
#' @return An object of class \code{"gflow_graph"} containing:
#'   \describe{
#'     \item{adjacency.list}{List of integer vectors representing the graph adjacency
#'       structure. Each element corresponds to a basin vertex and contains indices of
#'       adjacent basins in the gradient flow graph. Indices refer to positions in the
#'       combined vertex list (ascending basins followed by descending basins).}
#'     \item{weight.list}{List of numeric vectors representing edge weights. Each weight
#'       equals the number of shared vertices between the corresponding basin pair,
#'       normalized by the size of the smaller basin to give a similarity measure.}
#'     \item{intersection.matrix}{Symmetric matrix where entry (i,j) records the number
#'       of vertices shared between basin i and basin j. Diagonal entries give basin
#'       sizes. This matrix provides complete information about pairwise basin overlaps.}
#'     \item{basin.metadata}{Data frame with one row per basin vertex containing:
#'       \describe{
#'         \item{index}{Sequential index in the combined vertex list.}
#'         \item{label}{Basin label (e.g., "m1", "M2") from the merged basin structure.}
#'         \item{type}{Basin type, either "ascending" for minima or "descending" for maxima.}
#'         \item{size}{Number of vertices in the basin.}
#'         \item{extremum.vertex}{1-based vertex index of the basin's extremum.}
#'         \item{extremum.value}{Function value at the extremum.}
#'       }}
#'     \item{n.ascending}{Number of ascending basin vertices in the graph.}
#'     \item{n.descending}{Number of descending basin vertices in the graph.}
#'     \item{edge.type}{Character string recording which edge type was used.}
#'     \item{min.intersection}{Integer recording the minimum intersection threshold used.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compute basins and cluster extrema
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' edgelen.list <- compute.edge.lengths(adj.list, weight.list)
#'
#' # Cluster and merge maxima
#' max.clusters <- cluster.local.extrema(basins, edgelen.list,
#'                                       extrema.type = "max",
#'                                       overlap.threshold = 0.15)
#' merged.max.basins <- merge.clustered.basins(basins, max.clusters,
#'                                             extrema.type = "max")
#'
#' # Cluster and merge minima
#' min.clusters <- cluster.local.extrema(merged.max.basins, edgelen.list,
#'                                       extrema.type = "min",
#'                                       overlap.threshold = 0.15)
#' merged.basins <- merge.clustered.basins(merged.max.basins, min.clusters,
#'                                         extrema.type = "min")
#'
#' # Construct gradient flow graph
#' gflow.graph <- construct.gflow.graph(merged.basins)
#'
#' # Examine the structure
#' print(gflow.graph)
#' cat("Number of vertices:", length(gflow.graph$adjacency.list), "\n")
#' cat("Number of edges:", sum(sapply(gflow.graph$adjacency.list, length)) / 2, "\n")
#'
#' # Find which basins connect to a specific basin
#' basin.idx <- 1
#' neighbors <- gflow.graph$adjacency.list[[basin.idx]]
#' neighbor.labels <- gflow.graph$basin.metadata$label[neighbors]
#' cat("Basin", gflow.graph$basin.metadata$label[basin.idx],
#'     "connects to:", paste(neighbor.labels, collapse = ", "), "\n")
#'
#' # Construct only Morse-Smale edges
#' ms.graph <- construct.gflow.graph(merged.basins, edge.type = "ms_only")
#'
#' # Use stricter intersection threshold
#' filtered.graph <- construct.gflow.graph(merged.basins, min.intersection = 5)
#' }
#'
#' @seealso
#' \code{\link{merge.clustered.basins}} for creating merged basin structures,
#' \code{\link{cluster.local.extrema}} for clustering basins,
#' \code{\link{compute.basins.of.attraction}} for initial basin computation
#'
#' @export
construct.gflow.graph <- function(merged.basins,
                                  min.intersection = 1,
                                  edge.type = c("all", "ms_only", "same_type")) {
  
  ## Input validation
  if (!inherits(merged.basins, "basins_of_attraction")) {
    stop("merged.basins must be of class 'basins_of_attraction'")
  }
  
  edge.type <- match.arg(edge.type)
  
  if (!is.numeric(min.intersection) || min.intersection < 1) {
    stop("min.intersection must be a positive integer")
  }
  min.intersection <- as.integer(min.intersection)
  
  ## Extract basin structures
  ascending.basins <- merged.basins$lmin_basins
  descending.basins <- merged.basins$lmax_basins
  
  if (length(ascending.basins) == 0 && length(descending.basins) == 0) {
    stop("No basins found in the merged basin structure")
  }
  
  ## Extract vertex sets from basins
  ## For ascending basins (minima)
  ascending.vertices <- list()
  ascending.labels <- character(0)
  ascending.extrema <- integer(0)
  ascending.values <- numeric(0)
  
  if (length(ascending.basins) > 0) {
    for (i in seq_along(ascending.basins)) {
      basin <- ascending.basins[[i]]
      ## Extract vertices from basin_df (first column)
      ascending.vertices[[i]] <- as.integer(basin$basin_df[, 1])
      ascending.extrema[i] <- basin$vertex
      ascending.values[i] <- basin$value
    }
    ## Create labels m1, m2, ... for ascending basins
    ascending.labels <- paste0("m", seq_along(ascending.basins))
    names(ascending.vertices) <- ascending.labels
  }
  
  ## For descending basins (maxima)
  descending.vertices <- list()
  descending.labels <- character(0)
  descending.extrema <- integer(0)
  descending.values <- numeric(0)
  
  if (length(descending.basins) > 0) {
    for (i in seq_along(descending.basins)) {
      basin <- descending.basins[[i]]
      ## Extract vertices from basin_df (first column)
      descending.vertices[[i]] <- as.integer(basin$basin_df[, 1])
      descending.extrema[i] <- basin$vertex
      descending.values[i] <- basin$value
    }
    ## Create labels M1, M2, ... for descending basins
    descending.labels <- paste0("M", seq_along(descending.basins))
    names(descending.vertices) <- descending.labels
  }
  
  ## Combine all basins
  all.vertices <- c(ascending.vertices, descending.vertices)
  all.labels <- c(ascending.labels, descending.labels)
  all.extrema <- c(ascending.extrema, descending.extrema)
  all.values <- c(ascending.values, descending.values)
  
  n.ascending <- length(ascending.basins)
  n.descending <- length(descending.basins)
  n.total <- n.ascending + n.descending
  
  ## Compute basin sizes
  basin.sizes <- sapply(all.vertices, length)
  
  ## Compute intersection matrix
  ## intersection.matrix[i, j] = number of vertices shared by basins i and j
  intersection.matrix <- matrix(0, nrow = n.total, ncol = n.total)
  rownames(intersection.matrix) <- all.labels
  colnames(intersection.matrix) <- all.labels
  
  ## Fill diagonal with basin sizes
  diag(intersection.matrix) <- basin.sizes
  
  ## Compute pairwise intersections
  ## Only compute upper triangle since matrix is symmetric
  for (i in seq_len(n.total - 1)) {
    for (j in (i + 1):n.total) {
      ## Compute set intersection
      intersection.size <- length(intersect(all.vertices[[i]], all.vertices[[j]]))
      intersection.matrix[i, j] <- intersection.size
      intersection.matrix[j, i] <- intersection.size
    }
  }
  
  ## Build adjacency structure based on edge.type
  adjacency.list <- vector("list", n.total)
  weight.list <- vector("list", n.total)
  names(adjacency.list) <- all.labels
  names(weight.list) <- all.labels
  
  for (i in seq_len(n.total)) {
    adjacency.list[[i]] <- integer(0)
    weight.list[[i]] <- numeric(0)
  }
  
  ## Helper function to determine if edge should be included
  should.include.edge <- function(i, j, edge.type, n.ascending) {
    i.is.ascending <- i <= n.ascending
    j.is.ascending <- j <= n.ascending
    
    if (edge.type == "all") {
      return(TRUE)
    } else if (edge.type == "ms_only") {
      ## Only include edges between different types (Morse-Smale edges)
      return(i.is.ascending != j.is.ascending)
    } else if (edge.type == "same_type") {
      ## Only include edges between same types
      return(i.is.ascending == j.is.ascending)
    }
    return(FALSE)
  }
  
  ## Build edges based on intersections and edge.type
  for (i in seq_len(n.total - 1)) {
    for (j in (i + 1):n.total) {
      
      ## Check if edge should be included based on type
      if (!should.include.edge(i, j, edge.type, n.ascending)) {
        next
      }
      
      intersection.size <- intersection.matrix[i, j]
      
      ## Create edge if intersection is large enough
      if (intersection.size >= min.intersection) {
        ## Compute edge weight as Dice similarity coefficient
        ## This normalizes by basin sizes to give a similarity measure
        weight <- (2 * intersection.size) / (basin.sizes[i] + basin.sizes[j])
        
        ## Add edge in both directions (undirected graph)
        adjacency.list[[i]] <- c(adjacency.list[[i]], j)
        weight.list[[i]] <- c(weight.list[[i]], weight)
        
        adjacency.list[[j]] <- c(adjacency.list[[j]], i)
        weight.list[[j]] <- c(weight.list[[j]], weight)
      }
    }
  }
  
  ## Create basin metadata
  basin.metadata <- data.frame(
    index = seq_len(n.total),
    label = all.labels,
    type = c(rep("ascending", n.ascending), rep("descending", n.descending)),
    size = basin.sizes,
    extremum.vertex = all.extrema,
    extremum.value = all.values,
    stringsAsFactors = FALSE
  )
  
  ## Construct result object
  result <- list(
    adjacency.list = adjacency.list,
    weight.list = weight.list,
    intersection.matrix = intersection.matrix,
    basin.metadata = basin.metadata,
    n.ascending = n.ascending,
    n.descending = n.descending,
    edge.type = edge.type,
    min.intersection = min.intersection
  )
  
  class(result) <- "gflow_graph"
  return(result)
}


#' Print Method for Gradient Flow Graph Objects
#'
#' @description
#' Prints a summary of a gradient flow graph object.
#'
#' @param x An object of class \code{"gflow_graph"} as created by
#'   \code{\link{construct.gflow.graph}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}. This function is called
#'   for its side effect of printing a summary to the console.
#'
#' @export
print.gflow_graph <- function(x, ...) {
  cat("Gradient Flow Graph (1-skeleton of nerve complex)\n")
  cat("=================================================\n\n")
  
  cat("Vertices:\n")
  cat(sprintf("  Ascending basins (minima):  %d\n", x$n.ascending))
  cat(sprintf("  Descending basins (maxima): %d\n", x$n.descending))
  cat(sprintf("  Total:                      %d\n\n",
              x$n.ascending + x$n.descending))
  
  ## Count edges (each edge counted once)
  n.edges <- sum(sapply(x$adjacency.list, length)) / 2
  cat("Edges:\n")
  cat(sprintf("  Total edges:                %d\n", n.edges))
  cat(sprintf("  Edge type:                  %s\n", x$edge.type))
  cat(sprintf("  Min intersection threshold: %d\n\n", x$min.intersection))
  
  ## Compute degree distribution
  degrees <- sapply(x$adjacency.list, length)
  cat("Degree distribution:\n")
  cat(sprintf("  Mean degree:                %.2f\n", mean(degrees)))
  cat(sprintf("  Min degree:                 %d\n", min(degrees)))
  cat(sprintf("  Max degree:                 %d\n\n", max(degrees)))
  
  ## Report isolated vertices
  n.isolated <- sum(degrees == 0)
  if (n.isolated > 0) {
    cat(sprintf("  Warning: %d isolated vertices (degree 0)\n\n", n.isolated))
  }
  
  invisible(x)
}


#' Summary Method for Gradient Flow Graph Objects
#'
#' @description
#' Generates a detailed summary of a gradient flow graph, including vertex
#' information, edge statistics, and connectivity properties.
#'
#' @param object An object of class \code{"gflow_graph"} as created by
#'   \code{\link{construct.gflow.graph}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing summary statistics:
#'   \describe{
#'     \item{n.vertices}{Total number of vertices (basins).}
#'     \item{n.ascending}{Number of ascending basins.}
#'     \item{n.descending}{Number of descending basins.}
#'     \item{n.edges}{Total number of edges.}
#'     \item{edge.type}{Type of edges included.}
#'     \item{degree.stats}{Summary statistics of vertex degrees.}
#'     \item{basin.size.stats}{Summary statistics of basin sizes.}
#'     \item{intersection.stats}{Summary statistics of nonzero intersections.}
#'   }
#'
#' @export
summary.gflow_graph <- function(object, ...) {
  degrees <- sapply(object$adjacency.list, length)
  n.edges <- sum(degrees) / 2
  
  ## Get nonzero off-diagonal intersections
  n <- nrow(object$intersection.matrix)
  intersections <- object$intersection.matrix[upper.tri(object$intersection.matrix)]
  nonzero.intersections <- intersections[intersections > 0]
  
  result <- list(
    n.vertices = length(object$adjacency.list),
    n.ascending = object$n.ascending,
    n.descending = object$n.descending,
    n.edges = n.edges,
    edge.type = object$edge.type,
    degree.stats = summary(degrees),
    basin.size.stats = summary(object$basin.metadata$size),
    intersection.stats = if (length(nonzero.intersections) > 0) {
      summary(nonzero.intersections)
    } else {
      "No basin intersections"
    }
  )
  
  class(result) <- "summary.gflow_graph"
  return(result)
}


#' Print Summary for Gradient Flow Graph Objects
#'
#' @description
#' Prints the summary generated by \code{\link{summary.gflow_graph}}.
#'
#' @param x An object of class \code{"summary.gflow_graph"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.summary.gflow_graph <- function(x, ...) {
  cat("Gradient Flow Graph Summary\n")
  cat("===========================\n\n")
  
  cat("Graph Structure:\n")
  cat(sprintf("  Total vertices:             %d\n", x$n.vertices))
  cat(sprintf("  Ascending basins (minima):  %d\n", x$n.ascending))
  cat(sprintf("  Descending basins (maxima): %d\n", x$n.descending))
  cat(sprintf("  Total edges:                %d\n", x$n.edges))
  cat(sprintf("  Edge type:                  %s\n\n", x$edge.type))
  
  cat("Vertex Degrees:\n")
  print(x$degree.stats)
  cat("\n")
  
  cat("Basin Sizes:\n")
  print(x$basin.size.stats)
  cat("\n")
  
  cat("Basin Intersections (nonzero):\n")
  if (is.character(x$intersection.stats)) {
    cat(sprintf("  %s\n", x$intersection.stats))
  } else {
    print(x$intersection.stats)
  }
  cat("\n")
  
  invisible(x)
}

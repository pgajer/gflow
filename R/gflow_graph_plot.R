#' Plot Method for Gradient Flow Graph Objects
#'
#' @description
#' Creates a visualization of a gradient flow graph showing basins as vertices
#' and their intersections as edges. The plot uses different visual encoding
#' to distinguish ascending from descending basins and to represent edge weights.
#'
#' @details
#' The visualization problem for gradient flow graphs differs from standard
#' graph drawing because the vertices carry semantic meaning beyond mere
#' connectivity. We wish to represent not only which basins intersect but
#' also the orientation of gradient flow (ascending versus descending) and
#' the relative importance of different connections.
#'
#' We address this through a layered visual encoding. Vertex position reflects
#' the extremum function value, creating a natural vertical ordering that
#' mirrors the gradient flow structure. Ascending basins (minima) appear toward
#' the bottom of the plot while descending basins (maxima) appear toward the top.
#' This spatial arrangement helps viewers immediately grasp the relationship
#' between basins and the underlying function landscape.
#'
#' Vertex size encodes basin extent, with larger circles representing basins
#' that capture more vertices from the original graph. This allows quick
#' identification of major versus minor topological features. Color distinguishes
#' basin type, with one palette for minima and another for maxima, making the
#' bipartite structure evident when Morse-Smale edges predominate.
#'
#' Edge styling communicates connection strength through both width and
#' transparency. Edges representing larger intersections appear thicker and
#' more opaque, drawing attention to regions where basins substantially overlap.
#' Weak intersections fade into the background, preventing visual clutter while
#' preserving the complete connectivity structure for careful examination.
#'
#' @param x An object of class \code{"gflow_graph"} as created by
#'   \code{\link{construct.gflow.graph}}.
#' @param layout Character string specifying the layout algorithm. Options are:
#'   \describe{
#'     \item{\code{"extremum"}}{Position vertices by extremum value (default).
#'       X-coordinates are random jitter, Y-coordinates reflect function values.}
#'     \item{\code{"force"}}{Use force-directed layout algorithm (requires igraph).}
#'     \item{\code{"circular"}}{Arrange vertices in a circle.}
#'     \item{\code{"bipartite"}}{Arrange ascending and descending basins in two rows.}
#'   }
#' @param vertex.size.scale Numeric scaling factor for vertex sizes (default: 1.0).
#' @param edge.width.scale Numeric scaling factor for edge widths (default: 1.0).
#' @param vertex.color.ascending Color for ascending basin vertices (default: "lightblue").
#' @param vertex.color.descending Color for descending basin vertices (default: "lightcoral").
#' @param vertex.label Logical indicating whether to show vertex labels (default: TRUE).
#' @param vertex.label.cex Numeric character expansion for labels (default: 0.7).
#' @param edge.color Color for edges (default: "gray60").
#' @param edge.transparency Numeric in \eqn{[0,1]} controlling edge transparency,
#'   where 0 is fully transparent and 1 is fully opaque (default: 0.5).
#' @param add.legend Logical indicating whether to add a legend (default: TRUE).
#' @param main Character string for plot title (default: "Gradient Flow Graph").
#' @param ... Additional graphical parameters passed to plot().
#'
#' @return Invisibly returns the coordinates used for vertex placement as a
#'   two-column matrix with rownames matching basin labels.
#'
#' @examples
#' \dontrun{
#' # After constructing a gradient flow graph
#' gflow.graph <- construct.gflow.graph(merged.basins)
#'
#' # Default plot
#' plot(gflow.graph)
#'
#' # Customize appearance
#' plot(gflow.graph,
#'      layout = "bipartite",
#'      vertex.size.scale = 1.5,
#'      vertex.color.ascending = "steelblue",
#'      vertex.color.descending = "firebrick",
#'      main = "Basin Connectivity Structure")
#'
#' # Force-directed layout
#' plot(gflow.graph, layout = "force", edge.transparency = 0.3)
#'
#' # Focus on Morse-Smale edges only
#' ms.graph <- construct.gflow.graph(merged.basins, edge.type = "ms_only")
#' plot(ms.graph, layout = "bipartite")
#' }
#'
#' @seealso
#' \code{\link{construct.gflow.graph}} for creating gradient flow graphs
#'
#' @export
plot.gflow_graph <- function(x,
                             layout = c("extremum", "force", "circular", "bipartite"),
                             vertex.size.scale = 1.0,
                             edge.width.scale = 1.0,
                             vertex.color.ascending = "lightblue",
                             vertex.color.descending = "lightcoral",
                             vertex.label = TRUE,
                             vertex.label.cex = 0.7,
                             edge.color = "gray60",
                             edge.transparency = 0.5,
                             add.legend = TRUE,
                             main = "Gradient Flow Graph",
                             ...) {
  
  layout <- match.arg(layout)
  
  n.total <- length(x$adjacency.list)
  n.asc <- x$n.ascending
  n.desc <- x$n.descending
  
  ## Compute vertex coordinates based on layout
  coords <- compute.layout(x, layout)
  
  ## Set up plot
  plot.new()
  plot.window(xlim = range(coords[, 1]) + c(-0.1, 0.1),
              ylim = range(coords[, 2]) + c(-0.1, 0.1))
  title(main = main, ...)
  
  ## Compute vertex sizes based on basin sizes
  basin.sizes <- x$basin.metadata$size
  vertex.sizes <- sqrt(basin.sizes) * 0.02 * vertex.size.scale
  
  ## Assign vertex colors
  vertex.colors <- rep(vertex.color.ascending, n.total)
  vertex.colors[(n.asc + 1):n.total] <- vertex.color.descending
  
  ## Convert edge transparency to alpha
  edge.alpha <- edge.transparency
  edge.color.with.alpha <- adjustcolor(edge.color, alpha.f = edge.alpha)
  
  ## Draw edges first (so they appear behind vertices)
  for (i in seq_len(n.total)) {
    if (length(x$adjacency.list[[i]]) == 0) next
    
    for (k in seq_along(x$adjacency.list[[i]])) {
      j <- x$adjacency.list[[i]][k]
      
      ## Only draw each edge once (when i < j)
      if (i < j) {
        x0 <- coords[i, 1]
        y0 <- coords[i, 2]
        x1 <- coords[j, 1]
        y1 <- coords[j, 2]
        
        ## Edge width based on weight
        edge.weight <- x$weight.list[[i]][k]
        edge.width <- 0.5 + 3 * edge.weight * edge.width.scale
        
        segments(x0, y0, x1, y1,
                col = edge.color.with.alpha,
                lwd = edge.width)
      }
    }
  }
  
  ## Draw vertices
  for (i in seq_len(n.total)) {
    points(coords[i, 1], coords[i, 2],
           pch = 21,
           bg = vertex.colors[i],
           col = "black",
           cex = vertex.sizes[i] * 2)
  }
  
  ## Add vertex labels
  if (vertex.label) {
    text(coords[, 1], coords[, 2],
         labels = x$basin.metadata$label,
         cex = vertex.label.cex,
         font = 2)
  }
  
  ## Add legend
  if (add.legend) {
    legend("topright",
           legend = c("Ascending (min)", "Descending (max)"),
           pch = 21,
           pt.bg = c(vertex.color.ascending, vertex.color.descending),
           pt.cex = 1.5,
           bty = "n")
  }
  
  invisible(coords)
}


#' Compute Layout Coordinates for Gradient Flow Graph
#'
#' @description
#' Internal function to compute vertex coordinates for different layout algorithms.
#'
#' @param gflow.graph An object of class \code{"gflow_graph"}.
#' @param layout Character string specifying layout type.
#'
#' @return A two-column matrix of coordinates with rownames matching basin labels.
#'
#' @keywords internal
compute.layout <- function(gflow.graph, layout) {
  
  n.total <- length(gflow.graph$adjacency.list)
  n.asc <- gflow.graph$n.ascending
  n.desc <- gflow.graph$n.descending
  
  coords <- matrix(0, nrow = n.total, ncol = 2)
  rownames(coords) <- gflow.graph$basin.metadata$label
  
  if (layout == "extremum") {
    ## Position by extremum value with random jitter
    coords[, 2] <- gflow.graph$basin.metadata$extremum.value
    
    ## Scale y to [0, 1]
    y.range <- range(coords[, 2])
    if (diff(y.range) > 0) {
      coords[, 2] <- (coords[, 2] - y.range[1]) / diff(y.range)
    } else {
      coords[, 2] <- 0.5
    }
    
    ## Random x with slight separation for ascending vs descending
    set.seed(123)
    coords[1:n.asc, 1] <- runif(n.asc, 0, 0.45)
    if (n.desc > 0) {
      coords[(n.asc + 1):n.total, 1] <- runif(n.desc, 0.55, 1)
    }
    
  } else if (layout == "circular") {
    ## Arrange vertices in a circle
    angles <- seq(0, 2 * pi, length.out = n.total + 1)[1:n.total]
    coords[, 1] <- cos(angles)
    coords[, 2] <- sin(angles)
    
  } else if (layout == "bipartite") {
    ## Two rows: ascending on bottom, descending on top
    if (n.asc > 0) {
      coords[1:n.asc, 1] <- seq(0, 1, length.out = n.asc)
      coords[1:n.asc, 2] <- 0
    }
    if (n.desc > 0) {
      coords[(n.asc + 1):n.total, 1] <- seq(0, 1, length.out = n.desc)
      coords[(n.asc + 1):n.total, 2] <- 1
    }
    
  } else if (layout == "force") {
    ## Force-directed layout using igraph
    if (!requireNamespace("igraph", quietly = TRUE)) {
      warning("Package 'igraph' needed for force layout. Using extremum layout instead.")
      return(compute.layout(gflow.graph, "extremum"))
    }
    
    ## Build igraph object
    edge.list <- NULL
    for (i in seq_len(n.total - 1)) {
      if (length(gflow.graph$adjacency.list[[i]]) > 0) {
        for (j in gflow.graph$adjacency.list[[i]]) {
          if (i < j) {
            edge.list <- rbind(edge.list, c(i, j))
          }
        }
      }
    }
    
    if (is.null(edge.list)) {
      ## No edges, fall back to circular
      return(compute.layout(gflow.graph, "circular"))
    }
    
    g <- igraph::graph_from_edgelist(edge.list, directed = FALSE)
    igraph::V(g)$name <- gflow.graph$basin.metadata$label
    
    layout.coords <- igraph::layout_with_fr(g)
    coords <- layout.coords
    rownames(coords) <- gflow.graph$basin.metadata$label
  }
  
  return(coords)
}


#' Export Gradient Flow Graph to igraph Format
#'
#' @description
#' Converts a gradient flow graph object to an igraph graph object for
#' compatibility with igraph's extensive graph analysis and visualization tools.
#'
#' @param gflow.graph An object of class \code{"gflow_graph"} as created by
#'   \code{\link{construct.gflow.graph}}.
#' @param include.vertex.attrs Logical indicating whether to include vertex
#'   attributes (basin type, size, extremum info) in the igraph object
#'   (default: TRUE).
#' @param include.edge.attrs Logical indicating whether to include edge weights
#'   as attributes (default: TRUE).
#'
#' @return An igraph graph object with vertices corresponding to basins and
#'   edges corresponding to basin intersections. If requested, vertex attributes
#'   include: label, type, size, extremum.vertex, extremum.value. Edge attributes
#'   include: weight, intersection.size.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' gflow.graph <- construct.gflow.graph(merged.basins)
#' g <- as_igraph(gflow.graph)
#'
#' # Use igraph functions
#' plot(g, vertex.size = V(g)$size / 10,
#'      vertex.color = ifelse(V(g)$type == "ascending", "blue", "red"))
#'
#' # Compute graph properties
#' betweenness(g)
#' clustering_coefficient(g)
#' }
#'
#' @seealso
#' \code{\link{construct.gflow.graph}}
#'
#' @export
as_igraph <- function(gflow.graph,
                      include.vertex.attrs = TRUE,
                      include.edge.attrs = TRUE) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for this function.")
  }
  
  n.total <- length(gflow.graph$adjacency.list)
  
  ## Build edge list
  edge.list <- NULL
  edge.weights <- NULL
  edge.intersections <- NULL
  
  for (i in seq_len(n.total)) {
    if (length(gflow.graph$adjacency.list[[i]]) == 0) next
    
    for (k in seq_along(gflow.graph$adjacency.list[[i]])) {
      j <- gflow.graph$adjacency.list[[i]][k]
      
      ## Only add each edge once
      if (i < j) {
        edge.list <- rbind(edge.list, c(i, j))
        if (include.edge.attrs) {
          edge.weights <- c(edge.weights, gflow.graph$weight.list[[i]][k])
          edge.intersections <- c(edge.intersections,
                                 gflow.graph$intersection.matrix[i, j])
        }
      }
    }
  }
  
  ## Create igraph object
  if (is.null(edge.list)) {
    ## Graph with no edges
    g <- igraph::make_empty_graph(n = n.total, directed = FALSE)
  } else {
    g <- igraph::graph_from_edgelist(edge.list, directed = FALSE)
  }
  
  ## Add vertex attributes
  if (include.vertex.attrs) {
    igraph::V(g)$name <- gflow.graph$basin.metadata$label
    igraph::V(g)$type <- gflow.graph$basin.metadata$type
    igraph::V(g)$size <- gflow.graph$basin.metadata$size
    igraph::V(g)$extremum.vertex <- gflow.graph$basin.metadata$extremum.vertex
    igraph::V(g)$extremum.value <- gflow.graph$basin.metadata$extremum.value
  } else {
    igraph::V(g)$name <- gflow.graph$basin.metadata$label
  }
  
  ## Add edge attributes
  if (include.edge.attrs && !is.null(edge.weights)) {
    igraph::E(g)$weight <- edge.weights
    igraph::E(g)$intersection.size <- edge.intersections
  }
  
  return(g)
}

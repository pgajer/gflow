#' Plot a Single Component from a Network
#'
#' Extracts and visualizes a single connected component from an igraph object
#' with customizable vertex labels, sizes, and layout.
#'
#' @param x An igraph object
#' @param comp.id Integer specifying which component to plot (by component ID)
#' @param show.labels Logical indicating whether to display vertex labels (default: TRUE)
#' @param label.cex Numeric value controlling label text size (default: 0.7)
#' @param label.dist Numeric value controlling distance of labels from vertices (default: 0.5)
#' @param label.degree Numeric value controlling label placement angle in radians (default: -pi/4)
#' @param vertex.size.scale Numeric multiplier for vertex sizes (default: 2)
#' @param detect.communities Logical indicating whether to detect and color communities (default: TRUE)
#' @param layout.fun Function to compute graph layout (default: layout_with_fr)
#' @param layout.coords Layout coordinates.
#' @param main Character string for plot title (default: auto-generated)
#' @param ... Additional arguments passed to plot.igraph()
#'
#' @return Invisibly returns the induced subgraph of the component
#'
#' @examples
#' \dontrun{
#' ## Plot largest component with labels
#' plot.component(phylo.net, comp.id = 1, show.labels = TRUE, label.cex = 0.6)
#'
#' ## Plot without labels
#' plot.component(phylo.net, comp.id = 1, show.labels = FALSE)
#'
#' ## Plot with custom layout
#' plot.component(phylo.net, comp.id = 1, layout.fun = layout_with_kk)
#' }
#'
#' @export
plot.component <- function(x,
                           comp.id = 1,
                           show.labels = TRUE,
                           label.cex = 0.7,
                           label.dist = 0.5,
                           label.degree = -pi/4,
                           vertex.size.scale = 2,
                           detect.communities = TRUE,
                           layout.fun = igraph::layout_with_fr,
                           layout.coords = NULL,
                           main = NULL,
                           ...) {
    graph <- x

    ## Validate inputs
    if (!inherits(graph, "igraph")) {
        stop("graph must be an igraph object")
    }

    ## Compute connected components
    comp <- igraph::components(graph)

    ## Validate component ID
    if (comp.id < 1 || comp.id > comp$no) {
        stop("comp.id must be between 1 and ", comp$no)
    }

    ## Extract vertices in specified component
    vertices.in.comp <- which(comp$membership == comp.id)

    if (length(vertices.in.comp) == 0) {
        stop("No vertices found in component ", comp.id)
    }

    ## Create induced subgraph
    subgraph <- igraph::induced_subgraph(graph, vertices.in.comp)

    ## Set up vertex labels
    if (show.labels) {
        vertex.labels <- igraph::V(subgraph)$name
    } else {
        vertex.labels <- NA
    }

    ## Detect communities if requested
    if (detect.communities && igraph::vcount(subgraph) > 1) {
        communities <- igraph::cluster_louvain(subgraph)
        vertex.colors <- igraph::membership(communities)
    } else {
        vertex.colors <- "lightblue"
    }

    ## Compute vertex sizes based on degree
    vertex.sizes <- igraph::degree(subgraph) * vertex.size.scale

    ## Ensure minimum size for visibility
    vertex.sizes[vertex.sizes < 5] <- 5

    ## Compute layout
    if (is.null(layout.coords)) {
        layout.coords <- layout.fun(subgraph)
    }


    ## Create default title if not provided
    if (is.null(main)) {
        main <- paste0("Component ", comp.id, " (", comp$csize[comp.id], " nodes)")
    }

    ## Plot the component
    plot(subgraph,
         vertex.color = vertex.colors,
         vertex.size = vertex.sizes,
         vertex.label = vertex.labels,
         vertex.label.cex = label.cex,
         vertex.label.dist = label.dist,
         vertex.label.degree = label.degree,
         vertex.label.color = "black",
         layout = layout.coords,
         main = main,
         ...)

    ## Return subgraph invisibly for further analysis
    invisible(list(subgraph = subgraph,
                   layout.coords = layout.coords,
                   vertex.labels = vertex.labels,
                   communities = communities))
}


#' Plot Multiple Components from a Network
#'
#' Visualizes multiple connected components in a multi-panel layout.
#'
#' @param x An igraph object
#' @param min.size Integer specifying minimum component size to plot (default: 3)
#' @param max.components Integer specifying maximum number of components to plot (default: 12)
#' @param show.labels Logical indicating whether to display vertex labels (default: FALSE)
#' @param label.cex Numeric value controlling label text size (default: 0.6)
#' @param mfrow Numeric vector of length 2 specifying panel layout (default: auto-computed)
#' @param ... Additional arguments passed to plot.component()
#'
#' @return Invisibly returns a vector of component IDs that were plotted
#'
#' @examples
#' \dontrun{
#' ## Plot all components with at least 3 nodes
#' plot.components.multi(phylo.net, min.size = 3)
#'
#' ## Plot with labels and custom layout
#' plot.components.multi(phylo.net, min.size = 3, show.labels = TRUE, mfrow = c(3, 3))
#' }
#'
#' @export
plot.components.multi <- function(x,
                                  min.size = 3,
                                  max.components = 12,
                                  show.labels = FALSE,
                                  label.cex = 0.6,
                                  mfrow = NULL,
                                  ...) {
  graph <- x

  ## Compute connected components
  comp <- igraph::components(graph)

  ## Find components meeting size criterion
  comp.ids <- which(comp$csize >= min.size)

  if (length(comp.ids) == 0) {
    message("No components with size >= ", min.size)
    return(invisible(NULL))
  }

  ## Limit to max.components
  if (length(comp.ids) > max.components) {
    ## Sort by size and take largest
    comp.sizes <- comp$csize[comp.ids]
    comp.ids <- comp.ids[order(-comp.sizes)[1:max.components]]
    message("Showing largest ", max.components, " components out of ",
            length(which(comp$csize >= min.size)))
  }

  ## Compute panel layout if not provided
  if (is.null(mfrow)) {
    n.plots <- length(comp.ids)
    n.cols <- ceiling(sqrt(n.plots))
    n.rows <- ceiling(n.plots / n.cols)
    mfrow <- c(n.rows, n.cols)
  }

  ## Set up multi-panel plot
  old.par <- par(mfrow = mfrow, mar = c(2, 2, 3, 1))
  on.exit(par(old.par))

  ## Plot each component
  for (comp.id in comp.ids) {
    plot.component(graph,
                  comp.id = comp.id,
                  show.labels = show.labels,
                  label.cex = label.cex,
                  ...)
  }

  invisible(comp.ids)
}

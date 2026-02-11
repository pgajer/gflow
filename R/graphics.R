.ensure_custom_shapes_registered <- function() {
    if (!("triangle" %in% names(igraph::shapes()))) {
        igraph::add_shape("triangle",
                          clip = igraph::shapes(shape = "circle")$clip,
                          plot = triangle.plot)
    }
    if (!("itriangle" %in% names(igraph::shapes()))) {
        igraph::add_shape("itriangle",
                          clip = igraph::shapes(shape = "circle")$clip,
                          plot = itriangle.plot)
    }
}

#' Create a kh.matrix object
#'
#' @description
#' Constructor function for creating a kh.matrix object that stores kernel and bandwidth
#' information for model selection or cross-validation procedures.
#'
#' @param kh.mat Numeric matrix containing kernel values or cross-validation results
#' @param existing.k Integer vector specifying which k values were evaluated
#' @param h.values Numeric vector of bandwidth (h) values tested
#' @param id Character string or identifier for this kh.matrix object
#'
#' @return An object of class "kh.matrix" containing the input components
#'
#' @export
kh.matrix <- function(kh.mat, existing.k, h.values, id) {
  structure(
    list(
      kh.mat = kh.mat,
      existing.k = existing.k,
      h.values = h.values,
      id = id
    ),
    class = "kh.matrix"
  )
}

#' Create a prediction.errors object
#'
#' @description
#' Constructor function for creating a prediction.errors object that stores model
#' prediction errors along with optional x-values for plotting or analysis.
#'
#' @param errors Numeric vector of prediction errors
#' @param xvals Optional numeric vector of x-values corresponding to the errors.
#'        Default is NULL.
#'
#' @return An object of class "prediction.errors" containing the errors with
#'         xvals stored as an attribute
#'
#' @export
prediction.errors <- function(errors, xvals = NULL) {
  structure(
    errors,
    xvals = xvals,
    class = "prediction.errors"
  )
}

#' Create a chain.with.path object
#'
#' @description
#' Constructor function for creating a chain.with.path object that represents
#' a graph structure with adjacency and weight information, typically used in
#' path-based algorithms or network analysis.
#'
#' @param adj.list List representing the adjacency structure of the graph
#' @param weight.list List of weights corresponding to edges in the adjacency list
#' @param gpd.obj Object containing additional graph path data or metadata
#'
#' @return An object of class "chain.with.path" containing the graph structure
#'         and associated data
#'
#' @export
chain.with.path <- function(adj.list, weight.list, gpd.obj) {
  structure(
    list(
      adj.list = adj.list,
      weight.list = weight.list,
      gpd.obj = gpd.obj
    ),
    class = "chain.with.path"
  )
}

#' Create a model.errors object
#'
#' @description
#' Constructor function for creating a model.errors object that stores integrated
#' error measures and their bootstrap distributions for model evaluation.
#'
#' @param integrals Numeric vector or matrix of integrated error measures
#' @param bb.integrals Numeric matrix containing bootstrap samples of the integrated
#'        error measures, used for uncertainty quantification
#'
#' @return An object of class "model.errors" containing the error integrals and
#'         their bootstrap distributions
#'
#' @export
model.errors <- function(integrals, bb.integrals) {
  structure(
    list(
      integrals = integrals,
      bb.integrals = bb.integrals
    ),
    class = "model.errors"
  )
}

#' Create a gflow graph object
#'
#' @param adj.list List where element i contains indices of vertices adjacent to vertex i
#' @param weight.list List where element i contains weights of edges from vertex i
#' @return An object of class "ggraph"
#' @export
ggraph <- function(adj.list, weight.list = NULL) {
  # Validate inputs
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }

  # If weight.list is not provided, create unit weights
  if (is.null(weight.list)) {
    weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
  }

  # Check that lengths match
  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have the same length")
  }

  # Check that each element has matching lengths
  for (i in seq_along(adj.list)) {
    if (length(adj.list[[i]]) != length(weight.list[[i]])) {
      stop(sprintf("adj.list[[%d]] and weight.list[[%d]] must have the same length", i, i))
    }
  }

  structure(
    list(
      adj.list = adj.list,
      weight.list = weight.list
    ),
    class = "ggraph"
  )
}

#' Create a graph.3d object
#'
#' @param plot.result A list returned by plot.graph()/plot.ggraph(), containing
#'        the graph and layout elements
#' @param z A numeric vector containing z values for each vertex in the graph
#' @return An object of class "graph.3d"
#' @export
create.graph.3d.obj <- function(plot.result, z) {
  # Validate inputs
  if (!is.list(plot.result) || is.null(plot.result$graph) || is.null(plot.result$layout)) {
    stop("plot.result must be a list containing 'graph' and 'layout' elements, as returned by plot.graph()")
  }

  if (!is.numeric(z)) {
    stop("z must be a numeric vector")
  }

  if (length(z) != nrow(plot.result$layout)) {
    stop("Length of z must match the number of vertices in the graph")
  }

  structure(
    list(
      graph = plot.result$graph,
      layout = plot.result$layout,
      z = z,
      vertex.color = plot.result$vertex.color  # Preserve vertex colors if available
    ),
    class = "graph.3d"
  )
}

#' triangle.plot
#'
#' Custom plot function to draw vertices as triangles in an igraph plot
#'
#' This function is used to plot vertices as equilateral triangles in an igraph
#' graph. The triangles are sized to have approximately the same area as circles
#' of the same vertex size. It is intended to be used with the igraph `add_shape`
#' function to add the "triangle" shape to graph vertices.
#'
#' @param coords A matrix of vertex coordinates. Each row should contain the x
#'        and y coordinates of a vertex.
#' @param v An optional vector of vertex indices to be plotted. If NULL, all
#'        vertices will be plotted.
#' @param params A list of plotting parameters, typically obtained from the
#'        igraph `params` function.
#'
#' @details The function adjusts the size of the triangles to match the area of
#'          circles with the same vertex size. This ensures visual consistency
#'          in the plot. The function also handles vertex colors, frame colors,
#'          and frame widths.
#'
#' @examples
#' \dontrun{
#'   library(igraph)
#'
#'   # Define the triangle plot function
#'   triangle.plot <- function(coords, v = NULL, params) {
#'       vertex.color <- params("vertex", "color")
#'       if (length(vertex.color) != 1 && !is.null(v)) {
#'           vertex.color <- vertex.color\[v\]
#'       }
#'       vertex.frame.color <- params("vertex", "frame.color")
#'       if (length(vertex.frame.color) != 1 && !is.null(v)) {
#'           vertex.frame.color <- vertex.frame.color\[v\]
#'       }
#'       vertex.frame.width <- params("vertex", "frame.width")
#'       if (length(vertex.frame.width) != 1 && !is.null(v)) {
#'           vertex.frame.width <- vertex.frame.width\[v\]
#'       }
#'       vertex.size <- 1/200 * params("vertex", "size")
#'       if (length(vertex.size) != 1 && !is.null(v)) {
#'           vertex.size <- vertex.size\[v\]
#'       }
#'       vertex.size <- rep(vertex.size, length.out = nrow(coords))
#'
#'       # Adjust the size for the triangle
#'       side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2
#'
#'       vertex.frame.color\[vertex.frame.width <= 0\] <- NA
#'       vertex.frame.width\[vertex.frame.width <= 0\] <- 1
#'
#'       for (i in 1:nrow(coords)) {
#'           x <- coords\[i, 1\]
#'           y <- coords\[i, 2\]
#'           size <- side.length\[i\]
#'           polygon(x + size * c(cos(pi/2), cos(7*pi/6), cos(11*pi/6)),
#'                   y + size * c(sin(pi/2), sin(7*pi/6), sin(11*pi/6)),
#'                   col = vertex.color\[i\], border = vertex.frame.color\[i\],
#'                   lwd = vertex.frame.width\[i\])
#'       }
#'   }
#'
#'   # Add the triangle shape to igraph
#'   add_shape("triangle", clip = shapes(shape = "circle")$clip, plot = triangle.plot)
#'
#'   # Example graph
#'   g <- make_ring(10)
#'   V(g)$shape <- rep(c("circle", "triangle"), length.out = vcount(g))
#'   plot(g, vertex.size = 15, vertex.color = "skyblue")
#' }
#'
#' @importFrom graphics polygon
#' @export
triangle.plot <- function(coords, v = NULL, params) {

    .ensure_custom_shapes_registered()

    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    vertex.size <- rep(vertex.size, length.out = nrow(coords))

    ## Adjust the size for the triangle
    side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2

    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1

    for (i in 1:nrow(coords)) {
        x <- coords[i, 1]
        y <- coords[i, 2]
        size <- side.length[i]
        polygon(x + size * c(cos(pi/2), cos(7*pi/6), cos(11*pi/6)),
                y + size * c(sin(pi/2), sin(7*pi/6), sin(11*pi/6)),
                col = vertex.color[i], border = vertex.frame.color[i], lwd = vertex.frame.width[i])
    }
}

#' itriangle.plot
#'
#' Custom plot function to draw vertices as inverted triangles in an igraph plot
#'
#' This function is used to plot vertices as inverted equilateral triangles in an igraph
#' graph. The triangles are sized to have approximately the same area as circles
#' of the same vertex size. It is intended to be used with the igraph `add_shape`
#' function to add the "itriangle" (inverted triangle) shape to graph vertices.
#'
#' @param coords A matrix of vertex coordinates. Each row should contain the x
#'        and y coordinates of a vertex.
#' @param v An optional vector of vertex indices to be plotted. If NULL, all
#'        vertices will be plotted.
#' @param params A list of plotting parameters, typically obtained from the
#'        igraph `params` function.
#'
#' @details The function adjusts the size of the triangles to match the area of
#'          circles with the same vertex size. This ensures visual consistency
#'          in the plot. The function also handles vertex colors, frame colors,
#'          and frame widths.
#'
#' @examples
#' \dontrun{
#'   library(igraph)
#'
#'   # Define the inverted triangle plot function
#'   itriangle.plot <- function(coords, v = NULL, params) {
#'       vertex.color <- params("vertex", "color")
#'       if (length(vertex.color) != 1 && !is.null(v)) {
#'           vertex.color <- vertex.color\[v\]
#'       }
#'       vertex.frame.color <- params("vertex", "frame.color")
#'       if (length(vertex.frame.color) != 1 && !is.null(v)) {
#'           vertex.frame.color <- vertex.frame.color\[v\]
#'       }
#'       vertex.frame.width <- params("vertex", "frame.width")
#'       if (length(vertex.frame.width) != 1 && !is.null(v)) {
#'           vertex.frame.width <- vertex.frame.width\[v\]
#'       }
#'       vertex.size <- 1/200 * params("vertex", "size")
#'       if (length(vertex.size) != 1 && !is.null(v)) {
#'           vertex.size <- vertex.size\[v\]
#'       }
#'       vertex.size <- rep(vertex.size, length.out = nrow(coords))
#'
#'       # Adjust the size for the inverted triangle
#'       side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2
#'
#'       vertex.frame.color\[vertex.frame.width <= 0\] <- NA
#'       vertex.frame.width\[vertex.frame.width <= 0\] <- 1
#'
#'       for (i in 1:nrow(coords)) {
#'           x <- coords\[i, 1\]
#'           y <- coords\[i, 2\]
#'           size <- side.length\[i\]
#'           polygon(x + size * c(cos(3*pi/2), cos(5*pi/6), cos(pi/6)),
#'                   y + size * c(sin(3*pi/2), sin(5*pi/6), sin(pi/6)),
#'                   col = vertex.color\[i\],
#'                   border = vertex.frame.color\[i\],
#'                   lwd = vertex.frame.width\[i\])
#'       }
#'   }
#'
#'   # Add the inverted triangle shape to igraph
#'   add_shape("itriangle", clip = shapes(shape = "circle")$clip,
#'             plot = itriangle.plot)
#'
#'   # Example graph
#'   g <- make_ring(10)
#'   V(g)$shape <- rep(c("circle", "itriangle"), length.out = vcount(g))
#'   plot(g, vertex.size = 15, vertex.color = "skyblue")
#' }
#'
#' @importFrom graphics polygon
#' @export
itriangle.plot <- function(coords, v = NULL, params) {

    .ensure_custom_shapes_registered()

    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    vertex.size <- rep(vertex.size, length.out = nrow(coords))

                                        # Adjust the size for the inverted triangle
    side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2

    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1

    for (i in 1:nrow(coords)) {
        x <- coords[i, 1]
        y <- coords[i, 2]
        size <- side.length[i]
        polygon(x + size * c(cos(3*pi/2), cos(5*pi/6), cos(pi/6)),
                y + size * c(sin(3*pi/2), sin(5*pi/6), sin(pi/6)),
                col = vertex.color[i], border = vertex.frame.color[i],
                lwd = vertex.frame.width[i])
    }
}

#' Plot a ggraph Object
#'
#' Creates a visualization of a ggraph object in either 2D or 3D.
#'
#' @param x A ggraph object
#' @param y Optional numeric vector over vertices (used for color coding if provided)
#' @param dim Embedding dimension: 2 or 3 (default 2)
#' @param vertex.size Numeric scalar for 2D point size (igraph)
#' @param vertex.radius Numeric sphere radius multiplier for 3D (combined with vertex.size)
#' @param vertex.label Character vector of vertex labels
#' @param layout Layout algorithm name (e.g., "kk", "fr", "in_circle", "lgl", ...)
#' @param vertex.color Optional vector of colors for vertices (overrides y if provided)
#' @param draw.edges Logical; draw edges in 3D
#' @param legend.cex Legend size when y is not NULL
#' @param legend.position Legend position for 2D (e.g., "topleft")
#' @param quantize.method Method to quantize y: "uniform" or "quantile"
#' @param label.cex Size of vertex labels in 3D
#' @param label.adj Horizontal/vertical adjustment of labels (3D)
#' @param use.saved.layout Optional matrix of coordinates to reuse
#' @param ... Additional parameters passed to igraph plotting in 2D
#'
#' @return Invisibly returns a list: graph, layout, vertex.color, and y-color info if used
#'
#' @examples
#' \dontrun{
#' adj_list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight_list <- list(c(1, 1), c(1, 1), c(1, 1))
#' g <- ggraph(adj_list, weight_list)
#' plot.res <- plot(g)
#' plot(g, vertex.size = 5, vertex.label = c("A", "B", "C"))
#' plot(g, layout = "in_circle")
#' plot(g, vertex.size = 2, use.saved.layout = plot.res$layout)
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'   plot(g, dim = 3)
#' }
#' }
#' @method plot ggraph
#' @export
plot.ggraph <- function(x,
                        y = NULL,
                        dim = 2,
                        vertex.size = 1,
                        vertex.radius = 0.1,
                        vertex.label = NA,
                        layout = "kk",
                        vertex.color = NULL,
                        draw.edges = TRUE,
                        legend.cex = 1,
                        legend.position = "topleft",
                        quantize.method = "uniform",
                        label.cex = 1.2,
                        label.adj = c(0.5, 0.5),
                        use.saved.layout = NULL,
                        ...) {

  if (!dim %in% c(2, 3)) stop("'dim' must be 2 or 3")
  dots <- list(...)
  # Backward compatibility for historical misspelling.
  if (is.null(vertex.color) && !is.null(dots$vertex.colory)) {
      vertex.color <- dots$vertex.colory
  }
  dots$vertex.colory <- NULL

    ## Extract from ggraph object
    graph.adj.list    <- x$adj.list
    graph.edge.lengths <- x$weight.list

    graph.obj <- convert.adjacency.to.edge.matrix(graph.adj.list, graph.edge.lengths)
    g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
    if (!is.null(graph.edge.lengths)) igraph::E(g)$weight <- graph.obj$weights

    ## Default vertex colors by components (stable palette)
    if (is.null(vertex.color)) {
        comps <- igraph::components(g)
        pal <- grDevices::hcl.colors(comps$no, palette = "Set3")
        vertex.color <- pal[comps$membership]
    }

    .ensure_custom_shapes_registered()

    ## Resolve layout function
    if (is.character(layout)) {
        layout_func <- switch(layout,
                              "nicely"    = igraph::layout_nicely,
                              "randomly"  = igraph::layout_randomly,
                              "in_circle" = igraph::layout_in_circle,
                              "on_sphere" = igraph::layout_on_sphere,
                              "fr"        = igraph::layout_with_fr,
                              "kk"        = igraph::layout_with_kk,
                              "as_tree"   = igraph::layout_as_tree,
                              "lgl"       = igraph::layout_with_lgl,
                              "graphopt"  = igraph::layout_with_graphopt,
                              "drl"       = igraph::layout_with_drl,
                              "dh"        = igraph::layout_with_dh,
                              "mds"       = igraph::layout_with_mds,
                              stop("Invalid layout specified")
                              )
    } else {
        layout_func <- layout
    }

    ## Coordinates
    if (!is.null(use.saved.layout)) {
        layout_coords <- use.saved.layout
    } else {
        ## many igraph layouts accept dim=2 or 3; ignored by some
        layout_coords <- layout_func(g, dim = dim)
    }
    layout_coords <- as.matrix(layout_coords)
    if (dim == 3 && ncol(layout_coords) < 3) {
        ## be defensive: pad with zeros if a 2D layout was returned
        layout_coords <- cbind(layout_coords, z = 0)
    }

    ## Color via y if supplied (both 2D and 3D)
    y_info <- NULL
    if (!is.null(y)) {
        q <- quantize.cont.var(y, method = quantize.method)
        y.cat <- q$x.cat
        y.col.tbl <- q$x.col.tbl
        vertex.color <- y.col.tbl[y.cat]
        y_info <- list(y.color = vertex.color, y.col.tbl = y.col.tbl)
    }

    ## Return info
    result <- c(list(graph = g, layout = layout_coords, vertex.color = vertex.color), y_info)

    ## 2D branch: uses base graphics; safe on headless (R uses an off-screen device)
    if (dim == 2) {
        do.call(
            igraph::plot.igraph,
            c(
                list(
                    g,
                    layout = layout_coords,
                    vertex.color = vertex.color,
                    vertex.size = vertex.size,
                    vertex.label = vertex.label
                ),
                dots
            )
        )
        if (!is.null(y_info)) {
            graphics::legend(legend.position, xpd = NA,
                             legend = names(y_info$y.col.tbl),
                             fill   = y_info$y.col.tbl,
                             inset  = 0.01, cex = legend.cex
                             )
        }
        return(invisible(result))
    }

    ## 3D branch: rgl-gated & headless-safe
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("This function requires the optional package 'rgl' for 3D plots. ",
             "Install with install.packages('rgl').", call. = FALSE)
    }

    use_null <- (!interactive()) ||
        identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
        (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
    old_opt <- options(rgl.useNULL = use_null)
    on.exit(options(old_opt), add = TRUE)

    rgl::open3d()
    if (use_null) {
        ## Only close the device automatically if using null device
        on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
    }
    rgl::clear3d()

    ## Points (spheres); radius scales by vertex.size
    rgl::plot3d(layout_coords, axes = FALSE, xlab = "", ylab = "", zlab = "")
    rgl::spheres3d(layout_coords, radius = vertex.radius * vertex.size, col = vertex.color)

    ## Labels
    if (!all(is.na(vertex.label))) {
        text_offset <- 1.3 * (vertex.radius * vertex.size)
        label_idx <- which(!is.na(vertex.label))
        if (length(label_idx)) {
            for (i in label_idx) {
                rgl::text3d(layout_coords[i, 1], layout_coords[i, 2], layout_coords[i, 3] + text_offset,
                            texts = vertex.label[i], color = "black", cex = label.cex, adj = label.adj)
            }
        }
    }

    ## Edges
    if (isTRUE(draw.edges)) {
        edges <- igraph::as_edgelist(g, names = FALSE)
        if (NROW(edges)) {
            for (k in seq_len(NROW(edges))) {
                v1 <- edges[k, 1]; v2 <- edges[k, 2]
                rgl::lines3d(rbind(layout_coords[v1, ], layout_coords[v2, ]),
                             col = "gray50", lwd = 1)
            }
        }
    }

    invisible(result)
}


#' Plot Method for graphMScx Objects
#'
#' @description
#' Visualizes various aspects of a Morse-Smale complex computed on a graph, including
#' the complex itself, cell graphs, branches of attraction, and the original graph
#' with cell assignments. The graphMScx object is typically created by the
#' \code{graph.MS.cx} function.
#'
#' @param x A graphMScx object containing the results of Morse-Smale complex
#'   construction, as returned by \code{graph.MS.cx}
#' @param type Character string specifying the type of plot to generate:
#'   \describe{
#'     \item{\code{"MS_cx_graph"}}{Plot the Morse-Smale complex graph with vertices
#'       colored by type (local minimum/maximum) using a bipartite layout}
#'     \item{\code{"MS_cx_nerve_graph"}}{Plot the Cech (nerve) graph of the
#'       Morse-Smale complex with edges weighted by the number of shared vertices
#'       between cells}
#'     \item{\code{"lmin_BoA"}}{Plot the graph with vertices colored by their
#'       local minima branches of attraction}
#'     \item{\code{"lmax_BoA"}}{Plot the graph with vertices colored by their
#'       local maxima branches of attraction}
#'     \item{\code{"MS_graphs"}}{Plot a specific Morse-Smale cell subgraph
#'       (specify cell with parameter \code{i})}
#'     \item{\code{"MS_cell"}}{Plot the entire graph highlighting a specific cell
#'       with vertices colored by function values}
#'     \item{\code{"Ey"}}{Plot the entire graph with vertices colored by function
#'       values (Ey)}
#'   }
#' @param i Integer specifying the index of the Morse-Smale cell to plot when
#'   \code{type = "MS_graphs"} or \code{type = "MS_cell"}. Default is 1.
#' @param vertex.size Numeric scalar or vector specifying vertex sizes. For local
#'   extrema, this is scaled by \code{lext.sfactor}. Default is 2.
#' @param vertex.label Character vector of vertex labels. If NULL (default),
#'   labels are determined based on \code{show.only.lextr.labels}.
#' @param show.only.lextr.labels Logical. If TRUE (default), only show labels
#'   for local extrema vertices.
#' @param edge.label.cex Numeric scalar controlling the size of edge labels in
#'   the Morse-Smale complex graph. Default is 2.
#' @param legend.cex Numeric scalar controlling legend text size. Default is 1.
#' @param lmin.color Color for local minima vertices. Default is "blue".
#' @param lmax.color Color for local maxima vertices. Default is "red".
#' @param label.dist Numeric scalar controlling the distance of vertex labels
#'   from vertices. Default is 1.5.
#' @param legend.position Character string specifying legend position. Options
#'   include "topleft", "topright", "bottomleft", "bottomright". Default is "topleft".
#' @param lext.sfactor Numeric scalar for scaling the size of local extrema
#'   vertices relative to regular vertices. Default is 5.
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisible NULL. The function is called for its side effect of creating a plot.
#'
#' @details
#' The function provides multiple visualization options for understanding the
#' structure of the Morse-Smale complex:
#' \itemize{
#'   \item The complex graph shows connections between critical points
#'   \item Branch visualizations show which extrema influence each vertex
#'   \item Cell visualizations show the decomposition of the graph
#' }
#'
#' Local minima are displayed as downward-pointing triangles (blue by default),
#' while local maxima are shown as upward-pointing triangles (red by default).
#'
#' @note
#' This function requires the igraph package for graph visualization. The input
#' object must be properly formatted with all required components from the
#' Morse-Smale complex computation.
#'
#' @examples
#' \dontrun{
#' # Create a simple graph
#' adj_list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
#' weight_list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
#' g <- ggraph(adj_list, weight_list)
#'
#' # Define function values on vertices
#' Ey <- c(0.5, 0.3, 0.8, 0.2)
#'
#' # Compute Morse-Smale complex
#' mscx <- graph.MS.cx(g, Ey)
#'
#' # Plot the Morse-Smale complex graph
#' plot(mscx)
#'
#' # Show branches of attraction for local minima
#' plot(mscx, type = "lmin_BoA")
#'
#' # Plot a specific cell
#' plot(mscx, type = "MS_graphs", i = 2)
#'
#' # Show the entire graph with first cell highlighted
#' plot(mscx, type = "MS_cell", i = 1)
#' }
#'
#' @seealso
#' \code{\link{graph.MS.cx}} for computing the Morse-Smale complex,
#' \code{\link{ggraph}} for creating graph objects
#'
#' @importFrom igraph graph_from_edgelist E V layout_with_kk layout_as_bipartite vcount graph_from_adjacency_matrix
#' @export
plot.graphMScx <- function(x,
                           type = c("MS_cx_graph", "MS_cx_nerve_graph", "lmin_BoA",
                                   "lmax_BoA", "MS_graphs", "MS_cell", "Ey"),
                           i = 1,
                           vertex.size = 2,
                           vertex.label = NULL,
                           show.only.lextr.labels = TRUE,
                           edge.label.cex = 2,
                           legend.cex = 1,
                           lmin.color = "blue",
                           lmax.color = "red",
                           label.dist = 1.5,
                           legend.position = "topleft",
                           lext.sfactor = 5,
                           ...) {

    # Validate input
    if (!inherits(x, "graphMScx")) {
        stop("x must be a graphMScx object")
    }

    # Match argument
    type <- match.arg(type)

    # Use consistent naming
    res <- x

    if (type == "MS_cx_graph") {

        graph.m <- convert.adjacency.to.edge.matrix(res$MS_graph, res$MS_graph_weights)
        g <- igraph::graph_from_edgelist(graph.m$edge.matrix, directed = TRUE)

        if (!is.null(res$MS_graph_weights)) {
            igraph::E(g)$weight <- graph.m$weights
        }

        igraph::V(g)$type <- res$MS_graph_types

        if (length(res$MS_graph) == 1) {
            layout <- igraph::layout_with_kk
        } else {
            layout <- igraph::layout_as_bipartite(g, igraph::V(g)$type == res$MS_graph_types[1])
        }

        if (!is.null(res$MS_graph_edge_label)) {
            edge.labels <- convert.edge.label.list.to.edge.label.vector(res$MS_graph_edge_label, rm.duplicates = FALSE)
            igraph::E(g)$label <- edge.labels
        }

        # Identify upper and lower nodes
        upper_nodes <- which(igraph::V(g)$type == 1)
        lower_nodes <- which(igraph::V(g)$type == 2)

        # Create vectors for label degree and distance
        label_degree <- rep(0, igraph::vcount(g))
        label_degree[upper_nodes] <- pi/2    # 90 degrees (above)
        label_degree[lower_nodes] <- 3*pi/2  # 270 degrees (below)

        label_dist_vec <- rep(label.dist, igraph::vcount(g))

        vertex.color <- ifelse(res$MS_graph_types == 1, lmin.color, lmax.color)

        plot(g,
             layout = layout,
             vertex.label.degree = label_degree,
             vertex.label.dist = label_dist_vec,
             vertex.color = vertex.color,
             edge.label.cex = edge.label.cex,
             vertex.shape = ifelse(igraph::V(g)$type == res$MS_graph_types[1], "csquare", "square"),
             vertex.size = vertex.size,
             vertex.label = res$MS_graph_labels,
             ...)

    } else if (type == "MS_cx_nerve_graph") {

        # Creating a nerve graph from the covering of the graph vertices by trajectory cells
        nerve.res <- nerve.graph(res$MS_cell_cc_vertices)

        # Extracting the corresponding weights matrix
        nerve.weights <- convert.adjacency.list.to.adjacency.matrix(
            nerve.res$adjacency.list,
            nerve.res$weights.list
        )

        # Creating an undirected weighted graph
        g <- igraph::graph_from_adjacency_matrix(
            nerve.weights,
            mode = "undirected",
            weighted = TRUE,
            diag = FALSE
        )

        plot(g, edge.label = igraph::E(g)$weight, vertex.shape = "square", ...)

    } else if (type == "lmin_BoA") {

        # Get extrema
        lmin <- unique(res$MS_cx[,"local_min"])
        lmax <- unique(res$MS_cx[,"local_max"])
        n <- length(res$Ey)

        # Set vertex shapes
        vertex.shape <- rep("circle", n)
        vertex.shape[lmin] <- "csquare"
        vertex.shape[lmax] <- "square"

        # Set vertex sizes
        vertex.size.vec <- rep(vertex.size, n)
        vertex.size.vec[c(lmin, lmax)] <- lext.sfactor * vertex.size

        # Color by local minimum branch
        vertex.lmin <- res$MS_cx[,"local_min"]
        color.tbl <- rainbow(length(unique(vertex.lmin)))
        names(color.tbl) <- as.character(sort(unique(vertex.lmin)))
        vertex.color <- color.tbl[as.character(vertex.lmin)]

        # Handle labels
        if (show.only.lextr.labels && is.null(vertex.label)) {
            lextr <- c(lmin, lmax)
            vertex.label <- rep(NA, n)
            vertex.label[lextr] <- lextr
        }

        # Create and plot ggraph
        g <- ggraph(res$graph$adj.list, res$graph$weight.list)
        plot(g,
             vertex.color = vertex.color,
             vertex.size = vertex.size.vec,
             vertex.label = vertex.label,
             vertex.shape = vertex.shape,
             ...)

        legend(legend.position,
               xpd = NA,
               legend = names(color.tbl),
               fill = color.tbl,
               inset = 0.01,
               cex = legend.cex,
               title = "lmin BoA")

    } else if (type == "lmax_BoA") {

        # Get extrema
        lmin <- unique(res$MS_cx[,"local_min"])
        lmax <- unique(res$MS_cx[,"local_max"])
        n <- length(res$Ey)

        # Set vertex shapes
        vertex.shape <- rep("circle", n)
        vertex.shape[lmin] <- "csquare"
        vertex.shape[lmax] <- "square"

        # Set vertex sizes
        vertex.size.vec <- rep(vertex.size, n)
        vertex.size.vec[c(lmin, lmax)] <- lext.sfactor * vertex.size

        # Color by local maximum branch
        vertex.lmax <- res$MS_cx[,"local_max"]
        color.tbl <- rainbow(length(unique(vertex.lmax)))
        names(color.tbl) <- as.character(sort(unique(vertex.lmax)))
        vertex.color <- color.tbl[as.character(vertex.lmax)]

        # Handle labels
        if (show.only.lextr.labels && is.null(vertex.label)) {
            lextr <- c(lmin, lmax)
            vertex.label <- rep(NA, n)
            vertex.label[lextr] <- lextr
        }

        # Create and plot ggraph
        g <- ggraph(res$graph$adj.list, res$graph$weight.list)
        plot(g,
             vertex.color = vertex.color,
             vertex.size = vertex.size.vec,
             vertex.label = vertex.label,
             vertex.shape = vertex.shape,
             ...)

        legend(legend.position,
               xpd = NA,
               legend = names(color.tbl),
               fill = color.tbl,
               inset = 0.01,
               cex = legend.cex,
               title = "lmax BoA")

    } else if (type == "MS_graphs") {

        # Validate cell index
        if (i < 1 || i > length(res$MS_cell_cc_vertices)) {
            stop(sprintf("Invalid cell index i=%d. Must be between 1 and %d",
                        i, length(res$MS_cell_cc_vertices)))
        }

        # Get cell information
        cell.vertices <- res$MS_cell_cc_vertices[[i]]
        lmin <- unique(res$MS_cx[,"local_min"])
        lmax <- unique(res$MS_cx[,"local_max"])

        # Find extrema in this cell
        lmin.in.cell <- intersect(cell.vertices, lmin)
        lmax.in.cell <- intersect(cell.vertices, lmax)

        # Map to local indices
        lmin.local.idx <- which(cell.vertices %in% lmin.in.cell)
        lmax.local.idx <- which(cell.vertices %in% lmax.in.cell)

        n <- length(cell.vertices)

        # Set vertex shapes
        vertex.shape <- rep("circle", n)
        vertex.shape[lmin.local.idx] <- "csquare"
        vertex.shape[lmax.local.idx] <- "square"

        # Set vertex sizes
        vertex.size.vec <- rep(vertex.size, n)
        vertex.size.vec[c(lmin.local.idx, lmax.local.idx)] <- lext.sfactor * vertex.size

        # Color by function values
        q <- quantize.cont.var(res$Ey[cell.vertices])
        vertex.color <- q$x.col.tbl[q$x.cat]

        # Handle labels
        if (show.only.lextr.labels && is.null(vertex.label)) {
            lextr.local <- c(lmin.local.idx, lmax.local.idx)
            vertex.label <- rep(NA, n)
            vertex.label[lextr.local] <- cell.vertices[lextr.local]
        } else if (is.null(vertex.label)) {
            vertex.label <- cell.vertices
        }

        # Plot the cell subgraph
        plot(res$MS_cell_graphs[[i]],
             vertex.size = vertex.size.vec,
             vertex.label = vertex.label,
             vertex.shape = vertex.shape,
             vertex.color = vertex.color,
             ...)

        legend(legend.position,
               xpd = NA,
               legend = names(q$x.col.tbl),
               fill = q$x.col.tbl,
               inset = 0.01,
               cex = legend.cex,
               title = sprintf("Cell %d", i))

    } else if (type == "MS_cell") {

        # Validate cell index
        if (i < 1 || i > length(res$MS_cell_cc_vertices)) {
            stop(sprintf("Invalid cell index i=%d. Must be between 1 and %d",
                        i, length(res$MS_cell_cc_vertices)))
        }

        cell.vertices <- res$MS_cell_cc_vertices[[i]]
        n.all <- length(res$Ey)

        # Get extrema
        lmin <- unique(res$MS_cx[,"local_min"])
        lmax <- unique(res$MS_cx[,"local_max"])

        # Initialize shapes for all vertices
        all.vertex.shape <- rep("circle", n.all)
        all.vertex.shape[lmin] <- "csquare"
        all.vertex.shape[lmax] <- "square"

        # Initialize sizes for all vertices
        all.vertex.size <- rep(vertex.size, n.all)
        all.vertex.size[c(lmin, lmax)] <- lext.sfactor * vertex.size

        # Color only the cell vertices
        q <- quantize.cont.var(res$Ey)
        all.vertex.color <- rep(NA, n.all)
        all.vertex.color[cell.vertices] <- q$x.col.tbl[q$x.cat[cell.vertices]]

        # Handle labels
        if (show.only.lextr.labels && is.null(vertex.label)) {
            lextr <- c(lmin, lmax)
            vertex.label <- rep(NA, n.all)
            vertex.label[lextr] <- lextr
        }

        # Create and plot ggraph
        g <- ggraph(res$graph$adj.list, res$graph$weight.list)
        plot(g,
             vertex.color = all.vertex.color,
             vertex.size = all.vertex.size,
             vertex.label = vertex.label,
             vertex.shape = all.vertex.shape,
             ...)

        legend(legend.position,
               xpd = NA,
               legend = names(q$x.col.tbl),
               fill = q$x.col.tbl,
               inset = 0.01,
               cex = legend.cex,
               title = sprintf("Cell %d highlighted", i))

    } else if (type == "Ey") {

        # Color by function values
        q <- quantize.cont.var(res$Ey)
        vertex.color <- q$x.col.tbl[q$x.cat]

        # Get extrema
        lmin <- unique(res$MS_cx[,"local_min"])
        lmax <- unique(res$MS_cx[,"local_max"])
        n <- length(res$Ey)

        # Set vertex shapes
        vertex.shape <- rep("circle", n)
        vertex.shape[lmin] <- "csquare"
        vertex.shape[lmax] <- "square"

        # Set vertex sizes
        vertex.size.vec <- rep(vertex.size, n)
        vertex.size.vec[c(lmin, lmax)] <- lext.sfactor * vertex.size

        # Handle labels
        if (show.only.lextr.labels && is.null(vertex.label)) {
            lextr <- c(lmin, lmax)
            vertex.label <- rep(NA, n)
            vertex.label[lextr] <- lextr
        }

        # Create and plot ggraph
        g <- ggraph(res$graph$adj.list, res$graph$weight.list)
        plot(g,
             vertex.color = vertex.color,
             vertex.size = vertex.size.vec,
             vertex.label = vertex.label,
             vertex.shape = vertex.shape,
             ...)

        legend(legend.position,
               xpd = NA,
               legend = names(q$x.col.tbl),
               fill = q$x.col.tbl,
               inset = 0.01,
               cex = legend.cex,
               title = "Function values")
    }

    invisible(NULL)
}

#' Plot a k-h Matrix
#'
#' This function creates a heatmap-style plot of a k-h matrix, which represents
#' the presence or absence of a local maximum for different combinations of
#' k (number of nearest neighbors) and h (hop size) values.
#'
#' @param x Either a kh.matrix object (list with components kh.mat, existing.k, h.values, id)
#'   or a matrix where rows represent k values and columns represent h values.
#'   Values should be 0 (absent) or 1 (present).
#' @param existing.k A vector of k values corresponding to the rows (used when x is a matrix).
#' @param h.values A vector of h values corresponding to the columns (used when x is a matrix).
#' @param id A string identifier for the local maximum being plotted (used when x is a matrix).
#' @param color.palette A vector of two colors for the plot. Default is c("white", "black").
#' @param xlab Label for the x-axis. Default is "k (number of nearest neighbors)".
#' @param ylab Label for the y-axis. Default is "h (hop size)".
#' @param main Title of the plot. Default is "Presence of Local Maximum" followed by id.
#' @param with.legend Logical. If TRUE, a legend is added to the plot. Default is FALSE.
#' @param ... Additional arguments passed to image()
#'
#' @return A plot is created on the current graphics device.
#'
#' @details
#' The function creates a heatmap where each cell represents a combination of k and h values.
#' Black cells indicate the presence of a local maximum, while white cells indicate its absence.
#'
#' @examples
#' \dontrun{
#' # Create a kh.matrix object
#' kh.mat <- matrix(c(0,1,1,0,1,0,1,1,0), nrow = 3)
#' existing.k <- c(10, 20, 30)
#' h.values <- c(1, 2, 3)
#' khm <- kh.matrix(kh.mat, existing.k, h.values, id = "LM1")
#'
#' # Basic plot
#' plot(khm)
#'
#' # With custom colors and a legend
#' plot(khm,
#'      color.palette = c("lightblue", "darkred"),
#'      with.legend = TRUE)
#'
#' # You can also override parameters from the object
#' plot(khm, main = "Custom Title", xlab = "k values")
#' }
#'
#' @export
plot.kh.matrix <- function(x, existing.k = NULL, h.values = NULL, id = NULL,
                           color.palette = c("white", "black"),
                           xlab = "k (number of nearest neighbors)",
                           ylab = "h (hop size)",
                           main = NULL,
                           with.legend = FALSE,
                           ...) {
    ## Handle both list and matrix input formats
    if (is.list(x)) {
        kh.mat <- x$kh.mat
        if (is.null(existing.k)) existing.k <- x$existing.k
        if (is.null(h.values)) h.values <- x$h.values
        if (is.null(id)) id <- x$id
    } else {
        kh.mat <- x
        if (is.null(existing.k) || is.null(h.values) || is.null(id)) {
            stop("When x is a matrix, existing.k, h.values, and id must be provided")
        }
    }

    ## Set default main title if not provided
    if (is.null(main)) {
        main <- paste("Presence of Local Maximum", id)
    }

    x.vals <- 1:(max(existing.k) + 1)
    y.vals <- 1:(max(h.values) + 1)

    image(x = x.vals,
          y = y.vals,
          z = kh.mat,
          col = color.palette,  # when color.palette = c("white", "black"), then White is for 0, black is for 1
          xlab = xlab,
          ylab = ylab,
          main = main,
          las = 1,
          xlim = range(existing.k),
          ylim = range(h.values),
          ...)  # Pass additional arguments to image()

    if (with.legend) {
        legend("bottomright", legend = c("Present", "Absent"),
               fill = c("black", "white"), border = "black", bty = "n")
    }
}


#' Create Error Plot for Model Comparisons
#'
#' @description
#' Creates a plot showing total absolute true errors for different models,
#' with error bars representing the 95% credible intervals based on Bayesian
#' bootstrap samples. Supports both bar plot and point plot visualizations.
#'
#' @param x An object of class "model.errors" containing the following elements:
#'   \itemize{
#'     \item \code{integrals}: Named numeric vector of integral values
#'     \item \code{bb.integrals}: List of Bayesian bootstrap integral values for each model
#'   }
#' @param method Visualization method: "bars" or "points" (default: "bars")
#' @param title Plot title (optional)
#' @param y.lab Y-axis label (default: "Total Absolute True Error")
#' @param point.col Color of points when method="points" (default: "black")
#' @param bar.col Color of bars when method="bars" (default: "white")
#' @param border.col Color of bar borders (default: "black")
#' @param error.bar.col Color of error bars (default: "black")
#' @param error.bar.width Width of error bar ends (default: 0.1)
#' @param bar.width Width of the bars (default: 0.7)
#' @param point.size Size of points when method="points" (default: 1.2)
#' @param x.margin Margin to add on both sides of x-axis range (default: 0.5)
#' @param line1 Distance of x-axis labels from axis (default: 1)
#' @param line2 Distance of y-axis label from axis (default: 2.8)
#' @param line3 Distance of title from plot (default: 1)
#' @param with.vertical.lines Logical, whether to add vertical grid lines (default: TRUE)
#' @param vertical.line.col Color of vertical grid lines (default: "gray")
#' @param vertical.line.lty Line type for vertical grid lines (default: 1)
#' @param margin.bottom Extra space for x-axis labels (default: 8)
#' @param ... Additional graphical parameters passed to plotting functions
#'
#' @export
plot.model.errors <- function(x,
                              method = c("bars", "points"),
                              title = NULL,
                              y.lab = "Total Absolute True Error",
                              point.col = "black",
                              bar.col = "white",
                              border.col = "black",
                              error.bar.col = "black",
                              error.bar.width = 0.1,
                              bar.width = 0.7,
                              point.size = 1.2,
                              x.margin = 0.5,
                              line1 = 1,
                              line2 = 2.8,
                              line3 = 1,
                              with.vertical.lines = TRUE,
                              vertical.line.col = "gray",
                              vertical.line.lty = 1,
                              margin.bottom = 8,
                              ...) {

    # Extract components from x
    integrals <- x$integrals
    bb.integrals <- x$bb.integrals

    method <- match.arg(method)

    # Compute confidence intervals from bootstrap samples
    ci.intervals <- lapply(bb.integrals, function(x) {
        quantile(unlist(x), probs = c(0.025, 0.975))
    })

    # Convert to matrices for easier handling
    ci.mat <- do.call(rbind, ci.intervals)

    # Basic plot setup
    x <- seq_along(integrals)
    y <- unname(integrals)

    # Adjust y-axis range based on method
    if(method == "bars") {
        y.max <- max(ci.mat[,2]) * 1.05  # Add 5% padding at top
        y.min <- 0
    } else {
        y.range <- max(ci.mat[,2]) - min(ci.mat[,1])
        y.max <- max(ci.mat[,2]) + 0.05 * y.range
        y.min <- min(ci.mat[,1]) - 0.05 * y.range
    }

    # Create empty plot with margins
    plot(x, y, type = 'n', xlab = "", ylab = "",
         ylim = c(y.min, y.max),
         xlim = c(min(x) - x.margin, max(x) + x.margin),
         axes = FALSE)

    # Add visualization based on method
    if(method == "bars") {
        half_width <- bar.width/2
        for(i in seq_along(x)) {
            graphics::rect(x[i] - half_width, 0, x[i] + half_width, y[i],
                 col = bar.col, border = border.col)
            # Add error bars
            graphics::segments(x[i], ci.mat[i,1], x[i], ci.mat[i,2],
                    col = error.bar.col)
            # Add error bar ends
            graphics::segments(x[i] - error.bar.width, ci.mat[i,1],
                    x[i] + error.bar.width, ci.mat[i,1],
                    col = error.bar.col)
            graphics::segments(x[i] - error.bar.width, ci.mat[i,2],
                    x[i] + error.bar.width, ci.mat[i,2],
                    col = error.bar.col)
        }
    } else {  # points
        # First add error bars
        for(i in seq_along(x)) {
            if (with.vertical.lines) {
                graphics::abline(v = x[i], col = vertical.line.col, lty = vertical.line.lty)
            }
            graphics::segments(x[i], ci.mat[i,1], x[i], ci.mat[i,2],
                    col = error.bar.col)
            graphics::segments(x[i] - error.bar.width, ci.mat[i,1],
                    x[i] + error.bar.width, ci.mat[i,1],
                    col = error.bar.col)
            graphics::segments(x[i] - error.bar.width, ci.mat[i,2],
                    x[i] + error.bar.width, ci.mat[i,2],
                    col = error.bar.col)
        }
        # Then add points on top
        graphics::points(x, y, pch = 19, col = point.col, cex = point.size)
    }

    # Add axes and labels
    graphics::axis(2, las = 1)
    graphics::axis(1, at = x, labels = FALSE)

    # Add box
    graphics::box()

    # Add model names
    graphics::mtext(names(integrals), side = 1, at = x,
          line = line1, las = 2)

    # Add y-axis label
    graphics::mtext(y.lab, side = 2, line = line2)

    # Add title if provided
    if(!is.null(title)) {
        graphics::mtext(title, side = 3, line = line3)
    }
}

#' Plot Prediction Errors for Multiple Smoothing Models
#'
#' Creates a line plot comparing prediction errors across different smoothing models,
#' with customizable visual parameters. Each model's errors are represented by a line
#' with distinct color, line type, and point markers.
#'
#' @param x A prediction.errors object or a named list, data frame, or matrix containing
#'   prediction errors for different models. If a list, each element should be a numeric
#'   vector of errors and the names of the list elements are used in the plot legend.
#'   If a data frame or matrix, each column represents errors for a different model and
#'   column names are used in the legend.
#' @param xvals Numeric vector of x-axis values (typically independent variable values).
#'   Can be NULL if x contains an 'xvals' attribute or component. Default is NULL.
#' @param cols Character vector of colors for the lines. Defaults to a preset
#'   palette of 15 colors. Must be at least as long as the number of models.
#' @param ltys Numeric vector of line types. Defaults to a preset pattern of
#'   15 line types. Must be at least as long as the number of models.
#' @param main Character string for the plot title (default: "")
#' @param xlab Character string for x-axis label (default: "Gaussian Mean")
#' @param ylab Character string for y-axis label (default: "Mean Absolute True Error")
#' @param legend.pos Character string specifying legend position (default: "topleft")
#' @param legend.cex Numeric scaling factor for legend text size (default: 0.7)
#' @param legend.ncol Integer number of columns in legend (default: 2)
#' @param point.cex Numeric scaling factor for point size (default: 0.8)
#' @param ylab.line Numeric value specifying the distance of the y-axis label from the plot
#'   (in margin line units, default: 2.8)
#' @param legend.inset Numeric value specifying the inward shift of the legend from the
#'   plot border as a fraction of the plot region (default: 0.02)
#' @param offset.factor Numeric value controlling the spacing between points on
#'   different lines (default: floor(length(x) / 140)). Larger values create more
#'   separation between points across different models. Set to 0 to align points
#'   vertically across all lines.
#' @param n.points Integer number of points to plot per line (default: 10)
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns a list containing the colors and line types used in the plot
#'
#' @details The function plots multiple lines representing prediction errors for different
#' models. Each line is distinguished by color, line type, and point markers. Points are
#' added at regular intervals along each line with slight offsets to prevent overlapping.
#' If more than 15 models are provided, the function will throw an error unless custom
#' color and line type vectors of sufficient length are supplied.
#'
#' @examples
#' \dontrun{
#' # Using a named list
#' errors_list <- list(
#'   model1 = c(1,2,3,4,5),
#'   model2 = c(2,3,4,5,6)
#' )
#' xvals <- c(1,2,3,4,5)
#'
#' # Create prediction.errors object
#' pe <- prediction.errors(errors_list, xvals)
#'
#' # Basic plot
#' plot(pe)
#'
#' # Using a data frame
#' errors_df <- data.frame(
#'   model1 = c(1,2,3,4,5),
#'   model2 = c(2,3,4,5,6)
#' )
#'
#' # Create prediction.errors object from data frame
#' pe_df <- prediction.errors(errors_df, xvals = 1:5)
#'
#' # Plot with custom parameters
#' plot(pe_df,
#'      cols = c("red", "blue"),
#'      ltys = c(1, 2),
#'      main = "Prediction Errors",
#'      legend.pos = "topright")
#' }
#'
#' @export
plot.prediction.errors <- function(x,
                                   xvals = NULL,
                                   cols = NULL,
                                   ltys = NULL,
                                   main = "",
                                   xlab = "Gaussian Mean",
                                   ylab = "Mean Absolute True Error",
                                   legend.pos = "topleft",
                                   legend.cex = 0.7,
                                   legend.ncol = 2,
                                   point.cex = 0.8,
                                   ylab.line = 2.8,
                                   legend.inset = 0.02,
                                   offset.factor = NULL,
                                   n.points = 10,
                                   ...) {

    # Extract errors data and x-axis values
    if (inherits(x, "prediction.errors")) {
        # If it's a prediction.errors object, extract components
        errors_data <- x
        if (is.null(xvals)) {
            xvals <- attr(x, "xvals")
            if (is.null(xvals)) {
                stop("xvals must be provided either as a parameter or as an attribute of x")
            }
        }
    } else {
        # For backward compatibility with direct list/data.frame/matrix input
        errors_data <- x
    }

    # Input validation and conversion
    if (!is.list(errors_data) && !is.matrix(errors_data) && !is.data.frame(errors_data)) {
        stop("x must be a prediction.errors object, named list, data frame, or matrix with column names")
    }

    # Convert matrix to data frame
    if (is.matrix(errors_data)) {
        if (is.null(colnames(errors_data))) {
            stop("Matrix input must have column names")
        }
        errors_data <- as.data.frame(errors_data)
    }

    # Convert data frame to list if necessary
    if (is.data.frame(errors_data)) {
        errors_data <- as.list(errors_data)
    }

    # Validate xvals parameter
    if (is.null(xvals)) {
        stop("xvals must be provided either as a parameter or as an attribute of the prediction.errors object")
    }

    n.models <- length(errors_data)
    n.total <- length(xvals)

    # Set offset_factor if not provided
    if (is.null(offset.factor)) {
        offset.factor <- floor(n.total / 140)
    }

    # Define default colors and line types if not provided
    if (is.null(cols)) {
        cols <- c("red", "darkgreen", "purple", "blue", "chocolate4",
                 "darkturquoise", "darkorange3", "darkgrey", "violetred3",
                 "steelblue", "seagreen4", "salmon4", "slateblue4",
                 "indianred3", "olivedrab")
    }

    if (is.null(ltys)) {
        ltys <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 1)
    }

    # Check if we have enough colors and line types
    if (n.models > length(cols) || n.models > length(ltys)) {
        stop("Number of models exceeds available colors or line types. ",
             "Please supply longer cols and ltys vectors.")
    }

    # Subset colors and line types to match number of models
    cols <- cols[seq(n.models)]
    ltys <- ltys[seq(n.models)]

    # Calculate plot limits
    all.errors <- unlist(errors_data)
    all.errors <- all.errors[is.finite(all.errors)]
    ylim <- range(all.errors)

    # Create base plot
    plot(xvals, errors_data[[1]], type = "n",
         ylim = ylim, las = 1,
         main = main,
         xlab = xlab,
         ylab = "")
    mtext(ylab, side = 2, line = ylab.line, outer = FALSE)

    # Add lines and points for each model
    for (i in seq(errors_data)) {
        errors <- errors_data[[i]]

        # Draw line
        lines(xvals, errors, col = cols[i], lty = ltys[i])

        # Calculate indices for equally spaced points with offset
        offset <- offset.factor * (i-1)  # Different offset for each line
        indices <- seq(from = 1 + offset,
                      to = n.total,
                      length.out = n.points)
        indices <- round(indices)
        indices <- indices[indices <= length(xvals)]

        # Add points
        points(xvals[indices], errors[indices],
               col = cols[i],
               pch = i %% 25,
               cex = point.cex)
    }

    # Add legend
    legend(legend.pos,
           legend = names(errors_data),
           col = cols,
           lty = ltys,
           pch = seq_along(errors_data) %% 25,
           bg = "white",
           inset = legend.inset,
           cex = legend.cex,
           ncol = legend.ncol)

    # Return invisibly the plotting parameters used
    invisible(list(cols = cols, ltys = ltys))
}


#' Plot a Chain Graph with Optional Highlighted Paths
#'
#' @description
#' Creates a visualization of a chain graph where vertices are arranged horizontally.
#' The function can highlight specific paths and vertices based on path data,
#' and displays edge weights when they differ from 1.
#'
#' @param x An object of class "chain.with.path" containing the following elements:
#'   \itemize{
#'     \item \code{adj.list}: List of adjacency lists for each vertex, where each element
#'           contains the indices of neighboring vertices
#'     \item \code{weight.list}: List of corresponding weights for each adjacency,
#'           matching the structure of adj.list
#'     \item \code{gpd.obj}: List of path data objects, typically output from a
#'           get.path.data function. Each object should contain a 'vertices' element
#'           listing vertex indices in the path. The first object should have a
#'           'ref_vertex' element indicating a reference vertex.
#'   }
#' @param vertex.size Numeric value for the size of vertex points (default: 0.5)
#' @param margin Numeric value for the margin around the plot (default: 0.5)
#' @param title Character string for the plot title (default: "Chain Graph with Highlighted Path")
#' @param y.offset Numeric value for vertical offset of vertex labels (default: 0.1)
#' @param ... Additional arguments (currently unused)
#'
#' @details
#' The graph is plotted with vertices arranged horizontally at equal intervals.
#' Edge weights are displayed above the edges when they differ from 1.
#' Vertices in the paths specified by gpd.obj are highlighted with red circles,
#' and the reference vertex is highlighted with a larger blue circle.
#'
#' @return
#' No return value, called for side effects (plotting)
#'
#' @examples
#' # Create a simple chain graph
#' adj_list <- list(
#'   c(2),    # Vertex 1 connected to 2
#'   c(1, 3), # Vertex 2 connected to 1 and 3
#'   c(2)     # Vertex 3 connected to 2
#' )
#' weight_list <- list(
#'   c(1),    # Weight for edge 1-2
#'   c(1, 2), # Weights for edges 2-1 and 2-3
#'   c(2)     # Weight for edge 3-2
#' )
#' path_data <- list(
#'   list(vertices = c(1, 2, 3), ref_vertex = 1)
#' )
#'
#' chain.with.path.obj <- chain.with.path(adj_list, weight_list, path_data)
#' plot(chain.with.path.obj)
#'
#' @importFrom graphics plot points segments text
#' @export
plot.chain.with.path <- function(x,
                                 vertex.size = 0.5,
                                 margin = 0.5,
                                 title = "Chain Graph with Highlighted Path",
                                 y.offset = 0.1,
                                 ...) {
    # Extract components from x
    adj.list <- x$adj.list
    weight.list <- x$weight.list
    gpd.obj <- x$gpd.obj

    # Validate inputs
    if (length(adj.list) < 2) {
        stop("Chain graph must have at least 2 vertices")
    }

    # Set up plotting parameters
    n.vertices <- length(adj.list)
    x.coords <- seq_len(n.vertices)  # Horizontal positions
    y.coords <- rep(1, n.vertices)   # All vertices at same y-level

    # Set up plotting area
    plot(NA, NA,
         xlim = c(1 - margin, n.vertices + margin),
         ylim = c(1 - margin, 1 + margin),
         xlab = "", ylab = "",
         type = "n",
         xaxt = "n", yaxt = "n",
         main = title)

    # Draw edges
    for (i in seq_len(n.vertices)) {
        neighbors <- adj.list[[i]]
        for (j in neighbors) {
            if (j > i) {  # Draw each edge only once
                segments(x.coords[i], y.coords[i],
                        x.coords[j], y.coords[j],
                        lwd = 1)

                # Add weight label if not unit weight
                if (abs(weight.list[[i]][which(adj.list[[i]] == j)] - 1) > 1e-10) {
                    weight <- round(weight.list[[i]][which(adj.list[[i]] == j)], 2)
                    text(mean(c(x.coords[i], x.coords[j])),
                         mean(c(y.coords[i], y.coords[j])) + 0.1,
                         weight,
                         cex = 0.7)
                }
            }
        }
    }

    # Draw vertices
    points(x.coords, y.coords, pch = 19, cex = vertex.size)

    # Add vertex labels
    text(x.coords, y.coords - y.offset,
         labels = seq_len(n.vertices),
         cex = 0.8)

    # Highlight paths from get.path.data output
    if (!is.null(gpd.obj) && length(gpd.obj) > 0) {
        for (path in gpd.obj) {
            if (!is.null(path$vertices)) {
                # Draw red circles around vertices in the path
                points(x.coords[path$vertices],
                       y.coords[path$vertices],
                       col = "red",
                       pch = 1,
                       cex = 2 * vertex.size)
            }
        }
    }

    path <- gpd.obj[[1]]

    points(x.coords[path$ref_vertex], y.coords[path$ref_vertex], pch = 1, cex = 3 * vertex.size, col = "blue", lwd = 2)
}

#' 3D Graph Visualization with Height Values
#'
#' Visualizes a graph in 3D where vertices are positioned at different heights
#' based on a function's values. Uses an off-screen rgl device by default so it
#' works on headless/CI environments.
#'
#' @param x Graph input. Can be one of:
#'   \itemize{
#'     \item A \code{ggraph} object
#'     \item A plot result from \code{plot.ggraph()}
#'     \item A \code{graph.3d} object
#'     \item A list with components \code{graph} and \code{layout}
#'     \item A list with adjacency and weight list components (recognized names include
#'           \code{adj_list}, \code{adj.list}, \code{weight_list}, \code{weight.list}, etc.)
#' }
#' @param z Numeric vector of height values for each vertex (required unless x is graph.3d).
#' @param layout Layout algorithm name (character) or a numeric matrix of 3D coordinates (n x 2).
#' @param conn.points Logical; connect vertices with edges in 3D.
#' @param use.spheres Logical; use spheres (TRUE) or points (FALSE) for vertices.
#' @param graph.alpha Numeric in \eqn{[0,1]}; transparency for base graph.
#' @param z.point.size Numeric; size for z-value points/spheres.
#' @param z.color Vector of colors for z values or NULL for automatic coloring.
#' @param z.alpha Numeric in \eqn{[0,1]}; transparency for z points.
#' @param edge.color Color for edges.
#' @param edge.width Numeric; width of edges.
#' @param base.plane Logical; show base graph in x-y plane.
#' @param base.vertex.size Numeric; size of base plane vertices.
#' @param z.scale Numeric; scaling factor for z values.
#' @param show.axes Logical; show coordinate axes.
#' @param vertical.lines Logical; draw vertical lines from base to z points.
#' @param vertical.line.style "solid" or "dashed" for verticals.
#' @param dash.length, gap.length Numerics controlling dashed verticals.
#' @param vertical.line.color, vertical.line.alpha, vertical.line.width Aesthetics for verticals.
#' @param vertical.line.alpha Numeric in \eqn{[0,1]}; transparency for vertical guide lines.
#' @param vertical.line.width Numeric; line width for vertical guide lines.
#' @param basins.df Optional data.frame with columns \code{evertex}, \code{is_max}, and optional \code{label}.
#' @param evertex.sphere.radius Numeric; radius for extrema spheres.
#' @param evertex.min.color, evertex.max.color Colors for minima/maxima.
#' @param evertex.max.color Color used for maxima (extrema) spheres and labels.
#' @param evertex.cex Numeric; text size for extrema labels.
#' @param evertex.adj Numeric length-2; text adjustment for extrema labels.
#' @param evertex.label.offset Numeric; z-offset for extrema labels.
#' @param gap.length Numeric; length of gaps between dashes when \code{vertical.line.style = "dashed"}.
#' @param ... Reserved for future extensions.
#'
#' @return NULL (invisibly). Opens an rgl device, draws, and closes it on exit.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'
#'   # Example 1: Using ggraph object
#'   adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2), c(2))
#'   weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1), c(1))
#'   z.values <- c(0.5, 1.2, 0.8, 1.5)
#'   g <- ggraph(adj.list, weight.list)
#'   plot3D.graph(g, z.values)
#'
#'   # Example 2: Using list with adjacency and weight components directly
#'   my.graph <- list(adj_list = adj.list, weight_list = weight.list)
#'   plot3D.graph(my.graph, z.values)
#'
#'   # Example 3: Alternative naming convention
#'   my.graph2 <- list(adj.list = adj.list, weight.list = weight.list)
#'   plot3D.graph(my.graph2, z.values)
#' }
#' }
#' @export
plot3D.graph <- function(x,
                         z = NULL,
                         layout = "kk",
                         conn.points = FALSE,
                         use.spheres = TRUE,
                         graph.alpha = 0.7,
                         z.point.size = 0.5,
                         z.color = NULL,
                         z.alpha = 1,
                         edge.color = "gray70",
                         edge.width = 1,
                         base.plane = TRUE,
                         base.vertex.size = 0.5,
                         z.scale = 1,
                         show.axes = TRUE,
                         vertical.lines = TRUE,
                         vertical.line.style = "dashed",
                         dash.length = 0.05,
                         gap.length = 0.05,
                         vertical.line.color = "darkgray",
                         vertical.line.alpha = 0.5,
                         vertical.line.width = 0.5,
                         basins.df = NULL,
                         evertex.sphere.radius = 0.2,
                         evertex.min.color = "blue",
                         evertex.max.color = "red",
                         evertex.cex = 1,
                         evertex.adj = c(0.5, 0.5),
                         evertex.label.offset = 0.3,
                         ...) {

    ## rgl gating
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("This function requires the optional package 'rgl'. ",
             "Install with install.packages('rgl').", call. = FALSE)
    }

    ## Headless/CI-safe: render to the null device; harmless on desktops
    use_null <- (!interactive()) ||
        identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
        (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
    old_opt <- options(rgl.useNULL = use_null)
    on.exit(options(old_opt), add = TRUE)

    ## ---- open/clear device ----
    ## Check if an rgl device is already open
    if (rgl::cur3d() == 0) {
        ## No device open, create a new one
        if (use_null) {
            ## Null device for headless environments
            rgl::open3d()
            on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
        } else {
            ## Interactive: create a large square window
            ## Get screen dimensions
            screen_info <- try(rgl::par3d("windowRect"), silent = TRUE)

            ## Calculate square window size (use ~80% of screen height for safety)
            if (inherits(screen_info, "try-error")) {
                ## Fallback if we can't get screen info
                window_size <- 800
            } else {
                ## Estimate available screen space
                ## par3d("windowRect") returns current window, not screen size
                ## Use a reasonable maximum
                window_size <- min(1200, 800)  ## Conservative default
            }

            rgl::open3d(windowRect = c(50, 50, 50 + window_size, 50 + window_size))
        }
    } else {
        ## Device already open, use it (but don't close it on exit)
        rgl::set3d(rgl::cur3d())
    }

    rgl::clear3d()

    ## ---- extract/construct graph and 3D layout ----
    ## Extracting or constructing the graph object 'g'
    if (inherits(x, "graph.3d")) {
        g <- x$graph
        layout.2d <- x$layout
        if (is.null(z)) z <- x$z
        vertex.color <- if (!is.null(x$vertex.color)) x$vertex.color else rep("gray", igraph::vcount(g))

    } else if (inherits(x, "ggraph")) {
        graph.adj.list     <- x$adj.list
        graph.edge.lengths <- x$weight.list

        graph.obj <- convert.adjacency.to.edge.matrix(graph.adj.list, graph.edge.lengths)
        g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
        if (!is.null(graph.edge.lengths)) igraph::E(g)$weight <- graph.obj$weights

        layout.2d <- NULL  # Will be computed below
        vertex.color <- rep("gray", igraph::vcount(g))

    } else if (is.list(x) && !is.null(x$graph) && !is.null(x$layout)) {
        g <- x$graph
        layout.2d <- x$layout
        vertex.color <- if (!is.null(x$vertex.color)) x$vertex.color else rep("gray", igraph::vcount(g))

    } else if (is.list(x) && length(x) >= 2) {
        ## Try to identify adjacency list and weight list components
        ## Check for common naming patterns
        adj_names <- c("adj_list", "adj.list", "adjacency_list", "adjacency.list", "adj", "adjacency")
        weight_names <- c("weight_list", "weight.list", "weights", "weight", "edge_weights", "edge.weights")

        adj_idx <- which(names(x) %in% adj_names)[1]
        weight_idx <- which(names(x) %in% weight_names)[1]

        ## If named components not found, assume first two list components
        ## where both are lists themselves (heuristic for adjacency + weights)
        if (is.na(adj_idx) || is.na(weight_idx)) {
            list_components <- which(sapply(x, is.list))
            if (length(list_components) >= 2) {
                adj_idx <- list_components[1]
                weight_idx <- list_components[2]
                message("Using components '", names(x)[adj_idx], "' and '",
                        names(x)[weight_idx], "' as adjacency and weight lists")
            } else {
                stop("`x` appears to be a list but cannot identify adjacency and weight list components")
            }
        }

        graph.adj.list <- x[[adj_idx]]
        graph.edge.lengths <- x[[weight_idx]]

        ## Convert to igraph object
        graph.obj <- convert.adjacency.to.edge.matrix(graph.adj.list, graph.edge.lengths)
        g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
        if (!is.null(graph.edge.lengths)) igraph::E(g)$weight <- graph.obj$weights

        layout.2d <- NULL  # Will be computed below
        vertex.color <- rep("gray", igraph::vcount(g))

    } else {
        stop("`x` must be a graph.3d object, a ggraph object, ",
             "a list with components 'graph' and 'layout', ",
             "or a list with adjacency and weight list components")
    }

    ## Compute layout if not already provided
    if (is.null(layout.2d)) {
        if (is.character(layout)) {
            layout_func <- switch(layout,
                                  "nicely"    = igraph::layout_nicely,
                                  "randomly"  = igraph::layout_randomly,
                                  "in_circle" = igraph::layout_in_circle,
                                  "on_sphere" = igraph::layout_on_sphere,
                                  "fr"        = igraph::layout_with_fr,
                                  "kk"        = igraph::layout_with_kk,
                                  "as_tree"   = igraph::layout_as_tree,
                                  "lgl"       = igraph::layout_with_lgl,
                                  "graphopt"  = igraph::layout_with_graphopt,
                                  "drl"       = igraph::layout_with_drl,
                                  "dh"        = igraph::layout_with_dh,
                                  "mds"       = igraph::layout_with_mds,
                                  stop("Invalid layout specified")
                                  )
            layout.2d <- layout_func(g, dim = 2)
        } else if (is.matrix(layout) && ncol(layout) == 2) {
            layout.2d <- layout
        } else {
            stop("'layout' must be a known layout name or an n x 2 numeric matrix")
        }
    }

    ## ---- validate z and layout ----
    if (is.null(z)) {
        stop("`z` values must be provided (argument `z` or inside a graph.3d object)")
    }
    if (!is.numeric(z)) stop("`z` must be a numeric vector", call. = FALSE)
    layout.2d <- as.matrix(layout.2d)
    if (!is.numeric(layout.2d)) stop("`layout` must be numeric coordinates", call. = FALSE)
    if (ncol(layout.2d) != 2L) stop("`layout` must have exactly 2 columns", call. = FALSE)
    if (length(z) != nrow(layout.2d)) {
        stop(sprintf("Length of `z` (%d) must match number of vertices (%d).",
                     length(z), nrow(layout.2d)), call. = FALSE)
    }
    if (!is.numeric(graph.alpha) || graph.alpha < 0 || graph.alpha > 1) {
        stop("`graph.alpha` must be in [0, 1]", call. = FALSE)
    }
    if (!is.numeric(z.alpha) || z.alpha < 0 || z.alpha > 1) {
        stop("`z.alpha` must be in [0, 1]", call. = FALSE)
    }

    ## ---- scale & colors ----
    z <- as.numeric(z) * z.scale
    layout.3d <- cbind(layout.2d, z)

    if (is.null(z.color)) {
        if (diff(range(z)) > 0) {
            z.norm <- (z - min(z)) / diff(range(z))
            pal <- grDevices::hcl.colors(100, palette = "Spectral")
            idx <- pmin(100L, pmax(1L, 1L + floor(z.norm * 99)))
            z.color <- pal[idx]
        } else {
            z.color <- rep("red", length(z))
        }
    } else if (length(z.color) == 1L) {
        z.color <- rep(z.color, length(z))
    } else if (length(z.color) != length(z)) {
        stop("`z.color` must be a single color or have the same length as `z`", call. = FALSE)
    }

    ## ---- base plane (z=0) ----
    if (isTRUE(base.plane)) {
        base.layout <- cbind(layout.2d, 0)

        if (isTRUE(use.spheres)) {
            rgl::spheres3d(
                     base.layout,
                     radius = base.vertex.size * 0.02,
                     col = vertex.color,
                     alpha = graph.alpha
                 )
        } else {
            rgl::points3d(
                     base.layout,
                     size = base.vertex.size * 10,
                     col = vertex.color,
                     alpha = graph.alpha
                 )
        }

        edges <- igraph::as_edgelist(g, names = FALSE)
        if (NROW(edges) > 0) {
            edge.coords <- matrix(NA_real_, nrow = NROW(edges) * 2L, ncol = 3L)
            for (i in seq_len(NROW(edges))) {
                v1 <- edges[i, 1]; v2 <- edges[i, 2]
                edge.coords[(i - 1L) * 2L + 1L, ] <- base.layout[v1, ]
                edge.coords[(i - 1L) * 2L + 2L, ] <- base.layout[v2, ]
            }
            rgl::segments3d(edge.coords, col = "gray50", alpha = graph.alpha, lwd = edge.width)
        }
    }

    ## ---- z points ----
    if (isTRUE(use.spheres)) {
        rgl::spheres3d(
                 layout.3d,
                 radius = z.point.size * 0.02,
                 col = z.color,
                 alpha = z.alpha
             )
    } else {
        rgl::points3d(
                 layout.3d,
                 size = z.point.size * 10,
                 col = z.color,
                 alpha = z.alpha
             )
    }

    ## ---- 3D edges between z points ----
    if (isTRUE(conn.points)) {
        edges <- igraph::as_edgelist(g, names = FALSE)
        if (NROW(edges) > 0) {
            for (i in seq_len(NROW(edges))) {
                v1 <- edges[i, 1]; v2 <- edges[i, 2]
                rgl::lines3d(rbind(layout.3d[v1, ], layout.3d[v2, ]),
                             col = edge.color, lwd = edge.width)
            }
        }
    }

    ## ---- vertical lines from base to z points ----
    if (isTRUE(vertical.lines) && isTRUE(base.plane)) {
        if (identical(vertical.line.style, "dashed")) {
            for (i in seq_len(nrow(layout.2d))) {
                z_height <- z[i]
                if (z_height <= 0) next
                seg_len <- dash.length + gap.length
                n_segments <- floor(z_height / seg_len)
                if (n_segments < 1L) n_segments <- 1L
                for (j in 0:(n_segments - 1L)) {
                    start_z <- j * seg_len
                    end_z   <- min(start_z + dash.length, z_height)
                    rgl::lines3d(
                             rbind(c(layout.2d[i, ], start_z),
                                   c(layout.2d[i, ], end_z)),
                             col   = vertical.line.color,
                             alpha = vertical.line.alpha,
                             lwd   = vertical.line.width
                         )
                }
            }
        } else {
            base.layout <- cbind(layout.2d, 0)
            for (i in seq_len(nrow(layout.2d))) {
                rgl::lines3d(
                         rbind(base.layout[i, ], layout.3d[i, ]),
                         col   = vertical.line.color,
                         alpha = vertical.line.alpha,
                         lwd   = vertical.line.width
                     )
            }
        }
    }

    ## ---- extrema (basins.df) ----
    if (!is.null(basins.df)) {
        if (!all(c("evertex", "is_max") %in% names(basins.df))) {
            stop("`basins.df` must contain columns 'evertex' and 'is_max'", call. = FALSE)
        }
        labels <- if (!is.null(basins.df$label)) basins.df$label else paste0("E", seq_len(nrow(basins.df)))
        for (i in seq_len(nrow(basins.df))) {
            v <- basins.df$evertex[i]
            if (!is.finite(v) || v < 1L || v > nrow(layout.3d)) next
            colv <- if (isTRUE(basins.df$is_max[i] == 1)) evertex.max.color else evertex.min.color
            rgl::spheres3d(layout.3d[v, 1], layout.3d[v, 2], layout.3d[v, 3],
                           radius = evertex.sphere.radius * 0.02, col = colv, alpha = 0.8)
            rgl::text3d(layout.3d[v, 1], layout.3d[v, 2], layout.3d[v, 3] + evertex.label.offset,
                        texts = labels[i], col = colv, cex = evertex.cex, adj = evertex.adj)
        }
    }

    ## ---- axes ----
    if (isTRUE(show.axes)) {
        rgl::axes3d()
        rgl::title3d(xlab = "X", ylab = "Y", zlab = "Z")
    }

    invisible(NULL)
}

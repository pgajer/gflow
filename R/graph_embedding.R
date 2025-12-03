#' Embed Graph in 2D or 3D Space
#'
#' @description
#' Computes a 2D or 3D layout for a graph using force-directed algorithms.
#' The function takes an adjacency list representation of a graph and produces
#' an embedding using either the Fruchterman-Reingold or Kamada-Kawai algorithm.
#'
#' @param adj.list A list representing the adjacency list of the graph. Each element
#'   is an integer vector of neighbor indices.
#' @param weights.list An optional list of numeric vectors representing edge weights.
#'   If provided, each element should correspond to the weights of the edges in the
#'   same position in adj.list. Must have the same length as adj.list, with matching
#'   element lengths.
#' @param invert.weights Logical. If TRUE (default), inverts the weights (1/weight)
#'   before applying them to the Fruchterman-Reingold algorithm. Has no effect on
#'   Kamada-Kawai algorithm or when weights are not provided. See Details.
#' @param dim Integer, either 2 or 3, specifying the dimension of the embedding.
#'   Default is 2.
#' @param method Character string specifying the layout algorithm to use.
#'   Either "fr" for Fruchterman-Reingold or "kk" for Kamada-Kawai.
#'   Default is "fr".
#' @param verbose Logical. If TRUE, print progress and timing information.
#'   Default is FALSE.
#'
#' @return A numeric matrix with n.vertices rows and dim columns, where each row
#'   contains the coordinates of a vertex in the embedded space. For empty graphs,
#'   returns a matrix with 0 rows and dim columns.
#'
#' @details
#' This function implements force-directed graph layout algorithms to position
#' vertices in low-dimensional space. The goal is to position connected vertices
#' close to each other while maintaining appropriate spacing based on edge weights.
#'
#' The two layout algorithms interpret weights differently. When weights are provided,
#' they are handled as follows. The Fruchterman-Reingold algorithm interprets weights
#' as spring constants, where higher weights create stronger springs that pull vertices
#' closer together. By default, weights are inverted (1/weight) so that higher original
#' weights result in greater distances between vertices, matching the intuition of
#' weights as distances or dissimilarities. The Kamada-Kawai algorithm interprets
#' weights directly as distances, where higher weights mean greater distances between
#' vertices. The invert.weights parameter controls weight inversion for the
#' Fruchterman-Reingold algorithm. When TRUE (default), it inverts the weights to
#' match the distance interpretation. When FALSE, it uses the weights as-is, where
#' higher weights pull vertices closer together.
#'
#' For graphs with no edges, the function returns a matrix of random positions
#' uniformly distributed in the range \eqn{[-1, 1]} for each dimension.
#'
#' @examples
#' ## Unweighted graph
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' embedding <- graph.embedding(adj.list, dim = 2, method = "fr")
#'
#' ## Weighted graph with distance interpretation
#' weights.list <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
#' embedding.weighted <- graph.embedding(adj.list, weights.list, dim = 2, method = "kk")
#'
#' ## 3D embedding
#' embedding.3d <- graph.embedding(adj.list, dim = 3, method = "fr", verbose = TRUE)
#'
#' @importFrom igraph graph_from_edgelist E layout_with_fr layout_with_kk set_edge_attr
#'
#' @export
graph.embedding <- function(adj.list,
                            weights.list = NULL,
                            invert.weights = TRUE,
                            dim = 2,
                            method = c("fr", "kk"),
                            verbose = FALSE) {

  ## Check igraph availability
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it.")
  }

  ## Validate inputs
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }

  if (!is.numeric(dim) || length(dim) != 1 || dim %% 1 != 0 || !(dim %in% c(2, 3))) {
    stop("dim must be either 2 or 3")
  }

  method <- match.arg(method)

  if (!is.logical(invert.weights) || length(invert.weights) != 1) {
    stop("invert.weights must be a single logical value")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  ## Validate weights.list if provided
  if (!is.null(weights.list)) {
    if (!is.list(weights.list) || length(weights.list) != length(adj.list)) {
      stop("weights.list must be a list with the same length as adj.list")
    }
    if (!all(mapply(function(a, w) length(a) == length(w), adj.list, weights.list))) {
      stop("Each element of weights.list must have the same length as the corresponding element of adj.list")
    }
    ## Convert to numeric (handles both integer and numeric inputs)
    weights.list <- lapply(weights.list, as.numeric)
  }

  ## Number of vertices
  n.vertices <- length(adj.list)

  ## Handle empty graph
  if (n.vertices == 0) {
    return(matrix(0, nrow = 0, ncol = dim))
  }

  ## Start timing if verbose
  if (verbose) {
    routine.ptm <- proc.time()
    ptm <- proc.time()
    cat("Converting graph adjacency list to edge matrix ... ")
  }

  ## Convert adjacency list to edge matrix
  res <- convert.adjacency.to.edge.matrix(adj.list, weights.list)

  if (verbose) {
    elapsed.time(ptm)
  }

  ## Handle graph with no edges
  if (nrow(res$edge.matrix) == 0) {
    if (verbose) {
      cat("Graph has no edges. Returning random positions.\n")
    }
    return(matrix(runif(n.vertices * dim, -1, 1), nrow = n.vertices, ncol = dim))
  }

  ## Create igraph object
  if (verbose) {
    ptm <- proc.time()
    cat("Creating igraph object from edge matrix ... ")
  }

  g <- igraph::graph_from_edgelist(res$edge.matrix, directed = FALSE)

  if (verbose) {
    elapsed.time(ptm)
  }

  ## Add edge weights if provided
  if (!is.null(weights.list)) {
    if (method == "fr" && invert.weights) {
      ## For FR algorithm, invert weights so higher weights = greater distances
      igraph::E(g)$weight <- 1 / res$weights
    } else {
      ## For KK algorithm or FR without inversion, use weights as-is
      igraph::E(g)$weight <- res$weights
    }
  }

  ## Compute layout
  if (verbose) {
    ptm <- proc.time()
    cat(sprintf("Computing %s layout in %dD space ... ",
                toupper(method), dim))
  }

  layout.coords <- switch(method,
                         fr = igraph::layout_with_fr(g, dim = dim),
                         kk = igraph::layout_with_kk(g, dim = dim))

  if (verbose) {
    elapsed.time(ptm)
    txt <- "Total elapsed time"
    elapsed.time(routine.ptm, txt, with.brackets = FALSE)
  }

  return(layout.coords)
}

#' @title Plot a Graph with Colored Vertices
#' @description Creates a visualization of a graph where vertices are colored according
#'              to a numeric function value, using base R graphics. The function displays
#'              the graph structure with edges as lines and vertices as colored points,
#'              with an optional color scale legend.
#'
#' @param embedding A numeric matrix with dimensions n x 2, where n is the number of
#'                 vertices. Contains the 2D coordinates for each vertex, typically
#'                 generated by a graph embedding algorithm.
#' @param adj.list A list of length n, where each element i contains a numeric vector
#'                of indices representing the vertices adjacent to vertex i.
#' @param vertex.colors A numeric vector of length n containing the function values
#'                     used to color the vertices.
#' @param vertex.size Numeric scalar controlling the size of vertices (default: 1).
#'                   Uses the same scale as base R's cex parameter
#' @param edge.alpha Numeric scalar between 0 and 1 controlling edge transparency
#'                  (default: 0.2). Higher values make edges more opaque.
#' @param color.palette Optional vector of colors defining the color gradient for
#'                     vertices. If NULL (default), uses a blue-white-red gradient.
#' @param main Character string for the plot title (default: "").
#' @param add.legend Logical indicating whether to add a color scale legend
#'                  (default: TRUE).
#'
#' @return No return value; produces a plot as a side effect.
#'
#' @details The function creates a 2D visualization of a graph structure where:
#'          - Edges are drawn as semi-transparent lines
#'          - Vertices are drawn as colored points
#'          - Colors are mapped to vertex_colors values using a continuous gradient
#'          - The aspect ratio is maintained at 1:1
#'          - Axes and tick marks are suppressed
#'          An optional color scale legend can be added to interpret vertex colors.
#'
#' @examples
#' \dontrun{
#' # Basic usage with simulated data
#' n <- 100  # number of vertices
#' embedding <- matrix(rnorm(2*n), ncol=2)
#' adj.list <- lapply(1:n, function(i) sample(1:n, 5))
#' vertex.colors <- rnorm(n)
#' plot.colored.graph(embedding, adj.list, vertex.colors)
#'
#' # Custom styling
#' plot2D.colored.graph(
#'     embedding,
#'     adj.list,
#'     vertex.colors,
#'     vertex.size = 1.5,
#'     edge.alpha = 0.3,
#'     color.palette = colorRampPalette(c("navy", "white", "darkred"))(100),
#'     main = "My Graph",
#'     add.legend = TRUE
#' )
#' }
#'
#' @note The function temporarily modifies graphical parameters using par()
#'       but restores them before exiting.
#'
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom graphics strwidth text
#' @export
plot2D.colored.graph <- function(embedding, adj.list, vertex.colors,
                               vertex.size = 1,
                               edge.alpha = 0.2,
                               color.palette = NULL,
                               main = "",
                               add.legend = TRUE) {

    # If no color palette is provided, create a default one
    if (is.null(color.palette)) {
        # Create a blue to red color palette
        cols <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
        # Scale the vertex colors to 1-100 for color mapping
        color.indices <- round((vertex.colors - min(vertex.colors)) /
                             (max(vertex.colors) - min(vertex.colors)) * 99 + 1)
        point.colors <- cols[color.indices]
    } else {
        # Use provided color palette
        cols <- color.palette
        color.indices <- round((vertex.colors - min(vertex.colors)) /
                             (max(vertex.colors) - min(vertex.colors)) * (length(cols) - 1) + 1)
        point.colors <- cols[color.indices]
    }

    # Calculate plot margins
    mar.right <- if(add.legend) 4 else 1
    par(mar = c(1, 1, 2, mar.right))

    # Create empty plot
    plot(embedding[,1], embedding[,2],
         type = "n",
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         main = main,
         asp = 1)  # Keep aspect ratio 1:1

    # Draw edges
    # Convert edge alpha to hex color
    edge.col <- grDevices::rgb(0, 0, 0, edge.alpha)

    # Draw all edges
    for(i in seq_along(adj.list)) {
        if(length(adj.list[[i]]) > 0) {
            segments(embedding[i,1], embedding[i,2],
                    embedding[adj.list[[i]],1], embedding[adj.list[[i]],2],
                    col = edge.col)
        }
    }

    # Draw vertices
    points(embedding[,1], embedding[,2],
           pch = 19,
           cex = vertex.size,
           col = point.colors)

    # Add legend if requested
    if(add.legend) {
        # Create legend values
        legend.vals <- round(seq(min(vertex.colors), max(vertex.colors), length.out = 5), 2)
        legend.cols <- cols[round(seq(1, length(cols), length.out = 5))]

        # Add color scale legend
        par(xpd = TRUE)  # Allow plotting in margins
        legend.x <- par("usr")[2] * 1.02  # Just outside the plot
        legend.y <- mean(par("usr")[3:4])

        # Draw color bar
        gradient.bars <- length(cols)
        bar.height <- (par("usr")[4] - par("usr")[3]) / gradient.bars

        for(i in 1:gradient.bars) {
            rect(legend.x,
                 par("usr")[3] + (i-1) * bar.height,
                 legend.x + graphics::strwidth("M"),  # Width of one character
                 par("usr")[3] + i * bar.height,
                 col = cols[i],
                 border = NA)
        }

        # Add text labels
        graphics::text(legend.x + graphics::strwidth("M") * 1.5,
             seq(par("usr")[3], par("usr")[4], length.out = 5),
             labels = legend.vals,
             adj = 0,
             cex = 0.8)
    }

    # Reset par
    par(xpd = FALSE)
}

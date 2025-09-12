##
## Graph related functions
##

#' Create a 2D or 3D Embedding of a Graph
#'
#' This function takes an adjacency list representation of a graph and produces
#' a 2D or 3D embedding using either the Fruchterman-Reingold or
#' Kamada-Kawai algorithm. It can handle both weighted and unweighted graphs.
#'
#' @param adj.list A list representing the adjacency list of the graph.
#' @param weights.list An optional list of numeric vectors representing the weights
#'   of the edges. If provided, each element should correspond to the weights of
#'   the edges in the same position in adj.list.
#' @param invert.weights Logical, if TRUE (default), inverts the weights (1/weight)
#'   before applying them to the Fruchterman-Reingold algorithm. Has no effect on
#'   Kamada-Kawai algorithm or when weights are not provided.
#' @param dim Integer, either 2 or 3, specifying the dimension of the embedding.
#'   Default is 3.
#' @param method Character string specifying the layout algorithm to use.
#'   Either "fr" for Fruchterman-Reingold or "kk" for Kamada-Kawai.
#'   Default is "fr".
#' @param verbose Logical, if TRUE (default), print progress messages.
#'
#' @return A matrix representing the graph embedding, with each row
#'   corresponding to a vertex and containing its coordinates.
#'
#' @details
#' The function can handle both weighted and unweighted graphs. When weights are provided,
#' they are interpreted differently by the two layout algorithms:
#'
#' - Fruchterman-Reingold (FR): Weights are interpreted as spring constants. Higher
#'   weights mean stronger springs, pulling vertices closer together. By default,
#'   weights are inverted (1/weight) so that higher original weights result in
#'   greater distances between vertices, matching the intuition of weights as distances.
#'
#' - Kamada-Kawai (KK): Weights are interpreted as distances. Higher weights mean
#'   greater distances between vertices.
#'
#' The `invert.weights` parameter allows control over this behavior for the FR algorithm.
#' When TRUE (default), it inverts the weights for FR to match the distance interpretation.
#' When FALSE, it uses the weights as-is, where higher weights pull vertices closer.
#'
#' @importFrom igraph graph_from_edgelist E layout_with_fr layout_with_kk
#'
#' @examples
#' # Unweighted graph
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' embedding <- graph.embedding(adj.list, dim = 2, method = "fr")
#'
#' # Weighted graph
#' weights.list <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
#' embedding_weighted <- graph.embedding(adj.list, weights.list, dim = 2, method = "kk")
#'
#' @export
graph.embedding <- function(adj.list,
                            weights.list = NULL,
                            invert.weights = TRUE,
                            dim = 3,
                            method = c("fr", "kk"),
                            verbose = TRUE) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required. Please install it.")
    }

    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    if (!is.numeric(dim) || length(dim) != 1 || dim %% 1 != 0 || !(dim %in% c(2, 3))) {
        stop("dim must be either 2 or 3")
    }

    method <- match.arg(method)

    if (!is.logical(verbose)) {
        stop("verbose must be a logical value")
    }

    # Convert integer weights to numeric if necessary
    if (!is.null(weights.list)) {
        if (!is.list(weights.list) || length(weights.list) != length(adj.list)) {
            stop("weights.list must be a list with the same length as adj.list")
        }
        if (!all(mapply(function(a, w) length(a) == length(w), adj.list, weights.list))) {
            stop("Each element of weights.list must have the same length as the corresponding element of adj.list")
        }
        # Convert to numeric (this will work for both integer and numeric inputs)
        weights.list <- lapply(weights.list, as.numeric)
    }

    if (verbose) {
        routine.ptm <- proc.time()
        ptm <- proc.time()
        cat("Converting the graph adjacency list to edge matrix ... ")
    }

    res <- convert.adjacency.to.edge.matrix(adj.list, weights.list)


    if (verbose) {
        elapsed.time(ptm)
        ptm <- proc.time()
        cat("Creating an igraph graph from the edge matrix ... ")
    }

    g <- igraph::graph_from_edgelist(res$edge.matrix, directed = FALSE)

    if (!is.null(weights.list)) {

        if (invert.weights) {
            igraph::E(g)$weight = 1 / res$weights
        } else {
            igraph::E(g)$weight = res$weights
        }
    }

    if (verbose) {
        elapsed.time(ptm)
        ptm <- proc.time()
        cat("Creating a graph layout ... ")
    }

    graph.embedding <- switch(method,
                              fr = igraph::layout_with_fr(g, dim = dim),
                              kk = igraph::layout_with_kk(g, dim = dim)
                              )
    if (verbose) {
        elapsed.time(ptm)
        txt <- sprintf("Total elapsed time")
        elapsed.time(routine.ptm, txt, with.brackets = FALSE)
    }

    return(graph.embedding)
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



#' Create Pruned Intersection Size List
#'
#' This function generates a list of intersection sizes for the edges of a pruned adjacency list.
#' It takes the original adjacency list, intersection size list, and pruned adjacency list as inputs,
#' and returns a new list containing intersection sizes corresponding to the pruned graph structure.
#'
#' @param adj.list A list where each element represents a node and contains indices of its neighbors
#'   in the original graph.
#' @param isize.list A list of the same length as adj.list, where each element contains intersection
#'   sizes corresponding to the edges in adj.list.
#' @param pruned.adj.list A list representing the pruned graph, where each element contains indices
#'   of neighbors after removing redundant edges.
#'
#' @return A list of the same length as pruned.adj.list, where each element contains intersection
#'   sizes corresponding to the edges in the pruned graph.
#'
#' @details
#' The function iterates through each node in the pruned adjacency list. For each node, it identifies
#' the neighbors in the pruned graph and retrieves their corresponding intersection sizes from the
#' original isize.list. The resulting pruned.isize.list maintains the structure of pruned.adj.list
#' but contains intersection sizes instead of neighbor indices.
#'
#' @note
#' - The function assumes that all nodes present in pruned.adj.list are also present in adj.list
#'   and isize.list.
#' - If a neighbor in pruned.adj.list is not found in the original adj.list (which should not happen
#'   under normal circumstances), the corresponding intersection size will be set to NA.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' adj.list <- list(c(2,3,4), c(1,3), c(1,2,4), c(1,3))
#' isize.list <- list(c(2,1,3), c(2,1), c(1,1,2), c(3,2))
#' pruned.adj.list <- list(c(2,4), c(1), c(4), c(1,3))
#' pruned.isize.list <- create.pruned.isize.list(adj.list, isize.list, pruned.adj.list)
#' print(pruned.isize.list)
#' }
#'
#' @export
create.pruned.isize.list <- function(adj.list, isize.list, pruned.adj.list) {
    ## Initialize an empty list to store the pruned intersection sizes
    pruned.isize.list <- vector("list", length(pruned.adj.list))

    ## Iterate through each node in the pruned adjacency list
    for (i in seq_along(pruned.adj.list)) {
        ## Get the neighbors of the current node in the pruned graph
        pruned.neighbors <- pruned.adj.list[[i]]

        ## Initialize a vector to store intersection sizes for the current node
        pruned.isizes <- numeric(length(pruned.neighbors))

        ## Iterate through each neighbor in the pruned graph
        for (j in seq_along(pruned.neighbors)) {
            ## Get the index of the neighbor in the original adjacency list
            neighbor.index <- which(adj.list[[i]] == pruned.neighbors[j])

            ## If the neighbor is found, get its intersection size from isize.list
            if (length(neighbor.index) > 0) {
                pruned.isizes[j] <- isize.list[[i]][neighbor.index]
            } else {
                ## If not found (shouldn't happen), set to NA or 0
                pruned.isizes[j] <- NA
            }
        }

        ## Store the pruned intersection sizes for the current node
        pruned.isize.list[[i]] <- pruned.isizes
    }

    return(pruned.isize.list)
}

#' Prune Long Edges in a Weighted Graph
#'
#' This function prunes long edges in a weighted graph based on the existence of shorter alternative paths.
#' It uses a C++ implementation for efficiency and returns the pruned graph along with information about the pruning process.
#'
#' @param graph A list of integer vectors representing the adjacency list of the graph.
#'   Each element of the list corresponds to a vertex, and contains the indices of its neighboring vertices.
#'   The graph should use 1-based indexing (as is standard in R).
#'
#' @param edge.lengths A list of numeric vectors representing the lengths of edges.
#'   Each element corresponds to a vertex, and contains the lengths of edges to its neighbors.
#'   The structure should match that of the `graph` parameter.
#'
#' @param alt.path.len.ratio.thld A single numeric value representing the threshold for the
#'   alternative path length ratio. Edges are pruned if an alternative path is found with
#'   length less than this ratio times the original edge length.
#'
#' @param use.total.length.constraint A logical value. If TRUE, the total length of the alternative path
#'   must be less than the original edge length. If FALSE, each edge in the alternative path
#'   must be shorter than the original edge length. Default is TRUE.
#'
#' @param verbose A logical value. If TRUE, progress information will be printed during the pruning process.
#'   Default is FALSE.
#'
#' @return A list containing four elements:
#'   \item{adj_list}{A list of integer vectors representing the adjacency list of the pruned graph.}
#'   \item{edge_lengths_list}{A list of numeric vectors representing the edge lengths in the pruned graph.}
#'   \item{path_lengths}{A numeric vector of alternative path lengths found during pruning.}
#'   \item{edge_lengths}{A numeric vector of original edge lengths corresponding to path_lengths.}
#'
#' @details
#' The function first converts the input graph to 0-based indexing for the C++ function,
#' then calls the C++ implementation to perform the pruning. The pruned graph is converted
#' back to 1-based indexing before being returned.
#'
#' The pruning process iterates through edges from longest to shortest. For each edge,
#' it searches for an alternative path. If a path is found that is shorter than the threshold
#' ratio times the original edge length, the edge is removed from the graph.
#'
#' @note
#' - The function assumes the input graph is undirected.
#' - The input graph and edge lengths must be consistent (same length and structure).
#' - The function may modify the order of vertices in the adjacency lists.
#'
#' @examples
#' \dontrun{
#' # Create a simple weighted graph
#' graph <- list(c(2,3), c(1,3), c(1,2))
#' edge.lengths <- list(c(1,2), c(1,3), c(2,3))
#' threshold <- 0.9
#'
#' # Prune the graph
#' result <- wgraph.prune.long.edges(graph, edge.lengths, threshold,
#'                                   use.total.length.constraint = TRUE,
#'                                   verbose = TRUE)
#' # Examine the results
#' print(result$adj_list)
#' print(result$edge_lengths_list)
#' print(result$path_lengths)
#' print(result$edge_lengths)
#' }
#' @export
wgraph.prune.long.edges <- function(graph,
                                    edge.lengths,
                                    alt.path.len.ratio.thld,
                                    use.total.length.constraint = TRUE,
                                    verbose = FALSE) {

    if (!is.list(graph) || !all(sapply(graph, is.integer)))
        stop("graph must be a list of integer vectors")
    if (!is.list(edge.lengths) || !all(sapply(edge.lengths, is.numeric)))
        stop("edge.lengths must be a list of numeric vectors")
    if (length(graph) != length(edge.lengths))
        stop("graph and edge.lengths must have the same length")
    if (!is.numeric(alt.path.len.ratio.thld) || length(alt.path.len.ratio.thld) != 1)
        stop("alt.path.len.ratio.thld must be a single numeric value")
    if (!is.logical(use.total.length.constraint) || length(use.total.length.constraint) != 1)
        stop("use.total.length.constraint must be a single logical value")
    if (!is.logical(verbose) || length(verbose) != 1)
        stop("verbose must be a single logical value")

    ## Converting graph to 0-based indexing
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))

    result <- .Call("S_wgraph_prune_long_edges",
                    graph.0based,
                    edge.lengths,
                    as.numeric(alt.path.len.ratio.thld),
                    as.logical(use.total.length.constraint),
                    as.logical(verbose))

    ## Converting the pruned graph back to 1-based indexing
    result$adj_list <- lapply(result$adj_list, function(x) x + 1)

    return(result)
}


#' Compare Two Adjacency Lists
#'
#' This function compares two adjacency lists to check if they are equivalent,
#' meaning each vertex's list of neighbors in the first adjacency list contains
#' exactly the same set of neighbors as in the corresponding list in the second
#' adjacency list, regardless of order.
#'
#' @param adj.list1 A list of integer vectors, where each vector represents the
#'                  adjacency list for a vertex in the first graph.
#' @param adj.list2 A list of integer vectors, where each vector represents the
#'                  adjacency list for a vertex in the second graph.
#'
#' @return A logical value; `TRUE` if the two adjacency lists are equivalent
#'         for all corresponding vertices, otherwise `FALSE`.
#'
#' @details
#' The function performs a length check on the adjacency lists and then uses
#' `setequal` for each corresponding pair of sub-lists to verify equivalence
#' of the adjacency sets. This function is useful for verifying the equality of
#' graphs in terms of connectivity, ignoring the order of nodes in the adjacency
#' lists and potential multiple edges between nodes.
#'
#' @examples
#' adj.list1 <- list(c(2, 3), c(1, 3), c(1, 2))
#' adj.list2 <- list(c(3, 2), c(3, 1), c(2, 1))
#' compare.adj.lists(adj.list1, adj.list2) # returns TRUE
#'
#' adj.list1 <- list(c(2, 3, 4), c(1, 3), c(1, 2), c(1))
#' adj.list2 <- list(c(3, 2), c(3, 1), c(2, 1), c(1))
#' compare.adj.lists(adj.list1, adj.list2) # returns FALSE
#'
#' @export
compare.adj.lists <- function(adj.list1, adj.list2) {
  if (length(adj.list1) != length(adj.list2)) {
    return(FALSE)
  }

  vertices.with.different.neighbors <- c()
  for (i in seq_along(adj.list1)) {
      if (!setequal(adj.list1[[i]], adj.list2[[i]])) {
          vertices.with.different.neighbors <- c(vertices.with.different.neighbors,i)
    }
  }

  if (length(vertices.with.different.neighbors) == 0) {
      return(TRUE)
  } else {
      cat("FALSE: Vertices with neighbor sets not the same: ")
      cat(vertices.with.different.neighbors)
      cat("\n")
      return(FALSE)
  }
}

#' Converts a graph adjacency list to an adjacency matrix
#'
#' This function converts a given graph adjacency list into an adjacency matrix,
#' where each entry in the matrix indicates whether there is an edge between the corresponding nodes.
#'
#' @param adj.list A list representing the adjacency list of a graph. Each element of the list
#'        is a vector of node indices that the corresponding node is connected to.
#' @param weights.list (optional) A list of the same length as `adj.list`, where each element
#'        is a vector of weights for the edges connecting the corresponding node to the nodes in `adj.list`.
#'        If not provided, the adjacency matrix will be binary (1 for edges and 0 for no edges).
#' @return A matrix representing the adjacency matrix of the graph.
#' @examples
#' adj.list <- list(c(2, 4), c(1, 3), c(2, 4), c(1, 3))
#' weights.list <- list(c(2, 4), c(2, 3), c(3, 1), c(4, 1))
#' convert.adjacency.list.to.adjacency.matrix(adj.list)
#' convert.adjacency.list.to.adjacency.matrix(adj.list, weights.list)
#' @export
convert.adjacency.list.to.adjacency.matrix <- function(adj.list, weights.list = NULL) {

  n <- length(adj.list)
  adj.matrix <- matrix(0, nrow = n, ncol = n)

  for (i in seq_along(adj.list)) {
    for (j in seq_along(adj.list[[i]])) {
      node <- adj.list[[i]][j]
      if (!is.null(weights.list)) {
        weight <- weights.list[[i]][j]
        adj.matrix[i, node] <- weight
      } else {
        adj.matrix[i, node] <- 1
      }
    }
  }

  adj.matrix
}

#' Converts a graph adjacency matrix to an adjacency list
#'
#' This function converts a given graph adjacency matrix into an adjacency list,
#' where each entry in the list corresponds to a node and contains the indices
#' of nodes it is connected to.
#'
#' @param M A matrix representing the adjacency matrix of a graph.
#' @return A list representing the adjacency list of the graph.
#' @export
convert.adjacency.matrix.to.adjacency.list <- function(M) {

    if (!is.matrix(M)) {
        stop("M has to be a matrix.")
    }

    L <- list()
    for (i in seq(nrow(M))) {
        L[[i]] <- which(M[i,] != 0)
    }

    L
}

#' Converts a graph weighted adjacency matrix to an adjacency list and weights list
#'
#' This function converts a given weighted graph adjacency matrix into two lists:
#' 1. An adjacency list where each entry corresponds to a node and contains the indices
#'    of nodes it is connected to.
#' 2. A weights list where each entry corresponds to a node and contains the weights
#'    of the edges connecting it to other nodes.
#'
#' @param M A matrix representing the weighted adjacency matrix of a graph.
#' @return A list with two elements: `adjacency.list` and `weights.list`.
#'         `adjacency.list` is a list where each element is a vector of connected node indices.
#'         `weights.list` is a list where each element is a vector of weights corresponding to the edges in the adjacency list.
#' @export
convert.weighted.adjacency.matrix.to.adjacency.list <- function(M) {

    if (!is.matrix(M)) {
        stop("M has to be a matrix.")
    }

    adjacency.list <- list()
    weights.list <- list()
    for (i in seq(nrow(M))) {
        connections <- which(M[i,] != 0)
        adjacency.list[[i]] <- connections
        weights.list[[i]] <- M[i, connections]
    }

    list(adjacency.list = adjacency.list, weights.list = weights.list)
}


#' Convert an Adjacency List to an Edgelist
#'
#' This function takes an adjacency list representation of a graph and converts it
#' into an edgelist format suitable for use with plotting functions like `igraph`.
#'
#' @param adj.list An adjacency list where each element is a vector of indices
#'   representing the neighbors of a node.
#' @return A list where each element is a vector of length two, representing
#'   an edge in the graph (source node, target node).
#' @export
convert.adjacency.to.edgelist <- function(adj.list) {
    edges <- list()
    ## Initialize empty list to store edges
    j <- 1
    for (i in seq_along(adj.list)) {
        source.node <- i
        for (target.node in adj.list[[i]]) {
            edges[[j]] <- c(source.node, target.node)
            j <- j + 1
        }
    }
    return(edges)
}

#' Convert Adjacency List to Edge Matrix
#'
#' This function converts an adjacency list representation of a graph to an edge matrix.
#' It can handle both weighted and unweighted graphs.
#'
#' @param adj.list A list where each element is a vector of integers representing
#'   the nodes adjacent to the node indexed by the element's position.
#' @param weights.list An optional list of numeric vectors representing the weights
#'   of the edges. If provided, each element should correspond to the weights of
#'   the edges in the same position in adj.list.
#'
#' @return A list with two elements:
#'   \item{edge.matrix}{A matrix where each row represents an edge, given as a pair of node indices.}
#'   \item{weights}{A numeric vector of weights corresponding to the edges in edge.matrix.
#'   This is NULL if weights.list was not provided.}
#'
#' @details The function converts the input adjacency list to 0-based indexing before
#'   passing it to the C++ function. The returned edge matrix uses 1-based indexing.
#'   For undirected graphs, each edge is included only once, with the lower index always
#'   appearing first.
#'
#' @examples
#' adj.list <- list(c(2,3), c(1,3), c(1,2))
#' result <- convert.adjacency.to.edge.matrix(adj.list)
#'
#' # With weights
#' weights.list <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
#' result_weighted <- convert.adjacency.to.edge.matrix(adj.list, weights.list)
#'
#' @export
convert.adjacency.to.edge.matrix <- function(adj.list, weights.list = NULL) {

    if (!is.null(weights.list) && !all(mapply(function(g, e) length(g) == length(e), adj.list, weights.list))) {
        stop("Each element of weights.list must have the same length as the corresponding element of adj.list.")
    }

    adj.list <- lapply(adj.list, function(x) as.integer(x - 1))
    return(.Call("S_convert_adjacency_to_edge_matrix", adj.list, weights.list, PACKAGE = "gflow"))
}

#' Converts an Edge Label List to an Edge Vector
#'
#' This function takes a list of edge labels and turns it into a vector.
#'
#' @param edge.label.list A list of edge labels.
#' @param rm.duplicates Set to TRUE to allow for duplicate edges.
#' @return A vector of edge labels.
convert.edge.label.list.to.edge.label.vector <- function(edge.label.list, rm.duplicates = TRUE) {
    edge.label <- c()
    for (i in seq_along(edge.label.list)) {
        source.node <- i
        for (label in edge.label.list[[i]]) {
            edge.label <- c(edge.label, label)
        }
    }

    if (rm.duplicates) {
        edge.label <- unique(edge.label)
    }

    edge.label
}


#' Converts a directed graph to an undirected graph
#'
#' This function takes the adjacency list representation of a directed graph and
#' converts it into the adjacency list representation of the corresponding undirected graph.
#'
#' @param adj.list A named list representing the adjacency list of the directed graph.
#'                 Each element of the list is a character vector containing the names
#'                 of the adjacent vertices for a given vertex.
#'
#' @return A named list representing the adjacency list of the corresponding undirected graph.
#'         Each element of the list is a character vector containing the names of the adjacent
#'         vertices for a given vertex in the undirected graph.
#'
#' @examples
#' directed_graph <- list(
#'   "A" = c("B", "C"),
#'   "B" = c("C"),
#'   "C" = c("D"),
#'   "D" = c("A")
#' )
#'
#' undirected_graph <- convert.to.undirected(directed_graph)
#' print(undirected_graph)
#'
#' @export
convert.to.undirected <- function(adj.list) {
    ## Create a new adjacency list to store the undirected graph
    undirected.adj.list <- list()

    if (is.null(names(adj.list))) {
        names(adj.list) <- as.character(seq_along(adj.list))
    }

    ## Iterate over each vertex and its adjacent vertices in the directed graph
    for (vertex in seq_along(adj.list)) {
        vertex.name <- names(adj.list)[vertex]
        for (neighbor in adj.list[[vertex]]) {
            ## Add the edge (vertex, neighbor) to the undirected graph
            if (!(vertex.name %in% names(undirected.adj.list))) {
                undirected.adj.list[[vertex.name]] <- c()
            }

            neighbor.name <- ifelse(neighbor %in% seq_along(adj.list), names(adj.list)[neighbor], as.character(neighbor))
            undirected.adj.list[[vertex.name]] <- c(undirected.adj.list[[vertex.name]], neighbor)

            ## Add the reverse edge (neighbor, vertex) to the undirected graph
            if (!(neighbor.name %in% names(undirected.adj.list))) {
                undirected.adj.list[[neighbor.name]] <- c()
            }
            undirected.adj.list[[neighbor.name]] <- c(undirected.adj.list[[neighbor.name]], vertex)
        }
    }

    ## Remove duplicates from the adjacency list of the undirected graph
    undirected.adj.list <- lapply(undirected.adj.list, unique)

    return(undirected.adj.list)
}

#' Removes self-loops from a graph
#'
#' This function takes an adjacency list representation of a graph and removes all self-loops
#' from it. A self-loop is an edge that connects a vertex to itself.
#'
#' @param adj.list A named list representing the adjacency list of the graph. Each element of the
#'                 list is a numeric vector containing the indices of the adjacent vertices for a
#'                 given vertex. The names of the list elements correspond to the vertex names.
#'
#' @return A named list representing the adjacency list of the graph with self-loops removed.
#'         The format is the same as the input adjacency list, but without any self-loops.
#'
#' @examples
#' graph <- list(
#'   "A" = c(1, 2),
#'   "B" = c(2, 3),
#'   "C" = c(1, 3),
#'   "D" = c(4)
#' )
#'
#' graph.no.self.loops <- rm.self.loops(graph)
#' print(graph.no.self.loops)
#'
#' @export
rm.self.loops <- function(adj.list) {
    ## Iterate over each vertex and its adjacent vertices in the graph
    for (vertex in seq_along(adj.list)) {
        ## Remove self-loops by filtering out the vertex itself from its adjacency list
        adj.list[[vertex]] <- adj.list[[vertex]][adj.list[[vertex]] != vertex]
    }

    return(adj.list)
}

## -------------------------------------------------------------------------------------------
##
## Distance between graphs
##
## -------------------------------------------------------------------------------------------

#' Hungarian-Frobenius Graph Matching Algorithm
#'
#' This function implements the Hungarian-Frobenius Graph Matching Algorithm to measure the similarity
#' between two graphs. It uses the Hungarian algorithm to find the optimal permutation of vertices that
#' minimizes the cost of matching the vertices between the graphs, and then calculates the Frobenius norm
#' distance between the distance matrices of the graphs after applying the optimal permutation.
#'
#' @param graph1 A list representing the adjacency list of the first graph.
#' @param graph2 A list representing the adjacency list of the second graph.
#' @return The similarity score between the two graphs, ranging from 0 to 1, where 0 indicates perfect similarity.
#' @export
hungarian.frobenius.graph.matching <- function(graph1, graph2) {
  g1.m <- convert.adjacency.to.edge.matrix(graph1)$edge.matrix
  g2.m <- convert.adjacency.to.edge.matrix(graph2)$edge.matrix

  g1 <- igraph::graph_from_edgelist(g1.m, directed = FALSE)
  g2 <- igraph::graph_from_edgelist(g2.m, directed = FALSE)

  igraph::V(g1)$label <- as.character(1:igraph::vcount(g1))
  igraph::V(g2)$label <- as.character(1:igraph::vcount(g2))

  ## Calculate the minimal path distance matrices
  D1 <- igraph::distances(g1)
  D2 <- igraph::distances(g2)

  if (igraph::vcount(g1) == igraph::vcount(g2)) {
    ## Compute the cost matrix
    C <- matrix(0, nrow = igraph::vcount(g1), ncol = igraph::vcount(g2))
    for (i in 1:igraph::vcount(g1)) {
      for (j in 1:igraph::vcount(g2)) {
        C[i, j] <- sum((D1[i, ] - D2[j, ])^2)
      }
    }

    ## Find the optimal permutation using the Hungarian algorithm
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' is required for this function. Please install it with install.packages('clue')")
    }
    mapping <- clue::solve_LSAP(C)

    ## Permute the vertices of graph 2 according to the optimal mapping
    perm_g2 <- igraph::permute(g2, mapping)

    ## Recompute the distance matrix for the permuted graph 2
    perm_D2 <- igraph::distances(perm_g2)  # Corrected: Changed shortest.paths to distances

    ## Calculate the Frobenius distance between the original and permuted distance matrices
    deviation <- norm(D1 - perm_D2, type = "F")
  } else {
    ## Determine the maximum number of vertices between the two graphs
    max_vertices <- max(igraph::vcount(g1), igraph::vcount(g2))

    ## Pad the distance matrices with dummy vertices
    padded_D1 <- matrix(max_vertices, nrow = max_vertices, ncol = max_vertices)
    padded_D1[1:igraph::vcount(g1), 1:igraph::vcount(g1)] <- D1

    padded_D2 <- matrix(max_vertices, nrow = max_vertices, ncol = max_vertices)
    padded_D2[1:igraph::vcount(g2), 1:igraph::vcount(g2)] <- D2

    ## Compute the cost matrix
    C <- matrix(0, nrow = max_vertices, ncol = max_vertices)
    for (i in 1:max_vertices) {
      for (j in 1:max_vertices) {
        C[i, j] <- sum((padded_D1[i, ] - padded_D2[j, ])^2)
      }
    }

    ## Find the optimal permutation using the Hungarian algorithm
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' is required for this function. Please install it with install.packages('clue')")
    }
    mapping <- clue::solve_LSAP(C, maximum = FALSE)

      S <- intersect(mapping[1:igraph::vcount(g2)], seq(igraph::vcount(g2)))
      valid_mapping <- c(S, setdiff(seq(igraph::vcount(g2)), S))

    ## Ensure valid_mapping is within the correct range
    if (length(valid_mapping) != igraph::vcount(g2)) {
      stop("Invalid mapping length. length(valid_mapping):", length(valid_mapping), " vcount(g2): ", igraph::vcount(g2))
    }

    ## Permute the vertices of graph 2 according to the optimal mapping
    perm_g2 <- igraph::permute(g2, valid_mapping)

    ## Recompute the distance matrix for the permuted graph 2
    perm_D2 <- igraph::distances(perm_g2)  # Corrected: Changed shortest.paths to distances

    ## Calculate the Frobenius distance between the original and permuted distance matrices
    deviation <- norm(padded_D1 - padded_D2, type = "F")
  }

  ## Normalize the deviation
  n <- vcount(g1)
  max_distance <- n - 1
  max_deviation <- sqrt(n * max_distance^2)
  normalized_deviation <- deviation / max_deviation

  return(normalized_deviation)
}


#' Weighted Graph Matching between two graphs with identical vertex sets
#'
#' This function computes a similarity measure between two weighted graphs
#' derived from the same dataset. It assumes that both graphs have identical
#' vertex sets and indexing, eliminating the need for vertex permutation.
#'
#' @param graph1.adj.list A list representing the adjacency list of the first graph.
#' @param graph1.weights A list representing the weights of the first graph.
#' @param graph2.adj.list A list representing the adjacency list of the second graph.
#' @param graph2.weights A list representing the weights of the second graph.
#' @param calculate.normalized.deviation Logical indicating whether to calculate the
#'        normalized deviation (default: FALSE). If TRUE, returns normalized deviation;
#'        if FALSE, returns raw distance difference.
#'
#' @return A numeric value between 0 and 1 representing the normalized
#'   deviation between the two graphs. A value closer to 0 indicates higher
#'   similarity between the graphs.
#'
#' @details
#' The function performs the following steps:
#' 1. Converts the input lists to igraph objects.
#' 2. Calculates the distance matrices for both graphs using edge weights.
#' 3. Computes the L1 norm of the difference between these distance matrices.
#' 4. Normalizes the result using the maximum possible deviation, which is
#'    based on the largest distance found in either graph.
#'
#' @note
#' - This function assumes that the graphs are undirected.
#' - The weights should be stored in the input lists as a vector named 'weight'.
#' - For disconnected graphs, infinite distances are ignored in calculating
#'   the maximum distance.
#'
#' @examples
#' # Create adjacency lists (as lists, not matrices)
#' g1.adj.list <- list(c(2), c(3), c(1))  # vertex 1->2, vertex 2->3, vertex 3->1
#' g1.weights <- list(c(1), c(2), c(3))   # corresponding weights
#'
#' g2.adj.list <- list(c(2), c(3), c(1))  # same structure
#' g2.weights <- list(c(1), c(3), c(2))   # different weights
#'
#' similarity <- identical.vertex.set.weighted.graph.similarity(
#'   g1.adj.list, g1.weights, g2.adj.list, g2.weights
#' )
#' print(similarity)
#'
#' @importFrom igraph graph_from_edgelist E distances vcount
#'
#' @export
identical.vertex.set.weighted.graph.similarity <- function(graph1.adj.list,
                                                           graph1.weights,
                                                           graph2.adj.list,
                                                           graph2.weights,
                                                           calculate.normalized.deviation = FALSE) {

    graph1.obj <- convert.adjacency.to.edge.matrix(graph1.adj.list, graph1.weights)
    graph2.obj <- convert.adjacency.to.edge.matrix(graph2.adj.list, graph2.weights)

    graph1 <- igraph::graph_from_edgelist(graph1.obj$edge.matrix, directed = FALSE)
    graph2 <- igraph::graph_from_edgelist(graph2.obj$edge.matrix, directed = FALSE)

    igraph::E(graph1)$weight <- graph1.obj$weights
    igraph::E(graph2)$weight <- graph2.obj$weights

    ## Calculate the minimal path distance matrices
    D1 <- igraph::distances(graph1, weights = igraph::E(graph1)$weight)
    D2 <- igraph::distances(graph2, weights = igraph::E(graph2)$weight)

    ## Calculate the difference using L1 norm
    deviation <- sum(abs(D1 - D2))

    if (calculate.normalized.deviation) {
        ## Find the maximum distance in each graph
        max_distance_graph1 <- max(D1[is.finite(D1)])
        max_distance_graph2 <- max(D2[is.finite(D2)])

        ## Use the larger of the two maximum distances
        max_distance <- max(max_distance_graph1, max_distance_graph2)

        ## Normalize the deviation
        n <- igraph::vcount(graph1)
        max_deviation <- n * n * max_distance  ## Maximum possible L1 difference
        deviation <- deviation / max_deviation
    }

    return(deviation)
}

#' Assign Vertices to Connected Components
#'
#' This function takes a graph representation and assigns each vertex to a connected component.
#'
#' @param adj.list A list representing the adjacency list of the graph.
#'        Each element of the list corresponds to a vertex and contains the indices of its neighbors.
#'
#' @return An integer vector where each element represents the connected component
#'         ID for the corresponding vertex.
#'
#' @examples
#' adj.list <- list(c(2,3), c(1), c(1), c(5), c(4))
#' components <- graph.connected.components(adj.list)
#' print(components)
#'
#' @export
graph.connected.components <- function(adj.list) {

    ## Check if input is a list
    if (!is.list(adj.list)) {
        stop("Input must be a list representing the graph's adjacency list.")
    }

    ## Check if the list is empty
    if (length(adj.list) == 0) {
        stop("Input graph is empty.")
    }

    ## Check each element of the list
    for (i in seq_along(adj.list)) {
        ## Check if each element is a numeric vector
        if (!is.numeric(adj.list[[i]])) {
            stop(paste("Element", i, "of the input list is not a numeric vector."))
        }

        ## Check if each element contains only positive integers
        if (any(adj.list[[i]] != as.integer(adj.list[[i]])) || any(adj.list[[i]] <= 0)) {
            stop(paste("Element", i, "of the input list contains non-positive or non-integer values."))
        }

        ## Check if each element references valid vertices
        if (any(adj.list[[i]] > length(adj.list))) {
            stop(paste("Element", i, "of the input list contains invalid vertex indices."))
        }
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call("S_graph_connected_components",
                    adj.list.0based)

    ## Check the result
    if (!is.integer(result) || length(result) != length(adj.list)) {
        stop("Internal error: Invalid result from C++ function.")
    }

    return(result)
}


#' Compute Edge Difference Between Two Graphs
#'
#' This function takes two graphs represented as adjacency lists and computes
#' the edge difference from graph1 to graph2. It returns a list where each
#' component contains the edges present in graph1 but not in graph2 for each vertex.
#'
#' @param graph1 A list representing the adjacency list of the first graph.
#'        Each element of the list corresponds to a vertex and contains the indices of its neighbors.
#' @param graph2 A list representing the adjacency list of the second graph.
#'        It must have the same number of vertices as graph1.
#'
#' @return A list of the same length as the input graphs. Each element contains
#'         the neighbors of the corresponding vertex that are present in graph1 but not in graph2.
#'
#' @examples
#' graph1 <- list(c(2,3), c(1,3), c(1,2), c(5), c(4))
#' graph2 <- list(c(2), c(1), c(), c(5), c(4))
#' diff <- edge.diff(graph1, graph2)
#' print(diff)
#'
#' @export
edge.diff <- function(graph1, graph2) {

    ## Check if inputs are lists
    if (!is.list(graph1) || !is.list(graph2)) {
        stop("Both inputs must be lists representing graph adjacency lists.")
    }

    ## Check if the graphs have the same number of vertices
    if (length(graph1) != length(graph2)) {
        stop("The two graphs must have the same number of vertices.")
    }

    ## Initialize the result list
    result <- vector("list", length(graph1))

    ## Compute the edge difference for each vertex
    for (i in seq_along(graph1)) {
        result[[i]] <- setdiff(graph1[[i]], graph2[[i]])
    }

    return(result)
}

#' Create a subgraph from a given graph
#'
#' This function creates a subgraph from a given graph based on specified indices or IDs.
#'
#' @param S.graph A list containing the original graph structure with adjacency and distance lists.
#' @param id.indices A vector of indices to include in the subgraph. Default is NULL.
#' @param ids A vector of IDs to include in the subgraph. Default is NULL.
#' @param S A data frame or matrix where rownames correspond to node IDs. Required if `ids` is provided. Default is NULL.
#' @param use.sequential.indices Logical. If TRUE, renumber the indices in the subgraph to be 1:length(id.indices). Default is FALSE.
#'
#' @return A list containing the subgraph structure with adjacency and distance lists.
#'
#' @details
#' The function can create a subgraph based on either `id.indices` or `ids`. If `ids` is provided,
#' `S` must also be provided to map the IDs to indices. The `use.sequential.indices` parameter
#' allows for renumbering the indices in the subgraph to be sequential, which can be useful for
#' certain applications.
#'
#' @examples
#' X <- runif.sphere(20, 2)
#' graph <- create.single.iknn.graph(X, k = 3, compute.full = TRUE, verbose = FALSE)
#' graph$dist_list <- graph$weight_list
#' graph$weight_list <- NULL
#' subgraph <- create.subgraph(graph, id.indices = c(1:10))
#' subgraph_sequential <- create.subgraph(graph, id.indices = c(1:10, 15:16),
#'                                        use.sequential.indices = TRUE)
#' @export
create.subgraph <- function(S.graph, id.indices = NULL, ids = NULL, S = NULL, use.sequential.indices = FALSE) {
    ## Parameter checks
    if (!is.list(S.graph) || !all(c("adj_list", "dist_list") %in% names(S.graph))) {
        stop("S.graph must be a list containing 'adj_list' and 'dist_list'")
    }
    if (is.null(id.indices) && is.null(ids)) {
        stop("Either id.indices or ids must be provided")
    }
    if (!is.null(ids) && is.null(S)) {
        stop("If ids are provided, S must also be provided")
    }
    if (!is.null(S) && !is.null(ids)) {
        if (!all(ids %in% rownames(S))) {
            stop("All ids must be present in rownames(S)")
        }
    }
    if (!is.null(id.indices) && !all(id.indices %in% seq_along(S.graph$adj_list))) {
        stop("All id.indices must be valid indices in S.graph")
    }
    if (!is.logical(use.sequential.indices)) {
        stop("use.sequential.indices must be a logical value (TRUE or FALSE)")
    }

    if (!is.null(ids) && !is.null(S)) {
        ## Get the indices of the ids in S
        id.indices <- match(ids, rownames(S))
    }

    ## Initialize the subgraph object
    S.subgraph <- list(adj_list = list(), dist_list = list())

    ## Create a mapping from original indices to sequential indices if needed
    if (use.sequential.indices) {
        index.map <- setNames(seq_along(id.indices), id.indices)
    }

    ## For each id in the subgraph
    for (i in seq_along(id.indices)) {
        ## Get the original index
        orig.index <- id.indices[i]

        ## Get the adjacent nodes and distances
        adj.nodes <- S.graph$adj_list[[orig.index]]
        dist.nodes <- S.graph$dist_list[[orig.index]]

        ## Handle NULL or empty adjacency lists
        if (is.null(adj.nodes) || length(adj.nodes) == 0) {
            S.subgraph$adj_list[[i]] <- integer(0)
            S.subgraph$dist_list[[i]] <- numeric(0)
            next
        }

        ## Find which of these nodes are in the subgraph
        in.subgraph <- adj.nodes %in% id.indices

        ## Add to the subgraph
        if (use.sequential.indices) {
            # Map the adjacent nodes to their new sequential indices
            mapped_indices <- index.map[as.character(adj.nodes[in.subgraph])]
            S.subgraph$adj_list[[i]] <- as.integer(mapped_indices)
        } else {
            S.subgraph$adj_list[[i]] <- adj.nodes[in.subgraph]
        }
        S.subgraph$dist_list[[i]] <- dist.nodes[in.subgraph]
    }

    ## Set appropriate names for the lists
    if (use.sequential.indices) {
        names(S.subgraph$adj_list) <- seq_along(id.indices)
        names(S.subgraph$dist_list) <- seq_along(id.indices)
    } else {
        names(S.subgraph$adj_list) <- id.indices
        names(S.subgraph$dist_list) <- id.indices
    }

    return(S.subgraph)
}


#' Calculate Degree Distribution Properties for Random Points on a Sphere
#'
#' @description
#' Simulates points uniformly on a sphere and computes the degree distribution
#' properties of their k-nearest neighbor graph. For each simulation, points are
#' generated, a k-NN graph is constructed, and the proportion of vertices with
#' each degree is calculated. The function returns mean proportions and confidence
#' intervals across all simulations.
#'
#' @param n.pts numeric; Number of points to generate on the sphere for each simulation.
#' @param n.sims numeric; Number of simulations to run.
#' @param k numeric; Number of nearest neighbors to use in constructing the graph.
#' @param dim numeric; Dimension of the sphere (e.g., 2 for circle, 3 for sphere).
#'
#' @return A list containing:
#' \describe{
#'   \item{mean.props}{Vector of mean proportions for each degree}
#'   \item{ci.lower}{Vector of lower 95% confidence interval bounds}
#'   \item{ci.upper}{Vector of upper 95% confidence interval bounds}
#'   \item{degrees}{Vector of degree values corresponding to the proportions}
#' }
#'
#' @details
#' The function generates uniform random points on a sphere using \code{runif.sphere}
#' and constructs k-nearest neighbor graphs using \code{create.single.iknn.graph}. For each
#' simulation, it computes the proportion of vertices with each degree. The final
#' results include means and 95% confidence intervals for these proportions across
#' all simulations.
#'
#' The confidence intervals are computed using the normal approximation:
#' mean +/- 1.96 * (standard deviation / sqrt(n.sims))
#'
#' @examples
#' \dontrun{
#' # Calculate degree distribution properties for 1000 points on a circle
#' circle_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
#'
#' # Calculate for points on a sphere
#' sphere_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 3)
#' }
#' @export
get.sphere.degree.props <- function(n.pts = 1000, n.sims = 100, k = 10, dim = 2) {

    if(!all(sapply(list(n.pts, n.sims, k, dim), is.numeric))) {
        stop("All parameters must be numeric")
    }
    if(!all(sapply(list(n.pts, n.sims, k, dim), function(x) x > 0))) {
        stop("All parameters must be positive")
    }
    if(dim < 1) stop("dimension must be at least 1")
    if(k >= n.pts) stop("k must be less than n.pts")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        X <- runif.sphere(n.pts, dim = dim)
        graph <- create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)  # Added to show which degrees correspond to the proportions
}


#' Generate Degree Distribution Properties in Tubular Neighborhood of Unit Circle
#'
#' @description
#' Generates random samples in a tubular neighborhood of the unit circle and computes
#' degree distribution properties of the resulting k-nearest neighbor graph. The function
#' performs multiple simulations to estimate mean proportions and confidence intervals
#' for each degree.
#'
#' @param n.pts Positive integer. Number of points to generate in each simulation.
#' @param n.sims Positive integer. Number of simulations to run.
#' @param k Positive integer. Number of nearest neighbors for graph construction. Must be less than n.pts.
#' @param noise Non-negative numeric. Standard deviation of the noise added to the radius.
#' @param noise.type Character string. Type of noise distribution to use ("laplace" or "normal").
#'
#' @return A list containing:
#' \itemize{
#'   \item mean.props: Vector of mean proportions for each degree
#'   \item ci.lower: Vector of lower 95% confidence interval bounds
#'   \item ci.upper: Vector of upper 95% confidence interval bounds
#'   \item degrees: Vector of degrees corresponding to the proportions
#' }
#'
#' @details
#' The function uses `generate.circle.data()` to create points and `create.single.iknn.graph()`
#' to construct the k-nearest neighbor graph. It computes degree distributions for each
#' simulation and aggregates the results to estimate population parameters.
#'
#' @examples
#' # Generate degree distribution properties with default parameters
#' results <- get.TN.S1.degree.props()
#'
#' # Use custom parameters
#' results <- get.TN.S1.degree.props(
#'   n.pts = 500,
#'   n.sims = 50,
#'   k = 5,
#'   noise = 0.05,
#'   noise.type = "normal"
#' )
#'
#' @importFrom stats sd
#' @export
get.TN.S1.degree.props <- function(n.pts = 1000,
                                   n.sims = 100,
                                   k = 10,
                                   noise = 0.1,
                                   noise.type = "laplace") {

    noise.type <- match.arg(noise.type, c("laplace", "normal"))

    if (!all(sapply(list(n.pts, n.sims, k, noise), is.numeric))) {
        stop("The first four parameters must be numeric")
    }
    if (!all(sapply(list(n.pts, n.sims, k), function(x) x > 0))) {
        stop("The first three parameters must be positive")
    }
    if (k >= n.pts) stop("k must be less than n.pts")
    if (noise < 0) stop("noise has to be non-negative")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        cX <- generate.circle.data(n.pts, radius = 1, noise = noise, type = "random", noise.type = noise.type)
        X <- cX[,1:2]
        graph <- create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)
}

#' Calculate Degree Distribution Properties for Random Points on a Torus
#'
#' @description
#' Simulates points uniformly on a torus and computes the degree distribution
#' properties of their k-nearest neighbor graph. For each simulation, points are
#' generated, a k-NN graph is constructed, and the proportion of vertices with
#' each degree is calculated. The function returns mean proportions and confidence
#' intervals across all simulations.
#'
#' @param n.pts numeric; Number of points to generate on the torus for each simulation.
#' @param n.sims numeric; Number of simulations to run.
#' @param k numeric; Number of nearest neighbors to use in constructing the graph.
#' @param dim numeric; Dimension of the torus (e.g., 1 for circle, 2 for 2-torus).
#'
#' @return A list containing:
#' \describe{
#'   \item{mean.props}{Vector of mean proportions for each degree}
#'   \item{ci.lower}{Vector of lower 95% confidence interval bounds}
#'   \item{ci.upper}{Vector of upper 95% confidence interval bounds}
#'   \item{degrees}{Vector of degree values corresponding to the proportions}
#' }
#'
#' @details
#' The function generates uniform random points on a torus using \code{\link{runif.torus}}
#' and constructs k-nearest neighbor graphs using \code{create.single.iknn.graph}. For each
#' simulation, it computes the proportion of vertices with each degree. The final
#' results include means and 95% confidence intervals for these proportions across
#' all simulations.
#'
#' The confidence intervals are computed using the normal approximation:
#' \eqn{\text{mean} \pm 1.96 \times (\text{standard deviation} / \sqrt{n.sims})}
#'
#' @examples
#' \dontrun{
#' # Calculate degree distribution properties for 1000 points on a circle
#' circle_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 1)
#'
#' # Calculate for points on a 2-torus
#' torus_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
#' }
#' @export
get.torus.degree.props <- function(n.pts = 1000, n.sims = 100, k = 10, dim = 1) {
    if(!all(sapply(list(n.pts, n.sims, k, dim), is.numeric))) {
        stop("All parameters must be numeric")
    }
    if(!all(sapply(list(n.pts, n.sims, k, dim), function(x) x > 0))) {
        stop("All parameters must be positive")
    }
    if(round(dim) != dim) {
        stop("dimension must be a positive integer")
    }
    if(k >= n.pts) stop("k must be less than n.pts")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        X <- runif.torus(n.pts, dim = dim)
        graph <- create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd, na.rm = TRUE)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)
}

#' Create Labels for Vertices in State Space
#'
#' @description
#' Creates unique labels for specified vertices in state space by concatenating
#' two-letter shortcuts of the most abundant taxa. The function ensures uniqueness
#' of labels by iteratively including additional taxa if needed.
#'
#' @param vertices Vector of indices corresponding to rows in state.space
#' @param state.space Matrix of ASV counts/abundances. When taxonomy is NULL,
#'        column names should contain species names (e.g., "Genus_species")
#' @param taxonomy Taxonomy information for ASVs (default: NULL). If NULL,
#'        species names are taken directly from state.space column names
#' @param min.relAb.thld Minimum relative abundance threshold for including taxa in labels (default: 0.05)
#' @param profile.length Integer specifying the number of top species to keep in each profile (default: 5)
#'
#' @return A list with two components:
#'   \item{labels}{Named character vector of vertex labels}
#'   \item{profiles}{Named list of relative abundance profiles for vertices, where names are point indices and each profile is a matrix containing the top profile.length species}
#'
#' @details
#' The function processes the specified vertices, creating unique labels by:
#' 1. For each vertex, identifying taxa above the relative abundance threshold
#' 2. Creating initial labels using two-letter shortcuts
#' 3. Ensuring uniqueness by incorporating additional taxa if needed
#'
#' When taxonomy is NULL, the function uses column names of state.space directly as species names.
#' The profile.length parameter controls how many species are kept in each profile, keeping the most abundant ones.
#'
#' @examples
#' \dontrun{
#' # Keep only top 5 species in profiles
#' result <- create.vertex.labels(c(1,3,5), state.space, taxonomy,
#'                               min.relAb.thld = 0.05, profile.length = 5)
#' }
#' @export
create.vertex.labels <- function(vertices, state.space, taxonomy = NULL, min.relAb.thld = 0.05, profile.length = 5) {
    # Helper function to create two-letter shortcuts from species names
    create.shortcut <- function(sp.name) {
        parts <- strsplit(sp.name, "_")[[1]]
        if (length(parts) >= 2) {
            paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))
        } else {
            paste0(substr(parts[1], 1, 1), "")
        }
    }

    # Helper function to get profile, handling both taxonomy and direct column names cases
    get.profile <- function(id, state.space, taxonomy) {
        if (is.null(taxonomy)) {
            # Use state.space directly
            abundances <- state.space[id, ]
            profile <- cbind(
                species = colnames(state.space),
                abundance = as.numeric(abundances)
            )
            # Sort by abundance in decreasing order
            profile <- profile[order(as.numeric(profile[, 2]), decreasing = TRUE), , drop = FALSE]
        } else {
            # Use prof.fn when taxonomy is provided
            profile <- prof.fn(id, state.space, taxonomy)
        }

        # Limit profile length
        if (nrow(profile) > profile.length) {
            profile <- profile[1:profile.length, , drop = FALSE]
        }

        return(profile)
    }

    # Initialize labels vector and profiles list
    labels <- character(length(vertices))
    names(labels) <- vertices
    profiles <- list()

    # First pass: create initial labels and store profiles
    for (i in seq_along(vertices)) {
        point.i <- vertices[i]
        id <- if (is.null(taxonomy)) point.i else rownames(state.space)[point.i]

        # Get profile for this point
        profile <- get.profile(id, state.space, taxonomy)
        # Store profile
        profiles[[as.character(point.i)]] <- profile

        abundances <- as.numeric(profile[, 2])
        above.threshold <- which(abundances >= min.relAb.thld)

        if (length(above.threshold) > 0) {
            taxa.names <- profile[above.threshold, 1]
            shortcuts <- sapply(taxa.names, create.shortcut)
            labels[i] <- paste(shortcuts, collapse = "")
        } else {
            labels[i] <- create.shortcut(profile[1, 1])
        }
    }

    # Second pass: ensure uniqueness
    while (any(duplicated(labels))) {
        dups <- labels[duplicated(labels)]
        for (dup.label in unique(dups)) {
            dup.indices <- which(labels == dup.label)
            for (idx in dup.indices) {
                point.i <- vertices[idx]
                profile <- profiles[[as.character(point.i)]]

                current.taxa <- ceiling(nchar(labels[idx])/2)
                if (current.taxa < nrow(profile)) {
                    new.taxon <- create.shortcut(profile[current.taxa + 1, 1])
                    labels[idx] <- paste0(labels[idx], new.taxon)
                }
            }
        }
    }

    return(list(
        labels = labels,
        profiles = profiles
    ))
}


#' Create Label Tables and Indicator Vectors for Vertices
#'
#' @description
#' Creates a label table and indicator vector for specified vertices in the state space.
#' The label table maps point IDs to their corresponding labels, while the indicator
#' vector marks the presence/absence of vertices in the state space.
#'
#' @param vertex.labels Named character vector of labels for vertices
#' @param state.space State space matrix where rownames correspond to point IDs
#'
#' @return A list with two components:
#'   \item{lab.tbl}{Named character vector mapping point IDs to vertex labels}
#'   \item{ind}{Numeric vector indicating vertices (1) and non-vertices (0)}
#'
#' @details
#' The function processes vertex labels to create:
#' 1. A label table that maps point IDs to their corresponding labels
#' 2. A binary indicator vector marking the presence (1) or absence (0) of vertices
#' The vector maintains the same length as the number of rows in state.space and uses consistent naming.
#'
#' @examples
#' \dontrun{
#' # Given vertex.labels and state space matrix state.space
#' result <- create.vertex.label.indicators(vertex.labels, state.space)
#' print(head(result$lab.tbl))
#' print(sum(result$ind)) # Number of vertices
#' }
#' @export
create.vertex.label.indicators <- function(vertex.labels, state.space) {
    # Convert label names to integers
    point.indices <- as.integer(names(vertex.labels))

    # Get corresponding row names from state.space
    point.ids <- rownames(state.space)[point.indices]

    # Create label table
    lab.tbl <- c()
    lab.tbl[point.ids] <- vertex.labels[as.character(point.indices)]

    # Create indicator vector
    ind <- numeric(nrow(state.space))
    names(ind) <- rownames(state.space)
    ind[point.ids] <- 1

    # Return results
    return(list(
        lab.tbl = lab.tbl,
        ind = ind
    ))
}


#' Calculate Shortest Path Distance Between Two Vertices in a Graph
#'
#' @description
#' Computes the shortest path distance between two vertices in any weighted graph,
#' using Dijkstra's algorithm.
#'
#' @param star.obj Output from generate.star.dataset(). If NULL, adj.list and edge.lengths must be provided.
#' @param i Index of source vertex
#' @param j Index of target vertex
#' @param adj.list List where each element i contains indices of vertices adjacent to vertex i
#' @param edge.lengths List where element i contains lengths of edges from vertex i to its adjacent vertices
#'
#' @return Numeric value representing shortest path distance between vertices i and j
#'
#' @examples
#' \dontrun{
#' # Using star_obj
#' star.obj <- generate.star.dataset(n.points = 10, n.arms = 3)
#' dist <- compute.graph.distance(star.obj = star.obj, i = 1, j = 5)
#'
#' # Using adjacency list and edge lengths directly
#' # Assuming adj_list and edge_lengths are defined
#' dist <- compute.graph.distance(i = 1, j = 5,
#'                               adj.list = adj_list,
#'                               edge.lengths = edge_lengths)
#' }
#' @export
compute.graph.distance <- function(star.obj = NULL,
                                   i,
                                   j,
                                   adj.list = NULL,
                                   edge.lengths = NULL) {
    ## Algorithm outline:
    ## 1. Setup: Extract/validate graph structure
    ## 2. Initialize data structures for Dijkstra's algorithm
    ## 3. Main loop: Find shortest paths until target reached
    ## 4. Return final distance to target

    ## Step 1: Setup and validation
    ## Extract graph structure from star.obj if provided
    if (!is.null(star.obj)) {
        adj.list <- star.obj$adj.list      ## List of adjacent vertices for each vertex
        edge.lengths <- star.obj$edge.lengths  ## Corresponding edge lengths
    }

    ## Verify we have required graph information
    if (is.null(adj.list) || is.null(edge.lengths)) {
        stop("Must provide either star.obj or both adj.list and edge.lengths")
    }

    ## Handle trivial case - distance to self is 0
    if (i == j) return(0)

    ## Validate vertex indices
    n <- length(adj.list)  ## Total number of vertices
    if (i > n || j > n || i < 1 || j < 1) {
        stop("Invalid vertex indices")
    }

    ## Step 2: Initialize Dijkstra's algorithm data structures
    distances <- rep(Inf, n)    ## Track shortest known distance to each vertex
    distances[i] <- 0           ## Distance to start vertex is 0
    visited <- rep(FALSE, n)    ## Track which vertices have been processed

    ## Step 3: Main loop of Dijkstra's algorithm
    while (!visited[j] && !all(visited)) {  ## Continue until target found or all vertices visited
        ## Find unvisited vertex with smallest current distance
        curr <- which.min(ifelse(visited, Inf, distances)) ## set distance to visited vertices Inf, this gives the index of an unvisited vertex with the smallest distance
        if (distances[curr] == Inf) break  ## No path exists to remaining vertices

        visited[curr] <- TRUE  ## Mark current vertex as processed

        ## Update distances to all neighbors of current vertex
        for (k in seq_along(adj.list[[curr]])) {
            neighbor <- adj.list[[curr]][k]             ## Get index of neighboring vertex
            if (!visited[neighbor]) {                   ## Only update unvisited neighbors
                new_dist <- distances[curr] + edge.lengths[[curr]][k]  ## Calculate potential new distance
                distances[neighbor] <- min(distances[neighbor], new_dist)  ## Update if shorter path found
            }
        }
    }

    ## Step 4: Return final shortest distance to target vertex
    return(distances[j])
}


#' Count Edges in an Undirected Graph
#'
#' @description
#' Computes the total number of edges in an undirected graph represented by an adjacency list.
#' This function works by counting all neighbor connections across the graph and then
#' dividing by 2, since each edge is represented twice in an undirected graph (once
#' from each endpoint).
#'
#' @param adj.list A list where each element \code{adj.list[[i]]} contains the neighbors
#'   of vertex \code{i} as a vector of indices. For an undirected graph, if \code{j} is in
#'   \code{adj.list[[i]]}, then \code{i} should also be in \code{adj.list[[j]]}.
#'
#' @return A numeric value representing the number of edges in the graph.
#'   Returns 0 for an empty graph.
#'
#' @examples
#' # Create an adjacency list for a simple undirected graph
#' # with 3 vertices and 2 edges: (1-2) and (2-3)
#' adj <- list(
#'   c(2),    # Vertex 1 is connected to vertex 2
#'   c(1, 3), # Vertex 2 is connected to vertices 1 and 3
#'   c(2)     # Vertex 3 is connected to vertex 2
#' )
#' count.edges(adj)  # Should return 2
#'
#' @export
count.edges <- function(adj.list) {
    n.edges <- 0
    for (i in seq_along(adj.list)) {
        nbrs <- adj.list[[i]]
        n.edges <- n.edges + length(nbrs)
    }

    return(n.edges / 2)
}

#' Get Unique Edge Weights from a Weighted Graph in Parallel
#'
#' @description
#' Returns a vector containing the weights of all edges in a weighted graph
#' without duplication, using parallel processing for improved performance
#' on large graphs.
#'
#' @param adj.list A list where each element contains the adjacency list for a vertex.
#'        Element i contains a vector of vertices that are adjacent to vertex i.
#' @param weight.list A list with the same structure as \code{adj.list}, where \code{weight.list[[i]][j]}.
#'        contains the weight of the edge between vertex i and its j-th neighbor.
#' @param n.cores Number of cores to use for parallel processing. Default is 2.
#'
#' @return A numeric vector containing the weights of all unique edges in the graph.
#'
#' @details
#' The function parallelizes the processing by dividing the vertices among different cores.
#' It handles undirected graphs by only processing edges where i < neighbor to ensure
#' each edge is counted exactly once.
#' The function requires the foreach and doParallel packages.
#'
#' @examples
#' \dontrun{
#' # Create a simple undirected weighted graph
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight.list <- list(c(5, 10), c(5, 7), c(10, 7))
#'
#' # Get weights of all edges using 2 cores
#' edge_weights <- get.edge.weights(adj.list, weight.list, n.cores = safe_cores(1))
#' }
#'
#' @import foreach
#' @import doParallel
#' @export
get.edge.weights <- function(adj.list,
                             weight.list,
                             n.cores = 12) {

    ## Load required libraries
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE)) {
        stop("Packages 'foreach' and 'doParallel' are required for this function")
    }

    ## Set up parallel backend
    doParallel::registerDoParallel(cores = n.cores)

    ## Create chunks of vertices to process in parallel
    n.vertices <- length(adj.list)
    vertices.per.chunk <- ceiling(n.vertices / n.cores)
    vertex.chunks <- split(1:n.vertices,
                           ceiling(seq_along(1:n.vertices) / vertices.per.chunk))

    ## Process each chunk in parallel and combine results
    chunk <- NULL  # to avoid R CMD check NOTE about undefined global
    results <- foreach::foreach(chunk = vertex.chunks,
                                .combine = 'c',
                                .packages = c()) %dopar% {
                                    ## Initialize storage for this chunk's edge weights
                                    chunk.weights <- c()

                                    ## Process each vertex in the chunk
                                    for (i in chunk) {
                                        nbrs <- adj.list[[i]]

                                        for (j in seq_along(nbrs)) {
                                            neighbor <- nbrs[j]

                                            ## Only process edges where the current vertex has a smaller index
                                            ## than its neighbor to avoid duplication
                                            if (i < neighbor) {
                                                weight <- weight.list[[i]][j]
                                                chunk.weights <- c(chunk.weights, weight)
                                            }
                                        }
                                    }

                                    return(chunk.weights)
                                }

    ## Stop parallel backend
    doParallel::stopImplicitCluster()

    return(results)
}



#' Generate a Gaussian Mixture Function on a Graph
#'
#' @description
#' Creates a function on a graph that is a mixture of Gaussian-like components centered
#' at specified vertices. Each component follows a decay based on the shortest path
#' distance from its center.
#'
#' @param adj.list List of adjacency lists, where \code{adj.list[[i]]} contains indices of
#'        vertices adjacent to vertex \code{i}.
#' @param weight.list List of edge weights, where \code{weight.list[[i]][j]} is the weight
#'        of the edge from vertex \code{i} to \code{adj.list[[i]][j]}.
#' @param centers Vector of vertex indices to use as centers for the Gaussian components
#' @param amplitudes Vector of amplitudes for each Gaussian component
#' @param sigmas Vector of sigma values (spread) for each Gaussian component
#' @param normalize Logical; if TRUE, normalize the resulting function to \code{[0,1]}.
#'
#' @return A numeric vector of function values at each vertex of the graph
#'
#' @examples
#' \dontrun{
#' # Create a grid graph
#' grid <- create.grid.graph(10, 10)
#'
#' # Generate a mixture of two Gaussians
#' centers <- c(1, 100)  # Corner and center vertices
#' amplitudes <- c(1.0, 0.7)
#' sigmas <- c(2.0, 3.0)
#'
#' y <- generate.graph.gaussian.mixture(
#'   grid$adj.list,
#'   grid$weight.list,
#'   centers,
#'   amplitudes,
#'   sigmas
#' )
#' }
#'
#' @export
generate.graph.gaussian.mixture <- function(
  adj.list,
  weight.list,
  centers,
  amplitudes = rep(1.0, length(centers)),
  sigmas = rep(2.0, length(centers)),
  normalize = TRUE
) {
  # Verify inputs
  if (length(centers) != length(amplitudes) || length(centers) != length(sigmas)) {
    stop("centers, amplitudes, and sigmas must have the same length")
  }

  n_vertices <- length(adj.list)
  n_components <- length(centers)

  # Initialize result vector
  result <- numeric(n_vertices)

  # Compute distance matrices from each center
  distance_matrices <- list()
  for (center in centers) {
    # Initialize distances for Dijkstra's algorithm
    distances <- rep(Inf, n_vertices)
    distances[center] <- 0

    # Set of vertices whose shortest distance is determined
    visited <- logical(n_vertices)

    # Priority queue (implemented as a while loop with min selection)
    while (!all(visited)) {
      # Find unvisited vertex with minimum distance
      current <- which.min(ifelse(visited, Inf, distances))

      # If all remaining vertices are unreachable, break
      if (is.infinite(distances[current])) {
        break
      }

      # Mark as visited
      visited[current] <- TRUE

      # Update distances to neighbors
      for (i in seq_along(adj.list[[current]])) {
        neighbor <- adj.list[[current]][i]
        weight <- weight.list[[current]][i]

        # Calculate potential new distance
        new_dist <- distances[current] + weight

        # Update if better path found
        if (new_dist < distances[neighbor]) {
          distances[neighbor] <- new_dist
        }
      }
    }

    # Store distance matrix
    distance_matrices[[length(distance_matrices) + 1]] <- distances
  }

  # Compute mixture of Gaussians
  for (i in 1:n_components) {
    # Compute Gaussian function using distances
    gauss_values <- amplitudes[i] * exp(-(distance_matrices[[i]]^2) / (2 * sigmas[i]^2))

    # Add to mixture
    result <- result + gauss_values
  }

  # Normalize if requested
  if (normalize && max(result) > 0) {
    result <- result / max(result)
  }

  return(result)
}

#' Visualize a Function on a Grid Graph
#'
#' Creates multiple visualizations of a function defined on a grid graph including
#' heatmap with contours, 3D perspective plot, and optionally an interactive 3D plot.
#'
#' @param grid.size Integer; size of the square grid (grid.size x grid.size).
#' @param z Numeric vector; function values at each vertex in row-wise order (length = grid.size^2).
#' @param centers Optional integer vector; vertex indices to highlight (e.g., peaks).
#' @param title Character string; title for the plots (default: "Function on Grid Graph").
#'
#' @return Invisibly returns the input \code{z}.
#'
#' @details
#' Left panel: heatmap with contour lines (and optional centers).
#' Right panel: 3D perspective plot (base graphics).
#' If \pkg{rgl} is available, a separate off-screen rgl window is opened to render a 3D point view.
#'
#' @examples
#' \dontrun{
#' grid.size <- 20
#' x <- rep(1:grid.size, grid.size) / grid.size
#' y <- rep(1:grid.size, each = grid.size) / grid.size
#' z <- exp(-10*((x-0.3)^2 + (y-0.3)^2)) + 0.5*exp(-8*((x-0.7)^2 + (y-0.6)^2))
#' centers <- which(z > 0.9)
#' visualize.grid.function(grid.size, z, centers)
#'
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'   visualize.grid.function(grid.size, z, centers)
#' }
#' }
#' @export
visualize.grid.function <- function(grid.size, z, centers = NULL, title = "Function on Grid Graph") {
  ## ---- Validation -----------------------------------------------------------
  if (length(grid.size) != 1L || !is.numeric(grid.size) || grid.size < 2 || !is.finite(grid.size)) {
    stop("'grid.size' must be a single finite numeric value >= 2", call. = FALSE)
  }
  grid.size <- as.integer(grid.size)
  n_vertices <- grid.size^2L

  if (!is.numeric(z) || length(z) != n_vertices) {
    stop(sprintf("'z' must be a numeric vector of length grid.size^2 (= %d)", n_vertices), call. = FALSE)
  }
  if (any(!is.finite(z))) {
    stop("'z' contains non-finite values; please remove or impute NA/Inf/NaN", call. = FALSE)
  }

  if (!is.null(centers)) {
    centers <- as.integer(centers)
    centers <- centers[is.finite(centers)]
    centers <- centers[centers >= 1L & centers <= n_vertices]
    if (!length(centers)) centers <- NULL
  }

  ## ---- Coordinates & helpers -----------------------------------------------
  x_coords <- rep(seq_len(grid.size), grid.size) / grid.size
  y_coords <- rep(seq_len(grid.size), each = grid.size) / grid.size

  vertex_to_coords <- function(vertices) {
    xv <- (vertices - 1L) %% grid.size + 1L
    yv <- ceiling(vertices / grid.size)
    list(x = xv / grid.size, y = yv / grid.size)
  }

  center_coords <- if (!is.null(centers)) vertex_to_coords(centers) else NULL

  ## ---- Set layout and par safely -------------------------------------------
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)

  graphics::layout(matrix(c(1, 2), nrow = 1L, ncol = 2L))
  on.exit(try(graphics::layout(1), silent = TRUE), add = TRUE)

  ## ---- Plot 1: Heatmap + contours ------------------------------------------
  graphics::par(mar = c(4, 4, 2, 1))
  z_matrix <- matrix(z, nrow = grid.size, ncol = grid.size, byrow = FALSE)

  pal2d <- grDevices::hcl.colors(100, palette = "YlOrRd")
  graphics::image(
    x = seq_len(grid.size) / grid.size,
    y = seq_len(grid.size) / grid.size,
    z = z_matrix,
    col = pal2d,
    main = title,
    xlab = "X", ylab = "Y"
  )

  graphics::contour(
    x = seq_len(grid.size) / grid.size,
    y = seq_len(grid.size) / grid.size,
    z = z_matrix,
    add = TRUE,
    col = "black"
  )

  if (!is.null(centers)) {
    graphics::points(center_coords$x, center_coords$y, pch = 19, col = "blue", cex = 1.5)
    graphics::text(center_coords$x, center_coords$y,
                   labels = paste("Center", seq_along(centers)),
                   pos = 3, offset = 0.7, cex = 0.8)
  }

  ## ---- Plot 2: Base-graphics 3D perspective -------------------------------
  graphics::par(mar = c(4, 4, 2, 1))
  graphics::persp(
    x = seq_len(grid.size) / grid.size,
    y = seq_len(grid.size) / grid.size,
    z = z_matrix,
    theta = 30, phi = 30,
    expand = 0.7,
    col = "lightblue",
    shade = 0.5,
    main = "3D Perspective",
    xlab = "X", ylab = "Y", zlab = "Value"
  )

  ## ---- Optional rgl 3D view (headless-safe) --------------------------------
  if (requireNamespace("rgl", quietly = TRUE)) {
    old_rgl <- options(rgl.useNULL = TRUE)
    on.exit(options(old_rgl), add = TRUE)

    rgl::open3d()
    on.exit(try(rgl::rgl.close(), silent = TRUE), add = TRUE)

    pal3d <- grDevices::hcl.colors(length(z), palette = "Spectral")
    rgl::plot3d(
      x_coords, y_coords, z,
      col = pal3d[rank(z, ties.method = "average")],
      size = 3,
      xlab = "X", ylab = "Y", zlab = "Z",
      type = "p"
    )
    rgl::title3d(main = title)

    if (!is.null(centers)) {
      for (i in seq_along(centers)) {
        rgl::spheres3d(
          center_coords$x[i], center_coords$y[i], z[centers[i]],
          radius = 0.02,
          color = "blue"
        )
      }
    }

    rgl::axes3d()
  }

  invisible(z)
}

#' Create a Random Graph
#'
#' @description
#' Creates a random graph with the specified number of vertices and
#' average degree (neighbors per vertex).
#'
#' @param n_vertices Number of vertices in the graph
#' @param avg_degree Average number of neighbors per vertex
#' @param connected Logical; if TRUE, ensure the graph is connected
#'
#' @return A list with adjacency and weight lists
#'
#' @examples
#' graph <- create.random.graph(100, 4)
#'
#' @export
create.random.graph <- function(n_vertices, avg_degree, connected = TRUE) {
  # Initialize empty adjacency and weight lists
  adj.list <- vector("list", n_vertices)
  weight.list <- vector("list", n_vertices)

  # Total number of edges to create
  n_edges <- floor(n_vertices * avg_degree / 2)

  # Initialize all vertices with empty adjacency lists
  for (i in 1:n_vertices) {
    adj.list[[i]] <- integer(0)
    weight.list[[i]] <- numeric(0)
  }

  # First, ensure the graph is connected by creating a spanning tree
  if (connected) {
    # Start with vertex 1
    connected_vertices <- c(1)
    unconnected_vertices <- setdiff(1:n_vertices, connected_vertices)

    # Connect each vertex to one already in the connected set
    while (length(unconnected_vertices) > 0) {
      # Choose a random unconnected vertex
      v2 <- sample(unconnected_vertices, 1)

      # Choose a random connected vertex to link to
      v1 <- sample(connected_vertices, 1)

      # Add edge in both directions
      adj.list[[v1]] <- c(adj.list[[v1]], v2)
      weight.list[[v1]] <- c(weight.list[[v1]], 1.0)

      adj.list[[v2]] <- c(adj.list[[v2]], v1)
      weight.list[[v2]] <- c(weight.list[[v2]], 1.0)

      # Update sets
      connected_vertices <- c(connected_vertices, v2)
      unconnected_vertices <- setdiff(unconnected_vertices, v2)
    }

    # Adjust the remaining edges to add
    n_edges <- n_edges - (n_vertices - 1)
  }

  # Add random edges until we reach the desired average degree
  edges_added <- 0
  while (edges_added < n_edges) {
    # Choose two different random vertices
    vertices <- sample(n_vertices, 2)
    v1 <- vertices[1]
    v2 <- vertices[2]

    # Check if edge already exists
    if (!(v2 %in% adj.list[[v1]])) {
      # Add edge in both directions (undirected graph)
      adj.list[[v1]] <- c(adj.list[[v1]], v2)
      weight.list[[v1]] <- c(weight.list[[v1]], 1.0)

      adj.list[[v2]] <- c(adj.list[[v2]], v1)
      weight.list[[v2]] <- c(weight.list[[v2]], 1.0)

      edges_added <- edges_added + 1
    }
  }

  # Calculate a more meaningful weight for edges (random, but consistent)
  # Here we use a simple approach: random weights between 0.5 and 1.5
  set.seed(42)  # For reproducibility
  for (i in 1:n_vertices) {
    n_neighbors <- length(adj.list[[i]])
    if (n_neighbors > 0) {
      weight.list[[i]] <- runif(n_neighbors, 0.5, 1.5)
    }
  }

  return(list(adj.list = adj.list, weight.list = weight.list))
}


#' Assess Fidelity of Graph-Based Geodesic Distances to Euclidean Geometry
#'
#' This function compares local neighborhood structures defined by graph-based distances
#' to those defined in the original Euclidean space. It evaluates how well a given graph
#' preserves the local geometry of a dataset by computing two metrics for each data point:
#'
#' \itemize{
#'   \item \strong{Jaccard Index:} Measures the similarity between the sets of neighbors
#'         within a radius \eqn{\tau} in both the Euclidean and graph-based spaces.
#'   \item \strong{Mean Absolute Deviation (MAD):} Measures the distortion in local distances
#'         over the intersection of the two neighborhood sets.
#' }
#'
#' The graph is assumed to be provided as an adjacency list and a corresponding weight list.
#' Any duplicate edges (i.e., undirected edges appearing in both directions) are automatically deduplicated.
#'
#' @param X A numeric matrix of shape \code{[n, d]}, where each row represents a data point in \code{d}-dimensional space.
#' @param adj.list A list of integer vectors of length \code{n}, giving the graph adjacency list. Each entry contains the 1-based indices of neighbors for that vertex.
#' @param weight.list A list of numeric vectors of same length as \code{adj.list}, where each element contains edge weights for the corresponding neighbors.
#' @param taus A numeric vector of radius values \eqn{\tau} used to define local neighborhoods around each point.
#' @param max.k Integer; the number of nearest neighbors to compute in the Euclidean space (used to approximate local neighborhoods).
#'
#' @return A named list of length equal to \code{length(taus)}. Each element is a list with:
#' \describe{
#'   \item{\code{mean.jaccard}}{Mean Jaccard index over all data points for a given \eqn{\tau}.}
#'   \item{\code{mean.mad}}{Mean absolute deviation in distances between Euclidean and graph metrics over intersecting neighbors.}
#'   \item{\code{jaccard.vals}}{A numeric vector of Jaccard index values per vertex.}
#'   \item{\code{mad.vals}}{A numeric vector of MAD values per vertex.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(300 * 2), ncol = 2)
#' g <- build_graph(X)  # your function to generate adj.list and weight.list
#' taus <- seq(0.05, 0.2, by = 0.01)
#' fidelity <- compute.local.distance.fidelity(X, g$adj.list, g$weight.list, taus)
#' }
#'
#' @importFrom FNN get.knn
#' @importFrom igraph graph_from_adj_list as_edgelist graph_from_data_frame E distances
#'
#' @export
compute.local.distance.fidelity <- function(X, adj.list, weight.list, taus, max.k = 50) {
    n <- nrow(X)

    ## 1. Euclidean distances
    nn.X <- FNN::get.knn(X, k = max.k)
    D.X <- as.matrix(dist(X))  ## For true Euclidean distances

    ## 2. Graph distances
    G <- igraph::graph_from_adj_list(adj.list, mode="all")
    ## Get edge list and remove duplicates
    edgelist <- igraph::as_edgelist(G)
    edge.df <- data.frame(from = pmin(edgelist[,1], edgelist[,2]),
                          to   = pmax(edgelist[,1], edgelist[,2]),
                          weight = unlist(weight.list))
    ## Remove duplicates
    edge.df <- edge.df[!duplicated(edge.df), ]
    ## Rebuild the graph and assign deduplicated weights
    G <- igraph::graph_from_data_frame(edge.df, directed = FALSE)
    igraph::E(G)$weight <- edge.df$weight

    D.G <- igraph::distances(G, v=1:n, to=1:n, weights=E(G)$weight)

    results <- list()
    for (tau.idx in seq(taus)) {
        tau <- taus[tau.idx]
        mad_vals <- numeric(n)
        jaccard_vals <- numeric(n)

        for (i in 1:n) {
            A <- which(D.X[i, ] < tau)
            B <- which(D.G[i, ] < tau)

            I <- intersect(A, B)
            U <- union(A, B)

            if (length(U) > 0) {
                jaccard_vals[i] <- length(I) / length(U)
            }

            if (length(I) > 0) {
                d_x <- D.X[i, I]
                d_g <- D.G[i, I]
                mad_vals[i] <- mean(abs(d_x - d_g))
            }
        }

        results[[tau.idx]] <- list(
            mean.jaccard = mean(jaccard_vals, na.rm = TRUE),
            mean.mad = mean(mad_vals, na.rm = TRUE),
            jaccard.vals = jaccard_vals,
            mad.vals = mad_vals
        )
    }

    return(results)
}


compute.kernel.graph.laplacian.eigenfunctions.I.minus.L.powered <- function(x.min,
                                                                            x.max,
                                                                            grid.size,
                                                                            x.grid = NULL,
                                                                            random = TRUE,
                                                                            n.evects,
                                                                            tau,
                                                                            power = 1) {
    if (is.null(x.grid)) {
        if (random) {
            x.grid <- sort(runif(grid.size, min = x.min, max = x.max))
        } else {
            x.grid <- seq(x.min, x.max, length.out = grid.size) ## Uniform grid
        }
    }

    ## Compute pairwise squared distances
    distsq.mat <- outer(x.grid, x.grid, function(a, b) (a - b)^2)

    ## Construct kernel-based affinity matrix W
    W <- exp(-distsq.mat / tau^2)
    diag(W) <- 0  # No self-loops

    ## Degree matrix D
    D <- diag(rowSums(W))

    ## Unnormalized Graph Laplacian
    L <- D - tau * W

    ## Construct (I - L)
    I_minus_L <- diag(grid.size) - L

    ## Raise (I - L) to the given power
    L_mod <- I_minus_L
    if (power > 1) {
        for (i in 2:power) {
            L_mod <- L_mod %*% I_minus_L
        }
    }

    ## Eigen decomposition
    eig <- eigen(L_mod, symmetric = TRUE)

    list(
        x.grid = x.grid,
        eigenvectors = eig$vectors[, 1:n.evects, drop = FALSE],
        eigenvalues = eig$values[1:n.evects],
        L_mod = L_mod,
        W = W,
        D = D,
        power = power
    )
}

## compare.kernel.graph.laplacian.with.continuum.laplacian <- function(result, x.min, x.max, n.evects = NULL) {
##     x.grid <- result$x.grid

##     if (is.null(n.evects)) {
##         n.evects <- ncol(result$eigenvectors)
##     }

##     if (n.evects == 1) {
##         op <- par(mfrow=c(1,1))
##     } else {
##         op <- par(mfrow = c(ceiling(n.evects / 2), 2), mar = c(2, 2, 2, 1))
##     }

##     for (j in 1:n.evects) {
##         ## Classical continuum Laplacian eigenfunctions (sine waves)
##         phi.j <- sin(j * pi * (x.grid - x.min) / (x.max - x.min))

##         ## Normalize both to have unit norm (L2)
##         phi.j <- phi.j / sqrt(sum(phi.j^2))
##         evect.j <- result$eigenvectors[, j]
##         evect.j <- evect.j / sqrt(sum(evect.j^2))

##         plot(x.grid, evect.j, type = "l", col = "blue", lwd = 2, las = 1,
##              ylab = "", xlab = "", main = paste("Eigenfunction", j))
##         lines(x.grid, phi.j, col = "red", lty = 2, lwd = 2)
##         ## legend("topright", legend = c(paste0(\"(I-L)^\", result$power, \" eigen\"), \"Classical (sine)\"),
##         ##        col = c(\"blue\", \"red\"), lty = c(1, 2), cex = 0.8, inset = 0.025)
##     }

##     par(op)
## }


#' Create Threshold Distance Graph from Distance Matrix
#'
#' This function constructs a graph where vertices are connected by an edge
#' if and only if their distance is less than a specified threshold.
#' The weight of each edge is the corresponding distance.
#'
#' @param dist.matrix A symmetric distance matrix where rows and columns correspond to vertices
#' @param threshold A numeric threshold value; vertices i and j are connected if dist(i,j) < threshold
#' @param include.names Logical; whether to include vertex names in the output (default: TRUE)
#'
#' @return A list with two components:
#'   \item{adj_list}{A list of integer vectors. Each vector contains the indices of vertices
#'         adjacent to the corresponding vertex.}
#'   \item{weight_list}{A list of numeric vectors. Each vector contains weights of edges
#'         corresponding to adjacencies in adj_list.}
#'
#' @examples
#' # Example distance matrix
#' dist <- matrix(c(
#'   0.0000000, 0.02834008, 0.05050505, 0.12500000, 0.1086957,
#'   0.02834008, 0.00000000, 0.88888889, 0.54166667, 1.0000000,
#'   0.05050505, 0.88888889, 0.00000000, 0.04166667, 0.1086957,
#'   0.12500000, 0.54166667, 0.04166667, 0.00000000, 1.0000000,
#'   0.10869565, 1.00000000, 0.10869565, 1.00000000, 0.0000000
#' ), nrow=5, byrow=TRUE)
#' rownames(dist) <- colnames(dist) <- c("M1", "M2", "M3", "M4", "M5")
#'
#' # Force symmetry by averaging with transpose
#' dist <- (dist + t(dist)) / 2
#'
#' # Create graph with threshold 0.15
#' graph <- create.threshold.distance.graph(dist, 0.15)
#'
#' @export
create.threshold.distance.graph <- function(dist.matrix, threshold, include.names = TRUE) {
    ## Check if the matrix is symmetric
    if (!isSymmetric(unname(dist.matrix))) {
        stop("The distance matrix must be symmetric")
    }

    ## Get the number of vertices
    n.vertices <- nrow(dist.matrix)

    ## Get vertex names if available
    vertex.names <- rownames(dist.matrix)
    if (is.null(vertex.names)) {
        vertex.names <- 1:n.vertices
    }

    ## Initialize adjacency and weight lists
    adj.list <- vector("list", n.vertices)
    weight.list <- vector("list", n.vertices)

    ## Populate the adjacency and weight lists
    for (i in 1:n.vertices) {
        ## Find all vertices j where dist(i,j) < threshold and i != j
        neighbors <- which(dist.matrix[i, ] < threshold & (1:n.vertices != i))

        ## Add to adjacency list
        adj.list[[i]] <- neighbors

        ## Add corresponding weights
        weight.list[[i]] <- dist.matrix[i, neighbors]
    }

    ## Add names if requested
    if (include.names && !is.null(vertex.names)) {
        names(adj.list) <- vertex.names
        names(weight.list) <- vertex.names
    }

    return(list(
        adj_list = adj.list,
        weight_list = weight.list
    ))
}


#' Compute local cluster evenness with full cluster support
#'
#' Computes the evenness (normalized entropy) of cluster label distribution
#' in the expanded neighborhood of each vertex. The entropy is computed using
#' the full set of cluster labels, including those not present in the neighborhood
#' (i.e., with zero frequency), to ensure comparability across vertices.
#'
#' @param adj.list List of adjacency vectors (1-based indices) for each vertex.
#' @param cltr A vector of cluster labels (factor, character, or numeric), one per vertex.
#'
#' @return Numeric vector of evenness values for each vertex.
#' @export
graph.cltr.evenness <- function(adj.list, cltr) {
    n.vertices <- length(adj.list)
    uq.clusters <- unique(cltr)

    cltr.evenness <- numeric(n.vertices)
    for (i in seq_len(n.vertices)) {
        nbrs <- c(i, adj.list[[i]])
                                        # Include all clusters, even if count is 0
        cltr.labels <- factor(cltr[nbrs], levels = uq.clusters)
        nbrs.cltr.freq <- table(cltr.labels)
        cltr.evenness[i] <- evenness(as.numeric(nbrs.cltr.freq))
    }

    return(cltr.evenness)
}


#' Compute the Diameter of a Weighted Undirected Graph
#'
#' @description
#' This function calculates the diameter of a weighted undirected graph
#' represented as adjacency and weight lists. The diameter is the
#' length of the longest shortest path between any two vertices in the graph.
#'
#' @param adj.list A list where each element i contains the indices of vertices
#'        adjacent to vertex i
#' @param weight.list A list where each element i contains the weights of edges
#'        connecting vertex i to its adjacent vertices in adj.list
#'
#' @return A list containing:
#'   \item{diameter}{The diameter of the graph (numeric value, or Inf if disconnected)}
#'   \item{message}{A descriptive message about the diameter}
#'   \item{farthest_vertices}{The pair of vertices that are farthest apart}
#'   \item{diameter_path}{The shortest path between the farthest vertices}
#'
#' @details
#' The function handles several special cases:
#' - If the graph has no edges, returns NA as the diameter
#' - If the graph is disconnected, returns Inf as the diameter
#' - For undirected graphs, edges are only added once to avoid duplication
#'
#' @note
#' Requires the 'igraph' package to be installed
#'
#' @examples
#' # Example with a simple graph
#' adj.list <- list(c(2,3), c(1,3), c(1,2))
#' weight.list <- list(c(1,2), c(1,3), c(2,3))
#' result <- compute.graph.diameter(adj.list, weight.list)
#' print(result$message)
#'
#' @importFrom igraph graph_from_edgelist E diameter
#'             farthest_vertices shortest_paths is_connected
#'
#' @export
compute.graph.diameter <- function(adj.list, weight.list) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Please install the igraph package: install.packages('igraph')")
    }

    res <- convert.adjacency.to.edge.matrix(adj.list, weight.list)

    g <- igraph::graph_from_edgelist(res$edge.matrix, directed = FALSE)
    igraph::E(g)$weight <- res$weights

    ## Compute diameter
    diam <- igraph::diameter(g, weights = igraph::E(g)$weight, directed = FALSE)

    ## Find farthest vertices
    farthest <- igraph::farthest_vertices(g, weights = igraph::E(g)$weight, directed = FALSE)

    ## Get the path that gives the diameter
    path <- igraph::shortest_paths(g,
                           from = farthest$vertices[1],
                           to = farthest$vertices[2],
                           weights = igraph::E(g)$weight,
                           output = "both")

    ## Return results as a list
    return(list(
        diameter = diam,
        message = paste("The diameter of the graph is:", diam),
        farthest_vertices = farthest,
        diameter_path = path
    ))
}

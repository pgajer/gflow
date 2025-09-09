#' Parameterize a Circular Graph Structure for Biological Cycle Analysis
#'
#' @description
#' Applies spectral methods to parameterize a graph with a circular structure,
#' calculating the position of each vertex along a circle using eigenvectors of
#' the graph Laplacian. This function was motivated by the discovery that cell
#' cycle genes form a circular structure in high-dimensional expression space,
#' enabling cell cycle position assignment.
#'
#' @details
#' This function embeds a graph with circular topology onto a circle using the
#' second and third eigenvectors of the graph Laplacian matrix. The algorithm
#' computes an angle for each vertex, positioning it on the unit circle in a way
#' that preserves the graph structure.
#'
#' The method is particularly useful for biological applications where the
#' underlying process has circular dynamics, such as the cell cycle. In such
#' cases, genes or cells can be ordered along a circle representing progression
#' through the cyclic process. The parameterization allows for synthetic
#' time-within-cycle assignment, where the computed angles represent positions
#' along the circular trajectory.
#'
#' The graph Laplacian is constructed as L = D - W, where D is the degree matrix
#' and W is the weighted adjacency matrix. When \code{use.edge.lengths = TRUE},
#' edge weights are incorporated into the Laplacian construction. The method
#' works best for graphs with approximately circular connectivity patterns and
#' average vertex degree between 3 and 5.
#'
#' @param adj.list A list of integer vectors representing the adjacency structure.
#'   Each element \code{adj.list[[i]]} contains the indices of vertices adjacent
#'   to vertex \code{i}. Uses 1-based indexing (R convention). In biological
#'   applications, vertices might represent genes or cells with similar expression
#'   patterns.
#' @param weight.list A list of numeric vectors containing edge weights. Each
#'   element \code{weight.list[[i]]} contains the weights for edges from vertex
#'   \code{i} to the vertices listed in \code{adj.list[[i]]}. Must have the same
#'   structure as \code{adj.list}. Weights can represent similarity or distance
#'   metrics between expression profiles.
#' @param use.edge.lengths Logical. If \code{TRUE}, edge weights from
#'   \code{weight.list} are used in constructing the graph Laplacian. If
#'   \code{FALSE}, all edges are treated as having unit weight. Default is
#'   \code{TRUE}.
#'
#' @return A list of class \code{"circular_parameterization"} containing:
#'   \describe{
#'     \item{\code{angles}}{Numeric vector of angles in radians (range \eqn{[0, 2\pi]})
#'       for each vertex, representing their position on the unit circle. In
#'       cell cycle applications, these angles represent the inferred position
#'       within the cycle.}
#'     \item{\code{eig_vec2}}{Numeric vector containing the second eigenvector of
#'       the graph Laplacian.}
#'     \item{\code{eig_vec3}}{Numeric vector containing the third eigenvector of
#'       the graph Laplacian.}
#'   }
#'
#' @section Note:
#' The function internally converts R's 1-based indexing to 0-based indexing
#' before calling the underlying C++ implementation.
#'
#' @references
#' Zheng, S. C., Stein-O'Brien, G., Augustin, J. J., Slosberg, J., Carosso, G. A.,
#' Winer, B., ... & Hansen, K. D. (2022). Universal prediction of cell-cycle
#' position using transfer learning. Genome Biology, 23(1), 41.
#' \doi{10.1186/s13059-021-02581-y}
#'
#' @examples
#' # Example 1: Simple circular graph (like cell cycle progression)
#' n <- 8
#' adj.list <- lapply(seq_len(n), function(i) {
#'   c(if (i == n) 1 else i + 1,  # next vertex
#'     if (i == 1) n else i - 1)  # previous vertex
#' })
#' weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))
#'
#' # Get circular parameterization
#' result <- parameterize.circular.graph(adj.list, weight.list, TRUE)
#'
#' # Display angles (positions along the cycle)
#' print(round(result$angles, 2))
#'
#' # Example 2: Biological application - genes with circular expression pattern
#' # Simulate a gene similarity network where genes are connected if their
#' # expression patterns are similar (simplified example)
#' n_genes <- 12
#' # Create connections based on proximity in the cycle
#' adj.list <- lapply(seq_len(n_genes), function(i) {
#'   # Connect to 2 neighbors on each side to simulate local similarity
#'   neighbors <- c((i - 2):(i - 1), (i + 1):(i + 2))
#'   neighbors <- ((neighbors - 1) %% n_genes) + 1
#'   neighbors[neighbors != i]  # Remove self-connections
#' })
#' # Weights decrease with distance in the cycle
#' weight.list <- lapply(seq_len(n_genes), function(i) {
#'   neighbors <- adj.list[[i]]
#'   weights <- sapply(neighbors, function(j) {
#'     dist <- min(abs(i - j), n_genes - abs(i - j))
#'     exp(-dist/2)  # Exponential decay
#'   })
#'   weights
#' })
#'
#' result <- parameterize.circular.graph(adj.list, weight.list)
#'
#' # Plot genes positioned by their inferred cycle position
#' plot(cos(result$angles), sin(result$angles),
#'      xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
#'      pch = 19, xlab = "x", ylab = "y", asp = 1,
#'      main = "Gene Positions in Expression Cycle")
#'
#' # Add gene labels
#' text(1.1 * cos(result$angles), 1.1 * sin(result$angles),
#'      labels = paste0("G", seq_len(n_genes)), cex = 0.8)
#'
#' # The angles can be interpreted as positions in the biological cycle
#' # For cell cycle: 0 = G1/S, pi/2 = S, pi = G2, 3*pi/2 = M phase
#'
#' @export
parameterize.circular.graph <- function(adj.list,
                                        weight.list,
                                        use.edge.lengths = TRUE) {
  # Input validation
  if (!is.list(adj.list)) {
    stop("'adj.list' must be a list", call. = FALSE)
  }

  if (!is.list(weight.list)) {
    stop("'weight.list' must be a list", call. = FALSE)
  }

  if (length(adj.list) != length(weight.list)) {
    stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
  }

  n <- length(adj.list)
  if (n == 0) {
    stop("'adj.list' must not be empty", call. = FALSE)
  }

  # Validate use.edge.lengths
  if (!is.logical(use.edge.lengths) || length(use.edge.lengths) != 1) {
    stop("'use.edge.lengths' must be a single logical value", call. = FALSE)
  }

  if (is.na(use.edge.lengths)) {
    stop("'use.edge.lengths' must not be NA", call. = FALSE)
  }

  # Validate each element of adj.list and weight.list
  for (i in seq_along(adj.list)) {
    # Check adj.list[[i]]
    if (!is.numeric(adj.list[[i]])) {
      stop(sprintf("adj.list[[%d]] must be numeric", i), call. = FALSE)
    }

    if (any(is.na(adj.list[[i]]))) {
      stop(sprintf("adj.list[[%d]] contains NA values", i), call. = FALSE)
    }

    if (any(adj.list[[i]] != round(adj.list[[i]]))) {
      stop(sprintf("adj.list[[%d]] must contain only integer values", i),
           call. = FALSE)
    }

    # Check vertex indices are in valid range
    if (length(adj.list[[i]]) > 0) {
      if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n)) {
        stop(sprintf("adj.list[[%d]] contains vertex indices outside valid range [1, %d]",
                     i, n), call. = FALSE)
      }
    }

    # Check weight.list[[i]]
    if (!is.numeric(weight.list[[i]])) {
      stop(sprintf("weight.list[[%d]] must be numeric", i), call. = FALSE)
    }

    if (any(is.na(weight.list[[i]]))) {
      stop(sprintf("weight.list[[%d]] contains NA values", i), call. = FALSE)
    }

    if (any(!is.finite(weight.list[[i]]))) {
      stop(sprintf("weight.list[[%d]] contains non-finite values", i),
           call. = FALSE)
    }

    if (any(weight.list[[i]] < 0)) {
      stop(sprintf("weight.list[[%d]] contains negative weights", i),
           call. = FALSE)
    }

    # Check matching lengths
    if (length(adj.list[[i]]) != length(weight.list[[i]])) {
      stop(sprintf("Vertex %d has %d edges but %d weights",
                   i, length(adj.list[[i]]), length(weight.list[[i]])),
           call. = FALSE)
    }
  }

  # Convert to 0-based indexing for C++
  adj.list.0based <- lapply(adj.list, function(x) {
    if (length(x) == 0) {
      integer(0)
    } else {
      as.integer(x - 1L)
    }
  })

  # Call the C++ function
  result <- tryCatch({
    .Call("S_parameterize_circular_graph",
          adj.list.0based,
          weight.list,
          as.logical(use.edge.lengths))
  }, error = function(e) {
    stop("Error in C++ function call: ", e$message, call. = FALSE)
  })

  # Add class to result for potential S3 methods
  class(result) <- c("circular_parameterization", "list")

  return(result)
}

#' Print Method for circular_parameterization Objects
#'
#' @description
#' Prints a summary of a circular parameterization result.
#'
#' @param x An object of class \code{"circular_parameterization"} as returned by
#'   \code{\link{parameterize.circular.graph}}.
#' @param digits Number of digits to display for numeric values. Default is 4.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' # Create example graph
#' n <- 6
#' adj.list <- lapply(seq_len(n), function(i) {
#'   c(if (i == n) 1 else i + 1, if (i == 1) n else i - 1)
#' })
#' weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))
#' result <- parameterize.circular.graph(adj.list, weight.list)
#' print(result)
#'
#' @export
print.circular_parameterization <- function(x, digits = 4, ...) {
  cat("Circular Graph Parameterization\n")
  cat("-------------------------------\n")
  cat("Number of vertices:", length(x$angles), "\n")
  cat("\nVertex angles (radians):\n")
  print(round(x$angles, digits = digits))
  cat("\nAngle range: [",
      round(min(x$angles), digits), ", ",
      round(max(x$angles), digits), "]\n", sep = "")
  invisible(x)
}

#' Plot Method for circular_parameterization Objects
#'
#' @description
#' Creates a visualization of vertices positioned on a circle according to their
#' computed angles.
#'
#' @param x An object of class \code{"circular_parameterization"} as returned by
#'   \code{\link{parameterize.circular.graph}}.
#' @param adj.list Optional adjacency list to draw edges. If provided, edges will
#'   be drawn between connected vertices.
#' @param vertex.labels Labels for vertices. Default is \code{seq_len(n)} where
#'   \code{n} is the number of vertices. Use \code{NA} to suppress labels.
#' @param vertex.cex Character expansion factor for vertex points. Default is 1.5.
#' @param label.cex Character expansion factor for vertex labels. Default is 0.8.
#' @param edge.col Color for edges. Default is "gray70".
#' @param vertex.col Color for vertex points. Default is "black".
#' @param ... Further graphical parameters passed to \code{\link{plot}}.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' # Create example graph
#' n <- 8
#' adj.list <- lapply(seq_len(n), function(i) {
#'   c(if (i == n) 1 else i + 1, if (i == 1) n else i - 1)
#' })
#' weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))
#' result <- parameterize.circular.graph(adj.list, weight.list)
#'
#' # Basic plot
#' plot(result)
#'
#' # Plot with edges
#' plot(result, adj.list = adj.list,
#'      main = "Circular Graph Visualization")
#'
#' @export
plot.circular_parameterization <- function(x, adj.list = NULL,
                                           vertex.labels = seq_len(length(x$angles)),
                                           vertex.cex = 1.5,
                                           label.cex = 0.8,
                                           edge.col = "gray70",
                                           vertex.col = "black",
                                           ...) {
  n <- length(x$angles)

  # Compute vertex positions
  x_pos <- cos(x$angles)
  y_pos <- sin(x$angles)

  # Set up plot
  plot(x_pos, y_pos,
       xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
       pch = 19, cex = vertex.cex, col = vertex.col,
       asp = 1, xlab = "", ylab = "", ...)

  # Add reference circle
  theta <- seq(0, 2 * pi, length.out = 100)
  lines(cos(theta), sin(theta), col = "lightgray", lty = 2)

  # Draw edges if adjacency list is provided
  if (!is.null(adj.list)) {
    for (i in seq_len(n)) {
      for (j in adj.list[[i]]) {
        if (i < j) {  # Avoid drawing each edge twice
          segments(x_pos[i], y_pos[i], x_pos[j], y_pos[j],
                   col = edge.col, lwd = 1.5)
        }
      }
    }
  }

  # Redraw vertices on top of edges
  points(x_pos, y_pos, pch = 19, cex = vertex.cex, col = vertex.col)

  # Add labels if requested
  if (!is.null(vertex.labels) && !all(is.na(vertex.labels))) {
    text(1.15 * x_pos, 1.15 * y_pos,
         labels = vertex.labels, cex = label.cex)
  }

  invisible(x)
}

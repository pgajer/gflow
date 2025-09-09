#' Test Monotonic Reachability Map
#'
#' @description
#' Tests the monotonic reachability map functionality by finding paths that optimize
#' monotonicity of a function over a graph.
#'
#' @param adj.list List of adjacency lists, where adj.list[[i]] contains the indices
#'        of vertices adjacent to vertex i
#' @param weight.list List of edge weights, where weight.list[[i]][j] is the weight
#'        of the edge from vertex i to adj.list[[i]][j]
#' @param y Numeric vector of function values at each graph vertex
#' @param ref.vertex Index of the reference vertex (starting point)
#' @param radius Maximum distance to search from reference vertex
#' @param ascending Logical; if TRUE, find paths of consistent ascent; if FALSE, find
#'        paths of consistent descent
#' @param plot Logical; if TRUE, create plots of the results
#'
#' @return A list with the results, including vertices, monotonicity indices, paths,
#'         and the best gradient direction
#'
#' @examples
#' \dontrun{
#' # Create a grid graph with a peak function
#' n <- 10
#' grid <- create.grid.graph(n, n)
#' x <- rep(1:n, n) / n
#' y <- rep(1:n, each=n) / n
#' z <- exp(-((x-0.5)^2 + (y-0.5)^2) * 10)
#'
#' # Test monotonic paths from a corner to find paths up the peak
#' result <- test.monotonic.reachability(
#'   grid$adj.list,
#'   grid$weight.list,
#'   z,
#'   ref.vertex = 1,  # Corner vertex
#'   radius = n,
#'   ascending = TRUE
#' )
#' }
#'
#' @export
test.monotonic.reachability <- function(
                                        adj.list,
                                        weight.list,
                                        y,
                                        ref.vertex,
                                        radius,
                                        ascending = TRUE
                                        ) {

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call the C++ function
    result <- .Call(
        "S_test_monotonic_reachability_map",
        adj.list.0based,
        weight.list,
        y,
        as.integer(ref.vertex), ## it is being transformed to 0-based in S_test_monotonic_reachability_map()
        as.double(radius),
        as.logical(ascending)
    )

    ## Add y values to result for convenience
    result$y_values <- y[result$vertices]

    ## Add reference vertex information
    result$ref_vertex <- ref.vertex
    result$ref_y <- y[ref.vertex]

                                        # Process best vertex info
    if (result$best_vertex > 0) {
        result$best_path <- result$paths[[which(result$vertices == result$best_vertex)]]
        result$best_y_values <- y[result$best_path]
        result$best_y_change <- diff(result$best_y_values)
        result$best_cum_abs_change <- sum(abs(result$best_y_change))
        result$best_total_change <- result$best_y_values[length(result$best_y_values)] -
            result$best_y_values[1]
    }

    class(result) <- "monotonic_reachability"

    return(result)
}

plot.monotonic_reachability <- function(result) {

    ## Create visualization if requested
    if (plot) {
        par(mfrow = c(2, 2))

        ## Plot 1: Monotonicity indices
        plot(
            result$monotonicity,
            main = paste("Monotonicity Indices -",
                         ifelse(ascending, "Ascending", "Descending")),
            xlab = "Vertex Rank",
            ylab = "Monotonicity Index",
            type = "h",
            col = "blue"
        )
        abline(h = 1, lty = 2, col = "red") ## Perfect monotonicity reference

        ## Plot 2: Scatter of monotonicity vs. total change
        plot(
            result$total_change,
            result$monotonicity,
            main = "Monotonicity vs. Total Change",
            xlab = "Total y Change",
            ylab = "Monotonicity Index",
            pch = 19,
            col = "darkgreen"
        )
        if (result$best_vertex > 0) {
            idx <- which(result$vertices == result$best_vertex)
            points(result$total_change[idx], result$monotonicity[idx],
                   col = "red", pch = 19, cex = 1.5)
        }

        ## Plot 3: Best path y-values
        if (result$best_vertex > 0) {
            plot(
                result$best_y_values,
                type = "b",
                main = "Best Path y-values",
                xlab = "Path Position",
                ylab = "y-value",
                pch = 19,
                col = "purple"
            )
            points(1, result$best_y_values[1], pch = 17, col = "red", cex = 1.5) ## Reference vertex
            points(length(result$best_y_values),
                   result$best_y_values[length(result$best_y_values)],
                   pch = 18, col = "blue", cex = 1.5) ## Best vertex

            ## Add monotonicity information
            legend(
                "topleft",
                legend = c(
                    paste("MI =", round(result$best_monotonicity, 4)),
                    paste("Total Change =", round(result$best_total_change, 4)),
                    paste("Path Length =", length(result$best_path))
                ),
                bty = "n"
            )
        } else {
            plot(1, type = "n", main = "No Valid Path Found", xlab = "", ylab = "")
        }

        ## Plot 4: y-changes along best path
        if (result$best_vertex > 0 && length(result$best_y_values) > 1) {
            plot(
                result$best_y_change,
                type = "h",
                main = "y-changes Along Best Path",
                xlab = "Path Step",
                ylab = "Change in y",
                col = ifelse(result$best_y_change >= 0, "darkgreen", "darkred")
            )
            abline(h = 0, lty = 2, col = "gray")
        } else {
            plot(1, type = "n", main = "No Valid Path Found", xlab = "", ylab = "")
        }

        par(mfrow = c(1, 1))
    }

    ## Return the result
    invisible(result)
}

##' Create a 2D Grid Graph
##'
##' @description
##' Creates a graph representing a 2D grid with edges between adjacent cells.
##'
##' @param nx Number of vertices in x-direction
##' @param ny Number of vertices in y-direction
##' @param diagonals Logical; if TRUE, include diagonal connections
##'
##' @return A list with components:
##'   \item{adj.list}{List of adjacency lists}
##'   \item{weight.list}{List of weight lists}
##'
##' @examples
##' grid <- create.grid.graph(10, 10)
##'
##' @export
create.grid.graph <- function(nx, ny, diagonals = FALSE) {
    n_vertices <- nx * ny
    adj.list <- vector("list", n_vertices)
    weight.list <- vector("list", n_vertices)

    ## Helper to convert (x,y) coordinates to vertex index
    coord_to_index <- function(x, y) {
        return((y-1) * nx + x)
    }

    ## Helper to convert vertex index to (x,y) coordinates
    index_to_coord <- function(idx) {
        y <- ceiling(idx / nx)
        x <- idx - (y-1) * nx
        return(c(x, y))
    }

    ## Create edges
    for (i in 1:n_vertices) {
        coords <- index_to_coord(i)
        x <- coords[1]
        y <- coords[2]

        ## Initialize empty lists
        adj.list[[i]] <- integer(0)
        weight.list[[i]] <- numeric(0)

        ## Check all potential neighbors
        neighbor_offsets <- list(
            c(1, 0),   ## Right
            c(-1, 0),  ## Left
            c(0, 1),   ## Up
            c(0, -1)   ## Down
        )

        if (diagonals) {
            neighbor_offsets <- c(
                neighbor_offsets,
                list(
                    c(1, 1),    ## Upper-right
                    c(-1, 1),   ## Upper-left
                    c(1, -1),   ## Lower-right
                    c(-1, -1)   ## Lower-left
                )
            )
        }

        for (offset in neighbor_offsets) {
            nx_pos <- x + offset[1]
            ny_pos <- y + offset[2]

            ## Check if neighbor is within bounds
            if (nx_pos >= 1 && nx_pos <= nx && ny_pos >= 1 && ny_pos <= ny) {
                neighbor_idx <- coord_to_index(nx_pos, ny_pos)

                ## Calculate Euclidean distance for edge weight
                weight <- sqrt(offset[1]^2 + offset[2]^2)

                ## Add to adjacency and weight lists
                adj.list[[i]] <- c(adj.list[[i]], neighbor_idx)
                weight.list[[i]] <- c(weight.list[[i]], weight)
            }
        }
    }

    return(list(adj.list = adj.list, weight.list = weight.list))
}

##' Test Monotonic Reachability on a Gaussian Peak
##'
##' @description
##' Demonstrates the monotonic reachability algorithm on a 2D grid with a
##' Gaussian peak function.
##'
##' @param grid.size Size of the grid (grid.size x grid.size)
##' @param peak.center Center coordinates of the peak (c(x, y))
##' @param peak.scale Scale factor for the peak (controls steepness)
##' @param start.vertex Starting vertex for the monotonic paths (if NULL, uses corner)
##' @param radius Search radius
##' @param ascending Logical; if TRUE, find ascending paths; if FALSE, descending
##'
##' @return Results from test.monotonic.reachability
##'
##' @examples
##' result <- test_gaussian_peak_monotonicity()
##'
##' @export
test_gaussian_peak_monotonicity <- function(
                                            grid.size = 20,
                                            peak.center = c(0.5, 0.5),
                                            peak.scale = 10,
                                            start.vertex = NULL,
                                            radius = NULL,
                                            ascending = TRUE
                                            ) {
    ## Create grid graph
    grid <- create.grid.graph(grid.size, grid.size)

    ## Generate coordinates
    x_coords <- rep(1:grid.size, grid.size) / grid.size
    y_coords <- rep(1:grid.size, each=grid.size) / grid.size

    ## Generate Gaussian peak function
    z <- exp(
        -((x_coords - peak.center[1])^2 + (y_coords - peak.center[2])^2) * peak.scale
    )

    ## Use corner vertex if not specified
    if (is.null(start.vertex)) {
        start.vertex <- 1  ## Corner vertex
    }

    ## Use grid.size as default radius if not specified
    if (is.null(radius)) {
        radius <- grid.size
    }

    ## Run test
    result <- test.monotonic.reachability(
        grid$adj.list,
        grid$weight.list,
        z,
        start.vertex,
        radius,
        ascending
    )

    ## Create additional 3D visualization
    if (requireNamespace("rgl", quietly = TRUE)) {
        ## Grid of 3D points
        plot_gaussian_peak_3d(
            grid.size,
            z,
            result$best_path,
            peak.center,
            ascending
        )
    } else {
        message("Package 'rgl' is required for 3D visualization. Please install it with install.packages('rgl').")
    }

    return(result)
}

##' Plot Gaussian Peak with Best Path in 3D
##'
##' @param grid.size Size of the grid
##' @param z Vector of height values
##' @param path Path vertices to highlight
##' @param peak.center Center of peak
##' @param ascending Whether path is ascending or descending
##'
##' @importFrom rgl plot3d points3d lines3d
plot_gaussian_peak_3d <- function(
                                  grid.size,
                                  z,
                                  path = NULL,
                                  peak.center = c(0.5, 0.5),
                                  ascending = TRUE
                                  ) {
    ## Check if package is available
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package 'rgl' is required for this function")
    }

    ## Generate coordinates
    x_coords <- rep(1:grid.size, grid.size) / grid.size
    y_coords <- rep(1:grid.size, each=grid.size) / grid.size

    ## Create 3D plot
    rgl::open3d()

    ## Plot surface
    rgl::plot3d(
             x_coords, y_coords, z,
             col = heat.colors(length(z)),
             size = 3,
             xlab = "X", ylab = "Y", zlab = "Z",
             main = ifelse(ascending, "Ascending", "Descending"),
             type = "p"
         )

    ## Add peak center
    peak_idx <- which.min((x_coords - peak.center[1])^2 + (y_coords - peak.center[2])^2)
    rgl::points3d(
             x_coords[peak_idx], y_coords[peak_idx], z[peak_idx],
             col = "black", size = 10
         )

    ## Add best path if provided
    if (!is.null(path) && length(path) > 0) {
        ## Extract path coordinates
        path_x <- x_coords[path]
        path_y <- y_coords[path]
        path_z <- z[path]

        ## Plot path points
        rgl::points3d(
                 path_x, path_y, path_z,
                 col = ifelse(ascending, "blue", "red"),
                 size = 8
             )

        ## Plot path lines
        rgl::lines3d(
                 path_x, path_y, path_z,
                 col = ifelse(ascending, "blue", "red"),
                 lwd = 3
             )

        ##  Highlight start and end
        rgl::points3d(
                 path_x[1], path_y[1], path_z[1],
                 col = "green", size = 12
             )
        rgl::points3d(
                 path_x[length(path_x)], path_y[length(path_y)], path_z[length(path_z)],
                 col = "purple", size = 12
             )
    }
}

#' Test Local Extrema Detection with Graph Gaussian Mixture
#'
#' @description
#' Creates a graph Gaussian mixture function and tests local extrema detection on it.
#'
#' @param graph_type Type of graph to create: "grid" or "random"
#' @param n_vertices Number of vertices in the graph (for random) or grid size (for grid)
#' @param n_components Number of Gaussian components in the mixture
#' @param max.radius Maximum radius for extrema detection
#' @param min.neighborhood.size Minimum neighborhood size for extrema detection
#' @param detect.maxima Logical; if TRUE, detect maxima; if FALSE, detect minima
#'
#' @return A list containing the graph, function values, and extrema detection results
#'
#' @examples
#' \dontrun{
#' # Test on a grid graph
#' result <- test.extrema.on.gaussian.mixture()
#'
#' # Test on a random graph
#' result <- test.extrema.on.gaussian.mixture("random", 100, 3)
#' }
#'
#' @export
test.extrema.on.gaussian.mixture <- function(
  graph_type = "grid",
  n_vertices = 20,  # grid size for "grid" type
  n_components = 3,
  max.radius = 2.0,
  min.neighborhood.size = 5,
  detect.maxima = TRUE
) {
  # Create graph based on type
  if (graph_type == "grid") {
    grid.size <- n_vertices
    graph <- create.grid.graph(grid.size, grid.size)
    n_total_vertices <- grid.size * grid.size
  } else if (graph_type == "random") {
    graph <- create.random.graph(n_vertices, 3)
    n_total_vertices <- n_vertices
    grid.size <- NULL
  } else {
    stop("Invalid graph_type. Must be 'grid' or 'random'")
  }

  # Randomly select centers for Gaussian components
  set.seed(42)  # For reproducibility
  centers <- sample(n_total_vertices, n_components)

  # Generate random amplitudes and sigmas
  amplitudes <- runif(n_components, 0.5, 1.0)

  # For sigmas, use values proportional to graph size
  if (graph_type == "grid") {
    max_sigma <- grid.size / 4
  } else {
    max_sigma <- sqrt(n_vertices) / 2
  }
  sigmas <- runif(n_components, max_sigma/3, max_sigma)

  # Generate Gaussian mixture function
  z <- generate.graph.gaussian.mixture(
    graph$adj.list,
    graph$weight.list,
    centers,
    amplitudes,
    sigmas
  )

  # If looking for minima, invert the function
  if (!detect.maxima) {
    z <- 1 - z
  }

  # Visualize the function
  if (graph_type == "grid") {
    visualize.grid.function(
      grid.size,
      z,
      centers,
      paste("Gaussian Mixture on Grid -", n_components, "Components")
    )
  } else {
    # For random graphs, we'd need a different visualization
    cat("Function generated on random graph with", n_components, "components\n")
    cat("Centers at vertices:", centers, "\n")
  }

  # Detect extrema
  extrema <- detect.local.extrema(
    graph$adj.list,
    graph$weight.list,
    z,
    max.radius,
    min.neighborhood.size,
    detect.maxima
  )

  # Compare detected extrema with known centers
  if (length(extrema$vertices) > 0) {
    cat("\nDetected", length(extrema$vertices),
        ifelse(detect.maxima, "maxima", "minima"), "\n")

    # For maxima detection, centers should correspond to extrema
    if (detect.maxima) {
      # Check if centers were detected as extrema
      centers_detected <- intersect(centers, extrema$vertices)
      cat(length(centers_detected), "out of", n_components,
          "true centers detected as extrema\n")

      # Show which centers were detected
      for (i in seq_along(centers)) {
        if (centers[i] %in% extrema$vertices) {
          idx <- which(extrema$vertices == centers[i])
          cat("  Center", i, "(vertex", centers[i], ") detected with value",
              extrema$values[idx], "and neighborhood size",
              extrema$neighborhood_sizes[idx], "\n")
        } else {
          cat("  Center", i, "(vertex", centers[i], ") not detected\n")
        }
      }

      # Show additional detected extrema
      extra_extrema <- setdiff(extrema$vertices, centers)
      if (length(extra_extrema) > 0) {
        cat(length(extra_extrema), "additional extrema detected\n")
      }
    } else {
      # For minima detection, extrema should be away from centers
      cat("When detecting minima, extrema should be away from the component centers\n")
    }
  } else {
    cat("No extrema detected\n")
  }

  # Return results
  return(list(
    graph = graph,
    grid.size = grid.size,
    centers = centers,
    amplitudes = amplitudes,
    sigmas = sigmas,
    z = z,
    extrema = extrema
  ))
}

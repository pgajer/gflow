#' Visualize a Function on a Grid Graph
#'
#' Creates multiple visualizations of a function defined on a grid graph including
#' heatmap with contours, 3D perspective plot, and optionally an interactive 3D plot.
#'
#' @param grid.size Integer; The size of the square grid (grid.size x grid.size)
#' @param z Numeric vector; Function values at each vertex of the grid graph,
#'          ordered row-wise (length = grid.size^2)
#' @param centers Optional integer vector; Indices of special vertices to highlight
#'               (e.g., local maxima or centers)
#' @param title Character string; Title for the plots (default: "Function on Grid Graph")
#'
#' @return Invisibly returns the input z values
#'
#' @details
#' The function creates a side-by-side visualization with:
#' \itemize{
#'   \item Left panel: Heatmap with contour lines and optional center points
#'   \item Right panel: 3D perspective plot of the surface
#' }
#'
#' If the rgl package is available, an additional interactive 3D visualization
#' is created in a separate window.
#'
#' The grid vertices are numbered from 1 to grid.size^2, going row by row.
#' The function values in z should be ordered to match this numbering.
#'
#' @examples
#' \dontrun{
#' # Create a simple example with a Gaussian-like function on a 20x20 grid
#' grid.size <- 20
#' n_vertices <- grid.size^2
#'
#' # Generate coordinates
#' x <- rep(1:grid.size, grid.size) / grid.size
#' y <- rep(1:grid.size, each=grid.size) / grid.size
#'
#' # Create a function with two peaks
#' z <- exp(-10*((x-0.3)^2 + (y-0.3)^2)) + 0.5*exp(-8*((x-0.7)^2 + (y-0.6)^2))
#'
#' # Find local maxima (simple approach)
#' centers <- which(z > 0.9)
#'
#' # Visualize
#' visualize.grid.function(grid.size, z, centers)
#' }
#'
#' @importFrom graphics layout image contour points text persp par
#' @importFrom grDevices heat.colors
#' @export
visualize.grid.function <- function(grid.size, z, centers = NULL, title = "Function on Grid Graph") {
  # Generate coordinates
  x_coords <- rep(1:grid.size, grid.size) / grid.size
  y_coords <- rep(1:grid.size, each=grid.size) / grid.size
  
  # Convert vertex indices to 2D coordinates
  vertex_to_coords <- function(vertices) {
    x <- (vertices - 1) %% grid.size + 1
    y <- ceiling(vertices / grid.size)
    return(list(x = x / grid.size, y = y / grid.size))
  }
  
  # Get center coordinates if provided
  if (!is.null(centers)) {
    center_coords <- vertex_to_coords(centers)
  }
  
  # Create visualization layout
  layout(matrix(c(1, 2), 1, 2))
  
  # Plot 1: Heatmap with contours
  par(mar = c(4, 4, 2, 1))
  
  # Reshape z for image
  z_matrix <- matrix(z, nrow = grid.size, byrow = FALSE)
  
  # Create heatmap
  image(
    1:grid.size / grid.size,
    1:grid.size / grid.size,
    z_matrix,
    col = heat.colors(100),
    main = title,
    xlab = "X", ylab = "Y"
  )
  
  # Add contour lines
  contour(
    1:grid.size / grid.size,
    1:grid.size / grid.size,
    z_matrix,
    add = TRUE,
    col = "black"
  )
  
  # Add center points if provided
  if (!is.null(centers)) {
    points(
      center_coords$x,
      center_coords$y,
      pch = 19,
      col = "blue",
      cex = 1.5
    )
    
    # Add labels
    text(
      center_coords$x,
      center_coords$y,
      labels = paste("Center", seq_along(centers)),
      pos = 3,
      offset = 0.7,
      cex = 0.8
    )
  }
  
  # Plot 2: 3D perspective plot
  par(mar = c(4, 4, 2, 1))
  
  # Create perspective plot
  persp(
    1:grid.size / grid.size,
    1:grid.size / grid.size,
    z_matrix,
    theta = 30, phi = 30,
    expand = 0.7,
    col = "lightblue",
    shade = 0.5,
    main = "3D Perspective",
    xlab = "X", ylab = "Y", zlab = "Value"
  )
  
  # Reset layout
  layout(1)
  
  # 3D visualization if rgl is available
  if (requireNamespace("rgl", quietly = TRUE)) {
    # Open 3D device
    rgl::open3d()
    
    # Plot surface
    rgl::plot3d(
      x_coords, y_coords, z,
      col = heat.colors(length(z))[rank(z)],
      size = 3,
      xlab = "X", ylab = "Y", zlab = "Z",
      main = title,
      type = "p"
    )
    
    # Add center points if provided
    if (!is.null(centers)) {
      for (i in seq_along(centers)) {
        # Get center coordinates and value
        center_x <- center_coords$x[i]
        center_y <- center_coords$y[i]
        center_z <- z[centers[i]]
        
        # Draw sphere at center
        rgl::spheres3d(
          center_x, center_y, center_z,
          radius = 0.02,
          color = "blue"
        )
      }
    }
  }
  
  invisible(z)
}